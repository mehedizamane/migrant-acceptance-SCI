#!/usr/bin/env python3
"""Construct the country and ADM1 SCI measures used in the paper.

Country-to-country SCI is constructed separately for the August 2020 and 2026
releases. Regional SCI is constructed only from the 2026 GADM1-to-country
file. Meta scales releases and geographic products independently, so this
script never treats raw values from different products as commensurate.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path
from zipfile import ZipFile

import geopandas as gpd
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import DERIVED, RAW, ensure_directories, require_columns


# August 2020 release, from the Bailey et al. (2021) replication archive (Raw_data/SCI/).
COUNTRY_2020_PATH = RAW / "meta" / "country_country_aug2020.tsv"
COUNTRY_2026_PATH = RAW / "meta" / "country_2026.csv"
REGIONAL_2026_ARCHIVE = RAW / "meta" / "all_region_to_country.zip"
REGIONAL_2026_MEMBER = "gadm1_to_country.csv"
POPULATION_PATH = RAW / "population" / "gadm1_ghs_pop_2025_population.csv"
GADM_PATH = RAW / "geography" / "gadm1_metadata.gpkg"
CROSSWALK_PATH = RAW / "geography" / "v1" / "gallup_to_gadm1_crosswalk_v1.csv"
COUNTRY_METADATA_PATH = RAW / "controls" / "country_metadata_reference.csv"

EXPECTED = {
    "country_2020_origins": 185,
    "country_2020_foreign_rows": 34_040,
    "country_2020_partners": 184,
    "country_2026_origins": 178,
    "country_2026_foreign_rows": 31_506,
    "country_2026_partners": 177,
    "regional_rows": 563_904,
    "regional_origins": 3_168,
    "regional_destinations": 178,
    "regional_foreign_destinations": 177,
}


def load_country_release(release: int) -> pd.DataFrame:
    if release == 2020:
        frame = pd.read_csv(
            COUNTRY_2020_PATH,
            sep="\t",
            keep_default_na=False,
            dtype={"user_loc": str, "fr_loc": str},
        ).rename(columns={"user_loc": "origin", "fr_loc": "destination"})
    elif release == 2026:
        frame = pd.read_csv(
            COUNTRY_2026_PATH,
            keep_default_na=False,
            dtype={"user_country": str, "friend_country": str},
        ).rename(
            columns={"user_country": "origin", "friend_country": "destination"}
        )
    else:
        raise ValueError(release)

    require_columns(frame, ["origin", "destination", "scaled_sci"], f"{release} SCI")
    frame = frame[["origin", "destination", "scaled_sci"]].copy()
    frame["origin"] = frame["origin"].astype(str).str.strip()
    frame["destination"] = frame["destination"].astype(str).str.strip()
    frame["scaled_sci"] = pd.to_numeric(frame["scaled_sci"], errors="raise")
    if frame["scaled_sci"].le(0).any():
        raise RuntimeError(f"{release} SCI contains nonpositive scaled_sci values")
    if frame.duplicated(["origin", "destination"]).any():
        raise RuntimeError(f"{release} SCI contains duplicate directed dyads")
    frame["release"] = release
    frame["log_sci"] = np.log(frame["scaled_sci"].astype(float))
    return frame


def validate_country_release(frame: pd.DataFrame, release: int) -> pd.DataFrame:
    origins = frame["origin"].nunique()
    foreign = frame.loc[frame["origin"].ne(frame["destination"])].copy()
    self_links = frame.loc[frame["origin"].eq(frame["destination"])]
    partners = foreign.groupby("origin")["destination"].nunique()
    if origins != EXPECTED[f"country_{release}_origins"]:
        raise RuntimeError(f"{release}: expected origin count changed: {origins}")
    if len(foreign) != EXPECTED[f"country_{release}_foreign_rows"]:
        raise RuntimeError(f"{release}: expected foreign dyad count changed: {len(foreign)}")
    if len(self_links) != origins or self_links["origin"].nunique() != origins:
        raise RuntimeError(f"{release}: expected one self-link per origin")
    expected_partners = EXPECTED[f"country_{release}_partners"]
    if len(partners) != origins or not partners.eq(expected_partners).all():
        raise RuntimeError(f"{release}: not every origin has {expected_partners} partners")
    if "NA" not in set(frame["origin"]):
        raise RuntimeError(f"{release}: Namibia ISO-2 code NA was not preserved")
    return foreign


def build_country_measures() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    release_frames = {release: load_country_release(release) for release in (2020, 2026)}
    foreign = {
        release: validate_country_release(frame, release)
        for release, frame in release_frames.items()
    }

    summaries: dict[int, pd.DataFrame] = {}
    for release in (2020, 2026):
        summaries[release] = (
            foreign[release]
            .groupby("origin", as_index=False)
            .agg(
                **{
                    f"global_sci_{release}": ("log_sci", "mean"),
                    f"global_sci_partner_log_variance_{release}": (
                        "log_sci",
                        lambda values: float(np.var(values, ddof=0)),
                    ),
                    f"foreign_partner_count_{release}": ("destination", "nunique"),
                }
            )
            .rename(columns={"origin": "iso2"})
        )
        summaries[release][f"global_sci_rank_{release}"] = summaries[release][
            f"global_sci_{release}"
        ].rank(ascending=False, method="min")

    measures = summaries[2026].merge(
        summaries[2020], on="iso2", how="outer", validate="one_to_one"
    )
    common_codes = sorted(set(summaries[2020]["iso2"]).intersection(summaries[2026]["iso2"]))
    if len(common_codes) != 166:
        raise RuntimeError(f"Expected 166 common country origins; observed {len(common_codes)}")

    common_dyads = foreign[2020][["origin", "destination", "log_sci"]].rename(
        columns={"log_sci": "log_sci_2020"}
    ).merge(
        foreign[2026][["origin", "destination", "log_sci"]].rename(
            columns={"log_sci": "log_sci_2026"}
        ),
        on=["origin", "destination"],
        how="inner",
        validate="one_to_one",
    )
    counts = common_dyads.groupby("origin")["destination"].nunique()
    if len(common_dyads) != 27_390 or len(counts) != 166 or not counts.eq(165).all():
        raise RuntimeError("Common-release dyads are not a complete 166 x 165 directed panel")
    harmonized = (
        common_dyads.groupby("origin", as_index=False)[["log_sci_2020", "log_sci_2026"]]
        .mean()
        .rename(
            columns={
                "origin": "iso2",
                "log_sci_2020": "global_sci_2020_common_partners",
                "log_sci_2026": "global_sci_2026_common_partners",
            }
        )
    )
    measures = measures.merge(harmonized, on="iso2", how="left", validate="one_to_one")

    scaling_rows: list[dict[str, object]] = []
    for release in (2020, 2026):
        source = f"global_sci_{release}"
        target = f"global_sci_{release}_all_origins_z"
        reference = measures[source].notna()
        mean = float(measures.loc[reference, source].mean())
        sd = float(measures.loc[reference, source].std(ddof=0))
        measures[target] = (measures[source] - mean) / sd
        scaling_rows.append(
            {
                "release": release,
                "source_variable": source,
                "standardized_variable": target,
                "reference_countries": int(reference.sum()),
                "foreign_partners_per_origin": EXPECTED[f"country_{release}_partners"],
                "reference_mean_log_sci": mean,
                "reference_sd_log_sci": sd,
                "definition": (
                    "Equal-destination-country mean of natural-log scaled SCI over all "
                    "non-self partners in this release; standardized over every origin "
                    "available in that release"
                ),
            }
        )

    for release in (2020, 2026):
        source = f"global_sci_{release}"
        target = f"global_sci_{release}_z"
        reference = measures["iso2"].isin(common_codes) & measures[source].notna()
        mean = float(measures.loc[reference, source].mean())
        sd = float(measures.loc[reference, source].std(ddof=0))
        measures[target] = (measures[source] - mean) / sd
        scaling_rows.append(
            {
                "release": release,
                "source_variable": source,
                "standardized_variable": target,
                "reference_countries": int(reference.sum()),
                "foreign_partners_per_origin": EXPECTED[f"country_{release}_partners"],
                "reference_mean_log_sci": mean,
                "reference_sd_log_sci": sd,
                "definition": (
                    "Equal-destination-country mean of natural-log scaled SCI over all "
                    "non-self partners in this release; standardized over the 166 "
                    "origins common to the 2020 and 2026 releases"
                ),
            }
        )

    source = "global_sci_2026_common_partners"
    mean = float(measures.loc[measures["iso2"].isin(common_codes), source].mean())
    sd = float(measures.loc[measures["iso2"].isin(common_codes), source].std(ddof=0))
    measures[f"{source}_z"] = (measures[source] - mean) / sd
    scaling_rows.append(
        {
            "release": "2026_common_2020_partner_set",
            "source_variable": source,
            "standardized_variable": f"{source}_z",
            "reference_countries": 166,
            "foreign_partners_per_origin": 165,
            "reference_mean_log_sci": mean,
            "reference_sd_log_sci": sd,
            "definition": (
                "2026 equal-country mean restricted to the 165 destinations present in "
                "the 2020 release; standardized over the 166 common origins"
            ),
        }
    )

    metadata = pd.read_csv(COUNTRY_METADATA_PATH, keep_default_na=False, low_memory=False)
    metadata_cols = [column for column in ["iso2", "iso3", "country_name", "region"] if column in metadata]
    metadata = metadata[metadata_cols].drop_duplicates("iso2")
    measures = measures.merge(metadata, on="iso2", how="left", validate="one_to_one")
    measures["is_common_release_origin"] = measures["iso2"].isin(common_codes)
    measures["global_sci_rank_change_2026_minus_2020"] = (
        measures["global_sci_rank_2026"] - measures["global_sci_rank_2020"]
    )

    common = measures.loc[measures["is_common_release_origin"]].copy()
    country_pearson = pearsonr(common["global_sci_2020"], common["global_sci_2026"])
    country_spearman = spearmanr(common["global_sci_2020"], common["global_sci_2026"])
    dyad_pearson = pearsonr(common_dyads["log_sci_2020"], common_dyads["log_sci_2026"])
    dyad_spearman = spearmanr(common_dyads["log_sci_2020"], common_dyads["log_sci_2026"])
    stability = pd.DataFrame(
        [
            {
                "level": "Country GlobalSCI",
                "n": len(common),
                "pearson_r": country_pearson.statistic,
                "pearson_p": country_pearson.pvalue,
                "spearman_rho": country_spearman.statistic,
                "spearman_p": country_spearman.pvalue,
            },
            {
                "level": "Common directed country dyads",
                "n": len(common_dyads),
                "pearson_r": dyad_pearson.statistic,
                "pearson_p": dyad_pearson.pvalue,
                "spearman_rho": dyad_spearman.statistic,
                "spearman_p": dyad_spearman.pvalue,
            },
        ]
    )

    measures = measures.sort_values("iso2").reset_index(drop=True)
    measures.to_csv(DERIVED / "country_global_sci_2020_2026.csv", index=False)
    common_dyads.to_csv(DERIVED / "country_sci_2020_2026_common_dyads.csv", index=False)
    pd.DataFrame(scaling_rows).to_csv(DERIVED / "country_sci_release_scaling.csv", index=False)
    stability.to_csv(DERIVED / "country_sci_release_stability.csv", index=False)
    return measures, pd.DataFrame(scaling_rows), stability


def accepted_adm1_survey_country_map() -> pd.Series:
    crosswalk = pd.read_csv(CROSSWALK_PATH, keep_default_na=False, dtype=str)
    require_columns(crosswalk, ["iso2", "adm1_gid", "match_status"], "Gallup-GADM1 crosswalk")
    accepted = crosswalk.loc[crosswalk["match_status"].eq("accepted"), ["iso2", "adm1_gid"]].copy()
    accepted = accepted.loc[accepted["adm1_gid"].str.strip().ne("")]
    conflicts = accepted.groupby("adm1_gid")["iso2"].nunique()
    if conflicts.gt(1).any():
        raise RuntimeError("An accepted ADM1 region maps to multiple Gallup country codes")
    return accepted.drop_duplicates("adm1_gid").set_index("adm1_gid")["iso2"]


def build_regional_measures(country_measures: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    with ZipFile(REGIONAL_2026_ARCHIVE) as archive:
        if REGIONAL_2026_MEMBER not in archive.namelist():
            raise RuntimeError(f"Regional SCI archive lacks {REGIONAL_2026_MEMBER}")
        with archive.open(REGIONAL_2026_MEMBER) as handle:
            dyads = pd.read_csv(
                handle,
                keep_default_na=False,
                dtype={
                    "user_country": str,
                    "friend_country": str,
                    "user_region": str,
                    "friend_region": str,
                },
            )
    require_columns(
        dyads,
        ["user_country", "friend_country", "user_region", "friend_region", "scaled_sci"],
        "2026 GADM1-to-country SCI",
    )
    if len(dyads) != EXPECTED["regional_rows"]:
        raise RuntimeError(f"Expected 563,904 regional rows; observed {len(dyads)}")
    if dyads["user_region"].nunique() != EXPECTED["regional_origins"]:
        raise RuntimeError("Expected 3,168 regional origins")
    per_origin = dyads.groupby("user_region")["friend_country"].nunique()
    if not per_origin.eq(EXPECTED["regional_destinations"]).all():
        raise RuntimeError("Every regional origin must have 178 country destinations")
    dyads["scaled_sci"] = pd.to_numeric(dyads["scaled_sci"], errors="raise")
    if dyads["scaled_sci"].le(0).any():
        raise RuntimeError("Regional SCI includes nonpositive values")
    dyads["log_sci"] = np.log(dyads["scaled_sci"].astype(float))
    foreign = dyads.loc[dyads["user_country"].ne(dyads["friend_country"])].copy()
    own_country = dyads.loc[dyads["user_country"].eq(dyads["friend_country"])].copy()
    counts = foreign.groupby("user_region")["friend_country"].nunique()
    if len(own_country) != EXPECTED["regional_origins"] or not counts.eq(177).all():
        raise RuntimeError("Expected one own-country and 177 foreign-country rows per ADM1")

    region = (
        foreign.groupby(["user_region", "user_country"], as_index=False)
        .agg(
            regional_sci_2026=("log_sci", "mean"),
            regional_sci_partner_log_variance_2026=(
                "log_sci",
                lambda values: float(np.var(values, ddof=0)),
            ),
            regional_sci_foreign_mean_raw_2026=("scaled_sci", "mean"),
            foreign_destination_count=("friend_country", "nunique"),
        )
        .rename(columns={"user_region": "adm1_gid", "user_country": "sci_country_iso2"})
    )
    own_country = own_country[["user_region", "scaled_sci", "log_sci"]].rename(
        columns={
            "user_region": "adm1_gid",
            "scaled_sci": "own_country_raw_sci_2026",
            "log_sci": "own_country_log_sci_2026",
        }
    )
    region = region.merge(own_country, on="adm1_gid", how="left", validate="one_to_one")

    population = pd.read_csv(POPULATION_PATH, keep_default_na=False)
    population["ghs_pop_2025"] = pd.to_numeric(population["ghs_pop_2025"], errors="coerce")
    region = region.merge(population, on="adm1_gid", how="left", validate="one_to_one")

    gadm = gpd.read_file(GADM_PATH, layer="gadm1_metadata")
    require_columns(gadm, ["GID_1", "GID_0", "COUNTRY", "NAME_1", "geometry"], "GADM metadata")
    projected = gadm[["GID_1", "geometry"]].to_crs(8857)
    area = pd.DataFrame(
        {
            "adm1_gid": projected["GID_1"].astype(str),
            "adm1_area_km2": projected.geometry.area.to_numpy(float) / 1_000_000,
        }
    )
    metadata = gadm[["GID_1", "GID_0", "COUNTRY", "NAME_1"]].rename(
        columns={
            "GID_1": "adm1_gid",
            "GID_0": "gadm_country_code",
            "COUNTRY": "gadm_country_name",
            "NAME_1": "adm1_name",
        }
    )
    region = region.merge(metadata, on="adm1_gid", how="left", validate="one_to_one")
    region = region.merge(area, on="adm1_gid", how="left", validate="one_to_one")
    survey_map = accepted_adm1_survey_country_map()
    region["analysis_iso2"] = region["adm1_gid"].map(survey_map).fillna(region["sci_country_iso2"])
    region["population_valid"] = region["ghs_pop_2025"].notna() & region["ghs_pop_2025"].gt(0)
    region["population_weight_within_country"] = np.nan
    valid = region["population_valid"]
    region.loc[valid, "population_weight_within_country"] = (
        region.loc[valid, "ghs_pop_2025"]
        / region.loc[valid].groupby("analysis_iso2")["ghs_pop_2025"].transform("sum")
    )

    population_country = (
        region.loc[valid]
        .assign(
            weighted_sci=lambda frame: frame["regional_sci_2026"]
            * frame["population_weight_within_country"]
        )
        .groupby("analysis_iso2")["weighted_sci"]
        .sum()
    )
    equal_country = region.groupby("analysis_iso2")["regional_sci_2026"].mean()
    region["country_average_regional_sci_2026_population_weighted"] = region[
        "analysis_iso2"
    ].map(population_country)
    region["country_average_regional_sci_2026_equal_adm1"] = region["analysis_iso2"].map(
        equal_country
    )
    region["regional_sci_deviation_from_country_2026"] = (
        region["regional_sci_2026"]
        - region["country_average_regional_sci_2026_population_weighted"]
    )
    region["regional_sci_deviation_from_equal_adm1_country_2026"] = (
        region["regional_sci_2026"]
        - region["country_average_regional_sci_2026_equal_adm1"]
    )

    variance_frame = region.loc[valid].copy()
    country_count = variance_frame["analysis_iso2"].nunique()
    variance_frame["equal_country_population_weight"] = (
        variance_frame["population_weight_within_country"] / country_count
    )
    grand_mean = float(
        np.sum(
            variance_frame["equal_country_population_weight"]
            * variance_frame["regional_sci_2026"]
        )
    )
    within_variance = float(
        np.sum(
            variance_frame["equal_country_population_weight"]
            * np.square(variance_frame["regional_sci_deviation_from_country_2026"])
        )
    )
    between_variance = float(
        np.sum(
            variance_frame["equal_country_population_weight"]
            * np.square(
                variance_frame["country_average_regional_sci_2026_population_weighted"]
                - grand_mean
            )
        )
    )
    total_variance = float(
        np.sum(
            variance_frame["equal_country_population_weight"]
            * np.square(variance_frame["regional_sci_2026"] - grand_mean)
        )
    )
    if not np.isclose(total_variance, between_variance + within_variance, atol=1e-12):
        raise RuntimeError("Regional SCI variance decomposition is not additive")
    within_sd = math.sqrt(within_variance)
    total_sd = math.sqrt(total_variance)
    region["regional_sci_deviation_within_country_z"] = (
        region["regional_sci_deviation_from_country_2026"] / within_sd
    )
    region["regional_sci_deviation_total_sd_units"] = (
        region["regional_sci_deviation_from_country_2026"] / total_sd
    )

    country_rows = []
    for iso2, part in region.groupby("analysis_iso2", sort=True):
        country_rows.append(
            {
                "iso2": iso2,
                "country_average_regional_sci_2026_population_weighted": population_country.get(
                    iso2, np.nan
                ),
                "country_average_regional_sci_2026_equal_adm1": equal_country.get(iso2, np.nan),
                "n_adm1_regions": part["adm1_gid"].nunique(),
                "n_adm1_regions_with_population": int(part["population_valid"].sum()),
                "ghs_pop_2025_total": part.loc[part["population_valid"], "ghs_pop_2025"].sum(),
            }
        )
    regional_country = pd.DataFrame(country_rows)
    between_reference = regional_country[
        "country_average_regional_sci_2026_population_weighted"
    ].dropna()
    between_mean = float(between_reference.mean())
    between_sd = float(between_reference.std(ddof=0))
    regional_country["country_average_regional_sci_2026_population_weighted_z"] = (
        regional_country["country_average_regional_sci_2026_population_weighted"]
        - between_mean
    ) / between_sd
    equal_reference = regional_country["country_average_regional_sci_2026_equal_adm1"]
    equal_mean = float(equal_reference.mean())
    equal_sd = float(equal_reference.std(ddof=0))
    regional_country["country_average_regional_sci_2026_equal_adm1_z"] = (
        equal_reference - equal_mean
    ) / equal_sd

    variance = pd.DataFrame(
        [
            {
                "component": "Between countries",
                "variance": between_variance,
                "standard_deviation": math.sqrt(between_variance),
                "share_of_total_variance": between_variance / total_variance,
                "countries": country_count,
                "adm1_regions": len(variance_frame),
            },
            {
                "component": "Within countries",
                "variance": within_variance,
                "standard_deviation": within_sd,
                "share_of_total_variance": within_variance / total_variance,
                "countries": country_count,
                "adm1_regions": len(variance_frame),
            },
            {
                "component": "Total",
                "variance": total_variance,
                "standard_deviation": total_sd,
                "share_of_total_variance": 1.0,
                "countries": country_count,
                "adm1_regions": len(variance_frame),
            },
        ]
    )

    if not np.allclose(
        region.loc[valid, "regional_sci_2026"],
        region.loc[valid, "country_average_regional_sci_2026_population_weighted"]
        + region.loc[valid, "regional_sci_deviation_from_country_2026"],
        atol=1e-12,
        rtol=0,
    ):
        raise RuntimeError("Regional score does not reconstruct from country mean plus deviation")

    # Country-file concordance is standardized/rank based because Meta scales
    # the country and regional products independently.
    sci_origin_country = (
        region.loc[region["population_valid"]]
        .assign(
            origin_weight=lambda frame: frame["ghs_pop_2025"]
            / frame.groupby("sci_country_iso2")["ghs_pop_2025"].transform("sum")
        )
        .assign(weighted=lambda frame: frame["regional_sci_2026"] * frame["origin_weight"])
        .groupby("sci_country_iso2", as_index=False)
        .agg(
            regional_to_country_sci_population_weighted_2026=("weighted", "sum"),
            n_adm1_with_population=("adm1_gid", "nunique"),
        )
        .rename(columns={"sci_country_iso2": "iso2"})
    )
    equal_origin = (
        region.groupby("sci_country_iso2", as_index=False)
        .agg(
            regional_to_country_sci_equal_adm1_2026=("regional_sci_2026", "mean"),
            n_adm1=("adm1_gid", "nunique"),
        )
        .rename(columns={"sci_country_iso2": "iso2"})
    )
    concordance = equal_origin.merge(
        sci_origin_country, on="iso2", how="left", validate="one_to_one"
    ).merge(
        country_measures[["iso2", "global_sci_2026", "global_sci_2026_z"]],
        on="iso2",
        how="left",
        validate="one_to_one",
    )
    for source in [
        "regional_to_country_sci_equal_adm1_2026",
        "regional_to_country_sci_population_weighted_2026",
    ]:
        concordance[f"{source}_z"] = (
            concordance[source] - concordance[source].mean()
        ) / concordance[source].std(ddof=0)
        concordance[f"{source}_rank"] = concordance[source].rank(
            ascending=False, method="min"
        )
    concordance["global_sci_2026_rank"] = concordance["global_sci_2026"].rank(
        ascending=False, method="min"
    )

    region["adm1_population_density_per_km2"] = (
        region["ghs_pop_2025"] / region["adm1_area_km2"]
    )
    region["adm1_log_population"] = np.log1p(region["ghs_pop_2025"])
    region["adm1_log_population_density"] = np.log1p(
        region["adm1_population_density_per_km2"]
    )
    region = region.sort_values(["analysis_iso2", "adm1_gid"]).reset_index(drop=True)
    regional_country = regional_country.sort_values("iso2").reset_index(drop=True)
    region.to_csv(DERIVED / "regional_sci_adm1_2026.csv", index=False)
    regional_country.to_csv(DERIVED / "regional_sci_country_2026.csv", index=False)
    variance.to_csv(DERIVED / "regional_sci_variance_decomposition.csv", index=False)
    concordance.to_csv(DERIVED / "country_regional_sci_concordance.csv", index=False)
    return region, variance


def main() -> None:
    ensure_directories()
    country, _, stability = build_country_measures()
    region, variance = build_regional_measures(country)
    print(
        "Constructed country SCI for 2020/2026 and regional SCI for 2026: "
        f"{len(country)} unique country/territory rows across releases; "
        f"{len(region)} ADM1 rows."
    )
    print(stability.to_string(index=False))
    print(variance.to_string(index=False))


if __name__ == "__main__":
    main()
