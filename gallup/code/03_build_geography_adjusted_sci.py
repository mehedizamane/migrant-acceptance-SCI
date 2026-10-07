#!/usr/bin/env python3
"""Construct geography-adjusted country and ADM1 foreign-SCI measures.

The country measure is the origin-country mean residual from a dyadic model of
log SCI on distance, contiguity, shared languages, colonial history, and
destination-country fixed effects. Origin fixed effects are deliberately
omitted because they would absorb the country-level quantity being built.

The ADM1 measure uses the same geographic covariates plus origin-country and
destination-country fixed effects. It is therefore a within-country regional
sensitivity: the origin-country effects remove country-average differences.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path
from zipfile import ZipFile

import geopandas as gpd
import numpy as np
import pandas as pd
from pyproj import Geod
from scipy import sparse
from scipy.sparse.linalg import lsqr

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import DERIVED, RAW, ensure_directories, require_columns


# August 2020 release, from the Bailey et al. (2021) replication archive (Raw_data/SCI/).
COUNTRY_2020_PATH = RAW / "meta" / "country_country_aug2020.tsv"
COUNTRY_2026_PATH = RAW / "meta" / "country_2026.csv"
REGIONAL_2026_ARCHIVE = RAW / "meta" / "all_region_to_country.zip"
REGIONAL_2026_MEMBER = "gadm1_to_country.csv"
GADM_PATH = RAW / "geography" / "gadm1_metadata.gpkg"
CEPII_PATH = DERIVED / "cepii_country_dyadic_controls_2020.csv.gz"
COUNTRY_METADATA_PATH = RAW / "controls" / "country_metadata_reference.csv"
POPULATION_PATH = RAW / "population" / "gadm1_ghs_pop_2025_population.csv"

GEOGRAPHY_COLUMNS = [
    "log_distance_km",
    "contig",
    "comlang_off",
    "comlang_ethno",
    "col_dep_ever",
]


def sparse_residualize(
    frame: pd.DataFrame,
    outcome: str,
    continuous: list[str],
    fixed_effects: list[str],
) -> tuple[np.ndarray, pd.DataFrame, dict[str, float]]:
    """Fit a sparse least-squares model and return residuals and diagnostics."""
    values = frame[continuous].astype(float).to_numpy()
    means = values.mean(axis=0)
    scales = values.std(axis=0, ddof=0)
    scales[scales == 0] = 1.0
    standardized = (values - means) / scales
    pieces: list[sparse.spmatrix] = [
        sparse.csr_matrix(np.ones((len(frame), 1))),
        sparse.csr_matrix(standardized),
    ]
    parameter_names = ["intercept", *continuous]
    row_index = np.arange(len(frame))
    for fixed_effect in fixed_effects:
        categories = pd.Categorical(
            frame[fixed_effect], categories=sorted(frame[fixed_effect].unique())
        )
        matrix = sparse.csr_matrix(
            (np.ones(len(frame)), (row_index, categories.codes)),
            shape=(len(frame), len(categories.categories)),
        )[:, 1:]
        pieces.append(matrix)
        parameter_names.extend(
            f"{fixed_effect}[{level}]" for level in categories.categories[1:]
        )
    design = sparse.hstack(pieces, format="csr")
    y = frame[outcome].astype(float).to_numpy()
    solution = lsqr(design, y, atol=1e-10, btol=1e-10, iter_lim=5_000)
    beta = solution[0]
    fitted = np.asarray(design @ beta).ravel()
    residual = y - fitted
    sse = float(np.sum(residual**2))
    sst = float(np.sum((y - y.mean()) ** 2))
    coefficient_rows = []
    for offset, name in enumerate(continuous, start=1):
        coefficient_rows.append(
            {
                "term": name,
                "coefficient_original_scale": beta[offset] / scales[offset - 1],
                "coefficient_standardized_predictor": beta[offset],
                "predictor_mean": means[offset - 1],
                "predictor_sd": scales[offset - 1],
            }
        )
    diagnostics = {
        "n_observations": float(len(frame)),
        "n_parameters": float(design.shape[1]),
        "r_squared": float(1 - sse / sst),
        "rmse": math.sqrt(sse / len(frame)),
        "lsqr_iterations": float(solution[2]),
        "lsqr_condition_estimate": float(solution[6]),
    }
    if not np.isfinite(residual).all():
        raise RuntimeError("Geography-adjustment model produced nonfinite residuals")
    return residual, pd.DataFrame(coefficient_rows), diagnostics


def load_country_sci(release: int) -> pd.DataFrame:
    if release == 2020:
        frame = pd.read_csv(
            COUNTRY_2020_PATH,
            sep="\t",
            keep_default_na=False,
            dtype={"user_loc": str, "fr_loc": str},
        ).rename(columns={"user_loc": "origin_iso2", "fr_loc": "destination_iso2"})
    elif release == 2026:
        frame = pd.read_csv(
            COUNTRY_2026_PATH,
            keep_default_na=False,
            dtype={"user_country": str, "friend_country": str},
        ).rename(
            columns={
                "user_country": "origin_iso2",
                "friend_country": "destination_iso2",
            }
        )
    else:
        raise ValueError(release)
    frame["scaled_sci"] = pd.to_numeric(frame["scaled_sci"], errors="raise")
    frame = frame.loc[frame["origin_iso2"].ne(frame["destination_iso2"])].copy()
    frame["log_sci"] = np.log(frame["scaled_sci"])
    return frame


def load_country_metadata() -> pd.DataFrame:
    metadata = pd.read_csv(COUNTRY_METADATA_PATH, keep_default_na=False)
    require_columns(metadata, ["iso2", "iso3"], "country metadata")
    return metadata[["iso2", "iso3"]].drop_duplicates("iso2")


def build_country_adjusted() -> tuple[pd.DataFrame, pd.DataFrame]:
    metadata = load_country_metadata()
    iso3 = metadata.set_index("iso2")["iso3"].to_dict()
    cepii = pd.read_csv(CEPII_PATH, keep_default_na=False)
    require_columns(
        cepii,
        ["iso3_o", "iso3_d", *GEOGRAPHY_COLUMNS],
        "CEPII dyadic controls",
    )
    release_dyads = {release: load_country_sci(release) for release in (2020, 2026)}
    common_origins = set(release_dyads[2020]["origin_iso2"]).intersection(
        release_dyads[2026]["origin_iso2"]
    )
    cepii_origins = set(cepii["iso3_o"].astype(str))
    cepii_destinations = set(cepii["iso3_d"].astype(str))
    common_supported = sorted(
        code
        for code in common_origins
        if iso3.get(code, "") in cepii_origins
        and iso3.get(code, "") in cepii_destinations
    )
    if len(common_supported) != 165:
        raise RuntimeError(
            "Expected 165 common country origins supported by both SCI releases and CEPII; "
            f"observed {len(common_supported)}"
        )
    outputs = []
    fit_rows = []
    for release in (2020, 2026):
        dyads = release_dyads[release].loc[
            release_dyads[release]["origin_iso2"].isin(common_supported)
            & release_dyads[release]["destination_iso2"].isin(common_supported)
        ].copy()
        counts = dyads.groupby("origin_iso2")["destination_iso2"].nunique()
        if len(dyads) != 165 * 164 or len(counts) != 165 or not counts.eq(164).all():
            raise RuntimeError(
                f"{release} geography adjustment lacks the complete 165 x 164 common-country panel"
            )
        dyads["iso3_o"] = dyads["origin_iso2"].map(iso3)
        dyads["iso3_d"] = dyads["destination_iso2"].map(iso3)
        merged = dyads.merge(
            cepii,
            on=["iso3_o", "iso3_d"],
            how="left",
            indicator=True,
            validate="many_to_one",
        )
        matched = merged.loc[merged["_merge"].eq("both")].drop(columns="_merge").copy()
        match_rate = len(matched) / len(merged)
        if match_rate != 1:
            raise RuntimeError(
                f"CEPII failed to match the complete common-country panel for {release}: "
                f"{match_rate:.1%}"
            )
        residual, coefficients, diagnostics = sparse_residualize(
            matched,
            outcome="log_sci",
            continuous=GEOGRAPHY_COLUMNS,
            fixed_effects=["iso3_d"],
        )
        matched["geography_adjusted_log_sci_residual"] = residual
        summary = (
            matched.groupby("origin_iso2", as_index=False)
            .agg(
                **{
                    f"global_sci_geography_adjusted_{release}": (
                        "geography_adjusted_log_sci_residual",
                        "mean",
                    ),
                    f"geography_adjusted_destination_count_{release}": (
                        "destination_iso2",
                        "nunique",
                    ),
                }
            )
            .rename(columns={"origin_iso2": "iso2"})
        )
        outputs.append(summary)
        coefficients.insert(0, "release", release)
        coefficients["model"] = (
            "Country log SCI on geography/history and destination-country fixed effects"
        )
        for key, value in diagnostics.items():
            coefficients[key] = value
        coefficients["match_rate"] = match_rate
        fit_rows.append(coefficients)
    country = outputs[0].merge(outputs[1], on="iso2", how="outer", validate="one_to_one")
    common = country[
        ["global_sci_geography_adjusted_2020", "global_sci_geography_adjusted_2026"]
    ].notna().all(axis=1)
    if int(common.sum()) != 165:
        raise RuntimeError("Geography-adjusted release comparison lost common origins")
    for release in (2020, 2026):
        source = f"global_sci_geography_adjusted_{release}"
        mean = float(country.loc[common, source].mean())
        sd = float(country.loc[common, source].std(ddof=0))
        country[f"{source}_z"] = (country[source] - mean) / sd
    country.to_csv(DERIVED / "country_global_sci_geography_adjusted_2020_2026.csv", index=False)
    fits = pd.concat(fit_rows, ignore_index=True)
    fits.to_csv(DERIVED / "country_sci_geography_adjustment_models.csv", index=False)
    return country, fits


def compute_centroids() -> tuple[pd.DataFrame, pd.DataFrame]:
    gadm = gpd.read_file(GADM_PATH, layer="gadm1_metadata")
    require_columns(gadm, ["GID_1", "GID_0", "geometry"], "GADM metadata")
    projected = gadm[["GID_1", "GID_0", "geometry"]].to_crs(8857)
    region_points = gpd.GeoSeries(projected.geometry.centroid, crs=8857).to_crs(4326)
    regions = pd.DataFrame(
        {
            "adm1_gid": projected["GID_1"].astype(str),
            "origin_lon": region_points.x,
            "origin_lat": region_points.y,
        }
    )
    dissolved = projected.dissolve(by="GID_0", as_index=False)
    country_points = gpd.GeoSeries(dissolved.geometry.centroid, crs=8857).to_crs(4326)
    countries = pd.DataFrame(
        {
            "iso3_d": dissolved["GID_0"].astype(str),
            "destination_lon": country_points.x,
            "destination_lat": country_points.y,
        }
    )
    return regions, countries


def build_regional_adjusted() -> tuple[pd.DataFrame, pd.DataFrame]:
    with ZipFile(REGIONAL_2026_ARCHIVE) as archive, archive.open(
        REGIONAL_2026_MEMBER
    ) as handle:
        dyads = pd.read_csv(
            handle,
            keep_default_na=False,
            dtype={
                "user_country": str,
                "friend_country": str,
                "user_region": str,
            },
        )
    dyads = dyads.loc[dyads["user_country"].ne(dyads["friend_country"])].copy()
    if len(dyads) != 3_168 * 177:
        raise RuntimeError("Unexpected number of foreign ADM1-to-country dyads")
    dyads["scaled_sci"] = pd.to_numeric(dyads["scaled_sci"], errors="raise")
    dyads["log_sci"] = np.log(dyads["scaled_sci"])
    dyads = dyads.rename(
        columns={
            "user_region": "adm1_gid",
            "user_country": "origin_iso2",
            "friend_country": "destination_iso2",
        }
    )
    metadata = load_country_metadata()
    iso3 = metadata.set_index("iso2")["iso3"].to_dict()
    dyads["iso3_o"] = dyads["origin_iso2"].map(iso3)
    dyads["iso3_d"] = dyads["destination_iso2"].map(iso3)
    cepii = pd.read_csv(CEPII_PATH, keep_default_na=False).drop(
        columns=["log_distance_km"], errors="ignore"
    )
    merged = dyads.merge(
        cepii,
        on=["iso3_o", "iso3_d"],
        how="left",
        indicator=True,
        validate="many_to_one",
    )
    matched = merged.loc[merged["_merge"].eq("both")].drop(columns="_merge").copy()
    match_rate = len(matched) / len(merged)
    if match_rate < 0.98:
        raise RuntimeError(f"CEPII matched only {match_rate:.1%} of regional dyads")
    regions, countries = compute_centroids()
    matched = matched.merge(regions, on="adm1_gid", how="left", validate="many_to_one")
    matched = matched.merge(countries, on="iso3_d", how="left", validate="many_to_one")
    geod = Geod(ellps="WGS84")
    _, _, distance_metres = geod.inv(
        matched["origin_lon"].to_numpy(),
        matched["origin_lat"].to_numpy(),
        matched["destination_lon"].to_numpy(),
        matched["destination_lat"].to_numpy(),
    )
    matched["distance_km"] = np.asarray(distance_metres) / 1_000
    valid_distance = np.isfinite(matched["distance_km"]) & matched["distance_km"].gt(0)
    matched = matched.loc[valid_distance].copy()
    matched["log_distance_km"] = np.log(matched["distance_km"])
    residual, coefficients, diagnostics = sparse_residualize(
        matched,
        outcome="log_sci",
        continuous=GEOGRAPHY_COLUMNS,
        fixed_effects=["iso3_o", "iso3_d"],
    )
    matched["geography_adjusted_log_sci_residual"] = residual
    region = (
        matched.groupby(["adm1_gid", "origin_iso2"], as_index=False)
        .agg(
            regional_sci_geography_adjusted_2026=(
                "geography_adjusted_log_sci_residual",
                "mean",
            ),
            geography_adjusted_destination_count_2026=("destination_iso2", "nunique"),
        )
        .rename(columns={"origin_iso2": "sci_country_iso2"})
    )
    population = pd.read_csv(POPULATION_PATH, keep_default_na=False)
    population["ghs_pop_2025"] = pd.to_numeric(population["ghs_pop_2025"], errors="coerce")
    region = region.merge(population, on="adm1_gid", how="left", validate="one_to_one")
    valid_pop = region["ghs_pop_2025"].notna() & region["ghs_pop_2025"].gt(0)
    region["population_weight_within_country"] = np.nan
    region.loc[valid_pop, "population_weight_within_country"] = (
        region.loc[valid_pop, "ghs_pop_2025"]
        / region.loc[valid_pop].groupby("sci_country_iso2")["ghs_pop_2025"].transform("sum")
    )
    country_mean = (
        region.loc[valid_pop]
        .assign(
            weighted=lambda x: x["regional_sci_geography_adjusted_2026"]
            * x["population_weight_within_country"]
        )
        .groupby("sci_country_iso2")["weighted"]
        .sum()
    )
    region["country_average_regional_sci_geography_adjusted_2026"] = region[
        "sci_country_iso2"
    ].map(country_mean)
    region["regional_sci_geography_adjusted_deviation_2026"] = (
        region["regional_sci_geography_adjusted_2026"]
        - region["country_average_regional_sci_geography_adjusted_2026"]
    )
    region.to_csv(DERIVED / "regional_sci_geography_adjusted_adm1_2026.csv", index=False)
    coefficients["model"] = (
        "ADM1 log SCI on geography/history, origin-country FE, and destination-country FE"
    )
    for key, value in diagnostics.items():
        coefficients[key] = value
    coefficients["match_rate"] = match_rate
    coefficients.to_csv(DERIVED / "regional_sci_geography_adjustment_model.csv", index=False)
    return region, coefficients


def main() -> None:
    ensure_directories()
    country, _ = build_country_adjusted()
    region, _ = build_regional_adjusted()
    if "NA" not in set(country["iso2"].astype(str)):
        raise RuntimeError("Namibia was lost from geography-adjusted country SCI")
    print(
        f"Built geography-adjusted SCI for {country['iso2'].nunique():,} country origins "
        f"and {region['adm1_gid'].nunique():,} ADM1 regions."
    )


if __name__ == "__main__":
    main()
