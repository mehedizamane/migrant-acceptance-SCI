#!/usr/bin/env python3
"""Link Gallup respondents to GADM 4.1 ADM1 regions and build the regional model input.

For each country-wave, the base REGION_* variable with the most coded respondents
is selected before outcomes are inspected. A region code is linked only when the
candidate inventory from 08_match_gallup_regions.py assigns it to exactly one
GADM1 unit by a deterministic rule (official code, primary or official alternate
name, administrative-type normalization, verified prior exact match, ISO 3166-2
subdivision name, or an exact lower-level unit nested in one GADM1 parent).
Broader regions, fuzzy candidates, and countries without compatible subnational
SCI remain unlinked; no respondent receives a country-average or nearby exposure.

Outputs:
  outputs/regional_linkage/  aggregate inputs to SI Tables 1 and 15
  data/restricted/derived/geography_direct_adm1/  restricted model input and crosswalk
"""

from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import DERIVED, OUTPUTS, RESTRICTED, YEARS, weighted_mean, weighted_sd


BASE = RESTRICTED / "derived" / "gallup_four_wave_sci_analytic.parquet"
ARCHIVE = RESTRICTED / "raw" / "Gallup_World_Poll_022026_ALL_WAVES.zip"
ARCHIVE_MEMBER = "The_Gallup_022026.dat"
CANDIDATES = RESTRICTED / "derived" / "geography_v2" / "all_observed_region_code_candidates.csv"
CANONICAL = DERIVED / "regional_sci_adm1_2026.csv"
VARIANCE = DERIVED / "regional_sci_variance_decomposition.csv"
REGIONAL_GEO_SCI = DERIVED / "regional_sci_geography_adjusted_adm1_2026.csv"
OUT = OUTPUTS / "regional_linkage"
RESTRICTED_OUT = RESTRICTED / "derived" / "geography_direct_adm1"

BASE_REGION_RE = re.compile(r"^REGION_[A-Z0-9]+$")
ACCEPTED_METHODS = {
    "exact_official_adm1_code",
    "exact_gadm1_primary_name",
    "exact_gadm1_official_alternate_name",
    "exact_normalized_administrative_core_candidate",
    "verified_version1_exact_normalized_gadm1_name",
    "exact_iso_3166_2_subdivision_name_candidate",
    "exact_lower_gadm_unit_to_unique_adm1_parent",
}
# SI Table 14 drops administrative-core normalization, ISO subdivision names,
# and lower-level nesting.
EXACT_ONLY_METHODS = {
    "exact_official_adm1_code",
    "exact_gadm1_primary_name",
    "exact_gadm1_official_alternate_name",
    "verified_version1_exact_normalized_gadm1_name",
}
INDIVIDUAL = [
    "age", "female", "education_level", "income_quintile", "foreign_born", "urbanicity",
]
EXPECTED_TOTAL = 3_028_581
EXPECTED_WAVE_ROWS = {2016: 150_724, 2019: 176_253, 2022: 143_686, 2023: 145_702}


def canonical_code(series: pd.Series) -> pd.Series:
    numeric = pd.to_numeric(series.replace("", np.nan), errors="coerce")
    return numeric.astype("Int64").astype("string").fillna("")


def rejection_reason(row: pd.Series) -> str:
    method = row["match_method"]
    count = row["constituent_adm1_count_num"]
    if method == "no_compatible_subnational_sci":
        return "Country absent from the compatible Meta GADM1 SCI origin universe"
    if method == "unmatched":
        return "No direct GADM1 name or code match"
    if count > 1:
        return "Gallup region aggregates multiple GADM1 units"
    if "fuzzy" in method:
        return "Fuzzy candidate withheld pending documented manual review"
    if "aggregate" in method or "historical" in method:
        return "Aggregate or historical concordance requires authoritative composition"
    if method.startswith("version1_"):
        return "Legacy mapping lacks sufficient source evidence for automatic reuse"
    return "Mapping method is outside the prespecified direct-match hierarchy"


def direct_class(method: str) -> str:
    if method in {"exact_official_adm1_code", "exact_gadm1_primary_name"}:
        return "Exact official code or primary GADM1 name"
    if method in {
        "exact_gadm1_official_alternate_name",
        "verified_version1_exact_normalized_gadm1_name",
    }:
        return "Exact official alternate or previously verified normalized name"
    if method == "exact_normalized_administrative_core_candidate":
        return "Unique exact name after removing administrative-type words"
    if method == "exact_iso_3166_2_subdivision_name_candidate":
        return "Unique exact ISO 3166-2 subdivision name"
    if method == "exact_lower_gadm_unit_to_unique_adm1_parent":
        return "Exact lower-level GADM unit nested in one GADM1 parent"
    return "Other"


def valid_rows(frame: pd.DataFrame, columns: list[str]) -> pd.Series:
    return frame[columns].notna().all(axis=1)


def stream_selected_base_regions(selection: pd.DataFrame, chunksize: int = 25_000) -> pd.DataFrame:
    """Read one prespecified base REGION_* field per country-wave from Gallup."""
    selected_map = selection.set_index(["iso2", "year_wave"])["gallup_region_variable"].to_dict()
    variables = sorted(selection["gallup_region_variable"].unique())
    process = subprocess.Popen(
        ["unzip", "-p", str(ARCHIVE), ARCHIVE_MEMBER], stdout=subprocess.PIPE
    )
    if process.stdout is None:
        raise RuntimeError("Could not stream the Gallup archive")
    pieces: list[pd.DataFrame] = []
    total = 0
    wave_counts = {year: 0 for year in YEARS}
    try:
        reader = pd.read_csv(
            process.stdout, sep="\t", usecols=["COUNTRY_ISO2", "YEAR_WAVE", *variables],
            dtype="string", keep_default_na=False, na_values=[], chunksize=chunksize,
            low_memory=False,
        )
        for chunk in reader:
            total += len(chunk)
            years = pd.to_numeric(chunk["YEAR_WAVE"].replace("", np.nan), errors="coerce")
            keep = years.isin(YEARS)
            if not keep.any():
                continue
            use = chunk.loc[keep]
            year_values = years.loc[keep].astype(int)
            iso2 = use["COUNTRY_ISO2"].astype("string").str.strip()
            for year, count in year_values.value_counts().items():
                wave_counts[int(year)] += int(count)
            linked = pd.DataFrame({
                "iso2_stream": iso2,
                "year_wave_stream": year_values,
                "gallup_region_variable_direct": "",
                "gallup_region_code_direct": "",
            }, index=use.index)
            for (country, year), variable in selected_map.items():
                mask = iso2.eq(country) & year_values.eq(int(year))
                if not mask.any():
                    continue
                linked.loc[mask, "gallup_region_variable_direct"] = variable
                linked.loc[mask, "gallup_region_code_direct"] = canonical_code(use.loc[mask, variable])
            pieces.append(linked.reset_index(drop=True))
    finally:
        process.stdout.close()
        process.wait()
    if process.returncode != 0:
        raise RuntimeError(f"Gallup unzip stream exited with status {process.returncode}")
    if total != EXPECTED_TOTAL or wave_counts != EXPECTED_WAVE_ROWS:
        raise RuntimeError(f"Gallup row counts changed: total={total}; waves={wave_counts}")
    return pd.concat(pieces, ignore_index=True)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    RESTRICTED_OUT.mkdir(parents=True, exist_ok=True)
    if not CANDIDATES.exists():
        raise FileNotFoundError("Run 08_match_gallup_regions.py to generate the candidate inventory")

    candidates = pd.read_csv(CANDIDATES, keep_default_na=False, dtype=str)
    candidates = candidates.loc[
        candidates["gallup_region_variable"].map(lambda value: bool(BASE_REGION_RE.fullmatch(value)))
        & pd.to_numeric(candidates["year_wave"], errors="coerce").isin(YEARS)
    ].copy()
    candidates["year_wave"] = pd.to_numeric(candidates["year_wave"], errors="raise").astype(int)
    candidates["gallup_region_code"] = canonical_code(candidates["gallup_region_code"])
    candidates["respondents"] = pd.to_numeric(candidates["respondents"], errors="coerce").fillna(0).astype(int)
    candidates["constituent_adm1_count_num"] = pd.to_numeric(
        candidates["constituent_adm1_count"], errors="coerce"
    ).fillna(0).astype(int)
    key = ["iso2", "year_wave", "gallup_region_variable", "gallup_region_code"]
    if candidates.duplicated(key).any():
        raise RuntimeError("The base REGION_* candidate inventory is not unique by country-wave-variable-code")

    # One base REGION_* variable per country-wave: the one with the most coded
    # respondents, chosen before outcomes are inspected.
    variable_summary = candidates.groupby(
        ["iso2", "year_wave", "gallup_region_variable"], as_index=False
    ).agg(
        observed_codes=("gallup_region_code", "nunique"),
        respondents_with_code=("respondents", "sum"),
    )
    variable_summary = variable_summary.sort_values(
        ["iso2", "year_wave", "respondents_with_code", "gallup_region_variable"],
        ascending=[True, True, False, True],
    )
    selection = variable_summary.drop_duplicates(["iso2", "year_wave"], keep="first").copy()
    candidates = candidates.merge(
        selection[["iso2", "year_wave", "gallup_region_variable"]],
        on=["iso2", "year_wave", "gallup_region_variable"], how="inner", validate="many_to_one",
    )

    accepted = candidates.loc[
        candidates["match_method"].isin(ACCEPTED_METHODS)
        & candidates["constituent_adm1_count_num"].eq(1)
        & ~candidates["constituent_adm1_ids"].str.contains(r"[|;,]", regex=True)
        & candidates["constituent_adm1_ids"].ne("")
    ].copy()
    accepted["adm1_gid_direct"] = accepted["constituent_adm1_ids"]
    accepted["direct_match_class"] = accepted["match_method"].map(direct_class)
    accepted["exact_only_sample"] = accepted["match_method"].isin(EXACT_ONLY_METHODS)
    accepted["mapping_id_direct"] = (
        "direct|" + accepted["iso2"] + "|" + accepted["year_wave"].astype(str)
        + "|" + accepted["gallup_region_variable"] + "|" + accepted["gallup_region_code"]
    )

    canonical = pd.read_csv(CANONICAL, keep_default_na=False, dtype={"analysis_iso2": str})
    numeric_canonical = [
        "regional_sci_2026", "regional_sci_foreign_mean_raw_2026", "ghs_pop_2025",
        "country_average_regional_sci_2026_population_weighted",
        "regional_sci_deviation_from_country_2026", "adm1_log_population",
        "adm1_log_population_density", "population_weight_within_country",
    ]
    for column in numeric_canonical:
        canonical[column] = pd.to_numeric(canonical[column], errors="coerce")
    canonical_keep = [
        "adm1_gid", "analysis_iso2", "adm1_name", "regional_sci_2026",
        "regional_sci_foreign_mean_raw_2026", "ghs_pop_2025",
        "country_average_regional_sci_2026_population_weighted",
        "regional_sci_deviation_from_country_2026", "adm1_log_population",
        "adm1_log_population_density", "population_weight_within_country",
    ]
    accepted = accepted.merge(
        canonical[canonical_keep], left_on="adm1_gid_direct", right_on="adm1_gid",
        how="left", validate="many_to_one",
    )
    missing_meta = accepted["regional_sci_2026"].isna()
    if missing_meta.any():
        raise RuntimeError(f"{missing_meta.sum()} accepted direct mappings lack canonical 2026 regional SCI")
    country_mismatch = accepted["iso2"].ne(accepted["analysis_iso2"])
    if country_mismatch.any():
        examples = accepted.loc[country_mismatch, ["iso2", "analysis_iso2", "adm1_gid_direct"]].head()
        raise RuntimeError(f"Direct mappings cross national boundaries:\n{examples.to_string(index=False)}")
    accepted[[
        "mapping_id_direct", *key, "gallup_region_label", "respondents", "adm1_gid_direct",
        "adm1_name", "direct_match_class", "match_method", "exact_only_sample",
    ]].to_csv(RESTRICTED_OUT / "direct_adm1_crosswalk.csv", index=False)

    # Region codes that were not linked, with the reason (SI Table 1, Panel C).
    accepted_keys = set(map(tuple, accepted[key].astype(str).to_numpy()))
    rejected = candidates.loc[
        ~candidates[key].astype(str).apply(tuple, axis=1).isin(accepted_keys)
    ].copy()
    rejected["rejection_reason"] = rejected.apply(rejection_reason, axis=1)
    rejected.groupby("rejection_reason", as_index=False).agg(
        country_wave_codes=("rejection_reason", "size"),
        respondents=("respondents", "sum"),
        countries=("iso2", "nunique"),
    ).sort_values("respondents", ascending=False).to_csv(OUT / "exclusion_reasons.csv", index=False)

    # Accepted match classes (SI Table 1, Panel B).
    accepted.groupby(["direct_match_class", "match_method"], as_index=False).agg(
        mapped_country_wave_codes=("mapping_id_direct", "size"),
        candidate_respondents=("respondents", "sum"),
        countries=("iso2", "nunique"),
    ).sort_values("candidate_respondents", ascending=False).to_csv(OUT / "match_classes.csv", index=False)

    base = pd.read_parquet(BASE).copy().reset_index(drop=True)
    streamed = stream_selected_base_regions(selection)
    if len(base) != len(streamed):
        raise RuntimeError(
            f"Frozen analytic file has {len(base):,} rows but raw geography stream has {len(streamed):,}"
        )
    base_year = pd.to_numeric(base["year_wave"], errors="coerce").astype("Int64")
    stream_year = pd.to_numeric(streamed["year_wave_stream"], errors="coerce").astype("Int64")
    same_keys = (
        base["iso2"].astype(str).eq(streamed["iso2_stream"].astype(str))
        & base_year.eq(stream_year)
    )
    if not same_keys.all():
        first = int(np.flatnonzero(~same_keys.to_numpy())[0])
        raise RuntimeError(f"Frozen analytic and raw Gallup streams diverge at row {first}")
    base["gallup_region_variable"] = streamed["gallup_region_variable_direct"].to_numpy()
    base["gallup_region_code"] = streamed["gallup_region_code_direct"].to_numpy()
    base["year_wave"] = pd.to_numeric(base["year_wave"], errors="raise").astype(int)
    base["gallup_region_variable"] = base["gallup_region_variable"].fillna("").astype(str)
    base["gallup_region_code"] = canonical_code(base["gallup_region_code"].fillna(""))
    link_columns = [
        *key, "mapping_id_direct", "adm1_gid_direct", "direct_match_class", "match_method",
        "exact_only_sample", "regional_sci_2026",
        "regional_sci_foreign_mean_raw_2026", "ghs_pop_2025",
        "country_average_regional_sci_2026_population_weighted",
        "regional_sci_deviation_from_country_2026", "adm1_log_population",
        "adm1_log_population_density",
    ]
    linked = base.merge(accepted[link_columns], on=key, how="left", validate="many_to_one", suffixes=("", "_direct"))
    linked["direct_adm1_match"] = linked["mapping_id_direct"].notna()
    linked["direct_adm1_gid"] = linked["adm1_gid_direct"].fillna("")
    linked["direct_exact_only"] = linked["exact_only_sample"].eq(True)
    if linked.loc[linked["direct_adm1_match"], "regional_sci_2026_direct"].isna().any():
        raise RuntimeError("A linked respondent lacks the direct regional SCI exposure")
    if linked.loc[linked["iso2"].eq("NA")].empty:
        raise RuntimeError("Namibia's NA identifier was lost")

    candidate_keys = candidates[key + ["match_method"]].rename(
        columns={"match_method": "candidate_match_method"}
    )
    linked = linked.merge(candidate_keys, on=key, how="left", validate="many_to_one")
    linked["usable_base_region_code"] = linked["gallup_region_variable"].ne("") & linked["gallup_region_code"].ne("")
    linked["compatible_subnational_sci"] = (
        linked["candidate_match_method"].notna()
        & linked["candidate_match_method"].ne("no_compatible_subnational_sci")
    )

    # Linked and unlinked respondents, unadjusted weighted mean MAI (SI Table 15).
    selection_rows = []
    for year, wave in linked.groupby("year_wave"):
        for matched, label in [(True, "Directly linked"), (False, "Not directly linked")]:
            use = wave[wave["direct_adm1_match"].eq(matched)].copy()
            use["weight_raw"] = pd.to_numeric(use["weight_raw"], errors="coerce")
            use["migrant_acceptance_index"] = pd.to_numeric(use["migrant_acceptance_index"], errors="coerce")
            use = use[use["weight_raw"].gt(0) & use["migrant_acceptance_index"].notna()]
            mean = float(np.average(use["migrant_acceptance_index"], weights=use["weight_raw"]))
            selection_rows.append({
                "year_wave": int(year), "linkage": label, "respondents": len(use),
                "countries": use["iso2"].nunique(), "weighted_mean_mai": mean,
            })
    pd.DataFrame(selection_rows).to_csv(OUT / "linked_vs_unlinked.csv", index=False)

    # Express the country-average and within-country components in the same
    # total RegionalSCI SD units (0.773 log-SCI units).
    variance = pd.read_csv(VARIANCE, keep_default_na=False)
    total_sd = float(
        variance.loc[variance["component"].eq("Total"), "standard_deviation"].iloc[0]
    )
    reference = canonical.loc[
        canonical["population_weight_within_country"].gt(0)
        & canonical["regional_sci_2026"].notna()
    ].copy()
    country_count = reference["analysis_iso2"].nunique()
    reference_weight = reference["population_weight_within_country"] / country_count
    grand_mean = float(np.average(reference["regional_sci_2026"], weights=reference_weight))

    linked["country_id"] = linked["iso2"].astype(str)
    linked["exposure_region_id"] = linked["direct_adm1_gid"]
    linked["age_10"] = (pd.to_numeric(linked["age"], errors="coerce") - 40) / 10
    linked["survey_region_sci_total_z"] = (linked["regional_sci_2026_direct"] - grand_mean) / total_sd
    linked["survey_region_sci_between_total_sd"] = (
        linked["country_average_regional_sci_2026_population_weighted_direct"] - grand_mean
    ) / total_sd
    linked["survey_region_sci_within_total_sd"] = (
        linked["regional_sci_deviation_from_country_2026_direct"] / total_sd
    )
    reconstruction = linked["survey_region_sci_between_total_sd"] + linked["survey_region_sci_within_total_sd"]
    check = linked["direct_adm1_match"]
    if not np.allclose(reconstruction.loc[check], linked.loc[check, "survey_region_sci_total_z"], atol=1e-12):
        raise RuntimeError("Direct ADM1 between/within components do not reconstruct total regional SCI")

    linked_rows = linked.loc[
        linked["direct_adm1_match"] & linked["weight_equal_country_equal_wave"].gt(0)
    ].copy()
    for source, target in [
        ("adm1_log_population_direct", "survey_region_log_population_z"),
        ("adm1_log_population_density_direct", "survey_region_log_density_z"),
    ]:
        mean = weighted_mean(linked_rows[source], linked_rows["weight_equal_country_equal_wave"])
        sd = weighted_sd(linked_rows[source], linked_rows["weight_equal_country_equal_wave"])
        if not np.isfinite(sd) or sd <= 0:
            raise RuntimeError(f"Cannot standardize {source}")
        linked[target] = (linked[source] - mean) / sd

    # Geography-adjusted regional deviation (SI Table 12, bottom panel), in units of
    # its population-weighted within-country SD over all supported ADM1 regions.
    geo = pd.read_csv(REGIONAL_GEO_SCI, keep_default_na=False)
    for column in ["regional_sci_geography_adjusted_deviation_2026", "population_weight_within_country"]:
        geo[column] = pd.to_numeric(geo[column], errors="coerce")
    geo_reference = geo.loc[
        geo["population_weight_within_country"].gt(0)
        & geo["regional_sci_geography_adjusted_deviation_2026"].notna()
    ]
    geo_weight = geo_reference["population_weight_within_country"] / geo_reference["sci_country_iso2"].nunique()
    geo_sd = float(np.sqrt(np.average(
        np.square(geo_reference["regional_sci_geography_adjusted_deviation_2026"]), weights=geo_weight
    )))
    geo_z = (geo.set_index("adm1_gid")["regional_sci_geography_adjusted_deviation_2026"] / geo_sd)
    linked["regional_sci_geography_adjusted_within_z"] = linked["direct_adm1_gid"].map(geo_z)

    required_primary = [
        "migrant_acceptance_index", "weight_raw", "country_id", "exposure_region_id",
        "survey_region_sci_total_z", "survey_region_log_population_z",
        "survey_region_log_density_z", *INDIVIDUAL,
    ]
    linked["regional_primary_complete"] = (
        linked["direct_adm1_match"] & linked["weight_raw"].gt(0)
        & valid_rows(linked, required_primary)
    )
    linked["exact_only_sample"] = linked["direct_exact_only"]

    keep = [
        "country_id", "exposure_region_id", "year_wave", "weight_raw",
        "migrant_acceptance_index", "exact_only_sample",
        "survey_region_sci_between_total_sd", "survey_region_sci_within_total_sd",
        "regional_sci_geography_adjusted_within_z",
        "age_10", "female", "education_level", "income_quintile", "foreign_born", "urbanicity",
        "survey_region_log_population_z", "survey_region_log_density_z",
    ]
    model_input = linked.loc[linked["regional_primary_complete"], keep].copy()
    model_input_path = RESTRICTED_OUT / "gallup_regional_model_input_direct_adm1.csv"
    model_input.to_csv(model_input_path, index=False, na_rep="NA_REAL")

    # Linkage coverage and model sample by wave (SI Table 1, Panel A).
    coverage_rows = []
    for year in YEARS:
        wave = linked.loc[linked["year_wave"].eq(year)]
        direct = wave["direct_adm1_match"]
        eligible = wave["usable_base_region_code"] & wave["compatible_subnational_sci"]
        model = wave.loc[wave["regional_primary_complete"]]
        coverage_rows.append({
            "year_wave": year,
            "all_wave_respondents": len(wave),
            "respondents_with_base_region_code": int(wave["usable_base_region_code"].sum()),
            "eligible_geography_respondents": int(eligible.sum()),
            "linked_respondents": int(direct.sum()),
            "linked_pct_all": float(100 * direct.mean()),
            "linked_pct_eligible": float(100 * direct.sum() / max(1, eligible.sum())),
            "model_respondents": len(model),
            "model_countries": model["iso2"].nunique(),
            "model_adm1_regions": model["direct_adm1_gid"].nunique(),
        })
    coverage = pd.DataFrame(coverage_rows)
    coverage.to_csv(OUT / "linkage_coverage.csv", index=False)
    print(coverage.to_string(index=False))
    print(f"Wrote {len(model_input):,} regional model-input rows to {model_input_path}")


if __name__ == "__main__":
    main()
