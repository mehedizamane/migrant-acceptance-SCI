#!/usr/bin/env python3
"""Build the restricted four-wave Gallup analytic file.

Streams the licensed Gallup archive with ``unzip -p`` so the 17.5 GB text
member is never extracted to disk.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import DERIVED, RAW, RESTRICTED, YEARS, ensure_directories, require_columns


GALLUP_ARCHIVE = RESTRICTED / "raw" / "Gallup_World_Poll_022026_ALL_WAVES.zip"
GALLUP_MEMBER = "The_Gallup_022026.dat"
CODEBOOK = RESTRICTED / "raw" / "GWP_022026_Codebook_dta.txt"
CROSSWALK = RAW / "geography" / "v1" / "gallup_to_gadm1_crosswalk_v1.csv"
COUNTRY_SCI = DERIVED / "country_global_sci_2020_2026.csv"
REGIONAL_SCI = DERIVED / "regional_sci_adm1_2026.csv"
COUNTRY_CONTROLS = DERIVED / "country_controls.csv"
EMPLOYMENT_CONTROLS = DERIVED / "country_wave_employment_controls.csv"
OUTPUT = RESTRICTED / "derived" / "gallup_four_wave_sci_analytic.parquet"

EXPECTED_TOTAL = 3_028_581
EXPECTED_WAVE_ROWS = {2016: 150_724, 2019: 176_253, 2022: 143_686, 2023: 145_702}

BASE_COLUMNS = [
    "COUNTRY_ISO2",
    "COUNTRY_ISO3",
    "COUNTRYNEW",
    "WGT",
    "YEAR_WAVE",
    "INDEX_MAI",
    "WP17639",
    "WP17640",
    "WP17641",
    "SC_22F",
    "WP1219",
    "WP1220",
    "WP3117",
    "INCOME_5",
    "EMP_2010",
    "DEGURBA",
    "WP14",
    "WP7572",
    "WP4657",
]

CODEBOOK_LABELS = {
    "INDEX_MAI": "Migrant Acceptance Index",
    "WP17639": "Immigrants Living in (Country)",
    "WP17640": "An Immigrant Becoming Your Neighbor",
    "WP17641": "An Immigrant Marrying One of Your Close Relatives",
    "SC_22F": "How often did you interact with strangers",
    "WP1219": "Gender",
    "WP1220": "Age",
    "WP3117": "Education Level",
    "INCOME_5": "Per Capita Income Quintiles",
    "EMP_2010": "Employment Status",
    "WP14": "Urban/Rural",
    "WP7572": "Urban/Rural",
    "WP4657": "Born in Country",
}


def numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series.replace("", np.nan), errors="coerce")


def recode_binary(series: pd.Series, true_code: int, false_code: int) -> pd.Series:
    values = numeric(series)
    result = pd.Series(np.nan, index=series.index, dtype=float)
    result.loc[values.eq(true_code)] = 1.0
    result.loc[values.eq(false_code)] = 0.0
    return result


def recode_component(series: pd.Series) -> pd.Series:
    # Gallup coding: 1 = good thing, 2 = bad thing; depends, DK, and refused
    # are excluded from the binary component outcome.
    return recode_binary(series, true_code=1, false_code=2)


def recode_urbanicity(frame: pd.DataFrame) -> pd.Series:
    degurba = numeric(frame["DEGURBA"])
    legacy = numeric(frame["WP14"])
    older = numeric(frame["WP7572"])
    result = pd.Series(pd.NA, index=frame.index, dtype="string")
    result.loc[degurba.eq(1)] = "rural"
    result.loc[degurba.eq(2)] = "town_or_semi_dense"
    result.loc[degurba.eq(3)] = "city"
    missing = result.isna()
    result.loc[missing & legacy.eq(1)] = "rural"
    result.loc[missing & legacy.eq(2)] = "town_or_semi_dense"
    result.loc[missing & legacy.eq(3)] = "city"
    result.loc[missing & legacy.eq(6)] = "city"
    missing = result.isna()
    result.loc[missing & older.eq(1)] = "city"
    result.loc[missing & older.eq(2)] = "rural"
    return result


def validate_codebook() -> None:
    text = CODEBOOK.read_text(encoding="utf-8", errors="replace")
    missing = [code for code, label in CODEBOOK_LABELS.items() if code not in text or label not in text]
    if "DEGURBA" not in text:
        missing.append("DEGURBA")
    if missing:
        raise RuntimeError(f"Gallup codebook labels or codes changed: {sorted(set(missing))}")


def load_crosswalk_maps() -> tuple[list[str], dict[tuple[str, str, str], str], dict[tuple[str, str], int], dict[str, list[str]]]:
    crosswalk = pd.read_csv(CROSSWALK, keep_default_na=False, dtype=str)
    required = ["iso2", "gallup_region_variable", "gallup_region_code", "adm1_gid", "match_status"]
    require_columns(crosswalk, required, "Gallup-GADM1 crosswalk")
    accepted = crosswalk.loc[crosswalk["match_status"].eq("accepted")].copy()
    accepted = accepted.loc[accepted["adm1_gid"].str.strip().ne("")]
    accepted["gallup_region_code"] = (
        pd.to_numeric(accepted["gallup_region_code"], errors="coerce")
        .astype("Int64")
        .astype(str)
    )
    conflicts = accepted.groupby(
        ["iso2", "gallup_region_variable", "gallup_region_code"]
    )["adm1_gid"].nunique()
    if conflicts.gt(1).any():
        raise RuntimeError("A Gallup country/region code maps to multiple accepted ADM1 regions")
    map_gid = (
        accepted.drop_duplicates(["iso2", "gallup_region_variable", "gallup_region_code"])
        .set_index(["iso2", "gallup_region_variable", "gallup_region_code"])["adm1_gid"]
        .to_dict()
    )
    priority_table = (
        accepted.groupby(["iso2", "gallup_region_variable"])["gallup_region_code"]
        .nunique()
        .rename("accepted_codes")
        .reset_index()
        .sort_values(["iso2", "accepted_codes", "gallup_region_variable"], ascending=[True, False, True])
    )
    priority: dict[tuple[str, str], int] = {}
    by_country: dict[str, list[str]] = {}
    for iso2, part in priority_table.groupby("iso2", sort=False):
        variables = part["gallup_region_variable"].astype(str).tolist()
        by_country[str(iso2)] = variables
        for rank, variable in enumerate(variables):
            priority[(str(iso2), variable)] = rank
    variables = sorted(accepted["gallup_region_variable"].astype(str).unique())
    return variables, map_gid, priority, by_country


def harmonize_chunk(
    use: pd.DataFrame,
    region_map: dict[tuple[str, str, str], str],
    priority: dict[tuple[str, str], int],
    by_country: dict[str, list[str]],
) -> pd.DataFrame:
    iso2 = use["COUNTRY_ISO2"].astype("string").str.strip()
    selected_variable = pd.Series("", index=use.index, dtype="string")
    selected_code = pd.Series("", index=use.index, dtype="string")
    selected_gid = pd.Series("", index=use.index, dtype="string")
    selected_priority = pd.Series(np.inf, index=use.index, dtype=float)
    for country_code in iso2.unique():
        variables = by_country.get(str(country_code), [])
        if not variables:
            continue
        country_mask = iso2.eq(country_code)
        for variable in variables:
            codes = numeric(use.loc[country_mask, variable]).astype("Int64").astype(str)
            gids = pd.Series(
                [region_map.get((str(country_code), variable, code), "") for code in codes],
                index=codes.index,
                dtype="string",
            )
            rank = priority[(str(country_code), variable)]
            mask = gids.ne("") & selected_priority.loc[codes.index].gt(rank)
            target_index = codes.index[mask]
            selected_variable.loc[target_index] = variable
            selected_code.loc[target_index] = codes.loc[target_index]
            selected_gid.loc[target_index] = gids.loc[target_index]
            selected_priority.loc[target_index] = rank

    living = recode_component(use["WP17639"])
    neighbor = recode_component(use["WP17640"])
    relative = recode_component(use["WP17641"])
    strict_components = pd.concat([living, neighbor, relative], axis=1)
    strict_mai = strict_components.sum(axis=1, min_count=3) * 3
    year = numeric(use["YEAR_WAVE"]).astype("Int64")
    result = pd.DataFrame(
        {
            "iso2": iso2,
            "iso3": use["COUNTRY_ISO3"].astype("string").str.strip(),
            "country_name_gallup": use["COUNTRYNEW"].astype("string").str.strip(),
            "year_wave": year,
            "weight_raw": numeric(use["WGT"]),
            "gallup_region_variable": selected_variable,
            "gallup_region_code": selected_code,
            "adm1_gid": selected_gid,
            "migrant_acceptance_index": numeric(use["INDEX_MAI"]).where(
                numeric(use["INDEX_MAI"]).between(0, 9)
            ),
            "migrant_acceptance_index_strict": strict_mai,
            "immigrants_living_in_country_good": living,
            "immigrant_neighbor_good": neighbor,
            "immigrant_marry_close_relative_good": relative,
            "stranger_contact_frequency": numeric(use["SC_22F"]).where(
                numeric(use["SC_22F"]).isin([1, 2, 3, 4, 5])
            ),
            "female": recode_binary(use["WP1219"], true_code=2, false_code=1),
            "age": numeric(use["WP1220"]).where(numeric(use["WP1220"]).between(15, 100)),
            "education_level": numeric(use["WP3117"]).where(
                numeric(use["WP3117"]).isin([1, 2, 3])
            ),
            "income_quintile": numeric(use["INCOME_5"]).where(
                numeric(use["INCOME_5"]).isin([1, 2, 3, 4, 5])
            ),
            "employment_status": numeric(use["EMP_2010"]).where(
                numeric(use["EMP_2010"]).isin([1, 2, 3, 4, 5, 6])
            ),
            "urbanicity": recode_urbanicity(use),
            # WP4657 asks whether the respondent was born in the survey
            # country: 1=yes and 2=no. We code 1 for foreign born.
            "foreign_born": recode_binary(use["WP4657"], true_code=2, false_code=1),
        }
    )
    return result


def stream_raw(chunksize: int = 25_000) -> tuple[pd.DataFrame, dict[str, object]]:
    validate_codebook()
    region_variables, region_map, priority, by_country = load_crosswalk_maps()
    usecols = list(dict.fromkeys(BASE_COLUMNS + region_variables))
    process = subprocess.Popen(
        ["unzip", "-p", str(GALLUP_ARCHIVE), GALLUP_MEMBER], stdout=subprocess.PIPE
    )
    if process.stdout is None:
        raise RuntimeError("Could not stream the Gallup archive")
    total = 0
    wave_counts = {year: 0 for year in YEARS}
    pieces: list[pd.DataFrame] = []
    try:
        reader = pd.read_csv(
            process.stdout,
            sep="\t",
            usecols=usecols,
            dtype="string",
            keep_default_na=False,
            na_values=[],
            chunksize=chunksize,
            low_memory=False,
        )
        for chunk in reader:
            total += len(chunk)
            years = numeric(chunk["YEAR_WAVE"])
            keep = years.isin(YEARS)
            if not keep.any():
                continue
            use = chunk.loc[keep].copy()
            for year, count in years.loc[keep].value_counts().items():
                wave_counts[int(year)] += int(count)
            pieces.append(harmonize_chunk(use, region_map, priority, by_country))
    finally:
        process.stdout.close()
        process.wait()
    if process.returncode != 0:
        raise RuntimeError(f"Gallup unzip stream exited with status {process.returncode}")
    if total != EXPECTED_TOTAL or wave_counts != EXPECTED_WAVE_ROWS:
        raise RuntimeError(f"Gallup row counts changed: total={total}; waves={wave_counts}")
    return pd.concat(pieces, ignore_index=True), {
        "source_mode": "raw_stream",
        "streamed_total_rows": total,
        "wave_rows": wave_counts,
    }


def merge_exposures_and_controls(frame: pd.DataFrame) -> pd.DataFrame:
    country_sci = pd.read_csv(COUNTRY_SCI, keep_default_na=False)
    country_controls = pd.read_csv(COUNTRY_CONTROLS, keep_default_na=False)
    regional = pd.read_csv(REGIONAL_SCI, keep_default_na=False)
    employment = pd.read_csv(EMPLOYMENT_CONTROLS, keep_default_na=False)
    for column in [
        "global_sci_2020",
        "global_sci_2020_z",
        "global_sci_2020_all_origins_z",
        "global_sci_2026",
        "global_sci_2026_z",
        "global_sci_2026_all_origins_z",
    ]:
        country_sci[column] = pd.to_numeric(country_sci[column], errors="coerce")
    country_columns = [
        "iso2",
        "global_sci_2020",
        "global_sci_2020_z",
        "global_sci_2020_all_origins_z",
        "global_sci_2026",
        "global_sci_2026_z",
        "global_sci_2026_all_origins_z",
        "is_common_release_origin",
    ]
    frame = frame.merge(country_sci[country_columns], on="iso2", how="left", validate="many_to_one")
    control_columns = [
        "iso2",
        "country_name",
        "region",
        "un_m49_region",
        "un_m49_subregion",
        "un_region_9",
        "log_gdp_per_capita",
        "log_population",
        "internet_users_pct",
        "mobile_cellular_subscriptions_per_100",
        "fixed_broadband_subscriptions_per_100",
        "log_secure_internet_servers_per_million_people",
        "urban_population_pct",
        "trade_pct_gdp",
        "tourism_arrivals_per_capita",
        "migrant_stock_share",
        "refugee_stock_share",
    ]
    for column in control_columns:
        if column not in {
            "iso2",
            "country_name",
            "region",
            "un_m49_region",
            "un_m49_subregion",
            "un_region_9",
        }:
            country_controls[column] = pd.to_numeric(
                country_controls[column], errors="coerce"
            )
    frame = frame.merge(
        country_controls[control_columns], on="iso2", how="left", validate="many_to_one"
    )
    employment_columns = ["iso2", "year_wave", "employment_to_population_ratio_15plus"]
    employment["year_wave"] = pd.to_numeric(employment["year_wave"], errors="coerce").astype("Int64")
    employment["employment_to_population_ratio_15plus"] = pd.to_numeric(
        employment["employment_to_population_ratio_15plus"], errors="coerce"
    )
    employment = employment[employment_columns].drop_duplicates(["iso2", "year_wave"])
    frame = frame.merge(
        employment,
        on=["iso2", "year_wave"],
        how="left",
        validate="many_to_one",
    )
    regional_columns = [
        "adm1_gid",
        "analysis_iso2",
        "sci_country_iso2",
        "regional_sci_2026",
        "country_average_regional_sci_2026_population_weighted",
        "regional_sci_deviation_from_country_2026",
        "regional_sci_deviation_within_country_z",
        "regional_sci_deviation_total_sd_units",
        "ghs_pop_2025",
        "adm1_area_km2",
        "adm1_log_population",
        "adm1_log_population_density",
        "population_valid",
    ]
    for column in regional_columns:
        if column not in {"adm1_gid", "analysis_iso2", "sci_country_iso2"}:
            regional[column] = pd.to_numeric(regional[column], errors="coerce")
    frame = frame.merge(regional[regional_columns], on="adm1_gid", how="left", validate="many_to_one")
    frame["regional_country_match"] = (
        frame["analysis_iso2"].astype(str).eq(frame["iso2"].astype(str))
        | frame["analysis_iso2"].isna()
        | frame["analysis_iso2"].eq("")
    )
    invalid = frame["adm1_gid"].astype(str).str.strip().ne("") & ~frame["regional_country_match"]
    if invalid.any():
        examples = frame.loc[invalid, ["iso2", "analysis_iso2", "adm1_gid"]].drop_duplicates().head()
        raise RuntimeError(f"Gallup-to-ADM1 country mismatch:\n{examples}")
    namibia = frame.loc[frame["iso2"].eq("NA"), "country_name"]
    # An unmatched left merge yields NaN rather than "", so both must be checked.
    if len(namibia) and (namibia.isna() | namibia.eq("")).all():
        raise RuntimeError("Namibia NA was lost during Gallup merges")
    return frame


def add_weights(frame: pd.DataFrame) -> pd.DataFrame:
    frame = frame.copy()
    frame["weight_raw"] = pd.to_numeric(frame["weight_raw"], errors="coerce")
    valid = frame["weight_raw"].notna() & frame["weight_raw"].gt(0)
    frame["weight_equal_country_wave"] = np.nan
    totals = frame.loc[valid].groupby(["year_wave", "iso2"])["weight_raw"].transform("sum")
    frame.loc[valid, "weight_equal_country_wave"] = frame.loc[valid, "weight_raw"] / totals
    countries_per_wave = (
        frame.loc[valid].groupby("year_wave")["iso2"].nunique().to_dict()
    )
    frame["weight_equal_country_equal_wave"] = frame["weight_equal_country_wave"]
    for year, countries in countries_per_wave.items():
        mask = frame["year_wave"].eq(year)
        frame.loc[mask, "weight_equal_country_equal_wave"] = (
            frame.loc[mask, "weight_equal_country_wave"] / countries / len(YEARS)
        )
    return frame


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--chunksize", type=int, default=25_000)
    args = parser.parse_args()
    ensure_directories()
    frame, _ = stream_raw(chunksize=args.chunksize)
    frame = merge_exposures_and_controls(frame)
    frame = add_weights(frame)
    weight_sums = (
        frame.loc[frame["weight_equal_country_wave"].notna()]
        .groupby(["year_wave", "iso2"])["weight_equal_country_wave"]
        .sum()
    )
    if not np.allclose(weight_sums, 1.0, atol=1e-10):
        raise RuntimeError("Country-normalized Gallup weights do not sum to one")
    frame.to_parquet(OUTPUT, index=False, compression="zstd")
    print(
        f"Wrote {len(frame):,} restricted Gallup rows across {frame['iso2'].nunique()} "
        f"countries to {OUTPUT}."
    )


if __name__ == "__main__":
    main()
