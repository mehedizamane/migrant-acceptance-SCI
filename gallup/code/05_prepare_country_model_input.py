#!/usr/bin/env python3
"""Prepare the restricted, fixed-scaling input for 06_fit_country_models.R."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import DERIVED, MODELS, RESTRICTED, YEARS, ensure_directories, weighted_mean, weighted_sd


INPUT = RESTRICTED / "derived" / "gallup_four_wave_sci_analytic.parquet"
OUTPUT = RESTRICTED / "derived" / "gallup_model_input.csv"
COUNTRY_CONTROLS = DERIVED / "country_controls.csv"
COUNTRY_GEO_SCI = DERIVED / "country_global_sci_geography_adjusted_2020_2026.csv"

OUTCOMES = [
    "migrant_acceptance_index",
    "migrant_acceptance_index_strict",
    "immigrants_living_in_country_good",
    "immigrant_neighbor_good",
    "immigrant_marry_close_relative_good",
]
INDIVIDUAL = [
    "age",
    "female",
    "education_level",
    "income_quintile",
    "foreign_born",
    "urbanicity",
]
MACRO = [
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


def numeric(frame: pd.DataFrame, columns: list[str]) -> None:
    for column in columns:
        frame[column] = pd.to_numeric(frame[column], errors="coerce")


def unweighted_country_scaling(frame: pd.DataFrame, rows: list[dict[str, object]]) -> None:
    country = pd.read_csv(COUNTRY_CONTROLS, keep_default_na=False)
    for column in MACRO:
        values = pd.to_numeric(country[column], errors="coerce")
        mean = float(values.mean())
        sd = float(values.std(ddof=0))
        target = f"{column}_z"
        frame[target] = (pd.to_numeric(frame[column], errors="coerce") - mean) / sd
        rows.append(
            {
                "source_variable": column,
                "standardized_variable": target,
                "mean": mean,
                "standard_deviation": sd,
                "reference": "All 2026 country-SCI origins with a valid control value",
                "weighting": "Each country receives equal weight",
            }
        )


def add_geography_adjusted_scaling(
    frame: pd.DataFrame, rows: list[dict[str, object]]
) -> pd.DataFrame:
    country = pd.read_csv(COUNTRY_GEO_SCI, keep_default_na=False)
    country_columns = [
        "iso2",
        "global_sci_geography_adjusted_2020_z",
        "global_sci_geography_adjusted_2026_z",
    ]
    frame = frame.merge(
        country[country_columns], on="iso2", how="left", validate="many_to_one"
    )
    return frame


def valid_rows(frame: pd.DataFrame, required: list[str]) -> pd.Series:
    return frame[required].notna().all(axis=1)


def main() -> None:
    ensure_directories()
    frame = pd.read_parquet(INPUT)
    numeric(
        frame,
        [
            "weight_raw",
            "year_wave",
            *OUTCOMES,
            "stranger_contact_frequency",
            *[column for column in INDIVIDUAL if column != "urbanicity"],
            *MACRO,
            "employment_to_population_ratio_15plus",
            "global_sci_2020_z",
            "global_sci_2020_all_origins_z",
            "global_sci_2026_z",
            "global_sci_2026_all_origins_z",
        ],
    )
    frame["year_wave"] = frame["year_wave"].astype("Int64")
    frame["age_10"] = (frame["age"] - 40) / 10
    frame["country_id"] = frame["iso2"].astype(str)
    frame["adm1_id"] = frame["adm1_gid"].astype(str)
    frame.loc[frame["adm1_id"].eq(""), "adm1_id"] = np.nan
    frame["wb_region"] = frame["region"].replace("", np.nan)
    frame["un_region"] = frame["un_region_9"].replace("", np.nan)

    scaling_rows: list[dict[str, object]] = []
    primary_required = [
        "migrant_acceptance_index",
        "weight_raw",
        "global_sci_2026_all_origins_z",
        *INDIVIDUAL,
    ]
    primary_eligible_sets = []
    for year in YEARS:
        wave = frame.loc[frame["year_wave"].eq(year)]
        primary_eligible_sets.append(
            set(wave.loc[valid_rows(wave, primary_required), "iso2"])
        )
    primary_countries = set.intersection(*primary_eligible_sets)

    release_required = [
        "migrant_acceptance_index",
        "weight_raw",
        "global_sci_2020_z",
        "global_sci_2026_z",
        *INDIVIDUAL,
    ]
    release_eligible_sets = []
    for year in YEARS:
        wave = frame.loc[frame["year_wave"].eq(year)]
        release_eligible_sets.append(
            set(wave.loc[valid_rows(wave, release_required), "iso2"])
        )
    release_comparison_countries = set.intersection(*release_eligible_sets)

    unweighted_country_scaling(frame, scaling_rows)
    # Employment is time-varying, so scale it separately among countries in
    # each Gallup wave. The scaling does not use respondent duplication.
    frame["employment_to_population_ratio_15plus_z"] = np.nan
    for year in YEARS:
        country_values = (
            frame.loc[frame["year_wave"].eq(year), ["iso2", "employment_to_population_ratio_15plus"]]
            .drop_duplicates("iso2")
        )
        values = country_values["employment_to_population_ratio_15plus"].dropna()
        mean = float(values.mean())
        sd = float(values.std(ddof=0))
        mask = frame["year_wave"].eq(year)
        frame.loc[mask, "employment_to_population_ratio_15plus_z"] = (
            frame.loc[mask, "employment_to_population_ratio_15plus"] - mean
        ) / sd
        scaling_rows.append(
            {
                "source_variable": "employment_to_population_ratio_15plus",
                "standardized_variable": "employment_to_population_ratio_15plus_z",
                "mean": mean,
                "standard_deviation": sd,
                "reference": f"Countries with 2026 SCI and employment data in {year}",
                "weighting": "Each country receives equal weight",
            }
        )

    frame = add_geography_adjusted_scaling(frame, scaling_rows)

    contact_required = [
        "migrant_acceptance_index",
        "weight_raw",
        "global_sci_2026_all_origins_z",
        "stranger_contact_frequency",
        *INDIVIDUAL,
    ]
    contact_reference = frame.loc[
        frame["year_wave"].eq(2022) & valid_rows(frame, contact_required)
    ].copy()
    contact_mean = weighted_mean(
        contact_reference["stranger_contact_frequency"],
        contact_reference["weight_equal_country_wave"],
    )
    contact_sd = weighted_sd(
        contact_reference["stranger_contact_frequency"],
        contact_reference["weight_equal_country_wave"],
    )
    frame["stranger_contact_z"] = (
        frame["stranger_contact_frequency"] - contact_mean
    ) / contact_sd
    scaling_rows.append(
        {
            "source_variable": "stranger_contact_frequency",
            "standardized_variable": "stranger_contact_z",
            "mean": contact_mean,
            "standard_deviation": contact_sd,
            "reference": (
                f"2022 country-level contact sample: {len(contact_reference):,} respondents "
                f"in {contact_reference['iso2'].nunique()} countries"
            ),
            "weighting": "Gallup weights normalized to equal country totals",
        }
    )
    country_contact_means = contact_reference.groupby("iso2", group_keys=False).apply(
        lambda part: np.average(
            part["stranger_contact_frequency"],
            weights=part["weight_equal_country_wave"],
        ),
        include_groups=False,
    )
    frame["stranger_contact_country_mean"] = frame["iso2"].map(country_contact_means)
    frame["stranger_contact_within"] = (
        frame["stranger_contact_frequency"] - frame["stranger_contact_country_mean"]
    )
    contact_reference = frame.loc[contact_reference.index]
    for source, target in [
        ("stranger_contact_within", "stranger_contact_within_z"),
        ("stranger_contact_country_mean", "stranger_contact_country_mean_z"),
    ]:
        mean = weighted_mean(contact_reference[source], contact_reference["weight_equal_country_wave"])
        sd = weighted_sd(contact_reference[source], contact_reference["weight_equal_country_wave"])
        frame[target] = (frame[source] - mean) / sd
        scaling_rows.append(
            {
                "source_variable": source,
                "standardized_variable": target,
                "mean": mean,
                "standard_deviation": sd,
                "reference": "2022 country-level 2026-SCI contact sample",
                "weighting": "Gallup weights normalized to equal country totals",
            }
        )

    frame["balanced_primary_2026_country"] = frame["iso2"].isin(primary_countries)
    frame["balanced_release_comparison_country"] = frame["iso2"].isin(
        release_comparison_countries
    )
    frame["contact_eligible_2026"] = False
    frame.loc[contact_reference.index, "contact_eligible_2026"] = True
    keep = [
        "country_id",
        "adm1_id",
        "year_wave",
        "weight_raw",
        "balanced_primary_2026_country",
        "balanced_release_comparison_country",
        "contact_eligible_2026",
        *OUTCOMES,
        "global_sci_2020_z",
        "global_sci_2026_z",
        "global_sci_2026_all_origins_z",
        "global_sci_geography_adjusted_2020_z",
        "global_sci_geography_adjusted_2026_z",
        "stranger_contact_z",
        "stranger_contact_within_z",
        "stranger_contact_country_mean_z",
        "age_10",
        "female",
        "education_level",
        "income_quintile",
        "foreign_born",
        "urbanicity",
        "wb_region",
        "un_region",
        *[f"{column}_z" for column in MACRO],
        "employment_to_population_ratio_15plus_z",
    ]
    # Keep every respondent with valid design weights and demographics; each
    # model then applies its own outcome/exposure complete-case restriction.
    base_keep = valid_rows(frame, ["weight_raw", *INDIVIDUAL]) & frame["weight_raw"].gt(0)
    model_input = frame.loc[base_keep, keep].copy()
    model_input.to_csv(OUTPUT, index=False, na_rep="NA_REAL")
    MODELS.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(scaling_rows).to_csv(MODELS / "model_scaling_parameters.csv", index=False)
    print(
        f"Prepared {len(model_input):,} restricted model-input rows; "
        f"primary 2026 balanced set contains {len(primary_countries)} countries; "
        f"release-comparison set contains {len(release_comparison_countries)} countries."
    )


if __name__ == "__main__":
    main()
