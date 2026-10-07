#!/usr/bin/env python3
"""Build national, ADM1, and dyadic control datasets with provenance fields."""

from __future__ import annotations

import sys
from pathlib import Path
from zipfile import ZipFile

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import DERIVED, RAW, ensure_directories, require_columns


COUNTRY_METADATA = RAW / "controls" / "country_metadata_reference.csv"
WDI_LONG = RAW / "controls" / "wdi_2019_2022_long.csv"
EMPLOYMENT_LONG = RAW / "controls" / "wdi_employment_to_population_ratio_country_year.csv"
MIGRANT_WORKBOOK = RAW / "controls" / "undesa_pd_2024_ims.xlsx"
REFUGEE_2022 = RAW / "controls" / "unhcr_refugee_stock_2022.csv"
UN_M49 = RAW / "controls" / "un_m49_country_regions.csv"
CEPII_ARCHIVE = RAW / "controls" / "Gravity_csv_V202211.zip"
CEPII_MEMBER = "Gravity_V202211.csv"
REGIONAL_SCI = DERIVED / "regional_sci_adm1_2026.csv"

WDI_INDICATORS = {
    "gdp_per_capita_current_usd": "NY.GDP.PCAP.CD",
    "population_total": "SP.POP.TOTL",
    "internet_users_pct": "IT.NET.USER.ZS",
    "mobile_cellular_subscriptions_per_100": "IT.CEL.SETS.P2",
    "fixed_broadband_subscriptions_per_100": "IT.NET.BBND.P2",
    "secure_internet_servers_per_million_people": "IT.NET.SECR.P6",
    "urban_population_pct": "SP.URB.TOTL.IN.ZS",
    "trade_pct_gdp": "NE.TRD.GNFS.ZS",
    "tourism_arrivals": "ST.INT.ARVL",
}


def select_wdi_values(long: pd.DataFrame) -> pd.DataFrame:
    require_columns(long, ["variable", "indicator", "iso3", "year", "value"], "WDI long file")
    pieces: list[pd.DataFrame] = []
    for variable, indicator in WDI_INDICATORS.items():
        part = long.loc[
            long["variable"].eq(variable) & long["indicator"].eq(indicator)
        ].copy()
        part["year"] = pd.to_numeric(part["year"], errors="coerce")
        part["value"] = pd.to_numeric(part["value"], errors="coerce")
        if variable == "tourism_arrivals":
            # Pre-pandemic tourism exposure. This is intentionally fixed at
            # 2019 for every Gallup wave in the fully adjusted model.
            part = part.loc[part["year"].eq(2019)]
        else:
            # 2022 is preferred; 2021, 2020, and 2019 are ordered fallbacks.
            # Drop missing values before ranking so a present-but-empty 2022 row
            # cannot mask a valid earlier observation. The selected value is held
            # constant across every Gallup wave.
            part = part.loc[part["year"].between(2019, 2022) & part["value"].notna()]
            part = part.sort_values(["iso3", "year"], ascending=[True, False])
        part = part.drop_duplicates("iso3")
        part = part.rename(columns={"value": variable, "year": f"{variable}_year"})
        pieces.append(part[["iso3", variable, f"{variable}_year"]])
    result = pieces[0]
    for part in pieces[1:]:
        result = result.merge(part, on="iso3", how="outer", validate="one_to_one")

    population_rows = long.loc[long["variable"].eq("population_total")].copy()
    population_rows["year"] = pd.to_numeric(population_rows["year"], errors="coerce")
    population_rows["value"] = pd.to_numeric(population_rows["value"], errors="coerce")
    for year in (2019, 2020, 2022):
        annual = population_rows.loc[
            population_rows["year"].eq(year), ["iso3", "value"]
        ].rename(columns={"value": f"population_total_{year}"})
        annual = annual.drop_duplicates("iso3")
        result = result.merge(annual, on="iso3", how="left", validate="one_to_one")
    return result


def parse_migrant_stock(country: pd.DataFrame) -> pd.DataFrame:
    data = pd.read_excel(MIGRANT_WORKBOOK, sheet_name="Table 1", header=10)
    destination = "Location code of destination"
    origin = "Location code of origin"
    require_columns(data, [destination, origin, 2020], "UN DESA migrant stock table")
    numeric_map = (
        country.assign(un_numeric_int=pd.to_numeric(country["un_numeric"], errors="coerce"))
        .dropna(subset=["un_numeric_int"])
        .assign(un_numeric_int=lambda frame: frame["un_numeric_int"].astype(int))
        .set_index("un_numeric_int")["iso2"]
        .to_dict()
    )
    totals = data.loc[pd.to_numeric(data[origin], errors="coerce").eq(900), [destination, 2020]].copy()
    totals["destination_code"] = pd.to_numeric(totals[destination], errors="coerce")
    totals["iso2"] = totals["destination_code"].map(numeric_map)
    totals["total_migrant_stock_2020"] = pd.to_numeric(totals[2020], errors="coerce")
    return (
        totals.dropna(subset=["iso2"])
        .groupby("iso2", as_index=False)["total_migrant_stock_2020"]
        .sum()
    )


def build_country_controls() -> pd.DataFrame:
    country = pd.read_csv(COUNTRY_METADATA, keep_default_na=False, low_memory=False)
    required = ["iso2", "iso3", "un_numeric", "country_name", "region"]
    require_columns(country, required, "country metadata")
    country = country[required].drop_duplicates("iso2").copy()
    if len(country) != 178 or "NA" not in set(country["iso2"]):
        raise RuntimeError("Country metadata must contain 178 countries including Namibia NA")

    wdi_long = pd.read_csv(WDI_LONG, keep_default_na=False)
    selected = select_wdi_values(wdi_long)
    result = country.merge(selected, on="iso3", how="left", validate="one_to_one")
    result = result.merge(parse_migrant_stock(country), on="iso2", how="left", validate="one_to_one")
    refugees = pd.read_csv(REFUGEE_2022, keep_default_na=False)
    result = result.merge(refugees, on="iso3", how="left", validate="one_to_one")

    numeric = [
        *WDI_INDICATORS,
        "population_total_2019",
        "population_total_2020",
        "population_total_2022",
        "total_migrant_stock_2020",
        "total_refugee_stock_2022",
    ]
    for column in numeric:
        result[column] = pd.to_numeric(result[column], errors="coerce")
    result["log_gdp_per_capita"] = np.log(
        result["gdp_per_capita_current_usd"].where(result["gdp_per_capita_current_usd"].gt(0))
    )
    result["log_population"] = np.log(
        result["population_total"].where(result["population_total"].gt(0))
    )
    result["log_secure_internet_servers_per_million_people"] = np.log1p(
        result["secure_internet_servers_per_million_people"]
    )
    result["tourism_arrivals_per_capita"] = (
        result["tourism_arrivals"] / result["population_total_2019"]
    )
    result["migrant_stock_share"] = (
        result["total_migrant_stock_2020"] / result["population_total_2020"]
    )
    result["refugee_stock_share"] = (
        result["total_refugee_stock_2022"].fillna(0) / result["population_total_2022"]
    )

    un_regions = pd.read_csv(UN_M49, keep_default_na=False, dtype=str)
    un_columns = [
        "iso2",
        "un_m49_region",
        "un_m49_subregion",
        "un_m49_intermediate_region",
        "un_region_9",
    ]
    require_columns(un_regions, un_columns, "UN M49 country regions")
    result = result.merge(
        un_regions[un_columns].drop_duplicates("iso2"),
        on="iso2",
        how="left",
        validate="one_to_one",
    )
    result = result.sort_values("iso2").reset_index(drop=True)
    result.to_csv(DERIVED / "country_controls.csv", index=False)

    employment = pd.read_csv(EMPLOYMENT_LONG, keep_default_na=False)
    require_columns(
        employment,
        ["iso3", "indicator_name", "indicator_code", "year_wave", "employment_to_population_ratio_15plus"],
        "employment-to-population file",
    )
    employment = employment.merge(
        country[["iso2", "iso3", "country_name"]], on="iso3", how="inner", validate="many_to_one"
    )
    employment = employment.loc[
        pd.to_numeric(employment["year_wave"], errors="coerce").isin([2016, 2019, 2022, 2023])
    ].copy()
    employment["employment_to_population_ratio_15plus"] = pd.to_numeric(
        employment["employment_to_population_ratio_15plus"], errors="coerce"
    )
    employment.to_csv(DERIVED / "country_wave_employment_controls.csv", index=False)
    return result


def build_adm1_controls() -> pd.DataFrame:
    region = pd.read_csv(REGIONAL_SCI, keep_default_na=False)
    columns = [
        "analysis_iso2",
        "sci_country_iso2",
        "adm1_gid",
        "gadm_country_name",
        "adm1_name",
        "ghs_pop_2025",
        "adm1_area_km2",
        "adm1_population_density_per_km2",
        "adm1_log_population",
        "adm1_log_population_density",
        "population_valid",
    ]
    require_columns(region, columns, "regional SCI data")
    controls = region[columns].rename(columns={"analysis_iso2": "iso2"}).copy()
    controls.to_csv(DERIVED / "adm1_controls.csv", index=False)
    return controls


def build_cepii_dyadic_controls() -> pd.DataFrame:
    columns = [
        "year",
        "iso3_o",
        "iso3_d",
        "country_exists_o",
        "country_exists_d",
        "distw_harmonic",
        "dist",
        "contig",
        "comlang_off",
        "comlang_ethno",
        "comcol",
        "col45",
        "col_dep_ever",
    ]
    parts: list[pd.DataFrame] = []
    with ZipFile(CEPII_ARCHIVE) as archive:
        if CEPII_MEMBER not in archive.namelist():
            raise RuntimeError(f"CEPII archive lacks {CEPII_MEMBER}")
        with archive.open(CEPII_MEMBER) as handle:
            for chunk in pd.read_csv(handle, usecols=columns, chunksize=500_000, low_memory=False):
                use = chunk.loc[
                    chunk["year"].eq(2020)
                    & chunk["country_exists_o"].eq(1)
                    & chunk["country_exists_d"].eq(1)
                ].copy()
                if not use.empty:
                    parts.append(use)
    result = pd.concat(parts, ignore_index=True).drop_duplicates(["iso3_o", "iso3_d"])
    required = ["distw_harmonic", "contig", "comlang_off", "comlang_ethno", "col_dep_ever"]
    result = result.dropna(subset=required)
    for column in ["contig", "comlang_off", "comlang_ethno", "comcol", "col45", "col_dep_ever"]:
        result[column] = pd.to_numeric(result[column], errors="coerce")
    result["log_distance_km"] = np.log(
        pd.to_numeric(result["distw_harmonic"], errors="coerce").where(
            pd.to_numeric(result["distw_harmonic"], errors="coerce").gt(0)
        )
    )
    result.to_csv(DERIVED / "cepii_country_dyadic_controls_2020.csv.gz", index=False, compression="gzip")
    return result


def main() -> None:
    ensure_directories()
    country = build_country_controls()
    adm1 = build_adm1_controls()
    cepii = build_cepii_dyadic_controls()
    print(
        f"Built controls for {len(country)} countries, {len(adm1)} ADM1 regions, "
        f"and {len(cepii):,} CEPII country dyads."
    )


if __name__ == "__main__":
    main()
