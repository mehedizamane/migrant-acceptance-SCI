#!/usr/bin/env python3
"""Seed source-cited historical and aggregate Gallup geography compositions.

The archived submarine-cable crosswalk contains useful many-to-many work, but
it is not itself an authoritative geographic source. This script promotes only
five partitions whose full constituent lists can be checked against cited
official publications: Ghana's former ten regions, Kenya's former eight
provinces, Nigeria's six geopolitical zones, Sri Lanka's nine provinces, and
Poland's sixteen voivodeships.

The output is an input to 08_match_gallup_regions.py. A Gallup region that
resolves to more than one ADM1 unit is excluded from the direct ADM1 linkage,
so these compositions determine which regional codes are eligible for it.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd


PACKAGE = Path(__file__).resolve().parents[1]
RAW_V2 = PACKAGE / "data" / "raw" / "geography" / "v2"
PRIOR = RAW_V2 / "prior_submarine_geography_crosswalk.csv"
OUTPUT = RAW_V2 / "authoritative_composition_overrides.csv"

COLUMNS = [
    "iso2", "year_from", "year_to", "gallup_region_variable",
    "gallup_region_code", "gallup_region_label", "adm1_gid",
    "population_overlap_2025", "boundary_definition", "boundary_vintage",
    "source_id", "source_title", "source_url", "source_access_date",
    "mapping_method", "notes",
]

SOURCES = {
    "GH": {
        "variable": "REGION_GHA",
        "codes": 10,
        "constituents": 16,
        "source_id": "ghana_gss_10_to_16_regions",
        "source_title": "Ghana Statistical Service CPI Technical Guide",
        "source_url": (
            "https://www.statsghana.gov.gh/gssmain/fileUpload/Price%20Indices/"
            "CPI_Technical_Guide_v5_Published_14102020.pdf"
        ),
        "definition": (
            "Former ten-region Ghana administrative partition concorded to "
            "the sixteen GADM 4.1 regions"
        ),
        "vintage": "Former ten-region scheme through 2018; GADM 4.1 targets",
    },
    "KE": {
        "variable": "REGION2_KEN",
        "codes": 8,
        "constituents": 47,
        "source_id": "kenya_districts_and_provinces_act_1992",
        "source_title": "Kenya Districts and Provinces Act, Third Schedule",
        "source_url": "https://new.kenyalaw.org/akn/ke/act/1992/5/eng@2022-12-31",
        "definition": (
            "Former eight-province Kenya partition concorded to the forty-seven "
            "GADM 4.1 counties"
        ),
        "vintage": "1992 province/district schedule; GADM 4.1 county targets",
    },
    "NG": {
        "variable": "REGION2_NGA",
        "codes": 6,
        "constituents": 37,
        "source_id": "nigeria_nbs_six_geopolitical_zones",
        "source_title": "Nigeria National Nutrition and Health Survey 2018, Table 3",
        "source_url": "https://www.nigerianstat.gov.ng/download/839",
        "definition": (
            "Six Nigerian geopolitical zones composed from the thirty-six states "
            "and Federal Capital Territory in GADM 4.1"
        ),
        "vintage": "NBS six-zone classification; GADM 4.1 state targets",
    },
}

MANUAL_PARTITIONS = {
    "LK": {
        "variable": "REGION_LKA",
        "source_id": "sri_lanka_dcs_nine_provinces_25_districts_2024",
        "source_title": "Sri Lanka Census of Population and Housing 2024 Preliminary Report",
        "source_url": (
            "https://www.statistics.gov.lk/Resource/en/Population/CPH_2024/"
            "Preliminary_Report.pdf"
        ),
        "definition": (
            "Nine Sri Lankan provinces composed from the twenty-five administrative "
            "districts represented as GADM 4.1 ADM1 origins in the Meta file"
        ),
        "vintage": "Nine-province structure documented for the 2024 census; GADM 4.1 district targets",
        "regions": {
            "1": ("Western", ["LKA.5_1", "LKA.7_1", "LKA.10_1"]),
            "2": ("Central", ["LKA.11_1", "LKA.16_1", "LKA.20_1"]),
            "3": ("Southern", ["LKA.6_1", "LKA.17_1", "LKA.8_1"]),
            "4": ("Northern", ["LKA.9_1", "LKA.13_1", "LKA.15_1", "LKA.19_1", "LKA.25_1"]),
            "5": ("Eastern", ["LKA.24_1", "LKA.4_1", "LKA.1_1"]),
            "6": ("Northwest", ["LKA.14_1", "LKA.22_1"]),
            "7": ("North Central", ["LKA.2_1", "LKA.21_1"]),
            "8": ("Uva", ["LKA.3_1", "LKA.18_1"]),
            "9": ("Sabaragamuwa", ["LKA.12_1", "LKA.23_1"]),
        },
    },
    "PL": {
        "variable": "REGION_POL",
        "source_id": "poland_statistics_16_voivodeships",
        "source_title": "Statistics Poland: Administrative division of Poland",
        "source_url": (
            "https://stat.gov.pl/en/regional-statistics/classification-of-territorial-units/"
            "administrative-division-of-poland/"
        ),
        "definition": (
            "Sixteen Polish voivodeships concorded one-to-one to GADM 4.1 ADM1 origins"
        ),
        "vintage": "Sixteen-voivodeship structure in force since 1999; GADM 4.1 targets",
        "regions": {
            "1": ("Lodz Voivodeship", ["POL.3_1"]),
            "2": ("Masovian Voivodeship", ["POL.7_1"]),
            "3": ("Lesser Poland Voivodeship", ["POL.6_1"]),
            "4": ("Silesian Voivodeship", ["POL.12_1"]),
            "5": ("Lublin Voivodeship", ["POL.4_1"]),
            "6": ("Subcarpathian Voivodeship", ["POL.9_1"]),
            "7": ("Swietokrzyskie Voivodeship", ["POL.13_1"]),
            "8": ("Podlaskie Voivodeship", ["POL.10_1"]),
            "9": ("Greater Poland Voivodeship", ["POL.15_1"]),
            "10": ("West Pomeranian Voivodeship", ["POL.16_1"]),
            "11": ("Lubusz Voivodeship", ["POL.5_1"]),
            "12": ("Lower Silesian Voivodeship", ["POL.1_1"]),
            "13": ("Opole Voivodeship", ["POL.8_1"]),
            "14": ("Kuyavian-Pomeranian Voivodeship", ["POL.2_1"]),
            "15": ("Warmian-Masurian Voivodeship", ["POL.14_1"]),
            "16": ("Pomeranian Voivodeship", ["POL.11_1"]),
        },
    },
}


def main() -> None:
    if not PRIOR.exists():
        raise FileNotFoundError(f"Missing archived composition input: {PRIOR}")
    prior = pd.read_csv(PRIOR, keep_default_na=False, dtype=str)
    accepted = prior.loc[
        prior["accepted"].astype(str).str.casefold().isin(["true", "1"])
        & prior["source_gadm1_id"].str.strip().ne("")
        & prior["iso2"].isin(SOURCES),
    ].copy()
    accepted = accepted.drop_duplicates(
        ["iso2", "gallup_region_variable", "gallup_region_code", "source_gadm1_id"]
    )

    rows: list[dict[str, str]] = []
    for iso2, source in SOURCES.items():
        part = accepted.loc[
            accepted["iso2"].eq(iso2)
            & accepted["gallup_region_variable"].eq(source["variable"])
        ].copy()
        observed_codes = part["gallup_region_code"].nunique()
        observed_constituents = part["source_gadm1_id"].nunique()
        if observed_codes != source["codes"] or observed_constituents != source["constituents"]:
            raise RuntimeError(
                f"{iso2} partition failed validation: expected "
                f"{source['codes']} codes/{source['constituents']} constituents; "
                f"found {observed_codes}/{observed_constituents}"
            )
        if part["source_gadm1_id"].duplicated().any():
            raise RuntimeError(f"{iso2} partition assigns a GADM1 constituent more than once")
        for item in part.itertuples(index=False):
            rows.append(
                {
                    "iso2": iso2,
                    "year_from": "",
                    "year_to": "",
                    "gallup_region_variable": source["variable"],
                    "gallup_region_code": item.gallup_region_code,
                    "gallup_region_label": item.gallup_region_label,
                    "adm1_gid": item.source_gadm1_id,
                    "population_overlap_2025": "",
                    "boundary_definition": source["definition"],
                    "boundary_vintage": source["vintage"],
                    "source_id": source["source_id"],
                    "source_title": source["source_title"],
                    "source_url": source["source_url"],
                    "source_access_date": "2026-08-16",
                    "mapping_method": "authoritative_historical_or_aggregate_composition_candidate",
                    "notes": (
                        "Composition cross-checked against the cited official partition and "
                        "the archived project crosswalk; two independent approvals required."
                    ),
                }
            )

    for iso2, source in MANUAL_PARTITIONS.items():
        constituents = [
            adm1_gid
            for _, adm1_ids in source["regions"].values()
            for adm1_gid in adm1_ids
        ]
        if len(constituents) != len(set(constituents)):
            raise RuntimeError(f"{iso2} manual partition assigns a GADM1 constituent more than once")
        for code, (label, adm1_ids) in source["regions"].items():
            for adm1_gid in adm1_ids:
                rows.append(
                    {
                        "iso2": iso2,
                        "year_from": "",
                        "year_to": "",
                        "gallup_region_variable": source["variable"],
                        "gallup_region_code": code,
                        "gallup_region_label": label,
                        "adm1_gid": adm1_gid,
                        "population_overlap_2025": "",
                        "boundary_definition": source["definition"],
                        "boundary_vintage": source["vintage"],
                        "source_id": source["source_id"],
                        "source_title": source["source_title"],
                        "source_url": source["source_url"],
                        "source_access_date": "2026-08-17",
                        "mapping_method": (
                            "authoritative_historical_or_aggregate_composition_candidate"
                        ),
                        "notes": (
                            "Composition transcribed from the cited official partition; "
                            "two independent approvals required."
                        ),
                    }
                )

    seeded = pd.DataFrame(rows, columns=COLUMNS)
    existing = pd.DataFrame(columns=COLUMNS)
    if OUTPUT.exists():
        existing = pd.read_csv(OUTPUT, keep_default_na=False, dtype=str)
        for column in COLUMNS:
            if column not in existing:
                existing[column] = ""
        managed_ids = {
            source["source_id"]
            for source in [*SOURCES.values(), *MANUAL_PARTITIONS.values()]
        }
        existing = existing.loc[~existing["source_id"].isin(managed_ids), COLUMNS]
    combined = pd.concat([existing, seeded], ignore_index=True)
    keys = [
        "iso2", "year_from", "year_to", "gallup_region_variable",
        "gallup_region_code", "adm1_gid",
    ]
    if combined.duplicated(keys).any():
        raise RuntimeError("Combined authoritative composition input contains duplicate rows")
    combined = combined.sort_values(
        ["iso2", "gallup_region_variable", "gallup_region_code", "adm1_gid"]
    )
    combined.to_csv(OUTPUT, index=False)
    print(
        f"Wrote {len(combined):,} authoritative composition rows to {OUTPUT}; "
        f"seeded {seeded['iso2'].nunique()} country partitions."
    )


if __name__ == "__main__":
    main()
