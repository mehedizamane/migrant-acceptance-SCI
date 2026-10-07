#!/usr/bin/env python3
"""Match observed Gallup regional codes to GADM 4.1 ADM1 units.

For every country, wave, base REGION_* variable, and observed region code, this
script records the candidate GADM1 unit(s) and the matching rule that produced
them: official ADM1 code, primary or official alternate GADM1 name, a name after
removing administrative-type words, an ISO 3166-2 subdivision name, a previously
verified exact match, or an exact lower-level GADM unit nested in one GADM1
parent. Documented multi-region compositions (07_seed_region_compositions.py)
are recorded as such.

The output is a restricted candidate inventory. 09_link_gallup_regions_to_adm1.py
applies the deterministic one-to-one rules that define the paper's regional
linkage (SI Section 2.1, SI Table 1). Fuzzy-similarity candidates are recorded
so they can be counted as exclusions, but they are never accepted.
"""

from __future__ import annotations

import argparse
import difflib
import hashlib
import math
import re
import sqlite3
import subprocess
import sys
import unicodedata
from collections import defaultdict
from pathlib import Path
from zipfile import ZipFile

import geopandas as gpd
import numpy as np
import pandas as pd
import pycountry
import pyogrio
from unidecode import unidecode

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import RAW, RESTRICTED, YEARS, ensure_directories


ARCHIVE = RESTRICTED / "raw" / "Gallup_World_Poll_022026_ALL_WAVES.zip"
MEMBER = "The_Gallup_022026.dat"
SYSFILE = RESTRICTED / "raw" / "GWP_022026_SysfileInfo_sav.xlsx"
META_ARCHIVE = RAW / "meta" / "all_region_to_country.zip"
META_MEMBER = "gadm1_to_country.csv"
GADM_COMPLETE = RAW / "geography" / "gadm_410_complete.gpkg"
GADM_SUBSET_V1 = RAW / "geography" / "gadm1_metadata.gpkg"
V1_CROSSWALK = RAW / "geography" / "v1" / "gallup_to_gadm1_crosswalk_v1.csv"

RAW_V2 = RAW / "geography" / "v2"
PRIOR_COMPOSITION_CROSSWALK = RAW_V2 / "prior_submarine_geography_crosswalk.csv"
AUTHORITATIVE_COMPOSITIONS = RAW_V2 / "authoritative_composition_overrides.csv"
RESTRICTED_V2 = RESTRICTED / "derived" / "geography_v2"

EXPECTED_TOTAL = 3_028_581
EXPECTED_WAVE_ROWS = {2016: 150_724, 2019: 176_253, 2022: 143_686, 2023: 145_702}
EXPECTED_REGION_VARIABLES = 299
EXPECTED_META_ROWS = 563_904
EXPECTED_META_ORIGINS = 3_168
EXPECTED_META_DESTINATIONS = 178

LEGACY_GEOMETRY_IDS = {
    "IND.14_1": "Z01.14_1",
    "PAK.1_1": "Z06.1_1",
    "PAK.6_1": "Z06.6_1",
}

INVALID_LABEL_RE = re.compile(
    r"(?:^|\b)(?:dk|don.?t know|refused|missing|not applicable|no answer)(?:\b|$)", re.I
)
INVALID_REGION_CODES = {
    "98", "99", "997", "998", "999", "9997", "9998", "9999",
}
REGION_RE = re.compile(r"^REGION(?P<level>[2-5]?)_[A-Z0-9]+$")
ADMIN_WORDS = {
    "administrative", "area", "autonomous", "canton", "capital", "city", "county",
    "department", "departement", "departamento", "district", "division", "emirate",
    "federal", "governorate", "grad", "kanton", "krai", "kray", "municipality",
    "municipio", "oblast", "of", "okrug", "prefecture", "provincia", "province",
    "rayon", "region", "regional", "republic", "special", "state", "territory",
    "the", "union", "voivodeship", "wilaya", "zone", "zupanija",
}
GENERIC_BROAD_NAMES = {
    "central", "center", "coast", "east", "eastern", "interior", "islands",
    "metropolitan", "north", "northeast", "northern", "northwest", "other",
    "south", "southeast", "southern", "southwest", "west", "western",
}


def canonical_code(value: object) -> str:
    text = str(value).strip()
    if not text or text.lower() in {"nan", "none", "<na>", "."}:
        return ""
    try:
        number = float(text)
        if math.isfinite(number) and number.is_integer():
            return str(int(number))
    except ValueError:
        pass
    return text


def normalize_name(value: object) -> str:
    text = unicodedata.normalize("NFKD", str(value))
    text = "".join(ch for ch in text if not unicodedata.combining(ch))
    text = unidecode(text)
    text = text.casefold().replace("&", " and ")
    text = re.sub(r"[’'`]", "", text)
    text = re.sub(r"[^a-z0-9]+", " ", text)
    return " ".join(text.split())


def normalize_admin_core(value: object) -> str:
    text = re.sub(r"\([^)]*\)", " ", str(value))
    text = re.sub(r"^\s*[A-Za-z]{2,3}\d{1,3}[A-Za-z]?(?:\s*[-:]\s*|\s+)", " ", text)
    tokens = [token for token in normalize_name(text).split() if token not in ADMIN_WORDS]
    return " ".join(tokens)


def split_aliases(value: object) -> list[str]:
    if value is None or pd.isna(value):
        return []
    text = str(value).strip()
    if not text:
        return []
    values = re.split(r"\s*\|\s*|\s*;\s*|\s*/\s*", text)
    return [item.strip() for item in values if item.strip()]


def stable_id(*values: object, length: int = 16) -> str:
    payload = "|".join(str(value) for value in values)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()[:length]


def ensure_v2_directories() -> None:
    ensure_directories()
    for path in [RAW_V2, RESTRICTED_V2]:
        path.mkdir(parents=True, exist_ok=True)

    if not AUTHORITATIVE_COMPOSITIONS.exists():
        pd.DataFrame(
            columns=[
                "iso2", "year_from", "year_to", "gallup_region_variable",
                "gallup_region_code", "gallup_region_label", "adm1_gid",
                "population_overlap_2025", "boundary_definition",
                "boundary_vintage", "source_id", "source_title", "source_url",
                "source_access_date", "mapping_method", "notes",
            ]
        ).to_csv(AUTHORITATIVE_COMPOSITIONS, index=False)


def parse_gallup_dictionary() -> tuple[pd.DataFrame, pd.DataFrame]:
    variable_info = pd.read_excel(SYSFILE, sheet_name="Variable Info", header=1)
    variable_info = variable_info.rename(columns={"Name": "variable", "Label": "variable_label"})
    variable_info["variable"] = variable_info["variable"].astype(str).str.strip()
    inventory = variable_info.loc[
        variable_info["variable"].str.match(REGION_RE),
        ["variable", "variable_label", "Position", "Measurement Level", "Format"],
    ].copy()
    inventory["region_level"] = inventory["variable"].str.extract(REGION_RE)["level"].replace("", "1").astype(int)
    if len(inventory) != EXPECTED_REGION_VARIABLES:
        raise RuntimeError(
            f"Expected {EXPECTED_REGION_VARIABLES} REGION through REGION5 fields; found {len(inventory)}"
        )

    raw = pd.read_excel(SYSFILE, sheet_name="Value Labels", header=None)
    rows: list[dict[str, str]] = []
    current = ""
    region_set = set(inventory["variable"])
    for first, second, third in raw.iloc[2:, :3].itertuples(index=False, name=None):
        if pd.notna(first) and str(first).strip():
            current = str(first).strip()
        if current not in region_set or pd.isna(second):
            continue
        code = canonical_code(second)
        label = "" if pd.isna(third) else str(third).strip()
        if code:
            rows.append({"variable": current, "gallup_region_code": code, "gallup_region_label": label})
    labels = pd.DataFrame(rows).drop_duplicates(["variable", "gallup_region_code"], keep="last")
    inventory = inventory.sort_values(["region_level", "variable"]).reset_index(drop=True)
    return inventory, labels


def substantive_codebook_categories(labels: pd.DataFrame) -> pd.DataFrame:
    """Return declared, labeled Gallup geography categories that are substantive."""

    categories = labels[["variable", "gallup_region_code", "gallup_region_label"]].copy()
    categories["gallup_region_code"] = categories["gallup_region_code"].map(canonical_code)
    categories["gallup_region_label"] = categories["gallup_region_label"].fillna("").astype(str).str.strip()
    categories = categories.loc[
        categories["gallup_region_code"].ne("")
        & categories["gallup_region_label"].ne("")
        & ~categories["gallup_region_code"].isin(INVALID_REGION_CODES)
        & ~categories["gallup_region_label"].str.contains(INVALID_LABEL_RE, na=False)
    ].drop_duplicates(["variable", "gallup_region_code"], keep="last")
    return categories.reset_index(drop=True)


def augment_counts_with_codebook_categories(
    counts: pd.DataFrame,
    labels: pd.DataFrame,
) -> pd.DataFrame:
    """Add zero-count rows for substantive categories declared by each used field.

    Gallup's realized sample need not contain respondents in every valid survey
    region. Partition validation therefore uses the complete codebook category
    list, while respondent coverage continues to use observed counts only.
    """

    key_columns = [
        "iso2", "year_wave", "gallup_region_variable", "gallup_region_code",
    ]
    required = key_columns + ["gallup_region_label", "respondents", "weight_sum"]
    missing = set(required) - set(counts.columns)
    if missing:
        raise RuntimeError(f"Gallup geography count cache lacks columns: {sorted(missing)}")

    observed = counts[required].copy()
    observed["gallup_region_code"] = observed["gallup_region_code"].map(canonical_code)
    observed["respondents"] = pd.to_numeric(observed["respondents"], errors="coerce").fillna(0).astype(int)
    observed["weight_sum"] = pd.to_numeric(observed["weight_sum"], errors="coerce").fillna(0.0)
    observed = (
        observed.groupby(key_columns, as_index=False, dropna=False)
        .agg(
            gallup_region_label=(
                "gallup_region_label",
                lambda values: next((str(value).strip() for value in values if str(value).strip()), ""),
            ),
            respondents=("respondents", "sum"),
            weight_sum=("weight_sum", "sum"),
        )
    )

    field_keys = observed[["iso2", "year_wave", "gallup_region_variable"]].drop_duplicates()
    declared = substantive_codebook_categories(labels).rename(
        columns={"variable": "gallup_region_variable"}
    )
    declared = field_keys.merge(declared, on="gallup_region_variable", how="inner")
    declared["declared_in_codebook"] = True

    universe = pd.concat(
        [
            observed[key_columns],
            declared[key_columns],
        ],
        ignore_index=True,
    ).drop_duplicates(key_columns)
    augmented = universe.merge(observed, on=key_columns, how="left", validate="one_to_one")
    augmented = augmented.merge(
        declared[key_columns + ["gallup_region_label", "declared_in_codebook"]].rename(
            columns={"gallup_region_label": "codebook_region_label"}
        ),
        on=key_columns,
        how="left",
        validate="one_to_one",
    )
    augmented["respondents"] = augmented["respondents"].fillna(0).astype(int)
    augmented["weight_sum"] = augmented["weight_sum"].fillna(0.0)
    augmented["observed_in_wave"] = augmented["respondents"].gt(0)
    augmented["declared_in_codebook"] = augmented["declared_in_codebook"].fillna(False).astype(bool)
    augmented["gallup_region_label"] = augmented["codebook_region_label"].where(
        augmented["codebook_region_label"].fillna("").astype(str).str.strip().ne(""),
        augmented["gallup_region_label"],
    ).fillna("")
    augmented = augmented.drop(columns="codebook_region_label")
    return augmented.sort_values(key_columns).reset_index(drop=True)


def stream_region_counts(
    inventory: pd.DataFrame,
    labels: pd.DataFrame,
    chunksize: int,
    rebuild: bool,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    counts_path = RESTRICTED_V2 / "gallup_country_wave_region_code_counts.csv"
    totals_path = RESTRICTED_V2 / "gallup_country_wave_geography_totals.csv"
    country_path = RESTRICTED_V2 / "gallup_country_identifiers.csv"
    if not rebuild and counts_path.exists() and totals_path.exists() and country_path.exists():
        counts = pd.read_csv(counts_path, keep_default_na=False, dtype=str)
        counts = augment_counts_with_codebook_categories(counts, labels)
        counts.to_csv(counts_path, index=False)
        return (
            counts,
            pd.read_csv(totals_path, keep_default_na=False),
            pd.read_csv(country_path, keep_default_na=False, dtype=str),
        )

    variables = inventory["variable"].tolist()
    variables_by_iso3: defaultdict[str, list[str]] = defaultdict(list)
    for variable in variables:
        variables_by_iso3[variable.rsplit("_", 1)[-1]].append(variable)
    usecols = ["COUNTRY_ISO2", "COUNTRY_ISO3", "COUNTRYNEW", "YEAR_WAVE", "WGT"] + variables
    label_map = labels.set_index(["variable", "gallup_region_code"])["gallup_region_label"].to_dict()
    process = subprocess.Popen(["unzip", "-p", str(ARCHIVE), MEMBER], stdout=subprocess.PIPE)
    if process.stdout is None:
        raise RuntimeError("Could not stream Gallup archive")
    total_rows = 0
    wave_counts = {year: 0 for year in YEARS}
    code_accumulator: defaultdict[tuple[str, int, str, str], list[float]] = defaultdict(
        lambda: [0.0, 0.0]
    )
    total_accumulator: defaultdict[tuple[str, int], list[float]] = defaultdict(
        lambda: [0.0, 0.0, 0.0, 0.0]
    )
    country_rows: dict[tuple[str, str], str] = {}
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
        for chunk_number, chunk in enumerate(reader, start=1):
            total_rows += len(chunk)
            year = pd.to_numeric(chunk["YEAR_WAVE"], errors="coerce").astype("Int64")
            keep = year.isin(YEARS)
            if not keep.any():
                continue
            use = chunk.loc[keep].copy()
            use["_year"] = year.loc[keep].astype(int)
            use["_weight"] = pd.to_numeric(use["WGT"], errors="coerce")
            use["_valid_weight"] = use["_weight"].where(use["_weight"].gt(0), 0.0)
            for observed_year, count in use["_year"].value_counts().items():
                wave_counts[int(observed_year)] += int(count)
            for iso2, iso3, name in use[["COUNTRY_ISO2", "COUNTRY_ISO3", "COUNTRYNEW"]].drop_duplicates().itertuples(index=False, name=None):
                country_rows[(str(iso2).strip(), str(iso3).strip())] = str(name).strip()

            # Each geography field is country-specific and therefore sparse.
            # Looping over columns avoids materializing a 600k x 299 long table.
            any_valid = pd.Series(False, index=use.index)
            active_variables = sorted(
                {
                    variable
                    for iso3 in use["COUNTRY_ISO3"].astype(str).str.strip().unique()
                    for variable in variables_by_iso3.get(iso3, [])
                }
            )
            for variable in active_variables:
                codes = use[variable].map(canonical_code)
                valid = codes.ne("")
                if not valid.any():
                    continue
                labels_for_codes = codes.map(lambda code: label_map.get((variable, code), ""))
                invalid = labels_for_codes.str.contains(INVALID_LABEL_RE, na=False) | codes.isin(
                    INVALID_REGION_CODES
                )
                valid &= ~invalid
                if not valid.any():
                    continue
                any_valid |= valid
                grouped = (
                    pd.DataFrame(
                        {
                            "iso2": use.loc[valid, "COUNTRY_ISO2"].astype(str).str.strip(),
                            "year_wave": use.loc[valid, "_year"].astype(int),
                            "code": codes.loc[valid],
                            "weight": use.loc[valid, "_valid_weight"],
                        }
                    )
                    .groupby(["iso2", "year_wave", "code"], sort=False)
                    .agg(respondents=("code", "size"), weight_sum=("weight", "sum"))
                    .reset_index()
                )
                for iso2, observed_year, code, respondents, weight_sum in grouped.itertuples(index=False, name=None):
                    bucket = code_accumulator[(str(iso2), int(observed_year), variable, str(code))]
                    bucket[0] += int(respondents)
                    bucket[1] += float(weight_sum)

            grouped_totals = use.groupby(["COUNTRY_ISO2", "_year"], sort=False)
            for (iso2, observed_year), index in grouped_totals.groups.items():
                bucket = total_accumulator[(str(iso2).strip(), int(observed_year))]
                idx = list(index)
                bucket[0] += len(idx)
                bucket[1] += float(use.loc[idx, "_valid_weight"].sum())
                bucket[2] += int(any_valid.loc[idx].sum())
                bucket[3] += float(use.loc[idx[0:0], "_valid_weight"].sum()) if not idx else float(
                    use.loc[[i for i in idx if bool(any_valid.loc[i])], "_valid_weight"].sum()
                )
            if chunk_number % 10 == 0:
                print(f"  scanned {total_rows:,} Gallup rows", flush=True)
    finally:
        process.stdout.close()
        process.wait()
    if process.returncode != 0:
        raise RuntimeError(f"Gallup unzip stream exited with status {process.returncode}")
    if total_rows != EXPECTED_TOTAL or wave_counts != EXPECTED_WAVE_ROWS:
        raise RuntimeError(f"Gallup row counts changed: total={total_rows}; waves={wave_counts}")

    count_rows = []
    for (iso2, year_wave, variable, code), (respondents, weight_sum) in code_accumulator.items():
        count_rows.append(
            {
                "iso2": iso2,
                "year_wave": year_wave,
                "gallup_region_variable": variable,
                "gallup_region_code": code,
                "gallup_region_label": label_map.get((variable, code), ""),
                "respondents": int(respondents),
                "weight_sum": weight_sum,
            }
        )
    counts = augment_counts_with_codebook_categories(pd.DataFrame(count_rows), labels)
    totals = pd.DataFrame(
        [
            {
                "iso2": iso2,
                "year_wave": year_wave,
                "all_wave_respondents": int(values[0]),
                "all_wave_weight_sum": values[1],
                "respondents_with_any_region_code": int(values[2]),
                "weight_with_any_region_code": values[3],
            }
            for (iso2, year_wave), values in total_accumulator.items()
        ]
    )
    countries = pd.DataFrame(
        [
            {"iso2": iso2, "iso3": iso3, "country_name_gallup": name}
            for (iso2, iso3), name in country_rows.items()
        ]
    ).drop_duplicates("iso2")
    counts.to_csv(counts_path, index=False)
    totals.to_csv(totals_path, index=False)
    countries.to_csv(country_path, index=False)
    return counts, totals, countries


def read_meta_origins() -> tuple[pd.DataFrame, pd.DataFrame]:
    with ZipFile(META_ARCHIVE) as archive:
        with archive.open(META_MEMBER) as handle:
            dyads = pd.read_csv(
                handle,
                keep_default_na=False,
                dtype={"user_country": str, "friend_country": str, "user_region": str},
            )
    if len(dyads) != EXPECTED_META_ROWS:
        raise RuntimeError(f"Expected {EXPECTED_META_ROWS:,} Meta dyads; found {len(dyads):,}")
    origins = dyads[["user_region", "user_country"]].drop_duplicates().rename(
        columns={"user_region": "adm1_gid", "user_country": "sci_country_iso2"}
    )
    if len(origins) != EXPECTED_META_ORIGINS:
        raise RuntimeError(f"Expected {EXPECTED_META_ORIGINS:,} Meta origins; found {len(origins):,}")
    destination_counts = dyads.groupby("user_region")["friend_country"].nunique()
    if not destination_counts.eq(EXPECTED_META_DESTINATIONS).all():
        raise RuntimeError("Every Meta GADM1 origin must have 178 country destinations")
    return origins, dyads


def read_gadm_attributes(origin_ids: set[str]) -> pd.DataFrame:
    geometry_ids = {LEGACY_GEOMETRY_IDS.get(gid, gid) for gid in origin_ids}
    placeholders = ",".join("?" for _ in geometry_ids)
    columns = [
        "GID_0", "NAME_0", "GID_1", "NAME_1", "VARNAME_1", "NL_NAME_1",
        "ISO_1", "HASC_1", "CC_1", "TYPE_1", "ENGTYPE_1", "VALIDFR_1",
        "GID_2", "NAME_2", "VARNAME_2", "NL_NAME_2", "HASC_2", "CC_2",
        "GID_3", "NAME_3", "VARNAME_3", "NL_NAME_3", "HASC_3", "CC_3",
        "GID_4", "NAME_4", "VARNAME_4", "CC_4", "GID_5", "NAME_5", "CC_5",
        "REGION", "VARREGION", "COUNTRY", "CONTINENT", "SUBCONT",
    ]
    query = f"SELECT {','.join(columns)} FROM gadm_410 WHERE GID_1 IN ({placeholders})"
    connection = sqlite3.connect(GADM_COMPLETE)
    try:
        attributes = pd.read_sql_query(query, connection, params=sorted(geometry_ids), dtype=str)
    finally:
        connection.close()
    inverse = {geometry_gid: source_gid for source_gid, geometry_gid in LEGACY_GEOMETRY_IDS.items()}
    attributes["geometry_gid"] = attributes["GID_1"]
    attributes["source_adm1_gid"] = attributes["GID_1"].map(inverse).fillna(attributes["GID_1"])
    missing = origin_ids.difference(attributes["source_adm1_gid"].dropna().unique())
    if missing:
        raise RuntimeError(f"Direct Meta origins absent from complete GADM attributes: {sorted(missing)}")
    return attributes


def query_union_geometry(field: str, identifiers: list[str]) -> gpd.GeoDataFrame:
    if not identifiers:
        return gpd.GeoDataFrame(columns=["boundary_gid", "geometry"], geometry="geometry", crs="EPSG:4326")
    quoted = ",".join("'" + value.replace("'", "''") + "'" for value in identifiers)
    sql = (
        f"SELECT {field} AS boundary_gid, ST_Union(geom) AS geom "
        f"FROM gadm_410 WHERE {field} IN ({quoted}) GROUP BY {field}"
    )
    return pyogrio.read_dataframe(GADM_COMPLETE, sql=sql, sql_dialect="SQLITE")


def build_canonical_adm1(
    origins: pd.DataFrame,
    attributes: pd.DataFrame,
    countries: pd.DataFrame,
) -> gpd.GeoDataFrame:
    existing = gpd.read_file(GADM_SUBSET_V1, layer="gadm1_metadata")
    existing = existing.rename(columns={"GID_1": "geometry_gid"})
    inverse = {geometry_gid: source_gid for source_gid, geometry_gid in LEGACY_GEOMETRY_IDS.items()}
    existing["adm1_gid"] = existing["geometry_gid"].map(inverse).fillna(existing["geometry_gid"])
    existing = existing.loc[existing["adm1_gid"].isin(origins["adm1_gid"])].copy()
    missing = sorted(set(origins["adm1_gid"]).difference(existing["adm1_gid"]))
    missing_geometry_ids = [LEGACY_GEOMETRY_IDS.get(gid, gid) for gid in missing]
    additions = query_union_geometry("GID_1", missing_geometry_ids)
    if len(additions) != len(missing_geometry_ids):
        found = set(additions["boundary_gid"])
        raise RuntimeError(f"Could not build GADM1 geometries for {sorted(set(missing_geometry_ids)-found)}")
    additions = additions.rename(columns={"boundary_gid": "geometry_gid"})
    additions["adm1_gid"] = additions["geometry_gid"].map(inverse).fillna(additions["geometry_gid"])
    geometry = pd.concat(
        [existing[["adm1_gid", "geometry_gid", "geometry"]], additions[["adm1_gid", "geometry_gid", "geometry"]]],
        ignore_index=True,
    )
    geometry = gpd.GeoDataFrame(geometry, geometry="geometry", crs="EPSG:4326")
    if geometry["adm1_gid"].nunique() != EXPECTED_META_ORIGINS:
        raise RuntimeError("Canonical geometry does not cover all 3,168 Meta origins")

    adm1_fields = [
        "geometry_gid", "GID_0", "NAME_0", "NAME_1", "VARNAME_1", "NL_NAME_1",
        "ISO_1", "HASC_1", "CC_1", "TYPE_1", "ENGTYPE_1", "VALIDFR_1",
        "COUNTRY", "CONTINENT", "SUBCONT",
    ]

    def modal_nonblank(series: pd.Series) -> str:
        values = series.fillna("").astype(str).str.strip()
        values = values.loc[values.ne("")]
        if values.empty:
            return ""
        counts = values.value_counts()
        return sorted(counts.loc[counts.eq(counts.max())].index)[0]

    # GADM repeats ADM1 fields on every lower-level record. Use the modal
    # nonblank value rather than an arbitrary first row so rare denormalized
    # parent-label errors cannot redefine a Meta origin.
    primary = (
        attributes.groupby("source_adm1_gid", as_index=False)[adm1_fields]
        .agg(modal_nonblank)
    )
    metadata = origins.merge(
        primary[
            [
                "source_adm1_gid", "geometry_gid", "GID_0", "NAME_0", "NAME_1",
                "VARNAME_1", "NL_NAME_1", "ISO_1", "HASC_1", "CC_1", "TYPE_1",
                "ENGTYPE_1", "VALIDFR_1", "COUNTRY", "CONTINENT", "SUBCONT",
            ]
        ],
        left_on="adm1_gid",
        right_on="source_adm1_gid",
        how="left",
        validate="one_to_one",
    ).drop(columns="source_adm1_gid")

    # Survey-country identifiers follow Gallup when a GADM sovereign/territory
    # has its own Gallup country code. Remaining origins retain Meta's country
    # grouping. Northern Cyprus is the one documented legacy bridge in v1.
    iso3_to_iso2 = countries.set_index("iso3")["iso2"].to_dict()
    metadata["survey_iso2"] = metadata["GID_0"].map(iso3_to_iso2).fillna(metadata["sci_country_iso2"])
    if V1_CROSSWALK.exists():
        v1 = pd.read_csv(V1_CROSSWALK, keep_default_na=False, dtype=str)
        accepted = v1.loc[
            v1["match_status"].eq("accepted") & v1["adm1_gid"].str.strip().ne(""),
            ["adm1_gid", "iso2"],
        ].drop_duplicates()
        conflicts = accepted.groupby("adm1_gid")["iso2"].nunique()
        if conflicts.gt(1).any():
            raise RuntimeError("Version-1 territory bridge contains conflicting survey countries")
        bridge = accepted.drop_duplicates("adm1_gid").set_index("adm1_gid")["iso2"]
        metadata["survey_iso2"] = metadata["adm1_gid"].map(bridge).fillna(metadata["survey_iso2"])
    metadata["survey_country_source"] = np.where(
        metadata["survey_iso2"].eq(metadata["sci_country_iso2"]),
        "Meta origin-country identifier",
        "Gallup ISO3 or documented v1 territory bridge",
    )
    canonical = geometry.merge(metadata, on=["adm1_gid", "geometry_gid"], how="left", validate="one_to_one")
    return canonical


def build_alias_indexes(
    attributes: pd.DataFrame,
    canonical: gpd.GeoDataFrame,
) -> dict[str, object]:
    country_by_gid = canonical.set_index("adm1_gid")["survey_iso2"].to_dict()
    attributes = attributes.copy()
    attributes["survey_iso2"] = attributes["source_adm1_gid"].map(country_by_gid)

    adm1_aliases: defaultdict[tuple[str, str], set[str]] = defaultdict(set)
    adm1_core_aliases: defaultdict[tuple[str, str], set[str]] = defaultdict(set)
    adm1_primary: defaultdict[tuple[str, str], set[str]] = defaultdict(set)
    adm1_codes: defaultdict[tuple[str, str], set[str]] = defaultdict(set)
    iso_subdivision_aliases: defaultdict[tuple[str, str], set[tuple[str, str, str]]] = defaultdict(set)
    iso_subdivision_core_aliases: defaultdict[tuple[str, str], set[tuple[str, str, str]]] = defaultdict(set)
    lower_aliases: defaultdict[tuple[str, str], set[tuple[str, str, int]]] = defaultdict(set)
    aggregate_aliases: defaultdict[tuple[str, str], set[str]] = defaultdict(set)
    label_by_gid = canonical.set_index("adm1_gid")["NAME_1"].fillna("").astype(str).to_dict()

    # ADM1 aliases and identifiers must come from the canonical one-row-per-
    # Meta-origin table. Repeated lower-level GADM rows occasionally contain
    # inconsistent denormalized parent labels and must not create extra ADM1
    # aliases (for example, a small number of England rows carry "Wales" in
    # NAME_1 even though Wales is a separate ADM1 origin).
    for row in canonical.itertuples(index=False):
        iso2 = str(row.survey_iso2)
        adm1_gid = str(row.adm1_gid)
        if not iso2 or iso2 == "nan":
            continue
        primary_names = [getattr(row, "NAME_1", "")]
        aliases = primary_names + split_aliases(getattr(row, "VARNAME_1", "")) + split_aliases(getattr(row, "NL_NAME_1", ""))
        for name in primary_names:
            normalized = normalize_name(name)
            if normalized:
                adm1_primary[(iso2, normalized)].add(adm1_gid)
        for name in aliases:
            normalized = normalize_name(name)
            if normalized:
                adm1_aliases[(iso2, normalized)].add(adm1_gid)
                core = normalize_admin_core(name)
                if core:
                    adm1_core_aliases[(iso2, core)].add(adm1_gid)
        for code_field in ["GID_1", "ISO_1", "HASC_1", "CC_1"]:
            code = str(getattr(row, code_field, "") or "").strip()
            if code:
                adm1_codes[(iso2, normalize_name(code))].add(adm1_gid)

    # The complete hierarchy is used only for lower-level official names and
    # GADM's named aggregate fields. Each lower unit retains its unique ADM1
    # parent from the Meta origin universe.
    for row in attributes.itertuples(index=False):
        iso2 = str(row.survey_iso2)
        adm1_gid = str(row.source_adm1_gid)
        if not iso2 or iso2 == "nan":
            continue
        for level in range(2, 6):
            gid = str(getattr(row, f"GID_{level}", "") or "").strip()
            if not gid:
                continue
            names = [
                getattr(row, f"NAME_{level}", ""),
                *split_aliases(getattr(row, f"VARNAME_{level}", "")),
                *split_aliases(getattr(row, f"NL_NAME_{level}", "")),
            ]
            for name in names:
                normalized = normalize_name(name)
                if normalized:
                    lower_aliases[(iso2, normalized)].add((adm1_gid, gid, level))
        for field in ["REGION", "VARREGION"]:
            for name in split_aliases(getattr(row, field, "")):
                normalized = normalize_name(name)
                if normalized:
                    aggregate_aliases[(iso2, normalized)].add(adm1_gid)

    # ISO 3166-2 supplies an independent, official subdivision-name bridge.
    # A candidate is created only when the official ISO code maps to exactly
    # one ADM1 in the Meta origin universe through GADM's ISO_1 field. The
    # Gallup label must then match the ISO name exactly after deterministic
    # normalization. These mappings still require two human approvals because
    # Gallup's regional boundary vintage is not established by the name alone.
    iso_code_to_gids: defaultdict[tuple[str, str], set[str]] = defaultdict(set)
    for row in canonical.itertuples(index=False):
        iso2 = str(row.survey_iso2)
        iso_code = str(getattr(row, "ISO_1", "") or "").strip().upper()
        if iso2 and iso2 != "nan" and iso_code:
            iso_code_to_gids[(iso2, iso_code)].add(str(row.adm1_gid))
    for subdivision in pycountry.subdivisions:
        iso2 = str(getattr(subdivision, "country_code", "") or "").strip()
        iso_code = str(getattr(subdivision, "code", "") or "").strip().upper()
        official_name = str(getattr(subdivision, "name", "") or "").strip()
        gids = iso_code_to_gids.get((iso2, iso_code), set())
        if len(gids) != 1 or not official_name:
            continue
        candidate = (next(iter(gids)), iso_code, official_name)
        normalized = normalize_name(official_name)
        if normalized:
            iso_subdivision_aliases[(iso2, normalized)].add(candidate)
        core = normalize_admin_core(official_name)
        if core:
            iso_subdivision_core_aliases[(iso2, core)].add(candidate)

    return {
        "adm1_aliases": adm1_aliases,
        "adm1_core_aliases": adm1_core_aliases,
        "adm1_primary": adm1_primary,
        "adm1_codes": adm1_codes,
        "iso_subdivision_aliases": iso_subdivision_aliases,
        "iso_subdivision_core_aliases": iso_subdivision_core_aliases,
        "lower_aliases": lower_aliases,
        "aggregate_aliases": aggregate_aliases,
        "label_by_gid": label_by_gid,
        "country_by_gid": country_by_gid,
    }


def load_documented_compositions(
    canonical: gpd.GeoDataFrame,
) -> tuple[dict[tuple[str, str, str], list[dict[str, object]]], dict[tuple[str, str, str], list[dict[str, object]]]]:
    """Load explicit many-to-many compositions without treating them as approvals.

    Authoritative overrides must carry a source and boundary definition. The
    archived submarine-cable crosswalk is useful prior work, but it lacks an
    external geography citation for each composition and therefore remains a
    reviewer candidate only.
    """

    valid_gid_country = canonical.set_index("adm1_gid")["survey_iso2"].astype(str).to_dict()

    def index_rows(frame: pd.DataFrame) -> dict[tuple[str, str, str], list[dict[str, object]]]:
        indexed: defaultdict[tuple[str, str, str], list[dict[str, object]]] = defaultdict(list)
        for record in frame.to_dict(orient="records"):
            key = (
                str(record["iso2"]).strip(),
                str(record["gallup_region_variable"]).strip(),
                canonical_code(record["gallup_region_code"]),
            )
            indexed[key].append(record)
        return dict(indexed)

    authoritative_columns = [
        "iso2", "year_from", "year_to", "gallup_region_variable",
        "gallup_region_code", "gallup_region_label", "adm1_gid",
        "population_overlap_2025", "boundary_definition", "boundary_vintage",
        "source_id", "source_title", "source_url", "source_access_date",
        "mapping_method", "notes",
    ]
    authoritative = pd.read_csv(
        AUTHORITATIVE_COMPOSITIONS, keep_default_na=False, dtype=str
    )
    missing = set(authoritative_columns).difference(authoritative.columns)
    if missing:
        raise RuntimeError(
            f"Authoritative-composition file lacks columns: {sorted(missing)}"
        )
    authoritative = authoritative[authoritative_columns].copy()
    if len(authoritative):
        duplicate_keys = [
            "iso2", "year_from", "year_to", "gallup_region_variable",
            "gallup_region_code", "adm1_gid",
        ]
        duplicate_rows = authoritative.duplicated(duplicate_keys, keep=False)
        if duplicate_rows.any():
            duplicates = authoritative.loc[duplicate_rows, duplicate_keys]
            raise RuntimeError(
                "Authoritative compositions contain duplicate constituents:\n"
                + duplicates.to_string(index=False)
            )
        for column in ["iso2", "gallup_region_variable", "gallup_region_code", "adm1_gid"]:
            if authoritative[column].str.strip().eq("").any():
                raise RuntimeError(f"Authoritative compositions contain blank {column}")
        required_citation = [
            "boundary_definition", "boundary_vintage", "source_id",
            "source_title", "source_url",
        ]
        if authoritative[required_citation].apply(
            lambda column: column.astype(str).str.strip().eq("")
        ).any().any():
            raise RuntimeError(
                "Every authoritative composition must document its definition, "
                "vintage, source ID, title, and URL"
            )
        bad_gid = authoritative.loc[
            ~authoritative["adm1_gid"].isin(valid_gid_country), "adm1_gid"
        ].unique()
        if len(bad_gid):
            raise RuntimeError(
                f"Authoritative compositions reference non-Meta ADM1 IDs: {sorted(bad_gid)}"
            )
        wrong_country = authoritative.loc[
            authoritative.apply(
                lambda row: valid_gid_country.get(row["adm1_gid"]) != row["iso2"], axis=1
            )
        ]
        if len(wrong_country):
            raise RuntimeError("Authoritative compositions cross survey-country boundaries")

    prior = pd.DataFrame(columns=authoritative_columns)
    if PRIOR_COMPOSITION_CROSSWALK.exists():
        source = pd.read_csv(PRIOR_COMPOSITION_CROSSWALK, keep_default_na=False, dtype=str)
        source = source.loc[
            source["accepted"].astype(str).str.casefold().isin(["true", "1"])
            & source["source_gadm1_id"].str.strip().ne("")
        ].copy()
        prior = pd.DataFrame(
            {
                "iso2": source["iso2"],
                "year_from": "",
                "year_to": "",
                "gallup_region_variable": source["gallup_region_variable"],
                "gallup_region_code": source["gallup_region_code"].map(canonical_code),
                "gallup_region_label": source["gallup_region_label"],
                "adm1_gid": source["source_gadm1_id"],
                "population_overlap_2025": "",
                "boundary_definition": source["historical_region_id"],
                "boundary_vintage": "prior project derivation; verify",
                "source_id": "prior_submarine_geography_crosswalk",
                "source_title": "Archived submarine-cable Gallup geography crosswalk",
                "source_url": "",
                "source_access_date": "",
                "mapping_method": source["mapping_method"],
                "notes": source["notes"],
            }
        )
        prior = prior.loc[prior["adm1_gid"].isin(valid_gid_country)].copy()
        prior = prior.loc[
            prior.apply(
                lambda row: valid_gid_country.get(row["adm1_gid"]) == row["iso2"], axis=1
            )
        ]
    return index_rows(authoritative), index_rows(prior)


def match_observed_codes(
    counts: pd.DataFrame,
    inventory: pd.DataFrame,
    indexes: dict[str, object],
    canonical: gpd.GeoDataFrame,
    authoritative_compositions: dict[tuple[str, str, str], list[dict[str, object]]],
    prior_compositions: dict[tuple[str, str, str], list[dict[str, object]]],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    level_map = inventory.set_index("variable")["region_level"].to_dict()
    available_countries = set(canonical["survey_iso2"].astype(str))
    adm1_aliases = indexes["adm1_aliases"]
    adm1_core_aliases = indexes["adm1_core_aliases"]
    adm1_primary = indexes["adm1_primary"]
    adm1_codes = indexes["adm1_codes"]
    iso_subdivision_aliases = indexes["iso_subdivision_aliases"]
    iso_subdivision_core_aliases = indexes["iso_subdivision_core_aliases"]
    lower_aliases = indexes["lower_aliases"]
    aggregate_aliases = indexes["aggregate_aliases"]
    label_by_gid = indexes["label_by_gid"]

    unique_adm1_fuzzy: defaultdict[str, dict[str, str]] = defaultdict(dict)
    for (iso2, alias), gids in adm1_aliases.items():
        if len(gids) == 1:
            unique_adm1_fuzzy[iso2][alias] = next(iter(gids))
    aggregate_fuzzy: defaultdict[str, dict[str, set[str]]] = defaultdict(dict)
    for (iso2, alias), gids in aggregate_aliases.items():
        if gids:
            aggregate_fuzzy[iso2][alias] = set(gids)

    def best_fuzzy(query: str, choices: list[str]) -> tuple[str, float, float]:
        scored = sorted(
            ((choice, 100.0 * difflib.SequenceMatcher(None, query, choice).ratio()) for choice in choices),
            key=lambda item: item[1],
            reverse=True,
        )
        if not scored:
            return "", np.nan, np.nan
        runner_up = scored[1][1] if len(scored) > 1 else 0.0
        return scored[0][0], scored[0][1], scored[0][1] - runner_up

    v1_map: dict[tuple[str, str, str], pd.Series] = {}
    if V1_CROSSWALK.exists():
        v1 = pd.read_csv(V1_CROSSWALK, keep_default_na=False, dtype=str)
        for _, row in v1.iterrows():
            key = (row["iso2"], row["gallup_region_variable"], canonical_code(row["gallup_region_code"]))
            v1_map[key] = row

    mapping_rows: list[dict[str, object]] = []
    composition_rows: list[dict[str, object]] = []
    for row in counts.itertuples(index=False):
        iso2 = str(row.iso2)
        year = int(row.year_wave)
        variable = str(row.gallup_region_variable)
        code = canonical_code(row.gallup_region_code)
        label = str(row.gallup_region_label).strip()
        normalized = normalize_name(label)
        core = normalize_admin_core(label)
        method = "unmatched"
        review_required = True
        adm1_gids: list[str] = []
        boundary_gid = ""
        boundary_level = ""
        source_id = ""
        source_note = ""
        source_title = ""
        source_url = ""
        boundary_definition = ""
        boundary_vintage = ""
        evidence_tier = "unresolved"
        provisional_usable = False
        overlap_by_gid: dict[str, float] = {}

        def eligible_documented_rows(
            source: dict[tuple[str, str, str], list[dict[str, object]]]
        ) -> list[dict[str, object]]:
            records = source.get((iso2, variable, code), [])
            eligible: list[dict[str, object]] = []
            for record in records:
                label_in_source = str(record.get("gallup_region_label", "")).strip()
                if label_in_source and normalize_name(label_in_source) != normalized:
                    continue
                year_from = pd.to_numeric(record.get("year_from", ""), errors="coerce")
                year_to = pd.to_numeric(record.get("year_to", ""), errors="coerce")
                if np.isfinite(year_from) and year < int(year_from):
                    continue
                if np.isfinite(year_to) and year > int(year_to):
                    continue
                eligible.append(record)
            return eligible

        if iso2 not in available_countries:
            method = "no_compatible_subnational_sci"
            review_required = False
            source_note = "Country has no origin region in the 2026 GADM1-to-country SCI universe."
        elif not label or INVALID_LABEL_RE.search(label):
            method = "invalid_or_unlabeled_code"
            review_required = False
        else:
            exact_primary = adm1_primary.get((iso2, normalized), set())
            exact_alias = adm1_aliases.get((iso2, normalized), set())
            exact_code = adm1_codes.get((iso2, normalized), set())
            iso_subdivision = iso_subdivision_aliases.get((iso2, normalized), set())
            iso_subdivision_core = iso_subdivision_core_aliases.get((iso2, core), set())
            lower = lower_aliases.get((iso2, normalized), set())
            aggregate = aggregate_aliases.get((iso2, normalized), set())
            documented = eligible_documented_rows(authoritative_compositions)
            if documented:
                # A source-cited historical or aggregate definition takes
                # precedence over a coincidentally identical current GADM1
                # name. This prevents, for example, an old survey region that
                # later split from being assigned only to its same-named child.
                adm1_gids = sorted({str(item["adm1_gid"]) for item in documented})
                method = str(documented[0].get("mapping_method", "")).strip() or "authoritative_composition_candidate"
                review_required = True
                source_id = str(documented[0]["source_id"])
                source_title = str(documented[0]["source_title"])
                source_url = str(documented[0]["source_url"])
                boundary_definition = str(documented[0]["boundary_definition"])
                boundary_vintage = str(documented[0]["boundary_vintage"])
                source_note = str(documented[0].get("notes", ""))
                overlap_by_gid = {
                    str(item["adm1_gid"]): float(pd.to_numeric(item.get("population_overlap_2025", np.nan), errors="coerce"))
                    for item in documented
                    if np.isfinite(pd.to_numeric(item.get("population_overlap_2025", np.nan), errors="coerce"))
                }
                evidence_tier = "source_backed_candidate"
                provisional_usable = True
            elif len(exact_code) == 1:
                adm1_gids = sorted(exact_code)
                method = "exact_official_adm1_code"
                review_required = False
                boundary_gid = adm1_gids[0]
                boundary_level = "1"
                source_id = "gadm_4_1_complete_database"
                evidence_tier = "authoritative_exact"
                provisional_usable = True
            elif len(exact_primary) == 1:
                adm1_gids = sorted(exact_primary)
                method = "exact_gadm1_primary_name"
                review_required = False
                boundary_gid = adm1_gids[0]
                boundary_level = "1"
                source_id = "gadm_4_1_complete_database"
                evidence_tier = "authoritative_exact"
                provisional_usable = True
            elif len(exact_alias) == 1:
                adm1_gids = sorted(exact_alias)
                method = "exact_gadm1_official_alternate_name"
                review_required = False
                boundary_gid = adm1_gids[0]
                boundary_level = "1"
                source_id = "gadm_4_1_complete_database"
                evidence_tier = "authoritative_exact"
                provisional_usable = True
            elif len(iso_subdivision) == 1 or len(iso_subdivision_core) == 1:
                candidates = iso_subdivision if len(iso_subdivision) == 1 else iso_subdivision_core
                adm1_gid, iso_code, official_name = next(iter(candidates))
                adm1_gids = [adm1_gid]
                method = "exact_iso_3166_2_subdivision_name_candidate"
                review_required = True
                boundary_gid = adm1_gid
                boundary_level = "1"
                source_id = "iso_3166_2_pycountry_24_6_1_gadm_iso1"
                source_title = "ISO 3166-2 subdivision names linked to GADM 4.1 ISO_1 codes"
                source_url = "https://www.iso.org/obp/ui/#iso:code:3166"
                boundary_definition = (
                    f"Gallup label matches official ISO 3166-2 subdivision "
                    f"{official_name} ({iso_code}); that code uniquely matches GADM ISO_1."
                )
                boundary_vintage = "pycountry 24.6.1 / GADM 4.1"
                source_note = (
                    "Exact official-name and code bridge; Gallup boundary vintage still requires review."
                )
                evidence_tier = "source_backed_candidate"
                provisional_usable = True
            elif lower:
                parents = {item[0] for item in lower}
                lower_units = {(item[1], item[2]) for item in lower}
                if len(parents) == 1 and len(lower_units) == 1:
                    adm1_gids = sorted(parents)
                    boundary_gid, boundary_level_int = next(iter(lower_units))
                    boundary_level = str(boundary_level_int)
                    method = "exact_lower_gadm_unit_to_unique_adm1_parent"
                    review_required = True
                    source_id = "gadm_4_1_complete_database"
                    evidence_tier = "source_backed_candidate"
                    provisional_usable = True
            if not adm1_gids and aggregate:
                adm1_gids = sorted(aggregate)
                method = "exact_gadm_named_aggregate"
                review_required = True
                source_id = "gadm_4_1_complete_database"
                evidence_tier = "source_backed_candidate"
                provisional_usable = True
            if not adm1_gids and core and core != normalized:
                core_matches = adm1_core_aliases.get((iso2, core), set())
                if len(core_matches) == 1:
                    adm1_gids = sorted(core_matches)
                    method = "exact_normalized_administrative_core_candidate"
                    review_required = True
                    boundary_gid = adm1_gids[0]
                    boundary_level = "1"
                    source_id = "gadm_4_1_complete_database"
                    evidence_tier = "source_backed_candidate"
                    provisional_usable = True
            if not adm1_gids:
                old = v1_map.get((iso2, variable, code))
                old_gid = "" if old is None else str(old.get("adm1_gid", "")).strip()
                old_method = "" if old is None else str(old.get("match_method", ""))
                old_label = "" if old is None else str(old.get("gallup_region_label", ""))
                old_exact_still_valid = (
                    old_gid in indexes["country_by_gid"]
                    and indexes["country_by_gid"].get(old_gid) == iso2
                    and normalize_name(old_label) == normalized
                )
                if old_exact_still_valid and old_method == "exact_normalized_name":
                    adm1_gids = [old_gid]
                    method = "verified_version1_exact_normalized_gadm1_name"
                    review_required = False
                    boundary_gid = old_gid
                    boundary_level = "1"
                    source_id = "gallup_gadm1_crosswalk_v1|gadm_4_1_complete_database"
                    source_note = "Revalidated against the same Gallup code/label and current Meta origin universe."
                    evidence_tier = "authoritative_exact"
                    provisional_usable = True
                elif old_gid in indexes["country_by_gid"] and old_method in {
                    "high_confidence_fuzzy_name", "documented_geographic_override"
                }:
                    adm1_gids = [old_gid]
                    method = (
                        "version1_fuzzy_candidate"
                        if old_method == "high_confidence_fuzzy_name"
                        else "version1_documented_override_candidate"
                    )
                    review_required = True
                    boundary_gid = old_gid
                    boundary_level = "1"
                    source_id = "gallup_gadm1_crosswalk_v1"
                    source_note = "Candidate only; requires an authoritative source and two-person review."
                    evidence_tier = "prior_project_candidate"
                elif (
                    old_gid in indexes["country_by_gid"]
                    and old is not None
                    and str(old.get("match_status", "")) == "low_confidence_fuzzy"
                    and float(pd.to_numeric(old.get("match_score", np.nan), errors="coerce")) >= 85
                    and core not in GENERIC_BROAD_NAMES
                    and normalized not in GENERIC_BROAD_NAMES
                ):
                    adm1_gids = [old_gid]
                    method = "version1_high_score_fuzzy_candidate"
                    review_required = True
                    boundary_gid = old_gid
                    boundary_level = "1"
                    source_id = "gallup_gadm1_crosswalk_v1"
                    source_note = "Spelling/translation candidate only; two reviewers must verify against an official source."
                    evidence_tier = "fuzzy_candidate"
            prior_documented = eligible_documented_rows(prior_compositions)
            if not adm1_gids and prior_documented:
                adm1_gids = sorted({str(item["adm1_gid"]) for item in prior_documented})
                method = "prior_project_historical_composition_candidate"
                review_required = True
                source_id = "prior_submarine_geography_crosswalk"
                source_title = "Archived submarine-cable Gallup geography crosswalk"
                boundary_definition = str(prior_documented[0].get("boundary_definition", ""))
                boundary_vintage = str(prior_documented[0].get("boundary_vintage", ""))
                source_note = (
                    "Prior project composition only. Supply an authoritative geography source "
                    "and obtain two approvals before final use."
                )
                evidence_tier = "prior_project_candidate"
            if not adm1_gids and normalized not in GENERIC_BROAD_NAMES:
                aggregate_choice, aggregate_score, aggregate_margin = best_fuzzy(
                    normalized, list(aggregate_fuzzy.get(iso2, {}))
                )
                if aggregate_score >= 82 and aggregate_margin >= 5:
                    adm1_gids = sorted(aggregate_fuzzy[iso2][aggregate_choice])
                    method = "fuzzy_gadm_named_aggregate_candidate"
                    review_required = True
                    source_id = "gadm_4_1_complete_database"
                    source_note = (
                        f"Reviewer candidate: label similarity {aggregate_score:.1f}, "
                        f"margin {aggregate_margin:.1f}; never auto-accepted."
                    )
                    evidence_tier = "fuzzy_candidate"
            if not adm1_gids and normalized not in GENERIC_BROAD_NAMES:
                adm1_choice, adm1_score, adm1_margin = best_fuzzy(
                    normalized, list(unique_adm1_fuzzy.get(iso2, {}))
                )
                if adm1_score >= 88 and adm1_margin >= 5:
                    adm1_gids = [unique_adm1_fuzzy[iso2][adm1_choice]]
                    method = "fuzzy_gadm1_candidate"
                    review_required = True
                    boundary_gid = adm1_gids[0]
                    boundary_level = "1"
                    source_id = "gadm_4_1_complete_database"
                    source_note = (
                        f"Reviewer candidate: label similarity {adm1_score:.1f}, "
                        f"margin {adm1_margin:.1f}; never auto-accepted."
                    )
                    evidence_tier = "fuzzy_candidate"

        if source_id == "gadm_4_1_complete_database":
            source_title = source_title or "GADM 4.1 database"
            source_url = source_url or "https://gadm.org/"
            boundary_definition = boundary_definition or (
                f"Gallup label resolved to {method.replace('_', ' ')} in GADM 4.1"
            )
            boundary_vintage = boundary_vintage or "GADM 4.1"
        elif source_id == "gallup_gadm1_crosswalk_v1|gadm_4_1_complete_database":
            source_title = source_title or "Archived v1 crosswalk revalidated against GADM 4.1"
            source_url = source_url or "https://gadm.org/"
            boundary_definition = boundary_definition or "Unique normalized GADM1 name"
            boundary_vintage = boundary_vintage or "GADM 4.1"

        source_complete_for_final = bool(
            adm1_gids
            and source_id
            and source_title
            and source_url
            and boundary_definition
            and boundary_vintage
        )
        mapping_id = stable_id(iso2, year, variable, code)
        review_id = stable_id(
            iso2, variable, code, normalized, method, "|".join(adm1_gids),
            source_id, boundary_vintage,
        )
        if adm1_gids:
            candidate_status = "pending_two_person_review" if review_required else "auto_approved_exact"
            for gid in adm1_gids:
                composition_rows.append(
                    {
                        "mapping_id": mapping_id,
                        "iso2": iso2,
                        "year_wave": year,
                        "gallup_region_variable": variable,
                        "gallup_region_code": code,
                        "gallup_region_label": label,
                        "adm1_gid": gid,
                        "adm1_name": label_by_gid.get(gid, ""),
                        "boundary_gid": boundary_gid,
                        "boundary_level": boundary_level,
                        "match_method": method,
                        "source_id": source_id,
                        "population_overlap_2025": overlap_by_gid.get(gid, np.nan),
                        "boundary_definition": boundary_definition,
                        "boundary_vintage": boundary_vintage,
                        "source_title": source_title,
                        "source_url": source_url,
                    }
                )
        else:
            candidate_status = "unmatched"
        mapping_rows.append(
            {
                "mapping_id": mapping_id,
                "review_id": review_id,
                "iso2": iso2,
                "year_wave": year,
                "gallup_region_variable": variable,
                "region_level": level_map.get(variable, np.nan),
                "gallup_region_code": code,
                "gallup_region_label": label,
                "normalized_gallup_region_label": normalized,
                "respondents": int(row.respondents),
                "weight_sum": float(row.weight_sum),
                "observed_in_wave": bool(int(row.respondents) > 0),
                "declared_in_codebook": bool(getattr(row, "declared_in_codebook", False)),
                "candidate_status": candidate_status,
                "match_method": method,
                "constituent_adm1_count": len(adm1_gids),
                "constituent_adm1_ids": "|".join(adm1_gids),
                "boundary_gid": boundary_gid,
                "boundary_level": boundary_level,
                "review_required": review_required and bool(adm1_gids),
                "source_id": source_id,
                "source_note": source_note,
                "source_title": source_title,
                "source_url": source_url,
                "boundary_definition": boundary_definition,
                "boundary_vintage": boundary_vintage,
                "evidence_tier": evidence_tier,
                "provisional_usable": provisional_usable and bool(adm1_gids),
                "source_complete_for_final": source_complete_for_final,
            }
        )

    mappings = pd.DataFrame(mapping_rows)
    compositions = pd.DataFrame(composition_rows)

    # Fuzzy scores are reviewer aids only. They never change match status.
    choices_by_country: dict[str, dict[str, str]] = {}
    for iso2, part in canonical.groupby("survey_iso2"):
        choices: dict[str, str] = {}
        for item in part.itertuples(index=False):
            for alias in [item.NAME_1, *split_aliases(item.VARNAME_1), *split_aliases(item.NL_NAME_1)]:
                normalized = normalize_name(alias)
                if normalized:
                    choices[normalized] = str(item.adm1_gid)
        choices_by_country[str(iso2)] = choices
    mappings["fuzzy_candidate_adm1_gid"] = ""
    mappings["fuzzy_candidate_name"] = ""
    mappings["fuzzy_score"] = np.nan
    for index, row in mappings.loc[mappings["candidate_status"].eq("unmatched")].iterrows():
        choices = choices_by_country.get(str(row["iso2"]), {})
        query = str(row["normalized_gallup_region_label"])
        if not choices or not query:
            continue
        result = difflib.get_close_matches(query, list(choices), n=1, cutoff=0.0)
        if not result:
            continue
        candidate_name = result[0]
        score = 100.0 * difflib.SequenceMatcher(None, query, candidate_name).ratio()
        mappings.loc[index, "fuzzy_candidate_adm1_gid"] = choices[candidate_name]
        mappings.loc[index, "fuzzy_candidate_name"] = candidate_name
        mappings.loc[index, "fuzzy_score"] = float(score)
    return mappings, compositions


CANDIDATES = RESTRICTED_V2 / "all_observed_region_code_candidates.csv"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--chunksize", type=int, default=50_000)
    parser.add_argument(
        "--rebuild-gallup-counts",
        action="store_true",
        help="Re-stream the Gallup archive instead of reusing cached region-code counts.",
    )
    args = parser.parse_args()
    ensure_v2_directories()
    inventory, labels = parse_gallup_dictionary()
    counts, _, countries = stream_region_counts(
        inventory, labels, args.chunksize, args.rebuild_gallup_counts
    )
    origins, _ = read_meta_origins()
    attributes = read_gadm_attributes(set(origins["adm1_gid"]))
    canonical = build_canonical_adm1(origins, attributes, countries)
    indexes = build_alias_indexes(attributes, canonical)
    authoritative_compositions, prior_compositions = load_documented_compositions(canonical)
    mappings, _ = match_observed_codes(
        counts,
        inventory,
        indexes,
        canonical,
        authoritative_compositions,
        prior_compositions,
    )
    mappings.to_csv(CANDIDATES, index=False)
    print(f"Wrote {len(mappings):,} region-code match candidates to {CANDIDATES}")


if __name__ == "__main__":
    main()
