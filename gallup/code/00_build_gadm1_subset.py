#!/usr/bin/env python3
"""Build the GADM 4.1 ADM1 layer used by the pipeline from the full GADM GeoPackage.

GADM does not permit redistribution of its boundaries, so the repository ships
only the list of ADM1 identifiers (data/raw/geography/gadm1_subset_ids.csv).
This script dissolves gadm_410.gpkg to those ADM1 units. It takes about 10
minutes and is skipped if the output already exists.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import pyogrio

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import RAW

GADM = RAW / "geography" / "gadm_410_complete.gpkg"
IDS = RAW / "geography" / "gadm1_subset_ids.csv"
OUTPUT = RAW / "geography" / "gadm1_metadata.gpkg"


def main() -> None:
    if OUTPUT.exists():
        print(f"{OUTPUT} already exists; skipping.")
        return
    ids = pd.read_csv(IDS, dtype=str)["GID_1"].tolist()
    quoted = ",".join("'" + value.replace("'", "''") + "'" for value in ids)
    sql = (
        "SELECT GID_1, GID_0, COUNTRY, NAME_1, VARNAME_1, TYPE_1, ENGTYPE_1, ISO_1, "
        f"ST_Union(geom) AS geom FROM gadm_410 WHERE GID_1 IN ({quoted}) GROUP BY GID_1"
    )
    layer = pyogrio.read_dataframe(GADM, sql=sql, sql_dialect="SQLITE")
    if len(layer) != len(ids):
        raise RuntimeError(f"Expected {len(ids):,} ADM1 units; built {len(layer):,}")
    layer.to_file(OUTPUT, layer="gadm1_metadata", driver="GPKG")
    print(f"Wrote {len(layer):,} ADM1 units to {OUTPUT}")


if __name__ == "__main__":
    main()
