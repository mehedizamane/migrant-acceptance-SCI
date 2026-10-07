#!/usr/bin/env python3
"""Run the full Gallup/SCI pipeline in order.

Requires the input files described in README.md under data/. Each stage can also
be run on its own; they communicate only through files.
"""

from __future__ import annotations

from pathlib import Path
import subprocess
import sys


PACKAGE = Path(__file__).resolve().parent
CODE = PACKAGE / "code"
RESTRICTED = PACKAGE / "data" / "restricted" / "derived"
OUTPUTS = PACKAGE / "outputs"

STAGES = [
    [sys.executable, CODE / "00_build_gadm1_subset.py"],
    [sys.executable, CODE / "01_build_sci_measures.py"],
    [sys.executable, CODE / "02_build_controls.py"],
    [sys.executable, CODE / "03_build_geography_adjusted_sci.py"],
    [sys.executable, CODE / "04_build_gallup_analysis_data.py"],
    [sys.executable, CODE / "05_prepare_country_model_input.py"],
    ["Rscript", CODE / "06_fit_country_models.R",
     RESTRICTED / "gallup_model_input.csv", OUTPUTS / "models"],
    [sys.executable, CODE / "07_seed_region_compositions.py"],
    [sys.executable, CODE / "08_match_gallup_regions.py"],
    [sys.executable, CODE / "09_link_gallup_regions_to_adm1.py"],
    ["Rscript", CODE / "10_fit_regional_models.R",
     RESTRICTED / "geography_direct_adm1" / "gallup_regional_model_input_direct_adm1.csv",
     OUTPUTS / "regional_models"],
    [sys.executable, CODE / "11_make_figures_and_tables.py"],
]


def main() -> None:
    for command in STAGES:
        command = [str(part) for part in command]
        print(f"\n=== {Path(command[1]).name} ===", flush=True)
        subprocess.run(command, cwd=PACKAGE, check=True)


if __name__ == "__main__":
    main()
