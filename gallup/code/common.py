#!/usr/bin/env python3
"""Shared paths and helpers for the pipeline."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


PACKAGE = Path(__file__).resolve().parents[1]
RAW = PACKAGE / "data" / "raw"
RESTRICTED = PACKAGE / "data" / "restricted"
DERIVED = PACKAGE / "data" / "derived"
OUTPUTS = PACKAGE / "outputs"
MODELS = OUTPUTS / "models"
YEARS = (2016, 2019, 2022, 2023)


def ensure_directories() -> None:
    for path in [RAW, RESTRICTED / "derived", DERIVED, OUTPUTS]:
        path.mkdir(parents=True, exist_ok=True)


def weighted_mean(values: pd.Series, weights: pd.Series) -> float:
    x = pd.to_numeric(values, errors="coerce").to_numpy(float)
    w = pd.to_numeric(weights, errors="coerce").to_numpy(float)
    valid = np.isfinite(x) & np.isfinite(w) & (w > 0)
    if not valid.any():
        return float("nan")
    return float(np.average(x[valid], weights=w[valid]))


def weighted_sd(values: pd.Series, weights: pd.Series) -> float:
    mean = weighted_mean(values, weights)
    x = pd.to_numeric(values, errors="coerce").to_numpy(float)
    w = pd.to_numeric(weights, errors="coerce").to_numpy(float)
    valid = np.isfinite(x) & np.isfinite(w) & (w > 0)
    return float(np.sqrt(np.average(np.square(x[valid] - mean), weights=w[valid])))


def require_columns(frame: pd.DataFrame, columns: Iterable[str], label: str) -> None:
    missing = sorted(set(columns).difference(frame.columns))
    if missing:
        raise RuntimeError(f"{label} lacks required columns: {missing}")
