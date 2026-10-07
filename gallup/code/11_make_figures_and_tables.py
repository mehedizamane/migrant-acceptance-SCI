#!/usr/bin/env python3
"""Build the paper's Gallup/SCI figures and tables from the fitted models.

Figures 1, 2, 4, and 5 are written to outputs/figures/. Every Gallup number in
SI Tables 1 and 8-17, plus values quoted only in the text, is written as CSV to
outputs/tables/. Figure 3 and SI Tables 4-7 and 18-27 report the European
Social Survey and Eurobarometer analyses, which are not part of this package.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import DERIVED, MODELS, OUTPUTS

REGIONAL_MODELS = OUTPUTS / "regional_models"
LINKAGE = OUTPUTS / "regional_linkage"
FIGURES = OUTPUTS / "figures"
TABLES = OUTPUTS / "tables"

LMM = "weighted_lmm_ml"
OLS = "weighted_ols_country_cr1"
COUNTRY_FE = "weighted_ols_country_fe_cr1"
SCI = "global_sci_2026_all_origins_z"
WAVES = ["2016", "2019", "2022", "2023"]
POOLED = "2016+2019+2022+2023"

BG = "#ffffff"
TEXT = "#181c1b"
MUTED = "#6d7774"
RULE = "#c7cecb"
TEAL = "#087f76"
TEAL_DARK = "#174f4b"
CORAL = "#dc5d4e"


# ---------------------------------------------------------------------------
# Model results
# ---------------------------------------------------------------------------

def load_models(directory: Path) -> pd.DataFrame:
    """Coefficients joined to their model registry rows."""
    coefficients = pd.read_csv(directory / "model_coefficients.csv")
    registry = pd.read_csv(directory / "model_registry.csv")
    keep = [c for c in ["model_id", "estimator", "year_wave", "n_respondents", "n_countries",
                        "n_adm1_regions", "n_exposure_clusters", "country_icc",
                        "exposure_region_icc"] if c in registry.columns]
    return coefficients.merge(registry[keep], on=["model_id", "estimator"], how="left",
                              validate="many_to_one")


COUNTRY = load_models(MODELS)
REGIONAL = load_models(REGIONAL_MODELS)


def result(model_id: str, term: str, estimator: str = LMM, frame: pd.DataFrame = COUNTRY) -> pd.Series:
    rows = frame[(frame["model_id"] == model_id) & (frame["estimator"] == estimator) & (frame["term"] == term)]
    if len(rows) != 1:
        raise KeyError(f"Expected one result for {model_id} / {estimator} / {term}; found {len(rows)}")
    return rows.iloc[0]


def interaction_term(model_id: str, estimator: str, *parts: str) -> str:
    """Interaction coefficient name; lme4 and fixest order the operands differently."""
    terms = COUNTRY.loc[(COUNTRY["model_id"] == model_id) & (COUNTRY["estimator"] == estimator), "term"]
    for term in terms:
        if ":" in term and set(term.split(":")) == set(parts):
            return term
    raise KeyError(f"No {':'.join(parts)} interaction in {model_id} / {estimator}")


def row(label: dict, r: pd.Series, *, n: bool = True, icc: bool = False) -> dict:
    out = {**label, "estimate": r["estimate"], "ci_low": r["conf_low"], "ci_high": r["conf_high"],
           "p_value": r["p_value"]}
    if n:
        out["respondents"] = int(r["n_respondents"])
        out["countries"] = int(r["n_countries"])
    if icc:
        out["country_icc"] = r["country_icc"]
    return out


def holm_adjust(p_values: pd.Series) -> np.ndarray:
    p = np.asarray(p_values, dtype=float)
    order = np.argsort(p)
    ranked = p[order]
    m = len(ranked)
    adjusted_sorted = np.maximum.accumulate((m - np.arange(m)) * ranked)
    adjusted_sorted = np.minimum(adjusted_sorted, 1.0)
    adjusted = np.empty(m, dtype=float)
    adjusted[order] = adjusted_sorted
    return adjusted


def write(frame: pd.DataFrame, name: str) -> pd.DataFrame:
    frame.to_csv(TABLES / f"{name}.csv", index=False)
    return frame


# ---------------------------------------------------------------------------
# Tables
# ---------------------------------------------------------------------------

def table_1() -> None:
    """Gallup-to-GADM1 linkage coverage, accepted match classes, and exclusions."""
    coverage = pd.read_csv(LINKAGE / "linkage_coverage.csv")
    write(coverage[["year_wave", "all_wave_respondents", "eligible_geography_respondents",
                    "linked_respondents", "linked_pct_all", "linked_pct_eligible",
                    "model_respondents"]], "si_table_01a_linkage_coverage")
    classes = pd.read_csv(LINKAGE / "match_classes.csv")
    display_class = {
        "Exact official code or primary GADM1 name": "Exact official code or primary GADM1 name",
        "Unique exact name after removing administrative-type words": "Exact name after administrative-type normalization",
        "Exact official alternate or previously verified normalized name": "Exact alternate, verified normalized, or ISO subdivision name",
        "Unique exact ISO 3166-2 subdivision name": "Exact alternate, verified normalized, or ISO subdivision name",
        "Exact lower-level GADM unit nested in one GADM1 parent": "Exact lower-level unit nested in one GADM1 parent",
    }
    classes["match_class"] = classes["direct_match_class"].map(display_class)
    if classes["match_class"].isna().any():
        raise RuntimeError("Unmapped direct match class")
    write(classes.groupby("match_class", as_index=False)
          .agg(country_wave_codes=("mapped_country_wave_codes", "sum"),
               respondents=("candidate_respondents", "sum"))
          .sort_values("respondents", ascending=False), "si_table_01b_match_classes")
    exclusions = pd.read_csv(LINKAGE / "exclusion_reasons.csv")
    display_reason = {
        "No direct GADM1 name or code match": "No deterministic ADM1 match",
        "Fuzzy candidate withheld pending documented manual review": "No deterministic ADM1 match",
        "Legacy mapping lacks sufficient source evidence for automatic reuse": "No deterministic ADM1 match",
        "Country absent from the compatible Meta GADM1 SCI origin universe": "Country not covered by ADM1 SCI data",
        "Gallup region aggregates multiple GADM1 units": "Broader, aggregate, or historical region without an authoritative composition",
        "Aggregate or historical concordance requires authoritative composition": "Broader, aggregate, or historical region without an authoritative composition",
    }
    exclusions["reason"] = exclusions["rejection_reason"].map(display_reason)
    if exclusions["reason"].isna().any():
        raise RuntimeError("Unmapped exclusion reason")
    write(exclusions.groupby("reason", as_index=False)
          .agg(country_wave_codes=("country_wave_codes", "sum"), respondents=("respondents", "sum"))
          .sort_values("respondents", ascending=False), "si_table_01c_exclusions")


def table_8() -> pd.DataFrame:
    """Primary country GlobalSCI models (also Figure 1)."""
    rows = [row({"wave": w}, result(f"country_primary_2026_{w}", SCI), icc=True) for w in WAVES]
    rows.append(row({"wave": "Pooled, four waves"}, result("country_pooled_primary_2026", SCI), icc=True))
    return write(pd.DataFrame(rows), "si_table_08_primary_country_models")


def table_9() -> None:
    """2026 and 2020 SCI releases on identical respondents and countries."""
    rows = [row({"sci_release": rel, "wave": w}, result(f"country_release_matched_{rel}_{w}", f"global_sci_{rel}_z"))
            for rel in (2026, 2020) for w in WAVES]
    write(pd.DataFrame(rows), "si_table_09_release_comparison")


def table_10() -> None:
    """Weighted OLS with country-clustered CR1 errors for the primary models."""
    rows = [row({"wave": w}, result(f"country_primary_2026_{w}", SCI, OLS)) for w in WAVES]
    rows.append(row({"wave": "Pooled"}, result("country_pooled_primary_2026", SCI, OLS)))
    write(pd.DataFrame(rows), "si_table_10_ols_sensitivity")


def table_11() -> None:
    """MAI components with Holm adjustment, strict MAI, and component interactions."""
    outcomes = {
        "immigrants_living_in_country_good": "Immigrants living in country",
        "immigrant_neighbor_good": "Immigrant as neighbor",
        "immigrant_marry_close_relative_good": "Immigrant marrying close relative",
    }
    rows = [row({"wave": w, "outcome": label}, result(f"country_component_2026_{w}_{slug}", SCI))
            for slug, label in outcomes.items() for w in WAVES]
    components = pd.DataFrame(rows)
    # Holm adjustment across the three components within each wave.
    components["holm_p"] = components.groupby("wave")["p_value"].transform(holm_adjust)
    write(components, "si_table_11a_mai_components")

    rows = [row({"wave": w}, result(f"country_strict_mai_2026_{w}", SCI)) for w in WAVES]
    write(pd.DataFrame(rows), "si_table_11b_strict_mai")

    rows = []
    for slug, label in outcomes.items():
        model_id = f"country_contact_component_2026_{slug}"
        term = interaction_term(model_id, LMM, SCI, "stranger_contact_z")
        rows.append(row({"outcome": label}, result(model_id, term)))
    contact = pd.DataFrame(rows)
    # Holm adjustment across the three 2022 component interactions.
    contact["holm_p"] = holm_adjust(contact["p_value"])
    write(contact, "si_table_11c_contact_component_interactions")


def table_12() -> None:
    """Geography-adjusted country GlobalSCI and regional deviation."""
    rows = []
    for rel in (2026, 2020):
        for w in [*WAVES, "pooled"]:
            model_id = (f"country_geography_adjusted_pooled_{rel}" if w == "pooled"
                        else f"country_geography_adjusted_{rel}_{w}")
            r = result(model_id, f"global_sci_geography_adjusted_{rel}_z")
            rows.append({**row({"measure": f"{rel} release: geography-adjusted GlobalSCI",
                                "wave": "2016-2023" if w == "pooled" else w}, r),
                         "adm1_regions": int(r["n_adm1_regions"])})
    for w in [*WAVES, "pooled"]:
        r = result(f"regional_geography_adjusted_within_{w}", "regional_sci_geography_adjusted_within_z",
                   frame=REGIONAL)
        rows.append({**row({"measure": "2026 release: geography-adjusted regional SCI deviation",
                            "wave": "2016-2023" if w == "pooled" else w}, r),
                     "adm1_regions": int(r["n_exposure_clusters"])})
    write(pd.DataFrame(rows), "si_table_12_geography_adjusted_sci")


def table_13() -> pd.DataFrame:
    """Country-average and within-country regional SCI (also Figure 2C)."""
    components = {
        "survey_region_sci_between_total_sd": "Country average",
        "survey_region_sci_within_total_sd": "Regional deviation",
    }
    rows = []
    for w in [*WAVES, "pooled"]:
        for term, label in components.items():
            r = result(f"regional_within_between_{w}", term, frame=REGIONAL)
            rows.append({**row({"wave": "Pooled" if w == "pooled" else w, "component": label}, r),
                         "adm1_regions": int(r["n_exposure_clusters"]),
                         "country_icc": r["country_icc"], "adm1_icc": r["exposure_region_icc"]})
    return write(pd.DataFrame(rows), "si_table_13_regional_contextual")


def table_14() -> None:
    """Pooled regional model restricted to exact code- or name-matched regions."""
    rows = [row({"component": label}, result("regional_within_between_exact_only", term, frame=REGIONAL))
            for term, label in [("survey_region_sci_between_total_sd", "Country-average regional SCI"),
                                ("survey_region_sci_within_total_sd", "Regional deviation from country average")]]
    write(pd.DataFrame(rows), "si_table_14_exact_only_regional")


def table_15() -> None:
    """Directly linked versus unlinked respondents (descriptive)."""
    write(pd.read_csv(LINKAGE / "linked_vs_unlinked.csv"), "si_table_15_linked_vs_unlinked")


def table_16() -> pd.DataFrame:
    """Equal-wave GlobalSCI slope by UN-M49-derived region (also Figure 2D)."""
    contrasts = pd.read_csv(MODELS / "model_linear_contrasts.csv")
    contrasts = contrasts[(contrasts["estimator"] == LMM)
                          & (contrasts["contrast"] == "Equal-wave average SCI slope")]
    registry = pd.read_csv(MODELS / "model_registry.csv")
    countries = registry[registry["estimator"] == LMM].set_index("model_id")["n_countries"]
    table = pd.DataFrame({
        "un_region": contrasts["region"],
        "estimate": contrasts["estimate"],
        "ci_low": contrasts["conf_low"],
        "ci_high": contrasts["conf_high"],
        "p_value": contrasts["p_value"],
        "countries": contrasts["model_id"].map(countries).astype(int),
    }).sort_values("un_region")
    table["inference"] = np.where(table["countries"] < 10, "Descriptive only (<10 countries)", "Reported")
    return write(table, "si_table_16_un_regions")


def table_17() -> None:
    """GlobalSCI, stranger contact, and their interaction in 2022 (also Figure 5)."""
    prefix = "country_contact_2026_"
    specs = [
        ("Contact only", "contact_only", [("Contact", "stranger_contact_z")]),
        ("Additive SCI + contact", "additive", [("GlobalSCI", SCI), ("Contact", "stranger_contact_z")]),
        ("SCI by contact", "interaction",
         [("GlobalSCI", SCI), ("Contact", "stranger_contact_z"), ("GlobalSCI x contact", (SCI, "stranger_contact_z"))]),
        ("Contextual contact", "within_between_contact",
         [("GlobalSCI", SCI), ("Within-country contact", "stranger_contact_within_z"),
          ("Country-average contact", "stranger_contact_country_mean_z"),
          ("GlobalSCI x within-country contact", (SCI, "stranger_contact_within_z"))]),
    ]
    rows = []
    for spec_label, spec, terms in specs:
        for estimator, estimator_label in [(LMM, "Mixed model"), (OLS, "OLS, CR1")]:
            model_id = prefix + spec
            for term_label, term in terms:
                name = interaction_term(model_id, estimator, *term) if isinstance(term, tuple) else term
                r = result(model_id, name, estimator)
                rows.append({**row({"specification": spec_label, "estimator": estimator_label, "term": term_label}, r),
                             "country_icc": r["country_icc"] if estimator == LMM else np.nan})
    model_id = prefix + "interaction_country_fe"
    for term_label, term in [("Within-country contact", "stranger_contact_within_z"),
                             ("GlobalSCI x within-country contact", (SCI, "stranger_contact_within_z"))]:
        name = interaction_term(model_id, COUNTRY_FE, *term) if isinstance(term, tuple) else term
        rows.append(row({"specification": "Country fixed effects", "estimator": "OLS, CR1",
                         "term": term_label}, result(model_id, name, COUNTRY_FE)))
    write(pd.DataFrame(rows), "si_table_17_stranger_contact")


def in_text_values(primary: pd.DataFrame, regional: pd.DataFrame, concordance_r: float,
                   concordance_n: int) -> None:
    """Numbers quoted in the text that do not appear in an SI table."""
    values = []
    for w in WAVES:
        r = result(f"country_available_2026_{w}", SCI)
        values.append((f"All available countries, {w}: GlobalSCI estimate", r["estimate"]))
        values.append((f"All available countries, {w}: countries", r["n_countries"]))
    for rel in (2026, 2020):
        values.append((f"Matched-release pooled estimate, {rel} release",
                       result(f"country_release_matched_pooled_{rel}", f"global_sci_{rel}_z")["estimate"]))
    stability = pd.read_csv(DERIVED / "country_sci_release_stability.csv")
    country = stability[stability["level"] == "Country GlobalSCI"].iloc[0]
    values.append(("2020 vs 2026 country GlobalSCI: Pearson r", country["pearson_r"]))
    values.append(("2020 vs 2026 country GlobalSCI: Spearman rho", country["spearman_rho"]))
    values.append(("Population-aggregated RegionalSCI vs country GlobalSCI: Pearson r", concordance_r))
    values.append(("Population-aggregated RegionalSCI vs country GlobalSCI: countries", concordance_n))
    variance = pd.read_csv(DERIVED / "regional_sci_variance_decomposition.csv").set_index("component")
    values.append(("RegionalSCI variance between countries (%)", 100 * variance.loc["Between countries", "share_of_total_variance"]))
    values.append(("RegionalSCI variance within countries (%)", 100 * variance.loc["Within countries", "share_of_total_variance"]))
    values.append(("Total RegionalSCI SD (log-SCI units)", variance.loc["Total", "standard_deviation"]))
    pooled = regional[regional["wave"] == "Pooled"].set_index("component")["estimate"]
    values.append(("Pooled between-country / within-country regional estimate",
                   pooled["Country average"] / pooled["Regional deviation"]))
    # Halving every pairwise SCI lowers each log by ln 2, so GlobalSCI falls by
    # ln 2 / SD(GlobalSCI) standard deviations (SD over all 178 2026 origins).
    scaling = pd.read_csv(DERIVED / "country_sci_release_scaling.csv")
    sd = float(scaling.loc[scaling["standardized_variable"] == SCI, "reference_sd_log_sci"].iloc[0])
    shift = np.log(2) / sd
    pooled_primary = primary.loc[primary["wave"] == "Pooled, four waves", "estimate"].iloc[0]
    values.append(("Halving cross-border SCI: change in GlobalSCI (SD)", shift))
    values.append(("Halving cross-border SCI: implied MAI points (x pooled estimate)", shift * pooled_primary))
    write(pd.DataFrame(values, columns=["quantity", "value"]), "in_text_values")


# ---------------------------------------------------------------------------
# Figures (plotting code from the manuscript build)
# ---------------------------------------------------------------------------

def set_plot_style() -> None:
    plt.rcParams.update(
        {
            "figure.facecolor": BG,
            "axes.facecolor": BG,
            "savefig.facecolor": BG,
            "font.family": "serif",
            "font.serif": ["Palatino Linotype", "Palatino", "Book Antiqua", "DejaVu Serif"],
            "axes.labelcolor": TEXT,
            "axes.titlecolor": TEXT,
            "text.color": TEXT,
            "xtick.color": MUTED,
            "ytick.color": TEXT,
            "axes.edgecolor": RULE,
            "axes.linewidth": 0.8,
            "font.size": 10,
            "axes.titlesize": 11,
            "axes.titleweight": "bold",
            "axes.labelsize": 9,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "legend.frameon": False,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def clean_axis(ax: plt.Axes, zero_line: bool = False) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color(RULE)
    ax.spines["bottom"].set_color(RULE)
    ax.grid(False)
    if zero_line:
        ax.axvline(0, color=RULE, linewidth=0.9, zorder=0)


def save_figure(fig: plt.Figure, stem: str, *, tight: bool = True) -> None:
    bbox_inches = "tight" if tight else None
    fig.savefig(FIGURES / f"{stem}.pdf", bbox_inches=bbox_inches, facecolor="white")
    fig.savefig(FIGURES / f"{stem}.png", dpi=320, bbox_inches=bbox_inches, facecolor="white")
    plt.close(fig)


def add_panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(0, 1.025, label, transform=ax.transAxes, ha="left", va="bottom",
            fontsize=11, fontweight="bold", color=TEXT)


def figure_1(primary: pd.DataFrame) -> None:
    """Country GlobalSCI coefficient in each Gallup wave."""
    main = primary[primary["wave"].isin(WAVES)].rename(
        columns={"ci_low": "low", "ci_high": "high", "respondents": "Respondents", "wave": "Wave"})
    fig, ax = plt.subplots(figsize=(7.0, 3.8))
    waves = main["Wave"].astype(str).tolist()
    wave_labels = [f"{row['Wave']}  (n = {int(row['Respondents']):,})" for _, row in main.iterrows()]
    y = np.arange(len(waves))[::-1]
    ax.errorbar(
        main["estimate"], y,
        xerr=[main["estimate"] - main["low"], main["high"] - main["estimate"]],
        fmt="o", color=TEAL, ecolor=TEAL, markersize=8, linewidth=2.2, capsize=0,
    )
    for yi, (_, row) in zip(y, main.iterrows()):
        ax.text(row["high"] + 0.035, yi, f"{row['estimate']:.2f}  [{row['low']:.2f}, {row['high']:.2f}]",
                va="center", fontsize=12, color=TEAL_DARK)
    ax.set_yticks(y, wave_labels, fontsize=12.5)
    ax.set_xlabel("MAI points per 1 SD higher GlobalSCI", fontsize=13)
    ax.set_xticks(np.arange(0, 1.01, 0.2))
    ax.tick_params(axis="x", labelsize=11.5)
    ax.set_xlim(0, 1.52)
    ax.set_ylim(-0.32, len(waves) - 0.68)
    ax.spines["left"].set_bounds(y.min(), y.max())
    ax.spines["bottom"].set_bounds(0, 1.0)
    clean_axis(ax, zero_line=True)
    fig.subplots_adjust(left=0.27, right=0.99, top=0.96, bottom=0.21)
    save_figure(fig, "figure_1_country_association", tight=False)


def figure_4() -> None:
    """GlobalSCI before and after full national adjustment, identical samples."""
    rows = []
    for w in [*WAVES, "pooled"]:
        for spec, label in [("matched_baseline", "Matched baseline"), ("fully_adjusted", "Fully adjusted")]:
            model_id = (f"country_national_adjustment_pooled_2026_{spec}" if w == "pooled"
                        else f"country_national_adjustment_2026_{w}_{spec}")
            rows.append(row({"Wave": "Pooled, four waves" if w == "pooled" else w, "Specification": label},
                            result(model_id, SCI)))
    adjustment = write(pd.DataFrame(rows), "figure_4_national_adjustment").rename(
        columns={"ci_low": "low", "ci_high": "high"})

    fig, ax = plt.subplots(figsize=(8.6, 6.9))
    waves = ["2016", "2019", "2022", "2023", "Pooled, four waves"]
    y_base = np.array([4.25, 3.25, 2.25, 1.25, 0.0])
    baseline = adjustment[adjustment["Specification"] == "Matched baseline"].set_index("Wave").reindex(waves)
    adjusted = adjustment[adjustment["Specification"] == "Fully adjusted"].set_index("Wave").reindex(waves)
    y_baseline = y_base + 0.10
    y_adjusted = y_base - 0.10
    label_x = 1.38
    for yb, ya, wave in zip(y_baseline, y_adjusted, waves):
        ax.plot([baseline.loc[wave, "estimate"], adjusted.loc[wave, "estimate"]], [yb, ya],
                color=RULE, linewidth=1.2, zorder=1)
    for data, y, color, marker in [(baseline, y_baseline, MUTED, "o"), (adjusted, y_adjusted, TEAL, "s")]:
        ax.errorbar(
            data["estimate"], y,
            xerr=[data["estimate"] - data["low"], data["high"] - data["estimate"]],
            fmt=marker, color=color, ecolor=color, markersize=5, linewidth=1.4, capsize=0,
        )
        for yi, (_, r) in zip(y, data.iterrows()):
            ax.text(label_x, yi, f"{r['estimate']:.3f} [{r['low']:.3f}, {r['high']:.3f}]",
                    va="center", ha="left", fontsize=9.5, color=color)
    ax.set_yticks(y_base, ["2016", "2019", "2022", "2023", "Pooled\n(four waves)"])
    ax.set_xlabel("MAI points per 1 SD higher 2026 GlobalSCI")
    ax.set_xticks(np.arange(0, 1.41, 0.2))
    ax.set_xlim(0, 1.94)
    ax.set_ylim(-0.55, y_base.max() + 0.65)
    ax.spines["left"].set_bounds(y_base.min(), y_base.max())
    ax.spines["bottom"].set_bounds(0, 1.4)
    clean_axis(ax, zero_line=True)
    ax.annotate("Matched baseline", xy=(baseline.iloc[0]["estimate"], y_baseline[0]), xytext=(-12, 18),
                textcoords="offset points", ha="right", fontsize=9, color=MUTED,
                arrowprops={"arrowstyle": "-", "color": MUTED, "lw": 0.7})
    ax.annotate("Fully adjusted", xy=(adjusted.iloc[0]["estimate"], y_adjusted[0]), xytext=(14, 18),
                textcoords="offset points", ha="left", fontsize=9, color=TEAL,
                arrowprops={"arrowstyle": "-", "color": TEAL, "lw": 0.7})
    fig.subplots_adjust(left=0.20, right=0.95, top=0.96, bottom=0.14)
    save_figure(fig, "figure_4_national_adjustment")


def figure_2(regional_table: pd.DataFrame, un_table: pd.DataFrame) -> tuple[float, int]:
    """Geographic scale: concordance, variance decomposition, contextual and UN-region estimates."""
    concordance = pd.read_csv(DERIVED / "country_regional_sci_concordance.csv", keep_default_na=False)
    concordance["regional_z"] = pd.to_numeric(
        concordance["regional_to_country_sci_population_weighted_2026_z"], errors="coerce")
    concordance["global_sci_2026_z"] = pd.to_numeric(concordance["global_sci_2026_z"], errors="coerce")
    concordance = concordance.dropna(subset=["regional_z", "global_sci_2026_z"]).copy()
    concordance["discrepancy"] = concordance["regional_z"] - concordance["global_sci_2026_z"]
    r_value = stats.pearsonr(concordance["global_sci_2026_z"], concordance["regional_z"]).statistic

    variance = pd.read_csv(DERIVED / "regional_sci_variance_decomposition.csv")
    regional = regional_table.rename(columns={"ci_low": "low", "ci_high": "high"}).copy()
    regional["Wave"] = regional["wave"].replace({"Pooled": POOLED})
    regional["Component"] = regional["component"].map(
        {"Country average": "Country-average regional SCI", "Regional deviation": "Regional SCI deviation"})
    un_region = un_table.rename(columns={"un_region": "UN world region", "ci_low": "low",
                                         "ci_high": "high", "countries": "Countries"})

    fig = plt.figure(figsize=(7.0, 7.25))
    ax_a = fig.add_axes([0.12, 0.68, 0.36, 0.28])
    ax_b = fig.add_axes([0.60, 0.68, 0.37, 0.28])
    ax_c = fig.add_axes([0.12, 0.39, 0.85, 0.20])
    ax_d = fig.add_axes([0.27, 0.07, 0.70, 0.24])
    axes = np.array([[ax_a, ax_b], [ax_c, ax_d]], dtype=object)

    ax = axes[0, 0]
    ax.scatter(concordance["global_sci_2026_z"], concordance["regional_z"], s=15, color="#8f9996",
               alpha=0.75, edgecolors="none")
    lower = min(concordance["global_sci_2026_z"].min(), concordance["regional_z"].min())
    upper = max(concordance["global_sci_2026_z"].max(), concordance["regional_z"].max())
    ax.plot([lower, upper], [lower, upper], linestyle="--", linewidth=1, color=RULE)
    labels = concordance.nlargest(3, "discrepancy").index.tolist() + concordance.nsmallest(3, "discrepancy").index.tolist()
    for idx in labels:
        r = concordance.loc[idx]
        ax.annotate(r["iso2"], (r["global_sci_2026_z"], r["regional_z"]), xytext=(4, 4),
                    textcoords="offset points", fontsize=8, color=TEXT)
    ax.text(0.04, 0.94, f"Pearson r = {r_value:.3f}", transform=ax.transAxes, fontsize=9, color=TEAL_DARK)
    ax.set_xlabel("Country GlobalSCI (SD)", fontsize=9)
    ax.set_ylabel("Population-aggregated ADM1 SCI (SD)", fontsize=9)
    ax.tick_params(labelsize=8)
    add_panel_label(ax, "A")
    clean_axis(ax)

    ax = axes[0, 1]
    shares = variance[variance["component"].isin(["Between countries", "Within countries"])].copy()
    shares["percent"] = 100 * shares["share_of_total_variance"]
    shares = shares.set_index("component").loc[["Between countries", "Within countries"]].reset_index()
    colors = [TEAL, CORAL]
    y = [1, 0]
    ax.barh(y, shares["percent"], color=colors, height=0.48)
    for yi, (_, r), color in zip(y, shares.iterrows(), colors):
        ax.text(r["percent"] + 1.2, yi, f"{r['percent']:.0f}%", va="center", color=color,
                fontsize=10, fontweight="bold")
    ax.set_yticks(y, ["Between countries", "Within countries"], fontsize=9)
    ax.set_xlim(0, 100)
    ax.set_xlabel("Share of RegionalSCI variance (%)", fontsize=9)
    ax.tick_params(axis="x", labelsize=8)
    add_panel_label(ax, "B")
    clean_axis(ax)

    ax = axes[1, 0]
    waves = [*WAVES, POOLED]
    labels_wave = [*WAVES, "Pooled"]
    base_y = np.arange(len(waves))[::-1]
    offsets = {"Country-average regional SCI": 0.16, "Regional SCI deviation": -0.16}
    styles = {
        "Country-average regional SCI": (TEAL, "o", "Country mean"),
        "Regional SCI deviation": (CORAL, "s", "Within-country deviation"),
    }
    for component, (color, marker, direct_label) in styles.items():
        data = regional[regional["Component"] == component].set_index("Wave").loc[waves]
        yy = base_y + offsets[component]
        ax.errorbar(data["estimate"], yy, xerr=[data["estimate"] - data["low"], data["high"] - data["estimate"]],
                    fmt=marker, color=color, ecolor=color, linewidth=1.5, markersize=5.5, capsize=0)
        label_y = {"Country-average regional SCI": 0.31, "Regional SCI deviation": -0.31}[component]
        label_x = {"Country-average regional SCI": 0.69, "Regional SCI deviation": 0.31}[component]
        ax.plot([data["estimate"].iloc[-1], label_x - 0.015], [yy[-1], label_y], color=color, linewidth=0.8)
        ax.text(label_x, label_y, direct_label, ha="left", va="center", fontsize=8, color=color)
    ax.set_yticks(base_y, labels_wave, fontsize=8.5)
    ax.set_xlabel("MAI points per overall RegionalSCI SD", fontsize=9)
    ax.tick_params(axis="x", labelsize=8)
    add_panel_label(ax, "C")
    ax.set_xlim(-0.15, 1.12)
    ax.set_ylim(-0.55, 4.55)
    clean_axis(ax, zero_line=True)

    ax = axes[1, 1]
    plot_order = [
        "Europe", "Sub-Saharan Africa", "Eastern and South-Eastern Asia", "Latin America and Caribbean",
        "Northern Africa and Western Asia", "Australia and New Zealand", "Central and Southern Asia",
    ]
    data = un_region.set_index("UN world region").loc[plot_order].reset_index()
    display_names = {
        "Europe": "Europe",
        "Sub-Saharan Africa": "Sub-Saharan Africa",
        "Eastern and South-Eastern Asia": "E & SE Asia",
        "Latin America and Caribbean": "Latin America/Carib.",
        "Northern Africa and Western Asia": "N Africa & W Asia",
        "Australia and New Zealand": "Australia & NZ",
        "Central and Southern Asia": "Central & S Asia",
    }
    y = np.arange(len(data))[::-1]
    x_min, x_max = -1.05, 2.0
    for yi, (_, r) in zip(y, data.iterrows()):
        sparse = r["Countries"] < 10
        color = TEAL if r["UN world region"] == "Europe" else MUTED
        markerface = "none" if sparse else color
        estimate = np.clip(r["estimate"], x_min + 0.04, x_max - 0.04)
        if r["low"] >= x_min and r["high"] <= x_max and not sparse:
            ax.hlines(yi, r["low"], r["high"], color=color, linewidth=1.5)
        ax.plot(estimate, yi, marker="o", markersize=5.5, markerfacecolor=markerface, markeredgecolor=color,
                markeredgewidth=1.2, linestyle="none")
        if r["estimate"] < x_min:
            ax.text(x_min + 0.03, yi, f"< {r['estimate']:.2f}", va="center", ha="left", fontsize=8, color=MUTED)
    ax.set_yticks(y, [f"{display_names[name]}  (n={int(n)})" for name, n in zip(data["UN world region"], data["Countries"])],
                  fontsize=8)
    ax.set_xlim(x_min, x_max)
    ax.set_xlabel("MAI points per GlobalSCI SD", fontsize=9)
    ax.tick_params(axis="x", labelsize=8)
    add_panel_label(ax, "D")
    ax.text(0.99, 0.02, "Open points: <10 countries", transform=ax.transAxes, ha="right", fontsize=8, color=MUTED)
    clean_axis(ax, zero_line=True)
    for ax in axes.flat:
        for text in ax.texts:
            if text.get_text() in {"A", "B", "C", "D"}:
                text.set_fontsize(11)
    save_figure(fig, "figure_2_geographic_scale", tight=False)
    return float(r_value), len(concordance)


def figure_5() -> None:
    """Stranger contact: main associations, interaction estimators, and margins."""
    prefix = "country_contact_2026_"

    def point(model_id: str, estimator: str, *term: str) -> dict:
        name = interaction_term(model_id, estimator, *term) if len(term) > 1 else term[0]
        r = result(model_id, name, estimator)
        return {"estimate": r["estimate"], "low": r["conf_low"], "high": r["conf_high"]}

    main_rows = pd.DataFrame([
        {**point(prefix + "contact_only", LMM, "stranger_contact_z"), "label": "Contact only: stranger contact"},
        {**point(prefix + "additive", LMM, SCI), "label": "Additive model: GlobalSCI"},
        {**point(prefix + "additive", LMM, "stranger_contact_z"), "label": "Additive model: stranger contact"},
    ])
    interaction_rows = pd.DataFrame([
        {**point(prefix + "interaction", LMM, SCI, "stranger_contact_z"), "label": "Country random intercept"},
        {**point(prefix + "interaction", OLS, SCI, "stranger_contact_z"), "label": "Clustered OLS"},
        {**point(prefix + "within_between_contact", LMM, SCI, "stranger_contact_within_z"), "label": "Contextual mixed model"},
        {**point(prefix + "interaction_country_fe", COUNTRY_FE, SCI, "stranger_contact_within_z"), "label": "Country fixed effects"},
    ])
    margins = pd.read_csv(MODELS / "contact_interaction_margins.csv")
    margins = margins[margins["model_id"] == prefix + "interaction"]

    fig = plt.figure(figsize=(12.4, 8.7))
    grid = fig.add_gridspec(2, 2, height_ratios=[0.82, 1.25], hspace=0.45, wspace=0.38)
    ax_a = fig.add_subplot(grid[0, 0])
    ax_b = fig.add_subplot(grid[0, 1])
    ax_c = fig.add_subplot(grid[1, :])

    y = np.arange(len(main_rows))[::-1]
    ax_a.errorbar(main_rows["estimate"], y,
                  xerr=[main_rows["estimate"] - main_rows["low"], main_rows["high"] - main_rows["estimate"]],
                  fmt="o", color=TEAL, ecolor=TEAL, linewidth=1.4, markersize=5.5, capsize=0)
    ax_a.set_yticks(y, main_rows["label"])
    ax_a.set_xlabel("MAI points per 1 SD")
    add_panel_label(ax_a, "A")
    ax_a.set_xlim(0, 0.88)
    clean_axis(ax_a, zero_line=True)

    y = np.arange(len(interaction_rows))[::-1]
    colors = [TEAL, MUTED, TEAL, TEAL]
    for yi, (_, r), color in zip(y, interaction_rows.iterrows(), colors):
        ax_b.hlines(yi, r["low"], r["high"], color=color, linewidth=1.4)
        ax_b.plot(r["estimate"], yi, "o", color=color, markersize=5.5)
    ax_b.set_yticks(y, interaction_rows["label"])
    ax_b.set_xlabel("SCI × contact coefficient")
    add_panel_label(ax_b, "B")
    ax_b.set_xlim(-0.105, 0.09)
    clean_axis(ax_b, zero_line=True)

    levels = [(-1, CORAL, "Less frequent contact (-1 SD)"), (0, MUTED, "Average contact"),
              (1, TEAL, "More frequent contact (+1 SD)")]
    margin_label_y = {-1: 0.86, 0: 1.05, 1: 1.24}
    for level, color, label in levels:
        data = margins[margins["stranger_contact_z"] == level].sort_values("sci_z")
        ax_c.fill_between(data["sci_z"], data["conf_low"], data["conf_high"], color=color, alpha=0.10, linewidth=0)
        ax_c.plot(data["sci_z"], data["estimate"], color=color, linewidth=1.8)
        last = data.iloc[-1]
        label_x = last["sci_z"] + 0.14
        ax_c.plot([last["sci_z"], label_x - 0.02], [last["estimate"], margin_label_y[level]], color=color, linewidth=0.8)
        ax_c.text(label_x, margin_label_y[level], label, va="center", fontsize=8, color=color)
    ax_c.axhline(0, color=RULE, linewidth=0.8)
    ax_c.set_xlabel("2026 country GlobalSCI (SD units)")
    ax_c.set_ylabel("Model-implied MAI difference (points)")
    add_panel_label(ax_c, "C")
    ax_c.set_xlim(-2.05, 2.85)
    clean_axis(ax_c)
    save_figure(fig, "figure_5_stranger_contact")


def main() -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    TABLES.mkdir(parents=True, exist_ok=True)
    set_plot_style()
    table_1()
    primary = table_8()
    table_9()
    table_10()
    table_11()
    table_12()
    regional = table_13()
    table_14()
    table_15()
    un_regions = table_16()
    table_17()
    figure_1(primary)
    concordance_r, concordance_n = figure_2(regional, un_regions)
    figure_4()
    figure_5()
    in_text_values(primary, regional, concordance_r, concordance_n)
    print(f"Wrote figures to {FIGURES} and tables to {TABLES}")


if __name__ == "__main__":
    main()
