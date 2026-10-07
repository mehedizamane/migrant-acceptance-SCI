#!/usr/bin/env Rscript
# Regional models (Figure 2C; SI Table 12, bottom panel; SI Tables 13 and 14).
#
# Fits MAI on the country-average component of RegionalSCI and each ADM1 region's
# deviation from it, both in the same total-RegionalSCI SD units, with country and
# ADM1 random intercepts. Uses the deterministic direct ADM1 linkage built by
# 09_link_gallup_regions_to_adm1.py.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: 10_fit_regional_models.R INPUT_CSV OUTPUT_DIR")
}
input_path <- normalizePath(args[[1]], mustWork = TRUE)
output_dir <- args[[2]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(lme4)
  library(lmerTest)
})

set.seed(20260816)
options(contrasts = c("contr.treatment", "contr.poly"))
d <- fread(input_path, na.strings = c("NA_REAL", ""), showProgress = TRUE)
d[, country_id := factor(country_id)]
d[, exposure_region_id := factor(exposure_region_id)]
d[, year_f := relevel(factor(year_wave), ref = "2016")]
d[, female_f := factor(female, levels = c(0, 1), labels = c("Male", "Female"))]
d[, foreign_born_f := factor(
  foreign_born, levels = c(0, 1),
  labels = c("Born in survey country", "Born outside survey country")
)]
d[, education_f := factor(education_level, levels = c(1, 2, 3))]
d[, income_f := factor(income_quintile, levels = c(1, 2, 3, 4, 5))]
d[, urbanicity_f := factor(
  urbanicity, levels = c("rural", "town_or_semi_dense", "city")
)]

years <- c(2016L, 2019L, 2022L, 2023L)
individual <- c(
  "age_10", "female_f", "foreign_born_f", "education_f", "income_f",
  "urbanicity_f"
)
region_controls <- c("survey_region_log_population_z", "survey_region_log_density_z")
random_terms <- "(1 | country_id) + (1 | exposure_region_id)"

coefficient_rows <- list()
registry_rows <- list()
convergence_rows <- list()
skipped_rows <- list()

normalize_weights <- function(use, pooled = FALSE) {
  if (pooled) {
    use[, analysis_weight := weight_raw / sum(weight_raw), by = .(year_wave, country_id)]
    counts <- use[, .(country_count = uniqueN(country_id)), by = year_wave]
    use <- merge(use, counts, by = "year_wave", all.x = TRUE, sort = FALSE)
    use[, analysis_weight := analysis_weight / country_count / uniqueN(year_wave)]
    use[, country_count := NULL]
  } else {
    use[, analysis_weight := weight_raw / sum(weight_raw), by = country_id]
  }
  use[, mixed_weight := analysis_weight / mean(analysis_weight)]
  use
}

formula_variables <- function(outcome, rhs) {
  all.vars(as.formula(paste(outcome, "~", paste(rhs, collapse = " + "))))
}

prepare_use <- function(frame, outcome, rhs, pooled = FALSE) {
  required <- unique(c(
    formula_variables(outcome, rhs), "country_id", "exposure_region_id", "weight_raw"
  ))
  use <- copy(frame[
    complete.cases(frame[, ..required]) & is.finite(weight_raw) & weight_raw > 0
  ])
  use <- droplevels(use)
  normalize_weights(use, pooled = pooled)
}

fit_model <- function(frame, outcome, rhs, model_id, sample, year_wave, pooled = FALSE) {
  use <- prepare_use(frame, outcome, rhs, pooled)
  if (nrow(use) == 0 || uniqueN(use$country_id) < 2 || uniqueN(use$exposure_region_id) < 2) {
    skipped_rows[[length(skipped_rows) + 1L]] <<- data.table(
      model_id = model_id, reason = "No estimable complete-case multilevel sample"
    )
    return(NULL)
  }
  fixed_formula <- as.formula(paste(outcome, "~", paste(rhs, collapse = " + ")))
  mixed_formula <- as.formula(paste(
    outcome, "~", paste(rhs, collapse = " + "), "+", random_terms
  ))
  control <- lmerControl(
    optimizer = "bobyqa", optCtrl = list(maxfun = 200000),
    check.conv.singular = .makeCC(action = "message", tol = 1e-5)
  )
  mixed <- tryCatch(
    lmerTest::lmer(
      mixed_formula, data = use, weights = mixed_weight, REML = FALSE,
      control = control
    ),
    error = function(error) error
  )
  if (inherits(mixed, "error")) {
    skipped_rows[[length(skipped_rows) + 1L]] <<- data.table(
      model_id = model_id, reason = conditionMessage(mixed)
    )
    return(NULL)
  }
  summary_table <- coef(summary(mixed))
  beta <- fixef(mixed)
  covariance <- as.matrix(vcov(mixed))[names(beta), names(beta), drop = FALSE]
  se <- sqrt(diag(covariance))
  dfs <- setNames(as.numeric(summary_table[names(beta), "df"]), names(beta))
  critical <- ifelse(is.finite(dfs), qt(.975, pmax(dfs, 1)), qnorm(.975))
  coefficient_rows[[length(coefficient_rows) + 1L]] <<- data.table(
    model_id = model_id, estimator = "weighted_lmm_ml", term = names(beta),
    estimate = as.numeric(beta), std_error = as.numeric(se),
    degrees_of_freedom = as.numeric(dfs), statistic = as.numeric(beta / se),
    p_value = as.numeric(summary_table[names(beta), "Pr(>|t|)"]),
    conf_low = as.numeric(beta - critical * se),
    conf_high = as.numeric(beta + critical * se),
    inference = "Satterthwaite model-based inference for weighted lmer"
  )
  vc <- as.data.table(as.data.frame(VarCorr(mixed)))
  diagonal <- vc[is.na(var2)]
  total <- sum(diagonal$vcov)
  country_variance <- diagonal[grp == "country_id", vcov][1]
  exposure_variance <- diagonal[grp == "exposure_region_id", vcov][1]
  gradient <- mixed@optinfo$derivs$gradient
  messages <- mixed@optinfo$conv$lme4$messages
  convergence_rows[[length(convergence_rows) + 1L]] <<- data.table(
    model_id = model_id, estimator = "weighted_lmm_ml",
    singular_fit = isSingular(mixed, tol = 1e-5),
    max_absolute_gradient = if (is.null(gradient)) NA_real_ else max(abs(gradient)),
    convergence_message = if (is.null(messages)) "" else paste(messages, collapse = "; ")
  )
  registry_rows[[length(registry_rows) + 1L]] <<- data.table(
    model_id = model_id, estimator = "weighted_lmm_ml", estimation = "maximum likelihood",
    sample = sample, year_wave = as.character(year_wave), outcome = outcome,
    fixed_formula = paste(deparse(fixed_formula), collapse = " "),
    random_effects = random_terms,
    n_respondents = nrow(use),
    n_countries = uniqueN(use$country_id),
    n_exposure_clusters = uniqueN(use$exposure_region_id),
    country_icc = ifelse(is.finite(country_variance), country_variance / total, NA_real_),
    exposure_region_icc = ifelse(is.finite(exposure_variance), exposure_variance / total, NA_real_),
    weighting = if (pooled) {
      "Gallup weights normalized within country-wave, equal country totals within wave, equal wave totals; rescaled to mean one for lmer"
    } else {
      "Gallup weights normalized to equal country totals; rescaled to mean one for lmer"
    },
    pooled = pooled,
    log_likelihood = as.numeric(logLik(mixed)), aic = AIC(mixed), bic = BIC(mixed)
  )
  invisible(mixed)
}

decomposition_rhs <- c(
  "survey_region_sci_between_total_sd", "survey_region_sci_within_total_sd",
  individual, region_controls
)

# Figure 2C and SI Table 13: each wave and the pooled four-wave model.
for (year in years) {
  fit_model(
    d[year_wave == year], "migrant_acceptance_index", decomposition_rhs,
    sprintf("regional_within_between_%s", year), "direct ADM1 linkage", year
  )
}
fit_model(
  d, "migrant_acceptance_index", c(decomposition_rhs, "year_f"),
  "regional_within_between_pooled", "direct ADM1 linkage",
  "2016+2019+2022+2023", pooled = TRUE
)

# SI Table 14: pooled model restricted to exact code- or name-matched regions.
fit_model(
  d[exact_only_sample == TRUE], "migrant_acceptance_index", c(decomposition_rhs, "year_f"),
  "regional_within_between_exact_only", "exact code or name matches only",
  "2016+2019+2022+2023", pooled = TRUE
)

# SI Table 12, bottom panel: geography-adjusted regional deviation, same linkage.
geography_rhs <- c("regional_sci_geography_adjusted_within_z", individual, region_controls)
for (year in years) {
  fit_model(
    d[year_wave == year], "migrant_acceptance_index", geography_rhs,
    sprintf("regional_geography_adjusted_within_%s", year), "direct ADM1 linkage", year
  )
}
fit_model(
  d, "migrant_acceptance_index", c(geography_rhs, "year_f"),
  "regional_geography_adjusted_within_pooled", "direct ADM1 linkage",
  "2016+2019+2022+2023", pooled = TRUE
)

skipped <- if (length(skipped_rows)) {
  rbindlist(skipped_rows, fill = TRUE)
} else {
  data.table(model_id = character(), reason = character())
}
fwrite(rbindlist(coefficient_rows, fill = TRUE), file.path(output_dir, "model_coefficients.csv"))
fwrite(rbindlist(registry_rows, fill = TRUE), file.path(output_dir, "model_registry.csv"))
fwrite(rbindlist(convergence_rows, fill = TRUE), file.path(output_dir, "model_convergence.csv"))
fwrite(skipped, file.path(output_dir, "skipped_models.csv"))
message("Saved ", length(registry_rows), " regional model specifications to ", normalizePath(output_dir))
