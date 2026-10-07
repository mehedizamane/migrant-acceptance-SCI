#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2 || length(args) > 3) {
  stop("Usage: 06_fit_country_models.R INPUT_CSV OUTPUT_DIR [FAMILY]")
}

input_path <- normalizePath(args[[1]], mustWork = TRUE)
output_dir <- args[[2]]
family_filter <- if (length(args) == 3) strsplit(args[[3]], ",", fixed = TRUE)[[1]] else "all"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

`%||%` <- function(x, y) if (is.null(x)) y else x

suppressPackageStartupMessages({
  library(data.table)
  library(lme4)
  library(lmerTest)
  library(fixest)
})

set.seed(20260814)
options(contrasts = c("contr.treatment", "contr.poly"))

d <- fread(input_path, na.strings = c("NA_REAL", ""), showProgress = TRUE)
d[, country_id := factor(country_id)]
d[, adm1_id := factor(adm1_id)]
d[, year_f := relevel(factor(year_wave), ref = "2016")]
d[, female_f := factor(female, levels = c(0, 1), labels = c("Male", "Female"))]
d[, foreign_born_f := factor(
  foreign_born,
  levels = c(0, 1),
  labels = c("Born in survey country", "Born outside survey country")
)]
d[, education_f := factor(education_level, levels = c(1, 2, 3))]
d[, income_f := factor(income_quintile, levels = c(1, 2, 3, 4, 5))]
d[, urbanicity_f := factor(
  urbanicity,
  levels = c("rural", "town_or_semi_dense", "city")
)]
d[, wb_region_f := factor(wb_region)]
d[, un_region_f := factor(un_region)]

years <- c(2016L, 2019L, 2022L, 2023L)
individual <- c(
  "age_10", "female_f", "foreign_born_f", "education_f", "income_f",
  "urbanicity_f"
)
macro <- list(
  development = c("log_gdp_per_capita_z", "log_population_z"),
  region = "wb_region_f",
  digital = c(
    "internet_users_pct_z", "mobile_cellular_subscriptions_per_100_z",
    "fixed_broadband_subscriptions_per_100_z",
    "log_secure_internet_servers_per_million_people_z"
  ),
  openness = c(
    "urban_population_pct_z", "trade_pct_gdp_z",
    "tourism_arrivals_per_capita_z"
  ),
  migration = c("migrant_stock_share_z", "refugee_stock_share_z")
)
full_macro <- unlist(macro, use.names = FALSE)
full_national_controls <- c(
  full_macro, "employment_to_population_ratio_15plus_z"
)

coefficient_rows <- list()
registry_rows <- list()
convergence_rows <- list()
contrast_rows <- list()
margin_rows <- list()
skipped_rows <- list()

run_family <- function(name) "all" %in% family_filter || name %in% family_filter

slug <- function(x) gsub("(^_+|_+$)", "", gsub("[^a-z0-9]+", "_", tolower(x)))

normalize_model_weights <- function(use, pooled = FALSE) {
  if (pooled) {
    use[, analysis_weight := weight_raw / sum(weight_raw), by = .(year_wave, country_id)]
    wave_country_counts <- use[, uniqueN(country_id), by = year_wave]
    setnames(wave_country_counts, "V1", "country_count")
    use <- merge(use, wave_country_counts, by = "year_wave", all.x = TRUE, sort = FALSE)
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

prepare_use <- function(frame, outcome, rhs, random_terms, pooled = FALSE) {
  required <- unique(c(
    formula_variables(outcome, rhs), "country_id", "weight_raw",
    if (grepl("adm1_id", random_terms, fixed = TRUE)) "adm1_id" else character(0)
  ))
  use <- copy(frame[complete.cases(frame[, ..required]) & is.finite(weight_raw) & weight_raw > 0])
  use <- droplevels(use)
  normalize_model_weights(use, pooled = pooled)
}

record_coefficients <- function(model_id, estimator, beta, covariance, df_values,
                                p_values, inference) {
  covariance <- covariance[names(beta), names(beta), drop = FALSE]
  se <- sqrt(diag(covariance))
  df_values <- df_values[names(beta)]
  critical <- ifelse(is.finite(df_values), qt(0.975, pmax(df_values, 1)), qnorm(0.975))
  coefficient_rows[[length(coefficient_rows) + 1L]] <<- data.table(
    model_id = model_id,
    estimator = estimator,
    term = names(beta),
    estimate = as.numeric(beta),
    std_error = as.numeric(se),
    degrees_of_freedom = as.numeric(df_values),
    statistic = as.numeric(beta / se),
    p_value = as.numeric(p_values[names(beta)]),
    conf_low = as.numeric(beta - critical * se),
    conf_high = as.numeric(beta + critical * se),
    inference = inference
  )
}

fit_pair <- function(frame, outcome, rhs, random_terms, model_id, metadata,
                     pooled = FALSE, reml = FALSE, ols = FALSE) {
  use <- prepare_use(frame, outcome, rhs, random_terms, pooled = pooled)
  n_clusters <- uniqueN(use$country_id)
  if (nrow(use) == 0 || n_clusters < 2) {
    skipped_rows[[length(skipped_rows) + 1L]] <<- data.table(
      model_id = model_id, reason = "No estimable complete-case sample"
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
      mixed_formula, data = use, weights = mixed_weight,
      REML = reml, control = control
    ),
    error = function(error) error
  )
  # Weighted OLS with country-clustered CR1 errors is fitted only where the
  # paper reports it (primary models, SI Table 10; contact models, SI Table 17).
  fit_ols <- ols
  ols <- if (!fit_ols) NULL else tryCatch(
    feols(
      fixed_formula, data = use, weights = ~mixed_weight,
      cluster = ~country_id,
      ssc = ssc(adj = TRUE, cluster.adj = TRUE, fixef.K = "none")
    ),
    error = function(error) error
  )

  base_registry <- data.table(
    model_id = model_id,
    analysis_family = metadata$analysis_family,
    sample = metadata$sample,
    release = as.character(metadata$release),
    year_wave = as.character(metadata$year_wave),
    outcome = outcome,
    outcome_label = metadata$outcome_label,
    focal_terms = paste(metadata$focal_terms, collapse = " | "),
    fixed_formula = paste(deparse(fixed_formula), collapse = " "),
    random_effects = random_terms,
    n_respondents = nrow(use),
    n_countries = n_clusters,
    n_adm1_regions = uniqueN(use$adm1_id),
    weighting = if (pooled) {
      "Gallup weights normalized to equal country totals within wave and equal wave totals; rescaled to mean one"
    } else {
      "Gallup weights normalized to equal country totals; rescaled to mean one"
    },
    pooled = pooled
  )

  if (inherits(mixed, "error")) {
    skipped_rows[[length(skipped_rows) + 1L]] <<- data.table(
      model_id = model_id, estimator = "weighted_lmm_ml",
      reason = conditionMessage(mixed)
    )
  } else {
    summary_table <- coef(summary(mixed))
    beta <- fixef(mixed)
    covariance <- as.matrix(vcov(mixed))[names(beta), names(beta), drop = FALSE]
    df_values <- setNames(as.numeric(summary_table[names(beta), "df"]), names(beta))
    p_values <- setNames(as.numeric(summary_table[names(beta), "Pr(>|t|)"]), names(beta))
    record_coefficients(
      model_id, "weighted_lmm_ml", beta, covariance, df_values, p_values,
      "Satterthwaite model-based inference for weighted lmer"
    )
    vc <- as.data.table(as.data.frame(VarCorr(mixed)))
    diagonal <- vc[is.na(var2)]
    total_random_residual <- sum(diagonal$vcov)
    country_variance <- diagonal[grp == "country_id", vcov][1] %||% NA_real_
    adm1_variance <- diagonal[grp == "adm1_id", vcov][1] %||% NA_real_
    residual_variance <- diagonal[grp == "Residual", vcov][1] %||% NA_real_
    country_icc <- if (length(country_variance) == 0 || !is.finite(country_variance)) {
      NA_real_
    } else country_variance / total_random_residual
    adm1_icc <- if (length(adm1_variance) == 0 || !is.finite(adm1_variance)) {
      NA_real_
    } else adm1_variance / total_random_residual
    gradient <- mixed@optinfo$derivs$gradient
    max_gradient <- if (is.null(gradient)) NA_real_ else max(abs(gradient))
    messages <- mixed@optinfo$conv$lme4$messages %||% character(0)
    convergence_rows[[length(convergence_rows) + 1L]] <<- data.table(
      model_id = model_id,
      estimator = "weighted_lmm_ml",
      singular_fit = isSingular(mixed, tol = 1e-5),
      max_absolute_gradient = max_gradient,
      convergence_message = paste(messages, collapse = "; ")
    )
    registry_rows[[length(registry_rows) + 1L]] <<- copy(base_registry)[, `:=`(
      estimator = "weighted_lmm_ml",
      estimation = if (reml) "REML" else "maximum likelihood",
      country_icc = country_icc,
      adm1_icc = adm1_icc,
      residual_variance = residual_variance,
      log_likelihood = as.numeric(logLik(mixed)),
      aic = AIC(mixed),
      bic = BIC(mixed)
    )]
  }

  if (is.null(ols)) {
    # OLS not requested for this specification.
  } else if (inherits(ols, "error")) {
    skipped_rows[[length(skipped_rows) + 1L]] <<- data.table(
      model_id = model_id, estimator = "weighted_ols_country_cr1",
      reason = conditionMessage(ols)
    )
  } else {
    beta <- coef(ols)
    covariance <- as.matrix(vcov(ols))[names(beta), names(beta), drop = FALSE]
    df_value <- n_clusters - 1
    df_values <- setNames(rep(df_value, length(beta)), names(beta))
    se <- sqrt(diag(covariance))
    p_values <- setNames(2 * pt(-abs(beta / se), df = df_value), names(beta))
    record_coefficients(
      model_id, "weighted_ols_country_cr1", beta, covariance, df_values,
      p_values, "Country-clustered CR1 inference with cluster df"
    )
    registry_rows[[length(registry_rows) + 1L]] <<- copy(base_registry)[, `:=`(
      estimator = "weighted_ols_country_cr1",
      estimation = "weighted least squares",
      country_icc = NA_real_, adm1_icc = NA_real_, residual_variance = NA_real_,
      log_likelihood = as.numeric(logLik(ols)), aic = AIC(ols), bic = BIC(ols)
    )]
  }
  list(mixed = if (inherits(mixed, "error")) NULL else mixed,
       ols = if (inherits(ols, "error")) NULL else ols,
       use = use)
}

fit_country_fe_ols <- function(frame, outcome, rhs, model_id, metadata,
                               pooled = FALSE) {
  use <- prepare_use(frame, outcome, rhs, "", pooled = pooled)
  n_clusters <- uniqueN(use$country_id)
  if (nrow(use) == 0 || n_clusters < 2) {
    skipped_rows[[length(skipped_rows) + 1L]] <<- data.table(
      model_id = model_id, estimator = "weighted_ols_country_fe_cr1",
      reason = "No estimable complete-case sample"
    )
    return(NULL)
  }
  fixed_formula <- as.formula(paste(
    outcome, "~", paste(rhs, collapse = " + "), "| country_id"
  ))
  fit <- tryCatch(
    feols(
      fixed_formula, data = use, weights = ~mixed_weight,
      cluster = ~country_id,
      ssc = ssc(adj = TRUE, cluster.adj = TRUE, fixef.K = "none")
    ),
    error = function(error) error
  )
  if (inherits(fit, "error")) {
    skipped_rows[[length(skipped_rows) + 1L]] <<- data.table(
      model_id = model_id, estimator = "weighted_ols_country_fe_cr1",
      reason = conditionMessage(fit)
    )
    return(NULL)
  }
  beta <- coef(fit)
  covariance <- as.matrix(vcov(fit))[names(beta), names(beta), drop = FALSE]
  df_value <- n_clusters - 1
  df_values <- setNames(rep(df_value, length(beta)), names(beta))
  se <- sqrt(diag(covariance))
  p_values <- setNames(2 * pt(-abs(beta / se), df = df_value), names(beta))
  record_coefficients(
    model_id, "weighted_ols_country_fe_cr1", beta, covariance, df_values,
    p_values, "Country-fixed-effects OLS with country-clustered CR1 inference"
  )
  registry_rows[[length(registry_rows) + 1L]] <<- data.table(
    model_id = model_id,
    analysis_family = metadata$analysis_family,
    sample = metadata$sample,
    release = as.character(metadata$release),
    year_wave = as.character(metadata$year_wave),
    outcome = outcome,
    outcome_label = metadata$outcome_label,
    focal_terms = paste(metadata$focal_terms, collapse = " | "),
    fixed_formula = paste(deparse(fixed_formula), collapse = " "),
    random_effects = "None; country fixed effects absorbed",
    n_respondents = nrow(use),
    n_countries = n_clusters,
    n_adm1_regions = uniqueN(use$adm1_id),
    weighting = if (pooled) {
      "Gallup weights normalized to equal country totals within wave and equal wave totals; rescaled to mean one"
    } else {
      "Gallup weights normalized to equal country totals; rescaled to mean one"
    },
    pooled = pooled,
    estimator = "weighted_ols_country_fe_cr1",
    estimation = "weighted least squares with absorbed country fixed effects",
    country_icc = NA_real_, adm1_icc = NA_real_, residual_variance = NA_real_,
    log_likelihood = as.numeric(logLik(fit)), aic = AIC(fit), bic = BIC(fit)
  )
  fit
}

record_linear_contrast <- function(fit, model_id, estimator, contrast, label,
                                   family, release, region = "") {
  if (is.null(fit)) return(invisible(NULL))
  beta <- if (inherits(fit, "merMod")) fixef(fit) else coef(fit)
  covariance <- as.matrix(vcov(fit))[names(beta), names(beta), drop = FALSE]
  vector <- setNames(rep(0, length(beta)), names(beta))
  shared <- intersect(names(contrast), names(vector))
  vector[shared] <- contrast[shared]
  estimate <- sum(vector * beta)
  se <- sqrt(as.numeric(t(vector) %*% covariance %*% vector))
  registry <- rbindlist(registry_rows, fill = TRUE)
  current_model_id <- model_id
  current_estimator <- estimator
  clusters <- registry[
    model_id == current_model_id & estimator == current_estimator,
    n_countries
  ][1]
  df <- if (estimator == "weighted_ols_country_cr1") clusters - 1 else Inf
  critical <- if (is.finite(df)) qt(0.975, df) else qnorm(0.975)
  p <- if (is.finite(df)) 2 * pt(-abs(estimate / se), df) else 2 * pnorm(-abs(estimate / se))
  contrast_rows[[length(contrast_rows) + 1L]] <<- data.table(
    model_id = model_id, estimator = estimator, contrast = label,
    analysis_family = family, release = as.character(release), region = region,
    estimate = estimate, std_error = se, degrees_of_freedom = df,
    p_value = p, conf_low = estimate - critical * se,
    conf_high = estimate + critical * se
  )
}

message("Loaded ", format(nrow(d), big.mark = ","), " model-input rows.")

# -------------------------------------------------------------------------
# Country GlobalSCI: 2026-primary models and matched release sensitivities.
# -------------------------------------------------------------------------
if (run_family("country_main")) {
  primary_sci <- "global_sci_2026_all_origins_z"
  primary_rhs <- c(primary_sci, individual)
  for (year in years) {
    fit_pair(
      d[year_wave == year & balanced_primary_2026_country == TRUE],
      "migrant_acceptance_index", primary_rhs, "(1 | country_id)",
      sprintf("country_primary_2026_%s", year),
      list(
        analysis_family = "country_global_sci_primary_balanced",
        sample = "Same 2026-SCI-eligible countries in all four waves; no ADM1 requirement",
        release = 2026, year_wave = year, outcome_label = "Migrant Acceptance Index (0-9)",
        focal_terms = primary_sci
      ), ols = TRUE
    )
    fit_pair(
      d[year_wave == year], "migrant_acceptance_index", primary_rhs,
      "(1 | country_id)", sprintf("country_available_2026_%s", year),
      list(
        analysis_family = "country_global_sci_available",
        sample = "All 2026-SCI-eligible countries in this wave; no ADM1 requirement",
        release = 2026, year_wave = year, outcome_label = "Migrant Acceptance Index (0-9)",
        focal_terms = primary_sci
      )
    )
  }
  fit_pair(
    d[balanced_primary_2026_country == TRUE], "migrant_acceptance_index",
    c(primary_sci, "year_f", individual), "(1 | country_id)",
    "country_pooled_primary_2026",
    list(
      analysis_family = "country_global_sci_primary_pooled",
      sample = "Balanced 2026-SCI country sample across four waves; no ADM1 requirement",
      release = 2026, year_wave = "2016+2019+2022+2023",
      outcome_label = "Migrant Acceptance Index (0-9)", focal_terms = primary_sci
    ), pooled = TRUE, ols = TRUE
  )

  # The release sensitivity uses identical countries and respondents and the
  # common-origin standardization for the independently scaled 2020 and 2026 files.
  for (release in c(2026, 2020)) {
    sci <- paste0("global_sci_", release, "_z")
    for (year in years) {
      fit_pair(
        d[year_wave == year & balanced_release_comparison_country == TRUE],
        "migrant_acceptance_index", c(sci, individual), "(1 | country_id)",
        sprintf("country_release_matched_%s_%s", release, year),
        list(
          analysis_family = "country_global_sci_release_comparison",
          sample = "Same countries and respondents for the 2020 and 2026 SCI releases",
          release = release, year_wave = year,
          outcome_label = "Migrant Acceptance Index (0-9)", focal_terms = sci
        )
      )
    }
    fit_pair(
      d[balanced_release_comparison_country == TRUE], "migrant_acceptance_index",
      c(sci, "year_f", individual), "(1 | country_id)",
      sprintf("country_release_matched_pooled_%s", release),
      list(
        analysis_family = "country_global_sci_release_comparison_pooled",
        sample = "Matched 2020/2026 release-comparison countries across four waves",
        release = release, year_wave = "2016+2019+2022+2023",
        outcome_label = "Migrant Acceptance Index (0-9)", focal_terms = sci
      ), pooled = TRUE
    )
  }

}

# -------------------------------------------------------------------------
# Dyadic-geography-adjusted country SCI sensitivity.
# -------------------------------------------------------------------------
if (run_family("country_geography")) {
  for (release in c(2026, 2020)) {
    sci <- paste0("global_sci_geography_adjusted_", release, "_z")
    sample_flag <- if (release == 2026) {
      "balanced_primary_2026_country"
    } else {
      "balanced_release_comparison_country"
    }
    for (year in years) {
      fit_pair(
        d[year_wave == year & get(sample_flag) == TRUE],
        "migrant_acceptance_index", c(sci, individual),
        "(1 | country_id)",
        sprintf("country_geography_adjusted_%s_%s", release, year),
        list(
          analysis_family = "country_geography_adjusted_sensitivity",
          sample = "Balanced countries supported by the common 165-country CEPII panel",
          release = release, year_wave = year,
          outcome_label = "Migrant Acceptance Index (0-9)", focal_terms = sci
        )
      )
    }
    fit_pair(
      d[get(sample_flag) == TRUE], "migrant_acceptance_index",
      c(sci, "year_f", individual), "(1 | country_id)",
      sprintf("country_geography_adjusted_pooled_%s", release),
      list(
        analysis_family = "country_geography_adjusted_sensitivity_pooled",
        sample = "Balanced countries supported by the common 165-country CEPII panel",
        release = release, year_wave = "2016+2019+2022+2023",
        outcome_label = "Migrant Acceptance Index (0-9)", focal_terms = sci
      ), pooled = TRUE
    )
  }
}

# -------------------------------------------------------------------------
# Matched baseline and fully adjusted national-control specifications.
# -------------------------------------------------------------------------
if (run_family("country_adjustment")) {
  primary_sci <- "global_sci_2026_all_origins_z"
  for (year in years) {
    full_rhs <- c(primary_sci, individual, full_national_controls)
    full_required <- unique(c(
      "migrant_acceptance_index", "weight_raw", "country_id", full_rhs
    ))
    common <- copy(d[
      year_wave == year & complete.cases(d[, ..full_required]) &
        is.finite(weight_raw) & weight_raw > 0
    ])
    for (specification in c("matched_baseline", "fully_adjusted")) {
      rhs <- if (specification == "matched_baseline") {
        c(primary_sci, individual)
      } else {
        full_rhs
      }
      fit_pair(
        common, "migrant_acceptance_index", rhs, "(1 | country_id)",
        sprintf("country_national_adjustment_2026_%s_%s", year, specification),
        list(
          analysis_family = "country_national_adjustment",
          sample = "Wave-specific national-control complete cases; no ADM1 requirement",
          release = 2026, year_wave = year,
          outcome_label = "Migrant Acceptance Index (0-9)", focal_terms = primary_sci
        )
      )
    }
  }

  pooled_full_rhs <- c(primary_sci, "year_f", individual, full_national_controls)
  pooled_required <- unique(c(
    "migrant_acceptance_index", "weight_raw", "country_id", pooled_full_rhs
  ))
  pooled_common <- copy(d[
    complete.cases(d[, ..pooled_required]) & is.finite(weight_raw) & weight_raw > 0
  ])
  for (specification in c("matched_baseline", "fully_adjusted")) {
    rhs <- if (specification == "matched_baseline") {
      c(primary_sci, "year_f", individual)
    } else {
      pooled_full_rhs
    }
    fit_pair(
      pooled_common, "migrant_acceptance_index", rhs, "(1 | country_id)",
      sprintf("country_national_adjustment_pooled_2026_%s", specification),
      list(
        analysis_family = "country_national_adjustment_pooled",
        sample = "Pooled national-control complete cases; no ADM1 requirement",
        release = 2026, year_wave = "2016+2019+2022+2023",
        outcome_label = "Migrant Acceptance Index (0-9)", focal_terms = primary_sci
      ), pooled = TRUE
    )
  }

}

# -------------------------------------------------------------------------
# MAI components and strict-MAI sensitivity.
# -------------------------------------------------------------------------
if (run_family("country_outcomes")) {
  components <- c(
    immigrants_living_in_country_good = "Immigrants living in the country: good thing",
    immigrant_neighbor_good = "Immigrant neighbor: good thing",
    immigrant_marry_close_relative_good = "Immigrant marrying close relative: good thing"
  )
  # SI Table 11, Panels A and B: the 2026 release in the balanced primary sample.
  release <- 2026
  sci <- "global_sci_2026_all_origins_z"
  sample_flag <- "balanced_primary_2026_country"
  sample_label <- "Balanced primary 2026 country set; no ADM1 requirement"
  rhs <- c(sci, individual)
  for (outcome in names(components)) {
    for (year in years) {
      fit_pair(
        d[year_wave == year & get(sample_flag) == TRUE], outcome, rhs,
        "(1 | country_id)", sprintf("country_component_%s_%s_%s", release, year, slug(outcome)),
        list(
          analysis_family = "country_mai_components",
          sample = paste(sample_label, "outcome-specific valid responses"),
          release = release, year_wave = year, outcome_label = components[[outcome]], focal_terms = sci
        )
      )
    }
  }
  for (year in years) {
    fit_pair(
      d[year_wave == year & get(sample_flag) == TRUE],
      "migrant_acceptance_index_strict", rhs, "(1 | country_id)",
      sprintf("country_strict_mai_%s_%s", release, year),
      list(
        analysis_family = "strict_mai_sensitivity",
        sample = paste(sample_label, "all three substantive component responses required"),
        release = release, year_wave = year, outcome_label = "Strict Migrant Acceptance Index (0, 3, 6, or 9)",
        focal_terms = sci
      )
    )
  }
}

# -------------------------------------------------------------------------
# 2022 stranger contact, additive and interaction models.
# -------------------------------------------------------------------------
if (run_family("contact")) {
  components <- c(
    "immigrants_living_in_country_good",
    "immigrant_neighbor_good",
    "immigrant_marry_close_relative_good"
  )
  # SI Table 17 and Figure 5: the 2026 release in the 2022 contact sample.
  for (release in c(2026)) {
    sci <- "global_sci_2026_all_origins_z"
    sample_flag <- paste0("contact_eligible_", release)
    use <- d[get(sample_flag) == TRUE]
    contact_models <- list(
      contact_only = c("stranger_contact_z", individual),
      additive = c(sci, "stranger_contact_z", individual),
      interaction = c(paste0(sci, " * stranger_contact_z"), individual),
      within_between_contact = c(
        paste0(sci, " * stranger_contact_within_z"),
        "stranger_contact_country_mean_z", individual
      )
    )
    for (spec in names(contact_models)) {
      result <- fit_pair(
        use, "migrant_acceptance_index", contact_models[[spec]], "(1 | country_id)",
        sprintf("country_contact_%s_%s", release, spec),
        list(
          analysis_family = "country_sci_stranger_contact",
          sample = "2022 country-level contact sample; no ADM1 requirement",
          release = release, year_wave = 2022, outcome_label = "Migrant Acceptance Index (0-9)",
          focal_terms = c(sci, "stranger_contact_z", paste0(sci, ":stranger_contact_z"))
        ), ols = TRUE
      )
      if (spec == "interaction" && !is.null(result$mixed)) {
        beta <- fixef(result$mixed)
        covariance <- as.matrix(vcov(result$mixed))[names(beta), names(beta), drop = FALSE]
        interaction_term <- paste0(sci, ":stranger_contact_z")
        for (contact in c(-1, 0, 1)) {
          for (sci_value in seq(-2, 2, by = 0.1)) {
            contrast <- setNames(rep(0, length(beta)), names(beta))
            contrast[sci] <- sci_value
            contrast["stranger_contact_z"] <- contact
            if (interaction_term %in% names(contrast)) {
              contrast[interaction_term] <- sci_value * contact
            }
            estimate <- sum(contrast * beta)
            se <- sqrt(as.numeric(t(contrast) %*% covariance %*% contrast))
            margin_rows[[length(margin_rows) + 1L]] <- data.table(
              model_id = sprintf("country_contact_%s_interaction", release),
              release = release, sci_z = sci_value, stranger_contact_z = contact,
              estimate = estimate, std_error = se,
              conf_low = estimate - qnorm(0.975) * se,
              conf_high = estimate + qnorm(0.975) * se
            )
          }
        }
      }
    }
    fit_country_fe_ols(
      use, "migrant_acceptance_index",
      c(
        "stranger_contact_within_z",
        paste0(sci, ":stranger_contact_within_z"),
        individual
      ),
      sprintf("country_contact_%s_interaction_country_fe", release),
      list(
        analysis_family = "country_sci_stranger_contact_country_fe_sensitivity",
        sample = "2022 country-level contact sample; no ADM1 requirement",
        release = release, year_wave = 2022,
        outcome_label = "Migrant Acceptance Index (0-9)",
        focal_terms = c(
          "stranger_contact_within_z",
          paste0(sci, ":stranger_contact_within_z")
        )
      )
    )
    for (outcome in components) {
      fit_pair(
        use, outcome, c(paste0(sci, " * stranger_contact_z"), individual),
        "(1 | country_id)", sprintf("country_contact_component_%s_%s", release, slug(outcome)),
        list(
          analysis_family = "country_sci_contact_component_interactions", sample = "2022 release-specific contact sample; component valid responses",
          release = release, year_wave = 2022, outcome_label = outcome,
          focal_terms = c(sci, "stranger_contact_z", paste0(sci, ":stranger_contact_z"))
        )
      )
    }
  }
}

# -------------------------------------------------------------------------
# UN-region heterogeneity from stacked wave models. 2026 is main; 2020 is SI.
# -------------------------------------------------------------------------
if (run_family("un_region")) {
  # SI Table 16 and Figure 2D: equal-wave average GlobalSCI slope within each
  # UN-M49-derived region, from a stacked four-wave model with wave-specific slopes.
  release <- 2026
  sci <- "global_sci_2026_all_origins_z"
  for (region in sort(unique(na.omit(d$un_region)))) {
    use <- d[balanced_primary_2026_country == TRUE & un_region == region]
    model_id <- sprintf("un_region_%s_%s", release, slug(region))
    result <- fit_pair(
      use, "migrant_acceptance_index",
      c(paste0("year_f * ", sci), individual),
      "(1 | country_id)", model_id,
      list(
        analysis_family = "un_region_heterogeneity", sample = "Balanced-country sample within UN region",
        release = release, year_wave = "stacked four waves", outcome_label = "Migrant Acceptance Index (0-9)",
        focal_terms = sci
      ), pooled = TRUE
    )
    if (is.null(result) || is.null(result$mixed)) next
    # Average of the four wave-specific slopes: the 2016 slope is the main term,
    # and each later wave adds its interaction, so each interaction gets weight 1/4.
    average <- setNames(1, sci)
    for (year in setdiff(years, 2016L)) {
      average[paste0("year_f", year, ":", sci)] <- 1 / length(years)
    }
    record_linear_contrast(
      result$mixed, model_id, "weighted_lmm_ml", average, "Equal-wave average SCI slope",
      "un_region_heterogeneity", release, region
    )
  }
}

coefficients <- rbindlist(coefficient_rows, fill = TRUE)
registry <- rbindlist(registry_rows, fill = TRUE)
convergence <- rbindlist(convergence_rows, fill = TRUE)
contrasts <- rbindlist(contrast_rows, fill = TRUE)
margins <- rbindlist(margin_rows, fill = TRUE)
skipped <- rbindlist(skipped_rows, fill = TRUE)

fwrite(coefficients, file.path(output_dir, "model_coefficients.csv"))
fwrite(registry, file.path(output_dir, "model_registry.csv"))
fwrite(convergence, file.path(output_dir, "model_convergence.csv"))
fwrite(contrasts, file.path(output_dir, "model_linear_contrasts.csv"))
fwrite(margins, file.path(output_dir, "contact_interaction_margins.csv"))
fwrite(skipped, file.path(output_dir, "skipped_models.csv"))

message(
  "Saved ", uniqueN(registry$model_id), " fitted specifications and ",
  nrow(coefficients), " coefficient rows to ", normalizePath(output_dir)
)
