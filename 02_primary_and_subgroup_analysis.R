# 02_primary_and_subgroup_analysis.R
# Fit the main congenital-syphilis model, check subgroup results, and draw figures.
# Run from the v7 analysis root after reviewing the data checks.
# Sections: 1) primary model; 2) covariance-aware contrasts; 3) subgroup audit
# and Figure 3; 4) Figure 2 and country selection; 5) subgroup annual curves;
# 6) adjustment-graph check.

# 0. Shared data setup
# Common paths and input panels for scripts run from the analysis root.
# Edit the two data roots below when files are stored elsewhere.

# Raw exports and prepared analysis inputs are stored separately.
raw_data_dir <- Sys.getenv("BPG_RAW_DATA_DIR", unset = file.path("data", "raw_data"))
processed_data_dir <- Sys.getenv("BPG_PROCESSED_DATA_DIR", unset = file.path("data", "processed_data"))

if (!file.exists(file.path(processed_data_dir, "data_combined.Rdata")) ||
    !file.exists(file.path(processed_data_dir, "data_combined_hiv.Rdata"))) {
  stop("Prepared panels are required in data/processed_data/. Raw exports belong in data/raw_data/.")
}


results_dir <- "results"
figures_dir <- "figures"
tables_dir <- "tables"
logs_dir <- "logs"

dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(logs_dir, showWarnings = FALSE, recursive = TRUE)

cs_environment <- new.env(parent = emptyenv())
load(file.path(processed_data_dir, "data_combined.Rdata"), envir = cs_environment)
cs_panel <- get("data_combined", envir = cs_environment)

hiv_environment <- new.env(parent = emptyenv())
load(file.path(processed_data_dir, "data_combined_hiv.Rdata"), envir = hiv_environment)
hiv_panel <- get("data_combined_hiv", envir = hiv_environment)

study_years <- 2010:2021
break_year <- 2016L

# Shared figure style used by plots in this file.
publication_theme <- function(base_size = 11) {
  ggplot2::theme_classic(base_size = base_size, base_family = "Arial") +
    ggplot2::theme(
      text = ggplot2::element_text(colour = "black"),
      axis.line = ggplot2::element_line(colour = "black", linewidth = 0.3),
      axis.ticks = ggplot2::element_line(colour = "black", linewidth = 0.3),
      axis.title = ggplot2::element_text(face = "bold"),
      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 8)),
      axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 8)),
      plot.title = ggplot2::element_text(face = "bold", size = 12, hjust = 0),
      plot.subtitle = ggplot2::element_text(size = 10, hjust = 0),
      panel.grid = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"),
      legend.position = "bottom",
      legend.title = ggplot2::element_blank(),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA)
    )
}
shortage_colours <- c("Shortage reported" = "#C63D40",
  "No shortage reported" = "#176B9A", "Unclassified" = "#737A83")


# -----------------------------------------------------------------------------
# Primary weighted controlled interrupted time-series model
# -----------------------------------------------------------------------------
local({
  # Reproduce the original v3 weighted country-level congenital-syphilis LMM.
  # Run from the analysis root: Rscript the primary-analysis section of 02_primary_and_subgroup_analysis.R

  library(lme4)
  library(lmerTest)

  # 1. Keep the frozen 89-country sample and the 2010-2021 study years.
  analysis_data <- cs_panel[
    cs_panel$year %in% study_years & cs_panel$shortage %in% c(0, 1),
  ]

  analysis_data$location_id <- factor(analysis_data$location_id)
  analysis_data$gbd_superregion <- factor(analysis_data$gbd_superregion)
  analysis_data$year <- as.integer(analysis_data$year)
  analysis_data$time_c <- analysis_data$year - break_year
  analysis_data$post <- as.integer(analysis_data$year >= break_year)
  analysis_data$rate <- 1000 * analysis_data$val / analysis_data$population

  # Weights are proportional to the GBD Birth population, as in v3.
  analysis_data$birth_weight <- analysis_data$population / sum(analysis_data$population)

  core_covariates <- c(
    "hdi", "sdi", "gni", "syphilis_prev_female", "syphilis_prev_male",
    "abor_misc_rate", "gbd_superregion"
  )
  model_columns <- c(
    "rate", "time_c", "shortage", "post", "population", "birth_weight",
    "location_id", core_covariates
  )
  analysis_data <- analysis_data[complete.cases(analysis_data[, model_columns]), ]

  if (nrow(analysis_data) != 1068L || length(unique(analysis_data$location_id)) != 89L) {
    stop("The original model sample differs from the frozen 89-country, 1,068-row sample.")
  }
  if (anyDuplicated(analysis_data[c("location_id", "year")])) {
    stop("The original model data have duplicate country-year rows.")
  }

  # 2. Fit the v3 primary model. The triple interaction is the differential
  # post-2016 change in annual slope for shortage versus no-shortage countries.
  primary_model <- lmerTest::lmer(
    rate ~ time_c * shortage * post +
      hdi + sdi + gni + gbd_superregion +
      syphilis_prev_female + syphilis_prev_male + abor_misc_rate +
      (1 + time_c | location_id),
    data = analysis_data,
    weights = birth_weight,
    REML = FALSE
  )

  print(summary(primary_model))

  # 3. Save full-precision coefficients with matching t-based confidence limits.
  coefficient_table <- as.data.frame(summary(primary_model)$coefficients)
  coefficient_table$term <- rownames(coefficient_table)
  rownames(coefficient_table) <- NULL
  names(coefficient_table)[names(coefficient_table) == "Std. Error"] <- "std_error"
  names(coefficient_table)[names(coefficient_table) == "t value"] <- "t_value"
  names(coefficient_table)[names(coefficient_table) == "Pr(>|t|)"] <- "p_value"
  coefficient_table$ci_lower <- with(
    coefficient_table,
    Estimate - qt(0.975, df) * std_error
  )
  coefficient_table$ci_upper <- with(
    coefficient_table,
    Estimate + qt(0.975, df) * std_error
  )
  coefficient_table$formula <- paste(
    deparse(formula(primary_model)), collapse = " "
  )
  coefficient_table$n_country_years <- nrow(analysis_data)
  coefficient_table$n_countries <- length(unique(analysis_data$location_id))
  coefficient_table$weight_sum <- sum(analysis_data$birth_weight)
  coefficient_table$singular_fit <- lme4::isSingular(primary_model, tol = 1e-4)

  write.csv(
    coefficient_table,
    file.path(results_dir, "original_primary_lmm_coefficients.csv"),
    row.names = FALSE,
    na = ""
  )
  write.csv(
    analysis_data,
    file.path(results_dir, "original_primary_lmm_sample.csv"),
    row.names = FALSE,
    na = ""
  )
  saveRDS(primary_model, file.path(results_dir, "original_primary_lmm.rds"))
  writeLines(capture.output(sessionInfo()), file.path(logs_dir, "session_info.txt"))

  cat("Saved the original primary LMM and its analysis sample.\n")
})


# -----------------------------------------------------------------------------
# Covariance-aware primary contrasts and counterfactual trajectories
# -----------------------------------------------------------------------------
local({
  # This diagnostic does not overwrite the historical publication estimates.
  library(lme4)
  library(lmerTest)
  model <- readRDS("results/original_primary_lmm.rds")
  contrast <- setNames(rep(0, length(fixef(model))), names(fixef(model)))
  contrast["shortage:post"] <- 1
  contrast["time_c:shortage:post"] <- 5
  checks <- list()
  for (year in 2017:2021) {
    for (group in c("No shortage", "Shortage", "Controlled difference")) {
      contrast[] <- 0
      if (group != "Controlled difference") {
        contrast["post"] <- 1
        contrast["time_c:post"] <- year - 2016
      }
      if (group != "No shortage") {
        contrast["shortage:post"] <- 1
        contrast["time_c:shortage:post"] <- year - 2016
      }
      check <- as.data.frame(contest1D(model, contrast, confint = TRUE))
      check$year <- year
      check$group <- group
      checks[[length(checks) + 1L]] <- check
    }
  }
  checks <- do.call(rbind, checks)
  write.csv(checks, "results/counterfactual_inference_audit.csv", row.names = FALSE)
  print(checks[checks$year == 2021, ])


  # -----------------------------------------------------------------------------
  # Covariance-aware primary counterfactual estimates
  # -----------------------------------------------------------------------------

  # Reconstruct conditional trajectory and ratio uncertainty using joint covariance.
  # Variance parameters are held at their fitted ML values, as in a Wald analysis.
  library(lme4)
  library(Matrix)
  model <- readRDS("results/original_primary_lmm.rds")
  data <- read.csv("results/original_primary_lmm_sample.csv")
  stopifnot(nrow(data) == nrow(model.frame(model)))

  # 1. Mixed-model equations in the spherical random-effect parameterisation.
  X <- getME(model, "X")
  ZL <- as.matrix(getME(model, "Z") %*% getME(model, "Lambda"))
  design <- cbind(X, ZL)
  weights <- model.weights(model.frame(model))
  information <- crossprod(design * sqrt(weights))
  p <- ncol(X)
  q <- ncol(ZL)
  information[(p + 1):(p + q), (p + 1):(p + q)] <-
    information[(p + 1):(p + q), (p + 1):(p + q)] + diag(q)
  scale <- sqrt(diag(information))
  scaled_information <- information / outer(scale, scale)
  covariance <- sigma(model)^2 * solve(scaled_information) / outer(scale, scale)
  fixed_covariance <- as.matrix(vcov(model))
  stopifnot(max(abs(covariance[1:p, 1:p] - fixed_covariance)) < 1e-7)
  parameters <- c(fixef(model), getME(model, "u"))

  # 2. Counterfactual fixed-effect design; random effects and annual covariates stay.
  counterfactual <- model.frame(model)
  counterfactual$post <- 0
  Xcf <- model.matrix(delete.response(terms(model, fixed.only = TRUE)),
                      counterfactual, contrasts.arg = attr(X, "contrasts"))
  Xcf <- Xcf[, colnames(X), drop = FALSE]
  design_cf <- cbind(Xcf, ZL)
  stopifnot(max(abs(as.vector(design %*% parameters) - fitted(model))) < 1e-8)
  stopifnot(max(abs(as.vector(design_cf %*% parameters) -
    predict(model, newdata = counterfactual, re.form = NULL))) < 1e-8)

  # 3. Weighted group trajectories and their complete joint covariance.
  predictions <- list()
  ratios <- list()
  for (year in 2010:2021) {
    group_designs <- list()
    for (shortage in 0:1) {
      rows <- which(data$year == year & data$shortage == shortage)
      annual_weights <- data$population[rows] / sum(data$population[rows])
      for (scenario in c("Observed", "Counterfactual post=0")) {
        current <- if (scenario == "Observed") design else design_cf
        average <- colSums(current[rows, , drop = FALSE] * annual_weights)
        estimate <- sum(average * parameters)
        se <- sqrt(as.numeric(t(average) %*% covariance %*% average))
        name <- paste0(scenario, ": ", if (shortage == 0) "non-shortage" else "shortage")
        predictions[[length(predictions) + 1L]] <- data.frame(year = year,
          scenario = name, estimate = estimate, se = se,
          ci_lower = estimate - qnorm(.975) * se,
          ci_upper = estimate + qnorm(.975) * se)
        group_designs[[paste(shortage, scenario)]] <- average
      }
    }
    if (year >= 2017) {
      gradients <- list()
      log_ratios <- c()
      for (shortage in 0:1) {
        fitted_design <- group_designs[[paste(shortage, "Observed")]]
        cf_design <- group_designs[[paste(shortage, "Counterfactual post=0")]]
        fitted_rate <- sum(fitted_design * parameters)
        cf_rate <- sum(cf_design * parameters)
        stopifnot(fitted_rate > 0, cf_rate > 0)
        log_ratios[shortage + 1] <- log(fitted_rate / cf_rate)
        gradients[[shortage + 1]] <- fitted_design / fitted_rate - cf_design / cf_rate
      }
      for (index in 1:3) {
        gradient <- if (index <= 2) gradients[[index]] else gradients[[2]] - gradients[[1]]
        log_estimate <- if (index <= 2) log_ratios[index] else log_ratios[2] - log_ratios[1]
        log_se <- sqrt(as.numeric(t(gradient) %*% covariance %*% gradient))
        ratios[[length(ratios) + 1L]] <- data.frame(year = year,
          group = c("No shortage", "Shortage", "Controlled ratio")[index],
          estimate = exp(log_estimate), std_error = exp(log_estimate) * log_se,
          log_std_error = log_se,
          lower = exp(log_estimate - qnorm(.975) * log_se),
          upper = exp(log_estimate + qnorm(.975) * log_se))
      }
    }
  }
  predictions <- do.call(rbind, predictions)
  ratios <- do.call(rbind, ratios)
  write.csv(predictions, "results/counterfactual_covariance_predictions.csv", row.names = FALSE)
  write.csv(ratios, "results/counterfactual_covariance_ratios.csv", row.names = FALSE)
  print(ratios[ratios$year == 2021, ])
})


# -----------------------------------------------------------------------------
# Subgroup heterogeneity audit and Figure 3
# -----------------------------------------------------------------------------
local({
  # This audit retains the archived subgroup coefficients; it fits no new model.
  data <- read.csv("results/figure3_subgroup_data.csv", stringsAsFactors = FALSE)
  current_group <- ""
  for (row in seq_len(nrow(data))) {
    if (nzchar(data$Group[row])) current_group <- data$Group[row]
    data$analysis_group[row] <- current_group
  }
  results <- list()
  for (group in unique(data$analysis_group)) {
    rows <- data[data$analysis_group == group, ]
    for (contrast in c("level", "slope")) {
      estimate <- rows[[paste0("Estimate_", contrast)]]
      se <- rows[[paste0("SE_", contrast)]]
      weights <- 1 / se^2
      pooled <- sum(weights * estimate) / sum(weights)
      Q <- sum(weights * (estimate - pooled)^2)
      df <- length(estimate) - 1
      p <- pchisq(Q, df = df, lower.tail = FALSE)
      label <- rows[[paste0("interaction_", contrast)]][1]
      if (label == "<0.001") stopifnot(p < .001)
      else stopifnot(abs(p - as.numeric(label)) < .0005)
      results[[length(results) + 1L]] <- data.frame(group = group,
        contrast = contrast, strata = length(estimate), Q = Q, df = df,
        p_value = p, publication_label = label)
    }
  }
  write.csv(do.call(rbind, results), "logs/subgroup_heterogeneity_audit.csv", row.names = FALSE)
  cat("All Figure 3 Q-test labels agree with their archived coefficients.\n")


  # -----------------------------------------------------------------------------
  # Figure 3 subgroup forest plot
  # -----------------------------------------------------------------------------

  # Draw Figure 3 from the verified original subgroup estimates.
  # Cochran Q compares separately fitted strata; it is not a pooled interaction test.

  library(forestploter)
  library(grid)
  library(ggplot2)

  # 1. Read the published model coefficients and their matched t intervals.
  combined <- read.csv("results/figure3_subgroup_data.csv", stringsAsFactors = FALSE)
  combined$interaction_level[is.na(combined$interaction_level)] <- ""
  combined$interaction_slope[is.na(combined$interaction_slope)] <- ""

  # 2. Prepare labels while preserving the original row order and P values.
  plot_table <- data.frame(
    Group = combined$Group, Subgroup = combined$Subgroup,
    `Level change (95% CI)` = sprintf("%.3f (%.3f, %.3f)",
      combined$Estimate_level, combined$CI_lower_level, combined$CI_upper_level),
    `Level contrast` = paste(rep(" ", 32), collapse = ""),
    `P for heterogeneity\nLevel change` = combined$interaction_level,
    `Slope change (95% CI)` = sprintf("%.3f (%.3f, %.3f)",
      combined$Estimate_slope, combined$CI_lower_slope, combined$CI_upper_slope),
    `Slope contrast` = paste(rep(" ", 32), collapse = ""),
    `P for heterogeneity\nSlope change` = combined$interaction_slope,
    check.names = FALSE
  )

  # 3. Retain the original forest-plot layout, axes, colours and shading.
  colours <- ggsci::pal_lancet()(2)
  style <- forest_theme(base_size = 10,
    refline_gp = gpar(lty = "dashed", col = "#666666"),
    ci_pch = 15, ci_col = colours[1], ci_fill = colours[1],
    arrow_type = "closed", ci_Theight = .25)
  figure <- forest(plot_table,
    est = list(combined$Estimate_level, combined$Estimate_slope),
    lower = list(combined$CI_lower_level, combined$CI_lower_slope),
    upper = list(combined$CI_upper_level, combined$CI_upper_slope),
    ci_column = c(4, 7), xlim = list(c(-.42, .28), c(-.08, .34)),
    ticks_at = list(c(-.4, -.2, 0, .2), c(0, .1, .2, .3)),
    arrow_lab = c("Decreased rate", "Increased rate"), theme = style)
  figure <- edit_plot(figure, row = seq_len(nrow(combined)), col = 7,
    which = "ci", gp = gpar(pch = 21, fill = colours[2], col = colours[2]))
  figure <- edit_plot(figure, row = which(combined$Group == ""),
    which = "background", gp = gpar(fill = "white"))
  figure <- edit_plot(figure, row = which(combined$Group != ""),
    which = "background", gp = gpar(fill = "gray80"))

  # 4. Save publication images without re-estimating any model.
  ggsave("figures/Figure_03_subgroups.png", figure, width = 14, height = 9,
    dpi = 400, bg = "white", limitsize = FALSE)
  ggsave("figures/Figure_03_subgroups.pdf", figure, width = 14, height = 9,
    bg = "white", limitsize = FALSE)
})


# -----------------------------------------------------------------------------
# Primary and country-selection figures
# -----------------------------------------------------------------------------
local({
  # Recreate manuscript Figure 2 and Supplementary Figure S1 from saved model output.

  library(ggplot2)

  # Figure 2 retains the reproduced v3 point predictions, with corrected
  # covariance-aware uncertainty from the counterfactual and breakpoint sections of 02_primary_and_subgroup_analysis.R and 03_sensitivity_analyses.R.
  prediction_file <- file.path(results_dir, "counterfactual_covariance_predictions.csv")
  primary_data_file <- file.path(results_dir, "original_primary_lmm_sample.csv")
  if (!file.exists(prediction_file) || !file.exists(primary_data_file)) {
    stop("Run the primary-analysis section of 02_primary_and_subgroup_analysis.R and the counterfactual and breakpoint sections of 02_primary_and_subgroup_analysis.R and 03_sensitivity_analyses.R first.")
  }
  predictions <- read.csv(prediction_file, stringsAsFactors = FALSE)
  primary_data <- read.csv(primary_data_file, stringsAsFactors = FALSE)
  primary_data$group <- ifelse(primary_data$shortage == 1,
                               "Shortage reported", "No shortage reported")
  observed_groups <- split(primary_data, list(primary_data$year, primary_data$group),
                           drop = TRUE)
  observed <- do.call(rbind, lapply(observed_groups, function(x) data.frame(
    year = x$year[1], group = x$group[1],
    observed_rate = weighted.mean(x$rate, x$population),
    population = sum(x$population)
  )))

  predictions$group <- ifelse(grepl("non-shortage", predictions$scenario),
                              "No shortage reported", "Shortage reported")
  predictions$trajectory <- ifelse(grepl("^Observed", predictions$scenario),
                                   "Fitted trajectory", "Counterfactual trajectory")
  pre_fitted <- predictions[predictions$year < break_year &
                              predictions$trajectory == "Fitted trajectory", ]
  pre_counterfactual <- predictions[predictions$year < break_year &
                                      predictions$trajectory == "Counterfactual trajectory", ]
  if (nrow(pre_fitted) != nrow(pre_counterfactual) ||
      any(abs(pre_fitted$estimate - pre_counterfactual$estimate) > 1e-10)) {
    stop("Pre-break fitted and counterfactual paths do not match.")
  }

  counterfactual_pre <- predictions[predictions$trajectory ==
                                      "Counterfactual trajectory" & predictions$year <= break_year, ]
  counterfactual_pre$trajectory <- "Fitted trajectory"
  counterfactual_pre$segment <- paste(counterfactual_pre$group, "pre-break")
  post_break <- predictions[predictions$year >= break_year, ]
  post_break$segment <- paste(post_break$group, post_break$trajectory, "post-break")
  line_data <- rbind(counterfactual_pre, post_break)
  level_jumps <- merge(
    predictions[predictions$year == break_year &
                  predictions$trajectory == "Fitted trajectory", ],
    predictions[predictions$year == break_year &
                  predictions$trajectory == "Counterfactual trajectory", ],
    by = "group", suffixes = c("_fitted", "_counterfactual")
  )

  figure2 <- ggplot() +
    geom_ribbon(data = line_data,
                aes(x = year, ymin = ci_lower, ymax = ci_upper,
                    fill = group, group = segment), alpha = 0.14) +
    geom_vline(xintercept = break_year, linetype = 2, colour = "grey45") +
    geom_line(data = line_data,
              aes(x = year, y = estimate, colour = group,
                  linetype = trajectory, group = segment), linewidth = 0.8) +
    geom_segment(data = level_jumps,
                 aes(x = break_year, xend = break_year,
                     y = estimate_counterfactual, yend = estimate_fitted,
                     colour = group), linewidth = 0.8) +
    geom_point(data = observed,
               aes(x = year, y = observed_rate, colour = group, shape = group),
               size = 2) +
    scale_colour_manual(values = c("Shortage reported" = "#C63D40",
                                   "No shortage reported" = "#176B9A")) +
    scale_fill_manual(values = c("Shortage reported" = "#C63D40",
                                 "No shortage reported" = "#176B9A")) +
    scale_linetype_manual(values = c("Fitted trajectory" = "solid",
                                     "Counterfactual trajectory" = "longdash")) +
    scale_x_continuous(breaks = study_years) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
    labs(x = "Year", y = "Congenital syphilis cases per 1,000 live births",
         colour = NULL, fill = NULL, shape = NULL, linetype = NULL) +
    publication_theme(base_size = 11) +
    theme(legend.position = "bottom", legend.key.width = grid::unit(1.8, "cm"),
          axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(figures_dir, "Figure_02_v3_presentation.png"), figure2,
         width = 8.5, height = 5, dpi = 300)
  ggsave(file.path(figures_dir, "Figure_02_v3_presentation.pdf"), device = grDevices::cairo_pdf, figure2,
         width = 8.5, height = 5)

  # Supplementary Figure S5: display the exact sample counts used in the CITS.
  selection <- data.frame(
    group = factor(c("Shortage reported", "No shortage reported", "Unclassified"),
                   levels = c("Shortage reported", "No shortage reported", "Unclassified")),
    countries = c(38, 51, 115)
  )
  selection_plot <- ggplot(selection, aes(x = group, y = countries, fill = group)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = countries), vjust = -0.35, size = 4) +
    scale_fill_manual(values = c("#C63D40", "#176B9A", "#999999")) +
    labs(x = NULL, y = "Number of countries",
         title = "Country selection by BPG shortage-reporting status") +
    publication_theme(base_size = 11) + theme(legend.position = "none")
  ggsave(file.path(figures_dir, "Figure_S05_Country_selection.png"), selection_plot,
         width = 7, height = 4.8, dpi = 300)
  ggsave(file.path(figures_dir, "Figure_S05_Country_selection.pdf"), device = grDevices::cairo_pdf, selection_plot,
         width = 7, height = 4.8)

  # S15 reports both endpoints and both contrasts, with identical adjustment.


  # Supplementary Figure S15 is generated by the HIV figure section of 04_additional_analyses.R after HIV model results are available.
})


# -----------------------------------------------------------------------------
# Subgroup annual-rate curves, including final Supplementary Figure S4B
# -----------------------------------------------------------------------------
local({
  # Recreate subgroup curves, alternative-breakpoint figures, and Supplementary Figure S4B.
  # Plotting only: subgroup rates are sums of cases divided by sums of live births.
  library(ggplot2)
  library(patchwork)

  # S4B and S7–S9: group curves use strata defined in 2021.
  d <- cs_panel[cs_panel$year %in% study_years, ]
  d$group <- ifelse(is.na(d$shortage), "Unclassified",
    ifelse(d$shortage == 1, "Shortage reported", "No shortage reported"))
  strata <- cs_panel[cs_panel$year == 2021,
    c("location_id", "gbd_superregion", "income_level", "hdi_level", "sdi_level2021")]
  figure_numbers <- c("S4B", "S07", "S08", "S09")
  variables <- c("gbd_superregion", "income_level", "hdi_level", "sdi_level2021")
  for (index in seq_along(variables)) {
    variable <- variables[index]
    labels <- strata[c("location_id", variable)]
    names(labels)[2] <- "stratum"
    plot_panel <- merge(d[c("location_id", "year", "val", "population", "group")],
                        labels, by = "location_id")
    totals <- aggregate(cbind(val, population) ~ year + group + stratum,
                        data = plot_panel, FUN = sum)
    totals$rate <- 1000 * totals$val / totals$population
    write.csv(totals, file.path(results_dir,
      paste0("figure_", tolower(figure_numbers[index]), "_curve_data.csv")), row.names = FALSE)
    curve <- ggplot(totals, aes(x = year, y = rate, colour = group, group = group)) +
      geom_vline(xintercept = c(2014, 2016), linetype = "dashed", colour = "grey65") +
      geom_line(linewidth = 0.8) + geom_point(size = 1.6) +
      scale_colour_manual(values = shortage_colours) +
      scale_x_continuous(breaks = c(2010, 2012, 2014, 2016, 2018, 2020, 2021)) +
      scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
      facet_wrap(~ stratum, ncol = 2) +
      labs(x = "Year", y = "Congenital syphilis cases per 1,000 live births", colour = NULL) +
      publication_theme()
    ggsave(file.path(figures_dir, paste0("Figure_", figure_numbers[index], "_curves.png")),
           curve, width = 10, height = if (index == 1) 8.5 else 6.5, dpi = 300)
    ggsave(file.path(figures_dir, paste0("Figure_", figure_numbers[index], "_curves.pdf")),
           curve, width = 10, height = if (index == 1) 8.5 else 6.5,
           device = grDevices::cairo_pdf)
  }


  # -----------------------------------------------------------------------------
  # Final S4B panel with region title
  # -----------------------------------------------------------------------------

  # Supplementary Figure S4B uses saved v7 descriptive data; no model is refitted.
  library(ggplot2)

  # (1) Annual regional rates and shortage-report groups
  d <- read.csv("results/figure_s4b_curve_data.csv")
  d$group <- factor(d$group, levels = c("No shortage reported",
    "Shortage reported", "Unclassified"))
  d$stratum <- factor(d$stratum, levels = unique(d$stratum))

  # (2) Reproduce the trajectories and add the missing panel title
  figure <- ggplot(d, aes(year, rate, colour = group, group = group)) +
    geom_vline(xintercept = c(2014, 2016), linetype = "dashed",
      colour = "grey65", linewidth = 0.4) +
    geom_line(linewidth = 0.65) + geom_point(size = 1.9) +
    facet_wrap(~stratum, ncol = 2) +
    scale_colour_manual(values = shortage_colours, drop = FALSE) +
    scale_x_continuous(breaks = c(seq(2010, 2020, 2), 2021)) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
    labs(title = "(B) Annual congenital syphilis incidence by GBD super-region",
      x = "Year", y = "Congenital syphilis cases per 1,000 live births") +
    publication_theme(10) +
    theme(strip.text = element_text(size = 9.5),
      legend.position = "bottom", plot.title = element_text(size = 12))

  # (3) Export the publication figure with the existing rates and colours
  ggsave("figures/Figure_S4B_curves.png", figure, width = 10,
    height = 8.8, dpi = 300, bg = "white")
  ggsave("figures/Figure_S4B_curves.pdf", figure, width = 10,
    height = 8.8, device = grDevices::cairo_pdf, bg = "white")
})


# -----------------------------------------------------------------------------
# Congenital-syphilis adjustment-graph check and supplied rendering
# -----------------------------------------------------------------------------
local({
  # Display the original causal graph with the verified covariate measure label.
  library(dagitty)
  dag_text <- readLines("data/dag/figure_s1_dagitty.txt", warn = FALSE)
  graph <- dagitty(paste(dag_text, collapse = "\n"))
  adjustment <- c("GNI per capita", "HDI", "SDI",
    "Incidence of maternal abortion and miscarriage",
    "Prevalence of syphilis of men aged 15-49",
    "Prevalence of syphilis of women aged 15-49")
  stopifnot(isAdjustmentSet(graph, adjustment))
  writeLines(c("The six-covariate set is valid under the specified graph.",
    capture.output(adjustmentSets(graph, type = "minimal"))),
    "results/cs_dag_adjustment_check.txt")
  # The publication panels retain the official DAGitty web renderer's layout.
  # Adjusted covariates use red fill, as described in data/dag/README.md.
  # Keep the exported SVG, PNG and PDF together in figures/dagitty_web.
  # Paste data/dag/figure_s1a_web_dagitty.txt into dagitty.net/dags.html
  # to inspect the same graph and adjustment markings.
  file.copy("figures/dagitty_web/figure_s1a.png",
            "figures/Figure_S1A_CS_adjustment_DAG.png", overwrite = TRUE)
  file.copy("figures/dagitty_web/figure_s1a.pdf",
            "figures/Figure_S1A_CS_adjustment_DAG.pdf", overwrite = TRUE)
})
