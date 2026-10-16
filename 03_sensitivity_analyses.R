# 03_sensitivity_analyses.R
# Sensitivity analyses for the primary congenital-syphilis association.
# Run from the v7 analysis root after 02_primary_and_subgroup_analysis.R.
# Country-label permutation defaults to 500 iterations and may take time.
# Sections: 1) alternative models and exclusions; 2) WHO service adjustment;
# 3) quantitative bias analysis; 4) 2014/2015 interruption years; 5) curves;
# 6) country-label permutation; 7) QBA heatmap.

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
# Alternative congenital-syphilis model specifications and exclusions
# -----------------------------------------------------------------------------
local({
  # Refit the main congenital-syphilis CITS model under prespecified changes.
  # Run from the analysis root: Rscript the congenital-syphilis sensitivity section of 03_sensitivity_analyses.R

  library(lme4)
  library(lmerTest)
  library(ggplot2)

  # 1. Prepare the classified 89-country panel and the original adjustment set.
  d <- cs_panel[cs_panel$year %in% study_years & cs_panel$shortage %in% c(0, 1), ]
  d$location_id <- factor(d$location_id)
  d$gbd_superregion <- factor(d$gbd_superregion)
  d$rate <- 1000 * d$val / d$population
  d$log_rate <- log(d$rate)
  d$hdi_z <- as.numeric(scale(d$hdi))
  d$sdi_z <- as.numeric(scale(d$sdi))
  d$gni_z <- as.numeric(scale(d$gni))
  d$female_syphilis_z <- as.numeric(scale(d$syphilis_prev_female))
  d$male_syphilis_z <- as.numeric(scale(d$syphilis_prev_male))
  d$abortion_rate_z <- as.numeric(scale(d$abor_misc_rate))

  # 2. Use one model-fitting step so every sensitivity has the same weights,
  #    random effects, estimation method, and target coefficient.
  fit_one <- function(data, label, break_at = 2016L, outcome = "rate",
                      socioeconomic = "all", drop_female_syphilis = FALSE,
                      quadratic = FALSE, calendar_year = FALSE) {
    data$time_c <- data$year - break_at
    data$post <- as.integer(data$year >= break_at)
    data$location_id <- droplevels(factor(data$location_id))
    socioeconomic_terms <- switch(
      socioeconomic,
      all = c("hdi_z", "sdi_z", "gni_z"),
      sdi = "sdi_z",
      hdi = "hdi_z",
      gni = "gni_z",
      stop("Unknown socioeconomic specification: ", socioeconomic)
    )
    covariates <- c(socioeconomic_terms, "gbd_superregion",
                    "female_syphilis_z", "male_syphilis_z", "abortion_rate_z")
    if (drop_female_syphilis) covariates <- setdiff(covariates, "female_syphilis_z")
    model_columns <- c(outcome, "year", "shortage", "population", "location_id", covariates)
    data <- data[complete.cases(data[, model_columns]), ]
    data$birth_weight <- data$population / sum(data$population)
    if (calendar_year) {
      fixed_terms <- c("factor(year)", "shortage", "time_c:shortage",
                       "shortage:post", "time_c:shortage:post", covariates)
    } else {
      fixed_terms <- c("time_c * shortage * post", covariates)
    }
    if (quadratic) {
      fixed_terms <- c(fixed_terms, "I(time_c^2)", "I(time_c^2):shortage",
                       "I(time_c^2):post", "I(time_c^2):shortage:post")
    }
    model_formula <- as.formula(paste(outcome, "~", paste(fixed_terms, collapse = " + "),
                                      "+ (1 + time_c | location_id)"))
    model <- lmerTest::lmer(model_formula, data = data, weights = birth_weight,
                            REML = FALSE,
                            control = lme4::lmerControl(optimizer = "bobyqa"))
    table <- as.data.frame(summary(model)$coefficients)
    normalized_terms <- vapply(strsplit(rownames(table), ":", fixed = TRUE),
                               function(x) paste(sort(x), collapse = ":"), character(1))
    target <- rownames(table)[normalized_terms == "post:shortage:time_c"][1]
    if (is.na(target)) stop("Target three-way interaction is not estimable for: ", label)
    estimate <- table[target, "Estimate"]
    standard_error <- table[target, "Std. Error"]
    degrees_freedom <- table[target, "df"]
    data.frame(
      analysis = label, estimate = estimate, std_error = standard_error,
      df = degrees_freedom,
      ci_lower = estimate - qt(0.975, degrees_freedom) * standard_error,
      ci_upper = estimate + qt(0.975, degrees_freedom) * standard_error,
      p_value = table[target, "Pr(>|t|)"], country_years = nrow(data),
      countries = length(unique(data$location_id)), weight_sum = sum(data$birth_weight),
      singular = lme4::isSingular(model), formula = paste(deparse(model_formula), collapse = " "),
      stringsAsFactors = FALSE
    )
  }

  # 3. Main additional specification and alternative model specifications.
  results <- list()
  results[["full adjusted model"]] <- fit_one(d, "full adjusted model")
  results[["omit female syphilis"]] <- fit_one(
    d, "omit female syphilis", drop_female_syphilis = TRUE
  )
  results[["SDI only"]] <- fit_one(
    d, "SDI only", socioeconomic = "sdi"
  )
  results[["HDI replaces SDI"]] <- fit_one(d, "HDI replaces SDI", socioeconomic = "hdi")
  results[["GNI replaces SDI"]] <- fit_one(d, "GNI replaces SDI", socioeconomic = "gni")
  results[["through 2019"]] <- fit_one(d[d$year <= 2019, ], "through 2019")
  results[["placebo 2012"]] <- fit_one(d[d$year <= 2015, ], "placebo 2012", break_at = 2012L)
  results[["placebo 2013"]] <- fit_one(d[d$year <= 2015, ], "placebo 2013", break_at = 2013L)
  results[["break 2014"]] <- fit_one(d, "break 2014", break_at = 2014L)
  results[["break 2015"]] <- fit_one(d, "break 2015", break_at = 2015L)
  results[["quadratic time"]] <- fit_one(d, "quadratic time", quadratic = TRUE)
  results[["calendar-year effects"]] <- fit_one(d, "calendar-year effects", calendar_year = TRUE)
  results[["log rate"]] <- fit_one(d, "log rate", outcome = "log_rate")
  results[["exclude India"]] <- fit_one(
    d[!grepl("India", d$location_name, ignore.case = TRUE), ], "exclude India"
  )

  # 4. Repeat the model after excluding each GBD super-region in turn.
  regions <- levels(droplevels(d$gbd_superregion))
  for (region in regions) {
    remaining <- d[as.character(d$gbd_superregion) != region, ]
    results[[paste("exclude region", region)]] <- fit_one(
      remaining, paste("exclude region", region)
    )
  }

  # 5. Leave out each classified country once to check whether one country dominates.
  country_ids <- unique(as.character(d$location_id))
  influence_rows <- lapply(country_ids, function(country) {
    reduced <- d[as.character(d$location_id) != country, ]
    result <- fit_one(reduced, paste("leave out", country))
    result$excluded_location_id <- country
    result
  })
  leave_one_country_out <- do.call(rbind, influence_rows)
  full_estimate <- results[["full adjusted model"]]$estimate
  leave_one_country_out$change_from_full <- leave_one_country_out$estimate - full_estimate
  write.csv(leave_one_country_out,
            file.path(results_dir, "cs_leave_one_country_out.csv"),
            row.names = FALSE, na = "")

  # 6. Exclude the highest-baseline-incidence countries in two prespecified tiers.
  baseline_country_rate <- aggregate(
    rate ~ location_id, data = d[d$year %in% 2010:2015, ], FUN = mean
  )
  baseline_country_rate <- baseline_country_rate[order(-baseline_country_rate$rate), ]
  n_top_5 <- ceiling(0.05 * nrow(baseline_country_rate))
  n_top_10 <- ceiling(0.10 * nrow(baseline_country_rate))
  top_5_ids <- as.character(baseline_country_rate$location_id[seq_len(n_top_5)])
  top_10_ids <- as.character(baseline_country_rate$location_id[seq_len(n_top_10)])
  results[["exclude top baseline 5 percent"]] <- fit_one(
    d[!as.character(d$location_id) %in% top_5_ids, ],
    "exclude top baseline 5 percent"
  )
  results[["exclude top baseline 10 percent"]] <- fit_one(
    d[!as.character(d$location_id) %in% top_10_ids, ],
    "exclude top baseline 10 percent"
  )

  # 7. Save estimates, then draw the forest plot from the saved table.
  result_table <- do.call(rbind, results)
  row.names(result_table) <- NULL
  write.csv(result_table, file.path(results_dir, "cs_sensitivity_models.csv"),
            row.names = FALSE, na = "")

  plot_data <- result_table[order(result_table$estimate), ]
  plot_data$analysis <- factor(plot_data$analysis, levels = plot_data$analysis)
  sensitivity_plot <- ggplot(plot_data, aes(x = estimate, y = analysis)) +
    geom_vline(xintercept = 0, colour = "grey55", linetype = 2) +
    geom_errorbar(aes(xmin = ci_lower, xmax = ci_upper), width = 0.18,
                  orientation = "y") +
    geom_point(size = 2) +
    labs(x = "Shortage by time by post estimate (95% CI)", y = NULL,
         title = "Congenital syphilis CITS sensitivity analyses") +
    theme_bw(base_size = 11)
  ggsave(file.path(figures_dir, "cs_sensitivity_forest.png"), sensitivity_plot,
         width = 9, height = 6.5, dpi = 300)
  ggsave(file.path(figures_dir, "cs_sensitivity_forest.pdf"), sensitivity_plot,
         width = 9, height = 6.5)

  cat("Saved", nrow(result_table), "sensitivity estimates. Review convergence/singularity and samples.\n")
})
who_service_files <- file.path(raw_data_dir, "WHOData", c(
  "Women accessing antenatal care (ANC) services who were tested for syphilis(%), reported.csv",
  "Antenatal care attendees who were positive for syphilis(%), reported.csv",
  "Antenatal care attendees positive for syphilis who received treatment (%), reported.csv"
))


# -----------------------------------------------------------------------------
# WHO antenatal-service adjustment sensitivity
# -----------------------------------------------------------------------------
if (all(file.exists(who_service_files))) {


# -----------------------------------------------------------------------------
# WHO service-series audit and adjusted models
# -----------------------------------------------------------------------------
local({
  # Audit and describe WHO ANC service series; fit observed-data sensitivities.
  # Descriptive means use observed values in an unbalanced annual panel.
  # The adjusted models retain unbalanced country-year observations.

  library(lme4)
  library(lmerTest)
  library(readr)
  library(ggplot2)

  who_dir <- file.path(raw_data_dir, "WHOData")
  screen_file <- file.path(who_dir,
    "Women accessing antenatal care (ANC) services who were tested for syphilis(%), reported.csv")
  positive_file <- file.path(who_dir,
    "Antenatal care attendees who were positive for syphilis(%), reported.csv")
  treatment_file <- file.path(who_dir,
    "Antenatal care attendees positive for syphilis who received treatment (%), reported.csv")

  # 1. Read each WHO GHO extract and keep the latest country-year record.
  read_who_series <- function(file, expected_code, output_name) {
    raw <- readr::read_csv(file, show_col_types = FALSE,
                           locale = readr::locale(encoding = "UTF-8"),
                           name_repair = "minimal")
    raw <- raw[raw$IndicatorCode == expected_code & raw$`Location type` == "Country" &
                 as.integer(raw$Period) %in% study_years, ]
    out <- data.frame(
      location_iso = toupper(trimws(as.character(raw$SpatialDimValueCode))),
      year = as.integer(raw$Period),
      value = suppressWarnings(as.numeric(raw$FactValueNumeric)),
      modified = suppressWarnings(as.Date(substr(raw$DateModified, 1, 10))),
      stringsAsFactors = FALSE
    )
    out$modified[is.na(out$modified)] <- as.Date("1900-01-01")
    out <- out[order(out$location_iso, out$year, -as.numeric(out$modified)), ]
    out <- out[!duplicated(out[c("location_iso", "year")]), ]
    names(out)[names(out) == "value"] <- output_name
    out$modified <- NULL
    out
  }

  screening <- read_who_series(screen_file, "PerctestedANC", "screening_percent")
  positivity <- read_who_series(positive_file, "PercposANC", "positivity_percent")
  treatment <- read_who_series(treatment_file, "PerctreatedANC", "treatment_percent")

  # 2. Report missingness by year and describe the available observations.
  service_rows <- merge(screening, positivity, by = c("location_iso", "year"), all = TRUE)
  service_rows <- merge(service_rows, treatment, by = c("location_iso", "year"), all = TRUE)
  service_variables <- c("screening_percent", "positivity_percent", "treatment_percent")
  all_countries <- unique(toupper(trimws(as.character(cs_panel$location_iso))))
  missingness <- do.call(rbind, lapply(service_variables, function(variable) {
    observed_all_years <- service_rows[is.finite(service_rows[[variable]]),
                                      c("location_iso", "year", variable)]
    observed_counts <- table(observed_all_years$location_iso)
    complete_iso <- names(observed_counts)[observed_counts == length(study_years)]
    by_year <- lapply(study_years, function(year) {
      observed_year <- service_rows[service_rows$year == year,
                                   c("location_iso", variable)]
      d_year <- merge(data.frame(location_iso = all_countries), observed_year,
                      by = "location_iso", all.x = TRUE)
      d_complete <- d_year[d_year$location_iso %in% complete_iso, ]
      data.frame(variable = variable, year = year,
                 countries_in_GBD_panel = length(all_countries),
                 countries_observed = sum(is.finite(d_year[[variable]])),
                 countries_missing = length(all_countries) -
                   sum(is.finite(d_year[[variable]])),
                 countries_complete_series = length(complete_iso),
                 complete_series_mean_percent = if (nrow(d_complete))
                   mean(d_complete[[variable]], na.rm = TRUE) else NA_real_,
                 observed_mean_percent = mean(d_year[[variable]], na.rm = TRUE))
    })
    do.call(rbind, by_year)
  }))
  write.csv(missingness, file.path(results_dir, "who_anc_service_missingness.csv"),
            row.names = FALSE, na = "")

  plot_rows <- list()
  complete_country_summary <- list()
  for (variable in service_variables) {
    observed <- service_rows[is.finite(service_rows[[variable]]),
                             c("location_iso", "year", variable)]
    year_count <- table(observed$location_iso)
    complete_iso <- names(year_count)[year_count == length(study_years)]
    complete_country_summary[[variable]] <- data.frame(
      indicator = variable,
      countries_with_any_data = length(unique(observed$location_iso)),
      countries_complete_2010_2021 = length(complete_iso),
      complete_country_years = length(complete_iso) * length(study_years)
    )
    observed <- observed[observed$location_iso %in% all_countries, ]
    if (nrow(observed)) {
      names(observed)[names(observed) == variable] <- "percent"
      observed$indicator <- variable
      plot_rows[[variable]] <- observed
    }
  }
  write.csv(do.call(rbind, complete_country_summary),
            file.path(results_dir, "who_anc_complete_series_counts.csv"),
            row.names = FALSE, na = "")
  observed_plot_data <- do.call(rbind, plot_rows)
  if (!is.null(observed_plot_data) && nrow(observed_plot_data)) {
    observed_plot_data$indicator <- factor(observed_plot_data$indicator,
      levels = service_variables,
      labels = c("ANC syphilis screening", "ANC syphilis positivity",
                 "Syphilis treatment among ANC-positive attendees"))
    who_plot <- ggplot(observed_plot_data,
        aes(x = year, y = percent, group = location_iso)) +
      geom_line(colour = "grey70", alpha = 0.55) +
      stat_summary(aes(group = 1), fun = mean, geom = "line",
                   colour = "#0072B2", linewidth = 1.1) +
      facet_wrap(~ indicator, scales = "free_y") +
      labs(x = "Year", y = "Percent",
           title = "WHO ANC service indicators: available country-year observations",
           subtitle = "Blue line is the annual arithmetic mean; countries vary by year") +
      theme_bw(base_size = 10)
    ggsave(file.path(figures_dir, "who_anc_service_observed_series.png"),
           who_plot, width = 10, height = 7, dpi = 300)
  }

  # 3. Merge WHO screening and positivity with the GBD congenital-syphilis panel.
  d <- merge(cs_panel, screening, by = c("location_iso", "year"), all.x = TRUE)
  d <- merge(d, positivity, by = c("location_iso", "year"), all.x = TRUE)
  d <- d[d$year %in% study_years & d$shortage %in% c(0, 1), ]
  d$location_id <- factor(d$location_id)
  d$gbd_superregion <- factor(d$gbd_superregion)
  d$rate <- 1000 * d$val / d$population
  d$time_c <- d$year - break_year
  d$post <- as.integer(d$year >= break_year)
  d$birth_weight <- d$population / mean(d$population, na.rm = TRUE)
  d$hdi_z <- as.numeric(scale(d$hdi))
  d$sdi_z <- as.numeric(scale(d$sdi))
  d$gni_z <- as.numeric(scale(d$gni))
  d$female_syphilis_z <- as.numeric(scale(d$syphilis_prev_female))
  d$male_syphilis_z <- as.numeric(scale(d$syphilis_prev_male))
  d$abortion_rate_z <- as.numeric(scale(d$abor_misc_rate))
  d$anc4_z <- as.numeric(scale(d$pregnant_antenatal_care4))
  base_covariates <- c("hdi_z", "sdi_z", "gni_z", "gbd_superregion",
                       "female_syphilis_z", "male_syphilis_z", "abortion_rate_z")

  # 4. Fit base and service-adjusted LMMs on identical observed rows and weights.
  service_specs <- list(
    screening = "screening_percent",
    screening_and_positivity = c("screening_percent", "positivity_percent"),
    screening_and_anc4 = c("screening_percent", "anc4_z")
  )
  service_results <- list()
  for (spec_name in names(service_specs)) {
    service_covariates <- service_specs[[spec_name]]
    columns <- c("rate", "shortage", "time_c", "post", "population", "location_id",
                 base_covariates, service_covariates)
    model_data <- d[complete.cases(d[, columns]), ]
    model_data$location_id <- droplevels(model_data$location_id)
    model_data$gbd_superregion <- droplevels(model_data$gbd_superregion)
    model_data$birth_weight <- model_data$population / sum(model_data$population)
    base_formula <- as.formula(paste(
      "rate ~ time_c * shortage * post +", paste(base_covariates, collapse = " + "),
      "+ (1 + time_c | location_id)"
    ))
    adjusted_formula <- as.formula(paste(
      "rate ~ time_c * shortage * post +",
      paste(c(base_covariates, service_covariates), collapse = " + "),
      "+ (1 + time_c | location_id)"
    ))
    base_fit <- lmerTest::lmer(base_formula, data = model_data,
                               weights = birth_weight, REML = FALSE)
    adjusted_fit <- lmerTest::lmer(adjusted_formula, data = model_data,
                                   weights = birth_weight, REML = FALSE)
    for (fit_name in c("same_sample_base", "service_adjusted")) {
      fit <- if (fit_name == "same_sample_base") base_fit else adjusted_fit
      coefficient_table <- as.data.frame(summary(fit)$coefficients)
      for (target in c("shortage:post", "time_c:shortage:post")) {
      estimate <- coefficient_table[target, "Estimate"]
      se <- coefficient_table[target, "Std. Error"]
      df <- coefficient_table[target, "df"]
      service_results[[length(service_results) + 1L]] <- data.frame(
        sensitivity = spec_name, model = fit_name, term = target,
        contrast = if (target == "shortage:post") "Level change" else "Slope change/year",
        adjusted_variables = if (fit_name == "same_sample_base") "none" else
          paste(service_covariates, collapse = ";"),
        estimate = estimate, std_error = se, df = df,
        ci_lower = estimate - qt(0.975, df) * se,
        ci_upper = estimate + qt(0.975, df) * se,
        p_value = coefficient_table[target, "Pr(>|t|)"],
        country_years = nrow(model_data),
        countries = length(unique(model_data$location_id)),
        weight_sum = sum(model_data$birth_weight),
        singular = lme4::isSingular(fit), stringsAsFactors = FALSE
      )
      }
    }
  }
  service_results <- do.call(rbind, service_results)
  write.csv(service_results, file.path(results_dir, "who_service_sensitivity_lmm.csv"),
            row.names = FALSE, na = "")

  # 5. PMM sensitivity retaining all 89 classified countries. WHO screening and
  #    ANC4+ are imputed from pooled observed donors matched on year, region,
  #    shortage group, SDI, outcome rate and population. This is a separate MAR
  #    sensitivity because no country has a complete WHO series for all years.
  core_columns <- c("rate", "shortage", "time_c", "post", "population",
                    "location_id", base_covariates)
  full_pmm_data <- d[complete.cases(d[, core_columns]), ]
  full_pmm_data$location_id <- droplevels(full_pmm_data$location_id)
  full_pmm_data$gbd_superregion <- droplevels(full_pmm_data$gbd_superregion)
  full_pmm_data$log_population <- log(full_pmm_data$population)
  pmm_variables <- c("screening_percent", "pregnant_antenatal_care4")
  set.seed(20261011)
  number_of_imputations <- 50L
  service_imputations <- vector("list", number_of_imputations)
  imputed_cell_counts <- setNames(integer(length(pmm_variables)), pmm_variables)

  for (imputation in seq_len(number_of_imputations)) {
    completed <- full_pmm_data
    for (variable in pmm_variables) {
      observed_rows <- which(is.finite(completed[[variable]]))
      missing_rows <- which(!is.finite(completed[[variable]]))
      donor_data <- completed[observed_rows, ]
      donor_formula <- as.formula(paste(
        variable, "~ year + I(year^2) + shortage + sdi + rate + log_population +",
        "gbd_superregion"
      ))
      donor_model <- lm(donor_formula, data = donor_data)
      donor_prediction <- as.numeric(predict(donor_model, newdata = donor_data))
      target_prediction <- as.numeric(predict(donor_model,
                                               newdata = completed[missing_rows, ]))
      residual_sd <- sd(residuals(donor_model), na.rm = TRUE)
      if (!is.finite(residual_sd)) residual_sd <- 0
      for (row_index in seq_along(missing_rows)) {
        prediction_with_noise <- target_prediction[row_index] +
          rnorm(1, mean = 0, sd = residual_sd)
        nearest <- order(abs(donor_prediction - prediction_with_noise))
        donor_pool <- nearest[seq_len(min(5L, length(nearest)))]
        chosen <- sample(donor_pool, size = 1L)
        completed[[variable]][missing_rows[row_index]] <-
          completed[[variable]][observed_rows[chosen]]
        if (imputation == 1L) imputed_cell_counts[variable] <-
          imputed_cell_counts[variable] + 1L
      }
    }
    completed$screening_z <- as.numeric(scale(completed$screening_percent))
    completed$anc4_z <- as.numeric(scale(completed$pregnant_antenatal_care4))
    completed$birth_weight <- completed$population / sum(completed$population)
    service_imputations[[imputation]] <- completed
  }

  service_pmm_formulas <- list(
    same_sample_base = as.formula(paste(
      "rate ~ time_c * shortage * post +", paste(base_covariates, collapse = " + "),
      "+ (1 + time_c | location_id)"
    )),
    adjusted = as.formula(paste(
      "rate ~ time_c * shortage * post +",
      paste(c(base_covariates, "screening_z", "anc4_z"), collapse = " + "),
      "+ (1 + time_c | location_id)"
    ))
  )
  service_pmm_results <- list()
  for (model_label in names(service_pmm_formulas)) {
    for (target_term in c("post:shortage", "post:shortage:time_c")) {
    estimates <- standard_errors <- rep(NA_real_, number_of_imputations)
    for (imputation in seq_len(number_of_imputations)) {
      completed <- service_imputations[[imputation]]
      fit <- tryCatch(
        lmerTest::lmer(service_pmm_formulas[[model_label]], data = completed,
                       weights = birth_weight, REML = FALSE),
        error = function(e) e
      )
      if (inherits(fit, "error")) next
      coefficient_table <- as.data.frame(summary(fit)$coefficients)
      term_order <- vapply(strsplit(rownames(coefficient_table), ":", fixed = TRUE),
                           function(x) paste(sort(x), collapse = ":"), character(1))
      target <- rownames(coefficient_table)[term_order == target_term][1]
      estimates[imputation] <- coefficient_table[target, "Estimate"]
      standard_errors[imputation] <- coefficient_table[target, "Std. Error"]
    }
    valid <- is.finite(estimates) & is.finite(standard_errors)
    m <- sum(valid)
    if (m < 40L) {
      pooled_estimate <- pooled_se <- df <- lower <- upper <- p_value <- NA_real_
      status <- "not_estimable: fewer than 40 converged imputations"
    } else {
      pooled_estimate <- mean(estimates[valid])
      within_variance <- mean(standard_errors[valid]^2)
      between_variance <- var(estimates[valid])
      total_variance <- within_variance + (1 + 1 / m) * between_variance
      pooled_se <- sqrt(total_variance)
      if (between_variance < .Machine$double.eps) {
        df <- Inf
        critical <- qnorm(0.975)
        p_value <- 2 * pnorm(abs(pooled_estimate / pooled_se), lower.tail = FALSE)
      } else {
        relative_increase <- (1 + 1 / m) * between_variance / within_variance
        df <- (m - 1) * (1 + 1 / relative_increase)^2
        critical <- qt(0.975, df)
        p_value <- 2 * pt(abs(pooled_estimate / pooled_se), df,
                          lower.tail = FALSE)
      }
      lower <- pooled_estimate - critical * pooled_se
      upper <- pooled_estimate + critical * pooled_se
      status <- "estimated"
    }
    service_pmm_results[[paste(model_label, target_term)]] <- data.frame(
      specification = "screening_plus_ANC4_full_sample_PMM", model = model_label,
      term = if (target_term == "post:shortage") "shortage:post" else "time_c:shortage:post",
      contrast = if (target_term == "post:shortage") "Level change" else "Slope change/year",
      adjusted_variables = if (model_label == "adjusted")
        "WHO ANC syphilis screening; ANC4+" else "none",
      country_years = nrow(full_pmm_data),
      countries = length(unique(full_pmm_data$location_id)),
      shortage_countries = length(unique(full_pmm_data$location_id[
        full_pmm_data$shortage == 1
      ])),
      no_shortage_countries = length(unique(full_pmm_data$location_id[
        full_pmm_data$shortage == 0
      ])),
      imputations_converged = m, imputations_requested = number_of_imputations,
      estimate = pooled_estimate, std_error = pooled_se, df = df,
      ci_lower = lower, ci_upper = upper, p_value = p_value, status = status,
      imputation_method = "Pooled PMM; donor predictors include year, region, exposure, SDI, outcome and population",
      stringsAsFactors = FALSE
    )
    }
  }
  write.csv(do.call(rbind, service_pmm_results),
            file.path(results_dir, "who_service_full_sample_pmm_lmm.csv"),
            row.names = FALSE, na = "")
  write.csv(data.frame(
    variable = names(imputed_cell_counts),
    observed_cells = vapply(names(imputed_cell_counts), function(variable)
      sum(is.finite(full_pmm_data[[variable]])), integer(1)),
    imputed_cells_per_dataset = as.integer(imputed_cell_counts),
    country_years = nrow(full_pmm_data), countries = length(unique(full_pmm_data$location_id)),
    imputations = number_of_imputations,
    stringsAsFactors = FALSE
  ), file.path(results_dir, "who_service_pmm_imputation_summary.csv"),
  row.names = FALSE, na = "")

  cat("WHO missingness, available-year mean descriptions, and same-sample models saved.\n")
})

} else {
  message('WHO service-sensitivity section skipped because one or more WHO ANC source files are missing.')
}


# -----------------------------------------------------------------------------
# Deterministic quantitative-bias analysis
# -----------------------------------------------------------------------------
local({
  # Deterministic quantitative-bias calibration for the ordinary country LMM.
  # Annual covariates, live-birth precision weights, 2016 index; no bootstrap.
  # This file fits the covariate calibration models when explicitly run.

  suppressPackageStartupMessages(library(nlme))
  results_root <- results_dir
  output_dir <- results_dir
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # 1. Use the same eligible annual panel and primary adjustment variables.
  panel <- cs_panel
  panel$rate <- 1000 * panel$val / panel$population
  panel <- panel[!is.na(panel$shortage) & panel$year %in% 2010:2021 &
                   is.finite(panel$rate) & is.finite(panel$population) & panel$population > 0, ]
  annual_variables <- c("hdi", "sdi", "gni", "syphilis_prev_female",
                        "syphilis_prev_male", "abor_misc_rate")
  for (variable in annual_variables) {
    panel[[variable]] <- as.numeric(panel[[variable]])
    if (any(!is.finite(panel[[variable]]))) stop("Missing core covariate: ", variable)
    panel[[paste0(variable, "_z")]] <- as.numeric(scale(panel[[variable]]))
  }
  covariates <- paste0(annual_variables, "_z")
  panel$location_id <- factor(panel$location_id)
  panel$gbd_superregion <- factor(panel$gbd_superregion)
  panel$time_c <- panel$year - 2016L
  panel$post <- as.integer(panel$year >= 2016L)
  panel$population_relative <- panel$population / mean(panel$population)
  panel <- panel[order(panel$location_id, panel$year), ]

  # 2. Baseline country averages are used only to calibrate imbalance.
  # They never replace the annual values in the fitted regression models.
  baseline <- panel[panel$year %in% 2010:2013, ]
  country_means <- aggregate(baseline[, covariates],
                            baseline[, c("location_id", "shortage")], mean)
  imbalance <- data.frame(
    covariate = covariates,
    standardized_difference = vapply(covariates, function(variable) {
      mean(country_means[country_means$shortage == 1, variable]) -
        mean(country_means[country_means$shortage == 0, variable])
    }, numeric(1)))

  # 3. Fit one ordinary weighted association model for each calibration variable.
  association_rows <- fit_rows <- list()
  for (variable in covariates) {
    formula <- as.formula(paste(
      "rate ~ time_c * shortage * post +", paste(covariates, collapse = " + "),
      "+ gbd_superregion +", paste0(variable, ":post + ", variable, ":time_c:post")))
    warnings <- character()
    fit <- tryCatch(withCallingHandlers(
      nlme::lme(fixed = formula, random = ~ 1 + time_c | location_id,
                data = panel, correlation = NULL,
                weights = nlme::varFixed(~ I(1 / population_relative)),
                method = "ML", na.action = na.fail,
                control = nlme::lmeControl(opt = "nlminb", msMaxIter = 250,
                                          msMaxEval = 500, niterEM = 25)),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }), error = function(e) e)
    converged <- inherits(fit, "lme") && !any(grepl(
      "convergence error|iteration limit|singular convergence|false convergence",
      warnings, ignore.case = TRUE))
    fit_rows[[variable]] <- data.frame(
      covariate = variable, model = "LMM", converged = converged,
      warnings = paste(warnings, collapse = " | "),
      error = if (inherits(fit, "error")) conditionMessage(fit) else "")
    if (!converged) stop("Ordinary QBA association LMM did not converge: ", variable)
    coefficients <- summary(fit)$tTable
    normalized_terms <- vapply(strsplit(rownames(coefficients), ":", fixed = TRUE),
                               function(parts) paste(sort(parts), collapse = ":"), character(1))
    for (contrast in c("level", "slope")) {
      parts <- if (contrast == "level") c(variable, "post") else c(variable, "time_c", "post")
      position <- match(paste(sort(parts), collapse = ":"), normalized_terms)
      estimate <- coefficients[position, "Value"]
      se <- coefficients[position, "Std.Error"]
      df <- coefficients[position, "DF"]
      critical <- qt(0.975, df)
      association_rows[[length(association_rows) + 1L]] <- data.frame(
        covariate = variable, association = contrast, model = "LMM", estimate = estimate,
        se = se, df = df, ci_lower = estimate - critical * se,
        ci_upper = estimate + critical * se, p_value = coefficients[position, "p-value"])
    }

  }
  associations <- do.call(rbind, association_rows)

  # 4. Fit the ordinary nlme primary LMM to reproduce the designated calibration
  #    and its conditional t intervals. This does not overwrite the v3/lmer model.
  primary_fit <- nlme::lme(
    fixed = rate ~ time_c * shortage * post + hdi + sdi + gni +
      gbd_superregion + syphilis_prev_female + syphilis_prev_male + abor_misc_rate,
    random = ~ 1 + time_c | location_id, data = panel, correlation = NULL,
    weights = nlme::varFixed(~ I(1 / population_relative)), method = "ML",
    control = nlme::lmeControl(opt = "nlminb", msMaxIter = 250,
      msMaxEval = 500, niterEM = 25))
  primary <- as.data.frame(summary(primary_fit)$tTable)
  primary$term <- rownames(primary)
  primary$estimate <- primary$Value
  primary$ci_lower <- primary$Value - qt(0.975, primary$DF) * primary$Std.Error
  primary$ci_upper <- primary$Value + qt(0.975, primary$DF) * primary$Std.Error
  primary <- primary[primary$term %in% c("shortage:post", "time_c:shortage:post"), ]
  if (nrow(primary) != 2L) stop("Expected two ordinary primary LMM contrasts.")
  max_delta <- max(abs(imbalance$standardized_difference))
  bias_rows <- tipping_rows <- list()
  for (contrast in c("level", "slope")) {
    term <- if (contrast == "level") "shortage:post" else "time_c:shortage:post"
    effect <- primary[primary$term == term, ]
    max_gamma <- max(abs(associations$estimate[associations$association == contrast]))
    grid <- expand.grid(delta_multiple = seq(-2, 2, .5), gamma_multiple = seq(-2, 2, .5))
    grid$model <- "LMM"
    grid$association <- contrast
    grid$bias <- grid$delta_multiple * max_delta * grid$gamma_multiple * max_gamma
    grid$base_estimate <- effect$estimate
    grid$adjusted_estimate <- effect$estimate - grid$bias
    grid$ci_lower <- effect$ci_lower - grid$bias
    grid$ci_upper <- effect$ci_upper - grid$bias
    grid$interval_source <- "Conditional LMM interval shifted deterministically"
    bias_rows[[contrast]] <- grid
    tipping_rows[[contrast]] <- data.frame(
      model = "LMM", association = contrast, primary_estimate = effect$estimate,
      base_ci_lower = effect$ci_lower, base_ci_upper = effect$ci_upper,
      maximum_measured_imbalance = max_delta,
      maximum_measured_association = max_gamma,
      maximum_measured_bias = max_delta * max_gamma,
      multiple_for_point_null = abs(effect$estimate) / (max_delta * max_gamma),
      multiple_for_ci_include_zero = if (effect$ci_lower <= 0 & effect$ci_upper >= 0) 0 else min(abs(c(effect$ci_lower, effect$ci_upper))) / (max_delta * max_gamma))
  }

  # 5. Keep ordinary-LMM exports separate from historical paired-model outputs.
  write.csv(imbalance, file.path(output_dir, "qba_measured_imbalance.csv"), row.names = FALSE)
  write.csv(associations, file.path(output_dir, "qba_measured_associations.csv"), row.names = FALSE)
  write.csv(do.call(rbind, fit_rows), file.path(output_dir, "qba_fit_attempts.csv"), row.names = FALSE)
  write.csv(do.call(rbind, bias_rows), file.path(output_dir, "qba_tipping_point_grid.csv"), row.names = FALSE)
  write.csv(do.call(rbind, tipping_rows), file.path(output_dir, "qba_tipping_point_summary.csv"), row.names = FALSE)
})


# -----------------------------------------------------------------------------
# Alternative 2014 and 2015 interruption-year checks
# -----------------------------------------------------------------------------
local({
  # -----------------------------------------------------------------------------
  # Alternative 2014 and 2015 breakpoints
  # -----------------------------------------------------------------------------

  # Refit the 2014 and 2015 breakpoint models and calculate complete 95% CIs.
  # Run the primary-analysis section of 02_primary_and_subgroup_analysis.R first so the frozen sample is available.

  library(lme4)
  library(Matrix)

  data <- read.csv("results/original_primary_lmm_sample.csv")
  data$location_id <- factor(data$location_id)
  data$gbd_superregion <- factor(data$gbd_superregion)

  breakpoints <- c(2014L, 2015L)

  for (breakpoint in breakpoints) {
    data$time_c <- data$year - breakpoint
    data$post <- as.integer(data$year >= breakpoint)

    model <- lmer(
      rate ~ time_c * shortage * post + hdi + sdi + gni +
        gbd_superregion + syphilis_prev_female + syphilis_prev_male +
        abor_misc_rate + (1 + time_c | location_id),
      data = data,
      weights = birth_weight,
      REML = FALSE
    )

    # Include fixed effects and country-specific random effects in predictions.
    fixed_design <- getME(model, "X")
    random_design <- as.matrix(getME(model, "Z") %*% getME(model, "Lambda"))
    observed_design <- cbind(fixed_design, random_design)

    model_weights <- model.weights(model.frame(model))
    information <- crossprod(observed_design * sqrt(model_weights))
    n_fixed <- ncol(fixed_design)
    n_random <- ncol(random_design)
    random_rows <- (n_fixed + 1):(n_fixed + n_random)
    information[random_rows, random_rows] <-
      information[random_rows, random_rows] + diag(n_random)

    # Scale the information matrix before inversion to improve numerical stability.
    scale_terms <- sqrt(diag(information))
    scaled_information <- information / outer(scale_terms, scale_terms)
    full_covariance <- sigma(model)^2 * solve(scaled_information) /
      outer(scale_terms, scale_terms)

    if (max(abs(full_covariance[1:n_fixed, 1:n_fixed] - as.matrix(vcov(model)))) >
        1e-7) {
      stop("The reconstructed fixed-effect covariance does not match lme4.")
    }

    counterfactual_data <- model.frame(model)
    counterfactual_data$post <- 0L
    counterfactual_fixed_design <- model.matrix(
      delete.response(terms(model, fixed.only = TRUE)),
      counterfactual_data,
      contrasts.arg = attr(fixed_design, "contrasts")
    )
    counterfactual_fixed_design <- counterfactual_fixed_design[
      , colnames(fixed_design), drop = FALSE
    ]
    counterfactual_design <- cbind(counterfactual_fixed_design, random_design)

    model_parameters <- c(fixef(model), getME(model, "u"))
    prediction_rows <- list()

    for (year in 2010:2021) {
      for (shortage_status in 0:1) {
        rows <- which(data$year == year & data$shortage == shortage_status)
        if (length(rows) == 0L) next

        population_weights <- data$population[rows] / sum(data$population[rows])
        group_label <- if (shortage_status == 1L) {
          "Shortage reported"
        } else {
          "No shortage reported"
        }

        scenarios <- character()
        if (year <= breakpoint) {
          scenarios <- c(scenarios, "Pre-index continuation")
        }
        if (year >= breakpoint) {
          scenarios <- c(scenarios, "Counterfactual trajectory", "Fitted trajectory")
        }

        for (scenario_type in scenarios) {
          if (scenario_type %in% c("Pre-index continuation",
                                   "Counterfactual trajectory")) {
            prediction_design <- counterfactual_design
            scenario_label <- paste("Counterfactual post=0:",
                                    if (shortage_status == 1L) "shortage" else "non-shortage")
            segment_label <- scenario_type
          } else {
            prediction_design <- observed_design
            scenario_label <- paste("Observed:",
                                    if (shortage_status == 1L) "shortage" else "non-shortage")
            segment_label <- "Fitted trajectory"
          }

          average_design <- colSums(
            prediction_design[rows, , drop = FALSE] * population_weights
          )
          estimate <- sum(average_design * model_parameters)
          standard_error <- sqrt(as.numeric(
            t(average_design) %*% full_covariance %*% average_design
          ))

          prediction_rows[[length(prediction_rows) + 1L]] <- data.frame(
            year = year,
            scenario = scenario_label,
            estimate = estimate,
            se = standard_error,
            ci_lower = estimate - qnorm(0.975) * standard_error,
            ci_upper = estimate + qnorm(0.975) * standard_error,
            group = group_label,
            scenario_type = if (scenario_type == "Pre-index continuation") {
              "Fitted trajectory"
            } else {
              scenario_type
            },
            segment = segment_label,
            stringsAsFactors = FALSE
          )
        }
      }
    }

    segments <- do.call(rbind, prediction_rows)
    output_number <- if (breakpoint == 2014L) "s10" else "s11"
    output_file <- file.path("results", paste0("figure_", output_number,
                                                 "_covariance_segments.csv"))

    # Compare point estimates with the previous v7 run before replacing it.
    if (file.exists(output_file)) {
      previous_segments <- read.csv(output_file, stringsAsFactors = FALSE)
      comparison <- merge(
        segments[, c("year", "scenario", "segment", "estimate")],
        previous_segments[, c("year", "scenario", "segment", "estimate")],
        by = c("year", "scenario", "segment"),
        suffixes = c("_new", "_previous")
      )
      comparison$point_difference <- comparison$estimate_new -
        comparison$estimate_previous
      write.csv(comparison,
        file.path("logs", paste0("figure_", output_number,
                                 "_point_reproduction_check.csv")),
        row.names = FALSE
      )
      if (nrow(comparison) != nrow(segments) ||
          max(abs(comparison$point_difference)) > 1e-5) {
        print(comparison[order(abs(comparison$point_difference), decreasing = TRUE), ][1:8, ])
        stop("Alternative-breakpoint point estimates differ from the previous v7 output. See the point-reproduction check in logs/.")
      }
    }

    write.csv(segments, output_file, row.names = FALSE)
  }
})


# -----------------------------------------------------------------------------
# Sensitivity curves and alternative-breakpoint figures
# -----------------------------------------------------------------------------
local({
  # S10/S11: fit and plot the 2014/2015 sensitivity trajectories from local analysis outputs.
  primary_data <- read.csv(file.path(results_dir, "original_primary_lmm_sample.csv"))
  primary_data$group <- ifelse(primary_data$shortage == 1,
    "Shortage reported", "No shortage reported")
  observed <- aggregate(
    cbind(val, population) ~ year + group,
    data = primary_data,
    FUN = sum
  )
  observed$rate <- 1000 * observed$val / observed$population

  for (index in 1:2) {
    source_number <- if (index == 1) "S08" else "S09"
    current_number <- if (index == 1) "S10" else "S11"
    intervention <- if (index == 1) 2014 else 2015
    segments <- read.csv(file.path(results_dir,
      paste0("figure_", tolower(current_number), "_covariance_segments.csv")))
    counterfactual_at_break <- segments[
      segments$year == intervention & segments$segment == "Counterfactual trajectory",
      c("group", "estimate")
    ]
    fitted_at_break <- segments[
      segments$year == intervention & segments$segment == "Fitted trajectory",
      c("group", "estimate")
    ]
    names(counterfactual_at_break)[2] <- "before"
    names(fitted_at_break)[2] <- "after"
    jumps <- merge(counterfactual_at_break, fitted_at_break, by = "group")
    fitted <- ggplot() +
      geom_ribbon(data = segments, aes(x = year, ymin = ci_lower, ymax = ci_upper,
        fill = group, group = interaction(group, scenario_type, segment)), alpha = 0.14) +
      geom_vline(xintercept = intervention, linetype = "dashed", colour = "grey50") +
      geom_segment(data = jumps, aes(x = intervention, xend = intervention,
        y = before, yend = after, colour = group), linewidth = 0.8) +
      geom_line(data = segments, aes(x = year, y = estimate, colour = group,
        linetype = scenario_type, group = interaction(group, scenario_type, segment)), linewidth = 0.8) +
      geom_point(data = observed, aes(x = year, y = rate, colour = group), size = 1.8) +
      scale_colour_manual(values = shortage_colours) +
      scale_fill_manual(values = shortage_colours) +
      scale_linetype_manual(values = c("Fitted trajectory" = "solid", "Counterfactual trajectory" = "longdash")) +
      scale_x_continuous(breaks = study_years) +
      scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
      labs(x = "Year", y = "Congenital syphilis cases per 1,000 live births",
           colour = NULL, fill = NULL, linetype = NULL) + publication_theme()
    ggsave(file.path(figures_dir, paste0("Figure_", current_number, "_CITS.png")),
           fitted, width = 9, height = 5.5, dpi = 300)
    ggsave(file.path(figures_dir, paste0("Figure_", current_number, "_CITS.pdf")),
           fitted, width = 9, height = 5.5, device = grDevices::cairo_pdf)
  }
})


# -----------------------------------------------------------------------------
# Country-label permutation analysis and Figures S12–S13
# -----------------------------------------------------------------------------
local({
  # Permute shortage labels across countries while preserving the 38/51 group sizes.
  # Run from the analysis root after the original primary model has been fitted.

  library(lme4)
  library(lmerTest)
  library(ggplot2)

  arguments <- commandArgs(trailingOnly = TRUE)
  n_permutations <- if (length(arguments) >= 1L) as.integer(arguments[1]) else 500L
  random_seed <- if (length(arguments) >= 2L) as.integer(arguments[2]) else 1L
  if (!is.finite(n_permutations) || n_permutations < 1L) stop("Permutation count must be positive.")

  # Prepare the same 89-country, 1,068-row model sample as the primary model.
  d <- cs_panel[cs_panel$year %in% study_years & cs_panel$shortage %in% c(0, 1), ]
  d$location_id <- factor(d$location_id)
  d$gbd_superregion <- factor(d$gbd_superregion)
  d$time_c <- d$year - break_year
  d$post <- as.integer(d$year >= break_year)
  d$rate <- 1000 * d$val / d$population
  covariates <- c("hdi", "sdi", "gni", "gbd_superregion",
                  "syphilis_prev_female", "syphilis_prev_male", "abor_misc_rate")
  required <- c("rate", "population", "shortage", "location_id", "year", covariates)
  d <- d[complete.cases(d[, required]), ]
  d$location_id <- droplevels(d$location_id)
  d$gbd_superregion <- droplevels(d$gbd_superregion)
  d$birth_weight <- d$population / sum(d$population)

  primary_file <- file.path(results_dir, "original_primary_lmm_coefficients.csv")
  if (!file.exists(primary_file)) stop("Run the primary-analysis section of 02_primary_and_subgroup_analysis.R first.")
  primary <- read.csv(primary_file, stringsAsFactors = FALSE)
  observed_slope <- primary$Estimate[primary$term == "time_c:shortage:post"]
  observed_level <- primary$Estimate[primary$term == "shortage:post"]

  country_labels <- unique(d[c("location_id", "shortage")])
  if (nrow(country_labels) != 89L || sum(country_labels$shortage == 1) != 38L ||
      sum(country_labels$shortage == 0) != 51L) {
    stop("Permutation requires the frozen 38 shortage / 51 comparison country groups.")
  }

  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  permuted_rows <- vector("list", n_permutations)
  for (iteration in seq_len(n_permutations)) {
    # Match the author-designated fixed per-replicate seeds 1,...,500.
    set.seed(random_seed + iteration - 1L)
    shuffled_labels <- sample(country_labels$shortage, replace = FALSE)
    names(shuffled_labels) <- as.character(country_labels$location_id)
    d$permuted_shortage <- unname(shuffled_labels[as.character(d$location_id)])
    permutation_formula <- as.formula(paste(
      "rate ~ time_c * permuted_shortage * post +",
      paste(covariates, collapse = " + "), "+ (1 + time_c | location_id)"
    ))
    fit <- tryCatch(
      suppressWarnings(lmerTest::lmer(
        permutation_formula, data = d, weights = birth_weight, REML = FALSE,
        control = lme4::lmerControl(optimizer = "nloptwrap",
                                  optCtrl = list(maxeval = 200000L))
      )), error = function(e) e
    )
    if (inherits(fit, "error")) {
      permuted_rows[[iteration]] <- data.frame(
        iteration = iteration, slope_estimate = NA_real_, level_estimate = NA_real_,
        singular = NA, status = conditionMessage(fit)
      )
      next
    }
    coefficient_table <- as.data.frame(summary(fit)$coefficients)
    normalized <- vapply(strsplit(rownames(coefficient_table), ":", fixed = TRUE),
                         function(x) paste(sort(x), collapse = ":"), character(1))
    slope_term <- rownames(coefficient_table)[normalized ==
      "permuted_shortage:post:time_c"][1]
    level_term <- rownames(coefficient_table)[normalized ==
      "permuted_shortage:post"][1]
    permuted_rows[[iteration]] <- data.frame(
      iteration = iteration,
      slope_estimate = coefficient_table[slope_term, "Estimate"],
      level_estimate = coefficient_table[level_term, "Estimate"],
      singular = lme4::isSingular(fit), status = "estimated"
    )
  }
  permutation_results <- do.call(rbind, permuted_rows)
  write.csv(permutation_results,
            file.path(results_dir, "country_label_permutation_estimates.csv"),
            row.names = FALSE, na = "")

  successful <- permutation_results[permutation_results$status == "estimated", ]
  summary <- data.frame(
    permutations_requested = n_permutations,
    permutations_estimated = nrow(successful),
    random_seed = random_seed,
    observed_slope = observed_slope,
    slope_two_sided_empirical_p =
      (1 + sum(abs(successful$slope_estimate) >= abs(observed_slope))) /
        (nrow(successful) + 1),
    observed_level = observed_level,
    level_two_sided_empirical_p =
      (1 + sum(abs(successful$level_estimate) >= abs(observed_level))) /
        (nrow(successful) + 1),
    stringsAsFactors = FALSE
  )
  write.csv(summary, file.path(results_dir, "country_label_permutation_summary.csv"),
            row.names = FALSE, na = "")


  cat("Saved", nrow(successful), "successful country-label permutations out of",
      n_permutations, "requested.\n")


  # -----------------------------------------------------------------------------
  # Permutation figures
  # -----------------------------------------------------------------------------

  # Restore the original v3 figure style using the current saved permutations.
  # Run from the analysis root after the permutation section of 03_sensitivity_analyses.R; no model fitting here.
  library(ggplot2)
  draws <- read.csv("results/country_label_permutation_estimates.csv")
  draws <- draws[draws$status == "estimated", ]
  test <- read.csv("results/country_label_permutation_summary.csv")
  original_theme <- theme_minimal(base_family = "Arial") +
    theme(text = element_text(colour = "black", size = 10),
      line = element_line(colour = "black", linewidth = 0.25),
      axis.line = element_line(colour = "black"),
      axis.ticks = element_line(colour = "black"),
      axis.title = element_text(face = "bold"),
      axis.title.x = element_text(margin = margin(t = 10)),
      axis.title.y = element_text(margin = margin(r = 10)),
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold", size = 12, hjust = 0),
      plot.subtitle = element_text(face = "italic", hjust = 0),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA))

  # S12: level change. Red dashed = observed; blue dashed = permutation mean.
  level_mean <- mean(draws$level_estimate)
  level_plot <- ggplot(draws, aes(x = level_estimate)) +
    geom_histogram(bins = 50, fill = "grey80", colour = "black") +
    geom_vline(xintercept = test$observed_level, colour = "#CD534CFF",
               linetype = "dashed", linewidth = 1) +
    geom_vline(xintercept = level_mean, colour = "blue",
               linetype = "dashed", linewidth = 1) +
    labs(title = "Permutation Test Distribution",
         subtitle = sprintf("Two-sided empirical P = %.3f", test$level_two_sided_empirical_p),
         x = "Level change", y = "Frequency",
         caption = "Red dashed: estimated effect; blue dashed: permutation average effect") +
    original_theme
  level_top <- max(ggplot_build(level_plot)$data[[1]]$count) * 1.1
  level_plot <- level_plot +
    annotate("text", x = test$observed_level - 0.005, y = level_top,
             label = "Estimated effect", hjust = 1, colour = "#CD534CFF") +
    annotate("text", x = level_mean + 0.005, y = level_top,
             label = "Permutation average effect", hjust = 0, colour = "blue")
  ggsave("figures/Figure_S12_permutation_level.png", level_plot,
         width = 9, height = 6, dpi = 300)
  ggsave("figures/Figure_S12_permutation_level.pdf", device = grDevices::cairo_pdf, level_plot, width = 9, height = 6)

  # S13: slope change. Display follows the original v3 style; P values use the
  # corrected (exceedances + 1)/(valid permutations + 1) two-sided calculation.
  slope_mean <- mean(draws$slope_estimate)
  slope_plot <- ggplot(draws, aes(x = slope_estimate)) +
    geom_histogram(bins = 50, fill = "grey80", colour = "black") +
    geom_vline(xintercept = test$observed_slope, colour = "#CD534CFF",
               linetype = "dashed", linewidth = 1) +
    geom_vline(xintercept = slope_mean, colour = "blue",
               linetype = "dashed", linewidth = 1) +
    labs(title = "Permutation Test Distribution",
         subtitle = sprintf("Two-sided empirical P = %.3f", test$slope_two_sided_empirical_p),
         x = "Slope change", y = "Frequency",
         caption = "Red dashed: estimated effect; blue dashed: permutation average effect") +
    original_theme
  slope_top <- max(ggplot_build(slope_plot)$data[[1]]$count) * 1.1
  slope_plot <- slope_plot +
    annotate("text", x = test$observed_slope + 0.005, y = slope_top,
             label = "Estimated effect", hjust = 0, colour = "#CD534CFF") +
    annotate("text", x = slope_mean + 0.005, y = slope_top,
             label = "Permutation average effect", hjust = 0, colour = "blue")
  ggsave("figures/Figure_S13_permutation_slope.png", slope_plot,
         width = 9, height = 6, dpi = 300)
  ggsave("figures/Figure_S13_permutation_slope.pdf", device = grDevices::cairo_pdf, slope_plot, width = 9, height = 6)
})


# -----------------------------------------------------------------------------
# Quantitative-bias calibration heatmap
# -----------------------------------------------------------------------------
local({
  # Supplementary Figure S17: calibrated level and slope changes over QBA scenarios.
  # Run this file from the v7 analysis root after the QBA results are available.
  # The heatmap uses saved deterministic QBA estimates and does not refit models.

  # 1. Load packages ----------------------------------------------------------
  suppressPackageStartupMessages(library(ggplot2))
  suppressPackageStartupMessages(library(patchwork))

  # 2. Read the saved QBA calibration grid -----------------------------------
  qba_dir <- results_dir
  bias_grid <- read.csv(
      file.path(qba_dir, "qba_tipping_point_grid.csv"),
      stringsAsFactors = FALSE
  )

  level_data <- bias_grid[
      bias_grid$model == "LMM" & bias_grid$association == "level",
  ]
  slope_data <- bias_grid[
      bias_grid$model == "LMM" & bias_grid$association == "slope",
  ]

  stopifnot(
      nrow(level_data) == 81L,
      nrow(slope_data) == 81L,
      all(is.finite(level_data$adjusted_estimate)),
      all(is.finite(slope_data$adjusted_estimate)),
      all(level_data$delta_multiple %in% seq(-2, 2, by = 0.5)),
      all(level_data$gamma_multiple %in% seq(-2, 2, by = 0.5))
  )

  # Use the same colour range within each panel, centred on zero.
  level_limit <- max(abs(level_data$adjusted_estimate))
  slope_limit <- max(abs(slope_data$adjusted_estimate))

  # 3. Plot the calibrated level change --------------------------------------
  level_plot <- ggplot(
      level_data,
      aes(x = delta_multiple, y = gamma_multiple, fill = adjusted_estimate)
  ) +
      geom_tile(width = 0.5, height = 0.5) +
      geom_contour(
          data = level_data,
          aes(x = delta_multiple, y = gamma_multiple, z = adjusted_estimate),
          breaks = 0,
          colour = "#20252B",
          linewidth = 0.65,
          inherit.aes = FALSE
      ) +
      scale_x_continuous(breaks = -2:2, limits = c(-2.25, 2.25), expand = c(0, 0)) +
      scale_y_continuous(breaks = -2:2, limits = c(-2.25, 2.25), expand = c(0, 0)) +
      scale_fill_gradient2(
          low = "#2878A1",
          mid = "#F7F7F7",
          high = "#D94B4B",
          midpoint = 0,
          limits = c(-level_limit, level_limit),
          name = "Calibrated level change\n(cases per 1,000 live births)"
      ) +
      coord_equal() +
      labs(
          title = "A  Level change",
          x = "Unmeasured-factor imbalance\n(multiple of measured maximum)",
          y = "Outcome association\n(multiple of measured maximum)"
      ) +
      guides(fill = guide_colourbar(
          title.position = "top",
          title.hjust = 0,
          barwidth = grid::unit(4.2, "cm"),
          barheight = grid::unit(0.28, "cm")
      )) +
      publication_theme(base_size = 11) +
      theme(
          axis.text = element_text(colour = "black"),
          axis.title = element_text(colour = "black"),
          plot.title = element_text(face = "bold", size = 13),
          legend.position = "bottom",
          legend.title = element_text(size = 9, face = "plain"),
          legend.text = element_text(size = 9),
          plot.margin = margin(8, 8, 8, 8)
      )

  # 4. Plot the calibrated slope change --------------------------------------
  slope_plot <- ggplot(
      slope_data,
      aes(x = delta_multiple, y = gamma_multiple, fill = adjusted_estimate)
  ) +
      geom_tile(width = 0.5, height = 0.5) +
      scale_x_continuous(breaks = -2:2, limits = c(-2.25, 2.25), expand = c(0, 0)) +
      scale_y_continuous(breaks = -2:2, limits = c(-2.25, 2.25), expand = c(0, 0)) +
      scale_fill_gradient2(
          low = "#2878A1",
          mid = "#F7F7F7",
          high = "#D94B4B",
          midpoint = 0,
          limits = c(-slope_limit, slope_limit),
          name = "Calibrated slope change\n(cases per 1,000 live births/year)"
      ) +
      coord_equal() +
      labs(
          title = "B  Slope change",
          x = "Unmeasured-factor imbalance\n(multiple of measured maximum)",
          y = "Outcome association\n(multiple of measured maximum)"
      ) +
      guides(fill = guide_colourbar(
          title.position = "top",
          title.hjust = 0,
          barwidth = grid::unit(4.2, "cm"),
          barheight = grid::unit(0.28, "cm")
      )) +
      publication_theme(base_size = 11) +
      theme(
          axis.text = element_text(colour = "black"),
          axis.title = element_text(colour = "black"),
          plot.title = element_text(face = "bold", size = 13),
          legend.position = "bottom",
          legend.title = element_text(size = 9, face = "plain"),
          legend.text = element_text(size = 9),
          plot.margin = margin(8, 8, 8, 8)
      )

  # 5. Combine panels and save the figure ------------------------------------
  qba_heatmap <- level_plot + slope_plot + plot_layout(nrow = 1, widths = c(1, 1))

  output_dir <- figures_dir
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  ggsave(
      file.path(output_dir, "Figure_S17_QBA_calibration_heatmap.png"),
      qba_heatmap,
      width = 10.5,
      height = 6,
      dpi = 600,
      bg = "white"
  )
  ggsave(
      file.path(output_dir, "Figure_S17_QBA_calibration_heatmap.pdf"),
      qba_heatmap,
      width = 10.5,
      height = 6,
      device = grDevices::cairo_pdf,
      bg = "white"
  )

  write.csv(
      rbind(level_data, slope_data),
      file.path(qba_dir, "qba_heatmap_plot_data.csv"),
      row.names = FALSE
  )
})
