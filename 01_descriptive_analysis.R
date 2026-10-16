# 01_descriptive_analysis.R
# Country and annual summaries for the congenital-syphilis analysis.
# Run from the v7 analysis root. Review the numbered sections and outputs.
# Sections: 1) input audit; 2) annual rates by shortage-report status;
# 3) included versus unclassified country profile.

# Purpose: Check the source files and assembled panels before any model is run.
# Run the file from the v7 analysis root to check inputs and produce summaries.

# 1. Choose the folder that contains the assembled analysis data.
# Set BPG_RAW_DATA_DIR and BPG_PROCESSED_DATA_DIR for external locations.
# Raw exports and prepared analysis inputs are stored separately.
raw_data_dir <- Sys.getenv("BPG_RAW_DATA_DIR", unset = file.path("data", "raw_data"))
processed_data_dir <- Sys.getenv("BPG_PROCESSED_DATA_DIR", unset = file.path("data", "processed_data"))

if (!file.exists(file.path(processed_data_dir, "data_combined.Rdata")) ||
    !file.exists(file.path(processed_data_dir, "data_combined_hiv.Rdata"))) {
  stop("Prepared panels are required in data/processed_data/. Raw exports belong in data/raw_data/.")
}

log_dir <- file.path("logs")
dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)

# 2. Check core files separately from files used only by optional modules.
core_files <- c(
  file.path(processed_data_dir, "data_combined.Rdata"),
  file.path(processed_data_dir, "data_combined_hiv.Rdata")
)

iqvia_raw_dir <- file.path(raw_data_dir, "IQVIA")
iqvia_processed_dir <- file.path(processed_data_dir, "IQVIA")

optional_paths <- c(
  "GBD extract folder" = file.path(raw_data_dir, "IHME-GBD_2021_DATA"),
  "WHO source-data folder" = file.path(raw_data_dir, "WHOData"),
  "UNAIDS source-data folder" = file.path(raw_data_dir, "AIDSInfo"),
  "UN population source-data folder" = file.path(raw_data_dir, "UN_WPP"),
  "Licensed IQVIA raw exports" = iqvia_raw_dir,
  "IQVIA curated inputs" = iqvia_processed_dir,
  "Prepared HIV model panels" = file.path(processed_data_dir, "hiv")
)
optional_paths <- optional_paths[nzchar(optional_paths)]

input_paths <- c(core_files, optional_paths)
input_names <- c(basename(core_files), names(optional_paths))
file_status <- file.exists(input_paths) | dir.exists(input_paths)
is_core_required <- c(rep(TRUE, length(core_files)),
                      rep(FALSE, length(optional_paths)))
input_audit <- data.frame(
  input = input_names,
  input_path = input_paths,
  available = file_status,
  required_for_core_analyses = is_core_required,
  stringsAsFactors = FALSE
)
write.csv(input_audit, file.path(log_dir, "input_audit.csv"), row.names = FALSE, na = "")
print(input_audit, row.names = FALSE)

missing_core_files <- core_files[!file_status[seq_along(core_files)]]
if (length(missing_core_files) > 0) {
  stop("One or more assembled analysis panels are missing. See logs/input_audit.csv.")
}

# 3. Load the two assembled panels in separate environments.
cs_environment <- new.env(parent = emptyenv())
load(file.path(processed_data_dir, "data_combined.Rdata"), envir = cs_environment)
cs_panel <- get("data_combined", envir = cs_environment)

hiv_environment <- new.env(parent = emptyenv())
load(file.path(processed_data_dir, "data_combined_hiv.Rdata"), envir = hiv_environment)
hiv_panel <- get("data_combined_hiv", envir = hiv_environment)

# 4. Check required columns and panel keys.
cs_columns <- c(
  "location_id", "location_iso", "location_name", "year", "val",
  "population", "shortage", "hdi", "sdi", "gni", "gbd_superregion"
)
hiv_columns <- c(
  "location_id", "location_iso", "location_name", "year", "val",
  "population", "shortage", "hdi", "sdi", "gni", "gbd_superregion"
)

missing_cs_columns <- setdiff(cs_columns, names(cs_panel))
missing_hiv_columns <- setdiff(hiv_columns, names(hiv_panel))
if (length(missing_cs_columns) > 0) {
  stop("CS panel is missing columns: ", paste(missing_cs_columns, collapse = ", "))
}
if (length(missing_hiv_columns) > 0) {
  stop("HIV panel is missing columns: ", paste(missing_hiv_columns, collapse = ", "))
}

cs_panel$location_id <- as.character(cs_panel$location_id)
hiv_panel$location_id <- as.character(hiv_panel$location_id)

cs_duplicate_keys <- duplicated(cs_panel[c("location_id", "year")])
hiv_duplicate_keys <- duplicated(hiv_panel[c("location_id", "year")])

if (any(cs_duplicate_keys)) stop("Duplicate location-year keys in the CS panel.")
if (any(hiv_duplicate_keys)) stop("Duplicate location-year keys in the HIV panel.")

if (any(!is.finite(cs_panel$population) | cs_panel$population <= 0)) {
  stop("CS panel has a missing, non-finite, or non-positive Birth population.")
}
if (any(!is.finite(hiv_panel$population) | hiv_panel$population <= 0)) {
  stop("HIV panel has a missing, non-finite, or non-positive Birth-age population.")
}

# 5. Save a concise record of the assembled panels.
cs_shortage_countries <- unique(cs_panel$location_id[cs_panel$shortage %in% c(0, 1)])
cs_shortage_only <- cs_panel[cs_panel$shortage %in% c(0, 1), ]
hiv_shortage_only <- hiv_panel[hiv_panel$shortage %in% c(0, 1), ]

panel_audit <- data.frame(
  panel = c("congenital_syphilis", "hiv_birth_age"),
  rows = c(nrow(cs_panel), nrow(hiv_panel)),
  locations = c(length(unique(cs_panel$location_id)), length(unique(hiv_panel$location_id))),
  first_year = c(min(cs_panel$year), min(hiv_panel$year)),
  last_year = c(max(cs_panel$year), max(hiv_panel$year)),
  classified_countries = c(length(cs_shortage_countries),
                           length(unique(hiv_shortage_only$location_id))),
  shortage_countries = c(length(unique(cs_shortage_only$location_id[cs_shortage_only$shortage == 1])),
                         length(unique(hiv_shortage_only$location_id[hiv_shortage_only$shortage == 1]))),
  no_shortage_countries = c(length(unique(cs_shortage_only$location_id[cs_shortage_only$shortage == 0])),
                            length(unique(hiv_shortage_only$location_id[hiv_shortage_only$shortage == 0])))
)

write.csv(panel_audit, file.path(log_dir, "panel_audit.csv"), row.names = FALSE)
print(panel_audit, row.names = FALSE)

# 6. Record the R runtime and versions of the analysis packages.
analysis_packages <- c(
  "dagitty", "dplyr", "forestploter", "ggplot2", "ggsci", "lme4",
  "lmerTest", "Matrix", "nlme", "openxlsx", "patchwork", "readr",
  "scales", "tidyr"
)
package_versions <- character(length(analysis_packages))
for (i in seq_along(analysis_packages)) {
  package_name <- analysis_packages[i]
  if (requireNamespace(package_name, quietly = TRUE)) {
    package_versions[i] <- as.character(utils::packageVersion(package_name))
  } else {
    package_versions[i] <- "not installed"
  }
}
names(package_versions) <- analysis_packages
session_record <- c(
  capture.output(sessionInfo()),
  "",
  "Analysis package versions:",
  paste(names(package_versions), package_versions, sep = ": ")
)
writeLines(session_record, file.path(log_dir, "session_info.txt"))

cat("Data checks completed. Review logs/input_audit.csv and logs/panel_audit.csv.\n")

# Create the standard output folders and analysis constants.
results_dir <- "results"
figures_dir <- "figures"
tables_dir <- "tables"
logs_dir <- "logs"
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(logs_dir, showWarnings = FALSE, recursive = TRUE)
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
# Annual rates by shortage-response group and country-selection profile
# -----------------------------------------------------------------------------
local({
  # Descriptive trends and country-selection summaries
  # Describe annual congenital-syphilis incidence by the three shortage groups.

  library(ggplot2)

  d <- cs_panel[cs_panel$year %in% study_years, ]
  d$rate_per_1000_births <- 1000 * d$val / d$population
  d$shortage_group <- ifelse(is.na(d$shortage), "Unclassified",
                             ifelse(d$shortage == 1, "Shortage reported",
                                    "No shortage reported"))

  # Summarise counts and denominators first; the population-weighted rate is the
  # sum of cases divided by the sum of Birth population, not a mean of rates.
  group_years <- split(d, list(d$year, d$shortage_group), drop = TRUE)
  summary_rows <- lapply(group_years, function(x) {
    data.frame(
      year = x$year[1], group = x$shortage_group[1],
      countries = length(unique(x$location_id)),
      congenital_syphilis_cases = sum(x$val, na.rm = TRUE),
      live_birth_population = sum(x$population, na.rm = TRUE),
      population_weighted_rate = 1000 * sum(x$val, na.rm = TRUE) /
        sum(x$population, na.rm = TRUE),
      arithmetic_mean_country_rate = mean(x$rate_per_1000_births, na.rm = TRUE),
      sd_country_rate = sd(x$rate_per_1000_births, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  annual_summary <- do.call(rbind, summary_rows)
  annual_summary <- annual_summary[order(annual_summary$group, annual_summary$year), ]
  write.csv(annual_summary, file.path(results_dir, "cs_annual_group_summary.csv"),
            row.names = FALSE, na = "")

  trend_plot <- ggplot(annual_summary,
                       aes(x = year, y = population_weighted_rate,
                           colour = group, group = group)) +
    scale_colour_manual(values = shortage_colours) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.8) +
    geom_vline(xintercept = break_year, linetype = 2, colour = "grey45") +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
    labs(x = "Year", y = "Congenital syphilis cases per 1,000 live births",
         colour = NULL,
         title = "Congenital syphilis incidence by BPG shortage-reporting group",
         subtitle = "Annual rates weighted by live births. Dashed line marks 2016") +
    publication_theme(base_size = 11) +
    theme(legend.position = "bottom")
  ggsave(file.path(figures_dir, "cs_three_group_annual_rates.png"),
         trend_plot, width = 9, height = 5.5, dpi = 300)
  ggsave(file.path(figures_dir, "cs_three_group_annual_rates.pdf"), device = grDevices::cairo_pdf,
         trend_plot, width = 9, height = 5.5)
  ggsave(file.path(figures_dir, "Figure_S06_Three_country_groups.png"),
         trend_plot, width = 9, height = 5.5, dpi = 300)
  ggsave(file.path(figures_dir, "Figure_S06_Three_country_groups.pdf"), device = grDevices::cairo_pdf,
         trend_plot, width = 9, height = 5.5)

  # Compare 2010–2021 represented versus unclassified countries at the country level.
  country_summary <- aggregate(
    cbind(rate_per_1000_births, sdi, hdi, gni, syphilis_prev_female) ~
      location_id + location_name + shortage_group,
    data = d[d$year == 2010, ], FUN = function(x) x[1]
  )
  write.csv(country_summary,
            file.path(results_dir, "cs_country_2010_characteristics.csv"),
            row.names = FALSE, na = "")

  cat("Saved three-group annual summaries and the trend figure.\n")


  # -----------------------------------------------------------------------------
  # Country-selection profile
  # -----------------------------------------------------------------------------

  # Compare countries with classified and unclassified WHO shortage responses.


  d <- cs_panel[cs_panel$year == 2021, ]
  d$selection_group <- ifelse(d$shortage %in% c(0, 1), "Included in 89-country CITS",
                              "Unclassified / excluded from CITS")
  d$rate_per_1000_births <- 1000 * d$val / d$population

  # Country-level summaries; categorical distributions are also exported separately.
  groups <- split(d, d$selection_group)
  continuous <- do.call(rbind, lapply(groups, function(x) data.frame(
    group = x$selection_group[1], countries = length(unique(x$location_id)),
    mean_sdi = mean(x$sdi, na.rm = TRUE), sd_sdi = sd(x$sdi, na.rm = TRUE),
    mean_hdi = mean(x$hdi, na.rm = TRUE), sd_hdi = sd(x$hdi, na.rm = TRUE),
    mean_gni = mean(x$gni, na.rm = TRUE), sd_gni = sd(x$gni, na.rm = TRUE),
    mean_cs_rate = mean(x$rate_per_1000_births, na.rm = TRUE),
    sd_cs_rate = sd(x$rate_per_1000_births, na.rm = TRUE),
    mean_female_syphilis_prevalence = mean(x$syphilis_prev_female, na.rm = TRUE),
    sd_female_syphilis_prevalence = sd(x$syphilis_prev_female, na.rm = TRUE),
    stringsAsFactors = FALSE
  )))
  write.csv(continuous, file.path(results_dir, "included_vs_unclassified_summary.csv"),
            row.names = FALSE, na = "")

  categorical_variables <- c("income_level", "sdi_level2021", "gbd_superregion")
  categorical_rows <- list()
  for (variable in categorical_variables) {
    counts <- as.data.frame(table(d$selection_group, d[[variable]], useNA = "ifany"),
                            stringsAsFactors = FALSE)
    names(counts) <- c("group", "category", "countries")
    counts$variable <- variable
    group_total <- aggregate(countries ~ group, data = counts, FUN = sum)
    counts <- merge(counts, group_total, by = "group", suffixes = c("", "_total"))
    counts$percent <- 100 * counts$countries / counts$countries_total
    categorical_rows[[variable]] <- counts
  }
  categorical <- do.call(rbind, categorical_rows)
  write.csv(categorical, file.path(results_dir, "included_vs_unclassified_categories.csv"),
            row.names = FALSE, na = "")

  cat("Saved 2021 included versus unclassified country comparisons.\n")
})
