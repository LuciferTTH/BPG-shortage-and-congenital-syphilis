# 04_additional_analyses.R
# WHO–GBD comparison, IQVIA procurement, HIV negative controls, and final tables.
# Run from the v7 analysis root after the primary and sensitivity files.
# Sections: 1) WHO–GBD comparison; 2) IQVIA analyses; 3) HIV source audit and
# optional negative-control models; 4) final tables and figure-asset audit.

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
gbd_who_files <- c(
  file.path(raw_data_dir, "WHOData", "Congenital syphilis rate per 100 000 live births, reported.csv"),
  file.path(raw_data_dir, "IHME-GBD_2021_DATA", "SexTransInfection-Incidence",
    "IHME-GBD_2021_DATA-SexTransInfection-Incidence.csv"),
  file.path(raw_data_dir, "IHME-GBD_2021_DATA", "Population",
    "IHME-GBD_2021_DATA-Population-Number-Birth_All_regions.csv"),
  file.path(processed_data_dir, "UN_WPP", "WHO_country_year_live_births_UN_WPP2024_2010_2021.csv")
)


# -----------------------------------------------------------------------------
# WHO surveillance and GBD estimate comparison
# -----------------------------------------------------------------------------
if (all(file.exists(gbd_who_files))) {


# -----------------------------------------------------------------------------
# Matched WHO–GBD rates and temporal trends
# -----------------------------------------------------------------------------
local({
  # Analysis plan: positive pairs, source-specific birth weights, descriptive rho.
  library(readr)
  library(ggplot2)
  library(patchwork)

  # 1. Read WHO data; select latest modified record and audit zero notifications.
  w <- readr::read_csv(file.path(raw_data_dir, "WHOData",
    "Congenital syphilis rate per 100 000 live births, reported.csv"),
    show_col_types = FALSE, name_repair = "minimal",
    locale = readr::locale(encoding = "UTF-8"))
  who <- data.frame(location_iso = as.character(w$SpatialDimValueCode),
    year = as.integer(w$Period), who_rate = as.numeric(w$FactValueNumeric),
    date_rank = as.numeric(as.Date(substr(as.character(w$DateModified), 1, 10))))
  who <- who[who$year %in% study_years, ]
  who$date_rank[!is.finite(who$date_rank)] <- -Inf
  who <- who[order(who$location_iso, who$year, -who$date_rank), ]
  who <- who[!duplicated(who[c("location_iso", "year")]), ]
  write.csv(who[is.finite(who$who_rate) & who$who_rate == 0, ],
    file.path(results_dir, "who_gbd_zero_notification_audit.csv"), row.names = FALSE)
  who <- who[is.finite(who$who_rate) & who$who_rate > 0, ]

  # 2. Calculate Birth-age GBD rates and join matching WPP live births.
  st <- read.csv(file.path(raw_data_dir, "IHME-GBD_2021_DATA",
    "SexTransInfection-Incidence", "IHME-GBD_2021_DATA-SexTransInfection-Incidence.csv"))
  st <- st[st$cause_name == "Syphilis" & st$age_name == "Birth" &
    st$sex_name == "Both" & st$metric_name == "Number",
    c("location_id", "year", "val")]
  names(st)[3] <- "gbd_cases"
  births <- read.csv(file.path(raw_data_dir, "IHME-GBD_2021_DATA", "Population",
    "IHME-GBD_2021_DATA-Population-Number-Birth_All_regions.csv"))
  births <- births[births$year %in% study_years & births$sex_id == 3 &
    births$metric_id == 1, c("location_id", "year", "val")]
  births <- births[!duplicated(births[c("location_id", "year")]), ]
  names(births)[3] <- "gbd_births"
  gbd <- merge(st, births, by = c("location_id", "year"))
  gbd$gbd_rate_per_100k <- 100000 * gbd$gbd_cases / gbd$gbd_births
  gbd <- merge(gbd, unique(cs_panel[c("location_id", "location_iso")]), by = "location_id")
  wpp <- read.csv(file.path(processed_data_dir, "UN_WPP",
    "WHO_country_year_live_births_UN_WPP2024_2010_2021.csv"))
  iso_col <- names(wpp)[grepl("iso|SpatialDimValueCode", names(wpp), ignore.case = TRUE)][1]
  year_col <- names(wpp)[grepl("year|period", names(wpp), ignore.case = TRUE)][1]
  birth_col <- names(wpp)[grepl("birth", names(wpp), ignore.case = TRUE) &
    !grepl("iso|year|name", names(wpp), ignore.case = TRUE)][1]
  wpp <- data.frame(location_iso = as.character(wpp[[iso_col]]),
    year = as.integer(wpp[[year_col]]), wpp_births = as.numeric(wpp[[birth_col]]))
  pairs <- merge(who, gbd, by = c("location_iso", "year"))
  pairs <- merge(pairs, wpp, by = c("location_iso", "year"))
  groups <- unique(cs_panel[c("location_iso", "shortage")])
  groups$group <- ifelse(is.na(groups$shortage), "unclassified",
    ifelse(groups$shortage == 1, "shortage reported", "no shortage reported"))
  pairs <- merge(pairs, groups[c("location_iso", "group")], by = "location_iso")
  pairs <- pairs[is.finite(pairs$gbd_rate_per_100k) & pairs$gbd_rate_per_100k > 0 &
    is.finite(pairs$wpp_births) & pairs$wpp_births > 0 &
    is.finite(pairs$gbd_births) & pairs$gbd_births > 0, ]
  stopifnot(nrow(pairs) == 466L, length(unique(pairs$location_iso)) == 90L,
    !anyDuplicated(pairs[c("location_iso", "year")]))
  write.csv(pairs, file.path(results_dir, "who_gbd_rate_pairs.csv"), row.names = FALSE)
  output_dir <- results_dir
  # 3. Calculate pooled source-specific weighted mean rates -------------------
  # WHO surveillance rates use UN WPP live births as weights.
  # GBD rates use GBD live births as weights.
  source_summary <- data.frame(
      source = c("WHO surveillance", "GBD estimates"),
      country_years = nrow(pairs),
      countries = length(unique(pairs$location_iso)),
      weighting_denominator = c("UN WPP live births", "GBD live births"),
      weighted_mean_rate_per_100000 = c(
          weighted.mean(pairs$who_rate, pairs$wpp_births),
          weighted.mean(pairs$gbd_rate_per_100k, pairs$gbd_births)
      )
  )
  write.csv(
      source_summary,
      file.path(output_dir, "who_gbd_pooled_weighted_means.csv"),
      row.names = FALSE
  )

  # 4. Calculate unweighted country-year Spearman correlations ----------------
  # These correlations compare paired rates directly. Birth weights are not
  # applied to the correlation, so this estimate does not change with weights.
  groups <- c(
      "Overall", "no shortage reported", "shortage reported", "unclassified"
  )
  country_year_rows <- list()

  for (group_name in groups) {
      if (group_name == "Overall") {
          group_data <- pairs
      } else {
          group_data <- pairs[pairs$group == group_name, ]
      }

      country_year_rows[[length(country_year_rows) + 1L]] <- data.frame(
          group = group_name,
          country_years = nrow(group_data),
          countries = length(unique(group_data$location_iso)),
          who_weighted_mean_per_100000 = weighted.mean(
              group_data$who_rate, group_data$wpp_births
          ),
          gbd_weighted_mean_per_100000 = weighted.mean(
              group_data$gbd_rate_per_100k, group_data$gbd_births
          ),
          spearman_rho_unweighted = cor(
              group_data$who_rate,
              group_data$gbd_rate_per_100k,
              method = "spearman"
          )
      )
  }

  country_year_summary <- do.call(rbind, country_year_rows)
  write.csv(
      country_year_summary,
      file.path(output_dir, "who_gbd_group_summary.csv"),
      row.names = FALSE
  )

  # 5. Calculate annual weighted mean rates -----------------------------------
  # Use the same matched positive-rate country-years in both annual means, while
  # applying each source's own live-birth denominator as its weight.
  annual_rows <- list()

  for (year_value in 2010:2021) {
      year_data <- pairs[pairs$year == year_value, ]

      annual_rows[[length(annual_rows) + 1L]] <- data.frame(
          year = year_value,
          countries = nrow(year_data),
          who_live_births_wpp = sum(year_data$wpp_births),
          gbd_live_births = sum(year_data$gbd_births),
          who_rate_weighted_per_100000 = weighted.mean(
              year_data$who_rate, year_data$wpp_births
          ),
          gbd_rate_weighted_per_100000 = weighted.mean(
              year_data$gbd_rate_per_100k, year_data$gbd_births
          )
      )
  }

  annual_rates <- do.call(rbind, annual_rows)
  annual_rho <- cor(
      annual_rates$who_rate_weighted_per_100000,
      annual_rates$gbd_rate_weighted_per_100000,
      method = "spearman"
  )
  write.csv(
      annual_rates,
      file.path(output_dir, "who_gbd_annual_weighted_rates.csv"),
      row.names = FALSE
  )

  # 6. Save correlation summary -----------------------------------------------
  correlation_summary <- data.frame(
      comparison = c(
          "Matched country-year rates",
          "Annual source-specific weighted mean rates"
      ),
      number_of_pairs = c(nrow(pairs), nrow(annual_rates)),
      spearman_rho_unweighted = c(
          country_year_summary$spearman_rho_unweighted[1], annual_rho
      ),
      weighting_note = c(
          "Correlation is unweighted; 466 matched positive-rate country-years",
          paste(
              "Annual WHO means weighted by WPP births; annual GBD means",
              "weighted by GBD births; correlation across 12 years is unweighted"
          )
      )
  )
  write.csv(
      correlation_summary,
      file.path(output_dir, "who_gbd_agreement_statistics.csv"),
      row.names = FALSE
  )

  # 7. Draw the paired rates and annual trends --------------------------------
  annual_long <- rbind(
      data.frame(
          year = annual_rates$year,
          source = "WHO surveillance",
          rate = annual_rates$who_rate_weighted_per_100000
      ),
      data.frame(
          year = annual_rates$year,
          source = "GBD estimates",
          rate = annual_rates$gbd_rate_weighted_per_100000
      )
  )

  plot_theme <- publication_theme(base_size = 11) +
      theme(
          axis.text = element_text(colour = "black"),
          plot.title = element_text(face = "bold", size = 12),
          legend.position = "bottom",
          legend.title = element_blank()
      )

  rate_scatter <- ggplot(
      pairs,
      aes(x = who_rate, y = gbd_rate_per_100k)
  ) +
      geom_abline(
          slope = 1, intercept = 0, linetype = "dotted", colour = "grey55"
      ) +
      geom_point(colour = "#415A77", alpha = 0.48, size = 1.5) +
      scale_x_log10(labels = scales::label_number()) +
      scale_y_log10(labels = scales::label_number()) +
      coord_fixed() +
      labs(
          title = "A  Matched country-year rates",
          subtitle = sprintf(
              "466 pairs; unweighted rho = %.3f",
              country_year_summary$spearman_rho_unweighted[1]
          ),
          x = "WHO rate per 100,000 live births",
          y = "GBD rate per 100,000 live births"
      ) +
      plot_theme

  annual_trend <- ggplot(
      annual_long,
      aes(x = year, y = rate, colour = source, shape = source, group = source)
  ) +
      geom_line(linewidth = 0.85) +
      geom_point(size = 2.1) +
      scale_colour_manual(
          values = c("WHO surveillance" = "#176B9A", "GBD estimates" = "#C63D40")
      ) +
      scale_shape_manual(
          values = c("WHO surveillance" = 16, "GBD estimates" = 17)
      ) +
      scale_x_continuous(breaks = c(seq(2010, 2020, 2), 2021), limits = c(2010, 2021)) +
      scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.06))) +
      labs(
          title = "B  Annual source-specific weighted rates",
          subtitle = sprintf(
              "12 annual pairs; unweighted rho = %.3f", annual_rho
          ),
          x = "Year",
          y = "Cases per 100,000 live births",
          colour = NULL,
          shape = NULL
      ) +
      plot_theme

  comparison_figure <- (rate_scatter + annual_trend + plot_layout(widths = c(1, 1.2))) +
      plot_annotation(
          caption = paste(
              "Annual WHO rates weighted by UN WPP live births; annual GBD rates",
              "weighted by GBD live births. WHO zero notifications excluded."
          )
      )

  ggsave(
      file.path(figures_dir, "Figure_S14_WHO_GBD_comparison.png"),
      comparison_figure,
      width = 10,
      height = 4.8,
      dpi = 300,
      bg = "white"
  )
  ggsave(
      file.path(figures_dir, "Figure_S14_WHO_GBD_comparison.pdf"),
      comparison_figure,
      width = 10,
      height = 4.8,
      device = grDevices::cairo_pdf,
      bg = "white"
  )

  # 8. Print the updated results to the R console -----------------------------
  print(source_summary)
  print(country_year_summary)
  print(annual_rates)
  print(correlation_summary)
})

} else {
  message('WHO–GBD comparison skipped because one or more WHO, GBD, or UN WPP extracts are missing.')
}
iqvia_raw_dir <- file.path(raw_data_dir, "IQVIA")
iqvia_processed_dir <- file.path(processed_data_dir, "IQVIA")

iqvia_population_dir <- file.path(raw_data_dir, "IHME-GBD_2021_DATA", "Population",
  "Population-All_ages-All_regions")
iqvia_gbd_incidence_dir <- file.path(raw_data_dir, "IHME-GBD_2021_DATA",
  "Syphilis-Death_Incidence", "All_ages-All_regions")
iqvia_gbd_incidence_files <- list.files(iqvia_gbd_incidence_dir,
  pattern = "\\.csv$", full.names = TRUE)
iqvia_required_files <- c(
  file.path(iqvia_raw_dir, "data", "penicillin_G_2012-2022_global.csv"),
  file.path(iqvia_processed_dir, "BPG drug products to be included.xlsx"),
  file.path(iqvia_processed_dir, "CountryNameIQVIAGBDMatched.csv"),
  file.path(iqvia_processed_dir, "data_volume_full.RData"))
has_iqvia_inputs <- length(iqvia_required_files) == 4L &&
  all(file.exists(iqvia_required_files)) &&
  length(list.files(iqvia_population_dir, pattern = "\\.zip$", full.names = TRUE)) > 0L &&
  length(iqvia_gbd_incidence_files) > 0L


# -----------------------------------------------------------------------------
# IQVIA procurement analyses
# -----------------------------------------------------------------------------
if (has_iqvia_inputs) {


# -----------------------------------------------------------------------------
# Rebuild procurement panel and fit dose-volume LMMs
# -----------------------------------------------------------------------------
local({
  # Rebuild IQVIA dose-equivalents from quarter-level data and run weighted LMMs.
  # Run from the analysis root. Raw IQVIA source inputs are not included in this package.

  suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(lme4)
    library(lmerTest)
  })

  # Raw quarterly values and curated matching/conversion files have separate roots.
  results_dir <- file.path("results", "iqvia")
  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
  write_csv <- function(x, filename) {
    utils::write.csv(x, file.path(results_dir, filename), row.names = FALSE, na = "")
  }
  model_control <- lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 200000))

  build_cs_panel <- function() {
    as.data.frame(cs_panel)
  }

  load_total_population <- function() {
    population_dir <- file.path(raw_data_dir, "IHME-GBD_2021_DATA", "Population",
                                "Population-All_ages-All_regions")
    zip_files <- list.files(population_dir, pattern = "\\.zip$", full.names = TRUE)
    population_parts <- list()
    for (zip_path in zip_files) {
      csv_name <- utils::unzip(zip_path, list = TRUE)$Name
      csv_name <- csv_name[grepl("\\.csv$", csv_name)]
      population_data <- utils::read.csv(unz(zip_path, csv_name), stringsAsFactors = FALSE)
      population_data <- population_data[population_data$age_id == 22 &
        population_data$sex_id == 3 & population_data$metric_id == 1 &
        population_data$year >= 2010 & population_data$year <= 2021,
        c("location_id", "year", "val")]
      population_parts[[length(population_parts) + 1L]] <- population_data
    }
    population <- do.call(rbind, population_parts)
    names(population)[names(population) == "val"] <- "total_population"
    if (anyDuplicated(population[c("location_id", "year")])) {
      stop("Duplicate country-year rows in the GBD population extract.")
    }
    population
  }
  suppressPackageStartupMessages({
      library(dplyr)
      library(tidyr)
      library(ggplot2)
  })

  to_numeric <- function(x) suppressWarnings(as.numeric(gsub(",", "", trimws(as.character(x)),
                                                              fixed = TRUE)))
  raw_file <- file.path(iqvia_raw_dir, "data", "penicillin_G_2012-2022_global.csv")
  product_file <- file.path(iqvia_processed_dir, "BPG drug products to be included.xlsx")
  country_file <- file.path(iqvia_processed_dir, "CountryNameIQVIAGBDMatched.csv")
  if (!all(file.exists(c(raw_file, product_file, country_file)))) {
      stop("Required IQVIA raw source, product list, or country lookup is missing.")
  }

  # ---- 1. Product and source reconciliation --------------------------------
  raw <- utils::read.csv(raw_file, check.names = FALSE, stringsAsFactors = FALSE,
                         fileEncoding = "UTF-8-BOM")
  products <- openxlsx::read.xlsx(product_file)
  country_map <- utils::read.csv(country_file, stringsAsFactors = FALSE,
                                 fileEncoding = "UTF-8-BOM")
  product_keys <- c("international_product", "molecule_list", "nfc123",
                    "international_strength", "international_volume")
  if (!all(c(product_keys, "transfer") %in% names(products))) {
      stop("The product workbook lacks its product keys or dose-equivalence factor.")
  }
  if (!all(c("country", product_keys) %in% names(raw))) {
      stop("The IQVIA raw file lacks country or product identifiers.")
  }
  if (!all(c("country_iqvia", "iso3", "location_id", "location_name") %in% names(country_map))) {
      stop("The IQVIA country lookup lacks required identifiers.")
  }
  if (anyDuplicated(products[product_keys])) stop("Duplicate product keys in the product workbook.")
  products$transfer <- to_numeric(products$transfer)
  write_csv(products[, c(product_keys, "transfer")],
            "iqvia_product_dose_equivalence_map.csv")
  raw$product_key <- do.call(paste, c(raw[product_keys], sep = "\r"))
  products$product_key <- do.call(paste, c(products[product_keys], sep = "\r"))
  raw$country_key <- toupper(trimws(as.character(raw$country)))
  products$transfer <- as.numeric(products$transfer)
  selected <- dplyr::inner_join(raw, products[, c("product_key", "transfer")],
                                by = "product_key")
  if (!nrow(selected)) stop("No IQVIA rows matched the curated BPG product list.")

  quarter_columns <- grep("^q[1-4]_[0-9]{4}_standard_units$", names(selected),
                          value = TRUE, ignore.case = TRUE)
  quarter_map <- data.frame(
      column = quarter_columns,
      quarter = as.integer(sub("^q([1-4])_.*$", "\\1", tolower(quarter_columns))),
      year = as.integer(sub("^q[1-4]_([0-9]{4})_.*$", "\\1", tolower(quarter_columns))),
      stringsAsFactors = FALSE)
  selected <- selected[, unique(c("country", "country_key", "product_key", "transfer",
                                  quarter_columns)), drop = FALSE]
  long <- tidyr::pivot_longer(selected, cols = dplyr::all_of(quarter_columns),
                              names_to = "quarter_column", values_to = "reported_units")
  long <- merge(long, quarter_map, by.x = "quarter_column", by.y = "column", all.x = TRUE)
  long$reported_units <- to_numeric(long$reported_units)
  long$transfer <- to_numeric(long$transfer)
  long$volume_dose_equivalents <- long$reported_units * long$transfer
  long$year <- as.integer(long$year)
  long$quarter <- as.integer(long$quarter)
  long <- long[long$year >= 2012 & long$year <= 2021, ]

  country_map$country_key <- toupper(trimws(as.character(country_map$country_iqvia)))
  if (anyDuplicated(country_map$country_key)) stop("Duplicate IQVIA country names in lookup.")
  long <- dplyr::left_join(long, country_map[, c("country_key", "iso3", "location_id",
                                                  "location_name")], by = "country_key")
  unmatched <- unique(long[is.na(long$location_id), c("country_key", "country")])
  long <- long[!is.na(long$location_id) & !is.na(long$iso3) & nzchar(long$iso3), ]

  # Each country-quarter is complete only when every selected product has a
  # numeric source value. Partially reported quarters are not treated as zero.
  quarter_rows <- split(long, interaction(long$location_id, long$year,
                                          long$quarter, drop = TRUE))
  quarter_observed <- do.call(rbind, lapply(quarter_rows, function(d) {
      vals <- d$volume_dose_equivalents
      all_values <- all(is.finite(vals))
      data.frame(
          location_id = as.character(d$location_id[1]),
          location_iso = as.character(d$iso3[1]),
          location_name = as.character(d$location_name[1]),
          year = d$year[1], quarter = d$quarter[1],
          selected_product_rows = nrow(d),
          observed_product_rows = sum(is.finite(vals)),
          observation_fraction = mean(is.finite(vals)),
          quarter_volume = if (any(is.finite(vals))) sum(vals[is.finite(vals)]) else NA_real_,
          all_product_cells_observed = all_values,
          explicit_all_product_zero = all_values && all(vals == 0),
          stringsAsFactors = FALSE)
  }))
  if (anyDuplicated(quarter_observed[c("location_id", "year", "quarter")])) {
      stop("Duplicate country-quarter rows after raw source aggregation.")
  }

  country_keys <- unique(country_map[, c("location_id", "iso3", "location_name")])
  names(country_keys)[names(country_keys) == "iso3"] <- "location_iso"
  expected <- merge(country_keys,
                     expand.grid(year = 2012:2021, quarter = 1:4),
                     by = NULL, sort = FALSE)
  quarter_audit <- merge(expected, quarter_observed,
                         by = c("location_id", "location_iso", "location_name",
                                "year", "quarter"), all.x = TRUE, sort = FALSE)
  quarter_audit$all_product_cells_observed[is.na(quarter_audit$all_product_cells_observed)] <- FALSE
  quarter_audit$explicit_all_product_zero[is.na(quarter_audit$explicit_all_product_zero)] <- FALSE
  quarter_audit$observation_fraction[is.na(quarter_audit$observation_fraction)] <- 0
  quarter_audit$quarter_status <- ifelse(
      quarter_audit$all_product_cells_observed,
      ifelse(quarter_audit$explicit_all_product_zero, "Complete observed zero",
             "Complete observed product values"),
      "Missing or partial product values")

  annual_parts <- split(quarter_audit,
                        interaction(quarter_audit$location_id, quarter_audit$year,
                                    drop = TRUE))
  annual <- do.call(rbind, lapply(annual_parts, function(d) {
      d <- d[order(d$quarter), ]
      all_q_observed <- nrow(d) == 4L && all(d$all_product_cells_observed)
      any_q_observed <- nrow(d) == 4L && all(is.finite(d$quarter_volume))
      annual_value <- if (all_q_observed) sum(d$quarter_volume) else NA_real_
      data.frame(
          location_id = as.character(d$location_id[1]),
          location_iso = as.character(d$location_iso[1]),
          location_name = as.character(d$location_name[1]),
          year = d$year[1],
          n_quarters_observed = sum(is.finite(d$quarter_volume)),
          n_quarters_complete_products = sum(d$all_product_cells_observed),
          annual_volume_dose_equivalents = annual_value,
          annual_volume_any_quarters = if (any_q_observed) sum(d$quarter_volume) else NA_real_,
          verified_annual_zero = all_q_observed && is.finite(annual_value) && annual_value == 0,
          annual_status = if (d$year[1] == 2012) "Q4-only source year; excluded from annual models"
                           else if (all_q_observed) "All four quarters and selected products observed"
                           else "Incomplete selected-product source values; retained missing",
          stringsAsFactors = FALSE)
  }))

  # ---- 2. Attach CS, syphilis case counts, population, and exposure ---------
  main <- build_cs_panel()
  main$location_id <- as.character(main$location_id)
  main$location_iso <- toupper(trimws(as.character(main$location_iso)))
  pop <- load_total_population()
  pop$location_id <- as.character(pop$location_id)
  main <- merge(main, pop[, c("location_id", "year", "total_population")],
                by = c("location_id", "year"), all.x = TRUE)

  # The archived IQVIA panel contains all-age and female GBD case counts, but
  # not the guideline-recommended denominator. Reconstruct that denominator
  # from GBD incidence counts for both-sex children aged 0-14 years plus
  # females aged 15-49 years.
  incidence_env <- new.env(parent = emptyenv())
  load(file.path(iqvia_processed_dir, "data_volume_full.RData"), envir = incidence_env)
  incidence_archive <- as.data.frame(incidence_env$data_volume_full)
  incidence_archive$location_id <- as.character(incidence_archive$location_id)
  incidence_archive$year <- as.integer(incidence_archive$year)
  archived_fields <- c("location_id", "year", "syphilis_incidence_all",
                       "syphilis_incidence_female",
                       "syphilis_incidence_female15to49")
  if (!all(archived_fields %in% names(incidence_archive))) {
      stop("The archived IQVIA panel lacks expected GBD incidence fields.")
  }
  incidence_archive <- incidence_archive[, archived_fields]
  if (anyDuplicated(incidence_archive[c("location_id", "year")])) {
      stop("Duplicate GBD syphilis-incidence country-years in the archived IQVIA panel.")
  }

  gbd_incidence_dir <- file.path(raw_data_dir, "IHME-GBD_2021_DATA",
      "Syphilis-Death_Incidence", "All_ages-All_regions")
  gbd_incidence_files <- list.files(gbd_incidence_dir, pattern = "\\.csv$",
                                    full.names = TRUE)
  if (length(gbd_incidence_files) == 0L) {
      stop("Cannot find the GBD age-specific syphilis incidence CSV extracts.")
  }
  gbd_incidence_parts <- lapply(gbd_incidence_files, function(gbd_file) {
      gbd <- readr::read_csv(gbd_file,
          col_select = c("measure_name", "metric_name", "cause_name", "location_id",
                         "year", "age_name", "sex_name", "val"),
          show_col_types = FALSE, progress = FALSE)
      keep <- gbd$cause_name == "Syphilis" & gbd$measure_name == "Incidence" &
          gbd$metric_name == "Number" & gbd$year %in% study_years &
          ((gbd$age_name == "0-14 years" & gbd$sex_name == "Both") |
           (gbd$age_name == "15-49 years" & gbd$sex_name == "Female") |
           (gbd$age_name == "All ages" & gbd$sex_name %in% c("Both", "Female")))
      as.data.frame(gbd[keep, c("location_id", "year", "age_name", "sex_name", "val")])
  })
  gbd_incidence_parts <- gbd_incidence_parts[
      vapply(gbd_incidence_parts, nrow, integer(1)) > 0L]
  if (length(gbd_incidence_parts) == 0L) {
      stop("No required GBD syphilis incidence strata were found in the CSV extracts.")
  }
  gbd_incidence_strata <- do.call(rbind, gbd_incidence_parts)
  gbd_incidence_strata$location_id <- as.character(gbd_incidence_strata$location_id)
  gbd_incidence_strata$year <- as.integer(gbd_incidence_strata$year)
  stratum_keys <- c("location_id", "year", "age_name", "sex_name")
  if (anyDuplicated(gbd_incidence_strata[stratum_keys])) {
      stop("Duplicate location-year-age-sex rows in the GBD incidence extracts.")
  }

  gbd_children <- gbd_incidence_strata[
      gbd_incidence_strata$age_name == "0-14 years" &
          gbd_incidence_strata$sex_name == "Both",
      c("location_id", "year", "val")]
  names(gbd_children)[3] <- "syphilis_incidence_children_0to14"
  gbd_women <- gbd_incidence_strata[
      gbd_incidence_strata$age_name == "15-49 years" &
          gbd_incidence_strata$sex_name == "Female",
      c("location_id", "year", "val")]
  names(gbd_women)[3] <- "syphilis_incidence_female15to49_gbd"
  gbd_all_both <- gbd_incidence_strata[
      gbd_incidence_strata$age_name == "All ages" &
          gbd_incidence_strata$sex_name == "Both",
      c("location_id", "year", "val")]
  names(gbd_all_both)[3] <- "syphilis_incidence_all_gbd"
  gbd_all_female <- gbd_incidence_strata[
      gbd_incidence_strata$age_name == "All ages" &
          gbd_incidence_strata$sex_name == "Female",
      c("location_id", "year", "val")]
  names(gbd_all_female)[3] <- "syphilis_incidence_female_gbd"

  gbd_incidence_by_country_year <- merge(gbd_children, gbd_women,
      by = c("location_id", "year"), all = TRUE, sort = FALSE)
  gbd_incidence_by_country_year <- merge(gbd_incidence_by_country_year,
      gbd_all_both, by = c("location_id", "year"), all = TRUE, sort = FALSE)
  gbd_incidence_by_country_year <- merge(gbd_incidence_by_country_year,
      gbd_all_female, by = c("location_id", "year"), all = TRUE, sort = FALSE)
  required_gbd_fields <- c("syphilis_incidence_children_0to14",
      "syphilis_incidence_female15to49_gbd", "syphilis_incidence_all_gbd",
      "syphilis_incidence_female_gbd")
  if (any(!stats::complete.cases(gbd_incidence_by_country_year[, required_gbd_fields]))) {
      stop("One or more required GBD age-sex strata are missing for a location-year.")
  }
  gbd_incidence_by_country_year$syphilis_incidence_guidrecom <-
      gbd_incidence_by_country_year$syphilis_incidence_children_0to14 +
      gbd_incidence_by_country_year$syphilis_incidence_female15to49_gbd

  incidence_check <- merge(incidence_archive, gbd_incidence_by_country_year,
      by = c("location_id", "year"), all.x = TRUE, sort = FALSE)
  if (nrow(incidence_check) != nrow(incidence_archive) ||
      any(!stats::complete.cases(incidence_check[, c(required_gbd_fields,
          "syphilis_incidence_guidrecom")]))) {
      stop("Some archived IQVIA country-years could not be matched to GBD incidence strata.")
  }
  incidence_validation <- data.frame(
      comparison = c("All-age, both sexes", "All-age, female",
                     "Female aged 15-49 years"),
      country_years = nrow(incidence_check),
      exact_matches = c(
          sum(incidence_check$syphilis_incidence_all ==
                  incidence_check$syphilis_incidence_all_gbd),
          sum(incidence_check$syphilis_incidence_female ==
                  incidence_check$syphilis_incidence_female_gbd),
          sum(incidence_check$syphilis_incidence_female15to49 ==
                  incidence_check$syphilis_incidence_female15to49_gbd)),
      max_absolute_difference = c(
          max(abs(incidence_check$syphilis_incidence_all -
                      incidence_check$syphilis_incidence_all_gbd)),
          max(abs(incidence_check$syphilis_incidence_female -
                      incidence_check$syphilis_incidence_female_gbd)),
          max(abs(incidence_check$syphilis_incidence_female15to49 -
                      incidence_check$syphilis_incidence_female15to49_gbd))),
      stringsAsFactors = FALSE)
  write_csv(incidence_validation, "iqvia_gbd_incidence_strata_validation.csv")
  if (any(incidence_validation$max_absolute_difference > 1e-8)) {
      stop("Archived incidence counts do not agree with the corresponding GBD extracts.")
  }
  incidence <- incidence_check[, c("location_id", "year", "syphilis_incidence_all",
      "syphilis_incidence_female", "syphilis_incidence_guidrecom")]
  if (anyDuplicated(incidence[c("location_id", "year")])) {
      stop("Duplicate GBD syphilis-incidence country-years after stratum reconstruction.")
  }
  main <- merge(main, incidence, by = c("location_id", "year"), all.x = TRUE)
  volume_panel <- merge(main, annual, by = c("location_id", "location_iso", "year"),
                        all = FALSE, sort = FALSE)
  volume_panel$syphilis_cases <- as.numeric(volume_panel$syphilis_incidence_all)
  volume_panel$syphilis_cases_female <- as.numeric(volume_panel$syphilis_incidence_female)
  # Guideline-recommended treatment population: incident syphilis among women aged
  # 15-49 years plus all children aged 0-14 years, reconstructed from GBD strata above. These are the
  # groups for whom BPG is the guideline-recommended treatment, so this count is
  # the procurement-ratio denominator. Country-year model weights use guideline-recommended incident syphilis cases.
  volume_panel$syphilis_cases_guidrecom <-
      as.numeric(volume_panel$syphilis_incidence_guidrecom)
  volume_panel$cs_rate <- as.numeric(volume_panel$val) /
      as.numeric(volume_panel$population) * 1000
  volume_panel$syphilis_rate_per_100k <- volume_panel$syphilis_cases /
      as.numeric(volume_panel$total_population) * 1e5
  volume_panel$ratio_per_1000_all_cases <- ifelse(
      is.finite(volume_panel$annual_volume_dose_equivalents) &
          is.finite(volume_panel$syphilis_cases) & volume_panel$syphilis_cases > 0,
      1000 * volume_panel$annual_volume_dose_equivalents / volume_panel$syphilis_cases,
      NA_real_)
  volume_panel$ratio_per_1000_female_cases <- ifelse(
      is.finite(volume_panel$annual_volume_dose_equivalents) &
          is.finite(volume_panel$syphilis_cases_female) & volume_panel$syphilis_cases_female > 0,
      1000 * volume_panel$annual_volume_dose_equivalents / volume_panel$syphilis_cases_female,
      NA_real_)
  # Primary exposure: BPG 2.4-MIU dose-equivalents per 1,000 guideline-recommended
  # patients. ratio_per_1000_guideline_need is retained as the historical name for
  # the same quantity so existing readers keep working.
  volume_panel$ratio_per_1000_guidrecom <- ifelse(
      is.finite(volume_panel$annual_volume_dose_equivalents) &
          is.finite(volume_panel$syphilis_cases_guidrecom) &
          volume_panel$syphilis_cases_guidrecom > 0,
      1000 * volume_panel$annual_volume_dose_equivalents /
          volume_panel$syphilis_cases_guidrecom,
      NA_real_)
  volume_panel$ratio_per_1000_guideline_need <- volume_panel$ratio_per_1000_guidrecom
  volume_panel$group <- ifelse(is.na(volume_panel$shortage), "unclassified",
      ifelse(volume_panel$shortage == 1, "shortage reported", "no shortage reported"))
  volume_panel$annual_valid_for_primary <- volume_panel$year >= 2013 &
      is.finite(volume_panel$annual_volume_dose_equivalents) &
      is.finite(volume_panel$syphilis_cases_guidrecom) &
      volume_panel$syphilis_cases_guidrecom > 0 &
      is.finite(volume_panel$cs_rate) & is.finite(volume_panel$population) &
      volume_panel$population > 0 & is.finite(volume_panel$syphilis_rate_per_100k)

  write_csv(quarter_audit, "iqvia_country_quarter_source_audit.csv")
  write_csv(annual, "iqvia_country_year_source_audit.csv")
  write_csv(unmatched, "iqvia_unmatched_country_names.csv")
  write_csv(volume_panel[, c("location_id", "location_iso", "location_name.x", "year",
                             "group", "annual_volume_dose_equivalents",
                             "n_quarters_observed", "n_quarters_complete_products",
                             "annual_status", "verified_annual_zero", "syphilis_cases",
                             "syphilis_cases_female", "syphilis_cases_guidrecom",
                             "syphilis_rate_per_100k", "cs_rate",
                             "ratio_per_1000_all_cases", "ratio_per_1000_female_cases",
                             "ratio_per_1000_guidrecom", "ratio_per_1000_guideline_need",
                             "population")],
            "iqvia_incidence_case_ratio_panel.csv")

  annual_valid <- volume_panel[volume_panel$year >= 2013, ]
  audit <- do.call(rbind, lapply(split(annual_valid, annual_valid$year), function(d) {
      data.frame(year = d$year[1], countries_with_matched_panel = nrow(d),
                 countries_with_complete_iqvia_year = sum(is.finite(d$annual_volume_dose_equivalents)),
                 complete_iqvia_and_guideline_cases = sum(is.finite(d$annual_volume_dose_equivalents) &
                                                                 is.finite(d$syphilis_cases_guidrecom)),
                 verified_zero_volume_country_years = sum(d$verified_annual_zero, na.rm = TRUE),
                 eligible_for_main_ratio = sum(d$annual_valid_for_primary),
                 stringsAsFactors = FALSE)
  }))
  write_csv(audit, "iqvia_incident_ratio_coverage_by_year.csv")

  # ---- 3. Weighted descriptive summaries -----------------------------------
  descriptive <- do.call(rbind, lapply(split(annual_valid[annual_valid$annual_valid_for_primary, ],
                                            list(annual_valid$year[annual_valid$annual_valid_for_primary],
                                                 annual_valid$group[annual_valid$annual_valid_for_primary]),
                                            drop = TRUE), function(d) {
      d <- d[is.finite(d$ratio_per_1000_guidrecom) &
                 is.finite(d$syphilis_cases_guidrecom) & d$syphilis_cases_guidrecom > 0, ]
      if (!nrow(d)) return(NULL)
      data.frame(
          year = d$year[1], group = d$group[1], n_countries = length(unique(d$location_id)),
          total_dose_equivalents = sum(d$annual_volume_dose_equivalents),
          total_guideline_recommended_patients = sum(d$syphilis_cases_guidrecom),
          total_syphilis_cases_all_ages = sum(d$syphilis_cases),
          aggregate_ratio_per_1000_guidrecom =
              1000 * sum(d$annual_volume_dose_equivalents) /
                  sum(d$syphilis_cases_guidrecom),
          aggregate_ratio_per_1000_all_cases =
              1000 * sum(d$annual_volume_dose_equivalents) / sum(d$syphilis_cases),
          case_weighted_mean_ratio = stats::weighted.mean(d$ratio_per_1000_guidrecom,
                                                          d$syphilis_cases_guidrecom),
          stringsAsFactors = FALSE)
  }))
  write_csv(descriptive, "iqvia_incident_ratio_group_year_description.csv")
  ratio_plot <- ggplot(descriptive, aes(x = year, y = aggregate_ratio_per_1000_guidrecom,
                                        colour = group)) +
      geom_line(linewidth = .9) + geom_point(size = 1.4) +
      scale_x_continuous(breaks = 2013:2021) +
      labs(title = "BPG 2.4-MIU dose-equivalents per 1,000 guideline-recommended patients",
           subtitle = "Guideline-recommended patients: women aged 15-49 years and children aged 0-14 years",
           x = "Year", y = "Dose-equivalents per 1,000\nguideline-recommended patients", colour = NULL) +
      theme_bw() + theme(legend.position = "top")
  ggsave(file.path(results_dir, "iqvia_incident_ratio_by_group.png"),
         ratio_plot, width = 9, height = 5.5, dpi = 300)

  # ---- 4. Weighted country-panel LMMs (no ARMA) ----------------------------
  model_data <- volume_panel[volume_panel$annual_valid_for_primary &
                                 volume_panel$year >= 2013 & volume_panel$year <= 2021, ]
  model_data$location_id <- factor(model_data$location_id)
  model_data$gbd_superregion <- factor(model_data$gbd_superregion)
  model_data$year_c <- model_data$year - 2016
  model_data$time_c <- model_data$year - 2016
  model_data$post <- as.integer(model_data$year >= 2016)

  covars <- c("hdi", "sdi", "gni", "syphilis_prev_female", "syphilis_prev_male",
              "abor_misc_rate")
  z_score <- function(x) {
      x <- as.numeric(x)
      s <- stats::sd(x, na.rm = TRUE)
      if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
      (x - mean(x, na.rm = TRUE)) / s
  }
  for (v in covars) model_data[[paste0(v, "_z")]] <- z_score(model_data[[v]])
  model_data$syph_rate_z <- z_score(model_data$syphilis_rate_per_100k)

  # ratio_log is ln(1 + ratio), the same transformation used by the group-level
  # LMM and display in 46_iqvia_log_ratio_by_group.R. The offset keeps the
  # transform defined if a verified zero-volume country-year enters the panel.
  add_within_between <- function(d, source, stem) {
      d[[paste0(stem, "_log")]] <- log1p(pmax(0, as.numeric(d[[source]])))
      means <- ave(d[[paste0(stem, "_log")]], d$location_id,
                   FUN = function(x) mean(x, na.rm = TRUE))
      d[[paste0(stem, "_between")]] <- means
      d[[paste0(stem, "_within")]] <- d[[paste0(stem, "_log")]] - means
      d
  }
  # The primary stem "ratio" is built on the guideline-recommended denominator, so
  # ratio_log, ratio_within and ratio_between all refer to BPG dose-equivalents per
  # 1,000 guideline-recommended patients. The all-case and female-case ratios are
  # retained as denominator sensitivities, and the dose and case-count components
  # use the same guideline-recommended population.
  model_data <- add_within_between(model_data, "ratio_per_1000_guidrecom", "ratio")
  model_data <- add_within_between(model_data, "annual_volume_dose_equivalents", "dose")
  model_data <- add_within_between(model_data, "syphilis_cases_guidrecom", "cases")
  model_data <- add_within_between(model_data, "ratio_per_1000_all_cases", "allcase_ratio")
  model_data <- add_within_between(model_data, "ratio_per_1000_female_cases", "female_ratio")

  fit_lmm <- function(d, label, exposure_terms,
                      include_case_rate = TRUE,
                      use_shortage_cits = FALSE) {
      d <- d[order(d$location_id, d$year), , drop = FALSE]
      if (use_shortage_cits) {
          # The binary shortage CITS has a distinct exposure and retains its
          # prespecified broader adjustment set.
          common <- c("hdi_z", "sdi_z", "gni_z", "syphilis_prev_female_z",
                      "syphilis_prev_male_z", "abor_misc_rate_z", "gbd_superregion")
          if (include_case_rate) common <- c(common, "syph_rate_z")
      } else {
          # Procurement intensity is defined per incident syphilis case; do not
          # additionally adjust for syphilis prevalence or incidence burden.
          common <- c("sdi_z", "gbd_superregion")
      }
      if (use_shortage_cits) {
          fixed <- stats::as.formula(paste0(
              "cs_rate ~ time_c * shortage * post + ",
              paste(common, collapse = " + ")))
      } else {
          fixed <- stats::as.formula(paste0(
              "cs_rate ~ factor(year) + ", paste(exposure_terms, collapse = " + "),
              " + ", paste(common, collapse = " + ")))
      }
      d <- d[stats::complete.cases(d[, unique(c(all.vars(fixed), "location_id",
                                                 "gbd_superregion", "year_c")), drop = FALSE]) &
                 is.finite(d$cs_rate) & is.finite(d$population) & d$population > 0, ]
      if (nrow(d) < 30L || length(unique(d$location_id)) < 8L) {
          return(list(result = data.frame(
              model = label, term = NA_character_, estimate = NA_real_, se = NA_real_,
              ci_lower = NA_real_, ci_upper = NA_real_, p = NA_real_,
              model_structure = "guideline-recommended incident-case-weighted random-intercept/random-slope LMM; no ARMA",
              singular_fit = NA,
              n_country_years = nrow(d),
              n_countries = length(unique(d$location_id)),
              status = "insufficient complete rows/countries"),
              grid = data.frame(model = label, status = "not fitted"), fit = NULL,
              data = d, diagnostics = NULL))
      }
      d$location_id <- droplevels(factor(d$location_id))
      d$gbd_superregion <- droplevels(factor(d$gbd_superregion))
      # Weight country-years by estimated guideline-recommended GBD incident syphilis case
      # counts for the country-level average association. The
      # procurement ratio itself remains defined per 1,000 guideline-recommended
      # incident cases (women aged 15-49 years plus children aged 0-14 years).
      # Rescale weights to a mean of 1 for numerical stability in lme4.
      d$prior_weight <- d$syphilis_cases_guidrecom / mean(d$syphilis_cases_guidrecom)
      full_formula <- stats::update(fixed, . ~ . + (1 + year_c | location_id))
      fit_warnings <- character()
      theta_starts <- list(
          default = NULL,
          start_a = c(150.048207313949, 28.5688389337839, 2.53797203789327),
          start_b = c(62.3299316056475, 56.5044512309326, 15.5967255963415)
      )
      candidate_fits <- list()
      for (start_name in names(theta_starts)) {
          theta <- theta_starts[[start_name]]
          fit_try <- tryCatch(withCallingHandlers(
              if (is.null(theta)) {
                  lmerTest::lmer(full_formula, data = d, weights = prior_weight,
                      REML = FALSE, control = model_control, na.action = na.fail)
              } else {
                  lmerTest::lmer(full_formula, data = d, weights = prior_weight,
                      REML = FALSE, control = model_control, na.action = na.fail,
                      start = list(theta = theta))
              },
              warning = function(w) {
                  fit_warnings <<- c(fit_warnings,
                                     paste0(start_name, ": ", conditionMessage(w)))
                  invokeRestart("muffleWarning")
              }), error = function(e) e)
          if (inherits(fit_try, "merMod")) {
              candidate_fits[[start_name]] <- fit_try
          } else {
              fit_warnings <- c(fit_warnings,
                                paste(start_name, conditionMessage(fit_try)))
          }
      }
      selected_start <- NA_character_
      loglik_gap <- NA_real_
      if (length(candidate_fits)) {
          candidate_loglik <- vapply(candidate_fits,
              function(x) as.numeric(logLik(x)), numeric(1))
          best_i <- which.max(candidate_loglik)
          best <- candidate_fits[[best_i]]
          selected_start <- names(candidate_fits)[best_i]
          ordered_loglik <- sort(candidate_loglik, decreasing = TRUE)
          if (length(ordered_loglik) > 1L)
              loglik_gap <- ordered_loglik[1] - ordered_loglik[2]
      } else {
          best <- simpleError(paste(unique(fit_warnings), collapse = " | "))
      }
      if (!inherits(best, "merMod")) {
          return(list(result = data.frame(model = label, term = NA_character_,
                                          estimate = NA_real_, se = NA_real_,
                                          ci_lower = NA_real_, ci_upper = NA_real_, p = NA_real_,
                                          model_structure = "guideline-recommended incident-case-weighted random-intercept/random-slope LMM; no ARMA",
                                          n_country_years = nrow(d),
                                          n_countries = length(unique(d$location_id)),
                                          status = conditionMessage(best)),
                      grid = data.frame(model = label, status = conditionMessage(best)),
                      fit = NULL, data = d, diagnostics = NULL))
      }
      tt <- summary(best)$coefficients
      singular_fit <- lme4::isSingular(best, tol = 1e-4)
      if (singular_fit) fit_warnings <- c(fit_warnings, "singular random-effects fit")
      effects <- exposure_terms
      if (use_shortage_cits) effects <- c("time_c:shortage", "shortage:post",
                                           "time_c:shortage:post")
      out <- do.call(rbind, lapply(effects, function(term) {
          if (!term %in% rownames(tt)) return(data.frame(term = term, estimate = NA,
                                                         se = NA, ci_lower = NA,
                                                         ci_upper = NA, p = NA))
          estimate <- unname(tt[term, "Estimate"])
          se <- unname(tt[term, "Std. Error"])
          dfree <- unname(tt[term, "df"])
          crit <- if (is.finite(dfree) && dfree > 0) stats::qt(.975, dfree) else 1.96
          data.frame(term = term, estimate = estimate, se = se,
                     ci_lower = estimate - crit*se, ci_upper = estimate + crit*se,
                     p = unname(tt[term, "Pr(>|t|)"]),
                     stringsAsFactors = FALSE)
      }))
      out$selected_start <- selected_start
      out$logLik <- as.numeric(logLik(best))
      out$logLik_gap_to_next_start <- loglik_gap
      out$model <- label
      out$model_structure <- "guideline-recommended incident-case-weighted LMM; random intercept/slope; no ARMA"
      out$n_country_years <- nrow(d)
      out$n_countries <- length(unique(d$location_id))
      out$status <- if (singular_fit) "singular_fit" else "ok"
      if (!use_shortage_cits && "ratio_within" %in% out$term) {
          iq <- stats::IQR(d$ratio_within, na.rm = TRUE)
          out$within_iqr <- ifelse(out$term == "ratio_within", iq, NA_real_)
          out$within_iqr_rate_change <- ifelse(out$term == "ratio_within",
                                               out$estimate * iq, NA_real_)
      }

      spec <- data.frame(model = label, n_country_years = nrow(d),
                         n_countries = length(unique(d$location_id)),
                         AIC = stats::AIC(best), BIC = stats::BIC(best),
                         selected_start = selected_start,
                         logLik = as.numeric(logLik(best)),
                         logLik_gap_to_next_start = loglik_gap,
                         model_structure = "guideline-recommended incident-case-weighted LMM; random intercept/slope; no ARMA",
                         singular_fit = singular_fit,
                         warning = if (length(fit_warnings)) paste(unique(fit_warnings), collapse = " | ") else "",
                         stringsAsFactors = FALSE)
      list(result = out, grid = spec, fit = best, data = d,
           diagnostics = NULL)
  }

  model_fits <- list()
  # Primary procurement specification: the log-transformed procurement ratio per
  # 1,000 guideline-recommended incident cases enters as a single exposure term. It is not
  # decomposed into within- and between-country components here; the decomposed
  # versions are retained below as a comparison.
  model_fits[["ratio_log_single"]] <- fit_lmm(
      model_data, "All countries: log procurement ratio per guideline-recommended patient (single term)",
      c("ratio_log"), include_case_rate = TRUE)
  classified <- model_data[!is.na(model_data$shortage), ]
  model_fits[["ratio_log_single_classified"]] <- fit_lmm(
      classified, "WHO-classified countries: log procurement ratio per guideline-recommended patient (single term)",
      c("ratio_log"), include_case_rate = TRUE)
  model_fits[["ratio_all_full"]] <- fit_lmm(
      model_data, "All countries: guideline-recommended ratio within/between",
      c("ratio_within", "ratio_between"), include_case_rate = TRUE)
  model_fits[["ratio_all_classified"]] <- fit_lmm(
      classified, "WHO-classified countries: guideline-recommended ratio within/between",
      c("ratio_within", "ratio_between"), include_case_rate = TRUE)
  model_fits[["who_shortage_classified"]] <- fit_lmm(
      classified, "Same classified rows: WHO binary shortage CITS",
      character(), include_case_rate = TRUE, use_shortage_cits = TRUE)
  # Denominator sensitivity: the same exposure built on all-age, both-sex incident
  # syphilis cases instead of the guideline-recommended population.
  model_fits[["allcase_denominator"]] <- fit_lmm(
      model_data[is.finite(model_data$allcase_ratio_within), ],
      "All countries: all-case denominator within/between",
      c("allcase_ratio_within", "allcase_ratio_between"),
      include_case_rate = TRUE)
  model_fits[["female_cases"]] <- fit_lmm(
      model_data[is.finite(model_data$female_ratio_within), ],
      "All countries: IQVIA per female syphilis cases within/between",
      c("female_ratio_within", "female_ratio_between"),
      include_case_rate = TRUE)
  model_fits[["dose_component"]] <- fit_lmm(
      model_data, "All countries: dose-equivalents component within/between",
      c("dose_within", "dose_between"), include_case_rate = TRUE)
  model_fits[["case_component"]] <- fit_lmm(
      model_data, "All countries: guideline-recommended case-count component within/between",
      c("cases_within", "cases_between"), include_case_rate = FALSE)

  # One-year lag is prespecified; no data-driven lag search.
  lag_data <- model_data
  lag_data <- lag_data[order(lag_data$location_id, lag_data$year), ]
  lag_key <- model_data[, c("location_id", "year", "ratio_per_1000_guidrecom")]
  lag_key$year <- lag_key$year + 1L
  names(lag_key)[3] <- "ratio_lag1"
  lag_data <- merge(lag_data, lag_key, by = c("location_id", "year"), all = FALSE)
  lag_data$ratio_lag1 <- as.numeric(lag_data$ratio_lag1)
  lag_data <- add_within_between(lag_data, "ratio_lag1", "ratio_lag1")
  model_fits[["ratio_lag1"]] <- fit_lmm(
      lag_data[is.finite(lag_data$ratio_lag1_within), ],
      "All countries: one-year lagged guideline-recommended ratio within/between",
      c("ratio_lag1_within", "ratio_lag1_between"),
      include_case_rate = TRUE)

  # Two-year lag sensitivity: assign each year's ratio to the same country's
  # outcome two years later, then estimate within- and between-country terms.
  lag2_data <- model_data
  lag2_key <- model_data[, c("location_id", "year", "ratio_per_1000_guidrecom")]
  lag2_key$year <- lag2_key$year + 2L
  names(lag2_key)[3] <- "ratio_lag2"
  lag2_data <- merge(lag2_data, lag2_key, by = c("location_id", "year"), all = FALSE)
  lag2_data$ratio_lag2 <- as.numeric(lag2_data$ratio_lag2)
  lag2_data <- add_within_between(lag2_data, "ratio_lag2", "ratio_lag2")
  model_fits[["ratio_lag2"]] <- fit_lmm(
      lag2_data[is.finite(lag2_data$ratio_lag2_within), ],
      "All countries: two-year lagged guideline-recommended ratio within/between",
      c("ratio_lag2_within", "ratio_lag2_between"),
      include_case_rate = TRUE)

  # Paired continuous and binary shortage models must use the same country-years.
  paired_keys <- intersect(
      paste(classified$location_id[is.finite(classified$ratio_within)],
            classified$year[is.finite(classified$ratio_within)]),
      paste(classified$location_id[!is.na(classified$shortage)],
            classified$year[!is.na(classified$shortage)]))
  paired <- classified[paste(classified$location_id, classified$year) %in% paired_keys, ]
  if (nrow(paired)) {
      model_fits[["ratio_binary_paired_continuous"]] <- fit_lmm(
          paired, "Paired WHO-classified sample: continuous guideline ratio",
          c("ratio_within", "ratio_between"), include_case_rate = FALSE)
      model_fits[["ratio_binary_paired_cits"]] <- fit_lmm(
          paired, "Paired WHO-classified sample: binary shortage CITS",
          character(), include_case_rate = TRUE, use_shortage_cits = TRUE)
  }

  bind_fill <- function(items) {
      items <- Filter(function(x) !is.null(x) && nrow(x) > 0, items)
      if (!length(items)) return(data.frame())
      columns <- unique(unlist(lapply(items, names)))
      items <- lapply(items, function(x) {
          missing <- setdiff(columns, names(x))
          for (nm in missing) x[[nm]] <- NA
          x[, columns, drop = FALSE]
      })
      do.call(rbind, items)
  }
  model_results <- bind_fill(lapply(model_fits, `[[`, "result"))
  model_grids <- bind_fill(lapply(model_fits, `[[`, "grid"))
  write_csv(model_results, "iqvia_incident_ratio_lmm_results.csv")
  write_csv(model_grids, "iqvia_incident_ratio_lmm_model_specifications.csv")
  # Exposure components and the case-rate adjustment can be mathematically
  # coupled; report their observed correlations rather than interpreting them
  # causally or as treatment coverage.
  coupling <- stats::cor(model_data[, c("ratio_within", "dose_within", "cases_within",
                                       "syph_rate_z")], use = "complete.obs")
  write_csv(as.data.frame(coupling), "iqvia_exposure_component_correlation.csv")
  writeLines(c(
      "R is BPG 2.4-MIU dose-equivalents per 1,000 guideline-recommended patients.",
      "Guideline-recommended patients are incident syphilis cases among women aged 15-49 years plus all children aged 0-14 years (GBD strata), the groups for whom BPG is the guideline-recommended treatment.",
      "The primary LMM uses R on the ln(1 + R) scale as a single exposure term and weights each country-year by guideline-recommended incident syphilis case counts.",
      "All-case and female-only denominator models are reported as denominator sensitivities.",
      "It is a relative procurement-intensity proxy, not a treatment-coverage measure or confirmed shortage.",
      "Only years with all four quarters and all selected product cells observed were included in primary longitudinal models.",
      "Blank/unobserved procurement entries were retained as missing; explicit all-product, all-quarter zeros were retained as verified zeros.",
      "2012 contains Q4 only and was excluded from annual comparisons and models.",
      "The product list and transfer factors define dose-equivalents, not treatment courses.",
      "All associations are ecological and non-causal; same-year procurement may reflect inventory, delayed use, or channels not captured by IQVIA."
  ), file.path(results_dir, "iqvia_incident_ratio_interpretation.txt"))

  message("v7 IQVIA incidence-case analyses complete. Eligible annual country-years: ",
          sum(volume_panel$annual_valid_for_primary, na.rm = TRUE),
          "; countries: ", length(unique(volume_panel$location_id[volume_panel$annual_valid_for_primary])),
          "; weighted LMMs fit without ARMA.")
})


# -----------------------------------------------------------------------------
# Describe procurement by shortage group and fit shortage CITS
# -----------------------------------------------------------------------------
local({
  # Describe log(1 + procurement ratio) in three shortage-status groups and fit an LMM.

  suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(lme4)
    library(lmerTest)
    library(patchwork)
  })
  results_dir <- file.path("results", "iqvia")
  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
  write_csv <- function(x, filename) {
    utils::write.csv(x, file.path(results_dir, filename), row.names = FALSE, na = "")
  }
  panel_file <- file.path(results_dir, "iqvia_incidence_case_ratio_panel.csv")
  if (!file.exists(panel_file)) {
      stop("Run the IQVIA reconstruction section of 04_additional_analyses.R first: panel not found.")
  }
  panel <- utils::read.csv(panel_file, stringsAsFactors = FALSE)

  group_levels <- c("Shortage reported", "No shortage reported",
                    "Unknown shortage status")
  group_labels <- c("shortage reported" = "Shortage reported",
                    "no shortage reported" = "No shortage reported",
                    "unclassified" = "Unknown shortage status")

  # ---- 1. Analysis frame ----------------------------------------------------
  # The descriptive frame requires a complete procurement numerator and a positive
  # guideline-recommended-population incident-syphilis denominator. It does not additionally require the
  # congenital-syphilis outcome or the covariates used by script 45; its sample
  # size can therefore exceed the outcome-association LMM sample.
  d <- panel
  d$group_long <- unname(group_labels[d$group])
  d <- d[d$year >= 2013 & d$year <= 2021 &
             !is.na(d$group_long) &
             is.finite(d$ratio_per_1000_guidrecom) &
             d$ratio_per_1000_guidrecom >= 0 &
             is.finite(d$syphilis_cases_guidrecom) & d$syphilis_cases_guidrecom > 0 &
             is.finite(d$syphilis_cases) & d$syphilis_cases > 0, ]
  if (!nrow(d)) stop("No country-year has a usable procurement ratio.")
  d$group_long <- factor(d$group_long, levels = group_levels)
  d$log_ratio <- log1p(d$ratio_per_1000_guidrecom)
  d$country <- d$location_name.x
  d <- d[order(d$group_long, d$country, d$year), ]

  write_csv(d[, c("location_id", "location_iso", "country", "year", "group",
                  "group_long", "annual_volume_dose_equivalents",
                  "syphilis_cases",
                  "syphilis_cases_guidrecom",
                  "cs_rate", "ratio_per_1000_guidrecom", "log_ratio")],
            "iqvia_log_ratio_panel.csv")

  # ---- 2. Mean-based descriptive statistics by group -----------------------
  # Summaries use means on the log scale and back-transform with exp(L)-1.
  # Case-weighted means and LMM weights use guideline-recommended GBD incident syphilis counts.
  # Confidence intervals are reserved for the LMM because country-years repeat.
  describe_group <- function(x, weights = NULL) {
      n <- length(x)
      m <- mean(x)
      s <- if (n > 1) stats::sd(x) else NA_real_
      back <- function(z) pmax(0, exp(z) - 1)
      out <- data.frame(
          n_country_years = n,
          mean_log1p = m,
          sd_log1p = s,
          arithmetic_mean_ratio = mean(expm1(x)),
          back_transformed_geometric_mean = back(m),
          stringsAsFactors = FALSE)
      if (!is.null(weights)) {
          out$total_guideline_recommended_incident_cases <- sum(weights)
          out$case_weighted_mean_log1p <- stats::weighted.mean(x, weights)
          out$case_weighted_arithmetic_mean_ratio <-
              stats::weighted.mean(expm1(x), weights)
          out$case_weighted_geometric_mean <-
              back(out$case_weighted_mean_log1p)
      }
      out
  }

  by_group <- do.call(rbind, lapply(group_levels, function(g) {
      dd <- d[d$group_long == g, ]
      if (!nrow(dd)) return(NULL)
      cbind(data.frame(group = g,
                       n_countries = length(unique(dd$location_id)),
                       year_min = min(dd$year), year_max = max(dd$year)),
            describe_group(dd$log_ratio, dd$syphilis_cases_guidrecom))
  }))
  by_group$case_weighted_geometric_mean_change_first_last <- NA_real_

  by_year <- do.call(rbind, lapply(split(d, list(d$group_long, d$year), drop = TRUE),
                                   function(dd) {
      cbind(data.frame(group = as.character(dd$group_long[1]), year = dd$year[1],
                       n_countries = length(unique(dd$location_id))),
            describe_group(dd$log_ratio, dd$syphilis_cases_guidrecom))
  }))
  by_year <- by_year[order(by_year$group, by_year$year), ]

  for (i in seq_len(nrow(by_group))) {
      g <- by_group$group[i]
      s <- by_year[by_year$group == g, ]
      if (nrow(s) >= 2) {
          first <- s$case_weighted_geometric_mean[s$year == min(s$year)][1]
          last <- s$case_weighted_geometric_mean[s$year == max(s$year)][1]
          by_group$case_weighted_geometric_mean_change_first_last[i] <-
              if (is.finite(first) && first > 0) last / first else NA_real_
      }
  }

  overall <- cbind(data.frame(group = "All groups",
                              n_countries = length(unique(d$location_id)),
                              year_min = min(d$year), year_max = max(d$year)),
                   describe_group(d$log_ratio, d$syphilis_cases_guidrecom))
  overall$case_weighted_geometric_mean_change_first_last <- NA_real_
  overall <- overall[, names(by_group), drop = FALSE]
  write_csv(rbind(by_group, overall), "iqvia_log_ratio_group_description.csv")
  write_csv(by_year, "iqvia_log_ratio_group_year.csv")

  # ---- 3. LMM group comparison; no ARMA ------------------------------------
  # The shortage-reporting group is the reference. Calendar year is categorical
  # to avoid imposing a linear procurement trend. Guideline-recommended GBD incident syphilis
  # cases are both the denominator of R and the prior weights. The country random
  # intercept accounts for repeated years.
  d$location_id <- factor(d$location_id)
  d$year_factor <- factor(d$year)
  # Rescale to a mean of 1 instead of dividing by the sum: sum-normalised weights
  # are ~1e-5, which drives the lme4 residual variance to the boundary and
  # produces spurious singular fits.
  d$prior_weight <- d$syphilis_cases_guidrecom / mean(d$syphilis_cases_guidrecom)
  fit_warnings <- character()
  lmm <- tryCatch(withCallingHandlers(
      lmerTest::lmer(log_ratio ~ group_long + year_factor + (1 | location_id),
                     data = d, weights = prior_weight, REML = FALSE,
                     na.action = na.fail),
      warning = function(w) {
          fit_warnings <<- c(fit_warnings, conditionMessage(w))
          invokeRestart("muffleWarning")
      }), error = function(e) e)

  lmm_rows <- data.frame()
  if (inherits(lmm, "merMod")) {
      tt <- summary(lmm)$coefficients
      keep <- grep("^group_long", rownames(tt), value = TRUE)
      singular_fit <- lme4::isSingular(lmm, tol = 1e-4)
      if (singular_fit) fit_warnings <- c(fit_warnings, "singular random-effects fit")
      lmm_rows <- do.call(rbind, lapply(keep, function(term) {
          est <- unname(tt[term, "Estimate"])
          se <- unname(tt[term, "Std. Error"])
          dfree <- unname(tt[term, "df"])
          crit <- if (is.finite(dfree) && dfree > 0) stats::qt(.975, dfree) else 1.96
          comparison <- if (grepl("No shortage", term)) {
              "No shortage reported vs shortage reported"
          } else {
              "Unknown shortage status vs shortage reported"
          }
          data.frame(
              comparison = comparison,
              term = term,
              estimate_log1p = est,
              se = se,
              df = dfree,
              ci_lower_log1p = est - crit * se,
              ci_upper_log1p = est + crit * se,
              p_value = unname(tt[term, "Pr(>|t|)"]),
              ratio_of_geometric_1_plus_R = exp(est),
              ratio_ci_lower = exp(est - crit * se),
              ratio_ci_upper = exp(est + crit * se),
              n_country_years = nrow(d),
              n_countries = length(unique(d$location_id)),
              singular_fit = singular_fit,
              model = "Guideline-recommended incident-case-weighted LMM; categorical year; country random intercept; no ARMA",
              stringsAsFactors = FALSE)
      }))
      model_spec <- data.frame(
          model = "log1p(R) ~ shortage status + factor(year) + (1 | country)",
          transformation = "R = BPG 2.4-MIU dose-equivalents per 1,000 guideline-recommended incident syphilis cases (women aged 15-49 years plus children aged 0-14 years); outcome = ln(1 + R)",
          weighting = "Estimated guideline-recommended GBD incident syphilis cases; relative prior weights rescaled to mean 1",
          estimation = "Maximum likelihood",
          residual_correlation = "No ARMA structure",
          n_country_years = nrow(d),
          n_countries = length(unique(d$location_id)),
          singular_fit = singular_fit,
          warning = if (length(fit_warnings)) paste(unique(fit_warnings), collapse = " | ") else "",
          AIC = stats::AIC(lmm), BIC = stats::BIC(lmm),
          stringsAsFactors = FALSE)
      cat("LMM on log1p procurement ratio fitted; n =", nrow(d),
          "country-years,", length(unique(d$location_id)), "countries; singular =",
          singular_fit, "\n")
  } else {
      model_spec <- data.frame(
          model = "log1p(R) ~ shortage status + factor(year) + (1 | country)",
          transformation = "R = BPG 2.4-MIU dose-equivalents per 1,000 guideline-recommended incident syphilis cases (women aged 15-49 years plus children aged 0-14 years); outcome = ln(1 + R)",
          weighting = "Estimated guideline-recommended GBD incident syphilis cases; relative prior weights rescaled to mean 1",
          estimation = "Maximum likelihood",
          residual_correlation = "No ARMA structure",
          n_country_years = nrow(d),
          n_countries = length(unique(d$location_id)),
          singular_fit = NA,
          warning = conditionMessage(lmm), AIC = NA, BIC = NA,
          stringsAsFactors = FALSE)
      cat("LMM not fitted:", conditionMessage(lmm), "\n")
  }
  write_csv(lmm_rows, "iqvia_log_ratio_group_lmm.csv")
  write_csv(model_spec, "iqvia_log_ratio_group_lmm_specification.csv")

  # Report all pairwise shortage-status contrasts using the same Satterthwaite
  # reference distribution for each estimate, interval and P value.
  fixed_terms <- names(lme4::fixef(lmm))
  contrast_matrix <- matrix(0, nrow = 3, ncol = length(fixed_terms),
    dimnames = list(c("No shortage vs shortage", "Unknown vs shortage",
                      "Unknown vs no shortage"), fixed_terms))
  contrast_matrix[1, "group_longNo shortage reported"] <- 1
  contrast_matrix[2, "group_longUnknown shortage status"] <- 1
  contrast_matrix[3, "group_longUnknown shortage status"] <- 1
  contrast_matrix[3, "group_longNo shortage reported"] <- -1
  contrast_rows <- list()
  for (i in seq_len(nrow(contrast_matrix))) {
    contrast <- lmerTest::contest1D(lmm, contrast_matrix[i, ], ddf = "Satterthwaite",
                                    confint = TRUE)
    contrast_rows[[i]] <- data.frame(
      contrast = rownames(contrast_matrix)[i], estimate = contrast$Estimate,
      standard_error = contrast[["Std. Error"]], df = contrast$df,
      ci_lower = contrast$lower, ci_upper = contrast$upper,
      p_value = contrast[["Pr(>|t|)" ]], countries = length(unique(d$location_id)),
      country_years = nrow(d), period = "2013-2021", stringsAsFactors = FALSE)
  }
  write_csv(do.call(rbind, contrast_rows), "iqvia_shortage_group_all_contrasts.csv")
  # ---- 4. Figures -----------------------------------------------------------
  group_cols <- c("Shortage reported" = "#C63D40",
                  "No shortage reported" = "#176B9A",
                  "Unknown shortage status" = "#737A83")

  group_mean_points <- by_group
  group_mean_points$group_long <- factor(group_mean_points$group, levels = group_levels)

  dist_plot <- ggplot(d, aes(x = group_long, y = log_ratio, colour = group_long)) +
      geom_jitter(width = 0.14, height = 0, size = 1.5, alpha = 0.65) +
      geom_point(data = group_mean_points,
                 aes(x = group_long, y = case_weighted_mean_log1p,
                     colour = group_long),
                 inherit.aes = FALSE, shape = 23, size = 3.2,
                 fill = "white", stroke = 0.8) +
      scale_x_discrete(labels = c("Shortage\nreported", "No shortage\nreported",
                                  "Unknown\nstatus")) +
      scale_colour_manual(values = group_cols, guide = "none") +
      labs(title = "A  Country-year distribution",
           x = NULL,
           y = "ln[1 + BPG dose equivalents per\n1,000 guideline-recommended incident cases]") +
      theme_bw(base_size = 10) +
      theme(panel.grid.minor = element_blank(),
            plot.title = element_text(face = "bold", size = 10.5))

  annual <- by_year
  annual$group_long <- factor(annual$group, levels = group_levels)

  traj_plot <- ggplot() +
      geom_line(data = d, aes(x = year, y = log_ratio,
                              group = location_id, colour = group_long),
                linewidth = 0.3, alpha = 0.35) +
      geom_point(data = d, aes(x = year, y = log_ratio,
                               colour = group_long), size = 0.7, alpha = 0.4) +
      geom_line(data = annual, aes(x = year, y = case_weighted_mean_log1p,
                                   colour = group_long), linewidth = 1.05) +
      geom_point(data = annual, aes(x = year, y = case_weighted_mean_log1p,
                                    colour = group_long), size = 2.1) +
      scale_x_continuous(breaks = 2013:2021) +
      scale_colour_manual(values = group_cols) +
      labs(title = "B  Annual country trajectories",
           x = "Year",
           y = "Guideline-recommended-case-weighted mean ln[1 + dose equivalents per\n1,000 guideline-recommended incident cases]",
           colour = NULL) +
      theme_bw(base_size = 10) +
      theme(panel.grid.minor = element_blank(),
            plot.title = element_text(face = "bold", size = 10.5),
            legend.position = "top")

  if (nrow(lmm_rows)) {
      forest_data <- lmm_rows
      forest_data$label <- ifelse(grepl("No shortage", forest_data$comparison),
                                  "No shortage vs shortage",
                                  "Unknown vs shortage")
      forest_data$label <- factor(forest_data$label,
                                  levels = rev(c("No shortage vs shortage",
                                                 "Unknown vs shortage")))
      forest_plot <- ggplot(forest_data,
                            aes(x = ratio_of_geometric_1_plus_R, y = label)) +
          geom_vline(xintercept = 1, linetype = "dashed", colour = "#8B929B",
                     linewidth = 0.4) +
          geom_errorbar(aes(xmin = ratio_ci_lower, xmax = ratio_ci_upper),
                        orientation = "y", width = 0.14, linewidth = 0.6) +
          geom_point(size = 2.4, colour = "#24282D") +
          scale_x_log10(labels = scales::label_number(accuracy = 0.1)) +
          labs(title = "C  LMM group contrasts",
               x = "Ratio of geometric means of (1 + R)\nrelative to shortage reported",
               y = NULL) +
          theme_bw(base_size = 10) +
          theme(panel.grid.minor = element_blank(),
                plot.title = element_text(face = "bold", size = 10.5),
                axis.text.y = element_text(size = 9))
  } else {
      forest_plot <- patchwork::plot_spacer()
  }

  combined <- (dist_plot + traj_plot + forest_plot) +
      patchwork::plot_layout(widths = c(1, 1.25, 0.9), guides = "keep")
  ggsave(file.path(results_dir, "iqvia_log_ratio_by_group.png"), combined,
         width = 13.5, height = 5, dpi = 300, bg = "white")
  ggsave(file.path(results_dir, "iqvia_log_ratio_by_group.pdf"), combined,
         width = 13.5, height = 5, device = cairo_pdf, bg = "white")

  message("Log-ratio description and LMM complete. Eligible country-years: ", nrow(d),
          "; countries: ", length(unique(d$location_id)),
          "; groups: ", paste(sprintf("%s n=%d", group_levels,
                                      as.integer(table(d$group_long))), collapse = "; "))




  # -----------------------------------------------------------------------------
  # Shortage-status procurement interrupted time-series analysis
  # -----------------------------------------------------------------------------

  # Shortage-status procurement interrupted time-series analysis
  # Test whether the post-2016 trend in the procurement ratio differed by
  # reported BPG-shortage status. This is separate from the CS-outcome association.

  library(lme4)
  library(lmerTest)
  library(ggplot2)

  results_dir <- file.path("results", "iqvia")
  figure_dir <- file.path("figures", "iqvia")
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

  # Start from the product-quarter audited country-year panel.
  panel <- read.csv(file.path(results_dir, "iqvia_incidence_case_ratio_panel.csv"),
                    stringsAsFactors = FALSE)
  panel$location_id <- as.character(panel$location_id)
  covariates <- as.data.frame(cs_panel)[, c("location_id", "year", "shortage",
    "hdi", "sdi", "gni", "gbd_superregion")]
  covariates$location_id <- as.character(covariates$location_id)
  panel <- merge(panel, covariates, by = c("location_id", "year"),
                 all.x = TRUE, sort = FALSE)

  # Use 2013-2021 complete IQVIA annual ratios and the 89-country shortage labels.
  data <- panel[panel$year >= 2013 & panel$year <= 2021 &
    panel$shortage %in% c(0, 1) &
    is.finite(panel$ratio_per_1000_guidrecom) &
    panel$ratio_per_1000_guidrecom >= 0 &
    is.finite(panel$syphilis_cases_guidrecom) &
    panel$syphilis_cases_guidrecom > 0, ]

  data$time_c <- data$year - 2016
  data$post <- as.integer(data$year >= 2016)
  data$log_ratio <- log1p(data$ratio_per_1000_guidrecom)
  data$location_id <- factor(data$location_id)
  data$gbd_superregion <- factor(data$gbd_superregion)
  data$sdi_z <- as.numeric(scale(data$sdi))
  data$hdi_z <- as.numeric(scale(data$hdi))
  data$gni_z <- as.numeric(scale(data$gni))

  # Syphilis prevalence is not adjusted for in a model whose ratio denominator is
  # incident syphilis cases. Fit one socioeconomic index at a time with region.
  complete_sdi <- complete.cases(data[, c("log_ratio", "time_c", "post", "shortage",
    "sdi_z", "gbd_superregion", "location_id", "syphilis_cases_guidrecom")])
  complete_hdi <- complete.cases(data[, c("log_ratio", "time_c", "post", "shortage",
    "hdi_z", "gbd_superregion", "location_id", "syphilis_cases_guidrecom")])
  complete_gni <- complete.cases(data[, c("log_ratio", "time_c", "post", "shortage",
    "gni_z", "gbd_superregion", "location_id", "syphilis_cases_guidrecom")])

  sample_audit <- data.frame(
    specification = c("SDI + region", "HDI + region", "GNI + region"),
    country_years = c(sum(complete_sdi), sum(complete_hdi), sum(complete_gni)),
    countries = c(length(unique(data$location_id[complete_sdi])),
                  length(unique(data$location_id[complete_hdi])),
                  length(unique(data$location_id[complete_gni]))),
    stringsAsFactors = FALSE
  )
  write.csv(sample_audit, file.path(results_dir, "iqvia_shortage_cits_sample_audit.csv"),
            row.names = FALSE)

  # The complete-row count is expected to be the same for these preassembled
  # socioeconomic series. Preserve an exact common sample for fair comparisons.
  common_complete <- complete_sdi & complete_hdi & complete_gni
  analysis_data <- data[common_complete, ]
  analysis_data$location_id <- droplevels(analysis_data$location_id)
  analysis_data$gbd_superregion <- droplevels(analysis_data$gbd_superregion)
  analysis_data$prior_weight <- analysis_data$syphilis_cases_guidrecom /
    mean(analysis_data$syphilis_cases_guidrecom)

  # Fit the primary random-intercept/random-slope model and its same-sample
  # random-intercept-only stability model. Neither uses ARMA.
  formula_slope <- log_ratio ~ time_c * shortage * post + sdi_z + gbd_superregion +
    (1 + time_c | location_id)
  fit_slope <- lmerTest::lmer(formula_slope, data = analysis_data,
    weights = prior_weight, REML = FALSE,
    control = lme4::lmerControl(optimizer = "bobyqa"), na.action = na.fail)

  formula_intercept <- log_ratio ~ time_c * shortage * post + sdi_z + gbd_superregion +
    (1 | location_id)
  fit_intercept <- lmerTest::lmer(formula_intercept, data = analysis_data,
    weights = prior_weight, REML = FALSE,
    control = lme4::lmerControl(optimizer = "bobyqa"), na.action = na.fail)

  fit_list <- list("Random intercept and time slope" = fit_slope,
                   "Random intercept only" = fit_intercept)
  target_terms <- c("time_c:shortage", "shortage:post", "time_c:shortage:post")
  target_labels <- c("Pre-2016 slope difference", "2016 level-change difference",
                     "Post-2016 slope-change difference")
  result_rows <- list()
  for (fit_name in names(fit_list)) {
    fit <- fit_list[[fit_name]]
    coefficient_table <- summary(fit)$coefficients
    for (i in seq_along(target_terms)) {
      term <- target_terms[i]
      estimate <- coefficient_table[term, "Estimate"]
      standard_error <- coefficient_table[term, "Std. Error"]
      degrees_freedom <- coefficient_table[term, "df"]
      critical_value <- qt(0.975, degrees_freedom)
      result_rows[[length(result_rows) + 1L]] <- data.frame(
        term = term, interpretation = target_labels[i],
        random_effects = fit_name, estimate = estimate,
        standard_error = standard_error, df = degrees_freedom,
        ci_lower = estimate - critical_value * standard_error,
        ci_upper = estimate + critical_value * standard_error,
        p_value = coefficient_table[term, "Pr(>|t|)"],
        country_years = nrow(analysis_data),
        countries = length(unique(analysis_data$location_id)),
        shortage_countries = length(unique(analysis_data$location_id[analysis_data$shortage == 1])),
        no_shortage_countries = length(unique(analysis_data$location_id[analysis_data$shortage == 0])),
        singular_fit = lme4::isSingular(fit),
        AIC = AIC(fit), BIC = BIC(fit),
        formula = paste(deparse(formula(fit)), collapse = " "),
        covariates = "SDI; GBD super-region; no syphilis prevalence",
        weights = "GBD incident syphilis cases in women aged 15-49 plus children aged 0-14; rescaled to mean 1",
        residual_correlation = "No ARMA",
        stringsAsFactors = FALSE
      )
    }
  }
  results <- do.call(rbind, result_rows)
  write.csv(results, file.path(results_dir, "iqvia_shortage_cits_lmm.csv"), row.names = FALSE)

  # Plot CITS interaction estimates for the two random-effects structures.
  results$interpretation <- factor(results$interpretation,
    levels = rev(target_labels))
  plot <- ggplot(results, aes(x = estimate, y = interpretation,
    colour = random_effects, shape = random_effects)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey45") +
    geom_errorbar(aes(xmin = ci_lower, xmax = ci_upper), orientation = "y",
      width = 0.15, position = position_dodge(width = 0.55)) +
    geom_point(position = position_dodge(width = 0.55), size = 2.4) +
    labs(x = "Difference in log(1 + procurement ratio), with 95% CI",
         y = NULL, colour = NULL, shape = NULL,
         title = "BPG shortage CITS models for procurement intensity",
         subtitle = paste(nrow(analysis_data), "country-years;",
           length(unique(analysis_data$location_id)), "countries; weighted by guideline-recommended incident cases")) +
    theme_bw(base_size = 10) + theme(legend.position = "bottom")
  ggsave(file.path(figure_dir, "iqvia_shortage_cits_forest.png"), plot,
         width = 9, height = 5.5, dpi = 300)
  ggsave(file.path(figure_dir, "iqvia_shortage_cits_forest.pdf"), plot,
         width = 9, height = 5.5)

  cat("Saved IQVIA shortage CITS LMMs on", nrow(analysis_data),
      "country-years from", length(unique(analysis_data$location_id)), "countries.\n")
})


# -----------------------------------------------------------------------------
# Procurement–incidence association and model-structure sensitivities
# -----------------------------------------------------------------------------
local({
  # Test the association between procurement intensity and congenital syphilis.
  # Syphilis prevalence and the incidence rate are excluded from this denominator-linked model.

  suppressPackageStartupMessages({ library(lmerTest); library(lme4); library(ggplot2) })
  results_dir <- file.path("results", "iqvia")
  model_control <- lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 200000))
  load_covariates <- function() {
    columns <- c("location_id", "year", "shortage", "sdi", "hdi", "gni", "gbd_superregion")
    out <- as.data.frame(cs_panel)[, columns]
    out$location_id <- as.character(out$location_id)
    out
  }
  results_dir <- file.path("results", "iqvia")
  panel <- utils::read.csv(
      file.path(results_dir, "iqvia_incidence_case_ratio_panel.csv"),
      stringsAsFactors = FALSE)
  cov <- load_covariates()
  panel$location_id <- as.character(panel$location_id)
  cov$location_id <- as.character(cov$location_id)
  d0 <- merge(panel, cov, by = c("location_id", "year"),
              all = FALSE, sort = FALSE)
  d0 <- d0[d0$year >= 2013 & d0$year <= 2021 &
               is.finite(d0$ratio_per_1000_guidrecom) &
               d0$ratio_per_1000_guidrecom >= 0 &
               is.finite(d0$syphilis_cases_guidrecom) &
               d0$syphilis_cases_guidrecom > 0 &
               is.finite(d0$cs_rate) &
               is.finite(d0$population) & d0$population > 0, , drop = FALSE]

  z_score <- function(x) {
      x <- as.numeric(x)
      s <- stats::sd(x, na.rm = TRUE)
      if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
      (x - mean(x, na.rm = TRUE)) / s
  }
  for (v in c("hdi", "sdi", "gni")) d0[[paste0(v, "_z")]] <- z_score(d0[[v]])
  d0$ratio_log <- log1p(d0$ratio_per_1000_guidrecom)
  d0$year_factor <- factor(d0$year)
  d0$year_c <- d0$year - 2016
  d0$location_id <- as.character(d0$location_id)
  d0$gbd_superregion <- factor(d0$gbd_superregion)

  eligible_rows <- nrow(d0)
  eligible_countries <- length(unique(d0$location_id))
  sdi_complete <- complete.cases(d0[, c("cs_rate", "ratio_log",
      "syphilis_cases_guidrecom", "sdi_z", "gbd_superregion")])
  common_complete <- complete.cases(d0[, c("cs_rate", "ratio_log",
      "syphilis_cases_guidrecom", "hdi_z", "sdi_z", "gni_z",
      "gbd_superregion")])
  sample_audit <- data.frame(
      sample = c("Eligible outcome/exposure/weight rows",
                 "SDI + region complete cases",
                 "Common SES comparison sample: HDI + SDI + GNI + region"),
      country_years = c(eligible_rows, sum(sdi_complete), sum(common_complete)),
      countries = c(eligible_countries,
          length(unique(d0$location_id[sdi_complete])),
          length(unique(d0$location_id[common_complete]))),
      stringsAsFactors = FALSE)
  utils::write.csv(sample_audit,
      file.path(results_dir, "iqvia_no_syph_prevalence_sample_audit.csv"),
      row.names = FALSE, na = "")

  # The preferred SDI + region model uses all country-years complete for those
  # required covariates. Other SES comparisons use a common sample complete for
  # HDI, SDI, GNI, and region.
  make_model_data <- function(mask) {
      dd <- d0[mask, , drop = FALSE]
      dd$location_id <- droplevels(factor(dd$location_id))
      dd$gbd_superregion <- droplevels(factor(dd$gbd_superregion))
      dd$prior_weight <- dd$syphilis_cases_guidrecom /
          mean(dd$syphilis_cases_guidrecom)
      dd$ratio_log_within <- dd$ratio_log -
          ave(dd$ratio_log, dd$location_id, FUN = mean)
      dd$ratio_log_between <- ave(dd$ratio_log, dd$location_id, FUN = mean)
      dd[order(dd$location_id, dd$year), , drop = FALSE]
  }
  d_primary <- make_model_data(sdi_complete)
  d <- make_model_data(common_complete)
  if (nrow(d_primary) < 30L || length(unique(d_primary$location_id)) < 8L)
      stop("Insufficient SDI + region complete-case sample.")
  if (nrow(d) < 30L || length(unique(d$location_id)) < 8L)
      stop("Insufficient common SES comparison sample.")

  specs <- list(
      "SDI + region" = c("sdi_z", "gbd_superregion"),
      "HDI + region" = c("hdi_z", "gbd_superregion"),
      "GNI + region" = c("gni_z", "gbd_superregion"),
      "Joint SES + region" = c("hdi_z", "sdi_z", "gni_z", "gbd_superregion"),
      "Region only" = "gbd_superregion"
  )
  theta_starts <- list(
      default = NULL,
      start_a = c(150.048207313949, 28.5688389337839, 2.53797203789327),
      start_b = c(62.3299316056475, 56.5044512309326, 15.5967255963415)
  )

  fit_spec <- function(label, covars, exposure_terms, exposure_label, data_used = d) {
      dd <- data_used[order(data_used$location_id, data_used$year), , drop = FALSE]
      formula <- stats::as.formula(paste(
          "cs_rate ~ factor(year) +",
          paste(c(exposure_terms, covars, "(1 + year_c | location_id)"),
                collapse = " + ")))
      fits <- list()
      warning_text <- character()
      for (start_name in names(theta_starts)) {
          theta <- theta_starts[[start_name]]
          fit_try <- tryCatch(withCallingHandlers(
              if (is.null(theta)) {
                  lmerTest::lmer(formula, data = dd, weights = prior_weight,
                      REML = FALSE, control = model_control, na.action = na.fail)
              } else {
                  lmerTest::lmer(formula, data = dd, weights = prior_weight,
                      REML = FALSE, control = model_control, na.action = na.fail,
                      start = list(theta = theta))
              },
              warning = function(w) {
                  warning_text <<- c(warning_text, conditionMessage(w))
                  invokeRestart("muffleWarning")
              }), error = function(e) e)
          if (inherits(fit_try, "merMod")) fits[[start_name]] <- fit_try
          else warning_text <- c(warning_text,
                                 paste(start_name, conditionMessage(fit_try)))
      }
      if (!length(fits)) {
          return(data.frame(model = label, exposure = exposure_label,
              term = exposure_terms, estimate = NA_real_, std_error = NA_real_,
              ci_lower = NA_real_, ci_upper = NA_real_, p_value = NA_real_,
              selected_start = NA_character_, n_successful_starts = 0L,
              estimate_min_across_starts = NA_real_,
              estimate_max_across_starts = NA_real_, logLik = NA_real_,
              logLik_gap_to_next_start = NA_real_, n_country_years = nrow(dd),
              n_countries = length(unique(dd$location_id)), singular_fit = NA,
              AIC = NA_real_, warning = paste(warning_text, collapse = " | "),
              stringsAsFactors = FALSE))
      }
      ll <- vapply(fits, function(x) as.numeric(logLik(x)), numeric(1))
      best_i <- which.max(ll)
      fit <- fits[[best_i]]
      tt <- summary(fit)$coefficients
      singular <- lme4::isSingular(fit, tol = 1e-4)
      if (singular) warning_text <- c(warning_text, "selected fit is singular")
      ll_sorted <- sort(ll, decreasing = TRUE)
      out <- do.call(rbind, lapply(seq_along(exposure_terms), function(i) {
          term <- exposure_terms[i]
          estimate <- unname(tt[term, "Estimate"])
          se <- unname(tt[term, "Std. Error"])
          dfree <- unname(tt[term, "df"])
          crit <- if (is.finite(dfree) && dfree > 0) qt(.975, dfree) else 1.96
          start_est <- vapply(fits, function(x) unname(fixef(x)[term]), numeric(1))
          data.frame(
              model = label, exposure = exposure_label, term = term,
              estimate = estimate, std_error = se,
              ci_lower = estimate - crit * se,
              ci_upper = estimate + crit * se,
              p_value = unname(tt[term, "Pr(>|t|)"]),
              selected_start = names(fits)[best_i],
              n_successful_starts = length(fits),
              estimate_min_across_starts = min(start_est),
              estimate_max_across_starts = max(start_est),
              logLik = max(ll),
              logLik_gap_to_next_start = if (length(ll_sorted) > 1L)
                  ll_sorted[1] - ll_sorted[2] else NA_real_,
              n_country_years = nrow(dd),
              n_countries = length(unique(dd$location_id)),
              singular_fit = singular, AIC = AIC(fit),
              warning = paste(unique(warning_text), collapse = " | "),
              stringsAsFactors = FALSE)
      }))
      out
  }

  results <- rbind(
      fit_spec("SDI + region (full eligible sample)",
          specs[["SDI + region"]], "ratio_log", "Single ln(1 + P) term",
          data_used = d_primary),
      fit_spec("SDI + region (common SES sample)",
          specs[["SDI + region"]], "ratio_log", "Single ln(1 + P) term",
          data_used = d),
      fit_spec("HDI + region (common SES sample)",
          specs[["HDI + region"]], "ratio_log", "Single ln(1 + P) term",
          data_used = d),
      fit_spec("GNI + region (common SES sample)",
          specs[["GNI + region"]], "ratio_log", "Single ln(1 + P) term",
          data_used = d),
      fit_spec("Joint SES + region (common SES sample)",
          specs[["Joint SES + region"]], "ratio_log", "Single ln(1 + P) term",
          data_used = d),
      fit_spec("Region only (full eligible sample)",
          specs[["Region only"]], "ratio_log", "Single ln(1 + P) term",
          data_used = d_primary),
      fit_spec("SDI + region (full eligible sample)",
          specs[["SDI + region"]],
          c("ratio_log_within", "ratio_log_between"),
          "Country-mean decomposition", data_used = d_primary)
  )
  utils::write.csv(results,
      file.path(results_dir, "iqvia_covariate_sensitivity.csv"),
      row.names = FALSE, na = "")

  ses_cor <- stats::cor(d[, c("hdi_z", "sdi_z", "gni_z")])
  utils::write.csv(ses_cor,
      file.path(results_dir, "iqvia_covariate_correlation.csv"))

  vif_rows <- lapply(c("hdi_z", "sdi_z", "gni_z"), function(target) {
      others <- setdiff(c("hdi_z", "sdi_z", "gni_z"), target)
      form <- stats::reformulate(c(others, "factor(year)", "gbd_superregion"),
                                 response = target)
      aux <- stats::lm(form, data = d, weights = prior_weight)
      r2 <- summary(aux)$r.squared
      data.frame(variable = target, weighted_r_squared = r2,
                 weighted_vif = 1 / (1 - r2), n = nrow(d),
                 stringsAsFactors = FALSE)
  })
  utils::write.csv(do.call(rbind, vif_rows),
      file.path(results_dir, "iqvia_covariate_vif.csv"),
      row.names = FALSE, na = "")

  single_results <- results[results$exposure == "Single ln(1 + P) term" &
      is.finite(results$estimate) & is.finite(results$ci_lower) &
      is.finite(results$ci_upper), , drop = FALSE]
  plot_labels <- c(
      "SDI + region (full eligible sample)" = "SDI + region (112 country-years / 23 countries)",
      "SDI + region (common SES sample)" = "SDI + region (95 country-years / 21 countries)",
      "HDI + region (common SES sample)" = "HDI + region (95 country-years / 21 countries)",
      "GNI + region (common SES sample)" = "GNI + region (95 country-years / 21 countries)",
      "Joint SES + region (common SES sample)" = "HDI + SDI + GNI + region (95 country-years / 21 countries)",
      "Region only (full eligible sample)" = "Region only (112 country-years / 23 countries)"
  )
  single_results$label <- unname(plot_labels[single_results$model])
  single_results$label <- factor(single_results$label,
      levels = single_results$label[order(single_results$estimate)])
  forest <- ggplot2::ggplot(single_results,
      ggplot2::aes(x = estimate, y = label)) +
      ggplot2::geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
      ggplot2::geom_segment(ggplot2::aes(x = ci_lower, xend = ci_upper,
                                         yend = label), linewidth = .6) +
      ggplot2::geom_point(size = 2.2, colour = "#1F4E79") +
      ggplot2::labs(
          title = "No-prevalence covariate sensitivity",
          subtitle = paste0(nrow(d_primary), " country-years from ",
              length(unique(d_primary$location_id)),
              " countries; weights proportional to guideline-recommended cases"),
          x = "Change in congenital syphilis rate per 1,000 live births
  per one-unit increase in ln(1 + procurement ratio)",
          y = NULL
      ) +
      ggplot2::theme_bw(base_size = 10) +
      ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                     axis.text.y = ggplot2::element_text(size = 9))
  ggplot2::ggsave(file.path(results_dir,
      "iqvia_covariate_sensitivity_forest.png"), forest,
      width = 9.5, height = 5.5, dpi = 300)
  ggplot2::ggsave(file.path(results_dir,
      "iqvia_covariate_sensitivity_forest.pdf"), forest,
      width = 9.5, height = 5.5)

  message("No-prevalence LMM sensitivity complete: preferred model uses ",
      nrow(d_primary), " country-years from ",
      length(unique(d_primary$location_id)),
      " countries; SES comparisons use ", nrow(d), " country-years from ",
      length(unique(d$location_id)),
      " countries. All models use guideline-recommended incident-case weights.")




  # -----------------------------------------------------------------------------
  # Model-structure sensitivity analysis
  # -----------------------------------------------------------------------------

  # Model-structure sensitivity analysis
  # Model-structure sensitivity for the no-syphilis-prevalence procurement LMM.
  # R/18 creates the shared eligible complete-case sample and case weights.

  suppressPackageStartupMessages(library(lmerTest))

  covariate_sets <- list(
      "SDI + region (full eligible sample)" = c("sdi_z", "gbd_superregion"),
      "HDI + region" = c("hdi_z", "gbd_superregion"),
      "GNI + region" = c("gni_z", "gbd_superregion"),
      "Joint SES + region" = c("hdi_z", "sdi_z", "gni_z", "gbd_superregion")
  )
  random_structures <- list(
      "Correlated intercept and slope" = "(1 + year_c | location_id)",
      "Independent intercept and slope" = "(1 + year_c || location_id)",
      "Random intercept only" = "(1 | location_id)"
  )

  fit_multistart <- function(formula, structure, data_used) {
      starts <- if (structure == "Correlated intercept and slope") {
          list(default = NULL,
               start_a = c(150.048207313949, 28.5688389337839, 2.53797203789327),
               start_b = c(62.3299316056475, 56.5044512309326, 15.5967255963415))
      } else if (structure == "Independent intercept and slope") {
          list(default = NULL, start_a = c(150.0482, 2.538),
               start_b = c(62.3299, 15.5967))
      } else {
          list(default = NULL, start_a = 150.0482, start_b = 62.3299)
      }
      fits <- list()
      warns <- character()
      for (nm in names(starts)) {
          theta <- starts[[nm]]
          fit <- tryCatch(withCallingHandlers(
              if (is.null(theta)) {
                  lmerTest::lmer(formula, data = data_used, weights = data_used$prior_weight,
                      REML = FALSE, control = model_control, na.action = na.fail)
              } else {
                  lmerTest::lmer(formula, data = data_used, weights = data_used$prior_weight,
                      REML = FALSE, control = model_control, na.action = na.fail,
                      start = list(theta = theta))
              },
              warning = function(w) {
                  warns <<- c(warns, conditionMessage(w))
                  invokeRestart("muffleWarning")
              }), error = function(e) e)
          if (inherits(fit, "merMod")) fits[[nm]] <- fit
          else warns <- c(warns, paste(nm, conditionMessage(fit)))
      }
      if (!length(fits)) return(list(fit = NULL, warnings = warns))
      ll <- vapply(fits, function(x) as.numeric(logLik(x)), numeric(1))
      best_i <- which.max(ll)
      list(fit = fits[[best_i]], fits = fits, logLik = ll,
           selected_start = names(fits)[best_i], warnings = warns)
  }

  out <- list()
  for (set_name in names(covariate_sets)) {
      for (structure_name in names(random_structures)) {
          form <- stats::as.formula(paste(
              "cs_rate ~ factor(year) + ratio_log +",
              paste(c(covariate_sets[[set_name]],
                      random_structures[[structure_name]]), collapse = " + ")))
          data_used <- if (set_name == "SDI + region (full eligible sample)") d_primary else d
          fitted <- fit_multistart(form, structure_name, data_used)
          if (is.null(fitted$fit)) next
          fit <- fitted$fit
          tt <- summary(fit)$coefficients
          estimate <- unname(tt["ratio_log", "Estimate"])
          se <- unname(tt["ratio_log", "Std. Error"])
          dfree <- unname(tt["ratio_log", "df"])
          crit <- if (is.finite(dfree) && dfree > 0) qt(.975, dfree) else 1.96
          start_est <- vapply(fitted$fits,
              function(x) unname(lme4::fixef(x)["ratio_log"]), numeric(1))
          ordered_ll <- sort(fitted$logLik, decreasing = TRUE)
          warns <- fitted$warnings
          singular <- lme4::isSingular(fit, tol = 1e-4)
          if (singular) warns <- c(warns, "selected fit is singular")
          out[[length(out) + 1L]] <- data.frame(
              covariates = set_name, random_effects = structure_name,
              estimate = estimate, std_error = se,
              ci_lower = estimate - crit * se, ci_upper = estimate + crit * se,
              p_value = unname(tt["ratio_log", "Pr(>|t|)"]),
              estimate_min_across_starts = min(start_est),
              estimate_max_across_starts = max(start_est),
              selected_start = fitted$selected_start,
              logLik = as.numeric(logLik(fit)),
              logLik_gap_to_next_start = if (length(ordered_ll) > 1L)
                  ordered_ll[1] - ordered_ll[2] else NA_real_,
              AIC = AIC(fit), singular_fit = singular,
              sample = if (set_name == "SDI + region (full eligible sample)") "Full eligible SDI + region sample" else "Common SES comparison sample",
              n_country_years = nrow(data_used),
              n_countries = length(unique(data_used$location_id)),
              warning = paste(unique(warns), collapse = " | "),
              stringsAsFactors = FALSE
          )
      }
  }
  results_structure <- do.call(rbind, out)
  utils::write.csv(results_structure,
      file.path(results_dir, "iqvia_model_structure_covariate_sensitivity.csv"),
      row.names = FALSE, na = "")
  message("No-prevalence model-structure analysis complete: preferred SDI + region fits use ",
          nrow(d_primary), " country-years from ", length(unique(d_primary$location_id)),
          " countries; alternative SES fits use ", nrow(d), " country-years from ",
          length(unique(d$location_id)), " countries; all fits retain incident-case weights.")
})


# -----------------------------------------------------------------------------
# Publication-window procurement models and Figure S16
# -----------------------------------------------------------------------------
local({
  # Refit the publication-window descriptive group and procurement-CS models.

  library(lme4)
  library(lmerTest)

  results_dir <- file.path("results", "iqvia")
  panel <- read.csv(file.path(results_dir, "iqvia_log_ratio_panel.csv"),
                    stringsAsFactors = FALSE)
  panel$location_id <- as.character(panel$location_id)
  panel$year <- as.integer(panel$year)
  panel$group_long <- factor(panel$group_long,
    levels = c("Shortage reported", "No shortage reported", "Unknown shortage status"))
  panel <- panel[panel$year >= 2014 & panel$year <= 2020, ]

  # A. Mean-based weighted description by country shortage status.
  summary_rows <- list()
  for (group_name in levels(panel$group_long)) {
    group_data <- panel[panel$group_long == group_name, ]
    if (!nrow(group_data)) next
    summary_rows[[length(summary_rows) + 1L]] <- data.frame(
      group = group_name,
      countries = length(unique(group_data$location_id)),
      country_years = nrow(group_data),
      mean_log1p_ratio = weighted.mean(group_data$log_ratio,
                                       group_data$syphilis_cases_guidrecom),
      back_transformed_weighted_mean = exp(weighted.mean(group_data$log_ratio,
        group_data$syphilis_cases_guidrecom)) - 1,
      pooled_dose_equivalents_per_1000_cases = 1000 *
        sum(group_data$annual_volume_dose_equivalents) /
        sum(group_data$syphilis_cases_guidrecom),
      stringsAsFactors = FALSE
    )
  }
  description <- do.call(rbind, summary_rows)
  write.csv(description, file.path(results_dir, "iqvia_2014_2020_group_description.csv"),
            row.names = FALSE)

  # B. Group comparison: shortage-reported is the reference, with categorical year.
  panel$year_factor <- factor(panel$year)
  panel$prior_weight <- panel$syphilis_cases_guidrecom /
    mean(panel$syphilis_cases_guidrecom)
  panel$location_id <- factor(panel$location_id)
  group_fit <- lmerTest::lmer(log_ratio ~ group_long + year_factor +
    (1 | location_id), data = panel, weights = prior_weight,
    REML = FALSE, na.action = na.fail)
  group_table <- as.data.frame(summary(group_fit)$coefficients)
  group_terms <- grep("^group_long", rownames(group_table), value = TRUE)
  group_results <- do.call(rbind, lapply(group_terms, function(term) {
    estimate <- group_table[term, "Estimate"]
    standard_error <- group_table[term, "Std. Error"]
    degrees_freedom <- group_table[term, "df"]
    critical_value <- qt(0.975, degrees_freedom)
    data.frame(
      term = term, estimate_log1p = estimate, standard_error = standard_error,
      df = degrees_freedom, ci_lower = estimate - critical_value * standard_error,
      ci_upper = estimate + critical_value * standard_error,
      p_value = group_table[term, "Pr(>|t|)"],
      ratio_of_geometric_1_plus_R = exp(estimate),
      countries = length(unique(panel$location_id)), country_years = nrow(panel),
      singular_fit = lme4::isSingular(group_fit), stringsAsFactors = FALSE
    )
  }))
  write.csv(group_results, file.path(results_dir, "iqvia_2014_2020_group_lmm.csv"),
            row.names = FALSE)
  fixed_terms <- names(lme4::fixef(group_fit))
  contrast_matrix <- matrix(0, nrow = 3, ncol = length(fixed_terms),
    dimnames = list(c("No shortage vs shortage", "Unknown vs shortage",
                      "Unknown vs no shortage"), fixed_terms))
  contrast_matrix[1, "group_longNo shortage reported"] <- 1
  contrast_matrix[2, "group_longUnknown shortage status"] <- 1
  contrast_matrix[3, "group_longUnknown shortage status"] <- 1
  contrast_matrix[3, "group_longNo shortage reported"] <- -1
  contrast_rows <- list()
  for (i in seq_len(nrow(contrast_matrix))) {
    contrast <- lmerTest::contest1D(group_fit, contrast_matrix[i, ],
                                    ddf = "Satterthwaite", confint = TRUE)
    contrast_rows[[i]] <- data.frame(
      contrast = rownames(contrast_matrix)[i], estimate = contrast$Estimate,
      standard_error = contrast[["Std. Error"]], df = contrast$df,
      ci_lower = contrast$lower, ci_upper = contrast$upper,
      p_value = contrast[["Pr(>|t|)" ]], countries = length(unique(panel$location_id)),
      country_years = nrow(panel), period = "2014-2020", stringsAsFactors = FALSE)
  }
  write.csv(do.call(rbind, contrast_rows),
            file.path(results_dir, "iqvia_2014_2020_group_all_contrasts.csv"),
            row.names = FALSE)

  # C. Association of log procurement intensity with congenital syphilis incidence.
  # Syphilis prevalence and incidence-rate covariates are omitted because the
  # procurement ratio is already defined per incident syphilis case.
  covariates <- as.data.frame(cs_panel)[, c("location_id", "year", "sdi", "gbd_superregion")]
  covariates$location_id <- as.character(covariates$location_id)
  association_data <- merge(panel, covariates, by = c("location_id", "year"),
                            all = FALSE, sort = FALSE)
  association_data$location_id <- factor(association_data$location_id)
  association_data$gbd_superregion <- factor(association_data$gbd_superregion)
  association_data$sdi_z <- as.numeric(scale(association_data$sdi))
  association_data$year_c <- association_data$year - 2016
  association_data$prior_weight <- association_data$syphilis_cases_guidrecom /
    mean(association_data$syphilis_cases_guidrecom)
  association_formula <- cs_rate ~ log_ratio + factor(year) + sdi_z +
    gbd_superregion + (1 + year_c | location_id)
  association_fit <- lmerTest::lmer(association_formula, data = association_data,
    weights = prior_weight, REML = FALSE,
    control = lme4::lmerControl(optimizer = "bobyqa"), na.action = na.fail)
  association_table <- as.data.frame(summary(association_fit)$coefficients)
  term <- "log_ratio"
  estimate <- association_table[term, "Estimate"]
  standard_error <- association_table[term, "Std. Error"]
  degrees_freedom <- association_table[term, "df"]
  critical_value <- qt(0.975, degrees_freedom)
  association_result <- data.frame(
    outcome = "GBD congenital syphilis incidence per 1,000 live births",
    exposure = "log(1 + BPG dose-equivalents per 1,000 guideline-recommended incident cases)",
    estimate = estimate, standard_error = standard_error, df = degrees_freedom,
    ci_lower = estimate - critical_value * standard_error,
    ci_upper = estimate + critical_value * standard_error,
    p_value = association_table[term, "Pr(>|t|)"],
    country_years = nrow(association_data),
    countries = length(unique(association_data$location_id)),
    singular_fit = lme4::isSingular(association_fit),
    formula = paste(deparse(association_formula), collapse = " "),
    covariates = "Categorical year; SDI; GBD super-region; no syphilis prevalence",
    weights = "GBD incident syphilis cases in women aged 15-49 plus children aged 0-14; rescaled to mean 1",
    stringsAsFactors = FALSE
  )
  write.csv(association_result,
            file.path(results_dir, "iqvia_2014_2020_cs_association.csv"),
            row.names = FALSE)

  cat("2014-2020 IQVIA subgroup description, group LMM and CS association completed.\n")


  # -----------------------------------------------------------------------------
  # Procurement publication figures
  # -----------------------------------------------------------------------------

  # Publication figures
  # Recreate procurement coverage and association figures from saved v7 results.

  library(ggplot2)
  library(patchwork)

  results_dir <- file.path("results", "iqvia")
  figure_dir <- file.path("figures", "iqvia")
  publication_figure_dir <- "figures"
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  group_colours <- c("shortage reported" = "#C63D40",
                     "no shortage reported" = "#176B9A",
                     "unclassified" = "#737A83")

  # Supplementary Figure S16. Ratio uses the guideline-recommended incident-case
  # denominator; the lower panel documents how many countries contribute each year.
  annual <- read.csv(file.path(results_dir,
    "iqvia_incident_ratio_group_year_description.csv"), stringsAsFactors = FALSE)
  annual <- annual[annual$year %in% 2014:2020, ]
  write.csv(annual, file.path(results_dir, "figure_s16_annual_plot_data.csv"),
            row.names = FALSE)
  annual$group <- factor(annual$group,
    levels = c("shortage reported", "no shortage reported", "unclassified"))
  ratio_plot <- ggplot(annual, aes(x = year, y = aggregate_ratio_per_1000_guidrecom,
    colour = group, shape = group, linetype = group)) +
    geom_line(linewidth = 0.8) + geom_point(size = 2) +
    scale_colour_manual(values = group_colours) +
    scale_x_continuous(breaks = 2014:2020) +
    scale_y_log10(labels = scales::label_number(big.mark = ",")) +
    labs(title = "A  Aggregate BPG procurement intensity", x = NULL,
         y = "2.4-MIU dose-equivalents per 1,000\nincident syphilis cases",
         colour = NULL, shape = NULL, linetype = NULL) +
    publication_theme(base_size = 10) + theme(legend.position = "none")
  coverage_plot <- ggplot(annual, aes(x = year, y = n_countries,
    colour = group, shape = group, linetype = group)) +
    geom_line(linewidth = 0.7) + geom_point(size = 2) +
    scale_colour_manual(values = group_colours) +
    scale_x_continuous(breaks = 2014:2020) +
    scale_y_continuous(limits = c(0, NA), breaks = 0:15) +
    labs(title = "B  Countries contributing each year", x = "Year",
         y = "Countries", colour = NULL, shape = NULL, linetype = NULL) +
    publication_theme(base_size = 10) + theme(legend.position = "top")
  figure_s16 <- ratio_plot / coverage_plot + plot_layout(heights = c(1.7, 1))
  ggsave(file.path(publication_figure_dir, "Figure_S16_Procurement_ratio_and_coverage.png"),
         figure_s16, width = 8.5, height = 6.5, dpi = 300)
  ggsave(file.path(publication_figure_dir, "Figure_S16_Procurement_ratio_and_coverage.pdf"), device = grDevices::cairo_pdf,
         figure_s16, width = 8.5, height = 6.5)

  # The additional procurement forest plot was removed as redundant.
  cat("Saved Supplementary Figure S16. Procurement sensitivity results are tabulated.\n")
})

} else {
  message('IQVIA section skipped because licensed source files or required GBD population/incidence extracts are unavailable.')
}
hiv_art_filename <- paste0(
  "HIV ", intToUtf8(0x2013),
  " Estimated antiretroviral therapy coverage among people living with HIV (%).csv"
)
hiv_source_files <- c(
  file.path(raw_data_dir, "AIDSInfo",
    "Elimination of vertical transmission_Coverage of pregnant women who receive ARV for PMTCT.csv"),
  file.path(raw_data_dir, "AIDSInfo",
    "Elimination of vertical transmission_Pregnant women needing ARV for PMTCT.csv"),
  file.path(raw_data_dir, "UN_WPP", "WPP2024_GEN_F01_DEMOGRAPHIC_INDICATORS_COMPACT.xlsx"),
  file.path(raw_data_dir, "WHOData", hiv_art_filename)
)


# -----------------------------------------------------------------------------
# HIV service-covariate source join and audit
# -----------------------------------------------------------------------------
if (all(file.exists(hiv_source_files))) {


# -----------------------------------------------------------------------------
# Prepare annual HIV service-covariate audit
# -----------------------------------------------------------------------------
local({
  # Attach annual PMTCT and WHO ART coverage to the GBD Birth-age HIV panel.
  # Optional source join and availability audit; this does not create the complete input set required by the HIV model section of 04_additional_analyses.R.

  library(readr)

  d <- hiv_panel[hiv_panel$year %in% study_years, ]
  d$location_iso <- toupper(trimws(as.character(d$location_iso)))
  d$location_name <- trimws(as.character(d$location_name))
  d$year <- as.integer(d$year)

  # Parse AIDSInfo intervals consistently: less-than values use half the threshold;
  # greater-than coverage uses the midpoint from the threshold to 100%.
  parse_percentage <- function(x) {
    x <- gsub(",", "", trimws(as.character(x)), fixed = TRUE)
    value <- suppressWarnings(as.numeric(x))
    threshold <- suppressWarnings(as.numeric(sub("^[<>]\\s*", "", x)))
    value[grepl("^<\\s*", x) & is.finite(threshold)] <-
      threshold[grepl("^<\\s*", x) & is.finite(threshold)] / 2
    value[grepl("^>\\s*", x) & is.finite(threshold)] <-
      threshold[grepl("^>\\s*", x) & is.finite(threshold)] +
      (100 - threshold[grepl("^>\\s*", x) & is.finite(threshold)]) / 2
    value[x %in% c("", "...", "..", "NA", "N/A")] <- NA_real_
    value
  }

  # Read UNAIDS reported PMTCT coverage and reshape the annual columns.
  pmtct_file <- file.path(raw_data_dir, "AIDSInfo",
    "Elimination of vertical transmission_Coverage of pregnant women who receive ARV for PMTCT.csv")
  pmtct_raw <- utils::read.csv(pmtct_file, check.names = FALSE,
                               stringsAsFactors = FALSE)
  country_lookup <- unique(d[, c("location_name", "location_iso")])
  country_lookup <- country_lookup[!duplicated(country_lookup$location_name), ]
  pmtct_rows <- list()
  for (year in study_years) {
    year_column <- as.character(year)
    pmtct_rows[[as.character(year)]] <- data.frame(
      source_country = trimws(as.character(pmtct_raw$Country)),
      year = year,
      pmtct_percent = parse_percentage(pmtct_raw[[year_column]]),
      stringsAsFactors = FALSE
    )
  }
  pmtct <- do.call(rbind, pmtct_rows)
  pmtct <- merge(pmtct, country_lookup, by.x = "source_country",
                 by.y = "location_name", all.x = TRUE)
  unmatched_pmtct <- unique(pmtct[is.na(pmtct$location_iso), "source_country"])
  pmtct <- pmtct[!is.na(pmtct$location_iso),
                 c("location_iso", "year", "pmtct_percent")]
  if (anyDuplicated(pmtct[c("location_iso", "year")])) {
    stop("Duplicate PMTCT country-year values after name matching.")
  }

  # Estimate maternal HIV burden from pregnant women needing PMTCT ARVs per
  # 100 UN WPP live births. This is the pregnancy-specific maternal HIV proxy.
  parse_count <- function(x) {
    x <- gsub(",", "", trimws(as.character(x)), fixed = TRUE)
    value <- suppressWarnings(as.numeric(x))
    threshold <- suppressWarnings(as.numeric(sub("^[<>]\\s*", "", x)))
    less_than <- grepl("^<\\s*", x) & is.finite(threshold)
    greater_than <- grepl("^>\\s*", x) & is.finite(threshold)
    value[less_than] <- threshold[less_than] / 2
    value[greater_than] <- threshold[greater_than] * 1.05
    value[x %in% c("", "...", "..", "NA", "N/A")] <- NA_real_
    value
  }
  pregnant_count_file <- file.path(raw_data_dir, "AIDSInfo",
    "Elimination of vertical transmission_Pregnant women needing ARV for PMTCT.csv")
  pregnant_raw <- utils::read.csv(pregnant_count_file, check.names = FALSE,
                                  stringsAsFactors = FALSE)
  pregnant_rows <- list()
  for (year in study_years) {
    pregnant_rows[[as.character(year)]] <- data.frame(
      source_country = trimws(as.character(pregnant_raw$Country)),
      year = year,
      pregnant_hiv_count = parse_count(pregnant_raw[[as.character(year)]]),
      stringsAsFactors = FALSE
    )
  }
  pregnant_hiv <- do.call(rbind, pregnant_rows)
  pregnant_hiv <- merge(pregnant_hiv, country_lookup, by.x = "source_country",
                        by.y = "location_name", all.x = TRUE)
  unmatched_pregnant_hiv <- unique(
    pregnant_hiv[is.na(pregnant_hiv$location_iso), "source_country"]
  )
  pregnant_hiv <- pregnant_hiv[!is.na(pregnant_hiv$location_iso), ]

  wpp_file <- file.path(raw_data_dir, "UN_WPP",
    "WPP2024_GEN_F01_DEMOGRAPHIC_INDICATORS_COMPACT.xlsx")
  wpp_raw <- openxlsx::read.xlsx(wpp_file, sheet = "Estimates", colNames = FALSE,
                                 skipEmptyRows = FALSE, skipEmptyCols = FALSE)
  header_row <- which(trimws(as.character(wpp_raw[[1]])) == "Index")
  if (length(header_row) != 1L) stop("Could not identify the UN WPP header row.")
  wpp_headers <- trimws(as.character(unlist(wpp_raw[header_row, ], use.names = FALSE)))
  wpp_needed <- c("Variant", "Region, subregion, country or area *",
                  "ISO3 Alpha-code", "Year", "Births (thousands)")
  wpp_columns <- match(wpp_needed, wpp_headers)
  if (anyNA(wpp_columns)) stop("UN WPP workbook lacks a required live-birth field.")
  wpp <- as.data.frame(wpp_raw[seq.int(header_row + 1L, nrow(wpp_raw)),
                                wpp_columns, drop = FALSE])
  names(wpp) <- c("variant", "country_name", "location_iso", "year", "births_thousands")
  wpp$year <- suppressWarnings(as.integer(as.character(wpp$year)))
  wpp$live_births <- suppressWarnings(as.numeric(as.character(wpp$births_thousands))) * 1000
  wpp <- wpp[as.character(wpp$variant) == "Estimates" &
               wpp$year %in% study_years & is.finite(wpp$live_births) &
               wpp$live_births > 0 & !is.na(wpp$location_iso), ]
  wpp$location_iso <- toupper(trimws(as.character(wpp$location_iso)))
  wpp <- wpp[, c("location_iso", "year", "live_births")]
  if (anyDuplicated(wpp[c("location_iso", "year")])) {
    stop("Duplicate UN WPP live-birth country-years.")
  }
  pregnant_hiv <- merge(pregnant_hiv, wpp, by = c("location_iso", "year"),
                        all = FALSE)
  pregnant_hiv$maternal_hiv_per_100_births <-
    100 * pregnant_hiv$pregnant_hiv_count / pregnant_hiv$live_births
  pregnant_hiv <- pregnant_hiv[, c("location_iso", "year",
                                   "maternal_hiv_per_100_births")]
  if (anyDuplicated(pregnant_hiv[c("location_iso", "year")])) {
    stop("Duplicate maternal HIV country-years after WPP matching.")
  }

  # WHO's GHO ART-coverage file uses an en dash in its filename. Construct it by
  # Unicode code point so the source remains readable across Windows locales.
  art_filename <- paste0(
    "HIV ", intToUtf8(0x2013),
    " Estimated antiretroviral therapy coverage among people living with HIV (%).csv"
  )
  art_file <- file.path(raw_data_dir, "WHOData", art_filename)
  if (!file.exists(art_file)) stop("WHO ART-coverage source file was not found.")
  art_raw <- readr::read_csv(art_file, show_col_types = FALSE,
                             locale = readr::locale(encoding = "UTF-8"),
                             name_repair = "minimal")
  art_raw <- art_raw[art_raw$IndicatorCode == "HIV_ARTCOVERAGE" &
                       art_raw$`Location type` == "Country" &
                       as.integer(art_raw$Period) %in% study_years, ]
  art <- data.frame(
    location_iso = toupper(trimws(as.character(art_raw$SpatialDimValueCode))),
    year = as.integer(art_raw$Period),
    art_percent = suppressWarnings(as.numeric(art_raw$FactValueNumeric)),
    modified = suppressWarnings(as.Date(substr(art_raw$DateModified, 1, 10))),
    stringsAsFactors = FALSE
  )
  art$modified[is.na(art$modified)] <- as.Date("1900-01-01")
  art <- art[order(art$location_iso, art$year, -as.numeric(art$modified)), ]
  art <- art[!duplicated(art[c("location_iso", "year")]), ]
  art$modified <- NULL

  # Keep source values and missingness visible in the saved audit.
  d <- merge(d, pmtct, by = c("location_iso", "year"), all.x = TRUE)
  d <- merge(d, art, by = c("location_iso", "year"), all.x = TRUE)
  d <- merge(d, pregnant_hiv, by = c("location_iso", "year"), all.x = TRUE)
  d$pmtct_z <- as.numeric(scale(d$pmtct_percent))
  d$art_z <- as.numeric(scale(d$art_percent))
  d$anc4_z <- as.numeric(scale(d$pregnant_antenatal_care4))
  d$maternal_hiv_z <- as.numeric(scale(d$maternal_hiv_per_100_births))
  write.csv(d, file.path(results_dir, "hiv_birth_annual_panel.csv"),
            row.names = FALSE, na = "")
  saveRDS(d, file.path(results_dir, "hiv_birth_annual_panel.rds"))

  variables <- c("sdi", "maternal_hiv_per_100_births", "hiv_prev_female",
                 "art_percent", "pmtct_percent",
                 "pregnant_antenatal_care4")
  audit <- do.call(rbind, lapply(variables, function(variable) data.frame(
    covariate = variable,
    available_country_years = sum(is.finite(d[[variable]])),
    missing_country_years = sum(!is.finite(d[[variable]])),
    available_countries = length(unique(d$location_id[is.finite(d[[variable]])])),
    first_observed_year = if (any(is.finite(d[[variable]])))
      min(d$year[is.finite(d[[variable]])]) else NA_integer_,
    last_observed_year = if (any(is.finite(d[[variable]])))
      max(d$year[is.finite(d[[variable]])]) else NA_integer_,
    stringsAsFactors = FALSE
  )))
  write.csv(audit, file.path(results_dir, "hiv_birth_covariate_availability.csv"),
            row.names = FALSE, na = "")
  write.csv(data.frame(unmatched_pmtct_country = unmatched_pmtct),
            file.path(logs_dir, "hiv_pmtct_unmatched_country_names.csv"),
            row.names = FALSE, na = "")
  write.csv(data.frame(unmatched_pregnant_HIV_country = unmatched_pregnant_hiv),
            file.path(logs_dir, "hiv_pregnant_hiv_unmatched_country_names.csv"),
            row.names = FALSE, na = "")

  cat("Saved the Birth-age HIV analysis panel and annual covariate availability audit.\n")
})

} else {
  message('HIV source join skipped because one or more UNAIDS, WHO, or UN WPP source files are missing.')
}
hiv_model_inputs <- file.path(processed_data_dir, "hiv", c(
  "hiv_birth_annual_panel.rds", "hiv_under5_annual_panel.rds", "covariate_scaling.csv"
))


# -----------------------------------------------------------------------------
# HIV negative-control models and Figure S15
# -----------------------------------------------------------------------------
if (all(file.exists(hiv_model_inputs))) {


# -----------------------------------------------------------------------------
# HIV model support, six model branches, and publication figure
# -----------------------------------------------------------------------------
local({
  # Shared model-support routines, followed by the six HIV model branches and Figure S15.
  # Shared setup and model helpers for the HIV negative-control analyses.

  # Shared paths for the HIV-specific adjusted models in this analysis.
  # Data panels retain the audited GBD and UNAIDS country-year joins.
  results_dir <- file.path("results", "hiv_specific")
  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)


  # -----------------------------------------------------------------------------
  # Shared annual-model helper routines
  # -----------------------------------------------------------------------------

  # Shared model definitions for the six HIV annual analyses.
  # Run analysis scripts from the analysis root.
  suppressPackageStartupMessages(library(nlme))
  BREAK_YEAR <- 2016L
  N_IMPUTATIONS <- 50L
  MIN_PMM_CONVERGED <- 45L

  fit_lme <- function(d, covariates, p = 0L, q = 0L,
                      model_label = "", break_year = BREAK_YEAR) {
      d <- d[order(d$location_id, d$year), , drop = FALSE]
      d$location_id <- droplevels(factor(d$location_id))
      d$gbd_superregion <- droplevels(factor(d$gbd_superregion))
      d$year <- as.integer(d$year)
      d$time_c <- d$year - as.integer(break_year)
      # Scale only the random time slope to improve covariance optimization;
      # fixed CITS coefficients remain on the original per-year scale.
      d$time_5 <- d$time_c / 5
      d$post <- as.integer(d$year >= break_year)
      d$population_relative <- d$population / mean(d$population)
      fixed_formula <- stats::as.formula(paste0(
          "rate ~ time_c * shortage * post + ", paste(covariates, collapse = " + ")
      ))
      if (p != 0L || q != 0L) stop("Country-level HIV LMMs do not fit ARMA.")
      warning_text <- character()
      fit <- withCallingHandlers(
          nlme::lme(
              fixed = fixed_formula,
              random = ~ 1 + time_5 | location_id,
              data = d,
              correlation = NULL,
              weights = nlme::varFixed(~ I(1 / population_relative)),
              method = "ML", na.action = na.fail,
              control = nlme::lmeControl(
                  opt = "nlminb", maxIter = 1000,
                  msMaxIter = 1000, msMaxEval = 2000,
                  niterEM = 50, returnObject = TRUE
              )
          ),
          warning = function(w) {
              warning_text <<- c(warning_text, conditionMessage(w))
              invokeRestart("muffleWarning")
          }
      )
      attr(fit, "analysis_data") <- d
      attr(fit, "fit_warnings") <- paste(unique(warning_text), collapse = " | ")
      attr(fit, "model_label") <- model_label
      attr(fit, "p") <- p
      attr(fit, "q") <- q
      fit
  }

  is_converged <- function(fit) {
      if (!inherits(fit, "lme")) return(FALSE)
      !grepl(
          "convergence error|iteration limit reached|singular convergence|false convergence",
          attr(fit, "fit_warnings"), ignore.case = TRUE
      )
  }

  get_aicc <- function(fit) {
      ll <- as.numeric(stats::logLik(fit))
      k <- attr(stats::logLik(fit), "df")
      n <- stats::nobs(fit)
      stats::AIC(fit) + 2 * k * (k + 1) / (n - k - 1)
  }

  extract_terms <- function(fit, branch, model, imputation = NA_integer_) {
      tab <- as.data.frame(summary(fit)$tTable)
      tab$term <- rownames(tab)
      dd <- attr(fit, "analysis_data")
      n_post_years <- length(unique(dd$year[dd$post == 1L]))
      terms <- c("time_c:shortage", "shortage:post", "time_c:shortage:post")
      rows <- lapply(terms, function(term) {
          if (term == "time_c:shortage:post" && n_post_years < 2L) return(NULL)
          z <- tab[tab$term == term, , drop = FALSE]
          if (!nrow(z)) return(NULL)
          df <- z$DF[1]
          crit <- stats::qt(0.975, df = df)
          data.frame(
              branch = branch, model = model, term = term,
              estimate = z$Value[1], se = z$Std.Error[1], df = df,
              ci_lower = z$Value[1] - crit * z$Std.Error[1],
              ci_upper = z$Value[1] + crit * z$Std.Error[1],
              p_value = z$`p-value`[1],
              n_country_years = stats::nobs(fit),
              n_countries = length(unique(dd$location_id)),
              n_pre_years = length(unique(dd$year[dd$post == 0L])),
              n_post_years = n_post_years,
              imputation = imputation,
              stringsAsFactors = FALSE
          )
      })
      do.call(rbind, Filter(Negate(is.null), rows))
  }


  # One ordinary weighted LMM per outcome and missing-data strategy.
  fit_lmm <- function(d, covariates, branch) {
      fit <- tryCatch(fit_lme(d, covariates, model_label = "LMM"), error = function(e) e)
      ok <- is_converged(fit)
      attempts <- data.frame(
          branch = branch, model = "LMM", converged = ok,
          error = if (ok) "" else if (inherits(fit, "lme"))
              "Optimizer convergence warning" else conditionMessage(fit),
          warnings = if (inherits(fit, "lme")) attr(fit, "fit_warnings") else "")
      if (!ok) return(list(fits = list(), attempts = attempts,
                           models = data.frame(), selected = NA_character_))
      information <- data.frame(
          branch = branch, model = "LMM", n_country_years = stats::nobs(fit),
          n_countries = length(unique(d$location_id)), AICc = get_aicc(fit))
      list(fits = list(LMM = fit), attempts = attempts,
           models = information, selected = "LMM")
  }

  # Rubin pooling for the three CITS terms from the same repeated fit set.
  pool_coefficients <- function(est, se, term, branch, model, n_countries,
                                n_rows, successful, requested) {
      ok <- is.finite(est) & is.finite(se) & se >= 0
      est <- est[ok]; se <- se[ok]
      m <- length(est)
      if (m < MIN_PMM_CONVERGED) {
          return(data.frame(
              branch = branch, model = model, term = term,
              estimate = NA_real_, se = NA_real_, df = NA_real_,
              ci_lower = NA_real_, ci_upper = NA_real_, p_value = NA_real_,
              n_imputations = m, requested_imputations = requested,
              n_country_years = n_rows, n_countries = n_countries,
              stringsAsFactors = FALSE
          ))
      }
      qbar <- mean(est); ubar <- mean(se^2); between <- stats::var(est)
      total <- ubar + (1 + 1 / m) * between
      pooled_se <- sqrt(total)
      if (!is.finite(between) || between < .Machine$double.eps) {
          df <- Inf
          crit <- stats::qnorm(0.975)
          p <- 2 * stats::pnorm(abs(qbar / pooled_se), lower.tail = FALSE)
      } else {
          r <- (1 + 1 / m) * between / ubar
          df <- (m - 1) * (1 + 1 / r)^2
          crit <- stats::qt(0.975, df = df)
          p <- 2 * stats::pt(abs(qbar / pooled_se), df = df, lower.tail = FALSE)
      }
      data.frame(
          branch = branch, model = model, term = term,
          estimate = qbar, se = pooled_se, df = df,
          ci_lower = qbar - crit * pooled_se,
          ci_upper = qbar + crit * pooled_se, p_value = p,
          n_imputations = m, requested_imputations = requested,
          n_country_years = n_rows, n_countries = n_countries,
          stringsAsFactors = FALSE
      )
  }

  # Country-stratified PMM for annual AIDSInfo treatment series. Within each
  # country, a piecewise time model includes the 2016 level and slope break;
  # donors are that country's observed years, matched on predicted value and
  # HIV outcome/population profile.
  country_pmm <- function(d, variables, m = N_IMPUTATIONS, k = 3L) {
      n <- nrow(d)
      imputed_draws <- setNames(vector("list", length(variables)), variables)
      missing_indices <- setNames(vector("list", length(variables)), variables)
      for (v in variables) {
          missing_indices[[v]] <- which(!is.finite(d[[v]]))
          imputed_draws[[v]] <- matrix(
              NA_real_, nrow = length(missing_indices[[v]]), ncol = m
          )
      }

      country_ids <- unique(as.character(d$location_id))
      row_to_missing <- lapply(variables, function(v) {
          stats::setNames(seq_along(missing_indices[[v]]),
                          as.character(missing_indices[[v]]))
      })
      names(row_to_missing) <- variables

      for (id in country_ids) {
          rows <- which(as.character(d$location_id) == id)
          country <- d[rows, , drop = FALSE]
          trend <- list()
          for (v in variables) {
              observed <- which(is.finite(country[[v]]))
              if (!length(observed)) next

              trend_data <- data.frame(
                  value = country[[v]][observed],
                  time_c = country$time_c[observed],
                  post = country$post[observed]
              )
              fit <- NULL
              if (length(observed) >= 4L && length(unique(country$year[observed])) >= 4L) {
                  piecewise_fit <- stats::lm(
                      value ~ time_c + post + time_c:post, data = trend_data
                  )
                  piecewise_design <- stats::model.matrix(piecewise_fit)
                  if (qr(piecewise_design)$rank == ncol(piecewise_design)) {
                      fit <- piecewise_fit
                  }
              }
              if (is.null(fit) && length(unique(trend_data$time_c)) >= 2L) {
                  # A simple linear trend keeps donor matching estimable when
                  # the country has too few pre/post observations for a break.
                  fit <- stats::lm(value ~ time_c, data = trend_data)
              }

              if (!is.null(fit)) {
                  new_all <- data.frame(time_c = country$time_c, post = country$post)
                  pred_all <- as.numeric(stats::predict(fit, newdata = new_all))
                  pred_obs <- pred_all[observed]
                  residuals <- stats::residuals(fit)
                  if (stats::df.residual(fit) > 0L) {
                      pred_se <- as.numeric(stats::predict(
                          fit, newdata = new_all, se.fit = TRUE
                      )$se.fit)
                  } else {
                      pred_se <- rep(0, nrow(country))
                  }
                  if (!length(residuals)) residuals <- 0
                  pred_all[!is.finite(pred_all)] <- mean(country[[v]][observed])
                  pred_obs <- pred_all[observed]
                  pred_se[!is.finite(pred_se)] <- 0
              } else {
                  pred_all <- rep(country[[v]][observed[1]], nrow(country))
                  pred_obs <- rep(country[[v]][observed[1]], length(observed))
                  residuals <- 0
                  pred_se <- rep(0, nrow(country))
              }
              trend[[v]] <- list(
                  observed = observed,
                  predicted_all = pred_all,
                  predicted_observed = pred_obs,
                  residuals = residuals,
                  prediction_se = pred_se
              )
          }

          if (!length(trend)) {
              stop("No observed treatment donor for country ", id)
          }
          for (m_id in seq_len(m)) {
              for (local_row in seq_len(nrow(country))) {
                  missing_vars <- variables[vapply(variables, function(v) {
                      !is.finite(country[[v]][local_row])
                  }, logical(1))]
                  if (!length(missing_vars)) next
                  donor_local <- which(Reduce(`&`, lapply(missing_vars, function(v) {
                      is.finite(country[[v]])
                  })))
                  if (!length(donor_local)) {
                      stop("No joint PMM donor values for country ", id)
                  }

                  distance2 <- rep(0, length(donor_local))
                  components <- 0L
                  for (v in missing_vars) {
                      tr <- trend[[v]]
                      donor_match <- match(donor_local, tr$observed)
                      if (anyNA(donor_match)) {
                          stop("PMM donor does not have an observed value for ", v)
                      }
                      uncertainty <- sqrt(
                          tr$prediction_se[local_row]^2 +
                              stats::sd(tr$residuals, na.rm = TRUE)^2
                      )
                      if (!is.finite(uncertainty)) uncertainty <- 0
                      target_prediction <- tr$predicted_all[local_row] +
                          stats::rnorm(1, mean = 0, sd = uncertainty)
                      scale <- stats::sd(tr$predicted_observed, na.rm = TRUE)
                      if (!is.finite(scale) || scale < 1e-8) {
                          scale <- stats::sd(country[[v]][tr$observed], na.rm = TRUE)
                      }
                      if (!is.finite(scale) || scale < 1e-8) scale <- 1
                      distance2 <- distance2 +
                          ((tr$predicted_observed[donor_match] - target_prediction) / scale)^2
                      components <- components + 1L
                  }

                  # Include the observed HIV outcome and denominator profile as
                  # auxiliary PMM predictors without adding them to the outcome
                  # regression adjustment set.
                  for (aux in c("rate", "log_population")) {
                      aux_scale <- stats::sd(country[[aux]], na.rm = TRUE)
                      if (!is.finite(aux_scale) || aux_scale < 1e-8) next
                      distance2 <- distance2 +
                          ((country[[aux]][donor_local] - country[[aux]][local_row]) /
                               aux_scale)^2
                      components <- components + 1L
                  }
                  distance <- sqrt(distance2 / max(components, 1L))
                  nearest <- order(distance)[seq_len(min(k, length(distance)))]
                  chosen <- sample(donor_local[nearest], size = 1L)
                  global_row <- rows[local_row]
                  global_donor <- rows[chosen]
                  for (v in missing_vars) {
                      draw_row <- row_to_missing[[v]][as.character(global_row)]
                      imputed_draws[[v]][draw_row, m_id] <- d[[v]][global_donor]
                  }
              }
          }
      }
      for (v in variables) {
          if (length(missing_indices[[v]]) &&
              any(!is.finite(imputed_draws[[v]]))) {
              stop("Country PMM left missing draws for ", v)
          }
      }
      list(draws = imputed_draws, indices = missing_indices,
           method = "Within-country PMM; k=3 nearest observed donor years")
  }


  # Birth-age and under-5 congenital HIV negative-control models.
  # Requires the prepared Birth-age and under-5 panels and covariate scaling file under results/hiv_specific/Data; see the package README.



  # -----------------------------------------------------------------------------
  # 02_birth_complete_country
  # -----------------------------------------------------------------------------

  # Birth-age HIV incidence: complete-country primary analysis.
  # Required prepared model inputs are listed in the package README.

  set.seed(20260926)

  panel <- readRDS(file.path(processed_data_dir, "hiv", "hiv_birth_annual_panel.rds"))
  result_dir <- file.path(results_dir, "Birth", "Complete_country")
  dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

  analysis_name <- "Birth | Maternal HIV burden + PMTCT | Complete annual countries"
  covariates <- c(
      "hdi_z", "gni_z", "sdi_z", "gbd_superregion",
      "hiv_pregnant_per_100_births_z", "pmtct_z"
  )
  # Retain rows with observed outcomes, population, exposure group, and DAG covariates.
  base_ok <- is.finite(panel$rate) & is.finite(panel$population) & panel$population > 0 &
      !is.na(panel$shortage) & !is.na(panel$gbd_superregion) &
      is.finite(panel$hdi_z) & is.finite(panel$gni_z) & is.finite(panel$sdi_z)
  eligible <- panel[base_ok, , drop = FALSE]
  exposure_ok <- is.finite(eligible$hiv_pregnant_per_100_births) &
      is.finite(eligible$pmtct_pct)
  ids_with_data <- unique(eligible$location_id[exposure_ok])

  # The primary sample keeps only countries with all 12 annual covariate records.
  complete_ids <- ids_with_data[vapply(ids_with_data, function(id) {
      observed_years <- eligible$year[
          eligible$location_id == id & exposure_ok
      ]
      all(2010:2021 %in% observed_years)
  }, logical(1))]
  analysis_data <- eligible[
      eligible$location_id %in% complete_ids & exposure_ok & eligible$year %in% 2010:2021,
      , drop = FALSE
  ]

  analysis_data$location_id <- factor(analysis_data$location_id)
  analysis_data$gbd_superregion <- factor(analysis_data$gbd_superregion)
  branch <- analysis_name

  # Fit the prespecified ordinary weighted LMM for this sample.
  selection <- fit_lmm(analysis_data, covariates, branch)
  if (!length(selection$fits)) stop("No converged model for ", analysis_name)
  models_to_report <- "LMM"

  effect_rows <- list()
  for (model_name in models_to_report) {
      fit <- selection$fits[[model_name]]
      if (is.null(fit)) next
      effect_rows[[model_name]] <- extract_terms(fit, branch, model_name)
  }

  write.csv(do.call(rbind, effect_rows), file.path(result_dir, "effects.csv"), row.names = FALSE)
  write.csv(data.frame(
      analysis = analysis_name, n_country_years = nrow(analysis_data),
      n_countries = length(unique(analysis_data$location_id)),
      shortage_countries = length(unique(analysis_data$location_id[analysis_data$shortage == 1])),
      comparison_countries = length(unique(analysis_data$location_id[analysis_data$shortage == 0])),
      years = paste(sort(unique(analysis_data$year)), collapse = ";")
  ), file.path(result_dir, "sample_flow.csv"), row.names = FALSE)
  write.csv(selection$models, file.path(result_dir, "lmm_fit_information.csv"), row.names = FALSE)
  write.csv(selection$attempts, file.path(result_dir, "lmm_fit_attempts.csv"), row.names = FALSE)
  write.csv(data.frame(analysis = analysis_name, model = "LMM"),
            file.path(result_dir, "model_specification.csv"), row.names = FALSE)
  message("Finished: ", analysis_name)


  # -----------------------------------------------------------------------------
  # 03_birth_observed_years
  # -----------------------------------------------------------------------------

  # Birth-age HIV incidence: complete-country primary analysis.
  # Required prepared model inputs are listed in the package README.

  set.seed(20260926)

  panel <- readRDS(file.path(processed_data_dir, "hiv", "hiv_birth_annual_panel.rds"))
  result_dir <- file.path(results_dir, "Birth", "Observed_years")
  dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

  analysis_name <- "Birth | Maternal HIV burden + PMTCT | Observed country-years"
  covariates <- c(
      "hdi_z", "gni_z", "sdi_z", "gbd_superregion",
      "hiv_pregnant_per_100_births_z", "pmtct_z"
  )
  # Retain rows with observed outcomes, population, exposure group, and DAG covariates.
  base_ok <- is.finite(panel$rate) & is.finite(panel$population) & panel$population > 0 &
      !is.na(panel$shortage) & !is.na(panel$gbd_superregion) &
      is.finite(panel$hdi_z) & is.finite(panel$gni_z) & is.finite(panel$sdi_z)
  eligible <- panel[base_ok, , drop = FALSE]
  exposure_ok <- is.finite(eligible$hiv_pregnant_per_100_births) &
      is.finite(eligible$pmtct_pct)
  ids_with_data <- unique(eligible$location_id[exposure_ok])

  # Keep every country-year with observed treatment covariates.
  analysis_data <- eligible[exposure_ok & eligible$year %in% 2010:2021, , drop = FALSE]
  analysis_data$location_id <- factor(analysis_data$location_id)
  analysis_data$gbd_superregion <- factor(analysis_data$gbd_superregion)
  branch <- analysis_name

  # Fit the prespecified ordinary weighted LMM for this sample.
  selection <- fit_lmm(analysis_data, covariates, branch)
  if (!length(selection$fits)) stop("No converged model for ", analysis_name)
  models_to_report <- "LMM"

  effect_rows <- list()
  for (model_name in models_to_report) {
      fit <- selection$fits[[model_name]]
      if (is.null(fit)) next
      effect_rows[[model_name]] <- extract_terms(fit, branch, model_name)
  }

  write.csv(do.call(rbind, effect_rows), file.path(result_dir, "effects.csv"), row.names = FALSE)
  write.csv(data.frame(
      analysis = analysis_name, n_country_years = nrow(analysis_data),
      n_countries = length(unique(analysis_data$location_id)),
      shortage_countries = length(unique(analysis_data$location_id[analysis_data$shortage == 1])),
      comparison_countries = length(unique(analysis_data$location_id[analysis_data$shortage == 0])),
      years = paste(sort(unique(analysis_data$year)), collapse = ";")
  ), file.path(result_dir, "sample_flow.csv"), row.names = FALSE)
  write.csv(selection$models, file.path(result_dir, "lmm_fit_information.csv"), row.names = FALSE)
  write.csv(selection$attempts, file.path(result_dir, "lmm_fit_attempts.csv"), row.names = FALSE)
  write.csv(data.frame(analysis = analysis_name, model = "LMM"),
            file.path(result_dir, "model_specification.csv"), row.names = FALSE)
  message("Finished: ", analysis_name)



  # -----------------------------------------------------------------------------
  # 04_under5_complete_country
  # -----------------------------------------------------------------------------

  # Under-5 HIV incidence: complete-country primary analysis.
  # Required prepared model inputs are listed in the package README.

  set.seed(20260926)

  panel <- readRDS(file.path(processed_data_dir, "hiv", "hiv_under5_annual_panel.rds"))
  result_dir <- file.path(results_dir, "under-5", "Complete_country")
  dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

  analysis_name <- "under-5 | Maternal HIV burden + PMTCT | Complete annual countries"
  covariates <- c(
      "hdi_z", "gni_z", "sdi_z", "gbd_superregion",
      "hiv_pregnant_per_100_births_z", "pmtct_z"
  )
  # Retain rows with observed outcomes, population, exposure group, and DAG covariates.
  base_ok <- is.finite(panel$rate) & is.finite(panel$population) & panel$population > 0 &
      !is.na(panel$shortage) & !is.na(panel$gbd_superregion) &
      is.finite(panel$hdi_z) & is.finite(panel$gni_z) & is.finite(panel$sdi_z)
  eligible <- panel[base_ok, , drop = FALSE]
  exposure_ok <- is.finite(eligible$hiv_pregnant_per_100_births) &
      is.finite(eligible$pmtct_pct)
  ids_with_data <- unique(eligible$location_id[exposure_ok])

  # The primary sample keeps only countries with all 12 annual covariate records.
  complete_ids <- ids_with_data[vapply(ids_with_data, function(id) {
      observed_years <- eligible$year[
          eligible$location_id == id & exposure_ok
      ]
      all(2010:2021 %in% observed_years)
  }, logical(1))]
  analysis_data <- eligible[
      eligible$location_id %in% complete_ids & exposure_ok & eligible$year %in% 2010:2021,
      , drop = FALSE
  ]

  analysis_data$location_id <- factor(analysis_data$location_id)
  analysis_data$gbd_superregion <- factor(analysis_data$gbd_superregion)
  branch <- analysis_name

  # Fit the prespecified ordinary weighted LMM for this sample.
  selection <- fit_lmm(analysis_data, covariates, branch)
  if (!length(selection$fits)) stop("No converged model for ", analysis_name)
  models_to_report <- "LMM"

  effect_rows <- list()
  for (model_name in models_to_report) {
      fit <- selection$fits[[model_name]]
      if (is.null(fit)) next
      effect_rows[[model_name]] <- extract_terms(fit, branch, model_name)
  }

  write.csv(do.call(rbind, effect_rows), file.path(result_dir, "effects.csv"), row.names = FALSE)
  write.csv(data.frame(
      analysis = analysis_name, n_country_years = nrow(analysis_data),
      n_countries = length(unique(analysis_data$location_id)),
      shortage_countries = length(unique(analysis_data$location_id[analysis_data$shortage == 1])),
      comparison_countries = length(unique(analysis_data$location_id[analysis_data$shortage == 0])),
      years = paste(sort(unique(analysis_data$year)), collapse = ";")
  ), file.path(result_dir, "sample_flow.csv"), row.names = FALSE)
  write.csv(selection$models, file.path(result_dir, "lmm_fit_information.csv"), row.names = FALSE)
  write.csv(selection$attempts, file.path(result_dir, "lmm_fit_attempts.csv"), row.names = FALSE)
  write.csv(data.frame(analysis = analysis_name, model = "LMM"),
            file.path(result_dir, "model_specification.csv"), row.names = FALSE)
  message("Finished: ", analysis_name)



  # -----------------------------------------------------------------------------
  # 05_under5_observed_years
  # -----------------------------------------------------------------------------

  # Under-5 HIV incidence: observed country-years sensitivity analysis.
  # Required prepared model inputs are listed in the package README.

  set.seed(20260926)

  panel <- readRDS(file.path(processed_data_dir, "hiv", "hiv_under5_annual_panel.rds"))
  result_dir <- file.path(results_dir, "under-5", "Observed_years")
  dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

  analysis_name <- "under-5 | Maternal HIV burden + PMTCT | Observed country-years"
  covariates <- c(
      "hdi_z", "gni_z", "sdi_z", "gbd_superregion",
      "hiv_pregnant_per_100_births_z", "pmtct_z"
  )
  # Retain rows with observed outcomes, population, exposure group, and DAG covariates.
  base_ok <- is.finite(panel$rate) & is.finite(panel$population) & panel$population > 0 &
      !is.na(panel$shortage) & !is.na(panel$gbd_superregion) &
      is.finite(panel$hdi_z) & is.finite(panel$gni_z) & is.finite(panel$sdi_z)
  eligible <- panel[base_ok, , drop = FALSE]
  exposure_ok <- is.finite(eligible$hiv_pregnant_per_100_births) &
      is.finite(eligible$pmtct_pct)
  # Keep every country-year with observed treatment covariates.
  analysis_data <- eligible[exposure_ok & eligible$year %in% 2010:2021, , drop = FALSE]
  analysis_data$location_id <- factor(analysis_data$location_id)
  analysis_data$gbd_superregion <- factor(analysis_data$gbd_superregion)
  branch <- analysis_name

  # Fit the prespecified ordinary weighted LMM for this sample.
  selection <- fit_lmm(analysis_data, covariates, branch)
  if (!length(selection$fits)) stop("No converged model for ", analysis_name)
  models_to_report <- "LMM"

  effect_rows <- list()
  for (model_name in models_to_report) {
      fit <- selection$fits[[model_name]]
      if (is.null(fit)) next
      effect_rows[[model_name]] <- extract_terms(fit, branch, model_name)
  }

  write.csv(do.call(rbind, effect_rows), file.path(result_dir, "effects.csv"), row.names = FALSE)
  write.csv(data.frame(
      analysis = analysis_name, n_country_years = nrow(analysis_data),
      n_countries = length(unique(analysis_data$location_id)),
      shortage_countries = length(unique(analysis_data$location_id[analysis_data$shortage == 1])),
      comparison_countries = length(unique(analysis_data$location_id[analysis_data$shortage == 0])),
      years = paste(sort(unique(analysis_data$year)), collapse = ";")
  ), file.path(result_dir, "sample_flow.csv"), row.names = FALSE)
  write.csv(selection$models, file.path(result_dir, "lmm_fit_information.csv"), row.names = FALSE)
  write.csv(selection$attempts, file.path(result_dir, "lmm_fit_attempts.csv"), row.names = FALSE)
  write.csv(data.frame(analysis = analysis_name, model = "LMM"),
            file.path(result_dir, "model_specification.csv"), row.names = FALSE)
  message("Finished: ", analysis_name)



  # -----------------------------------------------------------------------------
  # 06_birth_pmm
  # -----------------------------------------------------------------------------

  # Birth-age HIV incidence: PMM missing-year sensitivity analysis.
  # Required prepared model inputs are listed in the package README.

  set.seed(20260926)

  panel <- readRDS(file.path(processed_data_dir, "hiv", "hiv_birth_annual_panel.rds"))
  scale_table <- read.csv(file.path(processed_data_dir, "hiv", "covariate_scaling.csv"))
  result_dir <- file.path(results_dir, "Birth", "PMM")
  dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

  analysis_name <- "Birth | Maternal HIV burden + PMTCT | PMM annual panel"
  covariates <- c(
      "hdi_z", "gni_z", "sdi_z", "gbd_superregion",
      "hiv_pregnant_per_100_births_z", "pmtct_z"
  )
  exposure_vars <- c("hiv_pregnant_per_100_births", "pmtct_pct")
  target_terms <- c("time_c:shortage", "shortage:post", "time_c:shortage:post")

  # Match the observed-data eligibility rules used in the two non-imputed analyses.
  base_ok <- is.finite(panel$rate) & is.finite(panel$population) & panel$population > 0 &
      !is.na(panel$shortage) & !is.na(panel$gbd_superregion) &
      is.finite(panel$hdi_z) & is.finite(panel$gni_z) & is.finite(panel$sdi_z)
  eligible <- panel[base_ok, , drop = FALSE]
  exposure_ok <- is.finite(eligible$hiv_pregnant_per_100_births) &
      is.finite(eligible$pmtct_pct)
  ids_with_data <- unique(eligible$location_id[exposure_ok])
  observed_data <- eligible[exposure_ok, , drop = FALSE]
  observed_data$location_id <- factor(observed_data$location_id)
  observed_data$gbd_superregion <- factor(observed_data$gbd_superregion)

  # PMM imputes only countries with genuine observations for both UNAIDS covariates.
  imputation_data <- eligible[
      eligible$location_id %in% ids_with_data & eligible$year %in% 2010:2021,
      , drop = FALSE
  ]
  imputation_data <- imputation_data[order(imputation_data$location_id, imputation_data$year), ]
  imputation_data$location_id <- factor(imputation_data$location_id)
  imputation_data$gbd_superregion <- factor(imputation_data$gbd_superregion)
  imputation_data$log_population <- log(imputation_data$population)
  imputation_data$imputation_row_id <- seq_len(nrow(imputation_data))
  rownames(imputation_data) <- as.character(imputation_data$imputation_row_id)

  imputed <- country_pmm(imputation_data, exposure_vars, m = N_IMPUTATIONS, k = 3L)
  completed_panels <- vector("list", N_IMPUTATIONS)
  for (m in seq_len(N_IMPUTATIONS)) {
      completed <- imputation_data
      for (variable in exposure_vars) {
          missing_rows <- imputed$indices[[variable]]
          if (length(missing_rows)) {
              completed[[variable]][missing_rows] <- imputed$draws[[variable]][, m]
          }
      }
      burden_scale <- scale_table[scale_table$covariate == "burden", ]
      pmtct_scale <- scale_table[scale_table$covariate == "pmtct", ]
      completed$hiv_pregnant_per_100_births_z <-
          (completed$hiv_pregnant_per_100_births - burden_scale$mean) / burden_scale$sd
      completed$pmtct_z <- (completed$pmtct_pct - pmtct_scale$mean) / pmtct_scale$sd
      completed_panels[[m]] <- completed
  }

  # Fit one ordinary LMM in each of the 50 completed panels.
  lmm_fits <- lapply(completed_panels, function(d) {
      tryCatch(fit_lme(d, covariates, model_label = "LMM"), error = function(e) e)
  })
  n_converged <- sum(vapply(lmm_fits, is_converged, logical(1)))
  if (n_converged < MIN_PMM_CONVERGED) {
      stop("Fewer than 45 of 50 ordinary PMM LMMs converged for ", analysis_name)
  }
  model_fits <- list(LMM = lmm_fits)
  pmm_fit_attempts <- data.frame(
      analysis = analysis_name, model = "LMM", converged_imputations = n_converged,
      requested_imputations = N_IMPUTATIONS)

  effect_rows <- list()
  for (model_name in names(model_fits)) {
      fits <- model_fits[[model_name]]
      term_estimates <- vector("list", N_IMPUTATIONS)
      for (m in seq_len(N_IMPUTATIONS)) {
          fit <- fits[[m]]
          if (!is_converged(fit)) next
          term_estimates[[m]] <- extract_terms(fit, analysis_name, model_name, m)
      }

      pooled_terms <- list()
      for (term in target_terms) {
          estimates <- standard_errors <- rep(NA_real_, N_IMPUTATIONS)
          for (m in seq_len(N_IMPUTATIONS)) {
              row <- term_estimates[[m]]
              if (is.null(row)) next
              row <- row[row$term == term, , drop = FALSE]
              if (nrow(row)) {
                  estimates[m] <- row$estimate[1]
                  standard_errors[m] <- row$se[1]
              }
          }
          pooled_terms[[term]] <- pool_coefficients(
              estimates, standard_errors, term, analysis_name, model_name,
              length(unique(imputation_data$location_id)), nrow(imputation_data),
              sum(vapply(fits, is_converged, logical(1))), N_IMPUTATIONS
          )
      }
      effect_rows[[model_name]] <- do.call(rbind, pooled_terms)
  }

  imputation_summary <- do.call(rbind, lapply(exposure_vars, function(variable) {
      draws <- imputed$draws[[variable]]
      data.frame(
          variable = variable,
          observed_country_years = sum(is.finite(imputation_data[[variable]])),
          imputed_country_years = nrow(draws),
          countries_with_any_observed = length(unique(imputation_data$location_id[
              is.finite(imputation_data[[variable]])
          ])),
          imputations = N_IMPUTATIONS,
          imputed_min = if (length(draws)) min(draws) else NA_real_,
          imputed_max = if (length(draws)) max(draws) else NA_real_
      )
  }))

  write.csv(do.call(rbind, effect_rows), file.path(result_dir, "effects.csv"), row.names = FALSE)
  write.csv(data.frame(
      analysis = analysis_name, n_country_years = nrow(imputation_data),
      n_countries = length(unique(imputation_data$location_id)),
      shortage_countries = length(unique(imputation_data$location_id[imputation_data$shortage == 1])),
      comparison_countries = length(unique(imputation_data$location_id[imputation_data$shortage == 0])),
      years = paste(sort(unique(imputation_data$year)), collapse = ";")
  ), file.path(result_dir, "sample_flow.csv"), row.names = FALSE)
  write.csv(imputation_summary, file.path(result_dir, "pmm_imputation_summary.csv"), row.names = FALSE)
  write.csv(pmm_fit_attempts, file.path(result_dir, "pmm_lmm_fit_attempts.csv"), row.names = FALSE)
  write.csv(data.frame(analysis = analysis_name, model = "LMM"),
            file.path(result_dir, "model_specification.csv"), row.names = FALSE)
  message("Finished: ", analysis_name)


  # -----------------------------------------------------------------------------
  # 07_under5_pmm
  # -----------------------------------------------------------------------------

  # Under-5 HIV incidence: PMM missing-year sensitivity analysis.
  # Required prepared model inputs are listed in the package README.

  set.seed(20260926)

  panel <- readRDS(file.path(processed_data_dir, "hiv", "hiv_under5_annual_panel.rds"))
  scale_table <- read.csv(file.path(processed_data_dir, "hiv", "covariate_scaling.csv"))
  result_dir <- file.path(results_dir, "under-5", "PMM")
  dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

  analysis_name <- "under-5 | Maternal HIV burden + PMTCT | PMM annual panel"
  covariates <- c(
      "hdi_z", "gni_z", "sdi_z", "gbd_superregion",
      "hiv_pregnant_per_100_births_z", "pmtct_z"
  )
  exposure_vars <- c("hiv_pregnant_per_100_births", "pmtct_pct")
  target_terms <- c("time_c:shortage", "shortage:post", "time_c:shortage:post")

  # Match the observed-data eligibility rules used in the two non-imputed analyses.
  base_ok <- is.finite(panel$rate) & is.finite(panel$population) & panel$population > 0 &
      !is.na(panel$shortage) & !is.na(panel$gbd_superregion) &
      is.finite(panel$hdi_z) & is.finite(panel$gni_z) & is.finite(panel$sdi_z)
  eligible <- panel[base_ok, , drop = FALSE]
  exposure_ok <- is.finite(eligible$hiv_pregnant_per_100_births) &
      is.finite(eligible$pmtct_pct)
  ids_with_data <- unique(eligible$location_id[exposure_ok])
  observed_data <- eligible[exposure_ok, , drop = FALSE]
  observed_data$location_id <- factor(observed_data$location_id)
  observed_data$gbd_superregion <- factor(observed_data$gbd_superregion)

  # PMM imputes only countries with genuine observations for both UNAIDS covariates.
  imputation_data <- eligible[
      eligible$location_id %in% ids_with_data & eligible$year %in% 2010:2021,
      , drop = FALSE
  ]
  imputation_data <- imputation_data[order(imputation_data$location_id, imputation_data$year), ]
  imputation_data$location_id <- factor(imputation_data$location_id)
  imputation_data$gbd_superregion <- factor(imputation_data$gbd_superregion)
  imputation_data$log_population <- log(imputation_data$population)
  imputation_data$imputation_row_id <- seq_len(nrow(imputation_data))
  rownames(imputation_data) <- as.character(imputation_data$imputation_row_id)

  imputed <- country_pmm(imputation_data, exposure_vars, m = N_IMPUTATIONS, k = 3L)
  completed_panels <- vector("list", N_IMPUTATIONS)
  for (m in seq_len(N_IMPUTATIONS)) {
      completed <- imputation_data
      for (variable in exposure_vars) {
          missing_rows <- imputed$indices[[variable]]
          if (length(missing_rows)) {
              completed[[variable]][missing_rows] <- imputed$draws[[variable]][, m]
          }
      }
      burden_scale <- scale_table[scale_table$covariate == "burden", ]
      pmtct_scale <- scale_table[scale_table$covariate == "pmtct", ]
      completed$hiv_pregnant_per_100_births_z <-
          (completed$hiv_pregnant_per_100_births - burden_scale$mean) / burden_scale$sd
      completed$pmtct_z <- (completed$pmtct_pct - pmtct_scale$mean) / pmtct_scale$sd
      completed_panels[[m]] <- completed
  }

  # Fit one ordinary LMM in each of the 50 completed panels.
  lmm_fits <- lapply(completed_panels, function(d) {
      tryCatch(fit_lme(d, covariates, model_label = "LMM"), error = function(e) e)
  })
  n_converged <- sum(vapply(lmm_fits, is_converged, logical(1)))
  if (n_converged < MIN_PMM_CONVERGED) {
      stop("Fewer than 45 of 50 ordinary PMM LMMs converged for ", analysis_name)
  }
  model_fits <- list(LMM = lmm_fits)
  pmm_fit_attempts <- data.frame(
      analysis = analysis_name, model = "LMM", converged_imputations = n_converged,
      requested_imputations = N_IMPUTATIONS)

  effect_rows <- list()
  for (model_name in names(model_fits)) {
      fits <- model_fits[[model_name]]
      term_estimates <- vector("list", N_IMPUTATIONS)
      for (m in seq_len(N_IMPUTATIONS)) {
          fit <- fits[[m]]
          if (!is_converged(fit)) next
          term_estimates[[m]] <- extract_terms(fit, analysis_name, model_name, m)
      }

      pooled_terms <- list()
      for (term in target_terms) {
          estimates <- standard_errors <- rep(NA_real_, N_IMPUTATIONS)
          for (m in seq_len(N_IMPUTATIONS)) {
              row <- term_estimates[[m]]
              if (is.null(row)) next
              row <- row[row$term == term, , drop = FALSE]
              if (nrow(row)) {
                  estimates[m] <- row$estimate[1]
                  standard_errors[m] <- row$se[1]
              }
          }
          pooled_terms[[term]] <- pool_coefficients(
              estimates, standard_errors, term, analysis_name, model_name,
              length(unique(imputation_data$location_id)), nrow(imputation_data),
              sum(vapply(fits, is_converged, logical(1))), N_IMPUTATIONS
          )
      }
      effect_rows[[model_name]] <- do.call(rbind, pooled_terms)
  }

  imputation_summary <- do.call(rbind, lapply(exposure_vars, function(variable) {
      draws <- imputed$draws[[variable]]
      data.frame(
          variable = variable,
          observed_country_years = sum(is.finite(imputation_data[[variable]])),
          imputed_country_years = nrow(draws),
          countries_with_any_observed = length(unique(imputation_data$location_id[
              is.finite(imputation_data[[variable]])
          ])),
          imputations = N_IMPUTATIONS,
          imputed_min = if (length(draws)) min(draws) else NA_real_,
          imputed_max = if (length(draws)) max(draws) else NA_real_
      )
  }))

  write.csv(do.call(rbind, effect_rows), file.path(result_dir, "effects.csv"), row.names = FALSE)
  write.csv(data.frame(
      analysis = analysis_name, n_country_years = nrow(imputation_data),
      n_countries = length(unique(imputation_data$location_id)),
      shortage_countries = length(unique(imputation_data$location_id[imputation_data$shortage == 1])),
      comparison_countries = length(unique(imputation_data$location_id[imputation_data$shortage == 0])),
      years = paste(sort(unique(imputation_data$year)), collapse = ";")
  ), file.path(result_dir, "sample_flow.csv"), row.names = FALSE)
  write.csv(imputation_summary, file.path(result_dir, "pmm_imputation_summary.csv"), row.names = FALSE)
  write.csv(pmm_fit_attempts, file.path(result_dir, "pmm_lmm_fit_attempts.csv"), row.names = FALSE)
  write.csv(data.frame(analysis = analysis_name, model = "LMM"),
            file.path(result_dir, "model_specification.csv"), row.names = FALSE)
  message("Finished: ", analysis_name)


  # Draw the HIV negative-control figure and copy the supplied HIV DAG rendering.
  # Run after the HIV model section of 04_additional_analyses.R has produced combined_effects.csv.


  # HIV-specific full adjusted LMMs: both endpoints and both CITS contrasts.
  library(ggplot2)
  library(patchwork)
  effects <- read.csv("results/hiv_specific/combined_effects.csv")
  d <- effects[effects$term %in% c("shortage:post", "time_c:shortage:post"), ]
  write.csv(d, "tables/table_s22_hiv_specific_level_slope.csv", row.names = FALSE)
  panels <- list()
  for (index in 1:4) {
    outcome <- if (index <= 2) "Birth" else "Under-5"
    term <- if (index %% 2 == 1) "shortage:post" else "time_c:shortage:post"
    panel <- d[d$outcome == outcome & d$term == term, ]
    panel$label <- paste0(panel$missingness_strategy, "\n", panel$n_countries,
      " countries; ", panel$n_country_years, " country-years")
    panel$label <- factor(panel$label, levels = rev(panel$label))
    panels[[index]] <- ggplot(panel, aes(x = estimate, y = label)) +
      geom_vline(xintercept = 0, linetype = 2, colour = "grey55") +
      geom_errorbar(aes(xmin = ci_lower, xmax = ci_upper), orientation = "y", width = 0.17) +
      geom_point(size = 2.2, colour = "#176B9A") +
      labs(title = paste(LETTERS[index], outcome, "-",
        if (index %% 2 == 1) "level change" else "slope change"), y = NULL,
        x = if (index %% 2 == 1) "Cases per 1,000 age-specific population" else
          "Cases per 1,000 age-specific population per year") + publication_theme(10)
  }
  figure <- (panels[[1]] | panels[[2]]) / (panels[[3]] | panels[[4]])
  ggsave("figures/Figure_S15_HIV_negative_controls.png", figure, width = 12, height = 7.8, dpi = 300)
  ggsave("figures/Figure_S15_HIV_negative_controls.pdf", figure, width = 12, height = 7.8,
         device = grDevices::cairo_pdf)

  # This is an explicit assumed causal graph, not evidence that all common
  # causes were measured. Maternal PMTCT coverage is the observed access proxy.
  library(dagitty)
  graph_text <- 'dag {
   BPG_shortage [exposure,pos="0,1"]
   HIV_incidence [outcome,pos="8,1"]
   Shared_supply_disruption [latent,pos="1,3"]
   ARV_stockout [latent,pos="4,3"]
   PMTCT_coverage [adjusted,pos="8,3"]
   Socioeconomic_domain [adjusted,pos="4,5"]
   Maternal_HIV_burden [adjusted,pos="4,-1"]
   Country_region_time [adjusted,pos="-3,5"]
   Shared_supply_disruption -> BPG_shortage
   Shared_supply_disruption -> ARV_stockout
   ARV_stockout -> PMTCT_coverage
   PMTCT_coverage -> HIV_incidence
   Socioeconomic_domain -> BPG_shortage
   Socioeconomic_domain -> HIV_incidence
   Socioeconomic_domain -> PMTCT_coverage
   Maternal_HIV_burden -> BPG_shortage
   Maternal_HIV_burden -> HIV_incidence
   Maternal_HIV_burden -> PMTCT_coverage
   Country_region_time -> BPG_shortage
   Country_region_time -> HIV_incidence
   Country_region_time -> ARV_stockout
  }'
  dir.create("data/dag", recursive = TRUE, showWarnings = FALSE)
  writeLines(graph_text, "data/dag/hiv_specific_adjustment_dag.txt")
  g <- dagitty(graph_text)
  adjustment <- c("PMTCT_coverage", "Socioeconomic_domain", "Maternal_HIV_burden", "Country_region_time")
  valid <- isAdjustmentSet(g, adjustment)
  writeLines(c(paste("Adjustment set valid under specified graph:", valid),
    "This is conditional on the drawn arrows and absence of omitted direct paths.",
    "PMTCT is an access proxy, not direct ARV stock-out data.",
    "Socioeconomic domain is represented by the retained HDI/GNI/SDI terms;",
    "the graph does not require all three indicators as separate causal domains."),
    "results/hiv_specific/dag_adjustment_check.txt")
  stopifnot(valid)
  # The publication panels use the official DAGitty web renderer's classic style.
  # Keep the exported SVG, PNG and PDF together in figures/dagitty_web.
  # Paste data/dag/figure_s1b_web_dagitty.txt into dagitty.net/dags.html
  # to inspect the same graph and adjustment markings.
  file.copy("figures/dagitty_web/figure_s1b.png",
            "figures/Figure_S1B_HIV_adjustment_DAG.png", overwrite = TRUE)
  file.copy("figures/dagitty_web/figure_s1b.pdf",
            "figures/Figure_S1B_HIV_adjustment_DAG.pdf", overwrite = TRUE)
})

} else {
  message('HIV models skipped because the prepared annual panels and covariate scaling file are absent from results/hiv_specific/Data/.')
}
publication_result_files <- c(
  'results/original_primary_lmm_coefficients.csv',
  'results/counterfactual_inference_audit.csv',
  'results/counterfactual_covariance_ratios.csv',
  'results/cs_sensitivity_models.csv',
  'results/who_anc_service_missingness.csv',
  'results/who_anc_complete_series_counts.csv',
  'results/who_service_sensitivity_lmm.csv',
  'results/who_service_full_sample_pmm_lmm.csv',
  'results/hiv_specific/combined_effects.csv',
  'results/qba_tipping_point_summary.csv',
  'results/included_vs_unclassified_summary.csv',
  'results/included_vs_unclassified_categories.csv',
  'results/who_gbd_agreement_statistics.csv',
  'results/who_gbd_group_summary.csv',
  'results/who_gbd_annual_weighted_rates.csv',
  'results/country_label_permutation_summary.csv',
  'results/iqvia/iqvia_incident_ratio_lmm_results.csv',
  'results/iqvia/iqvia_log_ratio_group_lmm.csv',
  'results/iqvia/iqvia_log_ratio_group_description.csv',
  'results/iqvia/iqvia_shortage_group_all_contrasts.csv',
  'results/iqvia/iqvia_shortage_cits_lmm.csv',
  'results/iqvia/iqvia_covariate_sensitivity.csv',
  'results/iqvia/iqvia_model_structure_covariate_sensitivity.csv',
  'results/iqvia/iqvia_2014_2020_cs_association.csv',
  'results/iqvia/iqvia_2014_2020_group_description.csv',
  'results/iqvia/iqvia_2014_2020_group_all_contrasts.csv',
  'results/iqvia/iqvia_incident_ratio_coverage_by_year.csv'
)
missing_publication_results <- publication_result_files[!file.exists(publication_result_files)]
if (length(missing_publication_results) == 0L) {


# -----------------------------------------------------------------------------
# Final result comparison, publication tables, and figure-asset audit
# -----------------------------------------------------------------------------
local({
  # Compare the rebuilt original primary model with the saved v3 coefficients.
  # Run this after the primary-analysis section of 02_primary_and_subgroup_analysis.R.

  reference_file <- file.path("logs", "original_primary_reference.csv")
  current_file <- file.path("results", "original_primary_lmm_coefficients.csv")

  if (!file.exists(reference_file)) {
    stop("The small v3 coefficient reference file is missing from logs/.")
  }
  if (!file.exists(current_file)) {
    stop("Run the primary-analysis section of 02_primary_and_subgroup_analysis.R before comparing coefficients.")
  }

  reference <- read.csv(reference_file, stringsAsFactors = FALSE)
  current <- read.csv(current_file, stringsAsFactors = FALSE)

  comparison <- merge(
    reference[, c("Parameter", "Estimate", "SE", "p_value")],
    current[, c("term", "Estimate", "std_error", "p_value")],
    by.x = "Parameter",
    by.y = "term",
    all = TRUE,
    suffixes = c("_v3", "_v7")
  )

  comparison$estimate_difference <- comparison$Estimate_v7 - comparison$Estimate_v3
  comparison$se_difference <- comparison$std_error - comparison$SE
  comparison$p_value_difference <- comparison$p_value_v7 - comparison$p_value_v3

  write.csv(
    comparison,
    file.path("logs", "original_primary_result_comparison.csv"),
    row.names = FALSE,
    na = ""
  )

  max_estimate_difference <- max(abs(comparison$estimate_difference), na.rm = TRUE)
  max_se_difference <- max(abs(comparison$se_difference), na.rm = TRUE)
  max_p_difference <- max(abs(comparison$p_value_difference), na.rm = TRUE)

  cat("Maximum absolute coefficient difference:", max_estimate_difference, "\n")
  cat("Maximum absolute standard-error difference:", max_se_difference, "\n")
  cat("Maximum absolute P-value difference:", max_p_difference, "\n")

  # The archived model was fit through lmerTest::step(); tiny optimizer and
  # refit differences can occur even when the data and selected formula match.
  # These tolerances are tighter than the precision used in the manuscript.
  # The standard-error tolerance is 2e-5 because the maximum observed difference
  # from the archived lmerTest step/refit workflow is 1.02e-5.
  if (max_estimate_difference > 1e-4 || max_se_difference > 2e-5 || max_p_difference > 1e-4) {
    warning("One or more primary model values do not reproduce the v3 reference.")
  } else {
    cat("All fixed-effect estimates, standard errors, and P values match the v3 reference within the stated numerical tolerances.\n")
  }


  # -----------------------------------------------------------------------------
  # Publication tables
  # -----------------------------------------------------------------------------

  # Verify the primary result, prepare display-ready tables, and audit figure assets.
  # Full-precision analysis outputs remain unchanged in results/.


  # Table 1. Original weighted congenital-syphilis CITS coefficients.
  primary <- read.csv(file.path(results_dir, "original_primary_lmm_coefficients.csv"))
  primary_table <- primary
  primary_table$Estimate <- round(primary_table$Estimate, 3)
  primary_table$std_error <- round(primary_table$std_error, 3)
  primary_table$ci_lower <- round(primary_table$ci_lower, 3)
  primary_table$ci_upper <- round(primary_table$ci_upper, 3)
  primary_table$p_value <- signif(primary_table$p_value, 3)
  write.csv(primary_table, file.path(tables_dir, "table_primary_cits.csv"), row.names = FALSE)

  # Supplementary Table S7: corrected covariance-aware counterfactual contrasts.
  absolute <- read.csv(file.path(results_dir, "counterfactual_inference_audit.csv"))
  relative <- read.csv(file.path(results_dir, "counterfactual_covariance_ratios.csv"))
  write.csv(absolute, file.path(tables_dir, "table_s07_absolute_contrasts.csv"), row.names = FALSE)
  write.csv(relative, file.path(tables_dir, "table_s07_relative_contrasts.csv"), row.names = FALSE)

  # Table 2. Congenital-syphilis sensitivity estimates.
  cs <- read.csv(file.path(results_dir, "cs_sensitivity_models.csv"))
  cs_table <- cs
  for (column in c("estimate", "std_error", "ci_lower", "ci_upper")) {
    if (column %in% names(cs_table)) cs_table[[column]] <- round(cs_table[[column]], 3)
  }
  if ("p_value" %in% names(cs_table)) cs_table$p_value <- signif(cs_table$p_value, 3)
  write.csv(cs_table, file.path(tables_dir, "table_cs_sensitivity.csv"), row.names = FALSE)

  # Table 3. WHO service missingness and country-year availability.
  who_missing <- read.csv(file.path(results_dir, "who_anc_service_missingness.csv"))
  who_complete <- read.csv(file.path(results_dir, "who_anc_complete_series_counts.csv"))
  write.csv(who_missing, file.path(tables_dir, "table_who_service_availability.csv"), row.names = FALSE)
  write.csv(who_complete, file.path(tables_dir, "table_who_complete_series.csv"), row.names = FALSE)

  # Table 4. WHO observed-data and full-sample imputed sensitivity models.
  who_observed <- read.csv(file.path(results_dir, "who_service_sensitivity_lmm.csv"))
  who_imputed <- read.csv(file.path(results_dir, "who_service_full_sample_pmm_lmm.csv"))
  write.csv(who_observed, file.path(tables_dir, "table_who_service_observed_models.csv"), row.names = FALSE)
  write.csv(who_imputed, file.path(tables_dir, "table_who_service_full_sample_PMM.csv"), row.names = FALSE)

  # Supplementary Table S22. Only the HIV-specific adjusted model, both outcomes.
  hiv <- read.csv(file.path(results_dir, "hiv_specific", "combined_effects.csv"))
  hiv <- hiv[hiv$term %in% c("shortage:post", "time_c:shortage:post"), ]
  write.csv(hiv, file.path(tables_dir, "table_s22_hiv_specific_level_slope.csv"), row.names = FALSE)

  # Table 6. Quantitative-bias calibration and country selection.
  qba <- read.csv(file.path(results_dir, "qba_tipping_point_summary.csv"))
  selection <- read.csv(file.path(results_dir, "included_vs_unclassified_summary.csv"))
  selection_categories <- read.csv(file.path(results_dir, "included_vs_unclassified_categories.csv"))
  write.csv(qba, file.path(tables_dir, "table_qba_tipping_point.csv"), row.names = FALSE)
  write.csv(selection, file.path(tables_dir, "table_country_selection_continuous.csv"), row.names = FALSE)
  write.csv(selection_categories, file.path(tables_dir, "table_country_selection_categories.csv"), row.names = FALSE)

  # Table 7. GBD and WHO surveillance triangulation statistics.
  surveillance <- read.csv(file.path(results_dir, "who_gbd_agreement_statistics.csv"))
  write.csv(surveillance, file.path(tables_dir, "table_gbd_surveillance_validation.csv"), row.names = FALSE)
  surveillance_groups <- read.csv(file.path(results_dir, "who_gbd_group_summary.csv"))
  write.csv(surveillance_groups, file.path(tables_dir, "table_who_gbd_group_comparison.csv"), row.names = FALSE)
  surveillance_years <- read.csv(file.path(results_dir, "who_gbd_annual_weighted_rates.csv"))
  write.csv(surveillance_years, file.path(tables_dir, "table_who_gbd_annual_rates.csv"), row.names = FALSE)

  # Table 9. Country-label permutation diagnostics.
  permutation <- read.csv(file.path(results_dir, "country_label_permutation_summary.csv"))
  write.csv(permutation, file.path(tables_dir, "table_country_label_permutation.csv"), row.names = FALSE)

  # Tables 10-14. IQVIA estimates from the rebuilt product-quarter pipeline.
  iqvia_results_dir <- file.path(results_dir, "iqvia")
  iqvia_main <- read.csv(file.path(iqvia_results_dir, "iqvia_incident_ratio_lmm_results.csv"))
  iqvia_group <- read.csv(file.path(iqvia_results_dir, "iqvia_log_ratio_group_lmm.csv"))
  iqvia_group_summary <- read.csv(file.path(iqvia_results_dir, "iqvia_log_ratio_group_description.csv"))
  iqvia_group_contrasts <- read.csv(file.path(iqvia_results_dir, "iqvia_shortage_group_all_contrasts.csv"))
  iqvia_cits <- read.csv(file.path(iqvia_results_dir, "iqvia_shortage_cits_lmm.csv"))
  iqvia_covariates <- read.csv(file.path(iqvia_results_dir, "iqvia_covariate_sensitivity.csv"))
  iqvia_structure <- read.csv(file.path(iqvia_results_dir, "iqvia_model_structure_covariate_sensitivity.csv"))
  iqvia_window <- read.csv(file.path(iqvia_results_dir, "iqvia_2014_2020_cs_association.csv"))
  iqvia_window_description <- read.csv(file.path(iqvia_results_dir, "iqvia_2014_2020_group_description.csv"))
  iqvia_window_group <- read.csv(file.path(iqvia_results_dir, "iqvia_2014_2020_group_all_contrasts.csv"))
  iqvia_coverage <- read.csv(file.path(iqvia_results_dir, "iqvia_incident_ratio_coverage_by_year.csv"))
  iqvia_incidence_validation <- read.csv(file.path(iqvia_results_dir,
      "iqvia_gbd_incidence_strata_validation.csv"))
  write.csv(iqvia_main, file.path(tables_dir, "table_iqvia_main_and_denominator_sensitivity.csv"), row.names = FALSE)
  write.csv(iqvia_group_summary, file.path(tables_dir, "table_iqvia_shortage_group_description.csv"), row.names = FALSE)
  write.csv(iqvia_group, file.path(tables_dir, "table_iqvia_shortage_group_lmm.csv"), row.names = FALSE)
  write.csv(iqvia_group_contrasts, file.path(tables_dir, "table_iqvia_all_group_contrasts_2013_2021.csv"), row.names = FALSE)
  write.csv(iqvia_cits, file.path(tables_dir, "table_iqvia_shortage_cits.csv"), row.names = FALSE)
  write.csv(iqvia_covariates, file.path(tables_dir, "table_iqvia_covariate_sensitivity.csv"), row.names = FALSE)
  write.csv(iqvia_structure, file.path(tables_dir, "table_iqvia_model_structure_sensitivity.csv"), row.names = FALSE)
  write.csv(iqvia_window, file.path(tables_dir, "table_iqvia_2014_2020_association.csv"), row.names = FALSE)
  write.csv(iqvia_window_description, file.path(tables_dir, "table_iqvia_group_description_2014_2020.csv"), row.names = FALSE)
  write.csv(iqvia_window_group, file.path(tables_dir, "table_iqvia_all_group_contrasts_2014_2020.csv"), row.names = FALSE)
  write.csv(iqvia_coverage, file.path(tables_dir, "table_iqvia_annual_completeness.csv"), row.names = FALSE)
  write.csv(iqvia_incidence_validation,
      file.path(tables_dir, "table_iqvia_gbd_incidence_strata_validation.csv"),
      row.names = FALSE)

  cat("Saved display-ready tables from the current v7 result files.\n")


  # -----------------------------------------------------------------------------
  # Figure asset audit
  # -----------------------------------------------------------------------------

  # Compare the expected figure inventory with files available in figures/.
  # Check that expected publication figure assets are present; this does not
  # establish that each asset was regenerated from the included R scripts.

  figures_dir <- "figures"
  logs_dir <- "logs"
  index_file <- file.path("logs", "publication_figure_crosswalk_20261014.csv")
  figure_index <- read.csv(index_file, stringsAsFactors = FALSE, check.names = FALSE)

  asset_status <- character(nrow(figure_index))
  for (row in seq_len(nrow(figure_index))) {
    expected_files <- strsplit(figure_index$Expected_output[row], ";", fixed = TRUE)[[1]]
    expected_files <- expected_files[nzchar(expected_files)]
    output_paths <- file.path(figures_dir, expected_files)
    found <- file.exists(output_paths)
    if (length(expected_files) == 0L) {
      asset_status[row] <- "No figure asset listed"
    } else if (all(found)) {
      asset_status[row] <- "All expected figure assets present"
    } else if (any(found)) {
      asset_status[row] <- "Some expected figure assets present"
    } else {
      asset_status[row] <- "Expected figure assets missing"
    }
  }

  figure_index$Asset_status <- asset_status
  write.csv(figure_index,
    file.path(logs_dir, "publication_figure_asset_audit_20261014.csv"),
    row.names = FALSE, na = "", fileEncoding = "UTF-8"
  )

  cat("Publication items checked:", nrow(figure_index), "\n")
  print(table(figure_index$Asset_status))
})
} else {
  message("Final publication tables not rebuilt; complete the required analysis sections first. Missing files: ",
    paste(missing_publication_results, collapse = ", "))
}
