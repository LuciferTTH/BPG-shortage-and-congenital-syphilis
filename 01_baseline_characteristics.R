# 0. File information --------------------------------------------------------
# File Name: 01_baseline_characteristics.R
# Description: Data extraction and analysis for congenital syphilis incidence
#            in relation to penicillin shortages
# Author: Huang Tao
# Date: 2025-11-17
# Version: 1.0
# R version: 4.3.3


# 1. Load required libraries and define functions ----------------------------
## Load libraries ------------------------------------------------------------
library(extrafont)
library(openxlsx)
library(dplyr)
library(ggplot2)
library(httpgd)
library(showtext)
library(ggsci)
library(scales)
library(survey)
library(its.analysis)
library(lmtest) # For robust standard errors
library(sandwich) # For clustered standard errors

## Set working directory -----------------------------------------------------
getwd()

## Define theme for ggplot2 --------------------------------------------------
theme_jama <- function() {
    theme_minimal(base_family = "sans") +
        theme(
            # Base elements
            text = element_text(color = "black", size = 10),
            line = element_line(color = "black", linewidth = 0.25),

            # Axes
            axis.line = element_line(color = "black"),
            axis.ticks = element_line(color = "black"),
            axis.title = element_text(face = "bold"),
            axis.title.x = element_text(margin = margin(t = 10)),
            axis.title.y = element_text(margin = margin(r = 10)),

            # Gridlines (removed)
            panel.grid = element_blank(),

            # Plot title/subtitle
            plot.title = element_text(face = "bold", size = 12, hjust = 0),
            plot.subtitle = element_text(face = "italic", hjust = 0),

            # Legend (top position to save space)
            legend.position = "top",
            legend.title = element_text(face = "bold"),

            # Panel and background
            panel.background = element_rect(fill = "white", color = NA),
            plot.background = element_rect(fill = "white", color = NA),
            strip.background = element_rect(fill = "gray95", color = "black")
        )
}

# 2. Load data and preprocess ------------------------------------------------
load("data_combined.Rdata")
nsummary <- data_combined %>%
    group_by(location_name) %>%
    summarise(n = n())

# 3. Data visualization and exploratory analysis -----------------------------
## 3.1 change ----------------------------------------------------------------
data_syphilis_summary_all <-
    data_combined %>%
    group_by(year) %>%
    summarise(
        total_number = sum(val, na.rm = TRUE),
        total_population = sum(population, na.rm = TRUE),
        proportion = total_number / total_population * 1000
    ) %>%
    ungroup()

str(data_combined)
x_factorvars <- c("all", "income_level", "hdi_level", "sdi_level2021", "gbd_superregion")

# all countries
results <- do.call(rbind, lapply(x_factorvars, function(var) {
    data <- data_combined
    data_sub <- data %>% filter(year == 2021)
    if (var == "all") {
        data_sub$subgroup <- "all"
    } else {
        data_sub$subgroup <- data_sub[[var]]
    }
    data <- data %>%
        left_join(
            data_sub %>% select("location_id", "subgroup"),
            by = "location_id"
        )

    results_sum <- data %>%
        group_by(year, subgroup) %>%
        summarise(
            n = n(),
            total_number = sum(val, na.rm = TRUE),
            total_population = sum(population, na.rm = TRUE),
            incidence_rate = total_number / total_population * 1000
        ) %>%
        filter(year == 2010 | year == 2021) %>%
        select(year, subgroup, n, incidence_rate) %>%
        tidyr::pivot_wider(
            names_from = year,
            values_from = incidence_rate,
            names_prefix = "year_"
        ) %>%
        mutate(
            subgroup = as.character(subgroup),
            change = (year_2021 - year_2010) / year_2010 * 100
        )
    results_sum <- rbind(
        c(var, NA, NA, NA, NA),
        results_sum
    )
    return(results_sum)
}))
View(results)
write.xlsx(results,
    file = "Results/syphilis_congenital_incidence_change_by_subgroups_2010_2021_all.xlsx",
    rowNames = FALSE
)

# countries with shortage reports
results_shortage <- do.call(rbind, lapply(x_factorvars, function(var) {
    data <- data_combined %>% filter(shortage == 1)
    data_sub <- data %>% filter(year == 2021)
    if (var == "all") {
        data_sub$subgroup <- "all"
    } else {
        data_sub$subgroup <- data_sub[[var]]
    }
    data <- data %>%
        left_join(
            data_sub %>% select("location_id", "subgroup"),
            by = "location_id"
        )

    results_sum <- data %>%
        group_by(year, subgroup) %>%
        summarise(
            n = n(),
            total_number = sum(val, na.rm = TRUE),
            total_population = sum(population, na.rm = TRUE),
            incidence_rate = total_number / total_population * 1000
        ) %>%
        filter(year == 2010 | year == 2021) %>%
        select(year, subgroup, n, incidence_rate) %>%
        tidyr::pivot_wider(
            names_from = year,
            values_from = incidence_rate,
            names_prefix = "year_"
        ) %>%
        mutate(
            subgroup = as.character(subgroup),
            change = (year_2021 - year_2010) / year_2010 * 100
        )
    results_sum <- rbind(
        c(var, NA, NA, NA, NA),
        results_sum
    )
    return(results_sum)
}))
View(results_shortage)
write.xlsx(results_shortage,
    file = "Results/syphilis_congenital_incidence_change_by_subgroups_2010_2021_shortage.xlsx",
    rowNames = FALSE
)

# countries without shortage reports
results_noshortage <- do.call(rbind, lapply(x_factorvars, function(var) {
    data <- data_combined %>% filter(shortage == 0)
    data_sub <- data %>% filter(year == 2021)
    if (var == "all") {
        data_sub$subgroup <- "all"
    } else {
        data_sub$subgroup <- data_sub[[var]]
    }

    data <- data %>%
        left_join(
            data_sub %>% select("location_id", "subgroup"),
            by = "location_id"
        )

    results_sum <- data %>%
        group_by(year, subgroup) %>%
        summarise(
            n = n(),
            total_number = sum(val, na.rm = TRUE),
            total_population = sum(population, na.rm = TRUE),
            incidence_rate = total_number / total_population * 1000
        ) %>%
        filter(year == 2010 | year == 2021) %>%
        select(year, subgroup, n, incidence_rate) %>%
        tidyr::pivot_wider(
            names_from = year,
            values_from = incidence_rate,
            names_prefix = "year_"
        ) %>%
        mutate(
            subgroup = as.character(subgroup),
            change = (year_2021 - year_2010) / year_2010 * 100
        )
    results_sum <- rbind(
        c(var, NA, NA, NA, NA),
        results_sum
    )
    return(results_sum)
}))
View(results_noshortage)
write.xlsx(results_noshortage,
    file = "Results/syphilis_congenital_incidence_change_by_subgroups_2010_2021_noshortage.xlsx",
    rowNames = FALSE
)

## 3.2 change 20110-2015 2016-2021 ----------------------------------------------------------------
data_syphilis_summary_all <-
    data_combined %>%
    group_by(year) %>%
    summarise(
        total_number = sum(val, na.rm = TRUE),
        total_population = sum(population, na.rm = TRUE),
        proportion = total_number / total_population * 1000
    ) %>%
    ungroup()

str(data_combined)
x_factorvars <- c("all", "income_level", "hdi_level", "sdi_level2021", "gbd_superregion")

# all countries
results <- do.call(rbind, lapply(x_factorvars, function(var) {
    data <- data_combined
    data_sub <- data %>% filter(year == 2021)
    if (var == "all") {
        data_sub$subgroup <- "all"
    } else {
        data_sub$subgroup <- data_sub[[var]]
    }
    data <- data %>%
        left_join(
            data_sub %>% select("location_id", "subgroup"),
            by = "location_id"
        )

    results_sum <- data %>%
        group_by(year, subgroup) %>%
        summarise(
            n = n(),
            total_number = sum(val, na.rm = TRUE),
            total_population = sum(population, na.rm = TRUE),
            incidence_rate = total_number / total_population * 1000
        ) %>%
        filter(year == 2010 | year == 2015 | year == 2016 | year == 2021) %>%
        select(year, subgroup, n, incidence_rate) %>%
        tidyr::pivot_wider(
            names_from = year,
            values_from = incidence_rate,
            names_prefix = "year_"
        ) %>%
        mutate(
            subgroup = as.character(subgroup),
            change = (year_2021 - year_2010) / year_2010 * 100,
            change10to15 = (year_2015 - year_2010) / year_2010 * 100,
            change16to21 = (year_2021 - year_2016) / year_2016 * 100
        )
    results_sum <- rbind(
        c(var, NA, NA, NA, NA),
        results_sum
    )
    return(results_sum)
}))
View(results)
write.xlsx(results,
    file = "Results/syphilis_congenital_incidence_change_by_subgroups_2010_2021_all.xlsx",
    rowNames = FALSE
)

# countries with shortage reports
results_shortage <- do.call(rbind, lapply(x_factorvars, function(var) {
    data <- data_combined %>% filter(shortage == 1)
    data_sub <- data %>% filter(year == 2021)
    if (var == "all") {
        data_sub$subgroup <- "all"
    } else {
        data_sub$subgroup <- data_sub[[var]]
    }
    data <- data %>%
        left_join(
            data_sub %>% select("location_id", "subgroup"),
            by = "location_id"
        )

    results_sum <- data %>%
        group_by(year, subgroup) %>%
        summarise(
            n = n(),
            total_number = sum(val, na.rm = TRUE),
            total_population = sum(population, na.rm = TRUE),
            incidence_rate = total_number / total_population * 1000
        ) %>%
        filter(year == 2010 | year == 2021) %>%
        select(year, subgroup, n, incidence_rate) %>%
        tidyr::pivot_wider(
            names_from = year,
            values_from = incidence_rate,
            names_prefix = "year_"
        ) %>%
        mutate(
            subgroup = as.character(subgroup),
            change = (year_2021 - year_2010) / year_2010 * 100
        )
    results_sum <- rbind(
        c(var, NA, NA, NA, NA),
        results_sum
    )
    return(results_sum)
}))
View(results_shortage)
write.xlsx(results_shortage,
    file = "Results/syphilis_congenital_incidence_change_by_subgroups_2010_2021_shortage.xlsx",
    rowNames = FALSE
)

# countries without shortage reports
results_noshortage <- do.call(rbind, lapply(x_factorvars, function(var) {
    data <- data_combined %>% filter(shortage == 0)
    data_sub <- data %>% filter(year == 2021)
    if (var == "all") {
        data_sub$subgroup <- "all"
    } else {
        data_sub$subgroup <- data_sub[[var]]
    }

    data <- data %>%
        left_join(
            data_sub %>% select("location_id", "subgroup"),
            by = "location_id"
        )

    results_sum <- data %>%
        group_by(year, subgroup) %>%
        summarise(
            n = n(),
            total_number = sum(val, na.rm = TRUE),
            total_population = sum(population, na.rm = TRUE),
            incidence_rate = total_number / total_population * 1000
        ) %>%
        filter(year == 2010 | year == 2021) %>%
        select(year, subgroup, n, incidence_rate) %>%
        tidyr::pivot_wider(
            names_from = year,
            values_from = incidence_rate,
            names_prefix = "year_"
        ) %>%
        mutate(
            subgroup = as.character(subgroup),
            change = (year_2021 - year_2010) / year_2010 * 100
        )
    results_sum <- rbind(
        c(var, NA, NA, NA, NA),
        results_sum
    )
    return(results_sum)
}))
View(results_noshortage)
write.xlsx(results_noshortage,
    file = "Results/syphilis_congenital_incidence_change_by_subgroups_2010_2021_noshortage.xlsx",
    rowNames = FALSE
)

## 3.2 Overall trend by shortage groups --------------------------------------
data_syphilis_summary <-
    data_combined %>%
    mutate(shortage = ifelse(is.na(shortage), 2, shortage)) %>%
    group_by(year, shortage) %>%
    summarise(
        total_number = sum(val, na.rm = TRUE),
        total_population = sum(population, na.rm = TRUE),
        proportion = total_number / total_population
    ) %>%
    ungroup()

View(data_syphilis_summary)
head(data_syphilis_summary)

incidence_plot <- ggplot(
    data = data_syphilis_summary,
    mapping = aes(x = year, y = proportion * 100, colour = factor(shortage))
) +
    geom_line() +
    geom_point(size = 4) +
    geom_vline(xintercept = c(2014, 2016), linetype = "dashed", color = "grey") +
    labs(
        title = "Congenital Syphilis Incidence Rate",
        x = "Year",
        y = "Incidence rate (cases per 1000 births)"
    ) +
    theme_jama() +
    scale_colour_manual(
        name = "Shortage between 2014 and 2016",
        values = c("blue", "red", "grey"),
        labels = c("No Shortage", "Shortage", "Unknown")
    ) +
    # scale_y_continuous(labels = scales::percent_format(scale = 1)) +
    scale_x_continuous(breaks = seq(2010, 2021, by = 1), labels = seq(2010, 2021, by = 1))
incidence_plot
ggsave("syphilis_congenital_incidence_plot.png",
    plot = incidence_plot,
    width = 10, height = 6, dpi = 300
)

### Subgroup -----------------------------------------------------------------
#### Prepare aggregated data by shortage groups ------------------------------
data_ts_analysis <- data_combined %>%
    filter(
        year >= 2010 & year <= 2021,
        !is.na(shortage)
    ) %>%
    group_by(year, shortage) %>%
    summarise(
        total_cases = sum(val, na.rm = TRUE),
        total_population = sum(population, na.rm = TRUE),
        proportion = (total_cases / total_population) * 1000,
        .groups = "drop"
    ) %>%
    mutate(
        time = year - 2010,
        intervention = ifelse(year >= 2016, 1, 0),
        postslope = ifelse(year >= 2016, year - 2016, 0)
    )

ggplot(
    data = data_ts_analysis,
    mapping = aes(x = year, y = proportion, colour = factor(shortage))
) +
    geom_line() +
    geom_point(size = 4) +
    geom_vline(xintercept = c(2014, 2016), linetype = "dashed", color = "grey") +
    labs(
        title = "Congenital Syphilis Incidence Rate",
        x = "Year",
        y = "Incidence rate (cases per 1000 births)"
    ) +
    theme_jama() +
    scale_colour_manual(
        name = "Shortage between 2014 and 2016",
        values = c("blue", "red", "grey"),
        labels = c("No Shortage", "Shortage", "Unknown")
    ) +
    scale_y_continuous(breaks = seq(0, 5, by = 1), limits = c(0, 5)) +
    scale_x_continuous(breaks = seq(1990, 2021, by = 1), labels = seq(1990, 2021, by = 1))
ggsave("syphilis_congenital_incidence_plot_its_all.png", width = 10, height = 6, dpi = 300)

#### By GBD regions ----------------------------------------------------------
data_syphilis_sum <- data_combined %>%
    filter(year >= 2010 & year <= 2021) %>%
    mutate(shortage = ifelse(is.na(shortage), 2, shortage)) %>%
    group_by(year, gbd_superregion, shortage) %>%
    summarise(
        total_number = sum(val, na.rm = TRUE),
        total_population = sum(population, na.rm = TRUE),
        proportion = total_number / total_population
    ) %>%
    ungroup()

ggplot(
    data = data_syphilis_sum,
    mapping = aes(x = year, y = proportion * 1000, colour = factor(shortage))
) +
    geom_line() +
    geom_point(size = 1) +
    geom_vline(xintercept = c(2014, 2016), linetype = "dashed", color = "grey") +
    labs(
        # title = "Congenital Syphilis Incidence Rate by Regions",
        x = "Year",
        y = "Incidence rate (cases per 1000 births)"
    ) +
    theme_jama() +
    scale_colour_manual(
        name = "Shortage between 2014 and 2016",
        values = c("blue", "red", "grey"),
        labels = c("No Shortage", "Shortage", "Unknown")
    ) +
    # scale_y_continuous(breaks = seq(0, 10, by = 1), limits = c(0, 10)) +
    scale_x_continuous(breaks = seq(1990, 2021, by = 1), labels = seq(1990, 2021, by = 1)) +
    facet_wrap(~gbd_superregion)

ggsave("syphilis_congenital_incidence_plot_its_gbdregions.png", width = 12, height = 8, dpi = 300)

#### SDI quintiles in 2021 ---------------------------------------------------
data_syphilis_summary <-
    data_combined %>%
    mutate(shortage = ifelse(is.na(shortage), 2, shortage)) %>%
    rename(sdi_level = sdi_level2021) %>%
    left_join(
        data_combined %>%
            filter(year == 2021) %>%
            select(location_id, sdi_level2021),
        by = c("location_id")
    ) %>%
    mutate(shortage = ifelse(is.na(shortage), 2, shortage)) %>%
    group_by(year, sdi_level2021, shortage) %>%
    summarise(
        total_number = sum(val, na.rm = TRUE),
        total_population = sum(population, na.rm = TRUE),
        proportion = total_number / total_population
    ) %>%
    ungroup()

data_syphilis_summary <-
    data_syphilis_summary %>%
    mutate(
        sdi_level2021 = factor(
            sdi_level2021,
            levels = c("Low", "Low-middle", "Middle", "High-middle", "High")
        )
    )
View(data_syphilis_summary)

incidence_plot <- ggplot(
    data = data_syphilis_summary,
    mapping = aes(x = year, y = proportion * 1000, colour = factor(shortage))
) +
    geom_line() +
    geom_point(size = 2) +
    geom_vline(xintercept = c(2014, 2016), linetype = "dashed", color = "grey") +
    labs(
        title = "Congenital Syphilis Incidence Rate across SDI Quintiles in 2021",
        x = "Year",
        y = "Incidence rate (cases per 1000 births)"
    ) +
    theme_jama() +
    scale_colour_manual(
        name = "Shortage between 2014 and 2016",
        values = c("blue", "red", "grey"),
        labels = c("No Shortage", "Shortage", "Unknown")
    ) +
    scale_y_continuous(breaks = seq(0, 10, by = 1), limits = c(0, 10)) +
    scale_x_continuous(breaks = seq(1990, 2021, by = 1), labels = seq(1990, 2021, by = 1)) +
    facet_wrap(~sdi_level2021)
plot(incidence_plot)
ggsave("syphilis_congenital_incidence_by_SDI_plot.png",
    plot = incidence_plot, width = 12, height = 6, dpi = 300
)

# check the data
data_syphilis_summary <-
    data_combined %>%
    mutate(shortage = ifelse(is.na(shortage), 2, shortage)) %>%
    rename(sdi_level = sdi_level2021) %>%
    left_join(
        data_combined %>%
            filter(year == 2021) %>%
            select(location_id, sdi_level2021),
        by = c("location_id")
    ) %>%
    filter(sdi_level2021 == "Middle" & shortage == 1) %>%
    group_by(year, location_name, shortage) %>%
    summarise(
        total_number = sum(val, na.rm = TRUE),
        total_population = sum(population, na.rm = TRUE),
        proportion = total_number / total_population
    ) %>%
    ungroup()
incidence_plot <- ggplot(
    data = data_syphilis_summary,
    mapping = aes(x = year, y = proportion * 1000, colour = factor(shortage))
) +
    geom_line() +
    geom_point(size = 2) +
    geom_vline(xintercept = c(2014, 2016), linetype = "dashed", color = "grey") +
    labs(
        title = "Congenital Syphilis Incidence Rate across SDI Quintiles in 2021",
        x = "Year",
        y = "Incidence rate (cases per 1000 births)"
    ) +
    theme_jama() +
    scale_colour_manual(
        name = "Shortage between 2014 and 2016",
        values = c("blue", "red", "grey"),
        labels = c("No Shortage", "Shortage", "Unknown")
    ) +
    # scale_y_continuous(breaks = seq(0, 10, by = 1), limits = c(0, 10)) +
    scale_x_continuous(breaks = seq(1990, 2021, by = 1), labels = seq(1990, 2021, by = 1)) +
    facet_wrap(~location_name)
plot(incidence_plot)

#### HDI levels in 2021 ------------------------------------------------------
data_syphilis_summary <-
    data_combined %>%
    mutate(shortage = ifelse(is.na(shortage), 2, shortage)) %>%
    left_join(
        data_combined %>%
            filter(year == 2021) %>%
            select(location_id, hdi_level) %>%
            rename(hdi_level2021 = hdi_level),
        by = c("location_id")
    ) %>%
    group_by(year, hdi_level2021, shortage) %>%
    summarise(
        total_number = sum(val, na.rm = TRUE),
        total_population = sum(population, na.rm = TRUE),
        proportion = total_number / total_population
    ) %>%
    ungroup()

unique(data_syphilis_summary$hdi_level2021)
unique(data_combined$location_name[is.na(data_combined$hdi_level) & data_combined$year == 2021])

incidence_plot <- ggplot(
    data = data_syphilis_summary %>% filter(!is.na(hdi_level2021)),
    mapping = aes(x = year, y = proportion * 1000, colour = factor(shortage))
) +
    geom_line() +
    geom_point(size = 2) +
    geom_vline(xintercept = c(2014, 2016), linetype = "dashed", color = "grey") +
    labs(
        title = "Congenital Syphilis Incidence Rate across HDI Levels in 2021",
        x = "Year",
        y = "Incidence rate (cases per 1000 births)"
    ) +
    theme_jama() +
    scale_colour_manual(
        name = "Shortage between 2014 and 2016",
        values = c("blue", "red", "grey"),
        labels = c("No Shortage", "Shortage", "Unknown")
    ) +
    scale_y_continuous(breaks = seq(0, 10, by = 1), limits = c(0, 10)) +
    # scale_y_continuous(labels = scales::percent_format(scale = 1)) +
    scale_x_continuous(breaks = seq(1990, 2021, by = 1), labels = seq(1990, 2021, by = 1)) +
    facet_wrap(~hdi_level2021, scales = "free")
incidence_plot
ggsave("syphilis_congenital_incidence_by_HDI_plot.png",
    plot = incidence_plot, width = 12, height = 6, dpi = 300
)

#### income levels in 2021 ---------------------------------------------------
data_syphilis_summary <-
    data_combined %>%
    mutate(shortage = ifelse(is.na(shortage), 2, shortage)) %>%
    left_join(
        data_combined %>%
            filter(year == 2021) %>%
            select(location_id, income_level) %>%
            rename(income_level2021 = income_level),
        by = c("location_id")
    ) %>%
    group_by(year, income_level2021, shortage) %>%
    summarise(
        total_number = sum(val, na.rm = TRUE),
        total_population = sum(population, na.rm = TRUE),
        proportion = total_number / total_population
    ) %>%
    ungroup()

unique(data_syphilis_summary$income_level2021)

incidence_plot <- ggplot(
    data = data_syphilis_summary %>% filter(!is.na(income_level2021)),
    mapping = aes(x = year, y = proportion * 1000, colour = factor(shortage))
) +
    geom_line() +
    geom_point(size = 2) +
    geom_vline(xintercept = c(2014, 2016), linetype = "dashed", color = "grey") +
    labs(
        title = "Congenital Syphilis Incidence Rate across Income Levels in 2021",
        x = "Year",
        y = "Incidence rate (cases per 1000 births)"
    ) +
    theme_jama() +
    scale_colour_manual(
        name = "Shortage between 2014 and 2016",
        values = c("blue", "red", "grey"),
        labels = c("No Shortage", "Shortage", "Unknown")
    ) +
    scale_y_continuous(breaks = seq(0, 10, by = 1), limits = c(0, 10)) +
    # scale_y_continuous(labels = scales::percent_format(scale = 1)) +
    scale_x_continuous(breaks = seq(1990, 2021, by = 1), labels = seq(1990, 2021, by = 1)) +
    facet_wrap(~income_level2021, scales = "free")
incidence_plot
ggsave("syphilis_congenital_incidence_by_income_plot.png",
    plot = incidence_plot, width = 12, height = 6, dpi = 300
)

### Country map --------------------------------------------------------------
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)

data_country_summary <- data_combined %>%
    mutate(shortage = ifelse(is.na(shortage), 2, shortage)) %>%
    group_by(year, shortage, location_iso) %>%
    summarise(
        total_number = sum(val, na.rm = TRUE),
        total_population = mean(population, na.rm = TRUE),
        proportion = total_number / total_population * 1000
    ) %>%
    ungroup()
View(data_country_summary)

data_country_map <- world %>%
    left_join(
        data_country_summary %>%
            filter(year == 2021) %>%
            select(location_iso, proportion),
        by = c("iso_a3_eh" = "location_iso")
    )
data_country_map <-
    data_country_map %>%
    mutate(proportion2 = cut(proportion,
        breaks = c(0, 0.5, 1, 5, 10, Inf),
        labels = c(
            "0~0.5", "0.5~1.0",
            "1.0~5.0", "5.0~10.0", "10.0+"
        ), include.lowest = TRUE, right = TRUE
    ))

incidence_rate_plot2021 <- ggplot(data = data_country_map) +
    geom_sf(aes(fill = proportion2)) +
    scale_fill_brewer(palette = "Spectral", direction = -1, name = "Incidence rate") +
    labs(title = "(B) Congenital Syphilis Incidence Rate (Cases per 1000 births) in 2021") +
    theme_void() +
    theme(
        # Base elements
        text = element_text(color = "black", size = 10),

        # Plot title/subtitle
        plot.title = element_text(face = "bold", size = 12, hjust = 0),
        plot.subtitle = element_text(face = "italic", hjust = 0),

        # Legend (top position to save space)
        legend.position = "inside",
        legend.position.inside = c(0.08, 0.4),
        legend.direction = "vertical",
        legend.title = element_text(face = "bold")
    )

data_country_map2 <- world %>%
    left_join(
        data_country_summary %>%
            filter(year == 2010) %>%
            select(location_iso, proportion),
        by = c("iso_a3_eh" = "location_iso")
    )
data_country_map2 <-
    data_country_map2 %>%
    mutate(proportion2 = cut(proportion,
        breaks = c(0, 0.5, 1, 5, 10, Inf),
        labels = c(
            "0~0.5", "0.5~1.0",
            "1.0~5.0", "5.0~10.0", "10.0+"
        ), include.lowest = TRUE, right = TRUE
    ))

incidence_rate_plot2010 <- ggplot(data = data_country_map2) +
    geom_sf(aes(fill = proportion2)) +
    scale_fill_brewer(palette = "Spectral", direction = -1, name = "Incidence rate") +
    labs(title = "(A) Congenital Syphilis Incidence Rate (Cases per 1000 births) in 2010") +
    theme_void() +
    theme(
        # Base elements
        text = element_text(color = "black", size = 10),

        # Plot title/subtitle
        plot.title = element_text(face = "bold", size = 12, hjust = 0),
        plot.subtitle = element_text(face = "italic", hjust = 0),

        # Legend (top position to save space)
        legend.position = "inside",
        legend.position.inside = c(0.08, 0.4),
        legend.direction = "vertical",
        legend.title = element_text(face = "bold")
    )
incidence_rate_plot2010

incidence_rate_plot2010_2021 <-
    ggpubr::ggarrange(
        incidence_rate_plot2010,
        incidence_rate_plot2021,
        ncol = 1, nrow = 2
    )

incidence_rate_plot2010_2021
ggsave("syphilis_congenital_incidence_rate_map_2010_2021.png",
    plot = incidence_rate_plot2010_2021,
    width = 20, height = 12, dpi = 300
)

#### Change between 2010 and 2021 --------------------------------------------
data_country_change <- data_country_summary %>%
    filter(year %in% c(2010, 2021)) %>%
    select(location_iso, shortage, year, proportion) %>%
    tidyr::pivot_wider(names_from = year, values_from = proportion, names_prefix = "year_") %>%
    mutate(
        change = (year_2021 - year_2010) / year_2010 * 100
    )

data_country_change %>% 
    filter(!is.na(shortage)) %>%
    mutate(change_cat = cut(change,
        breaks = c(-Inf, -25, 0, 25, Inf),
        include.lowest = TRUE, right = TRUE
    )) %>%
    group_by(shortage) %>%
    mutate(ntotal = n()) %>%
    group_by(shortage, change_cat) %>%
    summarise(
        n = n(),
        ntotal = first(ntotal),
        prop = mean(n/ntotal * 100),
        mean_change = mean(change, na.rm = TRUE),
        median_change = median(change, na.rm = TRUE),
        min_change = min(change, na.rm = TRUE),
        max_change = max(change, na.rm = TRUE),
        sd_change = sd(change, na.rm = TRUE),
        q25_change = quantile(change, 0.25, na.rm = TRUE),
        q75_change = quantile(change, 0.75, na.rm = TRUE),
        iqr_change = IQR(change, na.rm = TRUE)
    )

data_country_map3 <- world %>%
    left_join(
        data_country_change %>%
            select(location_iso, shortage, change),
        by = c("iso_a3_eh" = "location_iso")
    )
summary(data_country_map3$change)

data_country_map3 <-
    data_country_map3 %>%
    mutate(change2 = cut(change,
        breaks = c(-Inf, -25, 0, 25, 50, Inf),
        include.lowest = TRUE, right = TRUE
    ))
levels <- levels(data_country_map3$change2)

data_country_change %>%
    group_by(shortage) %>%
    summarise(
        n = n(),
        mean_change = mean(change, na.rm = TRUE),
        median_change = median(change, na.rm = TRUE),
        min_change = min(change, na.rm = TRUE),
        max_change = max(change, na.rm = TRUE)
    )

incidence_rate_change_plot <- ggplot(data = data_country_map3) +
    geom_sf(aes(fill = change2)) +
    scale_fill_brewer(
        palette = "Spectral", direction = -1,
        name = "Percent change (%)"
    ) +
    labs(title = "(A) Congenital Syphilis Incidence Rate Change without Shortage reports") +
    theme_void() +
    theme(
        # Base elements
        text = element_text(color = "black", size = 10),

        # Plot title/subtitle
        plot.title = element_text(face = "bold", size = 12, hjust = 0),
        plot.subtitle = element_text(face = "italic", hjust = 0),

        # Legend (top position to save space)
        legend.position = "inside",
        legend.position.inside = c(0.08, 0.4),
        legend.direction = "vertical",
        legend.title = element_text(face = "bold")
    )
incidence_rate_change_plot
ggsave("syphilis_congenital_incidence_rate_change.png",
    plot = incidence_rate_change_plot,
    width = 12, height = 6, dpi = 300
)

library(RColorBrewer)
brewer.pal(n = 7, name = "Spectral")

incidence_rate_change_noshortage_plot <- ggplot(data = data_country_map3 %>%
    mutate(
        change2 = ifelse(shortage == 0, as.character(change2), NA),
        change2 = factor(change2, levels = levels)
    )) +
    geom_sf(aes(fill = change2)) +
    scale_fill_brewer(palette = "Spectral", direction = -1, name = "Incidence rate change (%)") +
    labs(title = "(A) Congenital Syphilis Incidence Rate Change without Shortage reports") +
    theme_void() +
    theme(
        # Base elements
        text = element_text(color = "black", size = 10),

        # Plot title/subtitle
        plot.title = element_text(face = "bold", size = 12, hjust = 0),
        plot.subtitle = element_text(face = "italic", hjust = 0),

        # Legend (top position to save space)
        legend.position = "inside",
        legend.position.inside = c(0.08, 0.4),
        legend.direction = "vertical",
        legend.title = element_text(face = "bold")
    )
incidence_rate_change_noshortage_plot

incidence_rate_change_shortage_plot <- ggplot(
    data = data_country_map3 %>%
        mutate(
            change2 = ifelse(shortage == 1, as.character(change2), NA),
            change2 = factor(change2, levels = levels)
        )
) +
    geom_sf(aes(fill = change2)) +
    scale_fill_brewer(palette = "Spectral", direction = -1, name = "Incidence rate change (%)") +
    labs(title = "(B) Congenital Syphilis Incidence Rate Change with Shortage reports") +
    theme_void() +
    theme(
        # Base elements
        text = element_text(color = "black", size = 10),

        # Plot title/subtitle
        plot.title = element_text(face = "bold", size = 12, hjust = 0),
        plot.subtitle = element_text(face = "italic", hjust = 0),

        # Legend (top position to save space)
        legend.position = "inside",
        legend.position.inside = c(0.08, 0.4),
        legend.direction = "vertical",
        legend.title = element_text(face = "bold")
    )
incidence_rate_change_shortage_plot

incidence_rate_change_by_shortage_plot <-
    ggpubr::ggarrange(
        incidence_rate_change_noshortage_plot,
        incidence_rate_change_shortage_plot,
        ncol = 1, nrow = 2
    )

incidence_rate_change_by_shortage_plot
ggsave("syphilis_congenital_incidence_rate_change_by_shortage.png",
    plot = incidence_rate_change_by_shortage_plot,
    width = 12, height = 12, dpi = 300
)

#### For each country ####
plot_its <- function(
    data,
    location_name_para,
    measure_name_para,
    shortage_para = "shortage"
) {
    message(paste(measure_name_para, location_name_para))
    data$shortage <- data[[shortage_para]]

    data_temp <- data %>%
        mutate(
            shortage = factor(ifelse(is.na(shortage), 2, shortage),
                levels = c(0, 1, 2)
            )
        ) %>%
        filter(location_name == location_name_para)

    data_temp_summary <- data_temp %>%
        group_by(year, shortage) %>%
        summarise(
            total_number = sum(val, na.rm = TRUE),
            total_population = mean(population, na.rm = TRUE),
            proportion = total_number / total_population
        ) %>%
        ungroup()

    title <- paste0(
        "Syphilis ",
        measure_name_para,
        " Proportion of Births Over Years in ",
        location_name_para
    )

    filename <- paste0(
        "Results/Figures/syphilis_congenital_",
        tolower(gsub(" ", "_", measure_name_para)), "_",
        tolower(gsub(" ", "_", location_name_para)), "_plot.png"
    )

    shortage_colour <- ifelse(
        all(data_temp$shortage == 1),
        ggsci::pal_jama("default")(3)[1],
        ifelse(
            any(data_temp$shortage == 2),
            ggsci::pal_jama("default")(3)[3],
            ggsci::pal_jama("default")(3)[2]
        )
    )

    shortage_label <- ifelse(
        all(data_temp$shortage == 1),
        "Shortage",
        ifelse(
            any(data_temp$shortage == 2),
            "Unknown",
            "No shortage"
        )
    )

    temp_plot <- ggplot(
        data = data_temp_summary,
        mapping = aes(x = year, y = proportion * 100, colour = shortage_colour)
    ) +
        geom_line() +
        geom_point(size = 4) +
        geom_vline(xintercept = c(2014, 2016), linetype = "dashed", color = "grey") +
        labs(
            title = title,
            x = "Year",
            y = paste(measure_name_para, " Proportion (%)")
        ) +
        theme_jama() +
        scale_colour_manual(
            name = "Shortage between 2014 and 2016",
            values = shortage_colour,
            labels = shortage_label
        ) +
        # scale_y_continuous(labels = scales::percent_format(scale = 1)) +
        scale_x_continuous(breaks = seq(1990, 2021, by = 5), labels = seq(1990, 2021, by = 5))

    message(paste("Plot saved as", filename))
    ggsave(filename, temp_plot, width = 10, height = 6, dpi = 300)

    # Fit a model to the incidence data
    # We will use a linear model with an interaction term for the time and interruption variable
    data_temp_itsa <- data_temp %>%
        filter(year >= 2010 & year <= 2021) %>%
        mutate(
            interrupt_var = ifelse(year >= 2016, 1, 0),
            proportion = val / population * 100
        ) %>%
        as.data.frame()

    itsa_model <- nlme::gls(
        proportion ~ year + interrupt_var + year:interrupt_var,
        data = data_temp_itsa
    )

    ## generate the 95% CI
    itsa_model_confint <- confint(itsa_model)
    itsa_results <- cbind(
        data.frame(coef = itsa_model$coefficients),
        itsa_model_confint
    )
    itsa_results$shortage_status <- shortage_label
    return(itsa_results)
}

results_incidence_list <- lapply(
    unique(data_combined$location_name),
    function(location) {
        results <- plot_its(
            data = data_combined,
            location_name_para = location,
            measure_name_para = "Incidence"
        )
        results$location_name <- location
        return(results)
    }
)
results_incidence_list <- do.call(rbind, results_incidence_list)
write.xlsx(results_incidence_list,
    file = "Results/syphilis_congenital_incidence_results.xlsx", rowNames = TRUE
)

# 4. Baseline characteristics tables -----------------------------------------
# Summary table by shortage groups
library(tableone)
colnames(data_combined)
x_vars <- c(
    "hdi", "hdi_level", "sdi", "sdi_level2021",
    "gni", "income_level", "pregnant_antenatal_care",
    "syphilis_prev_female", "syphilis_prev_male", "syphilis_prev_pregnant",
    "abor_misc_rate", "gbd_superregion"
)
x_factorVars <- c("hdi_level", "sdi_level2021", "income_level", "gbd_superregion")
x_nonnormalVars <- c(
    "hdi", "sdi", "gni", "pregnant_antenatal_care",
    "syphilis_prev_female", "syphilis_prev_male", "syphilis_prev_pregnant",
    "abor_misc_rate"
)

# Function to calculate missing data by variable
na_by_var <- function(data, vars) {
    data %>%
        summarise(across(all_of(vars), ~ sum(is.na(.)))) %>%
        tidyr::pivot_longer(everything(),
            names_to = "variable",
            values_to = "n_missing"
        ) %>%
        mutate(
            n_rows = nrow(data),
            pct_missing = 100 * n_missing / n_rows
        ) %>%
        arrange(desc(pct_missing))
}

# Calculate missing data for 2010 and 2021
na_by_var(
    data = data_combined %>%
        mutate(shortage = ifelse(is.na(shortage), 2, shortage)) %>%
        filter(year == 2010 & shortage != 2),
    vars = x_vars
)

na_by_var(
    data = data_combined %>%
        mutate(shortage = ifelse(is.na(shortage), 2, shortage)) %>%
        filter(year == 2021 & shortage != 2),
    vars = x_vars
)

# Create TableOne for 2010
x_vars <- c(
    "gni", "income_level", "hdi", "hdi_level", "sdi", "sdi_level2021", "gbd_superregion",
    "syphilis_prev_female", "syphilis_prev_male", "abor_misc_rate",
    "syphilis_prev_pregnant", "pregnant_antenatal_care"
)
x_factorVars <- c("hdi_level", "sdi_level2021", "income_level", "gbd_superregion")
x_nonnormalVars <- c(
    "hdi", "sdi", "gni", "pregnant_antenatal_care",
    "syphilis_prev_female", "syphilis_prev_male", "syphilis_prev_pregnant",
    "abor_misc_rate"
)

tb1 <- tableone::CreateTableOne(
    data = data_combined %>%
        mutate(shortage = ifelse(is.na(shortage), 2, shortage)) %>%
        filter(year == 2021 & shortage != 2),
    vars = x_vars,
    factorVars = x_factorVars,
    strata = "shortage",
    addOverall = TRUE,
    includeNA = TRUE
)
tb1 <- print(
    tb1,
    contDigits = 2,
    nonnormal = x_nonnormalVars,
    exact = x_factorVars,
    smd = TRUE,
    missing = TRUE,
    quote = TRUE,
    noSpaces = TRUE,
    pDigits = 2,
)
tb1

write.xlsx(tb1,
    file = "Results/syphilis_congenital_baseline_characteristics_table.xlsx",
    rowNames = TRUE
)
