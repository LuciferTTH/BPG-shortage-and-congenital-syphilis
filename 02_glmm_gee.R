# File: 02_glmm_gee.R
# Purpose: Mixed (GLMM) and GEE analyses for syphilis outcomes
# Author: Tao Huang (restructured)
# Date: 2025-11-26
# Notes: tidy steps, save outputs to Results/


# 00 Setup ####
required_pkgs <- c(
    "openxlsx", "dplyr", "ggplot2", "httpgd", "showtext", "ggsci",
    "scales", "survey", "its.analysis", "lmtest", "sandwich", "lme4",
    "lmerTest", "jtools", "geepack", "glmmTMB", "performance", "nlme"
)
missing_pkgs <- required_pkgs[
    !sapply(required_pkgs, requireNamespace, quietly = TRUE)
]
if (length(missing_pkgs)) {
    stop(
        "Install missing packages: ",
        paste(missing_pkgs, collapse = ", ")
    )
}
lapply(required_pkgs, library, character.only = TRUE)

if (!dir.exists("Results")) {
    dir.create("Results", recursive = TRUE)
}

## Define theme for ggplot2 ####
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

# 01 Load data ####
load("data_combined.Rdata")

data_combined %>%
    filter(year == 2021) %>%
    filter(is.na(income_level)) %>%
    select(location_id, location_name)

data_combined %>%
    filter(year == 2021) %>%
    filter(is.na(hdi_level)) %>%
    select(location_name)


# 02 Prepare analysis dataset ####
data_analysis <- data_combined %>%
    filter(!is.na(shortage)) %>%
    mutate(
        location_id = factor(location_id),
        time = year,
        time_c = year - 2016,
        rate = (val / population) * 1000, # per 1000 births
        post = as.integer(year >= 2016),
        weights = population / sum(population)
    )
colnames(data_analysis)

data_combined %>%
    mutate(
        location_id = factor(location_id),
        time = year,
        time_c = year - 2016,
        rate = (val / population) * 1000, # per 1000 births
        post = as.integer(year >= 2016),
        weights = population / sum(population)
    ) %>%
    group_by(time) %>%
    summarise(
        n = n(),
        sum_weights = sum(weights, na.rm = TRUE),
        mean_rate = weighted.mean(rate, w = population, na.rm = TRUE),
        .groups = "drop"
    )

data_combined %>%
    filter(!is.na(shortage)) %>%
    mutate(
        location_id = factor(location_id),
        time = year,
        time_c = year - 2016,
        rate = (val / population) * 1000, # per 1000 births
        post = as.integer(year >= 2016),
        weights = population / sum(population)
    ) %>%
    filter(shortage == 1) %>%
    group_by(time, shortage) %>%
    summarise(
        n = n(),
        sum_weights = sum(weights, na.rm = TRUE),
        mean_rate = weighted.mean(rate, w = population, na.rm = TRUE),
        .groups = "drop"
    )


# volint plot with smoothed line
ggplot(
    data_analysis,
    aes(x = factor(time), y = rate)
) +
    geom_violin(
        aes(weight = population, fill = as.factor(shortage)),
        position = position_dodge(width = 0.1),
        alpha = 0.5,
        width = 0.5
    ) +
    geom_smooth(
        data = data_analysis,
        aes(
            x = as.numeric(factor(time)), y = rate,
            weight = population,
            colour = as.factor(shortage)
        ),
        method = "loess"
    ) +
    labs(
        title = "Violin plot with smoothed line",
        x = "Year",
        y = "Incidence rate (per 1,000 births)"
    ) +
    scale_color_discrete(
        name = "With shortage reports",
        breaks = c("0", "1"),
        palette = ggsci::pal_lancet()(2),
        labels = c("No", "Yes")
    ) +
    scale_fill_discrete(
        name = "With shortage reports",
        breaks = c("0", "1"),
        palette = ggsci::pal_lancet()(2),
        labels = c("No", "Yes")
    ) +
    scale_x_discrete(
        breaks = seq(2010, 2021, by = 1)
    ) +
    # scale_y_continuous(breaks = seq(0, 1, by = 1)) +
    facet_wrap(
        ~shortage,
        labeller = as_labeller(c(
            "0" = "Countries without shortage reports",
            "1" = "Countries with shortage reports"
        ))
    ) +
    theme_jama()
ggsave("Results/violin_plot_smoothed_line.png", width = 9, height = 6)

# Box plot with smoothed line
ggplot(
    data_analysis,
    aes(x = factor(time), y = rate)
) +
    geom_boxplot(
        aes(
            weight = weights,
            fill = as.factor(shortage)
        ),
        position = position_dodge(width = 0.1),
        width = 0.5,
        outlier.shape = NA,
        alpha = 0.5
    ) +
    geom_smooth(
        data = data_analysis,
        aes(
            x = as.numeric(factor(time)), y = rate,
            weight = weights,
            colour = as.factor(shortage)
        ),
        method = "loess"
    ) +
    labs(
        title = "Box plot with smoothed line",
        x = "Year",
        y = "Incidence rate (per 1,000 births)"
    ) +
    scale_color_discrete(
        name = "With shortage reports",
        breaks = c("0", "1"),
        palette = ggsci::pal_lancet()(2),
        labels = c("No", "Yes")
    ) +
    scale_fill_discrete(
        name = "With shortage reports",
        breaks = c("0", "1"),
        palette = ggsci::pal_lancet()(2),
        labels = c("No", "Yes")
    ) +
    scale_x_discrete(
        breaks = seq(2010, 2021, by = 1)
    ) +
    # scale_y_continuous(breaks = seq(0, 1, by = 1)) +
    facet_wrap(
        ~shortage,
        labeller = as_labeller(c(
            "0" = "Countries without shortage reports",
            "1" = "Countries with shortage reports"
        ))
    ) +
    theme_jama()
ggsave("Results/box_plot_smoothed_line.png", width = 9, height = 6)

# scatter plot with smoothed line
ggplot(
    data_analysis,
    aes(x = time, y = rate, colour = as.factor(shortage))
) +
    geom_jitter(width = 0.1, alpha = 0.3) +
    geom_smooth(
        aes(weight = population),
        method = "loess"
    ) +
    labs(
        title = "Scatter plot with smoothed line",
        x = "Year",
        y = "Incidence rate (per 1,000 births)"
    ) +
    scale_color_discrete(
        name = "With shortage reports",
        breaks = c("0", "1"),
        palette = ggsci::pal_lancet()(2),
        labels = c("No", "Yes")
    ) +
    scale_y_continuous(breaks = seq(0, 15, by = 1)) +
    scale_x_continuous(breaks = seq(2010, 2021, by = 1)) +
    facet_wrap(
        ~shortage,
        labeller = as_labeller(c(
            "0" = "Countries without shortage reports",
            "1" = "Countries with shortage reports"
        ))
    ) +
    theme_jama()
ggsave("Results/scatter_plot_smoothed_line.png", width = 9, height = 6)


# 03 Mixed effects models ####
##  03.1 Linear mixed model (rate) ####
data_analysis <- data_combined %>%
    filter(!is.na(shortage)) %>%
    mutate(
        location_id = factor(location_id),
        time = year,
        time_c = year - 2016,
        rate = (val / population) * 1000, # per 1000 births
        post = as.integer(year >= 2016),
        weights = population / sum(population)
    )
lmm_formula <- rate ~ time_c * shortage * post +
    hdi + sdi + gni + gbd_superregion +
    syphilis_prev_female + syphilis_prev_male + abor_misc_rate +
    (1 + time_c | location_id)
stepwise_lmm <- step(
    lmer(
        lmm_formula,
        data = data_analysis,
        weights = weights,
        REML = FALSE
    ),
    reduce.fixed = FALSE
)
stepwise_lmm

# get final model
stepwise_lmm <- get_model(stepwise_lmm)
stepwise_lmm
summary(stepwise_lmm)

# diagnostics
glmm_modelcheck <- performance::check_model(stepwise_lmm)
glmm_modelcheck
png(
    file = "Results/lmm_modelcheck1219.png",
    width = 1200, height = 2400, res = 100
)
glmm_modelcheck
dev.off()

# summary
lmm_sum <- summary(stepwise_lmm)
lmm_sum
fixed_lmm <- as.data.frame(lmm_sum$coefficients)
fixed_lmm
fixed_lmm_results <- data.frame(
    Parameter = rownames(fixed_lmm),
    Estimate = fixed_lmm[, "Estimate"],
    SE = fixed_lmm[, "Std. Error"],
    t_value = fixed_lmm[, "t value"], # t-tests use Satterthwaite's method
    p_value = fixed_lmm[, "Pr(>|t|)"]
) %>%
    mutate(
        CI_lower = Estimate - 1.96 * SE,
        CI_upper = Estimate + 1.96 * SE
    )
fixed_lmm_results
write.xlsx(
    fixed_lmm_results,
    "Results/linear_mixed_effects_results.xlsx"
)

jtools::plot_summs(stepwise_lmm,
    omit.coefs = c("(Intercept)", "hdi", "sdi", "gni", "gbd_superregion*"),
    model.names = "LMM (rate)",
    robust = TRUE
) +
    labs(
        title = "Linear Mixed Model: Estimates (95% CI)",
        x = "Estimate",
        y = "Parameter"
    )
ggsave("Results/lmm_coefplot.png", width = 9, height = 6)

## 03.2 Counterfactual estimates from LMM (Conditional Effects) ####
model_cf <- stepwise_lmm

years <- 2010:2021

# 1) Observed: shortage group (post as observed)
scenario_obs_shortage <- data_analysis %>%
  filter(shortage == 1) %>%
  mutate(scenario = "Observed: shortage")

# 2) Observed: non-shortage group (post as observed)
scenario_obs_noshortage <- data_analysis %>%
  filter(shortage == 0) %>%
  mutate(scenario = "Observed: non-shortage")

# 3) Counterfactual: shortage group with post = 0 for all years
scenario_cf_shortage <- data_analysis %>%
  filter(shortage == 1) %>%
  mutate(post = 0, scenario = "Counterfactual post=0: shortage")

# 4) Counterfactual: non-shortage group with post = 0 for all years
scenario_cf_noshortage <- data_analysis %>%
  filter(shortage == 0) %>%
  mutate(post = 0, scenario = "Counterfactual post=0: non-shortage")

# Combine all scenarios
all_scenarios <- bind_rows(
  scenario_obs_shortage,
  scenario_obs_noshortage,
  scenario_cf_shortage,
  scenario_cf_noshortage
)

# Make conditional predictions with standard errors
pred_obj <- predict(model_cf, newdata = all_scenarios, re.form = NULL, 
                    allow.new.levels = TRUE, se.fit = TRUE)

all_scenarios$estimate <- pred_obj$fit
all_scenarios$se_fit <- pred_obj$se.fit

# Aggregate by year and scenario
pred_4scenario_conditional <- all_scenarios %>%
  group_by(year, scenario) %>%
  summarise(
    n = sum(!is.na(estimate)),
    sum_weights = sum(population, na.rm = TRUE),
    # Population-weighted mean of predictions
    wm_estimate = weighted.mean(estimate, w = population, na.rm = TRUE),
    # Correct SE: population-weighted variance of model SEs
    # SE of weighted mean = sqrt(sum of weighted variances)
    se = sqrt(sum(population / sum(population, na.rm = TRUE) * se_fit^2, na.rm = TRUE)) ,
    .groups = "drop"
  ) %>%
  rename(estimate = wm_estimate) %>%
  mutate(
    ci_lower = estimate - 1.96 * se,
    ci_upper = estimate + 1.96 * se
  ) %>%
  select(year, scenario, estimate, se, ci_lower, ci_upper)

write.csv(pred_4scenario_conditional, "Results/lmm_4scenario_predictions_conditional.csv", row.names = FALSE)

# Calculate observed rates by year and shortage status
observed_rates <- data_analysis %>%
  group_by(year, shortage) %>%
  summarise(
    observed_rate = weighted.mean(rate, w = population, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Group = ifelse(shortage == 0, "Non-shortage", "Shortage"),
    Group = factor(Group, levels = c("Non-shortage", "Shortage"))
  )

# Prepare plot data with model predictions
plot_dat <- pred_4scenario_conditional %>%
  mutate(
    Group = ifelse(grepl("non-shortage", scenario), "Non-shortage", "Shortage"),
    Type = ifelse(grepl("^Observed", scenario), "Observed", "Counterfactual")
  ) %>%
  mutate(
    Group = factor(Group, levels = c("Non-shortage", "Shortage")),
    Type = factor(Type, levels = c("Observed", "Counterfactual"))
  ) %>%
  select(year, Group, Type, estimate, ci_lower, ci_upper)

plot_dat

# Create the plot with confidence bands and observed data
p4 <- ggplot() +
  geom_vline(xintercept = 2016, linetype = "dashed", color = "gray50") +
  geom_text(
    data = data.frame(x = 2016.2, y = 4, label = "Intervention Year = 2016"),
    aes(x = x, y = y, label = label),
    vjust = 0, hjust = 0,
    color = "gray50", size = 3
  ) +
  # Observed data points
  geom_point(
    data = observed_rates,
    aes(x = year, y = observed_rate, color = Group, shape = Group),
    size = 3, alpha = 0.6
  ) +
  # Model predictions with confidence bands
  geom_line(
    data = plot_dat,
    aes(x = year, y = estimate, color = Group, linetype = Type),
    linewidth = 1.1
  ) +
  # Confidence interval ribbons
  geom_ribbon(
    data = plot_dat %>% filter(year >= 2016 & Type == "Counterfactual"),  # Only add ribbons for counterfactual scenarios
    aes(x = year, ymin = ci_lower, ymax = ci_upper, fill = Group, linetype = Type),
    alpha = 0.2,
    color = NA
  ) +
  geom_ribbon(
    data = plot_dat %>% filter(Type == "Observed"),  # Only add ribbons for observed scenarios
    aes(x = year, ymin = ci_lower, ymax = ci_upper, fill = Group, linetype = Type),
    alpha = 0.2,
    color = NA
  ) +
  scale_color_lancet(name = "BPG Shortage") +
  scale_fill_lancet(name = "BPG Shortage") +
  scale_linetype_manual(
    name = "Scenario",
    values = c("Observed" = "solid", "Counterfactual" = "dashed")
  ) +
  scale_shape_manual(
    name = "BPG Shortage",
    values = c("Non-shortage" = 16, "Shortage" = 17)
  ) +
  scale_x_continuous(breaks = 2010:2021) +
  scale_y_continuous(breaks = seq(0, 5, by = 0.5)) +
  coord_cartesian(ylim = c(0, NA)) +
  labs(
    x = "Year",
    y = "Incidence rate per 1,000 births"
  ) +
  theme_jama() +
  theme(
    legend.position = "right",
    legend.key.width = unit(1.5, "cm")  # Increase legend key width
  ) +
  guides(linetype = guide_legend(override.aes = list(linewidth = 1.2)))
p4
ggsave("Results/lmm_observed_vs_counterfactual_4lines_conditional.png", p4, width = 12, height = 7, dpi = 300)

## 03.2 Poisson GLMM: counts with offset(log(population)) ####
glmm <- glmmTMB(
    val ~ time_c * shortage * post +
        hdi + sdi + gni + gbd_superregion +
        syphilis_prev_female + syphilis_prev_male + abor_misc_rate +
        (1 + time_c | location_id),
    data = data_analysis,
    offset = log(population),
    family = poisson
)
summary(glmm)
glmm

# diagnostics
glmm_modelcheck <- performance::check_model(glmm)
png(
    file = "Results/mixed_glmm_poisson_modelcheck.png",
    width = 1200, height = 2400, res = 100
)
glmm_modelcheck
dev.off()

glmm_sum <- summary(glmm)
fixed_glmm <- as.data.frame(glmm_sum$coefficients$cond)
fixed_glmm

fixed_glmm_results <- data.frame(
    Parameter = rownames(fixed_glmm),
    Estimate = fixed_glmm[, "Estimate"],
    SE = fixed_glmm[, "Std. Error"],
    z_value = fixed_glmm[, "z value"],
    p_value = fixed_glmm[, "Pr(>|z|)"]
) %>%
    mutate(
        CI_lower = Estimate - 1.96 * SE,
        CI_upper = Estimate + 1.96 * SE
    )
fixed_glmm_results

exp_tab <- fixed_glmm_results %>%
    mutate(
        Exp = exp(Estimate),
        Exp_CI_low = exp(CI_lower),
        Exp_CI_up = exp(CI_upper)
    )
exp_tab
write.xlsx(exp_tab, "Results/mixed_glmm_poisson_exp_estimates2016.xlsx")

jtools::plot_summs(glmm,
    exp = TRUE,
    omit.coefs = c("(Intercept)"),
    model.names = "GLMM Poisson"
) +
    labs(
        title = "Poisson Mixed Model: Exponentiated Estimates",
        x = "Estimate (exp)",
        y = "Parameter"
    )
ggsave("Results/mixed_glmm_poisson_plot.png", width = 9, height = 6)


# 04 Sensitive analysis ####
## 04.1 Sensitive analysis (index year = 2014) ####
data_analysis <- data_combined %>%
    filter(!is.na(shortage)) %>%
    mutate(
        location_id = factor(location_id),
        time_c = year - 2014,
        rate = (val / population) * 1000, # per 1000 births
        post = as.integer(year >= 2014),
        weights = population / sum(population)
    )
lmm_formula <- rate ~ time_c * shortage * post +
    hdi + sdi + gni + gbd_superregion +
    syphilis_prev_female + syphilis_prev_male + abor_misc_rate +
    (1 + time_c | location_id)
stepwise_lmm <- step(
    lmer(
        lmm_formula,
        data = data_analysis,
        weights = weights,
        REML = FALSE
    ),
    reduce.fixed = FALSE
)
stepwise_lmm

# get final model
stepwise_lmm <- get_model(stepwise_lmm)
stepwise_lmm
summary(stepwise_lmm)

# diagnostics
glmm_modelcheck <- performance::check_model(stepwise_lmm)
png(
    file = "Results/lmm_modelcheck_2014.png",
    width = 1200, height = 2400, res = 100
)
glmm_modelcheck
dev.off()

# summary
lmm_sum <- summary(stepwise_lmm)
lmm_sum
fixed_lmm <- as.data.frame(lmm_sum$coefficients)
fixed_lmm
fixed_lmm_results <- data.frame(
    Parameter = rownames(fixed_lmm),
    Estimate = fixed_lmm[, "Estimate"],
    SE = fixed_lmm[, "Std. Error"],
    t_value = fixed_lmm[, "t value"], # t-tests use Satterthwaite's method
    p_value = fixed_lmm[, "Pr(>|t|)"]
) %>%
    mutate(
        CI_lower = Estimate - 1.96 * SE,
        CI_upper = Estimate + 1.96 * SE
    )
fixed_lmm_results
write.xlsx(
    fixed_lmm_results,
    "Results/mixed_effects_linear_fixed_effects2014.xlsx"
)

jtools::plot_summs(stepwise_lmm,
    omit.coefs = c("(Intercept)", "hdi", "sdi", "gni"),
    model.names = "LMM (rate)",
    robust = TRUE
) +
    labs(
        title = "Linear Mixed Model: Estimates (95% CI)",
        x = "Estimate",
        y = "Parameter"
    )
ggsave("Results/mixed_lmm_plot2014.png", width = 9, height = 6)

## 04.2 Sensitive analysis (index year = 2015) ####
data_analysis <- data_combined %>%
    filter(!is.na(shortage)) %>%
    mutate(
        location_id = factor(location_id),
        time_c = year - 2015,
        rate = (val / population) * 1000, # per 1000 births
        post = as.integer(year >= 2015),
        weights = population / sum(population)
    )
lmm_formula <- rate ~ time_c * shortage * post +
    hdi + sdi + gni + gbd_superregion +
    syphilis_prev_female + syphilis_prev_male + abor_misc_rate +
    (1 + time_c | location_id)
stepwise_lmm <- step(
    lmer(
        lmm_formula,
        data = data_analysis,
        weights = weights,
        REML = FALSE
    ),
    reduce.fixed = FALSE
)
stepwise_lmm

# get final model
stepwise_lmm <- get_model(stepwise_lmm)
stepwise_lmm
summary(stepwise_lmm)

# diagnostics
glmm_modelcheck <- performance::check_model(stepwise_lmm)
glmm_modelcheck
png(
    file = "Results/lmm_modelcheck2015.png",
    width = 1200, height = 2400, res = 100
)
glmm_modelcheck
dev.off()

# summary
lmm_sum <- summary(stepwise_lmm)
lmm_sum
fixed_lmm <- as.data.frame(lmm_sum$coefficients)
fixed_lmm
fixed_lmm_results <- data.frame(
    Parameter = rownames(fixed_lmm),
    Estimate = fixed_lmm[, "Estimate"],
    SE = fixed_lmm[, "Std. Error"],
    t_value = fixed_lmm[, "t value"], # t-tests use Satterthwaite's method
    p_value = fixed_lmm[, "Pr(>|t|)"]
) %>%
    mutate(
        CI_lower = Estimate - 1.96 * SE,
        CI_upper = Estimate + 1.96 * SE
    )
fixed_lmm_results
write.xlsx(
    fixed_lmm_results,
    "Results/mixed_effects_linear_fixed_effects2015.xlsx"
)

jtools::plot_summs(stepwise_lmm,
    omit.coefs = c("(Intercept)", "hdi", "sdi", "gni"),
    model.names = "LMM (rate)",
    robust = TRUE
) +
    labs(
        title = "Linear Mixed Model: Estimates (95% CI)",
        x = "Estimate",
        y = "Parameter"
    )
ggsave("Results/mixed_lmm_plot2015.png", width = 9, height = 6)

## 04.3 Sensitive analysis (use the SstepcAIC to select model) ####
data_analysis <- data_combined %>%
    filter(!is.na(shortage)) %>%
    mutate(
        location_id = factor(location_id),
        time_c = year - 2016,
        rate = (val / population) * 1000, # per 1000 births
        post = as.integer(year >= 2016),
        weights = population / sum(population)
    )

lmm_formula <- rate ~ time_c * shortage * post +
    hdi + sdi + gni + gbd_superregion +
    syphilis_prev_female + syphilis_prev_male + abor_misc_rate +
    (1 + time_c | location_id)
stepwise_lmm <- cAIC4::stepcAIC(
    lmer(
        lmm_formula,
        data = data_analysis,
        weights = weights,
        REML = FALSE
    )
)

# get final model
stepwise_lmm <- stepwise_lmm$finalModel
stepwise_lmm
summary(stepwise_lmm)

# diagnostics
glmm_modelcheck <- performance::check_model(stepwise_lmm)
glmm_modelcheck
png(
    file = "Results/lmm_modelcheck_stepcAIC.png",
    width = 1200, height = 2400, res = 100
)
glmm_modelcheck
dev.off()

# summary
lmm_sum <- summary(stepwise_lmm)
lmm_sum
fixed_lmm <- as.data.frame(lmm_sum$coefficients)
fixed_lmm
fixed_lmm_results <- data.frame(
    Parameter = rownames(fixed_lmm),
    Estimate = fixed_lmm[, "Estimate"],
    SE = fixed_lmm[, "Std. Error"],
    t_value = fixed_lmm[, "t value"], # t-tests use Satterthwaite's method
    p_value = fixed_lmm[, "Pr(>|t|)"]
) %>%
    mutate(
        CI_lower = Estimate - 1.96 * SE,
        CI_upper = Estimate + 1.96 * SE
    )
fixed_lmm_results
write.xlsx(
    fixed_lmm_results,
    "Results/mixed_effects_linear_fixed_effects2016_stepAIC.xlsx"
)

jtools::plot_summs(stepwise_lmm,
    omit.coefs = c("(Intercept)", "hdi", "sdi", "gni"),
    model.names = "LMM (rate)",
    robust = TRUE
) +
    labs(
        title = "Linear Mixed Model: Estimates (95% CI)",
        x = "Estimate",
        y = "Parameter"
    )
ggsave("Results/mixed_lmm_plot2016_stepAIC.png", width = 9, height = 6)

## 04.5 Add covariates of prevalence of HIV ####
head(data_analysis)
summary(data_analysis$rate)

data_analysis <- data_combined %>%
    filter(!is.na(shortage)) %>%
    mutate(
        location_id = factor(location_id),
        time_c = year - 2016,
        rate = (val / population) * 1000, # per 1000 births
        post = as.integer(year >= 2016),
        weights = population / sum(population)
    )
lmm_formula <- rate ~ time_c * shortage * post +
    hdi + sdi + gni + gbd_superregion + hiv_prev_female + hiv_prev_male +
    syphilis_prev_female + syphilis_prev_male + abor_misc_rate +
    (1 + time_c | location_id)
stepwise_lmm <- step(
    lmer(
        lmm_formula,
        data = data_analysis,
        weights = weights,
        REML = FALSE
    ),
    reduce.fixed = FALSE
)
stepwise_lmm

# get final model
stepwise_lmm <- get_model(stepwise_lmm)
stepwise_lmm
summary(stepwise_lmm)

# diagnostics
glmm_modelcheck <- performance::check_model(stepwise_lmm)
png(
    file = "Results/AddHIVPrev_lmm_modelcheck.png",
    width = 1200, height = 2400, res = 100
)
glmm_modelcheck
dev.off()

# summary
lmm_sum <- summary(stepwise_lmm)
lmm_sum
fixed_lmm <- as.data.frame(lmm_sum$coefficients)
fixed_lmm
fixed_lmm_results <- data.frame(
    Parameter = rownames(fixed_lmm),
    Estimate = fixed_lmm[, "Estimate"],
    SE = fixed_lmm[, "Std. Error"],
    t_value = fixed_lmm[, "t value"], # t-tests use Satterthwaite's method
    p_value = fixed_lmm[, "Pr(>|t|)"]
) %>%
    mutate(
        CI_lower = Estimate - 1.96 * SE,
        CI_upper = Estimate + 1.96 * SE
    )
fixed_lmm_results
write.xlsx(
    fixed_lmm_results,
    "Results/AddHIVPrev_linear_mixed_effects_results.xlsx"
)

jtools::plot_summs(stepwise_lmm,
    omit.coefs = c("(Intercept)", "hdi", "sdi", "gni"),
    model.names = "LMM (rate)",
    robust = TRUE
) +
    labs(
        title = "Linear Mixed Model: Estimates (95% CI)",
        x = "Estimate",
        y = "Parameter"
    )
ggsave("Results/AddHIVPrev_lmm_coefplot.png",
    width = 9, height = 6
)

# 05 Subgroup analysis ####
## 05.1 GBD regions ####
gbd_superregion <- levels(data_combined$gbd_superregion)

data_combined <- data_combined %>%
    mutate(
        region = case_when(
            gbd_superregion == "Sub-Saharan Africa" ~ "Sub-Saharan Africa",
            gbd_superregion == "High-income" ~ "High-income",
            TRUE ~ "Other regions"
        ),
        region = factor(
            region,
            levels = c(
                "Sub-Saharan Africa",
                "High-income",                
                "Other regions"
            )
        )
    )
table(
    data_combined$region,
    data_combined$shortage
)
unique(data_combined$region)

results <- lapply(
    unique(data_combined$region),
    function(x) {
        message("Processing region: ", x)
        # data
        data_analysis <- data_combined %>%
            filter(!is.na(shortage)) %>%
            filter(region == x) %>%
            mutate(
                location_id = factor(location_id),
                time_c = year - 2016,
                rate = (val / population) * 1000, # per 1000 births
                post = as.integer(year >= 2016),
                weights = population / sum(population)
            )
        # step LMM
        lmm_formula <- rate ~ time_c * shortage * post +
            hdi + sdi + gni +
            syphilis_prev_female + syphilis_prev_male + abor_misc_rate +
            (1 + time_c | location_id)
        if (x == "Other regions") {
            lmm_formula <- rate ~ time_c * shortage * post +
                hdi + sdi + gni + gbd_superregion +
                syphilis_prev_female + syphilis_prev_male + abor_misc_rate +
                (1 + time_c | location_id)
        }
        stepwise_lmm <- step(
            lmer(
                lmm_formula,
                data = data_analysis,
                weights = weights,
                REML = FALSE
            ),
            reduce.fixed = FALSE
        )
        stepwise_lmm

        # get final model
        stepwise_lmm <- get_model(stepwise_lmm)
        stepwise_lmm
        summary(stepwise_lmm)

        # diagnostics
        glmm_modelcheck <- performance::check_model(stepwise_lmm)
        png(
            file = paste0("Results/lmm_modelcheck_", x, ".png"),
            width = 1200, height = 2400, res = 100
        )
        glmm_modelcheck
        dev.off()
        # results
        lmm_sum <- summary(stepwise_lmm)
        fixed_lmm <- as.data.frame(lmm_sum$coefficients)
        fixed_lmm_results <- data.frame(
            Subgroup = x,
            Parameter = rownames(fixed_lmm),
            Estimate = fixed_lmm[, "Estimate"],
            SE = fixed_lmm[, "Std. Error"],
            t_value = fixed_lmm[, "t value"], # t-tests use Satterthwaite's method
            p_value = fixed_lmm[, "Pr(>|t|)"]
        ) %>%
            mutate(
                CI_lower = Estimate - 1.96 * SE,
                CI_upper = Estimate + 1.96 * SE
            )
        return(fixed_lmm_results)
    }
)
gbd_results <- do.call(
    rbind,
    results
)
str(gbd_results)

write.xlsx(
    gbd_results,
    "Results/gbd_linear_mixed_effects_results.xlsx"
)

## 05.2 income level ####
data_income <- data_combined %>%
    filter(
        year >= 2010 & year <= 2021,
        !is.na(shortage)
    ) %>%
    left_join(
        data_combined %>%
            filter(year == 2021) %>%
            select(location_id, income_level) %>%
            rename(income_level2021 = income_level),
        by = c("location_id")
    )
income_level <- levels(data_income$income_level2021)
table(
    data_income$income_level2021,
    data_income$shortage
)

results <- lapply(
    income_level,
    function(x) {
        message("Processing income_level2021: ", x)
        # data
        data_analysis <- data_income %>%
            filter(!is.na(shortage)) %>%
            filter(income_level2021 == x) %>%
            mutate(
                location_id = factor(location_id),
                time_c = year - 2016,
                rate = (val / population) * 1000, # per 1000 births
                post = as.integer(year >= 2016),
                weights = population / sum(population)
            )
        # step LMM
        if (length(unique(data_analysis$gbd_superregion)) > 1) {
            gbd_term <- " + gbd_superregion"
        } else {
            gbd_term <- ""
        }

        lmm_formula <- as.formula(paste0(
            "rate ~ time_c * shortage * post +",
            " hdi + sdi + gni",
            gbd_term,
            " + syphilis_prev_female + syphilis_prev_male + abor_misc_rate +
            (1 + time_c | location_id)"
        ))
        # step LMM
        stepwise_lmm <- step(
            lmer(
                lmm_formula,
                data = data_analysis,
                weights = weights,
                REML = FALSE
            ),
            reduce.fixed = FALSE
        )
        stepwise_lmm <- get_model(stepwise_lmm)
        stepwise_lmm
        # diagnostics
        glmm_modelcheck <- performance::check_model(stepwise_lmm)
        png(
            file = paste0("Results/lmm_modelcheck_", x, ".png"),
            width = 1200, height = 2400, res = 100
        )
        glmm_modelcheck
        dev.off()
        # results
        lmm_sum <- summary(stepwise_lmm)
        fixed_lmm <- as.data.frame(lmm_sum$coefficients)
        fixed_lmm_results <- data.frame(
            Subgroup = x,
            Parameter = rownames(fixed_lmm),
            Estimate = fixed_lmm[, "Estimate"],
            SE = fixed_lmm[, "Std. Error"],
            t_value = fixed_lmm[, "t value"], # t-tests use Satterthwaite's method
            p_value = fixed_lmm[, "Pr(>|t|)"]
        ) %>%
            mutate(
                CI_lower = Estimate - 1.96 * SE,
                CI_upper = Estimate + 1.96 * SE
            )
        return(fixed_lmm_results)
    }
)
income_results <- do.call(
    rbind,
    results
)
str(income_results)

write.xlsx(
    income_results,
    "Results/income_linear_mixed_effects_results.xlsx"
)

## 05.5 CS Rate level ####
data_ratelevel <- data_combined %>%
    filter(
        year >= 2010 & year <= 2021,
        !is.na(shortage)
    ) %>%
    left_join(
        data_combined %>%
            filter(year == 2010) %>%
            mutate(
                rate = (val / population) * 1000, # per 1000 births
                ratelevel = case_when(
                    rate <= 0.5 ~ "Lower than 0.5",
                    rate > 0.5 ~ "Above 0.5"
                ),
                ratelevel = factor(
                    ratelevel,
                    levels = c("Lower than 0.5", "Above 0.5")
                )
            ) %>% 
            select(location_id, ratelevel),
        by = c("location_id")
    )
ratelevel <- levels(data_ratelevel$ratelevel)
table(
    data_ratelevel$ratelevel,
    data_ratelevel$shortage
)

results <- lapply(
    ratelevel,
    function(x) {
        message("Processing ratelevel: ", x)
        # data
        data_analysis <- data_ratelevel %>%
            filter(!is.na(shortage)) %>%
            filter(ratelevel == x) %>%
            mutate(
                location_id = factor(location_id),
                time_c = year - 2016,
                rate = (val / population) * 1000, # per 1000 births
                post = as.integer(year >= 2016),
                weights = population / sum(population)
            )
        # step LMM
        if (length(unique(data_analysis$gbd_superregion)) > 1) {
            gbd_term <- " + gbd_superregion"
        } else {
            gbd_term <- ""
        }
        lmm_formula <- as.formula(paste0(
            "rate ~ time_c * shortage * post +",
            " hdi + sdi + gni",
            gbd_term,
            " + syphilis_prev_female + syphilis_prev_male + abor_misc_rate +
            (1 + time_c | location_id)"
        ))
        stepwise_lmm <- step(
            lmer(
                lmm_formula,
                data = data_analysis,
                weights = weights,
                REML = FALSE
            ),
            reduce.fixed = FALSE
        )
        stepwise_lmm <- get_model(stepwise_lmm)
        stepwise_lmm
        # diagnostics
        glmm_modelcheck <- performance::check_model(stepwise_lmm)
        png(
            file = paste0("Results/lmm_modelcheck_", x, ".png"),
            width = 1200, height = 2400, res = 100
        )
        glmm_modelcheck
        dev.off()
        # results
        lmm_sum <- summary(stepwise_lmm)
        fixed_lmm <- as.data.frame(lmm_sum$coefficients)
        fixed_lmm_results <- data.frame(
            Subgroup = x,
            Parameter = rownames(fixed_lmm),
            Estimate = fixed_lmm[, "Estimate"],
            SE = fixed_lmm[, "Std. Error"],
            t_value = fixed_lmm[, "t value"], # t-tests use Satterthwaite's method
            p_value = fixed_lmm[, "Pr(>|t|)"]
        ) %>%
            mutate(
                CI_lower = Estimate - 1.96 * SE,
                CI_upper = Estimate + 1.96 * SE
            )
        return(fixed_lmm_results)
    }
)
ratelevel_results <- do.call(
    rbind,
    results
)
str(ratelevel_results)
write.xlsx(
    ratelevel_results,
    "Results/ratelevel_linear_mixed_effects_results.xlsx"
)

## 05.4 SDI qunintile ####
data_sdi <- data_combined %>%
    filter(
        year >= 2010 & year <= 2021,
        !is.na(shortage)
    ) %>%
    select(-sdi_level2021) %>%
    left_join(
        data_combined %>%
            filter(year == 2021) %>%
            select(location_id, sdi_level2021),
        by = c("location_id")
    )

sdi_level <- levels(data_sdi$sdi_level2021)
table(
    data_sdi$sdi_level2021,
    data_sdi$gbd_superregion
)

# results <- lapply(
#     sdi_level,
#     function(x) {
#         message("Processing sdi_level2021: ", x)
#         # data
#         data_analysis <- data_sdi %>%
#             filter(sdi_level2021 == x) %>%
#             mutate(
#                 location_id = factor(location_id),
#                 time_c = year - 2016,
#                 rate = (val / population) * 1000, # per 1000 births
#                 post = as.integer(year >= 2016),
#                 weights = population / sum(population)
#             )
#         # step LMM
#         if (length(unique(data_analysis$gbd_superregion)) > 1) {
#             gbd_term <- " + gbd_superregion"
#         } else {
#             gbd_term <- ""
#         }
#         lmm_formula <- as.formula(paste0(
#             "rate ~ time_c * shortage * post +",
#             " hdi + sdi + gni",
#             gbd_term,
#             " + syphilis_prev_female + syphilis_prev_male + abor_misc_rate +
#             (1 + time_c | location_id)"
#         ))
#         stepwise_lmm <- step(
#             lmer(
#                 lmm_formula,
#                 data = data_analysis,
#                 weights = weights,
#                 REML = FALSE
#             ),
#             reduce.fixed = FALSE
#         )
#         stepwise_lmm <- get_model(stepwise_lmm)
#         message("Stepwise LMM completed for sdi_level2021: ", x)

#         # diagnostics
#         glmm_modelcheck <- performance::check_model(stepwise_lmm)
#         png(
#             file = paste0("Results/lmm_modelcheck_", x, ".png"),
#             width = 1200, height = 2400, res = 100
#         )
#         glmm_modelcheck
#         dev.off()
#         # results
#         lmm_sum <- summary(stepwise_lmm)
#         fixed_lmm <- as.data.frame(lmm_sum$coefficients)
#         fixed_lmm_results <- data.frame(
#             Subgroup = x,
#             Parameter = rownames(fixed_lmm),
#             Estimate = fixed_lmm[, "Estimate"],
#             SE = fixed_lmm[, "Std. Error"],
#             t_value = fixed_lmm[, "t value"], # t-tests use Satterthwaite's method
#             p_value = fixed_lmm[, "Pr(>|t|)"]
#         ) %>%
#             mutate(
#                 CI_lower = Estimate - 1.96 * SE,
#                 CI_upper = Estimate + 1.96 * SE
#             )
#         return(fixed_lmm_results)
#     }
# )

# sdi_results <- do.call(
#     rbind,
#     results
# )
# str(sdi_results)

for (i in seq_along(sdi_level)) {
    x <- sdi_level[i]

    message("\n", paste0(rep("=", 50), collapse = ""))
    message("Processing sdi_level2021: ", x, " (", i, "/", length(sdi_level), ")")
    message(paste0(rep("=", 50), collapse = ""))

    # Data preparation
    tryCatch(
        {
            data_analysis1 <- data_sdi %>%
                filter(sdi_level2021 == x) %>%
                mutate(
                    location_id = factor(location_id),
                    time_c = year - 2016,
                    rate = (val / population) * 1000, # per 1000 births
                    post = as.integer(year >= 2016),
                    weights = population / sum(population, na.rm = TRUE)
                )

            # Check if we have enough data
            n_locations <- length(unique(data_analysis1$location_id))
            n_obs <- nrow(data_analysis1)
            message("Data dimensions: ", n_obs, " observations, ", n_locations, " locations")

            # Skip if insufficient data
            if (n_obs < 10 | n_locations < 3) {
                message("Insufficient data, skipping sdi_level: ", x)
                next
            }

            # Determine if gbd_superregion should be included
            gbd_term <- ifelse(length(unique(data_analysis1$gbd_superregion)) > 1,
                " + gbd_superregion", ""
            )

            # Construct formula
            lmm_formula <- as.formula(paste0(
                "rate ~ time_c * shortage * post +",
                " hdi + sdi + gni",
                gbd_term,
                " + syphilis_prev_female + syphilis_prev_male + abor_misc_rate +",
                " (1 + time_c | location_id)"
            ))

            message("Model formula: ", format(lmm_formula))

            # Fit initial model
            message("Fitting initial model...")
            initial_fit <- lmer(
                lmm_formula,
                data = data_analysis1,
                weights = weights,
                REML = FALSE,
                control = lmerControl(
                    optimizer = "bobyqa",
                    optCtrl = list(maxfun = 2e5),
                    check.nobs.vs.nRE = "ignore"
                )
            )

            # Stepwise selection (only for random effects if reduce.fixed = FALSE)
            message("Performing stepwise selection...")
            tryCatch(
                {
                    stepwise_result <- step(
                        initial_fit,
                        reduce.fixed = FALSE
                    )

                    stepwise_lmm <- get_model(stepwise_result)
                    message("Stepwise selection completed")
                },
                error = function(e) {
                    message("Stepwise selection failed: ", e$message)
                    message("Using initial model instead")
                    stepwise_lmm <- initial_fit
                }
            )

            # Diagnostics
            message("Generating diagnostics...")
            tryCatch(
                {
                    glmm_modelcheck <- performance::check_model(stepwise_lmm)
                    png(
                        file = paste0("Results/lmm_modelcheck_", gsub("[^[:alnum:]]", "_", x), ".png"),
                        width = 1200, height = 2400, res = 100
                    )
                    print(glmm_modelcheck)
                    dev.off()
                    message("Diagnostic plot saved")
                },
                error = function(e) {
                    message("Could not create diagnostic plot: ", e$message)
                }
            )

            # Extract results
            message("Extracting results...")
            tryCatch(
                {
                    lmm_sum <- summary(stepwise_lmm)
                    fixed_lmm <- as.data.frame(lmm_sum$coefficients)

                    fixed_lmm_results <- data.frame(
                        Subgroup = x,
                        Parameter = rownames(fixed_lmm),
                        Estimate = fixed_lmm[, "Estimate"],
                        SE = fixed_lmm[, "Std. Error"],
                        t_value = fixed_lmm[, "t value"],
                        p_value = fixed_lmm[, "Pr(>|t|)"]
                    ) %>%
                        mutate(
                            CI_lower = Estimate - 1.96 * SE,
                            CI_upper = Estimate + 1.96 * SE
                        )

                    # Store results
                    results[[as.character(x)]] <- list(
                        coefficients = fixed_lmm_results,
                        model_summary = lmm_sum,
                        model_object = stepwise_lmm,
                        data_info = c(N_obs = n_obs, N_locations = n_locations)
                    )

                    message("Successfully processed sdi_level: ", x)
                },
                error = function(e) {
                    message("Error extracting results: ", e$message)
                }
            )

            # Clean up to free memory
            rm(data_analysis1, initial_fit, stepwise_lmm)
            gc()
        },
        error = function(e) {
            message("Error processing sdi_level ", x, ": ", e$message)
        }
    )

    # Progress update
    progress <- round(i / length(sdi_level) * 100, 1)
    message("\nProgress: ", progress, "% complete")
    message(paste0(rep("-", 50), collapse = ""))
}

sdi_results <- do.call(rbind, lapply(results, function(x) x$coefficients))
rownames(sdi_results) <- NULL

write.xlsx(
    sdi_results,
    "Results/sdi_linear_mixed_effects_results.xlsx"
)

## 05.3 HDI level ####
data_hdi <- data_combined %>%
    filter(
        year >= 2010 & year <= 2021,
        !is.na(shortage)
    ) %>%
    left_join(
        data_combined %>%
            filter(year == 2021) %>%
            select(location_id, hdi_level) %>%
            rename(hdi_level2021 = hdi_level),
        by = c("location_id")
    )
hdi_level <- levels(data_hdi$hdi_level2021)
table(
    data_hdi$hdi_level2021,
    data_hdi$shortage
)

results <- lapply(
    hdi_level,
    function(x) {
        message("Processing hdi_level2021: ", x)
        # data
        data_analysis <- data_hdi %>%
            filter(!is.na(shortage)) %>%
            filter(hdi_level2021 == x) %>%
            mutate(
                location_id = factor(location_id),
                time_c = year - 2016,
                rate = (val / population) * 1000, # per 1000 births
                post = as.integer(year >= 2016),
                weights = population / sum(population)
            )
        # step LMM
        if (length(unique(data_analysis$gbd_superregion)) > 1) {
            gbd_term <- " + gbd_superregion"
        } else {
            gbd_term <- ""
        }
        lmm_formula <- as.formula(paste0(
            "rate ~ time_c * shortage * post +",
            " hdi + sdi + gni",
            gbd_term,
            " + syphilis_prev_female + syphilis_prev_male + abor_misc_rate +
            (1 + time_c | location_id)"
        ))
        stepwise_lmm <- step(
            lmer(
                lmm_formula,
                data = data_analysis,
                weights = weights,
                REML = FALSE
            ),
            reduce.fixed = FALSE
        )
        stepwise_lmm <- get_model(stepwise_lmm)
        stepwise_lmm
        # diagnostics
        glmm_modelcheck <- performance::check_model(stepwise_lmm)
        png(
            file = paste0("Results/lmm_modelcheck_", x, ".png"),
            width = 1200, height = 2400, res = 100
        )
        glmm_modelcheck
        dev.off()
        # results
        lmm_sum <- summary(stepwise_lmm)
        fixed_lmm <- as.data.frame(lmm_sum$coefficients)
        fixed_lmm_results <- data.frame(
            Subgroup = x,
            Parameter = rownames(fixed_lmm),
            Estimate = fixed_lmm[, "Estimate"],
            SE = fixed_lmm[, "Std. Error"],
            t_value = fixed_lmm[, "t value"], # t-tests use Satterthwaite's method
            p_value = fixed_lmm[, "Pr(>|t|)"]
        ) %>%
            mutate(
                CI_lower = Estimate - 1.96 * SE,
                CI_upper = Estimate + 1.96 * SE
            )
        return(fixed_lmm_results)
    }
)
hdi_results <- do.call(
    rbind,
    results
)
str(hdi_results)

write.xlsx(
    hdi_results,
    "Results/hdi_linear_mixed_effects_results.xlsx"
)

## 05.5 Forestplot (forestplot package) for subgroup effects -------
# P_for_interaction
p_for_interaction_fixed <- function(beta, se) {
    # beta: vector of subgroup effects (length K)
    # se:   vector of standard errors (same length)

    if (length(beta) != length(se)) {
        stop("beta and se must have the same length.")
    }

    K <- length(beta)
    w <- 1 / se^2 # inverse-variance weights
    beta_fixed <- sum(w * beta) / sum(w) # pooled fixed-effect estimate

    # Cochran's Q for subgroup differences
    Q <- sum(w * (beta - beta_fixed)^2)

    df <- K - 1
    p_value <- pchisq(Q, df = df, lower.tail = FALSE)

    return(list(
        Q = Q,
        df = df,
        p_interaction = p_value,
        beta_pooled = beta_fixed
    ))
}

p_for_interaction_random <- function(beta, se) {
    # beta: subgroup pooled effects from random-effects models
    # se:   subgroup standard errors (from each random-effects model)

    if (length(beta) != length(se)) {
        stop("beta and se must have the same length.")
    }

    K <- length(beta)
    var_k <- se^2
    w_k <- 1 / var_k

    # pooled effect across subgroups (random-effects subgroup test uses IV pooled effect)
    beta_pooled <- sum(w_k * beta) / sum(w_k)

    # Q statistic for subgroup differences
    Q_subgroup <- sum(w_k * (beta - beta_pooled)^2)

    df <- K - 1
    p_value <- pchisq(Q_subgroup, df = df, lower.tail = FALSE)

    return(list(
        Q = Q_subgroup,
        df = df,
        p_interaction = p_value,
        beta_pooled = beta_pooled
    ))
}

# combinde all subgroup effect results
subgroup_names <- c("Region", "Income Level", "HDI Level", "SDI Level", "Rate level")
subgroup_effects <- list(
    Region = gbd_results,
    Income_Level = income_results,
    HDI_Level = hdi_results,
    SDI_Level = sdi_results,
    Rate_level = ratelevel_results
)

sdi_results

# prepare level_difference
level_difference <- lapply(names(subgroup_effects), function(x) {
    df <- cbind(
        data.frame(Group = c(subgroup_names[which(names(subgroup_effects) == x)])),
        subgroup_effects[[x]] %>%
            filter(Parameter == "shortage:post") %>%
            select(Subgroup, Estimate, SE, CI_lower, CI_upper)
    )
    p_results <- p_for_interaction_random(beta = df$Estimate, se = df$SE)
    rf <- cbind(
        df,
        data.frame(
            Q = c(p_results$Q, rep(NA, nrow(df) - 1)),
            p_value = c(p_results$p_interaction, rep(NA, nrow(df) - 1)),
            Estimated_pooled = c(p_results$beta_pooled, rep(NA, nrow(df) - 1))
        )
    )
    rownames(rf) <- NULL
    rf
}) %>%
    do.call(rbind, .)
level_difference

# prepare slope_difference analogously
slope_difference <- lapply(names(subgroup_effects), function(x) {
    df <- cbind(
        data.frame(Group = c(subgroup_names[which(names(subgroup_effects) == x)])),
        subgroup_effects[[x]] %>%
            filter(grepl("time_c:shortage:post", Parameter)) %>%
            select(Subgroup, Estimate, SE, CI_lower, CI_upper)
    )
    p_results <- p_for_interaction_fixed(beta = df$Estimate, se = df$SE)
    rf <- cbind(
        df,
        data.frame(
            Q = c(p_results$Q, rep(NA, nrow(df) - 1)),
            p_value = c(p_results$p_interaction, rep(NA, nrow(df) - 1)),
            Estimated_pooled = c(p_results$beta_pooled, rep(NA, nrow(df) - 1))
        )
    )
    rownames(rf) <- NULL
    rf
}) %>%
    do.call(rbind, .)
slope_difference

library(forestploter)
library(grid)
combined <- inner_join(level_difference,
    slope_difference,
    by = c("Group", "Subgroup"),
    suffix = c("_level", "_slope")
)
View(combined)

# formatted CI text columns used as label columns by forestploter
combined$Group <- c(
    "Region", rep("", 2), "Income level", rep("", 3),
    "HDI level", rep("", 3), "SDI qunintile", rep("", 4), "Rate level", rep("", 1)
)

str(combined)
combined$`Level change (95%CI) ` <- sprintf(
    "%.3f (%.3f, %.3f)",
    combined$Estimate_level,
    combined$CI_lower_level,
    combined$CI_upper_level
)
combined$`Level change (95%CI)` <- paste(rep(" ", 40), collapse = "")
combined$`P for interaction\n(Level change)` <- ifelse(
    is.na(combined$p_value_level), "",
    ifelse(combined$p_value_level < 0.001, "<0.001",
        sprintf("%.3f", combined$p_value_level)
    )
)

combined$`Slope change (95%CI) ` <- sprintf(
    "%.3f (%.3f, %.3f)",
    combined$Estimate_slope,
    combined$CI_lower_slope,
    combined$CI_upper_slope
)
combined$`Slope change (95%CI)` <- paste(rep(" ", 40), collapse = "")
combined$`P for interaction\n(Slope change)` <- ifelse(
    is.na(combined$p_value_slope), "",
    ifelse(combined$p_value_slope < 0.001, "<0.001",
        sprintf("%.3f", combined$p_value_slope)
    )
)

pal <- ggsci::pal_lancet(palette = c("lanonc"), alpha = 1)

# theme for forest
tm <- forestploter::forest_theme(
    base_size = 10,
    refline_gp = grid::gpar(lty = "dashed", col = "#666666"),
    ci_pch = 15,
    ci_col = pal(2)[1],
    ci_fill = pal(2)[1],
    arrow_type = "closed",
    ci_Theight = 0.25
)

str(combined)
fp_grob <- forestploter::forest(
    combined[, c(1, 2, 17, 18, 19, 20, 21, 22)],
    labeltext = combined[, c("Group", "Subgroup", "Level change (95%CI) ", "Slope change (95%CI) ")],
    est = list(combined$Estimate_level, combined$Estimate_slope),
    lower = list(combined$CI_lower_level, combined$CI_lower_slope),
    upper = list(combined$CI_upper_level, combined$CI_upper_slope),
    ci_colum = c(4, 7),
    arrow_lab = c("Decreased rate", "Increased rate"),
    # ticks_at = list(
    #     seq(-1, 1, by = 0.1),
    #     seq(-0.1, 0.1, by = 0.01)
    # ),
    xlab = "",
    theme = tm
)
fp_grob <- edit_plot(fp_grob,
    row = seq_len(nrow(combined)),
    col = 7,
    which = "ci",
    gp = gpar(pch = 21, fill = pal(2)[2], col = pal(2)[2])
)
fp_grob <- edit_plot(fp_grob,
    row = which(combined$Group == ""),
    which = "background",
    gp = gpar(fill = "white")
)
fp_grob <- edit_plot(fp_grob,
    row = which(combined$Group != ""),
    which = "background",
    gp = gpar(fill = "gray80")
)
fp_grob <- edit_plot(fp_grob,
    row = which(combined$Group != ""),
    which = "background",
    gp = gpar(fill = "gray80")
)
fp_grob

png(paste0("Results/lmm_combined_forest_forestploter.png"),
    width = 14, height = 9, units = "in", res = 300
)
fp_grob
dev.off()
