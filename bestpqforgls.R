#' Select optimal ARMA orders for GLS model with DW test
#' @param data Data frame containing the time series data
#' @param max_p Maximum AR order to test
#' @param max_q Maximum MA order to test
#' @param lb_lag Lag for Ljung-Box test
#' @return Data frame with model selection results
select_optimal_arma <- function(
    data,
    formula = as.formula(paste0(
        "proportion ~ time + intervention + postslope + shortage +",
        "time:shortage + intervention:shortage + postslope:shortage"
    )),
    max_p = 3,
    max_q = 3,
    lb_lag = 4) {
    # Create results dataframe
    results <- expand.grid(p = 0:max_p, q = 0:max_q)
    results$AIC <- NA
    results$BIC <- NA
    results$LB_pvalue <- NA
    results$DW_stat <- NA
    results$DW_pvalue <- NA
    results$converged <- FALSE
    results$valid <- FALSE

    # Try different ARMA orders
    for (i in seq_len(nrow(results))) {
        p <- results$p[i]
        q <- results$q[i]

        # Fit model
        model <- try(nlme::gls(
            formula,
            data = data,
            correlation = if (p > 0 || q > 0) {
                nlme::corARMA(form = ~ time | shortage, p = p, q = q)
            } else {
                NULL
            },
            method = "ML"
        ), silent = TRUE)

        # Check model and diagnostics
        if (!inherits(model, "try-error")) {
            residuals <- residuals(model, type = "normalized")
            lb_test <- Box.test(residuals, lag = lb_lag, type = "Ljung-Box")
            dw_test <- lmtest::dwtest(residuals ~ 1)

            results$AIC[i] <- AIC(model)
            results$BIC[i] <- BIC(model)
            results$LB_pvalue[i] <- lb_test$p.value
            results$DW_stat[i] <- dw_test$statistic
            results$DW_pvalue[i] <- dw_test$p.value
            results$dw_delta[i] <- abs(dw_test$statistic - 2)
            results$converged[i] <- TRUE
            results$valid[i] <- lb_test$p.value > 0.05 & dw_test$p.value > 0.05
        }
    }
    return(results)
}

#' Find best ARMA model that passes all diagnostics
#' @param results Data frame from select_optimal_arma
#' @return List containing best model specifications
get_best_model <- function(results) {
    # Check if results is valid
    if (missing(results) || is.null(results) || nrow(results) == 0) {
        return(list(found_valid = FALSE))
    }

    # Ensure dw_delta column exists
    if (!"dw_delta" %in% names(results) ||
        !"converged" %in% names(results) ||
        !"valid" %in% names(results)) {
        results$dw_delta <- NA
    }

    valid_models <- results %>%
        filter(converged & valid) %>%
        arrange(dw_delta, AIC)

    if (nrow(valid_models) > 0) {
        best <- valid_models[1, ]
        return(list(
            best_p = best$p,
            best_q = best$q,
            AIC = best$AIC,
            BIC = best$BIC,
            LB_pvalue = best$LB_pvalue,
            DW_stat = best$DW_stat,
            DW_pvalue = best$DW_pvalue,
            found_valid = TRUE
        ))
    } else {
        return(list(
            found_valid = FALSE
        ))
    }
}


# ' Perform enhanced diagnostics for GLS model with ARMA correlation
# ' @param data Data frame containing the time series data
# ' @param subgroup Variable name for subgroup analysis
# ' @param value Value of the subgroup variable to filter data
# ' @param filepath File path to save diagnostic plots
perform_enhanced_diagnostics <- function(data,
                                         subgroup = NULL,
                                         value = NULL,
                                         filepath = "Results_1121/Diagnostics/",
                                         formula = as.formula(paste0(
                                             "proportion ~ time + intervention + postslope + shortage +",
                                             "time:shortage + intervention:shortage + postslope:shortage"
                                         )),
                                         xvar = "year",
                                         max_p = 3,
                                         max_q = 3,
                                         lb_lag = 4) {
    # Filter data for specific subgroup
    if (is.null(subgroup) & is.null(value)) {
        value <- "overall"
        subgroup_data <- data
        filepath <- paste0(filepath, "overall_diagnostics.pdf")
    } else if (!subgroup %in% names(data)) {
        stop("Subgroup variable not found in data.")
    } else {
        data$subgroup <- data[[subgroup]]
        subgroup_data <- data %>% filter(subgroup == value)
        filepath <- paste0(filepath, "subgroup_", value, "_diagnostics.pdf")
    }

    # Test different ARMA orders
    arma_results <- select_optimal_arma(
        data = subgroup_data,
        max_p = max_p,
        max_q = max_p,
        lb_lag = lb_lag
    )

    # Get best model specifications
    best_model <- get_best_model(arma_results)
    if (best_model$found_valid) {
        if ((best_model$best_p > 0 | best_model$best_q > 0)) {
            # Fit final model with optimal orders
            final_model <- gls(
                proportion ~ time + intervention + postslope + shortage +
                    time:shortage + intervention:shortage + postslope:shortage,
                data = subgroup_data,
                correlation = corARMA(
                    form = ~ time | shortage,
                    p = best_model$best_p,
                    q = best_model$best_q
                ),
                method = "ML"
            )

            # Print results
            cat("Optimal ARMA orders:\n")
            cat("p =", best_model$best_p, "\n")
            cat("q =", best_model$best_q, "\n")
            cat("AIC =", best_model$AIC, "\n")
            cat("Ljung-Box p-value =", best_model$LB_pvalue, "\n")
        } else if ((best_model$best_p == 0 & best_model$best_q == 0)) {
            # Fit final model with optimal orders
            final_model <- gls(
                proportion ~ time + intervention + postslope + shortage +
                    time:shortage + intervention:shortage + postslope:shortage,
                data = subgroup_data,
                method = "ML"
            )

            # Print results
            cat("Optimal ARMA orders:\n")
            cat("p =", 0, "\n")
            cat("q =", 0, "\n")
            cat("AIC =", best_model$AIC, "\n")
            cat("Ljung-Box p-value =", best_model$LB_pvalue, "\n")
        }
    } else {
        cat("No valid model found. Consider increasing max_p and max_q.\n")
    }

    # Extract residuals from best model
    residuals <- residuals(final_model, type = "normalized")

    # Perform diagnostic tests
    dw_test <- dwtest(residuals ~ 1)
    lb_tests <- lapply(c(1, 4, 8), function(lag) {
        Box.test(residuals, lag = lag, type = "Ljung-Box")
    })

    # Create diagnostic plots
    pdf(filepath, width = 8, height = 6)
    par(mfrow = c(2, 2))
    # Residuals vs Time
    x <- subgroup_data[[xvar]]
    plot(x, residuals,
        type = "o",
        main = paste("Residuals vs Time", value),
        xlab = xvar,
        ylab = "Normalized Residuals"
    )
    abline(h = 0, col = "red", lty = 2)
    # ACF
    acf(residuals, main = "ACF of Residuals")
    # PACF
    pacf(residuals, main = "PACF of Residuals")
    # QQ Plot
    qqnorm(residuals, main = "Normal Q-Q Plot")
    qqline(residuals)
    dev.off()

    # Return results
    list(
        Subgroup = value,
        AR = best_model$best_p,
        MR = best_model$best_q,
        DW_Stat = dw_test$statistic,
        DW_PValue = dw_test$p.value,
        LB_Results = data.frame(
            Lag = c(1, 4, 8),
            Statistic = sapply(lb_tests, function(x) x$statistic),
            PValue = sapply(lb_tests, function(x) x$p.value)
        ),
        AIC = AIC(final_model),
        BIC = BIC(final_model)
    )
}


#' Perform ARMA ITS analysis for specified subgroup
#' @param data Data frame containing the time series data
#' @param subgroup Variable name for subgroup analysis
#' @param value Value of the subgroup variable to filter data
#' @param filepath File path to save diagnostic plots
#' @export


armaits_analysis <- function(data,
                             subgroup = NULL,
                             value = NULL,
                             filepath = "Results_1121/",
                             formula = as.formula(paste0(
                                 "proportion ~ time + intervention + postslope + shortage +",
                                 "time:shortage + intervention:shortage + postslope:shortage"
                             )),
                             xvar = "year",
                             max_p = 3,
                             max_q = 3,
                             lb_lag = 4) {
    # Filter data for specific region
    if (is.null(subgroup) & is.null(value)) {
        value <- "overall"
        subgroup_data <- data
        filepath <- paste0(filepath, "GLS_overall_diagnostics.pdf")
    } else if (!subgroup %in% names(data)) {
        stop("Subgroup variable not found in data.")
    } else {
        data$subgroup <- data[[subgroup]]
        subgroup_data <- data %>% filter(subgroup == value)
        filepath <- paste0(filepath, "GLS_", value, "_diagnostics.pdf")
    }

    # Test different ARMA orders
    arma_results <- select_optimal_arma(
        data = subgroup_data,
        max_p = max_p,
        max_q = max_q,
        lb_lag = lb_lag
    )
    # Get best model specifications
    best_model <- get_best_model(arma_results)
    if (best_model$found_valid) {
        if ((best_model$best_p > 0 | best_model$best_q > 0)) {
            # Fit final model with optimal orders
            final_model <- gls(
                formula,
                data = subgroup_data,
                correlation = corARMA(
                    form = ~ time | shortage,
                    p = best_model$best_p,
                    q = best_model$best_q
                ),
                method = "ML"
            )

            # Print results
            cat("Optimal ARMA orders:\n")
            cat("p =", best_model$best_p, "\n")
            cat("q =", best_model$best_q, "\n")
            cat("AIC =", best_model$AIC, "\n")
            cat("Ljung-Box p-value =", best_model$LB_pvalue, "\n")
        } else if ((best_model$best_p == 0 & best_model$best_q == 0)) {
            # Fit final model with optimal orders
            final_model <- gls(
                proportion ~ time + intervention + postslope + shortage +
                    time:shortage + intervention:shortage + postslope:shortage,
                data = subgroup_data,
                method = "ML"
            )

            # Print results
            cat("Optimal ARMA orders:\n")
            cat("p =", 0, "\n")
            cat("q =", 0, "\n")
            cat("AIC =", best_model$AIC, "\n")
            cat("Ljung-Box p-value =", best_model$LB_pvalue, "\n")
        }
    } else {
        cat("No valid model found. Consider increasing max_p and max_q.\n")
    }

    # Extract model results
    summary_table <- summary(final_model)$tTable

    # Create results dataframe
    results <- data.frame(
        Subgroup = value,
        Parameter = rownames(summary_table),
        Estimate = summary_table[, "Value"],
        SE = summary_table[, "Std.Error"],
        t_value = summary_table[, "t-value"],
        p_value = summary_table[, "p-value"]
    ) %>%
        mutate(
            CI_Lower = Estimate - 1.96 * SE,
            CI_Upper = Estimate + 1.96 * SE
        )

    return(results)
}
