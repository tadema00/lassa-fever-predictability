# ============================================================================
# 04_spi_sarimax_zinb_pooled_wildboot.R
# Combines:04a_spi_sarimax_zinb_models.R        - SARIMAX grid, ZINB robustness
#   04b_pooled_nb_cluster_interaction.R  - pooled NB with state FE & interaction
#   04c_wild_cluster_bootstrap.R         - wild bootstrap for interaction term
# ============================================================================

source("000_setup.R")

library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(readr)
library(forecast)
library(lmtest)
library(pscl)
library(sandwich)   # for vcovCL

# ----------------------------------------------------------------------------
# 1. LOAD DATA
# ----------------------------------------------------------------------------

prepped <- readRDS(file.path(OUT_DIR, "prepared_data.rds"))
weekly  <- prepped$weekly

lag_results <- readRDS(file.path(OUT_DIR, "lag_tables_sensitivity.rds"))
lag_table   <- lag_results$interp$lag_table

SPI_SCALE <- 4
clim_col  <- paste0("SPI_", SPI_SCALE)

# ----------------------------------------------------------------------------
# 2. SARIMAX GRID DEFINITION
# ----------------------------------------------------------------------------
ORDER_GRID <- expand.grid(
  p = 0:1, d = 0:1, q = 0:1,
  P = 0:1, D = 0:1, Q = 0:1
)
cat("SARIMAX candidate models:", nrow(ORDER_GRID), "\n")

# ----------------------------------------------------------------------------
# 3. ROBUST SARIMAX FITTING FUNCTION
# ----------------------------------------------------------------------------
fit_best_sarimax <- function(y, xreg = NULL, period = 52, model_name = "SARIMAX") {
  best <- list(aicc = Inf, aic = Inf, bic = Inf, fit = NULL,
               order = NULL, seasonal = NULL, model_name = model_name)
  successful_models <- list()
  
  for (i in seq_len(nrow(ORDER_GRID))) {
    o <- ORDER_GRID[i, ]
    fit <- tryCatch(
      forecast::Arima(y,
                      order = c(o$p, o$d, o$q),
                      seasonal = list(order = c(o$P, o$D, o$Q), period = period),
                      xreg = xreg,
                      include.drift = FALSE,
                      method = "ML"),
      error = function(e) NULL,
      warning = function(w) invokeRestart("muffleWarning")
    )
    if (is.null(fit)) next
    
    aic_val <- tryCatch(AIC(fit), error = function(e) Inf)
    aicc_val <- tryCatch(fit$aicc, error = function(e) Inf)
    bic_val <- tryCatch(BIC(fit), error = function(e) Inf)
    
    if (!is.finite(aic_val) || !is.finite(aicc_val) || !is.finite(bic_val)) next
    
    successful_models[[length(successful_models) + 1]] <- tibble(
      model = paste0("(", o$p, ",", o$d, ",", o$q, ")x(", o$P, ",", o$D, ",", o$Q, ")[", period, "]"),
      p = o$p, d = o$d, q = o$q, P = o$P, D = o$D, Q = o$Q,
      AIC = aic_val, AICc = aicc_val, BIC = bic_val
    )
    
    if (aicc_val < best$aicc) {
      best <- list(aicc = aicc_val, aic = aic_val, bic = bic_val,
                   fit = fit,
                   order = c(o$p, o$d, o$q),
                   seasonal = c(o$P, o$D, o$Q),
                   model_name = model_name)
    }
  }
  
  best$candidate_table <- if (length(successful_models) > 0) bind_rows(successful_models) %>% arrange(AICc) else tibble()
  return(best)
}

# ============================================================================
# 4. FIT ALL STATES (SARIMAX + ZINB) 
# All model selection and coefficient testing now use TRAINING DATA (2018-2023)
# ============================================================================

results_list <- list()
candidate_results <- list()

# Define the start of the validation period to ensure consistency
VALIDATION_START_YEAR <- 2024

for (st in FOCAL_STATES) {
  cat("\n============================================\nProcessing state:", st, "\n============================================\n")
  
  # 1. DATA SPLITTING: Create training-only dataset to prevent leakage
  df_state_full <- weekly %>% filter(State == st) %>% arrange(week_idx)
  df_state_train <- df_state_full %>% filter(Year < VALIDATION_START_YEAR)
  
  # 2. RETRIEVE LAG: Use the lag calculated from training data (from 03a)
  cons_lag <- lag_table %>% filter(State == st) %>% pull(Consensus_lag_wk)
  if (length(cons_lag) != 1) stop(sprintf("Expected one consensus lag for %s, got %d.", st, length(cons_lag)))
  
  # 3. PREPARE TRAINING VECTORS
  x_lagged_train <- dplyr::lag(df_state_train[[clim_col]], cons_lag)
  valid_train    <- is.finite(x_lagged_train) & is.finite(df_state_train$Cases)
  
  y_train_ts <- ts(df_state_train$Cases[valid_train], frequency = 52)
  x_train_mat <- matrix(x_lagged_train[valid_train], ncol = 1, dimnames = list(NULL, "SPI_lag"))
  
  cat("  Consensus lag (Training-based):", cons_lag, "weeks\n")
  cat("  Training observations:", length(y_train_ts), "\n")
  
  # 4. SARIMAX MODEL SELECTION (AICc) ON TRAINING DATA ONLY
  # This finds the best p,d,q structure based only on historical data.
  informed <- fit_best_sarimax(y = y_train_ts, xreg = x_train_mat, period = 52, 
                               model_name = "Climate-informed SARIMAX")
  
  naive    <- fit_best_sarimax(y = y_train_ts, xreg = NULL, period = 52, 
                               model_name = "Climate-naive SARIMA")
  
  # 5. EXTRACT TRAINING COEFFICIENTS
  spi_beta <- NA_real_; spi_se <- NA_real_; spi_p <- NA_real_
  if (!is.null(informed$fit)) {
    coef_tab <- informed$fit$coef
    if ("SPI_lag" %in% names(coef_tab)) {
      spi_beta <- unname(coef_tab["SPI_lag"])
      se_tab <- tryCatch(sqrt(diag(informed$fit$var.coef)), error = function(e) NULL)
      if (!is.null(se_tab) && "SPI_lag" %in% names(se_tab)) {
        spi_se <- unname(se_tab["SPI_lag"])
        spi_p <- 2 * pnorm(-abs(spi_beta / spi_se))
      }
    }
  }
  
  sarimax_order_str <- if (!is.null(informed$fit)) paste0("(", paste(informed$order, collapse = ","), ")x(", paste(informed$seasonal, collapse = ","), ")52") else NA_character_
  naive_order_str   <- if (!is.null(naive$fit))    paste0("(", paste(naive$order, collapse = ","), ")x(", paste(naive$seasonal, collapse = ","), ")52") else NA_character_
  
  delta_aic  <- if (is.finite(informed$aic) && is.finite(naive$aic))  naive$aic - informed$aic  else NA_real_
  delta_aicc <- if (is.finite(informed$aicc) && is.finite(naive$aicc)) naive$aicc - informed$aicc else NA_real_
  delta_bic  <- if (is.finite(informed$bic) && is.finite(naive$bic))   naive$bic - informed$bic  else NA_real_
  
  # 6. RESIDUAL DIAGNOSTICS (ON TRAINING FIT)
  residual_mean <- residual_sd <- lb_12 <- lb_24 <- NA_real_
  if (!is.null(informed$fit)) {
    resids <- residuals(informed$fit)
    resids <- resids[is.finite(resids)]
    if (length(resids) > 0) {
      residual_mean <- mean(resids); residual_sd <- sd(resids)
      n_params <- length(informed$fit$coef)
      # Check if the chosen model adequately handles training autocorrelation
      lb_12 <- tryCatch(Box.test(resids, lag = 12, type = "Ljung-Box", fitdf = n_params)$p.value, error = function(e) NA_real_)
      lb_24 <- tryCatch(Box.test(resids, lag = 24, type = "Ljung-Box", fitdf = n_params)$p.value, error = function(e) NA_real_)
    }
  }
  
  # 7. ZINB ROBUSTNESS (ON TRAINING DATA ONLY)
  zinb_df <- tibble(cases = round(as.numeric(y_train_ts)), spi = as.numeric(x_train_mat))
  zinb_fit <- tryCatch(pscl::zeroinfl(cases ~ spi | 1, data = zinb_df, dist = "negbin"), error = function(e) NULL)
  zinb_naive_fit <- tryCatch(pscl::zeroinfl(cases ~ 1 | 1, data = zinb_df, dist = "negbin"), error = function(e) NULL)
  
  zinb_lr <- zinb_p <- zinb_coef <- zinb_coef_p <- NA_real_
  if (!is.null(zinb_fit) && !is.null(zinb_naive_fit)) {
    lr <- tryCatch(lmtest::lrtest(zinb_naive_fit, zinb_fit), error = function(e) NULL)
    if (!is.null(lr)) { zinb_lr <- lr$Chisq[2]; zinb_p <- lr$`Pr(>Chisq)`[2] }
    coef_tab <- tryCatch(summary(zinb_fit)$coefficients$count, error = function(e) NULL)
    if (!is.null(coef_tab) && "spi" %in% rownames(coef_tab)) {
      zinb_coef <- coef_tab["spi", "Estimate"]
      zinb_coef_p <- coef_tab["spi", "Pr(>|z|)"]
    }
  }
  
  # 8. CONSOLIDATE TRAINING-BASED RESULTS
  results_list[[st]] <- tibble(
    State = st,
    Consensus_lag_wk = cons_lag,
    N_train = length(y_train_ts),
    SPI_beta_SARIMAX_train = round(spi_beta, 4),
    SPI_SE_train = round(spi_se, 4),
    SPI_p_train = signif(spi_p, 4),
    SARIMAX_order = sarimax_order_str,
    SARIMAX_AICc_train = round(informed$aicc, 2),
    Naive_AICc_train = round(naive$aicc, 2),        # <--- ADDED HERE
    Delta_AICc_train = round(delta_aicc, 2),
    LjungBox_p_lag12_train = round(lb_12, 4),
    ZINB_LR_train = round(zinb_lr, 3),
    ZINB_coef_train = round(zinb_coef, 4),
    ZINB_coef_p_train = signif(zinb_coef_p, 4)
  )
  
  candidate_results[[st]] <- list(climate_informed = informed$candidate_table,
                                  climate_naive = naive$candidate_table)
  
  cat("  Best Order identified from training:", sarimax_order_str, "\n")
  cat("  AICc (train):", round(informed$aicc, 2), "\n")
}

# 9. SAVE OUTPUTS
results <- bind_rows(results_list)
write_csv(results, file.path(TAB_DIR, "spi_table3_training_models_diagnostics.csv"))
saveRDS(results, file.path(OUT_DIR, "spi_training_models_results.rds"))
saveRDS(candidate_results, file.path(OUT_DIR, "spi_training_candidate_models.rds"))

cat("\nTraining-only model selection and coefficient analysis complete.\n")





