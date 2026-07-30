# ============================================================================
# 04_spi_sarimax_zinb_models.R
# Stationarity assessment (ADF + KPSS, both series), SARIMAX baseline with
# SPI-4 at the consensus lag as exogenous regressor, nested comparison
# against a climate-naive SARIMA, ZINB robustness check (Table 3), and
# expanded residual diagnostics.
#
# Stationarity (ADF + KPSS, plus ndiffs()/nsdiffs() as a cross-check) is
# assessed formally, before any SARIMAX/ZINB fitting, on both the Cases
# series and the lagged SPI-4 exogenous series for every state - this is
# what supports the d = D = 0 choice described in manuscript Section 2.6.
# ============================================================================

source("00_setup.R")

prepped <- readRDS(file.path(OUT_DIR, "prepared_data.rds"))
weekly  <- prepped$weekly

lag_results <- readRDS(file.path(OUT_DIR, "lag_tables_sensitivity.rds"))
lag_table   <- lag_results$interp$lag_table

SPI_SCALE   <- 4
clim_col    <- paste0("SPI_", SPI_SCALE)

# ---------------------------------------------------------------------------
# 0. Stationarity assessment (Section 2.6) - run BEFORE any SARIMAX/ZINB
#    fitting, on the exact series each model actually consumes:
#      - Cases:   the raw weekly count fed to Arima() and (rounded) to
#                 pscl::zeroinfl(); this is what "was it OK to fix d=D=0"
#                 is really asking about, so testing a transformed series
#                 instead (e.g. log1p) would not match what was fit.
#      - SPI_lag: the state-specific-lag exogenous regressor; SPI is
#                 constructed to be approximately N(0,1) by design (McKee
#                 et al. 1993), so this mainly serves as a manuscript-facing
#                 confirmation rather than an expected source of surprise.
#
# Two complementary tests are used because they have OPPOSITE null
# hypotheses, which is the standard way to guard against either test's
# specific blind spots:
#   - ADF (tseries::adf.test):  H0 = unit root present (non-stationary).
#                               Rejecting H0 (p < 0.05) supports stationarity.
#   - KPSS (tseries::kpss.test): H0 = level-stationary.
#                               Rejecting H0 (p < 0.05) supports NON-stationarity.
# Agreement between the two (ADF rejects + KPSS does not reject) is the
# strongest evidence of stationarity; disagreement is flagged as "Ambiguous"
# rather than silently resolved one way.
#
# forecast::ndiffs()/nsdiffs() are also reported alongside as a second,
# independent cross-check on the specific d=0/D=0 decision used throughout
# this script - if either suggests a nonzero order, that is a signal the
# fixed-differencing assumption may not hold for that state and is worth
# flagging in the manuscript rather than fitting through silently.
# ---------------------------------------------------------------------------
run_stationarity_checks <- function(x, series_label, state_label, freq = 52) {
  adf  <- tryCatch(tseries::adf.test(x, alternative = "stationary"), error = function(e) NULL)
  kpss <- tryCatch(tseries::kpss.test(x, null = "Level"), error = function(e) NULL)
  
  adf_stat  <- if (!is.null(adf))  unname(adf$statistic) else NA_real_
  adf_p     <- if (!is.null(adf))  adf$p.value           else NA_real_
  kpss_stat <- if (!is.null(kpss)) unname(kpss$statistic) else NA_real_
  kpss_p    <- if (!is.null(kpss)) kpss$p.value            else NA_real_
  
  suggested_d <- tryCatch(forecast::ndiffs(x), error = function(e) NA_integer_)
  suggested_D <- tryCatch(forecast::nsdiffs(ts(x, frequency = freq)), error = function(e) NA_integer_)
  
  conclusion <- dplyr::case_when(
    is.na(adf_p) | is.na(kpss_p) ~ "Test failed",
    adf_p < 0.05 & kpss_p >= 0.05 ~ "Stationary (ADF rejects unit root; KPSS does not reject stationarity)",
    adf_p >= 0.05 & kpss_p < 0.05 ~ "Non-stationary (ADF fails to reject unit root; KPSS rejects stationarity)",
    adf_p < 0.05 & kpss_p < 0.05 ~ "Ambiguous (ADF and KPSS disagree)",
    TRUE ~ "Ambiguous (neither test conclusive)"
  )
  
  tibble(
    State = state_label, Series = series_label,
    ADF_stat = round(adf_stat, 3), ADF_p = round(adf_p, 3),
    KPSS_stat = round(kpss_stat, 3), KPSS_p = round(kpss_p, 3),
    Suggested_d = suggested_d, Suggested_D = suggested_D,
    Conclusion = conclusion
  )
}

stationarity_results <- map_dfr(FOCAL_STATES, function(st) {
  df_state <- weekly %>% filter(State == st) %>% arrange(week_idx)
  
  cons_lag_st <- lag_table %>% filter(State == st) %>% pull(Consensus_lag_wk)
  spi_lag_series <- dplyr::lag(df_state[[clim_col]], cons_lag_st) %>% na.omit() %>% as.numeric()
  
  bind_rows(
    run_stationarity_checks(df_state$Cases, "Cases", st),
    run_stationarity_checks(spi_lag_series, "SPI_lag", st)
  )
})

print(stationarity_results, n = Inf, width = Inf)
write_csv(stationarity_results, file.path(TAB_DIR, "spi_supp_table_stationarity_checks.csv"))
cat("\nStationarity checks complete (ADF + KPSS + ndiffs/nsdiffs), saved to\n",
    file.path(TAB_DIR, "spi_supp_table_stationarity_checks.csv"),
    "\nReview the Conclusion/Suggested_d/Suggested_D columns before trusting the fixed\n",
    "d = D = 0 SARIMAX specification below for every state; if any state shows\n",
    "Suggested_d/Suggested_D > 0 or a 'Non-stationary'/'Ambiguous' conclusion for\n",
    "Cases, that is worth reporting as a caveat (or re-fitting that state's models\n",
    "with differencing) rather than silently proceeding with d = D = 0 everywhere.\n\n", sep = "")

# Grid search order (p,q in {0,1}; P,Q in {0,1}; d=D=0 fixed, per the
# stationarity checks above) - fit EXHAUSTIVELY over all 16 combinations,
# not via a stepwise or auto.arima() search; see Section 2.6.
ORDER_GRID <- expand.grid(p = 0:1, q = 0:1, P = 0:1, Q = 0:1)

fit_best_sarimax <- function(y, xreg = NULL, period = 52) {
  best <- list(aic = Inf, fit = NULL, order = NULL, seasonal = NULL)
  for (i in seq_len(nrow(ORDER_GRID))) {
    o <- ORDER_GRID[i, ]
    fit <- tryCatch(
      forecast::Arima(y, order = c(o$p, 0, o$q),
                      seasonal = list(order = c(o$P, 0, o$Q), period = period),
                      xreg = xreg, method = "CSS-ML"),
      error = function(e) NULL
    )
    if (!is.null(fit) && AIC(fit) < best$aic) {
      best <- list(aic = AIC(fit), fit = fit, order = c(o$p, 0, o$q),
                   seasonal = c(o$P, 0, o$Q))
    }
  }
  best
}

results <- map_dfr(FOCAL_STATES, function(st) {
  df_state <- weekly %>% filter(State == st) %>% arrange(week_idx)
  
  cons_lag <- lag_table %>% filter(State == st) %>% pull(Consensus_lag_wk)
  if (length(cons_lag) != 1) {
    stop(sprintf(
      "Expected exactly one Consensus_lag_wk for state '%s', got %d. ",
      st, length(cons_lag)))
  }
  
  y <- ts(df_state$Cases, frequency = 52)
  x_lagged <- dplyr::lag(df_state[[clim_col]], cons_lag)
  
  valid <- !is.na(x_lagged)
  y_trim <- y[valid]
  x_trim <- matrix(x_lagged[valid], ncol = 1, dimnames = list(NULL, "SPI_lag"))
  
  informed <- fit_best_sarimax(y_trim, xreg = x_trim)
  naive    <- fit_best_sarimax(y_trim, xreg = NULL)
  
  spi_beta <- if (!is.null(informed$fit)) informed$fit$coef[["SPI_lag"]] else NA_real_
  
  sarimax_order_str <- if (!is.null(informed$fit)) {
    paste0("(", paste(informed$order, collapse = ","), ")x(",
           paste(informed$seasonal, collapse = ","), ")52")
  } else {
    NA_character_
  }
  
  # Residual Diagnostics
  if (!is.null(informed$fit)) {
    resids <- residuals(informed$fit)
    resids <- resids[!is.na(resids)]
    
    n_params <- length(informed$fit$coef)
    lb_12 <- tryCatch(Box.test(resids, lag = 12, type = "Ljung-Box", fitdf = n_params)$p.value, error = function(e) NA_real_)
    lb_24 <- tryCatch(Box.test(resids, lag = 24, type = "Ljung-Box", fitdf = n_params)$p.value, error = function(e) NA_real_)
    
    resid_mean <- mean(resids)
    resid_sd   <- sd(resids)
  } else {
    lb_12 <- NA_real_
    lb_24 <- NA_real_
    resid_mean <- NA_real_
    resid_sd   <- NA_real_
  }
  
  # ZINB robustness check
  zinb_df <- tibble(cases = round(y_trim), spi = as.numeric(x_trim))
  
  zinb_fit <- tryCatch(
    pscl::zeroinfl(cases ~ spi | 1, data = zinb_df, dist = "negbin"),
    error = function(e) NULL, warning = function(w) NULL
  )
  zinb_naive_fit <- tryCatch(
    pscl::zeroinfl(cases ~ 1 | 1, data = zinb_df, dist = "negbin"),
    error = function(e) NULL, warning = function(w) NULL
  )
  
  if (is.null(zinb_fit) || is.null(zinb_naive_fit)) {
    lr_stat <- NA_real_; lr_p <- NA_real_
    zinb_coef_est <- NA_real_; zinb_coef_p <- NA_real_
  } else {
    lr <- lmtest::lrtest(zinb_naive_fit, zinb_fit)
    lr_stat <- lr$Chisq[2]
    lr_p    <- lr$`Pr(>Chisq)`[2]
    
    zinb_coef <- summary(zinb_fit)$coefficients$count["spi", ]
    zinb_coef_est <- zinb_coef["Estimate"]
    zinb_coef_p   <- zinb_coef["Pr(>|z|)"]
  }
  
  tibble(
    State                 = st,
    Consensus_lag_wk      = cons_lag,
    SPI_beta_SARIMAX      = round(spi_beta, 4),
    SARIMAX_order         = sarimax_order_str,
    SARIMAX_AIC           = round(informed$aic, 1),
    Naive_SARIMA_AIC      = round(naive$aic, 1),
    Delta_AIC             = round(naive$aic - informed$aic, 1),
    Residual_Mean         = round(resid_mean, 4),
    Residual_SD           = round(resid_sd, 4),
    LjungBox_p_lag12      = round(lb_12, 3),
    LjungBox_p_lag24      = round(lb_24, 3),
    ZINB_LR               = round(lr_stat, 2),
    ZINB_p                = round(lr_p, 3),
    ZINB_coef             = round(zinb_coef_est, 4),
    ZINB_coef_p           = round(zinb_coef_p, 3)
  )
})

print(results)
write_csv(results, file.path(TAB_DIR, "spi_table3_sarimax_zinb_diagnostics.csv"))
saveRDS(results, file.path(OUT_DIR, "spi_sarimax_zinb_results.rds"))
cat("SPI-based SARIMAX/ZINB modeling and residual diagnostics complete (Table 3).\n")
