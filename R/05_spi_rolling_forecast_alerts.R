# ============================================================================
# 05_spi_rolling_forecast_alerts.R
# Expanding-window rolling-origin forecast evaluation using SPI-4 as the
# exogenous climatic driver (Section 2.5; Table 4, Figures 6-7) and the
# four-tier operational alert system with bootstrap sensitivity/PPV
# validation (Section 2.6; Table 5).
# ============================================================================

source("00_setup.R")

prepped <- readRDS(file.path(OUT_DIR, "prepared_data.rds"))
weekly  <- prepped$weekly

lag_results <- readRDS(file.path(OUT_DIR, "lag_tables_sensitivity.rds"))
lag_table   <- lag_results$interp$lag_table

SPI_SCALE           <- 4
clim_col            <- paste0("SPI_", SPI_SCALE)
EVAL_START_WEEK_IDX <- NULL
REFIT_EVERY         <- 12
OUTBREAK_PCTILE     <- 0.95
MAX_LEAD            <- 4

# ---------------------------------------------------------------------------
# Rolling-origin one-step evaluation for a single state using SPI
# ---------------------------------------------------------------------------
rolling_origin_eval <- function(df_state, cons_lag, eval_start_idx, season = 52,
                                outbreak_pctile = OUTBREAK_PCTILE) {
  
  x_lag <- dplyr::lag(df_state[[clim_col]], cons_lag)
  train_hist_threshold <- quantile(df_state$Cases[1:eval_start_idx], outbreak_pctile, na.rm = TRUE)
  
  eval_idx <- seq(eval_start_idx, nrow(df_state) - 1)
  preds_clim   <- numeric(length(eval_idx))
  preds_naiveS <- numeric(length(eval_idx))
  refit_fit <- NULL
  
  for (i in seq_along(eval_idx)) {
    t <- eval_idx[i]
    if ((i - 1) %% REFIT_EVERY == 0 || is.null(refit_fit)) {
      train_y <- ts(df_state$Cases[1:t], frequency = season)
      train_x <- matrix(x_lag[1:t], ncol = 1, dimnames = list(NULL, "SPI_lag"))
      valid <- !is.na(train_x[, 1])
      refit_fit <- tryCatch(
        forecast::Arima(train_y[valid], order = c(1, 0, 1),
                        seasonal = list(order = c(1, 0, 1), period = season),
                        xreg = train_x[valid, , drop = FALSE]),
        error = function(e) NULL
      )
    }
    if (!is.null(refit_fit) && !is.na(x_lag[t + 1])) {
      fc <- tryCatch(
        forecast::forecast(refit_fit, h = 1, xreg = matrix(x_lag[t + 1], ncol = 1,
                                                           dimnames = list(NULL, "SPI_lag"))),
        error = function(e) NULL
      )
      preds_clim[i] <- if (!is.null(fc)) as.numeric(fc$mean) else NA
    } else {
      preds_clim[i] <- NA
    }
    preds_naiveS[i] <- if (t - season + 1 > 0) df_state$Cases[t - season + 1] else NA
  }
  
  actual <- df_state$Cases[eval_idx + 1]
  train_actual <- df_state$Cases[1:eval_start_idx]
  
  mase_clim  <- mase(actual, preds_clim, train_actual, season)
  mase_naive <- mase(actual, preds_naiveS, train_actual, season)
  
  outbreak_weeks <- which(actual > train_hist_threshold)
  detected <- sapply(outbreak_weeks, function(w) {
    lookback <- max(1, w - MAX_LEAD):max(1, w - 1)
    any(preds_clim[lookback] > train_hist_threshold, na.rm = TRUE)
  })
  lead_times <- sapply(outbreak_weeks[detected], function(w) {
    lookback <- max(1, w - MAX_LEAD):max(1, w - 1)
    hit <- lookback[which(preds_clim[lookback] > train_hist_threshold)]
    if (length(hit) == 0) return(NA)
    w - max(hit)
  })
  
  list(
    mase_clim = mase_clim, mase_naive = mase_naive,
    n_outbreaks = length(outbreak_weeks), n_detected = sum(detected),
    detection_rate = mean(detected), mean_lead = mean(lead_times, na.rm = TRUE),
    actual = actual, preds_clim = preds_clim, preds_naiveS = preds_naiveS,
    eval_idx = eval_idx, hist_threshold = train_hist_threshold
  )
}

eval_results <- list()
table4 <- map_dfr(FOCAL_STATES, function(st) {
  df_state <- weekly %>% filter(State == st) %>% arrange(week_idx)
  
  cons_lag <- lag_table %>% filter(State == st) %>% pull(Consensus_lag_wk)
  if (length(cons_lag) != 1) {
    stop(sprintf(
      "Expected exactly one Consensus_lag_wk for state '%s', got %d. ",
      st, length(cons_lag)))
  }
  
  eval_start_idx <- df_state %>% filter(Year == 2024, Epi_Week == 1) %>% pull(week_idx)
  if (length(eval_start_idx) == 0) {
    eval_start_idx <- nrow(df_state) - 96
  } else if (length(eval_start_idx) > 1) {
    warning(sprintf(
      "State '%s': found %d rows matching Year==2024 & Epi_Week==1. Using the first.",
      st, length(eval_start_idx)))
    eval_start_idx <- eval_start_idx[1]
  }
  
  res <- rolling_origin_eval(df_state, cons_lag, eval_start_idx)
  eval_results[[st]] <<- res
  
  tibble(
    State = st,
    MASE_climate = round(res$mase_clim, 3),
    MASE_naive   = round(res$mase_naive, 3),
    Outbreaks_detected = sprintf("%d/%d", res$n_detected, res$n_outbreaks),
    Detection_rate = sprintf("%.0f%%", 100 * res$detection_rate),
    Lead_weeks = round(res$mean_lead, 1)
  )
})

print(table4)
write_csv(table4, file.path(TAB_DIR, "spi_table4_forecast_evaluation.csv"))
saveRDS(eval_results, file.path(OUT_DIR, "spi_rolling_eval_results.rds"))

# ---------------------------------------------------------------------------
# Figures 6 & 7: rolling-origin forecast plots by cluster
# ---------------------------------------------------------------------------
plot_forecast_panel <- function(states, panel_letters, filename) {
  panels <- lapply(seq_along(states), function(i) {
    st <- states[i]
    res <- eval_results[[st]]
    df <- tibble(
      week = res$eval_idx + 1,
      Observed = res$actual,
      Climate_informed = res$preds_clim,
      Seasonal_naive = res$preds_naiveS
    ) %>% pivot_longer(-week, names_to = "series", values_to = "value")
    
    ggplot(df, aes(x = week, y = value, color = series, linetype = series)) +
      geom_line(linewidth = 0.6) +
      scale_color_manual(values = c(Observed = "black", Climate_informed = "red",
                                    Seasonal_naive = "blue"), name = NULL) +
      scale_linetype_manual(values = c(Observed = "solid", Climate_informed = "solid",
                                       Seasonal_naive = "dashed"), name = NULL) +
      labs(title = sprintf("(%s) %s", panel_letters[i], st), x = "Week", y = "Weekly cases") +
      theme_minimal(base_size = 10) +
      theme(legend.position = if (i == length(states)) "right" else "none")
  })
  combined <- patchwork::wrap_plots(panels, ncol = length(states))
  ggsave(file.path(FIG_DIR, filename), combined, width = 13, height = 3.2, dpi = 220, bg = "white")
  cat("Saved", file.path(FIG_DIR, filename), "\n")
}

plot_forecast_panel(c("Ondo", "Ebonyi", "Edo"), c("A", "B", "C"),
                    "fig6_spi_rolling_forecast_southern.png")
plot_forecast_panel(c("Bauchi", "Taraba", "Plateau"), c("A", "B", "C"),
                    "fig7_spi_rolling_forecast_northern.png")

# ---------------------------------------------------------------------------
# Four-tier operational alert system (Section 2.6): Green (<50th),
# Yellow (50th-75th), Orange (75th-95th), Red (>95th), bootstrapped
# sensitivity/PPV for the Red tier (Table 5)
# ---------------------------------------------------------------------------
alert_bootstrap <- function(actual, predicted, hist_thresholds, n_boot = 1000) {
  ok <- !is.na(predicted) & !is.na(actual)
  actual <- actual[ok]; predicted <- predicted[ok]
  n <- length(actual)
  if (n == 0) return(tibble(sensitivity = NA, sens_lo = NA, sens_hi = NA,
                            ppv = NA, ppv_lo = NA, ppv_hi = NA))
  
  red_thresh <- hist_thresholds[["95%"]]
  
  boot_stat <- function(data, idx) {
    a <- data$actual[idx]; p <- data$predicted[idx]
    true_red  <- a > red_thresh
    pred_red  <- p > red_thresh
    sens <- if (sum(true_red) > 0) sum(pred_red & true_red) / sum(true_red) else NA
    ppv  <- if (sum(pred_red) > 0) sum(pred_red & true_red) / sum(pred_red) else NA
    c(sens = sens, ppv = ppv)
  }
  
  bd <- data.frame(actual = actual, predicted = predicted)
  b <- boot::boot(bd, boot_stat, R = n_boot)
  
  sens_ci <- tryCatch(boot::boot.ci(b, type = "perc", index = 1)$percent[4:5], error = function(e) c(NA, NA))
  ppv_ci  <- tryCatch(boot::boot.ci(b, type = "perc", index = 2)$percent[4:5], error = function(e) c(NA, NA))
  
  tibble(
    sensitivity = b$t0["sens"], sens_lo = sens_ci[1], sens_hi = sens_ci[2],
    ppv = b$t0["ppv"], ppv_lo = ppv_ci[1], ppv_hi = ppv_ci[2]
  )
}

# NOTE ON REPRODUCIBILITY: `boot::boot()` inside alert_bootstrap() is the only
# stochastic step in this script (rolling_origin_eval()/Arima/forecast are all
# deterministic given the data, so table4 needs no seed). The seed must be set
# immediately before the call that actually consumes randomness - table5's
# map_dfr() below - not before alert_bootstrap()'s definition a few lines up.
# Seeding before a function DEFINITION does nothing for calls made later,
# because any other code executed between the seed and the call advances the
# RNG stream first; this is a common and easy-to-miss source of "reproducible
# in theory, not reproducible in practice" results.
set.seed(SEED_GLOBAL)
table5 <- map_dfr(FOCAL_STATES, function(st) {
  df_state <- weekly %>% filter(State == st) %>% arrange(week_idx)
  res <- eval_results[[st]]
  hist_thresholds <- quantile(df_state$Cases[1:min(res$eval_idx)], c(0.5, 0.75, 0.95),
                              na.rm = TRUE) %>% setNames(c("50%", "75%", "95%"))
  
  bt <- alert_bootstrap(res$actual, res$preds_clim, hist_thresholds)
  bt %>% mutate(State = st, .before = 1)
})

print(table5)
write_csv(table5, file.path(TAB_DIR, "spi_table5_alert_validation.csv"))

cat("SPI-based rolling-origin forecast evaluation, figures, and alert validation complete.\n")

