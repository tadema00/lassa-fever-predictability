# ============================================================================
# 05_spi_rolling_forecast_alerts.R - CORRECTED FOR TEMPORAL INTEGRITY
# Primary rolling-origin evaluation and operational early-warning analysis.
# Evaluates 2024-2025 performance using pre-specified training-stage model 
# orders, consensus lags, and historical thresholds.
# ============================================================================

source("000_setup.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(readr)
  library(ggplot2)
  library(forecast)
  library(boot)
  library(patchwork)
})

# Conflict resolution
if (requireNamespace("conflicted", quietly = TRUE)) {
  suppressMessages({
    conflicted::conflicts_prefer(dplyr::filter, .quiet = TRUE)
    conflicted::conflicts_prefer(dplyr::select, .quiet = TRUE)
    conflicted::conflicts_prefer(dplyr::arrange, .quiet = TRUE)
    conflicted::conflicts_prefer(dplyr::lag, .quiet = TRUE)
  })
}

HAS_SURVEILLANCE <- requireNamespace("surveillance", quietly = TRUE)
if (!HAS_SURVEILLANCE) {
  warning("Package 'surveillance' not installed; surveillance algorithm comparators will be skipped.")
}

# ----------------------------------------------------------------------------
# 1. LOAD DATA AND TRAINING-STAGE DECISIONS
# ----------------------------------------------------------------------------
prepped <- readRDS(file.path(OUT_DIR, "prepared_data.rds"))
weekly  <- prepped$weekly

# CRITICAL: Load model orders and consensus lags identified during training phase (Script 04)
# Expected structure: tibble/data.frame with columns [State, Consensus_lag_wk, SARIMAX_order]
train_results <- readRDS(file.path(OUT_DIR, "spi_training_models_results.rds"))

# ----------------------------------------------------------------------------
# 2. GLOBAL PARAMETERS & SETUP
# ----------------------------------------------------------------------------
SPI_SCALE       <- 4
CLIM_COL        <- paste0("SPI_", SPI_SCALE)
SEASONAL_PERIOD <- 52
EVAL_YEAR       <- 2024
MAX_LEAD        <- 4
OUTBREAK_PCTILE <- 0.95
BOOT_R          <- 2000
MIN_VALID_OBS   <- 104  # Minimum valid observations required for fitting (2 years)

# CUSUM parameters (exploratory)
CUSUM_K <- 1
CUSUM_H <- 4

set.seed(SEED_GLOBAL)

# ----------------------------------------------------------------------------
# 3. HELPER FUNCTIONS
# ----------------------------------------------------------------------------

# Parse order string like "(1,0,1)(1,0,1)[52]" or "(1,0,0)x(1,1,1)52"
parse_order <- function(order_str) {
  if (is.na(order_str) || !is.character(order_str)) {
    return(list(order = c(1, 0, 1), seasonal = c(1, 0, 1)))
  }
  nums <- as.numeric(unlist(regmatches(order_str, gregexpr("[0-9]+", order_str))))
  if (length(nums) < 6) {
    stop(sprintf("Failed to parse SARIMAX order string: '%s'. Expected at least 6 integers.", order_str))
  }
  list(
    order    = c(nums[1], nums[2], nums[3]),
    seasonal = c(nums[4], nums[5], nums[6])
  )
}

calc_mae <- function(actual, predicted) {
  ok <- is.finite(actual) & is.finite(predicted)
  if (!any(ok)) return(NA_real_)
  mean(abs(actual[ok] - predicted[ok]))
}

calc_rmse <- function(actual, predicted) {
  ok <- is.finite(actual) & is.finite(predicted)
  if (!any(ok)) return(NA_real_)
  sqrt(mean((actual[ok] - predicted[ok])^2))
}

calc_mape <- function(actual, predicted) {
  ok <- is.finite(actual) & is.finite(predicted) & actual != 0
  if (!any(ok)) return(NA_real_)
  mean(abs((actual[ok] - predicted[ok]) / actual[ok])) * 100
}

calc_mase <- function(actual, predicted, train_y, season = 52) {
  ok <- is.finite(actual) & is.finite(predicted)
  if (!any(ok)) return(NA_real_)
  
  # Mean absolute error of in-sample seasonal naive forecast
  if (length(train_y) <= season) return(NA_real_)
  scale <- mean(abs(diff(train_y, lag = season)), na.rm = TRUE)
  if (is.na(scale) || scale == 0) return(NA_real_)
  
  mean(abs(actual[ok] - predicted[ok])) / scale
}

classify_alert <- function(predicted, thresholds) {
  out <- rep(NA_character_, length(predicted))
  ok  <- is.finite(predicted)
  out[ok & predicted < thresholds[["50%"]]] <- "Green"
  out[ok & predicted >= thresholds[["50%"]] & predicted < thresholds[["75%"]]] <- "Yellow"
  out[ok & predicted >= thresholds[["75%"]] & predicted <= thresholds[["95%"]]] <- "Orange"
  out[ok & predicted > thresholds[["95%"]]] <- "Red"
  out
}

binary_metrics <- function(actual_red, predicted_red) {
  ok <- !is.na(actual_red) & !is.na(predicted_red)
  actual_red <- actual_red[ok]
  predicted_red <- predicted_red[ok]
  
  TP <- sum(actual_red & predicted_red)
  FP <- sum(!actual_red & predicted_red)
  FN <- sum(actual_red & !predicted_red)
  TN <- sum(!actual_red & !predicted_red)
  
  sensitivity <- if ((TP + FN) > 0) TP / (TP + FN) else NA_real_
  specificity <- if ((TN + FP) > 0) TN / (TN + FP) else NA_real_
  ppv         <- if ((TP + FP) > 0) TP / (TP + FP) else NA_real_
  npv         <- if ((TN + FN) > 0) TN / (TN + FN) else NA_real_
  accuracy    <- if ((TP + FP + FN + TN) > 0) (TP + TN) / (TP + FP + FN + TN) else NA_real_
  
  tibble(TP = TP, FP = FP, FN = FN, TN = TN,
         sensitivity = sensitivity, specificity = specificity,
         ppv = ppv, npv = npv, accuracy = accuracy)
}

bootstrap_alert_metrics <- function(actual, predicted, red_threshold,
                                    R = BOOT_R, block_length = 8) {
  ok <- is.finite(actual) & is.finite(predicted)
  actual <- actual[ok]; predicted <- predicted[ok]
  n <- length(actual)
  if (n == 0) {
    return(tibble(sensitivity = NA_real_, sens_lo = NA_real_, sens_hi = NA_real_,
                  ppv = NA_real_, ppv_lo = NA_real_, ppv_hi = NA_real_,
                  specificity = NA_real_, spec_lo = NA_real_, spec_hi = NA_real_,
                  n_actual_red = 0L, n_predicted_red = 0L,
                  status = "No valid forecasts"))
  }
  
  actual_red <- actual > red_threshold
  predicted_red <- predicted > red_threshold
  n_actual_red <- sum(actual_red)
  n_predicted_red <- sum(predicted_red)
  
  point <- binary_metrics(actual_red, predicted_red)
  sensitivity <- point$sensitivity; ppv <- point$ppv; specificity <- point$specificity
  
  status <- dplyr::case_when(
    n_predicted_red == 0 & n_actual_red == 0 ~ "No actual or predicted Red alerts",
    n_predicted_red == 0 ~ "No predicted Red alerts; PPV undefined",
    n_actual_red == 0 ~ "No observed Red events; sensitivity undefined",
    TRUE ~ "Estimable"
  )
  
  block_length <- max(1L, min(as.integer(block_length), n))
  n_blocks <- ceiling(n / block_length)
  
  boot_statistics <- matrix(NA_real_, nrow = R, ncol = 3)
  colnames(boot_statistics) <- c("sensitivity", "ppv", "specificity")
  
  for (r in seq_len(R)) {
    starts <- sample(seq_len(n - block_length + 1), size = n_blocks, replace = TRUE)
    indices <- unlist(lapply(starts, function(s) s:(s + block_length - 1)))
    indices <- indices[seq_len(n)]
    
    a <- actual[indices]; p <- predicted[indices]
    a_red <- a > red_threshold; p_red <- p > red_threshold
    m <- binary_metrics(a_red, p_red)
    boot_statistics[r, ] <- c(m$sensitivity, m$ppv, m$specificity)
  }
  
  safe_percentile_ci <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) < 20) return(c(NA_real_, NA_real_))
    as.numeric(quantile(x, probs = c(0.025, 0.975), na.rm = TRUE, names = FALSE, type = 6))
  }
  
  sens_ci <- safe_percentile_ci(boot_statistics[, "sensitivity"])
  ppv_ci  <- safe_percentile_ci(boot_statistics[, "ppv"])
  spec_ci <- safe_percentile_ci(boot_statistics[, "specificity"])
  
  tibble(
    sensitivity = sensitivity, sens_lo = sens_ci[1], sens_hi = sens_ci[2],
    ppv = ppv, ppv_lo = ppv_ci[1], ppv_hi = ppv_ci[2],
    specificity = specificity, spec_lo = spec_ci[1], spec_hi = spec_ci[2],
    n_actual_red = n_actual_red, n_predicted_red = n_predicted_red,
    status = status
  )
}

identify_episodes <- function(x, threshold) {
  exceed <- is.finite(x) & x > threshold
  if (!any(exceed)) {
    return(tibble(
      episode = integer(),
      start = integer(),
      end = integer(),
      duration = integer()
    ))
  }
  starts <- which(exceed & !dplyr::lag(exceed, default = FALSE))
  ends   <- which(exceed & !dplyr::lead(exceed, default = FALSE))
  tibble(
    episode = seq_along(starts),
    start   = starts,
    end     = ends,
    duration = ends - starts + 1
  )
}

# ----------------------------------------------------------------------------
# 4. ROLLING-ORIGIN EVALUATION FUNCTION (STRICT TEMPORAL INTEGRITY)
# ----------------------------------------------------------------------------
rolling_origin_eval_robust <- function(df_state, st_name, training_info, 
                                       season = SEASONAL_PERIOD, 
                                       outbreak_pctile = OUTBREAK_PCTILE) {
  df_state <- df_state %>% dplyr::arrange(week_idx)
  
  # A. Identify Training/Validation Split Point
  eval_candidates <- which(df_state$Year == EVAL_YEAR & df_state$Epi_Week == 1)
  eval_start_idx  <- if (length(eval_candidates) > 0) eval_candidates[1] else max(1, nrow(df_state) - 96)
  
  if (eval_start_idx < 2) stop("eval_start_idx must be >= 2")
  if (eval_start_idx >= nrow(df_state)) stop(sprintf("Invalid evaluation start index for state %s.", st_name))
  
  # B. Get the pre-determined fixed Lag and SARIMA Order from training step
  st_train <- training_info %>% dplyr::filter(State == st_name)
  if (nrow(st_train) != 1) {
    stop(sprintf("Expected exactly 1 training result record for state '%s', found %d.", st_name, nrow(st_train)))
  }
  
  cons_lag <- as.integer(st_train$Consensus_lag_wk)
  ord_list <- parse_order(st_train$SARIMAX_order)
  
  # C. Outbreak Threshold fixed from training period only
  hist_cases     <- df_state$Cases[1:(eval_start_idx - 1)]
  hist_threshold <- as.numeric(quantile(hist_cases, probs = outbreak_pctile, na.rm = TRUE, names = FALSE, type = 7))
  
  # D. Construct climate covariate with appropriate lag
  x_lag <- dplyr::lag(df_state[[CLIM_COL]], n = cons_lag)
  
  eval_idx <- seq(eval_start_idx, nrow(df_state))
  n_eval   <- length(eval_idx)
  actual   <- df_state$Cases[eval_idx]
  
  preds_clim  <- rep(NA_real_, n_eval)
  preds_naive <- rep(NA_real_, n_eval)
  
  # E. Rolling Loop (Refitting Parameters, Holding Model Order Fixed)
  for (i in seq_along(eval_idx)) {
    t <- eval_idx[i]
    
    train_y <- df_state$Cases[1:(t - 1)]
    train_x <- x_lag[1:(t - 1)]
    valid   <- is.finite(train_y) & is.finite(train_x)
    
    if (sum(valid) >= MIN_VALID_OBS) {
      train_y_ts  <- ts(train_y[valid], frequency = season)
      train_x_mat <- matrix(train_x[valid], ncol = 1)
      colnames(train_x_mat) <- "SPI_lag"
      
      # Fit model using fixed order established during training
      fit <- tryCatch({
        forecast::Arima(
          y = train_y_ts,
          order = ord_list$order,
          seasonal = list(order = ord_list$seasonal, period = season),
          xreg = train_x_mat,
          method = "ML"
        )
      }, error = function(e) {
        # Fallback to CSS-ML if ML fails
        tryCatch({
          forecast::Arima(
            y = train_y_ts,
            order = ord_list$order,
            seasonal = list(order = ord_list$seasonal, period = season),
            xreg = train_x_mat,
            method = "CSS-ML"
          )
        }, error = function(e2) NULL)
      })
      
      if (!is.null(fit) && is.finite(x_lag[t])) {
        xreg_new <- matrix(x_lag[t], nrow = 1, ncol = 1)
        colnames(xreg_new) <- "SPI_lag"
        fc <- tryCatch(forecast::forecast(fit, h = 1, xreg = xreg_new), error = function(e) NULL)
        if (!is.null(fc)) preds_clim[i] <- as.numeric(fc$mean[1])
      }
    } else {
      warning(sprintf("%s: Insufficient training observations at t=%d.", st_name, t))
    }
    
    # Seasonal Naive benchmark (t - 52)
    if (t - season > 0) {
      preds_naive[i] <- df_state$Cases[t - season]
    }
  }
  
  # Error Metrics
  mae_clim   <- calc_mae(actual, preds_clim)
  mae_naive  <- calc_mae(actual, preds_naive)
  rmse_clim  <- calc_rmse(actual, preds_clim)
  rmse_naive <- calc_rmse(actual, preds_naive)
  mape_clim  <- calc_mape(actual, preds_clim)
  mape_naive <- calc_mape(actual, preds_naive)
  
  mase_clim  <- calc_mase(actual, preds_clim, hist_cases, season = season)
  mase_naive <- calc_mase(actual, preds_naive, hist_cases, season = season)
  
  # Episode-Based Early Warning Alerts Evaluation
  alert_ok     <- is.finite(actual) & is.finite(preds_clim)
  alert_actual <- actual[alert_ok]
  alert_pred   <- preds_clim[alert_ok]
  
  predicted_red <- alert_pred > hist_threshold
  episodes      <- identify_episodes(alert_actual, hist_threshold)
  n_episodes    <- nrow(episodes)
  
  detected   <- logical(n_episodes)
  lead_times <- numeric(0)
  
  if (n_episodes > 0) {
    for (j in seq_len(n_episodes)) {
      onset <- episodes$start[j]
      if (onset <= 1) {
        detected[j] <- FALSE
        next
      }
      start_w <- max(1, onset - MAX_LEAD)
      end_w   <- onset - 1
      if (start_w > end_w) {
        detected[j] <- FALSE
        next
      }
      window_alerts <- which(predicted_red[start_w:end_w]) + start_w - 1
      if (length(window_alerts) > 0) {
        detected[j] <- TRUE
        lead_times  <- c(lead_times, onset - max(window_alerts))
      } else {
        detected[j] <- FALSE
      }
    }
  }
  
  n_detected     <- sum(detected, na.rm = TRUE)
  detection_rate <- if (n_episodes > 0) n_detected / n_episodes else NA_real_
  mean_lead      <- if (length(lead_times) > 0) mean(lead_times, na.rm = TRUE) else NA_real_
  median_lead    <- if (length(lead_times) > 0) median(lead_times, na.rm = TRUE) else NA_real_
  q1_lead        <- if (length(lead_times) > 0) quantile(lead_times, 0.25, na.rm = TRUE) else NA_real_
  q3_lead        <- if (length(lead_times) > 0) quantile(lead_times, 0.75, na.rm = TRUE) else NA_real_
  
  list(
    mae_clim = mae_clim, mae_naive = mae_naive,
    rmse_clim = rmse_clim, rmse_naive = rmse_naive,
    mape_clim = mape_clim, mape_naive = mape_naive,
    mase_clim = mase_clim, mase_naive = mase_naive,
    n_episodes = n_episodes, n_detected = n_detected,
    detection_rate = detection_rate,
    mean_lead = mean_lead, median_lead = median_lead,
    q1_lead = q1_lead, q3_lead = q3_lead,
    actual = actual, preds_clim = preds_clim, preds_naiveS = preds_naive,
    alert_actual = alert_actual, alert_pred = alert_pred,
    predicted_red = predicted_red,
    eval_idx = eval_idx, hist_threshold = hist_threshold,
    lead_times = lead_times, episodes = episodes
  )
}

# ============================================================================
# 5. RUN FOR ALL STATES & GENERATE TABLE 4
# ============================================================================
eval_results <- list()

table4 <- purrr::map_dfr(FOCAL_STATES, function(st) {
  message("\nExecuting rolling-origin forecast for state: ", st)
  df_state <- weekly %>% dplyr::filter(State == st) %>% dplyr::arrange(week_idx)
  
  res <- rolling_origin_eval_robust(
    df_state = df_state, 
    st_name = st, 
    training_info = train_results
  )
  
  eval_results[[st]] <<- res
  
  tibble(
    State              = st,
    MAE_climate        = round(res$mae_clim, 2),
    MAE_naive          = round(res$mae_naive, 2),
    RMSE_climate       = round(res$rmse_clim, 2),
    RMSE_naive         = round(res$rmse_naive, 2),
    MASE_climate       = round(res$mase_clim, 3),
    MASE_naive         = round(res$mase_naive, 3),
    Outbreaks_detected = sprintf("%d/%d", res$n_detected, res$n_episodes),
    Detection_rate     = ifelse(is.na(res$detection_rate), "\u2014", sprintf("%.0f%%", 100 * res$detection_rate)),
    Lead_weeks         = ifelse(is.na(res$mean_lead), "\u2014", as.character(round(res$mean_lead, 1)))
  )
}) %>%
  dplyr::mutate(
    MAE_skill  = round((MAE_naive - MAE_climate) / MAE_naive * 100, 1),
    RMSE_skill = round((RMSE_naive - RMSE_climate) / RMSE_naive * 100, 1)
  ) %>%
  dplyr::select(
    State, MAE_climate, MAE_naive, MAE_skill,
    RMSE_climate, RMSE_naive, RMSE_skill,
    MASE_climate, MASE_naive,
    Outbreaks_detected, Detection_rate, Lead_weeks
  )

print(table4)
write_csv(table4, file.path(TAB_DIR, "spi_table4_forecast_evaluation.csv"))
saveRDS(eval_results, file.path(OUT_DIR, "spi_rolling_eval_results.rds"))

# ============================================================================
# 6. FORECAST PANEL FIGURES
# ============================================================================
plot_forecast_panel <- function(states, panel_letters, filename) {
  panels <- lapply(seq_along(states), function(i) {
    st  <- states[i]
    res <- eval_results[[st]]
    df  <- tibble(
      week             = res$eval_idx,
      Observed         = res$actual,
      Climate_informed = res$preds_clim,
      Seasonal_naive   = res$preds_naiveS
    ) %>% tidyr::pivot_longer(-week, names_to = "series", values_to = "value")
    
    ggplot(df, aes(x = week, y = value, color = series, linetype = series)) +
      geom_line(linewidth = 0.6, na.rm = TRUE) +
      labs(title = sprintf("(%s) %s", panel_letters[i], st),
           x = "Evaluation week", y = "Weekly cases") +
      theme_minimal(base_size = 10) +
      theme(legend.position = if (i == length(states)) "right" else "none")
  })
  
  combined <- patchwork::wrap_plots(panels, ncol = length(states))
  ggsave(file.path(FIG_DIR, filename), combined, width = 13, height = 3.2, dpi = 300, bg = "white")
  cat("Saved figure panel to:", file.path(FIG_DIR, filename), "\n")
}

plot_forecast_panel(c("Ondo", "Ebonyi", "Edo"), c("A", "B", "C"), "fig6_spi_rolling_forecast_southern.png")
plot_forecast_panel(c("Bauchi", "Taraba", "Plateau"), c("A", "B", "C"), "fig7_spi_rolling_forecast_northern.png")

# ============================================================================
# 7. ALERT VALIDATION (TABLE 5)
# ============================================================================
table5_spi <- purrr::map_dfr(FOCAL_STATES, function(st) {
  message("Calculating alert validation metrics for: ", st)
  df_state <- weekly %>% dplyr::filter(State == st) %>% dplyr::arrange(week_idx)
  res      <- eval_results[[st]]
  
  hist_end   <- min(res$eval_idx) - 1
  hist_cases <- df_state$Cases[seq_len(hist_end)]
  
  hist_thresholds <- quantile(hist_cases, probs = c(0.50, 0.75, 0.95), na.rm = TRUE, names = FALSE)
  names(hist_thresholds) <- c("50%", "75%", "95%")
  
  alert_class <- classify_alert(res$preds_clim, hist_thresholds)
  bt          <- bootstrap_alert_metrics(actual = res$actual, predicted = res$preds_clim,
                                         red_threshold = hist_thresholds[["95%"]], R = BOOT_R)
  
  lead_text <- if (is.na(res$median_lead)) {
    "\u2014"
  } else {
    sprintf("%.1f (%.1f\u2013%.1f)", res$median_lead, res$q1_lead, res$q3_lead)
  }
  
  tibble(
    State            = st,
    threshold_50     = hist_thresholds[["50%"]],
    threshold_75     = hist_thresholds[["75%"]],
    threshold_95     = hist_thresholds[["95%"]],
    sensitivity      = bt$sensitivity,
    sens_lo          = bt$sens_lo, sens_hi = bt$sens_hi,
    ppv              = bt$ppv, ppv_lo = bt$ppv_lo, ppv_hi = bt$ppv_hi,
    specificity      = bt$specificity,
    spec_lo          = bt$spec_lo, spec_hi = bt$spec_hi,
    n_actual_red     = bt$n_actual_red,
    n_predicted_red  = bt$n_predicted_red,
    alert_status     = bt$status,
    Lead_median_IQR  = lead_text,
    green_n          = sum(alert_class == "Green", na.rm = TRUE),
    yellow_n         = sum(alert_class == "Yellow", na.rm = TRUE),
    orange_n         = sum(alert_class == "Orange", na.rm = TRUE),
    red_n            = sum(alert_class == "Red", na.rm = TRUE)
  )
})

print(table5_spi)
write_csv(table5_spi, file.path(TAB_DIR, "spi_table5_alert_validation.csv"))

# ============================================================================
# 8. OPTIONAL SURVEILLANCE COMPARATORS
# ============================================================================
if (HAS_SURVEILLANCE) {
  message("\nRunning surveillance-package algorithmic comparators...")
  surveillance_results <- list()
  farrington_warnings  <- list()
  
  for (st in FOCAL_STATES) {
    message("Evaluating surveillance algorithms for: ", st)
    df_state   <- weekly %>% dplyr::filter(State == st) %>% dplyr::arrange(week_idx)
    start_year <- min(df_state$Year, na.rm = TRUE)
    start_week <- min(df_state$Epi_Week[df_state$Year == start_year], na.rm = TRUE)
    
    sts_data <- surveillance::sts(
      observed  = df_state$Cases,
      start     = c(start_year, start_week),
      frequency = 52
    )
    
    eval_start <- which(df_state$Year == EVAL_YEAR & df_state$Epi_Week == 1)[1]
    if (is.na(eval_start)) { 
      warning("No evaluation year week 1 found for ", st); next 
    }
    eval_range <- eval_start:nrow(df_state)
    
    # Farrington Flexible
    farrington_control <- list(range = eval_range, b = 3, w = 2, reweight = TRUE, trend = TRUE, alpha = 0.05)
    farrington_fit     <- tryCatch(
      withCallingHandlers(
        surveillance::farringtonFlexible(sts_data, control = farrington_control),
        warning = function(w) {
          farrington_warnings[[st]] <<- c(farrington_warnings[[st]], conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      ),
      error = function(e) {
        warning("Farrington failed for ", st, ": ", e$message)
        NULL
      }
    )
    farrington_alarm <- if (!is.null(farrington_fit)) surveillance::alarms(farrington_fit) else rep(NA, nrow(df_state))
    
    # EARS-C2
    ears_control <- list(range = eval_range, method = "C2", alpha = 0.05)
    ears_fit     <- tryCatch(surveillance::earsC(sts_data, control = ears_control),
                             error = function(e) { warning("EARS failed for ", st, ": ", e$message); NULL })
    ears_alarm   <- if (!is.null(ears_fit)) surveillance::alarms(ears_fit) else rep(NA, nrow(df_state))
    
    # CUSUM
    cusum_control <- list(range = eval_range, k = CUSUM_K, h = CUSUM_H, trans = "rossi")
    cusum_fit     <- tryCatch(surveillance::cusum(sts_data, control = cusum_control),
                              error = function(e) { warning("CUSUM failed for ", st, ": ", e$message); NULL })
    cusum_alarm   <- if (!is.null(cusum_fit)) surveillance::alarms(cusum_fit) else rep(NA, nrow(df_state))
    
    surveillance_results[[st]] <- list(
      farrington = farrington_alarm,
      ears_c2    = ears_alarm,
      cusum      = cusum_alarm
    )
  }
  
  saveRDS(surveillance_results, file.path(OUT_DIR, "surveillance_algorithm_results.rds"))
  
  comparator_table <- purrr::map_dfr(FOCAL_STATES, function(st) {
    x <- surveillance_results[[st]]
    tibble(
      State             = st,
      Farrington_alerts = sum(x$farrington, na.rm = TRUE),
      EARS_C2_alerts    = sum(x$ears_c2, na.rm = TRUE),
      CUSUM_alerts      = sum(x$cusum, na.rm = TRUE)
    )
  })
  
  print(comparator_table)
  write_csv(comparator_table, file.path(TAB_DIR, "spi_surveillance_algorithm_comparison.csv"))
}

# ============================================================================
# 9. DIAGNOSTIC REPORT & COMPLETE OUTPUT SAVE
# ============================================================================
diagnostic_table <- purrr::map_dfr(FOCAL_STATES, function(st) {
  res <- eval_results[[st]]
  tibble(
    State             = st,
    Evaluation_weeks  = length(res$actual),
    Valid_forecasts   = sum(is.finite(res$preds_clim)),
    Missing_forecasts = sum(!is.finite(res$preds_clim)),
    Observed_episodes = res$n_episodes,
    Detected_episodes = res$n_detected,
    Detection_rate    = ifelse(is.na(res$detection_rate), NA_real_, res$detection_rate),
    Median_lead       = res$median_lead,
    Hist_95_threshold = res$hist_threshold
  )
})

print(diagnostic_table)
write_csv(diagnostic_table, file.path(TAB_DIR, "spi_alert_diagnostic_summary.csv"))

saveRDS(list(table4 = table4, table5_spi = table5_spi,
             eval_results = eval_results, diagnostics = diagnostic_table),
        file.path(OUT_DIR, "spi_complete_forecast_alert_analysis.rds"))

cat("\n============================================================\n",
    "SPI-based rolling-origin forecast and operational alert\n",
    "analysis completed successfully with strict temporal integrity.\n",
    "============================================================\n")
