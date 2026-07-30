# ============================================================================
# 05b_spi_alert_sensitivity_PI.R
# SUPPLEMENTARY SENSITIVITY ANALYSIS - does the Table 5 conclusion depend on
# using the SARIMAX *point forecast mean* as the Red-tier alert trigger?
#
# Rationale: the primary analysis (05_spi_rolling_forecast_alerts.R) fires a
# Red alert only when the forecast MEAN exceeds the historical 95th-percentile
# outbreak threshold - a conservative test, since an outbreak is by
# definition an extreme deviation and a conditional-mean forecast is
# structurally bad at reaching extremes even when the model has captured
# real signal. This script repeats the identical rolling-origin evaluation
# but also triggers alerts off the forecast's 80% and 95% UPPER prediction
# limits, so forecast uncertainty is taken into account.
#
# This is a companion/supplementary analysis, not a replacement:
#   - Primary analysis (unchanged):  05_spi_rolling_forecast_alerts.R -> Table 5
#   - Sensitivity analysis (this):   05b_spi_alert_sensitivity_PI.R  -> Supp. Table/Figure
#
# EVALUATION METRICS (all bootstrapped, 1000 resamples, percentile CIs):
#   Sensitivity (TPR), Specificity (TNR), PPV (precision), NPV,
#   False Positive Rate (FPR = 1 - Specificity), F1 score, Balanced Accuracy.
# Reporting the full confusion-matrix set - not just sensitivity/PPV - is
# what lets a reviewer see exactly how much specificity/precision is traded
# away when the trigger is loosened from the mean to a prediction interval.
#
# Inputs:  outputs/prepared_data.rds, outputs/lag_tables_sensitivity.rds
#          (both already produced by 01a and 03a - no need to re-run 05 first)
# Outputs: outputs/tables/spi_supp_table5_alert_sensitivity_PI_long.csv
#          outputs/tables/spi_supp_table5_alert_sensitivity_PI_report.csv
#          outputs/figures/supp_alert_sensitivity_PI.png  (300 dpi, manuscript-ready)
#          outputs/figures/supp_alert_sensitivity_PI.pdf  (vector, for print submission)
# ============================================================================

source("00_setup.R")

prepped     <- readRDS(file.path(OUT_DIR, "prepared_data.rds"))
weekly      <- prepped$weekly

lag_results <- readRDS(file.path(OUT_DIR, "lag_tables_sensitivity.rds"))
lag_table   <- lag_results$interp$lag_table

SPI_SCALE       <- 4
clim_col        <- paste0("SPI_", SPI_SCALE)
REFIT_EVERY     <- 12
OUTBREAK_PCTILE <- 0.95
MAX_LEAD        <- 4
PI_LEVELS       <- c(80, 95)   # matches forecast::forecast(..., level = PI_LEVELS)
N_BOOT          <- 1000

# Full label (table/report) vs. short label (figure axis) for each rule
RULE_LABELS <- c(
  mean    = "Forecast mean (primary analysis)",
  upper80 = "80% upper prediction limit",
  upper95 = "95% upper prediction limit"
)
RULE_LABELS_SHORT <- c(
  mean    = "Mean",
  upper80 = "80% PI",
  upper95 = "95% PI"
)

# ---------------------------------------------------------------------------
# Rolling-origin evaluation that keeps the forecast mean AND both upper
# prediction limits at every step (the only change vs. 05's
# rolling_origin_eval(): three prediction series come out instead of one).
# ---------------------------------------------------------------------------
rolling_origin_eval_PI <- function(df_state, cons_lag, eval_start_idx, season = 52,
                                   outbreak_pctile = OUTBREAK_PCTILE,
                                   pi_levels = PI_LEVELS) {
  
  x_lag <- dplyr::lag(df_state[[clim_col]], cons_lag)
  train_hist_threshold <- quantile(df_state$Cases[1:eval_start_idx], outbreak_pctile, na.rm = TRUE)
  
  eval_idx <- seq(eval_start_idx, nrow(df_state) - 1)
  n_eval   <- length(eval_idx)
  
  preds_mean    <- numeric(n_eval)
  preds_upper80 <- numeric(n_eval)
  preds_upper95 <- numeric(n_eval)
  refit_fit <- NULL
  
  for (i in seq_len(n_eval)) {
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
        forecast::forecast(refit_fit, h = 1, level = pi_levels,
                           xreg = matrix(x_lag[t + 1], ncol = 1,
                                         dimnames = list(NULL, "SPI_lag"))),
        error = function(e) NULL
      )
      if (!is.null(fc)) {
        preds_mean[i]    <- as.numeric(fc$mean)
        preds_upper80[i] <- as.numeric(fc$upper[, 1])  # first level in pi_levels (80%)
        preds_upper95[i] <- as.numeric(fc$upper[, 2])  # second level in pi_levels (95%)
      } else {
        preds_mean[i] <- NA; preds_upper80[i] <- NA; preds_upper95[i] <- NA
      }
    } else {
      preds_mean[i] <- NA; preds_upper80[i] <- NA; preds_upper95[i] <- NA
    }
  }
  
  actual <- df_state$Cases[eval_idx + 1]
  
  list(
    actual = actual,
    preds_mean = preds_mean, preds_upper80 = preds_upper80, preds_upper95 = preds_upper95,
    eval_idx = eval_idx, hist_threshold = train_hist_threshold
  )
}

# ---------------------------------------------------------------------------
# Detection rate / lead time for ANY prediction series against the same
# outbreak-week definition and lookback window used in 05 (identical logic,
# generalized to take the prediction vector as an argument).
# ---------------------------------------------------------------------------
detection_summary <- function(actual, preds, threshold, max_lead = MAX_LEAD) {
  outbreak_weeks <- which(actual > threshold)
  detected <- sapply(outbreak_weeks, function(w) {
    lookback <- max(1, w - max_lead):max(1, w - 1)
    any(preds[lookback] > threshold, na.rm = TRUE)
  })
  lead_times <- sapply(outbreak_weeks[detected], function(w) {
    lookback <- max(1, w - max_lead):max(1, w - 1)
    hit <- lookback[which(preds[lookback] > threshold)]
    if (length(hit) == 0) return(NA_real_)
    w - max(hit)
  })
  list(
    n_outbreaks = length(outbreak_weeks),
    n_detected  = sum(detected),
    detection_rate = if (length(outbreak_weeks) > 0) mean(detected) else NA_real_,
    mean_lead = if (length(lead_times) > 0) mean(lead_times, na.rm = TRUE) else NA_real_
  )
}

# ---------------------------------------------------------------------------
# Bootstrap the full confusion-matrix metric set for a given prediction
# series & threshold: Sensitivity, Specificity, PPV, NPV, FPR, F1, Balanced
# Accuracy - each with a percentile bootstrap CI (1000 resamples by default).
# Returns one row PER METRIC (long format) so it can be bound straight into
# the state x rule table and reused directly for faceted plotting.
# ---------------------------------------------------------------------------
METRIC_NAMES <- c("sensitivity", "specificity", "ppv", "npv", "fpr", "f1", "balanced_accuracy")

alert_bootstrap_full <- function(actual, predicted, red_thresh, n_boot = N_BOOT) {
  ok <- !is.na(predicted) & !is.na(actual)
  actual <- actual[ok]; predicted <- predicted[ok]
  n <- length(actual)
  
  if (n == 0) {
    return(tibble(metric = METRIC_NAMES, estimate = NA_real_, lo = NA_real_, hi = NA_real_))
  }
  
  boot_stat <- function(data, idx) {
    a <- data$actual[idx]; p <- data$predicted[idx]
    true_red <- a > red_thresh
    pred_red <- p > red_thresh
    
    TP <- sum(pred_red & true_red)
    FP <- sum(pred_red & !true_red)
    FN <- sum(!pred_red & true_red)
    TN <- sum(!pred_red & !true_red)
    
    sens <- if ((TP + FN) > 0) TP / (TP + FN) else NA_real_
    spec <- if ((TN + FP) > 0) TN / (TN + FP) else NA_real_
    ppv  <- if ((TP + FP) > 0) TP / (TP + FP) else NA_real_
    npv  <- if ((TN + FN) > 0) TN / (TN + FN) else NA_real_
    fpr  <- if ((FP + TN) > 0) FP / (FP + TN) else NA_real_
    f1   <- if (!is.na(ppv) && !is.na(sens) && (ppv + sens) > 0) 2 * ppv * sens / (ppv + sens) else NA_real_
    bacc <- if (!is.na(sens) && !is.na(spec)) (sens + spec) / 2 else NA_real_
    
    c(sensitivity = sens, specificity = spec, ppv = ppv, npv = npv,
      fpr = fpr, f1 = f1, balanced_accuracy = bacc)
  }
  
  bd <- data.frame(actual = actual, predicted = predicted)
  b <- tryCatch(boot::boot(bd, boot_stat, R = n_boot), error = function(e) NULL)
  
  if (is.null(b)) {
    return(tibble(metric = METRIC_NAMES, estimate = NA_real_, lo = NA_real_, hi = NA_real_))
  }
  
  estimate <- as.numeric(b$t0)
  
  map_dfr(seq_along(METRIC_NAMES), function(i) {
    ci <- tryCatch(boot::boot.ci(b, type = "perc", index = i)$percent[4:5],
                   error = function(e) c(NA_real_, NA_real_))
    tibble(metric = METRIC_NAMES[i], estimate = estimate[i], lo = ci[1], hi = ci[2])
  })
}

# ---------------------------------------------------------------------------
# Run the PI-based rolling evaluation once per state, then score all three
# alert rules (Mean / 80% PI / 95% PI) against the identical outbreak
# threshold and lookback window, with the full bootstrapped metric set.
# ---------------------------------------------------------------------------
pi_eval_results <- list()

supp_table5_long <- map_dfr(FOCAL_STATES, function(st) {
  df_state <- weekly %>% filter(State == st) %>% arrange(week_idx)
  
  cons_lag <- lag_table %>% filter(State == st) %>% pull(Consensus_lag_wk)
  if (length(cons_lag) != 1) {
    stop(sprintf("Expected exactly one Consensus_lag_wk for state '%s', got %d.", st, length(cons_lag)))
  }
  
  eval_start_idx <- df_state %>% filter(Year == 2024, Epi_Week == 1) %>% pull(week_idx)
  if (length(eval_start_idx) == 0) {
    eval_start_idx <- nrow(df_state) - 96
  } else {
    eval_start_idx <- eval_start_idx[1]
  }
  
  res <- rolling_origin_eval_PI(df_state, cons_lag, eval_start_idx)
  pi_eval_results[[st]] <<- res
  
  pred_series <- list(
    mean    = res$preds_mean,
    upper80 = res$preds_upper80,
    upper95 = res$preds_upper95
  )
  
  map_dfr(names(pred_series), function(rule) {
    preds <- pred_series[[rule]]
    det <- detection_summary(res$actual, preds, res$hist_threshold)
    mt  <- alert_bootstrap_full(res$actual, preds, res$hist_threshold)
    
    mt %>% mutate(
      State  = st,
      Rule   = RULE_LABELS[[rule]],
      RuleShort = RULE_LABELS_SHORT[[rule]],
      Outbreaks_detected = sprintf("%d/%d", det$n_detected, det$n_outbreaks),
      Detection_rate = ifelse(is.na(det$detection_rate), NA_character_,
                              sprintf("%.0f%%", 100 * det$detection_rate)),
      Lead_weeks = round(det$mean_lead, 1),
      .before = 1
    )
  })
})

# Long format: one row per State x Rule x Metric - full reproducibility record
write_csv(supp_table5_long, file.path(TAB_DIR, "spi_supp_table5_alert_sensitivity_PI_long.csv"))

# Human-readable "report" table: one row per State x Rule, each metric shown
# as "estimate [lo, hi]" - ready to paste into a supplementary Word/PDF table.
fmt_est_ci <- function(est, lo, hi) {
  ifelse(is.na(est), "-", sprintf("%.2f [%.2f, %.2f]", est, lo, hi))
}

supp_table5_report <- supp_table5_long %>%
  mutate(cell = fmt_est_ci(estimate, lo, hi)) %>%
  select(State, Rule, Outbreaks_detected, Detection_rate, Lead_weeks, metric, cell) %>%
  pivot_wider(names_from = metric, values_from = cell) %>%
  rename(
    Sensitivity          = sensitivity,
    Specificity          = specificity,
    PPV                  = ppv,
    NPV                  = npv,
    `FPR`                = fpr,
    `F1 score`           = f1,
    `Balanced accuracy`  = balanced_accuracy
  ) %>%
  arrange(State, Rule)

write_csv(supp_table5_report, file.path(TAB_DIR, "spi_supp_table5_alert_sensitivity_PI_report.csv"))

print(supp_table5_report, n = Inf, width = Inf)

# ---------------------------------------------------------------------------
# Manuscript-ready figure: Sensitivity, Specificity, PPV, NPV by trigger
# rule, faceted by state. FPR/F1/Balanced Accuracy are reported in the table
# only (FPR = 1 - Specificity is redundant to plot; F1/Balanced Accuracy are
# composite indices better read as numbers than as bars with CIs).
#
# No title/subtitle baked into the image - captions belong in the manuscript
# text, not the artwork. Sized to a single manuscript page width, vector PDF
# included for print submission, PNG at 300 dpi for word processors/preprint
# servers.
# ---------------------------------------------------------------------------
plot_metrics <- c("sensitivity", "specificity", "ppv", "npv")
plot_metric_labels <- c(sensitivity = "Sensitivity", specificity = "Specificity",
                        ppv = "PPV", npv = "NPV")

plot_df <- supp_table5_long %>%
  filter(metric %in% plot_metrics) %>%
  mutate(
    Metric = factor(plot_metric_labels[metric], levels = unname(plot_metric_labels)),
    RuleShort = factor(RuleShort, levels = unname(RULE_LABELS_SHORT))
  )

supp_plot <- ggplot(plot_df, aes(x = RuleShort, y = estimate, fill = RuleShort)) +
  geom_col(width = 0.7, color = "grey20", linewidth = 0.2) +
  geom_errorbar(aes(ymin = pmax(lo, 0), ymax = pmin(hi, 1)), width = 0.15, linewidth = 0.35, na.rm = TRUE) +
  facet_grid(Metric ~ State) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1), expand = expansion(mult = c(0, 0.03))) +
  scale_fill_manual(values = c("Mean" = "#95A5A6", "80% PI" = "#F39C12", "95% PI" = "#C0392B"), guide = "none") +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 9, base_family = "") +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.background = element_rect(fill = "grey92", color = "grey40", linewidth = 0.3),
    strip.text = element_text(face = "bold", size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    axis.text.y = element_text(size = 7),
    panel.spacing = unit(0.35, "lines"),
    plot.margin = margin(4, 6, 4, 4)
  )

ggsave(file.path(FIG_DIR, "supp_alert_sensitivity_PI.png"), supp_plot,
       width = 7.0, height = 6.2, dpi = 300, bg = "white")
ggsave(file.path(FIG_DIR, "supp_alert_sensitivity_PI.pdf"), supp_plot,
       width = 7.0, height = 6.2, device = cairo_pdf, bg = "white")

cat("\nAlert-rule sensitivity analysis complete. Bauchi should remain the\n",
    "strongest-performing state under every rule if the primary conclusion is robust;\n",
    "sensitivity/NPV are expected to rise and specificity/PPV to fall moving\n",
    "Mean -> 80% PI -> 95% PI. Ebonyi returning near-zero under every rule indicates\n",
    "a genuine absence of climate-driven predictive signal there, not an overly\n",
    "strict threshold.\n\n", sep = "")

cat("Saved:\n",
    " -", file.path(TAB_DIR, "spi_supp_table5_alert_sensitivity_PI_long.csv"), "\n",
    " -", file.path(TAB_DIR, "spi_supp_table5_alert_sensitivity_PI_report.csv"), "\n",
    " -", file.path(FIG_DIR, "supp_alert_sensitivity_PI.png"), "\n",
    " -", file.path(FIG_DIR, "supp_alert_sensitivity_PI.pdf"), "\n")


################################################################################
# Define metric mappings (including NPV)
metrics_to_plot <- c(
  sensitivity = "Sensitivity",
  ppv         = "PPV",
  npv         = "NPV"
)

# Rule label mapping for x-axis
rule_label_map <- c(
  "Mean"   = "Forecast Mean",
  "80% PI" = "80% Upper PI",
  "95% PI" = "95% Upper PI"
)

# Iterate over each metric and output individual figures
for (m_key in names(metrics_to_plot)) {
  
  df_sub <- supp_table5_long %>%
    filter(metric == m_key) %>%
    mutate(
      # Map short rules to full clear titles
      RuleShort = factor(
        dplyr::recode(RuleShort, !!!rule_label_map), 
        levels = unname(rule_label_map)
      )
    )
  
  p_single <- ggplot(df_sub, aes(x = RuleShort, y = estimate, fill = RuleShort)) +
    geom_col(width = 0.65, color = "grey20", linewidth = 0.2) +
    geom_errorbar(
      aes(ymin = pmax(lo, 0), ymax = pmin(hi, 1)), 
      width = 0.15, linewidth = 0.35, na.rm = TRUE
    ) +
    facet_wrap(~ State, nrow = 1) +
    scale_y_continuous(
      limits = c(0, 1), 
      breaks = seq(0, 1, 0.25),
      labels = c("0", "0.25", "0.50", "0.75", "1.00"),
      expand = expansion(mult = c(0, 0.03))
    ) +
    scale_fill_manual(
      values = c(
        "Forecast Mean" = "#95A5A6", 
        "80% Upper PI"  = "#F39C12", 
        "95% Upper PI"  = "#C0392B"
      ), 
      guide = "none"
    ) +
    labs(x = NULL, y = metrics_to_plot[[m_key]]) +
    theme_bw(base_size = 9) +
    theme(
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_blank(),
      strip.background   = element_rect(fill = "grey92", color = "grey40", linewidth = 0.3),
      strip.text         = element_text(face = "bold", size = 8.5),
      axis.title.y       = element_text(size = 9, face = "bold"),
      axis.text.x        = element_text(angle = 45, hjust = 1, size = 7.5, color = "black"),
      axis.text.y        = element_text(size = 7.5, color = "black"),
      panel.spacing      = unit(0.4, "lines"),
      plot.margin        = margin(4, 6, 4, 4)
    )
  
  # Save PNG and PDF for each metric
  ggsave(file.path(FIG_DIR, paste0("supp_alert_", m_key, ".png")), p_single,
         width = 7.0, height = 3.2, dpi = 300, bg = "white")
}


library(ggplot2)
library(dplyr)

# Define metric mappings (including NPV)
metrics_to_plot <- c(
  sensitivity = "Sensitivity",
  ppv         = "PPV",
  npv         = "NPV"
)

# x-axis labels for the three alert-trigger rules being compared
rule_label_map <- c(
  "Mean"   = "Mean",
  "80% PI" = "80% PI",
  "95% PI" = "95% PI"
)

# Iterate over each metric and output individual figures
for (m_key in names(metrics_to_plot)) {
  
  df_sub <- supp_table5_long %>%
    filter(metric == m_key) %>%
    mutate(
      # Map short rules to consistent x-axis labels
      RuleShort = factor(
        dplyr::recode(RuleShort, !!!rule_label_map), 
        levels = unname(rule_label_map)
      )
    )
  
  p_single <- ggplot(df_sub, aes(x = RuleShort, y = estimate, fill = RuleShort)) +
    geom_col(width = 0.65, color = "grey20", linewidth = 0.2) +
    geom_errorbar(
      aes(ymin = pmax(lo, 0), ymax = pmin(hi, 1)), 
      width = 0.15, linewidth = 0.35, na.rm = TRUE
    ) +
    facet_wrap(~ State, nrow = 1) +
    scale_y_continuous(
      limits = c(0, 1), 
      breaks = seq(0, 1, 0.25),
      labels = c("0", "0.25", "0.50", "0.75", "1.00"),
      expand = expansion(mult = c(0, 0.03))
    ) +
    scale_fill_manual(
      values = c(
        "Mean"   = "#95A5A6", 
        "80% PI" = "#F39C12", 
        "95% PI" = "#C0392B"
      ), 
      guide = "none"
    ) +
    labs(x = NULL, y = metrics_to_plot[[m_key]]) +
    theme_bw(base_size = 9) +
    theme(
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_blank(),
      strip.background   = element_rect(fill = "grey92", color = "grey40", linewidth = 0.3),
      strip.text         = element_text(face = "bold", size = 8.5),
      axis.title.y       = element_text(size = 9, face = "bold"),
      axis.text.x        = element_text(angle = 45, hjust = 1, size = 7.5, color = "black"),
      axis.text.y        = element_text(size = 7.5, color = "black"),
      panel.spacing      = unit(0.4, "lines"),
      plot.margin        = margin(4, 6, 4, 4)
    )
  
  # Save PNG and PDF for each metric
  ggsave(file.path(FIG_DIR, paste0("supp_alert_", m_key, ".png")), p_single,
         width = 7.0, height = 3.2, dpi = 300, bg = "white")
}
