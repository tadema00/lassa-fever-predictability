# ============================================================================
# 03a_spi_wavelet_lag_analysis.R
# Wavelet coherence + AR(1) Monte Carlo surrogate significance test
# using Standardized Precipitation Index (SPI-4 / SPI-8) as the climatic driver,
# three-way lag triangulation -> consensus lag (Table 1),
# and peak-timing consistency (Table 2), run on BOTH the interpolated
# ("weekly") and original ("weekly_original") datasets.
#
# Loads `baseline_params` from prepared_data.rds (produced by 01a) rather
# than assuming it exists in the session, so this script can be run
# independently. Saves `lag_tables_sensitivity.rds` at the end, with
# `$interp$lag_table` / `$original$lag_table` - this is the exact object
# that 04_spi_sarimax_zinb_models.R and 05_spi_rolling_forecast_alerts.R
# both read back in.
# ============================================================================

source("000_setup.R")

prepped         <- readRDS(file.path(OUT_DIR, "prepared_data.rds"))
weekly          <- prepped$weekly
weekly_original <- prepped$weekly_original
baseline_params <- prepped$baseline_params

compute_spi_from_baseline <- function(df, baseline_df, scale_val) {
  acc_col <- paste0("Rain_acc_", scale_val)
  spi_col <- paste0("SPI_", scale_val)
  
  base_sub <- baseline_df %>% filter(Scale == scale_val)
  
  df_joined <- df %>%
    left_join(base_sub %>% select(State, Epi_Week, shape, rate, q, empirical_fn), by = c("State", "Epi_Week"))
  
  spi_vals <- numeric(nrow(df_joined))
  
  for (i in seq_len(nrow(df_joined))) {
    val <- df_joined[[acc_col]][i]
    q_val <- df_joined$q[i]
    shape <- df_joined$shape[i]
    rate <- df_joined$rate[i]
    
    emp_cell <- df_joined$empirical_fn[i]
    emp_fn <- if (is.list(emp_cell)) emp_cell[[1]] else NULL
    
    if (is.na(val) || is.na(q_val)) {
      spi_vals[i] <- NA_real_
      next
    }
    
    if (val == 0) {
      H_val <- q_val
    } else {
      if (is.na(shape) || is.na(rate) || shape <= 0 || rate <= 0 || !is.function(emp_fn)) {
        prob_nz <- if (is.function(emp_fn)) emp_fn(val) else 0.5
        H_val <- q_val + (1 - q_val) * prob_nz
      } else {
        prob_nz <- pgamma(val, shape = shape, rate = rate)
        H_val <- q_val + (1 - q_val) * prob_nz
      }
    }
    
    H_val <- pmin(pmax(H_val, 1e-10), 1.0 - 1e-10)
    s_val <- qnorm(H_val)
    spi_vals[i] <- if (is.finite(s_val)) s_val else 0.0
  }
  
  df_joined[[spi_col]] <- spi_vals
  df_joined %>% select(-shape, -rate, -q, -empirical_fn)
}

# Recompute SPI on weekly_original 
if (!"SPI_4" %in% names(weekly_original)) {
  weekly_original <- weekly_original %>%
    group_by(State) %>%
    arrange(Year, week_idx) %>%
    mutate(
      Rain_acc_1 = Rainfall_mm,
      Rain_acc_4 = rollapply(Rainfall_mm, width = 4, FUN = sum, align = "right", fill = NA, na.rm = FALSE),
      Rain_acc_8 = rollapply(Rainfall_mm, width = 8, FUN = sum, align = "right", fill = NA, na.rm = FALSE)
    ) %>%
    ungroup()
  
  weekly_original <- compute_spi_from_baseline(weekly_original, baseline_params, 1)
  weekly_original <- compute_spi_from_baseline(weekly_original, baseline_params, 4)
  weekly_original <- compute_spi_from_baseline(weekly_original, baseline_params, 8)
}

ANNUAL_BAND    <- c(48, 56)   # weeks
N_SURROGATES   <- 1000
N_MC_PVALUE    <- 300
RUN_MC_PVALUE  <- FALSE
MAX_CCF_LAG <- 30
MAX_AIC_LAG <- 30
MIN_LAG_OBS  <- 30
SPI_SCALE      <- 4           # SPI-4 as the primary meteorological drought index for hydrological/epidemiological coupling

clim_col       <- paste0("SPI_", SPI_SCALE)

# ---------------------------------------------------------------------------
# Helper: fit an AR(1) surrogate
# ---------------------------------------------------------------------------
ar1_surrogate <- function(x) {
  n <- length(x)
  mu <- mean(x, na.rm = TRUE)
  xc <- x - mu
  
  alpha <- cor(xc[-n], xc[-1], use = "complete.obs")
  if (is.na(alpha)) alpha <- 0.0
  alpha <- max(min(alpha, 0.95), -0.95)
  
  sigma_z <- sd(xc, na.rm = TRUE) * sqrt(1 - alpha^2)
  
  z <- numeric(n)
  z[1] <- rnorm(1, sd = sd(xc, na.rm = TRUE))
  eps <- rnorm(n - 1, sd = sigma_z)
  for (t in 2:n) z[t] <- alpha * z[t - 1] + eps[t - 1]
  z + mu
}

# ---------------------------------------------------------------------------
# Wavelet coherence for one state
# ---------------------------------------------------------------------------
run_wtc <- function(week_idx, clim_var, cases, nrands) {
  x <- cbind(week_idx, clim_var)
  y_log <- log1p(cases)
  y <- cbind(
    week_idx,
    scale(y_log)[, 1]
  )
  
  tryCatch({
    biwavelet::wtc(x, y, nrands = nrands, quiet = TRUE)
  }, error = function(e) {
    tryCatch({
      biwavelet::wtc(x, y, nrands = 0, quiet = TRUE)
    }, error = function(e2) {
      return(NULL)
    })
  })
}

summarize_band <- function(wtc_res) {
  period   <- wtc_res$period
  band_idx <- which(period >= ANNUAL_BAND[1] & period <= ANNUAL_BAND[2])
  valid_mask <- outer(period, wtc_res$coi, FUN = function(p, c) p <= c)
  
  band_rsq   <- wtc_res$rsq[band_idx, , drop = FALSE]
  band_phase <- wtc_res$phase[band_idx, , drop = FALSE]
  valid      <- valid_mask[band_idx, , drop = FALSE]
  
  mean_rsq   <- mean(band_rsq[valid], na.rm = TRUE)
  mean_phase <- mean(band_phase[valid], na.rm = TRUE)
  mean_period <- mean(period[band_idx])
  lag_weeks <- round((mean_phase / (2 * pi)) * mean_period)
  
  list(mean_rsq = mean_rsq, lag_weeks = lag_weeks, valid = valid, band_idx = band_idx)
}

prop_sig_annual_band <- function(wtc_res, band_summary) {
  signif_band <- wtc_res$signif[band_summary$band_idx, , drop = FALSE]
  band_rsq    <- wtc_res$rsq[band_summary$band_idx, , drop = FALSE]
  mean((band_rsq[band_summary$valid] / signif_band[band_summary$valid]) >= 1, na.rm = TRUE)
}

ccf_lag <- function(wrai, cases, max_lag = MAX_CCF_LAG) {
  cc <- ccf(wrai, cases, lag.max = max_lag, plot = FALSE)
  lags <- cc$lag[, 1, 1]
  vals <- cc$acf[, 1, 1]
  nonneg <- lags >= 0
  lags[nonneg][which.max(vals[nonneg])]
}

aic_lag <- function(wrai, cases, max_lag = MAX_AIC_LAG) {
  aics <- sapply(0:max_lag, function(l) {
    x_lag <- dplyr::lag(wrai, l)
    ok <- !is.na(x_lag)
    if (sum(ok) < MIN_LAG_OBS) return(Inf)
    fit <- glm(cases[ok] ~ x_lag[ok], family = poisson())
    AIC(fit)
  })
  list(best_lag = (0:max_lag)[which.min(aics)], best_aic = min(aics))
}
# ---------------------------------------------------------------------------
# Full analysis function using SPI
# ---------------------------------------------------------------------------
run_analysis_spi <- function(weekly_df, label) {
  
  # CRITICAL FIX: Eliminate Temporal Leakage
  # Estimate lags based ONLY on data available before the forecast period starts.
  weekly_df_train <- weekly_df %>% filter(Year < 2024)
  
  lag_table <- map_dfr(FOCAL_STATES, function(st) {
    # Use the TRAINING data for lag estimation
    df_state_raw <- weekly_df_train %>% filter(State == st) %>% arrange(week_idx)
    df_state     <- df_state_raw %>% tidyr::drop_na(all_of(clim_col), Cases)
    
    # Run Wavelet Coherence on training period
    wtc_res <- run_wtc(df_state$week_idx, df_state[[clim_col]], df_state$Cases, nrands = N_SURROGATES)
    
    if (is.null(wtc_res)) {
      return(tibble(Dataset = label, State = st, Consensus_lag_wk = NA_real_))
    }
    
    band_summary <- summarize_band(wtc_res)
    prop_sig     <- prop_sig_annual_band(wtc_res, band_summary)
    
    # Run CCF and AIC on training period
    ccf_l <- ccf_lag(df_state[[clim_col]], df_state$Cases)
    aic_l <- aic_lag(df_state[[clim_col]], df_state$Cases)
    
    consensus_lag <- median(c(band_summary$lag_weeks, ccf_l, aic_l$best_lag), na.rm = TRUE)
    
    tibble(
      Dataset              = label,
      State                = st,
      Wavelet_lag_wk       = band_summary$lag_weeks,
      CCF_lag_wk           = ccf_l,
      AIC_lag_wk           = aic_l$best_lag,
      Consensus_lag_wk     = consensus_lag,
      Annual_Rsq           = round(band_summary$mean_rsq, 3), # Note: Reviewer 2 wants this renamed in text
      Prop_Sig_Annual_Band = round(prop_sig, 3)
    )
  })
  return(lag_table)
}

# ---------------------------------------------------------------------------
# Run on both datasets and compare
# ---------------------------------------------------------------------------
set.seed(SEED_RESAMPLING)
interp_results   <- run_analysis_spi(weekly,          "Interpolated")
set.seed(SEED_RESAMPLING)
original_results <- run_analysis_spi(weekly_original, "Original")

spi_comparison <- bind_rows(interp_results, original_results) |>
  arrange(State, Dataset)

print(spi_comparison)
write_csv(spi_comparison, file.path(TAB_DIR, "spi_table1_lag_triangulation_sensitivity.csv"))

# Save the per-dataset lag tables in the exact shape 04_spi_sarimax_zinb_models.R
# and 05_spi_rolling_forecast_alerts.R expect: lag_results$interp$lag_table /
# lag_results$original$lag_table. This is what was missing before.
saveRDS(
  list(
    interp   = list(lag_table = interp_results),
    original = list(lag_table = original_results)
  ),
  file.path(OUT_DIR, "lag_tables_sensitivity.rds")
)

cat("\nSPI-based wavelet lag analysis and triangulation sensitivity complete.\n")


# ============================================================================
# 03b_spi_wavelet_figures.R (Refined Phase Wrapping & Plotting)
# ============================================================================

source("00_setup.R")

prepped <- readRDS(file.path(OUT_DIR, "prepared_data.rds"))
weekly  <- prepped$weekly

ANNUAL_BAND  <- c(48, 56)
PERIOD_TICKS <- c(2, 4, 8, 16, 32, 52, 104)
SPI_SCALE    <- 4
clim_col     <- paste0("SPI_", SPI_SCALE)

# ----------------------------------------------------------------------------
# 1. Compute Wavelet Coherence (WTC)
# ----------------------------------------------------------------------------
wtc_for_state <- function(state) {
  df_state <- weekly %>%
    filter(State == state) %>%
    arrange(week_idx) %>%
    tidyr::drop_na(all_of(clim_col), Cases)
  
  x <- cbind(df_state$week_idx, df_state[[clim_col]])
  y <- cbind(df_state$week_idx, scale(log1p(df_state$Cases))[, 1])
  
  wtc_obj <- biwavelet::wtc(x, y, nrands = 300, quiet = TRUE)
  wtc_obj$state_name <- state
  return(wtc_obj)
}

# ----------------------------------------------------------------------------
# 2. Extract Phase Lead with Unwrapped Phase Boundary [-90°, +270°]
# ----------------------------------------------------------------------------
extract_annual_phase <- function(wtc_obj, band = ANNUAL_BAND) {
  band_idx <- which(wtc_obj$period >= band[1] & wtc_obj$period <= band[2])
  if (length(band_idx) == 0) band_idx <- which.min(abs(wtc_obj$period - 52))
  
  phase_sub <- wtc_obj$phase[band_idx, , drop = FALSE]
  rsq_sub   <- wtc_obj$rsq[band_idx, , drop = FALSE]
  
  sin_mean <- colMeans(sin(phase_sub), na.rm = TRUE)
  cos_mean <- colMeans(cos(phase_sub), na.rm = TRUE)
  mean_phase <- atan2(sin_mean, cos_mean) # radians [-pi, pi]
  mean_rsq   <- colMeans(rsq_sub, na.rm = TRUE)
  
  # Shift domain from [-pi, pi] to [-pi/2, 3pi/2] (i.e. [-90°, 270°])
  # This prevents -180° / +180° anti-phase wrapping line snaps
  phase_shifted <- ifelse(mean_phase < -pi/2, mean_phase + 2*pi, mean_phase)
  
  # Lead time in weeks (0 to 52 week cycle)
  lead_time_weeks <- (phase_shifted / (2 * pi)) * 52
  
  data.frame(
    week_idx        = wtc_obj$t,
    State           = wtc_obj$state_name,
    phase_rad       = mean_phase,
    phase_deg       = mean_phase * (180 / pi),
    lead_time_weeks = lead_time_weeks,
    coherence_rsq   = mean_rsq
  )
}

add_period_axis <- function() {
  axis(2, at = log2(PERIOD_TICKS), labels = PERIOD_TICKS, las = 1)
}

# ----------------------------------------------------------------------------
# 3. Wavelet Coherence Panel Plotter
# ----------------------------------------------------------------------------
plot_cluster_coherence_panel <- function(states, panel_letters, filename, band = ANNUAL_BAND) {
  tiff(file.path(FIG_DIR, filename), width = 12, height = 9, units = "in",
       res = 600, compression = "lzw")
  on.exit(dev.off(), add = TRUE)
  
  layout(matrix(c(1, 2, 3, 3), nrow = 2, byrow = TRUE))
  par(mar = c(4, 4, 2.5, 1))
  
  wtc_list <- lapply(states, wtc_for_state)
  
  for (i in seq_along(states)) {
    wtc_obj <- wtc_list[[i]]
    plot(wtc_obj,
         plot.phase = TRUE,
         plot.coi   = TRUE, lwd.coi = 2, col.coi = "black",
         lwd.sig    = 2,
         plot.cb    = FALSE,
         yaxt       = "n",
         main       = sprintf("(%s) %s Wavelet Coherence (SPI-4 vs. Cases)", panel_letters[i], states[i]),
         xlab       = "Week Index", ylab = "Period (weeks)")
    
    add_period_axis()
    abline(h = log2(band), lty = 2, col = "red", lwd = 2)
  }
  return(wtc_list)
}

# ----------------------------------------------------------------------------
# 4. Phase & Lead-Time Panel Plotter 
# ----------------------------------------------------------------------------
plot_cluster_phase_panel <- function(wtc_list, panel_letters, filename, band = ANNUAL_BAND) {
  phase_dfs <- lapply(wtc_list, extract_annual_phase, band = band)
  combined_df <- bind_rows(phase_dfs)
  
  # Smooth out remaining boundary jumps by inserting NA across wrap transitions
  combined_df <- combined_df %>%
    group_by(State) %>%
    mutate(
      jump = c(0, abs(diff(lead_time_weeks))) > 15,
      lead_time_clean = if_else(jump, NA_real_, lead_time_weeks),
      sig_label = if_else(coherence_rsq >= 0.5, "High Coherence (R² ≥ 0.5)", "Low Coherence (R² < 0.5)")
    ) %>%
    ungroup()
  
  # Figure A: Continuous Phase Lead Time
  p1 <- ggplot(combined_df, aes(x = week_idx, color = State, group = State)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.6) +
    geom_line(aes(y = lead_time_clean), alpha = 0.85, linewidth = 0.9, na.rm = TRUE) +
    geom_point(aes(y = lead_time_weeks, shape = sig_label, alpha = coherence_rsq), size = 1.5) +
    scale_shape_manual(values = c("High Coherence (R² ≥ 0.5)" = 16, "Low Coherence (R² < 0.5)" = 1)) +
    scale_alpha_continuous(range = c(0.2, 0.95), guide = "none") +
    scale_y_continuous(breaks = seq(-12, 36, 6)) +
    coord_cartesian(ylim = c(-14, 40)) +  # Zooms without dropping data or throwing warnings
    labs(
      title = "",
      subtitle = "",
      x = "Week Index",
      y = "Phase Lead (+) / Lag (-) [Weeks]",
      shape = "Coherence Level",
      color = "State"
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )
  
  # Figure B: Density of Phase Angles
  p2 <- ggplot(combined_df %>% filter(coherence_rsq >= 0.40), aes(x = phase_deg, fill = State)) +
    geom_density(alpha = 0.45) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.6) +
    scale_x_continuous(limits = c(-180, 180), breaks = seq(-180, 180, 45)) +
    labs(
      title = "",
      subtitle = "",
      x = "Phase Angle (Degrees)",
      y = "Density",
      fill = "State"
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(face = "bold")
    )
  
  combined_plot <- gridExtra::grid.arrange(p1, p2, ncol = 1, heights = c(1.2, 1))
  
  ggsave(
    filename = file.path(FIG_DIR, filename),
    plot     = combined_plot,
    width    = 11,
    height   = 8,
    dpi      = 600,
    compression = "lzw"
  )
}

# ============================================================================
# EXECUTION
# ============================================================================

cat("Generating Southern Forest Cluster Wavelet Coherence & Phase...\n")
wtc_south <- plot_cluster_coherence_panel(
  c("Ondo", "Ebonyi", "Edo"), 
  c("A", "B", "C"),
  "fig4_spi_wavelet_coherence_southern.tiff"
)
plot_cluster_phase_panel(wtc_south, c("A", "B"), "fig4_spi_wavelet_phase_southern.tiff")

cat("Generating Northern Savannah Cluster Wavelet Coherence & Phase...\n")
wtc_north <- plot_cluster_coherence_panel(
  c("Bauchi", "Taraba", "Plateau"), 
  c("A", "B", "C"),
  "fig5_spi_wavelet_coherence_northern.tiff"
)
plot_cluster_phase_panel(wtc_north, c("A", "B"), "fig5_spi_wavelet_phase_northern.tiff")

cat("SPI wavelet coherence and phase figures regenerated successfully.\n")

# ============================================================================
# 03c_peak_month_rayleigh_test.R
# Table 2: Peak-timing consistency (Rayleigh test)
# Computes the peak-timing consistency table referenced alongside the
# consensus-lag results in 03a_spi_wavelet_lag_analysis.R.
#
# WHAT "PEAK-TIMING CONSISTENCY" MEANS HERE: for each state, does the
# within-year peak week of Lassa fever incidence fall at roughly the same
# time of year across the 2018-2025 study period, or does it drift/vary
# unpredictably? Week-of-year is a circular (not linear) quantity - week 52
# and week 1 are adjacent, not 51 weeks apart - so a circular (directional)
# statistic is the correct tool, not an ordinary mean/SD or t-test. The
# Rayleigh test (Fisher 1993; Mardia & Jupp 2000) tests the null hypothesis
# that peak weeks are uniformly scattered around the year against the
# alternative that they cluster around one preferred direction (i.e. a
# genuine peak "season"); the mean resultant length (r, 0-1) measures how
# tightly they cluster regardless of significance.
#
# TWO TABLES ARE PRODUCED:
#   Table 2  - peak week of Lassa CASE incidence per state, across years
#              (the primary, manuscript-facing table).
#   Table 2b - peak (driest) week of SPI-4 per state, across years - a
#              companion analysis showing whether the climatic driver itself
#              has a stable annual minimum. This supports the "annual band"
#              wavelet-coherence interpretation in Figs 4-5, but is an
#              addition beyond the original Table 2 spec - drop it if you
#              only want the case-incidence table.
#
# Inputs:  outputs/prepared_data.rds
# Outputs: outputs/tables/spi_table2_peak_timing_rayleigh_cases.csv
#          outputs/tables/spi_table2b_peak_timing_rayleigh_spi4.csv
#          outputs/figures/supp_peak_week_rose_diagrams.png
# ============================================================================

source("000_setup.R")
library(circular)   # rayleigh.test()

prepped <- readRDS(file.path(OUT_DIR, "prepared_data.rds"))
weekly  <- prepped$weekly

SPI_SCALE <- 4
clim_col  <- paste0("SPI_", SPI_SCALE)

# A year needs most of its 52 weeks actually observed to trust its "peak
# week" - otherwise a still-incomplete year (e.g. the current year, or one
# with a long reporting gap) could show a false early "peak" simply because
# the true peak hasn't happened yet within the observed window.
MIN_WEEKS_FOR_YEAR <- 40

# ---------------------------------------------------------------------------
# Core helper: Rayleigh test + mean direction for a set of week-of-year
# values (1-52), returned as a one-row tibble. `n_period = 52` maps week 52
# back to adjacent to week 1 on the circle.
# ---------------------------------------------------------------------------
peak_week_rayleigh <- function(peak_weeks, n_period = 52) {
  peak_weeks <- peak_weeks[!is.na(peak_weeks)]
  n <- length(peak_weeks)
  
  if (n < 2) {
    return(tibble(
      n_years = n, mean_peak_week = NA_real_, mean_resultant_length = NA_real_,
      rayleigh_stat = NA_real_, rayleigh_p = NA_real_
    ))
  }
  
  angles_deg <- (peak_weeks - 1) / n_period * 360
  circ <- circular::circular(angles_deg, units = "degrees", template = "none",
                             modulo = "2pi", zero = 0, rotation = "counter")
  
  rt <- circular::rayleigh.test(circ)
  r  <- circular::rho.circular(circ)
  mean_dir_deg   <- as.numeric(circular::mean.circular(circ)) %% 360
  mean_peak_week <- (mean_dir_deg / 360) * n_period + 1
  
  tibble(
    n_years = n,
    mean_peak_week = round(mean_peak_week, 1),
    mean_resultant_length = round(as.numeric(r), 4),
    rayleigh_stat = round(as.numeric(rt$statistic), 4),
    rayleigh_p = signif(as.numeric(rt$p.value), 4)
  )
}

# ---------------------------------------------------------------------------
# Table 2: peak week of Lassa CASE incidence, per state per year, then
# Rayleigh-tested across years within each state.
# ---------------------------------------------------------------------------
year_summary_cases <- weekly %>%
  group_by(State, Year) %>%
  summarise(n_weeks = n(), total_cases = sum(Cases, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    excluded_incomplete_year = n_weeks < MIN_WEEKS_FOR_YEAR,
    excluded_zero_cases      = !excluded_incomplete_year & total_cases == 0
  )

case_peaks <- weekly %>%
  semi_join(
    year_summary_cases %>% filter(!excluded_incomplete_year, !excluded_zero_cases),
    by = c("State", "Year")
  ) %>%
  group_by(State, Year) %>%
  # Ties (multiple weeks sharing the max case count) are resolved by
  # which.max()'s "first occurrence" rule - a simplification worth noting if
  # a state has many tied/low-count years, but immaterial for states with a
  # clear seasonal spike.
  summarise(peak_epi_week = Epi_Week[which.max(Cases)], .groups = "drop")

table2_cases <- case_peaks %>%
  group_by(State) %>%
  group_modify(~ peak_week_rayleigh(.x$peak_epi_week)) %>%
  ungroup()

exclusion_counts <- year_summary_cases %>%
  group_by(State) %>%
  summarise(
    n_years_total              = n(),
    n_years_excluded_incomplete = sum(excluded_incomplete_year),
    n_years_excluded_zero_cases = sum(excluded_zero_cases),
    .groups = "drop"
  )

table2_cases <- table2_cases %>%
  left_join(exclusion_counts, by = "State") %>%
  select(State, n_years_total, n_years_excluded_incomplete, n_years_excluded_zero_cases,
         n_years, mean_peak_week, mean_resultant_length, rayleigh_stat, rayleigh_p) %>%
  arrange(State)

print(table2_cases, n = Inf)
write_csv(table2_cases, file.path(TAB_DIR, "spi_table2_peak_timing_rayleigh_cases.csv"))

# ---------------------------------------------------------------------------
# Table 2b (companion): peak DRIEST week of SPI-4, per state per year - is
# the climatic driver's own annual minimum stable in timing? Uses the same
# MIN_WEEKS_FOR_YEAR completeness rule; years where SPI_4 is entirely NA
# (e.g. the first weeks of the observed record for that state, before a full
# 4-week accumulation window exists) are dropped automatically via drop_na().
# ---------------------------------------------------------------------------
spi_peaks <- weekly %>%
  filter(!is.na(.data[[clim_col]])) %>%
  semi_join(
    year_summary_cases %>% filter(!excluded_incomplete_year),
    by = c("State", "Year")
  ) %>%
  group_by(State, Year) %>%
  summarise(trough_epi_week = Epi_Week[which.min(.data[[clim_col]])], .groups = "drop")

table2b_spi <- spi_peaks %>%
  group_by(State) %>%
  group_modify(~ peak_week_rayleigh(.x$trough_epi_week)) %>%
  ungroup() %>%
  arrange(State)

print(table2b_spi, n = Inf)
write_csv(table2b_spi, file.path(TAB_DIR, "spi_table2b_peak_timing_rayleigh_spi4.csv"))

# ---------------------------------------------------------------------------
# Table 2 (manuscript display format): calendar MONTH of peak weekly Lassa
# case counts, one row per state, one column per year - the human-readable
# format for the main text, alongside the statistical Rayleigh-test version
# above. Uses every year with at least one reported case (no completeness
# filter), since this is a descriptive display table rather than the
# inferential Rayleigh test, which needs the stricter exclusions above.
#
# CAVEAT (documented, not hidden): epi-week -> calendar month is approximated
# by treating epi-week 1 as starting Jan 1 of that year and adding 7 days per
# subsequent week. True epidemiological week standards (e.g. MMWR/ISO 8601
# weeks) anchor slightly differently at year boundaries, so a peak reported
# here as, say, early January vs. late December could shift by a few days
# under a different week-numbering convention. This does not affect the
# Rayleigh test above (which works directly in week-space), only the month
# label shown in this display table.
# ---------------------------------------------------------------------------
epiweek_to_month_abbr <- function(year, epi_week) {
  week_start_date <- as.Date(sprintf("%d-01-01", year)) + (epi_week - 1) * 7
  format(week_start_date, "%b")
}

case_peaks_all_years <- weekly %>%
  group_by(State, Year) %>%
  summarise(
    total_cases   = sum(Cases, na.rm = TRUE),
    peak_epi_week = Epi_Week[which.max(Cases)],
    .groups = "drop"
  ) %>%
  mutate(
    peak_month = if_else(total_cases > 0,
                         epiweek_to_month_abbr(Year, peak_epi_week),
                         NA_character_)
  )

year_cols_sorted <- sort(unique(case_peaks_all_years$Year))

table2_calendar <- case_peaks_all_years %>%
  select(State, Year, peak_month) %>%
  pivot_wider(names_from = Year, values_from = peak_month) %>%
  arrange(State) %>%
  select(State, all_of(as.character(year_cols_sorted)))

print(table2_calendar, n = Inf, width = Inf)
write_csv(table2_calendar, file.path(TAB_DIR, "spi_table2_calendar_month_by_state_year.csv"))

# ---------------------------------------------------------------------------
# Supplementary figure: rose diagram per state showing each year's case-peak
# week as a point around the annual circle, with the mean direction arrow.
# ---------------------------------------------------------------------------
png(file.path(FIG_DIR, "supp_peak_week_rose_diagrams.png"),
    width = 10, height = 7, units = "in", res = 300)
par(mfrow = c(2, 3), mar = c(1, 1, 2.5, 1))

for (st in FOCAL_STATES) {
  yrs <- case_peaks %>% filter(State == st) %>% pull(peak_epi_week)
  if (length(yrs) < 2) next
  
  angles_deg <- (yrs - 1) / 52 * 360
  circ <- circular::circular(angles_deg, units = "degrees", template = "none",
                             modulo = "2pi", zero = 0, rotation = "counter")
  
  circular::plot.circular(circ, stack = TRUE, bins = 52, shrink = 1.3,
                          main = st, sep = 0.05, col = "#2874A6", cex = 1.1)
  circular::arrows.circular(circular::mean.circular(circ),
                            y = circular::rho.circular(circ),
                            col = "#C0392B", lwd = 2, length = 0.12)
}
dev.off()

cat("\nTable 2 (peak-timing consistency, Rayleigh test) complete.\n")
cat("Saved:\n",
    " -", file.path(TAB_DIR, "spi_table2_peak_timing_rayleigh_cases.csv"), "\n",
    " -", file.path(TAB_DIR, "spi_table2_calendar_month_by_state_year.csv"), "\n",
    " -", file.path(TAB_DIR, "spi_table2b_peak_timing_rayleigh_spi4.csv"), "\n",
    " -", file.path(FIG_DIR, "supp_peak_week_rose_diagrams.png"), "\n")


# ============================================================================
# 03d_spi_scale_comparison.R
# SUPPLEMENTARY ANALYSIS - does SPI-1, SPI-4, or SPI-8 associate most
# strongly with Lassa fever incidence?
#
# 01b_rainfall_interpolation_sensitivity.R checks whether SPI's own
# *distribution* (mean/SD) is sensitive to the rainfall-interpolation
# choice at each scale. It does not compare the three scales against each
# other in terms of association with incidence - this script fills that
# gap, using the same wavelet-coherence + lag-triangulation machinery as
# 03a_spi_wavelet_lag_analysis.R, generalized to loop over SPI_1/4/8.
#
# For each focal state and each SPI scale, this reports:
#   - Consensus_lag_wk    : median of wavelet-phase, CCF, and AIC-guided lag
#   - Annual_Rsq          : mean wavelet coherence in the annual band (48-56
#                           week period), i.e. coupling strength
#   - Prop_Sig_Annual_Band: proportion of that band exceeding the AR(1)
#                           surrogate significance threshold
#
# This directly supports (or challenges) the choice of SPI-4 as the primary
# meteorological driver used in 03a/04/05: if SPI-4's Annual_Rsq is not
# distinguishably higher than SPI-1's or SPI-8's, that should be stated as
# a limitation rather than implied by only ever showing SPI-4 results.
#
# Uses the interpolated ("weekly") series only - this is a scale
# comparison, not an interpolation-sensitivity check (see 01b for that).
#
# Inputs:  outputs/prepared_data.rds
# Outputs: outputs/tables/spi_table1b_scale_comparison.csv
#          outputs/figures/supp_spi_scale_comparison.png
# ============================================================================

source("000_setup.R")

prepped         <- readRDS(file.path(OUT_DIR, "prepared_data.rds"))
weekly          <- prepped$weekly
baseline_params <- prepped$baseline_params

ANNUAL_BAND   <- c(48, 56)   # weeks
N_SURROGATES  <- 1000
MAX_LAG       <- 30
MIN_LAG_OBS   <- 30
SPI_SCALES    <- c(1, 4, 8)

# ---------------------------------------------------------------------------
# Helper functions (mirrors 03a_spi_wavelet_lag_analysis.R; kept local here
# so this script has no hidden dependency on 03a's R session)
# ---------------------------------------------------------------------------
run_wtc <- function(week_idx, clim_var, cases, nrands) {
  x <- cbind(week_idx, clim_var)
  y_log <- log1p(cases)
  y <- cbind(week_idx, scale(y_log)[, 1])
  
  tryCatch({
    biwavelet::wtc(x, y, nrands = nrands, quiet = TRUE)
  }, error = function(e) {
    tryCatch({
      biwavelet::wtc(x, y, nrands = 0, quiet = TRUE)
    }, error = function(e2) NULL)
  })
}

summarize_band <- function(wtc_res) {
  period     <- wtc_res$period
  band_idx   <- which(period >= ANNUAL_BAND[1] & period <= ANNUAL_BAND[2])
  valid_mask <- outer(period, wtc_res$coi, FUN = function(p, c) p <= c)
  
  band_rsq   <- wtc_res$rsq[band_idx, , drop = FALSE]
  band_phase <- wtc_res$phase[band_idx, , drop = FALSE]
  valid      <- valid_mask[band_idx, , drop = FALSE]
  
  mean_rsq    <- mean(band_rsq[valid], na.rm = TRUE)
  mean_phase  <- mean(band_phase[valid], na.rm = TRUE)
  mean_period <- mean(period[band_idx])
  lag_weeks   <- round((mean_phase / (2 * pi)) * mean_period)
  
  list(mean_rsq = mean_rsq, lag_weeks = lag_weeks, valid = valid, band_idx = band_idx)
}

prop_sig_annual_band <- function(wtc_res, band_summary) {
  signif_band <- wtc_res$signif[band_summary$band_idx, , drop = FALSE]
  band_rsq    <- wtc_res$rsq[band_summary$band_idx, , drop = FALSE]
  mean((band_rsq[band_summary$valid] / signif_band[band_summary$valid]) >= 1, na.rm = TRUE)
}

ccf_lag <- function(clim_var, cases, max_lag = MAX_LAG) {
  cc     <- ccf(clim_var, cases, lag.max = max_lag, plot = FALSE)
  lags   <- cc$lag[, 1, 1]
  vals   <- cc$acf[, 1, 1]
  nonneg <- lags >= 0
  lags[nonneg][which.max(vals[nonneg])]
}

aic_lag <- function(clim_var, cases, max_lag = MAX_LAG) {
  aics <- sapply(0:max_lag, function(l) {
    x_lag <- dplyr::lag(clim_var, l)
    ok <- !is.na(x_lag)
    if (sum(ok) < MIN_LAG_OBS) return(Inf)
    fit <- glm(cases[ok] ~ x_lag[ok], family = poisson())
    AIC(fit)
  })
  list(best_lag = (0:max_lag)[which.min(aics)], best_aic = min(aics))
}

# ---------------------------------------------------------------------------
# Run the lag/coherence analysis for one SPI scale, across all focal states
# ---------------------------------------------------------------------------
run_analysis_for_scale <- function(weekly_df, spi_scale) {
  clim_col <- paste0("SPI_", spi_scale)
  
  map_dfr(FOCAL_STATES, function(st) {
    df_state_raw <- weekly_df |> dplyr::filter(State == st) |> dplyr::arrange(week_idx)
    df_state     <- df_state_raw |> tidyr::drop_na(all_of(clim_col), Cases)
    
    wtc_res <- run_wtc(df_state$week_idx, df_state[[clim_col]], df_state$Cases, nrands = N_SURROGATES)
    
    if (is.null(wtc_res)) {
      return(tibble(
        SPI_Scale = spi_scale, State = st,
        Consensus_lag_wk = NA_real_, Annual_Rsq = NA_real_, Prop_Sig_Annual_Band = NA_real_
      ))
    }
    
    band_summary <- summarize_band(wtc_res)
    prop_sig     <- prop_sig_annual_band(wtc_res, band_summary)
    
    ccf_l <- ccf_lag(df_state[[clim_col]], df_state$Cases)
    aic_l <- aic_lag(df_state[[clim_col]], df_state$Cases)
    consensus_lag <- median(c(band_summary$lag_weeks, ccf_l, aic_l$best_lag))
    
    tibble(
      SPI_Scale             = spi_scale,
      State                 = st,
      Consensus_lag_wk      = consensus_lag,
      Annual_Rsq            = round(band_summary$mean_rsq, 3),
      Prop_Sig_Annual_Band  = round(prop_sig, 3)
    )
  })
}

# ---------------------------------------------------------------------------
# Run across SPI-1, SPI-4, SPI-8 and combine
# ---------------------------------------------------------------------------
set.seed(SEED_RESAMPLING)
scale_comparison <- map_dfr(SPI_SCALES, function(sc) run_analysis_for_scale(weekly, sc)) |>
  arrange(State, SPI_Scale)

print(scale_comparison)
write_csv(scale_comparison, file.path(TAB_DIR, "spi_table1b_scale_comparison.csv"))

# ---------------------------------------------------------------------------
# Supplementary figure: coherence strength and consensus lag by SPI scale,
# faceted by state
# ---------------------------------------------------------------------------
library(ggplot2)
library(patchwork) # Required for the '/' syntax

plot_df <- scale_comparison |>
  mutate(SPI_Scale_lab = factor(paste0("SPI-", SPI_Scale), levels = c("SPI-1", "SPI-4", "SPI-8")))

p_rsq <- ggplot(plot_df, aes(x = SPI_Scale_lab, y = Annual_Rsq, fill = SPI_Scale_lab)) +
  geom_col(width = 0.65, show.legend = FALSE) +
  facet_wrap(~ State, nrow = 1) +
  scale_fill_viridis_d(option = "D") +
  labs(x = NULL, y = "Mean wavelet coherence\n(annual band)",
       title = "Coupling strength by SPI scale") +
  theme_minimal(base_size = 11) +
  theme(strip.text = element_text(face = "bold"))

p_lag <- ggplot(plot_df, aes(x = SPI_Scale_lab, y = Consensus_lag_wk, fill = SPI_Scale_lab)) +
  geom_col(width = 0.65, show.legend = FALSE) +
  facet_wrap(~ State, nrow = 1) +
  scale_fill_viridis_d(option = "D") +
  labs(x = NULL, y = "Consensus lag (weeks)",
       title = "Rainfall-to-incidence lag by SPI scale") +
  theme_minimal(base_size = 11) +
  theme(strip.text = element_text(face = "bold"))

# Added plot_annotation tag_levels if formal "A" / "B" labels are required
combined <- (p_rsq / p_lag) + 
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(face = "bold", size = 14))

ggsave(file.path(FIG_DIR, "supp_spi_scale_comparison.png"), combined,
       width = 11, height = 7, dpi = 300, bg = "white")


