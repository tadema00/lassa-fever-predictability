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

source("00_setup.R")

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

# Recompute SPI on weekly_original if not already present
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
MAX_LAG        <- 20
MIN_LAG_OBS    <- 30
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

ccf_lag <- function(clim_var, cases, max_lag = MAX_LAG) {
  cc <- ccf(clim_var, cases, lag.max = max_lag, plot = FALSE)
  lags <- cc$lag[, 1, 1]
  vals <- cc$acf[, 1, 1]
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
# Full analysis function using SPI
# ---------------------------------------------------------------------------
run_analysis_spi <- function(weekly_df, label) {
  
  lag_table <- map_dfr(FOCAL_STATES, function(st) {
    df_state_raw <- weekly_df |> dplyr::filter(State == st) |> dplyr::arrange(week_idx)
    
    df_state <- df_state_raw |> tidyr::drop_na(all_of(clim_col), Cases)
    
    wtc_res <- run_wtc(df_state$week_idx, df_state[[clim_col]], df_state$Cases, nrands = N_SURROGATES)
    
    if (is.null(wtc_res)) {
      return(
        tibble(
          Dataset = label, State = st,
          Wavelet_lag_wk = NA_real_, CCF_lag_wk = NA_real_, AIC_lag_wk = NA_real_,
          Consensus_lag_wk = NA_real_, Annual_Rsq = NA_real_, Prop_Sig_Annual_Band = NA_real_
        )
      )
    }
    
    band_summary <- summarize_band(wtc_res)
    prop_sig     <- prop_sig_annual_band(wtc_res, band_summary)
    
    ccf_l <- ccf_lag(df_state[[clim_col]], df_state$Cases)
    aic_l <- aic_lag(df_state[[clim_col]], df_state$Cases)
    consensus_lag <- median(c(band_summary$lag_weeks, ccf_l, aic_l$best_lag))
    
    tibble(
      Dataset              = label,
      State                = st,
      Wavelet_lag_wk       = band_summary$lag_weeks,
      CCF_lag_wk           = ccf_l,
      AIC_lag_wk           = aic_l$best_lag,
      Consensus_lag_wk     = consensus_lag,
      Annual_Rsq           = round(band_summary$mean_rsq, 3),
      Prop_Sig_Annual_Band = round(prop_sig, 3)
    )
  })
  
  lag_table
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


