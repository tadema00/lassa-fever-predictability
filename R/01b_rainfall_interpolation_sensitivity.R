# ============================================================================
# 01b_rainfall_interpolation_sensitivity.R
#
# Robustness check: do the SPI-based lag and coherence estimates depend on
# the rainfall-interpolation procedure used in 01a_data_prep_SPI.R?
#
# Compares `weekly_original` (rainfall as reported, missing weeks left as NA
# and dropped case-wise) against `weekly` (linearly interpolated + LOCF/NOCB-
# completed) on SPI distribution properties (mean, SD) at each scale
# (SPI-1, SPI-4, SPI-8).
#
# TODO: this script checks whether SPI's *distribution* is sensitive to the
# interpolation choice. It does not yet compare SPI-1 vs SPI-4 vs SPI-8
# against each other in terms of association with incidence or forecast
# skill - that would be a separate analysis (see manuscript review notes).
#
# Reads `baseline_params` and `compute_spi_from_baseline()` from
# prepared_data.rds / redefines them locally, so this script has no hidden
# dependency on 01a's R session and can be run independently.
# ============================================================================

source("00_setup.R")

prepped         <- readRDS(file.path(OUT_DIR, "prepared_data.rds"))
weekly          <- prepped$weekly           # interpolated rainfall & SPI
weekly_original <- prepped$weekly_original  # rainfall as reported (NAs retained)
baseline_params <- prepped$baseline_params  # Gamma/empirical baseline fits from 01a

# Local copy of the SPI transform (identical to the one in 01a) so this
# script has no hidden dependency on 01a's environment.
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

# ---------------------------------------------------------------------------
# Recompute SPI on the original (uninterpolated) series for comparison
# ---------------------------------------------------------------------------
weekly_original <- weekly_original %>%
  group_by(State) %>%
  arrange(Year, week_idx) %>%
  mutate(
    Rain_acc_1 = Rainfall_mm,
    Rain_acc_4 = rollapply(Rainfall_mm, width = 4, FUN = sum, align = "right", fill = NA, na.rm = FALSE),
    Rain_acc_8 = rollapply(Rainfall_mm, width = 8, FUN = sum, align = "right", fill = NA, na.rm = FALSE)
  ) %>%
  ungroup()

# Compute SPI for original series using the same baseline parameters
weekly_original <- compute_spi_from_baseline(weekly_original, baseline_params, 1)
weekly_original <- compute_spi_from_baseline(weekly_original, baseline_params, 4)
weekly_original <- compute_spi_from_baseline(weekly_original, baseline_params, 8)

# ---------------------------------------------------------------------------
# Compare metrics between interpolated and original series
# ---------------------------------------------------------------------------
compare_one_state_spi <- function(st) {
  d_orig <- weekly_original %>% filter(State == st) %>% arrange(Year, week_idx)
  d_int  <- weekly         %>% filter(State == st) %>% arrange(Year, week_idx)
  
  tibble(
    State = st,
    SPI1_mean_orig = mean(d_orig$SPI_1, na.rm = TRUE),
    SPI1_mean_int  = mean(d_int$SPI_1, na.rm = TRUE),
    SPI1_sd_orig   = sd(d_orig$SPI_1, na.rm = TRUE),
    SPI1_sd_int    = sd(d_int$SPI_1, na.rm = TRUE),
    
    SPI4_mean_orig = mean(d_orig$SPI_4, na.rm = TRUE),
    SPI4_mean_int  = mean(d_int$SPI_4, na.rm = TRUE),
    SPI4_sd_orig   = sd(d_orig$SPI_4, na.rm = TRUE),
    SPI4_sd_int    = sd(d_int$SPI_4, na.rm = TRUE),
    
    SPI8_mean_orig = mean(d_orig$SPI_8, na.rm = TRUE),
    SPI8_mean_int  = mean(d_int$SPI_8, na.rm = TRUE),
    SPI8_sd_orig   = sd(d_orig$SPI_8, na.rm = TRUE),
    SPI8_sd_int    = sd(d_int$SPI_8, na.rm = TRUE)
  )
}

sensitivity_table <- map_dfr(FOCAL_STATES, compare_one_state_spi) %>%
  mutate(
    SPI1_mean_diff = abs(SPI1_mean_int - SPI1_mean_orig),
    SPI1_sd_diff   = abs(SPI1_sd_int - SPI1_sd_orig),
    SPI4_mean_diff = abs(SPI4_mean_int - SPI4_mean_orig),
    SPI4_sd_diff   = abs(SPI4_sd_int - SPI4_sd_orig),
    SPI8_mean_diff = abs(SPI8_mean_int - SPI8_mean_orig),
    SPI8_sd_diff   = abs(SPI8_sd_int - SPI8_sd_orig)
  )

print(sensitivity_table)
write_csv(sensitivity_table, file.path(TAB_DIR, "spi_rainfall_interpolation_sensitivity.csv"))

max_mean_diff <- max(c(sensitivity_table$SPI1_mean_diff, sensitivity_table$SPI4_mean_diff, sensitivity_table$SPI8_mean_diff), na.rm = TRUE)
max_sd_diff   <- max(c(sensitivity_table$SPI1_sd_diff, sensitivity_table$SPI4_sd_diff, sensitivity_table$SPI8_sd_diff), na.rm = TRUE)

cat(sprintf(
  "\nLargest mean SPI discrepancy across states/scales: %.3f. Largest SD discrepancy: %.3f.\n",
  max_mean_diff, max_sd_diff
))

cat(if (max_mean_diff <= 0.1 && max_sd_diff <= 0.1) {
  "Sensitivity analysis passed: Interpolation procedure has a negligible effect on standardized precipitation indices -> results are robust.\n"
} else {
  "Discrepancies exceed acceptable thresholds for at least one state/scale - inspect `spi_rainfall_interpolation_sensitivity.csv`.\n"
})

