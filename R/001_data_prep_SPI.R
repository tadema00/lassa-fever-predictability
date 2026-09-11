# ============================================================================
# 01a_data_prep_SPI.R
# Multi-Scale Standardized Precipitation Index (SPI-1, SPI-4, SPI-8) pipeline
# SPI was computed following McKee et al. (1993), assuming positive rainfall
# totals are Gamma distributed, with maximum likelihood parameter estimation
# as recommended by Guttman (1999) and Stagge et al. (2015). The 2018-2023
# period served as the climatological reference period for SPI
# standardization. To increase the reference sample size to approximately
# 42 observations while preserving seasonal homogeneity, a +/-3
# epidemiological-week moving window was adopted.
# Saves `baseline_params` into prepared_data.rds so that downstream scripts
# (01b, 03a) can recompute SPI from the same fitted baseline without
# re-running this script or refitting the Gamma distributions.
# ============================================================================

source("000_setup.R")

# ---------------------------------------------------------------------------
# 1. Reshape wide weekly data into long State/Epi_Week/Year panel
# ---------------------------------------------------------------------------
weekly_data <- readxl::read_excel(WEEKLY_XLSX)
glimpse(weekly_data)

# Read annual data only once
annual_data <- readxl::read_excel(ANNUAL_XLSX)
glimpse(annual_data)

cases_long <- weekly_data %>%
  dplyr::select(years, weeks, ends_with("_cases")) %>%
  pivot_longer(cols = ends_with("_cases"), names_to = "State", values_to = "Cases") %>%
  mutate(State = str_to_title(str_remove(State, "_cases")))

rain_long <- weekly_data %>%
  dplyr::select(years, weeks, ends_with("_rain_mm")) %>%
  pivot_longer(cols = ends_with("_rain_mm"), names_to = "State", values_to = "Rainfall_mm") %>%
  mutate(State = str_to_title(str_remove(State, "_weekly_rain_mm")))

weekly_raw <- cases_long %>%
  mutate(across(c(years, weeks), as.integer)) %>%
  inner_join(
    rain_long %>% mutate(across(c(years, weeks), as.integer)),
    by = c("years", "weeks", "State")
  ) %>%
  rename(Year = years, week_idx = weeks) %>%
  mutate(Epi_Week = ((week_idx - 1) %% 52) + 1) %>%
  filter(State %in% FOCAL_STATES) %>%
  arrange(State, Year, week_idx)

# Check that all focal states are present
missing_states <- setdiff(FOCAL_STATES, unique(weekly_raw$State))
if (length(missing_states) > 0) {
  stop("The following focal states are missing from the weekly data: ",
       paste(missing_states, collapse = ", "))
}

# ---------------------------------------------------------------------------
# 2. Fill short case-reporting gaps (<=3 weeks, linear interpolation);
#    longer absences -> 0 ("true epidemiological silence")
# ---------------------------------------------------------------------------
fill_short_gaps <- function(x, max_gap = 3) {
  n <- length(x)
  is_na <- is.na(x)
  if (!any(is_na)) return(x)
  r <- rle(is_na)
  idx_end <- cumsum(r$lengths)
  idx_start <- idx_end - r$lengths + 1
  out <- x
  for (i in seq_along(r$lengths)) {
    if (r$values[i] && r$lengths[i] <= max_gap) {
      s <- idx_start[i]; e <- idx_end[i]
      lo <- if (s > 1) x[s - 1] else NA
      hi <- if (e < n) x[e + 1] else NA
      if (!is.na(lo) && !is.na(hi)) {
        out[s:e] <- approx(x = c(s - 1, e + 1), y = c(lo, hi), xout = s:e)$y
      }
    } else if (r$values[i]) {
      s <- idx_start[i]; e <- idx_end[i]
      out[s:e] <- 0
    }
  }
  out
}

weekly <- weekly_raw %>%
  group_by(State) %>%
  arrange(Year, week_idx) %>%
  mutate(Cases = fill_short_gaps(Cases), week_idx = row_number()) %>%
  ungroup()

# ---------------------------------------------------------------------------
# 3. Rainfall interpolation & multi-scale accumulated precipitation (1, 4, 8 weeks)
# ---------------------------------------------------------------------------
weekly <- weekly %>%
  group_by(State) %>%
  arrange(Year, week_idx) %>%
  mutate(Rainfall_mm = zoo::na.approx(Rainfall_mm, na.rm = FALSE)) %>%
  mutate(
    Rainfall_mm = zoo::na.locf(Rainfall_mm, na.rm = FALSE),
    Rainfall_mm = zoo::na.locf(Rainfall_mm, fromLast = TRUE, na.rm = FALSE)
  ) %>%
  mutate(
    Rain_acc_1 = Rainfall_mm,
    Rain_acc_4 = rollapply(Rainfall_mm, width = 4, FUN = sum, align = "right", fill = NA),
    Rain_acc_8 = rollapply(Rainfall_mm, width = 8, FUN = sum, align = "right", fill = NA)
  ) %>%
  ungroup()

# ---------------------------------------------------------------------------
# 4. Precomputed Baseline Parameter Table (+/-3 Weeks / ~42 Observations Window)
#    Includes Gamma MLE/MoM fitting and Anderson-Darling goodness-of-fit testing.
# ---------------------------------------------------------------------------
build_baseline_parameters <- function(df, scales = c(1, 4, 8), ref_start = 2018, ref_end = 2023) {
  ref_data <- df %>% filter(Year >= ref_start & Year <= ref_end)
  states <- unique(df$State)
  weeks <- 1:52
  
  baseline_list <- vector("list", length(states) * length(scales) * length(weeks))
  idx <- 1
  
  for (st in states) {
    df_st <- ref_data %>% filter(State == st)
    
    for (sc in scales) {
      acc_col <- paste0("Rain_acc_", sc)
      
      for (w in weeks) {
        raw_neigh <- seq(w - 3, w + 3)
        neighbor_weeks <- ((raw_neigh - 1) %% 52) + 1
        
        sample_vals <- df_st %>%
          filter(Epi_Week %in% neighbor_weeks) %>%
          pull(.data[[acc_col]])
        
        sample_vals <- sample_vals[!is.na(sample_vals)]
        n_total <- length(sample_vals)
        
        if (n_total == 0) {
          baseline_list[[idx]] <- tibble(
            State = st, Scale = sc, Epi_Week = w,
            shape = NA_real_, rate = NA_real_, q = NA_real_,
            n_total = 0, n_zero = 0, method = "none",
            ad_stat = NA_real_, ad_pval = NA_real_,
            empirical_fn = list(NULL)
          )
          idx <- idx + 1
          next
        }
        
        zeros <- sample_vals == 0
        n_zeros <- sum(zeros)
        q_val <- n_zeros / n_total
        non_zero <- sample_vals[!zeros]
        
        shape_est <- NA_real_
        rate_est <- NA_real_
        fit_method <- "none"
        ad_stat <- NA_real_
        ad_pval <- NA_real_
        emp_fn <- if (length(non_zero) > 0) ecdf(non_zero) else NULL
        
        if (length(non_zero) >= 2 && sd(non_zero) > 0) {
          fit_res <- tryCatch({
            fit <- MASS::fitdistr(non_zero, densfun = "gamma", lower = c(0.001, 0.001))
            list(shape = fit$estimate["shape"], rate = fit$estimate["rate"])
          }, error = function(e) NULL)
          
          if (!is.null(fit_res)) {
            shape_est <- fit_res$shape
            rate_est <- fit_res$rate
            fit_method <- "MLE"
          } else {
            mean_nz <- mean(non_zero)
            var_nz <- var(non_zero)
            if (var_nz > 0 && mean_nz > 0) {
              shape_est <- (mean_nz^2) / var_nz
              scale_mom <- var_nz / mean_nz
              rate_est <- 1.0 / scale_mom
              fit_method <- "MoM"
            }
          }
          
          if (!is.na(shape_est) && shape_est > 0 && rate_est > 0) {
            ad_test <- tryCatch({
              goftest::ad.test(non_zero, "pgamma", shape = shape_est, rate = rate_est)
            }, error = function(e) NULL)
            
            if (!is.null(ad_test)) {
              ad_stat <- ad_test$statistic
              ad_pval <- ad_test$p.value
            }
          }
        }
        
        baseline_list[[idx]] <- tibble(
          State = st, Scale = sc, Epi_Week = w,
          shape = shape_est, rate = rate_est, q = q_val,
          n_total = n_total, n_zero = n_zeros, method = fit_method,
          ad_stat = ad_stat, ad_pval = ad_pval,
          empirical_fn = list(emp_fn)
        )
        idx <- idx + 1
      }
    }
  }
  
  bind_rows(baseline_list)
}

cat("\nBuilding baseline parameters with +/-3 week window...\n")
baseline_params <- build_baseline_parameters(weekly)

# Export parameter table without list columns
baseline_export <- baseline_params %>% select(-empirical_fn)
write_csv(baseline_export, file.path(TAB_DIR, "spi_baseline_parameters.csv"))

# Fitting convergence & failure summary
fitting_convergence <- baseline_params %>%
  group_by(Scale, State, method) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(Scale, State) %>%
  mutate(pct = 100 * count / sum(count))

write_csv(fitting_convergence, file.path(TAB_DIR, "spi_fitting_convergence_by_state.csv"))

# ---------------------------------------------------------------------------
# 5. Vectorized SPI Calculation via McKee (1993) Formulation
#    Using smooth empirical CDF fallback via ecdf()
# ---------------------------------------------------------------------------
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

weekly <- compute_spi_from_baseline(weekly, baseline_params, 1)
weekly <- compute_spi_from_baseline(weekly, baseline_params, 4)
weekly <- compute_spi_from_baseline(weekly, baseline_params, 8)

# Export missingness report by state and scale
missing_spi_report <- weekly %>%
  group_by(State) %>%
  summarise(
    Missing_SPI1 = sum(is.na(SPI_1)),
    Missing_SPI4 = sum(is.na(SPI_4)),
    Missing_SPI8 = sum(is.na(SPI_8)),
    .groups = "drop"
  )
write_csv(missing_spi_report, file.path(TAB_DIR, "spi_missing_counts_by_state.csv"))

# ---------------------------------------------------------------------------
# 6. Comprehensive Validation & Target Benchmarking (2018-2023)
# ---------------------------------------------------------------------------
ref_subset <- weekly %>% filter(Year <= 2023)

state_spi_validation <- ref_subset %>%
  group_by(State) %>%
  summarise(
    SPI1_Mean = mean(SPI_1, na.rm = TRUE), SPI1_SD = sd(SPI_1, na.rm = TRUE), SPI1_Pass = abs(SPI1_Mean) < 0.1 && abs(SPI1_SD - 1.0) < 0.1,
    SPI4_Mean = mean(SPI_4, na.rm = TRUE), SPI4_SD = sd(SPI_4, na.rm = TRUE), SPI4_Pass = abs(SPI4_Mean) < 0.1 && abs(SPI4_SD - 1.0) < 0.1,
    SPI8_Mean = mean(SPI_8, na.rm = TRUE), SPI8_SD = sd(SPI_8, na.rm = TRUE), SPI8_Pass = abs(SPI8_Mean) < 0.1 && abs(SPI8_SD - 1.0) < 0.1,
    .groups = "drop"
  )

write_csv(state_spi_validation, file.path(TAB_DIR, "state_spi_validation_targets.csv"))

# ---------------------------------------------------------------------------
# 7. Global Diagnostic Summary Table (Median AD p-value & Convergence)
# ---------------------------------------------------------------------------
diagnostic_summary <- tibble(
  Scale = character(), Mean = numeric(), SD = numeric(),
  Skewness = numeric(), Kurtosis = numeric(),
  `Median AD p-value` = numeric(), `% AD Pass (p>0.05)` = numeric(), `MLE success (%)` = numeric()
)

for (sc in c(1, 4, 8)) {
  col_name <- paste0("SPI_", sc)
  vals <- ref_subset[[col_name]]
  vals <- vals[!is.na(vals)]
  
  base_sc <- baseline_params %>% filter(Scale == sc)
  med_ad_p <- median(base_sc$ad_pval, na.rm = TRUE)
  pct_pass <- 100 * sum(base_sc$ad_pval > 0.05, na.rm = TRUE) / sum(!is.na(base_sc$ad_pval))
  mle_succ <- 100 * sum(base_sc$method == "MLE", na.rm = TRUE) / nrow(base_sc)
  
  diagnostic_summary <- bind_rows(diagnostic_summary, tibble(
    Scale = paste0("SPI-", sc),
    Mean = round(mean(vals), 2),
    SD = round(sd(vals), 2),
    Skewness = round(e1071::skewness(vals), 2),
    Kurtosis = round(e1071::kurtosis(vals), 2),
    `Median AD p-value` = round(med_ad_p, 2),
    `% AD Pass (p>0.05)` = round(pct_pass, 1),
    `MLE success (%)` = round(mle_succ, 1)
  ))
}

write_csv(diagnostic_summary, file.path(TAB_DIR, "spi_global_diagnostic_summary.csv"))

# ---------------------------------------------------------------------------
# 8. ggplot2 Publication-Quality Q-Q Plots with Confidence Bands
# ---------------------------------------------------------------------------
dir.create(file.path(FIG_DIR, "supp_qq"), showWarnings = FALSE, recursive = TRUE)

for (st in FOCAL_STATES) {
  for (sc in c(1, 4, 8)) {
    col_name <- paste0("SPI_", sc)
    df_qq <- ref_subset %>% filter(State == st) %>% select(val = all_of(col_name)) %>% drop_na()
    
    p <- ggplot(df_qq, aes(sample = val)) +
      stat_qq(color = "#2c3e50", alpha = 0.6) +
      stat_qq_line(color = "#e74c3c", linetype = "dashed", linewidth = 1) +
      labs(
        title = sprintf("Normal Q-Q Plot: SPI-%d (%s)", sc, st),
        subtitle = "Climatological Reference Period (2018-2023)",
        x = "Theoretical Quantiles",
        y = "Sample Quantiles"
      ) +
      theme_bw(base_size = 12) +
      theme(
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(color = "gray40", size = 10)
      )
    
    ggsave(file.path(FIG_DIR, "supp_qq", sprintf("SPI_QQ_SPI%d_%s.pdf", sc, st)), p, width = 6, height = 5)
  }
}

# ---------------------------------------------------------------------------
# 9. Annual Panel & Integrity Checks
# ---------------------------------------------------------------------------
annual_raw <- readxl::read_excel(ANNUAL_XLSX) %>% rename_with(str_trim)

annual <- annual_raw %>%
  mutate(across(where(is.numeric), ~ replace_na(.x, 0))) %>%
  filter(State != "Total", State != "Jos") %>%
  mutate(State = if_else(State == "Fct", "Federal Capital Territory", State))

year_cols <- str_subset(names(annual), "^Cumulative Confirmed_")

annual_long <- annual %>%
  select(State, all_of(year_cols)) %>%
  pivot_longer(-State, names_to = "Year", values_to = "Cases") %>%
  mutate(
    Year  = as.integer(str_extract(Year, "\\d{4}")),
    Cases = replace_na(Cases, 0)
  )

annual_cumulative <- annual %>% select(State, Total)
total_row <- annual_raw %>% filter(`S/N` == "Total") %>% slice(1)

integrity_check <- annual_long %>%
  group_by(Year) %>%
  summarise(computed_total = sum(Cases), .groups = "drop") %>%
  mutate(
    reported_total = sapply(Year, function(y) {
      col <- paste0("Cumulative Confirmed_", y)
      if (col %in% names(total_row)) total_row[[col]] else NA_real_
    })
  )

write_csv(integrity_check, file.path(TAB_DIR, "national_total_integrity_check.csv"))
integrity_check
# ---------------------------------------------------------------------------
# Save prepared data (includes weekly_original for sensitivity analysis)
# ---------------------------------------------------------------------------
saveRDS(
  list(
    weekly           = weekly,
    weekly_original  = weekly_raw,   # uninterpolated rainfall, pre-SPI
    baseline_params  = baseline_params,
    annual_long      = annual_long,
    annual_cumulative = annual_cumulative
  ),
  file.path(OUT_DIR, "prepared_data.rds")
)
cat("\nData preparation completed\n")


# ============================================================================
# 01b_rainfall_interpolation_sensitivity.R
# Robustness check: do the SPI-based lag and coherence estimates depend on
# the rainfall-interpolation procedure used in 01a_data_prep_SPI.R?
#
# Compares `weekly_original` (rainfall as reported, missing weeks left as NA
# and dropped case-wise) against `weekly` (linearly interpolated + LOCF/NOCB-
# completed) on SPI distribution properties (mean, SD) at each scale
# (SPI-1, SPI-4, SPI-8).
#
# Reads `baseline_params` and `compute_spi_from_baseline()` from
# prepared_data.rds / redefines them locally, so this script has no hidden
# dependency on 01a's R session and can be run independently.
# ============================================================================

# (This section is already included above; for clarity, we run it here as a separate block.)
# If you split into two files, move everything below this line into 01b.R

source("00_setup.R")   # if run separately, but we are in the same script

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

