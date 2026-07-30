# ============================================================================
# 06_spi_bivariate_state_synchrony.R
# Bivariate wavelet coherence between STATE PAIRS' case-incidence series,
# within each ecological cluster (Table 6), adapted for the SPI workflow.
# ============================================================================

source("00_setup.R")

prepped <- readRDS(file.path(OUT_DIR, "prepared_data.rds"))
weekly  <- prepped$weekly

ANNUAL_BAND  <- c(48, 56)
N_SURROGATES <- 1000

get_state_series <- function(st) {
  weekly %>% dplyr::filter(State == st) %>% dplyr::arrange(week_idx)
}

run_wtc_states <- function(df_a, df_b, nrands = N_SURROGATES) {
  stopifnot(nrow(df_a) == nrow(df_b))
  x <- cbind(df_a$week_idx, scale(log1p(df_a$Cases))[, 1])
  y <- cbind(df_b$week_idx, scale(log1p(df_b$Cases))[, 1])
  biwavelet::wtc(x, y, nrands = nrands, quiet = TRUE)
}

summarize_band_pair <- function(wtc_res) {
  period      <- wtc_res$period
  band_idx    <- which(period >= ANNUAL_BAND[1] & period <= ANNUAL_BAND[2])
  valid_mask  <- outer(period, wtc_res$coi, FUN = function(p, c) p <= c)
  
  band_rsq    <- wtc_res$rsq[band_idx, , drop = FALSE]
  band_phase  <- wtc_res$phase[band_idx, , drop = FALSE]
  band_signif <- wtc_res$signif[band_idx, , drop = FALSE]
  valid       <- valid_mask[band_idx, , drop = FALSE]
  
  mean_rsq    <- mean(band_rsq[valid], na.rm = TRUE)
  mean_phase  <- mean(band_phase[valid], na.rm = TRUE)
  mean_period <- mean(period[band_idx])
  lag_weeks   <- round((mean_phase / (2 * pi)) * mean_period)
  prop_sig    <- mean((band_rsq[valid] / band_signif[valid]) >= 1, na.rm = TRUE)
  
  list(mean_rsq = mean_rsq, lag_weeks = lag_weeks, prop_sig = prop_sig)
}

cluster_pairs <- function(states) combn(states, 2, simplify = FALSE)

run_pair_synchrony <- function(pairs, label) {
  map_dfr(pairs, function(pair) {
    df_a <- get_state_series(pair[1])
    df_b <- get_state_series(pair[2])
    wtc_res <- run_wtc_states(df_a, df_b)
    summ <- summarize_band_pair(wtc_res)
    tibble(
      Cluster              = label,
      State_A              = pair[1],
      State_B              = pair[2],
      Annual_Rsq           = round(summ$mean_rsq, 3),
      Prop_Sig_Annual_Band = round(summ$prop_sig, 3),
      Lag_weeks_A_vs_B     = summ$lag_weeks
    )
  })
}

run_cluster_synchrony <- function(states, cluster_label) {
  run_pair_synchrony(cluster_pairs(states), cluster_label)
}

CROSS_CLUSTER_PAIRS <- list(
  c("Ondo", "Bauchi"),
  c("Edo", "Taraba"),
  c("Ebonyi", "Plateau")
)

set.seed(SEED_RESAMPLING)
southern_synchrony <- run_cluster_synchrony(SOUTHERN_CLUSTER, "Southern Forest")
set.seed(SEED_RESAMPLING)
northern_synchrony <- run_cluster_synchrony(NORTHERN_CLUSTER, "Northern Savannah")
set.seed(SEED_RESAMPLING)
cross_synchrony    <- run_pair_synchrony(CROSS_CLUSTER_PAIRS, "Cross-cluster (control)")

table6_synchrony <- bind_rows(southern_synchrony, northern_synchrony, cross_synchrony)
print(table6_synchrony)
write_csv(table6_synchrony, file.path(TAB_DIR, "spi_table6_bivariate_state_synchrony.csv"))

within_cluster_rsq <- c(southern_synchrony$Annual_Rsq, northern_synchrony$Annual_Rsq)
cross_cluster_rsq  <- cross_synchrony$Annual_Rsq

summary_comparison <- tibble(
  Group           = c("Within-cluster (n=6 pairs)", "Cross-cluster control (n=3 pairs)"),
  Mean_Annual_Rsq = round(c(mean(within_cluster_rsq), mean(cross_cluster_rsq)), 3),
  Min_Annual_Rsq  = round(c(min(within_cluster_rsq), min(cross_cluster_rsq)), 3),
  Max_Annual_Rsq  = round(c(max(within_cluster_rsq), max(cross_cluster_rsq)), 3)
)
print(summary_comparison)
write_csv(summary_comparison, file.path(TAB_DIR, "spi_table6b_within_vs_cross_cluster_summary.csv"))

gap <- mean(within_cluster_rsq) - mean(cross_cluster_rsq)
cat(sprintf(
  "\nMean within-cluster Annual_Rsq = %.3f vs. mean cross-cluster Annual_Rsq = %.3f (gap = %.3f).\n",
  mean(within_cluster_rsq), mean(cross_cluster_rsq), gap
))
if (gap > 0.10) {
  cat("Within-cluster synchrony is notably higher than cross-cluster -- consistent with a real,\n",
      "cluster-specific synchrony finding rather than just shared national seasonality.\n")
} else {
  cat("Within-cluster and cross-cluster synchrony are similar in magnitude -- the honest framing\n",
      "is likely nationwide seasonal coupling (already shown by the Rayleigh test), not something\n",
      "specific to ecological cluster membership. Report Table 6 alongside this comparison rather\n",
      "than presenting within-cluster coherence alone as evidence of cluster-specific synchrony.\n")
}

cat("\nSPI-based bivariate state-pair wavelet synchrony analysis complete (Table 6).\n")

