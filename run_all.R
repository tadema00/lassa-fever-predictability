# ============================================================================
# run_all.R
# Runs the full Lassa fever SPI early-warning pipeline end to end, in the
# order the stages depend on each other.
#
# Usage:
#   Rscript run_all.R
# or, interactively, source() this file from the project root (the folder
# that contains all the numbered .R scripts, plus data/ and outputs/).
#
# Each stage script sources 00_setup.R itself using a bare relative path
# ("00_setup.R"), so this file - and every stage script - must be run with
# the project root as the working directory.
# ============================================================================

stages <- c(
  "00_setup.R",
  "01a_data_prep_SPI.R",
  "01b_rainfall_interpolation_sensitivity.R",
  "02_spatial_choropleth_moran.R",
  "03a_spi_wavelet_lag_analysis.R",
  "03b_spi_wavelet_figures.R",
  "03c_peak_month_rayleigh_test.R",
  "04_spi_sarimax_zinb_models.R",
  "04b_pooled_nb_cluster_interaction.R",
  "05_spi_rolling_forecast_alerts.R",
  "05b_spi_alert_sensitivity_PI.R",
  "06_spi_bivariate_state_synchrony.R",
  "07_spi_bivariate_synchrony_figures.R"
)

for (stage in stages) {
  cat("\n============================================================\n")
  cat("Running:", stage, "\n")
  cat("============================================================\n")
  source(stage, chdir = FALSE)
}

cat("\nAll stages complete. Outputs written to outputs/tables and outputs/figures.\n")
