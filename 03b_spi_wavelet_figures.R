# ============================================================================
# 03b_spi_wavelet_figures.R
# Figures 4 & 5: wavelet coherence panels for the Southern Forest Cluster
# (Ondo, Ebonyi, Edo) and Northern Savannah Cluster (Bauchi, Taraba, Plateau)
# using Standardized Precipitation Index (SPI-4) as the climatic driver.
#
# `weekly` already carries SPI_1/4/8 from 01a, so this stage has no
# baseline_params dependency of its own.
# ============================================================================

source("00_setup.R")

prepped <- readRDS(file.path(OUT_DIR, "prepared_data.rds"))
weekly  <- prepped$weekly

ANNUAL_BAND  <- c(48, 56)
PERIOD_TICKS <- c(2, 4, 8, 16, 32, 52, 104)
SPI_SCALE    <- 4
clim_col     <- paste0("SPI_", SPI_SCALE)

wtc_for_state <- function(state) {
  # FIX: Rain_acc_4/Rain_acc_8 (and therefore SPI_4/SPI_8) are NA for each
  # state's first few weeks (rollapply(..., align = "right", fill = NA) has
  # no full window yet). Feeding those NAs into biwavelet::wtc() poisons the
  # entire FFT with NaN, so plot.wtc()'s range(rsq, na.rm = TRUE) collapses to
  # Inf/-Inf ("invalid z limits"). Stage 03a already drop_na()s before calling
  # wtc(); do the same here. Since the NAs are only a short leading run per
  # state, dropping them just starts the series a few weeks later - it stays
  # contiguous and evenly spaced, so the wavelet transform is unaffected.
  df_state <- weekly %>%
    filter(State == state) %>%
    arrange(week_idx) %>%
    tidyr::drop_na(all_of(clim_col), Cases)
  x <- cbind(df_state$week_idx, df_state[[clim_col]])
  y <- cbind(df_state$week_idx, scale(log1p(df_state$Cases))[, 1])
  biwavelet::wtc(x, y, nrands = 300, quiet = TRUE)
}

add_period_axis <- function() {
  axis(2, at = log2(PERIOD_TICKS), labels = PERIOD_TICKS, las = 1)
}

plot_cluster_panel <- function(states, panel_letters, filename, band = ANNUAL_BAND) {
  tiff(file.path(FIG_DIR, filename), width = 13, height = 9, units = "in",
       res = 600, compression = "lzw")
  on.exit(dev.off(), add = TRUE)
  
  layout(matrix(c(1, 2,
                  3, 3), nrow = 2, byrow = TRUE))
  par(mar = c(4, 4, 2, 1))
  
  wtc_list <- lapply(states, wtc_for_state)
  
  for (i in seq_along(states)) {
    wtc_obj <- wtc_list[[i]]
    
    plot(wtc_obj,
         plot.phase = TRUE,
         plot.coi = TRUE, lwd.coi = 2, col.coi = "black",
         lwd.sig = 2,
         plot.cb = FALSE,
         yaxt = "n",
         main = sprintf("(%s) %s", panel_letters[i], states[i]),
         xlab = "Week", ylab = "Period (weeks)")
    
    add_period_axis()
    abline(h = log2(band), lty = 2, col = "red", lwd = 2)
  }
}

# Figure 4 (Southern Forest Cluster)
plot_cluster_panel(c("Ondo", "Ebonyi", "Edo"), c("A", "B", "C"),
                   "fig4_spi_wavelet_coherence_southern.tiff")

# Figure 5 (Northern Savannah Cluster)
plot_cluster_panel(c("Bauchi", "Taraba", "Plateau"), c("A", "B", "C"),
                   "fig5_spi_wavelet_coherence_northern.tiff")

cat("SPI wavelet coherence cluster figures (Figs 4 & 5) generated successfully.\n")
