# ============================================================================
# 07_spi_bivariate_synchrony_figures.R
# Companion figures for 06_spi_bivariate_state_synchrony.R: cross-wavelet power
# spectrum + phase-difference panels per state pair, adapted for the SPI
# workflow consistency.
# ============================================================================
source("00_setup.R")

prepped <- readRDS(file.path(OUT_DIR, "prepared_data.rds"))
weekly  <- prepped$weekly

ANNUAL_BAND  <- c(48, 56)
PERIOD_TICKS <- c(2, 4, 8, 16, 32, 52, 104)   # Grinsted et al. (2004) convention; 52 = annual

get_state_series <- function(st) {
  weekly %>% dplyr::filter(State == st) %>% dplyr::arrange(week_idx)
}

xwt_for_pair <- function(st_a, st_b, sig_level = 0.95) {
  df_a <- get_state_series(st_a)
  df_b <- get_state_series(st_b)
  stopifnot(nrow(df_a) == nrow(df_b))
  x <- cbind(df_a$week_idx, scale(log1p(df_a$Cases))[, 1])
  y <- cbind(df_b$week_idx, scale(log1p(df_b$Cases))[, 1])
  biwavelet::xwt(x, y, sig.level = sig_level)
}

add_period_axis <- function() {
  axis(2, at = log2(PERIOD_TICKS), labels = PERIOD_TICKS, las = 1)
}

annual_band_phase_series <- function(xwt_res, band = ANNUAL_BAND) {
  period     <- xwt_res$period
  band_idx   <- which(period >= band[1] & period <= band[2])
  phase_band <- xwt_res$phase[band_idx, , drop = FALSE]
  apply(phase_band, 2, function(col) {
    atan2(mean(sin(col), na.rm = TRUE), mean(cos(col), na.rm = TRUE))
  }) * (180 / pi)
}

plot_synchrony_panel <- function(states, filename, band = ANNUAL_BAND) {
  pairs <- combn(states, 2, simplify = FALSE)
  n_pairs <- length(pairs)
  
  tiff(file.path(FIG_DIR, filename), width = 13, height = 4 * n_pairs, units = "in",
       res = 600, compression = "lzw")
  on.exit(dev.off(), add = TRUE)
  
  layout(matrix(seq_len(2 * n_pairs), nrow = n_pairs, byrow = TRUE),
         widths = c(1.3, 1))
  par(mar = c(4, 4, 2.5, 1))
  
  for (i in seq_along(pairs)) {
    st_a <- pairs[[i]][1]; st_b <- pairs[[i]][2]
    xwt_res <- xwt_for_pair(st_a, st_b)
    
    plot(xwt_res,
         plot.phase = TRUE,
         plot.coi = TRUE, lwd.coi = 2, col.coi = "black",
         lwd.sig = 2,
         plot.cb = (i == n_pairs),
         yaxt = "n",
         main = sprintf("%s vs %s: cross-wavelet power", st_a, st_b),
         xlab = "Week", ylab = "Period (weeks)")
    add_period_axis()
    abline(h = log2(band), lty = 2, col = "red", lwd = 2)
    
    phase_deg <- annual_band_phase_series(xwt_res, band)
    week_idx  <- get_state_series(st_a)$week_idx
    plot(week_idx, phase_deg, type = "l", lwd = 1.5, col = "#2874A6",
         ylim = c(-180, 180),
         main = sprintf("%s vs %s: annual-band phase difference", st_a, st_b),
         xlab = "Week", ylab = "Phase difference (deg)")
    abline(h = 0, lty = 2, col = "grey40")
    abline(h = c(-180, 180), lty = 3, col = "grey70")
  }
}

plot_synchrony_panel(SOUTHERN_CLUSTER, "fig8_spi_bivariate_synchrony_southern.tiff")
plot_synchrony_panel(NORTHERN_CLUSTER, "fig9_spi_bivariate_synchrony_northern.tiff")

cat("SPI-based bivariate synchrony figures complete: cross-wavelet power + phase-difference\n",
    "panels per state pair, Southern (fig8) and Northern (fig9) clusters.\n")



