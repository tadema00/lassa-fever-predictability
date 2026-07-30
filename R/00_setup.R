# ============================================================================
# 00_setup.R
# Lassa fever Early Warning System pipeline (SPI version)
#
# Packages, global constants, file paths, and shared helper functions used
# by every other script in the pipeline. Source this file at the top of
# every stage script; do not run it standalone.
# ============================================================================

pkgs <- c(
  "tidyverse",   # dplyr, tidyr, ggplot2, purrr, readr, stringr
  "readxl",      # read .xlsx surveillance/climate data
  "openxlsx",    # write .xlsx outputs
  "lubridate",   # date/epiweek handling
  "sf",          # spatial vector data (GADM shapefile)
  "spdep",       # Queen contiguity weights + Global Moran's I
  "biwavelet",   # continuous wavelet transform + wavelet coherence + AR(1) MC surrogate test
  "forecast",    # SARIMA/SARIMAX (Arima with xreg), MASE-style accuracy()
  "MASS",        # fitdistr (Gamma MLE), glm.nb
  "pscl",        # zeroinfl (Zero-Inflated Negative Binomial)
  "lmtest",      # lrtest, waldtest for ZINB vs NB nesting
  "signal",      # butter() + filtfilt() zero-phase Butterworth bandpass
  "zoo",         # na.approx (linear interpolation) + na.locf (LOCF/NOCB fill) for rainfall gaps
  "geosphere",   # great-circle (haversine) distance between state centroids
  "circular",    # circular statistics: Rayleigh test for peak-month consistency
  "boot",        # bootstrap resampling for alert sensitivity/PPV CIs
  "tseries",     # adf.test()/kpss.test() - stationarity checks before SARIMAX/ZINB fitting (stage 04)
  "viridis",     # perceptually uniform palettes for choropleths
  "patchwork",   # combine ggplot panels (Fig 3-7 style multi-panel layouts)
  "scico",       # additional diverging/sequential palettes
  "cowplot",     # cowplot::get_legend() for the shared choropleth legend (stage 02)
  "e1071",       # skewness()/kurtosis() for the SPI diagnostic summary (stage 01a)
  "goftest"      # Anderson-Darling goodness-of-fit test for the Gamma fit (stage 01a)
)

new_pkgs <- pkgs[!(pkgs %in% installed.packages()[, "Package"])]
if (length(new_pkgs) > 0) install.packages(new_pkgs, repos = "https://cloud.r-project.org")

invisible(lapply(pkgs, library, character.only = TRUE))

# ---------------------------------------------------------------------------
# Namespace-conflict guard: `signal::filter()` and `MASS::select()` both
# load *after* tidyverse in the list above and silently mask
# `dplyr::filter()`/`dplyr::select()`. Those replacements don't support
# non-standard evaluation on data-frame columns, so every `filter(State ==
# ...)` / `select(State, ...)` call downstream fails with errors like
# "argument x is missing" or "unused arguments" - the failure looks like a
# bug in the calling script, but it's actually the wrong function being
# called. Pin the tidyverse verbs explicitly so load order never matters.
# ---------------------------------------------------------------------------
filter <- dplyr::filter
select <- dplyr::select
lag    <- dplyr::lag   # MASS/stats also both define lag(); same guard

# ---------------------------------------------------------------------------
# Random seeds
#
# Two seeds are used across this pipeline, for different purposes:
#   SEED_GLOBAL     - general reproducibility (set once here; governs any
#                      randomness in data prep / model fitting stages).
#   SEED_RESAMPLING - Monte Carlo / bootstrap procedures (wavelet coherence
#                      AR(1) surrogates, bivariate synchrony surrogates,
#                      alert-sensitivity bootstrap CIs). Set locally inside
#                      each of those scripts immediately before the relevant
#                      resampling call.
#
# NOTE: the manuscript (Section 2.9) currently states "a fixed seed (42)
# ensured reproducibility." That should be revised to describe both seeds
# explicitly, since the Monte Carlo / bootstrap stages use SEED_RESAMPLING,
# not SEED_GLOBAL.
# ---------------------------------------------------------------------------
SEED_GLOBAL     <- 42
SEED_RESAMPLING <- 2026

set.seed(SEED_GLOBAL)

# ---------------------------------------------------------------------------
# Global constants (shared across all scripts in this pipeline)
# ---------------------------------------------------------------------------

FOCAL_STATES     <- c("Ondo", "Edo", "Ebonyi", "Bauchi", "Taraba", "Plateau")
SOUTHERN_CLUSTER <- c("Ondo", "Edo", "Ebonyi")       # Southern Forest Cluster
NORTHERN_CLUSTER <- c("Bauchi", "Taraba", "Plateau") # Northern Savannah Cluster

cluster_of <- function(state) {
  dplyr::case_when(
    state %in% SOUTHERN_CLUSTER ~ "Southern Forest",
    state %in% NORTHERN_CLUSTER ~ "Northern Savannah",
    TRUE ~ "Other (not in study)"
  )
}
# ---------------------------------------------------------------------------
# Project directories
#
# Raw input datasets should be placed in:
#   data/raw/
#
# Intermediate processed datasets may be written to:
#   data/processed/
#
# All generated figures, tables, and .rds files are written to:
#   outputs/
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# Project directories
# ---------------------------------------------------------------------------

RAW_DATA_DIR  <- file.path("data", "raw")
PROC_DATA_DIR <- file.path("data", "processed")

dir.create(RAW_DATA_DIR,  recursive = TRUE, showWarnings = FALSE)
dir.create(PROC_DATA_DIR, recursive = TRUE, showWarnings = FALSE)

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(RDS_DIR, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Input data
# ---------------------------------------------------------------------------

SHAPEFILE   <- file.path(RAW_DATA_DIR, "gadm41_NGA_1.shp")
ANNUAL_XLSX <- file.path(RAW_DATA_DIR, "annual_data.xlsx")
WEEKLY_XLSX <- file.path(RAW_DATA_DIR, "Cases_rainfal_data.xlsx")


# Choropleth binning (matches the Python script's BINS/BIN_LABELS/COLORS)
CASE_BREAKS <- c(-0.1, 0, 10, 50, 100, 300, Inf)
CASE_LABELS <- c("0", "1 - 10", "11 - 50", "51 - 100", "101 - 300", "> 300")
CASE_COLORS <- c("#F2F2F2", "#F9E79F", "#F3C57A", "#E67E22", "#C0392B", "#641E16")

DATA_NOTE <- paste(
  "National totals reported by NCDC exceed the summed state values in some",
  "years because several low-incidence states contained missing or",
  "suppressed counts in the compiled dataset."
)

# ---------------------------------------------------------------------------
# Small shared helpers
# ---------------------------------------------------------------------------

# Standard MASE against a seasonal-naive benchmark (Hyndman & Koehler 2006),
# computed manually so the seasonal period and training window are explicit
# rather than relying on forecast::accuracy()'s in-sample scaling default.
mase <- function(actual, predicted, train_actual, season = 52) {
  naive_errors <- abs(diff(train_actual, lag = season))
  scale <- mean(naive_errors, na.rm = TRUE)
  mean(abs(actual - predicted), na.rm = TRUE) / scale
}

# Great-circle distance in km between two lon/lat points
great_circle_km <- function(lon1, lat1, lon2, lat2) {
  geosphere::distHaversine(cbind(lon1, lat1), cbind(lon2, lat2)) / 1000
}

cat("Setup complete. Packages loaded, constants defined.\n")

