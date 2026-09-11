# ============================================================================
# 00_setup.R
# Lassa fever Early Warning System pipeline
# Packages, global constants, file paths, and shared helper functions used
# by every other script in the pipeline. Source this file at the top of
# every stage script; do not run it standalone.
# ============================================================================

# ---------------------------------------------------------------------------
# Use 'here' to locate files relative to project root
# ---------------------------------------------------------------------------
if (!requireNamespace("here", quietly = TRUE)) {
  install.packages("here", repos = "https://cloud.r-project.org")
}
library(here)

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

# Install missing packages if any
new_pkgs <- pkgs[!(pkgs %in% installed.packages()[, "Package"])]
if (length(new_pkgs) > 0) {
  message("Installing missing packages: ", paste(new_pkgs, collapse = ", "))
  install.packages(new_pkgs, repos = "https://cloud.r-project.org")
}
invisible(lapply(pkgs, library, character.only = TRUE))


filter <- dplyr::filter
select <- dplyr::select
lag    <- dplyr::lag   # MASS/stats also both define lag(); same guard

# ---------------------------------------------------------------------------
# Random seeds
# Two seeds are used across this pipeline, for different purposes:
#   SEED_GLOBAL     - general reproducibility (set once here; governs any
#                      randomness in data prep / model fitting stages).
#   SEED_RESAMPLING - Monte Carlo / bootstrap procedures (wavelet coherence
#                      AR(1) surrogates, bivariate synchrony surrogates,
#                      alert-sensitivity bootstrap CIs). Set locally inside
#                      each of those scripts immediately before the relevant
#                      resampling call.
# ---------------------------------------------------------------------------
SEED_GLOBAL     <- 42
SEED_RESAMPLING <- 2026

set.seed(SEED_GLOBAL)

# ---------------------------------------------------------------------------
# File paths – all relative to project root via here::here()
# ---------------------------------------------------------------------------
DATA_DIR <- here("data")
OUT_DIR  <- here("outputs")
FIG_DIR  <- file.path(OUT_DIR, "figures")
TAB_DIR  <- file.path(OUT_DIR, "tables")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)

SHAPEFILE   <- file.path(DATA_DIR, "gadm41_NGA_1.shp")
ANNUAL_XLSX <- file.path(DATA_DIR, "annual_data.xlsx")
WEEKLY_XLSX <- file.path(DATA_DIR, "Cases_rainfal_data.xlsx")   # state, epi_week, year, cases, rainfall_mm

# ---------------------------------------------------------------------------
# Check that required input data files exist
# ---------------------------------------------------------------------------
required_files <- c(SHAPEFILE, ANNUAL_XLSX, WEEKLY_XLSX)
missing <- required_files[!file.exists(required_files)]
if (length(missing) > 0) {
  stop(
    "The following required data files are missing:\n",
    paste("  -", missing, collapse = "\n"),
    "\nPlease place them in the 'data/' folder."
  )
}

# ---------------------------------------------------------------------------
# Focal states and clustering
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
# Choropleth binning (matches the Python script's BINS/BIN_LABELS/COLORS)
# ---------------------------------------------------------------------------
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

cat("Setup complete. Packages loaded, constants defined, data files verified.\n")

