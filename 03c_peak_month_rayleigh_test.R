# ============================================================================
# 03c_peak_month_rayleigh_test.R
# Table 2: Peak-timing consistency (Rayleigh test)
#
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

source("00_setup.R")
library(circular)   # rayleigh.test(), rho.circular(), mean.circular(), circular()

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
# Base-R `circular` plotting (no ggplot equivalent exists for this).
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

