# Lassa Fever SPI Early Warning Pipeline

R analysis pipeline supporting *"Predictability and Timing of Lassa Fever
Outbreaks in Nigeria: A Wavelet-Informed Early Warning Framework Using the
Standardized Precipitation Index."* The pipeline estimates climate-disease
lags using wavelet phase, cross-correlation, and AIC-guided triangulation;
tests wavelet coherence between the Standardized Precipitation Index (SPI-4)
and Lassa fever incidence; fits SARIMAX and Zero-Inflated Negative Binomial
(ZINB) models; and evaluates a four-tier operational alert system via
rolling-origin forecast validation.

## Requirements

- R >= 4.4.0
- Packages (installed automatically on first run if missing):
  `tidyverse`, `readxl`, `openxlsx`, `lubridate`, `sf`, `spdep`, `biwavelet`,
  `forecast`, `MASS`, `pscl`, `lmtest`, `signal`, `zoo`, `geosphere`,
  `circular`, `boot`, `tseries`, `viridis`, `patchwork`, `scico`, `cowplot`,
  `e1071`, `goftest`

## Setup

1. Clone this repository.
2. Place the required data files in `data/` (see `data/README.md` for the
   exact filenames expected).
3. From the project root:

```r
source("run_all.R")
```

or, from the command line:

```bash
Rscript run_all.R
```

All tables are written to `outputs/tables/` and all figures to
`outputs/figures/`. Intermediate `.rds` objects (e.g. `prepared_data.rds`,
`lag_tables_sensitivity.rds`) are written to `outputs/` and read back in by
downstream scripts, so each stage can also be run independently after the
stages it depends on have been run once.

## Pipeline stages

Scripts are numbered in the order they should run. Each stage sources
`00_setup.R` itself and is independently runnable, provided its listed
dependencies have already been produced.

| Script | Purpose | Depends on | Key outputs |
|---|---|---|---|
| `000_setup.R` | Packages, global constants, file paths, shared helpers (`mase()`, `great_circle_km()`). Sourced by every other script. | - | - |
| `001_data_prep_SPI.R` | Reshapes raw data, fills short case-reporting gaps, computes multi-scale SPI (SPI-1/4/8) via Gamma MLE with empirical-CDF fallback, Q-Q diagnostics. | `00_setup.R` | `prepared_data.rds`, `spi_baseline_parameters.csv`, `spi_global_diagnostic_summary.csv`, Q-Q plots |
| `002_spatial_choropleth.R` | National choropleth maps (annual + cumulative incidence)  | `01a` | `national_annual_choropleth.png`|
| `003_spi_wavelet_lag_analysis.R` | Continuous wavelet transform + AR(1) Monte Carlo surrogate coherence test; three-way lag triangulation (wavelet phase, CCF, AIC) → consensus lag (Table 1). Run on both interpolated and original rainfall series. | `01a` | `lag_tables_sensitivity.rds`, `spi_table1_lag_triangulation_sensitivity.csv` | `| Wavelet coherence panel figures (Figs 4-5) for the Southern Forest and Northern Savannah clusters. | `01a` | figures | `03c_peak_month_rayleigh_test.R` | Peak-timing consistency (Rayleigh test) across the study period (Table 2). | `01a` | `spi_table2_peak_timing_rayleigh_cases.csv`, rose diagrams |
| `004_spi_sarimax_zinb_models.R` | Formal stationarity assessment (ADF/KPSS), SARIMAX with SPI-4 at the consensus lag vs. climate-naive benchmark, ZINB robustness check (Table 3). | `01a`, `03a` | `spi_table3_sarimax_zinb_diagnostics.csv`, `spi_sarimax_zinb_results.rds` |
| `005_spi_rolling_forecast_alerts.R` | Expanding-window rolling-origin forecast evaluation (Table 4) and four-tier alert system validated by bootstrap sensitivity/PPV (Table 5). | `01a`, `03a` | `spi_table4_forecast_evaluation.csv`, `spi_table5_alert_validation.csv`, `spi_rolling_eval_results.rds` || `05b_spi_alert_sensitivity_PI.R` | Supplementary check: does the Table 5 conclusion depend on triggering alerts off the forecast mean vs. its 80%/95% upper prediction limits? Reports the full confusion-matrix metric set (sensitivity, specificity, PPV, NPV, F1, balanced accuracy), bootstrapped. | `01a`, `03a` | `spi_supp_table5_alert_sensitivity_PI_*.csv`, figures |


## Reproducibility notes

- **Random seeds.** Two named seeds are used, defined once in `00_setup.R`:
  `SEED_GLOBAL` (42) for general reproducibility, and `SEED_RESAMPLING`
  (2026) for Monte Carlo / bootstrap procedures specifically (the AR(1)
  wavelet-coherence surrogate test, the bivariate-synchrony surrogate test,
  and the alert-sensitivity bootstrap CIs). The manuscript's current
  Methods text ("a fixed seed (42) ensured reproducibility") should be
  updated to describe both seeds explicitly.
- **R version and session info.** Analyses were run under R 4.5.3. Run
  `sessionInfo()` after `source("run_all.R")` and save the output alongside
  any archived release (e.g. to Zenodo) so package versions are recorded.
- **Namespace conflicts.** `signal::filter()` and `MASS::select()` both mask
  `dplyr::filter()`/`dplyr::select()` when loaded after `tidyverse`.
  `00_setup.R` pins `filter`, `select`, and `lag` to their `dplyr` versions
  explicitly so stage scripts behave the same regardless of package load
  order.

## Known limitations (see manuscript Discussion / Limitations)


## License

Add a license file appropriate for your intended use (e.g. MIT for code;
consider CC-BY for any accompanying processed data you choose to release).
