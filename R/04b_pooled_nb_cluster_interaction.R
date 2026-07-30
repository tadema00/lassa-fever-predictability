# ============================================================================
# 04b_pooled_nb_cluster_interaction.R
# Manuscript Section 2.7 / 3.7: Pooled Cross-Cluster Climate-Sensitivity Test
#
#   log(mu_s,t) = alpha_s + beta1 * SPI_(s, t-l_s) + beta2 * (SPI_(s,t-l_s) x Cluster_s)
#
# All six focal states analysed simultaneously in one pooled Negative Binomial
# regression, with:
#   - alpha_s   : a state fixed effect (intercept per state)
#   - beta1     : the SPI-4 slope in the reference cluster (Northern Savannah)
#   - beta2     : the ADDITIONAL SPI-4 slope for the Southern Forest cluster
#                 (Cluster_s coded 1 = Southern Forest, 0 = Northern Savannah,
#                 so the Southern-cluster slope is beta1 + beta2)
#
# SPI_(s, t-l_s) uses each state's OWN consensus lag l_s from Table 1
# (lag_tables_sensitivity.rds), exactly as in the per-state SARIMAX (stage 04)
# and rolling-forecast (stage 05) analyses - so the pooled model is on the
# same climate-exposure definition as the rest of the manuscript, not a fresh
# re-derivation.
#
# NOTE ON WHAT'S NOT INCLUDED: a main effect for `Cluster` alone is NOT in the
# model. Every state belongs to exactly one cluster, so Cluster is a perfectly
# collinear function of the State fixed effects (alpha_s) - including both
# would make Cluster's main-effect coefficient inestimable (R would silently
# return NA for it). alpha_s already absorbs any baseline-incidence difference
# between clusters; beta2 isolates the DIFFERENTIAL CLIMATE SENSITIVITY the
# manuscript is testing, uncontaminated by that baseline difference.
#
# Inference is reported two ways:
#   1. Model-based SEs/Wald z/p from glm.nb() directly.
#   2. Cluster-robust (sandwich) SEs, clustered by State, via sandwich::vcovCL.
#      Weekly observations within a state are serially correlated, so the
#      naive glm.nb() SEs likely understate uncertainty; the robust version
#      is the more defensible one to lead with in the manuscript, with the
#      naive version shown alongside for transparency.
#   3. A likelihood-ratio test of beta2 = 0 (interaction vs. no-interaction
#      model), which is the direct test of "does climate sensitivity differ
#      by ecological cluster" and is invariant to the SE choice above.
#
# Inputs:  outputs/prepared_data.rds, outputs/lag_tables_sensitivity.rds
# Outputs: outputs/tables/spi_table7_pooled_nb_cluster_interaction.csv
#          outputs/tables/spi_table7b_pooled_nb_lr_test.csv
# ============================================================================

source("00_setup.R")

# `sandwich` isn't in 00_setup.R's central package list (it's only needed by
# this one supplementary script), so it's checked/installed locally here.
if (!requireNamespace("sandwich", quietly = TRUE)) {
  install.packages("sandwich", repos = "https://cloud.r-project.org")
}
library(sandwich)

prepped <- readRDS(file.path(OUT_DIR, "prepared_data.rds"))
weekly  <- prepped$weekly

lag_results <- readRDS(file.path(OUT_DIR, "lag_tables_sensitivity.rds"))
lag_table   <- lag_results$interp$lag_table

SPI_SCALE <- 4
clim_col  <- paste0("SPI_", SPI_SCALE)

# ---------------------------------------------------------------------------
# 1. Build the pooled analysis dataset: each state's SPI-4 lagged by ITS OWN
#    consensus lag (Table 1), stacked into one long panel with State + Cluster.
# ---------------------------------------------------------------------------
pooled_data <- map_dfr(FOCAL_STATES, function(st) {
  df_state <- weekly %>% filter(State == st) %>% arrange(week_idx)
  
  cons_lag <- lag_table %>% filter(State == st) %>% pull(Consensus_lag_wk)
  if (length(cons_lag) != 1) {
    stop(sprintf("Expected exactly one Consensus_lag_wk for state '%s', got %d.", st, length(cons_lag)))
  }
  
  df_state %>%
    mutate(SPI_lag = dplyr::lag(.data[[clim_col]], cons_lag)) %>%
    select(State, week_idx, Year, Epi_Week, Cases, SPI_lag)
}) %>%
  filter(!is.na(SPI_lag)) %>%
  mutate(
    State        = factor(State),
    Cluster      = factor(cluster_of(as.character(State)),
                          levels = c("Northern Savannah", "Southern Forest")),
    ClusterSouth = as.integer(Cluster == "Southern Forest")
  )

cat(sprintf("Pooled dataset: %d state-weeks across %d states (%d Southern, %d Northern).\n",
            nrow(pooled_data), n_distinct(pooled_data$State),
            sum(pooled_data$ClusterSouth == 1), sum(pooled_data$ClusterSouth == 0)))

# ---------------------------------------------------------------------------
# 2. Fit the full model (with interaction) and the reduced model (without),
#    both as pooled Negative Binomial regressions with state fixed effects.
# ---------------------------------------------------------------------------
full_fit <- MASS::glm.nb(
  Cases ~ State + SPI_lag + SPI_lag:ClusterSouth,
  data = pooled_data
)

reduced_fit <- MASS::glm.nb(
  Cases ~ State + SPI_lag,
  data = pooled_data
)

# ---------------------------------------------------------------------------
# 3. Model-based (naive) coefficient table for the two climate terms.
# ---------------------------------------------------------------------------
coef_naive <- summary(full_fit)$coefficients
beta1_naive <- coef_naive["SPI_lag", ]
beta2_naive <- coef_naive["SPI_lag:ClusterSouth", ]

# ---------------------------------------------------------------------------
# 4. Cluster-robust (sandwich) coefficient table, clustered by State - the
#    recommended version to report given within-state serial correlation.
# ---------------------------------------------------------------------------
vcov_cl <- sandwich::vcovCL(full_fit, cluster = pooled_data$State)
coef_robust <- lmtest::coeftest(full_fit, vcov. = vcov_cl)
beta1_robust <- coef_robust["SPI_lag", ]
beta2_robust <- coef_robust["SPI_lag:ClusterSouth", ]

# ---------------------------------------------------------------------------
# 5. Likelihood-ratio test of H0: beta2 = 0 (full vs. reduced model) - the
#    direct test of differential climate sensitivity between clusters.
# ---------------------------------------------------------------------------
lr_test <- lmtest::lrtest(reduced_fit, full_fit)
lr_stat <- lr_test$Chisq[2]
lr_df   <- lr_test$`#Df`[2] - lr_test$`#Df`[1]
lr_p    <- lr_test$`Pr(>Chisq)`[2]

# ---------------------------------------------------------------------------
# 6. Assemble the Section 3.7 results table.
# ---------------------------------------------------------------------------
table7 <- tibble(
  Term = c("beta1: SPI-4 (Northern Savannah, reference)",
           "beta2: SPI-4 x Cluster (Southern Forest, additional slope)"),
  Estimate         = c(beta1_naive["Estimate"], beta2_naive["Estimate"]),
  SE_model         = c(beta1_naive["Std. Error"], beta2_naive["Std. Error"]),
  z_model          = c(beta1_naive["z value"], beta2_naive["z value"]),
  p_model          = c(beta1_naive["Pr(>|z|)"], beta2_naive["Pr(>|z|)"]),
  SE_robust        = c(beta1_robust["Std. Error"], beta2_robust["Std. Error"]),
  z_robust         = c(beta1_robust["z value"], beta2_robust["z value"]),
  p_robust         = c(beta1_robust["Pr(>|z|)"], beta2_robust["Pr(>|z|)"])
) %>%
  mutate(across(where(is.numeric), ~ round(.x, 4)))

# Southern-cluster total slope (beta1 + beta2) with its own robust SE, via
# the delta method (Var(b1+b2) = Var(b1) + Var(b2) + 2*Cov(b1,b2)).
south_slope_est <- beta1_naive["Estimate"] + beta2_naive["Estimate"]
south_slope_var_robust <- vcov_cl["SPI_lag", "SPI_lag"] +
  vcov_cl["SPI_lag:ClusterSouth", "SPI_lag:ClusterSouth"] +
  2 * vcov_cl["SPI_lag", "SPI_lag:ClusterSouth"]
south_slope_se_robust <- sqrt(south_slope_var_robust)
south_slope_z_robust <- south_slope_est / south_slope_se_robust
south_slope_p_robust <- 2 * pnorm(-abs(south_slope_z_robust))

table7_south <- tibble(
  Term = "Southern Forest total SPI-4 slope (beta1 + beta2)",
  Estimate = round(south_slope_est, 4),
  SE_model = NA_real_, z_model = NA_real_, p_model = NA_real_,
  SE_robust = round(south_slope_se_robust, 4),
  z_robust  = round(south_slope_z_robust, 4),
  p_robust  = round(south_slope_p_robust, 4)
)

table7_full <- bind_rows(table7, table7_south)

model_summary <- tibble(
  Model = c("Reduced (no interaction)", "Full (with interaction)"),
  AIC   = round(c(AIC(reduced_fit), AIC(full_fit)), 1),
  theta = round(c(reduced_fit$theta, full_fit$theta), 4),
  SE_theta = round(c(reduced_fit$SE.theta, full_fit$SE.theta), 4)
)

lr_summary <- tibble(
  Test = "LR test, H0: beta2 = 0 (interaction term)",
  Chisq = round(lr_stat, 3),
  df = lr_df,
  p_value = signif(lr_p, 4)
)

print(table7_full)
cat("\n")
print(model_summary)
cat("\n")
print(lr_summary)

write_csv(table7_full, file.path(TAB_DIR, "spi_table7_pooled_nb_cluster_interaction.csv"))
write_csv(bind_rows(
  model_summary %>% mutate(Chisq = NA_real_, df = NA_integer_, p_value = NA_real_, Test = Model) %>% select(Test, AIC, theta, SE_theta, Chisq, df, p_value),
  lr_summary %>% mutate(AIC = NA_real_, theta = NA_real_, SE_theta = NA_real_) %>% select(Test, AIC, theta, SE_theta, Chisq, df, p_value)
),
file.path(TAB_DIR, "spi_table7b_pooled_nb_lr_test.csv")
)

cat("\nSection 2.7/3.7 pooled cross-cluster NB interaction model complete (Table 7).\n")
cat("Report the cluster-robust (SE_robust/p_robust) columns as primary; model-based\n",
    "columns (SE_model/p_model) are shown for transparency but likely understate\n",
    "uncertainty given within-state serial correlation across weeks.\n", sep = "")



