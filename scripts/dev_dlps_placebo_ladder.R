## DEV REFERENCE (SRV-2 "fixed" variant + SRV-3) -- to be ported into guarded
## `dlps_placebo` and `dlps_ladder` blocks in analyze_application.R; see
## docs/server_build_list_2026-07.md for output schemas. Delete after porting.
## Reproduces the 2026-07-03 numbers in docs/extend_preperiod_2005.md.
##
## (1) DLPS in-time placebo: hold BV's propensity-score weights fixed (they are
##     functions of pre-2014 data only), restrict the outcome regression to
##     pre-2014 windows, and estimate at simulated onsets.
## (2) DID -> DLPS bridge ladder: from pop-weighted TWFE DID to BV's preferred
##     spec, one ingredient per rung, 2009-2017.
## PS construction copied from scripts/analyze_application.R (DLPS block).
suppressPackageStartupMessages({library(dplyr); library(readr); library(tidyr)})
set.seed(20240101)

adf    <- read_csv("paper/data/analysis_data.csv", show_col_types = FALSE)
ext    <- read_csv("paper/data/analysis_data_2005.csv", show_col_types = FALSE)
bv_cov <- read_csv("paper/data/bv_covariates.csv", show_col_types = FALSE)
mort_yearly <- bind_rows(
  read_csv("paper/data/county_mortality_pre.csv", show_col_types = FALSE) %>%
    select(fips, year, crude_rate),
  adf %>% filter(year == 2009) %>% select(fips, year, crude_rate)
) %>% filter(year >= 2005, year <= 2009) %>%
  pivot_wider(names_from = year, values_from = crude_rate, names_prefix = "mort_")
mort_vars      <- paste0("mort_", 2005:2009)
panel_controls <- c("pct_white", "pct_55_64", "log_20_64", "log_35_44", "log_f_20_64", "unemp")
bv_cov_vars    <- c("pct_male", "pct_black", "pct_hispanic", "pct_2024", "pct_2534", "pct_3544",
                    "pct_4554", "pct_5564", "poverty_rate", "log_median_income", "log_pop_density",
                    "uninsured_rate", "dem_governor_2010", "obama_share_2008", "obama_share_2012")
dlps_base <- adf %>% filter(year <= 2013) %>%
  group_by(fips, state_fips = as.numeric(state_fips), expansion) %>%
  summarise(across(all_of(panel_controls), ~ mean(.x, na.rm = TRUE)),
            crude_rate_pre = mean(crude_rate), pop_2064_pre = mean(population), .groups = "drop") %>%
  left_join(bv_cov %>% select(fips, all_of(bv_cov_vars)), by = "fips") %>%
  left_join(mort_yearly, by = "fips") %>% drop_na()
control_vars <- c(panel_controls, bv_cov_vars, mort_vars)
X_mat <- scale(as.matrix(dlps_base[, control_vars])); X_mat[!is.finite(X_mat)] <- 0
D_vec <- dlps_base$expansion; Y_vec <- rowMeans(dlps_base[, mort_vars])
if (requireNamespace("RPtests", quietly = TRUE)) {
  sel_treat <- which(abs(RPtests::sqrt_lasso(X_mat, as.numeric(D_vec))) > 1e-8)
  sel_out   <- which(abs(RPtests::sqrt_lasso(X_mat, as.numeric(scale(Y_vec)))) > 1e-8)
} else {
  sel_treat <- which(hdm::rlassologit(X_mat, D_vec, post = TRUE)$index)
  sel_out   <- which(hdm::rlasso(X_mat, Y_vec, post = TRUE)$index)
}
sel_union <- sort(unique(c(sel_treat, sel_out)))
ps_df <- data.frame(D = D_vec, X_mat[, sel_union, drop = FALSE])
dlps_base$ps <- predict(glm(D ~ ., data = ps_df, family = binomial()), type = "response")
dlps_base <- dlps_base %>%
  mutate(ipw = if_else(expansion == 1, 1, ps / (1 - ps)))
dlps_trim <- dlps_base %>% filter(ps >= 0.038, ps <= 0.971)
cat(sprintf("PS rebuilt: %d complete-case, %d trimmed away, %d selected covars\n",
            nrow(dlps_base), nrow(dlps_base) - nrow(dlps_trim), length(sel_union)))

twfe <- function(dat, wcol, controls = NULL) {
  rhs <- paste(c("treated_post", controls), collapse = " + ")
  m <- fixest::feols(stats::as.formula(paste0("crude_rate ~ ", rhs, " | fips + year")),
                     data = dat, weights = stats::as.formula(paste0("~", wcol)),
                     cluster = ~ state_fips_n)
  ct <- fixest::coeftable(m)["treated_post", ]
  sprintf("%7.2f (%.2f)", ct[1], ct[2])
}

## ---- (1) DLPS in-time placebos -------------------------------------------
prep <- function(panel, trim_tbl) {
  panel %>%
    inner_join(trim_tbl %>% select(fips, ipw), by = "fips") %>%
    mutate(final_weight = population * ipw, state_fips_n = as.numeric(state_fips))
}
cat("\n== DLPS in-time placebos (BV weights held fixed; est (SE)) ==\n")
cat("-- window 2009-2013 (canonical) --\n")
p_can <- prep(adf %>% filter(year <= 2013), dlps_trim)
for (py in 2011:2012) {
  d <- p_can %>% mutate(treated_post = expansion * as.integer(year >= py))
  cat(sprintf("  onset %d: base %s   controls %s\n", py,
              twfe(d, "final_weight"), twfe(d, "final_weight", panel_controls)))
}
cat("-- window 2005-2013 (extended) --\n")
p_ext <- prep(ext %>% filter(year <= 2013), dlps_trim)
for (py in 2009:2012) {
  d <- p_ext %>% mutate(treated_post = expansion * as.integer(year >= py))
  cat(sprintf("  onset %d: base %s   controls %s\n", py,
              twfe(d, "final_weight"), twfe(d, "final_weight", panel_controls)))
}

## ---- (2) DID -> DLPS bridge ladder (2009-2017) ---------------------------
cat("\n== Bridge ladder, pop-weighted TWFE DID -> BV DLPS (2009-2017) ==\n")
full <- adf %>% mutate(treated_post = expansion * as.integer(year >= 2014),
                       state_fips_n = as.numeric(state_fips), final_weight = population)
cc    <- full %>% semi_join(dlps_base, by = "fips")
trimd <- full %>% semi_join(dlps_trim, by = "fips")
ipw_cc   <- full  %>% inner_join(dlps_base %>% select(fips, ipw), by = "fips") %>%
  mutate(final_weight = population * ipw)
ipw_trim <- full  %>% inner_join(dlps_trim %>% select(fips, ipw), by = "fips") %>%
  mutate(final_weight = population * ipw)
cat(sprintf("  L1 pop TWFE, all %d counties           : %s\n",
            n_distinct(full$fips),  twfe(full,  "final_weight")))
cat(sprintf("  L2 pop TWFE, complete-case (%d)        : %s\n",
            n_distinct(cc$fips),    twfe(cc,    "final_weight")))
cat(sprintf("  L3 pop TWFE, PS-trimmed sample (%d)    : %s\n",
            n_distinct(trimd$fips), twfe(trimd, "final_weight")))
cat(sprintf("  L4 pop x IPW, complete-case untrimmed   : %s\n", twfe(ipw_cc,   "final_weight")))
cat(sprintf("  L5 pop x IPW, trimmed  [= BV base]      : %s\n", twfe(ipw_trim, "final_weight")))
cat(sprintf("  L6 + panel controls    [= BV preferred] : %s\n",
            twfe(ipw_trim, "final_weight", panel_controls)))
cat("  [reference: package did_w -18.39; committed DLPS base -9.56 (5.05), controls -6.62 (4.14)]\n")
