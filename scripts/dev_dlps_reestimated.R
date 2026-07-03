## DEV REFERENCE (SRV-2, "reestimated" variant) -- to be ported into a
## guarded `dlps_placebo` block in analyze_application.R; see
## docs/server_build_list_2026-07.md for the output schema. Delete after
## porting. Reproduces the re-estimated-DLPS placebo numbers in
## docs/extend_preperiod_2005.md (2026-07-03): the full pipeline (double-
## lasso selection, logit PS, trim) is re-fit at each fake onset on
## pre-onset data only, with ALL pre-onset mortality years as predictors
## (deliberately more generous than BV's 4-year-holdout convention).
## Residual caveat: the ACS/SAHIE/political covariates keep their fixed
## 2009-2013 vintage (slow-moving; unavoidable).
suppressPackageStartupMessages({library(dplyr); library(readr); library(tidyr)})
set.seed(20240101)

ext    <- read_csv("paper/data/analysis_data_2005.csv", show_col_types = FALSE)
bv_cov <- read_csv("paper/data/bv_covariates.csv", show_col_types = FALSE)
bv_cov_vars <- c("pct_male", "pct_black", "pct_hispanic", "pct_2024", "pct_2534", "pct_3544",
                 "pct_4554", "pct_5564", "poverty_rate", "log_median_income", "log_pop_density",
                 "uninsured_rate", "dem_governor_2010", "obama_share_2008", "obama_share_2012")
panel_controls <- c("pct_white", "pct_55_64", "log_20_64", "log_35_44", "log_f_20_64", "unemp")
mort_wide <- ext %>% select(fips, year, crude_rate) %>%
  pivot_wider(names_from = year, values_from = crude_rate, names_prefix = "mort_")

dlps_at_onset <- function(py, w_start) {
  mort_vars_py <- paste0("mort_", w_start:(py - 1))   # ALL pre-onset years (best case)
  base <- ext %>% filter(year >= w_start, year < py) %>%
    group_by(fips, state_fips = as.numeric(state_fips), expansion) %>%
    summarise(across(all_of(panel_controls), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
    left_join(bv_cov %>% select(fips, all_of(bv_cov_vars)), by = "fips") %>%
    left_join(mort_wide %>% select(fips, all_of(mort_vars_py)), by = "fips") %>%
    drop_na()
  cv <- c(panel_controls, bv_cov_vars, mort_vars_py)
  X <- scale(as.matrix(base[, cv])); X[!is.finite(X)] <- 0
  D <- base$expansion; Yv <- rowMeans(base[, mort_vars_py])
  st <- which(abs(RPtests::sqrt_lasso(X, as.numeric(D))) > 1e-8)
  so <- which(abs(RPtests::sqrt_lasso(X, as.numeric(scale(Yv)))) > 1e-8)
  su <- sort(unique(c(st, so)))
  base$ps <- predict(glm(D ~ ., data = data.frame(D = D, X[, su, drop = FALSE]),
                         family = binomial()), type = "response")
  trim <- base %>% filter(ps >= 0.038, ps <= 0.971) %>%
    mutate(ipw = if_else(expansion == 1, 1, ps / (1 - ps)))
  dat <- ext %>% filter(year >= w_start, year <= 2013) %>%
    inner_join(trim %>% select(fips, ipw), by = "fips") %>%
    mutate(final_weight = population * ipw, state_fips_n = as.numeric(state_fips),
           treated_post = expansion * as.integer(year >= py))
  fe <- function(controls = NULL) {
    rhs <- paste(c("treated_post", controls), collapse = " + ")
    m <- fixest::feols(stats::as.formula(paste0("crude_rate ~ ", rhs, " | fips + year")),
                       data = dat, weights = ~final_weight, cluster = ~state_fips_n)
    ct <- fixest::coeftable(m)["treated_post", ]
    sprintf("%7.2f (%.2f)", ct[1], ct[2])
  }
  cat(sprintf("  onset %d (PS re-fit on %d-%d, %d selected, trim %d): base %s   controls %s\n",
              py, w_start, py - 1, length(su), nrow(base) - nrow(trim),
              fe(), fe(panel_controls)))
}
cat("== DLPS placebo, FULL pipeline re-estimated at each onset ==\n")
cat("-- window 2005-2013 (mortality predictors = all pre-onset years) --\n")
for (py in 2009:2012) dlps_at_onset(py, 2005)
cat("-- window 2009-2013 --\n")
for (py in 2011:2012) dlps_at_onset(py, 2009)
