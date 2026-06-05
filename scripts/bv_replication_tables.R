## Replicate the all-cause-mortality components of Borgschulte & Vogler (2020)
## Tables 1 and 2 with public data, and compare against the published values
## (docs/Borgschulte and Vogler (2020).pdf; text in docs/bv2020_text.txt).
##
## Procedure follows BV Sections 3.1-3.2:
##   (i)   collapse to one row per county, averaging 2009-2013 covariates;
##   (ii)  double lasso (square-root lasso per BV footnote 8, via RPtests;
##         falls back to hdm post-lasso when RPtests is absent) selecting
##         predictors of (a) pre-expansion mortality and (b) expansion;
##   (iii) logit propensity score on the union of selected predictors;
##   (iv)  trim p-hat outside [0.038, 0.971] (BV's published thresholds);
##   (v)   weighted DID with county + year FE, state-clustered SEs, weights =
##         county population 20-64 x [T + (1-T) p/(1-p)].
##
## Known remaining deviations from BV (see docs/bv_replication_todo.md):
##   - Public CDC WONDER mortality (suppression) vs restricted NVSS.
##   - Outcome-lasso DV uses avg 2009-2013 mortality until the panel is
##     rebuilt with 2005-2009 data (BV hold out 2010-2013 and use 2005-2009).
##   - Missing PS candidates: state health & welfare expenditures 2005-2013.
##   - log_20_64 / pct_55_64 are corrected here from existing columns
##     (log(WONDER population); ACS pct_5564) until analysis_data.csv is
##     regenerated with the fixed county_dataset_creation.R.
##
## Output:
##   paper/data/bv_table1_comparison.csv  (Table 1, ours vs published)
##   paper/data/bv_table2_comparison.csv  (Table 2 Panel A cols 1/2/10/11)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
})

USE_SQRT_LASSO <- requireNamespace("RPtests", quietly = TRUE)
TRIM_LO <- 0.038
TRIM_HI <- 0.971

# ---------------------------------------------------------------------------
# 1. Load data; corrected control variables
# ---------------------------------------------------------------------------

## Panel covariates are correct as of the SEER rebuild (fixed log_20_64 and
## pct_55_64 definitions, time-varying so they survive the county FE).
panel <- read_csv("paper/data/analysis_data.csv", show_col_types = FALSE)

bv_cov <- read_csv("paper/data/bv_covariates.csv", show_col_types = FALSE)

## BV's six panel-lasso-selected DID controls (their table notes)
panel_controls <- c("pct_white", "pct_55_64", "log_20_64",
                    "log_35_44", "log_f_20_64", "unemp")

## NOTE: avg_mortality_pre (2009-13 avg from the OLD panel) is deliberately
## absent -- it isn't used (the outcome lasso targets the 2005-09 yearly
## average) and joining it would drop_na() away every county recovered by the
## suppression imputation, since the old panel never contained them.
bv_cov_vars <- c("pct_male", "pct_black", "pct_hispanic",
                 "pct_2024", "pct_2534", "pct_3544", "pct_4554", "pct_5564",
                 "poverty_rate", "log_median_income", "log_pop_density",
                 "uninsured_rate", "dem_governor_2010",
                 "obama_share_2008", "obama_share_2012")
bv_cov_vars <- intersect(bv_cov_vars, names(bv_cov))  # tolerate missing pulls

# ---------------------------------------------------------------------------
# 2. County-level baseline frame (2009-2013 means)
# ---------------------------------------------------------------------------

## Year-by-year pre-expansion mortality 2005-2009 (PS candidates, per BV
## Appendix Table A.2; 2010-2013 outcomes stay held out). 2005-2008 come from
## the rebuilt county_mortality_pre.csv, 2009 from the panel itself.
mort_yearly <- bind_rows(
  read_csv("paper/data/county_mortality_pre.csv", show_col_types = FALSE) %>%
    select(fips, year, crude_rate),
  panel %>% filter(year == 2009) %>% select(fips, year, crude_rate)
) %>%
  filter(year >= 2005, year <= 2009) %>%
  tidyr::pivot_wider(names_from = year, values_from = crude_rate,
                     names_prefix = "mort_")

MORT_YEARLY_VARS <- paste0("mort_", 2005:2009)

baseline <- panel %>%
  filter(year >= 2009, year <= 2013) %>%
  group_by(fips, state_fips = as.numeric(state_fips), expansion) %>%
  summarise(across(c(pct_white, pct_55_64, log_20_64, log_35_44,
                     log_f_20_64, unemp),
                   ~ mean(.x, na.rm = TRUE)),
            crude_rate_pre = mean(crude_rate, na.rm = TRUE),
            pop_2064_pre   = mean(population, na.rm = TRUE),
            .groups = "drop") %>%
  left_join(bv_cov %>% select(fips, any_of(bv_cov_vars)), by = "fips") %>%
  left_join(mort_yearly, by = "fips") %>%
  drop_na()

message("Complete-case counties: ", nrow(baseline),
        " (", sum(baseline$expansion), " expansion)")

# ---------------------------------------------------------------------------
# 3. Double lasso (square-root lasso per BV fn. 8) and propensity score
# ---------------------------------------------------------------------------

## Candidate set now mirrors BV fn. 8: demographics/economics/politics plus
## 2005-2009 all-cause mortality year-by-year. The 2009-2013 average
## (avg_mortality_pre) is excluded -- BV hold out 2010-2013 outcomes.
control_vars <- unique(c(panel_controls,
                         setdiff(bv_cov_vars, "avg_mortality_pre"),
                         MORT_YEARLY_VARS))
X_mat <- scale(as.matrix(baseline[, control_vars]))
X_mat[!is.finite(X_mat)] <- 0
D_vec <- baseline$expansion
Y_vec <- rowMeans(baseline[, MORT_YEARLY_VARS])  # 2005-2009 avg, holdout intact

out_keep <- control_vars  # mortality target is no longer itself a candidate
out_idx  <- which(control_vars %in% out_keep)

if (USE_SQRT_LASSO) {
  message("Variable selection: square-root lasso (RPtests), as in BV fn. 8")
  fit_treat   <- RPtests::sqrt_lasso(X_mat, as.numeric(D_vec))
  fit_outcome <- RPtests::sqrt_lasso(X_mat[, out_idx, drop = FALSE],
                                     as.numeric(scale(Y_vec)))
  sel_treat   <- which(abs(fit_treat)   > 1e-8)
  sel_outcome <- out_idx[which(abs(fit_outcome) > 1e-8)]
} else {
  message("Variable selection: hdm post-lasso fallback (install RPtests for sqrt-lasso)")
  fit_treat   <- hdm::rlassologit(X_mat, D_vec, post = TRUE)
  fit_outcome <- hdm::rlasso(X_mat[, out_idx, drop = FALSE], Y_vec, post = TRUE)
  sel_treat   <- which(fit_treat$index)
  sel_outcome <- out_idx[which(fit_outcome$index)]
}

sel_union <- sort(unique(c(sel_treat, sel_outcome)))
message("Selected ", length(sel_union), " of ", length(control_vars),
        " candidates: ", paste(control_vars[sel_union], collapse = ", "))

ps_df <- data.frame(D = D_vec, X_mat[, sel_union, drop = FALSE])
baseline$ps <- predict(glm(D ~ ., data = ps_df, family = binomial()),
                       type = "response")

trimmed <- baseline %>% filter(ps >= TRIM_LO, ps <= TRIM_HI)
message("Trim removes ", nrow(baseline) - nrow(trimmed), " counties (",
        sum(baseline$expansion) - sum(trimmed$expansion), " expansion, ",
        sum(1 - baseline$expansion) - sum(1 - trimmed$expansion),
        " non-expansion); final sample ", nrow(trimmed), " counties (",
        sum(trimmed$expansion), " expansion / ", sum(1 - trimmed$expansion),
        " non-expansion). BV: 563 trimmed (78/485); final 2260 (1174/1086).")

## Diagnostic: trimming big counties materially changes the (weighted)
## estimand, so report what the trim removes.
trimmed_out <- baseline %>% filter(ps < TRIM_LO | ps > TRIM_HI)
cat("\n-- 10 largest trimmed counties (by 20-64 pop) --\n")
print(trimmed_out %>% arrange(desc(pop_2064_pre)) %>%
        transmute(fips, expansion, ps = round(ps, 3),
                  pop_2064 = round(pop_2064_pre)) %>% head(10),
      n = 10)
cat(sprintf("share of expansion 20-64 pop trimmed: %.1f%% | control: %.1f%%\n",
            100 * sum(trimmed_out$pop_2064_pre[trimmed_out$expansion == 1]) /
              sum(baseline$pop_2064_pre[baseline$expansion == 1]),
            100 * sum(trimmed_out$pop_2064_pre[trimmed_out$expansion == 0]) /
              sum(baseline$pop_2064_pre[baseline$expansion == 0])))

# ---------------------------------------------------------------------------
# 4. Table 1: population-weighted baseline means (SD), ours vs published
# ---------------------------------------------------------------------------

wmean <- function(x, w) sum(x * w, na.rm = TRUE) / sum(w[!is.na(x)])
wsd   <- function(x, w) {
  m <- wmean(x, w)
  sqrt(sum(w * (x - m)^2, na.rm = TRUE) / sum(w[!is.na(x)]))
}

table1_rows <- tribble(
  ~variable,                              ~col,                  ~bv_exp_mean, ~bv_exp_sd, ~bv_ctrl_mean, ~bv_ctrl_sd,
  "All cause mortality (per 100,000)",    "crude_rate_pre",      315.17, 94.69, 359.06, 117.06,
  "% Population ages 20-24",              "pct_2024",            0.12, 0.03, 0.12, 0.04,
  "% Population ages 25-34",              "pct_2534",            0.23, 0.03, 0.22, 0.03,
  "% Population ages 35-44",              "pct_3544",            0.22, 0.02, 0.22, 0.02,
  "% Population ages 45-54",              "pct_4554",            0.24, 0.02, 0.24, 0.02,
  "% Population ages 55-64",              "pct_5564",            0.20, 0.03, 0.20, 0.04,
  "% Male",                               "pct_male",            0.50, 0.01, 0.49, 0.02,
  "% Hispanic",                           "pct_hispanic",        0.19, 0.16, 0.16, 0.19,
  "% White",                              "pct_white",           0.78, 0.14, 0.79, 0.13,
  "% Black",                              "pct_black",           0.11, 0.11, 0.15, 0.12,
  "Unemployed rate",                      "unemp_frac",          0.09, 0.02, 0.08, 0.02,
  "Poverty rate",                         "poverty_rate",        0.15, 0.05, 0.16, 0.05,
  "Real median income ($10,000)",         "med_income_10k",      5.95, 1.46, 5.33, 1.49,
  "Uninsured rate",                       "uninsured_rate",      0.19, 0.06, 0.24, 0.07,
  "Obama 2008 vote share",                "obama_share_2008",    0.46, 0.12, 0.38, 0.13,
  "Obama 2012 vote share",                "obama_share_2012",    0.43, 0.13, 0.35, 0.14,
  "Democratic governor",                  "dem_governor_2010",   0.76, 0.43, 0.49, 0.50
)

tab1_frame <- trimmed %>%
  left_join(bv_cov %>% select(fips, med_income), by = "fips") %>%
  mutate(unemp_frac = unemp / 100,            # panel unemp is 0-100
         med_income_10k = med_income / 10000)

table1 <- table1_rows %>%
  rowwise() %>%
  mutate(
    our_exp_mean  = if (col %in% names(tab1_frame))
      wmean(tab1_frame[[col]][tab1_frame$expansion == 1],
            tab1_frame$pop_2064_pre[tab1_frame$expansion == 1]) else NA_real_,
    our_exp_sd    = if (col %in% names(tab1_frame))
      wsd(tab1_frame[[col]][tab1_frame$expansion == 1],
          tab1_frame$pop_2064_pre[tab1_frame$expansion == 1]) else NA_real_,
    our_ctrl_mean = if (col %in% names(tab1_frame))
      wmean(tab1_frame[[col]][tab1_frame$expansion == 0],
            tab1_frame$pop_2064_pre[tab1_frame$expansion == 0]) else NA_real_,
    our_ctrl_sd   = if (col %in% names(tab1_frame))
      wsd(tab1_frame[[col]][tab1_frame$expansion == 0],
          tab1_frame$pop_2064_pre[tab1_frame$expansion == 0]) else NA_real_
  ) %>%
  ungroup() %>%
  mutate(gap_exp = our_exp_mean - bv_exp_mean,
         gap_ctrl = our_ctrl_mean - bv_ctrl_mean)

write_csv(table1, "paper/data/bv_table1_comparison.csv")

cat("\n==== Table 1 comparison (population-weighted means, 2009-2013) ====\n")
print(table1 %>%
        transmute(variable,
                  ours_exp = round(our_exp_mean, 2), bv_exp = bv_exp_mean,
                  ours_ctrl = round(our_ctrl_mean, 2), bv_ctrl = bv_ctrl_mean),
      n = 30)

## Diagnostic: BV's expansion-group Dem-governor mean (0.76) equals our
## UNWEIGHTED share of expansion counties in Dem-governor states, suggesting
## the political rows of their Table 1 may be unweighted despite the caption.
cat("\n-- Political variables, unweighted county means (BV exp/ctrl: ",
    "Obama08 0.46/0.38, Obama12 0.43/0.35, DemGov 0.76/0.49) --\n", sep = "")
for (v in c("obama_share_2008", "obama_share_2012", "dem_governor_2010")) {
  cat(sprintf("%-18s exp %.2f  ctrl %.2f\n", v,
              mean(trimmed[[v]][trimmed$expansion == 1], na.rm = TRUE),
              mean(trimmed[[v]][trimmed$expansion == 0], na.rm = TRUE)))
}

# ---------------------------------------------------------------------------
# 5. Table 2 Panel A: weighted DID, ours vs published
# ---------------------------------------------------------------------------

stopifnot(requireNamespace("fixest", quietly = TRUE))

panel_ipw <- panel %>%
  inner_join(trimmed %>%
               transmute(fips, ps,
                         uninsured_rate_base = uninsured_rate,
                         ipw = ifelse(expansion == 1, 1, ps / (1 - ps))),
             by = "fips") %>%
  mutate(final_weight = population * ipw,
         treated_post = expansion * post)

baseline_mort_exp <- wmean(trimmed$crude_rate_pre[trimmed$expansion == 1],
                           trimmed$pop_2064_pre[trimmed$expansion == 1])

run_did <- function(dat, controls = NULL, label = "") {
  rhs <- if (is.null(controls)) "treated_post" else
    paste(c("treated_post", controls), collapse = " + ")
  mod <- fixest::feols(
    stats::as.formula(paste0("crude_rate ~ ", rhs, " | fips + year")),
    data = dat, weights = ~ final_weight, cluster = ~ state_fips)
  ct <- fixest::coeftable(mod)["treated_post", ]
  tibble(model = label, estimate = ct[1], se = ct[2],
         pct_effect = 100 * ct[1] / baseline_mort_exp,
         n_obs = stats::nobs(mod),
         n_clusters = length(unique(dat$state_fips)))
}

uninsured_median <- stats::median(trimmed$uninsured_rate, na.rm = TRUE)

table2 <- bind_rows(
  run_did(panel_ipw, NULL, "(1) Base"),
  run_did(panel_ipw, panel_controls, "(2) Controls"),
  run_did(panel_ipw %>% filter(uninsured_rate_base >  uninsured_median),
          panel_controls, "(10) High uninsured"),
  run_did(panel_ipw %>% filter(uninsured_rate_base <= uninsured_median),
          panel_controls, "(11) Low uninsured")
) %>%
  mutate(bv_estimate = c(-14.83, -11.36, -12.40, -3.96),
         bv_se       = c(6.12, 3.59, 4.22, 4.42),
         bv_pct      = c(-4.71, -3.60, -3.81, -1.30),
         bv_n_obs    = c(20340, 20340, 18882, 11214))

write_csv(table2, "paper/data/bv_table2_comparison.csv")

cat("\n==== Table 2 Panel A comparison (all-cause mortality) ====\n")
print(table2 %>%
        transmute(model,
                  ours = sprintf("%.2f (%.2f)", estimate, se),
                  bv   = sprintf("%.2f (%.2f)", bv_estimate, bv_se),
                  ours_pct = round(pct_effect, 2), bv_pct,
                  n_obs, bv_n_obs))

cat("\nNotes:\n",
    "- BV columns (10)/(11) split at the median baseline uninsured rate for\n",
    "  ages 19-64; we use SAHIE 18-64 and an unweighted county median.\n",
    "- Gender/age columns (3)-(9) require by-sex/by-age WONDER exports.\n",
    "- Rebuild analysis_data.csv (suppression imputation + covariate fixes)\n",
    "  before quoting these numbers; see docs/bv_replication_todo.md.\n")
