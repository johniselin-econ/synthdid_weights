# Plan: Thorough DLPS Replication of Borgschulte & Vogler (2020)

## Goal

Replicate the double-lasso propensity score (DLPS) weighted DiD estimate of Borgschulte & Vogler (2020) as closely as possible using publicly available data. Their headline finding: ACA Medicaid expansion reduced mortality by **11.36 deaths per 100,000** among working-age adults (ages 20-64). This serves as an independent benchmark for our population-weighted SDID estimate (~13.5).

## Current state

Our current DLPS implementation (Supplement Appendix F) uses:
- 6 covariates: pct_white, pct_55_64, log_20_64, log_35_44, log_f_20_64, unemp
- `glmnet::cv.glmnet` (single lasso, not double-lasso)
- Arbitrary trimming at 4th/96th percentiles
- Result: estimate diverges substantially from B&V's -11.36

## What B&V use that we're missing

### Covariates (from their Table 1 / Section III)

| Variable | Source | Status |
|----------|--------|--------|
| Population pct in 5 age groups (20-24, 25-34, 35-44, 45-54, 55-64) | Census ACS | **Partial** (we have pct_55_64, pop_20_64) |
| Pct male | Census ACS | **Missing** |
| Pct white | Census ACS | **Have** |
| Pct Black | Census ACS | **Missing** |
| Pct Hispanic | Census ACS | **Missing** |
| Unemployment rate | BLS LAUS | **Have** |
| Poverty rate | Census SAIPE | **Missing** |
| Log real median household income | Census SAIPE | **Missing** |
| Log population | Census | **Have** (log_20_64) |
| Population density | Census | **Missing** (computable from pop + land area) |
| Obama vote share 2008 | MIT Election Lab | **Missing** |
| Obama vote share 2012 | MIT Election Lab | **Missing** |
| Democratic governor 2010 | Manual coding | **Missing** |
| Pre-expansion uninsured rate (non-elderly adults) | Census SAHIE | **Missing** |
| Log avg state health/welfare expenditures 2005-2013 | Census of Governments | **Missing** |
| Pre-expansion mortality rates 2005-2009 (held out) | CDC WONDER | **Computable** from existing data |
| Avg all-cause/amenable/non-amenable mortality 2005-2013 | CDC WONDER | **Partial** (we have all-cause crude_rate) |

### Key missing categories:
1. **Demographics**: pct_male, pct_black, pct_hispanic, finer age groups
2. **Economic**: poverty rate, median income, population density
3. **Political**: Obama vote shares, Democratic governor
4. **Health access**: uninsured rate, state health expenditures

## Data pull plan

### Source 1: Census ACS 5-year (demographics + economics)
- **API**: `tidycensus::get_acs()` with year=2013 (2009-2013 5-year estimates)
- **Geography**: county level
- **Variables needed**:
  - B01001: Sex by age (for age group shares, pct male)
  - B03002: Hispanic or Latino origin by race (pct Black, pct Hispanic)
  - B17001: Poverty status (poverty rate)
  - B19013: Median household income
  - B01003: Total population (for density with land area)
- **Effort**: ~30 lines of R code using tidycensus

### Source 2: Census SAIPE (poverty — alternative to ACS)
- **API**: `censusapi` package or direct download
- **URL**: https://www.census.gov/data/datasets/time-series/demo/saipe/model-tables.html
- **Variables**: county poverty rate, median income
- **Years**: 2009-2013 average (pre-expansion baseline)

### Source 3: Census SAHIE (uninsured rate)
- **API**: https://www.census.gov/data/datasets/time-series/demo/sahie/estimates-acs.html
- **Variables**: pct uninsured, non-elderly adults, by county
- **Years**: 2013 (last pre-expansion year)
- **Effort**: CSV download + merge

### Source 4: MIT Election Lab (vote shares)
- **URL**: https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/VOQCHQ
- **Variables**: county-level presidential vote shares 2008, 2012
- **Compute**: Obama share = dem_votes / total_votes
- **Effort**: CSV download + merge by county FIPS

### Source 5: Democratic governor 2010
- **Source**: Manual coding (Wikipedia/Ballotpedia)
- **Variables**: binary indicator by state
- **Effort**: 50-state lookup table, ~15 minutes

### Source 6: Census of Governments (state health expenditures)
- **URL**: https://www.census.gov/programs-surveys/cog.html
- **Variables**: state-level health and welfare expenditures, 2005-2013 average
- **Effort**: Download + compute log average by state

### Source 7: Pre-expansion mortality (from existing data)
- **Already available**: crude_rate for 2009-2013
- **Compute**: average pre-expansion mortality by county (for PS model)
- **Hold-out**: exclude 2010-2013 from PS model, use for pre-trend verification

## Implementation plan

### Step 1: Data assembly script
**File**: `scripts/pull_bv_covariates.R`

```r
# Pull and merge all covariates needed for B&V DLPS replication
# Output: paper/data/bv_covariates.csv (one row per county)
#
# Sources:
#   1. Census ACS 5-year 2013 (tidycensus)
#   2. Census SAHIE 2013 (uninsured)
#   3. MIT Election Lab (presidential vote shares)
#   4. Democratic governor 2010 (manual)
#   5. Pre-expansion mortality averages (from analysis_data.csv)
```

Steps:
1. Pull ACS demographics + economics via tidycensus
2. Download SAHIE uninsured rates
3. Download MIT election data, compute Obama shares
4. Code Democratic governor indicator
5. Compute pre-expansion mortality averages from existing data
6. Merge all by county FIPS
7. Save as `paper/data/bv_covariates.csv`

### Step 2: Proper hdm double-lasso implementation
**File**: Update `paper/weighted_sdid_supplement.Rmd` Appendix F

```r
library(hdm)

# Step 1: Prepare baseline covariates (pre-treatment averages)
X_mat <- as.matrix(baseline[, all_covariates])
D_vec <- baseline$expansion
Y_vec <- baseline$avg_mortality_pre

# Step 2: Double-lasso variable selection
# Outcome model: Y ~ X
outcome_lasso <- rlasso(Y_vec ~ X_mat)
vars_outcome <- which(outcome_lasso$coefficients[-1] != 0)

# Treatment model: D ~ X
treatment_lasso <- rlassologit(D_vec ~ X_mat)
vars_treatment <- which(treatment_lasso$coefficients[-1] != 0)

# Union of selected variables
vars_selected <- union(vars_outcome, vars_treatment)

# Step 3: Propensity score from selected variables
ps_model <- glm(D_vec ~ X_mat[, vars_selected], family = binomial)
ps <- predict(ps_model, type = "response")

# Step 4: IPW weights (B&V formula)
# Treated: weight = pop
# Control: weight = pop * p/(1-p)
ipw <- ifelse(D_vec == 1, 1, ps / (1 - ps))
final_weight <- pop * ipw

# Step 5: Overlap trimming
# Drop counties outside [min(ps|D=1), max(ps|D=0)]
ps_overlap <- ps >= min(ps[D_vec==1]) & ps <= max(ps[D_vec==0])

# Step 6: Weighted DiD via fixest
fixest::feols(y ~ expansion:post | fips + year,
              data = panel_trimmed,
              weights = ~ final_weight,
              cluster = ~ state_fips)
```

### Step 3: Pre-trend verification
- Exclude 2010-2013 mortality from the propensity score model
- Estimate event-study with IPW weights
- Verify flat pre-trends in the held-out window (2010-2013)

### Step 4: Comparison table
Report side-by-side:
- B&V published estimate: -11.36
- Our DLPS (full covariates): [X]
- Our DLPS (limited covariates, current): [Y]
- Our population-weighted SDID: -13.5
- Our population-weighted DiD: -14.25

### Step 5: Sensitivity analysis
- Report with and without political variables
- Report with and without uninsured rate
- Show how the estimate changes as covariates are added sequentially

## Required R packages
- `tidycensus` (ACS pulls — needs Census API key)
- `hdm` (double-lasso)
- `fixest` (weighted DiD with clustering)
- `glmnet` (backup for PS estimation)

## Timeline estimate
- Data pulls (Sources 1-6): 2-3 hours (API setup + downloads + cleaning)
- Implementation (hdm + IPW + trimming): 1-2 hours
- Verification + comparison table: 1 hour
- Writing up results in supplement: 1 hour

## Files to create/modify
- `scripts/pull_bv_covariates.R` — new data pull script
- `paper/data/bv_covariates.csv` — new covariate file
- `paper/weighted_sdid_supplement.Rmd` — rewrite Appendix F with full replication
- `paper/references.bib` — add any new references

## Success criteria
- DLPS estimate within 3-4 deaths per 100,000 of B&V's -11.36
- Flat pre-trends in held-out window under IPW
- Clear documentation of which covariates we have vs. B&V
- Comparison table showing our SDID and DLPS estimates bracket B&V's result
