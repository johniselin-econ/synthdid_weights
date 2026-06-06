## Rebuild the analysis panel from the public CDC WONDER county exports in
## paper/data/raw/ (2005-2011 + 2012-2017, exported with "Show Suppressed" and
## "Show Zero Values"). This is the public-data rebuild that runs anywhere;
## scripts/county_dataset_creation.R remains the from-scratch builder that
## additionally constructs the SEER/LAUS covariates from raw files on the Mac.
##
## What it does:
##   1. Combine the split county exports; impute suppressed cells (1-9 deaths)
##      by allocating the exact state-year residual (state total minus
##      unsuppressed county sum, validated feasible for all 663 state-years)
##      proportional to population, clamped to [1, 9].
##   2. Build treatment variables with NUMERIC state fips (the character/
##      numeric %in% mismatch previously misclassified AZ/AR/CA/CO/CT as
##      control and let AK escape the drop filter).
##   3. Carry over the SEER/LAUS covariates from the previous analysis_data
##      by county-year (counties new to the panel get NA covariates until
##      county_dataset_creation.R is re-run against the raw SEER/LAUS files;
##      SDID main results never use those covariates).
##   4. Write paper/data/analysis_data.csv (2009-2017 balanced panel) and
##      paper/data/county_mortality_pre.csv (2005-2008, for the PS model's
##      year-by-year pre-expansion mortality with 2010-2013 held out).
##
## Sensitivity: set IMPUTE_METHOD to "lower"/"upper" for deaths=1/deaths=9
## bound runs, or "poisson" for the truncated-Poisson fallback.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

RAW_DIR       <- "paper/data/raw"
OUT_PANEL     <- "paper/data/analysis_data.csv"
OUT_PRE       <- "paper/data/county_mortality_pre.csv"
OLD_PANEL     <- "paper/data/analysis_data.csv"   # source of carried covariates
IMPUTE_METHOD <- "residual"

START_YEAR <- 2009
END_YEAR   <- 2017

## ACA expansion states (first half of 2014), incl. DC
aca  <- c(4, 5, 6, 8, 9, 10, 11, 15, 17, 19, 21, 24, 25, 26, 27, 32, 34,
          35, 36, 38, 39, 41, 44, 50, 53, 54)
## Drop AK, IN, LA, MT, NH, PA (later expansions)
drop <- c(2, 18, 22, 30, 33, 42)

# ---------------------------------------------------------------------------
# 1. Read county exports and the state-level totals
# ---------------------------------------------------------------------------

read_wonder <- function(f) {
  d <- read.csv(f)
  names(d) <- tolower(names(d))
  d %>% filter(county.code != "")   # drops the metadata footer
}

mort <- bind_rows(
  read_wonder(file.path(RAW_DIR, "cdc_wonder_data_2005-2011.csv")),
  read_wonder(file.path(RAW_DIR, "cdc_wonder_data_2012-2017.csv"))
) %>%
  transmute(
    fips       = as.numeric(county.code),
    state_fips = fips %/% 1000,
    year       = as.numeric(year),
    death_supp = as.integer(deaths == "Suppressed"),
    death_miss = as.integer(deaths == "Missing"),
    deaths     = suppressWarnings(as.numeric(ifelse(
                   deaths %in% c("Suppressed", "Missing"), NA, deaths))),
    population = suppressWarnings(as.numeric(ifelse(
                   population == "Missing", NA, population)))
  )

stopifnot(nrow(mort) > 0, !any(duplicated(mort[, c("fips", "year")])))

states <- read.csv(file.path(RAW_DIR, "cdc_wonder_states.csv"))
names(states) <- tolower(names(states))
states <- states %>%
  transmute(state_fips   = suppressWarnings(as.numeric(state.code)),
            year         = suppressWarnings(as.numeric(year)),
            state_deaths = suppressWarnings(as.numeric(deaths))) %>%
  filter(!is.na(state_fips))

# ---------------------------------------------------------------------------
# 2. Impute suppressed cells
# ---------------------------------------------------------------------------

if (IMPUTE_METHOD == "residual") {

  resid_sy <- mort %>%
    group_by(state_fips, year) %>%
    summarise(observed = sum(deaths[death_supp == 0 & death_miss == 0],
                             na.rm = TRUE),
              supp_pop = sum(population[death_supp == 1], na.rm = TRUE),
              n_supp   = sum(death_supp),
              .groups = "drop") %>%
    left_join(states, by = c("state_fips", "year")) %>%
    mutate(residual = state_deaths - observed)

  ## Feasibility guard: residual must lie in [n_supp, 9 * n_supp]
  infeasible <- resid_sy %>%
    filter(n_supp > 0,
           (residual < n_supp | residual > 9 * n_supp))
  if (nrow(infeasible) > 0) {
    print(infeasible)
    stop("Residual allocation infeasible for the state-years above; ",
         "check that county and state exports use the same query.")
  }

  mort <- mort %>%
    left_join(resid_sy %>% select(state_fips, year, residual, supp_pop),
              by = c("state_fips", "year")) %>%
    mutate(deaths = ifelse(
      death_supp == 1 & !is.na(population) & supp_pop > 0,
      pmin(pmax(residual * population / supp_pop, 1), 9),
      deaths)) %>%
    select(-residual, -supp_pop)

} else if (IMPUTE_METHOD == "poisson") {
  state_rates <- mort %>%
    filter(death_supp == 0, death_miss == 0, !is.na(population)) %>%
    group_by(state_fips, year) %>%
    summarise(state_rate = sum(deaths) / sum(population) * 100000,
              .groups = "drop")
  trunc_pois_mean <- function(lambda) {
    d <- 1:9
    sapply(lambda, function(l) {
      if (is.na(l) || l <= 0) return(4.5)
      p <- dpois(d, l)
      if (sum(p) == 0) return(ifelse(l > 9, 9, 1))
      sum(d * p) / sum(p)
    })
  }
  mort <- mort %>%
    left_join(state_rates, by = c("state_fips", "year")) %>%
    mutate(deaths = ifelse(
      death_supp == 1 & !is.na(population),
      trunc_pois_mean(population * state_rate / 100000),
      deaths)) %>%
    select(-state_rate)
} else if (IMPUTE_METHOD %in% c("lower", "upper")) {
  mort <- mort %>%
    mutate(deaths = ifelse(death_supp == 1,
                           ifelse(IMPUTE_METHOD == "lower", 1, 9), deaths))
} else stop("Unknown IMPUTE_METHOD: ", IMPUTE_METHOD)

mort <- mort %>%
  mutate(crude_rate = deaths / population * 100000)

message(sprintf("Imputed %d suppressed county-years (%s method)",
                sum(mort$death_supp == 1 & !is.na(mort$deaths)),
                IMPUTE_METHOD))

# ---------------------------------------------------------------------------
# 3. Treatment variables (numeric state fips!) and sample restrictions
# ---------------------------------------------------------------------------

mort <- mort %>%
  filter(!(state_fips %in% drop)) %>%
  mutate(expansion = as.integer(state_fips %in% aca),
         post      = as.integer(year >= 2014),
         treated   = post * expansion)

## Analysis years; require a complete balanced panel of valid rates
panel <- mort %>%
  filter(year >= START_YEAR, year <= END_YEAR) %>%
  group_by(fips) %>%
  mutate(ever_missing_flag = as.integer(any(is.na(crude_rate))),
         count_year        = n()) %>%
  ungroup() %>%
  filter(ever_missing_flag == 0,
         count_year == (END_YEAR - START_YEAR + 1))

# ---------------------------------------------------------------------------
# 4. Covariates: freshly built SEER/LAUS file when available (covers ALL
#    counties, fixed definitions; scripts/build_seer_laus_covariates.R),
#    otherwise carry over from the previous panel. unemp falls back to the
#    carryover where LAUS is missing (BLS blocks scripted downloads).
# ---------------------------------------------------------------------------

old <- read_csv(OLD_PANEL, show_col_types = FALSE) %>%
  select(any_of(c("fips", "year", "state_abb", "pop_total", "pct_white",
                  "pop_20_64", "pct_55_64", "log_20_64", "log_35_44",
                  "log_f_20_64", "unemp"))) %>%
  distinct(fips, year, .keep_all = TRUE)

state_abb_lookup <- old %>%
  filter(!is.na(state_abb)) %>%
  distinct(state_fips = floor(fips / 1000), state_abb)

SEER_COV <- "paper/data/seer_laus_covariates.csv"
if (file.exists(SEER_COV)) {
  seer <- read_csv(SEER_COV, show_col_types = FALSE)
  panel <- panel %>%
    left_join(seer %>% select(-any_of("state_abb")), by = c("fips", "year")) %>%
    left_join(old %>% select(fips, year, unemp_old = unemp),
              by = c("fips", "year")) %>%
    mutate(unemp = coalesce(unemp, unemp_old)) %>%
    select(-unemp_old) %>%
    left_join(state_abb_lookup, by = "state_fips")
  message("Covariates: fresh SEER build (+ LAUS/carryover unemp)")
} else {
  panel <- panel %>%
    left_join(old %>% select(-state_abb), by = c("fips", "year")) %>%
    left_join(state_abb_lookup, by = "state_fips")
  message("Covariates: carryover from previous panel (SEER build not found)")
}

n_no_cov <- panel %>% filter(is.na(pct_white)) %>% distinct(fips) %>% nrow()
message(sprintf(
  "Panel: %d counties (%d expansion / %d control); %d counties lack SEER/LAUS covariates (rebuild on Mac to fill)",
  n_distinct(panel$fips),
  n_distinct(panel$fips[panel$expansion == 1]),
  n_distinct(panel$fips[panel$expansion == 0]),
  n_no_cov))

## Drop the deaths count from the public file: WONDER's data-use agreement
## prohibits publishing sub-national counts of 1-9, and although imputed
## values are model estimates rather than actual counts, shipping only the
## rate (death_supp flags which cells are imputed) keeps us clearly inside
## the agreement. No downstream code uses the count.
panel <- panel %>% select(-deaths, -death_miss)

write_csv(panel, OUT_PANEL)
message("Wrote ", OUT_PANEL)

# ---------------------------------------------------------------------------
# 5. Pre-2009 mortality (PS model inputs; 2010-2013 stay held out)
# ---------------------------------------------------------------------------

pre <- mort %>%
  filter(year >= 2005, year <= 2008) %>%
  select(fips, year, crude_rate, population, death_supp)  # no counts (DUA)

write_csv(pre, OUT_PRE)
message("Wrote ", OUT_PRE, " (", n_distinct(pre$fips), " counties, 2005-2008)")
