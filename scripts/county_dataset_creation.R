# =============================================================================
# Author: John Iselin and Erica Ryan
# Date:   November 24th, 2025
# File:   dataset_construction
#
# Project:  Dataset Construction for Borgschulte and Vogler (2020) replication
#           Construct County-Year panel from 2009 through 2017.
#           The original paper used confidential mortality data, which we do not
#           have access to. Rather, we use public data via CDC wonder.
# =============================================================================

## CLEAR ALL
rm(list = ls())
gc()


# --- Libraries ---
library(here)       # project-root relative paths
library(tidyverse)
library(data.table)
library(readr)      # data import
library(stringr)    # string helpers

# --- Project Metadata ---
project_name <- "synth_weights"
date_run     <- format(Sys.Date(), "%Y-%m-%d")

# --- File Paths (relative to Rproj root) ---
dir_data     <- here("paper", "data")
dir_data_raw <- "/Users/johniselin/Library/CloudStorage/Dropbox/sdid_weights/data/"

# --- Random Seed for Reproducibility ---
set.seed(56403)

# --- Parameters ---
start_year       <- 2009
end_year         <- 2017

## ACA Expansion states (First half of 2014)
aca <- c(4, 5, 6, 8, 9, 10, 11,	15, 17, 19, 21, 24, 25, 26, 27, 32, 34,
         35, 36, 38, 39, 41, 44, 50, 53, 54)

## DROP AK IN, LA, MT, NH, and PA for later expansions (post first half of 2014)
drop <- c(2, 18, 22, 30, 33, 42)

## Mortality Data
## Imported via CDC Wonder on November 24, 2025
## All cause mortality, 20-64-year-olds by county and year

## Suppression handling (2026-06 rework). CDC WONDER reports "Suppressed" for
## county-year cells with 1-9 deaths and a numeric 0 for true zeros, so the two
## are distinguishable in the raw export. The original version zero-filled
## suppressed cells and then dropped any county without 9 years of data, which
## removed ~380 (mostly small, rural, high-mortality, non-expansion) counties
## relative to BV2020's 2,823-county pre-trim sample. We now impute suppressed
## cells instead of dropping those counties.
##
## impute_method options:
##   "residual" - allocate the state-year residual (exact state total minus the
##                sum of unsuppressed county deaths) across suppressed cells in
##                proportion to population, clamped to [1, 9]. Requires a
##                state-level WONDER export at dir_data_raw/cdc_wonder_states.csv
##                (state totals are never suppressed for all-cause 20-64).
##                Falls back to "poisson" if that file is absent.
##   "poisson"  - E[D | 1 <= D <= 9] under Poisson with lambda = population x
##                state-year rate estimated from unsuppressed counties.
##   "lower"/"upper" - bound runs (deaths = 1 / deaths = 9) for sensitivity.
impute_suppressed <- TRUE
impute_method     <- "residual"

## Load data
df_mortality <- read.csv(file = file.path(dir_data_raw,"cdc_wonder_data.csv"))

## Rename all to lower case
names(df_mortality) <- tolower(names(df_mortality))

## Clean Data
df_mortality <- df_mortality %>%

  ## Drop crude rate
  select(-crude.rate, -year.code, -notes) %>%

  ## Tag counties with suppression
  mutate(death_supp = ifelse(deaths == "Suppressed", 1,0),
         death_miss = ifelse(deaths == "Missing", 1, 0),
         pop_miss = ifelse(population == "Missing", 1, 0)) %>%

  ## Suppressed/missing deaths -> NA (imputed below); missing population -> NA
  mutate(deaths = ifelse(deaths %in% c("Suppressed", "Missing"), NA, deaths),
         population = ifelse(population == "Missing", NA, population)) %>%

  ## Rename geographic vars
  rename(fips = county.code) %>%

  ## Convert deaths and population to numeric
  mutate(deaths = as.numeric(deaths), population = as.numeric(population)) %>%
  mutate(state_fips_mort = fips %/% 1000)

## Impute suppressed cells (1-9 deaths by construction)
if (impute_suppressed) {

  ## State-year rates from unsuppressed counties (poisson fallback + weights)
  state_rates <- df_mortality %>%
    filter(death_supp == 0, death_miss == 0, !is.na(population)) %>%
    group_by(state_fips_mort, year) %>%
    summarise(state_rate = sum(deaths) / sum(population) * 100000,
              .groups = "drop")

  df_mortality <- df_mortality %>%
    left_join(state_rates, by = c("state_fips_mort", "year"))

  trunc_pois_mean <- function(lambda) {
    ## E[D | 1 <= D <= 9] under Poisson(lambda)
    d <- 1:9
    sapply(lambda, function(l) {
      if (is.na(l) || l <= 0) return(4.5)
      p <- dpois(d, l)
      if (sum(p) == 0) return(ifelse(l > 9, 9, 1))
      sum(d * p) / sum(p)
    })
  }

  state_file <- file.path(dir_data_raw, "cdc_wonder_states.csv")
  method_used <- impute_method
  if (impute_method == "residual" && !file.exists(state_file)) {
    message("cdc_wonder_states.csv not found; falling back to poisson imputation")
    method_used <- "poisson"
  }

  if (method_used == "residual") {
    ## Exact state-year totals from the state-level WONDER export
    df_states <- read.csv(state_file)
    names(df_states) <- tolower(names(df_states))
    df_states <- df_states %>%
      transmute(state_fips_mort = as.numeric(state.code),
                year = year,
                state_deaths = as.numeric(deaths))

    residuals_sy <- df_mortality %>%
      group_by(state_fips_mort, year) %>%
      summarise(observed = sum(deaths[death_supp == 0 & death_miss == 0],
                               na.rm = TRUE),
                supp_pop = sum(population[death_supp == 1], na.rm = TRUE),
                .groups = "drop") %>%
      left_join(df_states, by = c("state_fips_mort", "year")) %>%
      mutate(residual = state_deaths - observed)

    df_mortality <- df_mortality %>%
      left_join(residuals_sy %>% select(state_fips_mort, year, residual, supp_pop),
                by = c("state_fips_mort", "year")) %>%
      mutate(deaths_imp = ifelse(
        death_supp == 1 & !is.na(population) & supp_pop > 0,
        pmin(pmax(residual * population / supp_pop, 1), 9),
        NA_real_)) %>%
      select(-residual, -supp_pop)
  } else if (method_used == "poisson") {
    df_mortality <- df_mortality %>%
      mutate(deaths_imp = ifelse(
        death_supp == 1 & !is.na(population),
        trunc_pois_mean(population * state_rate / 100000),
        NA_real_))
  } else if (method_used == "lower") {
    df_mortality <- df_mortality %>%
      mutate(deaths_imp = ifelse(death_supp == 1, 1, NA_real_))
  } else if (method_used == "upper") {
    df_mortality <- df_mortality %>%
      mutate(deaths_imp = ifelse(death_supp == 1, 9, NA_real_))
  } else stop("Unknown impute_method: ", impute_method)

  df_mortality <- df_mortality %>%
    mutate(deaths = ifelse(death_supp == 1, deaths_imp, deaths)) %>%
    select(-deaths_imp, -state_rate)
}

df_mortality <- df_mortality %>%

  ## Create mortality rate
  mutate(crude_rate = (deaths / population) * 100000 ) %>%

  ## Select data (death_supp kept so downstream can flag/exclude imputed cells)
  select(year, fips, crude_rate, deaths, population, death_supp)


## Population data via https://seer.cancer.gov/popdata/download.html
## 1990-2020, 4 Expanded Races by Origin, All US
df_pop <- read_fwf(file = file.path(dir_data_raw,"covariates/us.1990_2020.19ages.adjusted.txt"),
                     fwf_widths( c(4, 2, 2, 3, 2, 1, 1, 1, 2, 8),
                                 c("year", "state_abb", "state_fips", "county_fips",
                                   "registry", "race", "hispanic", "sex", "age", "pop"))) %>%
    filter(year >= 2009) %>% filter(year <= 2017) %>%
    mutate(age = as.numeric(age), pop = as.numeric(pop)) %>%
    select(-registry, -hispanic) %>%
    group_by(year, state_abb, state_fips, county_fips, race, sex, age) %>%
    summarize(pop = sum(pop)) %>% group_by() %>%
    mutate(fips = as.numeric(paste0(state_fips, county_fips)))

## Create relevant population samples

## Percent of population that is white
df_pop_race <- df_pop %>%

  ## Keep required variables
  select(year, state_abb, state_fips, county_fips, fips, race, pop) %>%

  ## Group to get total pop by race
  group_by(year, state_abb, state_fips, county_fips, race, fips) %>%
  summarize(pop = sum(pop)) %>% group_by() %>%

  ## Pivot to county X year level
  pivot_wider(names_from = race, values_from = pop, values_fill = 0) %>%

  ## Rename by race
  rename(white = '1', black = '2', aian = '3', asian = '4') %>%

  ## Fill-in missing values
  replace_na(list(white = 0, black = 0, aian = 0, asian = 0)) %>%

  ## Create percent white
  mutate(pct_white = white / (white + black + aian + asian)) %>%
  select(year, fips, state_abb, state_fips, county_fips, pct_white) %>%
  ungroup()

## Percent of the population that is elderly
df_pop_age <- df_pop %>%

  ## Keep required variables
  select(year, state_abb, state_fips, county_fips, age, pop) %>%

  ## Group by to get totals by age
  group_by(year, state_abb, state_fips, county_fips, age) %>%
  summarize(pop = sum(pop)) %>% group_by() %>%

  ## Pivot to get county X year level
  pivot_wider(names_from = age, values_from = pop, values_fill = 0) %>%

  ## Rename variables
  rename(pop_00 = '0', pop_01_04 = '1', pop_05_09 = '2',
         pop_10_14 = '3', pop_15_19 = '4', pop_20_24 = '5',
         pop_25_29 = '6', pop_30_34 = '7', pop_35_39 = '8',
         pop_40_44 = '9', pop_45_49 = '10', pop_50_54 = '11',
         pop_55_59 = '12', pop_60_64 = '13', pop_65_69 = '14',
         pop_70_74 = '15', pop_75_79 = '16', pop_80_84 = '17',
         pop_85 = '18' ) %>%

  ## Create summary variables
  ## NOTE (2026-06): two bug fixes relative to the original version.
  ##   (1) pop_20_64 / log_20_64 previously used across(pop_20_24:pop_20_24),
  ##       a one-column range, so they captured only ages 20-24.
  ##   (2) pct_55_64 previously used pop_total as the denominator; BV2020's
  ##       five age shares sum to 1 within the 20-64 population (Table 1),
  ##       so the denominator is the working-age population.
  ## analysis_data.csv must be regenerated before these fixes reach the paper.
  mutate(pop_total = rowSums(across(pop_00:pop_85)),
         pop_20_64 = rowSums(across(pop_20_24:pop_60_64))) %>%
  mutate(pct_55_64 = (pop_55_59 + pop_60_64) / pop_20_64,
         log_20_64 = log(pop_20_64),
         log_35_44 = log(rowSums(across(pop_35_39:pop_40_44)))) %>%
  select(year, state_abb, state_fips, county_fips, pct_55_64, log_20_64, log_35_44, pop_20_64, pop_total) %>%
  ungroup()

## By-age and gender statistics
df_pop_f_age <- df_pop %>%

  ## Keep if sex == 2
  filter(sex == 2) %>%

  ## Keep required variables
  select(year, state_abb, state_fips, county_fips, age, pop) %>%

  ## Group by to get totals by age
  group_by(year, state_abb, state_fips, county_fips, age) %>%
  summarize(pop = sum(pop)) %>% group_by() %>%

  ## Pivot to get county X year level
  pivot_wider(names_from = age, values_from = pop, values_fill = 0) %>%
  rename(pop_00 = '0', pop_01_04 = '1', pop_05_09 = '2',
         pop_10_14 = '3',  pop_15_19 = '4', pop_20_24 = '5',
         pop_25_29 = '6', pop_30_34 = '7', pop_35_39 = '8',
         pop_40_44 = '9', pop_45_49 = '10', pop_50_54 = '11',
         pop_55_59 = '12', pop_60_64 = '13', pop_65_69 = '14',
         pop_70_74 = '15', pop_75_79 = '16', pop_80_84 = '17',
         pop_85 = '18' ) %>%

  ## Get summary variables
  mutate(pop_total = rowSums(across(pop_00:pop_85))) %>%
  mutate(log_f_20_64 = log(rowSums(across(pop_20_24:pop_60_64)))) %>%
  select(year, state_abb, state_fips, county_fips, log_f_20_64) %>%
  ungroup()

## Merge together population covariates
df_pop_cov <- df_pop_race %>%
  full_join(df_pop_age,
            by = join_by(year, state_abb, state_fips, county_fips)) %>%
    left_join(df_pop_f_age,
              by = join_by(year, state_abb, state_fips, county_fips)) %>%
    select(year, fips, state_abb, state_fips, county_fips, pop_total, pct_white,
           pop_20_64, pct_55_64, log_20_64, log_35_44, log_f_20_64)

## Unemployment data via https://download.bls.gov/pub/time.series/la/
## la.data.64.County
df_unempl <- read_tsv(file = file.path(dir_data_raw,"covariates/la.data.64.County")) %>%

  ## Keep year range
  filter(year >= 2009, year <= 2017) %>%
    select(-footnote_codes) %>%
    mutate(series = as.numeric(substr(series_id, 19, 20)),
           state_fips = substr(series_id, 6, 7),
           county_fips = substr(series_id, 8, 10)) %>%
    filter(series == 3) %>%
    group_by(state_fips, county_fips, series_id, year) %>%
    summarise(unemp = mean(value)) %>%
    group_by() %>%
    mutate(fips = as.numeric(paste0(state_fips, county_fips)),
           state_fips = as.numeric(state_fips)) %>%
    select(year, fips, unemp)

## Combine data
df <- df_pop_cov %>%

  ## Merge with Unemployment
  left_join(df_unempl, join_by(year, fips)) %>%

  ## Merge with Mortality
  left_join(df_mortality, join_by(year, fips)) %>%

  ## Assign ACA variables
  mutate(expansion = ifelse(state_fips %in% aca,1, 0),
         post = ifelse(year >= 2014,1,0)) %>%
  mutate(treated = post * expansion) %>%

  ## Drop if in select states with odd ACA timing
  filter(!(state_fips %in% drop)) %>%

  ## Keep years with mortality data
  filter(year >= 2009 & year <= 2017) %>%

  ## Drop if missing mortality data
  group_by(fips) %>%
  mutate(ever_missing_flag = as.numeric(any(is.na(crude_rate))),
         count_year = n()) %>%
  ungroup() %>%

  ## Keep fips with no missing data
  filter(ever_missing_flag == 0,
         count_year == 9 )


## Determine state assignment

## Export version of data
write.csv(df,file = file.path(dir_data,"analysis_data.csv"))

rm(list = ls())
