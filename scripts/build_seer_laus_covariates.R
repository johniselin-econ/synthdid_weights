## Build the SEER population covariates + LAUS unemployment for ALL counties,
## 2009-2017, from public raw files in paper/data/raw/ (downloaded from
## seer.cancer.gov/popdata and download.bls.gov; both gitignored for size).
##
## This replaces the covariate carryover in rebuild_mortality_panel.R: the 397
## counties recovered by the suppression imputation previously had NA
## covariates because the old panel never contained them. Definitions match
## the FIXED county_dataset_creation.R:
##   - pop_20_64 / log_20_64 use ages 20-64 (the original one-column across()
##     bug captured only 20-24)
##   - pct_55_64 uses the 20-64 denominator (BV2020 Table 1 convention)
##
## Output: paper/data/seer_laus_covariates.csv
##   (fips, year, state_abb, pop_total, pct_white, pop_20_64, pct_55_64,
##    log_20_64, log_35_44, log_f_20_64, unemp)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
})

RAW_DIR  <- "paper/data/raw"
OUT_CSV  <- "paper/data/seer_laus_covariates.csv"
## NOTE: SEER retired the 1990-2020 19ages vintage; current is 1990-2024 with
## 20 age groups. Age codes 5-13 (20-24 ... 60-64) are identical; only the
## 85+ tail is split differently, which we never use. Population estimates
## are Vintage 2024 and differ slightly from the Vintage-2020-based values in
## the original Mac pull -- document, don't be surprised.
SEER_GZ  <- file.path(RAW_DIR, "seer_us.1990_2024.20ages.adjusted.txt.gz")
LAUS_TXT <- file.path(RAW_DIR, "la.data.64.County")

START_YEAR <- 2009
END_YEAR   <- 2017

# ---------------------------------------------------------------------------
# 1. SEER population by county x year x race x sex x 19 age groups
# ---------------------------------------------------------------------------

message("Reading SEER population file (large; ~1-2 min)...")
df_pop <- read_fwf(
  SEER_GZ,
  fwf_widths(c(4, 2, 2, 3, 2, 1, 1, 1, 2, 8),
             c("year", "state_abb", "state_fips", "county_fips",
               "registry", "race", "hispanic", "sex", "age", "pop")),
  col_types = cols(year = col_double(), state_abb = col_character(),
                   state_fips = col_character(), county_fips = col_character(),
                   registry = col_character(), race = col_double(),
                   hispanic = col_double(), sex = col_double(),
                   age = col_double(), pop = col_double())
) %>%
  filter(year >= START_YEAR, year <= END_YEAR) %>%
  mutate(fips = as.numeric(paste0(state_fips, county_fips))) %>%
  select(year, state_abb, fips, race, sex, age, pop)

## Age codes: 0='00', 1='01-04', 2='05-09', ..., 4='15-19', 5='20-24',
## 6='25-29', 7='30-34', 8='35-39', 9='40-44', 10='45-49', 11='50-54',
## 12='55-59', 13='60-64', 14+='65+'

pop_age <- df_pop %>%
  group_by(year, state_abb, fips, age) %>%
  summarise(pop = sum(pop), .groups = "drop")

cov_age <- pop_age %>%
  group_by(year, state_abb, fips) %>%
  summarise(
    pop_total = sum(pop),
    pop_20_64 = sum(pop[age >= 5  & age <= 13]),
    pop_55_64 = sum(pop[age >= 12 & age <= 13]),
    pop_35_44 = sum(pop[age >= 8  & age <= 9]),
    .groups = "drop"
  ) %>%
  mutate(
    pct_55_64 = pop_55_64 / pop_20_64,        # 20-64 denominator (BV)
    log_20_64 = log(pop_20_64),
    log_35_44 = log(pop_35_44)
  )

## Race codes: 1=white, 2=black, 3=AI/AN, 4=API
cov_race <- df_pop %>%
  group_by(year, fips, race) %>%
  summarise(pop = sum(pop), .groups = "drop") %>%
  pivot_wider(names_from = race, values_from = pop, values_fill = 0,
              names_prefix = "race_") %>%
  mutate(pct_white = race_1 / (race_1 + race_2 + race_3 + race_4)) %>%
  select(year, fips, pct_white)

## Female 20-64 (sex code 2)
cov_f <- df_pop %>%
  filter(sex == 2, age >= 5, age <= 13) %>%
  group_by(year, fips) %>%
  summarise(log_f_20_64 = log(sum(pop)), .groups = "drop")

seer <- cov_age %>%
  left_join(cov_race, by = c("year", "fips")) %>%
  left_join(cov_f,    by = c("year", "fips")) %>%
  select(year, fips, state_abb, pop_total, pct_white, pop_20_64, pct_55_64,
         log_20_64, log_35_44, log_f_20_64)

message("SEER covariates: ", n_distinct(seer$fips), " counties x ",
        n_distinct(seer$year), " years")

# ---------------------------------------------------------------------------
# 2. LAUS county unemployment rate (series measure 03 = unemployment rate).
#    BLS blocks scripted downloads (403); fetch la.data.64.County manually in
#    a browser from download.bls.gov/pub/time.series/la/ and drop it in
#    paper/data/raw/. Until then unemp is left NA here and the panel rebuild
#    falls back to the old-panel carryover values.
# ---------------------------------------------------------------------------

if (file.exists(LAUS_TXT)) {
  message("Reading LAUS county file...")
  laus <- read_tsv(LAUS_TXT, col_types = cols(.default = col_character())) %>%
    mutate(year = as.numeric(year), value = as.numeric(value)) %>%
    filter(year >= START_YEAR, year <= END_YEAR) %>%
    mutate(measure     = substr(series_id, 19, 20),
           state_fips  = substr(series_id, 6, 7),
           county_fips = substr(series_id, 8, 10)) %>%
    filter(measure == "03") %>%
    group_by(state_fips, county_fips, year) %>%
    summarise(unemp = mean(value, na.rm = TRUE), .groups = "drop") %>%
    mutate(fips = as.numeric(paste0(state_fips, county_fips))) %>%
    select(year, fips, unemp)
  message("LAUS: ", n_distinct(laus$fips), " counties")
} else {
  message("LAUS file not found (", LAUS_TXT, "); unemp left NA -- ",
          "download la.data.64.County in a browser (BLS blocks scripts)")
  laus <- tibble(year = numeric(), fips = numeric(), unemp = numeric())
}

# ---------------------------------------------------------------------------
# 3. Merge and write
# ---------------------------------------------------------------------------

out <- seer %>%
  left_join(laus, by = c("year", "fips"))

stopifnot(nrow(out) == nrow(seer))
write_csv(out, OUT_CSV)
message("Wrote ", OUT_CSV, " (", nrow(out), " county-years; unemp missing for ",
        sum(is.na(out$unemp)), ")")
