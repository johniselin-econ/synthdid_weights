## Pull the Borgschulte & Vogler (2020) covariate menu via public sources for
## use as inputs to the DLPS replication. Output: paper/data/bv_covariates.csv.
##
## Sources covered:
##   - ACS 2009-2013 5-year (tidycensus): demographics, age groups, race,
##     poverty, median income, total population.
##   - tigris county shapefile 2010: land area (for population density).
##   - Census SAHIE 2013: uninsured rate, non-elderly adults (ages 18-64).
##   - Pre-expansion mortality 2009-2013: rolled up from paper/data/analysis_data.csv.
##   - Democratic governor 2010: state lookup table inline.
##
## Sources NOT yet covered (would-be additional improvements; flagged TODO):
##   - MIT Election Lab county vote shares 2008/2012.
##   - Census of Governments state health & welfare expenditures 2005-2013.

suppressPackageStartupMessages({
  library(tidycensus)
  library(tigris)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(httr)
  library(jsonlite)
})

options(tigris_use_cache = TRUE)

OUT_CSV <- "paper/data/bv_covariates.csv"

# ---------------------------------------------------------------------------
# 1. ACS 2009-2013 5-year demographics (variables match BV2020 covariate menu)
# ---------------------------------------------------------------------------

acs_vars <- c(
  total      = "B01001_001",
  male       = "B01001_002",
  # Male age cells we sum to BV's 5 groups
  m_20a      = "B01001_008",   # M 20
  m_20b      = "B01001_009",   # M 21
  m_20c      = "B01001_010",   # M 22-24
  m_2529     = "B01001_011",
  m_3034     = "B01001_012",
  m_3539     = "B01001_013",
  m_4044     = "B01001_014",
  m_4549     = "B01001_015",
  m_5054     = "B01001_016",
  m_5559     = "B01001_017",
  m_6061     = "B01001_018",
  m_6264     = "B01001_019",
  # Female age cells
  f_20a      = "B01001_032",   # F 20
  f_20b      = "B01001_033",   # F 21
  f_20c      = "B01001_034",   # F 22-24 — note: ACS 5-yr 2013 puts 22-24 here
  f_2529     = "B01001_035",
  f_3034     = "B01001_036",
  f_3539     = "B01001_037",
  f_4044     = "B01001_038",
  f_4549     = "B01001_039",
  f_5054     = "B01001_040",
  f_5559     = "B01001_041",
  f_6061     = "B01001_042",
  f_6264     = "B01001_043",
  # Race / Hispanic origin (B03002)
  hisp_total = "B03002_001",
  nh_white   = "B03002_003",
  nh_black   = "B03002_004",
  hispanic   = "B03002_012",
  # Poverty (S1701-style direct counts via B17001)
  pov_universe = "B17001_001",
  pov_below    = "B17001_002",
  # Median household income
  med_income = "B19013_001"
)

acs_raw <- get_acs(
  geography = "county",
  variables = acs_vars,
  year      = 2013,
  survey    = "acs5",
  output    = "wide"
)

acs <- acs_raw %>%
  transmute(
    fips         = as.numeric(GEOID),
    total_pop    = totalE,
    pct_male     = maleE / totalE,
    pop_2024     = m_20aE + m_20bE + m_20cE + f_20aE + f_20bE + f_20cE,
    pop_2534     = m_2529E + m_3034E + f_2529E + f_3034E,
    pop_3544     = m_3539E + m_4044E + f_3539E + f_4044E,
    pop_4554     = m_4549E + m_5054E + f_4549E + f_5054E,
    pop_5564     = m_5559E + m_6061E + m_6264E + f_5559E + f_6061E + f_6264E,
    pop_2064     = pop_2024 + pop_2534 + pop_3544 + pop_4554 + pop_5564,
    # BV2020 Table 1 reports age shares as fractions of working-age (20-64)
    # population, not total population (their five shares sum to ~1).
    pct_2024     = pop_2024 / pop_2064,
    pct_2534     = pop_2534 / pop_2064,
    pct_3544     = pop_3544 / pop_2064,
    pct_4554     = pop_4554 / pop_2064,
    pct_5564     = pop_5564 / pop_2064,
    pct_white    = nh_whiteE   / hisp_totalE,
    pct_black    = nh_blackE   / hisp_totalE,
    pct_hispanic = hispanicE   / hisp_totalE,
    poverty_rate = pov_belowE  / pov_universeE,
    med_income   = med_incomeE,
    log_median_income = log(pmax(med_incomeE, 1))
  )

# ---------------------------------------------------------------------------
# 2. Land area for population density (2010 county shapefile)
# ---------------------------------------------------------------------------

land <- counties(cb = TRUE, year = 2010, progress_bar = FALSE) %>%
  sf::st_drop_geometry() %>%
  transmute(
    fips      = as.numeric(paste0(STATE, COUNTY)),
    aland_sqmi = CENSUSAREA   # CB shapefile gives land area in sq miles
  )

# ---------------------------------------------------------------------------
# 3. SAHIE 2013: uninsured rate, ages 18-64, all races, both sexes
#    Census API: api.census.gov/data/timeseries/healthins/sahie
# ---------------------------------------------------------------------------

sahie_url <- paste0(
  "https://api.census.gov/data/timeseries/healthins/sahie?",
  "get=NAME,NUI_PT,NIPR_PT,PCTUI_PT&",
  "for=county:*&",
  # AGECAT=3 = ages 18-64 (closest match to BV's "non-elderly adults" 19-64)
  "AGECAT=3&RACECAT=0&SEXCAT=0&IPRCAT=0&",
  "time=2013"
)
sahie_resp <- httr::GET(sahie_url, httr::timeout(60))
stopifnot(httr::status_code(sahie_resp) == 200)
sahie_raw <- jsonlite::fromJSON(rawToChar(sahie_resp$content))
sahie <- as.data.frame(sahie_raw[-1, , drop = FALSE], stringsAsFactors = FALSE)
colnames(sahie) <- sahie_raw[1, ]
sahie <- sahie %>%
  mutate(
    fips           = as.numeric(paste0(state, county)),
    uninsured_rate = as.numeric(PCTUI_PT) / 100  # API returns 0-100, want 0-1
  ) %>%
  select(fips, uninsured_rate)

# ---------------------------------------------------------------------------
# 4. Pre-expansion mortality average 2009-2013 from existing analysis data
# ---------------------------------------------------------------------------

panel <- read_csv("paper/data/analysis_data.csv", show_col_types = FALSE)
mort_pre <- panel %>%
  filter(year >= 2009, year <= 2013) %>%
  group_by(fips) %>%
  summarise(avg_mortality_pre = mean(crude_rate, na.rm = TRUE), .groups = "drop")

# ---------------------------------------------------------------------------
# 5. Democratic governor 2010 indicator (lookup table at state level)
#    Source: National Governors Association / Wikipedia, cross-checked
# ---------------------------------------------------------------------------

dem_gov_2010 <- tribble(
  ~state_fips, ~dem_governor_2010,
   1, 0,  2, 0,  4, 0,  5, 1,  6, 0,  8, 1,  9, 1, 10, 1, 11, 1, 12, 0,
  13, 0, 15, 1, 16, 0, 17, 1, 18, 0, 19, 1, 20, 0, 21, 1, 22, 1, 23, 0,
  24, 1, 25, 1, 26, 0, 27, 0, 28, 0, 29, 1, 30, 1, 31, 0, 32, 0, 33, 1,
  34, 0, 35, 0, 36, 1, 37, 1, 38, 1, 39, 1, 40, 0, 41, 1, 42, 1, 44, 1,
  45, 0, 46, 0, 47, 1, 48, 0, 49, 0, 50, 0, 51, 1, 53, 1, 54, 1, 55, 0, 56, 0
)

# ---------------------------------------------------------------------------
# 6. Merge everything
# ---------------------------------------------------------------------------

bv_cov <- acs %>%
  left_join(land,     by = "fips") %>%
  left_join(sahie,    by = "fips") %>%
  left_join(mort_pre, by = "fips") %>%
  mutate(
    state_fips = floor(fips / 1000),
    pop_density   = total_pop / pmax(aland_sqmi, 1),
    log_pop_density = log(pmax(pop_density, 1))
  ) %>%
  left_join(dem_gov_2010, by = "state_fips") %>%
  select(
    fips, state_fips,
    pct_male, pct_white, pct_black, pct_hispanic,
    pct_2024, pct_2534, pct_3544, pct_4554, pct_5564,
    poverty_rate, med_income, log_median_income,
    pop_density, log_pop_density,
    uninsured_rate,
    avg_mortality_pre,
    dem_governor_2010
  )

write_csv(bv_cov, OUT_CSV)
message("Wrote ", OUT_CSV, " (", nrow(bv_cov), " counties; ",
        sum(complete.cases(bv_cov)), " complete-case)")
