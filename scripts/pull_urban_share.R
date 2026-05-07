## Pull 2010 Census urban-share by county for the ACA application's heterogeneity
## binscatter. We use 2010 Decennial SF1 because it predates the 2014 expansion
## and is the canonical source for county-level urban-population shares.
##
## Output: data/urban_share.csv with columns (fips, pct_urban).

suppressPackageStartupMessages({
  library(tidycensus)
  library(dplyr)
  library(tidyr)
  library(readr)
})

OUT_CSV <- "data/urban_share.csv"

raw <- get_decennial(
  geography = "county",
  variables = c(total = "P002001", urban = "P002002"),
  year      = 2010,
  sumfile   = "sf1",
  output    = "wide"
)

urban <- raw %>%
  transmute(
    fips      = as.numeric(GEOID),
    pct_urban = ifelse(total > 0, urban / total, NA_real_)
  )

write_csv(urban, OUT_CSV)
message("Wrote ", OUT_CSV, " (", nrow(urban), " counties)")
