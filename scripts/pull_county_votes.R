## Pull county-level Obama vote shares 2008/2012, a propensity-score candidate
## used by Borgschulte & Vogler (2020) in the DLPS replication
## (supplement Appendix E.5; pipeline order in scripts/README.md).
##
## Preferred source: MIT Election Data and Science Lab, "County Presidential
## Election Returns" (Harvard Dataverse, doi:10.7910/DVN/VOQCHQ). That file is
## guestbook-protected, so it cannot be downloaded programmatically: download
## countypres_2000-2024.csv manually from
##   https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/VOQCHQ
## and save it as paper/data/raw/countypres_2000-2024.csv. If that file is
## present this script uses it.
##
## Fallback (automatic): tonmcg/US_County_Level_Election_Results_08-24 on
## GitHub (compiled from Guardian/townhall.com results). Vote totals differ
## trivially from MIT for a handful of counties; fine for a PS covariate, but
## switch to the MIT file before publication for citability.
##
## Output: paper/data/county_votes.csv with columns
##   (fips, obama_share_2008, obama_share_2012, votes_source).

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
})

OUT_CSV  <- "paper/data/county_votes.csv"
MIT_CSV  <- "paper/data/raw/countypres_2000-2024.csv"
GH_URL   <- paste0("https://raw.githubusercontent.com/tonmcg/",
                   "US_County_Level_Election_Results_08-24/master/",
                   "US_County_Level_Presidential_Results_08-16.csv")

if (file.exists(MIT_CSV)) {
  message("Using MIT Election Lab file: ", MIT_CSV)
  votes_raw <- read_csv(MIT_CSV, show_col_types = FALSE)

  obama <- votes_raw %>%
    filter(year %in% c(2008, 2012), office == "US PRESIDENT") %>%
    mutate(fips = as.numeric(county_fips)) %>%
    filter(!is.na(fips)) %>%
    group_by(year, fips) %>%
    summarise(dem_votes   = sum(candidatevotes[party == "DEMOCRAT"], na.rm = TRUE),
              total_votes = max(totalvotes, na.rm = TRUE),
              .groups = "drop") %>%
    mutate(obama_share = ifelse(total_votes > 0, dem_votes / total_votes, NA_real_)) %>%
    select(year, fips, obama_share) %>%
    pivot_wider(names_from = year, values_from = obama_share,
                names_prefix = "obama_share_") %>%
    mutate(votes_source = "mit_dataverse")

} else {
  message("MIT file not found at ", MIT_CSV, "; using tonmcg GitHub fallback")
  votes_raw <- read_csv(GH_URL, show_col_types = FALSE)

  ## Recompute totals from components: the file's total_* column has at least
  ## one transcription error (Laclede County MO, total_2008 = 2024 vs
  ## dem+gop+oth = 16,323), so component sums are the safer denominator.
  obama <- votes_raw %>%
    mutate(tot08 = dem_2008 + gop_2008 + oth_2008,
           tot12 = dem_2012 + gop_2012 + oth_2012) %>%
    transmute(fips = as.numeric(fips_code),
              obama_share_2008 = ifelse(tot08 > 0, dem_2008 / tot08, NA_real_),
              obama_share_2012 = ifelse(tot12 > 0, dem_2012 / tot12, NA_real_),
              votes_source = "tonmcg_github")
}

stopifnot(nrow(obama) > 3000,
          all(obama$obama_share_2008 >= 0 & obama$obama_share_2008 <= 1,
              na.rm = TRUE))

write_csv(obama, OUT_CSV)
message("Wrote ", OUT_CSV, " (", nrow(obama), " counties, source: ",
        obama$votes_source[1], ")")
