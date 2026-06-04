## Count papers using DiD / SCM / SDID in title or abstract, by year × journal,
## from the AEA Replication Tracker SQLite produced by the
## `paulgp/replication-package-db` pipeline (extended to include QJE, JPE,
## ReStud, Econometrica). Output is `paper/data/lit_counts.csv`, which is read by the
## paper at render time so we don't re-query the SQLite on every knit.

suppressPackageStartupMessages({
  library(DBI)
  library(RSQLite)
  library(dplyr)
  library(stringr)
  library(readr)
  library(tidyr)
})

DB_PATH <- "../replication-package-db/data/replication_tracker.db"
OUT_CSV <- "paper/data/lit_counts.csv"

if (!file.exists(DB_PATH)) {
  stop("SQLite not found at ", DB_PATH,
       " — run scripts/02_fetch_papers.py in replication-package-db first.")
}

con <- dbConnect(SQLite(), DB_PATH)
on.exit(dbDisconnect(con), add = TRUE)

papers <- dbGetQuery(con, "
  SELECT openalex_id, journal_name, publication_year,
         COALESCE(title, '') AS title,
         COALESCE(abstract, '') AS abstract
  FROM papers
  WHERE publication_year BETWEEN 2014 AND 2024
") |> as_tibble()

## Tier journals into Top-5 vs AEJ Field for the figure.
top5_journals <- c(
  "American Economic Review",
  "The Quarterly Journal of Economics",
  "Journal of Political Economy",
  "The Review of Economic Studies",
  "Econometrica"
)
aej_journals <- c(
  "American Economic Journal Applied Economics",
  "American Economic Journal Economic Policy",
  "American Economic Journal Macroeconomics",
  "American Economic Journal Microeconomics"
)

papers <- papers |>
  mutate(
    text  = str_to_lower(paste(title, abstract, sep = " ")),
    tier  = case_when(
      journal_name %in% top5_journals ~ "Top-5",
      journal_name %in% aej_journals  ~ "AEJ field",
      TRUE                            ~ "Other"
    )
  )

## Method dictionaries. SDID is most specific, so we test it first and exclude
## its hits from the DiD bucket; SCM and DiD overlap freely. Patterns allow
## hyphen-or-space separators and singular/plural forms ("difference[s] in
## difference[s]", "synthetic[-]control[s]"). We deliberately do not search
## for the bare abbreviation "DiD" because of false positives ("Did the ...").
re_did  <- "differences?[-\\s]in[-\\s]differences?|diff[-\\s]?in[-\\s]?diff"
re_scm  <- "synthetic[-\\s]controls?(?:\\s+method)?"
re_sdid <- "synthetic\\s+differences?[-\\s]in[-\\s]differences?"

papers <- papers |>
  mutate(
    has_sdid = str_detect(text, re_sdid),
    has_scm  = str_detect(text, re_scm),
    has_did  = str_detect(text, re_did) & !has_sdid  # don't double-count SDID as DiD
  )

counts <- papers |>
  filter(tier %in% c("Top-5", "AEJ field")) |>
  group_by(publication_year, tier) |>
  summarise(
    n_papers = dplyr::n(),
    n_did    = sum(has_did),
    n_scm    = sum(has_scm),
    n_sdid   = sum(has_sdid),
    .groups  = "drop"
  ) |>
  pivot_longer(
    cols      = c(n_did, n_scm, n_sdid),
    names_to  = "method",
    values_to = "n_method"
  ) |>
  mutate(method = recode(method,
                         n_did  = "DiD",
                         n_scm  = "SCM",
                         n_sdid = "SDID"))

write_csv(counts, OUT_CSV)
message("Wrote ", OUT_CSV, " (", nrow(counts), " rows)")
