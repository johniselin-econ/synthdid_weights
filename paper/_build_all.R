#!/usr/bin/env Rscript
# Knit the paper and/or supplement from a clean state.
# Usage:
#   Rscript _build_all.R              # knit both
#   Rscript _build_all.R paper        # knit only the paper
#   Rscript _build_all.R supplement   # knit only the supplement
# Run from the paper/ directory regardless of where the script is invoked
args0    <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", grep("^--file=", args0, value = TRUE))
if (length(file_arg) == 1) setwd(dirname(normalizePath(file_arg)))
stopifnot(file.exists("weighted_sdid_paper.Rmd"),
          file.exists("weighted_sdid_supplement.Rmd"))

# Locate pandoc on any OS. rmarkdown::render() needs a pandoc binary; it is a
# system tool, not an R package, so renv cannot provide it. We try, in order:
#   1. RSTUDIO_PANDOC, if already set to a valid directory
#   2. pandoc on the system PATH (e.g. a static binary dropped in ~/bin)
#   3. quarto's bundled pandoc tools directory
#   4. known bundled locations on Windows / macOS
# and fail with an actionable message if none is found.
locate_pandoc_dir <- function() {
  has_pandoc <- function(d) nzchar(d) &&
    (file.exists(file.path(d, "pandoc")) ||
     file.exists(file.path(d, "pandoc.exe")))
  # 1. honour an existing RSTUDIO_PANDOC
  env_dir <- Sys.getenv("RSTUDIO_PANDOC")
  if (has_pandoc(env_dir)) return(env_dir)
  # 2. pandoc on PATH
  on_path <- Sys.which("pandoc")
  if (nzchar(on_path)) return(dirname(on_path))
  # 3. quarto bundles pandoc under <quarto>/../tools (or tools/<arch>)
  quarto <- Sys.which("quarto")
  if (nzchar(quarto)) {
    cand <- file.path(dirname(dirname(quarto)), "tools")
    hits <- list.files(cand, pattern = "^pandoc(\\.exe)?$",
                       recursive = TRUE, full.names = TRUE)
    if (length(hits)) return(dirname(hits[1]))
  }
  # 4. known bundled locations
  candidates <- c(
    "C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools",
    "/Applications/RStudio.app/Contents/Resources/app/bin/quarto/bin/tools",
    "/usr/local/bin"
  )
  for (d in candidates) if (has_pandoc(d)) return(d)
  NA_character_
}

pandoc_dir <- locate_pandoc_dir()
if (is.na(pandoc_dir)) {
  stop(
    "No pandoc binary found. rmarkdown::render() cannot run without pandoc.\n",
    "Fix this one of these ways:\n",
    "  * Install pandoc (https://pandoc.org/installing.html) or quarto, OR\n",
    "  * Set RSTUDIO_PANDOC to a directory containing the pandoc binary.\n",
    "On a cluster with no system pandoc (and where you lack root), download a\n",
    "static pandoc binary into ~/bin and ensure ~/bin is on your PATH.\n",
    call. = FALSE
  )
}
Sys.setenv(RSTUDIO_PANDOC = pandoc_dir)
rmarkdown::find_pandoc(dir = pandoc_dir)
cat(sprintf("pandoc: %s (version %s)\n",
            rmarkdown::pandoc_exec(),
            as.character(rmarkdown::pandoc_version())))

# Ensure a LaTeX engine exists; bootstrap TinyTeX if neither system TeX nor
# an existing TinyTeX install is present.
if (!nzchar(Sys.which("pdflatex")) && !tinytex::is_tinytex()) {
  cat("No LaTeX engine found; installing TinyTeX (one-time)...\n")
  tinytex::install_tinytex()
}

knit_one <- function(rmd) {
  cat(sprintf("\n==== Knitting %s ====\n", rmd))
  t0 <- Sys.time()
  ok <- tryCatch({
    rmarkdown::render(rmd, quiet = FALSE, envir = new.env())
    TRUE
  }, error = function(e) {
    cat(sprintf("[ERROR] %s failed: %s\n", rmd, conditionMessage(e)))
    FALSE
  })
  dt <- difftime(Sys.time(), t0, units = "mins")
  cat(sprintf("==== %s %s in %.1f min ====\n",
              rmd, if (ok) "OK" else "FAILED", as.numeric(dt)))
  invisible(ok)
}

args <- commandArgs(trailingOnly = TRUE)
target <- if (length(args) == 0) "both" else args[1]
stopifnot(target %in% c("both", "paper", "supplement"))

t_start <- Sys.time()
ok_paper      <- if (target %in% c("both", "paper"))      knit_one("weighted_sdid_paper.Rmd")      else NA
ok_supplement <- if (target %in% c("both", "supplement")) knit_one("weighted_sdid_supplement.Rmd") else NA
total <- difftime(Sys.time(), t_start, units = "mins")

cat(sprintf("\n==== Build summary: paper=%s, supplement=%s, total=%.1f min ====\n",
            if (is.na(ok_paper))      "skipped" else if (ok_paper)      "OK" else "FAIL",
            if (is.na(ok_supplement)) "skipped" else if (ok_supplement) "OK" else "FAIL",
            as.numeric(total)))

if (isFALSE(ok_paper) || isFALSE(ok_supplement)) quit(status = 1)
