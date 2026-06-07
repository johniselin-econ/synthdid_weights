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

# Point rmarkdown at RStudio's bundled pandoc when no system pandoc exists
# (the Windows build machine has no system-wide pandoc)
pandoc_dir <- "C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools"
if (file.exists(file.path(pandoc_dir, "pandoc.exe"))) {
  Sys.setenv(RSTUDIO_PANDOC = pandoc_dir)
  rmarkdown::find_pandoc(dir = pandoc_dir)
}
cat(sprintf("pandoc: %s (version %s)\n",
            rmarkdown::pandoc_exec(),
            as.character(rmarkdown::pandoc_version())))

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
