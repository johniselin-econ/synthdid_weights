#!/usr/bin/env Rscript
# Knit main paper then supplement from clean state.
# Usage: Rscript _build_all.R

setwd("C:/Users/ji252/Documents/GitHub/synthdid_weights/paper")
stopifnot(file.exists("weighted_sdid_paper.Rmd"),
          file.exists("weighted_sdid_supplement.Rmd"))

# Point rmarkdown at RStudio's bundled pandoc (no system-wide pandoc installed)
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

t_start <- Sys.time()
ok_paper      <- knit_one("weighted_sdid_paper.Rmd")
ok_supplement <- knit_one("weighted_sdid_supplement.Rmd")
total <- difftime(Sys.time(), t_start, units = "mins")

cat(sprintf("\n==== Build summary: paper=%s, supplement=%s, total=%.1f min ====\n",
            if (ok_paper) "OK" else "FAIL",
            if (ok_supplement) "OK" else "FAIL",
            as.numeric(total)))

if (!ok_paper || !ok_supplement) quit(status = 1)
