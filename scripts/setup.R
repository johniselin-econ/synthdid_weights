#!/usr/bin/env Rscript
# One-command environment bootstrap for reproducing the paper.
# Run from the repository root:
#   Rscript scripts/setup.R
#
# It restores the exact R package versions pinned in renv.lock, installs the
# local synthdid package into the project library, and checks for the system
# tools (pandoc, LaTeX) that renv cannot manage. After this succeeds, knit with
#   Rscript paper/_build_all.R

# The project .Rprofile sources renv/activate.R, which bootstraps renv itself
# on a fresh machine. Guard in case this script is run with --vanilla.
if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv", repos = "https://cloud.r-project.org")
}

cat("==== Restoring R packages from renv.lock ====\n")
renv::restore(prompt = FALSE)

cat("\n==== Installing local synthdid package ====\n")
renv::install(".")

# --- System-tool checks (not R packages; renv cannot provide these) ----------
cat("\n==== Checking system tools ====\n")

pandoc_ok <- nzchar(Sys.which("pandoc")) ||
  nzchar(Sys.getenv("RSTUDIO_PANDOC")) ||
  nzchar(Sys.which("quarto"))
if (pandoc_ok) {
  cat("pandoc: found\n")
} else {
  cat("pandoc: NOT found.\n",
      "  rmarkdown::render() needs pandoc. Options:\n",
      "  * Install pandoc (https://pandoc.org/installing.html) or quarto, OR\n",
      "  * set RSTUDIO_PANDOC to a directory containing the pandoc binary.\n",
      "  On a cluster without system pandoc and no root, download a static\n",
      "  pandoc binary into ~/bin and put ~/bin on your PATH.\n", sep = "")
}

latex_ok <- nzchar(Sys.which("pdflatex")) ||
  (requireNamespace("tinytex", quietly = TRUE) && tinytex::is_tinytex())
if (latex_ok) {
  cat("LaTeX: found\n")
} else {
  cat("LaTeX: NOT found. paper/_build_all.R will bootstrap TinyTeX on first run,\n",
      "  or install one manually with tinytex::install_tinytex().\n", sep = "")
}

cat("\n==== Setup complete ====\n")
if (pandoc_ok) {
  cat("Build the paper with: Rscript paper/_build_all.R\n")
} else {
  cat("Resolve pandoc (above), then build with: Rscript paper/_build_all.R\n")
}
