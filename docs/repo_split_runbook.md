# Repo split runbook (step 3) — 2026-07-10

Decision (settled): **current repo becomes the package** (`wsynthdid`), the
**paper moves to a new repo** that John creates on GitHub. `gh` here is authed
as `johniselin-budget-lab`, but `origin` is SSH under `johniselin-econ`, so the
remote repo must be created by John; pushing over SSH then works.

Recommended names: package R name **`wsynthdid`** (GitHub repo stays at
`johniselin-econ/synthdid_weights` until John renames it — a 1-click GitHub
rename with redirects, which I can't do via the budget-lab token); paper repo
**`weighted-sdid`** (adjust freely).

## File partition

**PACKAGE (stays in current repo):**
`R/`, `man/`, `tests/`, `vignettes/{synthdid.Rmd, synthdid_weighted.Rmd,
more-plotting.Rmd, output/, all-simulations.rds}`, `data/`, `DESCRIPTION`,
`NAMESPACE`, `LICENSE`, `_pkgdown.yml`, `.Rbuildignore`, `DEVELOPING.md`,
`synthdid.Rproj`, `.gitignore`, `README.md` (rewrite: package-only),
`CLAUDE.md` (rewrite: package-only).

**PAPER (moves to new repo):**
`paper/`, `scripts/`, `results/`, `results_2005/`, `docs/` (tracked planning
docs), `renv/`, `renv.lock`, `.Rprofile`, `LICENSE`, `.gitignore`,
`vignettes/paper-results.Rmd` (paper-specific vignette moves here), plus a new
paper `README.md` and paper-focused `CLAUDE.md`.

## Edits on each side

**Package (`DESCRIPTION`):** `Package: synthdid` -> `Package: wsynthdid`;
drop the "Currently provides methods only for the case that all treated units
adopt treatment at the same time" sentence (staggered now supported); point
`URL`/`BugReports` at the package repo. Change `library(synthdid)` ->
`library(wsynthdid)` in vignettes and any script/test that names it. Keep the
`upstream` remote (this is the fork).

**Paper:** manuscript footnote + data-availability URLs — package pointer ->
`wsynthdid` repo, replication pointer -> the paper repo (both currently point
to `synthdid_weights`). `scripts/*` `library(synthdid)` -> `library(wsynthdid)`.
Pin `wsynthdid` in `renv.lock` by git SHA (AFTER the package is pushed). Drop
the `upstream` remote.

## Execution order (remote steps gated on John)

1. [reversible, done now] Build the paper repo as an isolated clone
   `../weighted-sdid`, partitioned to the PAPER set above. Current repo and all
   remotes untouched. Undo = `rm -rf ../weighted-sdid`.
2. [John] Create empty `johniselin-econ/weighted-sdid` on GitHub.
3. Point the clone's `origin` at it, push. Verify the paper repo renders /
   scripts reference `wsynthdid`.
4. [John] Rename `johniselin-econ/synthdid_weights` -> `wsynthdid` on GitHub
   (optional; redirects preserve the old URL).
5. In the current repo: strip the PAPER set, apply the package `DESCRIPTION`
   edits, rewrite package README/CLAUDE, push. (Irreversible on the public repo,
   but recoverable from history + the paper repo.)
6. In the paper repo: `renv::install` the pushed `wsynthdid` at its SHA,
   `renv::snapshot()`, commit the lockfile, push.
7. Re-knit paper + supplement against the installed `wsynthdid`; run
   `R CMD check` on the package. (Both blocked locally by the intermittent R
   segfaults on this box — do on a healthy environment / the cluster.)

## Caveats

- R is segfaulting intermittently here; the package build (`R CMD check`) and
  the renv pin should be verified on a healthy environment before relying on the
  split.
- The renv SHA pin is chicken-and-egg: it requires the package pushed first
  (step 5 before step 6).
- `results/`, `results_2005/`, and `paper/data/` are large; they live only in
  the paper repo after the split.
