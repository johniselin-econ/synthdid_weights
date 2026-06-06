# Pre-submission review checklist (2026-06-06)

From a full read of both built PDFs (paper 19pp, supplement 22pp). Ordered by severity.

## Blockers — fix before anyone else sees it

1. **Wrong appendix letter, 4×.** The main paper cites "Appendix F" for the BV
   reproduction (Table 3 caption, Sec. 3.5, data statement ×2), but in the
   supplement it is **Section E.4**; F is "Alternative Approaches."
2. **Ghost journal reference.** Supplement D.2 runtime note refers to "the
   companion *Economics Letters* supplement (weighted_sdid_supplement.pdf,
   Tables E.1 and E.2 in that document)" — a self-referential leftover from
   the old dual-version setup, names the wrong journal, and says "for this
   draft we omit the rendered tables." Decide: render the MC tables or
   rewrite the note. Kill "for this draft" language everywhere.
3. **Abstract is ~370 words.** REStat expectation is ~100–150. Cut by two
   thirds; lead with the finding (−17.5 vs ≈0), not the method genealogy.
4. **highlights.txt and README still say the old story** ("estimates flip
   sign…"; README says "approximately 13.5"). Neither is a sign flip anymore
   (−0.71 vs −17.45, both negative). Update both before pushing anywhere
   public-facing.
5. **Numbered paragraph-heading artifacts**: "D.1.0.2", "D.2.0.1", "F.0.0.1"
   render as deep section numbers. Set `secnumdepth` / use unnumbered
   `\paragraph{}` styling.

## Consistency — things a referee will cross-check

6. **Small-county story across Figures 3 and 4.** Sec. 3.3 says bottom
   quartiles "oscillate around zero (−6 to +15)"; Sec. 3.3's urban figure
   says the most-rural quintile is **+21**; the Discussion says "near zero or
   positive." Pick one framing and explain the +21 (imputation noise in tiny
   counties?) or soften it.
7. **"Sign change / flip" sweep**: grep both Rmds + README + highlights for
   residual flip language (intro and abstract are fixed; verify nothing else).
8. **MC spec**: D reports 18 configs × **100** sims; confirm no stray claim
   of 500 sims, and decide whether 100 is the final spec (if so, the 50-hour
   runtime note can shrink to a sentence).
9. **Old coarse-review punch list** (docs/coarse_review_2026-04-17.md, items
   M1–M3 / A1–A7 / D1–D11 + docs/todo.txt): confirm each is addressed or
   consciously dropped. The MC-evidence gap (M-level) is now addressed.
10. **If any data re-pull happens**, re-run scripts/bv_replication_tables.R
    and re-knit BOTH docs (knitr caches don't track data changes). Appendix
    E.4 numbers were verified identical to the script on 2026-06-06.

## Submission hygiene — decide once

11. **WONDER DUA**: the public repo's analysis_data.csv now contains imputed
    1–9 death counts. They're model estimates, not actual counts, but
    consider shipping crude_rate only (drop the deaths column) to be safe.
12. **AI-use disclosure**: substantial AI assistance in code, replication, and
    draft editing this cycle — REStat/AEA policies want this disclosed at
    submission (not authorship; a disclosure line).
13. **Setup chunks still `source("../R/*.R")`** instead of
    `library(synthdid)`. Works, but the replication README must say so, or
    switch now that the fork installs cleanly.
14. **Affiliations/thanks**: header.tex hardcodes "Erica Ryan (Amazon)" — confirm
    she wants the affiliation + add any disclaimer her employer requires.
15. Cosmetics: float-specifier warnings (`!h`→`!ht`) are benign; lit-counts
    figure caption is very long (consider moving the OpenAlex/method detail
    to a footnote); check title page date renders as intended.

Verified clean: zero undefined references in both PDFs; supplement E.4 table
matches the standalone replication script exactly; LOO, placebo, event-study
inline numbers all recomputed on the corrected panel; political-rows footnote
renders; eq. (1) reference resolves.
