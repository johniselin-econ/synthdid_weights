# Pre-Submission Review: "Weighted Synthetic Difference-in-Differences"

**Authors:** John Iselin and Erica Ryan
**Target:** Economics Letters
**Date:** April 15, 2026
**Reviewer:** Automated pre-submission review (Claude)

---

## 1. Spelling, Grammar & Style

### Spelling/Typos
- No spelling errors detected in the main text or supplement.

### Grammar & Phrasing
- **Line 61 (main):** "This is an implicit design choice, not a requirement of the methodology." -- Clean and punchy. No issue.
- **Line 63 (main):** Long sentence ("The modification is minimal---one equation change---but it affects..."). Consider splitting for clarity. The parenthetical "---one equation change---" is effective but the sentence runs to 4 clauses.
- **Line 94 (main):** "Weighted DiD ($\hat{\omega} = 1/N_0$, $\hat{\lambda} = 1/T_0$) and weighted SC ($\hat{\lambda} = 0$) emerge as special cases." -- These should be $\hat{\omega} = \mathbf{1}/N_0$ (vector) and $\hat{\lambda} = \mathbf{1}/T_0$ (vector), not scalars. See Section 4 below.
- **Line 199 (main):** "consistent with greater baseline uninsurance rates and more scope for coverage gains in urban areas" -- This is a mechanism claim stated without citation or evidence. See Section 3 below.
- **Line 232 (main):** "The SDID pre-treatment fit is slightly tighter than DiD" -- "slightly tighter" is imprecise; consider quantifying (e.g., RMSPE of pre-treatment residuals).
- **Line 240 (main):** "Standard errors are modestly larger under population weighting (approximately 2.4 vs. 2.0 for SDID)" -- These specific numbers are hardcoded prose but should be computed inline (or verified against the computed results). If the code output changes, these will be stale.

### Style
- **Passive voice:** Generally well-controlled. A few instances ("We have proposed," line 247) use active voice appropriately.
- **Tense consistency:** Present tense used for methodology, past/present for results. Consistent throughout.
- **Filler words:** None detected. The writing is admirably tight for a short letter.

### Typographic Consistency
- Em-dashes (---) used consistently throughout.
- Numerical formatting: "13.5 deaths per 100,000" in abstract matches body text. Good.
- "per-capita" (line 63) -- should be "per capita" (no hyphen when used as a standalone adjective after the noun, but acceptable as compound modifier before noun). Since it precedes "treatment effect," the hyphenated form is defensible.

---

## 2. Internal Consistency & Cross-References

### Numerical Consistency
- **Abstract:** "near zero" and "approximately 13.5 deaths per 100,000" -- matches the body text (line 199).
- **Abstract vs. body:** The abstract says "reduction of approximately 13.5 deaths per 100,000." The body says "approximately $-13.5$ deaths per 100,000." These are consistent.
- **Hardcoded SE values (line 240):** "approximately 2.4 vs. 2.0 for SDID" -- These are not dynamically computed from the R code. If the data or code changes, these could become inconsistent with Table 1. **RECOMMENDATION:** Replace with inline R expressions.

### Cross-Reference Correctness
- **Main text:** References to Table 1 (line 199) and Figure 1 (line 232) are correct.
- **Main text line 94:** "proof in Supplementary Appendix C" -- The proof IS in Appendix C of the supplement. Correct.
- **Main text line 102:** "Supplementary Appendix B for details" -- Appendix B covers inference. Correct.
- **Main text line 240:** "Supplementary Appendix F" -- Appendix F covers the DLPS benchmark. Correct.
- **Main text line 242:** "Supplementary Appendix G" -- Appendix G covers extensions including clustered SEs. Correct.
- **Main text line 247:** "Supplementary Appendix E" -- Appendix E covers Monte Carlo simulations. Correct.
- **CRITICAL -- Supplement broken \eqref references:** The proof in Appendix C (line 152) references `\eqref{eq:weighted-collapsed}` and `\eqref{eq:tau-weighted}`, and Appendix G (line 717) references `\eqref{eq:weighted-collapsed}` and `\eqref{eq:lambda-w}`. **None of these labels are defined anywhere in the supplement.** The supplement defines `eq:collapsed`, `eq:sdid-regression`, `eq:tau-sdid`, `eq:omega`, `eq:lambda`, `eq:boot-weights`, `eq:jk-weights`, but NOT `eq:weighted-collapsed`, `eq:tau-weighted`, or `eq:lambda-w`. These will render as "(??" in the PDF. **FIX REQUIRED:** Add labeled weighted equations to the supplement, or change the references to point to the main text equations `eq:collapsed-w` and `eq:tau-w`.
- **Supplement Table E.1 (line 517):** HHI "= 0.10" mentioned in text but the configurations use HHI = hhi_mult/N1, i.e., {1.5, 3}/{5, 10, 20}. For N1=10, high concentration is 3/10 = 0.30. For N1=20, high is 3/20 = 0.15. The value 0.10 does not appear as a configuration. **Possible error:** should this be "HHI $= 3/N_1$" or "HHI $\leq 0.30$"?

### Terminology Consistency
- "Equally weighted" vs. "unweighted" -- The paper uses both interchangeably (e.g., Table 1 header says "Equally weighted" but body text says "unweighted SDID estimate"). **RECOMMENDATION:** Standardize. "Equally weighted" is more precise since $1/N_1$ IS a weighting scheme.
- "Population weighted" vs. "population-weighted" -- Hyphenation is inconsistent. Used as "population-weighted" (compound adjective, correct) in most places, but the table header and legend say "Population weighted" (no hyphen, as two-word label). Acceptable in table/figure labels.

### Citation Verification
- **Cited in main text:** `arkhangelsky2021`, `borgschulte2020`, `ciccia2024`, `belloni2014`. All present in references.bib.
- **Cited in supplement:** `arkhangelsky2021`, `borgschulte2020`, `belloni2014`, `ciccia2024`, `callaway2021`, `sun2021`, `dechaisemartin2020`. All present in references.bib.
- **Uncited entries in references.bib:** `sommers2017`, `miller2021`, `abadie2010`, `abadie2015`. These are in the .bib but never cited in either the main text or supplement. They will NOT appear in the bibliography (natbib only includes cited references), so this is harmless but untidy. If they are there for future use, leave them; otherwise clean up.
- **Sommers et al. citation:** The bib entry `sommers2017` has `year={2016}` but key says `2017`. Minor bib housekeeping issue (won't affect output since uncited).

### Sample Description
- Panel dimensions are computed dynamically from data (lines 121-126), so they will be internally consistent. Good practice.
- Treatment year 2014 is hardcoded in `mutate(post = as.integer(time >= 2014))` (line 119) and `T0 <- sum(sort(unique(panel$time)) < 2014)` (line 126). Consistent.

---

## 3. Unsupported Claims & Identification Integrity

### Causal Language
- **Line 199:** "Larger counties experience larger mortality reductions from Medicaid expansion, consistent with greater baseline uninsurance rates and more scope for coverage gains in urban areas." This is a causal mechanism claim. The paper provides no direct evidence on uninsurance rates or coverage gains. **RECOMMENDATION:** Soften to "which may reflect greater baseline uninsurance rates..." or cite supporting evidence (e.g., Sommers et al. 2016, Miller et al. 2021 -- both are in the .bib but uncited).
- **Line 65:** "a sign change driven entirely by the weighting choice" -- The word "entirely" is strong. This is appropriate since the same data and same estimator are used, with only the weights changing. Acceptable.

### Generalization Beyond Sample
- **Line 247:** "a finding that generalizes to any setting where treatment effects are heterogeneous and correlated with unit size" -- This is an overly broad claim. The theoretical result (weighted estimator targets a different estimand) is general, but the empirical finding of a sign change is specific to this application. **RECOMMENDATION:** Soften to "a finding relevant to any setting..." or "an issue that arises whenever..."

### Missing Caveats
- **State-level clustering (line 241-242):** The paper correctly flags that county-level resampling does not account for within-state correlation. This is a significant limitation for an application where treatment is assigned at the state level. The paper acknowledges it but then uses county-level SEs in the main results without stating how this might bias inference. **RECOMMENDATION:** Add a sentence noting the direction of likely bias (county-level SEs are likely too small, so the true confidence intervals are wider).
- **No covariates in main specification:** The main SDID estimates use no covariates (X = empty). This is standard for SDID (which handles level differences via the intercept), but for a mortality application with large county heterogeneity, it is worth noting. The DLPS benchmark in the supplement does use covariates, providing some reassurance.
- **Borgschulte and Vogler (2020) comparison:** The main text (line 240) references this but the supplement (Appendix F) notes the covariate set is "smaller than the one used by Borgschulte and Vogler." The main text does not flag this limitation.

### Statistical vs. Economic Significance
- The paper appropriately focuses on the sign change (qualitative difference) rather than precision. However, the hardcoded SEs "approximately 2.4 vs. 2.0" suggest both estimates are imprecise relative to their magnitudes. The unweighted estimate is "near zero" with SE ~2.0, so it is clearly insignificant. The weighted estimate at ~-13.5 with SE ~2.4 is significant at the 5% level. This distinction is implicit but never stated explicitly. **RECOMMENDATION:** State significance explicitly in the results discussion.

---

## 4. Mathematics, Equations & Notation

### Equation Correctness
- **Equation (1) (line 73-76):** Collapsed form. The indexing in the sum for the last row uses $j=1$ to $N_1$ with $Y_{N_0+j,t}$. Correct. The last column uses $1/T_1 \sum_{t=T_0+1}^{T} Y_{it}$. Correct.
- **Equation (2) (line 78-80):** The bilinear form for the SDID estimate. The left vector has $-\hat{\omega}$ (length $N_0$) and $\frac{1}{N_1}\mathbf{1}_{N_1}$ (length $N_1$), giving a vector of length $N$. The right vector has $-\hat{\lambda}$ (length $T_0$) and $\frac{1}{T_1}\mathbf{1}_{T_1}$ (length $T_1$), giving a vector of length $T$. The product $v^\top Y w$ is well-defined. Correct.
- **Equation (3) (line 86-89):** Weighted collapsed form. Only the last ROW uses $\tilde{\omega}_j$; the last COLUMN still uses $1/T_1$. This is consistent with the paper's claim that period weights are left uniform in the main text (the period weight extension is in the supplement). Correct.
- **Equation (4) (line 91-93):** Weighted bilinear form. The treated block uses $\tilde{\omega}$ instead of $\frac{1}{N_1}\mathbf{1}_{N_1}$. Correct. The time vector is unchanged. Correct.
- **Equation (5) (line 99-101):** Bootstrap weight renormalization. The formula correctly multiplies original weights by bootstrap multiplicities and renormalizes. Correct.

### Notation Consistency
- $\hat{\omega}$ for control unit weights, $\tilde{\omega}$ for researcher-specified treated weights. This convention is used consistently throughout.
- $\hat{\omega}^w$ for control weights estimated from the weighted collapsed form. Consistent.
- $\hat{\lambda}^w$ for time weights under the weighted estimator. Consistent.
- $\Delta^{N_0}$ for the unit simplex. Standard notation. Consistent.
- **ISSUE (line 94):** "Weighted DiD ($\hat{\omega} = 1/N_0$, $\hat{\lambda} = 1/T_0$)" -- Here $1/N_0$ and $1/T_0$ should be vectors: $\hat{\omega} = (1/N_0)\mathbf{1}_{N_0}$ and $\hat{\lambda} = (1/T_0)\mathbf{1}_{T_0}$. As written, they appear scalar. Similarly, "weighted SC ($\hat{\lambda} = 0$)" should be $\hat{\lambda} = \mathbf{0}_{T_0}$. This is a minor but noticeable imprecision for a methods paper.

### Undefined Symbols
- $\tilde{\pi}_s$ appears in the supplement Proposition 1 (line 149) as period weights but is NOT defined in the main text (which uses $1/T_1$ uniformly). The supplement should define $\tilde{\pi}$ before using it in Proposition 1, or the proposition should be stated using only $\tilde{\omega}$ (since the main text does not introduce period weights).
- $\bar{\tau}_j$ (line 85): defined as "the average treatment effect for unit $j$." Slightly ambiguous -- average over what? Presumably over post-treatment periods: $\bar{\tau}_j = T_1^{-1} \sum_{t=T_0+1}^{T} \tau_{jt}$. Worth making explicit.

### LaTeX Formatting
- All equation labels use `\label{eq:...}` and are referenced via `\eqref{eq:...}`. Standard and correct in the main text.
- The supplement has broken `\eqref` references (see Section 2 above).
- The `\ind` command (defined in header.tex as `\mathbf{1}`) is not used in the text (which uses `\mathbf{1}` directly). This is fine but the custom command is wasted.

---

## 5. Tables, Figures & Documentation

### Table 1
- **Caption (line 192):** "Effect of ACA Medicaid expansion on mortality (deaths per 100,000, ages 20--64). Jackknife standard errors." Clear and complete.
- **Column headers:** "Equally weighted" and "Population weighted" spanning headers with "Estimate" and "SE" sub-columns. Well-structured.
- **Notes:** No explicit table notes beyond the caption. **RECOMMENDATION:** Add a note specifying the sample (counties, years), outcome definition (CDC WONDER), and treatment definition (state expanded Medicaid by Jan 2014). This is especially important for a short letter where the data description is brief.
- **Digits:** Set to 2 decimal places. Appropriate for deaths per 100,000.
- **Missing from table:** No asterisks or significance indicators. This is actually fine for a letter focused on the sign change rather than significance testing, but consider adding them.

### Figure 1
- **Caption (line 201):** "Event-study decomposition under equal and population weighting. Left panel: DiD. Right panel: SDID. Shaded bands show 95% bootstrap confidence intervals (200 replications). The dashed vertical line marks treatment onset (2014)." Thorough.
- **Figure dimensions:** 8 x 4 inches, faceted by estimator. Appropriate for a two-panel figure.
- **Color scheme:** steelblue (equally weighted) vs. firebrick (population weighted). Colorblind-accessible? Blue/red is generally distinguishable for most forms of color blindness, and the line styles differ (both solid but with different colors + CI bands). **RECOMMENDATION:** Consider adding different point shapes for additional accessibility.
- **Y-axis label:** "Effect on mortality\n(deaths per 100,000)" -- Good. The line break is appropriate for the figure.
- **X-axis:** Unlabeled (NULL), relying on year tick marks. Acceptable for a time-series event study.
- **Bootstrap replications:** 200. This is on the low side for publication but acceptable for a short letter. Consider noting this choice in the text or increasing to 500.

### Data Availability Statement (line 253-256)
- Sources listed: CDC WONDER, SEER, BLS. Correct per the data loading code (crude_rate, population, unemployment).
- GitHub repo listed. Good.
- **Missing:** No mention of which specific CDC WONDER query parameters (age range, cause of death, etc.) were used. The paper says "all-cause mortality rate per 100,000 among adults aged 20--64" but the data availability statement does not repeat this specificity.

### Replication Readiness
- The Rmd file sources from `../R/` (lines 48-55), which contains the package code. The code is self-contained with a seed (line 57).
- Data is loaded from `../data/analysis_data.csv` (line 111). This file is presumably in the repo but is not verified here.
- `cache = TRUE` is used on computationally intensive chunks. Good for development but should be cleared and re-run from scratch before final submission.
- **RECOMMENDATION:** Include a `Makefile` or script that knits both the main paper and supplement cleanly from scratch.

---

## 6. Contribution Evaluation (Leading Field Journal Standard)

### Central Contribution Rating: **Significant (low end) / Incremental (high end)**

The paper makes a clear, well-executed methodological contribution: extending SDID to handle non-uniform treated-unit weights. The modification is simple (one equation change to the collapsed form), which is both a strength (easy to implement) and a weakness (may appear too simple for a standalone paper). The application is compelling: a sign change in the estimated treatment effect demonstrates that the extension is practically consequential.

However:
- The weighted extension is conceptually straightforward -- it amounts to recognizing that the collapsed form uses an implicit weighting scheme and allowing the researcher to change it.
- No new asymptotic theory is provided; the paper adapts existing variance estimators.
- The Monte Carlo evidence (in the supplement) supports the claim but does not reveal surprising properties.

For Economics Letters specifically, this is a good fit: the contribution is focused, the paper is short, and the result is striking. The bar for novelty is lower than AER/QJE, and "short notes" on consequential methodological modifications are exactly what EL publishes.

### Identification Credibility
- The SDID identification relies on parallel trends (conditional on the synthetic control weights), which is standard.
- The event-study figure shows pre-treatment coefficients near zero under both weighting schemes, supporting the identifying assumption.
- The state-level treatment / county-level observation mismatch is a real concern for inference (acknowledged in the paper).
- The application is well-chosen: the ACA Medicaid expansion is a canonical natural experiment with a plausible pre-treatment parallel-trends story.

### Required Analyses (for acceptance)
1. **Clustered standard errors:** The paper acknowledges this limitation but does not provide them. State-clustered SEs (or at minimum wild cluster bootstrap p-values) should be reported alongside county-level SEs, even if only in a robustness table. The main result could be overturned if state-clustered SEs are much larger.
2. **Robustness to weight specification:** Show results with alternative weight specifications (e.g., log population, population shares from a different base year, employment-based weights). This tests whether the sign change is specific to 2013 population weights.
3. **SC estimates:** The main text omits SC results (only DiD and SDID). The supplement includes them. At least report SC in the main table to complete the comparison.
4. **Formal significance testing of the divergence:** The paper claims the weighting "determines the sign" but does not formally test whether the equally-weighted and population-weighted estimates are statistically distinguishable from each other. A Hausman-type test or joint confidence region would strengthen the claim.

### Suggested Analyses (to strengthen)
1. **Decomposition by county size:** Show a binned scatter or quantile analysis of unit-level treatment effects vs. population. This would directly demonstrate the heterogeneity-size correlation that drives the result.
2. **Leave-one-out sensitivity:** What happens when the largest county (likely LA County or Cook County) is dropped? If the result is fragile to a single large county, this undermines the generalizability.
3. **Comparison to Borgschulte & Vogler (2020) headline numbers:** The supplement implements a DLPS benchmark but with fewer covariates. Directly compare to the published estimates from BV2020.
4. **Pre-treatment balance table:** Show covariate means for treated vs. (weighted) control under both weighting schemes.
5. **Placebo treatment dates:** Run the estimator with fake treatment dates (e.g., 2010, 2011) as a falsification test.

### Literature Positioning
- The paper correctly positions against Arkhangelsky et al. (2021) as the base method.
- The staggered-adoption literature (Callaway & Sant'Anna, de Chaisemartin & D'Haultfoeuille, Sun & Abraham) is cited in the supplement for the staggered extension discussion.
- **Missing:** Ben-Michael, Feller & Rothstein (2021, JASA) on augmented SCM with covariates; Cattaneo, Feng & Titiunik (2021) on prediction intervals for SC; Pang, Liu, Xu & Li (2022) on factor-model-based approaches. These are relevant competitors/complements for the weighted SC literature.
- **Missing:** The "aggregation weights" discussion in the heterogeneous-TE literature (e.g., Gibbons, Suarez Serrato & Urbancic 2019 on implicit regression weights; Sloczynski 2022 on OLS interpretation) is conceptually related and could strengthen the motivation.

### Journal Fit
Economics Letters is a strong fit. The paper is:
- Short (under 2,000 words including code output)
- Focused on a single, clean methodological point
- Accompanied by a compelling application
- Accompanied by a thorough supplement with proofs, simulations, and extensions

The main risk is a referee who views the contribution as "too obvious" (just change the weights in the collapsed form). The sign-change application is the key defense against this critique.

### Referee Questions (4-7 pointed questions a referee might ask)

1. **"Why not just estimate heterogeneous treatment effects directly?"** If the concern is that effects vary with county size, one could estimate $\tau_j$ for each county using unit-by-unit SDID (your Solution 2 in Appendix F) and then aggregate. How does this compare? Is the collapsed-form modification computationally necessary or merely convenient?

2. **"How sensitive is the sign change to the largest counties?"** If the top 10 counties (which account for a large share of weight) are dropped, does the population-weighted estimate remain negative? Is this a "Los Angeles County drives the result" story?

3. **"What are the state-clustered standard errors?"** You acknowledge this limitation but do not provide the numbers. Given that treatment is assigned at the state level, this is essential. If the state-clustered SE for the weighted estimate exceeds 7, the result is no longer significant at 5%.

4. **"The placebo variance estimator uses uniform weights for pseudo-treated units. How much does this matter?"** You note this discrepancy and offer alternatives in the software. Which one should practitioners use? The simulation evidence should directly compare the coverage of the uniform, size-match, and permute variants.

5. **"How does this relate to the choice of estimand in the heterogeneous treatment effects literature?"** Callaway & Sant'Anna (2021) emphasize that the aggregation scheme matters. Your contribution can be viewed as bringing this insight to the SDID framework. Is there a formal equivalence result?

6. **"Can you provide asymptotic theory?"** The paper adapts variance estimators without new asymptotics. Under what conditions does the weighted SDID converge to the population-weighted ATT? What is the asymptotic variance? Without this, the Monte Carlo evidence is the only basis for validity of inference.

7. **"Is the sign change real or an artifact of weight concentration?"** With an effective N1 of ~[computed value], the population-weighted estimate is driven by a small number of large counties. Could compositional differences between large and small counties (not treatment effect heterogeneity) explain the divergence?

---

## Priority Action Items

### CRITICAL (must fix before submission)

1. **Fix broken cross-references in the supplement.** `\eqref{eq:weighted-collapsed}`, `\eqref{eq:tau-weighted}`, and `\eqref{eq:lambda-w}` are undefined and will render as "??" in the PDF. Either add labeled equations for the weighted collapsed form and weighted estimator to the supplement, or redirect the references to the main-text labels `eq:collapsed-w` and `eq:tau-w`.

2. **Report state-clustered standard errors** (at minimum in a robustness table). The main results use county-level jackknife SEs, but treatment is at the state level. This is the single most predictable referee objection.

3. **Verify or dynamize the hardcoded SE values on line 240.** "approximately 2.4 vs. 2.0 for SDID" must match the computed output. Replace with inline R.

### MAJOR (strongly recommended before submission)

4. **Standardize "unweighted" vs. "equally weighted" terminology.** Pick one and use it consistently throughout main text, tables, and figures.

5. **Soften the unsupported mechanism claim (line 199)** about uninsurance rates and urban coverage gains. Either cite evidence or hedge with "may reflect."

6. **Soften the overgeneralization (line 247)** from "a finding that generalizes to any setting" to "an issue relevant to any setting."

7. **Fix the notation on line 94** -- DiD and SC special-case weights should be vector-valued, not scalar.

8. **Add SC estimates to the main Table 1.** Three estimators (DiD, SC, SDID) x two weighting schemes completes the comparison.

9. **Supplement Table E.1 text (line 517):** Verify the HHI = 0.10 claim against the actual simulation configurations. It appears to be incorrect.

10. **Define $\bar{\tau}_j$ explicitly** (line 85) as the average over post-treatment periods.

11. **Define $\tilde{\pi}_s$ in the supplement before Proposition 1** or remove it from the proposition statement (the main text does not introduce period weights).

### MINOR (polish items)

12. Clean up uncited .bib entries (`sommers2017`, `miller2021`, `abadie2010`, `abadie2015`) or cite them where appropriate (e.g., `sommers2017` and `miller2021` could support the mechanism claim on line 199).

13. Fix `sommers2017` bib entry: key says 2017 but `year={2016}`.

14. Add point shapes to Figure 1 for colorblind accessibility.

15. Add a brief table note to Table 1 specifying sample size, years, and outcome definition.

16. Consider increasing bootstrap replications from 200 to 500 for the published figure.

17. Add missing literature references (Ben-Michael et al. 2021; Gibbons, Suarez Serrato & Urbancic 2019 on implicit weights).

18. Clear all Rmd caches and verify a clean knit from scratch before submission.

---

## Preliminary Recommendation

**Revise and resubmit (minor-to-major).** The paper makes a clear, useful contribution well-suited to Economics Letters. The sign-change application is striking and the implementation is clean. However, the absence of state-clustered SEs is a near-certain rejection trigger for any referee familiar with the ACA literature. The broken supplement cross-references are a presentation-quality issue that must be fixed. With the critical and major items addressed, this paper has a realistic path to acceptance at Economics Letters.
