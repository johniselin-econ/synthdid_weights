# Pre-Submission Review: "Weighted Synthetic Difference-in-Differences"

**Authors:** John Iselin and Erica Ryan
**Target:** Economics Letters
**Date:** April 15, 2026
**Reviewer:** Automated pre-submission review (Claude)

---

## 1. Spelling, Grammar & Style

### Spelling/Typos
### Grammar & Phrasing
- **Line 232 (main):** "The SDID pre-treatment fit is slightly tighter than DiD" -- "slightly tighter" is imprecise; consider quantifying (e.g., RMSPE of pre-treatment residuals).

### Style
### Typographic Consistency

---

## 2. Internal Consistency & Cross-References

### Numerical Consistency
### Cross-Reference Correctness
### Terminology Consistency
### Citation Verification
### Sample Description

---

## 3. Unsupported Claims & Identification Integrity

### Causal Language
### Generalization Beyond Sample
### Missing Caveats
### Statistical vs. Economic Significance


---

## 4. Mathematics, Equations & Notation

### Equation Correctness
### Notation Consistency
### Undefined Symbols
### LaTeX Formatting

---

## 5. Tables, Figures & Documentation

### Table 1
- **Missing from table:** No asterisks or significance indicators. This is actually fine for a letter focused on the sign change rather than significance testing, but consider adding them.

### Figure 1

### Data Availability Statement (line 253-256)
- GitHub repo listed. Good.
- **Missing:** No mention of which specific CDC WONDER query parameters (age range, cause of death, etc.) were used. The paper says "all-cause mortality rate per 100,000 among adults aged 20--64" but the data availability statement does not repeat this specificity.

### Replication Readiness
- `cache = TRUE` is used on computationally intensive chunks. Good for development but should be cleared and re-run from scratch before final submission.
- **RECOMMENDATION:** Include a `Makefile` or script that knits both the main paper and supplement cleanly from scratch.

---

## 6. Contribution Evaluation (Leading Field Journal Standard)

### Required Analyses (for acceptance)
1. **Clustered standard errors:** The paper acknowledges this limitation but does not provide them. State-clustered SEs (or at minimum wild cluster bootstrap p-values) should be reported alongside county-level SEs, even if only in a robustness table. The main result could be overturned if state-clustered SEs are much larger.
2. **Robustness to weight specification:** Show results with alternative weight specifications (e.g., log population, population shares from a different base year, employment-based weights). This tests whether the sign change is specific to 2013 population weights.
3. **SC estimates:** The main text omits SC results (onlyDIDand SDID). The supplement includes them. At least report SC in the main table to complete the comparison.
4. **Formal significance testing of the divergence:** The paper claims the weighting "determines the sign" but does not formally test whether the equally-weighted and population-weighted estimates are statistically distinguishable from each other. A Hausman-type test or joint confidence region would strengthen the claim.

### Suggested Analyses (to strengthen)
1. **Decomposition by county size:** Show a binned scatter or quantile analysis of unit-level treatment effects vs. population. This would directly demonstrate the heterogeneity-size correlation that drives the result.
2. **Leave-one-out sensitivity:** What happens when the largest county (likely LA County or Cook County) is dropped? If the result is fragile to a single large county, this undermines the generalizability.
3. **Comparison to Borgschulte & Vogler (2020) headline numbers:** The supplement implements a DLPS benchmark but with fewer covariates. Directly compare to the published estimates from BV2020.
4. **Pre-treatment balance table:** Show covariate means for treated vs. (weighted) control under both weighting schemes.
5. **Placebo treatment dates:** Run the estimator with fake treatment dates (e.g., 2010, 2011) as a falsification test.

### Literature Positioning
### Journal Fit
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

### MAJOR (strongly recommended before submission)

8. **Add SC estimates to the main Table 1.** Three estimators (DiD, SC, SDID) x two weighting schemes completes the comparison.

### MINOR (polish items)

17. Add missing literature references (Ben-Michael et al. 2021; Gibbons, Suarez Serrato & Urbancic 2019 on implicit weights).
18. Clear all Rmd caches and verify a clean knit from scratch before submission.

---

## Preliminary Recommendation

**Revise and resubmit (minor-to-major).** The paper makes a clear, useful contribution well-suited to Economics Letters. The sign-change application is striking and the implementation is clean. However, the absence of state-clustered SEs is a near-certain rejection trigger for any referee familiar with the ACA literature. The broken supplement cross-references are a presentation-quality issue that must be fixed. With the critical and major items addressed, this paper has a realistic path to acceptance at Economics Letters.






REMAINING:

Pre-treatment balance table — new table
Cache clearing / clean knit — pre-submission, verify clean knit from scratch
