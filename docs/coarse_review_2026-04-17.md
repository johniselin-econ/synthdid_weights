# Weighted Synthetic Difference-in-Differences

**Date**: 04/17/2026
**Domain**: social_sciences/economics
**Taxonomy**: academic/working_paper
**Filter**: Active comments

---

## Overall Feedback

Here are some overall reactions to the document.

**Outline**

This paper proposes replacing the uniform 1/N₁ averaging of treated units in the SDID estimator of Arkhangelsky et al. (2021) with researcher-specified weights. The modification is applied to ACA Medicaid expansion using county-level mortality data, where equally-weighted and population-weighted estimates differ in sign. The idea is sound and practically relevant, but the paper as currently written lacks the theoretical and simulation infrastructure expected of a methods contribution targeting a top econometrics venue.

The paper identifies a genuine gap: standard SDID implicitly commits to an equal-weighted estimand, which is often not what applied researchers want. The ACA application is well-chosen because the dramatic sign reversal vividly illustrates the stakes. The exposition is clean and the modification itself is elegantly minimal. These are real strengths. But the paper stops well short of what is needed to establish the weighted estimator as a trustworthy tool.

**No formal identification or asymptotic theory in the main text**

The paper proposes a new estimator but states no formal assumptions under which it is consistent, derives no rate of convergence, and characterizes no asymptotic distribution. Arkhangelsky et al. (2021) establish their results under a latent factor model with specific regularity conditions. When the treated-unit averaging changes from uniform to arbitrary weights ω̃, the target of the synthetic control optimization shifts, and there is no guarantee that the original consistency argument carries over without modification. At minimum, the paper needs to state the weighted parallel trends or factor model condition under which τ̃ʷ converges to the target estimand τ^(ω̃), and verify that the regularity conditions in Arkhangelsky et al. still hold. Without this, a reader has no basis for trusting the estimator beyond the single application shown.

**Monte Carlo evidence is absent from the main text**

The conclusion mentions that finite-sample coverage simulations appear in Supplementary Appendix E, but these results are never summarized in the paper itself. For a methods paper, this is a serious gap. The reader needs to see, at minimum: (i) that the adapted bootstrap, jackknife, and placebo SEs achieve nominal coverage under known DGPs with heterogeneous treatment effects correlated with unit size; (ii) how coverage degrades as weight concentration increases (i.e., as effective N shrinks); and (iii) how the weighted estimator compares to standard SDID, SC, and DiD across DGPs where weighting should and should not matter. A table or figure in the main text showing coverage rates and RMSE across 3-4 DGPs would substantially strengthen the paper. Relegating all simulation evidence to an appendix leaves the paper's central methodological claim—that the adapted inference procedures work—unsubstantiated in the main narrative.

**Placebo variance estimator is inconsistent with the paper's own motivation**

Section 2.3 states that the placebo method uses uniform weights for pseudo-treated control units because the original ω̃ 'do not map naturally to control units.' But the entire premise of the paper is that unit heterogeneity makes equal weighting inappropriate. If the placebo SE estimator assumes homogeneous units while the actual estimand weights units heterogeneously, the resulting variance estimate may be biased—likely downward, since it misses the additional variance from weight concentration. The paper acknowledges alternative placebo weighting schemes exist (Appendix B) but does not resolve which one is valid or show that any of them achieves correct coverage. This needs a clear recommendation backed by simulation evidence showing that the chosen placebo procedure has correct size.

**Effective sample size collapse and inference reliability**

Population weighting reduces the effective number of treated units from 864 to 118 (Herfindahl index 0.0085), with the top ten counties holding 21% of the weight. The paper reports modestly larger SEs under population weighting (5.3 vs. 2.0) but does not investigate whether this modest increase is sufficient. With 46 state clusters and an effective N of 118, bootstrap confidence intervals may undercover. The paper should report empirical coverage of its bootstrap procedure under a DGP calibrated to the ACA application—specifically one with similar weight concentration and cluster structure. The phrase 'modestly larger' standard errors deserves scrutiny: a threefold increase (from 2 to 5.3) is not modest, and it would help to decompose how much comes from weight concentration versus clustering.

**Single application with no robustness to weight specification**

The paper's empirical evidence rests on one application with one set of weights (population shares). The domain-specific calibration notes correctly flag that sensitivity to weight choice matters. What happens with alternative weighting schemes—inverse variance weights, propensity-score-derived weights, or trimmed population weights that cap influence of the largest counties? If the -13.5 estimate is sensitive to, say, trimming the top 5% of counties by population, then the result may be driven by a handful of large counties rather than a general pattern. A robustness table showing estimates under 3-4 alternative weight specifications would directly address whether the finding is an artifact of extreme weight concentration.

**SUTVA concerns in the ACA application are unaddressed**

Treatment is assigned at the state level, but the paper analyzes county-level outcomes. Medicaid expansion in treated states may generate spillovers to neighboring counties in non-expansion states through provider networks, patient migration, or federal funding reallocation. These spillovers would contaminate the control group and bias the treatment effect estimate, plausibly toward zero for the equally-weighted estimator (if spillovers are diffuse) and away from zero for the population-weighted estimator (if spillovers concentrate in large border counties). The paper does not discuss SUTVA, test for border-county contamination, or cite the substantial literature on Medicaid spillovers. Even a brief discussion of why spillovers are unlikely to drive the sign reversal would be helpful.

**Staggered adoption is deferred but may affect the application**

The conclusion mentions staggered treatment adoption as a future extension, yet several states expanded Medicaid after 2014 (e.g., Montana in 2016, Louisiana in 2016). The paper appears to define treatment as a single 2014 event. If late-expanding states are classified as controls during 2014-2015 and then switch to treatment, they contaminate the control group. If they are dropped entirely, sample composition changes. The paper should clarify how late expanders are handled and whether results are sensitive to this choice. This is not a hypothetical concern—Clarke et al. (2023) develop staggered SDID specifically because collapsing heterogeneous adoption timing into a single cohort can bias estimates.

**No worked parametric example showing when weighting changes the estimand**

The paper's central claim is that weighting can flip the sign of the treatment effect. A transparent algebraic demonstration would make this concrete. Consider a two-treated-unit example under a simple interactive fixed effects model Y_it = alpha_i + gamma_t + delta_i * F_t + tau_i * W_it + epsilon_it, where unit 1 is large (population weight 0.9) with tau_1 = -15 and unit 2 is small (weight 0.1) with tau_2 = +10. Computing tau^sdid and tau^w explicitly under this DGP would show exactly how the sign reversal arises algebraically, and more importantly, would let readers see which features of the data-generating process make the equally-weighted estimator misleading. This kind of worked special case is standard in methods papers and is distinct from the Monte Carlo question (which concerns inference). Without it, readers must take the ACA application on faith rather than developing portable intuition.

**No diagnostic criteria for when practitioners should prefer weighted SDID**

The paper argues convincingly that weighting matters when treatment effects correlate with unit size, but never tells applied researchers how to detect this situation ex ante. A practical decision rule is missing. For instance, the paper could propose comparing pre-treatment outcome trajectories under equal and proposed weights—if these diverge substantially, the estimands differ and weighting is consequential. Or it could suggest a Hausman-type test comparing the two estimates. The evaluation standards for methods papers in this area call for decision rules that guide applied researchers, and the current draft offers none. The ACA application shows that weighting can matter dramatically, but a reader working on a different application has no guidance for assessing whether it will matter in their setting. Even an informal flowchart would help.

**Missing comparison to heterogeneity-robust DiD estimators**

The paper's motivation—treatment effect heterogeneity correlated with unit characteristics—is exactly the problem that Callaway and Sant'Anna (2021), Sun and Abraham (2021), and de Chaisemartin and D'Haultfoeuille (2020) address in the DiD context. The paper never discusses how weighted SDID relates to these estimators. This is a significant positioning gap. In the ACA application, where treatment timing is arguably uniform (2014 for the main expansion cohort), a Callaway-Sant'Anna estimate with population-weighted aggregation would provide a direct benchmark. Does weighted SDID offer something these estimators do not—better pre-treatment fit, smaller variance, robustness to violations of parallel trends that the time weights buy? The paper needs to make this case explicitly, or readers at a methods journal will wonder why they should adopt a new estimator when existing tools already handle heterogeneous effects.

**No decomposition of what drives the sign reversal in the application**

The paper attributes the sign change to treatment effects being heterogeneous and correlated with county size, but never directly demonstrates this correlation. Showing county-level or bin-level estimates plotted against population—even rough ones—would let readers verify the claimed mechanism. An alternative explanation is that a handful of very large counties with negative effects are driving the entire population-weighted result, which is a fragility concern rather than evidence of a general pattern. The paper reports that the top ten counties hold 21% of weight but does not show their individual contributions to the estimate. Computing the weighted SDID estimate with and without the ten largest counties, or showing how the estimate evolves as a function of a population cap, would distinguish between a broad pattern and a few influential outliers. This matters because the application is the paper's primary evidence that the extension is consequential.

**Recommendation**: Major revision. The paper identifies an important practical issue—the implicit equal-weighting choice in SDID—and proposes an elegant fix. But it currently reads as an extended application note rather than a methods contribution. The absence of formal identification conditions, the relegation of all simulation evidence to appendices, and the reliance on a single empirical application with no robustness checks leave the paper's central claims insufficiently supported for a methods journal.

**Key revision targets**:

1. State formal identification assumptions (e.g., a weighted factor model condition) under which the estimator is consistent for the target estimand, and verify that the regularity conditions of Arkhangelsky et al. (2021) carry through under non-uniform weights.
2. Include Monte Carlo simulations in the main text showing coverage rates and RMSE of the adapted variance estimators across at least 3 DGPs: one where weighting matters (heterogeneous effects correlated with size), one where it does not, and one calibrated to the ACA application's weight concentration and cluster structure.
3. Resolve the placebo inference inconsistency: recommend a specific placebo weighting scheme, justify it theoretically, and show via simulation that it achieves nominal coverage under weight concentration.
4. Add robustness checks for the ACA application: alternative weight specifications (trimmed population, inverse variance), sensitivity to treatment of late-expanding states, and a border-county placebo test addressing SUTVA.

**Status**: [Pending]

---

## Detailed Comments (11)

### 1. Asymmetric SE methods confound the central comparison

**Status**: [Pending]

**Quote**:
> Equally-weighted SEs: county-level jackknife. Population-weighted SEs: state-clustered bootstrap (200 replications).

**Feedback**:
The two columns of Table 1 use different inference procedures—county-level jackknife for equally-weighted, state-clustered bootstrap for population-weighted. The reported SE increase from ~2 to ~5.3 therefore conflates two effects: the impact of population weighting on precision and the switch from county-level to state-clustered inference. State-clustered bootstrap absorbs within-state correlation that the jackknife ignores, so part of the SE increase may simply reflect that the equally-weighted SEs are too small. Report state-clustered bootstrap SEs for the equally-weighted estimator as well (even as a supplementary column), so readers can isolate the precision cost of reweighting from the effect of clustering.

---

### 2. Cluster-robust inference presented as future extension despite being used in the application

**Status**: [Pending]

**Quote**:
> Important extensions include staggered treatment adoption, where the weighted estimator can be applied within each adoption cohort, and cluster-robust inference for settings where treatment assignment occurs at a higher level than the unit of observation (Supplementary Appendix G).

**Feedback**:
The conclusion lists cluster-robust inference as a future "important extension," but the ACA application already has exactly this structure—state-level treatment assignment with county-level outcomes—and Table 1 already reports state-clustered bootstrap SEs. A reader will be confused about whether clustering is a solved problem or an open one. Restructure to acknowledge that clustering is already operative in the reported results, then frame the extension as formal asymptotic theory for cluster-robust weighted SDID or generalization beyond the bootstrap approach used here.

---

### 3. DiD pre-trend failure undercuts the cross-estimator robustness claim

**Status**: [Pending]

**Quote**:
> The consistency of results across DiD and SDID confirms that the finding is robust to the choice of estimator.

**Feedback**:
Two paragraphs earlier, the paper reports that population-weighted DiD exhibits "non-trivial pre-treatment fluctuations, raising concerns about the parallel-trends assumption." If DiD's identifying assumption fails under population weighting, agreement between DiD and SDID does not constitute an informative robustness check—both estimators could be picking up the same pre-existing divergence. The sentence needs qualification. Something like: "though the DiD estimate under population weighting should be interpreted cautiously given the pre-treatment fluctuations noted above" would acknowledge the tension rather than let the robustness claim stand unqualified.

---

### 4. Hat/tilde notation inconsistency between introduction and formal definitions

**Status**: [Pending]

**Quote**:
> unit weights $\hat{\omega}$ and time weights $\tilde{\lambda}$ while allowing level differences between groups

**Feedback**:
The first paragraph of the introduction uses $\hat{\omega}$ for unit weights and $\tilde{\lambda}$ for time weights. Section 2.1 defines both with tildes ($\tilde{\omega} \in \Delta^{N_0}$, $\tilde{\lambda} \in \Delta^{T_0}$). The second introduction paragraph then writes $\tilde{\omega}$ and $\hat{\lambda}$, swapping the decorations again. Three different conventions across two pages will lose readers before they reach the formal setup. Pick one and apply it consistently.

---

### 5. Tilde-omega denotes both researcher-chosen and control weights

**Status**: [Pending]

**Quote**:
> ging with researcher-chosen weights $\tilde{\omega}_{j}$ in the collapsed form of the outcome matrix. The modification is minimal—one equation change—but it affects the estimation of both control unit weights $\tilde{\omega}$ and time weights $\hat{\lambda}$,

**Feedback**:
The researcher-specified treated-unit weights ($\tilde{\omega}_j$) and the standard SDID control weights ($\tilde{\omega}$) share the same letter and decorator, distinguished only by a subscript that could naturally index either set of units. These are fundamentally different objects—one is a user input, the other an optimization output. Section 2.2 introduces $\hat{\omega}^w$ for the weighted control weights, which helps, but the introduction has not yet made this distinction. Either preview the $\hat{\omega}^w$ notation here or use a different letter for the treated averaging weights.

---

### 6. Sign change claim needs the base estimate's sign stated explicitly

**Status**: [Pending]

**Quote**:
> The equally-weighted SDID estimate is near zero; the population-weighted estimate indicates a mortality reduction of approximately 13.5 deaths per 100,000

**Feedback**:
The introduction claims a "sign change" but describes the equally-weighted estimate only as "near zero" without stating its sign. Table 1 shows it is +0.35. For the sign-reversal claim to land, readers need the positive sign in the introduction—otherwise the narrative reads as a magnitude shift from approximately zero to -13.5, not a reversal. State the point estimate: "The equally-weighted SDID estimate is +0.35 deaths per 100,000" makes the sign change concrete.

---

### 7. "Modestly larger" understates a 2.65-fold SE increase

**Status**: [Pending]

**Quote**:
> Standard errors are modestly larger under population weighting (5.3 vs. 2 for SDID), reflecting greater heterogeneity in the reweighted sample.

**Feedback**:
An SE of 5.3 versus 2.0 is a factor of 2.65, meaning the confidence interval is 165% wider. That is not modest by any standard convention. The sentence also attributes the increase entirely to "greater heterogeneity in the reweighted sample," but as noted separately, part of the increase likely comes from switching to state-clustered inference. Replace "modestly larger" with "larger" or "substantially larger" and acknowledge the confound with the inference method change.

---

### 8. 90/10 ratio overstated as "order of magnitude"

**Status**: [Pending]

**Quote**:
> The population distribution is highly skewed: the 90th percentile county (176,243 residents) is an order of magnitude larger than the 10th percentile (5,928).

**Feedback**:
The ratio 176,243 / 5,928 is approximately 30—closer to 1.5 orders of magnitude than one. Correct to "roughly 30 times larger." The 90/10 ratio is also the wrong summary for weight concentration; the Herfindahl and top-share statistics reported two paragraphs later are more informative for the paper's argument. Consider moving those forward or replacing the 90/10 comparison.

---

### 9. Sample construction details needed for replicability

**Status**: [Pending]

**Quote**:
>  use a balanced county-year panel covering 2009–2017 with 2,443 counties (864 treated, 1,579 control), following *Borgschulte and Vogler (2020)*. The out

**Feedback**:
"Following" is ambiguous. Does this mean the sample is identical to Borgschulte and Vogler's, or that a similar panel was built from the same data sources? CDC WONDER suppresses counts below 10, so achieving a balanced panel necessarily excludes the smallest counties. Since the paper argues that weight concentration across heterogeneous units matters, readers need to know whether sample construction truncates the lower tail of the size distribution and by how much. One sentence stating the selection criterion (e.g., "counties with non-suppressed mortality in all nine years") would resolve this.

---

### 10. 46 states not identified

**Status**: [Pending]

**Quote**:
> Sample: 2,443 counties, 2009-2017. Outcome: all-cause mortality per 100,000, ages 20-64 (CDC WONDER). 46 states as clusters.

**Feedback**:
The US has 50 states plus DC. The table note says 46 but does not explain which are excluded or why. Likely candidates include DC and states dropped for data availability, but the reader should not have to guess. State the exclusion criterion so readers can assess whether the dropped jurisdictions matter for external validity.

---

### 11. 200 bootstrap replications may be too few for stable confidence intervals

**Status**: [Pending]

**Quote**:
> Population-weighted SEs: state-clustered bootstrap (200 replications).

**Feedback**:
While the simulation error in the bootstrap SE is tolerable at ~5%, percentile-based 95% confidence intervals from 200 draws depend on the 5th and 195th order statistics, which can be quite volatile. Standard practice for published work is at least 999 replications. If computation is binding, note that explicitly; otherwise increase to 999 and report whether the SEs and CI endpoints change.

---
