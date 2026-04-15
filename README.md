# synthdid: Weighted Synthetic Difference-in-Differences

This is a fork of the [`synthdid`](https://github.com/synth-inference/synthdid) R package that extends the synthetic difference-in-differences estimator of Arkhangelsky et al. (2021) with **user-specified weights on treated units**.

The standard SDID estimator averages treated units equally. When treated units differ in size---for example, counties ranging from 10,000 to 10 million residents---the equally-weighted and population-weighted estimands can diverge substantially. This package allows the researcher to specify weights that match their estimand of interest.

### What this fork adds

- **`synthdid_estimate_weighted()`**: weighted SDID with researcher-specified treated-unit weights
- **`sc_estimate_weighted()`**, **`did_estimate_weighted()`**: weighted SC and DiD variants
- **`synthdid_event_study()`**: event-study decomposition for both weighted and unweighted estimates
- **Adapted variance estimators**: bootstrap, jackknife, and placebo SEs with weight renormalization
- **Cluster-robust standard errors**: cluster bootstrap and jackknife for settings where treatment is assigned at a higher level than the unit of observation
- **Placebo weight options**: uniform, size-matched, and permuted weights for placebo inference
- **Effective sample size tuning**: optional Kish's N<sub>eff</sub> adjustment for the regularization parameter

### Installation

```R
devtools::install_github("johniselin-econ/synthdid_weights")
```

### Quick example

```R
library(synthdid)

# Standard SDID (equal weights)
data('california_prop99')
setup = panel.matrices(california_prop99)
tau.hat = synthdid_estimate(setup$Y, setup$N0, setup$T0)
sprintf('point estimate: %1.2f', tau.hat)

# Weighted SDID (custom weights on treated units)
N1 = nrow(setup$Y) - setup$N0
pop_weights = runif(N1)  # replace with real population weights
tau.w = synthdid_estimate_weighted(setup$Y, setup$N0, setup$T0,
          treated.weights = pop_weights)
sprintf('weighted estimate: %1.2f', tau.w)

# Cluster-robust standard errors
cluster = rep(1:5, length.out = nrow(setup$Y))
tau.cl = synthdid_estimate_weighted(setup$Y, setup$N0, setup$T0,
           treated.weights = pop_weights, cluster = cluster)
se = sqrt(vcov(tau.cl, method = "bootstrap", replications = 200))

# Event study
es = synthdid_event_study(tau.w)
plot_event_study(es)
```

### Paper

John Iselin and Erica Ryan. **Weighted Synthetic Difference-in-Differences**. Working paper, 2026.

The paper applies the weighted estimator to the ACA Medicaid expansion, following Borgschulte and Vogler (2020). The equally-weighted SDID estimate is near zero, while the population-weighted estimate indicates a mortality reduction of approximately 13.5 deaths per 100,000---demonstrating that the weighting choice can be economically consequential.

### References

Dmitry Arkhangelsky, Susan Athey, David A. Hirshberg, Guido W. Imbens, and Stefan Wager.
**Synthetic Difference in Differences.** *American Economic Review*, 111(12):4088-4118, 2021.

Damian Clarke, Daniel Pailanir, Susan Athey, and Guido Imbens.
**Synthetic Difference-in-Differences Estimation.** IZA Discussion Paper No. 15907, 2023.

Mark Borgschulte and Jacob Vogler.
**Did the ACA Medicaid Expansion Save Lives?** *Journal of Health Economics*, 72:102333, 2020.

### Upstream

This fork is based on [synth-inference/synthdid](https://github.com/synth-inference/synthdid). The original package API is fully preserved; all weighted functions are additive extensions.

---

**Note:** The weighted vignette (`vignettes/synthdid_weighted.Rmd`) should be updated to use `library(synthdid)` once the weighted functions are merged into the installed package. Currently it sources `../R/` files directly.
