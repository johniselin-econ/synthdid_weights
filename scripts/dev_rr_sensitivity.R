## DEV REFERENCE (press-side) -- this computation belongs in the supplement
## Rmd (E.6), reading the committed results/event_study_draws.csv; it is NOT
## a server block. Keep alongside the supplement chunk until E.6 is written.
##
## Rambachan-Roth-style (M = 0) trend sensitivity for the population-weighted
## SDID event study, using the committed state-clustered bootstrap draws.
## Per draw: fit an OLS line through the PRE-period coefficients (2009-2013),
## extrapolate it through 2014-2017, and compute the trend-adjusted average
## post effect. The draw distribution gives the bootstrap CI of the adjusted
## estimator (slope-estimation uncertainty included by construction).
suppressPackageStartupMessages({library(dplyr); library(readr)})

draws <- read_csv("results/event_study_draws.csv", show_col_types = FALSE)
pt    <- read_csv("results/event_studies.csv",     show_col_types = FALSE) %>%
  filter(Estimator == "SDID", Weighting == "Population weighted") %>%
  select(year, estimate)

adjust <- function(yr, cf) {
  pre <- yr <= 2013
  m   <- stats::lm(cf[pre] ~ yr[pre])
  b0  <- unname(stats::coef(m)[1]); b1 <- unname(stats::coef(m)[2])
  c(raw = mean(cf[!pre]), adj = mean(cf[!pre] - (b0 + b1 * yr[!pre])), slope = b1)
}

point <- adjust(pt$year, pt$estimate)
boot  <- draws %>% group_by(rep) %>%
  summarise(res = list(adjust(year, value)), .groups = "drop") %>%
  mutate(raw = sapply(res, `[[`, "raw"),
         adj = sapply(res, `[[`, "adj"),
         slope = sapply(res, `[[`, "slope"))

ci <- function(x) stats::quantile(x, c(.025, .975))
cat("Population-weighted SDID event study, 2009-2017 (500 state-clustered draws)\n\n")
cat(sprintf("raw average post effect      : %6.2f  (boot SE %4.2f, 95%% CI [%6.2f, %6.2f])\n",
            point["raw"], sd(boot$raw), ci(boot$raw)[1], ci(boot$raw)[2]))
cat(sprintf("pre-period differential slope: %6.2f per year (SE %4.2f)\n",
            point["slope"], sd(boot$slope)))
cat(sprintf("trend-adjusted (M = 0) effect: %6.2f  (boot SE %4.2f, 95%% CI [%6.2f, %6.2f])\n",
            point["adj"], sd(boot$adj), ci(boot$adj)[1], ci(boot$adj)[2]))
cat(sprintf("\nshare of CI draws: adj < 0: %.0f%%   adj < -11.36 (BV): %.0f%%   adj < -17.45 (raw): %.0f%%\n",
            100 * mean(boot$adj < 0), 100 * mean(boot$adj < -11.36), 100 * mean(boot$adj < -17.45)))
cat(sprintf("\n[point-estimate counterparts: detrended SDID -2.29 (SE 4.35); binned -4.77 (SE 3.93)]\n"))
