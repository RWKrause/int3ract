# int3ract 2.0.0

## Bug fix

* The two-way standard errors were wrong. Following the delta method, the
  variance of the conditional effect is
  `Var(b_focal) + 2m Cov(b_focal, b_interaction) + m^2 Var(b_interaction)`, but
  the covariance term used was `Cov(b_focal, b_other_main_effect)` instead of
  the covariance with the interaction. Confidence bands, `p` values and
  therefore the regions of significance of every two-way frequentist analysis
  produced by version 1.0.x are affected, and results from those versions
  should be regenerated. Three-way analyses and all Bayesian analyses were
  computed correctly and are unchanged.

## New interface

* New generic `JN()`, which dispatches on the fitted object. It replaces
  `JNK_freq()` and `JNK_bayes()`, which split the package by inference
  paradigm rather than by input. Whether an analysis is frequentist or Bayesian
  now follows from the object handed in: models carrying point estimates and a
  covariance matrix are analysed with Wald tests, objects carrying draws as
  conditional posteriors.

* New generic `jn_input()`, the extension point of the package, with methods
  for `lm`, `glm`, `merMod`, `sienaFit`, `matrix`, `mcmc`, `mcmc.list`,
  `data.frame`, `sienaBayesFit` and `multiSiena`. Supporting a further model
  class means writing one `jn_input()` method returning either `jn_wald()`
  (point estimates and covariance) or `jn_posterior()` (draws); nothing else
  needs to change. This covers, for instance, models fitted by variational
  inference, where the variational parameters approximate the posterior.

* `JN()` returns a classed object with `print()`, `summary()`, `plot()`,
  `autoplot()`, `as.data.frame()`, `coef()` and `vcov()` methods. Analyses
  producing one result per group return a `JN_list` carrying the same methods.
  Figures are built on demand by `plot()` rather than stored in the result.

* `JNK_freq()` and `JNK_bayes()` are deprecated. They still work, delegate to
  `JN()`, and return the 1.0.x layout, but warn once per session. They will be
  removed in a future release.

## New features

* `summary()` reports the regions of significance directly: the stretches of
  the moderator over which each conditional effect is distinguishable from
  zero, with the share of observations falling inside each. For two-way Wald
  analyses the boundaries are solved for exactly, so their precision no longer
  depends on `range_size`. Also available programmatically as `jn_regions()`.

* `plot()` attaches data-density panels: a histogram of the observed moderator
  values below two-way plots, and marginal histograms along both axes of the
  heatmaps. A region of significance covering moderator values that were
  hardly observed is the standard failure mode of this technique, and the
  panels make it visible. Observed values are taken from `lm`, `glm` and
  `lme4` models automatically, and can be supplied for any other input through
  the `support` argument. Turn the panels off with
  `jn_style(show_density = FALSE)`.

* New `jn_style()` collects the colour and pattern settings that were
  previously repeated across both function signatures.

* `plot()` draws a figure for every variable involved in the interaction,
  since each takes its turn as the focal one: two figures for a two-way
  analysis, three for a three-way one. Pass `which` for a single figure.
  It draws and returns invisibly, so nothing is rendered twice.

* New `jn_plots()` returns every figure of an analysis at once without
  drawing them, and `autoplot()` returns a single one; `jn_save()` writes
  them to disk, replacing the `save` and `folder` arguments.

* Two-way Bayesian analyses can be drawn as a posterior mean with credible
  band against the moderator (`plot(x, type = "band")`), directly comparable
  to the frequentist figure, in addition to the overlaid conditional densities.

* For `sienaBayes()` results, the level the draws were taken from is now
  recorded on the result and shown by `print()`, `summary()` and the figures:
  `Eta` when the parameters involved are shared across groups, `Mu` when they
  vary and the hyper-parameter is used. Version 1.0.x labelled the figures this
  way; the information had been dropped.

* `hyper_only = FALSE` no longer produces group-level analyses when none of the
  parameters involved varies between groups. It now reports that and returns
  the population-level analysis alone, since the group analyses would have been
  identical copies.

## Other changes

* The three-way extension is no longer called the Johnson-Neyman-Krause or JNK
  technique; it is referred to as JN3.

* The default `thresholds` for Bayesian analyses are now `c(alpha/2, 1 -
  alpha/2)`, that is `c(0.025, 0.975)`. The previous default of
  `c(0.49999999999999999, 0.5)` marked essentially every cell significant.
  `JNK_bayes()` keeps the old default.

* `control_fdr` is documented as what it is: the Benjamini-Hochberg step-up
  procedure controlling the false discovery rate across the moderator grid.
  Earlier documentation described it as Bonferroni-Holm.

* `dplyr`, `tibble`, `tidyr` and `scales` are no longer dependencies. `lme4`
  moved from Imports to Suggests, since it is needed only for mixed-model
  input. `patchwork` was added, for the density panels.

* The package now has a test suite.

# int3ract 1.0.7
* Fixing stray files in package for CRAN submission

# int3ract 1.0.6
* Fixing save directory for CRAN submission

# int3ract 1.0.5
* Fixing typos for CRAN submission

# int3ract 1.0.2

* Initial CRAN submission.
* Implements the Johnson-Neyman (JN - Johnson and Fay (1950)
<doi:10.1007/BF02288864>) technique for two-way interactions and its three-way
extension.
* `JNK_freq()` supports `lm`, `glm`, `lme4::lmer`/`glmer`, and
`RSiena::siena` objects, as well as raw coefficient vectors and
covariance matrices.
* `JNK_bayes()` supports `multiSiena`/`sienaBayesFit` objects and
raw posterior-draw matrices; produces density plots for two-way
interactions and posterior-mean heatmaps with Bayesian p-value
overlays for three-way interactions.
