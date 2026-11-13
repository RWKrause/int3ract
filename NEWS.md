# int3ract 2.2.0

## New features

* Simple-slopes figures: `plot(x, type = "slopes")` draws the regions of
  significance in outcome space, for every kind of analysis. The focal
  variable is on the x axis, the model's prediction of the outcome on the y
  axis, and simple slopes are drawn at the moderator values bounding each
  region, the fan between them shaded in the heatmap's positive (red) or
  negative (blue) colour. That gives one dashed slope at the moderator mean
  when the effect is nowhere significant, two slopes for a region of one sign,
  and four for regions of both signs. An effect significant throughout gets
  the slopes at the lowest and highest observed moderator values and a dotted
  one at the most common value. Regions are computed over the observed
  moderator values only. `jn_style(slope_shading = "density")` shades the fan
  by how many observations lie at each moderator value.

  Three-way analyses are drawn as slices: one panel per value of one
  moderator (its quartiles by default, or the new `at` argument of `plot()`,
  `autoplot()` and `jn_plots()`), each the two-way figure along the other.
  For Wald analyses the boundaries within each panel are solved exactly.

  The outcome axis is on the linear-predictor scale with other covariates at
  their means. The constant this needs is extracted from `lm`, `glm`, `lme4`
  and `lavaan` fits, and can be given to `jn_wald()` and `jn_posterior()`
  through their new `intercept` and `outcome_label` arguments; without it the
  axis is labelled as being up to a constant.

* Binary moderators. A variable observed at exactly two values is evaluated
  at those values only. Its two-way figure shows the conditional effect with
  its interval at each value; its simple-slopes figure shows both slopes,
  each marked significant or not. In a three-way analysis the heatmap is
  replaced by the two-way figure along the other moderator, once for each
  value of the binary one. `jn_regions()` and `summary()` report the values at
  which the effect is significant. The new `binary` argument of `JN()` names
  binary variables where the observed values are not available.

* `JN()` gains `extend` (default 0.1): every moderator range, whether taken
  from the data or supplied as a length-two `theta_*_vals`, is widened by this
  fraction on either side. The margin is flagged in the new `extrapolated`
  column of the result tables and drawn fainter in the figures
  (`jn_style(extrap_alpha = )`), with the observed range marked. `summary()`
  reports the observed range alongside the evaluated one.

## Changes in default output

* Because of `extend`, default grids now reach 10% beyond the observed (or
  supplied) range, so regions of significance that ran to the edge of the
  range may now extend further, and printed boundaries and grids differ from
  2.1.0. `extend = 0` reproduces the earlier results exactly. The deprecated
  `JNK_freq()` and `JNK_bayes()` keep their 1.x behaviour.


# int3ract 2.1.0

## New features

* Structural equation and path models fitted with `lavaan` are supported.
  A structural model usually has several endogenous variables, and the same
  interaction may predict more than one of them with different coefficients,
  so the new `outcome` argument names the dependent variable whose regression
  is wanted. It may be omitted when only one outcome is regressed on all of
  the focal variables; when several are, `JN()` stops and lists them rather
  than guessing.

  `lavaan` writes two-way products as `x:z` in its own syntax and those are
  found automatically, but it does not accept three-way products at all: those
  have to be built as variables in the data, and their names passed as
  `theta_int_12`, `theta_int_13`, `theta_int_23` and `theta_int_123`.

  Latent outcomes need nothing special. Multiple-group fits are declined
  rather than pooled, since each group carries its own structural
  coefficients, and `blavaan` fits are declined with a pointer to passing
  their draws, since a Wald test on posterior means would discard the
  posterior.


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

* The default `thresholds` for Bayesian analyses are now 
  `c(alpha/2, 1 - alpha/2)`, that is `c(0.025, 0.975)`. The previous default of
  `c(0.49999999999999999, 0.5)` marked essentially every cell significant.
  `JNK_bayes()` keeps the old default.

* `control_fdr` is documented as what it is: the Benjamini-Hochberg step-up
  procedure controlling the false discovery rate across the moderator grid.
  Earlier documentation described it as Bonferroni-Holm.

* `dplyr`, `tibble`, `tidyr` and `scales` are no longer dependencies. `lme4`
  moved from Imports to Suggests, since it is needed only for mixed-model
  input. `patchwork` was added, for the density panels.

* The package now has a test suite and a `Getting started with int3ract`
  vignette.

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
