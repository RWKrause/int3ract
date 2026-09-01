# int3ract

Johnson-Neyman analysis of two- and three-way interactions, for frequentist
and Bayesian models.

An interaction coefficient says little on its own. The classic Johnson-Neyman
(JN) technique (Johnson and Neyman 1936) reports instead the *region* of the
moderator's range over which the focal effect is distinguishable from zero.
**int3ract** implements that technique, extends it to three-way interactions
over a two-dimensional moderator grid (JN3), and applies it to Bayesian models
by working on posterior draws rather than point estimates.

The package has one entry point, `JN()`, which dispatches on the fitted
object. Models carrying point estimates and a covariance matrix are analysed
with Wald tests; objects carrying draws are analysed as conditional posteriors.
Results are classed objects with `print()`, `summary()` and `plot()` methods.

## Installation

```r
# CRAN:
install.packages("int3ract")

# Development version:
# install.packages("remotes")
remotes::install_github("RWKrause/int3ract")
```

## Quick example

```r
library(int3ract)

set.seed(1402)
dat <- data.frame(x = rnorm(100), z = rnorm(100), w = rnorm(100))
dat$y <- dat$x + 0.5 * dat$z - 0.5 * dat$w +
  0.5 * dat$x * dat$z * dat$w + rnorm(100, sd = 4)

jn <- JN(lm(y ~ x * z, data = dat), theta_1 = "x", theta_2 = "z")

summary(jn)   # the regions of significance
plot(jn, which = "x")

# three-way, over a two-dimensional moderator grid
jn3 <- JN(lm(y ~ x * z * w, data = dat),
          theta_1 = "x", theta_2 = "z", theta_3 = "w")
plot(jn3, which = "x")

# any matrix of posterior draws takes the same route
post <- as.matrix(MCMCpack::MCMCregress(y ~ x * z, data = dat))
JN(post, theta_1 = "x", theta_2 = "z",
   theta_1_vals = seq(-3, 3, 0.5), theta_2_vals = seq(-3, 3, 0.5))
```

## Stochastic actor-oriented models

`RSiena` and `multiSiena` results are handled directly. Effects are addressed
by integer position, because effect names are not unique within a model; if a
position is out of range the error lists the available effects.

```r
# RSiena: siena() / siena07() -> Wald tests
JN(saom_fit, theta_1 = 4, theta_2 = 7, theta_int_12 = 12,
   theta_1_vals = c(0, 6), theta_2_vals = c(-2, 2))

# multiSiena: sienaBayes() -> conditional posteriors
JN(bayes_fit, theta_1 = 4, theta_2 = 7, theta_int_12 = 12,
   theta_1_vals = seq(0, 6, 1), theta_2_vals = seq(-2, 2, 1))
```

For `sienaBayes()` results the package works out for itself whether each
parameter was estimated as shared across groups (`eta`) or as varying between
them (`mu`), takes the corresponding draws, and reports which it used.
`hyper_only = FALSE` adds one analysis per group, exactly as
`fixed_only = FALSE` does for `lme4` models. It is skipped, with a message,
when no parameter involved varies between groups — the group plots would
otherwise be identical copies.

## Data support

A region of significance covering moderator values that were hardly observed
is not evidence of much. Where the observed moderator values are known --
automatically for `lm`, `glm` and `lme4` models, and through the `support`
argument otherwise -- every figure carries a histogram of them, and
`summary()` reports the share of observations falling inside each region.

## Supporting further model classes

`JN()` works on any object with a `jn_input()` method. Writing one means
returning either point estimates with their covariance matrix, or a matrix of
draws:

```r
jn_input.myfit <- function(object, theta_1, theta_2, theta_3 = NULL, ...) {
  idx <- c(theta_1, theta_2, paste(theta_1, theta_2, sep = ":"))
  jn_wald(coefficients = object$estimates[idx],
          vcov         = object$covariance[idx, idx],
          labels       = c(theta_1, theta_2))
}
```

Nothing else needs to change: `print()`, `summary()` and `plot()` work on the
result immediately. See `?jn_input`.

## Migrating from 1.0.x

`JNK_freq()` and `JNK_bayes()` are deprecated in favour of `JN()`. They still
work and return the old layout, but warn once per session.

Note that version 2.0.0 **fixes a bug in the two-way standard errors**: 1.0.x
used the covariance between the two main effects where the delta method calls
for the covariance between the focal main effect and the interaction. Two-way
frequentist results from earlier versions should be regenerated. See
`NEWS.md`.

## Reference

Krause, R. W. (2026). *int3ract: Johnson-Neyman Technique and its Three-Way
Extension for Frequentist and Bayesian Models in R.* arXiv:2604.22051.
<https://doi.org/10.48550/arXiv.2604.22051>
