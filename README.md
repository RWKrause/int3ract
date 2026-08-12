# int3ract

Visualizing 3-way interactions and obtaining regions of significance.

Interaction effects are ubiquitous in applied statistical modelling, yet their meaningful interpretation remains challenging. The classic Johnson-Neyman (JN) technique (Johnson and Neyman 1936) addresses this challenge for two-way interactions by identifying the regions of a moderator's range over which a focal effect is and is not statistically significant. The int3ract package for R implements the JN technique and its three-way extension (the Johnson-Neyman-Krause, or JNK, technique) for both frequentist and Bayesian models. The function JNK_freq() auto-detects models fitted via lm()/glm(), RSiena's siena(), or lme4's lmer()/glmer(), but can also be applied to multiplicative interactions from (virtually) any model family by supplying a coefficient vector and covariance matrix directly. For Bayesian Stochastic Actor-Oriented Models (SAOMs) estimated with multiSiena, or any model producing posterior draws, the function JNK_bayes() produces conditional posterior distributions. For two-way interactions, classic shaded confidence-band plots are created that visually demarcate significant and non-significant regions along the moderator range; three-way interactions yield colour-gradient heatmaps with optional crosshatch overlays for non-significant regions. The package is designed to encourage richer, region-specific reporting of interaction effects in place of the conventional single-slope spotlight approach.

## Installation

```r
# CRAN (once available):
install.packages("int3ract")

# Development version:
# install.packages("remotes")
remotes::install_github("RWKrause/int3ract")
```

## Quick example

```r
library(int3ract)

dat <- data.frame(y = rnorm(100), x = rnorm(100),
                  z = rnorm(100), w = rnorm(100))
res <- lm(y ~ x * z * w, dat)

# two-way JN plot
JNK_freq(res, theta_1 = "x", theta_2 = "z")$plots$z

# three-way JNK heatmap
JNK_freq(res, theta_1 = "x", theta_2 = "z", theta_3 = "w")$plots$z
```

## Reference

Krause, R. W. (2026). *int3ract: Johnson-Neyman Technique and its Three-Way
Extension for Frequentist and Bayesian Models in R.* arXiv:2604.22051.
<https://doi.org/10.48550/arXiv.2604.22051>