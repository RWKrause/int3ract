# int3ract 1.0.8

# int3ract 1.0.7
* Fixing stray files in package for CRAN submission
# int3ract 1.0.6
* Fixing save directory for CRAN submission
# int3ract 1.0.5
* Fixing typos for CRAN submission
# int3ract 1.0.2

* Initial CRAN submission.
* Implements the Johnson–Neyman (JN - Johnson and Fay (1950) <doi:10.1007/BF02288864>) technique for two-way
interactions and its three-way extension (JNK technique).
* `JNK_freq()` supports `lm`, `glm`, `lme4::lmer`/`glmer`, and
`RSiena::siena` objects, as well as raw coefficient vectors and
covariance matrices.
* `JNK_bayes()` supports `multiSiena`/`sienaBayesFit` objects and
raw posterior-draw matrices; produces density plots for two-way
interactions and posterior-mean heatmaps with Bayesian p-value
overlays for three-way interactions.
