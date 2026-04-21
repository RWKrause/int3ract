# int3ract 1.0.1

* Initial CRAN submission.
* Implements the Johnson–Neyman (JN) technique for two-way
interactions and its three-way extension (JNK technique).
* `JNK_freq()` supports `lm`, `glm`, `lme4::lmer`/`glmer`, and
`RSiena::siena` objects, as well as raw coefficient vectors and
covariance matrices.
* `JNK_bayes()` supports `multiSiena`/`sienaBayesFit` objects and
raw posterior-draw matrices; produces density plots for two-way
interactions and posterior-mean heatmaps with Bayesian p-value
overlays for three-way interactions.