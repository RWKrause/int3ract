## R CMD check results

0 errors | 0 warnings | 0 notes

Checked on Windows (R 4.6.1), win-builder (R-devel and R-release) and the
macOS builder.

## Notes on this submission

* This is a major version. The user-facing interface is now the single generic
  `JN()`, which dispatches on the fitted object. The previous entry points
  `JNK_freq()` and `JNK_bayes()` are deprecated rather than removed: they still
  work, delegate to `JN()`, return the old layout and warn once per session.

* It also fixes a bug in the two-way standard errors. Versions 1.0.x used the
  covariance between the two main effects where the delta method calls for the
  covariance between the focal main effect and the interaction, so two-way
  frequentist results produced by those versions change. This is documented in
  NEWS.md, in the README and in the vignette.

* There are no reverse dependencies.
