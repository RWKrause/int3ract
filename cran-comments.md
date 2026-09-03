## Test environments

* Windows 11, R 4.6.1 (local)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Notes on this submission

* This is a feature release. It adds a `jn_input()` method for `lavaan`
  structural equation and path models, with a new `outcome` argument naming
  the dependent variable whose regression is to be analysed, since a
  structural model usually has more than one. Nothing in the existing
  interface changes and no existing behaviour is affected.

* `lavaan` is used only by that method, and is in Suggests.

* There are no reverse dependencies.
