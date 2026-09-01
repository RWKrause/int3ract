#' int3ract: Johnson-Neyman analysis of two- and three-way interactions
#'
#' Interaction coefficients say little on their own. \pkg{int3ract} reports and
#' draws the conditional effect of each variable involved in a multiplicative
#' interaction across the range of its moderators, together with the region
#' over which that effect is distinguishable from zero.
#'
#' The package has one entry point, \code{\link{JN}}, which dispatches on the
#' fitted object. Models carrying point estimates and a covariance matrix are
#' analysed with Wald tests; objects carrying draws are analysed as conditional
#' posteriors. Adding support for a further model class means writing one
#' \code{\link{jn_input}} method for it.
#'
#' @section Main functions:
#' \describe{
#'   \item{\code{\link{JN}}}{run the analysis.}
#'   \item{\code{\link{jn_input}}}{the extension point for new model classes.}
#'   \item{\code{\link{summary.JN}}}{the regions of significance.}
#'   \item{\code{\link{plot.JN}}}{the figures, with data-density panels.}
#'   \item{\code{\link{jn_style}}}{their appearance.}
#' }
#'
#' @importFrom stats coef model.frame pnorm qnorm quantile sd setNames vcov
#' @importFrom utils head modifyList
#' @importFrom grDevices dev.interactive devAskNewPage
#' @importFrom ggplot2 autoplot ggsave
#' @import ggplot2
#' @import ggpattern
#' @import patchwork
#' @keywords internal
"_PACKAGE"

utils::globalVariables(".data")
