# =============================================================================
# JN(): the user-facing generic.
#
# Dispatch is on the fitted object. Whether the analysis ends up frequentist or
# Bayesian is decided by which carrier the object's jn_input() method returns,
# not by which function the user called.
# =============================================================================

#' Johnson-Neyman analysis of a two- or three-way interaction
#'
#' Computes the conditional effect of each variable involved in a
#' multiplicative interaction across the range of its moderators, together with
#' the regions over which that effect is statistically distinguishable from
#' zero. Two-way interactions give the classic Johnson-Neyman analysis;
#' three-way interactions give its extension over a two-dimensional moderator
#' grid, referred to here as JN3.
#'
#' \code{JN()} dispatches on \code{object}. Fitted models carrying point
#' estimates and a covariance matrix (\code{lm}, \code{glm}, \code{lmerMod},
#' \code{sienaFit}, ...) are analysed with Wald \eqn{z} tests; objects carrying
#' draws (a posterior matrix, \code{mcmc}, \code{sienaBayesFit}, ...) are
#' analysed as conditional posterior distributions with Bayesian \eqn{p}
#' values. Both paths return an object of class \code{JN} with the same set of
#' methods.
#'
#' Support for a model class that is not listed above is added by writing a
#' single \code{\link{jn_input}} method for it; \code{JN()} itself needs no
#' change. See \code{\link{jn_input}}.
#'
#' @param object a fitted model, a matrix of draws, or a \code{jn_input}
#'   object built with \code{\link{jn_wald}} or \code{\link{jn_posterior}}.
#' @param theta_1,theta_2 name (or, for \code{sienaFit} and
#'   \code{sienaBayesFit}, integer position) of the first and second variable
#'   involved in the interaction.
#' @param theta_3 the third variable, or \code{NULL} (default) for a two-way
#'   interaction.
#' @param theta_1_vals,theta_2_vals,theta_3_vals the values each variable takes
#'   when it acts as a moderator. A vector of length two is read as a range and
#'   filled in with \code{range_size} points; any other length is used
#'   verbatim. Defaults to the observed range where the model carries its data,
#'   and is required otherwise.
#' @param support optional list of length 2 or 3 with the observed values of
#'   each variable, used for the data-density panels of \code{\link{plot.JN}}.
#'   Extracted automatically from \code{lm}, \code{glm} and \pkg{lme4} models.
#' @param alpha significance level for the Wald path, and the basis of the
#'   default \code{thresholds} for the posterior path. Default 0.05.
#' @param control_fdr control the false discovery rate across the moderator
#'   grid with the Benjamini-Hochberg procedure? Default \code{FALSE}.
#' @param range_size number of moderator values used when a range is expanded.
#'   Defaults to 1000 for two-way and 50 for three-way Wald analyses, and to 13
#'   and 50 respectively for posterior analyses, where each value becomes its
#'   own conditional posterior.
#' @param thresholds two Bayesian \eqn{p} values bounding the region treated as
#'   inconclusive; the conditional effect counts as significant outside them.
#'   Defaults to \code{c(alpha/2, 1 - alpha/2)}.
#' @param ... passed to the \code{\link{jn_input}} method for \code{object}.
#'
#' @returns An object of class \code{JN}, with methods
#'   \code{\link[=print.JN]{print}}, \code{\link[=summary.JN]{summary}},
#'   \code{\link[=plot.JN]{plot}} and
#'   \code{\link[=as.data.frame.JN]{as.data.frame}}. Analyses producing one
#'   result per group return a \code{JN_list}, which carries the same methods.
#'
#' @references
#' Johnson PO, Neyman J (1936). "Tests of Certain Linear Hypotheses and Their
#' Application to Some Educational Problems." \emph{Statistical Research
#' Memoirs}, 1, 57-93.
#'
#' Bauer DJ, Curran PJ (2005). "Probing Interactions in Fixed and Multilevel
#' Regression: Inferential and Graphical Techniques." \emph{Multivariate
#' Behavioral Research}, 40(3), 373-400.
#' \doi{10.1207/s15327906mbr4003_5}
#'
#' @examples
#' set.seed(1402)
#' dat <- data.frame(x = rnorm(100), z = rnorm(100), w = rnorm(100))
#' dat$y <- dat$x + 0.5 * dat$z - 0.5 * dat$w +
#'   0.5 * dat$x * dat$z * dat$w + rnorm(100, sd = 4)
#'
#' ## two-way
#' fit2 <- lm(y ~ x * z, data = dat)
#' jn2  <- JN(fit2, theta_1 = "x", theta_2 = "z")
#' jn2
#' summary(jn2)
#' plot(jn2, which = "x")
#'
#' ## three-way
#' fit3 <- lm(y ~ x * z * w, data = dat)
#' jn3  <- JN(fit3, theta_1 = "x", theta_2 = "z", theta_3 = "w",
#'            range_size = 20)
#' summary(jn3)
#'
#' ## a matrix of posterior draws takes the same route
#' post <- cbind(x = rnorm(500, 0.5, 0.2), z = rnorm(500, -0.3, 0.2),
#'               `x:z` = rnorm(500, 0.4, 0.2))
#' JN(post, theta_1 = "x", theta_2 = "z",
#'    theta_1_vals = seq(-2, 2, 1), theta_2_vals = seq(-2, 2, 1))
#'
#' @export
JN <- function(object, ...) UseMethod("JN")


#' @rdname JN
#' @export
JN.default <- function(object,
                       theta_1, theta_2, theta_3 = NULL,
                       theta_1_vals = NULL, theta_2_vals = NULL,
                       theta_3_vals = NULL,
                       support = NULL,
                       alpha = 0.05,
                       control_fdr = FALSE,
                       range_size = NULL,
                       thresholds = NULL,
                       ...) {
  input <- jn_input(object,
                    theta_1 = theta_1, theta_2 = theta_2, theta_3 = theta_3,
                    ...)
  if (!is.null(support))
    input$support <- .jn_check_support(support, length(input$labels))

  .jn_build(input,
            ranges      = .jn_ranges_arg(theta_1_vals, theta_2_vals,
                                         theta_3_vals, theta_3),
            range_size  = range_size,
            alpha       = alpha,
            control_fdr = control_fdr,
            thresholds  = thresholds %||% c(alpha / 2, 1 - alpha / 2))
}


#' @rdname JN
#' @export
JN.jn_input <- function(object,
                        theta_1_vals = NULL, theta_2_vals = NULL,
                        theta_3_vals = NULL,
                        alpha = 0.05, control_fdr = FALSE,
                        range_size = NULL, thresholds = NULL, ...) {
  .jn_build(object,
            ranges      = .jn_ranges_arg(theta_1_vals, theta_2_vals,
                                         theta_3_vals,
                                         if (length(object$labels) == 3L) TRUE
                                         else NULL),
            range_size  = range_size,
            alpha       = alpha,
            control_fdr = control_fdr,
            thresholds  = thresholds %||% c(alpha / 2, 1 - alpha / 2))
}


#' @rdname JN
#' @param fixed_only analyse only the fixed effects? When \code{FALSE}, a
#'   separate analysis is returned for every level of the grouping factor,
#'   using group-specific coefficients (fixed effects plus conditional modes).
#' @param group_var name of the grouping factor to use when
#'   \code{fixed_only = FALSE}. Defaults to the first one.
#' @export
JN.merMod <- function(object,
                      theta_1, theta_2, theta_3 = NULL,
                      theta_1_vals = NULL, theta_2_vals = NULL,
                      theta_3_vals = NULL,
                      support = NULL,
                      alpha = 0.05, control_fdr = FALSE,
                      range_size = NULL, thresholds = NULL,
                      fixed_only = TRUE, group_var = NULL, ...) {

  fixed <- JN.default(object,
                      theta_1 = theta_1, theta_2 = theta_2,
                      theta_3 = theta_3,
                      theta_1_vals = theta_1_vals,
                      theta_2_vals = theta_2_vals,
                      theta_3_vals = theta_3_vals,
                      support = support, alpha = alpha,
                      control_fdr = control_fdr, range_size = range_size,
                      thresholds = thresholds, ...)
  if (fixed_only) return(fixed)

  groups <- .jn_lme4_groups(object, fixed, group_var)
  structure(c(list(fixed = fixed), groups),
            class = c("JN_list", "list"),
            n_way = fixed$n_way)
}

#' @rdname JN
#' @export
JN.lmerModLmerTest <- JN.merMod


# ---- per-group analyses for lme4 --------------------------------------------

# The covariance matrix is the one of the fixed effects throughout: conditional
# modes are predictions, not estimates, and lme4 offers no sampling covariance
# for them. Group-specific coefficients shift the conditional effect but not
# its nominal uncertainty, which is what the plots then show.
.jn_lme4_groups <- function(object, fixed, group_var) {
  if (!requireNamespace("lme4", quietly = TRUE))
    stop("Package 'lme4' is required for mixed-model input.", call. = FALSE)

  if (is.null(group_var)) {
    group_var <- names(lme4::ranef(object))[1]
    message("Using grouping variable: ", group_var)
  }

  gc_all <- stats::coef(object)[[group_var]]
  cn_all <- colnames(gc_all)
  vars   <- fixed$labels
  wanted <- .jn_par_index(names(lme4::fixef(object)), vars)

  missing_cn <- wanted[!wanted %in% cn_all]
  if (length(missing_cn))
    warning("Not all interaction terms vary by ", group_var,
            "; using the fixed effect for: ",
            paste(missing_cn, collapse = ", "), call. = FALSE)

  stats::setNames(
    lapply(rownames(gc_all), function(g) {
      b <- vapply(seq_along(wanted), function(j) {
        if (wanted[j] %in% cn_all) gc_all[g, wanted[j]]
        else fixed$input$coefficients[j]
      }, numeric(1))

      inp <- jn_wald(coefficients = b,
                     vcov         = fixed$input$vcov,
                     labels       = vars,
                     ranges       = fixed$input$ranges,
                     support      = fixed$input$support)
      res <- .jn_build(inp,
                       ranges      = vector("list", length(vars)),
                       range_size  = length(fixed$grid[[1]]),
                       alpha       = fixed$alpha,
                       control_fdr = fixed$control_fdr,
                       thresholds  = fixed$thresholds,
                       group       = g)
      res
    }),
    rownames(gc_all))
}


# ---- helpers ----------------------------------------------------------------

.jn_ranges_arg <- function(v1, v2, v3, theta_3) {
  if (is.null(theta_3)) list(v1, v2) else list(v1, v2, v3)
}
