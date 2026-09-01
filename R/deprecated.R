# =============================================================================
# The 1.0.x interface, kept working.
#
# JNK_freq() and JNK_bayes() split the package by inference paradigm rather
# than by input object, were not generic, and returned unclassed lists. JN()
# replaces all three properties. These wrappers delegate to it and reshape the
# result into the old layout so that scripts written against 1.0.x -- including
# the replication material of the preprint -- keep running.
# =============================================================================

.jn_deprecation_warned <- new.env(parent = emptyenv())

.jn_deprecate <- function(old, new) {
  if (isTRUE(.jn_deprecation_warned[[old]])) return(invisible(NULL))
  .jn_deprecation_warned[[old]] <- TRUE
  warning(old, "() is deprecated and will be removed in a future release; ",
          "use ", new, "() instead. It dispatches on the fitted object and ",
          "returns a classed result with print(), summary() and plot() ",
          "methods. See ?JN.", call. = FALSE)
}


#' Johnson-Neyman plots for frequentist models (deprecated)
#'
#' Superseded by \code{\link{JN}}, which dispatches on the fitted object rather
#' than on the inference paradigm and returns a classed object with
#' \code{print()}, \code{summary()} and \code{plot()} methods. This wrapper
#' delegates to \code{JN()} and reshapes the result into the layout returned by
#' \pkg{int3ract} 1.0.x.
#'
#' Note that the two-way standard errors returned here differ from those of
#' \pkg{int3ract} 1.0.x, which used the covariance between the two main effects
#' where the delta method calls for the covariance between the focal main
#' effect and the interaction. The values are now correct; see \code{NEWS.md}.
#'
#' Figures are drawn without the data-density panels that \code{\link{plot.JN}}
#' adds, so that they remain plain \pkg{ggplot2} objects as before.
#'
#' @param x a fitted model, or \code{NULL} when supplying \code{covar} and
#'   \code{coefs} directly.
#' @param theta_1,theta_2,theta_3 the variables involved in the interaction.
#' @param theta_int_12,theta_int_13,theta_int_23,theta_int_123 interaction
#'   positions, for \code{sienaFit} and generic input.
#' @param theta_1_vals,theta_2_vals,theta_3_vals moderator ranges.
#' @param covar,coefs,name covariance matrix, coefficient vector and variable
#'   names for generic input.
#' @param group_var,fixed_only grouping factor and whether to restrict the
#'   analysis to the fixed effects (\pkg{lme4} only).
#' @param control_fdr,alpha,round_res,range_size analysis settings.
#' @param sig_color,non_sig_color,line_color,color_mid,color_low,color_high
#'   colour settings, passed to \code{\link{jn_style}}.
#' @param color_values,color_grid,grid_density,grid_spacing,crosshatch_non_sig
#'   further appearance settings, passed to \code{\link{jn_style}}.
#' @param save,folder write the figures to disk?
#'
#' @returns A list, in the layout used by \pkg{int3ract} 1.0.x.
#'
#' @seealso \code{\link{JN}}
#' @keywords internal
#' @export
JNK_freq <- function(x = NULL, theta_1, theta_2, theta_3 = NULL,
                     theta_int_12 = NULL, theta_int_13 = NULL,
                     theta_int_23 = NULL, theta_int_123 = NULL,
                     theta_1_vals = NULL, theta_2_vals = NULL,
                     theta_3_vals = NULL,
                     covar = NULL, coefs = NULL, name = NULL,
                     group_var = NULL, fixed_only = TRUE,
                     control_fdr = FALSE, alpha = 0.05, round_res = 3,
                     range_size = NULL,
                     sig_color = "seagreen3", non_sig_color = "chocolate",
                     line_color = "black", color_mid = "#EBCC2A",
                     color_low = "#3B9AB2", color_high = "#F21A00",
                     color_values = "grey40", color_grid = "black",
                     grid_density = 0.01, grid_spacing = 0.1,
                     crosshatch_non_sig = TRUE,
                     save = FALSE, folder = NULL) {
  .jn_deprecate("JNK_freq", "JN")

  threeWay <- !is.null(theta_3)
  style <- jn_style(sig_color = sig_color, non_sig_color = non_sig_color,
                    line_color = line_color, color_mid = color_mid,
                    color_low = color_low, color_high = color_high,
                    color_values = color_values, color_grid = color_grid,
                    grid_density = grid_density, grid_spacing = grid_spacing,
                    crosshatch_non_sig = crosshatch_non_sig,
                    show_density = FALSE)

  args <- list(theta_1 = theta_1, theta_2 = theta_2, theta_3 = theta_3,
               theta_1_vals = theta_1_vals, theta_2_vals = theta_2_vals,
               theta_3_vals = theta_3_vals,
               alpha = alpha, control_fdr = control_fdr,
               range_size = range_size)

  res <- if (is.null(x) && !is.null(covar) && !is.null(coefs)) {
    labels <- name %||% (if (threeWay) c(theta_1, theta_2, theta_3)
                         else c(theta_1, theta_2))
    inp <- jn_wald(coefficients = coefs, vcov = covar, labels = labels)
    do.call(JN, c(list(inp), args[!names(args) %in%
                                   c("theta_1", "theta_2", "theta_3")]))
  } else if (inherits(x, "sienaFit")) {
    do.call(JN, c(list(x), args,
                  list(theta_int_12 = theta_int_12,
                       theta_int_13 = theta_int_13,
                       theta_int_23 = theta_int_23,
                       theta_int_123 = theta_int_123)))
  } else if (inherits(x, c("lmerMod", "glmerMod", "lmerModLmerTest"))) {
    do.call(JN, c(list(x), args,
                  list(fixed_only = fixed_only, group_var = group_var)))
  } else {
    do.call(JN, c(list(x), args))
  }

  out <- if (inherits(res, "JN_list")) {
    list(fixed = .jn_legacy_freq(res$fixed, style, round_res, save, folder),
         random_groups = lapply(res[-1], .jn_legacy_freq,
                                style = style, round_res = round_res,
                                save = save, folder = folder))
  } else {
    .jn_legacy_freq(res, style, round_res, save, folder)
  }
  out
}


.jn_legacy_freq <- function(res, style, round_res, save, folder) {
  figs <- jn_plots(res, style = style)
  if (save) jn_save(res, folder = folder, style = style)

  if (res$n_way == 2L) {
    tabs <- lapply(seq_along(res$labels), function(i) {
      tb <- res$tables[[i]]
      data.frame(theta      = tb$focal,
                 moderator  = tb$moderator,
                 mod_vals   = round(tb$mod_value, round_res),
                 theta_vals = round(tb$estimate,  round_res),
                 theta_se   = round(tb$std.error, round_res),
                 theta_p    = round(tb$p.value,   3),
                 sig        = tb$significant,
                 stringsAsFactors = FALSE)
    })
    return(list(param_table = stats::setNames(tabs, res$labels),
                plots       = figs))
  }

  mats <- function(column) {
    stats::setNames(lapply(res$tables, .jn_long_to_matrix, column = column),
                    res$labels)
  }
  list(thetas          = mats("estimate"),
       standard_errors = mats("std.error"),
       p_values        = mats("p.value"),
       significance    = mats("significant"),
       plots           = figs)
}

# Long grid back to the rows-by-columns matrix the 1.0.x interface returned.
.jn_long_to_matrix <- function(tab, column) {
  r <- sort(unique(tab$mod1_value))
  c_ <- sort(unique(tab$mod2_value))
  m <- matrix(tab[[column]][order(match(tab$mod2_value, c_),
                                  match(tab$mod1_value, r))],
              nrow = length(r), ncol = length(c_),
              dimnames = list(r, c_))
  m
}


#' Johnson-Neyman plots for Bayesian models (deprecated)
#'
#' Superseded by \code{\link{JN}}, which dispatches on the fitted object rather
#' than on the inference paradigm and returns a classed object with
#' \code{print()}, \code{summary()} and \code{plot()} methods. This wrapper
#' delegates to \code{JN()} and reshapes the result into the layout returned by
#' \pkg{int3ract} 1.0.x.
#'
#' @param x a matrix of posterior draws, or a \code{multiSiena} object.
#' @param theta_1,theta_2,theta_3 the variables involved in the interaction.
#' @param theta_int_12,theta_int_13,theta_int_23,theta_int_123 interaction
#'   positions.
#' @param theta_1_vals,theta_2_vals,theta_3_vals moderator values.
#' @param burn_in,thin burn-in and thinning.
#' @param thresholds Bayesian \eqn{p} values bounding the inconclusive region.
#' @param hyper_only analyse only the hyper-parameter (\code{multiSiena})?
#' @param round_res rounding applied to the returned tables.
#' @param noTitle ignored; retained for compatibility.
#' @param color_mid,color_low,color_high,color_values,color_grid colour
#'   settings, passed to \code{\link{jn_style}}.
#' @param grid_density,grid_spacing crosshatch settings.
#' @param save,folder write the figures to disk?
#'
#' @returns A list, in the layout used by \pkg{int3ract} 1.0.x.
#'
#' @seealso \code{\link{JN}}
#' @keywords internal
#' @export
JNK_bayes <- function(x, theta_1, theta_2, theta_3 = NULL,
                      theta_int_12 = NULL, theta_int_13 = NULL,
                      theta_int_23 = NULL, theta_int_123 = NULL,
                      theta_1_vals, theta_2_vals, theta_3_vals = NULL,
                      burn_in = NULL, thin = 1, thresholds = NULL,
                      hyper_only = TRUE, round_res = 3, noTitle = NULL,
                      color_mid = "#EBCC2A", color_low = "#3B9AB2",
                      color_high = "#F21A00", color_values = "grey40",
                      color_grid = "black", grid_density = 0.01,
                      grid_spacing = 0.1, save = FALSE, folder = NULL) {
  .jn_deprecate("JNK_bayes", "JN")

  if (is.null(thresholds)) thresholds <- c(0.49999999999999999, 0.5)
  style <- jn_style(color_mid = color_mid, color_low = color_low,
                    color_high = color_high, color_values = color_values,
                    color_grid = color_grid, grid_density = grid_density,
                    grid_spacing = grid_spacing, show_density = FALSE)

  args <- list(theta_1 = theta_1, theta_2 = theta_2, theta_3 = theta_3,
               theta_1_vals = theta_1_vals, theta_2_vals = theta_2_vals,
               theta_3_vals = theta_3_vals,
               thresholds = thresholds, thin = thin)
  if (!is.null(burn_in)) args$burn_in <- burn_in

  res <- if (inherits(x, c("multiSiena", "sienaBayesFit"))) {
    do.call(JN, c(list(x), args,
                  list(theta_int_12 = theta_int_12,
                       theta_int_13 = theta_int_13,
                       theta_int_23 = theta_int_23,
                       theta_int_123 = theta_int_123,
                       hyper_only = hyper_only)))
  } else {
    do.call(JN, c(list(x), args,
                  list(theta_int_12 = theta_int_12,
                       theta_int_13 = theta_int_13,
                       theta_int_23 = theta_int_23,
                       theta_int_123 = theta_int_123)))
  }

  if (inherits(res, "JN_list"))
    return(stats::setNames(
      list(.jn_legacy_bayes(res[[1]], style, round_res, save, folder),
           lapply(res[-1], .jn_legacy_bayes, style = style,
                  round_res = round_res, save = save, folder = folder)),
      c("Mu", "random_groups_effects")))

  .jn_legacy_bayes(res, style, round_res, save, folder)
}


.jn_legacy_bayes <- function(res, style, round_res, save, folder) {
  if (save) jn_save(res, folder = folder, style = style)

  if (res$n_way == 2L) {
    tabs <- lapply(seq_along(res$labels), function(i) {
      tb <- res$tables[[i]]
      data.frame(theta         = tb$focal,
                 moderator     = tb$moderator,
                 modValue      = as.factor(round(tb$mod_value, 3)),
                 thetaPostMean = tb$post_mean,
                 thetaPostSD   = tb$post_sd,
                 bayes_p       = tb$bayes_p,
                 thetaPost2.5  = tb$q2.5,
                 thetaPost97.5 = tb$q97.5,
                 stringsAsFactors = FALSE)
    })
    figs <- lapply(seq_along(res$labels),
                   function(i) plot(res, which = i, type = "density",
                                    style = style))
    return(stats::setNames(c(tabs, figs),
                           c(paste0(res$labels, "_table"),
                             paste0(res$labels, "_plot"))))
  }

  tabs <- stats::setNames(lapply(res$tables, function(tb)
    data.frame(theta         = tb$focal,
               mod1          = tb$mod1,
               mod2          = tb$mod2,
               mod1Val       = tb$mod1_value,
               mod2Val       = tb$mod2_value,
               thetaPostMean = tb$post_mean,
               thetaPostSD   = tb$post_sd,
               bayes_p       = tb$bayes_p,
               thetaPost2.5  = tb$q2.5,
               thetaPost97.5 = tb$q97.5,
               sig           = tb$significant,
               pattern       = ifelse(!tb$significant, "crosshatch", "none"),
               stringsAsFactors = FALSE)),
    res$labels)

  figs <- jn_plots(res, style = style)
  list(result_tables   = tabs,
       post_mean_plots = figs$post_mean,
       bayes_p_plots   = figs$bayes_p)
}
