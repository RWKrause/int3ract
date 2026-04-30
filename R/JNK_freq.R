#' Johnson-Neyman(-Krause) plots for frequentist models
#'
#' Unified function that accepts \code{lm}, \code{glm}, \code{sienaFit}
#' ('RSiena'), or \code{lmerMod}/\code{glmerMod} ('lme4') objects, or raw
#' coefficient vectors and covariance matrices.
#'
#' @param x model object (\code{lm}, \code{glm}, \code{sienaFit},
#'   \code{lmerMod}, \code{glmerMod}) or \code{NULL} when supplying
#'   \code{covar}/\code{coefs} directly.
#' @param theta_1 character or numeric; name (lm/glm/lmerMod/glmerMod) or index (sienaFit)
#'   of the first variable.
#' @param theta_2 character or numeric; second variable.
#' @param theta_3 character or numeric; third variable. Default NULL (two-way).
#' @param theta_int_12 numeric; index of the theta_1:theta_2 interaction
#'   (sienaFit / generic only). For m/glm/lmerMod/glmerMod inputs the interaction name is
#'   resolved automatically.
#' @param theta_int_13 numeric; index of the theta_1:theta_3 interaction.
#'   Default NULL.
#' @param theta_int_23 numeric; index of the theta_2:theta_3 interaction.
#'   Default NULL.
#' @param theta_int_123 numeric; index of the three-way interaction.
#'   Default NULL.
#' @param theta_1_vals numeric; evaluation range for theta_1.
#'   Auto-derived from model data for lm/glm/lmerMod/glmerMod if NULL.
#' @param theta_2_vals numeric; evaluation range for theta_2.
#' @param theta_3_vals numeric; evaluation range for theta_3. Default NULL.
#' @param covar matrix; covariance matrix of the relevant parameters.
#'   Required only when \code{x} is NULL.
#' @param coefs numeric; coefficient vector.
#'   Required only when \code{x} is NULL.
#' @param name character; variable names.
#'   Required only when \code{x} is NULL.
#' @param group_var character; ('lme4' only) grouping variable for random
#'   effects. Defaults to the first grouping factor.
#' @param fixed_only logical; ('lme4' only) produce only fixed-effects plots?
#'   If FALSE, per-group plots are produced for groups with random interaction
#'   terms. Default TRUE.
#' @param control_fdr logical; apply Bonferroni-Holm correction? Default FALSE.
#' @param alpha numeric; significance level. Default 0.05.
#' @param round_res integer; rounding precision. Default 3.
#' @param range_size integer; number of moderator values. Default 1000 for
#'   two-way, 50 for three-way.
#' @param sig_color character; significant region colour (2-way). Default
#'   'seagreen3'.
#' @param non_sig_color character; non-significant region colour (2-way).
#'   Default 'chocolate'.
#' @param line_color character; line colour (2-way). Default 'black'.
#' @param color_mid character; midpoint colour (3-way heatmap). Default
#'   '#EBCC2A'.
#' @param color_low character; low-value colour. Default '#3B9AB2'.
#' @param color_high character; high-value colour. Default '#F21A00'.
#' @param color_values character; value label colour. Default 'grey40'.
#' @param color_grid character; crosshatch colour. Default 'black'.
#' @param grid_density numeric; crosshatch density. Default 0.01.
#' @param grid_spacing numeric; crosshatch spacing. Default 0.1.
#' @param crosshatch_non_sig logical; crosshatch non-significant cells?
#'   Default TRUE.
#' @param save logical; save plots via ggsave()? Default FALSE.
#' @param folder character; output folder for saved plots. Default NULL,
#'   which writes into a session-temporary directory
#'   (\code{file.path(tempdir(), 'int3ract JNKplots')}). Set explicitly to
#'   write elsewhere.
#'
#' @returns A list containing tables and plots. For two-way interactions:
#'   \code{param_table} and \code{plots}. For three-way: \code{thetas},
#'   \code{standard_errors}, \code{p_values}, \code{significance}, and
#'   \code{plots}. When \code{fixed_only = FALSE} ('lme4'), returns a list
#'   with \code{fixed} and \code{random_groups} elements.
#'
#' @import ggplot2
#' @import ggpattern
#' @import scales
#' @import tidyr
#' @import tibble
#'
#' @export
#'
#' @examples
#' # --- two-way lm ---
#' set.seed(1)
#' dat <- data.frame(y = rnorm(100), x = rnorm(100),
#'                   z = rnorm(100), w = rnorm(100))
#' res <- lm(y ~ x * z * w, dat)
#'
#' x2 <- JNK_freq(res, theta_1 = 'x', theta_2 = 'z',
#'                range_size = 100)
#' x2$plots$z
#'
#' # --- three-way lm (small grid for speed) ---
#' x3 <- JNK_freq(res, theta_1 = 'x', theta_2 = 'z', theta_3 = 'w',
#'                range_size = 10)
#' x3$plots$z
#'
#' # --- generic (covariance + coefficients) ---
#' JNK_freq(covar = vcov(res)[c('x','z','x:z'), c('x','z','x:z')],
#'          coefs = coef(res)[c('x','z','x:z')],
#'          name  = c('x', 'z'),
#'          theta_1 = 'x',
#'          theta_2 = 'z',
#'          theta_1_vals = c(-3, 3),
#'          theta_2_vals = c(-3, 3),
#'          range_size = 100)
JNK_freq <- function(x = NULL,
                     theta_1,
                     theta_2,
                     theta_3 = NULL,
                     theta_int_12 = NULL,
                     theta_int_13 = NULL,
                     theta_int_23 = NULL,
                     theta_int_123 = NULL,
                     theta_1_vals = NULL,
                     theta_2_vals = NULL,
                     theta_3_vals = NULL,
                     covar = NULL,
                     coefs = NULL,
                     name = NULL,
                     group_var = NULL,
                     fixed_only = TRUE,
                     control_fdr = FALSE,
                     alpha = 0.05,
                     round_res = 3,
                     range_size = NULL,
                     sig_color = 'seagreen3',
                     non_sig_color = 'chocolate',
                     line_color = 'black',
                     color_mid = '#EBCC2A',
                     color_low = '#3B9AB2',
                     color_high = '#F21A00',
                     color_values = 'grey40',
                     color_grid = 'black',
                     grid_density = 0.01,
                     grid_spacing = 0.1,
                     crosshatch_non_sig = TRUE,
                     save = FALSE,
                     folder = NULL) {

  threeWay <- !is.null(theta_3)
  folder   <- folder %||% file.path(tempdir(), 'int3ract JNKplots')

  if (is.null(range_size))
    range_size <- if (threeWay) 50 else 1000

  if (inherits(x, 'sienaFit')) {
    out <- .extract_siena(x, theta_1, theta_2, theta_3,
                          theta_int_12, theta_int_13,
                          theta_int_23, theta_int_123,
                          theta_1_vals, theta_2_vals, theta_3_vals)

  } else if (inherits(x, c('lm', 'glm'))) {
    out <- .extract_lm(x, theta_1, theta_2, theta_3,
                       theta_1_vals, theta_2_vals, theta_3_vals)

  } else if (inherits(x, c('lmerMod', 'glmerMod', 'lmerModLmerTest'))) {
    if (!requireNamespace("lme4", quietly = TRUE))
      stop("Package 'lme4' is required for lmerMod/glmerMod inputs.")
    out <- .extract_lme4(x, theta_1, theta_2, theta_3,
                         theta_1_vals, theta_2_vals, theta_3_vals)

  } else if (!is.null(covar) && !is.null(coefs)) {
    out <- list(
      name  = name %||% (if (threeWay) c(theta_1, theta_2, theta_3)
                         else c(theta_1, theta_2)),
      coefs = coefs,
      covar = covar,
      vals  = if (threeWay) list(theta_1_vals, theta_2_vals, theta_3_vals)
              else list(theta_1_vals, theta_2_vals)
    )

  } else {
    stop("x must be an lm, glm, sienaFit, lmerMod, or glmerMod object,\n",
         "or supply covar and coefs directly.")
  }

  plot_args <- list(alpha = alpha, control_fdr = control_fdr,
                    round_res = round_res, range_size = range_size,
                    sig_color = sig_color, non_sig_color = non_sig_color,
                    line_color = line_color,
                    color_mid = color_mid, color_low = color_low,
                    color_high = color_high, color_values = color_values,
                    color_grid = color_grid, grid_density = grid_density,
                    grid_spacing = grid_spacing,
                    crosshatch_non_sig = crosshatch_non_sig,
                    save = save, folder = folder)

  result <- do.call(.freq_core,
                    c(list(name = out$name, covar = out$covar,
                           coefs = out$coefs, vals = out$vals,
                           group = NULL),
                      plot_args))

  if (inherits(x, c('lmerMod', 'glmerMod', 'lmerModLmerTest')) && !fixed_only) {
    grp_result <- .lme4_groups(x, theta_1, theta_2, theta_3,
                               out, group_var, plot_args)
    result <- list(fixed = result, random_groups = grp_result)
  }

  result
}
