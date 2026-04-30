#' Create Johnson-Neyman plots for Bayesian models
#'
#' Accepts either a raw matrix of posterior draws (rows = iterations,
#' columns = parameters) or a \code{multiSiena} object produced by
#' \code{sienaBayes()}.
#'
#' @param x matrix or 'multiSiena'; posterior draws or sienaBayes() output.
#'   If a matrix, columns should be named or referenced by index.
#'   If a 'multiSiena' object, parameters are referenced by their position
#'   in the (rate-excluded) effects object.
#' @param theta_1 character or numeric; name/index of the first variable
#'   involved in the interaction. (For \code{multiSiena} input, this is the 
#'   position of the effect in the effects object ignoring rates -- \code{x[x$type != 'rate', ]}).
#' @param theta_2 character or numeric; name/index of the second variable. 
#' @param theta_3 character or numeric; name/index of the third variable.
#'   Default NULL (two-way interaction).
#' @param theta_int_12 numeric; index of the interaction between theta_1
#'   and theta_2. Not needed when theta_1 etc. are character names
#'   (matrix input and \code{multiSiena} only).
#' @param theta_int_13 numeric; index of the theta_1:theta_3 interaction.
#'   Default NULL.
#' @param theta_int_23 numeric; index of the theta_2:theta_3 interaction.
#'   Default NULL.
#' @param theta_int_123 numeric; index of the three-way interaction.
#'   Default NULL.
#' @param theta_1_vals numeric; range of the statistic theta_1 is
#'   multiplied with.
#' @param theta_2_vals numeric; range of the statistic theta_2 is
#'   multiplied with.
#' @param theta_3_vals numeric; range of the statistic theta_3 is
#'   multiplied with. Default NULL.
#' @param burn_in numeric; burn-in iterations to remove.
#'   For \code{multiSiena} input defaults to \code{max(x$nwarm, 1)};
#'   for matrix input defaults to 0.
#' @param thin numeric; thinning interval. Default 1.
#' @param thresholds numeric; threshold for significance hashing.
#'   Default \code{c(0.49999999999999999, 0.5)}.
#' @param hyper_only logical; (\code{multiSiena} only) use only the
#'   hyper-parameter, or also produce group-level plots? Default TRUE.
#' @param round_res numeric; rounding digits. Default 3.
#' @param noTitle character; optional plot title.
#' @param color_mid character; mid-point colour. Default '#EBCC2A'.
#' @param color_low character; low-value colour. Default '#3B9AB2'.
#' @param color_high character; high-value colour. Default '#F21A00'.
#' @param color_values character; number colour. Default 'grey40'.
#' @param color_grid character; grid colour. Default 'black'.
#' @param grid_density numeric; hash-grid density. Default 0.01.
#' @param grid_spacing numeric; hash-grid spacing. Default 0.1.
#' @param save logical; save plots with ggsave()? Default FALSE.
#' @param folder character; save folder. Default NULL, which writes into a
#'   session-temporary directory
#'   (\code{file.path(tempdir(), 'int3ract JNKplots')}). Set explicitly to
#'   write elsewhere.
#'
#' @returns A list containing tables and plots. For two-way interactions:
#'   \code{param_table} and \code{plots}. For three-way: \code{thetas},
#'   \code{standard_errors}, \code{p_values}, \code{significance}, and
#'   \code{plots}. When \code{hyper_only = FALSE} ('multiSiena'), also returns a
#'   list of group-level results under \code{random_groups_effects}.
#' @export
#'
#' @examples
#' # --- two-way: raw posterior matrix (fast, no extra packages) ---
#' set.seed(1)
#' n_iter <- 500
#' post2 <- cbind(x     = rnorm(n_iter,  0.5, 0.2),
#'                z     = rnorm(n_iter, -0.3, 0.2),
#'                `x:z` = rnorm(n_iter,  0.4, 0.2))
#'
#' jnk_bayes2 <- JNK_bayes(post2,
#'                         theta_1 = 'x', theta_2 = 'z',
#'                         theta_1_vals = seq(-2, 2, 1),
#'                         theta_2_vals = seq(-2, 2, 1))
#'
#' # --- three-way: raw posterior matrix ---
#' post3 <- cbind(x       = rnorm(n_iter,  0.5, 0.2),
#'                z       = rnorm(n_iter, -0.3, 0.2),
#'                w       = rnorm(n_iter,  0.2, 0.2),
#'                `x:z`   = rnorm(n_iter,  0.4, 0.2),
#'                `x:w`   = rnorm(n_iter,  0.1, 0.2),
#'                `z:w`   = rnorm(n_iter, -0.1, 0.2),
#'                `x:z:w` = rnorm(n_iter,  0.2, 0.2))
#'
#' jnk_bayes3 <- JNK_bayes(post3,
#'                         theta_1 = 'x', theta_2 = 'z', theta_3 = 'w',
#'                         theta_1_vals = seq(-2, 2, 1),
#'                         theta_2_vals = seq(-2, 2, 1),
#'                         theta_3_vals = seq(-2, 2, 1))
#'
#' # --- two-way: integration with MCMCpack (only if installed) ---
#' \donttest{
#' if (requireNamespace("MCMCpack", quietly = TRUE)) {
#'   set.seed(1402)
#'   dat   <- data.frame(x = rnorm(100), z = rnorm(100))
#'   dat$y <- dat$x + 0.5 * dat$x * dat$z - 0.5 * dat$z + rnorm(100, sd = 4)
#'   mod_bayes2 <- MCMCpack::MCMCregress(y ~ x * z, data = dat,
#'                                       burnin = 1000, mcmc = 10000,
#'                                       thin = 1, verbose = 0)
#'   jnk_bayes2_e <- JNK_bayes(mod_bayes2, theta_1 = 'x', theta_2 = 'z',
#'                             theta_1_vals = seq(-3, 3, 0.5),
#'                             theta_2_vals = seq(-3, 3, 0.5))
#' }
#' }
#'
JNK_bayes <- function(x,
                      theta_1,
                      theta_2,
                      theta_3 = NULL,
                      theta_int_12 = NULL,
                      theta_int_13 = NULL,
                      theta_int_23 = NULL,
                      theta_int_123 = NULL,
                      theta_1_vals,
                      theta_2_vals,
                      theta_3_vals = NULL,
                      burn_in = NULL,
                      thin = 1,
                      thresholds = NULL,
                      hyper_only = TRUE,
                      round_res = 3,
                      noTitle = NULL,
                      color_mid = '#EBCC2A',
                      color_low = '#3B9AB2',
                      color_high = '#F21A00',
                      color_values = 'grey40',
                      color_grid = 'black',
                      grid_density = 0.01,
                      grid_spacing = 0.1,
                      save = FALSE,
                      folder = NULL) {
  
  threeWay <- !is.null(theta_3)
  if (is.null(thresholds)) thresholds <- c(0.49999999999999999, 0.5)
  folder <- folder %||% file.path(tempdir(), 'int3ract JNKplots')
  
  call_support <- function(theta, group) {
    if (threeWay) {
      do.call(jnb_support3,
              list(theta = theta,
                   group = group,
                   theta_1 = theta_1,
                   theta_2 = theta_2,
                   theta_3 = theta_3,
                   theta_1n = theta_1n,
                   theta_2n = theta_2n,
                   theta_3n = theta_3n,
                   theta_int_12 = theta_int_12,
                   theta_int_13 = theta_int_13,
                   theta_int_23 = theta_int_23,
                   theta_int_123 = theta_int_123,
                   theta_1_vals = theta_1_vals,
                   theta_2_vals = theta_2_vals,
                   theta_3_vals = theta_3_vals,
                   thresholds = thresholds,
                   color_mid = color_mid,
                   color_low = color_low,
                   color_high = color_high,
                   color_values = color_values,
                   color_grid = color_grid,
                   grid_density = grid_density,
                   grid_spacing = grid_spacing,
                   save = save,
                   folder = folder))
    } else {
      do.call(jnb_support2,
              list(theta = theta,
                   group = group,
                   theta_1 = theta_1,
                   theta_2 = theta_2,
                   theta_1n = theta_1n,
                   theta_2n = theta_2n,
                   theta_int_12 = theta_int_12,
                   theta_1_vals = theta_1_vals,
                   theta_2_vals = theta_2_vals,
                   color_mid = color_mid,
                   color_low = color_low,
                   color_high = color_high,
                   save = save,
                   folder = folder))
    }
  }
  
  if (inherits(x, 'multiSiena') || inherits(x, 'sienaBayesFit')) {
    if (is.null(burn_in)) burn_in <- max(x$nwarm, 1)

    eff <-  x$effects[!x$basicRate, ]
    rownames(eff) <- seq_len(nrow(eff))
    
    theta_1n <- .clean_effect_name(eff$effectName[theta_1])
    theta_2n <- .clean_effect_name(eff$effectName[theta_2])
    theta_3n <- if (threeWay) .clean_effect_name(eff$effectName[theta_3]) else NULL
    
    idx <- seq(burn_in, nrow(x$ThinParameters), thin)
    
    thetas <- if (threeWay) {
      c(theta_1, theta_2, theta_3, theta_int_12,
        theta_int_13, theta_int_23, theta_int_123)
    } else {
      c(theta_1, theta_2, theta_int_12)
    }
    
    any_random <- any(eff$randomEffects[thetas])
    
    if (!any_random) {
      theta <- x$ThinParameters[idx, 1, thetas]
      return(call_support(theta, 'Eta'))
    }
    
    extract_theta <- function(thetas) {
      eff_ran <- eff[eff$randomEffects, ]
      t(sapply(thetas, function(i) {
        if (eff$randomEffects[i]) {
          x$ThinPosteriorMu[idx, which(as.numeric(rownames(eff_ran)) == i) + 1]
        } else {
          x$ThinParameters[idx, 1, i + 1]
        }
      })) |> t()
    }
    
    theta      <- extract_theta(thetas)
    returnList <- call_support(theta, 'Mu')
    
    if (!hyper_only) {
      randomResults <- setNames(
        lapply(seq_len(x$nGroup), function(i) {
          theta_g <- x$ThinParameters[idx, i, thetas]
          call_support(theta_g, paste0('group', i))
        }),
        paste0('group', seq_len(x$nGroup))
      )
      returnList <- setNames(list(returnList, randomResults),
                             c('Mu', 'random_groups_effects'))
    }
    
    return(returnList)
  }
  
  # ---- matrix path ----
  if (any(c(is.character(theta_1),
            is.character(theta_2),
            is.character(theta_3))) &
      any(c(is.numeric(theta_1),
            is.numeric(theta_2),
            is.numeric(theta_3)))) {
    stop("Please only provide either numbers or characters to theta_1[23] etc.")
  }
  
  if (is.character(theta_1)) {
    theta_int_12 <- paste0(theta_1, ':', theta_2)
    if (threeWay) {
      theta_int_13  <- paste0(theta_1, ':', theta_3)
      theta_int_23  <- paste0(theta_2, ':', theta_3)
      theta_int_123 <- paste0(theta_1, ':', theta_2, ':', theta_3)
    }
  }
  
  if (is.null(burn_in)) burn_in <- 0
  
  theta <- x[seq(burn_in + 1, nrow(x), thin),
             c(theta_1, theta_2, theta_3,
               theta_int_12, theta_int_13,
               theta_int_23, theta_int_123)]
  
  if (is.numeric(theta_1)) {
    colnames(theta) <- if (threeWay) {
      c('V1', 'V2', 'V3', 'V1:V2', 'V1:V3', 'V2:V3', 'V1:V2:V3')
    } else {
      c('V1', 'V2', 'V1:V2')
    }
  }
  
  theta_1n <- colnames(theta)[1]
  theta_2n <- colnames(theta)[2]
  theta_3n <- if (threeWay) colnames(theta)[3] else NULL
  
  call_support(theta, NULL)
}
