# =============================================================================
# Bayesian stochastic actor-oriented models (multiSiena / sienaBayesFit).
#
# Effects are addressed by their position in the rate-excluded effects object,
# x$effects[!x$basicRate, ]. Whether a parameter was estimated as fixed across
# groups (eta) or as varying between them (mu) decides which set of draws is
# the right one, and that is worked out here rather than left to the user.
# =============================================================================

.clean_effect_name <- function(n) {
  gsub(":", "", gsub("/", "o", gsub("\\^", "", n)))
}

.jn_siena_effects <- function(object) {
  eff <- object$effects[!object$basicRate, ]
  rownames(eff) <- seq_len(nrow(eff))
  eff
}

.jn_siena_positions <- function(object, theta_1, theta_2, theta_3,
                                theta_int_12, theta_int_13,
                                theta_int_23, theta_int_123) {
  eff      <- .jn_siena_effects(object)
  threeWay <- !is.null(theta_3)

  pos <- if (threeWay) {
    c(theta_1, theta_2, theta_3, theta_int_12,
      theta_int_13, theta_int_23, theta_int_123)
  } else {
    c(theta_1, theta_2, theta_int_12)
  }
  nms <- if (threeWay) {
    c("theta_1", "theta_2", "theta_3", "theta_int_12",
      "theta_int_13", "theta_int_23", "theta_int_123")
  } else {
    c("theta_1", "theta_2", "theta_int_12")
  }

  overview <- function(...) {
    stop(..., "\n",
         "Effects are addressed by their position in the rate-excluded ",
         "effects object, x$effects[!x$basicRate, ]. This model has ",
         nrow(eff), " such effects:\n",
         paste0("  ", format(seq_len(nrow(eff))), ": ", eff$effectName,
                collapse = "\n"),
         call. = FALSE)
  }

  if (length(pos) != length(nms))
    overview("Every interaction position must be supplied: ",
             paste(nms, collapse = ", "), ".")
  if (!is.numeric(pos))
    overview("For multiSiena and sienaBayesFit input, the parameters must be ",
             "given as numeric positions.")
  bad <- is.na(pos) | pos < 1 | pos > nrow(eff) | pos != round(pos)
  if (any(bad))
    overview("Effect position out of range: ",
             paste0(nms[bad], " = ", pos[bad], collapse = ", "), ".")

  list(pos = pos, eff = eff)
}


#' @rdname jn_input
#' @param group for \code{sienaBayesFit} input, the index of the group whose
#'   draws to extract. \code{NULL} (the default) uses the hyper-parameter, or
#'   the fixed parameter when none of the effects varies between groups.
#' @export
jn_input.sienaBayesFit <- function(object, theta_1, theta_2, theta_3 = NULL,
                                   theta_int_12 = NULL, theta_int_13 = NULL,
                                   theta_int_23 = NULL, theta_int_123 = NULL,
                                   burn_in = NULL, thin = 1, group = NULL,
                                   ranges = NULL, support = NULL, ...) {
  p   <- .jn_siena_positions(object, theta_1, theta_2, theta_3,
                             theta_int_12, theta_int_13,
                             theta_int_23, theta_int_123)
  pos <- p$pos
  eff <- p$eff

  if (is.null(burn_in)) burn_in <- max(object$nwarm, 1)

  # ThinParameters and ThinPosteriorMu carry each group's basic rate parameters
  # before the remaining effects, so effect i sits at i + n_rate.
  n_rate <- sum(object$basicRate) / object$nGroup
  if (n_rate != round(n_rate) ||
      dim(object$ThinParameters)[3] != n_rate + nrow(eff))
    stop("Cannot align the effects object with x$ThinParameters: expected ",
         "x$nGroup (", object$nGroup, ") blocks of ",
         dim(object$ThinParameters)[3], " parameters made up of the basic ",
         "rates (", sum(object$basicRate), " in total) plus ", nrow(eff),
         " further effects.", call. = FALSE)

  idx    <- seq(burn_in, nrow(object$ThinParameters), by = thin)
  labels <- .clean_effect_name(
    eff$effectName[if (is.null(theta_3)) c(theta_1, theta_2)
                   else c(theta_1, theta_2, theta_3)])

  # The carrier records which level of the model its draws came from, so that
  # the analysis can say so: eta is the parameter shared across groups, mu the
  # mean of a parameter that varies between them, and the rest are the
  # group-specific draws themselves.
  if (!is.null(group)) {
    draws <- object$ThinParameters[idx, group, pos + n_rate]
    return(structure(jn_posterior(draws, labels, ranges, support),
                     jn_level = paste0("group", group)))
  }

  if (!any(eff$randomEffects[pos])) {
    draws <- object$ThinParameters[idx, 1, pos + n_rate]
    return(structure(jn_posterior(draws, labels, ranges, support),
                     jn_level = "Eta"))
  }

  if (ncol(object$ThinPosteriorMu) != n_rate + sum(eff$randomEffects))
    stop("Cannot align the effects object with x$ThinPosteriorMu: expected ",
         n_rate, " basic rate column(s) plus ", sum(eff$randomEffects),
         " random effects, but found ", ncol(object$ThinPosteriorMu),
         " columns.", call. = FALSE)

  eff_ran <- eff[eff$randomEffects, ]
  draws   <- vapply(pos, function(i) {
    if (eff$randomEffects[i])
      object$ThinPosteriorMu[idx,
        which(as.numeric(rownames(eff_ran)) == i) + n_rate]
    else
      object$ThinParameters[idx, 1, i + n_rate]
  }, numeric(length(idx)))

  structure(jn_posterior(draws, labels, ranges, support), jn_level = "Mu")
}

#' @rdname jn_input
#' @export
jn_input.multiSiena <- jn_input.sienaBayesFit


#' @rdname JN
#' @param hyper_only analyse only the population-level parameter? When
#'   \code{FALSE}, a separate analysis is returned for each group as well, in
#'   the same way that \code{fixed_only = FALSE} works for \pkg{lme4} models.
#'   Group-level analyses are only produced when at least one of the parameters
#'   involved actually varies between groups; when none does, there is nothing
#'   for them to show and the population-level analysis is returned alone.
#' @export
JN.sienaBayesFit <- function(object,
                             theta_1, theta_2, theta_3 = NULL,
                             theta_1_vals = NULL, theta_2_vals = NULL,
                             theta_3_vals = NULL,
                             support = NULL,
                             alpha = 0.05, control_fdr = FALSE,
                             range_size = NULL, thresholds = NULL,
                             hyper_only = TRUE, ...) {

  ranges     <- .jn_ranges_arg(theta_1_vals, theta_2_vals, theta_3_vals,
                               theta_3)
  thresholds <- thresholds %||% c(alpha / 2, 1 - alpha / 2)

  build <- function(input, group) {
    if (!is.null(support))
      input$support <- .jn_check_support(support, length(input$labels))
    .jn_build(input, ranges = ranges, range_size = range_size,
              alpha = alpha, control_fdr = control_fdr,
              thresholds = thresholds, group = group)
  }

  input <- jn_input(object, theta_1 = theta_1, theta_2 = theta_2,
                    theta_3 = theta_3, ...)
  level <- attr(input, "jn_level")
  hyper <- build(input, level)

  if (hyper_only) return(hyper)

  # Every parameter involved is shared across groups, so the group-level
  # analyses would be nGroup copies of the same picture. Say so rather than
  # producing them.
  if (identical(level, "Eta")) {
    message("None of the parameters involved varies between groups; ",
            "returning the population-level (Eta) analysis alone.")
    return(hyper)
  }

  groups <- stats::setNames(
    lapply(seq_len(object$nGroup), function(g) {
      gi <- jn_input(object, theta_1 = theta_1, theta_2 = theta_2,
                     theta_3 = theta_3, group = g, ...)
      build(gi, attr(gi, "jn_level"))
    }),
    paste0("group", seq_len(object$nGroup)))

  structure(c(list(hyper = hyper), groups),
            class = c("JN_list", "list"),
            n_way = hyper$n_way)
}

#' @rdname JN
#' @export
JN.multiSiena <- JN.sienaBayesFit
