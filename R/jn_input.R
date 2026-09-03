# =============================================================================
# jn_input(): the extension point
#
# Every fitted object supported by int3ract is turned into one of two
# normalized carriers before any computation happens:
#
#   jn_wald()      point estimates + covariance matrix  -> Wald z tests
#   jn_posterior() a matrix of draws                -> conditional posteriors
#
# Supporting a new model class therefore means writing a single jn_input()
# method; nothing downstream needs to change.
# =============================================================================

# Coefficients and draws are always stored in the canonical order
#   two-way:   (theta_1, theta_2, theta_1:theta_2)
#   three-way: (theta_1, theta_2, theta_3, theta_1:theta_2, theta_1:theta_3,
#               theta_2:theta_3, theta_1:theta_2:theta_3)
.jn_n_par <- function(k) if (k == 2L) 3L else 7L


#' Carriers for Johnson-Neyman input
#'
#' \code{jn_input()} is the generic that turns a fitted model into the
#' normalized representation used by \code{\link{JN}}. It is the extension
#' point of the package: to support a model class that \pkg{int3ract} does not
#' know about, write a \code{jn_input()} method for it that returns either a
#' \code{jn_wald} or a \code{jn_posterior} object built with the constructors
#' documented here.
#'
#' \code{jn_wald()} carries point estimates and their covariance matrix, and
#' leads to Wald \eqn{z} tests of the conditional effect. \code{jn_posterior()}
#' carries a matrix of draws, and leads to conditional posterior distributions
#' and Bayesian \eqn{p} values. Which of the two a method returns is what makes
#' an analysis frequentist or Bayesian: the distinction is a property of the
#' input, not of the function the user calls.
#'
#' @param object a fitted model object.
#' @param theta_1,theta_2 name (or, for \code{sienaFit} and
#'   \code{sienaBayesFit}, integer position) of the first and second variable
#'   involved in the interaction.
#' @param theta_3 the third variable, or \code{NULL} (default) for a two-way
#'   interaction.
#' @param coefficients numeric vector of length 3 (two-way) or 7 (three-way),
#'   in the canonical order \code{(theta_1, theta_2, theta_1:theta_2)} or
#'   \code{(theta_1, theta_2, theta_3, theta_1:theta_2, theta_1:theta_3,
#'   theta_2:theta_3, theta_1:theta_2:theta_3)}.
#' @param vcov covariance matrix of \code{coefficients}, in the same order.
#' @param draws matrix of draws with one row per iteration and one column per
#'   parameter, in the same canonical order.
#' @param labels character vector of length 2 or 3 with display names for the
#'   focal variables.
#' @param ranges list of length 2 or 3; the range over which each variable is
#'   evaluated when it acts as a moderator. Elements may be \code{NULL} if they
#'   are to be derived from the data.
#' @param support list of length 2 or 3 holding the observed values of each
#'   variable, or \code{NULL}. Used to draw the data-density panels of
#'   \code{\link{plot.JN}}, which show how much empirical support each part of
#'   the moderator range actually has. Supplying it is optional but strongly
#'   recommended: a region of significance covering moderator values that were
#'   barely observed is not evidence of much.
#' @param theta_int_12,theta_int_13,theta_int_23,theta_int_123 the interaction
#'   parameters. Required for \code{sienaFit} and \code{sienaBayesFit} input,
#'   where effects are addressed by integer position rather than by name;
#'   optional for matrices of draws whose columns are named. For \code{lavaan}
#'   input they are term names, needed only when a product is not written in
#'   lavaan's own \code{x:z} form -- which is always the case for three-way
#'   interactions, since lavaan does not accept those in its model syntax.
#' @param outcome for \code{lavaan} input, the dependent variable whose
#'   regression is to be analysed. A structural model usually has several, and
#'   the same interaction may predict more than one of them with different
#'   coefficients. May be omitted when exactly one outcome is regressed on all
#'   of the focal variables; when several are, \code{jn_input()} stops and
#'   lists them rather than guessing.
#' @param burn_in number of initial draws to discard.
#' @param thin thinning interval applied after burn-in.
#' @param ... passed to methods; ignored by the constructors.
#'
#' @returns An object of class \code{c("jn_wald", "jn_input")} or
#'   \code{c("jn_posterior", "jn_input")}.
#'
#' @seealso \code{\link{JN}} for the analysis itself.
#'
#' @examples
#' set.seed(1)
#' dat <- data.frame(x = rnorm(100), z = rnorm(100))
#' dat$y <- dat$x + 0.5 * dat$x * dat$z + rnorm(100)
#' fit <- lm(y ~ x * z, data = dat)
#'
#' ## what JN() does internally
#' inp <- jn_input(fit, theta_1 = "x", theta_2 = "z")
#' inp
#'
#' ## the same thing built by hand, for a model class without a method
#' idx <- c("x", "z", "x:z")
#' jn_wald(coefficients = coef(fit)[idx],
#'         vcov         = vcov(fit)[idx, idx],
#'         labels       = c("x", "z"),
#'         ranges       = list(range(dat$x), range(dat$z)),
#'         support      = list(dat$x, dat$z))
#'
#' @export
jn_input <- function(object, ...) UseMethod("jn_input")


#' @rdname jn_input
#' @export
jn_wald <- function(coefficients, vcov, labels, ranges = NULL,
                    support = NULL) {
  labels <- as.character(labels)
  k      <- length(labels)
  .jn_check_k(k)

  coefficients <- as.numeric(coefficients)
  vcov         <- as.matrix(vcov)
  n_par        <- .jn_n_par(k)

  if (length(coefficients) != n_par)
    stop("`coefficients` must have ", n_par, " elements for a ", k,
         "-way interaction, not ", length(coefficients), ".", call. = FALSE)
  if (!identical(dim(vcov), c(n_par, n_par)))
    stop("`vcov` must be ", n_par, " by ", n_par, ", not ",
         paste(dim(vcov), collapse = " by "), ".", call. = FALSE)
  if (anyNA(coefficients) || anyNA(vcov))
    stop("`coefficients` and `vcov` must not contain missing values.",
         call. = FALSE)

  structure(
    list(coefficients = coefficients,
         vcov         = vcov,
         labels       = labels,
         ranges       = .jn_check_ranges(ranges, k),
         support      = .jn_check_support(support, k)),
    class = c("jn_wald", "jn_input"))
}


#' @rdname jn_input
#' @export
jn_posterior <- function(draws, labels, ranges = NULL, support = NULL) {
  labels <- as.character(labels)
  k      <- length(labels)
  .jn_check_k(k)

  draws <- as.matrix(draws)
  n_par <- .jn_n_par(k)

  if (ncol(draws) != n_par)
    stop("`draws` must have ", n_par, " columns for a ", k,
         "-way interaction, not ", ncol(draws), ".", call. = FALSE)
  if (nrow(draws) < 2L)
    stop("`draws` must have at least two rows.", call. = FALSE)
  if (anyNA(draws))
    stop("`draws` must not contain missing values.", call. = FALSE)

  structure(
    list(draws   = draws,
         labels  = labels,
         ranges  = .jn_check_ranges(ranges, k),
         support = .jn_check_support(support, k)),
    class = c("jn_posterior", "jn_input"))
}


.jn_check_k <- function(k) {
  if (!k %in% c(2L, 3L))
    stop("int3ract handles two- and three-way interactions, so `labels` must ",
         "have 2 or 3 elements, not ", k, ".", call. = FALSE)
  invisible(TRUE)
}

.jn_check_ranges <- function(ranges, k) {
  if (is.null(ranges)) return(vector("list", k))
  if (!is.list(ranges) || length(ranges) != k)
    stop("`ranges` must be a list of length ", k, ".", call. = FALSE)
  lapply(ranges, function(r) if (is.null(r)) NULL else as.numeric(r))
}

.jn_check_support <- function(support, k) {
  if (is.null(support)) return(NULL)
  if (!is.list(support) || length(support) != k)
    stop("`support` must be a list of length ", k,
         ", one element per focal variable.", call. = FALSE)
  support <- lapply(support, function(s) {
    if (is.null(s)) return(NULL)
    s <- suppressWarnings(as.numeric(s))
    s[is.finite(s)]
  })
  if (all(vapply(support, length, integer(1)) == 0L)) NULL else support
}


#' @export
print.jn_input <- function(x, ...) {
  type <- if (inherits(x, "jn_wald")) {
    "Wald (point estimates and covariance matrix)"
  } else {
    paste0("posterior (", nrow(x$draws), " draws)")
  }
  cat("<jn_input>\n")
  cat("  variables: ", paste(x$labels, collapse = ", "), "\n", sep = "")
  cat("  inference: ", type, "\n", sep = "")
  cat("  support:   ",
      if (is.null(x$support)) "not supplied" else "observed values available",
      "\n", sep = "")
  invisible(x)
}


#' @rdname jn_input
#' @export
jn_input.default <- function(object, ...) {
  stop("int3ract does not know how to extract Johnson-Neyman input from an ",
       "object of class ",
       paste(sQuote(class(object)), collapse = "/"), ".\n",
       "Supported out of the box: lm, glm, lmerMod/glmerMod, lavaan, ",
       "sienaFit, matrix/mcmc of draws, and sienaBayesFit/multiSiena.\n",
       "For anything else, either build the input by hand with jn_wald() or ",
       "jn_posterior(), or write a jn_input() method for the class. ",
       "See ?jn_input.", call. = FALSE)
}


# ---- lm / glm ---------------------------------------------------------------

#' @rdname jn_input
#' @export
jn_input.lm <- function(object, theta_1, theta_2, theta_3 = NULL, ...) {
  vars <- .jn_vars(theta_1, theta_2, theta_3)
  cf   <- stats::coef(object)
  idx  <- .jn_par_index(names(cf), vars)
  mf   <- object$model

  jn_wald(coefficients = cf[idx],
          vcov         = stats::vcov(object)[idx, idx, drop = FALSE],
          labels       = vars,
          ranges       = .jn_ranges_from_frame(mf, vars),
          support      = .jn_support_from_frame(mf, vars))
}


# ---- lme4 -------------------------------------------------------------------

#' @rdname jn_input
#' @export
jn_input.merMod <- function(object, theta_1, theta_2, theta_3 = NULL, ...) {
  if (!requireNamespace("lme4", quietly = TRUE))
    stop("Package 'lme4' is required for mixed-model input.", call. = FALSE)

  vars <- .jn_vars(theta_1, theta_2, theta_3)
  fe   <- lme4::fixef(object)
  idx  <- .jn_par_index(names(fe), vars)
  mf   <- stats::model.frame(object)

  jn_wald(coefficients = fe[idx],
          vcov         = as.matrix(stats::vcov(object))[idx, idx, drop = FALSE],
          labels       = vars,
          ranges       = .jn_ranges_from_frame(mf, vars),
          support      = .jn_support_from_frame(mf, vars))
}

#' @rdname jn_input
#' @export
jn_input.lmerModLmerTest <- jn_input.merMod


# ---- lavaan -----------------------------------------------------------------
#
# A structural model normally has several endogenous variables, so unlike a
# regression there is no single "the" interaction to analyse: `x:z` may predict
# two different outcomes with different coefficients. `outcome` names the
# dependent variable whose regression is wanted. It may be omitted when only
# one outcome is regressed on all of the focal variables, which is the case an
# lm user would recognize.

#' @rdname jn_input
#' @export
jn_input.lavaan <- function(object, theta_1, theta_2, theta_3 = NULL,
                            outcome = NULL,
                            theta_int_12 = NULL, theta_int_13 = NULL,
                            theta_int_23 = NULL, theta_int_123 = NULL, ...) {
  if (!requireNamespace("lavaan", quietly = TRUE))
    stop("Package 'lavaan' is required for structural equation model input.",
         call. = FALSE)
  if (inherits(object, "blavaan"))
    stop("A blavaan fit carries posterior draws, and analysing its posterior ",
         "means with Wald tests would throw that away.\nExtract the draws ",
         "(blavaan::blavInspect(object, \"draws\")) and pass the matrix to ",
         "JN() instead.", call. = FALSE)

  vars <- .jn_vars(theta_1, theta_2, theta_3)

  if (lavaan::lavInspect(object, "ngroups") > 1L)
    stop("Multiple-group lavaan fits are not supported: each group has its ",
         "own structural coefficients, so there is no single regression to ",
         "analyse.\nFit the groups separately and call JN() on each.",
         call. = FALSE)

  reg <- lavaan::parameterEstimates(object)
  reg <- reg[reg$op == "~", , drop = FALSE]
  if (!nrow(reg))
    stop("This model has no regression paths (`~`), so it contains no ",
         "interaction to analyse.", call. = FALSE)

  outcome <- .jn_lavaan_outcome(reg, vars, outcome)
  rhs     <- reg$rhs[reg$lhs == outcome]

  # Main effects first: a missing one is a plainer problem than a missing
  # product, and saying so beats a message about the interaction.
  absent <- setdiff(vars, rhs)
  if (length(absent))
    stop("Outcome ", sQuote(outcome), " is not regressed on ",
         paste(sQuote(absent), collapse = ", "),
         ".\nIts predictors are: ", paste(sQuote(rhs), collapse = ", "), ".",
         call. = FALSE)

  terms <- c(vars, .jn_lavaan_ints(rhs, vars, theta_int_12, theta_int_13,
                                   theta_int_23, theta_int_123))

  named <- setdiff(terms, rhs)
  if (length(named))
    stop("Outcome ", sQuote(outcome), " is not regressed on ",
         paste(sQuote(named), collapse = ", "),
         ".\nIts predictors are: ", paste(sQuote(rhs), collapse = ", "), ".",
         call. = FALSE)

  # lavaan is S4, and stats::coef()/stats::vcov() called by name bypass its
  # methods, so the estimates and their covariance are taken from lavaan's own
  # accessors. Both label structural paths as "<outcome>~<term>".
  idx <- paste0(outcome, "~", terms)
  vc  <- as.matrix(lavaan::lavInspect(object, "vcov"))

  fixed <- setdiff(idx, colnames(vc))
  if (length(fixed))
    stop("These paths are not freely estimated, so they carry no sampling ",
         "variance: ", paste(sQuote(fixed), collapse = ", "),
         ".\nA Johnson-Neyman analysis needs the covariance matrix of every ",
         "term involved.", call. = FALSE)

  rows <- match(terms, reg$rhs[reg$lhs == outcome])
  cf   <- reg$est[reg$lhs == outcome][rows]
  mf   <- .jn_lavaan_frame(object)

  jn_wald(coefficients = cf,
          vcov         = vc[idx, idx, drop = FALSE],
          labels       = vars,
          ranges       = .jn_ranges_from_frame(mf, vars),
          support      = .jn_support_from_frame(mf, vars))
}


# Which dependent variable's regression to analyse. One candidate is
# unambiguous; several is the normal SEM case and the caller has to choose.
.jn_lavaan_outcome <- function(reg, vars, outcome) {
  lhs <- unique(reg$lhs)

  if (!is.null(outcome)) {
    outcome <- as.character(outcome)
    if (length(outcome) != 1L)
      stop("`outcome` must name a single dependent variable, not ",
           length(outcome), ".", call. = FALSE)
    if (!outcome %in% lhs)
      stop("Nothing is regressed on ", sQuote(outcome),
           " in this model.\nIts outcomes are: ",
           paste(sQuote(lhs), collapse = ", "), ".", call. = FALSE)
    return(outcome)
  }

  cand <- lhs[vapply(lhs,
                     function(o) all(vars %in% reg$rhs[reg$lhs == o]),
                     logical(1))]

  if (length(cand) == 1L) return(cand)
  if (!length(cand))
    stop("No outcome in this model is regressed on all of ",
         paste(sQuote(vars), collapse = ", "), ".\nIts outcomes are: ",
         paste(sQuote(lhs), collapse = ", "), ".", call. = FALSE)
  stop(length(cand), " outcomes are regressed on ",
       paste(sQuote(vars), collapse = ", "), ": ",
       paste(sQuote(cand), collapse = ", "),
       ".\nUse `outcome` to say which of them to analyse.", call. = FALSE)
}


# lavaan writes two-way products as "x:z" in its own syntax, so those are found
# the usual way. It refuses three-way products outright -- they have to be
# built as a variable in the data, under a name of the user's choosing -- so
# for those the term name has to be supplied.
.jn_lavaan_ints <- function(rhs, vars, i12, i13, i23, i123) {
  pick <- function(given, arg, parts) {
    if (!is.null(given)) return(as.character(given))
    hit <- .find_int_quiet(rhs, parts)
    if (!is.null(hit)) return(hit)
    stop("No predictor matches the interaction of ",
         paste(sQuote(parts), collapse = " and "), ".\n",
         if (length(parts) == 3L)
           paste0("lavaan does not accept three-way products in its model ",
                  "syntax, so this term has to be built as a variable in the ",
                  "data. Pass its name as `", arg, "`.\n")
         else
           paste0("lavaan writes this product as \"",
                  paste(parts, collapse = ":"),
                  "\"; if it is named differently, pass that name as `",
                  arg, "`.\n"),
         "Predictors available: ", paste(sQuote(rhs), collapse = ", "), ".",
         call. = FALSE)
  }

  if (length(vars) == 2L)
    return(pick(i12, "theta_int_12", vars[1:2]))

  c(pick(i12,  "theta_int_12",  vars[c(1, 2)]),
    pick(i13,  "theta_int_13",  vars[c(1, 3)]),
    pick(i23,  "theta_int_23",  vars[c(2, 3)]),
    pick(i123, "theta_int_123", vars[c(1, 2, 3)]))
}


# The observed data, when the fit carries any. Models estimated from a sample
# covariance matrix carry none, and then the density panels are simply absent.
.jn_lavaan_frame <- function(object) {
  d <- try(lavaan::lavInspect(object, "data"), silent = TRUE)
  if (inherits(d, "try-error") || is.null(d) || !length(d)) return(NULL)
  as.data.frame(d, check.names = FALSE)
}


# ---- RSiena -----------------------------------------------------------------

#' @rdname jn_input
#' @export
jn_input.sienaFit <- function(object, theta_1, theta_2, theta_3 = NULL,
                              theta_int_12 = NULL, theta_int_13 = NULL,
                              theta_int_23 = NULL, theta_int_123 = NULL,
                              ranges = NULL, support = NULL, ...) {
  sn       <- object$effects$effectName
  threeWay <- !is.null(theta_3)

  idx <- if (threeWay) {
    c(theta_1, theta_2, theta_3,
      theta_int_12, theta_int_13, theta_int_23, theta_int_123)
  } else {
    c(theta_1, theta_2, theta_int_12)
  }
  nms <- if (threeWay) {
    c("theta_1", "theta_2", "theta_3", "theta_int_12",
      "theta_int_13", "theta_int_23", "theta_int_123")
  } else {
    c("theta_1", "theta_2", "theta_int_12")
  }

  if (length(idx) != length(nms))
    stop("For sienaFit input every interaction position must be supplied: ",
         paste(nms, collapse = ", "), ".", call. = FALSE)
  if (!is.numeric(idx))
    stop("For sienaFit input, effects are addressed by their integer position ",
         "in x$effects, not by name.", call. = FALSE)
  bad <- idx > length(sn) | idx < 1
  if (any(bad))
    stop("Parameter positions out of range: ",
         paste(nms[bad], collapse = ", "),
         "\nSee x$effects$effectName for valid positions.", call. = FALSE)

  labels <- if (threeWay) {
    sn[c(theta_1, theta_2, theta_3)]
  } else {
    sn[c(theta_1, theta_2)]
  }

  jn_wald(coefficients = object$theta[idx],
          vcov         = object$covtheta[idx, idx, drop = FALSE],
          labels       = labels,
          ranges       = ranges,
          support      = support)
}


# ---- posterior draws --------------------------------------------------------

#' @rdname jn_input
#' @export
jn_input.matrix <- function(object, theta_1, theta_2, theta_3 = NULL,
                            theta_int_12 = NULL, theta_int_13 = NULL,
                            theta_int_23 = NULL, theta_int_123 = NULL,
                            burn_in = 0, thin = 1,
                            ranges = NULL, support = NULL, ...) {
  threeWay <- !is.null(theta_3)
  by_name  <- is.character(theta_1)

  if (by_name != is.character(theta_2) ||
      (threeWay && by_name != is.character(theta_3)))
    stop("Address the parameters either all by name or all by column ",
         "position, not a mix of the two.", call. = FALSE)

  if (by_name) {
    cn <- colnames(object)
    if (is.null(cn))
      stop("The matrix of draws has no column names, so the parameters ",
           "cannot be addressed by name.", call. = FALSE)
    theta_int_12 <- theta_int_12 %||% .find_int(cn, theta_1, theta_2)
    if (threeWay) {
      theta_int_13  <- theta_int_13  %||% .find_int(cn, theta_1, theta_3)
      theta_int_23  <- theta_int_23  %||% .find_int(cn, theta_2, theta_3)
      theta_int_123 <- theta_int_123 %||% .find_int(cn, theta_1, theta_2,
                                                    theta_3)
    }
  }

  idx <- if (threeWay) {
    c(theta_1, theta_2, theta_3,
      theta_int_12, theta_int_13, theta_int_23, theta_int_123)
  } else {
    c(theta_1, theta_2, theta_int_12)
  }
  if (length(idx) != .jn_n_par(if (threeWay) 3L else 2L))
    stop("Every interaction parameter must be identified; supply the ",
         "theta_int_* arguments when the columns are addressed by position.",
         call. = FALSE)

  keep  <- seq(burn_in + 1, nrow(object), by = thin)
  draws <- object[keep, idx, drop = FALSE]

  labels <- if (by_name) {
    if (threeWay) c(theta_1, theta_2, theta_3) else c(theta_1, theta_2)
  } else {
    if (threeWay) c("V1", "V2", "V3") else c("V1", "V2")
  }

  jn_posterior(draws = draws, labels = labels,
               ranges = ranges, support = support)
}

#' @rdname jn_input
#' @export
jn_input.mcmc <- function(object, ...) {
  jn_input.matrix(as.matrix(object), ...)
}

#' @rdname jn_input
#' @export
jn_input.mcmc.list <- function(object, ...) {
  jn_input.matrix(do.call(rbind, lapply(object, as.matrix)), ...)
}

#' @rdname jn_input
#' @export
jn_input.data.frame <- function(object, ...) {
  jn_input.matrix(as.matrix(object), ...)
}


# ---- helpers ----------------------------------------------------------------

.jn_vars <- function(theta_1, theta_2, theta_3 = NULL) {
  vars <- c(theta_1, theta_2, theta_3)
  if (!is.character(vars))
    stop("For this model class the variables must be given by name.",
         call. = FALSE)
  vars
}

# Locate the main effects and every interaction among a coefficient vector.
.jn_par_index <- function(nms, vars) {
  if (length(vars) == 2L) {
    c(vars, .find_int(nms, vars[1], vars[2]))
  } else {
    c(vars,
      .find_int(nms, vars[1], vars[2]),
      .find_int(nms, vars[1], vars[3]),
      .find_int(nms, vars[2], vars[3]),
      .find_int(nms, vars[1], vars[2], vars[3]))
  }
}

.jn_ranges_from_frame <- function(mf, vars) {
  lapply(vars, function(v) {
    if (is.null(mf) || !v %in% names(mf)) return(NULL)
    r <- suppressWarnings(range(as.numeric(mf[[v]]), na.rm = TRUE))
    if (all(is.finite(r))) r else NULL
  })
}

.jn_support_from_frame <- function(mf, vars) {
  if (is.null(mf)) return(NULL)
  .jn_check_support(
    lapply(vars, function(v) {
      if (!v %in% names(mf)) return(NULL)
      suppressWarnings(as.numeric(mf[[v]]))
    }),
    length(vars))
}
