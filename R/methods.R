# =============================================================================
# Methods on the JN class.
# =============================================================================

.jn_inference_label <- function(x) {
  if (x$inference == "wald")
    paste0("Wald z tests, alpha = ", format(x$alpha))
  else
    paste0("conditional posteriors, thresholds = ",
           paste(format(x$thresholds), collapse = " and "))
}

#' Print a Johnson-Neyman analysis
#'
#' @param x a \code{JN} object.
#' @param ... ignored.
#' @returns \code{x}, invisibly.
#' @export
print.JN <- function(x, ...) {
  cat("Johnson-Neyman analysis (",
      if (x$n_way == 2L) "two-way" else "three-way", ", ",
      .jn_inference_label(x), ")\n", sep = "")
  cat("Variables: ", paste(x$labels, collapse = ", "),
      .jn_group_suffix(x), "\n", sep = "")
  cat("Grid:      ",
      paste(vapply(x$grid, length, integer(1)), collapse = " x "),
      " moderator values\n", sep = "")
  if (x$control_fdr)
    cat("           false discovery rate controlled across the grid\n")
  cat("\n")
  .jn_print_regions(jn_regions(x), x)
  cat("\nUse summary() for the full regions, plot() for the figures.\n")
  invisible(x)
}


#' Summarize a Johnson-Neyman analysis
#'
#' Reports the regions of significance: the stretches of the moderator over
#' which the conditional effect of each focal variable is distinguishable from
#' zero. Where the observed moderator values are known, each region is
#' annotated with the share of observations falling inside it, so that a region
#' resting on almost no data is visible as such.
#'
#' Regions that run to the edge of the evaluated range are marked with an
#' asterisk: their boundary is where the moderator stopped being evaluated, not
#' a Johnson-Neyman boundary.
#'
#' @param object a \code{JN} object.
#' @param at for three-way analyses, the values of the first moderator at which
#'   to report the region along the second. Defaults to the quartiles.
#' @param ... ignored.
#'
#' The range each variable was evaluated over is reported alongside the
#' regions, since a boundary can only ever be found inside it and the regions
#' are not interpretable without it.
#'
#' @returns An object of class \code{summary.JN}, holding the region table in
#'   its \code{regions} element and the evaluated range of each variable in its
#'   \code{ranges} element.
#'
#' @examples
#' set.seed(1402)
#' dat <- data.frame(x = rnorm(200), z = rnorm(200))
#' dat$y <- dat$x + 0.6 * dat$x * dat$z + rnorm(200)
#' summary(JN(lm(y ~ x * z, data = dat), theta_1 = "x", theta_2 = "z"))
#'
#' @export
summary.JN <- function(object, at = NULL, ...) {
  structure(list(regions = jn_regions(object, at = at),
                 ranges  = .jn_ranges_table(object),
                 object  = object),
            class = "summary.JN")
}

# The evaluated range of every variable, as a data frame so that it is
# available to callers and not only to the printed output.
.jn_ranges_table <- function(object) {
  data.frame(variable = object$labels,
             from     = vapply(object$grid, min,    numeric(1)),
             to       = vapply(object$grid, max,    numeric(1)),
             n        = vapply(object$grid, length, integer(1)),
             row.names        = NULL,
             stringsAsFactors = FALSE)
}


#' @rdname summary.JN
#' @param x a \code{summary.JN} object.
#' @export
print.summary.JN <- function(x, ...) {
  obj <- x$object
  cat("Johnson-Neyman analysis (",
      if (obj$n_way == 2L) "two-way" else "three-way", ", ",
      .jn_inference_label(obj), ")\n", sep = "")
  cat("Variables: ", paste(obj$labels, collapse = ", "),
      .jn_group_suffix(obj), "\n\n", sep = "")
  .jn_print_ranges(x$ranges)
  cat("\n")
  .jn_print_regions(x$regions, obj)
  invisible(x)
}


.jn_print_ranges <- function(rng) {
  cat("Evaluated over\n")
  w <- max(nchar(rng$variable))
  for (i in seq_len(nrow(rng)))
    cat("  ", formatC(rng$variable[i], width = -w),
        "   [", trimws(.jn_fmt(rng$from[i])),
        ", ",   trimws(.jn_fmt(rng$to[i])), "]",
        "   ", rng$n[i], " values\n", sep = "")
  invisible(rng)
}


# ---- printing the regions ---------------------------------------------------

.jn_fmt <- function(v) formatC(v, format = "f", digits = 3, width = 8)

.jn_print_regions <- function(reg, obj) {
  cat("Regions of significance\n")
  has_share <- any(!is.na(reg$data_share))
  clipped   <- FALSE

  for (f in unique(reg$focal)) {
    sub <- reg[reg$focal == f, , drop = FALSE]

    if (obj$n_way == 2L) {
      cat("\n  ", f, ", moderated by ", sub$moderator[1], "\n", sep = "")
      clipped <- .jn_print_block(sub, obj, has_share, indent = "    ",
                                 share_of = "observations") || clipped
    } else {
      cat("\n  ", f, ", moderated by ", sub$mod1[1], " and ",
          sub$moderator[1], "\n", sep = "")
      for (a in unique(sub$mod1_value)) {
        cat("    at ", sub$mod1[1], " = ", .jn_fmt(a), "\n", sep = "")
        clipped <- .jn_print_block(sub[sub$mod1_value == a, , drop = FALSE],
                                   obj, has_share, indent = "      ",
                                   share_of = paste(sub$moderator[1],
                                                    "values")) || clipped
      }
    }
  }

  if (clipped)
    cat("\n* reaches the edge of the evaluated range.\n", sep = "")

  invisible(reg)
}

# Prints one focal variable's regions. Returns TRUE if any of them was bounded
# by the evaluated range rather than by a genuine boundary.
.jn_print_block <- function(sub, obj, has_share, indent, share_of) {
  if (all(is.na(sub$from))) {
    cat(indent, "no region of significance over the evaluated range\n",
        sep = "")
    return(FALSE)
  }
  mod     <- sub$moderator[1]
  full    <- range(obj$grid[[match(mod, obj$labels)]])
  clipped <- FALSE

  for (r in seq_len(nrow(sub))) {
    if (is.na(sub$from[r])) next
    open_lo <- isTRUE(all.equal(sub$from[r], full[1]))
    open_hi <- isTRUE(all.equal(sub$to[r],   full[2]))
    clipped <- clipped || open_lo || open_hi

    where <- if (open_lo && open_hi) {
      paste0("the whole evaluated range of ", mod, " *")
    } else if (open_lo) {
      paste0(mod, " up to ", .jn_fmt(sub$to[r]), " *")
    } else if (open_hi) {
      paste0(mod, " from ", .jn_fmt(sub$from[r]), " upwards *")
    } else {
      paste0(mod, " from ", .jn_fmt(sub$from[r]),
             " to ", .jn_fmt(sub$to[r]))
    }

    cat(indent, where,
        "   effect ", if (identical(sub$sign[r], "+")) "positive"
                      else "negative", sep = "")
    if (has_share && !is.na(sub$data_share[r]))
      cat("   (", formatC(100 * sub$data_share[r], format = "f", digits = 1),
          "% of ", share_of, ")", sep = "")
    cat("\n")
  }
  clipped
}


# ---- extracting the numbers -------------------------------------------------

#' Extract the conditional effects of a Johnson-Neyman analysis
#'
#' \code{as.data.frame()} returns the full grid: one row per focal variable and
#' moderator value, with the conditional estimate, its uncertainty, and whether
#' it is significant. \code{coef()} and \code{vcov()} return the model
#' parameters the analysis was built from.
#'
#' @param x a \code{JN} object.
#' @param object a \code{JN} object.
#' @param row.names,optional ignored, for compatibility with the generic.
#' @param ... ignored.
#'
#' @returns A data frame, respectively a named numeric vector and a matrix.
#'
#' @examples
#' set.seed(1)
#' dat <- data.frame(x = rnorm(100), z = rnorm(100))
#' dat$y <- dat$x + 0.5 * dat$x * dat$z + rnorm(100)
#' jn <- JN(lm(y ~ x * z, data = dat), theta_1 = "x", theta_2 = "z",
#'          range_size = 5)
#' head(as.data.frame(jn))
#' coef(jn)
#'
#' @export
as.data.frame.JN <- function(x, row.names = NULL, optional = FALSE, ...) {
  out <- do.call(rbind, x$tables)
  rownames(out) <- NULL
  if (!is.null(x$group)) out <- cbind(group = x$group, out)
  out
}


#' @rdname as.data.frame.JN
#' @export
coef.JN <- function(object, ...) {
  nms <- .jn_par_names(object$labels)
  if (object$inference == "wald")
    stats::setNames(object$input$coefficients, nms)
  else
    stats::setNames(colMeans(object$input$draws), nms)
}

.jn_par_names <- function(labels) {
  if (length(labels) == 2L)
    c(labels, paste(labels[1], labels[2], sep = ":"))
  else
    c(labels,
      paste(labels[1], labels[2], sep = ":"),
      paste(labels[1], labels[3], sep = ":"),
      paste(labels[2], labels[3], sep = ":"),
      paste(labels, collapse = ":"))
}


#' @rdname as.data.frame.JN
#' @export
vcov.JN <- function(object, ...) {
  nms <- .jn_par_names(object$labels)
  V   <- if (object$inference == "wald") object$input$vcov
         else stats::var(object$input$draws)
  dimnames(V) <- list(nms, nms)
  V
}


# ---- grouped analyses -------------------------------------------------------

#' Methods for grouped Johnson-Neyman analyses
#'
#' A \code{JN_list} holds one analysis per group, as produced by
#' \code{\link{JN}} with \code{fixed_only = FALSE} for \pkg{lme4} models or
#' \code{hyper_only = FALSE} for \pkg{multiSiena} models. Its first element is
#' the fixed-effect, respectively hyper-parameter, analysis.
#'
#' @param x,object a \code{JN_list} object.
#' @param ... passed to the corresponding \code{JN} method.
#'
#' @returns \code{print()} returns its argument invisibly; \code{summary()} a
#'   list of summaries; \code{as.data.frame()} the grids of every group stacked
#'   into one data frame.
#'
#' @export
print.JN_list <- function(x, ...) {
  cat("Johnson-Neyman analyses for ", length(x), " groups: ",
      paste(utils::head(names(x), 6), collapse = ", "),
      if (length(x) > 6) ", ...", "\n\n", sep = "")
  print(x[[1]])
  if (length(x) > 1)
    cat("\n[", length(x) - 1, " further group-level analyses; ",
        "use x[[i]], summary(x) or plot(x)]\n", sep = "")
  invisible(x)
}

#' @rdname print.JN_list
#' @export
summary.JN_list <- function(object, ...) lapply(object, summary, ...)

#' @rdname print.JN_list
#' @export
as.data.frame.JN_list <- function(x, ...) {
  out <- do.call(rbind, lapply(names(x), function(g) {
    d <- as.data.frame(x[[g]])
    if (!"group" %in% names(d)) d <- cbind(group = g, d)
    d
  }))
  rownames(out) <- NULL
  out
}
