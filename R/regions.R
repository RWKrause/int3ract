# =============================================================================
# Regions of significance.
#
# The point of a Johnson-Neyman analysis is the *region*: the stretch of the
# moderator over which the focal effect is distinguishable from zero. The grid
# in the result tables approximates it; the functions here report it directly,
# and for the two-way Wald case solve for the boundaries exactly.
# =============================================================================

# Share of the observed moderator values falling inside [from, to]. This is the
# guard against the classic failure mode of JN analyses: a region of
# significance sitting where hardly any data were observed.
.jn_data_share <- function(obs, from, to) {
  if (is.null(obs) || !length(obs)) return(NA_real_)
  mean(obs >= from & obs <= to)
}


# Refine a run boundary by linear interpolation on a score that is positive
# where the effect is significant.
.jn_runs <- function(x, score) {
  sig <- score > 0
  sig[is.na(sig)] <- FALSE
  if (!any(sig)) return(data.frame(from = numeric(0), to = numeric(0)))

  r      <- rle(sig)
  ends   <- cumsum(r$lengths)
  starts <- c(1L, utils::head(ends, -1L) + 1L)
  keep   <- which(r$values)

  interp <- function(i, j) {
    # zero crossing of `score` between grid points i and j
    d <- score[j] - score[i]
    if (!is.finite(d) || d == 0) return(x[j])
    x[i] + (x[j] - x[i]) * (0 - score[i]) / d
  }

  do.call(rbind, lapply(keep, function(r_i) {
    s <- starts[r_i]; e <- ends[r_i]
    data.frame(from = if (s > 1L)          interp(s - 1L, s) else x[s],
               to   = if (e < length(x))   interp(e + 1L, e) else x[e])
  }))
}


# ---- exact two-way Wald boundaries ------------------------------------------

# theta(m) is significant where theta(m)^2 > z_crit^2 Var[theta(m)]. Both sides
# are quadratic in m, so the boundaries are the roots of
#   A m^2 + B m + C = 0,  A = b3^2 - z^2 V33,
#                         B = 2 (b1 b3 - z^2 V13),
#                         C = b1^2 - z^2 V11,
# and the significant set is where that quadratic is positive.
.jn_exact_2way <- function(b_main, b_int, V_main, V_int, V_cross,
                           alpha, lo, hi) {
  z2 <- stats::qnorm(alpha / 2)^2
  A  <- b_int^2 - z2 * V_int
  B  <- 2 * (b_main * b_int - z2 * V_cross)
  C  <- b_main^2 - z2 * V_main

  positive <- function(m) A * m^2 + B * m + C > 0

  if (isTRUE(all.equal(A, 0))) {
    if (isTRUE(all.equal(B, 0))) {
      iv <- if (C > 0) data.frame(from = -Inf, to = Inf)
            else       data.frame(from = numeric(0), to = numeric(0))
    } else {
      root <- -C / B
      iv <- if (B > 0) data.frame(from = root, to = Inf)
            else       data.frame(from = -Inf, to = root)
    }
  } else {
    disc <- B^2 - 4 * A * C
    if (disc <= 0) {
      iv <- if (positive(0)) data.frame(from = -Inf, to = Inf)
            else             data.frame(from = numeric(0), to = numeric(0))
    } else {
      r <- sort(c((-B - sqrt(disc)) / (2 * A), (-B + sqrt(disc)) / (2 * A)))
      iv <- if (A > 0) {
        data.frame(from = c(-Inf, r[2]), to = c(r[1], Inf))
      } else {
        data.frame(from = r[1], to = r[2])
      }
    }
  }

  if (!nrow(iv)) return(iv)
  # clip to the range actually evaluated, and drop what falls outside it
  iv$from <- pmax(iv$from, lo)
  iv$to   <- pmin(iv$to,   hi)
  iv[iv$from < iv$to, , drop = FALSE]
}


# ---- dispatcher -------------------------------------------------------------

#' Regions of significance of a Johnson-Neyman analysis
#'
#' Reports, for each focal variable, the stretch of the moderator over which
#' the conditional effect is significant. This is what
#' \code{\link{summary.JN}} prints.
#'
#' For two-way Wald analyses the boundaries are solved for exactly rather than
#' read off the grid, so their precision does not depend on \code{range_size}.
#' For posterior analyses, and for the three-way case, the boundaries are
#' obtained by linear interpolation between adjacent grid points.
#'
#' @param object a \code{JN} object.
#' @param at for three-way analyses, the values of the first moderator at which
#'   the region along the second moderator is reported. Defaults to its
#'   quartiles, which keeps the printed output readable; the full grid is
#'   always available through \code{\link[=as.data.frame.JN]{as.data.frame}}.
#' @param ... ignored.
#'
#' @returns A data frame with one row per region, holding the focal variable,
#'   the moderator, the interval, the sign of the effect inside it, and -- when
#'   the observed moderator values are known -- the share of observations that
#'   fall inside it.
#'
#' @seealso \code{\link{summary.JN}}
#' @export
jn_regions <- function(object, ...) UseMethod("jn_regions")


#' @rdname jn_regions
#' @export
jn_regions.JN_2way <- function(object, ...) {
  exact <- object$inference == "wald"

  out <- lapply(seq_along(object$labels), function(i) {
    tab <- object$tables[[i]]
    mod <- 3L - i
    obs <- object$support[[mod]]
    m   <- tab$mod_value

    iv <- if (exact) {
      V <- object$input$vcov; b <- object$input$coefficients
      .jn_exact_2way(b_main = b[i], b_int = b[3],
                     V_main = V[i, i], V_int = V[3, 3], V_cross = V[i, 3],
                     alpha = object$alpha, lo = min(m), hi = max(m))
    } else {
      # positive where the Bayesian p value is outside the threshold band
      score <- pmax(min(object$thresholds) - tab$bayes_p,
                    tab$bayes_p - max(object$thresholds))
      .jn_runs(m, score)
    }

    if (!nrow(iv)) {
      return(data.frame(focal = object$labels[i],
                        moderator = object$labels[mod],
                        from = NA_real_, to = NA_real_, sign = NA_character_,
                        data_share = NA_real_, stringsAsFactors = FALSE))
    }

    mid <- (iv$from + iv$to) / 2
    est <- if (exact) {
      object$input$coefficients[i] + object$input$coefficients[3] * mid
    } else {
      stats::approx(m, tab$post_mean, xout = mid, rule = 2)$y
    }

    data.frame(focal      = object$labels[i],
               moderator  = object$labels[mod],
               from       = iv$from,
               to         = iv$to,
               sign       = ifelse(est > 0, "+", "-"),
               data_share = mapply(.jn_data_share,
                                   list(obs), iv$from, iv$to),
               stringsAsFactors = FALSE)
  })

  res <- do.call(rbind, out)
  rownames(res) <- NULL
  res
}


#' @rdname jn_regions
#' @export
jn_regions.JN_3way <- function(object, at = NULL, ...) {
  out <- lapply(seq_along(object$labels), function(i) {
    tab  <- object$tables[[i]]
    m1   <- sort(unique(tab$mod1_value))
    mod2 <- tab$mod2[1]
    obs2 <- object$support[[match(mod2, object$labels)]]

    at_i <- at %||% stats::quantile(m1, c(0.25, 0.5, 0.75), names = FALSE)
    at_i <- vapply(at_i, function(a) m1[which.min(abs(m1 - a))], numeric(1))
    at_i <- unique(at_i)

    do.call(rbind, lapply(at_i, function(a) {
      slice <- tab[tab$mod1_value == a, , drop = FALSE]
      slice <- slice[order(slice$mod2_value), , drop = FALSE]

      score <- if (object$inference == "wald") {
        abs(slice$statistic) - abs(stats::qnorm(object$alpha / 2))
      } else {
        pmax(min(object$thresholds) - slice$bayes_p,
             slice$bayes_p - max(object$thresholds))
      }
      iv <- .jn_runs(slice$mod2_value, score)

      if (!nrow(iv))
        return(data.frame(focal = object$labels[i],
                          mod1 = slice$mod1[1], mod1_value = a,
                          moderator = mod2,
                          from = NA_real_, to = NA_real_,
                          sign = NA_character_, data_share = NA_real_,
                          stringsAsFactors = FALSE))

      mid <- (iv$from + iv$to) / 2
      est <- stats::approx(slice$mod2_value,
                           if (object$inference == "wald") slice$estimate
                           else slice$post_mean,
                           xout = mid, rule = 2)$y

      data.frame(focal      = object$labels[i],
                 mod1       = slice$mod1[1],
                 mod1_value = a,
                 moderator  = mod2,
                 from       = iv$from,
                 to         = iv$to,
                 sign       = ifelse(est > 0, "+", "-"),
                 data_share = mapply(.jn_data_share,
                                     list(obs2), iv$from, iv$to),
                 stringsAsFactors = FALSE)
    }))
  })

  res <- do.call(rbind, out)
  rownames(res) <- NULL
  res
}
