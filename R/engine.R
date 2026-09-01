# =============================================================================
# Computation. Nothing in this file draws anything: every function here turns a
# jn_input carrier into tidy data frames. Plotting lives in R/plot.R and is
# only ever triggered by plot()/autoplot() on the finished object.
# =============================================================================

# Conditional effect of a focal variable given two moderator values, and its
# delta-method standard error (Equations 4 and 5 of the manuscript).
.jn_theta3 <- function(m1, m2, bx, bm1x, bm2x, bm1m2x) {
  bx + m1 * bm1x + m2 * bm2x + m1 * m2 * bm1m2x
}

.jn_se3 <- function(m1, m2,
                    Vx, Vm1x, Vm2x, Vm1m2x,
                    covX_m1x, covX_m2x, covX_m1m2x,
                    covm1x_m2x, covm1x_m1m2x, covm2x_m1m2x) {
  sqrt(Vx +
         m1^2 * Vm1x +
         m2^2 * Vm2x +
         (m1 * m2)^2 * Vm1m2x +
         2 * m1 * covX_m1x +
         2 * m2 * covX_m2x +
         2 * m1 * m2 * covX_m1m2x +
         2 * m1 * m2 * covm1x_m2x +
         2 * m1^2 * m2 * covm1x_m1m2x +
         2 * m1 * m2^2 * covm2x_m1m2x)
}


# ---- moderator grids --------------------------------------------------------

# A supplied vector of length two is read as a range and expanded to
# `range_size` points; any other length is used verbatim, so that users can ask
# for exactly the moderator values they care about.
.jn_expand <- function(vals, range_size) {
  vals <- as.numeric(vals)
  if (length(vals) == 2L) seq(vals[1], vals[2], length.out = range_size)
  else sort(unique(vals))
}

.jn_grid <- function(input, ranges, range_size) {
  k     <- length(input$labels)
  given <- input$ranges

  lapply(seq_len(k), function(i) {
    v <- ranges[[i]] %||% given[[i]]
    if (is.null(v))
      stop("No values to evaluate ", sQuote(input$labels[i]), " over. ",
           "This model class does not carry its data, so supply the ",
           "moderator ranges explicitly, e.g. theta_", i,
           "_vals = c(-3, 3).", call. = FALSE)
    .jn_expand(v, range_size)
  })
}


# ---- significance flags -----------------------------------------------------

# The `control_fdr` argument applies the Benjamini-Hochberg step-up procedure
# across the grid, which is what is wanted here: the tests along a moderator
# range are many and highly dependent, and controlling the false discovery rate
# is far less brutal than controlling the family-wise error rate.
.jn_significant <- function(p, alpha, control_fdr) {
  if (!control_fdr) return(p < alpha)
  ord     <- order(p)
  sig_ord <- p[ord] < alpha * seq_along(p) / length(p)
  if (!all(sig_ord))
    sig_ord[which.max(!sig_ord):length(p)] <- FALSE
  out      <- logical(length(p))
  out[ord] <- sig_ord
  out
}

.jn_two_sided_p <- function(z) 2 * pmin(stats::pnorm(z), 1 - stats::pnorm(z))


# ---- two-way, Wald ----------------------------------------------------------

.jn_wald_2way <- function(input, grid, alpha, control_fdr) {
  b      <- input$coefficients
  V      <- input$vcov
  labels <- input$labels
  z_crit <- abs(stats::qnorm(alpha / 2))

  tables <- lapply(1:2, function(i) {
    mod <- 3L - i                       # the other main effect
    m   <- grid[[mod]]

    est <- b[i] + b[3] * m
    # Var[theta_i(m)] = Var(b_i) + 2m Cov(b_i, b_int) + m^2 Var(b_int).
    # Note the covariance is between the focal main effect and the
    # *interaction*, not between the two main effects.
    se  <- sqrt(V[i, i] + 2 * m * V[i, 3] + m^2 * V[3, 3])
    z   <- est / se
    p   <- .jn_two_sided_p(z)

    data.frame(focal       = labels[i],
               moderator   = labels[mod],
               mod_value   = m,
               estimate    = est,
               std.error   = se,
               statistic   = z,
               p.value     = p,
               conf.low    = est - z_crit * se,
               conf.high   = est + z_crit * se,
               significant = .jn_significant(p, alpha, control_fdr),
               stringsAsFactors = FALSE)
  })
  stats::setNames(tables, labels)
}


# ---- three-way, Wald --------------------------------------------------------

.jn_wald_3way <- function(input, grid, alpha, control_fdr) {
  b      <- input$coefficients
  V      <- input$vcov
  labels <- input$labels

  x_idx  <- c(1L, 2L, 3L)
  m1_idx <- c(4L, 4L, 5L)
  m2_idx <- c(5L, 6L, 6L)
  v_idx  <- list(c(2L, 3L), c(1L, 3L), c(1L, 2L))

  tables <- lapply(1:3, function(i) {
    xi <- x_idx[i]; m1i <- m1_idx[i]; m2i <- m2_idx[i]
    vi <- v_idx[[i]]
    g  <- expand.grid(mod1_value = grid[[vi[1]]],
                      mod2_value = grid[[vi[2]]])

    est <- .jn_theta3(g$mod1_value, g$mod2_value,
                      bx = b[xi], bm1x = b[m1i],
                      bm2x = b[m2i], bm1m2x = b[7])
    se  <- .jn_se3(g$mod1_value, g$mod2_value,
                   Vx           = V[xi,  xi],
                   Vm1x         = V[m1i, m1i],
                   Vm2x         = V[m2i, m2i],
                   Vm1m2x       = V[7,   7],
                   covX_m1x     = V[xi,  m1i],
                   covX_m2x     = V[xi,  m2i],
                   covX_m1m2x   = V[xi,  7],
                   covm1x_m2x   = V[m1i, m2i],
                   covm1x_m1m2x = V[m1i, 7],
                   covm2x_m1m2x = V[m2i, 7])
    z <- est / se
    p <- .jn_two_sided_p(z)

    data.frame(focal       = labels[i],
               mod1        = labels[vi[1]],
               mod2        = labels[vi[2]],
               mod1_value  = g$mod1_value,
               mod2_value  = g$mod2_value,
               estimate    = est,
               std.error   = se,
               statistic   = z,
               p.value     = p,
               significant = .jn_significant(p, alpha, control_fdr),
               stringsAsFactors = FALSE)
  })
  stats::setNames(tables, labels)
}


# ---- two-way, posterior -----------------------------------------------------

# Returns both the summary table and the conditional draws themselves; the
# two-way Bayesian figure is a set of overlaid densities, so plot() needs the
# draws, not just their summaries.
.jn_post_2way <- function(input, grid, thresholds) {
  draws  <- input$draws
  labels <- input$labels
  n_it   <- nrow(draws)

  out <- lapply(1:2, function(i) {
    mod <- 3L - i
    m   <- grid[[mod]]

    # n_it x length(m): draw s at moderator value m equals b_i^(s) + b_int^(s) m
    cond <- outer(draws[, 3], m, "*") + draws[, i]

    tab <- data.frame(
      focal       = labels[i],
      moderator   = labels[mod],
      mod_value   = m,
      post_mean   = colMeans(cond),
      post_sd     = apply(cond, 2, stats::sd),
      bayes_p     = colSums(cond > 0) / n_it,
      q2.5        = apply(cond, 2, stats::quantile, probs = 0.025),
      q97.5       = apply(cond, 2, stats::quantile, probs = 0.975),
      stringsAsFactors = FALSE)
    tab$significant <- tab$bayes_p <= min(thresholds) |
                       tab$bayes_p >= max(thresholds)
    rownames(tab) <- NULL

    list(table = tab, draws = cond)
  })

  list(tables = stats::setNames(lapply(out, `[[`, "table"), labels),
       draws  = stats::setNames(lapply(out, `[[`, "draws"), labels))
}


# ---- three-way, posterior ---------------------------------------------------

.jn_post_3way <- function(input, grid, thresholds) {
  draws  <- input$draws
  labels <- input$labels
  n_it   <- nrow(draws)

  # columns of `draws` entering theta_i, in the order (main, m1, m2, m1:m2)
  cols  <- list(c(1L, 4L, 5L, 7L), c(2L, 4L, 6L, 7L), c(3L, 5L, 6L, 7L))
  v_idx <- list(c(2L, 3L), c(1L, 3L), c(1L, 2L))

  tables <- lapply(1:3, function(i) {
    vi <- v_idx[[i]]
    g  <- expand.grid(mod1_value = grid[[vi[1]]],
                      mod2_value = grid[[vi[2]]])

    design <- cbind(1, g$mod1_value, g$mod2_value,
                    g$mod1_value * g$mod2_value)
    cond   <- draws[, cols[[i]], drop = FALSE] %*% t(design)

    tab <- data.frame(
      focal      = labels[i],
      mod1       = labels[vi[1]],
      mod2       = labels[vi[2]],
      mod1_value = g$mod1_value,
      mod2_value = g$mod2_value,
      post_mean  = colMeans(cond),
      post_sd    = apply(cond, 2, stats::sd),
      bayes_p    = colSums(cond > 0) / n_it,
      q2.5       = apply(cond, 2, stats::quantile, probs = 0.025),
      q97.5      = apply(cond, 2, stats::quantile, probs = 0.975),
      stringsAsFactors = FALSE)
    tab$significant <- tab$bayes_p <= min(thresholds) |
                       tab$bayes_p >= max(thresholds)
    rownames(tab) <- NULL
    tab
  })
  stats::setNames(tables, labels)
}


# ---- assembly ---------------------------------------------------------------

# The single place where a jn_input becomes a JN object.
.jn_build <- function(input, ranges, range_size, alpha, control_fdr,
                      thresholds, group = NULL) {
  k        <- length(input$labels)
  is_wald  <- inherits(input, "jn_wald")

  if (is.null(range_size))
    range_size <- if (is_wald) {
      if (k == 2L) 1000L else 50L
    } else {
      if (k == 2L) 13L else 50L
    }

  grid <- .jn_grid(input, ranges, range_size)

  cond <- NULL
  if (is_wald) {
    tables <- if (k == 2L) .jn_wald_2way(input, grid, alpha, control_fdr)
              else         .jn_wald_3way(input, grid, alpha, control_fdr)
  } else {
    if (k == 2L) {
      res    <- .jn_post_2way(input, grid, thresholds)
      tables <- res$tables
      cond   <- res$draws
    } else {
      tables <- .jn_post_3way(input, grid, thresholds)
    }
  }

  structure(
    list(labels     = input$labels,
         n_way      = k,
         inference  = if (is_wald) "wald" else "posterior",
         alpha      = alpha,
         thresholds = thresholds,
         control_fdr = control_fdr,
         grid       = stats::setNames(grid, input$labels),
         support    = input$support,
         tables     = tables,
         cond_draws = cond,
         group      = group,
         input      = input),
    class = c(if (is_wald) "JN_wald" else "JN_posterior",
              if (k == 2L) "JN_2way" else "JN_3way",
              "JN"))
}
