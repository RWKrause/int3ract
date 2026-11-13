# =============================================================================
# Simple-slopes figures: the regions of significance drawn in outcome space.
#
# The Johnson-Neyman figure puts the moderator on the x axis and the focal
# effect on the y axis. This one puts the focal variable on the x axis and the
# outcome on the y axis, and draws the simple slopes at the moderator values
# bounding each region of significance, shading the fan between them. Every
# slope inside that fan is significant with the sign of its colour.
#
# The regions are recomputed here over the moderator's *observed* range: the
# fan should cover the people in the data, not the margin JN() adds around them.
# =============================================================================

# ---- algebra ----------------------------------------------------------------

# Point estimates of the canonical coefficients: the estimates themselves for
# Wald input, the posterior means for draws.
.jn_beta <- function(x) {
  if (x$inference == "wald") x$input$coefficients
  else colMeans(x$input$draws)
}

# The linear predictor at the rows of `v` (one column per variable), up to the
# constant when the model did not supply one.
.jn_linpred <- function(x, v) {
  (x$input$intercept %||% 0) + .jn_terms_sum(.jn_beta(x), v)
}

# The focal effect along moderator r, with the remaining variable of a
# three-way interaction held at vals, is linear in the moderator value m:
# theta(m) = L1'b + m L2'b. Both contrasts are returned.
.jn_slice_contrasts <- function(x, i, r, vals) {
  v0 <- vals; v0[r] <- 0
  v1 <- vals; v1[r] <- 1
  L1 <- .jn_contrast(x$n_way, i, v0)
  list(L1 = L1, L2 = .jn_contrast(x$n_way, i, v1) - L1)
}

# Regions of significance of that slice over [lo, hi], with their sign. For
# Wald input the slice is a two-way JN problem in its own right, so the
# boundaries are solved for exactly, in the three-way case as well. For draws
# they are found on a fine grid of the moderator.
.jn_slice_regions <- function(x, i, r, vals, lo, hi) {
  L <- .jn_slice_contrasts(x, i, r, vals)
  b <- .jn_beta(x)

  if (x$inference == "wald") {
    V  <- x$input$vcov
    iv <- .jn_exact_2way(b_main  = sum(L$L1 * b),
                         b_int   = sum(L$L2 * b),
                         V_main  = drop(L$L1 %*% V %*% L$L1),
                         V_int   = drop(L$L2 %*% V %*% L$L2),
                         V_cross = drop(L$L1 %*% V %*% L$L2),
                         alpha = x$alpha, lo = lo, hi = hi)
  } else {
    D  <- x$input$draws
    a  <- drop(D %*% L$L1)
    s  <- drop(D %*% L$L2)
    m  <- seq(lo, hi, length.out = 201)
    bp <- vapply(m, function(mm) mean(a + s * mm > 0), numeric(1))
    th <- x$thresholds
    iv <- .jn_runs(m, pmax(min(th) - bp, bp - max(th)))
  }

  iv <- iv[iv$to - iv$from > 1e-8 * max(1, abs(hi - lo)), , drop = FALSE]
  mid <- (iv$from + iv$to) / 2
  iv$sign <- if (nrow(iv))
    ifelse(sum(L$L1 * b) + mid * sum(L$L2 * b) > 0, "+", "-")
  else character(0)
  iv
}

# Estimate and test of the focal effect at one set of values of the others.
.jn_cond_at <- function(x, i, vals) {
  L <- .jn_contrast(x$n_way, i, vals)
  if (x$inference == "wald") {
    est <- sum(L * x$input$coefficients)
    se  <- sqrt(drop(L %*% x$input$vcov %*% L))
    p   <- .jn_two_sided_p(est / se)
    list(estimate = est, p = p, significant = p < x$alpha)
  } else {
    d  <- drop(x$input$draws %*% L)
    bp <- mean(d > 0)
    list(estimate = mean(d), p = bp,
         significant = bp <= min(x$thresholds) || bp >= max(x$thresholds))
  }
}


# ---- which slopes to draw ---------------------------------------------------

# The simple slopes for one panel: focal variable i, moderator r, the other
# variable of a three-way interaction held at vals. Returns the lines (one row
# per moderator value, with its category) and the bands to shade between them.
.jn_slopes_choice <- function(x, i, r, vals) {
  rng <- x$observed_range[[r]]
  obs <- x$support[[r]]

  line <- function(m, category, band = NA_integer_)
    data.frame(mod_value = m, category = category, band = band,
               stringsAsFactors = FALSE)
  no_bands <- data.frame(band = integer(0), from = numeric(0),
                         to = numeric(0), sign = character(0))

  # A binary moderator has two slopes and nothing in between them.
  if (isTRUE(x$binary[[r]])) {
    lines <- do.call(rbind, lapply(x$grid[[r]], function(m) {
      vv <- vals; vv[r] <- m
      ca <- .jn_cond_at(x, i, vv)
      line(m, if (!ca$significant) "ns"
              else if (ca$estimate > 0) "pos" else "neg")
    }))
    return(list(lines = lines, bands = no_bands))
  }

  iv <- .jn_slice_regions(x, i, r, vals, rng[1], rng[2])

  # Nothing significant anywhere in the data: one slope, at the mean.
  if (!nrow(iv)) {
    m0 <- if (!is.null(obs)) mean(obs) else mean(rng)
    return(list(lines = line(m0, "ns"), bands = no_bands))
  }

  bands <- data.frame(band = seq_len(nrow(iv)), from = iv$from, to = iv$to,
                      sign = iv$sign, stringsAsFactors = FALSE)
  cat   <- ifelse(iv$sign == "+", "pos", "neg")
  lines <- rbind(line(iv$from, cat, bands$band),
                 line(iv$to,   cat, bands$band))

  # Significant throughout: the two slopes are the ends of the data, and a
  # third at the densest moderator value shows where most people are.
  tol   <- 1e-8 * max(1, diff(rng))
  whole <- nrow(iv) == 1L && iv$from <= rng[1] + tol && iv$to >= rng[2] - tol
  if (whole && !is.null(obs) && length(unique(obs)) > 2L) {
    dn    <- stats::density(obs, from = rng[1], to = rng[2])
    lines <- rbind(lines, line(dn$x[which.max(dn$y)], "mode"))
  }
  list(lines = lines, bands = bands)
}


# ---- evaluating them over the focal variable --------------------------------

# Values of the focal variable to draw over: its two values if binary, its
# observed range otherwise.
.jn_focal_values <- function(x, i) {
  if (isTRUE(x$binary[[i]])) return(x$grid[[i]])
  rng <- x$observed_range[[i]]
  seq(rng[1], rng[2], length.out = 101)
}

.jn_slopes_frames <- function(x, i, r, vals, style, panel = "") {
  ch <- .jn_slopes_choice(x, i, r, vals)
  xf <- .jn_focal_values(x, i)

  y_at <- function(m) {
    v <- matrix(rep(vals, each = length(xf)), ncol = x$n_way)
    v[, i] <- xf
    v[, r] <- m
    .jn_linpred(x, v)
  }

  lines <- do.call(rbind, lapply(seq_len(nrow(ch$lines)), function(k) {
    data.frame(panel = panel, line = paste(panel, k),
               focal = xf, y = y_at(ch$lines$mod_value[k]),
               category = ch$lines$category[k],
               mod_value = ch$lines$mod_value[k],
               stringsAsFactors = FALSE)
  }))

  # The slopes at the ends of the observed moderator range fix the vertical
  # scale, drawn or not: a lone near-flat slope must not be zoomed into
  # looking steep.
  rng  <- x$observed_range[[r]]
  span <- data.frame(panel = panel, focal = rep(xf, 2),
                     y = c(y_at(rng[1]), y_at(rng[2])),
                     stringsAsFactors = FALSE)

  bands <- NULL
  if (nrow(ch$bands)) {
    obs  <- x$support[[r]]
    dens <- style$slope_shading == "density" && !is.null(obs) &&
            length(unique(obs)) > 2L
    if (dens) dn <- stats::density(obs, from = rng[1], to = rng[2])
    bands <- do.call(rbind, lapply(seq_len(nrow(ch$bands)), function(k) {
      bd <- ch$bands[k, ]
      # flat: one wedge between the bounding slopes. density: many thin
      # wedges, each as opaque as the data are dense at its moderator value.
      cuts  <- if (dens) seq(bd$from, bd$to, length.out = 61)
               else c(bd$from, bd$to)
      mids  <- (utils::head(cuts, -1L) + cuts[-1L]) / 2
      alpha <- if (dens) {
        dd <- stats::approx(dn$x, dn$y, xout = mids, rule = 2)$y
        pmin(0.85, 2 * style$shade_alpha * dd / max(dn$y))
      } else {
        style$shade_alpha
      }
      do.call(rbind, lapply(seq_along(mids), function(j) {
        y1 <- y_at(cuts[j]); y2 <- y_at(cuts[j + 1L])
        data.frame(panel = panel, group = paste(panel, k, j),
                   focal = xf, ymin = pmin(y1, y2), ymax = pmax(y1, y2),
                   sign = bd$sign, alpha = alpha[j],
                   stringsAsFactors = FALSE)
      }))
    }))
  }
  list(lines = lines, bands = bands, span = span)
}


# ---- drawing ----------------------------------------------------------------

.jn_slope_categories <- function(style) {
  list(values   = c(pos = .jn_pos_color(style), neg = .jn_neg_color(style),
                    ns = "grey45", mode = style$mode_color),
       linetype = c(pos = "solid", neg = "solid", ns = "dashed",
                    mode = "dotted"),
       labels   = c(pos = "significantly positive",
                    neg = "significantly negative",
                    ns = "not significant", mode = "most observations"))
}

.jn_outcome_axis <- function(x) {
  out <- x$input$outcome
  if (is.null(x$input$intercept)) {
    paste0(out %||% "Linear predictor", " (up to a constant)")
  } else {
    out %||% "Outcome"
  }
}

.jn_fmt_value <- function(v)
  vapply(v, function(u) format(signif(u, 3)), character(1))

.jn_draw_slopes <- function(fr, x, i, mod, style, xlim) {
  cats <- .jn_slope_categories(style)
  used <- intersect(names(cats$values), unique(fr$lines$category))

  p <- ggplot2::ggplot() + ggplot2::theme_bw() +
    ggplot2::geom_blank(data = fr$span,
                        ggplot2::aes(.data$focal, .data$y))

  if (!is.null(fr$bands))
    p <- p +
      ggplot2::geom_ribbon(
        data = fr$bands,
        ggplot2::aes(x = .data$focal, ymin = .data$ymin, ymax = .data$ymax,
                     group = .data$group, fill = .data$sign,
                     alpha = .data$alpha),
        colour = NA) +
      ggplot2::scale_fill_manual(values = c(`+` = .jn_pos_color(style),
                                            `-` = .jn_neg_color(style)),
                                 guide = "none") +
      ggplot2::scale_alpha_identity()

  # each slope is labelled with its moderator value at its right-hand end,
  # pushed apart where slopes meet there
  ends <- do.call(rbind, lapply(split(fr$lines, fr$lines$line),
                                function(d) d[nrow(d), , drop = FALSE]))
  ends$label <- paste(mod, "=", .jn_fmt_value(ends$mod_value))
  ends <- do.call(rbind, lapply(split(ends, ends$panel), function(e) {
    yr <- range(c(fr$span$y, fr$lines$y))
    e$y <- .jn_spread(e$y, 0.06 * diff(yr))
    e
  }))

  p <- p +
    ggplot2::geom_line(data = fr$lines,
                       ggplot2::aes(.data$focal, .data$y, group = .data$line,
                                    colour = .data$category,
                                    linetype = .data$category),
                       linewidth = 0.9) +
    ggplot2::geom_text(data = ends,
                       ggplot2::aes(.data$focal, .data$y, label = .data$label,
                                    colour = .data$category),
                       hjust = -0.08, size = 2.8, show.legend = FALSE) +
    ggplot2::scale_colour_manual(values = cats$values[used],
                                 labels = cats$labels[used],
                                 breaks = used, name = NULL) +
    ggplot2::scale_linetype_manual(values = cats$linetype[used],
                                   labels = cats$labels[used],
                                   breaks = used, name = NULL) +
    ggplot2::coord_cartesian(xlim = xlim, clip = "off") +
    ggplot2::labs(x = x$labels[i], y = .jn_outcome_axis(x)) +
    ggplot2::theme(legend.position = "bottom",
                   plot.subtitle = ggplot2::element_text(size = 9))

  if (isTRUE(x$binary[[i]]))
    p <- p + ggplot2::scale_x_continuous(breaks = x$grid[[i]])
  p
}

# Move label positions apart until neighbours are at least `gap` apart, keeping
# their centre where it was.
.jn_spread <- function(y, gap) {
  if (length(y) < 2L || !is.finite(gap) || gap <= 0) return(y)
  o <- order(y)
  s <- y[o]
  for (k in seq_along(s)[-1L]) s[k] <- max(s[k], s[k - 1L] + gap)
  s <- s - (mean(s) - mean(y))
  y[o] <- s
  y
}

# x limits leaving room on the right for the slope labels
.jn_slopes_xlim <- function(x, i) {
  xf  <- .jn_focal_values(x, i)
  rng <- range(xf)
  w   <- diff(rng)
  c(rng[1] - 0.03 * w, rng[2] + 0.3 * w)
}

.jn_slopes_subtitle <- function(lines, focal, mod) {
  cats <- unique(lines$category)
  txt  <- if (identical(cats, "ns"))
    paste0("Not significant for any observed ", mod,
           "; slope drawn at its mean")
  else
    paste0("Shaded where the effect of ", focal, " is significant (observed ",
           mod, " only)")
  .jn_wrap(txt)
}

# Subtitles are wrapped so that they fit a figure a column wide.
.jn_wrap <- function(txt, width = 48) paste(strwrap(txt, width), collapse = "\n")


# ---- the figures ------------------------------------------------------------

.jn_plot_2way_slopes <- function(x, i, style) {
  r    <- 3L - i
  mod  <- x$labels[r]
  fr   <- .jn_slopes_frames(x, i, r, c(0, 0), style)
  xlim <- .jn_slopes_xlim(x, i)

  p <- .jn_draw_slopes(fr, x, i, mod, style, xlim) +
    ggplot2::labs(
      title = paste0("Simple slopes of ", x$labels[i], " by ", mod,
                     .jn_group_suffix(x)),
      subtitle = if (isTRUE(x$binary[[r]]))
        paste0("Slopes at both values of ", mod)
      else .jn_slopes_subtitle(fr$lines, x$labels[i], mod))

  # The panel below shows the focal variable this time: it is the x axis.
  obs <- .jn_support_for(x, x$labels[i])
  if (!style$show_density || is.null(obs)) return(p)
  binary <- isTRUE(x$binary[[i]])
  dens <- .jn_density_bottom(obs, xlim, x$labels[i], style, binary = binary)
  p <- p + ggplot2::labs(x = NULL) +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank())
  # the legend goes underneath everything, not between the two panels
  patchwork::wrap_plots(p, dens, ncol = 1,
                        heights = c(1, style$density_height),
                        guides = "collect") &
    ggplot2::theme(legend.position = "bottom")
}


# Three-way regions are two-dimensional, so the figure shows slices: one panel
# per value of one moderator, each the two-way figure along the other. A binary
# moderator is always the one sliced on; otherwise the first moderator is, at
# `at` or by default at its quartiles.
.jn_plot_3way_slopes <- function(x, i, style, at = NULL) {
  tab  <- x$tables[[i]]
  role <- .jn_3way_roles(x, i)
  cnm  <- tab[[role$cond]][1]
  rnm  <- tab[[role$run]][1]
  ci   <- match(cnm, x$labels)
  r    <- match(rnm, x$labels)

  if (isTRUE(x$binary[[ci]])) {
    at   <- x$grid[[ci]]
    tags <- rep("", length(at))
  } else if (is.null(at)) {
    src  <- x$support[[ci]] %||% x$observed_range[[ci]]
    at   <- stats::quantile(src, c(0.25, 0.5, 0.75), names = FALSE)
    tags <- c(" (Q1)", " (median)", " (Q3)")
  } else {
    tags <- rep("", length(at))
  }
  panels <- paste0(cnm, " = ", .jn_fmt_value(at), tags)

  frs <- lapply(seq_along(at), function(k) {
    vals <- numeric(3); vals[ci] <- at[k]
    .jn_slopes_frames(x, i, r, vals, style, panel = panels[k])
  })
  fr <- list(lines = do.call(rbind, lapply(frs, `[[`, "lines")),
             bands = do.call(rbind, lapply(frs, `[[`, "bands")),
             span  = do.call(rbind, lapply(frs, `[[`, "span")))
  fr$lines$panel <- factor(fr$lines$panel, levels = panels)
  fr$span$panel  <- factor(fr$span$panel,  levels = panels)
  if (!is.null(fr$bands))
    fr$bands$panel <- factor(fr$bands$panel, levels = panels)

  .jn_draw_slopes(fr, x, i, rnm, style, .jn_slopes_xlim(x, i)) +
    ggplot2::facet_wrap(ggplot2::vars(.data$panel), nrow = 1) +
    ggplot2::labs(
      title = paste0("Simple slopes of ", x$labels[i], " by ", rnm,
                     ", at values of ", cnm, .jn_group_suffix(x)),
      subtitle = if (isTRUE(x$binary[[r]]))
        paste0("Slopes at both values of ", rnm)
      else paste0("Shaded where the effect of ", x$labels[i],
                  " is significant (observed ", rnm, " only)"))
}
