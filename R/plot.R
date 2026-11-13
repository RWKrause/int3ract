# =============================================================================
# Figures. Built on demand from a finished JN object, never stored inside it.
#
# Every figure that has a moderator on an axis also carries a panel showing how
# many observations actually lie at each moderator value. A region of
# significance stretching over moderator values that were barely observed is an
# extrapolation, and the reader should be able to see that at a glance.
# =============================================================================

#' Appearance of Johnson-Neyman figures
#'
#' Collects the colour and pattern settings shared by every figure the package
#' draws, so that they can be set once and reused.
#'
#' @param sig_color,non_sig_color ribbon fill for the significant and
#'   non-significant parts of a two-way plot.
#' @param line_color colour of the conditional-effect line.
#' @param color_low,color_mid,color_high the diverging fill scale used by the
#'   heatmaps and by the two-way posterior densities. The defaults are taken
#'   from the \code{Zissou1} palette.
#' @param color_values colour of the within-cell value labels drawn on small
#'   heatmaps.
#' @param color_grid,grid_density,grid_spacing,grid_linewidth colour, density,
#'   spacing and line weight of the crosshatch drawn over non-significant
#'   heatmap cells. The defaults are chosen to read as texture over the fill
#'   rather than to hide it; widen \code{grid_spacing} for a coarser hatch.
#' @param crosshatch_non_sig crosshatch the non-significant cells
#'   (\code{TRUE}, the default) or the significant ones?
#' @param show_density draw the data-density panels? Defaults to \code{TRUE},
#'   and is ignored when the analysis carries no observed values.
#' @param density_fill,density_color fill and outline of the density panels.
#' @param density_height size of the density panels relative to the main panel.
#' @param density_bins number of histogram bins in the density panels.
#' @param extrap_alpha opacity of the parts of a figure lying beyond the
#'   observed range of the moderator, the margin added by the \code{extend}
#'   argument of \code{\link{JN}}. Lower is fainter; the observed part is drawn
#'   at 0.7 (bands) or 1 (heatmap tiles).
#' @param pos_color,neg_color line and shading colour of the simple-slopes
#'   figure (\code{plot(type = "slopes")}) where the focal effect is
#'   significantly positive, respectively negative. Default to
#'   \code{color_high} and \code{color_low}, the ends of the heatmap scale.
#' @param shade_alpha opacity of the shading between the simple slopes.
#' @param slope_shading \code{"flat"} (the default) shades the region between
#'   the simple slopes uniformly; \code{"density"} makes the shading follow how
#'   many observations lie at each moderator value, so that the part of a
#'   region where most of the data sit stands out. Needs the observed values.
#' @param mode_color colour of the simple slope drawn at the most common
#'   moderator value when the focal effect is significant throughout.
#'
#' @returns A list of settings, of class \code{jn_style}.
#'
#' @examples
#' set.seed(1)
#' dat <- data.frame(x = rnorm(100), z = rnorm(100))
#' dat$y <- dat$x + 0.5 * dat$x * dat$z + rnorm(100)
#' jn <- JN(lm(y ~ x * z, data = dat), theta_1 = "x", theta_2 = "z")
#'
#' plot(jn, style = jn_style(sig_color = "steelblue", show_density = FALSE))
#'
#' @export
jn_style <- function(sig_color = "seagreen3",
                     non_sig_color = "chocolate",
                     line_color = "black",
                     color_low = "#3B9AB2",
                     color_mid = "#EBCC2A",
                     color_high = "#F21A00",
                     color_values = "grey40",
                     color_grid = "black",
                     grid_density = 0.03,
                     grid_spacing = 0.05,
                     grid_linewidth = 0,
                     crosshatch_non_sig = TRUE,
                     show_density = TRUE,
                     density_fill = "grey70",
                     density_color = "grey30",
                     density_height = 0.28,
                     density_bins = 30,
                     extrap_alpha = 0.25,
                     pos_color = "seagreen3",
                     neg_color = "orangered4",
                     shade_alpha = 0.35,
                     slope_shading = c("flat", "density"),
                     mode_color = "grey30") {
  slope_shading <- match.arg(slope_shading)
  structure(as.list(environment()), class = "jn_style")
}

.jn_pos_color <- function(style) style$pos_color %||% style$color_high
.jn_neg_color <- function(style) style$neg_color %||% style$color_low

.jn_as_style <- function(style) {
  if (inherits(style, "jn_style")) return(style)
  if (is.list(style)) return(utils::modifyList(jn_style(), style))
  stop("`style` must be built with jn_style().", call. = FALSE)
}


# ---- density panels ---------------------------------------------------------

.jn_theme_density <- function() {
  ggplot2::theme_bw() +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   axis.title.y = ggplot2::element_text(size = 7),
                   axis.text.y  = ggplot2::element_text(size = 6))
}

# A histogram hanging downwards from the top of its panel, so that it reads as
# the base the plot above stands on rather than as a second plot.
.jn_density_bottom <- function(obs, rng, label, style, binary = FALSE) {
  bars <- if (binary) {
    ggplot2::geom_bar(width = 0.1 * diff(range(obs)),
                      fill = style$density_fill,
                      colour = style$density_color,
                      linewidth = 0.2)
  } else {
    ggplot2::geom_histogram(bins = style$density_bins,
                            fill = style$density_fill,
                            colour = style$density_color,
                            linewidth = 0.2)
  }
  ggplot2::ggplot(data.frame(value = obs), ggplot2::aes(x = .data$value)) +
    bars +
    ggplot2::scale_y_reverse(n.breaks = 3) +
    ggplot2::coord_cartesian(xlim = rng) +
    ggplot2::labs(x = label, y = "n") +
    .jn_theme_density()
}

# The same, standing upright; used above a heatmap.
.jn_density_top <- function(obs, rng, style) {
  ggplot2::ggplot(data.frame(value = obs), ggplot2::aes(x = .data$value)) +
    ggplot2::geom_histogram(bins = style$density_bins,
                            fill = style$density_fill,
                            colour = style$density_color,
                            linewidth = 0.2) +
    ggplot2::scale_y_continuous(n.breaks = 3) +
    ggplot2::coord_cartesian(xlim = rng) +
    ggplot2::labs(x = NULL, y = "n") +
    .jn_theme_density() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank())
}

# ... and lying on its side, used to the right of a heatmap.
.jn_density_right <- function(obs, rng, style) {
  ggplot2::ggplot(data.frame(value = obs), ggplot2::aes(y = .data$value)) +
    ggplot2::geom_histogram(bins = style$density_bins,
                            fill = style$density_fill,
                            colour = style$density_color,
                            linewidth = 0.2) +
    ggplot2::scale_x_continuous(n.breaks = 3) +
    ggplot2::coord_cartesian(ylim = rng) +
    ggplot2::labs(x = "n", y = NULL) +
    .jn_theme_density() +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank())
}

.jn_support_for <- function(x, label) {
  if (is.null(x$support)) return(NULL)
  i <- match(label, x$labels)
  if (is.na(i)) return(NULL)
  obs <- x$support[[i]]
  if (is.null(obs) || !length(obs)) NULL else obs
}


# ---- two-way ----------------------------------------------------------------

.jn_plot_2way_wald <- function(x, i, style) {
  tab <- x$tables[[i]]
  mod <- x$labels[3L - i]
  d   <- data.frame(mod_value = tab$mod_value, est = tab$estimate,
                    lo = tab$conf.low, hi = tab$conf.high,
                    significant = tab$significant,
                    extrapolated = tab$extrapolated)

  .jn_plot_effect(d, style, binary = isTRUE(x$binary[[mod]]),
                  obs_rng = x$observed_range[[mod]]) +
    ggplot2::labs(x = mod,
                  y = paste0("Conditional effect of ", x$labels[i]),
                  title = paste0("Moderation analysis for ", x$labels[i],
                                 .jn_group_suffix(x)))
}


.jn_plot_2way_post_band <- function(x, i, style) {
  tab <- x$tables[[i]]
  mod <- x$labels[3L - i]
  d   <- data.frame(mod_value = tab$mod_value, est = tab$post_mean,
                    lo = tab$q2.5, hi = tab$q97.5,
                    significant = tab$significant,
                    extrapolated = tab$extrapolated)

  .jn_plot_effect(d, style, binary = isTRUE(x$binary[[mod]]),
                  obs_rng = x$observed_range[[mod]]) +
    ggplot2::labs(x = mod,
                  y = paste0("Conditional effect of ", x$labels[i]),
                  title = paste0("Conditional posterior for ", x$labels[i],
                                 .jn_group_suffix(x)))
}


# The conditional effect against its moderator, from a data frame holding the
# moderator value, the estimate, its interval, and the significance and
# extrapolation flags; an optional `panel` column facets it. A binary moderator
# gets the estimate and its interval at each of its two values; anything else a
# band shaded by significance and drawn fainter beyond the observed values.
.jn_plot_effect <- function(d, style, binary, obs_rng = NULL) {
  p <- ggplot2::ggplot(d, ggplot2::aes(.data$mod_value, .data$est)) +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::theme_bw()
  if (!is.null(d$panel))
    p <- p + ggplot2::facet_wrap(ggplot2::vars(.data$panel), nrow = 1)

  if (binary) {
    d$sig_label <- factor(ifelse(d$significant, "significant",
                                 "not significant"),
                          levels = c("significant", "not significant"))
    lv <- sort(unique(d$mod_value))
    return(p +
      ggplot2::geom_errorbar(data = d,
                             ggplot2::aes(ymin = .data$lo, ymax = .data$hi,
                                          colour = .data$sig_label),
                             width = 0.08 * diff(range(lv)), linewidth = 0.8) +
      ggplot2::geom_point(data = d, ggplot2::aes(colour = .data$sig_label),
                          size = 2.5) +
      ggplot2::scale_colour_manual(
        values = c(significant = style$sig_color,
                   `not significant` = style$non_sig_color),
        drop = FALSE, name = NULL) +
      ggplot2::scale_x_continuous(breaks = lv) +
      ggplot2::expand_limits(x = lv + c(-0.25, 0.25) * diff(lv)) +
      ggplot2::theme(legend.position = "bottom"))
  }

  # One ribbon per run of constant significance and extrapolation within a
  # panel, so that colour and opacity change exactly at the boundary rather
  # than being blended across it. Each run reaches to the first point of the
  # next, so that coarse grids leave no gap.
  panel  <- if (is.null(d$panel)) rep("", nrow(d)) else as.character(d$panel)
  runs   <- rle(paste(panel, d$significant, d$extrapolated))
  ends   <- cumsum(runs$lengths)
  starts <- c(1L, utils::head(ends, -1L) + 1L)

  for (r in seq_along(runs$values)) {
    s <- starts[r]; e <- ends[r]
    ext <- d$extrapolated[s]
    sig <- d$significant[s]
    if (e < nrow(d) && panel[e + 1L] == panel[e]) e <- e + 1L
    # a single point at the end of a panel reaches back instead
    if (e == s && s > 1L && panel[s - 1L] == panel[s]) s <- s - 1L
    run <- d[s:e, , drop = FALSE]
    p <- p +
      ggplot2::geom_ribbon(
        data = run,
        ggplot2::aes(ymin = .data$lo, ymax = .data$hi),
        fill  = if (sig) style$sig_color else style$non_sig_color,
        alpha = if (ext) style$extrap_alpha else 0.7) +
      ggplot2::geom_path(data = run, linewidth = 1,
                         colour = style$line_color,
                         linetype = if (ext) "22" else "solid")
  }

  # mark where the data end, when the figure goes beyond them
  if (any(d$extrapolated) && !is.null(obs_rng))
    p <- p + ggplot2::geom_vline(xintercept = obs_rng, linetype = "dotted",
                                 colour = "grey40")
  p
}


.jn_plot_2way_post_density <- function(x, i, style) {
  cond <- x$cond_draws[[i]]
  m    <- x$tables[[i]]$mod_value
  mod  <- x$labels[3L - i]

  d <- data.frame(
    parameter = as.vector(cond),
    modValue  = rep(m, each = nrow(cond)))

  leg <- if (nchar(mod) < 26) mod else "Moderator"

  ggplot2::ggplot(d, ggplot2::aes(x = .data$parameter,
                                  group = .data$modValue,
                                  fill  = .data$modValue,
                                  colour = .data$modValue)) +
    ggplot2::geom_vline(xintercept = 0) +
    ggplot2::geom_density(alpha = 0.1) +
    ggplot2::scale_fill_gradient2(low = style$color_low, mid = style$color_mid,
                                  high = style$color_high, midpoint = 0) +
    ggplot2::scale_color_gradient2(low = style$color_low, mid = style$color_mid,
                                   high = style$color_high, midpoint = 0) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title  = paste0("Posterior density for ", x$labels[i],
                      "\nmoderated by ", mod, .jn_group_suffix(x)),
      x = x$labels[i], y = "Posterior density",
      colour = leg, fill = leg)
}


# ---- three-way --------------------------------------------------------------

.jn_plot_3way <- function(x, i, fill_var, midpoint, fill_label, title, style) {
  tab <- x$tables[[i]]
  pat <- tab[tab$significant != style$crosshatch_non_sig, , drop = FALSE]
  pat$pattern <- "crosshatch"

  p <- ggplot2::ggplot(tab, ggplot2::aes(.data$mod1_value, .data$mod2_value)) +
    ggplot2::geom_tile(ggplot2::aes(fill = .data[[fill_var]],
                                    alpha = .data$extrapolated)) +
    ggplot2::scale_alpha_manual(values = c(`FALSE` = 1,
                                           `TRUE` = style$extrap_alpha),
                                guide = "none") +
    ggplot2::scale_color_identity() +
    ggplot2::scale_fill_gradient2(low = style$color_low, mid = style$color_mid,
                                  high = style$color_high,
                                  midpoint = midpoint) +
    ggplot2::theme_bw() +
    ggplot2::guides(pattern = "none") +
    ggplot2::labs(title = title, fill = fill_label,
                  x = tab$mod1[1], y = tab$mod2[1])

  if (nrow(pat))
    p <- p + ggpattern::geom_tile_pattern(
      data = pat, ggplot2::aes(pattern = .data$pattern),
      pattern_density = style$grid_density,
      pattern_spacing = style$grid_spacing,
      # ggpattern strokes the outline of each hatch polygon in pattern_color
      # and fills its interior with pattern_fill. Setting only the former
      # leaves pale interiors ringed in black; both are set to the same colour
      # so the hatch reads as a plain line of the requested weight.
      pattern_color   = style$color_grid,
      pattern_fill    = style$color_grid,
      pattern_size    = style$grid_linewidth,
      alpha = 0)

  # outline the observed part of the grid when the figure goes beyond it
  if (any(tab$extrapolated)) {
    r1 <- x$observed_range[[tab$mod1[1]]]
    r2 <- x$observed_range[[tab$mod2[1]]]
    p  <- p + ggplot2::annotate("rect", xmin = r1[1], xmax = r1[2],
                                ymin = r2[1], ymax = r2[2], fill = NA,
                                colour = "grey20", linetype = "dashed",
                                linewidth = 0.4)
  }

  # print the value in each cell only when the grid is coarse enough to read
  if (length(unique(tab$mod1_value)) < 8 &&
      length(unique(tab$mod2_value)) < 8)
    p <- p + ggplot2::geom_text(
      ggplot2::aes(label = round(.data[[fill_var]], 3)),
      colour = style$color_values, size = 2.5)

  p
}


# A three-way analysis in which a moderator is binary. A heatmap along an axis
# with two values is no heatmap, so the two-way figure along the other
# moderator is drawn once for each level of the binary one.
.jn_plot_3way_binary <- function(x, i, style) {
  tab  <- x$tables[[i]]
  role <- .jn_3way_roles(x, i)
  cv   <- paste0(role$cond, "_value")
  rv   <- paste0(role$run,  "_value")
  cnm  <- tab[[role$cond]][1]
  rnm  <- tab[[role$run]][1]

  if (x$inference == "wald") {
    z   <- abs(stats::qnorm(x$alpha / 2))
    est <- tab$estimate
    lo  <- est - z * tab$std.error
    hi  <- est + z * tab$std.error
  } else {
    est <- tab$post_mean; lo <- tab$q2.5; hi <- tab$q97.5
  }
  d <- data.frame(mod_value = tab[[rv]], est = est, lo = lo, hi = hi,
                  significant = tab$significant,
                  extrapolated = tab$extrapolated,
                  level = tab[[cv]])
  d  <- d[order(d$level, d$mod_value), , drop = FALSE]
  lv <- sort(unique(d$level))
  d$panel <- factor(paste(cnm, "=", .jn_fmt_value(d$level)),
                    levels = paste(cnm, "=", .jn_fmt_value(lv)))

  .jn_plot_effect(d, style, binary = isTRUE(x$binary[[rnm]]),
                  obs_rng = x$observed_range[[rnm]]) +
    ggplot2::labs(x = rnm,
                  y = paste0(if (x$inference == "wald") "Conditional effect"
                             else "Posterior mean", " of ", x$labels[i]),
                  title = paste0(x$labels[i], .jn_group_suffix(x)))
}


# ---- composition ------------------------------------------------------------

.jn_group_suffix <- function(x) {
  if (is.null(x$group)) "" else paste0(" (", x$group, ")")
}

# Attach the density panel(s) if there is anything to draw, otherwise hand back
# the bare ggplot. Both print and ggsave() the same way.
.jn_compose_2way <- function(p, x, i, style, aligned) {
  mod <- x$labels[3L - i]
  obs <- .jn_support_for(x, mod)
  if (!style$show_density || is.null(obs)) return(p)

  binary <- isTRUE(x$binary[[mod]])
  rng    <- range(x$tables[[i]]$mod_value)
  if (binary) rng <- rng + c(-0.25, 0.25) * diff(rng)
  dens   <- .jn_density_bottom(obs, rng, mod, style, binary = binary)

  if (aligned) {
    # the main panel already carries the moderator on its x axis, so drop its
    # axis and let the density panel below label the shared scale
    p <- p + ggplot2::coord_cartesian(xlim = rng) +
      ggplot2::labs(x = NULL) +
      ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                     axis.ticks.x = ggplot2::element_blank())
  } else {
    # the main panel shows the parameter, not the moderator; the panel below is
    # a separate chart and says so
    dens <- dens + ggplot2::labs(x = paste0("Distribution of ", mod))
  }

  patchwork::wrap_plots(p, dens, ncol = 1,
                        heights = c(1, style$density_height))
}

.jn_compose_3way <- function(p, x, i, style) {
  tab  <- x$tables[[i]]
  obs1 <- .jn_support_for(x, tab$mod1[1])
  obs2 <- .jn_support_for(x, tab$mod2[1])
  if (!style$show_density || (is.null(obs1) && is.null(obs2))) return(p)

  rng1 <- range(tab$mod1_value)
  rng2 <- range(tab$mod2_value)

  # The title moves to the composition, otherwise it is stranded between the
  # top histogram and the heatmap; the legend moves underneath, otherwise it
  # wedges itself between the heatmap and the right-hand histogram.
  title <- p$labels$title
  p <- p +
    ggplot2::coord_cartesian(xlim = rng1, ylim = rng2) +
    ggplot2::labs(title = NULL) +
    ggplot2::theme(legend.position = "bottom")

  top   <- if (is.null(obs1)) patchwork::plot_spacer() else
             .jn_density_top(obs1, rng1, style)
  right <- if (is.null(obs2)) patchwork::plot_spacer() else
             .jn_density_right(obs2, rng2, style)

  patchwork::wrap_plots(top, patchwork::plot_spacer(), p, right,
                        ncol = 2,
                        widths  = c(1, style$density_height),
                        heights = c(style$density_height, 1)) +
    patchwork::plot_annotation(title = title)
}


# ---- the plot methods -------------------------------------------------------

#' Plot a Johnson-Neyman analysis
#'
#' Draws the figure for one focal variable. Two-way Wald analyses give the
#' classic plot of the conditional effect against the moderator, with the
#' confidence band shaded by significance; three-way analyses give a heatmap
#' over the two-dimensional moderator grid, crosshatched where the effect is
#' not significant; two-way posterior analyses give either overlaid conditional
#' posterior densities or, with \code{type = "band"}, the posterior mean and
#' credible band against the moderator. Where a moderator is binary the
#' two-way figure shows the effect with its interval at the two values, and a
#' three-way heatmap becomes the two-way figure along the other moderator, once
#' for each value of the binary one.
#'
#' With \code{type = "slopes"}, any analysis is drawn instead in outcome space:
#' the focal variable on the x axis, the model's prediction of the outcome on
#' the y axis, and simple slopes at the moderator values bounding each region
#' of significance, the fan between them shaded in \code{pos_color} where the
#' effect is significantly positive and \code{neg_color} where it is
#' significantly negative. A focal effect that is nowhere significant gets one
#' dashed slope at the moderator mean; one that is significant throughout gets
#' the slopes at the lowest and highest observed moderator values plus a
#' dotted one at the most common value. A binary moderator gets its two slopes,
#' each marked as significant or not. The regions are computed over the
#' observed moderator values only, so the fan never reaches into the margin
#' added by \code{extend}. For three-way analyses the figure is faceted: one
#' panel per value of one moderator (a binary one if there is one, otherwise
#' the first, at \code{at}), each showing the two-way figure along the other.
#' For Wald analyses the boundaries of every panel are solved for exactly. The
#' outcome axis is on the model's scale -- the link scale for \code{glm} --
#' with other covariates at their means, and is drawn up to an additive
#' constant when the model does not supply its intercept (see
#' \code{\link{jn_input}}).
#'
#' When the observed values of the moderators are known -- automatically for
#' \code{lm}, \code{glm} and \pkg{lme4} models, and otherwise from the
#' \code{support} argument of \code{\link{JN}} -- a histogram of them is
#' attached to the figure: below the panel for two-way plots, along the top and
#' right edges for heatmaps. Regions of significance covering moderator values
#' that were hardly observed are the standard failure mode of this technique,
#' and these panels make them visible. For heatmaps the histograms are
#' marginal, so a cell can look well supported on both axes even though few
#' observations lie near that combination of moderator values; the joint
#' distribution is not shown.
#'
#' @param x a \code{JN} object.
#' @param which the focal variable(s) to plot, by name or position. Defaults to
#'   all of them: every variable involved in the interaction takes its turn as
#'   the focal one, so a two-way analysis has two figures and a three-way
#'   analysis three.
#' @param type for two-way posterior analyses, \code{"density"} (the default)
#'   for overlaid conditional posteriors or \code{"band"} for the posterior
#'   mean with its credible band; for three-way posterior analyses,
#'   \code{"estimate"} (the default) for the posterior mean or \code{"bayes_p"}
#'   for the Bayesian \eqn{p} value. \code{"slopes"} draws the simple-slopes
#'   figure for any analysis.
#' @param style appearance settings from \code{\link{jn_style}}.
#' @param at for three-way simple-slopes figures, the values of the moderator
#'   to draw one panel at. Defaults to its quartiles; ignored when a moderator
#'   is binary, whose two values are used.
#' @param ... ignored.
#'
#' @returns Invisibly, the figure when a single \code{which} was given, and a
#'   named list of them otherwise. Each figure is a \pkg{ggplot2} object, or a
#'   \pkg{patchwork} composition of one where density panels are attached; both
#'   print, extend with \code{+} and save with \code{ggsave()} in the usual way.
#'
#'   \code{plot()} draws; \code{\link{jn_plots}} and \code{autoplot()} return
#'   figures without drawing them.
#'
#' @examples
#' set.seed(1)
#' dat <- data.frame(x = rnorm(200), z = rnorm(200))
#' dat$y <- dat$x + 0.6 * dat$x * dat$z + rnorm(200)
#' jn <- JN(lm(y ~ x * z, data = dat), theta_1 = "x", theta_2 = "z")
#'
#' plot(jn)                 # both figures: x moderated by z, and z by x
#' plot(jn, which = "x")     # just the one
#' plot(jn, style = jn_style(show_density = FALSE))
#'
#' \donttest{
#' # the regions of significance as simple slopes of the outcome
#' plot(jn, which = "x", type = "slopes")
#' plot(jn, which = "x", type = "slopes",
#'      style = jn_style(slope_shading = "density"))
#' }
#'
#' @export
plot.JN <- function(x, which = NULL, type = NULL, style = jn_style(),
                    at = NULL, ...) {
  style <- .jn_as_style(style)
  idx   <- if (is.null(which)) seq_along(x$labels)
           else vapply(which, function(w) .jn_which(x, w), integer(1))

  figs <- stats::setNames(
    lapply(idx, function(i) .jn_figure(x, i, type, style, at)),
    x$labels[idx])

  # More than one figure and a screen device: page through them, as plot.lm
  # does for its diagnostic plots.
  if (length(figs) > 1L && grDevices::dev.interactive()) {
    ask <- grDevices::devAskNewPage(TRUE)
    on.exit(grDevices::devAskNewPage(ask), add = TRUE)
  }
  for (p in figs) print(p)

  invisible(if (length(figs) == 1L) figs[[1L]] else figs)
}


# Builds one figure. Everything about which figure to build lives here; whether
# to draw it is the caller's business.
.jn_figure <- function(x, i, type, style, at = NULL) {
  if (identical(type, "slopes"))
    return(if (x$n_way == 2L) .jn_plot_2way_slopes(x, i, style)
           else .jn_plot_3way_slopes(x, i, style, at))

  if (x$n_way == 3L) {
    tab <- x$tables[[i]]
    if (isTRUE(x$binary[[tab$mod1[1]]]) || isTRUE(x$binary[[tab$mod2[1]]]))
      return(.jn_plot_3way_binary(x, i, style))
  }

  if (x$n_way == 2L) {
    if (x$inference == "wald") {
      .jn_compose_2way(.jn_plot_2way_wald(x, i, style), x, i, style,
                       aligned = TRUE)
    } else {
      type <- match.arg(type %||% "density", c("density", "band"))
      if (type == "band")
        .jn_compose_2way(.jn_plot_2way_post_band(x, i, style), x, i, style,
                         aligned = TRUE)
      else
        .jn_compose_2way(.jn_plot_2way_post_density(x, i, style), x, i, style,
                         aligned = FALSE)
    }
  } else {
    p <- if (x$inference == "wald") {
      .jn_plot_3way(x, i, "estimate", 0, "Conditional effect",
                    paste0(x$labels[i], .jn_group_suffix(x)), style)
    } else {
      type <- match.arg(type %||% "estimate", c("estimate", "bayes_p"))
      if (type == "bayes_p")
        .jn_plot_3way(x, i, "bayes_p", 0.5,
                      expression("Bayesian" ~ italic("p") * " value"),
                      paste0("Bayesian p value for ", x$labels[i],
                             .jn_group_suffix(x)), style)
      else
        .jn_plot_3way(x, i, "post_mean", 0, "Posterior mean",
                      paste0("Posterior mean for ", x$labels[i],
                             .jn_group_suffix(x)), style)
    }
    .jn_compose_3way(p, x, i, style)
  }
}


# A package that supplies an autoplot() method should re-export the generic;
# otherwise autoplot() is only reachable after attaching ggplot2 as well.
#' @importFrom ggplot2 autoplot
#' @export
ggplot2::autoplot


#' @rdname plot.JN
#' @param object a \code{JN} object.
#' @export
autoplot.JN <- function(object, which = 1L, type = NULL,
                        style = jn_style(), at = NULL, ...) {
  # autoplot() returns a figure rather than drawing one, so it takes a single
  # focal variable; use jn_plots() for all of them at once.
  .jn_figure(object, .jn_which(object, which), type, .jn_as_style(style), at)
}


#' @rdname plot.JN
#' @export
plot.JN_list <- function(x, which = NULL, ...) {
  # A grouped analysis holds one result per group, so drawing everything can
  # mean dozens of figures. Draw the population-level one and say where the
  # rest are.
  if (length(x) > 1L)
    message("Drawing the ", names(x)[1], " analysis only. ",
            "Use plot(x[[i]]) for one group, or jn_save(x) to write every ",
            "figure to disk.")
  invisible(plot(x[[1L]], which = which, ...))
}


#' Every figure of a Johnson-Neyman analysis
#'
#' Returns all of the figures at once, named by focal variable, without
#' drawing any of them. \code{\link{plot.JN}} draws instead; use this when the
#' figures are wanted as objects, to arrange, modify or save.
#'
#' @param x a \code{JN} object.
#' @param type which figure to draw; see \code{\link{plot.JN}}. For three-way
#'   posterior analyses both quantities are returned unless
#'   \code{type = "slopes"}.
#' @param style appearance settings from \code{\link{jn_style}}.
#' @param at passed to \code{\link{plot.JN}}.
#' @param ... ignored.
#'
#' @returns A named list of figures. For three-way posterior analyses the list
#'   has two elements, \code{post_mean} and \code{bayes_p}, each holding one
#'   figure per focal variable.
#'
#' @examples
#' set.seed(1)
#' dat <- data.frame(x = rnorm(100), z = rnorm(100))
#' dat$y <- dat$x + 0.5 * dat$x * dat$z + rnorm(100)
#' figs <- jn_plots(JN(lm(y ~ x * z, data = dat),
#'                     theta_1 = "x", theta_2 = "z"))
#' names(figs)
#'
#' # the simple-slopes figures instead
#' slopes <- jn_plots(JN(lm(y ~ x * z, data = dat),
#'                       theta_1 = "x", theta_2 = "z"), type = "slopes")
#'
#' @export
jn_plots <- function(x, type = NULL, style = jn_style(), at = NULL, ...) {
  stopifnot(inherits(x, "JN"))
  style <- .jn_as_style(style)
  one   <- function(tp)
    stats::setNames(lapply(seq_along(x$labels),
                           function(i) .jn_figure(x, i, tp, style, at)),
                    x$labels)

  if (x$n_way == 3L && x$inference == "posterior" &&
      !identical(type, "slopes"))
    list(post_mean = one("estimate"), bayes_p = one("bayes_p"))
  else
    one(type)
}


#' Save the figures of a Johnson-Neyman analysis
#'
#' Writes every figure to disk with \code{ggplot2::ggsave()}. Useful for
#' grouped analyses, where inspecting dozens of figures interactively is
#' impractical.
#'
#' @param x a \code{JN} or \code{JN_list} object.
#' @param folder directory to write into; created if it does not exist.
#'   Defaults to a session-temporary directory.
#' @param device file extension passed to \code{ggsave()}. Default
#'   \code{"png"}.
#' @param width,height,dpi passed to \code{ggsave()}.
#' @param ... passed to \code{\link{jn_plots}}.
#'
#' @returns The paths written, invisibly.
#'
#' @examples
#' set.seed(1)
#' dat <- data.frame(x = rnorm(100), z = rnorm(100))
#' dat$y <- dat$x + 0.5 * dat$x * dat$z + rnorm(100)
#' jn <- JN(lm(y ~ x * z, data = dat), theta_1 = "x", theta_2 = "z")
#' jn_save(jn, folder = tempfile("jn"))
#'
#' @export
jn_save <- function(x, folder = NULL, device = "png",
                    width = 10, height = 7, dpi = 600, ...) {
  folder <- folder %||% file.path(tempdir(), "int3ract-figures")
  if (!dir.exists(folder)) dir.create(folder, recursive = TRUE)

  figs <- if (inherits(x, "JN_list")) {
    unlist(lapply(names(x), function(g) {
      f <- jn_plots(x[[g]], ...)
      stats::setNames(f, paste0(g, "-", names(f)))
    }), recursive = FALSE)
  } else {
    jn_plots(x, ...)
  }
  figs <- .jn_flatten(figs)

  paths <- vapply(names(figs), function(nm) {
    path <- file.path(folder, paste0(.jn_safe_name(nm), ".", device))
    ggplot2::ggsave(path, plot = figs[[nm]],
                    width = width, height = height, dpi = dpi)
    path
  }, character(1))

  invisible(unname(paths))
}

.jn_flatten <- function(x) {
  if (!is.list(x) || inherits(x, "ggplot") || inherits(x, "patchwork"))
    return(list(x))
  out <- list()
  for (nm in names(x)) {
    part <- .jn_flatten(x[[nm]])
    names(part) <- if (length(part) == 1L && is.null(names(x[[nm]]))) nm
                   else paste0(nm, "-", names(part))
    out <- c(out, part)
  }
  out
}

.jn_safe_name <- function(x) gsub("[^A-Za-z0-9._-]+", "_", x)

.jn_which <- function(x, which) {
  if (is.character(which)) {
    i <- match(which, x$labels)
    if (is.na(i))
      stop("`which` must be one of ",
           paste(sQuote(x$labels), collapse = ", "), ".", call. = FALSE)
    return(i)
  }
  i <- as.integer(which)
  if (is.na(i) || i < 1L || i > length(x$labels))
    stop("`which` must be between 1 and ", length(x$labels), ".",
         call. = FALSE)
  i
}
