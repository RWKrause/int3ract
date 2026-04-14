# =============================================================================
# Shared utilities
# =============================================================================

# ---- used by: JNK_bayes (multiSiena path) ----

.clean_effect_name <- function(n) {
  gsub(':', '', gsub('/', 'o', gsub('\\^', '', n)))
}


# ---- used by: .freq_core (three-way) ----

parFun <- function(m1,m2,
                   bx,
                   bm1x,
                   bm2x,
                   bm1m2x) {
  bx + m1 * bm1x  + m2 * bm2x  + m1 * m2 * bm1m2x
}

seFun <- function(m1,m2,
                  Vx,
                  Vm1x,
                  Vm2x,
                  Vm1m2x,
                  covX_m1x,
                  covX_m2x,
                  covX_m1m2x,
                  covm1x_m2x,
                  covm1x_m1m2x,
                  covm2x_m1m2x) {
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


# =============================================================================
# Internal support for JNK_freq
# =============================================================================

# ---- used by: JNK_freq, .extract_lm, .extract_lme4, .lme4_groups ----

.find_int <- function(nms, ...) {
  parts <- c(...)
  if (length(parts) == 2) {
    candidates <- c(paste(parts, collapse = ':'),
                    paste(rev(parts), collapse = ':'))
  } else {
    p <- parts
    candidates <- c(
      paste(p[1], p[2], p[3], sep = ':'),
      paste(p[1], p[3], p[2], sep = ':'),
      paste(p[2], p[1], p[3], sep = ':'),
      paste(p[2], p[3], p[1], sep = ':'),
      paste(p[3], p[1], p[2], sep = ':'),
      paste(p[3], p[2], p[1], sep = ':'))
  }
  found <- candidates[candidates %in% nms]
  if (length(found) == 0)
    stop("Interaction term not found for: ", paste(parts, collapse = " * "),
         "\nAvailable: ", paste(nms, collapse = ", "))
  found[1]
}


# ---- used by: JNK_freq (lm/glm path) ----

.extract_lm <- function(x, theta_1, theta_2, theta_3,
                         theta_1_vals, theta_2_vals, theta_3_vals) {
  cn   <- names(coef(x))
  covT <- vcov(x)
  threeWay <- !is.null(theta_3)

  if (!threeWay) {
    int_12 <- .find_int(cn, theta_1, theta_2)
    idx    <- c(theta_1, theta_2, int_12)
    list(name  = c(theta_1, theta_2),
         coefs = coef(x)[idx],
         covar = covT[idx, idx],
         vals  = list(theta_1_vals %||% range(x$model[[theta_1]], na.rm = TRUE),
                      theta_2_vals %||% range(x$model[[theta_2]], na.rm = TRUE)))
  } else {
    int_12  <- .find_int(cn, theta_1, theta_2)
    int_13  <- .find_int(cn, theta_1, theta_3)
    int_23  <- .find_int(cn, theta_2, theta_3)
    int_123 <- .find_int(cn, theta_1, theta_2, theta_3)
    idx     <- c(theta_1, theta_2, theta_3, int_12, int_13, int_23, int_123)
    list(name  = c(theta_1, theta_2, theta_3),
         coefs = coef(x)[idx],
         covar = covT[idx, idx],
         vals  = list(theta_1_vals %||% range(x$model[[theta_1]], na.rm = TRUE),
                      theta_2_vals %||% range(x$model[[theta_2]], na.rm = TRUE),
                      theta_3_vals %||% range(x$model[[theta_3]], na.rm = TRUE)))
  }
}


# ---- used by: JNK_freq (sienaFit path) ----

.extract_siena <- function(x, theta_1, theta_2, theta_3,
                            theta_int_12, theta_int_13,
                            theta_int_23, theta_int_123,
                            theta_1_vals, theta_2_vals, theta_3_vals) {
  sn   <- x$effects$effectName
  covT <- x$covtheta
  threeWay <- !is.null(theta_3)

  all_idx <- c(theta_1, theta_2, theta_3,
               theta_int_12, theta_int_13, theta_int_23, theta_int_123)
  if (any(all_idx > length(sn))) {
    nms <- c('theta_1','theta_2','theta_3','theta_int_12',
             'theta_int_13','theta_int_23','theta_int_123')
    bad <- nms[seq_along(all_idx)][all_idx > length(sn)]
    stop("Parameter indices out of range: ", paste(bad, collapse = ", "),
         "\nSee x$effects$effectName for valid indices.")
  }

  if (!threeWay) {
    idx <- c(theta_1, theta_2, theta_int_12)
    list(name  = c(sn[theta_1], sn[theta_2]),
         coefs = x$theta[idx],
         covar = covT[idx, idx],
         vals  = list(theta_1_vals, theta_2_vals))
  } else {
    idx <- c(theta_1, theta_2, theta_3,
             theta_int_12, theta_int_13, theta_int_23, theta_int_123)
    list(name  = c(sn[theta_1], sn[theta_2], sn[theta_3]),
         coefs = x$theta[idx],
         covar = covT[idx, idx],
         vals  = list(theta_1_vals, theta_2_vals, theta_3_vals))
  }
}


# ---- used by: JNK_freq (lme4 path) ----

.extract_lme4 <- function(x, theta_1, theta_2, theta_3,
                           theta_1_vals, theta_2_vals, theta_3_vals) {
  fe   <- lme4::fixef(x)
  covT <- as.matrix(vcov(x))
  mf   <- model.frame(x)
  cn   <- names(fe)
  threeWay <- !is.null(theta_3)

  if (!threeWay) {
    int_12 <- .find_int(cn, theta_1, theta_2)
    idx    <- c(theta_1, theta_2, int_12)
    list(name  = c(theta_1, theta_2),
         coefs = fe[idx],
         covar = covT[idx, idx],
         vals  = list(theta_1_vals %||% range(mf[[theta_1]], na.rm = TRUE),
                      theta_2_vals %||% range(mf[[theta_2]], na.rm = TRUE)))
  } else {
    int_12  <- .find_int(cn, theta_1, theta_2)
    int_13  <- .find_int(cn, theta_1, theta_3)
    int_23  <- .find_int(cn, theta_2, theta_3)
    int_123 <- .find_int(cn, theta_1, theta_2, theta_3)
    idx     <- c(theta_1, theta_2, theta_3, int_12, int_13, int_23, int_123)
    list(name  = c(theta_1, theta_2, theta_3),
         coefs = fe[idx],
         covar = covT[idx, idx],
         vals  = list(theta_1_vals %||% range(mf[[theta_1]], na.rm = TRUE),
                      theta_2_vals %||% range(mf[[theta_2]], na.rm = TRUE),
                      theta_3_vals %||% range(mf[[theta_3]], na.rm = TRUE)))
  }
}


# ---- used by: JNK_freq (lme4 group-level path) ----

.lme4_groups <- function(x, theta_1, theta_2, theta_3,
                          out, group_var, plot_args) {
  if (is.null(group_var)) {
    group_var <- names(lme4::ranef(x))[1]
    message("Using grouping variable: ", group_var)
  }

  group_coefs <- coef(x)[[group_var]]
  group_names <- rownames(group_coefs)
  cn_all      <- colnames(group_coefs)
  threeWay    <- !is.null(theta_3)

  if (!threeWay) {
    int_12 <- .find_int(cn_all, theta_1, theta_2)
    cn <- c(theta_1, theta_2, int_12)
  } else {
    int_12  <- .find_int(cn_all, theta_1, theta_2)
    int_13  <- .find_int(cn_all, theta_1, theta_3)
    int_23  <- .find_int(cn_all, theta_2, theta_3)
    int_123 <- .find_int(cn_all, theta_1, theta_2, theta_3)
    cn <- c(theta_1, theta_2, theta_3, int_12, int_13, int_23, int_123)
  }

  missing_cn <- cn[!cn %in% cn_all]
  if (length(missing_cn) > 0) {
    warning("Not all interaction terms have random effects. ",
            "Using fixed effects for: ", paste(missing_cn, collapse = ", "))
  }

  setNames(
    lapply(group_names, function(g) {
      g_coefs <- numeric(length(cn))
      names(g_coefs) <- cn
      for (j in seq_along(cn)) {
        g_coefs[j] <- if (cn[j] %in% cn_all) {
          group_coefs[g, cn[j]]
        } else {
          out$coefs[j]
        }
      }
      do.call(.freq_core,
              c(list(name = out$name, covar = out$covar,
                     coefs = g_coefs, vals = out$vals,
                     group = g),
                plot_args))
    }),
    group_names
  )
}


# ---- used by: JNK_freq (via .lme4_groups and directly) ----

.freq_core <- function(name, covar, coefs, vals,
                       group = NULL,
                       range_size,
                       alpha = 0.05,
                       round_res = 3,
                       control_fdr = FALSE,
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
                       folder = 'int3ract JNKplots') {

  k      <- length(name)
  z_crit <- abs(qnorm(alpha / 2))

  grp_label <- if (!is.null(group)) paste0(" (", group, ")") else ""

  for (var in seq_len(k))
    vals[[var]] <- seq(min(vals[[var]]),
                       max(vals[[var]]),
                       length.out = range_size)

  # ---- two-way ----
  if (k == 2) {

    add_ribbons <- function(p, df) {
      runs   <- rle(df$sig)
      ends   <- cumsum(runs$lengths)
      starts <- c(1, head(ends, -1) + 1)
      for (r in seq_along(runs$values)) {
        p <- p + ggplot2::geom_ribbon(
          data = df[starts[r]:ends[r], ],
          ggplot2::aes(ymin = theta_vals - z_crit * theta_se,
                       ymax = theta_vals + z_crit * theta_se),
          fill  = if (runs$values[r]) sig_color else non_sig_color,
          alpha = 0.7)
      }
      p
    }

    td <- lapply(1:k, function(var) {
      mv      <- c(1:k)[-var]
      theta1s <- coefs[var] + coefs[3] * vals[[mv]]
      seT1    <- sqrt(covar[var, var] +
                        vals[[mv]] * 2 * covar[mv, var] +
                        vals[[mv]]^2 * covar[3, 3])
      p1      <- 2 * pmin(pnorm(theta1s / seT1), 1 - pnorm(theta1s / seT1))

      sig11 <- if (control_fdr) {
        ord     <- order(p1)
        sig_ord <- p1[ord] < alpha * seq_along(p1) / length(p1)
        if (!all(sig_ord)) sig_ord[which.max(!sig_ord):length(p1)] <- FALSE
        out        <- logical(length(p1))
        out[ord]   <- sig_ord
        out
      } else {
        p1 < alpha
      }

      data.frame(theta      = name[var],
                 moderator  = name[mv],
                 mod_vals   = round(vals[[mv]],   round_res),
                 theta_vals = round(theta1s,       round_res),
                 theta_se   = round(seT1,          round_res),
                 theta_p    = round(p1, 3),
                 sig        = sig11)
    })
    names(td) <- name

    plots <- lapply(1:2, function(i) {
      p <- ggplot2::ggplot(td[[i]], ggplot2::aes(mod_vals, theta_vals)) +
        ggplot2::geom_hline(ggplot2::aes(yintercept = 0)) +
        ggplot2::theme_bw() +
        ggplot2::ylab(paste0("Theta for ", name[i])) +
        ggplot2::xlab(name[-i]) +
        ggplot2::ggtitle(paste0("Moderation analysis for ", name[i],
                                grp_label))

      p <- add_ribbons(p, td[[i]]) +
        ggplot2::geom_path(linewidth = 1, color = line_color)

      if (save) {
        if (!dir.exists(folder)) dir.create(folder, recursive = TRUE)
        ggplot2::ggsave(
          file.path(folder, paste0("Moderation analysis for ", name[i],
                                   grp_label, ".png")),
          plot = p, dpi = 600, width = 10, height = 7)
      }
      p
    })
    names(plots) <- name

    return(list(param_table = td, plots = plots))
  }

  # ---- three-way ----
  x_idx  <- c(1L, 2L, 3L)
  m1_idx <- c(4L, 4L, 5L)
  m2_idx <- c(5L, 6L, 6L)
  v_idx  <- list(c(2L, 3L), c(1L, 3L), c(1L, 2L))

  mat_to_long <- function(mat, val_name) {
    mat |>
      as.data.frame() |>
      tibble::rownames_to_column("Var1") |>
      tidyr::pivot_longer(-Var1, names_to = "Var2", values_to = val_name) |>
      dplyr::mutate(Var1 = as.numeric(Var1), Var2 = as.numeric(Var2))
  }

  mats <- lapply(1:3, function(i) {
    xi  <- x_idx[i];  m1i <- m1_idx[i];  m2i <- m2_idx[i]
    vi  <- v_idx[[i]]
    th  <- outer(vals[[vi[1]]], vals[[vi[2]]], parFun,
                 bx = coefs[xi], bm1x = coefs[m1i],
                 bm2x = coefs[m2i], bm1m2x = coefs[7])
    se  <- outer(vals[[vi[1]]], vals[[vi[2]]], seFun,
                 Vx           = covar[xi,  xi],
                 Vm1x         = covar[m1i, m1i],
                 Vm2x         = covar[m2i, m2i],
                 Vm1m2x       = covar[7,   7],
                 covX_m1x     = covar[xi,  m1i],
                 covX_m2x     = covar[xi,  m2i],
                 covX_m1m2x   = covar[xi,  7],
                 covm1x_m2x   = covar[m1i, m2i],
                 covm1x_m1m2x = covar[m1i, 7],
                 covm2x_m1m2x = covar[m2i, 7])
    row.names(th) <- row.names(se) <- vals[[vi[1]]]
    colnames(th)  <- colnames(se)  <- vals[[vi[2]]]
    list(theta = th, se = se)
  })

  thetaMat <- setNames(lapply(mats, `[[`, "theta"), name)
  seMat    <- setNames(lapply(mats, `[[`, "se"),    name)
  ZMat     <- setNames(lapply(1:3, \(i) thetaMat[[i]] / seMat[[i]]), name)
  pMat     <- setNames(lapply(ZMat, \(z) 2 * pmin(pnorm(z), 1 - pnorm(z))), name)
  sigMat   <- setNames(lapply(pMat, \(p) p < alpha), name)

  small_vals <- all(sapply(vals, length) < 8)

  figures <- lapply(1:3, function(i) {
    dat2 <- mat_to_long(thetaMat[[i]], "Parameter Value")
    sig2 <- mat_to_long(sigMat[[i]],  "value")

    sig2$pattern <- ifelse(sig2$value == crosshatch_non_sig,
                           "none", "crosshatch")
    sig3 <- sig2[!sig2$value, ]

    p <- ggplot2::ggplot(dat2, ggplot2::aes(Var1, Var2)) +
      ggplot2::geom_tile(ggplot2::aes(fill = `Parameter Value`)) +
      ggplot2::scale_color_identity() +
      ggplot2::scale_fill_gradient2(low = color_low, high = color_high,
                                    mid = color_mid, midpoint = 0) +
      ggpattern::geom_tile_pattern(data = sig3, ggplot2::aes(pattern = pattern),
                                   pattern_density = grid_density,
                                   pattern_spacing = grid_spacing,
                                   pattern_color   = color_grid, alpha = 0) +
      ggplot2::ggtitle(paste0(name[i], grp_label)) +
      ggplot2::theme_bw() +
      ggplot2::guides(pattern = "none") +
      ggplot2::xlab(name[-i][1]) +
      ggplot2::ylab(name[-i][2])

    if (small_vals)
      p <- p + ggplot2::geom_text(
        ggplot2::aes(label = round(`Parameter Value`, round_res)),
        color = color_values)

    if (save) {
      if (!dir.exists(folder)) dir.create(folder, recursive = TRUE)
      ggplot2::ggsave(
        file.path(folder, paste0(name[i], grp_label, ".png")),
        plot = p, dpi = 600, width = 10, height = 10)
    }
    p
  })
  names(figures) <- name

  list(thetas = thetaMat,
       standard_errors = seMat,
       p_values = pMat,
       significance = sigMat,
       plots = figures)
}


# =============================================================================
# Internal support for JNK_bayes
# =============================================================================

# ---- used by: JNK_bayes (via call_support, two-way) ----
#' @import ggplot2
#' @import ggpattern
#' @import scales
#' @import tidyr
#' @import tibble
#' @import dplyr
jnb_support2 <- function(theta,
                         group,
                         theta_1,
                         theta_2,
                         theta_1n,
                         theta_2n,
                         theta_int_12,
                         theta_1_vals,
                         theta_2_vals,
                         color_mid = color_mid,
                         color_low = color_low,
                         color_high = color_high,
                         save, folder) {


  make_td <- function(param_vec, mod_vals, theta_name, mod_name) {
    d <- data.frame(
      modValue  = as.factor(rep(round(mod_vals, 3), each = nrow(theta))),
      parameter = param_vec
    )
    if (!is.null(group)) d$group <-  group
    s <- d |>
      dplyr::summarize(
        thetaPostMean = mean(parameter),
        thetaPostSD   = sd(parameter),
        bayes_p       = sum(parameter > 0) / dplyr::n(),
        thetaPost2.5  = quantile(parameter, 0.025),
        thetaPost97.5 = quantile(parameter, 0.975),
        .by = modValue
      )
    cbind(data.frame(theta     = theta_name,
                     moderator = mod_name), s)
  }


  make_plot <- function(plot_data, x_name, mod_name,
                        color_mid = color_mid,
                        color_low = color_low,
                        color_high = color_high) {
    leg_label <- if (nchar(mod_name) < 26) {mod_name} else {"Moderator"}

    plot_data$modValue <- as.numeric(as.character(plot_data$modValue))

    ggplot2::ggplot(plot_data,
                    ggplot2::aes(x      = parameter,
                                 group  = modValue,
                                 fill   = modValue,
                                 color  = modValue)) +
      ggplot2::geom_vline(xintercept = 0) +
      ggplot2::geom_density(alpha = 0.1) +
      ggplot2::scale_fill_gradient2(low  = color_low,  mid = color_mid,
                                    high = color_high, midpoint = 0) +
      ggplot2::scale_color_gradient2(low  = color_low,  mid = color_mid,
                                     high = color_high, midpoint = 0) +
      ggplot2::theme_bw() +
      ggplot2::labs(
        title = if (is.null(group)) {
            paste0("Posterior density for ", x_name,
                 "\n moderated by ", mod_name)
          } else {
            paste0("Posterior density for ", x_name,
                   "\n moderated by ", mod_name, " (", group, ")")
            },
        x     = x_name,
        y     = "Posterior Density",
        color = leg_label,
        fill  = leg_label
      )
  }


  theta_1x <- as.vector(outer(theta[, 3], theta_2_vals, `*`)) + theta[, 1]
  theta_2x <- as.vector(outer(theta[, 3], theta_1_vals, `*`)) + theta[, 2]

  t1plotData <- data.frame(modValue  = as.factor(rep(round(theta_2_vals, 3),
                                                     each = nrow(theta))),
                           parameter = theta_1x)
  t2plotData <- data.frame(modValue  = as.factor(rep(round(theta_1_vals, 3),
                                                     each = nrow(theta))),
                           parameter = theta_2x)
  if (!is.null(group)) {
    t1plotData$group <- group
    t2plotData$group <- group
  }


  t1d <- make_td(theta_1x, theta_2_vals, theta_1n, theta_2n)
  t2d <- make_td(theta_2x, theta_1_vals, theta_2n, theta_1n)

  g1  <- make_plot(t1plotData, theta_1n, theta_2n,
                   color_mid = color_mid,
                   color_low = color_low,
                   color_high = color_high)
  g2  <- make_plot(t2plotData, theta_2n, theta_1n,
                   color_mid = color_mid,
                   color_low = color_low,
                   color_high = color_high)

  if (save) {
    make_path <- function(tn, mn) {
      if (is.null(group)) {
        file.path(folder, paste0("Posterior density for ", tn,
                                 " moderated by ", mn,".png"))
      } else {
        file.path(folder, paste0("Posterior density for ", tn,
                                 " moderated by ", mn, " (", group, ").png"))
      }
    }
    ggplot2::ggsave(make_path(theta_1n, theta_2n), plot = g1, dpi = 600,
                    width = 10, height = 10)
    ggplot2::ggsave(make_path(theta_2n, theta_1n), plot = g2, dpi = 600,
                    width = 10, height = 10)
  }

  setNames(list(t1d, t2d, g1, g2),
           c(paste(theta_1n, "_table"),
             paste(theta_2n, "_table"),
             paste(theta_1n, "_plot"),
             paste(theta_2n, "_plot")))
}


# ---- used by: JNK_bayes (via call_support, three-way) ----

jnb_support3 <- function(theta,
                         group,
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
                         folder = folder) {


  make_heatmap <- function(data,
                           fill_var,
                           midpoint,
                           title,
                           xlab,
                           ylab,
                           fill_label,
                           pat,
                           color_low,
                           color_mid,
                           color_high,
                           grid_density,
                           grid_spacing,
                           color_grid) {
    ggplot2::ggplot(data, ggplot2::aes(mod1Val, mod2Val)) +
      ggplot2::geom_tile(ggplot2::aes(fill = .data[[fill_var]])) +
      ggplot2::scale_color_identity() +
      ggplot2::scale_fill_gradient2(low = color_low,
                                    high = color_high,
                                    mid = color_mid,
                                    midpoint = midpoint) +
      ggpattern::geom_tile_pattern(data = pat,
                                 ggplot2::aes(pattern = pattern),
                                 pattern_density = grid_density,
                                 pattern_spacing = grid_spacing,
                                 pattern_color   = color_grid,
                                 alpha = 0) +
      ggplot2::theme_bw() +
      ggplot2::guides(pattern = "none") +
      ggplot2::ggtitle(title) +
      ggplot2::xlab(xlab) + ggplot2::ylab(ylab) +
      ggplot2::labs(fill = fill_label)
  }

  make_plot_data <- function(theta, vals1, vals2,
                             theta_name, mod1_name, mod2_name,
                             theta_cols,
                             thresholds) {

    grid <- expand.grid(mod1Val = vals1, mod2Val = vals2)

    coef_mat  <- cbind(1, grid$mod1Val, grid$mod2Val,
                       grid$mod1Val * grid$mod2Val)
    theta_sub <- theta[, theta_cols]
    param_mat <- theta_sub %*% t(coef_mat)

    plotData <- data.frame(
      mod1Val   = rep(grid$mod1Val, each = nrow(theta)),
      mod2Val   = rep(grid$mod2Val, each = nrow(theta)),
      parameter = as.vector(param_mat)
    )

    d <- plotData |>
      dplyr::summarize(
        thetaPostMean = mean(parameter),
        thetaPostSD   = sd(parameter),
        bayes_p       = sum(parameter > 0) / dplyr::n(),
        thetaPost2.5  = quantile(parameter, 0.025),
        thetaPost97.5 = quantile(parameter, 0.975),
        .by = c(mod1Val, mod2Val)
      )

    d$theta   <- theta_name
    d$mod1    <- mod1_name
    d$mod2    <- mod2_name
    d$sig     <- d$bayes_p <= min(thresholds) | d$bayes_p >= max(thresholds)
    d$pattern <- ifelse(!d$sig, "crosshatch", "none")

    d[, c("theta", "mod1", "mod2", "mod1Val", "mod2Val",
          "thetaPostMean", "thetaPostSD", "bayes_p",
          "thetaPost2.5", "thetaPost97.5", "sig", "pattern")]
  }
  nt <- nrow(theta)

  t1d <- make_plot_data(theta, theta_2_vals, theta_3_vals,
                        theta_1n, theta_2n, theta_3n,
                        theta_cols = c(1, 4, 5, 7), thresholds)

  t2d <- make_plot_data(theta, theta_1_vals, theta_3_vals,
                        theta_2n, theta_1n, theta_3n,
                        theta_cols = c(2, 4, 6, 7), thresholds)

  t3d <- make_plot_data(theta, theta_1_vals, theta_2_vals,
                        theta_3n, theta_1n, theta_2n,
                        theta_cols = c(3, 5, 6, 7), thresholds)

  ns               <- c(theta_1n, theta_2n, theta_3n)
  tables           <- list(t1d, t2d, t3d)
  names(tables)    <- ns
  plotsMean        <- vector(mode = "list", length = 3)
  names(plotsMean) <- ns
  plotsP           <- plotsMean


  for (i in 1:3) {
    pat      <- tables[[i]][!tables[[i]]$sig, ]
    mods_out <- ns[-i]

    plotsMean[[i]] <- make_heatmap(
      data       = tables[[i]],
      fill_var   = "thetaPostMean",
      midpoint   = 0,
      title      = if (is.null(group)) {
          paste0("Average posterior density for ",ns[i])
        } else {
          paste0("Average posterior density for ", ns[i], " (", group, ")")
          },
      xlab       = mods_out[[1]],
      ylab       = mods_out[[2]],
      fill_label = "Posterior Mean",
      pat = pat,
      color_low = color_low,
      color_mid = color_mid,
      color_high = color_high,
      grid_density = grid_density,
      grid_spacing = grid_spacing,
      color_grid = color_grid
    )

    plotsP[[i]] <- make_heatmap(
      data       = tables[[i]],
      fill_var   = "bayes_p",
      midpoint   = 0.5,
      title      = if (is.null(group)) {
          paste0("Bayesian p-value for ", ns[i])
        } else {
          paste0("Bayesian p-value for ", ns[i], " (", group, ")")
        },
      xlab       = mods_out[[1]],
      ylab       = mods_out[[2]],
      fill_label = expression("Bayesian" ~ italic("p") * "-value"),
      pat = pat,
      color_low = color_low,
      color_mid = color_mid,
      color_high = color_high,
      grid_density = grid_density,
      grid_spacing = grid_spacing,
      color_grid = color_grid
    )

    if (save) {
      if (is.null(group)) {
        mod_str  <- paste0(mods_out[[1]], " and ", mods_out[[2]])
      } else {
        mod_str  <- paste0(mods_out[[1]], " and ", mods_out[[2]],
                           " (", group, ")")
      }

      base_path <- file.path(folder,
                             paste0("%s for ",
                                    ns[i],
                                    " moderated by ",
                                    mod_str, ".png"))
      ggplot2::ggsave(sprintf(base_path, "Posterior density"),
                      plot = plotsMean[[i]], dpi = 600, width = 10, height = 10)
      ggplot2::ggsave(sprintf(base_path, "Bayesian p-value"),
                      plot = plotsP[[i]],   dpi = 600, width = 10, height = 10)
    }
  }

  returnList <- list(tables,plotsMean,plotsP)
  names(returnList) <- c('result_tables',
                         'post_mean_plots',
                         'bayes_p_plots')
  return(returnList)
}
