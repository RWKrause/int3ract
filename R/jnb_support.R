#' Support function for 2-way JN plots for Bayesian SAOMs
#'
#' Computes conditional posterior distributions and produces Johnson-Neyman
#' density plots for a two-way interaction between two SAOM parameters.
#'
#' @param theta matrix; posterior draws (rows = iterations) with three columns:
#'   \[,1\] main effect of theta_1, \[,2\] main effect of theta_2,
#'   \[,3\] interaction between theta_1 and theta_2. Extracted from either
#'   \code{ThinParameters} (fixed/Eta) or \code{ThinPosteriorMu} (random/Mu).
#' @param group character; label identifying the parameter source, either
#'   \code{'Eta'} for fixed effects, \code{'Mu'} for the hyper-mean of random
#'   effects, or \code{'group<i>'} for group-specific random effects.

#' @param theta_1 numeric; index of the first parameter in the (rate-excluded)
#'   effects object. See \code{modelOut$effects[modelOut$effects$type != 'rate',]}.
#' @param theta_2 numeric; index of the second parameter in the (rate-excluded)
#'   effects object. See \code{modelOut$effects[modelOut$effects$type != 'rate',]}.
#' @param theta_1n character; display name for theta_1, derived from
#'   \code{effectName} with special characters stripped.
#' @param theta_2n character; display name for theta_2, derived from
#'   \code{effectName} with special characters stripped.
#' @param theta_int_12 numeric; index of the interaction between theta_1 and theta_2
#'   in the (rate-excluded) effects object.
#' @param theta_1_vals numeric vector; values of the theta_1 statistic at which
#'   to evaluate the conditional effect of theta_2 (i.e., the moderator range
#'   for the second plot).
#' @param theta_2_vals numeric vector; values of the theta_2 statistic at which
#'   to evaluate the conditional effect of theta_1 (i.e., the moderator range
#'   for the first plot).
#'
#' @returns A named list with four elements:
#'   \describe{
#'     \item{\code{<theta_1n>_table}}{data.frame of posterior summaries
#'       (mean, SD, Bayesian p, 2.5th and 97.5th percentiles) for the
#'       conditional effect of theta_1 across \code{theta_2_vals}.}
#'     \item{\code{<theta_2n>_table}}{As above for the conditional effect
#'       of theta_2 across \code{theta_1_vals}.}
#'     \item{\code{<theta_1n>_ggplot}}{ggplot2 density plot for the
#'       conditional effect of theta_1.}
#'     \item{\code{<theta_2n>_ggplot}}{ggplot2 density plot for the
#'       conditional effect of theta_2.}
#'   }
#'
#' @import ggplot2
#' @import ggpattern
#' @import scales
#' @import tidyr
#' @import tibble
#' @import wesanderson
#'
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
      parameter = param_vec,
      group     = group
    )
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
  
  
  make_plot <- function(plot_data, x_name, mod_name) {
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
        title = paste0("Posterior density for ", x_name,
                       "\n moderated by ", mod_name, " (", group, ")"),
        x     = x_name,
        y     = "Posterior Density",
        color = leg_label,
        fill  = leg_label
      )
  }
  
  make_path <- function(tn, mn) {
    file.path(folder, paste0("Posterior density for ", tn,
                             " moderated by ", mn, " (", group, ").png"))
  }
  theta_1x <- as.vector(outer(theta[, 3], theta_2_vals, `*`)) + theta[, 1]
  theta_2x <- as.vector(outer(theta[, 3], theta_1_vals, `*`)) + theta[, 2]
  
  t1plotData <- data.frame(modValue  = as.factor(rep(round(theta_2_vals, 3),
                                                     each = nrow(theta))),
                           parameter = theta_1x, 
                           group = group)
  t2plotData <- data.frame(modValue  = as.factor(rep(round(theta_1_vals, 3),
                                                     each = nrow(theta))),
                           parameter = theta_2x, group = group)
  
  t1d <- make_td(theta_1x, theta_2_vals, theta_1n, theta_2n)
  t2d <- make_td(theta_2x, theta_1_vals, theta_2n, theta_1n)
  
  g1  <- make_plot(t1plotData, theta_1n, theta_2n)
  g2  <- make_plot(t2plotData, theta_2n, theta_1n)
  
  if (save) {
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

#' Support function for 3-way JN plots for Bayesian SAOMs
#'
#' Computes conditional posterior distributions and produces Johnson-Neyman
#' heatmap plots for a three-way interaction among three SAOM parameters.
#' For each parameter, its conditional effect is evaluated across a grid of
#' values for the other two, yielding both posterior mean and Bayesian p-value
#' heatmaps. Cells outside the significance threshold are overlaid with a
#' crosshatch pattern.
#'
#' @param theta matrix; posterior draws (rows = iterations) with seven columns
#'   in the order: \[,1\] theta_1, \[,2\] theta_2, \[,3\] theta_3,
#'   \[,4\] Int12, \[,5\] Int13, \[,6\] Int23, \[,7\] Int123.
#'   Extracted from either \code{ThinParameters} or \code{ThinPosteriorMu}.
#' @param group character; label identifying the parameter source, either
#'   \code{'Eta'} for fixed effects, \code{'Mu'} for the hyper-mean of random
#'   effects, or \code{'group<i>'} for group-specific random effects.
#' @param theta_1 numeric; index of the first parameter in the (rate-excluded)
#'   effects object. See \code{modelOut$effects[modelOut$effects$type != 'rate',]}.
#' @param theta_2 numeric; index of the second parameter in the (rate-excluded)
#'   effects object.
#' @param theta_3 numeric; index of the third parameter in the (rate-excluded)
#'   effects object.
#' @param theta_1n character; display name for theta_1, derived from
#'   \code{effectName} with special characters stripped.
#' @param theta_2n character; display name for theta_2, derived from
#'   \code{effectName} with special characters stripped.
#' @param theta_3n character; display name for theta_3, derived from
#'   \code{effectName} with special characters stripped.
#' @param theta_int_12 numeric; index of the theta_1 × theta_2 interaction term.
#' @param theta_int_13 numeric; index of the theta_1 × theta_3 interaction term.
#' @param theta_int_23 numeric; index of the theta_2 × theta_3 interaction term.
#' @param theta_int_123 numeric; index of the three-way theta_1 × theta_2 × theta_3
#'   interaction term.
#' @param theta_1_vals numeric vector; grid of theta_1 statistic values over which
#'   to evaluate conditional effects.
#' @param theta_2_vals numeric vector; grid of theta_2 statistic values over which
#'   to evaluate conditional effects.
#' @param theta_3_vals numeric vector; grid of theta_3 statistic values over which
#'   to evaluate conditional effects.
#' @param thresholds numeric vector of length 2; lower and upper Bayesian
#'   p-value bounds used to determine significance. Cells outside this range
#'   are considered significant and are not crosshatched. E.g. \code{c(0.05, 0.95)}.
#' @param color_mid character; color for the midpoint of the fill scale
#'   (posterior mean plots: midpoint = 0; p-value plots: midpoint = 0.5).
#'   Default \code{'white'}.
#' @param color_low character; color for the low end of the fill scale.
#'   Default \code{'#F05039'} (red).
#' @param color_high character; color for the high end of the fill scale.
#'   Default \code{'#000066'} (dark blue).
#' @param color_values character; color for displayed numeric values.
#'   Default \code{'grey40'}.
#' @param color_grid character; color of the crosshatch pattern overlaid on
#'   non-significant cells. Default \code{'black'}.
#' @param grid_density numeric; density of the crosshatch pattern lines.
#'   Default \code{0.01}.
#' @param grid_spacing numeric; spacing between crosshatch pattern lines.
#'   Default \code{0.1}.
#'
#' @returns A named list with three elements:
#'   \describe{
#'     \item{\code{result_tables}}{Named list of three data.frames (one per
#'       parameter) containing posterior mean, SD, Bayesian p, and 2.5th/97.5th
#'       percentiles across all moderator value combinations.}
#'     \item{\code{post_mean_plots}}{Named list of three ggplot2 heatmaps
#'       showing the posterior mean of each parameter's conditional effect.}
#'     \item{\code{bayes_p_plots}}{Named list of three ggplot2 heatmaps
#'       showing the Bayesian p-value of each parameter's conditional effect.}
#'   }
#'
#' @import ggplot2
#' @import ggpattern
#' @import scales
#' @import tibble
#' @import dplyr
#'
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
      ggplot2::geom_tile_pattern(data = pat, 
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
    
    mod_str  <- paste0(mods_out[[1]], 
                       " and ", 
                       mods_out[[2]], 
                       " (", group, ")")
    base_path <- file.path(folder,
                           paste0("%s for ", 
                                  ns[i], 
                                  " moderated by ", 
                                  mod_str, ".png"))
    
    plotsMean[[i]] <- make_heatmap(
      data       = tables[[i]],
      fill_var   = "thetaPostMean",
      midpoint   = 0,
      title      = paste0("Average posterior density for ", 
                          ns[i], 
                          " (", group, ")"),
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
      title      = paste0("Bayesian p-value for ", 
                          ns[i], 
                          " (", group, ")"),
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