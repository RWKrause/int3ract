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
#' @param sbo multiSiena; the \code{multiSiena()} output object. Used to access
#'   model metadata such as the number of groups.
#' @param theta_1 numeric; index of the first parameter in the (rate-excluded)
#'   effects object. See \code{sbo$effects[sbo$effects$type != 'rate',]}.
#' @param theta_2 numeric; index of the second parameter in the (rate-excluded)
#'   effects object. See \code{sbo$effects[sbo$effects$type != 'rate',]}.
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
#'
jnb_support2 <- function(theta,
                         group,
                         sbo = sbo,
                         theta_1 = theta_1, 
                         theta_2 = theta_2, 
                         theta_1n = theta_1n, 
                         theta_2n = theta_2n, 
                         theta_int_12 = theta_int_12, 
                         theta_1_vals = theta_1_vals, 
                         theta_2_vals = theta_2_vals, 
                         save = save,
                         folder = folder) { 
  
  
  theta_1x <- as.vector(outer(theta[,3] , theta_2_vals, FUN = '*')) + theta[,1]
  theta_2x <- as.vector(outer(theta[,3] , theta_1_vals, FUN = '*')) + theta[,2]
  
  t1plotData <- data.frame(modValue = as.factor(rep(round(theta_2_vals,3),
                                                    each = nrow(theta))),
                           parameter = theta_1x,
                           group = group)
  
  t2plotData <- data.frame(modValue = as.factor(rep(round(theta_1_vals,3),
                                                    each = nrow(theta))),
                           parameter = theta_2x,
                           group = group)
  
  t1d <- t1plotData |>
    summarize(thetaPostMean = mean(parameter),
              thetaPostSD = sd(parameter),
              bayes_p = sum(parameter > 0) / length(parameter),
              thetaPost2.5    = quantile(parameter,0.025),
              thetaPost97.5   = quantile(parameter,0.975),
              .by = modValue
    )
  
  
  t2d <- t2plotData |>
    summarize(thetaPostMean = mean(parameter),
              thetaPostSD = sd(parameter),
              bayes_p = sum(parameter > 0) / length(parameter),
              thetaPost2.5    = quantile(parameter,0.025),
              thetaPost97.5   = quantile(parameter,0.975),
              .by = modValue
    )
  
  t1d <- cbind(data.frame(theta     = rep(theta_1n,length(theta_2_vals)),
                          moderator = rep(theta_2n,length(theta_2_vals))),
               t1d)
  
  t2d <- cbind(data.frame(theta     = rep(theta_2n,length(theta_1_vals)),
                          moderator = rep(theta_1n,length(theta_1_vals))),
               t2d)
  
  
  g1 <- ggplot(t1plotData, aes(x = parameter)) +
    geom_vline(xintercept = 0) +
    geom_density(alpha = 0.1, aes(fill = modValue, color = modValue)) +
    theme_bw()
  
  if (nchar(theta_2n) < 26) {
    g1 <- g1 + 
      labs(title = paste0('Posterior density for ',
                          theta_1n,
                          '\n moderated by ',
                          theta_2n,
                          ' (',
                          group,
                          ')'),
           x = theta_1n,
           y = 'Posterior Density',
           color = theta_2n,
           fill = theta_2n)
  } else {
    g1 <- g1 + 
      labs(title = paste0('Posterior density for ',
                          theta_1n,
                          '\n moderated by ',
                          theta_2n,
                          ' (',
                          group,
                          ')'),
           x = theta_1n,
           y = 'Posterior Density',
           color = 'Moderator',
           fill = 'Moderator')
  }
  
  
  g2 <- ggplot(t2plotData, aes(x = parameter)) +
    geom_vline(xintercept = 0) +
    geom_density(alpha = 0.1, aes(fill = modValue, color = modValue))  +
    theme_bw()
  
  if (nchar(theta_1n) < 26) {
    g2 <- g2 + 
      labs(title = paste0('Posterior density for ',
                          theta_2n,
                          '\n moderated by ',
                          theta_1n,
                          ' (',
                          group,
                          ')'),
           x = theta_2n,
           y = 'Posterior Density',
           color = theta_1n,
           fill = theta_1n)
  } else {
    g2 <- g2 + 
      labs(title = paste0('Posterior density for ',
                          theta_2n,
                          '\n moderated by ',
                          theta_1n,
                          ' (',
                          group,
                          ')'),
           x = theta_2n,
           y = 'Posterior Density',
           color = 'Moderator',
           fill = 'Moderator') 
  }
  
  if (grepl('group',group)) {
    g1n <- paste0(folder,
                  '/Posterior density for ',
                  theta_1n,
                  ' moderated by ',
                  theta_2n,
                  ' (',
                  group,
                  ')',
                  '.png')
    
    g2n <- paste0(folder,
                  '/Posterior density for ',
                  theta_2n,
                  ' moderated by ',
                  theta_1n,
                  ' (',
                  group,
                  ')',
                  '.png')
  } else {
    g1n <- paste0(folder,
                  '/Posterior density for ',
                  theta_1n,
                  ' moderated by ',
                  theta_2n,
                  ' (',
                  group,
                  ')',
                  '.png')
    
    g2n <- paste0(folder,
                  '/Posterior density for ',
                  theta_2n,
                  ' moderated by ',
                  theta_1n,
                  ' (',
                  group,
                  ')',
                  '.png')
    
    plot(g1)
    plot(g2)
  }
  
  if (save) {
    ggsave(g1n, plot = g1, dpi = 600, width = 10, height = 10)
    ggsave(g2n, plot = g2, dpi = 600, width = 10, height = 10)
  }
 
  
  returnList <- list(t1d,t2d,g1,g2)
  names(returnList) <- c(paste(theta_1n, '_table'),
                         paste(theta_2n, '_table'),
                         paste(theta_1n, '_ggplot'),
                         paste(theta_2n, '_ggplot'))
  return(returnList)
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
#' @param sbo multiSiena; the \code{multiSiena()} output object.
#' @param theta_1 numeric; index of the first parameter in the (rate-excluded)
#'   effects object. See \code{sbo$effects[sbo$effects$type != 'rate',]}.
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
                         sbo = sbo,
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
                         grid_spacing = grid_spacing) { 
  
  nt <- nrow(theta)
  n1 <- nt * length(theta_2_vals) * length(theta_3_vals)
  
  t1plotData <- data.frame(mod1      = rep(NA,n1),
                           mod1Val   = rep(NA,n1),
                           mod2      = rep(NA,n1),
                           mod2Val   = rep(NA,n1),
                           parameter = rep(NA,n1))  
  
  count <- 0
  for (i in 1:length(theta_2_vals)) {
    for (j in 1:length(theta_3_vals)) {
      count <- count + 1
      thetaX <- theta[,1] + theta_2_vals[i] * theta[,4] + 
        theta_3_vals[j] * theta[,5] +
        theta_2_vals[i] * theta_3_vals[j] * theta[,7]
      t1plotData$mod1[(-nt + count * nt + 1):(count * nt)]      <- theta_2n
      t1plotData$mod1Val[(-nt + count * nt + 1):(count * nt)]   <- theta_2_vals[i]
      t1plotData$mod2[(-nt + count * nt + 1):(count * nt)]      <- theta_3n
      t1plotData$mod2Val[(-nt + count * nt + 1):(count * nt)]   <- theta_3_vals[j]
      t1plotData$parameter[(-nt + count * nt + 1):(count * nt)] <- thetaX
    }
  }
  
  n2 <- nt * length(theta_1_vals) * length(theta_3_vals)
  
  t2plotData <- data.frame(mod1      = rep(NA,n2),
                           mod1Val   = rep(NA,n2),
                           mod2      = rep(NA,n2),
                           mod2Val   = rep(NA,n2),
                           parameter = rep(NA,n2))  
  
  count <- 0
  for (i in 1:length(theta_1_vals)) {
    for (j in 1:length(theta_3_vals)) {
      count <- count + 1
      thetaX <- theta[,2] + theta_1_vals[i] * theta[,4] + 
        theta_3_vals[j] * theta[,6] +
        theta_1_vals[i] * theta_3_vals[j] * theta[,7]
      t2plotData$mod1[(-nt + count * nt + 1):(count * nt)]      <- theta_1n
      t2plotData$mod1Val[(-nt + count * nt + 1):(count * nt)]   <- theta_1_vals[i]
      t2plotData$mod2[(-nt + count * nt + 1):(count * nt)]      <- theta_3n
      t2plotData$mod2Val[(-nt + count * nt + 1):(count * nt)]   <- theta_3_vals[j]
      t2plotData$parameter[(-nt + count * nt + 1):(count * nt)] <- thetaX
    }
  }
  
  n3 <- nt * length(theta_1_vals) * length(theta_2_vals)
  
  t3plotData <- data.frame(mod1      = rep(NA,n3),
                           mod1Val   = rep(NA,n3),
                           mod2      = rep(NA,n3),
                           mod2Val   = rep(NA,n3),
                           parameter = rep(NA,n3))  
  
  count <- 0
  for (i in 1:length(theta_1_vals)) {
    for (j in 1:length(theta_2_vals)) {
      count <- count + 1
      thetaX <- theta[,3] + theta_1_vals[i] * theta[,5] + 
        theta_2_vals[j] * theta[,6] +
        theta_1_vals[i] * theta_2_vals[j] * theta[,7]
      t3plotData$mod1[(-nt + count * nt + 1):(count * nt)]      <- theta_1n
      t3plotData$mod1Val[(-nt + count * nt + 1):(count * nt)]   <- theta_1_vals[i]
      t3plotData$mod2[(-nt + count * nt + 1):(count * nt)]      <- theta_2n
      t3plotData$mod2Val[(-nt + count * nt + 1):(count * nt)]   <- theta_2_vals[j]
      t3plotData$parameter[(-nt + count * nt + 1):(count * nt)] <- thetaX
    }
  }
  
  
  t1d <- t1plotData |>
    dplyr::summarize(thetaPostMean = mean(parameter),
              thetaPostSD = sd(parameter),
              bayes_p = sum(parameter > 0) / length(parameter),
              thetaPost2.5    = quantile(parameter,0.025),
              thetaPost97.5   = quantile(parameter,0.975),
              .by = c(mod1Val,mod2Val)
    )
  
  
  t2d <- t2plotData |>
    dplyr::summarize(thetaPostMean = mean(parameter),
              thetaPostSD = sd(parameter),
              bayes_p = sum(parameter > 0) / length(parameter),
              thetaPost2.5    = quantile(parameter,0.025),
              thetaPost97.5   = quantile(parameter,0.975),
              .by = c(mod1Val,mod2Val)
    )  
  
  t3d <- t3plotData |>
    dplyr::summarize(thetaPostMean = mean(parameter),
              thetaPostSD = sd(parameter),
              bayes_p = sum(parameter > 0) / length(parameter),
              thetaPost2.5    = quantile(parameter,0.025),
              thetaPost97.5   = quantile(parameter,0.975),
              .by = c(mod1Val,mod2Val)
    )
  
  
  t1d <- cbind(data.frame(theta = rep(theta_1n,nrow(t1d)),
                          mod1  = rep(theta_2n,nrow(t1d)),
                          mod2  = rep(theta_3n,nrow(t1d))),
               t1d)
  
  t2d <- cbind(data.frame(theta = rep(theta_2n,nrow(t2d)),
                          mod1  = rep(theta_1n,nrow(t2d)),
                          mod2  = rep(theta_3n,nrow(t2d))),
               t2d)
  
  t3d <- cbind(data.frame(theta = rep(theta_3n,nrow(t3d)),
                          mod1  = rep(theta_1n,nrow(t3d)),
                          mod2  = rep(theta_2n,nrow(t3d))),
               t3d)
  
  
  t1d$sig <- t1d$bayes_p <= min(thresholds) |
    t1d$bayes_p >= max(thresholds)
  
  t2d$sig <- t2d$bayes_p <= min(thresholds) |
    t2d$bayes_p >= max(thresholds)
  
  t3d$sig <- t3d$bayes_p <= min(thresholds) |
    t3d$bayes_p >= max(thresholds)
  
  t1d$pattern <- ifelse(!t1d$sig, "crosshatch",'none')
  t2d$pattern <- ifelse(!t2d$sig, "crosshatch",'none')
  t3d$pattern <- ifelse(!t3d$sig, "crosshatch",'none')
  
  ns <- c(theta_1n,theta_2n,theta_3n)
  plotsMean <- vector(mode = 'list', length = 3)
  tables <- list(t1d,t2d,t3d)
  
  names(plotsMean) <- names(tables) <- ns
  
  plotsP <- plotsMean
  
  for (i in 1:3) {
    pat <- tables[[i]][!tables[[i]]$sig,]
    
    plotsMean[[i]] <- ggplot(tables[[i]], aes(mod1Val, mod2Val)) +
      geom_tile(aes(fill = thetaPostMean )) +
      scale_color_identity() +
      scale_fill_gradient2(low = color_low, 
                           high = color_high,
                           mid = color_mid,
                           midpoint = 0) +
      geom_tile_pattern(data = pat, aes(pattern = pattern),
                        pattern_density = grid_density,
                        pattern_spacing = grid_spacing,
                        pattern_color = color_grid,
                        alpha = 0) +
      theme_bw() +
      guides(pattern = "none") + 
      ggtitle(paste0('Average posterior density for ',
                     ns[i],' (',group,')')) +
      xlab(ns[-i][[1]]) +
      ylab(ns[-i][[2]]) +
      labs(fill = "Posterior Mean")
    
    
    plotsP[[i]] <- ggplot(tables[[i]], aes(mod1Val, mod2Val)) +
      geom_tile(aes(fill = bayes_p)) +
      scale_color_identity() +
      scale_fill_gradient2(low = color_low, 
                           high = color_high,
                           mid = color_mid,
                           midpoint = 0.5) +
      geom_tile_pattern(data = pat, aes(pattern = pattern),
                        pattern_density = grid_density,
                        pattern_spacing = grid_spacing,
                        pattern_color = color_grid,
                        alpha = 0) +
      theme_bw() +
      guides(pattern = "none") + 
      ggtitle(paste0('Bayesian p-value for ',
                     ns[i],' (',group,')')) +
      xlab(ns[-i][[1]]) +
      ylab(ns[-i][[2]]) +
      labs(fill = expression("Bayesian" ~ italic('p') * "-value"))
    
    if (grepl('group',group)) {
      plotName <- paste0(folder,
                         '/Posterior density for ',
                         ns[i],
                         ' moderated by ',
                         ns[-i][[1]], 
                         ' and ',
                         ns[-i][[2]],
                         ' (',
                         group,
                         ')',
                         '.png')
      plotNameP <- paste0(folder,
                          '/Bayesian p-value for ',
                          ns[i],
                          ' moderated by ',
                          ns[-i][[1]], 
                          ' and ',
                          ns[-i][[2]],
                          ' (',
                          group,
                          ')',
                          '.png')
    } else {
      plotName <- paste0(folder,
                         '/Posterior density for ',
                         ns[i],
                         ' moderated by ',
                         ns[-i][[1]], 
                         ' and ',
                         ns[-i][[2]],
                         ' (',
                         group,
                         ')',
                         '.png')
      plotNameP <- paste0(folder,
                          '/Bayesian p-value for ',
                          ns[i],
                          ' moderated by ',
                          ns[-i][[1]], 
                          ' and ',
                          ns[-i][[2]],
                          ' (',
                          group,
                          ')',
                          '.png')
    }
    
    if (save) {
      ggsave(plotName,  plot = plotsMean[[i]], dpi = 600, 
             width = 10, height = 10)
      ggsave(plotNameP, plot = plotsP[[i]],    dpi = 600, 
             width = 10, height = 10)
    }

  }
  
  
  
  returnList <- list(tables,plotsMean,plotsP)
  names(returnList) <- c('result_tables',
                         'post_mean_plots',
                         'bayes_p_plots')
  return(returnList)
}