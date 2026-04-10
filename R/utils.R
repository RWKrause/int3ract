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
  ggplot2::ggplot(plot_data, ggplot2::aes(x = parameter)) +
    ggplot2::geom_vline(xintercept = 0) +
    ggplot2::geom_density(alpha = 0.1,
                          ggplot2::aes(fill = modValue, color = modValue)) +
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

extract_theta <- function(thetas) {
  eff_ran <- eff[eff$randomEffects, ]
  t(sapply(thetas, function(i) {
    if (eff$randomEffects[i]) {
      sbo$ThinPosteriorMu[idx, which(as.numeric(rownames(eff_ran)) == i)]
    } else {
      sbo$ThinParameters[idx, 1, i]
    }
  })) |> t()   
}