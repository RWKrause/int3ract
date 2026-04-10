#' Function for creating 2- and 3-way plots and tables based on any frequentist multiplicative model.
#'
#' @param name character; names of the variables involved in interaction
#' @param covar matrix; covariance of relevant parameters only
#' @param coefs numeric; vector of coefficient values
#' @param vals list; values of the variables
#' @param alpha see JNK_lm() or JNK_siena()
#' @param round_res see JNK_lm() or JNK_siena()
#' @param control_fdr see JNK_lm() or JNK_siena()
#' @param range_size see JNK_lm() or JNK_siena()
#' @param sig_color see JNK_lm() or JNK_siena()
#' @param non_sig_color see JNK_lm() or JNK_siena()
#' @param line_color see JNK_lm() or JNK_siena()
#' @param color_mid see JNK_lm() or JNK_siena()
#' @param color_low see JNK_lm() or JNK_siena()
#' @param color_high see JNK_lm() or JNK_siena()
#' @param color_values see JNK_lm() or JNK_siena()
#' @param color_grid see JNK_lm() or JNK_siena()
#' @param grid_density see JNK_lm() or JNK_siena()
#' @param grid_spacing see JNK_lm() or JNK_siena()
#' @param crosshatch_non_sig see JNK_lm() or JNK_siena()
#'
#' @returns Returns tables and plots for interactions
#' 
#' @import ggplot2
#' @import ggpattern
#' @import scales
#' @import tidyr
#' @import tibble
#' 
#' @export
JNK_int <- function(name,
                     covar,
                     coefs,
                     vals,
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
                     crosshatch_non_sig = TRUE) {
  
  k      <- length(name)
  z_crit <- abs(qnorm(alpha / 2)) 
  
  for (var in 1:k)
    vals[[var]] <- seq(min(vals[[var]]), 
                       max(vals[[var]]), 
                       length.out = range_size)
  
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
        ggplot2::ggtitle(paste0("Moderation analysis for ", name[i]))
      
      add_ribbons(p, td[[i]]) +
        ggplot2::geom_path(linewidth = 1, color = line_color)
    })
    names(plots) <- name
    
    return(list(param_table = td, plots = plots))
  }
  
  x_idx  <- c(1L, 2L, 3L)
  m1_idx <- c(4L, 4L, 5L)
  m2_idx <- c(5L, 6L, 6L)
  v_idx  <- list(c(2L,3L), c(1L,3L), c(1L,2L))
  
  mat_to_long <- function(mat, val_name) {
    mat |>
      as.data.frame() |>
      tibble::rownames_to_column("Var1") |>
      tidyr::pivot_longer(-Var1, names_to = "Var2", values_to = val_name) |>
      dplyr::mutate(Var1 = as.numeric(Var1), Var2 = as.numeric(Var2))
  }
  
  mats <- lapply(1:3, function(i) {
    xi  <- x_idx[i]
    m1i <- m1_idx[i]
    m2i <- m2_idx[i]
    vi <- v_idx[[i]]
    th  <- outer(vals[[vi[1]]], 
                 vals[[vi[2]]], 
                 parFun,
                 bx = coefs[xi], 
                 bm1x = coefs[m1i],
                 bm2x = coefs[m2i],
                 bm1m2x = coefs[7])
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
  pMat     <- setNames(lapply(ZMat, \(z) 2 * pmin(pnorm(z), 1 - pnorm(z))),name)
  sigMat   <- setNames(lapply(pMat, \(p) p < alpha),name)
  
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
      ggplot2::ggtitle(name[i]) +
      ggplot2::theme_bw() +
      ggplot2::guides(pattern = "none") +
      ggplot2::xlab(name[-i][1]) +
      ggplot2::ylab(name[-i][2])
    
    if (small_vals)
      p <- p + ggplot2::geom_text(
        ggplot2::aes(label = round(`Parameter Value`, round_res)),
        color = color_values)
    p
  })
  names(figures) <- name
  
  list(thetas = thetaMat, 
       standard_errors = seMat,
       p_values = pMat, 
       significance = sigMat, 
       plots = figures)
}