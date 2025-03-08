#' Internal frequentist interaction finction
#'
#' @param name charachter; names of the variables involved in interaction
#' @param covar matrix; covariance of relevant parameters
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
freq_int <- function(name,
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
                     color_mid = 'white',
                     color_low = '#F05039',
                     color_high = '#000066',
                     color_values = 'grey40',
                     color_grid = 'black',
                     grid_density = 0.01,
                     grid_spacing = 0.1,
                     crosshatch_non_sig = TRUE
                     ) {
  
  
  library(ggplot2)
  library(ggpattern)
  library(tidyverse)
  library(scales)
  
  k <- length(name)
  for (var in 1:k) {
    vals[[var]] <- seq(min(vals[[var]], nam = TRUE), 
                       max(vals[[var]], nam = TRUE),
                       length.out = range_size)
  }

  td <- list()
  if (length(vals) == 2) {
    for (var in 1:k) {
      theta1s <- coefs[var] + coefs[3] * vals[[c(1:k)[-var]]]
      seT1 <- sqrt(covar[var,var] + 
                     vals[[c(1:k)[-var]]] * 2 * covar[c(1:k)[-var],var] + 
                     vals[[c(1:k)[-var]]] ^ 2 * covar[3,3])
      z1 <- theta1s/seT1
      p1 <- c()
      for (i in 1:length(vals[[c(1:k)[-var]]])) {
        p1[i] <-  2 * pmin(pnorm(z1[[i]]), (1 - pnorm(z1[[i]])))
      }
      if (control_fdr) {
        p1O <- p1[order(p1)]
        bhT1 <- alpha * c(1:length(vals[[c(1:k)[-var]]])) / 
          length(vals[[c(1:k)[-var]]])
        sig1 <- p1O < bhT1
        
        if (!all(sig1)) {
          first_false <- which.max(!sig1)
          sig1[first_false:length(vals[[c(1:k)[-var]]])] <- FALSE
        }
        
        sig11 <- vector(length = length(vals[[c(1:k)[-var]]]))
        for (i in 1:length(vals[[c(1:k)[-var]]])) {
          sig11[order(p1)[i]] <- sig1[i]
        }
      } else {
        sig11 <- p1 < alpha
      }
      
      td[[var]] <- data.frame(theta      = rep(name[var],
                                               length(vals[[c(1:k)[-var]]])),
                              moderator  = rep(name[var],
                                               length(vals[[c(1:k)[-var]]])),
                              mod_vals   = round(vals[[c(1:k)[-var]]],
                                                 round_res),
                              theta_vals = round(theta1s,round_res),
                              theta_se   = round(seT1,round_res),
                              theta_p    = round(p1,3),
                              sig = sig11)
    }
    
    plots <- vector(mode = 'list', length = 2)
    names(plots) <- names(td) <- name
    
    for (i in 1:2) {
      plots[[i]] <- ggplot(data = td[[i]], 
                           aes(x = mod_vals, 
                               y = theta_vals)) + 
        geom_hline(aes(yintercept = 0)) +
        theme_bw() + 
        ylab(paste0('Theta for ',names(plots)[i])) +
        xlab(names(plots)[-i]) +
        ggtitle(paste0('Moderation analysis for ',names(plots)[i]))
      
      if (!any(td[[i]]$sig)) {
        plots[[i]] <- plots[[i]] + 
          geom_ribbon(aes(ymin = -abs(qnorm(alpha / 2)) * theta_se + 
                            theta_vals,
                          ymax =  abs(qnorm(alpha / 2)) * theta_se + 
                            theta_vals),
                      fill = non_sig_color,
                      alpha = 0.7)
      } else if (all(td[[i]]$sig)) {
        plots[[i]] <- plots[[i]] + 
          geom_ribbon(aes(ymin = -abs(qnorm(alpha / 2)) * theta_se +
                            theta_vals,
                          ymax =  abs(qnorm(alpha / 2)) * theta_se +
                            theta_vals),
                      fill = sig_color,
                      alpha = 0.7)
      } else {
        nvals <- nrow(td[[i]])
        if (td[[i]]$sig[1]) {
          first_false <- which.max(!td[[i]]$sig)
          last_false  <- nrow(td[[i]]) - which.max(rev(!td[[i]]$sig)) + 1
          
          
          plots[[i]] <- plots[[i]] + 
            geom_ribbon(data = td[[i]][1:(first_false - 1),],
                        aes(ymin = -abs(qnorm(alpha / 2)) * theta_se +
                              theta_vals,
                            ymax =  abs(qnorm(alpha / 2)) * theta_se +
                              theta_vals),
                        fill = sig_color,
                        alpha = 0.7) + 
            geom_ribbon(data = td[[i]][first_false:last_false,],
                        aes(ymin = -abs(qnorm(alpha / 2)) * theta_se +
                              theta_vals,
                            ymax =  abs(qnorm(alpha / 2)) * theta_se +
                              theta_vals),
                        fill = non_sig_color,
                        alpha = 0.7)
          if (td[[i]]$sig[nvals]) {
            plots[[i]] <- plots[[i]] + 
              geom_ribbon(data = td[[i]][(last_false + 1):nvals,],
                          aes(ymin = -abs(qnorm(alpha / 2)) * theta_se +
                                theta_vals,
                              ymax =  abs(qnorm(alpha / 2)) * theta_se +
                                theta_vals),
                          fill = sig_color,
                          alpha = 0.7)
          }
        } else {
          first_true <- which.max(td[[i]]$sig)
          last_true  <- nrow(td[[i]]) - which.max(rev(td[[i]]$sig)) + 1
          
          
          plots[[i]] <- plots[[i]] + 
            geom_ribbon(data = td[[i]][1:(first_true - 1),],
                        aes(ymin = -abs(qnorm(alpha / 2)) * theta_se +
                            theta_vals,
                            ymax =  abs(qnorm(alpha / 2)) * theta_se +
                              theta_vals),
                        fill = non_sig_color,
                        alpha = 0.7) + 
            geom_ribbon(data = td[[i]][first_true:last_true,],
                        aes(ymin = -abs(qnorm(alpha / 2)) * theta_se +
                              theta_vals,
                            ymax =  abs(qnorm(alpha / 2)) * theta_se +
                              theta_vals),
                        fill = sig_color,
                        alpha = 0.7)
          if (!td[[i]]$sig[nvals]) {
            plots[[i]] <- plots[[i]] + 
              geom_ribbon(data = td[[i]][(last_true + 1):nvals,],
                          aes(ymin = -abs(qnorm(alpha / 2)) * theta_se +
                                theta_vals,
                              ymax =  abs(qnorm(alpha / 2)) * theta_se +
                                theta_vals),
                          fill = non_sig_color,
                          alpha = 0.7)
          }
        }
      }
      plots[[i]] <- plots[[i]] +  
        geom_path(linewidth = 1, color = line_color) 
    }
    
    
    
    return_list <- list(param_table = td,  plots = plots)
    
  } else {
    thetaMat <- vector(mode = 'list', length = 3)
    seMat    <- vector(mode = 'list', length = 3)
    ZMat     <- vector(mode = 'list', length = 3)
    pMat     <- vector(mode = 'list', length = 3)
    sigMat   <- vector(mode = 'list', length = 3)
    figures  <- vector(mode = 'list', length = 3)
    names(thetaMat) <- names(seMat) <- names(figures) <- name
    names(ZMat) <- names(pMat) <- names(sigMat) <- name
    
    thetaMat[[1]] <- outer(vals[[2]],vals[[3]],
                           parFun,
                           bx     = coefs[1],
                           bm1x   = coefs[4],
                           bm2x   = coefs[5],
                           bm1m2x = coefs[7])
    
    thetaMat[[2]] <- outer(vals[[1]],vals[[3]],
                           parFun,
                           bx     = coefs[2],
                           bm1x   = coefs[4],
                           bm2x   = coefs[6],
                           bm1m2x = coefs[7])
    
    thetaMat[[3]] <- outer(vals[[1]],vals[[3]],
                           parFun,
                           bx     = coefs[3],
                           bm1x   = coefs[5],
                           bm2x   = coefs[6],
                           bm1m2x = coefs[7])
    
    
    seMat[[1]] <-  outer(vals[[2]],vals[[3]],
                         seFun,
                         Vx           = covar[1,1],
                         Vm1x         = covar[4,4],
                         Vm2x         = covar[5,5],
                         Vm1m2x       = covar[7,7],
                         covX_m1x     = covar[1,4],
                         covX_m2x     = covar[1,5],
                         covX_m1m2x   = covar[1,7],
                         covm1x_m2x   = covar[4,5],
                         covm1x_m1m2x = covar[4,7],
                         covm2x_m1m2x = covar[5,7])
    
    seMat[[2]] <- outer(vals[[1]],vals[[3]],
                        seFun,
                        Vx           = covar[2,2],
                        Vm1x         = covar[4,4],
                        Vm2x         = covar[6,6],
                        Vm1m2x       = covar[7,7],
                        covX_m1x     = covar[2,4],
                        covX_m2x     = covar[2,6],
                        covX_m1m2x   = covar[2,7],
                        covm1x_m2x   = covar[4,6],
                        covm1x_m1m2x = covar[4,7],
                        covm2x_m1m2x = covar[6,7])
    
    seMat[[3]] <- outer(vals[[1]],vals[[2]],
                        seFun,
                        Vx           = covar[3,3],
                        Vm1x         = covar[5,5],
                        Vm2x         = covar[6,6],
                        Vm1m2x       = covar[7,7],
                        covX_m1x     = covar[3,5],
                        covX_m2x     = covar[3,6],
                        covX_m1m2x   = covar[3,7],
                        covm1x_m2x   = covar[5,6],
                        covm1x_m1m2x = covar[5,7],
                        covm2x_m1m2x = covar[6,7])
    
    for (i in 1:3) {
      row.names(thetaMat[[i]]) <- row.names(seMat[[i]]) <- vals[-i][[1]]
      colnames(thetaMat[[i]])  <- colnames(seMat[[i]])  <- vals[-i][[2]]
      ZMat[[i]]   <- thetaMat[[i]] / seMat[[i]]
      pMat[[i]]   <- 2 * pmin(pnorm(ZMat[[i]]), (1 - pnorm(ZMat[[i]])))
      sigMat[[i]] <- pMat[[i]] < alpha
      
      
      
      dat2 <- thetaMat[[i]] |>
        as.data.frame() |>
        rownames_to_column("Var1") |>
        pivot_longer(-Var1, names_to = "Var2", values_to = "Parameter Value") 
      
      sig2 <- sigMat[[i]] |>
        as.data.frame() |>
        rownames_to_column("Var1") |>
        pivot_longer(-Var1, names_to = "Var2", values_to = "value") 
      
      dat2$Var1 <- as.numeric(dat2$Var1)
      dat2$Var2 <- as.numeric(dat2$Var2)
      sig2$Var1 <- as.numeric(dat2$Var1)
      sig2$Var2 <- as.numeric(dat2$Var2)
      
      if (crosshatch_non_sig) {
        sig2$pattern <- ifelse(sig2$value, "none",'crosshatch')
      } else {
        sig2$pattern <- ifelse(sig2$value, "crosshatch",'none')
      }
      
      sig3 <- sig2[!sig2$value,]
      
      
      figures[[i]] <- ggplot(dat2, aes(Var1, Var2)) +
        geom_tile(aes(fill = `Parameter Value`)) +
        scale_color_identity() +
        scale_fill_gradient2(low = color_low, 
                             high = color_high,
                             mid = color_mid,
                             midpoint = 0) +
        geom_tile_pattern(data = sig3, aes(pattern = pattern),
                          pattern_density = grid_density,
                          pattern_spacing = grid_spacing,
                          pattern_color = color_grid,
                          alpha = 0) +
        ggtitle(name[i]) +
        theme_bw() +
        guides(pattern = "none") + 
        xlab(name[-i][1]) +
        ylab(name[-i][2])
      
      if (all(c(length(vals[[1]]) < 8,
                length(vals[[2]]) < 8,
                length(vals[[3]]) < 8))) {
        
        figures[[i]] <- figures[[i]] + geom_text(aes(
          label = round(`Parameter Value`, round_res)),
          color = color_values) 
      } 
      
    }
    return_list <- list(thetas = thetaMat,
                        standard_errors = seMat,
                        p_values = pMat,
                        significance = sigMat,
                        plots = figures)
  }
  return(return_list)
}