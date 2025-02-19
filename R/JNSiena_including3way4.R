#with3way
JNSiena <- function(siena07out, 
                    theta1,
                    theta2, 
                    theta3 = NULL,
                    thetaInt12, 
                    thetaInt13 = NULL,
                    thetaInt23 = NULL,
                    thetaInt123 = NULL, 
                    theta1Vals, 
                    theta2Vals,  
                    theta3Vals = NULL, 
                    use_range_only = TRUE,
                    range_size = 50,
                    control_fdr = FALSE,
                    alpha = 0.05,
                    round_res = 3,
                    sig_color_2 = 'seagreen3',
                    non_sig_color_2 = 'chocolate',
                    line_color_2 = 'black',
                    color_mid = 'white',
                    color_low = '#F05039',
                    color_high = '#000066',
                    color_values = 'grey40',
                    color_grid = 'black',
                    grid_density = 0.01,
                    grid_spacing = 0.1,
                    crosshatch_non_sig = TRUE) { 
  
  
  library(ggplot2)
  library(ggpattern)
  library(tidyverse)
  library(scales)
  
  if (class(siena07out) != 'sienaFit') {
    stop('sienaout needs to be a sienaFit object created by siena07().')
  }
  
  sn <- siena07out$effects$effectName
  
  if (any(c(theta1 ,theta2, theta3, thetaInt12, thetaInt13, thetaInt23, 
            thetaInt123) > length(sn))) {
    cat('The following parameter numbers are incorrect:\n',
        if (theta1 >length(sn)) {'theta1\n'},
        if (theta2 > length(sn)) {'theta2\n'},
        if (theta3 > length(sn)) {'theta3\n'},
        if (thetaInt12 > length(sn)) {'thetaInt12\n'},
        if (thetaInt13 > length(sn)) {'thetaInt13\n'},
        if (thetaInt23 > length(sn)) {'thetaInt23\n'},
        if (thetaInt123 > length(sn)) {'thetaInt123\n'},
        'Numbers need to be the number of the effect in the sienaFit object.\n',
        'See most left column of your siena07() result call.\n\n')
    stop()
  }
  
  theta1n <- sn[theta1]
  theta2n <- sn[theta2]
  
  covT <- siena07out$covtheta
  thetas <- siena07out$theta
  
  if (use_range_only) {
    theta1Vals <- seq(min(theta1Vals, na.rm = TRUE), 
                      max(theta1Vals, na.rm = TRUE),
                      length.out = range_size)
    theta2Vals <- seq(min(theta2Vals, na.rm = TRUE), 
                      max(theta2Vals, na.rm = TRUE),
                      length.out = range_size)
    if (!is.null(theta3Vals)) {
      theta3Vals <- seq(min(theta3Vals, na.rm = TRUE), 
                        max(theta3Vals, na.rm = TRUE),
                        length.out = range_size)
    }
  }
  
  if (is.null(theta3)) {
    theta1s <- thetas[theta1] + thetas[thetaInt12] * theta2Vals
    seT1 <- sqrt(covT[theta1,theta1] + 
                   theta2Vals * 2 * covT[thetaInt12,theta1] + 
                   theta2Vals ^ 2 * covT[thetaInt12,thetaInt12])
    
    z1 <- theta1s/seT1
    p1 <- c()
    for (i in 1:length(theta2Vals)) {
      p1[i] <-  2 * pmin(pnorm(z1[[i]]), (1 - pnorm(z1[[i]])))
    }
    if (control_fdr) {
      p1O <- p1[order(p1)]
      bhT1 <- alpha * c(1:length(theta2Vals))/length(theta2Vals)
      sig1 <- p1O < bhT1
      
      if (!all(sig1)) {
        first_false <- which.max(!sig1)
        sig1[first_false:length(theta2Vals)] <- FALSE
      }

      
      
      
      sig11 <- vector(length = length(theta2Vals))
      for (i in 1:length(theta2Vals)) {
        sig11[order(p1)[i]] <- sig1[i]
      }
    } else {
      sig11 <- p1 < alpha
    }
    
    
    
    t1d <- data.frame(theta      = rep(theta1n,length(theta2Vals)),
                      moderator  = rep(theta2n,length(theta2Vals)),
                      mod_vals   = round(theta2Vals,round_res),
                      theta_vals = round(theta1s,round_res),
                      theta_se   = round(seT1,round_res),
                      theta_p    = round(p1,3),
                      sig = sig11)
    
    
    theta2s <- thetas[theta2] + thetas[thetaInt12] * theta1Vals
    seT2 <- sqrt(covT[theta2,theta2] + 
                   theta1Vals * 2 * covT[thetaInt12,theta2] + 
                   theta1Vals ^ 2 * covT[thetaInt12,thetaInt12])
    
    z2 <- theta2s/seT2
    p2 <- c()
    for (i in 1:length(theta1Vals)) {
      p2[i] <-  2 * pmin(pnorm(z2[[i]]), (1 - pnorm(z2[[i]])))
    }
    
    
    
    for (i in 1:length(theta2Vals)) {
      p1[i] <-  2 * pmin(pnorm(z1[[i]]), (1 - pnorm(z1[[i]])))
    }
    if (control_fdr) {
      p2O <- p2[order(p2)]
      bhT2 <- alpha * c(1:length(theta1Vals))/length(theta1Vals)
      sig2 <- p2O < bhT2
      
      if (!all(sig2)) {
        first_false <- which.max(!sig2)
        sig2[first_false:length(theta1Vals)] <- FALSE
      }
      
      
      
      sig22 <- vector(length = length(theta1Vals))
      for (i in 1:length(theta1Vals)) {
        sig22[order(p2)[i]] <- sig2[i]
      }
    } else {
      sig22 <- p2 < alpha
    }
    
    t2d <- data.frame(theta      = rep(theta2n,length(theta1Vals)),
                      moderator  = rep(theta1n,length(theta1Vals)),
                      mod_vals   = round(theta1Vals,round_res),
                      theta_vals = round(theta2s,round_res),
                      theta_se   = round(seT2,round_res),
                      theta_p    = round(p2,3),
                      sig = sig22)
    
    tables <- list(t1d,t2d)
    plots <- vector(mode = 'list', length = 2)
    names(plots) <- names(tables) <- c(theta1n,theta2n)
    
    for (i in 1:2) {
      plots[[i]] <- ggplot(data = tables[[i]], 
                           aes(x = mod_vals, 
                               y = theta_vals)) + 
        geom_hline(aes(yintercept = 0)) +
        theme_bw() + 
        ylab(paste0('Theta for ',names(plots)[i])) +
        xlab(names(plots)[-i]) +
        ggtitle(paste0('Moderation analysis for ',names(plots)[i]))
      
      if (!any(tables[[i]]$sig)) {
        plots[[i]] <- plots[[i]] + 
          geom_ribbon(aes(ymin = -abs(qnorm(alpha / 2)) * theta_se + theta_vals,
                          ymax =  abs(qnorm(alpha / 2)) * theta_se + theta_vals),
                      fill = non_sig_color_2,
                      alpha = 0.7)
      } else if (all(tables[[i]]$sig)) {
        plots[[i]] <- plots[[i]] + 
          geom_ribbon(aes(ymin = -abs(qnorm(alpha / 2)) * theta_se + theta_vals,
                          ymax =  abs(qnorm(alpha / 2)) * theta_se + theta_vals),
                      fill = sig_color_2,
                      alpha = 0.7)
      } else {
        nvals <- nrow(tables[[i]])
        if (tables[[i]]$sig[1]) {
          first_false <- which.max(!tables[[i]]$sig)
          last_false  <- nrow(tables[[i]]) - which.max(rev(!tables[[i]]$sig)) + 1
          
          
          plots[[i]] <- plots[[i]] + 
            geom_ribbon(data = tables[[i]][1:(first_false - 1),],
                        aes(ymin = -abs(qnorm(alpha / 2)) * theta_se + theta_vals,
                            ymax =  abs(qnorm(alpha / 2)) * theta_se + theta_vals),
                        fill = sig_color_2,
                        alpha = 0.7) + 
            geom_ribbon(data = tables[[i]][first_false:last_false,],
                        aes(ymin = -abs(qnorm(alpha / 2)) * theta_se + theta_vals,
                            ymax =  abs(qnorm(alpha / 2)) * theta_se + theta_vals),
                        fill = non_sig_color_2,
                        alpha = 0.7)
          if (tables[[i]]$sig[nvals]) {
            plots[[i]] <- plots[[i]] + 
              geom_ribbon(data = tables[[i]][(last_false + 1):nvals,],
                          aes(ymin = -abs(qnorm(alpha / 2)) * theta_se + theta_vals,
                              ymax =  abs(qnorm(alpha / 2)) * theta_se + theta_vals),
                          fill = sig_color_2,
                          alpha = 0.7)
          }
        } else {
          first_true <- which.max(tables[[i]]$sig)
          last_true  <- nrow(tables[[i]]) - which.max(rev(tables[[i]]$sig)) + 1
          
          
          plots[[i]] <- plots[[i]] + 
            geom_ribbon(data = tables[[i]][1:(first_true - 1),],
                        aes(ymin = -abs(qnorm(alpha / 2)) * theta_se + theta_vals,
                            ymax =  abs(qnorm(alpha / 2)) * theta_se + theta_vals),
                        fill = non_sig_color_2,
                        alpha = 0.7) + 
            geom_ribbon(data = tables[[i]][first_true:last_true,],
                        aes(ymin = -abs(qnorm(alpha / 2)) * theta_se + theta_vals,
                            ymax =  abs(qnorm(alpha / 2)) * theta_se + theta_vals),
                        fill = sig_color_2,
                        alpha = 0.7)
          if (!tables[[i]]$sig[nvals]) {
            plots[[i]] <- plots[[i]] + 
              geom_ribbon(data = tables[[i]][(last_true + 1):nvals,],
                          aes(ymin = -abs(qnorm(alpha / 2)) * theta_se + theta_vals,
                              ymax =  abs(qnorm(alpha / 2)) * theta_se + theta_vals),
                          fill = non_sig_color_2,
                          alpha = 0.7)
          }
        }
      }
      plots[[i]] <- plots[[i]] +  
        geom_path(linewidth = 1, color = line_color_2) 
    }
    
    
    
    return_list <- list(tables = tables,
                        plots = plots)
    
  } else {
    theta3n <- sn[theta3]
    
    thetaMat <- vector(mode = 'list', length = 3)
    seMat    <- vector(mode = 'list', length = 3)
    ZMat     <- vector(mode = 'list', length = 3)
    pMat     <- vector(mode = 'list', length = 3)
    sigMat   <- vector(mode = 'list', length = 3)
    figures  <- vector(mode = 'list', length = 3)
    names(thetaMat) <- names(seMat) <- names(figures) <- c(theta1n,theta2n,
                                                           theta3n)
    names(ZMat) <- names(pMat) <- names(sigMat) <- c(theta1n,theta2n,theta3n)
    
    
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
    
    
    
    thetaMat[[1]] <- outer(theta2Vals,theta3Vals,
                           parFun,
                           bx     = thetas[theta1],
                           bm1x   = thetas[thetaInt12],
                           bm2x   = thetas[thetaInt13],
                           bm1m2x = thetas[thetaInt123])
    
    thetaMat[[2]] <- outer(theta1Vals,theta3Vals,
                           parFun,
                           bx     = thetas[theta2],
                           bm1x   = thetas[thetaInt12],
                           bm2x   = thetas[thetaInt23],
                           bm1m2x = thetas[thetaInt123])
    
    thetaMat[[3]] <- outer(theta1Vals,theta2Vals,
                           parFun,
                           bx     = thetas[theta3],
                           bm1x   = thetas[thetaInt13],
                           bm2x   = thetas[thetaInt23],
                           bm1m2x = thetas[thetaInt123])
    
    
    seMat[[1]] <-  outer(theta2Vals,theta3Vals,
                         seFun,
                         Vx           = covT[theta1,theta1],
                         Vm1x         = covT[thetaInt12,thetaInt12],
                         Vm2x         = covT[thetaInt13,thetaInt13],
                         Vm1m2x       = covT[thetaInt123,thetaInt123],
                         covX_m1x     = covT[theta1,thetaInt12],
                         covX_m2x     = covT[theta1,thetaInt13],
                         covX_m1m2x   = covT[theta1,thetaInt123],
                         covm1x_m2x   = covT[thetaInt12,thetaInt13],
                         covm1x_m1m2x = covT[thetaInt12,thetaInt123],
                         covm2x_m1m2x = covT[thetaInt13,thetaInt123])
    
    seMat[[2]] <- outer(theta1Vals,theta3Vals,
                        seFun,
                        Vx           = covT[theta2,theta2],
                        Vm1x         = covT[thetaInt12,thetaInt12],
                        Vm2x         = covT[thetaInt23,thetaInt23],
                        Vm1m2x       = covT[thetaInt123,thetaInt123],
                        covX_m1x     = covT[theta2,thetaInt12],
                        covX_m2x     = covT[theta2,thetaInt23],
                        covX_m1m2x   = covT[theta2,thetaInt123],
                        covm1x_m2x   = covT[thetaInt12,thetaInt23],
                        covm1x_m1m2x = covT[thetaInt12,thetaInt123],
                        covm2x_m1m2x = covT[thetaInt23,thetaInt123])
    
    seMat[[3]] <- outer(theta1Vals,theta2Vals,
                        seFun,
                        Vx           = covT[theta3,theta3],
                        Vm1x         = covT[thetaInt13,thetaInt13],
                        Vm2x         = covT[thetaInt23,thetaInt23],
                        Vm1m2x       = covT[thetaInt123,thetaInt123],
                        covX_m1x     = covT[theta3,thetaInt13],
                        covX_m2x     = covT[theta3,thetaInt23],
                        covX_m1m2x   = covT[theta3,thetaInt123],
                        covm1x_m2x   = covT[thetaInt13,thetaInt23],
                        covm1x_m1m2x = covT[thetaInt13,thetaInt123],
                        covm2x_m1m2x = covT[thetaInt23,thetaInt123])
    
    vals <- list(theta1Vals,theta2Vals,theta3Vals)
    ns <- list(theta1n,theta2n,theta3n)
    for (i in 1:3) {
      #naming with parameter name might interfere with plot
      # thus only values so far, otherwise:paste(ns[-i][[1]],vals[-i][[1]]) 
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
        ggtitle(ns[i]) +
        theme_bw() +
        guides(pattern = "none") + 
        xlab(ns[-i][[1]]) +
        ylab(ns[-i][[2]])
      
      if (all(c(length(theta1Vals) < 8,
                length(theta3Vals) < 8,
                length(theta2Vals) < 8))) {
        
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
