JNK_sienaFit <- function(sienaFit, 
                        theta1,
                        theta2, 
                        theta3 = NULL,
                        thetaInt12 = NULL, 
                        thetaInt13 = NULL,
                        thetaInt23 = NULL,
                        thetaInt123 = NULL, 
                        theta1vals = NULL, 
                        theta2vals = NULL,  
                        theta3vals = NULL, 
                        use_range_only = TRUE,
                        range_size = 500,
                        control_fdr = FALSE,
                        alpha = 0.05,
                        round_res = 3,
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
                        crosshatch_non_sig = TRUE) { 
  
  sn <- sienaFit$effects$effectName
  
  if (any(c(theta1 ,theta2, theta3, thetaInt12, thetaInt13, thetaInt23, 
            thetaInt123) > length(sn))) {
    cat('The following parameter numbers are incorrect:\n',
        if (theta1 > length(sn)) {'theta1\n'},
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
  
  
  covT <- sienaFit$covtheta
  
  if (is.null(theta3)) {
    vals <- list(theta1vals, theta2vals)
    coefs <- c(sienaFit$theta[theta1],
               sienaFit$theta[theta2],
               sienaFit$theta[thetaInt12])
    name <- c(sn[theta1],sn[theta2])
    covar <- covT[c(theta1,theta2,thetaInt12),
                  c(theta1,theta2,thetaInt12)]
  } else {
    vals <- list(theta1vals, theta2vals, theta3vals)
    coefs <- c(sienaFit$theta[theta1],
               sienaFit$theta[theta2],
               sienaFit$theta[theta3],
               sienaFit$theta[thetaInt12],
               sienaFit$theta[thetaInt13],
               sienaFit$theta[thetaInt23],
               sienaFit$theta[thetaInt123])
    name <- c(sn[theta1],sn[theta2],sn[theta3])
    covar <- covT[c(theta1,theta2,theta3,thetaInt12,
                    thetaInt13,thetaInt23,thetaInt123),
                  c(theta1,theta2,theta3,thetaInt12,
                    thetaInt13,thetaInt23,thetaInt123)]
  }
  
  
  freq_int(name = name,
           covar = covar,
           coefs = coefs,
           vals = vals,
           range_size = range_size,
           alpha = alpha,
           control_fdr = control_fdr,
           round_res = round_res,
           sig_color = sig_color,
           non_sig_color = non_sig_color,
           line_color = line_color,
           color_mid = color_mid,
           color_low = color_low,
           color_high = color_high,
           color_values = color_values,
           color_grid = color_grid,
           grid_density = grid_density,
           grid_spacing = grid_spacing,
           crosshatch_non_sig = crosshatch_non_sig)
}
