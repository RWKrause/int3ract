JN.lm <- function(modelOut, 
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
  
  covT <- vcov(modelOut)
  
  if (is.null(theta3)) {
    vals <- list(range(modelOut$model[[theta1]], na.rm = TRUE),
                 range(modelOut$model[[theta2]], na.rm = TRUE))
    coefs <- modelOut$coefficients[c(theta1,theta2)]
    name <- c(theta1,theta2)
    covar <- covT[c(theta1,theta2,paste0(theta1,':',theta2)),
                  c(theta1,theta2,paste0(theta1,':',theta2))]
  } else {
    vals <- list(range(modelOut$model[[theta1]], na.rm = TRUE),
                 range(modelOut$model[[theta2]], na.rm = TRUE),
                 range(modelOut$model[[theta3]], na.rm = TRUE))
    coefs <- modelOut$coefficients[c(theta1,theta2,theta3)]
    name <- c(theta1,theta2,theta3)
    covar <- covT[c(theta1,theta2,theta3,
                    paste0(theta1,':',theta2),
                    paste0(theta1,':',theta3),
                    paste0(theta2,':',theta3),
                    paste(theta1,theta2,theta3, sep = ':')),
                  c(theta1,theta2,theta3,
                    paste0(theta1,':',theta2),
                    paste0(theta1,':',theta3),
                    paste0(theta2,':',theta3),
                    paste(theta1,theta2,theta3, sep = ':'))]
  }
  
  
  freq_int(name = name,
           covar = covar,
           coefs = coefs,
           vals = vals,
           use_range_only = use_range_only,
           range_size = range_size,
           alpha = alpha,
           control_fdr = control_fdr,
           round_res = round_res,
           sig_color_2 = sig_color_2,
           non_sig_color_2 = non_sig_color_2,
           line_color_2 = line_color_2,
           color_mid = color_mid,
           color_low = color_low,
           color_high = color_high,
           color_values = color_values,
           color_grid = color_grid,
           grid_density = grid_density,
           grid_spacing = grid_spacing,
           crosshatch_non_sig = crosshatch_non_sig)
}
    
    
