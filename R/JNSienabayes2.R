# all numbers ignore rates!
# see sbo$effects[sbo$effects$type != 'rate',]

JNSienaBayes <- function(sbo, # sienaBayes output
                         theta1, #number of the first parameter involved in the 
                         # interaction - check the model results
                         theta2, # number of the 2nd parameter  
                         theta3 = NULL, # number of the 3rd parameter  
                         thetaInt12, # number of the interaction
                         thetaInt13 = NULL, # number of the interaction
                         thetaInt23 = NULL, # number of the interaction
                         thetaInt123 = NULL, # number of the interaction
                         theta1vals, # range of the statistics for the 1st parameter
                         theta2vals,  # range of the statistics for the 2nd parameter 
                         theta3vals = NULL,  # range of the statistics for the 3rd parameter 
                         burnIn = NULL, # how many iterations should be cut as burn in, NULL = nwarm
                         thin = 1, #should iterations be thinned
                         hyperOnly = TRUE,
                         round_res = 3,
                         noTitle = NULL,
                         color_mid = 'white',
                         color_low = '#F05039',
                         color_high = '#000066',
                         color_values = 'grey40',
                         color_grid = 'black',
                         grid_density = 0.01,
                         grid_spacing = 0.1,
                         thresholds = NULL) { # if any of the parameters is random should only the hyper parameter mu be used (otherwise plots and tables for all groups included!)
  eff <- sbo$effects
  
  deps <- unique(eff$name)
  
  if (is.null(burnIn)) {
    burnIn <- sbo$nwarm
    if (burnIn == 0) {
      burnIn <- 1
    }
  }
  
  first_rates <- c()
  for (i in 1:length(deps)) {
    first_rates <- c(first_rates,
                     which(eff$name == deps[i] & 
                             eff$functionName == 'Amount of network change in period 1'))
  }
  
  rates <- sbo$basicRate
  rates[first_rates] <- FALSE
  
  eff <- eff[!rates,]
  
  
  theta1n <- eff$effectName[theta1]
  theta2n <- eff$effectName[theta2]
  theta1n <- gsub('\\^','',theta1n)
  theta1n <- gsub('\\/','o',theta1n)
  theta1n <- gsub(':','',theta1n)
  theta2n <- gsub('\\^','',theta2n)
  theta2n <- gsub(':','',theta2n)
  theta2n <- gsub('\\/','o',theta2n)
  
  
  if (is.null(burnIn)) {
    burnIn <- sbo$nwarm + 1
  } 
  
  source('jnb_support.R')
  
  if (is.null(theta3)) {
    if (!any(eff$randomEffects[c(theta1, theta2, thetaInt12)])) {
      theta <- sbo$ThinParameters[seq(burnIn,
                                      nrow(sbo$ThinParameters),
                                      thin),1,c(theta1, theta2, thetaInt12)]
      returnList <- jnb_support2(theta = theta,
                                 group = 'Eta',
                                 sbo = sbo,
                                 theta1 = theta1, 
                                 theta2 = theta2, 
                                 theta1n = theta1n, 
                                 theta2n = theta2n, 
                                 thetaInt12 = thetaInt12, 
                                 theta1vals = theta1vals, 
                                 theta2vals = theta2vals)
    } else {
      rownames(eff) <- 1:nrow(eff)
      eff_ran <- eff[eff$randomEffects,]
      
      thetas <- c(theta1,theta2,
                  thetaInt12)
      
      theta <- matrix(NA, 
                      length(seq(burnIn,
                                 nrow(sbo$ThinPosteriorMu),
                                 thin)),3)
      
      
      collumn <- 0
      for (i in thetas) {
        collumn <- collumn + 1
        if (eff$randomEffects[i]) {
          theta[,collumn] <- sbo$ThinPosteriorMu[seq(burnIn,
                                                     nrow(sbo$ThinPosteriorMu),
                                                     thin),
                                                 which(as.numeric(
                                                   rownames(eff_ran)) == i)]
        } else {
          theta[,collumn] <- sbo$ThinParameters[seq(burnIn,
                                                    nrow(sbo$ThinParameters),
                                                    thin),1,i]
        }
      }
      
      returnList <- jnb_support2(theta = theta,
                                 group = 'Mu',
                                 sbo = sbo,
                                 theta1 = theta1, 
                                 theta2 = theta2, 
                                 theta1n = theta1n, 
                                 theta2n = theta2n, 
                                 thetaInt12 = thetaInt12, 
                                 theta1vals = theta1vals, 
                                 theta2vals = theta2vals)
      
      if (!hyperOnly) {
        randomResults <- vector(mode = 'list', length = sbo$nGroup)
        names(randomResults) <- paste0('group',1:sbo$nGroup)
        for (i in 1:sbo$nGroup) {
          theta <- sbo$ThinParameters[seq(burnIn,
                                          nrow(sbo$ThinParameters),
                                          thin),
                                      i,
                                      c(theta1, theta2, thetaInt12)]
          randomResults[[i]] <- jnb_support2(theta = theta,
                                             group = paste0('group',i),
                                             sbo = sbo,
                                             theta1 = theta1, 
                                             theta2 = theta2, 
                                             theta1n = theta1n, 
                                             theta2n = theta2n, 
                                             thetaInt12 = thetaInt12, 
                                             theta1vals = theta1vals, 
                                             theta2vals = theta2vals)
        }
        returnList <- list(returnList,randomResults)
        names(returnList) <- c('Mu','random_groups_effects')
      } 
    }
    
  } else {
    theta3n <- eff$effectName[theta3]
    theta3n <- gsub('\\^','',theta3n)
    theta3n <- gsub('\\/','o',theta3n)
    theta3n <- gsub(':','',theta3n)
    if (!any(eff$randomEffects[c(theta1, theta2, theta3, 
                                 thetaInt12, thetaInt13, thetaInt23,
                                 thetaInt123)])) {
      
      theta <- sbo$ThinParameters[seq(burnIn,
                                      nrow(sbo$ThinParameters),
                                      thin),1,c(theta1, theta2, theta3, 
                                                thetaInt12, thetaInt13, 
                                                thetaInt23, thetaInt123)]
      returnList <- jnb_support3(theta = theta,
                                 group = 'Eta',
                                 sbo = sbo,
                                 theta1 = theta1, 
                                 theta2 = theta2, 
                                 theta3 = theta3, 
                                 theta1n = theta1n, 
                                 theta2n = theta2n, 
                                 theta3n = theta3n, 
                                 thetaInt12 = thetaInt12, 
                                 thetaInt13 = thetaInt13, 
                                 thetaInt23 = thetaInt23, 
                                 thetaInt123 = thetaInt123, 
                                 theta1vals = theta1vals, 
                                 theta2vals = theta2vals, 
                                 theta3vals = theta3vals,
                                 thresholds = thresholds,
                                 color_mid = color_mid,
                                 color_low = color_low,
                                 color_high = color_high,
                                 color_values = color_values,
                                 color_grid = color_grid,
                                 grid_density = grid_density,
                                 grid_spacing = grid_spacing)
    } else {
      rownames(eff) <- 1:nrow(eff)
      eff_ran <- eff[eff$randomEffects,]
      
      thetas <- c(theta1,theta2,theta3,
                  thetaInt12,thetaInt13,thetaInt23,
                  thetaInt123)
      
      theta <- matrix(NA, 
                      length(seq(burnIn,
                                 nrow(sbo$ThinPosteriorMu),
                                 thin)),
                      7)
      
      
      collumn <- 0
      for (i in thetas) {
        collumn <- collumn + 1
        if (eff$randomEffects[i]) {
          theta[,collumn] <- sbo$ThinPosteriorMu[seq(burnIn,
                                                     nrow(sbo$ThinPosteriorMu),
                                                     thin),
                                                 which(as.numeric(
                                                   rownames(eff_ran)) == i)]
        } else {
          theta[,collumn] <- sbo$ThinParameters[seq(burnIn,
                                                    nrow(sbo$ThinParameters),
                                                    thin),1,i]
        }
      }
      
      returnList <- jnb_support3(theta = theta,
                                 group = 'Eta',
                                 sbo = sbo,
                                 theta1 = theta1, 
                                 theta2 = theta2, 
                                 theta3 = theta3, 
                                 theta1n = theta1n, 
                                 theta2n = theta2n, 
                                 theta3n = theta3n, 
                                 thetaInt12 = thetaInt12, 
                                 thetaInt13 = thetaInt13, 
                                 thetaInt23 = thetaInt23, 
                                 thetaInt123 = thetaInt123, 
                                 theta1vals = theta1vals, 
                                 theta2vals = theta2vals, 
                                 theta3vals = theta3vals,
                                 thresholds = thresholds,
                                 color_mid = color_mid,
                                 color_low = color_low,
                                 color_high = color_high,
                                 color_values = color_values,
                                 color_grid = color_grid,
                                 grid_density = grid_density,
                                 grid_spacing = grid_spacing)
      
      if (!hyperOnly) {
        randomResults <- vector(mode = 'list', length = sbo$nGroup)
        names(randomResults) <- paste0('group',1:sbo$nGroup)
        for (i in 1:sbo$nGroup) {
          theta <- sbo$ThinParameters[seq(burnIn,
                                          nrow(sbo$ThinParameters),
                                          thin),
                                      i,thetas]
          randomResults[[i]] <- jnb_support3(theta = theta,
                                             group = 'Eta',
                                             sbo = sbo,
                                             theta1 = theta1, 
                                             theta2 = theta2, 
                                             theta3 = theta3, 
                                             theta1n = theta1n, 
                                             theta2n = theta2n, 
                                             theta3n = theta3n, 
                                             thetaInt12 = thetaInt12, 
                                             thetaInt13 = thetaInt13, 
                                             thetaInt23 = thetaInt23, 
                                             thetaInt123 = thetaInt123, 
                                             theta1vals = theta1vals, 
                                             theta2vals = theta2vals, 
                                             theta3vals = theta3vals,
                                             thresholds = thresholds,
                                             color_mid = color_mid,
                                             color_low = color_low,
                                             color_high = color_high,
                                             color_values = color_values,
                                             color_grid = color_grid,
                                             grid_density = grid_density,
                                             grid_spacing = grid_spacing)
        }
        returnList <- list(returnList,randomResults)
        names(returnList) <- c('Mu','random_groups_effects')
      } 
    }
  }
  
  
  return(returnList)
}
