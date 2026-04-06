#' Create Johnson-Neyman plots for Bayesian SAOMs
#'
#' @param sbo multiSiena; multiSiena() output
#' @param theta_1 numeric; number of the parameter in the effects object. All numbers ignore rates! see 'sbo$effects[sbo$effects$type != 'rate',]'
#' @param theta_2 numeric; number of the parameter in the effects object. All numbers ignore rates! see 'sbo$effects[sbo$effects$type != 'rate',]'
#' @param theta_3 numeric; number of the parameter in the effects object. All numbers ignore rates! see 'sbo$effects[sbo$effects$type != 'rate',]'. Default NULL
#' @param theta_int_12 numeric; number of the interaction between theta_1 and theta_2 in the effects object. All numbers ignore rates! see 'sbo$effects[sbo$effects$type != 'rate',]'
#' @param theta_int_13 numeric; number of the interaction between theta_1 and theta_2 in the effects object. All numbers ignore rates! see 'sbo$effects[sbo$effects$type != 'rate',]'. Default NULL
#' @param theta_int_23 numeric; number of the interaction between theta_1 and theta_2 in the effects object. All numbers ignore rates! see 'sbo$effects[sbo$effects$type != 'rate',]'. Default NULL
#' @param theta_int_123 numeric; number of the interaction between theta_1, theta_2, and theta_3 in the effects object. All numbers ignore rates! see 'sbo$effects[sbo$effects$type != 'rate',]'. Default NULL
#' @param theta_1_vals numeric; possible range of the statistic the parameter is multiplied with.
#' @param theta_2_vals numeric; possible range of the statistic the parameter is multiplied with.
#' @param theta_3_vals numeric; possible range of the statistic the parameter is multiplied with. Default NULL
#' @param burnIn numeric; burn-in iteration to be removed. Default NULL
#' @param thin numeric; Should the posteriors be thinned? Default = 1, i.e., every iteration is used. Thin = 2 means every 2nd iteration is used.
#' @param thresholds numeric; threshold for 'significance' hashing. (Note that all concerns with regards to 'significance' in Bayesian models apply.) Default = NULL
#' @param hyper_only logical; should only the hyper parameter be used or plots be made for the random parameters?
#' @param round_res numeric; to how many digits should results be rounded? Default = 3
#' @param noTitle character; You can provide a title to be used for the plots. Titles can always be changed later by using regular ggplot2 coding.
#' @param color_mid character; Color chosen for the mid-point of the scale. Default = 'white'
#' @param color_low character; color chosen for low values. Default = '#F05039', (red)
#' @param color_high character; color chosen for high values. Default = '#000066', (blue)
#' @param color_values character; color chosen for the numbers. Default = 'grey40' 
#' @param color_grid character; color chosen for the grid. Default = 'black'
#' @param grid_density numeric; density of the hash-grid. Default = 0.01
#' @param grid_spacing numeric; spacing of the hash-grid. Default = 0.1
#' @param save logical; should the plots be saved directly with ggsave(). This is advisable if you have many plots and want to look at them with more ease. See also the folder parameter below. Default = FALSE
#' @param folder character; Where should the plots be saved. If NULL, a new folder for the plots - int3ract JNplots - will be created. If that folder does not exist, RStudio will ask if it should be created.

#'
#' @returns list of plots.
#' @export
#'
#' @import ggplot2
#' @import ggpattern
#' @import scales
#' @import tidyr
#' @import tibble
#'
JNSienaBayes <- function(sbo,
                         theta_1,
                         theta_2, 
                         theta_3 = NULL,
                         theta_int_12,
                         theta_int_13 = NULL,
                         theta_int_23 = NULL,
                         theta_int_123 = NULL,
                         theta_1_vals,
                         theta_2_vals,
                         theta_3_vals = NULL,
                         burnIn = NULL,
                         thin = 1,
                         thresholds = NULL,
                         hyper_only = TRUE,
                         round_res = 3,
                         noTitle = NULL,
                         color_mid = 'white',
                         color_low = '#F05039',
                         color_high = '#000066',
                         color_values = 'grey40',
                         color_grid = 'black',
                         grid_density = 0.01,
                         grid_spacing = 0.1,
                         save = FALSE,
                         folder = NULL) {
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
  
  if (is.null(folder)) {
    folder <- 'int3ract JNplots'
  }
  theta_1n <- eff$effectName[theta_1]
  theta_2n <- eff$effectName[theta_2]
  theta_1n <- gsub('\\^','',theta_1n)
  theta_1n <- gsub('\\/','o',theta_1n)
  theta_1n <- gsub(':','',theta_1n)
  theta_2n <- gsub('\\^','',theta_2n)
  theta_2n <- gsub(':','',theta_2n)
  theta_2n <- gsub('\\/','o',theta_2n)
  
  if (is.null(theta_3)) {
    if (!any(eff$randomEffects[c(theta_1, theta_2, theta_int_12)])) {
      theta <- sbo$ThinParameters[seq(burnIn,
                                      nrow(sbo$ThinParameters),
                                      thin),1,c(theta_1, theta_2, theta_int_12)]
      returnList <- jnb_support2(theta = theta,
                                 group = 'Eta',
                                 sbo = sbo,
                                 theta_1 = theta_1, 
                                 theta_2 = theta_2, 
                                 theta_1n = theta_1n, 
                                 theta_2n = theta_2n, 
                                 theta_int_12 = theta_int_12, 
                                 theta_1_vals = theta_1_vals, 
                                 theta_2_vals = theta_2_vals,
                                 save = save,
                                 folder = folder)
    } else {
      rownames(eff) <- 1:nrow(eff)
      eff_ran <- eff[eff$randomEffects,]
      
      thetas <- c(theta_1,theta_2,
                  theta_int_12)
      
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
                                 theta_1 = theta_1, 
                                 theta_2 = theta_2, 
                                 theta_1n = theta_1n, 
                                 theta_2n = theta_2n, 
                                 theta_int_12 = theta_int_12, 
                                 theta_1_vals = theta_1_vals, 
                                 theta_2_vals = theta_2_vals,
                                 save = save,
                                 folder = folder)
      
      if (!hyper_only) {
        randomResults <- vector(mode = 'list', length = sbo$nGroup)
        names(randomResults) <- paste0('group',1:sbo$nGroup)
        for (i in 1:sbo$nGroup) {
          theta <- sbo$ThinParameters[seq(burnIn,
                                          nrow(sbo$ThinParameters),
                                          thin),
                                      i,
                                      c(theta_1, theta_2, theta_int_12)]
          randomResults[[i]] <- jnb_support2(theta = theta,
                                             group = paste0('group',i),
                                             sbo = sbo,
                                             theta_1 = theta_1, 
                                             theta_2 = theta_2, 
                                             theta_1n = theta_1n, 
                                             theta_2n = theta_2n, 
                                             theta_int_12 = theta_int_12, 
                                             theta_1_vals = theta_1_vals, 
                                             theta_2_vals = theta_2_vals, 
                                             save = save,
                                             folder = folder)
        }
        returnList <- list(returnList,randomResults)
        names(returnList) <- c('Mu','random_groups_effects')
      } 
    }
    
  } else {
    theta_3n <- eff$effectName[theta_3]
    theta_3n <- gsub('\\^','',theta_3n)
    theta_3n <- gsub('\\/','o',theta_3n)
    theta_3n <- gsub(':','',theta_3n)
    if (!any(eff$randomEffects[c(theta_1, theta_2, theta_3, 
                                 theta_int_12, theta_int_13, theta_int_23,
                                 theta_int_123)])) {
      
      theta <- sbo$ThinParameters[seq(burnIn,
                                      nrow(sbo$ThinParameters),
                                      thin),1,c(theta_1, theta_2, theta_3, 
                                                theta_int_12, theta_int_13, 
                                                theta_int_23, theta_int_123)]
      returnList <- jnb_support3(theta = theta,
                                 group = 'Eta',
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
                                 grid_spacing = grid_spacing, 
                                 save = save,
                                 folder = folder)
    } else {
      rownames(eff) <- 1:nrow(eff)
      eff_ran <- eff[eff$randomEffects,]
      
      thetas <- c(theta_1,theta_2,theta_3,
                  theta_int_12,theta_int_13,theta_int_23,
                  theta_int_123)
      
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
                                 folder = folder)
      
      if (!hyper_only) {
        randomResults <- vector(mode = 'list', length = sbo$nGroup)
        names(randomResults) <- paste0('group',1:sbo$nGroup)
        for (i in 1:sbo$nGroup) {
          theta <- sbo$ThinParameters[seq(burnIn,
                                          nrow(sbo$ThinParameters),
                                          thin),
                                      i,thetas]
          randomResults[[i]] <- jnb_support3(theta = theta,
                                             group = paste0('group',i),
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
                                             grid_spacing = grid_spacing, 
                                             save = save,
                                             folder = folder)
        }
        returnList <- list(returnList,randomResults)
        names(returnList) <- c('Mu','random_groups_effects')
      } 
    }
  }
  
  
  return(returnList)
}