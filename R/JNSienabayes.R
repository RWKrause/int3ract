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
                         color_mid = '#EBCC2A', 
                         color_low = '#3B9AB2',
                         color_high = '#F21A00', 
                         color_values = 'grey40',
                         color_grid = 'black', 
                         grid_density = 0.01,
                         grid_spacing = 0.1, 
                         save = FALSE, 
                         folder = NULL) {
  
  eff    <- sbo$effects
  deps   <- unique(eff$name)
  folder <- folder %||% 'int3ract JNplots'
  
  if (is.null(burnIn)) burnIn <- max(sbo$nwarm, 1)
  
  first_rates <- sapply(deps, function(d) {
    which(eff$name == d & 
            eff$functionName == 'Amount of network change in period 1')})
  rates <- sbo$basicRate
  rates[unlist(first_rates)] <- FALSE
  eff <- eff[!rates, ]
  rownames(eff) <- seq_len(nrow(eff))
  
  clean_name <- function(n) gsub(':', '', gsub('/', 'o', gsub('\\^', '', n)))
  
  theta_1n <- clean_name(eff$effectName[theta_1])
  theta_2n <- clean_name(eff$effectName[theta_2])
  theta_3n <- if (!is.null(theta_3)) {
    clean_name(eff$effectName[theta_3])
  } else {
    NULL
  }
  
  idx <- seq(burnIn, nrow(sbo$ThinParameters), thin)
  
  color_args <- list(color_mid = color_mid, 
                     color_low = color_low,
                     color_high = color_high,
                     color_values = color_values,
                     color_grid = color_grid, 
                     grid_density = grid_density,
                     grid_spacing = grid_spacing)
  
  call_support <- function(theta, group, threeWay) {
    if (threeWay) {
      do.call(jnb_support3,
              c(list(theta = theta, 
                     group = group, 
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
                     save = save, 
                     folder = folder),
                color_args))
    } else {
      do.call(jnb_support2,
              list(theta = theta, 
                   group = group, 
                   theta_1 = theta_1, 
                   theta_2 = theta_2,
                   theta_1n = theta_1n,
                   theta_2n = theta_2n,
                   theta_int_12 = theta_int_12,
                   theta_1_vals = theta_1_vals, 
                   theta_2_vals = theta_2_vals,
                   color_mid = color_mid, 
                   color_low = color_low,
                   color_high = color_high,
                   save = save, 
                   folder = folder))
    }
  }
  
  threeWay <- !is.null(theta_3)
  thetas   <- if (threeWay) {
    c(theta_1, theta_2, theta_3, theta_int_12, 
      theta_int_13, theta_int_23, theta_int_123)
  } else { 
    c(theta_1, theta_2, theta_int_12)
  }
  
  any_random <- any(eff$randomEffects[thetas])
  
  if (!any_random) {
    theta <- sbo$ThinParameters[idx, 1, thetas]
    return(call_support(theta, 'Eta', threeWay))
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
  theta      <- extract_theta(thetas)
  returnList <- call_support(theta, 'Mu', threeWay)
  
  if (!hyper_only) {
    randomResults <- setNames(
      lapply(seq_len(sbo$nGroup), function(i) {
        theta_g <- sbo$ThinParameters[idx, i, thetas]
        call_support(theta_g, paste0('group', i), threeWay)
      }),
      paste0('group', seq_len(sbo$nGroup))
    )
    returnList <- setNames(list(returnList, randomResults),
                           c('Mu', 'random_groups_effects'))
  }
  
  returnList
}