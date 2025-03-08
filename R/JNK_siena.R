#' Title
#'
#' @param modelOut 
#' @param theta1 integer; number of the first parameter involved in the interaction
#' @param theta2 integer; number of the second parameter involved in the interaction
#' @param theta3 integer; number of the third parameter involved in the interaction
#' @param thetaInt12 integer; number of the interaction between the first and second parameter
#' @param thetaInt13 integer; number of the interaction between the first and third parameter
#' @param thetaInt23 integer; number of the interaction between the second and third parameter
#' @param thetaInt123 integer; number of the 3-way interaction 
#' @param theta1vals numeric; change statistic values for the first parameter
#' @param theta2vals numeric; change statistic values for the second parameter
#' @param theta3vals numeric; change statistic values for the third parameter
#' @param control_fdr logical; should Bonferroni-Holms correction be used
#' @param alpha numeric; what is the alpha-level
#' @param round_res integer; to which level should results be rounded
#' @param range_size integer; number of calculated parameter combinations. Default for 2-way is 1000, default for 3-way is 50. Note that for 3-way \code{range_size^2} values will be calculated.
#' @param sig_color character; sets color for significant area in 2-way plot
#' @param non_sig_color character; sets color for insignificant area in 2-way plot
#' @param line_color character; sets color for line in 2-way plot
#' @param color_mid character; sets color for \code{0} value in 3-way plot
#' @param color_low character; sets color for negative value in 3-way plot
#' @param color_high character; sets color for positive values in 3-way plot
#' @param color_grid character; sets color of grid 
#' @param color_values character; sets color of moderated parameter values, only active when range_size is 7 or less
#' @param grid_density numeric; sets density of grid 
#' @param grid_spacing numeric; sets spacing of grid 
#' @param crosshatch_non_sig logical; should insignificant area be hashed out? If FALSE, significant area will be hashed instead.
#'
#' @returns Returns a list with tables for each of the moderated parameters and plots for each.
#' @export
#'
#' @examples
#' \dontrun{
#' 
#' library(RSiena)
#' 
#' friendship <- sienaDependent(array(
#'   c(s501,s502,s503),
#'   dim = c(50,50,3)), allowOnly = FALSE)
#' 
#' smoke <- varCovar(s50s[,1:2],centered = TRUE)
#' 
#' netDynamics <- sienaDataCreate(friendship,smoke)
#' 
#' first.model <- getEffects(netDynamics)
#' first.model <- includeEffects(first.model,
#'                               transTrip1)
#' first.model <- includeEffects(first.model,
#'                               egoX,altX,egoXaltX,
#'                               interaction1 = 'smoke')
#' first.model <- includeInteraction(first.model,egoX,transTrip1,
#'                                   interaction1 = c('smoke',''))
#' first.model <- includeInteraction(first.model,altX,transTrip1,
#'                                   interaction1 = c('smoke',''))
#' first.model <- includeInteraction(first.model,egoXaltX,transTrip1,
#'                                   interaction1 = c('smoke',''))
#' 
#' estimation.options <- sienaAlgorithmCreate(useStdInits = FALSE,
#'                                            cond = FALSE,
#'                                            seed = 1402, 
#'                                            n3 = 50000,
#'                                            lessMem = TRUE)
#' 
#' res <- siena07(estimation.options,
#'                data = netDynamics,
#'                effects = first.model,
#'                batch = FALSE,verbose = FALSE,
#'                returnDeps = TRUE,
#'                nbrNodes = 10,
#'                useCluster = TRUE)
#' 
#' x2 <- JNK_siena(res,
#'                 theta1 = 6,
#'                 theta2 = 7, 
#'                 thetaInt12 = 8,
#'                 theta1vals = c(-10:10), 
#'                 theta2vals = c(-10:10))
#' 
#' x2$plots$`smoke alter`
#' 
#' x3 <- JNK_siena(res,
#'                 theta1 = 6,
#'                 theta2 = 7, 
#'                 theta3 = 5, 
#'                 thetaInt12 = 8,
#'                 thetaInt13 = 10,
#'                 thetaInt23 = 9,
#'                 thetaInt123 = 11,
#'                 theta1vals = c(-10:10), 
#'                 theta2vals = c(-10:10), 
#'                 theta3vals = c(0:6),
#'                 range_size = 10)
#' 
#' x3$plots$`smoke alter`
#' }
#' 
JNK_siena <- function(modelOut, 
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
  
  sn <- modelOut$effects$effectName
  
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
        'Numbers need to be the number of the effect in the modelOut object.\n',
        'See most left column of your siena07() result call.\n\n')
    stop()
  }
  
  
  covT <- modelOut$covtheta
  
  if (is.null(theta3)) {
    vals <- list(theta1vals, theta2vals)
    coefs <- c(modelOut$theta[theta1],
               modelOut$theta[theta2],
               modelOut$theta[thetaInt12])
    name <- c(sn[theta1],sn[theta2])
    covar <- covT[c(theta1,theta2,thetaInt12),
                  c(theta1,theta2,thetaInt12)]
  } else {
    vals <- list(theta1vals, theta2vals, theta3vals)
    coefs <- c(modelOut$theta[theta1],
               modelOut$theta[theta2],
               modelOut$theta[theta3],
               modelOut$theta[thetaInt12],
               modelOut$theta[thetaInt13],
               modelOut$theta[thetaInt23],
               modelOut$theta[thetaInt123])
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
