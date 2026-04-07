#' Title
#'
#' @param modelOut 
#' @param theta_1 integer; number of the first parameter involved in the interaction
#' @param theta_2 integer; number of the second parameter involved in the interaction
#' @param theta_3 integer; number of the third parameter involved in the interaction
#' @param theta_int_12 integer; number of the interaction between the first and second parameter
#' @param theta_int_13 integer; number of the interaction between the first and third parameter
#' @param theta_int_23 integer; number of the interaction between the second and third parameter
#' @param theta_int_123 integer; number of the 3-way interaction 
#' @param theta_1vals numeric; change statistic values for the first parameter
#' @param theta_2vals numeric; change statistic values for the second parameter
#' @param theta_3vals numeric; change statistic values for the third parameter
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
#'                 theta_1 = 6,
#'                 theta_2 = 7, 
#'                 theta_int_12 = 8,
#'                 theta_1vals = c(-10:10), 
#'                 theta_2vals = c(-10:10))
#' 
#' x2$plots$`smoke alter`
#' 
#' x3 <- JNK_siena(res,
#'                 theta_1 = 6,
#'                 theta_2 = 7, 
#'                 theta_3 = 5, 
#'                 theta_int_12 = 8,
#'                 theta_int_13 = 10,
#'                 theta_int_23 = 9,
#'                 theta_int_123 = 11,
#'                 theta_1vals = c(-10:10), 
#'                 theta_2vals = c(-10:10), 
#'                 theta_3vals = c(0:6),
#'                 range_size = 10)
#' 
#' x3$plots$`smoke alter`
#' }
#' 
JNK_siena <- function(modelOut, 
                      theta_1,
                      theta_2, 
                      theta_3 = NULL,
                      theta_int_12 = NULL, 
                      theta_int_13 = NULL,
                      theta_int_23 = NULL,
                      theta_int_123 = NULL, 
                      theta_1vals = NULL, 
                      theta_2vals = NULL,  
                      theta_3vals = NULL, 
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
  
  if (any(c(theta_1 ,theta_2, theta_3, theta_int_12, theta_int_13, theta_int_23, 
            theta_int_123) > length(sn))) {
    cat('The following parameter numbers are incorrect:\n',
        if (theta_1 > length(sn)) {'theta_1\n'},
        if (theta_2 > length(sn)) {'theta_2\n'},
        if (theta_3 > length(sn)) {'theta_3\n'},
        if (theta_int_12 > length(sn)) {'theta_int_12\n'},
        if (theta_int_13 > length(sn)) {'theta_int_13\n'},
        if (theta_int_23 > length(sn)) {'theta_int_23\n'},
        if (theta_int_123 > length(sn)) {'theta_int_123\n'},
        'Numbers need to be the number of the effect in the modelOut object.\n',
        'See most left column of your siena07() result call.\n\n')
    stop()
  }
  
  
  covT <- modelOut$covtheta
  
  if (is.null(theta_3)) {
    vals <- list(theta_1vals, theta_2vals)
    coefs <- c(modelOut$theta[theta_1],
               modelOut$theta[theta_2],
               modelOut$theta[theta_int_12])
    name <- c(sn[theta_1],sn[theta_2])
    covar <- covT[c(theta_1,theta_2,theta_int_12),
                  c(theta_1,theta_2,theta_int_12)]
  } else {
    vals <- list(theta_1vals, theta_2vals, theta_3vals)
    coefs <- c(modelOut$theta[theta_1],
               modelOut$theta[theta_2],
               modelOut$theta[theta_3],
               modelOut$theta[theta_int_12],
               modelOut$theta[theta_int_13],
               modelOut$theta[theta_int_23],
               modelOut$theta[theta_int_123])
    name <- c(sn[theta_1],sn[theta_2],sn[theta_3])
    covar <- covT[c(theta_1,theta_2,theta_3,theta_int_12,
                    theta_int_13,theta_int_23,theta_int_123),
                  c(theta_1,theta_2,theta_3,theta_int_12,
                    theta_int_13,theta_int_23,theta_int_123)]
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
