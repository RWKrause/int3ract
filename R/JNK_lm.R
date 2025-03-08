#' Johnson-Neyman or Johnson-Neyman-Krause plot for lm or glm outputs
#'
#' @param modelOut lm or glm output
#' @param theta1 character; name of the first variable involved in the interaction
#' @param theta2 character; name of the second variable involved in the interaction
#' @param theta3 character; name of the third variable involved in the interaction
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
#' @param grid_density numeric; sets density of grid 
#' @param grid_spacing numeric; sets spacing of grid 
#' @param crosshatch_non_sig logical; should insignificant area be hashed out? If FALSE, significant area will be hashed instead.
#'
#' @returns Returns a list with tables for each of the moderated variables and plots for each.
#' @export
#'
#' @examples
#' \dontrun{
#' dat <- data.frame(y = rnorm(100),
#'                   x = rnorm(100),
#'                   z = rnorm(100),
#'                   w = rnorm(100))
#' 
#' res <- lm(y ~ x * z * w,dat)
#' 
#' x2  <- JN.lm(res,
#'              theta1 = 'x',
#'              theta2 = 'z')
#' x2$plots$z
#' 
#' x3  <- JN.lm(res,
#'              theta1 = 'x',
#'              theta2 = 'z',
#'              theta3 = 'w')
#' x3$plots$z
#' }
JNK_lm <- function(modelOut, 
                  theta1,
                  theta2, 
                  theta3 = NULL,
                  control_fdr = FALSE,
                  alpha = 0.05,
                  round_res = 3,
                  range_size = NULL,
                  sig_color = 'seagreen3',
                  non_sig_color = 'chocolate',
                  line_color = 'black',
                  color_mid = 'white',
                  color_low = '#F05039',
                  color_high = '#000066',
                  color_grid = 'black',
                  grid_density = 0.01,
                  grid_spacing = 0.1,
                  crosshatch_non_sig = TRUE) { 
  
  covT <- vcov(modelOut)
  
  if (is.null(range_size)) {
    if (is.null(theta3)) {
      range_size <- 1000
    } else {
      range_size <- 50
    }
  }
  
  if (is.null(theta3)) {
    vals <- list(range(modelOut$model[[theta1]], na.rm = TRUE),
                 range(modelOut$model[[theta2]], na.rm = TRUE))
    coefs <- modelOut$coefficients[c(theta1,theta2,paste0(theta1,':',theta2))]
    name <- c(theta1,theta2)
    covar <- covT[c(theta1,theta2,paste0(theta1,':',theta2)),
                  c(theta1,theta2,paste0(theta1,':',theta2))]
  } else {
    vals <- list(range(modelOut$model[[theta1]], na.rm = TRUE),
                 range(modelOut$model[[theta2]], na.rm = TRUE),
                 range(modelOut$model[[theta3]], na.rm = TRUE))
    coefs <- modelOut$coefficients[c(theta1,theta2,theta3,
                                     paste0(theta1,':',theta2),
                                     paste0(theta1,':',theta3),
                                     paste0(theta2,':',theta3),
                                     paste(theta1,theta2,theta3, sep = ':'))]
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
           alpha = alpha,
           control_fdr = control_fdr,
           round_res = round_res,
           range_size = range_size,
           sig_color = sig_color,
           non_sig_color = non_sig_color,
           line_color = line_color,
           color_mid = color_mid,
           color_low = color_low,
           color_high = color_high,
           color_grid = color_grid,
           grid_density = grid_density,
           grid_spacing = grid_spacing,
           crosshatch_non_sig = crosshatch_non_sig)
}