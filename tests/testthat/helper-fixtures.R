# Shared fixtures. The same seed and data-generating process as the manuscript,
# so that a failure here is directly comparable with what the paper reports.

jn_test_data <- function(n = 100, seed = 1402) {
  set.seed(seed)
  d <- data.frame(x = stats::rnorm(n), z = stats::rnorm(n), w = stats::rnorm(n))
  d$y <- d$x + 0.5 * d$z - 0.5 * d$w +
    0.5 * d$x * d$z * d$w + stats::rnorm(n, sd = 4)
  d
}

# Stand-ins for RSiena and multiSiena results. Both packages are heavyweight
# and their fitting routines are far too slow to run in a test suite, so these
# carry just the structure the extractors read.

jn_mock_sienafit <- function(n_eff = 5, seed = 3) {
  set.seed(seed)
  nm <- c("altX", "egoX", "altX x egoX", "recip", "transTrip")[seq_len(n_eff)]
  V  <- crossprod(matrix(stats::rnorm(n_eff * n_eff), n_eff))
  structure(list(effects  = data.frame(effectName = nm,
                                       stringsAsFactors = FALSE),
                 theta    = stats::rnorm(n_eff, 0.3, 0.5),
                 covtheta = V),
            class = "sienaFit")
}

jn_mock_sbf <- function(n_iter = 400, nGroup = 3, n_rate = 2, n_eff = 5,
                        random = c(TRUE, FALSE, TRUE, FALSE, FALSE),
                        seed = 9) {
  set.seed(seed)
  effects <- data.frame(
    effectName = c(paste0("rate ", seq_len(n_rate * nGroup)),
                   c("altX", "egoX", "altX x egoX", "recip",
                     "transTrip")[seq_len(n_eff)]),
    randomEffects = c(rep(FALSE, n_rate * nGroup), random),
    stringsAsFactors = FALSE)

  n_par <- n_rate + n_eff
  structure(
    list(effects         = effects,
         basicRate       = c(rep(TRUE, n_rate * nGroup), rep(FALSE, n_eff)),
         nGroup          = nGroup,
         nwarm           = 50,
         ThinParameters  = array(stats::rnorm(n_iter * nGroup * n_par, 0.3, 0.2),
                                 dim = c(n_iter, nGroup, n_par)),
         ThinPosteriorMu = matrix(stats::rnorm(n_iter * (n_rate + sum(random)),
                                               0.4, 0.2),
                                  nrow = n_iter)),
    class = "sienaBayesFit")
}

jn_test_draws <- function(n = 600, three_way = FALSE, seed = 7) {
  set.seed(seed)
  if (!three_way)
    return(cbind(x = stats::rnorm(n, 0.5, 0.2),
                 z = stats::rnorm(n, -0.3, 0.2),
                 "x:z" = stats::rnorm(n, 0.4, 0.2)))
  cbind(x = stats::rnorm(n, 0.5, 0.2),
        z = stats::rnorm(n, -0.3, 0.2),
        w = stats::rnorm(n, 0.2, 0.2),
        "x:z" = stats::rnorm(n, 0.4, 0.2),
        "x:w" = stats::rnorm(n, 0.1, 0.2),
        "z:w" = stats::rnorm(n, -0.1, 0.2),
        "x:z:w" = stats::rnorm(n, 0.2, 0.2))
}
