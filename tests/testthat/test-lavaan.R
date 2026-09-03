# A structural model has several endogenous variables, so unlike a regression
# there is no single interaction to analyse. Most of what is worth testing here
# is how that ambiguity is resolved and what happens when it cannot be.

jn_lav_data <- function(n = 300, seed = 1402) {
  set.seed(seed)
  d <- data.frame(x = stats::rnorm(n), z = stats::rnorm(n),
                  w = stats::rnorm(n))
  d$y1 <- d$x + 0.5 * d$z + 0.6 * d$x * d$z + stats::rnorm(n)
  d$y2 <- -d$x + 0.3 * d$z - 0.4 * d$x * d$z + stats::rnorm(n)
  d
}

jn_lav_two_outcomes <- function(d = jn_lav_data()) {
  lavaan::sem("y1 ~ x + z + x:z\n y2 ~ x + z + x:z", data = d)
}


test_that("a lavaan fit is addressed by name and analysed with Wald tests", {
  skip_if_not_installed("lavaan")
  d   <- jn_lav_data()
  fit <- lavaan::sem("y1 ~ x + z + x:z", data = d)

  inp <- jn_input(fit, theta_1 = "x", theta_2 = "z")
  expect_s3_class(inp, "jn_wald")
  expect_equal(inp$labels, c("x", "z"))

  res <- JN(fit, theta_1 = "x", theta_2 = "z")
  expect_s3_class(res, "JN_wald")
  expect_equal(res$inference, "wald")
})


test_that("the estimates match what lm() gives for the same equation", {
  skip_if_not_installed("lavaan")
  d <- jn_lav_data()

  lav <- JN(lavaan::sem("y1 ~ x + z + x:z", data = d),
            theta_1 = "x", theta_2 = "z")
  ols <- JN(lm(y1 ~ x * z, data = d), theta_1 = "x", theta_2 = "z")

  # identical point estimates; the covariance differs only by lavaan's ML
  # divisor of n against OLS's n - k
  expect_equal(unname(coef(lav)), unname(coef(ols)), tolerance = 1e-8)
  expect_equal(unname(vcov(lav)), unname(vcov(ols)), tolerance = 1e-2)
  expect_equal(jn_regions(lav)$focal, jn_regions(ols)$focal)
})


test_that("several outcomes without `outcome` stop and name the candidates", {
  skip_if_not_installed("lavaan")
  expect_error(JN(jn_lav_two_outcomes(), theta_1 = "x", theta_2 = "z"),
               "2 outcomes are regressed on")
  expect_error(JN(jn_lav_two_outcomes(), theta_1 = "x", theta_2 = "z"),
               "outcome")
})


test_that("`outcome` picks the equation, and the equations differ", {
  skip_if_not_installed("lavaan")
  fit <- jn_lav_two_outcomes()

  j1 <- JN(fit, theta_1 = "x", theta_2 = "z", outcome = "y1")
  j2 <- JN(fit, theta_1 = "x", theta_2 = "z", outcome = "y2")

  # y1 was generated with a positive interaction and y2 with a negative one,
  # so picking the wrong outcome is not a harmless mistake
  expect_gt(coef(j1)[3], 0)
  expect_lt(coef(j2)[3], 0)
  expect_false(isTRUE(all.equal(coef(j1), coef(j2))))
})


test_that("a single candidate outcome needs no `outcome` argument", {
  skip_if_not_installed("lavaan")
  d <- jn_lav_data()
  # y2 is an outcome too, but only y1 is regressed on both focal variables
  fit <- lavaan::sem("y1 ~ x + z + x:z\n y2 ~ x", data = d)

  expect_equal(unname(coef(JN(fit, theta_1 = "x", theta_2 = "z"))),
               unname(coef(JN(fit, theta_1 = "x", theta_2 = "z",
                              outcome = "y1"))))
})


test_that("an unusable `outcome` is refused with the available ones listed", {
  skip_if_not_installed("lavaan")
  fit <- jn_lav_two_outcomes()

  expect_error(JN(fit, theta_1 = "x", theta_2 = "z", outcome = "nope"),
               "Nothing is regressed on")
  expect_error(JN(fit, theta_1 = "x", theta_2 = "z",
                  outcome = c("y1", "y2")),
               "single dependent variable")
})


test_that("an outcome missing a main effect says so, not the interaction", {
  skip_if_not_installed("lavaan")
  d   <- jn_lav_data()
  fit <- lavaan::sem("y1 ~ x + z + x:z\n y2 ~ x", data = d)

  expect_error(JN(fit, theta_1 = "x", theta_2 = "z", outcome = "y2"),
               "not regressed on 'z'")
})


test_that("three-way products are addressed by name, since lavaan needs them
           built by hand", {
  skip_if_not_installed("lavaan")
  d <- jn_lav_data()
  d$xz <- d$x * d$z; d$xw <- d$x * d$w
  d$zw <- d$z * d$w; d$xzw <- d$x * d$z * d$w
  d$y3 <- d$x + 0.5 * d$z - 0.3 * d$w + 0.6 * d$xzw + stats::rnorm(nrow(d))

  fit <- lavaan::sem("y3 ~ x + z + w + xz + xw + zw + xzw", data = d)
  res <- JN(fit, theta_1 = "x", theta_2 = "z", theta_3 = "w",
            theta_int_12 = "xz", theta_int_13 = "xw",
            theta_int_23 = "zw", theta_int_123 = "xzw",
            range_size = 10)

  expect_s3_class(res, "JN_3way")
  expect_length(coef(res), 7L)
  expect_equal(unname(coef(res)),
               unname(coef(JN(lm(y3 ~ x * z * w, data = d),
                              theta_1 = "x", theta_2 = "z", theta_3 = "w",
                              range_size = 10))),
               tolerance = 1e-6)
})


test_that("a product that cannot be found names the argument that supplies it", {
  skip_if_not_installed("lavaan")
  d <- jn_lav_data()
  d$xz <- d$x * d$z
  fit <- lavaan::sem("y1 ~ x + z + xz", data = d)

  expect_error(JN(fit, theta_1 = "x", theta_2 = "z"), "theta_int_12")
  expect_no_error(JN(fit, theta_1 = "x", theta_2 = "z",
                     theta_int_12 = "xz"))
})


test_that("observed values are picked up for the density panels", {
  skip_if_not_installed("lavaan")
  d   <- jn_lav_data()
  fit <- lavaan::sem("y1 ~ x + z + x:z", data = d)

  inp <- jn_input(fit, theta_1 = "x", theta_2 = "z")
  expect_length(inp$support, 2L)
  expect_equal(vapply(inp$support, length, integer(1)),
               c(nrow(d), nrow(d)))
  expect_false(is.na(jn_regions(JN(fit, theta_1 = "x",
                                   theta_2 = "z"))$data_share[1]))
})


test_that("a fit carrying no raw data still works once ranges are supplied", {
  skip_if_not_installed("lavaan")
  d  <- jn_lav_data()
  d$xz <- d$x * d$z
  v  <- c("y1", "x", "z", "xz")
  fit <- lavaan::sem("y1 ~ x + z + xz",
                     sample.cov = stats::cov(d[, v]),
                     sample.mean = colMeans(d[, v]),
                     sample.nobs = nrow(d))

  inp <- jn_input(fit, theta_1 = "x", theta_2 = "z", theta_int_12 = "xz")
  expect_null(inp$support)

  res <- JN(fit, theta_1 = "x", theta_2 = "z", theta_int_12 = "xz",
            theta_1_vals = c(-3, 3), theta_2_vals = c(-3, 3))
  expect_s3_class(res, "JN_wald")
  expect_true(is.na(jn_regions(res)$data_share[1]))
})


test_that("a latent outcome is analysed like any other", {
  skip_if_not_installed("lavaan")
  d <- jn_lav_data()
  d$i1 <- d$y1 + stats::rnorm(nrow(d))
  d$i2 <- d$y1 + stats::rnorm(nrow(d))
  d$i3 <- d$y1 + stats::rnorm(nrow(d))

  fit <- lavaan::sem("F =~ i1 + i2 + i3\n F ~ x + z + x:z", data = d)
  res <- JN(fit, theta_1 = "x", theta_2 = "z")

  expect_s3_class(res, "JN_wald")
  # the moderators are observed even though the outcome is not
  expect_length(jn_input(fit, theta_1 = "x", theta_2 = "z")$support, 2L)
})


test_that("multiple-group fits are refused rather than silently pooled", {
  skip_if_not_installed("lavaan")
  d <- jn_lav_data()
  d$g <- rep(c("a", "b"), length.out = nrow(d))
  fit <- lavaan::sem("y1 ~ x + z + x:z", data = d, group = "g")

  expect_error(JN(fit, theta_1 = "x", theta_2 = "z"),
               "Multiple-group")
})


test_that("a model without regression paths is refused", {
  skip_if_not_installed("lavaan")
  d <- jn_lav_data()
  d$i1 <- d$y1 + stats::rnorm(nrow(d))
  d$i2 <- d$y1 + stats::rnorm(nrow(d))
  d$i3 <- d$y1 + stats::rnorm(nrow(d))

  expect_error(JN(lavaan::sem("F =~ i1 + i2 + i3", data = d),
                  theta_1 = "x", theta_2 = "z"),
               "no regression paths")
})
