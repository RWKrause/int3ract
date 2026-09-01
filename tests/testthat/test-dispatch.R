test_that("JN() dispatches to the right carrier for each supported class", {
  d <- jn_test_data()

  expect_s3_class(jn_input(lm(y ~ x * z, data = d),
                           theta_1 = "x", theta_2 = "z"), "jn_wald")
  expect_s3_class(jn_input(glm(I(y > 0) ~ x * z, family = binomial, data = d),
                           theta_1 = "x", theta_2 = "z"), "jn_wald")
  expect_s3_class(jn_input(jn_test_draws(),
                           theta_1 = "x", theta_2 = "z"), "jn_posterior")

  expect_equal(JN(lm(y ~ x * z, data = d),
                  theta_1 = "x", theta_2 = "z")$inference, "wald")
  expect_equal(JN(jn_test_draws(), theta_1 = "x", theta_2 = "z",
                  theta_1_vals = seq(-2, 2, 1),
                  theta_2_vals = seq(-2, 2, 1))$inference, "posterior")
})


test_that("an unsupported class is refused with a message pointing at the fix", {
  expect_error(JN(structure(list(), class = "nonsense"),
                  theta_1 = "x", theta_2 = "z"),
               "jn_input")
})


test_that("a user-supplied jn_input method is enough to support a new class", {
  # exactly the case the JSS editor raised: a model fitted by variational
  # inference, whose variational parameters approximate the posterior
  d   <- jn_test_data()
  fit <- lm(y ~ x * z, data = d)
  idx <- c("x", "z", "x:z")
  vb  <- structure(list(mu = coef(fit)[idx], Sigma = vcov(fit)[idx, idx]),
                   class = "fakevbfit")

  jn_input.fakevbfit <- function(object, theta_1, theta_2, theta_3 = NULL, ...)
    jn_wald(coefficients = object$mu, vcov = object$Sigma,
            labels = c(theta_1, theta_2),
            ranges = list(c(-3, 3), c(-3, 3)))
  registerS3method("jn_input", "fakevbfit", jn_input.fakevbfit)

  res <- JN(vb, theta_1 = "x", theta_2 = "z", range_size = 50)

  expect_s3_class(res, "JN")
  expect_equal(res$inference, "wald")
  # and every method on the class works on it, with no further code
  expect_output(print(res), "Regions of significance")
  expect_output(print(summary(res)), "Regions of significance")
  expect_s3_class(plot(res, which = 1), "ggplot")
  expect_s3_class(as.data.frame(res), "data.frame")

  # the answer matches what the equivalent lm gives over the same grid
  ref <- JN(fit, theta_1 = "x", theta_2 = "z",
            theta_1_vals = c(-3, 3), theta_2_vals = c(-3, 3),
            range_size = 50)
  expect_equal(res$tables$x$estimate,  ref$tables$x$estimate)
  expect_equal(res$tables$x$std.error, ref$tables$x$std.error)
})


test_that("a hand-built carrier can be handed straight to JN()", {
  d   <- jn_test_data()
  fit <- lm(y ~ x * z, data = d)
  idx <- c("x", "z", "x:z")

  inp <- jn_wald(coefficients = coef(fit)[idx], vcov = vcov(fit)[idx, idx],
                 labels = c("x", "z"),
                 ranges = list(range(d$x), range(d$z)),
                 support = list(d$x, d$z))
  expect_s3_class(inp, "jn_input")
  expect_output(print(inp), "observed values available")
  expect_s3_class(JN(inp, range_size = 20), "JN")
})


test_that("the carriers reject malformed input", {
  expect_error(jn_wald(1:2, diag(2), labels = c("a", "b")), "3 elements")
  expect_error(jn_wald(1:3, diag(2), labels = c("a", "b")), "3 by 3")
  expect_error(jn_wald(1:3, diag(3), labels = c("a", "b", "c", "d")),
               "2 or 3 elements")
  expect_error(jn_posterior(matrix(rnorm(20), ncol = 2),
                            labels = c("a", "b")), "3 columns")
  expect_error(jn_wald(1:3, diag(3), labels = c("a", "b"),
                       support = list(1:5)), "length 2")
})


test_that("interaction terms are found whichever order the formula produced", {
  d <- jn_test_data()
  # z * x rather than x * z, so the coefficient is named "x:z" or "z:x"
  fit <- lm(y ~ z * x, data = d)
  expect_s3_class(JN(fit, theta_1 = "x", theta_2 = "z"), "JN")
})
