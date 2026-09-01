test_that("JNK_freq() warns once and returns the 1.0.x two-way layout", {
  d <- jn_test_data()
  expect_warning(
    old <- JNK_freq(lm(y ~ x * z, data = d), theta_1 = "x", theta_2 = "z",
                    range_size = 20),
    "deprecated")

  expect_named(old, c("param_table", "plots"))
  expect_named(old$param_table, c("x", "z"))
  expect_named(old$plots, c("x", "z"))
  expect_true(all(c("theta", "moderator", "mod_vals", "theta_vals",
                    "theta_se", "theta_p", "sig") %in%
                    names(old$param_table$x)))
  expect_s3_class(old$plots$x, "ggplot")
  expect_false(inherits(old$plots$x, "patchwork"))
})


test_that("the deprecated wrapper agrees with JN() on the numbers", {
  d   <- jn_test_data()
  fit <- lm(y ~ x * z, data = d)
  suppressWarnings(
    old <- JNK_freq(fit, theta_1 = "x", theta_2 = "z", range_size = 20))
  new <- JN(fit, theta_1 = "x", theta_2 = "z", range_size = 20)

  expect_equal(old$param_table$x$theta_vals, round(new$tables$x$estimate, 3))
  expect_equal(old$param_table$x$theta_se,   round(new$tables$x$std.error, 3))
  expect_equal(old$param_table$x$sig,        new$tables$x$significant)
})


test_that("JNK_freq() returns the 1.0.x three-way layout of matrices", {
  d <- jn_test_data()
  suppressWarnings(
    old <- JNK_freq(lm(y ~ x * z * w, data = d), theta_1 = "x",
                    theta_2 = "z", theta_3 = "w", range_size = 8))

  expect_named(old, c("thetas", "standard_errors", "p_values",
                      "significance", "plots"))
  expect_true(is.matrix(old$thetas$x))
  expect_equal(dim(old$thetas$x), c(8L, 8L))

  new <- JN(lm(y ~ x * z * w, data = d), theta_1 = "x", theta_2 = "z",
            theta_3 = "w", range_size = 8)
  # the matrix is the long grid folded back, row = mod1, column = mod2
  tab <- new$tables$x
  expect_equal(old$thetas$x[1, 1],
               tab$estimate[tab$mod1_value == min(tab$mod1_value) &
                            tab$mod2_value == min(tab$mod2_value)])
  expect_equal(as.numeric(rownames(old$thetas$x)),
               sort(unique(tab$mod1_value)))
})


test_that("JNK_bayes() returns the 1.0.x layout for both widths", {
  suppressWarnings(
    old2 <- JNK_bayes(jn_test_draws(), theta_1 = "x", theta_2 = "z",
                      theta_1_vals = seq(-2, 2, 1),
                      theta_2_vals = seq(-2, 2, 1)))
  expect_named(old2, c("x_table", "z_table", "x_plot", "z_plot"))
  expect_true(all(c("thetaPostMean", "thetaPostSD", "bayes_p",
                    "thetaPost2.5", "thetaPost97.5") %in%
                    names(old2$x_table)))
  expect_s3_class(old2$x_plot, "ggplot")

  suppressWarnings(
    old3 <- JNK_bayes(jn_test_draws(three_way = TRUE),
                      theta_1 = "x", theta_2 = "z", theta_3 = "w",
                      theta_1_vals = seq(-1, 1, 1),
                      theta_2_vals = seq(-1, 1, 1),
                      theta_3_vals = seq(-1, 1, 1),
                      thresholds = c(0.05, 0.95)))
  expect_named(old3, c("result_tables", "post_mean_plots", "bayes_p_plots"))
  expect_named(old3$result_tables, c("x", "z", "w"))
  expect_true("pattern" %in% names(old3$result_tables$x))
})


test_that("JNK_freq() still accepts a bare coefficient vector and covariance", {
  d   <- jn_test_data()
  fit <- lm(y ~ x * z, data = d)
  idx <- c("x", "z", "x:z")

  suppressWarnings(
    old <- JNK_freq(covar = vcov(fit)[idx, idx], coefs = coef(fit)[idx],
                    name = c("x", "z"), theta_1 = "x", theta_2 = "z",
                    theta_1_vals = c(-3, 3), theta_2_vals = c(-3, 3),
                    range_size = 20))
  expect_named(old, c("param_table", "plots"))
  expect_equal(nrow(old$param_table$x), 20)
})
