# RSiena and multiSiena objects cannot be constructed cheaply, so these tests
# work against stand-ins carrying the structure the extractors rely on.

test_that("sienaFit input is addressed by position and analysed with Wald tests", {
  fit <- jn_mock_sienafit()

  inp <- jn_input(fit, theta_1 = 1, theta_2 = 2, theta_int_12 = 3)
  expect_s3_class(inp, "jn_wald")
  expect_equal(inp$labels, c("altX", "egoX"))

  res <- JN(fit, theta_1 = 1, theta_2 = 2, theta_int_12 = 3,
            theta_1_vals = c(-2, 2), theta_2_vals = c(-2, 2),
            range_size = 20)
  expect_s3_class(res, "JN_wald")
  expect_equal(res$inference, "wald")
  expect_equal(coef(res), stats::setNames(fit$theta[1:3],
                                          c("altX", "egoX", "altX:egoX")))
})


test_that("a sienaFit position out of range names the available effects", {
  fit <- jn_mock_sienafit()
  expect_error(
    JN(fit, theta_1 = 1, theta_2 = 99, theta_int_12 = 3,
       theta_1_vals = c(-2, 2), theta_2_vals = c(-2, 2)),
    "out of range")
})


test_that("sienaBayesFit input dispatches to the posterior path", {
  x   <- jn_mock_sbf()
  res <- JN(x, theta_1 = 1, theta_2 = 2, theta_int_12 = 3,
            theta_1_vals = seq(-2, 2, 1), theta_2_vals = seq(-2, 2, 1))

  expect_s3_class(res, "JN_posterior")
  expect_equal(res$inference, "posterior")
  expect_equal(res$labels, c("altX", "egoX"))
  # not a Wald analysis, whatever else happens
  expect_false(inherits(res, "JN_wald"))
  expect_output(print(res), "conditional posteriors")
})


test_that("multiSiena objects reach the same method", {
  x <- jn_mock_sbf()
  class(x) <- "multiSiena"
  expect_s3_class(jn_input(x, theta_1 = 1, theta_2 = 2, theta_int_12 = 3),
                  "jn_posterior")
})


test_that("the level the draws came from is detected and reported", {
  # position 1 and 3 vary between groups -> the hyper-parameter mu is used
  random <- JN(jn_mock_sbf(), theta_1 = 1, theta_2 = 2, theta_int_12 = 3,
               theta_1_vals = seq(-2, 2, 1), theta_2_vals = seq(-2, 2, 1))
  expect_equal(random$group, "Mu")

  # nothing varies between groups -> the shared parameter eta is used
  fixed <- JN(jn_mock_sbf(random = rep(FALSE, 5)),
              theta_1 = 1, theta_2 = 2, theta_int_12 = 3,
              theta_1_vals = seq(-2, 2, 1), theta_2_vals = seq(-2, 2, 1))
  expect_equal(fixed$group, "Eta")

  # and it shows up where the user can see it
  expect_output(print(random), "Mu")
})


test_that("hyper_only = FALSE adds one analysis per group", {
  x   <- jn_mock_sbf()
  res <- JN(x, theta_1 = 1, theta_2 = 2, theta_int_12 = 3,
            theta_1_vals = seq(-2, 2, 1), theta_2_vals = seq(-2, 2, 1),
            hyper_only = FALSE)

  expect_s3_class(res, "JN_list")
  expect_named(res, c("hyper", "group1", "group2", "group3"))
  expect_equal(res$hyper$group, "Mu")
  expect_equal(res$group2$group, "group2")
  # the group analyses are genuinely different draws, not copies
  expect_false(isTRUE(all.equal(res$group1$tables[[1]]$post_mean,
                                res$group2$tables[[1]]$post_mean)))
  # and the JN_list methods work on the result
  expect_output(print(res), "further group-level analyses")
  expect_length(summary(res), 4)
  expect_s3_class(as.data.frame(res), "data.frame")
})


test_that("no group-level analyses are produced when nothing varies by group", {
  x <- jn_mock_sbf(random = rep(FALSE, 5))
  expect_message(
    res <- JN(x, theta_1 = 1, theta_2 = 2, theta_int_12 = 3,
              theta_1_vals = seq(-2, 2, 1), theta_2_vals = seq(-2, 2, 1),
              hyper_only = FALSE),
    "does not vary|varies between groups")
  # a JN_list of identical pictures would be worse than nothing
  expect_s3_class(res, "JN_posterior")
  expect_false(inherits(res, "JN_list"))
})


test_that("three-way sienaBayesFit analyses work", {
  x   <- jn_mock_sbf()
  res <- JN(x, theta_1 = 1, theta_2 = 2, theta_3 = 4,
            theta_int_12 = 3, theta_int_13 = 5,
            theta_int_23 = 3, theta_int_123 = 5,
            theta_1_vals = seq(-2, 2, 1), theta_2_vals = seq(-2, 2, 1),
            theta_3_vals = seq(-2, 2, 1))
  expect_s3_class(res, "JN_3way")
  expect_s3_class(res, "JN_posterior")
  expect_length(res$labels, 3)
})


test_that("burn_in and thin reach the draws", {
  x    <- jn_mock_sbf(n_iter = 400)
  full <- jn_input(x, theta_1 = 1, theta_2 = 2, theta_int_12 = 3)
  thin <- jn_input(x, theta_1 = 1, theta_2 = 2, theta_int_12 = 3,
                   burn_in = 100, thin = 4)
  expect_lt(nrow(thin$draws), nrow(full$draws))
})


test_that("the deprecated JNK_bayes layout still comes back for multiSiena", {
  x <- jn_mock_sbf()
  suppressWarnings(
    old <- JNK_bayes(x, theta_1 = 1, theta_2 = 2, theta_int_12 = 3,
                     theta_1_vals = seq(-2, 2, 1),
                     theta_2_vals = seq(-2, 2, 1),
                     hyper_only = FALSE))
  expect_named(old, c("Mu", "random_groups_effects"))
  expect_named(old$random_groups_effects, c("group1", "group2", "group3"))
  expect_true(all(c("altX_table", "egoX_table") %in% names(old$Mu)))
})


test_that("a structurally inconsistent sienaBayesFit is rejected, not guessed at", {
  x <- jn_mock_sbf()
  x$ThinParameters <- x$ThinParameters[, , -1, drop = FALSE]
  expect_error(
    JN(x, theta_1 = 1, theta_2 = 2, theta_int_12 = 3,
       theta_1_vals = seq(-2, 2, 1), theta_2_vals = seq(-2, 2, 1)),
    "ThinParameters")
})
