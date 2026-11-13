test_that("the two-way conditional effect and its standard error follow the delta method", {
  d   <- jn_test_data()
  fit <- lm(y ~ x * z, data = d)
  idx <- c("x", "z", "x:z")
  b   <- coef(fit)[idx]
  V   <- vcov(fit)[idx, idx]

  jn  <- JN(fit, theta_1 = "x", theta_2 = "z", theta_1_vals = c(-2, 2),
            theta_2_vals = c(-2, 2), range_size = 5)
  tab <- jn$tables$x
  m   <- tab$mod_value

  # theta_x(m) = b_x + b_xz m
  expect_equal(tab$estimate, unname(b["x"] + b["x:z"] * m))

  # se from a' V a with a = (1, 0, m), the general delta-method form
  se_ref <- vapply(m, function(mm) {
    a <- c(1, 0, mm)
    sqrt(drop(t(a) %*% V %*% a))
  }, numeric(1))
  expect_equal(tab$std.error, se_ref)

  # and specifically: the covariance term is with the interaction, not with
  # the other main effect (the bug present in int3ract 1.0.x)
  wrong <- sqrt(V[1, 1] + 2 * m * V[1, 2] + m^2 * V[3, 3])
  expect_false(isTRUE(all.equal(tab$std.error, wrong)))
})


test_that("the three-way conditional effect and standard error follow the delta method", {
  d   <- jn_test_data()
  fit <- lm(y ~ x * z * w, data = d)
  idx <- c("x", "z", "w", "x:z", "x:w", "z:w", "x:z:w")
  b   <- coef(fit)[idx]
  V   <- vcov(fit)[idx, idx]

  jn  <- JN(fit, theta_1 = "x", theta_2 = "z", theta_3 = "w",
            theta_1_vals = c(-1, 1), theta_2_vals = c(-1, 1),
            theta_3_vals = c(-1, 1), range_size = 3)
  tab <- jn$tables$x

  # theta_x(m, w) = b_x + b_xz m + b_xw w + b_xzw m w, so a = d theta / d b
  ref <- t(apply(tab[, c("mod1_value", "mod2_value")], 1, function(r) {
    a <- c(1, 0, 0, r[[1]], r[[2]], 0, r[[1]] * r[[2]])
    c(drop(a %*% b), sqrt(drop(t(a) %*% V %*% a)))
  }))
  expect_equal(tab$estimate,  ref[, 1])
  expect_equal(tab$std.error, ref[, 2])
})


test_that("the posterior path summarizes the conditional draws it says it does", {
  draws <- jn_test_draws()
  vals  <- seq(-2, 2, 1)
  jn    <- JN(draws, theta_1 = "x", theta_2 = "z",
              theta_1_vals = vals, theta_2_vals = vals)
  tab   <- jn$tables$x

  # theta_x^(s)(m) = b_x^(s) + b_xz^(s) m, one conditional posterior per m
  for (i in seq_along(vals)) {
    cond <- draws[, "x"] + draws[, "x:z"] * vals[i]
    expect_equal(tab$post_mean[i], mean(cond))
    expect_equal(tab$post_sd[i],   sd(cond))
    expect_equal(tab$bayes_p[i],   mean(cond > 0))
  }
})


test_that("a length-two range is expanded and any other length is used verbatim", {
  d  <- jn_test_data()
  jn <- JN(lm(y ~ x * z, data = d), theta_1 = "x", theta_2 = "z",
           theta_1_vals = c(-3, 3), theta_2_vals = c(-3, 3), range_size = 11,
           extend = 0)
  expect_length(jn$grid$x, 11)
  expect_equal(range(jn$grid$x), c(-3, 3))

  jn2 <- JN(lm(y ~ x * z, data = d), theta_1 = "x", theta_2 = "z",
            theta_1_vals = c(-1, 0, 1), theta_2_vals = c(-1, 0, 1))
  expect_equal(jn2$grid$x, c(-1, 0, 1))
})


test_that("ranges are widened by `extend`, and the margin is flagged", {
  d  <- jn_test_data()
  jn <- JN(lm(y ~ x * z, data = d), theta_1 = "x", theta_2 = "z",
           range_size = 101)
  r  <- range(d$z)
  expect_equal(range(jn$grid$z), r + c(-0.1, 0.1) * diff(r))
  expect_equal(jn$observed_range$z, r)
  tab <- jn$tables$x
  expect_equal(tab$extrapolated, tab$mod_value < r[1] | tab$mod_value > r[2])
  expect_true(any(tab$extrapolated) && !all(tab$extrapolated))

  # supplied ranges are widened too; supplied values are not
  jn2 <- JN(lm(y ~ x * z, data = d), theta_1 = "x", theta_2 = "z",
            theta_1_vals = c(-2, 2), theta_2_vals = c(-1, 0, 1),
            extend = 0.25)
  expect_equal(range(jn2$grid$x), c(-3, 3))
  expect_equal(jn2$grid$z, c(-1, 0, 1))

  # extend = 0 gives the old grid exactly
  jn0 <- JN(lm(y ~ x * z, data = d), theta_1 = "x", theta_2 = "z",
            range_size = 101, extend = 0)
  expect_equal(range(jn0$grid$z), r)
  expect_false(any(jn0$tables$x$extrapolated))

  expect_error(JN(lm(y ~ x * z, data = d), theta_1 = "x", theta_2 = "z",
                  extend = -1), "non-negative")
})


test_that("binary variables are detected, evaluated at their levels only", {
  d   <- jn_test_data()
  d$b <- rep(c(0, 1), length.out = nrow(d))
  jn  <- JN(lm(y ~ x * b, data = d), theta_1 = "x", theta_2 = "b")
  expect_equal(unname(jn$binary), c(FALSE, TRUE))
  expect_equal(jn$grid$b, c(0, 1))
  expect_equal(nrow(jn$tables$x), 2L)
  expect_false(any(jn$tables$x$extrapolated))

  # without observed values a range c(0, 1) is a range unless declared binary
  draws <- jn_test_draws()
  jr <- JN(draws, theta_1 = "x", theta_2 = "z",
           theta_1_vals = c(-1, 1), theta_2_vals = c(0, 1))
  expect_gt(length(jr$grid$z), 2L)
  jb <- JN(draws, theta_1 = "x", theta_2 = "z",
           theta_1_vals = c(-1, 1), theta_2_vals = c(0, 1), binary = "z")
  expect_equal(jb$grid$z, c(0, 1))

  expect_error(JN(draws, theta_1 = "x", theta_2 = "z",
                  theta_1_vals = c(-1, 1), theta_2_vals = c(0, 1),
                  binary = "nope"), "not in the interaction")
})


test_that("regions along a binary moderator are its significant levels", {
  set.seed(5)
  d <- data.frame(x = rnorm(300), b = rbinom(300, 1, 0.5))
  d$y <- 0.05 * d$x + 0.8 * d$x * d$b + rnorm(300)
  jn  <- JN(lm(y ~ x * b, data = d), theta_1 = "x", theta_2 = "b")
  reg <- jn_regions(jn)
  rx  <- reg[reg$focal == "x", ]
  expect_equal(rx$from, 1)
  expect_equal(rx$to, 1)
  expect_equal(rx$sign, "+")
  expect_equal(rx$data_share, mean(d$b == 1))
  expect_output(print(jn), "b = 1.000")
})


test_that("false discovery rate control only ever removes significance", {
  d   <- jn_test_data()
  fit <- lm(y ~ x * z, data = d)
  plain <- JN(fit, theta_1 = "x", theta_2 = "z", range_size = 200)
  fdr   <- JN(fit, theta_1 = "x", theta_2 = "z", range_size = 200,
              control_fdr = TRUE)
  expect_true(all(fdr$tables$x$significant <= plain$tables$x$significant))
})
