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
           theta_1_vals = c(-3, 3), theta_2_vals = c(-3, 3), range_size = 11)
  expect_length(jn$grid$x, 11)
  expect_equal(range(jn$grid$x), c(-3, 3))

  jn2 <- JN(lm(y ~ x * z, data = d), theta_1 = "x", theta_2 = "z",
            theta_1_vals = c(-1, 0, 1), theta_2_vals = c(-1, 0, 1))
  expect_equal(jn2$grid$x, c(-1, 0, 1))
})


test_that("false discovery rate control only ever removes significance", {
  d   <- jn_test_data()
  fit <- lm(y ~ x * z, data = d)
  plain <- JN(fit, theta_1 = "x", theta_2 = "z", range_size = 200)
  fdr   <- JN(fit, theta_1 = "x", theta_2 = "z", range_size = 200,
              control_fdr = TRUE)
  expect_true(all(fdr$tables$x$significant <= plain$tables$x$significant))
})
