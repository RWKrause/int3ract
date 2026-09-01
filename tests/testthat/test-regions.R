test_that("the exact two-way boundaries are the roots of the JN quadratic", {
  d   <- jn_test_data(n = 300)
  fit <- lm(y ~ x * z, data = d)
  jn  <- JN(fit, theta_1 = "x", theta_2 = "z",
            theta_1_vals = c(-6, 6), theta_2_vals = c(-6, 6))
  reg <- jn_regions(jn)
  reg <- reg[reg$focal == "x" & !is.na(reg$from), , drop = FALSE]
  skip_if(nrow(reg) == 0)

  idx <- c("x", "z", "x:z")
  b   <- coef(fit)[idx]
  V   <- vcov(fit)[idx, idx]
  z2  <- qnorm(0.025)^2

  # at a boundary the test statistic sits exactly on the critical value
  for (m in c(reg$from, reg$to)) {
    if (abs(m) >= 6) next          # clipped at the edge, not a real boundary
    est <- b[["x"]] + b[["x:z"]] * m
    se  <- sqrt(V[1, 1] + 2 * m * V[1, 3] + m^2 * V[3, 3])
    expect_equal(abs(est / se), abs(qnorm(0.025)), tolerance = 1e-6)
    expect_equal(est^2 - z2 * se^2, 0, tolerance = 1e-6)
  }
})


test_that("the exact boundaries agree with the significance flags on the grid", {
  d   <- jn_test_data(n = 300)
  jn  <- JN(lm(y ~ x * z, data = d), theta_1 = "x", theta_2 = "z",
            theta_1_vals = c(-6, 6), theta_2_vals = c(-6, 6),
            range_size = 4001)
  reg <- jn_regions(jn)
  reg <- reg[reg$focal == "x" & !is.na(reg$from), , drop = FALSE]
  tab <- jn$tables$x

  inside <- Reduce(`|`, Map(function(a, b) tab$mod_value >= a &
                                           tab$mod_value <= b,
                            reg$from, reg$to),
                   init = rep(FALSE, nrow(tab)))
  # allow the two grid points straddling each boundary to disagree
  expect_lt(sum(inside != tab$significant), 2 * nrow(reg) + 1)
})


test_that("regions of significance are reported per moderator and signed", {
  d   <- jn_test_data(n = 300)
  jn  <- JN(lm(y ~ x * z, data = d), theta_1 = "x", theta_2 = "z")
  reg <- jn_regions(jn)

  expect_true(all(c("focal", "moderator", "from", "to", "sign",
                    "data_share") %in% names(reg)))
  expect_setequal(reg$focal, c("x", "z"))
  ok <- !is.na(reg$sign)
  expect_true(all(reg$sign[ok] %in% c("+", "-")))
  # the share of observations inside a region is a proportion
  sh <- reg$data_share[!is.na(reg$data_share)]
  expect_true(all(sh >= 0 & sh <= 1))
})


test_that("summary() prints the regions and flags range-bounded edges", {
  d  <- jn_test_data(n = 300)
  # a narrow evaluation range guarantees the region runs off the edge
  jn <- JN(lm(y ~ x * z, data = d), theta_1 = "x", theta_2 = "z",
           theta_1_vals = c(-0.2, 0.2), theta_2_vals = c(-0.2, 0.2))
  out <- capture.output(print(summary(jn)))
  expect_true(any(grepl("Regions of significance", out)))
  expect_true(any(grepl("\\*", out)))
})


test_that("three-way regions are reported along the second moderator", {
  d  <- jn_test_data(n = 300)
  jn <- JN(lm(y ~ x * z * w, data = d), theta_1 = "x", theta_2 = "z",
           theta_3 = "w", range_size = 30)
  reg <- jn_regions(jn)

  expect_true(all(c("focal", "mod1", "mod1_value", "moderator",
                    "from", "to") %in% names(reg)))
  # three conditioning values per focal variable by default (the quartiles)
  expect_lte(length(unique(reg$mod1_value[reg$focal == "x"])), 3)

  reg2 <- jn_regions(jn, at = 0)
  expect_length(unique(reg2$mod1_value[reg2$focal == "x"]), 1)
})


test_that("an analysis with no significant region says so", {
  set.seed(99)
  d <- data.frame(x = rnorm(60), z = rnorm(60))
  d$y <- rnorm(60, sd = 10)          # nothing to find
  jn  <- JN(lm(y ~ x * z, data = d), theta_1 = "x", theta_2 = "z")
  expect_output(print(summary(jn)), "no region of significance")
})
