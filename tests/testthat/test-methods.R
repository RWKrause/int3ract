test_that("print() and summary() describe the analysis they ran", {
  d  <- jn_test_data()
  jn <- JN(lm(y ~ x * z, data = d), theta_1 = "x", theta_2 = "z",
           range_size = 50)

  out <- capture.output(print(jn))
  expect_true(any(grepl("two-way", out)))
  expect_true(any(grepl("Wald z tests", out)))
  expect_true(any(grepl("Variables: x, z", out)))

  s <- summary(jn)
  expect_s3_class(s, "summary.JN")
  expect_s3_class(s$regions, "data.frame")
})


test_that("summary() reports the range each variable was evaluated over", {
  d  <- jn_test_data()
  jn <- JN(lm(y ~ x * z, data = d), theta_1 = "x", theta_2 = "z",
           range_size = 40)
  s  <- summary(jn)

  expect_s3_class(s$ranges, "data.frame")
  expect_equal(s$ranges$variable, c("x", "z"))
  expect_equal(s$ranges$n, c(40L, 40L))
  # the reported range is the grid actually used, not the raw data range
  expect_equal(s$ranges$from, vapply(jn$grid, min, numeric(1)),
               ignore_attr = TRUE)
  expect_equal(s$ranges$to, vapply(jn$grid, max, numeric(1)),
               ignore_attr = TRUE)

  out <- capture.output(print(s))
  expect_true(any(grepl("Evaluated over", out)))
  expect_true(any(grepl("\\[.*,.*\\]", out)))
})


test_that("the evaluated range follows explicitly supplied moderator values", {
  jn <- JN(jn_test_draws(), theta_1 = "x", theta_2 = "z",
           theta_1_vals = seq(-3, 3, 0.5), theta_2_vals = c(-2, 2))
  s  <- summary(jn)
  expect_equal(s$ranges$from, c(-3, -2))
  expect_equal(s$ranges$to,   c(3, 2))
  expect_equal(s$ranges$n[1], 13L)     # used verbatim
})


test_that("three-way summaries report all three ranges", {
  d  <- jn_test_data()
  s  <- summary(JN(lm(y ~ x * z * w, data = d), theta_1 = "x",
                   theta_2 = "z", theta_3 = "w", range_size = 12))
  expect_equal(nrow(s$ranges), 3L)
  expect_equal(s$ranges$variable, c("x", "z", "w"))
  expect_true(all(s$ranges$n == 12L))
})


test_that("posterior analyses are described as such", {
  jn <- JN(jn_test_draws(), theta_1 = "x", theta_2 = "z",
           theta_1_vals = seq(-2, 2, 1), theta_2_vals = seq(-2, 2, 1))
  out <- capture.output(print(jn))
  expect_true(any(grepl("conditional posteriors", out)))
  expect_false(any(grepl("Wald", out)))
})


test_that("as.data.frame() returns the full grid and coef()/vcov() the parameters", {
  d  <- jn_test_data()
  jn <- JN(lm(y ~ x * z, data = d), theta_1 = "x", theta_2 = "z",
           range_size = 25)

  df <- as.data.frame(jn)
  expect_s3_class(df, "data.frame")
  expect_equal(nrow(df), 50)              # 25 grid points for each of x and z
  expect_setequal(unique(df$focal), c("x", "z"))

  expect_named(coef(jn), c("x", "z", "x:z"))
  expect_equal(dim(vcov(jn)), c(3L, 3L))
  expect_equal(colnames(vcov(jn)), c("x", "z", "x:z"))

  jn3 <- JN(lm(y ~ x * z * w, data = d), theta_1 = "x", theta_2 = "z",
            theta_3 = "w", range_size = 5)
  expect_named(coef(jn3),
               c("x", "z", "w", "x:z", "x:w", "z:w", "x:z:w"))
  expect_equal(nrow(as.data.frame(jn3)), 3 * 25)
})


test_that("the default Bayesian thresholds follow alpha", {
  draws <- jn_test_draws()
  jn <- JN(draws, theta_1 = "x", theta_2 = "z",
           theta_1_vals = seq(-2, 2, 1), theta_2_vals = seq(-2, 2, 1))
  expect_equal(jn$thresholds, c(0.025, 0.975))

  jn10 <- JN(draws, theta_1 = "x", theta_2 = "z", alpha = 0.10,
             theta_1_vals = seq(-2, 2, 1), theta_2_vals = seq(-2, 2, 1))
  expect_equal(jn10$thresholds, c(0.05, 0.95))
})


test_that("moderator ranges must be supplied when the input carries no data", {
  expect_error(JN(jn_test_draws(), theta_1 = "x", theta_2 = "z"),
               "supply the moderator ranges")
})


test_that("grouped lme4 analyses return a JN_list with the same methods", {
  skip_if_not_installed("lme4")
  set.seed(11)
  n  <- 400
  gp <- factor(rep(seq_len(20), each = n / 20))
  d  <- data.frame(g = gp, x = rnorm(n), z = rnorm(n))
  d$y <- d$x + 0.4 * d$x * d$z + rnorm(20)[as.integer(gp)] + rnorm(n)

  fit <- lme4::lmer(y ~ x * z + (1 | g), data = d)
  res <- suppressWarnings(
    JN(fit, theta_1 = "x", theta_2 = "z", range_size = 20,
       fixed_only = FALSE, group_var = "g"))

  expect_s3_class(res, "JN_list")
  expect_equal(names(res)[1], "fixed")
  expect_length(res, 21)
  expect_output(print(res), "further group-level analyses")
  expect_s3_class(as.data.frame(res), "data.frame")
  expect_length(summary(res), 21)
})
