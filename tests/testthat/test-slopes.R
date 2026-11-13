# Fixtures with a known shape, large enough that the regions are not in doubt.
slopes_data <- function(n = 400, seed = 11) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n), z = rnorm(n), w = rnorm(n),
                  b = rbinom(n, 1, 0.4))
  d$y_cross <- 0.2 * d$x + 0.5 * d$z + 0.8 * d$x * d$z + rnorm(n)
  d$y_all   <- 2 * d$x + 0.3 * d$x * d$z + rnorm(n)
  d$y_none  <- rnorm(n)
  d$y_bin   <- 0.5 * d$b + 0.7 * d$x * d$b + rnorm(n)
  d$y_3     <- d$x + 0.5 * d$z - 0.5 * d$w + 0.6 * d$x * d$z * d$w +
               0.4 * d$x * d$z + rnorm(n, sd = 2)
  d$y_3b    <- d$x + 0.5 * d$z + 0.5 * d$b + 0.6 * d$x * d$z * d$b +
               rnorm(n, sd = 1.5)
  d
}

slope_lines <- function(jn, i = 1L) {
  .jn_slopes_choice(jn, i, 3L - i, c(0, 0))$lines
}


test_that("a crossover interaction gives four slopes, two of each sign", {
  d  <- slopes_data()
  jn <- JN(lm(y_cross ~ x * z, data = d), theta_1 = "x", theta_2 = "z")
  ln <- slope_lines(jn)
  expect_equal(nrow(ln), 4L)
  expect_equal(sort(table(ln$category)), sort(c(neg = 2L, pos = 2L)),
               ignore_attr = TRUE)
  # the bounding slopes are the JN boundaries, within the observed data
  reg <- jn_regions(JN(lm(y_cross ~ x * z, data = d), theta_1 = "x",
                       theta_2 = "z", extend = 0))
  reg <- reg[reg$focal == "x", ]
  expect_equal(sort(ln$mod_value), sort(c(reg$from, reg$to)))
})


test_that("an effect significant throughout gives the data limits and the mode", {
  d  <- slopes_data()
  jn <- JN(lm(y_all ~ x * z, data = d), theta_1 = "x", theta_2 = "z")
  ln <- slope_lines(jn)
  expect_equal(sort(ln$category), c("mode", "pos", "pos"))
  expect_equal(sort(ln$mod_value[ln$category == "pos"]), range(d$z))
  m <- ln$mod_value[ln$category == "mode"]
  expect_true(m > min(d$z) && m < max(d$z))
})


test_that("an effect significant nowhere gives one slope at the mean", {
  d  <- slopes_data()
  jn <- JN(lm(y_none ~ x * z, data = d), theta_1 = "x", theta_2 = "z")
  ln <- slope_lines(jn)
  expect_equal(nrow(ln), 1L)
  expect_equal(ln$category, "ns")
  expect_equal(ln$mod_value, mean(d$z))
})


test_that("a binary moderator gives one slope per level, each tested", {
  d  <- slopes_data()
  jn <- JN(lm(y_bin ~ x * b, data = d), theta_1 = "x", theta_2 = "b")
  ln <- slope_lines(jn)
  expect_equal(ln$mod_value, c(0, 1))
  expect_equal(ln$category, c("ns", "pos"))
})


test_that("the drawn slopes are the model's predictions", {
  d   <- slopes_data()
  fit <- lm(y_all ~ x * z + w, data = d)
  jn  <- JN(fit, theta_1 = "x", theta_2 = "z")
  nd  <- data.frame(x = c(-1, 0.5, 2), z = c(0.7, -0.3, 1.1), w = mean(d$w))
  expect_equal(.jn_linpred(jn, cbind(nd$x, nd$z)),
               unname(predict(fit, nd)))
  expect_identical(jn$input$outcome, "y_all")
})


test_that("glm slopes are on the link scale and say so", {
  d   <- slopes_data()
  d$s <- rbinom(nrow(d), 1, plogis(d$x + d$x * d$z))
  fit <- glm(s ~ x * z, data = d, family = binomial)
  jn  <- JN(fit, theta_1 = "x", theta_2 = "z")
  expect_identical(jn$input$outcome, "s (logit scale)")
  nd  <- data.frame(x = c(-1, 1), z = c(0.5, -0.5))
  expect_equal(.jn_linpred(jn, as.matrix(nd)),
               unname(predict(fit, nd, type = "link")))
})


test_that("without an intercept the outcome axis says it is up to a constant", {
  jn <- JN(jn_test_draws(), theta_1 = "x", theta_2 = "z",
           theta_1_vals = c(-2, 2), theta_2_vals = c(-2, 2))
  expect_null(jn$input$intercept)
  expect_match(.jn_outcome_axis(jn), "up to a constant")
})


test_that("three-way slice boundaries agree with the grid-based regions", {
  d  <- slopes_data()
  jn <- JN(lm(y_3 ~ x * z * w, data = d), theta_1 = "x", theta_2 = "z",
           theta_3 = "w", range_size = 400, extend = 0)
  at  <- sort(unique(jn$tables$x$mod1_value))[200]
  reg <- jn_regions(jn, at = at)
  reg <- reg[reg$focal == "x" & !is.na(reg$from), ]
  rng <- range(jn$grid$w)
  ex  <- .jn_slice_regions(jn, 1L, 3L, c(0, at, 0), rng[1], rng[2])
  expect_equal(nrow(ex), nrow(reg))
  step <- diff(rng) / 399
  expect_true(all(abs(ex$from - reg$from) < step))
  expect_true(all(abs(ex$to   - reg$to)   < step))
  expect_equal(ex$sign, reg$sign)
})


test_that("posterior slices find their regions from the draws", {
  jn <- JN(jn_test_draws(), theta_1 = "x", theta_2 = "z",
           theta_1_vals = c(-2, 2), theta_2_vals = c(-2, 2))
  iv <- .jn_slice_regions(jn, 1L, 2L, c(0, 0), -2, 2)
  expect_true(nrow(iv) >= 1L)
  expect_true(all(iv$sign %in% c("+", "-")))
})


test_that("slopes figures build for every kind of analysis", {
  d <- slopes_data()
  builds <- function(p) {
    expect_true(inherits(p, "ggplot") || inherits(p, "patchwork"))
    expect_silent(ggplot2::ggplot_build(if (inherits(p, "patchwork"))
                                          p[[1]] else p))
  }
  jn2 <- JN(lm(y_cross ~ x * z, data = d), theta_1 = "x", theta_2 = "z")
  builds(autoplot(jn2, "x", type = "slopes"))
  builds(autoplot(jn2, "x", type = "slopes",
                  style = jn_style(slope_shading = "density")))
  builds(autoplot(JN(lm(y_bin ~ x * b, data = d), theta_1 = "x",
                     theta_2 = "b"), "b", type = "slopes"))
  builds(autoplot(JN(jn_test_draws(), theta_1 = "x", theta_2 = "z",
                     theta_1_vals = c(-2, 2), theta_2_vals = c(-2, 2)),
                  "x", type = "slopes"))

  jn3 <- JN(lm(y_3 ~ x * z * w, data = d), theta_1 = "x", theta_2 = "z",
            theta_3 = "w", range_size = 15)
  p3 <- autoplot(jn3, "x", type = "slopes")
  builds(p3)
  expect_equal(nrow(ggplot2::ggplot_build(p3)$layout$layout), 3L)
  builds(autoplot(jn3, "x", type = "slopes", at = c(-1, 1)))
  builds(autoplot(JN(jn_test_draws(three_way = TRUE), theta_1 = "x",
                     theta_2 = "z", theta_3 = "w",
                     theta_1_vals = c(-1, 1), theta_2_vals = c(-1, 1),
                     theta_3_vals = c(-1, 1), range_size = 8),
                  "x", type = "slopes"))
})


test_that("a binary moderator in a three-way analysis becomes the facet", {
  d   <- slopes_data()
  jn  <- JN(lm(y_3b ~ x * z * b, data = d), theta_1 = "x", theta_2 = "z",
            theta_3 = "b", range_size = 15)
  lay <- function(p) ggplot2::ggplot_build(p)$layout$layout
  # the regular figure is the two-way one per level, not a heatmap
  expect_equal(nrow(lay(autoplot(jn, "x"))), 2L)
  expect_equal(nrow(lay(autoplot(jn, "x", type = "slopes"))), 2L)
  # with the binary variable focal, the heatmap stays
  expect_s3_class(autoplot(jn, "b"), "patchwork")
})


test_that("jn_plots() returns the slopes figures on request", {
  d  <- slopes_data()
  jn <- JN(lm(y_cross ~ x * z, data = d), theta_1 = "x", theta_2 = "z")
  figs <- jn_plots(jn, type = "slopes")
  expect_named(figs, c("x", "z"))

  jb <- JN(jn_test_draws(three_way = TRUE), theta_1 = "x", theta_2 = "z",
           theta_3 = "w", theta_1_vals = c(-1, 1), theta_2_vals = c(-1, 1),
           theta_3_vals = c(-1, 1), range_size = 6)
  expect_named(jn_plots(jb, type = "slopes"), c("x", "z", "w"))
  expect_named(jn_plots(jb), c("post_mean", "bayes_p"))
})
