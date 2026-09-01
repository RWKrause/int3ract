test_that("density panels are attached when the observed values are known", {
  d  <- jn_test_data()
  jn <- JN(lm(y ~ x * z, data = d), theta_1 = "x", theta_2 = "z",
           range_size = 50)

  expect_false(is.null(jn$support))
  expect_s3_class(plot(jn, which = "x"), "patchwork")
  # and are dropped on request, leaving a plain ggplot
  bare <- plot(jn, which = "x", style = jn_style(show_density = FALSE))
  expect_s3_class(bare, "ggplot")
  expect_false(inherits(bare, "patchwork"))
})


test_that("without observed values the figure is a plain ggplot", {
  jn <- JN(jn_test_draws(), theta_1 = "x", theta_2 = "z",
           theta_1_vals = seq(-2, 2, 1), theta_2_vals = seq(-2, 2, 1))
  expect_null(jn$support)
  expect_s3_class(plot(jn, which = 1), "ggplot")
  expect_false(inherits(plot(jn, which = 1), "patchwork"))
})


test_that("support supplied by hand drives the density panels", {
  draws <- jn_test_draws()
  jn <- JN(draws, theta_1 = "x", theta_2 = "z",
           theta_1_vals = seq(-2, 2, 1), theta_2_vals = seq(-2, 2, 1),
           support = list(rnorm(200), rnorm(200)))
  expect_false(is.null(jn$support))
  expect_s3_class(plot(jn, which = 1, type = "band"), "patchwork")
})


test_that("three-way heatmaps carry marginal histograms on both axes", {
  d  <- jn_test_data()
  jn <- JN(lm(y ~ x * z * w, data = d), theta_1 = "x", theta_2 = "z",
           theta_3 = "w", range_size = 12)
  p  <- plot(jn, which = "x")
  expect_s3_class(p, "patchwork")
  # main panel plus two histograms plus a spacer
  expect_equal(length(p$patches$plots) + 1L, 4L)
})


test_that("which= accepts a name or a position and refuses anything else", {
  d  <- jn_test_data()
  jn <- JN(lm(y ~ x * z, data = d), theta_1 = "x", theta_2 = "z",
           range_size = 20)
  expect_s3_class(plot(jn, which = 2), "ggplot")
  expect_s3_class(plot(jn, which = "z"), "ggplot")
  expect_error(plot(jn, which = "nope"), "must be one of")
  expect_error(plot(jn, which = 5), "between 1 and 2")
})


test_that("two-way posterior analyses can be drawn as densities or as a band", {
  jn <- JN(jn_test_draws(), theta_1 = "x", theta_2 = "z",
           theta_1_vals = seq(-2, 2, 1), theta_2_vals = seq(-2, 2, 1))
  expect_s3_class(plot(jn, which = 1, type = "density"), "ggplot")
  expect_s3_class(plot(jn, which = 1, type = "band"), "ggplot")
  expect_error(plot(jn, which = 1, type = "nonsense"))
})


test_that("plot() draws every focal variable by default", {
  d  <- jn_test_data()
  jn <- JN(lm(y ~ x * z, data = d), theta_1 = "x", theta_2 = "z",
           range_size = 20)

  # both variables take their turn as the focal one
  all_figs <- plot(jn)
  expect_type(all_figs, "list")
  expect_named(all_figs, c("x", "z"))

  # a single `which` still gives back the figure itself, not a list of one
  one <- plot(jn, which = "z")
  expect_s3_class(one, "ggplot")

  # and an explicit subset is honoured
  expect_named(plot(jn, which = c("z", "x")), c("z", "x"))
})


test_that("three-way analyses draw all three figures by default", {
  d  <- jn_test_data()
  jn <- JN(lm(y ~ x * z * w, data = d), theta_1 = "x", theta_2 = "z",
           theta_3 = "w", range_size = 8)
  expect_named(plot(jn), c("x", "z", "w"))
})


test_that("plot() returns invisibly so nothing is drawn twice", {
  d  <- jn_test_data()
  jn <- JN(lm(y ~ x * z, data = d), theta_1 = "x", theta_2 = "z",
           range_size = 20)
  expect_invisible(plot(jn, which = 1))
  expect_invisible(plot(jn))
})


test_that("autoplot() returns a figure without drawing it", {
  d  <- jn_test_data()
  jn <- JN(lm(y ~ x * z, data = d), theta_1 = "x", theta_2 = "z",
           range_size = 20)
  # bare autoplot(), not ggplot2::autoplot(): the generic must be reachable
  # from int3ract alone, without the user attaching ggplot2
  expect_true("autoplot" %in% getNamespaceExports("int3ract"))
  p <- autoplot(jn, which = "z")
  expect_s3_class(p, "ggplot")
  expect_false(is.list(p) && !inherits(p, "ggplot"))
})


test_that("plotting a grouped analysis draws the population-level one", {
  skip_if_not_installed("lme4")
  set.seed(11)
  n  <- 300
  gp <- factor(rep(seq_len(15), each = n / 15))
  d  <- data.frame(g = gp, x = rnorm(n), z = rnorm(n))
  d$y <- d$x + 0.4 * d$x * d$z + rnorm(15)[as.integer(gp)] + rnorm(n)
  res <- suppressWarnings(
    JN(lme4::lmer(y ~ x * z + (1 | g), data = d), theta_1 = "x",
       theta_2 = "z", range_size = 12, fixed_only = FALSE, group_var = "g"))

  # 15 groups times two focal variables is not what anyone wants at the console
  expect_message(plot(res), "jn_save")
})


test_that("jn_plots() returns one figure per focal variable", {
  d  <- jn_test_data()
  f2 <- jn_plots(JN(lm(y ~ x * z, data = d), theta_1 = "x", theta_2 = "z",
                    range_size = 20))
  expect_named(f2, c("x", "z"))

  jn3 <- JN(jn_test_draws(three_way = TRUE),
            theta_1 = "x", theta_2 = "z", theta_3 = "w",
            theta_1_vals = seq(-1, 1, 1), theta_2_vals = seq(-1, 1, 1),
            theta_3_vals = seq(-1, 1, 1))
  f3 <- jn_plots(jn3)
  expect_named(f3, c("post_mean", "bayes_p"))
  expect_named(f3$post_mean, c("x", "z", "w"))
})


test_that("jn_save() writes one file per figure", {
  d   <- jn_test_data()
  jn  <- JN(lm(y ~ x * z, data = d), theta_1 = "x", theta_2 = "z",
            range_size = 20)
  dir <- tempfile("jn-save")
  paths <- jn_save(jn, folder = dir, width = 4, height = 3, dpi = 72)
  expect_length(paths, 2)
  expect_true(all(file.exists(paths)))
})


test_that("figures stay extendable with the usual ggplot2 layer syntax", {
  d  <- jn_test_data()
  jn <- JN(lm(y ~ x * z, data = d), theta_1 = "x", theta_2 = "z",
           range_size = 20)
  p  <- plot(jn, which = 1, style = jn_style(show_density = FALSE)) +
    ggplot2::labs(title = "custom")
  expect_s3_class(p, "ggplot")
})
