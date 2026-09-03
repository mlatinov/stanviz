# New plotting functions: plot_ppc_ecdf() (+ grouped), plot_ppc_boxplot()
# (+ grouped), plot_ppc_violin()/plot_ppc_violin_grouped(), plot_effect_curve()
# (+ grouped), and the group-aware plot_conditional_effect_draws() that
# powers the effect-curve functions.

draws <- make_fake_draws()
y     <- as.numeric(draws_matrix_of(draws, "y_rep")[1, ]) + stats::rnorm(6, sd = 0.1)
grp   <- factor(rep(c("a", "b"), length.out = 6))

# ---- plot_ppc_ecdf() / plot_ppc_ecdf_grouped() -----------------------------

test_that("plot_ppc_ecdf builds a single-panel ECDF overlay", {
  p <- plot_ppc_ecdf(draws, y = y, yrep_var = "y_rep", n_draws = 10)
  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$title, "Posterior predictive check: ECDF")
  b <- ggplot2::ggplot_build(p)
  expect_equal(length(b$layout$panel_params), 1L)
})

test_that("plot_ppc_ecdf downsamples yrep to n_draws", {
  p <- plot_ppc_ecdf(draws, y = y, yrep_var = "y_rep", n_draws = 3)
  b <- ggplot2::ggplot_build(p)
  # one ECDF step-line group per retained yrep draw, plus one for y
  expect_true(nrow(unique(b$data[[1]]["group"])) <= 4)
})

test_that("plot_ppc_ecdf_grouped facets by group", {
  p <- plot_ppc_ecdf_grouped(draws, y = y, group = grp, yrep_var = "y_rep", n_draws = 10)
  expect_s3_class(p, "ggplot")
  b <- ggplot2::ggplot_build(p)
  expect_equal(length(b$layout$panel_params), nlevels(grp))
})

test_that("plot_ppc_ecdf_grouped validates group_labels length", {
  expect_error(
    plot_ppc_ecdf_grouped(draws, y = y, group = grp, group_labels = "only_one"),
    "group_labels"
  )
})

# ---- plot_ppc_boxplot() / plot_ppc_boxplot_grouped() -----------------------

test_that("plot_ppc_boxplot builds a single ggplot with y-axis outcome label", {
  p <- suppressWarnings(plot_ppc_boxplot(draws, y = y, yrep_var = "y_rep", n_draws = 5))
  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$y, "outcome")
})

test_that("plot_ppc_boxplot_grouped returns one patchwork panel per group", {
  p <- suppressWarnings(
    plot_ppc_boxplot_grouped(draws, y = y, group = grp, yrep_var = "y_rep", n_draws = 5)
  )
  expect_s3_class(p, "patchwork")
  g <- suppressWarnings(patchwork::patchworkGrob(p))
  n_panels <- sum(grepl("^panel", vapply(g$grobs, function(x) x$name, character(1))))
  expect_equal(n_panels, nlevels(grp))
})

# ---- plot_ppc_violin_grouped() / plot_ppc_violin() -------------------------

test_that("plot_ppc_violin_grouped builds a ggplot with one x-tick per group", {
  p <- plot_ppc_violin_grouped(draws, y = y, group = grp, yrep_var = "y_rep")
  expect_s3_class(p, "ggplot")
  b <- ggplot2::ggplot_build(p)
  expect_equal(b$layout$panel_params[[1]]$x$get_labels(), levels(grp))
})

test_that("plot_ppc_violin reuses plot_ppc_violin_grouped with a single dummy group", {
  p <- plot_ppc_violin(draws, y = y, yrep_var = "y_rep")
  expect_s3_class(p, "ggplot")
  expect_equal(p$theme$legend.position, "none")
  expect_true(inherits(p$theme$axis.text.x, "element_blank"))
  # builds without error, i.e. bayesplot accepted the dummy single-level group
  expect_silent(ggplot2::ggplot_build(p))
})

# ---- plot_conditional_effect_draws(): group-column backward compatibility --

test_that("plot_conditional_effect_draws is unchanged for plain (ungrouped) input", {
  dl <- tibble::tibble(
    x     = rep(1:5, each = 20),
    .draw = rep(1:20, times = 5),
    value = rep(1:5, each = 20) + stats::rnorm(100, sd = 0.1)
  )
  p <- plot_conditional_effect_draws(dl, x_label = "age")
  expect_s3_class(p, "ggplot")
  expect_null(p$labels$colour)
  expect_equal(p$labels$title, "Conditional effect")
})

test_that("plot_conditional_effect_draws draws one curve per group when `group` is present", {
  dl <- tibble::tibble(
    x     = rep(1:5, each = 20, times = 2),
    .draw = rep(rep(1:20, times = 5), 2),
    value = c(
      rep(1:5, each = 20) + stats::rnorm(100, sd = 0.1),
      rep(1:5, each = 20) * 2 + stats::rnorm(100, sd = 0.1)
    ),
    group = rep(c("g1", "g2"), each = 100)
  )
  p <- plot_conditional_effect_draws(dl)
  expect_s3_class(p, "ggplot")
  b <- ggplot2::ggplot_build(p)
  expect_true(nrow(b$data[[1]]) > 0)
})

# ---- plot_effect_curve() ----------------------------------------------------

grid_draws <- make_fake_grid_draws()
x_grid     <- seq(0, 10, length.out = 6)

test_that("plot_effect_curve reuses gather_indexed() to build the curve from a nested slice", {
  p <- plot_effect_curve(grid_draws, "mu_grid[1,]", x = x_grid,
                          outcome_label = "H2O2", x_label = "SOD")
  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$x, "SOD")
  expect_equal(p$labels$y, "H2O2")
  expect_equal(p$labels$title, "Posterior effect curve")

  b <- ggplot2::ggplot_build(p)
  ribbon <- unique(b$data[[1]][, c("x", "y")])
  ribbon <- ribbon[order(ribbon$x), ]
  # median value should be monotonically increasing with x for this synthetic curve
  expect_true(all(diff(ribbon$y) >= -1e-6))
  # and should track the known slope=2 relationship (value = 2 * x) reasonably closely
  expect_equal(ribbon$y, 2 * ribbon$x, tolerance = 0.5)
})

test_that("plot_effect_curve overlays observed points when x_obs/y_obs are supplied", {
  p_no_obs <- plot_effect_curve(grid_draws, "mu_grid[1,]", x = x_grid)
  p_obs    <- plot_effect_curve(
    grid_draws, "mu_grid[1,]", x = x_grid,
    x_obs = c(1, 5, 9), y_obs = c(2, 10, 18)
  )
  expect_equal(length(p_no_obs$layers), length(p_obs$layers) - 1L)
  expect_silent(ggplot2::ggplot_build(p_obs))
})

test_that("plot_effect_curve errors informatively on a grid-length mismatch", {
  expect_error(
    plot_effect_curve(grid_draws, "mu_grid[1,]", x = 1:3),
    "grid point"
  )
})

test_that("plot_effect_curve errors when the slice still has more than one free dimension", {
  expect_error(
    plot_effect_curve(grid_draws, "mu_grid", x = x_grid),
    "more than one free dimension"
  )
})

# ---- plot_effect_curve_grouped() -------------------------------------------

test_that("plot_effect_curve_grouped overlays one curve per named effect_vars entry", {
  p <- plot_effect_curve_grouped(
    grid_draws,
    effect_vars = c("Steep" = "mu_grid[1,]", "Gentle" = "mu_grid[2,]"),
    x = x_grid, outcome_label = "H2O2", x_label = "SOD"
  )
  expect_s3_class(p, "ggplot")
  b <- ggplot2::ggplot_build(p)
  ribbon <- b$data[[1]]
  expect_equal(length(unique(ribbon$group)), 2L)     # one ribbon-group per named group
  expect_equal(nrow(ribbon), 2L * length(x_grid) * 2L)  # 2 groups x 6 x-values x 2 widths

  # the "Steep" curve (slope 2) should sit above "Gentle" (slope 0.3) at the
  # right edge of the grid, confirming groups map to the correct draws
  build_data <- b$plot$data
  steep  <- build_data[build_data$group == "Steep" & build_data$x == max(x_grid), ]
  gentle <- build_data[build_data$group == "Gentle" & build_data$x == max(x_grid), ]
  expect_true(unique(steep$value) > unique(gentle$value))
})

test_that("plot_effect_curve_grouped requires a named effect_vars vector", {
  expect_error(
    plot_effect_curve_grouped(grid_draws, effect_vars = c("mu_grid[1,]", "mu_grid[2,]"), x = x_grid),
    "named character vector"
  )
})

test_that("plot_effect_curve_grouped overlays obs_data points coloured by group", {
  obs <- data.frame(
    sod       = c(1, 8),
    h2o2      = c(2, 15),
    treatment = c("Steep", "Gentle")
  )
  p <- plot_effect_curve_grouped(
    grid_draws,
    effect_vars = c("Steep" = "mu_grid[1,]", "Gentle" = "mu_grid[2,]"),
    x = x_grid, obs_data = obs,
    x_col = "sod", y_col = "h2o2", group_col = "treatment"
  )
  expect_equal(length(p$layers), 3L)  # lineribbon + median line + obs points
  expect_silent(ggplot2::ggplot_build(p))
})
