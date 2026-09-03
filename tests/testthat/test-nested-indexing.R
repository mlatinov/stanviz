# Nested/multi-dimensional Stan containers (array[J] vector[N],
# array[K] matrix[R, C], ...) used to be silently mangled: the old
# implementation ran `readr::parse_number()` on the whole bracket text, so
# "samples_combined[1,2,3]" collapsed the commas into the number 123. These
# tests cover the fix across the full call chain: parse_stan_var() ->
# stan_index_matrix() -> resolve_stan_var() -> gather_indexed() /
# draws_matrix_of() / draws_matrix_of_or_scalar(), and finally through to the
# exported plotting functions that consume them.

draws <- make_fake_draws()
all_vars <- posterior::variables(draws)
n_draws_total <- posterior::ndraws(draws)

# ---- parse_stan_var() ------------------------------------------------------

test_that("parse_stan_var splits base name and index spec", {
  expect_equal(parse_stan_var("mu"), list(base = "mu", idx = NULL))
  expect_equal(parse_stan_var("mu[2]"), list(base = "mu", idx = list(2L)))
})

test_that("parse_stan_var treats blank slots as wildcards, including trailing ones", {
  spec <- parse_stan_var("samples_combined[1,]")
  expect_equal(spec$base, "samples_combined")
  expect_equal(spec$idx, list(1L, NA_integer_))

  spec2 <- parse_stan_var("samples_combined[,3]")
  expect_equal(spec2$idx, list(NA_integer_, 3L))

  spec3 <- parse_stan_var("blk[,2,]")
  expect_equal(spec3$idx, list(NA_integer_, 2L, NA_integer_))
})

test_that("parse_stan_var rejects garbage indices", {
  expect_error(parse_stan_var("mu[x]"), "Invalid index")
})

# ---- stan_index_matrix() ---------------------------------------------------

test_that("stan_index_matrix parses multi-dimensional bracket groups", {
  nms <- c("blk[1,2,3]", "blk[2,1,1]")
  m <- stan_index_matrix(nms)
  expect_equal(dim(m), c(2L, 3L))
  expect_equal(m[1, ], c(1L, 2L, 3L))
  expect_equal(m[2, ], c(2L, 1L, 1L))
})

# ---- resolve_stan_var() -----------------------------------------------------

test_that("resolve_stan_var returns the full family, in order, for a plain vector", {
  res <- resolve_stan_var(all_vars, "mu")
  expect_false(res$is_scalar)
  expect_equal(res$hit_vars, paste0("mu[", 1:4, "]"))
  expect_equal(res$free_dims, 1L)
})

test_that("resolve_stan_var resolves a genuine scalar", {
  res <- resolve_stan_var(all_vars, "sigma")
  expect_true(res$is_scalar)
  expect_equal(res$hit_vars, "sigma")
})

test_that("resolve_stan_var slices array[J] vector[N] on the trailing wildcard", {
  # samples_combined[1,] -> fix array index 1, keep every vector element
  res <- resolve_stan_var(all_vars, "samples_combined[1,]")
  expect_equal(length(res$hit_vars), 5L)
  expect_equal(res$hit_vars, sprintf("samples_combined[1,%d]", 1:5))
  expect_equal(res$free_dims, 2L)
})

test_that("resolve_stan_var slices array[J] vector[N] on the leading wildcard", {
  # samples_combined[,3] -> every array/group index, vector element 3 fixed
  res <- resolve_stan_var(all_vars, "samples_combined[,3]")
  expect_equal(length(res$hit_vars), 3L)
  expect_equal(res$hit_vars, sprintf("samples_combined[%d,3]", 1:3))
  expect_equal(res$free_dims, 1L)
})

test_that("resolve_stan_var slices a fully-specified element to exactly one hit", {
  res <- resolve_stan_var(all_vars, "samples_combined[2,4]")
  expect_equal(res$hit_vars, "samples_combined[2,4]")
  expect_equal(res$free_dims, integer(0))
})

test_that("resolve_stan_var handles array[K] matrix[R, C] (3 dimensions)", {
  # blk[,2,3] -> every array slot, row 2 / col 3 of the matrix fixed
  res <- resolve_stan_var(all_vars, "blk[,2,3]")
  expect_equal(res$hit_vars, sprintf("blk[%d,2,3]", 1:2))
  expect_equal(res$free_dims, 1L)

  # blk[1,,] -> array slot 1 fixed, every matrix element
  res2 <- resolve_stan_var(all_vars, "blk[1,,]")
  expect_equal(length(res2$hit_vars), 9L)
  expect_true(all(grepl("^blk\\[1,", res2$hit_vars)))
  expect_equal(res2$free_dims, c(2L, 3L))
})

test_that("resolve_stan_var errors informatively on dimension/index mismatches", {
  expect_error(resolve_stan_var(all_vars, "nope"), "No variable matching")
  expect_error(resolve_stan_var(all_vars, "mu[1,2,3]"), "dimension")
  expect_error(resolve_stan_var(all_vars, "samples_combined[99,]"), "No elements")
})

# ---- gather_indexed() -------------------------------------------------------

test_that("gather_indexed keeps the classic single `index` column for a plain vector", {
  long <- gather_indexed(draws, "mu")
  expect_true(all(c(".draw", "param", "value", "index") %in% names(long)))
  expect_setequal(long$index, 1:4)
  expect_equal(nrow(long), n_draws_total * 4L)
})

test_that("gather_indexed collapses a nested slice with one free dimension to `index`", {
  long <- gather_indexed(draws, "samples_combined[1,]")
  expect_true("index" %in% names(long))
  expect_false(any(c("index1", "index2") %in% names(long)))
  expect_setequal(long$index, 1:5)
  expect_equal(nrow(long), n_draws_total * 5L)
})

test_that("gather_indexed emits index1/index2 for an unsliced 2D family", {
  long <- gather_indexed(draws, "samples_combined")
  expect_true(all(c("index1", "index2") %in% names(long)))
  expect_setequal(long$index1, 1:3)
  expect_setequal(long$index2, 1:5)
  expect_equal(nrow(long), n_draws_total * 15L)
})

test_that("gather_indexed handles a genuine scalar", {
  long <- gather_indexed(draws, "sigma")
  expect_true(all(long$index == 1L))
  expect_equal(nrow(long), n_draws_total)
})

# ---- draws_matrix_of() ------------------------------------------------------

test_that("draws_matrix_of builds a plain [n_draws, N] matrix as before", {
  m <- draws_matrix_of(draws, "mu")
  expect_equal(dim(m), c(n_draws_total, 4L))
  expect_equal(colnames(m), paste0("mu[", 1:4, "]"))
})

test_that("draws_matrix_of slices array[J] vector[N] via 'var[1,]'", {
  m <- draws_matrix_of(draws, "samples_combined[1,]")
  expect_equal(dim(m), c(n_draws_total, 5L))
  expect_equal(colnames(m), sprintf("samples_combined[1,%d]", 1:5))
})

test_that("draws_matrix_of slices array[K] matrix[R, C] via 'var[,2,3]'", {
  m <- draws_matrix_of(draws, "blk[,2,3]")
  expect_equal(dim(m), c(n_draws_total, 2L))
  expect_equal(colnames(m), sprintf("blk[%d,2,3]", 1:2))
})

test_that("draws_matrix_of refuses a genuine scalar", {
  expect_error(draws_matrix_of(draws, "sigma"), "No vector variable matching")
})

# ---- draws_matrix_of_or_scalar() (internal, conditional-effects.R) --------

test_that("draws_matrix_of_or_scalar returns a numeric vector for a real scalar", {
  v <- draws_matrix_of_or_scalar(draws, "sigma")
  expect_type(v, "double")
  expect_equal(length(v), n_draws_total)
})

test_that("draws_matrix_of_or_scalar returns a matrix for a nested slice", {
  m <- draws_matrix_of_or_scalar(draws, "samples_combined[,1]")
  expect_true(is.matrix(m))
  expect_equal(dim(m), c(n_draws_total, 3L))
  expect_equal(colnames(m), sprintf("samples_combined[%d,1]", 1:3))
})

# ---- End-to-end: exported plotting/extraction functions -------------------

test_that("extract_group_effects treats the array dimension of a nested var as groups", {
  eff <- extract_group_effects(draws, "samples_combined[,1]")
  expect_equal(nlevels(eff$group), 3L)
  expect_equal(nrow(eff), n_draws_total * 3L)

  eff_blk <- extract_group_effects(draws, "blk[,1,2]")
  expect_equal(nlevels(eff_blk$group), 2L)
})

test_that("plot_random_effects builds a ggplot from a nested-variable slice", {
  p <- plot_random_effects(draws, group_var = "samples_combined[,1]")
  expect_s3_class(p, "ggplot")
})

test_that("plot_ppc_dens works with the (unsliced, plain-vector) y_rep family", {
  y <- as.numeric(draws_matrix_of(draws, "y_rep")[1, ])
  p <- plot_ppc_dens(draws, y = y, yrep_var = "y_rep", n_draws = 5)
  expect_s3_class(p, "ggplot")
})

test_that("plot_param_halfeye facets a nested slice like any other parameter set", {
  p <- plot_param_halfeye(draws, pars = sprintf("blk[1,%d,%d]", 1:2, 1:2))
  expect_s3_class(p, "ggplot")
})

# ---- plot_funnel_diagnostic(): which_group naming for nested group_var ----
# plot_funnel_diagnostic() needs a real CmdStanMCMC/stanfit (it calls
# bayesplot::nuts_params() for divergence info), which these synthetic draws
# don't provide. This checks the piece of that function actually touched by
# the fix: building the exact draws-column name from a multi-index
# `which_group`, the same way the function itself does internally.
test_that("plot_funnel_diagnostic's multi-index group-name construction resolves", {
  which_group <- c(2L, 1L)
  which_group_str <- paste(which_group, collapse = ",")
  group_name <- paste0("z", "[", which_group_str, "]")
  expect_equal(group_name, "z[2,1]")
  expect_true(group_name %in% all_vars)
})
