# Synthetic posterior draws exercising every shape the plotting helpers need
# to cope with: a plain vector, a genuine scalar, an `array[J] vector[N]`
# (2 dimensions), an `array[K] matrix[R, C]` (3 dimensions), and a small
# `array[J] vector[N]` sized to double as a hierarchical "group_var".
make_fake_draws <- function(n_draws = 20, n_chains = 2, seed = 42) {
  set.seed(seed)

  mu_names <- paste0("mu[", 1:4, "]")                       # vector[4]

  sc_idx    <- expand.grid(n = 1:5, a = 1:3)                # array[3] vector[5]
  sc_names  <- sprintf("samples_combined[%d,%d]", sc_idx$a, sc_idx$n)

  blk_idx   <- expand.grid(c = 1:3, r = 1:3, a = 1:2)       # array[2] matrix[3,3]
  blk_names <- sprintf("blk[%d,%d,%d]", blk_idx$a, blk_idx$r, blk_idx$c)

  z_idx     <- expand.grid(n = 1:2, a = 1:3)                # array[3] vector[2]
  z_names   <- sprintf("z[%d,%d]", z_idx$a, z_idx$n)

  y_rep_names <- paste0("y_rep[", 1:6, "]")                 # vector[6]

  scalar_names <- "sigma"

  all_names <- c(mu_names, sc_names, blk_names, z_names, y_rep_names, scalar_names)
  n_total   <- n_draws * n_chains

  arr <- array(
    stats::rnorm(n_total * length(all_names)),
    dim      = c(n_draws, n_chains, length(all_names)),
    dimnames = list(NULL, NULL, all_names)
  )
  posterior::as_draws_array(arr)
}

# Synthetic draws for a nested predictor-grid quantity, e.g.
# `array[K] vector[G] mu_grid` in `generated quantities` - the shape
# plot_effect_curve()/plot_effect_curve_grouped() are built to consume.
# Group `a=1` rises steeply with the grid index, `a=2` rises gently, so the
# two groups' posterior curves are visibly distinguishable in tests.
make_fake_grid_draws <- function(n_draws = 30, n_chains = 1, n_group = 2, n_grid = 6, seed = 7) {
  set.seed(seed)

  idx  <- expand.grid(g = seq_len(n_grid), a = seq_len(n_group))
  nms  <- sprintf("mu_grid[%d,%d]", idx$a, idx$g)
  n_total <- n_draws * n_chains

  slope    <- ifelse(idx$a == 1, 2, 0.3)
  base_val <- slope * idx$g

  arr <- array(
    rep(base_val, each = n_total) + stats::rnorm(n_total * length(nms), sd = 0.2),
    dim      = c(n_draws, n_chains, length(nms)),
    dimnames = list(NULL, NULL, nms)
  )
  posterior::as_draws_array(arr)
}
