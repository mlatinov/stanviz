# =============================================================================
# combined-predictive-check.R
# -----------------------------------------------------------------------------
# Runs the prior predictive check and the posterior predictive check against
# the same observed sample in one call, and stacks them into a single
# patchwork figure so you can see, at a glance, how much the data updated
# the model relative to what the priors alone implied.
#
# Internals reused from the existing library: `draws_matrix_of()` to pull a
# generated-quantities matrix off either fit, `theme_bayes()`/`bayes_palette()`
# for consistent styling, and the existing `plot_prior_predictive_dens()` /
# `plot_ppc_dens()` single-panel functions - this is a thin composition layer
# on top of those, not a reimplementation.
# =============================================================================


#' Prior + posterior predictive check, stacked
#'
#' Convenience wrapper that produces the prior predictive density overlay
#' (prior-only fit vs. observed sample) and the posterior predictive density
#' overlay (posterior fit vs. observed sample), and stacks them into one
#' patchwork figure. This is the "run everything at once" companion to
#' `plot_prior_predictive_dens()` and `plot_ppc_dens()` - use those directly
#' if you only want one panel or want to customize each independently.
#'
#' @param prior_model A prior-only Stan fit (likelihood switched off), i.e.
#'   the same object you'd pass to `plot_prior_predictive_dens()`.
#' @param model A fitted (posterior) Stan model, i.e. the same object you'd
#'   pass to `plot_ppc_dens()`.
#' @param y Numeric vector of observed outcomes of length `N`. Used as the
#'   comparison sample for both panels.
#' @param prior_yrep_var Name of the simulated-outcome vector family in the
#'   prior fit's `generated quantities`. Default `"y_rep"`.
#' @param post_yrep_var Name of the replicated-outcome vector family in the
#'   posterior fit's `generated quantities`. Default `"y_rep"`. Kept as a
#'   separate argument from `prior_yrep_var` in case the two Stan programs
#'   name it differently.
#' @param n_draws Number of draws to overlay per panel (passed through to
#'   both underlying single-panel functions). Default 100.
#' @param outcome_label X-axis label, shared by both panels so they line up.
#' @param ncol Number of columns in the stacked layout. Default 1 (prior on
#'   top, posterior below) - the natural reading order for "before/after
#'   seeing the data." Pass 2 to place them side by side instead.
#'
#' @return A `patchwork` composite `ggplot` object with two panels: prior
#'   predictive density overlay (top/left) and posterior predictive density
#'   overlay (bottom/right), sharing an overall title.
#' @export
#'
#' @examples
#' \dontrun{
#' plot_prior_posterior_check(prior_fit, fit, y = df$y)
#'
#' # different generated-quantity names between the two Stan programs
#' plot_prior_posterior_check(
#'   prior_fit, fit, y = df$y,
#'   prior_yrep_var = "y_sim", post_yrep_var = "y_rep"
#' )
#' }
plot_prior_posterior_check <- function(
    prior_model,
    model,
    y,
    prior_yrep_var = "y_rep",
    post_yrep_var  = "y_rep",
    n_draws        = 100,
    outcome_label  = "outcome",
    ncol           = 1
) {
  p_prior <- plot_prior_predictive_dens(
    model         = prior_model,
    yrep_var      = prior_yrep_var,
    y             = y,
    n_draws       = n_draws,
    outcome_label = outcome_label
  ) +
    ggplot2::labs(subtitle = "Before seeing the data: observed (dark) vs. prior-implied datasets (light)")

  p_post <- plot_ppc_dens(
    model         = model,
    y             = y,
    yrep_var      = post_yrep_var,
    n_draws       = n_draws,
    outcome_label = outcome_label
  ) +
    ggplot2::labs(subtitle = "After seeing the data: observed (dark) vs. replicated datasets (light)")

  patchwork::wrap_plots(list(p_prior, p_post), ncol = ncol) +
    patchwork::plot_annotation(
      title    = "Prior vs. posterior predictive check",
      subtitle = "Same observed sample against what the model generates before and after fitting",
      caption  = "A big narrowing/shift from top to bottom means the data did real work; little change may mean the priors already dominated or the likelihood is weakly identified.",
      theme    = theme_bayes()
    )
}