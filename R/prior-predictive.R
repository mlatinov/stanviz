# =============================================================================
# prior-predictive.R
# -----------------------------------------------------------------------------
# Prior predictive checks. The question these answer: "Before seeing the
# data, does my model generate outcomes that are even remotely plausible?"
#
# Workflow assumption: you have run the SAME Stan program with the likelihood
# switched off (e.g. a `prior_only` flag in `data`, or a separate
# *_prior.stan file that only samples from the priors and produces a
# `*_rep` / `y_sim` generated quantity). You then pass that prior-only fit
# here exactly as you would a posterior fit.
# =============================================================================


#' Prior predictive density overlay
#'
#' Overlay many prior-predictive simulated datasets against the (optional)
#' observed outcome. If `y` is supplied it is drawn on top so you can see
#' whether the priors at least bracket reality; if not, you just see the
#' prior-implied spread of outcomes.
#'
#' @param model A prior-only Stan fit (anything `as_draws_safe()` accepts).
#' @param yrep_var Name of the simulated-outcome vector family in
#'   `generated quantities`.
#' @param y Optional numeric vector of observed outcomes of length `N`.
#' @param n_draws How many prior draws to overlay (keeps the plot readable).
#' @param outcome_label X-axis label.
#'
#' @return A `ggplot` object.
#' @export
#'
#' @examples
#' \dontrun{
#' plot_prior_predictive_dens(prior_fit, yrep_var = "y_rep", y = df$y)
#' }
plot_prior_predictive_dens <- function(
    model,
    yrep_var      = "y_rep",
    y             = NULL,
    n_draws       = 100,
    outcome_label = "outcome"
) {
  yrep <- draws_matrix_of(model, yrep_var)

  if (nrow(yrep) > n_draws) {
    idx  <- sample.int(nrow(yrep), n_draws)
    yrep <- yrep[idx, , drop = FALSE]
  }

  p <-
    if (is.null(y)) {
      bayesplot::ppc_dens_overlay(
        y    = apply(yrep, 2, stats::median),
        yrep = yrep
      )
    } else {
      bayesplot::ppc_dens_overlay(y = y, yrep = yrep)
    }

  p +
    ggplot2::labs(
      x        = outcome_label,
      y        = NULL,
      title    = "Prior predictive distribution",
      subtitle = if (is.null(y))
        "Outcomes implied by the priors alone"
      else
        "Observed outcome (dark) against prior-implied datasets (light)",
      caption  = "If the dark line sits far in the tails, the priors are too tight or mis-scaled."
    ) +
    theme_bayes()
}


#' Prior predictive summary statistic
#'
#' Distribution of a summary statistic (mean, sd, max, ...) across
#' prior-predictive datasets, with the observed statistic optionally marked.
#' Often more informative than the full density: it tells you whether the
#' priors imply, say, a mean outcome of 0 +/- 2 or 0 +/- 2000.
#'
#' @param model A prior-only Stan fit.
#' @param yrep_var Name of the simulated-outcome vector family in
#'   `generated quantities`.
#' @param y Optional numeric vector of observed outcomes; if supplied, the
#'   observed statistic is overlaid.
#' @param stat A function name as a string (`"mean"`, `"sd"`, `"max"`,
#'   `"min"`) or a function.
#' @param outcome_label Used to label the x-axis.
#'
#' @return A `ggplot` object.
#' @export
#'
#' @examples
#' \dontrun{
#' plot_prior_predictive_stat(prior_fit, stat = "sd", y = df$y)
#' }
plot_prior_predictive_stat <- function(
    model,
    yrep_var      = "y_rep",
    y             = NULL,
    stat          = "mean",
    outcome_label = "outcome"
) {
  yrep     <- draws_matrix_of(model, yrep_var)
  stat_fun <- match.fun(stat)
  stat_lab <- if (is.character(stat)) stat else "statistic"

  sim_stat <- apply(yrep, 1, stat_fun)
  df       <- tibble::tibble(value = sim_stat)
  pal      <- bayes_palette()

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$value)) +
    ggdist::stat_halfeye(
      .width         = c(0.66, 0.95),
      point_interval = "median_qi",
      fill           = pal[["effect"]],
      alpha          = 0.7
    ) +
    ggplot2::labs(
      x        = paste0(stat_lab, " of ", outcome_label),
      y        = NULL,
      title    = paste0("Prior predictive ", stat_lab),
      subtitle = "Range of this statistic implied by the priors"
    ) +
    theme_bayes()

  if (!is.null(y)) {
    obs <- stat_fun(y)
    p <- p +
      ggplot2::geom_vline(xintercept = obs, colour = pal[["null"]], linewidth = 1) +
      ggplot2::annotate(
        "text", x = obs, y = 0, label = "observed",
        hjust = -0.1, vjust = -0.5, colour = pal[["null"]], size = 3
      )
  }
  p
}


#' Realised parameter priors
#'
#' Visualise the marginal priors actually realised for a set of *parameters*
#' (not outcomes). Useful sanity check that, e.g., a `normal(0, 0.5)` prior on
#' a slope really does concentrate where you think.
#'
#' @param model A prior-only Stan fit.
#' @param pars Optional character vector of parameter names to plot. `NULL`
#'   plots all non-bookkeeping parameters in the draws.
#' @param facet_ncol Number of facet columns.
#'
#' @return A `ggplot` object.
#' @export
#'
#' @examples
#' \dontrun{
#' plot_prior_draws(prior_fit, pars = c("alpha", "beta", "sigma"))
#' }
plot_prior_draws <- function(
    model,
    pars       = NULL,
    facet_ncol = 3
) {
  draws <- as_draws_safe(model, variables = pars)
  long  <- draws |>
    posterior::as_draws_df() |>
    tibble::as_tibble() |>
    dplyr::select(-dplyr::any_of(c(".chain", ".iteration", ".draw"))) |>
    tidyr::pivot_longer(
      dplyr::everything(),
      names_to  = "param",
      values_to = "value"
    )

  pal <- bayes_palette()

  ggplot2::ggplot(long, ggplot2::aes(x = .data$value)) +
    ggdist::stat_halfeye(
      .width         = c(0.66, 0.95),
      point_interval = "median_qi",
      fill           = pal[["effect"]],
      alpha          = 0.7
    ) +
    ggplot2::facet_wrap(~ .data$param, scales = "free", ncol = facet_ncol) +
    ggplot2::labs(
      x = NULL, y = NULL,
      title    = "Realised parameter priors",
      subtitle = "What the priors imply for each parameter before seeing data"
    ) +
    theme_bayes()
}