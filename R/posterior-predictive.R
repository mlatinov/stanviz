# =============================================================================
# posterior-predictive.R
# -----------------------------------------------------------------------------
# Posterior predictive checks (PPCs). The question: "After fitting, does the
# model regenerate data that looks like what we actually observed?"
# =============================================================================


#' Posterior predictive density overlay
#'
#' Observed density vs. many posterior-predictive densities. The workhorse PPC.
#'
#' @param model A fitted Stan model (anything `as_draws_safe()` accepts).
#' @param y Numeric vector of observed outcomes of length `N`.
#' @param yrep_var Name of the replicated-outcome vector family in
#'   `generated quantities`.
#' @param n_draws Number of posterior draws to overlay.
#' @param outcome_label X-axis label.
#'
#' @return A `ggplot` object.
#' @export
#'
#' @examples
#' \dontrun{
#' plot_ppc_dens(fit, y = df$y, yrep_var = "y_rep")
#' }
plot_ppc_dens <- function(
    model,
    y,
    yrep_var      = "y_rep",
    n_draws       = 100,
    outcome_label = "outcome"
) {
  yrep <- draws_matrix_of(model, yrep_var)
  if (nrow(yrep) > n_draws) {
    yrep <- yrep[sample.int(nrow(yrep), n_draws), , drop = FALSE]
  }

  bayesplot::ppc_dens_overlay(y = y, yrep = yrep) +
    ggplot2::labs(
      x        = outcome_label,
      title    = "Posterior predictive check",
      subtitle = "Observed (dark) vs. replicated datasets (light)",
      caption  = "Systematic gaps indicate the likelihood is missing structure."
    ) +
    theme_bayes()
}


#' Posterior predictive test statistic
#'
#' Compare a single test statistic between observed and replicated data.
#'
#' @param model A fitted Stan model.
#' @param y Numeric vector of observed outcomes.
#' @param yrep_var Name of the replicated-outcome vector family.
#' @param stat A function name string (`"mean"`, `"sd"`, `"max"`, `"min"`)
#'   or a function.
#'
#' @return A `ggplot` object.
#' @export
#'
#' @examples
#' \dontrun{
#' plot_ppc_stat(fit, y = df$y, stat = "sd")
#' }
plot_ppc_stat <- function(
    model,
    y,
    yrep_var = "y_rep",
    stat     = "mean"
) {
  yrep <- draws_matrix_of(model, yrep_var)
  bayesplot::ppc_stat(y = y, yrep = yrep, stat = stat) +
    ggplot2::labs(
      title    = paste0("PPC test statistic: ",
                        if (is.character(stat)) stat else "custom"),
      subtitle = "Observed statistic (dark) within the replicated distribution"
    ) +
    theme_bayes()
}


#' Posterior predictive test-statistic grid
#'
#' Several test statistics stitched together with `patchwork`. Covers centre,
#' spread, and tails at a glance.
#'
#' @param model A fitted Stan model.
#' @param y Numeric vector of observed outcomes.
#' @param yrep_var Name of the replicated-outcome vector family.
#' @param stats Character vector of statistic names.
#'
#' @return A `patchwork` composite ggplot.
#' @export
#'
#' @examples
#' \dontrun{
#' plot_ppc_stat_grid(fit, y = df$y, stats = c("mean", "sd", "max"))
#' }
plot_ppc_stat_grid <- function(
    model,
    y,
    yrep_var = "y_rep",
    stats    = c("mean", "sd", "max", "min")
) {
  yrep <- draws_matrix_of(model, yrep_var)

  plots <- lapply(stats, function(s) {
    bayesplot::ppc_stat(y = y, yrep = yrep, stat = s) +
      ggplot2::labs(title = s, x = NULL) +
      theme_bayes(base_size = 10) +
      ggplot2::theme(legend.position = "none")
  })

  patchwork::wrap_plots(plots) +
    patchwork::plot_annotation(
      title    = "Posterior predictive test statistics",
      subtitle = "Each panel: observed statistic vs. replicated distribution",
      theme    = theme_bayes()
    )
}


#' Per-observation predictive intervals
#'
#' Per-observation predictive intervals against observed points. Helpful for
#' identifying which observations the model fits poorly. Optionally orders
#' the x-axis by the observed value to reveal trends.
#'
#' @param model A fitted Stan model.
#' @param y Numeric vector of observed outcomes.
#' @param yrep_var Name of the replicated-outcome vector family.
#' @param x Optional numeric vector for the x-axis. Defaults to observation
#'   index.
#' @param order_by_y Logical; sort observations by `y` for trend inspection.
#' @param outcome_label Y-axis label.
#'
#' @return A `ggplot` object.
#' @export
#'
#' @examples
#' \dontrun{
#' plot_ppc_intervals(fit, y = df$y, yrep_var = "y_rep")
#' }
plot_ppc_intervals <- function(
    model,
    y,
    yrep_var      = "y_rep",
    x             = NULL,
    order_by_y    = TRUE,
    outcome_label = "outcome"
) {
  yrep <- draws_matrix_of(model, yrep_var)

  if (is.null(x)) x <- seq_along(y)
  if (order_by_y) {
    ord  <- order(y)
    y    <- y[ord]
    yrep <- yrep[, ord, drop = FALSE]
    x    <- seq_along(y)
  }

  bayesplot::ppc_intervals(y = y, yrep = yrep, x = x) +
    ggplot2::labs(
      x = if (order_by_y) "observation (ordered by observed value)" else "observation",
      y = outcome_label,
      title    = "Per-observation predictive intervals",
      subtitle = "Points = observed; intervals = posterior predictive (50% / 90%)"
    ) +
    theme_bayes()
}


#' Predictive error scatter
#'
#' Average predictive error (`y - mean(yrep)`) against the predictor or
#' fitted value. Surfaces heteroskedasticity and nonlinearity that the
#' density plot hides.
#'
#' @param model A fitted Stan model.
#' @param y Numeric vector of observed outcomes.
#' @param yrep_var Name of the replicated-outcome vector family.
#' @param x Optional predictor or fitted vector; defaults to `colMeans(yrep)`.
#' @param x_label X-axis label.
#'
#' @return A `ggplot` object.
#' @export
#'
#' @examples
#' \dontrun{
#' plot_ppc_error_scatter(fit, y = df$y, x = df$predictor)
#' }
plot_ppc_error_scatter <- function(
    model,
    y,
    yrep_var = "y_rep",
    x        = NULL,
    x_label  = "fitted / predictor"
) {
  yrep <- draws_matrix_of(model, yrep_var)
  yhat <- colMeans(yrep)
  err  <- y - yhat
  if (is.null(x)) x <- yhat

  df  <- tibble::tibble(x = x, err = err)
  pal <- bayes_palette()

  ggplot2::ggplot(df, ggplot2::aes(x = .data$x, y = .data$err)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey30") +
    ggplot2::geom_point(alpha = 0.5, colour = pal[["effect"]]) +
    ggplot2::geom_smooth(method = "loess", formula = y ~ x,
                         colour = pal[["control"]], se = TRUE) +
    ggplot2::labs(
      x = x_label, y = "predictive error (y - E[y_rep])",
      title    = "Predictive error structure",
      subtitle = "Flat band around zero = well calibrated; trends = missing structure"
    ) +
    theme_bayes()
}


#' Probability-integral-transform (PIT) histogram
#'
#' Calibration check. If the model is calibrated, PIT values are
#' `Uniform(0,1)` and the histogram is flat. U-shapes mean the predictive is
#' too narrow; dome shapes mean too wide. Pass the `pit` vector that the Stan
#' template computes in `generated quantities`.
#'
#' @param model A fitted Stan model.
#' @param pit_var Name of the PIT vector family in `generated quantities`.
#' @param bins Number of histogram bins.
#'
#' @return A `ggplot` object.
#' @export
#'
#' @examples
#' \dontrun{
#' plot_ppc_pit(fit, pit_var = "pit")
#' }
plot_ppc_pit <- function(
    model,
    pit_var = "pit",
    bins    = 20
) {
  pit_mat <- draws_matrix_of(model, pit_var)
  pit     <- colMeans(pit_mat)
  df      <- tibble::tibble(pit = pit)
  pal     <- bayes_palette()

  n_obs    <- length(pit)
  expected <- n_obs / bins

  ggplot2::ggplot(df, ggplot2::aes(x = .data$pit)) +
    ggplot2::geom_histogram(bins = bins, fill = pal[["effect"]],
                            alpha = 0.7, colour = "white") +
    ggplot2::geom_hline(yintercept = expected, linetype = "dashed",
                        colour = "grey30") +
    ggplot2::labs(
      x = "PIT value", y = "count",
      title    = "Calibration: PIT histogram",
      subtitle = "Flat = calibrated; U-shape = over-confident; dome = under-confident"
    ) +
    theme_bayes()
}