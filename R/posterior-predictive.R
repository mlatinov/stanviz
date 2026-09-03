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


#' Posterior predictive ECDF overlay
#'
#' Empirical-CDF counterpart to [plot_ppc_dens()]. Often shows tail and skew
#' mismatches more clearly than an overlaid density, especially for
#' bounded, skewed, or heavy-tailed outcomes.
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
#' plot_ppc_ecdf(fit, y = df$y, yrep_var = "y_rep")
#' }
plot_ppc_ecdf <- function(
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

  bayesplot::ppc_ecdf_overlay(y = y, yrep = yrep) +
    ggplot2::labs(
      x        = outcome_label,
      title    = "Posterior predictive check: ECDF",
      subtitle = "Observed (dark) vs. replicated (light) empirical CDFs",
      caption  = "Gaps concentrated in the tails point to skew or shape the model isn't capturing."
    ) +
    theme_bayes()
}


#' Posterior predictive ECDF overlay, by group
#'
#' Grouped counterpart to [plot_ppc_ecdf()]: one panel per group, so you can
#' see whether the model reproduces the observed distribution's shape
#' equally well across treatment arms, sites, days, etc.
#'
#' @param model A fitted Stan model.
#' @param y Numeric vector of observed outcomes, length `N`.
#' @param group Integer or factor vector of group memberships, length `N`.
#' @param yrep_var Name of the replicated-outcome vector family.
#' @param n_draws Number of posterior draws to overlay per facet.
#' @param group_labels Optional character vector of group display labels.
#' @param outcome_label X-axis label.
#'
#' @return A `ggplot` object.
#' @export
#'
#' @examples
#' \dontrun{
#' plot_ppc_ecdf_grouped(fit, y = df$y, group = df$treatment)
#' }
plot_ppc_ecdf_grouped <- function(
    model,
    y,
    group,
    yrep_var      = "y_rep",
    n_draws       = 50,
    group_labels  = NULL,
    outcome_label = "outcome"
) {
  yrep <- draws_matrix_of(model, yrep_var)
  if (nrow(yrep) > n_draws) {
    yrep <- yrep[sample.int(nrow(yrep), n_draws), , drop = FALSE]
  }
  group <- as_group_factor(group, group_labels)

  bayesplot::ppc_ecdf_overlay_grouped(y = y, yrep = yrep, group = group) +
    ggplot2::labs(
      x        = outcome_label,
      title    = "Posterior predictive check by group: ECDF",
      subtitle = "Observed (dark) vs. replicated (light) empirical CDFs, within each group"
    ) +
    theme_bayes()
}


#' Posterior predictive boxplots
#'
#' Observed outcome and a handful of replicated datasets, each shown as a
#' box-and-whisker plot. Complements [plot_ppc_dens()] / [plot_ppc_ecdf()]:
#' coarser, but makes median/IQR/outlier mismatches easy to read at a
#' glance. `yrep` is downsampled hard (`n_draws` defaults to `8`) because,
#' unlike the overlay plots, every retained draw gets its own box.
#'
#' @param model A fitted Stan model.
#' @param y Numeric vector of observed outcomes.
#' @param yrep_var Name of the replicated-outcome vector family.
#' @param n_draws Number of replicated datasets to draw as boxes.
#' @param outcome_label Y-axis label.
#'
#' @return A `ggplot` object.
#' @export
#'
#' @examples
#' \dontrun{
#' plot_ppc_boxplot(fit, y = df$y, yrep_var = "y_rep")
#' }
plot_ppc_boxplot <- function(
    model,
    y,
    yrep_var      = "y_rep",
    n_draws       = 8,
    outcome_label = "outcome"
) {
  yrep <- draws_matrix_of(model, yrep_var)
  if (nrow(yrep) > n_draws) {
    yrep <- yrep[sample.int(nrow(yrep), n_draws), , drop = FALSE]
  }

  bayesplot::ppc_boxplot(y = y, yrep = yrep) +
    ggplot2::labs(
      y        = outcome_label,
      title    = "Posterior predictive check: boxplots",
      subtitle = "Observed (y) vs. a handful of replicated datasets (yrep)",
      caption  = "yrep is downsampled to n_draws replicated datasets so each box stays readable."
    ) +
    theme_bayes()
}


#' Posterior predictive boxplots, by group
#'
#' Grouped counterpart to [plot_ppc_boxplot()]. bayesplot has no built-in
#' grouped boxplot, so this stitches one [plot_ppc_boxplot()]-style panel
#' per group together with `patchwork`, the same small-multiples approach
#' [plot_ppc_stat_grid()] uses.
#'
#' @param model A fitted Stan model.
#' @param y Numeric vector of observed outcomes, length `N`.
#' @param group Integer or factor vector of group memberships, length `N`.
#' @param yrep_var Name of the replicated-outcome vector family.
#' @param n_draws Number of replicated datasets to draw as boxes, per group.
#' @param group_labels Optional character vector of group display labels.
#' @param outcome_label Y-axis label, shared across panels.
#' @param ncol Number of facet columns; `NULL` lets `patchwork` choose.
#'
#' @return A `patchwork` composite `ggplot` object.
#' @export
#'
#' @examples
#' \dontrun{
#' plot_ppc_boxplot_grouped(fit, y = df$y, group = df$treatment)
#' }
plot_ppc_boxplot_grouped <- function(
    model,
    y,
    group,
    yrep_var      = "y_rep",
    n_draws       = 8,
    group_labels  = NULL,
    outcome_label = "outcome",
    ncol          = NULL
) {
  yrep  <- draws_matrix_of(model, yrep_var)
  group <- as_group_factor(group, group_labels)

  plots <- lapply(levels(group), function(g) {
    keep <- group == g
    yrep_g <- yrep[, keep, drop = FALSE]
    if (nrow(yrep_g) > n_draws) {
      yrep_g <- yrep_g[sample.int(nrow(yrep_g), n_draws), , drop = FALSE]
    }
    bayesplot::ppc_boxplot(y = y[keep], yrep = yrep_g) +
      ggplot2::labs(title = g, y = NULL) +
      theme_bayes(base_size = 10)
  })

  patchwork::wrap_plots(plots, ncol = ncol) +
    patchwork::plot_annotation(
      title    = "Posterior predictive check by group: boxplots",
      subtitle = paste0(outcome_label, ": observed vs. replicated, within each group"),
      theme    = theme_bayes()
    )
}


#' Posterior predictive violins, by group
#'
#' Density estimate of the replicated draws within each group, shown as a
#' violin, with the observed outcome overlaid. Reveals whether the spread
#' *and* shape of the predictive distribution track the data equally well
#' across groups, not just the central tendency.
#'
#' @param model A fitted Stan model.
#' @param y Numeric vector of observed outcomes, length `N`.
#' @param group Integer or factor vector of group memberships, length `N`.
#' @param yrep_var Name of the replicated-outcome vector family.
#' @param group_labels Optional character vector of group display labels.
#' @param outcome_label Y-axis label.
#'
#' @return A `ggplot` object.
#' @export
#'
#' @examples
#' \dontrun{
#' plot_ppc_violin_grouped(fit, y = df$y, group = df$treatment)
#' }
plot_ppc_violin_grouped <- function(
    model,
    y,
    group,
    yrep_var      = "y_rep",
    group_labels  = NULL,
    outcome_label = "outcome"
) {
  yrep  <- draws_matrix_of(model, yrep_var)
  group <- as_group_factor(group, group_labels)

  bayesplot::ppc_violin_grouped(y = y, yrep = yrep, group = group) +
    ggplot2::labs(
      y        = outcome_label,
      title    = "Posterior predictive check by group: violins",
      subtitle = "Replicated distribution (violin) vs. observed outcome, within each group"
    ) +
    theme_bayes()
}


#' Posterior predictive violin
#'
#' Single-panel counterpart to [plot_ppc_violin_grouped()]: the replicated
#' draws' distribution as one violin, with the observed outcome overlaid.
#' Implemented as a thin wrapper - every observation is put in one dummy
#' group - so it shares its rendering with the grouped version rather than
#' duplicating it.
#'
#' @param model A fitted Stan model.
#' @param y Numeric vector of observed outcomes.
#' @param yrep_var Name of the replicated-outcome vector family.
#' @param outcome_label Y-axis label.
#'
#' @return A `ggplot` object.
#' @export
#'
#' @examples
#' \dontrun{
#' plot_ppc_violin(fit, y = df$y, yrep_var = "y_rep")
#' }
plot_ppc_violin <- function(
    model,
    y,
    yrep_var      = "y_rep",
    outcome_label = "outcome"
) {
  dummy_group <- factor(rep("", length(y)))

  plot_ppc_violin_grouped(
    model, y = y, group = dummy_group,
    yrep_var = yrep_var, outcome_label = outcome_label
  ) +
    ggplot2::labs(
      title    = "Posterior predictive check: violin",
      subtitle = "Replicated distribution (violin) vs. observed outcome"
    ) +
    ggplot2::theme(
      legend.position = "none",
      axis.text.x     = ggplot2::element_blank(),
      axis.ticks.x    = ggplot2::element_blank()
    )
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