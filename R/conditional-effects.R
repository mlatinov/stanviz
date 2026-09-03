# =============================================================================
# conditional-effects.R
# -----------------------------------------------------------------------------
# Effects on the OUTCOME scale: counterfactual marginal means, average
# treatment effects, conditional/partial effects of a predictor.
#
# This file works entirely off raw Stan draws (no marginaleffects dependency).
# The typical workflow: compute arm-specific predictions inside Stan's
# `generated quantities`, then pass their names here.
# =============================================================================


# ---- internal helper ----
# Tolerates either a vector family (including nested/multi-dimensional
# containers and R-array-style slices, e.g. "samples_combined[1,]") or a
# scalar parameter. Not exported.
draws_matrix_of_or_scalar <- function(model, var) {
  draws    <- as_draws_safe(model)
  all_vars <- posterior::variables(draws)
  res      <- resolve_stan_var(all_vars, var)
  if (res$is_scalar) {
    as.numeric(posterior::extract_variable(draws, res$hit_vars))
  } else {
    m <- posterior::subset_draws(draws, variable = res$hit_vars) |>
      posterior::as_draws_matrix()
    m[, res$hit_vars, drop = FALSE]
  }
}


#' Extract ATE draws from arm-specific Stan predictions
#'
#' Given two generated-quantity variables representing `E[Y | treated]` and
#' `E[Y | control]`, returns a tibble of per-draw counterfactual means and
#' their difference (the ATE).
#'
#' If either variable is per-observation (length `N`), it is averaged within
#' each draw to give the marginal mean for that arm.
#'
#' Typical Stan-side setup:
#' \preformatted{
#' generated quantities \{
#'   vector[N] mu_treated = alpha + beta * 1 + ...;
#'   vector[N] mu_control = alpha + beta * 0 + ...;
#' \}
#' }
#'
#' @param model A fitted Stan model.
#' @param mu_treated_var Name of the treated-arm expected outcome.
#' @param mu_control_var Name of the control-arm expected outcome.
#'
#' @return A [tibble::tibble] with columns `.draw`, `mu_control`,
#'   `mu_treated`, `ATE`.
#' @export
#'
#' @examples
#' \dontrun{
#' ate <- extract_ate_from_draws(fit, "mu_treated", "mu_control")
#' plot_ate_posterior(ate)
#' }
extract_ate_from_draws <- function(
    model,
    mu_treated_var,
    mu_control_var
) {
  to_marginal <- function(var) {
    m <- draws_matrix_of_or_scalar(model, var)
    if (is.matrix(m) && ncol(m) > 1) rowMeans(m) else as.numeric(m)
  }
  mu1 <- to_marginal(mu_treated_var)
  mu0 <- to_marginal(mu_control_var)

  tibble::tibble(
    .draw      = seq_along(mu1),
    mu_control = mu0,
    mu_treated = mu1,
    ATE        = mu1 - mu0
  )
}


#' ATE posterior halfeye
#'
#' Visualises the posterior of the average treatment effect.
#'
#' @param ate_draws A tibble with an `ATE` column (typically from
#'   [extract_ate_from_draws()]).
#' @param outcome_label Label for the title; describes the outcome measured.
#' @param effect_color Fill colour for the halfeye.
#'
#' @return A `ggplot` object.
#' @export
#'
#' @examples
#' \dontrun{
#' ate <- extract_ate_from_draws(fit, "mu_treated", "mu_control")
#' plot_ate_posterior(ate, outcome_label = "mass change")
#' }
plot_ate_posterior <- function(
    ate_draws,
    outcome_label = "outcome",
    effect_color  = "red4"
) {
  pal <- bayes_palette(effect = effect_color, treatment = effect_color)

  ggplot2::ggplot(ate_draws, ggplot2::aes(x = .data$ATE)) +
    ggdist::stat_halfeye(
      .width = c(0.66, 0.95), point_interval = "median_qi",
      fill = pal[["effect"]], alpha = 0.7
    ) +
    geom_null_line(0) +
    ggplot2::labs(
      x = "ATE", y = NULL,
      title    = paste0("Posterior average treatment effect: ", outcome_label),
      subtitle = "Median, 66% and 95% credible intervals",
      caption  = "Dashed line = no effect."
    ) +
    theme_bayes()
}


#' Counterfactual marginal means
#'
#' Side-by-side counterfactual marginal means under each arm.
#'
#' Expects a long tibble of predictive draws with a treatment indicator
#' column (`0`/`1`) and a value column called `draw` (the per-draw mean for
#' that arm). You can build this from [extract_ate_from_draws()] output with
#' `tidyr::pivot_longer()`.
#'
#' @param pred_draws Long tibble of predictive draws.
#' @param treatment_var Name of the treatment indicator column.
#' @param outcome_label Y-axis label.
#' @param control_color Colour for the control arm.
#' @param treatment_color Colour for the treated arm.
#'
#' @return A `ggplot` object.
#' @export
#'
#' @examples
#' \dontrun{
#' ate_long <- extract_ate_from_draws(fit, "mu_treated", "mu_control") |>
#'   tidyr::pivot_longer(c(mu_control, mu_treated),
#'                       names_to  = "treatment",
#'                       values_to = "draw") |>
#'   dplyr::mutate(treatment = ifelse(treatment == "mu_treated", 1, 0))
#' plot_counterfactual_means(ate_long)
#' }
plot_counterfactual_means <- function(
    pred_draws,
    treatment_var   = "treatment",
    outcome_label   = "outcome",
    control_color   = "lightblue4",
    treatment_color = "red4"
) {
  pal <- bayes_palette(control = control_color, treatment = treatment_color)

  ggplot2::ggplot(
    pred_draws,
    ggplot2::aes(
      x    = factor(.data[[treatment_var]], labels = c("Untreated", "Treated")),
      y    = .data$draw,
      fill = factor(.data[[treatment_var]])
    )
  ) +
    ggdist::stat_halfeye(.width = c(0.66, 0.95), point_interval = "median_qi") +
    ggplot2::scale_fill_manual(values = unname(pal[c("control", "treatment")])) +
    ggplot2::labs(
      x = NULL, y = outcome_label,
      title    = "Counterfactual marginal means",
      subtitle = "Average outcome under each treatment arm"
    ) +
    theme_bayes() +
    ggplot2::theme(legend.position = "none")
}


#' Conditional effect from a grid of generated-quantity draws
#'
#' Posterior median + credible ribbon(s) for an expected outcome computed
#' across a predictor grid inside Stan's `generated quantities`. If
#' `draws_long` carries a `group` column, one coloured curve/ribbon is drawn
#' per group instead of a single one - this is what powers
#' [plot_effect_curve_grouped()], and lets you build custom grouped effect
#' plots without going through a Stan variable name at all.
#'
#' Typical Stan-side setup: declare a grid of predictor values in `data`,
#' then a corresponding `vector[G] mu_grid` in `generated quantities`.
#' Extract those draws in R, tidy to long format with columns `x` (the grid
#' value), `.draw`, `value` (and optionally `group`), and pass to this
#' function. [plot_effect_curve()] / [plot_effect_curve_grouped()] do this
#' extraction for you directly off a fitted model.
#'
#' @param draws_long A tibble with columns `x`, `.draw`, `value`, and
#'   optionally `group` (any type; coerced to factor).
#' @param outcome_label Y-axis label.
#' @param x_label X-axis label.
#' @param .width Numeric vector of credible-interval widths.
#'
#' @return A `ggplot` object.
#' @export
#'
#' @examples
#' \dontrun{
#' plot_conditional_effect_draws(draws_long, x_label = "age")
#' }
plot_conditional_effect_draws <- function(
    draws_long,
    outcome_label = "outcome",
    x_label       = "predictor",
    .width        = c(0.66, 0.95)
) {
  pal       <- bayes_palette()
  has_group <- "group" %in% names(draws_long)

  if (has_group) {
    draws_long <- dplyr::mutate(draws_long, group = as.factor(.data$group))
    summ <- draws_long |>
      dplyr::group_by(.data$group, .data$x) |>
      ggdist::median_qi(.data$value, .width = .width)
  } else {
    summ <- draws_long |>
      dplyr::group_by(.data$x) |>
      ggdist::median_qi(.data$value, .width = .width)
  }

  base_aes <- if (has_group) {
    ggplot2::aes(
      x = .data$x, y = .data$value,
      colour = .data$group, fill = .data$group
    )
  } else {
    ggplot2::aes(x = .data$x, y = .data$value)
  }

  ggplot2::ggplot(summ, base_aes) +
    {
      if (has_group) {
        ggdist::geom_lineribbon(
          ggplot2::aes(ymin = .data$.lower, ymax = .data$.upper, group = interaction(.data$group, .data$.width)),
          alpha = 0.25
        )
      } else {
        ggdist::geom_lineribbon(
          ggplot2::aes(ymin = .data$.lower, ymax = .data$.upper, group = .data$.width),
          fill = pal[["effect"]], alpha = 0.25
        )
      }
    } +
    {
      if (has_group) {
        ggplot2::geom_line(linewidth = 1)
      } else {
        ggplot2::geom_line(colour = pal[["effect"]], linewidth = 1)
      }
    } +
    ggplot2::labs(
      x = x_label, y = outcome_label, colour = NULL, fill = NULL,
      title    = "Conditional effect",
      subtitle = "Posterior median with credible bands"
    ) +
    theme_bayes()
}


# ---- internal helper ----
# Builds the (.draw, x, value[, group]) long tibble that both
# plot_effect_curve() and plot_effect_curve_grouped() feed to
# plot_conditional_effect_draws(). Reuses gather_indexed() rather than
# re-implementing Stan variable extraction. Not exported.
effect_curve_draws <- function(model, var, x, group = NULL) {
  dl <- gather_indexed(model, var)
  if (!("index" %in% names(dl))) {
    stop(sprintf(
      "'%s' resolves to more than one free dimension; slice it down to a single grid dimension (e.g. '%s[1,]').",
      var, sub("\\[.*", "", var)
    ), call. = FALSE)
  }
  n_grid <- length(unique(dl$index))
  if (n_grid != length(x)) {
    stop(sprintf(
      "'%s' has %d grid point(s) but `x` has length %d.",
      var, n_grid, length(x)
    ), call. = FALSE)
  }
  dl$x <- x[dl$index]
  if (!is.null(group)) dl$group <- group
  dl
}


#' Posterior effect / response curve
#'
#' General-purpose posterior effect curve: how a model-implied outcome
#' changes across a predictor grid, with posterior uncertainty. Unlike a
#' smoothed fit through the raw observations, the curve and ribbon(s) come
#' straight from posterior draws of a Stan-computed prediction grid - e.g. a
#' `vector[G] mu_grid` evaluated at `G` predictor values in `generated
#' quantities`.
#'
#' Built entirely on existing stanviz machinery: [gather_indexed()] pulls
#' the grid draws off the model (nested/multi-dimensional grids can be
#' reached with an index slice, e.g. `"mu_grid[1,]"` - see
#' [gather_indexed()]), and [plot_conditional_effect_draws()] does the
#' summarising and drawing.
#'
#' @param model A fitted Stan model (anything [as_draws_safe()] accepts).
#' @param effect_var Name of the predictor-grid vector in `generated
#'   quantities`, e.g. `"mu_grid"`. May be an R-array-style slice of a
#'   nested/multi-dimensional container, e.g. `"mu_grid[1,]"`.
#' @param x Numeric vector of predictor grid values, the same length and in
#'   the same order as `effect_var`'s indices.
#' @param x_obs,y_obs Optional numeric vectors of observed predictor/outcome
#'   pairs, overlaid as points for reference. `NULL` (the default) omits
#'   them. Note these are raw data, not part of the posterior curve itself.
#' @param .width Numeric vector of credible-interval widths for the
#'   ribbon(s); pass two values (the default) for a narrower and a wider
#'   band.
#' @param outcome_label Y-axis label.
#' @param x_label X-axis label.
#'
#' @return A `ggplot` object.
#' @export
#'
#' @examples
#' \dontrun{
#' # e.g. an SOD -> H2O2 dose-response curve
#' plot_effect_curve(
#'   fit, "mu_grid", x = sod_grid,
#'   x_obs = df$sod, y_obs = df$h2o2,
#'   outcome_label = "H2O2", x_label = "SOD"
#' )
#' }
plot_effect_curve <- function(
    model,
    effect_var,
    x,
    x_obs         = NULL,
    y_obs         = NULL,
    .width        = c(0.66, 0.95),
    outcome_label = "outcome",
    x_label       = "predictor"
) {
  draws_long <- effect_curve_draws(model, effect_var, x)

  p <- plot_conditional_effect_draws(
    draws_long,
    outcome_label = outcome_label,
    x_label       = x_label,
    .width        = .width
  ) +
    ggplot2::labs(
      title    = "Posterior effect curve",
      subtitle = "Posterior median with credible band(s), from model draws"
    )

  if (!is.null(x_obs) && !is.null(y_obs)) {
    obs <- tibble::tibble(x = x_obs, y = y_obs)
    p <- p + ggplot2::geom_point(
      data = obs, ggplot2::aes(x = .data$x, y = .data$y),
      inherit.aes = FALSE, colour = "grey20", alpha = 0.5, size = 1.5
    )
  }
  p
}


#' Posterior effect / response curves, by group
#'
#' Grouped counterpart to [plot_effect_curve()]: overlays one posterior
#' effect curve per group in a single panel - e.g. separate treatment arms
#' or biological conditions - so their shapes can be compared directly.
#' Built on the exact same machinery as [plot_effect_curve()]: one
#' [gather_indexed()] pull per group, stacked into a single long tibble and
#' handed to [plot_conditional_effect_draws()], which draws one
#' colour-coded curve/ribbon per group.
#'
#' @param model A fitted Stan model.
#' @param effect_vars Named character vector: names become group labels,
#'   values are the Stan grid-variable name (or index slice) for that
#'   group's predictions. Example for a nested `array[2] vector[G] mu_grid`:
#'   `c("Treated" = "mu_grid[1,]", "Control" = "mu_grid[2,]")`. Or, for two
#'   separate generated quantities: `c("Treated" = "mu_grid_t", "Control" =
#'   "mu_grid_c")`.
#' @param x Numeric vector of predictor grid values, shared across groups,
#'   the same length and order as each variable's indices.
#' @param obs_data Optional data frame of observed points, overlaid as
#'   points coloured by group. `NULL` (the default) omits them.
#' @param x_col,y_col,group_col Column names to use within `obs_data`.
#' @param .width Numeric vector of credible-interval widths for the
#'   ribbon(s).
#' @param outcome_label Y-axis label.
#' @param x_label X-axis label.
#'
#' @return A `ggplot` object.
#' @export
#'
#' @examples
#' \dontrun{
#' plot_effect_curve_grouped(
#'   fit,
#'   effect_vars = c("Treated" = "mu_grid[1,]", "Control" = "mu_grid[2,]"),
#'   x = sod_grid,
#'   obs_data = df, x_col = "sod", y_col = "h2o2", group_col = "treatment"
#' )
#' }
plot_effect_curve_grouped <- function(
    model,
    effect_vars,
    x,
    obs_data      = NULL,
    x_col         = "x",
    y_col         = "y",
    group_col     = "group",
    .width        = c(0.66, 0.95),
    outcome_label = "outcome",
    x_label       = "predictor"
) {
  if (is.null(names(effect_vars)) || any(!nzchar(names(effect_vars)))) {
    stop("`effect_vars` must be a named character vector; names become group labels.", call. = FALSE)
  }

  draws_long <- dplyr::bind_rows(lapply(names(effect_vars), function(g) {
    effect_curve_draws(model, effect_vars[[g]], x, group = g)
  }))
  draws_long <- dplyr::mutate(
    draws_long,
    group = factor(.data$group, levels = names(effect_vars))
  )

  p <- plot_conditional_effect_draws(
    draws_long,
    outcome_label = outcome_label,
    x_label       = x_label,
    .width        = .width
  ) +
    ggplot2::labs(
      title    = "Posterior effect curves by group",
      subtitle = "Posterior median with credible band(s), from model draws"
    )

  if (!is.null(obs_data)) {
    obs <- tibble::tibble(
      x     = obs_data[[x_col]],
      y     = obs_data[[y_col]],
      group = factor(obs_data[[group_col]], levels = names(effect_vars))
    )
    p <- p + ggplot2::geom_point(
      data = obs,
      ggplot2::aes(x = .data$x, y = .data$y, colour = .data$group),
      inherit.aes = FALSE, alpha = 0.5, size = 1.5
    )
  }
  p
}