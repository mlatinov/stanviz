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
# Tolerates either a vector family or a scalar parameter.
# Not exported.
draws_matrix_of_or_scalar <- function(model, var) {
  draws    <- as_draws_safe(model)
  all_vars <- posterior::variables(draws)
  if (any(grepl(paste0("^", var, "\\["), all_vars))) {
    draws_matrix_of(model, var)
  } else if (var %in% all_vars) {
    as.numeric(posterior::extract_variable(draws, var))
  } else {
    stop(sprintf("Variable '%s' not found.", var))
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
#' Posterior median + credible ribbon for an expected outcome computed across
#' a predictor grid inside Stan's `generated quantities`.
#'
#' Typical Stan-side setup: declare a grid of predictor values in `data`,
#' then a corresponding `vector[G] mu_grid` in `generated quantities`.
#' Extract those draws in R, tidy to long format with columns `x` (the grid
#' value), `.draw`, `value`, and pass to this function.
#'
#' @param draws_long A tibble with columns `x`, `.draw`, and `value`.
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
  pal <- bayes_palette()

  summ <- draws_long |>
    dplyr::group_by(.data$x) |>
    ggdist::median_qi(.data$value, .width = .width)

  ggplot2::ggplot(summ, ggplot2::aes(x = .data$x, y = .data$value)) +
    ggdist::geom_lineribbon(
      ggplot2::aes(ymin = .data$.lower, ymax = .data$.upper),
      fill = pal[["effect"]], alpha = 0.25
    ) +
    ggplot2::geom_line(colour = pal[["effect"]], linewidth = 1) +
    ggplot2::labs(
      x = x_label, y = outcome_label,
      title    = "Conditional effect",
      subtitle = "Posterior median with credible bands"
    ) +
    theme_bayes()
}