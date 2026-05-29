# =============================================================================
# mediation.R
# -----------------------------------------------------------------------------
# Mediation-specific figures, designed to consume the generated quantities the
# Stan mediation templates emit: direct effect, per-mediator natural indirect
# effects (NIEs), total indirect, total effect, and proportion mediated.
# =============================================================================


#' Effect decomposition: direct, indirect, total
#'
#' The headline mediation figure. Direct, total-indirect, and total effects
#' shown as stacked halfeyes so the reader sees how the total splits.
#'
#' @param model A fitted Stan model.
#' @param direct_var Name of the direct-effect scalar in `generated
#'   quantities`.
#' @param indirect_var Name of the total indirect-effect scalar.
#' @param total_var Name of the total-effect scalar.
#' @param effect_units String used in the x-axis label.
#'
#' @return A `ggplot` object.
#' @export
#'
#' @examples
#' \dontrun{
#' plot_effect_decomposition(fit)
#' }
plot_effect_decomposition <- function(
    model,
    direct_var   = "direct_effect",
    indirect_var = "total_indirect_effect",
    total_var    = "total_effect",
    effect_units = "std. units"
) {
  pars   <- c(direct_var, indirect_var, total_var)
  labels <- stats::setNames(
    c("Direct", "Total indirect", "Total"),
    pars
  )
  long <- tidy_pars(model, pars, labels)
  long <- dplyr::mutate(
    long,
    param = factor(.data$param, levels = c("Direct", "Total indirect", "Total"))
  )
  pal <- bayes_palette()

  ggplot2::ggplot(long, ggplot2::aes(x = .data$value, y = .data$param)) +
    geom_null_line(0) +
    ggdist::stat_halfeye(
      .width = c(0.66, 0.95), point_interval = "median_qi",
      fill = pal[["effect"]], alpha = 0.7
    ) +
    ggplot2::labs(
      x = paste0("effect (", effect_units, ")"), y = NULL,
      title    = "Effect decomposition",
      subtitle = "Total effect = direct + total indirect"
    ) +
    theme_bayes()
}


#' Indirect effects by mediator
#'
#' Per-mediator natural indirect effects on one forest plot.
#'
#' @param model A fitted Stan model.
#' @param nie_vars Named character vector: names become axis labels, values
#'   are Stan variable names. Example:
#'   `c("Diet" = "NIE_diet", "Activity" = "NIE_activity")`.
#' @param effect_units String used in the x-axis label.
#' @param sort Logical; sort mediators by posterior median.
#'
#' @return A `ggplot` object.
#' @export
#'
#' @examples
#' \dontrun{
#' plot_indirect_effects(
#'   fit,
#'   nie_vars = c("Diet" = "NIE_diet", "Activity" = "NIE_activity")
#' )
#' }
plot_indirect_effects <- function(
    model,
    nie_vars,
    effect_units = "std. units",
    sort         = TRUE
) {
  labels <- stats::setNames(names(nie_vars), unname(nie_vars))
  long   <- tidy_pars(model, unname(nie_vars), labels)

  if (sort) {
    ord <- long |>
      dplyr::group_by(.data$param) |>
      dplyr::summarise(m = stats::median(.data$value), .groups = "drop") |>
      dplyr::arrange(.data$m) |>
      dplyr::pull(.data$param)
    long <- dplyr::mutate(long, param = factor(.data$param, levels = ord))
  }
  pal <- bayes_palette()

  ggplot2::ggplot(long, ggplot2::aes(x = .data$value, y = .data$param)) +
    geom_null_line(0) +
    ggdist::stat_pointinterval(
      .width = c(0.66, 0.95), point_interval = "median_qi",
      colour = pal[["effect"]]
    ) +
    ggplot2::labs(
      x = paste0("indirect effect (", effect_units, ")"), y = NULL,
      title    = "Natural indirect effects by mediator",
      subtitle = "Treatment -> mediator -> outcome paths"
    ) +
    theme_bayes()
}


#' Proportion-mediated posteriors
#'
#' Proportion of the total effect transmitted through each mediator.
#' This quantity is notoriously unstable when the total effect is near zero,
#' so the function clips to a sensible display range and warns in the caption.
#'
#' @param model A fitted Stan model.
#' @param pm_vars Named character vector: names become axis labels, values
#'   are Stan variable names for proportion-mediated quantities.
#' @param clip Length-2 numeric vector giving display clipping bounds.
#'
#' @return A `ggplot` object.
#' @export
#'
#' @examples
#' \dontrun{
#' plot_proportion_mediated(
#'   fit,
#'   pm_vars = c("Diet" = "proportion_mediated_diet")
#' )
#' }
plot_proportion_mediated <- function(
    model,
    pm_vars,
    clip = c(-1, 2)
) {
  labels <- stats::setNames(names(pm_vars), unname(pm_vars))
  long   <- tidy_pars(model, unname(pm_vars), labels) |>
    dplyr::mutate(value = pmin(pmax(.data$value, clip[1]), clip[2]))
  pal <- bayes_palette()

  ggplot2::ggplot(long, ggplot2::aes(x = .data$value, y = .data$param)) +
    geom_null_line(0) +
    geom_null_line(1) +
    ggdist::stat_pointinterval(
      .width = c(0.66, 0.95), point_interval = "median_qi",
      colour = pal[["effect"]]
    ) +
    ggplot2::labs(
      x = "proportion mediated", y = NULL,
      title    = "Proportion of the total effect mediated",
      subtitle = "Reference lines at 0 (none) and 1 (fully mediated)",
      caption  = sprintf(
        "Display clipped to [%g, %g]; unstable when the total effect is near zero.",
        clip[1], clip[2]
      )
    ) +
    theme_bayes()
}


#' Mediation path diagram
#'
#' A simple generated path diagram (treatment -> mediators -> outcome) with
#' posterior-median path coefficients as edge labels. Communication aid that
#' keeps the diagram and the numbers in one place.
#'
#' @param model A fitted Stan model.
#' @param a_vars Named character vector of treatment->mediator coefficient
#'   names. Names become mediator labels.
#' @param b_vars Named character vector of mediator->outcome coefficient
#'   names, in the same order as `a_vars`.
#' @param c_var Name of the direct-effect scalar.
#' @param treatment_name Display name for the treatment node.
#' @param outcome_name Display name for the outcome node.
#'
#' @return A `ggplot` object.
#' @export
#'
#' @examples
#' \dontrun{
#' plot_mediation_paths(
#'   fit,
#'   a_vars = c("Diet" = "beta_treatment_diet"),
#'   b_vars = c("Diet" = "beta_diet_mass")
#' )
#' }
plot_mediation_paths <- function(
    model,
    a_vars,
    b_vars,
    c_var          = "direct_effect",
    treatment_name = "Treatment",
    outcome_name   = "Outcome"
) {
  med_names <- names(a_vars)
  n_med     <- length(med_names)

  med_summary <- function(v) stats::median(draws_matrix_of_or_scalar(model, v))

  a_vals <- vapply(unname(a_vars), med_summary, numeric(1))
  b_vals <- vapply(unname(b_vars), med_summary, numeric(1))
  c_val  <- med_summary(c_var)

  ys    <- if (n_med == 1) 0 else seq(1, -1, length.out = n_med)
  nodes <- tibble::tibble(
    name = c(treatment_name, med_names, outcome_name),
    x    = c(0, rep(1, n_med), 2),
    y    = c(0, ys, 0)
  )

  edges <- dplyr::bind_rows(
    tibble::tibble(
      x = 0, y = 0, xend = 1, yend = ys,
      label = sprintf("a=%.2f", a_vals)
    ),
    tibble::tibble(
      x = 1, y = ys, xend = 2, yend = 0,
      label = sprintf("b=%.2f", b_vals)
    ),
    tibble::tibble(
      x = 0, y = 0, xend = 2, yend = 0,
      label = sprintf("c'=%.2f", c_val)
    )
  )

  pal <- bayes_palette()

  ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = edges,
      ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend),
      arrow  = ggplot2::arrow(length = ggplot2::unit(0.15, "cm")),
      colour = "grey50"
    ) +
    ggplot2::geom_label(
      data = edges,
      ggplot2::aes(
        x = (.data$x + .data$xend) / 2,
        y = (.data$y + .data$yend) / 2,
        label = .data$label
      ),
      size = 3, label.size = 0, fill = "white"
    ) +
    ggplot2::geom_label(
      data = nodes,
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$name),
      fill = pal[["effect"]], colour = "white", fontface = "bold"
    ) +
    ggplot2::labs(
      title    = "Mediation path diagram",
      subtitle = "Edge labels = posterior median path coefficients"
    ) +
    ggplot2::theme_void(base_size = 12) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
}