# =============================================================================
# parameter-posteriors.R
# -----------------------------------------------------------------------------
# Visualising the posterior of model PARAMETERS (coefficients, sigmas, derived
# scalars). These are the "what did we learn" plots: coefficient forests,
# halfeye densities, ROPE comparisons.
# =============================================================================


# ---- internal helper ----
# Pull named scalar parameters into a long tibble (.draw, param, value).
# Not exported.
tidy_pars <- function(model, pars, labels = NULL) {
  draws <- as_draws_safe(model, variables = pars)
  long  <- draws |>
    posterior::as_draws_df() |>
    tibble::as_tibble() |>
    dplyr::select(".draw", dplyr::all_of(pars)) |>
    tidyr::pivot_longer(-".draw", names_to = "param", values_to = "value")

  if (!is.null(labels)) {
    long <- dplyr::mutate(long, param = dplyr::recode(.data$param, !!!labels))
  }
  long
}


#' Coefficient forest plot
#'
#' Forest plot of several coefficients at once - the standard "regression
#' table as a picture". Sorted by posterior median by default.
#'
#' @param model A fitted Stan model.
#' @param pars Character vector of parameter names.
#' @param labels Optional named character vector mapping parameter names
#'   (values) to display labels (names).
#' @param sort Logical; sort by posterior median.
#' @param null_line Numeric reference line; pass `NULL` to omit.
#'
#' @return A `ggplot` object.
#' @export
#'
#' @examples
#' \dontrun{
#' plot_coef_forest(
#'   fit,
#'   pars   = c("beta_1", "beta_2", "beta_3"),
#'   labels = c("Age" = "beta_1", "Sex" = "beta_2", "Treatment" = "beta_3")
#' )
#' }
plot_coef_forest <- function(
    model,
    pars,
    labels    = NULL,
    sort      = TRUE,
    null_line = 0
) {
  long <- tidy_pars(model, pars, labels)

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
    {
      if (!is.null(null_line))
        ggplot2::geom_vline(xintercept = null_line, linetype = "dashed", colour = "grey30")
    } +
    ggdist::stat_pointinterval(
      .width = c(0.66, 0.95), point_interval = "median_qi",
      colour = pal[["effect"]]
    ) +
    ggplot2::labs(
      x = "posterior value", y = NULL,
      title    = "Coefficient posteriors",
      subtitle = "Point = median; lines = 66% and 95% credible intervals"
    ) +
    theme_bayes()
}


#' Halfeye density for one or several parameters
#'
#' Richer than a forest plot when you want to see the shape (skew, bimodality)
#' rather than just the credible interval.
#'
#' @param model A fitted Stan model.
#' @param pars Character vector of parameter names.
#' @param labels Optional named character vector mapping parameter names to
#'   display labels.
#' @param null_line Numeric reference line; pass `NULL` to omit.
#' @param facet_ncol Number of facet columns.
#'
#' @return A `ggplot` object.
#' @export
#'
#' @examples
#' \dontrun{
#' plot_param_halfeye(fit, pars = c("alpha", "beta"))
#' }
plot_param_halfeye <- function(
    model,
    pars,
    labels     = NULL,
    null_line  = 0,
    facet_ncol = 2
) {
  long <- tidy_pars(model, pars, labels)
  pal  <- bayes_palette()

  ggplot2::ggplot(long, ggplot2::aes(x = .data$value)) +
    {
      if (!is.null(null_line))
        ggplot2::geom_vline(xintercept = null_line, linetype = "dashed", colour = "grey30")
    } +
    ggdist::stat_halfeye(
      .width = c(0.66, 0.95), point_interval = "median_qi",
      fill = pal[["effect"]], alpha = 0.7
    ) +
    ggplot2::facet_wrap(~ .data$param, scales = "free", ncol = facet_ncol) +
    ggplot2::labs(
      x = "posterior value", y = NULL,
      title    = "Parameter posteriors",
      subtitle = "Median, 66% and 95% credible intervals"
    ) +
    theme_bayes()
}


#' Effect posterior vs. Region Of Practical Equivalence
#'
#' Plots a single effect's posterior with a ROPE band shaded, and reports the
#' share of posterior mass inside the ROPE in the subtitle.
#'
#' @param model A fitted Stan model.
#' @param par Single parameter name.
#' @param rope Length-2 numeric vector giving the ROPE bounds.
#' @param label Optional display label for the title.
#' @param effect_units String used in the x-axis label.
#'
#' @return A `ggplot` object.
#' @export
#'
#' @examples
#' \dontrun{
#' plot_rope(fit, par = "treatment_effect", rope = c(-0.1, 0.1))
#' }
plot_rope <- function(
    model,
    par,
    rope         = c(-0.1, 0.1),
    label        = NULL,
    effect_units = "std. units"
) {
  long <- tidy_pars(model, par)
  v    <- long$value
  p_in <- mean(v >= rope[1] & v <= rope[2])
  pal  <- bayes_palette()
  ttl  <- label %||% par

  ggplot2::ggplot(long, ggplot2::aes(x = .data$value)) +
    ggplot2::annotate(
      "rect", xmin = rope[1], xmax = rope[2],
      ymin = -Inf, ymax = Inf, fill = pal[["null"]], alpha = 0.15
    ) +
    ggplot2::geom_vline(xintercept = rope, linetype = "dotted", colour = pal[["null"]]) +
    ggdist::stat_halfeye(
      .width = c(0.66, 0.95), point_interval = "median_qi",
      fill = pal[["effect"]], alpha = 0.7
    ) +
    ggplot2::labs(
      x = paste0("effect (", effect_units, ")"), y = NULL,
      title    = paste0("Effect vs. ROPE: ", ttl),
      subtitle = sprintf(
        "%.1f%% of posterior mass inside the ROPE [%.2f, %.2f]",
        100 * p_in, rope[1], rope[2]
      ),
      caption  = "Shaded band = region of practical equivalence."
    ) +
    theme_bayes()
}


#' Probability the effect exceeds a threshold
#'
#' Replaces the binary "is it significant" with a continuous "how sure are we
#' the effect is at least this big". `favours_neg = TRUE` plots `P(effect < c)`
#' (useful when reduction is the desired direction).
#'
#' @param model A fitted Stan model.
#' @param par Single parameter name.
#' @param favours_neg Logical; if `TRUE`, plot `P(effect < c)`, else
#'   `P(effect > c)`.
#' @param label Optional display label for the title.
#' @param effect_units String used in the x-axis label.
#'
#' @return A `ggplot` object.
#' @export
#'
#' @examples
#' \dontrun{
#' plot_threshold_curve(fit, par = "treatment_effect", favours_neg = TRUE)
#' }
plot_threshold_curve <- function(
    model,
    par,
    favours_neg  = TRUE,
    label        = NULL,
    effect_units = "std. units"
) {
  v   <- tidy_pars(model, par)$value
  thr <- seq(min(v), max(v), length.out = 250)
  prob_fun <- if (favours_neg) function(c) mean(v < c) else function(c) mean(v > c)
  df  <- tibble::tibble(c = thr, prob = vapply(thr, prob_fun, numeric(1)))
  pal <- bayes_palette()
  ttl <- label %||% par

  ggplot2::ggplot(df, ggplot2::aes(x = .data$c, y = .data$prob)) +
    ggplot2::geom_area(fill = pal[["effect"]], alpha = 0.25) +
    ggplot2::geom_line(colour = pal[["effect"]], linewidth = 1) +
    ggplot2::geom_hline(yintercept = c(0.5, 0.95), linetype = "dotted", colour = "grey40") +
    ggplot2::annotate(
      "text", x = max(thr), y = 0.95, label = "95% certainty",
      vjust = -0.4, hjust = 1, size = 3, colour = "grey30"
    ) +
    ggplot2::annotate(
      "text", x = max(thr), y = 0.50, label = "50% certainty",
      vjust = -0.4, hjust = 1, size = 3, colour = "grey30"
    ) +
    ggplot2::labs(
      x = paste0("threshold c (", effect_units, ")"), y = NULL,
      title    = paste0("P(effect beyond threshold): ", ttl),
      subtitle = if (favours_neg) "P(effect < c)" else "P(effect > c)"
    ) +
    theme_bayes()
}