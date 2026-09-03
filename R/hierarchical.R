# =============================================================================
# hierarchical.R
# -----------------------------------------------------------------------------
# Plots specific to hierarchical / multilevel models. Designed for the
# non-centered parameterization with 5-50 groups at a single grouping level.
#
# Typical Stan structure these functions assume:
#
#   parameters {
#     real mu;                        // grand mean
#     real<lower=0> tau;              // group-level SD
#     vector[J] z;                    // standardized group offsets
#   }
#   transformed parameters {
#     vector[J] alpha_j = mu + tau * z;   // group-level intercepts
#   }
#
# Functions take the NAMES of these quantities so they work regardless of
# your specific naming conventions.
# =============================================================================


#' Extract group-level effects into a tidy tibble
#'
#' Foundation helper for the hierarchical plotting functions. Pulls a
#' vector-valued group parameter (e.g. `alpha_j` or `z`) into a long tibble
#' with one row per group per draw.
#'
#' @param model A fitted Stan model.
#' @param group_var Name of the group-level vector in `parameters` or
#'   `transformed parameters` (e.g. `"alpha_j"`). For a nested/
#'   multi-dimensional container (e.g. `array[J] vector[K] u` or
#'   `array[K] matrix[J, R] u`), pass an R-array-style slice that leaves the
#'   group dimension free and fixes the rest, e.g. `"u[,1]"` for the first
#'   component of every group.
#' @param group_labels Optional character vector of group display names. Must
#'   have length `J` (the number of groups). If `NULL`, groups are labelled
#'   `"1"`, `"2"`, ....
#'
#' @return A [tibble::tibble] with columns `.draw`, `group` (factor),
#'   `group_id` (integer), `value`.
#' @export
#'
#' @examples
#' \dontrun{
#' extract_group_effects(fit, "alpha_j", group_labels = school_names)
#'
#' # array[J] vector[K] u -> every group's first component
#' extract_group_effects(fit, "u[,1]", group_labels = school_names)
#' }
extract_group_effects <- function(
    model,
    group_var,
    group_labels = NULL
) {
  mat <- draws_matrix_of(model, group_var)   # [n_draws, J]
  J   <- ncol(mat)

  if (is.null(group_labels)) {
    group_labels <- as.character(seq_len(J))
  } else if (length(group_labels) != J) {
    stop(sprintf(
      "`group_labels` has length %d but `%s` has %d groups.",
      length(group_labels), group_var, J
    ))
  }

  tibble::tibble(
    .draw    = rep(seq_len(nrow(mat)), times = J),
    group_id = rep(seq_len(J),         each  = nrow(mat)),
    group    = factor(rep(group_labels, each = nrow(mat)),
                      levels = group_labels),
    value    = as.numeric(mat)
  )
}


#' Caterpillar plot of group-level effects
#'
#' Forest / caterpillar plot of every group's posterior, sorted by median.
#' Two display modes:
#'
#' * `mode = "absolute"`: shows the group-level estimates themselves
#'   (e.g. `alpha_j`), with the grand mean as a reference line.
#' * `mode = "deviation"`: shows deviations from the grand mean
#'   (`alpha_j - mu`), centred on zero.
#'
#' If `mode = "deviation"` and `group_var` already contains deviations (e.g.
#' you pass `"z"`), leave `grand_mean_var = NULL`.
#'
#' @param model A fitted Stan model.
#' @param group_var Name of the group-level vector.
#' @param grand_mean_var Name of the grand-mean scalar (e.g. `"mu"`). Required
#'   when `mode = "deviation"` and `group_var` is on the absolute scale.
#' @param mode Either `"absolute"` or `"deviation"`.
#' @param group_labels Optional character vector of length `J`.
#' @param sort Logical; sort groups by posterior median.
#' @param outcome_label String used in the x-axis label.
#'
#' @return A `ggplot` object.
#' @export
#'
#' @examples
#' \dontrun{
#' # Absolute group-level intercepts
#' plot_random_effects(fit, "alpha_j", mode = "absolute",
#'                     grand_mean_var = "mu")
#'
#' # Deviations from the grand mean
#' plot_random_effects(fit, "alpha_j", grand_mean_var = "mu",
#'                     mode = "deviation")
#' }
plot_random_effects <- function(
    model,
    group_var,
    grand_mean_var = NULL,
    mode           = c("absolute", "deviation"),
    group_labels   = NULL,
    sort           = TRUE,
    outcome_label  = "group effect"
) {
  mode <- match.arg(mode)
  long <- extract_group_effects(model, group_var, group_labels)
  pal  <- bayes_palette()

  grand_mean_draws <- NULL
  if (!is.null(grand_mean_var)) {
    grand_mean_draws <- as.numeric(
      posterior::extract_variable(as_draws_safe(model), grand_mean_var)
    )
    grand_mean_median <- stats::median(grand_mean_draws)
  }

  if (mode == "deviation") {
    if (is.null(grand_mean_var)) {
      # Assume group_var already contains deviations
      ref_line <- 0
    } else {
      # Subtract grand mean draw-by-draw
      gm_per_draw <- grand_mean_draws[long$.draw]
      long <- dplyr::mutate(long, value = .data$value - gm_per_draw)
      ref_line <- 0
    }
    xlab <- paste0(outcome_label, " (deviation from grand mean)")
  } else {
    ref_line <- if (!is.null(grand_mean_var)) grand_mean_median else NULL
    xlab <- outcome_label
  }

  if (sort) {
    ord <- long |>
      dplyr::group_by(.data$group) |>
      dplyr::summarise(m = stats::median(.data$value), .groups = "drop") |>
      dplyr::arrange(.data$m) |>
      dplyr::pull(.data$group)
    long <- dplyr::mutate(long, group = factor(.data$group, levels = ord))
  }

  p <- ggplot2::ggplot(long, ggplot2::aes(x = .data$value, y = .data$group))
  if (!is.null(ref_line)) {
    p <- p + ggplot2::geom_vline(
      xintercept = ref_line, linetype = "dashed", colour = "grey30"
    )
  }
  p +
    ggdist::stat_pointinterval(
      .width = c(0.66, 0.95), point_interval = "median_qi",
      colour = pal[["effect"]]
    ) +
    ggplot2::labs(
      x = xlab, y = NULL,
      title    = "Group-level effects",
      subtitle = if (mode == "absolute")
        "Posterior medians and credible intervals; dashed line = grand mean"
      else
        "Deviations from grand mean; dashed line = zero"
    ) +
    theme_bayes()
}


#' Shrinkage plot: raw vs. partially-pooled group means
#'
#' Visualises what partial pooling does. Plots raw per-group means against
#' partially-pooled posterior medians, with arrows showing how much each group
#' shrinks toward the grand mean.
#'
#' @param model A fitted Stan model.
#' @param group_var Name of the group-level vector (absolute scale).
#' @param raw_means Numeric vector of length `J`: the raw per-group sample
#'   means computed from your data.
#' @param grand_mean_var Optional name of the grand-mean scalar; if supplied,
#'   the grand mean is drawn as a horizontal reference line.
#' @param group_labels Optional character vector of length `J`.
#' @param outcome_label Axis label.
#'
#' @return A `ggplot` object.
#' @export
#'
#' @examples
#' \dontrun{
#' raw <- tapply(df$y, df$group, mean)
#' plot_shrinkage(fit, "alpha_j", raw_means = raw, grand_mean_var = "mu")
#' }
plot_shrinkage <- function(
    model,
    group_var,
    raw_means,
    grand_mean_var = NULL,
    group_labels   = NULL,
    outcome_label  = "outcome"
) {
  long <- extract_group_effects(model, group_var, group_labels)
  J    <- length(raw_means)
  if (J != length(unique(long$group))) {
    stop(sprintf(
      "`raw_means` has length %d but `%s` has %d groups.",
      J, group_var, length(unique(long$group))
    ))
  }

  pooled <- long |>
    dplyr::group_by(.data$group_id, .data$group) |>
    dplyr::summarise(pooled = stats::median(.data$value), .groups = "drop")

  df <- pooled |>
    dplyr::mutate(raw = raw_means[.data$group_id])

  pal <- bayes_palette()
  gm  <- if (!is.null(grand_mean_var)) {
    stats::median(as.numeric(
      posterior::extract_variable(as_draws_safe(model), grand_mean_var)
    ))
  } else NULL

  p <- ggplot2::ggplot(df) +
    ggplot2::geom_segment(
      ggplot2::aes(
        x = .data$raw, xend = .data$pooled,
        y = .data$group, yend = .data$group
      ),
      arrow  = ggplot2::arrow(length = ggplot2::unit(0.12, "cm"), type = "closed"),
      colour = "grey60"
    ) +
    ggplot2::geom_point(
      ggplot2::aes(x = .data$raw, y = .data$group),
      colour = pal[["control"]], size = 2.5
    ) +
    ggplot2::geom_point(
      ggplot2::aes(x = .data$pooled, y = .data$group),
      colour = pal[["effect"]], size = 2.5
    )

  if (!is.null(gm)) {
    p <- p + ggplot2::geom_vline(xintercept = gm, linetype = "dashed", colour = "grey30")
  }

  p +
    ggplot2::labs(
      x = outcome_label, y = NULL,
      title    = "Shrinkage from partial pooling",
      subtitle = "Open dot = raw group mean; filled dot = partial-pooled estimate",
      caption  = "Arrow length = amount of shrinkage toward the grand mean."
    ) +
    theme_bayes()
}


#' Variance decomposition across hierarchical levels
#'
#' Posterior of the proportion of variance attributable to each level.
#' For a typical two-level model with `tau` (group-level SD) and `sigma`
#' (residual SD), this shows what fraction of variance lives between vs.
#' within groups.
#'
#' @param model A fitted Stan model.
#' @param sd_vars Named character vector of standard-deviation parameters.
#'   Names become level labels in the plot. Example:
#'   `c("Between groups" = "tau", "Within groups" = "sigma")`.
#'
#' @return A `ggplot` object.
#' @export
#'
#' @examples
#' \dontrun{
#' plot_variance_decomposition(
#'   fit,
#'   sd_vars = c("Between schools" = "tau", "Within schools" = "sigma")
#' )
#' }
plot_variance_decomposition <- function(
    model,
    sd_vars
) {
  draws <- as_draws_safe(model, variables = unname(sd_vars))
  df    <- posterior::as_draws_df(draws) |>
    tibble::as_tibble() |>
    dplyr::select(dplyr::all_of(unname(sd_vars)))

  # Variances per draw, then proportions per draw
  variances  <- df^2
  totals     <- rowSums(variances)
  props      <- as.data.frame(variances / totals)
  names(props) <- names(sd_vars)

  long <- props |>
    tibble::as_tibble() |>
    tidyr::pivot_longer(
      dplyr::everything(),
      names_to  = "level",
      values_to = "proportion"
    ) |>
    dplyr::mutate(level = factor(.data$level, levels = names(sd_vars)))

  pal <- bayes_palette()

  ggplot2::ggplot(long, ggplot2::aes(x = .data$proportion, y = .data$level)) +
    ggdist::stat_halfeye(
      .width = c(0.66, 0.95), point_interval = "median_qi",
      fill = pal[["effect"]], alpha = 0.7
    ) +
    ggplot2::scale_x_continuous(limits = c(0, 1), labels = scales::percent) +
    ggplot2::labs(
      x = "proportion of total variance", y = NULL,
      title    = "Variance decomposition",
      subtitle = "Posterior share of variance at each level"
    ) +
    theme_bayes()
}


#' Intraclass correlation (ICC) posterior
#'
#' Specialised version of variance decomposition: the posterior of the
#' between-group share of variance — the standard two-level ICC.
#' Defined as `tau^2 / (tau^2 + sigma^2)`.
#'
#' @param model A fitted Stan model.
#' @param tau_var Name of the group-level SD.
#' @param sigma_var Name of the residual / within-group SD.
#'
#' @return A `ggplot` object.
#' @export
#'
#' @examples
#' \dontrun{
#' plot_icc(fit, tau_var = "tau", sigma_var = "sigma")
#' }
plot_icc <- function(
    model,
    tau_var   = "tau",
    sigma_var = "sigma"
) {
  draws <- as_draws_safe(model, variables = c(tau_var, sigma_var))
  tau   <- as.numeric(posterior::extract_variable(draws, tau_var))
  sigma <- as.numeric(posterior::extract_variable(draws, sigma_var))
  icc   <- tau^2 / (tau^2 + sigma^2)

  df  <- tibble::tibble(icc = icc)
  pal <- bayes_palette()

  ggplot2::ggplot(df, ggplot2::aes(x = .data$icc)) +
    ggdist::stat_halfeye(
      .width = c(0.66, 0.95), point_interval = "median_qi",
      fill = pal[["effect"]], alpha = 0.7
    ) +
    ggplot2::scale_x_continuous(limits = c(0, 1), labels = scales::percent) +
    ggplot2::labs(
      x = "ICC", y = NULL,
      title    = "Intraclass correlation",
      subtitle = "Posterior share of variance attributable to groups",
      caption  = "ICC = tau^2 / (tau^2 + sigma^2)"
    ) +
    theme_bayes()
}


#' Posterior predictive check faceted by group
#'
#' Faceted PPC density: one panel per group. Reveals whether the model fits
#' equally well across groups, or only well on average.
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
#' plot_ppc_by_group(fit, y = df$y, group = df$school_id)
#' }
plot_ppc_by_group <- function(
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

  bayesplot::ppc_dens_overlay_grouped(
    y     = y,
    yrep  = yrep,
    group = group
  ) +
    ggplot2::labs(
      x        = outcome_label,
      title    = "Posterior predictive check by group",
      subtitle = "Observed (dark) vs. replicated (light) within each group"
    ) +
    theme_bayes()
}


#' Funnel diagnostic for hierarchical parameterization
#'
#' Scatter of group-level estimates against `log(tau)`, with divergent
#' transitions highlighted. The classic "funnel" pattern, with divergences
#' clustered at small `tau`, signals that you should switch from centered to
#' non-centered parameterization (or vice versa, in rare cases).
#'
#' Requires a CmdStanR or rstan fit (needs sampler diagnostics).
#'
#' @param model A `CmdStanMCMC` or `stanfit` object.
#' @param group_var Name of the group-level vector (or `z` if non-centered).
#'   For a nested/multi-dimensional container (e.g. `array[J] vector[K] z`),
#'   pass the base name and use `which_group` to pick one element.
#' @param tau_var Name of the group-level SD.
#' @param which_group Which group index to plot on the y-axis. A single
#'   integer for a plain vector (the default, `1`), or an integer vector
#'   giving one index per dimension for a nested/multi-dimensional
#'   `group_var` (e.g. `c(1, 2)` for `group_var[1,2]`).
#'
#' @return A `ggplot` object.
#' @export
#'
#' @examples
#' \dontrun{
#' plot_funnel_diagnostic(fit, group_var = "z", tau_var = "tau")
#'
#' # array[J] vector[K] z -> group 1, component 2
#' plot_funnel_diagnostic(fit, group_var = "z", which_group = c(1, 2))
#' }
plot_funnel_diagnostic <- function(
    model,
    group_var,
    tau_var      = "tau",
    which_group  = 1
) {
  if (!inherits(model, "CmdStanMCMC") && !inherits(model, "stanfit")) {
    stop("plot_funnel_diagnostic() needs a CmdStanR or rstan fit.")
  }

  draws <- as_draws_safe(model)
  all_vars <- posterior::variables(draws)
  which_group_str <- paste(which_group, collapse = ",")
  group_name <- paste0(group_var, "[", which_group_str, "]")
  if (!(group_name %in% all_vars)) {
    stop(sprintf("Could not find '%s' in the draws.", group_name))
  }

  group_draws <- as.numeric(posterior::extract_variable(draws, group_name))
  tau_draws   <- as.numeric(posterior::extract_variable(draws, tau_var))

  np <- bayesplot::nuts_params(model)
  div <- np |>
    dplyr::filter(.data$Parameter == "divergent__") |>
    dplyr::pull(.data$Value)

  df <- tibble::tibble(
    group_value = group_draws,
    log_tau     = log(tau_draws),
    divergent   = as.logical(div)
  )
  pal <- bayes_palette()

  ggplot2::ggplot(df, ggplot2::aes(x = .data$group_value, y = .data$log_tau)) +
    ggplot2::geom_point(
      data = subset(df, !df$divergent),
      colour = pal[["control"]], alpha = 0.4, size = 0.8
    ) +
    ggplot2::geom_point(
      data = subset(df, df$divergent),
      colour = pal[["null"]], alpha = 0.9, size = 1.2
    ) +
    ggplot2::labs(
      x = sprintf("%s[%s]", group_var, which_group_str),
      y = sprintf("log(%s)", tau_var),
      title    = "Funnel diagnostic",
      subtitle = "Red points = divergent transitions",
      caption  = "Divergences concentrated at small tau suggest the parameterization is fighting the geometry."
    ) +
    theme_bayes()
}


#' Random intercept-slope correlation
#'
#' For models with correlated random intercepts and slopes. Shows the joint
#' distribution of per-group `(intercept, slope)` posterior medians plus the
#' posterior of the correlation parameter as a marginal halfeye.
#'
#' @param model A fitted Stan model.
#' @param intercept_var Name of the random-intercept vector.
#' @param slope_var Name of the random-slope vector.
#' @param rho_var Name of the correlation scalar (e.g. `"rho"` or
#'   `"L_Omega[2,1]"` — for the latter, extract first and use the scalar form).
#' @param group_labels Optional character vector of length `J`.
#'
#' @return A patchwork composite ggplot.
#' @export
#'
#' @examples
#' \dontrun{
#' plot_random_correlation(fit, "alpha_j", "beta_j", rho_var = "rho")
#' }
plot_random_correlation <- function(
    model,
    intercept_var,
    slope_var,
    rho_var,
    group_labels = NULL
) {
  a <- extract_group_effects(model, intercept_var, group_labels) |>
    dplyr::group_by(.data$group) |>
    dplyr::summarise(intercept = stats::median(.data$value), .groups = "drop")
  b <- extract_group_effects(model, slope_var, group_labels) |>
    dplyr::group_by(.data$group) |>
    dplyr::summarise(slope = stats::median(.data$value), .groups = "drop")
  joint <- dplyr::inner_join(a, b, by = "group")

  rho_draws <- as.numeric(
    posterior::extract_variable(as_draws_safe(model), rho_var)
  )
  rho_df <- tibble::tibble(rho = rho_draws)

  pal <- bayes_palette()

  p_scatter <- ggplot2::ggplot(joint, ggplot2::aes(x = .data$intercept, y = .data$slope)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
    ggplot2::geom_point(colour = pal[["effect"]], size = 2.5, alpha = 0.8) +
    ggplot2::labs(
      x = "random intercept (posterior median)",
      y = "random slope (posterior median)",
      title = "Group-level intercepts vs. slopes"
    ) +
    theme_bayes()

  p_rho <- ggplot2::ggplot(rho_df, ggplot2::aes(x = .data$rho)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey30") +
    ggdist::stat_halfeye(
      .width = c(0.66, 0.95), point_interval = "median_qi",
      fill = pal[["effect"]], alpha = 0.7
    ) +
    ggplot2::scale_x_continuous(limits = c(-1, 1)) +
    ggplot2::labs(
      x = "correlation (rho)", y = NULL,
      title = "Correlation posterior"
    ) +
    theme_bayes()

  patchwork::wrap_plots(p_scatter, p_rho, ncol = 1, heights = c(2, 1))
}