# =============================================================================
# prior-sensitivity.R
# -----------------------------------------------------------------------------
# Prior / likelihood sensitivity analysis via power-scaling (priorsense).
# These functions require the priorsense package (Suggests). The Stan model
# must expose pointwise `log_lik` and total `lprior` in generated quantities
# for power-scaling to work.
# =============================================================================


#' Power-scaling sensitivity summary
#'
#' Tidy table of prior/likelihood sensitivity diagnostics for selected
#' parameters.
#'
#' @param model A fitted Stan model (CmdStanR, rstan, or brms).
#' @param variable Optional character vector restricting which parameters to
#'   diagnose.
#'
#' @return A [tibble::tibble] with one row per parameter and columns for
#'   prior sensitivity, likelihood sensitivity, and a diagnosis flag.
#' @export
#'
#' @examples
#' \dontrun{
#' powerscale_summary(fit, variable = c("alpha", "beta"))
#' }
powerscale_summary <- function(model, variable = NULL) {
  rlang::check_installed("priorsense", reason = "for power-scaling diagnostics.")
  s <- priorsense::powerscale_sensitivity(model, variable = variable)
  tibble::as_tibble(s$sensitivity)
}


#' Power-scaling density plot
#'
#' Shows how the posterior of each parameter shifts as the prior and
#' likelihood are scaled up/down. A near-flat posterior across scaling means
#' the conclusion is robust; strong movement under prior-scaling means the
#' conclusion is prior-driven.
#'
#' @param model A fitted Stan model.
#' @param variable Character vector of parameter names to plot.
#' @param plot_type Currently only `"dens"` is wired through; reserved for
#'   future variants.
#'
#' @return A `ggplot` object.
#' @export
#'
#' @examples
#' \dontrun{
#' plot_powerscale(fit, variable = c("alpha", "beta"))
#' }
plot_powerscale <- function(
    model,
    variable,
    plot_type = c("dens", "quantities", "ecdf")
) {
  rlang::check_installed("priorsense", reason = "for power-scaling diagnostics.")
  plot_type <- match.arg(plot_type)

  p <- priorsense::powerscale_plot_dens(
    priorsense::powerscale_sequence(model),
    variable = variable
  )

  if (inherits(p, "ggplot")) {
    p <- p +
      ggplot2::labs(
        title    = "Prior / likelihood power-scaling",
        subtitle = "Posterior movement as prior and likelihood are perturbed"
      ) +
      theme_bayes()
  }
  p
}


#' Power-scaling summaries plot
#'
#' Compact alternative to overlaid densities: plots how posterior summaries
#' (mean, sd) of chosen quantities track the power-scaling factor.
#'
#' @param model A fitted Stan model.
#' @param variable Character vector of parameter names.
#'
#' @return A `ggplot` object.
#' @export
#'
#' @examples
#' \dontrun{
#' plot_powerscale_quantities(fit, variable = c("alpha", "beta"))
#' }
plot_powerscale_quantities <- function(
    model,
    variable
) {
  rlang::check_installed("priorsense", reason = "for power-scaling diagnostics.")
  seq_obj <- priorsense::powerscale_sequence(model)
  p <- priorsense::powerscale_plot_quantities(seq_obj, variable = variable)
  if (inherits(p, "ggplot")) {
    p <- p +
      ggplot2::labs(
        title    = "Power-scaling: posterior quantities",
        subtitle = "How summaries respond to perturbing prior vs. likelihood"
      ) +
      theme_bayes()
  }
  p
}