# =============================================================================
# model-comparison.R
# -----------------------------------------------------------------------------
# Out-of-sample fit and model comparison via PSIS-LOO. These rely on the
# pointwise `log_lik` vector emitted by the Stan templates. The `loo` package
# is in Suggests.
# =============================================================================


#' Compute PSIS-LOO from a fitted Stan model
#'
#' Convenience wrapper: pulls `log_lik` out of a fit and runs PSIS-LOO with
#' the `r_eff` correction (which uses MCMC effective sample sizes for honest
#' standard errors).
#'
#' @param model A fitted Stan model.
#' @param loglik_var Name of the pointwise log-likelihood vector family in
#'   `generated quantities`.
#'
#' @return A `loo` object.
#' @export
#'
#' @examples
#' \dontrun{
#' loo_obj <- compute_loo(fit)
#' plot_loo_pareto_k(loo_obj)
#' }
compute_loo <- function(model, loglik_var = "log_lik") {
  rlang::check_installed("loo", reason = "for PSIS-LOO computation.")
  ll <- draws_matrix_of(model, loglik_var)
  draws <- as_draws_safe(model)
  r_eff <- loo::relative_eff(
    exp(ll),
    chain_id = rep(
      seq_len(posterior::nchains(draws)),
      each = posterior::ndraws(draws) / posterior::nchains(draws)
    )
  )
  loo::loo(ll, r_eff = r_eff)
}


#' Pareto-k diagnostic plot
#'
#' Each point is one observation; high `k` (above 0.7) means the importance
#' sampling approximation is unreliable there - a flag for outliers or model
#' misspecification.
#'
#' @param loo_obj A `loo` object (from [compute_loo()]).
#'
#' @return A `ggplot` object.
#' @export
#'
#' @examples
#' \dontrun{
#' plot_loo_pareto_k(compute_loo(fit))
#' }
plot_loo_pareto_k <- function(loo_obj) {
  k  <- loo_obj$diagnostics$pareto_k
  df <- tibble::tibble(obs = seq_along(k), k = k)
  pal <- bayes_palette()

  ggplot2::ggplot(df, ggplot2::aes(x = .data$obs, y = .data$k)) +
    ggplot2::geom_hline(
      yintercept = c(0.5, 0.7), linetype = "dashed",
      colour = c("grey50", pal[["null"]])
    ) +
    ggplot2::geom_point(ggplot2::aes(colour = .data$k > 0.7), alpha = 0.7) +
    ggplot2::scale_colour_manual(
      values = c(`FALSE` = pal[["control"]], `TRUE` = pal[["null"]]),
      guide  = "none"
    ) +
    ggplot2::labs(
      x = "observation", y = "Pareto k",
      title    = "PSIS-LOO Pareto-k diagnostic",
      subtitle = "Points above 0.7 (highlighted) make LOO unreliable there"
    ) +
    theme_bayes()
}


#' LOO-PIT calibration
#'
#' Leave-one-out predictive PIT values compared to uniform. A more honest
#' calibration check than in-sample PIT because it uses out-of-sample
#' predictions.
#'
#' @param model A fitted Stan model.
#' @param y Numeric vector of observed outcomes.
#' @param yrep_var Name of the replicated-outcome vector family.
#' @param loglik_var Name of the pointwise log-likelihood vector family.
#'
#' @return A `ggplot` object.
#' @export
#'
#' @examples
#' \dontrun{
#' plot_loo_pit(fit, y = df$y, yrep_var = "y_rep")
#' }
plot_loo_pit <- function(
    model,
    y,
    yrep_var   = "y_rep",
    loglik_var = "log_lik"
) {
  rlang::check_installed("loo", reason = "for LOO-PIT calibration.")
  yrep  <- draws_matrix_of(model, yrep_var)
  ll    <- draws_matrix_of(model, loglik_var)
  r_eff <- loo::relative_eff(exp(ll), chain_id = rep(1, nrow(ll)))
  psis  <- loo::psis(-ll, r_eff = r_eff)
  wts   <- weights(psis)

  bayesplot::ppc_loo_pit_overlay(y = y, yrep = yrep, lw = wts) +
    ggplot2::labs(
      title    = "LOO-PIT calibration",
      subtitle = "Empirical LOO-PIT (dark) vs. uniform reference"
    ) +
    theme_bayes()
}


#' Compare several models on ELPD
#'
#' Pass a NAMED list of `loo` objects. Plots `elpd_diff` ± 2 SE relative to
#' the best model; bars whose interval excludes zero are meaningfully worse.
#'
#' @param loo_list A named list of `loo` objects.
#'
#' @return A `ggplot` object.
#' @export
#'
#' @examples
#' \dontrun{
#' plot_loo_compare(list(M1 = loo_m1, M2 = loo_m2))
#' }
plot_loo_compare <- function(loo_list) {
  rlang::check_installed("loo", reason = "for model comparison.")
  stopifnot(is.list(loo_list), !is.null(names(loo_list)))

  cmp <- loo::loo_compare(loo_list)
  df  <- tibble::as_tibble(cmp, rownames = "model") |>
    dplyr::mutate(model = factor(.data$model, levels = rev(.data$model)))
  pal <- bayes_palette()

  ggplot2::ggplot(df, ggplot2::aes(x = .data$elpd_diff, y = .data$model)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey30") +
    ggplot2::geom_pointrange(
      ggplot2::aes(
        xmin = .data$elpd_diff - 2 * .data$se_diff,
        xmax = .data$elpd_diff + 2 * .data$se_diff
      ),
      colour = pal[["effect"]]
    ) +
    ggplot2::labs(
      x = "ELPD difference vs. best model (+/- 2 SE)", y = NULL,
      title    = "Model comparison (PSIS-LOO)",
      subtitle = "Higher = better predictive fit; best model fixed at 0"
    ) +
    theme_bayes()
}