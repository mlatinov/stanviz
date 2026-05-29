# =============================================================================
# diagnostics.R
# -----------------------------------------------------------------------------
# Sampler / convergence diagnostics. Before trusting ANY estimate, confirm the
# chains actually explored the posterior: R-hat near 1, healthy ESS, no
# divergences, sane energy.
# =============================================================================


#' Trace plots
#'
#' Classic trace plots. Well-mixed chains look like fuzzy caterpillars that
#' overlap. With no `pars` argument, picks the first few non-bookkeeping
#' scalars to avoid plotting hundreds of parameters by accident.
#'
#' @param model A fitted Stan model.
#' @param pars Optional character vector of parameter names.
#' @param n_auto If `pars` is `NULL`, how many parameters to pick
#'   automatically.
#'
#' @return A `ggplot` object.
#' @export
#'
#' @examples
#' \dontrun{
#' plot_trace(fit, pars = c("alpha", "beta", "sigma"))
#' }
plot_trace <- function(
    model,
    pars   = NULL,
    n_auto = 8
) {
  draws <- as_draws_safe(model, variables = pars)
  if (is.null(pars)) {
    vars  <- setdiff(posterior::variables(draws), c("lp__"))
    pars  <- utils::head(vars, n_auto)
    draws <- posterior::subset_draws(draws, variable = pars)
  }
  bayesplot::mcmc_trace(draws) +
    ggplot2::labs(
      title    = "Trace plots",
      subtitle = "Overlapping, stationary 'caterpillars' = good mixing"
    ) +
    theme_bayes()
}


#' Rank (trank) plots
#'
#' A more sensitive successor to trace plots. Under good mixing every chain's
#' ranks are uniform, so the bars sit flat across chains.
#'
#' @param model A fitted Stan model.
#' @param pars Optional character vector of parameter names.
#' @param n_auto If `pars` is `NULL`, how many parameters to pick.
#'
#' @return A `ggplot` object.
#' @export
#'
#' @examples
#' \dontrun{
#' plot_rank(fit, pars = c("alpha", "beta"))
#' }
plot_rank <- function(
    model,
    pars   = NULL,
    n_auto = 8
) {
  draws <- as_draws_safe(model, variables = pars)
  if (is.null(pars)) {
    vars  <- setdiff(posterior::variables(draws), c("lp__"))
    pars  <- utils::head(vars, n_auto)
    draws <- posterior::subset_draws(draws, variable = pars)
  }
  bayesplot::mcmc_rank_overlay(draws) +
    ggplot2::labs(
      title    = "Rank plots (trank)",
      subtitle = "Uniform, overlapping histograms across chains = converged"
    ) +
    theme_bayes()
}


#' R-hat across all parameters
#'
#' Dotplot of R-hat for every parameter. The reference at 1.01 is the modern
#' threshold; anything to the right warrants concern.
#'
#' @param model A fitted Stan model.
#'
#' @return A `ggplot` object.
#' @export
#'
#' @examples
#' \dontrun{
#' plot_rhat(fit)
#' }
plot_rhat <- function(model) {
  draws <- as_draws_safe(model)
  summ  <- posterior::summarise_draws(draws, rhat = posterior::rhat)
  bayesplot::mcmc_rhat(summ$rhat) +
    ggplot2::labs(
      title    = "R-hat across all parameters",
      subtitle = "Values should sit below 1.01"
    ) +
    bayesplot::yaxis_text(hjust = 1) +
    theme_bayes()
}


#' Effective sample size ratio
#'
#' ESS divided by total draws. Low ratios mean high autocorrelation: your
#' nominal draw count overstates the real information.
#'
#' @param model A fitted Stan model.
#' @param type Either `"bulk"` (default) or `"tail"`.
#'
#' @return A `ggplot` object.
#' @export
#'
#' @examples
#' \dontrun{
#' plot_ess(fit, type = "bulk")
#' }
plot_ess <- function(model, type = c("bulk", "tail")) {
  type    <- match.arg(type)
  draws   <- as_draws_safe(model)
  ess_fun <- if (type == "bulk") posterior::ess_bulk else posterior::ess_tail
  summ    <- posterior::summarise_draws(draws, ess = ess_fun)
  n_draws <- posterior::ndraws(draws)

  bayesplot::mcmc_neff(summ$ess / n_draws) +
    ggplot2::labs(
      title    = paste0(tools::toTitleCase(type), " effective sample size ratio"),
      subtitle = "Ratios below 0.1 flag heavy autocorrelation"
    ) +
    bayesplot::yaxis_text(hjust = 1) +
    theme_bayes()
}


#' NUTS energy diagnostic
#'
#' The two overlaid histograms - marginal energy and energy transitions -
#' should look alike. A much narrower transition distribution signals the
#' sampler struggling with the geometry (often fixed by reparameterisation).
#' Requires a CmdStanR or rstan fit (needs sampler diagnostics).
#'
#' @param model A `CmdStanMCMC` or `stanfit` object.
#'
#' @return A `ggplot` object.
#' @export
#'
#' @examples
#' \dontrun{
#' plot_energy(fit)
#' }
plot_energy <- function(model) {
  np <-
    if (inherits(model, "CmdStanMCMC")) {
      bayesplot::nuts_params(model)
    } else if (inherits(model, "stanfit")) {
      bayesplot::nuts_params(model)
    } else {
      stop("plot_energy() needs a CmdStanR or rstan fit with sampler diagnostics.")
    }
  bayesplot::mcmc_nuts_energy(np) +
    ggplot2::labs(
      title    = "NUTS energy diagnostic",
      subtitle = "Overlapping histograms = well-behaved geometry"
    ) +
    theme_bayes()
}


#' Pairs plot with divergent transitions highlighted
#'
#' Divergences clustering in one region of parameter space pinpoint where the
#' geometry is pathological - the classic funnel. The single most useful plot
#' for deciding what to reparameterise.
#'
#' @param model A `CmdStanMCMC` or `stanfit` object.
#' @param pars Character vector of parameter names. At least two required.
#'
#' @return A `bayesplot` pairs grid (a `gtable` / ggplot composite).
#' @export
#'
#' @examples
#' \dontrun{
#' plot_pairs_divergences(fit, pars = c("alpha", "tau"))
#' }
plot_pairs_divergences <- function(
    model,
    pars
) {
  if (missing(pars) || length(pars) < 2) {
    stop("Supply at least two parameter names in `pars`.")
  }
  np <-
    if (inherits(model, "CmdStanMCMC") || inherits(model, "stanfit")) {
      bayesplot::nuts_params(model)
    } else {
      stop("plot_pairs_divergences() needs a CmdStanR or rstan fit.")
    }
  draws <- as_draws_safe(model, variables = pars)
  bayesplot::mcmc_pairs(
    draws,
    np            = np,
    pars          = pars,
    off_diag_args = list(size = 0.8, alpha = 0.4)
  )
}


#' Diagnostic summary table
#'
#' Tidy tibble of mean, sd, quantiles, R-hat, and both ESS measures. Not a
#' plot, but the companion table you want next to the plots in a report.
#'
#' @param model A fitted Stan model.
#' @param pars Optional character vector of parameter names; `NULL` for all.
#'
#' @return A [tibble::tibble].
#' @export
#'
#' @examples
#' \dontrun{
#' diagnostic_summary(fit, pars = c("alpha", "beta"))
#' }
diagnostic_summary <- function(model, pars = NULL) {
  draws <- as_draws_safe(model, variables = pars)
  posterior::summarise_draws(
    draws,
    mean, sd,
    ~ posterior::quantile2(.x, probs = c(0.025, 0.5, 0.975)),
    rhat      = posterior::rhat,
    ess_bulk  = posterior::ess_bulk,
    ess_tail  = posterior::ess_tail
  ) |>
    tibble::as_tibble()
}