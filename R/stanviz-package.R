#' stanviz: Plotting Helpers for Stan Models
#'
#' Reusable ggplot2-based plotting functions for fitted Stan models —
#' prior/posterior predictive checks, MCMC diagnostics, parameter posteriors,
#' conditional effects, mediation decomposition, prior sensitivity, and
#' PSIS-LOO model comparison.
#'
#' @keywords internal
#' @importFrom rlang .data %||%
"_PACKAGE"

utils::globalVariables(c(".draw", ".chain", ".iteration"))