# =============================================================================
# helpers.R
# -----------------------------------------------------------------------------
# Shared theme, palette, and low-level helpers used across stanviz.
# =============================================================================


#' Bayesian colour palette
#'
#' A named vector of four colours used throughout stanviz so every plot speaks
#' the same visual language.
#'
#' @param control Colour for the reference / untreated arm.
#' @param treatment Colour for the treated arm.
#' @param effect Colour for anything representing an estimated effect or
#'   posterior distribution.
#' @param null Colour for "no effect" references (zero lines, ROPE bands).
#'
#' @return A named character vector of length four with names
#'   `control`, `treatment`, `effect`, `null`.
#' @export
#'
#' @examples
#' bayes_palette()
#' bayes_palette(effect = "steelblue")
bayes_palette <- function(
    control   = "lightblue4",
    treatment = "red4",
    effect    = "red4",
    null      = "#E76F51"
) {
  c(control = control, treatment = treatment, effect = effect, null = null)
}


#' A clean, publication-friendly ggplot theme
#'
#' Light horizontal grid only, bold titles, legend at the bottom.
#' All sizes scale with `base_size`.
#'
#' @param base_size Base font size in points.
#'
#' @return A ggplot2 theme object.
#' @export
#'
#' @examples
#' library(ggplot2)
#' ggplot(mtcars, aes(wt, mpg)) + geom_point() + theme_bayes()
theme_bayes <- function(base_size = 12) {
  ggplot2::theme_minimal(base_size = base_size, base_family = "sans") +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(face = "bold", size = ggplot2::rel(1.15)),
      plot.subtitle    = ggplot2::element_text(colour = "grey30", size = ggplot2::rel(0.95)),
      plot.caption     = ggplot2::element_text(colour = "grey40", size = ggplot2::rel(0.8), hjust = 0),
      axis.title       = ggplot2::element_text(face = "plain"),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(colour = "grey92"),
      panel.grid.major.x = ggplot2::element_blank(),
      strip.text       = ggplot2::element_text(face = "bold"),
      legend.position  = "bottom"
    )
}


#' Apply the stanviz theme to bayesplot
#'
#' Convenience wrapper that sets bayesplot's colour scheme and theme so that
#' [bayesplot::ppc_dens_overlay()] and friends match the rest of your figures.
#' Call once at the top of a script.
#'
#' @param scheme A bayesplot colour scheme name; see
#'   [bayesplot::color_scheme_set()].
#'
#' @return Invisibly `TRUE`.
#' @export
#'
#' @examples
#' \dontrun{
#' set_bayesplot_theme("red")
#' }
set_bayesplot_theme <- function(scheme = "red") {
  bayesplot::color_scheme_set(scheme)
  bayesplot::bayesplot_theme_set(theme_bayes())
  invisible(TRUE)
}


#' Robustly extract posterior draws from a fitted model
#'
#' Accepts a CmdStanR fit, an rstan `stanfit`, an already-extracted draws
#' object, or anything `posterior::as_draws_df()` understands. Returns a
#' `draws_df` for downstream consumption.
#'
#' @param model A fitted model object.
#' @param variables Optional character vector of parameter names to keep.
#'   `NULL` (the default) keeps all.
#'
#' @return A [posterior::draws_df] object.
#' @export
#'
#' @examples
#' \dontrun{
#' draws <- as_draws_safe(fit, variables = c("alpha", "beta"))
#' }
as_draws_safe <- function(model, variables = NULL) {
  draws <-
    if (inherits(model, "draws")) {
      posterior::as_draws_df(model)
    } else if (inherits(model, "CmdStanMCMC") || inherits(model, "CmdStanFit")) {
      model$draws(variables = variables, format = "draws_df")
    } else if (inherits(model, "stanfit")) {
      posterior::as_draws_df(model)
    } else {
      posterior::as_draws_df(model)
    }

  if (!is.null(variables) && !inherits(model, "CmdStanMCMC")) {
    keep <- intersect(variables, posterior::variables(draws))
    if (length(keep)) draws <- posterior::subset_draws(draws, variable = keep)
  }
  draws
}


#' Tidy an indexed parameter family into long format
#'
#' Stan vector parameters arrive as `name[1]`, `name[2]`, .... This pulls a
#' single family out and returns a long tibble with one row per draw per index.
#'
#' @param model A fitted model (anything [as_draws_safe()] accepts).
#' @param var Base name of the parameter family, e.g. `"mu"` or `"y_rep"`.
#'
#' @return A [tibble::tibble] with columns `.draw`, `param`, `value`, `index`.
#' @export
#'
#' @examples
#' \dontrun{
#' gather_indexed(fit, "y_rep")
#' }
gather_indexed <- function(model, var) {
  draws <- as_draws_safe(model)
  all_vars <- posterior::variables(draws)
  pattern  <- paste0("^", var, "\\[")
  hit_vars <- all_vars[grepl(pattern, all_vars)]
  if (!length(hit_vars)) {
    hit_vars <- all_vars[all_vars == var]
  }
  if (!length(hit_vars)) {
    stop(sprintf("No variable matching '%s' found in the model draws.", var))
  }

  posterior::subset_draws(draws, variable = hit_vars) |>
    posterior::as_draws_df() |>
    tibble::as_tibble() |>
    dplyr::select(".draw", dplyr::all_of(hit_vars)) |>
    tidyr::pivot_longer(
      -".draw",
      names_to  = "param",
      values_to = "value"
    ) |>
    dplyr::mutate(
      index = readr::parse_number(sub(".*\\[", "[", .data$param)),
      index = ifelse(is.na(.data$index), 1L, as.integer(.data$index))
    )
}


#' Build a draws-by-observation matrix for a vector generated quantity
#'
#' Many bayesplot `ppc_*` functions want a matrix of dimension
#' `[n_draws, N]`. This builds that matrix from a Stan variable family such
#' as `"y_rep"`, with columns ordered by observation index.
#'
#' @param model A fitted model (anything [as_draws_safe()] accepts).
#' @param var Base name of a vector generated quantity.
#'
#' @return A numeric matrix with `n_draws` rows and `N` columns.
#' @export
#'
#' @examples
#' \dontrun{
#' yrep <- draws_matrix_of(fit, "y_rep")
#' dim(yrep)
#' }
draws_matrix_of <- function(model, var) {
  draws <- as_draws_safe(model)
  all_vars <- posterior::variables(draws)
  hit_vars <- all_vars[grepl(paste0("^", var, "\\["), all_vars)]
  if (!length(hit_vars)) {
    stop(sprintf("No vector variable matching '%s' found.", var))
  }
  ord <- order(readr::parse_number(sub(".*\\[", "[", hit_vars)))
  hit_vars <- hit_vars[ord]

  m <- posterior::subset_draws(draws, variable = hit_vars) |>
    posterior::as_draws_matrix()
  m[, hit_vars, drop = FALSE]
}


#' Dashed "no effect" reference line
#'
#' Small wrapper around [ggplot2::geom_vline()] used by several stanviz
#' plots to mark a null reference.
#'
#' @param xintercept Numeric position of the line.
#' @param colour Line colour.
#'
#' @return A ggplot2 layer.
#' @export
#'
#' @examples
#' library(ggplot2)
#' ggplot(mtcars, aes(wt, mpg)) + geom_point() + geom_null_line(3)
geom_null_line <- function(xintercept = 0, colour = "grey30") {
  ggplot2::geom_vline(
    xintercept = xintercept,
    linetype   = "dashed",
    colour     = colour
  )
}