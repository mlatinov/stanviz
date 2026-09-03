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


#' Parse a Stan variable name/index spec into base name + index filter
#'
#' Not exported. Handles plain scalar names (`"mu"`), simple vector families
#' (`"mu[2]"`), and R-array-style slices of nested/multi-dimensional
#' containers, e.g. `"samples_combined[1,]"` (fix the first index, keep every
#' element of the rest) or `"samples_combined[,2,3]"`. A blank slot means
#' "every value of this dimension"; an omitted trailing slot is treated the
#' same as a blank one.
#'
#' @param var A single Stan variable name, optionally with a trailing
#'   `[i,j,...]` index/slice spec.
#'
#' @return A list with `base` (character) and `idx` (`NULL`, or a list with
#'   one element per dimension: `NA_integer_` for "all", or an integer for a
#'   fixed index).
#' @keywords internal
#' @noRd
parse_stan_var <- function(var) {
  m <- regmatches(var, regexec("^([A-Za-z_][A-Za-z0-9_]*)(?:\\[(.*)\\])?$", var))[[1]]
  if (!length(m)) {
    stop(sprintf("'%s' is not a valid Stan variable name.", var), call. = FALSE)
  }
  base     <- m[2]
  spec_str <- m[3]

  idx <- NULL
  if (nzchar(spec_str)) {
    n_commas <- lengths(regmatches(spec_str, gregexpr(",", spec_str, fixed = TRUE)))
    n_tok    <- n_commas + 1L
    tok      <- strsplit(spec_str, ",", fixed = TRUE)[[1]]
    if (length(tok) < n_tok) tok <- c(tok, rep("", n_tok - length(tok)))
    tok <- trimws(tok)
    idx <- lapply(tok, function(t) {
      if (!nzchar(t)) return(NA_integer_)
      n <- suppressWarnings(as.integer(t))
      if (is.na(n)) {
        stop(sprintf(
          "Invalid index '%s' in '%s'; use integers, or leave a slot blank for 'all'.",
          t, var
        ), call. = FALSE)
      }
      n
    })
  }
  list(base = base, idx = idx)
}


#' Parse the trailing `[i,j,...]` of Stan draws column names into a matrix
#'
#' Not exported. `names_vec` must already be a family sharing the same
#' number of dimensions (e.g. everything matched by `^base\\[`).
#'
#' @param names_vec Character vector of draws variable names such as
#'   `"samples_combined[1,2,3]"`.
#'
#' @return An integer matrix, one row per name, one column per dimension.
#' @keywords internal
#' @noRd
stan_index_matrix <- function(names_vec) {
  bracket <- regmatches(names_vec, regexpr("\\[[^][]*\\]$", names_vec))
  inner   <- substring(bracket, 2, nchar(bracket) - 1)
  split   <- strsplit(inner, ",", fixed = TRUE)
  n_dim   <- max(lengths(split))
  mat     <- do.call(rbind, lapply(split, as.integer))
  mat[, seq_len(n_dim), drop = FALSE]
}


#' Resolve a (possibly indexed/sliced) Stan variable spec against draws
#'
#' Not exported. The workhorse behind [gather_indexed()] and
#' [draws_matrix_of()]: turns a name like `"mu"`, `"mu[2]"`, or
#' `"samples_combined[1,]"` into the exact, correctly-ordered set of draws
#' column names it refers to, so nested containers (`array[J] vector[N]`,
#' `array[K] matrix[R, C]`, ...) can be sliced the same way you'd slice an R
#' array.
#'
#' @param all_vars Character vector, typically `posterior::variables(draws)`.
#' @param var A Stan variable name, optionally with a trailing index/slice.
#'
#' @return A list with `hit_vars` (character, in row-major dimension order),
#'   `index_mat` (integer matrix of their parsed indices, or `NULL` for a
#'   true scalar), `free_dims` (integer positions, into the columns of
#'   `index_mat`, of the dimensions left unfixed - i.e. the ones that still
#'   vary across `hit_vars`), and `is_scalar` (logical).
#' @keywords internal
#' @noRd
resolve_stan_var <- function(all_vars, var) {
  spec <- parse_stan_var(var)
  base <- spec$base
  idx  <- spec$idx

  family <- all_vars[grepl(paste0("^", base, "\\["), all_vars)]

  if (!length(family)) {
    if (is.null(idx) && base %in% all_vars) {
      return(list(hit_vars = base, index_mat = NULL, free_dims = integer(0), is_scalar = TRUE))
    }
    stop(sprintf("No variable matching '%s' found in the model draws.", var), call. = FALSE)
  }

  mat       <- stan_index_matrix(family)
  keep      <- seq_len(nrow(mat))
  free_dims <- seq_len(ncol(mat))

  if (!is.null(idx)) {
    if (length(idx) > ncol(mat)) {
      stop(sprintf(
        "'%s' has %d dimension(s), but '%s' supplies %d index position(s).",
        base, ncol(mat), var, length(idx)
      ), call. = FALSE)
    }
    if (length(idx) < ncol(mat)) idx <- c(idx, rep(list(NA_integer_), ncol(mat) - length(idx)))
    free_dims <- which(vapply(idx, is.na, logical(1)))
    for (d in seq_along(idx)) {
      if (!is.na(idx[[d]])) keep <- keep[mat[keep, d] == idx[[d]]]
    }
    if (!length(keep)) {
      stop(sprintf("No elements of '%s' match the index '%s'.", base, var), call. = FALSE)
    }
  }

  ord  <- do.call(order, as.data.frame(mat[keep, , drop = FALSE]))
  keep <- keep[ord]

  list(
    hit_vars  = family[keep],
    index_mat = mat[keep, , drop = FALSE],
    free_dims = free_dims,
    is_scalar = FALSE
  )
}


#' Tidy an indexed parameter family into long format
#'
#' Stan vector parameters arrive as `name[1]`, `name[2]`, .... This pulls a
#' single family out and returns a long tibble with one row per draw per
#' index. Nested/multi-dimensional families (`array[J] vector[N]`,
#' `array[K] matrix[R, C]`, ...) work too: pass the base name for every
#' element, or an R-array-style slice such as `"samples_combined[1,]"` to
#' fix one dimension and keep the rest.
#'
#' @param model A fitted model (anything [as_draws_safe()] accepts).
#' @param var Base name of the parameter family, e.g. `"mu"` or `"y_rep"`,
#'   optionally with a trailing `[i,j,...]` index/slice (blank slots mean
#'   "keep every value of this dimension").
#'
#' @return A [tibble::tibble] with columns `.draw`, `param`, `value`, and
#'   either `index` (whenever exactly one dimension actually varies - a
#'   plain vector, or a slice of a nested container that fixes every
#'   dimension but one) or `index1`, `index2`, ... (one per dimension still
#'   varying, for an unsliced multi-dimensional family).
#' @export
#'
#' @examples
#' \dontrun{
#' gather_indexed(fit, "y_rep")
#'
#' # array[3] vector[N] samples_combined -> first array slot, every N
#' gather_indexed(fit, "samples_combined[1,]")
#' }
gather_indexed <- function(model, var) {
  draws    <- as_draws_safe(model)
  all_vars <- posterior::variables(draws)
  res      <- resolve_stan_var(all_vars, var)

  long <- posterior::subset_draws(draws, variable = res$hit_vars) |>
    posterior::as_draws_df() |>
    tibble::as_tibble() |>
    dplyr::select(".draw", dplyr::all_of(res$hit_vars)) |>
    tidyr::pivot_longer(
      -".draw",
      names_to  = "param",
      values_to = "value"
    )

  if (res$is_scalar) {
    long$index <- 1L
    return(long)
  }

  row_idx <- match(long$param, res$hit_vars)
  free    <- res$free_dims
  if (length(free) <= 1L) {
    # A single (or fully-fixed) free dimension collapses to the classic,
    # backward-compatible `index` column - including the plain-vector case.
    long$index <- if (length(free) == 1L) res$index_mat[row_idx, free] else 1L
  } else {
    idx_cols <- as.data.frame(res$index_mat[row_idx, free, drop = FALSE])
    names(idx_cols) <- paste0("index", seq_along(free))
    long <- dplyr::bind_cols(long, idx_cols)
  }
  long
}


#' Build a draws-by-observation matrix for a vector generated quantity
#'
#' Many bayesplot `ppc_*` functions want a matrix of dimension
#' `[n_draws, N]`. This builds that matrix from a Stan variable family such
#' as `"y_rep"`, with columns ordered by index. Nested/multi-dimensional
#' families (`array[J] vector[N]`, `array[K] matrix[R, C]`, ...) work too:
#' pass an R-array-style slice such as `"samples_combined[1,]"` to fix one
#' dimension and keep every value of the rest as columns.
#'
#' @param model A fitted model (anything [as_draws_safe()] accepts).
#' @param var Base name of a vector generated quantity, optionally with a
#'   trailing `[i,j,...]` index/slice (blank slots mean "keep every value of
#'   this dimension").
#'
#' @return A numeric matrix with `n_draws` rows and one column per matched
#'   element.
#' @export
#'
#' @examples
#' \dontrun{
#' yrep <- draws_matrix_of(fit, "y_rep")
#' dim(yrep)
#'
#' # array[3] vector[N] samples_combined -> first array slot, every N
#' draws_matrix_of(fit, "samples_combined[1,]")
#' }
draws_matrix_of <- function(model, var) {
  draws    <- as_draws_safe(model)
  all_vars <- posterior::variables(draws)
  res      <- resolve_stan_var(all_vars, var)
  if (res$is_scalar) {
    stop(sprintf("No vector variable matching '%s' found.", var), call. = FALSE)
  }

  m <- posterior::subset_draws(draws, variable = res$hit_vars) |>
    posterior::as_draws_matrix()
  m[, res$hit_vars, drop = FALSE]
}


#' Coerce a grouping vector to a labelled factor
#'
#' Not exported. Shared by every "grouped" plotting function in stanviz
#' (posterior predictive checks and hierarchical plots alike) so grouping
#' and labelling behave identically everywhere.
#'
#' @param group A vector of group memberships, one element per observation.
#' @param group_labels Optional character vector of display labels, one per
#'   distinct level of `group` (in `sort()` order). `NULL` keeps whatever
#'   levels `factor(group)` produces.
#'
#' @return A factor the same length as `group`.
#' @keywords internal
#' @noRd
as_group_factor <- function(group, group_labels = NULL) {
  group <- as.factor(group)
  if (!is.null(group_labels)) {
    if (length(group_labels) != nlevels(group)) {
      stop(sprintf(
        "`group_labels` has length %d but `group` has %d levels.",
        length(group_labels), nlevels(group)
      ), call. = FALSE)
    }
    levels(group) <- group_labels
  }
  group
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