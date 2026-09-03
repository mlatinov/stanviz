
<!-- README.md is generated from README.Rmd. Please edit that file. -->

# stanviz

> Reusable ggplot2-based plotting functions for fitted Stan models.

`stanviz` standardises the figures you produce when working with
Bayesian models in Stan: prior and posterior predictive checks, MCMC
diagnostics, coefficient inference, conditional and counterfactual
effects, mediation decomposition, hierarchical / multilevel results,
prior sensitivity, and PSIS-LOO model comparison. Every function takes a
fitted Stan model (CmdStanR first-class; rstan and `draws_df` also
accepted) and returns a `ggplot` object that’s ready to ship into a
report.

The design assumes the typical Stan workflow: declare the quantities you
care about (`y_rep`, `log_lik`, `pit`, per-arm predictions, mediation
scalars, group-level effects) in `generated quantities`, then point the
plotting functions at them by name. The plots layer comes ready; the
modelling layer stays in your hands.

## Installation

``` r
# install.packages("remotes")
remotes::install_github("mlatinov/stanviz")
```

The heavier dependencies (`priorsense`, `loo`, `cmdstanr`, `rstan`) are
declared as **Suggests**, so they aren’t installed automatically.
Install them as needed:

``` r
install.packages(c("loo", "priorsense"))
# cmdstanr separately:
install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev",
                                       getOption("repos")))
```

## Quick start

``` r
library(stanviz)

# Convergence first
plot_trace(fit, pars = c("alpha", "beta", "sigma"))
plot_rhat(fit)
plot_energy(fit)

# Does the model regenerate the data?
plot_ppc_dens(fit, y = df$y, yrep_var = "y_rep")
plot_ppc_stat_grid(fit, y = df$y, yrep_var = "y_rep")
plot_ppc_pit(fit, pit_var = "pit")

# Coefficient inference
plot_coef_forest(fit, pars = c("beta_1", "beta_2", "beta_3"))
plot_rope(fit, par = "treatment_effect", rope = c(-0.1, 0.1))

# Out-of-sample fit
loo_obj <- compute_loo(fit)
plot_loo_pareto_k(loo_obj)
```

## Function map

The library is organised by what question the plot answers. Use this as
a lookup table when you’re not sure which function to reach for.

### Prior predictive — *“are my priors sane before fitting?”*

| Function | Purpose |
|----|----|
| `plot_prior_predictive_dens()` | Outcome densities implied by the priors alone. |
| `plot_prior_predictive_stat()` | Distribution of a summary statistic (mean, sd, max). |
| `plot_prior_draws()` | Realised marginal priors per parameter. |

### Posterior predictive — *“does the fitted model regenerate the data?”*

| Function | Purpose |
|----|----|
| `plot_ppc_dens()` | Observed density vs. many replicates. |
| `plot_ppc_ecdf()` | Observed vs. replicated empirical CDFs - tails/skew often show up here first. |
| `plot_ppc_ecdf_grouped()` | `plot_ppc_ecdf()`, one panel per group. |
| `plot_ppc_boxplot()` | Observed vs. a handful of replicated datasets, as boxplots. |
| `plot_ppc_boxplot_grouped()` | `plot_ppc_boxplot()`, one panel per group. |
| `plot_ppc_violin_grouped()` | Replicated distribution (violin) vs. observed, per group. |
| `plot_ppc_violin()` | Single-panel `plot_ppc_violin_grouped()`. |
| `plot_ppc_stat()` | Single test statistic, observed vs. replicates. |
| `plot_ppc_stat_grid()` | Several test statistics at once (centre, spread, tails). |
| `plot_ppc_intervals()` | Per-observation predictive intervals. |
| `plot_ppc_error_scatter()` | Predictive error structure across fitted values. |
| `plot_ppc_pit()` | PIT calibration histogram. |

### MCMC diagnostics — *“did the sampler actually work?”*

| Function | Purpose |
|----|----|
| `plot_trace()` | Trace plots (overlapping caterpillars = good). |
| `plot_rank()` | Rank/trank plots — more sensitive than traces. |
| `plot_rhat()` | R-hat across every parameter. |
| `plot_ess()` | Effective sample size ratio. |
| `plot_energy()` | NUTS energy diagnostic. |
| `plot_pairs_divergences()` | Pairs plot with divergent transitions. |
| `diagnostic_summary()` | Tidy table of mean / sd / quantiles / Rhat / ESS. |

### Parameter posteriors — *“what did we learn?”*

| Function | Purpose |
|----|----|
| `plot_coef_forest()` | Forest plot of several coefficients. |
| `plot_param_halfeye()` | Halfeye densities (shape, not just intervals). |
| `plot_rope()` | Effect posterior vs. Region Of Practical Equivalence. |
| `plot_threshold_curve()` | P(effect beyond threshold) as a continuous curve. |

### Conditional & causal effects — *“on the outcome scale, what’s the effect?”*

| Function | Purpose |
|----|----|
| `extract_ate_from_draws()` | Build ATE draws from per-arm Stan predictions. |
| `plot_ate_posterior()` | Halfeye of the ATE. |
| `plot_counterfactual_means()` | Marginal means under each arm. |
| `plot_conditional_effect_draws()` | Partial-effect curve from a tidy `(x, .draw, value[, group])` tibble. |
| `plot_effect_curve()` | Posterior response curve straight off a Stan predictor grid (e.g. an SOD -> H2O2 dose-response). |
| `plot_effect_curve_grouped()` | `plot_effect_curve()`, one curve per group, overlaid for direct comparison. |

### Mediation — *“how does the effect transmit?”*

| Function | Purpose |
|----|----|
| `plot_effect_decomposition()` | Direct, total-indirect, and total effects. |
| `plot_indirect_effects()` | Per-mediator natural indirect effects. |
| `plot_proportion_mediated()` | Posterior proportion mediated, with stability caveats. |
| `plot_mediation_paths()` | Path diagram with posterior-median coefficients. |

### Hierarchical / multilevel — *“how do the groups differ, and what does pooling do?”*

| Function | Purpose |
|----|----|
| `extract_group_effects()` | Tidy tibble of per-group posterior draws. |
| `plot_random_effects()` | Caterpillar plot of group-level effects (absolute or deviation). |
| `plot_shrinkage()` | Raw vs. partial-pooled group means, with shrinkage arrows. |
| `plot_variance_decomposition()` | Share of variance at each hierarchical level. |
| `plot_icc()` | Intraclass correlation posterior. |
| `plot_ppc_by_group()` | Faceted PPC: does the model fit equally across groups? |
| `plot_funnel_diagnostic()` | Group-level estimate vs. log(tau), with divergences. |
| `plot_random_correlation()` | Joint random intercepts/slopes + correlation posterior. |

### Prior sensitivity — *“are my conclusions driven by data or priors?”*

| Function | Purpose |
|----|----|
| `powerscale_summary()` | Tidy table of prior/likelihood sensitivity. |
| `plot_powerscale()` | Posterior movement under prior/likelihood scaling. |
| `plot_powerscale_quantities()` | Compact alternative summarising mean/sd response. |

### Model comparison — *“which model predicts out-of-sample better?”*

| Function | Purpose |
|----|----|
| `compute_loo()` | PSIS-LOO from a fitted Stan model’s `log_lik`. |
| `plot_loo_pareto_k()` | Pareto-k diagnostic (flags influential observations). |
| `plot_loo_pit()` | LOO-PIT calibration. |
| `plot_loo_compare()` | ELPD differences between models, with SE. |

### Helpers — *theming and infrastructure*

| Function | Purpose |
|----|----|
| `theme_bayes()` | The shared ggplot theme. |
| `bayes_palette()` | The shared four-colour palette. |
| `set_bayesplot_theme()` | Make bayesplot match. |
| `as_draws_safe()` | Robust draws extractor (CmdStanR / rstan / draws). |
| `draws_matrix_of()` | Build an `[n_draws, N]` matrix from a vector parameter. |
| `gather_indexed()` | Tidy an indexed parameter family to long format. |
| `geom_null_line()` | Dashed reference line at a given x. |

## What `stanviz` expects in your Stan code

To get the most out of the library, have your `generated quantities`
emit these standard names. Functions accept the variable name as an
argument, so you can deviate — but the defaults assume these:

| Quantity | Used by |
|----|----|
| `vector[N] y_rep` | Every PPC plot. |
| `vector[N] log_lik` | `compute_loo()` and all LOO plots. |
| `vector[N] pit` | `plot_ppc_pit()`. |
| `real lprior` (alongside `log_lik`) | `priorsense` power-scaling functions. |
| Per-arm vectors, e.g. `mu_treated`, `mu_control` | `extract_ate_from_draws()` and counterfactual plots. |
| A predictor-grid vector, e.g. `vector[G] mu_grid` (or `array[K] vector[G] mu_grid` for `K` groups) | `plot_effect_curve()` / `plot_effect_curve_grouped()`. |
| Mediation scalars: `direct_effect`, `total_effect`, `total_indirect_effect`, `NIE_*`, `proportion_mediated_*` | Mediation plots. |
| Hierarchical: a vector of group effects (`alpha_j`), the grand mean (`mu`), the group SD (`tau`) | Hierarchical plots. |

## Design conventions

- **Every plot returns a `ggplot`.** You can keep customising
  (`+ ggplot2::labs(...)`) or composing with patchwork (`p1 / p2`).
- **One palette, one theme.** All figures use `bayes_palette()` and
  `theme_bayes()` so multi-panel reports look coherent.
- **CmdStanR first-class.** All non-marginaleffects functions work with
  raw CmdStanR fits; nothing needs a brms wrapper.
- **You provide the variable names.** No magic auto-detection; the
  function takes what to plot as an explicit argument. Trades a few
  keystrokes for total transparency about what’s being shown.

## Development

``` r
# After editing R/ source:
devtools::document()      # regenerate help and NAMESPACE
devtools::load_all()      # reload without installing
devtools::check()         # R CMD check
devtools::install()       # install into your library
```

## License

MIT © Metodi Latinov
