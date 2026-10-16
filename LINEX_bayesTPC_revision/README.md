# LINEX revision: bayesTPC analysis

This folder replaces the earlier revision bundle. Models and samplers are built with bayesTPC, and its compiled NIMBLE sampler runs the chains. The two CSVs are the same inputs as before. The analysis and the RH90 figure script are both included here.

## Run the main analysis

Install `bayesTPC`, `nimble`, `coda`, `posterior` and `loo`, with the C++ toolchain needed by NIMBLE. From the directory containing this folder, run:

```r
source("LINEX_bayesTPC_revision/run_LINEX_revision.R")
```

Or in a terminal:

```sh
Rscript LINEX_bayesTPC_revision/run_LINEX_revision.R
```

This runs the three stephensi RH datasets, the Aedes example, and the four simulated examples. There are six chains per fit. It reports diagnostics, fitted curves and predictive checks. Prior-sensitivity fits and exact leave-one-temperature-out refits are now **off by default**: these were adding a large amount of work to what should have been a straightforward revision run. They remain available by setting `run_prior_sensitivity = TRUE` and `run_exact_loo_if_needed = TRUE` in the `linex_revision` options. If PSIS-LOO is unreliable with exact refits switched off, the script leaves the ranking unwritten.

For a first check of RH90 only, use this in RStudio:

```r
options(linex_revision = list(
  dataset_filter = "stephensi_RH90",
  run_simulations = FALSE,
  run_prior_sensitivity = FALSE,
  run_exact_loo_if_needed = FALSE
))
source("LINEX_bayesTPC_revision/run_LINEX_revision.R")
```

Reset `options(linex_revision = NULL)` before running the full analysis.

The main analysis fits the unrestricted LINEX with starts on both sides of zero, then explicit `a > 0` and `a < 0` versions, plus the quadratic for the real datasets. The sign-specific versions use the original normal prior truncated to the stated sign. We can compare them as conditional models; their results do not establish the relative posterior weights of the two modes in the unrestricted model.

The model is now built **once per fit**, and its compiled NIMBLE sampler is reused for the six chains. Each chain gets separate starting values, a separate seed and its full burn-in; adaptation resets for each chain. The two setup draws used when building the bayesTPC object are discarded. One separate R process handles each fit, then exits to release its compilation libraries. Each completed chain is still cached individually, so an interrupted fit can resume with only the missing chains. The fit's `_worker.log` file records compilation and sampling progress.

## Resuming after the DLL-limit error

Keep the existing `results/` directory. Replace `run_LINEX_revision.R` in the same analysis folder with this version, then restart R once to clear the libraries left by the failed run. To restore the settings used in that run and continue:

```r
cfg <- readRDS("LINEX_bayesTPC_revision/results/bayesTPC/settings.rds")
cfg$run_prior_sensitivity <- FALSE
cfg$run_exact_loo_if_needed <- FALSE
options(linex_revision = cfg)
source("LINEX_bayesTPC_revision/run_LINEX_revision.R")
```

Use the actual analysis-folder path if yours has a different name. Completed matching chain caches are loaded; a chain interrupted before it was saved must run again. This update retains the existing chain-cache format, model definitions, priors and chain lengths. The two expensive additional analyses are disabled by the commands above; their already saved fits are retained. The DLL fix does not require raising R's library limit.

The WAIC variance and Pareto-k warnings are separate from the crash. They are still relevant to predictive comparison. Exact refits are no longer triggered automatically under the new default settings.

## What to look at

The main script creates `results/bayesTPC/`. This keeps the new results apart from the previous analysis. No fits or final result figures are supplied in this ZIP.

- `diagnostics.csv` and `parameter_diagnostics.csv`: Rhat and effective sample sizes, including the log posterior. Split chains should fail these checks.
- `figures/`: a 2-by-2 trace PDF for every fit, with all six chains overlaid. Curves are also plotted for fits that pass the checks.
- `posterior_predictive_checks.csv` and `predictive_checks_by_temperature.csv`: how well the fitted models reproduce the observations.
- `predictive_scores.csv`: WAIC and PSIS-LOO for fits that pass diagnostics, with warnings recorded. WAIC alone is not a convergence check.
- `model_comparisons.csv`: LOO comparisons, only when every candidate has passed the required checks. With exact refits switched off, an unreliable PSIS result leaves this comparison unwritten.
- `prior_sensitivity.csv` and `simulation_recovery.csv`: produced when the corresponding fits pass diagnostics.

These are generated outputs, so some will be absent if a run is unfinished or the relevant fits fail. The chain caches are in `results/bayesTPC/fits/`.

## RH90 figure and unrestricted starting-point check

The second script is included:

```r
source("LINEX_bayesTPC_revision/plot_RH90_bayesTPC_six_chains.R")
```

This uses the original unrestricted LINEX model for both sets of chains: six original/default-start chains and six negative-start chains. Unlike the sign-specific main-analysis fits, **the negative-start chains are free to cross zero**. It also fits the quadratic with six chains.

After a successful run, `results/figures/` contains:

- `RH90_bayesTPC_LINEX_positive_six_chain_trace.pdf`
- `RH90_bayesTPC_LINEX_negative_six_chain_trace.pdf`
- `RH90_bayesTPC_quadratic_six_chain_trace.pdf`
- `RH90_bayesTPC_four_panel_diagnostic.pdf`

Each trace figure is a single 2-by-2 page with six coloured chains. The curve figure shows the two LINEX modes, a 50:50 combination and the quadratic. That equal weighting is a sensitivity illustration, not an estimate of mode probabilities or stacking weights. The CSV `results/RH90_bayesTPC_chain_modes.csv` records which regions the chains sampled.

If we already have RH90 fit caches, copy the three six-chain RDS files from the old `results/` directory into this folder's `results/` directory. The figure script reuses them. The old single-chain negative-start RDS cannot replace the six-chain negative-start fit.

We should inspect the traces and numerical diagnostics before deciding what replaces Figure 5. Finding two posterior modes is not the same as showing that the biological temperature response has two peaks.
