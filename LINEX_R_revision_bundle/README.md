# LINEX revision

The reviewer found a real problem with the RH90 fit: the LINEX model has two very different solutions, and an individual chain can sit in one without ever visiting the other. In our first bayesTPC run, all six default-start chains ended up near positive `a` and `Topt = 30.3 C`. A separate chain started near negative `a` stayed around `Topt = 21.8 C`. We now need to see how **six** negative-start chains behave, alongside the six original chains and the quadratic fit. This is still the same unrestricted LINEX model; starting near negative `a` does not force it to remain there.

## Folder contents:

- `run_LINEX_revision.R` is the broader RStan check: six-chain fits across the example datasets, diagnostics, predictions, prior sensitivity and model comparisons.
- `plot_RH90_bayesTPC_six_chains.R` is the focused bayesTPC check for RH90 and the proposed Figure 5 layout.
- `data/` contains the two input CSVs.

Each script creates a `results/` directory beside itself when run. If you already have cached fits from an earlier run, keep that directory in place; the plotting script will reuse the original LINEX and quadratic fits. It must fit the *new* six negative-start chains once, saving them as `results/RH90_bayesTPC_LINEX_negative_start_6chains.rds`. The previous one-chain negative-start fit has a different filename and is not reused.

## Run the RH90 figure check

From the repository root, with `bayesTPC`, `nimble` and `coda` installed:

```sh
Rscript LINEX_R_revision_bundle/plot_RH90_bayesTPC_six_chains.R
```

Or in RStudio, from the repository root:

```r
source("LINEX_R_revision_bundle/plot_RH90_bayesTPC_six_chains.R")
```

Once it finishes, look in `LINEX_R_revision_bundle/results/figures/` for:

| PDF created by the script | What it shows |
| --- | --- |
| `RH90_bayesTPC_LINEX_positive_six_chain_trace.pdf` | Four overlaid trace plots for the original six LINEX chains: `a`, `Topt`, `cc`, and survival at 30 C. |
| `RH90_bayesTPC_LINEX_negative_six_chain_trace.pdf` | The same four plots for the six new negative-start LINEX chains. The LINEX plots share vertical scales. |
| `RH90_bayesTPC_quadratic_six_chain_trace.pdf` | Six quadratic chains overlaid for `B0`, `B1`, `B2`, and survival at 30 C. |
| `RH90_bayesTPC_four_panel_diagnostic.pdf` | RH90 observations with the positive-mode LINEX, negative-mode LINEX, an illustrative combined curve, and the quadratic curve. |

The script also creates `results/RH90_bayesTPC_chain_modes.csv` so we can check where each LINEX chain ended up. The combined curve gives equal weight to draws from each mode. **That 50:50 split is a visual sensitivity check, not an estimated posterior weight or stacking analysis.** We should not present it as a resolved full-posterior fit.

## Broader check

Run `Rscript LINEX_R_revision_bundle/run_LINEX_revision.R` from the repository root (or `source()` it in RStudio). It needs `rstan`, `loo`, `posterior` and a working C++ toolchain. The default run includes the other datasets and can take a long time, particularly when leave-one-temperature-out refits are needed. Its results also go under `LINEX_R_revision_bundle/results/`; start with `diagnostics.csv`, `posterior_predictive_checks.csv` and `predictive_scores.csv`. Model comparisons are only written when the required fits and validation pass. If you alter the data, priors or model code, remove the affected cached fit files before rerunning.

For the manuscript, I think the immediate question is whether both modes give defensible fits to the observations, and whether the simpler LINEX curve is adequate at RH90 given the dip in observed survival around 28 C. That dip does not by itself establish two biological temperature optima. We should decide what to put in Figure 5 after inspecting the new traces and diagnostics, and keep the LINEX contribution centred on a transparent, usable method with explicit checks for problematic fits.
