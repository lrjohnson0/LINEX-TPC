# Validation of the current bayesTPC workflow

Checked on 16 October 2026 using R 4.3.3, bayesTPC 1.1.6, NIMBLE 1.4.3, coda 0.19-4.1, posterior 1.5.0 and loo 2.6.0.

## Current execution and defaults

bayesTPC builds one model and sampler per fit. Its compiled NIMBLE sampler is reused for six chains, each with its own initialization, seed, full burn-in and reset adaptation. Two setup draws are discarded. Each fit runs in a separate R process so generated DLLs are released when the process exits. Each completed chain is cached immediately.

Prior-sensitivity fits and exact leave-one-temperature-out refits are off by default. They remain explicit options. Unreliable PSIS diagnostics do not trigger extra fits under the new defaults; the model ranking is left unwritten instead. Chain lengths, priors, likelihoods and existing cache signatures are unchanged.

## Completed checks

- Both scripts parse, and both input CSVs match the supplied data byte for byte.
- All eight example datasets pass the binomial input checks.
- All four model definitions build in bayesTPC/NIMBLE, and their log densities match independent R calculations, including the truncated-prior normalising constants.
- The current shared-compilation worker was tested on positive-conditional LINEX and quadratic fits: six chains per model, 40 iterations and 20 burn-in per chain. Each fit performed exactly one bayesTPC build, returned six 20-draw chains, and did not accumulate DLLs between chains. The parent DLL count did not increase. These are execution tests, not scientific results.
- A partial-resume test preserved five old chain caches byte for byte and fitted only the missing sixth chain using the shared-compilation worker. Existing cache signatures remain compatible.
- Earlier short actual sampling tests also covered unrestricted and negative-conditional LINEX. A complete earlier RH90 short run produced 24 chain caches. Those caches remain readable by the updated analysis.
- Controlled six-chain tests detect separated posterior regions and withhold scores from failed fits. Synthetic-draw tests exercise WAIC/LOO, posterior predictive checks, prior-sensitivity summaries and the exact leave-one-row-out comparison path.
- The RH90 plotting script runs from synthetic caches and creates all four named PDFs. Every trace PDF is a single 2-by-2 page. Negative and quadratic traces and the combined-curve panel have been rendered and visually checked.
- Archive contents, README filenames and checksum manifest are checked.

The full-length scientific analysis across all datasets has not been completed here. No short-run or synthetic result figures are included. Sampling convergence and predictive adequacy still need to be assessed from the full runs.

Interfaces checked against:
https://github.com/johnwilliamsmithjr/bayesTPC/blob/master/R/b_TPC.R
https://r-nimble.org/manual/cha-mcmc.html
