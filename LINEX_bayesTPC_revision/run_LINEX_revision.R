#!/usr/bin/env Rscript
# bayesTPC builds each model/sampler once; its NIMBLE engine runs six chains.
# Run from the folder above: Rscript LINEX_bayesTPC_revision/run_LINEX_revision.R
# Or source the same file in R/RStudio. See README.md for the RH90 figure script.

linex_root <- function() {
  args <- sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value = TRUE))
  frames <- unlist(lapply(sys.frames(), function(x) x$ofile), use.names = FALSE)
  paths <- c(args, frames, file.path(getwd(), 'run_LINEX_revision.R'),
             file.path(getwd(), 'LINEX_bayesTPC_revision', 'run_LINEX_revision.R'))
  for (p in unique(paths[nzchar(paths)])) {
    if (file.exists(p)) {
      root <- dirname(normalizePath(p))
      if (file.exists(file.path(root, 'data', 'headlee_data.csv'))) return(root)
    }
  }
  stop('Cannot find the script beside its data/ directory; source or run the full script.')
}

linex_config <- function() {
  defaults <- list(chains = 6L, niter = 25000L, burn = 15000L,
    sampler = 'RW', quadratic_sampler = 'AF_slice', seed = 20261016L,
    rhat_max = 1.01, ess_min = 400, pareto_k_max = 0.7,
    run_simulations = TRUE, run_prior_sensitivity = FALSE,
    run_exact_loo_if_needed = FALSE, dataset_filter = NULL,
    predictive_draws = 2000L)
  cfg <- utils::modifyList(defaults, getOption('linex_revision', list()))
  if (cfg$chains != 6L) stop('This revision requires six chains per fit.')
  if (cfg$burn >= cfg$niter) stop('burn must be less than niter.')
  cfg
}

linex_read_data <- function(root, cfg) {
  mosquito <- read.csv(file.path(root, 'data', 'binomial-juvenile-survival-data.csv'))
  headlee <- read.csv(file.path(root, 'data', 'headlee_data.csv'))
  sets <- lapply(c(45, 75, 90), function(rh) {
    z <- mosquito[mosquito$rh == rh, ]
    data.frame(Temp = z$Temp, Trait = z$Trait, n = z$n)
  })
  names(sets) <- paste0('stephensi_RH', c(45, 75, 90))
  sets$Aedes_Headlee <- data.frame(Temp = headlee$interactor1temp,
                                 Trait = headlee$interactor1number, n = headlee$N)
  if (cfg$run_simulations) {
    simulate <- function(a, Topt, seed) {
      set.seed(seed)
      temp <- rep(c(5, 10, 15, 20, 25, 30), each = 10)
      d <- a * (temp - Topt)
      z <- data.frame(Temp = temp, Trait = rbinom(length(temp), 20,
                          plogis(-exp(d) + d + 2)), n = 20)
      attr(z, 'truth') <- c(a = a, Topt = Topt, cc = 2)
      z
    }
    sets$sim_left <- simulate(0.3, 21, 1983)
    sets$sim_right <- simulate(-0.3, 14, 1984)
    sets$sim_left_missing10 <- sets$sim_left[sets$sim_left$Temp != 10, ]
    sets$sim_right_missing10 <- sets$sim_right[sets$sim_right$Temp != 10, ]
    attr(sets$sim_left_missing10, 'truth') <- attr(sets$sim_left, 'truth')
    attr(sets$sim_right_missing10, 'truth') <- attr(sets$sim_right, 'truth')
  }
  if (!is.null(cfg$dataset_filter)) {
    unknown <- setdiff(cfg$dataset_filter, names(sets))
    if (length(unknown)) stop('Unknown dataset: ', paste(unknown, collapse = ', '))
    sets <- sets[cfg$dataset_filter]
  }
  for (nm in names(sets)) {
    z <- sets[[nm]]
    if (anyNA(z) || any(!is.finite(as.matrix(z))) || any(z$n < 1) ||
        any(z$Trait < 0 | z$Trait > z$n) || any(z$Trait != round(z$Trait)) ||
        any(z$n != round(z$n))) stop('Invalid binomial data in ', nm)
  }
  sets
}

# Conditional models retain the original normal prior, truncated to a sign.
# They are different analysis targets from unrestricted chains started at a sign.
linex_model <- function(kind, prior, name) {
  if (kind == 'quadratic') {
    return(bayesTPC::specify_binomial_model(name,
      parameters = c(B0 = 'dunif(-100,100)', B1 = 'dunif(-10,10)', B2 = 'dunif(-1,1)'),
      formula = expression(B0 + B1 * Temp + B2 * Temp^2)))
  }
  a_prior <- sprintf('dnorm(0, %.17g)', 1 / prior$a_sd^2)
  if (kind == 'positive') a_prior <- paste0('T(', a_prior, ', 0, Inf)')
  if (kind == 'negative') a_prior <- paste0('T(', a_prior, ', -Inf, 0)')
  bayesTPC::specify_binomial_model(name,
    parameters = c(a = a_prior,
      Topt = sprintf('dnorm(%.17g, %.17g)', prior$Topt_mu, 1 / prior$Topt_sd^2),
      cc = 'dexp(0.1)'),
    formula = expression(-exp(a * (Temp - Topt)) + a * (Temp - Topt) + cc))
}

linex_initials <- function(kind, z, chain) {
  if (kind == 'quadratic') {
    b <- tryCatch(suppressWarnings(coef(glm(cbind(Trait, n - Trait) ~
                       Temp + I(Temp^2), family = binomial(), data = z))),
                  error = function(e) rep(NA_real_, 3L))
    if (any(!is.finite(b))) b <- c(-20, 1.5, -0.025)
    b <- pmax(pmin(unname(b), c(95, 9, 0.95)), c(-95, -9, -0.95))
    return(list(B0 = b[1], B1 = b[2], B2 = b[3]))
  }
  positive <- kind == 'positive' || (kind == 'unrestricted' && chain %% 2L == 1L)
  peak <- z$Temp[which.max(z$Trait / z$n)]
  list(a = (if (positive) 1 else -1) * (0.15 + 0.025 * chain),
       Topt = peak + (if (positive) 1 else -1) * diff(range(z$Temp)) / 6 +
                     (chain - 3.5) * 0.25,
       cc = 2.5 + chain * 0.1)
}

linex_find_samples <- function(fit) {
  x <- fit$samples
  if (inherits(x, 'mcmc.list')) {
    if (length(x) != 1L) stop('Expected one chain from each separately initialized b_TPC call.')
    return(x[[1L]])
  }
  if (inherits(x, 'mcmc')) return(x)
  stop('Cannot find a coda mcmc sample in bayesTPC output; inspect str(fit$samples).')
}

linex_atomic_save <- function(x, path) {
  tmp <- tempfile('linex-cache-', tmpdir = dirname(path), fileext = '.pending')
  on.exit(unlink(tmp))
  saveRDS(x, tmp)
  if (!file.rename(tmp, path)) stop('Cannot finish writing cache: ', path)
}

# Each model's compilation lives in a fresh R process. Ending it releases
# NIMBLE's generated DLLs; rm()/gc() alone does not reliably release them.
linex_chain_worker <- function(request_path) {
  request <- readRDS(request_path)
  .libPaths(request$libpaths)
  suppressPackageStartupMessages(library(nimble))
  suppressPackageStartupMessages(library(bayesTPC))
  for (pkg in c('bayesTPC', 'nimble')) {
    if (!identical(as.character(packageVersion(pkg)), request$signature[[pkg]]))
      stop('Worker package version differs from the parent: ', pkg)
  }
  spec <- linex_model(request$kind, request$prior,
                       paste0('linex_revision_', request$kind, '_', request$key))
  if (!is.null(request$jobs)) {
    # bayesTPC builds the model and configures the SAME sampler once. These
    # two setup draws are discarded; they are not any of the retained chains.
    first <- request$jobs[[1]]
    set.seed(first$seed)
    fit <- bayesTPC::b_TPC(data = request$data, model = spec, nchains = 1L,
      inits = first$init, samplerType = request$sampler,
      burn = 0L, niter = 2L, store_nimble_model = TRUE)
    completed <- character()
    for (job in request$jobs) {
      message('Running chain ', job$chain, ' using the compiled bayesTPC sampler')
      # runMCMC resets adaptive samplers by default; each chain receives its
      # own initial values, seed, full warmup and retained-draw schedule.
      samples <- nimble::runMCMC(fit$mcmc, niter = request$niter,
        nburnin = request$burn, nchains = 1L, inits = job$init,
        setSeed = job$seed, samplesAsCodaMCMC = TRUE,
        summary = FALSE, WAIC = FALSE)
      x <- linex_find_samples(list(samples = samples))
      params <- if (request$kind == 'quadratic') c('B0', 'B1', 'B2') else c('a', 'Topt', 'cc')
      x <- coda::mcmc(as.matrix(x)[, params, drop = FALSE],
                      start = start(x), thin = coda::thin(x))
      if (any(!is.finite(as.matrix(x)))) stop('Nonfinite samples in worker.')
      linex_atomic_save(list(samples = x, signature = job$signature,
        worker_pid = Sys.getpid(), worker_DLLs = length(getLoadedDLLs())), job$path)
      completed <- c(completed, job$path)
    }
    linex_atomic_save(list(signature = request$signature, completed = completed,
                          worker_pid = Sys.getpid()), request$result)
    return(invisible(NULL))
  }
  set.seed(request$seed)
  fit <- bayesTPC::b_TPC(data = request$data, model = spec, nchains = 1L,
    inits = request$init, samplerType = request$sampler,
    burn = request$burn, niter = request$niter, store_nimble_model = FALSE)
  x <- linex_find_samples(fit)
  params <- if (request$kind == 'quadratic') c('B0', 'B1', 'B2') else c('a', 'Topt', 'cc')
  x <- coda::mcmc(as.matrix(x)[, params, drop = FALSE],
                  start = start(x), thin = coda::thin(x))
  if (any(!is.finite(as.matrix(x)))) stop('Nonfinite samples in worker.')
  linex_atomic_save(list(samples = x, signature = request$signature,
    worker_pid = Sys.getpid(), worker_DLLs = length(getLoadedDLLs())), request$result)
}

linex_run_worker <- function(request, script, log_path) {
  if (is.null(script) || !file.exists(script))
    stop('Cannot find run_LINEX_revision.R for the independent chain worker.')
  request_file <- tempfile('linex-request-', fileext = '.rds')
  result_file <- tempfile('linex-result-', fileext = '.rds')
  on.exit(unlink(c(request_file, result_file)))
  request$result <- result_file
  request$libpaths <- .libPaths()
  saveRDS(request, request_file)
  status <- system2(file.path(R.home('bin'), 'Rscript'),
    args = c('--vanilla', shQuote(normalizePath(script)),
             '--linex-worker', shQuote(request_file)),
    stdout = log_path, stderr = log_path)
  if (!identical(as.integer(status), 0L) || !file.exists(result_file)) {
    details <- if (file.exists(log_path)) tail(readLines(log_path, warn = FALSE), 12) else character()
    stop('Chain worker failed. Completed caches are retained. Log: ', log_path,
         '\n', paste(details, collapse = '\n'))
  }
  saved <- readRDS(result_file)
  if (!identical(saved$signature, request$signature)) stop('Worker result settings do not match.')
  saved
}

linex_fit <- function(z, kind, key, prior, cfg, out, seed_offset = 0L) {
  # Compile once per model; keep independent starts and per-chain checkpoints.
  settings <- list(data = z[, c('Temp', 'Trait', 'n')], kind = kind, prior = prior,
    niter = cfg$niter, burn = cfg$burn,
    sampler = if (kind == 'quadratic') cfg$quadratic_sampler else cfg$sampler,
    bayesTPC = as.character(packageVersion('bayesTPC')),
    nimble = as.character(packageVersion('nimble')), schema = 1L)
  chains <- vector('list', cfg$chains)
  jobs <- list()
  for (j in seq_len(cfg$chains)) {
    path <- file.path(out, 'fits', paste0(key, '_', kind, '_chain', j, '.rds'))
    init <- linex_initials(kind, z, j)
    seed <- as.integer(cfg$seed + seed_offset + j)
    signature <- c(settings, list(init = init, seed = seed))
    if (file.exists(path)) {
      saved <- readRDS(path)
      if (!identical(saved$signature, signature))
        stop('Cached fit settings changed: ', path,
             '. Move that cache file aside before rerunning.')
      chains[[j]] <- saved$samples
      message('Loaded ', key, ' / ', kind, ' / chain ', j)
    } else {
      jobs[[length(jobs) + 1L]] <- list(chain = j, path = file.path(normalizePath(dirname(path)), basename(path)),
                                        init = init, seed = seed, signature = signature)
    }
  }
  if (length(jobs)) {
    log_path <- file.path(out, 'fits', paste0(key, '_', kind, '_worker.log'))
    message('Fitting ', key, ' / ', kind, ': ', length(jobs),
            ' missing chains, one compilation setup. Log: ', log_path)
    linex_run_worker(list(data = z, kind = kind, key = key, prior = prior,
      jobs = jobs, sampler = settings$sampler, burn = cfg$burn,
      niter = cfg$niter, signature = settings), cfg$worker_script, log_path)
    for (job in jobs) {
      saved <- readRDS(job$path)
      if (!identical(saved$signature, job$signature)) stop('Worker cache settings do not match.')
      chains[[job$chain]] <- saved$samples
    }
  }
  do.call(coda::mcmc.list, chains)
}

linex_eta <- function(x, temp, kind) {
  if (kind == 'quadratic') return(outer(x[, 'B0'], rep(1, length(temp))) +
                    outer(x[, 'B1'], temp) + outer(x[, 'B2'], temp^2))
  d <- outer(x[, 'a'], temp) - x[, 'a'] * x[, 'Topt']
  -exp(pmin(d, 700)) + d + x[, 'cc']
}
linex_softplus <- function(x) pmax(x, 0) + log1p(exp(-abs(x)))
linex_loglik <- function(x, z, kind) {
  eta <- linex_eta(x, z$Temp, kind)
  lp <- -linex_softplus(-eta)
  lq <- -linex_softplus(eta)
  ll <- matrix(0, nrow(x), nrow(z))
  for (i in seq_len(nrow(z))) ll[, i] <- lchoose(z$n[i], z$Trait[i]) +
    (if (z$Trait[i] == 0) 0 else z$Trait[i] * lp[, i]) +
    (if (z$Trait[i] == z$n[i]) 0 else (z$n[i] - z$Trait[i]) * lq[, i])
  ll
}
linex_logprior <- function(x, kind, prior) {
  if (kind == 'quadratic') return(rep(-log(200) - log(20) - log(2), nrow(x)))
  dnorm(x[, 'a'], 0, prior$a_sd, log = TRUE) +
    dnorm(x[, 'Topt'], prior$Topt_mu, prior$Topt_sd, log = TRUE) +
    dexp(x[, 'cc'], 0.1, log = TRUE)
}
linex_diagnostics <- function(chains, z, kind, prior, dataset, cfg) {
  mats <- lapply(chains, as.matrix)
  pars <- colnames(mats[[1]])
  a <- array(NA_real_, c(nrow(mats[[1]]), length(mats), length(pars) + 1L),
             dimnames = list(NULL, NULL, c(pars, 'log_posterior')))
  for (j in seq_along(mats)) {
    a[, j, pars] <- mats[[j]]
    a[, j, 'log_posterior'] <- rowSums(linex_loglik(mats[[j]], z, kind)) +
                                linex_logprior(mats[[j]], kind, prior)
  }
  tab <- as.data.frame(posterior::summarise_draws(posterior::as_draws_array(a),
              'rhat', 'ess_bulk', 'ess_tail'))
  valid <- all(is.finite(as.matrix(tab[, c('rhat', 'ess_bulk', 'ess_tail')]))) &&
    all(tab$rhat <= cfg$rhat_max) && all(tab$ess_bulk >= cfg$ess_min) &&
    all(tab$ess_tail >= cfg$ess_min)
  summary <- data.frame(dataset = dataset, model = kind,
    max_Rhat = max(tab$rhat), min_bulk_ESS = min(tab$ess_bulk),
    min_tail_ESS = min(tab$ess_tail), diagnostics_pass = valid)
  list(valid = valid, summary = summary,
       parameters = cbind(dataset = dataset, model = kind, tab))
}
linex_array_loglik <- function(chains, z, kind) {
  out <- array(NA_real_, c(nrow(chains[[1]]), length(chains), nrow(z)))
  for (j in seq_along(chains)) out[, j, ] <- linex_loglik(as.matrix(chains[[j]]), z, kind)
  out
}
linex_score <- function(chains, z, kind, key, valid, cfg) {
  row <- data.frame(dataset = key, model = kind, diagnostics_pass = valid,
    WAIC = NA_real_, WAIC_SE = NA_real_, p_WAIC = NA_real_,
    WAIC_variance_warning = NA, elpd_loo = NA_real_, loo_SE = NA_real_,
    max_pareto_k = NA_real_, psis_pass = FALSE)
  if (!valid) return(list(summary = row, loo = NULL))
  ll <- linex_array_loglik(chains, z, kind)
  result <- tryCatch({
    w <- loo::waic(ll)
    # Scaling per row prevents underflow without changing relative ESS.
    scaled <- ll
    for (i in seq_len(nrow(z))) scaled[, , i] <- exp(ll[, , i] - max(ll[, , i]))
    r_eff <- loo::relative_eff(scaled)
    l <- loo::loo(ll, r_eff = r_eff)
    row$WAIC <- w$estimates['waic', 'Estimate']
    row$WAIC_SE <- w$estimates['waic', 'SE']
    row$p_WAIC <- w$estimates['p_waic', 'Estimate']
    row$WAIC_variance_warning <- any(w$pointwise[, 'p_waic'] > 0.4)
    row$elpd_loo <- l$estimates['elpd_loo', 'Estimate']
    row$loo_SE <- l$estimates['elpd_loo', 'SE']
    row$max_pareto_k <- max(loo::pareto_k_values(l))
    row$psis_pass <- all(is.finite(loo::pareto_k_values(l))) &&
                        all(loo::pareto_k_values(l) <= cfg$pareto_k_max)
    list(summary = row, loo = l)
  }, error = function(e) {
    warning('Score unavailable for ', key, '/', kind, ': ', conditionMessage(e))
    list(summary = row, loo = NULL)
  })
  result
}
linex_pool <- function(chains, n = NULL) {
  x <- do.call(rbind, lapply(chains, as.matrix))
  if (!is.null(n) && nrow(x) > n) x <- x[unique(round(seq(1, nrow(x), length.out = n))), , drop = FALSE]
  x
}
linex_pdf <- function(path, width = 10, height = 8) {
  if (isTRUE(capabilities('cairo'))) grDevices::cairo_pdf(path, width, height)
  else grDevices::pdf(path, width, height, useDingbats = FALSE)
}
linex_trace <- function(chains, kind, label, path) {
  linex_pdf(path)
  on.exit(dev.off())
  par(mfrow = c(2, 2), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))
  pars <- colnames(chains[[1]])
  cols <- grDevices::hcl.colors(length(chains), 'Dark 3')
  for (v in c(pars, 'p30')) {
    ys <- lapply(chains, function(x) if (v == 'p30')
      as.vector(plogis(linex_eta(as.matrix(x), 30, kind))) else as.matrix(x)[, v])
    plot(NA, xlim = c(1, max(lengths(ys))), ylim = extendrange(range(unlist(ys))),
      xlab = 'Stored iteration', ylab = if (v == 'p30') 'Survival at 30 C' else v,
      main = if (v == 'p30') 'Predicted survival at 30 C' else v)
    for (j in seq_along(ys)) {
      ix <- unique(round(seq(1, length(ys[[j]]), length.out = min(1500, length(ys[[j]])))))
      lines(ix, ys[[j]][ix], col = adjustcolor(cols[j], alpha.f = 0.75), lwd = 0.6)
    }
    if (v == pars[1]) legend('topright', paste('Chain', seq_along(chains)),
                             col = cols, lty = 1, cex = 0.65, bg = 'white')
  }
  mtext(label, outer = TRUE, cex = 0.85)
}
linex_prediction <- function(chains, z, kind, key, cfg) {
  x <- linex_pool(chains, cfg$predictive_draws)
  grid <- seq(min(z$Temp), max(z$Temp), length.out = 201)
  p <- plogis(linex_eta(x, grid, kind))
  b <- apply(p, 2, quantile, probs = c(0.025, 0.5, 0.975))
  obs_p <- plogis(linex_eta(x, z$Temp, kind))
  replicated <- vapply(seq_len(nrow(z)), function(i)
    rbinom(nrow(x), z$n[i], obs_p[, i]) / z$n[i], numeric(nrow(x)))
  rb <- apply(replicated, 2, quantile, probs = c(0.025, 0.5, 0.975))
  # Compare discrepancy in replicate binomial counts with observed counts.
  denom <- pmax(sweep(obs_p * (1 - obs_p), 2, z$n, '*'), 1e-10)
  expected <- sweep(obs_p, 2, z$n, '*')
  residual <- sweep(expected, 2, z$Trait, '-')
  rep_residual <- sweep(replicated, 2, z$n, '*') - expected
  pval <- mean(rowSums(rep_residual^2 / denom) >= rowSums(residual^2 / denom))
  list(grid = data.frame(dataset = key, model = kind, Temp = grid,
                         lower = b[1, ], median = b[2, ], upper = b[3, ]),
       rows = data.frame(dataset = key, model = kind, Temp = z$Temp,
         observed = z$Trait / z$n, lower_predictive = rb[1, ],
         upper_predictive = rb[3, ], inside_95 = z$Trait / z$n >= rb[1, ] &
                                  z$Trait / z$n <= rb[3, ]),
       ppc = data.frame(dataset = key, model = kind, bayesian_p = pval,
         rows_inside_95 = sum(z$Trait / z$n >= rb[1, ] & z$Trait / z$n <= rb[3, ]),
         rows_total = nrow(z)))
}
linex_exact_loo <- function(z, kind, key, prior, cfg, out) {
  values <- rep(NA_real_, nrow(z))
  for (i in seq_len(nrow(z))) {
    fold_key <- paste0(key, '_loo', i)
    fit <- linex_fit(z[-i, , drop = FALSE], kind, fold_key, prior, cfg, out, 50000L + 100L * i)
    diag <- linex_diagnostics(fit, z[-i, , drop = FALSE], kind, prior, fold_key, cfg)
    write.csv(diag$parameters, file.path(out, 'fits', paste0(fold_key, '_', kind, '_diagnostics.csv')), row.names = FALSE)
    if (!diag$valid) next
    ll <- as.vector(linex_loglik(linex_pool(fit), z[i, , drop = FALSE], kind))
    m <- max(ll)
    values[i] <- m + log(mean(exp(ll - m)))
  }
  values
}
linex_write_rows <- function(rows, path) {
  if (length(rows)) write.csv(do.call(rbind, rows), path, row.names = FALSE)
}

linex_main <- function(root = linex_root(), cfg = linex_config()) {
  if (is.null(cfg$worker_script)) cfg$worker_script <- file.path(root, 'run_LINEX_revision.R')
  if (!requireNamespace('nimble', quietly = TRUE)) stop('Install the R package nimble.')
  suppressPackageStartupMessages(library(nimble))
  needed <- c('bayesTPC', 'nimble', 'coda', 'posterior', 'loo')
  missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop('Install missing R packages: ', paste(missing, collapse = ', '))
  suppressPackageStartupMessages(library(bayesTPC))
  # Keep these outputs separate from earlier results made with other software.
  out <- file.path(root, 'results', 'bayesTPC')
  for (d in c(out, file.path(out, 'fits'), file.path(out, 'figures')))
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
  writeLines(capture.output(sessionInfo()), file.path(out, 'session_info.txt'))
  saveRDS(cfg, file.path(out, 'settings.rds'))
  sets <- linex_read_data(root, cfg)
  # Remove derived tables from an earlier run; leave resumable fit caches.
  # Otherwise a failed diagnostic could leave an old prediction/ranking visible.
  unlink(file.path(out, c('data_audit.csv', 'diagnostics.csv',
    'parameter_diagnostics.csv', 'predictive_scores.csv', 'predictions.csv',
    'predictive_checks_by_temperature.csv', 'posterior_predictive_checks.csv',
    'prior_sensitivity.csv', 'simulation_recovery.csv', 'model_comparisons.csv')))
  prior <- list(a_sd = 1, Topt_mu = 25, Topt_sd = sqrt(1000))
  sens_prior <- list(a_sd = 0.5, Topt_mu = 28, Topt_sd = 10)
  audit <- lapply(names(sets), function(nm) {
    z <- sets[[nm]]
    write.csv(z, file.path(out, paste0(nm, '_input.csv')), row.names = FALSE)
    data.frame(dataset = nm, rows = nrow(z), temperatures = length(unique(z$Temp)),
               trials = sum(z$n), successes = sum(z$Trait))
  })
  linex_write_rows(audit, file.path(out, 'data_audit.csv'))
  all_diag <- all_parameters <- scores <- predictions <- predictive_rows <- ppcs <- sensitivity <- recovery <- comparisons <- list()
  for (key in names(sets)) {
    z <- sets[[key]]
    unlink(file.path(out, 'figures', paste0(key, '_validated_curves.pdf')))
    simulated <- grepl('^sim_', key)
    kinds <- c('unrestricted', 'positive', 'negative', if (!simulated) 'quadratic')
    fits <- diags <- scored <- list()
    for (kind in kinds) {
      fits[[kind]] <- linex_fit(z, kind, key, prior, cfg, out,
                                seed_offset = 1000L * match(key, names(sets)))
      d <- linex_diagnostics(fits[[kind]], z, kind, prior, key, cfg)
      diags[[kind]] <- d
      all_diag[[paste(key, kind)]] <- d$summary
      all_parameters[[paste(key, kind)]] <- d$parameters
      linex_trace(fits[[kind]], kind, paste(key, kind, '(bayesTPC)'),
                   file.path(out, 'figures', paste0(key, '_', kind, '_traces.pdf')))
      if (kind != 'quadratic') {
        mode_rows <- lapply(seq_along(fits[[kind]]), function(j) {
          a <- as.matrix(fits[[kind]][[j]])[, 'a']
          data.frame(dataset = key, model = kind, chain = j, mean_a = mean(a),
                      fraction_a_positive = mean(a > 0), crosses_sign = min(a) < 0 && max(a) > 0)
        })
        linex_write_rows(mode_rows, file.path(out, paste0(key, '_', kind, '_chain_modes.csv')))
      }
      if (!simulated) {
        scored[[kind]] <- linex_score(fits[[kind]], z, kind, key, d$valid, cfg)
        scores[[paste(key, kind)]] <- scored[[kind]]$summary
      }
      if (d$valid) {
        set.seed(cfg$seed + match(kind, kinds) + match(key, names(sets)))
        pred <- linex_prediction(fits[[kind]], z, kind, key, cfg)
        predictions[[paste(key, kind)]] <- pred$grid
        predictive_rows[[paste(key, kind)]] <- pred$rows
        ppcs[[paste(key, kind)]] <- pred$ppc
        if (simulated) {
          x <- linex_pool(fits[[kind]])
          truth <- attr(z, 'truth')
          for (v in names(truth)) {
            q <- quantile(x[, v], c(0.025, 0.5, 0.975))
            recovery[[paste(key, kind, v)]] <- data.frame(dataset = key, model = kind,
              parameter = v, truth = truth[v], lower = q[1], median = q[2], upper = q[3],
              covers_truth = truth[v] >= q[1] && truth[v] <= q[3])
          }
        }
      }
      linex_write_rows(all_diag, file.path(out, 'diagnostics.csv'))
      linex_write_rows(all_parameters, file.path(out, 'parameter_diagnostics.csv'))
      linex_write_rows(scores, file.path(out, 'predictive_scores.csv'))
      linex_write_rows(predictions, file.path(out, 'predictions.csv'))
      linex_write_rows(predictive_rows, file.path(out, 'predictive_checks_by_temperature.csv'))
      linex_write_rows(ppcs, file.path(out, 'posterior_predictive_checks.csv'))
    }
    if (!simulated && cfg$run_prior_sensitivity) {
      for (kind in c('positive', 'negative')) {
        sens_key <- paste0(key, '_alternative_prior')
        sf <- linex_fit(z, kind, sens_key, sens_prior, cfg, out, 30000L + 100L * match(key, names(sets)))
        sd <- linex_diagnostics(sf, z, kind, sens_prior, sens_key, cfg)
        all_diag[[paste(sens_key, kind)]] <- sd$summary
        all_parameters[[paste(sens_key, kind)]] <- sd$parameters
        linex_trace(sf, kind, paste(key, kind, 'alternative prior'),
                     file.path(out, 'figures', paste0(sens_key, '_', kind, '_traces.pdf')))
        if (sd$valid && diags[[kind]]$valid) {
          for (which_prior in c('original', 'alternative')) {
            x <- linex_pool(if (which_prior == 'original') fits[[kind]] else sf)
            for (v in colnames(x)) {
              q <- quantile(x[, v], c(0.025, 0.5, 0.975))
              sensitivity[[paste(key, kind, which_prior, v)]] <- data.frame(dataset = key,
                model = kind, prior = which_prior, parameter = v,
                lower = q[1], median = q[2], upper = q[3])
            }
          }
        }
      }
    }
    # Comparisons are conditional sign models vs quadratic. The unrestricted
    # fit is retained as a mixing check; it is not ranked if its chains split.
    if (!simulated) {
      candidates <- c('positive', 'negative', 'quadratic')
      if (all(vapply(diags[candidates], function(d) d$valid, logical(1)))) {
        pointwise <- NULL
        if (all(vapply(scored[candidates], function(s) s$summary$psis_pass, logical(1)))) {
          pointwise <- lapply(scored[candidates], function(s) s$loo$pointwise[, 'elpd_loo'])
          method <- 'PSIS-LOO'
        } else if (cfg$run_exact_loo_if_needed) {
          pointwise <- setNames(lapply(candidates, function(k)
            linex_exact_loo(z, k, key, prior, cfg, out)), candidates)
          method <- 'exact leave-one-temperature-row-out'
        }
        if (!is.null(pointwise) && all(vapply(pointwise, function(p) all(is.finite(p)), logical(1)))) {
          elpd <- vapply(pointwise, sum, numeric(1))
          best <- names(which.max(elpd))
          for (k in candidates) {
            difference <- pointwise[[k]] - pointwise[[best]]
            comparisons[[paste(key, k)]] <- data.frame(dataset = key, model = k,
              method = method, elpd = elpd[k], best_model = best,
              elpd_difference = sum(difference),
              difference_SE = sqrt(length(difference) * var(difference)))
          }
        } else message('No validated LOO ranking for ', key, '; inspect diagnostics and scores.')
      }
    }
    valid_names <- names(diags)[vapply(diags, function(d) d$valid, logical(1))]
    if (length(valid_names)) {
      linex_pdf(file.path(out, 'figures', paste0(key, '_validated_curves.pdf')), 8, 6)
      plot(z$Temp, z$Trait / z$n, ylim = c(0, 1), xlab = 'Temperature (C)',
           ylab = 'Survival probability', main = key, pch = 19)
      cols <- c(unrestricted = 'black', positive = '#1B9E77', negative = '#D95F02', quadratic = '#7570B3')
      for (k in valid_names) {
        p <- predictions[[paste(key, k)]]
        polygon(c(p$Temp, rev(p$Temp)), c(p$lower, rev(p$upper)),
                 col = adjustcolor(cols[k], alpha.f = 0.15), border = NA)
        lines(p$Temp, p$median, col = cols[k], lwd = 2)
      }
      points(z$Temp, z$Trait / z$n, pch = 19)
      legend('topright', valid_names, col = cols[valid_names], lty = 1, cex = 0.8)
      dev.off()
    }
    linex_write_rows(all_diag, file.path(out, 'diagnostics.csv'))
    linex_write_rows(all_parameters, file.path(out, 'parameter_diagnostics.csv'))
    linex_write_rows(sensitivity, file.path(out, 'prior_sensitivity.csv'))
    linex_write_rows(recovery, file.path(out, 'simulation_recovery.csv'))
    linex_write_rows(comparisons, file.path(out, 'model_comparisons.csv'))
  }
  message('Completed bayesTPC workflow. Outputs: ', out)
  invisible(list(output = out, diagnostics = all_diag, config = cfg))
}

linex_worker_args <- commandArgs(TRUE)
if (length(linex_worker_args) >= 2L && identical(linex_worker_args[1], '--linex-worker')) {
  linex_chain_worker(linex_worker_args[2])
} else if (!identical(Sys.getenv('LINEX_REVISION_NO_RUN'), '1')) linex_main()
