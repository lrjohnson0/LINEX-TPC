#!/usr/bin/env Rscript
# Full reviewer-response workflow for the submitted LINEX examples.
# Run from a shell with Rscript, or source() this file in R/RStudio.
# Input CSVs live in data/ beside this script. Outputs go to results/.
# Requires rstan, loo, posterior; ggplot2 is not needed.

required <- c("rstan", "loo", "posterior")
missing_packages <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop("Install missing R packages before running: ", paste(missing_packages, collapse = ", "))
}

resolve_root <- function() {
  script_name <- "run_LINEX_revision.R"
  cli <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  candidates <- if (length(cli)) sub("^--file=", "", cli) else character()
  sourced <- lapply(sys.frames(), function(frame) frame$ofile)
  candidates <- c(candidates, unlist(sourced, use.names = FALSE))
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    active <- tryCatch(rstudioapi::getActiveDocumentContext()$path,
                       error = function(e) "")
    candidates <- c(candidates, active)
  }
  candidates <- c(candidates,
                  file.path(getwd(), script_name),
                  file.path(getwd(), "LINEX_R_revision_bundle", script_name))
  for (candidate in unique(candidates[nzchar(candidates)])) {
    if (!file.exists(candidate)) next
    folder <- dirname(normalizePath(candidate))
    if (file.exists(file.path(folder, "data",
                              "binomial-juvenile-survival-data.csv"))) return(folder)
  }
  stop("Cannot locate the script beside its data/ folder. From the repository root, ",
       "run: Rscript LINEX_R_revision_bundle/run_LINEX_revision.R; ",
       "or in R: source('LINEX_R_revision_bundle/run_LINEX_revision.R').")
}
root <- resolve_root()
out <- file.path(root, "results")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "models"), showWarnings = FALSE)
dir.create(file.path(out, "fits"), showWarnings = FALSE)
dir.create(file.path(out, "figures"), showWarnings = FALSE)

# Larger values are advisable for final publication fits; cached fits can be removed
# and the script rerun. The workflow never treats an unvalidated fit as accepted.
cfg <- list(chains = 6L, unrestricted_chains = 6L,
            iter = 4000L, warmup = 2000L, adapt_delta = 0.98,
            max_treedepth = 13L, cores = min(4L, parallel::detectCores()),
            rhat_max = 1.01, ess_min = 400,
            pareto_k_max = 0.7, seed = 20260923L,
            run_exact_loo_if_needed = TRUE,
            run_simulations = TRUE, run_prior_sensitivity = TRUE,
            dataset_filter = NULL) # e.g. "stephensi_RH90" for a focused first run
rstan::rstan_options(auto_write = TRUE)
options(mc.cores = cfg$cores)
writeLines(capture.output(sessionInfo()), file.path(out, "session_info.txt"))

# The same logit-binomial likelihood and priors as the submitted code.
# For quadratic: coefficients retain the exact submitted uniform priors.
stan_common_data <- '
data {
  int<lower=1> N;
  vector[N] temp;
  int<lower=0> n[N];
  int<lower=0> y[N];
  int<lower=1> G;
  vector[G] grid;
  int<lower=-1,upper=1> a_sign;
  real<lower=0> prior_a_sd;
  real prior_Topt_mu;
  real<lower=0> prior_Topt_sd;
}
'
stan_gq <- '
generated quantities {
  vector[N] log_lik;
  vector[N] fitted_p;
  vector[G] grid_p;
  int y_rep[N];
  for (i in 1:N) {
    fitted_p[i] = inv_logit(eta[i]);
    log_lik[i] = binomial_logit_lpmf(y[i] | n[i], eta[i]);
    y_rep[i] = binomial_rng(n[i], fitted_p[i]);
  }
  for (g in 1:G) {
    grid_p[g] = inv_logit(PREDICTOR);
  }
}
'
linex_pred <- '-exp(a * (TEMP - Topt)) + a * (TEMP - Topt) + cc'
quad_pred <- 'B0 + B1 * TEMP + B2 * square(TEMP)'
make_stan <- function(kind) {
  if (kind == "unrestricted") {
    pars <- 'parameters { real a; real Topt; real<lower=0> cc; }'
    transformed_decl <- ''
    transformed_assign <- ''
    pred <- linex_pred
    prior <- 'a ~ normal(0, prior_a_sd); Topt ~ normal(prior_Topt_mu, prior_Topt_sd); cc ~ exponential(0.1);'
  } else if (kind %in% c("left", "right")) {
    # Conditional half-normal prior: original N(0,1) restricted to a sign.
    pars <- 'parameters { real<lower=0> a_mag; real Topt; real<lower=0> cc; }'
    transformed_decl <- 'real a;'
    transformed_assign <- 'a = a_sign * a_mag;'
    pred <- linex_pred
    prior <- 'a_mag ~ normal(0, prior_a_sd); Topt ~ normal(prior_Topt_mu, prior_Topt_sd); cc ~ exponential(0.1);'
  } else if (kind == "quadratic") {
    pars <- paste('parameters { real<lower=-100,upper=100> B0;',
                  'real<lower=-10,upper=10> B1;',
                  'real<lower=-1,upper=1> B2; }')
    pred <- quad_pred
    transformed_decl <- ''
    transformed_assign <- ''
    prior <- paste('B0 ~ uniform(-100,100); B1 ~ uniform(-10,10);',
                   'B2 ~ uniform(-1,1);')
  } else stop("Unknown model: ", kind)
  eta <- gsub("TEMP", "temp[i]", pred, fixed = TRUE)
  grid_eta <- gsub("TEMP", "grid[g]", pred, fixed = TRUE)
  gq <- gsub("PREDICTOR", grid_eta, stan_gq, fixed = TRUE)
  paste(stan_common_data, pars,
        'transformed parameters {', transformed_decl,
        'vector[N] eta;', transformed_assign,
        'for (i in 1:N) eta[i] =',
        eta, '; }',
        'model {', prior, 'y ~ binomial_logit(n, eta); }', gq, sep = "\n")
}

model_kinds <- c("unrestricted", "left", "right", "quadratic")
compiled <- setNames(vector("list", length(model_kinds)), model_kinds)
for (kind in model_kinds) {
  filename <- file.path(out, "models", paste0(kind, ".stan"))
  writeLines(make_stan(kind), filename)
  compiled[[kind]] <- rstan::stan_model(file = filename)
}

read_data <- function() {
  mosquito <- read.csv(file.path(root, "data", "binomial-juvenile-survival-data.csv"))
  headlee <- read.csv(file.path(root, "data", "headlee_data.csv"))
  out_data <- lapply(c(45, 75, 90), function(rh) {
    z <- mosquito[mosquito$rh == rh, ]
    data.frame(temp = z$Temp, y = z$Trait, n = z$n)
  })
  names(out_data) <- paste0("stephensi_RH", c(45, 75, 90))
  out_data$Aedes_Headlee <- data.frame(temp = headlee$interactor1temp,
                                       y = headlee$interactor1number,
                                       n = headlee$N)
  if (cfg$run_simulations) {
    simulate <- function(a, Topt, seed) {
      set.seed(seed)
      temperature <- rep(c(5, 10, 15, 20, 25, 30), each = 10)
      d <- a * (temperature - Topt)
      eta <- -exp(d) + d + 2
      data.frame(temp = temperature,
                 y = rbinom(length(temperature), 20, plogis(eta)), n = 20)
    }
    # The paper generates 10 observations per temperature at each of six settings.
    out_data$sim_left <- simulate(0.3, 21, 1983)
    out_data$sim_right <- simulate(-0.3, 14, 1984)
    out_data$sim_left_missing10 <- subset(out_data$sim_left, temp != 10)
    out_data$sim_right_missing10 <- subset(out_data$sim_right, temp != 10)
  }
  for (nm in names(out_data)) {
    z <- out_data[[nm]]
    if (anyNA(z) || any(z$n < 1) || any(z$y < 0) || any(z$y > z$n) ||
        any(z$y != round(z$y)) || any(z$n != round(z$n))) {
      stop("Invalid binomial input in ", nm)
    }
    write.csv(z, file.path(out, paste0(nm, "_input.csv")), row.names = FALSE)
  }
  out_data
}

dataset <- read_data()
if (!is.null(cfg$dataset_filter)) {
  absent <- setdiff(cfg$dataset_filter, names(dataset))
  if (length(absent)) stop("Unknown dataset(s): ", paste(absent, collapse = ", "))
  dataset <- dataset[cfg$dataset_filter]
}
data_audit <- do.call(rbind, lapply(names(dataset), function(nm) {
  z <- dataset[[nm]]
  data.frame(dataset = nm, rows = nrow(z), temperatures = length(unique(z$temp)),
             trials = sum(z$n), successes = sum(z$y),
             zero_rows = sum(z$y == 0))
}))
write.csv(data_audit, file.path(out, "data_audit.csv"), row.names = FALSE)
if ("stephensi_RH90" %in% names(dataset) &&
    !identical(as.integer(dataset$stephensi_RH90$n[
      dataset$stephensi_RH90$temp == 38]), 900L)) {
  warning("RH90 at 38 C differs from the submitted 0/900 row; verify the data")
}

original_prior <- list(a_sd = 1, Topt_mu = 25, Topt_sd = sqrt(1000))
sensitivity_prior <- list(a_sd = 0.5, Topt_mu = 28, Topt_sd = 10)
stan_data <- function(z, kind, grid, prior = original_prior) {
  list(N = nrow(z), temp = as.vector(z$temp),
       n = as.integer(z$n), y = as.integer(z$y),
       G = length(grid), grid = as.vector(grid),
       a_sign = if (kind == "left") 1L else if (kind == "right") -1L else 0L,
       prior_a_sd = prior$a_sd,
       prior_Topt_mu = prior$Topt_mu,
       prior_Topt_sd = prior$Topt_sd)
}

initials <- function(kind, chains, z) {
  if (kind == "unrestricted") {
    return(lapply(seq_len(chains), function(i) {
      positive <- i %% 2L == 1L
      list(a = if (positive) 0.25 else -0.25,
           Topt = if (positive) 30 else 23, cc = 3)
    }))
  }
  if (kind %in% c("left", "right")) {
    return(lapply(seq_len(chains), function(i)
      list(a_mag = 0.15 + i * 0.03,
           Topt = if (kind == "left") 29 + i / 3 else 21 + i / 3,
           cc = 2.5 + i / 5)))
  }
  glm_fit <- suppressWarnings(glm(cbind(y, n - y) ~ temp + I(temp^2),
                                   family = binomial(), data = z))
  b <- coef(glm_fit)
  if (any(!is.finite(b))) b <- c(-20, 1.5, -0.025)
  b <- pmax(pmin(unname(b), c(95, 9, 0.95)), c(-95, -9, -0.95))
  lapply(seq_len(chains), function(i)
    list(B0 = b[1], B1 = b[2], B2 = b[3]))
}

fit_one <- function(z, kind, key, grid, seed_offset = 0L,
                    prior = original_prior) {
  chains <- if (kind == "unrestricted") cfg$unrestricted_chains else cfg$chains
  # Encode the chain count to avoid silently reusing earlier four-chain fits.
  cache <- file.path(out, "fits", paste0(key, "_", kind, "_", chains, "chains.rds"))
  if (file.exists(cache)) return(readRDS(cache))
  fit <- rstan::sampling(compiled[[kind]], data = stan_data(z, kind, grid, prior),
                         chains = chains, iter = cfg$iter, warmup = cfg$warmup,
                         init = initials(kind, chains, z),
                         seed = cfg$seed + as.integer(seed_offset),
                         cores = min(chains, cfg$cores), refresh = 250,
                         control = list(adapt_delta = cfg$adapt_delta,
                                        max_treedepth = cfg$max_treedepth))
  saveRDS(fit, cache)
  fit
}

fit_diagnostics <- function(fit, kind, dataset_name) {
  variables <- if (kind == "quadratic") c("B0", "B1", "B2", "lp__")
               else c("a", "Topt", "cc", "lp__")
  draws <- posterior::as_draws_array(as.array(fit, pars = variables))
  summary <- as.data.frame(posterior::summarise_draws(draws))
  sampler <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
  divergences <- sum(vapply(sampler, function(x) sum(x[, "divergent__"]), numeric(1)))
  max_depth <- sum(vapply(sampler, function(x)
    sum(x[, "treedepth__"] >= cfg$max_treedepth), numeric(1)))
  valid <- all(is.finite(summary$rhat)) &&
    all(summary$rhat <= cfg$rhat_max) &&
    all(summary$ess_bulk >= cfg$ess_min) &&
    all(summary$ess_tail >= cfg$ess_min) &&
    divergences == 0 && max_depth == 0
  summary$dataset <- dataset_name
  summary$model <- kind
  summary$divergences <- divergences
  summary$max_treedepth_hits <- max_depth
  summary$diagnostics_pass <- valid
  list(summary = summary, valid = valid)
}

trace_pdf <- function(fit, kind, name) {
  vars <- if (kind == "quadratic") c("B0", "B1", "B2") else c("a", "Topt", "cc")
  pdf(file.path(out, "figures", paste0(name, "_", kind, "_traces.pdf")),
      width = 10, height = 7)
  print(rstan::traceplot(fit, pars = vars, inc_warmup = FALSE))
  dev.off()
}

predict_summary <- function(fit, name, kind, grid) {
  p <- rstan::extract(fit, pars = "grid_p")$grid_p
  q <- t(apply(p, 2, quantile, c(0.025, 0.5, 0.975)))
  data.frame(dataset = name, model = kind, temp = grid,
             lower = q[, 1], median = q[, 2], upper = q[, 3])
}

ppc_summary <- function(fit, z, name, kind) {
  rep <- rstan::extract(fit, pars = "y_rep")$y_rep
  bounds <- t(apply(rep, 2, quantile, c(0.025, 0.5, 0.975)))
  data.frame(dataset = name, model = kind, temp = z$temp,
             observed = z$y, trials = z$n,
             rep_lower = bounds[, 1], rep_median = bounds[, 2],
             rep_upper = bounds[, 3],
             outside_95 = z$y < bounds[, 1] | z$y > bounds[, 3])
}

sensitivity_summary <- function(fit, name, kind) {
  x <- rstan::extract(fit, pars = c("a", "Topt", "cc"))
  quantities <- list(a = x$a, Topt = x$Topt, Bpk = plogis(x$cc - 1))
  do.call(rbind, lapply(names(quantities), function(q) {
    v <- quantile(quantities[[q]], c(0.025, 0.5, 0.975))
    data.frame(dataset = name, model = kind, quantity = q,
               lower = unname(v[1]), median = unname(v[2]), upper = unname(v[3]))
  }))
}

score_fit <- function(fit, kind, name, diagnostics_pass) {
  if (!diagnostics_pass) {
    return(data.frame(dataset = name, model = kind, diagnostics_pass = FALSE,
                      waic = NA_real_, looic = NA_real_, elpd_loo = NA_real_,
                      max_pareto_k = NA_real_, psis_pass = FALSE))
  }
  log_lik <- loo::extract_log_lik(fit, parameter_name = "log_lik", merge_chains = FALSE)
  r_eff <- loo::relative_eff(exp(log_lik))
  w <- loo::waic(log_lik)
  l <- loo::loo(log_lik, r_eff = r_eff)
  k <- loo::pareto_k_values(l)
  data.frame(dataset = name, model = kind, diagnostics_pass = TRUE,
             waic = unname(w$estimates["waic", "Estimate"]),
             looic = unname(l$estimates["looic", "Estimate"]),
             elpd_loo = unname(l$estimates["elpd_loo", "Estimate"]),
             max_pareto_k = max(k),
             psis_pass = all(is.finite(k) & k <= cfg$pareto_k_max))
}

log_mean_exp <- function(x) {
  m <- max(x)
  m + log(mean(exp(x - m)))
}

# Exact leave-one-observation-out refits for the small real datasets if PSIS
# is unreliable. This can be computationally expensive (up to 3*N refits).
exact_loo <- function(z, kind, name, grid) {
  pointwise <- numeric(nrow(z))
  for (i in seq_len(nrow(z))) {
    fold_key <- paste0(name, "_holdout", i)
    training <- z[-i, , drop = FALSE]
    fit <- fit_one(training, kind, fold_key, grid, seed_offset = 1000L + i)
    diag <- fit_diagnostics(fit, kind, fold_key)
    if (!diag$valid) {
      warning("Exact LOO refit failed diagnostics: ", fold_key, ", ", kind)
      return(rep(NA_real_, nrow(z)))
    }
    held <- z[i, ]
    if (kind == "quadratic") {
      x <- rstan::extract(fit, pars = c("B0", "B1", "B2"))
      eta <- x$B0 + x$B1 * held$temp + x$B2 * held$temp^2
    } else {
      x <- rstan::extract(fit, pars = c("a", "Topt", "cc"))
      d <- x$a * (held$temp - x$Topt)
      eta <- -exp(d) + d + x$cc
    }
    pointwise[i] <- log_mean_exp(dbinom(held$y, held$n,
                                       plogis(eta), log = TRUE))
  }
  pointwise
}

all_diag <- list()
all_scores <- list()
all_predictions <- list()
all_ppc <- list()
prior_sensitivity <- list()
comparison <- list()
index <- 0L
for (name in names(dataset)) {
  z <- dataset[[name]]
  simulated <- grepl("^sim_", name)
  grid <- seq(min(z$temp), max(z$temp), length.out = 201)
  fits <- list()
  diags <- list()
  scores <- list()
  kinds <- if (simulated) c("unrestricted", "left", "right") else model_kinds
  message("Processing ", name, " (", nrow(z), " rows)")
  for (kind in kinds) {
    fit <- fit_one(z, kind, name, grid, seed_offset = match(name, names(dataset)))
    fits[[kind]] <- fit
    diags[[kind]] <- fit_diagnostics(fit, kind, name)
    index <- index + 1L
    all_diag[[index]] <- diags[[kind]]$summary
    write.csv(do.call(rbind, all_diag), file.path(out, "diagnostics.csv"),
              row.names = FALSE)
    trace_pdf(fit, kind, name)
    if (diags[[kind]]$valid) {
      all_predictions[[paste(name, kind)]] <- predict_summary(fit, name, kind, grid)
      all_ppc[[paste(name, kind)]] <- ppc_summary(fit, z, name, kind)
    }
    if (!simulated) {
      scores[[kind]] <- score_fit(fit, kind, name, diags[[kind]]$valid)
      all_scores[[paste(name, kind)]] <- scores[[kind]]
      write.csv(do.call(rbind, all_scores),
                file.path(out, "predictive_scores.csv"), row.names = FALSE)
    }
  }

  if (!simulated && cfg$run_prior_sensitivity) {
    for (kind in c("left", "right")) {
      sens_key <- paste0(name, "_alternative_prior")
      sens_fit <- fit_one(z, kind, sens_key, grid,
                          seed_offset = 400L + match(name, names(dataset)),
                          prior = sensitivity_prior)
      sens_diag <- fit_diagnostics(sens_fit, kind, sens_key)
      index <- index + 1L
      all_diag[[index]] <- sens_diag$summary
      write.csv(do.call(rbind, all_diag), file.path(out, "diagnostics.csv"),
                row.names = FALSE)
      trace_pdf(sens_fit, kind, sens_key)
      if (sens_diag$valid && diags[[kind]]$valid) {
        prior_sensitivity[[paste(name, kind, "original")]] <-
          transform(sensitivity_summary(fits[[kind]], name, kind), prior = "original")
        prior_sensitivity[[paste(name, kind, "alternative")]] <-
          transform(sensitivity_summary(sens_fit, name, kind), prior = "alternative")
      }
    }
  }

  if (!simulated) {
    score_table <- do.call(rbind, scores)
    candidates <- c("left", "right", "quadratic")
    candidate_scores <- score_table[match(candidates, score_table$model), ]
    # The unrestricted model is only shown diagnostically unless it actually mixes.
    if (diags$unrestricted$valid && scores$unrestricted$psis_pass) {
      candidates <- c("unrestricted", candidates)
      candidate_scores <- score_table[match(candidates, score_table$model), ]
    }
    if (all(candidate_scores$diagnostics_pass) &&
        all(candidate_scores$psis_pass)) {
      for (kind in candidates) {
        llik <- loo::extract_log_lik(fits[[kind]], parameter_name = "log_lik",
                                      merge_chains = FALSE)
        l <- loo::loo(llik, r_eff = loo::relative_eff(exp(llik)))
        comparison[[paste(name, kind)]] <- list(dataset = name,
                                                model = kind,
                                                elpd = unname(l$estimates["elpd_loo", "Estimate"]),
                                                pointwise = l$pointwise[, "elpd_loo"],
                                                method = "PSIS-LOO")
      }
    } else if (cfg$run_exact_loo_if_needed && !simulated &&
               all(candidate_scores$diagnostics_pass)) {
      message("PSIS unreliable for ", name, "; refitting every held-out row")
      exact_group <- list()
      for (kind in candidates) {
        pointwise <- exact_loo(z, kind, name, grid)
        if (all(is.finite(pointwise))) {
          exact_group[[kind]] <- list(dataset = name,
                                      model = kind,
                                      elpd = sum(pointwise),
                                      pointwise = pointwise,
                                      method = "exact row-wise LOO")
        }
      }
      if (length(exact_group) == length(candidates)) {
        for (kind in candidates)
          comparison[[paste(name, kind)]] <- exact_group[[kind]]
      } else warning("Incomplete exact LOO: no model ranking for ", name)
    } else {
      warning("No defensible model ranking for ", name,
              "; see failed diagnostics in results/diagnostics.csv")
    }
  }
}

diagnostics <- do.call(rbind, all_diag)
write.csv(diagnostics, file.path(out, "diagnostics.csv"), row.names = FALSE)
if (length(all_scores)) {
  write.csv(do.call(rbind, all_scores), file.path(out, "predictive_scores.csv"),
            row.names = FALSE)
}
if (length(all_predictions)) {
  preds <- do.call(rbind, all_predictions)
  write.csv(preds, file.path(out, "predictions.csv"), row.names = FALSE)
  recovery <- list()
  for (name in unique(preds$dataset)) {
    sub <- preds[preds$dataset == name, ]
    z <- dataset[[name]]
    pdf(file.path(out, "figures", paste0(name, "_validated_curves.pdf")),
        width = 7, height = 5)
    plot(z$temp, z$y / z$n, pch = 19, xlab = "Temperature (C)",
         ylab = "Observed success fraction / model probability",
         xlim = range(sub$temp), ylim = c(0, 1), main = name)
    colors <- c(left = "#1B9E77", right = "#D95F02",
                quadratic = "#7570B3", unrestricted = "#333333")
    for (kind in unique(sub$model)) {
      v <- sub[sub$model == kind, ]
      polygon(c(v$temp, rev(v$temp)), c(v$lower, rev(v$upper)),
              border = NA, col = grDevices::adjustcolor(colors[[kind]], alpha.f = 0.13))
      lines(v$temp, v$median, col = colors[[kind]], lwd = 2)
      if (grepl("^sim_left", name) || grepl("^sim_right", name)) {
        true_a <- if (grepl("^sim_left", name)) 0.3 else -0.3
        true_t <- if (grepl("^sim_left", name)) 21 else 14
        d <- true_a * (v$temp - true_t)
        truth <- plogis(-exp(d) + d + 2)
        recovery[[paste(name, kind)]] <- data.frame(
          dataset = name, model = kind,
          mean_absolute_error = mean(abs(v$median - truth)),
          maximum_absolute_error = max(abs(v$median - truth)))
      }
    }
    if (grepl("^sim_left", name) || grepl("^sim_right", name)) {
      true_a <- if (grepl("^sim_left", name)) 0.3 else -0.3
      true_t <- if (grepl("^sim_left", name)) 21 else 14
      d <- true_a * (unique(sub$temp) - true_t)
      lines(unique(sub$temp), plogis(-exp(d) + d + 2),
            col = "black", lty = 2, lwd = 2)
    }
    if (grepl("^sim_", name)) {
      legend("topright", legend = c(unique(sub$model), "generating curve"),
             col = c(colors[unique(sub$model)], "black"),
             lty = c(rep(1, length(unique(sub$model))), 2), lwd = 2, bty = "n")
    } else {
      legend("topright", legend = unique(sub$model),
             col = colors[unique(sub$model)], lwd = 2, bty = "n")
    }
    dev.off()
  }
  if (length(recovery)) {
    write.csv(do.call(rbind, recovery),
              file.path(out, "simulation_recovery.csv"), row.names = FALSE)
  }
}
if (length(all_ppc)) {
  write.csv(do.call(rbind, all_ppc), file.path(out, "posterior_predictive_checks.csv"),
            row.names = FALSE)
}
if (length(prior_sensitivity)) {
  write.csv(do.call(rbind, prior_sensitivity),
            file.path(out, "prior_sensitivity.csv"), row.names = FALSE)
}

if (length(comparison)) {
  comparisons <- list()
  for (name in unique(vapply(comparison, `[[`, character(1), "dataset"))) {
    group <- comparison[vapply(comparison, function(x) x$dataset == name, logical(1))]
    best <- which.max(vapply(group, `[[`, numeric(1), "elpd"))
    best_points <- group[[best]]$pointwise
    for (item in group) {
      diff <- item$pointwise - best_points
      comparisons[[paste(name, item$model)]] <- data.frame(
        dataset = name, model = item$model, method = item$method,
        elpd = item$elpd, elpd_difference_from_best = sum(diff),
        se_difference = sqrt(length(diff) * var(diff)))
    }
  }
  write.csv(do.call(rbind, comparisons), file.path(out, "model_comparisons.csv"),
            row.names = FALSE)
}
message("Complete. Start with results/diagnostics.csv; never rank a failed fit.")
