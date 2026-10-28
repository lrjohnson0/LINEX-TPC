#!/usr/bin/env Rscript
# RH90 revision check using bayesTPC/NIMBLE. The equal-mode prediction below
# assigns 50% to each separately sampled LINEX mode for illustration only.
# Its weights are not estimated posterior mode probabilities or LOO stacking.

needed <- c("nimble", "coda")
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "))

# NIMBLE-generated model methods can resolve helpers (including
# getNimbleOption) through the attached package search path. Loading its
# namespace with requireNamespace() alone is insufficient for some releases.
library(nimble)
if (!requireNamespace("bayesTPC", quietly = TRUE)) stop("Missing R package: bayesTPC")
library(bayesTPC)

script_paths <- c(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)),
                  unlist(lapply(sys.frames(), function(x) x$ofile), use.names = FALSE),
                  file.path(getwd(), "plot_RH90_bayesTPC_six_chains.R"),
                  file.path(getwd(), "LINEX_R_revision_bundle",
                            "plot_RH90_bayesTPC_six_chains.R"))
root <- NULL
for (path in unique(script_paths[nzchar(script_paths)])) {
  if (file.exists(path) &&
      file.exists(file.path(dirname(normalizePath(path)), "data",
                            "binomial-juvenile-survival-data.csv"))) {
    root <- dirname(normalizePath(path))
    break
  }
}
if (is.null(root)) stop("Run the script from the bundle directory or its repository root")
output <- file.path(root, "results", "figures")
dir.create(output, recursive = TRUE, showWarnings = FALSE)
# Remove only obsolete plots produced by earlier versions of this script.
unlink(file.path(output, c("RH90_bayesTPC_standard_per_chain_panels.pdf",
                           "RH90_bayesTPC_combined_plot_DIAGNOSTIC_ONLY.pdf")))

dat <- read.csv(file.path(root, "data", "binomial-juvenile-survival-data.csv"))
dat <- dat[dat$rh == 90, , drop = FALSE]
stopifnot(nrow(dat) > 0, all(dat$Trait <= dat$n))

bayesTPC::reset_models()
spec <- bayesTPC::specify_binomial_model(
  "linex_binom_six_chain_figure",
  parameters = c(cc = "dexp(0.1)", a = "dnorm(0, 1)",
                 Topt = "dnorm(25, 1/1000)"),
  formula = expression(-exp(-a * (Topt - Temp)) + a * (Temp - Topt) + cc)
)

# Use the same priors and sampler choices as the submitted bayesTPC code.
# The six chains may settle in either mode; their allocation is not evidence
# about the modes' relative posterior probabilities if they fail to mix.
set.seed(20260923)
cached_fit <- function(filename, fit_fun) {
  path <- file.path(root, "results", filename)
  if (file.exists(path)) return(readRDS(path))
  fit <- fit_fun()
  saveRDS(fit, path)
  fit
}
linex_fit <- cached_fit("RH90_bayesTPC_LINEX_6chains.rds", function()
  bayesTPC::b_TPC(data = dat, model = spec, nchains = 6L,
                  burn = 15000L, niter = 25000L))
quadratic_fit <- cached_fit("RH90_bayesTPC_quadratic_6chains.rds", function()
  bayesTPC::b_TPC(
    data = dat, model = "binomial_glm_quad", nchains = 6L,
    samplerType = "AF_slice",
    priors = list(B0 = "dunif(-100,100)", B1 = "dunif(-10,10)",
                  B2 = "dunif(-1,1)"),
    inits = list(B0 = 10, B1 = 1.2, B2 = -0.02),
    burn = 5000L, niter = 15000L))

# bayesTPC stores chains as a coda mcmc.list; locate it without assuming a
# particular element name, which may vary across bayesTPC versions.
find_chains <- function(x) {
  if (inherits(x, "mcmc.list")) return(x)
  if (inherits(x, "mcmc")) return(coda::mcmc.list(x))
  if (is.list(x)) {
    for (part in x) {
      found <- find_chains(part)
      if (!is.null(found)) return(found)
    }
  }
  NULL
}
chains <- find_chains(linex_fit)
if (is.null(chains) || length(chains) != 6L)
  stop("Expected six LINEX chains as a coda mcmc.list in the bayesTPC fit")
draws <- lapply(chains, as.matrix)
stopifnot(all(vapply(draws, function(x)
  all(c("a", "Topt", "cc") %in% colnames(x)), logical(1))))

chain_info <- data.frame(
  chain = seq_along(draws),
  mean_a = vapply(draws, function(x) mean(x[, "a"]), numeric(1)),
  mean_Topt = vapply(draws, function(x) mean(x[, "Topt"]), numeric(1)),
  source = "six-chain fit"
)
high <- which(chain_info$mean_a > 0)
# A NEW independent six-chain fit of the SAME unrestricted model. Starting
# all six near the negative mode is diagnostic, not a constraint on a < 0.
# Its distinct cache name prevents reuse of the older one-chain fit.
set.seed(20260924)
low_fit <- cached_fit("RH90_bayesTPC_LINEX_negative_start_6chains.rds", function()
  bayesTPC::b_TPC(data = dat, model = spec, nchains = 6L,
                  inits = list(a = -0.32, Topt = 22, cc = 3.5),
                  burn = 15000L, niter = 25000L))
low_chains <- find_chains(low_fit)
if (is.null(low_chains) || length(low_chains) != 6L)
  stop("Expected six independently sampled negative-start LINEX chains")
low_start_draws <- lapply(low_chains, as.matrix)
stopifnot(all(vapply(low_start_draws, function(x)
  all(c("a", "Topt", "cc") %in% colnames(x)), logical(1))))
draws <- c(draws, low_start_draws)
chain_info <- rbind(chain_info, data.frame(
  chain = 7:12,
  mean_a = vapply(low_start_draws, function(x) mean(x[, "a"]), numeric(1)),
  mean_Topt = vapply(low_start_draws, function(x) mean(x[, "Topt"]), numeric(1)),
  source = "six-chain negative-start fit"))
low <- which(chain_info$mean_a < 0 & chain_info$chain > 6L)
if (length(low) != 6L)
  warning("Some negative-start chains did not stay in the negative mode; ",
          "inspect all twelve traces and chain modes before interpretation")
write.csv(chain_info, file.path(root, "results", "RH90_bayesTPC_chain_modes.csv"),
          row.names = FALSE)

grid <- seq(min(dat$Temp), max(dat$Temp), length.out = 201)
thin_draws <- function(x, n = 1000L) {
  ix <- unique(round(seq(1, nrow(x), length.out = min(nrow(x), n))))
  x[ix, , drop = FALSE]
}
mode_draws <- function(indices, n = 1000L) {
  if (!length(indices)) return(NULL)
  x <- do.call(rbind, lapply(draws[indices], as.matrix))
  thin_draws(x, n)
}
high_draws <- mode_draws(high)
low_draws <- mode_draws(low)
curve_probs <- function(x) {
  p <- vapply(grid, function(t) {
    d <- x[, "a"] * (t - x[, "Topt"])
    plogis(-exp(pmin(d, 700)) + d + x[, "cc"])
  }, numeric(nrow(x)))
  if (is.null(dim(p))) p <- matrix(p, nrow = nrow(x))
  p
}
curve_summary <- function(p, center = c("median", "mean")) {
  center <- match.arg(center)
  rbind(apply(p, 2, quantile, probs = 0.025),
        if (center == "mean") colMeans(p) else apply(p, 2, median),
        apply(p, 2, quantile, probs = 0.975))
}
high_p <- if (!is.null(high_draws)) curve_probs(high_draws) else NULL
low_p <- if (!is.null(low_draws)) curve_probs(low_draws) else NULL
high_band <- if (!is.null(high_p)) curve_summary(high_p) else NULL
low_band <- if (!is.null(low_p)) curve_summary(low_p) else NULL
# Each mode supplies the same number of draws, irrespective of how many
# initial chains landed there. These are predictive scenario weights only.
if (!is.null(high_p) && !is.null(low_p)) {
  n_mix <- min(nrow(high_p), nrow(low_p))
  pick <- function(n) unique(round(seq(1, n, length.out = n_mix)))
  mixed_p <- rbind(high_p[pick(nrow(high_p)), , drop = FALSE],
                   low_p[pick(nrow(low_p)), , drop = FALSE])
  mixed_band <- curve_summary(mixed_p, "mean")
  write.csv(data.frame(Temp = grid, mean_prediction = mixed_band[2, ],
                       lower_95 = mixed_band[1, ], upper_95 = mixed_band[3, ]),
            file.path(root, "results", "RH90_bayesTPC_equal_mode_prediction.csv"),
            row.names = FALSE)
} else {
  mixed_band <- NULL
}
quadratic_chains <- find_chains(quadratic_fit)
if (is.null(quadratic_chains) || length(quadratic_chains) != 6L)
  stop("Expected six quadratic chains as a coda mcmc.list")
quad_draws <- do.call(rbind, lapply(quadratic_chains, as.matrix))
stopifnot(all(c("B0", "B1", "B2") %in% colnames(quad_draws)))
qix <- unique(round(seq(1, nrow(quad_draws),
                        length.out = min(nrow(quad_draws), 1500L))))
qp <- vapply(grid, function(t)
  plogis(quad_draws[qix, "B0"] + quad_draws[qix, "B1"] * t +
           quad_draws[qix, "B2"] * t^2), numeric(length(qix)))
if (is.null(dim(qp))) qp <- matrix(qp, nrow = length(qix))
quad_band <- apply(qp, 2, quantile, probs = c(0.025, 0.5, 0.975))
base_panel <- function(title) {
  plot(dat$Temp, dat$Trait / dat$n, xlim = range(grid), ylim = c(0, 1),
       xlab = "Temperature (C)", ylab = "Juvenile survival", main = title,
       pch = 19, col = "black")
}
draw_band <- function(b, title, col) {
  base_panel(title)
  if (is.null(b)) {
    text(mean(grid), 0.48, "Mode not sampled")
    return(invisible(NULL))
  }
  polygon(c(grid, rev(grid)), c(b[1, ], rev(b[3, ])),
          col = grDevices::adjustcolor(col, alpha.f = 0.18), border = NA)
  lines(grid, b[2, ], col = col, lwd = 2)
  points(dat$Temp, dat$Trait / dat$n, pch = 19)
}

# Cairo embeds fonts where available and avoids spaced/overprinted text in
# PDF viewers that substitute fonts for base R's unembedded Helvetica.
open_pdf <- function(path, width, height) {
  if (isTRUE(capabilities("cairo"))) {
    grDevices::cairo_pdf(path, width = width, height = height, family = "sans")
  } else {
    grDevices::pdf(path, width = width, height = height,
                   useDingbats = FALSE)
  }
}
open_pdf(file.path(output, "RH90_bayesTPC_four_panel_diagnostic.pdf"),
         width = 11, height = 8)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1), oma = c(2, 0, 0, 0))
draw_band(high_band, "(A) LINEX: high mode (six-chain fit)", "#1B9E77")
draw_band(low_band, "(B) LINEX: low mode (six low starts)", "#D95F02")
draw_band(mixed_band, "(C) LINEX: 50:50 mode scenario", "#4D4D4D")
base_panel("(D) Quadratic: six chains pooled")
polygon(c(grid, rev(grid)), c(quad_band[1, ], rev(quad_band[3, ])),
        col = grDevices::adjustcolor("#7570B3", alpha.f = 0.18), border = NA)
lines(grid, quad_band[2, ], col = "#7570B3", lwd = 2)
points(dat$Temp, dat$Trait / dat$n, pch = 19)
mtext("50:50 weighting is illustrative; these chains do not estimate mode probabilities or stacking weights.",
      side = 1, outer = TRUE, line = 0.5, cex = 0.75)
dev.off()

cols <- grDevices::hcl.colors(6L, "Dark 3")
linex_trace_value <- function(x, name) {
  if (name != "p30") return(x[, name])
  d <- x[, "a"] * (30 - x[, "Topt"])
  plogis(-exp(pmin(d, 700)) + d + x[, "cc"])
}
quad_trace_value <- function(x, name) {
  if (name != "p30") return(x[, name])
  plogis(x[, "B0"] + 30 * x[, "B1"] + 30^2 * x[, "B2"])
}
trace_panel <- function(selected_chains, comparison_chains, name, title,
                        value_fun, labels) {
  selected <- lapply(selected_chains, value_fun, name = name)
  # Both LINEX files share the same parameter scales for direct comparison.
  vals <- unlist(lapply(comparison_chains, value_fun, name = name),
                 use.names = FALSE)
  plot(NA, xlim = c(1, max(vapply(selected, length, integer(1)))),
       ylim = extendrange(range(vals), f = 0.08), xlab = "Stored iteration",
       ylab = title, main = title)
  if (name == "a") abline(h = 0, col = "#777777", lty = 3)
  for (j in seq_len(6L)) {
    y <- selected[[j]]
    ix <- unique(round(seq(1, length(y), length.out = min(length(y), 1500L))))
    lines(ix, y[ix], col = grDevices::adjustcolor(cols[j], alpha.f = 0.75),
          lwd = 0.65)
  }
  if (name == "a") legend("topright", inset = 0.01,
    legend = paste("Chain", labels), col = cols,
    lty = 1,
    lwd = 1, cex = 0.65, bg = "white")
  if (name == "B0") legend("topright", inset = 0.01,
    legend = paste("Chain", labels), col = cols,
    lty = 1, lwd = 1, cex = 0.65, bg = "white")
}
write_trace <- function(filename, selected_chains, comparison_chains,
                        variables, titles, value_fun, label, labels) {
  stopifnot(length(selected_chains) == 6L,
            all(vapply(selected_chains, function(x)
              all(variables[variables != "p30"] %in% colnames(x)), logical(1))))
  open_pdf(file.path(output, filename), width = 10, height = 8)
  par(mfrow = c(2, 2), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))
  for (i in seq_along(variables))
    trace_panel(selected_chains, comparison_chains,
                variables[i], titles[i], value_fun, labels)
  mtext(label, outer = TRUE, cex = 0.85)
  dev.off()
}
linex_titles <- c("LINEX a", "Optimal temperature", "LINEX cc",
                  "Predicted survival at 30 C")
linex_vars <- c("a", "Topt", "cc", "p30")
write_trace("RH90_bayesTPC_LINEX_positive_six_chain_trace.pdf",
            draws[1:6], draws, linex_vars, linex_titles, linex_trace_value,
            "Unrestricted LINEX: six original chains (default starts)", 1:6)
write_trace("RH90_bayesTPC_LINEX_negative_six_chain_trace.pdf",
            draws[7:12], draws, linex_vars, linex_titles, linex_trace_value,
            "Unrestricted LINEX: six separate negative-start chains", 1:6)
quad_chain_draws <- lapply(quadratic_chains, as.matrix)
write_trace("RH90_bayesTPC_quadratic_six_chain_trace.pdf",
            quad_chain_draws, quad_chain_draws,
            c("B0", "B1", "B2", "p30"),
            c("Quadratic intercept B0", "Quadratic slope B1",
              "Quadratic curvature B2", "Predicted survival at 30 C"),
            quad_trace_value, "Quadratic: six chains", 1:6)
# Avoid leaving an older, ambiguous two-page LINEX trace next to the three
# newly named one-page files.
unlink(file.path(output, "RH90_bayesTPC_LINEX_six_chain_trace.pdf"))

message("RH90 bayesTPC diagnostic plots written to ", output)
