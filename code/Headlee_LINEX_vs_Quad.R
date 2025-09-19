# Aedes Traits exploration: Model fitting and plotting
# Author: Paul Huxley 
# Date: 2025-09-19

rm(list = ls())
graphics.off()

# Load libraries
library(nimble)
library(bayesTPC)
library(tidyverse)
library(purrr)
library(patchwork)
library(cowplot)
library(EnvStats)
library(parallel)
library(hdi)
library(HDInterval)
library(grid)
library(gridExtra) 
library(vroom)

# Load and process data 
headlee <- read_csv("data/headlee_data.csv")

# Quick checks
head(headlee)
unique(headlee$doi)    # should be one dataset (or you want the whole file as single dataset)
unique(headlee$originaltraitname)

#----------------------------------------------------------
# Specify LINEX model exactly as you did (unchanged)
#----------------------------------------------------------
reset_models()
my_name<- "LinExp_binom"
LinExp_formula <- expression(-exp( -a*(Topt-Temp) ) + a*(Temp-Topt) + cc)
LinExp_priors <- c(
  cc = "dexp(0.1)",
  a = "dnorm(0, 1)",
  Topt  = "dnorm(25, 1/1000)"
)
LinExp_binom <- specify_binomial_model(
  my_name,
  parameters = LinExp_priors,
  formula = LinExp_formula
)

#----------------------------------------------------------
# Prepare data for a single dataset fit (avoid 1:N warning)
#----------------------------------------------------------
# Treat the entire headlee dataframe as one dataset (no per-doi loop)
success <- headlee$interactor1number       # number of successes per observation (vector)
trials_per_obs <- headlee$N                # trials per observation (vector)
Temp <- headlee$interactor1temp
Nobs <- length(success)                    # scalar: number of observations

# Build the data list including both the scalar Nobs (for model loops) and
# the per-observation trials vector n (for dbin calls)
fit_dat <- list(Trait = success, Temp = Temp, N = Nobs, n = trials_per_obs)

#----------------------------------------------------------
# Fit LINEX (binomial) and binomial quadratic models
#----------------------------------------------------------
set.seed(2025)

linexp_fit <- b_TPC(
  data = fit_dat,
  model = LinExp_binom,
  inits = list(a = 0.01, Topt = 28, cc = 0.1),
  burn = 15000,
  niter = 25000
)

quad_fit <- bayesTPC::b_TPC(
  data = fit_dat,
  model = "binomial_glm_quad",
  samplerType = "AF_slice",
  priors = list(B0 = "dunif(-100,100)", 
                B1 = "dunif(-10,10)", 
                B2 = "dunif(-1,1)"),
  inits = list(B0 = 10, B1 = 1.2, B2 = -0.02),
  burn = 5000, niter = 15000
)

#----------------------------------------------------------
# Plot fits side-by-side
#----------------------------------------------------------
# Single PNG with two panels: LINEX | Quadratic
png("figures/Headlee_fits.png", width = 800, height = 400)
par(mfrow = c(1,2), mar = c(4,4,3,1))
plot(linexp_fit, ylim = c(0,1), legend = FALSE, main = "LINEX")
plot(quad_fit, ylim = c(0,1), legend = FALSE, main = "Quadratic")
dev.off()

# Also show on screen/device
par(mfrow = c(1,2))
plot(linexp_fit, ylim = c(0,1), legend = FALSE, main = "LINEX")
plot(quad_fit, ylim = c(0,1), legend = FALSE, main = "Quadratic")

#----------------------------------------------------------
# Compute and save WAIC
#----------------------------------------------------------
WAIC_df <- data.frame(
  Model = c("LINEX", "QUAD"),
  WAIC = c(get_WAIC(linexp_fit)[1], get_WAIC(quad_fit)[1])
)

print(WAIC_df)

png("figures/headlee_WAIC_table.png", width = 400, height = 200)  
grid.table(WAIC_df)
dev.off()

# End of script