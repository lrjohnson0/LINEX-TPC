#--------------------------
# Aedes Trsits exploration: Model fitting and plotting
# Author: Paul Huxley
# Date: Sept 2025
#--------------------------
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

unique(headlee$N)
unique(headlee$notes)
unique(headlee$originaltraitname)
unique(headlee$interactor1number)
unique(headlee$citation)

ggplot(data = headlee) +
  geom_jitter(aes(interactor1temp, originaltraitvalue, group = doi),
              width = 0.2, alpha=0.75) +
  facet_wrap(~doi, ncol = 5, scales = "free_y") +
  theme(legend.position = "none")

#----------------------------------------------------------
# Specify LINEX model exactly as you did
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
# Fit LINEX and binomial quadratic models for each citation
#----------------------------------------------------------

dois <- unique(headlee$doi)
n_dois <- length(dois)

linexp_fits <- vector("list", n_dois)
quad_fits <- vector("list", n_dois)
WAIC_df <- data.frame(doi = dois, LINEX = NA, QUAD = NA)

par(mfrow = c(2, 2), bty = "n")

for (i in seq_along(dois)) {
  subdat <- headlee %>% filter(doi == dois[i])
  # Prepare binomial data
  success <- subdat$interactor1number
  N <- subdat$N
  Temp <- subdat$interactor1temp
  
  fit_dat <- list(Trait = success, Temp = Temp, n = N)
  
  # Fit LINEX binomial using your specified model
  linexp_fits[[i]] <- b_TPC(
    data = fit_dat,
    model = LinExp_binom,
    inits = list(a = 0.01, Topt = 28, cc = 0.1),
    burn = 15000,
    niter = 25000
  )
  
  # Fit binomial quadratic GLM
  quad_fits[[i]] <- bayesTPC::b_TPC(
    data = fit_dat,
    "binomial_glm_quad",
    samplerType = "AF_slice",
    priors = list(B0 = "dunif(-100,100)", 
                  B1 = "dunif(-10,10)", 
                  B2 = "dunif(-1,1)"),
    inits = list(B0 = 10, B1 = 1.2, 
                 B2 = -0.02),
    burn = 5000, niter = 15000
  )
  
  # Plot fits
  par(mfrow = c(1,2))
  plot(linexp_fits[[i]], ylim = c(0,1), legend = FALSE,
       main = paste("LINEX:", dois[i]))
  plot(quad_fits[[i]], ylim = c(0,1), legend = FALSE,
       main = paste("Quadratic:", dois[i]))
  
  # Compute WAIC
  WAIC_df[i, "LINEX"] <- get_WAIC(linexp_fits[[i]])[1]
  WAIC_df[i, "QUAD"] <- get_WAIC(quad_fits[[i]])[1]
}

print(WAIC_df)

library(gridExtra)

# Save WAIC table as PNG
png("figures/headlee_WAIC_table.png", width = 750, height = 400)  
grid.table(WAIC_df)
dev.off()

# Save all fits to one PNG (1 rows, 2 columns)
png("figures/headlee.png", width = 800, height = 400)  

par(mfrow = c(1,2)) # 2 rows, 7 columns

# Top row: LINEX fits
for (i in seq_along(dois)) {
  plot(linexp_fits[[i]], ylim = c(0,1), legend = FALSE,
       main = paste("LINEX:", dois[i]))
}

# Bottom row: Quadratic fits
for (i in seq_along(dois)) {
  plot(quad_fits[[i]], ylim = c(0,1), legend = FALSE,
       main = paste("Quadratic:", dois[i]))
}

dev.off()

