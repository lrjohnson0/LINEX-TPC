#--------------------------
# Doeurk et al Aedes eggs: Model fitting and plotting
# Author: Paul Huxley
# Date: September 2025
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
df <- read_csv("data/Doeurk_et_al_2025.csv")

df_surv <- df %>%
  group_by(Species, Temperature, Replicate) %>%
  filter(Hatch.rate != 0) %>%
  mutate(Trait =  Hatching / eggs.hatching, 
         Temp = Temperature
         ) %>%
  ungroup()

# Filter out some that don't seem suitable for fitting (NOT CONCLUSIVE!)
#model_dat <- df_max %>%
#  filter(!citation_short %in% c("Reisen_1993", "TurellandLundstrom_1990", "Cornel_1993"))

# Plot 
ggplot(data = df_surv) +
  geom_jitter(aes(Temp, Trait, group = Temp, colour = Species),
              width = 0.2, alpha=0.75) +
  #facet_wrap(~citation_short, ncol = 3, scales = "free_y") +
  theme(legend.position = "none")



#----------------------------------------------------------
# Specify LINEX model
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

species <- unique(df_surv$Species)
n_species <- length(species)

linexp_fits <- vector("list", n_species)
quad_fits <- vector("list", n_species)
WAIC_df <- data.frame(species = species, LINEX = NA, QUAD = NA)

par(mfrow = c(2, 2), bty = "n")

for (i in seq_along(species)) {
  #subdat <- model_dat %>% filter(citation_short == citations[i])
  # Prepare binomial data
  success <- df_surv$Hatching
  N <- df_surv$eggs.hatching
  Temp <- df_surv$Temp
  
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
       main = paste("LINEX:", species[i]))
  plot(quad_fits[[i]], ylim = c(0,1), legend = FALSE,
       main = paste("Quadratic:", species[i]))
  
  # Compute WAIC
  WAIC_df[i, "LINEX"] <- get_WAIC(linexp_fits[[i]])[1]
  WAIC_df[i, "QUAD"] <- get_WAIC(quad_fits[[i]])[1]
}

print(WAIC_df)

library(gridExtra)

# Save WAIC table as PNG
png("figures/Doeurk_WAIC_table.png", width = 750, height = 400)  
grid.table(WAIC_df)
dev.off()

# Save all fits to one PNG (2 rows, 7 columns)
png("figures/Doeurk_egg_fits.png", width = 2800, height = 800)  

par(mfrow = c(2, 7), mar = c(5,5,4,2) + 0.1) # 2 rows, 7 columns

# Top row: LINEX fits
for (i in seq_along(species)) {
  plot(linexp_fits[[i]], ylim = c(0,1), legend = FALSE,
       main = paste("LINEX:", species[i]))
}

# Bottom row: Quadratic fits
for (i in seq_along(species)) {
  plot(quad_fits[[i]], ylim = c(0,1), legend = FALSE,
       main = paste("Quadratic:", species [i]))
}

dev.off()
