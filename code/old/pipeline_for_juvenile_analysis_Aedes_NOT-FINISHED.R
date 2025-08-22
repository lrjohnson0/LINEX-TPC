#--------------------------
# Project: Evaluation of Linex vs Binomial quadratic glm using 
# Aedes data 
# Author: Paul Huxley
# Date: August 2025
#--------------------------
rm(list = ls())
graphics.off()

# Load libraries and data
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

# Trait data
df <- readRDS("~/Desktop/LINEX-TPC/data/zj-backup-of-all-pseudo-distributions-inc-global.rds")
ggplot() +
  geom_point(aes(Temp, Trait /N), df, shape=21, fill="transparent", col="blue", size=1) +
  facet_wrap(~curve_id, ncol = 5, scales = "free_y") +
  #scale_y_continuous(expression(plain(paste("trait value")))) +
  theme_bw() +
  theme(text = element_text(size = 10, family = "Palatino")) +
  labs(x = expression(plain(paste(" Temperature, ", degree, "C")))) +
  theme(legend.position = "none")


juvenile_survival_data <- df

#_______________________________________________________________________________
#                               TPC fitting 
#_______________________________________________________________________________
data <- split(juvenile_survival_data, juvenile_survival_data$curve_id) 
temps <- seq(0, 45, length = 2000)
library(foreach)
library(doParallel)

# Parallel model fitting for quadratic and LINEX models, with egg constraint at 45°C

# Setup parallel backend
cores <- parallel::detectCores()
cl <- parallel::makeCluster(min(10, cores))
doParallel::registerDoParallel(cl)

iterations <- 15000

# --- Quadratic GLM fit ---
output_models <- foreach(i=1:length(data)) %dopar% {
  library(nimble)
  library(bayesTPC)
  dat <- data[[i]]
  data.list <- list(
    Temp = dat$Temp,
    Trait = dat$Trait,
    n = dat$N
  )
  
  
  quad_priors <- 
    # Strong negative prior on B2 ensures curve drops by 45°C
    list(B0 = "dunif(-100,45)",
         B1 = "dunif(-10,10)",
         B2 = "dunif(-1, -0.001)")
  
  
  # Fit model
  fit <- bayesTPC::b_TPC(
    data = data.list,
    model = "binomial_glm_quad",
    samplerType = "AF_slice",
    priors = quad_priors,
    inits = list(B0 = 20, B1 = -5, B2 = -0.02),
    burn = 7500,
    niter = iterations
  )
  gc()
  list(
    samp = fit$samples,
    MAP = fit$MAP_parameters,
    summary = predict(fit, temp_interval = temps),
    waic = bayesTPC::get_WAIC(fit),
    data = dat
  )
}

# --- LINEX fit ---
linex_models <- foreach(i=1:length(data)) %dopar% {
  library(nimble)
  library(bayesTPC)
  reset_models()
  dat <- data[[i]]
  data.list <- list(
    Temp = dat$Temp,
    Trait = dat$Trait,
    n = dat$N
  )

  # Define formula for LINEX model
  LinExp_formula <- expression(-exp(-a * (Topt - Temp)) + a * (Temp - Topt) + cc)
  
  LinExp_priors <- 
    # Constrain Topt below 45 and use steep a prior 
    c(
      cc = "dexp(0.1)",
      a = "dunif(0.3, 2)",
      Topt = "dunif(10, 45)"
    )
  
  my_name <- "LinExp_binom"
  LinExp_binom <- specify_binomial_model(
    my_name,
    parameters = LinExp_priors,
    formula = LinExp_formula
  )
  inits <- list(a = 0.5, Topt = median(dat$Temp), cc = 1)
  
  fit <- bayesTPC::b_TPC(
    data = data.list,
    model = "LinExp_binom",
    inits = inits,
    samplerType = "AF_slice",
    burn = 7500,
    niter = iterations
  )
  gc()
  list(
    fit = fit,
    samp = fit$samples,
    MAP = fit$MAP_parameters,
    summary = predict(fit, temp_interval = temps),
    waic = bayesTPC::get_WAIC(fit),
    data = dat
  )
}

parallel::stopCluster(cl)

## Check convergence of all models
#par(mfrow = c(5,3), mar = c(4,3,3,1) +.1)
#for (fit in output_models) coda::traceplot(fit$samp)
#for (fit in linex_models) coda::traceplot(fit$samp)

#_______________________________________________________________________________
# Extract/estimate TPC params
pars <- list()
fits <- list()
param_hpds <- list()
zj.post <- list()  
zj.HPD <- list()    

for (i in 1:length(output_models)) {
  print(i)
  
  library(coda)
  library(HDInterval)
  
  samples  <- output_models[[i]]$samp
  map      <- output_models[[i]]$MAP
  fit_sum  <- output_models[[i]]$summary
  hpds     <- HPDinterval(samples)
  
  mybinomial <- get_model_function("binomial_glm_quad")
  l <- dim(output_models[[1]]$samp)[1]
  
  mydat1 <- matrix(NA, nrow = l, ncol = length(temps))
  
  for (j in 1:l) {
    mydat1[j,] <-  mybinomial(Temp = temps,
                              B0 = samples[j, 1],
                              B1 = samples[j, 2],
                              B2 = samples[j, 3])
  }
  
  ss <- l
  ts <- length(temps)
  
  zj.post[[i]] <- mydat1  # Assign predictions to zj.post[[i]]
  
  zj.HPD[[i]] <- matrix(NA, nrow = 2, ncol = ts)  # HPD interval matrix
  
  for (j in 1:ts) {
    mytemp <- mcmc(zj.post[[i]][, j])
    hpd_values <- coda::HPDinterval(mytemp)
    zj.HPD[[i]][, j] <- c(hpd_values[1], hpd_values[2])
  }
  
  names(zj.HPD)[i] <- paste0("curve_", i)
  
  get_Tpk_indices <- apply(mydat1, MARGIN = 1, FUN = which.max)
  get_Bpks        <- apply(mydat1, MARGIN = 1, FUN = max)
  get_Tpk_values  <- temps[get_Tpk_indices]
  
  Tpk_hpd <- HPDinterval(as.mcmc(get_Tpk_values), prob = .95)
  Bpk_hpd <- HPDinterval(as.mcmc(get_Bpks), prob = .95)
  
  survival_tolerance <- 0.01 #Sigmoids asymptote at 1 and 0, so we need to set a cut off for min/max temps
  low_cutoffs <- apply(mydat1, MARGIN = 1, FUN = function(x) min(which(x >= survival_tolerance)))
  high_cutoffs <- apply(mydat1, MARGIN = 1, FUN = function(x) max(which(x >= survival_tolerance)))
  
  get_min_temps <- temps[low_cutoffs]
  get_max_temps <- temps[high_cutoffs]
  
  Tmin_hpd <- HPDinterval(as.mcmc(get_min_temps), prob = .95)
  Tmax_hpd <- HPDinterval(as.mcmc(get_max_temps), prob = .95)
  
  pars[[i]] <- list(
    Tmin_mean = mean(get_min_temps),
    Tmin_median = median(get_min_temps),
    Tmin_lwr = Tmin_hpd[1],
    Tmin_upr = Tmin_hpd[2],
    Tmax_mean = mean(get_max_temps),
    Tmax_median = median(get_max_temps),
    Tmax_lwr = Tmax_hpd[1],
    Tmax_upr = Tmax_hpd[2],
    Tpk_mean = mean(get_Tpk_values),
    Tpk_median = median(get_Tpk_values),
    Tpk_lwr = Tpk_hpd[1],
    Tpk_upr = Tpk_hpd[2],
    Bpk_mean = mean(get_Bpks),
    Bpk_median = median(get_Bpks),
    Bpk_lwr = Bpk_hpd[1],
    Bpk_upr = Bpk_hpd[2],
    traitname = "survival", 
    curve_id = unique(data[[i]]$curve_id)
  )
  
  fit_sum$temp_interval <- NULL
  
  fits[[i]] <- bind_cols(Temp = temps,fit_sum,
                         curve_id = unique(data[[i]]$curve_id)
  )
  
  param_hpds[[i]] <- data.frame(
    Tmin_lwr = Tmin_hpd[1],
    Tmin_upr = Tmin_hpd[2],
    Tmax_lwr = Tmax_hpd[1],
    Tmax_upr = Tmax_hpd[2],
    traitname = "survival", 
    curve_id = unique(data[[i]]$curve_id)
  )
}

param_hpds_combined <- do.call(rbind, param_hpds)
param_hpds_combined

all_pars <- bind_rows(pars)
all_pars

zjHPD_df <- lapply(1:length(zj.HPD), function(i) {
  data.frame(
    Temp = temps,
    Trait_lwr = zj.HPD[[i]][1,],
    Trait_upr = zj.HPD[[i]][2,],
    curve_id = unique(data[[i]]$curve_id), length.out = length(temps))
  
})

zjHPD_combined <- do.call(rbind, zjHPD_df)

all_fits  <- bind_rows(fits) %>%
  dplyr::rename(Trait = medians, Trait_lwr = lower_bounds, Trait_upr = upper_bounds) %>% 
  dplyr::select(-temp_interval) 

linex_all_fits  <- bind_rows(lapply(linex_models, function(x) {
  fit_sum <- x$summary
  fit_sum$temp_interval <- NULL
  curve_id_rep <- rep(unique(x$data$curve_id), length.out = length(temps))
  bind_cols(Temp = temps, fit_sum,
            curve_id = curve_id_rep)
})) %>%
  dplyr::rename(Trait = medians, Trait_lwr = lower_bounds, Trait_upr = upper_bounds) %>%
  dplyr::select(-temp_interval)

# For plotting, update plotdata/plotfits with all six curves/labels

label_levels <- c(
  "90", "75", "60", "45", "30"
)
my_cols <- c(
  "90" = "#0078B2",
  "75" = "#E2A233",
  "60" = "#4BACC6",
  "45" = "#AE5650",
  "30" = "#7FBFF5"
)

plotdata <- juvenile_survival_data %>%
  mutate(
    curve_id = factor(curve_id,  levels = label_levels))


plotfits <- all_fits %>%
  mutate(
    curve_id = factor(curve_id, levels = levels(plotdata$curve_id)) 
  )

linex_plotfits <- linex_all_fits %>%
  mutate(
    curve_id = factor(curve_id, levels = levels(plotdata$curve_id)) 
)


plotdata %>% distinct(curve_id)
plotfits %>% distinct(curve_id)

linex_plotfits %>% distinct(curve_id)

# --- Quadratic fit plot ---
output_graph <- 
  plotdata %>% 
  ggplot() +
  geom_ribbon(aes(Temp, ymin = Trait_lwr, ymax = Trait_upr, 
                  group = curve_id, 
                  fill = curve_id,
                  col = curve_id), plotfits,
              alpha = 0.5,
              inherit.aes = T) +
  geom_line(aes(Temp, Trait, group = curve_id, col = curve_id), plotfits) + 
  geom_point(aes(Temp, 
                 Trait / n, 
                 group = curve_id, 
                 shape = curve_id,
                 fill =  curve_id),
             col = "white", position = position_jitter(width = .25, height = 0),
             size = 2, 
             stroke=0.1) +
  scale_fill_manual(name=expression(bold("")),
                    values = my_cols,
                    guide = guide_legend(nrow=5,ncol=1,
                                         direction = "vertical",
                                         title.position = "top",
                                         title.hjust=0.5))+
  scale_colour_manual(name=expression(bold("")),
                      values = my_cols,
                      guide = guide_legend(nrow=5,ncol=1,
                                           direction = "vertical",
                                           title.position = "top",
                                           title.hjust=0.5))+
  scale_shape_manual(
    name = expression(bold("")),
    values = rep(23, length(label_levels)),
    guide = guide_legend(
      nrow = 5, ncol = 1, direction = "vertical",
      title.position = "top", title.hjust = 0.5
    )
  ) +
  scale_y_continuous(expression("survival probability"),
                     limits =c(-0.01,1.01),
                     expand = c(0, 0),
                     breaks=seq(0,1, by=0.25)) +
  theme_bw() +
  theme(text = element_text(size = 12, family = "Palatino")) +
  labs(x = expression(plain(paste(" Temperature, ", degree, "C")))) +
  theme(legend.position = c(.15, .65),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_rect(colour = "white")
  ) +
  ggtitle("(a) Quadratic")

output_graph 

# --- LINEX fit plot ---
linex_output_graph <- 
  plotdata %>%
  ggplot() +
  geom_ribbon(aes(Temp, ymin = Trait_lwr, ymax = Trait_upr, 
                  group = curve_id, 
                  fill = curve_id,
                  col = curve_id), linex_plotfits,
              alpha = 0.5,
              inherit.aes = TRUE) +
  geom_line(aes(Temp, Trait, group = curve_id, col = curve_id), linex_plotfits) +
  geom_point(aes(Temp, 
                 Trait / n, 
                 group = curve_id, 
                 shape = curve_id,
                 fill = curve_id),
             col = "white", position = position_jitter(width = .25, height = 0),
             size = 2,
             stroke = 0.1) +
  scale_fill_manual(name = expression(bold("")),
                    values = my_cols,
                    guide = guide_legend(nrow = 5, ncol = 1,
                                         direction = "vertical",
                                         title.position = "top",
                                         title.hjust = 0.5)) +
  scale_colour_manual(name = expression(bold("")),
                      values = my_cols,
                      guide = guide_legend(nrow = 5, ncol = 1,
                                           direction = "vertical",
                                           title.position = "top",
                                           title.hjust = 0.5)) +
  scale_shape_manual(
    name = expression(bold("")),
    values = rep(23, length(label_levels)),
    guide = guide_legend(
      nrow = 5, ncol = 1, direction = "vertical",
      title.position = "top", title.hjust = 0.5
    )
  ) +
  scale_y_continuous(expression("survival probability"),
                     limits =c(-0.01,1.01),
                     expand = c(0, 0),
                     breaks=seq(0,1, by=0.25)) +
  theme_bw() +
  theme(text = element_text(size = 12, family = "Palatino")) +
  labs(x = expression(plain(paste(" Temperature, ", degree, "C")))) +
  theme(legend.position = c(.15, .65),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_rect(colour = "white")
  ) +
  ggtitle("(b) LINEX")

linex_output_graph 

output_graph + linex_output_graph

# ---------------- WAIC extraction and summary -----------------

extract_waic <- function(x) {
  if (!is.null(x$waic) && is.numeric(x$waic) && !is.null(names(x$waic)) && "WAIC" %in% names(x$waic)) {
    return(x$waic["WAIC"])
  }
  return(NA_real_)
}

quad_WAIC <- map2_dfr(output_models, levels(plotfits$curve_id), ~tibble(
  curve_id = .y,
  waic_quadratic = extract_waic(.x)
))

linex_WAIC <- map2_dfr(linex_models, levels(linex_plotfits$curve_id), ~tibble(
  curve_id = .y,
  waic_linex = extract_waic(.x)
))

all_WAIC <- full_join(quad_WAIC, linex_WAIC, by = "curve_id") %>%
  mutate(
    waic_quadratic = sprintf("%.1f", waic_quadratic),
    waic_linex = sprintf("%.1f", waic_linex),
    curve_id = factor(curve_id, levels = levels(plotfits$curve_id))
  ) %>%
  arrange(curve_id)

print(all_WAIC)

waic_sum_quadratic <- sum(as.numeric(all_WAIC$waic_quadratic), na.rm = TRUE)
waic_sum_linex     <- sum(as.numeric(all_WAIC$waic_linex), na.rm = TRUE)

cat("Quadratic WAIC sum:", waic_sum_quadratic, "\n")
cat("LINEX WAIC sum:", waic_sum_linex, "\n")



# ---------- WAIC table as 3rd panel and save to PDF ----------

library(grid)
library(gridExtra)
library(cowplot)
library(patchwork)

# Prepare a printable WAIC table (include sums as a final row)
waic_tbl_print <- all_WAIC %>%
  mutate(
    across(c(waic_quadratic, waic_linex), ~ ifelse(is.na(.), "", .))
  ) %>%
  rename(
    curve_id = curve_id,
    `WAIC (Quadratic)` = waic_quadratic,
    `WAIC (LINEX)` = waic_linex
  ) %>%
  bind_rows(tibble(
    curve_id = "Sum",
    `WAIC (Quadratic)` = sprintf("%.1f", waic_sum_quadratic),
    `WAIC (LINEX)` = sprintf("%.1f", waic_sum_linex)
  ))

# Build a table grob
tab_theme <- ttheme_minimal(
  base_size = 11,
  core = list(fg_params = list(fontfamily = "Palatino", hjust = 0.5, x = 0.5)),
  colhead = list(fg_params = list(fontfamily = "Palatino", fontface = "bold"))
)

tab_grob <- tableGrob(waic_tbl_print, rows = NULL, theme = tab_theme)

# Optionally tweak column widths to keep the table narrower
# Adjust "in" sizes as you prefer
tab_grob$widths <- unit(c(0.7, 1.3, 1.3), "in")

# Wrap table grob in a ggplot-compatible object and add a panel label
waic_table_plot <- ggdraw() +
  draw_grob(tab_grob, x = 0, y = 0, width = 1, height = 1)  # +
  #draw_label("(c) WAIC", x = 0.02, y = 0.98, hjust = 0, vjust = 1,
  #           fontface = "bold", 
  #           size = 12)

# Combine: (a) Quadratic, (b) LINEX, (c) WAIC table
combined_plot <- output_graph + linex_output_graph + waic_table_plot +
  plot_layout(nrow = 1, ncol = 3, widths = c(1, 1, 1.2))

# Save as PDF (use cairo_pdf for better font support)
#out_pdf <- "~/Desktop/LINEX-TPC/figures/survival_quad_linex_waic.pdf"
#ggsave(
#  filename = out_pdf,
#  plot = combined_plot,
#  width = 12, height = 3.75, units = "in",
#  device = cairo_pdf
#)

#message("Saved 1x3 plot with WAIC table to: ", out_pdf)

