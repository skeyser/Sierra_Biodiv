## -------------------------------------------------------------
##
## Script name: GDM Ferrier Z Matrix for the Fire Only Sites
##
## Script purpose: GDM for Sierra Bioacoustics Latent Z using
## traditional Ferrier GDM approach
##
## Author: Spencer R Keyser
##
## Date Created: 2024-10-22
##
## Email: srk252@cornell.edu
##
## Github: https://github.com/skeyser
##
## -------------------------------------------------------------
##
## Notes:
##
##
## -------------------------------------------------------------

## Defaults
options(scipen = 6, digits = 4)

## -------------------------------------------------------------

## Package Loading
library(dplyr)
library(ggplot2)
library(here)
library(purrr)

## GDM
library(gdm)
library(betapart)

## -------------------------------------------------------------

## Load in the Z-matrices
load(here("./Data/JAGS_Data/Occ2GDM_Data_95thresh_Fires.Rdata"))

## Reformat data
Z <- OccData$Z.posterior
dimnames(Z) <- list(NULL, OccData$ZcolNames, NULL)

## Read in the meta data
# aru_meta <- read.csv(here("./Data/ARU_Meta_120_FilteredOcc.csv")) 
# aru_sub <- aru_meta |> 
#   select(Cell_Unit, 
#          utme, utmn,
#          topo_elev,
#          standage_f3_mn,
#          cpycovr_f3_mn,
#          contains("high_prop"),
#          contains("lowmod_prop")
#   ) |> 
#   mutate(fire1_5yr_high_prop = fire1yr_high_prop + fire2_5yr_high_prop,
#          fire1_5yr_lowmod_prop = fire1yr_high_prop + fire2_5yr_lowmod_prop) |> 
#   select(Cell_Unit:cpycovr_f3_mn,
#          fire1_5yr_high_prop,
#          fire6_10yr_high_prop,
#          fire11_35yr_high_prop,
#          fire1_5yr_lowmod_prop,
#          fire6_10yr_lowmod_prop,
#          fire11_35yr_lowmod_prop)

aru_meta <- OccData$SiteMeta
aru_meta <- aru_meta |> select(-X, -Y)

## Combine fire classes
aru_meta <- aru_meta |> 
  mutate(HSF = fire1_5yr_high_prop + fire6_10yr_high_prop + fire11_35yr_high_prop,
         LMSF = fire1_5yr_lowmod_prop + fire6_10yr_lowmod_prop + fire11_35yr_lowmod_prop) |> 
  select(Cell_Unit:cancov, HSF, LMSF)


## -------------------------------------------------------------
##
## Begin Section: Format for GDM
##
## -------------------------------------------------------------

## Species data
mr <- vector(mode = "list", length = 100)
for(i in 1:dim(Z)[3]){
  Z.tmp <- Z[,,i]
  missRows <- which(rowSums(Z.tmp) == 0)
  mr[[i]] <- missRows
  #Z.tmp <- as.data.frame(Z.tmp)
  #Z.tmp <- Z.tmp[-missRows, ]
}

missRows <- unique(unlist(mr))
Z <- Z[-missRows,,]

## Environmental data columns
envTab <- aru_meta |>
  as.data.frame()
envTab <- envTab[-missRows, ]

## Format the data for GDM
npost <- dim(Z)[3]

gdm.fit.list <- replicate(3, vector(mode = "list", length = npost), simplify = FALSE)
gdm.fit.list <- tibble(
  beta_metric = c("Brep", "Brich", "Btotal"),
  models = vector("list", 3)) |>
  mutate(models = map(1:n(), ~vector("list", npost)))


# Create a progress bar
pb <- txtProgressBar(min = 0, max = npost, style = 3)

for(i in 1:npost){
  
  Z.samp <- as.data.frame(Z[,,i])
  #Z.samp$Cell_Unit <- envTab$Cell_Unit
  rownames(Z.samp) <- envTab$Cell_Unit
  envTab$Cell_Unit <- 1:length(unique(envTab$Cell_Unit))
  
  ## Partition the beta diversity based on Jaccard-based Podani metrics
  bpart <- adespatial::beta.div.comp(Z.samp, coef = "J")
  
  ## Replacement
  brep <- as.matrix(bpart$repl)
  dimnames(brep) <- list(envTab$Cell_Unit, envTab$Cell_Unit)
  brep <- cbind(Cell_Unit = envTab$Cell_Unit, brep)
  
  ## Richness
  brich <- as.matrix(bpart$rich)
  dimnames(brich) <- list(envTab$Cell_Unit, envTab$Cell_Unit)
  brich <- cbind(Cell_Unit = envTab$Cell_Unit, brich)
  
  ## Total Beta
  btotal <- as.matrix(bpart$D)
  dimnames(btotal) <- list(envTab$Cell_Unit, envTab$Cell_Unit)
  btotal <- cbind(Cell_Unit = envTab$Cell_Unit, btotal)
  
  # gdmTab <- formatsitepair(bioData = Z.samp,
  #                          bioFormat = 1,
  #                          dist = "jaccard",
  #                          XColumn = "utme",
  #                          YColumn = "utmn",
  #                          siteColumn = "Cell_Unit",
  #                          predData = envTab)
  gdmTab.brep <- formatsitepair(bioData = brep,
                           bioFormat = 3,
                           XColumn = "utme",
                           YColumn = "utmn",
                           siteColumn = "Cell_Unit",
                           predData = envTab)
  
  gdmTab.brich <- formatsitepair(bioData = brich,
                                bioFormat = 3,
                                XColumn = "utme",
                                YColumn = "utmn",
                                siteColumn = "Cell_Unit",
                                predData = envTab)
  
  gdmTab.btotal <- formatsitepair(bioData = btotal,
                                bioFormat = 3,
                                XColumn = "utme",
                                YColumn = "utmn",
                                siteColumn = "Cell_Unit",
                                predData = envTab)

  
  ## Fit GDM for one Z matrix slice
  gdm.fit.brep <- gdm(gdmTab.brep, geo = TRUE)
  gdm.fit.brich <- gdm(gdmTab.brich, geo = TRUE)
  gdm.fit.btotal <- gdm(gdmTab.btotal, geo = TRUE)
  
  ## GDM fit list
  gdm.fit.list[gdm.fit.list$beta_metric == "Brep",]$models[[1]][[i]] <- gdm.fit.brep
  gdm.fit.list[gdm.fit.list$beta_metric == "Brich",]$models[[1]][[i]] <- gdm.fit.brich
  gdm.fit.list[gdm.fit.list$beta_metric == "Btotal",]$models[[1]][[i]] <- gdm.fit.btotal
  
  ## Update progress
  setTxtProgressBar(pb, i)
  
}

close(pb)

gc()

## Next we need to deal with the big lists
fit.brep <- gdm.fit.list |> 
  filter(beta_metric == "Brep") |> 
  pull(models) |> 
  unlist(recursive = F)

fit.brich <- gdm.fit.list |> 
  filter(beta_metric == "Brich") |> 
  pull(models) |> 
  unlist(recursive = F)

fit.btotal <- gdm.fit.list |> 
  filter(beta_metric == "Btotal") |> 
  pull(models) |> 
  unlist(recursive = F)


## Summary stats on the fits
summary_stats <- gdm.fit.list %>%
  mutate(
    # First unlist and extract the explained deviance
    deviance_explained = map(models, function(model_list) {
      sapply(model_list, function(x) x$explained)
    })
  ) %>%
  group_by(beta_metric) %>%
  summarise(
    n = n(),
    mean_deviance = mean(unlist(deviance_explained)),
    lower_ci = quantile(unlist(deviance_explained), 0.025),
    upper_ci = quantile(unlist(deviance_explained), 0.975),
    sd_deviance = sd(unlist(deviance_explained))
  )

## Get the effects of each predictor
coeffs <- split(matrix(gdm.fit.list[[1]]$coefficients, ncol = 3, byrow = T), gdm.fit.list[[1]]$predictors)
coeffSum <- lapply(coeffs, sum)
coeffDF <- data.frame(Pred = names(coeffSum), Bsum = unlist(coeffSum))
rownames(coeffDF) <- NULL
coeffDF <- coeffDF |> mutate(Pred = case_when(Pred == "ele" ~ "Elevation",
                                              Pred == "stage" ~ "Stand Age",
                                              Pred == "cancov" ~ "Canopy Cover",
                                              Pred == "Geographic" ~ "Geog. Dist.",
                                              Pred == "fire1_5yr_high_prop" ~ "High Sev. Fire: 1-5yr",
                                              Pred == "fire6_10yr_high_prop" ~ "High Sev. Fire: 6-10yr",
                                              Pred == "fire11_35yr_high_prop" ~ "High Sev. Fire: 11-35yr",
                                              Pred == "fire1_5yr_lowmod_prop" ~ "Low/Mod Sev. Fire: 1-5yr",
                                              Pred == "fire6_10yr_lowmod_prop" ~ "Low/Mod Sev. Fire: 6-10yr",
                                              Pred == "fire11_35yr_lowmod_prop" ~ "Low/Mod Sev. Fire: 11-35yr")) |> 
  mutate(Color = case_when(Pred %in% c("Elevation", "Geog. Dist.") ~ "#7570b3",
                           Pred %in% c("Canopy Cover", "Stand Age") ~ "#1b9e77",
                           Pred %in% c("High Sev. Fire: 1-5yr",
                                       "High Sev. Fire: 6-10yr",
                                       "High Sev. Fire: 11-35yr",
                                       "Low/Mod Sev. Fire: 1-5yr",
                                       "Low/Mod Sev. Fire: 6-10yr",
                                       "Low/Mod Sev. Fire: 11-35yr") ~ "#d95f02"))

## Lollipop for the summed coefficients
coeff_lolli <- ggplot(coeffDF, aes(x = reorder(Pred, Bsum), y = Bsum, color = Pred)) +
  geom_segment(aes(x = reorder(Pred, Bsum),
                   xend = reorder(Pred, Bsum),
                   y = 0, 
                   yend = Bsum),
               size = 1) +
  geom_point(aes(color = Pred), size = 3) +
  scale_color_manual(values = setNames(coeffDF$Color, coeffDF$Pred)) +
  coord_flip() +
  theme_bw() +
  labs(x = "",
       y = "Effect Size") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.position = "none")


## Significance testing via varImp## Significance testing via varImp
gdmVI <- gdm::gdm.varImp(gdmTab, geo = T)

## GDM partition variance
varSet <- vector("list", 3)

names(varSet) <- c("topo", "forest", "fire")

varSet$topo <- c("topo_elev") 
varSet$forest <- c("standage_f3_mn") 
varSet$fire <- c("fire1yr_high_prop",
                 "fire2_5yr_high_prop",
                 "fire6_10yr_high_prop",
                 "fire11_35yr_high_prop")
gdm.part.dev <- gdm::gdm.partition.deviance(gdmTab, varSets = varSet, partSpace = FALSE)


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Plotting
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Plot topo
gdm.spline <- isplineExtract(gdm.fit.list[[1]])

#plot(gdm.spline$x[,"topo_elev"], gdm.spline$y[,"topo_elev"], lwd=3, type="l", xlab="Elevation", ylab="Partial ecological distance")

## Elevation
topo.plot <- gdm.spline |> 
  map(~ as.data.frame(.x)) |>
  map(~ select(.x, ele)) |> 
  bind_cols() |> 
  rename(x = `ele...1`, y = `ele...2`) |> 
  ggplot() + 
  geom_line(aes(x=x, y=y), color = "black", size = 2) +
  theme_bw() + 
  theme(axis.title = element_text(family = "sans", size = 16),
        axis.text = element_text(family = "sans", size = 14)) + 
  scale_y_continuous(limits = c(0, 0.4)) + 
  xlab("Elevation (meters)") + 
  ylab("")

# ggsave(plot = topo.plot,
#        filename = here("Figures/Exploration/GDM_Topo_pareff.jpg"),
#        height = 8, width = 12, dpi = 600)

## Stand Age
sage.plot <- gdm.spline |> 
  map(~ as.data.frame(.x)) |>
  map(~ select(.x, stage)) |> 
  bind_cols() |> 
  rename(x = `stage...1`, y = `stage...2`) |> 
  ggplot() + 
  geom_line(aes(x=x, y=y), color = "black", size = 2) +
  theme_bw() + 
  theme(axis.title = element_text(family = "sans", size = 16),
        axis.text = element_text(family = "sans", size = 14)) +
  scale_y_continuous(limits = c(0, 0.4)) +
  xlab("Stand Age") +   
  ylab("")

## Geographic
geo.plot <- gdm.spline |> 
  map(~ as.data.frame(.x)) |>
  map(~ select(.x, Geographic)) |> 
  bind_cols() |> 
  rename(x = `Geographic...1`, y = `Geographic...2`) |> 
  ggplot() + 
  geom_line(aes(x=x, y=y), color = "black", size = 2) +
  theme_bw() + 
  theme(axis.title = element_text(family = "sans", size = 16),
        axis.text = element_text(family = "sans", size = 14)) +
  scale_y_continuous(limits = c(0, 0.4)) +
  xlab("Distance") +   
  ylab("")

# ggsave(plot = sage.plot,
#        filename = here("Figures/Exploration/GDM_StandAge_pareff.jpg"),
#        height = 8, width = 12, dpi = 600)

## High Severity Fire

fire1_5hs.plot <- gdm.spline |> 
  map(~ as.data.frame(.x)) |>
  map(~ select(.x, fire1_5yr_high_prop)) |> 
  bind_cols() |> 
  rename(x = `fire1_5yr_high_prop...1`, y = `fire1_5yr_high_prop...2`) |> 
  ggplot() + 
  geom_line(aes(x=x, y=y), color = "black", size = 2) +
  theme_bw() + 
  theme(axis.title = element_text(family = "sans", size = 16),
        axis.text = element_text(family = "sans", size = 14)) +
  scale_y_continuous(limits = c(0, 0.4)) +
  xlab("Proportion High Severity Fire: 1-5 years") +  
  ylab("Partial Ecological Distance")

fire6_10hs.plot <- gdm.spline |> 
  map(~ as.data.frame(.x)) |>
  map(~ select(.x, fire6_10yr_high_prop)) |> 
  bind_cols() |> 
  rename(x = `fire6_10yr_high_prop...1`, y = `fire6_10yr_high_prop...2`) |> 
  ggplot() + 
  geom_line(aes(x=x, y=y), color = "black", size = 2) +
  theme_bw() + 
  theme(axis.title = element_text(family = "sans", size = 16),
        axis.text = element_text(family = "sans", size = 14)) +
  scale_y_continuous(limits = c(0, 0.4)) +
  xlab("Proportion High Severity Fire: 6-10 years") +  
  ylab("")

fire11_35hs.plot <- gdm.spline |> 
  map(~ as.data.frame(.x)) |>
  map(~ select(.x, fire11_35yr_high_prop)) |> 
  bind_cols() |> 
  rename(x = `fire11_35yr_high_prop...1`, y = `fire11_35yr_high_prop...2`) |> 
  ggplot() + 
  geom_line(aes(x=x, y=y), color = "black", size = 2) +
  theme_bw() + 
  theme(axis.title = element_text(family = "sans", size = 16),
        axis.text = element_text(family = "sans", size = 14)) +
  scale_y_continuous(limits = c(0, 0.4)) +
  xlab("Proportion High Severity Fire: 11-35 years") +  
  ylab("")

## Low Mod Fire
fire1_5lm.plot <- gdm.spline |> 
  map(~ as.data.frame(.x)) |>
  map(~ select(.x, fire1_5yr_lowmod_prop)) |> 
  bind_cols() |> 
  rename(x = `fire1_5yr_lowmod_prop...1`, y = `fire1_5yr_lowmod_prop...2`) |> 
  ggplot() + 
  geom_line(aes(x=x, y=y), color = "black", size = 2) +
  theme_bw() + 
  theme(axis.title = element_text(family = "sans", size = 16),
        axis.text = element_text(family = "sans", size = 14)) +
  scale_y_continuous(limits = c(0, 0.4)) +
  xlab("Proportion Low/Mod Severity Fire: 1-5 years") +  
  ylab("")

fire6_10lm.plot <- gdm.spline |> 
  map(~ as.data.frame(.x)) |>
  map(~ select(.x, fire6_10yr_lowmod_prop)) |> 
  bind_cols() |> 
  rename(x = `fire6_10yr_lowmod_prop...1`, y = `fire6_10yr_lowmod_prop...2`) |> 
  ggplot() + 
  geom_line(aes(x=x, y=y), color = "black", size = 2) +
  theme_bw() + 
  theme(axis.title = element_text(family = "sans", size = 16),
        axis.text = element_text(family = "sans", size = 14)) +
  scale_y_continuous(limits = c(0, 0.4)) +
  xlab("Proportion Low/Mod Severity Fire: 6-10 years") +  
  ylab("")

fire11_35lm.plot <- gdm.spline |> 
  map(~ as.data.frame(.x)) |>
  map(~ select(.x, fire11_35yr_lowmod_prop)) |> 
  bind_cols() |> 
  rename(x = `fire11_35yr_lowmod_prop...1`, y = `fire11_35yr_lowmod_prop...2`) |> 
  ggplot() + 
  geom_line(aes(x=x, y=y), color = "black", size = 2) +
  theme_bw() + 
  theme(axis.title = element_text(family = "sans", size = 16),
        axis.text = element_text(family = "sans", size = 14)) +
  scale_y_continuous(limits = c(0, 0.4)) +
  xlab("Proportion Low/Mod Severity Fire: 11-35 years") +  
  ylab("")


ggsave(plot = fire.plot,
       filename = here("Figures/Preliminary/GDM_Fire6-10_pareff.jpg"),
       height = 8, width = 12, dpi = 600)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Plotting all 
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Plot this all so that we can visualize each model fit
# Extract splines for all models in gdm.fit.list
all_splines <- map(gdm.fit.list, isplineExtract)

# Create a data frame with all splines and their mean
spline_df <- map(all_splines, function(spline) {
  data.frame(
    x.ele = spline$x[,"ele"],
    y.ele = spline$y[,"ele"],
    x.geo = spline$x[,"Geographic"],
    y.geo = spline$y[,"Geographic"],
    x.sage = spline$x[,"stage"],
    y.sage = spline$y[,"stage"],
    x.cc = spline$x[,"cancov"],
    y.cc = spline$y[,"cancov"],
    x.f15hp = spline$x[,"fire1_5yr_high_prop"],
    y.f15hp = spline$y[,"fire1_5yr_high_prop"],
    x.f610hp = spline$x[,"fire6_10yr_high_prop"],
    y.f610hp = spline$y[,"fire6_10yr_high_prop"],
    x.f1135hp = spline$x[,"fire11_35yr_high_prop"],
    y.f1135hp = spline$y[,"fire11_35yr_high_prop"],
    x.f15lm = spline$x[,"fire1_5yr_lowmod_prop"],
    y.f15lm = spline$y[,"fire1_5yr_lowmod_prop"],
    x.f610lm = spline$x[,"fire6_10yr_lowmod_prop"],
    y.f610lm = spline$y[,"fire6_10yr_lowmod_prop"],
    x.f1135lm = spline$x[,"fire11_35yr_lowmod_prop"],
    y.f1135lm = spline$y[,"fire11_35yr_lowmod_prop"]
  )
}) |> 
  bind_rows(.id = "model")

# Calculate mean spline
mean_spline_ele <- spline_df |> 
  group_by(x.ele) |> 
  summarize(mean_y = mean(y.ele))

# Create the plot
topo.plot <- ggplot() +
  # Individual splines
  geom_line(data = spline_df, 
            aes(x = x.ele, y = y.ele, group = model),
            color = "#7570b3", 
            alpha = 0.5) +
  # Mean spline
  geom_line(data = mean_spline_ele,
            aes(x = x.ele, y = mean_y),
            color = "black",
            size = 1) +
  theme_bw() +
  theme(axis.title = element_text(family = "sans", size = 16),
        axis.text = element_text(family = "sans", size = 14)) +
  scale_y_continuous(limits = c(0, 0.4)) +
  xlab("Elevation (meters)") +
  ylab("Ecological Distance")

topo.plot


# Calculate mean spline
mean_spline_geo <- spline_df |> 
  group_by(x.geo) |> 
  summarize(mean_y = mean(y.geo))

# Create the plot
geo.plot <- ggplot() +
  # Individual splines
  geom_line(data = spline_df, 
            aes(x = x.geo, y = y.geo, group = model),
            color = "#7570b3", 
            alpha = 0.5) +
  # Mean spline
  geom_line(data = mean_spline_geo,
            aes(x = x.geo, y = mean_y),
            color = "black",
            size = 1) +
  theme_bw() +
  theme(axis.title = element_text(family = "sans", size = 16),
        axis.text = element_text(family = "sans", size = 14)) +
  scale_y_continuous(limits = c(0, 0.4)) +
  xlab("Geographic Distance") +
  ylab("Ecological Distance")

geo.plot

# Calculate mean spline
mean_spline_stage <- spline_df |> 
  group_by(x.sage) |> 
  summarize(mean_y = mean(y.sage))

# Create the plot
stage.plot <- ggplot() +
  # Individual splines
  geom_line(data = spline_df, 
            aes(x = x.sage, y = y.sage, group = model),
            color = "#1b9e77", 
            alpha = 0.5) +
  # Mean spline
  geom_line(data = mean_spline_stage,
            aes(x = x.sage, y = mean_y),
            color = "black",
            size = 1) +
  theme_bw() +
  theme(axis.title = element_text(family = "sans", size = 16),
        axis.text = element_text(family = "sans", size = 14)) +
  scale_y_continuous(limits = c(0, 0.4)) +
  xlab("Stand Age") +
  ylab("Ecological Distance")

stage.plot


# Calculate mean spline
mean_spline_cc <- spline_df |> 
  group_by(x.cc) |> 
  summarize(mean_y = mean(y.cc))

# Create the plot
cancov.plot <- ggplot() +
  # Individual splines
  geom_line(data = spline_df, 
            aes(x = x.cc, y = y.cc, group = model),
            color = "#1b9e77", 
            alpha = 0.5) +
  # Mean spline
  geom_line(data = mean_spline_cc,
            aes(x = x.cc, y = mean_y),
            color = "black",
            size = 1) +
  theme_bw() +
  theme(axis.title = element_text(family = "sans", size = 16),
        axis.text = element_text(family = "sans", size = 14)) +
  scale_y_continuous(limits = c(0, 0.4)) +
  xlab("Canopy Cover (%)") +
  ylab("Ecological Distance")

cancov.plot

# Calculate mean spline
mean_spline_f15hp <- spline_df |> 
  group_by(x.f15hp) |> 
  summarize(mean_y = mean(y.f15hp))

# Create the plot
f15hp.plot <- ggplot() +
  # Individual splines
  geom_line(data = spline_df, 
            aes(x = x.f15hp, y = y.f15hp, group = model),
            color = "#d95f02", 
            alpha = 0.5) +
  # Mean spline
  geom_line(data = mean_spline_f15hp,
            aes(x = x.f15hp, y = mean_y),
            color = "black",
            size = 1) +
  theme_bw() +
  theme(axis.title = element_text(family = "sans", size = 16),
        axis.text = element_text(family = "sans", size = 14)) +
  scale_y_continuous(limits = c(0, 0.4)) +
  xlab("Prop. High Sev. Fire 1-5yr") +
  ylab("Ecological Distance")

f15hp.plot

# Calculate mean spline
mean_spline_f610hp <- spline_df |> 
  group_by(x.f610hp) |> 
  summarize(mean_y = mean(y.f610hp))

# Create the plot
f610hp.plot <- ggplot() +
  # Individual splines
  geom_line(data = spline_df, 
            aes(x = x.f610hp, y = y.f610hp, group = model),
            color = "#d95f02", 
            alpha = 0.5) +
  # Mean spline
  geom_line(data = mean_spline_f610hp,
            aes(x = x.f610hp, y = mean_y),
            color = "black",
            size = 1) +
  theme_bw() +
  theme(axis.title = element_text(family = "sans", size = 16),
        axis.text = element_text(family = "sans", size = 14)) +
  scale_y_continuous(limits = c(0, 0.4)) +
  xlab("Prop. High Sev. Fire 6-10yr") +
  ylab("Ecological Distance")

f610hp.plot

# Calculate mean spline
mean_spline_f1135hp <- spline_df |> 
  group_by(x.f1135hp) |> 
  summarize(mean_y = mean(y.f1135hp))

# Create the plot
f1135hp.plot <- ggplot() +
  # Individual splines
  geom_line(data = spline_df, 
            aes(x = x.f1135hp, y = y.f1135hp, group = model),
            color = "#d95f02", 
            alpha = 0.5) +
  # Mean spline
  geom_line(data = mean_spline_f1135hp,
            aes(x = x.f1135hp, y = mean_y),
            color = "black",
            size = 1) +
  theme_bw() +
  theme(axis.title = element_text(family = "sans", size = 16),
        axis.text = element_text(family = "sans", size = 14)) +
  scale_y_continuous(limits = c(0, 0.4)) +
  xlab("Prop. High Sev. Fire 11-35yr") +
  ylab("Ecological Distance")

f1135hp.plot

# Calculate mean spline
mean_spline_f15lm <- spline_df |> 
  group_by(x.f15lm) |> 
  summarize(mean_y = mean(y.f15lm))

# Create the plot
f15lm.plot <- ggplot() +
  # Individual splines
  geom_line(data = spline_df, 
            aes(x = x.f15lm, y = y.f15lm, group = model),
            color = "#d95f02", 
            alpha = 0.5) +
  # Mean spline
  geom_line(data = mean_spline_f15lm,
            aes(x = x.f15lm, y = mean_y),
            color = "black",
            size = 1) +
  theme_bw() +
  theme(axis.title = element_text(family = "sans", size = 16),
        axis.text = element_text(family = "sans", size = 14)) +
  scale_y_continuous(limits = c(0, 0.4)) +
  xlab("Prop. Low/Mod Sev. Fire 1-5yr") +
  ylab("Ecological Distance")

f15lm.plot

# Calculate mean spline
mean_spline_f610lm <- spline_df |> 
  group_by(x.f610lm) |> 
  summarize(mean_y = mean(y.f610lm))

# Create the plot
f610lm.plot <- ggplot() +
  # Individual splines
  geom_line(data = spline_df, 
            aes(x = x.f610lm, y = y.f610lm, group = model),
            color = "#d95f02", 
            alpha = 0.5) +
  # Mean spline
  geom_line(data = mean_spline_f610lm,
            aes(x = x.f610lm, y = mean_y),
            color = "black",
            size = 1) +
  theme_bw() +
  theme(axis.title = element_text(family = "sans", size = 16),
        axis.text = element_text(family = "sans", size = 14)) +
  scale_y_continuous(limits = c(0, 0.4)) +
  xlab("Prop. Low/Mod Sev. Fire 6-10yr") +
  ylab("Ecological Distance")

f610lm.plot

# Calculate mean spline
mean_spline_f1135lm <- spline_df |> 
  group_by(x.f1135lm) |> 
  summarize(mean_y = mean(y.f1135lm))

# Create the plot
f1135lm.plot <- ggplot() +
  # Individual splines
  geom_line(data = spline_df, 
            aes(x = x.f1135lm, y = y.f1135lm, group = model),
            color = "#d95f02", 
            alpha = 0.5) +
  # Mean spline
  geom_line(data = mean_spline_f1135lm,
            aes(x = x.f1135lm, y = mean_y),
            color = "black",
            size = 1) +
  theme_bw() +
  theme(axis.title = element_text(family = "sans", size = 16),
        axis.text = element_text(family = "sans", size = 14)) +
  scale_y_continuous(limits = c(0, 0.4)) +
  xlab("Prop. Low/Mod Sev. Fire 11-35yr") +
  ylab("Ecological Distance")

f1135lm.plot

library(patchwork)

pgrid_lyt <- "
AAB
CDE
FHI
JKL
"
pgrid <- free(coeff_lolli) + topo.plot + stage.plot + cancov.plot + geo.plot + f15hp.plot + f610hp.plot + f1135hp.plot + f15lm.plot + f610lm.plot + f1135lm.plot + plot_layout(design = pgrid_lyt)

ggsave(filename = here("./Figures/Preliminary/GridLolli_GDM_AllSamples_FireOnly.jpg"),
       plot = pgrid, height = 12, width = 12,
       dpi = 600)



pgrid <- topo.plot + sage.plot + geo.plot + fire1_5hs.plot + fire6_10hs.plot + fire11_35hs.plot + fire1_5lm.plot + fire6_10lm.plot + fire11_35lm.plot + plot_layout(guides = "collect")

ggsave(plot = pgrid,
       filename = here("Figures/Preliminary/GDM_All_Responses_Zsample1.jpg"),
       height = 8, width = 16, dpi = 600)

gdm.pred <- predict(object=gdm.fit1, data=gdmTab)

head(gdm.pred)

plot(gdmTab$distance, 
     gdm.pred, 
     xlab="Observed dissimilarity", 
     ylab="Predicted dissimilarity", 
     xlim=c(0,1), 
     ylim=c(0,1), 
     pch=20, 
     col=rgb(0,0,1,0.5))
lines(c(-1,2), c(-1,2))
