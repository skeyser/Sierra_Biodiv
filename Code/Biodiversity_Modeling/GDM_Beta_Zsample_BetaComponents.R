## -------------------------------------------------------------
##
## Script name: GDM Ferrier Z Matrix
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
library(tidyr)
library(sf)

## GDM
library(gdm)
library(adespatial)

## -------------------------------------------------------------

## Load in the Z-matrices
#load(here("./Data/JAGS_Data/Occ2GDM_Data_SpThresh_975minMaxPrec.Rdata"))
load(here("./Data/SpOccupancy_Data/Occ2GDM_Data_SpThresh_975minMaxPrec_spOcc.Rdata"))

## Reformat data
Z <- OccData$Z.posterior
hist(apply(apply(Z, c(1,3), sum), 1, mean))
hist(apply(apply(Z, c(1,3), sum), 1, sd))
hist(apply(apply(Z, c(1,3), sum), 1, function(x) quantile(x, probs = 0.975)) - apply(apply(Z, c(1,3), sum), 1, function(x) quantile(x, probs = 0.025)))
mean(apply(apply(Z, c(1,3), sum), 1, function(x) quantile(x, probs = 0.975)) - apply(apply(Z, c(1,3), sum), 1, function(x) quantile(x, probs = 0.025)))


dimnames(Z) <- list(NULL, OccData$ZcolNames, NULL)

## Take the environmental data
aru_meta <- OccData$SiteMeta
aru_meta <- aru_meta |> 
  select(-utme, -utmn) |> 
  st_as_sf(coords = c("X", "Y"), crs = 4326) |> 
  st_transform(crs = 3310) |> 
  mutate(Lat = st_coordinates(geometry)[,2],
         Long = st_coordinates(geometry)[,1]) |> 
  st_drop_geometry()

hist(aru_meta$Long)
hist(aru_meta$Lat)

## Vars we want
v.keep <- c("Cell_Unit", "Lat", "Long", 
            "ele", "ppt", 
            "fire1_5yr_cbi_mn", "fire6_10yr_cbi_mn", "fire11_35yr_cbi_mn",
            "cancov", "ch_res")

aru_meta <- aru_meta |> 
  select(all_of(v.keep))


## -------------------------------------------------------------
##
## Begin Section: Format for GDM
##
## -------------------------------------------------------------

## Posterior samples
n.samples <- OccData$Z.draws

## Species data
mr <- vector(mode = "list", length = n.samples)
for(i in 1:n.samples){
  Z.tmp <- Z[,,i]
  missRows <- which(rowSums(Z.tmp) == 0)
  mr[[i]] <- missRows
}

## Remove rows with 0 species
missRows <- unique(unlist(mr))
Z <- Z[-missRows,,]

## Environmental data columns
envTab <- aru_meta |>
  as.data.frame()
envTab <- envTab[-missRows, ]

## Format the data for GDM
## !! This section is for testing !! ##
## !! Need >=1000 runs for full on different comp !! ##
samp <- sort(ceiling(runif(100, 0, 1000)))
Z <- Z[,,samp]
npost <- dim(Z)[3]

## Make the list to hold model fits
gdm.fit.list <- tibble(
  beta_metric = c("Brep", "Brich", "Btotal"),
  models = vector("list", 3)) |>
  mutate(models = map(1:n(), ~vector("list", npost)))

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Pairwise Correlations from Mokany et al., 2022
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Set up site-pair table, environmental tabular data 
Z.samp <- as.data.frame(Z[,,1])
rownames(Z.samp) <- envTab$Cell_Unit
envTab$Cell_Unit <- 1:length(unique(envTab$Cell_Unit))

## Partition the beta diversity based on Jaccard-based Podani metrics
bpart <- adespatial::beta.div.comp(Z.samp, coef = "J")
bpart.bas <- adespatial::beta.div.comp(Z.samp, coef = "BJ")

## Total Beta
btotal <- as.matrix(bpart$D)
dimnames(btotal) <- list(envTab$Cell_Unit, envTab$Cell_Unit)
btotal <- cbind(Cell_Unit = envTab$Cell_Unit, btotal)

## Format the data for pairs - bioFormat = 3 accepts preformed D-mat
gdmTab.samp <- formatsitepair(bioData = btotal,
                              bioFormat = 3,
                              XColumn = "Long",
                              YColumn = "Lat",
                              siteColumn = "Cell_Unit",
                              predData = envTab)

# Calculate correlation coefficient for predictor variables  
env.preds <- c("ele","ppt","cancov","ch_res","fire1_5yr_cbi_mn","fire6_10yr_cbi_mn","fire11_35yr_cbi_mn") 
env.cor <- cor(envTab[,which(colnames(envTab) %in% env.preds)], method = "pearson")  
# Calculate correlation coefficients for site-pairs  
env.pair.dif <- matrix(0, ncol = length(env.preds), nrow = nrow(gdmTab.samp)) 
colnames(env.pair.dif) <- paste0(env.preds,"_diff") 
for(i in 1:length(env.preds)) {  
  env.pair.dif[,i] <- abs(gdmTab.samp[,which(colnames(gdmTab.samp) == paste0("s1.",env.preds[i]))] - 
                          gdmTab.samp[,which(colnames(gdmTab.samp) == paste0("s2.",env.preds[i]))]) 
  } # end for i  
env.dif.cor <- cor(env.pair.dif, method = "pearson")  
# compare the correlation coefficients for the predictors vs the differnece between the predictors  
plot(abs(env.cor),
     abs(env.dif.cor), 
     xlab="Predictor correlation", 
     ylab="Predictor-pair correlation")
lines(c(0,1),c(0,1),lty=3)

## -------------------------------------------------------------
##
## End Section: Pairwise correlation checks
##
## -------------------------------------------------------------

## Remove the OccData object prior to loop
rm(OccData)
gc()

## Progress bar
pb <- txtProgressBar(min = 0, max = npost, style = 3)

## Loop to fit GDMs for components of beta diversity
for(i in 1:npost){
  
  ## Pluck one iteration of the posterior
  Z.samp <- as.data.frame(Z[,,i])
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
  
  ## Format the data for pairs - bioFormat = 3 accepts preformed D-mat
  gdmTab.brep <- formatsitepair(bioData = brep,
                                bioFormat = 3,
                                XColumn = "Long",
                                YColumn = "Lat",
                                siteColumn = "Cell_Unit",
                                predData = envTab)
  
  gdmTab.brich <- formatsitepair(bioData = brich,
                                 bioFormat = 3,
                                 XColumn = "Long",
                                 YColumn = "Lat",
                                 siteColumn = "Cell_Unit",
                                 predData = envTab)
  
  gdmTab.btotal <- formatsitepair(bioData = btotal,
                                  bioFormat = 3,
                                  XColumn = "Long",
                                  YColumn = "Lat",
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

## Clean up
gc()

## Write the model list to file
#saveRDS(gdm.fit.list, file = here("./Data/GDM_Out/GDMPosteriorFitsBetaComps_SpVarThresh_975minMaxPrec.RDS"))
gdm.fit.list <- readRDS(file = here("./Data/GDM_Out/GDMPosteriorFitsBetaComps_SpVarThresh_975minMaxPrec.RDS"))

AICFxn <- function(model){ 
  mod <- glm((1-model$observed)~model$ecological, family=binomial(link=log))  
  k <- length(which(model$coefficients>0))+1 
  # Number of coefficients, plus 1 for the model intercept 
  AIC <- (2*k)-(2*logLik(mod)) 
  dev <- ((mod$null.deviance-mod$deviance)/mod$null.deviance)*100  
  return(list(AIC,dev)) 
} 

test <- readRDS("R:/Users/skeyser/Postdoc/GDM_Output/Output_GDM_1.RDS")

test <- test$models[[1]]

AICFxn(test)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Variable Importance
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class(gdm.fit.list$models[[3]][[1]])
VI <- gdm.varImp(gdm.fit.list$models[[3]][[1]],
           predSelect = T,
           nPerm = 50,
           parallel = TRUE,
           cores = 4)

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
summary_stats

## Get the effects of each predictor
## Function to extract coefficients from one model
extract_coeffs <- function(model) {
  # Get the number of predictors
  n_pred <- length(model$predictors)
  
  # Extract coefficients matrix (each row is a predictor, columns are splines)
  coef_mat <- matrix(model$coefficients, nrow = n_pred, ncol = 3, byrow = TRUE)
  
  # Create data frame ensuring predictors match their coefficients
  data.frame(
    Pred = model$predictors,
    spline1 = coef_mat[,1],
    spline2 = coef_mat[,2],
    spline3 = coef_mat[,3]
  )
}

# Let's verify with a single model first
test_extract <- extract_coeffs(gdm.fit.list$models[[1]][[1]])
print(test_extract) 

# Now process all models
coef_summary <- gdm.fit.list %>%
  rowwise() %>%
  mutate(
    # Extract coefficients from all models in the list
    all_coeffs = list(
      map(models, extract_coeffs) %>%
        bind_rows(.id = "model_id")
    )
  ) %>%
  ungroup() %>%
  unnest(all_coeffs) %>%
  # Group and summarize
  group_by(beta_metric, Pred) %>%
  summarise(
    spline1_mean = mean(spline1),
    spline1_high = quantile(spline1, probs = 0.975),
    spline1_low = quantile(spline1, probs = 0.025),
    spline2_mean = mean(spline2),
    spline2_high = quantile(spline2, probs = 0.975),
    spline2_low = quantile(spline2, probs = 0.025),
    spline3_mean = mean(spline3),
    spline3_high = quantile(spline3, probs = 0.975),
    spline3_low = quantile(spline3, probs = 0.025),
    .groups = "drop"
  ) |> mutate(PredPretty = case_when(Pred == "ele" ~ "Elevation",
                            #Pred == "stage" ~ "Stand Age",
                            Pred == "cancov" ~ "Canopy Cover",
                            Pred == "ch_res" ~ "Canopy Height Res.",
                            Pred == "Geographic" ~ "Geog. Dist.",
                            #Pred == "fire1_5yr_high_prop" ~ "High Sev. Fire: 1-5yr",
                            #Pred == "fire6_10yr_high_prop" ~ "High Sev. Fire: 6-10yr",
                            #Pred == "fire11_35yr_high_prop" ~ "High Sev. Fire: 11-35yr",
                            #Pred == "fire1_5yr_lowmod_prop" ~ "Low/Mod Sev. Fire: 1-5yr",
                            #Pred == "fire6_10yr_lowmod_prop" ~ "Low/Mod Sev. Fire: 6-10yr",
                            #Pred == "fire11_35yr_lowmod_prop" ~ "Low/Mod Sev. Fire: 11-35yr"),
                            Pred == "ppt" ~ "Precipitation",
                            Pred == "fire1_5yr_cbi_mn" ~ "Mean Fire Severity: 1-5yr",
                            Pred == "fire6_10yr_cbi_mn" ~ "Mean Fire Severity: 6-10yr",
                            Pred == "fire11_35yr_cbi_mn" ~ "Mean Fire Severity: 11-35yr"
                            )) |> 
  mutate(Color = case_when(PredPretty %in% c("Elevation", 
                                             "Geog. Dist.",
                                             "Precipitation") ~ "#7570b3",
                           PredPretty %in% c("Canopy Cover", 
                                             #"Stand Age"
                                             "Canopy Height Res."
                                             ) ~ "#1b9e77",
                           PredPretty %in% c(
                             "Mean Fire Severity: 1-5yr",
                             "Mean Fire Severity: 6-10yr",
                             "Mean Fire Severity: 11-35yr"
                                       # "High Sev. Fire: 1-5yr",
                                       # "High Sev. Fire: 6-10yr",
                                       # "High Sev. Fire: 11-35yr",
                                       # "Low/Mod Sev. Fire: 1-5yr",
                                       # "Low/Mod Sev. Fire: 6-10yr",
                                       # "Low/Mod Sev. Fire: 11-35yr"
                                       ) ~ "#d95f02"))


## Total spline influence
coef_totals <- coef_summary %>%
  mutate(
    total_mean = spline1_mean + spline2_mean + spline3_mean,
    total_low = spline1_low + spline2_low + spline3_low,
    total_high = spline1_high + spline2_high + spline3_high
  ) |> 
  select(beta_metric, Pred, PredPretty, Color, contains("total")) |> 
  mutate(beta_metric = factor(case_when(
    beta_metric == "Brep" ~ "Turnover",
    beta_metric == "Brich" ~ "Nestedness",
    beta_metric == "Btotal" ~ "Total Beta"
  ), levels = c("Total Beta", "Turnover", "Nestedness")))

## Lollipop for the summed coefficients
new_labels <- c(
  'Brep' = 'Turnover',
  'Brich' = 'Nestedness',
  'Btotal' = 'Total Beta'
)

## Coeff totals difference in effects
coef_totals |> 
  group_by(beta_metric) |> 
  arrange(desc(total_mean))

coeff_lolli <- ggplot(coef_totals, aes(x = reorder(PredPretty, total_mean), y = total_mean, color = PredPretty)) +
  geom_segment(aes(x = reorder(PredPretty, total_mean),
                   xend = reorder(PredPretty, total_mean),
                   y = 0, 
                   yend = total_mean),
               size = 1) +
  geom_point(aes(color = PredPretty), size = 3) +
  # geom_errorbar(aes(ymin = total_low,
  #                   ymax = total_high),
  #               width = 0.2) +
  scale_color_manual(values = setNames(coef_totals$Color, coef_totals$PredPretty)) +
  coord_flip() +
  facet_wrap(~beta_metric, labeller = labeller(beta_metric = new_labels)) +
  scale_y_continuous(breaks = seq(0, 1, by = .2),  # Adjust these numbers as needed
                     limits = c(0, 1)) +  # Adjust limits as needed
  theme_bw() +
  labs(x = "",
       y = "Effect Size") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.position = "none",
        strip.background = element_rect(fill = "lightgray"),
        panel.spacing = unit(0.75, "lines")) 

coeff_lolli

coeff_lolli2 <- ggplot(coef_totals, aes(x = reorder(PredPretty, total_mean), 
                                                       y = total_mean)) +
  # Draw segments
  geom_segment(aes(x = reorder(PredPretty, total_mean),
                   xend = reorder(PredPretty, total_mean),
                   y = 0, 
                   yend = total_mean),
               color = "gray50",
               size = 1) +
  # Draw white background for points to "break" the line
  geom_point(size = 5, 
             color = "white",
             fill = "white") +
  # Draw actual points
  geom_point(aes(color = Color,
                 shape = beta_metric),
             size = 4) +
  scale_color_identity() +
  scale_shape_manual(values = c("Total Beta" = 16,  # Filled circle
                                "Turnover" = 18,      # Open circle
                                "Nestedness" = 17),   # Open triangle
                     name = "Metric") +
  coord_flip() +
  scale_y_continuous(breaks = seq(0, 1, by = .2),
                     limits = c(0, 1)) +
  theme_bw() +
  labs(x = "",
       y = "Effect Size") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.position = "right",
        panel.spacing = unit(0.75, "lines"))

coeff_lolli2

# ggsave(plot = coeff_lolli, filename = here("./Figures/Preliminary/BComponentsLolli.jpg"),
#        height = 8, width = 12, dpi = 600)

## Significance testing via varImp
gdmVI <- gdm::gdm.varImp(gdmTab, geo = T)

## GDM partition variance
varSet <- vector("list", 3)

names(varSet) <- c("topo", "forest", "fire")

varSet$topo <- c("ele") 
varSet$forest <- c("stage", "cancov") 
varSet$fire <- c("fire1_5yr_high_prop",
                 "fire6_10yr_high_prop",
                 "fire11_35yr_high_prop",
                 "fire1_5yr_lowmod_prop",
                 "fire6_10yr_lowmod_prop",
                 "fire11_35yr_lowmod_prop")
gdm.part.dev <- gdm::gdm.partition.deviance(gdmTab, varSets = varSet, partSpace = FALSE)

saveRDS(gdmVI, file = here("./Data/GDM_Out/GDMVarImpSingleFullData.RData"))
saveRDS(gdm.fit.list, file = here("./Data/GDM_Out/GDMPosteriorFits.RData"))

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Automating plotting functions for I-splines 
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ## Extract splines for all models from the list of posterior models
# all_splines <- gdm.fit.list %>%
#   rowwise() %>%
#   mutate(
#     splines = list(map(models, ~isplineExtract(.)))
#   ) %>%
#   ungroup()
# 
# # For a single variable (example with elevation):
# spline_df <- all_splines %>%
#   filter(beta_metric == "Brep") %>%  # or whichever metric you want
#   pull(splines) %>%
#   .[[1]] %>%
#   map_dfr(function(spline) {
#     data.frame(
#       x = spline$x[,"ele"],
#       y = spline$y[,"ele"]
#     )
#   }, .id = "model")

## Function to extract I-spline matrix for a given variable
make_spline_df <- function(model_list, var_name){
  spline_df <- gdm.fit.list |>
  rowwise() |>
  mutate(
    spline_data = list(
      map_dfr(models, function(model) {
        spline <- isplineExtract(model)
        data.frame(
          x = spline$x[,var_name],
          y = spline$y[,var_name]
        )
      }, .id = "model")
    )
  ) |>
  unnest(spline_data)
  
  ## Find the covariate name
  cov_rename <- data.frame(Pred = c("ele", 
                                    "cancov", 
                                    "Geographic", 
                                    "ch_res", 
                                    "ppt", 
                                    #"stage", 
                                    #"fire1_5yr_high_prop", "fire6_10yr_high_prop", "fire11_35yr_high_prop", 
                                    #"fire1_5yr_lowmod_prop", "fire6_10yr_lowmod_prop", "fire11_35yr_lowmod_prop"
                                    "fire1_5yr_cbi_mn",
                                    "fire6_10yr_cbi_mn",
                                    "fire11_35yr_cbi_mn"
                                    ), 
                           PredPretty = c("Elevation", 
                                          "Canopy Cover",
                                          "Geog. Dist.", 
                                          "Canopy Height Res.",
                                          "Precipitation",
                                          # "Stand Age",
                                          # "High Sev. Fire: 1-5yr",
                                          # "High Sev. Fire: 6-10yr",
                                          # "High Sev. Fire: 11-35yr",
                                          # "Low/Mod Sev. Fire: 1-5yr",
                                          # "Low/Mod Sev. Fire: 6-10yr",
                                          # "Low/Mod Sev. Fire: 11-35yr"
                                          "Mean Fire Severity: 1-5yr",
                                          "Mean Fire Severity: 6-10yr",
                                          "Mean Fire Severity: 11-35yr"
                                          )) |> 
    mutate(Color = case_when(PredPretty %in% c("Elevation", 
                                               "Geog. Dist.",
                                               "Precipitation"
                                               ) ~ "#7570b3",
                             PredPretty %in% c("Canopy Cover", 
                                               #"Stand Age"
                                               "Canopy Height Res."
                                               ) ~ "#1b9e77",
                             PredPretty %in% c("Mean Fire Severity: 1-5yr",
                                               "Mean Fire Severity: 6-10yr",
                                               "Mean Fire Severity: 11-35yr"
                                               # "High Sev. Fire: 1-5yr",
                                               # "High Sev. Fire: 6-10yr",
                                               # "High Sev. Fire: 11-35yr",
                                               # "Low/Mod Sev. Fire: 1-5yr",
                                               # "Low/Mod Sev. Fire: 6-10yr",
                                               # "Low/Mod Sev. Fire: 11-35yr"
                                               ) ~ "#d95f02"))
  
  ## Return a list with the necessary info to automate plotting
  spline_df <- list(spline_df = spline_df,
                    plot_name = cov_rename$PredPretty[cov_rename$Pred == var_name],
                    color = cov_rename$Color[cov_rename$Pred == var_name])
  
  return(spline_df)
}

## Function for plotting
spline_plot <- function(spline_list){
  
  ## Take the DF for the plot
  spline_df <- spline_list$spline_df
  
  ## Convert m to km for geog dist
  if(spline_list$plot_name == "Geog. Dist.") {
    spline_df$x <- spline_df$x/1000
  }
  
  # Calculate mean spline
  mean_spline <- spline_df |> 
    group_by(beta_metric, x) |> 
    summarize(mean_y = mean(y)) |> 
    mutate(beta_metric = factor(case_when(
      beta_metric == "Brep" ~ "Turnover",
      beta_metric == "Brich" ~ "Nestedness",
      beta_metric == "Btotal" ~ "Total Beta"
    ), levels = c("Total Beta", "Turnover", "Nestedness")))
  
  ## Calculate the 95% CI for splines
  spline_ci <- spline_df |> 
    group_by(beta_metric, x) |> 
    summarize(hi_y = quantile(y, probs = 0.975),
              lo_y = quantile(y, probs = 0.025)) |> 
    mutate(beta_metric = factor(case_when(
      beta_metric == "Brep" ~ "Turnover",
      beta_metric == "Brich" ~ "Nestedness",
      beta_metric == "Btotal" ~ "Total Beta"
    ), levels = c("Total Beta", "Turnover", "Nestedness")))
  
  splot <- ggplot() +
    # Individual splines
    # geom_line(data = spline_df, 
    #           aes(x = x, y = y, group = model),
    #           color = color, 
    #           alpha = 0.5) +
    # Mean spline
    geom_line(data = mean_spline,
              aes(x = x, 
                  y = mean_y,
                  linetype = beta_metric),
              color = spline_list$color,
              size = 1) +
    scale_linetype_manual(values = c("Total Beta" = "solid",
                                     "Turnover" = "dashed",
                                     "Nestedness" = "dotted"),
                          name = "Component") +
    scale_fill_manual(values = c("Total Beta" = "gray70",
                                 "Turnover" = "gray70", 
                                 "Nestedness" = "gray70"),
                      name = "Component") +
    geom_ribbon(data = spline_ci,
                aes(x = x, 
                    ymin = lo_y,
                    ymax = hi_y,
                    group = beta_metric),
                alpha = 0.3,
                fill = spline_list$color) +
    theme_bw() +
    theme(axis.title = element_text(family = "sans", size = 16),
          axis.text = element_text(family = "sans", size = 14)) +
    scale_y_continuous(limits = c(0, 1)) +
    xlab(spline_list$plot_name) +
    ylab("Ecological Distance") +
    labs(linetype = NULL) +
    guides(linetype = guide_legend(
      override.aes = list(color = "black")  # Make legend lines black
    ))

  ## Return the plot
return(splot)
}

## Plot each of the variable for patchworking later
## Make the DFs
ele_df <- make_spline_df(model_list = gdm.fit.list, var = "ele")
ppt_df <- make_spline_df(model_list = gdm.fit.list, var = "ppt")
cancov_df <- make_spline_df(model_list = gdm.fit.list, var = "cancov")
canht_df <- make_spline_df(model_list = gdm.fit.list, var = "ch_res")
geo_df <- make_spline_df(model_list = gdm.fit.list, var = "Geographic")
f5h_df <- make_spline_df(model_list = gdm.fit.list, var = "fire1_5yr_cbi_mn")
f10h_df <- make_spline_df(model_list = gdm.fit.list, var = "fire6_10yr_cbi_mn")
f35h_df <- make_spline_df(model_list = gdm.fit.list, var = "fire11_35yr_cbi_mn")
# f5h_df <- make_spline_df(model_list = gdm.fit.list, var = "fire1_5yr_high_prop")
# f10h_df <- make_spline_df(model_list = gdm.fit.list, var = "fire6_10yr_high_prop")
# f35h_df <- make_spline_df(model_list = gdm.fit.list, var = "fire11_35yr_high_prop")
# f5lm_df <- make_spline_df(model_list = gdm.fit.list, var = "fire1_5yr_lowmod_prop")
# f10lm_df <- make_spline_df(model_list = gdm.fit.list, var = "fire6_10yr_lowmod_prop")
# f35lm_df <- make_spline_df(model_list = gdm.fit.list, var = "fire11_35yr_lowmod_prop")

## Feed into the plotting function
p.ele <- spline_plot(ele_df)
p.ppt <- spline_plot(ppt_df)
p.cc <- spline_plot(cancov_df)
p.ch <- spline_plot(canht_df)
#p.sage <- spline_plot(stage_df)
p.geo <- spline_plot(geo_df)
p.f5h <- spline_plot(f5h_df)
p.f10h <- spline_plot(f10h_df)
p.f35h <- spline_plot(f35h_df)
# p.f5lm <- spline_plot(f5lm_df)
# p.f10lm <- spline_plot(f10lm_df)
# p.f35lm <- spline_plot(f35lm_df)


library(patchwork)

pgrid_lyt <- "
AAB
CDE
FHI
J##
"
pgrid <- free(coeff_lolli2) + 
  p.geo + 
  p.ele +
  p.ppt +
  #p.sage + 
  p.cc +
  p.ch +
  p.f5h + 
  p.f10h + 
  p.f35h + 
  #p.f5lm + 
  #p.f10lm + 
  #p.f35lm + 
  plot_layout(design = pgrid_lyt, guides = "collect")

ggsave(filename = here("./Figures/Preliminary/GridLolli_GDM_AllSamples_BetaComponents_SpVarThresh.jpg"),
       plot = pgrid, height = 13, width = 13,
       dpi = 600)

