## -------------------------------------------------------------
##
## Script name: JAGS Output Processing
##
## Script purpose:
##
## Author: Spencer R Keyser
##
## Date Created: 2025-01-14
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
library(stringr)
library(ggplot2)
library(here)
library(MCMCvis)

## -------------------------------------------------------------

logit2prob <- function(logit){
  odds <- exp(logit)
  prob <- odds / (1 + odds)
  return(prob)
}

## -------------------------------------------------------------
##
## Begin Section: Load in model output from JAGS runs
##
## -------------------------------------------------------------

## Load runs
load("R:/Users/skeyser/Postdoc/MSOM_Ragged_JAGS_Summaries_95thresh.Rdata")
load("R:/Users/skeyser/Postdoc/MSOM_Ragged_JAGS_Zout_95thresh.Rdata")

## Model summary
MCMCvis::MCMCsummary(out2, round = 2)



MCMCvis::MCMCtrace(out2, params = c("mu.beta1"), n.eff = T)
MCMCvis::MCMCtrace(out2, params = c("mu.beta2"), n.eff = T)

## Tree plot
out2$parameters

## Community-wide effects
## Occupancy model
MCMCvis::MCMCplot(out2,
         params = c(paste0("mu.beta", seq(1,9))),
         ci = c(50,90),
         ref_ovl = TRUE,
         labels = c("Latitude", 
                    "Elevation", 
                    expression("Elevation"^2), 
                    "Mean CBI: 1 yr post",
                    "Mean CBI: 2-5 yr post",
                    "Mean CBI: 6-10 yr post",
                    "Mean CBI: 11-35 yr post",
                    "Stand Age",
                    "Canopy Cover"),
         rank = T)

## Detection model
MCMCvis::MCMCplot(out2,
                  params = c(paste0("mu.alpha", seq(1,3))),
                  ci = c(50,90),
                  ref_ovl = TRUE,
                  labels = c("Efforts Hours",
                             "JDate",
                             expression("JDate"^2)),
                  rank = T)

par(mfrow = c(1,2))
psi.sample <- plogis(rnorm(10^6, mean = out2$mean$mu.lpsi, sd = out2$mean$sd.lpsi))
p.sample <- plogis(rnorm(10^6, mean = out2$mean$mu.lp, sd = out2$mean$sd.lp))

hist(psi.sample, freq = F, breaks = 50, col = "grey", xlab = "Species occupancy
probability", ylab = "Density", main = "")

hist(p.sample, freq = F, breaks = 50, col = "grey", xlab = "Species detection probability",
     ylab = "Density", main = "")
summary(psi.sample) ; summary(p.sample)

## Species-specific responses
sp.index <- read.csv(here("Code/Occupancy_Modeling/SpeciesIndex_Filtered.csv"))
sp.index$Index <- as.character(sp.index$Index)

str(out3)
all3 <- as.matrix(out3)
rm(out3)
gc()


pm <- apply(all3, 2, mean)
cri <- apply(all3, 2, function(x) quantile(x, prob = c(0.025, 0.975)))

nspec <- 95
npar <- 14
N <- nspec*npar
sp.resp <- data.frame(Par = names(pm[1:N]),
                      Mean = pm[1:N],
                      UCI = cri[2,1:N],
                      LCI = cri[1,1:N])
sp.resp <- sp.resp |> 
  mutate(Species = str_extract(Par, "\\[\\d+\\]")) |> 
  mutate(Species = gsub("[[:punct:]]", "", Species)) |> 
  left_join(sp.index, by = c("Species" = "Index")) |> 
  select(Par, Mean, UCI, LCI, Species = Species.y) |> 
  mutate(Par = gsub("\\[\\d+\\]", "", Par))


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Detection Covariates
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Detecton Effort
labels_det <- c(
  "Survey Hours" = "Survey~Hours",
  "JDate" = "JDate",
  "JDate2" = "JDate^2"
)

sp.resp |> 
  filter(str_detect(Par, "alpha")) |> 
  mutate(ParPretty = case_when(Par == "alpha1" ~ "Survey Hours",
                               Par == "alpha2" ~ "JDate",
                               Par == "alpha3" ~ "JDate2")) |>
  ggplot(aes(y = Mean, x = reorder(Species, -Mean))) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 1.2) +
  geom_pointrange(aes(ymin = LCI, ymax = UCI,
                      color = case_when(
                        LCI > 0 ~ "Above",
                        UCI < 0 ~ "Below",
                        TRUE ~ "Intersects")), 
                  size = 0.8,
                  show.legend = F) +
  scale_color_manual(values = c("Above" = "#99CCFF",
                                "Below" = "#FF6600",
                                "Intersect" = "gray")) + 
  coord_flip() + 
  theme_bw() +
  xlab("Species") + 
  ylab("Parameter Estimate") + 
  facet_wrap(~ParPretty)

## Overall detection estimates by species
sp.resp |> 
  filter(Par == "lp") |> 
  mutate(Mean = plogis(Mean), UCI = plogis(UCI), LCI = plogis(LCI)) |> 
  ggplot(aes(y = Mean, x = reorder(Species, -Mean))) +
  #geom_hline(yintercept = mean(Mean), linetype = "dashed", color = "black", size = 1.2) +
  geom_pointrange(aes(ymin = LCI, ymax = UCI,
                      color = case_when(
                        LCI > 0 ~ "Above",
                        UCI < 0 ~ "Below",
                        TRUE ~ "Intersects")), 
                  size = 0.8,
                  show.legend = F) +
  scale_color_manual(values = c("Above" = "#99CCFF",
                                "Below" = "#FF6600",
                                "Intersect" = "gray")) + 
  coord_flip() + 
  theme_bw() +
  xlab("Species") + 
  ylab("Parameter Estimate")

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Occupancy Covariates
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Overall occupancy estimates by species
sp.resp |> 
  filter(Par == "lpsi") |> 
  mutate(Mean = plogis(Mean), UCI = plogis(UCI), LCI = plogis(LCI)) |> 
  ggplot(aes(y = Mean, x = reorder(Species, -Mean))) +
  #geom_hline(yintercept = mean(Mean), linetype = "dashed", color = "black", size = 1.2) +
  geom_pointrange(aes(ymin = LCI, ymax = UCI,
                      color = case_when(
                        LCI > 0 ~ "Above",
                        UCI < 0 ~ "Below",
                        TRUE ~ "Intersects")), 
                  size = 0.8,
                  show.legend = F) +
  scale_color_manual(values = c("Above" = "#99CCFF",
                                "Below" = "#FF6600",
                                "Intersect" = "gray")) + 
  coord_flip() + 
  theme_bw() +
  xlab("Species") + 
  ylab("Parameter Estimate")

## Geo/Topo
labels <- c(
  "Latitude" = "Latitude",
  "Elevation" = "Elevation",
  "Elevation2" = "Elevation^2"
)

geo.eff.plot <- sp.resp |> 
  filter(str_detect(Par, "beta")) |> 
  mutate(ParPretty = case_when(Par == "beta1" ~ "Latitude",
                               Par == "beta2" ~ "Elevation",
                               Par == "beta3" ~ "Elevation2",
                               Par == "beta4" ~ "Mean CBI: 1 yr post",
                               Par == "beta5" ~ "Mean CBI: 2-5 yr post",
                               Par == "beta6" ~ "Mean CBI: 6-10 yr post",
                               Par == "beta7" ~ "Mean CBI: 11-35 yr post",
                               Par == "beta8" ~ "Stand Age",
                               Par == "beta9" ~ "Canopy Cover")) |>
  filter(ParPretty %in% c("Latitude", "Elevation", "Elevation2")) |> 
  ggplot(aes(y = Mean, x = reorder(Species, -Mean))) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 1.2) +
  geom_pointrange(aes(ymin = LCI, ymax = UCI,
                      color = case_when(
                        LCI > 0 ~ "Above",
                        UCI < 0 ~ "Below",
                        TRUE ~ "Intersects")), 
                  size = 0.8,
                  show.legend = F) +
  scale_color_manual(values = c("Above" = "#FF6600",
                                "Below" = "#99CCFF",
                                "Intersect" = "gray")) + 
  coord_flip() + 
  theme_bw() +
  xlab("Species") + 
  ylab("Parameter Estimate") + 
  facet_wrap(~ParPretty, labeller = as_labeller(labels, default = label_parsed))

## Habitat
hab.eff.plot <- sp.resp |> 
  filter(str_detect(Par, "beta")) |> 
  mutate(ParPretty = case_when(Par == "beta1" ~ "Latitude",
                               Par == "beta2" ~ "Elevation",
                               Par == "beta3" ~ "Elevation2",
                               Par == "beta4" ~ "Mean CBI: 1 yr post",
                               Par == "beta5" ~ "Mean CBI: 2-5 yr post",
                               Par == "beta6" ~ "Mean CBI: 6-10 yr post",
                               Par == "beta7" ~ "Mean CBI: 11-35 yr post",
                               Par == "beta8" ~ "Stand Age",
                               Par == "beta9" ~ "Canopy Cover")) |>
  filter(ParPretty %in% c("Stand Age", "Canopy Cover")) |> 
  ggplot(aes(y = Mean, x = reorder(Species, -Mean))) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 1.2) +
  geom_pointrange(aes(ymin = LCI, ymax = UCI,
                      color = case_when(
                        LCI > 0 ~ "Above",
                        UCI < 0 ~ "Below",
                        TRUE ~ "Intersects")), 
                  size = 0.8,
                  show.legend = F) +
  scale_color_manual(values = c("Above" = "#FF6600",
                                "Below" = "#99CCFF",
                                "Intersect" = "gray")) +
  coord_flip() + 
  theme_bw() +
  xlab("Species") + 
  ylab("Parameter Estimate") + 
  facet_wrap(~ParPretty, nrow = 1)

## Fire
fire.eff.plot <- sp.resp |> 
  filter(str_detect(Par, "beta")) |> 
  mutate(ParPretty = case_when(Par == "beta1" ~ "Latitude",
                               Par == "beta2" ~ "Elevation",
                               Par == "beta3" ~ "Elevation2",
                               Par == "beta4" ~ "Mean CBI: 1 yr post",
                               Par == "beta5" ~ "Mean CBI: 2-5 yr post",
                               Par == "beta6" ~ "Mean CBI: 6-10 yr post",
                               Par == "beta7" ~ "Mean CBI: 11-35 yr post",
                               Par == "beta8" ~ "Stand Age",
                               Par == "beta9" ~ "Canopy Cover")) |>
  filter(str_detect(ParPretty, "Mean CBI")) |> 
  mutate(ParPretty = factor(ParPretty,
                               levels = c(
                                 "Mean CBI: 1 yr post",
                               "Mean CBI: 2-5 yr post",
                               "Mean CBI: 6-10 yr post",
                               "Mean CBI: 11-35 yr post"))) |> 
  ggplot(aes(y = Mean, x = reorder(Species, -Mean))) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 1.2) +
  geom_pointrange(aes(ymin = LCI, ymax = UCI,
                      color = case_when(
                        LCI > 0 ~ "Above",
                        UCI < 0 ~ "Below",
                        TRUE ~ "Intersects")), 
                  size = 0.8,
                  show.legend = F) +
  scale_color_manual(values = c("Above" = "#FF6600",
                                "Below" = "#99CCFF",
                                "Intersect" = "gray")) + 
  coord_flip() + 
  theme_bw() +
  xlab("Species") + 
  ylab("Parameter Estimate") + 
  facet_wrap(~ParPretty, nrow = 1)

ggsave(plot = fire.eff.plot, filename = here("./Figures/Preliminary/MSOM_EffPlot_FireCovariates.jpg"),
       height = 10, width = 10, dpi = 600)


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Beta-diversity with MSOM
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
aru_meta <- readr::read_csv(here("Data/ARU_120m.csv"))
aru_meta$Cell_Unit <- paste0(aru_meta$cell_id, "_", aru_meta$unit_numbe)

## Select the variables we are interested in for detection and occupancy
aru_meta <- aru_meta |> select(Cell_Unit, 
                               utme, utmn,
                               X, Y,
                               topo_elev, topo_tpi,
                               tmx_bcm_mn, ppt_bcm_mn,
                               fire1yr_cbi_mn, fire2_5yr_cbi_mn,
                               fire6_10yr_cbi_mn, fire11_35yr_cbi_mn,
                               standage_f3_mn, cpycovr_f3_mn) |> 
  filter(Cell_Unit %in% cu.map$Cell_Unit) |> 
  arrange(match(Cell_Unit, cu.map$Cell_Unit))

load(here("./Data/JAGS_Data/MSOM_Ragged_2021_95cut.RData"))
str(all3)
## Z matrix
nsite <- win.data.rag$nsite
nspec <- win.data.rag$nspec
nsamp <- dim(all3)[1]
z <- array(NA, dim = c(nsite, nspec, nsamp))
Jacc <- array(NA, dim = c(nsite, nspec, nsamp))

## For the 99 cut the range of Z is: 1205:(1204 + nspec*nsite)
for(j in 1:nsamp){
  cat(paste("\nMCMC sample", j, "\n"))
  z[,,j] <- all3[j, 1447:(1446 + nspec*nsite)]
}

#saveRDS(z, file = here("./Data/JAGS_Data/Z_Array_Assembled_95thresh.RData"))

# Restrict computations to observed species
zobs <- z[,1:91,] # Species 1 to 145
# Compute Jaccard index for sites and for species
Jsite <- array(NA, dim = c(nsite, nsamp))
Jspec <- array(NA, dim = c(91, nsamp))
Jbetap <- array(NA, dim = c(nsite, nsamp))
# Choose reference site and species for Jaccard indices
ref.site <- 95 # Just choose first site
ref.species <- 1 # Acorn Woodpecker (check object 'obs.occ')
# Get posterior distributions for Jsite and Jspec (for references)
for(k in 1:nsamp){
  for(i in 1:nsite){ # Jaccard index for sites (in terms of shared species)
    Jsite[i,k] <- sum(zobs[ref.site,,k] * zobs[i,,k]) / (sum(zobs[ref.site,,k]) +
                                                           sum(zobs[i,,k]) - sum(zobs[ref.site,,k] * zobs[i,,k]))
  }
  for(i in 1:nspec){ # Jacc. index for species (in terms of shared sites)
    Jspec[i,k] <- sum(zobs[,ref.species,k] * zobs[,i,k]) / (sum(zobs[,ref.species,k]) +
                                                              sum(zobs[,i,k]) - sum(zobs[,ref.species,k] * zobs[,i,k]))
  }
}

for(k in 1:nsamp){
  test <- betapart::beta.pair(zobs[,,k], index.family = "sor")
  sim <- test$beta.sim
}

# NA's arise when a site has no species or a species no sites
# Get posterior means, standard deviations and 95% CRI
# Jaccard index for sites, compared to reference site 1
pm <- apply(Jsite, 1, mean, na.rm = TRUE) # Post. mean of Jsite wrt. site 1
psd <- apply(Jsite, 1, sd, na.rm = TRUE) # Post. sd of Jsite wrt. site 1
cri <- apply(Jsite, 1, function(x) quantile(x, prob = c(0.025, 0.975), na.rm =
                                              TRUE)) # CRI
cbind('post. mean' = pm, 'post. sd' = psd, '2.5%' = cri[1,], '97.5%' = cri[2,])

x <- 3 # size setting for plotting symbol

JZ <- data.frame(Y = aru_meta$Y, X = aru_meta$X,
                 Elevation = aru_meta$topo_elev,
                 Beta = pm, Beta.sd = psd, Beta.cril = cri[1,], Beta.criu = cri[2,])

JZsf <- JZ |> st_as_sf(coords = c("X", "Y"), crs = 4326)

ggplot(JZsf) + 
  geom_sf(aes(size = Beta)) + 
  scale_size_continuous() + 
  theme_bw()

library(mapview)

mapView(JZsf, zcol = "Beta")

## Typical vegan based adonis2
test <- zobs[,,1]
test <- test[-rmsites,]
aru_meta <- as.data.frame(aru_meta)
aru_meta <- aru_meta[-rmsites,]
colnames(test) <- paste0("Species", 1:91)
vegmod <- vegan::adonis2(test ~ utmn + topo_elev + standage_f3_mn + cpycovr_f3_mn + fire6_10yr_cbi_mn + fire11_35yr_cbi_mn, data = aru_meta, method = "jaccard", by = "terms")
vegmod

# Jaccard index for species, compared with a reference species
# (species 13, European Sparrowhawk)
pm <- apply(Jspec, 1, mean, na.rm = TRUE) # Post. mean of Jspec wrt. species 1
psd <- apply(Jspec, 1, sd, na.rm = TRUE) # Post. sd of Jspec wrt. species 1
cri <- apply(Jspec, 1, function(x) quantile(x, prob = c(0.025, 0.975), na.rm =
                                              TRUE)) # CRI
tmp <- cbind('post. mean' = pm, 'post. sd' = psd, '2.5%' = cri[1,], '97.5%' = cri[2,])
rownames(tmp) <- sp.names[[3]]
print(tmp) # print in systematic order
print(tmp[rev(order(tmp[,1])),]) # print in order of decreasing Jacc. values
