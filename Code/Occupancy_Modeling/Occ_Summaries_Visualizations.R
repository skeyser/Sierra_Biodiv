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
load("R:/Users/skeyser/Postdoc/MSOM_Ragged_JAGS_Summaries.Rdata")
load("R:/Users/skeyser/Postdoc/MSOM_Ragged_JAGS_Zout.Rdata")

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
                    "High Sev. Burn: 1 yr post",
                    "High Sev. Burn: 2-5 yr post",
                    "High Sev. Burn: 6-10 yr post",
                    "High Sev. Burn: 11-35 yr post",
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
sp.index <- read.csv(here("Code/Occupancy_Modeling/SpeciesIndex.csv"))
sp.index$Index <- as.character(sp.index$Index)

str(out3)
all3 <- as.matrix(out3)
rm(out3)


pm <- apply(all3, 2, mean)
cri <- apply(all3, 2, function(x) quantile(x, prob = c(0.025, 0.975)))

nspec <- 91
npar <- 13
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
  "Survey Hours" = "Survey Hours",
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
  facet_wrap(~ParPretty, as_labeller(labels_det, default = label_parsed))

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

sp.resp |> 
  filter(str_detect(Par, "beta")) |> 
  mutate(ParPretty = case_when(Par == "beta1" ~ "Latitude",
                               Par == "beta2" ~ "Elevation",
                               Par == "beta3" ~ "Elevation2",
                               Par == "beta4" ~ "High Sev. Burn: 1 yr post",
                               Par == "beta5" ~ "High Sev. Burn: 2-5 yr post",
                               Par == "beta6" ~ "High Sev. Burn: 6-10 yr post",
                               Par == "beta7" ~ "High Sev. Burn: 11-35 yr post",
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
  scale_color_manual(values = c("Above" = "#99CCFF",
                                "Below" = "#FF6600",
                                "Intersect" = "gray")) + 
  coord_flip() + 
  theme_bw() +
  xlab("Species") + 
  ylab("Parameter Estimate") + 
  facet_wrap(~ParPretty, labeller = as_labeller(labels, default = label_parsed))

## Habitat
sp.resp |> 
  filter(str_detect(Par, "beta")) |> 
  mutate(ParPretty = case_when(Par == "beta1" ~ "Latitude",
                               Par == "beta2" ~ "Elevation",
                               Par == "beta3" ~ "Elevation2",
                               Par == "beta4" ~ "High Sev. Burn: 1 yr post",
                               Par == "beta5" ~ "High Sev. Burn: 2-5 yr post",
                               Par == "beta6" ~ "High Sev. Burn: 6-10 yr post",
                               Par == "beta7" ~ "High Sev. Burn: 11-35 yr post",
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
  scale_color_manual(values = c("Above" = "#99CCFF",
                                "Below" = "#FF6600",
                                "Intersect" = "gray")) + 
  coord_flip() + 
  theme_bw() +
  xlab("Species") + 
  ylab("Parameter Estimate") + 
  facet_wrap(~ParPretty, nrow = 1)

## Fire
sp.resp |> 
  filter(str_detect(Par, "beta")) |> 
  mutate(ParPretty = case_when(Par == "beta1" ~ "Latitude",
                               Par == "beta2" ~ "Elevation",
                               Par == "beta3" ~ "Elevation2",
                               Par == "beta4" ~ "High Sev. Burn: 1 yr post",
                               Par == "beta5" ~ "High Sev. Burn: 2-5 yr post",
                               Par == "beta6" ~ "High Sev. Burn: 6-10 yr post",
                               Par == "beta7" ~ "High Sev. Burn: 11-35 yr post",
                               Par == "beta8" ~ "Stand Age",
                               Par == "beta9" ~ "Canopy Cover")) |>
  filter(str_detect(ParPretty, "High Sev.")) |> 
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
  facet_wrap(~ParPretty, nrow = 1)


