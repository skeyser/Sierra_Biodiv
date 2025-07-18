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
library(jagsUI)

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
#load("R:/Users/skeyser/Postdoc/MSOM_Ragged_JAGS_Summaries_95thresh.Rdata")
#load("R:/Users/skeyser/Postdoc/MSOM_Ragged_JAGS_Zout_95thresh.Rdata")

load(here("./Data/JAGS_Output/MSOM_Ragged_JAGS_Summaries_VarThresh_975min.Rdata"))
load(here("./Data/JAGS_Output/MSOM_Ragged_JAGS_Zout_VarThresh_975min.Rdata"))

## -------------------------------------------------------------
##
## Begin Section: Posterior Predictive Checks
##
## -------------------------------------------------------------
out2
summary(out2)
pp.check(out2, observed = 'fitY', simulated = 'fitY.new')
#pp.check(out2, observed = 'fitZ', simulated = 'fitZ.new')

## -------------------------------------------------------------
##
## End Section: Posterior Predictive Checks
##
## -------------------------------------------------------------


## Model summary
MCMCvis::MCMCsummary(out2, round = 2)



MCMCvis::MCMCtrace(out2, params = c("mu.beta1"), n.eff = T, pdf = F)
MCMCvis::MCMCtrace(out2, params = c("mu.beta2"), n.eff = T, pdf = F)

## Tree plot
out2$parameters

## Community-wide effects
## Occupancy model
com_eff <- MCMCvis::MCMCplot(out2,
         params = c(paste0("mu.beta", seq(1,9))),
         ci = c(50,89),
         ref_ovl = TRUE,
         labels = c("Latitude", 
                    "Elevation", 
                    expression("Elevation"^2),
                    "Precipitation",
                    #"Mean CBI: 1 yr post",
                    "Mean CBI: 1-5 yr post",
                    "Mean CBI: 6-10 yr post",
                    "Mean CBI: 11-35 yr post",
                    #"Stand Age",
                    "Canopy Cover",
                    "Canopy Height"),
         rank = T)

## Detection model
MCMCvis::MCMCplot(out2,
                  params = c(paste0("mu.alpha", seq(1,3))),
                  ci = c(50,89),
                  ref_ovl = TRUE,
                  labels = c("Efforts Hours",
                             "JDate",
                             expression("JDate"^2)),
                  rank = T)

par(mfrow = c(1,2))
psi.sample <- plogis(rnorm(10^6, mean = out2$mean$mu.lpsi, sd = out2$mean$sd.lpsi))
p.sample <- plogis(rnorm(10^6, mean = out2$mean$mu.lp, sd = out2$mean$sd.lp))
quantile(p.sample, probs = c(0.025, 0.975))
hist(psi.sample, freq = F, breaks = 50, col = "grey", xlab = "Species occupancy
probability", ylab = "Density", main = "")

hist(p.sample, freq = F, breaks = 50, col = "grey", xlab = "Species detection probability",
     ylab = "Density", main = "")
summary(psi.sample) ; summary(p.sample)

# Density plot
dens <- density(p.sample)
plot(dens, main="Community Detection Probability Distribution",
     xlab="Detection Probability", ylab="Density")
abline(v=c(0.277, 0.67, 0.93), col=c("red", "blue", "red"), lty=c(2,1,2))
legend("topright", legend=c("95% CI", "Mean"), 
       col=c("red", "blue"), lty=c(2,1))


## Histograms of species richness
hist(out2$mean$Nsite, main = "Species Richness", breaks = 20)
mean(out2$mean$Nsite)
mean(out2$q2.5$Nsite)
mean(out2$q97.5$Nsite)
var(out2$mean$Nsite)
sd(out2$mean$Nsite)


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Tidybayes Prettier Plots
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(tidybayes)
library(tidyr)
# First, gather all mean parameters
mcmc_draws <- spread_draws(out2, 
                           mu.beta1, 
                           mu.beta2, 
                           mu.beta3, 
                           mu.beta4, 
                           mu.beta5, 
                           mu.beta6, 
                           mu.beta7, 
                           mu.beta8, 
                           mu.beta9) |>
  pivot_longer(cols = contains("mu"),
               names_to = "Parameter",
               values_to = "Draws")

# Create more informative parameter labels (optional)
mcmc_draws <- mcmc_draws |>
  mutate(Parameter = factor(Parameter, 
                            levels = unique(Parameter),
                            labels = str_remove(unique(Parameter), "mu\\."))) |> 
  mutate(Pred_Pretty = case_when(
    Parameter == "beta1" ~ "Latitude",
    Parameter == "beta2" ~ "Elevation",
    Parameter == "beta3" ~ "Elevation^2",
    Parameter == "beta4" ~ "Precipitation",
    Parameter == "beta5" ~ "Fire Sev. 1-5yr",
    Parameter == "beta6" ~ "Fire Sev. 6-10yr",
    Parameter == "beta7" ~ "Fire Sev. 11-35yr",
    Parameter == "beta8" ~ "Canopy Cover",
    Parameter == "beta9" ~ "Canopy Hght."
  )) |> 
  mutate(Pred_Cat = case_when(str_detect(Pred_Pretty, "Latitude|Elevation") ~ "Geographic",
                              Pred_Pretty == "Precipitation" ~ "Climate",
                              str_detect(Pred_Pretty, "Fire") ~ "Fire",
                              str_detect(Pred_Pretty, "Canopy") ~ "Forest")) |> 
  mutate(Pred_Pretty = factor(Pred_Pretty,
                              levels = c(
                                "Latitude",
                                "Elevation",
                                "Elevation^2",
                                "Precipitation",
                                "Fire Sev. 1-5yr",
                                "Fire Sev. 6-10yr",
                                "Fire Sev. 11-35yr",
                                "Canopy Cover",
                                "Canopy Hght." 
                              ))) |> 
  mutate(Pred_Cat = factor(Pred_Cat,
                           levels = rev(c(
                             "Geographic",
                             "Climate",
                             "Fire",
                             "Forest"
                           ))))

# Create the half-eye plot
# Create the half-eye plot with different colors
pretty_eff <- ggplot(mcmc_draws, aes(y = Pred_Pretty, x = Draws, 
                       fill = Pred_Cat)) +
  stat_halfeye(alpha = 0.7,
               .width = c(0.95, 0.89, 0.5)) +
  geom_vline(xintercept = 0, linetype = "dashed", 
             color = "black", alpha = 0.6, size = 0.8) +
  scale_fill_manual(values = c(
    "Fire" = "#CC3311",           # red
    "Forest" = "#009988", # green
    "Climate" = "#0077BB",         # blue
    "Geographic" = "#7570b3"        # purple
  )) +
  scale_y_discrete(labels = function(x) {
    x <- gsub("\\^2", "²", x)  # Replace ^2 with ²
    return(x)
  },
  expand = c(0,0.1)) +
  scale_x_continuous(limits = c(-0.7, 0.2),
                     labels = scales::number_format(accuracy = 0.1, drop0trailing = T)) +
  theme_minimal() +
  labs(x = "Effect Size (Assemblage)",
       y = "") +
  theme(axis.text.y = element_text(hjust = 0, family = "sans", size = 12),
        axis.text.x = element_text(family = "sans", size = 12),
        axis.title = element_text(family = "sans", size = 12),
        plot.title = element_text(hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(family = "sans", size = 12),
        legend.position = "bottom")


ggsave(plot = pretty_eff, filename = here("./Figures/Preliminary/CommunityEff_HalfEye_SpVarThresh.jpg"),
       height = 8, width = 8, dpi = 600)


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Species-level effects plots
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sp.index <- read.csv(here("Code/Occupancy_Modeling/SpeciesIndex_Filtered_VarThresh.csv"))
sp.index$Index <- as.character(sp.index$Index)

str(out3)
all3 <- as.matrix(out3)
rm(out3)
gc()


pm <- apply(all3, 2, mean)
cri <- apply(all3, 2, function(x) quantile(x, prob = c(0.025, 0.975)))

nspec <- nrow(sp.index)
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
  "Elevation2" = "Elevation^2",
  "Precipitation" = "Precip"
)

geo.eff.plot <- sp.resp |> 
  filter(str_detect(Par, "beta")) |> 
  mutate(ParPretty = case_when(Par == "beta1" ~ "Latitude",
                               Par == "beta2" ~ "Elevation",
                               Par == "beta3" ~ "Elevation2",
                               Par == "beta4" ~ "Precipitation",
                               Par == "beta5" ~ "Mean CBI: 1-5 yr post",
                               Par == "beta6" ~ "Mean CBI: 6-10 yr post",
                               Par == "beta7" ~ "Mean CBI: 11-35 yr post",
                               Par == "beta8" ~ "Canopy Cover",
                               Par == "beta9" ~ "Canopy Height Res.")) |>
  filter(ParPretty %in% c("Latitude", "Elevation", "Elevation2", "Precipitation")) |> 
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

ggsave(plot = geo.eff.plot, filename = here("./Figures/Preliminary/MSOM_EffPlot_GeoCovariates_SpVarThresh.jpg"),
       height = 10, width = 10, dpi = 600)

## Habitat
hab.eff.plot <- sp.resp |> 
  filter(str_detect(Par, "beta")) |> 
  mutate(ParPretty = case_when(Par == "beta1" ~ "Latitude",
                               Par == "beta2" ~ "Elevation",
                               Par == "beta3" ~ "Elevation2",
                               Par == "beta4" ~ "Precipitation",
                               Par == "beta5" ~ "Mean CBI: 1-5 yr post",
                               Par == "beta6" ~ "Mean CBI: 6-10 yr post",
                               Par == "beta7" ~ "Mean CBI: 11-35 yr post",
                               Par == "beta8" ~ "Canopy Cover",
                               Par == "beta9" ~ "Canopy Height Res.")) |>
  filter(ParPretty %in% c("Canopy Cover", "Canopy Height Res.")) |> 
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

ggsave(plot = hab.eff.plot, filename = here("./Figures/Preliminary/MSOM_EffPlot_HabCovariates_SpVarThresh.jpg"),
       height = 10, width = 10, dpi = 600)

## Fire
fire.eff.plot <- sp.resp |> 
  filter(str_detect(Par, "beta")) |> 
  mutate(ParPretty = case_when(Par == "beta1" ~ "Latitude",
                               Par == "beta2" ~ "Elevation",
                               Par == "beta3" ~ "Elevation2",
                               Par == "beta4" ~ "Precipitation",
                               Par == "beta5" ~ "Fire 1-5yr",
                               Par == "beta6" ~ "Fire 6-10yr",
                               Par == "beta7" ~ "Fire 11-35yr",
                               Par == "beta8" ~ "Canopy Cover",
                               Par == "beta9" ~ "Canopy Height Res.")) |>
  filter(str_detect(ParPretty, "Fire")) |> 
  mutate(ParPretty = factor(ParPretty,
                               levels = c(
                               "Fire 1-5yr",
                               "Fire 6-10yr",
                               "Fire 11-35yr"))) |> 
  ggplot(aes(y = Mean, x = reorder(Species, -Mean))) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 1.2) +
  geom_pointrange(aes(ymin = LCI, ymax = UCI,
                      color = case_when(
                        LCI > 0 ~ "Above",
                        UCI < 0 ~ "Below",
                        TRUE ~ "Intersects")), 
                  size = 0.5,
                  show.legend = F) +
  scale_color_manual(values = c("Above" = "#FF6600",
                                "Below" = "#99CCFF",
                                "Intersect" = "gray")) + 
  scale_y_continuous(labels = scales::number_format(drop0trailing = T)) +
  coord_flip() + 
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text = element_text(family = "sans", size = 12),
        axis.text.y = element_text(family = "sans", size = 9),
        axis.text.x = element_text(family = "sans", size = 12),
        axis.title = element_text(family = "sans", size = 12),
        plot.title = element_text(hjust = 0.5)) + 
  xlab("Species") + 
  ylab("Effect Size (Species)") + 
  facet_wrap(~ParPretty, nrow = 1)

ggsave(plot = fire.eff.plot, filename = here("./Figures/Preliminary/MSOM_EffPlot_FireCovariates_SpVarThresh.jpg"),
       height = 10, width = 12, dpi = 600)


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Create Fig 2
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Fire plus the overall effects
fig2 <- pretty_eff + free(fire.eff.plot) + plot_layout(widths = c(0.5, 1)) + plot_annotation(tag_levels = 'A')
ggsave(plot = fig2, filename = here("./Figures/Final/Fig2.jpg"),
       height = 10, width = 10, dpi = 800)


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Stats on species responses to fire
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Significant response to either 
sp.resp |> 
  filter(Par %in% c("beta5", "beta6")) |> 
  mutate(Sig = ifelse(UCI < 0 | LCI > 0, "Sig", "NonSig")) |> 
  group_by(Par) |> 
  count(Sig)

sp.resp |> 
  filter(Par %in% c("beta5", "beta6")) |> 
  mutate(Sig = ifelse(UCI < 0 | LCI > 0, "Sig", "NonSig")) |> 
  filter(Sig == "Sig") |>
  distinct(Species) |> 
  tally()


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: # of sites with atleast 50% high severity fire
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
load(here("./Data/JAGS_Data/Occ2GDM_Data_SpThresh_975minMaxPrec.Rdata"))

## Take the environmental data
aru_meta <- OccData$SiteMeta
aru_meta <- aru_meta |> 
  select(-utme, -utmn) |> 
  st_as_sf(coords = c("X", "Y"), crs = 4326) |> 
  st_transform(crs = 3310) |> 
  mutate(Lat = st_coordinates(geometry)[,2],
         Long = st_coordinates(geometry)[,1]) |> 
  st_drop_geometry()

## Prop Fire
aru_meta |>
  select(Cell_Unit, fire1_5yr_high_prop, fire6_10yr_high_prop) |> 
  filter(fire1_5yr_high_prop > 0.5 | fire6_10yr_high_prop > 0.5) |> 
  distinct(Cell_Unit) |> 
  tally()
