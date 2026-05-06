## -------------------------------------------------------------
##
## Script name: spOccupancy Output Processing
##
## Script purpose:
##
## Author: Spencer R Keyser
##
## Date Created: 2025-07-30
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
library(patchwork)
library(here)
library(MCMCvis)
library(spOccupancy)

## -------------------------------------------------------------

logit2prob <- function(logit){
  odds <- exp(logit)
  prob <- odds / (1 + odds)
  return(prob)
}

## -------------------------------------------------------------
##
## Begin Section: Load in model from spOccupancy
##
## -------------------------------------------------------------

## Load spOccupancy model object
out.fo <- readRDS(here("./Data/spOccupancy_Data/Model_sfMsPGOccFO_2slf_R2.rds"))
out <- readRDS(here("./Data/spOccupancy_Data/Model_sfMsPGOcc_2slf_R2.rds"))

## -------------------------------------------------------------
##
## Begin Section: Posterior Predictive Checks
##
## -------------------------------------------------------------
out
summary(out.fo)

# Takes a few minutes to run. 
ppc.sfms.out.ft <- ppcOcc(out, 'freeman-tukey', group = 1)
ppc.sfms.out.rep.ft <- ppcOcc(out, 'freeman-tukey', group = 2)
ppc.sfms.out.cs <- ppcOcc(out, 'chi-squared', group = 1)
ppc.sfms.out.rep.cs <- ppcOcc(out, 'chi-squared', group = 2)

## Summarize the output
summary(ppc.sfms.out.ft)
summary(ppc.sfms.out.rep.ft)
summary(ppc.sfms.out.cs)
summary(ppc.sfms.out.rep.cs)

## -------------------------------------------------------------
##
## End Section: Posterior Predictive Checks
##
## -------------------------------------------------------------

## Model summary - community
summary(out, level = "community")

## Model summary - species
summary(out, level = "species", param = "beta")

## Trace and density plot for the community level coefficients
plot(out, level = "community", param = "beta.comm")

## Old code for JAGS models
# ## Tree plot
# out2$parameters
# 
# ## Community-wide effects
# ## Occupancy model
# com_eff <- MCMCvis::MCMCplot(out,
#                              params = c(paste0("mu.beta", seq(1,9))),
#                              ci = c(50,89),
#                              ref_ovl = TRUE,
#                              labels = c("Latitude", 
#                                         "Elevation", 
#                                         expression("Elevation"^2),
#                                         "Precipitation",
#                                         #"Mean CBI: 1 yr post",
#                                         "Mean CBI: 1-5 yr post",
#                                         "Mean CBI: 6-10 yr post",
#                                         "Mean CBI: 11-35 yr post",
#                                         #"Stand Age",
#                                         "Canopy Cover",
#                                         "Canopy Height"),
#                              rank = T)
# 
# ## Detection model
# MCMCvis::MCMCplot(out2,
#                   params = c(paste0("mu.alpha", seq(1,3))),
#                   ci = c(50,89),
#                   ref_ovl = TRUE,
#                   labels = c("Efforts Hours",
#                              "JDate",
#                              expression("JDate"^2)),
#                   rank = T)
# 
# par(mfrow = c(1,2))
# psi.sample <- plogis(rnorm(10^6, mean = out2$mean$mu.lpsi, sd = out2$mean$sd.lpsi))
# p.sample <- plogis(rnorm(10^6, mean = out2$mean$mu.lp, sd = out2$mean$sd.lp))
# quantile(p.sample, probs = c(0.025, 0.975))
# hist(psi.sample, freq = F, breaks = 50, col = "grey", xlab = "Species occupancy
# probability", ylab = "Density", main = "")
# 
# hist(p.sample, freq = F, breaks = 50, col = "grey", xlab = "Species detection probability",
#      ylab = "Density", main = "")
# summary(psi.sample) ; summary(p.sample)
# 
# # Density plot
# dens <- density(p.sample)
# plot(dens, main="Community Detection Probability Distribution",
#      xlab="Detection Probability", ylab="Density")
# abline(v=c(0.277, 0.67, 0.93), col=c("red", "blue", "red"), lty=c(2,1,2))
# legend("topright", legend=c("95% CI", "Mean"), 
#        col=c("red", "blue"), lty=c(2,1))
# 
# 
# ## Histograms of species richness
# hist(out2$mean$Nsite, main = "Species Richness", breaks = 20)
# mean(out2$mean$Nsite)
# mean(out2$q2.5$Nsite)
# mean(out2$q97.5$Nsite)
# var(out2$mean$Nsite)
# sd(out2$mean$Nsite)


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Tidybayes Prettier Plots
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(tidybayes)
library(tidyr)

# First, gather all mean parameters
mcmc_draws <- out$beta.comm.samples |>
  as.data.frame() |> 
  pivot_longer(cols = everything(),
               names_to = "Parameter",
               values_to = "Draws") |> 
  mutate(Parameter = gsub("scale\\(([^)]*)\\)", "\\1", Parameter)) |> 
  mutate(Parameter = case_when(Parameter == "(Intercept)" ~ "Intercept",
                               Parameter == "I(ele.res^2)" ~ "ele_sq",
                               Parameter == "I(ppt^2)" ~ "ppt_sq",
                               .default = Parameter))
glimpse(mcmc_draws)
unique(mcmc_draws$Parameter)
# Create more informative parameter labels (optional)
mcmc_draws <- mcmc_draws |>
  filter(Parameter != "Intercept") |> 
  mutate(Parameter = factor(Parameter, 
                            levels = unique(Parameter))) |> 
  mutate(Pred_Pretty = case_when(
    #Parameter == "Intercept" ~ "Intercept",
    Parameter == "ele.res" ~ "Elevation",
    Parameter == "ele_sq" ~ "Elevation^2",
    Parameter == "ppt" ~ "Precipitation",
    Parameter == "ppt_sq" ~ "Precipitation^2",
    Parameter == "cbi1_5" ~ "Fire Sev. 1-5yr",
    Parameter == "cbi6_10" ~ "Fire Sev. 6-10yr",
    Parameter == "cbi11_35" ~ "Fire Sev. 11-35yr",
    Parameter == "cc" ~ "Canopy Cover",
    Parameter == "ch.res" ~ "Canopy Hght."
  )) |> 
  mutate(Pred_Cat = case_when(str_detect(Pred_Pretty, "Latitude|Elevation|Precipitation") ~ "Climate/Geographic",
                              #Pred_Pretty == "Precipitation" ~ "Climate",
                              str_detect(Pred_Pretty, "Fire") ~ "Fire",
                              str_detect(Pred_Pretty, "Canopy") ~ "Forest")) |> 
  mutate(Pred_Pretty = factor(Pred_Pretty,
                              levels = rev(c(
                                "Latitude",
                                "Elevation",
                                "Elevation^2",
                                "Precipitation",
                                "Precipitation^2",
                                "Fire Sev. 1-5yr",
                                "Fire Sev. 6-10yr",
                                "Fire Sev. 11-35yr",
                                "Canopy Cover",
                                "Canopy Hght." 
                              )))) |> 
  mutate(Pred_Cat = factor(Pred_Cat,
                           levels = rev(c(
                             "Climate/Geographic",
                             #"Climate",
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
    #"Climate/Geographic" = "#0077BB",         # blue
    "Climate/Geographic" = "#7570b3"        # purple
  )) +
  scale_y_discrete(labels = function(x) {
    x <- gsub("\\^2", "²", x)  # Replace ^2 with ²
    return(x)
  },
  expand = c(0,0.1)) +
  scale_x_continuous(limits = c(-0.8, 0.5),
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
        legend.position = "bottom") + 
  ggtitle("Full MSOM")


ggsave(plot = pretty_eff, filename = here("./Figures/Preliminary/CommunityEff_HalfEye_SpVarThresh_SpOcc.jpg"),
       height = 8, width = 8, dpi = 600)


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Fire Only Model
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# First, gather all mean parameters
mcmc_draws <- out.fo$beta.comm.samples |>
  as.data.frame() |> 
  pivot_longer(cols = everything(),
               names_to = "Parameter",
               values_to = "Draws") |> 
  mutate(Parameter = gsub("scale\\(([^)]*)\\)", "\\1", Parameter)) |> 
  mutate(Parameter = case_when(Parameter == "(Intercept)" ~ "Intercept",
                               Parameter == "I(ele.res^2)" ~ "ele_sq",
                               Parameter == "I(ppt^2)" ~ "ppt_sq",
                               .default = Parameter))
glimpse(mcmc_draws)
unique(mcmc_draws$Parameter)
# Create more informative parameter labels (optional)
mcmc_draws <- mcmc_draws |>
  filter(Parameter != "Intercept") |> 
  mutate(Parameter = factor(Parameter, 
                            levels = unique(Parameter))) |> 
  mutate(Pred_Pretty = case_when(
    #Parameter == "Intercept" ~ "Intercept",
    Parameter == "ele.res" ~ "Elevation",
    Parameter == "ele_sq" ~ "Elevation^2",
    Parameter == "ppt" ~ "Precipitation",
    Parameter == "ppt_sq" ~ "Precipitation^2",
    Parameter == "cbi1_5" ~ "Fire Sev. 1-5yr",
    Parameter == "cbi6_10" ~ "Fire Sev. 6-10yr",
    Parameter == "cbi11_35" ~ "Fire Sev. 11-35yr",
    Parameter == "cc" ~ "Canopy Cover",
    Parameter == "ch.res" ~ "Canopy Hght."
  )) |> 
  mutate(Pred_Cat = case_when(str_detect(Pred_Pretty, "Latitude|Elevation|Precipitation") ~ "Climate/Geographic",
                              #Pred_Pretty == "Precipitation" ~ "Climate",
                              str_detect(Pred_Pretty, "Fire") ~ "Fire",
                              str_detect(Pred_Pretty, "Canopy") ~ "Forest")) |> 
  mutate(Pred_Pretty = factor(Pred_Pretty,
                              levels = rev(c(
                                "Latitude",
                                "Elevation",
                                "Elevation^2",
                                "Precipitation",
                                "Precipitation^2",
                                "Fire Sev. 1-5yr",
                                "Fire Sev. 6-10yr",
                                "Fire Sev. 11-35yr",
                                "Canopy Cover",
                                "Canopy Hght." 
                              )))) |> 
  mutate(Pred_Cat = factor(Pred_Cat,
                           levels = rev(c(
                             "Climate/Geographic",
                             #"Climate",
                             "Fire",
                             "Forest"
                           ))))

# Create the half-eye plot
# Create the half-eye plot with different colors
pretty_eff_fo <- ggplot(mcmc_draws, aes(y = Pred_Pretty, x = Draws, 
                                     fill = Pred_Cat)) +
  stat_halfeye(alpha = 0.7,
               .width = c(0.95, 0.89, 0.5)) +
  geom_vline(xintercept = 0, linetype = "dashed", 
             color = "black", alpha = 0.6, size = 0.8) +
  scale_fill_manual(values = c(
    "Fire" = "#CC3311",           # red
    "Forest" = "#009988", # green
    #"Climate/Geographic" = "#0077BB",         # blue
    "Climate/Geographic" = "#7570b3"        # purple
  )) +
  scale_y_discrete(labels = function(x) {
    x <- gsub("\\^2", "²", x)  # Replace ^2 with ²
    return(x)
  },
  expand = c(0,0.1)) +
  scale_x_continuous(limits = c(-0.8, 0.5),
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
        legend.position = "bottom") + 
  ggtitle("Full MSOM")

pretty_eff_fo


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Combined plot
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mcmc_draws_full <- out$beta.comm.samples |>
  as.data.frame() |> 
  pivot_longer(cols = everything(),
               names_to = "Parameter",
               values_to = "Draws") |> 
  mutate(Parameter = gsub("scale\\(([^)]*)\\)", "\\1", Parameter)) |> 
  mutate(Parameter = case_when(Parameter == "(Intercept)" ~ "Intercept",
                               Parameter == "I(ele.res^2)" ~ "ele_sq",
                               Parameter == "I(ppt^2)" ~ "ppt_sq",
                               .default = Parameter)) |> 
  mutate(model = "Full")

mcmc_draws_fo <- out.fo$beta.comm.samples |>
  as.data.frame() |> 
  pivot_longer(cols = everything(),
               names_to = "Parameter",
               values_to = "Draws") |> 
  mutate(Parameter = gsub("scale\\(([^)]*)\\)", "\\1", Parameter)) |> 
  mutate(Parameter = case_when(Parameter == "(Intercept)" ~ "Intercept",
                               Parameter == "I(ele.res^2)" ~ "ele_sq",
                               Parameter == "I(ppt^2)" ~ "ppt_sq",
                               .default = Parameter)) |> 
  mutate(model = "Fire Only")

## MCMC Draws both
mcmc_draws <- bind_rows(mcmc_draws_full, mcmc_draws_fo)
glimpse(mcmc_draws)
# Create more informative parameter labels (optional)
mcmc_draws <- mcmc_draws |>
  filter(Parameter != "Intercept") |> 
  mutate(Parameter = factor(Parameter, 
                            levels = unique(Parameter))) |> 
  mutate(Pred_Pretty = case_when(
    #Parameter == "Intercept" ~ "Intercept",
    Parameter == "ele.res" ~ "Elevation",
    Parameter == "ele_sq" ~ "Elevation^2",
    Parameter == "ppt" ~ "Precipitation",
    Parameter == "ppt_sq" ~ "Precipitation^2",
    Parameter == "cbi1_5" ~ "Fire Sev. 1-5yr",
    Parameter == "cbi6_10" ~ "Fire Sev. 6-10yr",
    Parameter == "cbi11_35" ~ "Fire Sev. 11-35yr",
    Parameter == "cc" ~ "Canopy Cover",
    Parameter == "ch.res" ~ "Canopy Hght."
  )) |> 
  mutate(Pred_Cat = case_when(str_detect(Pred_Pretty, "Latitude|Elevation|Precipitation") ~ "Climate/Geographic",
                              #Pred_Pretty == "Precipitation" ~ "Climate",
                              str_detect(Pred_Pretty, "Fire") ~ "Fire",
                              str_detect(Pred_Pretty, "Canopy") ~ "Forest")) |> 
  mutate(Pred_Pretty = factor(Pred_Pretty,
                              levels = rev(c(
                                "Latitude",
                                "Elevation",
                                "Elevation^2",
                                "Precipitation",
                                "Precipitation^2",
                                "Fire Sev. 1-5yr",
                                "Fire Sev. 6-10yr",
                                "Fire Sev. 11-35yr",
                                "Canopy Cover",
                                "Canopy Hght." 
                              )))) |> 
  mutate(Pred_Cat = factor(Pred_Cat,
                           levels = rev(c(
                             "Climate/Geographic",
                             #"Climate",
                             "Fire",
                             "Forest"
                           ))))

## Combined
pretty_eff_both <- ggplot(mcmc_draws, aes(y = Pred_Pretty, 
                                          x = Draws, 
                                          #fill = Pred_Cat, 
                                          group = factor(model, 
                                                         levels = c("Full", "Fire Only")),
                                          #linetype = model,
                                          fill = model)) +
  stat_halfeye(
    data = ~ filter(.x, model == "Full"),
    alpha = 0.7,
    .width = c(0.95, 0.89, 0.5),
    side = "top",
    position = position_nudge(y = 0.05),
    scale = 0.5) +
  stat_halfeye(
    data = ~ filter(.x, model == "Fire Only"),
    alpha = 0.7,
    .width = c(0.95, 0.89, 0.5),
    side = "bottom",
    position = position_nudge(y = -0.05),
    scale = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", 
             color = "black", alpha = 0.6, size = 0.8) +
  # scale_fill_manual(values = c(
  #   "Fire" = "#CC3311",           # red
  #   "Forest" = "#009988", # green
  #   #"Climate/Geographic" = "#0077BB",         # blue
  #   "Climate/Geographic" = "#7570b3"        # purple
  # )) +
  scale_fill_manual(values = c(
    "Fire Only" = "#CC3311", # red
    "Full" = "#009988" # green
  ),
  labels = c("Fire Only", "Full")) +
  scale_y_discrete(labels = function(x) {
    x <- gsub("\\^2", "²", x)  # Replace ^2 with ²
    return(x)
  },
  expand = c(0,0.1)) +
  # scale_linetype_manual(values = c("Full" = "solid", "Fire Only" = "dashed"),
  #                       name = "Model Type") +
  scale_x_continuous(limits = c(-0.8, 0.5),
                     labels = scales::number_format(accuracy = 0.1, drop0trailing = T)) +
  theme_minimal() +
  labs(x = "Effect Size (Assemblage)",
       y = "Predictor",
       fill = NULL) +
  # theme(axis.text.y = element_text(hjust = 0, family = "sans", size = 12),
  #       axis.text.x = element_text(family = "sans", size = 12),
  #       axis.title = element_text(family = "sans", size = 12),
  #       plot.title = element_text(hjust = 0.5),
  #       legend.title = element_blank(),
  #       legend.text = element_text(family = "sans", size = 12),
  #       legend.position = "bottom") + 
  theme(strip.background = element_blank(),
        strip.text = element_text(family = "sans", size = 12),
        axis.text.y = element_text(family = "sans", size = 9),
        axis.text.x = element_text(family = "sans", size = 12),
        axis.title = element_text(family = "sans", size = 12),
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom")  

pretty_eff_both

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Species-level effects plots
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sp.resp <- out$beta.samples |> 
  as.data.frame() |> 
  pivot_longer(cols = everything(),
               names_to = "Parameter",
               values_to = "Draws") |> 
  separate(col = Parameter, into = c("Parameter", "Species"), sep = "(?:-)", extra = "merge") |> 
  mutate(Parameter = gsub("scale\\(([^)]*)\\)", "\\1", Parameter)) |> 
  mutate(Parameter = case_when(Parameter == "(Intercept)" ~ "Intercept",
                               Parameter == "I(ele.res^2)" ~ "ele_sq",
                               .default = Parameter)) |> 
  group_by(Parameter, Species) |> 
  summarise(Mean = mean(Draws),
            UCI = quantile(Draws, probs = 0.975),
            LCI = quantile(Draws, probs = 0.025)) |> 
  mutate(model = "full")

sp.resp.fo <- out.fo$beta.samples |> 
  as.data.frame() |> 
  pivot_longer(cols = everything(),
               names_to = "Parameter",
               values_to = "Draws") |> 
  separate(col = Parameter, into = c("Parameter", "Species"), sep = "(?:-)", extra = "merge") |> 
  mutate(Parameter = gsub("scale\\(([^)]*)\\)", "\\1", Parameter)) |> 
  mutate(Parameter = case_when(Parameter == "(Intercept)" ~ "Intercept",
                               Parameter == "I(ele.res^2)" ~ "ele_sq",
                               .default = Parameter)) |> 
  group_by(Parameter, Species) |> 
  summarise(Mean = mean(Draws),
            UCI = quantile(Draws, probs = 0.975),
            LCI = quantile(Draws, probs = 0.025)) |> 
  mutate(model = "fire only")

glimpse(sp.resp)
glimpse(sp.resp.fo)

sp.resp <- bind_rows(sp.resp, sp.resp.fo)

# pm <- apply(all3, 2, mean)
# cri <- apply(all3, 2, function(x) quantile(x, prob = c(0.025, 0.975)))
# 
# nspec <- nrow(sp.index)
# npar <- 14
# N <- nspec*npar
# sp.resp <- data.frame(Par = names(pm[1:N]),
#                       Mean = pm[1:N],
#                       UCI = cri[2,1:N],
#                       LCI = cri[1,1:N])
# sp.resp <- sp.resp |> 
#   mutate(Species = str_extract(Par, "\\[\\d+\\]")) |> 
#   mutate(Species = gsub("[[:punct:]]", "", Species)) |> 
#   left_join(sp.index, by = c("Species" = "Index")) |> 
#   select(Par, Mean, UCI, LCI, Species = Species.y) |> 
#   mutate(Par = gsub("\\[\\d+\\]", "", Par))


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
  filter(Parameter == "lp") |> 
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
  ylab("Parameter Estimate") +
  facet_wrap(~model)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Occupancy Covariates
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Overall occupancy estimates by species
sp.resp |> 
  filter(Parameter == "Intercept") |> 
  mutate(Mean = plogis(Mean), UCI = plogis(UCI), LCI = plogis(LCI)) |> 
  ggplot(aes(y = Mean, x = reorder(Species, -Mean), color = model, group = model)) +
  #geom_hline(yintercept = mean(Mean), linetype = "dashed", color = "black", size = 1.2) +
  geom_pointrange(aes(ymin = LCI, ymax = UCI), 
                  size = 0.8,
                  show.legend = T) + 
  coord_flip() + 
  theme_bw() +
  xlab("Species") + 
  ylab(expression(psi))

## Geo/Topo
labels <- c(
  "Elevation" = "Elevation",
  "Elevation^2" = "Elevation^2",
  "Precipitation" = "Precip",
  "Precipitation^2" = "Precip^2"
)

geo.eff.plot <- sp.resp |> 
  ungroup() |> 
  filter(Parameter != "Intercept") |> 
  mutate(Pred_Pretty = case_when(
    #Parameter == "Intercept" ~ "Intercept",
    Parameter == "ele.res" ~ "Elevation",
    Parameter == "ele_sq" ~ "Elevation^2",
    Parameter == "ppt" ~ "Precipitation",
    Parameter == "I(ppt^2)" ~ "Precipitation^2",
    Parameter == "cbi1_5" ~ "Fire Sev. 1-5yr",
    Parameter == "cbi6_10" ~ "Fire Sev. 6-10yr",
    Parameter == "cbi11_35" ~ "Fire Sev. 11-35yr",
    Parameter == "cc_cfo" ~ "Canopy Cover",
    Parameter == "ch_res" ~ "Canopy Hght."
  )) |> 
  filter(Pred_Pretty %in% c("Elevation", "Elevation^2", "Precipitation", "Precipitation^2")) |> 
  ggplot(aes(y = Mean, x = reorder(Species, -Mean), 
             shape = factor(model), group = factor(model))) +
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
  scale_shape_manual(values = c(16, 21)) +
  coord_flip() + 
  theme_bw() +
  xlab("Species") + 
  ylab("Parameter Estimate") + 
  facet_wrap(~Pred_Pretty, labeller = as_labeller(labels, default = label_parsed))

ggsave(plot = geo.eff.plot, filename = here("./Figures/Preliminary/MSOM_EffPlot_GeoCovariates_SpVarThresh.jpg"),
       height = 10, width = 10, dpi = 600)

## Habitat
hab.eff.plot <- sp.resp |> 
  ungroup() |> 
  filter(Parameter != "Intercept") |> 
  mutate(Pred_Pretty = case_when(
    #Parameter == "Intercept" ~ "Intercept",
    Parameter == "ele" ~ "Elevation",
    Parameter == "ele_sq" ~ "Elevation^2",
    Parameter == "ppt" ~ "Precipitation",
    Parameter == "cbi1_5" ~ "Fire Sev. 1-5yr",
    Parameter == "cbi6_10" ~ "Fire Sev. 6-10yr",
    Parameter == "cbi11_35" ~ "Fire Sev. 11-35yr",
    Parameter == "cc_cfo" ~ "Canopy Cover",
    Parameter == "ch_res" ~ "Canopy Hght."
  )) |> 
  filter(Pred_Pretty %in% c("Canopy Cover", "Canopy Hght.")) |> 
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
  facet_wrap(~Pred_Pretty, nrow = 1)

ggsave(plot = hab.eff.plot, filename = here("./Figures/Preliminary/MSOM_EffPlot_HabCovariates_SpVarThresh.jpg"),
       height = 10, width = 10, dpi = 600)

## All besides fire for SI
## Geo/Topo
labels <- c(
  "Elevation" = "Elevation",
  "Elevation^2" = "Elevation^2",
  "Precipitation" = "Precip",
  "Canopy Cover" = "'Canopy Cover'",
  "Canopy Hght." = "'Canopy Hght.'"
)

geohab.eff.plot <- sp.resp |> 
  ungroup() |> 
  filter(Parameter != "Intercept") |> 
  mutate(Pred_Pretty = case_when(
    #Parameter == "Intercept" ~ "Intercept",
    Parameter == "ele" ~ "Elevation",
    Parameter == "ele_sq" ~ "Elevation^2",
    Parameter == "ppt" ~ "Precipitation",
    Parameter == "cbi1_5" ~ "Fire Sev. 1-5yr",
    Parameter == "cbi6_10" ~ "Fire Sev. 6-10yr",
    Parameter == "cbi11_35" ~ "Fire Sev. 11-35yr",
    Parameter == "cc_cfo" ~ "Canopy Cover",
    Parameter == "ch_res" ~ "Canopy Hght."
  )) |> 
  filter(Pred_Pretty %in% c("Elevation", "Elevation^2", "Precipitation", "Canopy Cover", "Canopy Hght.")) |> 
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
  facet_wrap(~Pred_Pretty, nrow = 1, labeller = as_labeller(labels, default = label_parsed))

ggsave(plot = geohab.eff.plot, filename = here("./Figures/Preliminary/MSOM_EffPlot_nonFireCovariates_SpVarThresh.jpg"),
       height = 10, width = 10, dpi = 600)


## Fire
fire.eff.plot <- sp.resp |> 
  ungroup() |> 
  filter(Parameter != "Intercept") |> 
  mutate(Pred_Pretty = case_when(
    #Parameter == "Intercept" ~ "Intercept",
    Parameter == "ele" ~ "Elevation",
    Parameter == "ele_sq" ~ "Elevation^2",
    Parameter == "ppt" ~ "Precipitation",
    Parameter == "cbi1_5" ~ "Fire Sev. 1-5yr",
    Parameter == "cbi6_10" ~ "Fire Sev. 6-10yr",
    Parameter == "cbi11_35" ~ "Fire Sev. 11-35yr",
    Parameter == "cc_cfo" ~ "Canopy Cover",
    Parameter == "ch_res" ~ "Canopy Hght."
  )) |> 
  filter(str_detect(Pred_Pretty, "Fire")) |> 
  mutate(Pred_Pretty = factor(Pred_Pretty,
                              levels = c(
                                "Fire Sev. 1-5yr",
                                "Fire Sev. 6-10yr",
                                "Fire Sev. 11-35yr"))) |> 
  ggplot(aes(y = Mean, x = reorder(Species, -Mean), 
             shape = factor(model))) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 1.2) +
  geom_pointrange(aes(ymin = LCI, ymax = UCI,
                      color = case_when(
                        LCI > 0 ~ "Above",
                        UCI < 0 ~ "Below",
                        TRUE ~ "Intersects")), 
                  size = 0.5,
                  show.legend = T,
                  position = position_dodge(width = 0.7)) +
  scale_color_manual(values = c("Above" = "#FF6600",
                                "Below" = "#99CCFF",
                                "Intersect" = "gray"),
                     guide = "none") + 
  scale_shape_manual(values = c("full" = 16, 
                                "fire only" = 17),
                     labels = c("Fire Only", "Full")) +
  scale_y_continuous(labels = scales::number_format(drop0trailing = T)) +
  coord_flip() + 
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text = element_text(family = "sans", size = 12),
        axis.text.y = element_text(family = "sans", size = 9),
        axis.text.x = element_text(family = "sans", size = 12),
        axis.title = element_text(family = "sans", size = 12),
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom") + 
  xlab("Species") + 
  ylab("Effect Size (Species)") + 
  labs(shape = NULL) +
  facet_wrap(~Pred_Pretty, nrow = 1)

ggsave(plot = fire.eff.plot, filename = here("./Figures/Final/MSOM_EffPlot_FireCovariates_SpVarThresh_FireOnlyFull.jpg"),
       height = 10, width = 12, dpi = 600)


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Create Fig 2
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Fire plus the overall effects
fig2 <- pretty_eff_both + free(fire.eff.plot) + plot_layout(widths = c(0.5, 1.5)) + plot_annotation(tag_levels = 'A')
ggsave(plot = fig2, filename = here("./Figures/Final/Fig2_spOcc_sfm2_BothModels.jpg"),
       height = 10, width = 10, dpi = 800)

ggsave(plot = fire.eff.plot, filename = here("./Figures/Final/Fig2_spOcc_FireOnly.jpg"),
       height = 10, width = 8, dpi = 800)


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Stats on species responses to fire
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Significant response to either 
sp.resp |> 
  filter(Parameter %in% c("cbi1_5", "cbi6_10", "cbi11_35")) |> 
  mutate(Sig = ifelse(UCI < 0 | LCI > 0, "Sig", "NonSig")) |>
  filter(Sig == "Sig") |> 
  mutate(Dir = ifelse(UCI < 0, "Neg", "Pos")) #|>
  group_by(Parameter) |> 
  count(Sig)

sp.resp |> 
  filter(Parameter %in% c("cbi1_5", "cbi6_10", "cbi11_35")) |>
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