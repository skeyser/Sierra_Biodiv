## -------------------------------------------------------------
##
## Script name: LCDB Sierra
##
## Script purpose:
##
## Author: Spencer R Keyser
##
## Date Created: 2025-05-29
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
options(scipen = 10, digits = 10)

## -------------------------------------------------------------

## Package Loading
library(dplyr)
library(ggplot2)
library(here)
library(adespatial)
library(adegraphics)
library(ade4)
library(performance)
library(DHARMa)
library(sf)
## -------------------------------------------------------------

## Load in the Z-matrices
load(here("./Data/SpOccupancy_Data/Occ2GDM_Data_SpThresh_975minMaxPrec_spOcc_2SLF.Rdata"))

## Reformat data
Z <- OccData$Z.posterior
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

# First create the dataset with FireP
aru_meta <- aru_meta |> 
  select(all_of(v.keep)) |>
  mutate(FireP = ifelse(fire1_5yr_cbi_mn > 0 | fire6_10yr_cbi_mn > 0 | fire11_35yr_cbi_mn > 0, "F", "NF"))

# Get number of F cases
n_fire <- sum(aru_meta$FireP == "F")
hist(aru_meta$fire1_5yr_cbi_mn)

## -------------------------------------------------------------
##
## Begin Section: Format for LCBD esimtation
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
# samp <- sort(ceiling(runif(100, 0, 1000)))
# Z <- Z[,,samp]
# npost <- dim(Z)[3]

## Calculate the LCBD values for each posterior draw
results <- list(
  BetaComp = list(),
  LCBD.D = list(),
  LCBD.Repl = list(),
  LCBD.Rich = list()
)

bdiv_list <- apply(Z, MARGIN = 3, function(x){
  #beta.div(x, "hellinger", nperm = 999)
  b.comp <- beta.div.comp(x, coef = "J", quant = F)
  out.d <- LCBD.comp(b.comp$D, sqrt.D = T)
  out.Repl <- LCBD.comp(b.comp$repl, sqrt.D = T)
  out.Diff <- LCBD.comp(b.comp$rich, sqrt.D = T)
  
  list(
    D = list(
      LCBD = out.d$LCBD
    ),
    Repl = list(
      LCBD = out.Repl$LCBD
    ),
    Rich = list(
      LCBD = out.Diff$LCBD
    )
  )
  
})

## Remove the OccData object prior to loop
rm(OccData)
gc()

## Pull the LCBD Values from the list
lcbd <- sapply(bdiv_list, function(x) x$D$LCBD)
lcbd.mean <- apply(lcbd, 1, mean)
lcbd.sd <- apply(lcbd, 1, sd)

## Plot species contribution to beta
scbd <- as.data.frame(bdiv$SCBD)
scbd$Species <- rownames(scbd)
colnames(scbd)[colnames(scbd) == "bdiv$SCBD"] <- "SCBD"

scbd |> 
  as_tibble() |> 
  arrange(desc(SCBD)) |> 
  tibble::remove_rownames() |> 
  mutate(Species = factor(Species, levels = Species)) |>  # This creates the ordered factor
  ggplot(aes(x = Species, y = SCBD)) +
  coord_flip() + 
  geom_point() + 
  theme_bw() + 
  xlab("Species") + 
  ylab("SCBD")


s.value(envTab[,c(3,2)], bdiv$LCBD, symbol = "circle", col = c("white", "brown"), main="Map of Sierra LCBD")
sig.lcbd <- which(bdiv$p.LCBD <= 0.05)
nonsig.lcbd <- which(bdiv$p.LCBD > 0.05)

g1 <- s.value(envTab[sig.lcbd,c(3,2)], bdiv$LCBD[sig.lcbd], 
        ppoint.alpha = 1,
        plegend.drawKey = FALSE,
        symbol = "circle",
        col = c("white", "brown"), 
        main="Map of Sierra LCBD")

g2 <- s.value(envTab[nonsig.lcbd,c(3,2)], bdiv$LCBD[nonsig.lcbd], 
              ppoint.alpha = 0.1,
              plegend.drawKey = FALSE,
              symbol = "circle",
              col = c("white", "blue"))

g1 + g2

## Extract the LCBD
lcbd <- bdiv$LCBD

## Pair with env
envTab$LCBD <- lcbd
hist(envTab$LCBD)

## Consider taking a random sample now
envTab <- envTab |> 
  mutate(FireP = case_when(
    fire1_5yr_cbi_mn == 0 & fire6_10yr_cbi_mn == 0 & fire11_35yr_cbi_mn == 0 ~ "Unburned",
    fire1_5yr_cbi_mn > 0 & fire1_5yr_cbi_mn < 2.25 ~ "LM Fire 1-5",
    fire1_5yr_cbi_mn >= 2.25 ~ "H Fire 1-5",
    fire6_10yr_cbi_mn > 0 & fire6_10yr_cbi_mn < 2.25 ~ "LM Fire 6-10",
    fire6_10yr_cbi_mn >= 2.25 ~ "H Fire 6-10",
    fire11_35yr_cbi_mn > 0 & fire11_35yr_cbi_mn < 2.25 ~ "LM Fire 11-35",
    fire11_35yr_cbi_mn >= 2.25 ~ "H Fire 11-35"
  )) |> 
  mutate(FireP = factor(FireP, 
                        levels = c("Unburned",
                                   "LM Fire 1-5",
                                   "LM Fire 6-10",
                                   "LM Fire 11-35",
                                   "H Fire 1-5",
                                   "H Fire 6-10",
                                   "H Fire 11-35")))


## -------------------------------------------------------------
##
## Begin Section: BRMS Modelling
##
## -------------------------------------------------------------
dat <- envTab
dat$LCBD <- lcbd.mean
dat$LCBDsd <- lcbd.sd
summary(dat)

## Transform value for Gamma
dat$LCBD_scaled <- dat$LCBD * nrow(dat)
dat$LCBDsd_scaled <- dat$LCBDsd * nrow(dat)
hist(dat$LCBD_scaled)
hist(dat$LCBDsd_scaled)
sum(dat$LCBD_scaled)

## Scale the predictors of interest
dat <- dat |> 
  mutate(across(c(Lat:ch_res),
                .fns = scale,
                .names = "{.col}_scaled"))

## fitting the model
levels(dat$FireP)
m2.me <- brm(
  LCBD_scaled ~ ele_scaled + I(ele_scaled^2) + Lat_scaled + ppt_scaled + cancov_scaled + ch_res_scaled + fire1_5yr_cbi_mn_scaled + fire6_10yr_cbi_mn_scaled + fire11_35yr_cbi_mn_scaled,
  data = dat,
  family = Gamma(link = "log"),
  control = list(adapt_delta = 0.95),
  cores = 3
)

m1.me <- brm(
  LCBD_scaled | mi(LCBDsd_scaled) ~ ele_scaled + I(ele_scaled^2) + Lat_scaled + ppt_scaled + cancov_scaled + ch_res_scaled + FireP,
  data = dat,
  family = Gamma(link = "log"),
  control = list(adapt_delta = 0.95),
  cores = 3
)

summary(m2.me)
summary(m1.me)

compare_models <- loo_compare(
  loo(m2.me),
  loo(m1.me)
)

## Check simulated values
par(mfrow = c(2,1))
pp_check(m2.me, ndraws = 1000)
pp_check(m1.me, ndraws = 1000)

## Plot the effects
ce = conditional_effects(m1.me)
ce$FireP

## Prediction df
pred.dat <- tibble(FireP = unique(dat$FireP),
                   ele_scaled = mean(dat$ele_scaled))

## Posterior draws
post.epred <- epred_draws(m1.me, newdata = pred.dat)

## Condtional effects plots
pce = plot(ce, ask=FALSE, plot=FALSE)

pce[["FireP"]] + 
  stat_halfeye(data = post.epred,
               inherit.aes = FALSE,
               aes(y = FireP, x = .epred),
               .width = c(0.95, 0.66),
               point_interval = "median_qi",
               fill = "red",
               alpha = 0.3)

## Get posterior draws for the categorical effect
# Extract posterior samples for fire categories
fire_effects <- m1.me |> 
  gather_draws(b_FirePLMFire1M5, b_FirePLMFire6M10, b_FirePLMFire11M35,
               b_FirePHFire1M5, b_FirePHFire6M10, b_FirePHFire11M35,
               b_ele_scaled, b_Iele_scaledE2, b_Lat_scaled, b_ppt_scaled, b_cancov_scaled,
               b_ch_res_scaled) |>
  # Clean up names for plotting
  mutate(.variable = str_remove(.variable, "b_FireP|b_")) |> 
  mutate(PredPretty = case_when(.variable == "LMFire1M5" ~ "LMSF: 1-5yr",
                                .variable == "LMFire6M10" ~ "LMSF: 6-10yr",
                                .variable == "LMFire11M35" ~ "LMSF: 11-35yr",
                                .variable == "HFire1M5" ~ "HSF: 1-5yr",
                                .variable == "HFire6M10" ~ "HSF: 6-10yr",
                                .variable == "HFire11M35" ~ "HSF: 11-35yr",
                                .variable == "ele_scaled" ~ "Elevation",
                                .variable == "Iele_scaledE2" ~ "Elevation^2",
                                .variable == "Lat_scaled" ~ "Latitude",
                                .variable == "ppt_scaled" ~ "Precipitation",
                                .variable == "cancov_scaled" ~ "Canopy Cover",
                                .variable == "ch_res_scaled" ~ "Canopy Height")) |> 
  mutate(PredCat = case_when(str_detect(.variable, "Fire") ~ "Fire",
                             str_detect(.variable, "ele|Lat") ~ "Topo/Geo",
                             str_detect(.variable, "ppt|can|ch_") ~ "Habitat"))

# Create halfeye plot
ggplot(fire_effects, aes(y = fct_reorder(PredPretty, .value), x = .value)) +
  stat_halfeye(
    .width = c(0.95),
    point_interval = "median_qi",
    fill = "lightblue",
    alpha = 0.7
  ) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw() +
  labs(x = "LCBD", 
       y = "") +
  theme(axis.text.y = element_text(hjust = 0),
        axis.text.x = element_text(angle = 90, vjust = 1)) +
  coord_flip() +
  facet_wrap(~PredCat, scales = "free_x")# Left-align category labels


## Patchwork
library(patchwork)

## Categorical effect plot
cat_plot <- ggplot(fire_effects |> filter(str_detect(.variable, "Fire")), 
                   aes(y = fct_reorder(PredPretty, .value), x = .value)) +
  stat_halfeye(
    .width = c(0.95),
    point_interval = "median_qi",
    fill = "lightblue",
    alpha = 0.7
  ) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw() +
  labs(x = "LCBD (x 1593)", 
       y = "",
       title = "Fire Effects") +
  theme(axis.text.y = element_text(hjust = 0)) +
  coord_flip()

# Get conditional effects data
ce_data <- conditional_effects(m1.me)

plot_theme_cat <- theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90),
    legend.position = "none")

plot_theme_lin <- theme_bw()


## Categorical effect plot
cat_plot <- 
  ce_data$FireP |> 
  mutate(PredPretty = case_when(FireP == "Unburned" ~ "Unburned",
    FireP == "LM Fire 1-5" ~ "LMSF: 1-5yr",
                   FireP == "LM Fire 6-10" ~ "LMSF: 6-10yr",
                   FireP == "LM Fire 11-35" ~ "LMSF: 11-35yr",
                   FireP == "H Fire 1-5" ~ "HSF: 1-5yr",
                   FireP == "H Fire 6-10" ~ "HSF: 6-10yr",
                   FireP == "H Fire 11-35" ~ "HSF: 11-35yr")) |>
  mutate(PredPretty = factor(PredPretty, 
                             levels = c("Unburned",
                                        "LMSF: 1-5yr",
                                        "LMSF: 6-10yr",
                                        "LMSF: 11-35yr",
                                        "HSF: 1-5yr",
                                        "HSF: 6-10yr",
                                        "HSF: 11-35yr"))) |> 
  ggplot(aes(y = PredPretty, x = estimate__, color = PredPretty)) +
  geom_linerange(aes(xmin = lower__, xmax = upper__), 
                 linewidth = 1) +
  geom_point(size = 2) + 
  scale_color_manual(values = c("#0571b0",
                                "#f4a582",
                                "#f4a582",
                                "#f4a582",
                                "#ca0020",
                                "#ca0020",
                                "#ca0020")) +
  theme_bw() +
  labs(x = "LCBD (x 1593)", 
       y = "") +
  coord_flip() + 
  plot_theme_cat

# Create list of linear effect plots with custom ggplot syntax
linear_plots <- list(
  # Elevation plot
  ggplot(ce_data$ele_scaled, aes(x = ele_scaled, y = estimate__)) +
    geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.2) +
    geom_line() +
    theme_bw() +
    labs(x = "Elevation", y = "LCBD (x 1593)") +
    plot_theme_lin,
    #theme(axis.text.x = element_text(angle = 0, hjust = 1)),
  
  # Latitude plot
  ggplot(ce_data$Lat_scaled, aes(x = Lat_scaled, y = estimate__)) +
    geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.2) +
    geom_line() +
    theme_bw() +
    labs(x = "Latitude", y = "LCBD (x 1593)") +
    plot_theme_lin,
  
  ## Precipitation 
  ggplot(ce_data$ppt_scaled, aes(x = ppt_scaled, y = estimate__)) +
    geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.2) +
    geom_line() +
    theme_bw() +
    labs(x = "Precipitation", y = "LCBD (x 1593)") +
    plot_theme_lin,
  
  ## Canopy Cover
  ggplot(ce_data$cancov_scaled, aes(x = cancov_scaled, y = estimate__)) +
    geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.2) +
    geom_line() +
    theme_bw() +
    labs(x = "Canopy Cover", y = "LCBD (x 1593)") +
    plot_theme_lin,
  
  ## Canopy Height
  ggplot(ce_data$ch_res_scaled, aes(x = ch_res_scaled, y = estimate__)) +
    geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.2) +
    geom_line() +
    theme_bw() +
    labs(x = "Canopy Height", y = "LCBD (x 1593)") +
    plot_theme_lin
)

# Then combine with your categorical plot as before using patchwork
lcbd_pred_plot <- cat_plot + linear_plots[[1]] + linear_plots[[2]] + linear_plots[[3]] + linear_plots[[4]] + linear_plots[[5]] + plot_layout(guides = "collect")

lcbd_pred_plot <- cat_plot + linear_plots[[1]] + linear_plots[[2]] + 
  linear_plots[[3]] + linear_plots[[4]] + linear_plots[[5]] + 
  plot_layout(guides = "collect") &
  theme(plot.margin = margin(b = 30, t = 5, l = 5, r = 5))

lcbd_pred_plot <- cowplot::plot_grid(cat_plot, linear_plots[[1]], linear_plots[[2]], 
                   linear_plots[[3]], linear_plots[[4]], linear_plots[[5]], align = "v", axis = "t")

ggsave(plot = lcbd_pred_plot, filename = here("./Figures/Preliminary/LCBD_CEplots_Patchwork.jpg"),
       height = 6, width = 8, dpi = 600)
## -------------------------------------------------------------
##
## End Section:
##
## -------------------------------------------------------------






m0 <- glm(LCBD_scale ~ 1, 
          family = Gamma(link = "inverse"),
          data = envTab)

## Fitting a full model for LCBD
m1 <- glmmTMB(log(LCBD) ~ ele + I(ele^2) + Lat + ppt + cancov + ch_res + FireP,
              family = ,
              data = envTab)

res = simulateResiduals(m1, plot = T)

summary(m1)

m2 <- glm(LCBD_scale ~ ele + I(ele^2) + Lat + ppt + cancov + ch_res + fire1_5yr_cbi_mn + fire6_10yr_cbi_mn + fire11_35yr_cbi_mn + FireP,
          family = Gamma(link = "inverse"),
          data = envTab)

1/coef(m1) 

AIC(m0, m1.geo, m1.hab, m1.habclim, m1.habfire, m1.fireonly)

library(ggeffects)
plot(ggpredict(m2, terms = "fire1_5yr_cbi_mn"), show_data = F)
plot(ggpredict(m2, terms = "fire6_10yr_cbi_mn"), show_data = F)
plot(ggpredict(m2, terms = "fire11_35yr_cbi_mn"), show_data = F)
plot(ggpredict(m1, terms = "FireP"), show_data = F)
plot(ggpredict(m1, terms = "cancov"), show_data = F)
plot(ggpredict(m1, terms = "ch_res"), show_data = F)
plot(ggpredict(m1, terms = "ppt"), show_data = F)
plot(ggpredict(m1, terms = "Lat"), show_data = F)
plot(ggpredict(m1, terms = "ele"), show_data = F)

modsel.m1 <- MuMIn::dredge(m1, rank = "BIC")
head(modsel.m1)
plot(modsel.m1)

modsel.m1

# Try with scaled response and explicit initial values
envTab$LCBD_scale <- envTab$LCBD * 1000

betafit <- brm(bf(LCBD_scale ~ ele + Lat + ppt + cancov + ch_res + FireP), 
               data = envTab,
               family = Gamma(link = "inverse"),
               chains = 3,
               control = list(adapt_delta = 0.99))

DT_fit1




## Checking model fit
plot(check_model(m1, panel = F))

# shape: 1 divided by dispersion parameter
breg_shape <- 1/0.1014104

# scale: mean/shape
breg_scale <- as.numeric(exp(coef(breg)))/breg_shape

hist(envTab$LCBD, breaks = 40, freq = FALSE)
curve(dgamma(x, shape = breg_shape, scale = breg_scale[1]), 
      from = 0, to = 3, col = "red", add = TRUE, n = 100)


# Get predicted mean at average predictor values
newdata <- data.frame(
  fire1_5yr_cbi_mn = mean(envTab$fire1_5yr_cbi_mn),
  fire6_10yr_cbi_mn = mean(envTab$fire6_10yr_cbi_mn),
  fire11_35yr_cbi_mn = mean(envTab$fire11_35yr_cbi_mn)
)

# Get predicted mean
pred_mean <- predict(breg, newdata = newdata, type = "response")

# Calculate shape and scale parameters for this prediction
breg_shape <- 1/summary(breg)$dispersion
breg_scale <- pred_mean/breg_shape

# Plot
hist(envTab$LCBD, breaks = 40, freq = FALSE)
curve(dgamma(x, shape = breg_shape, scale = breg_scale), 
      from = 0, to = 1.5, col = "red", add = TRUE, n = 1000)


# Assuming your model is called 'mod'
library(DHARMa)  # Great for GLM diagnostics
library(performance)  # Additional diagnostic tools

# 1. Basic diagnostic plots
par(mfrow = c(2,2))
plot(breg)  # Shows 4 standard diagnostic plots

# 2. DHARMa residuals (highly recommended for GLMs)
simulationOutput <- simulateResiduals(breg)
plot(simulationOutput)  # Comprehensive diagnostic plots

# 3. Specific tests
# Check residuals vs predicted
plot(fitted(breg), residuals(breg, type = "deviance"),
     xlab = "Fitted values", 
     ylab = "Deviance residuals")
abline(h = 0, lty = 2)

# 4. Check for influential observations
plot(cooks.distance(breg))
abline(h = 4/length(fitted(breg)), col = "red")  # Rule of thumb threshold

# 5. Additional specific tests
# Check model assumptions
check_model(mod)  # from performance package

# Check overdispersion
check_overdispersion(mod)

# Summary of various model checks
model_diagnostics <- check_distribution(mod)
print(model_diagnostics)
