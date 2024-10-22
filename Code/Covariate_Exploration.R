## -------------------------------------------------------------
##
## Script name: Covariate Exploration
##
## Script purpose: Explore the availability of covariate information
## for the Sierra ARU sites
##
## Author: Spencer R Keyser
##
## Date Created: 2024-10-15
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
library(sf)
library(terra)
library(patchwork)
## -------------------------------------------------------------

## Load in the ARU metadata
aru.meta <- readr::read_csv(here("Data/ARU_120m.csv"))
glimpse(aru.meta)
names(aru.meta)

## Make the dataframe an SF object
aru.meta <- aru.meta |> 
  st_as_sf(coords = c("X", "Y"), crs = 4326)


## -------------------------------------------------------------
##
## Begin Section: Histograms and Maps
##
## -------------------------------------------------------------

## Land Cover Variables

## Relative Height GEDI
hist(aru.meta$rh100mean_gedi_mn)
hist(aru.meta$rh100sd_gedi_mn)

## Climate Variables
hist(aru.meta$pet_bcm_mn)
hist(aru.meta$pet_bcm_sd)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Fire Variables
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Fire 1 yr
## CBI
h.f1_cbi <- ggplot(aru.meta) + 
  geom_histogram(aes(x = fire1yr_cbi_mn)) + 
  theme_bw() + 
  xlab("Mean CBI Fire 1yr Post")

m.f1_cbi <- ggplot(aru.meta) + 
  geom_sf(aes(color = fire1yr_cbi_mn)) + 
  scale_color_viridis_c() + 
  labs(color = "Mean CBI \nFire 1yr Post") + 
  theme_bw()

m.f1_cbi + h.f1_cbi

## Proportion Any Burn
h.f1_anyprop <- ggplot(aru.meta) + 
  geom_histogram(aes(x = fire1yr_all_prop)) + 
  theme_bw() + 
  xlab("Proportion Burned Any Severity")

m.f1_anyprop <- ggplot(aru.meta) + 
  geom_sf(aes(color = fire1yr_all_prop)) + 
  scale_color_viridis_c() + 
  labs(color = "Proportion Burned \nAny Severity") + 
  theme_bw()

## Figure for any proportion burned over severity
(m.f1_anyprop + h.f1_anyprop) / (m.f1_cbi + h.f1_cbi)

# ggsave(filename = here("Figures/Exploration/PropAny_MeanCBI.jpg"),
#        height = 8, width = 10, dpi = 600)


## Proportion Low-Moderate
h.f1_lowmod <- ggplot(aru.meta) + 
  geom_histogram(aes(x = fire1yr_lowmod_prop)) + 
  theme_bw() + 
  xlab("Proportion Burned Low-Moderate Severity")

m.f1_lowmod <- ggplot(aru.meta) + 
  geom_sf(aes(color = fire1yr_lowmod_prop)) + 
  scale_color_viridis_c() + 
  labs(color = "Proportion Burned \nLow-Moderate Severity") + 
  theme_bw()

## Figure for low-mod proportion burned over severity
(m.f1_lowmod + h.f1_lowmod) / (m.f1_cbi + h.f1_cbi)

ggsave(filename = here("Figures/Exploration/PropLowMod_MeanCBI.jpg"),
       height = 8, width = 10, dpi = 600)


## Proportion High Severity
h.f1_high <- ggplot(aru.meta) + 
  geom_histogram(aes(x = fire1yr_high_prop)) + 
  theme_bw() + 
  xlab("Proportion Burned High Severity")

m.f1_high <- ggplot(aru.meta) + 
  geom_sf(aes(color = fire1yr_high_prop)) + 
  scale_color_viridis_c() + 
  labs(color = "Proportion Burned \nHigh Severity") + 
  theme_bw()

## Figure for low-mod proportion burned over severity
(m.f1_high + h.f1_high) / (m.f1_cbi + h.f1_cbi)

ggsave(filename = here("Figures/Exploration/PropHi_MeanCBI.jpg"),
       height = 8, width = 10, dpi = 600)


## Fire 2-5 yr
## CBI
h.f2_5_cbi <- ggplot(aru.meta) + 
  geom_histogram(aes(x = fire2_5yr_cbi_mn)) + 
  theme_bw() + 
  xlab("Mean CBI Fire 2-5yr Post")

m.f2_5_cbi <- ggplot(aru.meta) + 
  geom_sf(aes(color = fire2_5yr_cbi_mn)) + 
  scale_color_viridis_c() + 
  labs(color = "Mean CBI \nFire 2-5yr Post") + 
  theme_bw()

m.f2_5_cbi + h.f2_5_cbi

## Proportion Any Burn
h.f2_5_anyprop <- ggplot(aru.meta) + 
  geom_histogram(aes(x = fire2_5yr_all_prop)) + 
  theme_bw() + 
  xlab("Proportion Burned Any Severity")

m.f2_5_anyprop <- ggplot(aru.meta) + 
  geom_sf(aes(color = fire2_5yr_all_prop)) + 
  scale_color_viridis_c() + 
  labs(color = "Proportion Burned \nAny Severity") + 
  theme_bw()

## Figure for any proportion burned over severity
(m.f2_5_anyprop + h.f2_5_anyprop) / (m.f2_5_cbi + h.f2_5_cbi)

# ggsave(filename = here("Figures/Exploration/PropAnyMeanCBI_2_5.jpg"),
#        height = 8, width = 10, dpi = 600)


## Proportion Low-Moderate
h.f2_5_lowmod <- ggplot(aru.meta) + 
  geom_histogram(aes(x = fire2_5yr_lowmod_prop)) + 
  theme_bw() + 
  xlab("Proportion Burned Low-Moderate Severity")

m.f2_5_lowmod <- ggplot(aru.meta) + 
  geom_sf(aes(color = fire2_5yr_lowmod_prop)) + 
  scale_color_viridis_c() + 
  labs(color = "Proportion Burned \nLow-Moderate Severity") + 
  theme_bw()

## Figure for low-mod proportion burned over severity
(m.f2_5_lowmod + h.f2_5_lowmod) / (m.f2_5_cbi + h.f2_5_cbi)

# ggsave(filename = here("Figures/Exploration/PropLowModMeanCBI_2_5.jpg"),
#        height = 8, width = 10, dpi = 600)


## Proportion High Severity
h.f2_5_high <- ggplot(aru.meta) + 
  geom_histogram(aes(x = fire2_5yr_high_prop)) + 
  theme_bw() + 
  xlab("Proportion Burned High Severity")

m.f2_5_high <- ggplot(aru.meta) + 
  geom_sf(aes(color = fire2_5yr_high_prop)) + 
  scale_color_viridis_c() + 
  labs(color = "Proportion Burned \nHigh Severity") + 
  theme_bw()

## Figure for low-mod proportion burned over severity
(m.f2_5_high + h.f2_5_high) / (m.f2_5_cbi + h.f2_5_cbi)

## Fire 6-10 yr
## CBI
h.f6_10_cbi <- ggplot(aru.meta) + 
  geom_histogram(aes(x = fire6_10yr_cbi_mn)) + 
  theme_bw() + 
  xlab("Mean CBI Fire 6-10yr Post")

m.f6_10_cbi <- ggplot(aru.meta) + 
  geom_sf(aes(color = fire6_10yr_cbi_mn)) + 
  scale_color_viridis_c() + 
  labs(color = "Mean CBI \nFire 6-10yr Post") + 
  theme_bw()

m.f6_10_cbi + h.f6_10_cbi

## Proportion Any Burn
h.f6_10_anyprop <- ggplot(aru.meta) + 
  geom_histogram(aes(x = fire6_10yr_all_prop)) + 
  theme_bw() + 
  xlab("Proportion Burned Any Severity")

m.f6_10_anyprop <- ggplot(aru.meta) + 
  geom_sf(aes(color = fire6_10yr_all_prop)) + 
  scale_color_viridis_c() + 
  labs(color = "Proportion Burned \nAny Severity") + 
  theme_bw()

## Figure for any proportion burned over severity
(m.f6_10_anyprop + h.f6_10_anyprop) / (m.f6_10_cbi + h.f6_10_cbi)

## Proportion Low-Moderate
h.f6_10_lowmod <- ggplot(aru.meta) + 
  geom_histogram(aes(x = fire6_10yr_lowmod_prop)) + 
  theme_bw() + 
  xlab("Proportion Burned Low-Moderate Severity")

m.f6_10_lowmod <- ggplot(aru.meta) + 
  geom_sf(aes(color = fire6_10yr_lowmod_prop)) + 
  scale_color_viridis_c() + 
  labs(color = "Proportion Burned \nLow-Moderate Severity") + 
  theme_bw()

## Figure for low-mod proportion burned over severity
(m.f6_10_lowmod + h.f6_10_lowmod) / (m.f6_10_cbi + h.f6_10_cbi)

## Proportion High Severity
h.f6_10_high <- ggplot(aru.meta) + 
  geom_histogram(aes(x = fire6_10yr_high_prop)) + 
  theme_bw() + 
  xlab("Proportion Burned High Severity")

m.f6_10_high <- ggplot(aru.meta) + 
  geom_sf(aes(color = fire6_10yr_high_prop)) + 
  scale_color_viridis_c() + 
  labs(color = "Proportion Burned \nHigh Severity") + 
  theme_bw()

## Figure for low-mod proportion burned over severity
(m.f6_10_high + h.f6_10_high) / (m.f6_10_cbi + h.f6_10_cbi)


## Fire 11-35 yr
## CBI
h.f11_35_cbi <- ggplot(aru.meta) + 
  geom_histogram(aes(x = fire11_35yr_cbi_mn)) + 
  theme_bw() + 
  xlab("Mean CBI Fire 11-35yr Post")

m.f11_35_cbi <- ggplot(aru.meta) + 
  geom_sf(aes(color = fire11_35yr_cbi_mn)) + 
  scale_color_viridis_c() + 
  labs(color = "Mean CBI \nFire 11-35yr Post") + 
  theme_bw()

m.f11_35_cbi + h.f11_35_cbi

## Proportion Any Burn
h.f11_35_anyprop <- ggplot(aru.meta) + 
  geom_histogram(aes(x = fire11_35yr_all_prop)) + 
  theme_bw() + 
  xlab("Proportion Burned Any Severity")

m.f11_35_anyprop <- ggplot(aru.meta) + 
  geom_sf(aes(color = fire11_35yr_all_prop)) + 
  scale_color_viridis_c() + 
  labs(color = "Proportion Burned \nAny Severity") + 
  theme_bw()

## Figure for any proportion burned over severity
(m.f11_35_anyprop + h.f11_35_anyprop) / (m.f11_35_cbi + h.f11_35_cbi)

## Proportion Low-Moderate
h.f11_35_lowmod <- ggplot(aru.meta) + 
  geom_histogram(aes(x = fire11_35yr_lowmod_prop)) + 
  theme_bw() + 
  xlab("Proportion Burned Low-Moderate Severity")

m.f11_35_lowmod <- ggplot(aru.meta) + 
  geom_sf(aes(color = fire11_35yr_lowmod_prop)) + 
  scale_color_viridis_c() + 
  labs(color = "Proportion Burned \nLow-Moderate Severity") + 
  theme_bw()

## Figure for low-mod proportion burned over severity
(m.f11_35_lowmod + h.f11_35_lowmod) / (m.f11_35_cbi + h.f11_35_cbi)

## Proportion High Severity
h.f11_35_high <- ggplot(aru.meta) + 
  geom_histogram(aes(x = fire11_35yr_high_prop)) + 
  theme_bw() + 
  xlab("Proportion Burned High Severity")

m.f11_35_high <- ggplot(aru.meta) + 
  geom_sf(aes(color = fire11_35yr_high_prop)) + 
  scale_color_viridis_c() + 
  labs(color = "Proportion Burned \nHigh Severity") + 
  theme_bw()

## Figure for low-mod proportion burned over severity
(m.f11_35_high + h.f11_35_high) / (m.f11_35_cbi + h.f11_35_cbi)


## Look at the proportion burned at high severity across
## time since fire

## High severity
m.f1_high + m.f2_5_high + m.f6_10_high + m.f11_35_high + plot_layout(guides = "collect")
h.f1_high + h.f2_5_high + h.f6_10_high + h.f11_35_high + plot_layout(guides = "collect")

## Low-moderate severity
m.f1_lowmod + m.f2_5_lowmod + m.f6_10_lowmod + m.f11_35_lowmod + plot_layout(guides = "collect")
h.f1_lowmod + h.f2_5_lowmod + h.f6_10_lowmod + h.f11_35_lowmod + plot_layout(guides = "collect")

## All fires
m.f1_anyprop + m.f2_5_anyprop + m.f6_10_anyprop + m.f11_35_anyprop + plot_layout(guides = "collect")
h.f1_anyprop + h.f2_5_anyprop + h.f6_10_anyprop + h.f11_35_anyprop + plot_layout(guides = "collect")
mapview::mapview(aru.meta,
                 zcol = "fire1yr_high_prop")

mapview::mapview(aru.meta,
                 zcol = "fire1yr_lowmod_prop")

mapview::mapview(aru.meta,
                 zcol = "fire1yr_all_prop")

mapview::mapview(aru.meta,
                 zcol = "fire2_5yr_all_prop")

mapview::mapview(aru.meta,
                 zcol = "fire6_10yr_all_prop")

mapview::mapview(aru.meta,
                 zcol = "fire11_35yr_all_prop")

## Mean CBI across the classes
yr_cbi <- m.f1_cbi + m.f2_5_cbi + m.f6_10_cbi + m.f11_35_cbi
yr_cbi_hist <- h.f1_cbi + h.f2_5_cbi + h.f6_10_cbi + h.f11_35_cbi
ggsave(here("Figures/Exploration/CBI_by_PostFireYr.jpg"),
       plot = yr_cbi, height = 8, width = 12, dpi = 600)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Treatments
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Change this for ANUs data

## Treatment 36 yr
h.t36_prop <- ggplot(aru.meta) + 
  geom_histogram(aes(x = all_treat_36yr_prop)) + 
  theme_bw() + 
  xlab("Treatment proportion: 36 year")

m.t36_prop <- ggplot(aru.meta) + 
  geom_sf(aes(color = all_treat_36yr_prop)) + 
  scale_color_viridis_c() + 
  labs(color = "Treatment proportion: 36 year") + 
  theme_bw()

## Treatment 5 yr
h.t5_prop <- ggplot(aru.meta) + 
  geom_histogram(aes(x = all_treat_5yr_prop)) + 
  theme_bw() + 
  xlab("Treatment proportion: 5 year")

m.t5_prop <- ggplot(aru.meta) + 
  geom_sf(aes(color = all_treat_5yr_prop)) + 
  scale_color_viridis_c() + 
  labs(color = "Treatment proportion: 5 year") + 
  theme_bw()

## -------------------------------------------------------------
##
## End Section: Histograms and Mapping
##
## -------------------------------------------------------------

## -------------------------------------------------------------
##
## Begin Section: Correlation Matrices
##
## -------------------------------------------------------------

## Correlation Matrices
aru.corr <- aru.meta |> 
  select(contains(c("fire", "aet", "ppt", "tmn_", "tmax_", "topo", "_treat_"))) |> 
  sf::st_drop_geometry()

cor.mat <- cor(aru.corr, use = "complete.obs")

corrplot(cor.mat, method = "circle")

plot(aru.corr$all_treat_36yr_prop, aru.corr$fire1yr_cbi_mn)
plot(aru.corr$all_treat_36yr_prop, aru.corr$fire2_5yr_cbi_mn)
plot(aru.corr$all_treat_36yr_prop, aru.corr$fire6_10yr_cbi_mn)
plot(aru.corr$all_treat_36yr_prop, aru.corr$fire11_35yr_cbi_mn)

## Number of ARUs 