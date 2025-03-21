## -------------------------------------------------------------
##
## Script name: ARU alpha diversity models
##
## Script purpose:
##
## Author: Spencer R Keyser
##
## Date Created: 2024-11-07
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
library(MASS)
## -------------------------------------------------------------

## Reading in data for GDM
sierra <- readr::read_csv(here::here("Data/Generated_DFs/2021_SeasonalSpeciesMat.csv"))

## Species Richness df
rich <- sierra |> 
  select(Cell_Unit, Long, Lat, 5:105) |>
  mutate(SR = rowSums(across(-c(Cell_Unit, Long, Lat)))) |> 
  select(Cell_Unit, Long, Lat, SR) |> 
  as.data.frame()

## SF for viz
rich.sf <- rich |> 
  st_as_sf(coords = c("Long", "Lat"), crs = 4326)

ggplot(data = rich.sf) +
  geom_sf(aes(color = SR)) + 
  scale_color_viridis_c()

## Merge with the ARU covariate data
aru_meta <- readr::read_csv(here("Data/ARU_120m.csv"))
names(aru_meta)
glimpse(aru_meta)
aru_meta$Cell_Unit <- paste0(aru_meta$cell_id, "_", aru_meta$unit_numbe)

## Merge
rich.meta <- rich |> 
  left_join(aru_meta, by = "Cell_Unit")

## Preliminary species richness models
## Reduce the dataset to some covariate of interest
rich.mod <- rich.meta |> 
  dplyr::select(Cell_Unit, Long, Lat, 
         SR, #species richness
         rh100sd_gedi_mn, #structural het
         cc_cfo_mn, #habitat
         tmx_bcm_mn, #thermal
         ppt_bcm_mn, #aridity
         topo_tpi, #topo het
         topo_elev #elev
         )


## Simple species richness model
## Standardize the covariates
rich.mod <- rich.mod |> 
  mutate(rh_sd_scale = scale(rh100sd_gedi_mn),
         cc_cfo_mn_scale = scale(cc_cfo_mn),
         tmx_bcm_mn_scale = scale(tmx_bcm_mn),
         ppt_bcm_mn_scale = scale(ppt_bcm_mn),
         topo_tpi_scale = scale(topo_tpi),
         topo_elev_scale = scale(topo_elev),
         Long_scale = scale(Long),
         Lat_scale = scale(Lat)) |> 
  tidyr::drop_na()

## Histograms
hist(rich.mod$SR)
hist(rich.mod$rh100sd_gedi_mn)
hist(rich.mod$cc_cfo_mn_scale)
hist(rich.mod$tmx_bcm_mn_scale)
hist(rich.mod$ppt_bcm_mn_scale)
hist(rich.mod$topo_tpi_scale)
hist(rich.mod$topo_elev_scale)

## Early predictor exploration
hab_het <- ggplot(data = rich.mod, mapping = aes(x = rh100sd_gedi_mn, y = SR)) + 
  geom_point() + 
  geom_smooth(method="loess") + 
  xlab('SD Canopy Height') + 
  ylab('Species richness')

can_cov <- ggplot(data = rich.mod, mapping = aes(x = cc_cfo_mn, y = SR)) + 
  geom_point() + 
  geom_smooth(method="loess") + 
  xlab('% Canopy Cover') + 
  ylab('Species richness')

therm_stress <- ggplot(data = rich.mod, mapping = aes(x = tmx_bcm_mn, y = SR)) + 
  geom_point() + 
  geom_smooth(method="loess") + 
  xlab('Max. Temp') + 
  ylab('Species richness')

precip <- ggplot(data = rich.mod, mapping = aes(x = ppt_bcm_mn, y = SR)) + 
  geom_point() + 
  geom_smooth(method="loess") + 
  xlab('Precipitation') + 
  ylab('Species richness')

elev_complex <- ggplot(data = rich.mod, mapping = aes(x = topo_tpi, y = SR)) + 
  geom_point() + 
  geom_smooth(method="loess") + 
  xlab('Topo. Pos. Index') + 
  ylab('Species richness')

elev_mean <- ggplot(data = rich.mod, mapping = aes(x = topo_elev, y = SR)) + 
  geom_point() + 
  geom_smooth(method="loess") + 
  xlab('Mean Elevation') + 
  ylab('Species richness')


## Predictor collinearity
# whether all variables are normally distributed.
cor_mat <- cor(rich.mod |> dplyr::select(contains("_scale")), method='spearman')

# We can visualise this correlation matrix. For better visibility, 
# we plot the correlation coefficients as percentages.
corrplot(cor_mat, method = "number")

## GLM poisson link for species richness
sr.mod1 <- glm(SR ~ Lat_scale*topo_tpi_scale + rh_sd_scale + cc_cfo_mn_scale + tmx_bcm_mn_scale + ppt_bcm_mn_scale,
               family = "poisson", data = rich.mod)


sr.mod1 <- glm.nb(SR ~ Lat_scale*topo_tpi_scale + rh_sd_scale + cc_cfo_mn_scale + tmx_bcm_mn_scale + ppt_bcm_mn_scale,
               data = rich.mod)


summary(sr.mod1)

## Deviance explained
## Pretty bad model here only 4% of the deviance explained
dev = 1- deviance(sr.mod1) / sr.mod1$null.deviance

## Fitt Mod1 
plot(sr.mod1$fitted.values, rich.mod$SR, xlab = "Fitted values", ylab = "Observed claims")
abline(lm(sr.mod1$fitted ~ rich.mod$SR), col="light blue", lwd=2)
abline(0, 1, col = "dark blue", lwd=2)

## Deviance
anova(sr.mod1, test="F")
1 - pchisq(deviance(sr.mod1), df = sr.mod1$df.residual)

## Pearson gof
# Pearson's goodness-of-fit
Pearson <- sum((rich.mod$SR - sr.mod1$fitted.values)^2 
               / sr.mod1$fitted.values)
1 - pchisq(Pearson, df = sr.mod1$df.residual)

## Overdispersion test
phi <- sum(residuals(sr.mod1, "pearson")^2)/sr.mod1$df.residual
library(AER)
dispersiontest(sr.mod1, trafo = 1)

lambdahat <-fitted(sr.mod1)
par(mfrow=c(1,2), pty="s")
plot(lambdahat,(rich.mod$SR-lambdahat)^2,
     xlab=expression(hat(lambda)), ylab=expression((y-hat(lambda))^2 ))
plot(lambdahat, resid(sr.mod1,type="pearson"), 
     xlab=expression(hat(lambda)), ylab="Pearson Residuals") 
