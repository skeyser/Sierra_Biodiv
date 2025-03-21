## -------------------------------------------------------------
##
## Script name: Multivariate Biodiversity Analyses
##
## Script purpose:
##
## Author: Spencer R Keyser
##
## Date Created: 2025-02-20
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
library(purrr)
library(ggplot2)
library(ggrepel)
library(here)
library(vegan)
library(adespatial)
library(ade4)
library(adegraphics)
library(sf)
library(CAbioacoustics)
## -------------------------------------------------------------

## Read in corrected species composition
load(here("./Data/JAGS_Data/Occ2GDM_Data_95thresh.Rdata"))

dat <- OccData
rm(OccData)

## Take a look
str(dat)

## Grab the posterior
zmat <- dat$Z.posterior
zmat <- zmat[,,1]
nrow(zmat)

## Test influence of preds using RDA
ind.miss <- which(rowSums(zmat) == 0)
zmat <- zmat[-ind.miss,]
nrow(zmat)

## Add environmental data
env <- dat$SiteMeta

## Env test
env.test <- env |> 
  slice(-ind.miss) |> 
  select(Cell_Unit, utmn, utme,
         ele:fire11_35yr_high_prop) |>
  mutate(across(where(is.numeric), 
                list(scaled = ~scale(., center = TRUE, scale = TRUE)))) |> 
  as.data.frame()

## Bring in the pyrodiversity metric
pyro <- read.csv(here("./Data/FireSeverity2021_MeanStDev_AllARUs.csv"))
interval_naming <- c(
  "2016_2020" = "1_5",
  "2011_2015" = "6_10",
  "1986_2010" = "11_35"
)
pyro <- pyro |> 
  select(deployment_name, contains("stdev")) |> 
  rename_with(
    ~str_replace_all(., interval_naming),
    contains("Fire_Sev")
  ) |> 
  mutate(Cell_Unit = gsub("G\\d{3}_V\\d{1}_", "", deployment_name)) |> 
  select(Cell_Unit, contains("Fire_Sev")) #|> 
  filter(Cell_Unit %in% env.test$Cell_Unit) |> 
  distinct(Cell_Unit, .keep_all = T)

## Add the pyrodiversity metric to the envTest
env.test <- left_join(env.test, pyro)


## Subset the data down to the points that have
## fire at all
env.fire <- env |>
  slice(-ind.miss) |> 
  filter(if_any(contains("fire"), ~. > 0)) |> 
  select(Cell_Unit, utmn, utme,
         ele:fire11_35yr_high_prop) |>
  mutate(across(where(is.numeric), 
                list(scaled = ~scale(., center = TRUE, scale = TRUE)))) |> 
  as.data.frame() 


fire.rows <- which(env.test$Cell_Unit %in% env.fire$Cell_Unit)

zmat.fire <- zmat[fire.rows,]
zmat <- zmat.fire

zmat_com <- adespatial::beta.div.comp(zmat, "J")
zmat_com$part
zmat_nest <- zmat_com$rich
zmat_turn <- zmat_com$repl
zmat_tot <- zmat_com$D

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Create dbMEMs for Coordinates
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Load the data with the coords and UTM zones
aru_meta <- read.csv(here("./Data/ARU_Meta_120_FilteredOcc.csv"))
aru_meta <- aru_meta |> 
  select(Cell_Unit, Year = survey_yea, utmn, utme, utm_zone) |> 
  filter(Cell_Unit %in% env.test$Cell_Unit) |> 
  group_split(utm_zone) |>
  map_dfr(cb_make_aru_sf) |> 
  mutate(x = st_coordinates(geometry)[,1],
         y = st_coordinates(geometry)[,2]) |> 
  arrange(match(Cell_Unit, env.test$Cell_Unit))

coords <- cbind(aru_meta$x, aru_meta$y)
hist(coords[,1])

mapview::mapview(aru_meta)

## Make the MEMs
MEMs <- dbmem(coords, MEM.autocor = "non-null", silent = F)
summary(MEMs)

s.label(coords, neig = attr(MEMs, "listw"))

s.value(coords, MEMs[,1:3])

## Association with Moran's I
test <- moran.randtest(MEMs, nrepet = 999)
plot(test$obs, attr(MEMs, "values"), xlab = "Moran's I", ylab = "Eigenvalues")

sign_idx <- which(test$pvalue < 0.05)
sign_MEMs <- MEMs[, sign_idx]
sign_MEMs <- as.data.frame(sign_MEMs)

## Visualize the MEMs in space
coord.sf <- aru_meta |> 
  bind_cols(sign_MEMs) |> 
  tidyr::pivot_longer(cols = starts_with("MEM"),
                      names_to = "MEM",
                      values_to = "MEMvals")
  

coord.sf |> 
  filter(MEM %in% ) |> 
  ggplot() + 
  geom_sf(aes(color = MEMvals)) + 
  facet_wrap(~MEM)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Fit dbRDAs
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#env.test <- env.test[-ind.miss,]
mout1 <- dbrda(zmat_nest ~ utmn_scaled + ele_scaled + stage_scaled + cancov_scaled + fire1_5yr_high_prop_scaled + fire6_10yr_high_prop_scaled 
               + fire11_35yr_high_prop_scaled + fire1_5yr_lowmod_prop_scaled + fire6_10yr_lowmod_prop_scaled 
               + fire11_35yr_lowmod_prop_scaled, data = env.test)
summary(mout1)
screeplot(mout1)

mout2 <- dbrda(zmat_turn ~ utmn_scaled + ele_scaled + stage_scaled + cancov_scaled + fire1_5yr_high_prop_scaled + fire6_10yr_high_prop_scaled 
               + fire11_35yr_high_prop_scaled + fire1_5yr_lowmod_prop_scaled + fire6_10yr_lowmod_prop_scaled 
               + fire11_35yr_lowmod_prop_scaled, data = env.test)
summary(mout2)
screeplot(mout2)

mout3 <- dbrda(zmat_tot ~ utmn_scaled + ele_scaled + stage_scaled + cancov_scaled + fire1_5yr_high_prop_scaled + fire6_10yr_high_prop_scaled 
               + fire11_35yr_high_prop_scaled + fire1_5yr_lowmod_prop_scaled + fire6_10yr_lowmod_prop_scaled 
               + fire11_35yr_lowmod_prop_scaled, data = env.test)
summary(mout3)
screeplot(mout3)

R2adjNest <- RsquareAdj(mout1)
R2adjRepl <- RsquareAdj(mout2)
R2adjBtot <- RsquareAdj(mout3)
vif.cca(mout1)

## Fire points only
env.fire.test <- env.fire |> 
  select(contains("scaled"))

mout1 <- dbrda(zmat_nest ~ . + Condition(MEM1, MEM2), data = env.fire.test)
summary(mout1)
mout2 <- dbrda(zmat_turn ~ ele_scaled + stage_scaled + cancov_scaled + fire1_5yr_high_prop_scaled + fire6_10yr_high_prop_scaled 
               + fire11_35yr_high_prop_scaled + fire1_5yr_lowmod_prop_scaled + fire6_10yr_lowmod_prop_scaled 
               + fire11_35yr_lowmod_prop_scaled + Condition(MEM1, MEM2), data = env.fire.test)

vp <- vegan::varpart(zmat_turn, env.test[, c(2,4)], env.test[,5:6], env.test[, 7:12])
plot(vp)
summary(mout2)

vp <- vegan::varpart(zmat_tot, env.fire.test[, 1:3], env.fire.test[,4:9], env.fire.test[, 10:11])
plot(vp)
mout3 <- dbrda(zmat_tot ~ utmn_scaled + utme_scaled + ele_scaled + stage_scaled + cancov_scaled + fire1_5yr_high_prop_scaled + fire6_10yr_high_prop_scaled 
               + fire11_35yr_high_prop_scaled + fire1_5yr_lowmod_prop_scaled + fire6_10yr_lowmod_prop_scaled 
               + fire11_35yr_lowmod_prop_scaled + Condition(sign_MEMs$MEM1, sign_MEMs$MEM2), data = env.fire.test)

summary(mout3)



mout.aov1 <- anova(mout1, by = "terms")
mout.aov1
mout.aov2 <- anova(mout2, by = "terms")

## Total beta diversity
## By terms
mout.aov3 <- anova(mout3, by = "terms")
## Overall
mout.aov3global <- anova(mout3)

mout <- dbrda(zmat ~ ., data = env.test, dist = "raup")
summary(mout)

mout.aov <- anova(mout, by = "terms")

R2_dbrda <- RsquareAdj(mout3)$r.squared

# Forward selection of variables:
back.sel <- ordiR2step(rda(zmat ~ 1, data = env.test), # lower model limit (simple!)
                      scope = formula(mout3), # upper model limit (the "full" model)
                      direction = "backward",
                      R2scope = TRUE, # can't surpass the "full" model's R2
                      pstep = 1000,
                      trace = FALSE) # change to TRUE to see the selection process!


plot_dbrda_density <- function(dbrda_result) {
  # Extract site scores
  site_scores <- as.data.frame(scores(dbrda_result, display = "sites"))
  env_scores <- as.data.frame(scores(dbrda_result, display = "bp"))
  
  rownames(env_scores) <- c("Latitude", 
                            "Elevation", 
                            "Stand Age",
                            "Canopy Cover", 
                            "Fire: Low/Mod 1-5yr",
                            "Fire: Low/Mod 6-10yr",
                            "Fire: Low/Mod 11-35yr",
                            "Fire: High 1-5yr",
                            "Fire: High 6-10yr",
                            "Fire: High 11-35yr ")
  
  # Create density plot
  p <- ggplot() +
    geom_hex(data = site_scores,
             aes(x = dbRDA1, y = dbRDA2),
             bins = 20,
             color = "white",
             alpha = 0.8) +
    # Add density contours
    # stat_density_2d(data = site_scores,
    #                 aes(x = dbRDA1, y = dbRDA2, fill = ..level..),
    #                 geom = "polygon") +
    # Add environmental vectors
    geom_segment(data = env_scores,
                 aes(x = 0, y = 0, xend = dbRDA1, yend = dbRDA2),
                 arrow = arrow(length = unit(0.1, "cm"),
                               type = "closed"),
                 linewidth = 1,
                 color = "#8856a7") +
    # # Add variable labels
    # geom_text(data = env_scores,
    #           aes(x = dbRDA1, y = dbRDA2, label = rownames(env_scores)),
    #           hjust = -0.2) +
    geom_text_repel(data = env_scores,
                    aes(x = dbRDA1, y = dbRDA2, label = rownames(env_scores)),
                    box.padding = 0.8,
                    point.padding = 0.8,
                    force = 1,
                    force_pull = 0,
                    fontface = "bold",
                    segment.color = "black",
                    max.overlaps = Inf) + 
    #scale_fill_viridis_c(name = "Count", option = "C") +
    scale_fill_gradientn(name = "Count",
                         colors = RColorBrewer::brewer.pal(9, "PuBuGn")[2:9]) +
    theme_bw() +
    coord_equal()
  
  return(p)
}

dbrdaHex <- plot_dbrda_density(dbrda_result = mout1)
dbrdaHex <- plot_dbrda_density(dbrda_result = mout2)
dbrdaHex <- plot_dbrda_density(dbrda_result = mout3)

ggsave(plot = dbrdaHex, filename = here("./Figures/Preliminary/dbRDA_Hex_Biplot_RC.jpg"),
       height = 6, width = 6, dpi = 600)

# Check the new model with forward-selected variables
fwd.sel$call

## Test for significant drivers
mout.sig <- anova.cca(mout2, step = 1000, by = "term")
