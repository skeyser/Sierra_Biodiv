## -------------------------------------------------------------
##
## Script name: GDM Sandbox
##
## Script purpose: Explore the use of GDM for Sierra Bioacoustics
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

## GDM
library(gdm)
library(betapart)

## -------------------------------------------------------------

## Prepackaged GDM data
str(southwest)

sppTab <- southwest[, c("species", "site", "Long", "Lat")]

envTab <- southwest[, c(2:ncol(southwest))]

gdmTab <- formatsitepair(bioData = sppTab, 
                         bioFormat = 2,
                         XColumn = "Long",
                         YColumn = "Lat",
                         sppColumn = "species",
                         siteColumn = "site",
                         predData = envTab)

## Reading in data for GDM
sierra <- readr::read_csv(here::here("Data/Generated_DFs/2021_SeasonalSpeciesMat.csv"))
sierra <- readr::read_csv(here::here("Data/Generated_DFs/Seasonal_Summaries/2021_SeasonalSpeciesMat.csv"))
sierra <- sierra |> rename(Long = lon, Lat = lat)

## Aru meta
aru_meta <- readr::read_csv(here("Data/ARU_120m.csv"))
names(aru_meta)
glimpse(aru_meta)
aru_meta$Cell_Unit <- paste0(aru_meta$cell_id, "_", aru_meta$unit_numbe)

## Take a couple of things for testing
aru_meta <- aru_meta |>
  dplyr::select(Cell_Unit, survey_yea, X, Y, topo_elev, standage_f3_mn, contains("fire"))

## Join wit the sierra
sierra <- left_join(sierra, aru_meta)


sppTab <- sierra |> 
  select(Cell_Unit, Long, Lat, 6:106) |>
  #tidyr::separate(Cell_Unit, into = c("Cell", "Unit"), sep = "_", remove = T) |>
  #group_by(Cell) |> 
  #summarise(across(matches("Lat|Long"), mean, na.rm = T),
  #          across(!matches("Cell|Lat|Long"), max, na.rm = T)) |>
  #select(-Unit) |> 
  mutate(SR = rowSums(across(-c(Cell_Unit, Long, Lat)))) |> 
  filter(SR > 0) |> 
  select(-SR) |> 
  distinct(Cell_Unit, .keep_all = T) |> 
  as.data.frame()

## Calculate pairwise beta diversity across the study sites
sierra.beta <- sierra |> 
  select(Cell_Unit, Long, Lat, 5:105) |>
  mutate(SR = rowSums(across(-c(Cell_Unit, Long, Lat)))) |> 
  filter(SR > 0) |> 
  select(-c(SR, Long, Lat)) |>
  distinct(Cell_Unit, .keep_all = T) |> 
  #slice_sample(prop = 0.2) |> 
  as.data.frame()

rownames(sierra.beta) <- sierra.beta$Cell_Unit 
sierra.beta <- sierra.beta[, -1]
sierra.dist <- vegan::vegdist(sierra.beta, method = "bray")


sierra.pwbeta <- beta.pair(sierra.beta)
names(sierra.pwbeta)

nmds <- vegan::metaMDS(sierra.dist,
                       distance = "bray",
                       k = 3, 
                       maxit = 999,
                       trymax = 100,
                       wascores = TRUE)
vegan::stressplot(nmds)
plot(nmds)

par(mfrow = c(1,3))
hist(sierra.pwbeta[[1]], main = "Beta: Replacement", xlab = "Replacement")
hist(sierra.pwbeta[[2]], main = "Beta: Nestedness", xlab = "Nestedness")
hist(sierra.pwbeta[[3]], main = "Beta: Total", xlab = "Beta Diversity")

repl <- apply(as.matrix(sierra.pwbeta[[1]]), 1, mean)
nest <- apply(as.matrix(sierra.pwbeta[[2]]), 1, mean)
total <- apply(as.matrix(sierra.pwbeta[[3]]), 1, mean)

beta.df <- data.frame(repl = repl, nest = nest, total = total)
beta.df$Cell_Unit <- rownames(beta.df)
rownames(beta.df) <- NULL

sierra.geo <- sierra |> select(Cell_Unit, Long, Lat)

beta.df <- beta.df |> tidyr::pivot_longer(cols = c(repl, nest, total),
                                          names_to = "Component",
                                          values_to = "Value") |> 
  left_join(sierra.geo, by = "Cell_Unit")

ca <- read_sf(here("./Data/US_Map/states.shp"))
ca <- ca |> 
  filter(STATE_NAME == "California") |> 
  st_transform(crs = 4326)

beta.sf.repl <- beta.df |>
  filter(Component == "repl") |> 
  rename(Lat = Long, Long = Lat) |> 
  st_as_sf(coords = c("Long", "Lat"), crs = 4326) |> 
  ggplot() + 
  geom_sf(data = ca, color = "black", fill = "lightgrey") + 
  geom_sf(aes(color = Value)) + 
  scale_color_gradient(low = "#fee8c8",
                        high = "#e34a33") + 
  theme_void() +
  labs(color = "Replacement") +
  theme(
    legend.position = "bottom",        # Places the legend at the bottom
    legend.title = element_text(hjust = 0.5, size = 16),  # Centers the legend title
    legend.key.width = unit(1, "cm"),  # Increases the width of legend keys (optional)
    legend.key.height = unit(0.5, "cm"),  # Increases the height of legend keys (optional)
    plot.title = element_text(hjust = 0.5, size = 16)  # Centers the title
  ) + 
  guides(
    color = guide_colorbar(title.position = "top", title.hjust = 0.5)  # Title above the color bar
  )

beta.sf.nest <- beta.df |> 
  filter(Component == "nest") |> 
  rename(Lat = Long, Long = Lat) |> 
  st_as_sf(coords = c("Long", "Lat"), crs = 4326) |> 
  ggplot() + 
  geom_sf(data = ca, color = "black", fill = "lightgrey") + 
  geom_sf(aes(color = Value)) + 
  scale_color_gradient(low = "#e5f5f9",
                       high = "#2ca25f") +
  theme_void() + 
  labs(color = "Nestedness") +
  theme(
    legend.position = "bottom",        # Places the legend at the bottom
    legend.title = element_text(hjust = 0.5, size = 16),  # Centers the legend title
    legend.key.width = unit(1, "cm"),  # Increases the width of legend keys (optional)
    legend.key.height = unit(0.5, "cm"),  # Increases the height of legend keys (optional)
    plot.title = element_text(hjust = 0.5, size = 16)  # Centers the title
  ) + 
  guides(
    color = guide_colorbar(title.position = "top", title.hjust = 0.5)  # Title above the color bar
  )

beta.sf.total <- beta.df |>
  filter(Component == "total") |> 
  rename(Lat = Long, Long = Lat) |> 
  st_as_sf(coords = c("Long", "Lat"), crs = 4326) |> 
  ggplot() + 
  geom_sf(data = ca, color = "black", fill = "lightgrey") + 
  geom_sf(aes(color = Value)) + 
  scale_color_gradient(low = "#e0ecf4",
                       high = "#8856a7") +
  theme_void() +
  labs(color = "Total Beta") +
  theme(
    legend.position = "bottom",        # Places the legend at the bottom
    legend.title = element_text(hjust = 0.5, size = 16),  # Centers the legend title
    legend.key.width = unit(1, "cm"),  # Increases the width of legend keys (optional)
    legend.key.height = unit(0.5, "cm"),  # Increases the height of legend keys (optional)
    plot.title = element_text(hjust = 0.5, size = 16)  # Centers the title
  ) + 
  guides(
    color = guide_colorbar(title.position = "top", title.hjust = 0.5)  # Title above the color bar
  )

map.grid <- beta.sf.total / (beta.sf.nest | beta.sf.repl)
ggsave(plot = map.grid, filename = here("./Figures/Exploration/BetaMapComponents.jpg"),
       width = 10, height = 12, dpi = 600)
ggsave(plot = beta.sf.total, filename = here("./Figures/Exploration/BetaMap.jpg"),
       width = 10, height = 12, dpi = 600)
ggsave(plot = beta.sf.nest, filename = here("./Figures/Exploration/BetaMapNest.jpg"),
       width = 8, height = 12, dpi = 600)
ggsave(plot = beta.sf.repl, filename = here("./Figures/Exploration/BetaMapRepl.jpg"),
       width = 8, height = 12, dpi = 600)


beta.hist <- ggplot(data = beta.df) +
  geom_histogram(aes(x = Value, fill = Component), color = "black", alpha = 0.6) +
  scale_fill_manual(values = c("#1b9e77", "#d95f02", "#7570b3"),
                    labels = c("Nestedness", "Replacement", "Beta")) + 
  theme_bw() + 
  theme(axis.text = element_text(family = "sans", size = 12),
        axis.title = element_text(family = "sans", size = 16),
        legend.text = element_text(family = "sans", size = 16),
        legend.title = element_text(family = "sans", size = 16)) +
  xlab("Mean Dissimilarity") + 
  ylab("Count")

ggsave(plot = beta.hist,
       filename = here("./Figures/Exploration/BetaComponentHisto.jpg"),
       dpi = 600,
       height = 8,
       width = 10)
  

vegan::meandist(sierra.pwbeta[[3]])
simper.sierra <- vegan::simper(sppTab[,4:ncol(sppTab)],
                               permutations = 999)

which.max(simper.sierra$total$average)
max(simper.sierra$total$average)

## Environmental data columns
envTab <- sierra |> 
  select(Cell_Unit, Long, Lat, 
         topo_elev,
         standage_f3_mn,
         fire1yr_high_prop,
         fire2_5yr_high_prop,
         fire6_10yr_high_prop,
         fire11_35yr_high_prop
         ) |> 
  filter(Cell_Unit %in% sppTab$Cell_Unit) |>
  distinct(Cell_Unit, .keep_all = T) |> 
  as.data.frame()

## Format the data for GDM
gdmTab <- formatsitepair(bioData = sppTab,
                         bioFormat = 1,
                         XColumn = "Long",
                         YColumn = "Lat",
                         siteColumn = "Cell_Unit",
                         predData = envTab)

## Preformatted distance matrix
beta.rep <- as.array(as.matrix(sierra.pwbeta[[2]]))
str(beta.rep)
beta.rep <- cbind(rownames(beta.rep), beta.rep)
colnames(beta.rep)[colnames(beta.rep) == ""] <- "Cell_Unit"
beta.rep.order <- beta.rep[,1]
beta.rep[,1] <- seq(from = 1, to = nrow(beta.rep), by = 1)
beta.rep <- apply(beta.rep, 2, as.numeric)
class(beta.rep[,1])

envTab <- envTab[match(beta.rep.order, envTab$Cell_Unit),]
envTab$Cell_Unit <- seq(from = 1, to = nrow(envTab), by = 1)
envTab$Cell_Unit <- as.numeric(envTab$Cell_Unit)
class(envTab$Cell_Unit)

gdmTab <- formatsitepair(bioData = beta.rep,
                         bioFormat = 3,
                         XColumn = "Long",
                         YColumn = "Lat",
                         siteColumn = "Cell_Unit",
                         predData = envTab)

class(gdmTab)
hist(gdmTab$distance)
hist(gdmTab$s1.topo_elev)
hist(gdmTab$s1.fire1yr_high_prop)
hist(gdmTab$s1.fire2_5yr_high_prop)
hist(gdmTab$s1.fire6_10yr_high_prop)
hist(gdmTab$s1.fire11_35yr_high_prop)

## Fit GDM
gdm.fit1 <- gdm(gdmTab, geo = TRUE)
summary(gdm.fit1)

length(gdm.fit1$predictors)

plot(gdm.fit1, plot.layout = c(3,3))

## Significance testing via varImp
gdm::gdm.varImp(gdmTab, geo = T)

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

## Plot topo
gdm.spline <- isplineExtract(gdm.fit1)

plot(gdm.spline$x[,"topo_elev"], gdm.spline$y[,"topo_elev"], lwd=3, type="l", xlab="Elevation", ylab="Partial ecological distance")

topo.plot <- gdm.spline |> 
  map(~ as.data.frame(.x)) |>
  map(~ select(.x, topo_elev)) |> 
  bind_cols() |> 
  rename(x = `topo_elev...1`, y = `topo_elev...2`) |> 
  ggplot() + 
  geom_line(aes(x=x, y=y), color = "black", size = 2) +
  theme_bw() + 
  theme(axis.title = element_text(family = "sans", size = 16),
        axis.text = element_text(family = "sans", size = 14)) + 
  xlab("Elevation (meters)") + 
  ylab("Partial Ecological Distance")

ggsave(plot = topo.plot,
       filename = here("Figures/Exploration/GDM_Topo_pareff.jpg"),
       height = 8, width = 12, dpi = 600)

sage.plot <- gdm.spline |> 
  map(~ as.data.frame(.x)) |>
  map(~ select(.x, standage_f3_mn)) |> 
  bind_cols() |> 
  rename(x = `standage_f3_mn...1`, y = `standage_f3_mn...2`) |> 
  ggplot() + 
  geom_line(aes(x=x, y=y), color = "black", size = 2) +
  theme_bw() + 
  theme(axis.title = element_text(family = "sans", size = 16),
        axis.text = element_text(family = "sans", size = 14)) +
  xlab("Stand Age") +   
  ylab("Partial Ecological Distance")

ggsave(plot = sage.plot,
       filename = here("Figures/Exploration/GDM_StandAge_pareff.jpg"),
       height = 8, width = 12, dpi = 600)


fire.plot <- gdm.spline |> 
  map(~ as.data.frame(.x)) |>
  map(~ select(.x, fire6_10yr_high_prop)) |> 
  bind_cols() |> 
  rename(x = `fire6_10yr_high_prop...1`, y = `fire6_10yr_high_prop...2`) |> 
  ggplot() + 
  geom_line(aes(x=x, y=y), color = "black", size = 2) +
  theme_bw() + 
  theme(axis.title = element_text(family = "sans", size = 16),
        axis.text = element_text(family = "sans", size = 14)) +
  xlab("Proportion High Severity Fire: 6-10 years") +  
  ylab("Partial Ecological Distance")

ggsave(plot = fire.plot,
       filename = here("Figures/Exploration/GDM_Fire6-10_pareff.jpg"),
       height = 8, width = 12, dpi = 600)



library(patchwork)

pgrid <- topo.plot + sage.plot + fire.plot + plot_layout(guides = "collect")

ggsave(plot = pgrid,
       filename = here("Figures/Exploration/GDM_Mod_3_Resp.jpg"),
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
