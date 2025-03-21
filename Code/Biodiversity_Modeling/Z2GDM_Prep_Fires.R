## -------------------------------------------------------------
##
## Script name: Data for GDM - 2020 Fire Anaylsis
##
## Script purpose: Prepping output from the MSOM for Bayesian GDM
## for the smaller fire-specific analysis.
##
## Author: Spencer R Keyser
##
## Date Created: 2025-01-31
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
library(sf)
library(rjags)
library(MCMCvis)
library(coda)

## -------------------------------------------------------------

## -------------------------------------------------------------
##
## Begin Section: Load in data
##
## -------------------------------------------------------------

## MCMC Samples
load("R:/Users/skeyser/Postdoc/MSOM_Ragged_JAGS_Zout_95thresh.Rdata")
str(out3)
all3 <- as.matrix(out3)
rm(out3)
gc()

## Occupancy model data
load(here("./Data/JAGS_Data/MSOM_Ragged_2021_95cut.RData"))

## Package the Zs up
## Z matrix
nsite <- win.data.rag$nsite
nspec <- win.data.rag$nspec
nsamp <- dim(all3)[1]
z <- array(NA, dim = c(nsite, nspec, nsamp))
Jacc <- array(NA, dim = c(nsite, nspec, nsamp))

## Find indices for the colnames of the Z posterior
z.ind <- which(str_detect(colnames(all3), "z"))
z.min <- min(z.ind)
z.max <- max(z.ind)
for(j in 1:nsamp){
  cat(paste("\nMCMC sample", j, "\n"))
  z[,,j] <- all3[j, z.min:z.max]
}

str(z)

gc()

## For our analysis we can sample 1000 posterior draws
sample.draws <- sort(sample(1:dim(z)[3], 100, replace = F))
z.samp <- z[,,sample.draws]

## Species-specific responses
sp.index <- read.csv(here("Code/Occupancy_Modeling/SpeciesIndex_Filtered.csv"))
sp.index$Index <- as.character(sp.index$Index)

## Add species names to the mix
dimnames(z.samp) <- list(NULL, sp.index$Species, NULL)

## -------------------------------------------------------------
##
## Begin Section: ARU Metadata
##
## -------------------------------------------------------------

## Select the variables we are interested in for detection and occupancy
## Load in the bird data for 2021
load(here("./Data/Generated_DFs/Occ_Mod_Data/95_Thresh_Cutoff/2021_OccSppList.RData"))

## ADDED removing species and name change for PACFLY
sp.det.list <- sp.det.list[which(!names(sp.det.list) %in% c("Red-tailed Hawk",
                                                            "Osprey",
                                                            "Red-shouldered Hawk",
                                                            "Clark's Nutcracker"))]

names(sp.det.list)[which(names(sp.det.list) == "Pacific-slope Flycatcher")] <- "Western Flycatcher"

## Cell Unit mapping file
cu.map <- data.frame(ID = 1:length(sp.det.list[[1]]$Cell_Unit), Cell_Unit = sp.det.list[[1]]$Cell_Unit)
cu.map <- cu.map |>
  mutate(Cell_Unit = ifelse(
    stringr::str_detect(string = Cell_Unit, pattern = "C[0-9]{3}"),
    gsub(pattern = "(C)([0-9]{3})(_U[0-9]+)$", replacement = "\\10\\2\\3", x = Cell_Unit),
    Cell_Unit
  ))


## Load in the ARU meta data
aru_meta <- readr::read_csv(here("Data/ARU_120m.csv"))
aru_meta$Cell_Unit <- paste0(aru_meta$cell_id, "_", aru_meta$unit_numbe)

aru_filter <- aru_meta |> select(Cell_Unit, 
                                 utme, utmn,
                                 X, Y,
                                 ele = topo_elev, tpi = topo_tpi,
                                 tmax = tmx_bcm_mn, ppt = ppt_bcm_mn,
                                 fire1yr_lowmod_prop, fire1yr_high_prop, 
                                 fire2_5yr_lowmod_prop, fire2_5yr_high_prop,
                                 fire6_10yr_lowmod_prop, fire6_10yr_high_prop, 
                                 fire11_35yr_lowmod_prop, fire11_35yr_high_prop,
                                 stage = standage_f3_mn, cancov = cpycovr_f3_mn) |> 
  filter(Cell_Unit %in% cu.map$Cell_Unit) |> 
  arrange(match(Cell_Unit, cu.map$Cell_Unit)) |> 
  mutate(fire1_5yr_lowmod_prop = fire1yr_lowmod_prop + fire2_5yr_lowmod_prop,
         fire1_5yr_high_prop = fire1yr_high_prop + fire2_5yr_high_prop) |> 
  select(-fire1yr_high_prop, -fire1yr_lowmod_prop, -fire2_5yr_lowmod_prop, -fire2_5yr_high_prop) |> 
  select(Cell_Unit:ele, stage, cancov, 
         fire1_5yr_lowmod_prop, fire6_10yr_lowmod_prop, fire11_35yr_lowmod_prop,
         fire1_5yr_high_prop, fire6_10yr_high_prop, fire11_35yr_high_prop)


## Add site names to the Zmatrices
site.ids <- aru_filter$Cell_Unit
## Add species names to the mix
dimnames(z.samp) <- list(site.ids, sp.index$Species, NULL)

## Load the files for the fire only sites
fires <- read.csv(here("./Data/Generated_DFs/ARU2021_Creek_NCP_Filtered.csv"))

## Filter the Aru data by the sites in the fires DF
aru_fire <- aru_filter |> 
  filter(Cell_Unit %in% fires$cell_unit)

## Filter the Zsamp array by the matching cell_unit names
fire_match <- match(aru_fire$Cell_Unit, rownames(z.samp)) 
z.samp.fire <- z.samp[fire_match,,]

## Set the names back to the null
dimnames(z.samp.fire) <- NULL

## Package all the data as a list
OccData <- list(Z.posterior = z.samp.fire,
                Nsite = dim(z.samp.fire)[1],
                Nspec = dim(z.samp.fire)[2],
                Z.draws = dim(z.samp.fire)[3],
                SiteMeta = aru_fire,
                ZcolNames = sp.index$Species,
                ZrowNames = aru_fire$Cell_Unit)

var.to.keep <- "OccData"

rm(list = setdiff(ls(), var.to.keep))

save(OccData, file = here("./Data/JAGS_Data/Occ2GDM_Data_95thresh_Fires.Rdata"))
