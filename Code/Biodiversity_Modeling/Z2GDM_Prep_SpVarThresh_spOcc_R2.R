## -------------------------------------------------------------
##
## Script name: Data for GDM from spOccupancy R2
##
## Script purpose: Prepping output from the MSOM for GDM
## using variable species thresholds with min 97.5 Max prec.
## Editted for revisions for FEE.
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
library(spOccupancy)
library(MCMCvis)
library(coda)

## -------------------------------------------------------------

## -------------------------------------------------------------
##
## Begin Section: Load in data
##
## -------------------------------------------------------------

## Load in the spOccupancy model output
## We will use the full model given better WAIC
out <- readRDS(here("./Data/spOccupancy_Data/Model_sfMsPGOcc_2slf_R2.rds"))

## MCMC Samples
#z <- readRDS("R:/Users/skeyser/Postdoc/spOccupancy/Z_Matrix_sfMsPGOcc_2slf.rds")
z <- out$z.samples
str(z)
z <- aperm(z, perm = c(3,2,1))

## Occupancy model data
OccData <- readRDS(here("./Data/SpOccupancy_Data/SpOccData_SpeciesThresh_975minMaxPrex_NewVars_2024Revalidation.rds"))
str(OccData)
## Package the Zs up
## Z matrix
nsite <- ncol(OccData$y)
nspec <- nrow(OccData$y)
nsamp <- dim(z)[3]

## For our analysis we can sample 1000 posterior draws
sample.draws <- sort(sample(1:dim(z)[3], 500, replace = F))
z.samp <- z[,,sample.draws]

## Species-specific responses
## Set a common species first
sp.trait <- readRDS(here("./Data/TrialSpeciesOrder.rds"))

sp.add <- sp.trait[sp.trait %in% unlist(dimnames(OccData$y)[1])]


sp.reorder <- unique(c("Brown Creeper",
                       "Acorn Woodpecker",
                       "Lazuli Bunting",
                       sp.add,
                       "Sphyrapicus spp.",
                       "Cassin's Vireo",
                       "Pacific-slope Flycatcher"))

sp.index <- data.frame(Species = sp.reorder, Index = 1:length(sp.reorder))
sp.index$Index <- as.character(sp.index$Index)

## Add species names to the mix
for(i in 1:dim(z.samp)[3]){
  colnames(z.samp[,,i]) <- sp.index$Species
}
## -------------------------------------------------------------
##
## Begin Section: ARU Metadata
##
## -------------------------------------------------------------

## Select the variables we are interested in for detection and occupancy
## Load in the bird data for 2021
load(here("./Data/Generated_DFs/Occ_Mod_Data/Thresh_By_Species/2021_99Conf_OccSppList.RData"))

## ADDED removing species and name change for PACFLY
sp.det.list <- sp.det.list[which(!names(sp.det.list) %in% c("Red-tailed Hawk",
                                                            "Osprey",
                                                            "Red-shouldered Hawk",
                                                            "Clark's Nutcracker",
                                                            "American Kestrel"
                                                            #"Wild Turkey",
                                                            #"Sooty Grouse",
                                                            #"Lawrence's Goldfinch",
                                                            #"California Thrasher"
))]

names(sp.det.list)[which(names(sp.det.list) == "Pacific-slope Flycatcher")] <- "Western Flycatcher"

combine_species <- function(species_list, species_to_combine, new_name, date_cols) {
  # Start with first species
  combined_df <- species_list[[species_to_combine[1]]]
  
  # Sum values across species
  for(sp in species_to_combine[-1]) {
    combined_df[, date_cols] <- combined_df[, date_cols] + species_list[[sp]][, date_cols]
  }
  
  # Binarize: convert all values > 0 to 1
  combined_df <- combined_df %>%
    mutate(across(all_of(date_cols), ~as.integer(. > 0)))
  
  # Remove original species and add combined
  species_list[species_to_combine] <- NULL
  species_list[[new_name]] <- combined_df
  
  return(species_list)
}

## Combine Sapsuckers
sp.det.list <- combine_species(species_list = sp.det.list, 
                               species_to_combine = c("Red-breasted Sapsucker", "Red-naped Sapsucker", "Williamson's Sapsucker"), 
                               new_name = "Sphyrapicus spp.",
                               date_cols = colnames(sp.det.list[[1]])[str_detect(colnames(sp.det.list[[1]]), "Cell_Unit", negate = T)])

## Combine vireos
sp.det.list <- combine_species(species_list = sp.det.list, 
                               species_to_combine = c("Cassin's Vireo", "Plumbeous Vireo"), 
                               new_name = "Vireo spp.",
                               date_cols = colnames(sp.det.list[[1]])[str_detect(colnames(sp.det.list[[1]]), "Cell_Unit", negate = T)])

names(sp.det.list)
str(sp.det.list)

## Cell Unit mapping file
cu.map <- data.frame(ID = 1:length(sp.det.list[[1]]$Cell_Unit), Cell_Unit = sp.det.list[[1]]$Cell_Unit)
cu.map <- cu.map |>
  mutate(Cell_Unit = ifelse(
    stringr::str_detect(string = Cell_Unit, pattern = "C[0-9]{3}"),
    gsub(pattern = "(C)([0-9]{3})(_U[0-9]+)$", replacement = "\\10\\2\\3", x = Cell_Unit),
    Cell_Unit
  ))


## Load in the ARU meta data
aru_meta <- readr::read_csv(here("Data/ARU_120m_New.csv"))
aru_meta$Cell_Unit <- paste0(aru_meta$cell_id, "_", aru_meta$unit_numbe)

aru_filter <- aru_meta |> select(Cell_Unit, 
                                 utme, utmn,
                                 X, Y,
                                 ele = topo_elev, 
                                 tpi = topo_tpi,
                                 tmax = tmx_bcm_mn, 
                                 ppt = ppt_bcm_mn,
                                 fire1yr_lowmod_prop, 
                                 fire1yr_high_prop, 
                                 fire2_5yr_lowmod_prop, 
                                 fire2_5yr_high_prop,
                                 fire1_5yr_lowmod_prop = fire5yr_lowmod_prop,
                                 fire1_5yr_high_prop = fire5yr_high_prop,
                                 fire6_10yr_lowmod_prop, 
                                 fire6_10yr_high_prop, 
                                 fire11_35yr_lowmod_prop, 
                                 fire11_35yr_high_prop,
                                 fire1yr_cbi_mn, 
                                 fire2_5yr_cbi_mn,
                                 fire1_5yr_cbi_mn = fire5yr_cbi_mn,
                                 fire6_10yr_cbi_mn, 
                                 fire11_35yr_cbi_mn,
                                 stage = standage_f3_mn, 
                                 cancov = cc_cfo_mn,
                                 canht = ch_cfo_mn
) |> 
  filter(Cell_Unit %in% cu.map$Cell_Unit) |> 
  arrange(match(Cell_Unit, cu.map$Cell_Unit)) |> 
  # mutate(fire1_5yr_lowmod_prop = fire1yr_lowmod_prop + fire2_5yr_lowmod_prop,
  #        fire1_5yr_high_prop = fire1yr_high_prop + fire2_5yr_high_prop,
  #        fire1_5yr_cbi_mn = (fire1yr_cbi_mn + fire2_5yr_cbi_mn) / 2) |> 
  select(-fire1yr_high_prop, 
         -fire1yr_lowmod_prop, 
         -fire2_5yr_lowmod_prop, 
         -fire2_5yr_high_prop,
         -fire1yr_cbi_mn,
         -fire2_5yr_cbi_mn) |> 
  select(Cell_Unit:ele, ppt, stage, cancov, canht,
         fire1_5yr_lowmod_prop, fire6_10yr_lowmod_prop, fire11_35yr_lowmod_prop,
         fire1_5yr_high_prop, fire6_10yr_high_prop, fire11_35yr_high_prop,
         fire1_5yr_cbi_mn, fire6_10yr_cbi_mn, fire11_35yr_cbi_mn)

## Create orthogonal predictor by lm
ch_res <- lm(canht ~ cancov, aru_filter)$residual

aru_filter <- aru_filter |> mutate(ch_res = ch_res)

## Package all the data as a list
OccData <- list(Z.posterior = z.samp,
                Nsite = dim(z.samp)[1],
                Nspec = dim(z.samp)[2],
                Z.draws = dim(z.samp)[3],
                SiteMeta = aru_filter,
                ZcolNames = sp.index$Species,
                ZrowNames = aru_filter$Cell_Unit)

var.to.keep <- "OccData"

rm(list = setdiff(ls(), var.to.keep))

save(OccData, file = here("./Data/SpOccupancy_Data/Occ2GDM_Data_SpThresh_975minMaxPrec_spOcc_2SLF.Rdata"))