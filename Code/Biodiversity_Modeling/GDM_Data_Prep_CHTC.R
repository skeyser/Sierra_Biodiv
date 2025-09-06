## -------------------------------------------------------------
##
## Script name: Prep Data for GDM runs on CHTC
##
## Script purpose: Create individual matrices for running on
## CHTC using the 1000 posterior samples of the Z-matrix.
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
library(purrr)
library(tidyr)
library(sf)

## GDM
library(gdm)
library(adespatial)

## -------------------------------------------------------------

## Load in the Z-matrices
load(here("./Data/JAGS_Data/Occ2GDM_Data_SpThresh_975minMaxPrec.Rdata"))
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

aru_meta <- aru_meta |> 
  select(all_of(v.keep))


## -------------------------------------------------------------
##
## Begin Section: Format for GDM
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

## Generate two input types
## 1. envTab which is the same
data.table::fwrite(envTab, file = here("./Data/Z_spOcc_Post_Samples_CHTC/envTab.csv"))

## 2. 1000 Z matrices
dim(Z)
for(i in 1:dim(Z)[3]){
  z.tmp <- Z[,,i]
  data.table::fwrite(z.tmp, file = paste0(here("./Data/Z_spOcc_Post_Samples_CHTC/"), "/Z_sample_", i, ".csv"))
}

## Generate a quick reference file for the looping with input and output
ZoutSort <- data.frame(Input = paste0(rep("Z_sample", times = dim(Z)[3]), "_", seq(1, dim(Z)[3]), ".csv"),
                       Output = paste0(rep("Output_GDM", times = dim(Z)[3]), "_", seq(1, dim(Z)[3]), ".rds"))
colnames(ZoutSort) <- NULL

data.table::fwrite(ZoutSort, file = here("./Data/Z_spOcc_Post_Samples_CHTC/GDM_Sort.csv"))
