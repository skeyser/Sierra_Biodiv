## -------------------------------------------------------------
##
## Script name: ARU Species Detection Pre-processing
##
## Script purpose: Format and process species detection histories
## from ARUs and link together with spatial and environmental metadata.
##
## Author: Spencer R Keyser
##
## Date Created: 2024-10-09
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
# Package Loading
library(here)
library(dplyr)
library(purrr)
library(ggplot2)
library(stringr)
library(CAbioacoustics)
library(sf)
library(mapview)

## -------------------------------------------------------------

## -------------------------------------------------------------
##
## Begin Section: Thresholds file exploration
##
## -------------------------------------------------------------
## *************************************************************
##
## Section Notes:
## Variable types
## 1. Species {character} = Species Common Name
## 2. keep {numeric} = Binary indicator variable
## 3. BirdNET_code {numeric} = Species' code names from BirdNET
## 4. intercept.c {numeric}
## 5. beta.c {numeric}
## 6. intercept.r {numeric}
## 7. beta.r {numeric}
## 8. cutoff85.r {numeric}
## 9. cutoff90.r {numeric}
## 10. cutoff95.r {numeric}
## 11. cutoff975.r {numeric}
## 12. cutoff99.r {numeric}
## 13. cutoff90.r_conf {numeric}
## 14. cutoff95.r_conf {numeric}
## 15. cutoff99.r_conf {numeric}
##
## *************************************************************

## Read thresholds in
thresh <- readr::read_csv(here("Data/Thresholds_2021_20230309.csv"))
glimpse(thresh)

## DF Dimensions
dim(thresh)

## How many species
length(unique(thresh$species))
print(thresh$species)

## Names of the columns
colnames(thresh)

## Breakdown for "keep"
table(thresh$keep)

## Overall summarization
summary(thresh)

## -------------------------------------------------------------
##
## Begin Section: Processing species detections
##
## -------------------------------------------------------------
## *************************************************************
##
## Section Notes: Source ac_det_filter() from Sierra_functions.R
## 
## To-do:
## 1. Streamline the multi-year species data processing
## 2. Unify with spatial coordinates from CAbioacoustics
## 3. Write yearly files for community data (summarized; naive biodiv. models)
## 4. Write yearly files for community data (date-detection; occ models)
##
## *************************************************************


## Source the ac_det_filter function
source(here("./Code/Acoustic_Data_Prep/Sierra_functions.R"))

## threshold file
thresh <- readr::read_csv("./Data/Thresholds_2021_20230309.csv")
thresh <- readr::read_csv("./Data/Thresholds_2021_20230309_BestThreshold_975min.csv")

## Run the function
aru_det_file_gen(det_dir = "C:/Users/srk252/Documents/Rprojs/Sierra_Biodiv/Data/Detections_By_Species/",
                 det_years = c("2021"),
                 seas_format = F,
                 seas_outdir = "C:/Users/srk252/Documents/Rprojs/Sierra_Biodiv/Data/Generated_DFs/Seasonal_Summaries/",
                 occ_format = T,
                 occ_outdir = "C:/Users/srk252/Documents/Rprojs/Sierra_Biodiv/Data/Generated_DFs/Occ_Mod_Data/Thresh_By_Species/",
                 eff_file = T,
                 coord_link = T,
                 d_thresh = thresh,
                 thresh_scale = "Conf",
                 thresh_transform = TRUE,
                 thresh_trans_dir = "conf2logit",
                 thresh_cut = "99",
                 species_thresh_cut = "BestThresh",
                 time_format = "ymd",
                 no_dets = 2,
                 binary = T,
                 date_range = c("2021-06-01", "2021-06-30"),
                 eff_site_name = "Cell_U",
                 eff_filter = 10,
                 verbose = F)

