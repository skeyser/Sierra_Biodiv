## -------------------------------------------------------------
##
## Script name: CBI Fire Extraction
##
## Script purpose: Extract and summarize fire variables for
## ARU locations across the Sierra Nevada
##
## Author: Spencer R Keyser
##
## Date Created: 2024-10-17
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
library(terra)
library(sf)
## -------------------------------------------------------------

## Pull in the ARU meta data for point locations

## Take a look at CBI data
cbi <- rast("c:/Users/srk252/Documents/data_for_spencer/cbi_sierra_cat_rasters/cbi_cat_2021.tif")
plot(cbi)
