## -------------------------------------------------------------
##
## Script name: CAbioacoustics Data Exploration
##
## Script purpose: Connect to and explore CAbioacoustics package
##
## Author: Spencer R Keyser
##
## Date Created: 2024-09-30
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
library(dbplyr)
library(CAbioacoustics)
library(sf)
library(terra)

## -------------------------------------------------------------

## -------------------------------------------------------------
##
## Begin Section: Conncec to the DB
##
## -------------------------------------------------------------
## *************************************************************
##
## Section Notes: This should work because we have already
## conncected a stashed the credentials in the keyring of the 
## computer.
##
## *************************************************************
cb_connect_db()

## -------------------------------------------------------------
##
## Begin Section: Explore the database
##
## -------------------------------------------------------------

## View the tables housed in the DB
DBI::dbListTables(conn = conn)

## View a single file
tbl(conn, "acoustic_efforts")

## Retrieve the Hex grid
hex <- cb_get_spatial(layer_name = "sierra_hexes")

## Retrieve the study area
region <- cb_get_spatial(layer_name = "sierra_study_area")

## Loaded hexes as an SF object
class(hex)

## What is the information within hex
names(hex)
head(hex)

## Plot the hexes
## Color code by elevation
ggplot(data = hex) +
  geom_sf(aes(color = elev_mean))
ggplot(data = hex) +
  geom_sf(aes(color = ownership))

## Plot the study region for the Sierra
ggplot(data = region) +
  geom_sf()

## How many unique ARU codes are there
length(unique(hex$cell_id))
length(unique(hex$ownership))

## Find the area of a single hex
st_area(hex[1, "geometry"])
st_area(region)

## Acoustic field visits are used for the spatial reference
