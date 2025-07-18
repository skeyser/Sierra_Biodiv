## -------------------------------------------------------------
##
## Script name: GDM Ferrier Z Matrix for CHTC
##
## Script purpose: GDM for Sierra Bioacoustics Latent Z using
## traditional Ferrier GDM approach running on CHTC
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
library(data.table)
library(dplyr)
library(purrr)
library(tidyr)
library(sf)

## GDM
library(gdm)
library(adespatial)

## -------------------------------------------------------------

## Command line argument
args <- commandArgs(trailingOnly = TRUE)

## Load in a single Z-matrix
Z <- fread(file = args[1])
Z <- as.data.frame(Z)

## Take the environmental data
envTab <- fread(file = args[2])
envTab <- as.data.frame(envTab)
## -------------------------------------------------------------
##
## Begin Section: Format for GDM
##
## -------------------------------------------------------------

## Pluck one iteration of the posterior
rownames(Z) <- envTab$Cell_Unit
envTab$Cell_Unit <- 1:length(unique(envTab$Cell_Unit))

## Partition the beta diversity based on Jaccard-based Podani metrics
bpart <- adespatial::beta.div.comp(Z, coef = "J")

## Replacement
brep <- as.matrix(bpart$repl)
dimnames(brep) <- list(envTab$Cell_Unit, envTab$Cell_Unit)
brep <- cbind(Cell_Unit = envTab$Cell_Unit, brep)

## Richness
brich <- as.matrix(bpart$rich)
dimnames(brich) <- list(envTab$Cell_Unit, envTab$Cell_Unit)
brich <- cbind(Cell_Unit = envTab$Cell_Unit, brich)

## Total Beta
btotal <- as.matrix(bpart$D)
dimnames(btotal) <- list(envTab$Cell_Unit, envTab$Cell_Unit)
btotal <- cbind(Cell_Unit = envTab$Cell_Unit, btotal)

## Format the data for pairs - bioFormat = 3 accepts preformed D-mat
gdmTab.brep <- formatsitepair(bioData = brep,
                              bioFormat = 3,
                              XColumn = "Long",
                              YColumn = "Lat",
                              siteColumn = "Cell_Unit",
                              predData = envTab)

gdmTab.brich <- formatsitepair(bioData = brich,
                               bioFormat = 3,
                               XColumn = "Long",
                               YColumn = "Lat",
                               siteColumn = "Cell_Unit",
                               predData = envTab)

gdmTab.btotal <- formatsitepair(bioData = btotal,
                                bioFormat = 3,
                                XColumn = "Long",
                                YColumn = "Lat",
                                siteColumn = "Cell_Unit",
                                predData = envTab)


## Fit GDM for one Z matrix slice
gdm.fit.brep <- gdm(gdmTab.brep, geo = TRUE)
gdm.fit.brich <- gdm(gdmTab.brich, geo = TRUE)
gdm.fit.btotal <- gdm(gdmTab.btotal, geo = TRUE)

## Add in the variable importance
vi.brep <- gdm.varImp(gdmTab.brep,
                 predSelect = FALSE,
                 geo = T,
                 nPerm = 100,
                 parallel = FALSE)

vi.brich <- gdm.varImp(gdmTab.brich,
                      predSelect = FALSE,
                      geo = T,
                      nPerm = 100,
                      parallel = FALSE)

vi.btotal <- gdm.varImp(gdmTab.btotal,
                      predSelect = FALSE,
                      geo = T,
                      nPerm = 100,
                      parallel = FALSE)


## Package up the data in a list
## Make the list to hold model fits
gdm.fit.list <- tibble(
  beta_metric = c("Brep", "Brich", "Btotal"),
  models = vector("list", 3),
  varImp = vector("list", 3)) 


## Store variables in list for output
gdm.fit.list[gdm.fit.list$beta_metric == "Brep",]$models[[1]] <- gdm.fit.brep
gdm.fit.list[gdm.fit.list$beta_metric == "Brep",]$varImp[[1]] <- vi.brep

gdm.fit.list[gdm.fit.list$beta_metric == "Brich",]$models[[1]] <- gdm.fit.brich
gdm.fit.list[gdm.fit.list$beta_metric == "Brich",]$varImp[[1]] <- vi.brich

gdm.fit.list[gdm.fit.list$beta_metric == "Btotal",]$models[[1]] <- gdm.fit.btotal
gdm.fit.list[gdm.fit.list$beta_metric == "Btotal",]$varImp[[1]] <- vi.btotal

## Save the local output as an RDS object
saveRDS(gdm.fit.list, file = args[3])