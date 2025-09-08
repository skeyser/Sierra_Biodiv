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

## If we want to remove the forest preds we just need to strip them from the data
envTab <- envTab |> select(-cancov, -ch_res)

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
## Fit GDM for one Z matrix slice with error handling
gdm.fit.brep <- tryCatch({
  fit <- gdm(gdmTab.brep, geo = TRUE)
  if(all(fit$splines == 0)) {
    message("Null model fit for brep")
    NULL
  } else {
    fit
  }
}, error = function(e) {
  warning(paste("Error in brep GDM fit:", e))
  NULL
})

gdm.fit.brich <- tryCatch({
  fit <- gdm(gdmTab.brich, geo = TRUE)
  if(all(fit$splines == 0)) {
    message("Null model fit for brich")
    NULL
  } else {
    fit
  }
}, error = function(e) {
  warning(paste("Error in brich GDM fit:", e))
  NULL
})

gdm.fit.btotal <- tryCatch({
  fit <- gdm(gdmTab.btotal, geo = TRUE)
  if(all(fit$splines == 0)) {
    message("Null model fit for btotal")
    NULL
  } else {
    fit
  }
}, error = function(e) {
  warning(paste("Error in btotal GDM fit:", e))
  NULL
})

## Variable importance only if model isn't null
vi.brep <- if(!is.null(gdm.fit.brep)) {
  gdm.varImp(gdmTab.brep,
             predSelect = FALSE,
             geo = T,
             nPerm = 50,
             sampleSites = 1,
             sampleSitePairs = 0.2,
             parallel = FALSE)
} else NULL

vi.brich <- if(!is.null(gdm.fit.brich)) {
  gdm.varImp(gdmTab.brich,
             predSelect = FALSE,
             geo = T,
             nPerm = 50,
             sampleSites = 1,
             sampleSitePairs = 0.2,
             parallel = FALSE)
} else NULL

vi.btotal <- if(!is.null(gdm.fit.btotal)) {
  gdm.varImp(gdmTab.btotal,
             predSelect = FALSE,
             geo = T,
             nPerm = 50,
             sampleSites = 1,
             sampleSitePairs = 0.2,
             parallel = FALSE)
} else NULL

## Package up the data in a list
## Make the list to hold model fits
gdm.fit.list <- tibble(
  beta_metric = c("Brep", "Brich", "Btotal"),
  models = vector("list", 3),
  varImp = vector("list", 3)) 


## Store variables in list for output
if(!is.null(gdm.fit.brep)){
  gdm.fit.list[gdm.fit.list$beta_metric == "Brep",]$models[[1]] <- gdm.fit.brep
  gdm.fit.list[gdm.fit.list$beta_metric == "Brep",]$varImp[[1]] <- vi.brep
}

if(!is.null(gdm.fit.brich)){
  gdm.fit.list[gdm.fit.list$beta_metric == "Brich",]$models[[1]] <- gdm.fit.brich
  gdm.fit.list[gdm.fit.list$beta_metric == "Brich",]$varImp[[1]] <- vi.brich
}

if(!is.null(gdm.fit.btotal)){
  gdm.fit.list[gdm.fit.list$beta_metric == "Btotal",]$models[[1]] <- gdm.fit.btotal
  gdm.fit.list[gdm.fit.list$beta_metric == "Btotal",]$varImp[[1]] <- vi.btotal
}

## Save the local output as an RDS object
saveRDS(gdm.fit.list, file = args[3])