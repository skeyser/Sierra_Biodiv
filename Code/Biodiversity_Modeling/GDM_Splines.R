## -------------------------------------------------------------
##
## Script name: GDM_Prep Script
##
## Script purpose:
##
## Author: Spencer R Keyser
##
## Date Created: 2025-02-05
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
library(splines2)
## -------------------------------------------------------------

## Load the predictor data for the I-spline generation
aru_meta <- readr::read_csv(here("./Data/ARU_Meta_120_FilteredOcc.csv"))

## Select the predictors for the GDM
## *************************************************************
##
## Section Notes: Predictor Selection
## 1. Elevation
## 2. Latitude
## 3. Prop High Sev Fire: 1-5, 6-10, 11-35
## 4. Prop Mod/Low Sev Fire: 1-5, 6-10, 11-35
## 5. Stand Age
## 6. Canopy Cover
## 
## *************************************************************

aru_sub <- aru_meta |> 
  select(Cell_Unit, 
         utme, utmn,
         topo_elev,
         standage_f3_mn,
         cpycovr_f3_mn,
         contains("high_prop"),
         contains("lowmod_prop")
         ) |> 
  mutate(fire1_5yr_high_prop = fire1yr_high_prop + fire2_5yr_high_prop,
         fire1_5yr_lowmod_prop = fire1yr_high_prop + fire2_5yr_lowmod_prop) |> 
  select(Cell_Unit:cpycovr_f3_mn,
         fire1_5yr_high_prop,
         fire6_10yr_high_prop,
         fire11_35yr_high_prop,
         fire1_5yr_lowmod_prop,
         fire6_10yr_lowmod_prop,
         fire11_35yr_lowmod_prop)

# Function to create I-splines
create.isplines <- function(x, n_knots = 3, degree = 3) {
  library(splines2)
  
  # Create knots
  knots <- quantile(x, probs = seq(0, 1, length.out = n_knots + 2))[-c(1, n_knots + 2)]
  
  # Generate I-splines
  iSplines <- iSpline(x, knots = knots, 
                      degree = degree, 
                      intercept = TRUE)
  
  return(iSplines)
}


## Apply the iSplines to the dataframe


x <- unlist(as.vector(as.data.frame(aru_sub[,3])))
str(x)
isplines <- create.isplines(x)

isplines <- create.isplines(x = x)

matplot(x, isplines, type = "l", 
        main = "I-spline Basis Functions",
        xlab = "x", ylab = "I(x)")



