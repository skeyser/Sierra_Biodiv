## -------------------------------------------------------------
##
## Script name: Sierra Function Lib
##
## Script purpose: Sourcing custom functions written for processing
## and analyzing Sierra Bioacoustics data.
##
## Author: Spencer R Keyser
##
## Date Created: 2024-10-18
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

## -------------------------------------------------------------


## -------------------------------------------------------------
##
## Begin Section: ARU species detection function
##
## -------------------------------------------------------------

## Add function to format the species data according to thresholds and other parameters of
## interest.
ac_det_filter <- function(d,
                          d_thresh,
                          thresh_scale = "Conf",
                          thresh_cut = "99",
                          time_format = "ymd",
                          species,
                          no_dets = 1,
                          binary = T,
                          date_range = NULL) {
  
  ## Retrieve the threshold from the species of interest
  sp.thresh <- d_thresh[d_thresh$species == species, ]
  
  ## Check the threshold type
  if (thresh_scale == "Conf") {
    tmp.cols <- which(stringr::str_detect(colnames(sp.thresh), pattern = "conf"))
    tmp.cols <- sp.thresh[, tmp.cols]
    if (all(stringr::str_detect(colnames(tmp.cols), thresh_cut))) {
      print("Threshold not named in threshold dataframe. Check column names and 'thresh_cut'.")
      tmp.thresh <- NA
    } else {
      tmp.thresh <- pull(tmp.cols[which(stringr::str_detect(colnames(tmp.cols), thresh_cut))])
    }
  } else if (thresh_scale == "Raw") {
    tmp.cols <- which(
      stringr::str_detect(colnames(sp.thresh), pattern = "(cutoff)([0-9]+)(.r$)|(cutoff)([0-9]+)(raw$)")
    )
    tmp.cols <- sp.thresh[, tmp.cols]
    if (all(stringr::str_detect(colnames(tmp.cols), thresh_cut))) {
      print("Threshold not named in threshold dataframe. Check column names and 'thresh_cut'.")
      tmp.thresh <- NA
    } else {
      tmp.thresh <- pull(tmp.cols[which(stringr::str_detect(colnames(tmp.cols), thresh_cut))])
    }
  }
  
  ## Check the acoustic data to make sure the threshold is applied
  ## Retrieve the columns names
  colnames(d)
  
  ## Acoustic data will likely contain date specific capture histories...these are the primary interest
  if(time_format == "ymd"){
    date.cols <- colnames(d)[which(
      stringr::str_detect(colnames(d), pattern = "(\\d{4})[[:punct:]]?(\\d{2})[[:punct:]]?(\\d{2})")
    )]
    other.cols <- colnames(d)[which(
      !stringr::str_detect(colnames(d), pattern = "(\\d{4})[[:punct:]]?(\\d{2})[[:punct:]]?(\\d{2})")
    )] 
  } else if(time_format %in% c("dmy", "mdy")){
    date.cols <- colnames(d)[which(
      stringr::str_detect(colnames(d), pattern = "(\\d{2})[[:punct:]]?(\\d{2})[[:punct:]]?(\\d{4})")
    )]
    other.cols <- colnames(d)[which(
      !stringr::str_detect(colnames(d), pattern = "(\\d{2})[[:punct:]]?(\\d{2})[[:punct:]]?(\\d{4})")
    )] 
  } else if(is.na(time_format)){
    stop("Time format not specified.")
  }
  
  ## We want to operate on the data within the date columns
  ## Find the minimum value for the df
  min.val <- suppressWarnings(min(apply(as.matrix(d[,date.cols]), 1, function(x) {min(x[x>0], na.rm = T)})))
  
  ## Replace the "." notation with a "0"
  d[d == "."] <- "0"
  
  ## Date range restriction
  if (!is.null(date_range)) {
    if (time_format == "ymd") {
      ## Take the start and end date
      s.date <- lubridate::ymd(date_range[1])
      e.date <- lubridate::ymd(date_range[2])
      
      # Convert the date columns to a format
      samp.dates <- lubridate::ymd(date.cols)
      
    } else if (time_format == "mdy") {
      ## Take the start and end date
      s.date <- lubridate::mdy(date_range[1])
      e.date <- lubridate::mdy(date_range[2])
      
      # Convert the date columns to a format
      samp.dates <- lubridate::mdy(date.cols)
      
    } else if (time_format == "dmy") {
      ## Take the start and end date
      s.date <- lubridate::dmy(date_range[1])
      e.date <- lubridate::dmy(date_range[2])
      
      # Convert the date columns to a format
      samp.dates <- lubridate::dmy(date.cols)
    }
    
    ## Check successful
    if (class(samp.dates) != "Date") {
      stop("Failed to convert to date format, but date range is provided.")
    }
    
    ## Subset the dates from the entire date list
    date.cols <- as.character(samp.dates[samp.dates >= s.date &
                                           samp.dates <= e.date])
    
    ## subset the columns of interest and return 'd'
    d <- d |>
      select(all_of(c(other.cols, date.cols)))
  }
  
  ## If binary
  if(binary) {
    d.thresh <- d |>
      mutate(across(all_of(date.cols), ~ as.numeric(.))) |> 
      mutate(across(all_of(date.cols), ~ ifelse(. >= tmp.thresh, 1, 0)))
  } else {
    d.thresh <- d |>
      mutate(across(all_of(date.cols), ~ as.numeric(.))) |> 
      mutate(across(all_of(date.cols), ~ ifelse(. >= tmp.thresh, ., 0)))
  }
  
  ## If the number of detections in the season matters
  if (is.numeric(no_dets)) {
    if (!binary) {
      tmp.detect <- d.thresh |>
        mutate(NDets = rowSums(across(all_of(date.cols)) > 0)) |> #sum logical > 0
        mutate(across(all_of(date.cols), ~ if_else(NDets >= no_dets, ., 0))) |>
        mutate(NDetsNew = rowSums(across(all_of(date.cols)) > 0))
    }
    
    if (binary) {
      tmp.detect <- d.thresh |>
        mutate(NDets = rowSums(across(all_of(date.cols)))) |>
        mutate(across(all_of(date.cols), ~ if_else(NDets >= no_dets, ., 0))) |>
        mutate(NDetsNew = rowSums(across(all_of(date.cols))))
    }
    
    ## Check the columns that changed
    ## If the statement below is true...this means the sum of all the scores below were
    ## less than the number of detections we need at the site over the season and these
    ## rows were made to be entirely nondetections
    all(dplyr::setdiff(tmp.detect$NDets, tmp.detect$NDetsNew) < (no_dets))
    
    ## Drop the helper columns and return d.thresh as the fixed df
    d.thresh <- tmp.detect |>
      dplyr::select(all_of(c(other.cols, date.cols)))
    
    ## Remove superfluous DFs made
    rm(tmp.detect)
  }
  return(d.thresh)
  
} #function

## -------------------------------------------------------------
##
## End Section: ARU species detection function
##
## -------------------------------------------------------------