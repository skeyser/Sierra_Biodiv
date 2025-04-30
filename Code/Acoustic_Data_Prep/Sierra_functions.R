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

## *************************************************************
##
## Section Notes:
## The following data files correspond to Step3 in the BirdNET.
## The purpose the following block should be to load the Step3 data and create 
## the detection history across the secondary sampling periods. We can write a small
## function to accomplish this based on some modifications we want.
##
##
## To-do: Build in functionality for:
## 1. Filtering by effort
## 2. Handling multiple years of data
##
## *************************************************************

## Add function to format the species data according to thresholds and other parameters of
## interest.
ac_det_filter <- function(d,
                          d_thresh,
                          thresh_scale = "Conf",
                          species_thresh_cut = NA,
                          thresh_transform = FALSE,
                          thresh_trans_dir = "conf2logit",
                          thresh_cut = "99",
                          time_format = "ymd",
                          species,
                          no_dets = 1,
                          binary = T,
                          date_range = NULL,
                          eff = NULL,
                          eff_site_name = NULL,
                          eff_summary_cols = NULL,
                          eff_filter = 1,
                          verbose = T){
  
  ## Retrieve the threshold from the species of interest
  sp.thresh <- d_thresh[d_thresh$species == species, ]
  
  ## Add a component for species-specific thresholds
  ## species_thresh_cut should be a column in a DF with the best threshold
  ## matching d_thresh and thresh_scale
  # if(!is.null(species_thresh_cut)){
  #   if(any(str_detect(colnames(sp.thresh), species_thresh_cut))){
  #     best_thresh <- sp.thresh[, str_detect(colnames(sp.thresh), species_thresh_cut)]
  #     
  #     if(thresh_scale == "Conf" | (thresh_transform & thresh_trans_dir == "logit2conf")){
  #       new_thresh_cut <- as.character(pull(best_thresh |> select(contains("conf"))))
  #       # Only update thresh_cut if new value is not NA
  #       if(!is.na(new_thresh_cut)) thresh_cut <- new_thresh_cut
  #     } else if(thresh_scale == "Raw" | (thresh_transform & thresh_trans_dir == "conf2logit")){
  #       new_thresh_cut <- as.character(pull(best_thresh[,str_detect(colnames(best_thresh), "\\.r$")]))
  #     }
  #     
  #     if(!is.na(new_thresh_cut)) {
  #       thresh_cut <- new_thresh_cut
  #       print(paste("Best threshold option selected. Optimal threshold for", species, "is", thresh_cut))
  #     } else {
  #       print(paste("Warning: Threshold calculation produced NA for", species, "using function default value:", thresh_cut))
  #     }
  #   }
  #   thresh_cut_meta <- data.frame(Default_Thresh_Scale = thresh_scale,
  #                                 New_Thresh_Cut = thresh_cut,
  #                                 Thresh_Transform = thresh_transform,
  #                                 Thresh_Tran_Direction = thresh_trans_dir)
  # }
  
  if(!is.na(species_thresh_cut)){
    if(any(str_detect(colnames(sp.thresh), species_thresh_cut))){
      best_thresh <- sp.thresh[, str_detect(colnames(sp.thresh), species_thresh_cut)]
      
      # Scenario 1: thresh_scale is "Conf" and no transformation
      if(thresh_scale == "Conf" && !thresh_transform) {
        new_thresh_cut <- as.character(pull(best_thresh |> select(contains("conf"))))
        flag <- "No transformation. Taking confidence score threshold."
      }
      # Scenario 2: thresh_scale is "Raw" and no transformation
      else if(thresh_scale == "Raw" && !thresh_transform) {
        new_thresh_cut <- as.character(pull(best_thresh[,str_detect(colnames(best_thresh), "\\.r$")]))
        flag <- "No transformation. Taking raw score threshold."
      }
      # Scenario 4: thresh_scale is "Conf" and transformation from conf to logit
      else if(thresh_scale == "Conf" && thresh_transform && thresh_trans_dir == "conf2logit") {
        new_thresh_cut <- as.character(pull(best_thresh[,str_detect(colnames(best_thresh), "\\.r$")]))
        flag <- "Transformation conf2logit. Taking raw score threshold."
      }
      # Scenario 4: thresh_scale is "Raw" and transformation from logit to conf
      else if(thresh_scale == "Raw" && thresh_transform && thresh_trans_dir == "logit2conf") {
        new_thresh_cut <- as.character(pull(best_thresh |> select(contains("conf"))))
        flag <- "Transformation logit2conf. Taking confidence score threshold."
      }
      
      if(!is.na(new_thresh_cut)) {
        thresh_cut <- new_thresh_cut
        message(paste("Best threshold option selected. Optimal threshold for", species, "is", thresh_cut))
        flag2 <- "Species Thresh"
      } else {
        message(paste("Warning: Threshold calculation produced NA for", species, "using function default value:", thresh_cut))
        flag2 <- "Default Thresh"
      }
    }
    thresh_cut_meta <- data.frame(Species = species,
                                  Default_Thresh_Scale = thresh_scale,
                                  Thresh_Used = thresh_cut,
                                  Thresh_Transform = thresh_transform,
                                  Thresh_Tran_Direction = thresh_trans_dir,
                                  Flag = flag,
                                  Flag2 = flag2)
  } else {thresh_cut_meta <- NA}
  
  ## Check the threshold type
  ## Add in a line for direct transformation of the confidence scores
  ## to logit
  if (thresh_scale == "Conf" & isFALSE(thresh_transform)) {
    tmp.cols <- which(stringr::str_detect(colnames(sp.thresh), pattern = "conf"))
    tmp.cols <- sp.thresh[, tmp.cols]
    if (all(stringr::str_detect(colnames(tmp.cols), thresh_cut))) {
      warning("Threshold not named in threshold dataframe. Check column names and 'thresh_cut'.")
      tmp.thresh <- NA
    } else {
      tmp.thresh <- pull(tmp.cols[which(stringr::str_detect(colnames(tmp.cols), thresh_cut))])
    }
  } else if (thresh_scale == "Raw" & isFALSE(thresh_transform)) {
    tmp.cols <- which(
      stringr::str_detect(colnames(sp.thresh), pattern = "(cutoff)([0-9]+)(.r$)|(cutoff)([0-9]+)(raw$)")
    )
    tmp.cols <- sp.thresh[, tmp.cols]
    if (all(stringr::str_detect(colnames(tmp.cols), thresh_cut))) {
      warning("Threshold not named in threshold dataframe. Check column names and 'thresh_cut'.")
      tmp.thresh <- NA
    } else {
      tmp.thresh <- pull(tmp.cols[which(stringr::str_detect(colnames(tmp.cols), thresh_cut))])
    }
  } else if (thresh_scale == "Conf" & isTRUE(thresh_transform)) {
    tmp.cols <- which(
      stringr::str_detect(colnames(sp.thresh), pattern = "(cutoff)")
    )
    tmp.cols <- sp.thresh[, tmp.cols]
    if (all(stringr::str_detect(colnames(tmp.cols), thresh_cut))) {
      warning("Threshold not named in threshold dataframe. Check column names and 'thresh_cut'.")
      tmp.thresh <- NA
    } else {
      tmp.thresh <- tmp.cols[which(stringr::str_detect(colnames(tmp.cols), thresh_cut))]
    }
  } else if (thresh_scale == "Raw" & isTRUE(thresh_transform)) {
    tmp.cols <- which(
      stringr::str_detect(colnames(sp.thresh), pattern = "(cutoff)")
    )
    tmp.cols <- sp.thresh[, tmp.cols]
    if (all(stringr::str_detect(colnames(tmp.cols), thresh_cut))) {
      warning("Threshold not named in threshold dataframe. Check column names and 'thresh_cut'.")
      tmp.thresh <- NA
    } else {
      tmp.thresh <- tmp.cols[which(stringr::str_detect(colnames(tmp.cols), thresh_cut))]
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
      
    } else if(!time_format %in% c("ymd", "mdy", "dmy")){
      stop("Date format does not match accepted. 
           Please use 'ymd', 'mdy', or 'dmy'")
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
  
  ## Add in the transformation
  if(thresh_transform){
    conf2logit <- function(x) { log(x / (1-x)) }
    logit2conf <- function(x) { 1/(1 + exp(-x)) }
    if(thresh_trans_dir == "conf2logit"){
      d <- d |>
        mutate(across(.cols = all_of(date.cols), as.numeric)) |>
        mutate(across(.cols = all_of(date.cols), conf2logit))
      
      tmp.thresh <- tmp.thresh |> select(!contains("_conf")) |> pull()
      
    } else if(thresh_trans_dir == "logit2conf"){
      d <- d |>
        mutate(across(.cols = all_of(date.cols), as.numeric)) |>
        mutate(across(.cols = all_of(date.cols), logit2conf))
      
      tmp.thresh <- tmp.thresh |> select(contains("_conf")) |> pull()
    }
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
  }
  
  ## Bring in the effort variables
  if(is.null(eff)){
    if(verbose == T){cat("Missing effort covariates.\nData is not filtered for ARU effort.")}
  } else {
    if(verbose == T){cat("Effort file provided.\nChecking effort covariates.")}
  }
  
  if(!is.null(eff) & !is.null(eff_summary_cols)){
    if(verbose == T){cat("Effort file provided.\nEffort summary columns provided.")}
    if(as.numeric(e.date - s.date) < eff_filter){
      stop("Effort day filter is greater than time interval.")
    } else {
      if(verbose == T){cat("Seasonal window:", 
                           as.numeric(e.date - s.date), 
                           "days. Effort filter:", 
                           eff_filter, "days")}
    }
    
    eff <- eff |> 
      select(all_of(c(eff_site_name, date.cols, eff_summary_cols)))|> 
      rename(Cell_Unit := !!eff_site_name) |> 
      filter(effort_days >= eff_filter)
    
    ## Filter the detection data by effort sites
    d.thresh <- d.thresh |> 
      filter(Cell_Unit %in% eff$Cell_Unit)
  }
  
  if(!is.null(eff) & is.null(eff_summary_cols)){
    
    if(verbose == T){cat("Effort file provided.\nSummary columns missing.
                         \nAttempting to generate summary cols.")}
    
    if(as.numeric(e.date - s.date) < eff_filter){
      stop("Effort day filter is greater than time interval.")
    } else {
      if(verbose == T){cat("Seasonal window:", 
                           as.numeric(e.date - s.date), 
                           "days. Effort filter:", 
                           eff_filter, "days")}
    }
    
    eff <- eff |> 
      select(all_of(c(eff_site_name, date.cols))) |> 
      rename(Cell_Unit := !!eff_site_name)
    
    # Apply function to each row to find the first non-zero value
    first_non_zero_per_row <- apply(eff[, date.cols], 1, function(row) {
      # Find the index of the first non-zero value
      non_zero_cols <- which(row != 0)
      if(length(non_zero_cols) > 0){
        return(names(eff[,date.cols])[non_zero_cols[1]])
      } else {
        return(NA)
      }
    })
    
    # Show the result
    eff$firstday <- purrr::map_vec(first_non_zero_per_row, ~ if (!is.na(.x)) as.Date(.x) else NA)
    eff$firstday_julian <- purrr::map_vec(eff$firstday, ~ if (!is.na(.x)) as.numeric(format(.x, "%j")) else NA)
    eff$effort_days <- rowSums(eff[, date.cols] > 0, na.rm = T)
    
    ## Filter ARUs by effort_days
    eff <- eff |> filter(effort_days >= eff_filter)
    
    ## Filter the ARU detections by this
    d.thresh <- d.thresh |> 
      filter(Cell_Unit %in% eff$Cell_Unit)
    
    ## Conditionally changing 0 to NA for unsurveyed sites
    ## Find the same column names
    eff.sub <- eff[,which(colnames(eff) %in% colnames(d.thresh))]
    d.thresh[eff.sub == 0] <- NA
    
  }
  
  return(list(sp.dat = d.thresh,
              eff.dat = eff,
              species = species,
              thresh_scale = thresh_scale,
              thresh_transform = thresh_transform,
              thresh_trans_dir = ifelse(thresh_transform, thresh_trans_dir, NA),
              thresh_cut = thresh_cut,
              species_thresh_cut = species_thresh_cut,
              date_range = date_range,
              no_dets = no_dets,
              binary = binary,
              effort_filter = eff_filter,
              thresh_cut_meta = thresh_cut_meta))
  
} #function

## -------------------------------------------------------------
##
## End Section: ARU species detection function
##
## -------------------------------------------------------------

## -------------------------------------------------------------
##
## Begin Section: ARU Data Generation function
##
## -------------------------------------------------------------

## Function for loading and writing these files piece meal
## limit the amount of files needed to be read it at once
## Function for loading and writing these files piece meal
## limit the amount of files needed to be read it at once

aru_det_file_gen <- function(det_dir = c("C:/Users/srk252/Documents/Rprojs/Sierra_Biodiv/Data/Detections_By_Species/"),
                             det_years = c("2021", "2022", "2023"),
                             occ_format = T,
                             occ_outdir = c("C:/Users/srk252/Documents/Rprojs/Sierra_Biodiv/Data/Generated_DFs/Occ_Mod_Data/"),
                             seas_format = F,
                             seas_outdir = NULL,
                             eff_file = T,
                             coord_link = F,
                             
                             ## Passing to ac_det_filter()
                             d_thresh = thresh,
                             thresh_scale = "Conf",
                             thresh_transform = FALSE,
                             thresh_trans_dir = "conf2logit",
                             thresh_cut = "99",
                             species_thresh_cut = NA,
                             time_format = "ymd",
                             no_dets = 2,
                             binary = T,
                             date_range = c("2021-06-01", "2021-07-30",
                                            "2022-06-01", "2022-07-30",
                                            "2023-06-01", "2023-07-30"),
                             eff_site_name = "Cell_U",
                             eff_filter = 10,
                             verbose = F,
                             ...){
  
  ## This is written as a for loop for now
  ## Control flow for ditching large files
  ## as they are processed.
  
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##
  ## Subsection: Control Flow for files and directories
  ##
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  if(!is.character(det_dir)){stop("det_dir expects directory path.")}
  tmp.dir <- sapply(det_years, function(x) paste(det_dir, x, sep = "/"))
  tmp.dir <- gsub("//", "/", tmp.dir)
  if(!all(dir.exists(tmp.dir))){
    stop("One or more input directories does not exist.")
  }
  
  ## Check multiple years
  if(length(det_years) > 1){
    message("More than 1 year of data detected.")
  }
  
  ## User prompts for directory management
  if(occ_format == T & !dir.exists(occ_outdir)){
    message("Occupancy output generation selected, but output dir does not exist.")
    
    ## User choice...where to create directory...or kill the fxn
    choice <- menu(c("Create an output directory at WD root.", 
                     "Create an output directory at specific dir.",
                     "Stop function."),
                   title = "Select the following:")
    if(choice == 1){
      message("Dir created here: ")
      print(paste0(getwd(), "/Occ_Mod_Data/", det_years))
      occ_outdir <- sapply(det_years, function(x) paste(getwd(), "Occ_Mod_Data", x, sep = "/"))
      dir.create(paste0(getwd(), "Occ_Mod_Data", det_years))
    } else if (choice == 2){
      occ_outdir <- readline(prompt = "Please enter occ_out dir: ")
      dir.create(occ_outdir)
    } else if(choice == 3){
      stop("Function stopped.")
    } else {
      stop("Invalid selection. Function stopped.")
    }
  }
  
  if(occ_format == T & length(list.files(occ_outdir)) > 0){
    message("Occupancy summary output generation select, but output dir has files.")
    message("Files in dir:", list.files(occ_outdir, full.names = T))
    
    ## User choice...potentially overwrite existing files...
    choice <- menu(c("Do NOT overwrite. Stop function.",
                     "Overwrite files if need be."))
    if(choice == 1){
      stop("Outdir has files and user stopped function.")
    } else if(choice == 2){
      message("User selected to progress with potential file overwriting.")
    } else {
      stop("Invalid selection. Function stopped.")
    }
  }
  
  ## User prompts for directory management
  if(seas_format == T){
    if(!dir.exists(seas_outdir)){
      message("Seasonal summary output generation selected, but output dir does not exist.")
      
      ## User choice...where to create directory...or kill the fxn
      choice <- menu(c("Create an output directory at WD root.", 
                       "Create an output directory at specific dir.",
                       "Stop function."),
                     title = "Select the following:")
      if(choice == 1){
        message("Dir created here: ", paste0(getwd(), "Seasonal_Summary_Data", det_years, collapse = "/"))
        seas_outdir <- paste0(getwd(), "Seasonal_Summary_Data", det_years)
        dir.create(paste0(getwd(), "Seasonal_Summary_Data", det_years))
      } else if (choice == 2){
        occ_outdir <- readline(prompt = "Please enter occ_out dir: ")
        dir.create(occ_outdir)
      } else if(choice == 3){
        stop("Function stopped.")
      } else {
        stop("Invalid selection. Function stopped.")
      }
    }
  }
  
  if(seas_format == T){
    if(length(list.files(seas_outdir)) > 0){
      message("Seasonal summary output generation select, but output dir has files.")
      message("Files in dir:", list.files(seas_outdir, full.names = T))
      
      ## User choice...potentially overwrite existing files...
      choice <- menu(c("Do NOT overwrite. Stop function.",
                       "Overwrite files if need be."))
      if(choice == 1){
        stop("Outdir has files and user stopped function.")
      } else if(choice == 2){
        message("User selected to progress with potential file overwriting.")
      } else {
        stop("Invalid selection. Function stopped.")
      }
    }
  }
  
  if(occ_format == F & seas_format == F){
    stop("No files will be generated. Function stopped.")
  }
  
  
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##
  ## Subsection: For loop for acoustic surveys
  ##
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  ## i-loop for years
  for(i in 1:length(det_years)){
    
    ## Report the year
    print(paste0("Working on year: ", det_years[i]))
    
    ## House keeping
    ## Stache ith year, tmp directory
    yr.tmp <- det_years[i]
    in.dir.tmp <- tmp.dir[i]
    
    ## Set dirs for output
    ## Right not this is only for a root dir (not broken into subdirs for years)
    if(occ_format == T){occ.dir.tmp <- occ_outdir}
    if(seas_format == T){seas.dir.tmp <- seas_outdir} 
    
    ## Pull the date range for the corresponding year
    drange.tmp <- date_range[stringr::str_detect(date_range, pattern = yr.tmp)]
    
    ## Read in multiple species
    sp.det.paths <- list.files(paste0(det_dir, yr.tmp, collapse = "/"), full.names = T)
    sp.det.paths <- sp.det.paths[!str_detect(sp.det.paths, "Effort|flac")]
    
    ## Load the files
    sp.det.files <- lapply(sp.det.paths, function(x) readr::read_csv(x, show_col_types = FALSE))
    names(sp.det.files) <- gsub(".*([0-9]{4}[[:punct:]])", "", gsub("_", " ", gsub(pattern = "_Gt.*", replacement = "", x = sp.det.paths)))
    
    ## Effort file
    if(eff_file == T){
      eff <- readr::read_csv(list.files(paste0(det_dir, yr.tmp, collapse = "/"), full.names = T, pattern = "Effort"))
      eff <- as.data.frame(eff)
    } else {
      eff <- NULL
    }
    
    ## Looping the ac_det_filter fxn over species
    ## Make a list to hold individual species
    sp.det.list <- vector(mode = "list", length = length(sp.det.files))
    thresh_cut_meta_list <- vector(mode = "list", length = length(sp.det.files))
    
    ## Start j loop
    for(j in 1:length(sp.det.files)){
      
      ## Apply the loop to each species in list sp.det.files
      ## args passed from top fxn...need to make this more dynamic 
      ## but ('...') aren't playing nicely
      filter.tmp <- ac_det_filter(d = sp.det.files[[j]],
                                  d_thresh = d_thresh,
                                  thresh_scale = thresh_scale,
                                  thresh_transform = thresh_transform,
                                  thresh_trans_dir = thresh_trans_dir,
                                  thresh_cut = thresh_cut,
                                  species_thresh_cut = species_thresh_cut,
                                  time_format = time_format,
                                  species = names(sp.det.files)[j],
                                  no_dets = no_dets,
                                  binary = binary,
                                  date_range = drange.tmp,
                                  eff = eff,
                                  eff_site_name = eff_site_name,
                                  eff_filter = eff_filter,
                                  verbose = verbose
      )
      
      ## Take the output from the first list element returned
      ## This is the species DF
      sp.det.list[[j]] <- filter.tmp$sp.dat
      
      ## filter temp staches the metadata for the function
      names(sp.det.list)[j] <- filter.tmp$species
      
      ## Report species
      message("Processed and adding:", names(sp.det.list)[j])
      
      ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ##
      ## Subsection: Save metadata for species thresholds if needed
      ##
      ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if(class(filter.tmp$thresh_cut_meta) == "data.frame"){
        thresh_cut_meta_list[[j]] <- filter.tmp$thresh_cut_meta
        if(j == length(sp.det.files)){
          thresh_cut_meta_df <- do.call(rbind, thresh_cut_meta_list)
          write.csv(thresh_cut_meta_df, file = paste0(occ.dir.tmp, yr.tmp, "_", thresh_cut, thresh_scale, "_OccSpeciesThreshMetaData.csv"))
        }
      }
      
      ## Create a metadata DF for the user so they can keep track 
      ## of fxn settings
      if(j == length(sp.det.files)){
        meta.df <- data.frame(Thresh_scale = filter.tmp$thresh_scale,
                              Thresh_cut = filter.tmp$thresh_cut,
                              Species_Thresh = filter.tmp$species_thresh_cut,
                              Thresh_transform = filter.tmp$thresh_transform,
                              Thresh_trans_dir = filter.tmp$thresh_trans_dir,
                              Start_date = filter.tmp$date_range[1],
                              End_date = filter.tmp$date_range[2],
                              Survey_Year = lubridate::year(filter.tmp$date_range[1]),
                              NumberDetections = filter.tmp$no_dets,
                              Binarized = filter.tmp$binary,
                              Effort_Filter = filter.tmp$effort_filter)
      }
      
      ## Create metadata file for the subsetted effort file
      if(j == length(sp.det.files) & eff_file){
        write.csv(filter.tmp$eff.dat, file = paste0(occ.dir.tmp, yr.tmp, "_", thresh_cut, thresh_scale, "_OccEffortFileSubset.csv"))
      }
      
    } #j closure
    
    names(sp.det.list)
    
    ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ##
    ## Subsection: File saving for occupancy model format
    ##
    ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    ## Create the data arrays for occupancy models
    ## Write the list as .RDA files + metadata
    if(occ_format == T){
      cat("Saving list of species filtered detections DFs to .RData file for year:", yr.tmp)
      save(sp.det.list, file = paste0(occ.dir.tmp,
                                      yr.tmp,
                                      "_", 
                                      thresh_cut, 
                                      thresh_scale,
                                      "_OccSppList.RData"))
      write.csv(meta.df, file = paste0(occ.dir.tmp, yr.tmp, "_", thresh_cut, thresh_scale, "_OccMetaFilterParams.csv"))
    }
    
    ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ##
    ## Subsection: Creating the seasonal community matrices
    ##
    ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    if(seas_format == T){
      
      ## Sum the detections in species richness DF
      cols_to_sum <- 2:ncol(sp.det.list[[1]])
      spec.rich <- Reduce(function(x,y){
        x[,cols_to_sum] <- x[,cols_to_sum] + y[,cols_to_sum]
        x
      }, sp.det.list)
      
      ## If we want to retrieve species richness for each site across the entire season
      season_pool <- function(x){
        season_occ <- x |> 
          mutate(row_sum = rowSums(across(where(is.numeric)))) |> 
          select(where(~ !is.numeric(.)), row_sum) |> 
          mutate(row_sum = ifelse(row_sum > 0, 1, row_sum)) |> 
          rename(Season_Occ = row_sum)
        
        return(season_occ)
      }
      
      ## Apply this function to the list of species
      seasonal_occ <- lapply(sp.det.list, season_pool)
      
      ## The join column should always be this...
      ## but perhaps this can change
      join_column <- "Cell_Unit"
      
      # Join the data frames in the list and rename columns
      seasonal_df <- reduce(names(seasonal_occ), function(acc, name) {
        df <- seasonal_occ[[name]]
        
        # Rename columns by appending the data frame name
        colnames(df)[-which(colnames(df) == join_column)] <- paste(name, colnames(df)[-which(colnames(df) == join_column)], sep = "_")
        colnames(df) <- gsub("_Season_Occ", "", colnames(df))
        
        # Perform the left join
        if (is.null(acc)) {
          return(df)
        } else {
          return(left_join(acc, df, by = join_column))
        }
      }, .init = NULL)
      
      ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ##
      ## Subsection: Link with the GPS coordinates from deployments
      ##
      ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if(coord_link == T){
        ## Connect to DB
        cb_connect_db()
        
        ## SF needs to be loaded to execute
        ownership <- c('any')
        cell_list <- cb_cells_by_ownership(ownership)
        
        ## Study type
        study <- c('Sierra Monitoring')
        
        ## Get the deployment information
        deployments_df <- 
          conn |> 
          dplyr::tbl('acoustic_field_visits') |> 
          dplyr::filter(
            is_invalid == 0,
            cell_id %in% cell_list,
            study_type %in% study,
            survey_year == yr.tmp
          ) |> 
          collect() |> 
          dplyr::select(deployment_name, survey_year, matches("utm"))
        
        ## Disconnect DB
        cb_disconnect_db()
        
        ## Create SF object
        deployments_geo <- 
          deployments_df |> 
          group_split(utm_zone) |> 
          map_dfr(cb_make_aru_sf) |> 
          mutate(Cell_Unit = stringr::str_remove(deployment_name, "G[0-9]+_V[0-9]+_")) |> 
          st_transform(crs = 4326) |> 
          dplyr::mutate(lon = sf::st_coordinates(geometry)[,2],
                        lat = sf::st_coordinates(geometry)[,1]) |> 
          st_drop_geometry() |>
          ## Remove
          select(Cell_Unit, deployment_name, survey_year, lon, lat)
        
        ## Link with the seasonal_df
        seasonal_df_geo <- seasonal_df |>
          mutate(Cell_Unit = ifelse(
            stringr::str_detect(string = Cell_Unit, pattern = "C[0-9]{3}"),
            gsub(pattern = "(C)([0-9]{3})(_U[0-9]+)$", replacement = "\\10\\2\\3", x = Cell_Unit),
            Cell_Unit
          )) |> 
          left_join(deployments_geo) |> 
          select(Cell_Unit, deployment_name, survey_year, lon, lat, everything())
        
        ## Rename the new df to save
        seasonal_df <- seasonal_df_geo
        
      }
      
      ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ##
      ## Subsection: Write files for the seasonal dfs
      ##
      ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ## Write the DF
      data.table::fwrite(seasonal_df,
                         paste0(seas_outdir, "/", yr.tmp, "_SeasonalSpeciesMat.csv"))
      
      ## Write the metadata
      write.csv(meta.df, file = paste0(seas_outdir, "/", yr.tmp, "_SeasonalMetaFilterParams.csv"))
      
    } #seasonal format
  } #for loop i
} #function end

## -------------------------------------------------------------
##
## End Section: ARU Data Generation function
##
## -------------------------------------------------------------

## -------------------------------------------------------------
##
## Begin Section: Query Spatial Locations from CAbioacoustics
##
## -------------------------------------------------------------

aru_sf_query <- function(years = c(2021, 2022, 2023, 2024),
                         .study_type = "Sierra Monitoring",
                         .is_invalid = 0){
# connect to database -----------------------------------------------------
require(CAbioacoustics)

# create connection; credentials stored in keyring
cb_connect_db()

# see that connection worked and list database tables
DBI::dbListTables(conn)

# just sierra monitoring deployments
deployments_df <-
  tbl(conn, "acoustic_field_visits") |>
  filter(survey_year %in% years & study_type == .study_type & is_invalid == .is_invalid) |>
  select(
    id,
    study_type,
    group_id,
    visit_id,
    cell_id,
    unit_number,
    deployment_name,
    survey_year,
    deploy_date,
    recover_date,
    utm_zone,
    utme,
    utmn
  ) |>
  # pull into memory
  collect()

# disconnect now
cb_disconnect_db()


# convert ARU data to simple feature point object -------------------------

# make sf, convert to WGS 84
deployments_sf <-
  deployments_df |>
  group_split(utm_zone) |>
  map_dfr(cb_make_aru_sf) |> 
  mutate(Long = st_coordinates(geometry)[,1],
         Lat = st_coordinates(geometry)[,2]) |>
  mutate(Cell_Unit = paste0(cell_id, "_", unit_number))

## Return the df
return(deployments_sf)

}

