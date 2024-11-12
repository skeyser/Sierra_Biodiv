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
## Begin Section: ARU and Species Detection Data
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
## *************************************************************

## -------------------------------------------------------------
##
## Begin Section: ARU species detection function
##
## -------------------------------------------------------------

## *************************************************************
##
## Section Notes:
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
                          thresh_cut = "99",
                          time_format = "ymd",
                          species,
                          no_dets = 1,
                          binary = T,
                          date_range = NULL,
                          eff,
                          eff_site_name,
                          eff_summary_cols,
                          eff_filter = 1) {
  
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
    
    ## Bring in the effort variables
    if(!missing(eff)){
      cat("Missing effort covariates.\nData is not filtered for ARU effort.")
    } else {
      cat("Effort file provided.\nChecking effort covariates.")
    }
    
    if(!missing(eff) & !missing(eff_summary_cols)){
      cat("Effort file provided.\nEffort summary columns provided.")
      if(as.numeric(e.date - s.date) < eff_filter){
        stop("Effort day filter is greater than time interval.")
      } else {
          cat("Seasonal window:", as.numeric(e.date - s.date), "days. Effort filter:", eff_filter, "days")
      }
      
      cat()
      
      eff <- eff |> 
        select(all_of(c(eff_site_name, date.cols, eff_summary_cols)))|> 
        rename(Cell_Unit := !!eff_site_name) |> 
        filter(effort_days >= eff_filter)
      
      ## Filter the detection data by effort sites
      d.thresh <- d.thresh |> 
        filter(Cell_Unit %in% eff$Cell_Unit)
    }
    
    if(!missing(eff) & missing(eff_summary_cols)){
      
      cat("Effort file provided.\nSummary columns missing.\nAttempting to generate summary cols.")
      
      if(as.numeric(e.date - s.date) < eff_filter){
        stop("Effort day filter is greater than time interval.")
      } else {
        cat("Seasonal window:", as.numeric(e.date - s.date), "days. Effort filter:", eff_filter, "days")
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
      eff$firstday_julian <- as.numeric(format(eff$firstday, "%j"))
      eff$effort_days <- rowSums(eff[, date.cols] > 0, na.rm = T)
      
      ## Filter ARUs by effort_days
      eff <- eff |> filter(effort_days >= eff_filter)
      
      ## Filter the ARU detections by this
      d.thresh <- d.thresh |> 
        filter(Cell_Unit %in% eff$Cell_Unit)
      
    }
    
    ## Remove superfluous DFs made
    rm(tmp.detect)
  }
  return(d.thresh)
  
} #function

fun.check <- ac_det_filter(d = test,
              d_thresh = thresh,
              thresh_scale = "Conf",
              thresh_cut = "99",
              time_format = "ymd",
              species = names(sp.det.files)[i],
              no_dets = 2,
              binary = T,
              date_range = c("2021-06-01", "2021-07-30"),
              eff = effort,
              eff_site_name = "Cell_U",
              eff_filter = 10)


## -------------------------------------------------------------
##
## End Section: ARU species detection function
##
## -------------------------------------------------------------

## -------------------------------------------------------------
##
## Begin Section: Processing species detections
##
## -------------------------------------------------------------

## Read in multiple species
sp.det.paths <- list.files(here("./Data/Detections_By_Species/2021/"), full.names = T)
sp.det.paths <- sp.det.paths[!str_detect(sp.det.paths, "Effort")]

sp.det.files <- lapply(sp.det.paths, readr::read_csv)
names(sp.det.files) <- gsub(".*([0-9]{4}[[:punct:]])", "", gsub("_", " ", gsub(pattern = "_Gt.*", replacement = "", x = sp.det.paths)))

## TEMPORARY
sp.det.paths22 <- list.files(here("./Data/Detections_By_Species/2022/"), full.names = T)
sp.det.paths22 <- sp.det.paths22[!str_detect(sp.det.paths22, "Effort")]

sp.det.files22 <- lapply(sp.det.paths22, readr::read_csv)
names(sp.det.files22) <- gsub(".*([0-9]{4}[[:punct:]])", "", gsub("_", " ", gsub(pattern = "_Gt.*", replacement = "", x = sp.det.paths22)))

## Effort file
effort <- readr::read_csv(here("./Data/Detections_By_Species/2021/2021_Effort_perUnit_0400-0859_1800-1959_split0.csv"))


## Looping the function over species
sp.det.list <- vector(mode = "list", length = length(sp.det.files))
for(i in 1:length(sp.det.files)){
  sp.det.list[[i]] <- ac_det_filter(d = sp.det.files[[i]],
                              d_thresh = thresh,
                              thresh_scale = "Conf",
                              thresh_cut = "99",
                              time_format = "ymd",
                              species = names(sp.det.files)[i],
                              no_dets = 2,
                              binary = T,
                              date_range = c("2021-06-01", "2021-07-30"))
  names(sp.det.list)[i] <- names(sp.det.files)[i]
  print(paste("Processed and adding:", names(sp.det.list)[i]))
  
}

names(sp.det.list)

## # of detections
lapply(sp.det.list, function(x) paste(names(x), "detections:", sum(x[,2:ncol(x)])))
paste(names(sp.det.list[1]), "detections:", sum(sp.det.list[[1]][2:ncol(sp.det.list[[1]])]))
paste(names(sp.det.list[2]), "detections:", sum(sp.det.list[[2]][2:ncol(sp.det.list[[2]])]))
paste(names(sp.det.list[3]), "detections:", sum(sp.det.list[[3]][2:ncol(sp.det.list[[3]])]))

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

seasonal_occ <- lapply(sp.det.list, season_pool)

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

# Print the result
print(seasonal_df)

## `seasonal_df` now captures the naive occupancy for each species
## across the Sierra ARUs and can be used as a site x species matrix
colnames(seasonal_df)

## -------------------------------------------------------------
##
## Begin Section: ARU spatial data
##
## -------------------------------------------------------------

## *************************************************************
##
## Section Notes: We now have the ARU data formatted into an ARU x Date matrix 
## with cells corresponding to the BirdNET confidence scores 
## (thresholded to the maximum detection score per day >= the estimated BirdNET thresholds). 
## The cells now represent naive occupancy for the species of interest. 
## Next, we will want to assign the ARUs a spatial location using spatial references for each ARU.
# 
# Spatial data can be accessed via the `CAbioacoustics` package under the table "acoustic_field_visits".
##
## *************************************************************

## Connect to the database
## Need to be connected to VPN for this to connect
cb_connect_db()

## Set some values for filters
syear <- 2021
eyear <- 2023

## SF needs to be loaded to execute
ownership <- c('any')
cell_list <- cb_cells_by_ownership(ownership)

## Study type
study <- c('Sierra Monitoring')

## Query and pull the table we want
deployments_df <- 
  conn |> 
  dplyr::tbl('acoustic_field_visits') |> 
  dplyr::filter(
    is_invalid == 0,
    cell_id %in% cell_list,
    study_type %in% study,
    survey_year >= syear,
    survey_year <= eyear
  ) |> 
  collect() |> 
  dplyr::select(deployment_name, survey_year, matches("utm"))

## Disconnect from the DB
cb_disconnect_db()

## SF creation
deployments_sf <- 
  deployments_df |> 
  group_split(utm_zone) |> 
  map_dfr(cb_make_aru_sf)

## Plot with MV
## hexes
hexes_map <-
  cb_get_spatial('sierra_hexes') |>
  mapview(layer.name = 'Hexes')

## deployments
deployments_map <-
  deployments_sf |>
  mutate(survey_year = as.factor(survey_year)) |>
  mapview(zcol = 'survey_year', layer.name = 'Survey year')

## combine
hexes_map + deployments_map


## Link Species Detections with ARU locations - From Jay's DB
glimpse(seasonal_df)
head(seasonal_df$Cell_Unit)
length(unique(seasonal_df$Cell_Unit))


glimpse(deployments_sf)

## Take only the 2021 spatial data
dep21 <- deployments_sf |> 
  filter(survey_year == 2021) |> 
  mutate(Cell_Unit = stringr::str_remove(deployment_name, "G[0-9]+_V[0-9]+_")) |> 
  st_transform(crs = 4326) |> 
  dplyr::mutate(lon = sf::st_coordinates(geometry)[,2],
                lat = sf::st_coordinates(geometry)[,1]) |> 
  st_drop_geometry() |>
  ## Remove
  select(Cell_Unit, survey_year, lon, lat) |> 
  distinct()

## *************************************************************
##
## Section Notes:
## Check the number of unit cell_unit combinations
## Check duplications in the cell_units
## Duplicates arise via different deployment_ids but not cell_unit ids
## V1 vs V2 in the deployment_id
## Step3 data only contains the cell_unit ID
## Deployments from owl DB have 1730 unique deployment_ids
## ATFL example species from 2021 has 1652 individual ARU deployments
##
## *************************************************************

nrow(dep21)
length(unique(dep21$deployment_name))
length(unique(dep21$Cell_Unit))
dup_unit <- dep21[which(duplicated(dep21$Cell_Unit)),]$Cell_Unit
dup21 <- dep21[dep21$Cell_Unit %in% dup_unit, ]
dep21 <- distinct(dep21)

## Match the ATFL data with the 2021 survey data
seasonal_df_meta <- seasonal_df |>
  mutate(Cell_Unit = ifelse(
    stringr::str_detect(string = Cell_Unit, pattern = "C[0-9]{3}"),
    gsub(pattern = "(C)([0-9]{3})(_U[0-9]+)$", replacement = "\\10\\2\\3", x = Cell_Unit),
    Cell_Unit
  )) |> 
  left_join(dep21) |> 
  select(Cell_Unit, survey_year, lon, lat, everything()) 
  
missing <- dep21[!dep21$Cell_Unit %in% seasonal_df_meta$Cell_Unit,]
# 
# atfl.sf.ok <- atfl.sf[which(!is.na(atfl.jay$survey_year)),]
# atfl.sf.ok <- atfl.sf.ok |> 
#   dplyr::select(Cell_Unit, deployment_name, survey_year, lon, lat, everything()) |> 
#   rowwise() |> 
#   mutate(sum = sum(c_across("2021-03-01":"2021-08-30"))) |> 
#   ungroup() |> 
#   dplyr::select(Cell_Unit, deployment_name, survey_year, lon, lat, sum) |> 
#   st_as_sf(coords = c("lon", "lat"), crs = 4326)
# 
# ggplot(data = atfl.sf.ok) + 
#   geom_sf(aes(color = sum))



## -------------------------------------------------------------
##
## Begin Section: 2021 ARU Metadata
##
## -------------------------------------------------------------
## *************************************************************
##
## Section Notes:
## Data file with locations from Connor. 
## It seems like there are some inconsistencies with the leading zero 
## in the cell ID name. 
## My guess is that importing the files into excel dropped the leading zeros 
## since all of the issues are occurring with cells that have 3 numbers when they should have 4. 
## Add in a leading zero if the string only has 3 numbers after the "C" and let's see what happens.
##
## The 2021 ARU meta data file from Connor is specific to 2021, where the DB from Jay has the location
## information for each year. 
## 
## For analyses, I plan to extract/create variables on interest so I really just need to spatial
## locations and the extractions can be done on my part.
##
## *************************************************************

## Load in the ARU metadata
aru_meta <- readr::read_csv(here("Data/ARU_120m.csv"))
names(aru_meta)
glimpse(aru_meta)
aru_meta$Cell_Unit <- paste0(aru_meta$cell_id, "_", aru_meta$unit_numbe)

## Take a couple of things for testing
aru_meta <- aru_meta |>
  dplyr::select(Cell_Unit, survey_yea, X, Y, topo_elev, standage_f3_mn, contains("fire"))

## Fixing the leading zero issue
seasonal_df <- seasonal_df |>
  mutate(Cell_Unit = ifelse(
    stringr::str_detect(string = Cell_Unit, pattern = "C[0-9]{3}"),
    gsub(pattern = "(C)([0-9]{3})(_U[0-9]+)$", replacement = "\\10\\2\\3", x = Cell_Unit),
    Cell_Unit
  ))

## Attempt to rejoin now
seasonal_df_meta <- seasonal_df |> 
  left_join(aru_meta) |> 
  select(Cell_Unit, 
         Year = survey_yea, 
         Long = X, 
         Lat = Y,
         everything())

## ARU Meta
aru_filter <- aru_meta |> 
  filter(Cell_Unit %in% seasonal_df_meta$Cell_Unit)

## Write the DF
data.table::fwrite(seasonal_df_meta,
                   here::here("Data/Generated_DFs/2021_SeasonalSpeciesMat.csv"))

