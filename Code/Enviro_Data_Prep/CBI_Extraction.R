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
library(stringr)
library(ggplot2)
library(here)
library(terra)
library(sf)
library(purrr)
library(CAbioacoustics)
## -------------------------------------------------------------

## Pull in the ARU meta data for point locations

## Take a look at CBI data
cbi <- rast("c:/Users/srk252/Documents/data_for_spencer/cbi_sierra_cat_rasters/cbi_cat_2021.tif")
plot(cbi)

## ARU metadata
## Take from Jays package?
## For now just use the one from Connor
aru_meta <- readr::read_csv(here("Data/ARU_120m.csv"))
names(aru_meta)
glimpse(aru_meta)
aru_meta$Cell_Unit <- paste0(aru_meta$cell_id, "_", aru_meta$unit_numbe)
aru_locs <- aru_meta |> 
  select(Cell_Unit, survey_yea, Long = X, Lat = Y)

aru_locs <- st_as_sf(aru_locs, coords = c("Long", "Lat"), crs = 4326)

## -------------------------------------------------------------
##
## Begin Section: Extraction function body
##
## -------------------------------------------------------------
## Environmental extraction function
## Acceptable fire products
## fire_interval, most_recent_fire, fire_frequency, number of fires,
## and landscape metrics
aru_fire_prep <- function(fire_prod = NULL, # character vector of desired fire output
                          locs_from_cabio = T, #override flag for using CAbioacoustic locations and metadata
                          aru_locs = NULL, # Data.frame with coordinates for custom locations                          
                          survey_years = NULL, # Survey year is only applicable for locs_from_cabio = TRUE
                          id_col = "deployment_name",
                          year_col = NULL,
                          x_col = "Long", # chr for x coordinate col name
                          y_col = "Lat", # chr for y coordinate col name
                          .crs = 4326, # default WGS84 for coordinates
                          des_out, # desired output
                          spat_ex, # chr for type (point, buff, hex)
                          buff_size = NULL, # vector of buffer sizes
                          #hex_size = NULL,
                          intervals = c("1-5", "6-10", "11-35"),
                          landscape_metrics = T,
                          ...
){
  
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##
  ## Subsection: CAbioacoustics Spatial
  ##
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  if(locs_from_cabio){aru_locs <- cabio_loc_query(years = survey_years)}
  
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##
  ## Subsection: Create SF objects from the ARU coordinates
  ##
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  ## Create SF objects from points
  if(!any(class(aru_locs) %in% "sf")){
    
    aru_locs <- st_as_sf(aru_locs, 
                         coords = c(x_col, 
                                    y_col), 
                         crs = .crs)
  }
  
  
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##
  ## Subsection: CBI Categorical Extraction and prep
  ##
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  ## Find the data for extraction
  #if(env_prod == "Fire_CBI"){
  
  ## Function bundle
  
  if(file.exists("C:/Users/srk252/Documents/GIS_Data/CBI_Sierra/CBI_1985_2024_ZeroFilling_Stack.tif")){
    
    print("CBI rasters are fixed and exist in directory. Loading the fixed stack.")
    cbi_stack <- rast("C:/Users/srk252/Documents/GIS_Data/CBI_Sierra/CBI_1985_2024_ZeroFilling_Stack.tif")
  
    } else {
    
    ## Find files
    cbi_path <- "C:/Users/srk252/Documents/data_for_spencer/cbi_sierra_cat_rasters/"
    cbi_files <- list.files(cbi_path, full.names = T, pattern = "(cbi_cat_)(\\d{4})(*.tif$)")
    
    ## Call function to stack if need be
    cbi_stack <- cbi_stack_fun(cbi_files)
    
    ## Add in zeros for the raster for downstream processing
    temp <- rast(nrows = nrow(cbi_stack[[1]]),
                 ncols = ncol(cbi_stack[[1]]),
                 xmin = xmin(cbi_stack[[1]]),
                 xmax = xmax(cbi_stack[[1]]),
                 ymin = ymin(cbi_stack[[1]]),
                 ymax = ymax(cbi_stack[[1]]),
                 crs = crs(cbi_stack[[1]]),
                 resolution = res(cbi_stack[[1]]),
                 vals = 0,
                 names = "template"
    )
    
    ## Reclassify the NAs to zero to fill in the rasters
    cbi_stack <- classify(cbi_stack, cbind(NA, 0))
    if(nlyr(cbi_stack) == length(cbi_files)){
      names(cbi_stack) <- str_extract(cbi_files, "\\d{4}")
    }
  }
  
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##
  ## Subsection: Fire year interval stacking
  ##
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  cbi_int <- int_ras_fun(ras_stack = cbi_stack,
                         intervals = intervals,
                         sum_int = "max")
  
  ## Project the points
  aru_locs <- st_transform(aru_locs,
                           crs = crs(cbi_stack))
  
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##
  ## Subsection: Fire Severity
  ##
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if("fire_severity" %in% fire_prod){
    if(!is.null(intervals)){
      message("Intervals set. Fire severity calculated from summarized fire years.")
      cbi_sev <- cbi_int
    } else {cbi_sev <- cbi_stack}
    if(!is.null(buff_size)){
      aru_buffer <- st_buffer(aru_locs,
                              dist = buff_size)
      fire_sev <- exactextractr::exact_extract(cbi_sev,
                                               aru_buffer,
                                               fun = c("mean", "stdev"),
                                               append_cols = id_col)
    } else {
      fire_sev <- terra::extract(cbi_sev,
                                 aru_locs,
                                 bind = T) |> st_as_sf()
    }
    colnames(fire_sev)[colnames(fire_sev) != id_col] <- paste0("Fire_Sev_", gsub("[[:punct:]]", "_", colnames(fire_sev)[colnames(fire_sev) != id_col]))
    
  }
  
  
  
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##
  ## Subsection: Calculate fire variables
  ##
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  years <- as.numeric(names(cbi_stack))
  
  if("time_since_fire" %in% fire_prod){
    if(!is.null(intervals)){
      message("Intervals are set. Time since fire doesn't accept intervals...output will be for single years.")
    } else {
      ## Time to most recent fire
      time_since_fire <- terra::app(cbi_stack, 
                                    fun = function(x) time_to_most_recent_fire(cell_values = x,
                                                                               years = years))
    }}
  
  if("fire_freq" %in% fire_prod){
    if(!is.null(intervals)){
      message("Intervals are set. Fire frequency doesn't accept intervals.")
    } else {
      ## Fire frequency (num fires/total record length)
      fire_freq <- terra::app(cbi_stack, 
                              fun = function(x) fire_freq_calc(cell_values = x,
                                                               years = years))
    }}
  
  if("fire_ret_int" %in% fire_prod){
    if(!is.null(intervals)){
      message("Intervals are set. Fire return interval doesn't accept intervals.")
    } else {
      ## Fire return interval (mean of time between successive fires)
      fire_return_int <- terra::app(cbi_stack, 
                                    fun = function(x) fire_return_int(cell_values = x, 
                                                                      years = years))
    }}
  
  ## Report the status of s2 geometry for convenience
  if(!is.null(buffer_size)){
    if(sf_use_s2() == T){
      print("s2 geometry enabled. Buff_size interpreted as meters.")
    } else {
      print("s2 disabled, if locations are geodetic (lat/lon) units interpretted as degrees.")
    }
  }
  
  
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##
  ## Subsection: Buffer Extraction from Custom Fire Variables
  ##
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  ## Buffer extraction code from Jay (tidy code)
  ## Retrieve variable names from the
  ## WIP
  if(any(fire_prod %in% c("time_since_fire", "fire_freq", "fire_ret_int"))){
    variable_name <- fire_prod
    fire_buff_out <- vector(mode = "list", length = length(fire_prod))
    
    for(var in 1:length(variable_name)){
      
      var.tmp <- variable_name[var]
      
      if(var.tmp == "time_since_fire"){
        fire_buff_extract <-
          buff_size %>%
          # create column names
          str_c(var.tmp, ., sep = '_') |>
          set_names() |>
          map_dfr(
            \(x)
            exactextractr::exact_extract(
              time_since_fire,
              # buffer points
              aru_locs |> st_buffer(as.numeric(str_extract(x, "\\d"))),
              fun = 'mean'
            )
          ) |> 
          bind_cols(st_drop_geometry(aru_locs[,id_col]))
        
        fire_buff_out[[var]] <- fire_buff_extract
        
      }
      
      if(var.tmp == "fire_freq"){
        fire_buff_extract <-
          buff_size %>%
          # create column names
          str_c(var.tmp, ., sep = '_') |>
          set_names() |>
          map_dfr(
            \(x)
            exactextractr::exact_extract(
              fire_freq,
              # buffer points
              aru_locs |> st_buffer(as.numeric(str_extract(x, "\\d"))),
              fun = 'mean'
            )
          ) |> 
          bind_cols(st_drop_geometry(aru_locs[,id_col]))
        
        fire_buff_out[[var]] <- fire_buff_extract
      }
      
      if(var.tmp == "fire_ret_int"){
        fire_buff_extract <-
          buff_size %>%
          # create column names
          str_c(var.tmp, ., sep = '_') |>
          set_names() |>
          map_dfr(
            \(x)
            exactextractr::exact_extract(
              fire_return_int,
              # buffer points
              aru_locs |> st_buffer(as.numeric(str_extract(x, "\\d"))),
              fun = 'mean'
            )
          ) |> 
          bind_cols(st_drop_geometry(aru_locs[,id_col]))
        
        fire_buff_out[[var]] <- fire_buff_extract
      }
      
    }
    fire_buff_merge <- Reduce(function(x, y) merge(x, y, by = id_col), fire_buff_out) 
  }
  
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##
  ## Subsection: Landscape Metrics
  ##
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## Function should estimate specific landscape metrics for the
  ## interval-based raster data for CBI
  fire_lscp_out <- fire_lscp_fun(ras_int = cbi_int,
                ras_stack = cbi_stack,
                locs = aru_locs,
                buff_size = buff_size,
                id_col = id_col,
                metrics = c("lsm_c_pland"))
  
#   fire_lscp_fun <- function(cbi_stack, aru_locs){
#     for(i in 1:length(cbi_int)){
#       
#       print(paste("Year:", i))
#       
#       ## Pull the ARUs for a specific year
#       aru_tmp <- aru_locs |> 
#         filter(survey_year == aru_years[i])
#       
#       ## Check that we intend to calculate landscape metrics
#       if(landscape_metrics == TRUE){
#         
#         ## Is there a buffer size?
#         if(!is.null(buff_size)){ #| is.null(hex_size)){
#           
#           ## Check to make sure sf
#           if(!is.null(id_col) & any(str_detect(class(aru_tmp), "sf"))){
#             id.v <- as.vector(unlist(st_drop_geometry(aru_tmp[, id_col])))
#           } else {
#             id.v <- NULL
#           }
#           
#           ## Individual years versus intervals for the fire data
#           if(is.null(intervals)){ 
#             message("Landscapemetrics done for individual years.")
#             cbi_lcpc <- cbi_stack 
#           } else { 
#             cbi_lcpc <- rast(cbi_int[i])
#             message("Landscapemetrics done for intervals.")
#           }
#           
#           ## Calculate LSCP
#           cbi_lsm <- landscapemetrics::sample_lsm(landscape = cbi_lcpc,
#                                                   y = aru_tmp,
#                                                   plot_id = id.v,
#                                                   shape = "circle",
#                                                   size = buff_size,
#                                                   what = c("lsm_c_pland"
#                                                            #"lsm_c_ed",
#                                                            #"lsm_c_lpi",
#                                                            #"lsm_l_shdi",
#                                                            #"lsm_l_msidi",
#                                                            #"lsm_l_msiei"
#                                                   ),
#                                                   #metric = "core",
#                                                   #level = "patch",
#                                                   return_raster = F,
#                                                   #type = "diversity metric",
#                                                   #classes_max = 3,
#                                                   verbose = FALSE,
#                                                   progress = T)
#           
#           ## Fix layer names
#           lyr_map <- data.frame(layer = 1:nlyr(cbi_lcpc), lyr_name = intervals)
#           
#           ## Merge these
#           cbi_lsm <- left_join(cbi_lsm, lyr_map)
#           
#           ## level of lsm
#           lsm_lev <- unique(cbi_lsm$level)
#           
#           if(any(lsm_lev == "landscape")){
#             cbi_lsm_l <- cbi_lsm |> 
#               filter(level == "landscape") |> 
#               select(plot_id, level, lyr_name, metric, value) |> 
#               tidyr::pivot_wider(names_from = c(metric, lyr_name), 
#                                  values_from = value,
#                                  names_glue = "{lyr_name}_{metric}_l") |> 
#               select(!level) |>
#               mutate(Year = aru_years[i])
#             
#           }
#           
#           if(any(lsm_lev == "class")){
#             
#             fire_class <- c("Unburned", "Low_Sev", "Mod_Sev", "High_Sev")
#             
#             cbi_lsm_c <- cbi_lsm |> 
#               filter(level == "class") |> 
#               mutate(FireClass = case_when(class == 0 ~ "Unburned",
#                                            class == 1 ~ "Low_Sev",
#                                            class == 2 ~ "Mod_Sev",
#                                            class == 3 ~ "High_Sev")) |> 
#               select(plot_id, level, lyr_name, FireClass, metric, value) |> 
#               tidyr::pivot_wider(names_from = c(metric, FireClass, lyr_name), 
#                                  values_from = value,
#                                  names_glue = "{FireClass}_{lyr_name}_{metric}_c") |> 
#               select(!level) |> 
#               mutate(Year = aru_years[i])
#             
#           }
#           
#           if(exists("cbi_lsm_l") & exists("cbi_lsm_c")){ cbi_lsm <- full_join(cbi_lsm_l, cbi_lsm_c) } 
#           if (exists("cbi_lsm_l") & !exists("cbi_lsm_c")) { cbi_lsm <- cbi_lsm_l }  
#           if (!exists("cbi_lsm_l") & exists("cbi_lsm_c")) { cbi_lsm <- cbi_lsm_c }
#         }
#       }
#       #} #lscpmet
#       
#       ## Grow the DF for multiple years
#       if(i == 1){
#         cbi_lsm_out <- cbi_lsm
#       } else {
#         cbi_lsm_out <- rbind(cbi_lsm_out, cbi_lsm)
#       }
#       
#     }
#     
#     ## Block for what to do with CBI Fire Data
#     if(landscape_metrics == TRUE){
#       if(!is.null(buff_size)){ #| is.null(hex_size)){
#         if(!is.null(id_col) & any(str_detect(class(aru_locs), "sf"))){
#           id.v <- as.vector(unlist(st_drop_geometry(aru_locs[, id_col])))
#         } else {
#           id.v <- NULL
#         }
#         if(is.null(intervals)){ 
#           message("Landscapemetrics done for individual years.")
#           cbi_lcpc <- cbi_stack 
#         } else { 
#           cbi_lcpc <- cbi_int
#           message("Landscapemetrics done for intervals.")
#         }
#         cbi_lsm <- landscapemetrics::sample_lsm(landscape = cbi_lcpc,
#                                                 y = aru_locs,
#                                                 plot_id = id.v,
#                                                 shape = "circle",
#                                                 size = buff_size,
#                                                 what = c("lsm_c_pland",
#                                                          "lsm_c_ed",
#                                                          "lsm_c_lpi",
#                                                          "lsm_l_shdi",
#                                                          "lsm_l_msidi",
#                                                          "lsm_l_msiei"),
#                                                 #metric = "core",
#                                                 #level = "patch",
#                                                 return_raster = F,
#                                                 #type = "diversity metric",
#                                                 #classes_max = 3,
#                                                 verbose = FALSE,
#                                                 progress = T)
#         
#         ## Fix layer names
#         lyr_map <- data.frame(layer = 1:nlyr(cbi_lcpc), lyr_name = names(cbi_lcpc))
#         
#         ## Merge these
#         cbi_lsm <- left_join(cbi_lsm, lyr_map)
#         
#         ## level of lsm
#         lsm_lev <- unique(cbi_lsm$level)
#         
#         if(any(lsm_lev == "landscape")){
#           cbi_lsm_l <- cbi_lsm |> 
#             filter(level == "landscape") |> 
#             select(plot_id, level, lyr_name, metric, value) |> 
#             tidyr::pivot_wider(names_from = c(metric, lyr_name), 
#                                values_from = value,
#                                names_glue = "{lyr_name}_{metric}_l") |> 
#             select(!level)
#           
#         }
#         
#         if(any(lsm_lev == "class")){
#           
#           fire_class <- c("Unburned", "Low_Sev", "Mod_Sev", "High_Sev")
#           
#           cbi_lsm_c <- cbi_lsm |> 
#             filter(level == "class") |> 
#             mutate(FireClass = case_when(class == 0 ~ "Unburned",
#                                          class == 1 ~ "Low_Sev",
#                                          class == 2 ~ "Mod_Sev",
#                                          class == 3 ~ "High_Sev")) |> 
#             select(plot_id, level, lyr_name, FireClass, metric, value) |> 
#             tidyr::pivot_wider(names_from = c(metric, FireClass, lyr_name), 
#                                values_from = value,
#                                names_glue = "{FireClass}_{lyr_name}_{metric}_c") |> 
#             select(!level)
#           
#         }
#         
#         if(exists("cbi_lsm_l") & exists("cbi_lsm_c")){
#           cbi_lsm <- full_join(cbi_lsm_l, cbi_lsm_c)
#         }
#       }
#     } #lscpmet
#   }
#   
#   #} #Env product CBI Fire
#   
#   return(list(FireMetrics = fire_buff_merge,
#               FireLscp_Class = cbi_lsm_c,
#               FireLscp_Landscape = cbi_lsm_l))
#   
# } ## function closure

## Test it

## -------------------------------------------------------------
##
## Begin Section: Chunking functions
##
## -------------------------------------------------------------

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: CAbioacoustics function
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Load the ROI
cabio_loc_query <- function(years = survey_years){
  
  ## Connect to CAbioacoustics DB
  CAbioacoustics::cb_connect_db()
  
  ## View the tables housed in the DB
  DBI::dbListTables(conn = conn)
  
  ## Retrieve the study area
  roi <- CAbioacoustics::cb_get_spatial(layer_name = "sierra_study_area")
  
  ## Grab the ARU location from the CAbioacoustics package
  ## Pull the spatial coordinates
  aru_df <-
    tbl(conn, "acoustic_field_visits") |>
    filter(survey_year %in% years & study_type == 'Sierra Monitoring' ) |>
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
      utm_zone,
      utme,
      utmn
    ) |>
    collect()
  
  aru_locs <- aru_df |>
    filter(!is.na(utm_zone)) |> 
    group_split(utm_zone) |>
    map_dfr(cb_make_aru_sf) |> 
    mutate(Cell_Unit = stringr::str_remove(deployment_name, "G[0-9]+_V[0-9]+_"))
  
  return(aru_locs)
  
  ## Remove connection
  CAbioacoustics::cb_disconnect_db()
}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: CBI Stacking Function
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cbi_stack_fun <- function(cbi_files){
  tryCatch({
    res <- terra::rast(cbi_files)
    return(res)
  }, error = function(e) {
    ## Error message
    message("Error:", e$message)
    
    ## if error is different spatial extents
    if(stringr::str_detect(e$message, "extents do not match")){
      
      ## Print a message
      cat("Extents do not match. Attempting to resolve. \n This could take a minute...")
      
      ## Fix extents
      cbi_list <- lapply(cbi_files, rast)
      ## Extents for all rasters
      l_extent <- lapply(cbi_list, ext)
      ## Make all extents bboxes
      cbi_bbox <- lapply(l_extent, function(x) vect(x))
      ## Union to find the largest extent
      big.ext <- Reduce(terra::union, cbi_bbox)
      ## Set the CRS
      crs(big.ext) <- crs(cbi_list[[1]])
      ## Extend all rasters to the same extent
      cbi_list <- lapply(cbi_list, function(x) extend(x, big.ext))
      ## Restack the rasters
      cbi_stack <- try({
        res <- do.call(c, cbi_list)
      }, silent = T)
      
      if(inherits(cbi_stack, "try-error")){
        stop("Check raster data. Could not stack.")
      } else {
        cat("Successful raster matching.")
        return(res)
      } #if error after alignment
      
    } #if error
    
  })}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Fire year interval stacking
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
int_ras_fun <- function(ras_stack,
                        intervals,
                        sum_int,
                        aru_locs = aru_locs){
  
  if(!is.null(intervals) & class(intervals) == "character"){
    int_len <- length(intervals)
    seq_list <- vector(mode = "list", length = int_len)
    int <- 0
    while(int < int_len){
      int <- int + 1
      int.t <- strsplit(intervals[int], "-")[[1]]
      if(length(int.t) > 1){
        s_seq <- int.t[[1]]
        e_seq <- int.t[[2]]
        seq_list[[int]] <- seq(s_seq, e_seq, by = 1)
      } else {seq_list[[int]] <- as.numeric(int.t)}
    }
    
    ## Check is the seq list is given in dates or integer years
    if(!any(sort(unlist(seq_list)) %in% names(ras_stack))){
      aru_years <- unique(aru_locs$survey_year)
      if(length(aru_years) > 1){
        
        ## Make sure all the years are integers
        years <- as.integer(aru_years)
        ## Find the number of years
        num_years <- length(years)
        ## Replicate the sequence list n times
        seq_list <- rep(list(seq_list), num_years)
        ## for each year replace the values to match the names of the fire data
        for(i in 1:num_years){
          seq_tmp <- seq_list[[i]]
          ## Back calculate the expected dates from the sequence of numbers
          ## Sort the unlisted seq_list vec
          seq_vec <- sort(unlist(seq_tmp))
          ## Find the year of interest
          years <- as.integer(aru_years[i])
          ## Find all years before the year of interest
          fire_years <- years - seq_vec
          ## Update the seq_list to the dates of the fires
          seq_list[[i]] <- lapply(seq_list[[i]], function(x) fire_years[x])
        }
      }
      
      if(length(aru_years) == 1){
        ## Back calculate the expected dates from the sequence of numbers
        ## Sort the unlisted seq_list vec
        seq_vec <- sort(unlist(seq_list))
        ## Find the year of interest
        years <- as.integer(aru_years)
        ## Find all years before the year of interest
        fire_years <- years - seq_vec
        ## Update the seq_list to the dates of the fires
        seq_list <- lapply(seq_list, function(x) fire_years[x])
      }
    } # fix year naming
    
    if(length(aru_years) == 1){
      ## Create the spatial products with the sequences for 1 year
      ras_list <- vector(mode = "list", length = int_len)
      for(i in 1:length(seq_list)){
        message(paste("Finding", sum_int, "value for Year:", aru_years[i], "For interval", j, "of", int_len))
        message("This is a single year workflow. Set by survey_years at the top-level.")
        seq_temp <- seq_list[[i]]
        ras_tmp <- ras_stack[[which(names(ras_stack) %in% seq_list[[i]])]]
        if(nlyr(ras_tmp) < 2){
          ras_list[[i]] <- ras_tmp
          names(ras_list[[i]]) <- names(ras_tmp)
        } else {
          ras_list[[i]] <- app(ras_tmp, fun = sum_int)
          names(ras_list[[i]]) <- paste0(names(ras_tmp)[1], "-", names(ras_tmp)[length(names(ras_tmp))])
          ras_int <- list(do.call(c, ras_list))
        }
      } ## Layer extraction
    }
    
    if(length(aru_years) > 1){
      ## Create the spatial products with the sequences for 1 year
      ras_list <- vector(mode = "list", length = int_len)
      ras_list <- rep(list(ras_list), num_years)
      for(i in 1:num_years){
        for(j in 1:int_len){
          print(paste("Finding", sum_int, "value for Year:", aru_years[i], "For interval", j, "of", int_len))
          seq_temp <- unlist(seq_list[[i]][j])
          ras_tmp <- ras_stack[[which(names(ras_stack) %in% seq_temp)]]
          if(nlyr(ras_tmp) < 2){
            ras_list[[i]][[j]] <- ras_tmp
            names(ras_list[[i]][[j]]) <- names(ras_tmp)
          } else {
            ras_list[[i]][[j]] <- app(ras_tmp, fun = sum_int)
            names(ras_list[[i]][[j]]) <- paste0(names(ras_tmp)[1], "-", names(ras_tmp)[length(names(ras_tmp))])
          }
        }
        
        ras_int <- vector(length(ras_list), mode = "list")
        for(i in 1:length(ras_list)){
          ras_int[[i]] <- do.call(c, ras_list[[i]])
        }
        
        
      } ## Layer extraction
    }
    
  } #FULL IF STATEMENT
  
  ## Return the stacked intervals
  return(ras_int)
}


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Landscape metrics function
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Landscape Metrics
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

calculate_landscape_metrics(ras_int = cbi_int,
                            ras_stack = ras_stack,
                            locations = current_locs,
                            years = years,
                            buff_size = buff_size,
                            id_col = id_col,
                            intervals = intervals)

## Check with consolidated version
calculate_landscape_metrics <- function(ras_int,
                                        ras_stack,
                                        locations,
                                        years,
                                        buff_size,
                                        id_col,
                                        intervals = NULL,
                                        metrics = c("lsm_c_pland")) {
  
  # Input validation
  if (!inherits(locations, "sf")) {
    stop("locations must be an sf object")
  }
  
  if (is.null(buff_size)) {
    stop("buff_size must be specified")
  }
  
  # Initialize empty list to store results for each year
  yearly_results <- vector("list", length(years))
  
  # Outer loop for years
  for(y in seq_along(years)) {
    current_year <- years[y]
    message(paste("Processing year:", current_year))
    
    # Filter locations for current year
    current_locs <- locations[locations$survey_year == current_year, ]
    current_locs <- current_locs[1:20,]
    
    if(nrow(current_locs) == 0) {
      message(paste("No locations found for year", current_year))
      next
    }
    
    # Inner loop for intervals within each year
    interval_results <- vector("list", nlyr(ras_int))
    
    for(i in seq_along(interval_results)) {
      # Get location IDs if specified
      loc_ids <- if (!is.null(id_col)) {
        as.vector(unlist(st_drop_geometry(current_locs[, id_col])))
      } else {
        NULL
      }
      
      # Select appropriate raster data
      current_raster <- if (is.null(intervals)) {
        message("Using individual years for landscape metrics")
        ras_stack
      } else {
        message("Using intervals for landscape metrics")
        ras_int[[i]]
      }
      
      # Calculate landscape metrics
      lsm <- landscapemetrics::sample_lsm(
        landscape = current_raster,
        y = current_locs,
        plot_id = loc_ids,
        shape = "circle",
        size = buff_size,
        what = metrics,
        return_raster = FALSE,
        verbose = FALSE,
        progress = TRUE
      )
      
      # Process results
      if (nrow(lsm) > 0) {
        # Add layer names
        lyr_map <- data.frame(
          layer = i,
          lyr_name = intervals[i]
        )
        lsm$lyr_name <- lyr_map$lyr_name
        
        # Process landscape-level metrics
        if ("landscape" %in% unique(lsm$level)) {
          lsm_landscape <- process_landscape_metrics(lsm, current_year)
        } else {lsm_landscape <- NULL}
        
        # Process class-level metrics
        if ("class" %in% unique(lsm$level)) {
          lsm_class <- process_class_metrics(lsm, current_year)
        } else {lsm_class <- NULL}
        
        # Combine results
        interval_results[[i]] <- combine_metrics(lsm_landscape, lsm_class)
      }
    }
    
    # Combine all intervals for this year
    yearly_results[[y]] <- do.call(cbind, interval_results)
  }
  
  # Combine all years
  final_results <- do.call(rbind, yearly_results)
  return(final_results)
}

# Helper functions remain the same
process_landscape_metrics <- function(lsm, year) {
  lsm %>%
    filter(level == "landscape") %>%
    select(plot_id, level, lyr_name, metric, value) %>%
    tidyr::pivot_wider(
      names_from = c(metric, lyr_name),
      values_from = value,
      names_glue = "{lyr_name}_{metric}_l"
    ) %>%
    select(!level) %>%
    mutate(Year = year)
}

process_class_metrics <- function(lsm, year) {
  fire_class <- c("Unburned", "Low_Sev", "Mod_Sev", "High_Sev")
  
  lsm %>%
    filter(level == "class") %>%
    mutate(
      FireClass = case_when(
        class == 0 ~ "Unburned",
        class == 1 ~ "Low_Sev",
        class == 2 ~ "Mod_Sev",
        class == 3 ~ "High_Sev"
      )
    ) %>%
    select(plot_id, level, lyr_name, FireClass, metric, value) %>%
    tidyr::pivot_wider(
      names_from = c(metric, FireClass, lyr_name),
      values_from = value,
      names_glue = "{FireClass}_{lyr_name}_{metric}_c"
    ) %>%
    select(!level) %>%
    mutate(Year = year)
}

combine_metrics <- function(landscape_metrics = NULL, class_metrics = NULL) {
  if (!is.null(landscape_metrics) && !is.null(class_metrics)) {
    full_join(landscape_metrics, class_metrics)
  } else if (!is.null(landscape_metrics)) {
    landscape_metrics
  } else if (!is.null(class_metrics)) {
    class_metrics
  }
}





## Function should estimate specific landscape metrics for the
## interval-based raster data for CBI
fire_lscp_fun <- function(ras_int = cbi_int, 
                          ras_stack = cbi_stack, 
                          locs = aru_locs,
                          years = survey_years,
                          buff_size = buff_size,
                          id_col = id_col,
                          metrics = c("lsm_c_pland")){
  
  for(i in 1:length(ras_int)){
    
    print(paste("Year:", years[i]))
    
    ## Pull the ARUs for a specific year
    locs_tmp <- locs |> 
      filter(survey_year == years[i])
    
    ## Check that we intend to calculate landscape metrics
    if(landscape_metrics == TRUE){
      
      ## Is there a buffer size?
      if(!is.null(buff_size)){ #| is.null(hex_size)){
        
        ## Check to make sure sf
        if(!is.null(id_col) & any(str_detect(class(locs_tmp), "sf"))){
          id.v <- as.vector(unlist(st_drop_geometry(locs_tmp[, id_col])))
        } else {
          id.v <- NULL
        }
        
        ## Individual years versus intervals for the fire data
        if(is.null(intervals)){ 
          message("Landscapemetrics done for individual years.")
          ras_lcpc <- ras_stack 
        } else { 
          ras_lcpc <- ras_int[[i]]
          message("Landscapemetrics done for intervals.")
        }
        
        ## Calculate LSCP
        lsm <- landscapemetrics::sample_lsm(landscape = ras_lcpc,
                                            y = locs_tmp,
                                            plot_id = id.v,
                                            shape = "circle",
                                            size = buff_size,
                                            what = metrics,
                                            #metric = "core",
                                            #level = "patch",
                                            return_raster = F,
                                            #type = "diversity metric",
                                            #classes_max = 3,
                                            verbose = FALSE,
                                            progress = T)
        
        ## Fix layer names
        lyr_map <- data.frame(layer = 1:nlyr(ras_lcpc), lyr_name = intervals)
        
        ## Merge these
        lsm <- left_join(lsm, lyr_map)
        
        ## level of lsm
        lsm_lev <- unique(lsm$level)
        
        if(any(lsm_lev == "landscape")){
          lsm_l <- lsm |> 
            filter(level == "landscape") |> 
            select(plot_id, level, lyr_name, metric, value) |> 
            tidyr::pivot_wider(names_from = c(metric, lyr_name), 
                               values_from = value,
                               names_glue = "{lyr_name}_{metric}_l") |> 
            select(!level) |>
            mutate(Year = years[i])
          
        }
        
        if(any(lsm_lev == "class")){
          
          fire_class <- c("Unburned", "Low_Sev", "Mod_Sev", "High_Sev")
          
          lsm_c <- lsm |> 
            filter(level == "class") |> 
            mutate(FireClass = case_when(class == 0 ~ "Unburned",
                                         class == 1 ~ "Low_Sev",
                                         class == 2 ~ "Mod_Sev",
                                         class == 3 ~ "High_Sev")) |> 
            select(plot_id, level, lyr_name, FireClass, metric, value) |> 
            tidyr::pivot_wider(names_from = c(metric, FireClass, lyr_name), 
                               values_from = value,
                               names_glue = "{FireClass}_{lyr_name}_{metric}_c") |> 
            select(!level) |> 
            mutate(Year = aru_years[i])
          
        }
        
        if(exists("lsm_l") & exists("lsm_c")){ lsm <- full_join(lsm_l, lsm_c) } 
        if (exists("lsm_l") & !exists("lsm_c")) { lsm <- lsm_l }  
        if (!exists("lsm_l") & exists("lsm_c")) { lsm <- lsm_c }
      }
    } else {stop("LandscapeMetrics calculation not selected.")}
    #} #lscpmet
    
    ## Grow the DF for multiple years
    if(i == 1){
      lsm_out <- lsm
    } else {
      lsm_out <- rbind(lsm_out, lsm)
    }
    
  }
  
  # ## Block for what to do with CBI Fire Data
  # if(landscape_metrics == TRUE){
  #   if(!is.null(buff_size)){ #| is.null(hex_size)){
  #     if(!is.null(id_col) & any(str_detect(class(aru_locs), "sf"))){
  #       id.v <- as.vector(unlist(st_drop_geometry(aru_locs[, id_col])))
  #     } else {
  #       id.v <- NULL
  #     }
  #     if(is.null(intervals)){ 
  #       message("Landscapemetrics done for individual years.")
  #       lcpc <- ras_stack 
  #     } else { 
  #       lcpc <- ras_int
  #       message("Landscapemetrics done for intervals.")
  #     }
  #     lsm <- landscapemetrics::sample_lsm(landscape = lcpc,
  #                                             y = locs,
  #                                             plot_id = id.v,
  #                                             shape = "circle",
  #                                             size = buff_size,
  #                                             what = c("lsm_c_pland",
  #                                                      "lsm_c_ed",
  #                                                      "lsm_c_lpi",
  #                                                      "lsm_l_shdi",
  #                                                      "lsm_l_msidi",
  #                                                      "lsm_l_msiei"),
  #                                             #metric = "core",
  #                                             #level = "patch",
  #                                             return_raster = F,
  #                                             #type = "diversity metric",
  #                                             #classes_max = 3,
  #                                             verbose = FALSE,
  #                                             progress = T)
  #     
  #     ## Fix layer names
  #     lyr_map <- data.frame(layer = 1:nlyr(lcpc), lyr_name = names(lcpc))
  #     
  #     ## Merge these
  #     lsm <- left_join(lsm, lyr_map)
  #     
  #     ## level of lsm
  #     lsm_lev <- unique(lsm$level)
  #     
  #     if(any(lsm_lev == "landscape")){
  #       lsm_l <- lsm |> 
  #         filter(level == "landscape") |> 
  #         select(plot_id, level, lyr_name, metric, value) |> 
  #         tidyr::pivot_wider(names_from = c(metric, lyr_name), 
  #                            values_from = value,
  #                            names_glue = "{lyr_name}_{metric}_l") |> 
  #         select(!level)
  #       
  #     }
  #     
  #     if(any(lsm_lev == "class")){
  #       
  #       fire_class <- c("Unburned", "Low_Sev", "Mod_Sev", "High_Sev")
  #       
  #       lsm_c <- lsm |> 
  #         filter(level == "class") |> 
  #         mutate(FireClass = case_when(class == 0 ~ "Unburned",
  #                                      class == 1 ~ "Low_Sev",
  #                                      class == 2 ~ "Mod_Sev",
  #                                      class == 3 ~ "High_Sev")) |> 
  #         select(plot_id, level, lyr_name, FireClass, metric, value) |> 
  #         tidyr::pivot_wider(names_from = c(metric, FireClass, lyr_name), 
  #                            values_from = value,
  #                            names_glue = "{FireClass}_{lyr_name}_{metric}_c") |> 
  #         select(!level)
  #       
  #     }
  #     
  #     if(exists("lsm_l") & exists("lsm_c")){
  #       lsm <- full_join(lsm_l, lsm_c)
  #     }
  #   }
  # } #lscpmet
  
  ## Return
  return(lsm_out)
  
}

fire_lscp_fun()

## -------------------------------------------------------------
##
## Begin Section: Fire Functions
##
## -------------------------------------------------------------
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Year of most recent fire
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Fire variables
## Time since fire
time_to_most_recent_fire <- function(cell_values, years) {
  # Check if there are any 1s in the cell
  if (any(cell_values > 0)) {
    # Find the most recent occurrence of 1 (event) - we look for the last occurrence of 1
    recent_event_index <- max(which(cell_values > 0), na.rm = TRUE)
    
    # Calculate the time since the most recent event (in years)
    recent_event_year <- years[recent_event_index]
    time_since_event_value <- as.numeric(max(years) - recent_event_year)  # Get current year and subtract
  } else {
    # If no event occurred, return NA
    time_since_event_value <- max(years) - min(years)
  }
  
  return(time_since_event_value)
}


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Fire frequency
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Fire freq
# Create the function to calculate time since the most recent event for each cell
fire_freq_calc <- function(cell_values, years) {
  
  # Check if there are any 1s in the cell
  if (any(cell_values > 0)) {
    
    # Find the most recent occurrence of 1 (event) - we look for the last occurrence of 1
    fire_years <- years[which(cell_values > 0)]
    num_fires <- length(fire_years)
    num_years <- length(years)
    
    fire_freq <- num_fires / num_years
    
  } else {
    # If no event occurred, return NA
    fire_freq <- 0
  }
  
  return(fire_freq)
}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Fire Return Interval
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Create the function to calculate fire return interval
## Return mean difference in time between fires (assumes yearly data)
## Returns NA under two conditions: Less than 2 fires produce NA
fire_return_int <- function(cell_values, years){
  # Check if there are any 1s in the cell
  if (any(cell_values > 0)) {
    
    # Find the most recent occurrence of 1 (event) - we look for the last occurrence of 1
    fire_years <- as.numeric(years[which(cell_values > 0)])
    
    if(length(fire_years) > 1){
      # Calculate the average time in-between fire events
      f_int <- diff(fire_years)
      ret_int <- mean(f_int)
    } else {
      ret_int <- NA
    }
  } else {
    # If no event occurred, return NA
    #ret_int <- length(years)
    ret_int <- NA
  }
  
  return(ret_int)
}


# ## -------------------------------------------------------------
# ##
# ## Begin Section: For Luca
# ##
# ## -------------------------------------------------------------
# ## *************************************************************
# ##
# ## Section Notes: Luca needs the following:
# ## Proportion of high-severity fire 2 years prior
# ## Proportion of high-severity fire 3-10 years prior
# ## Proportion of high-severity fire 11-35 years prior
# ## 
# ## For 2021 -> 2019-2020; 2010-2018; 1985-2009
# ##
# ## *************************************************************
# 
# ## 2021
# cbi_luca_2021 <- cbi_lsm |> 
#   select(plot_id, `High_Sev_1985-2009_pland_c`, `High_Sev_2010-2018_pland_c`, `High_Sev_2019-2020_pland_c`) |> 
#   mutate(across(2:4, ~if_else(is.na(.), 0, .)))
# 
# cbi_luca_2021 <- full_join(aru_locs, cbi_luca_2021, by = c("deployment_name" = "plot_id"))
# 
# cbi_luca_2021 <- cbi_luca_2021 |> 
#   rename(Prop_High_Sev_11_35yr = `High_Sev_1985-2009_pland_c`, Prop_High_Sev_3_10yr = `High_Sev_2010-2018_pland_c`, Prop_High_Sev_1_2yr = `High_Sev_2019-2020_pland_c`)
# 
# ## 2022
# cbi_luca_2022 <- cbi_lsm |> 
#   select(plot_id, `High_Sev_1986-2010_pland_c`, `High_Sev_2011-2019_pland_c`, `High_Sev_2020-2021_pland_c`) |> 
#   mutate(across(2:4, ~if_else(is.na(.), 0, .)))
# 
# cbi_luca_2022 <- full_join(aru_locs, cbi_luca_2022, by = c("deployment_name" = "plot_id"))
# 
# cbi_luca_2022 <- cbi_luca_2022 |> 
#   rename(Prop_High_Sev_11_35yr = `High_Sev_1986-2010_pland_c`, Prop_High_Sev_3_10yr = `High_Sev_2011-2019_pland_c`, Prop_High_Sev_1_2yr = `High_Sev_2020-2021_pland_c`)
# 
# ## 2023
# cbi_luca_2023 <- cbi_lsm |> 
#   select(plot_id, `High_Sev_1987-2011_pland_c`, `High_Sev_2012-2020_pland_c`, `High_Sev_2021-2022_pland_c`) |> 
#   mutate(across(2:4, ~if_else(is.na(.), 0, .)))
# 
# cbi_luca_2023 <- full_join(aru_locs, cbi_luca_2023, by = c("deployment_name" = "plot_id"))
# 
# cbi_luca_2023 <- cbi_luca_2023 |> 
#   rename(Prop_High_Sev_11_35yr = `High_Sev_1987-2011_pland_c`, Prop_High_Sev_3_10yr = `High_Sev_2012-2020_pland_c`, Prop_High_Sev_1_2yr = `High_Sev_2021-2022_pland_c`)
# 
# ## 2024
# cbi_luca_2024 <- cbi_lsm |> 
#   select(plot_id, `High_Sev_1988-2012_pland_c`, `High_Sev_2013-2021_pland_c`, `High_Sev_2022-2023_pland_c`) |> 
#   mutate(across(2:4, ~if_else(is.na(.), 0, .)))
# 
# cbi_luca_2024 <- cbi_luca_2024 |> 
#   rename(Prop_High_Sev_11_35yr = `High_Sev_1988-2012_pland_c`, Prop_High_Sev_3_10yr = `High_Sev_2013-2021_pland_c`, Prop_High_Sev_1_2yr = `High_Sev_2022-2023_pland_c`)
# 
# cbi_luca_2024 <- full_join(aru_locs, cbi_luca_2024, by = c("deployment_name" = "plot_id"))
# 
# ## Bind
# cbi_luca <- rbind(cbi_luca_2021, cbi_luca_2022, cbi_luca_2023, cbi_luca_2024)
# 
# cbi_luca |>
#   filter(survey_year == 2022) |> 
#   ggplot() + 
#   geom_sf(aes(color = Prop_High_Sev_1_2yr)) + 
#   scale_color_viridis_c()
# 
# cbi_luca <- cbi_luca |> 
#   st_transform(crs = 4326) |> 
#   mutate(lon = st_coordinates(cbi_luca)[,1],
#          lat = st_coordinates(cbi_luca)[,2]) |> 
#   st_drop_geometry() |> 
#   select(id:deploy_date, lat, lon, everything())
# 
# 
# #data.table::fwrite(cbi_luca, file = "C:/Users/srk252/Documents/ARU_2021_2024_PropHighSevFire_Luca.csv")
# 
# rthet_luca <- data.table::fread("C:/Users/srk252/Documents/ARU_2021_2024_ExtractedDataRelHet_Luca.csv")
# rthet_luca <- rthet_luca |> tidyr::as_tibble()
# 
# luca_da <- cbi_luca |> 
#   select(id, Prop_High_Sev_11_35yr:Prop_High_Sev_1_2yr) |> 
#   full_join(rthet_luca, by = "id") |> 
#   select(id, study_type:deploy_date, lon, lat, Prop_High_Sev_1_2yr, Prop_High_Sev_3_10yr, Prop_High_Sev_11_35yr, RelativeTempC, ThermalHet, MeanElevation)
# 
# luca_da |> 
#   filter(survey_year == 2022) |> 
#   st_as_sf(coords = c("lon", "lat"), crs = 4326) |> 
#   mapview::mapview(zcol = "Prop_High_Sev_1_2yr")
# 
# data.table::fwrite(luca_da, file = "C:/Users/srk252/Documents/ARU_2021_2024_FireTempElevationVars_Luca.csv")
# 
# 
# ## -------------------------------------------------------------
# ##
# ## Begin Section: Garrett
# ##
# ## -------------------------------------------------------------
# datgar <- cbi_lsm_out
# 
# write.csv(datgar, here("C:/Users/srk252/Documents/ARU_2021_2024_Garrett_PropFire.csv"))

