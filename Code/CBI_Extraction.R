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

## Environmental extraction function
aru_env_prep <- function(env_prod, # character vector of prod names
                         aru_locs, # Data.frame with coordinates
                         id_col = "Cell_Unit",
                         year_col = NULL,
                         x_col = "Long", # chr for x coordinate col name
                         y_col = "Lat", # chr for y coordinate col name
                         .crs = 4326, # default WGS84 for coordinates
                         des_out, # desired output
                         spat_ex, # chr for type (point, buff, hex)
                         buff_size = NULL, # vector of buffer sizes
                         hex_size = NULL,
                         intervals = NULL,
                         landscape_metrics = F,
                         ...
                         ){
  
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##
  ## Subsection: Load ROI
  ##
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## Load the ROI
  CAbioacoustics::cb_connect_db()
  
  ## View the tables housed in the DB
  DBI::dbListTables(conn = conn)
  
  ## Retrieve the study area
  roi <- CAbioacoustics::cb_get_spatial(layer_name = "sierra_study_area")
  
  ## Remove connection
  CAbioacoustics::cb_disconnect_db()
  
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##
  ## Subsection: Create SF object from the ARU coordinates
  ##
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  ## Create SF objects from points
  if(!any(class(aru_locs) %in% "sf")){
    
    aru_locs <- st_as_sf(aru_locs, 
                         coords = c(x_col, 
                                    y_col), 
                         crs = .crs,
                         ...)
  }
  
  
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##
  ## Subsection: CBI Categorical Extraction and prep
  ##
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  ## Find the data for extraction
  if(env_prod == "Fire_CBI"){
    
    ## Find files
    cbi_path <- "C:/Users/srk252/Documents/data_for_spencer/cbi_sierra_cat_rasters/"
    cbi_files <- list.files(cbi_path, full.names = T, pattern = "(cbi_cat_)(\\d{4})(*.tif$)")
    
    ## Read stack of CBI files
    cbi_stack <- tryCatch({
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
      
    })
    
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
    
    ## Reclassify the NAs to zero
    ## Is this okay?
    cbi_stack <- classify(cbi_stack, cbind(NA, 0))
    if(nlyr(cbi_stack) == length(cbi_files)){
      names(cbi_stack) <- str_extract(cbi_files, "\\d{4}")
    }
    
    ## If interval grouping for fire
    ## End goal here is a combine raster stack
    ## with the intervals aggregated for the CBI data
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
      
      ## Create the spatial products with the sequences
      cbi_list <- vector(mode = "list", length = int_len)
      for(i in 1:length(seq_list)){
        seq_temp <- seq_list[[i]]
        cbi_tmp <- cbi_stack[[seq_temp]]
        if(nlyr(cbi_tmp) < 2){
          cbi_list[[i]] <- cbi_tmp
          names(cbi_list[[i]]) <- names(cbi_tmp)
          } else {
          cbi_list[[i]] <- app(cbi_tmp, fun = max)
          names(cbi_list[[i]]) <- paste0(names(cbi_tmp)[1], "-", names(cbi_tmp)[length(names(cbi_tmp))])
        }
      }
    }
    
    ## Project the points
    aru_locs <- st_transform(aru_locs,
                             crs = crs(cbi_stack))
    
    ## Block for what to do with CBI Fire Data
    if(landscape_metrics == TRUE){
      if(is.null(buff_size) | is.null(hex_size)){
        if(!is.null(id_col) & any(str_detect(class(aru_locs), "sf"))){
          id.v <- as.vector(unlist(st_drop_geometry(aru_test[, id_col])))
        } else {
          id.v <- NULL
        }
      cbi_lsm <- landscapemetrics::sample_lsm(landscape = cbi_stack,
                                              y = aru_test,
                                              plot_id = id.v,
                                              shape = "circle",
                                              size = buff_size,
                                              metric = "pland",
                                              level = "class",
                                              return_raster = T,
                                              #type = "diversity metric",
                                              #classes_max = 3,
                                              verbose = FALSE,
                                              progress = T)
      
      ## Fix layer names
      lyr_map <- data.frame(layer = 1:nlyr(cbi_stack), lyr_name = names(cbi_stack))
      
      ## Merge these
      cbi_lsm <- left_join(cbi_lsm, lyr_map)
      
      ## DF manipulation - wide
      test <- cbi_lsm |> 
        mutate(FireClass = case_when(class == 0 ~ "Unburned",
                                     class == 1 ~ "Low_Sev",
                                     class == 2 ~ "Mod_Sev",
                                     class == 3 ~ "High_Sev")) |> 
        select(plot_id, level, lyr_name, FireClass, value) |> 
        tidyr::pivot_wider(names_from = FireClass, values_from = value) |>
        rowwise() |> 
        mutate(SumFlag = sum(c_across(Unburned:Low_Sev), na.rm = T)) |> 
        ungroup() |> 
        mutate(Unburned = if_else(is.na(Unburned) & ceiling(SumFlag) == 100, 0, Unburned),
               Low_Sev = if_else(is.na(Low_Sev) & ceiling(SumFlag) == 100, 0, Low_Sev),
               Mod_Sev = if_else(is.na(Mod_Sev) & ceiling(SumFlag) == 100, 0, Mod_Sev),
               High_Sev = if_else(is.na(High_Sev) & ceiling(SumFlag) == 100, 0, High_Sev))
        
      
      }
    }
    
    
  } #Env product CBI Fire
  
  
  
} ## function closure


## Scratch code
## Delete afterwards
r.hold <- c()
for(i in 1:length(cbi_files)){
  print(cbi_files[i])
  r <- rast(cbi_files[i])
  ext.r <- ext(r)
  r.hold <- c(r.hold, ext.r)
  
}
