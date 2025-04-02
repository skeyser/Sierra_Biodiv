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
## Begin Section: Processing species detections
##
## -------------------------------------------------------------
## *************************************************************
##
## Section Notes: Source ac_det_filter() from Sierra_functions.R
## 
## To-do:
## 1. Streamline the multi-year species data processing
## 2. Unify with spatial coordinates from CAbioacoustics
## 3. Write yearly files for community data (summarized; naive biodiv. models)
## 4. Write yearly files for community data (date-detection; occ models)
##
## *************************************************************


## Source the ac_det_filter function
source(here("./Code/Acoustic_Data_Prep/Sierra_functions.R"))

## threshold file
thresh <- readr::read_csv("./Data/Thresholds_2021_20230309.csv")

## Run the function
aru_det_file_gen(det_dir = "C:/Users/srk252/Documents/Rprojs/Sierra_Biodiv/Data/Detections_By_Species/",
                 det_years = c("2021"),
                 seas_format = F,
                 seas_outdir = "C:/Users/srk252/Documents/Rprojs/Sierra_Biodiv/Data/Generated_DFs/Seasonal_Summaries/",
                 occ_format = T,
                 occ_outdir = "C:/Users/srk252/Documents/Rprojs/Sierra_Biodiv/Data/Generated_DFs/Occ_Mod_Data/Flocker/",
                 eff_file = T,
                 coord_link = T,
                 d_thresh = thresh,
                 thresh_scale = "Conf",
                 thresh_cut = "90",
                 time_format = "ymd",
                 no_dets = 2,
                 binary = F,
                 date_range = c("2021-06-01", "2021-06-30"),
                 eff_site_name = "Cell_U",
                 eff_filter = 10,
                 verbose = F)


## Commented out 11/13
## Revisit to see if it's okay to delete


# ## # of detections
# lapply(sp.det.list, function(x) paste(names(x), "detections:", sum(x[,2:ncol(x)])))
# paste(names(sp.det.list[1]), "detections:", sum(sp.det.list[[1]][2:ncol(sp.det.list[[1]])]))
# paste(names(sp.det.list[2]), "detections:", sum(sp.det.list[[2]][2:ncol(sp.det.list[[2]])]))
# paste(names(sp.det.list[3]), "detections:", sum(sp.det.list[[3]][2:ncol(sp.det.list[[3]])]))
# 
# ## Sum the detections in species richness DF
# cols_to_sum <- 2:ncol(sp.det.list[[1]])
# spec.rich <- Reduce(function(x,y){
#   x[,cols_to_sum] <- x[,cols_to_sum] + y[,cols_to_sum]
#   x
# }, sp.det.list)
# 
# ## If we want to retrieve species richness for each site across the entire season
# season_pool <- function(x){
# season_occ <- x |> 
#   mutate(row_sum = rowSums(across(where(is.numeric)))) |> 
#   select(where(~ !is.numeric(.)), row_sum) |> 
#   mutate(row_sum = ifelse(row_sum > 0, 1, row_sum)) |> 
#   rename(Season_Occ = row_sum) 
# 
# return(season_occ)
# }
# 
# seasonal_occ <- lapply(sp.det.list, season_pool)
# 
# join_column <- "Cell_Unit"
# 
# # Join the data frames in the list and rename columns
# seasonal_df <- reduce(names(seasonal_occ), function(acc, name) {
#   df <- seasonal_occ[[name]]
#   
#   # Rename columns by appending the data frame name
#   colnames(df)[-which(colnames(df) == join_column)] <- paste(name, colnames(df)[-which(colnames(df) == join_column)], sep = "_")
#   colnames(df) <- gsub("_Season_Occ", "", colnames(df))
#   
#   # Perform the left join
#   if (is.null(acc)) {
#     return(df)
#   } else {
#     return(left_join(acc, df, by = join_column))
#   }
# }, .init = NULL)
# 
# # Print the result
# print(seasonal_df)
# 
# ## `seasonal_df` now captures the naive occupancy for each species
# ## across the Sierra ARUs and can be used as a site x species matrix
# seasonal_df$Survey_Year <- 2022
# seasonal_df <- seasonal_df |> 
#   select(Cell_Unit, Survey_Year, everything())
# colnames(seasonal_df)
# 
# ## -------------------------------------------------------------
# ##
# ## Begin Section: ARU spatial data
# ##
# ## -------------------------------------------------------------
# 
# ## *************************************************************
# ##
# ## Section Notes: We now have the ARU data formatted into an ARU x Date matrix 
# ## with cells corresponding to the BirdNET confidence scores 
# ## (thresholded to the maximum detection score per day >= the estimated BirdNET thresholds). 
# ## The cells now represent naive occupancy for the species of interest. 
# ## Next, we will want to assign the ARUs a spatial location using spatial references for each ARU.
# # 
# # Spatial data can be accessed via the `CAbioacoustics` package under the table "acoustic_field_visits".
# ##
# ## *************************************************************
# 
# ## Connect to the database
# ## Need to be connected to VPN for this to connect
# cb_connect_db()
# 
# ## Set some values for filters
# syear <- 2021
# eyear <- 2024
# 
# ## SF needs to be loaded to execute
# ownership <- c('any')
# cell_list <- cb_cells_by_ownership(ownership)
# 
# ## Study type
# study <- c('Sierra Monitoring')
# 
# ## Query and pull the table we want
# deployments_df <- 
#   conn |> 
#   dplyr::tbl('acoustic_field_visits') |> 
#   dplyr::filter(
#     is_invalid == 0,
#     cell_id %in% cell_list,
#     study_type %in% study,
#     survey_year >= syear,
#     survey_year <= eyear
#   ) |> 
#   collect() |> 
#   dplyr::select(deployment_name, survey_year, matches("utm"))
# 
# ## Disconnect from the DB
# cb_disconnect_db()
# 
# ## SF creation
# deployments_sf <- 
#   deployments_df |> 
#   group_split(utm_zone) |> 
#   map_dfr(cb_make_aru_sf)
# 
# ## Plot with MV
# ## hexes
# hexes_map <-
#   cb_get_spatial('sierra_hexes') |>
#   mapview(layer.name = 'Hexes')
# 
# ## deployments
# deployments_map <-
#   deployments_sf |>
#   mutate(survey_year = as.factor(survey_year)) |>
#   mapview(zcol = 'survey_year', layer.name = 'Survey year')
# 
# ## combine
# hexes_map + deployments_map
# 
# 
# ## Link Species Detections with ARU locations - From Jay's DB
# glimpse(seasonal_df)
# head(seasonal_df$Cell_Unit)
# length(unique(seasonal_df$Cell_Unit))
# 
# 
# glimpse(deployments_sf)
# 
# ## Take only the 2021 spatial data
# dep21 <- deployments_sf |> 
#   filter(survey_year == 2022) |> 
#   mutate(Cell_Unit = stringr::str_remove(deployment_name, "G[0-9]+_V[0-9]+_")) |> 
#   st_transform(crs = 4326) |> 
#   dplyr::mutate(lon = sf::st_coordinates(geometry)[,2],
#                 lat = sf::st_coordinates(geometry)[,1]) |> 
#   st_drop_geometry() |>
#   ## Remove
#   select(Cell_Unit, deployment_name, survey_year, lon, lat) |> 
#   distinct()
# 
# ## *************************************************************
# ##
# ## Section Notes:
# ## Check the number of unit cell_unit combinations
# ## Check duplications in the cell_units
# ## Duplicates arise via different deployment_ids but not cell_unit ids
# ## V1 vs V2 in the deployment_id
# ## Step3 data only contains the cell_unit ID
# ## Deployments from owl DB have 1730 unique deployment_ids
# ## ATFL example species from 2021 has 1652 individual ARU deployments
# ##
# ## *************************************************************
# 
# nrow(dep21)
# length(unique(dep21$deployment_name))
# length(unique(dep21$Cell_Unit))
# dup_unit <- dep21[which(duplicated(dep21$Cell_Unit)),]$Cell_Unit
# dup21 <- dep21[dep21$Cell_Unit %in% dup_unit, ]
# dep21 <- distinct(dep21)
# 
# ## Match the ATFL data with the 2021 survey data
# seasonal_df_meta <- seasonal_df |>
#   mutate(Cell_Unit = ifelse(
#     stringr::str_detect(string = Cell_Unit, pattern = "C[0-9]{3}"),
#     gsub(pattern = "(C)([0-9]{3})(_U[0-9]+)$", replacement = "\\10\\2\\3", x = Cell_Unit),
#     Cell_Unit
#   )) |> 
#   left_join(dep21) |> 
#   select(Cell_Unit, deployment_name, survey_year, lon, lat, everything()) 
#   
# missing <- dep21[!dep21$Cell_Unit %in% seasonal_df_meta$Cell_Unit,]
# # 
# # atfl.sf.ok <- atfl.sf[which(!is.na(atfl.jay$survey_year)),]
# # atfl.sf.ok <- atfl.sf.ok |> 
# #   dplyr::select(Cell_Unit, deployment_name, survey_year, lon, lat, everything()) |> 
# #   rowwise() |> 
# #   mutate(sum = sum(c_across("2021-03-01":"2021-08-30"))) |> 
# #   ungroup() |> 
# #   dplyr::select(Cell_Unit, deployment_name, survey_year, lon, lat, sum) |> 
# #   st_as_sf(coords = c("lon", "lat"), crs = 4326)
# # 
# # ggplot(data = atfl.sf.ok) + 
# #   geom_sf(aes(color = sum))
# 
# 
# 
# ## -------------------------------------------------------------
# ##
# ## Begin Section: 2021 ARU Metadata
# ##
# ## -------------------------------------------------------------
# ## *************************************************************
# ##
# ## Section Notes:
# ## Data file with locations from Connor. 
# ## It seems like there are some inconsistencies with the leading zero 
# ## in the cell ID name. 
# ## My guess is that importing the files into excel dropped the leading zeros 
# ## since all of the issues are occurring with cells that have 3 numbers when they should have 4. 
# ## Add in a leading zero if the string only has 3 numbers after the "C" and let's see what happens.
# ##
# ## The 2021 ARU meta data file from Connor is specific to 2021, where the DB from Jay has the location
# ## information for each year. 
# ## 
# ## For analyses, I plan to extract/create variables on interest so I really just need to spatial
# ## locations and the extractions can be done on my part.
# ##
# ## *************************************************************
# 
# ## Load in the ARU metadata
# aru_meta <- readr::read_csv(here("Data/ARU_120m.csv"))
# names(aru_meta)
# glimpse(aru_meta)
# aru_meta$Cell_Unit <- paste0(aru_meta$cell_id, "_", aru_meta$unit_numbe)
# 
# ## Take a couple of things for testing
# aru_meta <- aru_meta |>
#   dplyr::select(Cell_Unit, survey_yea, X, Y, topo_elev, standage_f3_mn, contains("fire"))
# 
# ## Fixing the leading zero issue
# seasonal_df <- seasonal_df |>
#   mutate(Cell_Unit = ifelse(
#     stringr::str_detect(string = Cell_Unit, pattern = "C[0-9]{3}"),
#     gsub(pattern = "(C)([0-9]{3})(_U[0-9]+)$", replacement = "\\10\\2\\3", x = Cell_Unit),
#     Cell_Unit
#   ))
# 
# ## Attempt to rejoin now
# seasonal_df_meta <- seasonal_df |> 
#   left_join(aru_meta) |> 
#   select(Cell_Unit, 
#          Year = survey_yea, 
#          Long = X, 
#          Lat = Y,
#          everything())
# 
# ## ARU Meta
# aru_filter <- aru_meta |> 
#   filter(Cell_Unit %in% seasonal_df_meta$Cell_Unit)
# 
# ## Write the DF
# data.table::fwrite(seasonal_df_meta,
#                    here::here("Data/Generated_DFs/2021_SeasonalSpeciesMat.csv"))

