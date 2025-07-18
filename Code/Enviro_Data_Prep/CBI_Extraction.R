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

# ## Pull in the ARU meta data for point locations
# 
# ## Take a look at CBI data
# cbi <- rast("c:/Users/srk252/Documents/data_for_spencer/cbi_sierra_cat_rasters/cbi_cat_2021.tif")
# plot(cbi)
# 
# ## ARU metadata
# ## Take from Jays package?
# ## For now just use the one from Connor
# aru_meta <- readr::read_csv(here("Data/ARU_120m.csv"))
# names(aru_meta)
# glimpse(aru_meta)
# aru_meta$Cell_Unit <- paste0(aru_meta$cell_id, "_", aru_meta$unit_numbe)
# aru_locs <- aru_meta |> 
#   select(Cell_Unit, survey_yea, Long = X, Lat = Y)
# 
# aru_locs <- st_as_sf(aru_locs, coords = c("Long", "Lat"), crs = 4326)


source(here("./Code/Enviro_Data_Prep/CBI_Extract_Funs.R"))

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
                          custom_locs = NULL, # Data.frame with coordinates for custom locations                          
                          survey_years = c(2021
                                           #2022,
                                           #2023,
                                           #2024
                                           ), # Survey year is only applicable for locs_from_cabio = TRUE
                          id_col = "deployment_name",
                          year_col = NULL,
                          x_col = "Long", # chr for x coordinate col name
                          y_col = "Lat", # chr for y coordinate col name
                          .crs = 4326, # default WGS84 for coordinates
                          des_out, # desired output
                          spat_ex, # chr for type (point, buff, hex)
                          buff_size = 120, # vector of buffer sizes
                          intervals = c("1-5", "6-10", "11-35"),
                          landscape_metrics = T
){
  
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##
  ## Subsection: CAbioacoustics Spatial
  ##
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  if(locs_from_cabio){
    aru_locs <- cabio_loc_query(years = survey_years)
  } else if (!is.null(custom_locs)) {
    aru_locs <- custom_locs
  } else {
    stop("No locations provided. Either set locs_from_cabio = TRUE or provide locations for custom_locs.")
  }
  
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
    
    message("CBI rasters are fixed and exist in directory. Loading the fixed stack.")
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
  if(!is.null(intervals)){
  cbi_int <- int_ras_fun(ras_stack = cbi_stack,
                         intervals = intervals,
                         sum_int = "max",
                         locations = aru_locs)
  } else {
    message("No intervals specified. Skipping interval binning.")
  }
  
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
      fire_sev <- exactextractr::exact_extract(cbi_sev[[1]],
                                               aru_buffer,
                                               fun = c("mean", "stdev"),
                                               append_cols = id_col)
    } else {
      fire_sev <- terra::extract(cbi_sev,
                                 aru_locs,
                                 bind = T) |> st_as_sf()
    }
    colnames(fire_sev)[colnames(fire_sev) != id_col] <- paste0("Fire_Sev_", gsub("[[:punct:]]", "_", colnames(fire_sev)[colnames(fire_sev) != id_col]))
    
  } else {
    fire_sev <- NULL
  }
  
  
  
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##
  ## Subsection: Calculate fire variables
  ##
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ras_years <- as.numeric(names(cbi_stack))
  
  if("time_since_fire" %in% fire_prod){
    if(!is.null(intervals)){
      message("Intervals are set. Time since fire doesn't accept intervals...output will be generated from single years.")
    }
    
    ## Time to most recent fire
    time_since_fire <- terra::app(cbi_stack, 
                                  fun = function(x) time_to_most_recent_fire(cell_values = x,
                                                                             years = ras_years))
  }
  
  if("fire_freq" %in% fire_prod){
    if(!is.null(intervals)){
      message("Intervals are set. Fire frequency doesn't accept intervals...output will be generated from single years.")
    } 
    ## Fire frequency (num fires/total record length)
    fire_freq <- terra::app(cbi_stack, 
                            fun = function(x) fire_freq_calc(cell_values = x,
                                                             years = ras_years))
  }
  
  if("fire_ret_int" %in% fire_prod){
    if(!is.null(intervals)){
      message("Intervals are set. Fire return interval doesn't accept intervals.,,output will be generated from single years.")
    } 
    ## Fire return interval (mean of time between successive fires)
    fire_return_int <- terra::app(cbi_stack, 
                                  fun = function(x) fire_return_int(cell_values = x, 
                                                                    years = ras_years))
  }
  
  ## Report the status of s2 geometry for convenience
  if(!is.null(buff_size)){
    if(sf_use_s2() == T){
      message("s2 geometry enabled. Buff_size interpreted as meters.")
    } else {
      message("s2 disabled, if locations are geodetic (lat/lon) units interpretted as degrees.")
    }
  }
  
  
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##
  ## Subsection: Buffer Extraction from Custom Fire Variables
  ##
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  if(any(fire_prod %in% c("time_since_fire", "fire_freq", "fire_ret_int"))){
    variable_name <- fire_prod[!fire_prod == "fire_severity"]
    fire_buff_out <- vector(mode = "list", length = length(variable_name))
    
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
  } else {
    fire_buff_merge <- NULL
  }
  
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##
  ## Subsection: Landscape Metrics
  ##
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## Function should estimate specific landscape metrics for the
  ## interval-based raster data for CBI
  if(landscape_metrics){
    fire_lscp_out <- fire_lscp_fun(ras_int = cbi_int,
                                   ras_stack = cbi_stack,
                                   locs = aru_locs,
                                   buff_size = buff_size,
                                   id_col = id_col,
                                   metrics = c("lsm_c_pland"))
  } else {
    fire_lscp_out <- NULL
  }

  return(list(FireSeverity = fire_sev,
              FireMetrics = fire_buff_merge,
              FireLscp = fire_lscp_out))

} ## function closure


buffSize <- 120

fire_sev21 <- aru_fire_prep(fire_prod = c("fire_severity"),
                            locs_from_cabio = TRUE,
                            survey_years = 2021,
                            intervals = c("1-5", "6-10", "11-35"),
                            id_col = "deployment_name",
                            buff_size = buffSize,
                            landscape_metrics = F
)

fire_sev21 <- fire_sev21$FireSeverity


## Histogram for the fire severity
nrow(fire_sev21$FireSeverity)
hist(fire_sev21$FireSeverity$Fire_Sev_mean_2016_2020)
hist(fire_sev21$FireSeverity$Fire_Sev_mean_2011_2015)
hist(fire_sev21$FireSeverity$Fire_Sev_mean_1986_2010)

## Bring in the dataframe from the ARU stuff
aru_meta <- read.csv(here("./Data/ARU_120m_New.csv")) |> 
  mutate(deployment_name = paste(group_id, visit_id, cell_id, unit_numbe, sep = "_"))

aru_meta <- aru_meta |> 
  filter(deployment_name %in% fire_sev21$FireSeverity$deployment_name) |> 
  left_join(fire_sev21$FireSeverity)

write.csv(fire_sev21, file = here("./Data/FireSeverity2021_MeanStDev_AllARUs.csv"))

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Checking Sensitivity to multiple buffer sizes
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## 120 meter buffer size used in our analyses
buffSize <- 120

fire_sev21_120 <- aru_fire_prep(fire_prod = c("fire_severity"),
                                 locs_from_cabio = TRUE,
                                 survey_years = 2021,
                                 intervals = c("1-5", "6-10", "11-35"),
                                 id_col = "deployment_name",
                                 buff_size = buffSize,
                                 landscape_metrics = F
)
fire_sev21_120 <- fire_sev21_120$FireSeverity


## 500 meter buffer size on fire severity
buffSize <- 500

fire_sev21_500 <- aru_fire_prep(fire_prod = c("fire_severity"),
                            locs_from_cabio = TRUE,
                            survey_years = 2021,
                            intervals = c("1-5", "6-10", "11-35"),
                            id_col = "deployment_name",
                            buff_size = buffSize,
                            landscape_metrics = F
)
fire_sev21_500 <- fire_sev21_500$FireSeverity

## 1200 meter buffer size on fire severity
buffSize <- 1200

fire_sev21_1200 <- aru_fire_prep(fire_prod = c("fire_severity"),
                                 locs_from_cabio = TRUE,
                                 survey_years = 2021,
                                 intervals = c("1-5", "6-10", "11-35"),
                                 id_col = "deployment_name",
                                 buff_size = buffSize,
                                 landscape_metrics = F
)
fire_sev21_1200 <- fire_sev21_1200$FireSeverity


## Check correlation between the three sizes
cor(fire_sev21_120$Fire_Sev_mean_2016_2020, fire_sev21_500$Fire_Sev_mean_2016_2020)
cor(fire_sev21_120$Fire_Sev_mean_2016_2020, fire_sev21_1200$Fire_Sev_mean_2016_2020)
cor(fire_sev21_500$Fire_Sev_mean_2016_2020, fire_sev21_1200$Fire_Sev_mean_2016_2020)

cor(fire_sev21_120$Fire_Sev_mean_2011_2015, fire_sev21_500$Fire_Sev_mean_2011_2015)
cor(fire_sev21_120$Fire_Sev_mean_2011_2015, fire_sev21_1200$Fire_Sev_mean_2011_2015)
cor(fire_sev21_500$Fire_Sev_mean_2011_2015, fire_sev21_1200$Fire_Sev_mean_2011_2015)

cor(fire_sev21_120$Fire_Sev_mean_1986_2010, fire_sev21_500$Fire_Sev_mean_1986_2010)
cor(fire_sev21_120$Fire_Sev_mean_1986_2010, fire_sev21_1200$Fire_Sev_mean_1986_2010)
cor(fire_sev21_500$Fire_Sev_mean_1986_2010, fire_sev21_1200$Fire_Sev_mean_1986_2010)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Fire trends
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fire_ras <- rast("C:/Users/srk252/Documents/GIS_Data/CBI_Sierra/CBI_1985_2024_ZeroFilling_Stack.tif")
#fire_ras <- fire_ras[[!names(fire_ras) %in% c("2022", "2023")]]

## Get the frequency table
freq_table <- freq(fire_ras, bylayer = TRUE)

## Cell area (30mx30m)
cell_area <- prod(res(fire_ras))

## Area calculation
## Total count * area
freq_table$Area <- freq_table$count * cell_area

## Figure to illustrate the change in total area burned and burned
## across frequency classes
fire_freq <- freq_table |> 
  filter(value != 0) |>
  mutate(Year = rep(1985:2023, each = 3)) |> 
  mutate(value = case_when(value == 1 ~ "Low",
                           value == 2 ~ "Moderate",
                           value == 3 ~ "High")) |> 
  mutate(value = factor(value, levels = c("Low",
                                             "Moderate",
                                             "High"))) |> 
  mutate(Year = as.factor(Year)) |> 
  mutate(AreaKm = Area/1000) 

fire_year_plot <- ggplot(data = fire_freq) + 
  geom_segment(aes(x = Year,
                   y = 0, yend = AreaKm,
                   color = value)) +
  geom_point(aes(x = Year, y = AreaKm, 
                 color = value),
             size = 2) +
  geom_smooth(aes(x = Year, y = AreaKm, group = value, color = value, fill = value)) + 
  scale_color_brewer(type = "qual", palette = "Dark2") +
  scale_fill_brewer(type = "qual", palette = "Dark2") +
  theme_bw() +
  theme(axis.text = element_text(angle = 90, vjust = 0)) +
  scale_y_continuous(expand = c(0, 0)) + 
  xlab("Year") + 
  ylab(expression("Burned Area ln(KM"^2*")"))


fire_loglin <- lm(log(AreaKm) ~ value + as.numeric(Year), data = fire_freq)
fire_lin <- lm(AreaKm ~ value + as.numeric(Year), data = fire_freq)
summary(fire_loglin)
(exp(fire_loglin$coefficients[4]) - 1) * 100  # Percentage change

plot(fire_loglin, which = 3)
plot(fire_lin, which = 3)

library(ggeffects)

# Get predicted values
model_pred <- ggpredict(fire_loglin, terms = c("Year", "value"))

# Plot with back-transformed values
ggplot(model_pred, aes(x = x, y = predicted, color = group)) +
  geom_line(linetype = "dashed") +
  geom_ribbon(aes(ymin = conf.low, 
                  ymax = conf.high, 
                  fill = group), alpha = 0.1) +
  labs(x = "Year",
       y = expression("Area (km"^2*")"),
       color = "Severity",
       fill = "Severity") +
  theme_bw()

# With raw data
ggplot(model_pred, aes(x = x, y = predicted, color = group)) +
  geom_point(data = fire_freq, aes(x = Year, y = AreaKm, color = value), alpha = 0.3) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, 
                  ymax = conf.high, 
                  fill = group), alpha = 0.1) +
  labs(x = "Year",
       y = expression("Area (km"^2*")"),
       color = "Severity",
       fill = "Severity") +
  theme_bw()
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

