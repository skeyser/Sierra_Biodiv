## ----------------------------------------------------------
##
## Script name: Cummulative Fire Patterns
##
## Script purpose: Estimate the total area burned from 1980 - 2020
## across the Sierra using CBI. Used for SI Fig X. in FEE paper.
##
## Author: Spencer R Keyser
##
## Date Created: 2026-05-05
##
## Email: skeyser@wisc.edu
##
## Github: https://github.com/skeyser
##
## -----------------------------------------------------------
##
## Notes:
##
##
## -----------------------------------------------------------

## Defaults
options(scipen = 6, digits = 4)

## -----------------------------------------------------------

## Package Loading
library(dplyr)
library(ggplot2)
library(data.table)
library(terra)
library(sf)
library(tidyr)
## -----------------------------------------------------------

## Create a function to estimate the cummulative proportion of the Sierra that has burned
cum_fire <- function(roi_path = "D:/GIS_Data/Sierra_ROI.shp", 
                     cbi_path = "D:/GIS_Data/CBI_Sierra/CBI_1985_2024_ZeroFilling_Stack_New.tif"){
  
  ## Load CBI from the raster
  cbi <- rast(cbi_path)
  
  ## Transform ROI to the matching CRS
  roi <- st_read(roi_path)
  roi <- st_transform(roi, crs = crs(cbi))
  
  ## Crop to start
  print("Cropping raster to ROI.")
  cbi <- crop(cbi, vect(roi), mask = T)
  
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##
  ## Subsection: Max CBI for total burn
  ##
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #cbi_max <- app(cbi, function(x) max(x, na.rm = T))
  print("Calculating max CBI for the whole series.")
  cbi_max <- max(cbi, na.rm = T)
  
  max_ex <- exactextractr::exact_extract(cbi_max, roi, 'frac', stack_apply = T, progress = T)
  
  max_df <- setNames(stack(max_ex), c("Value", "Variable")) |>
    as_tibble() |> 
    separate(col = Variable,
             into = c("Var", "Class"),
             sep = "_") |>
    mutate(Fire = case_when(Class == 0 ~ "Unburned",
                            Class == 1 | Class == 2 ~ "Low/Mod",
                            Class == 3 ~ "High",
                            TRUE ~ "Other")) |>
    group_by(Fire) |> 
    summarise(Value = sum(Value))
  
  ## Clean
  rm(cbi_max)
  gc()
  
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##
  ## Subsection: Cummulative total burn x severity
  ##
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  print("Moving window for cummulative burns.")
  
  ## Create temp dir
  temp_dir <- file.path(tempdir(), "cbi_cumulative")
  dir.create(temp_dir, showWarnings = FALSE)
  
  nyear <- nlyr(cbi)
  
  ## Stash the current max cum file
  current_max_file <- file.path(temp_dir, "current_max.tif")
  
  ## Initialization
  writeRaster(cbi[[1]], current_max_file, overwrite = TRUE)
  
  ## Create a list of extractions
  cum_extractions <- list()
  
  ## Make a cummulative stack
  cum_extractions[[1]] <- exactextractr::exact_extract(rast(current_max_file), roi,
                                                       'frac', stack_apply = T, progress = F)
  
  ## Loop the years
  ## Different than the total max bc we get a moving total of burns
  for(i in 2:nyear){
    print(paste("On year", i, "of", nyear))
    
    ## Load cumulative raster from disk
    current_max <- rast(current_max_file)
    
    ## New cumulative max file
    temp_max_file <- file.path(temp_dir, "temp_max.tif")
    
    ## Compute directly to disk
    max(current_max, cbi[[i]],
        filename = temp_max_file,
        overwrite = TRUE,
        na.rm = TRUE)
    
    file.rename(temp_max_file, current_max_file)
    
    # Extract and save
    cum_extractions[[i]] <- exactextractr::exact_extract(rast(current_max_file), roi, 'frac', 
                                                         stack_apply = T, progress = F)
    
    ## Free up memory
    rm(current_max)
    print("Garbage collection")
    gc()
    
  }
  
  ## Clean temp directory
  unlink(temp_dir, recursive = TRUE)
  
  ## Process extractions
  for(i in 1:nyear) {
    year_name <- names(cbi)[i]
    temp_df <- setNames(stack(cum_extractions[[i]]), c("Value", "Variable"))
    temp_df$Year <- str_extract(year_name, "\\d{4}")  # Extract year from layer name
    cum_extractions[[i]] <- temp_df
  }
  
  # Combine all extractions
  cum_df <- bind_rows(cum_extractions) |>
    as_tibble() |>
    mutate(Class = as.numeric(str_extract(Variable, "\\d+"))) |>
    filter(!is.na(Class)) |>
    mutate(Fire = case_when(Class == 0 ~ "Unburned",
                            Class == 1 | Class == 2 ~ "Low/Mod",
                            Class == 3 ~ "High",
                            TRUE ~ "Other")) |>
    group_by(Fire, Year) |>
    summarise(Value = sum(Value, na.rm = T), .groups = "drop") |>
    group_by(Year) |>
    mutate(Check = sum(Value, na.rm = T)) |>
    ungroup()
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##
  ## Subsection: Annual Burn by severity
  ##
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  print("Annual proportion burned.")
  ## Fraction of cells that are in each class
  ## Find total cells that are not NA
  prop_ex <- exactextractr::exact_extract(cbi, roi, 'frac', stack_apply = T, progress = T)
  
  ## Package as a dataframe
  prop_df <- setNames(stack(prop_ex), c("Value", "Variable")) |>
    as_tibble() |> 
    separate(col = Variable,
             into = c("Var", "Class", "Year"),
             sep = "_|\\.") |>
    mutate(Fire = case_when(Class == 0 ~ "Unburned",
                            Class == 1 | Class == 2 ~ "Low/Mod",
                            Class == 3 ~ "High",
                            TRUE ~ "Other")) |> 
    group_by(Fire, Year) |> 
    select(-Class) |> 
    summarise(Value = sum(Value, na.rm = T)) |> 
    ungroup() |> 
    group_by(Year) |> 
    mutate(Check = sum(Value, na.rm = T)) |> 
    ungroup()
  
  ## Return output
  print("Returning output.")
  return(list(Max = max_df,
              Cum = cum_df,
              Annual = prop_df))
  
}


## Run the function
fire_sum <- cum_fire()
str(fire_sum)
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Summaries and Plots
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Plot changes in region-wide fire data
ann.p <- ggplot(data = fire_sum$Annual |> mutate(Fire = factor(Fire, levels = c("Unburned",
                                                                       "Low/Mod",
                                                                       "High"))), 
       aes(x = Year, y = Value*100, 
           group = Fire, color = Fire)) + 
  geom_point(size = 2) + 
  geom_line() +
  scale_color_manual(values = c("Unburned" = "#313695",
                                "Low/Mod" = "#fee090",
                                "High" = "#a50026"),
                     labels = c("Unburned",
                                "Low/Moderate",
                                "High")) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90),
        axis.text = element_text(family = "sans", size = 12),
        axis.title = element_text(family = "sans", size = 12),
        legend.text = element_text(family = "sans", size = 12),
        legend.title = element_text(family = "sans", size = 12)) + 
  ylab("Proportion (Annual)") + 
  xlab("Year") + 
  labs(color = "Fire Severity")

cum.p <- ggplot(data = fire_sum$Cum |> mutate(Fire = factor(Fire, levels = c("Unburned",
                                                                    "Low/Mod",
                                                                    "High"))), 
       aes(x = Year, y = Value*100, 
           group = Fire, color = Fire)) + 
  geom_point(size = 2) + 
  geom_line() +
  scale_color_manual(values = c("Unburned" = "#313695",
                                "Low/Mod" = "#fee090",
                                "High" = "#a50026"),
                     labels = c("Unburned",
                                "Low/Moderate",
                                "High")) + 
  # geom_hline(data = fire_sum$Max, aes(color = Fire, yintercept = Value * 100),
  #            linetype = "dashed", size = 1) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90),
        axis.text = element_text(family = "sans", size = 12),
        axis.title = element_text(family = "sans", size = 12),
        legend.text = element_text(family = "sans", size = 12),
        legend.title = element_text(family = "sans", size = 12)) + 
  ylab("Proportion (Cumulative)") + 
  xlab("Year") + 
  labs(color = "Fire Severity")

p <- ann.p / cum.p + plot_layout(guides = "collect") + plot_annotation(tag_levels = "A")

ggsave(plot = p, filename = here("./Figures/R1/PropSierraBurn.png"),
       height = 8, width = 8, dpi = 600)
