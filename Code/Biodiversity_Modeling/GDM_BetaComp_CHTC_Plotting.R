## -------------------------------------------------------------
##
## Script name: GDM Summaries
##
## Script purpose: GDM for Sierra Bioacoustics Latent Z using
## traditional Ferrier GDM approach with model fits and VI
## from CHTC runs
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
library(dplyr)
library(ggplot2)
library(here)
library(purrr)
library(tidyr)
library(sf)

## GDM
library(gdm)
library(adespatial)

## -------------------------------------------------------------

## Load the model fits
gdm.fits <- list.files(path = "R:/Users/skeyser/Postdoc/GDM_Output/", full.names = TRUE)
test <- readRDS(gdm.fits[1])

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Variable Importance
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
str(test)
VI <- test$varImp[[3]]
varImp <- VI$`Predictor Importance`
varImp$Predictor <- rownames(varImp)
rownames(varImp) <- NULL
colnames(varImp) <- c("Importance", "Predictor")
pval <- VI$`Predictor p-values`
pval$Predictor <- rownames(pval)
rownames(pval) <- NULL
colnames(pval) <- c("pval", "Predictor")

varImp <- left_join(varImp, pval)

# Clean predictor names (optional)
varImp$Predictor <- factor(varImp$Predictor,
                           levels = c("Geographic", "ele", "ppt", "fire1_5yr_cbi_mn", 
                                      "fire6_10yr_cbi_mn", "fire11_35yr_cbi_mn", "cancov", "ch_res"),
                           labels = c("Geographic", "Elevation", "Precipitation", "Fire 1-5 yr", 
                                      "Fire 6-10 yr", "Fire 11-35 yr", "Canopy Cover", "Canopy Height"))

# Add significance stars
varImp$sig <- ifelse(varImp$pval < 0.05, "*", "")

# Create lollipop chart with significance stars
ggplot(varImp, aes(x = reorder(Predictor, Importance), y = Importance)) +
  geom_segment(aes(x = reorder(Predictor, Importance), 
                   xend = reorder(Predictor, Importance),
                   y = 0, 
                   yend = Importance),
               color = "grey50") +
  geom_point(size = 3, color = "steelblue") +
  geom_text(aes(label = sig), 
            vjust = -0.5,
            size = 5) +  # Adjust size as needed
  coord_flip() +
  theme_minimal() +
  labs(x = "",
       y = "Relative Importance (%)",
       title = "Variable Importance") +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14)
  )


# Read and process all files
gdm.fits <- list.files(path = "R:/Users/skeyser/Postdoc/GDM_Output/", full.names = TRUE)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Variable Importance
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Function to extract variable importance from a single model
extract_var_imp <- function(model) {
  
  map_dfr(1:nrow(model), function(i) {
    
    ## Variable Importance
    VI <- model$varImp[[i]]
    varImp <- VI$`Predictor Importance`
    varImp$Predictor <- rownames(varImp)
    rownames(varImp) <- NULL
    colnames(varImp) <- c("Importance", "Predictor")
    
    ## P-value for predictors
    pval <- VI$`Predictor p-values`
    pval$Predictor <- rownames(pval)
    rownames(pval) <- NULL
    colnames(pval) <- c("pval", "Predictor")
    
    result <- left_join(varImp, pval)
    result$Model <- model$beta_metric[i]
    return(result)
  })
}

# Process all models
all_results <- map_dfr(gdm.fits, function(x) {
  model <- readRDS(x)
  extract_var_imp(model)
}, .id = "model_id")

# Calculate summary statistics
summary_stats <- all_results %>%
  group_by(Predictor, Model) %>%
  summarise(
    mean_importance = mean(Importance),
    lower_ci = quantile(Importance, 0.025),
    upper_ci = quantile(Importance, 0.975),
    sig_prop = mean(pval < 0.05)  # proportion of significant results
  )

# Clean predictor names
summary_stats$Predictor <- factor(summary_stats$Predictor,
                                  levels = c("Geographic", "ele", "ppt", "fire1_5yr_cbi_mn", 
                                             "fire6_10yr_cbi_mn", "fire11_35yr_cbi_mn", "cancov", "ch_res"))

# Add significance indicator
summary_stats$sig <- ifelse(summary_stats$sig_prop > 0.95, "*", "") 

summary_stats <- summary_stats |> mutate(PredPretty = case_when(Predictor == "ele" ~ "Elevation",
                                                                Predictor == "cancov" ~ "Canopy Cover",
                                                                Predictor == "ch_res" ~ "Canopy Height Res.",
                                                                Predictor == "Geographic" ~ "Geog. Dist.",
                                                                Predictor == "ppt" ~ "Precipitation",
                                                                Predictor == "fire1_5yr_cbi_mn" ~ "Fire Severity: 1-5yr",
                                                                Predictor == "fire6_10yr_cbi_mn" ~ "Fire Severity: 6-10yr",
                                                                Predictor == "fire11_35yr_cbi_mn" ~ "Fire Severity: 11-35yr"
)) |> 
  mutate(Color = case_when(PredPretty %in% c("Elevation", 
                                             "Geog. Dist.",
                                             "Precipitation") ~ "#7570b3",
                           PredPretty %in% c("Canopy Cover",
                                             "Canopy Height Res."
                           ) ~ "#1b9e77",
                           PredPretty %in% c(
                             "Fire Severity: 1-5yr",
                             "Fire Severity: 6-10yr",
                             "Fire Severity: 11-35yr"
                           ) ~ "#d95f02")) |> 
  mutate(Model = factor(case_when(
    Model == "Brep" ~ "Turnover",
    Model == "Brich" ~ "Nestedness",
    Model == "Btotal" ~ "Total Beta"
  ), levels = c("Total Beta", "Turnover", "Nestedness")))

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Plotting the lollipop plot for VarImp
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ## Lollipop plot for the VI
# vi_lolli <- ggplot(summary_stats, aes(x = reorder(PredPretty, mean_importance), 
#                           y = mean_importance)) +
#   # Draw segments
#   geom_segment(aes(x = reorder(PredPretty, mean_importance),
#                    xend = reorder(PredPretty, mean_importance),
#                    y = 0, 
#                    yend = mean_importance),
#                color = "gray50",
#                size = 1) +
#   # Draw white background for points to "break" the line
#   geom_point(size = 5, 
#              color = "white",
#              fill = "white") +
#   # Draw actual points
#   geom_point(aes(color = Color,
#                  shape = Model),
#              size = 3) +
#   geom_text(aes(label = sig), 
#             vjust = -0.5,
#             size = 5) +
#   scale_color_identity() +
#   scale_shape_manual(values = c("Total Beta" = 16,  # Filled circle
#                                 "Turnover" = 18,      # Open circle
#                                 "Nestedness" = 17),   # Open triangle
#                      name = "Component") +
#   coord_flip() +
#   scale_y_continuous(breaks = seq(0, 50, by = 10),
#                      limits = c(0, 50)) +
#   theme_bw() +
#   labs(x = "",
#        y = "Predictor Importance") +
#   theme(axis.text = element_text(size = 12),
#         axis.title = element_text(size = 14),
#         legend.position = "right",
#         panel.spacing = unit(0.75, "lines"),
#         legend.background = element_blank())
# 
# 
# vi_lolli <- ggplot(summary_stats, aes(x = reorder(PredPretty, mean_importance), 
#                                       y = mean_importance)) +
#   # Draw segments
#   geom_segment(aes(x = reorder(PredPretty, mean_importance),
#                    xend = reorder(PredPretty, mean_importance),
#                    y = 0, 
#                    yend = mean_importance),
#                color = "gray50",
#                size = 1) +
#   # Draw white background for points to "break" the line
#   geom_point(size = 5, 
#              color = "white",
#              fill = "white") +
#   # Draw actual points
#   geom_point(aes(color = Color,
#                  shape = Model),
#              size = 3) +
#   # Add dummy geom_line for legend
#   geom_line(aes(linetype = Model),
#             color = "black",
#             show.legend = TRUE) +
#   geom_text(aes(label = sig), 
#             vjust = -0.5,
#             size = 5) +
#   scale_color_identity() +
#   scale_shape_manual(values = c("Total Beta" = 16,
#                                 "Turnover" = 18,
#                                 "Nestedness" = 17),
#                      name = "Component") +
#   scale_linetype_manual(values = c("Total Beta" = "solid",
#                                    "Turnover" = "dashed",
#                                    "Nestedness" = "dotted"),
#                         name = "Component") +
#   guides(shape = guide_legend(override.aes = list(color = "black"))) +
#   coord_flip() +
#   scale_y_continuous(breaks = seq(0, 50, by = 10),
#                      limits = c(0, 50)) +
#   theme_bw() +
#   labs(x = "",
#        y = "Predictor Importance") +
#   theme(axis.text = element_text(size = 12),
#         axis.title = element_text(size = 14),
#         legend.position = "right",
#         panel.spacing = unit(0.75, "lines"),
#         legend.background = element_blank())


# Create segment data with preserved ordering
# Set the desired order
plot_order <- rev(c( 
  "Geog. Dist.",
  "Elevation",
  "Precipitation",
  "Canopy Cover",
  "Canopy Height Res.",
  "Fire Severity: 1-5yr",
  "Fire Severity: 6-10yr",
  "Fire Severity: 11-35yr"))

vi_lolli2 <- summary_stats |> 
  mutate(PredPretty = factor(PredPretty, levels = plot_order)) |>  
  ggplot() +
  # Draw segments using separate data frame
  geom_segment(aes(x = PredPretty,  
                   y = 0, 
                   yend = mean_importance,
                   group = Model),
               color = "gray50",
               size = 0.7,
               position = position_dodge(width = 0.5)) +
  # Draw points using original data
  geom_point(aes(x = PredPretty,
                 y = mean_importance,
                 color = Color,
                 shape = Model),
             size = 3,
             position = position_dodge(width = 0.5)) +
  # Add dummy geom_line for legend
  geom_line(aes(x = PredPretty,
                y = mean_importance,
                linetype = Model),
            color = "black",
            show.legend = TRUE) +
  geom_hline(yintercept = 0,
             linetype = "dashed",
             color = "black",
             size = 0.5) +
  geom_text(aes(x = PredPretty,
                y = mean_importance,
                label = sig,
                group = Model), 
            vjust = 0.7,
            hjust = -1.5,
            size = 5,
            position = position_dodge(width = 0.5)) +
  scale_color_identity() +
  scale_shape_manual(values = c("Total Beta" = 16,
                                "Turnover" = 18,
                                "Nestedness" = 17),
                     name = "Component") +
  scale_linetype_manual(values = c("Total Beta" = "solid",
                                   "Turnover" = "dashed",
                                   "Nestedness" = "dotted"),
                        name = "Component") +
  guides(shape = guide_legend(override.aes = list(color = "black"))) +
  coord_flip() +
  scale_y_continuous(breaks = seq(0, 50, by = 10),
                     limits = c(0, 50)) +
  theme_bw() +
  labs(x = "",
       y = "Predictor Importance") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.position = "right",
        panel.spacing = unit(0.75, "lines"),
        legend.background = element_blank())

vi_lolli2

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: VI Lollipop using total spline and perm significance
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
total_spline_df <- function(model.paths) {
  spline_df <- map_dfr(seq_along(model.paths), function(i){ 
    message(sprintf("Processing model %d of %d", i, length(model.paths)))
    path <- model.paths[i]
    model <- readRDS(path)
    
    model |> 
      rowwise() |> 
      mutate(
        spline_data = list({
          # Get number of predictors and coefficients
          n_pred <- length(models$predictors)
          coef_mat <- matrix(models$coefficients, nrow = n_pred, ncol = 3, byrow = TRUE)
          
          # Create data frame with predictor names and spline sums
          data.frame(
            Pred = models$predictors,
            spline_sum = rowSums(abs(coef_mat))
          )
        })
      ) |>
      unnest(spline_data)
  }, .id = "model_id")
  
  # Calculate summary statistics by beta_metric and predictor
  summary_df <- spline_df %>%
    group_by(beta_metric, Pred) %>%
    summarise(
      mean_spline = mean(spline_sum),
      lower_ci = quantile(spline_sum, 0.025),
      upper_ci = quantile(spline_sum, 0.975),
      .groups = 'drop'
    )
  
  # Add pretty names
  cov_rename <- data.frame(
    Pred = c("ele", "cancov", "Geographic", "ch_res", "ppt",
             "fire1_5yr_cbi_mn", "fire6_10yr_cbi_mn", "fire11_35yr_cbi_mn"), 
    PredPretty = c("Elevation", "Canopy Cover", "Geog. Dist.", 
                   "Canopy Height Res.", "Precipitation",
                   "Fire Severity: 1-5yr", "Fire Severity: 6-10yr", 
                   "Fire Severity: 11-35yr")
  )
  
  summary_df <- summary_df %>%
    left_join(cov_rename, by = "Pred")
  
  return(summary_df)
}

total_spline <- total_spline_df(model.paths = gdm.fits)
glimpse(total_spline)
gc()


# Create segment data with preserved ordering
# Set the desired order
plot_order <- rev(c( 
  "Geog. Dist.",
  "Elevation",
  "Precipitation",
  "Canopy Cover",
  "Canopy Height Res.",
  "Fire Severity: 1-5yr",
  "Fire Severity: 6-10yr",
  "Fire Severity: 11-35yr"))

spline_lolli <- total_spline |>
  mutate(Model = factor(case_when(
    beta_metric == "Brep" ~ "Turnover",
    beta_metric == "Brich" ~ "Nestedness",
    beta_metric == "Btotal" ~ "Total Beta"
  ), levels = c("Total Beta", "Turnover", "Nestedness"))) |> 
  left_join(summary_stats |> dplyr::select(PredPretty, Model, sig_prop, sig, Color)) |> 
  mutate(PredPretty = factor(PredPretty, levels = plot_order)) |>  
  ggplot() +
  # Draw segments using separate data frame
  geom_segment(aes(x = PredPretty,  
                   y = 0, 
                   yend = mean_spline,
                   group = Model),
               color = "gray50",
               size = 0.7,
               position = position_dodge(width = 0.5)) +
  # Draw points using original data
  geom_point(aes(x = PredPretty,
                 y = mean_spline,
                 color = Color,
                 shape = Model),
             size = 3,
             position = position_dodge(width = 0.5)) +
  # Add dummy geom_line for legend
  geom_line(aes(x = PredPretty,
                y = mean_spline,
                linetype = Model),
            color = "black",
            show.legend = TRUE) +
  geom_hline(yintercept = 0,
             linetype = "dashed",
             color = "black",
             size = 0.5) +
  geom_text(aes(x = PredPretty,
                y = mean_spline,
                label = sig,
                group = Model), 
            vjust = 0.7,
            hjust = -1.5,
            size = 5,
            position = position_dodge(width = 0.5)) +
  scale_color_identity() +
  scale_shape_manual(values = c("Total Beta" = 16,
                                "Turnover" = 18,
                                "Nestedness" = 17),
                     name = "Component") +
  scale_linetype_manual(values = c("Total Beta" = "solid",
                                   "Turnover" = "dashed",
                                   "Nestedness" = "dotted"),
                        name = "Component") +
  guides(shape = guide_legend(override.aes = list(color = "black"))) +
  coord_flip() +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1),
                     limits = c(0, 1)) +
  theme_bw() +
  labs(x = "",
       y = "Effect Size") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.position = "right",
        panel.spacing = unit(0.75, "lines"),
        legend.background = element_blank())

spline_lolli


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Spline functions
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Summary stats on the fits
dev_exp <- gdm.fit %>%
  mutate(
    # First unlist and extract the explained deviance
    deviance_explained = map(models, function(model_list) {
      sapply(model_list, function(x) x$explained)
    })
  ) %>%
  group_by(beta_metric) %>%
  summarise(
    n = n(),
    mean_deviance = mean(unlist(deviance_explained)),
    lower_ci = quantile(unlist(deviance_explained), 0.025),
    upper_ci = quantile(unlist(deviance_explained), 0.975),
    sd_deviance = sd(unlist(deviance_explained))
  )
summary_stats

## Get the effects of each predictor
## Function to extract coefficients from one model
make_spline_df <- function(model.paths, var_name) {
  spline_df <- map_dfr(seq_along(model.paths), function(i){ 
    message(sprintf("Processing model %d of %d", i, length(model.paths)))
    path <- model.paths[i]
    model <- readRDS(path)
    model |> 
      rowwise() |> 
      mutate(
        spline_data = list({
          spline = isplineExtract(models)
          data.frame(
            x = spline$x[,var_name],
            y = spline$y[,var_name]
          )
        })
      ) |> 
      unnest(spline_data) |>
      select(beta_metric, x, y)
  }, .id = "model_id")
  ## Find the covariate name
  cov_rename <- data.frame(Pred = c("ele", 
                                    "cancov", 
                                    "Geographic", 
                                    "ch_res", 
                                    "ppt",
                                    "fire1_5yr_cbi_mn",
                                    "fire6_10yr_cbi_mn",
                                    "fire11_35yr_cbi_mn"
  ), 
  PredPretty = c("Elevation", 
                 "Canopy Cover",
                 "Geog. Dist.", 
                 "Canopy Height Res.",
                 "Precipitation",
                 "Fire Severity: 1-5yr",
                 "Fire Severity: 6-10yr",
                 "Fire Severity: 11-35yr"
  )) |> 
    mutate(Color = case_when(PredPretty %in% c("Elevation", 
                                               "Geog. Dist.",
                                               "Precipitation"
    ) ~ "#7570b3",
    PredPretty %in% c("Canopy Cover",
                      "Canopy Height Res."
    ) ~ "#1b9e77",
    PredPretty %in% c("Fire Severity: 1-5yr",
                      "Fire Severity: 6-10yr",
                      "Fire Severity: 11-35yr"
    ) ~ "#d95f02"))
  
  ## Return a list with the necessary info to automate plotting
  spline_df <- list(spline_df = spline_df,
                    plot_name = cov_rename$PredPretty[cov_rename$Pred == var_name],
                    color = cov_rename$Color[cov_rename$Pred == var_name])
  
  return(spline_df)
}

## Function for plotting
spline_plot <- function(spline_list){
  
  ## Take the DF for the plot
  spline_df <- spline_list$spline_df
  
  ## Convert m to km for geographic distance
  if(spline_list$plot_name == "Geog. Dist.") {
    spline_df$x <- spline_df$x/1000
    spline_list$plot_name <- "Geog. Dist. (km)"
  }
  
  # Calculate mean spline
  mean_spline <- spline_df |> 
    group_by(beta_metric, x) |> 
    summarize(mean_y = mean(y)) |> 
    mutate(beta_metric = factor(case_when(
      beta_metric == "Brep" ~ "Turnover",
      beta_metric == "Brich" ~ "Nestedness",
      beta_metric == "Btotal" ~ "Total Beta"
    ), levels = c("Total Beta", "Turnover", "Nestedness")))
  
  ## Calculate the 95% CI for splines
  spline_ci <- spline_df |> 
    group_by(beta_metric, x) |> 
    summarize(hi_y = quantile(y, probs = 0.975),
              lo_y = quantile(y, probs = 0.025)) |> 
    mutate(beta_metric = factor(case_when(
      beta_metric == "Brep" ~ "Turnover",
      beta_metric == "Brich" ~ "Nestedness",
      beta_metric == "Btotal" ~ "Total Beta"
    ), levels = c("Total Beta", "Turnover", "Nestedness")))
  
  splot <- ggplot() +
    # Individual splines
    # geom_line(data = spline_df, 
    #           aes(x = x, y = y, group = model),
    #           color = color, 
    #           alpha = 0.5) +
    # Mean spline
    geom_line(data = mean_spline,
              aes(x = x, 
                  y = mean_y,
                  linetype = beta_metric),
              color = spline_list$color,
              size = 1) +
    scale_linetype_manual(values = c("Total Beta" = "solid",
                                     "Turnover" = "dashed",
                                     "Nestedness" = "dotted"),
                          #name = "Component"
                          guide = "none"
    ) +
    scale_fill_manual(values = c("Total Beta" = "gray70",
                                 "Turnover" = "gray70",
                                 "Nestedness" = "gray70"),
                      name = "Component") +
    geom_ribbon(data = spline_ci,
                aes(x = x, 
                    ymin = lo_y,
                    ymax = hi_y,
                    group = beta_metric),
                alpha = 0.3,
                fill = spline_list$color) +
    theme_bw() +
    theme(axis.title = element_text(family = "sans", size = 14),
          axis.text = element_text(family = "sans", size = 12),
          legend.position = "none"
          # legend.direction = "horizontal",
          # legend.box = "horizontal",
          # legend.margin = margin(t = 5, b = 5), # Adjusts top and bottom margins
          # legend.spacing.x = unit(0.5, 'cm')
    ) +
    scale_y_continuous(limits = c(0, 1)) +
    xlab(gsub("Mean ", "", spline_list$plot_name)) +
    ylab("Ecological Distance") +
    labs(linetype = NULL)
  
  ## Return the plot
  return(splot)
}

## Plot each of the variable for patchworking later
## Make the DFs
ele_df <- make_spline_df(gdm.fits, var_name = "ele")
ppt_df <- make_spline_df(gdm.fits, var_name = "ppt")
cancov_df <- make_spline_df(gdm.fits, var_name = "cancov")
canht_df <- make_spline_df(gdm.fits, var_name = "ch_res")
geo_df <- make_spline_df(gdm.fits, var_name = "Geographic")
f5h_df <- make_spline_df(gdm.fits, var_name = "fire1_5yr_cbi_mn")
f10h_df <- make_spline_df(gdm.fits, var_name = "fire6_10yr_cbi_mn")
f35h_df <- make_spline_df(gdm.fits, var_name = "fire11_35yr_cbi_mn")


## Feed into the plotting function
p.ele <- spline_plot(ele_df)
p.ppt <- spline_plot(ppt_df)
p.cc <- spline_plot(cancov_df)
p.ch <- spline_plot(canht_df)
p.geo <- spline_plot(geo_df)
p.f5h <- spline_plot(f5h_df)
p.f10h <- spline_plot(f10h_df)
p.f35h <- spline_plot(f35h_df)


library(patchwork)

pgrid_lyt <- "
AAB
AAC
DEF
HIJ
"
pgrid <- free(spline_lolli) + 
  p.geo + 
  p.ele +
  p.ppt +
  p.cc +
  p.ch +
  p.f5h + 
  p.f10h + 
  p.f35h + 
  plot_layout(design = pgrid_lyt, guides = "collect") &
  theme(legend.position = "bottom",
        legend.key.size = unit(1, "cm"),  # Make legend symbols bigger
        legend.key.width = unit(2, "cm"),  # Make legend lines longer
        legend.text = element_text(size = 12),  # Adjust text size if needed
        plot.margin = margin(t = 5, r = 15, b = 5, l = 5))

ggsave(filename = here("./Figures/Final/GridLolliSpline_GDM_AllSamples_BetaComponents_SpVarThresh_Fig3.jpg"),
       plot = pgrid, height = 12, width = 12,
       dpi = 800)
ggsave(filename = here("./Figures/Final/Fig3.jpg"),
       plot = pgrid, height = 12, width = 12,
       dpi = 800)