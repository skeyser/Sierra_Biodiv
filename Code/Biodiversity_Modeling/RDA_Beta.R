## -------------------------------------------------------------
##
## Script name: Multivariate Biodiversity Analyses
##
## Script purpose:
##
## Author: Spencer R Keyser
##
## Date Created: 2025-02-20
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
library(purrr)
library(ggplot2)
library(ggrepel)
library(here)
library(vegan)
library(adespatial)
library(ade4)
library(adegraphics)
library(sf)
library(spdep)
library(CAbioacoustics)
## -------------------------------------------------------------

## Read in corrected species composition
#load(here("./Data/JAGS_Data/Occ2GDM_Data_95thresh.Rdata"))
load(here("./Data/JAGS_Data/Occ2GDM_Data_SpThresh_975minMaxPrec.Rdata"))

dat <- OccData
rm(OccData)

## Take a look
str(dat)

## Grab the posterior
zmat <- dat$Z.posterior
zmat <- zmat[,,1]
nrow(zmat)

## Test influence of preds using RDA
ind.miss <- which(rowSums(zmat) == 0)
zmat <- zmat[-ind.miss,]
nrow(zmat)
colnames(zmat) <- dat$ZcolNames

## Add environmental data
env <- dat$SiteMeta

## Env test
env.test <- env |> 
  slice(-ind.miss) |> 
  select(Cell_Unit, Y, X,
         ele, ppt, cancov, ch_res,
         fire1_5yr_cbi_mn,
         fire6_10yr_cbi_mn,
         fire11_35yr_cbi_mn) |>
  mutate(across(where(is.numeric), 
                list(scaled = ~scale(., center = TRUE, scale = TRUE)))) |> 
  as.data.frame()

## Bring in the pyrodiversity metric
pyro <- read.csv(here("./Data/FireSeverity2021_MeanStDev_AllARUs.csv"))
interval_naming <- c(
  "2016_2020" = "1_5",
  "2011_2015" = "6_10",
  "1986_2010" = "11_35"
)
pyro <- pyro |> 
  select(deployment_name, contains("stdev")) |> 
  rename_with(
    ~str_replace_all(., interval_naming),
    contains("Fire_Sev")
  ) |> 
  mutate(Cell_Unit = gsub("G\\d{3}_V\\d{1}_", "", deployment_name)) |> 
  select(Cell_Unit, contains("Fire_Sev")) #|> 
  filter(Cell_Unit %in% env.test$Cell_Unit) |> 
  distinct(Cell_Unit, .keep_all = T)

## Add the pyrodiversity metric to the envTest
env.test <- left_join(env.test, pyro)


## Subset the data down to the points that have
## fire at all
env.fire <- env |>
  slice(-ind.miss) |> 
  filter(if_any(contains("fire"), ~. > 0)) |> 
  select(Cell_Unit, utmn, utme,
         ele:ch_res) |>
  mutate(across(where(is.numeric), 
                list(scaled = ~scale(., center = TRUE, scale = TRUE)))) |> 
  as.data.frame() 


fire.rows <- which(env.test$Cell_Unit %in% env.fire$Cell_Unit)

zmat.fire <- zmat[fire.rows,]
zmat <- zmat.fire

zmat_com <- adespatial::beta.div.comp(zmat, "BS", quant = FALSE)
zmat_com$part
zmat_nest <- zmat_com$rich
zmat_turn <- zmat_com$repl
zmat_tot <- zmat_com$D

ad4::is.euclidean

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Create dbMEMs for Coordinates
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Load the data with the coords and UTM zones
aru_meta <- read.csv(here("./Data/ARU_120m_New.csv"))
aru_meta <- aru_meta |>
  rename(Cell_Unit = cell_unit) |> 
  select(Cell_Unit, Year = survey_yea, X, Y) |> 
  filter(Cell_Unit %in% env.test$Cell_Unit) |>  
  arrange(match(Cell_Unit, env.test$Cell_Unit))

## Take the xy data
coords <- cbind(aru_meta$X, aru_meta$Y)
hist(coords[,1], breaks = 50)
hist(coords[,2], breaks = 50)

## Make a neighborhood list w. Gabriel neighborhoods
nbgab <- graph2nb(gabrielneigh(coords), sym = T)
gab.p <- s.label(coords, nb = nbgab, pnb.edge.col = "red", main = "Gabriel", plot = T)

## Weight by distance
distgab <- nbdists(nbgab, coords)
fdist <- lapply(distgab, function(x) 1 - x/max(dist(coords)))

## Create spatial weighting matrix
listwgab <- nb2listw(nbgab, glist = fdist)

## Make the MEMs
mem.sierra <- mem(listwgab)
mem.sierra

## Plot MEMs
barplot(attr(mem.sierra, "values"), 
        main = "Eigenvalues of the spatial weighting matrix", cex.main = 0.7)

plot(mem.sierra[,c(1, 5, 10, 20, 30, 40, 50, 60, 70)], SpORcoords = coords, symbol = "circle")

## Calcualte Moran'sI
moranI <- moran.randtest(mem.sierra, listwgab, 99)

## Examine the spatial patterns in the predictors
MC.env <- moran.randtest(env.test[,4:10], listwgab, nrepet = 999)

mc.bounds <- moran.bounds(listwgab)
mc.bounds

env.maps <- s1d.barchart(MC.env$obs, labels = MC.env$names, plot = FALSE, xlim = 1.1 * mc.bounds, paxes.draw = TRUE, pgrid.draw = FALSE)
addline(env.maps, v = mc.bounds, plot = TRUE, pline.col = 'red', pline.lty = 3)

## Testing for spatial autocorrelation in species data
pca.bird <- dudi.pca(zmat, scale = FALSE, scannf = F, nf = 2)

moran.randtest(pca.bird$li, listw = listwgab)
s.value(coords, pca.bird$li, Sp = sierra.hull, symbol = "circle", col = c("white", "palegreen4"), ppoint.cex = 0.6)

ms.bird <- multispati(pca.bird, listw = listwgab, scannf = F)
summary(ms.bird)

g.ms.maps <- s.value(coords, ms.bird$li, 
                     symbol = "circle", 
                     col = c("white", "palegreen4"), 
                     ppoint.cex = 0.6,
                     xlim = c(min(coords[,1]), max(coords[,1])),
                     ylim = c(min(coords[,2]), max(coords[,2])))


## Species-specific spatial patterning
g.ms.spe <- s.arrow(ms.bird$c1, plot = FALSE)
g.abund <- s.value(coords, zmat[, c(1,50,20,10)],
                   xlim = c(min(coords[,1]), max(coords[,1])),
                   ylim = c(min(coords[,2]), max(coords[,2])), 
                   symbol = "circle", 
                   col = c("black", "palegreen4"), 
                   plegend.drawKey = FALSE, 
                   ppoint.cex = 0.4, 
                   plot = FALSE)
p1 <- list(c(0.05, 0.65), c(0.01, 0.25), c(0.74, 0.58), c(0.55, 0.05))
for (i in 1:4)
  g.ms.spe <- insert(g.abund[[i]], g.ms.spe, posi = p1[[i]], ratio = 0.25, plot = FALSE)
g.ms.spe

## Scalogram
scalo <- scalogram(zmat[,1], mem.sierra, nblocks = 50)
plot(scalo)

## Choosing the best MEMs for the spatial component of the analysis
mem.sierra.sel <- mem.select(bird.pca$tab, listw = listwgab)
dim(mem.sierra.sel$MEM.select)
## Take some 115 MEMs

rda.sierra <- pcaiv(bird.pca, mem.sierra.sel$MEM.select, scannf = FALSE)
test.rda <- randtest(rda.sierra)
test.rda
plot(test.rda)

s.value(coords, rda.sierra$li, 
        symbol = "circle", 
        col = c("white", "palegreen4"), 
        ppoint.cex = 0.6,
        xlim = c(min(coords[,1]), max(coords[,1])),
        ylim = c(min(coords[,2]), max(coords[,2])))

## Variance partitioning
vp1 <- varpart(bird.pca$tab, env.test.ade[,-1], mem.sierra.sel$MEM.select)
vp1
plot(vp1, bg = c(3, 5), Xnames = c("environment", "spatial"))

## Significance of fractions
# Test fraction [a] - Pure environmental effects
mem.df <- as.data.frame(mem.sierra.sel$MEM.select)
env.mem <- cbind(env.test.ade[,-1], mem.df)
rda.a <- rda(bird.pca$tab, env.mem[,1:7])
anova(rda.a)

# Test fraction [b] - Pure spatial effects
rda.b <- rda(bird.pca$tab, env.mem[,8:ncol(env.mem)])
anova(rda.b)

# Test full model [a+b+c]
rda.full <- rda(bird.pca$tab ~ ., 
                data = env.mem)
anova(rda.full)

# Note: Fraction [c] (shared variation) cannot be tested directly

## Model selection criteria for the variables of interest
# Stepwise selection
step.sel <- ordistep(rda(bird.pca$tab ~ 1, data=env.mem),  # null model
                     scope = formula(rda(bird.pca$tab ~ ., data=env.mem)),  # full model
                     direction = "both",
                     permutations = 999)

# View results
step.sel$anova

## -------------------------------------------------------------
##
## End Section: Moran's Eigen Maps
##
## -------------------------------------------------------------


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Fit dbRDAs
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#env.test <- env.test[-ind.miss,]
mout1 <- dbrda(zmat_nest ~ Y_scaled + ele_scaled + ppt_scaled + cancov_scaled + ch_res_scaled + fire1_5yr_cbi_mn_scaled + fire6_10yr_cbi_mn_scaled + fire11_35yr_cbi_mn_scaled, 
               data = env.test)
summary(mout1)
screeplot(mout1)

mout2 <- dbrda(zmat_turn ~ Y_scaled + ele_scaled + ppt_scaled + cancov_scaled + ch_res_scaled + fire1_5yr_cbi_mn_scaled + fire6_10yr_cbi_mn_scaled + fire11_35yr_cbi_mn_scaled, 
               data = env.test)
summary(mout2)
screeplot(mout2)

mout3 <- dbrda(zmat_tot ~ ppt_scaled + cancov_scaled + ch_res_scaled + fire1_5yr_cbi_mn_scaled + fire6_10yr_cbi_mn_scaled + fire11_35yr_cbi_mn_scaled + Condition(Y_scaled + ele_scaled), 
               data = env.test)
summary(mout3)
screeplot(mout3)

R2adjNest <- RsquareAdj(mout1)
R2adjRepl <- RsquareAdj(mout2)
R2adjBtot <- RsquareAdj(mout3)
vif.cca(mout1)

## Fire points only
env.fire.test <- env.fire |> 
  select(contains("scaled"))

mout1 <- dbrda(zmat_nest ~ . + Condition(MEM1, MEM2), data = env.fire.test)
summary(mout1)
mout2 <- dbrda(zmat_turn ~ ele_scaled + stage_scaled + cancov_scaled + fire1_5yr_high_prop_scaled + fire6_10yr_high_prop_scaled 
               + fire11_35yr_high_prop_scaled + fire1_5yr_lowmod_prop_scaled + fire6_10yr_lowmod_prop_scaled 
               + fire11_35yr_lowmod_prop_scaled + Condition(MEM1, MEM2), data = env.fire.test)

vp <- vegan::varpart(zmat_turn, env.test[, c(2,4)], env.test[,5:6], env.test[, 7:12])
plot(vp)
summary(mout2)

vp <- vegan::varpart(zmat_tot, env.fire.test[, 1:3], env.fire.test[,4:9], env.fire.test[, 10:11])
plot(vp)
mout3 <- dbrda(zmat_tot ~ utmn_scaled + utme_scaled + ele_scaled + stage_scaled + cancov_scaled + fire1_5yr_high_prop_scaled + fire6_10yr_high_prop_scaled 
               + fire11_35yr_high_prop_scaled + fire1_5yr_lowmod_prop_scaled + fire6_10yr_lowmod_prop_scaled 
               + fire11_35yr_lowmod_prop_scaled + Condition(sign_MEMs$MEM1, sign_MEMs$MEM2), data = env.fire.test)

summary(mout3)



mout.aov1 <- anova.cca(mout1, by = "terms")
mout.aov1
mout.aov2 <- anova.cca(mout2, by = "terms")
mout.aov2
## Total beta diversity
## By terms
mout.aov3 <- anova.cca(mout3, by = "terms")
mout.aov3

## Overall
mout.aov1global <- anova(mout1)
mout.aov2global <- anova(mout2)
mout.aov3global <- anova(mout3)

mout <- dbrda(zmat ~ ., data = env.test, dist = "raup")
summary(mout)

mout.aov <- anova(mout, by = "terms")

R2_dbrda <- RsquareAdj(mout3)$r.squared

# Forward selection of variables:
back.sel <- ordiR2step(rda(zmat ~ 1, data = env.test), # lower model limit (simple!)
                      scope = formula(mout3), # upper model limit (the "full" model)
                      direction = "backward",
                      R2scope = TRUE, # can't surpass the "full" model's R2
                      pstep = 1000,
                      trace = FALSE) # change to TRUE to see the selection process!


plot_dbrda_density <- function(dbrda_result) {
  # Extract site scores
  site_scores <- as.data.frame(scores(dbrda_result, display = "sites"))
  env_scores <- as.data.frame(scores(dbrda_result, display = "bp"))
  
  rownames(env_scores) <- c("Latitude", 
                            "Elevation", 
                            "Precipitation",
                            "Canopy Cover",
                            "Canopy Height",
                            "Fire Severity 1-5yr",
                            "Fire Severity 6-10yr",
                            "Fire Severity 11-35yr")
  
  # Create density plot
  p <- ggplot() +
    geom_hex(data = site_scores,
             aes(x = dbRDA1, y = dbRDA2),
             bins = 20,
             color = "white",
             alpha = 0.8) +
    # Add density contours
    # stat_density_2d(data = site_scores,
    #                 aes(x = dbRDA1, y = dbRDA2, fill = ..level..),
    #                 geom = "polygon") +
    # Add environmental vectors
    geom_segment(data = env_scores,
                 aes(x = 0, y = 0, xend = dbRDA1, yend = dbRDA2),
                 arrow = arrow(length = unit(0.1, "cm"),
                               type = "closed"),
                 linewidth = 1,
                 color = "#8856a7") +
    # # Add variable labels
    # geom_text(data = env_scores,
    #           aes(x = dbRDA1, y = dbRDA2, label = rownames(env_scores)),
    #           hjust = -0.2) +
    geom_text_repel(data = env_scores,
                    aes(x = dbRDA1, y = dbRDA2, label = rownames(env_scores)),
                    box.padding = 0.8,
                    point.padding = 0.8,
                    force = 1,
                    force_pull = 0,
                    fontface = "bold",
                    segment.color = "black",
                    max.overlaps = Inf) + 
    #scale_fill_viridis_c(name = "Count", option = "C") +
    scale_fill_gradientn(name = "Count",
                         colors = RColorBrewer::brewer.pal(9, "PuBuGn")[2:9]) +
    theme_bw() +
    coord_equal()
  
  return(p)
}

dbrdaHex <- plot_dbrda_density(dbrda_result = mout1)
dbrdaHex <- plot_dbrda_density(dbrda_result = mout2)
dbrdaHex <- plot_dbrda_density(dbrda_result = mout3)

ggsave(plot = dbrdaHex, filename = here("./Figures/Preliminary/dbRDA_Hex_Biplot_RC_SpVar_Rdiff.jpg"),
       height = 6, width = 6, dpi = 600)

# Check the new model with forward-selected variables
fwd.sel$call

## Test for significant drivers
mout.sig <- anova.cca(mout2, step = 1000, by = "term")

# Function for dbRDA sensitivity analysis with distance matrix
dbrda_sensitivity <- function(dist_matrix, env_data, 
                              size_seq=NULL, n_iter=100,
                              sqrt_dist=FALSE, add=NULL) {
  
  require(vegan)
  
  # If size sequence not specified, create default
  if(is.null(size_seq)) {
    n <- nrow(as.matrix(dist_matrix))
    size_seq <- round(seq(n/4, n, length.out=10))
  }
  
  # Store results
  results <- data.frame(
    sample_size = numeric(),
    iteration = numeric(),
    R2 = numeric(),
    adj_R2 = numeric(),
    F_stat = numeric(),
    p_value = numeric()
  )
  
  # Run sensitivity analysis
  for(size in size_seq) {
    for(i in 1:n_iter) {
      # Random subsample
      samp <- sample(1:nrow(as.matrix(dist_matrix)), size=size)
      
      # Subset distance matrix and environmental data
      dist_sub <- as.dist(as.matrix(dist_matrix)[samp, samp])
      env_sub <- env_data[samp,]
      env_sub <- env_sub |> select(contains("scaled"))
      
      # Run dbRDA
      dbrda_result <- dbrda(dist_sub ~ ., data=env_sub, 
                            sqrt.dist=sqrt_dist, add=add)
      
      # Extract results
      results <- rbind(results, data.frame(
        sample_size = size,
        iteration = i,
        R2 = RsquareAdj(dbrda_result)$r.squared,
        adj_R2 = RsquareAdj(dbrda_result)$adj.r.squared,
        F_stat = anova(dbrda_result)$F[1],
        p_value = anova(dbrda_result)$Pr[1]
      ))
    }
  }
  
  return(results)
}

sensRDA <- dbrda_sensitivity(dist_matrix = zmat_turn, 
                             env_data = env.test, 
                             size_seq = nrow(zmat_turn)*seq(0.25,1,0.25),
                             n_iter = 5)

# Plot function remains the same
plot_sensitivity <- function(sensitivity_results) {
  require(ggplot2)
  
  # Calculate means and CI for each sample size
  summary_stats <- aggregate(
    cbind(R2, adj_R2) ~ sample_size, 
    data=sensitivity_results,
    FUN=function(x) c(mean=mean(x), 
                      ci_lower=quantile(x, 0.025),
                      ci_upper=quantile(x, 0.975))
  )
  
  # Reshape for plotting
  summary_stats <- data.frame(
    sample_size = summary_stats$sample_size,
    R2_mean = summary_stats$R2[,1],
    R2_lower = summary_stats$R2[,2],
    R2_upper = summary_stats$R2[,3],
    adjR2_mean = summary_stats$adj_R2[,1],
    adjR2_lower = summary_stats$adj_R2[,2],
    adjR2_upper = summary_stats$adj_R2[,3]
  )
  
  # Create plot
  p <- ggplot(summary_stats, aes(x=sample_size)) +
    geom_ribbon(aes(ymin=R2_lower, ymax=R2_upper), alpha=0.2) +
    geom_line(aes(y=R2_mean, color="R2")) +
    geom_ribbon(aes(ymin=adjR2_lower, ymax=adjR2_upper), alpha=0.2) +
    geom_line(aes(y=adjR2_mean, color="Adjusted R2")) +
    labs(x="Sample Size", y="R-squared", color="Measure") +
    theme_bw()
  
  return(p)
}

plot_sensitivity(sensRDA)

# Example usage:
# sens_results <- dbrda_sensitivity(your_dist_matrix, env_data, 
#                                 sqrt_dist=TRUE, add="lingoes")
# plot_sensitivity(sens_results)

## -------------------------------------------------------------
##
## Begin Section: Ade4 Approach
##
## -------------------------------------------------------------

library(ade4)

bird.pca <- dudi.pca(zmat, scannf = F)
env.test.ade <- env.test |> select(contains("scaled")) |> select(-X_scaled)
rda.sierra <- pcaiv(bird.pca, env.test.ade, scannf = F, nf = 2)
plot(rda.sierra)
summary(rda.sierra)
randtest(rda.sierra)

## ratio
sum(rda.sierra$ls[, 1]^2 * rda.sierra$lw) / bird.pca$eig[1]
