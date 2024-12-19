## -------------------------------------------------------------
##
## Script name: MSOM Data Prep
##
## Script purpose: Prepare data for Bayesian MSOM model fit.
##
## Author: Spencer R Keyser
##
## Date Created: 2024-11-27
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
library(abind)

## -------------------------------------------------------------

## -------------------------------------------------------------
##
## Begin Section: Species Data
##
## -------------------------------------------------------------

## Load in the bird data for 2021
load(here("./Data/Generated_DFs/Occ_Mod_Data/2021_OccSppList.RData"))

## Cell Unit mapping file
cu.map <- data.frame(ID = 1:length(sp.det.list[[1]]$Cell_Unit), Cell_Unit = sp.det.list[[1]]$Cell_Unit)
cu.map <- cu.map |>
  mutate(Cell_Unit = ifelse(
    stringr::str_detect(string = Cell_Unit, pattern = "C[0-9]{3}"),
    gsub(pattern = "(C)([0-9]{3})(_U[0-9]+)$", replacement = "\\10\\2\\3", x = Cell_Unit),
    Cell_Unit
  ))

## Reduce the number of dates for sampling
## Check the duration of the sampling period
## 60 sampling periods
length(colnames(sp.det.list[[1]])) - 1

## Array for the data
## D1 (i) = Site, D2 (j) = Sampling Date, D3 (k) = species
samp.cols <- colnames(sp.det.list[[1]])[str_detect(colnames(sp.det.list[[1]]), "\\d")]

## Create different sampling periods
second_samp <- function(DAT, interval, id_col, eff = F, e.var = "Days"){
  tmp.cols <- colnames(DAT)[str_detect(colnames(DAT), "\\d")]
  tmp.splits <- split(tmp.cols, ceiling(seq_along(tmp.cols) / interval))
  tmp.new <- as.data.frame(cbind(DAT[, id_col], matrix(data = NA, nrow = nrow(DAT), ncol = length(tmp.splits), dimnames = list(NULL, paste0("J", seq(1:length(tmp.splits)))))))
  for(i in 1:length(tmp.splits)){
    if(!eff){
      Jsum <- rowSums(DAT[, tmp.splits[[i]]], na.rm = T)
      Jsum <- ifelse(Jsum > 0, 1, 0)
    } 
    if(eff & e.var == "Days"){
      Jsum <- ifelse(DAT[, tmp.splits[[i]]] > 0, 1, 0)
      Jsum <- rowSums(Jsum, na.rm = T)
    }
    if(eff & e.var == "Hrs"){
      Jsum <- rowSums(DAT[, tmp.splits[[i]]], na.rm = T)
    }
    if(eff & e.var == "FirstJDay"){
      mjd <- as.Date(gsub("[[:punct:]]", "-", tmp.splits[[i]]), format = "%Y-%m-%d")
      mjd <- median(lubridate::yday(mjd))
      Jsum <- mjd
    }
    
    tmp.new[,i+1] <- Jsum
  }
  return(tmp.new)
}

## Apply the function
sp.det.list.r <- lapply(sp.det.list, function(x) second_samp(DAT = x, interval = 6, id_col = "Cell_Unit", eff = F))

samp.cols <- colnames(sp.det.list.r[[1]])[str_detect(colnames(sp.det.list.r[[1]]), "J")]
nsite <- nrow(sp.det.list.r[[1]])
nrep <- ncol(sp.det.list.r[[1]]) - 1
nspec <- length(names(sp.det.list.r))
sp.det <- lapply(sp.det.list.r, function(x) x |> select(all_of(samp.cols)))
y <- array(unlist(lapply(sp.det, as.matrix)), dim = c(nsite, nrep, nspec))
dimnames(y) <- list(NULL, NULL, names(sp.det.list))

## No NAs present in the data
table(nsurveys <- apply(y[,,1], 1, function(x) sum(!is.na(x))))

## Species with 0 occurrences
## 9 species
tmp <- apply(y, c(1,3), max, na.rm = TRUE)
tmp[tmp == -Inf] <- NA
sort(obs.occ <- apply(tmp, 2, sum, na.rm = TRUE))

drop.sp <- which(obs.occ == 0)
y <- y[,,-drop.sp]
sp.names <- dimnames(y)

## Redefine nspec
nspec <- dim(y)[3]
# Get observed number of species per site
tmp <- apply(y, c(1,3), max, na.rm = TRUE)
tmp[tmp == "-Inf"] <- NA
sort(C <- apply(tmp, 1, sum)) # Compute and print sorted species counts



## -------------------------------------------------------------
##
## Begin Section: Detection Covariates
##
## -------------------------------------------------------------

## Get the number of hours per survey for the detection covariate
eff.dat <- read.csv(here("./Data/Generated_DFs/Occ_Mod_Data/2021_OccEffortFileSubset.csv"))
eff.dat <- eff.dat[,-1]
colnames(eff.dat) <- gsub("[[:punct:]]", "_", gsub("X", "", colnames(eff.dat)))

## Summarize the effort data at the same temporal interval
eff.days <- second_samp(DAT = eff.dat, interval = 6, id_col = "Cell_Unit", eff = T, e.var = "Days")
colnames(eff.days)[colnames(eff.days) == "V1"] <- "Cell_Unit"
eff.days <- eff.days[,-1]
colnames(eff.days) <- NULL

## Summarize the total number of hours surveyed per sampling unit
eff.hrs <- second_samp(DAT = eff.dat, interval = 6, id_col = "Cell_Unit", eff = T, e.var = "Hrs")
colnames(eff.hrs)[colnames(eff.hrs) == "V1"] <- "Cell_Unit"
#eff.hrs[eff.hrs == 0] <- NA
eff.hrs <- eff.hrs[,-1]
colnames(eff.hrs) <- NULL


## Summarize the median Jdate
eff.jday <- second_samp(DAT = eff.dat, interval = 6, id_col = "Cell_Unit", eff = T, e.var = "FirstJDay")
colnames(eff.jday)[colnames(eff.jday) == "V1"] <- "Cell_Unit"
#eff.jday[is.na(eff.days)] <- NA
eff.jday <- eff.jday[,-1]
colnames(eff.jday) <- NULL

## Set y to NA for days without sampling
# Iterate through the matrix and set corresponding values in the 3D array to NA
# Iterate through the matrix and set corresponding values in the 3D array to NA for all layers
for (i in 1:nrow(eff.days)) {
  for (j in 1:ncol(eff.days)) {
    if (as.numeric(eff.days[i, j]) == 0) {
      # Set all layers (third dimension) of the array at position (i, j) to NA
      y[i, j, ] <- NA
    }
  }
}

# Print the modified 3D array
print(y)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: "Ragged array" data input for skipping NAs
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## *************************************************************
##
## Section Notes:
## We want a way to skip the pesky NAs that comprise ~30% of the
## data and will inevitably slow the model down. The solution
## is to link the data only the the sites with sampling effort
## and skip the NAs. Because the ARUs don't sample at a given site
## this creates a similar scheme for all species and can be handled
## for one species and attributed to all remaining species.
##
## *************************************************************

## One slice of the species matrix
y_sub <- y[,,1]

## Find the indices with NAs
has_data <- which(
  !is.na(y_sub),
  arr.ind = T
)

## Make the covariates into long format
obs_cov_long <- matrix(
  NA,
  nrow(has_data),
  ncol = 3
)

## Place the correct data into the new long format
eff.hrs <- as.matrix(eff.hrs)
eff.days <- as.matrix(eff.days)
eff.jday <- as.matrix(eff.jday)

## Make these in long format
for(i in 1:nrow(has_data)){
  obs_cov_long[i,1] <- eff.hrs[
    has_data[i,1], # site
    #1,
    has_data[i,2] # rep
  ]
  
  obs_cov_long[i,2] <- eff.days[
    has_data[i,1], # site
    #1,
    has_data[i,2] # rep
  ]
  
  obs_cov_long[i,3] <- eff.jday[
    has_data[i,1], # site
    #1,
    has_data[i,2] # rep
  ]
  
}

head(obs_cov_long)

## Using body mass as a detection covariate



## Scale the detection covariates
eff.hrs <- obs_cov_long[,1]
mean.eff.hrs <- mean(eff.hrs)
sd.eff.hrs <- sd(eff.hrs)
eff.hrs.scale <- (eff.hrs - mean.eff.hrs) / sd.eff.hrs

eff.days <- obs_cov_long[,2]
mean.eff.days <- mean(eff.days)
sd.eff.days <- sd(eff.days)
eff.days.scale <- (eff.days - mean.eff.days) / sd.eff.days

eff.jday <- obs_cov_long[,3]
mean.eff.jday <- mean(eff.jday)
sd.eff.jday <- sd(eff.jday)
eff.jday.scale <- (eff.jday - mean.eff.jday) / sd.eff.jday

# {
# eff.hrs <- scale(eff.hrs)
# attr(eff.hrs, "scaled:center") <- NULL
# attr(eff.hrs, "scaled:scale") <- NULL
# 
# eff.days <- scale(eff.days)
# attr(eff.days, "scaled:center") <- NULL
# attr(eff.days, "scaled:scale") <- NULL
# }
#eff.jday <- as.matrix(scale(unlist(eff.))

## Make the response variable in long format
y_long <- matrix(data = NA, nrow = length(y_sub[!is.na(y_sub)]), ncol = dim(y)[3])
for(i in 1:dim(y)[3]){
  y.tmp <- y[,,i]
  y_long[,i] <- y.tmp[!is.na(y.tmp)]
}

## -------------------------------------------------------------
##
## Begin Section: Occupancy Covariates
##
## -------------------------------------------------------------

## Right now lets pull in the ARU meta data
aru_meta <- readr::read_csv(here("Data/ARU_120m.csv"))
aru_meta$Cell_Unit <- paste0(aru_meta$cell_id, "_", aru_meta$unit_numbe)

## Select the variables we are interested in for detection and occupancy
aru_meta <- aru_meta |> select(Cell_Unit, 
                               utme, utmn,
                               X, Y,
                               topo_elev, topo_tpi,
                               tmx_bcm_mn, ppt_bcm_mn,
                               fire1yr_cbi_mn, fire2_5yr_cbi_mn,
                               fire6_10yr_cbi_mn, fire11_35yr_cbi_mn,
                               standage_f3_mn, cpycovr_f3_mn) |> 
  filter(Cell_Unit %in% cu.map$Cell_Unit) |> 
  arrange(match(Cell_Unit, cu.map$Cell_Unit))

## Scale the preds of interest
{
utme <- as.vector(scale(aru_meta$utme))
utmn <- as.vector(scale(aru_meta$utmn))
ele <- as.vector(scale(aru_meta$topo_elev))
ppt <- as.vector(scale(aru_meta$ppt_bcm_mn))
tmx <- as.vector(scale(aru_meta$tmx_bcm_mn))
cbi1 <- as.vector(scale(aru_meta$fire1yr_cbi_mn))
cbi2_5 <- as.vector(scale(aru_meta$fire2_5yr_cbi_mn))
cbi6_10 <- as.vector(scale(aru_meta$fire6_10yr_cbi_mn))
cbi11_35 <- as.vector(scale(aru_meta$fire11_35yr_cbi_mn))
stage <- as.vector(scale(aru_meta$standage_f3_mn)) 
cc <- as.vector(scale(aru_meta$cpycovr_f3_mn))
}

## -------------------------------------------------------------
##
## Begin Section: Prep Data for JAGS
##
## -------------------------------------------------------------

dimnames(y) <- NULL
win.data <- list(y = y, 
             nsite = dim(y)[1],
             nrep = dim(y)[2],
             nspec = dim(y)[3],
             eff.days = eff.days,
             eff.hrs = eff.hrs,
             utmn = utmn,
             ele = ele,
             ppt = ppt,
             tmx = tmx,
             cbi1 = cbi1,
             cbi2_5 = cbi2_5,
             cbi6_10 = cbi6_10,
             cbi11_35 = cbi11_35,
             stage = stage,
             cc = cc
             )
str(win.data)

## Win data new
win.data.new <- list(y = y_long,
                     nsite = dim(y)[1],
                     N = nrow(y_long),
                     nspec = ncol(y_long),
                     site_id = has_data[,1],
                     eff.hrs = eff.hrs.scale,
                     eff.jday = eff.jday.scale,
                     utmn = utmn,
                     ele = ele,
                     ppt = ppt,
                     tmx = tmx,
                     cbi1 = cbi1,
                     cbi2_5 = cbi2_5,
                     cbi6_10 = cbi6_10,
                     cbi11_35 = cbi11_35,
                     stage = stage,
                     cc = cc
                     )
str(win.data.new)
#                  nrep = dim(y)[2],
#                  nspec = dim(y)[3],
#                  eff.days = eff.days,
#                  eff.hrs = eff.hrs,
#                  utmn = utmn,
#                  ele = ele,
#                  ppt = ppt,
#                  tmx = tmx,
#                  cbi1 = cbi1,
#                  cbi2_5 = cbi2_5,
#                  cbi6_10 = cbi6_10,
#                  cbi11_35 = cbi11_35,
#                  stage = stage,
#                  cc = cc
# )


## -------------------------------------------------------------
##
## End Section:
##
## -------------------------------------------------------------


sp.det <- lapply(sp.det.list.r, function(x) x |> select(all_of(samp.cols)))
y <- array(unlist(lapply(sp.det, as.matrix)), dim = c(nsite, nrep, nspec))
dimnames(y) <- list(NULL, NULL, names(sp.det.list))

## No NAs present in the data
table(nsurveys <- apply(y[,,1], 1, function(x) sum(!is.na(x))))

## Species with 0 occurrences
## 9 species
tmp <- apply(y, c(1,3), max, na.rm = TRUE)
tmp[tmp == -Inf] <- NA
sort(obs.occ <- apply(tmp, 2, sum, na.rm = TRUE))

drop.sp <- which(obs.occ == 0)
y <- y[,,-drop.sp]
dimnames(y)

## Redefine nspec
nspec <- dim(y)[3]

# Get observed number of species per site
tmp <- apply(y, c(1,3), max, na.rm = TRUE)
tmp[tmp == "-Inf"] <- NA
sort(C <- apply(tmp, 1, sum)) # Compute and print sorted species counts

## We have an average of 20 species detected across the sites
plot(table(C), xlim = c(0, 60), xlab = "Observed number of species", ylab = "Number of
     quadrats", frame = F)
abline(v = mean(C, na.rm = TRUE), col = "blue", lwd = 3)

## Generate the NA data
ysum <- apply(y, c(1,3), sum, na.rm = T)

## Let's try a really simple model
win.data <- list(ysum = ysum, M = nrow(ysum), 
                 J = ,
                 nspec = dim(ysum)[2])

# Specify model in BUGS language
sink("model5.txt")
cat("
model {
# Priors
for(k in 1:nspec){ # Loop over species
psi[k] ~ dunif(0, 1)
p[k] ~ dunif(0, 1)
}
# Ecological model for latent occurrence z (process model)
for(k in 1:nspec){ # Loop over species
for (i in 1:M) { # Loop over sites
z[i,k] ~ dbern(psi[k])
}
}
# Observation model for observed data ysum
for(k in 1:nspec){ # Loop over species
for (i in 1:M) {
mup[i,k] <- z[i,k] * p[k]
ysum[i,k] ~ dbin(mup[i,k], J[i])
}
}
# Derived quantities
for(k in 1:nspec){ # Loop over species
Nocc.fs[k] <- sum(z[,k]) # Add up number of occupied sites among the 267
}
for (i in 1:M) { # Loop over sites
Nsite[i] <- sum(z[i,]) # Add up number of occurring species at each site
}
}
",fill = TRUE)
sink()

# Initial values
zst <- apply(y, c(1,3), max) # Observed occurrence as inits for z
zst[is.na(zst)] <- 1
inits <- function() list(z = zst, psi = rep(0.4, nspec), p = rep(0.4, nspec))
# Parameters monitored
params <- c("psi", "p", "Nsite", "Nocc.fs")

# MCMC settings
ni <- 2500 ; nt <- 2 ; nb <- 500 ; nc <- 3
# Call JAGS from R (ART 2.1 min)
out5 <- jags(win.data, inits, params, "model5.txt", n.chains = nc, n.thin = nt,
             n.iter = ni, n.burnin = nb, parallel = TRUE)
par(mfrow = c(4,4)) ; traceplot(out5) ; print(out5, dig = 3)


# and in a plot
cbind(obs.occu = obs.occ, out5$summary[which(str_detect(rownames(out5$summary), "Nsite")), c(1,3,7)])
plot(obs.occ, out5$summary[558:702, 1], xlab = "Observed number of occupied sites",
     ylab = "Estimated version of quantity", ylim = c(0, 267), frame = F, pch = 16)
abline(0,1)
segments(obs.occ, out5$summary[558:702,3], obs.occ, out5$summary[558:702,7], col
         = "grey", lwd = 2)


plot(out5$summary[1:92,1], out5$summary[93:184,1], xlab = "Occupancy estimate",
     ylab = "Detection estimate", xlim = c(0,1), ylim = c(0,1), frame = F, pch = 16)
segments(out5$summary[1:92,3], out5$summary[93:184,1], out5$summary[1:92,7],
         out5$summary[93:184,1], col = "grey", lwd = 2)
segments(out5$summary[1:92,1], out5$summary[93:184,3], out5$summary[1:92,1],
         out5$summary[93:184,7], col = "grey", lwd = 2)
## Looking at the posterior
rich.df <- data.frame(mnN = out5$mean$Nsite,
                      lN = out5$q2.5$Nsite,
                      uN = out5$q97.5$Nsite,
                      obsN = apply(apply(y, c(1,3), max), 1, sum))

ggplot(data = rich.df) + 
  geom_point(aes(x = obsN, y = mnN)) +
  geom_errorbar(aes(ymin = lN, ymax = uN)) + 
  geom_abline(slope = 1) + 
  theme_bw()
