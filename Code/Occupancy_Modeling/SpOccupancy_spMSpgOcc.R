## -------------------------------------------------------------
##
## Script name: spOccupancy Model Fitting
##
## Script purpose:
##
## Author: Spencer R Keyser
##
## Date Created: 2025-05-02
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
options(scipen = 10, digits = 10)

## -------------------------------------------------------------

## Package Loading
library(dplyr)
library(ggplot2)
library(here)
library(spOccupancy)

## -------------------------------------------------------------

## Load the data for spOccupancy Multispecies Spatial Models
spOcc.data <- readRDS(here("./Data/SpOccupancy_Data/SpOccData_SpeciesThresh_975minMaxPrex_NewVars.rds"))
spOcc.data <- readRDS("F:/skeyser/PostDoc/SpOccupancy_Data/SpOccData_SpeciesThresh_975minMaxPrex_NewVars.rds")
spOcc.data$occ.covs$lat <- spOcc.data$coords$Y

ncol(spOcc.data$det.covs$eff.hrs)
cc_cfo_dt <- spOcc.data$occ.covs$cc_cfo
colnames(cc_cfo_dt) <- NULL
str(spOcc.data$det.covs$eff.hrs)
cc_cfo_df <- as.matrix(cc_cfo_dt)
spOcc.data$det.covs$cc_cfo <- cc_cfo_dt

## Site_ID numeric map
site.ids <- data.frame(Sites = dimnames(spOcc.data$y)$sites,
                       Site.id = 1:length(unique(dimnames(spOcc.data$y)$sites)))

## Random effect for sites
str(site.ids)
spOcc.data$det.covs$site.id <- site.ids$Site.id

## Specify the fomula
occ.ms.formula <- ~ scale(lat) + ele + I(ele^2) + ppt + cbi1_5 + cbi6_10 + cbi11_35 + cc_cfo + ch_res
det.ms.formula <- ~ scale(eff.jday) + I(scale(eff.jday)^2) + scale(eff.hrs) + cc_cfo + (1|site.id) 

## Set some initial values
N <- dim(spOcc.data$y)[1]

## Non-spatial multi-species model
ms.inits <- list(alpha.comm = 0, 
                 beta.comm = 0, 
                 beta = 0, 
                 alpha = 0,
                 tau.sq.beta = 1, 
                 tau.sq.alpha = 1, 
                 z = apply(spOcc.data$y, c(1, 2), max, na.rm = TRUE))

ms.priors <- list(beta.comm.normal = list(mean = 0, var = 2.72),
                  alpha.comm.normal = list(mean = 0, var = 2.72), 
                  tau.sq.beta.ig = list(a = 0.1, b = 0.1), 
                  tau.sq.alpha.ig = list(a = 0.1, b = 0.1))

## Run multi-species model
out.ms <- msPGOcc(occ.formula = occ.ms.formula, 
                  det.formula = det.ms.formula, 
                  data = spOcc.data, 
                  inits = ms.inits, 
                  n.samples = 30000, 
                  priors = ms.priors, 
                  n.omp.threads = 1, 
                  verbose = TRUE, 
                  n.report = 6000, 
                  n.burn = 10000,
                  n.thin = 50, 
                  n.chains = 3)

summary(out.ms, level = "community")
summary(out.ms, level = "species")

## Posterior predictive checking
ppc.ms.ft <- ppcOcc(out.ms, fit.stat = "freeman-tukey", group = 1)
ppc.ms.cs <- ppcOcc(out.ms, fit.stat = "chi-squared", group = 1)

ppc.ms.ft.rep <- ppcOcc(out.ms, fit.stat = "freeman-tukey", group = 2)
ppc.ms.cs.rep <- ppcOcc(out.ms, fit.stat = "chi-squared", group = 2)


summary(ppc.ms.ft)
summary(ppc.ms.ft.rep)

## Plotting fit statistics for the two different groups
ppc.df <- data.frame(fit = ppc.ms.ft.rep$fit.y, 
                     fit.rep = ppc.ms.ft.rep$fit.y.rep, 
                     color = 'lightskyblue1')
ppc.df$color[ppc.df$fit.rep.1 > ppc.df$fit.1] <- 'lightsalmon'
plot(ppc.df$fit.1, ppc.df$fit.rep.1, bg = ppc.df$color, pch = 21, 
     ylab = 'Fit', xlab = 'True')
lines(ppc.df$fit.1, ppc.df$fit.1, col = 'black')

## Look at the fitted values for detection probability
fit.ms <- fitted(out.ms)
str(fit.ms)
apply(fit.ms$p.samples, 2, mean)

# Assuming fit.ms is your fitted object
# Dimensions: [posterior samples, species, sites, visits]

# Calculate mean detection probability per species
n.species <- dim(fit.ms$p.samples)[2]
det.prob.species <- numeric(n.species)

# Mean across all posterior samples, sites, and visits
for(i in 1:n.species) {
  det.prob.species[i] <- mean(fit.ms$p.samples[,i,,], na.rm=TRUE)
}

# With credible intervals
det.prob.CI <- matrix(NA, nrow=n.species, ncol=3)
for(i in 1:n.species) {
  # Flatten arrays for each species
  p.flat <- as.vector(fit.ms$p.samples[,i,,])
  det.prob.CI[i,] <- quantile(p.flat, probs=c(0.025, 0.5, 0.975), na.rm=TRUE)
}

# Create summary dataframe
det.summary <- data.frame(
  Species = dimnames(spOcc.data$y)[1],
  Mean_Det = det.prob.species,
  Lower_CI = det.prob.CI[,1],
  Median = det.prob.CI[,2],
  Upper_CI = det.prob.CI[,3]
)

hist(det.summary$Mean_Det, breaks = 20)

## Bayesian p-value
bp <- ppc.ms.ft$fit.y.rep > ppc.ms.ft$fit.y
bp <- apply(bp, 2, FUN = function(x) sum(x)/length(x))

## Put this with the probability of detection
det.summary$BP <- bp
cor(det.summary$Mean_Det, det.summary$BP, method = "spearman")
plot(det.summary$Mean_Det, det.summary$BP)
abline(h = 0.1, lty = "dashed", col = "red")

## Naive occupancy
n.occ <- apply(spOcc.data$y, 1, mean, na.rm = T)
str(as.vector(n.occ))
## Add naive occ
det.summary$nOcc <- as.vector(n.occ)
plot(det.summary$nOcc, det.summary$BP)
abline(h = 0.1, lty = "dashed", col = "red")

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Spatial Multispecies Occ Models
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Specify the fomula
occ.ms.sp.formula <- ~ ele + I(ele^2) + ppt + cbi1_5 + cbi6_10 + cbi11_35 + cc_cfo + ch_res
det.ms.sp.formula <- ~ scale(eff.jday) + I(scale(eff.jday)^2) + scale(eff.hrs) + (1|site.id) 

# Number of species
N <- dim(spOcc.data$y)[1]
# Distances between sites
spOcc.data$coords <- spOcc.data$coords 
dist.sierra <- dist(spOcc.data$coords)/1000
# Exponential covariance model
cov.model <- "exponential"
ms.inits <- list(alpha.comm = 0, 
                 beta.comm = 0, 
                 beta = 0, 
                 alpha = 0,
                 tau.sq.beta = 1, 
                 tau.sq.alpha = 1, 
                 z = apply(spOcc.data$y, c(1, 2), max, na.rm = TRUE), 
                 sigma.sq = 2, 
                 phi = 3 / mean(dist.sierra), 
                 w = matrix(0, N, dim(spOcc.data$y)[2]))

# Minimum value is 0, so need to grab second element.
min.dist <- sort(unique(dist.sierra))[2]
max.dist <- max(dist.sierra)
ms.priors <- list(beta.comm.normal = list(mean = 0, var = 2.72),
                  alpha.comm.normal = list(mean = 0, var = 2.72), 
                  tau.sq.beta.ig = list(a = 0.1, b = 0.1), 
                  tau.sq.alpha.ig = list(a = 0.1, b = 0.1),
                  sigma.sq.ig = list(a = 2, b = 2), 
                  phi.unif = list(a = 3 / max.dist, b = 3 / min.dist))

batch.length <- 25
n.batch <- 400
n.burn <- 2000
n.thin <- 20
n.chains <- 3
ms.tuning <- list(phi = 0.5)
n.omp.threads <- 1
# Values for reporting
verbose <- TRUE
n.report <- 100

# Approx. run time: 10 min
out.sp.ms <- spMsPGOcc(occ.formula = occ.ms.sp.formula, 
                       det.formula = det.ms.sp.formula, 
                       data = spOcc.data, 
                       inits = ms.inits, 
                       n.batch = n.batch, 
                       batch.length = batch.length, 
                       accept.rate = 0.43, 
                       priors = ms.priors, 
                       cov.model = cov.model, 
                       tuning = ms.tuning, 
                       n.omp.threads = n.omp.threads, 
                       verbose = TRUE, 
                       NNGP = TRUE, 
                       n.neighbors = 5, 
                       n.report = n.report, 
                       n.burn = n.burn, 
                       n.thin = n.thin, 
                       n.chains = n.chains)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Spatial Latent Factor MSOM
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ***********************************************************
##
## Section Notes:
## The spatial MSOM takes a prohibitively long time to fit.
## Perhaps shifting to the faster fitting latent spatial factor
## MSOM is a fruitful endeavor.
##
## ***********************************************************
occ.splf.formula <- ~ ele + I(ele^2) + ppt + cbi1_5 + cbi6_10 + cbi11_35 + cc_cfo + ch_res
det.splf.formula <- ~ scale(eff.jday) + I(scale(eff.jday)^2) + scale(eff.hrs) + (1|site.id)

## Reorder the data according to the factors of interest
## That is...
## Species 1 - 
sort(apply(spOcc.data$y, 1, mean, na.rm = TRUE))

## Set the factors
n.factors <- 3

## Set a common species first
sp.trait <- readRDS("F:/skeyser/PostDoc/SpOccupancy_Data/TrialSpeciesOrder.rds")

sp.add <- sp.trait[sp.trait %in% unlist(dimnames(spOcc.data$y)[1])]


sp.reorder <- c("Red-breasted Nuthatch", 
                sp.add,
                "Sphyrapicus spp.",
                "Vireo spp.")

# Create new detection-nondetection data matrix in the new order
y.new <- spOcc.data$y[sp.reorder, , ]
# Create a new data array
spOcc.data.ordered <- spOcc.data
# Change the data to the new ordered data
spOcc.data.ordered$y <- y.new
str(spOcc.data.ordered)

# Pair-wise distance between all sites
dist.sierra <- dist(spOcc.data.ordered$coords)/1000
# Exponential correlation model
cov.model <- "exponential"
# Specify all other initial values identical to lfMsPGOcc() from before
# Number of species
N <- nrow(spOcc.data.ordered$y)
# Initiate all lambda initial values to 0. 
lambda.inits <- matrix(0, N, n.factors)
# Set diagonal elements to 1
diag(lambda.inits) <- 1
# Set lower triangular elements to random values from a standard normal dist
lambda.inits[lower.tri(lambda.inits)] <- rnorm(sum(lower.tri(lambda.inits)))
# Check it out
lambda.inits

# Create list of initial values. 
inits <- list(alpha.comm = 0,
              beta.comm = 0,
              beta = 0,
              alpha = 0,
              tau.sq.beta = 1,
              tau.sq.alpha = 1,
              lambda = lambda.inits, 
              phi = 3 / mean(dist.sierra),
              z = apply(spOcc.data.ordered$y, c(1, 2), max, na.rm = TRUE))

priors <- list(beta.comm.normal = list(mean = 0, var = 2.72),
               alpha.comm.normal = list(mean = 0, var = 2.72),
               tau.sq.beta.ig = list(a = 0.1, b = 0.1),
               tau.sq.alpha.ig = list(a = 0.1, b = 0.1))

batch.length <- 25
n.batch <- 400
n.burn <- 5000
n.thin <- 10
n.chains <- 3

(((n.batch*batch.length-n.burn)/n.thin))

min.dist <- min(dist.sierra)
max.dist <- max(dist.sierra)
priors <- list(beta.comm.normal = list(mean = 0, var = 2.72),
               alpha.comm.normal = list(mean = 0, var = 2.72),
               tau.sq.beta.ig = list(a = 0.1, b = 0.1),
               tau.sq.alpha.ig = list(a = 0.1, b = 0.1), 
               phi.unif = list(3 / max.dist, 3 / min.dist))

tuning <- list(phi = 1)

n.omp.threads <- 1
verbose <- TRUE
n.report <- 50 # Report progress at every 50th batch.

# Approx run time: 2 min
out.sfMsPGOcc <- sfMsPGOcc(occ.formula = occ.splf.formula, 
                           det.formula = det.splf.formula, 
                           data = spOcc.data.ordered, 
                           inits = inits, 
                           n.batch = n.batch, 
                           batch.length = batch.length, 
                           accept.rate = 0.43, 
                           priors = priors, 
                           n.factors = n.factors,
                           cov.model = cov.model, 
                           tuning = tuning, 
                           n.omp.threads = n.omp.threads, 
                           verbose = TRUE, 
                           NNGP = TRUE, 
                           n.neighbors = 5, 
                           n.report = n.report, 
                           n.burn = n.burn, 
                           n.thin = n.thin, 
                           n.chains = n.chains)

summary(out.sfMsPGOcc)
summary(out.sfMsPGOcc$lambda.samples)
# Takes a few seconds to run. 
ppc.sfms.out <- ppcOcc(out.sfMsPGOcc, 'freeman-tukey', group = 1)
ppc.sfms.out.rep <- ppcOcc(out.sfMsPGOcc, 'freeman-tukey', group = 2)
ppc.sfms.out.cs <- ppcOcc(out.sfMsPGOcc, 'chi-squared', group = 1)
ppc.sfms.out.rep.cs <- ppcOcc(out.sfMsPGOcc, 'chi-squared', group = 2)

summary(ppc.sfms.out)
summary(ppc.sfms.out.rep)

## Compare to the non-spatial MSOM
waicOcc(out.sfMsPGOcc)
waicOcc(out.ms)

## Psi samples
str(out.sfMsPGOcc$psi.samples)
mean.psi <- apply(out.sfMsPGOcc$psi.samples, 2, mean, na.rm = T)
psi.low <- apply(out.sfMsPGOcc$psi.samples, 2, function(x) quantile(x, probs = 0.025, na.rm = T))
psi.hi <- apply(out.sfMsPGOcc$psi.samples, 2, function(x) quantile(x, probs = 0.975, na.rm = T))

sp.reorder
bp <- ppc.sfms.out$fit.y.rep > ppc.sfms.out$fit.y
bp <- apply(bp, 2, FUN = function(x) sum(x)/length(x))

psi.df <- data.frame(sp = sp.reorder,
                     mean.psi = mean.psi,
                     lo.psi = psi.low,
                     hi.psi = psi.hi,
                     naive.psi = as.vector(apply(spOcc.data.ordered$y, 1, mean, na.rm = T)),
                     bp = bp) |> 
  mutate(psi.range = hi.psi - lo.psi)

plot(psi.df$naive.psi, psi.df$mean.psi)
plot(psi.df$naive.psi, psi.df$psi.range)
plot(psi.df$mean.psi, psi.df$bp)

## Save the model object
print(object.size(out.sfMsPGOcc), units = "GB", standard = "SI") # For gigabytes (10^9 bytes)

## Sample the Z-matrix
z.mat <- out.sfMsPGOcc$z.samples
str(z.mat)

## Rearrange the Zmatrix to be more condusive the earlier code
## Sites x Species x Draw
z.mat <- aperm(z.mat, c(3,2,1))

#saveRDS(z.mat, file = "F:/skeyser/PostDoc/SpOccupancy_Data/Z_Matrix_sfMsPGOcc.rds")

## Save the model object
saveRDS(out.sfMsPGOcc, file = "F:/skeyser/PostDoc/SpOccupancy_Data/Model_sfMsPGOcc.rds")

## -----------------------------------------------------------
##
## Begin Section: Spatial Single Species Occupancy Models
##
## -----------------------------------------------------------

## ***********************************************************
##
## Section Notes: Our traditional MSOM is failing to pass PPC
## for specific species. Consistently underfitting models
## suggests that we are missing some critical predictors or
## we are inadequately capturing spatial processes structure
## site level differences in occupancy or detection. Currently,
## fit statistics are worse between binned data across sites.
## To try to assess these issues...a spatial component might be
## important. Additionally, integrating more information in either
## the occupancy or detection submodel could be useful.
##
## Below: Attempting to fit a spatial model for one species to 
## start. Acorn Woodpecker has been a persistent problem species
## so we will start with that. 
##
## ***********************************************************
spOcc.data.ssom <- spOcc.data

## Let's subset just the ACWO
sp.test <- "Hermit Warbler"
y <- spOcc.data$y
str(y)

y.ssom <- y[sp.test, , ]

spOcc.data.ssom$y <- y.ssom
sum(apply(y.ssom, 1, max, na.rm = T))

## Drop the erroneous sites to see how much of an effect they have
spOcc.data.ssom$occ.covs <- spOcc.data.ssom$occ.covs[-erroneous.sites,]
spOcc.data.ssom$coords <- spOcc.data.ssom$coords[-erroneous.sites,]
spOcc.data.ssom$y <- spOcc.data.ssom$y[-erroneous.sites,]
spOcc.data.ssom$det.covs$eff.hrs <- spOcc.data.ssom$det.covs$eff.hrs[-erroneous.sites,]
spOcc.data.ssom$det.covs$eff.jday <- spOcc.data.ssom$det.covs$eff.jday[-erroneous.sites,]
spOcc.data.ssom$det.covs$cc_cfo <- spOcc.data.ssom$det.covs$cc_cfo[-erroneous.sites]
spOcc.data.ssom$det.covs$site.id <- spOcc.data.ssom$det.covs$site.id[-erroneous.sites]

## Distances
avg.disp <- 70 * 1000
dist.sites <- dist(spOcc.data.ssom$coords)
cov.model <- "exponential"
ms.inits <- list(alpha.comm = 0,
                 beta.comm = 0,
                 beta = 0,
                 alpha = 0,
                 tau.sq.beta = 1,
                 tau.sq.alpha = 1,
                 z = apply(spOcc.data.ssom$y, 1, max, na.rm = T),
                 sigma.sq = 2,
                 phi = 3 / avg.disp, #6 km is average dispersal for ACWO
                 w = rep(0, nrow(spOcc.data.ssom$y)))
str(ms.inits)

# Minimum value is 0, so need to grab second element.
# min.dist <- sort(unique(dist.sites))[2]
# max.dist <- max(dist.sites)
# ms.priors <- list(beta.comm.normal = list(mean = 0, var = 2.72),
#                   alpha.comm.normal = list(mean = 0, var = 2.72), 
#                   tau.sq.beta.ig = list(a = 0.1, b = 0.1), 
#                   tau.sq.alpha.ig = list(a = 0.1, b = 0.1),
#                   sigma.sq.ig = list(a = 2, b = 2), 
#                   phi.unif = list(a = 3 / max.dist, b = 3 / min.dist))

batch.length <- 25
n.batch <- 400
n.samples <- n.batch * batch.length
n.burn <- 2000
n.thin <- 20 #50
n.chains <- 3
ms.tuning <- list(phi = 1)
n.omp.threads <- 1
# Values for reporting
verbose <- TRUE
n.report <- 100

((n.samples - n.burn) / n.thin) * n.chains

## ACWO priors
## Restrict the max distance to the appx. upper recorded dispersal (200 km for ACWO)
## See if this fixes issues with phi convergence
min.dist <- min(dist.sites)
max.dist <- max(dist.sites)/5
ssom.priors <- list(beta.normal = list(mean = 0, var = 2.72), 
                    alpha.normal = list(mean = 0, var = 2.72), 
                    sigma.sq.ig = c(2, 1), 
                    phi.unif = c(3/max.dist, 3/min.dist))

## Formula
occ.ss.formula <- ~ ele + I(ele^2) + ppt + cbi1_5 + cbi6_10 + cbi11_35 + cc_cfo + ch_res
det.ss.formula <- ~ scale(eff.jday) + I(scale(eff.jday)^2) + scale(eff.hrs) + (1|site.id)

# Approx. run time: 10 min
out.sp.ss <- spPGOcc(occ.formula = occ.ss.formula, 
                     det.formula = det.ss.formula, 
                     data = spOcc.data.ssom, 
                     inits = ms.inits, 
                     n.batch = n.batch, 
                     batch.length = batch.length, 
                     accept.rate = 0.43, 
                     priors = ssom.priors, 
                     cov.model = cov.model, 
                     tuning = ms.tuning, 
                     n.omp.threads = n.omp.threads, 
                     verbose = TRUE, 
                     NNGP = TRUE, 
                     n.neighbors = 5, 
                     n.report = n.report, 
                     n.burn = n.burn, 
                     n.thin = n.thin, 
                     n.chains = n.chains)

## Summary
summary(out.sp.ss)

str(out.sp.ss$psi.samples)
exp.psi <- mean(out.sp.ss$psi.samples)
naive.psi <- sum(apply(spOcc.data.ssom$y, 1, max, na.rm = T))/nrow(spOcc.data.ssom$y)
exp.psi - naive.psi
## PPC
ppc.sp.out <- ppcOcc(out.sp.ss, fit.stat = 'freeman-tukey', group = 1)
summary(ppc.sp.out)

ppc.sp.out <- ppcOcc(out.sp.ss, fit.stat = 'freeman-tukey', group = 2)
summary(ppc.sp.out)

ppc.sp.out <- ppcOcc(out.sp.ss, fit.stat = 'chi-squared', group = 1)
summary(ppc.sp.out)

ppc.sp.out <- ppcOcc(out.sp.ss, fit.stat = 'chi-squared', group = 2)
summary(ppc.sp.out)


## PPC Visualization
## Plotting fit statistics for the two different groups
ppc.df <- data.frame(fit = ppc.sp.out$fit.y, 
                     fit.rep = ppc.sp.out$fit.y.rep, 
                     color = 'lightskyblue1')
ppc.df$color[ppc.df$fit.rep > ppc.df$fit] <- 'lightsalmon'
plot(ppc.df$fit, ppc.df$fit.rep, bg = ppc.df$color, pch = 21, 
     ylab = 'Fit', xlab = 'True')
lines(ppc.df$fit, ppc.df$fit, col = 'black')

diff.fit <- ppc.sp.out$fit.y.rep.group.quants[3, ] - ppc.sp.out$fit.y.group.quants[3, ]
plot(diff.fit, pch = 19, xlab = 'Site ID', ylab = 'Replicate - True Discrepancy')

mean(diff.fit)
quantile(diff.fit, probs = seq(0,1,0.1))

erroneous.sites <- which(diff.fit < quantile(diff.fit, probs = 0.05))

site.issue <- rownames(y.acwo[erroneous.sites,])

spOcc.data.ssom$occ.covs[erroneous.sites,]

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Acorn Woodpecker Non-spatial
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
s.inits <- list(alpha.comm = 0,
                beta.comm = 0,
                beta = 0,
                alpha = 0,
                tau.sq.beta = 1,
                tau.sq.alpha = 1,
                z = apply(spOcc.data.ssom$y, 1, max, na.rm = T),
                sigma.sq = 2)
str(s.inits)

n.samples <- 30000
n.burn <- 10000 
n.thin <- 50 #50
n.chains <- 3
n.omp.threads <- 1
# Values for reporting
verbose <- TRUE
n.report <- 10000

((n.samples - n.burn) / n.thin) * n.chains

acwo.priors.ns <- list(beta.normal = list(mean = 0, var = 2.72), 
                       alpha.normal = list(mean = 0, var = 2.72), 
                       sigma.sq.ig = c(2, 1))

## Formula
occ.ss.formula <- ~ ele + I(ele^2) + ppt + cbi1_5 + cbi6_10 + cbi11_35 + cc_cfo + ch_res
det.ss.formula <- ~ scale(eff.jday) + I(scale(eff.jday)^2) + scale(eff.hrs) + (1|site.id)

# Approx. run time: 10 min
out.ns.ss <- PGOcc(occ.formula = occ.ss.formula, 
                   det.formula = det.ss.formula, 
                   data = spOcc.data.ssom, 
                   inits = s.inits, 
                   priors = acwo.priors.ns, 
                   n.omp.threads = n.omp.threads, 
                   verbose = TRUE,
                   n.samples = n.samples,
                   n.report = n.report, 
                   n.burn = n.burn, 
                   n.thin = n.thin, 
                   n.chains = n.chains)

## Summary
summary(out.ns.ss)

## PPC
ppc.ns.out <- ppcOcc(out.ns.ss, fit.stat = 'freeman-tukey', group = 1)
summary(ppc.ns.out)

ppc.ns.out <- ppcOcc(out.ns.ss, fit.stat = 'freeman-tukey', group = 2)
summary(ppc.ns.out)

ppc.ns.out <- ppcOcc(out.ns.ss, fit.stat = 'chi-squared', group = 1)
summary(ppc.ns.out)

ppc.ns.out <- ppcOcc(out.ns.ss, fit.stat = 'chi-squared', group = 2)
summary(ppc.ns.out)

## WAIC between the two model types
waicOcc(out.sp.ss)
waicOcc(out.ns.ss)


## Load in the data
acwo.naive <- apply(y.acwo,1,max, na.rm = T)

coords <- spOcc.data.ssom$coords |> 
  sf::st_as_sf(coords = c("X", "Y"), crs = 3310) |> 
  mutate(Cell_Unit = rownames(y.acwo)) |> 
  mutate(Status = ifelse(Cell_Unit %in% site.issue, "Bad", "Good")) |> 
  mutate(NOcc = acwo.naive) |> 
  bind_cols(spOcc.data.ssom$occ.covs)


ggplot(coords) + 
  geom_sf(aes(color = Status))

mapview::mapview(coords, zcol = "Status")

## Predictions for detection
out.sp.fit <- fitted(out.sp.ss)
str(out.sp.fit$p.samples)
which(is.na(out.sp.fit$p.samples))

## Site-level detection probability
mean.det <- mean(apply(out.sp.fit$p.samples, 2, mean, na.rm = T))
hist(apply(out.sp.fit$p.samples, 2, mean, na.rm = T))

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Testing with Hubbard Brook
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls())

data(hbef2015)

occ.ms.formula <- ~ scale(Elevation) + I(scale(Elevation)^2)
det.ms.formula <- ~ scale(day) + scale(tod) + I(scale(day)^2)
str(hbef2015)

N <- dim(hbef2015$y)[1]
ms.inits <- list(alpha.comm = 0, 
                 beta.comm = 0, 
                 beta = 0, 
                 alpha = 0,
                 tau.sq.beta = 1, 
                 tau.sq.alpha = 1, 
                 z = apply(hbef2015$y, c(1, 2), max, na.rm = TRUE))

ms.priors <- list(beta.comm.normal = list(mean = 0, var = 2.72),
                  alpha.comm.normal = list(mean = 0, var = 2.72), 
                  tau.sq.beta.ig = list(a = 0.1, b = 0.1), 
                  tau.sq.alpha.ig = list(a = 0.1, b = 0.1))

# Approx. run time:  6 min
out.ms <- msPGOcc(occ.formula = occ.ms.formula, 
                  det.formula = det.ms.formula, 
                  data = hbef2015, 
                  inits = ms.inits, 
                  n.samples = 30000, 
                  priors = ms.priors, 
                  n.omp.threads = 1, 
                  verbose = TRUE, 
                  n.report = 6000, 
                  n.burn = 10000,
                  n.thin = 50, 
                  n.chains = 3)
