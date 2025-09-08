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
n.factors <- 2

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


batch.length <- 25
n.batch <- 800
n.burn <- 5000
n.thin <- 20
n.chains <- 3

n.batch*batch.length

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
n.report <- 200 # Report progress at every 50th batch.

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
any(out.sfMsPGOcc$rhat$beta.comm > 1.1)
any(out.sfMsPGOcc$rhat$alpha.comm > 1.1)
any(out.sfMsPGOcc$rhat$tau.sq.beta > 1.1)
any(out.sfMsPGOcc$rhat$tau.sq.alpha > 1.1)
any(out.sfMsPGOcc$rhat$beta > 1.1)
which(out.sfMsPGOcc$rhat$beta > 1.1)
out.sfMsPGOcc$rhat$beta[which(out.sfMsPGOcc$rhat$beta > 1.1)]
any(out.sfMsPGOcc$rhat$alpha > 1.1)
which(out.sfMsPGOcc$rhat$alpha > 1.1)
out.sfMsPGOcc$rhat$alpha[which(out.sfMsPGOcc$rhat$alpha > 1.1)]
any(out.sfMsPGOcc$rhat$theta > 1.1)
any(out.sfMsPGOcc$rhat$lambda.lower.tri > 1.1)
which(out.sfMsPGOcc$rhat$lambda.lower.tri > 1.1)
any(out.sfMsPGOcc$rhat$sigma.sq.p > 1.1)

## ***********************************************************
##
## Section Notes:
## Running the model with 5 latent factors did not improve model
## fit...
## Going to try with 2 LF next.
## If this fails, need to consider reordering species
## or running with 1 chain and assessing convergence this way
##
## ***********************************************************

# Takes a few seconds to run. 
ppc.sfms.out <- ppcOcc(out.sfMsPGOcc, 'freeman-tukey', group = 1)
ppc.sfms.out.rep <- ppcOcc(out.sfMsPGOcc, 'freeman-tukey', group = 2)
ppc.sfms.out.cs <- ppcOcc(out.sfMsPGOcc, 'chi-squared', group = 1)
ppc.sfms.out.rep.cs <- ppcOcc(out.sfMsPGOcc, 'chi-squared', group = 2)

summary(ppc.sfms.out)
summary(ppc.sfms.out.rep)

## Compare to the non-spatial MSOM
waicOcc(out.sfMsPGOcc)

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
abline(h = 0.1, lty = "dashed", col = "red")
## Save the model object
print(object.size(out.sfMsPGOcc), units = "GB", standard = "SI") # For gigabytes (10^9 bytes)

## Sample the Z-matrix
z.mat <- out.sfMsPGOcc$z.samples
str(z.mat)

## Rearrange the Zmatrix to be more condusive the earlier code
## Sites x Species x Draw
z.mat <- aperm(z.mat, c(3,2,1))

saveRDS(z.mat, file = "F:/skeyser/PostDoc/SpOccupancy_Data/Z_Matrix_sfMsPGOcc_2slf.rds")

## Save the model object
saveRDS(out.sfMsPGOcc, file = "F:/skeyser/PostDoc/SpOccupancy_Data/Model_sfMsPGOcc_2slf.rds")
