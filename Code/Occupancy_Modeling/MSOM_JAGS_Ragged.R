## -------------------------------------------------------------
##
## Script name: MSOM JAGS Code
##
## Script purpose:
##
## Author: Spencer R Keyser
##
## Date Created: 2024-12-03
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
library(jagsUI)
library(rjags)
library(MCMCvis)
## -------------------------------------------------------------

## Load the data
load(here("./Data/JAGS_Data/MSOM_Ragged_2021_95cut.RData"))
load(here("./Data/JAGS_Data/MSOM_Ragged_2021_SpeciesThresh_975minMaxPrex_NewVars.RData"))

## -------------------------------------------------------------
##
## Begin Section: JAGS Model Code
##
## -------------------------------------------------------------

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Model 1 w/ "ragged array" for NAs
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sink("Code/JAGS_Models/Sierra_MSOM_Covs_NA_95_BM.txt")
cat("
    model {
    
    # Book-keeping: Ragged Array
    # i = site 
    # N = rep x site 
    # k = species
    # beta = coefficients for occupancy
    # alpha = coefficient for detection
    
    ############
    ## Priors ##
    ############
    
    for(k in 1:nspec){
      # Psi
      # Species random psi intercept
      lpsi[k] ~ dnorm(mu.lpsi, tau.lpsi) #uniform prior for species occ
      
      # Species random psi slopes
      beta1[k] ~ dnorm(mu.beta1, tau.beta1)
      beta2[k] ~ dnorm(mu.beta2, tau.beta2)
      beta3[k] ~ dnorm(mu.beta3, tau.beta3)
      beta4[k] ~ dnorm(mu.beta4, tau.beta4)
      beta5[k] ~ dnorm(mu.beta5, tau.beta5)
      beta6[k] ~ dnorm(mu.beta6, tau.beta6)
      beta7[k] ~ dnorm(mu.beta7, tau.beta7)
      beta8[k] ~ dnorm(mu.beta8, tau.beta8)
      beta9[k] ~ dnorm(mu.beta9, tau.beta9)
      #beta10[k] ~ dnorm(mu.beta10, tau.beta10)
      #beta11[k] ~ dnorm(mu.beta11, tau.beta11)
      
      # p
      # Species random p intercept
      lp[k] ~ dnorm(mu.lp, tau.lp)
      
      # Species random p slope
      alpha1[k] ~ dnorm(mu.alpha1, tau.alpha1)
      alpha2[k] ~ dnorm(mu.alpha2, tau.alpha2)
      alpha3[k] ~ dnorm(mu.alpha3, tau.alpha3)
      #alpha4[k] ~ dnorm(mu.alpha4, tau.alpha4)
    }
    
    #################
    ## Hyperpriors ##
    #################
    
    # Occ hypers
    # Intercept
    mu.lpsi ~ dnorm(0,0.01)
    tau.lpsi <- pow(sd.lpsi, -2)
    sd.lpsi ~ dunif(0,4) 
    
    # Latitude
    mu.beta1 ~ dnorm(0,0.1)
    tau.beta1 <- pow(sd.beta1, -2)
    sd.beta1 ~ dunif(0, 4)
    
    # Elevation 
    mu.beta2 ~ dnorm(0,0.1)
    tau.beta2 <- pow(sd.beta2, -2)
    sd.beta2 ~ dunif(0,4)
    
    # Elevation Poly
    mu.beta3 ~ dnorm(0,0.1)
    tau.beta3 <- pow(sd.beta3, -2)
    sd.beta3 ~ dunif(0,4)
    
    # Fire Sev 1yr
    mu.beta4 ~ dnorm(0,0.1)
    tau.beta4 <- pow(sd.beta4, -2)
    sd.beta4 ~ dunif(0,4)
    
    # Fire Sev 2-5 years
    mu.beta5 ~ dnorm(0,0.1)
    tau.beta5 <- pow(sd.beta5, -2)
    sd.beta5 ~ dunif(0,4)
    
    # Fire Sev 6-10 yrs
    mu.beta6 ~ dnorm(0,0.1)
    tau.beta6 <- pow(sd.beta6, -2)
    sd.beta6 ~ dunif(0,4)
    
    # Fire Sev 11-35 yrs
    mu.beta7 ~ dnorm(0,0.1)
    tau.beta7 <- pow(sd.beta7, -2)
    sd.beta7 ~ dunif(0,4)
    
    # Stand Age
    mu.beta8 ~ dnorm(0,0.1)
    tau.beta8 <- pow(sd.beta8, -2)
    sd.beta8 ~ dunif(0,4)
    
    # Canopy Cover
    mu.beta9 ~ dnorm(0,0.1)
    tau.beta9 <- pow(sd.beta9, -2)
    sd.beta9 ~ dunif(0,4)
    
    # Detection hypers
    
    # Intercept
    mu.lp ~ dnorm(0, 0.1)
    tau.lp <- pow(sd.lp, -2)
    sd.lp ~ dunif(0, 2)
    
    # Number of hrs sampled (per secondary sample)
    mu.alpha1 ~ dnorm(0, 0.1)
    tau.alpha1 <- pow(sd.alpha1, -2)
    sd.alpha1 ~ dunif(0, 2)
    
    # Number of hours sampled (per secondary sample)
    mu.alpha2 ~ dnorm(0, 0.1)
    tau.alpha2 <- pow(sd.alpha2, -2)
    sd.alpha2 ~ dunif(0, 2)
    
    # JDay for sampling
    mu.alpha3 ~ dnorm(0, 0.1)
    tau.alpha3 <- pow(sd.alpha3, -2)
    sd.alpha3 ~ dunif(0, 2)
    
    # Body mass for detection
    #mu.alpha4 ~ dnorm(0, 0.1)
    #tau.alpha4 <- pow(sd.alpha4, -2)
    #sd.alpha4 ~ dunif(0, 2)
    
    #################################################
    ## Ecological model for the latent process (z) ##
    #################################################
    # Add covariates for occupancy
    for(k in 1:nspec){ #species loop
      for(i in 1:nsite){ #site loop
      # Occupancy model w/ covs
      # Covs: Lat, Elevation, Elevation^2, Prop Burn Sev (PBS) 1yr, PBS 2-5, PBD 6-10, PBS 11-35, Stand Age, % Canopy Cover  
      logit(psi[i,k]) <- lpsi[k] + beta1[k] * utmn[i] + beta2[k] * ele[i] + beta3[k] * pow(ele[i], 2) + beta4[k] * cbi1[i] + beta5[k] * cbi2_5[i] + beta6[k] * cbi6_10[i] +
        beta7[k] * cbi11_35[i] + beta8[k] * stage[i] + beta9[k] * cc[i] 
      
      # True latent state (Z-matrix)
      z[i,k] ~ dbern(psi[i,k])
      
      # Model Assement via Chi-squared GoF
      evalZ[i,k] <- psi[i,k]
      EZ[i,k] <- pow((z[i,k] - evalZ[i,k]), 2) / (evalZ[i,k] + 0.5)
        
      # Replicated data for new comparison
      z.new[i,k] ~ dbern(psi[i,k])
      EZ.new[i,k] <- pow((z.new[i,k] - evalZ[i,k]), 2) / (evalZ[i,k] + 0.5)
      
      }
    }
    
    ###################################################
    ## Observation submodel for replicate det/nondet ##
    ###################################################
    
    ## Observation Sub-model with nested indexing on sites
    
    # Significant NAs in the response for variable ARU start dates
    # Detection heterogeneity is built into this model
    # Detection Covs: Survey Hrs (sum across 6 day secondary sampling), Jdate (median per 6 day interval), Jdate^2
    for(k in 1:nspec){ #columns y 
      for(j in 1:N){ # rows y
        # Detection model on logit scale
        logit(p[j,k]) <- lp[k] + alpha1[k] * eff.hrs[j] + alpha2[k] * eff.jday[j] + alpha3[k] * pow(eff.jday[j], 2) #+ alpha4[k] * bmass[k]
          
        # Latent state and detection
        # Site_id nested index
        mup[j,k] <- z[site_id[j], k] * p[j,k]
          
        # Observation Bernoulli draw from z*p
        y[j,k] ~ dbern(mup[j,k])
        
        # Model Assement via Chi-squared GoF
        eval[j,k] <- mup[j,k]
        E[j,k] <- pow((y[j,k] - eval[j,k]), 2) / (eval[j,k] + 0.5)
        
        # Replicated data for new comparison
        y.new[j,k] ~ dbern(mup[j,k])
        E.new[j,k] <- pow((y.new[j,k] - eval[j,k]), 2) / (eval[j,k] + 0.5)

      }
    }
    
    ########################################
    ## Derived estimates of the community ##
    ########################################
    
    # Species-specific # of occupied sites
    for(k in 1:nspec){
      Nocc.fs[k] <- sum(z[,k])
    }
    # Species Richness
    for(i in 1:nsite){
      Nsite[i] <- sum(z[i,])
    }
    
    fitZ <- sum(EZ[,])
    fitZ.new <- sum(EZ.new[,])
    fitY <- sum(E[,])
    fitY.new <- sum(E.new[,])
    
    }
    ", fill = TRUE)
sink()

## -------------------------------------------------------------
##
## End Section: JAGS Model Code
##
## -------------------------------------------------------------

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Model 2 w/ "ragged array" for NAs
## New covariates subbing in for now
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sink("Code/JAGS_Models/Sierra_MSOM_Covs_NA_VarThresh_975min.txt")
cat("
    model {
    
    # Book-keeping: Ragged Array
    # i = site 
    # N = rep x site 
    # k = species
    # beta = coefficients for occupancy
    # alpha = coefficient for detection
    
    ############
    ## Priors ##
    ############
    
    for(k in 1:nspec){
      # Psi
      # Species random psi intercept
      lpsi[k] ~ dnorm(mu.lpsi, tau.lpsi) #uniform prior for species occ
      
      # Species random psi slopes
      beta1[k] ~ dnorm(mu.beta1, tau.beta1)
      beta2[k] ~ dnorm(mu.beta2, tau.beta2)
      beta3[k] ~ dnorm(mu.beta3, tau.beta3)
      beta4[k] ~ dnorm(mu.beta4, tau.beta4)
      beta5[k] ~ dnorm(mu.beta5, tau.beta5)
      beta6[k] ~ dnorm(mu.beta6, tau.beta6)
      beta7[k] ~ dnorm(mu.beta7, tau.beta7)
      beta8[k] ~ dnorm(mu.beta8, tau.beta8)
      beta9[k] ~ dnorm(mu.beta9, tau.beta9)
      #beta10[k] ~ dnorm(mu.beta10, tau.beta10)
      #beta11[k] ~ dnorm(mu.beta11, tau.beta11)
      
      # p
      # Species random p intercept
      lp[k] ~ dnorm(mu.lp, tau.lp)
      
      # Species random p slope
      alpha1[k] ~ dnorm(mu.alpha1, tau.alpha1)
      alpha2[k] ~ dnorm(mu.alpha2, tau.alpha2)
      alpha3[k] ~ dnorm(mu.alpha3, tau.alpha3)
      #alpha4[k] ~ dnorm(mu.alpha4, tau.alpha4)
    }
    
    #################
    ## Hyperpriors ##
    #################
    
    # Occ hypers
    # Intercept
    mu.lpsi ~ dnorm(0,0.01)
    tau.lpsi <- pow(sd.lpsi, -2)
    sd.lpsi ~ dunif(0,4) 
    
    # Latitude
    mu.beta1 ~ dnorm(0,0.1)
    tau.beta1 <- pow(sd.beta1, -2)
    sd.beta1 ~ dunif(0, 4)
    
    # Elevation 
    mu.beta2 ~ dnorm(0,0.1)
    tau.beta2 <- pow(sd.beta2, -2)
    sd.beta2 ~ dunif(0,4)
    
    # Elevation Poly
    mu.beta3 ~ dnorm(0,0.1)
    tau.beta3 <- pow(sd.beta3, -2)
    sd.beta3 ~ dunif(0,4)
    
    # Precipitation
    mu.beta4 ~ dnorm(0,0.1)
    tau.beta4 <- pow(sd.beta4, -2)
    sd.beta4 ~ dunif(0,4)
    
    # Fire Sev 1-5 years
    mu.beta5 ~ dnorm(0,0.1)
    tau.beta5 <- pow(sd.beta5, -2)
    sd.beta5 ~ dunif(0,4)
    
    # Fire Sev 6-10 yrs
    mu.beta6 ~ dnorm(0,0.1)
    tau.beta6 <- pow(sd.beta6, -2)
    sd.beta6 ~ dunif(0,4)
    
    # Fire Sev 11-35 yrs
    mu.beta7 ~ dnorm(0,0.1)
    tau.beta7 <- pow(sd.beta7, -2)
    sd.beta7 ~ dunif(0,4)
    
    # Stand Age
    mu.beta8 ~ dnorm(0,0.1)
    tau.beta8 <- pow(sd.beta8, -2)
    sd.beta8 ~ dunif(0,4)
    
    # Canopy Cover
    mu.beta9 ~ dnorm(0,0.1)
    tau.beta9 <- pow(sd.beta9, -2)
    sd.beta9 ~ dunif(0,4)
    
    # Detection hypers
    
    # Intercept
    mu.lp ~ dnorm(0, 0.1)
    tau.lp <- pow(sd.lp, -2)
    sd.lp ~ dunif(0, 2)
    
    # Number of hrs sampled (per secondary sample)
    mu.alpha1 ~ dnorm(0, 0.1)
    tau.alpha1 <- pow(sd.alpha1, -2)
    sd.alpha1 ~ dunif(0, 2)
    
    # Number of hours sampled (per secondary sample)
    mu.alpha2 ~ dnorm(0, 0.1)
    tau.alpha2 <- pow(sd.alpha2, -2)
    sd.alpha2 ~ dunif(0, 2)
    
    # JDay for sampling
    mu.alpha3 ~ dnorm(0, 0.1)
    tau.alpha3 <- pow(sd.alpha3, -2)
    sd.alpha3 ~ dunif(0, 2)
    
    # Body mass for detection
    #mu.alpha4 ~ dnorm(0, 0.1)
    #tau.alpha4 <- pow(sd.alpha4, -2)
    #sd.alpha4 ~ dunif(0, 2)
    
    #################################################
    ## Ecological model for the latent process (z) ##
    #################################################
    # Add covariates for occupancy
    for(k in 1:nspec){ #species loop
      for(i in 1:nsite){ #site loop
      # Occupancy model w/ covs
      # Covs: Lat, Elevation, Elevation^2, Prop Burn Sev (PBS) 1yr, PBS 2-5, PBD 6-10, PBS 11-35, Stand Age, % Canopy Cover  
      logit(psi[i,k]) <- lpsi[k] + beta1[k] * lat[i] + beta2[k] * ele[i] + beta3[k] * pow(ele[i], 2) + beta4[k] * ppt[i] + beta5[k] * cbi1_5[i] + beta6[k] * cbi6_10[i] +
        beta7[k] * cbi11_35[i] + beta8[k] * cc_cfo[i] + beta9[k] * ch_cfo[i]
      
      # True latent state (Z-matrix)
      z[i,k] ~ dbern(psi[i,k])
      
      # Model Assement via Chi-squared GoF
      evalZ[i,k] <- psi[i,k]
      EZ[i,k] <- pow((z[i,k] - evalZ[i,k]), 2) / (evalZ[i,k] + 0.5)
        
      # Replicated data for new comparison
      z.new[i,k] ~ dbern(psi[i,k])
      EZ.new[i,k] <- pow((z.new[i,k] - evalZ[i,k]), 2) / (evalZ[i,k] + 0.5)
      
      }
    }
    
    ###################################################
    ## Observation submodel for replicate det/nondet ##
    ###################################################
    
    ## Observation Sub-model with nested indexing on sites
    
    # Significant NAs in the response for variable ARU start dates
    # Detection heterogeneity is built into this model
    # Detection Covs: Survey Hrs (sum across 6 day secondary sampling), Jdate (median per 6 day interval), Jdate^2
    for(k in 1:nspec){ #columns y 
      for(j in 1:N){ # rows y
        # Detection model on logit scale
        logit(p[j,k]) <- lp[k] + alpha1[k] * eff.hrs[j] + alpha2[k] * eff.jday[j] + alpha3[k] * pow(eff.jday[j], 2)
          
        # Latent state and detection
        # Site_id nested index
        mup[j,k] <- z[site_id[j], k] * p[j,k]
          
        # Observation Bernoulli draw from z*p
        y[j,k] ~ dbern(mup[j,k])
        
        # Model Assement via Chi-squared GoF
        eval[j,k] <- mup[j,k]
        E[j,k] <- pow((y[j,k] - eval[j,k]), 2) / (eval[j,k] + 0.5)
        
        # Replicated data for new comparison
        y.new[j,k] ~ dbern(mup[j,k])
        E.new[j,k] <- pow((y.new[j,k] - eval[j,k]), 2) / (eval[j,k] + 0.5)

      }
    }
    
    ########################################
    ## Derived estimates of the community ##
    ########################################
    
    # Species-specific # of occupied sites
    for(k in 1:nspec){
      Nocc.fs[k] <- sum(z[,k])
    }
    # Species Richness
    for(i in 1:nsite){
      Nsite[i] <- sum(z[i,])
    }
    
    fitZ <- sum(EZ[,])
    fitZ.new <- sum(EZ.new[,])
    fitY <- sum(E[,])
    fitY.new <- sum(E.new[,])
    
    }
    ", fill = TRUE)
sink()

## -------------------------------------------------------------
##
## End Section: JAGS Model Code
##
## -------------------------------------------------------------

## -------------------------------------------------------------
##
## Begin Section: Running JAGS Models
##
## -------------------------------------------------------------

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Initialization
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

nspec <- win.data.rag$nspec

## Initial values for model
## Zst 1 if even recorded, 0 otherwise
zst <- apply(y, c(1,3), function(x) max(x, na.rm = T))
colnames(zst) <- NULL
inits <- function() list(z = zst, 
                         lpsi = rep(0.4, nspec), 
                         beta1 = rep(0, nspec),
                         beta2 = rep(0, nspec),
                         beta3 = rep(0, nspec),
                         beta4 = rep(0, nspec),
                         beta5 = rep(0, nspec),
                         beta6 = rep(0, nspec),
                         beta7 = rep(0, nspec),
                         beta8 = rep(0, nspec),
                         beta9 = rep(0, nspec),
                         lp = rep(0.5, nspec),
                         alpha1 = rep(0, nspec),
                         alpha2 = rep(0, nspec),
                         alpha3 = rep(0, nspec),
                         alpha3 = rep(0, nspec))

## Parameters to monitor
params1 <- c("mu.lpsi",
             "sd.lpsi",
             "mu.beta1",
             "sd.beta1",
             "mu.beta2",
             "sd.beta2",
             "mu.beta3",
             "sd.beta3",
             "mu.beta4",
             "sd.beta4",
             "mu.beta5",
             "sd.beta5",
             "mu.beta6",
             "sd.beta6",
             "mu.beta7",
             "sd.beta7",
             "mu.beta8",
             "sd.beta8",
             "mu.beta9",
             "sd.beta9",
             "mu.lp",
             "sd.lp",
             "mu.alpha1",
             "sd.alpha1",
             "mu.alpha2",
             "sd.alpha2",
             "mu.alpha3",
             "sd.alpha3",
             #"mu.alpha4",
             #"sd.alpha4",
             "Ntotal",
             "Nsite",
             "fitZ",
             "fitZ.new",
             "fitY",
             "fitY.new")
# params2 <- c("mu.lpsi",
#              "sd.lpsi",
#              "mu.beta1",
#              "sd.beta1",
#              "mu.beta2",
#              "sd.beta2",
#              "mu.beta3",
#              "sd.beta3",
#              "mu.beta4",
#              "sd.beta4",
#              "mu.beta5",
#              "sd.beta5",
#              "mu.beta6",
#              "sd.beta6",
#              "mu.beta7",
#              "sd.beta7",
#              "mu.beta8",
#              "sd.beta8",
#              "mu.beta9",
#              "sd.beta9",
#              "lpsi",
#              "beta1",
#              "beta2",
#              "beta3",
#              "beta4",
#              "beta5",
#              "beta6",
#              "beta7",
#              "beta8",
#              "beta9",
#              "lp",
#              "alpha1",
#              "alpha2",
#              "alpha3",
#              "z")

## Settings for the MCMC
ni <- 10000
nt <- 10
nb <- 1000
nc <- 3

## to keep
# to_keep <- c("out2", "out3", "nsite", "nspec", "aru_meta", "sp.names", "win.data", 
#              "cbi1", "cbi2_5", "cbi6_10", "cbi11_35", "ele", "stage", "cc", "utmn", "eff.hrs", "eff.days")
# to_remove <- setdiff(ls(), to_keep)
# 
# rm(list = to_remove)

## Running the model
out2 <- jags(data = win.data.rag,
             inits = inits,
             params1,
             "Code/JAGS_Models/Sierra_MSOM_Covs_NA_95_BM.txt",
             n.chains = nc,
             n.thin = nt,
             n.iter = ni,
             n.burnin = nb,
             parallel = TRUE)

# out3 <- jags.basic(win.data.new, 
#                    inits, 
#                    params2, 
#                    "Sierra_MSOM_Covs_NA.txt", 
#                    n.chains = nc,
#                    n.thin = nt, 
#                    n.iter = ni, 
#                    n.burnin = nb, 
#                    parallel = TRUE)


## -------------------------------------------------------------
##
## End Section: JAGS Model Runs
##
## -------------------------------------------------------------

## -------------------------------------------------------------
##
## Begin Section: Prior Predictive Checks
##
## -------------------------------------------------------------

sd.lpsi <- runif(1000, 0, 4) 
hist(sd.lpsi, breaks = 50)
tau.lpsi <- sd.lpsi^-2
hist(tau.lpsi, breaks = 10)
mu.lpsi <- rnorm(1000, 0, 0.1)
hist(mu.lpsi, breaks = 50)
lpsi <- rnorm(1000, mu.lpsi, 1/tau.lpsi)
quantile(lpsi, probs = seq(0,1, by = .05))
hist(lpsi, breaks = 50)

# Detection submodel
# Number of species
K <- 91

# J number of rep
# Covariates
eff.hrs <- win.data.rag$eff.hrs
eff.jday <- win.data.rag$eff.jday
J <- length(eff.hrs)

# Species random intercept
sd.lp <- runif(K, 0, 2)
tau.lp <- sd.lp^-2
mu.lp <- rnorm(K, 0, 1/0.01) #convert precision to SD 

## Random intercept per species
lp <- rnorm(mu.lp, 1/tau.lp)

# Species random slope
# Alpha 1
sd.alpha1 <- runif(K, 0, 2)
tau.alpha1 <- sd.alpha1^-2
mu.alpha1 <- rnorm(K, 0, 1/0.01)
hist(mu.alpha1)

# Random slope
alpha1 <- rnorm(mu.alpha1, 1/tau.alpha1)

# Alpha 2
sd.alpha2 <- runif(K, 0, 2)
tau.alpha2 <- sd.alpha2^-2
mu.alpha2 <- rnorm(K, 0, 1/0.01)
hist(mu.alpha2)

# Random slope
alpha2 <- rnorm(mu.alpha2, 1/tau.alpha2)

# Alpha 3
sd.alpha3 <- runif(K, 0, 2)
tau.alpha3 <- sd.alpha3^-2
mu.alpha3 <- rnorm(K, 0, 1/0.01)
hist(mu.alpha3)

# Random slope
alpha3 <- rnorm(mu.alpha3, 1/tau.alpha3)

## Estimate prior distribution of detection probability from random draws
p <- matrix(NA, nrow = J, ncol = K)

for(k in 1:K){
  for(j in 1:J){
    p[j,k] <- plogis(lp[k] + alpha1[k] * eff.hrs[j] + alpha2[k] * eff.jday[j] + alpha3[k] * eff.jday[j]^2)
  }
}

## Prior estimates of detection probability per species 
hist(p)

## -------------------------------------------------------------
##
## Begin Section: Posterior Predictive Checks
##
## -------------------------------------------------------------
out2
summary(out2)
pp.check(out2, observed = 'fitY', simulated = 'fitY.new')
pp.check(out2, observed = 'fitZ', simulated = 'fitZ.new')



## *************************************************************
##
## Section Notes: Code below has not been updated!!!
##
## *************************************************************

## -------------------------------------------------------------
##
## Begin Section: JAGS Output Processing and Viz
##
## -------------------------------------------------------------

all3 <- as.matrix(out3)
pm <- apply(all3, 2, mean)
cri <- apply(all3, 2,function(x) quantile(x, prob = c(0.025, 0.0975)))
#save(out2, file = here("./Data/MSOMchainOut.Rdata"))
#load(file = here::here("./Data/MSOMchainOut.Rdata"))

## Z matrix
nsite <- nsite
nspec <- nspec
nsamp <- dim(all3)[1]
z <- array(NA, dim = c(nsite, nspec, nsamp))
Jacc <- array(NA, dim = c(nsite, nspec, nsamp))

for(j in 1:nsamp){
  cat(paste("\nMCMC sample", j, "\n"))
  z[,,j] <- all3[j, 1391:149174]
}

# Restrict computations to observed species
zobs <- z[,1:91,] # Species 1 to 145
# Compute Jaccard index for sites and for species
Jsite <- array(NA, dim = c(nsite, nsamp))
Jspec <- array(NA, dim = c(91, nsamp))
Jbetap <- array(NA, dim = c(nsite, nsamp))
# Choose reference site and species for Jaccard indices
ref.site <- 95 # Just choose first site
ref.species <- 1 # Acorn Woodpecker (check object 'obs.occ')
# Get posterior distributions for Jsite and Jspec (for references)
for(k in 1:nsamp){
  for(i in 1:nsite){ # Jaccard index for sites (in terms of shared species)
    Jsite[i,k] <- sum(zobs[ref.site,,k] * zobs[i,,k]) / (sum(zobs[ref.site,,k]) +
                                                           sum(zobs[i,,k]) - sum(zobs[ref.site,,k] * zobs[i,,k]))
  }
  for(i in 1:nspec){ # Jacc. index for species (in terms of shared sites)
    Jspec[i,k] <- sum(zobs[,ref.species,k] * zobs[,i,k]) / (sum(zobs[,ref.species,k]) +
                                                              sum(zobs[,i,k]) - sum(zobs[,ref.species,k] * zobs[,i,k]))
  }
}

for(k in 1:nsamp){
  test <- betapart::beta.pair(zobs[,,k], index.family = "sor")
  sim <- test$beta.sim
}

# NA's arise when a site has no species or a species no sites
# Get posterior means, standard deviations and 95% CRI
# Jaccard index for sites, compared to reference site 1
pm <- apply(Jsite, 1, mean, na.rm = TRUE) # Post. mean of Jsite wrt. site 1
psd <- apply(Jsite, 1, sd, na.rm = TRUE) # Post. sd of Jsite wrt. site 1
cri <- apply(Jsite, 1, function(x) quantile(x, prob = c(0.025, 0.975), na.rm =
                                              TRUE)) # CRI
cbind('post. mean' = pm, 'post. sd' = psd, '2.5%' = cri[1,], '97.5%' = cri[2,])

x <- 3 # size setting for plotting symbol

JZ <- data.frame(Y = aru_meta$Y, X = aru_meta$X,
                 Beta = pm, Beta.sd = psd, Beta.cril = cri[1,], Beta.criu = cri[2,])

JZsf <- JZ |> st_as_sf(coords = c("X", "Y"), crs = 4326)

ggplot(JZsf) + 
  geom_sf(aes(color = Beta)) + 
  scale_color_viridis_c(option = "H") + 
  theme_bw()

library(mapview)

mapView(JZsf, zcol = "Beta")


# Jaccard index for species, compared with a reference species
# (species 13, European Sparrowhawk)
pm <- apply(Jspec, 1, mean, na.rm = TRUE) # Post. mean of Jspec wrt. species 1
psd <- apply(Jspec, 1, sd, na.rm = TRUE) # Post. sd of Jspec wrt. species 1
cri <- apply(Jspec, 1, function(x) quantile(x, prob = c(0.025, 0.975), na.rm =
                                              TRUE)) # CRI
tmp <- cbind('post. mean' = pm, 'post. sd' = psd, '2.5%' = cri[1,], '97.5%' = cri[2,])
rownames(tmp) <- sp.names[[3]]
print(tmp) # print in systematic order
print(tmp[rev(order(tmp[,1])),]) # print in order of decreasing Jacc. values


psi.sample <- plogis(rnorm(10^6, mean = out2$mean$mu.lpsi, sd = out2$mean$sd.lpsi))
p.sample <- plogis(rnorm(10^6, mean = out2$mean$mu.lp, sd = out2$mean$sd.lp))
ppsi.dist <- data.frame(psi = psi.sample, p = p.sample)
psid <- ggplot(ppsi.dist) +
  geom_histogram(aes(x=psi), fill = "deepskyblue2", alpha = 0.8, color = "black") +
  geom_vline(xintercept = mean(psi.sample), color = "black", size = 1, linetype = "dashed") + 
  theme_bw() + 
  xlab(expression(psi)) + 
  ylab("Density")

detd <- ggplot(ppsi.dist) +
  geom_histogram(aes(x=p), fill = "darkseagreen", alpha = 0.8, color = "black") +
  geom_vline(xintercept = mean(p.sample), color = "black", size = 1, linetype = "dashed") + 
  theme_bw() + 
  xlab(expression(p)) + 
  ylab("Density")

library(patchwork)
psid + detd

ggsave(filename = here::here("./Figures/Preliminary/MSOMmeanpsip.jpg"),
       height = 6, width = 8, dpi = 600)


summary(psi.sample) ; summary(p.sample)

LaplacesDemon::invlogit(outTest$mean$beta0)
outTest$mean$beta1
outTest$mean$Nocc.fs

## Plotting some of the responses of the community
mean.ele <- mean(aru_meta$topo_elev)
sd.ele <- sd(aru_meta$topo_elev)
o.ele <- seq(320, 2900,,500) # Get covariate values for prediction

# Sage
mean.sage <- mean(aru_meta$standage_f3_mn)
sd.sage <- sd(aru_meta$standage_f3_mn)
o.sage <- seq(13, 252,,500)


o.dat <- seq(15, 120,,500)
o.dur <- seq(100, 420,,500)
ele.pred <- (o.ele - mean.ele) / sd.ele
sage.pred <- (o.sage - mean.sage) / sd.sage
dat.pred <- (o.dat - mean.date) / sd.date
dur.pred <- (o.dur - mean.dur) / sd.dur

str( tmp <- out2$sims.list ) # grab MCMC samples
nsamp <- length(tmp[[1]]) # number of mcmc samples
predC <- array(NA, dim = c(500, nsamp, 2)) # "C" for 'community mean'

for(i in 1:nsamp){
  predC[,i,1] <- plogis(tmp$mu.lpsi[i] + tmp$mu.beta2[i] * ele.pred +
                          tmp$mu.beta3[i] * ele.pred^2 )
  predC[,i,2] <- plogis(tmp$mu.lpsi[i] + tmp$mu.beta10[i] * sage.pred)
  
}

pmC <- apply(predC, c(1,3), mean)
criC <- apply(predC, c(1,3), function(x) quantile(x, prob = c(0.025, 0.975)))
plot(o.ele, pmC[,1], col = "blue", lwd = 3, type = 'l', lty = 1, frame = F,
     ylim = c(0, 0.2), xlab = "Elevation (m a.s.l)", ylab = "Community mean occupancy")
matlines(o.ele, t(criC[,,1]), col = "grey", lty = 1)
plot(o.sage, pmC[,2], col = "blue", lwd = 3, type = 'l', lty = 1, frame = F,
     ylim = c(0, 0.2), xlab = "Mean Stand Age (years)", ylab = "Community mean occupancy")
matlines(o.sage, t(criC[,,2]), col = "grey", lty = 1)

# Species-level responses
predS <- array(NA, dim = c(500, nspec, 1))
p.coef <- cbind(lp=pm[1292:1436], alpha1 = pm[1:145], alpha2 = pm[216:360])
psi.coef <- cbind(lpsi=pm[1507:1651], beta1 = pm[646:790], beta2 = pm[861:1005],
                  beta3 = pm[1076:1220])
for(i in 1:nspec){ # Loop over 145 observed species
  predS[,i,1] <- plogis(p.coef[i,1] + p.coef[i,2] * dat.pred +
                          p.coef[i,3] * dat.pred^2 ) # p ~ date
  predS[,i,2] <- plogis(p.coef[i,1] + p.coef[i,4] * dur.pred) # p ~ duration
  predS[,i,3] <- plogis(psi.coef[i,1] + psi.coef[i,2] * ele.pred +
                          psi.coef[i,3] * ele.pred^2 ) # psi ~ elevation
  predS[,i,4] <- plogis(psi.coef[i,1] + psi.coef[i,4] * for.pred) # psi ~ forest
}
# Plots for detection probability and survey date and duration (Fig. 11-24)
par(mfrow = c(1,2), cex.lab = 1.3, cex.axis = 1.3)
plot(o.dat, predS[,1,1], lwd = 3, type = 'l', lty = 1, frame = F,
     ylim = c(0, 1), xlab = "Survey date (1 = 1 April)",
     ylab = "Detection probability")
for(i in 2:145){
  lines(o.dat, predS[,i,1], col = i, lwd = 3)
}