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

## -------------------------------------------------------------

## -------------------------------------------------------------
##
## Begin Section: JAGS Models
##
## -------------------------------------------------------------

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Model 1
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## *************************************************************
##
## Section Notes: MSOM with modeled detection heterogeneity
## species-species random effects for intercepts and slopes
## and detection heterogeneity is modeled. Occupancy is
## is a species species random effect, but not modeled with 
## covariates.
##
## *************************************************************

sink("Sierra_MSOM.txt")
cat("
    model {
    
    ############
    ## Priors ##
    ############
    
    for(k in 1:nspec){
      # Psi
      psi[k] ~ dunif(0, 1) #uniform prior for species occ
      
      # p
      beta0[k] ~ dnorm(mu.beta0, tau.beta0)
      beta1[k] ~ dnorm(mu.beta1, tau.beta1)
      beta2[k] ~ dnorm(mu.beta2, tau.beta2)
    }
    
    #################
    ## Hyperpriors ##
    #################
    
    # Occ hypers
    
    # Detection hypers
    mu.beta0 ~ dnorm(0, 0.1)
    tau.beta0 <- pow(sd.beta0, -2)
    sd.beta0 ~ dunif(0, 2)
    mu.beta1 ~ dnorm(0, 0.1)
    tau.beta1 <- pow(sd.beta1, -2)
    sd.beta1 ~ dunif(0, 2)
    mu.beta2 ~ dnorm(0, 0.1)
    tau.beta2 <- pow(sd.beta2, -2)
    sd.beta2 ~ dunif(0, 2)
    
    #################################################
    ## Ecological model for the latent process (z) ##
    #################################################
    # Add covariates for occupancy
    for(k in 1:nspec){ #species loop
      for(i in 1:nsite){ #site loop
      ## insert occupancy covariates in the future
      ## for now latent process only
      z[i,k] ~ dbern(psi[k]) #true latent state estimation from psi
      }
    }
    
    ###################################################
    ## Observation submodel for replicate det/nondet ##
    ###################################################
    
    # Significant NAs in the response for variable ARU start dates
    # Detection heterogeneity is built into this model
    for(k in 1:nspec){
      for(i in 1:nsite){
        for(j in 1:nrep){
          logit(p[i,j,k]) <- beta0[k] + beta1[k] * eff.days[i,j] + beta2[k] * eff.hrs[i,j] #+ beta3[k] * eff.jday[i,j] + beta4[k] * pow(eff.jday[i,j], 2)
          mup[i,j,k] <- z[i,k] * p[i,j,k]
          y[i,j,k] ~ dbern(mup[i,j,k])
        }
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
    
    }
    ", fill = TRUE)
sink()

## Initial values for model
zst <- apply(y, c(1,3), max)
colnames(zst) <- NULL
inits <- function() list(z = zst, 
                         psi = rep(0.4, nspec), 
                         beta0 = rep(0.4, nspec),
                         beta1 = rep(0, nspec),
                         beta2 = rep(0, nspec))

## Parameters to monitor
params <- c("psi", "beta0", "beta1", "beta2", "Nocc.fs", "Nsite", "z")

## Settings for the MCMC
ni <- 2500
nt <- 2
nb <- 500
nc <- 3

## Running the model
outTest <- jags(data = win.data,
                inits = inits,
                params,
                "Sierra_MSOM.txt",
                n.chains = nc,
                n.thin = nt,
                n.iter = ni,
                n.burnin = nb,
                parallel = T)

LaplacesDemon::invlogit(outTest$mean$beta0)
outTest$mean$beta1
outTest$mean$Nocc.fs

MCMCsummary(outTest, params = 'beta2', round = 2)
mu.beta0 <- rnorm(2500, 0, 100)
sd.beta0 <- runif(2500, 0, 2)
tau.beta0 <- I(sd.beta0^-2)
PR <- rnorm(mu.beta0, sd.beta0)

MCMCvis::MCMCtrace(outTest, params = 'beta1[1]', 
                   ISB = FALSE, exact = TRUE, priors = PR, ind = TRUE,
                   Rhat = TRUE, n.eff = TRUE, pdf = FALSE)
summary(outTest)

## Generate the Z matrix
allOut <- as.matrix(outTest)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Model 2 
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## *************************************************************
##
## Section Notes: Ideally we can represent species detection
## and (maybe?) occupancy RE as arising from a multivariate
## normal conditional on phylogenetic relatedness. This would 
## involve building in the structure for linking priors of species
## intercepts for p and psi from the tree itself.
##
## *************************************************************

sink("Sierra_MSOM_Phylo.txt")
cat("
    model {
    
    ############
    ## Priors ##
    ############
    
    for(k in 1:nspec){
      # Psi (phylo)
      psi[k] <- eta[k,1] #uniform prior for species occ
      p[k] <- eta[k,2]
      eta[k, 1:2] ~ dmnorm(mu.eta[], Omega[,]) #mu.eta would be mean psi and p & omega is the variance covariance (based on phylogeny)
      
      # p
      beta0[k] ~ dnorm(mu.beta0, tau.beta0)
      beta1[k] ~ dnorm(mu.beta1, tau.beta1)
      beta2[k] ~ dnorm(mu.beta2, tau.beta2)
    }
    
    #
    
    #################
    ## Hyperpriors ##
    #################
    
    for(v in 1:2){
      mu.eta[v] <- log(probs[v] / (1-probs[v])) #logit scale mean det and occ
      probs[v] ~ dunif(0,1)
      
    }
    
    mu.beta0 ~ dnorm(0, 0.1)
    tau.beta0 <- pow(sd.beta0, -2)
    sd.beta0 ~ dunif(0, 2)
    mu.beta1 ~ dnorm(0, 0.1)
    tau.beta1 <- pow(sd.beta1, -2)
    sd.beta1 ~ dunif(0, 2)
    mu.beta2 ~ dnorm(0, 0.1)
    tau.beta2 <- pow(sd.beta2, -2)
    sd.beta2 ~ dunif(0, 2)
    
    # Prior for variance-covariance matrix
    Omega[1:nspec, 1:nspec] ~ dwish(PhyCov[,], df)
    Sigma[1:nspec, 1:nspec] <- inverse(Omega[,])
    
    #################################################
    ## Ecological model for the latent process (z) ##
    #################################################
    # Add covariates for occupancy
    for(k in 1:nspec){ #species loop
      for(i in 1:nsite){ #site loop
      ## insert occupancy covariates in the future
      ## for now latent process only
      z[i,k] ~ dbern(psi[k]) #true latent state estimation from psi
      }
    }
    
    ###################################################
    ## Observation submodel for replicate det/nondet ##
    ###################################################
    
    # Significant NAs in the response for variable ARU start dates
    # Detection heterogeneity is built into this model
    for(k in 1:nspec){
      for(i in 1:nsite){
        for(j in 1:nrep){
          logit(p[i,j,k]) <- beta0[k] + beta1[k] * eff.days[i,j] + beta2[k] * eff.hrs[i,j] #+ beta3[k] * eff.jday[i,j] + beta4[k] * pow(eff.jday[i,j], 2)
          mup[i,j,k] <- z[i,k] * p[i,j,k]
          y[i,j,k] ~ dbern(mup[i,j,k])
        }
      }
    }
    
    ########################################
    ## Derived estimates of the community ##
    ########################################
    
    rho <- Sigma[1,2] / sqrt(Sigma[1,1]*Sigma[2,2])
    
    # Species-specific # of occupied sites
    for(k in 1:nspec){
      Nocc.fs[k] <- sum(z[,k])
    }
    # Species Richness
    for(i in 1:nsite){
      Nsite[i] <- sum(z[i,]) 
    }
    
    }
    ", fill = TRUE)
sink()

## Initial values for model
zst <- apply(y, c(1,3), max)
colnames(zst) <- NULL
inits <- function() list(z = zst,
                         Omega = diag(nspec),
                         eta = matrix(0, nrow = nspec, ncol = 2), 
                         beta0 = rep(0.4, nspec),
                         beta1 = rep(0, nspec),
                         beta2 = rep(0, nspec))

## Parameters to monitor
params <- c("psi", "beta0", "beta1", "beta2", "Nocc.fs", "Nsite", "Sigma", "rho")

## Settings for the MCMC
ni <- 100
nt <- 2
nb <- 50
nc <- 3

## Windata
win.data$PhyCov <- phy.vcv

## Running the model
outTest <- jags(data = win.data,
                inits = inits,
                params,
                "Sierra_MSOM_Phylo.txt",
                n.chains = nc,
                n.thin = nt,
                n.iter = ni,
                n.burnin = nb,
                parallel = T)

LaplacesDemon::invlogit(outTest$mean$beta0)
outTest$mean$beta1
outTest$mean$Nocc.fs

MCMCsummary(outTest, params = 'beta2', round = 2)
mu.beta0 <- rnorm(2500, 0, 100)
sd.beta0 <- runif(2500, 0, 2)
tau.beta0 <- I(sd.beta0^-2)
PR <- rnorm(mu.beta0, sd.beta0)

MCMCvis::MCMCtrace(outTest, params = 'beta1[1]', 
                   ISB = FALSE, exact = TRUE, priors = PR, ind = TRUE,
                   Rhat = TRUE, n.eff = TRUE, pdf = FALSE)
summary(outTest)

## Generate the Z matrix
allOut <- as.matrix(outTest)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Model 3 Detection and Occupancy Covariates
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sink("Sierra_MSOM_Covs.txt")
cat("
    model {
    
    # Book-keeping:
    # i = site; j = rep; k = species
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
      beta10[k] ~ dnorm(mu.beta10, tau.beta10)
      beta11[k] ~ dnorm(mu.beta11, tau.beta11)
      
      # p
      # Species random p intercept
      lp[k] ~ dnorm(mu.lp, tau.lp)
      
      # Species random p slope
      alpha1[k] ~ dnorm(mu.alpha1, tau.alpha1)
      alpha2[k] ~ dnorm(mu.alpha2, tau.alpha2)
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
    # Precip
    mu.beta4 ~ dnorm(0,0.1)
    tau.beta4 <- pow(sd.beta4, -2)
    sd.beta4 ~ dunif(0,4)
    # Max T
    mu.beta5 ~ dnorm(0,0.1)
    tau.beta5 <- pow(sd.beta5, -2)
    sd.beta5 ~ dunif(0,4)
    # Fire Sev 1yr
    mu.beta6 ~ dnorm(0,0.1)
    tau.beta6 <- pow(sd.beta6, -2)
    sd.beta6 ~ dunif(0,4)
    # Fire Sev 2-5 years
    mu.beta7 ~ dnorm(0,0.1)
    tau.beta7 <- pow(sd.beta7, -2)
    sd.beta7 ~ dunif(0,4)
    # Fire Sev 6-10 yrs
    mu.beta8 ~ dnorm(0,0.1)
    tau.beta8 <- pow(sd.beta8, -2)
    sd.beta8 ~ dunif(0,4)
    # Fire Sev 11-35 yrs
    mu.beta9 ~ dnorm(0,0.1)
    tau.beta9 <- pow(sd.beta9, -2)
    sd.beta9 ~ dunif(0,4)
    # Stand Age
    mu.beta10 ~ dnorm(0,0.1)
    tau.beta10 <- pow(sd.beta10, -2)
    sd.beta10 ~ dunif(0,4)
    # Canopy Cover
    mu.beta11 ~ dnorm(0,0.1)
    tau.beta11 <- pow(sd.beta11, -2)
    sd.beta11 ~ dunif(0,4)
    
    # Detection hypers
    # Intercept
    mu.lp ~ dnorm(0, 0.1)
    tau.lp <- pow(sd.lp, -2)
    sd.lp ~ dunif(0, 2)
    # Number of days sampled (per secondary sample)
    mu.alpha1 ~ dnorm(0, 0.1)
    tau.alpha1 <- pow(sd.alpha1, -2)
    sd.alpha1 ~ dunif(0, 2)
    # Number of hours sampled (per secondary sample)
    mu.alpha2 ~ dnorm(0, 0.1)
    tau.alpha2 <- pow(sd.alpha2, -2)
    sd.alpha2 ~ dunif(0, 2)
    
    #################################################
    ## Ecological model for the latent process (z) ##
    #################################################
    # Add covariates for occupancy
    for(k in 1:nspec){ #species loop
      for(i in 1:nsite){ #site loop
      # Occupancy model w/ covs
      logit(psi[i,k]) <- lpsi[k] + beta1[k] * utmn[i] + beta2[k] * ele[i] + beta3[k] * pow(ele[i], 2) + 
        beta4[k] * ppt[i] + beta5[k] * tmx[i] + beta6[k] * cbi1[i] + beta7[k] * cbi2_5[i] + beta8[k] * cbi6_10[i] +
        beta9[k] * cbi11_35[i] + beta10[k] * stage[i] + beta11[k] * cc[i]
      
      # True latent state (Z-matrix)
      z[i,k] ~ dbern(psi[i,k]) 
      }
    }
    
    ###################################################
    ## Observation submodel for replicate det/nondet ##
    ###################################################
    
    # Significant NAs in the response for variable ARU start dates
    # Detection heterogeneity is built into this model
    for(k in 1:nspec){
      for(i in 1:nsite){
        for(j in 1:nrep){
          # Detection model on logit scale
          logit(p[i,j,k]) <- lp[k] + alpha1[k] * eff.days[i,j] + alpha2[k] * eff.hrs[i,j] #+ beta3[k] * eff.jday[i,j] + beta4[k] * pow(eff.jday[i,j], 2)
          
          # Latent state and detection 
          mup[i,j,k] <- z[i,k] * p[i,j,k]
          
          # Observation Bernoulli draw from z*p
          y[i,j,k] ~ dbern(mup[i,j,k])
        }
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
    
    }
    ", fill = TRUE)
sink()

## Initial values for model
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
                         beta10 = rep(0, nspec),
                         beta11 = rep(0, nspec),
                         lp = rep(0.5, nspec),
                         alpha1 = rep(0, nspec),
                         alpha2 = rep(0, nspec))

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
             "mu.beta10",
             "sd.beta10",
             "mu.beta11",
             "sd.beta11",
             "mu.lp",
             "sd.lp",
             "mu.alpha1",
             "sd.alpha1",
             "mu.alpha2",
             "sd.alpha2",
             "Ntotal",
             "Nsite")
params2 <- c("mu.lpsi",
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
             "mu.beta10",
             "sd.beta10",
             "mu.beta11",
             "sd.beta11",
             "lpsi",
             "beta1",
             "beta2",
             "beta3",
             "beta4",
             "beta5",
             "beta6",
             "beta7",
             "beta8",
             "beta9",
             "beta10",
             "beta11",
             "lp",
             "alpha1",
             "alpha2",
             "z")

## Settings for the MCMC
ni <- 1000
nt <- 2
nb <- 200
nc <- 3

## to keep
to_keep <- c("out2", "out3", "nsite", "nspec", "aru_meta", "sp.names", "win.data", 
             "cbi1", "cbi2_5", "cbi6_10", "cbi11_35", "ele", "stage", "cc", "tmx", "ppt", "utmn", "eff.hrs", "eff.days")
to_remove <- setdiff(ls(), to_keep)

rm(list = to_remove)

## Running the model
out2 <- jags(data = win.data,
             inits = inits,
             params1,
             "Sierra_MSOM_Covs.txt",
             n.chains = nc,
             n.thin = nt,
             n.iter = ni,
             n.burnin = nb,
             parallel = T)

out3 <- jags.basic(win.data, 
                   inits, 
                   params2, 
                   "Sierra_MSOM_Covs.txt", 
                   n.chains = nc,
                   n.thin = nt, 
                   n.iter = ni, 
                   n.burnin = nb, 
                   parallel = TRUE)

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