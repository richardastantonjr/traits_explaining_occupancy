## Fit a models to 48 species with main effects of shrub cover and 3 land uses, without interactions
##------------------------------------------------------------------------------

## libraries
library(jagsUI)
##------------------------------------------------------------------------------

## load survey data, sampling covariates and site covariates
load("./GenerateAbunHists.RData")
rm(list = grep("^AbunH", ls(), value = TRUE, invert = TRUE)) ## keep only AbunHists
samplingCovs <- read.csv("./data/SamplingCovs.csv")
siteCovs <- read.csv("./data/siteCovs_numeric_grids.csv")

## Modify site and sampling covariates to work with JAGS/BUGS
#site-level covariates CANNOT have NAs
point <- as.numeric(factor(c(1:nrow(siteCovs))))
siteCovs$LandUse <-  relevel(siteCovs$LandUse, ref = "Protected")
land <- as.numeric(factor(siteCovs$LandUse))
nland <- length(levels(siteCovs$LandUse))
shrub <- scale (siteCovs$ShrubMean)
grid = as.numeric(factor(siteCovs$grid))

#sampling-level covariates CANNOT have NAs
date = scale(samplingCovs[,c(2:5)])
date[is.na(date)] <- 0      ## fill missing values with mean imputation
time = scale(samplingCovs[,c(12:15)])
time[is.na(time)] <- 0
timeSq <- time*time
wind = samplingCovs[,c(17:20)]
temp = scale(samplingCovs[,c(22:25)]) 
cloud = scale(samplingCovs[,c(27:30)])

##--------------------------------------------------------------------------------------------------------
## Next, modify data to suit the models with interaction terms, i.e. drop sites in sugar estate where there is negligible shrub cover
siteCovs <- siteCovs[which(siteCovs$LandUse != "SugarEstate"),]
siteCovs <- droplevels(siteCovs)
SamplingCovs <- samplingCovs[which(siteCovs$LandUse != "SugarEstate"),]
## remove sites in sugar cane and the 5th visits; nearly all 5th visits are NA
AbunHists.subset <- AbunHists.subset[which(siteCovs$LandUse != "SugarEstate"), 1:4, ]
##------------------------------------------------------------------------------------
## Redefine in terms of subsetted data
point<-as.numeric(factor(c(1:nrow(siteCovs))))
land<-as.numeric(factor(siteCovs$LandUse))
nland<-length(levels(siteCovs$LandUse))
shrub<-scale(siteCovs$ShrubMean)
grid = as.numeric(factor(siteCovs$grid))
date = scale(samplingCovs[,c(2:5)])
date[is.na(date)] <- 0 
time = scale(samplingCovs[,c(12:15)])
time[is.na(time)] <- 0
timeSq <- time*time
wind = samplingCovs[,c(17:20)]
temp = scale(samplingCovs[,c(22:25)]) 
cloud = scale(samplingCovs[,c(27:30)])

## create an empty vector to fill with each species' WAIC
occ.WAIC <- rep(NA, 48)
##---------------------------------------------------------------------------------

##---------------------------------------------------------------------------------------
## fit Occupancy models for each spp.
##---------------------------------------------------------------------------------------
##---------------------------------------------------------------------------------------
# 11.6.1 Simplest community occupancy model: n-fold single species 
#        occupancy model with species treated as fixed effects
# ---------------------------------------------------------------------------------------

##---------------------------------------------------------------------------------------
## (psi~ shrub + land)(p~ date + time + timeSq)
# Bundle and summarize data set
for (i in c(1:48)){         
  detect.i <- AbunHists.subset[,,i]
  y <- as.matrix(detect.i)
  # Bundle data
  jags.data <- list(y = y, 
                    R = length(point), 
                    J = ncol(y),
                    ngrid = length(unique(siteCovs$grid)),
                    grid = grid,
                    date = date,
                    timeSq = timeSq,
                    time = time,
                    shrub = shrub,
                    land = land)
  
  # Initial values
  zst <- apply(cbind(y, rep(0, nrow(y))), 1, max, na.rm = TRUE)  # use max of each row as starting value for z (0 / 1)
  
  inits <- function(){list(z = zst,
                           alpha1 = runif(1),
                           alpha2 = runif(1),
                           alpha3 = runif(1),
                           alpha4 = runif(1),
                           beta = runif(1),
                           beta1  = runif(1),
                           beta.l = c(NA, rnorm(2, mean = mean(y, na.rm = T), 3)), ## Cannot supply initial value for 1st land use [fixed at 0]
                           grid.effect = rnorm(length(unique(grid)), mean = mean(y, na.rm = T)), sd.grid = runif(1,0,5)
  ) }
  
  # Parameters monitored
  params <- c("pd", "log_pd") ## each species pd and lpd are ~ 1.1 gB
  #params <- c("beta", "beta1", "beta.l") ## run this if you want effects; 
  
  # MCMC settings
  ni <- 180000
  nt <- 50
  nb <- 15000 ##
  nc <- 4
  
  # Call JAGS from R 
  out.rep <- jagsUI(jags.data, inits, params, "./code/jagsmodel_3_land_uses.R", n.thin = nt,
                    n.chains = nc, n.burnin = nb, n.adapt = 5000, 
                    n.iter = ni, parallel = TRUE)
} 

## Save results
#saveRDS(occ.WAIC, "./outputs/WAIC_interaction_models_updates") ## Just WAIC for the models
                                                               
save.image("./Occupancy_analysis_Swazi_birds_Bayesian_quadTime_simplified_20180505.RData")
