
## Plot community response to shrub cover in different land uses, according to three models:
## 1) A model with a shrub cover by land use interaction and 3 land uses (sugar cane plantation was excluded)
## 2) A model with main effects only, and 3 land uses
## 3) A main effects only model, with 4 land uses (sugar can included)

# libraries ----------------------------------------------------------------------------
library(tidyverse)
library(reshape2) 
library(MCMCvis)

# data -----------------------------------------------------------------------------
load("./outputs/Occupancy_analysis_Swazi_birds_Bayesian_quadTime_WAIC_20180528.RData")
#load("./Occupancy_analysis_Swazi_birds_Bayesian_int_models_quadTime_fixed.RData")
#load("./Occupancy_analysis_Swazi_birds_Bayesian_quadTime_WAIC_20180528.RData")
# processing  -----------------------------------------------------------------------
## Put output from jagsUI outputs in a list- each corresponds to a single species
Outputs <- list(IntOutput1,IntOutput2,IntOutput3,IntOutput4,IntOutput5,
                IntOutput6,IntOutput7,IntOutput8,
                IntOutput9,IntOutput10,IntOutput11,IntOutput12,
                IntOutput13,IntOutput14,IntOutput15,IntOutput16,
                IntOutput17,IntOutput18,IntOutput19,IntOutput20,
                IntOutput21,IntOutput22,IntOutput23,IntOutput24,
                IntOutput25,IntOutput26,IntOutput27,IntOutput28,
                IntOutput29,IntOutput30,IntOutput31,IntOutput32,
                IntOutput33,IntOutput34,IntOutput35,IntOutput36,
                IntOutput37,IntOutput38,IntOutput39,IntOutput40,
                IntOutput41,IntOutput42,IntOutput43,IntOutput44,
                IntOutput45,IntOutput46,IntOutput47,IntOutput48)

## The objects Output1... Output48 are Large jagsUI; Only the posterior MCMC chains are needed.
Outputs <- ((map(Outputs, MCMCchains))) ## MCMCchains extracts the chains, map applies to each item in the list 

## take a random draw of ~2000 iterations from each posterior to get equal sample sizes
to_keep <- matrix(nrow = 2000, ncol = 48)
for(i in 1:48){
  to_keep[,i] <- sample(1:dim(Outputs[[i]])[1], size = 2000)                  
  Outputs[[i]] <- subset(Outputs[[i]][c(to_keep[,i]), c(1:9)])   
}                  ## keep row numbers that are in column i of to_keep, for params of interest

## Assemble outputs from each species in a dataframe in long format
Outputs_long <- as.data.frame (do.call(rbind, Outputs))

## create a column for species name for downstream processing
sppName <- vector(mode = "character", length = dim(Outputs[[1]])[1])  
sppName <- rep(dimnames(AbunHists.subset)[[3]], dim(Outputs[[1]])[1])  ## TraitData$Species
sppName <- sort(sppName)
Outputs_long$Species <- sppName 
Outputs_long$X <- rep(1:dim(Outputs[[1]])[1], 48)   ## add a row number for dcast to convert long to wide format

# subset by betas
shrub.post <- subset(Outputs_long, select = c("X", "Species","beta1"))
homestead.post <- subset(Outputs_long, select = c("X", "Species","beta.l[3]"))
pasture.post <- subset(Outputs_long, select = c("X", "Species","beta.l[2]"))  
protected.post <- subset(Outputs_long, select = c("X", "Species","beta")) ## intercept for this parameterization
pas_by_shrub.post <- subset(Outputs_long, select = c("X", "Species","beta.l.int[2]"))  
home_by_shrub.post <- subset(Outputs_long, select = c("X", "Species","beta.l.int[3]"))  

# reformat to make species columns and posterior samples rows with dcast
# then remove the 'X' column
shrub.post <- dcast(shrub.post, X~Species)
shrub.post <- shrub.post[, 2:49]
homestead.post <- dcast(homestead.post, X~Species)
homestead.post <- homestead.post[, 2:49]
pasture.post <- dcast(pasture.post, X~Species)
pasture.post <- pasture.post[, 2:49]
intercept.post <- dcast(protected.post, X~Species)
intercept.post <- intercept.post[, 2:49]
pas_by_shrub.post <- dcast(pas_by_shrub.post, X~Species)
pas_by_shrub.post <- pas_by_shrub.post[, 2:49]
home_by_shrub.post <- dcast(home_by_shrub.post, X~Species)
home_by_shrub.post <- home_by_shrub.post[, 2:49]

intercept <- rowMeans(intercept.post)
shrub <- rowMeans(shrub.post) 
pasture <- rowMeans(pasture.post)
homestead <- rowMeans(homestead.post)
pas_x_shrub <-rowMeans(pas_by_shrub.post)
home_x_shrub <- rowMeans(home_by_shrub.post)

o.shrub <- seq(0,100,,334)  ## get covariate values for prediction
shrub.pred <- (o.shrub - mean(siteCovs$ShrubMean)) / sd(siteCovs$ShrubMean)

predC_home <- array(NA, dim = c(334, 2000, 1))
predC_pasture <- array(NA, dim = c(334, 2000, 1))
predC_not_home <- array(NA, dim = c(334, 2000, 1))

for (i in 1:2000) {
    predC_home[, i, 1] <- plogis(intercept[i] + shrub[i] * shrub.pred + homestead[i] +
                              home_x_shrub[i] * shrub.pred)
    predC_pasture[, i, 1] <- plogis(intercept[i] + shrub[i] * shrub.pred + pasture[i] + 
                                   pas_x_shrub[i] * shrub.pred)
    predC_not_home[, i, 1] <- plogis(intercept[i] + shrub[i] * shrub.pred)
}

pmC_home <- apply(predC_home, 1, mean)
criC_home <- apply(predC_home, 1, function(x) quantile(x, prob = c(0.025, 0.975)))
pmC_pasture <- apply(predC_pasture, 1, mean)
criC_pasture <- apply(predC_pasture, 1, function(x) quantile(x, prob = c(0.025, 0.975)))
pmC_not_home <- apply(predC_not_home, 1, mean)
criC_not_home <- apply(predC_not_home, 1, function(x) quantile(x, prob = c(0.025, 0.975)))

# plotting --------------------------------------------
par(mfrow = c(3,1))
plot(o.shrub[1:119], pmC_home[1:119], col = "turquoise1", lwd = 3, type = 'l', lty = 1, frame = FALSE, ylim = c(0, 1),
     ylab = "Mean occupancy", las = 1, xlim = c(0, 100),
     xlab = "Shrub cover, %", main = "3 land uses, with land use by shrub cover interaction")
## matlines(o.shrub[1:119], criC_home[1,1:119], col = "turquoise1", lty = 1)
## matlines(o.shrub[1:119], criC_home[2,1:119], col = "turquoise1", lty = 1)
legend("bottomright", legend = c("protected areas","pastures", "homesteads", "plantation"),
       col = c("blue","green","turquoise1","goldenrod"), bty = "n", pch = 19, cex = 1.3
       )

lines(o.shrub, pmC_pasture, col = "green", lwd = 3, type = 'l', lty = 2)
## matlines(o.shrub, criC_pasture[1,], col = "green", lty = 2)
## matlines(o.shrub, criC_pasture[2,], col = "green", lty = 2)

lines(o.shrub, pmC_not_home, col = "blue", lwd = 3, type = 'l', lty = 3)
## matlines(o.shrub, criC_not_home[1,], col = "blue", lty = 3)
## matlines(o.shrub, criC_not_home[2,], col = "blue", lty = 3)

##----------------------
## Get species-specific predictions of occupancy at mean shrub cover to compare to model below
shrub_mean_spp_wise <- apply(shrub.post, 2, mean)
intercept_mean_spp_wise <- apply(intercept.post, 2, mean)
pasture_mean_spp_wise <- apply(pasture.post, 2, mean)
homestead_mean_spp_wise <- apply(homestead.post, 2, mean)
pas_by_shrub_spp_wise <- apply(pas_by_shrub.post, 2, mean)
home_by_shrub_spp_wise <- apply(home_by_shrub.post, 2, mean)

pred_home_int <- matrix(NA, 48, 1)
pred_pasture_int <- matrix(NA, 48, 1)

for (i in 1:48) {
  pred_home_int[i,] <- plogis(intercept_mean_spp_wise[i] + shrub_mean_spp_wise[i] * mean(shrub.pred) + 
                                homestead_mean_spp_wise[i] + home_by_shrub_spp_wise[i])
  pred_pasture_int[i,] <- plogis(intercept_mean_spp_wise[i] + shrub_mean_spp_wise[i] * mean(shrub.pred) + 
                                   pasture_mean_spp_wise[i] + pas_by_shrub_spp_wise[i])
}
##----------------------

## repeat for reduced model w/o interaction terms
load("./Occupancy_analysis_Swazi_birds_Bayesian_quadTime_simplified_20180505.RData")

## Put output from jagsUI in a list- each corresponds to a single species
Outputs <- list(Output1,Output2,Output3,Output4,Output5,
                Output6,Output7,Output8,
                Output9,Output10,Output11,Output12,
                Output13,Output14,Output15,Output16,
                Output17,Output18,Output19,Output20,
                Output21,Output22,Output23,Output24,
                Output25,Output26,Output27,Output28,
                Output29,Output30,Output31,Output32,
                Output33,Output34,Output35,Output36,
                Output37,Output38,Output39,Output40,
                Output41,Output42,Output43,Output44,
                Output45,Output46,Output47,Output48)

Outputs <- ((map(Outputs, MCMCchains)))    

## take a random draw of ~2000 iterations from each posterior to get equal sample sizes
to_keep <- matrix(nrow = 2000, ncol = 48)
for(i in 1:48){
  to_keep[,i] <- sample(1:dim(Outputs[[i]])[1], size = 2000)                  
  Outputs[[i]] <- subset(Outputs[[i]][c(to_keep[,i]), c(4, 5, 7, 8, 9)])  
}                  ## keep row numbers that are in column i of to_keep, for params of interest

## Assemble outputs from each species in a dataframe in long format
Outputs_long <- as.data.frame (do.call(rbind, Outputs))

## create a column for species name for downstream processing
sppName <- vector(mode = "character", length = dim(Outputs[[1]])[1])  
sppName <- rep(dimnames(AbunHists.subset)[[3]], dim(Outputs[[1]])[1])  ## TraitData$Species
sppName <- sort(sppName)
Outputs_long$Species <- sppName 
Outputs_long$X <- rep(1:dim(Outputs[[1]])[1], 48)   ## add a row number for dcast to convert long to wide format

#subset by effect size or beta, keeping posterior draw, species beta:
shrub.post<-subset(Outputs_long, select=c("X", "Species","beta1"))
homestead.post<-subset(Outputs_long, select=c("X", "Species","beta.l[3]"))
pasture.post<-subset(Outputs_long, select=c("X", "Species","beta.l[2]"))  
intercept.post<-subset(Outputs_long, select=c("X", "Species","beta")) ## intercept for this parameterization

## reformat to make species columns and posterior samples rows with dcast
## then drop the 'X' column
shrub.post <- dcast(shrub.post, X~Species)
shrub.post <- shrub.post[,2:49]
homestead.post <- dcast(homestead.post, X~Species)
homestead.post <- homestead.post[,2:49]
pasture.post <- dcast(pasture.post, X~Species)
pasture.post <- pasture.post[,2:49]
intercept.post <- dcast(intercept.post, X~Species)
intercept.post <- intercept.post[,2:49]

intercept <- rowMeans(intercept.post)
shrub <- rowMeans(shrub.post) 
pasture <- rowMeans(pasture.post)
homestead <- rowMeans(homestead.post)

predC_home <- array(NA, dim = c(334, 2000, 1))
predC_pasture <- array(NA, dim = c(334, 2000, 1))
predC_not_home <- array(NA, dim = c(334, 2000, 1))

for (i in 1:2000) {
  predC_home[, i, 1] <- plogis(intercept[i] + shrub[i] * shrub.pred + homestead[i])
  predC_pasture[, i, 1] <- plogis(intercept[i] + shrub[i] * shrub.pred + pasture[i])
  predC_not_home[, i, 1] <- plogis(intercept[i] + shrub[i] * shrub.pred)
}

pmC_home <- apply(predC_home, 1, mean)
criC_home <- apply(predC_home, 1, function(x) quantile(x, prob = c(0.025, 0.975)))
pmC_pasture <- apply(predC_pasture, 1, mean)
criC_pasture <- apply(predC_pasture, 1, function(x) quantile(x, prob = c(0.025, 0.975)))
pmC_not_home <- apply(predC_not_home, 1, mean)
criC_not_home <- apply(predC_not_home, 1, function(x) quantile(x, prob = c(0.025, 0.975)))

# plotting --------------------------------------------
plot(o.shrub[1:119], pmC_home[1:119], col = "turquoise1", lwd = 3, type = 'l', lty = 1, frame = FALSE, ylim = c(0, 1),
     ylab = "Mean occupancy", las = 1, xlim = c(0, 100),
     xlab = "Shrub cover, %", main = "3 land uses, without land use by shrub cover interaction")
 ##matlines(o.shrub[1:119], criC_home[1,1:119], col = "turquoise1", lty = 1)
 ##matlines(o.shrub[1:119], criC_home[2,1:119], col = "turquoise1", lty = 1)

lines(o.shrub, pmC_pasture, col = "green", lwd = 3, type = 'l', lty = 2)
 ##matlines(o.shrub, criC_pasture[1,], col = "green", lty = 2)
 ##matlines(o.shrub, criC_pasture[2,], col = "green", lty = 2)

lines(o.shrub, pmC_not_home, col = "blue", lwd = 3, type = 'l', lty = 3)
 ##matlines(o.shrub, criC_not_home[1,], col = "blue", lty = 3)
 ##matlines(o.shrub, criC_not_home[2,], col = "blue", lty = 3)

##---------------
##---------------------------------------------------------------
## Get predicted ocupancy at mean shrub cover for each species to compare to the model above
shrub_mean_spp_wise <- apply(shrub.post, 2, mean)
intercept_mean_spp_wise <- apply(intercept.post, 2, mean)
pasture_mean_spp_wise <- apply(pasture.post, 2, mean)
homestead_mean_spp_wise <- apply(homestead.post, 2, mean)

pred_home <- matrix(NA, 48, 1)
pred_pasture <- matrix(NA, 48, 1)
pred_not_home <- matrix(NA, 48, 1)

for (i in 1:48) {
  pred_home[i,] <- plogis(intercept_mean_spp_wise[i] + shrub_mean_spp_wise[i] * mean(shrub.pred) + 
                            homestead_mean_spp_wise[i])
  pred_pasture[i,] <- plogis(intercept_mean_spp_wise[i] + shrub_mean_spp_wise[i] * mean(shrub.pred) + 
                               pasture_mean_spp_wise[i])
}

##----------------------------------------------------------------
## For how many species is predicted occupancy higher with the interaction in
## the model?
sum(pred_home_int > pred_home) 
sum(pred_pasture_int > pred_pasture)

## And of the subset of species where the interaction is supported by WAIC?
## Species numbers copied from compare_WAIC.R
sum(pred_home_int[c(1,8,12,20,22,25,29,31,35,36,37,39,41,43,44,45,47)] > 
      pred_home[c(1,8,12,20,22,25,29,31,35,36,37,39,41,43,44,45,47)]) 
sum(pred_pasture_int[c(1,8,12,20,22,25,29,31,35,36,37,39,41,43,44,45,47)] > 
      pred_pasture[c(1,8,12,20,22,25,29,31,35,36,37,39,41,43,44,45,47)])
##---------------

## models with no interaction and sugar cane included
load(("./Occupancy_analysis_Swazi_birds_Bayesian_quadTime_fixed.RData"))

Outputs <- list(Output1,Output2,Output3,Output4,Output5,
                Output6,Output7,Output8,
                Output9,Output10,Output11,Output12,
                Output13,Output14,Output15,Output16,
                Output17,Output18,Output19,Output20,
                Output21,Output22,Output23,Output24,
                Output25,Output26,Output27,Output28,
                Output29,Output30,Output31,Output32,
                Output33,Output34,Output35,Output36,
                Output37,Output38,Output39,Output40,
                Output41,Output42,Output43,Output44,
                Output45,Output46,Output47,Output48)

## The objects Output1... Output48 are Large jagsUI; Only the posterior MCMC chains are needed.
Outputs <- ((map(Outputs, MCMCchains)))     ## MCMCchains extracts the chains, map applies 
## extraction to each item in the list 


## take a random draw of ~2000 iterations from each posterior to get equal sample sizes
to_keep <- matrix(nrow = 2000, ncol = 48)
for(i in 1:48){
  to_keep[,i] <- sample(1:dim(Outputs[[i]])[1], size = 2000)                  
  Outputs[[i]] <- subset(Outputs[[i]][c(to_keep[,i]), c(4, 5, 7, 8, 9)])  
}                  ## keep row numbers that are in column i of to_keep, for params of interest

## Assemble outputs from each species in a dataframe in long format
Outputs_long <- as.data.frame (do.call(rbind, Outputs))

## create a column for species name for downstream processing
sppName <- vector(mode = "character", length = dim(Outputs[[1]])[1])  
sppName <- rep(dimnames(AbunHists.subset)[[3]], dim(Outputs[[1]])[1])  ## TraitData$Species
sppName <- sort(sppName)
Outputs_long$Species <- sppName 
Outputs_long$X <- rep(1:dim(Outputs[[1]])[1], 48)   ## add a row number for dcast to convert long to wide format

#subset by beta
shrub.post<-subset(Outputs_long, select=c("X", "Species","beta1"))
homestead.post<-subset(Outputs_long, select=c("X", "Species","beta.l[3]"))
pasture.post<-subset(Outputs_long, select=c("X", "Species","beta.l[2]"))  
intercept.post<-subset(Outputs_long, select=c("X", "Species","beta")) ## intercept for this parameterization
sugarEstate.post<-subset(Outputs_long, select=c("X", "Species","beta.l[4]"))

## reformat to make species columns and posterior samples rows with dcast
## then drop the 'X' column
shrub.post <- dcast(shrub.post, X~Species)
shrub.post <- shrub.post[,2:49]
homestead.post <- dcast(homestead.post, X~Species)
homestead.post <- homestead.post[,2:49]
pasture.post <- dcast(pasture.post, X~Species)
pasture.post <- pasture.post[,2:49]
intercept.post <- dcast(intercept.post, X~Species)
intercept.post <- intercept.post[,2:49]
sugarEstate.post <- dcast(sugarEstate.post, X~Species)
sugarEstate.post <- sugarEstate.post[,2:49]

intercept <- rowMeans(intercept.post)
shrub <- rowMeans(shrub.post) 
pasture <- rowMeans(pasture.post)
homestead <- rowMeans(homestead.post)
plantation <- rowMeans(sugarEstate.post)

predC_home <- array(NA, dim = c(334, 2000, 1))
predC_pasture <- array(NA, dim = c(334, 2000, 1))
predC_not_home <- array(NA, dim = c(334, 2000, 1))
predC_plantation <- array(NA, dim = c(334, 2000, 1))

for (i in 1:2000) {
  predC_home[, i, 1] <- plogis(intercept[i] + shrub[i] * shrub.pred + homestead[i])
  predC_pasture[, i, 1] <- plogis(intercept[i] + shrub[i] * shrub.pred + pasture[i])
  predC_not_home[, i, 1] <- plogis(intercept[i] + shrub[i] * shrub.pred)
  predC_plantation[, i, 1] <- plogis(intercept[i] + shrub[i]*shrub.pred + plantation[i])
}

pmC_home <- apply(predC_home, 1, mean)
criC_home <- apply(predC_home, 1, function(x) quantile(x, prob = c(0.025, 0.975)))
pmC_pasture <- apply(predC_pasture, 1, mean)
criC_pasture <- apply(predC_pasture, 1, function(x) quantile(x, prob = c(0.025, 0.975)))
pmC_not_home <- apply(predC_not_home, 1, mean)
criC_not_home <- apply(predC_not_home, 1, function(x) quantile(x, prob = c(0.025, 0.975)))
pmC_plantation <- apply(predC_plantation, 1, mean)
criC_plantation <- apply(predC_plantation, 1, function(x) quantile(x, prob = c(0.025, 0.975)))
    
# plotting --------------------------------------------
## Restrict to observed values of shrub cover
plot(o.shrub[1:119], pmC_home[1:119], col = "turquoise1", lwd = 3, type = 'l', lty = 1, frame = FALSE, ylim = c(0, 1),
     ylab = "Mean occupancy", las = 1, xlim = c(0,100),
     xlab = "Shrub cover, %", main = "4 land uses, no land use by shrub cover interaction")
 ##matlines(o.shrub[1:119], criC_home[1,1:119], col = "turquoise1", lty = 1)
 ##matlines(o.shrub[1:119], criC_home[2,1:119], col = "turquoise1", lty = 1)

lines(o.shrub, pmC_pasture, col = "green", lwd = 3, type = 'l', lty = 2)
 ##matlines(o.shrub, criC_pasture[1,], col = "green", lty = 2)
 ##matlines(o.shrub, criC_pasture[2,], col = "green", lty = 2)

lines(o.shrub, pmC_not_home, col = "blue", lwd = 3, type = 'l', lty = 3)
 ##matlines(o.shrub, criC_not_home[1,], col = "blue", lty = 3)
 ##matlines(o.shrub, criC_not_home[2,], col = "blue", lty = 3)

## Restrict to observed values
lines(o.shrub[1:82], pmC_plantation[1:82], col = "goldenrod", lwd = 3, type = 'l', lty = 4)
## matlines(o.shrub[1:82], criC_plantation[1,1:82], col = "goldenrod", lty = 4)
## matlines(o.shrub[1:82], criC_plantation[2,1:82], col = "goldenrod", lty = 4)

## Plot predictions for each spp, One panel per species
shrub_mean_spp_wise <- apply(shrub.post, 2, mean)
intercept_mean_spp_wise <- apply(intercept.post, 2, mean)
pasture_mean_spp_wise <- apply(pasture.post, 2, mean)
homestead_mean_spp_wise <- apply(homestead.post, 2, mean)
plantation_mean_spp_wise <- apply(sugarEstate.post, 2, mean)

shrub_cri_spp_wise <- apply(shrub.post, 2, function(x) quantile (x, prob = c(0.025, 0.975)))
intercept_cri_spp_wise <- apply(intercept.post, 2, function(x) quantile (x, prob = c(0.025, 0.975)))
pasture_cri_spp_wise <- apply(pasture.post, 2, function(x) quantile (x, prob = c(0.025, 0.975)))
homestead_cri_spp_wise <- apply(homestead.post, 2, function(x) quantile (x, prob = c(0.025, 0.975)))
plantation_cri_spp_wise <- apply(sugarEstate.post, 2, function(x) quantile (x, prob = c(0.025, 0.975)))

pred_home <- matrix(NA, 48, 334)
pred_pasture <- matrix(NA, 48, 334)
pred_not_home <- matrix(NA, 48, 334)
pred_plantation <- matrix(NA, 48, 334)

for (i in 1:48) {
  pred_home[i,] <- plogis(intercept_mean_spp_wise[i] + shrub_mean_spp_wise[i] * shrub.pred + 
                           homestead_mean_spp_wise[i])
  pred_pasture[i,] <- plogis(intercept_mean_spp_wise[i] + shrub_mean_spp_wise[i] * shrub.pred + 
                              pasture_mean_spp_wise[i])
  pred_not_home[i,] <- plogis(intercept_mean_spp_wise[i] + shrub_mean_spp_wise[i] * shrub.pred)
  pred_plantation[i,] <- plogis(intercept_mean_spp_wise[i] + shrub_mean_spp_wise[i] * shrub.pred +
                                plantation_mean_spp_wise[i])
}

tiff("./outputs/Fig_S1.tiff", width = 18, height = 18, units = 'in', res = 300)
par(mfrow = c(6, 8))
for (j in 1:48) {
  plot(o.shrub[1:119], pred_home[j,1:119], col = "turquoise1", lwd = 3, type = 'l', lty = 1, frame = FALSE, 
       ylim = c(0, 1), ylab = "Mean occupancy", las = 1, xlim = c(0, 100),
       xlab = "Shrub cover, %", main = dimnames(AbunHists.subset)[[3]][j])
  lines(o.shrub, pred_pasture[j,], col = "green", lwd = 3, type = 'l', lty = 2)
  lines(o.shrub, pred_not_home[j,], col = "blue", lwd = 3, type = 'l', lty = 3)
  lines(o.shrub[1:82], pred_plantation[j,1:82], col = "goldenrod", lwd = 3, type = 'l', lty = 4)
}
dev.off()

## Plot all land use effects, one panel per land use, one line per species
tiff("./outputs/Fig_SXX.tiff", width = 8, height = 8, units = 'in', res = 300)
par(mfrow = c(2,2))
## Protected areas
j=1
plot(o.shrub, pred_not_home[j,], col = "blue", lwd = 2, type = 'l', lty = 1, frame = FALSE, 
     ylim = c(0, 1), ylab = "Mean occupancy", las = 1, xlim = c(0,100),
     xlab = "", main = "Protected areas")

for (j in 2:48) {
  lines(o.shrub, pred_not_home[j,], col = "blue", lwd = 2, type = 'l', lty = 1)
}

## Pastures
j=1
plot(o.shrub, pred_pasture[j,], col = "green", lwd = 2, type = 'l', lty = 1, frame = FALSE, 
     ylim = c(0, 1), ylab = "", las = 1, xlab = "", main = "Pastures", xlim = c(0,100))

for (j in 2:48) {
  lines(o.shrub, pred_pasture[j,], col = "green", lwd = 2, type = 'l', lty = 1)
}

## Homesteads
j=1
plot(o.shrub[1:119], pred_home[j,1:119], col = "turquoise1", lwd = 2, type = 'l', lty = 1, frame = FALSE, 
     ylim = c(0, 1), ylab = "Mean occupancy", las = 1, xlim = c(0, 100),
     xlab = "Shrub cover, %", main = "Homesteads")

for (j in 2:48) {
  lines(o.shrub[1:119], pred_home[j,1:119], col = "turquoise1", lwd = 2, type = 'l', lty = 1)
}

## Plantation
j=1
plot(o.shrub[1:82], pred_plantation[j,1:82], col = "goldenrod", lwd = 2, type = 'l', lty = 1, frame = FALSE, 
     ylim = c(0, 1), ylab = "", las = 1, xlim = c(0, 100),
     xlab = "Shrub cover, %", main = "Plantation")

for (j in 2:48) {
  lines(o.shrub[1:82], pred_plantation[j,1:82], col = "goldenrod", lwd = 2, type = 'l', lty = 1)
}
dev.off()
