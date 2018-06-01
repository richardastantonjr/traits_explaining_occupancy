# Summary of what this script does ----------------------------------------
## Extract posterior MCMC chains and process to get posterior 95% CRIs for suites of species with
## common traits, regress trait values on shrub cover and land-use effects, and display results
## in figures

# functions ---------------------------------------------------------------
# return mean and 95% CRI
meanAndCRI <- function(x) {
  c(mean <- mean(x),quantile(x,c(0.025, 0.975))) 
}

meanAnd90CRI <- function(x) {
  c(mean <- mean(x),quantile(x,c(0.05, 0.95))) 
}

# libraries ----------------------------------------------------------------------------
library(tidyverse)
library(reshape2) 
library(gridExtra)
library(MCMCvis)

# load data -----------------------------------------------------------------------------
TraitData <- read.csv("./data/TraitData.csv")
TraitData[, 1] <- as.factor(TraitData[,1])

## load models with quadratic time of day effects and effects parameterization 
load("./outputs/Occupancy_analysis_Swazi_birds_Bayesian_quadTime_WAIC_20180528.RData")

# subset and process data  -----------------------------------------------------------------------
TraitData <- filter(TraitData,TraitData$Species %in% dimnames(AbunHists.subset)[[3]]) 
TraitData <- TraitData[,c(1:6, 9:14, 23, 24, 26)]
TraitData[is.na(TraitData)] <- 0 
TraitData$Species <- droplevels(TraitData$Species)

## determine if each species is predatory and include that as a trait
predator<-rep(NA, nrow(TraitData))
for (i in 1:nrow(TraitData)){
  if((TraitData[i, 13] > 0) || (TraitData[i, 14] > 0) || (TraitData[i, 15] > 0)){
    predator[i] = 1
  } else {
    predator[i] = 0} 
}
TraitData$predator <- predator

## Calculate and add pseudo wing loading as a trait
TraitData$pseudo_loading <- TraitData$Mass / TraitData$Wing

## Process the model Output
## Put output from jagsUI in a list- each corresponds to a single species
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

## Extracts and assemble DIC and deviance for each model
dic_vals_int_models <- unlist(map(Outputs, "DIC"))

deviance_vals_int_models <- unlist(map(Outputs, function(x)
{MCMCsummary(x, params = 'deviance')[1]}))

## The objects Output1... Output48 are Large jagsUI; Only the posterior MCMC chains are needed.
Outputs <- ((map(Outputs, MCMCchains))) ## MCMCchains extracts the chains, map applies 
                                        ## extraction to each item in the list 

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

## remove extraneous material
rm(list=ls()[! ls() %in% c("Outputs_long","TraitData", "meanAndCRI", "meanAnd90CRI",
                           "dic_vals_int_models", "deviance_vals_int_models" )])
##----------------------------------------------------------------------------------
##----------------------------------------------------------------------------------
#create vectors of species by trait
## nesting substrates
nest.shrub <- subset(TraitData, ShrubNest == 1, select=c(Species))
nest.cavity <- subset(TraitData, CavityNest == 1, select=c(Species))
nest.tree <- subset(TraitData, TreeNest == 1, select=c(Species))
nest.grass <- subset(TraitData, GrassNest == 1, select=c(Species))
 
## diets
predators <- subset(TraitData, predator == 1, select=c(Species))
diet.invert <- subset(TraitData, Diet == "invertebrates", select=c(Species))
diet.fruit <- subset(TraitData, Diet == "fruit", select=c(Species))
diet.seed <- subset(TraitData, Diet == "seeds", select=c(Species))
diet.nectar <- subset(TraitData, Diet == "nectar", select=c(Species))

## merge posteriors and species traits 
Outputs_long <- merge(Outputs_long, TraitData)

#subset by effect size or beta, keeping posterior draw, species beta:
shrub.post<-subset(Outputs_long, select = c("X", "Species","beta1"))
homestead.post<-subset(Outputs_long, select = c("X", "Species","beta.l[3]"))
pasture.post<-subset(Outputs_long, select = c("X", "Species","beta.l[2]"))  
protected.post<-subset(Outputs_long, select = c("X", "Species","beta")) ## intercept for this parameterization
pas_by_shrub.post <- subset(Outputs_long, select = c("X", "Species","beta.l.int[2]"))  
home_by_shrub.post <- subset(Outputs_long, select = c("X", "Species","beta.l.int[3]"))  

#reformat to make species columns and posterior samples rows with dcast
shrub.post <- dcast(shrub.post, X~Species)
shrub.post <- shrub.post[,2:49]
homestead.post <- dcast(homestead.post, X~Species)
homestead.post <- homestead.post[,2:49]
pasture.post <- dcast(pasture.post, X~Species)
pasture.post <- pasture.post[,2:49]
intercept.post <- dcast(protected.post, X~Species)
intercept.post <- intercept.post[,2:49]
pas_by_shrub.post <- dcast(pas_by_shrub.post, X~Species)
pas_by_shrub.post <- pas_by_shrub.post[,2:49]
home_by_shrub.post <- dcast(home_by_shrub.post, X~Species)
home_by_shrub.post <- home_by_shrub.post[,2:49]

## summarize results for the community of common species
pooledInterceptCRI <- meanAndCRI(rowMeans(intercept.post))
pooledHomesteadCRI <-  meanAndCRI(rowMeans(homestead.post))
pooledPastureCRI <-  meanAndCRI(rowMeans(pasture.post))
pooledShrubCRI <- meanAndCRI(rowMeans(shrub.post))
pas_by_shrubCRI <- meanAndCRI(rowMeans(pas_by_shrub.post))
home_by_shrubCRI <- meanAndCRI(rowMeans(home_by_shrub.post))

pooledIntercept90CRI <- meanAnd90CRI(rowMeans(intercept.post))
pooledHomestead90CRI <- meanAnd90CRI(rowMeans(homestead.post))
pooledPasture90CRI <- meanAnd90CRI(rowMeans(pasture.post))
pooledShrub90CRI <- meanAnd90CRI(rowMeans(shrub.post))  
pas_by_shrub90CRI <- meanAnd90CRI(rowMeans(pas_by_shrub.post))
home_by_shrub90CRI <- meanAnd90CRI(rowMeans(home_by_shrub.post))

res_pooled_95 <- rbind(pooledInterceptCRI,
                       pooledShrubCRI,
                       pooledPastureCRI, 
                       pooledHomesteadCRI,
                       pas_by_shrubCRI,
                       home_by_shrubCRI)

res_pooled_90 <- rbind(pooledIntercept90CRI,
                       pooledShrub90CRI,
                       pooledPasture90CRI, 
                       pooledHomestead90CRI,
                       pas_by_shrub90CRI,
                       home_by_shrub90CRI)

res_pooled <- cbind(res_pooled_95, res_pooled_90)
colnames(res_pooled) <- c("Mean", "95% LCRI", "95% UCRI", "Mean", "90% LCRI", "90% UCRI")
rownames(res_pooled) <- c("Intercept", "Shrub cover", "Pasture","Homestead",
                          "Pasture x shrub", "Homestead x shrub")
res_pooled <- res_pooled[,-4] ## Remove duplicated mean
write.csv(res_pooled, "./outputs/AppendixS3/AppendixS3.2.csv")

## plot Community-level interaction posterior densities 
plot(density(rowMeans(pas_by_shrub.post)), xlab = "Land use x shrub cover effect size (beta coeff.)", main = "",
     xlim = c(-10, 10), col = "blue", lwd = 3) 
lines(density(rowMeans(home_by_shrub.post)), col = "salmon", lty = 2, lwd = 3) 
abline(v=0)
legend("topleft", c("Pasture x shrub","Homestead x shrub"),
       fill = c("blue", "salmon"), bty = "n")

## Community-level main effects of each covariate
plot(density(rowMeans(shrub.post)), xlim = c(-8, 8), xlab = "Beta estimate", 
     main = "", lwd = 3, las = 1, ylim = c(0, 1))
lines(density(rowMeans(pasture.post)), lty = 2, col = "blue", lwd = 3)
lines(density(rowMeans(homestead.post)), lty = 3, col = "orange", lwd = 3)
abline(v=0)
legend("topleft", c("shrub","pasture","homestead"),
       fill = c("black", "blue","orange"), bty = "n")

#take rowMeans based on species list
shrub.post$nest.shrub <- rowMeans(shrub.post[,names(shrub.post) %in% nest.shrub$Species])
shrub.post$nest.tree <- rowMeans(shrub.post[,names(shrub.post) %in% nest.tree$Species])
shrub.post$nest.cavity <- rowMeans(shrub.post[,names(shrub.post) %in% nest.cavity$Species])
shrub.post$nest.grass <- rowMeans(shrub.post[,names(shrub.post) %in% nest.grass$Species])

shrub.post$predators <- rowMeans(shrub.post[,names(shrub.post) %in% predators$Species])
shrub.post$diet.invert <- rowMeans(shrub.post[,names(shrub.post) %in% diet.invert$Species])
shrub.post$diet.fruit <- rowMeans(shrub.post[,names(shrub.post) %in% diet.fruit$Species])
shrub.post$diet.seeds <- rowMeans(shrub.post[,names(shrub.post) %in% diet.seed$Species])
shrub.post$diet.nectar <- rowMeans(shrub.post[,names(shrub.post) %in% diet.nectar$Species])

intercept.post$nest.shrub <- rowMeans(intercept.post[,names(intercept.post) %in% nest.shrub$Species])
intercept.post$nest.tree <- rowMeans(intercept.post[,names(intercept.post) %in% nest.tree$Species])
intercept.post$nest.cavity <- rowMeans(intercept.post[,names(intercept.post) %in% nest.cavity$Species])
intercept.post$nest.grass <- rowMeans(intercept.post[,names(intercept.post) %in% nest.grass$Species])

intercept.post$predators <- rowMeans(intercept.post[,names(intercept.post) %in% predators$Species])
intercept.post$diet.invert <- rowMeans(intercept.post[,names(intercept.post) %in% diet.invert$Species])
intercept.post$diet.fruit <- rowMeans(intercept.post[,names(intercept.post) %in% diet.fruit$Species])
intercept.post$diet.seeds <- rowMeans(intercept.post[,names(intercept.post) %in% diet.seed$Species])
intercept.post$diet.nectar <- rowMeans(intercept.post[,names(intercept.post) %in% diet.nectar$Species])

pasture.post$nest.shrub <- rowMeans(pasture.post[,names(pasture.post) %in% nest.shrub$Species])
pasture.post$nest.tree <- rowMeans(pasture.post[,names(pasture.post) %in% nest.tree$Species])
pasture.post$nest.cavity <- rowMeans(pasture.post[,names(pasture.post) %in% nest.cavity$Species])
pasture.post$nest.grass <- rowMeans(pasture.post[,names(pasture.post) %in% nest.grass$Species])

pasture.post$predators <- rowMeans(pasture.post[,names(pasture.post) %in% predators$Species])
pasture.post$diet.invert <- rowMeans(pasture.post[,names(pasture.post) %in% diet.invert$Species])
pasture.post$diet.fruit <- rowMeans(pasture.post[,names(pasture.post) %in% diet.fruit$Species])
pasture.post$diet.seeds <- rowMeans(pasture.post[,names(pasture.post) %in% diet.seed$Species])
pasture.post$diet.nectar <- rowMeans(pasture.post[,names(pasture.post) %in% diet.nectar$Species])

homestead.post$nest.shrub <- rowMeans(homestead.post[,names(homestead.post) %in% nest.shrub$Species])
homestead.post$nest.tree <- rowMeans(homestead.post[,names(homestead.post) %in% nest.tree$Species])
homestead.post$nest.cavity <- rowMeans(homestead.post[,names(homestead.post) %in% nest.cavity$Species])
homestead.post$nest.grass <- rowMeans(homestead.post[,names(homestead.post) %in% nest.grass$Species])

homestead.post$predators <- rowMeans(homestead.post[,names(homestead.post) %in% predators$Species])
homestead.post$diet.invert <- rowMeans(homestead.post[,names(homestead.post) %in% diet.invert$Species])
homestead.post$diet.fruit <- rowMeans(homestead.post[,names(homestead.post) %in% diet.fruit$Species])
homestead.post$diet.seeds <- rowMeans(homestead.post[,names(homestead.post) %in% diet.seed$Species])
homestead.post$diet.nectar <- rowMeans(homestead.post[,names(homestead.post) %in% diet.nectar$Species])

## interaction terms
pas_by_shrub.post$nest.shrub <- rowMeans(pas_by_shrub.post[,names(pas_by_shrub.post) %in% nest.shrub$Species])
pas_by_shrub.post$nest.tree <- rowMeans(pas_by_shrub.post[,names(pas_by_shrub.post) %in% nest.tree$Species])
pas_by_shrub.post$nest.cavity <- rowMeans(pas_by_shrub.post[,names(pas_by_shrub.post) %in% nest.cavity$Species])
pas_by_shrub.post$nest.grass <- rowMeans(pas_by_shrub.post[,names(pas_by_shrub.post) %in% nest.grass$Species])

pas_by_shrub.post$predators <- rowMeans(pas_by_shrub.post[,names(pas_by_shrub.post) %in% predators$Species])
pas_by_shrub.post$diet.invert <- rowMeans(pas_by_shrub.post[,names(pas_by_shrub.post) %in% diet.invert$Species])
pas_by_shrub.post$diet.fruit <- rowMeans(pas_by_shrub.post[,names(pas_by_shrub.post) %in% diet.fruit$Species])
pas_by_shrub.post$diet.seeds <- rowMeans(pas_by_shrub.post[,names(pas_by_shrub.post) %in% diet.seed$Species])
pas_by_shrub.post$diet.nectar <- rowMeans(pas_by_shrub.post[,names(pas_by_shrub.post) %in% diet.nectar$Species])

home_by_shrub.post$nest.shrub <- rowMeans(home_by_shrub.post[,names(home_by_shrub.post) %in% nest.shrub$Species])
home_by_shrub.post$nest.tree <- rowMeans(home_by_shrub.post[,names(home_by_shrub.post) %in% nest.tree$Species])
home_by_shrub.post$nest.cavity <- rowMeans(home_by_shrub.post[,names(home_by_shrub.post) %in% nest.cavity$Species])
home_by_shrub.post$nest.grass <- rowMeans(home_by_shrub.post[,names(home_by_shrub.post) %in% nest.grass$Species])

home_by_shrub.post$predators <- rowMeans(home_by_shrub.post[,names(home_by_shrub.post) %in% predators$Species])
home_by_shrub.post$diet.invert <- rowMeans(home_by_shrub.post[,names(home_by_shrub.post) %in% diet.invert$Species])
home_by_shrub.post$diet.fruit <- rowMeans(home_by_shrub.post[,names(home_by_shrub.post) %in% diet.fruit$Species])
home_by_shrub.post$diet.seeds <- rowMeans(home_by_shrub.post[,names(home_by_shrub.post) %in% diet.seed$Species])
home_by_shrub.post$diet.nectar <- rowMeans(home_by_shrub.post[,names(home_by_shrub.post) %in% diet.nectar$Species])

## create a data frame of effect sizes and CRIs for each species for plotting species names
Species <- levels(as.factor(Outputs_long$Species))

## create a data frame of one effect and CRIs for all species
shrubSppEfx <- apply(shrub.post, 2, meanAndCRI )
shrubSppEfx <- as.data.frame(t(shrubSppEfx[, 2:49]))
shrubSppEfx <- cbind(Species, shrubSppEfx)
colnames(shrubSppEfx) <- c("Species", "effect", "LCL", "UCL")
shrubSppEfx$betaName <- "Shrub cover"
rownames(shrubSppEfx)<- NULL

## repeat for the remaining effects
homesteadSppEfx <- apply(homestead.post[,1:48], 2, meanAndCRI )
homesteadSppEfx <- as.data.frame(t(homesteadSppEfx))
homesteadSppEfx <- cbind(Species, homesteadSppEfx)
colnames(homesteadSppEfx) <- c("Species", "effect", "LCL", "UCL")
homesteadSppEfx$betaName <- "Homestead"
rownames(homesteadSppEfx)<- NULL

pastureSppEfx <- apply(pasture.post[,1:48], 2, meanAndCRI )
pastureSppEfx <- as.data.frame(t(pastureSppEfx))
pastureSppEfx <- cbind(Species, pastureSppEfx)
colnames(pastureSppEfx) <- c("Species", "effect", "LCL", "UCL")
pastureSppEfx$betaName <- "Pasture"
rownames(pastureSppEfx)<- NULL

interceptSppEfx <- apply(intercept.post[,1:48], 2, meanAndCRI )
interceptSppEfx <- as.data.frame(t(interceptSppEfx))
interceptSppEfx <- cbind(Species, interceptSppEfx)
colnames(interceptSppEfx) <- c("Species", "effect", "LCL", "UCL")
interceptSppEfx$betaName <- "Intercept"
rownames(interceptSppEfx)<- NULL

pas_by_shrubSppEfx <- apply(pas_by_shrub.post[,1:48], 2, meanAndCRI )
pas_by_shrubSppEfx <- as.data.frame(t(pas_by_shrubSppEfx))
pas_by_shrubSppEfx <- cbind(Species, pas_by_shrubSppEfx)
colnames(pas_by_shrubSppEfx) <- c("Species", "effect", "LCL", "UCL")
pas_by_shrubSppEfx$betaName <- "pas_by_shrub"
rownames(pas_by_shrubSppEfx)<- NULL

home_by_shrubSppEfx <- apply(home_by_shrub.post[,1:48], 2, meanAndCRI )
home_by_shrubSppEfx <- as.data.frame(t(home_by_shrubSppEfx))
home_by_shrubSppEfx <- cbind(Species, home_by_shrubSppEfx)
colnames(home_by_shrubSppEfx) <- c("Species", "effect", "LCL", "UCL")
home_by_shrubSppEfx$betaName <- "home_by_shrub"
rownames(home_by_shrubSppEfx)<- NULL

## bind rows to create one frame for plotting 
sppEfx <- rbind.data.frame(shrubSppEfx, interceptSppEfx, pastureSppEfx, homesteadSppEfx,
                           pas_by_shrubSppEfx, home_by_shrubSppEfx)
sppEfx <- merge(sppEfx, TraitData[,1:2]) ## sort phylogenetically
sppEfx<- sppEfx[order(sppEfx$RobertsNum), ]

## add an indicator variable to distinguish significant positive and negative effects from inconclusive ones
for(i in 1:dim(sppEfx)[1]){
  sppEfx$POSNEG[i] <- if_else (all(c(sppEfx$LCL[i] < 0, sppEfx$UCL[i] < 0) == TRUE), "Neg", 
                               if_else(all(c(sppEfx$LCL[i] > 0, sppEfx$UCL[i] > 0) == TRUE), "Pos", "Inconclusive"))
}

## mass/wing chord against mean betas
shrubEfxAndTraits <- merge(sppEfx, TraitData) %>% 
  filter(betaName == "Shrub cover")
protEfxAndTraits <- merge(sppEfx, TraitData) %>% 
  filter(betaName == "Protected")
pastEfxAndTraits <- merge(sppEfx, TraitData) %>% 
  filter(betaName == "Pasture")
homeEfxAndTraits <- merge(sppEfx, TraitData) %>% 
  filter(betaName == "Homestead")
pas_by_shrubEfxAndTraits <- merge(sppEfx, TraitData) %>% 
  filter(betaName == "pas_by_shrub")
home_by_shrubEfxAndTraits <- merge(sppEfx, TraitData) %>% 
  filter(betaName == "home_by_shrub")

####
## full posterior samples. Break into 2000 data sets for each beta,
## each representing a posterior draw and associated continuous trait values.
## Should be condensed into a function
##-------------------------
## Shrub cover effect
beta1_samp <- list()
pseudo_loadings <- list()
masses <- list()
for (m in 1:2000) {
  beta1_samp[[m]] <- Outputs_long$beta1[c(seq(m,96000,2000))]
  pseudo_loadings[[m]] <- Outputs_long$pseudo_loading[c(seq(m,96000,2000))]  
  masses[[m]] <- Outputs_long$Mass[c(seq(m,96000,2000))]  
}

beta1_slopes_pseudo_load <- NULL
beta1_slopes_mass <- NULL
for (mm in 1:2000) {
  beta1_slopes_pseudo_load[mm] <- coef(lm(beta1_samp[[mm]] ~ pseudo_loadings[[mm]]))[2]
  beta1_slopes_mass[mm] <- coef(lm(beta1_samp[[mm]] ~ masses[[mm]]))[2]
} 
##-----------------

## Protection effect
beta_samp <- list()
for (m in 1:2000) {
  beta_samp[[m]] <- Outputs_long$beta[c(seq(m,96000,2000))]
}

beta_slopes_pseudo_load <- NULL
beta_slopes_mass <- NULL
for (mm in 1:2000) {
  beta_slopes_pseudo_load[mm] <- coef(lm(beta_samp[[mm]] ~ pseudo_loadings[[mm]]))[2]
  beta_slopes_mass[mm] <- coef(lm(beta_samp[[mm]] ~ masses[[mm]]))[2]
} 
##-------------------

## Pasture effect
betaL2_samp <- list()
for (m in 1:2000) {
  betaL2_samp[[m]] <- Outputs_long$`beta.l[2]`[c(seq(m,96000,2000))]
}

betaL2_slopes_pseudo_load <- NULL
betaL2_slopes_mass <- NULL
for (mm in 1:2000) {
  betaL2_slopes_pseudo_load[mm] <- coef(lm(betaL2_samp[[mm]] ~ pseudo_loadings[[mm]]))[2]
  betaL2_slopes_mass[mm] <- coef(lm(betaL2_samp[[mm]] ~ masses[[mm]]))[2]
} 
##-----------------------

## Homestead effect
betaL3_samp <- list()
for (m in 1:2000) {
  betaL3_samp[[m]] <- Outputs_long$`beta.l[3]`[c(seq(m,96000,2000))]
}

betaL3_slopes_pseudo_load <- NULL
betaL3_slopes_mass <- NULL
for (mm in 1:2000) {
  betaL3_slopes_pseudo_load[mm] <- coef(lm(betaL3_samp[[mm]] ~ pseudo_loadings[[mm]]))[2]
  betaL3_slopes_mass[mm] <- coef(lm(betaL3_samp[[mm]] ~ masses[[mm]]))[2]
} 

## pas_by_shrub effect
betaLI2_samp <- list()
for (m in 1:2000) {
  betaLI2_samp[[m]] <- Outputs_long$`beta.l.int[2]`[c(seq(m,96000,2000))]
}

betaLI2_slopes_pseudo_load <- NULL
betaLI2_slopes_mass <- NULL
for (mm in 1:2000) {
  betaLI2_slopes_pseudo_load[mm] <- coef(lm(betaLI2_samp[[mm]] ~ pseudo_loadings[[mm]]))[2]
  betaLI2_slopes_mass[mm] <- coef(lm(betaLI2_samp[[mm]] ~ masses[[mm]]))[2]
} 

## home_by_shrub effect
betaLI3_samp <- list()
for (m in 1:2000) {
  betaLI3_samp[[m]] <- Outputs_long$`beta.l.int[3]`[c(seq(m,96000,2000))]
}

betaLI3_slopes_pseudo_load <- NULL
betaLI3_slopes_mass <- NULL
for (mm in 1:2000) {
  betaLI3_slopes_pseudo_load[mm] <- coef(lm(betaLI3_samp[[mm]] ~ pseudo_loadings[[mm]]))[2]
  betaLI3_slopes_mass[mm] <- coef(lm(betaLI3_samp[[mm]] ~ masses[[mm]]))[2]
} 
##--------------------------

## bundle the results
slopes_95 <- rbind.data.frame(meanAndCRI(beta1_slopes_mass),
                              meanAndCRI(beta_slopes_mass),
                              meanAndCRI(betaL2_slopes_mass),
                              meanAndCRI(betaL3_slopes_mass),
                              meanAndCRI(betaLI2_slopes_mass),
                              meanAndCRI(betaLI3_slopes_mass),
                              
                              meanAndCRI(beta1_slopes_pseudo_load),
                              meanAndCRI(beta_slopes_pseudo_load),
                              meanAndCRI(betaL2_slopes_pseudo_load),
                              meanAndCRI(betaL3_slopes_pseudo_load),
                              meanAndCRI(betaLI2_slopes_pseudo_load),
                              meanAndCRI(betaLI3_slopes_pseudo_load)
)
colnames(slopes_95) <- c("effect", "LCL", "UCL")
slopes_95$betaName <- c("Shrub","Intercept","Pasture","Homestead",
                      "Pas_by_shrub", "Home_x_shrub", 
                      "Shrub","Intercept","Pasture","Homestead",
                      "Pas_by_shrub","Home_x_shrub")
slopes_95$measure <- c("mass","mass","mass","mass","mass","mass", 
                       "loading","loading","loading",
                       "loading","loading", "loading")



slopes_90 <- rbind.data.frame(meanAnd90CRI(beta1_slopes_mass),
                              meanAnd90CRI(beta_slopes_mass),
                              meanAnd90CRI(betaL2_slopes_mass),
                              meanAnd90CRI(betaL3_slopes_mass),
                              meanAnd90CRI(betaLI2_slopes_mass),
                              meanAnd90CRI(betaLI3_slopes_mass),
                              
                              meanAnd90CRI(beta1_slopes_pseudo_load),
                              meanAnd90CRI(beta_slopes_pseudo_load),
                              meanAnd90CRI(betaL2_slopes_pseudo_load),
                              meanAnd90CRI(betaL3_slopes_pseudo_load),
                              meanAnd90CRI(betaLI2_slopes_pseudo_load),
                              meanAnd90CRI(betaLI3_slopes_pseudo_load)
)
colnames(slopes_90) <- c("effect", "LCL", "UCL")
slopes_90$betaName <- c("Shrub","Intercept","Pasture","Homestead",
                        "Pas_by_shrub", "Home_x_shrub", 
                        "Shrub","Intercept","Pasture","Homestead",
                        "Pas_by_shrub","Home_x_shrub")
slopes_90$measure <- c("mass","mass","mass","mass","mass","mass",
                       "loading","loading","loading",
                       "loading","loading", "loading")

cont_trait_slopes <- cbind.data.frame(slopes_95, slopes_90)
cont_trait_slopes <- cont_trait_slopes[,c(1:3, 7:10)]
colnames(cont_trait_slopes)[2:5] <- c("95% LCL", "95% UCL", "90%LCL", "90% UCL")
write.csv(cont_trait_slopes, "./outputs/AppendixS4/AppendixS4cont_trait_slopes.2.csv")
##-----------------------------------------------------------------------------------

##-----------------------------------------------------------------------------------
## build data frames with nest and diet summary effects for plotting
## Nest substrates
r1 <- meanAndCRI(shrub.post$nest.shrub)
r2 <- meanAndCRI(shrub.post$nest.tree)
r3 <- meanAndCRI(shrub.post$nest.cavity)
r4 <- meanAndCRI(shrub.post$nest.grass)
shrub.post.nest <- rbind.data.frame(r1,r2,r3,r4)
colnames(shrub.post.nest) <- c("effect", "LCL", "UCL")
shrub.post.nest$betaName <- "Shrub cover"
shrub.post.nest$nest <- c("shrubs", "trees", "cavities", "grass")

r5 <- meanAndCRI (intercept.post$nest.shrub)
r6 <- meanAndCRI (intercept.post$nest.tree)
r7 <- meanAndCRI (intercept.post$nest.cavity)
r8 <- meanAndCRI (intercept.post$nest.grass)
intercept.post.nest <- rbind.data.frame(r5,r6,r7,r8)
colnames(intercept.post.nest) <- c("effect", "LCL","UCL")
intercept.post.nest$betaName <- "Intercept"
intercept.post.nest$nest <- c("shrubs", "trees", "cavities","grass")

r9 <- meanAndCRI (pasture.post$nest.shrub)
r10 <- meanAndCRI (pasture.post$nest.tree)
r11 <- meanAndCRI (pasture.post$nest.cavity)
r12<- meanAndCRI (pasture.post$nest.grass)
pasture.post.nest <- rbind.data.frame(r9,r10,r11,r12)
colnames(pasture.post.nest) <- c("effect", "LCL","UCL")
pasture.post.nest$betaName <- "Pasture"
pasture.post.nest$nest <- c("shrubs", "trees", "cavities","grass")

r13 <- meanAndCRI (homestead.post$nest.shrub)
r14 <- meanAndCRI (homestead.post$nest.tree)
r15 <- meanAndCRI (homestead.post$nest.cavity)
r16 <- meanAndCRI (homestead.post$nest.grass)
homestead.post.nest <- rbind.data.frame(r13,r14,r15,r16)
colnames(homestead.post.nest) <- c("effect", "LCL","UCL")
homestead.post.nest$betaName <- "Homestead"
homestead.post.nest$nest <- c("shrubs", "trees", "cavities","grass")

r21 <- meanAndCRI (pas_by_shrub.post$nest.shrub)
r22 <- meanAndCRI (pas_by_shrub.post$nest.tree)
r23 <- meanAndCRI (pas_by_shrub.post$nest.cavity)
r24 <- meanAndCRI (pas_by_shrub.post$nest.grass)
pas_by_shrub.post.nest <- rbind.data.frame(r21,r22,r23,r24)
colnames(pas_by_shrub.post.nest) <- c("effect", "LCL","UCL")
pas_by_shrub.post.nest$betaName <- "pas_by_shrub"
pas_by_shrub.post.nest$nest <- c("shrubs", "trees", "cavities","grass")

r25 <- meanAndCRI (home_by_shrub.post$nest.shrub)
r26 <- meanAndCRI (home_by_shrub.post$nest.tree)
r27 <- meanAndCRI (home_by_shrub.post$nest.cavity)
r28 <- meanAndCRI (home_by_shrub.post$nest.grass)
home_by_shrub.post.nest <- rbind.data.frame(r25,r26,r27,r28)
colnames(home_by_shrub.post.nest) <- c("effect", "LCL","UCL")
home_by_shrub.post.nest$betaName <- "home_by_shrub"
home_by_shrub.post.nest$nest <- c("shrubs", "trees", "cavities","grass")

## Diets
r29 <- meanAndCRI(shrub.post$predators)
r30 <- meanAndCRI(shrub.post$diet.invert)
r31 <- meanAndCRI(shrub.post$diet.fruit)
r32 <- meanAndCRI(shrub.post$diet.seeds)
r33 <- meanAndCRI(shrub.post$diet.nectar)
shrub.post.diet <- rbind.data.frame(r29,r30,r31,r32,r33)
colnames(shrub.post.diet) <- c("effect", "LCL","UCL")
shrub.post.diet$betaName <- "Shrub cover"
shrub.post.diet$diet <- c("predatory","invertebrates","fruit","seeds","nectar")

r34 <- meanAndCRI (intercept.post$predators)
r35 <- meanAndCRI (intercept.post$diet.invert)
r36 <- meanAndCRI (intercept.post$diet.fruit)
r37 <- meanAndCRI (intercept.post$diet.seeds)
r38 <- meanAndCRI (intercept.post$diet.nectar)
intercept.post.diet <- rbind.data.frame(r34,r35,r36,r37,r38)
colnames(intercept.post.diet) <- c("effect", "LCL","UCL")
intercept.post.diet$betaName <- "Intercept"
intercept.post.diet$diet <- c("predatory","invertebrates","fruit","seeds","nectar")

r39 <- meanAndCRI (pasture.post$predators)
r40 <- meanAndCRI (pasture.post$diet.invert)
r41 <- meanAndCRI (pasture.post$diet.fruit)
r42 <- meanAndCRI (pasture.post$diet.seeds)
r43 <- meanAndCRI (pasture.post$diet.nectar)
pasture.post.diet <- rbind.data.frame(r39,r40,r41,r42,r43)
colnames(pasture.post.diet) <- c("effect", "LCL","UCL")
pasture.post.diet$betaName <- "Pasture"
pasture.post.diet$diet <- c("predatory","invertebrates","fruit","seeds","nectar")

r44 <- meanAndCRI (homestead.post$predators)
r45 <- meanAndCRI (homestead.post$diet.invert)
r46 <- meanAndCRI (homestead.post$diet.fruit)
r47 <- meanAndCRI (homestead.post$diet.seeds)
r48 <- meanAndCRI (homestead.post$diet.nectar)
homestead.post.diet <- rbind.data.frame(r44,r45,r46,r47,r48)
colnames(homestead.post.diet) <- c("effect", "LCL","UCL")
homestead.post.diet$betaName <- "Homestead"
homestead.post.diet$diet <- c("predatory","invertebrates","fruit","seeds","nectar")

r54 <- meanAndCRI (pas_by_shrub.post$predators)
r55 <- meanAndCRI (pas_by_shrub.post$diet.invert)
r56 <- meanAndCRI (pas_by_shrub.post$diet.fruit)
r57 <- meanAndCRI (pas_by_shrub.post$diet.seeds)
r58 <- meanAndCRI (pas_by_shrub.post$diet.nectar)
pas_by_shrub.post.diet <- rbind.data.frame(r54,r55,r56,r57,r58)
colnames(pas_by_shrub.post.diet) <- c("effect", "LCL","UCL")
pas_by_shrub.post.diet$betaName <- "pas_by_shrub"
pas_by_shrub.post.diet$diet <- c("predatory","invertebrates","fruit","seeds","nectar")

r59 <- meanAndCRI (home_by_shrub.post$predators)
r60 <- meanAndCRI (home_by_shrub.post$diet.invert)
r61 <- meanAndCRI (home_by_shrub.post$diet.fruit)
r62 <- meanAndCRI (home_by_shrub.post$diet.seeds)
r63 <- meanAndCRI (home_by_shrub.post$diet.nectar)
home_by_shrub.post.diet <- rbind.data.frame(r59,r60,r61,r62,r63)
colnames(home_by_shrub.post.diet) <- c("effect", "LCL","UCL")
home_by_shrub.post.diet$betaName <- "home_by_shrub"
home_by_shrub.post.diet$diet <- c("predatory","invertebrates","fruit","seeds","nectar")

## collect
nest.summary <- rbind.data.frame(shrub.post.nest, intercept.post.nest, pasture.post.nest, homestead.post.nest,
                                 pas_by_shrub.post.nest, home_by_shrub.post.nest,
                                 stringsAsFactors = FALSE)
diet.summary <- rbind.data.frame(shrub.post.diet, intercept.post.diet, pasture.post.diet ,homestead.post.diet,
                                 pas_by_shrub.post.diet, home_by_shrub.post.diet,
                                 stringsAsFactors = FALSE)
 
write.csv(nest.summary, "./outputs/AppendixS4/AppendixS4nest_trait_effects.2.csv")
write.csv(diet.summary, "./outputs/AppendixS4/AppendixS4diet_trait_effects.2.csv")
###------------------------------------------------------------------------------------
#some graph code from Isabel's paper
#Fig 4: species-specific forest plots with multiple effect sizes
sppEfx_minus_int <- filter(sppEfx, betaName != "Intercept")
forest_plot_spp <- 
ggplot(data = sppEfx_minus_int, aes(x = Species, y = effect)) + 
  geom_errorbar(aes(ymin = LCL, ymax = UCL), width = 0) +
  geom_point(colour="gray20", shape = 21, size = 4, 
             aes(fill = factor(POSNEG))) +                     
  scale_fill_manual(values = c("gray50","orange","blue")) +
  geom_hline(aes(yintercept = 0)) +
  ylab("Median" ~beta~ ", 95% CRI") +         
  scale_x_discrete(limits = rev(sppEfx[order(sppEfx$RobertsNum),]$Species), expand = c(0.02,0.01)) +
  coord_flip() +  
  facet_wrap(~factor(betaName,levels = c("Pasture", "Homestead","Plantation","Shrub cover",
                                         "pas_by_shrub","home_by_shrub")), 
             scales = "fixed", ncol = 5) + ## modify to change column order
  theme_bw()+
  theme(axis.title.x = element_text(vjust = 1.5, size = 18, colour = "black"),
        axis.text.x  = element_text(vjust = 0.5, size = 10, colour = "black", margin = unit(c(0.1,0.1,0.1,0.1), "cm")),
        axis.title.y = element_blank(),
        axis.text.y  = element_text(size=12, colour = "black", margin = unit(c(0.3,0.3,0.3,0.3), "cm")),
        axis.ticks.length = unit(-0.05, "cm"),
        strip.background = element_rect(fill = "gray95"),
        strip.text.x = element_text(size = 12, face = "bold"),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(), 
        legend.position = "none")

#nesting substrate and diet fig
#nest effect
pd <- position_dodge(0.5)
#nest substrates
nest.summary_minus_int <- filter(nest.summary, betaName != "Intercept")
nest.fig <- 
ggplot(data = nest.summary, aes(x = betaName, y = effect, group = nest)) +
  geom_errorbar(aes(ymin = LCL, ymax = UCL), width = 0, size = 0.3, position = pd) +
  geom_point(colour = "gray20", shape = 21, size = 3, position = pd, aes(fill = factor(nest))) + 
  scale_fill_manual(values = c("white", "yellow", "gray50", "blue"), name = "Nesting substrate")+
  geom_hline(aes(yintercept = 0)) +
  ylab("Median" ~beta~ ", 95% CRI") +  
  ylim(c(-15, 15)) +
  scale_x_discrete(limits = rev(levels(factor(nest.summary$betaName, 
                  levels = c("Pasture", "Homestead","Shrub cover",
                             "pas_by_shrub","home_by_shrub") )) )) +
  coord_flip() +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y  = element_text(size = 9, colour = "black", 
                                    margin=unit(c(0.2,0.2,0.2,0.2), "cm")),
        axis.ticks.length=unit(-0.1, "cm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.text=element_text(size = 9), legend.title = element_text(face = "bold",
                                size = 8.6))

diet.summary_minus_int <- filter(diet.summary, betaName != "Intercept")
diet.fig <- 
ggplot(data = diet.summary, aes(x = betaName, y = effect, group = diet)) + 
  geom_errorbar(aes(ymin = LCL, ymax = UCL), width = 0, size = 0.3, position = pd) +
  geom_point(colour = "gray20", shape = 21, size = 3, position = pd, aes(fill = factor(diet))) + 
  scale_fill_manual(values = c("white", "yellow", "gray50","blue","black"), name = "Diet") +
  geom_hline(aes(yintercept = 0)) +
  ylab("Median" ~beta~ ", 95% CRI") +  
  ylim(c(-15, 15)) +
  scale_x_discrete(limits = rev(levels(factor(nest.summary$betaName, 
                   levels = c("Pasture", "Homestead", "Shrub cover",
                              "pas_by_shrub","home_by_shrub") )) )) +
  coord_flip() +
  theme_bw() +
  theme(axis.title.x = element_text(face = "bold", vjust = 0.3, size = 11,
                                    colour = "black"),
        axis.text.x  = element_text(vjust = 0.5, size = 9, colour = "black",
                                    margin = unit(c(0.2,0.2,0.2,0.2), "cm")),
        axis.title.y = element_blank(), 
        axis.text.y  = element_text(size = 9, colour = "black", 
                                    margin = unit(c(0.2,0.2,0.2,0.2), "cm")),
        axis.ticks.length = unit(-0.15, "cm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.text = element_text(size = 9), legend.title = element_text(face = "bold", 
                      size = 10))

##--------------------------------------------------------------------
##               Effects by mass / "wing loading"
##--------------------------------------------------------------------
## filter sppEfx by type of beta
pooled_mass <-  merge(sppEfx, TraitData)
shrub_mass <- merge(sppEfx, TraitData) %>% 
  filter(betaName == "Shrub cover")
intercept_mass <- merge(sppEfx, TraitData) %>% 
  filter(betaName == "Intercept")
past_mass <- merge(sppEfx, TraitData) %>% 
  filter(betaName == "Pasture")
home_mass <- merge(sppEfx, TraitData) %>% 
  filter(betaName == "Homestead")
pas_by_shrub_mass <- merge(sppEfx, TraitData) %>% 
  filter(betaName == "pas_by_shrub")
home_by_shrub_mass <- merge(sppEfx, TraitData) %>% 
  filter(betaName == "home_by_shrub")


## Plot species responses in order of log(mass)
##mass.fig.shrub <- 
ggplot(data = shrub_mass, aes(x = log(Mass), y = effect)) + 
  geom_errorbar(aes(ymin = LCL, ymax = UCL), width = 0) +
  geom_point(colour="gray20", shape = 21, size = 4, aes(fill = factor(POSNEG))) +                     
  scale_fill_manual(values = c("gray50","orange","blue"))+
  scale_x_continuous(breaks = c(2,3,4,5), labels = trunc (c(exp(2), 
                                                            exp(3), exp(4), exp(5) ))) +
  theme_bw()+
  ggtitle("Shrub cover") +
  theme(axis.title.x = element_blank(),
        axis.text.x  = element_text(vjust = 0.5,size = 10,colour = "black", 
                                    margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm")),
        axis.title.y = element_blank(),
        axis.text.y  = element_text(size=12, colour = "black", margin = unit(c(0.3, 0.3, 0.3, 0.3),
                                                                             "cm")),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(), 
        legend.position = "none")

##mass.fig.past <- 
ggplot(data = past_mass, aes(x = log(Mass), y = effect)) + 
  geom_errorbar(aes(ymin = LCL, ymax = UCL), width = 0) +
  geom_point(colour="gray20", shape=21, size = 4, 
             aes(fill = factor(POSNEG))) +                     
  scale_fill_manual(values = c("gray50","orange","blue"))+
  ylim(-30, 30)+
  theme_bw()+
  ggtitle("Pastures") +
  theme(
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(), 
    legend.position = "none")


##mass.fig.home <- 
ggplot(data = home_mass, aes(x = log(Mass), y = effect)) + 
  geom_errorbar(aes(ymin = LCL, ymax = UCL), width = 0) +
  geom_point(colour="gray20", shape=21, size = 4, 
             aes(fill = factor(POSNEG))) +                     
  scale_fill_manual(values = c("gray50","orange","blue"))+
  ylim(-30, 30)+
  theme_bw()+
  ggtitle("Homesteads") +
  theme(element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(), 
        axis.text.y  = element_text(size=12, colour = "black", 
                                    margin = unit(c(0.3, 0.3, 0.3, 0.3), "cm")),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(), 
        legend.position = "none")

##mass.fig.pas_by_shrub <- 
ggplot(data = pas_by_shrub_mass, aes(x = log(Mass), y = effect)) + 
  geom_errorbar(aes(ymin = LCL, ymax = UCL), width = 0) +
  geom_point(colour="gray20", shape=21, size = 4, 
             aes(fill = factor(POSNEG))) +                     
  scale_fill_manual(values = c("gray50","orange","blue"))+
  ylim(-30, 30)+
  theme_bw()+
  ggtitle("pas_by_shrub") +
  theme(element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(), 
        axis.text.y  = element_text(size=12, colour = "black", 
                                    margin = unit(c(0.3, 0.3, 0.3, 0.3), "cm")),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(), 
        legend.position = "none")

##mass.fig.home_by_shrub <- 
ggplot(data = home_by_shrub_mass, aes(x = log(Mass), y = effect)) + 
  geom_errorbar(aes(ymin = LCL, ymax = UCL), width = 0) +
  geom_point(colour="gray20", shape=21, size = 4, 
             aes(fill = factor(POSNEG))) +                     
  scale_fill_manual(values = c("gray50","orange","blue"))+
  ylim(-30, 30)+
  theme_bw()+
  ggtitle("home_by_shrub") +
  theme(element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(), 
        axis.text.y  = element_text(size=12, colour = "black", 
                                    margin = unit(c(0.3, 0.3, 0.3, 0.3), "cm")),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(), 
        legend.position = "none")

## plot species responses in order of pseudo_loading
##ps.load.fig.shrub <- 
ggplot(data = shrub_mass, aes(x = pseudo_loading, y = effect)) + 
  geom_errorbar(aes(ymin = LCL, ymax = UCL), width = 0) +
  geom_point(colour="gray20", shape = 21, size = 4, aes(fill = factor(POSNEG))) +                     
  scale_fill_manual(values = c("gray50","orange","blue"))+
  ylim(-30, 30)+
  theme_bw()+
  ggtitle("Shrub cover") +
  theme(axis.title.x = element_blank(),
        axis.text.x  = element_text(vjust = 0.5,size = 10,colour = "black", 
                                    margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm")),
        axis.title.y = element_blank(),
        axis.text.y  = element_text(size=12, colour = "black", 
                                    margin = unit(c(0.3, 0.3, 0.3, 0.3), "cm")),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(), 
        legend.position = "none")


##ps.load.fig.past <- 
ggplot(data = past_mass, aes(x = pseudo_loading, y = effect)) + 
  geom_errorbar(aes(ymin = LCL, ymax = UCL), width = 0) +
  geom_point(colour="gray20", shape=21, size = 4, 
             aes(fill = factor(POSNEG))) +                     
  scale_fill_manual(values = c("gray50","orange","blue"))+
  ylim(-30, 30)+
  theme_bw()+
  ggtitle("Pastures") +
  theme(
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(), 
    legend.position = "none")


##ps.load.fig.home <- 
ggplot(data = home_mass, aes(x = pseudo_loading, y = effect)) + 
  geom_errorbar(aes(ymin = LCL, ymax = UCL), width = 0) +
  geom_point(colour="gray20", shape=21, size = 4, 
             aes(fill = factor(POSNEG))) +                     
  scale_fill_manual(values = c("gray50","orange","blue"))+
  ylim(-30, 30)+
  theme_bw()+
  ggtitle("Homesteads") +
  theme(element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(), 
        axis.text.y  = element_text(size=12, colour = "black", 
                                    margin = unit(c(0.3, 0.3, 0.3, 0.3), "cm")),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(), 
        legend.position = "none")

##ps.load.fig.pas_by_shrub <- 
ggplot(data = pas_by_shrub_mass, aes(x = pseudo_loading, y = effect)) + 
  geom_errorbar(aes(ymin = LCL, ymax = UCL), width = 0) +
  geom_point(colour="gray20", shape=21, size = 4, 
             aes(fill = factor(POSNEG))) +                     
  scale_fill_manual(values = c("gray50","orange","blue"))+
  ylim(-30, 30)+
  theme_bw()+
  ggtitle("pas_by_shrub") +
  theme(element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(), 
        axis.text.y  = element_text(size=12, colour = "black", 
                                    margin = unit(c(0.3, 0.3, 0.3, 0.3), "cm")),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(), 
        legend.position = "none")

##ps.load.fig.home_by_shrub <- 
ggplot(data = home_by_shrub_mass, aes(x = pseudo_loading, y = effect)) + 
  geom_errorbar(aes(ymin = LCL, ymax = UCL), width = 0) +
  geom_point(colour="gray20", shape=21, size = 4, 
             aes(fill = factor(POSNEG))) +                     
  scale_fill_manual(values = c("gray50","orange","blue"))+
  ylim(-30, 30)+
  theme_bw()+
  ggtitle("home_by_shrub") +
  theme(element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(), 
        axis.text.y  = element_text(size=12, colour = "black", 
                                    margin = unit(c(0.3, 0.3, 0.3, 0.3), "cm")),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(), 
        legend.position = "none")

## compare traits of species exhibiting interactions to species pool 
TraitDataSubset <- TraitData[c(1,8,12,20,22,25,29,31,35,36,37,39,41,43,44,45,47),]
sum(TraitData$Diet == "invertebrates")/48
sum(TraitDataSubset$Diet == "invertebrates")/dim(TraitDataSubset)[[2]]
fisher.test(rbind(c(31,48-31), c(12,17-12)))

TraitDataSubset <- TraitData[c(1,8,12,20,22,25,29,31,35,36,37,39,41,43,44,45,47),]
sum(TraitData$Diet == "seeds")/48
sum(TraitDataSubset$Diet == "seeds")/dim(TraitDataSubset)[[2]]
fisher.test(rbind(c(8,48-8), c(3,17-3)))

TraitDataSubset <- TraitData[c(1,8,12,20,22,25,29,31,35,36,37,39,41,43,44,45,47),]
sum(TraitData$Diet == "fruit")/48
sum(TraitDataSubset$Diet == "fruit")/dim(TraitDataSubset)[[2]]
fisher.test(rbind(c(7,48-7), c(1,17-1)))

TraitDataSubset <- TraitData[c(1,8,12,20,22,25,29,31,35,36,37,39,41,43,44,45,47),]
sum(TraitData$Diet == "nectar")/48
sum(TraitDataSubset$Diet == "nectar")/dim(TraitDataSubset)[[2]]
fisher.test(rbind(c(2,48-2), c(1,17-1)))

sum(TraitData$ShrubNest == 1)/48
sum(TraitDataSubset$ShrubNest == 1)/dim(TraitDataSubset)[[2]]
fisher.test(rbind(c(29,48-29), c(10,17-10)))

sum(TraitData$TreeNest == 1)/48
sum(TraitDataSubset$TreeNest == 1)/dim(TraitDataSubset)[[2]]
fisher.test(rbind(c(30,48-30), c(11,17-11)))

sum(TraitData$CavityNest == 1)/48
sum(TraitDataSubset$CavityNest == 1)/dim(TraitDataSubset)[[2]]
fisher.test(rbind(c(1,48-1), c(0,17-0)))

sum(TraitData$GroundNest == 1)/48
sum(TraitDataSubset$GroundNest == 1)/dim(TraitDataSubset)[[2]]
fisher.test(rbind(c(9,48-9), c(3,17-3)))

sum(TraitData$GrassNest == 1)/48
sum(TraitDataSubset$GrassNest == 1)/dim(TraitDataSubset)[[2]]
fisher.test(rbind(c(6,48-6), c(2,17-2)))

t.test(TraitData$pseudo_loading, TraitDataSubset$pseudo_loading)
t.test(TraitData$Mass, TraitDataSubset$Mass)
