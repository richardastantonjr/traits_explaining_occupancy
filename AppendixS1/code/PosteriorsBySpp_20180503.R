
# Summary of what this script does ----------------------------------------
## Extract posterior MCMC chains and process to get posterior 95% CRIs for suites of species with
## common traits, regress trait values on shrub cover and land-use effects, and display results
## in figures

# functions ---------------------------------------------------------------
# return mean and 90% CRI
meanAnd90CRI <- function(x) {
  c(mean <- mean(x), quantile(x, c(0.05, 0.95))) 
}

# return mean and 95% CRI
meanAndCRI <- function(x) {
  c(mean <- mean(x), quantile(x, c(0.025, 0.975))) 
}

# libraries ----------------------------------------------------------------------------
library(tidyverse)
library(reshape2) 
library(cowplot)
library(gridExtra)
library(MCMCvis)

# load data -----------------------------------------------------------------------------
TraitData <- read.csv("./data/TraitData.csv")
TraitData[, 1] <- as.factor(TraitData[,1])

## load models with quadratic time of day effects and effects parameterization 
## Prior version used means parameterization and linear time of day effects.
load(("./Occupancy_analysis_Swazi_birds_Bayesian_quadTime_fixed.RData"))

# subset and process data  -----------------------------------------------------------------------
TraitData <- filter(TraitData,TraitData$Species %in% dimnames(AbunHists.subset)[[3]]) 
TraitData <- TraitData[,c(1:6, 9:14, 23, 24, 26)]
TraitData[is.na(TraitData)] <- 0 
TraitData <- TraitData[order(TraitData$RobertsNum), ] ## This doesn't accomplish anything
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

## Determine the # of spp observed in each land use
observed <- matrix(0, nrow = 371, ncol = 48)
for (i in 1:371) {
  for (j in 1:48) {
    observed[i,j] <- max(AbunHists.subset[i,,j])
  }
}
observed[is.na(observed)] <- 0 ## some points have 0 species observed

## Substract any species not observed in a land use from the total
48 - sum(colSums(observed[which(siteCovs$LandUse == "Protected"),]) == 0)
48 - sum(colSums(observed[which(siteCovs$LandUse == "ComPasture"),]) == 0)
48 - sum(colSums(observed[which(siteCovs$LandUse == "Homestead"),]) == 0)
48 - sum(colSums(observed[which(siteCovs$LandUse == "SugarEstate"),]) == 0)

## Number of species observed at each point by land use
summary(rowSums(observed[which(siteCovs$LandUse == "Protected"),]))
summary(rowSums(observed[which(siteCovs$LandUse == "ComPasture"),]))
summary(rowSums(observed[which(siteCovs$LandUse == "Homestead"),]))
summary(rowSums(observed[which(siteCovs$LandUse == "SugarEstate"),]))

## mean # of species observed at a point by land use
mean(rowSums(observed[which(siteCovs$LandUse == "Protected"),]))   ## 13.65 +/- 0.40 SE
mean(rowSums(observed[which(siteCovs$LandUse == "ComPasture"),]))  ## 13.05 +/- 0.44 SE
mean(rowSums(observed[which(siteCovs$LandUse == "Homestead"),]))   ## 10.77 +/- 0.52 SE
mean(rowSums(observed[which(siteCovs$LandUse == "SugarEstate"),])) ##  3.92 +/- 0.61 SE

## Standard deviations
sd(rowSums(observed[which(siteCovs$LandUse == "Protected"),]))
sd(rowSums(observed[which(siteCovs$LandUse == "ComPasture"),]))
sd(rowSums(observed[which(siteCovs$LandUse == "Homestead"),]))
sd(rowSums(observed[which(siteCovs$LandUse == "SugarEstate"),]))

## Standard errors
sd(rowSums(observed[which(siteCovs$LandUse == "Protected"),])) / 
    sqrt(length(rowSums(observed[which(siteCovs$LandUse == "Protected"),])))
sd(rowSums(observed[which(siteCovs$LandUse == "ComPasture"),])) /
    sqrt(length(rowSums(observed[which(siteCovs$LandUse == "ComPasture"),])))
sd(rowSums(observed[which(siteCovs$LandUse == "Homestead"),])) /
    sqrt(length(rowSums(observed[which(siteCovs$LandUse == "Homestead"),])))
sd(rowSums(observed[which(siteCovs$LandUse == "SugarEstate"),])) /
    sqrt(length(rowSums(observed[which(siteCovs$LandUse == "SugarEstate"),])))

## remove extraneous material
rm(list=ls()[! ls() %in% c("Outputs_long","TraitData", "meanAndCRI", "meanAnd90CRI")])
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

## summarize results for the community of common species
pooledInterceptCRI <- meanAndCRI(rowMeans(intercept.post))
pooledHomesteadCRI <-  meanAndCRI(rowMeans(homestead.post))
pooledPastureCRI <-  meanAndCRI(rowMeans(pasture.post))
pooledSugarEstateCRI <- meanAndCRI(rowMeans(sugarEstate.post))
pooledShrubCRI <- meanAndCRI(rowMeans(shrub.post))

pooledIntercept90CRI <- meanAnd90CRI(rowMeans(intercept.post))
pooledHomestead90CRI <- meanAnd90CRI(rowMeans(homestead.post))
pooledPasture90CRI <- meanAnd90CRI(rowMeans(pasture.post))
pooledSugarEstate90CRI <- meanAnd90CRI(rowMeans(sugarEstate.post)) 
pooledShrub90CRI <- meanAnd90CRI(rowMeans(shrub.post))  

res_pooled_95 <- rbind(pooledInterceptCRI,
                       pooledShrubCRI,
                       pooledPastureCRI, 
                       pooledHomesteadCRI,
                       pooledSugarEstateCRI)

res_pooled_90 <- rbind(pooledIntercept90CRI,
                       pooledShrub90CRI,
                       pooledPasture90CRI, 
                       pooledHomestead90CRI,
                       pooledSugarEstate90CRI)


res_pooled <- cbind(res_pooled_95, res_pooled_90)
rownames(res_pooled) <- c("Intercept", "Shrub cover", "Pasture","Homestead","Plantation")
colnames(res_pooled) <- c("Mean", "95% LCRI", "95% UCRI", "Mean", "90% LCRI", "90% UCRI")
res_pooled <- res_pooled[,-4] ## Remove duplicated mean
write.csv(res_pooled, "./outputs/AppendixS3/AppendixS3.1.csv")

## Overall effects of each covariate
plot(density(rowMeans(shrub.post)), xlim = c(-10, 10), xlab = "Beta estimate", 
             main = "", lwd = 2, las = 1)
lines(density(rowMeans(pasture.post)), lty = 2, col = "blue", lwd = 2)
lines(density(rowMeans(homestead.post)), lty = 3, col = "grey", lwd = 2)
lines(density(rowMeans(sugarEstate.post)), lty = 4, col = "orange", lwd = 2)
abline(v=0)
legend("topleft", c("shrub","pasture","homestead","plantation"),
       fill = c("black", "blue", "grey","orange"), bty = "n")

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

sugarEstate.post$nest.shrub <- rowMeans(sugarEstate.post[,names(sugarEstate.post) %in% nest.shrub$Species])
sugarEstate.post$nest.tree <- rowMeans(sugarEstate.post[,names(sugarEstate.post) %in% nest.tree$Species])
sugarEstate.post$nest.cavity <- rowMeans(sugarEstate.post[,names(sugarEstate.post) %in% nest.cavity$Species])
sugarEstate.post$nest.grass <- rowMeans(sugarEstate.post[,names(sugarEstate.post) %in% nest.grass$Species])

sugarEstate.post$predators <- rowMeans(sugarEstate.post[,names(sugarEstate.post) %in% predators$Species])
sugarEstate.post$diet.invert <- rowMeans(sugarEstate.post[,names(sugarEstate.post) %in% diet.invert$Species])
sugarEstate.post$diet.fruit <- rowMeans(sugarEstate.post[,names(sugarEstate.post) %in% diet.fruit$Species])
sugarEstate.post$diet.seeds <- rowMeans(sugarEstate.post[,names(sugarEstate.post) %in% diet.seed$Species])
sugarEstate.post$diet.nectar <- rowMeans(sugarEstate.post[,names(sugarEstate.post) %in% diet.nectar$Species])

## two examples of density plots; Forest plots below.
plot(density(shrub.post$predators), xlim = c(-15,5), col = "black", lwd= 3, main = "",
     las = 1)
lines(density(pasture.post$predators), col = "blue", lty = 2, lwd = 3)
lines(density(homestead.post$predators), col = "grey", lty = 3, lwd = 3)
lines(density(sugarEstate.post$predators), col = "orange", lty = 4, lwd = 3)
abline(v = 0)
legend("topleft", c("shrub","pasture","homestead","plantation"),
       fill = c("black", "blue", "grey","orange"), bty = "n")

plot(density(shrub.post$diet.invert), xlim = c(-15,5), col = "black", lwd= 3, main = "",
     las = 1)
lines(density(pasture.post$diet.invert), col = "blue", lty = 2, lwd = 3)
lines(density(homestead.post$diet.invert), col = "grey", lty = 3, lwd = 3)
lines(density(sugarEstate.post$diet.invert), col = "orange", lty = 4, lwd = 3)
abline(v = 0)
legend("topleft", c("shrub","pasture","homestead","plantation"),
       fill = c("black", "blue", "grey","orange"), bty = "n")

## create a data frame of effect sizes and CRIs for each species for plotting species names
Species <- levels(as.factor(Outputs_long$Species))

## create a data frame of one effect and CRIs for all species  
shrubSppEfx <- apply(shrub.post, 2, meanAndCRI )
shrubSppEfx <- as.data.frame(t(shrubSppEfx[, 1:48]))
shrubSppEfx <- cbind(Species, shrubSppEfx)
colnames(shrubSppEfx) <- c("Species", "effect", "LCL", "UCL")
shrubSppEfx$betaName <- "Shrub cover"
rownames(shrubSppEfx)<- NULL

## repeat for the remaining effects
homesteadSppEfx <- apply(homestead.post, 2, meanAndCRI )
homesteadSppEfx <- as.data.frame(t(homesteadSppEfx[, 1:48]))
homesteadSppEfx <- cbind(Species, homesteadSppEfx)
colnames(homesteadSppEfx) <- c("Species", "effect", "LCL", "UCL")
homesteadSppEfx$betaName <- "Homestead"
rownames(homesteadSppEfx)<- NULL

pastureSppEfx <- apply(pasture.post, 2, meanAndCRI )
pastureSppEfx <- as.data.frame(t(pastureSppEfx[, 1:48]))
pastureSppEfx <- cbind(Species, pastureSppEfx)
colnames(pastureSppEfx) <- c("Species", "effect", "LCL", "UCL")
pastureSppEfx$betaName <- "Pasture"
rownames(pastureSppEfx)<- NULL

interceptSppEfx <- apply(intercept.post, 2, meanAndCRI )
interceptSppEfx <- as.data.frame(t(interceptSppEfx[, 1:48]))
interceptSppEfx <- cbind(Species, interceptSppEfx)
colnames(interceptSppEfx) <- c("Species", "effect", "LCL", "UCL")
interceptSppEfx$betaName <- "Intercept"
rownames(interceptSppEfx)<- NULL

sugarSppEfx <- apply(sugarEstate.post, 2, meanAndCRI )
sugarSppEfx <- as.data.frame(t(sugarSppEfx[, 1:48]))
sugarSppEfx <- cbind(Species, sugarSppEfx)
colnames(sugarSppEfx) <- c("Species", "effect", "LCL", "UCL")
sugarSppEfx$betaName <- "Plantation"
rownames(sugarSppEfx)<- NULL

## bind rows to create one frame for plotting 
sppEfx <- rbind.data.frame(shrubSppEfx, interceptSppEfx, pastureSppEfx, homesteadSppEfx, sugarSppEfx )
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
intEfxAndTraits <- merge(sppEfx, TraitData) %>% 
  filter(betaName == "Protected")
pastEfxAndTraits <- merge(sppEfx, TraitData) %>% 
  filter(betaName == "Pasture")
homeEfxAndTraits <- merge(sppEfx, TraitData) %>% 
  filter(betaName == "Homestead")
sugarEfxAndTraits <- merge(sppEfx, TraitData) %>% 
  filter(betaName == "Plantation")

## full posterior samples. Break into 2000 data sets for each beta,
## each representing a posterior draw and associated continuous trait values.
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

## Intercepts
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
##------------------

## Plantation effect
betaL4_samp <- list()
for (m in 1:2000) {
  betaL4_samp[[m]] <- Outputs_long$`beta.l[4]`[c(seq(m,96000,2000))]
}

betaL4_slopes_pseudo_load <- NULL
betaL4_slopes_mass <- NULL
for (mm in 1:2000) {
  betaL4_slopes_pseudo_load[mm] <- coef(lm(betaL4_samp[[mm]] ~ pseudo_loadings[[mm]]))[2]
  betaL4_slopes_mass[mm] <- coef(lm(betaL4_samp[[mm]] ~ masses[[mm]]))[2]
} 
##--------------------------

## bundle the results
slopes_95 <- rbind.data.frame(meanAndCRI(beta1_slopes_mass),
                                   meanAndCRI(beta_slopes_mass),
                                   meanAndCRI(betaL2_slopes_mass),
                                   meanAndCRI(betaL3_slopes_mass),
                                   meanAndCRI(betaL4_slopes_mass),
                                              
                                   meanAndCRI(beta1_slopes_pseudo_load),
                                   meanAndCRI(beta_slopes_pseudo_load),
                                   meanAndCRI(betaL2_slopes_pseudo_load),
                                   meanAndCRI(betaL3_slopes_pseudo_load),
                                   meanAndCRI(betaL4_slopes_pseudo_load)
                                   )
colnames(slopes_95) <- c("effect", "LCL", "UCL")
slopes_95$betaName <- c("Shrub cover","Intercept","Pasture","Homestead","Plantation",
                        "Shrub cover","Intercept","Pasture","Homestead","Plantation")
slopes_95$measure <- c("mass","mass","mass","mass","mass", "loading","loading","loading",
                      "loading","loading")


slopes_90 <- rbind.data.frame(meanAnd90CRI(beta1_slopes_mass),
                              meanAnd90CRI(beta_slopes_mass),
                              meanAnd90CRI(betaL2_slopes_mass),
                              meanAnd90CRI(betaL3_slopes_mass),
                              meanAnd90CRI(betaL4_slopes_mass),
                              
                              meanAnd90CRI(beta1_slopes_pseudo_load),
                              meanAnd90CRI(beta_slopes_pseudo_load),
                              meanAnd90CRI(betaL2_slopes_pseudo_load),
                              meanAnd90CRI(betaL3_slopes_pseudo_load),
                              meanAnd90CRI(betaL4_slopes_pseudo_load)
)
colnames(slopes_90) <- c("median", "LCL", "UCL")
slopes_90$betaName <- c("Shrub","Intercept","Pasture","Homestead","Plantation",
                        "Shrub","Intercept","Pasture","Homestead","Plantation")
slopes_90$measure <- c("mass","mass","mass","mass","mass", "loading","loading","loading",
                      "loading","loading")

cont_trait_slopes <- cbind.data.frame(slopes_95, slopes_90)
cont_trait_slopes <- cont_trait_slopes[,c(1:3, 7:10)]
colnames(cont_trait_slopes)[2:5] <- c("95% LCL", "95% UCL", "90%LCL", "90% UCL")
write.csv(cont_trait_slopes, "./outputs/AppendixS4/AppendixS4cont_trait_slopes.1.csv")
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
intercept.post.nest$betaName <- "intercept"
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

r17 <- meanAndCRI (sugarEstate.post$nest.shrub)
r18 <- meanAndCRI (sugarEstate.post$nest.tree)
r19 <- meanAndCRI (sugarEstate.post$nest.cavity)
r20 <- meanAndCRI (sugarEstate.post$nest.grass)
sugarEstate.post.nest <- rbind.data.frame(r17,r18,r19,r20)
colnames(sugarEstate.post.nest) <- c("effect", "LCL","UCL")
sugarEstate.post.nest$betaName <- "Plantation"
sugarEstate.post.nest$nest <- c("shrubs", "trees", "cavities","grass")


## Diets
r21 <- meanAndCRI(shrub.post$predators)
r22 <- meanAndCRI(shrub.post$diet.invert)
r23 <- meanAndCRI(shrub.post$diet.fruit)
r24 <- meanAndCRI(shrub.post$diet.seeds)
r25 <- meanAndCRI(shrub.post$diet.nectar)
shrub.post.diet <- rbind.data.frame(r21,r22,r23,r24,r25)
colnames(shrub.post.diet) <- c("effect", "LCL","UCL")
shrub.post.diet$betaName <- "Shrub cover"
shrub.post.diet$diet <- c("predatory","invertebrates","fruit","seeds","nectar")

r26 <- meanAndCRI (intercept.post$predators)
r27 <- meanAndCRI (intercept.post$diet.invert)
r28 <- meanAndCRI (intercept.post$diet.fruit)
r29 <- meanAndCRI (intercept.post$diet.seeds)
r30 <- meanAndCRI (intercept.post$diet.nectar)
intercept.post.diet <- rbind.data.frame(r26,r27,r28,r29,r30)
colnames(intercept.post.diet) <- c("effect", "LCL","UCL")
intercept.post.diet$betaName <- "intercept"
intercept.post.diet$diet <- c("predatory","invertebrates","fruit","seeds","nectar")

r31 <- meanAndCRI (pasture.post$predators)
r32 <- meanAndCRI (pasture.post$diet.invert)
r33 <- meanAndCRI (pasture.post$diet.fruit)
r34 <- meanAndCRI (pasture.post$diet.seeds)
r35 <- meanAndCRI (pasture.post$diet.nectar)
pasture.post.diet <- rbind.data.frame(r31,r32,r33,r34,r35)
colnames(pasture.post.diet) <- c("effect", "LCL","UCL")
pasture.post.diet$betaName <- "Pasture"
pasture.post.diet$diet <- c("predatory","invertebrates","fruit","seeds","nectar")

r36 <- meanAndCRI (homestead.post$predators)
r37 <- meanAndCRI (homestead.post$diet.invert)
r38 <- meanAndCRI (homestead.post$diet.fruit)
r39 <- meanAndCRI (homestead.post$diet.seeds)
r40 <- meanAndCRI (homestead.post$diet.nectar)
homestead.post.diet <- rbind.data.frame(r36,r37,r38,r39,r40)
colnames(homestead.post.diet) <- c("effect", "LCL","UCL")
homestead.post.diet$betaName <- "Homestead"
homestead.post.diet$diet <- c("predatory","invertebrates","fruit","seeds","nectar")

r41 <- meanAndCRI (sugarEstate.post$predators)
r42 <- meanAndCRI (sugarEstate.post$diet.invert)
r43 <- meanAndCRI (sugarEstate.post$diet.fruit)
r44 <- meanAndCRI (sugarEstate.post$diet.seeds)
r45 <- meanAndCRI (sugarEstate.post$diet.nectar)
sugarEstate.post.diet <- rbind.data.frame(r41,r42,r42,r44,r45)
colnames(sugarEstate.post.diet) <- c("effect", "LCL","UCL")
sugarEstate.post.diet$betaName <- "Plantation"
sugarEstate.post.diet$diet <- c("predatory","invertebrates","fruit","seeds","nectar")

nest.summary <- rbind.data.frame(shrub.post.nest, intercept.post.nest, pasture.post.nest, homestead.post.nest,
                                 sugarEstate.post.nest, stringsAsFactors = FALSE)
diet.summary <- rbind.data.frame(shrub.post.diet, intercept.post.diet, pasture.post.diet ,homestead.post.diet,
                                 sugarEstate.post.diet, stringsAsFactors = FALSE)

write.csv(nest.summary, "./outputs/AppendixS4/AppendixS4nest_trait_effects.1.csv")
write.csv(diet.summary, "./outputs/AppendixS4/AppendixS4diet_trait_effects.1.csv")
####-----------------------------------------------------------------------------------
###------------------------------------------------------------------------------------
#some graph code from Isabel's paper
#Fig 4: species-specific forest plots with multiple effect sizes
sppEfx_minus_int <- filter(sppEfx, betaName != "Intercept")
forest_plot_spp <- ggplot(data = sppEfx_minus_int, aes(x = Species, y = effect)) + 
  geom_errorbar(aes(ymin = LCL, ymax = UCL), width = 0) +
  geom_point(colour="gray20", shape = 21, size = 4, 
             aes(fill = factor(POSNEG))) +                     
  scale_fill_manual(values = c("gray50","orange","blue")) +
  geom_hline(aes(yintercept = 0)) +
  ylab("Median" ~beta~ ", 95% CRI") +         
  scale_x_discrete(limits = rev(sppEfx[order(sppEfx$RobertsNum),]$Species), expand = c(0.02,0.01)) +
  coord_flip() +  
  facet_wrap(~factor(betaName,levels = c("Pasture", "Homestead","Plantation","Shrub cover")), 
        scales = "fixed", ncol = 4) + ## modify to change column order
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

sppEfx_subset <- filter(sppEfx_minus_int, Species == "Golden-tailed Woodpecker" | Species == "Black-headed Oriole"
                        | Species == "Yellow-breasted Apalis" | Species == "Spotted Flycatcher" | 
                          Species == "Black-backed Puffback")

sppEfx_subset <- filter(sppEfx_subset, betaName == "Homestead" | betaName == "Shrub cover")

forest_plot_mitigation_cands <- ggplot(data = sppEfx_subset, aes(x = Species, y = effect)) + 
  geom_errorbar(aes(ymin = LCL, ymax = UCL), width = 0) +
  geom_point(colour = "gray20", shape = 21, size = 4, 
             aes(fill = factor(POSNEG))) +                     
  scale_fill_manual(values = c("blue", "orange")) +
  geom_hline(aes(yintercept = 0)) +
  ylab("Median" ~beta~ ", 95% CRI") +         
  scale_x_discrete(limits = rev(sppEfx_subset[order(sppEfx_subset$RobertsNum),]$Species), expand = c(0.02,0.01)) +
  coord_flip() +  
  facet_wrap(~factor(betaName,levels = c("Protected", "Pasture", "Homestead","Plantation","Shrub cover")), 
             scales = "fixed", ncol = 5) +
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

sppEfx_subset2 <- filter(sppEfx_minus_int, Species == "Barn Swallow" | 
                         Species == "Southern Fiscal" | Species == "Southern Grey-headed Sparrow" | 
                         Species == "Black-crowned Tchagra")

sppEfx_subset2 <- filter(sppEfx_subset2, betaName == "Homestead" | betaName == "Shrub cover")

forest_plot_shr_neg <- ggplot(data = sppEfx_subset2, aes(x = Species, y = effect)) + 
  geom_errorbar(aes(ymin = LCL, ymax = UCL), width = 0) +
  geom_point(colour = "gray20", shape = 21, size = 4, 
             aes(fill = factor(POSNEG))) +                     
  scale_fill_manual(values = c("grey", "orange", "blue")) +
  geom_hline(aes(yintercept = 0)) +
  ylab("Median" ~beta~ ", 95% CRI") +         
  coord_flip() +  
  facet_wrap(~factor(betaName,levels = c("Pasture", "Homestead","Plantation","Shrub cover")), 
             scales = "fixed", ncol = 5) +
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

## plotting detail
pd <- position_dodge(0.5)
#nesting substrate and diet fig
#nest substrates
nest.summary_minus_int <- filter(nest.summary, betaName != "Intercept")
nest.fig <- ggplot(data = nest.summary_minus_int, aes(x = betaName, y = effect, group = nest)) +
  geom_errorbar(aes(ymin = LCL, ymax = UCL), width = 0, size = 0.3, position = pd) +
  geom_point(colour = "gray20", shape = 21, size = 3, position = pd, aes(fill = factor(nest))) + 
  scale_fill_manual(values = c("white", "yellow", "gray50", "blue"), name = "Nesting substrate")+
  geom_hline(aes(yintercept = 0)) +
  ggtitle("Nest substrate effects")+
  ylab("Median" ~beta~ ", 95% CRI") +  
  ylim(c(-15, 15)) +
  scale_x_discrete(limits = rev(levels(factor(nest.summary$betaName, 
          levels = c("Pasture", "Homestead","Plantation","Shrub cover") )) )) +
  coord_flip() +
  theme_bw() +
  theme(axis.title.x = element_text(face = "bold", vjust = 0.3, size = 11,
                                    colour = "black"),
        axis.text.x  = element_text(vjust = 0.5, size = 9, colour = "black",
                                    margin = unit(c(0.2,0.2,0.2,0.2), "cm")),
        axis.title.y = element_blank(), 
        axis.text.y  = element_text(size = 9, colour = "black", 
                                    margin = unit(c(0.2,0.2,0.2,0.2), "cm")),
        axis.ticks.length=unit(-0.1, "cm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.text=element_text(size = 9), legend.title = element_text(face = "bold",
            size = 8.6))

diet.summary_minus_int <- filter(diet.summary, betaName != "Intercept")
diet.fig <- ggplot(data = diet.summary_minus_int, aes(x = betaName, y = effect, group = diet)) + 
  geom_errorbar(aes(ymin = LCL, ymax = UCL), width = 0, size = 0.3, position = pd) +
  geom_point(colour = "gray20", shape = 21, size = 3, position = pd, aes(fill = factor(diet))) + 
  scale_fill_manual(values = c("white", "yellow", "gray50","blue","black"), name = "Diet") +
  geom_hline(aes(yintercept = 0)) +
  ggtitle("Diet effects")+
  ylab("Median" ~beta~ ", 95% CRI") +  
  ylim(c(-15, 15)) +
  scale_x_discrete(limits = rev(levels(factor(nest.summary$betaName, 
        levels = c("Pasture", "Homestead","Plantation","Shrub cover") )) )) +
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
res_pooled_95 <- as.data.frame(res_pooled_95)
colnames(res_pooled_95) <- c("effect", "LCL", "UCL")
res_pooled_95$betaName <- slopes_95$betaName[c(2,1,3,4,5)]
main_efx <- ggplot(data = res_pooled_95[2:5,], aes(x = betaName, y = effect)) + 
  ggtitle("Main effects")+
  geom_errorbar(aes(ymin = LCL, ymax = UCL), width = 0, size = 0.3, position = pd) +
  geom_point(shape = 21, size = 3, position = pd, fill = "salmon")+
  geom_hline(aes(yintercept = 0)) +
  ylab("Median"~beta~ "95% CRI") +  
  ylim(c(-15, 15)) +
  scale_x_discrete(limits = rev(levels(factor(res_pooled_95$betaName, 
                                              levels = c("Pasture", "Homestead", "Plantation", "Shrub cover") )) )) +
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

reg_on_loading <- ggplot(data = slopes_95[c(6,8,9,10),], aes(x = betaName, y = effect)) + 
  ggtitle(~beta~ "regressed on maneuverability")+
  geom_errorbar(aes(ymin = LCL, ymax = UCL), width = 0, size = 0.3, position = pd) +
  geom_point(shape = 21, size = 3, position = pd, fill = "salmon")+
  geom_hline(aes(yintercept = 0)) +
  ylab("Median slope, 95% CRI") +  
  ylim(c(-12, 12)) +
  scale_x_discrete(limits = rev(levels(factor(slopes_95$betaName, 
              levels = c("Pasture", "Homestead", "Plantation", "Shrub cover") )) )) +
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

reg_on_mass <- ggplot(data = slopes_95[c(1,3,4,5),], aes(x = betaName, y = effect)) + 
  ggtitle(~beta~ "regressed on mass")+
  geom_errorbar(aes(ymin = LCL, ymax = UCL), width = 0, size = 0.3, position = pd) +
  geom_point(shape = 21, size = 3, position = pd, fill = "salmon")+
  geom_hline(aes(yintercept = 0)) +
  ylab("") +  
  ylim(c(-0.1, 0.1)) +
  scale_x_discrete(limits = rev(levels(factor(slopes_95$betaName, 
                            levels = c("Pasture", "Homestead", "Plantation", "Shrub cover") )) )) +
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

## filter sppEfx by type of beta
pooled_mass <-  merge(sppEfx, TraitData)
shrub_mass <- merge(sppEfx, TraitData) %>% 
  filter(betaName == "Shrub cover")
int_mass <- merge(sppEfx, TraitData) %>% 
  filter(betaName == "Intercept")
past_mass <- merge(sppEfx, TraitData) %>% 
  filter(betaName == "Pasture")
home_mass <- merge(sppEfx, TraitData) %>% 
  filter(betaName == "Homestead")
plant_mass <- merge(sppEfx, TraitData) %>% 
  filter(betaName == "Plantation")

## Plot species responses in order of log(mass)
mass.fig.shrub <- ggplot(data = shrub_mass, aes(x = log(Mass), y = effect)) + 
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
        axis.text.y  = element_blank(),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(), 
        legend.position = "none")

mass.fig.past <- ggplot(data = past_mass, aes(x = log(Mass), y = effect)) + 
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
    axis.text.y  = element_text(vjust = 0.5, size = 9, colour = "black",
                                margin = unit(c(0.2,0.2,0.2,0.2), "cm")),
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(), 
    legend.position = "none")


mass.fig.home <- ggplot(data = home_mass, aes(x = log(Mass), y = effect)) + 
  geom_errorbar(aes(ymin = LCL, ymax = UCL), width = 0) +
  geom_point(colour="gray20", shape = 21, size = 4, 
             aes(fill = factor(POSNEG))) +                     
  scale_fill_manual(values = c("gray50","orange","blue"))+
  scale_x_continuous(breaks = c(2,3,4,5), labels = trunc (c(exp(2), 
                                                            exp(3), exp(4), exp(5) ))) +
  ylim(-30, 30)+
  theme_bw()+
  ggtitle("Homesteads") +
  theme(element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.y = element_blank(), 
    axis.text.y  = element_blank(), 
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(), 
    legend.position = "none")


mass.fig.plant <- ggplot(data = plant_mass, aes(x = log(Mass), y = effect)) + 
  geom_errorbar(aes(ymin = LCL, ymax = UCL), width = 0) +
  geom_point(colour = "gray20", shape = 21, size = 4, 
             aes(fill = factor(POSNEG))) +                     
  scale_fill_manual(values = c("gray50","orange","blue")) +
  scale_x_continuous(breaks = c(2,3,4,5), labels = trunc (c(exp(2), 
                    exp(3), exp(4), exp(5) ))) +
  ylim(-30, 30)+
  theme_bw()+
  ggtitle("Sugar cane plantation") +
  theme( axis.text.x  = element_text(vjust = 0.5,size = 10,colour = "black", 
        margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm")),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y  = element_text(vjust = 0.5,size = 10,colour = "black", 
                                margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm")),
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(), 
    legend.position = "none")

## Plot species responses in order of wing loading proxy
ps.load.fig.shrub <- ggplot(data = shrub_mass, aes(x = pseudo_loading, y = effect)) + 
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
        axis.text.y  = element_blank(),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(), 
        legend.position = "none")

ps.load.fig.past <- ggplot(data = past_mass, aes(x = pseudo_loading, y = effect)) + 
  geom_errorbar(aes(ymin = LCL, ymax = UCL), width = 0) +
  geom_point(colour="gray20", shape=21, size = 4, 
             aes(fill = factor(POSNEG))) +                     
  scale_fill_manual(values = c("gray50","orange","blue"))+
  ylim(-30, 30)+
  theme_bw()+
  ggtitle("Pastures") +
  theme(
    axis.title.x = element_blank(),
    axis.text.x  =element_blank(),
    axis.title.y = element_blank(),
    axis.text.y  = element_text(vjust = 0.5,size = 10,colour = "black", 
                                margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm")),
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(), 
    legend.position = "none")


ps.load.fig.home <- ggplot(data = home_mass, aes(x = pseudo_loading, y = effect)) + 
  geom_errorbar(aes(ymin = LCL, ymax = UCL), width = 0) +
  geom_point(colour="gray20", shape=21, size = 4, 
             aes(fill = factor(POSNEG))) +                     
  scale_fill_manual(values = c("gray50","orange","blue"))+
  ylim(-30, 30)+
  theme_bw()+
  ggtitle("Homesteads") +
  theme(element_blank(),
        axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.title.y = element_blank(), 
        axis.text.y  = element_blank(),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(), 
        legend.position = "none")


ps.load.fig.plant <- ggplot(data = plant_mass, aes(x = pseudo_loading, y = effect)) + 
  geom_errorbar(aes(ymin = LCL, ymax = UCL), width = 0) +
  geom_point(colour = "gray20", shape = 21, size = 4, 
             aes(fill = factor(POSNEG))) +                     
  scale_fill_manual(values = c("gray50","orange","blue")) +
  ylim(-30, 30)+
  theme_bw()+
  ggtitle("Sugar cane plantation") +
  theme( axis.text.x  = element_text(vjust = 0.5,size = 10, colour = "black", 
                                     margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm")),
         axis.title.x = element_blank(),
         axis.title.y = element_blank(),
         axis.text.y  = element_text(vjust = 0.5,size = 10, colour = "black", 
                                     margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm")),
         panel.grid.major = element_blank(),  
         panel.grid.minor = element_blank(), 
         legend.position = "none")

# save plots --------------------------------------------------------------
## fig 1 is the study area map with shrub cover inlay, produced in another script
fig2 <- grid.arrange(forest_plot_spp)
fig3 <- plot_grid(main_efx, nest.fig, diet.fig, labels = c("a)", "b)", "c)"), 
                  ncol = 3, rel_widths = c(0.72, 1, 1))
figX <- plot_grid(nest.fig, diet.fig, ncol = 1, labels = c("a)","b)"))
figXX <- grid.arrange(forest_plot_mitigation_cands)
figXXX <- grid.arrange(forest_plot_shr_neg)
figS4_1 <- plot_grid(reg_on_mass, reg_on_loading, ncol = 1, labels = c("a)", "b)"))
figS7_1 <- grid.arrange(mass.fig.past, mass.fig.home, mass.fig.plant, mass.fig.shrub,
                  ncol = 2, bottom = "Mass, g", 
                  #left = textGrob (label = expression(paste("Median", beta, ", 95% CRI")))
                  left = "Median beta, 95% CRI") 
figS7_2 <- grid.arrange(ps.load.fig.past, ps.load.fig.home, ps.load.fig.plant, ps.load.fig.shrub,
         ncol = 2, bottom = "Wing loading (higher is less maneuverable)", left = "Median beta, 95% CRI") ## Edit to Greek

ggsave(file = "./outputs/figure2.tiff", fig2, width = 12, height = 10, units ="in", dpi = 600)
ggsave(file = "./outputs/figure3.tiff", fig3, width = 12, height = 10, units ="in", dpi = 300)
ggsave(file = "./outputs/diet_nest_efx.tiff", figx, width = 5, height = 8, units ="in", dpi = 300)
ggsave(file = "./outputs/continuous_efx.tiff", figS4_1, width = 5, height = 8, units ="in", dpi = 300)
ggsave(file = "./outputs/forest_plot_mitigation_cands.tiff", figXX, width = 10, height = 6, units ="in", dpi = 600)
ggsave(file = "./outputs/forest_plot_shr_neg.tiff", figXXX, width = 10, height = 6, units ="in", dpi = 600)
ggsave(file = "./outputs/mass_efx.tiff", figS4_2_upper, width = 10, height = 8, units ="in", dpi = 300)
ggsave(file = "./outputs/wing_efx.tiff", figS4_2_lower, width = 10, height = 8, units ="in", dpi = 300)

## collect trait results into a sorted csv file, omitting the 90% CRIs
colnames(diet.summary)[5] <- "trait"
colnames(nest.summary)[5] <- "trait"
colnames(slopes_95)[5] <- "trait"
diet.summary <- arrange(diet.summary, trait)
nest.summary <- arrange(nest.summary, trait)
trait_res_95_add_models <- rbind.data.frame(nest.summary, diet.summary, slopes_95)
write.csv(trait_res_95_add_models, "./outputs/trait_res_95_add_models.csv")
