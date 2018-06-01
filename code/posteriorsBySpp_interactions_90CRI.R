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
## Prior version used means parameterization and linear time of day effects.
load("./Outputs/Occupancy_analysis_Swazi_birds_Bayesian_quadTime_WAIC_20180528.RData")

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

## The objects Output1... Output48 are Large jagsUI; Only the posterior MCMC chains are needed.
Outputs <- ((map(Outputs, MCMCchains)))     ## MCMCchains extracts the chains, map applies 
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
rm(list=ls()[! ls() %in% c("Outputs_long","TraitData", "meanAndCRI", "meanAnd90CRI")])
##----------------------------------------------------------------------------------

##----------------------------------------------------------------------------------
## summarize results for the community of common species

results <- matrix(NA, nrow = 3, ncol = 8)
for (x in c(1:8)) {
  results[,x] <- meanAnd90CRI(Outputs_long[,x])
}
colnames(results) <- colnames(Outputs_long)[1:8]
rownames(results) <- c("Mean", "95% LCRI", "95% UCRI")

results90 <- matrix(NA, nrow = 3, ncol = 8)
for (x in c(1:8)) {
  results90[,x] <- meanAnd90CRI(Outputs_long[,x])
}

colnames(results90) <- colnames(Outputs_long)[1:8]
rownames(results90) <- c("Mean", "90% LCRI", "90% UCRI")
res_pooled <- t(rbind(results, results90))
res_pooled <- res_pooled[-3,2:6]
rownames(res_pooled) <- c("Protected", "Shrub cover", "Pasture","Homestead",
                          "Protected x shrub","Pasture x shrub", "Homestead x shrub")

#write.csv(res_pooled, "./outputs/Appendix_4_2.csv")

## No overall effect of any covariate
par(mfrow = c(3,3))
for(i in c(1,2,4:8)) {
  plot(density(Outputs_long[,i]), xlab = "Effect size (beta coeff.)", main = "",
       xlim = c(-20, 20)) 
  abline(v = quantile(Outputs_long[,i], probs = 0.025), col = "blue" )
  abline(v = quantile(Outputs_long[,i], probs = 0.975), col= "blue" )
}

## plot interaction posterior densities together
par(mfrow = c(1, 1))
plot(density(Outputs_long[,6]), xlab = "Land use x shrub cover effect size (beta coeff.)", main = "",
     xlim = c(-30, 30), col = "blue", lwd = 3) 
lines(density(Outputs_long[,7]), col = "salmon", lty = 2, lwd = 3) 
lines(density(Outputs_long[,8]), col = "black", lty = 3, lwd = 3 ) 

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
pro_by_shrub.post <- subset(Outputs_long, select = c("X", "Species","beta.l.int[1]"))  
pas_by_shrub.post <- subset(Outputs_long, select = c("X", "Species","beta.l.int[2]"))  
home_by_shrub.post <- subset(Outputs_long, select = c("X", "Species","beta.l.int[3]"))  

#reformat to make species columns and posterior samples rows with dcast
shrub.post <- dcast(shrub.post, X~Species)
homestead.post <- dcast(homestead.post, X~Species)
pasture.post <- dcast(pasture.post, X~Species)
protected.post <- dcast(protected.post, X~Species)
pro_by_shrub.post <- dcast(pro_by_shrub.post, X~Species)  
pas_by_shrub.post <- dcast(pas_by_shrub.post, X~Species)  
home_by_shrub.post <- dcast(home_by_shrub.post, X~Species)


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

protected.post$nest.shrub <- rowMeans(protected.post[,names(protected.post) %in% nest.shrub$Species])
protected.post$nest.tree <- rowMeans(protected.post[,names(protected.post) %in% nest.tree$Species])
protected.post$nest.cavity <- rowMeans(protected.post[,names(protected.post) %in% nest.cavity$Species])
protected.post$nest.grass <- rowMeans(protected.post[,names(protected.post) %in% nest.grass$Species])

protected.post$predators <- rowMeans(protected.post[,names(protected.post) %in% predators$Species])
protected.post$diet.invert <- rowMeans(protected.post[,names(protected.post) %in% diet.invert$Species])
protected.post$diet.fruit <- rowMeans(protected.post[,names(protected.post) %in% diet.fruit$Species])
protected.post$diet.seeds <- rowMeans(protected.post[,names(protected.post) %in% diet.seed$Species])
protected.post$diet.nectar <- rowMeans(protected.post[,names(protected.post) %in% diet.nectar$Species])

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
pro_by_shrub.post$nest.shrub <- rowMeans(pro_by_shrub.post[,names(pro_by_shrub.post) %in% nest.shrub$Species])
pro_by_shrub.post$nest.tree <- rowMeans(pro_by_shrub.post[,names(pro_by_shrub.post) %in% nest.tree$Species])
pro_by_shrub.post$nest.cavity <- rowMeans(pro_by_shrub.post[,names(pro_by_shrub.post) %in% nest.cavity$Species])
pro_by_shrub.post$nest.grass <- rowMeans(pro_by_shrub.post[,names(pro_by_shrub.post) %in% nest.grass$Species])

pro_by_shrub.post$predators <- rowMeans(pro_by_shrub.post[,names(pro_by_shrub.post) %in% predators$Species])
pro_by_shrub.post$diet.invert <- rowMeans(pro_by_shrub.post[,names(pro_by_shrub.post) %in% diet.invert$Species])
pro_by_shrub.post$diet.fruit <- rowMeans(pro_by_shrub.post[,names(pro_by_shrub.post) %in% diet.fruit$Species])
pro_by_shrub.post$diet.seeds <- rowMeans(pro_by_shrub.post[,names(pro_by_shrub.post) %in% diet.seed$Species])
pro_by_shrub.post$diet.nectar <- rowMeans(pro_by_shrub.post[,names(pro_by_shrub.post) %in% diet.nectar$Species])

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
shrubSppEfx <- apply(shrub.post, 2, meanAnd90CRI )
shrubSppEfx <- as.data.frame(t(shrubSppEfx[, 2:49]))
shrubSppEfx <- cbind(Species, shrubSppEfx)
colnames(shrubSppEfx) <- c("Species", "effect", "LCL", "UCL")
shrubSppEfx$betaName <- "Shrub cover"
rownames(shrubSppEfx)<- NULL

## repeat for the remaining effects
homesteadSppEfx <- apply(homestead.post, 2, meanAnd90CRI )
homesteadSppEfx <- as.data.frame(t(homesteadSppEfx[,2:49]))
homesteadSppEfx <- cbind(Species, homesteadSppEfx)
colnames(homesteadSppEfx) <- c("Species", "effect", "LCL", "UCL")
homesteadSppEfx$betaName <- "Homestead"
rownames(homesteadSppEfx)<- NULL

pastureSppEfx <- apply(pasture.post, 2, meanAnd90CRI )
pastureSppEfx <- as.data.frame(t(pastureSppEfx[,2:49]))
pastureSppEfx <- cbind(Species, pastureSppEfx)
colnames(pastureSppEfx) <- c("Species", "effect", "LCL", "UCL")
pastureSppEfx$betaName <- "Pasture"
rownames(pastureSppEfx)<- NULL

protectedSppEfx <- apply(protected.post, 2, meanAnd90CRI )
protectedSppEfx <- as.data.frame(t(protectedSppEfx[,2:49]))
protectedSppEfx <- cbind(Species, protectedSppEfx)
colnames(protectedSppEfx) <- c("Species", "effect", "LCL", "UCL")
protectedSppEfx$betaName <- "Protected"
rownames(protectedSppEfx)<- NULL

pro_by_shrubSppEfx <- apply(pro_by_shrub.post, 2, meanAnd90CRI )
pro_by_shrubSppEfx <- as.data.frame(t(pro_by_shrubSppEfx[,2:49]))
pro_by_shrubSppEfx <- cbind(Species, pro_by_shrubSppEfx)
colnames(pro_by_shrubSppEfx) <- c("Species", "effect", "LCL", "UCL")
pro_by_shrubSppEfx$betaName <- "pro_by_shrub"
rownames(pro_by_shrubSppEfx)<- NULL

pas_by_shrubSppEfx <- apply(pas_by_shrub.post, 2, meanAnd90CRI )
pas_by_shrubSppEfx <- as.data.frame(t(pas_by_shrubSppEfx[,2:49]))
pas_by_shrubSppEfx <- cbind(Species, pas_by_shrubSppEfx)
colnames(pas_by_shrubSppEfx) <- c("Species", "effect", "LCL", "UCL")
pas_by_shrubSppEfx$betaName <- "pas_by_shrub"
rownames(pas_by_shrubSppEfx)<- NULL

home_by_shrubSppEfx <- apply(home_by_shrub.post, 2, meanAnd90CRI )
home_by_shrubSppEfx <- as.data.frame(t(home_by_shrubSppEfx[,2:49]))
home_by_shrubSppEfx <- cbind(Species, home_by_shrubSppEfx)
colnames(home_by_shrubSppEfx) <- c("Species", "effect", "LCL", "UCL")
home_by_shrubSppEfx$betaName <- "home_by_shrub"
rownames(home_by_shrubSppEfx)<- NULL


## bind rows to create one frame for plotting 
sppEfx <- rbind.data.frame(shrubSppEfx, protectedSppEfx, pastureSppEfx, homesteadSppEfx,
                           pro_by_shrubSppEfx, pas_by_shrubSppEfx, home_by_shrubSppEfx)
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
pro_by_shrubEfxAndTraits <- merge(sppEfx, TraitData) %>% 
  filter(betaName == "pro_by_shrub")
pas_by_shrubEfxAndTraits <- merge(sppEfx, TraitData) %>% 
  filter(betaName == "pas_by_shrub")
home_by_shrubEfxAndTraits <- merge(sppEfx, TraitData) %>% 
  filter(betaName == "home_by_shrub")

##-----------------------------------------------------------------------------------
## build data frames with nest and diet summary effects for plotting
## Nest substrates
r1 <- meanAnd90CRI(shrub.post$nest.shrub)
r2 <- meanAnd90CRI(shrub.post$nest.tree)
r3 <- meanAnd90CRI(shrub.post$nest.cavity)
r4 <- meanAnd90CRI(shrub.post$nest.grass)
shrub.post.nest <- rbind.data.frame(r1,r2,r3,r4)
colnames(shrub.post.nest) <- c("effect", "LCL", "UCL")
shrub.post.nest$betaName <- "Shrub cover"
shrub.post.nest$nest <- c("shrubs", "trees", "cavities", "grass")

r5 <- meanAnd90CRI (protected.post$nest.shrub)
r6 <- meanAnd90CRI (protected.post$nest.tree)
r7 <- meanAnd90CRI (protected.post$nest.cavity)
r8 <- meanAnd90CRI (protected.post$nest.grass)
protected.post.nest <- rbind.data.frame(r5,r6,r7,r8)
colnames(protected.post.nest) <- c("effect", "LCL","UCL")
protected.post.nest$betaName <- "Protected"
protected.post.nest$nest <- c("shrubs", "trees", "cavities","grass")

r9 <- meanAnd90CRI (pasture.post$nest.shrub)
r10 <- meanAnd90CRI (pasture.post$nest.tree)
r11 <- meanAnd90CRI (pasture.post$nest.cavity)
r12<- meanAnd90CRI (pasture.post$nest.grass)
pasture.post.nest <- rbind.data.frame(r9,r10,r11,r12)
colnames(pasture.post.nest) <- c("effect", "LCL","UCL")
pasture.post.nest$betaName <- "Pasture"
pasture.post.nest$nest <- c("shrubs", "trees", "cavities","grass")

r13 <- meanAnd90CRI (homestead.post$nest.shrub)
r14 <- meanAnd90CRI (homestead.post$nest.tree)
r15 <- meanAnd90CRI (homestead.post$nest.cavity)
r16 <- meanAnd90CRI (homestead.post$nest.grass)
homestead.post.nest <- rbind.data.frame(r13,r14,r15,r16)
colnames(homestead.post.nest) <- c("effect", "LCL","UCL")
homestead.post.nest$betaName <- "Homestead"
homestead.post.nest$nest <- c("shrubs", "trees", "cavities","grass")

r17 <- meanAnd90CRI (pro_by_shrub.post$nest.shrub)
r18 <- meanAnd90CRI (pro_by_shrub.post$nest.tree)
r19 <- meanAnd90CRI (pro_by_shrub.post$nest.cavity)
r20 <- meanAnd90CRI (pro_by_shrub.post$nest.grass)
pro_by_shrub.post.nest <- rbind.data.frame(r17,r18,r19,r20)
colnames(pro_by_shrub.post.nest) <- c("effect", "LCL","UCL")
pro_by_shrub.post.nest$betaName <- "pro_by_shrub"
pro_by_shrub.post.nest$nest <- c("shrubs", "trees", "cavities","grass")

r21 <- meanAnd90CRI (pas_by_shrub.post$nest.shrub)
r22 <- meanAnd90CRI (pas_by_shrub.post$nest.tree)
r23 <- meanAnd90CRI (pas_by_shrub.post$nest.cavity)
r24 <- meanAnd90CRI (pas_by_shrub.post$nest.grass)
pas_by_shrub.post.nest <- rbind.data.frame(r21,r22,r23,r24)
colnames(pas_by_shrub.post.nest) <- c("effect", "LCL","UCL")
pas_by_shrub.post.nest$betaName <- "pas_by_shrub"
pas_by_shrub.post.nest$nest <- c("shrubs", "trees", "cavities","grass")

r25 <- meanAnd90CRI (home_by_shrub.post$nest.shrub)
r26 <- meanAnd90CRI (home_by_shrub.post$nest.tree)
r27 <- meanAnd90CRI (home_by_shrub.post$nest.cavity)
r28 <- meanAnd90CRI (home_by_shrub.post$nest.grass)
home_by_shrub.post.nest <- rbind.data.frame(r25,r26,r27,r28)
colnames(home_by_shrub.post.nest) <- c("effect", "LCL","UCL")
home_by_shrub.post.nest$betaName <- "home_by_shrub"
home_by_shrub.post.nest$nest <- c("shrubs", "trees", "cavities","grass")

## Diets
r29 <- meanAnd90CRI(shrub.post$predators)
r30 <- meanAnd90CRI(shrub.post$diet.invert)
r31 <- meanAnd90CRI(shrub.post$diet.fruit)
r32 <- meanAnd90CRI(shrub.post$diet.seeds)
r33 <- meanAnd90CRI(shrub.post$diet.nectar)
shrub.post.diet <- rbind.data.frame(r29,r30,r31,r32,r33)
colnames(shrub.post.diet) <- c("effect", "LCL","UCL")
shrub.post.diet$betaName <- "Shrub cover"
shrub.post.diet$diet <- c("predatory","invertebrates","fruit","seeds","nectar")

r34 <- meanAnd90CRI (protected.post$predators)
r35 <- meanAnd90CRI (protected.post$diet.invert)
r36 <- meanAnd90CRI (protected.post$diet.fruit)
r37 <- meanAnd90CRI (protected.post$diet.seeds)
r38 <- meanAnd90CRI (protected.post$diet.nectar)
protected.post.diet <- rbind.data.frame(r34,r35,r36,r37,r38)
colnames(protected.post.diet) <- c("effect", "LCL","UCL")
protected.post.diet$betaName <- "Protected"
protected.post.diet$diet <- c("predatory","invertebrates","fruit","seeds","nectar")

r39 <- meanAnd90CRI (pasture.post$predators)
r40 <- meanAnd90CRI (pasture.post$diet.invert)
r41 <- meanAnd90CRI (pasture.post$diet.fruit)
r42 <- meanAnd90CRI (pasture.post$diet.seeds)
r43 <- meanAnd90CRI (pasture.post$diet.nectar)
pasture.post.diet <- rbind.data.frame(r39,r40,r41,r42,r43)
colnames(pasture.post.diet) <- c("effect", "LCL","UCL")
pasture.post.diet$betaName <- "Pasture"
pasture.post.diet$diet <- c("predatory","invertebrates","fruit","seeds","nectar")

r44 <- meanAnd90CRI (homestead.post$predators)
r45 <- meanAnd90CRI (homestead.post$diet.invert)
r46 <- meanAnd90CRI (homestead.post$diet.fruit)
r47 <- meanAnd90CRI (homestead.post$diet.seeds)
r48 <- meanAnd90CRI (homestead.post$diet.nectar)
homestead.post.diet <- rbind.data.frame(r44,r45,r46,r47,r48)
colnames(homestead.post.diet) <- c("effect", "LCL","UCL")
homestead.post.diet$betaName <- "Homestead"
homestead.post.diet$diet <- c("predatory","invertebrates","fruit","seeds","nectar")

r49 <- meanAnd90CRI (pro_by_shrub.post$predators)
r50 <- meanAnd90CRI (pro_by_shrub.post$diet.invert)
r51 <- meanAnd90CRI (pro_by_shrub.post$diet.fruit)
r52 <- meanAnd90CRI (pro_by_shrub.post$diet.seeds)
r53 <- meanAnd90CRI (pro_by_shrub.post$diet.nectar)
pro_by_shrub.post.diet <- rbind.data.frame(r49,r50,r51,r52,r53)
colnames(pro_by_shrub.post.diet) <- c("effect", "LCL","UCL")
pro_by_shrub.post.diet$betaName <- "pro_by_shrub"
pro_by_shrub.post.diet$diet <- c("predatory","invertebrates","fruit","seeds","nectar")

r54 <- meanAnd90CRI (pas_by_shrub.post$predators)
r55 <- meanAnd90CRI (pas_by_shrub.post$diet.invert)
r56 <- meanAnd90CRI (pas_by_shrub.post$diet.fruit)
r57 <- meanAnd90CRI (pas_by_shrub.post$diet.seeds)
r58 <- meanAnd90CRI (pas_by_shrub.post$diet.nectar)
pas_by_shrub.post.diet <- rbind.data.frame(r54,r55,r56,r57,r58)
colnames(pas_by_shrub.post.diet) <- c("effect", "LCL","UCL")
pas_by_shrub.post.diet$betaName <- "pas_by_shrub"
pas_by_shrub.post.diet$diet <- c("predatory","invertebrates","fruit","seeds","nectar")

r59 <- meanAnd90CRI (home_by_shrub.post$predators)
r60 <- meanAnd90CRI (home_by_shrub.post$diet.invert)
r61 <- meanAnd90CRI (home_by_shrub.post$diet.fruit)
r62 <- meanAnd90CRI (home_by_shrub.post$diet.seeds)
r63 <- meanAnd90CRI (home_by_shrub.post$diet.nectar)
home_by_shrub.post.diet <- rbind.data.frame(r59,r60,r61,r62,r63)
colnames(home_by_shrub.post.diet) <- c("effect", "LCL","UCL")
home_by_shrub.post.diet$betaName <- "home_by_shrub"
home_by_shrub.post.diet$diet <- c("predatory","invertebrates","fruit","seeds","nectar")

## collect
nest.summary <- rbind.data.frame(shrub.post.nest, protected.post.nest, pasture.post.nest, homestead.post.nest,
                                 pro_by_shrub.post.nest, pas_by_shrub.post.nest, home_by_shrub.post.nest,
                                 stringsAsFactors = FALSE)
diet.summary <- rbind.data.frame(shrub.post.diet, protected.post.diet, pasture.post.diet ,homestead.post.diet,
                                 pro_by_shrub.post.diet, pas_by_shrub.post.diet, home_by_shrub.post.diet,
                                 stringsAsFactors = FALSE)

###------------------------------------------------------------------------------------
#some graph code from Isabel's paper
#Fig 4: species-specific forest plots with multiple effect sizes
forest_plot_spp <- ggplot(data = sppEfx, aes(x = Species, y = effect)) + 
  geom_errorbar(aes(ymin = LCL, ymax = UCL), width = 0) +
  geom_point(colour="gray20", shape = 21, size = 4, 
             aes(fill = factor(POSNEG))) +                     
  scale_fill_manual(values = c("gray50","orange","blue")) +
  geom_hline(aes(yintercept = 0)) +
  ylab("Median" ~beta~ ", 95% CRI") +         
  scale_x_discrete(limits = rev(sppEfx[order(sppEfx$RobertsNum),]$Species), expand = c(0.02,0.01)) +
  coord_flip() +  
  facet_wrap(~factor(betaName,levels = c("Protected", "Pasture", "Homestead","Plantation","Shrub cover",
                                         "pro_by_shrub","pas_by_shrub","home_by_shrub")), 
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
pd <- position_dodge(0.5) # move them .05 to the left and right
#levels = c("Protected", "Pasture", "Homestead","Plantation","Shrub cover"))

#nest substrates
nest.fig <- ggplot(data = nest.summary, aes(x = betaName, y = effect, group = nest)) +
  geom_errorbar(aes(ymin = LCL, ymax = UCL), width = 0, size = 0.3, position = pd) +
  geom_point(colour = "gray20", shape = 21, size = 3, position = pd, aes(fill = factor(nest))) + 
  scale_fill_manual(values = c("white", "yellow", "gray50", "blue"), name = "Nesting substrate")+
  geom_hline(aes(yintercept = 0)) +
  ylab("Median" ~beta~ ", 95% CRI") +  
  ylim(c(-15, 15)) +
  scale_x_discrete(limits = rev(levels(factor(nest.summary$betaName, 
                  levels = c("Protected", "Pasture", "Homestead","Shrub cover",
                             "pro_by_shrub","pas_by_shrub","home_by_shrub") )) )) +
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

diet.fig <- ggplot(data = diet.summary, aes(x = betaName, y = effect, group = diet)) + 
  geom_errorbar(aes(ymin = LCL, ymax = UCL), width = 0, size = 0.3, position = pd) +
  geom_point(colour = "gray20", shape = 21, size = 3, position = pd, aes(fill = factor(diet))) + 
  scale_fill_manual(values = c("white", "yellow", "gray50","blue","black"), name = "Diet") +
  geom_hline(aes(yintercept = 0)) +
  ylab("Median" ~beta~ ", 95% CRI") +  
  ylim(c(-15, 15)) +
  scale_x_discrete(limits = rev(levels(factor(nest.summary$betaName, 
                   levels = c("Protected", "Pasture", "Homestead", "Shrub cover",
                              "pro_by_shrub","pas_by_shrub","home_by_shrub") )) )) +
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
prot_mass <- merge(sppEfx, TraitData) %>% 
  filter(betaName == "Protected")
past_mass <- merge(sppEfx, TraitData) %>% 
  filter(betaName == "Pasture")
home_mass <- merge(sppEfx, TraitData) %>% 
  filter(betaName == "Homestead")
pro_by_shrub_mass <- merge(sppEfx, TraitData) %>% 
  filter(betaName == "pro_by_shrub")
pas_by_shrub_mass <- merge(sppEfx, TraitData) %>% 
  filter(betaName == "pas_by_shrub")
home_by_shrub_mass <- merge(sppEfx, TraitData) %>% 
  filter(betaName == "home_by_shrub")


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
        axis.text.y  = element_text(size=12, colour = "black", margin = unit(c(0.3, 0.3, 0.3, 0.3),
                                                                             "cm")),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(), 
        legend.position = "none")


mass.fig.prot <- ggplot(data = prot_mass, aes(x = log(Mass), y = effect)) + 
  geom_errorbar(aes(ymin = LCL, ymax = UCL), width = 0) +
  geom_point(colour="gray20", shape=21, size = 4, 
             aes(fill = factor(POSNEG))) +                     
  scale_fill_manual(values = c("gray50","orange","blue"))+
  ylim(-30, 30)+
  theme_bw()+
  ggtitle("Protected areas") +
  theme(
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.title.y = element_blank(), 
    axis.text.y  = element_text(size = 12, colour = "black", 
                                margin = unit(c(0.3,0.3,0.3,0.3), "cm")),
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
    axis.text.y  = element_blank(),
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(), 
    legend.position = "none")


mass.fig.home <- ggplot(data = home_mass, aes(x = log(Mass), y = effect)) + 
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

mass.fig.pro_by_shrub <- ggplot(data = pro_by_shrub_mass, aes(x = log(Mass), y = effect)) + 
  geom_errorbar(aes(ymin = LCL, ymax = UCL), width = 0) +
  geom_point(colour="gray20", shape=21, size = 4, 
             aes(fill = factor(POSNEG))) +                     
  scale_fill_manual(values = c("gray50","orange","blue"))+
  ylim(-30, 30)+
  theme_bw()+
  ggtitle("pro_by_shrub") +
  theme(element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(), 
        axis.text.y  = element_text(size=12, colour = "black", 
                                    margin = unit(c(0.3, 0.3, 0.3, 0.3), "cm")),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(), 
        legend.position = "none")

mass.fig.pas_by_shrub <- ggplot(data = pas_by_shrub_mass, aes(x = log(Mass), y = effect)) + 
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

mass.fig.home_by_shrub <- ggplot(data = home_by_shrub_mass, aes(x = log(Mass), y = effect)) + 
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
        axis.text.y  = element_text(size=12, colour = "black", 
                                    margin = unit(c(0.3, 0.3, 0.3, 0.3), "cm")),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(), 
        legend.position = "none")


ps.load.fig.prot <- ggplot(data = prot_mass, aes(x = pseudo_loading, y = effect)) + 
  geom_errorbar(aes(ymin = LCL, ymax = UCL), width = 0) +
  geom_point(colour="gray20", shape=21, size = 4, 
             aes(fill = factor(POSNEG))) +                     
  scale_fill_manual(values = c("gray50","orange","blue"))+
  ylim(-30, 30)+
  theme_bw()+
  ggtitle("Protected areas") +
  theme(
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.title.y = element_blank(), 
    axis.text.y  = element_text(size = 12, colour = "black", 
                                margin = unit(c(0.3,0.3,0.3,0.3), "cm")),
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
    axis.text.x  = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
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
        axis.text.x = element_blank(),
        axis.title.y = element_blank(), 
        axis.text.y  = element_text(size=12, colour = "black", 
                                    margin = unit(c(0.3, 0.3, 0.3, 0.3), "cm")),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(), 
        legend.position = "none")

ps.load.fig.pro_by_shrub <- ggplot(data = pro_by_shrub_mass, aes(x = pseudo_loading, y = effect)) + 
  geom_errorbar(aes(ymin = LCL, ymax = UCL), width = 0) +
  geom_point(colour="gray20", shape=21, size = 4, 
             aes(fill = factor(POSNEG))) +                     
  scale_fill_manual(values = c("gray50","orange","blue"))+
  ylim(-30, 30)+
  theme_bw()+
  ggtitle("pro_by_shrubs") +
  theme(element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(), 
        axis.text.y  = element_text(size=12, colour = "black", 
                                    margin = unit(c(0.3, 0.3, 0.3, 0.3), "cm")),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(), 
        legend.position = "none")

ps.load.fig.pas_by_shrub <- ggplot(data = pas_by_shrub_mass, aes(x = pseudo_loading, y = effect)) + 
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

ps.load.fig.home_by_shrub <- ggplot(data = home_by_shrub_mass, aes(x = pseudo_loading, y = effect)) + 
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

# Violin and box plots ----------------------------------------------
## By diet, in PAs
xxd1 <- ggplot(data = Outputs_long, aes(factor(Diet), beta)) + 
  geom_violin(aes(fill = factor(Diet)), show.legend = FALSE) + 
  geom_boxplot(fill = "white", width = 0.1, outlier.size = 0.2) +
  geom_hline(aes(yintercept = 0)) +
  ggtitle("Protected areas") +
  scale_y_continuous(limits = c(-40, 40)) + 
  ylab("") +  
  coord_flip() +
  theme_bw() +
  theme(axis.text.x  = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y  = element_text(size = 9, colour = "black", 
                                    margin=unit(c(0.2, 0.2, 0.2, 0.2), "cm")),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.ticks.length=unit(-0.1, "cm") )

## By diet, in pastures
xxd2 <- ggplot(data = Outputs_long, aes(factor(Diet), `beta.l[2]`)) + 
  geom_violin(aes(fill = factor(Diet)), show.legend = FALSE) + 
  geom_boxplot(fill = "white", width = 0.1, outlier.size = 1) +
  geom_hline(aes(yintercept = 0)) +
  ggtitle("Pastures") +
  scale_y_continuous(limits = c(-40, 40)) + 
  coord_flip() +
  ylab("") +  
  theme_bw() +
  theme(axis.text.x  = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y  = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.ticks.length=unit(-0.1, "cm") )

## By diet, in homesteads
xxd3 <- ggplot(data = Outputs_long, aes(factor(Diet), `beta.l[3]`)) + 
  geom_violin(aes(fill = factor(Diet)), show.legend = FALSE) + 
  geom_boxplot(fill = "white", width = 0.1, outlier.size = 1) +
  geom_hline(aes(yintercept = 0)) +
  ggtitle("Homesteads") +
  scale_y_continuous(limits = c(-40, 40)) + 
  ylab("") +  
  coord_flip() +
  theme_bw() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y  = element_text(size = 9, colour = "black", 
                                margin=unit(c(0.2, 0.2, 0.2, 0.2), "cm")),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    axis.ticks.length=unit(-0.1, "cm") )

## by diet and z-transformed shrub cover
xxd5 <- ggplot(data = Outputs_long, aes(factor(Diet), beta1)) + 
  geom_violin(aes(fill = factor(Diet)), show.legend = FALSE) + 
  geom_boxplot(fill = "white", width = 0.1, outlier.size = 1) +
  geom_hline(aes(yintercept = 0)) +
  ggtitle("Shrub cover") +
  scale_y_continuous(limits = c(-20, 30)) + 
  ylab("Median" ~beta~ ", 95% CRI") +  
  coord_flip() +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y  = element_text(size = 9, colour = "black", 
                                    margin=unit(c(0.2, 0.2, 0.2, 0.2), "cm")),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.ticks.length=unit(-0.1, "cm") )


xxd5 <- ggplot(data = Outputs_long, aes(factor(Diet), `beta.l.int[1]`)) + 
  geom_violin(aes(fill = factor(Diet)), show.legend = FALSE) + 
  geom_boxplot(fill = "white", width = 0.1, outlier.size = 1) +
  geom_hline(aes(yintercept = 0)) +
  ggtitle("pro_byshrub") +
  scale_y_continuous(limits = c(-40, 40)) + 
  ylab("") +  
  coord_flip() +
  theme_bw() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y  = element_text(size = 9, colour = "black", 
                                margin=unit(c(0.2, 0.2, 0.2, 0.2), "cm")),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    axis.ticks.length=unit(-0.1, "cm") )


xxd6 <- ggplot(data = Outputs_long, aes(factor(Diet), `beta.l.int[2]`)) + 
  geom_violin(aes(fill = factor(Diet)), show.legend = FALSE) + 
  geom_boxplot(fill = "white", width = 0.1, outlier.size = 1) +
  geom_hline(aes(yintercept = 0)) +
  ggtitle("Homesteads") +
  scale_y_continuous(limits = c(-40, 40)) + 
  ylab("") +  
  coord_flip() +
  theme_bw() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y  = element_text(size = 9, colour = "black", 
                                margin=unit(c(0.2, 0.2, 0.2, 0.2), "cm")),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    axis.ticks.length=unit(-0.1, "cm") )


xxd7 <- ggplot(data = Outputs_long, aes(factor(Diet), `beta.l.int[3]`)) + 
  geom_violin(aes(fill = factor(Diet)), show.legend = FALSE) + 
  geom_boxplot(fill = "white", width = 0.1, outlier.size = 1) +
  geom_hline(aes(yintercept = 0)) +
  ggtitle("Homesteads") +
  scale_y_continuous(limits = c(-40, 40)) + 
  ylab("") +  
  coord_flip() +
  theme_bw() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y  = element_text(size = 9, colour = "black", 
                                margin=unit(c(0.2, 0.2, 0.2, 0.2), "cm")),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    axis.ticks.length=unit(-0.1, "cm") )

colnames(diet.summary)[5] <- "trait"
colnames(nest.summary)[5] <- "trait"
colnames(diet.summary)[2:3] <- c("90% LCL", "90% UCL")
colnames(nest.summary)[2:3] <-  c("90% LCL", "90% UCL")
cat_traits_90 <- rbind.data.frame(nest.summary, diet.summary)
write.csv(cat_traits_90, "./outputs/AppendixS4/cat_traits_90.2.csv")
