## Assemble MCMC results from the single-species occupancy models fitted in FILENAME to:
##   1) view which models have R-hat >= 1.105
##   2) Write a .csv file containing tuning parameters used for all models, by species
##   3) View saved beta parameters for all species, and write to two csv files for manual
##      combination in Excel
##   4) view traceplots for all models
##   5) view forest plots of interaction-terms models for all species

# Packages ----------------------------------------------------------------
library(MCMCvis)
library(tidyverse)

# Data --------------------------------------------------------------------
## load models without interaction terms
load(("./Occupancy_analysis_Swazi_birds_Bayesian_quadTime_fixed.RData"))

## load models with interaction terms
load("./outputs/Occupancy_analysis_Swazi_birds_Bayesian_quadTime_WAIC_20180528.RData")

## get Roberts number, Latin name, and common name to sort tables 
TraitData <- read.csv("./data/TraitData.csv")
spp_name_info <- TraitData[,1:3]


# Processing --------------------------------------------------------------
## assemble model outputs into lists
outputs <- list(Output1,Output2,Output3,Output4,Output5,
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

outputs_int <- list(IntOutput1,IntOutput2,IntOutput3,IntOutput4,IntOutput5,
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

## create rows for species common and Latin names for downstream processing
sppName <- vector(mode = "character", length = 48)  
sppName <- dimnames(AbunHists.subset)[[3]] 
spp_name_info <- spp_name_info[which(spp_name_info$Species %in% as.factor(sppName) ),]
spp_name_info <- spp_name_info[order(spp_name_info$Species),] 

## Find any Rhat >= 1.105
which (unlist (lapply(outputs, function(x){max(unlist(x$Rhat), na.rm = TRUE) >= 1.105 }) ) )
which (unlist (lapply(outputs_int, function(x){max(unlist(x$Rhat), na.rm = TRUE) >= 1.105 }) ) )

## Get summaries of model parameters and assemble into a table
output_summary <- lapply(outputs, function(x){
MCMCsummary(x, 
            excl = c('deviance', 'beta.l\\[1\\]'),
            ISB = FALSE,
            Rhat = FALSE,
            digits = 3,
            func = function(x) quantile(x, probs = c(0.05, 0.95)),
            func_name = c('5%', '95%'))}
)

output_int_summary <- lapply(outputs_int, function(x){
  MCMCsummary(x, 
              excl = c('deviance', 'beta.l\\[1\\]'),
              ISB = FALSE,
              Rhat = FALSE,
              digits = 3,
              func = function(x) quantile(x, probs = c(0.05, 0.95)),
              func_name = c('5%', '95%'))}
)

## drop superfluous columns
output_summary <- lapply(output_summary, function(x) {x[,3:7] })
output_int_summary <- lapply(output_int_summary, function(x) {x[c(1:4,6,7),3:7] })

names(output_summary) <- sppName 
names(output_int_summary) <- sppName 

## re-arrange into a 48 spp by 40 item dataframe
output_summary <- data.frame(t(sapply(output_summary,c)))
output_int_summary <- data.frame(t(sapply(output_int_summary,c)))

## Change column names
colnames(output_summary) <- c("Date, 2.5%", "50%", "97.5%", "5%", "95%",
                              "Time, 2.5%", "50%", "97.5%", "5%", "95%",
                              "Time, quadratic, 2.5%", "50%", "97.5%", "5%", "95%",
                              "Intercept, 2.5%", "50%", "97.5%", "5%", "95%",
                              "Shrub cover, 2.5%", "50%", "97.5%", "5%", "95%",
                              "Pasture,2.5%", "50%", "97.5%", "5%", "95%",
                              "Homestead, 2.5%", "50%", "97.5%", "5%", "95%",
                              "Plantation, 2.5%", "50%", "97.5%", "5%", "95%") 

colnames(output_int_summary) <- c("Intercept, 2.5%", "50%", "97.5%", "5%", "95%", 
                                  "Shrub cover, 2.5%", "50%", "97.5%", "5%", "95%",
                                  "Pasture, 2.5%", "50%", "97.5%", "5%", "95%",
                                  "Homestead, 2.5%", "50%", "97.5%", "5%", "95%",
                                  "Pasture x shrub, 2.5%", "50%", "97.5%", "5%", "95%",
                                  "Homestead x shrub, 2.5%", "50%", "97.5%", "5%", "95%")

## Get WAIC 
## Data. WAIC values for the two models
occ.WAIC_int <- readRDS("./outputs/WAIC_interaction_models_2")
occ.WAIC_updates <- readRDS("./outputs/WAIC_interaction_models_updates")
## replace WAIC_int with updated values for species that required extra iterations
for (i in c(10, 18, 27, 29, 38, 42)) {
  occ.WAIC_int[i] <- occ.WAIC_updates[i]
}

## Add WAIC to tables
output_summary$WAIC <- rep("NA", 48)
output_int_summary$WAIC <- occ.WAIC_int

write.csv(output_summary, "./outputs/AppendixS2.1.csv")
write.csv(output_int_summary, "./outputs/AppendixS2.2.csv")

##--------------------------
## Repeat for 3 land use, additive effexts model 
## remove everything since names are repeated
rm(list = ls())
load("./Occupancy_analysis_Swazi_birds_Bayesian_quadTime_simplified_20180505.RData")

## get Roberts number, Latin name, and common name to sort tables 
TraitData <- read.csv("./data/TraitData.csv")
spp_name_info <- TraitData[,1:3]


# Processing --------------------------------------------------------------
## assemble model outputs into lists
outputs <- list(Output1,Output2,Output3,Output4,Output5,
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

## create rows for species common and Latin names for downstream processing
sppName <- vector(mode = "character", length = 48)  
sppName <- dimnames(AbunHists.subset)[[3]] 
spp_name_info <- spp_name_info[which(spp_name_info$Species %in% as.factor(sppName) ),]
spp_name_info <- spp_name_info[order(spp_name_info$Species),] 

## Get summaries of model parameters and assemble into a table
output_summary <- lapply(outputs, function(x){
  MCMCsummary(x, 
              excl = c('deviance', 'beta.l\\[1\\]',"alpha2","alpha3","alpha4"),
              ISB = FALSE,
              Rhat = FALSE,
              digits = 3,
              func = function(x) quantile(x, probs = c(0.05, 0.95)),
              func_name = c('5%', '95%'))}
)

## drop superfluous columns and add species names
output_summary <- lapply(output_summary, function(x) {x[,3:7] })
names(output_summary) <- sppName 

## re-arrange into a 48 spp by 40 item dataframe
output_summary <- data.frame(t(sapply(output_summary,c)))

## Change column names
colnames(output_summary) <- c("Intercept, 2.5%", "50%", "97.5%", "5%", "95%",
                              "Shrub cover, 2.5%", "50%", "97.5%", "5%", "95%",
                              "Pasture,2.5%", "50%", "97.5%", "5%", "95%",
                              "Homestead, 2.5%", "50%", "97.5%", "5%", "95%") 
## Get WAIC 
## Data. WAIC values 
occ.WAIC_add <- readRDS("./outputs/WAIC_additive_models")
output_summary$WAIC <- occ.WAIC_add
write.csv(output_summary, "./outputs/AppendixS2.3.csv")

