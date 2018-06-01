## Compare WAIC for two models 

## Data. WAIC values for the two models
occ.WAIC_add <- readRDS("./outputs/WAIC_additive_models")
occ.WAIC_int <- readRDS("./outputs/WAIC_interaction_models_2")
occ.WAIC_updates <- readRDS("./outputs/WAIC_interaction_models_updates")
## replace WAIC_int with updated values for species that required extra iterations
for (i in c(10, 18, 27, 29, 38, 42)) {
  occ.WAIC_int[i] <- occ.WAIC_updates[i]
}

## Differences in WAIC
delta_WAIC <- abs(occ.WAIC_add - occ.WAIC_int)
summary(delta_WAIC)

## Number of species where the model with an interaction is supported by WAIC
## 17 of 48, 35.4% of species
sum(occ.WAIC_int < occ.WAIC_add)
sum(occ.WAIC_int < occ.WAIC_add)/48

## For which species is an interaction supported? 
int_supported <- which(occ.WAIC_int < occ.WAIC_add )
int_supported

##
summary(occ.WAIC_int[int_supported] - occ.WAIC_add[int_supported])

## Which and how many are "strongly" supported?
abs(int_supported) >= 7
sum(abs(int_supported) >= 7)    ## 16
sum(abs(int_supported) >= 7)/48 ## 33.3% of species
