# traits_explaining_occupancy
data and code supporting the manuscript, Traits explaining shrub encroachment and land use effects on bird occupancy in an African savanna

Files:
## updated_study_area.R
     Produces: Fig. 1 Study area maps and shrub cover strip plots. These need to be manually added to study area photos and arranged in
               PowerPoint, then exported as .tiff at suitable resolution.
              ./GenerateAbunHists.RData  ## Array of detection-nondetection events for each survey station
              ./data/siteCovs_numeric_grids.csv

## n_fold_single_spp_occ_Swazi.R, n_fold_single_spp_occ_int_models.R, n_fold_single_spp_occ_3_land_uses_add_efx.R,
     Requires:./GenerateAbunHists.RData          ## Array of detection-nondetection events for each survey station
              ./data/SamplingCovs.csv            ## Sampling covariate data
              ./data/siteCovs_numeric_grids.csv  ## Site covariate data
              ./code/jagsmodel1.R                ## jags model with additive effects only
              ./code/jagsmodel2.R                ## jags model with land use x shrub cover interaction terms added, plantation surveys excluded
              ./code/jagsmodel_3_land_uses.R     ## Additive effects only, plantation surveys excluded
	Produces: 
              ./Occupancy_analysis_Swazi_birds_Bayesian_quadTime_fixed.RData
              ./Occupancy_analysis_Swazi_birds_Bayesian_int_models_quadTime_fixed.RData
              ./Occupancy_analysis_Swazi_birds_Bayesian_quadTime_simplified_20180505.RData
        *These files contain jagsUI objects used in downstream scripts, mainly the posterior samples. 
        *Note that running the script does not implement different tuning parameters for the handful of species 
         needing additional iterations. That needs to be implemented manually by modifying the "in" argument in the for loop.

## n_fold_single_spp_occ_get_WAIC.R, n_fold_single_spp_occ_get_WAIC_2.R
              Run models and collect posterior pd and lpd, then compute WAIC for the additive and interaction models. This takes 
              >1gB of memory, so we ran it independently of saving posterior parameter estimates from the same models.

## compare_WAIC.R
   Requires:./outputs/WAIC_additive_models
            ./outputs/WAIC_interaction_models_2"
            ./outputs/WAIC_interaction_models_updates
   Produces: When called by other scripts, this introduces into the global environment a 
             1)vector of species for which WAIC indicated support for the model with a land use by shrub cover interaction, and;
             2) summery stats on how many species had WIAC support at different thresholds for deltaWAIC.

## PosteriorsBySPP_20180717.R    ## updated figures from 20180503 script
     Requires: ./Occupancy_analysis_Swazi_birds_Bayesian_quadTime_fixed.RData
               ./data/TraitData.csv
     Produces:
               1) Figs. 2-3
               2) 1st part of Appendix 4
               3) Appendix 7
               4) Appendix 8

## posteriorsBySpp_interactions20180601
     Requires: ./Occupancy_analysis_Swazi_birds_Bayesian_int_models_quadTime_fixed.RData
               ./data/TraitData.csv
     Produces:
      1) Model DIC and deviance ## not used
      2) Appendix S3 part 2 ## mean, 90, and 95% CRIs for main efffects and interactions 
      3) Posterior density plot of pooled interaction term effects and main effects  ## not used
      4) Appendix S4 mean, 90, and 95% CRIs of continuous trait effects on each land-cover and land-use response (slopes) 
      5) AppendixS4nest_trait_effects.2.csv  ## nest substrate effects under this model
      6) AppendixS4diet_trait_effects.2.csv  ## diet and predatory status effects under this model
      7) View Figs. 2-4 produced using the main effects from the interaction effects terms models
      8) View results of Fisher tests and t-tests comparing traits of spp. for which interaction is supported or not

## posteriorsBySpp_interactions_90CRI.R
     Produces the same materials as above but with 90% CRIs, with some extraneous code producing violin plots that was not used. Figures with 90% CRIs
     can be viewed but were not saved as graphics files.
     1) ./outputs/Appendix_4_2.csv
     2) ./outputs/AppendixS4/cat_traits_90.2.csv
  
## plot_pred_psi.R
      Requires: ./Occupancy_analysis_Swazi_birds_Bayesian_int_models_quadTime_fixed.RData
                ./Occupancy_analysis_Swazi_birds_Bayesian_quadTime_simplified_20180505.RData
                ./Occupancy_analysis_Swazi_birds_Bayesian_quadTime_fixed.RData
                ./code/code/compare_WAIC.R  ## sourced in script
      Produces: 
               1) Fig 4  ## Community-level effects of shrub cover for each land use
               2) Fig S1    ## One panel for each species, showing predicting occupancy for each land use across observed range of shrub cover,
                            ## all 4 land uses, no interactions
               3) Fig S_spp ## Figure with panels only for species with interaction supported by WAIC
               4) Fig SXX   ## Plot all species responses to shrub cover together, with a panel for each land use ## Not used

## write_model_results_table.R
    ## Assemble MCMC results from the single-species occupancy models
    ##   1) view which models have R-hat >= 1.105 and up iterations if this occurs
    ##   2) Write a .csv file containing tuning parameters used for all models, by species
    ##   3) View saved beta parameters for all species, and write to two csv files for manual combination in Excel
    ##   4) view traceplots for all models
    ##   5) view forest plots of interaction-terms models for all species

    Requires: ./Occupancy_analysis_Swazi_birds_Bayesian_quadTime_simplified_20180505.RData
              ./Occupancy_analysis_Swazi_birds_Bayesian_quadTime_fixed.RData
              ./outputs/Occupancy_analysis_Swazi_birds_Bayesian_quadTime_WAIC_20180528.RData
    Produces:
           ./outputs/AppendixS2.1.csv
           ./outputs/AppendixS2.2.csv
           ./outputs/AppendixS2.3.csv




