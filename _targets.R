#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#### SCRIPT INTRODUCTION ####
#
#' @name _targets.R  
#' @description R script to launch the target pipeline
#' @author Julien BARRERE
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Options and packages ----------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Load targets
library(targets)
# Load functions
lapply(grep("R$", list.files("R"), value = TRUE), function(x) source(file.path("R", x)))
# install if needed and load packages
packages.in <- c("dplyr", "ggplot2", "matreex", "tidyr", "data.table", 
                 "factoextra", "modi", "sf", "rnaturalearth", "scales", 
                 "cowplot", "multcomp", "piecewiseSEM", "future", "FD", "GGally", 
                 "statmod")
for(i in 1:length(packages.in)) if(!(packages.in[i] %in% rownames(installed.packages()))) install.packages(packages.in[i])
# Targets options
options(tidyverse.quiet = TRUE, clustermq.scheduler = "multiprocess")
tar_option_set(packages = packages.in,
               memory = "transient")
future::plan(future::multisession, workers = 6)
set.seed(2)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Targets workflow --------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

list(
  
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # -- Load data -----
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  # File with data of climate and species present per FUNDIV plot
  tar_target(FUNDIV_climate_species_file, "data/FUNDIV_climate_species.csv",
             format = "file"), 
  
  # Read the file
  tar_target(FUNDIV_climate_species, fread(FUNDIV_climate_species_file)),
  
  # List of species for which we have data
  tar_target(all.species.name, 
             colnames(FUNDIV_climate_species)[grep("_", colnames(FUNDIV_climate_species))]), 
  
  # Get demographic parameters for all species
  tar_target(fit.list.allspecies, load_param_demo(all.species.name)),
  
  # Generate some harvest rules
  tar_target(harv_rules.ref, c(Pmax = 0.25, dBAmin = 3, freq = 10, alpha = 1)),
  
  # Data.frame containing storm disturbance to apply
  tar_target(disturbance.df_storm, data.frame(
    type = rep("storm", 3), intensity = rep(0.5, 3), 
    IsSurv = rep(FALSE, 3), t = c(500:502))),
  
  
  
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # -- Make simulations with storm disturbances -----
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  # Generate some climates
  # -- iterations along all climates that will be created (one iteration per climate)
  tar_target(ID.climate_storm, c(1:10)),
  # -- list of climates
  tar_target(climate_list_storm, create_climate_list(length(ID.climate_storm), 
                                                     quantile.range = c(0, 0.8))),
  # -- generate one climate object per iteration with branching
  tar_target(climate_storm, make_climate(
      FUNDIV_climate_species, quantiles.in = climate_list_storm[[ID.climate_storm]], 
      "storm", 10, exclude.in = c("Carpinus_betulus", "Quercus_ilex"), 
      method = "frequency"), pattern = map(ID.climate_storm), iteration = "list"), 
  
  # Make species objects
  # -- First: list all species object to make
  tar_target(species_list_storm, make_species_list(climate_storm, "storm")),
  # -- Make a vector of ID for each species to make
  tar_target(ID.species_storm, species_list_storm$ID.species), 
  # -- Make the species via branching over ID.species
  tar_target(species_storm, make_species_rds(
    fit.list.allspecies, climate_storm, species_list_storm, 
    ID.species.in = ID.species_storm), 
    pattern = map(ID.species_storm), iteration = "vector", format = "file"),
  
  # Make simulations 
  # -- Start with a list of forest to simulate
  tar_target(forest_list_storm, make_forest_list(climate_storm, "storm")), 
  # -- Make a vector of ID for each forest to simulate
  tar_target(ID.forest_storm, forest_list_storm$ID.forest),
  # -- Make simulations till equilibrium
  tar_target(sim_equilibrium_storm, make_simulations_equilibrium(
   climate_storm, harv_rules.ref, species_list_storm, forest_list_storm, species_storm, ID.forest_storm), 
   pattern = map(ID.forest_storm), iteration = "vector", format = "file"),
  # -- Make simulations with disturbance
  tar_target(sim_disturbance_storm, make_simulations_disturbance(
   climate_storm, harv_rules.ref, species_list_storm, forest_list_storm, 
   species_storm, sim_equilibrium_storm, ID.forest.in = ID.forest_storm, disturbance.df_storm), 
   pattern = map(ID.forest_storm), iteration = "vector", format = "file"),
   
  # Extract results
  # -- Get functional diversity
  tar_target(FD_storm, get_FD(forest_list_storm, sim_disturbance_storm, pc1_per_species)),
  # -- Get resilience metrics
  tar_target(resilience_storm, get_resilience_metrics(
   sim_disturbance_storm, disturbance.df_storm, forest_list_storm)),
  # -- Format data together
  tar_target(data_model_storm, get_data_model(climate_storm, resilience_storm, FD_storm)),
  
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # -- Make simulations with storm disturbances and random selection -----
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  # -- generate one climate object per iteration with branching
  tar_target(climate_storm_random, make_climate(
    FUNDIV_climate_species, quantiles.in = climate_list_storm[[ID.climate_storm]], 
    "storm", 10, exclude.in = c("Carpinus_betulus", "Quercus_ilex"), 
    method = "random"), pattern = map(ID.climate_storm), iteration = "list"), 
  
  # Make species objects
  # -- First: list all species object to make
  tar_target(species_list_storm_random, make_species_list(climate_storm_random, "storm_random")),
  # -- Make a vector of ID for each species to make
  tar_target(ID.species_storm_random, species_list_storm_random$ID.species), 
  # -- Make the species via branching over ID.species
  tar_target(species_storm_random, make_species_rds(
    fit.list.allspecies, climate_storm_random, species_list_storm_random, 
    ID.species.in = ID.species_storm_random), 
    pattern = map(ID.species_storm_random), iteration = "vector", format = "file"),
  
  # Make simulations 
  # -- Start with a list of forest to simulate
  tar_target(forest_list_storm_random, make_forest_list(climate_storm_random, "storm_random")), 
  # -- Make a vector of ID for each forest to simulate
  tar_target(ID.forest_storm_random, forest_list_storm_random$ID.forest),
  # -- Make simulations till equilibrium
  tar_target(sim_equilibrium_storm_random, make_simulations_equilibrium(
    climate_storm_random, harv_rules.ref, species_list_storm_random, 
    forest_list_storm_random, species_storm_random, ID.forest_storm_random), 
    pattern = map(ID.forest_storm_random), iteration = "vector", format = "file"),
  # -- Make simulations with disturbance
  tar_target(sim_disturbance_storm_random, make_simulations_disturbance(
    climate_storm_random, harv_rules.ref, species_list_storm_random, 
    forest_list_storm_random, species_storm_random, sim_equilibrium_storm_random, 
    ID.forest.in = ID.forest_storm_random, disturbance.df_storm), 
    pattern = map(ID.forest_storm_random), iteration = "vector", format = "file"),
  
  # Extract results
  # -- Get functional diversity
  tar_target(FD_storm_random, get_FD(
    forest_list_storm_random, sim_disturbance_storm_random, pc1_per_species)),
  # -- Get resilience metrics
  tar_target(resilience_storm_random, get_resilience_metrics(
    sim_disturbance_storm_random, disturbance.df_storm, forest_list_storm_random)),
  # -- Format data together
  tar_target(data_model_storm_random, get_data_model(
    climate_storm_random, resilience_storm_random, FD_storm_random)),
  
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # -- Traits -----
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  # Traits data file
  tar_target(traits_file, "data/traits.csv", format = "file"),
  
  # Read traits data
  tar_target(traits, fread(traits_file)),
  
  # Get coordinates in the first pca axis per species
  tar_target(pc1_per_species, get_pc1_per_species(traits)),
  
  
  
   
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # -- Informative plots -----
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
  # Map of the different climates used for the analysis
  tar_target(fig_map_climates, map_climates(
   FUNDIV_climate_species, climate.in = climate_list_storm, 
   "output/fig_informative/map_climates.jpg"), format = "file"), 
  
  # Plot the pca with species traits
  tar_target(fig_traits_pca, plot_traits_pca(
    traits, "output/fig_informative/traits.jpg"), format = "file"),
  
  
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # -- Plots for analyses -----
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  # Plot the effect of FD on resilience
  tar_target(fig_FD_effect_resilience, plot_FD_effect_resilience(
    data_model_storm, "output/fig_analyses/fig_FD_effect_resilience_storm.jpg"), 
    format = "file"),
  
  # Plot the effect of FD and climate on resilience
  tar_target(fig_FD_and_climate_effect_resilience, plot_FD_and_climate_effect_resilience(
    data_model_storm, "output/fig_analyses/fig_FD_and_climate_effect_resilience_storm.jpg"), 
    format = "file"),
  
  # Plot predictions of resilience with FD metrics and climate
  tar_target(fig_predictions, plot_predictions(
    data_model_storm, "output/fig_analyses/fig_predictions_storm.jpg"), format = "file"),
  
  
  # Make a network analysis with peacewise SEM
  tar_target(fig_sem_storm_FD, plot_sem(
    data_model_storm, "FD", "output/fig_analyses/sem_storm_FD.jpg"), format = "file"), 
  
  # Plot how the FD effect changes with climate
  tar_target(fig_FD_effect_vs_climate, plot_FD_effect_vs_climate(
    data_model_storm, "output/fig_analyses/fig_fd_effect_vs_climate.jpg"), format = "file"),
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # -- Exploratory plots -----
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
  # changes in FW and CWM over time
  tar_target(fig_cwm_and_fd_overtime_storm, plot_cwm_fd_overtime(
   sim_equilibrium_storm, forest_list_storm, pc1_per_species, 
   "output/fig_exploratory/cwm_and_fd_over_time_storm.jpg"), format = "file"), 

  # Proportion of species with sensitivity estimation per climate
  tar_target(fig_prop.species_per_climate, plot_prop.species_per_climate(
    FUNDIV_climate_species, "storm", "output/fig_exploratory/fig_prop_per_clim_storm.jpg"), 
    format = "file"), 
  
  # Co-variation between resilience metrics
  tar_target(fig_covariation_FD_storm, plot_covariation_FD(
    data_model_storm, "storm", "output/fig_exploratory/fig_covar_FD_storm.jpg"), 
    format = "file"),

  # Functional density distribution of random vs selected communities
  tar_target(fig_pca1_selection_vs_random_storm, plot_pca1_selection_vs_random(
   climate_storm, pc1_per_species, 
   "output/fig_exploratory/fig_pca1_selection_vs_random_storm.jpg"), format = "file")
  
  
  
  
  
  
)

