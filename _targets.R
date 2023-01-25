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
                 "cowplot", "multcomp", "piecewiseSEM")
for(i in 1:length(packages.in)) if(!(packages.in[i] %in% rownames(installed.packages()))) install.packages(packages.in[i])
# Targets options
options(tidyverse.quiet = TRUE, clustermq.scheduler = "multiprocess")
tar_option_set(packages = packages.in,
               memory = "transient")
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
  
  # Data.frame containing disturbances to apply
  tar_target(disturbance.df, data.frame(
    type = rep("storm", 3), intensity = rep(0.5, 3), 
    IsSurv = rep(FALSE, 3), t = c(500:502))),
  
  
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # -- Make simulations -----
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  # Generate some climates
  # -- iterations along all climates that will be created (one iteration per climate)
  tar_target(ID.climate, c(1:10)),
  # -- list of climates
  tar_target(climate_list, create_climate_list(length(ID.climate))),
  # -- generate one climate object per iteration with branching
  tar_target(climate, make_climate(
      FUNDIV_climate_species, quantiles.in = climate_list[[ID.climate]], 
      nsp_per_richness = 10), pattern = map(ID.climate), iteration = "list"), 
  
  # Make species objects
  # -- First: list all species object to make
  tar_target(species_list, make_species_list(climate)),
  # -- Make a vector of ID for each species to make
  tar_target(ID.species, species_list$ID.species), 
  # -- Make the species via branching over ID.species
  tar_target(species, make_species_rds(
    fit.list.allspecies, climate, species_list, ID.species.in = ID.species), 
    pattern = map(ID.species), iteration = "vector", format = "file"),
  
  # Make simulations 
  # -- Start with a list of forest to simulate
  tar_target(forest_list, make_forest_list(climate)), 
  # -- Make a vector of ID for each forest to simulate
  tar_target(ID.forest, forest_list$ID.forest),
  # -- Make simulations till equilibrium
  tar_target(sim_equilibrium, make_simulations_equilibrium(
    climate, harv_rules.ref, species_list, forest_list, species, ID.forest), 
    pattern = map(ID.forest), iteration = "vector", format = "file"),
  # -- Make simulations with disturbance
  tar_target(sim_disturbance, make_simulations_disturbance(
    climate, harv_rules.ref, species_list, forest_list, species, sim_equilibrium, 
    ID.forest.in = ID.forest, disturbance.df), 
    pattern = map(ID.forest), iteration = "vector", format = "file"),
  
  # Extract results
  # -- Get functional diversity
  tar_target(FD, get_FD(forest_list, sim_disturbance, pc1_per_species)),
  # -- Get resilience metrics
  tar_target(resilience, get_resilience_metrics(
    sim_disturbance, disturbance.df, forest_list)),
  # -- Format data together
  tar_target(data_model, get_data_model(climate, resilience, FD)),
 
  
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # -- Traits -----
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  # Traits data file
  tar_target(traits_file, "data/traits.csv", format = "file"),
  
  # Read traits data
  tar_target(traits, fread(traits_file)),
   
  # Get coordinates in the first pca axis per species
  tar_target(pc1_per_species, get_pc1_per_species(traits)),
  
  # Plot the pca with species traits
  tar_target(fig_traits_pca, plot_traits_pca(traits, "output/traits.jpg"), 
             format = "file"),
  
  
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # -- Plots -----
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  # Map of the different climates used for the analysis
  tar_target(fig_map_climate, map_climates(
    FUNDIV_climate_species, climate_list, "output/map_climate.jpg"), 
    format = "file"), 
  
  # Plot the position of each plot and climate along the sgdd - wai pca
  tar_target(fig_pca_climate, plot_pca_climate(
    FUNDIV_climate_species, climate_list, "output/pca_climate.jpg"), 
    format = "file"), 
  
  # Plot resilience vs FD anc CWM
  tar_target(fig_resilience_vs_cmw_and_fd, plot_resilience_vs_CMW_and_FD(
    data_model, "output/fig_resilience_vs_CMW_and_FD.jpg"), format = "file"),
  
  # Plot effect of FD and CWM on resilience along climatic gradient
  tar_target(fig_CMW_and_FD_effect_climate, plot_CMW_and_FD_effect_climate(
    data_model, "output/fig_CMW_and_FD_effect_climate.jpg"), format = "file"),
  
  # Plot resilience vs climate
  tar_target(fig_resilience_vs_climate, plot_resilience_vs_climate(
    data_model, "output/fig_resilience_vs_climate.jpg"), format = "file"), 
  
  # Plot FD and CWM along the climatic gradient
  tar_target(fig_FD_and_CWM_vs_climate, plot_FD_and_CWM_vs_climate(
    data_model, "output/fig_FD_and_CWM_vs_climate.jpg"), 
    format = "file"), 
  
  # Make a network analysis with peacewise SEM
  tar_target(fig_sem, plot_sem(data_model, "output/fig_sem.jpg"), 
             format = "file")
  
  

 
  
)

