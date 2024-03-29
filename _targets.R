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
                 "statmod", "xtable", "car", "modi", "grid", "gridExtra", "semEff")
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
  
  # Disturbance coefficients
  tar_target(disturb_coef.in, fread("data/disturb_coef.csv")),
  
  # Raw data from FUNDIV
  tar_target(FUNDIV_tree_file, "data/FunDiv_trees_Nadja.csv", format = "file"),
  tar_target(FUNDIV_plot_file, "data/FunDiv_plots_Nadja.csv", format = "file"),
  tar_target(FUNDIV_species_file, "data/FunDiv_species_Nadja.csv", format = "file"),
  tar_target(FUNDIV_climate_file, "data/moreno_chelsa_fundiv_clim.csv", format = "file"),
  
  # Read and format FUNDIV data
  tar_target(FUNDIV_data, read_FUNDIV(
    FUNDIV_tree_file, FUNDIV_plot_file, FUNDIV_climate_file, FUNDIV_species_file)),
  tar_target(FUNDIV_climate_species, get_FUNDIV_species_per_climate(FUNDIV_data)),
  
  # List of species for which we have data
  tar_target(all.species.name, 
             colnames(FUNDIV_climate_species)[grep("_", colnames(FUNDIV_climate_species))]), 
  
  # Get demographic parameters for all species
  tar_target(fit.list.allspecies, load_param_demo(all.species.name)),
  
  # Generate some harvest rules
  tar_target(harv_rules.ref, c(Pmax = 0.25, dBAmin = 3, freq = 5, alpha = 1)),
  
  # Data.frame containing storm disturbance to apply
  tar_target(disturbance.df_storm, data.frame(
    type = rep("storm", 3), intensity = rep(0.5, 3), 
    IsSurv = rep(FALSE, 3), t = c(500:502))),
  
  
  
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # -- Prepare data for simulations -----
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  # Generate some climates
  # -- iterations along all climates that will be created (one iteration per climate)
  tar_target(ID.climate_storm, c(1:10)),
  # -- list of climates
  tar_target(climate_list_storm, create_climate_list(length(ID.climate_storm), 
                                                     quantile.range = c(0.2, 1))),
  # -- generate one climate object per iteration with branching
  tar_target(climate_storm, make_climate(
      FUNDIV_climate_species, quantiles.in = climate_list_storm[[ID.climate_storm]], 
      "storm", 10, exclude.in = c("Carpinus_betulus", "Quercus_ilex", "Salix_caprea"), 
      method = "frequency", disturb_coef.in, traits), 
      pattern = map(ID.climate_storm), iteration = "list"), 
  
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
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # -- Make simulations -----
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
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
   climate_storm, harv_rules.ref, species_list_storm, forest_list_storm, disturb_coef.in,
   species_storm, sim_equilibrium_storm, ID.forest.in = ID.forest_storm, disturbance.df_storm), 
   pattern = map(ID.forest_storm), iteration = "vector", format = "file"),
   
  
  
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # -- Prepare traits data -----
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  # -- Traits data file
  tar_target(traits_file, "data/traits.csv", format = "file"),
  # -- Read traits data
  tar_target(traits, fread(traits_file)),
  # -- Get recruitment traits
  tar_target(traits_rec, get_recruitment_traits(
    fit.list.allspecies, FUNDIV_data, "specific")), 
  # -- Get all species included in the simulations for filtering
  tar_target(sp.in.sim, unique(species_list_storm$species)),
  # -- Get coordinates on the pca axis per species
  tar_target(pca12_per_species, get_pc12_per_species(
    subset(traits, species %in% sp.in.sim), 
    subset(traits_rec, species %in% sp.in.sim))), 
  
  
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # -- Extract results form the simulations -----
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  # -- Get functional diversity
  tar_target(FD_multivar, get_FD_multivar(
    forest_list_storm, sim_disturbance_storm, pca12_per_species)),
  # -- Get resilience metrics
  tar_target(resilience_storm, get_resilience_metrics(
    sim_disturbance_storm, disturbance.df_storm, forest_list_storm)),
  # -- Format data together
  tar_target(data_model_all, get_data_model(climate_storm, resilience_storm, FD_multivar)),
  # -- Only keep simulations that reached equilibrium
  tar_target(data_model, subset(data_model_all, (SD < 0.15))),
  
  
   
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # -- Plots and tables for the methods -----
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
  # Map of the different climates used for the analysis
  tar_target(fig_map_climates, map_climates(
   FUNDIV_climate_species, climate.in = climate_list_storm, 
   "output/methods/map_climates.jpg"), format = "file"), 
  
  # Plot the pca with species traits
  tar_target(fic_pca12, plot_traits_pca12(
    subset(traits, species %in% sp.in.sim), subset(traits_rec, species %in% sp.in.sim), 
    species_list_storm, "output/methods/pca_traits.jpg"), format = "file"),
  
  # Plot of the resilience metrics
  tar_target(fig_metrics, plot_metrics(
    sim_disturbance_storm[64], "output/methods/metrics.jpg"), format = "file"),
  
  # Export species information in a table
  tar_target(table_species, export_species_table(
    FUNDIV_data, pca12_per_species, species_list_storm, 
    "output/methods/species_table.tex"), format = "file"), 
  
  
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # -- Main analyses -----
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  # Plot the effect of FD on resilience
  tar_target(fig_FD_effect_resilience_multivar, plot_FD_effect_resilience_multivar(
    data_model, "output/analyses/H1"), format = "file"),
  
  
  # Plot how the FD effect changes with climate
  tar_target(fig_FD_effect_vs_climate_multivar, plot_FD_effect_vs_climate_multivar(
    data_model, "H", "AIC", "output/analyses/H2"), format = "file"), 
  
  # Make a network analysis with peacewise SEM
  # -- Fit sem model
  tar_target(mod_sem, fit_sem(data_model)),
  # -- Make the plot
  tar_target(fig_sem_multivar, plot_sem_multivar(
    mod_sem, "output/analyses/H3"), format = "file"), 
  
  
  
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # -- Plots for supplementary -----
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
  # Co-variation between resilience metrics
  tar_target(fig_covariation_FD_storm, plot_covariation(
    data_model, "output/supplementary/fig_covar_storm.png")),
  
  # Relation between diversity and density of trees
  tar_target(fig_tree_packing, plot_tree_packing(
    data_model, "output/supplementary/fig_treepacking.jpg"), format = "file"), 
  
  # Relation between climate, diversity and structure
  tar_target(fig_climate_diversity_structure, plot_climate_vs_diversity_and_structure(
    data_model, "output/supplementary/fig_climate_vs_diversity_and_structure.jpg"), 
    format = "file"), 
  
  # SEM with additional relations
  # -- Fit sem model
  tar_target(mod_sem_supp, fit_sem_supp(data_model)),
  # -- Make the plot
  tar_target(fig_sem_multivar_supp, plot_sem_multivar(
    mod_sem_supp, "output/supplementary"), format = "file"), 
  
  # Plot some simulations
  tar_target(fig_simulations, plot_sim_dist(
    sim_disturbance_storm[data_model$ID.forest[c(401:412)]], 
    "output/supplementary/simulations.jpg"), format = "file"), 
  tar_target(fig_simulations_eq, plot_sim_dist(
    sim_equilibrium_storm[data_model$ID.forest[c(401:412)]], 
    "output/supplementary/simulations_eq.jpg"), format = "file"), 
  
  # Plot diversity and structure of raw data
  tar_target(fig_div_str_data, plot_clim_vs_div_and_str_data(
    FUNDIV_data, pca12_per_species, climate_list_storm, 
    "output/supplementary/div_str.jpg"), format = "file"), 

  
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # -- Plots and analyses for revision -----
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  
  # Extract species composition at equilibrium, post-dist and post-recov. 
  tar_target(sp_composition_dist, get_sp_composition(
    sim_disturbance_storm, sp.in.sim, c(1, max(disturbance.df_storm$t)+1, 4000))), 
  
  # Plot changes in species composition
  tar_target(fig_spcomposition, plot_spcomposition(
    sp_composition_dist, data_model, "output/supplementary"), format = "file"), 
  
  # Plot the timing of disturbance and recovery
  tar_target(fig_timing_dist, plot_timing_dist(
    sim_disturbance_storm[64], disturbance.df_storm, 
    "output/supplementary/fig_timing.jpg"), format = "file")
  
  
  
  
)

