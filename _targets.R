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
packages.in <- c("dplyr", "ggplot2", "treeforce", "tidyr", "data.table", 
                 "factoextra", "modi")
for(i in 1:length(packages.in)) if(!(packages.in[i] %in% rownames(installed.packages()))) install.packages(packages.in[i])
# Targets options
options(tidyverse.quiet = TRUE)
tar_option_set(packages = packages.in)
set.seed(2)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Targets workflow --------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

list(
  
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # -- Load and format traits data -----
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  # Traits data file
  tar_target(traits_file, "data/traits.csv", format = "file"),
  
  # Read traits data
  tar_target(traits, fread(traits_file)),
  
  # Get coordinates in the first pca axis per species
  tar_target(pc1_per_species, get_pc1_per_species(traits)),
  
  # Plot the co-variation between traits
  tar_target(fig_traits_pca, plot_traits_pca(traits, "output/fig_traits.jpg"), 
             format = "file"),
  
  
  
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # -- Load demographic data and parameters -----
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  # Species for which to build IPM
  tar_target(species.names, c("Fagus_sylvatica", "Picea_abies", "Abies_alba", 
                              "Quercus_robur", "Pinus_sylvestris")), 
  
  # Load demographic parameters for each species
  tar_target(fit.list, load_param_demo(species.names)), 
  
  # Get the optimum climate of the first species as reference
  tar_target(climate.ref, load_and_format_climate_ipm(species.names[1])), 
  
  # Generate some harvest rules
  tar_target(harv_rules.ref, c(Pmax = 0.25, dBAmin = 3, freq = 10, alpha = 1)),
  
  # Data.frame containing disturbances to apply
  tar_target(disturbance.df, data.frame(
    type = rep("storm", 3), intensity = rep(0.5, 3), 
    IsSurv = rep(FALSE, 3), t = c(500:502))), 
  
  
  
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # -- Make IPMs -----
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  # Fit IPM for all species in the list
  tar_target(IPM.list, make_IPM_multispecies(fit.list, climate.ref, clim_lab.in = "climate_Fag.syl")), 
  
  # Fit IPM for all species in the list and three different climates
  tar_target(IPM.list_multiclim, make_IPM_multisp_multiclim(fit.list)), 
  
  
  
  
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # -- Make simulations -----
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  ## -- With one reference climate
  
  # Create a list of species object to make simulations
  tar_target(species.list, generate_species.list(IPM.list)), 
  
  # Generate a list of forest with one, two or three species
  tar_target(forest.list, generate_forest_list(species.list, harv_rules.ref)), 
  
  # Run simulations from the list of forests generated
  tar_target(sim.list, run_simulations_from_list(forest.list, tlim = 3000)), 
  
  # Make simulations with disturbance starting at equilibrium
  tar_target(sim.list.disturbed, disturb_forest.list(
    sim.list, forest.list, disturbance.df)),
  
  # Extract resilience and functional diversity from simulations
  tar_target(FD_and_resilience, get_resilience_and_FD(
    sim.list.disturbed, pc1_per_species)),
  
  
  
  ## -- With several climate per species
  
  # Create a list of species object to make simulations
  tar_target(species.list_multiclim, generate_species.list_multiclim(IPM.list_multiclim)), 
  
  # Generate a list of forest with one species but multiple climates
  tar_target(forest.list_multiclim, generate_forest_list_multiclim(
    species.list_multiclim, maxsp = 1, harv_rules.ref)), 
  
  # Run simulations from the list of species generated
  tar_target(sim.list_multiclim, run_simulations_from_list_multiclim(
    forest.list_multiclim, tlim = 5000)), 
  
  # Make simulations with disturbance starting at equilibrium
  tar_target(sim.list.disturbed_multiclim, disturb_forest.list_multiclim(
    sim.list_multiclim, forest.list_multiclim, disturbance.df)), 
  
  
  
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # -- Plot results -----
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  ## -- With one reference climate
  
  # Plot the list of simulations generated (no disturbance so far)
  tar_target(fig_sim.list, plot_sim.list(sim.list, "output/fig_sim_nodist.jpg"), 
             format = "file"), 
  
  # Plot the list of simulations disturbed
  tar_target(fig_sim.list.disturbed, plot_sim.list(sim.list.disturbed, 
                                                   "output/fig_sim_dist.jpg"), 
             format = "file"),
  
  # FD vs resilience
  tar_target(fig_FD_vs_resilience, plot_FD_vs_resilience(
    FD_and_resilience, "output/fig_FD_vs_resilience.jpg"), format = "file"),
  
  
  
  ## -- With several climate per species
  
  # Create a list of species object to make simulations
  tar_target(fig_sim.list_multiclim, plot_sim.list_multiclim(
    sim.list_multiclim, "output/fig_simulticlim_nodist.jpg"), format = "file"), 
  
  # Plot the list of simulations disturbed
  tar_target(fig_sim.list_multiclim.disturbed, plot_sim.list_multiclim(
    sim.list.disturbed_multiclim, "output/fig_simulticlim_dist.jpg"), format = "file")
  
)

