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
packages.in <- c("dplyr", "ggplot2", "treeforce", "tidyr")
for(i in 1:length(packages.in)) if(!(packages.in[i] %in% rownames(installed.packages()))) install.packages(packages.in[i])
# Targets options
options(tidyverse.quiet = TRUE)
tar_option_set(packages = packages.in)
set.seed(2)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Targets workflow --------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

list(
  
  # Species for which to build IPM
  tar_target(species.names, c("Fagus_sylvatica", "Picea_abies", "Abies_alba", 
                              "Quercus_robur", "Pinus_sylvestris")), 
  
  # Load demographic parameters for each species
  tar_target(fit.list, load_param_demo(species.names)), 
  
  # Get the optimum climate of the first species as reference
  tar_target(climate.ref, load_and_format_climate_ipm(species.names[1])), 
  
  # Fit IPM for all species in the list
  tar_target(IPM.list, make_IPM_multispecies(fit.list, climate.ref, clim_lab.in = "climate_Fag.syl")), 
  
  # Generate some harvest rules
  tar_target(harv_rules.ref, c(Pmax = 0.25, dBAmin = 3, freq = 10, alpha = 1)),
  
  # Create a list of species object to make simulations
  tar_target(species.list, generate_species.list(IPM.list)), 
  
  # Generate a list of forest with one, two or three species
  tar_target(forest.list, generate_forest_list(species.list, harv_rules.ref)), 
  
  # Run simulations from the list of forests generated
  tar_target(sim.list, run_simulations_from_list(forest.list, tlim = 3000)), 
  
  # Data.frame containing disturbances to apply
  tar_target(disturbance.df, data.frame(type = "storm", intensity = 0.5, 
                                        IsSurv = FALSE, t = 100)), 
  
  # Make simulations with disturbance starting at equilibrium
  tar_target(sim.list.disturbed, disturb_forest.list(
    sim.list, forest.list, disturbance.df)),
  
  
  # Plot the list of simulations generated (no disturbance so far)
  tar_target(fig_sim.list, plot_sim.list(sim.list, "output/fig_sim_nodist.jpg"), 
             format = "file"), 
  
  # Plot the list of simulations disturbed
  tar_target(fig_sim.list.disturbed, plot_sim.list(sim.list.disturbed, 
                                                   "output/fig_sim_dist.jpg"), 
             format = "file")
  
)

