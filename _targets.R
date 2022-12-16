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
  
  # Disturbance parameters to load and format 
  tar_target(disturbance_parameters_file, "data/parameters_disturbance.Rdata", format = "file"),
  tar_target(disturbance_parameters, load_and_format_dist.param(disturbance_parameters_file)),
  tar_target(disturbance_parameters_alliter, load_and_format_dist.param(disturbance_parameters_file, 
                                                                        calc.type = "all")),
  
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
  
  # -- Generate a list of forest with one, two or three species
  tar_target(forest.list, generate_forest_list(species.list, harv_rules.ref)), 
  
  # Run simulations from the list of forests generated
  tar_target(sim.list, run_simulations_from_list(forest.list, tlim = 3000)), 
  
  
  
  # Plot comparison between disturbance based on mean parameter or all iteration
  # tar_target(fig_meanParam_vs_meanProba, plot_meanProba_vs_meanParam(
  #   species.names, disturbance.time = 1500, disturbance.in = "storm", intensity.in = 0.8, IPM.list, 
  #   duration.disturbance = 5, disturbance_parameters, disturbance_parameters_alliter, 
  #   "output/fig_meanProba_vs_meanParam.jpg"), format = "file"), 
  
  # Plot disturbance applied to monospecific stands
  tar_target(fig_monosp_disturbed, plot.monosp.disturbances(
    species.list, disturbance_parameters, disturbance.in = "storm", intensity.in = 0.8, 
    duration.disturbance = 5, sim.time = 5000, file.in = "output/fig_disturb_monosp.jpg"), 
    format = "file"), 
  
  # Plot the list of simulations generated (no disturbance so far)
  tar_target(fig_sim.list, plot_sim.list(sim.list, "output/fig_sim_nodist.jpg"), 
             format = "file")
  
)

