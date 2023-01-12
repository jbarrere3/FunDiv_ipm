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
                 "cowplot", "multcomp")
for(i in 1:length(packages.in)) if(!(packages.in[i] %in% rownames(installed.packages()))) install.packages(packages.in[i])
# Targets options
options(tidyverse.quiet = TRUE, clustermq.scheduler = "multiprocess")
tar_option_set(packages = packages.in)
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
  tar_target(clim.index, c(1:5)),
  # -- list of climates
  tar_target(climate_list, create_climate_list(length(clim.index))),
  # -- generate one climate object per iteration with branching
  tar_target(
    climate,
    make_climate(FUNDIV_climate_species, 
                 quantiles.in = climate_list[[clim.index]], 
                 nsp_per_richness = 10),
    pattern = map(clim.index), 
    iteration = "list"
  ), 
  
  # Fit IPMs, one list per climate
  tar_target(IPM.list, make_IPM_multispecies(
    fit.list.allspecies[climate[[clim.index]]$species], 
    climate[[clim.index]]$climate, 
    clim_lab.in = names(climate_list)[clim.index]),
  pattern = map(clim.index), 
  iteration = "list"), 
  
  # Create species object
  tar_target(species.list, generate_species.list(IPM.list[[clim.index]]), 
             pattern = map(clim.index), iteration = "list"), 
  
  # Create forest object
  tar_target(forest.list, generate_forest_from_combinations(
    species.list[[clim.index]], harv_rules.ref, climate[[clim.index]]$combinations), 
    pattern = map(clim.index), iteration = "list"), 
  
  # Run simulation from random population to get equilibrium
  tar_target(sim.list, run_simulations_from_list(forest.list[[clim.index]], tlim = 4000), 
             pattern = map(clim.index), iteration = "list"),
  
  # Run simulations with disturbance starting at equilibrium
  tar_target(sim.list.disturbed_cold, disturb_forest.list(
    sim.list[[clim.index]], forest.list[[clim.index]], disturbance.df), 
    pattern = map(clim.index), iteration = "list")
  
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # -- Make IPMs -----
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  # # Make different climates before fitting IPM
  # tar_target(climate.cold, make_climate(FUNDIV_climate_species, 
  #                                       quantiles.in = c(0.075, 0.125))), 
  # tar_target(climate.mid, make_climate(FUNDIV_climate_species, 
  #                                      quantiles.in = c(0.475, 0.525))), 
  # tar_target(climate.hot, make_climate(FUNDIV_climate_species, 
  #                                      quantiles.in = c(0.875, 0.925))), 
  # 
  # 
  # # Fit IPM for species list from cold climate
  # tar_target(IPM.list_cold, make_IPM_multispecies(
  #   fit.list.allspecies[climate.cold$species], climate.cold$climate, 
  #   clim_lab.in = "cold")), 
  # 
  # # Fit IPM for species list from mid climate
  # tar_target(IPM.list_mid, make_IPM_multispecies(
  #   fit.list.allspecies[climate.mid$species], climate.mid$climate, 
  #   clim_lab.in = "mid")), 
  # 
  # # Fit IPM for species list from hot climate
  # tar_target(IPM.list_hot, make_IPM_multispecies(
  #   fit.list.allspecies[climate.hot$species], climate.hot$climate, 
  #   clim_lab.in = "hot")), 
  # 
  # 
  # 
  # 
  # 
  # #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # # -- Make simulations and extract outputs -----
  # #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # 
  # 
  # ## - With cold climate
  # 
  # # Make species list
  # tar_target(species.list_cold, generate_species.list(IPM.list_cold)), 
  # # Generate forest
  # tar_target(forest.list_cold, generate_forest_from_combinations(
  #   species.list_cold, harv_rules.ref, climate.cold$combinations)), 
  # # Run simulations
  # tar_target(sim.list_cold, run_simulations_from_list(
  #   forest.list_cold, tlim = 4000)), 
  # # Make simulations with disturbance starting at equilibrium
  # tar_target(sim.list.disturbed_cold, disturb_forest.list(
  #   sim.list_cold, forest.list_cold, disturbance.df)),
  # # Get resilience
  # tar_target(resilience_cold, get_resilience_metrics(
  #   sim.list.disturbed_cold, disturbance.df)),
  # 
  # 
  # ## - With mid climate
  # 
  # # Make species list
  # tar_target(species.list_mid, generate_species.list(IPM.list_mid)), 
  # # Generate forest
  # tar_target(forest.list_mid, generate_forest_from_combinations(
  #   species.list_mid, harv_rules.ref, climate.mid$combinations)), 
  # # Run simulations
  # tar_target(sim.list_mid, run_simulations_from_list(
  #   forest.list_mid, tlim = 4000)), 
  # # Make simulations with disturbance starting at equilibrium
  # tar_target(sim.list.disturbed_mid, disturb_forest.list(
  #   sim.list_mid, forest.list_mid, disturbance.df)),
  # # Get resilience
  # tar_target(resilience_mid, get_resilience_metrics(
  #   sim.list.disturbed_mid, disturbance.df)),
  # 
  # 
  # ## - With hot climate
  # 
  # # Make species list
  # tar_target(species.list_hot, generate_species.list(IPM.list_hot)), 
  # # Generate forest
  # tar_target(forest.list_hot, generate_forest_from_combinations(
  #   species.list_hot, harv_rules.ref, climate.hot$combinations)), 
  # # Run simulations
  # tar_target(sim.list_hot, run_simulations_from_list(
  #   forest.list_hot, tlim = 4000)), 
  # # Make simulations with disturbance starting at equilibrium
  # tar_target(sim.list.disturbed_hot, disturb_forest.list(
  #   sim.list_hot, forest.list_hot, disturbance.df)),
  # # Get resilience
  # tar_target(resilience_hot, get_resilience_metrics(
  #   sim.list.disturbed_hot, disturbance.df)),
  # 
  # 
  # 
  # 
  # #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # # -- Traits -----
  # #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # 
  # # Traits data file
  # tar_target(traits_file, "data/traits.csv", format = "file"),
  # 
  # # Read traits data
  # tar_target(traits, fread(traits_file)),
  # 
  # # Get coordinates in the first pca axis per species
  # tar_target(pc1_per_species, get_pc1_per_species(traits)),
  # 
  # # Extract functional diversity of species combinations in each climate
  # tar_target(FD_cold, get_FD(sim.list.disturbed_cold, pc1_per_species)),
  # tar_target(FD_mid, get_FD(sim.list.disturbed_mid, pc1_per_species)),
  # tar_target(FD_hot, get_FD(sim.list.disturbed_hot, pc1_per_species)),
  # 
  # # Format resilience and functional diversity data from different climates
  # tar_target(data_models, format_resilience_FD(
  #   list.in = list(cold = list(resilience = resilience_cold, FD = FD_cold), 
  #                  mid = list(resilience = resilience_mid, FD = FD_mid), 
  #                  hot = list(resilience = resilience_hot, FD = FD_hot))
  # )), 
  # 
  # 
  # 
  # 
  # #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # # -- Plots -----
  # #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # 
  # 
  # # Plot the co-variation between traits
  # tar_target(fig_traits_pca, plot_traits_pca(traits, "output/fig_traits.jpg"), 
  #            format = "file"),
  # 
  # # Plot the pca of sgdd and wai and the different climate selected
  # tar_target(fig_pca_climate, plot_pca_climate(FUNDIV_climate_species, climate.list = list(
  #   cold = c(0.075, 0.125),  mid = c(0.475, 0.525), hot = c(0.875, 0.925)), 
  #   file.in = "output/pca_climate.jpg"), format = "file"),
  # 
  # # Plot the effect of community weighted mean on resilience metrics
  # tar_target(fig_resilience_vs_CWM, plot_resilience_vs_CMW(
  #   data_models, "output/analyses/fig_resilience_vs_CWM.jpg"), format = "file"),
  # 
  # # Plot the effect of climate on resilience
  # tar_target(fig_resilience_vs_climate, plot_resilience_vs_climate(
  #   data_models, "output/analyses/fig_resilience_vs_climate.jpg"), format = "file"),
  # 
  # # Plot the effect of functional diversity on resilience
  # tar_target(fig_resilience_vs_FD, plot_resilience_vs_FD(
  #   data_models, "output/analyses/fig_resilience_vs_FD.jpg"), format = "file"),
  # 
  # # Plot combined effect of CWM and FD on resilience metrics
  # # -- all climates
  # tar_target(fig_res_vs_fd_and_cwm_all, plot_resilience_vs_CMW_and_FD(
  #   data_models, "output/analyses/fig_res_vs_fd_and_cwm_all.jpg"), 
  #   format = "file"),
  # # -- cold climate
  # tar_target(fig_res_vs_fd_and_cwm_cold, plot_resilience_vs_CMW_and_FD(
  #   subset(data_models, climate == "cold"), "output/analyses/fig_res_vs_fd_and_cwm_cold.jpg"), 
  #   format = "file"),
  # # -- mid climate
  # tar_target(fig_res_vs_fd_and_cwm_mid, plot_resilience_vs_CMW_and_FD(
  #   subset(data_models, climate == "mid"), "output/analyses/fig_res_vs_fd_and_cwm_mid.jpg"), 
  #   format = "file"),
  # # -- hot climate
  # tar_target(fig_res_vs_fd_and_cwm_hot, plot_resilience_vs_CMW_and_FD(
  #   subset(data_models, climate == "hot"), "output/analyses/fig_res_vs_fd_and_cwm_hot.jpg"), 
  #   format = "file"),
  # 
  # ## -- Simulations with cold climate
  # 
  # # Plot the list of simulations generated (no disturbance so far)
  # tar_target(fig_sim.list_cold, plot_sim.list(
  #   sim.list_cold, "output/cold/fig_sim_nodist_cold.jpg"), 
  #   format = "file"), 
  # # Plot the list of simulations disturbed
  # tar_target(fig_sim.list.disturbed_cold, plot_sim.list(
  #   sim.list.disturbed_cold, "output/cold/fig_sim_dist_cold.jpg"), 
  #   format = "file"),
  # 
  # 
  # 
  # ## -- Simulations with mid climate
  # 
  # # Plot the list of simulations generated (no disturbance so far)
  # tar_target(fig_sim.list_mid, plot_sim.list(
  #   sim.list_mid, "output/mid/fig_sim_nodist_mid.jpg"), 
  #   format = "file"), 
  # # Plot the list of simulations disturbed
  # tar_target(fig_sim.list.disturbed_mid, plot_sim.list(
  #   sim.list.disturbed_mid, "output/mid/fig_sim_dist_mid.jpg"), 
  #   format = "file"),
  # 
  # 
  # 
  # ## -- Simulations with hot climate
  # 
  # # Plot the list of simulations generated (no disturbance so far)
  # tar_target(fig_sim.list_hot, plot_sim.list(
  #   sim.list_hot, "output/hot/fig_sim_nodist_hot.jpg"), 
  #   format = "file"), 
  # # Plot the list of simulations disturbed
  # tar_target(fig_sim.list.disturbed_hot, plot_sim.list(
  #   sim.list.disturbed_hot, "output/hot/fig_sim_dist_hot.jpg"), 
  #   format = "file")
  
)

