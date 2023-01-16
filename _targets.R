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
  tar_target(ID.climate, c(1:5)),
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
  
  # Make simulations to reach equilibrium
  # -- Start with a list of forest to simulate
  tar_target(forest_list, make_forest_list(climate)), 
  # -- Make a vector of ID for each forest to simulate
  tar_target(ID.forest, forest_list$ID.forest),
  # # Make IPM
  # # -- generate a list of IPM
  # tar_target(IPM_list, make_IPM_list(climate)),
  # # -- vector of IPM id
  # tar_target(ID.IPM, IPM_list$ID.IPM),
  # # -- Fit IPMs, one per element of IPM_list
  # tar_target(IPMs,make_IPM(
  #   species = IPM_list$species[ID.IPM], 
  #   climate = climate[[IPM_list$ID.climate[ID.IPM]]]$climate, 
  #   fit =  fit.list.allspecies[[IPM_list$species[ID.IPM]]],
  #   clim_lab = as.character(IPM_list$ID.climate[ID.IPM]),
  #   mesh = c(m = 700, L = 100, U = as.numeric(
  #     fit.list.allspecies[[IPM_list$species[ID.IPM]]]$info[["max_dbh"]]) * 1.1),
  #   BA = 0:200, verbose = TRUE
  # ),
  # pattern = map(ID.IPM), 
  # iteration = "list"),
  # 
  # # Create species object
  # tar_target(species, species(
  #   IPMs[[ID.IPM]], init_pop = def_initBA(20), harvest_fun = def_harv, 
  #   disturb_fun = def_disturb), pattern = map(ID.IPM), iteration = "list"), 
  # 
  # # Create the forests
  # # -- First, make a table indicating species combination and climate per forest
  # tar_target(forest_list, make_forest_list(climate, IPM_list)), 
  # # -- Then create a vector with forest ID for branching
  # tar_target(ID.forest, forest_list$ID.forest), 
  # # -- Finally, make one forest per ID
  # tar_target(forests, new_forest(
  #   species = species[as.numeric(unlist(strsplit(forest_list$ID.IPM_code[ID.forest], "\\.")))], 
  #   harv_rules = harv_rules.ref), pattern = map(ID.forest), iteration = "list"), 
  # 
  # # Run simulations without disturbance to reach equilibrium 
  # tar_target(sim.equilibrium, sim_deter_forest(
  #   forests[[ID.forest]], tlim = 4000, equil_time = 10000, equil_dist = 250, 
  #   equil_diff = 1, harvest = "default", SurfEch = 0.03, verbose = TRUE), 
  #   pattern = map(ID.forest), iteration = "list"),
  # 
  # # Run simulations with disturbances
  # # -- First, get the equilibrium for each simulation and verify that it was reached
  # tar_target(equil, get_equil(sim.equilibrium[[ID.forest]]), 
  #            pattern = map(ID.forest), iteration = "list"), 
  # # -- Then run simulations with disturbance from that equilibrium
  # tar_target(sim.disturbance, make_simulation_disturbance(
  #   forests[[ID.forest]], equil[[ID.forest]], disturbance.df, time.sim = 3000), 
  #   pattern = map(ID.forest), iteration = "list"), 
  # 
  # # Extract resilience metrics
  # tar_target(resilience_metrics, get_resilience_metrics(
  #   sim.disturbance, disturbance.df, forest_list)),

  
   
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # -- Traits -----
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  # Traits data file
  tar_target(traits_file, "data/traits.csv", format = "file"),
  
  # Read traits data
  tar_target(traits, fread(traits_file)),
   
  # Get coordinates in the first pca axis per species
  tar_target(pc1_per_species, get_pc1_per_species(traits))
   
  # Extract functional diversity and cwm per forest simulated
  # tar_target(FD, get_FD(ID.forest, sim.disturbance, pc1_per_species)), 
  # 
  # # Format data before running the models
  # tar_target(data_models, get_data_model(climate, resilience_metrics, FD))
  
  

 
  
)

