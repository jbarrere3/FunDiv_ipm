#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#### SCRIPT INTRODUCTION ####
#
#' @name functions_data.R  
#' @description R script containing all functions relative to data
#               importation and formatting
#' @author Julien Barrere, Marianne Bernard
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



#' Function to get the path of a file, and create directories if they don't exist
#' @param file.in character: path of the file, filename included (ex: "plot/plot.png")
create_dir_if_needed <- function(file.in){
  
  path.in <- strsplit(file.in, "/")[[1]]
  if(length(path.in) > 1){
    for(i in 1:(length(path.in)-1)){
      if(i == 1) path.in_i <- path.in[i]
      else path.in_i <- paste(path.in_i, path.in[i], sep = "/")
      if(!dir.exists(path.in_i)) dir.create(path.in_i)
    }
  }
}


#' Function to load and format disturbance parameters per species
#' @param disturbance_parameters_file file containing parameters value 
#' @param calc.type character indicating whether to use mean parameter value ("mean")
#'                  or all posterior iterations ("all")
load_and_format_dist.param <- function(disturbance_parameters_file, calc.type = "mean"){
  
  # Load the rdata file
  load(disturbance_parameters_file)
  
  # Remove unused disturbance types and species
  out <- param_per_iteration %>%
    filter(!(species %in% c("Other broadleaf", "Other conifer"))) %>%
    filter(disturbance %in% c("storm", "fire", "biotic"))
  
  # If the calculation type is mean, average parameters over all mcmc iteration
  if(calc.type == "mean"){
    out <- out %>%
      gather(key = "parameter", value = "value", "a0", "a1", "b", "c", "dbh.intercept", 
             "dbh.slope", "logratio.intercept", "logratio.slope") %>%
      group_by(disturbance, species, parameter) %>%
      summarize(value.mean = mean(value)) %>%
      spread(key = "parameter", value = "value.mean")
  }
  
  # Return output
  return(out)
}

#' Function to load demographic parameters for several species
#' @param species.names character vector of species ("Genus_species" format)
load_param_demo <- function(species.names){
  
  # Initialize list of fits (one element per species)
  fit.list <- list()
  
  # Loop on all species to load and store demographic parameters
  for(i in 1:length(species.names)){
    eval(parse(text=paste0("fit.list$", species.names[i], " <- fit_", species.names[i])))
  }
  
  # Return the list
  return(fit.list)
  
}

#' Function to load species climate data and format it as input for the IPM
#' @param species.name Name of the species foor which to load climate (Genus_species format)
#' @param N.ref Integer specifying whether to keep species cold margin (1), optimum (2) or hot margin (3)
load_and_format_climate_ipm <- function(species.name, N.ref = 2){
  
  # Load the data
  data("climate_species")
  
  # Format to fit IPM 
  out  <- drop(as.matrix(
    climate_species %>%
      filter(N == N.ref & sp == species.name) %>%
      dplyr::select(-sp)
  )) 
  
  # Return formatted vector
  return(out)
}



#' Fit IPM for several species and one climate
#' @param fit.list List containing the demographic parameters of each species
#' @param climate.ref Climate to use for each fit of the IPM
#' @param clim_lab.in Name of the reference climate
#' @param mesh.m numeric indicating the mesh size
#' @param mesh.L numeric indicating the lower bound of size
#' @param BA.max maximum basal area for the integration
make_IPM_multispecies <- function(fit.list, climate.ref, clim_lab.in, 
                                  mesh.m = 700, mesh.L = 100, BA.max = 200){
  
  # Identify the species names
  species.names <- names(fit.list)
  
  # Initialize the list of IPM
  IPM.list <- list()
  
  # Loop on all species
  for(i in 1:length(species.names)){
    # Print the species 
    print(paste0("fit IPM for species ", i, "/", length(species.names), " : ", gsub("\\_", "\\ ", species.names[i])))
    # Make the IPM
    ipm.i <- make_IPM(
      species = species.names[i], climate = climate.ref, fit =  fit.list[[i]],
      clim_lab = clim_lab.in,
      mesh = c(m = mesh.m, L = mesh.L, U = as.numeric(fit.list[[i]]$info[["max_dbh"]]) * 1.1),
      BA = 0:BA.max, verbose = TRUE
    )
    # Add to the list
    eval(parse(text=paste0("IPM.list$", species.names[i], " <- ipm.i")))
  }
  
  # Return the list generated
  return(IPM.list)
}




#' Generate a list of species object from the list of IPM
#' @param IPM.list list of fitted IPM
#' @param f.init Function to initialize basal area
generate_species.list <- function(IPM.list, f.init = def_initBA(20)){
  
  # Names of the species
  species.names <- names(IPM.list)
  
  # Initialize the list of species
  species.list <- list()
  
  # Loop on all IPM (and thus on all species)
  for(i in 1:length(species.names)){
    # Generate species object for species i
    species.i <- new_species(IPM.list[[i]], init_pop = f.init,
                             harvest_fun = Uneven_harv)
    # Add it to the list
    eval(parse(text=paste0("species.list$", species.names[i], " <- species.i")))
  }
  
  # Return the output list
  return(species.list)
}


#' Function to run simulations for a list of forests
#' @param forest.list List of forest objects
#' @param tlim Number of simulation iterations (in years)
#' @param targetBA basal area targetted by the harvest module (in m2)
#' @param SurfEch Sampled area in ha
run_simulations_from_list <- function(forest.list, tlim = 2000, targetBA = 40, 
                                      SurfEch = 0.03){
  
  # Initialize output list
  list.out <- list()
  
  # Loop on all forests in the list
  for(i in 1:length(names(forest.list))){
    
    # Printer
    print(paste0("Running simulation ", i, "/", length(names(forest.list)), 
                 " - : ", gsub("\\.", "\\ and\\ ", names(forest.list)[i])))
    
    # Run simulation without harvest
    sim.i <- sim_deter_forest(
      forest.list[[i]], tlim = tlim, equil_time = tlim, equil_dist = 1,
      harvest = "default", targetBA = targetBA, SurfEch = SurfEch, verbose = FALSE
    )
    
    # Add simulation to the output list
    eval(parse(text=paste0("list.out$", names(forest.list)[i], " <- sim.i")))
    
  }
  
  # Return the output list
  return(list.out)
}


#' Generat a list oof forest, based on all combinations of one, two or three 
#' species from a list of species objects
#' @param species.list List of species objects
#' @param harv_rules.ref List of rules for harvesting
generate_forest_list <- function(species.list, harv_rules.ref){
  
  # Identify the list of species
  species.in <- names(species.list)
  
  # All combinations of mono specific stands
  onesp <- species.in
  # All combinations of two-species stands
  twosp <- (expand.grid(sp1 = species.in, sp2 = species.in) %>%
              filter(sp1 != sp2) %>%
              mutate(sp = paste(sp1, sp2, sep = ".")))$sp
  # All combinations of three-species stands
  threesp <- (expand.grid(sp1 = species.in, sp2 = species.in, sp3 = species.in) %>%
                filter(sp1 != sp2 & sp2 != sp3 & sp1 != sp3) %>%
                mutate(sp = paste(sp1, sp2, sp3, sep = ".")))$sp
  # All species combinations
  sp.combinations <- c(onesp, twosp, threesp)
  
  # Final list of combinations (without same combination but different order)
  sp.combinations.final <- c()
  for(i in 1:length(sp.combinations)){
    combi.i <- paste(
      unlist(strsplit(sp.combinations[i], "\\."))[order(unlist(strsplit(sp.combinations[i], "\\.")))], 
      collapse = "."
    )
    if(!(combi.i %in% sp.combinations.final)) sp.combinations.final <- c(sp.combinations.final, combi.i)
  }
  # Initialize forest list
  list.out <- list()
  
  # Loop on all species combinations
  for(i in 1:length(sp.combinations.final)){
    
    # Vector of species combination i 
    species.i <- unlist(strsplit(sp.combinations.final[i], "\\."))
    
    # Generate forest i
    forest.i <- new_forest(species = species.list[species.i], harv_rules = harv_rules.ref)
    
    # Add forest i to the output list
    eval(parse(text=paste0("list.out$", sp.combinations.final[i], " <- forest.i")))
  }
  
  # Return the output list
  return(list.out)
}

#' Function to disturb the population at the last iteration of a simulation
#' @param sim simulation object
#' @param disturbance.in which disturbance to apply
#' @param intensity.in intensity (between 0 and 1) of the disturbance
#' @param duration.disturbance how long should the disturbance kill trees 
#' @param disturbance_parameters Parameters averaged over all iterations
disturb.population <- function(sim, disturbance.in, intensity.in, 
                               duration.disturbance, disturbance_parameters){
  
  # Identify the last year of simulation
  disturbance.time <- dim(sim)[2] - 1
  
  # Identify the species in the simulation
  species.in <- unique(gsub("\\..+", "", rownames(sim)))
  
  # Extract vector to disturb 
  pop.to.disturb <- tree_format(sim) %>%
    # Restrict to last year post disturbance and to tree number per size class
    filter(time == disturbance.time & var == "m" & species %in% species.in) %>%
    dplyr::select(-equil) %>%
    distinct()
  
  # Calculate quadratic diameter at time of disturbance
  dqm.disturbance.time <- (pop.to.disturb %>%
                             group_by(size) %>%
                             summarize(value = sum(value)) %>%
                             ungroup() %>%
                             summarize(dqm = sqrt(sum(size*size*value)/sum(value))))$dqm
  
  # Disturb the population with the mean of parameter values
  disturbed.proba <- pop.to.disturb %>%
    # Add disturbance parameters
    mutate(species = gsub("\\_", "\\ ", species), 
           disturbance = disturbance.in) %>%
    left_join(disturbance_parameters, by = c("species", "disturbance")) %>%
    # Calculate variables used to calculate probability
    mutate(dqm = dqm.disturbance.time, 
           logratio = log(size/dqm), 
           dbh.scaled = dbh.intercept + size*dbh.slope, 
           logratio.scaled = logratio.intercept + logratio*logratio.slope, 
           I = intensity.in, 
           t = duration.disturbance) %>%
    # Calculate probability of survival from disturbance per size class
    mutate(p = (1 - plogis(a0 + a1*logratio.scaled + b*I^(c*dbh.scaled)))^t) %>%
    dplyr::select(size, value, p) %>%
    mutate(distribution.disturbed = value*p)
  
  # Length of the mesh
  mesh <- (tree_format(sim) %>% summarize(mesh.max = max(mesh, na.rm = TRUE)))$mesh.max
  
  # Species distribution disturbed
  distr.before.disturbance <- sim[grep(paste0(species.in,".m"), rownames(sim)),
                                  grep("eqt", colnames(sim))]
  distr.disturbed <- distr.before.disturbance
  distr.disturbed[c(1:mesh)] <- disturbed.proba$distribution.disturbed
  
  # Return disturbed population
  return(distr.disturbed)
}

