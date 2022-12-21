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


#' Get coordinate in first pca traits axis per species
#' @param traits dataframe containing trait value per species
get_pc1_per_species <- function(traits){
  
  # - Make PCA 
  pca <- prcomp((traits %>% dplyr::select(-species)), 
                center = T, scale = T)
  
  # - Extract the coordinates of the ndividuals on pca axis
  out <- data.frame(species = traits$species, 
                    pca1 = get_pca_ind(pca)[[1]][, 1]) 
  
  # return the output
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



#' Fit IPM for several species and one climate
#' @param fit.list List containing the demographic parameters of each species
#' @param mesh.m numeric indicating the mesh size
#' @param mesh.L numeric indicating the lower bound of size
#' @param BA.max maximum basal area for the integration
make_IPM_multisp_multiclim <- function(fit.list, mesh.m = 700, mesh.L = 100, BA.max = 200){
  
  # Identify the species names
  species.names <- names(fit.list)
  
  # Initialize the list of IPM
  IPM.list <- list()
  
  # Loop on all species
  for(i in 1:length(species.names)){
    # Print the species 
    print(paste0("fit IPM for species ", i, "/", length(species.names), " : ", gsub("\\_", "\\ ", species.names[i])))
    
    # For cold margin
    print("--- cold edge")
    climate.cold.i <- load_and_format_climate_ipm(species.names[i], N.ref = 1)
    ipm.cold.i <- make_IPM(
      species = species.names[i], climate = climate.cold.i, fit =  fit.list[[i]],
      clim_lab = paste0(species.names[i], "_cold"),
      mesh = c(m = mesh.m, L = mesh.L, U = as.numeric(fit.list[[i]]$info[["max_dbh"]]) * 1.1),
      BA = 0:BA.max, verbose = FALSE
    )
    
    # For climate optimum
    print("--- optimum")
    climate.optimum.i <- load_and_format_climate_ipm(species.names[i], N.ref = 2)
    ipm.optimum.i <- make_IPM(
      species = species.names[i], climate = climate.optimum.i, fit =  fit.list[[i]],
      clim_lab = paste0(species.names[i], "_optimum"),
      mesh = c(m = mesh.m, L = mesh.L, U = as.numeric(fit.list[[i]]$info[["max_dbh"]]) * 1.1),
      BA = 0:BA.max, verbose = FALSE
    )
    
    # For hot margin
    print("--- hot edge")
    climate.hot.i <- load_and_format_climate_ipm(species.names[i], N.ref = 3)
    ipm.hot.i <- make_IPM(
      species = species.names[i], climate = climate.hot.i, fit =  fit.list[[i]],
      clim_lab = paste0(species.names[i], "_hot"),
      mesh = c(m = mesh.m, L = mesh.L, U = as.numeric(fit.list[[i]]$info[["max_dbh"]]) * 1.1),
      BA = 0:BA.max, verbose = FALSE
    )
    
    # Add to the list
    eval(parse(text=paste0("IPM.list$", species.names[i], "$cold <- ipm.cold.i")))
    eval(parse(text=paste0("IPM.list$", species.names[i], "$optimum <- ipm.optimum.i")))
    eval(parse(text=paste0("IPM.list$", species.names[i], "$hot <- ipm.hot.i")))
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
    species.i <- species(IPM.list[[i]], init_pop = f.init,
                         harvest_fun = def_harv, disturb_fun = def_disturb)
    # Add it to the list
    eval(parse(text=paste0("species.list$", species.names[i], " <- species.i")))
  }
  
  # Return the output list
  return(species.list)
}



#' Generate a list of species object from the list of IPM
#' @param IPM.list_multiclim list of fitted IPM with different climate
#' @param f.init Function to initialize basal area
generate_species.list_multiclim <- function(IPM.list_multiclim, f.init = def_initBA(20)){
  
  # Names of the species
  species.names <- names(IPM.list_multiclim)
  
  # Initialize the list of species
  species.list <- list()
  
  # Loop on all IPM (and thus on all species)
  for(i in 1:length(species.names)){
    # Generate species object for each climate of species i
    species.cold.i <- species(IPM.list_multiclim[[i]]$cold, init_pop = f.init,
                              harvest_fun = def_harv, disturb_fun = def_disturb)
    species.optimum.i <- species(IPM.list_multiclim[[i]]$optimum, init_pop = f.init,
                                 harvest_fun = def_harv, disturb_fun = def_disturb)
    species.hot.i <- species(IPM.list_multiclim[[i]]$hot, init_pop = f.init,
                             harvest_fun = def_harv, disturb_fun = def_disturb)
    
    # Add it to the list
    eval(parse(text=paste0("species.list$", species.names[i], "$cold <- species.cold.i")))
    eval(parse(text=paste0("species.list$", species.names[i], "$optimum <- species.optimum.i")))
    eval(parse(text=paste0("species.list$", species.names[i], "$hot <- species.hot.i")))
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



#' Function to run simulations for a list of forests
#' @param forest.list_multiclim List of forest or species objects with multiple climates
#' @param tlim Number of simulation iterations (in years)
#' @param targetBA basal area targetted by the harvest module (in m2)
#' @param SurfEch Sampled area in ha
run_simulations_from_list_multiclim <- function(
  forest.list, tlim = 2000, targetBA = 40, SurfEch = 0.03){
  
  # Initialize output list
  list.out <- list()
  
  # Loop on all forests in the list
  for(i in 1:length(names(forest.list))){
    
    # Species printer
    print(paste0("Running simulations for species ", i, "/", length(names(forest.list)), 
                 " : ", gsub("\\.", "\\ and\\ ", names(forest.list)[i])))
    
    # Run simulations for cold edge
    print("--- cold edge")
    sim.cold.i <- sim_deter_forest(
      forest.list[[i]]$cold, tlim = tlim, equil_time = tlim, equil_dist = 1,
      harvest = "default", targetBA = targetBA, SurfEch = SurfEch, verbose = FALSE
    )
    eval(parse(text=paste0("list.out$", names(forest.list)[i], "$cold <- sim.cold.i")))
    
    # Run simulations for optimum
    print("--- optimum")
    sim.optimum.i <- sim_deter_forest(
      forest.list[[i]]$optimum, tlim = tlim, equil_time = tlim, equil_dist = 1,
      harvest = "default", targetBA = targetBA, SurfEch = SurfEch, verbose = FALSE
    )
    eval(parse(text=paste0("list.out$", names(forest.list)[i], "$optimum <- sim.optimum.i")))
    
    # Run simulations for hot edge
    print("--- hot edge")
    sim.hot.i <- sim_deter_forest(
      forest.list[[i]]$hot, tlim = tlim, equil_time = tlim, equil_dist = 1,
      harvest = "default", targetBA = targetBA, SurfEch = SurfEch, verbose = FALSE
    )
    eval(parse(text=paste0("list.out$", names(forest.list)[i], "$hot <- sim.hot.i")))
    
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



#' Generat a list oof forest, based on all combinations of one, two or three 
#' species from a list of species objects
#' @param species.list List of species objects
#' @param maxsp numeric: up to how many species per stand should we combine ?
#' @param harv_rules.ref List of rules for harvesting
generate_forest_list_multiclim <- function(species.list, maxsp = 1, harv_rules.ref){
  
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
  if(maxsp == 1) sp.combinations <- onesp
  if(maxsp == 2) sp.combinations <- c(onesp, twosp)
  if(maxsp == 3) sp.combinations <- c(onesp, twosp, threesp)
  
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
    
    # Loop on all climate
    for(j in 1:length(names(species.list[[1]]))){
      
      # Initialize list of species to include in forest i, climate j
      species.list.ij <- sapply(species.i, function(x) NULL)
      
      # Include species objects with the right climate 
      for(k in species.i) species.list.ij[[k]] <- species.list[[k]][[j]]
      
      # Generate forest i j
      forest.ij <- new_forest(species = species.list.ij, harv_rules = harv_rules.ref)
      
      # Add forest ij to the output list
      eval(parse(text=paste0("list.out$", sp.combinations.final[i], "$", 
                             names(species.list[[1]])[j]," <- forest.ij")))
      
    }
    
  }
  
  # Return the output list
  return(list.out)
}





#' Function to simulate disturbances from forest list and equilibrium
#' @param sim.list List of simulations to get population at equilibrium
#' @param forest.list List of forest used to get sim.list
#' @param disturbance.df disturbance dataframe
disturb_forest.list <- function(sim.list, forest.list, disturbance.df){
  # Initiate the output list
  list.out <- sim.list
  
  # Timing for the simulations with disturbance
  time.sim  <- max(tree_format(sim.list[[1]])$time, na.rm = TRUE) + max(disturbance.df$t)
  
  # Loop on all simulations
  for(i in 1:length(names(sim.list))){
    
    # Printer
    print(paste0("Running simulation ", i, "/", length(names(forest.list)), 
                 " - : ", gsub("\\.", "\\ and\\ ", names(forest.list)[i])))
    
    # Identify the species present in simulation i
    species.in.i <- unlist(strsplit(names(sim.list)[i], "\\."))
    
    # Initiate population 
    forest.in.i <- forest.list[[i]]
    
    # Formatted output of the simulation i
    memor.i <- tree_format(sim.list[[i]])
    
    # Loop on all species to update the forest
    for(j in 1:length(species.in.i)){
      
      # Get equilibrium for species i
      equil.j <- memor.i %>%
        filter(var == "m", equil, species == species.in.i[j]) %>% 
        pull(value)
      
      # Initiate the population at equilibrium
      forest.in.i$species[[j]]$init_pop <- def_init_k(equil.j*0.03)
      
      # Update disturbance function
      forest.in.i$species[[j]]$disturb_fun <- disturb_fun
      
      # Add disturbance coefficients
      forest.in.i$species[[j]]$disturb_coef <- filter(treeforce::disturb_coef, 
                                                      species == species.in.i[j])
      
    }
    
    
    # Simulate disturbance and add to the output list
    list.out[[i]] <- sim_deter_forest(forest.in.i, tlim = time.sim,
                                      equil_dist = time.sim, equil_time = time.sim,
                                      disturbance  = disturbance.df, verbose = FALSE) 
    
  }
  
  # Return the output list
  return(list.out)
  
}


#' Function to simulate disturbances from forest list and equilibrium
#' @param sim.list_multiclim List of simulations to get population at equilibrium
#' @param forest.list_multiclim List of forest used to get sim.list
#' @param disturbance.df disturbance dataframe
disturb_forest.list_multiclim <- function(sim.list, forest.list, disturbance.df){
  
  # Initiate the output list
  list.out <- sim.list
  
  # Timing for the simulations with disturbance
  time.sim  <- max(tree_format(sim.list[[1]][[1]])$time, na.rm = TRUE) + max(disturbance.df$t)
  
  # Loop on all simulations
  for(i in 1:length(names(sim.list))){
    
    # Printer
    print(paste0("Running simulation for forest ", i, "/", length(names(forest.list)), 
                 " : ", gsub("\\.", "\\ and\\ ", names(forest.list)[i])))
    
    # Identify the species present in simulation i
    species.in.i <- unlist(strsplit(names(sim.list)[i], "\\."))
    
    # Loop on all climates
    for(k in 1:length(names(sim.list[[i]]))){
      
      # Initiate population 
      forest.in.ik <- forest.list[[i]][[k]]
      
      # Formatted output of the simulation i
      memor.ik <- tree_format(sim.list[[i]][[k]])
      
      # Loop on all species to update the forest
      for(j in 1:length(species.in.i)){
        
        # Get equilibrium for species i
        equil.ijk <- memor.ik %>%
          filter(var == "m", equil, species == species.in.i[j]) %>% 
          pull(value)
        
        # Initiate the population at equilibrium
        forest.in.ik$species[[j]]$init_pop <- def_init_k(equil.ijk*0.03)
        
        # Update disturbance function
        forest.in.ik$species[[j]]$disturb_fun <- disturb_fun
        
        # Add disturbance coefficients
        forest.in.ik$species[[j]]$disturb_coef <- filter(treeforce::disturb_coef, 
                                                        species == species.in.i[j])
        
      }
      
      
      # Simulate disturbance and add to the output list
      list.out[[i]][[k]] <- sim_deter_forest(forest.in.ik, tlim = time.sim,
                                        equil_dist = time.sim, equil_time = time.sim,
                                        disturbance  = disturbance.df, verbose = FALSE) 
      
    }

  }
  
  # Return the output list
  return(list.out)
  
}

#' Get resilience and FD from simulations
#' @param sim.list.disturbed List of simulations where a disturbance occured
#' @param pc1_per_species Position of each species along the growth-mortality trade-off
get_resilience_and_FD <- function(sim.list.disturbed, pc1_per_species){
  
  # Initialize the output
  out <- data.frame(
    ID = c(1:length(names(sim.list.disturbed))),
    sp.combination = names(sim.list.disturbed), 
    FD = NA_real_, CWM = NA_real_, Resilience = NA_real_
  )
  
  
  # Loop on all species combination
  for(i in 1:dim(out)[1]){
    
    # Format the output
    data.i <- tree_format(sim.list.disturbed[[i]]) %>%
      filter(var == "BAsp") %>%
      filter(!equil)
    
    # Calculate resilience
    out$Resilience[i] <- 1/sum((data.i %>%
                                  group_by(time) %>%
                                  summarize(BA = sum(value)) %>%
                                  mutate(BA0 = .[which(.$time == 1), "BA"]) %>%
                                  mutate(diff = abs(BA - BA0)))$diff)
    
    # Calculate functional diversity and CWM at equilibrium
    FD_CWM_i <- data.i %>%
      filter(time == 1) %>%
      left_join(pc1_per_species, by = "species") %>%
      summarise(CWM = weighted.mean(pca1, w = value), 
                FD = weighted.var(pca1, w = value))
    out$CWM[i] <- FD_CWM_i$CWM
    out$FD[i] <- FD_CWM_i$FD
  }
  
  # Replace NA by 0
  out <- out %>% mutate(FD = ifelse(is.na(FD), 0, FD))
  
  # Return output
  return(out)
}




#' Disturbance function
#'
#' @param x population state distribution at time t
#' @param species The species class object of interest to get mesh and RDIcoef
#' values from. RDIcoef is a one line dataframe with RDI coefficient for one
#' species.
#' @param disturb Disturbance parameters. Highly depend on the disturbance
#' impact parameters given to the species.
#' @param ... Not used in this case.
#' \describe{
#' \item{qmd}{Forest Quadratic Mean Diameter}
#' }
#' @author Maxime Jeaunatre
#'
disturb_fun <- function(x, species, disturb = NULL, ...){
  dots <- list(...)
  qmd <- dots$qmd
  size <- species$IPM$mesh
  coef <- species$disturb_coef
  if(any(disturb$type %in% coef$disturbance)){
    coef <- subset(coef, disturbance == disturb$type)
  } else {
    stop(sprintf("The species %s miss this disturbance type (%s) parameters",
                 sp_name(species), disturb$type))
  }
  logratio <- log(size / qmd)
  dbh.scaled = coef$dbh.intercept + size * coef$dbh.slope
  logratio.scaled = coef$logratio.intercept + logratio * coef$logratio.slope
  Pkill <- plogis(coef$a0 + coef$a1 * logratio.scaled +
                    coef$b * disturb$intensity^(coef$c * dbh.scaled))
  return(x* Pkill) # always return the mortality distribution
}
