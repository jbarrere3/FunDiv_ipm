#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#### SCRIPT INTRODUCTION ####
#
#' @name functions_data.R  
#' @description R script containing all functions relative to data
#               importation and formatting
#' @author Julien Barrere
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


#' Function to create a list of IPM to run
#' @param climate list of climate objects
make_species_list = function(climate){
  
  # Loop on all climates
  for(i in 1:length(names(climate))){
    # Dataframe for climate i
    out.i = data.frame(
      ID.climate = i, 
      species = climate[[i]]$species
    ) %>%
      mutate(
        file = paste0("rds/climate_", ID.climate, "/species/", species, ".rds"))
    
    # Add to final dataset
    if(i == 1) out = out.i
    else out = rbind(out, out.i)
    
  }
  
  # Final formatting
  out = out %>%
    mutate(ID.species = c(1:dim(.)[1])) %>%
    dplyr::select(ID.species, ID.climate, species, file)
  
  # Return output
  return(out)
  
}


#' Function to make a species object, save it as rds and return filename
#' @param fit.list.allspecies demographic parameter for all species
#' @param climate.in object created by the function make climate
#' @param species_list table containing all species to create
#' @param ID.species.in ID of the species to make in species_list
make_species_rds = function(
  fit.list.allspecies, climate, species_list, ID.species.in){
  
  # Make IPM
  IPM.in = make_IPM(
    species = species_list$species[ID.species.in], 
    climate = climate[[species_list$ID.climate[ID.species.in]]]$climate, 
    fit =  fit.list.allspecies[[species_list$species[ID.species.in]]],
    clim_lab = paste0("climate_", species_list$ID.climate[ID.species.in]),
    mesh = c(m = 700, L = 100, U = as.numeric(
      fit.list.allspecies[[species_list$species[ID.species.in]]]$info[["max_dbh"]]) * 1.1),
    BA = 0:200, verbose = TRUE
  )
  
  # Create species object 
  species.in = species(
    IPM.in, init_pop = def_initBA(20), harvest_fun = def_harv, disturb_fun = def_disturb)
  
  # Save species object in a rdata
  create_dir_if_needed(species_list$file[ID.species.in])
  saveRDS(species.in, species_list$file[ID.species.in])
  
  # Return output list
  return(species_list$file[ID.species.in])
  
}


#' Function to create a list of forest for the simulations
#' @param climate list of climate objects
make_forest_list = function(climate){
  
  # Loop on all climates
  for(i in 1:length(names(climate))){
    # Dataframe for climate i
    out.i = data.frame(
      ID.climate = i, 
      combination = climate[[i]]$combinations
    ) %>%
      mutate(
        file.sim.equil = paste0(
          "rds/climate_", ID.climate, "/sim_equilibrium/", combination, ".rds"), 
        file.sim.dist = paste0(
          "rds/climate_", ID.climate, "/sim_disturbance/", combination, ".rds"))
    
    # Add to final dataset
    if(i == 1) out = out.i
    else out = rbind(out, out.i)
    
  }
  
  # Final formatting
  out = out %>%
    mutate(ID.forest = c(1:dim(.)[1])) %>%
    dplyr::select(ID.forest, ID.climate, combination, file.sim.equil, file.sim.dist)
  
  # Return output
  return(out)
  
}


#' Function to get for each simulation the equilibrium 
#' @param sim.equilibrium list of simulations until the equilibrium
get_equil = function(sim.equilibrium){
  
  # Initialize output
  list.out <- list()
  
  # First element of output: boolean to indicate if equilibrium is reached
  list.out$reached_equil = ifelse(
    is.na(sum((tree_format(sim.equilibrium) %>%
                 filter(var == "BAsp") %>%
                 filter(time == max(.$time) - 1))$value)), 
    FALSE, TRUE
  )
  
  # Second element of output: equilibrium state
  list.out$equil = tree_format(sim.equilibrium) %>%
    filter(var == "m", equil)
  
  # return output
  return(list.out)
}



#' Function to simulate disturbances from forest list and equilibrium
#' @param sim.list List of simulations to get population at equilibrium
#' @param forest.list List of forest used to get sim.list
#' @param disturbance.df disturbance dataframe
make_simulation_disturbance <- function(forest, equil, disturbance.df, time.sim = 3000){
  
  # Verify that equilibrium is reached
  if(equil$reached_equil){
    
    # Initialize the forest that will be used for the simulation
    forest.in = forest
    
    # Loop on all species present in the forest
    for(i in 1:length(names(forest$species))){
      
      # Name of species i
      species.i = names(forest$species)[i]
      
      # Get equilibrium for species i
      equil.i <- equil$equil %>%
        filter(species == species.i) %>% 
        pull(value)
      
      # Initiate the population at equilibrium
      forest.in$species[[i]]$init_pop <- def_init_k(equil.i*0.03)
      
      # Update disturbance function
      forest.in$species[[i]]$disturb_fun <- disturb_fun
      
      # Add disturbance coefficients
      forest.in$species[[i]]$disturb_coef <- filter(matreex::disturb_coef, 
                                                    species == species.i)
      
    }
    
    # Run simulation
    out = sim_deter_forest(forest.in, tlim = time.sim, equil_time = time.sim,
                           disturbance  = disturbance.df, verbose = FALSE) 
  } else { 
    # If equilibrium is not reached, return an empty matrix
    out = matrix()
  }
  
  
  # Return the output list
  return(out)
  
}



#' Get resilience, resistance and recovery from simulations with disturbance
#' @param sim.list.disturbed List of simulations where a disturbance occurred
#' @param disturbance.df disturbance dataset used to generate the disturbance
#' @param forest_list Table giving the information on each forest generated
get_resilience_metrics <- function(sim.list.disturbed, disturbance.df, 
                                   forest_list){
  
  # Initialize the output
  out <- forest_list %>%
    mutate(resistance = NA_real_, recovery = NA_real_, resilience = NA_real_, 
           t0 = NA_real_, thalf = NA_real_)
  
  # Loop on all species combination
  for(i in 1:dim(out)[1]){
    
    # First, verify that equilibrium was reached
    if(!is.na(sim.list.disturbed[[i]][1, 1])){
      
      # Format the output
      data.i <- tree_format(sim.list.disturbed[[i]]) %>%
        filter(var == "BAsp") %>%
        filter(!equil) %>%
        group_by(time) %>%
        summarize(BA = sum(value))
      
      
      ## Calculate resistance
      #  - Basal area at equilibrium
      Beq.i = mean((data.i %>% filter(time < min(disturbance.df$t)))$BA)
      # - Basal area after disturbance
      Bdist.i = (data.i %>% filter(time == max(disturbance.df$t)+1))$BA
      # - Resistance
      out$resistance[i] = Beq.i/(Beq.i - Bdist.i)
      
      ## Calculate recovery
      #  - Time at which population recovered fully
      Rec.time.i = min((data.i %>% 
                          filter(time > max(disturbance.df$t)) %>%
                          filter(BA > Beq.i))$time)
      # - Basal area 20 years after disturbance
      Bdist20.i = (data.i %>% filter(time == max(disturbance.df$t)+21))$BA
      # - Recovery = slope of BA increase in teh 20 years after disturbance
      out$recovery[i] = abs(Bdist20.i - Bdist.i)/20
      
      ## Calculate resilience
      out$resilience[i] <- 1/sum((data.i %>%
                                    mutate(BA0 = .[which(.$time == 1), "BA"]) %>%
                                    mutate(diff = abs(BA - BA0)))$diff)
      
      ## Calculate t0
      #  - Time at which population recovered to 5% of the basal area lost
      Rec.0.time.i = min((data.i %>% 
                            filter(time > max(disturbance.df$t)) %>%
                            filter(BA > (Beq.i + 19*Bdist.i)/20))$time)
      # - Recovery = time to recover minus time of disturbance
      out$t0[i] = Rec.0.time.i - max(disturbance.df$t)
      
      ## Calculate thalf
      #  - Time at which population recovered to 50% of the basal area lost
      Rec.half.time.i = min((data.i %>% 
                               filter(time > max(disturbance.df$t)) %>%
                               filter(BA > (Beq.i + Bdist.i)/2))$time)
      # - Recovery = time to recover minus time of disturbance
      out$thalf[i] = Rec.half.time.i - max(disturbance.df$t)
      
    }
    
  }
  
  # Return output
  return(out)
}

#' Get functional diversity from simulations
#' @param ID.forest Vector containing the ID of each forest simulated
#' @param sim.list.disturbed List of simulations where a disturbance occured
#' @param pc1_per_species Position of each species along the growth-mortality trade-off
get_FD <- function(ID.forest, sim.list.disturbed, pc1_per_species){
  
  # Initialize the output
  out <- data.frame(
    ID.forest = ID.forest,
    FD = NA_real_, CWM = NA_real_)
  
  
  # Loop on all species combination
  for(i in 1:dim(out)[1]){
    
    # First, verify that equilibrium was reached
    if(!is.na(sim.list.disturbed[[i]][1, 1])){
      # Format the output
      data.i <- tree_format(sim.list.disturbed[[i]]) %>%
        filter(var == "BAsp") %>%
        filter(time == 1) %>%
        left_join(pc1_per_species, by = "species") %>%
        summarise(CWM = weighted.mean(pca1, w = value), 
                  FD = weighted.var(pca1, w = value))
      # Add to the final dataframe
      out$CWM[i] <- data.i$CWM
      out$FD[i] <- data.i$FD
    }
    
  }
  
  # Replace NA by 0
  out <- out %>% mutate(FD = ifelse(is.na(FD), 0, FD))
  
  # Return output
  return(out)
}


#' Function to generate a list with climate and all possible sp combinations
#' @param FUNDIV_climate_species data with climate and sp presence per plot
#' @param quantiles.in range between 0 and 1 of pca1 value to select
#' @param disturbance.in name of the disturbance we plan to apply to filter 
#'                       species combinations compatible
#' @param nsp_per_richness number of sp combinations to select per sp richness
make_climate <- function(FUNDIV_climate_species, quantiles.in, 
                         disturbance.in = "storm", 
                         nsp_per_richness = 10){
  
  # Initialize output list
  out = list()
  
  # Get the range of pca_values based on quantiles.in
  climate.in = as.numeric(quantile(FUNDIV_climate_species$pca1, 
                                   probs = quantiles.in))
  
  # climate vector
  out$climate = setNames(
    object = as.numeric(
      FUNDIV_climate_species %>%
        filter(pca1 > climate.in[1] & pca1 < climate.in[2]) %>%
        summarize(sgdd = mean(sgdd), wai = mean(wai), PC1 = mean(pca1), PC2 = mean(pca2)) %>%
        mutate(sgdd2 = sgdd^2, wai2 = wai^2, sgddb = 1/sgdd, 
               waib = 1/(1 + wai), N = 2, SDM = 0) %>%
        dplyr::select(sgdd, wai, sgddb, waib, wai2, sgdd2, PC1, PC2, N, SDM)
    ), c("sgdd", "wai", "sgddb", "waib", "wai2", "sgdd2", "PC1", "PC2", "N", "SDM")
  )
  
  # Vector of all species
  species_vec = colnames(FUNDIV_climate_species)[grep("_", colnames(FUNDIV_climate_species))]
  # Remove hornbeam from this vector as too unstable in the simulations
  species_vec = species_vec[species_vec != "Carpinus_betulus"]
  
  # Adjust species combinations to the disturbance if one is specified
  if(disturbance.in %in% c("storm", "fire", "biotic")){
    # Vector of all species for which we have disturbance parameters
    data("disturb_coef")
    species_vec_dist = (disturb_coef %>%
                          filter(disturbance %in% disturbance.in))$species
    # Restrict the species vector to these species
    species_vec = species_vec[which(species_vec %in% species_vec_dist)]
  }
  
  # species combinations for this climate
  # -- Paste all species presence absence to create binary code
  eval(parse(text = paste0(
    "data.in <- FUNDIV_climate_species %>% mutate(combi = paste0(", 
    paste(species_vec, collapse = ", "), "))")))
  # -- Calculate the number of species per species combination
  eval(parse(text = paste0(
    "data.in <- data.in %>% mutate(n.sp = ", 
    paste(species_vec, collapse = " + "), ")")))
  # -- Restrict to the climate specified
  data.in <- data.in %>%
    filter(pca1 > climate.in[1] & pca1 < climate.in[2])
  # -- Select the main combinations
  data_codes <- data.in %>%
    group_by(combi, n.sp) %>%
    summarize(n = n()) %>%
    filter(n.sp > 0) %>%
    arrange(desc(n.sp), desc(n))
  codes = c()
  for(j in 1:length(unique(data_codes$n.sp))){
    codes.j = (data_codes %>%
                 filter(n.sp == unique(data_codes$n.sp)[j]))$combi
    if(length(codes.j) > nsp_per_richness) codes = c(codes, codes.j[c(1:nsp_per_richness)])
    else codes = c(codes, codes.j)
  }
  # -- initialize combinations and species vector
  combinations.in = c(); species.in = c()
  # -- Loop on all codes
  for(i in 1:length(codes)){
    # Decode to get a vector of species
    vec.i = decode_species(codes[i], species_vec)
    # Add combination to the vector
    combinations.in = c(combinations.in, paste(vec.i, collapse = "."))
    # Store all species in species vector
    species.in = c(species.in, vec.i)
  } 
  # Add to the final list the combinations and the list of all species
  out$combinations = combinations.in
  out$species = unique(species.in)
  
  # Return output
  return(out)
}


#' Function that creates climate list based on quantiles
#' @param n.clim integer: number of climate to create in the list
create_climate_list = function(n.clim){
  
  # Vector that contains a sequence of quantiles value from 0 to 1
  vec.in = seq(from = 0, to = 1, length.out = n.clim+1)
  
  # Initialize the output list
  list.out = vector(mode = "list", length = n.clim)
  
  # Loop on all climates
  for(i in 1:n.clim){
    # Attribute a name to climate i
    names(list.out)[i] = paste("quantile", vec.in[i], vec.in[i+1], sep = "_")
    # Add vector of quantile for climate i
    list.out[[i]] = c(vec.in[i], vec.in[i+1])
  }
  
  # Return final list
  return(list.out)
}






#' Function to convert a binary code in a vector of species
#' @param code binary code where one indicates presence, 0 absence
#' @param species_vec vector of species, same length as code
decode_species <- function(code, species_vec){
  (data.frame(present = as.numeric(strsplit(code, split = "")[[1]]), 
              species = species_vec) %>%
     filter(present == 1))$species
}





#' Function to format resilience and fd for different climates
#' @param list.in list with one element per climate, each element being a list 
#'                with two elements: resilience and FD
format_resilience_FD = function(list.in){
  # Loop on all climates
  for(i in 1:length(names(list.in))){
    
    # Format climate i as dataframe
    data.i = list.in[[i]]$resilience %>%
      left_join(list.in[[i]]$FD, by = c("ID", "sp.combination")) %>%
      mutate(climate = names(list.in)[i])
    
    # Add to final output
    if(i == 1) data = data.i
    else data = rbind(data, data.i)
  }
  
  # Return final data frame
  return(data)
}


#' Format data for the models
#' @param climate list of climate objects used for the simulations
#' @param resilience_metrics df with resilience per forest ID
#' @param FD df with FD and CWM per forest ID
get_data_model = function(climate, resilience_metrics, FD){
  
  # Build a dataset to associate ID climate with sgdd, wai and pca1
  data.climate = data.frame(
    ID.climate = c(1:length(names(climate))), 
    pca1 = NA_real_, 
    sgdd = NA_real_, 
    wai = NA_real_
  )
  for(i in 1:dim(data.climate)[1]){
    data.climate$pca1[i] = climate[[i]]$climate[7]
    data.climate$sgdd[i] = climate[[i]]$climate[1]
    data.climate$wai[i] = climate[[i]]$climate[2]
  }
  
  # Format final dataset
  out = resilience_metrics %>%
    left_join(FD, by = "ID.forest") %>%
    left_join(data.climate, by = "ID.climate") %>%
    dplyr::select(ID.forest, ID.climate, forest.composition = combinations, 
                  sgdd, wai, pca1, FD, CWM, resistance, recovery, resilience)
  
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
