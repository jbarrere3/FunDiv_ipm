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


#' Function to generate a list with climate with species combinations
#' @param FUNDIV_climate_species data with climate and sp presence per plot
#' @param quantiles.in range between 0 and 1 of pca1 value to select
#' @param disturbance.in name of the disturbance we plan to apply to filter 
#'                       species combinations compatible
#' @param nsp_per_richness number of sp combinations to select per sp richness
#' @param exclude.in vector of species to exclude if bad estimation or IPM fit
#' @param method way to select species combination: most frequent ("frequency") or "random" 
make_climate <- function(FUNDIV_climate_species, quantiles.in, disturbance.in, 
                         nsp_per_richness = 10, exclude.in = c("Carpinus_betulus"), 
                         method = "frequency"){
  
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
  # Remove species with bad estimation or unstable in simulations
  species_vec = species_vec[!(species_vec %in% exclude.in)]
  
  # Adjust species combinations to the disturbance if one is specified
  if(disturbance.in %in% c("storm", "fire", "biotic")){
    # Vector of all species for which we have disturbance parameters
    data("disturb_coef")
    species_vec_dist = (disturb_coef %>%
                          filter(disturbance %in% disturbance.in))$species
    # Restrict the species vector to these species
    species_vec = species_vec[which(species_vec %in% species_vec_dist)]
  }
  
  
  
  # species combinations for this climate based on data (frequency method)
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
  # -- Identify all species for this climate
  species.in = unique(species.in)
  
  
  
  # Make a random selection of species with same number of forest per richness
  # -- Identify the number of combinations per species richness
  combi.per.richness = data_codes %>%
    filter(combi %in% codes) %>%
    group_by(n.sp) %>%
    summarize(n = n()) %>%
    arrange(n.sp)
  # -- Initialize the vector of random combinations
  combinations.random.in = c()
  # -- Loop on all levels of richness
  for(r in 1:dim(combi.per.richness)[1]){
    # All existing combinations for richness r
    combinations.r.all = apply(as.data.frame(t(combn(species.in, r))), 1, 
                               paste, collapse = "." )
    # Randomly sample the same number of forest per richness as with the data approach
    combinations.r = sample(combinations.r.all, size = combi.per.richness$n[r], 
                            replace = FALSE)
    # Add to the vector of random species combinations
    combinations.random.in = c(combinations.random.in, combinations.r)
  }
  
  
  # Add to the final list the combinations and the list of all species
  if(method == "frequency") out$combinations = combinations.in
  if(method == "random") out$combinations = combinations.random.in
  out$species = species.in
  
  # Return output
  return(out)
}


#' Function that creates climate list based on quantiles
#' @param n.clim integer: number of climate to create in the list
#' @param quantile.range numeric vector of length two indicating the climatic range 
create_climate_list = function(n.clim, quantile.range = c(0, 1)){
  
  # Vector that contains a sequence of quantiles value from 0 to 1
  vec.in = seq(from = quantile.range[1], to = quantile.range[2], 
               length.out = n.clim+1)
  
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



#' Function to create a list of IPM to run
#' @param climate list of climate objects
#' @param disturbance.in character: name of the disturbance (for file archiving)
make_species_list = function(climate, disturbance.in){
  
  # Loop on all climates
  for(i in 1:length(names(climate))){
    # Dataframe for climate i
    out.i = data.frame(
      ID.climate = i, 
      species = climate[[i]]$species
    ) %>%
      mutate(
        file = paste0("rds/", disturbance.in, "/climate_", ID.climate, "/species/", species, ".rds"))
    
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
    BA = 0:200, verbose = TRUE, correction = "none"
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
#' @param disturbance.in character: name of the disturbance (for file archiving)
make_forest_list = function(climate, disturbance.in){
  
  # Loop on all climates
  for(i in 1:length(names(climate))){
    # Dataframe for climate i
    out.i = data.frame(
      ID.climate = i, 
      combination = climate[[i]]$combinations
    ) %>%
      mutate(
        file.sim.equil = paste0(
          "rds/", disturbance.in, "/climate_", ID.climate, "/sim_equilibrium/", combination, ".rds"), 
        file.sim.dist = paste0(
          "rds/", disturbance.in, "/climate_", ID.climate, "/sim_disturbance/", combination, ".rds"))
    
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

#' Function to make a list of simulations till equilibrium
#' @param climate object generated by make_climate
#' @param harv_rules.ref rules for harvesting (needed to generate forest)
#' @param species_list df with information on all species object
#' @param forest_list df with information on all forest to simulate 
#' @param species vector containing all species rds files created
#' @param ID.forest.in ID of the forest to simulate in forest_list
make_simulations_equilibrium = function(climate, harv_rules.ref, species_list, 
                                        forest_list, species, ID.forest.in){
  
  # Identify ID of the climate
  ID.climate.in = forest_list$ID.climate[ID.forest.in]
  
  # Identify species combination i
  combination.in = forest_list$combination[ID.forest.in]
  
  # vector of species in forest i
  species.in = unlist(strsplit(combination.in, "\\."))
  
  # Initialize list of species to create
  list.species <- vector("list", length(species.in))
  names(list.species) = species.in
  
  # Loop on all species
  for(i in 1:length(species.in)){
    
    # Identify the file in species containing species i
    species.file.i = species[(species_list %>%
                                filter(species == species.in[i]) %>%
                                filter(ID.climate == ID.climate.in))$ID.species]
    
    # Store the file in the list
    list.species[[i]] = readRDS(species.file.i)
    
  }
  
  
  # Make forest
  forest.in = new_forest(species = list.species, harv_rules = harv_rules.ref)
  
  # Run simulation till equilibrium
  sim.in = sim_deter_forest(
    forest.in, tlim = 4000, equil_time = 10000, equil_dist = 250, 
    equil_diff = 1, harvest = "default", SurfEch = 0.03, verbose = TRUE)
  
  # Save simulation in a rdata
  create_dir_if_needed(forest_list$file.sim.equil[ID.forest.in])
  saveRDS(sim.in, forest_list$file.sim.equil[ID.forest.in])
  
  # Return output list
  return(forest_list$file.sim.equil[ID.forest.in])
}




#' Function to make a list of simulations with disturbance
#' @param climate object generated by make_climate
#' @param harv_rules.ref rules for harvesting (needed to generate forest)
#' @param species_list df with information on all species object
#' @param forest_list df with information on all forest to simulate 
#' @param species vector containing all species rds files created
#' @param sim_equilibrium Vector containing file names of simulations till equil
#' @param ID.forest.in ID of the forest to simulate in forest_list
#' @param disturbance.df disturbance dataframe
make_simulations_disturbance = function(
  climate, harv_rules.ref, species_list, forest_list, species, sim_equilibrium, 
  ID.forest.in, disturbance.df){
  
  # Identify ID of the climate
  ID.climate.in = forest_list$ID.climate[ID.forest.in]
  
  # Identify species combination i
  combination.in = forest_list$combination[ID.forest.in]
  
  # vector of species in forest i
  species.in = unlist(strsplit(combination.in, "\\."))
  
  # Initialize list of species to create
  list.species <- vector("list", length(species.in))
  names(list.species) = species.in
  
  # Read the simulation at equilibrium
  sim_equilibrium.in = readRDS(sim_equilibrium[ID.forest.in])
  
  # Checked that the population reached equilibrium
  reached_equil = ifelse(
    is.na(sum((sim_equilibrium.in %>%
                 filter(var == "BAsp") %>%
                 filter(time == max(.$time) - 1))$value)), 
    FALSE, TRUE
  )
  
  # Only make the simulation id ==f population reached an equilibrium
  if(reached_equil){
    # Loop on all species
    for(i in 1:length(species.in)){
      
      # Identify the file in species containing species i
      species.file.i = species[(species_list %>%
                                  filter(species == species.in[i]) %>%
                                  filter(ID.climate == ID.climate.in))$ID.species]
      
      # Store the file in the list
      list.species[[i]] = readRDS(species.file.i)
      
      # Extract the equilibrium for species i
      equil.i = sim_equilibrium.in %>%
        filter(var == "n", equil, species == species.in[i]) %>% 
        pull(value)
      
      # Initiate the population at equilibrium
      list.species[[i]]$init_pop <- def_init_k(equil.i*0.03)
      
      # Update disturbance function
      list.species[[i]]$disturb_fun <- disturb_fun
      
      # Add disturbance coefficients
      list.species[[i]]$disturb_coef <- filter(matreex::disturb_coef, 
                                               species == species.in[i])
    }
    
    
    # Make forest
    forest.in = new_forest(species = list.species, harv_rules = harv_rules.ref)
    
    # Run simulation till equilibrium
    sim.in = sim_deter_forest(
      forest.in, tlim = 4000, equil_time = 4000, disturbance = disturbance.df, 
      SurfEch = 0.03, verbose = TRUE)
  } else {
    sim.in = matrix()
  }
  
  
  # Save simulation in a rdata
  create_dir_if_needed(forest_list$file.sim.dist[ID.forest.in])
  saveRDS(sim.in, forest_list$file.sim.dist[ID.forest.in])
  
  # Return output list
  return(forest_list$file.sim.dist[ID.forest.in])
}




#' Get resilience, resistance and recovery from simulations with disturbance
#' @param sim_disturbance vector containing the file names of sim with disturbance
#' @param disturbance.df disturbance dataset used to generate the disturbance
#' @param forest_list Table giving the information on each forest generated
get_resilience_metrics <- function(sim_disturbance, disturbance.df, 
                                   forest_list){
  
  # Initialize the output
  out <- forest_list %>%
    dplyr::select(ID.forest, ID.climate, combination) %>%
    mutate(resistance = NA_real_, recovery = NA_real_, resilience = NA_real_, 
           t0 = NA_real_, thalf = NA_real_)
  
  # Loop on all species combination
  for(i in 1:length(sim_disturbance)){
    
    # Printer
    print(paste0("Reading simulation ", i, "/", length(sim_disturbance)))
    
    # Read simulation i
    sim.i = readRDS(sim_disturbance[i])
    
    # First, verify that equilibrium was reached
    if(!is.na(sim.i[1, 1])){
      
      # Format the output
      data.i <- sim.i %>%
        filter(var == "BAsp") %>%
        filter(!equil) %>%
        group_by(time) %>%
        summarize(BA = sum(value))
      
      
      ## Calculate resistance
      #  - Basal area at equilibrium
      Beq.i = mean((data.i %>% filter(time < min(disturbance.df$t)))$BA)
      # - Basal area after disturbance
      Bdist.i = (data.i %>% filter(time == max(disturbance.df$t)+1))$BA
      # - Resistance : logit of the percentage of basal area that survived 
      #out$resistance[i] = Beq.i/(Beq.i - Bdist.i)
      out$resistance[i] = log((Bdist.i/Beq.i)/(1 - (Bdist.i/Beq.i)))
      
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
#' @param forest_list df containing info on each forest simulated
#' @param sim_disturbance vector containing the file names of sim with disturbance
#' @param pc1_per_species Position of each species along the growth-mortality trade-off
get_FD_original <- function(forest_list, sim_disturbance, pc1_per_species){
  
  # Initialize the output
  out <- forest_list %>%
    dplyr::select(ID.forest, ID.climate, combination) %>%
    mutate(nsp = NA_real_, FD = NA_real_)
  
  
  # Loop on all species combination
  for(i in 1:length(sim_disturbance)){
    
    # Fill the number of species
    out$nsp[i] = length(unlist(strsplit(out$combination[i], "\\.")))
    
    # Read simulation i
    sim.i = readRDS(sim_disturbance[i])
    
    # First, verify that equilibrium was reached
    if(!is.na(sim.i[1, 1])){
      
      # Format the output
      data.i <- sim.i %>%
        filter(var == "BAsp") %>%
        filter(time == 1) %>%
        left_join(pc1_per_species, by = "species") %>%
        summarise(FD = weighted.var(pca1, w = value))
      # Add to the final dataframe
      out$FD[i] <- data.i$FD
    }
    
  }
  
  # Replace NA by 0
  out <- out %>% mutate(FD = ifelse(is.na(FD), 0, FD))
  
  # Return output
  return(out)
}

#' Get functional diversity from simulations with FD package
#' @param forest_list df containing info on each forest simulated
#' @param sim_disturbance vector containing the file names of sim with disturbance
#' @param pc1_per_species Position of each species along the growth-mortality trade-off
get_FD_dimension <- function(forest_list, sim_disturbance, pc1_per_species){
  
  # Initialize vector of all successful simulations
  vec.sim = c()
  
  # Loop on all species combination
  for(i in 1:length(sim_disturbance)){
    
    # Read simulation i
    sim.i = readRDS(sim_disturbance[i])
    
    # First, verify that equilibrium was reached
    if(!is.na(sim.i[1, 1])){
      
      # Format the output
      data.i <- sim.i %>%
        mutate(ID.climate = forest_list$ID.climate[i], 
               ID.forest = forest_list$ID.forest[i], 
               ID.community = paste(ID.climate, ID.forest, sep = ".")) %>%
        filter(var == "BAsp") %>%
        filter(time == 1) %>%
        dplyr::select(ID.climate, ID.forest, ID.community, species, value)
      
      # Also check we have trait value for all species in the community
      if(all(data.i$species %in% pc1_per_species$species)){
        # Add to the final dataframe
        if(length(vec.sim) == 0) data = data.i
        else data = rbind(data, data.i)
        
        # Increment the counter
        vec.sim = c(vec.sim, i)
      }
    }
    
  }
  
  # Abundance dataframe
  abun.df = data %>%
    spread(key = "species", value = "value") %>% 
    replace(is.na(.), 0)
  
  # Abundance matrix
  abun.matrix = as.matrix(
    abun.df %>% dplyr::select(-ID.climate, -ID.forest, -ID.community))
  rownames(abun.matrix) = abun.df$ID.community
  
  # Trait df
  trait.df = data.frame(species = colnames(abun.matrix)) %>%
    left_join(pc1_per_species, by = "species") %>%
    dplyr::select(-species)
  rownames(trait.df) = colnames(abun.matrix)
  
  # Calculate FD per community
  fd.raw <- dbFD(trait.df, abun.matrix, w.abun = TRUE)
  
  # Final dataframe
  out <- data.frame(
    ID.forest = forest_list$ID.forest[vec.sim], 
    ID.climate = forest_list$ID.climate[vec.sim], 
    combination = forest_list$combination[vec.sim], 
    FRic = fd.raw$FRic, 
    FDis = fd.raw$FDis, 
    CWM = fd.raw$CWM$pca1
  )
  
  # Return output
  return(out)
}

#' Get functional diversity from simulations 
#' @param forest_list df containing info on each forest simulated
#' @param sim_disturbance vector containing the file names of sim with disturbance
#' @param pc1_per_species Position of each species along the growth-mortality trade-off
get_FD <- function(forest_list, sim_disturbance, pc1_per_species){
  
  # Get the original FD index
  FD_original = get_FD_original(forest_list, sim_disturbance, pc1_per_species)
  
  # Get the original FD index
  FD_dimension = get_FD_dimension(forest_list, sim_disturbance, pc1_per_species)
  
  # Join the two datasets
  out = left_join(FD_original, FD_dimension, 
                  by = c("ID.forest", "ID.climate", "combination")) %>%
    filter(!is.na(FRic))
  
  # Return output
  return(out)
}



#' Format data for the models
#' @param climate list of climate objects used for the simulations
#' @param resilience df with resilience per forest ID
#' @param FD df with FD and CWM per forest ID
get_data_model = function(climate, resilience, FD){
  
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
  out = resilience %>%
    left_join(FD, by = c("ID.forest", "ID.climate", "combination")) %>%
    left_join(data.climate, by = "ID.climate") %>%
    rename(forest.composition = combination)
  
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
  
  # edits for delay
  size[size == 0] <- min(size[size !=0])
  
  logratio <-  log(size / qmd)
  dbh.scaled = coef$dbh.intercept + size * coef$dbh.slope
  logratio.scaled = coef$logratio.intercept + logratio * coef$logratio.slope
  Pkill <- plogis(coef$a0 + coef$a1 * logratio.scaled + 
                    coef$b * disturb$intensity ^(coef$c * dbh.scaled))
  
  return(x* Pkill) # always return the mortality distribution
}
