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
                    pca1 = get_pca_ind(pca)[[1]][, 1]) %>%
    mutate(species = ifelse(species == "Betula_pendula", "Betula", species))
  
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
#' @param disturb_coef.in disturbance coefficients
make_climate <- function(FUNDIV_climate_species, quantiles.in, disturbance.in, 
                         nsp_per_richness = 10, exclude.in = c("Carpinus_betulus"), 
                         method = "frequency", disturb_coef.in, pc1_per_species){
  
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
  # Remove species for which we have no trait estimation
  species_vec = species_vec[(species_vec %in% pc1_per_species$species)]
  
  # Adjust species combinations to the disturbance if one is specified
  if(disturbance.in %in% c("storm", "fire", "biotic")){
    # Vector of all species for which we have disturbance parameters
    species_vec_dist = (disturb_coef.in %>%
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
    forest.in, tlim = 4000, equil_time = 50000, equil_dist = 2000, 
    equil_diff = 0.5, harvest = "default", SurfEch = 0.03, verbose = TRUE)
  
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
#' @param disturb_coef.in table containing the disturbance coefficients per species
#' @param species vector containing all species rds files created
#' @param sim_equilibrium Vector containing file names of simulations till equil
#' @param ID.forest.in ID of the forest to simulate in forest_list
#' @param disturbance.df disturbance dataframe
make_simulations_disturbance = function(
  climate, harv_rules.ref, species_list, forest_list, disturb_coef.in, species, 
  sim_equilibrium, ID.forest.in, disturbance.df){
  
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
      list.species[[i]]$disturb_coef <- filter(disturb_coef.in, 
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
           t0 = NA_real_, thalf = NA_real_, SD = NA_real_, BA_diff = NA_real_, 
           BA_eq = NA_real_, dbh_mean = NA_real_, dbh_q10 = NA_real_, 
           dbh_q90 = NA_real_, dbh_mean_postdist = NA_real_, 
           dbh_q10_postdist = NA_real_, dbh_q90_postdist = NA_real_)
  
  # Identify disturbance time
  tdist = min(disturbance.df$t)
  
  # Loop on all species combination
  for(i in 1:length(sim_disturbance)){
    
    # Printer
    print(paste0("Reading simulation ", i, "/", length(sim_disturbance)))
    
    # Read simulation i
    sim.i = readRDS(sim_disturbance[i])
    
    # First, verify that equilibrium was reached
    if(!is.na(sim.i[1, 1])){
      
      # mean dbh at equilibrium and after disturbance
      dbh_i = sim.i %>%
        filter(var == "n") %>%
        filter(time %in% c(1, (max(disturbance.df$t)+1))) %>%
        group_by(size, time) %>%
        summarize(ntot = sum(value)) %>%
        ungroup() %>% group_by(time) %>%
        filter(size > 0) %>%
        mutate(ntot_size = ntot*size) %>%
        summarize(mean_dbh = weighted.mean(size, w = ntot), 
                  q10_dbh = weighted.quantile(size, w = ntot, prob = 0.1), 
                  q90_dbh = weighted.quantile(size, w = ntot, prob = 0.9))
      out$dbh_mean[i] <- subset(dbh_i, time == 1)$mean_dbh
      out$dbh_q10[i] <- subset(dbh_i, time == 1)$q10_dbh
      out$dbh_q90[i] <- subset(dbh_i, time == 1)$q90_dbh
      out$dbh_mean_postdist[i] <- subset(dbh_i, time != 1)$mean_dbh
      out$dbh_q10_postdist[i] <- subset(dbh_i, time != 1)$q10_dbh
      out$dbh_q90_postdist[i] <- subset(dbh_i, time != 1)$q90_dbh
      
      # Format the output
      data.i <- sim.i %>%
        filter(var == "BAsp") %>%
        filter(!equil) %>%
        group_by(time) %>%
        summarize(BA = sum(value))
      
      ## Calculate stability before disturbance (to check equilibrium)
      out$SD[i] = sd(subset(data.i, time %in% c(1:(tdist-1)))$BA)
      out$BA_diff[i] = diff(range(subset(data.i, time %in% c(1:(tdist-1)))$BA))
      
      ## Calculate resistance
      #  - Basal area at equilibrium
      Beq.i = mean((data.i %>% filter(time < min(disturbance.df$t)))$BA)
      out$BA_eq[i] = Beq.i
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



#' Get functional diversity from simulations with FD package
#' @param forest_list df containing info on each forest simulated
#' @param sim_disturbance vector containing the file names of sim with disturbance
#' @param pc1_per_species Position of each species along the growth-mortality trade-off
get_FD <- function(forest_list, sim_disturbance, pc1_per_species){
  
  # Initialize vector of all successful simulations
  vec.sim = c()
  
  # Initialize the original fd data
  data.fd.original <- forest_list %>%
    dplyr::select(ID.forest, ID.climate, combination) %>%
    mutate(nsp = NA_real_, FD = NA_real_, H = NA_real_, D = NA_real_, Nha = NA_real_)
  
  # Loop on all species combination
  for(i in 1:length(sim_disturbance)){
    
    # Printer
    print(paste0(i, "/", length(sim_disturbance)))
    
    # Fill the number of species
    data.fd.original$nsp[i] = length(unlist(strsplit(data.fd.original$combination[i], "\\.")))
    
    # Read simulation i
    sim.i = readRDS(sim_disturbance[i])
    
    # First, verify that equilibrium was reached
    if(!is.na(sim.i[1, 1])){
      
      ## - Calculate the number of trees per ha at equilibrium
      data.fd.original$Nha[i] = sum((sim.i %>%
                                       filter(var == "N") %>%
                                       filter(time == 1))$value)
      
      ## - FD with the original approach
      data.fd.original.i <- sim.i %>%
        filter(var == "BAsp") %>%
        filter(time == 1) %>%
        left_join(pc1_per_species, by = "species") %>%
        mutate(p = value/sum(.$value), 
               plnp = p*log(p), 
               p2 = p^2) %>%
        summarise(FD = weighted.var(pca1, w = value), 
                  H = -sum(plnp), 
                  D = 1/sum(p2)) 
      data.fd.original$FD[i] <- data.fd.original.i$FD
      data.fd.original$H[i] <- data.fd.original.i$H
      data.fd.original$D[i] <- data.fd.original.i$D
      
      ## - FD with the FD package (dimension approach)
      data.fd.dimension.i <- sim.i %>%
        mutate(ID.climate = forest_list$ID.climate[i], 
               ID.forest = forest_list$ID.forest[i], 
               ID.community = paste(ID.climate, ID.forest, sep = ".")) %>%
        filter(var == "BAsp") %>%
        filter(time == 1) %>%
        dplyr::select(ID.climate, ID.forest, ID.community, species, value)
      
      # Also check we have trait value for all species in the community
      if(all(data.fd.dimension.i$species %in% pc1_per_species$species)){
        # Add to the final dataframe
        if(length(vec.sim) == 0) data.fd.dimension = data.fd.dimension.i
        else data.fd.dimension = rbind(data.fd.dimension, data.fd.dimension.i)
        
        # Increment the counter
        vec.sim = c(vec.sim, i)
      }
    }
    
  }
  
  # Replace NA by 0 in original approach
  data.fd.original <- data.fd.original %>% mutate(FD = ifelse(is.na(FD), 0, FD))
  
  # Abundance dataframe
  abun.df = data.fd.dimension %>%
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
  data.fd.dimension.final <- data.frame(
    ID.forest = forest_list$ID.forest[vec.sim], 
    ID.climate = forest_list$ID.climate[vec.sim], 
    combination = forest_list$combination[vec.sim], 
    FRic = fd.raw$FRic, 
    FDis = fd.raw$FDis, 
    CWM = fd.raw$CWM$pca1
  )
  
  # Join the two datasets
  out = left_join(data.fd.original, data.fd.dimension.final, 
                  by = c("ID.forest", "ID.climate", "combination")) %>%
    filter(!is.na(FRic))
  
  # Return output
  return(out)}


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



#' Function to extract the list of species which are included in the simulations
#' @param climate.in object generated by the function make_climate
get_species_sim = function(climate.in){
  
  # Initialize output
  out = c()
  
  # Loop on all climates to extract the species present
  for(i in 1:length(names(climate.in))) out = c(out, climate.in[[i]]$species)
  
  # Remove duplicates
  out = unique(out)
  
  # Return vector generated
  return(out)
}


#' Function to get a new forest_list with only simulations that did not reach equilibrium
#' @param data_model_all formatted simulation output, with all simul 
#' @param forest_list data frame with all information on forest simulated
#' @param sim_equilibrium vector containing name of rds files with equilibrium simulations
get_forest_list_noneq = function(data_model_all, forest_list, sim_equilibrium){
  
  # Type of disturbance
  disturbance.in = gsub("\\/.+", "", gsub("rds\\/", "", forest_list$file.sim.equil[1]))
  
  # Format the output based on initial forest list
  out = forest_list %>%
    # Restrict to the simulations that did not reach equilibrium
    filter(ID.forest %in% subset(data_model_all, SD > 0.15)$ID.forest) %>%
    # Change the name of the files to save
    mutate(file.sim.dist = paste0("rds/", disturbance.in, "/noneq/dist_climate", 
                                  ID.climate, "_", combination, ".rds"))
  
  # Return output
  return(out)
}


#' Function to make a list of simulations with disturbance for pop not at equilibrium
#' @param climate object generated by make_climate
#' @param harv_rules.ref rules for harvesting (needed to generate forest)
#' @param species_list df with information on all species object
#' @param forest_list df with information on all forest to simulate 
#' @param species vector containing all species rds files created
#' @param sim_equilibrium Vector containing file names of simulations till equil
#' @param ID.forest.in ID of the forest to simulate in forest_list
#' @param disturbance.df disturbance dataframe
make_simulations_disturbance_noneq = function(
  climate, harv_rules.ref, species_list, forest_list_noneq, species, sim_equilibrium, 
  ID.forest.in, disturbance.df){
  
  # Identify ID of the climate
  ID.climate.in = subset(forest_list_noneq, ID.forest == ID.forest.in)$ID.climate
  
  # Identify species combination i
  combination.in = subset(forest_list_noneq, ID.forest == ID.forest.in)$combination
  
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
    
    
    # First loop on all species to build the forest to simulate
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
      
    }
    
    
    # Make forest
    forest.in = new_forest(species = list.species, harv_rules = harv_rules.ref)
    
    # Make a simulation from initial equilibrium to find a better one
    sim.equil = sim_deter_forest(
      forest.in, tlim = 2500, equil_time = 100000, equil_dist = 2000, 
      equil_diff = 0.5, harvest = "default", SurfEch = 0.03, verbose = TRUE)
    
    # Second loop on all species to build the forest to disturb
    for(j in 1:length(species.in)){
      
      # Extract the equilibrium for species i
      equil.j = sim.equil %>%
        filter(var == "n", equil, species == species.in[j]) %>% 
        pull(value)
      
      # Initiate the population at equilibrium
      list.species[[j]]$init_pop <- def_init_k(equil.j*0.03)
      
      # Update disturbance function
      list.species[[j]]$disturb_fun <- disturb_fun
      
      # Add disturbance coefficients
      list.species[[j]]$disturb_coef <- filter(matreex::disturb_coef, 
                                               species == species.in[j])
    }
    
    # Make forest
    forest.in = new_forest(species = list.species, harv_rules = harv_rules.ref)
    
    # Run simulation till equilibrium
    sim.out = sim_deter_forest(
      forest.in, tlim = 4000, equil_time = 4000, disturbance = disturbance.df, 
      SurfEch = 0.03, verbose = TRUE)
    
  } else {
    sim.out = matrix()
  }
  
  
  # Save simulation in a rdata
  create_dir_if_needed(subset(forest_list_noneq, ID.forest == ID.forest.in)$file.sim.dist)
  saveRDS(sim.out, subset(forest_list_noneq, ID.forest == ID.forest.in)$file.sim.dist)
  
  # Return output list
  return(subset(forest_list_noneq, ID.forest == ID.forest.in)$file.sim.dist)
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




#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Below: FUNCTIONS FOR THE REVISION TO FUN ECOL
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



#' Function to load and pre-format data from FUNDIV files
#' @param FUNDIV_tree_file Location of the file containing FUNDIV tree data
#' @param FUNDIV_plot_file Location of the file containing FUNDIV plot data
#' @param FUNDIV_climate_file Location of the file containing FUNDIV climate data
#' @param FUNDIV_species_file Location of the file containing FUNDIV species data
read_FUNDIV = function(FUNDIV_tree_file, FUNDIV_plot_file, 
                       FUNDIV_climate_file, FUNDIV_species_file){
  
  # Read climatic data
  FUNDIV_climate = fread(FUNDIV_climate_file) %>%
    # Remove Belgium
    filter(country != "WA") %>%
    # Only keep the mean sgdd and wai
    dplyr::select(plotcode, sgdd = sgdd_moreno_m, wai = wai_moreno_m) %>%
    # Remove null sgdd (points too close from the sea)
    filter(sgdd != 0)
  
  # Read tree, species and plot data
  FUNDIV_plot = fread(FUNDIV_plot_file)
  FUNDIV_tree = fread(FUNDIV_tree_file)
  FUNDIV_species = fread(FUNDIV_species_file)
  
  # Correct the weight in FUNDIV_tree
  FUNDIV_tree = FUNDIV_tree %>%
    mutate(weight1 = case_when(country == "DE" ~ weight1, 
                               country == "FG" ~ ba_ha1/ba1,
                               country %in% c("ES", "FI", "SW", "WA") ~ 10000/(pi*weight1^2)), 
           weight2 = case_when(country == "DE" ~ weight2, 
                               country == "FG" ~ ba_ha2/ba2,
                               country %in% c("ES", "FG", "FI", "SW", "WA") ~ 10000/(pi*weight2^2))) %>% 
    mutate_if(is.numeric, list(~na_if(., Inf)))
  
  
  # Correct  species anomalies 
  FUNDIV_tree$speciesid[FUNDIV_tree$speciesid %in% c(46,47)] <- 48
  FUNDIV_species$species[FUNDIV_species$id ==277] <- "pubescens"
  FUNDIV_species <- FUNDIV_species %>% mutate(sp = paste(genus, species)) %>%
    dplyr::select(c(id, sp))
  FUNDIV_species[FUNDIV_species$id == 277, "sp"] <- "Quercus pubescens"
  FUNDIV_species[FUNDIV_species$id == 48, "sp"] <- "Betula"
  
  ## - Make pca with sgdd and wai
  pca <- prcomp((FUNDIV_climate %>% dplyr::select(sgdd, wai) %>% na.omit()), 
                center = TRUE, scale = TRUE)
  
  ## - Format dataset
  out <- FUNDIV_tree %>%
    # Remove belgium and points in the sea
    filter(plotcode %in% FUNDIV_climate$plotcode) %>%
    # Add species name
    left_join((FUNDIV_species %>% 
                 na.omit() %>%
                 dplyr::select(speciesid = id, species = sp)), 
              by = "speciesid") %>% 
    # Add plot coordinates and dates
    left_join((FUNDIV_plot %>% 
                 dplyr::select(plotcode, longitude, latitude, yearsbetweensurveys)), 
              by = "plotcode") %>%
    # Add climate
    left_join((FUNDIV_climate %>% 
                 dplyr::select(plotcode, sgdd, wai) %>%
                 na.omit() %>%
                 mutate(pca1 = get_pca_ind(pca)[[1]][, 1], 
                        pca2 = get_pca_ind(pca)[[1]][, 2])), 
              by = "plotcode") %>%
    # Calculate competition
    group_by(plotcode) %>% mutate(BAtot = sum(ba_ha1) - ba_ha1) %>% 
    ungroup %>% group_by(plotcode, species) %>% 
    mutate(BAtotSP = sum(ba_ha1) - ba_ha1, BAtotNONSP = BAtot - BAtotSP) %>%
    ungroup()
  
  # Return output
  return(out)
  
}

#' Function to create a dataset with only climate and presence of species
#' @param FUNDIV_data Pre-formatted FUNDIV data
get_FUNDIV_species_per_climate = function(FUNDIV_data){
  
  FUNDIV_data %>%
    # Keep columns of interest
    dplyr::select(plotcode, longitude, latitude, sgdd, wai, pca1, pca2, species) %>%
    # Keep only species in the IPM
    filter(species %in% gsub("\\_", "\\ ", unique(climate_species$sp))) %>%
    # Remove duplicates and na
    distinct() %>%
    na.omit() %>%
    # Transform species presence in a binary variable per column
    mutate(value = 1, species = gsub("\\ ", "\\_", species)) %>%
    tidyr::spread(key = "species", value = "value") %>% 
    replace(is.na(.), 0)
  
}


#' Function to extract recruitment traits from demographic parameters
#' @param fit.list.allspecies demographic parameters of all species
#' @param FUNDIV_data pre-formatted tree data from FUNDIV
#' @param comp.ref character indicating how to calculate competition: 
#'                 "same" = mean competition across all species in dataset
#'                 "specific" = mean competition per species
get_recruitment_traits = function(fit.list.allspecies, FUNDIV_data, comp.ref){
  
  # Remove in-growth
  FUNDIV_data = subset(FUNDIV_data, treestatus_th != 1)
  
  # Calculate the mean competition in the entire dataset
  BASP.mean = mean(FUNDIV_data$BAtotSP, na.rm = TRUE)
  BANONSP.mean = mean(FUNDIV_data$BAtotNONSP, na.rm = TRUE)
  
  # If we choose to attribute the same competition for all species: 
  if(comp.ref == "same"){
    # Provide mean value of competition across all species
    data_species = data.frame(species = unique(climate_species$sp), 
                              BATOTSP = BASP.mean, 
                              BATOTNonSP = BANONSP.mean, 
                              logBATOTSP = log(BASP.mean), 
                              intercept = 1) %>%
      # Add climatic data
      left_join((climate_species %>%
                   filter(N == 2) %>%
                   dplyr::select(species = sp, wai, wai2, 
                                 waib, sgdd, sgdd2, sgddb)), 
                by = "species")
  }
  
  # If we choose to attribute different competition per species
  if(comp.ref == "specific"){
    # Calculate mean competition per species
    data_species = FUNDIV_data %>%
      mutate(species = gsub("\\ ", "\\_", species)) %>%
      filter(species %in% unique(climate_species$sp)) %>%
      group_by(species) %>%
      summarize(BATOTSP = mean(BAtotSP, na.rm = TRUE), 
                BATOTNonSP = mean(BAtotNONSP, na.rm = TRUE)) %>%
      mutate(logBATOTSP = log(BATOTSP), 
             intercept = 1) %>%
      # Add climatic data
      left_join((climate_species %>%
                   filter(N == 2) %>%
                   dplyr::select(species = sp, wai, wai2, 
                                 waib, sgdd, sgdd2, sgddb)), 
                by = "species")
  }
  
  # Initialize the output dataset
  traits_rec = data.frame(species = names(fit.list.allspecies), 
                          recruitment = NA_real_, 
                          delay = NA_real_)
  
  # Loop on all species to gather traits
  for(i in 1:dim(traits_rec)[1]){
    
    # Species i
    sp.i = traits_rec$species[i]
    
    # Give the delay from fit list
    traits_rec$delay[i] = as.numeric(fit.list.allspecies[[i]]$info["delay"])
    
    # Recruitment parameters for species i
    vec.rec.i = fit.list.allspecies[[i]]$rec$params_m
    
    # Associated vector of variables
    vec.var.i = as.vector(subset(data_species, species == sp.i)[, names(vec.rec.i)])
    
    # Get the recruitment
    traits_rec$recruitment[i] = exp(sum(vec.rec.i*vec.var.i))
    
  }
  
  # Return the trait dataset generated
  return(traits_rec)
}


#' Get coordinate in first pca traits axis per species
#' @param traits dataframe containing trait value per species
#' @param traits_rec recruitment traits per species
get_pc12_per_species <- function(traits, traits_rec){
  
  # Compile the traits data
  data_traits = left_join(traits, traits_rec, by = "species") %>%
    dplyr::select(-delay) %>%
    drop_na()
  
  # Make the PCA
  pca <- prcomp((data_traits %>% dplyr::select(-species)), 
                center = T, scale = T)
  
  # Extract data for individuals
  out = data.frame(species = data_traits$species, 
                   pca1 = get_pca_ind(pca)[[1]][, 1], 
                   pca2 = get_pca_ind(pca)[[1]][, 2])
  # return the output
  return(out)
}



#' Plot traits PCA
#' @param traits dataframe containing trait value per species
#' @param traits_rec recruitment traits per species
#' @param species_list dataset generated by get_species_list() function
#' @param file.in name including path of the file to save
plot_traits_pca12 <- function(traits, traits_rec, species_list, file.in){
  
  # Create directory if needed
  create_dir_if_needed(file.in)
  
  # Compile the traits data
  data_traits = left_join(traits, traits_rec, by = "species") %>%
    drop_na() %>%
    dplyr::select(-delay)
  
  # Make the PCA
  pca <- prcomp((data_traits %>% dplyr::select(-species)), 
                center = T, scale = T)
  
  # Extract data for individuals
  data.ind = data.frame(species = data_traits$species, 
                        pca1 = get_pca_ind(pca)[[1]][, 1], 
                        pca2 = get_pca_ind(pca)[[1]][, 2])  %>%
    filter(species %in% unique(species_list$species)) %>%
    mutate(species = gsub("\\_", "\\ ", species))
  
  # Range of pca1 and pca2
  range.pca1 = diff(range(data.ind$pca1))
  range.pca2 = diff(range(data.ind$pca2))
  
  # Extract data for variables
  data.var = data.frame(var = rownames(get_pca_var(pca)[[1]]), 
                        pca1 = get_pca_var(pca)[[1]][, 1], 
                        pca2 = get_pca_var(pca)[[1]][, 2]) %>%
    mutate(var = gsub("\\.", "\\ ", var), 
           pca1 = pca1*(max(abs(data.ind$pca1))/max(abs(pca1)))*0.9, 
           pca2 = pca2*(max(abs(data.ind$pca2))/max(abs(pca2)))*0.9, 
           pca1.txt = pca1*1.2, pca2.txt = pca2*1.2) 
  
  # Space fraction between axis and text
  space.frac = 0.025
  
  # Range of x and y axis for plotting
  range.x = range(data.var$pca1.txt)*1.5
  range.y = range(data.var$pca2.txt)*1.1
  
  # Make the plot
  plot.out = data.ind %>%
    # Add coordinates for species along x axis
    left_join((data.ind %>%
                 arrange(pca1) %>%
                 mutate(pca1_x = seq(from = range.x[1]+0.1, to = range.x[2]-0.1, 
                                     length.out = dim(.)[1]), 
                        pca1_y = range.y[2] + space.frac*diff(range.y)) %>%
                 dplyr::select(species, pca1_x, pca1_y)), by = "species") %>%
    # Add coordinates for species along y axis
    left_join((data.ind %>%
                 arrange(pca2) %>%
                 mutate(pca2_y = seq(from = range.y[1]+0.1, to = range.y[2]-0.1, 
                                     length.out = dim(.)[1]), 
                        pca2_x = range.x[2] + space.frac*diff(range.x)) %>%
                 dplyr::select(species, pca2_x, pca2_y)), by = "species") %>%
    ggplot(aes(x = pca1, y = pca2)) + 
    geom_segment(data = data.var, aes(x = 0, xend = pca1, y = 0, yend = pca2), 
                 arrow = arrow(length = unit(0.3, "cm"))) + 
    geom_point(size = 2, shape = 21, fill = "#22333B", color = "black", alpha = 0.5) +
    geom_text(data = data.var, aes(label = var, x = pca1.txt, y = pca2.txt)) +
    # Frame of the graph
    geom_rect(aes(xmin = range.x[1], xmax = range.x[2], ymin = range.y[1], ymax = range.y[2]), 
              color = "black", fill = NA) +
    # Zero horizontal and vertical lines
    geom_segment(x = range.x[1], xend = range.x[2], y=0, yend=0, linetype = "dashed") +
    geom_segment(x = 0, xend = 0, y = range.y[1], yend = range.y[2], linetype = "dashed") +
    # Segment to connect points to species name
    geom_segment(aes(xend = pca1_x), yend=range.y[2], linetype = "dotted", 
                 inherit.aes = TRUE, color = "#ADB5BD", size = 0.3) +
    geom_segment(aes(yend = pca2_y), xend=range.x[2], linetype = "dotted", 
                 inherit.aes = TRUE, color = "#ADB5BD", size = 0.3) +
    # Axis label
    xlab(paste0("PCA1 (", round(summary(pca)$importance[2, 1]*100, digits = 2), "%)",
                "\nHigh survival <--> High growth")) +
    ylab(paste0("PCA2 (", round(summary(pca)$importance[2, 2]*100, digits = 2), "%)",
                "\nHigh recruitment <--> Low recruitment")) +
    # Text axis
    # -- x axis
    geom_text(data = (data.frame(pca1 = c(-10:10), pca2 = range.y[1] - space.frac*diff(range.y)) %>%
                        filter(pca1 > range.x[1] & pca1 < range.x[2])), 
              aes(label = pca1), color = "#6C757D", size = 3) +
    # -- y axis
    geom_text(data = (data.frame(pca1 = range.x[1] - space.frac*diff(range.x), pca2 = c(-10:10)) %>%
                        filter(pca2 > range.y[1] & pca2 < range.y[2])), 
              aes(label = pca2), color = "#6C757D", size = 3) +
    # text of the x top axis
    geom_text(angle = 90, size = 2, 
              aes(x = pca1_x, y = pca1_y, label = species), color = "#6C757D", 
              inherit.aes = TRUE, vjust = 0, hjust = 0, fontface = "italic") + 
    # text of the y top axis
    geom_text(size = 2, 
              aes(x = pca2_x, y = pca2_y, label = species), color = "#6C757D", 
              inherit.aes = TRUE, vjust = 0, hjust = 0, fontface = "italic") + 
    # Theme
    theme(panel.background = element_blank(), 
          panel.grid = element_blank(), 
          axis.ticks = element_blank(), 
          axis.text = element_blank(), 
          axis.title = element_blank()) + 
    # Limits of the plot
    xlim(c(1.2, 1.5)*range.x) + ylim(c(1.2, 1.5)*range.y) + 
    # X title manual
    geom_text(x = mean(range.x), y = range.y[1] - 3*space.frac*diff(range.y), 
              label = paste0("PCA1 (", round(summary(pca)$importance[2, 1]*100, digits = 2), "%)"), 
              size = 3) + 
    geom_text(x = mean(range.x), y = range.y[1] - 5*space.frac*diff(range.y), 
              label = "High survival <--> High growth", fontface = "bold", size = 3) + 
    # Y title manual
    geom_text(y = mean(range.y), x = range.x[1] - 5*space.frac*diff(range.x), 
              label = paste0("PCA2 (", round(summary(pca)$importance[2, 2]*100, digits = 2), "%)"), 
              size = 3, angle = 90) + 
    geom_text(y = mean(range.y), x = range.x[1] - 3*space.frac*diff(range.x), angle = 90,
              label = "High recruitment <--> Low recruitment", fontface = "bold", size = 3)
  
  # - Save the plot
  ggsave(file.in, plot.out, width = 12, height = 12, 
         units = "cm", dpi = 600, bg = "white")
  
  
  # return the name of all the plots made
  return(file.in)
}



#' Get functional diversity from simulations with FD package
#' @param forest_list df containing info on each forest simulated
#' @param sim_disturbance vector containing the file names of sim with disturbance
#' @param pc12_per_species Position of each species along the fecundity and gr-surv axes
get_FD_multivar <- function(forest_list, sim_disturbance, pc12_per_species){
  
  # Initialize vector of all successful simulations
  vec.sim = c()
  
  # Initialize the original fd data
  data.fd.original <- forest_list %>%
    dplyr::select(ID.forest, ID.climate, combination) %>%
    mutate(nsp = NA_real_, FD1 = NA_real_, FD2 = NA_real_, H = NA_real_, 
           D = NA_real_, Nha = NA_real_)
  
  # Loop on all species combination
  for(i in 1:length(sim_disturbance)){
    
    # Printer
    print(paste0(i, "/", length(sim_disturbance)))
    
    # Fill the number of species
    data.fd.original$nsp[i] = length(unlist(strsplit(data.fd.original$combination[i], "\\.")))
    
    # Read simulation i
    sim.i = readRDS(sim_disturbance[i])
    
    # First, verify that equilibrium was reached
    if(!is.na(sim.i[1, 1])){
      
      ## - Calculate the number of trees per ha at equilibrium
      data.fd.original$Nha[i] = sum((sim.i %>%
                                       filter(var == "N") %>%
                                       filter(time == 1))$value)
      
      ## - FD with the original approach
      data.fd.original.i <- sim.i %>%
        filter(var == "BAsp") %>%
        filter(time == 1) %>%
        left_join(pc12_per_species, by = "species") %>%
        mutate(p = value/sum(.$value), 
               plnp = p*log(p), 
               p2 = p^2) %>%
        summarise(FD1 = weighted.var(pca1, w = value), 
                  FD2 = weighted.var(pca2, w = value), 
                  H = -sum(plnp), 
                  D = 1/sum(p2)) 
      data.fd.original$FD1[i] = data.fd.original.i$FD1
      data.fd.original$FD2[i] <- data.fd.original.i$FD2
      data.fd.original$H[i] <- data.fd.original.i$H
      data.fd.original$D[i] <- data.fd.original.i$D
      
      ## - FD with the FD package (dimension approach)
      data.fd.dimension.i <- sim.i %>%
        mutate(ID.climate = forest_list$ID.climate[i], 
               ID.forest = forest_list$ID.forest[i], 
               ID.community = paste(ID.climate, ID.forest, sep = ".")) %>%
        filter(var == "BAsp") %>%
        filter(time == 1) %>%
        dplyr::select(ID.climate, ID.forest, ID.community, species, value)
      
      # Also check we have trait value for all species in the community
      if(all(data.fd.dimension.i$species %in% pc12_per_species$species)){
        # Add to the final dataframe
        if(length(vec.sim) == 0) data.fd.dimension = data.fd.dimension.i
        else data.fd.dimension = rbind(data.fd.dimension, data.fd.dimension.i)
        
        # Increment the counter
        vec.sim = c(vec.sim, i)
      }
    }
    
  }
  
  # Replace NA by 0 in original approach
  data.fd.original <- data.fd.original %>% 
    mutate(FD1 = ifelse(is.na(FD1), 0, FD1), 
           FD2 = ifelse(is.na(FD2), 0, FD2))
  
  # Abundance dataframe
  abun.df = data.fd.dimension %>%
    spread(key = "species", value = "value") %>% 
    replace(is.na(.), 0)
  
  # Abundance matrix
  abun.matrix = as.matrix(
    abun.df %>% dplyr::select(-ID.climate, -ID.forest, -ID.community))
  rownames(abun.matrix) = abun.df$ID.community
  
  # Trait df
  trait.df = data.frame(species = colnames(abun.matrix)) %>%
    left_join(pc12_per_species, by = "species") %>%
    dplyr::select(-species)
  rownames(trait.df) = colnames(abun.matrix)
  
  # Calculate FD per community
  fd.raw <- dbFD(trait.df, abun.matrix, w.abun = TRUE)
  
  # Final dataframe
  data.fd.dimension.final <- data.frame(
    ID.forest = forest_list$ID.forest[vec.sim], 
    ID.climate = forest_list$ID.climate[vec.sim], 
    combination = forest_list$combination[vec.sim], 
    FRic = fd.raw$FRic, 
    FDis = fd.raw$FDis, 
    FEve = fd.raw$FEve, 
    FDiv = fd.raw$FDiv, 
    CWM1 = fd.raw$CWM$pca1, 
    CWM2 = fd.raw$CWM$pca2
  )
  
  # Join the two datasets
  out = left_join(data.fd.original, data.fd.dimension.final, 
                  by = c("ID.forest", "ID.climate", "combination")) 
  
  # Return output
  return(out)
}



#' Function to plot the structure and composition of calibration data
#' @param FUNDIV_data Pre-formatted FUNDIV dataset
#' @param pc1_per_species position of each species on the gr-surv avis
#' @param climate_list list of the quantile of each climate selected
#' @param file.in Name of the file to save, including path
plot_clim_vs_div_and_str_data = function(
  FUNDIV_data, pc1_per_species, climate_list, file.in){
  
  # Create output directory if needed
  create_dir_if_needed(file.in)
  
  # Species composition data
  data.sp = FUNDIV_data %>%
    filter(treestatus_th != 1) %>%
    mutate(species = gsub("\\ ", "\\_", species)) %>%
    group_by(plotcode, species) %>%
    summarise(ba_ha = sum(ba_ha1, na.rm = TRUE)) %>%
    left_join((pc1_per_species %>% rename(pca1_trait = pca1)), 
              by = "species") %>%
    ungroup() %>% group_by(plotcode) %>%
    mutate(p = ba_ha/sum(ba_ha, na.rm = TRUE), 
           plnp = p*log(p)) %>%
    summarise(FD = weighted.var(pca1_trait, w = ba_ha, na.rm = TRUE),
              CWM = weighted.mean(pca1_trait, w = ba_ha, na.rm = TRUE),
              H = -sum(plnp)) %>%
    mutate(FD = ifelse(H == 0, 0, FD))
  
  # Forest structure data
  data.str = FUNDIV_data %>%
    filter(treestatus_th != 1) %>%
    group_by(plotcode) %>%
    mutate(n_ha1 = ifelse(ba1 == 0, NA_real_, ba_ha1/ba1)) %>%
    summarise(pca1 = mean(pca1, na.rm = TRUE), 
              BA = sum(ba_ha1, na.rm = TRUE), 
              dbh_mean = weighted.mean(dbh1, w = weight1), 
              Nha = sum(n_ha1))
  
  # Merge data together
  data.out = data.str %>%
    left_join(data.sp, by = "plotcode") %>%
    mutate(clim = NA_real_)
  
  # Loop on all climates to fill the climate ID
  for(i in 1:length(names(climate_list))){
    range.i = quantile(data.out$pca1, probs = climate_list[[i]])
    data.out = data.out %>% mutate(clim = ifelse(
      (pca1 >= range.i[1] & pca1 < range.i[2]), paste0("clim", i), clim))
  }
  
  # Make final plot
  plot.out = data.out %>%
    filter(!is.na(clim)) %>%
    mutate(clim = factor(clim, levels = paste0("clim", c(1:10)))) %>%
    dplyr::select(clim, H, FD, CWM, BA, dbh_mean, Nha) %>%
    gather(key = "variable", value = "value", "H", "FD", "CWM", "BA", "Nha", 
           "dbh_mean") %>%
    # Only keep variables of interest
    filter(variable %in% c("BA", "dbh_mean", "Nha")) %>%
    mutate(variable = factor(variable, levels = c("BA", "dbh_mean", "Nha"))) %>%
    ggplot(aes(x = clim, y = value, fill = clim)) + 
    geom_boxplot(alpha = 0.7) + 
    scale_fill_manual(values = colorRampPalette(c("orange", "blue"))(10)) +
    facet_wrap(~ variable, scales = "free") + 
    theme(legend.position = "none", 
          panel.grid = element_blank(), 
          panel.background = element_rect(fill = "white", color = "black"), 
          strip.background = element_blank(), 
          axis.text.x = element_text(angle = 60, hjust = 1), 
          axis.title = element_blank()) 
  
  # Save the plot
  ggsave(file.in, plot.out, width = 17, height = 8, 
         units = "cm", dpi = 600, bg = "white")
  
  # Return the name of the file
  return(file.in)
}



#' Function to estimate of FD metrics on resilience
#' @param data_model df formatted to fit model
#' @param dir.in Name of the directory where to save files
plot_FD_effect_resilience_multivar = function(data_model, dir.in){
  
  # Name of the files
  fig.file.in = paste0(dir.in, "/fig_FD_effect_resilience_storm.jpg")
  fig.file.predictions = paste0(dir.in, "/observed_vs_predicted_H1.jpg")
  fig.file.residuals = paste0(dir.in, "/residuals.jpg")
  fig.file.observations = paste0(dir.in, "/observations.jpg")
  table.file.stats = paste0(dir.in, "/stats_H1.tex")
  
  # create output directory if it doesn't exist
  create_dir_if_needed(fig.file.in)
  
  # Vector of response variables for which to run models
  response.vec = c("resistance", "recovery", "resilience")
  
  # Remove from data model data with na
  data_model = data_model %>% filter(!is.na(CWM1))
  
  # Data to fit the models
  data.in = cbind(data_model[, response.vec], 
                  scale((data_model %>% dplyr::select(
                    "H.scaled" = "H", "FD.scaled" = "FDis", 
                    "CWM1.scaled" = "CWM1", "CWM2.scaled" = "CWM2")), 
                    center = TRUE, scale = TRUE)) %>%
    mutate(resilience = log(resilience), 
           recovery = log(recovery)) %>%
    drop_na()
  
  # Initialize the table with statistics
  table.stats = data.frame(col1 = "", col2 = "Est", col3 = "se", col4 = "F", 
                           col5 = "p", col6 = "VIF")
  
  # Initialize the list that will contain graphs composing the final plot
  plotlist.out = list()
  
  # Loop on all response variables
  for(j in 1:length(response.vec)){
    
    # Fit model
    eval(parse(text = paste0(
      "model.j = lm(", response.vec[j], 
      " ~ H.scaled + FD.scaled + CWM1.scaled + CWM2.scaled, data = data.in)"
    )))
    
    # Data for predictions vs observations
    data.predict.j = data.in %>%
      mutate(predicted = predict(model.j, newdata = .)) %>%
      rename("observed" = response.vec[j]) %>%
      mutate(response.var = response.vec[j], 
             residuals = model.j$residuals, 
             H = data_model$H, FD = data_model$FDis, CWM1 = data_model$CWM1,
             CWM2 = data_model$CWM2) %>%
      dplyr::select(response.var, H, H.scaled, FD, FD.scaled, CWM1, CWM1.scaled, 
                    CWM2, CWM2.scaled, predicted, observed, residuals)
    
    # Models to unscale explanatory variables
    unscale.H = lm(H ~ H.scaled, data.predict.j)
    unscale.FD = lm(FD ~ FD.scaled, data.predict.j)
    unscale.CWM1 = lm(CWM1 ~ CWM1.scaled, data.predict.j)
    unscale.CWM2 = lm(CWM2 ~ CWM2.scaled, data.predict.j)
    
    # Data for fit
    for(i in 1:4){
      var.exp.i = c("H", "FD", "CWM1", "CWM2")[i]
      data.fit.i = data.frame(
        H.scaled = seq(from = min(data.in$H.scaled), to = max(data.in$H.scaled), length.out = 100), 
        FD.scaled = seq(from = min(data.in$FD.scaled), to = max(data.in$FD.scaled), length.out = 100), 
        CWM1.scaled = seq(from = min(data.in$CWM1.scaled), to = max(data.in$CWM1.scaled), length.out = 100), 
        CWM2.scaled = seq(from = min(data.in$CWM2.scaled), to = max(data.in$CWM2.scaled), length.out = 100), 
        exp.var = var.exp.i, 
        response.var = response.vec[j]) %>%
        # Set to the mean for variables that are not var i
        mutate(H.scaled = ifelse(exp.var == "H", H.scaled, mean(data.predict.j$H.scaled)), 
               FD.scaled = ifelse(exp.var == "FD", FD.scaled, mean(data.predict.j$FD.scaled)), 
               CWM1.scaled = ifelse(exp.var == "CWM1", CWM1.scaled, mean(data.predict.j$CWM1.scaled)), 
               CWM2.scaled = ifelse(exp.var == "CWM2", CWM2.scaled, mean(data.predict.j$CWM2.scaled))) %>%
        # Add unscaled values
        mutate(H = predict(unscale.H, newdata = .), 
               FD = predict(unscale.FD, newdata = .), 
               CWM1 = predict(unscale.CWM1, newdata = .), 
               CWM2 = predict(unscale.CWM2, newdata = .)) %>%
        cbind(predict(model.j, newdata = ., interval = "confidence")) %>%
        # Modify prediction
        mutate(fit = ifelse(response.var == "resistance", plogis(fit), exp(fit)), 
               lwr = ifelse(response.var == "resistance", plogis(lwr), exp(lwr)), 
               upr = ifelse(response.var == "resistance", plogis(upr), exp(upr))) %>%
        # Keep only columns needed
        dplyr::select("response.var", "exp.var", "exp.var.value" = var.exp.i, 
                      "observed" = "fit", "lwr", "upr")
      if(i == 1 & j == 1) data.fit = data.fit.i
      else data.fit = rbind(data.fit, data.fit.i)
      
    }
    
    # Output data set for model i j 
    data.out.j = data.frame(
      var.resp = response.vec[j], 
      var.exp = c("H", "FD", "CWM1", "CWM2"), 
      var.pos = c(1:4),
      est = as.numeric(coef(model.j)[-1]), 
      est.low = as.numeric(confint(model.j)[-1, 1]), 
      est.high = as.numeric(confint(model.j)[-1, 2])
    )
    
    # Complete the statistics table
    table.stats = table.stats %>%
      rbind(data.frame(
        col1 = c("", response.vec[j], data.out.j$var.exp), 
        col2 = c("", "(R2: ", round(as.numeric(coef(model.j)[-1]), digits = 2)), 
        col3 = c("", paste0(round(summary(model.j)$r.squared, digits = 2), ")"), 
                 round(summary(model.j)$coefficients[-1, 2], digits = 2)), 
        col4 = c("", "", round(anova(model.j)[-dim(anova(model.j))[1], 4], digits = 2)), 
        col5 = c("", "", scales::pvalue(anova(model.j)[-dim(anova(model.j))[1], 5])), 
        col6 = c("", "", round(vif(model.j), digits = 2))
      ))
    
    # Range of estimate value for model j
    range.j = max(data.out.j$est.high) - min(data.out.j$est.low)
    
    # Plot for response variable j
    plot.j = data.out.j %>%
      mutate(significance = ifelse(est.low > 0 | est.high < 0, "yes", "no")) %>%
      mutate(var.exp = case_when(var.exp == "CWM1" ~ "CWM1\nSu->Gr", 
                                 var.exp == "CWM2" ~ "CWM2\nR+->R-", 
                                 TRUE ~ var.exp)) %>%
      ggplot(aes(x = var.exp, y = est, color = significance)) + 
      geom_errorbar(aes(ymin = est.low, ymax = est.high),
                    width = 0.1) + 
      geom_point(shape = 21, size = 3, fill = c("#9B2226", "#386641", "#33658A")[j]) +
      xlab("") + ylab(paste0("Effect on ", response.vec[j])) +
      scale_color_manual(values = c(`no` = "gray", `yes` = "black")) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      theme(panel.background = element_rect(color = "black", fill = "white"), 
            panel.grid = element_blank(), 
            strip.background = element_blank(), 
            strip.text = element_text(face = "bold"), 
            legend.position = "none", 
            axis.ticks.length=unit(-0.1, "cm")) + 
      coord_flip()
    
    # Add plot j to final list
    eval(parse(text = paste0("plotlist.out$", response.vec[j], " = plot.j")))
    
    # Add to the final output dataset
    if(j == 1){
      data.predict = data.predict.j
    } else {
      data.predict = rbind(data.predict, data.predict.j)
    } 
  }
  
  
  # Plot predicted vs observed
  plot.predictions = data.predict %>%
    mutate(observed = ifelse(response.var == "resistance", 
                             plogis(observed), exp(observed)), 
           predicted = ifelse(response.var == "resistance", 
                              plogis(predicted), exp(predicted))) %>%
    ggplot(aes(x = observed, y = predicted)) + 
    geom_point(shape = 21, color = "black", fill = "grey", alpha = 0.5) + 
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") + 
    facet_wrap(~ response.var, scales = "free") + 
    xlab("Observed value") + ylab("Predicted value") + 
    theme(panel.background = element_rect(fill = "white", color = "black"), 
          panel.grid = element_blank(), 
          strip.background = element_blank())
  
  # Plot residuals vs explanatory variables
  plot.residuals = data.predict %>%
    dplyr::select(response.var, residuals, H, FD, CWM1, CWM2) %>%
    gather(key = "exp.var", value = "exp.var.value", "H", "FD", "CWM1", "CWM2") %>%
    ggplot(aes(x = exp.var.value, y = residuals, fill = response.var, color = response.var)) + 
    scale_fill_manual(values = c("#9B2226", "#386641", "#33658A")) +
    scale_color_manual(values = c("#9B2226", "#386641", "#33658A")) +
    geom_point(shape = 21, color = "black", alpha = 0.5) + 
    facet_grid(response.var ~ exp.var, scales = "free") + 
    geom_hline(yintercept = 0, linetype = "dashed") + 
    geom_smooth(method = "loess") +
    xlab("Value of species composition metric") +
    theme(panel.background = element_rect(color = "black", fill = "white"), 
          panel.grid = element_blank(), 
          strip.background = element_blank(), 
          strip.text = element_text(face = "bold"), 
          legend.position = "none")
  
  # Plot observations along with prediction
  plot.observed = data.predict %>%
    dplyr::select(response.var, observed, predicted, H, FD, CWM1, CWM2) %>%
    gather(key = "exp.var", value = "exp.var.value", "H", "FD", "CWM1", "CWM2") %>%
    mutate(observed = ifelse(response.var == "resistance", 
                             plogis(observed), exp(observed))) %>%
    ggplot(aes(x = exp.var.value, y = observed, fill = response.var, 
               color = response.var, group = interaction(response.var, exp.var))) + 
    geom_ribbon(data = data.fit, aes(ymin = lwr, ymax = upr), inherit.aes = TRUE, 
                alpha = 0.2, color = NA) +
    geom_line(data = data.fit, inherit.aes = TRUE, size = 1) +
    geom_smooth(method = "loess", size = 1, color = "black") +
    scale_fill_manual(values = c("#9B2226", "#386641", "#33658A")) +
    scale_color_manual(values = c("#9B2226", "#386641", "#33658A")) +
    geom_point(shape = 21, color = "black", alpha = 0.4) + 
    facet_grid(response.var ~ exp.var, scales = "free") + 
    xlab("Value of species composition metric") +
    theme(panel.background = element_rect(color = "black", fill = "white"), 
          panel.grid = element_blank(), 
          strip.background = element_blank(), 
          strip.text = element_text(face = "bold"), 
          legend.position = "none")
  
  # Plot the estimates
  plot.out = plot_grid(plotlist = plotlist.out, nrow = 1, scale = 0.9, 
                       labels = c("(a)", "(b)", "(c)"))
  
  # Save plot 
  ggsave(fig.file.in, plot.out, width = 22, height = 8, units = "cm", 
         dpi = 600, bg = "white")
  ggsave(fig.file.predictions, plot.predictions, width = 14, height = 5, 
         units = "cm", dpi = 600, bg = "white")
  ggsave(fig.file.residuals, plot.residuals, width = 21, height = 15, 
         units = "cm", dpi = 600, bg = "white")
  ggsave(fig.file.observations, plot.observed, width = 21, height = 15, 
         units = "cm", dpi = 600, bg = "white")
  
  # Save tex file
  print(xtable(table.stats, type = "latex", 
               caption = "Statistics of the models used to test the effect of species composition on forest response to disturbances (H1)", 
               label = "table_stat_H1"), 
        include.rownames=FALSE, hline.after = c(0, 1, dim(table.stats)[1]), 
        include.colnames = FALSE, caption.placement = "top", 
        file = table.file.stats)
  
  
  # Return name of the file saved
  return(c(fig.file.in, fig.file.predictions, table.file.stats))
}




#' Plot FD effect on resilienc along climate with quadratic term
#' @param data_model df formatted to fit model
#' @param R_metric Richness metric to choose ("nsp", "H" or "D")
#' @param mod_selection character to specify if model should be selected using 
#'                      AIC ("AIC") or p-values ("pvalue")
#' @param file.in Name of the file to save, inlcuding path
plot_FD_effect_vs_climate_multivar = function(
  data_model, R_metric = "H", mod_selection = "pvalue", dir.in){
  
  # create output directory if it doesn't exist
  create_dir_if_needed(paste0(dir.in, "/test"))
  
  # Vector of response variables for which to run models
  response.vec = c("resistance", "recovery", "resilience")
  
  # Vector of explanatory variables that are not climate
  vec.exp = c("H", "FD", "CWM1", "CWM2")
  
  # Initialize data with the effect pvalue
  data.text = expand.grid(var.resp = response.vec, 
                          var.exp = vec.exp, 
                          text = NA)
  
  # Data to fit the models for disturbance i
  data.in = cbind(data_model[, response.vec], 
                  scale((data_model %>% 
                           dplyr::select("H" = R_metric, "FD" = "FDis", "CWM1", 
                                         "CWM2", "Clim" = "pca1")), 
                        center = TRUE, scale = TRUE)) %>%
    # Convert strictly positive variables to log scale
    mutate(resilience = log(resilience), 
           recovery = log(recovery), 
           Clim2 = Clim^2)
  
  # Model to scale climate
  scale.clim = lm(pca1 ~ Clim, 
                  data = data.frame(pca1 = data_model$pca1, 
                                    Clim = data.in$Clim))
  
  # Initialize data with climatic gradient
  data.clim = data.frame(
    Clim = seq(from = min(data.in$Clim), to = max(data.in$Clim), length.out = 100)) %>%
    mutate(pca1 = predict(scale.clim, newdata = .), 
           pca1.scaled.2 = Clim^2) %>%
    rename(pca1.scaled = Clim)
  
  # Vector that contains interaction terms
  vec.interaction = c("", "*Clim", "*Clim2")
  
  
  # Initialize list of models
  model.list = list()
  
  # Initialize dataset contianing prediction vs data
  data.predict = data.frame(resp.var = character(0), 
                            observed = numeric(0), 
                            fitted = numeric(0))
  
  # Initialize dataset to plot model predictions against climate
  data.fit.clim = data.frame(
    var.resp = character(0), pca1 = numeric(0), fit = numeric(0), 
    lwr = numeric(0), upr = numeric(0)
  )
  
  # All possible factor combination to include in model
  # -- All possible interactions to include in the model
  possible.interactions = c(paste0(vec.exp, "*Clim"), paste0(vec.exp, "*Clim2"))
  # -- All possible combinations of these interactions
  combi.interactions = expand.grid(rep(list(0:1), length(possible.interactions)))
  # -- Create vector of all factor combinations
  combi.formula = data.frame(model = paste0("model", c(1:dim(combi.interactions)[1])), 
                             formula = array(paste(vec.exp, collapse = " + "), 
                                             dim = dim(combi.interactions)[1]), 
                             nvar = NA_real_, AIC = NA_real_)
  # -- Loop on all possible combination of interactions
  for(k in 1:dim(combi.interactions)[1]){
    if(sum(combi.interactions[k, ]) > 0){
      combi.formula$formula[k] = paste(
        c(combi.formula$formula[k], possible.interactions[which(combi.interactions[k, ] == 1)]), 
        collapse = " + ")
    }
  }
  
  # Initialize the list of all models
  list.models = list()
  
  # Initialize dataset containing metric vs climate data + fit
  # Loop on all response variables
  for(j in 1:length(response.vec)){
    
    # Initialize the table with formulas
    combi.formula.j = combi.formula
    
    # Initialize the list of models for response var j
    eval(parse(text = paste0("list.models$", response.vec[j], " = list()")))
    
    # Loop on all possible formulas
    for(i in 1:dim(combi.formula.j)[1]){
      
      # Run model with response var j and formula i
      eval(parse(text = paste0("model.ij = lm(", response.vec[j], " ~ ", 
                               combi.formula.j$formula[i], ", data = data.in)")))
      
      # Add AIC in the table
      combi.formula.j$AIC[i] = AIC(model.ij)
      # Add the number of variables in the table
      combi.formula.j$nvar[i] = length(coefficients(model.ij)) - 1
      # Store model in model list
      eval(parse(text = paste0("list.models$", response.vec[j], "$", 
                               combi.formula.j$model[i], " = model.ij")))
    }
    
    # Spot the best model based on the AIC
    id.model.j = which(combi.formula.j$AIC == min(combi.formula.j$AIC))
    
    # Assign the reference model 
    model.j = list.models[[j]][[id.model.j]]
    
    # Name of the stat tables to export 
    files.stat.table = paste0(dir.in, "/stat_H2_", response.vec, ".tex")
    
    ## Export stats of the model
    # -- Create table
    table.stat.j = data.frame(
      col1 = c("", names(coef(model.j)[-1])), 
      col2 = c("Est", round(as.numeric(coef(model.j)[-1]), digits = 2)), 
      col3 = c("se", round(summary(model.j)$coefficients[-1, 2], digits = 2)), 
      col4 = c("F", round(anova(model.j)[-dim(anova(model.j))[1], 4], digits = 2)), 
      col5 = c("p", scales::pvalue(anova(model.j)[-dim(anova(model.j))[1], 5])), 
      col6 = c("VIF", round(vif(model.j), digits = 2))
    )
    # -- Export the table
    print(xtable(table.stat.j, type = "latex",
                 caption = paste0("Statistics of the models used to test the ",
                                  "effect of climate
                                  and species composition on ",
                                  response.vec[j], " (H2)"),
                 label = paste0("table_stat_H1_", response.vec[j])),
          include.rownames=FALSE, hline.after = c(0, 1, dim(table.stat.j)[1]),
          include.colnames = FALSE, caption.placement = "top",
          file = files.stat.table[j])
    
    # Coefficients of the model
    beta.j = coef(model.j)
    
    # Variance covariance matrix of the model
    vcov.j = vcov(model.j)
    
    # Loop on all explanatory variables
    for(i2 in 1:length(vec.exp)){
      
      # Data for predictions
      # -- Initialize matrix with 0
      newdata.j = matrix(0, nrow = dim(data.clim)[1], ncol = length(beta.j), 
                         dimnames = list(c(), names(beta.j)))
      # -- Add response variable (set to 1) and wirte pvalue
      newdata.j[, vec.exp[i2]] = 1
      text_ij = paste0(vec.exp[i2], ": ", pvalue(
        summary(model.j)$coefficients[vec.exp[i2], 4], add_p = TRUE, accuracy = 0.01))
      # -- Add interaction with climate (equals 1 * climate scaled) if included in the model
      if(paste0(vec.exp[i2], ":Clim") %in% names(beta.j) | 
         paste0("Clim:", vec.exp[i2]) %in% names(beta.j)){
        if(paste0(vec.exp[i2], ":Clim") %in% names(beta.j)) int.name.ij = paste0(vec.exp[i2], ":Clim")
        if(paste0("Clim:", vec.exp[i2]) %in% names(beta.j)) int.name.ij = paste0("Clim:", vec.exp[i2])
        newdata.j[, int.name.ij] = data.clim$pca1.scaled
        # -- Also add text to the text vector
        text_ij = c(text_ij, paste0(vec.exp[i2], "*Clim: ", pvalue(
          summary(model.j)$coefficients[int.name.ij, 4], add_p = TRUE, accuracy = 0.01)))
      }
      # -- Add interaction with climate quadratic (equals 1 * climate scaled^2)
      if(paste0(vec.exp[i2], ":Clim2") %in% names(beta.j) | 
         paste0("Clim2:", vec.exp[i2]) %in% names(beta.j)){
        if(paste0(vec.exp[i2], ":Clim2") %in% names(beta.j)) int2.name.ij = paste0(vec.exp[i2], ":Clim2")
        if(paste0("Clim2:", vec.exp[i2]) %in% names(beta.j)) int2.name.ij = paste0("Clim2:", vec.exp[i2])
        newdata.j[, int2.name.ij] = data.clim$pca1.scaled.2
        # -- Also add text to the text vector
        text_ij = c(text_ij, paste0(vec.exp[i2], "*Clim2: ", pvalue(
          summary(model.j)$coefficients[int2.name.ij, 4], add_p = TRUE, accuracy = 0.01)))
      }
      
      
      # Predictions
      # -- Mean effect
      pred.ij = newdata.j %*% beta.j
      # -- Standard error
      pred.ij.se <- sqrt(diag(newdata.j %*% vcov.j %*% t(newdata.j)))
      # -- Statistical criteria
      alpha <- 0.05
      crit <- -qnorm(alpha/2)
      # -- Upper and lower interval
      lwr.ij <- pred.ij-crit*pred.ij.se 
      upr.ij <- pred.ij+crit*pred.ij.se 
      
      # Data ij
      data.ij = data.frame(pca1 = data.clim$pca1, 
                           effect = paste0(vec.exp[i2], "_effect"), 
                           mean = pred.ij, 
                           lwr = lwr.ij, 
                           upr = upr.ij, 
                           var.resp = response.vec[j])
      
      # Fill the text dataset
      id_text_ij = which(data.text$var.exp == vec.exp[i2] & data.text$var.resp == response.vec[j])
      data.text[id_text_ij, "text"] = paste(text_ij, collapse = "\n")
      
      
      # Add to the final output dataset
      if(j == 1 & i2 == 1) data.out = data.ij
      else data.out = rbind(data.out, data.ij)
    }
    
    # Prepare data to plot predicted resilience vs climate
    newdata.clim.j = matrix(0, nrow = dim(data.clim)[1], ncol = length(beta.j), 
                            dimnames = list(c(), names(beta.j)))
    # -- Add response variable (set to 1) and write pvalue
    if("Clim" %in% names(beta.j)) newdata.clim.j[, "Clim"] = data.clim$pca1.scaled
    if("Clim2" %in% names(beta.j)) newdata.clim.j[, "Clim2"] = data.clim$pca1.scaled.2
    newdata.clim.j[, "(Intercept)"] = 1
    pred.clim.j = newdata.clim.j %*% beta.j
    pred.se.clim.j <- sqrt(diag(newdata.clim.j %*% vcov.j %*% t(newdata.clim.j)))
    data.fit.clim.j = data.frame(
      var.resp = response.vec[j], 
      pca1 = data.clim$pca1,
      fit = pred.clim.j, 
      lwr = pred.clim.j-1.96*pred.se.clim.j, 
      upr = pred.clim.j+1.96*pred.se.clim.j 
    )
    data.fit.clim = rbind(data.fit.clim, data.fit.clim.j)
    
    # Add model to the output model list
    eval(parse(text = paste0("model.list$", response.vec[j], " = model.j")))
    
    # Add prediction vs observed
    data.predict = rbind(
      data.predict, 
      data.frame(resp.var = response.vec[j],
                 observed = data.in[, response.vec[j]], 
                 fitted = predict(model.j, newdata = data.in)))
  }
  
  # Extract coefficients of the simpler models
  for(j2 in 1:length(list.models)){
    data.simple.j2 = data.frame(
      var.resp = names(list.models)[j2], 
      pca1 = max(data.out$pca1) + 0.05*diff(range(data.out$pca1)), 
      effect = paste0(names(coefficients(list.models[[j2]][[1]]))[-1], " effect"), 
      mean = as.numeric(coefficients(list.models[[j2]][[1]])[-1]), 
      lwr = as.numeric(confint(list.models[[j2]][[1]])[-1, 1]), 
      upr = as.numeric(confint(list.models[[j2]][[1]])[-1, 2]), 
      clim = NA_character_
    )
    if(j2 == 1) data.simple = data.simple.j2
    else data.simple = rbind(data.simple, data.simple.j2)
  }
  data.simple = data.simple  %>%
    mutate(var.resp = factor(var.resp, levels = c("resistance", "recovery", "resilience"))) %>%
    mutate(effect = case_when(effect == "CWM1 effect" ~ "CWM1 effect\ngrowth <--> surv", 
                              effect == "CWM2 effect" ~ "CWM2 effect\nhigh R <--> low R", 
                              TRUE ~ effect)) %>%
    mutate(effect = factor(effect, levels = c("H effect", "FD effect", 
                                              "CWM1 effect\ngrowth <--> surv", 
                                              "CWM2 effect\nhigh R <--> low R")))
  
  
  
  # Plot the prediction and data of resilience vs climate
  plot.fit.clim = data.fit.clim %>%
    mutate(fit = ifelse(var.resp == "resistance", fit, exp(fit)), 
           lwr = ifelse(var.resp == "resistance", lwr, exp(lwr)), 
           upr = ifelse(var.resp == "resistance", upr, exp(upr))) %>%
    ggplot(aes(y = fit, x = pca1, group = 1)) + 
    geom_line() + 
    geom_point(data = (data_model %>%
                         gather(key = "var.resp", value = "fit", response.vec)), 
               inherit.aes = TRUE, 
               shape = 21, color = "black", fill = "grey", alpha = 0.5) +
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.5) + 
    facet_wrap(~ var.resp, nrow = 1, scales = "free") + 
    xlab("Coordinate on the wai-sgdd pca") + ylab("Value of\nresilience metric") + 
    theme(panel.background = element_rect(fill = "white", color = "black"), 
          panel.grid = element_blank(), 
          strip.background = element_blank())
  
  # Categorize climates
  data.clim = data.frame(pca1_min = quantile(data.ij$pca1, c(0:9)/10), 
                         pca1_max = quantile(data.ij$pca1, c(1:10)/10)) %>%
    mutate(clim = paste0("clim_", letters[c(1:dim(.)[1])])) %>%
    merge(data.frame(pca1 = data.ij[, "pca1"])) %>%
    filter(pca1 >= pca1_min & pca1 <= pca1_max) %>%
    dplyr::select(pca1, clim)
  
  # Add position on x and y axis for data.text
  data.text = data.text %>%
    mutate(pca1 = min(data.out$pca1) + 0.05*diff(range(data.out$pca1)), 
           mean = max(data.out$mean) + 0.14*diff(range(data.out$mean))) %>%
    mutate(effect = paste0(var.exp, " effect")) %>%
    mutate(effect = case_when(effect == "CWM1 effect" ~ "CWM1 effect\ngrowth <--> surv", 
                              effect == "CWM2 effect" ~ "CWM2 effect\nhigh R <--> low R", 
                              TRUE ~ effect)) %>%
    mutate(effect = factor(effect, levels = c("H effect", "FD effect", 
                                              "CWM1 effect\ngrowth <--> surv", 
                                              "CWM2 effect\nhigh R <--> low R"))) %>%
    mutate(clim = NA_character_)
  
  # Plot predicted vs observed
  plot.predictions = data.predict %>%
    mutate(observed = ifelse(resp.var == "resistance", 
                             plogis(observed), exp(observed)), 
           fitted = ifelse(resp.var == "resistance", 
                           plogis(fitted), exp(fitted))) %>%
    ggplot(aes(x = observed, y = fitted)) + 
    geom_point(shape = 21, color = "black", fill = "grey", alpha = 0.5) + 
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") + 
    facet_wrap(~ resp.var, scales = "free") + 
    xlab("Observed value") + ylab("Predicted value") + 
    theme(panel.background = element_rect(fill = "white", color = "black"), 
          panel.grid = element_blank(), 
          strip.background = element_blank())
  
  # Re-format data.out
  data.out = data.out %>%
    mutate(effect = gsub("\\_", " ", effect)) %>%
    mutate(var.resp = factor(var.resp, levels = c("resistance", "recovery", "resilience"))) %>%
    left_join(data.clim, by = "pca1")
  data.out = data.out  %>%
    rbind.data.frame((
      data.out %>%
        # Keep one line per climate and effect
        group_by(clim) %>% mutate(clim.min = min(pca1)) %>%
        filter(pca1 == min(pca1)) %>%
        # Remove one line per effect (min one)
        ungroup() %>% group_by(effect, var.resp) %>%
        mutate(effect.min = min(pca1)) %>% filter(pca1 != effect.min) %>%
        ungroup() %>%
        # Add new climate
        left_join(data.frame(clim = paste0("clim_", letters[c(2:10)]), 
                             clim2 = paste0("clim_", letters[c(1:9)])), 
                  by = "clim") %>%
        dplyr::select("pca1", "effect", "mean", "lwr", "upr", "var.resp", "clim" = "clim2")
    ))
  
  # --- NEW APPROACH FOR THE GRAPH
  # Last formatting of data.out
  data.out = data.out %>%
    mutate(effect = case_when(effect == "CWM1 effect" ~ "CWM1 effect\ngrowth <--> surv", 
                              effect == "CWM2 effect" ~ "CWM2 effect\nhigh R <--> low R", 
                              TRUE ~ effect)) %>%
    mutate(var.resp = as.character(var.resp))
  # Initialize plotlist
  plotlist.out = list()
  y1 = 1
  y2 = 1
  # Loop on response vesponse vec
  for(y1 in 1:length(unique(data.out$var.resp))){
    # Initialize plot list for this variable
    plotlist.y1 = list()
    # Name of the response var
    resp.y1 = as.character(unique(data.out$var.resp)[y1])
    # Range for y axis
    range.y1 = c(min(subset(data.out, var.resp == resp.y1)$lwr), 
                 max(subset(data.out, var.resp == resp.y1)$upr))
    # Loop on all effects
    for(y2 in 1:length(unique(data.out$effect))){
      # Effect y2
      effect.y2 = unique(data.out$effect)[y2]
      # Make plot y2
      plot.y1.y2 = data.out %>%
        filter(var.resp == resp.y1 & effect == effect.y2) %>%
        ggplot(aes(x = pca1, y = mean, group = clim, fill = clim)) + 
        geom_ribbon(aes(ymin = lwr, ymax = upr), color = NA, alpha = 0.5) +
        geom_line(color = "#001524") + 
        xlab("") + 
        ylab(effect.y2) +
        geom_hline(yintercept = 0, linetype = "dashed") +
        geom_text(data = (data.text %>% 
                            filter(var.resp == resp.y1 & effect == effect.y2) %>%
                            mutate(mean = range.y1[2])), 
                  aes(label = text), inherit.aes = TRUE, 
                  hjust = "inward", size = 2.5, alpha = 0.8) +
        scale_fill_manual(values = colorRampPalette(c("orange", "blue"))(10)) +
        ylim(range.y1*c(1, 1.4)) +
        theme(panel.background = element_rect(color = "black", fill = "white"), 
              panel.grid = element_blank(), 
              strip.background = element_blank(), 
              strip.text = element_text(face = "bold"), 
              legend.position = "none", 
              axis.ticks.length=unit(-0.1, "cm")) 
      # Add to plotlist of y1
      eval(parse(text = paste0("plotlist.y1$plot", y2, " = plot.y1.y2")))
    }
    # Build plot y1
    plot.y1 = plot_grid(plotlist = plotlist.y1, nrow = 1, align = "v")
    # Add plot to the final plot list
    eval(parse(text = paste0("plotlist.out$void", y1, " = ggplot() + theme_void()")))
    eval(parse(text = paste0("plotlist.out$plot", y1, " = plot.y1")))
  }
  # Compile the new plot
  plot.out = plot_grid(plotlist = plotlist.out, ncol = 1, align = "hv", 
                       labels = c("(a) resistance", "", "(b) recovery", "", "(c) resilience", ""), 
                       rel_heights = c(0.1, 1, 0.1, 1, 0.1, 1))
  plot.out = ggdraw(add_sub(
    plot.out, "Coordinate on the sgdd-wai PCA\n(Hot-dry to cold-wet climatic gradient)", 
    vpadding=grid::unit(0,"lines"),y=6, x=0.5, vjust=5.5))
  # --- END OF NEW APPROACH
  
  
  # Plot the estimates
  # plot.out = data.out %>%
  #   ggplot(aes(x = pca1, y = mean, group = clim, fill = clim)) + 
  #   geom_ribbon(aes(ymin = lwr, ymax = upr), color = NA, alpha = 0.5) +
  #   geom_line(color = "#001524") + 
  #   xlab("Coordinate on the sgdd-wai PCA\n(Hot-dry to cold-wet climatic gradient)") + 
  #   ylab("Effect on response to\ndisturbance metric") +
  #   facet_grid(effect ~ var.resp, scales = "free") +
  #   geom_hline(yintercept = 0, linetype = "dashed") +
  #   geom_text(aes(label = text), data = data.text, inherit.aes = TRUE, 
  #             hjust = "inward", size = 2.5, alpha = 0.8) +
  #   scale_fill_manual(values = colorRampPalette(c("orange", "blue"))(10)) +
  #   ylim(c(min(data.out$lwr), max(data.out$upr) + 0.2*diff(range(data.out$mean)))) +
  #   theme(panel.background = element_rect(color = "black", fill = "white"), 
  #         panel.grid = element_blank(), 
  #         strip.background = element_blank(), 
  #         strip.text = element_text(face = "bold"), 
  #         legend.position = "none") 
  
  
  # Name of the plot to save
  file.plot = paste0(dir.in, "/fd_effect_climate.jpg")
  file.plot.predictions = paste0(dir.in, "/observed_vs_predicted.jpg")
  file.plot.fitclim = paste0(dir.in, "/data_and_fit_clim.jpg")
  
  # Save plots
  ggsave(file.plot, plot.out, width = 28, height = 20 , units = "cm", 
         dpi = 600, bg = "white")
  ggsave(file.plot.predictions, plot.predictions, width = 14, height = 5, 
         units = "cm", dpi = 600, bg = "white")
  ggsave(file.plot.fitclim, plot.fit.clim, width = 14, height = 5, 
         units = "cm", dpi = 600, bg = "white")
  
  # Return name of the file saved
  return(c(file.plot, file.plot.predictions, file.plot.fitclim, files.stat.table))
}


#' Analyse the data with a structural equation model approach
#' @param data_model formatted model output
#' @param FD_metric Functional diversity metric to choose ("FDis", "FRic or "FD")
#' @param recovery_metric Recovery metric to choose ("recovery", "thalf)
#' @param R_metric Richness metric to choose ("nsp", "H" or "D")
#' @param dir.in name of the directory where to save outputs
plot_sem_multivar = function(data_model, FD_metric = "FDis", R_metric = "H", 
                             recovery_metric = "recovery", dir.in){
  
  # File for the figure
  fig.file.in = paste0(dir.in, "/sem_storm.jpg")
  table.file.in = paste0(dir.in, "/stat_sem.tex")
  
  # Create directory if needed
  create_dir_if_needed(fig.file.in)
  
  # -- Start by formatting data before fitting the model
  data.in = data_model %>%
    # Choose the right recovery metric
    rename("recov" = recovery_metric) %>%
    # Log transform resilience metrics to fit normality assumption
    mutate(resistance.log = resistance, 
           recovery.log = log(recov), 
           resilience.log = log(resilience)) %>%
    # Select the right FD metric
    rename("FD_chosen" = FD_metric, "R_chosen" = R_metric) %>%
    # scale all variables used in models
    mutate(climate_scaled = as.numeric(scale(pca1, center = TRUE, scale = TRUE)), 
           FD_scaled = as.numeric(scale(FD_chosen, center = FALSE, scale = TRUE)), 
           CWM1_scaled = as.numeric(scale(CWM1, center = TRUE, scale = TRUE)), 
           CWM2_scaled = as.numeric(scale(CWM2, center = TRUE, scale = TRUE)), 
           H_scaled = as.numeric(scale(R_chosen, center = FALSE, scale = TRUE)),
           resistance.log_scaled = as.numeric(scale(resistance.log, center = TRUE, scale = TRUE)), 
           recovery.log_scaled = as.numeric(scale(recovery.log, center = TRUE, scale = TRUE)), 
           resilience.log_scaled = as.numeric(scale(resilience.log, center = TRUE, scale = TRUE)))
  
  
  # -- Make model
  mod_sem = psem(
    glm(FD_scaled ~ climate_scaled + H_scaled, family = tweedie(var.power = 1), 
        data = data.in), 
    glm(H_scaled ~ climate_scaled, family = tweedie(var.power = 1), 
        data = data.in), 
    lm(CWM1_scaled ~ climate_scaled, data = data.in), 
    lm(CWM2_scaled ~ climate_scaled, data = data.in), 
    lm(resistance.log_scaled ~ FD_scaled + CWM1_scaled + CWM2_scaled + 
         H_scaled + climate_scaled, data = data.in), 
    lm(recovery.log_scaled ~ FD_scaled + CWM1_scaled + CWM2_scaled + 
         climate_scaled + H_scaled, data = data.in), 
    lm(resilience.log_scaled ~ resistance.log_scaled + recovery.log_scaled + 
         FD_scaled + CWM1_scaled + CWM2_scaled + climate_scaled + H_scaled, data = data.in)
  )
  
  
  
  
  
  # Plot the results
  
  # -- Prepare some parameters about the box to drax
  box.height = 2.5
  box.width = 5
  box.height.spacing = 2
  
  # -- Make a dataset to plot the box with the text
  data.plot.box = data.frame(
    text = c("climate", "FD", "CWM1", "CWM2", R_metric, "recovery", "resistance", "resilience"), 
    center.x = c(0, -7, 4.5, 12, -13, 2, -7, -1), 
    height.level = c(4, 3, 3, 3, 3, 2, 2, 1)) %>%
    mutate(ymin = (height.level - 0.5*box.height)*(box.height + box.height.spacing), 
           ymax = ymin + box.height,
           xmin = center.x - box.width/2, 
           xmax = center.x + box.width/2,
           center.y = 0.5*(ymin + ymax))
  
  
  # Loop on all models to extract output
  for(i in 1:(length(names(mod_sem)) - 1)){
    
    # model output for model i
    data.plot.arrow.i = data.frame(
      var.resp = as.character(summary(mod_sem[[i]])$terms[[2]]), 
      var.exp = as.character(rownames(summary(mod_sem[[i]])$coefficients)[-1]), 
      est = as.numeric(summary(mod_sem[[i]])$coefficients[-1, 1]), 
      se = as.numeric(summary(mod_sem[[i]])$coefficients[-1, 2]), 
      p = as.numeric(summary(mod_sem[[i]])$coefficients[-1, 4]), 
      class = ifelse("glm" %in% class(mod_sem[[i]]), "glm", "lm")
    ) %>%
      mutate(R2 = ifelse(
        class == "glm", 
        with(summary(mod_sem[[i]]), 1 - deviance/null.deviance), 
        summary(mod_sem[[i]])$r.squared))
    # Add to final output
    if(i == 1) data.plot.arrow = data.plot.arrow.i
    else data.plot.arrow = rbind(data.plot.arrow, data.plot.arrow.i)
  }
  
  # -- Prepare model outputs for plotting
  data.plot.arrow = data.plot.arrow %>%
    # simplify the name of variables
    mutate(var.exp = gsub("\\..+", "", gsub("\\_.+", "", var.exp)), 
           var.resp = gsub("\\..+", "", gsub("\\_.+", "", var.resp))) %>%
    # Add plot information on the resp variable
    left_join((data.plot.box %>%
                 dplyr::select(var.resp = text, level.resp = height.level, 
                               center.x.resp = center.x, center.y.resp = center.y)), 
              by = "var.resp") %>%
    # Add plot information on the exp variable
    left_join((data.plot.box %>%
                 dplyr::select(var.exp = text, level.exp = height.level, 
                               center.x.exp = center.x, center.y.exp = center.y)), 
              by = "var.exp") %>%
    # Add start and end of each arrow
    mutate(
      arrow.beg.y = case_when(
        (var.exp %in% c(R_metric, "CWM1", "CWM2") & var.resp == "resilience") ~ center.y.resp, 
        (var.exp == "H" & var.resp == "FD") ~ center.y.resp, 
        TRUE ~ (center.y.exp - 0.5*box.height)), 
      arrow.beg.x = ifelse((var.exp == "H" & var.resp == "FD"),
                           (center.x.exp + 0.5*box.width), center.x.exp), 
      arrow.end.y = case_when(
        (var.exp %in% c(R_metric, "CWM1", "CWM2") & var.resp == "resilience") ~ center.y.resp, 
        (var.exp == "H" & var.resp == "FD") ~ center.y.resp,
        TRUE ~ (center.y.resp + 0.5*box.height)), 
      arrow.end.x = case_when(
        (var.exp == R_metric & var.resp == "resilience") ~ (center.x.resp - 0.5*box.width), 
        (var.exp %in% c("CWM1", "CWM2") & var.resp == "resilience") ~ (center.x.resp + 0.5*box.width),
        (var.exp == "H" & var.resp == "FD") ~ (center.x.resp - 0.5*box.width),
        TRUE ~ center.x.resp)) %>%
    # Report significance, sign and absolute value of estimate
    mutate(significance = ifelse(p <= 0.05, "yes", "no"), 
           sign = ifelse(est < 0, "negative", "positive"), 
           est.abs = abs(est), 
           magnitude = case_when(
             est.abs <= quantile(est.abs, 0.33) ~ "low", 
             est.abs > quantile(est.abs, 0.33) & 
               est.abs <= quantile(est.abs, 0.66) ~ "mid", 
             TRUE ~ "high"
           )) 
  
  # Statistics table
  # -- Make the table
  table.stat = data.frame(
    col1 = c("Resp. var", data.plot.arrow$var.resp),
    col2 = c("Exp. var", data.plot.arrow$var.exp),
    col3 = c("Est", round(data.plot.arrow$est, digits = 2)),
    col4 = c("se", round(data.plot.arrow$se, digits = 2)),
    col5 = c("p", scales::pvalue(data.plot.arrow$p))
  )
  # -- Export the table
  print(xtable(table.stat, type = "latex", 
               caption = "Estimate, standard error and pvalue of the models included in the structural equation model", 
               label = "table_stat_sem"), 
        include.rownames=FALSE, hline.after = c(0, 1, dim(table.stat)[1]), 
        include.colnames = FALSE, caption.placement = "top", 
        file = table.file.in)
  
  # Add stat to the plot box data
  data.plot.box = data.plot.box %>%
    left_join((data.plot.arrow %>% 
                 dplyr::select(text = var.resp, class, R2) %>% 
                 mutate(R2 = round(R2, digits = 3)) %>% 
                 distinct() %>% 
                 mutate(text_stat = case_when(
                   text == "climate" ~ "Hot-dry to\ncold-wet", 
                   text == "CWM1" ~ paste0("CWM1\nsurv->growth\nR2 = ", R2),
                   text == "CWM2" ~ paste0("CWM2\nhigh->lowRec\nR2 = ", R2),
                   TRUE ~ paste0(text, "\n", "R2 = ", R2)
                 ))), 
              by = "text") %>%
    mutate(text_stat = ifelse(text == "climate", "Hot-dry to\ncold-wet", text_stat))
  
  # Add rectangles
  data.box.cat = data.plot.box %>%
    mutate(cat = case_when(
      text == "climate" ~ "Climate", 
      text %in% c(R_metric, "FD", "CWM1", "CWM2") ~ "Species\ncomposition", 
      text %in% c("resistance", "recovery", "resilience") ~ "Response to\ndisturbance"
    )) %>%
    group_by(cat) %>%
    summarise(xmin = min(xmin) - 0.05*diff(range(data.plot.box$center.x)), 
              xmax = max(xmax) + 0.05*diff(range(data.plot.box$center.x)), 
              ymin = min(ymin) - 0.05*diff(range(data.plot.box$center.y)), 
              ymax = max(ymax) + 0.05*diff(range(data.plot.box$center.y))) %>%
    mutate(center.x = xmax + 0.02*diff(range(data.plot.box$center.y)), 
           center.y = (ymin+ymax)/2)
  
  # -- Plot the boxes
  plot.box = data.plot.box %>%
    ggplot(aes(x = center.x, y = center.y))  + 
    geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = cat),
              data = data.box.cat,  color = "black", inherit.aes = TRUE, 
              show.legend = FALSE, alpha = 0.5)+ 
    geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), 
              color = "black", fill = "#2B2D42")+ 
    geom_text(aes(label = text_stat), color = "white", size = 4) + 
    geom_text(aes(label = cat), data = data.box.cat,
              hjust = 0, size = 4, inherit.aes = TRUE) +
    scale_fill_manual(values = c("#F28F3B", "#A9D6E5", "#38B000")) +
    xlim(min(data.box.cat$xmin) - 0.05*diff(range(data.plot.box$center.x)), 
         max(data.box.cat$center.x) + 0.3*diff(range(data.plot.box$center.x))) +
    theme(axis.title = element_blank(), 
          axis.ticks = element_blank(), 
          axis.text = element_blank(), 
          panel.background = element_rect(fill = "white", color = "black"), 
          panel.grid = element_blank(), 
          legend.key = element_blank())
  
  
  
  
  # -- plot the box with arrows
  plot.out = plot.box + 
    # Vertical segment
    geom_segment(data = subset(data.plot.arrow, var.exp %in% c(R_metric, "CWM1", "CWM2") & 
                                 var.resp == "resilience"), 
                 aes(x = center.x.exp, xend = center.x.exp, 
                     y = center.y.exp - 0.5*box.height, yend = center.y.resp, 
                     size = magnitude, linetype = significance, color = sign)) + 
    # Arrows
    geom_segment(data = data.plot.arrow, 
                 aes(x = arrow.beg.x, xend = arrow.end.x, 
                     y = arrow.beg.y, yend = arrow.end.y, 
                     size = magnitude, linetype = significance, color = sign), 
                 arrow = arrow(length = unit(0.15, "cm")), 
                 type = "closed") + 
    # -- scale linetype, color and size
    scale_color_manual(values = c(`negative` = "#0081A7", `positive` = "#F07167")) + 
    scale_linetype_manual(values = c(`yes` = "solid", `no` = "dashed")) + 
    scale_size_manual(values = c(`low` = 0.3, `mid` = 0.6, `high` = 1.1))
  
  
    # - Save the plot
  ggsave(fig.file.in, plot.out, width = 27, height = 14, units = "cm", dpi = 600, bg = "white")
  
  # return the name of all the plots made
  return(c(fig.file.in))
  
}

#' Analyse the data with a structural equation model approach
#' @param data_model formatted model output
#' @param FD_metric Functional diversity metric to choose ("FDis", "FRic or "FD")
#' @param recovery_metric Recovery metric to choose ("recovery", "thalf)
#' @param R_metric Richness metric to choose ("nsp", "H" or "D")
#' @param dir.in name of the directory where to save outputs
plot_sem_multivar_supp = function(data_model, FD_metric = "FDis", R_metric = "H",
                                  recovery_metric = "recovery", dir.in){
  
  # File for the figure
  fig.file.in = paste0(dir.in, "/sem_storm_supp.jpg")
  table.file.in = paste0(dir.in, "/stat_sem_supp.tex")
  
  # Create directory if needed
  create_dir_if_needed(fig.file.in)
  
  # -- Start by formatting data before fitting the model
  data.in = data_model %>%
    # Choose the right recovery metric
    rename("recov" = recovery_metric) %>%
    # Log transform resilience metrics to fit normality assumption
    mutate(resistance.log = resistance, 
           recovery.log = log(recov), 
           resilience.log = log(resilience)) %>%
    # Select the right FD metric
    rename("FD_chosen" = FD_metric, "R_chosen" = R_metric) %>%
    # scale all variables used in models
    mutate(climate_scaled = as.numeric(scale(pca1, center = TRUE, scale = TRUE)), 
           FD_scaled = as.numeric(scale(FD_chosen, center = FALSE, scale = TRUE)), 
           CWM1_scaled = as.numeric(scale(CWM1, center = TRUE, scale = TRUE)), 
           CWM2_scaled = as.numeric(scale(CWM2, center = TRUE, scale = TRUE)), 
           H_scaled = as.numeric(scale(R_chosen, center = FALSE, scale = TRUE)),
           resistance.log_scaled = as.numeric(scale(resistance.log, center = TRUE, scale = TRUE)), 
           recovery.log_scaled = as.numeric(scale(recovery.log, center = TRUE, scale = TRUE)), 
           resilience.log_scaled = as.numeric(scale(resilience.log, center = TRUE, scale = TRUE)))
  
  
  # -- Make model
  mod_sem = psem(
    glm(FD_scaled ~ climate_scaled + H_scaled, family = tweedie(var.power = 1), 
        data = data.in), 
    glm(H_scaled ~ climate_scaled, family = tweedie(var.power = 1), 
        data = data.in), 
    lm(CWM1_scaled ~ climate_scaled + FD_scaled + H_scaled, data = data.in), 
    lm(CWM2_scaled ~ climate_scaled + FD_scaled + H_scaled + CWM1_scaled, data = data.in), 
    lm(resistance.log_scaled ~ FD_scaled + CWM1_scaled + CWM2_scaled + 
         H_scaled + climate_scaled, data = data.in), 
    lm(recovery.log_scaled ~ FD_scaled + CWM1_scaled + CWM2_scaled + 
         climate_scaled + H_scaled + resistance.log_scaled, data = data.in), 
    lm(resilience.log_scaled ~ resistance.log_scaled + recovery.log_scaled + 
         FD_scaled + CWM1_scaled + CWM2_scaled + climate_scaled + H_scaled, data = data.in)
  )
  
  
  
  
  
  # Plot the results
  
  # -- Prepare some parameters about the box to drax
  box.height = 2.5
  box.width = 5
  box.height.spacing = 2
  
  # -- Make a dataset to plot the box with the text
  data.plot.box = data.frame(
    text = c("climate", "FD", "CWM1", "CWM2", R_metric, "recovery", "resistance", "resilience"), 
    center.x = c(0, -7, 4.5, 12, -13, 2, -7, -1), 
    height.level = c(4, 3, 3, 3, 3, 2, 2, 1)) %>%
    mutate(ymin = (height.level - 0.5*box.height)*(box.height + box.height.spacing), 
           ymax = ymin + box.height,
           xmin = center.x - box.width/2, 
           xmax = center.x + box.width/2,
           center.y = 0.5*(ymin + ymax))
  
  
  # Loop on all models to extract output
  for(i in 1:(length(names(mod_sem)) - 1)){
    
    # model output for model i
    data.plot.arrow.i = data.frame(
      var.resp = as.character(summary(mod_sem[[i]])$terms[[2]]), 
      var.exp = as.character(rownames(summary(mod_sem[[i]])$coefficients)[-1]), 
      est = as.numeric(summary(mod_sem[[i]])$coefficients[-1, 1]), 
      se = as.numeric(summary(mod_sem[[i]])$coefficients[-1, 2]), 
      p = as.numeric(summary(mod_sem[[i]])$coefficients[-1, 4]), 
      class = ifelse("glm" %in% class(mod_sem[[i]]), "glm", "lm")
    ) %>%
      mutate(R2 = ifelse(
        class == "glm", 
        with(summary(mod_sem[[i]]), 1 - deviance/null.deviance), 
        summary(mod_sem[[i]])$r.squared))
    # Add to final output
    if(i == 1) data.plot.arrow = data.plot.arrow.i
    else data.plot.arrow = rbind(data.plot.arrow, data.plot.arrow.i)
  }
  
  # -- Prepare model outputs for plotting
  data.plot.arrow = data.plot.arrow %>%
    # simplify the name of variables
    mutate(var.exp = gsub("\\..+", "", gsub("\\_.+", "", var.exp)), 
           var.resp = gsub("\\..+", "", gsub("\\_.+", "", var.resp))) %>%
    # Add plot information on the resp variable
    left_join((data.plot.box %>%
                 dplyr::select(var.resp = text, level.resp = height.level, 
                               center.x.resp = center.x, center.y.resp = center.y)), 
              by = "var.resp") %>%
    # Add plot information on the exp variable
    left_join((data.plot.box %>%
                 dplyr::select(var.exp = text, level.exp = height.level, 
                               center.x.exp = center.x, center.y.exp = center.y)), 
              by = "var.exp") %>%
    # Add start and end of each arrow
    mutate(
      arrow.beg.y = case_when(
        (var.exp %in% c(R_metric, "CWM1", "CWM2") & var.resp == "resilience") ~ center.y.resp, 
        (var.exp == "H" & var.resp == "FD") ~ center.y.resp, 
        TRUE ~ (center.y.exp - 0.5*box.height)), 
      arrow.beg.x = ifelse((var.exp == "H" & var.resp == "FD"),
                           (center.x.exp + 0.5*box.width), center.x.exp), 
      arrow.end.y = case_when(
        (var.exp %in% c(R_metric, "CWM1", "CWM2") & var.resp == "resilience") ~ center.y.resp, 
        (var.exp == "H" & var.resp == "FD") ~ center.y.resp,
        TRUE ~ (center.y.resp + 0.5*box.height)), 
      arrow.end.x = case_when(
        (var.exp == R_metric & var.resp == "resilience") ~ (center.x.resp - 0.5*box.width), 
        (var.exp %in% c("CWM1", "CWM2") & var.resp == "resilience") ~ (center.x.resp + 0.5*box.width),
        (var.exp == "H" & var.resp == "FD") ~ (center.x.resp - 0.5*box.width),
        TRUE ~ center.x.resp)) %>%
    # Report significance, sign and absolute value of estimate
    mutate(significance = ifelse(p <= 0.05, "yes", "no"), 
           sign = ifelse(est < 0, "negative", "positive"), 
           est.abs = abs(est), 
           magnitude = case_when(
             est.abs <= quantile(est.abs, 0.33) ~ "low", 
             est.abs > quantile(est.abs, 0.33) & 
               est.abs <= quantile(est.abs, 0.66) ~ "mid", 
             TRUE ~ "high"
           )) 
  
  # Statistics table
  # -- Make the table
  table.stat = data.frame(
    col1 = c("Resp. var", data.plot.arrow$var.resp),
    col2 = c("Exp. var", data.plot.arrow$var.exp),
    col3 = c("Est", round(data.plot.arrow$est, digits = 2)),
    col4 = c("se", round(data.plot.arrow$se, digits = 2)),
    col5 = c("p", scales::pvalue(data.plot.arrow$p))
  )
  # -- Export the table
  print(xtable(table.stat, type = "latex", 
               caption = "Estimate, standard error and pvalue of the models included in the structural equation model", 
               label = "table_stat_sem_supp"), 
        include.rownames=FALSE, hline.after = c(0, 1, dim(table.stat)[1]), 
        include.colnames = FALSE, caption.placement = "top", 
        file = table.file.in)
  
  # Add stat to the plot box data
  data.plot.box = data.plot.box %>%
    left_join((data.plot.arrow %>% 
                 dplyr::select(text = var.resp, class, R2) %>% 
                 mutate(R2 = round(R2, digits = 3)) %>% 
                 distinct() %>% 
                 mutate(text_stat = case_when(
                   text == "climate" ~ "Hot-dry to\ncold-wet", 
                   text == "CWM1" ~ paste0("CWM1\nsurv->growth\nR2 = ", R2),
                   text == "CWM2" ~ paste0("CWM2\nhigh->lowRec\nR2 = ", R2),
                   TRUE ~ paste0(text, "\n", "R2 = ", R2)
                 ))), 
              by = "text") %>%
    mutate(text_stat = ifelse(text == "climate", "Hot-dry to\ncold-wet", text_stat))
  
  # Add rectangles
  data.box.cat = data.plot.box %>%
    mutate(cat = case_when(
      text == "climate" ~ "Climate", 
      text %in% c(R_metric, "FD", "CWM1", "CWM2") ~ "Species\ncomposition", 
      text %in% c("resistance", "recovery", "resilience") ~ "Response to\ndisturbance"
    )) %>%
    group_by(cat) %>%
    summarise(xmin = min(xmin) - 0.05*diff(range(data.plot.box$center.x)), 
              xmax = max(xmax) + 0.05*diff(range(data.plot.box$center.x)), 
              ymin = min(ymin) - 0.05*diff(range(data.plot.box$center.y)), 
              ymax = max(ymax) + 0.05*diff(range(data.plot.box$center.y))) %>%
    mutate(center.x = xmax + 0.02*diff(range(data.plot.box$center.y)), 
           center.y = (ymin+ymax)/2)
  
  # -- Adjust data plot arrow to remove relations that are not of interest
  data.plot.arrow <- data.plot.arrow %>%
    filter(!(var.resp == "CWM1" & var.exp %in% c("H", "FD"))) %>%
    filter(!(var.resp == "CWM2" & var.exp %in% c("H", "FD", "CWM1"))) %>%
    filter(!(var.resp == "recovery" & var.exp == "resistance")) 
  
  # -- Plot the boxes
  plot.box = data.plot.box %>%
    ggplot(aes(x = center.x, y = center.y))  + 
    geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = cat),
              data = data.box.cat,  color = "black", inherit.aes = TRUE, 
              show.legend = FALSE, alpha = 0.5)+ 
    geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), 
              color = "black", fill = "#2B2D42")+ 
    geom_text(aes(label = text_stat), color = "white", size = 4) + 
    geom_text(aes(label = cat), data = data.box.cat,
              hjust = 0, size = 4, inherit.aes = TRUE) +
    scale_fill_manual(values = c("#F28F3B", "#A9D6E5", "#38B000")) +
    xlim(min(data.box.cat$xmin) - 0.05*diff(range(data.plot.box$center.x)), 
         max(data.box.cat$center.x) + 0.3*diff(range(data.plot.box$center.x))) +
    theme(axis.title = element_blank(), 
          axis.ticks = element_blank(), 
          axis.text = element_blank(), 
          panel.background = element_rect(fill = "white", color = "black"), 
          panel.grid = element_blank(), 
          legend.key = element_blank())
  
  
  
  
  # -- plot the box with arrows
  plot.out = plot.box + 
    # Vertical segment
    geom_segment(data = subset(data.plot.arrow, var.exp %in% c(R_metric, "CWM1", "CWM2") & 
                                 var.resp == "resilience"), 
                 aes(x = center.x.exp, xend = center.x.exp, 
                     y = center.y.exp - 0.5*box.height, yend = center.y.resp, 
                     size = magnitude, linetype = significance, color = sign)) + 
    # Arrows
    geom_segment(data = data.plot.arrow, 
                 aes(x = arrow.beg.x, xend = arrow.end.x, 
                     y = arrow.beg.y, yend = arrow.end.y, 
                     size = magnitude, linetype = significance, color = sign), 
                 arrow = arrow(length = unit(0.15, "cm")), 
                 type = "closed") + 
    # -- scale linetype, color and size
    scale_color_manual(values = c(`negative` = "#0081A7", `positive` = "#F07167")) + 
    scale_linetype_manual(values = c(`yes` = "solid", `no` = "dashed")) + 
    scale_size_manual(values = c(`low` = 0.3, `mid` = 0.6, `high` = 1.1))
  
  
  # - Save the plot
  ggsave(fig.file.in, plot.out, width = 27, height = 14, units = "cm", dpi = 600, bg = "white")
  
  # return the name of all the plots made
  return(c(fig.file.in))
  
}





#' Export species information in a table
#' @param FUNDIV_data Pre-formatted data from FUNDIV
#' @param pca12_per_species Coordinates of each species on the trait PCA
#' @param species_list list of species for which IPM is to be computed
#' @param file.out Character: name including path of the file to export (end by .tex)
export_species_table = function(FUNDIV_data, pca12_per_species, species_list, file.out){
  
  # Create output directory if needed
  create_dir_if_needed(file.out)
  
  # Get proportion per species
  species_table = FUNDIV_data %>%
    # Remove in-growth trees
    filter(treestatus_th != 1) %>%
    # Calculate basal area per species
    group_by(species) %>%
    summarize(ba_ha = sum(ba_ha1, na.rm = TRUE), 
              dbh_mean = weighted.mean(dbh1, w = weight1, na.rm = TRUE)) %>%
    # Convert in proportion across dataset
    ungroup() %>%
    mutate(ba_percent = round(ba_ha/sum(ba_ha, na.rm = TRUE), digits = 3)*100) %>%
    # Keep only species in the simulations
    mutate(species = gsub("\\ ", "\\_", species)) %>%
    filter(species %in% unique(species_list$species)) %>%
    # Add climate niche
    # -- cold margin
    left_join((climate_species %>%
                 filter(N == 3) %>%
                 mutate(cold = paste0(round(wai, 2), "/", round(sgdd, 0))) %>%
                 dplyr::select(species = sp, cold)), by = "species") %>%
    # -- Optimum
    left_join((climate_species %>%
                 filter(N == 2) %>%
                 mutate(optimum = paste0(round(wai, 2), "/", round(sgdd, 0))) %>%
                 dplyr::select(species = sp, optimum)), by = "species") %>%
    # -- hot margin
    left_join((climate_species %>%
                 filter(N == 1) %>%
                 mutate(hot = paste0(round(wai, 2), "/", round(sgdd, 0))) %>%
                 dplyr::select(species = sp, hot)), by = "species") %>%
    # Add traits value
    left_join(pca12_per_species, by = "species") %>%
    # Final formatting
    mutate(species = gsub("\\_", "\\ ", species), 
           dbh_mean = as.character(round(dbh_mean, digits = 0)), 
           ba_percent = as.character(ba_percent), 
           pca1 = as.character(round(pca1, digits = 2)), 
           pca2 = as.character(round(pca2, digits = 2))) %>%
    dplyr::select(-ba_ha)
  
  # Finalize the table for latex
  tex.table = data.frame(
    species = "", dbh_mean = c("mean dbh", "(mm)"), ba_percent = c("share in ", "basal area (%)"), 
    cold = c("", "cold mar."), optimum = c("wai/sgdd", "optimum"), hot = c("", "hot mar."), 
    pca1 = c("trait", "PC1"), pca2 = c("value", "PC2")) %>%
    rbind(species_table) 
  
  # Save tex file
  print(xtable(tex.table, type = "latex", label = "table_species",
               caption = paste0("Mean dbh, share in basal area and climatic niche (cold ", 
                                "margin, optimum, hot margin) within the NFI dataset", 
                                " and trait value of the 15 species included in the simulations"), 
               align = c("c", "c", "c", "c", "r", "c", "l", "r", "l")), 
        include.rownames=FALSE, hline.after = c(0, 2, dim(tex.table)[1]), 
        include.colnames = FALSE, caption.placement = "top", size = "\\scriptsize",
        file = file.out)
  
  # Return the file exported
  return(file.out)
  
  
}