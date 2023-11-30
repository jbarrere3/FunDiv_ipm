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
  
  # Also export as csv file
  # -- Change filename
  file.csv.out = gsub("\\.tex", "\\.csv", file.out)
  write.table(tex.table, file = file.csv.out, row.names = FALSE, 
              col.names = FALSE, sep = ",")
  
  # Return the file exported
  return(c(file.out, file.csv.out))
  
  
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




#' Extract for all simulations the basal area of species present at a given timestep
#' @param sim_disturbance character vector containing name of rds files  
#'                        containing disturbance simulations
#' @param sp.in.sim character vector of all species present in the simulations
#' @param times.in vector of time steps for which to extract sp compo
get_sp_composition = function(sim_disturbance, sp.in.sim, times.in){
  
  # Initialize output dataframe for each time step
  df.out = data.frame(ID.forest = c(1:length(sim_disturbance)))
  for(sp in sp.in.sim) eval(parse(text = paste0("df.out$", sp, " = 0")))
  
  # Initialize output list
  out = list()
  for(t in 1:length(times.in)) eval(parse(text = paste0(
    "out$time", times.in[t], " = df.out")))
  
  # Loop on all species combination
  for(i in 1:length(sim_disturbance)){
    
    # Printer
    print(paste0("Reading simulation ", i, "/", length(sim_disturbance)))
    
    # Read simulation i
    sim.i = readRDS(sim_disturbance[i])
    
    # First, verify that equilibrium was reached
    if(!is.na(sim.i[1, 1])){
      
      # If equilibrium is reached, then extract species composition at each time
      for(k in 1:length(times.in)){
        
        # Restrict sim output to the year of interest
        data.ik = sim.i %>% 
          filter(time == times.in[k] & var == "BAsp") %>%
          filter(!equil)
        
        # Loop on species present to complete output
        for(j in 1:dim(data.ik)[1]) out[[k]][i, data.ik$species[j]] = data.ik$value[j]
      }
      
    }
    
  }
  
  # Return output
  return(out)
  
}






