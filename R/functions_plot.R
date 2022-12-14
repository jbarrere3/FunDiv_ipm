#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#### SCRIPT INTRODUCTION ####
#
#' @name functions_plot.R  
#' @description R script containing all functions relative to data
#               importation and formatting
#' @author Julien Barrere, Marianne Bernard
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



#' Function to plot the difference between probability calculated with mean 
#' parameter values or with all iterations plus uncertainty
#' @param species.in vector of species to analyze
#' @param disturbance.time time at which disturbance will be applied
#' @param disturbance.in which disturbance to apply
#' @param intensity.in intensity (between 0 and 1) of the disturbance
#' @param duration.disturbance how long should the disturbance kill trees 
#' @param disturbance_parameters Parameters averaged over all iterations
#' @param disturbance_parameters_alliter Parameters non averaged
#' @param file.in name including path of the file to save
plot_meanProba_vs_meanParam <- function(
  species.in, disturbance.time = 1500, disturbance.in = "storm", intensity.in = 0.8, IPM.list,
  duration.disturbance = 5, disturbance_parameters, disturbance_parameters_alliter, file.in){
  
  # Create directory if needed
  create_dir_if_needed(file.in)
  
  # List of species objects
  species.list.init <- generate_species.list(IPM.list[species.in])
  
  # Loop on all species
  for(i in 1:length(species.in)){
    
    # Species object for species i
    species.init <- species.list.init[[species.in[i]]]
    
    # Make simulation for species i
    sim <- sim_deter_forest(
      species.init, tlim = disturbance.time, equil_time = disturbance.time, equil_dist = 1,
      harvest = "default", targetBA = 40, SurfEch = 0.03, correction = "cut", verbose = TRUE
    )
    
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
      mutate(p.mean.param = (1 - plogis(a0 + a1*logratio.scaled + b*I^(c*dbh.scaled)))^t) %>%
      dplyr::select(size, value, p.mean.param)
    
    
    
    # Disturb population with all iterations
    disturbed.proba.alliter <- pop.to.disturb %>%
      # Add disturbance parameters
      mutate(species = gsub("\\_", "\\ ", species), 
             disturbance = disturbance.in) %>%
      merge(expand.grid(Chain = unique(disturbance_parameters_alliter$Chain), 
                        Iteration = unique(disturbance_parameters_alliter$Iteration))) %>%
      left_join(disturbance_parameters_alliter, by = c("species", "disturbance", "Chain", "Iteration")) %>%
      # Calculate variables used to calculate probability
      mutate(dqm = dqm.disturbance.time, 
             logratio = log(size/dqm), 
             dbh.scaled = dbh.intercept + size*dbh.slope, 
             logratio.scaled = logratio.intercept + logratio*logratio.slope, 
             I = intensity.in, 
             t = duration.disturbance) %>%
      # Calculate probability of survival from disturbance per size class
      mutate(p = (1 - plogis(a0 + a1*logratio.scaled + b*I^(c*dbh.scaled)))^t) %>%
      # Summarize over each size class
      group_by(size) %>%
      summarize(p.alliter.mean = mean(p), 
                p.alliter.inf = quantile(p, probs = 0.025), 
                p.alliter.sup = quantile(p, probs = 0.975)) 
    
    # Gather datasets together
    data.i <- disturbed.proba %>%
      left_join(disturbed.proba.alliter, by = "size") %>%
      rename(undisturbed = value, meanParam = p.mean.param, meanProba = p.alliter.mean) %>%
      mutate(meanParam = meanParam*undisturbed, meanProba = meanProba*undisturbed, 
             p.alliter.inf = p.alliter.inf*undisturbed, p.alliter.sup = p.alliter.sup*undisturbed) %>%
      gather(key = "method", value = "Ntrees", "undisturbed", "meanParam", "meanProba") %>%
      mutate(species = gsub("\\_", "\\ ", species.in[i]))
    
    # Add to the final dataset
    if(i == 1) data <- data.i
    else data <- rbind(data, data.i)
    
  }
  
  # Make the plot
  plot.out <- data %>%
    ggplot(aes(x = size, y = Ntrees, group = method, color = method)) + 
    geom_ribbon(aes(ymin = p.alliter.inf, ymax = p.alliter.sup), 
                fill = "#E9ECEF", color = NA) +
    geom_line() + 
    scale_color_manual(values = c("#15616D", "#FF7D00", "#001524")) +
    facet_wrap(~ species) +
    theme(panel.background = element_rect(fill = "white", color = "black"), 
          panel.grid = element_blank(), 
          legend.title = element_blank(), 
          legend.key = element_blank(), 
          strip.background = element_blank())
  
  # - Save the plot
  ggsave(file.in, plot.out, width = 20, height = 12, units = "cm", dpi = 600, bg = "white")
  
  # return the name of all the plots made
  return(file.in)
  
}


#' Function to plot disturbance effect on a monospecific stand
#' @param species.list List of species objects
#' @param disturbance_parameters Parameters for the disturbance mortality
#' @param disturbance.in character: disturbance to apply (e.g. "storm")
#' @param intensity.in numeric: intensity of the disturbance btw 0 and 1
#' @param duration.disturbance numeric: how long should the disturbance effect last (in years)
#' @param sim.time How long should the simulation last (disturbance occur in the middle)
#' @param SurfEch.in Area in ha to sample
#' @param file.in Name of the file to save, including path
plot.monosp.disturbances <- function(species.list, disturbance_parameters, disturbance.in, 
                                     intensity.in, duration.disturbance, sim.time, 
                                     SurfEch.in = 0.03, file.in){
  
  # Create directory for file if needed
  create_dir_if_needed(file.in)
  
  # Timing of the disturbance
  dist.time <- round(sim.time/2, digits = 0)
  
  # function to initialize basal area
  default_init <- function(x){
    force(x)
    function(mesh, SurfEch = SurfEch.in) {
      return(x)
    }
  }
  
  
  # Loop on all species in the list
  for(i in 1:length(names(species.list))){
    
    # Species printer
    print(paste0("Starting simulations for ", 
                 gsub("\\_", "\\ ", names(species.list)[i])))
    
    # Simulation printer
    print("---- Simulation before disturbance")
    
    # Run simulation for the first species before disturbance
    sim.i.before <- sim_deter_forest(
      species.list[[i]], tlim = dist.time, equil_time = dist.time, equil_dist = 1,
      harvest = "default", targetBA = 40, SurfEch = SurfEch.in, verbose = FALSE
    )
    
    # Disturb the population
    dist.pop.i <- disturb.population(sim.i.before, disturbance.in, intensity.in, 
                                     duration.disturbance, disturbance_parameters)
    
    # Create species object disturbed
    species.i.disturbed <- species.list[[i]]
    species.i.disturbed$init_pop <- default_init(dist.pop.i*SurfEch.in)
    
    # Simulation printer
    print("---- Simulation after disturbance")
    
    # Run simulation for the first species after disturbance
    sim.i.after <- sim_deter_forest(
      species.i.disturbed, tlim = (sim.time - dist.time), equil_time = (sim.time - dist.time), 
      equil_dist = 1, harvest = "default", targetBA = 40, SurfEch = SurfEch.in, verbose = FALSE
    )
    
    # Data formatting
    data.i <- tree_format(sim.i.after) %>%
      mutate(time = time + dist.time) %>%
      rbind(tree_format(sim.i.before)) %>%
      filter(var == "BAsp") %>%
      dplyr::select(-equil) %>%
      distinct()
    
    # Add to final dataset
    if(i == 1) data <- data.i
    else data <- rbind(data, data.i)
  }
  
  # Make the plot
  plot.out <- data %>%
    mutate(species = gsub("\\_", "\\ ", species)) %>%
    ggplot(aes(x = time, y = value, group = 1)) + 
    geom_line() + 
    facet_wrap(~ species) + 
    geom_vline(xintercept = dist.time, color = "red", linetype = "dashed") + 
    theme(panel.background = element_rect(fill = "white", color = "black"), 
          panel.grid = element_blank(), 
          strip.background = element_blank())
  
  # - Save the plot
  ggsave(file.in, plot.out, width = 20, height = 12, units = "cm", dpi = 600, bg = "white")
  
  # return the name of all the plots made
  return(file.in)
  
  
  
}





#' Plot a list of simulations
#' @param sim.list list of simulations
#' @param file.in Name including path of the file to save
plot_sim.list <- function(sim.list, file.in){
  
  # Create directory if needed
  create_dir_if_needed(file.in)
  
  # Loop on all species combinations
  for(i in 1:length(names(sim.list))){
    # Vector of species for simulation i
    sp.vec.i <- unlist(strsplit(names(sim.list)[i], "\\."))
    
    # Make a nice name from species vector
    name.i <- paste(paste0(substr(sp.vec.i, 1, 2), ". ",
                           substr(gsub(".+\\_", "", sp.vec.i), 1, 2)), 
                    collapse = " - ")
    
    # Format data for plotting simulation i
    data.i <- tree_format(sim.list[[i]]) %>%
      # Remove unused lines and columns
      filter(var == "BAsp") %>%
      dplyr::select(species, time, value) %>%
      distinct() %>%
      # Calculate the sum of basal area
      spread(key = "species", value = "value") %>%
      mutate(all = rowSums(.[, c(2:dim(.)[2])])) %>%
      gather(key = "species", value = "value", c(sp.vec.i, "all")) %>%
      # Add species combination
      mutate(sp.combination = name.i)
    
    # Add to final dataset
    if(i == 1) data = data.i
    else data <- rbind(data, data.i)
  }
  
  # Make the plot 
  plot.out <- data %>%
    mutate(species = factor(
      species, levels = c("all", unique(unlist(strsplit(names(sim.list), "\\.")))))) %>%
    ggplot(aes(x = time, y = value, group = species, color = species)) + 
    scale_color_manual(values = c("#D90429", "#1E6091", "#168AAD", "#34A0A4", "#76C893", 
                                  "#B5E48C", "#00AFB9")[c(1:length(unique(data$species)))]) +
    geom_line() + 
    facet_wrap(~ sp.combination) + 
    theme(panel.background = element_rect(fill = "white", color = "black"), 
          panel.grid = element_blank(), 
          strip.background = element_blank(), 
          legend.title = element_blank(), 
          legend.key = element_blank()) + 
    xlab("Time (years)") + ylab("Basal area")
  
  # - Save the plot
  ggsave(file.in, plot.out, width = 30, height = 20, units = "cm", dpi = 600, bg = "white")
  
  # return the name of all the plots made
  return(file.in)
}




