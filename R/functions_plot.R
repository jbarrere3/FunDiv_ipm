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


#' Plot a list of simulations with multiple climates
#' @param sim.list_multiclim list of simulations with multiple climates
#' @param file.in Name including path of the file to save
plot_sim.list_multiclim <- function(sim.list, file.in){
  
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
    
    # Loop on all climates
    for(j in 1:length(names(sim.list[[i]]))){
      
      # Format data for plotting forest i and climate j
      data.ij <- tree_format(sim.list[[i]][[j]]) %>%
        # Remove unused lines and columns
        filter(var == "BAsp") %>%
        dplyr::select(species, time, value) %>%
        distinct() %>%
        # Calculate the sum of basal area
        spread(key = "species", value = "value") %>%
        mutate(all = rowSums(.[, c(2:dim(.)[2])])) %>%
        gather(key = "species", value = "value", c(sp.vec.i, "all")) %>%
        # Add species combination
        mutate(sp.combination = name.i) %>%
        # Add climate
        mutate(climate = names(sim.list[[i]])[j])
      
      # Add to final dataset
      if(i == 1 & j == 1) data = data.ij
      else data <- rbind(data, data.ij)
    }
    
  }
  
  # Make the plot 
  plot.out <- data %>%
    filter(species == "all") %>%
    mutate(climate = factor(climate, levels = c("cold", "optimum", "hot"))) %>%
    ggplot(aes(x = time, y = value, group = climate, color = climate)) + 
    scale_color_manual(values = c("#0466C8", "#02040F", "#C20A16")) +
    geom_line() + 
    facet_wrap(~ sp.combination) + 
    theme(panel.background = element_rect(fill = "white", color = "black"), 
          panel.grid = element_blank(), 
          strip.background = element_blank(), 
          legend.title = element_blank(), 
          legend.key = element_blank()) + 
    xlab("Time (years)") + ylab("Basal area")
  
  # - Save the plot
  ggsave(file.in, plot.out, width = 20, height = 12, units = "cm", dpi = 600, bg = "white")
  
  # return the name of all the plots made
  return(file.in)
}


#' Plot traits PCA
#' @param traits dataframe containing trait value per species
#' @param file.in name including path of the file to save
plot_traits_pca <- function(traits, file.in){
  
  # Create directory if needed
  create_dir_if_needed(file.in)
  
  # - Make PCA 
  pca <- prcomp((traits %>% dplyr::select(-species)), 
                center = T, scale = T)
  
  # - Extract the coordinates of the ndividuals on pca axis
  res.ind <- data.frame(species = traits$species, 
                        pca1 = get_pca_ind(pca)[[1]][, 1], 
                        pca2 = get_pca_ind(pca)[[1]][, 2]) 
  
  # - Extract the coordinates of the variables on pca axis
  res.var <- data.frame(var = rownames(get_pca_var(pca)[[1]]), 
                        pca1 = get_pca_var(pca)[[1]][, 1], 
                        pca2 = get_pca_var(pca)[[1]][, 2]) %>%
    mutate(var = gsub("\\.", "\\ ", var))
  
  # - Minimum and maximum in each pca axis
  pca.xmin <- -max(abs(res.ind$pca1))
  pca.xmax <- max(abs(res.ind$pca1))
  pca.ymin <- -max(abs(res.ind$pca2))
  pca.ymax <- max(abs(res.ind$pca2))
  
  # Make the plot
  plot.out <- res.ind %>%
    ggplot(aes(x = pca1, y = pca2)) + 
    geom_point(fill = "#023E8A", color = "black", shape = 21, size = 3) +
    geom_segment(data = (res.var %>% mutate(pca1 = pca1*1.5, pca2 = pca2*1.5)), 
                 aes(x = 0, xend = pca1, y = 0, yend = pca2), 
                 arrow = arrow(length = unit(0.1, "cm")), 
                 type = "closed", color = "#D90429") + 
    geom_text(data = (res.var %>% mutate(pca1 = pca1*2, pca2 = pca2*2)), 
              aes(label = var), color = "#D90429", size = 5, 
              nudge_x = ifelse(res.var$pca1 < 0, pca.xmin/12, pca.xmax/12)) +
    geom_hline(size = 0.2, yintercept = 0, color = "#6C757D", linetype = "dashed") + 
    geom_vline(size = 0.2, xintercept = 0, color = "#6C757D", linetype = "dashed") + 
    xlim((pca.xmin-0.2), (pca.xmax+0.2)) + 
    ylim((pca.ymin-0.2), (pca.ymax+0.2)) +
    xlab(paste0("PCA1 (", round(summary(pca)$importance[2, 1]*100, digits = 2), "%)")) +
    ylab(paste0("PCA2 (", round(summary(pca)$importance[2, 2]*100, digits = 2), "%)")) +
    theme(panel.background = element_rect(color = "black", fill = "white"), 
          panel.grid = element_blank(), 
          axis.title = element_text(size = 15))
  
  # - Save the plot
  ggsave(file.in, plot.out, width = 14, height = 14, units = "cm", dpi = 600, bg = "white")
  
  # return the name of all the plots made
  return(file.in)
}



#' Plot the effect of functional diversity on resilience
#' @param FD_and_resilience Dataset contaiing fd and resilience per species
#' @param file.in Name including path of the file to save
plot_FD_vs_resilience <- function(FD_and_resilience, file.in){
  
  # Create directory if needed
  create_dir_if_needed(file.in)
  
  # Fit model
  mod <- lm(log(Resilience) ~ FD*CWM, data = FD_and_resilience)
  
  # Data with model predictions for plotting
  data.fit <- expand.grid(FD = seq(from = min(FD_and_resilience$FD), 
                                   to = max(FD_and_resilience$FD), 
                                   length.out = 30), 
                          CWM = quantile(FD_and_resilience$CWM, 
                                         probs = c(0.1, 0.5, 0.9)))
  data.fit$Resilience <- exp(predict(mod, newdata = data.fit))
  
  # Make the plot
  plot.out <- FD_and_resilience %>%
    ggplot(aes(x = FD, y = Resilience, color = CWM, group = CWM)) +
    geom_point() + 
    geom_line(data = data.fit, inherit.aes = TRUE) +
    theme(panel.background = element_rect(color = "black", fill = "white"), 
          panel.grid = element_blank(), 
          axis.title = element_text(size = 15))
  
  # Save the plot
  ggsave(file.in, plot.out, width = 17, height = 14, units = "cm", dpi = 600, bg = "white")
  
  # return the name of all the plots made
  return(file.in)
}
