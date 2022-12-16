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




