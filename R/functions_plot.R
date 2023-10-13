#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#### SCRIPT INTRODUCTION ####
#
#' @name functions_plot.R  
#' @description R script containing all functions relative to data
#               importation and formatting
#' @author Julien Barrere
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# -- Informative plots -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



#' Plot traits PCA
#' @param traits dataframe containing trait value per species
#' @param species.in vector of the species present in the simulations
#' @param file.in name including path of the file to save
plot_traits_pca <- function(traits, species.in, file.in){
  
  # Create directory if needed
  create_dir_if_needed(file.in)
  
  # - Load the list of species included in IPM
  data("fit_species")
  
  # - Make PCA 
  pca <- prcomp((traits %>% dplyr::select(-species)), 
                center = T, scale = T)
  
  # - Extract the coordinates of the variables on pca axis
  res.var <- data.frame(var = rownames(get_pca_var(pca)[[1]]), 
                        pca1 = get_pca_var(pca)[[1]][, 1]) %>%
    mutate(var = gsub("\\.", "\\ ", var)) %>%
    mutate(var.pos = c(1:dim(.)[1]))
  
  # - Extract the coordinates of the ndividuals on pca axis
  res.ind <- data.frame(species = traits$species, 
                        pca1 = get_pca_ind(pca)[[1]][, 1], 
                        var.pos = max(res.var$var.pos) + 1) %>%
    filter(species %in% fit_species) %>%
    mutate(species = gsub("\\_", "\\ ", species))  %>%
    arrange(pca1) %>%
    mutate(pca1_seq = seq(from = min(.$pca1), to = max(.$pca1), 
                          length.out = dim(.)[1]))
  
  # -Extend arrow length to match max of individuals
  res.var$pca1 = res.var$pca1*(max(abs(res.ind$pca1))/max(abs(res.var$pca1)))*0.9
  
  # - Vector of the pca axis scale 
  vec.pca.axis = eval(parse(text = paste0(
    "c(", paste(round(range(res.ind$pca1), digits = 0), collapse = ":"), ")")))
  
  # - Remove species not present in the simulations
  res.ind = res.ind %>% filter(species %in% gsub("\\_", "\\ ", species.in))
  
  # Make the plot
  plot.out <- res.var %>%
    ggplot(aes(x = var.pos, xend = var.pos, y = 0, yend = pca1)) + 
    geom_segment(arrow = arrow(length = unit(0.1, "cm"))) + 
    scale_x_continuous(breaks = res.var$var.pos, 
                       labels = res.var$var, 
                       limits = c(0.5, max(res.ind$var.pos)+9))+ 
    scale_y_continuous(breaks = vec.pca.axis)+ 
    xlab("") + 
    ylab(paste0("PCA1 (", round(summary(pca)$importance[2, 1]*100, digits = 2), "%)"))+
    # vertical 0 line
    geom_segment(x=0.235, xend=max(res.ind$var.pos), y=0, yend=0, 
                 linetype = "dashed") +
    # x bottom axis
    geom_segment(x=0, xend=0, y=min(vec.pca.axis), yend=max(vec.pca.axis)) + 
    # x top axis
    geom_segment(x=max(res.ind$var.pos), xend=max(res.ind$var.pos), 
                 y=min(res.ind$pca1), yend=max(res.ind$pca1)) + 
    # ticks of the x top axis
    geom_segment(xend=max(res.ind$var.pos)+0.95, data = res.ind, 
                 aes(x = var.pos, y = pca1, yend = pca1_seq), 
                 size = 0.2) + 
    # ticks of the x top axis bis
    geom_segment(xend=max(res.ind$var.pos)-0.3, data = res.ind, 
                 aes(x = var.pos, y = pca1, yend = pca1), 
                 size = 0.2) + 
    # text of the top axis
    annotate(geom = "text", fontface = 'italic', x = res.ind$var.pos+1,
             y = res.ind$pca1_seq, label = res.ind$species,
             angle = 90, hjust = 0, size = 3, color = "black") + 
    coord_flip() + 
    theme(panel.background = element_rect(fill = "white", color = "white"), 
          panel.grid = element_blank(), 
          axis.ticks.y = element_blank(), 
          axis.text.y = element_text(size = 10), 
          legend.position = "none") 
  
  # - Save the plot
  ggsave(file.in, plot.out, width = 14, height = 6, units = "cm", dpi = 600, bg = "white")
  
  # return the name of all the plots made
  return(file.in)
  
}







#' function to map the spatial distribution of storm and fire disturbed pop
#' @param FUNDIV_climate_species FUNDIV data
#' @param climate.in list of climate list objects, where names are disturbances
#' @param file.in Name of the file to save including path
map_climates = function(FUNDIV_climate_species, climate.in, file.in){
  
  # Create directory of file.in if it doesn't exist
  create_dir_if_needed(file.in)
  
  # -- 
  # Start by making the climate pca
  # -- 
  
  # - Make PCA 
  pca <- prcomp(FUNDIV_climate_species[, c("sgdd", "wai")], 
                center = TRUE, scale = TRUE)
  
  # Number of categories for plotting
  n.cat = 15
  
  
  # - Extract the coordinates of the ndividuals on pca axis
  res.ind <- data.frame(plot = FUNDIV_climate_species$plotcode, 
                        pca1 = FUNDIV_climate_species$pca1)  %>%
    mutate(pca1_cut = cut(pca1, breaks = seq(min(pca1), max(pca1), length.out = n.cat)), 
           pca1_min = as.numeric(gsub("\\(", "", gsub("\\,.+", "", pca1_cut))), 
           pca1_max = as.numeric(gsub(".+\\,", "", gsub("\\]", "", pca1_cut))), 
           pca1_median = (pca1_min + pca1_max)/2) %>%
    filter(!is.na(pca1_cut))
  
  # - Extract the coordinates of the variables on pca axis and classify by category
  res.var <- data.frame(var = rownames(get_pca_var(pca)[[1]]), 
                        # Negative because inverse of pca in original data
                        pca1 = get_pca_var(pca)[[1]][, 1]) %>%
    mutate(var.pos = c(1:dim(.)[1]))
  
  # Plot the arrows of the first PCA axis
  plot.var = res.var %>%
    ggplot(aes(x = var.pos, xend = var.pos, y = 0, yend = pca1)) + 
    geom_segment(arrow = arrow(length = unit(0.1, "cm"))) + 
    scale_x_continuous(breaks = res.var$var.pos, 
                       labels = res.var$var, 
                       limits = c(0.5, 2.5)) + 
    xlab("") + geom_hline(yintercept = 0, linetype = "dashed") +
    ylab(paste0("PCA1 (", round(summary(pca)$importance[2, 1]*100, digits = 2), "%)")) + 
    coord_flip() + 
    theme(panel.background = element_rect(fill = "white", color = "black"), 
          panel.grid = element_blank()) 
  
  # Plot the distribution of plots along the pca axis
  plot.ind = res.ind %>%
    group_by(pca1_median) %>%
    summarize(n = n()) %>%
    ggplot(aes(x = pca1_median, y = n, fill = pca1_median)) +
    geom_bar(color = "black", stat = "identity") +
    scale_fill_gradientn(colors = colorRampPalette(c("orange", "blue"))(n.cat)) +
    ylab("Number of\nNFI plots") + 
    theme(panel.background = element_rect(color = "black", fill = "white"), 
          panel.grid = element_blank(), 
          axis.text.x = element_blank(), 
          axis.title.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          legend.position = "none") 
  
  # Plot the pca climate
  plot.pca = plot_grid(plot.ind, plot.var, ncol = 1, align = "v", 
                       rel_heights = c(1, 0.5))
  
  # Memorize the color vector 
  color.data = res.ind %>%
    dplyr::select(pca1_min, pca1_max, pca1_median) %>%
    arrange(pca1_min) %>%
    distinct() %>%
    mutate(color = colorRampPalette(c("orange", "blue"))(dim(.)[1]))
  
  
  
  # -- 
  # Make map
  # -- 
  
  
  
  # Initialize the data 
  data = res.ind %>%
    left_join((FUNDIV_climate_species %>% 
                 dplyr::select(plot = plotcode, longitude, latitude)), 
              by = "plot") %>%
    mutate(climate = 0) %>%
    dplyr::select(plot, longitude, latitude, pca1, pca1_min, pca1_max, 
                  pca1_median, climate) %>%
    st_as_sf(coords = c("longitude", "latitude"), crs = 4326, agr = "constant")
  
  # Loop on all climates to specify which plot is in which climate
  for(j in 1:length(names(climate.in))){
    
    # Identify the range of pca1 associated with climate i
    range.j = as.numeric(quantile(data$pca1, probs = climate.in[[j]]))
    
    # Modify data.i for plots that are in the range
    data = data %>%
      mutate(climate = ifelse(pca1 > range.j[1] & pca1 < range.j[2], j, climate))
  }
  
  # restrict color dataset to the climatic range studied
  color.data.map = color.data %>%
    filter(pca1_median %in% subset(data, climate != 0)$pca1_median)
  
  # Convert climate in a factor
  datatest = data %>%
    mutate(climate = ifelse(climate == 0, "beyond\nselected\nrange", paste0("clim", climate))) %>%
    mutate(climate = factor(climate, levels = c(
      "beyond\nselected\nrange", paste0("clim", c(0:max(data$climate)))
    )))
  
  # Histogram adapted to the climatic gradient
  hist = res.ind %>%
    group_by(pca1_median) %>%
    summarize(n = n()) %>%
    ggplot(aes(x = pca1_median, y = n, fill = pca1_median)) +
    geom_bar(color = NA, stat = "identity") +
    scale_fill_gradientn(colors = (
      color.data %>% mutate(
        color = ifelse(pca1_median %in% color.data.map$pca1_median, color, "gray"))
    )$color) +
    ylab("Number of NFI plots") + 
    theme(panel.background = element_rect(color = "black", fill = "white"), 
          panel.grid = element_blank(), 
          axis.text = element_blank(), 
          axis.title = element_blank(), 
          axis.ticks = element_blank(), 
          legend.position = "none") 
  
  # Map for disturbance i
  plot.map <- ne_countries(scale = "medium", returnclass = "sf") %>%
    ggplot(aes(geometry = geometry)) +
    geom_sf(fill = "#343A40", color = "gray", show.legend = F, size = 0.2) + 
    geom_sf(data = datatest, aes(color = climate), size = 0.0001, shape = 20) +
    scale_color_manual(
      values = c("gray", colorRampPalette(
        c(color.data.map$color[1], color.data.map$color[dim(color.data.map)[1]]))(
          length(names(climate.in))))) +
    coord_sf(xlim = c(-10, 32), ylim = c(36, 71)) +
    theme(panel.background = element_rect(color = 'black', fill = 'white'), 
          panel.grid = element_blank(), 
          legend.title = element_blank(), 
          legend.key = element_blank()) + 
    guides(color = guide_legend(override.aes = list(size=5, alpha = 0.85))) +
    annotation_custom(ggplotGrob(hist), 
                      xmin = -11, xmax = 5, ymin = 61, ymax = 72)
  
  
  
  
  # -- 
  # Make histogram with diversity
  # -- 
  
  
  # Vector of all species
  species_vec = colnames(FUNDIV_climate_species)[grep("_", colnames(FUNDIV_climate_species))]
  # Remove species with bad estimation or unstable in simulations
  species_vec = species_vec[!(species_vec %in% c("Quercus_ilex", "Carpinus_betulus"))]
  # Vector of all species for which we have disturbance parameters
  data("disturb_coef")
  species_vec_dist = (disturb_coef %>%
                        filter(disturbance %in% "storm"))$species
  # Restrict the species vector to these species
  species_vec = species_vec[which(species_vec %in% species_vec_dist)]
  
  
  
  
  # species combinations for this climate based on data (frequency method)
  # -- Paste all species presence absence to create binary code
  eval(parse(text = paste0(
    "data.sp <- FUNDIV_climate_species %>% mutate(combi = paste0(", 
    paste(species_vec, collapse = ", "), "))")))
  # -- Calculate the number of species per species combination
  eval(parse(text = paste0(
    "data.sp <- data.sp %>% mutate(n.sp = ", 
    paste(species_vec, collapse = " + "), ")")))
  # -- Restrict to columns of interest
  data.sp = data.sp %>%
    dplyr::select(plot = plotcode, combi, n.sp)
  # -- count combinations per richness and per climate
  data.plot.sp = datatest %>%
    # Add combinations and richness
    left_join(data.sp, by = "plot") %>%
    filter(climate != "beyond\nselected\nrange" & n.sp > 0) %>%
    st_drop_geometry() %>%
    # Count combinations per richness and climate
    dplyr::select(climate, combi, n.sp) %>%
    distinct() %>%
    group_by(climate, n.sp) %>%
    summarise(n.combi = n()) %>%
    # Divide into combinations selected and not selected
    mutate(included = ifelse(n.combi <= 10, n.combi, 10), 
           not_included = ifelse(n.combi >= 10, n.combi - 10, 0)) %>%
    dplyr::select(-n.combi) %>%
    gather(key = "inclusion", value = "n.combi", "included", "not_included")
  
  # Plot the combinations per richness and per climate
  plot.sp = data.plot.sp %>%
    mutate(inclusion = factor(inclusion, levels = c("not_included", "included"))) %>%
    ggplot(aes(x = n.sp, y = n.combi, fill = climate, color = inclusion)) + 
    geom_bar(stat = "identity", aes(alpha = inclusion)) + 
    facet_wrap(~ climate, nrow = 2) + 
    scale_alpha_manual(values = c(0.8, 1)) +
    scale_color_manual(values = c("gray", "black")) +
    scale_fill_manual(values = colorRampPalette(
      c(color.data.map$color[1], color.data.map$color[dim(color.data.map)[1]]))(
        length(names(climate.in)))) +
    geom_hline(yintercept = 10, linetype = "dashed", color = "red") +
    xlab("Species richness") + ylab("Number of observed\ncombinations") +
    theme(panel.background = element_rect(fill = "white", color = "black"), 
          panel.grid = element_blank(), 
          strip.background = element_blank(), 
          legend.position = "none")
  
  
  # Final plot
  plot.out = plot_grid(plot_grid(plot.pca, plot.map, nrow = 1, rel_widths = c(0.55, 1), 
                                 labels = c("(a)", "(b)"), scale = c(0.85, 0.85)), 
                       plot.sp, nrow = 2, labels = c("", "(c)"), 
                       scale = c(1, 0.85))
  
  
  # Save plot i
  ggsave(file.in, plot.out, width = 17, height = 17, units = "cm", 
         dpi = 1400, bg = "white")
  
  # Return name of the file
  return(file.in)
  
}



#' Function to plot the three metric of forest response
#' @param sim_disturbance_file Name of rds file of disturbance simulation
#' @param file.in Name of the file to save, including path
plot_metrics = function(sim_disturbance_file, file.in){
  
  # Create directory if needed
  create_dir_if_needed(file.in)
  
  # Read simulation i
  sim.i = readRDS(sim_disturbance_file)
  
  # Format the output
  data.i <- sim.i %>%
    # Calculate total basal area per year
    filter(var == "BAsp") %>%
    filter(!equil) %>%
    group_by(time) %>%
    summarize(BA = sum(value)) %>%
    mutate(time = round(time/10)*10) %>%
    ungroup() %>% group_by(time) %>%
    summarize(BA = mean(BA)) %>%
    filter(time <= 3500)
  
  # Time of disturbance
  Tdist = 500
  # Time used to calculate recovery
  recov.time = 50
  # Identify Basal area at equilibrium and after disturbance
  BAeq.i = mean((data.i %>% filter(time < Tdist))$BA)
  BAdist.i = min(data.i$BA)
  # Identify slope and intercept after disturbance
  coef.i = abs(diff((data.i %>% filter(time %in% c(Tdist, Tdist+recov.time)))$BA))/recov.time
  intercept.i = BAdist.i - coef.i*Tdist
  # Span in Basal area (for plot scaling)
  BA.span = diff(range(data.i$BA))
  
  # Colors for the different response metrics
  resistance.color = "#9B2226"
  recovery.color = "#386641"
  resilience.color = "#33658A"
  
  # Make the final plot
  plot.out = data.i %>%
    mutate(lwr = ifelse(BA <= BAeq.i, BA, BAeq.i), 
           upr = ifelse(BA <= BAeq.i, BAeq.i, BA)) %>%
    ggplot(aes(x = time, y = BA, group = 1)) + 
    geom_line() + 
    # RESISTANCE
    annotate('text', x = 100, y = BAeq.i + 0.25*BA.span, parse = TRUE, size=5, hjust = 0,
             label = "Resistance==frac(BA[dist],BA[eq])", color = resistance.color) + 
    annotate('text', x = 10, y = BAeq.i+0.08*BA.span, size=5, hjust = 0, color = resistance.color,
             label = "BA[eq]", parse = TRUE) +
    annotate('text', x = 10, y = BAdist.i+0.08*BA.span, size=5, hjust = 0, color = resistance.color,
             label = "BA[dist]", parse = TRUE) +
    geom_segment(x=0, xend=Tdist, y=BAdist.i, yend=BAdist.i, 
                 color = resistance.color, linetype = "dashed") +
    geom_segment(x=0, xend=3000, y=BAeq.i, yend=BAeq.i, 
                 color = resistance.color, linetype = "dashed") +
    # RESILIENCE
    annotate('text', x = 1200, y = BAeq.i - 0.2*BA.span, parse = TRUE, size=5, 
             hjust = 0, color = resilience.color,
             label = "Resilience==frac(1, integral(sqrt((BA(t) - BA[eq])^2)*dt, t[dist], t[dist+3000]))") +
    geom_ribbon(aes(ymin = lwr, ymax = upr), fill = resilience.color, alpha = 0.4, color = NA) +
    # RECOVERY
    geom_segment(x=Tdist, xend=(Tdist+75), y=BAdist.i, yend=intercept.i+coef.i*(Tdist+75), 
                 color = recovery.color, linetype = "solid", 
                 arrow = arrow(length = unit(0.2, "cm")), size = 1) +
    annotate('text', x = (Tdist + 100), y = (BAdist.i + 0.05*BA.span), parse = TRUE, size=5, 
             color = recovery.color, hjust = 0,
             label = "Recovery==frac((BA(t[dist]+20) - BA[dist]), 20)") +
    # Finish plotting
    geom_segment(x=Tdist, xend=Tdist, y=min(data.i$BA) - 0.3*BA.span, yend=BAdist.i, 
                 color = "black", linetype = "dashed") +
    annotate('text', x = Tdist+25, y = min(data.i$BA) - 0.25*BA.span, size=5, 
             hjust = 0, label = "t[dist]", parse = TRUE) +
    ylim(min(data.i$BA) - 0.3*BA.span, 
         max(data.i$BA) + 0.3*BA.span) + 
    xlab("Time (years)") + ylab("Stand basal area (m2)") + 
    theme(panel.grid = element_blank(), 
          panel.background = element_rect(fill = "white", color = "black"), 
          axis.text = element_text(size = 12), 
          axis.title = element_text(size = 15))
  
  # - Save the plot
  ggsave(file.in, plot.out, width = 17, height = 12, units = "cm", dpi = 600, bg = "white")
  
  # return the name of all the plots made
  return(file.in)
}



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# -- Plots for analyses -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






#' Analyse the data with a structural equation model approach
#' @param data_model formatted model output
#' @param FD_metric Functional diversity metric to choose ("FDis", "FRic or "FD")
#' @param recovery_metric Recovery metric to choose ("recovery", "thalf)
#' @param R_metric Richness metric to choose ("nsp", "H" or "D")
#' @param dir.in name of the directory where to save outputs
plot_sem = function(data_model, FD_metric = "FD", R_metric = "nsp", 
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
           CWM_scaled = as.numeric(scale(CWM, center = TRUE, scale = TRUE)), 
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
    lm(CWM_scaled ~ climate_scaled, data = data.in), 
    lm(resistance.log_scaled ~ FD_scaled + CWM_scaled + H_scaled + climate_scaled, 
       data = data.in), 
    lm(recovery.log_scaled ~ FD_scaled + CWM_scaled + climate_scaled + H_scaled, 
       data = data.in), 
    lm(resilience.log_scaled ~ resistance.log_scaled + recovery.log_scaled + 
         FD_scaled + CWM_scaled + climate_scaled + H_scaled, data = data.in)
  )
  
  
  
  
  
  # Plot the results
  
  # -- Prepare some parameters about the box to drax
  box.height = 2
  box.width = 4
  box.height.spacing = 2
  
  # -- Make a dataset to plot the box with the text
  data.plot.box = data.frame(
    text = c("climate", "FD", "CWM", R_metric, "recovery", "resistance", "resilience"), 
    center.x = c(0, -7, 4.5, -13, 2, -7, -1), 
    height.level = c(4, 3, 3, 3, 2, 2, 1)) %>%
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
        (var.exp %in% c(R_metric, "CWM") & var.resp == "resilience") ~ center.y.resp, 
        (var.exp == "H" & var.resp == "FD") ~ center.y.resp, 
        TRUE ~ (center.y.exp - 0.5*box.height)), 
      arrow.beg.x = ifelse((var.exp == "H" & var.resp == "FD"),
                           (center.x.exp + 0.5*box.width), center.x.exp), 
      arrow.end.y = case_when(
        (var.exp %in% c(R_metric, "CWM") & var.resp == "resilience") ~ center.y.resp, 
        (var.exp == "H" & var.resp == "FD") ~ center.y.resp,
        TRUE ~ (center.y.resp + 0.5*box.height)), 
      arrow.end.x = case_when(
        (var.exp == R_metric & var.resp == "resilience") ~ (center.x.resp - 0.5*box.width), 
        (var.exp == "CWM" & var.resp == "resilience") ~ (center.x.resp + 0.5*box.width),
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
                 mutate(text_stat = paste0(text, "\n", "R2 = ", R2))), 
              by = "text") %>%
    mutate(text_stat = ifelse(text == "climate", "Hot-dry to\ncold-wet", text_stat))
  
  # Add rectangles
  data.box.cat = data.plot.box %>%
    mutate(cat = case_when(
      text == "climate" ~ "Climate", 
      text %in% c(R_metric, "FD", "CWM") ~ "Species\ncomposition", 
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
    geom_segment(data = subset(data.plot.arrow, var.exp %in% c(R_metric, "CWM") & 
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
  
  
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # Add additional plot with boxplots
  
  fig.file.bxplt.in = paste0(dir.in, "/fig_sem_full.jpg")
  #' Sub-function to plot changes in a given variable with climate
  #' @param data_model data formatted 
  #' @param var.in variable to plot (resistance, resilience or recovery)
  #' @param lab.in name of the label on the plot
  plot.bxplt = function(data_model, var.in, lab.in){
    
    data_model %>%
      rename("recov" = recovery_metric) %>%
      mutate(climate = paste0("clim", ID.climate)) %>%
      mutate(climate = factor(
        climate, levels = paste0("clim", c(1:length(unique(.$ID.climate)))))) %>%
      dplyr::select(climate, pca1, resistance, recovery = recov, resilience) %>%
      gather(key = "variable", value = "value", "resistance", "recovery", "resilience") %>%
      mutate(value = ifelse(variable == "resistance", value, log(value))) %>%
      filter(variable == var.in) %>%
      ggplot(aes(x = climate, y = value, fill = climate)) + 
      geom_boxplot(alpha = 0.7, outlier.shape = 21) + 
      scale_fill_manual(values = colorRampPalette(c("orange", "blue"))(10)) +
      ylab(lab.in) + 
      theme(legend.position = "none", 
            panel.grid = element_blank(), 
            panel.background = element_rect(fill = "white", color = "black"), 
            strip.background = element_blank(), 
            axis.text.x = element_text(angle = 60, hjust = 1), 
            axis.title.x = element_blank()) 
  }
  
  # Build the boxplots
  plot.resistance = plot.bxplt(data_model, "resistance", "logit(Resistance)")
  plot.recovery = plot.bxplt(data_model, "recovery", "log(Recovery)")
  plot.resilience = plot.bxplt(data_model, "resilience", "log(Resilience)")
  plot.bx = plot_grid(plot.resistance, plot.recovery, plot.resilience, 
                      align = "h", scale = 0.95, labels = c("(b)", "(c)", "(d)"), 
                      nrow = 1)
  plot.out.full = plot_grid(plot.out, plot.bx, ncol = 1, labels = c("(a)", ""), 
                            rel_heights = c(1, 0.4), scale = 0.9)
  ggsave(fig.file.bxplt.in, plot.out.full, width = 20, height = 18, units = "cm", 
         dpi = 600, bg = "white")
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  
  # - Save the plot
  ggsave(fig.file.in, plot.out, width = 20, height = 12, units = "cm", dpi = 600, bg = "white")
  
  # return the name of all the plots made
  return(c(fig.file.in))
  
}



#' Function to estimate of FD metrics on resilience
#' @param data_model df formatted to fit model
#' @param dir.in Name of the directory where to save files
plot_FD_effect_resilience = function(data_model, dir.in){
  # Name of the files
  fig.file.in = paste0(dir.in, "/fig_FD_effect_resilience_storm.jpg")
  fig.file.predictions = paste0(dir.in, "/observed_vs_predicted_H1.jpg")
  fig.file.residuals = paste0(dir.in, "/residuals.jpg")
  fig.file.observations = paste0(dir.in, "/observations.jpg")
  table.file.stats = paste0(dir.in, "/stats_H1.tex")
  
  # create output directory if it doesn't exist
  create_dir_if_needed(fig.file.in)
  
  # Vector of response variables for which to run models
  response.vec = c("resilience", "resistance", "recovery")
  
  # Remove from data model data with na
  data_model = data_model %>% filter(!is.na(CWM))
  
  # Data to fit the models
  data.in = cbind(data_model[, response.vec], 
                  scale((data_model %>% dplyr::select(
                    "H.scaled" = "H", "FD.scaled" = "FD", "CWM.scaled" = "CWM")), 
                    center = TRUE, scale = TRUE)) %>%
    mutate(resilience = log(resilience), 
           recovery = log(recovery))
  
  # Initialize the table with statistics
  table.stats = data.frame(col1 = "", col2 = "Est", col3 = "se", col4 = "F", col5 = "p", col6 = "VIF")
  
  # Loop on all response variables
  for(j in 1:length(response.vec)){
    
    # Fit model
    eval(parse(text = paste0(
      "model.j = lm(", response.vec[j], " ~ H.scaled + FD.scaled + CWM.scaled, data = data.in)"
    )))
    
    # Data for predictions vs observations
    data.predict.j = data.in %>%
      mutate(predicted = predict(model.j, newdata = .)) %>%
      rename("observed" = response.vec[j]) %>%
      mutate(response.var = response.vec[j], 
             residuals = model.j$residuals, 
             H = data_model$H, FD = data_model$FD, CWM = data_model$CWM) %>%
      dplyr::select(response.var, H, H.scaled, FD, FD.scaled, CWM, CWM.scaled, 
                    predicted, observed, residuals)
    
    # Models to unscale explanatory variables
    unscale.H = lm(H ~ H.scaled, data.predict.j)
    unscale.FD = lm(FD ~ FD.scaled, data.predict.j)
    unscale.CWM = lm(CWM ~ CWM.scaled, data.predict.j)
    
    # Data for fit
    for(i in 1:3){
      var.exp.i = c("H", "FD", "CWM")[i]
      data.fit.i = data.frame(
        H.scaled = seq(from = min(data.in$H.scaled), to = max(data.in$H.scaled), length.out = 100), 
        FD.scaled = seq(from = min(data.in$FD.scaled), to = max(data.in$FD.scaled), length.out = 100), 
        CWM.scaled = seq(from = min(data.in$CWM.scaled), to = max(data.in$CWM.scaled), length.out = 100), 
        exp.var = var.exp.i, 
        response.var = response.vec[j]) %>%
        # Set to the mean for variables that are not var i
        mutate(H.scaled = ifelse(exp.var == "H", H.scaled, mean(data.predict.j$H.scaled)), 
               FD.scaled = ifelse(exp.var == "FD", FD.scaled, mean(data.predict.j$FD.scaled)), 
               CWM.scaled = ifelse(exp.var == "CWM", CWM.scaled, mean(data.predict.j$CWM.scaled))) %>%
        # Add unscaled values
        mutate(H = predict(unscale.H, newdata = .), 
               FD = predict(unscale.FD, newdata = .), 
               CWM = predict(unscale.CWM, newdata = .)) %>%
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
      var.exp = c("H", "FD", "CWM"), 
      var.pos = c(1:3),
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
    
    # Add to the final output dataset
    if(j == 1){
      data.out = data.out.j
      data.predict = data.predict.j
    } else {
      data.out = rbind(data.out, data.out.j)
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
    dplyr::select(response.var, residuals, H, FD, CWM) %>%
    gather(key = "exp.var", value = "exp.var.value", "H", "FD", "CWM") %>%
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
    dplyr::select(response.var, observed, predicted, H, FD, CWM) %>%
    gather(key = "exp.var", value = "exp.var.value", "H", "FD", "CWM") %>%
    mutate(observed = ifelse(response.var == "resistance", 
                             plogis(observed), exp(observed))) %>%
    ggplot(aes(x = exp.var.value, y = observed, fill = exp.var, 
               color = exp.var, group = interaction(response.var, exp.var))) + 
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
  plot.out = data.out %>%
    mutate(significance = ifelse(est.low > 0 | est.high < 0, "yes", "no")) %>%
    mutate(var.resp = factor(var.resp, levels = c("resistance", "recovery", "resilience"))) %>%
    ggplot(aes(x = var.exp, y = est, color = significance, fill = var.resp)) + 
    geom_errorbar(aes(ymin = est.low, ymax = est.high),
                  width = 0) + 
    geom_point(shape = 21, size = 3) +
    xlab("") + ylab("Effect on response to\ndisturbance metric") +
    scale_color_manual(values = c(`no` = "gray", `yes` = "black")) +
    scale_fill_manual(values = c("#9B2226", "#386641", "#33658A")) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    facet_grid(. ~ var.resp, scales = "free") +
    theme(panel.background = element_rect(color = "black", fill = "white"), 
          panel.grid = element_blank(), 
          strip.background = element_blank(), 
          strip.text = element_text(face = "bold"), 
          legend.position = "none") + 
    coord_flip()
  
  
  
  # Save plot 
  ggsave(fig.file.in, plot.out, width = 13, height = 5, units = "cm", 
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
plot_FD_effect_vs_climate_quadra = function(
  data_model, R_metric = "nsp", mod_selection = "pvalue", dir.in){
  
  # create output directory if it doesn't exist
  create_dir_if_needed(paste0(dir.in, "/test"))
  
  # Vector of response variables for which to run models
  response.vec = c("resilience", "resistance", "recovery")
  
  # Vector of explanatory variables that are not climate
  vec.exp = c("H", "FD", "CWM")
  
  # Initialize data with the effect pvalue
  data.text = expand.grid(var.resp = response.vec, 
                          var.exp = vec.exp, 
                          text = NA)
  
  # Data to fit the models for disturbance i
  data.in = cbind(data_model[, response.vec], 
                  scale((data_model %>% 
                           dplyr::select("H" = R_metric, "FD", "CWM", "Clim" = "pca1")), 
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
    mutate(var.resp = factor(var.resp, levels = c("resistance", "recovery", "resilience")))
  
  
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
    mutate(effect = factor(effect, levels = c("H effect", "FD effect", "CWM effect"))) %>%
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
  
  # Plot the estimates
  plot.out = data.out %>%
    mutate(effect = factor(effect, levels = c("H effect", "FD effect", "CWM effect"))) %>%
    ggplot(aes(x = pca1, y = mean, group = clim, fill = clim)) + 
    geom_ribbon(aes(ymin = lwr, ymax = upr), color = NA, alpha = 0.5) +
    geom_line(color = "#001524") + 
    xlab("Coordinate on the sgdd-wai PCA\n(Hot-dry to cold-wet climatic gradient)") + 
    ylab("Effect on response to\ndisturbance metric") +
    facet_grid(effect ~ var.resp, scales = "free") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_text(aes(label = text), data = data.text, inherit.aes = TRUE, 
              hjust = "inward", size = 2.5, alpha = 0.8) +
    scale_fill_manual(values = colorRampPalette(c("orange", "blue"))(10)) +
    ylim(c(min(data.out$lwr), max(data.out$upr) + 0.2*diff(range(data.out$mean)))) +
    theme(panel.background = element_rect(color = "black", fill = "white"), 
          panel.grid = element_blank(), 
          strip.background = element_blank(), 
          strip.text = element_text(face = "bold"), 
          legend.position = "none") + 
    # Add point for the effect without interactions
    geom_point(data = data.simple, aes(color = var.resp), 
               inherit.aes = TRUE, shape = 16) +
    geom_errorbar(data = data.simple, inherit.aes = TRUE, width = 0, 
                  aes(color = var.resp, ymin = lwr, ymax = upr)) +
    scale_color_manual(values = c("#9B2226", "#386641", "#33658A")) + 
    geom_vline(xintercept = max(data.out$pca1), linetype = "solid", color = "grey")
  
  # Name of the plot to save
  file.plot = paste0(dir.in, "/fd_effect_climate.jpg")
  file.plot.predictions = paste0(dir.in, "/observed_vs_predicted.jpg")
  file.plot.fitclim = paste0(dir.in, "/data_and_fit_clim.jpg")
  
  # Save plots
  ggsave(file.plot, plot.out, width = 16, height = 12 , units = "cm", 
         dpi = 600, bg = "white")
  ggsave(file.plot.predictions, plot.predictions, width = 14, height = 5, 
         units = "cm", dpi = 600, bg = "white")
  ggsave(file.plot.fitclim, plot.fit.clim, width = 14, height = 5, 
         units = "cm", dpi = 600, bg = "white")
  
  # Return name of the file saved
  return(c(file.plot, file.plot.predictions, file.plot.fitclim, files.stat.table))
}

  
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# -- Plots for supplementary -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



#' Function to show co-variations between resilience metrics
#' @param data_model Data from simulations formatted
#' @param var.in vector of variables to include in the pairwise analysis
#' @param file.in Name of the file to save, including path
plot_covariation = function(data_model, var.in, file.in){
  
  # Create directory if needed
  create_dir_if_needed(file.in)
  
  # Save the plot
  png(file = file.in, width = 20, height = 18,
       unit = "cm", res = 600)
  
  # Make the plot
  pairs((data_model %>%
           mutate(recovery = log(recovery), 
                  resilience = log(resilience))%>% 
           dplyr::select(var.in) ), pch = 19, lower.panel = NULL)
  
  
  # Return name of the file generated 
  return(file.in)
  
}


#' Plot the relation between tree diversity and tree density at equilibrium
#' @param data_model dataset formatted for the analyses
#' @param file.in Name of the file to save, including path
plot_tree_packing = function(data_model, file.in){
  
  # Create directory if needed
  create_dir_if_needed(file.in)
  
  # Plot the relation between diversity and density
  plot.out = data_model %>%
    mutate(climate = paste0("clim", ID.climate)) %>%
    mutate(climate = factor(
      climate, levels = paste0("clim", c(1:length(unique(.$ID.climate)))))) %>%
    ggplot(aes(x = H, y = Nha, fill = climate)) + 
    geom_point(shape = 21, color = "black", alpha = 0.5) + 
    scale_fill_manual(values = colorRampPalette(c("orange", "blue"))(10)) +
    facet_wrap(~ climate, nrow = 2) + 
    geom_smooth(method = "lm", color = "red", fill = NA) + 
    theme(legend.position = "none", 
          panel.grid = element_blank(), 
          panel.background = element_rect(fill = "white", color = "black"), 
          strip.background = element_blank())
  
  # Save the plot
  ggsave(file.in, plot.out, width = 17, height = 8, 
         units = "cm", dpi = 600, bg = "white")
  
  # Return the name of the file
  return(file.in)
}


#' Plot the relation between climate, diversity and structure
#' @param data_model dataset formatted for the analyses
#' @param file.in Name of the file to save, including path
plot_climate_vs_diversity_and_structure = function(data_model, file.in){
  
  # Create directory if needed
  create_dir_if_needed(file.in)
  
  # Make the plot
  plot.out = data_model %>%
    mutate(climate = paste0("clim", ID.climate)) %>%
    mutate(climate = factor(
      climate, levels = paste0("clim", c(1:length(unique(.$ID.climate)))))) %>%
    mutate(dbh_mean_diff = dbh_mean_postdist - dbh_mean, 
           dbh_q10_diff = dbh_q10_postdist - dbh_q10, 
           dbh_q90_diff = dbh_q90_postdist - dbh_q90) %>%
    dplyr::select(climate, H, FD, CWM, BA_diff, BA_eq, dbh_mean, Nha,
                  dbh_mean_diff, dbh_q10_diff, dbh_q90_diff) %>%
    drop_na() %>%
    gather(key = "variable", value = "value", "H", "FD", "CWM", "BA_eq", "Nha", 
           "dbh_mean", "dbh_mean_diff", "dbh_q10_diff", "dbh_q90_diff") %>%
    mutate(variable = factor(variable, levels = c(
      "H", "FD", "CWM", "BA_eq", "dbh_mean", "Nha", "dbh_mean_diff", 
      "dbh_q10_diff", "dbh_q90_diff"))) %>%
    ggplot(aes(x = climate, y = value, fill = climate)) + 
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
  ggsave(file.in, plot.out, width = 17, height = 18, 
         units = "cm", dpi = 600, bg = "white")
  
  # Return the name of the file
  return(file.in)
}


#' Analyse the data with a structural equation model approach
#' @param data_model formatted model output
#' @param FD_metric Functional diversity metric to choose ("FDis", "FRic or "FD")
#' @param recovery_metric Recovery metric to choose ("recovery", "thalf)
#' @param R_metric Richness metric to choose ("nsp", "H" or "D")
#' @param dir.in name of the directory where to save outputs
plot_sem_supp = function(data_model, FD_metric = "FD", R_metric = "H", 
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
           CWM_scaled = as.numeric(scale(CWM, center = TRUE, scale = TRUE)), 
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
    lm(CWM_scaled ~ climate_scaled + FD_scaled + H_scaled, data = data.in), 
    lm(resistance.log_scaled ~ FD_scaled + CWM_scaled + H_scaled + climate_scaled, 
       data = data.in), 
    lm(recovery.log_scaled ~ FD_scaled + CWM_scaled + climate_scaled + H_scaled + resistance.log_scaled, 
       data = data.in), 
    lm(resilience.log_scaled ~ resistance.log_scaled + recovery.log_scaled + 
         FD_scaled + CWM_scaled + climate_scaled + H_scaled, data = data.in)
  )
  
  
  # Plot the results
  
  # -- Prepare some parameters about the box to drax
  box.height = 2
  box.width = 4
  box.height.spacing = 2
  
  # -- Make a dataset to plot the box with the text
  data.plot.box = data.frame(
    text = c("climate", "FD", "CWM", R_metric, "recovery", "resistance", "resilience"), 
    center.x = c(0, -7, 4.5, -13, 2, -7, -1), 
    height.level = c(4, 3, 3, 3, 2, 2, 1)) %>%
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
        (var.exp %in% c(R_metric, "CWM") & var.resp == "resilience") ~ center.y.resp, 
        (var.exp == "H" & var.resp == "FD") ~ center.y.resp, 
        TRUE ~ (center.y.exp - 0.5*box.height)), 
      arrow.beg.x = ifelse((var.exp == "H" & var.resp == "FD"),
                           (center.x.exp + 0.5*box.width), center.x.exp), 
      arrow.end.y = case_when(
        (var.exp %in% c(R_metric, "CWM") & var.resp == "resilience") ~ center.y.resp, 
        (var.exp == "H" & var.resp == "FD") ~ center.y.resp,
        TRUE ~ (center.y.resp + 0.5*box.height)), 
      arrow.end.x = case_when(
        (var.exp == R_metric & var.resp == "resilience") ~ (center.x.resp - 0.5*box.width), 
        (var.exp == "CWM" & var.resp == "resilience") ~ (center.x.resp + 0.5*box.width),
        (var.exp == "H" & var.resp == "FD") ~ (center.x.resp - 0.5*box.width),
        TRUE ~ center.x.resp)) %>%
    # Report significance, sign and absolute value of estimate
    mutate(signif = ifelse(p <= 0.05, "yes", "no"), 
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
                 mutate(text_stat = paste0(text, "\n", "R2 = ", R2))), 
              by = "text") %>%
    mutate(text_stat = ifelse(text == "climate", "Hot-dry to\ncold-wet", text_stat))
  
  # Add rectangles
  data.box.cat = data.plot.box %>%
    mutate(cat = case_when(
      text == "climate" ~ "Climate", 
      text %in% c(R_metric, "FD", "CWM") ~ "Species\ncomposition", 
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
    filter(!(var.resp == "CWM" & var.exp == "H")) %>%
    filter(!(var.resp == "CWM" & var.exp == "FD")) %>%
    filter(!(var.resp == "recovery" & var.exp == "resistance")) 
  
  # -- Plot the boxes
  plot.box = data.plot.box %>%
    ggplot(aes(x = center.x, y = center.y))  + 
    geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = cat),
              data = data.box.cat,  color = "black", inherit.aes = TRUE, 
              show.legend = FALSE, alpha = 0.5)+ 
    geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), 
              color = "black", fill = "#415A77", alpha = 0.7)+ 
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
    geom_segment(data = subset(data.plot.arrow, var.exp %in% c(R_metric, "CWM") & 
                                 var.resp == "resilience"), 
                 aes(x = center.x.exp, xend = center.x.exp, 
                     y = center.y.exp - 0.5*box.height, yend = center.y.resp, 
                     size = magnitude, linetype = signif, color = sign)) + 
    # Arrows
    geom_segment(data = data.plot.arrow, 
                 aes(x = arrow.beg.x, xend = arrow.end.x, 
                     y = arrow.beg.y, yend = arrow.end.y, 
                     size = magnitude, linetype = signif, color = sign), 
                 arrow = arrow(length = unit(0.15, "cm")), 
                 type = "closed") + 
    # -- scale linetype, color and size
    scale_color_manual(values = c(`negative` = "#0081A7", `positive` = "#F07167")) + 
    scale_linetype_manual(values = c(`yes` = "solid", `no` = "dashed")) + 
    scale_size_manual(values = c(`low` = 0.3, `mid` = 0.6, `high` = 1.1))
  
  
  # - Save the plot
  ggsave(fig.file.in, plot.out, width = 20, height = 12, units = "cm", dpi = 600, bg = "white")
  
  # return the name of all the plots made
  return(c(fig.file.in))
}




#' Plot a list of simulations
#' @param sim_disturbance.in list of simulations
#' @param file.in Name including path of the file to save
plot_sim_dist <- function(sim_disturbance.in, file.in){
  
  # Create directory if needed
  create_dir_if_needed(file.in)
  
  # Loop on all species combinations
  for(i in 1:length(sim_disturbance.in)){
    
    # Printer
    print(paste0(i, "/", length(sim_disturbance.in)))
    
    # Read simulation
    sim.i = readRDS(sim_disturbance.in[i])
    
    # Vector of species for simulation i
    sp.vec.i <- unique(sim.i$species)
    
    # Format data for plotting simulation i
    data.i <- sim.i %>%
      # Remove unused lines and columns
      filter(var == "BAsp") %>%
      dplyr::select(species, time, value) %>%
      distinct() %>%
      # Calculate the sum of basal area
      spread(key = "species", value = "value") %>%
      mutate(all = rowSums(.[, c(2:dim(.)[2])])) %>%
      gather(key = "species", value = "value", c(sp.vec.i, "all")) %>%
      # Add simulation ID
      mutate(ID = i)
    
    # Add to final dataset
    if(i == 1) data = data.i
    else data <- rbind(data, data.i)
  }
  
  # Make the plot 
  plot.out <- data %>%
    mutate(species = factor(
      species, levels = c("all", unique((filter(data, species != "all"))$species)))) %>%
    ggplot(aes(x = time, y = value, group = species, color = species)) + 
    geom_line() + 
    facet_wrap(~ ID) + 
    theme(panel.background = element_rect(fill = "white", color = "black"), 
          panel.grid = element_blank(), 
          strip.background = element_blank(), 
          strip.text = element_blank(),
          legend.title = element_blank(), 
          legend.key = element_blank()) + 
    xlab("Time (years)") + ylab("Basal area")
  
  # - Save the plot
  ggsave(file.in, plot.out, width = 21, height = 14, units = "cm", dpi = 600, bg = "white")
  
  # return the name of all the plots made
  return(file.in)
}


