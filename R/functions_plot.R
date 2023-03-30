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
#' @param file.in name including path of the file to save
plot_traits_pca <- function(traits, file.in){
  
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
  
  # Make the plot
  plot.out <- res.var %>%
    ggplot(aes(x = var.pos, xend = var.pos, y = 0, yend = pca1)) + 
    geom_segment(arrow = arrow(length = unit(0.1, "cm"))) + 
    scale_x_continuous(breaks = res.var$var.pos, 
                       labels = res.var$var, 
                       limits = c(0.5, max(res.ind$var.pos)+3))+ 
    scale_y_continuous(breaks = vec.pca.axis)+ 
    xlab("") + 
    ylab(paste0("PCA1 (", round(summary(pca)$importance[2, 1]*100, digits = 2), "%)"))+
    # vertical 0 line
    geom_segment(x=0.235, xend=max(res.ind$var.pos), y=0, yend=0, 
                 linetype = "dashed") +
    # x bottom axis
    geom_segment(x=0.22, xend=0.22, y=min(vec.pca.axis), yend=max(vec.pca.axis)) + 
    # x top axis
    geom_segment(x=max(res.ind$var.pos), xend=max(res.ind$var.pos), 
                 y=min(res.ind$pca1), yend=max(res.ind$pca1)) + 
    # ticks of the x top axis
    geom_segment(xend=max(res.ind$var.pos)+0.25, data = res.ind, 
                 aes(x = var.pos, y = pca1, yend = pca1_seq), 
                 size = 0.2) + 
    # ticks of the x top axis bis
    geom_segment(xend=max(res.ind$var.pos)-0.1, data = res.ind, 
                 aes(x = var.pos, y = pca1, yend = pca1), 
                 size = 0.2) + 
    # text of the top axis
    annotate(geom = "text", fontface = 'italic',
             x = res.ind$var.pos+0.30, y = res.ind$pca1_seq, label = res.ind$species,
             angle = 90, hjust = 0, size = 3) + 
    coord_flip() + 
    theme(panel.background = element_rect(fill = "white", color = "white"), 
          panel.grid = element_blank(), 
          axis.ticks.y = element_blank()) 
  
  # - Save the plot
  ggsave(file.in, plot.out, width = 14, height = 8, units = "cm", dpi = 600, bg = "white")
  
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
                        pca1 = -get_pca_var(pca)[[1]][, 1]) %>%
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
    scale_fill_gradientn(colors = colorRampPalette(c("blue", "orange"))(n.cat)) +
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
    mutate(color = colorRampPalette(c("blue", "orange"))(dim(.)[1]))
  
  
  
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
  data = data %>%
    mutate(climate = factor(climate, levels = as.character(c(0:max(.$climate)))))
  
  # Histogram adapted to the climatic gradient
  hist = res.ind %>%
    group_by(pca1_median) %>%
    summarize(n = n()) %>%
    ggplot(aes(x = pca1_median, y = n, fill = pca1_median)) +
    geom_bar(color = "black", stat = "identity") +
    scale_fill_gradientn(colors = (
      color.data %>% mutate(
        color = ifelse(pca1_median %in% color.data.map$pca1_median, color, "gray"))
    )$color) +
    ylab("Number of\nNFI plots") + 
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
    geom_sf(data = data, aes(color = climate), size = 0.01, alpha = 0.7) +
    scale_color_manual(
      values = c("gray", colorRampPalette(
        c(color.data.map$color[1], color.data.map$color[dim(color.data.map)[1]]))(
          length(names(climate.in))+1))) +
    coord_sf(xlim = c(-10, 32), ylim = c(36, 71)) +
    theme(panel.background = element_rect(color = 'black', fill = 'white'), 
          panel.grid = element_blank(), 
          legend.position = "none") +
    annotation_custom(ggplotGrob(hist), 
                      xmin = -10, xmax = 4, ymin = 62, ymax = 72)
  
  
  # Final plot
  plot.out = plot_grid(plot.pca, plot.map, nrow = 1, rel_widths = c(0.7, 1), 
                       labels = c("(a)", "(b)"), scale = c(0.85, 0.85))
  
  
  # Save plot i
  ggsave(file.in, plot.out, width = 18, height = 13, units = "cm", 
         dpi = 600, bg = "white")
  
  # Return name of the file
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
#' @param file.in name of the file to save, including path
plot_sem = function(data_model, FD_metric = "FD", R_metric = "nsp", 
                    recovery_metric = "recovery", file.in){
  
  # Create directory if needed
  create_dir_if_needed(file.in)
  
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
           SR_scaled = as.numeric(scale(R_chosen, center = TRUE, scale = TRUE)),
           resistance.log_scaled = as.numeric(scale(resistance.log, center = TRUE, scale = TRUE)), 
           recovery.log_scaled = as.numeric(scale(recovery.log, center = TRUE, scale = TRUE)), 
           resilience.log_scaled = as.numeric(scale(resilience.log, center = TRUE, scale = TRUE)))
  
  
  # -- Make model
  mod_sem = psem(
    glm(FD_scaled ~ climate_scaled, family = tweedie(var.power = 2), 
        data = data.in), 
    lm(CWM_scaled ~ climate_scaled, data = data.in), 
    lm(resistance.log_scaled ~ FD_scaled + CWM_scaled + SR_scaled, 
       data = data.in), 
    lm(recovery.log_scaled ~ FD_scaled + CWM_scaled + climate_scaled + SR_scaled, 
       data = data.in), 
    lm(resilience.log_scaled ~ resistance.log_scaled + recovery.log_scaled + 
         FD_scaled + CWM_scaled + climate_scaled + SR_scaled, data = data.in)
  )
  
  
  
  
  
  # Plot the results
  
  # -- Prepare some parameters about the box to drax
  box.height = 2
  box.width = 4
  box.height.spacing = 2
  
  # -- Make a dataset to plot the box with the text
  data.plot.box = data.frame(
    text = c("climate", "FD", "CWM", "SR", "recovery", "resistance", "resilience"), 
    center.x = c(0, -4, 4, -10, 1, -7, -2), 
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
      arrow.beg.y = ifelse((var.exp %in% c("SR", "CWM") & var.resp == "resilience"),
                           center.y.resp, (center.y.exp - 0.5*box.height)), 
      arrow.beg.x = center.x.exp, 
      arrow.end.y = ifelse((var.exp %in% c("SR", "CWM") & var.resp == "resilience"),
                           center.y.resp, (center.y.resp + 0.5*box.height)), 
      arrow.end.x = case_when(
        (var.exp == "SR" & var.resp == "resilience") ~ (center.x.resp - 0.5*box.width), 
        (var.exp == "CWM" & var.resp == "resilience") ~ (center.x.resp + 0.5*box.width),
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
  
  # Add stat to the plot box data
  data.plot.box = data.plot.box %>%
    left_join((data.plot.arrow %>% 
                 dplyr::select(text = var.resp, class, R2) %>% 
                 mutate(R2 = round(R2, digits = 3)) %>% 
                 distinct() %>% 
                 mutate(text_stat = paste0(text, "\n", "R2 = ", R2))), 
              by = "text") %>%
    mutate(text_stat = ifelse(is.na(text_stat), text, text_stat))
  
  # Add rectangles
  data.box.cat = data.plot.box %>%
    mutate(cat = case_when(
      text == "climate" ~ "Climate", 
      text %in% c("SR", "FD", "CWM") ~ "Species\ncomposition", 
      text %in% c("resistance", "recovery", "resilience") ~ "Resilience"
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
              color = "black", fill = "#415A77", alpha = 0.7)+ 
    geom_text(aes(label = text_stat), color = "white", size = 4) + 
    geom_text(aes(label = cat), data = data.box.cat,
              hjust = 0, size = 4, inherit.aes = TRUE) +
    scale_fill_manual(values = c("#FF6F59", "#7678ED", "#43AA8B")) +
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
    geom_segment(data = subset(data.plot.arrow, var.exp %in% c("SR", "CWM") & 
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
  ggsave(file.in, plot.out, width = 20, height = 12, units = "cm", dpi = 600, bg = "white")
  
  # return the name of all the plots made
  return(file.in)
  
}



#' Function to estimate of FD metrics on resilience
#' @param data_model df formatted to fit model
#' @param R_metric Richness metric to choose ("nsp", "H" or "D")
#' @param file.in Name of the file to save, inlcuding path
plot_FD_effect_resilience = function(data_model, R_metric = "nsp", file.in){
  
  # create output directory if it doesn't exist
  create_dir_if_needed(file.in)
  
  # Vector of response variables for which to run models
  response.vec = c("resilience", "resistance", "recovery")
  
  
  
  # Data to fit the models
  data.in = cbind(data_model[, response.vec], 
                  scale((data_model %>% dplyr::select("R" = R_metric, "FD", "CWM")), 
                        center = TRUE, scale = TRUE)) %>%
    mutate(resilience = log(resilience), 
           recovery = log(recovery))
  
  # Loop on all response variables
  for(j in 1:length(response.vec)){
    
    # Fit model
    eval(parse(text = paste0(
      "model.j = lm(", response.vec[j], " ~ R + FD + CWM, data = data.in)"
    )))
    
    # Output data set for model i j 
    data.out.j = data.frame(
      var.resp = response.vec[j], 
      var.exp = c("R", "FD", "CWM"), 
      var.pos = c(1:3),
      est = as.numeric(coef(model.j)[-1]), 
      est.low = as.numeric(confint(model.j)[-1, 1]), 
      est.high = as.numeric(confint(model.j)[-1, 2])
    )
    
    # Add to the final output dataset
    if(j == 1) data.out = data.out.j
    else data.out = rbind(data.out, data.out.j)
  }
  
  
  # Plot the estimates
  plot.out = data.out %>%
    mutate(significance = ifelse(est.low > 0 | est.high < 0, "yes", "no")) %>%
    mutate(var.resp = factor(var.resp, levels = c("resistance", "recovery", "resilience"))) %>%
    ggplot(aes(x = var.exp, y = est, color = significance)) + 
    geom_errorbar(aes(ymin = est.low, ymax = est.high),
                  width = 0) + 
    geom_point(shape = 21, fill = "#2A9D8F") +
    xlab("") + ylab("Effect on resilience metric") +
    scale_color_manual(values = c(`no` = "gray", `yes` = "black")) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    facet_grid(. ~ var.resp, scales = "free") +
    theme(panel.background = element_rect(color = "black", fill = "white"), 
          panel.grid = element_blank(), 
          strip.background = element_blank(), 
          strip.text = element_text(face = "bold"), 
          legend.position = "none") + 
    coord_flip()
  
  # Save plot 
  ggsave(file.in, plot.out, width = 13, height = 5, units = "cm", 
         dpi = 600, bg = "white")
  
  # Return name of the file saved
  return(file.in)
  
}


#' Function to estimate of FD metrics on resilience
#' @param data_model df formatted to fit model
#' @param R_metric Richness metric to choose ("nsp", "H" or "D")
#' @param file.in Name of the file to save, inlcuding path
plot_FD_and_climate_effect_resilience = function(data_model, R_metric = "nsp", file.in){
  
  # create output directory if it doesn't exist
  create_dir_if_needed(file.in)
  
  # Vector of response variables for which to run models
  response.vec = c("resilience", "resistance", "recovery")
  
  
  
  # Data to fit the models for disturbance i
  data.in = cbind(data_model[, response.vec], 
                  scale((data_model %>% 
                           dplyr::select("R" = R_metric, "FD", "CWM", "Clim" = "pca1")), 
                        center = TRUE, scale = TRUE)) %>%
    mutate(resilience = log(resilience), 
           recovery = log(recovery))
  
  
  
  # Data to fit the models for disturbance i
  data.in = cbind(data_model[, response.vec], 
                  scale((data_model %>% 
                           dplyr::select(R = nsp, FD, CWM, Clim = pca1)), 
                        center = TRUE, scale = TRUE)) %>%
    mutate(resilience = log(resilience), 
           recovery = log(recovery))
  
  # Loop on all response variables
  for(j in 1:length(response.vec)){
    
    # Fit model
    eval(parse(text = paste0(
      "model.j = lm(", response.vec[j], " ~ R*Clim + FD*Clim + CWM*Clim, data = data.in)"
    )))
    
    # Output data set for model i j 
    data.out.j = data.frame(
      var.resp = response.vec[j], 
      var.exp = names(coef(model.j))[-1], 
      est = as.numeric(coef(model.j)[-1]), 
      est.low = as.numeric(confint(model.j)[-1, 1]), 
      est.high = as.numeric(confint(model.j)[-1, 2])
    )
    
    # Add to the final output dataset
    if(j == 1) data.out = data.out.j
    else data.out = rbind(data.out, data.out.j)
  }
  
  
  
  # Plot the estimates
  plot.out = data.out %>%
    mutate(significance = ifelse(est.low > 0 | est.high < 0, "yes", "no")) %>%
    mutate(var.resp = factor(var.resp, levels = c("resistance", "recovery", "resilience"))) %>%
    ggplot(aes(x = var.exp, y = est, color = significance)) + 
    geom_errorbar(aes(ymin = est.low, ymax = est.high),
                  width = 0) + 
    geom_point(shape = 21, fill = "#2A9D8F") +
    xlab("") + ylab("Effect on resilience metric") +
    scale_color_manual(values = c(`no` = "gray", `yes` = "black")) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    facet_grid(. ~ var.resp, scales = "free") +
    theme(panel.background = element_rect(color = "black", fill = "white"), 
          panel.grid = element_blank(), 
          strip.background = element_blank(), 
          strip.text = element_text(face = "bold"), 
          legend.position = "none") + 
    coord_flip()
  
  # Save plot i
  ggsave(file.in, plot.out, width = 14, height = 7, units = "cm", 
         dpi = 600, bg = "white")
  
  # Return name of the file saved
  return(file.in)
  
}




#' Function to estimate of FD metrics on resilience
#' @param data_model df formatted to fit model
#' @param R_metric Richness metric to choose ("nsp", "H" or "D")
#' @param file.in Name of the file to save, inlcuding path
plot_predictions = function(data_model, R_metric = "nsp", file.in){
  
  # create output directory if it doesn't exist
  create_dir_if_needed(file.in)
  
  # Vector of response variables for which to run models
  response.vec = c("resilience", "resistance", "recovery")
  
  
  
  # Data to fit the models for disturbance i
  data.in = cbind(data_model[, response.vec], 
                  scale((data_model %>% 
                           dplyr::select("R" = R_metric, "FD", "CWM", "Clim" = "pca1")), 
                        center = TRUE, scale = TRUE)) %>%
    # Convert strictly positive variables to log scale
    mutate(resilience = log(resilience), 
           recovery = log(recovery))
  
  # Model to scale climate
  scale.clim = lm(pca1 ~ Clim, 
                  data = data.frame(pca1 = data_model$pca1, 
                                    Clim = data.in$Clim))
  
  # Loop on all response variables
  for(j in 1:length(response.vec)){
    
    # Fit model
    eval(parse(text = paste0(
      "model.j = lm(", response.vec[j], " ~ R*Clim + FD*Clim + CWM*Clim, data = data.in)"
    )))
    
    # Output data set for model i j 
    data0.j = expand.grid(
      Clim = seq(from = min(data.in$Clim), to = max(data.in$Clim), length.out = 100),
      var.resp = response.vec[j], 
      exp.quantile = c(0.2, 0.8))
    
    # Rescale climate
    data0.j$pca1 = predict(scale.clim, newdata = data0.j)
    
    # Vector of explanatory variables that are not climate
    vec.exp = c("R", "FD", "CWM")
    
    # Loop on the three explanatory variables
    for(k in 1:length(vec.exp)){
      
      # Add name of the explanatory variable
      data.jk = data0.j %>% mutate(var.exp = vec.exp[k])
      
      # Write variable k with different value depending on quantile
      eval(parse(text = paste0(
        "data.jk = data.jk %>% mutate(", vec.exp[k] ,"= quantile(data.in$", 
        vec.exp[k], ", probs = exp.quantile, na.rm = TRUE))"
      )))
      
      # Other variables have their mean value
      eval(parse(text = paste0(
        "data.jk = data.jk %>% mutate(", 
        paste(paste0(vec.exp[-k], " = mean(data.in$", vec.exp[-k], ", na.rm = TRUE)"), 
              collapse = ", "), ")"
      )))
      
      # Predict resilience metric and confidence interval
      data.jk = cbind.data.frame(data.jk, predict(
        model.j, newdata = data.jk, type = "response", 
        interval = "confidence", level = 0.95)) %>%
        gather(key = "variable", value = "value", "fit", "lwr", "upr") %>%
        mutate(value = ifelse(var.resp == "resistance", value, exp(value))) %>%
        spread(key = "variable", value = "value")
      
      # Add to the final output dataset
      if(j == 1 & k == 1) data.out = data.jk
      else data.out = rbind(data.out, data.jk)
      
    }
    
  }
  
  
  # Plot the estimates
  plot.out = data.out %>%
    mutate(FD_metric = paste0(100*exp.quantile, "% quantile")) %>%
    mutate(var.resp = factor(var.resp, levels = c("resistance", "recovery", "resilience"))) %>%
    ggplot(aes(x = pca1, y = fit, group = FD_metric, 
               alpha = FD_metric, linetype = FD_metric)) + 
    geom_line(color = "#001524") + 
    geom_ribbon(aes(ymin = lwr, ymax = upr), color = NA, fill = "#2A9D8F") +
    xlab("Coordinate on the sgdd-wai PCA") + ylab("Resilience metric") +
    scale_alpha_manual(values = c(0.3, 0.6)) +
    scale_linetype_manual(values = c("dashed", "solid")) +
    facet_grid(var.resp ~ var.exp, scales = "free") +
    theme(panel.background = element_rect(color = "black", fill = "white"), 
          panel.grid = element_blank(), 
          strip.background = element_blank(), 
          strip.text = element_text(face = "bold"), 
          legend.key = element_blank())
  
  # Save plot i
  ggsave(file.in, plot.out, width = 18, height = 9 , units = "cm", 
         dpi = 600, bg = "white")
  
  # Return name of the file saved
  return(file.in)
  
}





#' Function to estimate of FD metrics on resilience
#' @param data_model df formatted to fit model
#' @param R_metric Richness metric to choose ("nsp", "H" or "D")
#' @param file.in Name of the file to save, inlcuding path
plot_FD_effect_vs_climate = function(data_model, R_metric = "nsp", file.in){
  
  # create output directory if it doesn't exist
  create_dir_if_needed(file.in)
  
  # Vector of response variables for which to run models
  response.vec = c("resilience", "resistance", "recovery")
  
  # Vector of explanatory variables that are not climate
  vec.exp = c("R", "FD", "CWM")
  
  # Initialize data with the effect pvalue
  data.text = expand.grid(var.resp = response.vec, 
                          var.exp = vec.exp, 
                          text = NA)
  
  # Data to fit the models for disturbance i
  data.in = cbind(data_model[, response.vec], 
                  scale((data_model %>% 
                           dplyr::select("R" = R_metric, "FD", "CWM", "Clim" = "pca1")), 
                        center = TRUE, scale = TRUE)) %>%
    # Convert strictly positive variables to log scale
    mutate(resilience = log(resilience), 
           recovery = log(recovery))
  
  # Model to scale climate
  scale.clim = lm(pca1 ~ Clim, 
                  data = data.frame(pca1 = data_model$pca1, 
                                    Clim = data.in$Clim))
  
  # Initialize data with climatic gradient
  data.clim = data.frame(
    Clim = seq(from = min(data.in$Clim), to = max(data.in$Clim), length.out = 100)) %>%
    mutate(pca1 = predict(scale.clim, newdata = .)) %>%
    rename(pca1.scaled = Clim)
  
  
  # Loop on all response variables
  for(j in 1:length(response.vec)){
    
    # Fit model
    eval(parse(text = paste0(
      "model.j = lm(", response.vec[j], " ~ R*Clim + FD*Clim + CWM*Clim, data = data.in)"
    )))
    
    # Initialize table with coefficients value
    data.coef = as.data.frame(matrix(
      NA, nrow = 10000, ncol = length(model.j$coefficients), dimnames = list(
        c(), gsub("\\:", "\\.", names(model.j$coefficients)))))
    
    # Loop on all coefficients 
    for(i in 2:dim(data.coef)[2]){
      # Estimate
      i.mean = summary(model.j)$coefficients[i, 1]
      # Standard deviation
      i.sd = summary(model.j)$coefficients[i, 2]
      # Generate distribution of coefficient (normal distribution)
      data.coef[, i] = rnorm(n = dim(data.coef)[1], mean = i.mean, sd = i.sd)
    }
    # Remove intercept column
    data.coef = data.coef[, -1]
    
    # Merge with climate data
    data.j = expand_grid(data.clim, data.coef) %>%
      mutate(var.resp = response.vec[j]) %>%
      mutate(R_effect = R + R.Clim*pca1.scaled, 
             FD_effect = FD + Clim.FD*pca1.scaled, 
             CWM_effect = CWM + Clim.CWM*pca1.scaled) %>%
      dplyr::select(pca1, R_effect, FD_effect, CWM_effect) %>%
      gather(key = "effect", value = "value", "R_effect", "FD_effect", "CWM_effect") %>%
      group_by(pca1, effect) %>%
      summarise(mean = mean(value, na.rm = TRUE), 
                lwr = quantile(value, 0.025, na.rm = TRUE), 
                upr = quantile(value, 0.975, na.rm = TRUE)) %>%
      mutate(var.resp = response.vec[j])
    
    # Fill the text dataset
    for(k in vec.exp){
      id_text_jk = which(data.text$var.exp == k & data.text$var.resp == response.vec[j])
      data.text[id_text_jk, "text"] = paste0(
        k, ": ", pvalue(summary(model.j)$coefficients[k, 4], add_p = TRUE, accuracy = 0.01), 
        "\n", k, "*Clim: ", pvalue(summary(model.j)$coefficients[
          intersect(grep(k, names(model.j$coefficients)), 
                    grep("Clim", names(model.j$coefficients))), 4], add_p = TRUE, accuracy = 0.01))
    }
    
    # Add to the final output dataset
    if(j == 1) data.out = data.j
    else data.out = rbind(data.out, data.j)
    
  }
  
  # Add position on x and y axis for data.text
  data.text = data.text %>%
    mutate(pca1 = min(data.out$pca1) + 0.05*diff(range(data.out$pca1)), 
           mean = max(data.out$mean) - 0.03*diff(range(data.out$mean))) %>%
    mutate(effect = paste0(var.exp, " effect"))
  
  # Plot the estimates
  plot.out = data.out %>%
    mutate(effect = gsub("\\_", " ", effect)) %>%
    mutate(var.resp = factor(var.resp, levels = c("resistance", "recovery", "resilience"))) %>%
    ggplot(aes(x = pca1, y = mean, group = 1)) + 
    geom_ribbon(aes(ymin = lwr, ymax = upr), color = NA, fill = "#2A9D8F", alpha = 0.5) +
    geom_line(color = "#001524") + 
    xlab("Coordinate on the sgdd-wai PCA") + ylab("Effect on resilience metric") +
    facet_grid(effect ~ var.resp) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_text(aes(label = text), data = data.text, inherit.aes = TRUE, 
              hjust = "inward", size = 2.5, alpha = 0.8) +
    theme(panel.background = element_rect(color = "black", fill = "white"), 
          panel.grid = element_blank(), 
          strip.background = element_blank(), 
          strip.text = element_text(face = "bold"), 
          legend.key = element_blank(), 
          legend.title = element_blank())
  
  # Save plot i
  ggsave(file.in, plot.out, width = 13, height = 12 , units = "cm", 
         dpi = 600, bg = "white")
  
  # Return name of the file saved
  return(file.in)
  
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# -- Exploratory plots -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Function to plot changes in CWM and FD along time till equilibrium
#' @param sim_equilibrium vector of the filenames of simulations saved as rds
#' @param forest_list dataset with information on all forests simulated
#' @param pc1_per_species df with coordinates of the climate pca per species
#' @param file.in Name of the file to save, including path
plot_cwm_fd_overtime = function(sim_equilibrium, forest_list, pc1_per_species, 
                                file.in){
  # Create the directories of file.in if needed
  create_dir_if_needed(file.in)
  
  # Initialize a counter of successful simulations
  k = 0
  
  # Loop on all simulations
  for(i in 1:length(sim_equilibrium)){
    
    # Printer
    print(paste0("Getting data for forest ", i, "/", length(sim_equilibrium)))
    
    # Extract the data for simulation i
    data.i = readRDS(sim_equilibrium[i]) %>%
      filter(equil == FALSE, var == "BAsp") %>%
      left_join(pc1_per_species, by = "species")
    
    # Check that the simulations reached equilibrium
    if(!is.na(sum(data.i$value))){
      
      # If equilibrium is reached, calculate cwm and fd for each time step
      data.i = data.i  %>%
        group_by(time) %>%
        summarize(CWM = weighted.mean(pca1, w = value, na.rm = TRUE), 
                  FD = weighted.var(pca1, w = value, na.rm = TRUE)) %>%
        mutate(sim.number = i, 
               ID.climate = forest_list$ID.climate[i], 
               sp.richness = length(unlist(strsplit(forest_list$combination[i], "\\."))))
      
      # Add to the final dataset
      if(k == 0) data = data.i
      else data = rbind(data, data.i)
      
      # The simulation was successful: increase the counter
      k = k+1
    } 
    
  }
  
  # Plot the data for CWM
  plot.cwm = data %>%
    mutate(climate = paste0("climate.", ID.climate)) %>%
    mutate(climate = factor(
      climate, levels = paste0("climate.", c(1:length(unique(data$ID.climate))))
    )) %>%
    ggplot(aes(x = time, y = CWM, group = sim.number, color = sp.richness)) + 
    geom_line() + 
    facet_wrap(~ climate, nrow = 2) + 
    scale_color_gradient(low = "blue", high = "red") + 
    theme(panel.background = element_rect(fill = "white", color = "black"), 
          panel.grid = element_blank(), 
          strip.background = element_blank())
  
  # Plot the data for functional diversity
  plot.fd = data %>%
    mutate(climate = paste0("climate.", ID.climate)) %>%
    mutate(climate = factor(
      climate, levels = paste0("climate.", c(1:length(unique(data$ID.climate))))
    )) %>%
    ggplot(aes(x = time, y = FD, group = sim.number, color = sp.richness)) + 
    geom_line() + 
    facet_wrap(~ climate, nrow = 2) + 
    scale_color_gradient(low = "blue", high = "red") + 
    theme(panel.background = element_rect(fill = "white", color = "black"), 
          panel.grid = element_blank(), 
          strip.background = element_blank())
  
  # Final plot
  plot.out = plot_grid(plot.cwm, plot.fd, nrow = 2, labels = c("(a)", "(b)"))
  
  # Save the plot
  ggsave(file.in, plot.out, width = 30, height = 20, units = "cm", 
         dpi = 600, bg = "white")
  
  # Return the name of the file saved
  return(file.in)
  
}




#' Function to show the distribution of mean and var pca values for
#' combination of species selected randomly or based on occurence
#' @param climate list of climate objects
#' @param pc1_per_species df containing climate pca coordinates per species
#' @param file.in name of the file to save, including path
plot_pca1_selection_vs_random = function(climate, pc1_per_species, file.in){
  
  # Create directory if needed
  create_dir_if_needed(file.in)
  
  # Loop on all climate
  for(i in 1:length(names(climate))){
    
    # Vector of available species for climate i
    sp.vec.i = climate[[i]]$species
    
    # df with the data for selected combinations
    data.selected.i = data.frame(combination = climate[[i]]$combinations) %>%
      mutate(source = "selected") %>%
      mutate(sp.richness = unlist(lapply(
        .[, "combination"], function(x) length(unlist(strsplit(x, "\\."))))))
    
    
    # Loop to identify all possible combinations for each level of species richness
    for(j in 1:max(data.selected.i$sp.richness)){
      # combinations for richness j
      data.random.ij = data.frame(
        combination = apply(as.data.frame(t(combn(sp.vec.i, j))), 1, paste, collapse = "." )
      ) %>%
        mutate(source = "random", sp.richness = j)
      
      # Add to the final df with random data
      if(j == 1) data.random.i = data.random.ij
      else data.random.i = rbind(data.random.i, data.random.ij)
    }
    
    # Final data for climate i: bind the two df
    data.i = rbind(data.selected.i, data.random.i) %>%
      # Create one line per species for each forest
      mutate(ID.climate = i, ID.forest = c(1:dim(.)[1])) %>%
      cbind(as.data.frame(matrix(0, nrow = dim(.)[1], ncol = length(sp.vec.i), 
                                 dimnames = list(NULL, sp.vec.i)))) %>%
      gather(key = "species", value = "present", sp.vec.i)
    # Only keep the species present in the forest 
    for(k in 1:dim(data.i)[1]) data.i$present[k] = ifelse(
      grepl(data.i$species[k], data.i$combination[k]), 1, 0)
    data.i = data.i %>%
      filter(present == 1) %>%
      dplyr::select(ID.forest, ID.climate, source, combination, sp.richness, species) %>%
      arrange(ID.forest)
    
    # Add to the final output data
    if(i == 1) data = data.i
    else data = rbind(data, data.i)
    
  }
  
  # Final formatting
  data.out = data %>%
    left_join(pc1_per_species, by = "species") %>%
    group_by(ID.forest, ID.climate, combination, sp.richness, source) %>%
    summarize(mean = mean(pca1), 
              var = var(pca1))  %>%
    mutate(climate = paste0("climate.", ID.climate)) %>%
    mutate(climate = factor(
      climate, levels = paste0("climate.", c(1:length(unique(data$ID.climate))))
    ))
  
  # Plot for mean pca value
  plot.mean = data.out %>%
    ggplot(aes(x = mean, fill = source)) + 
    geom_density(color = "black", aes(y = stat(density)), alpha = 0.5) + 
    facet_wrap(~ climate, nrow = 2) + 
    xlab("Mean PCA value") + 
    scale_fill_manual(values = c("black", "red")) +
    theme(panel.background = element_rect(fill = "white", color = "black"), 
          panel.grid = element_blank(), 
          strip.background = element_blank(), 
          legend.position = "none")
  
  # Plot for mean pca value
  plot.var = data.out %>%
    ggplot(aes(x = var, fill = source)) + 
    geom_density(color = "black", aes(y = stat(density)), alpha = 0.5) + 
    facet_wrap(~ climate, nrow = 2) + 
    xlab("Mean PCA variance") + 
    scale_fill_manual(values = c("black", "red")) +
    theme(panel.background = element_rect(fill = "white", color = "black"), 
          panel.grid = element_blank(), 
          strip.background = element_blank(),
          legend.title = element_blank(),
          legend.position = "none")
  
  # Final plot
  plot.out = plot_grid(
    plot_grid(plot.mean, plot.var, nrow = 2, labels = c("(a)", "(b)"), scale = 0.95), 
    get_legend(plot.var + theme(legend.position = "right")), 
    nrow = 1, rel_widths = c(1, 0.2)
  )
  
  # Save the plot
  ggsave(file.in, plot.out, width = 18, height = 12, units = "cm", 
         dpi = 600, bg = "white")
  
  # Return the name of the file saved
  return(file.in)
  
}


#' Function to plot the proportion of species and of trees for which we have 
#' disturbance data per climate.
#' @param FUNDIV_climate_species df containing species presence and climate per FUNDIV plot 
#' @param disturbance.in "storm", "fire" or "biotic"
#' @param file.in Name of the file to save, including path
plot_prop.species_per_climate = function(
  FUNDIV_climate_species, disturbance.in = "storm", file.in){
  
  
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
  
  
  # - Extract the coordinates of the individuals on pca axis
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
                        pca1 = -get_pca_var(pca)[[1]][, 1]) %>%
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
    scale_fill_gradientn(colors = colorRampPalette(c("blue", "orange"))(n.cat)) +
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
    mutate(color = colorRampPalette(c("blue", "orange"))(dim(.)[1]))
  # Vector of all species
  species_vec = colnames(FUNDIV_climate_species)[grep("_", colnames(FUNDIV_climate_species))]
  
  
  
  
  
  
  # -- 
  # Plots of the proportion of trees per climate
  # -- 
  
  
  # Load vector with disturbance parameters per species and disturbance
  data("disturb_coef")
  
  # Create a climate list with 15 quantile
  climate_list.in = create_climate_list(15, quantile.range = c(0, 1))
  
  # Threshold of tree percentage above which we exclude disturbance
  tree.percent.max = 85
  
  # Vector of all species for which we have disturbance parameters
  species_vec_dist = (disturb_coef %>%
                        filter(disturbance == disturbance.in))$species
  
  # Loop on all climate
  for(i in 1:length(names(climate_list.in))){
    
    # Format data for climate i
    data.ij = FUNDIV_climate_species %>%
      filter(pca1 > quantile(FUNDIV_climate_species$pca1, climate_list.in[[i]][1])) %>%
      filter(pca1 < quantile(FUNDIV_climate_species$pca1, climate_list.in[[i]][2])) %>%
      gather(key = "species", value = "present", species_vec) %>%
      group_by(species) %>%
      summarize(n = sum(present)) %>%
      mutate(present.in.data = ifelse(n == 0, 0, 1),
             present.in.dist = ifelse(species %in% species_vec_dist, 1, 0),
             n.present.in.data = present.in.data*n,
             n.present.in.dist = n*present.in.dist,
             pca1.mean = quantile(FUNDIV_climate_species$pca1, sum(climate_list.in[[i]])/2)) %>%
      ungroup() %>% group_by(pca1.mean) %>%
      summarise(percent.trees.present = sum(n.present.in.dist)/sum(n.present.in.data)*100) %>%
      mutate(disturbance = disturbance.in)
    
    # Add color
    data.ij = data.ij %>%
      mutate(color = ifelse(percent.trees.present <= tree.percent.max, "gray", 
                            (color.data %>% 
                               filter(pca1_min < data.ij$pca1.mean) %>%
                               filter(pca1_max > data.ij$pca1.mean))$color))
    
    # Add to the final dataset
    if(i == 1) data.j = data.ij
    else data.j = rbind(data.j, data.ij)
  }
  
  # Make the plot for disturbance j
  plot.prop.trees = data.j %>%
    mutate(pca1.factor = as.character(pca1.mean)) %>% 
    ggplot(aes(x = pca1.mean, y = percent.trees.present, fill = pca1.factor)) + 
    geom_point(shape = 21, color = "black") + 
    scale_fill_manual(values = setNames(data.j$color, as.character(data.j$pca1.mean))) + 
    ylim(0, 100) + 
    ylab(paste0("Percentage of trees with\nan estimation of ", 
                disturbance.in, " sensitivity")) + 
    xlab(paste0("PCA1 (", round(summary(pca)$importance[2, 1]*100, digits = 2), "%)")) + 
    geom_hline(yintercept = tree.percent.max, linetype = "dashed") + 
    theme(panel.background = element_rect(fill = "white", color = "black"), 
          panel.grid = element_blank(), 
          legend.position = "none")
  
  
  
  
  
  # Final plot 
  plot.out = plot_grid(plot.pca, plot.prop.trees, nrow = 1, align = "h", 
                       scale = c(1, 0.9))
  
  # - Save the plot
  ggsave(file.in, plot.out, width = 14, height = 9, units = "cm", dpi = 600, bg = "white")
  
  # return the name of all the plots made
  return(file.in)
  
}


#' Function to show co-variations between resilience metrics
#' @param data_model Data from simulations formatted
#' @param var.in vector of variables to include in the pairwise analysis
#' @param file.in Name of the file to save, including path
plot_covariation_FD = function(
  data_model, var.in = c("FD", "FRic", "FDis", "CWM", "nsp", "H", "D"), file.in){
  
  # Create directory if needed
  create_dir_if_needed(file.in)
  
  # Make the plot
  plot.out =  data_model %>%
    dplyr::select(var.in) %>%
    ggpairs(aes(alpha = 0.4)) + 
    theme(panel.background = element_rect(color = "black", fill = "white"), 
          panel.grid = element_blank(), 
          legend.key = element_blank(), 
          strip.background = element_blank(), 
          strip.text = element_text(face = "bold"))
  
  # Save the plot
  ggsave(file.in, plot.out, width = 20, height = 18, 
         units = "cm", dpi = 600, bg = "white")
  
  # Return name of the file generated 
  return(file.in)
  
}