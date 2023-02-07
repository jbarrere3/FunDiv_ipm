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


#' Function to plot the sgdd wai pca and the climate defined
#' @param FUNDIV_climate_species data with species presence and climate per FUNDIV plot
#' @param climate.list list of climates to show in the pca
#' @param file.in name of the fileto save, including path
plot_pca_climate = function(FUNDIV_climate_species, climate.list, file.in){
  
  # Create directory if needed
  create_dir_if_needed(file.in)
  
  # Remake pca to plot the arrows
  pca = prcomp(FUNDIV_climate_species[, c("sgdd", "wai")], 
               center = TRUE, scale = TRUE)
  
  # Extract results for individuals
  data_pca_ind = data.frame(plotcode = FUNDIV_climate_species$plotcode) %>%
    cbind(as.data.frame(get_pca_ind(pca)[[1]])) %>%
    rename("pca1" = "Dim.1", "pca2" = "Dim.2") %>%
    mutate(climate = "none")
  
  # Extract results for variables
  data_var_pca = get_pca_var(pca)[[1]] %>%
    as.data.frame() %>%
    mutate(var = rownames(.)) %>%
    rename("pca1" = "Dim.1", "pca2" = "Dim.2") %>%
    mutate(pca1 = pca1*2, pca2 = pca2*2)
  
  # Loop on all climates to specify which plot is in which climate
  for(i in 1:length(names(climate.list))){
    
    # Identify the range of pca1 associated with climate i
    range.i = as.numeric(quantile(data_pca_ind$pca1, probs = climate.list[[i]]))
    
    # Modify data.in for plots that are in the range
    data_pca_ind = data_pca_ind %>%
      mutate(climate = ifelse(pca1 > range.i[1] & pca1 < range.i[2], 
                              names(climate.list)[i], climate))
  }
  
  
  # Plot the pca
  plot.out = data_pca_ind %>%
    filter(climate != "none") %>%
    mutate(climate = factor(climate, levels = names(climate.list))) %>%
    ggplot(aes(x = pca1, y = pca2, color = climate)) + 
    geom_point(alpha = 0.5) + 
    geom_segment(aes(x = 0, y = 0, xend = pca1, yend = pca2), 
                 data = data_var_pca, type = "closed", color = "black",
                 arrow = arrow(length = unit(0.1, "cm"))) + 
    geom_text(data = data_var_pca, aes(label = var), size = 5, color = "black", 
              nudge_y = 0.1) +
    scale_color_manual(values = colorRampPalette(c("blue", "red"))(length(names(climate.list)))) +
    geom_hline(size = 0.2, yintercept = 0, color = "#6C757D", linetype = "dashed") + 
    geom_vline(size = 0.2, xintercept = 0, color = "#6C757D", linetype = "dashed") + 
    xlab(paste0("PCA1 (", round(summary(pca)$importance[2, 1]*100, digits = 2), "%)")) +
    ylab(paste0("PCA2 (", round(summary(pca)$importance[2, 2]*100, digits = 2), "%)")) +
    theme(panel.background = element_rect(color = "black", fill = "white"), 
          panel.grid = element_blank(), 
          axis.title = element_text(size = 15), 
          legend.key = element_blank())
  
  # - Save the plot
  ggsave(file.in, plot.out, width = 16, height = 13, 
         units = "cm", dpi = 600, bg = "white")
  
  # return the name of all the plots made
  return(file.in)
}






#' Function to make a map of Fundiv plots with different climates
#' @param FUNDIV_climate_species Data with species presence and climate per plot
#' @param climate.list list of range (0 to 1) along the pca axis hot/dry to cold/wet
#' @param file.in name including path of the file to save
map_climates = function(FUNDIV_climate_species, climate.list, file.in){
  
  # Create directory if needed 
  create_dir_if_needed(file.in)
  
  # Initialize the data for plotting
  data.in = FUNDIV_climate_species %>%
    mutate(climate = "none") %>%
    dplyr::select(plotcode, longitude, latitude, pca1, climate) %>%
    st_as_sf(coords = c("longitude", "latitude"), crs = 4326, agr = "constant")
  
  # Loop on all climates to specify which plot is in which climate
  for(i in 1:length(names(climate.list))){
    
    # Identify the range of pca1 associated with climate i
    range.i = as.numeric(quantile(data.in$pca1, probs = climate.list[[i]]))
    
    # Modify data.in for plots that are in the range
    data.in = data.in %>%
      mutate(climate = ifelse(pca1 > range.i[1] & pca1 < range.i[2], 
                              names(climate.list)[i], climate))
  }
  
  # Final formatting
  data.in = data.in  %>% 
    filter(climate != "none") %>%
    # Factorize climate list for plotting
    mutate(climate = factor(climate, levels = names(climate.list)))
  
  ## - Final plot
  plot.out <- ne_countries(scale = "medium", returnclass = "sf") %>%
    ggplot(aes(geometry = geometry)) +
    geom_sf(fill = "#343A40", color = "gray", show.legend = F, size = 0.2) + 
    geom_sf(data = (data.in %>% filter(climate != "none")), 
            aes(color = climate), size = 0.01, alpha = 0.7) +
    scale_color_manual(
      values = colorRampPalette(c("blue", "red"))(length(names(climate.list)))) +
    coord_sf(xlim = c(-10, 32), ylim = c(36, 71)) +
    theme(panel.background = element_rect(color = 'black', fill = 'white'), 
          panel.grid = element_blank(),
          legend.title = element_blank(),
          legend.key = element_blank())
  
  # Save the plot
  ggsave(file.in, plot.out, width = 14, height = 14, 
         units = "cm", dpi = 600, bg = "white")
  
  # Return the name of the file
  return(file.in)
  
}




#' Plot the effect of functional strategy on resilience metrics (H1)
#' @param data_models data with functional diversity, resilience metrics
#' @param file.in name including path of the file to save
plot_resilience_vs_CMW_and_FD = function(data_models, file.in){
  
  # Create directory if needed
  create_dir_if_needed(file.in)
  
  # Format data
  data.in = data_models %>%
    filter(resistance > 0) %>%
    mutate(CWM.scaled = scale(.$CWM), 
           FD.scaled = scale(.$FD)) %>%
    gather(key = "variable", value = "value",  
           "resistance", "resilience", "recovery") %>%
    filter(is.finite(value))
  
  # model to unscale CWM
  mod.unscale.CWM = lm(CWM ~ CWM.scaled, data = subset(data.in, var = "resistance"))
  mod.unscale.FD = lm(FD ~ FD.scaled, data = subset(data.in, var = "resistance"))
  
  # Initialize list of plots and of models
  list_plots = list()
  list_plots_effect = list()
  
  # Loop on all variables tested
  for(i in 1:length(unique(data.in$variable))){
    
    # Variable i
    var.i = unique(data.in$variable)[i]
    
    # Make model
    mod.i = lm(log(value) ~ CWM.scaled*FD.scaled, 
               data = subset(data.in, variable == var.i))
    
    # Data fit
    data_fit.i = expand.grid(CWM.scaled = seq(from = min(data.in$CWM.scaled), 
                                              to = max(data.in$CWM.scaled), 
                                              length.out = 1000), 
                             FD.scaled = quantile(
                               subset(data.in, variable == var.i)$FD.scaled, 
                               probs = c(0.1, 0.5, 0.9)
                             )) %>%
      mutate(CWM = predict(mod.unscale.CWM, newdata = .), 
             FD = predict(mod.unscale.FD, newdata = .)) %>%
      cbind(exp(predict(mod.i, newdata = ., interval = "confidence", 
                        level = 0.95))) %>%
      rename(value = fit)
    
    # Maximum value that can be plotted (so that variation are visible)
    max.i = 2*max(c(max(data_fit.i$value, na.rm = TRUE), 
                    max(subset(data.in, variable == var.i)$value, na.rm = TRUE)))
    data_fit.i = data_fit.i %>%
      mutate(upr = ifelse(upr > max.i, max.i, upr))
    
    # Plot values
    plot.i = data_fit.i %>%
      ggplot(aes(x = CWM, y = value, group = FD, color = FD, fill = FD)) + 
      geom_line(show.legend = FALSE) + 
      geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.3, color = NA) + 
      geom_point(data = (subset(data.in, variable == var.i) %>%
                           filter(value <= max(data_fit.i))), shape = 21, 
                 color = "black", inherit.aes = TRUE, size = 2, alpha = 0.7) + 
      scale_fill_gradientn(colors = colorRampPalette(c("#DCE1DE", "#3E8914"))(5), 
                           values = quantile(data.in$FD, c(0, 0.2, 0.5, 0.8, 1))) +
      scale_color_gradientn(colors = colorRampPalette(c("#DCE1DE", "#3E8914"))(5), 
                            values = quantile(data.in$FD, c(0, 0.2, 0.5, 0.8, 1))) +
      ylab(toupper(var.i)) + 
      theme(panel.background = element_rect(fill = "white", color = "black"), 
            panel.grid = element_blank(), 
            legend.key = element_blank(), 
            legend.position = "none", 
            legend.text = element_text(size = 14), 
            legend.title = element_text(size = 16))
    
    # Plot confidence interval around estimates
    plot.effect.i = data.frame(
      var = c("Int", "CWM", "FD", "CWM:FD"), 
      est = coef(summary(mod.i))[, 1], 
      low = confint(mod.i)[, 1], 
      high = confint(mod.i)[, 2]
    ) %>%
      filter(var != "Int") %>%
      mutate(signif = ifelse(low > 0 | high < 0, "yes", "no")) %>%
      mutate(signif = factor(signif, levels = c("no", "yes"))) %>%
      ggplot(aes(x = var, y = est, color = signif)) + 
      geom_point() + 
      geom_errorbar(aes(ymin = low, ymax = high), width = 0) + 
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray") + 
      scale_color_manual(values = c(`no` = "gray", `yes` = "black")) + 
      theme(panel.background = element_rect(fill = "white", color = "black"), 
            panel.grid = element_blank(), 
            legend.position = "none") + 
      xlab("") + ylab("") +
      coord_flip()
    
    # Add to the output lists
    eval(parse(text = paste0("list_plots$", var.i, " = plot.i")))
    eval(parse(text = paste0("list_plots_effect$", var.i, " = plot.effect.i")))
  }
  
  # Final plot
  plot.out <- plot_grid(
    plot_grid(
      plot_grid(plotlist = list_plots_effect, align = "h", nrow = 1, scale = 0.9, 
                labels = paste0("(", letters[c(1:length(unique(data.in$variable)))], ")")), 
      plot_grid(plotlist = list_plots, align = "h", nrow = 1, scale = 0.9), 
      nrow = 2, rel_heights = c(0.5, 1), align = "v"
    ), 
    get_legend(plot.i + theme(legend.position = "left")), 
    nrow = 1, rel_widths = c(1, 0.1)
  )
  
  # Save the plot
  ggsave(file.in, plot.out, width = 24, height = 10, units = "cm", 
         dpi = 600, bg = "white")
  
  # return the name of the file
  return(file.in)
  
}


#' Plot the effect of functional strategy on resilience metrics (H1)
#' @param data_models data with functional diversity, resilience metrics
#' @param file.in name including path of the file to save
plot_resilience_vs_climate = function(data_models, file.in){
  
  # Create directory if needed
  create_dir_if_needed(file.in)
  
  # Format data
  data.in = data_models %>%
    filter(resistance > 0) %>%
    gather(key = "variable", value = "value",  
           "resistance", "resilience", "recovery") %>%
    mutate(ID.climate = factor(ID.climate, levels = unique(.$ID.climate))) %>%
    filter(is.finite(value))
  
  
  
  # Initialize list of plots
  list_plots = list()
  
  # Loop on all variables tested
  for(i in 1:length(unique(data.in$variable))){
    
    # Variable i
    var.i = unique(data.in$variable)[i]
    
    # Make model
    mod.i = lm(log(value) ~ pca1, 
               data = subset(data.in, variable == var.i))
    
    # Data fit
    data_fit.i = expand.grid(pca1 = seq(from = min(subset(data.in, variable == var.i)$pca1), 
                                        to = max(subset(data.in, variable == var.i)$pca1), 
                                        length.out = 1000), 
                             climate = NA_character_) %>%
      cbind(exp(predict(mod.i, newdata = ., interval = "confidence", 
                        level = 0.95))) %>%
      rename(value = fit)
    
    
    # Make a color vector
    color.vec = colorRampPalette(c("blue", "red"))(length(unique(data.in$ID.climate)))
    names(color.vec) = paste0("climate.", c(1:length(unique(data.in$ID.climate))))
    
    # Make plot i
    plot.i = subset(data.in, variable == var.i) %>%
      group_by(ID.climate, pca1) %>%
      summarize(lwr = quantile(value, 0.025), 
                upr = quantile(value, 0.975), 
                value = mean(value, na.rm = TRUE)) %>%
      mutate(climate = paste0("climate.", ID.climate)) %>%
      ggplot(aes(x = pca1, y = value, ymin = lwr, ymax = upr, 
                 fill = climate, group = 1)) +
      geom_errorbar(width = 0) +
      geom_point(shape = 21, color = "black") + 
      geom_line(data = data_fit.i, inherit.aes = TRUE) + 
      geom_ribbon(data = data_fit.i, inherit.aes = TRUE, alpha = 0.3, color = NA) + 
      scale_fill_manual(values = color.vec) + 
      xlab("CLIMATE") + ylab(toupper(var.i)) + 
      ggtitle(paste0(
        "F = ", round(anova(mod.i)[1, 4], digits = 1), 
        ", ", scales::pvalue(anova(mod.i)[1, 5], accuracy = 0.01, add_p = TRUE)
      )) +
      theme(panel.background = element_rect(fill = "white", color = "black"), 
            panel.grid = element_blank(), 
            legend.position = "none", 
            legend.title = element_blank(), 
            axis.text.x = element_blank(), 
            axis.ticks.x = element_blank())
    
    # Add to the output lists
    eval(parse(text = paste0("list_plots$", var.i, " = plot.i")))
  }
  
  
  # Final plot
  plot.out = plot_grid(plotlist = list_plots, align = "hv", scale = 0.9, nrow = 1,
                       labels = letters[c(1:(length(names(list_plots)) - 1))])
  
  # Save the plot
  ggsave(file.in, plot.out, width = 21, height = 7, units = "cm", 
         dpi = 600, bg = "white")
  
  # return the name of the file
  return(file.in)
}


#' Plot the effect of climate on functional strategy and diversity
#' @param data_models data with functional diversity, resilience metrics
#' @param file.in name including path of the file to save
plot_FD_and_CWM_vs_climate = function(data_models, file.in){
  
  # Create directory if needed
  create_dir_if_needed(file.in)
  
  # Format data
  data.in = data_models %>%
    filter(resistance > 0) %>%
    gather(key = "variable", value = "value",  
           "FD", "CWM") %>%
    mutate(ID.climate = factor(ID.climate, levels = unique(.$ID.climate))) %>%
    filter(is.finite(value))
  
  # Initialize list of plots
  list_plots = list()
  
  # Loop on all variables tested
  for(i in 1:length(unique(data.in$variable))){
    
    # Variable i
    var.i = unique(data.in$variable)[i]
    
    # Make model
    mod.i = aov(value ~ ID.climate, 
                data = subset(data.in, variable == var.i))
    
    # Data for plotting
    data.i = data.in %>%
      # Only keep the right variable
      filter(variable == var.i)
    
    # Data with the significance letters
    data.labels.i = data.i %>%
      # Calculate the 75% quantile
      group_by(pca1) %>%
      summarize(label.pos = max(value)) %>%
      # Add 10% of max value 
      mutate(label.pos = label.pos + 0.2*quantile(data.i$value, 0.975)) %>%
      left_join(
        data.frame(
          pca1 = unique(data.in$pca1), 
          label = as.character(cld(glht(mod.i, linfct=mcp(ID.climate="Tukey")))$mcletters$Letters)
        ), 
        by = "pca1") 
    
    # Make plot i
    plot.i = data.i %>%
      group_by(ID.climate, pca1) %>%
      summarize(mean = mean(value, na.rm = TRUE), 
                lwr = quantile(value, 0.025), 
                upr = quantile(value, 0.975)) %>%
      mutate(climate = paste0("climate.", ID.climate)) %>%
      ggplot(aes(x = pca1, y = mean, fill = climate)) +
      geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0) +
      geom_point(shape = 21, color = "black", size = 2) + 
      scale_fill_manual(
        values = colorRampPalette(c("blue", "red"))(length(unique(data.in$ID.climate)))) + 
      geom_text(data = (data.labels.i %>% mutate(climate = NA)),
                aes(label = label, y = label.pos), 
                size = 6, color = "black") + 
      xlab("CLIMATE") + ylab(toupper(var.i)) + 
      ggtitle(paste0(
        "F = ", round(as.numeric(summary(mod.i)[[1]][["F value"]][1]), digits = 1), 
        ", ", scales::pvalue(as.numeric(summary(mod.i)[[1]][["Pr(>F)"]][1]), 
                             accuracy = 0.01, add_p = TRUE)
      )) +
      theme(panel.background = element_rect(fill = "white", color = "black"), 
            panel.grid = element_blank(), 
            legend.position = "none", 
            legend.title = element_blank(), 
            axis.text.x = element_blank(), 
            axis.ticks.x = element_blank(), 
            legend.key = element_blank())
    
    # Add to the output lists
    eval(parse(text = paste0("list_plots$", var.i, " = plot.i")))
  }
  
  # Add legend to the output list
  list_plots$legend = get_legend(plot.i + theme(legend.position = "left", 
                                                legend.text = element_text(size = 16)))
  
  # Final plot
  plot.out = plot_grid(plotlist = list_plots, align = "hv", scale = 0.9, nrow = 1,
                       labels = c(letters[c(1:(length(names(list_plots)) - 1))], ""))
  
  # Save the plot
  ggsave(file.in, plot.out, width = 18, height = 6, units = "cm", 
         dpi = 600, bg = "white")
  
  # return the name of the file
  return(file.in)
  
}


#' Plot the effects of CWM and FD along climatic gradient
#' @param data_models data with functional diversity, resilience metrics
#' @param file.in name including path of the file to save
plot_CMW_and_FD_effect_climate = function(data_models, file.in){
  
  # Create directory if needed
  create_dir_if_needed(file.in)
  
  # Format data
  data.in = data_models %>%
    filter(resistance > 0) %>%
    mutate(CWM.scaled = scale(.$CWM), 
           FD.scaled = scale(.$FD)) %>%
    gather(key = "variable", value = "value",  
           "resistance", "resilience", "recovery") %>%
    filter(is.finite(value))
  
  # Initialize output dataset
  data.out = expand.grid(ID.climate = unique(data.in$ID.climate), 
                         variable = unique(data.in$variable), 
                         var.exp = c("CWM", "FD", "CWM:FD")) %>%
    mutate(climate = paste0("climate_", ID.climate)) %>%
    mutate(est = NA_real_, lwr = NA_real_, upr = NA_real_)
  
  # model to unscale CWM
  mod.unscale.CWM = lm(CWM ~ CWM.scaled, data = subset(data.in, var = "resistance"))
  mod.unscale.FD = lm(FD ~ FD.scaled, data = subset(data.in, var = "resistance"))
  
  
  # Loop on all variables tested
  for(i in 1:length(unique(data.in$variable))){
    
    # Variable i
    var.i = unique(data.in$variable)[i]
    
    # Loop on all climate tested
    for(j in 1:length(unique(data.in$ID.climate))){
      
      # Make model for variable i and climate j
      mod.ij = lm(log(value) ~ CWM.scaled*FD.scaled, 
                  data = subset(data.in, variable == var.i & ID.climate == j))
      
      # Store model output for each explanatory variable
      for(k in 1:3){
        
        # explanatory variable k
        var.exp.k = c("CWM", "FD", "CWM:FD")[k]
        
        # ID of climate j, var i and var exp k in data.out
        ID.ijk = which(data.out$ID.climate == j & 
                         data.out$variable == var.i & 
                         data.out$var.exp == var.exp.k)
        
        # Fill data.out with model output
        data.out$est[ID.ijk] = coef(summary(mod.ij))[(k+1), 1]
        data.out$lwr[ID.ijk] = confint(mod.ij)[(k+1), 1]
        data.out$upr[ID.ijk] = confint(mod.ij)[(k+1), 2]
      }
    }
  }
  
  # Final plot
  plot.out <- data.out %>%
    mutate(signif = ifelse(lwr > 0 | upr < 0, "yes", "no")) %>%
    left_join((data.in %>% dplyr::select(ID.climate, pca1) %>% distinct()), 
              by = "ID.climate") %>%
    ggplot(aes(x = pca1, y = est, color = signif, group = 1)) + 
    geom_point() + 
    geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0) + 
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray") + 
    scale_color_manual(values = c(`no` = "gray", `yes` = "black")) + 
    theme(panel.background = element_rect(fill = "white", color = "black"), 
          panel.grid = element_blank(), 
          legend.position = "none", 
          strip.background = element_blank()) + 
    xlab("PCA1 climate") + ylab("Effect") +
    facet_grid(var.exp ~ variable, scales = "free")
  
  # Save the plot
  ggsave(file.in, plot.out, width = 12, height = 9, units = "cm", 
         dpi = 600, bg = "white")
  
  # return the name of the file
  return(file.in)
  
}



#' Analyse the data with a structural equation model approach
#' @param data_model formatted model output
#' @param file.in name of the file to save, including path
plot_sem = function(data_model, file.in){
  
  # Create directory if needed
  create_dir_if_needed(file.in)
  
  # Make the sem model
  
  # -- Start by formatting data before fitting the model
  data.in = data_model %>%
    # Remove invalid data
    filter(resistance > 0) %>%
    # Log transform resilience metrics to fit normality assumption
    mutate(resistance.log = log(resistance), 
           recovery.log = log(recovery), 
           resilience.log = log(resilience)) %>%
    # scale all variables used in models
    mutate(climate_scaled = as.numeric(scale(pca1, center = TRUE, scale = TRUE)), 
           FD_scaled = as.numeric(scale(FD, center = TRUE, scale = TRUE)), 
           CWM_scaled = as.numeric(scale(CWM, center = TRUE, scale = TRUE)), 
           resistance.log_scaled = as.numeric(scale(resistance.log, center = TRUE, scale = TRUE)), 
           recovery.log_scaled = as.numeric(scale(recovery.log, center = TRUE, scale = TRUE)), 
           resilience.log_scaled = as.numeric(scale(resilience.log, center = TRUE, scale = TRUE)))
  
  # -- Make model
  mod_sem = psem(
    lm(FD_scaled ~ climate_scaled, data = data.in), 
    lm(CWM_scaled ~ climate_scaled, data = data.in), 
    lm(resistance.log_scaled ~ FD_scaled + CWM_scaled, data = data.in), 
    lm(recovery.log_scaled ~ FD_scaled + CWM_scaled + climate_scaled, data = data.in), 
    lm(resilience.log_scaled ~ resistance.log_scaled + recovery.log_scaled + 
         FD_scaled + CWM_scaled + climate_scaled, data = data.in)
  )
  
  
  
  
  
  # Plot the results
  
  # -- Prepare some parameters about the box to drax
  box.height = 2
  box.width = 4
  box.width.spacing = 2
  box.height.spacing = 2
  
  # -- Make a dataset to plot the box with the text
  data.plot.box = data.frame(
    text = c("climate", "FD", "CWM", "recovery", "resistance", "resilience"), 
    position = c("c", "l", "r", "l", "r", "c"), 
    height.level = c(4, 3, 3, 2, 2, 1)) %>%
    mutate(ymin = (height.level - 0.5*box.height)*(box.height + box.height.spacing), 
           ymax = ymin + box.height,
           xmin = case_when(position == "c" ~ -box.width/2, 
                            position == "l" ~ -(box.width + 0.5*box.width.spacing), 
                            position == "r" ~ box.width.spacing*0.5), 
           xmax = xmin + box.width,
           center.y = 0.5*(ymin + ymax), 
           center.x = 0.5*(xmin + xmax))
  
  # -- Plot the boxes
  plot.box = data.plot.box %>%
    ggplot(aes(x = center.x, y = center.y))  + 
    geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), 
              color = "black", fill = "#415A77")+ 
    geom_text(aes(label = text), color = "white", size = 7) + 
    theme(axis.title = element_blank(), 
          axis.ticks = element_blank(), 
          axis.text = element_blank(), 
          panel.background = element_rect(fill = "white", color = "black"), 
          legend.position = "none",
          panel.grid = element_blank())
  
  # -- Extract model results and prepare for plotting
  data.plot.arrow = as.data.frame(summary(mod_sem, .progressBar = F)$coefficients) %>%
    # Extract variable of interest
    dplyr::select(var.resp = Response, var.exp = Predictor, 
                  est = Estimate, p = P.Value) %>%
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
    # Difference between the two boxes connected by the arrow
    mutate(level.diff = level.exp - level.resp) %>%
    # Calculate intermediate height when level difference is 2
    mutate(height.mid = case_when(
      level.diff == 2 & level.exp == 4 ~ center.y.resp + 0.5*box.height + 0.5*box.height.spacing, 
      level.diff == 2 & level.exp == 3 ~ center.y.exp - 0.5*box.height - 0.7*box.height.spacing, 
      TRUE ~ NA_real_
    )) %>%
    # Slightly move the center x of response to make arrows visible
    mutate(center.x.resp = case_when(
      center.x.resp != center.x.exp ~ center.x.resp + 0.05*center.x.exp, 
      TRUE ~ center.x.resp
    )) %>%
    # Slightly move the center x of explanatory when level diff == 2 to be more visible
    mutate(center.x.exp = case_when(
      level.diff == 2 & level.exp == 3 ~ center.x.exp*0.97, 
      level.diff == 1 & level.exp == 3 ~ center.x.exp*1.05,
      TRUE ~ center.x.exp
    )) %>%
    # Report significance, sign and absolute value of estimate
    mutate(signif = ifelse(p <= 0.05, "yes", "no"), 
           sign = ifelse(est < 0, "negative", "positive"), 
           est.abs = abs(est), 
           est.abs.cat = case_when(
             est.abs <= quantile(est.abs, 0.33) ~ "low", 
             est.abs > quantile(est.abs, 0.33) & 
               est.abs <= quantile(est.abs, 0.66) ~ "mid", 
             TRUE ~ "high"
           )) 
  
  # -- plot the box with arrows
  plot.out = plot.box + 
    # direct arrows for difference in levels lower or higher than 2
    geom_segment(data = subset(data.plot.arrow, level.diff != 2), 
                 aes(x = center.x.exp, xend = center.x.resp, 
                     y = center.y.exp - 0.5*box.height, 
                     yend = center.y.resp + 0.5*box.height, 
                     size = est.abs.cat, linetype = signif, color = sign), 
                 arrow = arrow(length = unit(0.15, "cm")), 
                 type = "closed") + 
    # Segments and arrow for difference in levels of 2
    # -- first vertical segment
    geom_segment(data = subset(data.plot.arrow, level.diff == 2), 
                 aes(x = center.x.exp, xend = center.x.exp, 
                     y = center.y.exp - 0.5*box.height, 
                     yend = height.mid, 
                     size = est.abs.cat, linetype = signif, color = sign)) + 
    # -- second segment, horizontal
    geom_segment(data = subset(data.plot.arrow, level.diff == 2), 
                 aes(x = center.x.exp, xend = center.x.resp, 
                     y = height.mid, yend = height.mid, 
                     size = est.abs.cat, linetype = signif, color = sign)) + 
    # -- third: vertical arrow
    geom_segment(data = subset(data.plot.arrow, level.diff == 2), 
                 aes(x = center.x.resp, xend = center.x.resp, 
                     y = height.mid, 
                     yend = center.y.resp + 0.5*box.height, 
                     size = est.abs.cat, linetype = signif, color = sign), 
                 arrow = arrow(length = unit(0.15, "cm")), 
                 type = "closed") + 
    # -- scale linetype, color and size
    scale_color_manual(values = c(`negative` = "#0081A7", `positive` = "#F07167")) + 
    scale_linetype_manual(values = c(`yes` = "solid", `no` = "dashed")) + 
    scale_size_manual(values = c(`low` = 0.3, `mid` = 0.6, `high` = 1.1))
  
  # - Save the plot
  ggsave(file.in, plot.out, width = 10, height = 16, units = "cm", dpi = 600, bg = "white")
  
  # return the name of all the plots made
  return(file.in)
  
}







## -- Exploratory plots ----

#' Function to plot changes in CWM and FD along time till equilibrium
#' @param sim_equilibrium vector of the filenames of simulations saved as rds
#' @param forest_list dataset with information on all forests simulated
#' @param pc1_per_species df with coordinates of the climate pca per species
#' @param file.in Name of the file to save, including path
plot_cwm_fd_overtime = function(sim_equilibrium, forest_list, pc1_per_species, 
                                file.in){
  # Create the directories of file.in if needed
  create_dir_if_needed(file.in)
  
  # Loop on all simulations
  for(i in 1:length(sim_equilibrium)){
    
    # Printer
    print(paste0("Getting data for forest ", i, "/", length(sim_equilibrium)))
    
    # Extract the data for simulation i
    data.i = readRDS(sim_equilibrium[i]) %>%
      tree_format() %>%
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
      if(i == 1) data = data.i
      else data = rbind(data, data.i)
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
#' @param climate_list List containing quantiles for each climate
#' @param FUNDIV_climate_species df containing species presence and climate per FUNDIV plot 
#' @param disturbance.in disturbance to focus on ("storm", "biotic", "fire", or "snow")
#' @param exclude.in vector of additional species to exclude 
#' @param file.in Name of the file to save, including path
plot_prop.species_per_climate = function(
  climate_list, FUNDIV_climate_species, disturbance.in, exclude.in = c(), file.in){
  
  # Vector of all species
  species_vec = colnames(FUNDIV_climate_species)[grep("_", colnames(FUNDIV_climate_species))]
  
  # Vector of all species for which we have disturbance parameters
  data("disturb_coef")
  species_vec_dist = (disturb_coef %>%
                        filter(disturbance %in% disturbance.in))$species
  species_vec_dist = species_vec_dist[!(species_vec_dist %in% exclude.in)]
  
  # Loop on all climate
  for(i in 1:length(names(climate_list))){
    
    # Format data for climate i
    data.i = FUNDIV_climate_species %>%
      filter(pca1 > quantile(FUNDIV_climate_species$pca1, climate_list[[i]][1])) %>%
      filter(pca1 < quantile(FUNDIV_climate_species$pca1, climate_list[[i]][2])) %>%
      gather(key = "species", value = "present", species_vec) %>%
      group_by(species) %>%
      summarize(n = sum(present)) %>%
      mutate(present.in.data = ifelse(n == 0, 0, 1),
             present.in.dist = ifelse(species %in% species_vec_dist, 1, 0),
             climate = names(climate_list)[i], 
             pca1.mean = quantile(FUNDIV_climate_species$pca1, sum(climate_list[[i]])/2)) %>%
      dplyr::select(climate, pca1.mean, species, present.in.data, present.in.dist, n)
    
    # Add to the final dataset
    if(i == 1) data = data.i
    else data = rbind(data, data.i)
  }
  
  # Final plot 
  plot.out = data %>%
    ungroup() %>%
    group_by(climate, pca1.mean) %>%
    summarize(species.estimated = sum(present.in.data*present.in.dist)/sum(present.in.data), 
              trees.estimated = sum(present.in.dist*n)/sum(n)) %>%
    gather(key = "variable", value = "Proportion", "species.estimated", "trees.estimated") %>%
    ggplot(aes(x = pca1.mean, y = Proportion, fill = variable)) + 
    geom_point(shape = 21, color = "black") + 
    xlab("Position in the PCA climate") + ylim(0, 1) +
    ggtitle(disturbance.in) +
    theme(panel.background = element_rect(fill = "white", color = "black"), 
          panel.grid = element_blank(), 
          legend.title = element_blank(), 
          legend.key = element_blank())
  
  # - Save the plot
  ggsave(file.in, plot.out, width = 12, height = 8, units = "cm", dpi = 600, bg = "white")
  
  # return the name of all the plots made
  return(file.in)
  
}