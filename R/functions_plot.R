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
    mod.i = aov(log(value) ~ ID.climate, 
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
      mutate(value.log = log(value)) %>%
      group_by(ID.climate, pca1) %>%
      summarize(mean.log = mean(value.log, na.rm = TRUE), 
                lwr.log = quantile(value.log, 0.025), 
                upr.log = quantile(value.log, 0.975)) %>%
      mutate(climate = paste0("climate.", ID.climate), 
             value = exp(mean.log), 
             lwr = exp(lwr.log), 
             upr = exp(upr.log)) %>%
      ggplot(aes(x = pca1, y = value, fill = climate)) +
      geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0) +
      geom_point(shape = 21, color = "black") + 
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
            axis.ticks.x = element_blank())
    
    # Add to the output lists
    eval(parse(text = paste0("list_plots$", var.i, " = plot.i")))
  }
  
  # Add legend to the output list
  list_plots$legend = get_legend(plot.i + theme(legend.position = "left", 
                                                legend.text = element_text(size = 16)))
  
  # Final plot
  plot.out = plot_grid(plotlist = list_plots, align = "hv", scale = 0.9,
                       labels = c(letters[c(1:(length(names(list_plots)) - 1))], ""))
  
  # Save the plot
  ggsave(file.in, plot.out, width = 16, height = 16, units = "cm", 
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

