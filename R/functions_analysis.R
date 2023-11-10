#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#### SCRIPT INTRODUCTION ####
#
#' @name functions_analysis.R  
#' @description R script containing all functions relative to data
#               analysis
#' @author Julien Barrere
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



#' Fit SEM model and bootstrap for direct and indirect effect 
#' @param data_model formatted model output
#' @param FD_metric Functional diversity metric to choose ("FDis", "FRic or "FD")
#' @param recovery_metric Recovery metric to choose ("recovery", "thalf)
#' @param R_metric Richness metric to choose ("nsp", "H" or "D")
fit_sem = function(data_model, FD_metric = "FDis", R_metric = "H", 
                   recovery_metric = "recovery"){
  
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
  mod_sem = list(
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
  
  
  
  # Extract the effects
  boot_sem = bootEff(mod_sem, R = 1000, seed = 13, parallel = "no")
  
  # Output list
  out = list()
  out$sem = mod_sem
  out$boot = boot_sem
  
  # Return output
  return(out)
}



#' Fit SEM model and bootstrap for direct and indirect effect (supp info sem)
#' @param data_model formatted model output
#' @param FD_metric Functional diversity metric to choose ("FDis", "FRic or "FD")
#' @param recovery_metric Recovery metric to choose ("recovery", "thalf)
#' @param R_metric Richness metric to choose ("nsp", "H" or "D")
fit_sem_supp = function(data_model, FD_metric = "FDis", R_metric = "H", 
                   recovery_metric = "recovery"){
  
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
  mod_sem = list(
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
  
  
  
  # Extract the effects
  boot_sem = bootEff(mod_sem, R = 1000, seed = 13, parallel = "no")
  
  # Output list
  out = list()
  out$sem = mod_sem
  out$boot = boot_sem
  
  # Return output
  return(out)
}


#' Table to make a statistics table from boot of sem
#' @param boot.in bootstrapping made by bootEff function of semEff package
make_sem_stat = function(boot.in){
  
  # Direct and indirect effects interpreted by semEff function from bootstrap
  eff.in = semEff(boot.in)
  
  # Vector of response variables
  var.resp.in = gsub("\\_", "\\.", names(boot.in))
  
  # Number of bootstrapping iterations
  R.in = boot.in[[1]]$R
  
  # Loop on all response variables
  for(i in 1:length(var.resp.in)){
    
    # Variable i
    var.i = var.resp.in[i]
    
    # Start by direct effects
    data.i = data.frame(
      response = var.i, type = "direct",
      explanatory = rep(names(eff.in[[2]][[var.i]]$Direct), each = R.in),
      effect = as.numeric(eff.in[[3]][[var.i]]$Direct)) %>%
      # Add total effects
      rbind(data.frame(
        response = var.i, type = "total",
        explanatory = rep(names(eff.in[[2]][[var.i]]$Total), each = R.in),
        effect = as.numeric(eff.in[[3]][[var.i]]$Total)))
    
    # If there are indirect effects, add indirect effects
    if(length(names(eff.in[[2]][[var.i]]$Indirect)) > 0){
      data.i = data.i %>%
        rbind(data.frame(
          response = var.i, type = "indirect",
          explanatory = rep(names(eff.in[[2]][[var.i]]$Indirect), each = R.in),
          effect = as.numeric(eff.in[[3]][[var.i]]$Indirect)))
    }
    
    if(i == 1) data = data.i
    else data = rbind(data, data.i)
    
  }
  
  
  # Finish data formatting
  data.out = data %>%
    # Remove intercepts
    filter(explanatory != "(Intercept)") %>%
    # Calculate summary statistics
    group_by(response, explanatory, type) %>%
    summarize(eff = mean(effect, na.rm = TRUE), 
              lwr = quantile(effect, 0.025, na.rm = TRUE), 
              upr = quantile(effect, 0.975, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(signif = ifelse(lwr*upr > 0, "*", "")) %>%
    # Format variable names
    mutate(response = gsub("\\..+", "", response), 
           explanatory = gsub("\\..+", "", explanatory)) %>%
    # Round numbers
    mutate(eff = round(eff, digits = 3), lwr = round(lwr, digits = 3), 
           upr = round(upr, digits = 3)) %>%
    # Text to print in the final table
    mutate(text = paste0(eff, " (", lwr, "-", upr, ")", signif)) %>%
    # Final formatting
    dplyr::select(response, explanatory, type, text) %>%
    spread(key = "type", value = "text") %>%
    mutate(indirect = ifelse(is.na(indirect), "", indirect)) %>%
    mutate(response = factor(response, levels = c(
      "H", "FD", "CWM1", "CWM2", "resistance", "recovery", "resilience"))) %>%
    arrange(response) %>%
    mutate(response = as.character(response))
  
  # Make the table ready to export as latex file
  out = data.frame(
    col1 = c("Response", "variable", data.out$response), 
    col2 = c("Explanatory", "variable", data.out$explanatory), 
    col11 = "",
    col3 = c("Direct effect", "Est. (95% CI)", data.out$direct),
    col4 = "",
    col5 = c("Indirect effect", "Est. (95% CI)", data.out$indirect),
    col6 = "",
    col7 = c("Total effect", "Est. (95% CI)", data.out$total)
  )
  
  # Return output
  return(out)
}
