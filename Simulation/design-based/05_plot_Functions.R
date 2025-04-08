all_measures_TMLE <- function(results_list) {
  
  input_name <- deparse(substitute(results_list))
  
  calculate_measures <- function(TMLE_list) {
    
    TrueATE <-  0.1101056
    
    estimate_vector <- vector(length = (length(TMLE_list)))
    variance_values <- vector(length = (length(TMLE_list)))
    CI_upper_values <- vector(length = (length(TMLE_list)))
    CI_lower_values <- vector(length = (length(TMLE_list)))
    measure_vector <- vector(length = 10)
    
    for (i in 1:(length(TMLE_list))) {
      estimate_vector[i] <- TMLE_list[[i]][[1]]
      variance_values[i] <- TMLE_list[[i]][[2]]
      CI_upper_values[i] <- TMLE_list[[i]][[1]] + 1.96*(sqrt(TMLE_list[[i]][[2]]))
      CI_lower_values[i] <- TMLE_list[[i]][[1]] - 1.96*(sqrt(TMLE_list[[i]][[2]]))
    }
    
    measure_vector[1] <- mean(estimate_vector) #mean
    measure_vector[2] <- measure_vector[1]-TrueATE #bias
    measure_vector[3] <- (measure_vector[2]/TrueATE)*100 # rel bias in %
    measure_vector[4] <- sqrt((1/(length(estimate_vector)-1))*(sum((estimate_vector-measure_vector[1])^2))) #emp.SE
    measure_vector[5] <- sqrt((1/(length(estimate_vector)))*(sum((estimate_vector-TrueATE)^2))) #RMSE
    measure_vector[6] <- sqrt(mean(variance_values)) #Mod.SE
    measure_vector[7] <- 100*((measure_vector[6]/measure_vector[4] - 1)) #Relative error in Mod.SE in %
    measure_vector[8] <- mean(ifelse((CI_upper_values >= TrueATE) & 
                                       (CI_lower_values <= TrueATE), 1, 0 )) #Coverage
    measure_vector[9] <- mean(ifelse((CI_upper_values >= measure_vector[1]) & 
                                       (CI_lower_values <= measure_vector[1]), 1, 0 )) #Bias.eliminated.Coverage
    measure_vector[10] <- mean(abs(CI_upper_values-CI_lower_values)) # Mean.CI.Length
    
    return(measure_vector)
  }
  
  calculate_CIProp <- function(Result_table) {
    # create vector to store calculated Proportions
    CILengthProportion <- vector(length = nrow(Result_table))
    # for-loop for every method in results dataframe
    for(i in 1:nrow(Result_table)) {
      # divides respective average CI-length through maximal occurred CI-length
      CILengthProportion[i] <- (Result_table$MeanCILength[i])/(max(Result_table$MeanCILength)) 
    }
    Result_table <- cbind(Result_table, CILengthProportion)
    return(Result_table)
  }
  
  results <- list()
  
  for (method in names(results_list)) {
    method_results <- results_list[[method]]
    measure_vector <- calculate_measures(method_results)
    results[[method]] <- measure_vector
  }
  
  results_df <- do.call(rbind, lapply(results, function(x) as.data.frame(t(x))))
  rownames(results_df) <- names(results_list)
  colnames(results_df) <- c("mean", "Bias", "RelBias", "empSE", "RMSE", "ModSE", "RelErrorModSE", "Coverage", "BiasElimCoverage", "MeanCILength")
  
  results_df <- calculate_CIProp(results_df)
  
  
  # Extract additional columns from input_name
  Scenario <- if (grepl("pos1", input_name)) {
    "Scenario2"
  } else if (grepl("pos2", input_name)) {
    "Scenario3"
  } else if (grepl("pos3", input_name)) {
    "Scenario4"
  } else {
    "Scenario1"
  }
  mDag <- substring(input_name, 5, 5)
  recov <- ifelse(grepl("mDag[AB]", input_name), "yes", "no")
  
  
  results_df$Method <- gsub("_", " ", sub("_res$", "", rownames(results_df)))
  results_df$Scenario <- Scenario
  results_df$mDag <- mDag
  results_df$recov <- recov
  
  # Update rownames of results_df
  new_rownames <- sapply(rownames(results_df), function(name) {
    new_name <- sub("_res$", paste0("_sce", substring(Scenario, 9, 9)), name)
    paste0("mDag", mDag, "_", new_name)
  })
  rownames(results_df) <- new_rownames
  
  return(results_df)
}

# Bias Plot 
OverviewPlotBias <- function(Result_table) {
  
  #title <- paste0("m-DAG ", Dag) 
  title <- "Relative Bias in %"
  
  Limits <- c(-80, 80)
  Breaks <- c(-80, -40, 0, 40, 80) 
  
  
  Result_table <- Result_table[(Result_table$Method == "CC" | Result_table$Method == "Ext" | Result_table$Method == "Ext MCMI" |
                                  Result_table$Method == "MI PMM" | Result_table$Method == "MI Int" |
                                  Result_table$Method == "MI CART" | Result_table$Method == "MI RF" | Result_table$Method == "MI Amelia"),]
  ylimits <- c("CC", "Ext", "Ext MCMI", "MI PMM", "MI Int", "MI CART", "MI RF", "MI Amelia")
  ylabels <- c("CC", "Ext", "Ext MCMI", "MI PMM", "MI Int", "MI CART", "MI RF", "MI Amelia")
  
  options(ggplot2.continuous.fill = scale_fill_distiller) # to recieve palette
  Result_table$RelBias <- round(Result_table$RelBias, digits = 0)
  
  custom_palette <-c('#d73027', '#fdae61', '#a6d96a',  '#006837', 
                     '#a6d96a', '#fdae61', '#d73027')
  custom_palette <- rev(custom_palette)
  
  # Modify the facet label
  Result_table$FacetLabel <- paste0("m-DAG ", Result_table$mDag)
  Result_table$FacetLabel <- factor(Result_table$FacetLabel, levels = c("m-DAG A", "m-DAG B", "m-DAG C", "m-DAG D", "m-DAG E"))
  
  plot <- ggplot(data = Result_table, aes(Scenario, Method, fill = RelBias, label=RelBias)) + geom_tile(size=1.5, stat="identity", height=1, width=1) +
    scale_fill_gradientn(colours = custom_palette, limits = Limits, breaks = Breaks, oob = scales::squish) + geom_text(size= 3.0,hjust=0.5, vjust=0.2) +
    theme(plot.title = element_text(hjust = 0.5), legend.title = element_blank(), 
          panel.background = element_rect(fill = "white", colour = "white"),
          strip.background = element_rect(fill = 'grey')) + 
    ggtitle(title) +
    scale_y_discrete(limits = ylimits, labels = ylabels) +
    scale_x_discrete(labels = function(x) sub("Scenario", "", x)) + 
    facet_wrap(~FacetLabel, nrow = 1)
  
  return(plot)
}

# Coverage Plot
OverviewPlotCoverage <- function(Result_table) {
  
  #title <- paste0("m-DAG ", Dag) 
  title <- "Nominal Coverage in %"
  
  Limits <- c(60, 100)
  Breaks <-c(60, 70, 80, 90, 100)
  
  
    Result_table <- Result_table[(Result_table$Method == "CC" | Result_table$Method == "Ext" | Result_table$Method == "Ext MCMI" |
                                    Result_table$Method == "MI PMM" | Result_table$Method == "MI Int" |
                                    Result_table$Method == "MI CART" | Result_table$Method == "MI RF" | Result_table$Method == "MI Amelia"),]
    ylimits <- c("CC", "Ext", "Ext MCMI", "MI PMM", "MI Int", "MI CART", "MI RF", "MI Amelia")
    ylabels <- c("CC", "Ext", "Ext MCMI", "MI PMM", "MI Int", "MI CART", "MI RF", "MI Amelia")
  
  # type = "viridis"
  options(ggplot2.continuous.fill = scale_fill_distiller) # to recieve palette
  Result_table$Coverage <- round(Result_table$Coverage*100, digits = 0)
  
  value_range <- c(70, 100) # Example range from 40 to 100
  colors <- c('#d73027', '#f46d43', '#fdae61', '#fee08b', '#d9ef8b', '#66bd63', '#006837', '#66bd63')
  values <- c(60, 70, 75, 80, 85, 90, 95, 100)
  normalized_values <- (values - min(values)) / (max(values) - min(values))
  
  custom_palette <-c('#d73027', '#fdae61', '#a6d96a',  '#006837', 
                     '#a6d96a', '#fdae61', '#d73027')
  custom_palette <- rev(custom_palette)
  
  # Modify the facet label
  Result_table$FacetLabel <- paste0("m-DAG ", Result_table$mDag)
  Result_table$FacetLabel <- factor(Result_table$FacetLabel, levels = c("m-DAG A", "m-DAG B", "m-DAG C", "m-DAG D", "m-DAG E"))
  
  plot <- ggplot(data = Result_table, aes(Scenario, Method, fill = Coverage, label= Coverage)) + geom_tile(size=1.5, stat="identity", height=1, width=1) +
    scale_fill_gradientn(colors = colors, values = scales::rescale(normalized_values), limits = Limits, breaks = Breaks, oob = scales::squish) + geom_text(size= 3.0,hjust=0.5, vjust=0.2) +
    theme(plot.title = element_text(hjust = 0.5), legend.title = element_blank(),
          panel.background = element_rect(fill = "white", colour = "white"),
          strip.background = element_rect(fill = 'grey')) + ggtitle(title) +
    ggtitle(title) +
    scale_y_discrete(limits = ylimits, labels = ylabels) +
    scale_x_discrete(labels = function(x) sub("Scenario", "", x)) + 
    facet_wrap(~FacetLabel, nrow = 1)
  
  return(plot)
}

# RMSE Plot
OverviewPlotRMSE <- function(Result_table) {
  
  #title <- paste0("m-DAG ", Dag) 
  title <- "Root Mean Squared Error"
  
  Limits <- c(0.001, 0.110)
  Breaks <- c(0.001, 0.05, 0.10)
  
    Result_table <- Result_table[(Result_table$Method == "CC" | Result_table$Method == "Ext" | Result_table$Method == "Ext MCMI" |
                                    Result_table$Method == "MI PMM" | Result_table$Method == "MI Int" |
                                    Result_table$Method == "MI CART" | Result_table$Method == "MI RF" | Result_table$Method == "MI Amelia"),]
    ylimits <- c("CC", "Ext", "Ext MCMI", "MI PMM", "MI Int", "MI CART", "MI RF", "MI Amelia")
    ylabels <- c("CC", "Ext", "Ext MCMI", "MI PMM", "MI Int", "MI CART", "MI RF", "MI Amelia")
  
  options(ggplot2.continuous.fill = scale_fill_distiller) # to recieve palette
  Result_table$RMSE <- round(Result_table$RMSE, digits = 2)
  
  custom_palette <-c('#d73027', '#fdae61', '#a6d96a',  '#006837', 
                     '#a6d96a', '#fdae61', '#d73027')
  custom_palette <- rev(custom_palette)
  
  # Modify the facet label
  Result_table$FacetLabel <- paste0("m-DAG ", Result_table$mDag)
  Result_table$FacetLabel <- factor(Result_table$FacetLabel, levels = c("m-DAG A", "m-DAG B", "m-DAG C", "m-DAG D", "m-DAG E"))
  
  plot <- ggplot(data = Result_table, aes(Scenario, Method, fill = RMSE, label=RMSE)) + geom_tile(size=1.5, stat="identity", height=1, width=1) +
    scale_fill_gradientn(colours = custom_palette, limits = Limits, breaks = Breaks, oob = scales::squish) + geom_text(size= 3.0,hjust=0.5, vjust=0.2) +
    theme(plot.title = element_text(hjust = 0.5), legend.title = element_blank(), 
          panel.background = element_rect(fill = "white", colour = "white"),
          strip.background = element_rect(fill = 'grey')) + 
    ggtitle(title) +
    scale_y_discrete(limits = ylimits, labels = ylabels) +
    scale_x_discrete(labels = function(x) sub("Scenario", "", x)) + 
    facet_wrap(~FacetLabel, nrow = 1)
  
  return(plot)
}


