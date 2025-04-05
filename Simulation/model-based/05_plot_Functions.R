

##### Bias ####
scenarioOverviewPlot <- function(Result_table, Dag = "A") {
  
  title <- paste0("m-DAG ", Dag) 
  
  Limits <- c(-180, 180)
  Breaks <- c(-150, -75, 0, 75, 150) 
  
  Result_table <- Result_table[(Result_table$mDag == Dag),]
  
  ylimits <- c("CC", "Ext", "Ext MCMI", "MI PMM", "MI Int", "MI CART", "MI RF", "MI Amelia")
  ylabels <- c("CC", "Ext", "Ext MCMI", "MI PMM", "MI Int", "MI CART", "MI RF", "MI Amelia")
  
  options(ggplot2.continuous.fill = scale_fill_distiller) # to recieve palette
  Result_table$RelBias <- round(Result_table$RelBias, digits = 0)
  
  custom_palette <-c('#d73027', '#fdae61', '#a6d96a',  '#006837', 
                     '#a6d96a', '#fdae61', '#d73027')
  custom_palette <- rev(custom_palette)
  
  plot <- ggplot(data = Result_table, aes(DGP, Method, fill = RelBias, label=RelBias)) + geom_tile(size=1.5, stat="identity", height=1, width=1) +
    scale_fill_gradientn(colours = custom_palette, limits = Limits, breaks = Breaks, oob = scales::squish) + geom_text(size= 3.0,hjust=0.5, vjust=0.2) +
    theme(plot.title = element_text(hjust = 0.5), legend.title = element_blank(), 
          panel.background = element_rect(fill = "white", colour = "white"),
          strip.background = element_rect(fill = 'grey')) + 
    ggtitle(title) +
    scale_y_discrete(limits = ylimits, labels = ylabels) +
    facet_wrap(~Scenario, nrow = 1)
  
  return(plot)
}



##### Coverage ####
scenarioOverviewPlotC <- function(Result_table, Dag = "A") {
  
  title <- paste0("m-DAG ", Dag) 
  
  Limits <- c(45, 100)
  Breaks <-c(45, 70, 95)
  
  Result_table <- Result_table[(Result_table$mDag == Dag),]
  ylimits <- c("CC", "Ext", "Ext MCMI", "MI PMM", "MI Int", "MI CART", "MI RF", "MI Amelia")
  ylabels <- c("CC", "Ext", "Ext MCMI", "MI PMM", "MI Int", "MI CART", "MI RF", "MI Amelia")
  
  # type = "viridis"
  options(ggplot2.continuous.fill = scale_fill_distiller) # to recieve palette
  Result_table$Coverage <- round(Result_table$Coverage*100, digits = 0)
  
  value_range <- c(37, 100) # Example range from 40 to 100
  colors <- c('#d73027', '#f46d43', '#fdae61', '#fee08b', '#d9ef8b', '#66bd63', '#006837', '#66bd63')
  values <- c(37, 50, 60, 70, 80, 90, 95, 100)
  normalized_values <- (values - min(values)) / (max(values) - min(values))
  
  plot <- ggplot(data = Result_table, aes(DGP, Method, fill = Coverage, label= Coverage)) + geom_tile(size=1.5, stat="identity", height=1, width=1) +
    scale_fill_gradientn(colors = colors, values = scales::rescale(normalized_values), limits = Limits, breaks = Breaks, oob = scales::squish) + geom_text(size= 3.0,hjust=0.5, vjust=0.2) +
    theme(plot.title = element_text(hjust = 0.5), legend.title = element_blank(),
          panel.background = element_rect(fill = "white", colour = "white"),
          strip.background = element_rect(fill = 'grey')) + ggtitle(title) +
    scale_y_discrete(limits = ylimits, labels = ylabels) +
    facet_wrap(~Scenario, nrow = 1)
  
  return(plot)
}


##### RMSE ####
scenarioOverviewPlotRMSE <- function(Result_table, Dag = "A") {
  
  title <- paste0("m-DAG ", Dag) 
  
  Limits <- c(0.05, 0.48)
  Breaks <- c(0.10, 0.20, 0.30, 0.40)
  
  
  Result_table <- Result_table[(Result_table$mDag == Dag),]
  
  ylimits <- c("CC", "Ext", "Ext MCMI", "MI PMM", "MI Int", "MI CART", "MI RF", "MI Amelia")
  ylabels <- c("CC", "Ext", "Ext MCMI", "MI PMM", "MI Int", "MI CART", "MI RF", "MI Amelia")
  
  options(ggplot2.continuous.fill = scale_fill_distiller) # to recieve palette
  Result_table$RMSE <- round(Result_table$RMSE, digits = 2)
  
  plot <- ggplot(data = Result_table, aes(DGP, Method, fill = RMSE, label=RMSE)) + geom_tile(size=1.5, stat="identity", height=1, width=1) +
    scale_fill_continuous(palette = "RdYlGn", limits = Limits, breaks = Breaks) + geom_text(size= 3.0,hjust=0.50, vjust=0.2) +
    theme(plot.title = element_text(hjust = 0.5), legend.title = element_blank(), 
          panel.background = element_rect(fill = "white", colour = "white"),
          strip.background = element_rect(fill = 'grey')) + 
    ggtitle(title) +
    scale_y_discrete(limits = ylimits, labels = ylabels) +
    facet_wrap(~Scenario, nrow = 1)
  
  
  
  return(plot)
}
