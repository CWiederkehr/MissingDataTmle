



check_data <- function(data_list) {
  # Initialize lists to store results
  Z1_means <- Z2_means <- Z3_props <- Z3_means <- Z4_means <- Z4_props <- Z4_means_sd <- Z5_means <- Z5_props <- Z5_means_sd <- NULL
  Z6_moments <- A_means <- Y_means <- Y_means_sd <- NULL
  
  # Check for the presence of each variable and calculate the required statistics
  if ("Z1" %in% names(data_list[[1]])) {
    Z1_means <- sapply(data_list, function(x) mean(x$Z1))
  }
  
  if ("Z2" %in% names(data_list[[1]])) {
    Z2_means <- sapply(data_list, function(x) mean(x$Z2))
  }
  
  if ("Z3" %in% names(data_list[[1]])) {
    Z3_props <- sapply(data_list, function(x) {
      if (is.factor(x$Z3)) {
        table(x$Z3) / nrow(x)
      } else {
        mean(x$Z3)
      }
    })
  }
  
  if ("Z4" %in% names(data_list[[1]])) {
    Z4_means <- sapply(data_list, function(x) {
      if (all(x$Z4 %in% c(0, 1))) {
        mean(x$Z4)
      } else {
        mean(x$Z4)
      }
    })
    Z4_means_sd <- sapply(data_list, function(x) {
      if (!all(x$Z4 %in% c(0, 1))) {
        sd(x$Z4)
      } else {
        NA
      }
    })
  }  
  
  if ("Z5" %in% names(data_list[[1]])) {
    Z5_means <- sapply(data_list, function(x) {
      if (all(x$Z5 %in% c(0, 1))) {
        mean(x$Z5)
      } else {
        mean(x$Z5)
      }
    })
    Z5_means_sd <- sapply(data_list, function(x) {
      if (!all(x$Z5 %in% c(0, 1))) {
        sd(x$Z5)
      } else {
        NA
      }
    })
  }
  
  if ("Z6" %in% names(data_list[[1]])) {
    Z6_moments <- sapply(data_list, function(x) {
      mean_Z6 <- mean(x$Z6)
      var_Z6 <- var(x$Z6)
      shape_est <- (mean_Z6^2) / var_Z6
      rate_est <- mean_Z6 / var_Z6
      c(shape_est, rate_est)
    })
  }
  
  if ("A" %in% names(data_list[[1]])) {
    A_means <- sapply(data_list, function(x) mean(x$A))
  }
  
  if ("Y" %in% names(data_list[[1]])) {
    Y_means <- sapply(data_list, function(x) mean(x$Y))
    Y_means_sd <- sapply(data_list, function(x) sd(x$Y))
  }
  
  # Combine results into a proportion vector
  Proportion <- c(
    if (!is.null(Z1_means)) mean(Z1_means) else NA,
    if (!is.null(Z2_means)) mean(Z2_means) else NA,
    if (!is.null(Z3_props) && is.matrix(Z3_props)) rowMeans(Z3_props) else if (is.numeric(Z3_props)) mean(Z3_props) else NA,
    if (!is.null(Z4_means)) mean(Z4_means) else NA,
    if (!is.na(Z4_means_sd[1])) mean(Z4_means_sd, na.rm = TRUE) else NA,
    if (!is.null(Z5_means)) mean(Z5_means) else NA,
    if (!is.na(Z5_means_sd[1])) mean(Z5_means_sd, na.rm = TRUE) else NA,
    if (!is.null(Z6_moments)) mean(Z6_moments[1,], na.rm = TRUE) else NA,
    if (!is.null(Z6_moments)) mean(Z6_moments[2,], na.rm = TRUE) else NA,
    if (!is.null(A_means)) mean(A_means) else NA,
    if (!is.null(Y_means)) mean(Y_means) else NA,
    if (!is.null(Y_means_sd)) mean(Y_means_sd) else NA
  )
  
  # Name the proportion vector
  names(Proportion) <- c(
    "Z1", "Z2", if (!is.null(Z3_props) && is.matrix(Z3_props)) paste0("Z3.", 1:nrow(Z3_props)) else "Z3",
    "Z4", "Z4.SD","Z5", "Z5.SD","Z6_shape", "Z6_rate", "A", "Y", "Y.SD"
  )
  
  return(Proportion)
}

check_positivity <- function(data, include_Z6 = FALSE) {
  # Define the propensity score model formula based on include_Z6
  if (include_Z6) {
    formula <- A ~ Z1 + Z2 + Z3 + Z4 + Z5 + Z6 + B + 
      Z1 : Z3 + Z1 : Z4 + Z1 : Z5 + Z3 : Z4 + Z3 : Z5 + Z4 : Z5 + 
      Z1 : Z6 + Z4 : Z6 + Z5 : Z6
  } else {
    formula <- A ~ Z1 + Z2 + Z3 + Z4 + Z5 + B + 
      Z1 : Z3 + Z1 : Z4 + Z1 : Z5 + Z3 : Z4 + Z3 : Z5 + Z4 : Z5
  }
  
  # Start timing
  start_time <- Sys.time()
  
  total_cores <- 4
  cores_per_task <- 2
  num_clusters <- total_cores / cores_per_task
  
  # Set up parallel backend
  cl <- makeCluster(num_clusters)
  registerDoParallel(cl, cores_per_task)
  
  # Create respective lists for models and weights using foreach
  results <- foreach(i = 1:length(data), .packages = c('glmnet')) %dopar% {
    psm <- glm(formula, family = binomial, data = data[[i]]) # Propensity-score model
    gW <- predict(psm, type = "response") # Propensity score values
    count <- ifelse(any(gW < 0.002 | (1 - gW) < 0.002), 1, 0) # Check if probability in dataset is smaller than 0.002
    sum_obs <- sum(ifelse(gW < 0.002 | (1 - gW) < 0.002, 1, 0)) # Count number of observations with smaller probability
    return(list(gW = gW, count = count, sum_obs = sum_obs))
    gc()
  }
  
  # Stop parallel backend
  stopCluster(cl)
  
  # End timing
  end_time <- Sys.time()
  print(end_time - start_time)
  
  # Extract results
  gW_list <- lapply(results, `[[`, 'gW')
  count_vec <- sapply(results, `[[`, 'count')
  sum_vec <- sapply(results, `[[`, 'sum_obs')
  
  # Random sample model
  sample_nr <- sample(x = 1:length(data), size = 1)
  psm <- glm(formula, family = binomial, data = data[[sample_nr]])
  
  # Create title based on pattern in the name of the input object "data"
  if (grepl("sce1", deparse(substitute(data)))) {
    title <- "Scenario 1"
  } else if (grepl("sce2", deparse(substitute(data)))) {
    title <- "Scenario 2"
  } else if (grepl("sce3", deparse(substitute(data)))) {
    title <- "Scenario 3"
  } else {
    title <- "Unknown Scenario"
  }
  
  gW <- predict(psm, type = "response") # Propensity score values
  df <- cbind(data[[sample_nr]], gW) # Create df for ggplot2
  df$A <- as.factor(ifelse(df$A == 0, "Control", "Treatment"))
  
  positivity <- vector("list", length = 5)
  positivity[[1]] <- paste("Proportion of datasets with violations:", sum(count_vec) / length(data))
  positivity[[2]] <- paste("Percentage of violations in each dataset:", (mean(sum_vec)/nrow(data[[1]]))*100, "%" )
  positivity[[3]] <- summary(as.vector(unlist(gW_list)))
  positivity[[4]] <- ggplot(data = df, aes(x = gW, color = A)) + 
    geom_density() + 
    scale_x_continuous(name = "Propensity Score P(A=1|W) and P(A=0|W)", limits = c(0, 1), breaks = c(0, 0.25, 0.50, 0.75, 1), position = 'bottom') +
    scale_y_continuous(name = "Probability density function") + 
    theme(legend.position = "right", plot.title = element_text(hjust = 0.5)) +
    ggtitle(title)
  
  names(positivity) <- c("DatasetViolations", "ProportionPerDataset", "OverallDistribution", "RandomPSDensitiy", "ComputationTime")
  
  positivity[[5]] <- end_time - start_time
  
  return(positivity)
}

# Generate function to analyze proportion of missingness
MissingProportion <- function(data_list) {
  
  # Check if Z6 is present in the first dataset
  has_Z6 <- "Z6" %in% colnames(data_list[[1]])
  
  # Initialize vectors for missingness proportions
  MissingZ2 <- vector(length = length(data_list))
  MissingZ3 <- vector(length = length(data_list))
  MissingZ4 <- vector(length = length(data_list))
  MissingA <- vector(length = length(data_list))
  MissingY <- vector(length = length(data_list))
  MissingAorY <- vector(length = length(data_list))
  MissingAny <- vector(length = length(data_list))
  
  if (has_Z6) {
    MissingZ6 <- vector(length = length(data_list))
  }
  
  for (i in 1:(length(data_list))) {
    MissingZ2[i] <- sum(is.na(data_list[[i]]$Z2))/nrow(data_list[[i]])
    MissingZ3[i] <- sum(is.na(data_list[[i]]$Z3))/nrow(data_list[[i]])
    MissingZ4[i] <- sum(is.na(data_list[[i]]$Z4))/nrow(data_list[[i]])
    if (has_Z6) {
      MissingZ6[i] <- sum(is.na(data_list[[i]]$Z6))/nrow(data_list[[i]])
    }
    MissingA[i] <- sum(is.na(data_list[[i]]$A))/nrow(data_list[[i]])
    MissingY[i] <- sum(is.na(data_list[[i]]$Y))/nrow(data_list[[i]])
    MissingAorY[i] <- sum(is.na(data_list[[i]]$A) | is.na(data_list[[i]]$Y))/nrow(data_list[[i]])
    if (has_Z6) {
      MissingAny[i] <- sum(is.na(data_list[[i]]$Z2) | is.na(data_list[[i]]$Z3) | is.na(data_list[[i]]$Z4) |
                             is.na(data_list[[i]]$Z6) |
                             is.na(data_list[[i]]$A) | is.na(data_list[[i]]$Y))/nrow(data_list[[i]])
    } else {
      MissingAny[i] <- sum(is.na(data_list[[i]]$Z2) | is.na(data_list[[i]]$Z3) | is.na(data_list[[i]]$Z4) |
                             is.na(data_list[[i]]$A) | is.na(data_list[[i]]$Y))/nrow(data_list[[i]])
    }
  }
  
  # Calculate mean missingness
  if (has_Z6) {
    mean_missing <- c(mean(MissingZ2), mean(MissingZ3), mean(MissingZ4), mean(MissingZ6),
                      mean(MissingA), mean(MissingY), mean(MissingAorY), mean(MissingAny))
    row.names <- c("Z2", "Z3", "Z4", "Z6", "A", "Y", "AorY", "Any")
  } else {
    mean_missing <- c(mean(MissingZ2), mean(MissingZ3), mean(MissingZ4),
                      mean(MissingA), mean(MissingY), mean(MissingAorY), mean(MissingAny))
    row.names <- c("Z2", "Z3", "Z4", "A", "Y", "AorY", "Any")
  }
  
  mean_missing <- as.data.frame(mean_missing)
  row.names(mean_missing) <- row.names
  return(t(mean_missing))
}

calculate_CIProp <- function(Result_table) {
  # create vector to store calculated Proportions
  CILengthProportion <- vector(length = nrow(Result_table))
  # for-loop for every method in results dataframe
  for(i in 1:nrow(Result_table)) {
    # divedes respective average CI-length through maximal occured CI-length
    CILengthProportion[i] <- (Result_table$MeanCILength[i])/(max(Result_table$MeanCILength)) 
  }
  Result_table <- cbind(Result_table, CILengthProportion)
  return(Result_table)
}

all_measures_TMLE <- function(results_list) {
  
  input_name <- deparse(substitute(results_list))
  
  calculate_measures <- function(TMLE_list, True_ATE) {
    
    estimate_vector <- vector(length = length(TMLE_list))
    variance_values <- vector(length = length(TMLE_list))
    CI_upper_values <- vector(length = length(TMLE_list))
    CI_lower_values <- vector(length = length(TMLE_list))
    measure_vector <- vector(length = 10)
    
    for (i in seq_along(TMLE_list)) {
      estimate_vector[i] <- TMLE_list[[i]][[1]]
      variance_values[i] <- TMLE_list[[i]][[2]]
      CI_upper_values[i] <- estimate_vector[i] + 1.96 * sqrt(variance_values[i])
      CI_lower_values[i] <- estimate_vector[i] - 1.96 * sqrt(variance_values[i])
    }
    
    measure_vector[1] <- mean(estimate_vector)  # Mean
    measure_vector[2] <- measure_vector[1] - True_ATE  # Bias
    measure_vector[3] <- (measure_vector[2] / True_ATE) * 100  # Relative Bias (%)
    measure_vector[4] <- sd(estimate_vector)  # Empirical SE
    measure_vector[5] <- sqrt(mean((estimate_vector - True_ATE)^2))  # RMSE
    measure_vector[6] <- sqrt(mean(variance_values))  # Model-based SE
    measure_vector[7] <- 100 * ((measure_vector[6] / measure_vector[4]) - 1)  # Relative Error in SE (%)
    measure_vector[8] <- mean((CI_lower_values <= True_ATE) & (CI_upper_values >= True_ATE))  # Coverage
    measure_vector[9] <- mean((CI_lower_values <= measure_vector[1]) & (CI_upper_values >= measure_vector[1]))  # Bias-eliminated Coverage
    measure_vector[10] <- mean(CI_upper_values - CI_lower_values)  # Mean CI Length
    
    return(measure_vector)
  }
  
  calculate_CIProp <- function(Result_table) {
    CILengthProportion <- Result_table$MeanCILength / max(Result_table$MeanCILength)
    Result_table <- cbind(Result_table, CILengthProportion)
    return(Result_table)
  }
  
  # Extract DGP and Scenario numbers from input_name
  DGP <- as.numeric(sub(".*_DGP(\\d+).*", "\\1", input_name))
  Scenario <- as.numeric(sub(".*_sce(\\d+).*", "\\1", input_name))
  
  # Assign True_ATE based on DGP and Scenario
  # Replace the following with your actual True_ATE values
  if (DGP == 1 && Scenario == 1) {
    True_ATE <- 0.20
  } else if (DGP == 1 && Scenario == 2) {
    True_ATE <- 0.20
  } else if (DGP == 1 && Scenario == 3) {
    True_ATE <- 0.24
  } else if (DGP == 2 && Scenario == 1) {
    True_ATE <- 0.18
  } else if (DGP == 2 && Scenario == 2) {
    True_ATE <- 0.20
  } else if (DGP == 2 && Scenario == 3) {
    True_ATE <- 0.23
  } else if (DGP == 3 && Scenario == 1) {
    True_ATE <- 0.18
  } else if (DGP == 3 && Scenario == 2) {
    True_ATE <- 0.19
  } else if (DGP == 3 && Scenario == 3) {
    True_ATE <- 0.22
  } else if (DGP == 4 && Scenario == 1) {
    True_ATE <- 0.18
  } else if (DGP == 4 && Scenario == 2) {
    True_ATE <- 0.20
  } else if (DGP == 4 && Scenario == 3) {
    True_ATE <- 0.22
  } else if (DGP == 5 && Scenario == 1) {
    True_ATE <- 0.18
  } else if (DGP == 5 && Scenario == 2) {
    True_ATE <- 0.20
  } else if (DGP == 5 && Scenario == 3) {
    True_ATE <- 0.22
  } else {
    stop("Unknown DGP and Scenario combination. Please define True_ATE for this combination.")
  }
  
  results <- list()
  
  for (method in names(results_list)) {
    method_results <- results_list[[method]]
    measure_vector <- calculate_measures(method_results, True_ATE)
    results[[method]] <- measure_vector
  }
  
  results_df <- do.call(rbind, lapply(results, function(x) as.data.frame(t(x))))
  rownames(results_df) <- names(results_list)
  colnames(results_df) <- c("Mean", "Bias", "RelBias", "EmpSE", "RMSE", "ModSE", "RelErrorModSE", "Coverage", "BiasElimCoverage", "MeanCILength")
  
  results_df <- calculate_CIProp(results_df)
  
  # Extract additional columns from input_name
  mDag <- substring(input_name, 5, 5)
  recov <- ifelse(grepl("mDag[AET]", input_name), "yes", "no")
  
  results_df$Method <- gsub("_", " ", sub("_res$", "", rownames(results_df)))
  results_df$Scenario <- paste0("Scenario", Scenario)
  results_df$mDag <- mDag
  results_df$recov <- recov
  results_df$DGP <- paste0("DGP", DGP)
  
  # Update rownames of results_df
  new_rownames <- sapply(rownames(results_df), function(name) {
    scenario <- sub(".*_sce(\\d+).*", "\\1", input_name)
    dgp <- sub(".*_DGP(\\d+).*", "\\1", input_name)
    new_name <- sub("_res$", paste0("_sce", scenario, "_DGP", dgp), name)
    paste0("mDag", mDag, "_", new_name)
  })
  rownames(results_df) <- new_rownames
  
  return(results_df)
}