
source("~/GitHub/MissingDataTmle/Simulation/model-based/01_used_libraries.R")
source("~/GitHub/MissingDataTmle/Simulation/model-based/Data_and_Missingness.R")
source("~/GitHub/MissingDataTmle/Simulation/model-based/general_Functions.R")
source("~/GitHub/MissingDataTmle/Simulation/model-based/analysis_Functions.R")

# Simulate datasets for DGP1, DGP2, DGP3, DGP4, and DGP5
DGP1_data <- list(
  sce1_DGP1 = simulateDatasets("sce1", "DGP1"),
  sce2_DGP1 = simulateDatasets("sce2", "DGP1"),
  sce3_DGP1 = simulateDatasets("sce3", "DGP1")
)

DGP2_data <- list(
  sce1_DGP2 = simulateDatasets("sce1", "DGP2"),
  sce2_DGP2 = simulateDatasets("sce2", "DGP2"),
  sce3_DGP2 = simulateDatasets("sce3", "DGP2")
)

DGP3_data <- list(
  sce1_DGP3 = simulateDatasets("sce1", "DGP3"),
  sce2_DGP3 = simulateDatasets("sce2", "DGP3"),
  sce3_DGP3 = simulateDatasets("sce3", "DGP3")
)

DGP4_data <- list(
  sce1_DGP4 = simulateDatasets("sce1", "DGP4"),
  sce2_DGP4 = simulateDatasets("sce2", "DGP4"),
  sce3_DGP4 = simulateDatasets("sce3", "DGP4")
)

DGP5_data <- list(
  sce1_DGP5 = simulateDatasets("sce1", "DGP5"),
  sce2_DGP5 = simulateDatasets("sce2", "DGP5"),
  sce3_DGP5 = simulateDatasets("sce3", "DGP5")
)


check_data(DGP1_data$sce1_DGP1) 
check_data(DGP1_data$sce2_DGP1) 
check_data(DGP1_data$sce3_DGP1) 

check_data(DGP2_data$sce1_DGP2) 
check_data(DGP2_data$sce2_DGP2) 
check_data(DGP2_data$sce3_DGP2) 

check_data(DGP3_data$sce1_DGP3)
check_data(DGP3_data$sce2_DGP3)
check_data(DGP3_data$sce3_DGP3)

check_data(DGP4_data$sce1_DGP4)
check_data(DGP4_data$sce2_DGP4)
check_data(DGP4_data$sce3_DGP4)

check_data(DGP5_data$sce1_DGP5)
check_data(DGP5_data$sce2_DGP5)
check_data(DGP5_data$sce3_DGP5)

#### Full Data Assesment ####

calculate_measures <- function(TMLE_list, True_ATE = 0.20) {
  
  estimate_vector <- vector(length = (length(TMLE_list)))
  variance_values <- vector(length = (length(TMLE_list)))
  CI_upper_values <- vector(length = (length(TMLE_list)))
  CI_lower_values <- vector(length = (length(TMLE_list)))
  measure_vector <- vector(length = 11)
  
  for (i in 1:(length(TMLE_list))) {
    estimate_vector[i] <- TMLE_list[[i]][[1]]
    variance_values[i] <- TMLE_list[[i]][[2]]
    CI_upper_values[i] <- TMLE_list[[i]][[1]] + 1.96*(sqrt(TMLE_list[[i]][[2]]))
    CI_lower_values[i] <- TMLE_list[[i]][[1]] - 1.96*(sqrt(TMLE_list[[i]][[2]]))
  }
  
  measure_vector[1] <- mean(estimate_vector) # mean
  measure_vector[2] <- measure_vector[1] - True_ATE # bias
  measure_vector[3] <- (measure_vector[2] / True_ATE) * 100 # rel bias in %
  measure_vector[4] <- sqrt((1 / (length(estimate_vector) - 1)) * (sum((estimate_vector - measure_vector[1])^2))) # emp.SE
  measure_vector[5] <- sqrt((1 / (length(estimate_vector))) * (sum((estimate_vector - True_ATE)^2))) # RMSE
  measure_vector[6] <- sqrt(mean(variance_values)) # Mod.SE
  measure_vector[7] <- 100 * ((measure_vector[6] / measure_vector[4] - 1)) # Relative error in Mod.SE in %
  measure_vector[8] <- mean(ifelse((CI_upper_values >= True_ATE) & (CI_lower_values <= True_ATE), 1, 0)) # Coverage
  measure_vector[9] <- mean(ifelse((CI_upper_values >= measure_vector[1]) & (CI_lower_values <= measure_vector[1]), 1, 0)) # Bias.eliminated.Coverage
  measure_vector[10] <- mean(abs(CI_upper_values - CI_lower_values)) # Mean.CI.Length
  measure_vector[11] <- mean(ifelse((CI_upper_values >= 0) & (CI_lower_values <= 0), 0, 1)) # effect.Power
  
  names(measure_vector) <- c("mean", "Bias", "relBias.%", "empSE", "RMSE", "ModSE", "ErrorModSE.%", "Coverage", "BiasElimCoverage", "MeanCILength", "effect.Power")
  
  return(measure_vector)
}

estimate_models <- function(DGP_type, data_list, cores = 10) {
  
  # Register parallel backend
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  
  # Define the formula based on DGP type
  formula <- switch(DGP_type,
                    "DGP1" = Y ~ Z1 + Z2 + Z3 + Z4 + Z5 + A +
                      Z1:Z3 + Z1:Z4 + Z1:Z5 + 
                      Z3:Z4 + Z3:Z5 + Z4:Z5 + 
                      Z1:Z2:Z4 + Z1:Z2:Z5 + Z1:Z4:Z5 + Z2:Z4:Z5 + Z1:Z2:Z4:Z5,
                    "DGP4" = Y ~ Z1 + Z2 + Z3 + Z4 + Z5 + A +
                      Z1:Z3 + Z1:Z4 + Z1:Z5 +
                      Z3:Z4 + Z3:Z5 + Z4:Z5 +
                      Z1:Z4:Z2 + Z1:Z5:Z2 + Z1:Z4:Z5 + Z4:Z5:Z2 + Z1:Z4:Z5:Z2,
                    "DGP5" = Y ~ Z1 + Z2 + Z3 + Z4 + Z5 + Z6 + A +
                      Z1:Z3 + Z1:Z4 + Z1:Z5 +
                      Z3:Z4 + Z3:Z5 + Z4:Z5 +
                      Z1:Z6 + Z4:Z6 + Z5:Z6 +
                      Z1:Z4:Z6 + Z1:Z5:Z6 + Z1:Z4:Z5 + Z4:Z5:Z6 + Z1:Z4:Z5:Z6,
                    stop("Invalid DGP type")
  )
  
  # Parallel model estimation
  results <- foreach(i = 1:length(data_list), .packages = 'stats') %dopar% {
    data <- data_list[[i]]
    model <- lm(formula, data = data)
    coef_A <- coef(model)["A"]
    var_A <- vcov(model)["A", "A"]
    return(list(coef_A = coef_A, var_A = var_A))
  }
  
  # Stop the cluster
  stopCluster(cl)
  
  return(results)
}

# Example usage:
# Set global options for number formatting
options(scipen = 1, digits = 4)

models_sce1_DGP1 <- estimate_models("DGP1", DGP1_data$sce1_DGP1)
models_sce2_DGP1 <- estimate_models("DGP1", DGP1_data$sce2_DGP1)
models_sce3_DGP1 <- estimate_models("DGP1", DGP1_data$sce3_DGP1)
calculate_measures(models_sce1_DGP1, True_ATE = 0.20)
calculate_measures(models_sce2_DGP1, True_ATE = 0.20)
calculate_measures(models_sce3_DGP1, True_ATE = 0.24)


models_sce1_DGP2 <- estimate_models("DGP1", DGP2_data$sce1_DGP2)
models_sce2_DGP2 <- estimate_models("DGP1", DGP2_data$sce2_DGP2)
models_sce3_DGP2 <- estimate_models("DGP1", DGP2_data$sce3_DGP2)
calculate_measures(models_sce1_DGP2, True_ATE = 0.18)
calculate_measures(models_sce2_DGP2, True_ATE = 0.20)
calculate_measures(models_sce3_DGP2, True_ATE = 0.23)


models_sce1_DGP3 <- estimate_models("DGP1", DGP3_data$sce1_DGP3)
models_sce2_DGP3 <- estimate_models("DGP1", DGP3_data$sce2_DGP3)
models_sce3_DGP3 <- estimate_models("DGP1", DGP3_data$sce3_DGP3)
calculate_measures(models_sce1_DGP3, True_ATE = 0.18)
calculate_measures(models_sce2_DGP3, True_ATE = 0.19)
calculate_measures(models_sce3_DGP3, True_ATE = 0.22)


models_sce1_DGP4 <- estimate_models("DGP4", DGP4_data$sce1_DGP4)
models_sce2_DGP4 <- estimate_models("DGP4", DGP4_data$sce2_DGP4)
models_sce3_DGP4 <- estimate_models("DGP4", DGP4_data$sce3_DGP4)
calculate_measures(models_sce1_DGP4, True_ATE = 0.18)
calculate_measures(models_sce2_DGP4, True_ATE = 0.20)
calculate_measures(models_sce3_DGP4, True_ATE = 0.22)


models_sce1_DGP5 <- estimate_models("DGP5", DGP5_data$sce1_DGP5)
models_sce2_DGP5 <- estimate_models("DGP5", DGP5_data$sce2_DGP5)
models_sce3_DGP5 <- estimate_models("DGP5", DGP5_data$sce3_DGP5)
calculate_measures(models_sce1_DGP5, True_ATE = 0.18)
calculate_measures(models_sce2_DGP5, True_ATE = 0.20)
calculate_measures(models_sce3_DGP5, True_ATE = 0.22)









