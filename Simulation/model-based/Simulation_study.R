
source("~/GitHub/MissingDataTmle/Simulation/model-based/01_used_libraries.R")
source("~/GitHub/MissingDataTmle/Simulation/model-based/Data_and_Missingness.R")
source("~/GitHub/MissingDataTmle/Simulation/model-based/general_Functions.R")
source("~/GitHub/MissingDataTmle/Simulation/model-based/analysis_Functions.R")

#### DGP ####

# Set simulation parameters
Sim <- 10   # Number of repetitions
n <- 2000     # Sample size per dataset

dgp_list <- c("DGP1", "DGP2", "DGP3", "DGP4", "DGP5")
scenarios <- c("sce1", "sce2", "sce3")

# Create data
big_data_list <- list()

for (dgp in dgp_list) {
  for (sce in scenarios) {
    list_name <- paste0(sce, "_", dgp)
    big_data_list[[list_name]] <- simulateDatasets(scenario = sce, dgp = dgp, n_datasets = Sim, n = n)
  }
}

save(big_data_list, file = "../model-based/Results/data_list.RData")

#### Full data assessment #####

# Check data proportions
full_data_proportions <- list()

for (data_name in names(big_data_list)) {
  full_data_proportions[[data_name]] <- check_data(big_data_list[[data_name]])
}

full_data_proportions
save(full_data_proportions, file = "../model-based/Results/data_proportions.RData")

# Check positivity violation in data
full_positivity_results <- list()

for (data_name in names(big_data_list)) {
  full_positivity_results[[data_name]] <- check_positivity(big_data_list[[data_name]])
}

full_positivity_results
save(full_positivity_results, file = "../model-based/Results/positivity_results.RData")

# Check power of effect
effective_power_results <- check_effective_power_bigdata(big_data_list, cores = 10)
effective_power_results
save(effective_power_results, file = "../model-based/Results/effective_power_results.RData")


#### Induce Missingness ####

missingness_types <- c("T", "A", "E", "I", "J")

all_missingnes_data <- list()      # Will store the modified (missingness-induced) datasets
all_missing_proportions <- list()  # Will store the missing proportions

for (m_type in missingness_types) {
  cat("Processing missingness type:", m_type, "\n")
  
  # Apply missingness function for this type on the entire big_data_list
  result <- apply_missingness_bigdata(big_data_list, missingness_type = m_type, coef_list = coef_list)
  
  # Iterate over each element of the result to rename and store outputs
  for (orig_key in names(result)) {
    # Extract scenario and DGP from the original key (expected format "sceX_DGPY")
    matches <- regmatches(orig_key, regexec("sce(\\d+)_DGP(\\d+)", orig_key))[[1]]
    if (length(matches) < 3) {
      stop(paste("Name", orig_key, "does not match the expected pattern."))
    }
    scenario <- matches[2]  
    dgp <- matches[3]     
    
    # Create the new key according to the desired pattern:
    new_key <- paste0("mDag", m_type, "_DGP", dgp, "_sce", scenario)
    
    # Store the outputs under the new key
    all_missingnes_data[[new_key]] <- result[[orig_key]]$modified_data
    all_missing_proportions[[new_key]] <- result[[orig_key]]$missing_summary
  }
}

save(all_missingnes_data, file = "../model-based/Results/all_missingnes_data.RData")


# Create a list to store the missingness data frames for each DGP
missing_proportion_list <- list()

# Loop over DGP numbers 1 to 5
for(d in 1:5) {
  dgp_pattern <- paste0("DGP", d)
  
  # Extract the keys from all_missing_proportions that contain the DGP pattern.
  keys_dgp <- names(all_missing_proportions)[grepl(dgp_pattern, names(all_missing_proportions))]
  df_dgp <- do.call(rbind, all_missing_proportions[keys_dgp])
  
  # Assign the row names as the keys (i.e. the element names from the original list)
  rownames(df_dgp) <- keys_dgp
  
  # Store the data frame in the list with the name of the DGP (e.g., "DGP1")
  missing_proportion_list[[dgp_pattern]] <- as.data.frame(df_dgp)
}

# Check missingness for the DGPs
missing_proportion_list 
save(missing_proportion_list, file = "../model-based/Results/missing_proportion_list.RData")



#### Analysis ####


all_results <- list()

for (orig_key in names(big_data_list)) {
  
  matches <- regmatches(orig_key, regexec("DGP(\\d+)_sce(\\d+)", orig_key))[[1]]
  if (length(matches) < 3) {
    stop(paste("Name", orig_key, "does not match the expected pattern 'DGP<dgp>_sce<scenario>'."))
  }
  
  DGP <- paste0("DGP", matches[2])
  scenario <- matches[3]  # as a string
  
  # Call the apply_all_methods function for this element, passing the proper DGP
  method_res <- apply_all_methods(
    data_list = big_data_list[[orig_key]], 
    m = 5, 
    cores = 10, 
    DGP = DGP, 
    truncation = TRUE, 
    maxit = 5
  )
  
  new_key <- paste0(orig_key, "_res")
  
  # Store the result under the new key
  all_results[[new_key]] <- method_res
}


# Step 1. Make sure each element from all_results is assigned to the global environment.
for(nm in names(all_results)) {
  assign(nm, all_results[[nm]])
}

# Step 2. Create a list of names and use these to call all_measures_TMLE, preserving the object names.
all_obj_names <- names(all_results)
all_measures_list <- lapply(all_obj_names, function(nm) {
  expr <- substitute(all_measures_TMLE(x), list(x = as.name(nm)))
  eval(expr, envir = .GlobalEnv)
})
names(all_measures_list) <- all_obj_names

# Step 3. Combine all results in one big data frame.
big_table <- do.call(rbind, all_measures_list)

# Post-processing: adjusting factor levels, etc.
big_table$Scenario <- factor(big_table$Scenario, levels = c("Scenario1", "Scenario2", "Scenario3"))
big_table$Method <- factor(big_table$Method, levels = c("CC", "Ext", "Ext MCMI", "MI PMM", "MI Int", "MI CART", "MI CART MI", "MI RF", "MI RF MI", "MI Amelia"))
big_table$DGP <- as.numeric(gsub("DGP", "", big_table$DGP))
rownames(big_table) <- sub(".*\\.", "", rownames(big_table))

# Inspect the resulting big_table
str(big_table)
head(big_table)
