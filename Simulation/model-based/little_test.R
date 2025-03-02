test <- all_missingnes_data[1:9]




test_method_results <- list()

for (orig_key in names(test)) {
  # Use regex to extract DGP and scenario from keys like "mDagT_DGP1_sce1"
  matches <- regmatches(orig_key, regexec("DGP(\\d+)_sce(\\d+)", orig_key))[[1]]
  if (length(matches) < 3) {
    stop(paste("Name", orig_key, "does not match the expected pattern 'DGP<dgp>_sce<scenario>'."))
  }
  
  # Construct the DGP string (e.g., "DGP1") and extract scenario number
  DGP <- paste0("DGP", matches[2])
  scenario <- matches[3]
  
  # Call apply_all_methods for the current element of 'test'
  method_res <- apply_all_methods(
    data_list = test[[orig_key]], 
    m = 5, 
    cores = 10, 
    DGP = DGP, 
    truncation = TRUE, 
    maxit = 5
  )
  
  new_key <- paste0(orig_key, "_res")
  # Store the result under the new key
  test_method_results[[new_key]] <- method_res
}


# Step 1. Make sure each element from test_method_results is assigned to the global environment.
for(nm in names(test_method_results)) {
  assign(nm, test_method_results[[nm]])
}

# Step 2. Create a list of names and use these to call all_measures_TMLE, preserving the object names.
all_obj_names <- names(test_method_results)
all_measures_list <- lapply(all_obj_names, function(nm) {
  expr <- substitute(all_measures_TMLE(x), list(x = as.name(nm)))
  eval(expr, envir = .GlobalEnv)
})
names(all_measures_list) <- all_obj_names

# Step 3. Combine all results in one big data frame.
big_table <- do.call(rbind, all_measures_list)

# (Optional) Post-processing: adjusting factor levels, etc.
big_table$Scenario <- factor(big_table$Scenario, levels = c("Scenario1", "Scenario2", "Scenario3"))
big_table$Method <- factor(big_table$Method, levels = c("CC", "Ext", "Ext MCMI", "MI PMM", "MI Int", "MI CART", "MI CART MI", "MI RF", "MI RF MI", "MI Amelia"))
big_table$DGP <- as.numeric(gsub("DGP", "", big_table$DGP))
rownames(big_table) <- sub(".*\\.", "", rownames(big_table))

# Inspect the resulting big_table
str(big_table)
head(big_table)
