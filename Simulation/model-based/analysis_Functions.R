
merged_CC_TMLE <- function(data_list, cores = 4, truncation = TRUE) {
  
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  
  start <- Sys.time()
  
  SL_lib <- c("SL.glm", "SL.glm.interaction", "SL.earth", "SL.mean")
  
  required_packages <- c("tmle", "SuperLearner", "ranger", "earth", "arm", "glmnet")
  
  results <- foreach(i = 1:length(data_list), .packages = required_packages) %dopar% {
    data <- data_list[[i]]
    data <- na.omit(data)
    
    cov_name <- setdiff(names(data), c("A", "Y", "B"))
    
    gbound_value <- if (truncation) 0.01 else 0.0001
    
    result <- tmle(Y = data$Y, A = data$A, 
                   W = data[, cov_name],
                   family = "gaussian", 
                   Q.SL.library = SL_lib, 
                   g.SL.library = SL_lib, 
                   cvQinit = FALSE, 
                   gbound = gbound_value, 
                   V.Q = 5, 
                   V.g = 5)
    
    res <- list(
      result$estimates$ATE$psi,
      result$estimates$ATE$var.psi
    )
    
    res
  }
  
  stopCluster(cl)
  end <- Sys.time()
  duration <- end - start
  print(duration)
  gc()
  
  return(results)
}


merged_Ext_TMLE <- function(data_list, cores = 4, DGP = "DGP1", truncation = TRUE) {
  
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  
  # Start timing
  start <- Sys.time()
  
  # Super Learner library
  SL_lib <- c("SL.glm", "SL.glm.interaction", "SL.earth", "SL.mean")
  
  # Required packages for foreach
  required_packages <- c("tmle", "SuperLearner", "ranger", "earth", "arm", "glmnet")
  
  # Apply TMLE in parallel using foreach
  results <- foreach(i = 1:length(data_list), .packages = required_packages) %dopar% {
    data <- data_list[[i]]
    cov_name <- setdiff(names(data), c("A", "Y", "B"))
    
    if (DGP == "DGP5") {
      # Exclude records with missing exposure, confounder data, or Z6
      data <- data[-(which(is.na(data$A)
                           |is.na(data$Z2)
                           |is.na(data$Z3)
                           |is.na(data$Z4)
                           |is.na(data$Z6))), ]
    } else {
      # Exclude records with missing exposure or confounder data
      data <- data[-(which(is.na(data$A)
                           |is.na(data$Z2)
                           |is.na(data$Z3)
                           |is.na(data$Z4))), ]
    }
    
    # Create deltaY variable to estimate the probability of having observed outcome conditional on the exposure and confounders
    data$deltaY <- ifelse(is.na(data$Y), 0, 1)
    
    gbound_value <- if (truncation) 0.01 else 0.0001
    
    result <- tmle(Y = data$Y, A = data$A, 
                   W = data[, cov_name],
                   family = "gaussian", 
                   Q.SL.library = SL_lib, 
                   g.SL.library = SL_lib, 
                   g.Delta.SL.library = SL_lib, 
                   Delta = data$deltaY, 
                   cvQinit = FALSE, 
                   gbound = gbound_value, 
                   V.Q = 5, 
                   V.g = 5)
    
    res <- list(
      result$estimates$ATE$psi,
      result$estimates$ATE$var.psi
    )
    
    res
  }
  
  # Stop the cluster
  stopCluster(cl)
  
  # End timing
  end <- Sys.time()
  duration <- end - start
  print(duration)
  
  # Garbage collection
  gc()
  
  return(results)
}


merged_Ext_MCMI_TMLE <- function(data_list, cores = 4, DGP = "DGP1", truncation = TRUE) {
  
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  
  # Start timing
  start <- Sys.time()
  
  # Super Learner library
  SL_lib <- c("SL.glm", "SL.glm.interaction", "SL.earth", "SL.mean")
  
  # Required packages for foreach
  required_packages <- c("tmle", "SuperLearner", "ranger", "earth", "arm", "glmnet")
  
  # Apply TMLE in parallel using foreach
  results <- foreach(i = 1:length(data_list), .packages = required_packages) %dopar% {
    data <- data_list[[i]]
    
    # Convert Confounders with missing observations to factors
    data$Z2 <- as.factor(data$Z2)
    data$Z3 <- as.factor(data$Z3)
    
    if (DGP == "DGP1") {
      data$Z4 <- as.factor(data$Z4)
      levels(data$Z4) <- c(levels(data$Z4), "Missing")
      data$Z4[is.na(data$Z4)] <- "Missing"
    } else if (DGP == "DGP5") {
      data$MZ4 <- ifelse(is.na(data$Z4), 0, 1) # Missing Indicator for Z4
      data$MZ6 <- ifelse(is.na(data$Z6), 0, 1) # Missing Indicator for Z6
      data$Z4[is.na(data$Z4)] <- 0  # Impute missing with 0
      data$Z6[is.na(data$Z6)] <- 0  # Impute missing with 0
    } else {
      data$MZ4 <- ifelse(is.na(data$Z4), 0, 1) # Missing Indicator
      data$Z4[is.na(data$Z4)] <- 0  # Impute missing with 0
    }
    
    # Create new category for 'NAs'
    levels(data$Z2) <- c(levels(data$Z2), "Missing")
    levels(data$Z3) <- c(levels(data$Z3), "Missing")
    
    # Change 'NAs' to 'Missing'
    data$Z2[is.na(data$Z2)] <- "Missing"
    data$Z3[is.na(data$Z3)] <- "Missing"
    
    # Exclude only observations with unobserved exposures
    data <- data[!is.na(data$A), ]
    cov_name <- setdiff(names(data), c("A", "Y", "B"))
    
    # Add deltaY for tmle (zero means y-value is missing)
    data$deltaY <- ifelse(is.na(data$Y), 0, 1)
    
    gbound_value <- if (truncation) 0.01 else 0.0001
    
    result <- tmle(Y = data$Y, A = data$A, 
                   W = data[, cov_name],
                   family = "gaussian", 
                   Q.SL.library = SL_lib, 
                   g.SL.library = SL_lib, 
                   g.Delta.SL.library = SL_lib, 
                   Delta = data$deltaY, 
                   cvQinit = FALSE, 
                   gbound = gbound_value, 
                   V.Q = 5, 
                   V.g = 5)
    
    res <- list(
      result$estimates$ATE$psi,
      result$estimates$ATE$var.psi
    )
    
    res
  }
  
  # Stop the cluster
  stopCluster(cl)
  
  # End timing
  end <- Sys.time()
  duration <- end - start
  print(duration)
  
  # Garbage collection
  gc()
  
  return(results)
}


merged_MI_PMM_TMLE <- function(data_list, m = 10, cores = 4, DGP = "DGP1", truncation = TRUE, maxit = 10) {
  
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  
  # Start timing
  start <- Sys.time()
  
  # Super Learner library
  SL_lib <- c("SL.glm", "SL.glm.interaction", "SL.earth", "SL.mean")
  
  # Required packages for foreach
  required_packages <- c("tmle", "mice", "SuperLearner", "ranger", "earth", "arm", "glmnet")
  
  # Apply TMLE with multiple imputations in parallel using foreach
  results <- foreach(i = 1:length(data_list), .packages = required_packages) %dopar% {
    
    data <- data_list[[i]]
    cov_name <- setdiff(names(data), c("A", "Y", "B"))
    gbound_value <- if (truncation) 0.01 else 0.0001
    
    # Convert Confounders and Exposure to factor
    data$Z2 <- as.factor(data$Z2)
    data$A <- as.factor(data$A)
    
    if (DGP == "DGP1") {
      data$Z3 <- as.factor(data$Z3)
      data$Z4 <- as.factor(data$Z4)
      imp_method <- c("", "", "logreg", "logreg", "logreg", "", "logreg", "pmm")
    } else if (DGP == "DGP4") {
      imp_method <- c("", "", "logreg", "polyreg", "norm", "", "logreg", "pmm")
    } else if (DGP == "DGP5") {
      imp_method <- c("", "", "logreg", "polyreg", "norm", "", "pmm", "logreg", "pmm")
    } else {
      data$Z3 <- as.factor(data$Z3)
      imp_method <- c("", "", "logreg", "logreg", "norm", "", "logreg", "pmm")
    }
    
    # Specify imputation methods and conduct imputation
    PMM_imp <- mice(data, method = imp_method, m = m, maxit = maxit)
    
    # Create lists to store the different imputed datasets and respective analysis
    PMM_imp_data_list <- vector("list", length = m)
    PMM_imp_TMLE <- vector("list", length = m)
    ModVar <- vector(length = m)
    Estimate <- vector(length = m)
    
    for(l in 1:m) {
      # Estimate for each imputation the ATE and store it
      PMM_imp_data_list[[l]] <- mice::complete(PMM_imp, l)
      PMM_imp_data_list[[l]]$Z2 <- as.integer(ifelse(PMM_imp_data_list[[l]]$Z2 == "0", 0, 1))
      PMM_imp_data_list[[l]]$A <- as.integer(ifelse(PMM_imp_data_list[[l]]$A == "0", 0, 1))
      
      if (DGP == "DGP1") {
        PMM_imp_data_list[[l]]$Z3 <- as.integer(ifelse(PMM_imp_data_list[[l]]$Z3 == "0", 0, 1))
        PMM_imp_data_list[[l]]$Z4 <- as.integer(ifelse(PMM_imp_data_list[[l]]$Z4 == "0", 0, 1))
      } else if (DGP == "DGP4" || DGP == "DGP5") {
        # do nothing
      } else {
        PMM_imp_data_list[[l]]$Z3 <- as.integer(ifelse(PMM_imp_data_list[[l]]$Z3 == "0", 0, 1))
      }
      
      PMM_imp_TMLE[[l]] <- tmle(Y = PMM_imp_data_list[[l]]$Y, A = PMM_imp_data_list[[l]]$A, 
                                W = PMM_imp_data_list[[l]][, cov_name],
                                family = "gaussian", 
                                Q.SL.library = SL_lib, 
                                g.SL.library = SL_lib,
                                cvQinit = FALSE, 
                                gbound = gbound_value, 
                                V.Q = 5, 
                                V.g = 5)
      
      ModVar[l] <- PMM_imp_TMLE[[l]]$estimates$ATE$var.psi
      Estimate[l] <- PMM_imp_TMLE[[l]]$estimates$ATE$psi
    }
    
    imp_tmle_results <- vector("list", length = 2)
    imp_tmle_results[[1]] <- mean(Estimate)
    # Total variance = average of the complete-data variances + variance between the m complete-data estimates
    imp_tmle_results[[2]] <- mean(ModVar) + (1 + 1/m) * ((1/(m-1)) * (sum((Estimate - mean(Estimate))^2)))
    
    imp_tmle_results
    
  }
  
  # Stop the cluster
  stopCluster(cl)
  
  # End timing
  end <- Sys.time()
  duration <- end - start
  print(duration)
  
  # Garbage collection
  gc()
  
  return(results)
}


merged_MI_CART_TMLE <- function(data_list, m = 10, cores = 4, DGP = "DGP1", truncation = TRUE, maxit = 10) {
  
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  
  # Start timing
  start <- Sys.time()
  
  # Super Learner library
  SL_lib <- c("SL.glm", "SL.glm.interaction", "SL.earth", "SL.mean")
  
  # Required packages for foreach
  required_packages <- c("tmle", "mice", "SuperLearner", "ranger", "earth", "arm", "glmnet")
  
  # Apply TMLE with multiple imputations in parallel using foreach
  results <- foreach(i = 1:length(data_list), .packages = required_packages) %dopar% {
    
    data <- data_list[[i]]
    cov_name <- setdiff(names(data), c("A", "Y", "B"))
    gbound_value <- if (truncation) 0.01 else 0.0001
    
    # Convert Confounders and Exposure to factor
    data$Z2 <- as.factor(data$Z2)
    data$A <- as.factor(data$A)
    
    if (DGP == "DGP1") {
      data$Z3 <- as.factor(data$Z3)
      data$Z4 <- as.factor(data$Z4)
      imp_method <- c("", "", "cart", "cart", "cart", "", "cart", "cart")
    } else if (DGP == "DGP4") {
      imp_method <-  c("", "", "cart", "cart", "cart", "", "cart", "cart")
    } else if (DGP == "DGP5") {
      imp_method <-  c("", "", "cart", "cart", "cart", "", "cart", "cart", "cart")
    } else if (DGP == "DGP2" || DGP == "DGP3") {
      data$Z3 <- as.factor(data$Z3)
      imp_method <- c("", "", "cart", "cart", "cart", "", "cart", "cart")
    }
    
    # Specify imputation methods and conduct imputation
    Cart_imp <- mice(data, method = imp_method, m = m, maxit = maxit)
    
    # Create lists to store the different imputed datasets and respective analysis
    Cart_imp_data_list <- vector("list", length = m)
    Cart_imp_TMLE <- vector("list", length = m)
    ModVar <- vector(length = m)
    Estimate <- vector(length = m)
    
    for(l in 1:m) {
      # Estimate for each imputation the ATE and store it
      Cart_imp_data_list[[l]] <- mice::complete(Cart_imp, l)
      Cart_imp_data_list[[l]]$Z2 <- as.integer(ifelse(Cart_imp_data_list[[l]]$Z2 == "0", 0, 1))
      Cart_imp_data_list[[l]]$A <- as.integer(ifelse(Cart_imp_data_list[[l]]$A == "0", 0, 1))
      
      if (DGP == "DGP1") {
        Cart_imp_data_list[[l]]$Z3 <- as.integer(ifelse(Cart_imp_data_list[[l]]$Z3 == "0", 0, 1))
        Cart_imp_data_list[[l]]$Z4 <- as.integer(ifelse(Cart_imp_data_list[[l]]$Z4 == "0", 0, 1))
      } else if (DGP == "DGP4" || DGP == "DGP5") {
        # do nothing
      } else if (DGP == "DGP2" || DGP == "DGP3") {
        Cart_imp_data_list[[l]]$Z3 <- as.integer(ifelse(Cart_imp_data_list[[l]]$Z3 == "0", 0, 1))
      }
      
      Cart_imp_TMLE[[l]] <- tmle(Y = Cart_imp_data_list[[l]]$Y, A = Cart_imp_data_list[[l]]$A, 
                                 W = Cart_imp_data_list[[l]][, cov_name],
                                 family = "gaussian", 
                                 Q.SL.library = SL_lib, 
                                 g.SL.library = SL_lib,
                                 cvQinit = FALSE, gbound = gbound_value, V.Q = 5, V.g = 5)
      
      ModVar[l] <- Cart_imp_TMLE[[l]]$estimates$ATE$var.psi
      Estimate[l] <- Cart_imp_TMLE[[l]]$estimates$ATE$psi
    }
    
    imp_tmle_results <- vector("list", length = 2)
    imp_tmle_results[[1]] <- mean(Estimate)
    # Total variance = average of the complete-data variances + variance between the m complete-data estimates
    imp_tmle_results[[2]] <- mean(ModVar) + (1 + 1/m) * ((1/(m-1)) * (sum((Estimate - mean(Estimate))^2)))
    
    imp_tmle_results
    
  }
  
  # Stop the cluster
  stopCluster(cl)
  
  # End timing
  end <- Sys.time()
  duration <- end - start
  print(duration)
  
  # Garbage collection
  gc()
  
  return(results)
}


merged_MI_RF_TMLE <- function(data_list, m = 10, cores = 4, DGP = "DGP1", truncation = TRUE, maxit = 10) {
  
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  
  # Start timing
  start <- Sys.time()
  
  # Super Learner library
  SL_lib <- c("SL.glm", "SL.glm.interaction", "SL.earth", "SL.mean")
  
  # Required packages for foreach
  required_packages <- c("tmle", "mice", "SuperLearner", "ranger", "earth", "arm", "glmnet")
  
  # Apply TMLE with multiple imputations in parallel using foreach
  results <- foreach(i = 1:length(data_list), .packages = required_packages) %dopar% {
    data <- data_list[[i]]
    cov_name <- setdiff(names(data), c("A", "Y", "B"))
    gbound_value <- if (truncation) 0.01 else 0.0001
    
    # Convert Confounders and Exposure to factor
    data$Z2 <- as.factor(data$Z2)
    data$A <- as.factor(data$A)
    
    if (DGP == "DGP1") {
      data$Z3 <- as.factor(data$Z3)
      data$Z4 <- as.factor(data$Z4)
      meth <- c("", "", "rf", "rf", "rf", "", "rf", "rf")
    } else if (DGP == "DGP4") {
      meth <- c("", "", "rf", "rf", "rf", "", "rf", "rf")
      # Include MZ4 and MZ6 if necessary
    } else if (DGP == "DGP5") {
      meth <- c("", "", "rf", "rf", "rf", "", "rf", "rf", "rf")
      # Include MZ4 and MZ6 if necessary
    } else if (DGP == "DGP2" || DGP == "DGP3") {
      data$Z3 <- as.factor(data$Z3)
      meth <- c("", "", "rf", "rf", "rf", "", "rf", "rf")
    }
    
    # Specify imputation methods and conduct imputation
    RF_imp <- mice(data, meth = meth, m = m, maxit = maxit)
    
    # Create lists to store the different imputed datasets and respective analysis
    RF_imp_data_list <- vector("list", length = m)
    RF_imp_TMLE <- vector("list", length = m)
    ModVar <- vector(length = m)
    Estimate <- vector(length = m)
    
    for(l in 1:m) {
      # Estimate for each imputation the ATE and store it
      RF_imp_data_list[[l]] <- mice::complete(RF_imp, l)
      RF_imp_data_list[[l]]$Z2 <- as.integer(ifelse(RF_imp_data_list[[l]]$Z2 == "0", 0, 1))
      RF_imp_data_list[[l]]$A <- as.integer(ifelse(RF_imp_data_list[[l]]$A == "0", 0, 1))
      
      if (DGP == "DGP1") {
        RF_imp_data_list[[l]]$Z3 <- as.integer(ifelse(RF_imp_data_list[[l]]$Z3 == "0", 0, 1))
        RF_imp_data_list[[l]]$Z4 <- as.integer(ifelse(RF_imp_data_list[[l]]$Z4 == "0", 0, 1))
      } else if (DGP == "DGP4" || DGP == "DGP5") {
        # do nothing
      } else if (DGP == "DGP2" || DGP == "DGP3") {
        RF_imp_data_list[[l]]$Z3 <- as.integer(ifelse(RF_imp_data_list[[l]]$Z3 == "0", 0, 1))
      }
      
      RF_imp_TMLE[[l]] <- tmle(Y = RF_imp_data_list[[l]]$Y, A = RF_imp_data_list[[l]]$A, 
                               W = RF_imp_data_list[[l]][, cov_name],
                               family = "gaussian", 
                               Q.SL.library = SL_lib, 
                               g.SL.library = SL_lib,
                               cvQinit = FALSE, gbound = gbound_value, V.Q = 5, V.g = 5)
      
      ModVar[l] <- RF_imp_TMLE[[l]]$estimates$ATE$var.psi
      Estimate[l] <- RF_imp_TMLE[[l]]$estimates$ATE$psi
    }
    
    imp_tmle_results <- vector("list", length = 2)
    imp_tmle_results[[1]] <- mean(Estimate)
    # Total variance = average of the complete-data variances + variance between the m complete-data estimates
    imp_tmle_results[[2]] <- mean(ModVar) + (1 + 1/m) * ((1/(m-1)) * (sum((Estimate - mean(Estimate))^2)))
    
    imp_tmle_results
  }
  
  # Stop the cluster
  stopCluster(cl)
  
  # End timing
  end <- Sys.time()
  duration <- end - start
  print(duration)
  
  # Garbage collection
  gc()
  
  return(results)
}


merged_MI_Amelia_TMLE <- function(data_list, m = 10, cores = 4, DGP = "DGP1", truncation = TRUE) {
  
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  
  # Start timing
  start <- Sys.time()
  
  # Super Learner library
  SL_lib <- c("SL.glm", "SL.glm.interaction", "SL.earth", "SL.mean")
  
  # Required packages for foreach
  required_packages <- c("tmle", "Amelia", "SuperLearner", "ranger", "earth", "arm", "glmnet")
  
  # Apply TMLE with multiple imputations in parallel using foreach
  results <- foreach(i = 1:length(data_list), .packages = required_packages) %dopar% {
    data <- data_list[[i]]
    cov_name <- setdiff(names(data), c("A", "Y", "B"))
    gbound_value <- if (truncation) 0.01 else 0.0001
    
    # Convert Confounders and Exposure to factor
    data$Z1 <- as.factor(data$Z1)
    data$Z2 <- as.factor(data$Z2)
    data$A <- as.factor(data$A)
    
    if (DGP == "DGP1") {
      data$Z3 <- as.factor(data$Z3)
      data$Z4 <- as.factor(data$Z4)
      data$Z5 <- as.factor(data$Z5)
      noms <- c("Z1", "Z2", "Z3", "Z4", "Z5", "A")
    } else if (DGP == "DGP4") {
      noms <- c("Z1", "Z2", "Z3", "A")
    } else if (DGP == "DGP5") {
      data$Z6 <- log(data$Z6)
      noms <- c("Z1", "Z2", "Z3", "A")
    } else if (DGP == "DGP2" || DGP == "DGP3") {
      data$Z3 <- as.factor(data$Z3)
      noms <- c("Z1", "Z2", "Z3", "A")
    }
    
    # Specify factors and conduct imputation
    Amelia_imp <- amelia(data, m = m, noms = noms)
    
    # Create lists to store the different imputed datasets and respective analysis
    Amelia_imp_data_list <- vector("list", length = m)
    Amelia_imp_TMLE <- vector("list", length = m)
    ModVar <- vector(length = m)
    Estimate <- vector(length = m)
    
    for(l in 1:m) {
      # Estimate for each imputation the ATE and store it
      Amelia_imp_data_list[[l]] <- Amelia_imp$imputations[[l]]
      Amelia_imp_data_list[[l]]$Z1 <- as.integer(ifelse(Amelia_imp_data_list[[l]]$Z1 == "0", 0, 1))
      Amelia_imp_data_list[[l]]$Z2 <- as.integer(ifelse(Amelia_imp_data_list[[l]]$Z2 == "0", 0, 1))
      Amelia_imp_data_list[[l]]$A <- as.integer(ifelse(Amelia_imp_data_list[[l]]$A == "0", 0, 1))
      
      if (DGP == "DGP1") {
        Amelia_imp_data_list[[l]]$Z3 <- as.integer(ifelse(Amelia_imp_data_list[[l]]$Z3 == "0", 0, 1))
        Amelia_imp_data_list[[l]]$Z4 <- as.integer(ifelse(Amelia_imp_data_list[[l]]$Z4 == "0", 0, 1))
        Amelia_imp_data_list[[l]]$Z5 <- as.integer(ifelse(Amelia_imp_data_list[[l]]$Z5 == "0", 0, 1))
      } else if (DGP == "DGP2" || DGP == "DGP3") {
        Amelia_imp_data_list[[l]]$Z3 <- as.integer(ifelse(Amelia_imp_data_list[[l]]$Z3 == "0", 0, 1))
      } else if (DGP == "DGP5") {
        Amelia_imp_data_list[[l]]$Z6 <- exp(Amelia_imp_data_list[[l]]$Z6)
      }
      
      Amelia_imp_TMLE[[l]] <- tmle(Y = Amelia_imp_data_list[[l]]$Y, A = Amelia_imp_data_list[[l]]$A, 
                                   W = Amelia_imp_data_list[[l]][, cov_name],
                                   family = "gaussian", 
                                   Q.SL.library = SL_lib, 
                                   g.SL.library = SL_lib,
                                   cvQinit = FALSE, gbound = gbound_value, V.Q = 5, V.g = 5)
      
      ModVar[l] <- Amelia_imp_TMLE[[l]]$estimates$ATE$var.psi
      Estimate[l] <- Amelia_imp_TMLE[[l]]$estimates$ATE$psi
    }
    
    imp_tmle_results <- vector("list", length = 2)
    imp_tmle_results[[1]] <- mean(Estimate)
    # Total variance = average of the complete-data variances + variance between the m complete-data estimates
    imp_tmle_results[[2]] <- mean(ModVar) + (1 + 1/m) * ((1/(m-1)) * (sum((Estimate - mean(Estimate))^2)))
    
    imp_tmle_results
  }
  
  # Stop the cluster
  stopCluster(cl)
  
  # End timing
  end <- Sys.time()
  duration <- end - start
  print(duration)
  
  # Garbage collection
  gc()
  
  return(results)
}

run_MI_Int <- function(data_list, m = 10, cores = 4, DGP = "DGP1", truncation = TRUE, verbose = FALSE, maxit = 10) {
  
  # Initialize parallel backend
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  
  if (verbose) cat("Parallel backend initialized with", cores, "cores.\n")
  
  # Start timing
  start_time <- Sys.time()
  
  # Define required packages for foreach
  required_packages <- c("mice", "tmle", "SuperLearner", "ranger", "earth", "arm", "glmnet")
  
  # Super Learner library
  SL_lib <- c("SL.glm", "SL.glm.interaction", "SL.earth", "SL.mean")
  
  # Parallel processing using foreach
  tmle_results <- foreach(i = 1:length(data_list), .packages = required_packages) %dopar% {
    # Initialize result for this iteration
    result <- list(success = FALSE, estimate = NULL, variance = NULL, error_msg = NULL)
    
    # Try-catch block to handle errors
    tryCatch({
      data <- data_list[[i]]
      cov_name <- names(data)
      gbound_value <- if (truncation) 0.01 else 0.0001
      
      if (DGP == "DGP1") {
        # Create interaction terms and additional columns for DGP1
        data <- cbind(Z2.fac = as.factor(data$Z2), Z3.fac = as.factor(data$Z3), Z4.fac = as.factor(data$Z4),
                      A.fac = as.factor(data$A), data, 
                      A.Y = (data$A * data$Y), A.Z1 = as.integer(data$A * data$Z1), 
                      Y.Z1 = as.integer(data$Y * data$Z1), A.Z2 = as.integer(data$A * data$Z2),
                      Y.Z2 = as.integer(data$Y * data$Z2), A.Z3 = as.integer(data$A * data$Z3), 
                      Y.Z3 = (data$Y * data$Z3), A.Z4 = as.integer(data$A * data$Z4), 
                      Y.Z4 = (data$Y * data$Z4), A.Z5 = as.integer(data$A * data$Z5),
                      Y.Z5 = (data$Y * data$Z5), Z1.Z3 = as.integer(data$Z1 * data$Z3), 
                      Z1.Z4 = as.integer(data$Z1 * data$Z4), Z1.Z5 = as.integer(data$Z1 * data$Z5), 
                      Z3.Z4 = as.integer(data$Z3 * data$Z4), Z3.Z5 = as.integer(data$Z3 * data$Z5),
                      Z4.Z5 = as.integer(data$Z4 * data$Z5), Z2.Z1 = as.integer(data$Z2 * data$Z1),
                      Z2.Z3 = as.integer(data$Z2 * data$Z3), Z2.Z4 = as.integer(data$Z2 * data$Z4),
                      Z2.Z5 = as.integer(data$Z2 * data$Z5), 
                      Z1.Z2.Z4 = as.integer(data$Z1 * data$Z2 * data$Z4), 
                      Z1.Z2.Z5 = as.integer(data$Z1 * data$Z2 * data$Z5),
                      Z1.Z4.Z5 = as.integer(data$Z1 * data$Z4 * data$Z5), 
                      Z2.Z4.Z5 = as.integer(data$Z2 * data$Z4 * data$Z5), 
                      Z1.Z2.Z4.Z5 = as.integer(data$Z1 * data$Z2 * data$Z4 * data$Z5))
        
        # Initializing mice to create a template
        template_mice <- mice(data, max = 0, print = FALSE)
        
        methods <- template_mice$method
        methods[c("Z2.fac", "Z3.fac", "Z4.fac", "A.fac", "Y")] <- c("logreg", "logreg", "logreg", "logreg", "pmm")
        methods[c("Z2", "Z3", "Z4", "A")] <- c("~I((as.integer(Z2.fac) -1))", "~I((as.integer(Z3.fac) -1))", "~I((as.integer(Z4.fac) -1))", "~I((as.integer(A.fac) -1))")
        methods[colnames(data)[13:38]] <- c("~I(A*Y)", "~I(A*Z1)", "~I(Y*Z1)", "~I(A*Z2)", "~I(Y*Z2)", "~I(A*Z3)", "~I(Y*Z3)", "~I(A*Z4)", "~I(Y*Z4)",
                                            "~I(A*Z5)", "~I(Y*Z5)", "~I(Z1*Z3)", "~I(Z1*Z4)", "~I(Z1*Z5)", "~I(Z3*Z4)", "~I(Z3*Z5)",
                                            "~I(Z4*Z5)", "~I(Z2*Z1)", "~I(Z2*Z3)", "~I(Z2*Z4)", "~I(Z2*Z5)", "~I(Z1*Z2*Z4)", "~I(Z1*Z2*Z5)", "~I(Z1*Z4*Z5)", "~I(Z2*Z4*Z5)",
                                            "~I(Z1*Z2*Z4*Z5)")
        
        predictions <- template_mice$predictorMatrix
        predictions[13:38, ] <- 0
        predictions["Z2.fac", c("Z2", "Z3", "Z4", "A", "A.Z2", "Y.Z2", "Z2.Z1", "Z2.Z3", "Z2.Z4", "Z2.Z5", "Z1.Z2.Z4", "Z1.Z2.Z5", "Z2.Z4.Z5", "Z1.Z2.Z4.Z5")] <- 0
        predictions["Z3.fac", c("Z2", "Z3", "Z4", "A", "A.Z3", "Y.Z3", "Z1.Z3", "Z3.Z4", "Z3.Z5", "Z2.Z3")] <- 0
        predictions["Z4.fac", c("Z2", "Z3", "Z4", "A", "A.Z4", "Y.Z4", "Z1.Z4", "Z3.Z4", "Z4.Z5", "Z2.Z4", "Z1.Z2.Z4", "Z1.Z4.Z5", "Z2.Z4.Z5", "Z1.Z2.Z4.Z5")] <- 0
        predictions["A.fac", c("Z2", "Z3", "Z4", "A", "A.Y", "A.Z1", "A.Z2", "A.Z3", "A.Z4", "A.Z5")] <- 0
        predictions["Y", c("Z2", "Z3", "Z4", "A", "A.Y", "Y.Z1", "Y.Z2", "Y.Z3", "Y.Z4", "Y.Z5")] <- 0
      } else if (DGP == "DGP4") {
        # Additional manipulation for DGP4
        data$Z3_no_NA <- ifelse(is.na(data$Z3), "Missing", as.character(data$Z3))
        data$Z3_no_NA <- as.factor(data$Z3_no_NA)
        Z3_dummies <- model.matrix(~ Z3_no_NA - 1, data = data)
        missing_vector <- as.numeric(data$Z3_no_NA == "Missing")
        Z3_dummies <- cbind(Z3_dummies, Missing = missing_vector)
        Z3_dummies <- apply(Z3_dummies, 2, function(x) ifelse(x == 0 & missing_vector == 1, NA, x))
        Z3_dummies <- Z3_dummies[, c(1:3)]
        colnames(Z3_dummies) <- c("Z3_1", "Z3_2", "Z3_3")
        data <- cbind(data[,-9], Z3_dummies)
        
        data <- cbind(Z2.fac = as.factor(data$Z2), A.fac = as.factor(data$A), data, 
                      A.Y = data$A * data$Y, A.Z1 = data$A * data$Z1, 
                      Y.Z1 = data$Y * data$Z1, A.Z2 = data$A * data$Z2,
                      Y.Z2 = data$Y * data$Z2, A.Z31 = data$A * data$Z3_1, 
                      Y.Z31 = data$Y * data$Z3_1, A.Z32 = data$A * data$Z3_2,
                      Y.Z32 = data$Y * data$Z3_2, A.Z33 = data$A * data$Z3_3,
                      Y.Z33 = data$Y * data$Z3_3, A.Z4 = data$A * data$Z4, 
                      Y.Z4 = data$Y * data$Z4, A.Z5 = data$A * data$Z5,
                      Y.Z5 = data$Y * data$Z5, Z31.Z1 = data$Z1 * data$Z3_1, 
                      Z32.Z1 = data$Z1 * data$Z3_2, Z33.Z1 = data$Z1 * data$Z3_3,
                      Z1.Z4 = data$Z1 * data$Z4, Z1.Z5 = data$Z1 * data$Z5, 
                      Z31.Z4 = data$Z3_1 * data$Z4, Z32.Z4 = data$Z3_2 * data$Z4,
                      Z33.Z4 = data$Z3_3 * data$Z4, Z31.Z5 = data$Z3_1 * data$Z5,
                      Z32.Z5 = data$Z3_2 * data$Z5, Z33.Z5 = data$Z3_3 * data$Z5,
                      Z4.Z5 = data$Z4 * data$Z5, Z2.Z1 = data$Z2 * data$Z1,
                      Z2.Z4 = data$Z2 * data$Z4,
                      Z2.Z5 = data$Z2 * data$Z5, 
                      Z1.Z4.Z2 = data$Z1 * data$Z4 * data$Z2, Z1.Z5.Z2 = data$Z1 * data$Z5 * data$Z2,
                      Z1.Z4.Z5 = data$Z1 * data$Z4 * data$Z5, Z4.Z5.Z2 = data$Z4 * data$Z5 * data$Z2,
                      Z1.Z4.Z5.Z2 = data$Z1 * data$Z4 * data$Z5 * data$Z2)
        
        template_mice <- vector("list", length = 1)
        template_mice <- mice(data, max = 0, print = FALSE)
        
        methods <- template_mice$method
        methods[c("Z2.fac", "A.fac", "Z3", "Z4", "Y")] <- c("logreg", "logreg", "polyreg", "norm", "pmm") 
        methods[c("Z2", "Z3_1", "Z3_2", "Z3_3", "A")] <- c("~I((as.integer(Z2.fac) -1))", "~I(as.integer(model.matrix(~ Z3 - 1)[,1]))",
                                                           "~I(as.integer(model.matrix(~ Z3 - 1)[,2]))", "~I(as.integer(model.matrix(~ Z3 - 1)[,3]))",
                                                           "~I((as.integer(A.fac) -1))")
        methods[colnames(data)[14:48]] <- c("~I(A*Y)", "~I(A*Z1)", "~I(Y*Z1)", "~I(A*Z2)", "~I(Y*Z2)", "~I(A*Z3_1)", "~I(Y*Z3_1)", "~I(A*Z3_2)", "~I(Y*Z3_2)", "~I(A*Z3_3)", "~I(Y*Z3_3)", 
                                            "~I(A*Z4)", "~I(Y*Z4)", "~I(A*Z5)", "~I(Y*Z5)",
                                            "~I(Z1*Z3_1)", "~I(Z1*Z3_2)", "~I(Z1*Z3_3)", "~I(Z1*Z4)", "~I(Z1*Z5)", "~I(Z3_1*Z4)", "~I(Z3_2*Z4)", "~I(Z3_3*Z4)", 
                                            "~I(Z3_1*Z5)", "~I(Z3_2*Z5)", "~I(Z3_3*Z5)", "~I(Z4*Z5)", 
                                            "~I(Z2*Z1)", "~I(Z2*Z4)", "~I(Z2*Z5)",
                                            "~I(Z1*Z4*Z2)", "~I(Z1*Z5*Z2)", "~I(Z1*Z4*Z5)", "~I(Z4*Z5*Z2)", "~I(Z1*Z4*Z5*Z2)")
        predictions <- template_mice$predictorMatrix
        predictions[14:48, ] <- 0
        predictions["Z2.fac",  c("Z2", "Z3_1", "Z3_2", "Z3_3", "A", "A.Z2", "Y.Z2", "Z2.Z1", "Z2.Z4", "Z2.Z5",
                                 "Z1.Z4.Z2", "Z1.Z5.Z2", "Z4.Z5.Z2", "Z1.Z4.Z5.Z2")] <- 0
        predictions["A.fac", c("Z2", "Z3_1", "Z3_2", "Z3_3", "A",
                               "A.Y", "A.Z1", "A.Z2", "A.Z31", "A.Z32", "A.Z33", "A.Z4", "A.Z5")] <- 0
        predictions["Z3", c("Z2", "Z3_1", "Z3_2", "Z3_3", "A", "A.Z31", "Y.Z31", "A.Z32", "Y.Z32", "A.Z33", "Y.Z33",
                            "Z31.Z4", "Z32.Z4", "Z33.Z4", "Z31.Z5", "Z32.Z5", "Z33.Z5")] <- 0
        predictions["Z4", c("Z2", "Z3_1", "Z3_2", "Z3_3", "A", "Z2.Z4",
                            "A.Z4", "Y.Z4", "Z1.Z4", "Z31.Z4", "Z32.Z4", "Z33.Z4", "Z4.Z5", 
                            "Z1.Z4.Z2", "Z1.Z4.Z5", "Z4.Z5.Z2", "Z1.Z4.Z5.Z2")] <- 0 
        predictions["Y", c("Z2", "Z3_1", "Z3_2", "Z3_3", "A",
                           "A.Y", "Y.Z1", "Y.Z2", "Y.Z31", "Y.Z32", "Y.Z33", "Y.Z4", "Y.Z5")] <- 0
      } else if (DGP == "DGP5") { 
        data$Z3_no_NA <- ifelse(is.na(data$Z3), "Missing", as.character(data$Z3))
        data$Z3_no_NA <- as.factor(data$Z3_no_NA)
        Z3_dummies <- model.matrix(~ Z3_no_NA - 1, data = data)
        missing_vector <- as.numeric(data$Z3_no_NA == "Missing")
        Z3_dummies <- cbind(Z3_dummies, Missing = missing_vector)
        Z3_dummies <- apply(Z3_dummies, 2, function(x) ifelse(x == 0 & missing_vector == 1, NA, x))
        Z3_dummies <- Z3_dummies[, c(1:3)]
        colnames(Z3_dummies) <- c("Z3_1", "Z3_2", "Z3_3")
        data <- cbind(data[,-10], Z3_dummies)
        
        data <- cbind(Z2.fac = as.factor(data$Z2), A.fac = as.factor(data$A), data, 
                      A.Y = data$A * data$Y, A.Z1 = data$A * data$Z1, 
                      Y.Z1 = data$Y * data$Z1, A.Z31 = data$A * data$Z3_1, 
                      Y.Z31 = data$Y * data$Z3_1, A.Z32 = data$A * data$Z3_2,
                      Y.Z32 = data$Y * data$Z3_2, A.Z33 = data$A * data$Z3_3,
                      Y.Z33 = data$Y * data$Z3_3, A.Z4 = data$A * data$Z4, 
                      Y.Z4 = data$Y * data$Z4, A.Z5 = data$A * data$Z5,
                      Y.Z5 = data$Y * data$Z5, A.Z6 = data$A * data$Z6,
                      Y.Z6 = data$Y * data$Z6, Z31.Z1 = data$Z1 * data$Z3_1, 
                      Z32.Z1 = data$Z1 * data$Z3_2, Z33.Z1 = data$Z1 * data$Z3_3,
                      Z1.Z4 = data$Z1 * data$Z4, Z1.Z5 = data$Z1 * data$Z5, 
                      Z31.Z4 = data$Z3_1 * data$Z4, Z32.Z4 = data$Z3_2 * data$Z4,
                      Z33.Z4 = data$Z3_3 * data$Z4, Z31.Z5 = data$Z3_1 * data$Z5,
                      Z32.Z5 = data$Z3_2 * data$Z5, Z33.Z5 = data$Z3_3 * data$Z5,
                      Z4.Z5 = data$Z4 * data$Z5, Z6.Z1 = data$Z6 * data$Z1, 
                      Z6.Z4 = data$Z6 * data$Z4, Z6.Z5 = data$Z6 * data$Z5,
                      Z1.Z4.Z6 = data$Z1 * data$Z4 * data$Z6, Z1.Z5.Z6 = data$Z1 * data$Z5 * data$Z6,
                      Z1.Z4.Z5 = data$Z1 * data$Z4 * data$Z5, Z4.Z5.Z6 = data$Z4 * data$Z5 * data$Z6,
                      Z1.Z4.Z5.Z6 = data$Z1 * data$Z4 * data$Z5 * data$Z6)
        
        template_mice <- mice(data, max = 0, print = FALSE)
        
        methods <- template_mice$method
        methods[c("Z2.fac", "A.fac", "Z3", "Z4", "Z6", "Y")] <- c("logreg", "logreg", "polyreg", "norm", "pmm", "pmm") 
        methods[c("Z2", "Z3_1", "Z3_2", "Z3_3", "A")] <- c("~I((as.integer(Z2.fac) -1))", "~I(as.integer(model.matrix(~ Z3 - 1)[,1]))",
                                                           "~I(as.integer(model.matrix(~ Z3 - 1)[,2]))", "~I(as.integer(model.matrix(~ Z3 - 1)[,3]))",
                                                           "~I((as.integer(A.fac) -1))")
        methods[colnames(data)[14:49]] <- c("~I(A*Y)", "~I(A*Z1)", "~I(Y*Z1)", "~I(A*Z3_1)", "~I(Y*Z3_1)", "~I(A*Z3_2)", "~I(Y*Z3_2)", "~I(A*Z3_3)", "~I(Y*Z3_3)", 
                                            "~I(A*Z4)", "~I(Y*Z4)", "~I(A*Z5)", "~I(Y*Z5)", "~I(A*Z6)", "~I(Y*Z6)", 
                                            "~I(Z1*Z3_1)", "~I(Z1*Z3_2)", "~I(Z1*Z3_3)", "~I(Z1*Z4)", "~I(Z1*Z5)", "~I(Z3_1*Z4)", "~I(Z3_2*Z4)", "~I(Z3_3*Z4)", 
                                            "~I(Z3_1*Z5)", "~I(Z3_2*Z5)", "~I(Z3_3*Z5)", "~I(Z4*Z5)", "~I(Z6*Z1)", "~I(Z6*Z2)", "~I(Z6*Z4)", "~I(Z6*Z5)",
                                            "~I(Z1*Z4*Z6)", "~I(Z1*Z5*Z6)", "~I(Z1*Z4*Z5)", "~I(Z4*Z5*Z6)", "~I(Z1*Z4*Z5*Z6)")
        predictions <- template_mice$predictorMatrix
        predictions[14:49,] <- 0
        predictions["Z2.fac", c("Z2", "Z3_1", "Z3_2", "Z3_3", "A")] <- 0
        predictions["A.fac", c("Z2", "Z3_1", "Z3_2", "Z3_3", "A", 
                               "A.Y", "A.Z1", "A.Z31", "A.Z32", "A.Z33", "A.Z4", "A.Z5", "A.Z6")] <- 0
        predictions["Z3", c("Z2", "Z3_1", "Z3_2", "Z3_3", "A", "A.Z31", "Y.Z31", "A.Z32", "Y.Z32", "A.Z33", "Y.Z33",
                            "Z31.Z4", "Z32.Z4", "Z33.Z4", "Z31.Z5", "Z32.Z5", "Z33.Z5")] <- 0
        predictions["Z4", c("Z2", "Z3_1", "Z3_2", "Z3_3", "A", "A.Z4", "Y.Z4", "Z1.Z4", "Z31.Z4", "Z32.Z4", "Z33.Z4", "Z4.Z5", "Z6.Z4",
                            "Z1.Z4.Z6", "Z1.Z4.Z5", "Z4.Z5.Z6", "Z1.Z4.Z5.Z6")] <- 0
        predictions["Z6", c("Z2", "Z3_1", "Z3_2", "Z3_3", "A", "A.Z6", "Y.Z6", "Z6.Z1", "Z6.Z4", "Z6.Z5",
                            "Z1.Z4.Z6", "Z1.Z5.Z6", "Z4.Z5.Z6", "Z1.Z4.Z5.Z6")] <- 0
        predictions["Y", c("Z2", "Z3_1", "Z3_2", "Z3_3", "A", 
                           "A.Y", "Y.Z1", "Y.Z31", "Y.Z32", "Y.Z33", "Y.Z4", "Y.Z5", "Y.Z6")] <- 0
        
      } else {
        # Create interaction terms and additional columns for other DGPs
        data <- cbind(Z2.fac = as.factor(data$Z2), Z3.fac = as.factor(data$Z3), 
                      A.fac = as.factor(data$A), data, 
                      A.Y = data$A*data$Y, A.Z1 = data$A*data$Z1, 
                      Y.Z1 = data$Y*data$Z1, A.Z2 = data$A*data$Z2,
                      Y.Z2 = data$Y*data$Z2, A.Z3 = data$A*data$Z3,
                      Y.Z3 = data$Y*data$Z3, A.Z4 = data$A*data$Z4, 
                      Y.Z4 = data$Y*data$Z4, A.Z5 = data$A*data$Z5,
                      Y.Z5 = data$Y*data$Z5, Z1.Z3 = data$Z1*data$Z3, 
                      Z1.Z4 = data$Z1*data$Z4, Z1.Z5 = data$Z1*data$Z5, 
                      Z3.Z4 = data$Z3*data$Z4, Z3.Z5 = data$Z3*data$Z5,
                      Z4.Z5 = data$Z4*data$Z5, Z2.Z1 = data$Z2*data$Z1,
                      Z2.Z3 = data$Z2*data$Z3, Z2.Z4 = data$Z2*data$Z4,
                      Z2.Z5 = data$Z2*data$Z5, 
                      Z1.Z4.Z2 = data$Z1*data$Z4*data$Z2, Z1.Z5.Z2 = data$Z1*data$Z5*data$Z2,
                      Z1.Z4.Z5 = data$Z1*data$Z4*data$Z5, Z4.Z5.Z2 = data$Z4*data$Z5*data$Z2,
                      Z1.Z4.Z5.Z2 = data$Z1*data$Z4*data$Z5*data$Z2)
        template_mice <- vector("list", length = 1)
        template_mice <- mice(data, max = 0, print = FALSE)
        
        methods <- template_mice$method
        methods[c("Z2.fac", "Z3.fac", "A.fac", "Z4", "Y")] <- c("logreg", "logreg", "logreg", "norm", "pmm") 
        methods[c("Z2", "Z3", "A")] <- c("~I((as.integer(Z2.fac) -1))", "~I((as.integer(Z3.fac) -1))", "~I((as.integer(A.fac) -1))")
        
        methods[colnames(data)[12:37]] <- c("~I(A*Y)", "~I(A*Z1)", "~I(Y*Z1)", "~I(A*Z2)", "~I(Y*Z2)", "~I(A*Z3)", "~I(Y*Z3)", 
                                            "~I(A*Z4)", "~I(Y*Z4)", "~I(A*Z5)", "~I(Y*Z5)",
                                            "~I(Z1*Z3)", "~I(Z1*Z4)", "~I(Z1*Z5)", "~I(Z3*Z4)", "~I(Z3*Z5)", "~I(Z4*Z5)", 
                                            "~I(Z2*Z1)", "~I(Z2*Z3)", "~I(Z2*Z4)", "~I(Z2*Z5)",
                                            "~I(Z1*Z4*Z2)", "~I(Z1*Z5*Z2)", "~I(Z1*Z4*Z5)", "~I(Z4*Z5*Z2)", "~I(Z1*Z4*Z5*Z2)")
        predictions <- template_mice$predictorMatrix
        predictions[12:37,] <- 0
        predictions["Z2.fac",  c("Z2", "Z3", "A", "A.Z2", "Y.Z2", "Z2.Z1", "Z2.Z3", "Z2.Z4", "Z2.Z5", "Z1.Z4.Z2", "Z1.Z5.Z2", "Z4.Z5.Z2", "Z1.Z4.Z5.Z2")] <- 0
        predictions["A.fac", c("Z2", "Z3", "A", "A.Y", "A.Z1", "A.Z2", "A.Z3", "A.Z4", "A.Z5")] <- 0
        predictions["Z3.fac", c("Z2", "Z3", "A", "A.Z3", "Y.Z3", "Z1.Z3", "Z3.Z4", "Z3.Z5", "Z2.Z3")] <- 0
        predictions["Z4", c("Z2", "Z3", "A", "A.Z4", "Y.Z4", "Z1.Z4", "Z3.Z4", "Z4.Z5", "Z2.Z4",
                            "Z1.Z4.Z2", "Z1.Z4.Z5", "Z4.Z5.Z2", "Z1.Z4.Z5.Z2")] <- 0 
        predictions["Y", c("Z2", "Z3", "A", "A.Y", "Y.Z1", "Y.Z2", "Y.Z3", "Y.Z4", "Y.Z5")] <- 0
      }
      
      # Run mice with specified methods and predictor matrix
      imp_data <- mice(data, meth = methods, pred = predictions, m = m, maxit = maxit, print = FALSE)
      
      # Create lists to store the different imputed datasets and respective analysis
      ImpDataList <- vector("list", length = m)
      TMLE_List <- vector("list", length = m)
      ModVar <- vector(length = m)
      Estimate <- vector(length = m)
      
      # Update covariate names (exclude certain variables)
      if (is.null(unique(imp_data$loggedEvents$out))) {
        cov_name <- cov_name[!(cov_name %in% c("A", "Y", "B"))]
      } else {
        cov_name <- cov_name[!(cov_name %in% c("A", "Y", "B", unique(imp_data$loggedEvents$out)))]
      }
      
      for (l in 1:m) {
        # Estimate for each imputation the ATE and store it
        ImpDataList[[l]] <- mice::complete(imp_data, l)
        TMLE_List[[l]] <- tmle(Y = ImpDataList[[l]]$Y, A = ImpDataList[[l]]$A, 
                               W = ImpDataList[[l]][, cov_name],
                               family = "gaussian", 
                               Q.SL.library = SL_lib, 
                               g.SL.library = SL_lib,
                               cvQinit = FALSE, gbound = gbound_value, V.Q = 5, V.g = 5)
        ModVar[l] <- TMLE_List[[l]]$estimates$ATE$var.psi
        Estimate[l] <- TMLE_List[[l]]$estimates$ATE$psi
      }
      
      # Compute combined estimate and variance using Rubin's rules
      combined_estimate <- mean(Estimate)
      combined_variance <- mean(ModVar) + (1 + 1/m) * ((1/(m-1)) * sum((Estimate - mean(Estimate))^2))
      
      # If successful, store the estimates
      result$success <- TRUE
      result$estimate <- combined_estimate
      result$variance <- combined_variance
      
    }, error = function(e) {
      # On error, record failure and store the error message
      result$success <- FALSE
      result$estimate <- NULL
      result$variance <- NULL
      result$error_msg <- conditionMessage(e)  # Capture the full error message
    })
    
    return(result)
  }
  
  # Stop the cluster
  stopCluster(cl)
  
  if (verbose) cat("Parallel backend stopped.\n")
  
  # End timing
  end_time <- Sys.time()
  duration <- end_time - start_time
  print(paste("TMLE analysis completed in", duration))
  
  # Process results to separate successful computations and failed iterations
  estimates <- numeric(length(data_list))
  variances <- numeric(length(data_list))
  error_messages <- list()  # Initialize a list to store error messages
  failed_iterations <- integer()  # Initialize as integer vector
  
  for (i in 1:length(tmle_results)) {
    if (tmle_results[[i]]$success) {
      estimates[i] <- tmle_results[[i]]$estimate
      variances[i] <- tmle_results[[i]]$variance
    } else {
      failed_iterations <- c(failed_iterations, i)
      error_messages[[as.character(i)]] <- tmle_results[[i]]$error_msg  # Store error message with iteration index as name
      estimates[i] <- NA  # Or you can choose to exclude it
      variances[i] <- NA
    }
  }
  
  # Garbage collection
  gc()
  
  return(list(
    estimates = estimates, 
    variances = variances,
    failed = failed_iterations, 
    error_messages = error_messages
  ))
}


retry_failed_MI_Int <- function(original_result, data_list, m = 10,
                                          DGP = "DGP1", truncation = TRUE, verbose = TRUE, maxit = 10) {
  # Check if there are any failed iterations
  if (length(original_result$failed) == 0) {
    if (verbose) cat("No failed iterations to retry.\n")
    return(original_result)
  }
  
  # Initialize
  failed_indices <- original_result$failed
  
  while (length(failed_indices) > 0) {
    num_failed <- length(failed_indices)
    if (verbose) cat("Retrying", num_failed, "failed iterations.\n")
    
    # Sample replacement datasets for each failed iteration
    sampled_data_list <- lapply(1:num_failed, function(x) {
      sample_index <- sample(1:length(data_list), 1)
      data_list[[sample_index]]
    })
    cores <- floor(detectCores(logical = FALSE)/2)
    # Run TMLE imputations on the sampled datasets
    retry_result <- run_MI_Int(
      data_list = sampled_data_list,
      m = m,
      cores = cores,
      DGP = DGP,
      truncation = truncation,
      verbose = FALSE,
      maxit = maxit
    )
    
    # Update the original_result with the retry results
    for (i in 1:num_failed) {
      iter_index <- failed_indices[i]
      if (!is.na(retry_result$estimates[i])) {
        # Successful retry: Update the estimates and variances
        original_result$estimates[iter_index] <- retry_result$estimates[i]
        original_result$variances[iter_index] <- retry_result$variances[i]
        
        # Remove from failed
        original_result$failed <- setdiff(original_result$failed, iter_index)
        
        # Remove the corresponding error message
        original_result$error_messages[[as.character(iter_index)]] <- NULL
      } else {
        # Retry also failed: Update the error message
        original_result$error_messages[[as.character(iter_index)]] <- retry_result$error_messages[[as.character(i)]]
        # The iteration remains in the failed list
      }
    }
    
    # Update failed_indices
    failed_indices <- original_result$failed
    
    # Inform the user about the retry status
    if (verbose) {
      num_remaining_failed <- length(failed_indices)
      if (num_remaining_failed > 0) {
        cat(num_remaining_failed, "iterations still failed after retry.\n")
      } else {
        cat("All failed iterations succeeded after retry.\n")
      }
    }
  }
  
  # Perform garbage collection to free up memory
  gc()
  
  return(original_result)
}


merged_MI_Int_TMLE <- function(data_list, m = 10, cores = 4, DGP = "DGP1", truncation = TRUE,
                               verbose = TRUE, maxit = 10) {
  
  first_stage <- run_MI_Int(data_list = data_list, m = m, cores = cores,
                                           DGP = DGP, truncation = truncation, verbose = verbose, maxit = maxit)
  
  second_stage <- retry_failed_MI_Int(original_result = first_stage, data_list = data_list, m = m, DGP = DGP,
                                                truncation = truncation, verbose = verbose, maxit = maxit)
  
  imp_tmle_results <- vector("list", length = length(data_list))
  for (i in 1:length(imp_tmle_results)) {
    imp_tmle_results[[i]] <- list(
      estimate = second_stage$estimates[i],
      variance = second_stage$variances[i])
  }
  
  return(imp_tmle_results) 
}  

apply_all_methods <- function(data_list, m = 10, cores = 4, DGP = "DGP1", truncation = TRUE, maxit = 10) {
  print("CC")
  CC_res <- merged_CC_TMLE(data_list, cores, truncation = truncation)
  
  print("Ext")
  Ext_res <- merged_Ext_TMLE(data_list, cores, DGP, truncation = truncation)
  
  print("Ext_MCMI")
  Ext_MCMI_res <- merged_Ext_MCMI_TMLE(data_list, cores, DGP, truncation = truncation)
  
  print("PMM")
  MI_PMM_res <- merged_MI_PMM_TMLE(data_list, m, cores, DGP, truncation = truncation, maxit = maxit)
  
  print("CART")
  MI_CART_res <- merged_MI_CART_TMLE(data_list, m, cores, DGP, truncation = truncation, maxit = maxit)
  
  print("RF")
  MI_RF_res <- merged_MI_RF_TMLE(data_list, m, cores, DGP, truncation = truncation, maxit = maxit)
  
  print("Amelia")
  MI_Amelia_res <- merged_MI_Amelia_TMLE(data_list, m, cores, DGP, truncation = truncation, maxit = maxit)
  
  print("Int")
  MI_Int_res <- merged_MI_Int_TMLE(data_list, m, cores, DGP, truncation = truncation, maxit = maxit)
  
  result_list <- list(
    CC_res = CC_res,
    Ext_res = Ext_res,
    Ext_MCMI_res = Ext_MCMI_res,
    MI_PMM_res = MI_PMM_res,
    MI_CART_res = MI_CART_res,
    MI_RF_res = MI_RF_res,
    MI_Amelia_res = MI_Amelia_res,
    MI_Int_res = MI_Int_res
  )
  return(result_list)
}
