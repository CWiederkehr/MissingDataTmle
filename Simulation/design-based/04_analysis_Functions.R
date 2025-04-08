CC <- function(data_list, cores = 4) {
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  
  start <- Sys.time()
  
  required_packages <- c("tmle", "SuperLearner", "ranger", "earth", "arm", "glmnet")
  
  if (grepl("500", deparse(substitute(data_list)))) {
    SL_lib <- c("SL.mean","SL.glm","SL.glm.interaction", "SL.bayesglm", "SL.gam", 
                "SL.glmnet", "SL.earth","SL.rpart", "SL.rpartPrune", "SL.ranger","SL.nnet")
  } else {
    SL_lib <- c("SL.glm", "SL.glm.interaction", "SL.earth", "SL.mean")
  }
  
  results <- foreach(i = 1:length(data_list), .packages = required_packages) %dopar% {
    
    data <- data_list[[i]]
    data <- na.omit(data)
    
    if (nrow(data) == 0) {
      stop("Data frame is empty after removing rows with NAs")
    }
    cov_name <- setdiff(names(data), c("a", "haz"))
    
    result <- tmle(Y = data$haz, A = data$a, W = data[, cov_name], family = "gaussian", 
                   Q.SL.library = SL_lib, g.SL.library = SL_lib, cvQinit = FALSE, gbound = 0.01, V.Q=5, V.g=5
    )
    if (!is.list(result$estimates) || is.null(result$estimates$ATE)) {
      stop("tmle output does not have the expected structure")
    }
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

Ext <- function(data_list, cores = 4) {
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  
  start <- Sys.time()
  
  required_packages <- c("tmle", "SuperLearner", "ranger", "earth", "arm", "glmnet")
  
  if (grepl("500", deparse(substitute(data_list)))) {
    SL_lib <- c("SL.mean","SL.glm","SL.glm.interaction", "SL.bayesglm", "SL.gam", 
                "SL.glmnet", "SL.earth","SL.rpart", "SL.rpartPrune", "SL.ranger","SL.nnet")
  } else {
    SL_lib <- c("SL.glm", "SL.glm.interaction", "SL.earth", "SL.mean")
  }
  
  results <- foreach(i = 1:length(data_list), .packages = required_packages) %dopar% {
    data <- data_list[[i]]
    
    # records with missing exposure or confounder data being excluded except 'haz' (-> outcome)
    data <- data[complete.cases(data[, c("a", "enwast", "W_mhtcm", "W_mwtkg", "W_parity")]), ]
    if (nrow(data) == 0) {
      stop("Data frame is empty after removing rows with NAs")
    }
    # in Paper 'tmle: An R Package for Targeted Maximum Likelihood Estimation' 0 is declared to be missing 
    #-> create variable deltaY to estimate probability of having observed outcome conditional on the exposure and confounders
    data$deltaY <- ifelse(is.na(data$haz), 0, 1)
    
    cov_name <- setdiff(names(data), c("a", "haz", "deltaY"))
    
    result <- tmle(Y = data$haz, A = data$a, W = data[, cov_name], family = "gaussian", 
                   Q.SL.library = SL_lib, g.SL.library = SL_lib, 
                   g.Delta.SL.library = SL_lib, Delta = data$deltaY, 
                   cvQinit = FALSE, gbound = 0.01, V.Q=5, V.g=5)
    if (!is.list(result$estimates) || is.null(result$estimates$ATE)) {
      stop("tmle output does not have the expected structure")
    }
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

Ext_MCMI <- function(data_list, cores = 4) {
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  
  start <- Sys.time()
  # Parallel execution using foreach
  required_packages <- c("tmle", "SuperLearner", "ranger", 
                         "earth", "arm", "glmnet")
  
  if (grepl("500", deparse(substitute(data_list)))) {
    SL_lib <- c("SL.mean","SL.glm","SL.glm.interaction", "SL.bayesglm", "SL.gam", 
                "SL.glmnet", "SL.earth","SL.rpart", "SL.rpartPrune", "SL.ranger","SL.nnet")
  } else {
    SL_lib <- c("SL.glm", "SL.glm.interaction", "SL.earth", "SL.mean")
  }
  
  
  results <- foreach(i = 1:(length(data_list)), .packages = required_packages) %dopar% {
    data <- data_list[[i]]
    
    # convert Confounders with missing observations to factors
    data$enwast <- as.factor(data$enwast)
    
    # create new category for 'NAs'
    levels(data$enwast) <- c(levels(data$enwast), "Missing") 
    
    # change 'NAs' to 'Missing'
    data$enwast[is.na(data$enwast)] <- c("Missing")
    
    data$delta_W_mhtcm <- ifelse(is.na(data$W_mhtcm), 0, 1) # In Paper MissingIndicator = 0 -> for missing values
    data$delta_W_mwtkg <- ifelse(is.na(data$W_mwtkg), 0, 1) # In Paper MissingIndicator = 0 -> for missing values
    data$delta_W_parity <- ifelse(is.na(data$W_parity), 0, 1) # In Paper MissingIndicator = 0 -> for missing values
    
    data$W_mhtcm[is.na(data$W_mhtcm)] <- mean(data$W_mhtcm, na.rm = TRUE)  
    data$W_mwtkg[is.na(data$W_mwtkg)] <- mean(data$W_mwtkg, na.rm = TRUE)  
    data$W_parity[is.na(data$W_parity)] <- 0  
    # exclude only observations with unobserved exposures
    data <- data[!is.na(data$a), ]
    # add deltaY for tmle (zero means y-value is missing)
    data$deltaY <- ifelse(is.na(data$haz), 0, 1) 
    
    cov_name <- names(data)
    cov_name <- cov_name[!(cov_name %in% c("a", "haz", "deltaY"))]
    
    result <- tmle(Y = data$haz, A = data$a, W = data[, cov_name], family = "gaussian", 
                   Q.SL.library = SL_lib, g.SL.library = SL_lib, g.Delta.SL.library = SL_lib, 
                   Delta = data$deltaY, cvQinit = FALSE, gbound = 0.01, V.Q=5, V.g=5)
    
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

MI_PMM <- function(data_list, cores = 4, m = 10, maxit = 10) {
  
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  
  start <- Sys.time()
  
  # Specify the packages required for each worker
  required_packages <- c("tmle", "SuperLearner", "mice", "ranger", "earth", "arm", "glmnet")
  
  if (grepl("500", deparse(substitute(data_list)))) {
    SL_lib <- c("SL.mean","SL.glm","SL.glm.interaction", "SL.bayesglm", "SL.gam", 
                "SL.glmnet", "SL.earth","SL.rpart", "SL.rpartPrune", "SL.ranger","SL.nnet")
  } else {
    SL_lib <- c("SL.glm", "SL.glm.interaction", "SL.earth", "SL.mean")
  }
  
  results <- foreach(i = seq_along(data_list), .packages = required_packages) %dopar% {
    data <- data_list[[i]]
    cov_name <- setdiff(names(data), c("a", "haz"))
    
    # Convert only necessary columns to factors
    data[["month"]] <- as.factor(data[["month"]])
    data[["brthmon"]] <- as.factor(data[["brthmon"]])
    data[["enwast"]] <- as.factor(data[["enwast"]])
    data[["a"]] <- as.factor(data[["a"]])
    
    # specify imputation methods and conduct imputation
    PMM_imp <- mice(data, method = c("", "", "", "", "", "", "", "", "logreg", "pmm", "",
                                     "pmm", "pmm", "", "logreg", "", "pmm"), m = m, maxit = maxit)
    
    Estimates <- numeric(m) # m = 10
    ModVar <- numeric(m)
    
    for(l in 1:m) {
      imp_data <- mice::complete(PMM_imp, l)
      # convert back to type factor
      imp_data[["month"]] <- as.integer(imp_data[["month"]])
      imp_data[["brthmon"]] <- as.integer(imp_data[["brthmon"]])
      imp_data[["enwast"]] <- as.integer(imp_data[["enwast"]]) - 1
      imp_data[["a"]] <- as.integer(imp_data[["a"]]) - 1
      
      tmle_results <- tmle(Y = imp_data$haz, A = imp_data$a, 
                           W = imp_data[, cov_name],
                           family = "gaussian", Q.SL.library = SL_lib,
                           g.SL.library = SL_lib,
                           cvQinit = FALSE, gbound = 0.01, V.Q=5, V.g=5)
      Estimates[l] <- tmle_results$estimates$ATE$psi
      ModVar[l] <- tmle_results$estimates$ATE$var.psi
      rm(imp_data, tmle_results)
      gc()
    }
    res <- list(
      mean(Estimates),
      mean(ModVar) + (1 + 1/m) * ((1/(m-1)) * sum((Estimates - mean(Estimates))^2))
    )
    res
  }
  stopCluster(cl)
  print(Sys.time() - start)
  gc()
  return(results)
}

MI_CART <- function(data_list, cores = 4, m = 10, maxit = 10) {
  
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  
  start <- Sys.time()
  
  # Specify the packages required for each worker
  required_packages <- c("tmle", "SuperLearner", "mice", "ranger", "earth", "arm", "glmnet")
  
  if (grepl("500", deparse(substitute(data_list)))) {
    SL_lib <- c("SL.mean","SL.glm","SL.glm.interaction", "SL.bayesglm", "SL.gam", 
                "SL.glmnet", "SL.earth","SL.rpart", "SL.rpartPrune", "SL.ranger","SL.nnet")
  } else {
    SL_lib <- c("SL.glm", "SL.glm.interaction", "SL.earth", "SL.mean")
  }
  
  results <- foreach(i = seq_along(data_list), .packages = required_packages) %dopar% {
    data <- data_list[[i]]
    cov_name <- setdiff(names(data), c("a", "haz"))
    
    # Convert only necessary columns to factors
    data[["month"]] <- as.factor(data[["month"]])
    data[["brthmon"]] <- as.factor(data[["brthmon"]])
    data[["enwast"]] <- as.factor(data[["enwast"]])
    data[["a"]] <- as.factor(data[["a"]])
    
    # specify imputation methods and conduct imputation
    Cart_imp <- mice(data, method = c("", "", "", "", "", "", "", "", "cart", "cart", "",
                                      "cart", "cart", "", "cart", "", "cart"), m = m, maxit = maxit)
    
    Estimates <- numeric(m) # m = 10
    ModVar <- numeric(m)
    
    for(l in 1:m) {
      imp_data <- mice::complete(Cart_imp, l)
      # convert back to type factor
      imp_data[["month"]] <- as.integer(imp_data[["month"]])
      imp_data[["brthmon"]] <- as.integer(imp_data[["brthmon"]])
      imp_data[["enwast"]] <- as.integer(imp_data[["enwast"]]) - 1
      imp_data[["a"]] <- as.integer(imp_data[["a"]]) - 1
      
      tmle_results <- tmle(Y = imp_data$haz, A = imp_data$a, 
                           W = imp_data[, cov_name],
                           family = "gaussian", Q.SL.library = SL_lib,
                           g.SL.library = SL_lib,
                           cvQinit = FALSE, gbound = 0.01, V.Q=5, V.g=5)
      Estimates[l] <- tmle_results$estimates$ATE$psi
      ModVar[l] <- tmle_results$estimates$ATE$var.psi
      rm(imp_data, tmle_results)
      gc()
    }
    res <- list(
      mean(Estimates),
      mean(ModVar) + (1 + 1/m) * ((1/(m-1)) * sum((Estimates - mean(Estimates))^2))
    )
    res
  }
  stopCluster(cl)
  print(Sys.time() - start)
  gc()
  return(results)
}

MI_RF <- function(data_list, cores = 4, m = 10, maxit = 10) {
  
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  
  start <- Sys.time()
  
  # Specify the packages required for each worker
  required_packages <- c("tmle", "SuperLearner", "mice", "ranger", "earth", "arm", "glmnet")
  
  if (grepl("500", deparse(substitute(data_list)))) {
    SL_lib <- c("SL.mean","SL.glm","SL.glm.interaction", "SL.bayesglm", "SL.gam", 
                "SL.glmnet", "SL.earth","SL.rpart", "SL.rpartPrune", "SL.ranger","SL.nnet")
  } else {
    SL_lib <- c("SL.glm", "SL.glm.interaction", "SL.earth", "SL.mean")
  }
  
  results <- foreach(i = seq_along(data_list), .packages = required_packages) %dopar% {
    data <- data_list[[i]]
    cov_name <- setdiff(names(data), c("a", "haz"))
    
    # Convert only necessary columns to factors
    data[["month"]] <- as.factor(data[["month"]])
    data[["brthmon"]] <- as.factor(data[["brthmon"]])
    data[["enwast"]] <- as.factor(data[["enwast"]])
    data[["a"]] <- as.factor(data[["a"]])
    
    # specify imputation methods and conduct imputation
    RF_imp <- mice(data, method = c("", "", "", "", "", "", "", "", "rf", "rf", "",
                                    "rf", "rf", "", "rf", "", "rf"), m = m, maxit = maxit)
    
    Estimates <- numeric(m) # m = 10
    ModVar <- numeric(m)
    
    for(l in 1:m) {
      imp_data <- mice::complete(RF_imp, l)
      # convert back to type factor
      imp_data[["month"]] <- as.integer(imp_data[["month"]])
      imp_data[["brthmon"]] <- as.integer(imp_data[["brthmon"]])
      imp_data[["enwast"]] <- as.integer(imp_data[["enwast"]]) - 1
      imp_data[["a"]] <- as.integer(imp_data[["a"]]) - 1
      
      tmle_results <- tmle(Y = imp_data$haz, A = imp_data$a, 
                           W = imp_data[, cov_name],
                           family = "gaussian", Q.SL.library = SL_lib,
                           g.SL.library = SL_lib,
                           cvQinit = FALSE, gbound = 0.01, V.Q=5, V.g=5)
      Estimates[l] <- tmle_results$estimates$ATE$psi
      ModVar[l] <- tmle_results$estimates$ATE$var.psi
      rm(imp_data, tmle_results)
      gc()
    }
    res <- list(
      mean(Estimates),
      mean(ModVar) + (1 + 1/m) * ((1/(m-1)) * sum((Estimates - mean(Estimates))^2))
    )
    res
  }
  stopCluster(cl)
  print(Sys.time() - start)
  gc()
  return(results)
}

MI_Int <- function(data_list, cores = 4, m = 10, maxit = 10) {
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  
  start <- Sys.time()
  # Parallel execution using foreach
  required_packages <- c("tmle", "SuperLearner", "mice", "ranger", 
                         "earth", "arm", "glmnet")
  
  results <- foreach(i = 1:(length(data_list)), .packages = required_packages) %dopar% {
    data <- data_list[[i]]
    cov_name <- names(data)
    cov_name <- cov_name[!(cov_name %in% c("a", "haz"))]
    # Add new columns with all 2-way interactions including W_feducyrs
    data <- cbind(data,
                  days_cm = data$agedays * data$W_mhtcm,
                  days_kg = data$agedays * data$W_mwtkg,
                  cm_kg = data$W_mhtcm * data$W_mwtkg,
                  days_mage = data$agedays * data$W_mage,
                  mage_cm = data$W_mage * data$W_mhtcm,
                  mage_kg = data$W_mage * data$W_mwtkg,
                  days_feducyrs = data$agedays * data$W_feducyrs,
                  cm_feducyrs = data$W_mhtcm * data$W_feducyrs,
                  kg_feducyrs = data$W_mwtkg * data$W_feducyrs,
                  mage_feducyrs = data$W_mage * data$W_feducyrs)
    
    # Convert specific columns to factor
    data$enwast <- as.factor(data$enwast)
    data$a <- as.factor(data$a)
    
    # Initialize mice without imputation to get methods and predictorMatrix
    template_mice <- mice(data, maxit = 0, print = FALSE)
    
    # Setup specific methods for imputation
    methods <- template_mice$method
    methods[c("a", "enwast")] <- "logreg"
    methods[colnames(data)[18:27]] <- c("~I(agedays*W_mhtcm)", 
                                        "~I(agedays*W_mwtkg)",
                                        "~I(W_mhtcm*W_mwtkg)",
                                        "~I(agedays*W_mage)",
                                        "~I(W_mage*W_mhtcm)",
                                        "~I(W_mage*W_mwtkg)",
                                        "~I(agedays*W_feducyrs)",
                                        "~I(W_mhtcm*W_feducyrs)",
                                        "~I(W_mwtkg*W_feducyrs)",
                                        "~I(W_mage*W_feducyrs)")
    
    # Setup predictor matrix
    predictions <- template_mice$predictorMatrix
    # Example adjustments, replace with actual requirements
    predictions[unique(template_mice$loggedEvents$out), ] <- 1
    predictions[18:27, ] <- 0  # Assuming the new columns are at these positions
    
    # Additional example adjustments
    predictions[c("W_mhtcm", "W_mwtkg", "enwast", "W_parity"), c("days_cm", "days_kg", "cm_kg", "days_mage",
                                                                 "mage_cm", "mage_kg", "days_feducyrs", "cm_feducyrs",
                                                                 "kg_feducyrs", "mage_feducyrs")] <- 0
    
    m <- 10
    # Perform imputation with the specified methods and predictor matrix
    Int_imp <- mice(data, meth = methods, pred = predictions, m = m, maxit = maxit)
    
    # Create a list for imputed data and TMLE estimations 
    Int_imp_data_list <- vector("list", length = m)
    Int_imp_TMLE <- vector("list", length = m)
    ModVar <- vector(length = m)
    Estimate <- vector(length = m)
    
    for(l in 1:m) {
      Int_imp_data_list[[l]] <- mice::complete(Int_imp, l)[,1:17] # mice::complete()
      # Convert factors to integers, adjust as necessary
      Int_imp_data_list[[l]]$enwast <- as.integer(Int_imp_data_list[[l]]$enwast) -1
      Int_imp_data_list[[l]]$a <- as.integer(Int_imp_data_list[[l]]$a) -1
      
      Int_imp_TMLE[[l]] <- tmle(Y = Int_imp_data_list[[l]]$haz, A = Int_imp_data_list[[l]]$a, 
                                W = Int_imp_data_list[[l]][, cov_name],
                                family = "gaussian", Q.SL.library = c("SL.glm", "SL.glm.interaction", "SL.earth", "SL.mean"),
                                g.SL.library = c("SL.glm", "SL.glm.interaction", "SL.earth", "SL.mean"),
                                cvQinit = FALSE, gbound = 0.01)
      ModVar[l] <- Int_imp_TMLE[[l]]$estimates$ATE$var.psi
      Estimate[l] <- (Int_imp_TMLE[[l]]$estimates$ATE$psi)
    }
    
    imp_tmle_results <- vector("list", length = 2)
    imp_tmle_results[[1]] <- mean(Estimate)
    # total variance = average of the complete-data variances + the variance between the m complete-data estimates
    # + extra simulation variance  -> https://stefvanbuuren.name/fimd/sec-whyandwhen.html    Chapter 2.3.2
    imp_tmle_results[[2]] <- mean(ModVar) + (1 + 1/m)*((1/(m-1))*(sum((Estimate-mean(Estimate))^2)))
    
    imp_tmle_results
  }
  stopCluster(cl)
  end <- Sys.time()
  duration <- end - start
  print(duration)
  return(results)
}

MI_Amelia <- function(data_list, cores = 4, m = 10, maxit = 10) {
  
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  
  start <- Sys.time()
  
  # Specify the packages required for each worker
  required_packages <- c("tmle", "SuperLearner", "Amelia", "ranger", "earth", "arm", "glmnet")
  
  if (grepl("500", deparse(substitute(data_list)))) {
    SL_lib <- c("SL.mean","SL.glm","SL.glm.interaction", "SL.bayesglm", "SL.gam", 
                "SL.glmnet", "SL.earth","SL.rpart", "SL.rpartPrune", "SL.ranger","SL.nnet")
  } else {
    SL_lib <- c("SL.glm", "SL.glm.interaction", "SL.earth", "SL.mean")
  }
  
  results <- foreach(i = seq_along(data_list), .packages = required_packages) %dopar% {
    data <- data_list[[i]]
    cov_name <- setdiff(names(data), c("a", "haz"))
    
    # Convert only necessary columns to factors
    data[["month"]] <- as.factor(data[["month"]])
    data[["brthmon"]] <- as.factor(data[["brthmon"]])
    data[["W_parity"]] <- as.factor(data[["W_parity"]])
    data[["enwast"]] <- as.factor(data[["enwast"]])
    data[["a"]] <- as.factor(data[["a"]])
    
    # specify factors and conduct imputation
    Amelia_imp <- amelia(data, m = m, noms = c("sex", "enwast", "a"), 
                         ords = c("month", "brthmon", "hfoodsec", "W_parity"), 
                         logs = c("W_mwtkg", "W_mbmi"))
    # -> month and brthmon are ordered because otherwise MI with Amelia has convergence issues ... W_mage will cause collinearity issues
    
    Estimates <- numeric(m) # m = 10
    ModVar <- numeric(m)
    
    for(l in 1:m) {
      imp_data <- Amelia_imp$imputations[[l]]
      # convert back to type factor
      imp_data[["month"]] <- as.integer(imp_data[["month"]])
      imp_data[["brthmon"]] <- as.integer(imp_data[["brthmon"]])
      imp_data[["W_parity"]] <- as.integer(imp_data[["W_parity"]])
      imp_data[["enwast"]] <- as.integer(imp_data[["enwast"]]) - 1
      imp_data[["a"]] <- as.integer(imp_data[["a"]]) - 1
      
      tmle_results <- tmle(Y = imp_data$haz, A = imp_data$a, 
                           W = imp_data[, cov_name],
                           family = "gaussian", Q.SL.library = SL_lib,
                           g.SL.library = SL_lib,
                           cvQinit = FALSE, gbound = 0.01, V.Q=5, V.g=5)
      Estimates[l] <- tmle_results$estimates$ATE$psi
      ModVar[l] <- tmle_results$estimates$ATE$var.psi
      rm(imp_data, tmle_results)
      gc()
    }
    res <- list(
      mean(Estimates),
      mean(ModVar) + (1 + 1/m) * ((1/(m-1)) * sum((Estimates - mean(Estimates))^2))
    )
    res
  }
  stopCluster(cl)
  print(Sys.time() - start)
  gc()
  return(results)
}


apply_all_methods <- function(data_list, m = 10, cores = 4, maxit = 10) {
  print("CC")
  CC_res <- CC(data_list, cores)
  
  print("Ext")
  Ext_res <- Ext(data_list, cores)
  
  print("Ext_MCMI")
  Ext_MCMI_res <- Ext_MCMI(data_list, cores)
  
  print("PMM")
  MI_PMM_res <- MI_PMM(data_list, cores, m = m, maxit = 10)
  
  print("CART")
  MI_CART_res <- MI_CART(data_list, cores, m = m, maxit = 10)
  
  print("RF")
  MI_RF_res <- MI_RF(data_list, cores, m = m, maxit = 10)
  
  print("Amelia")
  MI_Amelia_res <- MI_Amelia(data_list, cores, m = m)
  
  print("Int")
  MI_Int_res <- MI_Int(data_list, cores, m = m, maxit = 10)
  
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