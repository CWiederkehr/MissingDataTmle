
## Acknowledgement: Haodong Li
# original functions of 'under_HAL', 'get_simu_data', 'ss_estimator', 'get_maxscore', 'calc_ate', 
# 'num_knots_generator', 'tune_knots' 
# can be found in: https://github.com/HaodongL/realistic_simu/blob/master/R/2_simulation_HAL.R

## corresponding Paper: https://doi.org/10.1002/sim.9348


# subetting helper function
"%w/o%" <- function(x, y) x[!x %in% y]

# -----------------------------------------------------------------------------
# Undersoomthed HAL function
# -----------------------------------------------------------------------------
# inputs: 
#       df is a dataframe
#       yname is a string
#       xname is a string
#       Nlam is a scalar, number of candidates lambda
#       type is a string, type of y
under_HAL <- function(df, yname, xname, Nlam, type = "gaussian", smooth = 1, base_knots_0 = 500, degree_crit = 20){
  
  # variables
  y <- as.numeric(as.matrix(df %>% select(all_of(yname))))
  x <- df %>% select(all_of(xname)) %>% 
    mutate_if(sapply(., is.factor), as.numeric)
  n <- nrow(df)
  
  # get initial fit
  print("---initial fitting---")
  tic()
  fit <- fit_hal(X = x,
                 Y = y, 
                 return_x_basis = TRUE,
                 family = type,
                 num_knots = num_knots_generator(
                   max_degree = ifelse(ncol(x) >= degree_crit, 2, 3),
                   smoothness_orders = smooth,
                   base_num_knots_0 = base_knots_0,
                   base_num_knots_1 = max(100, ceiling(sqrt(n)))
                 )
  )
  toc()
  
  # only non-zero direction
  init_coef <-fit$coefs[-1]
  nonzero_col <- which(init_coef != 0)
  
  # if all coefs are zero, skip undersmooth and use the initial fit
  if (length(nonzero_col) == 0){
    res <- list("lambda_init" = fit$lambda_star,
                "lambda_under"= fit$lambda_star)
  }else{
    # refit on new lambda sequence
    print("---candidates fitting---")
    tic()
    us_lambda <- fit$lambda_star*10^seq(from=0, to=-3, length=Nlam)
    us_fit <- glmnet(fit$x_basis, y, lambda=us_lambda, family = type, standardize = FALSE)
    toc()
    
    # evaluate refits
    if (type == "gaussian"){
      pred_mat <- predict(us_fit, fit$x_basis)
      preds_init <- predict(fit, new_data = x)
    }else if (type == "binomial"){
      pred_mat <- predict(us_fit, fit$x_basis, type = "response")
      preds_init <- predict(fit, new_data = x) 
    }
    
    resid_mat <- pred_mat-y
    basis_mat <- as.matrix(fit$x_basis)
    basis_mat <- as.matrix(basis_mat[, nonzero_col])
    
    # estimates of sd in each direction using initial fit
    resid_init <- preds_init-y
    sd_est <- rep(NA, ncol(basis_mat))
    
    for (i in 1:ncol(basis_mat)){
      u <- basis_mat[,i]
      score_init <- resid_init * u
      sd_est[i] <- sd(score_init)
    }
    
    # check the criterion 
    print("---checking criterion---")
    tic()
    max_score <- get_maxscore(basis_mat = basis_mat, 
                              resid_mat = resid_mat,
                              sd_est = sd_est, 
                              Nlam = Nlam, us_fit = us_fit)
    toc()
    
    # get the first lambda that satisfies the criteria
    lambda_under <- us_lambda[max_score <= 1/(sqrt(n)*log(n))][1]
    
    # collect and save results
    coef_mat <- as.matrix(us_fit$beta)
    
    df_under <- data.frame("lambda" = NA,
                           "l1_norm" = NA,
                           "n_coef" = NA)
    
    for (j in 1:Nlam){
      df_under[j,1] = us_lambda[j]
      df_under[j,2] = sum(abs(coef_mat[,j]))
      df_under[j,3] = sum(coef_mat[,j] != 0)
    }
    
    res <- list("lambda_init" = fit$lambda_star,
                "lambda_under" = lambda_under,
                "df_under" = df_under)
  }
  return(res)
}

# # make plot
# plot(x = df_under$lambda, y = df_under$l1_norm)
# abline(v = lambda_under, col="red", lwd=3, lty=2)


# -----------------------------------------------------------------------------
# Simulate data
# -----------------------------------------------------------------------------
# inputs: 
#       df is a dataframe
#       w is a set of covariates' names 
#       a is a string, the name of treatment variable
#       y is a string, the name of outcome variable
#       type is a string, outcome type
#       g_fit is the undersmoothed HAL fit of g model
#       Q_fit is the undersmoothed HAL fit of Q model
#       rv is a scalar, the residual variance

get_simu_data <- function(df, w, a, y, g_fit, Q_fit, rv, weights = NULL, Size){
  
  df$id <- 1:(nrow(df)*Size)
  
  # If weights are NULL, create a uniform weights vector
  if (is.null(weights)) {
    weights <- rep(1 / (nrow(df)*Size), nrow(df)*Size)
  }
  
  simu_data <- df[sample(1:nrow(df), size = nrow(df)*Size, prob = weights, replace = TRUE), ]
  
  # Filter for observations drawn more than 3 times
  counts <- table(simu_data$id)
  more_than_four <- as.data.frame(counts)
  colnames(more_than_four) <- c("id", "count")
  more_than_four$id <- as.numeric(more_than_four$id)
  simu_data <- simu_data %>% select(-id)
  
  # A Bootstrapping: undersmoothed HAL
  
  # generate predictions for A using g fit on original data and binomial sampling
  covariates <- simu_data %>% select(all_of(w)) %>% 
    mutate_if(sapply(., is.factor), as.numeric)
  a_preds <- predict(g_fit, new_data = covariates)
  
  # Save bootstrapped intervention A to dataframe
  simu_data$a <- rbinom(length(a_preds), 1, prob = a_preds)
  
  # Y Bootstrapping: Super Learner
  # generate Y predictions using Q fit on original data and normal error
  covariates <- simu_data %>% select(all_of(c(w, a))) %>% 
    mutate_if(sapply(., is.factor), as.numeric)
  y_preds <- predict(Q_fit, new_data = covariates)
  y_preds_error <- rnorm(length(y_preds), mean = 0, sd = sqrt(rv))
  y_preds <- y_preds + y_preds_error
  
  # Save bootstrapped outcome Y to dataframe
  simu_data$haz <- y_preds
  
  
  # Fit SL for g and Q on bootstrap data (for influence curves)
  
  #0630, drop constant cols
  drop_cov = 0
  for (covs in names(simu_data) %w/o% c('haz', 'a')) {
    if (length(unique(simu_data[[covs]])) == 1) {
      simu_data = simu_data %>% select(-all_of(covs))
      drop_cov = drop_cov + 1
    }
  }
  w =  colnames(simu_data) %w/o% c('haz', 'a')
  
  output <- list('data' = simu_data, 'drop_cov'=drop_cov, 'more_than_four' = more_than_four) 
  return(output)
}


# -----------------------------------------------------------------------------
# Simple substitution estimator for true ATE
# -----------------------------------------------------------------------------
ss_estimator <- function(df, w, a, y){
  # Fit undersmoothed HAL for g on original data
  res_g <- under_HAL(df = df, yname = a, xname = w, Nlam = 100, type = "binomial", smooth = 0, base_knots_0 = 300)
  
  # Preparing covariates for g
  covariates <- df %>% select(all_of(w)) %>% 
    mutate_if(sapply(., is.factor), as.numeric)
  n <- nrow(df)
  
  # Fitting HAL model for g
  g_fit <- fit_hal(X = covariates,
                   Y = as.numeric(as.matrix(df %>% select(all_of(a)))),
                   family = "binomial",
                   num_knots = num_knots_generator(
                     max_degree = ifelse(ncol(covariates) >= 20, 2, 3),
                     smoothness_orders = 0,
                     base_num_knots_0 = 300,
                     base_num_knots_1 = 100
                   ),
                   fit_control = list(
                     cv_select = FALSE,
                     use_min = TRUE,
                     lambda.min.ratio = 1e-4,
                     prediction_bounds = "default"
                   ),
                   lambda = res_g$lambda_under,
                   return_lasso = FALSE
  )
  
  # Fit undersmoothed HAL for Q on original data
  res_Q <- under_HAL(df = df, yname = y, xname = c(w, a), Nlam = 100, type = "gaussian", smooth = 0, base_knots_0 = 300)
  
  # Preparing covariates for Q and fitting HAL model for Q
  covariates <- df %>% select(all_of(c(w, a))) %>% 
    mutate_if(sapply(., is.factor), as.numeric)
  Q_fit <- fit_hal(X = covariates,
                   Y = as.numeric(as.matrix(df %>% select(all_of(y)))),
                   family = "gaussian",
                   num_knots = num_knots_generator(
                     max_degree = ifelse(ncol(covariates) >= 20, 2, 3),
                     smoothness_orders = 0,
                     base_num_knots_0 = 300,
                     base_num_knots_1 = 100
                   ),
                   fit_control = list(
                     cv_select = FALSE,
                     use_min = TRUE,
                     lambda.min.ratio = 1e-4,
                     prediction_bounds = "default"
                   ),
                   lambda = res_Q$lambda_under,
                   return_lasso = FALSE
  )
  
  # Calculate true ATE by generating a large sample
  print("---calculate the truth---")
  
  num_iterations <- if (n == 500) {
    500
  } else if (n == 1000) {
    250
  } else if (n == 2000) {
    125
  } else {
    250  # Default if not 500, 1000, or 2000
  }
  
  list_ate <- rep(NA, num_iterations)
  list_rss <- rep(NA, num_iterations)
  
  # Divide and conquer
  tic()
  for (m in 1:num_iterations){
    print(paste0("calc ate: ", m, "/", num_iterations))
    
    res_ate <- calc_ate(df = df, 
                        g_fit = g_fit, 
                        Q_fit = Q_fit, 
                        n_sample = n,  # Updated as per row count
                        w = w,
                        a = a,
                        y = y)
    
    list_ate[m] <- res_ate$psi_ss
    list_rss[m] <- res_ate$rss
  }
  toc()
  
  # Aggregate results
  psi_ss <- mean(list_ate)
  rv <- mean(list_rss)
  
  # Return all models and estimates
  ss <- c('psi' = psi_ss, 'g_fit' = list(g_fit), 'Q_fit' = list(Q_fit), 'rv' = rv)
  return(ss)
}


# -----------------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------------
# undersoomthed HAL helper function
get_maxscore <- function(basis_mat, resid_mat, sd_est, Nlam, us_fit){
  
  score_all <- matrix(NA, nrow = Nlam, 
                      ncol = ncol(basis_mat))
  
  for (i in 1:ncol(basis_mat)){
    u <- basis_mat[,i]
    score_mat <- resid_mat * u / sd_est[i]
    score_mean <- apply(score_mat, 2, mean)
    score_all[,i] <- score_mean
  }
  
  # absolute value
  max_score <- apply(abs(score_all), 1, max)
  return(max_score)
}

# helper function to calculate the true ATE
calc_ate <- function(df, g_fit, Q_fit, n_sample = 5*10^4, w, a, y){
  # sample W from emp
  large_data <- dplyr::sample_n(df, size = n_sample, replace = TRUE)
  
  # generate predictions for A using g HAL fit
  covariates <- large_data %>% select(all_of(w)) %>% 
    mutate_if(sapply(., is.factor), as.numeric)
  
  a_preds <- predict(g_fit, new_data = covariates)
  
  # generate A
  large_data$a <- rbinom(length(a_preds), 1, prob = a_preds)
  
  covariates <- large_data %>% select(all_of(c(w, a))) %>% 
    mutate_if(sapply(., is.factor), as.numeric)
  # generate Y
  df0 <- data.frame(covariates)
  df1 <- data.frame(covariates)
  df0$a = 0
  df1$a = 1
  
  # QbarAW
  y_preds <- predict(Q_fit, new_data = covariates)
  
  # Qbar1W
  covariates1 <- df1 %>% select(all_of(c(w, a))) %>% 
    mutate_if(sapply(., is.factor), as.numeric)
  y1_preds <- predict(Q_fit, new_data = covariates1)
  
  # Qbar0W
  covariates0 <- df0 %>% select(all_of(c(w, a))) %>% 
    mutate_if(sapply(., is.factor), as.numeric)
  y0_preds <- predict(Q_fit, new_data = covariates0)
  
  # ate
  psi_ss <- mean(y1_preds - y0_preds)
  rss <- (sum((y_preds - large_data %>% select(all_of(y)))^2))/n_sample
  
  return(list("psi_ss" = psi_ss,
              "rss" = rss))
}



num_knots_generator <- function(max_degree, smoothness_orders, base_num_knots_0 = 500,
                                base_num_knots_1 = 200) {
  if (all(smoothness_orders > 0)) {
    return(sapply(seq_len(max_degree), function(d) {
      round(base_num_knots_1 / 2^(d - 1))
    }))
  }
  else {
    return(sapply(seq_len(max_degree), function(d) {
      round(base_num_knots_0 / 2^(d - 1))
    }))
  }
}

# tuning the num_knots (very slow)
tune_knots <- function(df, 
                       yname, 
                       xname, 
                       num_knots_candi = c(150, 200, 250, 300, 350),
                       type = "gaussian"){
  # variables
  y <- as.numeric(as.matrix(df %>% select(all_of(yname))))
  x <- df %>% select(all_of(xname)) %>% 
    mutate_if(sapply(., is.factor), as.numeric)
  
  df_tune <- data.frame("num_knots" = NA,
                        "mse" = NA,
                        "n_coef" = NA,
                        "runtime" = NA)
  
  # tuning
  for(i in 1:length(num_knots_candi)){
    # fit
    start_time <- Sys.time()
    fit <- fit_hal(x,
                   y, 
                   return_x_basis = TRUE,
                   family = type,
                   num_knots = num_knots_generator(
                     max_degree = ifelse(ncol(x) >= 20, 2, 3),
                     smoothness_orders = 1,
                     base_num_knots_0 = 500,
                     base_num_knots_1 = num_knots_candi[i]
                   )
    )
    time_elapse <- as.numeric(Sys.time() - start_time, units="mins")
    
    # non-zero coef
    nonzero_col <- which(fit$coefs[-1] != 0)
    
    # performance
    preds <- predict(fit, new_data = x)
    
    df_tune[i, 1] <- num_knots_candi[i]
    df_tune[i, 2] <- mean((preds - y)^2)
    df_tune[i, 3] <- length(nonzero_col)
    df_tune[i, 4] <- time_elapse
  }
  return(df_tune)
}


#### new helper functions ####

count_violations <- function(probs, threshold = 0.01) {
  
  threshold_high <- 1-threshold
  threshold_low <- threshold
  prop_violation <- (sum(probs > threshold_high) +  sum(probs < threshold_low))/length(probs)
  return(prop_violation)
}

calculate_custom_weights <- function(pred_A, real_A = NULL, top_percentage = 0.10, bottom_percentage = 0.05, weight = 5) {
  # Sort the probabilities in ascending order
  sorted_pred <- sort(pred_A)
  
  # Determine the cutoff index for the top 10% highest probabilities
  top_cutoff_index <- ceiling((1 - top_percentage) * length(pred_A))
  top_threshold <- sorted_pred[top_cutoff_index]
  
  # Determine the cutoff index for the bottom 5% lowest probabilities
  bottom_cutoff_index <- ceiling(bottom_percentage * length(pred_A))
  bottom_threshold <- sorted_pred[bottom_cutoff_index]
  
  # Initialize weights to 1
  weights <- rep(1, length(pred_A))
  
  # If real_A is provided, adjust weights based on alignment with pred_A
  if (!is.null(real_A)) {
    for (i in 1:length(pred_A)) {
      if (pred_A[i] >= top_threshold) {
        if (real_A[i] == 1) {
          weights[i] <- weight
        } else {
          weights[i] <- 0
        }
      } else if (pred_A[i] <= bottom_threshold) {
        if (real_A[i] == 0) {
          weights[i] <- weight
        } else {
          weights[i] <- 0
        }
      }
    }
  } else {
    # Assign weights of 5 to probabilities in the top 10% or bottom 5%
    weights <- ifelse(pred_A >= top_threshold | pred_A <= bottom_threshold, weight, 1)
  }
  
  # Normalize the weights
  weights <- weights / sum(weights)
  
  return(weights)
}

# Draw Bootstrap samples
run_DGP_via_HAL <- function(df, w_cols, aname, yname, g_fit, Q_fit, rv, weights = NULL, Nlam = 100, 
                           Sim = 1000, Size = 1, vio = FALSE, total_cores = 20, degree_crit = 10) {
  # Validate inputs
  validate_inputs <- function() {
    if (!is.character(aname) || length(aname) != 1) stop("Error: 'aname' should be a single string.")
    if (!is.character(yname) || length(yname) != 1) stop("Error: 'yname' should be a single string.")
    if (!all(sapply(w_cols, is.character))) stop("Error: 'w_cols' should be a vector of strings.")
    if (!all(w_cols %in% names(df))) stop("Error: Some covariates specified in 'w_cols' are not present in the dataframe.")
    if (!(aname %in% names(df))) stop("Error: Target variable specified in 'aname' is not present in the dataframe.")
    if (!(yname %in% names(df))) stop("Error: Response variable specified in 'yname' is not present in the dataframe.")
  }
  
  validate_inputs()
  
  # Initialize storage vectors and lists
  simu_data_list <- vector("list", length = Sim)
  drop_cov_vector <- numeric(Sim)
  more_than_four <- vector("list", length = Sim)
  
  if (vio) {
    propScore_preds <- vector("list", length = Sim)
    violation <- numeric(Sim)
    error_messages <- vector("list", length = Sim)
  }
  
  # Set up parallel backend
  cores_per_task <- 2
  num_clusters <- total_cores / cores_per_task
  cl <- makeCluster(num_clusters)
  registerDoParallel(cl, cores_per_task)
  
  # Parallel processing
  simulation_results <- foreach(i = 1:Sim, .packages = c("foreach", "doParallel", "glmnet", "hal9001", "dplyr"),
                                .export = c("get_simu_data", "count_violations", "num_knots_generator", "%w/o%")) %dopar% {
                                  simu_output <- get_simu_data(df = df, w = w_cols, a = aname, y = yname, g_fit = g_fit, Q_fit = Q_fit, 
                                                               rv = rv, weights = weights, Size = Size)
                                  
                                  result <- list(data = simu_output$data, drop_cov = simu_output$drop_cov, 
                                                 more_than_four = simu_output$more_than_four)
                                  
                                  if (vio) {
                                    covariates <- simu_output$data[, w_cols]
                                    for (col in w_cols) {
                                      if (is.factor(covariates[[col]])) {
                                        covariates[[col]] <- as.numeric(covariates[[col]])
                                      }
                                    }
                                    fit_result <- tryCatch({
                                      fit_model <- fit_hal(
                                        X = covariates,
                                        Y = as.numeric(simu_output$data[, aname]),
                                        family = "binomial",
                                        smoothness_orders = 1,
                                        num_knots = num_knots_generator(
                                          max_degree = ifelse(ncol(covariates) >= degree_crit, 2, 3),
                                          smoothness_orders = 0,
                                          base_num_knots_0 = 300,
                                          base_num_knots_1 = 100
                                        ),
                                        return_lasso = FALSE
                                      )
                                      prop_scores <- predict(fit_model, new_data = covariates, type = "response")
                                      violation_count <- count_violations(prop_scores)
                                      list(propScore_pred = prop_scores, violation = violation_count, error = NULL)
                                    }, error = function(e) {
                                      list(propScore_pred = NA, violation = NA, error = conditionMessage(e))
                                    })
                                    
                                    result$propScore_pred <- fit_result$propScore_pred
                                    result$violation <- fit_result$violation
                                    result$error <- fit_result$error
                                  }
                                  
                                  result
                                }
  
  stopCluster(cl)
  
  # Unpack the parallel results into the original list structure
  for (i in 1:Sim) {
    simu_data_list[[i]] <- simulation_results[[i]]$data
    drop_cov_vector[i] <- simulation_results[[i]]$drop_cov
    more_than_four[[i]] <- simulation_results[[i]]$more_than_four
    
    if (vio) {
      propScore_preds[[i]] <- simulation_results[[i]]$propScore_pred
      violation[i] <- simulation_results[[i]]$violation
      error_messages[[i]] <- simulation_results[[i]]$error
    }
  }
  
  # Detect presence of numbers in the 'df' name for naming purposes
  df_name <- deparse(substitute(df))
  suffix <- ""
  if (grepl("500", df_name)) {
    suffix <- "500"
  }
  if (grepl("1000", df_name)) {
    suffix <- "1000"
  }
  if (grepl("2000", df_name)) {
    suffix <- "2000"
  }
  
  result_elements <- list(
    simu_data_list = simu_data_list,
    drop_cov_vector = drop_cov_vector,
    more_than_four = more_than_four
  )
  
  names_result_elements <- c(
    paste0("simu_data_list", suffix),
    paste0("drop_cov_vector", suffix),
    paste0("more_than_four", suffix)
  )
  
  if (vio) {
    result_elements$propScore_preds <- propScore_preds
    result_elements$violation <- violation
    result_elements$error_messages <- error_messages
    names_result_elements <- c(names_result_elements, paste0("propScore_preds", suffix), paste0("violation", suffix), paste0("error_messages", suffix))
  }
  
  names(result_elements) <- names_result_elements
  
  return(result_elements)
}
