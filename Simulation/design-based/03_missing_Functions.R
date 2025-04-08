generate_Missing_Data <- function(data_list, MissingIndicator_list) {
  
  for (i in 1:length(data_list)) {
    data_list[[i]]$haz[MissingIndicator_list[[i]]$MY == 1] <- NA
    data_list[[i]]$a[MissingIndicator_list[[i]]$MA == 1] <- NA
    data_list[[i]]$enwast[MissingIndicator_list[[i]]$MZ1 == 1] <- NA
    data_list[[i]]$W_mhtcm[MissingIndicator_list[[i]]$MZ2 == 1] <- NA
    data_list[[i]]$W_mwtkg[MissingIndicator_list[[i]]$MZ3 == 1] <- NA
    data_list[[i]]$W_parity[MissingIndicator_list[[i]]$MZ4 == 1] <- NA
  }
  return(data_list)
}

generate_mDag_Data <- function(data_list, coefs, m_dag) {
  
  # Fixed coefficients for different m-DAGs with real variable names
  fixed_coefs <- list(
    E = list(
      Z_enwast = 0.6,
      Z_W_mhtcm = 0.1,
      Z_W_mwtkg = 0.1,
      Z_W_parity = 0.1,
      Z_sex = -0.6,
      Z_agedays = 0.1,
      Z_W_meducyrs = -0.1,
      Z_a = -0.6,
      Z_haz = 0.1
    ),
    D = list(
      Z_enwast = 0.6,
      Z_W_mhtcm = 0.1,
      Z_W_mwtkg = 0.1,
      Z_W_parity = 0.1,
      Z_sex = -0.6,
      Z_agedays = 0.1,
      Z_W_meducyrs = -0.1,
      Z_a = -0.6,
      Z_haz = 0.1
    ),
    C = list(
      Z_enwast = 0.6,
      Z_W_mhtcm = 0.1,
      Z_W_mwtkg = 0.1,
      Z_W_parity = 0.1,
      Z_sex = -0.6,
      Z_agedays = 0.1,
      Z_W_meducyrs = -0.1,
      Z_a = -0.6,
      Z_haz = 0
    ),
    B = list(
      Z_enwast = 0,
      Z_W_mhtcm = 0,
      Z_W_mwtkg = 0,
      Z_W_parity = 0,
      Z_sex = -0.6,
      Z_agedays = 0.1,
      Z_W_meducyrs = -0.1,
      Z_a = 0,
      Z_haz = 0
    ),
    A = list(
      Z_enwast = 0,
      Z_W_mhtcm = 0,
      Z_W_mwtkg = 0,
      Z_W_parity = 0,
      Z_sex = 0,
      Z_agedays = 0,
      Z_W_meducyrs = 0,
      Z_a = 0,
      Z_haz = 0
    )
  )
  
  # Center the variables in data_list
  center_variables <- function(df) {
    df$agedays <- scale(df$agedays, center = TRUE, scale = TRUE)
    df$W_meducyrs <- scale(df$W_meducyrs, center = TRUE, scale = FALSE)
    df$W_mhtcm <- scale(df$W_mhtcm, center = TRUE, scale = FALSE)
    df$W_mwtkg <- scale(df$W_mwtkg, center = TRUE, scale = FALSE)
    df$W_parity <- scale(df$W_parity, center = TRUE, scale = FALSE)
    df$haz <- scale(df$haz, center = TRUE, scale = FALSE)
    df$sex <- ifelse(as.integer(df$sex) == 1, 0, 1)
    return(df)
  }
  
  indicator_list <- lapply(data_list, function(df) {
    df <- center_variables(df)
    
    MZ1 <- rbinom(nrow(df), size = 1, prob = plogis(
      coefs$MZ1$intercept +
        fixed_coefs[[m_dag]]$Z_enwast * df$enwast +
        fixed_coefs[[m_dag]]$Z_sex * df$sex +
        fixed_coefs[[m_dag]]$Z_agedays * df$agedays +
        fixed_coefs[[m_dag]]$Z_W_meducyrs * df$W_meducyrs +
        fixed_coefs[[m_dag]]$Z_a * df$a +
        fixed_coefs[[m_dag]]$Z_haz * df$haz
    ))
    df$MZ1 <- MZ1
    
    MZ2 <- rbinom(nrow(df), size = 1, prob = plogis(
      coefs$MZ2$intercept +
        fixed_coefs[[m_dag]]$Z_W_mhtcm * df$W_mhtcm +
        fixed_coefs[[m_dag]]$Z_sex * df$sex +
        fixed_coefs[[m_dag]]$Z_agedays * df$agedays +
        fixed_coefs[[m_dag]]$Z_W_meducyrs * df$W_meducyrs +
        fixed_coefs[[m_dag]]$Z_a * df$a +
        fixed_coefs[[m_dag]]$Z_haz * df$haz +
        coefs$MZ2$MZ1 * df$MZ1
    ))
    df$MZ2 <- MZ2
    
    MZ3 <- rbinom(nrow(df), size = 1, prob = plogis(
      coefs$MZ3$intercept +
        fixed_coefs[[m_dag]]$Z_W_mwtkg * df$W_mwtkg +
        fixed_coefs[[m_dag]]$Z_sex * df$sex +
        fixed_coefs[[m_dag]]$Z_agedays * df$agedays +
        fixed_coefs[[m_dag]]$Z_W_meducyrs * df$W_meducyrs +
        fixed_coefs[[m_dag]]$Z_a * df$a +
        fixed_coefs[[m_dag]]$Z_haz * df$haz +
        coefs$MZ3$MZ1 * df$MZ1 +
        coefs$MZ3$MZ2 * df$MZ2
    ))
    df$MZ3 <- MZ3
    
    MZ4 <- rbinom(nrow(df), size = 1, prob = plogis(
      coefs$MZ4$intercept +
        fixed_coefs[[m_dag]]$Z_W_parity * df$W_parity +
        fixed_coefs[[m_dag]]$Z_sex * df$sex +
        fixed_coefs[[m_dag]]$Z_agedays * df$agedays +
        fixed_coefs[[m_dag]]$Z_W_meducyrs * df$W_meducyrs +
        fixed_coefs[[m_dag]]$Z_a * df$a +
        fixed_coefs[[m_dag]]$Z_haz * df$haz +
        coefs$MZ4$MZ1 * df$MZ1 +
        coefs$MZ4$MZ2 * df$MZ2 +
        coefs$MZ4$MZ3 * df$MZ3
    ))
    df$MZ4 <- MZ4
    
    MA <- rbinom(nrow(df), size = 1, prob = plogis(
      coefs$MA$intercept +
        fixed_coefs[[m_dag]]$Z_enwast * df$enwast +
        fixed_coefs[[m_dag]]$Z_W_mhtcm * df$W_mhtcm +
        fixed_coefs[[m_dag]]$Z_W_mwtkg * df$W_mwtkg +
        fixed_coefs[[m_dag]]$Z_W_parity * df$W_parity +
        fixed_coefs[[m_dag]]$Z_sex * df$sex +
        fixed_coefs[[m_dag]]$Z_agedays * df$agedays +
        fixed_coefs[[m_dag]]$Z_W_meducyrs * df$W_meducyrs +
        fixed_coefs[[m_dag]]$Z_a * df$a +
        fixed_coefs[[m_dag]]$Z_haz * df$haz +
        coefs$MA$MZ1 * df$MZ1 +
        coefs$MA$MZ2 * df$MZ2 +
        coefs$MA$MZ3 * df$MZ3 +
        coefs$MA$MZ4 * df$MZ4 
    ))
    df$MA <- MA
    
    MY <- rbinom(nrow(df), size = 1, prob = plogis(
      coefs$MY$intercept +
        fixed_coefs[[m_dag]]$Z_enwast * df$enwast +
        fixed_coefs[[m_dag]]$Z_W_mhtcm * df$W_mhtcm +
        fixed_coefs[[m_dag]]$Z_W_mwtkg * df$W_mwtkg +
        fixed_coefs[[m_dag]]$Z_W_parity * df$W_parity +
        fixed_coefs[[m_dag]]$Z_sex * df$sex +
        fixed_coefs[[m_dag]]$Z_agedays * df$agedays +
        fixed_coefs[[m_dag]]$Z_W_meducyrs * df$W_meducyrs +
        fixed_coefs[[m_dag]]$Z_a * df$a +
        (if (m_dag == "I") 0 else fixed_coefs[[m_dag]]$Z_haz * df$haz) +
        coefs$MY$MZ1 * df$MZ1 +
        coefs$MY$MZ2 * df$MZ2 +
        coefs$MY$MZ3 * df$MZ3 +
        coefs$MY$MZ4 * df$MZ4 +
        coefs$MY$MA * df$MA 
    ))
    df$MY <- MY
    
    return(data.frame(MZ1, MZ2, MZ3, MZ4, MA, MY))
  })
  
  missing_data_list <- generate_Missing_Data(data_list, indicator_list)
  
  return(missing_data_list)
}

# wrapper function to generate all m-DAGs
generate_All_mDags <- function(simu_data_list, pos = c("no_pos", "pos1", "pos2", "pos3")) {
  # Choose the setting: "pos2" or "pos3"
  pos <- match.arg(pos)
  
  # Define coefficient lists based on the chosen 'pos'
  if (pos == "no_pos") {
    # Coefficient sets for "no pos"
    coefsA <- list(
      MZ1 = list(intercept = -1.75),
      MZ2 = list(intercept = -1.80, MZ1 = 2),
      MZ3 = list(intercept = -1.65, MZ1 = 2, MZ2 = 2),
      MZ4 = list(intercept = -3.55, MZ1 = 2.1, MZ2 = 2.1, MZ3 = 2.1),
      MA  = list(intercept = -2.75, MZ1 = 2, MZ2 = 2, MZ3 = 2, MZ4 = 2),
      MY  = list(intercept = -2.00, MZ1 = 0.6, MZ2 = 0.6, MZ3 = 0.6, MZ4 = 0.6, MA = -0.2)
    )
    coefsB <- list(
      MZ1 = list(intercept = -1.50),
      MZ2 = list(intercept = -1.60, MZ1 = 2),
      MZ3 = list(intercept = -1.40, MZ1 = 2, MZ2 = 2),
      MZ4 = list(intercept = -3.40, MZ1 = 2.1, MZ2 = 2.1, MZ3 = 2.1),
      MA  = list(intercept = -2.45, MZ1 = 2, MZ2 = 2, MZ3 = 2, MZ4 = 2),
      MY  = list(intercept = -1.65, MZ1 = 0.5, MZ2 = 0.5, MZ3 = 0.5, MZ4 = 0.5, MA = -0.2)
    )
    coefsC <- list(
      MZ1 = list(intercept = -1.15),
      MZ2 = list(intercept = -1.30, MZ1 = 2.2),
      MZ3 = list(intercept = -1.15, MZ1 = 2.2, MZ2 = 2.2),
      MZ4 = list(intercept = -3.05, MZ1 = 2.2, MZ2 = 2.2, MZ3 = 2.2),
      MA  = list(intercept = -2.05, MZ1 = 1.7, MZ2 = 1.7, MZ3 = 1.7, MZ4 = 1.7),
      MY  = list(intercept = -1.20, MZ1 = 0.4, MZ2 = 0.4, MZ3 = 0.4, MZ4 = 0.4, MA = -0.65)
    )
    coefsD <- list(
      MZ1 = list(intercept = -1.15),
      MZ2 = list(intercept = -1.30, MZ1 = 2.2),
      MZ3 = list(intercept = -1.15, MZ1 = 2.2, MZ2 = 2.2),
      MZ4 = list(intercept = -3.05, MZ1 = 2.2, MZ2 = 2.2, MZ3 = 2.2),
      MA  = list(intercept = -2.05, MZ1 = 1.7, MZ2 = 1.7, MZ3 = 1.7, MZ4 = 1.7),
      MY  = list(intercept = -1.20, MZ1 = 0.4, MZ2 = 0.4, MZ3 = 0.4, MZ4 = 0.4, MA = -0.65)
    )
    coefsE <- list(
      MZ1 = list(intercept = -1.15),
      MZ2 = list(intercept = -1.30, MZ1 = 2.2),
      MZ3 = list(intercept = -1.15, MZ1 = 2.2, MZ2 = 2.2),
      MZ4 = list(intercept = -3.05, MZ1 = 2.2, MZ2 = 2.2, MZ3 = 2.2),
      MA  = list(intercept = -2.05, MZ1 = 1.7, MZ2 = 1.7, MZ3 = 1.7, MZ4 = 1.7),
      MY  = list(intercept = -1.20, MZ1 = 0.4, MZ2 = 0.4, MZ3 = 0.4, MZ4 = 0.4, MA = -0.65)
    )
    
  } else if (pos == "pos1") {
    # Coefficient sets for pos1
    coefsA <- list(
      MZ1 = list(intercept = -1.75),
      MZ2 = list(intercept = -1.80, MZ1 = 2),
      MZ3 = list(intercept = -1.65, MZ1 = 2, MZ2 = 2),
      MZ4 = list(intercept = -3.55, MZ1 = 2.1, MZ2 = 2.1, MZ3 = 2.1),
      MA  = list(intercept = -2.75, MZ1 = 2, MZ2 = 2, MZ3 = 2, MZ4 = 2),
      MY  = list(intercept = -2.00, MZ1 = 0.6, MZ2 = 0.6, MZ3 = 0.6, MZ4 = 0.6, MA = -0.2)
    )
    coefsB <- list(
      MZ1 = list(intercept = -1.50),
      MZ2 = list(intercept = -1.60, MZ1 = 2),
      MZ3 = list(intercept = -1.40, MZ1 = 2, MZ2 = 2),
      MZ4 = list(intercept = -3.40, MZ1 = 2.1, MZ2 = 2.1, MZ3 = 2.1),
      MA  = list(intercept = -2.50, MZ1 = 2, MZ2 = 2, MZ3 = 2, MZ4 = 2),
      MY  = list(intercept = -1.65, MZ1 = 0.5, MZ2 = 0.5, MZ3 = 0.5, MZ4 = 0.5, MA = -0.2)
    )
    coefsC <- list(
      MZ1 = list(intercept = -1.15),
      MZ2 = list(intercept = -1.35, MZ1 = 2.2),
      MZ3 = list(intercept = -1.20, MZ1 = 2.2, MZ2 = 2.2),
      MZ4 = list(intercept = -3.10, MZ1 = 2.2, MZ2 = 2.2, MZ3 = 2.2),
      MA  = list(intercept = -2.25, MZ1 = 1.8, MZ2 = 1.8, MZ3 = 1.8, MZ4 = 1.8),
      MY  = list(intercept = -1.30, MZ1 = 0.4, MZ2 = 0.4, MZ3 = 0.4, MZ4 = 0.4, MA = -0.65)
    )
    coefsD <- list(
      MZ1 = list(intercept = -1.15),
      MZ2 = list(intercept = -1.30, MZ1 = 2.2),
      MZ3 = list(intercept = -1.20, MZ1 = 2.2, MZ2 = 2.2),
      MZ4 = list(intercept = -3.10, MZ1 = 2.2, MZ2 = 2.2, MZ3 = 2.2),
      MA  = list(intercept = -2.25, MZ1 = 1.8, MZ2 = 1.8, MZ3 = 1.8, MZ4 = 1.8),
      MY  = list(intercept = -1.35, MZ1 = 0.4, MZ2 = 0.4, MZ3 = 0.4, MZ4 = 0.4, MA = -0.65)
    )
    coefsE <- list(
      MZ1 = list(intercept = -1.15),
      MZ2 = list(intercept = -1.30, MZ1 = 2.2),
      MZ3 = list(intercept = -1.20, MZ1 = 2.2, MZ2 = 2.2),
      MZ4 = list(intercept = -3.10, MZ1 = 2.2, MZ2 = 2.2, MZ3 = 2.2),
      MA  = list(intercept = -2.25, MZ1 = 1.8, MZ2 = 1.8, MZ3 = 1.8, MZ4 = 1.8),
      MY  = list(intercept = -1.35, MZ1 = 0.4, MZ2 = 0.4, MZ3 = 0.4, MZ4 = 0.4, MA = -0.65)
    )
 } else if (pos == "pos2") {
    coefsA <- list(
      MZ1 = list(intercept = -1.75),
      MZ2 = list(intercept = -1.80, MZ1 = 2),
      MZ3 = list(intercept = -1.65, MZ1 = 2, MZ2 = 2),
      MZ4 = list(intercept = -3.55, MZ1 = 2.1, MZ2 = 2.1, MZ3 = 2.1),
      MA = list(intercept = -2.75, MZ1 = 2, MZ2 = 2, MZ3 = 2, MZ4 = 2),
      MY = list(intercept = -2.00, MZ1 = 0.6, MZ2 = 0.6, MZ3 = 0.6, MZ4 = 0.6, MA = -0.2)
    )
    coefsB <- list(
      MZ1 = list(intercept = -1.50),
      MZ2 = list(intercept = -1.60, MZ1 = 2),
      MZ3 = list(intercept = -1.40, MZ1 = 2, MZ2 = 2),
      MZ4 = list(intercept = -3.40, MZ1 = 2.1, MZ2 = 2.1, MZ3 = 2.1),
      MA = list(intercept = -2.50, MZ1 = 2, MZ2 = 2, MZ3 = 2, MZ4 = 2),
      MY = list(intercept = -1.65, MZ1 = 0.5, MZ2 = 0.5, MZ3 = 0.5, MZ4 = 0.5, MA = -0.2)
    )
    coefsC <- list(
      MZ1 = list(intercept = -1.15),
      MZ2 = list(intercept = -1.35, MZ1 = 2.2),
      MZ3 = list(intercept = -1.20, MZ1 = 2.2, MZ2 = 2.2),
      MZ4 = list(intercept = -3.10, MZ1 = 2.2, MZ2 = 2.2, MZ3 = 2.2),
      MA = list(intercept = -2.15, MZ1 = 1.8, MZ2 = 1.8, MZ3 = 1.8, MZ4 = 1.8),
      MY = list(intercept = -1.30, MZ1 = 0.4, MZ2 = 0.4, MZ3 = 0.4, MZ4 = 0.4, MA = -0.65)
    )
    coefsD <- list(
      MZ1 = list(intercept = -1.15),
      MZ2 = list(intercept = -1.30, MZ1 = 2.2),
      MZ3 = list(intercept = -1.20, MZ1 = 2.2, MZ2 = 2.2),
      MZ4 = list(intercept = -3.10, MZ1 = 2.2, MZ2 = 2.2, MZ3 = 2.2),
      MA = list(intercept = -2.20, MZ1 = 1.8, MZ2 = 1.8, MZ3 = 1.8, MZ4 = 1.8),
      MY = list(intercept = -1.30, MZ1 = 0.4, MZ2 = 0.4, MZ3 = 0.4, MZ4 = 0.4, MA = -0.65)
    )
    coefsE <- list(
      MZ1 = list(intercept = -1.15),
      MZ2 = list(intercept = -1.30, MZ1 = 2.2),
      MZ3 = list(intercept = -1.20, MZ1 = 2.2, MZ2 = 2.2),
      MZ4 = list(intercept = -3.10, MZ1 = 2.2, MZ2 = 2.2, MZ3 = 2.2),
      MA = list(intercept = -2.20, MZ1 = 1.8, MZ2 = 1.8, MZ3 = 1.8, MZ4 = 1.8),
      MY = list(intercept = -1.30, MZ1 = 0.4, MZ2 = 0.4, MZ3 = 0.4, MZ4 = 0.4, MA = -0.65)
    )
  } else if (pos == "pos3") {
    coefsA <- list(
      MZ1 = list(intercept = -1.75),
      MZ2 = list(intercept = -1.80, MZ1 = 2),
      MZ3 = list(intercept = -1.65, MZ1 = 2, MZ2 = 2),
      MZ4 = list(intercept = -3.55, MZ1 = 2.1, MZ2 = 2.1, MZ3 = 2.1),
      MA = list(intercept = -2.75, MZ1 = 2, MZ2 = 2, MZ3 = 2, MZ4 = 2),
      MY = list(intercept = -2.00, MZ1 = 0.6, MZ2 = 0.6, MZ3 = 0.6, MZ4 = 0.6, MA = -0.2)
    )
    coefsB <- list(
      MZ1 = list(intercept = -1.50),
      MZ2 = list(intercept = -1.60, MZ1 = 2),
      MZ3 = list(intercept = -1.40, MZ1 = 2, MZ2 = 2),
      MZ4 = list(intercept = -3.40, MZ1 = 2.1, MZ2 = 2.1, MZ3 = 2.1),
      MA = list(intercept = -2.50, MZ1 = 2, MZ2 = 2, MZ3 = 2, MZ4 = 2),
      MY = list(intercept = -1.65, MZ1 = 0.5, MZ2 = 0.5, MZ3 = 0.5, MZ4 = 0.5, MA = -0.2)
    )
    coefsC <- list(
      MZ1 = list(intercept = -1.15),
      MZ2 = list(intercept = -1.30, MZ1 = 2.2),
      MZ3 = list(intercept = -1.15, MZ1 = 2.2, MZ2 = 2.2),
      MZ4 = list(intercept = -3.05, MZ1 = 2.2, MZ2 = 2.2, MZ3 = 2.2),
      MA = list(intercept = -2.15, MZ1 = 1.8, MZ2 = 1.8, MZ3 = 1.8, MZ4 = 1.8),
      MY = list(intercept = -1.25, MZ1 = 0.4, MZ2 = 0.4, MZ3 = 0.4, MZ4 = 0.4, MA = -0.65)
    )
    coefsD <- list(
      MZ1 = list(intercept = -1.15),
      MZ2 = list(intercept = -1.30, MZ1 = 2.2),
      MZ3 = list(intercept = -1.15, MZ1 = 2.2, MZ2 = 2.2),
      MZ4 = list(intercept = -3.10, MZ1 = 2.2, MZ2 = 2.2, MZ3 = 2.2),
      MA = list(intercept = -2.15, MZ1 = 1.8, MZ2 = 1.8, MZ3 = 1.8, MZ4 = 1.8),
      MY = list(intercept = -1.25, MZ1 = 0.4, MZ2 = 0.4, MZ3 = 0.4, MZ4 = 0.4, MA = -0.65)
    )
    coefsE <- list(
      MZ1 = list(intercept = -1.15),
      MZ2 = list(intercept = -1.30, MZ1 = 2.2),
      MZ3 = list(intercept = -1.15, MZ1 = 2.2, MZ2 = 2.2),
      MZ4 = list(intercept = -3.10, MZ1 = 2.2, MZ2 = 2.2, MZ3 = 2.2),
      MA = list(intercept = -2.15, MZ1 = 1.8, MZ2 = 1.8, MZ3 = 1.8, MZ4 = 1.8),
      MY = list(intercept = -1.25, MZ1 = 0.4, MZ2 = 0.4, MZ3 = 0.4, MZ4 = 0.4, MA = -0.65)
    )
  }
  
  # Generate each m-DAG using your existing generateDag_Data() function.
  mDagA <- generate_mDag_Data(simu_data_list, coefsA, "A")
  mDagB <- generate_mDag_Data(simu_data_list, coefsB, "B")
  mDagC <- generate_mDag_Data(simu_data_list, coefsC, "C")
  mDagD <- generate_mDag_Data(simu_data_list, coefsD, "D")
  mDagE <- generate_mDag_Data(simu_data_list, coefsE, "E")
  
  # Return a list containing all five m-DAG objects.
  # Create a suffix for names if pos is not "no_pos"
  suffix <- if (pos == "no_pos") "" else paste0("_", pos)
  result_list <- setNames(list(mDagA, mDagB, mDagC, mDagD, mDagE),
                          c(paste0("mDagA", suffix),
                            paste0("mDagB", suffix),
                            paste0("mDagC", suffix),
                            paste0("mDagD", suffix),
                            paste0("mDagE", suffix)))
  return(result_list)
}

# help function for missingness proportions
MissingProportion <- function(data_list) {
  
  Missing_enwast <- vector(length = length(data_list))
  Missing_W_mhtcm <- vector(length = length(data_list))
  Missing_W_mwtkg <- vector(length = length(data_list))
  Missing_W_parity <- vector(length = length(data_list))
  Missing_a <- vector(length = length(data_list))
  Missing_haz <- vector(length = length(data_list))
  Missing_a_haz <- vector(length = length(data_list))
  MissingAny <- vector(length = length(data_list))
  
  for (i in 1:(length(data_list))) {
    Missing_enwast[i] <- sum(is.na(data_list[[i]]$enwast))/nrow(data_list[[i]])
    Missing_W_mhtcm[i] <- sum(is.na(data_list[[i]]$W_mhtcm))/nrow(data_list[[i]])
    Missing_W_mwtkg[i] <- sum(is.na(data_list[[i]]$W_mwtkg))/nrow(data_list[[i]])
    Missing_W_parity[i] <- sum(is.na(data_list[[i]]$W_parity))/nrow(data_list[[i]])
    Missing_a[i] <- sum(is.na(data_list[[i]]$a))/nrow(data_list[[i]])
    Missing_haz[i] <- sum(is.na(data_list[[i]]$haz))/nrow(data_list[[i]])
    Missing_a_haz[i] <- sum(is.na(data_list[[i]]$haz)|is.na(data_list[[i]]$a))/nrow(data_list[[i]])
    MissingAny[i] <- sum(is.na(data_list[[i]]$enwast)|is.na(data_list[[i]]$W_mhtcm)|
                           is.na(data_list[[i]]$W_mwtkg)|is.na(data_list[[i]]$W_parity)|
                           is.na(data_list[[i]]$a)|is.na(data_list[[i]]$haz))/nrow(data_list[[i]])
  }
  mean_missing <- c(mean(Missing_enwast), mean(Missing_W_mhtcm), mean(Missing_W_mwtkg),
                    mean(Missing_W_parity), mean(Missing_a), mean(Missing_haz), mean(Missing_a_haz), mean(MissingAny))
  #mean_missing <- as.data.frame(t(mean_missing), row.names = deparse(substitute(data_list))) -> cut after first '_' for better naming
  names(mean_missing) <- c("enwast", "W_mhtcm", "W_mwtkg", "W_parity", "a", "haz", "aorHaz", "Any")
  return(mean_missing)
} 

create_missing_df <- function(list_object) {
  # Apply MissingProportion to each element of the list
  missing_list <- lapply(list_object, MissingProportion)
  df <- do.call(rbind, missing_list)
  df <- as.data.frame(df)
  # include the names of each list element as a column in the data frame.
  # df <- cbind(list_name = names(list_object), df)
  
  return(df)
}


