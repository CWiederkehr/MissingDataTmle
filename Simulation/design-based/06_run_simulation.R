source("../MisingDataTmle/Simulation/design-based/01_used_libraries.R")
source("../MisingDataTmle/Simulation/design-based/02_DGP_HAL_Functions.R")
source("../MisingDataTmle/Simulation/design-based/03_missing_Functions.R")
source("../MisingDataTmle/Simulation/design-based/04_analysis_Functions.R")
source("../MisingDataTmle/Simulation/design-based/05_plot_Functions.R")
load("../MisingDataTmle/Simulation/design-based/Data/wash_b.RData")  

set.seed(5)
Sim <- 1000

#### Fit underHAL for DGP ####
# 2000 sample size -> 125 runs for 'true ATE'
ss_estimate2000 <-  ss_estimator(wash_b, colnames(wash_b)[-c(9,10)], "a", "haz")

ss_estimate2000$psi  # simulated ATE  0.1101056)
ss_estimate2000$rv  # simulated residuals: 0.4263114
sum(ss_estimate2000$g_fit$coefs !=0)  # basis inicators for exposure-model:
sum(ss_estimate2000$Q_fit$coefs !=0)  # basis indicators for Outcome-model

# convert 'factor' variables 
wash_b_cov_2000 <- wash_b %>% select(all_of(names(wash_b)[-c(9,10)])) %>% 
  mutate_if(sapply(., is.factor), as.numeric)

# generate prob_predictions
A_pred_2000 <- predict(ss_estimate2000$g_fit, new_data = wash_b_cov_2000)

save(ss_estimate2000, A_pred_2000, wash_b, file="../MisingDataTmle/Simulation/design-based/Results/ss_estimate2000.RData")

load("../MisingDataTmle/Simulation/design-based/Preliminary_Results/ss_estimate2000.RData")

#### DGP ####
## 1% Positvity
sim_data_list <- run_DGP_via_HAL(df = wash_b, w_cols = colnames(wash_b)[-c(9, 10)], aname = "a",
                                 yname = "haz", g_fit = ss_estimate2000$g_fit, Q_fit = ss_estimate2000$Q_fit, rv = ss_estimate2000$rv,
                                 weights = NULL, Nlam = 100, Sim = Sim, Size = 1, vio = TRUE, total_cores = 20)

mean(sim_data_list$violation) # already 1%
save(sim_data_list, file = "../MisingDataTmle/Simulation/design-based/Results/sim_data_list.RData")


## 5.5% Positivity
new_A_weights2 <- calculate_custom_weights(pred_A = A_pred_2000, real_A = wash_b$a, weight = 2)

sim_data_list_pos1 <- run_DGP_via_HAL(df = wash_b, w_cols = colnames(wash_b)[-c(9, 10)], aname = "a",
                                      yname = "haz", g_fit = ss_estimate2000$g_fit, Q_fit = ss_estimate2000$Q_fit, rv = ss_estimate2000$rv,
                                      weights = new_A_weights2, Nlam = 100, Sim = Sim, Size = 1, vio = TRUE, total_cores = 20)

mean(sim_data_list_pos1$violation2000) # 5.5%
save(sim_data_list_pos1, file="../MisingDataTmle/Simulation/design-based/Results/sim_data_list_pos1.RData")


## 12.5% Positivity
new_A_weights5 <- calculate_custom_weights(pred_A = A_pred_2000, real_A = wash_b$a, weight = 5)

sim_data_list_pos2 <- run_DGP_via_HAL(df = wash_b, w_cols = colnames(wash_b)[-c(9, 10)], aname = "a",
                                      yname = "haz", g_fit = ss_estimate2000$g_fit, Q_fit = ss_estimate2000$Q_fit, rv = ss_estimate2000$rv,
                                      weights = new_A_weights5, Nlam = 100, Sim = Sim, Size = 1, vio = TRUE, total_cores = 20)


mean(sim_data_list_pos2$violation2000) # 12.6%
save(sim_data_list_pos2, file="../MisingDataTmle/Simulation/design-based/Results/sim_data_list_pos2.RData")


## 20% Positivity
new_A_weights10 <- calculate_custom_weights(pred_A = A_pred_2000, real_A = wash_b$a, weight = 10)

sim_data_list_pos3 <- run_DGP_via_HAL(df = wash_b, w_cols = colnames(wash_b)[-c(9, 10)], aname = "a",
                                      yname = "haz", g_fit = ss_estimate2000$g_fit, Q_fit = ss_estimate2000$Q_fit, rv = ss_estimate2000$rv,
                                      weights = new_A_weights10, Nlam = 100, Sim = Sim, Size = 1, vio = TRUE, total_cores = 20)


mean(sim_data_list_pos3$violation2000) # 19.3%
save(sim_data_list_pos3, file="../MisingDataTmle/Simulation/design-based/Results/sim_data_list_pos3.RData")



#### Impose Missingness ####
load("../MisingDataTmle/Simulation/design-based/Preliminary_Results/all_sim_data.RData")

#sim_data_list$simu_data_list2000 <- sim_data_list$simu_data_list2000[1:100] # reduce computational burden
#sim_data_list_pos1$simu_data_list2000 <- sim_data_list_pos1$simu_data_list2000[1:100]
#sim_data_list_pos2$simu_data_list2000 <- sim_data_list_pos2$simu_data_list2000[1:100]
#sim_data_list_pos3$simu_data_list2000 <- sim_data_list_pos3$simu_data_list2000[1:100]

mDags <- generate_All_mDags(simu_data_list = sim_data_list$simu_data_list2000, pos = "no_pos")
mDags_pos1 <- generate_All_mDags(simu_data_list = sim_data_list_pos1$simu_data_list2000, pos = "pos1")
mDags_pos2 <- generate_All_mDags(simu_data_list = sim_data_list_pos2$simu_data_list2000, pos = "pos2")
mDags_pos3 <- generate_All_mDags(simu_data_list = sim_data_list_pos3$simu_data_list2000, pos = "pos3")

create_missing_df(mDags); create_missing_df(mDags_pos1); create_missing_df(mDags_pos2); create_missing_df(mDags_pos3)

# big list mDags
all_missingness_data <- c(mDags, mDags_pos1, mDags_pos2, mDags_pos3)


all_results <- list()

for (orig_key in names(all_missingness_data)) {
  
  # Call the apply_all_methods function for this element, passing the proper DGP
  method_res <- apply_all_methods(
    data_list = all_missingness_data[[orig_key]], 
    cores = 4, 
    m = 100, 
    maxit = 10
  )
  
  new_key <- paste0(orig_key, "_res")
  
  # Store the result under the new key
  all_results[[new_key]] <- method_res
}

save(all_results, file="../MisingDataTmle/Simulation/design-based/Results/all_results.RData")
load("../MisingDataTmle/Simulation/design-based/Preliminary_Results/all_results.RData")

#### Result-Plots ####

# Make sure each element from all_results is assigned to the global environment.
for(nm in names(all_results)) {
  assign(nm, all_results[[nm]])
}

# Create a list of names and use these to call all_measures_TMLE, preserving the object names.
all_obj_names <- names(all_results)
all_measures_list <- lapply(all_obj_names, function(nm) {
  expr <- substitute(all_measures_TMLE(x), list(x = as.name(nm)))
  eval(expr, envir = .GlobalEnv)
})
names(all_measures_list) <- all_obj_names

big_table <- do.call(rbind, all_measures_list)

big_table$Scenario <- factor(big_table$Scenario, levels = c("Scenario1", "Scenario2", "Scenario3", "Scenario4"))
big_table$Method <- factor(big_table$Method, levels = c("CC", "Ext", "Ext MCMI", "MI PMM", "MI Int", "MI CART",
                                                        "MI RF", "MI RF MI", "MI Amelia"))
big_table$mDag <- factor(big_table$mDag, levels = c("A", "B", "C", "D", "E"))
rownames(big_table) <- sub(".*\\.", "", rownames(big_table))

# Inspect the resulting big_table
str(big_table)
head(big_table)

save(big_table, file="../MisingDataTmle/Simulation/design-based/Results/big_table_results.RData")

#### Result Plots ####

# Bias
plot_Bias <- OverviewPlotBias(big_table)

ggsave("../MisingDataTmle/Simulation/design-based/Results/Plot_Bias.png",
       plot = plot_Bias,
       width = 8.88, height = 4.80)

# Coverage
plot_Coverage <- OverviewPlotCoverage(big_table)

ggsave("../MisingDataTmle/Simulation/design-based/Results/Plot_Coverage.png",
       plot = plot_Coverage,
       width = 8.88, height = 4.80)

# RMSE
plot_RMSE <- OverviewPlotRMSE(big_table)
ggsave("../MisingDataTmle/Simulation/design-based/Results/Plot_RMSE.png",
       plot = plot_RMSE,
       width = 8.88, height = 4.80)



