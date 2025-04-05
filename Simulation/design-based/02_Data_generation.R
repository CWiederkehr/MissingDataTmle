source("../MissingDataTmle/Simulation/design-based/01_used_libraries.R")
source("../MissingDataTmle/Simulation/design-based/UnderHAL_Functions.R")
load("../MissingDataTmle/Simulation/design-based/Data/wash_b.RData")

set.seed(5)

# Sample 2000 rows with replacement
wash_b_2000 <- wash_b[sample(nrow(wash_b), 2000, replace = TRUE), ]
str(wash_b_2000)


# 2000 sample size -> 125 runs for 'true ATE'
ss_estimate2000 <-  ss_estimator(wash_b_2000, colnames(wash_b_2000)[-c(9,10)], "a", "haz")

ss_estimate2000$psi  # In original paper: 0.0507 (vs here  0.1101056)
ss_estimate2000$rv  # here: 0.4263114
ss_estimate2000$g_fit$lambda_star
ss_estimate2000$Q_fit$lambda_star
sum(ss_estimate2000$g_fit$coefs !=0)  # In original paper: 124 (vs here 131)
sum(ss_estimate2000$Q_fit$coefs !=0)  # In original paper: 1496 (vs here 2210)

# convert 'factor' variables 
wash_b_cov_2000 <- wash_b_2000 %>% select(all_of(names(wash_b_2000)[-c(9,10)])) %>% 
  mutate_if(sapply(., is.factor), as.numeric)

# generate prob_predictions
A_pred_2000 <- predict(ss_estimate2000$g_fit, new_data = wash_b_cov_2000)

save(ss_estimate2000, A_pred_2000, wash_b_2000, file="../MissingDataTmle/Simulation/design-based/Results/ss_estimate2000.RData")



## n=2000 
simu_2000 <- run_simulation(df = wash_b_2000, w_cols = colnames(wash_b_2000)[-c(9, 10)], aname = "a",
                            yname = "haz", g_fit = ss_estimate2000$g_fit, Q_fit = ss_estimate2000$Q_fit, rv = ss_estimate2000$rv,
                            weights = NULL, Nlam = 100, Sim = 1000, Size = 1)
save(simu_2000, file = "simu_2000_data.RData")

simu_2000 <- calculate_vio_metrics2(simu_2000, w_cols = colnames(simu_2000$simu_data_list2000[[1]])[-c(9, 10)], aname = "a", start_index = 1, end_index = 1000, 
                                    degree_crit = 10, total_cores = 20)
save(simu_2000, file = "simu_2000_data.RData")

mean(simu_2000$violatio) # already 1%

new_A_pred_2000 <- calculate_custom_weights(pred_A = A_pred_2000, real_A = wash_b_2000$a, weight = 10)

simu_2000_pos1 <- run_simulation(df = wash_b_2000, w_cols = colnames(wash_b_2000)[-c(9, 10)], aname = "a",
                                 yname = "haz", g_fit = ss_estimate2000$g_fit, Q_fit = ss_estimate2000$Q_fit, rv = ss_estimate2000$rv,
                                 weights = new_A_pred_2000, Nlam = 100, Sim = 1000, Size = 1, vio = TRUE, total_cores = 20)


simu_2000_pos1$violation
mean(unlist(lapply(simu_2000_pos1$propScore_preds2000, count_violations))) # 19.2%
save(simu_2000_pos1, file="simu_2000_pos1.RData")



new_A_pred_2000 <- calculate_custom_weights(pred_A = A_pred_2000, real_A = wash_b_2000$a, weight = 5)

simu_2000_pos2 <- run_simulation(df = wash_b_2000, w_cols = colnames(wash_b_2000)[-c(9, 10)], aname = "a",
                                 yname = "haz", g_fit = ss_estimate2000$g_fit, Q_fit = ss_estimate2000$Q_fit, rv = ss_estimate2000$rv,
                                 weights = new_A_pred_2000, Nlam = 100, Sim = 1000, Size = 1, vio = TRUE, total_cores = 20)


simu_2000_pos2$violation
mean(unlist(lapply(simu_2000_pos2$propScore_preds2000, count_violations))) # 12.6%
save(simu_2000_pos2, file="simu_2000_pos2.RData")


new_A_pred_2000 <- calculate_custom_weights(pred_A = A_pred_2000, real_A = wash_b_2000$a, weight = 2)

simu_2000_pos3 <- run_simulation(df = wash_b_2000, w_cols = colnames(wash_b_2000)[-c(9, 10)], aname = "a",
                                 yname = "haz", g_fit = ss_estimate2000$g_fit, Q_fit = ss_estimate2000$Q_fit, rv = ss_estimate2000$rv,
                                 weights = new_A_pred_2000, Nlam = 100, Sim = 1000, Size = 1, vio = TRUE, total_cores = 20)


simu_2000_pos3$violation
mean(unlist(lapply(simu_2000_pos3$propScore_preds2000, count_violations))) # 5.5%
save(simu_2000_pos3, file="simu_2000_pos3.RData")