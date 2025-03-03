

source("01_used_libraries.R")
source("UnderHalSimulationFcs.R")

load("~/GitHub/MissingDataTmle/Simulation/design-based/Data/wash_b.RData")

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

save(ss_estimate2000, A_pred_2000, wash_b_2000, file="ss_estimate2000.RData")