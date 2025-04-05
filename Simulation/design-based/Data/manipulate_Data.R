
wash_b <- read.csv("../MissingDataTmle/Simulation/design-based/Data/original_wash_b.csv")

# Convert characters to factors
wash_b$sex <- as.factor(wash_b$sex)
wash_b$hfoodsec <- as.factor(wash_b$hfoodsec)
wash_b$hfoodsec <- factor(wash_b$hfoodsec, levels = c("Food Insecure","Mildly Food Insecure","Food Secure"))

# Remove redundant indicator variables
wash_b <- wash_b[, !(names(wash_b) %in% c("delta_W_mage", "delta_W_mhtcm", "delta_W_mwtkg", 
                                          "delta_W_mbmi", "delta_enwast", "delta_impsan", 
                                          "delta_W_feducyrs", "delta_W_parity", "impsan"))]

save(wash_b, file = "../MissingDataTmle/Simulation/design-based/Data/wash_b.RData")
