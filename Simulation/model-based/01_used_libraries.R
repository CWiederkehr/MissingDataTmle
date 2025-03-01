
packages <- c("tmle", "SuperLearner", "mice", "Amelia", "ggplot2", 
              "gridExtra", "scales", "foreach", "doParallel", 
              "copula", "ranger", "earth", "arm", "glmnet")

new.packages <- packages[!(packages %in% installed.packages()[, "Package"])]

if(length(new.packages)) {
  install.packages(new.packages, dependencies = TRUE)
}

for(pkg in packages) {
  library(pkg, character.only = TRUE)
}
