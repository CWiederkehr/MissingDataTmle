
# Define a vector of required packages (removing duplicates)
packages <- c("data.table", "dplyr", "foreach", "stringr", 
              "glmnet", "here", "tidyverse", "forcats", 
              "origami", "hal9001", "tictoc", "doParallel",
              "tmle", "SuperLearner", "mice", "Amelia", 
              "ranger", "earth", "arm", "rpart")

packages <- unique(packages)
new.packages <- packages[!(packages %in% installed.packages()[, "Package"])]

if(length(new.packages)) {
  install.packages(new.packages, dependencies = TRUE)
}

for(pkg in packages) {
  library(pkg, character.only = TRUE)
}
