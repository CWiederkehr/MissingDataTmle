
# Define a vector of required packages (removing duplicates)
packages <- c("data.table", "dplyr", "foreach", "stringr", 
              "glmnet", "here", "tidyverse", "forcats", 
              "origami", "hal9001", "tictoc", "doParallel")

# Identify packages that are not yet installed
new.packages <- packages[!(packages %in% installed.packages()[, "Package"])]

# Install any packages that are missing, along with their dependencies
if(length(new.packages)) {
  install.packages(new.packages, dependencies = TRUE)
}

# Load each package into the current session
for(pkg in packages) {
  library(pkg, character.only = TRUE)
}