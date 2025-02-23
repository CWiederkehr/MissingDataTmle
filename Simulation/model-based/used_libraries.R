
# Define a vector with the names of the required packages
packages <- c("tmle", "SuperLearner", "mice", "Amelia", "ggplot2", 
              "gridExtra", "scales", "foreach", "doParallel", 
              "copula", "ranger", "earth", "arm", "glmnet")

# Check which packages are not installed.
# The installed.packages() function returns a matrix where the column "Package" contains the names of installed packages.
new.packages <- packages[!(packages %in% installed.packages()[, "Package"])]

# If there are packages that are not installed, install them.
if(length(new.packages)) {
  install.packages(new.packages, dependencies = TRUE)
}

# Load each package using a loop. The argument 'character.only = TRUE' ensures that the package name is taken as a string.
for(pkg in packages) {
  library(pkg, character.only = TRUE)
}
