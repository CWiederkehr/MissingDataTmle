
source("GitHub/MissingDataTmle/Simulation/model-based/used_libraries.R")
source("GitHub/MissingDataTmle/Simulation/model-based/Data_and_Missingness.R")
source("GitHub/MissingDataTmle/Simulation/model-based/general_Functions.R")
source("GitHub/MissingDataTmle/Simulation/model-based/analysis_Functions.R")

# Simulate datasets for DGP1, DGP2, DGP3, DGP4, and DGP5
DGP1_data <- list(
  sce1_DGP1 = simulateDatasets("sce1", "DGP1"),
  sce2_DGP1 = simulateDatasets("sce2", "DGP1"),
  sce3_DGP1 = simulateDatasets("sce3", "DGP1")
)

DGP2_data <- list(
  sce1_DGP2 = simulateDatasets("sce1", "DGP2"),
  sce2_DGP2 = simulateDatasets("sce2", "DGP2"),
  sce3_DGP2 = simulateDatasets("sce3", "DGP2")
)

DGP3_data <- list(
  sce1_DGP3 = simulateDatasets("sce1", "DGP3"),
  sce2_DGP3 = simulateDatasets("sce2", "DGP3"),
  sce3_DGP3 = simulateDatasets("sce3", "DGP3")
)

DGP4_data <- list(
  sce1_DGP4 = simulateDatasets("sce1", "DGP4"),
  sce2_DGP4 = simulateDatasets("sce2", "DGP4"),
  sce3_DGP4 = simulateDatasets("sce3", "DGP4")
)

DGP5_data <- list(
  sce1_DGP5 = simulateDatasets("sce1", "DGP5"),
  sce2_DGP5 = simulateDatasets("sce2", "DGP5"),
  sce3_DGP5 = simulateDatasets("sce3", "DGP5")
)