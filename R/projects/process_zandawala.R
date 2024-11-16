# Load libraries
source("R/startup/packages.R")
source("R/startup/functions.R")

# Save mip files
neuronbridger::nrrd_to_mip(fiji.path = neuronbridger:::fiji())
