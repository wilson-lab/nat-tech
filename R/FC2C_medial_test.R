source("R/startup/packages.R")
source("R/startup/functions.R")

#hemibrain_to_JRC2018U_nrrd("FC2C")

# challenge: 
# 1: make a new function to read a .nrrd file (read.nrrd) file and plot with the right hemibrain neuron and template brain (let's just assume read.nrrd)
# 2: read a saved swc file (read.neuron) and plot with the right hemibrain neuron and template brain (let's just assume read.nrrd)
neuron_to_hemibrain("/Users/WilsonLab/Desktop/SNTsnapshots/fc2c_test-000.swc","FC2C")
