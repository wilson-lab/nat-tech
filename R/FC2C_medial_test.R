source("R/startup/packages.R")
source("R/startup/functions.R")

#hemibrain_to_JRC2018U_nrrd("FC2C")

# 1: make a new function to read a .nrrd file (read.nrrd) file and plot with the right hemibrain neuron and template brain (let's just assume read.nrrd)
nrrd_to_hemibrain("/Users/WilsonLab/Desktop/Registration/Reformatted/FC2Cmedial/JRC2018_FC2C_R50B05AD_R64A11GDBD_02_warp_m0g80c8e1e-1x26r4.nrrd","FC2C")

# 2: read a saved swc file (read.neuron) and plot with the right hemibrain neuron and template brain (let's just assume read.nrrd)
neuron_to_hemibrain("/Users/WilsonLab/Desktop/SNTsnapshots/fc2c_test-000.swc","FC2C")

#test function to save a hemibrain neuron to any template space
hemibrain_to_nrrd("FC2C","JRC2018F")

#turn flywire neuron into a nrrd file
flywireid_to_nrrd("720575940619059120","idk")

#get all FB neurons
hemibrain_to_nrrd("FB8A")
hemibrain_to_nrrd("FB8B")
hemibrain_to_nrrd("FB8C")
hemibrain_to_nrrd("FB8D")
hemibrain_to_nrrd("FB8E")
hemibrain_to_nrrd("FB8F")
hemibrain_to_nrrd("FB8G")
hemibrain_to_nrrd("FB8H")