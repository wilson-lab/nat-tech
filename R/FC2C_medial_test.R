source("R/startup/packages.R")

#get all FC1 neurons
FC2C.info <- neuprint_search("type:FC2C.*")
#FC1_split <- split(FC1.info, FC1.info$type)

#get the neuron skeletons
FC2C_skel = neuprint_read_neurons(FC2C.info$bodyid)


#FC1A
FC2C.jrc2018f = xform_brain(FC2C_skel*8/1000, reference="JRC2018F", sample="JRCFIB2018F")
FC2C.jrc2018u = xform_brain(FC2C.jrc2018f, reference='JRC2018U', sample='JRC2018F')

nopen3d()
#plot FC1 neurons with different colors
plot3d(FC2C.jrc2018u,lwd=3,col='black',WithNodes=FALSE,soma=FALSE)
plot3d(JRC2018U)

#save the FC2C to a .nrrd file 
write.nrrd(as.im3d(FC2C.jrc2018u),'data/FC2C_JRC2018.nrrd')
