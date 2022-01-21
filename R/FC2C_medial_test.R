library(nat)
library(natverse)
library(neuprintr)
library(nat.templatebrains)
library(nat.jrcbrains)
register_saalfeldlab_registrations()


#establish neuprint database
conn = neuprint_login(server= "https://neuprint.janelia.org/",
                      token= "eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJlbWFpbCI6ImtlbGxvMjJlQG10aG9seW9rZS5lZHUiLCJsZXZlbCI6Im5vYXV0aCIsImltYWdlLXVybCI6Imh0dHBzOi8vbGgzLmdvb2dsZXVzZXJjb250ZW50LmNvbS9hL0FBVFhBSnoxcWZuYXZJTHRHd1FyWGY2OEZGYVdTcFdhU09IMTIwTlRsSFU9czk2LWM_c3o9NTA_c3o9NTAiLCJleHAiOjE4MTA1MTA0MDd9.7IWz4boFJ2iDUmr0zNIjFXGNhSa1KWESateQovBExEE")
options(neuprint_dataset="hemibrain:v1.2.1")

#get all FC1 neurons
FC1.info <- neuprint_search("type:FC2C.*")
#FC1_split <- split(FC1.info, FC1.info$type)

#get the neuron skeletons
FC2C_skel = neuprint_read_neurons(FC1.info$bodyid)


#FC1A
FC2C.jrc2018f = xform_brain(FC2C_skel*8/1000, reference="JRC2018F", sample="JRCFIB2018F")
FC2C.jrc2018u = xform_brain(FC2C.jrc2018f, reference='JRC2018U', sample='JRC2018F')

#nopen3d()
#plot FC1 neurons with different colors
#plot3d(FC2C.jrc2018u,lwd=3,col='black',WithNodes=FALSE,soma=FALSE)

#save the FC2C to a .nrrd file 
write.nrrd(FC2C.jrc2018u,'Desktop/FC2C_JRC2018')
