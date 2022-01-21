library(nat)
library(natverse)
library(neuprintr)
library(nat.templatebrains)
library(nat.jrcbrains)
register_saalfeldlab_registrations()

#ask user what neurons they want to plot 
neuron_type <- readline(prompt="Which FC1 subytpe do you want to plot? ")

#establish neuprint database
conn = neuprint_login(server= "https://neuprint.janelia.org/",
                      token= "eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJlbWFpbCI6ImtlbGxvMjJlQG10aG9seW9rZS5lZHUiLCJsZXZlbCI6Im5vYXV0aCIsImltYWdlLXVybCI6Imh0dHBzOi8vbGgzLmdvb2dsZXVzZXJjb250ZW50LmNvbS9hL0FBVFhBSnoxcWZuYXZJTHRHd1FyWGY2OEZGYVdTcFdhU09IMTIwTlRsSFU9czk2LWM_c3o9NTA_c3o9NTAiLCJleHAiOjE4MTA1MTA0MDd9.7IWz4boFJ2iDUmr0zNIjFXGNhSa1KWESateQovBExEE")
options(neuprint_dataset="hemibrain:v1.2.1")

#neuron we think is FC1
#import traced neuron (from light level image)
neuronTracing1 <- read.neuron('/Users/WilsonLab/Desktop/SNTsnapshots/R20E08AD_R64A11_15-000.swc', class="neuron")
neuronTracing2 <- read.neuron('/Users/WilsonLab/Desktop/SNTsnapshots/R20E08AD_R64A11_15-001.swc', class="neuron")
nopen3d()
plot3d(neuronTracing1,lwd=3,col='black',WithNodes=FALSE)
plot3d(neuronTracing2,lwd=3,col='black',WithNodes=FALSE)

#import registration file and apply registration to the traced neuron file
reg <- as.cmtkreg('/Users/WilsonLab/Desktop/Registration/Registration_text_files/registration_15_3')
registeredNeuron1 <- xform(neuronTracing1, reg)
registeredNeuron2 <- xform(neuronTracing2, reg)

#get all FC1 neurons
FC1.info <- neuprint_search("type:FC1.*")
FC1_split <- split(FC1.info, FC1.info$type)

#get the neuron skeletons
FC1_A = neuprint_read_neurons((FC1_split$FC1A)$bodyid)
FC1_B = neuprint_read_neurons((FC1_split$FC1B)$bodyid)
FC1_C = neuprint_read_neurons((FC1_split$FC1D)$bodyid)
FC1_D = neuprint_read_neurons((FC1_split$FC1E)$bodyid)
FC1_E = neuprint_read_neurons((FC1_split$FC1E)$bodyid)
FC1_F = neuprint_read_neurons((FC1_split$FC1F)$bodyid)

#FC1A
FC1A.jrc2018f = xform_brain(FC1_A*8/1000, reference="JRC2018F", sample="JRCFIB2018F")
FC1A.jrc2018u = xform_brain(FC1A.jrc2018f, reference='JRC2018U', sample='JRC2018F')

#FC1B
FC1B.jrc2018f = xform_brain(FC1_B*8/1000, reference="JRC2018F", sample="JRCFIB2018F")
FC1B.jrc2018u = xform_brain(FC1B.jrc2018f, reference='JRC2018U', sample='JRC2018F')

#FC1C
FC1C.jrc2018f = xform_brain(FC1_C*8/1000, reference="JRC2018F", sample="JRCFIB2018F")
FC1C.jrc2018u = xform_brain(FC1C.jrc2018f, reference='JRC2018U', sample='JRC2018F')

#FC1D
FC1D.jrc2018f = xform_brain(FC1_D*8/1000, reference="JRC2018F", sample="JRCFIB2018F")
FC1D.jrc2018u = xform_brain(FC1D.jrc2018f, reference='JRC2018U', sample='JRC2018F')

#FC1E
FC1E.jrc2018f = xform_brain(FC1_E*8/1000, reference="JRC2018F", sample="JRCFIB2018F")
FC1E.jrc2018u = xform_brain(FC1E.jrc2018f, reference='JRC2018U', sample='JRC2018F')

#FC1F
FC1F.jrc2018f = xform_brain(FC1_F*8/1000, reference="JRC2018F", sample="JRCFIB2018F")
FC1F.jrc2018u = xform_brain(FC1F.jrc2018f, reference='JRC2018U', sample='JRC2018F')

nopen3d()
#plot FC1 neurons with different colors
plot3d(FC1A.jrc2018u,lwd=3,col='red',WithNodes=FALSE,soma=FALSE)
plot3d(FC1B.jrc2018u,lwd=3,col='green',WithNodes=FALSE,soma=FALSE)
plot3d(FC1C.jrc2018u,lwd=3,col='yellow',WithNodes=FALSE,soma=FALSE)
plot3d(FC1D.jrc2018u,lwd=3,col='orange',WithNodes=FALSE,soma=FALSE)
plot3d(FC1E.jrc2018u,lwd=3,col='purple',WithNodes=FALSE,soma=FALSE)
plot3d(FC1F.jrc2018u,lwd=3,col='pink',WithNodes=FALSE,soma=FALSE)

#plot trace in black
plot3d(registeredNeuron1,lwd=3,col='black',WithNodes=FALSE)
plot3d(registeredNeuron2,lwd=3,col='black',WithNodes=FALSE)
plot3d(JRC2018U)