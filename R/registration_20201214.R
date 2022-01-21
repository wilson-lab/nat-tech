library(nat)
library(natverse)
library(nat.templatebrains)
library(nat.jrcbrains)
library(neuprintr)
#download_saalfeldlab_registrations()
register_saalfeldlab_registrations()

#plot neuprint neuron FC1 to compare to confocal


#import traced neuron (from light level image)
neuronTracing1 <- read.neuron('/Users/WilsonLab/Desktop/SNTsnapshots/R20E08AD_R64A11_15-000.swc', class="neuron")

nopen3d()
plot3d(neuronTracing1,lwd=3,col='black',WithNodes=FALSE)


#import registration file and apply registration to the traced neuron file
#reg <- as.cmtkreg('/Users/WilsonLab/Desktop/Registration/Registration_text_files/registration_15_3')
#registeredNeuron1 <- xform(neuronTracing1, reg)
#registeredNeuron2 <- xform(neuronTracing2, reg)

#Plot neuron along with registered brain template 
#nopen3d()
#plot3d(registeredNeuron1,lwd=3,col='black',WithNodes=FALSE)
#plot3d(registeredNeuron2,lwd=3,col='black',WithNodes=FALSE)
plot3d(JRC2018U)

#save the registered neuron to a .swc file 
#write.neuron(registeredNeuron,'data/registeredNeuron_JFRC2013',format='swc')

#transform neuron from JFC2018U template brain space to FAFB14 template brain space 
neuron_FAFB=xform_brain(registeredNeuron,sample=JRC2018U,reference=FAFB14)

#plot the neuron with its new template brain
nopen3d()
plot3d(neuron_FAFB,lwd=3,col='black',WithNodes=FALSE)
plot3d(FAFB14)
