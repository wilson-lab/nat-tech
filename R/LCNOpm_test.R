# write flywire neuron to nrrd, sample = "FAFB14", xyzmatrix
flywire_to_nrrd <- function(flywire_id, cell_type, ref="JRC2018U", savefolder = "data", plot3D =TRUE){
  
  #read in flywire ID
  flywire_neuron = read_cloudvolume_meshes(flywire_id)
  
  #sanity check - check to see you successfully read in flywire neuron
  if(plot3D){
    plot3d(flywire_neuron)
    plot3d(FAFB14)
  }
  
  
  #transform neuron into the correct template space
  flywire.reg = xform_brain(flywire_neuron*8/1000, reference=ref, sample="FAFB14")
  
  #sanity check - check to see you successfully transformed flywire neuron
  if(plot3D){
    nopen3d()
    plot3d(flywire.reg,lwd=3,col='black',WithNodes=FALSE,soma=FALSE)
    plot3d(JRC2018U)
  }
  
  #make im3d of flywire neuron
  x <- get(ref)
  points = xyzmatrix(flywire.reg)
  y = as.im3d(points,x)
  
  #save the flywire neuron as a .nrrd file
  dir.create(savefolder)
  write.nrrd(y, file.path(savefolder, sprintf("%s_%s.nrrd",cell_type,ref)))
}


#make im3d of flywire neuron
x <- get(ref)
points = xyzmatrix(flywire_neuron)
y = as.im3d(points,x)

#sanity check - check to see you successfully transformed flywire neuron
if(plot3D){
  nopen3d()
  plot3d(flywire.reg,lwd=3,col='red',WithNodes=FALSE,soma=FALSE)
  plot3d(x)
}