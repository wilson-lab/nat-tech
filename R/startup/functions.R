#this doesn't work, read in a .nrrd file of a confocal image and plot together with a hemibrain neuron
nrrd_to_hemibrain <- function (file, cell_type){
  #read in nrrd file (i'm just assuming that this doesn't have to be transformed)
  i <- read.nrrd(file)
  df <- data.frame(x=attr(i,'x'),y=attr(i,'y'),z=attr(i,'z'))
  
  #get hemibrain neuron and skeletonize
  hbn <- neuprint_search(sprintf("type:%s",cell_type))
  hbn_skel <- neuprint_read_neurons(hbn.$bodyid)
  
  #transform from hemibrain to template space (JRC2018U)
  hbn.jrc2018f = xform_brain(hbn_skel*8/1000, reference="JRC2018F", sample="JFRCIB2018F")
  hbn.jrc2018u = xform_brain(hbn_skel*8/1000, reference="JRC2018U", sample="JRC2018F")
  
  #plot nrrd and hemibrain neuron
  nopen3d()
  plot3d(neuron())
  plot3d(hbn.jrc2018u, lwd=3, col='black',WithNodes=FALSE,some=FALSE)
  plot3d(JRC2018U)
  
}

#take a traced neuron from a confocal stack and plot together with a neuron from the hemibrain
neuron_to_hemibrain <- function(trace, cell_type){
  # 2: read a saved swc file (read.neuron) and plot with the right hemibrain neuron and template brain (let's just assume read.nrrd)
  #read in swc file
  neuronTracing <- read.neuron(sprintf('%s', trace), class="neuron")
  
  #get hemibrain neuron
  hbn.info <- neuprint_search(sprintf("type:%s.*",cell_type))
  
  #get the neuron skeletons
  hbn_skel = neuprint_read_neurons(hbn.info$bodyid)
  
  #transform hemibrain neuron into correct template space
  hbn.jrc2018f = xform_brain(hbn_skel*8/1000, reference="JRC2018F", sample="JRCFIB2018F")
  hbn.jrc2018u = xform_brain(hbn.jrc2018f, reference='JRC2018U', sample='JRC2018F')
  
  #plot .swc file and hemibrain neuron in correct template space
  nopen3d()
  plot3d(neuronTracing,lwd=3,col='green',WithNodes=FALSE)
  plot3d(hbn.jrc2018u, lwd=3, col='black',WithNodes=FALSE,some=FALSE)
  plot3d(JRC2018U)
}

#transform a hemibrain neuron into a .nrrd file
hemibrain_to_nrrd <- function(cell_type, ref="JRC2018U", savefolder = "data", plot3D = TRUE){
  #get hemibrain neuron
  hbn.info <- neuprint_search(sprintf("type:%s.*",cell_type))
  
  #get the neuron skeletons
  hbn_skel = neuprint_read_neurons(hbn.info$bodyid)
  class(hbn_skel)
  
  #transform hemibrain neuron to template space
  hbn.reg = xform_brain(hbn_skel*8/1000, reference=ref, sample="JRCFIB2018F")
  
  # make im3d
  x <- get(ref)
  points=xyzmatrix(hbn.reg)
  I=as.im3d(points,x)
  
  if (plot3D){
    nopen3d()
    #plot FC1 neurons with different colors
    plot3d(hbn.reg,lwd=3,col='black',WithNodes=FALSE,soma=FALSE)
    plot3d(x)
    points3d(points,col="green")
  }
  
  #save the hbn to a .nrrd file
  dir.create(savefolder)
  write.nrrd(I, file.path(savefolder, sprintf("%s_%s.nrrd",cell_type,ref)))
}


# write flywire neuron to nrrd, sample = "FAFB14", xyzmatrix
flywire_to_nrrd <- function(flywire_id, cell_type, ref="JRC2018U", savefolder = "data", plot3D =TRUE){
  #read in flywire ID
  flywire_neuron <- read_cloudvolume_meshes(flywire_id)
  
  #transform neuron into the correct template space
  flywire.reg = xform_brain(flywire_neuron*8/1000, reference=ref, sample="FAFB14")
  
  #plot transformed neuron
  if(plot3D){
    nopen3d()
    #plot neuron)
    plot3d(flywire.reg, col = 'black')
    plot3d(get(ref))
  }
  
  # make im3d
  x <- get(ref)
  points = xyzmatrix(flywire.reg)
  I = as.im3d(points,x)
  
  #save the flywire neuron as a .nrrd file
  dir.create(savefolder)
  write.nrrd(I, file.path(savefolder, sprintf("%s_%s.nrrd",cell_type,ref)))
}
