hemibrain_to_JRC2018U_nrrd <- function(cell_type, savefolder = "data", plot3D = TRUE){
  
  #get all FC1 neurons
  hbn.info <- neuprint_search(sprintf("type:%s.*",cell_type))

  #get the neuron skeletons
  hbn_skel = neuprint_read_neurons(hbn.info$bodyid)
  
  #FC1A
  hbn.jrc2018f = xform_brain(hbn_skel*8/1000, reference="JRC2018F", sample="JRCFIB2018F")
  hbn.jrc2018u = xform_brain(hbn.jrc2018f, reference='JRC2018U', sample='JRC2018F')
  
  # make im3d
  points=xyzmatrix(hbn.jrc2018u)
  I=as.im3d(points,JRC2018U)
  
  if (plot3D){
    nopen3d()
    #plot FC1 neurons with different colors
    plot3d(hbn.jrc2018u,lwd=3,col='black',WithNodes=FALSE,soma=FALSE)
    plot3d(JRC2018U)
    points3d(points,col="green")
  }
  
  #save the hbn to a .nrrd file
  dir.create(savefolder)
  write.nrrd(I, file.path(savefolder, sprintf("%s_JRC2018.nrrd",cell_type)))
  
}