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
flywireid_to_nrrd <- function(flywire_id, cell_type, ref="JRC2018U", savefolder = "data", plot3D =TRUE){
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

#file format: 20220408(1)_JRCU2018U(2)_FC1(3)_AD(4)_GDBD(5)_01(6).tif(7)
#take the full file name and returns just the cell type
get_cell_type <- function(file_name){
  file_name = basename(file_name)
  #split the name at "_"
  name_arr = strsplit(file_name, "_")
  return(name_arr[[1]][3])
}

#file format: 20220408(1)_JRCU2018U(2)_FC1(3)_AD(4)_GDBD(5)_01(6).tif(7)
#returns the name of the folder where the .nrrd files are saved in the Registration folder (should be in the format of AD_GDBD_num)
get_image_folder <- function(file_name){
  file_name = basename(file_name)
  #removes the .tif from the number at the end of the file name
  no_tif = strsplit(file_name,"[.]")
  
  #splits up the name string based on underscore
  name_arr = strsplit(no_tif[[1]][1],"_")
  
  #returns the folder name in the format "AD_GDBD_num"
  return(sprintf("%s_%s_%s",name_arr[[1]][4],name_arr[[1]][5],name_arr[[1]][6]))
}

#correctly formats the date and time based on how the FIJI plugin CMTK gui formats it
time_date_format <- function(){
  time_date = as.character(ymd_hms(Sys.time()))
  i = strsplit(time_date, " ")
  time = strsplit(i[[1]][2], ":")
  return(sprintf("%s_%s.%s.%s",i[[1]][1],time[[1]][1],time[[1]][2],time[[1]][3]))
}

#file format: 20220408(1)_JRCU2018U(2)_FC1(3)_AD(4)_GDBD(5)_01(6).tif(7)
#get the template brain from the file name
get_registration_brain <- function(file_name){
  file_name = basename(file_name)
  #split the name at "_" and set the template name
  name_arr = strsplit(file_name, "_")
  template = name_arr[[1]][2]
  
  #get all of the template brains in the Refbrain folder
  template_folder = list.files('/Users/WilsonLab/Desktop/Registration/Refbrain')
  #loop through all of the template brains in the folder and compare to template from the file name
  for(var in template_folder){
    x = strsplit(var, "_")
    y = x[[1]][1]
  
    #if the refbrain from the refbrain folder is equal to the template from the file name return the full name of the template
    if(y == template){
     return(var)
    }else if(template == "FCWB" & y == "FCWB.nrrd"){
      return(var)
    }
  }
}


#need a function to make the cmtkreg to run in terminal 
#only works for a very specific file name date_celltype_AD_GDBD_num
write_cmtkreg <- function(file_name,
                          template_path = "JRC2018U_38um_iso_16bit.nrrd",
                          registration_folder = "/Users/wilsonlab/Desktop/Registration"){
  folder = get_image_folder(file_name)
  date_time = time_date_format()
  save_file_name = sprintf("munger_%s.command", date_time)
  
  
  #for o2 will need to change the save file path
  #file path for munger file -> for o2 should be /Volumes/Neurobio/Wilson Lab/Emily/unprocessed/Registration/commands
  save_path = file.path("/Users/wilsonlab/Desktop/Registration/Commands", save_file_name)
    
  #creates the array of the commands for the cmtk registration
  array = c("#!/bin/bash", 
            sprintf("# %s",date_time), 
            #for 02 gotta change this path to point to the unprocessed folder \"/Volumes/Neurobio/Wilson Lab/Emily/unprocessed/Registration\""
            sprintf("cd \"%s\"",registration_folder), 
            sprintf("\"/Applications/Fiji.app/bin/cmtk/munger\" -b \"/Applications/Fiji.app/bin/cmtk\" -a -w -r 0102  -X 26 -C 8 -G 80 -R 4 -A \"--accuracy 0.4\" -W \"--accuracy 0.4\"  -T 8 -s \"Refbrain/%s\" images/%s", template_path, folder)
  )
  writeLines(array,con=save_path)
  #paste0("sh ", save_path)
  save_path
}

#runs the macro in FIJI
runMacro <- function (macro = "", macroArg = "", headless = FALSE, batch = TRUE, 
                      MinMem = MaxMem, MaxMem = "2500m", IncrementalGC = TRUE, 
                      Threads = NULL, fijiArgs = NULL, javaArgs = NULL, ijArgs = NULL, 
                      fijiPath = fiji(), DryRun = FALSE) {
  os <- neuronbridger:::get_os()
  if (os == "windows") {
    fijiPath <- paste0("\"", fijiPath, "\"")
  }
  if (is.null(Threads) & macroArg == "") {
    simple = TRUE
  }
  else {
    simple = FALSE
  }
  if (headless) 
    fijiArgs = c(fijiArgs, "--headless")
  fijiArgs = paste(fijiArgs, collapse = " ")
  javaArgs = c(paste("-Xms", MinMem, sep = ""), paste("-Xmx", 
                                                      MaxMem, sep = ""), javaArgs)
  if (IncrementalGC) 
    javaArgs = c(javaArgs, "-Xincgc")
  javaArgs = paste(javaArgs, collapse = " ")
  threadAdjust = ifelse(is.null(Threads), "", paste("run(\"Memory & Threads...\", \"parallel=", 
                                                    Threads, "\");", sep = ""))
  ijArgs = paste(c(ijArgs, ifelse(batch, "-batch", "")), collapse = " ")
  if (simple) {
    macroCall = sprintf(" -macro \"%s\"", macro)
    cmd <- paste(fijiPath, javaArgs, fijiArgs, macroCall, 
                 ijArgs)
  }
  else {
    macro
    macroCall = paste(" -eval '", threadAdjust, "runMacro(\"", 
                      macro, "\",\"", macroArg, "\");' ", sep = "")
    cmd <- paste(fijiPath, javaArgs, fijiArgs, "--", macroCall, 
                 ijArgs)
  }
  if (DryRun) 
    return(cmd)
  return(0 == system(cmd))
}