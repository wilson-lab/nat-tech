#this specific function doesn't work, the purpose is to read in a .nrrd file of a confocal image and plot together with a hemibrain neuron
nrrd_to_hemibrain <- function (file, cell_type){
  print("Still being worked on")
  # #read in nrrd file (i'm just assuming that this doesn't have to be transformed)
  # i <- read.nrrd(file)
  # df <- data.frame(x=attr(i,'x'),y=attr(i,'y'),z=attr(i,'z'))
  # open
  # #get hemibrain neuron and skeletonize
  # hbn <- neuprint_search(sprintf("type:%s",cell_type))
  # hbn_skel <- neuprint_read_neurons(hbn.$bodyid)
  # 
  # #transform from hemibrain to template space (JRC2018U)
  # hbn.jrc2018f = xform_brain(hbn_skel*8/1000, reference="JRC2018F", sample="JFRCIB2018F")
  # hbn.jrc2018u = xform_brain(hbn_skel*8/1000, reference="JRC2018U", sample="JRC2018F")
  # 
  # #plot nrrd and hemibrain neuron
  # nopen3d()
  # plot3d(neuron())
  # plot3d(hbn.jrc2018u, lwd=3, col='black',WithNodes=FALSE,some=FALSE)
  # plot3d(JRC2018U)
  
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
hemibrain_to_nrrd <- function(cell_type, ref="JRC2018U", savefolder = "data", plot3D = FALSE, mesh = TRUE){
  #get hemibrain neuron
  hbn.info <- neuprint_search(sprintf("type:%s.*",cell_type))
  
  #get the neuron skeletons
  message("reading hemibrain skeletons")
  hbn_skels = neuprintr::neuprint_read_neurons(hbn.info$bodyid)
  message('resampling hemibrain skeletons')
  hbn_skels = nat::nlapply(hbn_skels, nat::resample, stepsize = 1)
  if(mesh){
    message('reading hemibrain meshes')
    hbn_neurons = hemibrainr:::hemibrain_neuron_meshes(hbn.info$bodyid)
    points<-rbind(nat::xyzmatrix(hbn_neurons)*(1000/8)/1000,
                  nat::xyzmatrix(hbn_skels))
    
  }else{
    points<-nat::xyzmatrix(hbn_skels)
  }

  #transform hemibrain neuron to template space
  message('transforming hemibrain points to JRCFIB2018F')
  points.reg <- xform_brain(points*8/1000, reference=ref, sample="JRCFIB2018F")
  
  # make im3d
  x <- get(ref)
  I <- nat::as.im3d(points.reg,x)
  
  if (plot3D){
    nopen3d()
    #plot3d(points.reg, lwd=3,col='black',WithNodes=FALSE,soma=FALSE)
    plot3d(x)
    points3d(points,col="black")
  }
  
  #save the hbn to a .nrrd file
  message('writing .nrrd file')
  dir.create(savefolder, showWarnings = FALSE)
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
#returns the name of the folder where the .nrrd files are saved in the Registration folder (should be in the format of AD_GDBD_num)
get_image_folder <- function(file_name){
  name_arr = get_name_array(basename(file_name))
  #removes the .tif from the number at the end of the file name
  no_tif = strsplit(name_arr[(length(name_arr))],"[.]")
  
  if(length(name_arr) == 6){
    return(sprintf("%s_%s_%s",name_arr[4],name_arr[5],no_tif[[1]][1]))
  }else if(length(name_arr) == 5){
    #returns the folder name in the format "Gal4_num"
    return(sprintf("%s_%s",name_arr[4],no_tif[[1]][1]))
  }
}

#correctly formats the date and time based on how the FIJI plugin CMTK gui formats it
time_date_format <- function(){
  time_date = as.character(ymd_hms(Sys.time()))
  i = strsplit(time_date, " ")
  time = strsplit(i[[1]][2], ":")
  return(sprintf("%s_%s.%s.%s",i[[1]][1],time[[1]][1],time[[1]][2],time[[1]][3]))
}

#gets the correct name of the template brain
#if full is true will return the full name of the file, otherwise will just return shortened version
get_registration_brain <- function(file_name, full = FALSE){
  #file_name = basename(file_name)
  #split the name at "_" and set the template name
  name_arr = strsplit(file_name, "_")
  
  #get all of the template brains in the Refbrain folder
  template_folder = list.files('/Users/WilsonLab/Desktop/Registration/Refbrain')
  
  #loop through all of the template brains in the folder and compare to template from the file name
  for(var in template_folder){
    #split template brains by "_" to just get the name
    temp_name = strsplit(var, "_")
    
    #set a temp variable to the first element of the template brain name array (i.e JRC2018U_38um_iso_16bit.nrrd would be just JRC2018U)
    y = temp_name[[1]][1]
    
    #set temp to be if the JRC2018U string exists in the file name, will return true in the location where string matches
    temp = grepl(y, name_arr[[1]], ignore.case = TRUE)
    
    #if true is returned, find the index of true, since i is an array use the length to return correct template
    i = grep("TRUE",temp)
    if(length(i) > 0){
      if(full == TRUE){
        return(var)
      }else{
        return(y)
      }
    }else if(var == "FCWB.nrrd" && length(i) > 0){
      if(full == TRUE){
        return(var)
      }else{
        return("FCWB")
      }
    }
  }
}

#need a function to make the cmtkreg to run in terminal 
#only works for a very specific file name date_celltype_AD_GDBD_num
write_cmtkreg <- function(file_name,
                          template_path = "JRC2018U_38um_iso_16bit.nrrd",
                          registration_folder = "~/Desktop/Registration"){
  #might need to change this
  folder = get_image_folder(file_name)
  date_time = time_date_format()
  save_file_name = sprintf("munger_%s.command", date_time)
  
  
  #for o2 will need to change the save file path
  #file path for munger file -> for o2 should be /Volumes/Neurobio/Wilson Lab/Emily/unprocessed/Registration/commands
  command_folder <- file.path(registration_folder,"Commands")
  image_folder <- file.path(registration_folder, "Images", folder)
  contents <- list.files(image_folder, full.names = TRUE)

  #sets the number of channels to reformat in the munger file, the max is 3 
  if(length(contents) == 2){
    channel_num <- "0102"
  }else if(length(contents) == 3){
    channel_num <- "010203"
  }
  
  dir.create(command_folder, showWarnings = FALSE)
  save_path = file.path(command_folder, save_file_name)
    
  #creates the array of the commands for the cmtk registration
  array = c("#!/bin/bash", 
            sprintf("# %s",date_time), 
            #for 02 gotta change this path to point to the unprocessed folder \"/Volumes/Neurobio/Wilson Lab/Emily/unprocessed/Registration\""
            sprintf("cd \"%s\"",registration_folder), 
            sprintf("\"/Applications/Fiji.app/bin/cmtk/munger\" -b \"/Applications/Fiji.app/bin/cmtk\" -a -w -r %s  -X 26 -C 8 -G 80 -R 4 -A \"--accuracy 0.4\" -W \"--accuracy 0.4\"  -T 4 -s \"Refbrain/%s\" images/%s", channel_num, template_path, folder))
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

#reconnect to the server
server.connect <- function(server="research.files.med.harvard.edu/Neurobio/") {
  user=Sys.getenv('neurobio_user')
  pw=Sys.getenv('neurobio_password')
  system(sprintf('/usr/bin/osascript  -e \'mount volume \"smb://%s:%s@%s/\"\'' ,user,pw,server))
  return(TRUE)  # success
}

#runs the server.connect function multiple times if it fails 
server.connect.cycle <- function(server="research.files.med.harvard.edu/Neurobio/",
                                 attempts = 3, 
                                 sleep.seconds = 60){
  for (i in 1:attempts) {
    res <- try(server.connect(server), silent = TRUE)
    if (!("try-error" %in% class(res))) {
      print("Connected...")
      return(res)
    }
    print(paste0("Attempt #", i, " failed"))
    Sys.sleep(sleep.seconds)
  }
  stop("Maximum number of connection attempts exceeded")
}

#takes file name and returns an array of elements in the correct order
#if cell is true will return the cell type and not the array
get_name_array <- function(file_name, cell = FALSE){
  #name_arr = strsplit(file_name, "_")
  #evaluate how the file name is separated with the correct characters, wil not accept a mix
  if(grepl("_", file_name) & !(grepl("-",file_name))){
    name_arr = strsplit(file_name, "_")
  }else if(grepl("-", file_name) & !(grepl("_",file_name))){
    name_arr = strsplit(file_name, "-")
  }else{
    stop("wrong file name format: separate file name with \'_\' or \'-\' only")
  }
  name_arr = name_arr[[1]]
  
  #get index of image number
  img_num <- grep(".tif", name_arr,ignore.case = TRUE)
  #get index of date
  date <- grep("TRUE",(grepl("^[0-9]+$", name_arr)))
  #get template of template brain
  template <- get_registration_brain(file_name)
  template_loc <- grep(template, name_arr,ignore.case = TRUE)
  
  if(length(name_arr) == 6){
    #get index of AD and DBD 
    AD <- grep("AD", name_arr, ignore.case = TRUE)
    GDBD <- grep("GDBD", name_arr, ignore.case = TRUE)
    
    #get cell type from remaining fragments
    cell_type = name_arr[!name_arr %in% c(name_arr[AD], name_arr[GDBD],name_arr[img_num],name_arr[template_loc],name_arr[date])]
    
    if(cell == TRUE){
      return(cell_type)
    }else{
      #return array in correct order
      c(name_arr[date],template,cell_type,toupper(name_arr[AD]),toupper(name_arr[GDBD]),name_arr[img_num])
    }
  }else if(length(name_arr) == 5){
    #if the name of the file is shorter, assume gal4 line
    gal4 <- grep("Gal4",ignore.case = TRUE, name_arr)
    
    #get cell type from remaining fragments
    cell_type = name_arr[!name_arr %in% c(name_arr[gal4],name_arr[img_num],name_arr[template_loc],name_arr[date])]
    
    if(cell == TRUE){
      return(cell_type)
    }else{
      #return array in correct order
      c(name_arr[date],template,cell_type,name_arr[gal4],name_arr[img_num])
    }
  }else{
    #checks if the file is too short, otherwise will finally throw an error that the file name is incorrectly formatted
    if(length(name_arr) > 6 | length(name_arr) < 5){
      stop("wrong file name format: name too long or too short, name should be in format date_template_celltype_genotype_expnum")
    }else{
      stop("wrong file name format")
    }
  }
}

#changes file name if its not correct already
correct_file_name <- function(full_file, full = FALSE){
  #gets the path
  path = dirname(full_file)
  #gets the array of the name in the correct order
  temp = get_name_array(basename(full_file))
  #gets the correct file name
  new = make_file_name(temp)
  #changes the file name to the correct format
  file.rename(full_file,sprintf("%s/%s",path,new))
  if(full){
    return(new())
  }else{
    return(file.path(path,new))
  }
}

#correctly formats name of the file so that no errors are thrown when sent to FIJI macros
#defaults to using underscores
make_file_name <- function(arr){
  if(length(arr) == 6){
    sprintf("%s_%s_%s_%s_%s_%s", arr[1],arr[2],arr[3],arr[4],arr[5],arr[6])
  }else if(length(arr) == 5){
    sprintf("%s_%s_%s_%s_%s", arr[1],arr[2],arr[3],arr[4],arr[5])
  }
  
}
