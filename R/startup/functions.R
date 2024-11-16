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
hemibrain_to_nrrd <- function(cell_type, ref="JRC2018U", 
                              savefolder = "data", savemaxfolder = savefolder,
                              plot3D = FALSE, mesh = TRUE, max_projection = FALSE,
                              overwrite = FALSE){
  
  # Skip?
  if(!overwrite){
    nrrd.file <- file.path(savefolder, sprintf("%s_%s.nrrd",cell_type,ref))
    if(file.exists(nrrd.file)){
      message("File alreadg exists: ", nrrd.file)
      return(NULL)
    }
  }
  if(!overwrite&max_projection){
    max.file <- file.path(savemaxfolder, sprintf("max_projection_%s_%s.tiff",cell_type,ref))
    if(file.exists(max.file)){
      message("File alreadg exists: ", max.file)
      return(NULL)
    }
  }
  
  #get hemibrain neuron
  hbn.info <- neuprint_search(sprintf("type:%s.*",cell_type))
  
  #get the neuron skeletons
  message("reading hemibrain skeletons")
  hbn_skels <- neuprintr::neuprint_read_neurons(hbn.info$bodyid)
  message('resampling hemibrain skeletons')
  hbn_skels <- nat::nlapply(hbn_skels, nat::resample, stepsize = 1)
  if(mesh){
    message('reading hemibrain meshes')
    hbn_neurons <- hemibrainr:::hemibrain_neuron_meshes(hbn.info$bodyid)
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
  write.nrrd(I, nrrd.file)
  if(max_projection){
    # Maximal projection in X and Y dimensions
    max_projection <- apply(I, c(1, 2), max)
    #write.nrrd(max_projection, file.path(savefolder, sprintf("max_projection_%s_%s.nrrd",cell_type,ref)))
    raster::writeRaster(raster::raster(t(max_projection)), 
                        filename = max.file, 
                        format = "GTiff", 
                        overwrite = overwrite)
  }
}


# write flywire neuron to nrrd, sample = "FAFB14", xyzmatrix
flywireid_to_nrrd <- function(flywire_id, cell_type, ref="JRC2018U", 
                              savefolder = "data", savemaxfolder = savefolder, 
                              plot3D =TRUE, compressed = FALSE, max_projection = TRUE,
                              overwrite = FALSE){
  
  # Skip?
  if(!overwrite){
    nrrd.file <- file.path(savefolder, sprintf("%s_%s.nrrd",cell_type,ref))
    if(file.exists(nrrd.file)){
      message("File alreadg exists: ", nrrd.file)
      return(NULL)
    }
  }
  if(!overwrite&max_projection){
    max.file <- file.path(savemaxfolder, sprintf("max_projection_%s_%s.tiff",cell_type,ref))
    if(file.exists(max.file)){
      message("File alreadg exists: ", max.file)
      return(NULL)
    }
  }
  
  #read in flywire ID
  flywire_neuron <- fafbseg::read_cloudvolume_meshes(flywire_id)
  
  #transform neuron into the correct template space
  flywire.reg <- nat.templatebrains::xform_brain(flywire_neuron, reference=ref, sample="FAFB14")
  
  #plot transformed neuron
  if(plot3D){
    nopen3d()
    #plot neuron)
    plot3d(flywire.reg, col = 'black')
    plot3d(get(ref))
  }
  
  # make im3d
  x <- get(ref)
  points <- xyzmatrix(flywire.reg)
  I <- nat::as.im3d(points,x)
  if(compressed){
    I <- nat::as.im3d(points, voxdims = c(0.519,0.52,1), BoundingBox = boundingbox(x))
  }
  
  #save the flywire neuron as a .nrrd file
  dir.create(savefolder, showWarnings = FALSE)
  write.nrrd(I, nrrd.file)
  
  if(max_projection){
    # Maximal projection in X and Y dimensions
    max_projection <- apply(I, c(1, 2), max)
    #write.nrrd(max_projection, file.path(savefolder, sprintf("max_projection_%s_%s.nrrd",cell_type,ref)))
    raster::writeRaster(raster::raster(t(max_projection)), 
                        filename = max.file, 
                        format = "GTiff",
                        overwrite = overwrite)
  }
  
}

#file format: 20220408(1)_JRCU2018U(2)_FC1(3)_AD(4)_GDBD(5)_01(6).tif(7)
#returns the name of the folder where the .nrrd files are saved in the Registration folder (should be in the format of AD_GDBD_num)
get_image_folder <- function(file_name, reg_folder){
  name_arr = get_name_array(basename(file_name), reg_folder=reg_folder)
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
get_registration_brain <- function(file_name, full = FALSE, reg_folder){
  #file_name = basename(file_name)
  #split the name at "_" and set the template name
  name_arr = strsplit(file_name, "_")
  
  #get all of the template brains in the Refbrain folder
  template_folder = list.files(file.path(reg_folder,"RefBrain")) #list.files('/Users/WilsonLab/Desktop/Registration/Refbrain')
  
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
                          template_path = "JRC2018U_38um_iso_16bit.nrrd", #the default is JRC2018U
                          registration_folder = "~/Desktop/Registration"){ #the default is that your registration folder is on the desktop
  #might need to change this
  folder = get_image_folder(file_name, reg_folder=registration_folder)
  date_time = time_date_format()
  save_file_name = sprintf("munger_%s.command", date_time)
  
  
  #for o2 will need to change the save file path
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
            sprintf("cd \"%s\"",registration_folder), 
            sprintf("\"/Applications/Fiji.app/bin/cmtk/munger\" -b \"/Applications/Fiji.app/bin/cmtk\" -a -w -r %s  -X 26 -C 8 -G 80 -R 4 -A \"--accuracy 0.4\" -W \"--accuracy 0.4\"  -T 8 -s \"Refbrain/%s\" images/%s", channel_num, template_path, folder))
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
get_name_array <- function(file_name, cell = FALSE, reg_folder){
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
  img_num <- grep(".tif", name_arr)
  #get index of date
  date <- grep("TRUE",(grepl("^[0-9]+$", name_arr)))
  #get template of template brain
  template <- get_registration_brain(file_name, reg_folder=reg_folder) # problem - 
  template_loc <- grep(template, name_arr,ignore.case = TRUE)
  
  if(length(name_arr) == 6){
    #get index of AD and DBD 
    AD <- grep("AD", name_arr)
    GDBD <- grep("GDBD", name_arr, ignore.case = TRUE)
    
    #get cell type from remaining fragments
    cell_type = name_arr[!name_arr %in% c(name_arr[AD], name_arr[GDBD],name_arr[img_num],name_arr[template_loc],name_arr[date])]
    
    if(cell == TRUE){
      return(cell_type)
    }else{
      #return array in correct order
      x = c(name_arr[date],template,cell_type,toupper(name_arr[AD]),toupper(name_arr[GDBD]),name_arr[img_num])
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
correct_file_name <- function(full_file, full = FALSE, reg_folder){
  #gets the path
  path = dirname(full_file)
  #gets the array of the name in the correct order
  temp = get_name_array(basename(full_file), reg_folder=reg_folder)
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

# Downlad hemibrain mesh
download_hemibrain_obj <- function (segments, save.obj = getwd(), ratio = 1, ...){
  # Mesh data source
  trimmed_dataset = stringr::str_match("hemibrain:v1.2.1", ".+:v[^.]+\\.[^.]+")
  cloudvolume.url = sprintf("precomputed://gs://neuroglancer-janelia-flyem-%s/segmentation", 
                            sub(":", "/", trimmed_dataset, fixed = T))
  fafbseg:::py_skel_imports()
  fafbseg:::py_cloudvolume(cloudvolume.url, ...)
  neurons = nat::neuronlist()
  message("Saving .obj files in: ", save.obj)
  pb <- progress::progress_bar$new(format = "  downloading [:bar] :current/:total eta: :eta", 
                                   total = length(segments), clear = FALSE, show_after = 1)
  for (id in segments) {
    pb$tick()
    reticulate::py_run_string(sprintf("id=%s", id), ...)
    counter = 3
    while (counter > 0) {
      res = try(reticulate::py_run_string("m = vol.mesh.get(id)[id]", 
                                          ...), silent = TRUE)
      if (class(res)[1] == "try-error") {
        if (grepl("HTTPError|Server", res)) {
          counter = counter - 1
          Sys.sleep(1)
        }
        else {
          break
        }
      }
      else {
        break
      }
    }
    if (class(res)[1] == "try-error") {
      warning(as.character(res))
    }
    else {
      reticulate::py_run_string("m = tm.Trimesh(m.vertices, m.faces)", 
                                ...)
      if (ratio != 1) {
        reticulate::py_run_string(sprintf("m = sk.pre.simplify(m, ratio=%s)", 
                                          ratio), ...)
      }
      ff = file.path(save.obj, paste0(id, ".obj"))
      reticulate::py_run_string(sprintf("s = m.export('%s')", 
                                        ff), ...)
    }
  }
}

combine_max_projection_tifs <- function(folder1, 
                                        folder2, 
                                        savefolder,
                                        scores=NULL){
  if(dir.exists(folder1)){
    files1 <- list.files(folder1, full.names = TRUE, pattern = "^max|^MAX")
  }else{
    files1 <- folder1
  }
  if(dir.exists(folder1)){
    files2 <- list.files(folder2, full.names = TRUE, pattern = "^max|^MAX")
  }else{
    files2 <- folder2
  }
  for(i in 1:length(files1)){
      file1 = files1[i]
      if (is.null(scores)) {
        score=""
      } else {
        score=scores[i]      
      }
    for(file2 in files2){
      message("file 1: ", file1)
      message("file 2: ", file2)
      data1 <- raster::raster(file1)
      data2 <- raster::raster(file2)
      raster::extent(data1) <- raster::extent(data2)
      if (!raster::compareRaster(data1, data2)) {
        stop("The two .tiff files must have the same dimensions.")
      }
      combined <- raster::brick(data1, data2, nl = 2)
      # data1 <- tiff::readTIFF(file1, native = FALSE)
      # data2 <- tiff::readTIFF(file2, native = FALSE)
      # combined <- array(c(data1/max(data1), data2/max(data2)), dim = c(dim(data1),2))
      message("combined")
      # save
      savefile <- file.path(savefolder,paste0("MERGED_",score,gsub("\\.tiff|\\.tif","",basename(file1)),"____",gsub("\\.tiff|\\.tif","",basename(file2)),".tiff"))
      raster::writeRaster(combined, filename = savefile, format = "GTiff")
    }
  }
  return(invisible())
}

combine_nrrds <- function(folder1, 
                          folder2, 
                          savefolder,
                          scores=NULL){
  if(dir.exists(folder1)){
    files1 <- list.files(folder1, full.names = TRUE, pattern = "^max|^MAX")
  }else{
    files1 <- folder1
  }
  if(dir.exists(folder1)){
    files2 <- list.files(folder2, full.names = TRUE, pattern = "^max|^MAX")
  }else{
    files2 <- folder2
  }
  
  for(i in 1:length(files1)){
    file1 = files1[i]
    if (is.null(scores)) {
      score=""
    } else {
      score=scores[i]      
    }
    for(file2 in files2){
      message("file 1: ", file1)
      message("file 2: ", file2)
      runMacro(macro = "/Users/sophiarenauld/Documents/GitHub/nat-tech/R/macros/combine_nrrds.ijm", 
               macroArg = paste(c(file1, file2, savefolder,score), collapse ="|"), 
               headless = FALSE,
               batch = FALSE,
               MinMem = "100m",
               MaxMem = "25000m",
               IncrementalGC = TRUE,
               Threads = NULL,
               fijiArgs = NULL,
               javaArgs = NULL, 
               ijArgs = NULL,
               fijiPath = neuronbridger:::fiji(),
               DryRun = FALSE)
  
    }
  }
  return(invisible())
}

nrrd_bridge_nrrd <- function(nrrd,
                             savefile = NULL,
                             sample = "IS2",
                             reference = "JRC2018U",
                             threshold = 100,
                             plot3D = TRUE,
                             compress = FALSE){
  
  # path to sample .nrrd in IS2 space
  nrrd <- "/Users/abates/Desktop/observations/IS2_TyrRII_no1_02_warp_m0g40c4e1e-1x16r3.nrrd"
  if(is.null(savefile)){
    nf <- basename(nrrd)
    nf <- gsub(sample,"",nf)
    save.file <- file.path(dirname(nrrd), paste0("JRC2018U_", nf))
  }
  array_3d <- read.im3d(nrrd)
  points_df <- data.frame(
    x = rep(1:dim(array_3d)[1], times = dim(array_3d)[2] * dim(array_3d)[3]),
    y = rep(rep(1:dim(array_3d)[2], each = dim(array_3d)[1]), times = dim(array_3d)[3]),
    z = rep(1:dim(array_3d)[3], each = dim(array_3d)[1] * dim(array_3d)[2]),
    value = as.vector(array_3d)
  )
  points_df <- points_df[points_df$value>threshold,]
  
  # Use voxel dimensions
  vd <- nat::voxdims(array_3d)
  xyz <- nat::xyzmatrix(points_df)
  xyz[, 1] = xyz[, 1] * vd[1]
  xyz[, 2] = xyz[, 2] * vd[2]
  xyz[, 3] = xyz[, 3] * vd[3]
  nat::xyzmatrix(points_df) <- xyz
  
  # Move to JRC2018U
  if(reference=="BANC"){
    voxdims <- banc_voxdims()
    ref.brain <- banc_decapitate(banc.surf, invert = TRUE)
    points_df_JRC2018F <- transform_points_chunked(points_df,
                                                   sample = sample,
                                                   reference = "JRC2018F",
                                                   chunk_size = 1000,
                                                   parallel = TRUE)    
    points_df <- transform_points_chunked(points_df_JRC2018F, 
                                          sample = "JRC2018F",
                                          reference = "BANC",
                                          inverse=TRUE, 
                                          chunk_size = 1000,
                                          parallel = TRUE)
  }else{
    ref.brain <- get(reference)
    voxdims <- voxdims(ref.brain)
    points_df <- transform_points_chunked(points_df,
                                                   sample = sample,
                                                   reference = reference,
                                                   chunk_size = 10000)  
  }
  points_df <- points_df %>%
    dplyr::filter(!is.na(x),!is.na(y),!is.na(z))
  
  # Plot
  if(plot3D){
    points3d(xyzmatrix(points_df));plot3d(ref.brain, alpha = 0.1)
  }
  
  # Make im3d, then .nrrd
  points <- nat::xyzmatrix(points_df)
  if(compress & reference=="JRC2018U"){
    I <- nat::as.im3d(points, 
                      voxdims = c(0.519,0.52,1), 
                      origin = nat::origin(nat.flybrains::JRC2018U),
                      BoundingBox = nat::boundingbox(nat.flybrains::JRC2018U))
  }else{
    points <- round(points)
    I <- nat::as.im3d(points, 
                      voxdims = dim(array_3d), 
                      BoundingBox = bbx,
                      origin = nat::origin(ref.brain))
  }
  bbx <- nat::boundingbox(ref.brain)
  inds <- nat::coord2ind(coords = nat::xyzmatrix(points_df),
                         imdims = dim(I),
                         voxdims = voxdims(I),
                         CheckRanges = FALSE,
                         linear.indices = FALSE)
  valid_points <- inds[,1] > 0 & inds[,1] <= dim(I)[1] &
    inds[,2] > 0 & inds[,2] <= dim(I)[2] &
    inds[,3] > 0 & inds[,3] <= bbx[3]
  I[inds[valid_points,]] <- points_df$value[valid_points]

  # Write each slice as a page in the TIFF/Nrrd
  if(grepl("nrrd$",savefile)){
    nat::write.nrrd(I, save.file)
  }else{
    volume <- I/max(I, na.rm = TRUE)
    img <- EBImage::Image(volume)
    EBImage::writeImage(img, save.file) 
  }
  
  # return nothing
  invisible()
  
}

transform_points_chunked <- function(points_df, 
                                     sample, 
                                     reference = "JRC2018U", 
                                     chunk_size = 10000,
                                     parallel = FALSE,
                                     n_workers = 6,
                                     ...) {
  
  # Convert to data.table for faster operations
  points_dt <- data.table::as.data.table(points_df)
  if(reference=="BANC"&!sample%in%c("JRC2018F","JRCVNC2018F")){
    stop("sample must be JRC2018F or JRCVNC2018F")
  }
  
  # Get total number of points
  total_points <- nrow(points_dt)
  n_chunks <- ceiling(total_points / chunk_size)
  
  # Set up parallel processing if requested
  if (parallel) {
    if (is.null(n_workers)) {
      n_workers <- max(1, parallel::detectCores() - 1)
    }
    cl <- parallel::makeCluster(n_workers)
    doParallel::registerDoParallel(cl)
    on.exit(parallel::stopCluster(cl))
    message(sprintf("Using parallel processing with %d workers", n_workers))
  } else {
    foreach::registerDoSEQ()
  }
  
  # Create progress bar
  pb <- progress::progress_bar$new(
    format = "Processing [:bar] :percent eta: :eta",
    total = n_chunks,
    clear = FALSE
  )
  
  # Extract xyz coordinates
  xyz <- nat::xyzmatrix(points_dt)
  
  # Process chunks
  results <- foreach::foreach(i = 1:n_chunks, 
                              .combine = 'c', 
                              .packages = c("nat", "nat.templatebrains", "nat.jrcbrains", "data.table")) %dopar% {
    start_idx <- ((i-1) * chunk_size) + 1
    end_idx <- min(i * chunk_size, nrow(xyz))
    
    chunk_data <- xyz[start_idx:end_idx, ]
    
    if(reference=="BANC"){
      chunk_result <- suppressMessages(
        bancr::banc_to_JRC2018F(
          x = chunk_data, 
          ...
        ))
    }else{
      chunk_result <- suppressMessages(
        nat.templatebrains::xform_brain(
          chunk_data,
          sample = sample,
          reference = reference
        )) 
    }
    
    pb$tick()
    
    list(data.table::as.data.table(chunk_result))
  }
  
  # Combine results
  combined_results <- data.table::rbindlist(results)
  
  # Update points_dt
  if (ncol(combined_results) == 3) {
    data.table::setnames(combined_results, c("x", "y", "z"))
    points_dt[, c("x", "y", "z") := combined_results]
  } else {
    warning("Transformed data does not have exactly 3 columns. Original data not updated.")
  }
  
  return(as.data.frame(points_dt))
}

# Function to process points in chunks
banc_transform_chunked <- function(points_df, 
                                   chunk_size = 10000,
                                   ...) {
  
  # Get total number of points
  total_points <- nrow(points_df)
  n_chunks <- ceiling(total_points / chunk_size)
  
  # Create progress bar
  pb <- progress::progress_bar$new(
    format = "Processing [:bar] :percent eta: :eta",
    total = n_chunks,
    clear = FALSE
  )
  
  # Initialize matrix to store results
  xyz <- nat::xyzmatrix(points_df)
  result_xyz <- matrix(0, nrow = nrow(xyz), ncol = ncol(xyz))
  
  # Process in chunks
  for(i in 1:n_chunks) {
    # Calculate chunk indices
    start_idx <- ((i-1) * chunk_size) + 1
    end_idx <- min(i * chunk_size, total_points)
    
    # Transform chunk
    chunk_result <- suppressMessages(bancr::banc_to_JRC2018F(
      xyz[start_idx:end_idx, ],
      ...
    ))
    
    # Store results
    result_xyz[start_idx:end_idx, ] <- chunk_result
    
    # Update progress bar
    pb$tick()
  }
  
  # Update the points_df with transformed coordinates
  nat::xyzmatrix(points_df) <- result_xyz
  return(points_df)
}


# Function to inflate bounding box by a percentage
inflate_bbox <- function(bbox, inflation_percent = 20) {
  # Convert percentage to multiplier
  multiplier <- inflation_percent / 100
  
  # Calculate center points
  center_x <- (bbox[1,1] + bbox[2,1]) / 2
  center_y <- (bbox[1,2] + bbox[2,2]) / 2
  center_z <- (bbox[1,3] + bbox[2,3]) / 2
  
  # Calculate current dimensions
  width <- bbox[2,1] - bbox[1,1]
  height <- bbox[2,2] - bbox[1,2]
  depth <- bbox[2,3] - bbox[1,3]
  
  # Calculate inflation amounts
  inflate_x <- width * multiplier / 2
  inflate_y <- height * multiplier / 2
  inflate_z <- depth * multiplier / 2
  
  # Create new bounding box
  new_bbox <- structure(c(
    bbox[1,1] - inflate_x, bbox[2,1] + inflate_x,  # x coordinates
    bbox[1,2] - inflate_y, bbox[2,2] + inflate_y,  # y coordinates
    bbox[1,3] - inflate_z, bbox[2,3] + inflate_z   # z coordinates
  ), dim = 2:3, class = "boundingbox")
  
  return(new_bbox)
}

