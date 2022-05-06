source("R/startup/functions.R")
source("R/startup/packages.R")

#' @title Get unprocessed confocal images and register them to a template brain
#' @description 
#' @
#' format of the file names date_celltype_templatespace_AD_DBD_num (have to think about gal4s)

# crontab - e
# [specify time] Rscript /my/file/nat-tech.R
# link: https://www.hostinger.com/tutorials/cron-job

# 1 -- find if there is a unprocessed file in your reg folder
#i'm gonna cut corners and just make a folder called "unprocessed"
#at the end of running registration move the file to procesed folder using -  move_files(files, destinations, overwrite = FALSE)
raw_data = "/Volumes/Neurobio/Wilson Lab/Emily/unprocessed"
processed_data = c(file.path(raw_data,"processed"), file.path(raw_data,"Registration"))
macro = "macros/create_registration_images.ijm"
to_register <- list.files(raw_data, full.names = TRUE)

#check to see if there are unprocessed files, if not quit
if(length(to_register) > 1){
  to_register = to_register[!to_register %in% processed_data]
}else{
  stop("No files to process")
}

# 2 -- for each file, split into two .nrrd files, one for each channel
# do this by runnng a Fiji macros in this repo using 'system'
# neuronbridger:::runFijiMacro
 #iterate through image files, call fiji macro on each to create file
for (var in to_register) {
  fiji.path = neuronbridger:::fiji()
  runMacro(macro = macro, 
           macroArg = var, 
           headless = FALSE,
           batch = FALSE,
           MinMem = "100m",
           MaxMem = "25000m",
           IncrementalGC = TRUE,
           Threads = NULL,
           fijiArgs = NULL,
           javaArgs = NULL, 
           ijArgs = NULL,
           fijiPath = fiji.path,
           DryRun = FALSE)
  #how to get it to stop running because
  
  # 3 -- for each file, set up a registration, inc. creating the CMTK command file, which is a .sh
  # run the .sh file using system
  file = basename(var)
  template = get_registration_brain(file)
  
  munger_name = write_cmtkreg(var,template)
  system2(command = "sh",
          args = c(munger_name))
  #system2(munger_name)
  #system2(paste("chmod u+r+x ", munger_name))
  #system2(paste("sh ", munger_name))
  
  # 4 -- read the hemibrain cell type from the file name
  # make a .nrrd file in the same space and save with the registered data
  #how to incorporate registered brain of choice? must be hard coded?
  hemibrain_to_nrrd(get_cell_type(file))
  
  #when done, loop through all of the unprocessed files and move them to the processed folder
  move_files(sprintf("/Volumes/Neurobio/Wilson Lab/Emily/unprocessed%s", file), "/Volumes/Neurobio/Wilson Lab/Emily/unprocessed/processed", overwrite = FALSE)
  
}
