setwd("/Users/wilsonlab/Documents/GitHub/nat-tech/R/")
source("startup/functions.R")
source("startup/packages.R")

# run R script at specified time 
# runs at 8pm Mon-Fri (supposedly)
# 0 20 * * 1-5  Rscript /Users/wilsonlab/Documents/GitHub/nat-tech/R/pipeline.R

# 1 -- find if there is a unprocessed file in the unprocessed folder on the server
# at the end of running registration move the file to processed folder using -  move_files(files, destinations, overwrite = FALSE)
raw_data = "/Volumes/Neurobio/Wilson Lab/Emily/unprocessed"
processed_data = c(file.path(raw_data,"processed"), file.path(raw_data,"Registration"))
macro1 = "macros/create_registration_images.ijm"
macro2 = "macros/create_composite.ijm"
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
  runMacro(macro = macro1, 
           macroArg = var, 
           headless = TRUE,
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
  
  # 3 -- for each confocal file, set up a registration,
  file = basename(var)
  #get the template brain from the file name
  template = get_registration_brain(file)
  
  #create the CMTK registration .sh file
  munger_name = write_cmtkreg(var,template)
  
  # run the .sh file using system
  system2(command = "sh", args = c(munger_name))
  
  # 4 -- read in the hemibrain cell type from the file name
  # saves the hemibrain .nrrd file in the same space as the registered confocal file (specifically the warp)
  temp = get_image_folder(file)
  hemibrain_to_nrrd(get_cell_type(file), sprintf(savefolder="/Users/wilsonlab/Desktop/Registration/Reformatted/%s", temp))
  
  #move all of the unprocessed files into processed folder
  move_files(sprintf("/Volumes/Neurobio/Wilson Lab/Emily/unprocessed/%s", file), "/Volumes/Neurobio/Wilson Lab/Emily/unprocessed/processed", overwrite = FALSE)
  
  #gets the the contents of the reformatted folder, should include both reformatted confocal channels + the hemibrain .nrrd file
  contents = list.files(sprintf("/Users/wilsonlab/Desktop/Registration/Reformatted/%s",temp), full.names = TRUE)
  
  #creates a composite of registered image and saves to the desktop
  #the problem here is that the final registered images probably won't exist by the time this tries to execute
  runMacro(macro = macro2, 
           macroArg = contents[2], 
           headless = TRUE,
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
}
