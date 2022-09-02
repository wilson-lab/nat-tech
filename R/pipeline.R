#setwd("/Users/wilsonlab/Documents/GitHub/nat-tech/")
# 4 -- read in the hemibrain cell type from the file name
# saves the hemibrain .nrrd file in the same space as the registered confocal file (specifically the warp)
#temp = get_image_folder(file)
source("/Users/wilsonlab/Documents/GitHub/nat-tech/R/startup/packages.R")
source("/Users/wilsonlab/Documents/GitHub/nat-tech/R/startup/functions.R")
message("imported functions")

# 0 -- set up crontab on your computer to ru the script (crontab -e)
# run at 5pm Mon-Fri
# 0 17 * * *  /usr/local/bin/Rscript /Users/wilsonlab/Documents/GitHub/nat-tech/R/pipeline.R > /Users/wilsonlab/Documents/GitHub/nat-tech/R/jobs/day.log 2>&1
# 50 16 * * 1-5  /usr/local/bin/Rscript /Users/wilsonlab/Documents/GitHub/nat-tech/R/pipeline.R > /Users/wilsonlab/Documents/GitHub/nat-tech/R/jobs/day.log 2>&1

# 1 -- find if there is a unprocessed file in the unprocessed folder on the server
# at the end of running registration move the file to processed folder using -  move_files(files, destinations, overwrite = FALSE)
raw_data = "/Volumes/Neurobio/Wilson\ Lab/Emily/unprocessed"
processed_data = c(file.path(raw_data,"processed"), file.path(raw_data,"Registration"))
macro1 = "/Users/wilsonlab/Documents/GitHub/nat-tech/R/macros/create_registration_images.ijm"
macro2 = "/Users/wilsonlab/Documents/GitHub/nat-tech/R/macros/create_composite.ijm"
to_register <- list.files(raw_data, full.names = TRUE)
message(length(to_register))
# check to see if there are unprocessed files, if not quit
if(length(to_register) > 1){
  to_register = to_register[!to_register %in% processed_data]
}else{
  stop("No files to process")
}
message(length(to_register))
# 2 -- for each file, split into two .nrrd files, one for each channel and save in the correct folder
# iterate through each unprocessed image, register, pull hemibrain neuron, move raw image to processed folder and create composite
for (var in to_register) {
  print(var)
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
  
  # create the CMTK registration .sh file
  munger_name = write_cmtkreg(var,template)
  
  # run the .sh file using system
  system(paste0("sh ", munger_name))
  
  # 4 -- move all of the unprocessed files into processed folder
  temp = get_image_folder(file)
  server.exists <- dir.exists("/Volumes/Neurobio/Wilson Lab/Emily/")
  message("Can I see the server? ", server.exists)
  if(!server.exists){
    server.connect.cycle()
    server.exists <- dir.exists("/Volumes/Neurobio/Wilson Lab/Emily/")
    message("Can I see the server now? ", server.exists)
  }
  move_files(sprintf("/Volumes/Neurobio/Wilson Lab/Emily/unprocessed/%s", file), 
             "/Volumes/Neurobio/Wilson Lab/Emily/unprocessed/processed", overwrite = FALSE)
  
  # 4.5 -- gets the the contents of the reformatted folder, should include both reformatted confocal channels
  contents = list.files(sprintf("/Users/wilsonlab/Desktop/Registration/Reformatted/%s",temp), full.names = TRUE)
  
  # 5 -- read in the hemibrain cell type from the file name
  # puts the hemibrain .nrrd file in the same space registered confocal files AFTER fetching both confocal files
  hemibrain_to_nrrd(cell_type = get_cell_type(file), savefolder=sprintf("/Users/wilsonlab/Desktop/Registration/Reformatted/%s", temp), plot3D = FALSE)
  message("ran hemibrain_to_nrrd funtion")
  
  #6 -- create a composite of the hemibrain neuron and the confocal image and save in the Registration/Reformatted folder
  #creates a composite of registered image and hemibrain neuron
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
  message("ran macro2")
}
