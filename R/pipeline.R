source("/Users/wilsonlab/Documents/GitHub/nat-tech/R/startup/packages.R")
source("/Users/wilsonlab/Documents/GitHub/nat-tech/R/startup/functions.R")

# 0 -- set up crontab on your computer to ru the script (crontab -e)
# example below: run scripte at 4:55pm Mon-Fri, create a log to review, point to where you have put this repository
# 55 16 * * 1-5  /usr/local/bin/Rscript /Users/[user name]/Documents/GitHub/nat-tech/R/pipeline.R > /Users/[user name]/Documents/GitHub/nat-tech/R/jobs/day.log 2>&1

# 1 -- find if there is a unprocessed file in the unprocessed folder on the server
# at the end of running registration move the file to processed folder using -  move_files(files, destinations, overwrite = FALSE)
raw_data = "/Volumes/Neurobio/Wilson\ Lab/Emily/unprocessed"
processed_data = c(file.path(raw_data,"processed"), file.path(raw_data,"Registration"))
macro1 = "/Users/wilsonlab/Documents/GitHub/nat-tech/R/macros/create_registration_images.ijm"
macro2 = "/Users/wilsonlab/Documents/GitHub/nat-tech/R/macros/create_composite.ijm"
to_register <- list.files(raw_data, full.names = TRUE)

# check to see if there are unprocessed files, if not quit
if(length(to_register) > 1){
  to_register = to_register[!to_register %in% processed_data]
}else{
  stop("No files to process")
}

# 2 -- for each file, split into two .nrrd files, one for each channel and save in the correct folder
# iterate through each unprocessed image, register, pull hemibrain neuron, move raw image to processed folder and create composite
for (var in to_register) {
  
  #changes the name of the file to the correct format
  var = correct_file_name(var)
  
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
  # get the template brain from the file name
  template = get_registration_brain(var, full = TRUE)
  
  # create the CMTK registration .sh file
  munger_name = write_cmtkreg(var,template_path = template)
  
  # run the .sh file using system
  system(paste0("sh ", munger_name))
  
  # 4 -- move all of the unprocessed files into processed folder
  temp = get_image_folder(var)
  
  # 5 -- check to see if the server is still connected, if not, reconnects to the server
  server.exists <- dir.exists("/Volumes/Neurobio/Wilson Lab/Emily/")
  message("Can I see the server? ", server.exists)
  if(!server.exists){
    server.connect.cycle()
    server.exists <- dir.exists("/Volumes/Neurobio/Wilson Lab/Emily/")
    message("Can I see the server now? ", server.exists)
  }
  
  #get the file name so you can move the file to the processed folder
  file = basename(var)
  move_files(sprintf("/Volumes/Neurobio/Wilson Lab/Emily/unprocessed/%s", file), 
             "/Volumes/Neurobio/Wilson Lab/Emily/unprocessed/processed", overwrite = FALSE)
  
  # 6 -- gets the the contents of the reformatted folder, should include both reformatted confocal channels
  contents = list.files(sprintf("/Users/wilsonlab/Desktop/Registration/Reformatted/%s",temp), full.names = TRUE)
  
  # 7 -- read in the hemibrain cell type from the file name
  # puts the hemibrain .nrrd file in the same space registered confocal files AFTER fetching both confocal files
  hemibrain_to_nrrd(cell_type = get_name_array(var,cell = TRUE), savefolder=sprintf("/Users/wilsonlab/Desktop/Registration/Reformatted/%s", temp), plot3D = FALSE)
  
  # 8 -- create a composite of the hemibrain neuron and the confocal image and save in the Registration/Reformatted folder
  # creates a composite of registered image and hemibrain neuron
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
