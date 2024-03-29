# You can't run this code in terminal without full path to source these files since the terminal can't see these files
# change these lines to your own user: Example - "/Users/[INSERT USER]/Documents/GitHub/nat-tech/R/parameters.R"

#### CHANGE THIS LINE #####
source("~/nat-tech/R/parameters.R")

#sources the packages and functions needed to run this code, the paths can be found in parameters.R
source(packages)
source(functs)

# 0 -- set up crontab on your computer to run the script automatically(crontab -e)
# example below: run script at 4:55pm Mon-Fri, create a log to review, point to where you have put this repository
# 02 17 * * 1-5  /usr/local/bin/Rscript /Users/[INSERT USER]/Documents/GitHub/nat-tech/R/pipeline.R > /Users/[INSERT USER]/Documents/GitHub/nat-tech/R/jobs/day.log 2>&1

# 1 -- find if there is a unprocessed file in the unprocessed folder on the server
# at the end of running registration move the file to processed folder using -  move_files(files, destinations, overwrite = FALSE)
to_register <- list.files(raw_data, full.names = TRUE)
if(length(to_register) == 0){
  stop("No files to process")
}

# 2 -- for each file, split into two .nrrd files, one for each channel and save in the correct folder
# iterate through each unprocessed image, register, pull hemibrain neuron, move raw image to processed folder and create composite
for (var in to_register) {
  
  #changes the name of the file to the correct format
  var = correct_file_name(var, reg_folder = registration_folder)
  
  #sends the fiji macro the location of the Registration folder 
  var_temp = paste(registration_folder, "/Images", sep="")
  var_temp = paste(var_temp, var,sep = "-")
  
  # 2.5 -- runs a fiji macro to split your .tif file into the correct .nrrd files and saves them in the Registration folder under Images
  # saves the images in a folder - AD_GDBD_expnum
  fiji.path = neuronbridger:::fiji()
  runMacro(macro = macro1, 
           macroArg = var_temp, 
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
  template = get_registration_brain(var, full =TRUE, reg_folder=registration_folder)
  
  # create the CMTK registration .sh file
  munger_name = write_cmtkreg(var,template_path = template,registration_folder=registration_folder)
  
  # run the .sh file using system
  system(paste0("sh ", munger_name))
  
  # 4 -- move all of the unprocessed files into processed folder
  temp = get_image_folder(var)
  
  # 5 -- check to see if the server is still connected, if not, reconnects to the server
  server.exists <- dir.exists(raw_data)
  message("Can I see the server? ", server.exists)
  if(!server.exists){
    server.connect.cycle()
    server.exists <- dir.exists(raw_data)
    message("Can I see the server now? ", server.exists)
  }
  
  #get the file name so you can move the file to the processed folder
  file = basename(var)
  move_files(file.path(raw_data, file), 
             processed_data, overwrite = FALSE)
  
  # 6 -- gets the the contents of the reformatted folder, should include both reformatted confocal channels
  reformatted_folder <- file.path(registration_folder,"Reformatted")
  dir.create(reformatted_folder, showWarnings = FALSE)
  reformatted_exp_folder <- file.path(reformatted_folder,temp)
  contents <- list.files(reformatted_exp_folder, full.names = TRUE)
  
  # 7 -- read in the hemibrain cell type from the file name
  # puts the hemibrain .nrrd file in the same space registered confocal files AFTER fetching both confocal files
  hemibrain_to_nrrd(cell_type = get_name_array(file,cell = TRUE), savefolder=reformatted_exp_folder, plot3D = FALSE)
  
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
  
  # 8 -- create a max z projection of the image you just made and save in the Registration/Reformatted folder
  runMacro(macro = macro3, 
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
  message("got to macro3")
}
