source("R/startup/functions.R")
source("R/startup/packages.R")

# crontab - e
# [specify time] Rscript /my/file/nat-tech.R
# link: https://www.hostinger.com/tutorials/cron-job

# 1 -- find if there is a unprocessed file in your reg folder
#i'm gonna cut corners and just make a folder called "unprocessed"
#at the end of running registration move the file to procesed folder using -  move_files(files, destinations, overwrite = FALSE)
to_register <- as.vector(list.files("/Volumes/Neurobio/Wilson Lab/Emily/unprocessed"))

#check to see if there are unprocessed files, if not quit
if(length(to_register) > 1){
  elements_to_remove = c("processed")
  to_register = to_register[!(list %in% elements_to_remove)]
}else{
  quit(save = "yes")
}

# 2 -- for each file, split into two .nrrd files, one for each channel
# do this by runnng a Fiji macros in this repo using 'system'
# neuronbridger:::runFijiMacro
 #iterate through image files, call fiji macro on each to create file
for (var in to_register) {
  neuronbridger:::runFijiMacro(macro = macro, macroArg = var, headless = FALSE,batch = FALSE,MinMem = MaxMem,MaxMem = "2500m",IncrementalGC = TRUE,Threads = NULL,fijiArgs = NULL,javaArgs = NULL,ijArgs = NULL,fijiPath = fiji.path,DryRun = FALSE)
}

# 3 -- for each file, set up a registration, inc. creating the CMTK command file, which is a .sh
# run the .sh file using system
for(var in to_register){
  munger = write_cmtkreg(var)
  system(munger)
}


# 4 -- read the hemibrain cell type from the file name
# make a .nrrd file in the same space and save with the registered data
#how to incorporate registered brain of choice? must be hard coded?
for(var in to_register){
  celltype = get_cell_type()
  hemibrain_to_nrrd(var)
}

#when done, loop through all of the unprocessed files and move them to the processed folder
for(var in to_register){
  move_files(sprintf("/Volumes/Neurobio/Wilson Lab/Emily/unprocessed%s", var), "/Volumes/Neurobio/Wilson Lab/Emily/unprocessed/processed",overwrite = FALSE)
}
