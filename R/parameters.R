# USER enter the path to your unregistered data here:

#path to where your registration folder is
registration_folder = "~/Desktop/Registration"

#the data you want the pipeline to process, folder can be named anything
data_folder = "~/Desktop/to_register" 

#folder where unprocessed tif files are, located in the data_folder
raw_data = file.path(data_folder,"unprocessed")

#folder that the registered tif files will be moved to, lcoated in the data_folder
processed_data = file.path(data_folder,"processed")

#paths to FIJI macros used to register and create composite images 
macro1 = "~/nat-tech/R/macros/create_registration_images.ijm"
macro2 = "~/nat-tech/R/macros/create_composite.ijm"
macro3 = "~/nat-tech/R/macros/create_max_projection.ijm"

#paths to source the right startup code
packages = "~/nat-tech/R/startup/packages.R"
functs = "~/nat-tech/R/startup/functions.R"

