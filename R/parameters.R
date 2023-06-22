# USER enter the path to your unregistered data here:

#path to where your registration folder is
registration_folder = "/Users/wilsonlab/Desktop/Registration"

#the data you want the pipeline to process
data_folder = "/Volumes/Neurobio/Wilson\ Lab/Emily"

#folder where unprocessed tif files are, located in the data_folder
raw_data = file.path(data_folder,"unprocessed")

#folder that the registered tif files will be moved to, lcoated in the data_folder
processed_data = file.path(data_folder,"processed")

#paths to FIJI macros used to register and create composite images 
macro1 = "/Users/wilsonlab/Documents/GitHub/nat-tech/R/macros/create_registration_images.ijm"
macro2 = "/Users/wilsonlab/Documents/GitHub/nat-tech/R/macros/create_composite.ijm"
macro3 = "/Users/wilsonlab/Documents/GitHub/nat-tech/R/macros/create_max_projection.ijm"
