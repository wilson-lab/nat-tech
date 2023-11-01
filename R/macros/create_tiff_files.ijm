path = getArgument();
name = File.getName(path);

if(name=="") exit ("No argument");

//path = "/Users/sophiarenauld/Documents/janelia_acending_neuron_stacks/brain/SS52147-20181026_47_C4-m-20x-brain-Split_GAL4-JRC2018_Unisex_20x_HR-aligned_stack.h5j";
//name = File.getName(arg);

//get folder path
parent_path = File.getParent(path); 

//gets the file list in the parent folder
file_list = getFileList(parent_path);

//open file
open(path);
selectWindow(name);


//split channels
run("Split Channels");

//select red channels
selectWindow("C1-" + name);

//make naming easier
x = split(name,".");
no_tif = x[0];
name_array = split(no_tif,"_");

// save tiff of red channel only
//new_name = no_tif + ".tif" ;
//saveAs("Tiff", parent_path+ "/" + new_name);

new_name = no_tif;
//run("Nrrd ...", "nrrd="+new_name);
run("32-bit");
run("Nrrd ... ", "nrrd=[" + parent_path + "/" + new_name + ".nrrd" + "]");


//select original window and close
//selectWindow(new_name);
close();

run("Quit");
//My macro is not ending
//exit();