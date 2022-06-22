path = getArgument();
name = File.getName(path);

//path = "/Volumes/Neurobio/wilsonlab/Emily/unprocessed/20220216_JRCU2018_vDeltaA_R38D01AD_R70H05GDBD_1.tif";


if(name=="") exit ("No argument");

open(path);
selectWindow(name);
run("Split Channels");

//this will mean that naming correctly is essential 20220322(0)_JRCU2018(1)_FC1(2)_VT033647AD(3)_R64A11GDBD(4)_1(5).tif
x = split(name,".");
no_tif = x[0];
name_array = split(no_tif,"_");
celltype = name_array[2];
genotype = name_array[3] + "_" + name_array[4] + "_";

folder_path = "/" + genotype + name_array[5] + "/";

//make a new directory for each image
File.makeDirectory("/Users/WilsonLab/Desktop/Registration/Images/" + genotype + name_array[5]);

selectWindow("C2-" + name);
//for o2 folder should be "/Volumes/Neurobio/Wilson Lab/Emily/unprocessed/Registration/Images"
run("Nrrd ... ", "nrrd=/Users/WilsonLab/Desktop/Registration/Images" + folder_path + celltype + "_" + genotype + "01.nrrd");

selectWindow("C1-" + name);
//for o2 folder should be "/Volumes/Neurobio/Wilson Lab/Emily/unprocessed/Registration/Images"
run("Nrrd ... ", "nrrd=/Users/WilsonLab/Desktop/Registration/Images" + folder_path + celltype + "_" + genotype + "02.nrrd");
run("Close All");
exit();
