arg = getArgument();
//new argument combines the registration folder location and the file location. you need to split at the "-" to get both paths
// for testing: arg = "/Users/wilsonlab/Desktop/Registration/Images-/Volumes/Neurobio/wilsonlab/Emily/unprocessed/20211203_JRC2018U_LT51_VT049899AD_VT034811GDBD_1.tif"
arg_split = split(arg,"-");

//assign registration folder location to reg_folder and .tif location to path
reg_folder = arg_split[0];
path = arg_split[1];
name = File.getName(path);

if(name=="") exit ("No argument");

open(path);
selectWindow(name);
run("Split Channels");
//add functionality for Gal4
//this will mean that naming correctly is essential 20220322(0)_JRCU2018(1)_FC1(2)_VT033647AD(3)_R64A11GDBD(4)_1(5).tif
//20220322(0)_JRCU2018(1)_FC1(2)_R64A11Gal4(3)_1(4).tif
x = split(name,".");
no_tif = x[0];
name_array = split(no_tif,"_");

if(name_array.length == 6){
	celltype = name_array[2];
	genotype = name_array[3] + "_" + name_array[4] + "_";
	folder_path = "/" + genotype + name_array[5] + "/";
	File.makeDirectory(reg_folder + "/" + genotype + name_array[5]);
	
}else if(name_array.lenth == 5){
	celltype = name_array[2];
	genotype = name_array[3];
	folder_path = "/" + genotype + name_array[4] + "/";
	File.makeDirectory(reg_folder + "/" + genotype + name_array[5]);
}

//handle 3 channels
list = getList("image.titles");

if(list.length==2){
	selectWindow("C2-" + name);
	//for o2 folder should be "/Volumes/Neurobio/Wilson Lab/Emily/unprocessed/Registration/Images"
	run("Nrrd ... ", "nrrd=" + reg_folder + "/" + folder_path + celltype + "_" + genotype + "01.nrrd");

	selectWindow("C1-" + name);
	//for o2 folder should be "/Volumes/Neurobio/Wilson Lab/Emily/unprocessed/Registration/Images"
	run("Nrrd ... ", "nrrd=" + reg_folder + "/" + folder_path + celltype + "_" + genotype + "02.nrrd");
	
}else if(list.length == 3){
	selectWindow("C1-" + name);
	//for o2 folder should be "/Volumes/Neurobio/Wilson Lab/Emily/unprocessed/Registration/Images"
	run("Nrrd ... ", "nrrd=" + reg_folder + "/" + folder_path + celltype + "_" + genotype + "01.nrrd");

	selectWindow("C2-" + name);
	//for o2 folder should be "/Volumes/Neurobio/Wilson Lab/Emily/unprocessed/Registration/Images"
	run("Nrrd ... ", "nrrd=" + reg_folder + "/" + folder_path + celltype + "_" + genotype + "02.nrrd");
	
	selectWindow("C3-" + name);
	//for o2 folder should be "/Volumes/Neurobio/Wilson Lab/Emily/unprocessed/Registration/Images"
	run("Nrrd ... ", "nrrd=" + reg_folder + "/" + folder_path + celltype + "_" + genotype + "03.nrrd");
}

run("Close All");
exit();
