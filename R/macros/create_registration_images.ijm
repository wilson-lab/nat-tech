//path = getArgument();
//name = File.getName(path);

path = "/Volumes/Neurobio/Wilson Lab/Emily/unprocessed/20220216_JRC2018U_vDeltaA_R11E08AD_R70H05GDBD_3.tif";
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
	File.makeDirectory("/Users/WilsonLab/Desktop/Registration/Images/" + genotype + name_array[5]);
	
}else if(name_array.lenth == 5){
	celltype = name_array[2];
	genotype = name_array[3];
	folder_path = "/" + genotype + name_array[4] + "/";
	File.makeDirectory("/Users/WilsonLab/Desktop/Registration/Images/" + genotype + name_array[4]);
}

//handle 3 channels
list = getList("image.titles");

if(list.length==2){
	selectWindow("C2-" + name);
	//for o2 folder should be "/Volumes/Neurobio/Wilson Lab/Emily/unprocessed/Registration/Images"
	run("Nrrd ... ", "nrrd=/Users/WilsonLab/Desktop/Registration/Images" + folder_path + celltype + "_" + genotype + "01.nrrd");

	selectWindow("C1-" + name);
	//for o2 folder should be "/Volumes/Neurobio/Wilson Lab/Emily/unprocessed/Registration/Images"
	run("Nrrd ... ", "nrrd=/Users/WilsonLab/Desktop/Registration/Images" + folder_path + celltype + "_" + genotype + "02.nrrd");
	
}else if(list.length == 3){
	selectWindow("C1-" + name);
	//for o2 folder should be "/Volumes/Neurobio/Wilson Lab/Emily/unprocessed/Registration/Images"
	run("Nrrd ... ", "nrrd=/Users/WilsonLab/Desktop/Registration/Images" + folder_path + celltype + "_" + genotype + "01.nrrd");

	selectWindow("C2-" + name);
	//for o2 folder should be "/Volumes/Neurobio/Wilson Lab/Emily/unprocessed/Registration/Images"
	run("Nrrd ... ", "nrrd=/Users/WilsonLab/Desktop/Registration/Images" + folder_path + celltype + "_" + genotype + "02.nrrd");
	
	selectWindow("C3-" + name);
	//for o2 folder should be "/Volumes/Neurobio/Wilson Lab/Emily/unprocessed/Registration/Images"
	run("Nrrd ... ", "nrrd=/Users/WilsonLab/Desktop/Registration/Images" + folder_path + celltype + "_" + genotype + "03.nrrd");
}

run("Close All");
exit();
