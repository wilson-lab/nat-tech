path = getArgument();
name = File.getName(path);

if(name=="") exit ("No argument");

//path = "/Users/wilsonlab/Desktop/Registration/Reformatted/VT008452AD_R35F03GDBD_4/JRC2018U_DNa03_VT008452AD_R35F03GDBD_02_warp_m0g80c8e1e-1x26r4.nrrd";
//name = File.getName(path);

//file name is JRC2018U(0)_vDeltaA(1)_R38D01AD(2)_R70H05GDBD(3)_02(4)_warp(5)_m0g80c8e1e-1x26r4(6).nrrd
//file name is JRC2018U(0)_vDeltaA(1)_R64A11Gal4(2)_02(3)_warp(4)_m0g80c8e1e-1x26r4(5).nrrd
x = File.getParent(path);

//counts the number of files in the folder, use to check for 3 channels
numfiles = getFileList(x);
print(numfiles.length);

//split up the file name in order to find the neuron type and reference brain used
name_array = split(name,"_");
refbrain = name_array[0];
celltype = name_array[1];

//checks if image is gal4 or AD/DBD
//genotype looks like R381D01AD_R70H05GDBD
ch3 = false; 
if(name_array.length > 6){
	genotype = name_array[2] + "_" + name_array[3];
	
	//check if there are 3 channels and then get the file name accordingly
	if(numfiles.length > 3){
		ch_name = name_array[0] + "_" + name_array[1] + "_" + name_array[2] + "_" + name_array[3] + "_03_" + name_array[5] + "_" + name_array[6];
		ch3 = true;
	}
}else if(name_array.length < 7){
	genotype = name_array[2];
	
	//check if there are 3 channels and then get the file name accordingly
	if(numfiles.length > 3){
		ch_name = name_array[0] + "_" + name_array[1] + "_" + name_array[2] + "_03_" + name_array[4] + name_array[5];
		ch3 = true;
	}
}
final = celltype + "_" + genotype + "_" + "composite.tif" ;


//nrrd_name is the name of the hemibrain file which looks like FB2A_JRC2018U.nrrd
nrrd_name = celltype + "_" + refbrain + ".nrrd";

//open registered image
open(path);

//open hemibrain image
open(x + "/" + nrrd_name);


//select hemibrain image
selectWindow(nrrd_name);
setMinAndMax(0, 4);
run("8-bit");

//select reformatted image
selectWindow(name);
bit = bitDepth(); 


//if the reformatted image is 8 bit, change nothing and just merge the two images 
if(bit == 8){
	//check if there are 3 channels and merge the 3 channels if true
	if(ch3){
		open(x + "/" + ch_name);
		selectWindow(ch_name);

		//combine the 3 channels
		run("Merge Channels...", "c2=" + name + " c4=" + ch_name + " c5=" + nrrd_name + " create");
	}else{
		run("Merge Channels...", "c2=" + name + " c5=" + nrrd_name + " create");
	}
//if the image is 16bit, change the saturation, and make the image 8 bit
}else{
	selectWindow(name);
	setMinAndMax(100, 3000);
	run("8-bit");
	
	//check if there are 3 channels and merge the 3 channels if true
	if(ch3){
		//open the 3rd channel and then adjust the exposure
		open(x + "/" + ch_name);
		selectWindow(ch_name);
		setMinAndMax(100, 3000);
		run("8-bit");
		//combine the 3 channels
		run("Merge Channels...", "c2=" + name + " c4=" + ch_name + " c5=" + nrrd_name + " create");
	}else{
		run("Merge Channels...", "c2=" + name + " c5=" + nrrd_name + " create");
	}
}

saveAs("Tiff", x + "/" + final);
run("Close All");
exit();