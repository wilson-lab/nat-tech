path = getArgument();
name = File.getName(path);

if(name=="") exit ("No argument");

//path = "/Users/wilsonlab/Desktop/Registration/Reformatted/VT057292AD_R53G06GDBD_1/JRC2018U_FC2C_VT057292AD_R53G06GDBD_02_warp_m0g80c8e1e-1x26r4.nrrd";
//name = File.getName(path);

//file name is JRC2018U(0)_vDeltaA(1)_R38D01AD(2)_R70H05GDBD(3)_02(4)_warp(5)_m0g80c8e1e-1x26r4.nrrd
x = File.getParent(path);

//split up the file name in order to find the neuron type and reference brain used
name_array = split(name,"_");
refbrain = name_array[0];
celltype = name_array[1];

//genotype looks like R381D01AD_R70H05GDBD
genotype = name_array[2] + "_" + name_array[3];
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
	run("Merge Channels...", "c2=" + name + " c5=" + nrrd_name + " create");
	saveAs("Tiff", x + "/" + final);
	run("Close All");
	exit();
//if the image is 16bit, change the saturation, and make the image 8 bit
}else{
	selectWindow(name);
	setMinAndMax(100, 3000);
	run("8-bit");
}

run("Merge Channels...", "c2=" + name + " c5=" + nrrd_name + " create");

saveAs("Tiff", x + "/" + final);

run("Close All");
exit();