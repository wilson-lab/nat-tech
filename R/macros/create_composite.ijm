path = getArgument();
name = File.getName(path);

if(name=="") exit ("No argument");

//path = "/Users/wilsonlab/Desktop/Registration/Reformatted/R38D01AD_R70H05GDBD_1/JRC2018U_vDeltaA_R38D01AD_R70H05GDBD_02_warp_m0g80c8e1e-1x26r4.nrrd";
//name = File.getName(path);

//file name is JRC2018U(0)_vDeltaA(1)_R38D01AD(2)_R70H05GDBD(3)_01(4)_warp(5)_m0g80c8e1e-1x26r4.nrrd
x = File.getParent(path);

name_array = split(name,"_");
refbrain = name_array[0];
celltype = name_array[1];
genotype = name_array[2] + "_" + name_array[3];
final = celltype + "_" + genotype + "_" + "composite.tif" ;
nrrd_name = celltype + "_" + refbrain + ".nrrd";


open(path);
open(x + "/" + nrrd_name);
selectWindow(nrrd_name);
run("8-bit");

run("Merge Channels...", "c2=" + name + " c5=" + nrrd_name + " create");

saveAs("Tiff", "/Users/wilsonlab/Desktop/" + final);

run("Close All");
exit();