path = getArgument();
name = File.getName(path);

if(name=="") exit ("No argument");

//path = "/Users/wilsonlab/Desktop/Registration/Reformatted/R38H08AD_VT045265GDBD_2/FB1J_R38H08AD_VT045265GDBD_composite.tif";
//name = File.getName(arg);

//get folder path
parent_path = File.getParent(path); 

//gets the file list in the parent folder
file_list = getFileList(parent_path);

//initialize a place to store the name of the .tif file you want to open
to_composite = ""; 

//goes through the file list and sees what file ends with .tif and opens that file
for (i = 0; i < file_list.length; i++) {
	//checks if given file has .tif extension or jpg or png
	tif_exists = matches(file_list[i], ".*\\.(jpg|png|tif)$");
	if(tif_exists == 1){
		//create the full path that will point to correct .tif file
		full_path = parent_path + "/" + file_list[i];
		//opens .tif file
		open(full_path);
		
		//run max z projection on the to_composite file
    run("Z Project...", "projection=[Max Intensity]");
    run("32-bit");

    to_composite = file_list[i];
    save_name = "MAX_" + to_composite;
    //save as a .tif file with MAX_ and then original name
    saveAs("Tiff", parent_path + "/" + save_name);
    close();
    close();
	}else if ((i == file_list.length - 1) && to_composite == ""){
		//will exit the macro if the for loop has iterated through all files and has not assigned the to_composite variable - meaning there is no .tif file
		exit("no .tif files");
	}
}
