args = getArgument();

argArray = split(args, "|");
file1 = argArray[0];
file2 = argArray[1];
savefolder = argArray[2];
score = argArray[3];

name1 = File.getName(file1);
name2 = File.getName(file2);

print("Arguments:");
print(name1);
print(name2);
print(savefolder);
print(score);

//open file
open(file1);
open(file2)

//combine images
// line then em neuron
run("Merge Channels...", "c1=" + name1 +" c2=" + name2 + " create");

//select combined image
selectWindow("Composite");

x = split(name1,".");
y = split(name2,".");
new_name = "MERGE_" +score + "____"+ name1 + "____"+ name2 + ".tif";
//run("Nrrd ...", "nrrd="+new_name);
run("32-bit");
saveAs("Tiff", savefolder + "/" + new_name);


//select original window and close
close();

run("Quit");