# Load libraries
source("R/startup/packages.R")
source("R/startup/functions.R")

#paths we will need
imagepath = "/Users/sophiarenauld/Documents/janelia_acending_neuron_stacks/lines_brain"
lineimages = "/Users/sophiarenauld/Documents/janelia_acending_neuron_stacks/rgl_plots"
flywirepath ="/Users/sophiarenauld/Documents/janelia_acending_neuron_stacks/flywire_brain"
processed_ans = "/Users/sophiarenauld/Documents/janelia_acending_neuron_stacks/processed_ans"
matching_folder = "/Users/sophiarenauld/Documents/janelia_acending_neuron_stacks/matched_ans"

#read in AN data & turn into nblast objects
# Get meta data, if this does not work for you obtain data here: https://github.com/flyconnectome/flywire_annotations
ft <-  read.table("data/Supplemental_file1_annotations.tsv", header=TRUE, sep="\t", quote = "", colClasses = "character")
ft.an <- subset(ft, ft$cell_class=="AN")
#ft.an$root_id <- as.character(ft.an$root_id)

#read in AN as mesh file and turn into vector cloud (can do lots of things with meshes)
#REMEMBER that the light level data is not actually in JRC2018U  
ans.dps = nat::neuronlist()
for (id in ft.an$root_id) {
  mesh = fafbseg::read_cloudvolume_meshes(id)
  #transform neuron into the correct template space
  flywire.reg = nat.templatebrains::xform_brain(mesh, reference="JRC2018U", sample="FAFB14")
  flywire.reg[[1]] = Rvcg::vcgQEdecim(flywire.reg[[1]], percent = 0.001)
  dps = nat::dotprops(flywire.reg)
  ans.dps = nat::union(ans.dps, dps)
}
# save all em ans as dots
savefile = file.path(processed_ans, "ans_dps.rda")
save(ans.dps, file = savefile)

#step 2 - get dot masks of light level data
#get image data
lm.images = list.files(imagepath, full = TRUE, pattern = "\\.h5j")
#regex - match strings, have to remove contexts of characters (\\)
#no longer need contents <- list.files(lineimages, full.names = TRUE)
for (pic in lm.images[1:25]) {
  runMacro(macro = "/Users/sophiarenauld/Documents/GitHub/nat-tech/R/macros/create_tiff_files.ijm", 
           macroArg = pic,
           headless = FALSE,
           batch = FALSE,
           MinMem = "100m",
           MaxMem = "25000m",
           IncrementalGC = TRUE,
           Threads = NULL,
           fijiArgs = NULL,
           javaArgs = NULL, 
           ijArgs = NULL,
           fijiPath = neuronbridger:::fiji(),
           DryRun = FALSE)
}

#create new max projections of light data
contents <- list.files(imagepath, full.names = TRUE)
runMacro(macro = "/Users/sophiarenauld/Documents/GitHub/nat-tech/R/macros/create_max_projection_serial.ijm", 
         macroArg = contents[2], 
         headless = FALSE,
         batch = FALSE,
         MinMem = "100m",
         MaxMem = "25000m",
         IncrementalGC = TRUE,
         Threads = NULL,
         fijiArgs = NULL,
         javaArgs = NULL, 
         ijArgs = NULL,
         fijiPath = neuronbridger:::fiji(),
         DryRun = FALSE)

lm.ans.dps = nat::neuronlist()

nrrd.images = unique(list.files(imagepath, full = TRUE, pattern = "\\.nrrd"))
for (file in nrrd.images) {
  message("working on ", file)
  pc = nat::read.im3d(file)
  xlookup = attr(pc,"x")
  ylookup = attr(pc,"y")
  zlookup = attr(pc,"z")
  m = reshape2::melt(pc, varnames = c("Row", "Column", "Dimension"), value.name = "brightness")
  points = m %>%
    dplyr::filter(brightness > stats::quantile(pc, 0.9999)) %>%
    dplyr::rename(x=Row, y=Column, z=Dimension) %>%
    dplyr::mutate(x=xlookup[x], y=ylookup[y], z=zlookup[z]) %>%
    dplyr::select(-brightness)
  ## Sanity check
  #I <- nat::as.im3d(nat::xyzmatrix(points), voxdims = c(0.519,0.52,1), BoundingBox = nat::boundingbox(nat.flybrains::JRC2018U))
  #nat::write.nrrd(I, "~/Desktop/test.nrrd")
  dps = nat::dotprops(points)
  entry = nat::neuronlist(dps, DATAFRAME = data.frame(file=basename(file)))
  names(entry) = basename(file)
  lm.ans.dps = nat::union(lm.ans.dps, entry)
  open3d(zoom=0.4581118, userMatrix = structure(c(0.99932849407196, -0.00221055885776877, -0.0365745350718498, 
                                                  0, -0.0118458131328225, -0.964066743850708, -0.265395909547806, 
                                                  0, -0.0346736125648022, 0.2656509578228, -0.963445603847504, 
                                                  0, 0, 0, 0, 1), .Dim = c(4L, 4L)), windowRect = c(0L, 45L, 1425L, 812L))
  points3d(points)
  plot3d(JRC2018U, alpha =0.1, color="gray")
  rgl.snapshot(file=file.path(lineimages, gsub("\\.nrrd", ".png", basename(file))))
  close3d()
}

#save all light level ans as dots
savefile = file.path(processed_ans, "lm_ans_dps.rda")
save(lm.ans.dps, file = savefile)

#step 3 - run Nblast
#uncomment in new page
load(file.path(processed_ans, "lm_ans_dps.rda"))
load(file.path(processed_ans, "ans_dps.rda"))

#run nblast
nb = nat.nblast::nblast(lm.ans.dps, ans.dps, normalised = TRUE)
colnames(nb) <- paste0("id",as.character(colnames(nb)))
rownames(nb) <- paste0("id",as.character(rownames(nb)))
nb.df <- reshape2::melt(nb, varnames = c("Row", "Column"), value.name = "nblast_score")
nb.hits.df <- nb.df %>%
  dplyr::rename(flywire = Row,
         line_name = Column) %>%
  dplyr::mutate(flywire = gsub("^id","",flywire),
          line_name = gsub("\\-.*", "", gsub("^id","",line_name)),
         proceed = nblast_score>0.1) %>%
  dplyr::filter(flywire!=line_name) %>%
  dplyr::arrange(line_name, dplyr::desc(nblast_score), proceed) %>%
  as.data.frame()
readr::write_csv(nb.hits.df, file = file.path(processed_ans, "ans_nblast_results.csv"))

#save just the an em images with a hit
top.ids = subset(nb.hits.df, nblast_score > 0.1)$flywire
for (an.id in top.ids){
  flywireid_to_nrrd(flywire_id = an.id, 
                    cell_type = an.id, 
                    ref="JRC2018U", 
                    savefolder = file.path(flywirepath,"image_stack"), 
                    savemaxfolder = file.path(flywirepath,"max_projection"),
                    plot3D = FALSE, 
                    compressed = TRUE,
                    max_projection = TRUE)
}

#done up until here !!!
#step 4 - save subset of nblast light level images
top.lines = subset(nb.hits.df, nblast_score > 0.1)$line_name
for (ln in top.lines){
  line.folder=file.path(matching_folder, ln)
  dir.create(line.folder, showWarnings = FALSE)
  dir.create(file.path(line.folder,"max_projection"),  showWarnings = FALSE)
  dir.create(file.path(line.folder,"image_stack"),  showWarnings = FALSE)
  # thing 1 move the em max projections and stacks into folder
  top.ids.lines = subset(nb.hits.df, nblast_score > 0.1 & line_name == ln)$flywire
  found.files = list.files(file.path(flywirepath,"max_projection"), full = TRUE, pattern = paste(top.ids.lines, collapse = "|"))
  found.files.3d = list.files(file.path(flywirepath,"image_stack"), full = TRUE, pattern = paste(top.ids.lines, collapse = "|"))
  # thing 2 make max projections of light data
  ln.max = list.files(imagepath, full = TRUE, pattern = paste0("^MAX.*", ln, ".*\\.tif"))
  ln.3d = list.files(imagepath, full = TRUE, pattern = paste0(ln, ".*\\.nrrd"))
  scores = subset(nb.hits.df, nblast_score > 0.1 & line_name == ln)$nblast_score
  scores = round(scores*100, 0)
  # thing 3 combine max projections of em & light data
  combine_max_projection_tifs(found.files,
                              ln.max,
                              savefolder = file.path(line.folder,"max_projection"),
                              scores = scores)
  combine_nrrds(found.files.3d,
                ln.3d,
                savefolder = file.path(line.folder,"image_stack"),
                scores = scores)
}

#turn into colorMIP file
neuronbridger::nrrd_to_mip()


#get transparent brain and plot neuron of interest
#plot3d(elmr::FAFB.surf, alpha = 0.1, col = "gray")
#plot3d(JRC2018U) for switching brain space
#> plot3d(mesh)
#> 
