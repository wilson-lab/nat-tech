# Load libraries
source("R/startup/packages.R")
source("R/startup/functions.R")

# # Get meta data, if this does not work for you obtain data here: https://github.com/flyconnectome/flywire_annotations
# ft <-  read.table("data/Supplemental_file1_annotations.tsv", header=TRUE, sep="\t", quote = "")
# ft.an <- subset(ft, ft$cell_class=="AN")

# Get AN IDs
dna02.proj <- googlesheets4::read_sheet("1URPvS5ysiEX4acKqt3Tx6CiPNkf_F9PtZvKsiv9rXKY", sheet = "flywire")
an.df <- subset(dna02.proj, class=="AN")
an.ids <- unlist(an.df$root_id)
an.ids <- fafbseg::flywire_updateids(an.ids)

# Create .nrrd file and .tiff max projection files
for(an.id in an.ids){
  flywireid_to_nrrd(flywire_id = an.id, 
                    cell_type = an.id, 
                    ref="JRC2018U", 
                    savefolder = "/Volumes/Neurobio/wilsonlab/Fernanda/Ascending Neurons/flywire/", 
                    plot3D = FALSE, 
                    compressed = TRUE,
                    max_projection = TRUE)
}

# Convert my images into maximal projections in 32 bit
contents <- list.files("/Volumes/wilsonlab/Fernanda/Ascending Neurons/red_channel/", full.names = TRUE)
runMacro(macro = "/Users/wilsonlab/Documents/GitHub/nat-tech/R/macros/create_max_projection_serial.ijm", 
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

# Combine our max projections in a single folder
combine_max_projection_tifs(folder1="/Volumes/Neurobio/wilsonlab/Fernanda/Ascending Neurons/flywire/", 
                            folder2="/Volumes/wilsonlab/Fernanda/Ascending Neurons/red_channel/", 
                            savefolder="/Volumes/wilsonlab/Fernanda/Ascending Neurons/an_match_screen")
