# Load libraries
source("R/startup/packages.R")
source("R/startup/functions.R")

# Get directories
save.dir = ""

# Get meta data, if this does not work for you obtain data here: https://github.com/flyconnectome/flywire_annotations
ft <-  read.table("data/Supplemental_file1_annotations.tsv", header=TRUE, sep="\t", quote = "")
ft.an <- subset(ft, ft$cell_class=="AN")

# temporary
ft.an <- ft.an[1:10,]

# Get AN IDs
dna02.proj <- googlesheets4::read_sheet("1URPvS5ysiEX4acKqt3Tx6CiPNkf_F9PtZvKsiv9rXKY", sheet = "flywire")
an.df <- subset(dna02.proj, class=="AN")
an.ids <- an.df$root_id
an.ids <- fafbseg::flywire_updateids(an.ids)

# Create .nrrds file
flywireid_to_nrrd(an.ids, 
                  cell_type = "AN", 
                  ref="JRC2018U", 
                  savefolder = "~/Desktop/Ascending Neurons/", 
                  plot3D =TRUE)
  
