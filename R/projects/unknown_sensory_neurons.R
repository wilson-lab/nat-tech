# Load libraries
source("R/startup/packages.R")
source("R/startup/functions.R")

# Working directory
flywirenrrdfolder = "/Volumes/neurobio/wilsonlab/asbates/neuroanat/flywire_nrrd"
flywiremaxfolder = "/Volumes/neurobio/wilsonlab/asbates/neuroanat/flywire_mp"
flywiremipfolder = "/Volumes/neurobio/wilsonlab/asbates/neuroanat/flywire_mip"

# Create folders
dir.create(flywirenrrdfolder, showWarnings = FALSE)
dir.create(flywiremaxfolder, showWarnings = FALSE)
dir.create(flywiremaxfolder, showWarnings = FALSE)

# Get flywire data
ft <- fafbseg::flytable_query("select _id, root_id, root_630, supervoxel_id, proofread, status, pos_x, pos_y, pos_z, nucleus_id, side, ito_lee_hemilineage, hartenstein_hemilineage, top_nt, top_nt_conf, flow, super_class, cell_class, cell_sub_class, cell_type, hemibrain_type, root_duplicated from info")
ft.us <- ft %>%
  dplyr::filter(!super_class %in% 
                  c("not_neuron","glia","Kenyon_cell","trachea","tadpole",
                    "large_fragment","fragment")) %>%
  dplyr::filter(super_class =="sensory",
                cell_sub_class %in% c("unknown", 
                                     "accessory_pharyngeal_nerve_sensory_group2", 
                                     "pharyngeal_nerve_sensory_group1", 
                                     "pharyngeal_nerve_sensory_group2"),
                !duplicated(root_id)) %>%
  dplyr::mutate(cell_type = ifelse(is.na(cell_type),"",cell_type),
                cell_class = ifelse(is.na(cell_class),"",cell_class),
                cell_sub_class = ifelse(is.na(cell_sub_class),"",cell_sub_class)) %>%
  as.data.frame()

for(csc in unique(ft.us$cell_sub_class)){
  csc.df <- subset(ft.us, ft.us$cell_sub_class %in% csc)
  for(cc in unique(ft.us$cell_class)){
    cc.df <- subset(csc.df, ft.us$cell_class %in% cc)
    ids <- unique(csc.df$root_630)
    
    # save flywire .nrrd files
    flywireid_to_nrrd(flywire_id = ids, 
                      cell_type = paste0(csc,"_",cc), 
                      ref="JRC2018U", 
                      savefolder = flywirenrrdfolder,
                      savemaxfolder = flywiremaxfolder,
                      plot3D = FALSE, 
                      compressed = FALSE,
                      max_projection = TRUE)
    
    # save MIP files
  }
}

# And those that were specifically identified
flange.inputs.1 <- unique(c("720575940645493283", "720575940626396803", "720575940625871181",
                     "720575940624910373", "720575940620910716", "720575940638357813",
                     "720575940614641202", "720575940630434616", "720575940624905572",
                     "720575940611959011", "720575940632040749", "720575940625785731",
                     "720575940621106465", "720575940623951591", "720575940606163842",
                     "720575940622575754", "720575940605682790", "720575940604355872",
                     "720575940617034713", "720575940629887695", "720575940633759840",
                     "720575940653233569", "720575940639786830", "720575940626311049"))
flange.inputs.2 <- unique(c("720575940634085786", "720575940644676900", "720575940626767152", 
                     "720575940645108936", "720575940606239666", "720575940617095005", 
                     "720575940626522910", "720575940611643033", "720575940645108936", 
                     "720575940634085786", "720575940645108936", "720575940644676900", 
                     "720575940606239666", "720575940617095005", "720575940617095005", 
                     "720575940633579360", "720575940644676900", 
                     "720575940606239666", "720575940626767152", "720575940634085786", 
                     "720575940617095005"))

# Save flywire .nrrd files
flywireid_to_nrrd(flywire_id = flange.inputs.1, 
                  cell_type = "flange_inputs_1", 
                  ref="JRC2018U", 
                  savefolder = flywirenrrdfolder,
                  savemaxfolder = flywiremaxfolder,
                  plot3D = FALSE, 
                  compressed = TRUE,
                  max_projection = TRUE)
flywireid_to_nrrd(flywire_id = flange.inputs.2, 
                  cell_type = "flange_inputs_2", 
                  ref="JRC2018U", 
                  savefolder = flywirenrrdfolder,
                  savemaxfolder = flywiremaxfolder,
                  plot3D = FALSE, 
                  compressed = TRUE,
                  max_projection = TRUE)

# Save mip files
neuronbridger::nrrd_to_mip(fiji.path = neuronbridger:::fiji())

