# Load libraries
source("R/startup/packages.R")
source("R/startup/functions.R")
library(tiff)

# Function to detect Greek letters
contains_greek <- function(s) {
  grepl("[\u0370-\u03FF]", s)
}

# Function to replace Greek letters with Latin letters
replace_greek_with_latin <- function(s) {
  greek_to_latin <- list(
    "\u03B1" = "a",  # α -> a
    "\u03B2" = "b",  # β -> b
    "\u03A9" = "o"   # Ω -> O
    # Add more mappings as needed
  )
  
  for (greek in names(greek_to_latin)) {
    latin <- greek_to_latin[[greek]]
    s <- gsub(greek, latin, s)
  }
  return(s)
}


# Download the kondo et al. files
# Download all the Kondo et al., 2020 files
options(timeout = max(2000, getOption("timeout")))
all.receptors <- c("5-HT1B",
               "5-HT2A",
               "5-HT2B",
               "5-HT7",
               "AkhR",
               "AstA-R1",
               'AstA-R2',
               "AstC-R1",
               "AstC-R2",
               'CCAP-R',
               "CCHa1-R",
               'CCHa2-R',
               "CCKLR-17D1",
               "CCKLR-17D3",
               'CG12344',
               'CG13229',
               "CG13575",
               "CG13995",
               "CG32547",
               "CG33639",
               'CG43795',
               "CG7589",
               "CapaR",
               "CrzR",
               "Dh-31R",
               "Dh44-R2",
               "Dop1R1",
               'Dop1R2',
               "Dop2R",
               "DopEcR",
               'GABA-B-R1',
               "GABA-B-R2",
               "GluCLα",
               "GluRIB",
               "GluRIIA",
               "GluRIIB",
               "GluRIIC",
               "GluRIID",
               'GluRIIE',
               "Grd",
               'Gyc76C', # not registered to IS2
               "KaiR1D",
               "Lcch3",
               "Lgr1",
               "Lgr3",
               "Lgr4",
               "NPFR",
               "Nmdar1",
               "Nmdar2",
               "Oamb",
               "Octβ2R",
               "Octβ3R",
               "PK1-R",
               'PK2-R2',
               "Pdfr",
               "RYa-R",
               "SIFaR",
               "SPR",
               'TkR86C',
               "TyrRII",
               "hec",
               "mAChR-A",	
               "mAChR-B",	
               'mGluR',	
               "mtt",	
               "nAChR%ce%b11", # "nAChRα1",	# nAChR%ce%b11
               "nAChR%ce%b12", # "nAChRα2",	#  nAChR%ce%b12
               "nAChR%ce%b13", # "nAChRα3",	# nAChR%ce%b13
               "nAChR%ce%b16",  # "nAChRα6",	# nAChR%ce%b16
               "nAChR%ce%b17", # "nAChRα7", # nAChR%ce%b17
               "nAChR%ce%b21", # "nAChRβ1", # nAChR%ce%b21
               "nAChR%ce%b22", # "nAChRβ2", # nAChR%ce%b22
               "nAChR%ce%b23", # "nAChRβ3",	# nAChR%ce%b23
               "ort",	
               "rk")
all.receptors.simple <- c("5-HT1B",
                   "5-HT2A",
                   "5-HT2B",
                   "5-HT7",
                   "AkhR",
                   "AstA-R1",
                   'AstA-R2',
                   "AstC-R1",
                   "AstC-R2",
                   'CCAP-R',
                   "CCHa1-R",
                   'CCHa2-R',
                   "CCKLR-17D1",
                   "CCKLR-17D3",
                   'CG12344',
                   'CG13229',
                   "CG13575",
                   "CG13995",
                   "CG32547",
                   "33639",
                   'CG43795',
                   "CG7589",
                   "CapaR",
                   "CrzR",
                   "Dh31-R",
                   "Dh44-R2",
                   "DopR1",
                   'DopR2',
                   "D2R",
                   "DopEcR",
                   'GABA-B-R1',
                   "GABA-B-R2",
                   "GluCla",
                   "GluRIB",
                   "GluRIIA",
                   "GluRIIB",
                   "GluRIIC",
                   "GluRIID",
                   'GluRIIE',
                   "Grd",
                   'Gyc76C',
                   "KaiR1D",
                   "Lcch3",
                   "Lgr1",
                   "Lgr3",
                   "Lgr4",
                   "NPFR",
                   "Nmdar1",
                   "Nmdar2",
                   "Oamb",
                   "C1-OctB2R",
                   "OctB3R",
                   "PK1-R",
                   'PK2-R2',
                   "Pdfr",
                   "RYa-R",
                   "SIFaR",
                   "SPR",
                   'TkR86C',
                   "TyrRII",
                   "hec",
                   "mAChR-A",	
                   "mAChR-B",	
                   'mGluR',	
                   "mtt",	
                   "nAChRa1",	
                   "nAChRa2",	
                   "nAChRa3",	
                   "nAChRa6",	
                   "nAChRa7",
                   "nAChRb1",
                   "nAChRb2",
                   "nAChRb3",	
                   "ort",	
                   "rk")
names(all.receptors) <- all.receptors.simple
receptorfolder <- "data/kondo_et_al_2020/nrrd/"
sources <- c("raw") # "src",
templates <- c("IS2","","templete","template")
receptorfolder.files <- dir(receptorfolder)
receptors.done <- sort(unique(gsub("IS2_|templete_|template_|_no.*","",receptorfolder.files)))
receptors.missing <- all.receptors[!all.receptors%in%receptors.done]
receptors.missing <- receptors.missing[!replace_greek_with_latin(receptors.missing)%in%receptors.done]
receptors.missing <- receptors.missing[!tolower(receptors.missing)%in%receptors.done]
receptors.missing <- receptors.missing[!tolower(replace_greek_with_latin(receptors.missing))%in%receptors.done]
receptors.missing <- receptors.missing[!names(receptors.missing)%in%receptors.done]
for(src in sources){
  for(template in templates){
    for(r in 1:length(receptors.missing)){
      receptor <- receptors.missing[r]
      for(num in 1:2){
        for(brain in 1:2){
          is.greek <- contains_greek(receptor)
          receptors <- c(receptor,names(receptor),tolower(receptor))
          if(is.greek){
            receptor3 <- replace_greek_with_latin(receptor)
            receptors <- c(receptors,receptor3)
          }
          for(receptor2 in receptors){
            try({
              # raw or src?
              url<-sprintf("https://gin.g-node.org/Tanimoto_lab/T2A-GAL4_Library/%s/master/%s/UAS-mCD8-GFP_brp-SNAP_IS2_Registration/%s_%s_no%s_0%s_warp_m0g40c4e1e-1x16r3.nrrd",
                           src, receptor, template, receptor2, brain, num)
              if(!RCurl::url.exists(url)){
                url<-sprintf("https://gin.g-node.org/Tanimoto_lab/T2A-GAL4_Library/%s/master/%s/UAS-mCD8-GFP_brp-SNAP_IS2_Registration/%s_%s_no%s_m2_0%s_warp_m0g40c4e1e-1x16r3.nrrd",
                             src, receptor, template, receptor2, brain, num)
                if(!RCurl::url.exists(url)){
                  message("URL not found: ", url)
                  next
                }
              }
              nam <- basename(url)
              file <- file.path(receptorfolder,nam)
              if(file.exists(file)){
                message("already downloaded: ", file)
                next
              }
              download.file(url = url, destfile = file, quiet = FALSE) 
              message("written: ", file)
            })  
          }
        }
      }
    }
  }
}
options(timeout = min(60, getOption("timeout")))

# Get the flywire neurons with the most DCVS
# Save these neurons as tiff files.

# Read tangential neurons
fw.meta <- readr::read_csv("/Users/abates/projects/wilson-lab/nat-tech/data/flywire/flywire_meta.csv", col_types = hemibrainr:::sql_col_types)
fw.tangential.meta <- subset(fw.meta, cell_class == "CX") # cell_sub_class!="tangential" & 
cts <- sort(unique(fw.tangential.meta$cell_type))

# Choose save folders
receptorfolder <- "data/kondo_et_al_2020/nrrd/"
flywirenrrdfolder <- "neuroanat/flywire_fb_nrrd"
flywiremaxfolder <- "neuroanat/flywire_fb_mp"
dir.create(flywirenrrdfolder, recursive = TRUE)
dir.create(flywiremaxfolder, recursive = TRUE)

# Cycle through tangential cell types and save as .tiff files.
for(ct in cts){
  
  # get IDs
  ids <- unique(subset(fw.tangential.meta, cell_type==ct)$root_id)
  nt <- unique(subset(fw.tangential.meta, cell_type==ct)$top_nt)[1]
  hb <- unique(subset(fw.tangential.meta, cell_type==ct)$hemibrain_type)[1]
  
  # Save flywire .nrrd files
  flywireid_to_nrrd(flywire_id = ids, 
                    cell_type = paste(ct, hb, nt, collapse = "_", sep = "_"), 
                    ref="IS2", 
                    savefolder = flywirenrrdfolder,
                    savemaxfolder = flywiremaxfolder,
                    plot3D = FALSE, 
                    compressed = FALSE,
                    max_projection = TRUE,
                    overwrite = FALSE)
}

# Same with hemibrain
# Read tangential neurons
hb.meta <- readr::read_csv("/Users/abates/projects/wilson-lab/nat-tech/data/hemibrain/Supplemental_file4_hemibrain_meta.csv", col_types = hemibrainr:::sql_col_types)
hb.tangential.meta <- subset(hb.meta, grepl("^FB|^FC|^hDelta|^vDelta|^PFL|^FS|^SA1|^SA2|^SA3",type))
hb.tangential.meta$cell_type <- gsub("_.*","",hb.tangential.meta$type)
hb.cts <- sort(unique(hb.tangential.meta$cell_type))

# Choose save folders
hemibrainnrrdfolder <- "neuroanat/hemibrain_fb_nrrd"
hemibrainmaxfolder <- "neuroanat/hemibrain_fb_mp"
dir.create(hemibrainnrrdfolder, recursive = TRUE)
dir.create(hemibrainmaxfolder, recursive = TRUE)

# Cycle through tangential cell types and save as .tiff files.
for(hct in hb.cts){
  
  # get IDs
  ids <- unique(subset(fw.tangential.meta, cell_type==ct)$root_id)
  nt <- unique(subset(fw.tangential.meta, cell_type==ct)$top_nt)[1]
  
  # Save flywire .nrrd files
  hemibrain_to_nrrd(cell_type = hct, 
                    savemaxfolder = hemibrainmaxfolder,
                    ref="IS2", 
                    savefolder = hemibrainnrrdfolder,
                    plot3D = FALSE, 
                    mesh = TRUE,
                    max_projection = TRUE,
                    overwrite = FALSE)
  
}

# Function to normalize an image
normalize_image <- function(img, 
                            thresh = FALSE, 
                            lower = 0.01, 
                            upper = 0.99) {
  if (thresh){
    img[img<quantile(img,lower)] <- quantile(img,lower)
    img[img>quantile(img,upper)] <- quantile(img,upper) 
  }
  min_val <- min(img)
  max_val <- max(img)
  normalized_img <- (img - min_val) / (max_val - min_val)
  return(normalized_img)
}

# Read and combine maximal projection .tiff file
receptorfolder <- "data/kondo_et_al_2020/nrrd/"
resultfolder <- "neuroanat/flywire_receptor_mps"
maxfolder.files <- rev(list.files(c(hemibrainmaxfolder), full.names = TRUE)) # flywiremaxfolder
receptorfolder.files <- list.files(receptorfolder, full.names = TRUE)
receptorfolder.files <- receptorfolder.files[grepl("02_warp",receptorfolder.files)]
only.fb <- FALSE
for(f in 1:length(maxfolder.files)){
  
  # Read neuron file
  tiff1 <- maxfolder.files[f]
  img1 <- readTIFF(tiff1)
  if(only.fb){
    img1 <- img1[100:214,243:519]
  }
  normalized_img1 <- normalize_image(img1, thresh = TRUE, lower = 0)
  
  # Flywire or hemibrain file?
  if(grepl("flywire",tiff1)){
    resultfolder <- "neuroanat/flywire_receptor_mps"
    dataset <- "flywire"
  }else{
    resultfolder <- "neuroanat/hemibrain_receptor_mps"
    dataset <- "hemibrain"
  }
  dir.create(resultfolder, showWarnings = FALSE)
  
  # Iterate over receptor files
  for(ff in 1:length(receptorfolder.files)){
    
    # buffer against errors
    # try({
      # overwrite?
      nrrd2 <- receptorfolder.files[ff]
      celltype <- gsub(".*max_projection_(.*)_IS2.*", "\\1", tiff1)
      celltype <- paste(unlist(strsplit(celltype,split = " ")),collapse="_")
      recept <- gsub("\\_warp.*", "", basename(nrrd2))
      if(only.fb){
        receptor.folder <- file.path(resultfolder,paste0("fb_",recept))
      }else{
        receptor.folder <- file.path(resultfolder,recept)
      }
      dir.create(receptor.folder, showWarnings = FALSE)
      save.file <- file.path(receptor.folder, paste0(dataset,"_",celltype,"_",recept,".tiff"))
      if(file.exists(save.file) && (f!=1 && f!=1)){
        message("File already exists: ", save.file)
        next
      }
      
      # Read receptor expression image
      m <-  nat::read.nrrd(nrrd2)
      
      # 74:133 planes in IS2 capture the FB
      if (only.fb){
        img2 <- t(apply(m[100:214,72:386,86:130],c(1,2),max))
      }else{
        img2 <- t(apply(m[,,86:130],c(1,2),max))
      }
      
      # Normalize the TIFF images
      normalized_img2 <- normalize_image(img2, thresh = TRUE)
      image_array <- array(0, dim = c(nrow(img2), ncol(img2), 3))
      image_array[,,1] <- normalized_img2   # Red channel
      image_array[,,2] <- normalized_img1 # Green channel
      image_array[,,3] <- normalized_img2 # Blue channel
      
      # Write
      writeTIFF(image_array, save.file)
      message("Wrote: ", save.file)
      # })
      
      # save max proj of FB
      if(f==1 && f==1){
        save.file2 <- file.path(receptor.folder, paste0(recept,".tiff"))
        writeTIFF(normalized_img2, save.file2)
        message("Wrote: ", save.file2)
      }
      
  }
}
















