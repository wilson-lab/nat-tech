library(tidyverse)
library(edgeR)
library(ggplot2)
library(Matrix)
library(cowplot)
library(Seurat)
library(plyr)
library(SingleCellExperiment)
library(matrixStats)
library(umap)
library(foreach)
library(DoubletFinder)
library(ggplot2)

# Glutamate receptor associates: https://string-db.org/network/7227.FBpp0076691
# GluRs:
glurs <- c("stg1","PICK1","CASK","Neto", # GluR correlates
          "GluRIA", "GluRIB", # AMPAR brain
          "GluRIIC", "GluRIID", "GluRIIE", # Muscle specific 
          "GluRIIA", "GluRIIB", # Muscle specific 
          "mGluR", # metabotropic 
          "GluClalpha", # GluCLalpha
          "clumsy", "CG11155", # Kainate
          "CG3822", # DKaiR1D, Kainate
          "CG5621", # DKaiR1C, Kainate
          "Nmdar1", "Nmdar2" # NMDA
          )
typeII <- c("Imp")
dl1 <- c("gcm")
dm <- c("Rx")

##### Croset et al. 2018 #####
# Read-in expression data
files <- list.files("/Users/abates/projects/wilson-lab/nat-tech/data/croset_et_al_2018/", full.names = TRUE)
df <- edgeR::readDGE(files)

