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

##### Croset et al. 2018 #####
# # Read-in expression data
# files <- list.files("/Users/abates/projects/wilson-lab/nat-tech/data/croset_et_al_2018/", full.names = TRUE)
# df <- edgeR::readDGE(files)

##### Step 1 #####
# Read in the data
# Data downloaded here: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE207799
sat00_rep1.raw <- Read10X_h5("data/park_et_al_2022/GSM6318432_2019JUN_water_G4_rep1_filtered_feature_bc_matrix.h5")
sat00_rep2.raw <- Read10X_h5("data/park_et_al_2022/GSM6318433_2019JUN_water_G4_rep2_filtered_feature_bc_matrix.h5")
sat00_rep1 <- CreateSeuratObject(counts = sat00_rep1.raw, project = "sat00_rep1")
sat00_rep2 <- CreateSeuratObject(counts = sat00_rep2.raw, project = "sat00_rep2")

##### Step 2 #####
# Calculate and plot percent of mitochondrial and ribosomal RNA, and ribosomal proteins.
# Cells with high level of mitochondrial or ribosomal RNA, 
# or ribosomal proteins are probably unhealthy, stressed, or dying.
sat00_rep1[["percent.mt"]] <- PercentageFeatureSet(sat00_rep1, pattern = "^mt:")
sat00_rep2[["percent.mt"]] <- PercentageFeatureSet(sat00_rep2, pattern = "^mt:")
sat00_rep1[["percent.rRNA"]] <- PercentageFeatureSet(sat00_rep1, pattern = "(S|-)rRNA")
sat00_rep2[["percent.rRNA"]] <- PercentageFeatureSet(sat00_rep2, pattern = "(S|-)rRNA")
sat00_rep1[["percent.rProt"]] <- PercentageFeatureSet(sat00_rep1, pattern = "Rp(L|S)")
sat00_rep2[["percent.rProt"]] <- PercentageFeatureSet(sat00_rep2, pattern = "Rp(L|S)")
seurat_list <- list(sat00_rep1, sat00_rep2)
mt_table <- list()
rr_table <- list()
rp_table <- list()
hsp.table <- list()
for(i in seq_along(seurat_list)){
  mt_table[[seurat_list[[i]]@project.name]] <- as.numeric(seurat_list[[i]]$percent.mt)
  rr_table[[seurat_list[[i]]@project.name]] <- as.numeric(seurat_list[[i]]$percent.rRNA)
  rp_table[[seurat_list[[i]]@project.name]] <- as.numeric(seurat_list[[i]]$percent.rProt)
}
mt_table <- data.frame(lapply(mt_table, `length<-`, max(lengths(mt_table))))
rr_table <- data.frame(lapply(rr_table, `length<-`, max(lengths(rr_table))))
rp_table <- data.frame(lapply(rp_table, `length<-`, max(lengths(rp_table))))
perc.remove.mito <- 15
perc.remove.rrna <- 10
perc.remove.rprot <- 15

# Plot percentage mitochondrial DNA
col_vector <- scale_color_manual(values=rep(c("#4e85c5", "#7e93a2", "#bb9e69", "#b373a6"), each=2))
mt_table %>% 
  gather(sample, perc_mito, na.rm=TRUE, factor_key=TRUE) %>% 
  ggplot(aes(x=sample, y=perc_mito)) + labs(x="Samples", y="Percentage of mitochondrial RNA") + 
  geom_jitter(size=.1, aes(colour=sample)) + 
  theme(axis.text.x = element_text(angle=45, hjust=1), legend.position = "none") + 
  geom_boxplot(outlier.size=-1, alpha=0.5) + geom_hline(yintercept=perc.remove.mito, linetype="dashed") + 
  col_vector

##### Step 3 #####
# Calculate and plot number of features and number of UMIs.
# Barcodes with too few features (genes) or UMIs don’t contain enough information. 
# Barcodes with too many UMIs or genes may be doublets. These barcodes should be removed.
num_table <- list()
umi_table <- list()
for(i in seq_along(seurat_list)){
  num_table[[seurat_list[[i]]@project.name]] <- as.numeric(seurat_list[[i]]$nFeature_RNA)
  umi_table[[seurat_list[[i]]@project.name]] <- as.numeric(seurat_list[[i]]$nCount_RNA)
}
num_table <- data.frame(lapply(num_table, `length<-`, max(lengths(num_table)))) %>% as_tibble %>% gather(sample, n_features, na.rm=TRUE, factor_key=TRUE)
umi_table <- data.frame(lapply(umi_table, `length<-`, max(lengths(umi_table)))) %>% as_tibble %>% gather(sample, n_umis, na.rm=TRUE, factor_key=TRUE)
min_features <- 300
max_features <- 4500
max_UMIs <- 25000

##### Step 4 #####
# Trim data down
trim_data <- function(object){
  data_trim <- object %>% subset(subset = nFeature_RNA >= min_features & nFeature_RNA <= max_features & percent.mt < perc.remove.mito & percent.rRNA < perc.remove.rrna & nCount_RNA <= max_UMIs & percent.rProt <= perc.remove.rprot)
  removed <- dim(object)[2] - dim(data_trim)[2]
  perc_removed <- removed/dim(object)[2]*100
  print(paste0("Number of cells removed in ", object@project.name, ": ", removed, " (", round(perc_removed,2), "%)"))
  return(data_trim)
}
sat00_rep1_trim <- sat00_rep1 %>% 
  trim_data()
sat00_rep2_trim <- sat00_rep2 %>% 
  trim_data()
thirst_list <- list(sat00_rep1_trim, sat00_rep2_trim)

# We normalise the data and add the expression level of roX1 to the metadata 
stim <- rep("sated", each=2)
for(i in 1:length(thirst_list)){
  thirst_list[[i]]$stim <- stim[i]
  thirst_list[[i]] <- NormalizeData(thirst_list[[i]])
  thirst_list[[i]] <- AddModuleScore(thirst_list[[i]], "lncRNA:roX1", name="sex", seed=123)
}
sat00 <- merge(thirst_list[[1]], thirst_list[[2]], add.cell.ids = c("sat00_rep1", "sat00_rep2"))
stim_list <- list(SCTransform(sat00, vars.to.regress=c("orig.ident", "sex1"), verbose=FALSE)) 

# Use the ‘anchors’ method from Seurat V3 to integrate all samples together
IntegrationPCs <- 30
thirst_features <- SelectIntegrationFeatures(object.list=stim_list, nfeatures=2000)
options(future.globals.maxSize=8*1024*1024^2)
stim_list <- PrepSCTIntegration(object.list=stim_list, anchor.features=thirst_features, verbose=FALSE)
thirst_anchors <- FindIntegrationAnchors(object.list=stim_list, dims = 1:IntegrationPCs, scale=F, normalization.method="SCT", reduction = "cca", anchor.features=thirst_features, verbose=FALSE)
Thirst2_SCT <- IntegrateData(anchorset=thirst_anchors, normalization.method="SCT", verbose=FALSE)
Idents(Thirst2_SCT) <- Thirst2_SCT$seurat_clusters <- Thirst2_SCT$integrated_snn_res.2
DimPlot(Thirst2_SCT, label=T, pt.size=.2) + NoAxes() + theme(legend.position="none")

# Seurat pipeline
# Scale
object_list <- list(sat00)
object_list[[1]] <- ScaleData(object_list[[1]], vars.to.regress = c("orig.ident", "nCount_RNA", "sex1"))
# Integrate
Thirst2.anchors <- FindIntegrationAnchors(object.list=object_list, dims = 1:IntegrationPCs, scale = F, reduction = "cca", anchor.features = 2000)
Thirst2 <- IntegrateData(anchorset = Thirst2.anchors, dims = 1:IntegrationPCs)
# Cluster
DefaultAssay(Thirst2) <- "integrated"
Thirst2 <- ScaleData(Thirst2, vars.to.regress="sex1")
Thirst2 <- RunPCA(Thirst2, npcs = 20)
Thirst2 <- RunUMAP(Thirst2, reduction="pca", dims=1:20, n.neighbors=30)
Thirst2 <- FindNeighbors(Thirst2, reduction="pca", dims=1:20, k.param=30, force.recalc=TRUE)
Thirst2 <- FindClusters(Thirst2, resolution=2, graph.name="integrated_snn")

##### Step 6 #####
# Remove doublets
seurat_list <- SplitObject(Thirst2, split.by="orig.ident")
seurat_list <- foreach(s=seurat_list) %do% {
  s <- NormalizeData(s)
  s <- ScaleData(s, vars.to.regress=c("nCount_RNA", "sex1"))
  s <- FindVariableFeatures(s)
  RunPCA(s, npcs=PCs, verbose=F)
}
bcmvn_tables <- foreach(s=seurat_list) %do% {
  sweep.res.list <- paramSweep_v3(s, PCs = 1:20)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- as.tibble(sweep.stats) %>% group_by(pK) %>% dplyr::summarize(MeanBC=mean(BCreal), VarBC=sd(BCreal)^2, BCmetric=mean(BCreal)/(sd(BCreal)^2))
}
Thirst2$DF_classification <- "Singlet"
foreach(i=unique(Thirst2$orig.ident), b=bcmvn_tables) %do% {
  object.i <- subset(Thirst2, subset=orig.ident==i)
  object.i <- NormalizeData(object.i)
  object.i <- FindVariableFeatures(object.i)
  object.i <- ScaleData(object.i, vars.to.regress=c("nCount_RNA", "sex1"))
  object.i <- RunPCA(object.i, npcs=PCs, verbose=F)
  homotypic.prop <- modelHomotypic(Idents(object.i))
  nExp_poi <- round(0.05*dim(object.i)[2])
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  pK.i <- b %>% filter(BCmetric==max(BCmetric)) %>% pull(pK) %>% as.character() %>% as.double()
  object.i <- doubletFinder_v3(object.i, PCs = 1:20, pN = 0.25, pK = pK.i, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
  doublet_found <- object.i@meta.data[,grepl("DF.classifications", names(object.i@meta.data))]
  Thirst2$DF_classification[Thirst2$orig.ident==i] <- doublet_found
}
Thirst2@misc$bcmvn_tables <- bcmvn_tables

# Filter down by Dale's law
pairs <- list(c("VAChT", "nrv2"), c("Gad1", "nrv2"), c("VGlut", "nrv2"), c("VAChT", "Gad1"), c("VAChT", "VGlut"), c("Gad1", "VGlut"))
lims <- data.frame(Gene=unlist(pairs), Thresh=c(1.5,3,2.3,3,2.1,3,1.5,2.3,1.5,2.1,2.3,2.1)) %>% mutate(Gene=factor(Gene, levels=unique(unlist(pairs))))
Thirst2_SCT$doublet_ID <- "Singlet"
tick=1
for(j in pairs){
  lim.j <- lims[c(tick,tick+1),]
  doublets_toflag <- names(Thirst2_SCT[,Thirst2_SCT@assays$RNA@data[j[1],] > lim.j[lim.j$Gene==j[1],]$Thresh & Thirst2_SCT@assays$RNA@data[j[2],] > lim.j[lim.j$Gene==j[2],]$Thresh]$nCount_RNA)
  Thirst2_SCT$doublet_ID[doublets_toflag][Thirst2_SCT$doublet_ID[doublets_toflag]!="Singlet"] <- paste(paste(j, collapse=" & "), Thirst2_SCT$doublet_ID[doublets_toflag][Thirst2_SCT$doublet_ID[doublets_toflag]!="Singlet"], sep=" + ")
  Thirst2_SCT$doublet_ID[doublets_toflag][Thirst2_SCT$doublet_ID[doublets_toflag]=="Singlet"] <- paste(j, collapse=" & ")
  tick=tick+2
}

# And by KC identity
kcmarkers <- AverageExpression(Thirst2_SCT, features=c("ey", "Dop1R2", "mub"), assays="RNA")
Thirst2_KC_SCT <- subset(Thirst2_SCT, idents=kcmarkers$RNA %>% t() %>% as.data.frame() %>% rownames_to_column("cluster") %>% as_tibble() %>% filter(ey>2 & Dop1R2>2 & mub>50) %>% pull(cluster))
DefaultAssay(Thirst2_KC_SCT) <- "RNA"
Thirst2_KC_SCT <- SCTransform(Thirst2_KC_SCT, vars.to.regress=c("sex1", "orig.ident"), verbose=FALSE)
Thirst2_KC_SCT <- RunPCA(Thirst2_KC_SCT, npcs=PCs)
Thirst2_KC_SCT <- FindNeighbors(Thirst2_KC_SCT, reduction = "pca", dims = 1:10)
Thirst2_KC_SCT <- FindClusters(object = Thirst2_KC_SCT, resolution = 1, save.ssn = T)
Thirst2_KC_SCT <- RunUMAP(Thirst2_KC_SCT, reduction = "pca", dims = 1:10)
pairsKC <- list(c("Ca-alpha1T", "ab"), c("Ca-alpha1T", "CG8641"), c("ab", "CG8641"))
limsKC <- data.frame(Gene=unlist(pairsKC), Thresh=c(1.5,1.5,1.5,2.2,1.5,2.2)) %>% mutate(Gene=factor(Gene, levels=unique(unlist(pairsKC))))
Thirst2_KC_SCT$KC_doublets <- "Singlet"
tick=1
for(j in pairsKC){
  lim.j <- limsKC[c(tick,tick+1),]
  doublets_toflag <- names(Thirst2_KC_SCT[,Thirst2_KC_SCT@assays$RNA@data[j[1],] > lim.j[lim.j$Gene==j[1],]$Thresh & Thirst2_KC_SCT@assays$RNA@data[j[2],] > lim.j[lim.j$Gene==j[2],]$Thresh]$nCount_RNA)
  Thirst2_KC_SCT$KC_doublets[doublets_toflag][Thirst2_KC_SCT$KC_doublets[doublets_toflag]!="Singlet"] <- paste(paste(j, collapse=" & "), Thirst2_KC_SCT$KC_doublets[doublets_toflag][Thirst2_KC_SCT$KC_doublets[doublets_toflag]!="Singlet"], sep=" + ")
  Thirst2_KC_SCT$KC_doublets[doublets_toflag][Thirst2_KC_SCT$KC_doublets[doublets_toflag]=="Singlet"] <- paste(j, collapse=" & ")
  tick=tick+2
}
Thirst2_KC_SCT$KC_doublets[WhichCells(Thirst2_KC_SCT, idents=c(14,15,17))] <- "Mixed_cluster"
Thirst2_SCT$KC_doublets <- "Singlet"
Thirst2_SCT$KC_doublets[colnames(Thirst2_SCT) %in% colnames(Thirst2_KC_SCT)] <- Thirst2_KC_SCT$KC_doublets

# visualise and remove
doublet_table <- as_tibble(Thirst2_SCT@reductions$umap@cell.embeddings) %>% mutate(Doublets=paste(Thirst2_SCT$doublet_ID, Thirst2_SCT$DF_classification_original, Thirst2_SCT$KC_doublets, sep=" + "), Cluster=Idents(Thirst2_SCT), orig.ident=as.factor(Thirst2_SCT$orig.ident)) %>% mutate(Doublets=str_replace(Doublets, "(Singlet \\+ )|( \\+ Singlet)", ""), Doublets=str_replace(Doublets, "(Singlet \\+ )|( \\+ Singlet)", "")) %>% mutate(Doublets=str_replace(Doublets, "Doublet", "DoubletFinder")) %>% mutate(Doublets=factor(Doublets, levels=unique(Doublets)))
doublet_table %>% mutate(Doublets=case_when(grepl("(+ DoubletFinder)|(DoubletFinder +)", Doublets) ~ "Coexpression +\nDoubletFinder", grepl("&", Doublets)|Doublets=="Mixed_cluster" ~ "Coexpression", TRUE ~ as.character(Doublets))) %>% mutate(Doublets=factor(Doublets, levels=c("Singlet", "DoubletFinder", "Coexpression", "Coexpression +\nDoubletFinder"))) %>% arrange(Doublets) %>% ggplot(aes(x=UMAP_1, y=UMAP_2, color=Doublets)) + geom_point(size=.2) + theme_void() + theme(legend.title=element_blank()) + guides(colour=guide_legend(override.aes=list(size=3), ncol=1)) + scale_color_manual(values=c("lightsteelblue2", "orange", "purple", "black"))
Thirst2_SCT_trimPlus <- Thirst2_SCT[, doublet_table$Doublets=="Singlet" & Thirst2_SCT$doublet_ratio<.5]
cat("Number of doublets removed:", dim(Thirst2_SCT)[2]-dim(Thirst2_SCT_trimPlus)[2])

##### Step 7 #####
# Assign major cel lclass groups
# PCA, Clustering & UMAP
DefaultAssay(Thirst2_SCT_trimPlus) <- "integrated"
Thirst2_SCT_trimPlus <- RunPCA(Thirst2_SCT_trimPlus, npcs = 20)
Thirst2_SCT_trimPlus <- FindNeighbors(Thirst2_SCT_trimPlus, reduction="pca", dims=1:20, k.param=30)
Thirst2_SCT_trimPlus <- FindClusters(Thirst2_SCT_trimPlus, resolution=2, graph.name="integrated_snn")
Thirst2_SCT_trimPlus <- RunUMAP(Thirst2_SCT_trimPlus, reduction="pca", dims=1:20, n.neighbors=30)
markers <- c("nSyb", "elav", "VAChT", "VGlut", "Gad1", "ey", "Dop1R2", "Pka-C1", "Vmat", "Ilp2", "Ilp3", "CG10433")
markers_avex <- AverageExpression(Thirst2_SCT_trimPlus, assays="RNA", features=markers)$RNA %>% rownames_to_column("Gene") %>% gather(Clusters, Avex, -Gene)
thresholds <- data.frame(Gene=factor(markers, levels=markers), thresh=c(8,2.5,2,12,10,3,6.5,30,40,500,250,5))
markers_avex %>% arrange(Gene, Avex) %>% unite("gene_clusters", Gene, Clusters, sep = "_", remove = FALSE) %>% mutate(gene_clusters=factor(gene_clusters, levels=gene_clusters), Gene=factor(Gene, levels=markers)) %>% ggplot(aes(x=gene_clusters,y=Avex)) + geom_point(size=.2) + facet_wrap(~Gene, scales="free", nrow=3) + theme(strip.text=element_text(face="italic"), axis.ticks.x=element_blank(), axis.text.x = element_blank()) + labs(x="Clusters", y="Average expression") + coord_cartesian(clip = "off") + geom_hline(data=thresholds, aes(yintercept=thresh), color="darkred")
thresholds <- thresholds %>% spread(Gene, thresh)
cellgroups <- list()
markers_values <- markers_avex %>% spread(Gene, Avex)
cellgroups$neurons <- markers_values %>% filter(nSyb>=thresholds$nSyb & elav>=thresholds$elav) %>% pull(Clusters)
cellgroups$cholinergic <- markers_values %>% filter(VAChT>=thresholds$VAChT) %>% pull(Clusters)
cellgroups$glutamatergic <- markers_values %>% filter(VGlut>=thresholds$VGlut) %>% pull(Clusters)
cellgroups$gabaergic <- markers_values %>% filter(Gad1>=thresholds$Gad1) %>% pull(Clusters)
cellgroups$kenyon <- markers_values %>% filter(ey>=thresholds$ey & Dop1R2>=thresholds$Dop1R2 & `Pka-C1`>=thresholds$`Pka-C1`) %>% pull(Clusters)
cellgroups$monoaminergic <- markers_values %>% filter(Vmat>=thresholds$Vmat) %>% pull(Clusters)
cellgroups$IPCs <- markers_values %>% filter(Ilp2>=thresholds$Ilp2 & Ilp3>=thresholds$Ilp3) %>% pull(Clusters)
cellgroups$glia <- markers_values %>% filter(CG10433>=thresholds$CG10433) %>% pull(Clusters)
ann_clusters <- as.character(Idents(Thirst2_SCT_trimPlus))
CT <- rep("Other", length(ann_clusters))
for(i in unique(ann_clusters)){
  if(i %in% cellgroups$glia & !i %in% cellgroups$neurons){
    CT[ann_clusters==i] <- "Glia"
  } else if(i %in% cellgroups$neurons & !i %in% cellgroups$glia){
    if(i %in% cellgroups$cholinergic & !i %in% c(cellgroups$glutamatergic, cellgroups$gabaergic, cellgroups$monoaminergic, cellgroups$IPCs)){
      if(i %in% cellgroups$kenyon){
        CT[ann_clusters==i] <- "KCs"
      } else {
        CT[ann_clusters==i] <- "ACh"
      }
    } else if(i %in% cellgroups$glutamatergic & !i %in% c(cellgroups$cholinergic, cellgroups$gabaergic, cellgroups$monoaminergic, cellgroups$IPCs)){
      CT[ann_clusters==i] <- "Glut"
    } else if(i %in% cellgroups$gabaergic & !i %in% c(cellgroups$cholinergic, cellgroups$glutamatergic, cellgroups$monoaminergic, cellgroups$IPCs)){
      CT[ann_clusters==i] <- "GABA"
    } else if(i %in% cellgroups$monoaminergic & !i %in% c(cellgroups$cholinergic, cellgroups$gabaergic, cellgroups$glutamatergic, cellgroups$IPCs)){
      CT[ann_clusters==i] <- "Monoamines"
    } 
  }
}
Thirst2_SCT_trimPlus$celltype <- CT

# Re-cluster major cell groups
celltypes <- SplitObject(Thirst2_SCT_trimPlus, split.by="celltype")[c(5,1,3,4,7,2,6)]
PCs <- c(19,17,16,12,16,14,12)
resolutions <- c(1,1,1,1,4,1,2)
Thirst2_celltypes <- foreach(c=celltypes, n=names(celltypes), p=PCs, r=resolutions) %do% {
  c@project.name <- n
  c <- SCTransform(c, vars.to.regress=c("orig.ident", "sex1"), verbose=FALSE)
  c <- RunPCA(c, npcs = p)
  c <- FindNeighbors(c, reduction = "pca", dims = 1:p)
  c <- FindClusters(c, resolution = r, save.ssn=T)
  c <- RunUMAP(c, reduction = "pca", dims = 1:p)
}
Thirst2_celltypes <- foreach(c=Thirst2_celltypes, p=PCs) %do% {
  
  # Build and process tree
  c <- BuildClusterTree(c, dims=1:p, assay="integrated")
  phylo_table <- as_tibble(c@tools$BuildClusterTree$edge) %>% mutate(lengths=c@tools$BuildClusterTree$edge.length, Cluster=V2-1)
  
  # Merge clusters
  c$merged_clusters <- Idents(c)
  nclust_orig <- length(levels(c))
  lim_n_sig <- 10
  full <- FALSE
  clear <- c()
  round <- 1
  while(full==FALSE){
    full <- TRUE
    i <- 1
    while(i<nrow(phylo_table)){
      if(phylo_table[i,]$V2<=nclust_orig & phylo_table[i+1,]$V2<=nclust_orig & phylo_table[i,]$V1==phylo_table[i+1,]$V1 & !phylo_table[i,]$Cluster %in% clear){
        markers_i <- FindMarkers(c, phylo_table[i,]$Cluster, phylo_table[i+1,]$Cluster, assay="RNA", min.pct=.25)
        n_sig <- markers_i %>% rownames_to_column("Gene") %>% dplyr::filter(abs(avg_logFC)>=1, p_val_adj<.05, !grepl("(lncRNA|tRNA|^mt:|rRNA)", Gene)) %>% dplyr::summarise(n()) %>% pull()
        if(n_sig<lim_n_sig){
          min_clus <- min(c(phylo_table[i,]$Cluster, phylo_table[i+1,]$Cluster))
          other_clus <- max(c(phylo_table[i,]$Cluster, phylo_table[i+1,]$Cluster))
          c$merged_clusters[c$merged_clusters==other_clus] <- min_clus
          phylo_table[i-1,] <- t(c(phylo_table[i-1,]$V1, min_clus+1, sum(phylo_table[c(i,i-1),]$lengths), min_clus))
          phylo_table <- phylo_table[-c(i,i+1),]
          full <- FALSE
          i=i-1
        } else{
          clear <- c(clear, phylo_table[i,]$Cluster, phylo_table[i+1,]$Cluster)
        }
      }
      i=i+1
    }
    round=round+1
  }
  c@tools$phylo_table <- phylo_table
  c
}
Thirst2_celltypes <- foreach(c=Thirst2_celltypes) %do% {
  Idents(c) <- c$merged_clusters
  c
}

# Save
save(Thirst2_SCT_trimPlus, file="Thirst2_SCT_trimPlus.Robj")
save(Thirst2_celltypes, file="Thirst2_celltypes.Robj")

######################
##### GluR plots #####
######################

# CX marker: retinal homeobox (rx) genetic neural lineage, earmuff (ermR09D11), Deadpan 
# DL1 marker: gcm?
# Primary marker: Imp

GluR_markers <- c("ct", "acj6", "Lim1")
GluR_thresholds <- data.frame(Gene=GluR_markers, thresh=c(7,5,7))
for(c in Thirst2_celltypes){if(c@project.name=="GluR"){GluR_object <- c; Idents(GluR_object) <- GluR_object$merged_clusters}}
GluR_markers_table <- as.data.frame(GluR_object@assays$RNA@data[GluR_markers,]) %>% rownames_to_column("Gene") %>% gather(Cells, Values, -Gene) %>% mutate(Gene=factor(Gene, levels=unique(Gene)))
GluR_markers_avex <- AverageExpression(GluR_object, assays="RNA", features=GluR_markers)$RNA %>% rownames_to_column("Gene") %>% gather(Clusters, Avex, -Gene)
GluR_markers_avex %>% arrange(Gene, Avex) %>% unite("gene_clusters", Gene, Clusters, sep = "_", remove = FALSE) %>% mutate(gene_clusters=factor(gene_clusters, levels=gene_clusters), Gene=factor(Gene, levels=GluR_markers)) %>% ggplot(aes(x=gene_clusters,y=Avex)) + geom_point(size=.2) + facet_wrap(~Gene, scales="free") + theme(strip.text=element_text(face="italic"), axis.ticks.x=element_blank(), axis.text.x = element_blank()) + labs(x="Clusters", y="Average expression") + coord_cartesian(clip = "off") + geom_hline(data=GluR_thresholds, aes(yintercept=thresh), color="darkred")

Idents(GluR_object) <- GluR_object$merged_clusters
GluR_object$ann_clusters_sub <- as.character(Idents(GluR_object))
for(i in unique(GluR_object$ann_clusters_sub)){
  if(min(GluR_markers_avex[GluR_markers_avex$Clusters==i,]$Avex>GluR_thresholds$thresh)==1){
    GluR_object$ann_clusters_sub[GluR_object$ann_clusters_sub==i] <- paste0(i, "_OlfactoryPNs")
  } else{
    GluR_object$ann_clusters_sub[GluR_object$ann_clusters_sub==i] <- paste0(i, "_GluR")
  }
}
Idents(GluR_object) <- GluR_object$ann_clusters_sub

p1 <- FeaturePlot(GluR_object, "VGluRT", order=T) + 
  NoAxes() + 
  NoLegend() + 
  ggtitle("VGluRT") + 
  theme(plot.title=element_text(hjust=.5), panel.border=element_rect(colour="black", size=2))




