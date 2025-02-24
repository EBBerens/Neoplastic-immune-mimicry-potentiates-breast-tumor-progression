#####################################################################################################################
#                       Swarbrucj (Wu et al) Dataset Breast Tumor Analysis Steps                                    #
#####################################################################################################################
# Step 3: Merge Tumor Singlets into Combined Seurat Object                                                          #
# Step 4: Extract Epithelial Cells                                                                                  #
#####################################################################################################################



############################################################
# Step 3: Merge Tumor Singlets into Combined Seurat Object # 
############################################################
# Load Libraries 
library(Seurat)
library(patchwork)
library(dplyr)
library(ggplot2)

#########################
# Load Individual Files #
#########################
# HR-pos Tumors
Swarbrick_CID3941_ER_Total_Singlets <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Total_Singlets/Swarbrick_CID3941_ER_Total_Singlets.rds")
Swarbrick_CID3948_ER_Total_Singlets <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Total_Singlets/Swarbrick_CID3948_ER_Total_Singlets.rds")
Swarbrick_CID4040_ER_Total_Singlets <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Total_Singlets/Swarbrick_CID4040_ER_Total_Singlets.rds")
Swarbrick_CID4067_ER_Total_Singlets <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Total_Singlets/Swarbrick_CID4067_ER_Total_Singlets.rds")
Swarbrick_CID4290A_ER_Total_Singlets <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Total_Singlets/Swarbrick_CID4290A_ER_Total_Singlets.rds")
Swarbrick_CID4398_ER_Total_Singlets <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Total_Singlets/Swarbrick_CID4398_ER_Total_Singlets.rds")
Swarbrick_CID4461_ER_Total_Singlets <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Total_Singlets/Swarbrick_CID4461_ER_Total_Singlets.rds")
Swarbrick_CID4463_ER_Total_Singlets <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Total_Singlets/Swarbrick_CID4463_ER_Total_Singlets.rds")
Swarbrick_CID4471_ER_Total_Singlets <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Total_Singlets/Swarbrick_CID4471_ER_Total_Singlets.rds")
Swarbrick_CID4530N_ER_Total_Singlets <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Total_Singlets/Swarbrick_CID4530N_ER_Total_Singlets.rds")
Swarbrick_CID4535_ER_Total_Singlets <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Total_Singlets/Swarbrick_CID4535_ER_Total_Singlets.rds")

# HER-pos Tumors
Swarbrick_CID3586_HER2_Total_Singlets <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Total_Singlets/Swarbrick_CID3586_HER2_Total_Singlets.rds")
Swarbrick_CID3838_HER2_Total_Singlets <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Total_Singlets/Swarbrick_CID3838_HER2_Total_Singlets.rds")
Swarbrick_CID3921_HER2_Total_Singlets <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Total_Singlets/Swarbrick_CID3921_HER2_Total_Singlets.rds")
Swarbrick_CID4066_HER2_Total_Singlets <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Total_Singlets/Swarbrick_CID4066_HER2_Total_Singlets.rds")
Swarbrick_CID45171_HER2_Total_Singlets <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Total_Singlets/Swarbrick_CID45171_HER2_Total_Singlets.rds")

# TNBC Tumors
Swarbrick_CID3946_TNBC_Total_Singlets <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Total_Singlets/Swarbrick_CID3946_TNBC_Total_Singlets.rds")
Swarbrick_CID3963_TNBC_Total_Singlets <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Total_Singlets/Swarbrick_CID3963_TNBC_Total_Singlets.rds")
Swarbrick_CID4465_TNBC_Total_Singlets <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Total_Singlets/Swarbrick_CID4465_TNBC_Total_Singlets.rds")
Swarbrick_CID4495_TNBC_Total_Singlets <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Total_Singlets/Swarbrick_CID4495_TNBC_Total_Singlets.rds")
Swarbrick_CID4513_TNBC_Total_Singlets <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Total_Singlets/Swarbrick_CID4513_TNBC_Total_Singlets.rds")
Swarbrick_CID4515_TNBC_Total_Singlets <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Total_Singlets/Swarbrick_CID4515_TNBC_Total_Singlets.rds")
Swarbrick_CID4523_TNBC_Total_Singlets <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Total_Singlets/Swarbrick_CID4523_TNBC_Total_Singlets.rds")
Swarbrick_CID44041_TNBC_Total_Singlets <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Total_Singlets/Swarbrick_CID44041_TNBC_Total_Singlets.rds")
Swarbrick_CID44971_TNBC_Total_Singlets <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Total_Singlets/Swarbrick_CID44971_TNBC_Total_Singlets.rds")
Swarbrick_CID44991_TNBC_Total_Singlets <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Total_Singlets/Swarbrick_CID44991_TNBC_Total_Singlets.rds")


##############################################
# Merge Individual Files into Single Dataset #
##############################################
# Merge Datasets
all_Swarbrick_breast_tumors_Total_Singlets <- merge(Swarbrick_CID3941_ER_Total_Singlets, y = c(Swarbrick_CID3948_ER_Total_Singlets, Swarbrick_CID4040_ER_Total_Singlets, Swarbrick_CID4067_ER_Total_Singlets, Swarbrick_CID4290A_ER_Total_Singlets, Swarbrick_CID4398_ER_Total_Singlets, 
                                                                                               Swarbrick_CID4461_ER_Total_Singlets, Swarbrick_CID4463_ER_Total_Singlets, Swarbrick_CID4471_ER_Total_Singlets, Swarbrick_CID4530N_ER_Total_Singlets, Swarbrick_CID4535_ER_Total_Singlets,
                                                                                               Swarbrick_CID3586_HER2_Total_Singlets, Swarbrick_CID3838_HER2_Total_Singlets, Swarbrick_CID3921_HER2_Total_Singlets, Swarbrick_CID4066_HER2_Total_Singlets, Swarbrick_CID45171_HER2_Total_Singlets, 
                                                                                               Swarbrick_CID3946_TNBC_Total_Singlets, Swarbrick_CID3963_TNBC_Total_Singlets, Swarbrick_CID4465_TNBC_Total_Singlets, Swarbrick_CID4495_TNBC_Total_Singlets, Swarbrick_CID4513_TNBC_Total_Singlets,  
                                                                                               Swarbrick_CID4515_TNBC_Total_Singlets, Swarbrick_CID4523_TNBC_Total_Singlets, Swarbrick_CID44041_TNBC_Total_Singlets, Swarbrick_CID44971_TNBC_Total_Singlets, Swarbrick_CID44991_TNBC_Total_Singlets))
# Join Layers
all_Swarbrick_breast_tumors_Total_Singlets <- JoinLayers(all_Swarbrick_breast_tumors_Total_Singlets)

# Remove individual objects
rm(Swarbrick_CID3941_ER_Total_Singlets)
rm(Swarbrick_CID3948_ER_Total_Singlets)
rm(Swarbrick_CID4040_ER_Total_Singlets)
rm(Swarbrick_CID4067_ER_Total_Singlets)
rm(Swarbrick_CID4290A_ER_Total_Singlets)
rm(Swarbrick_CID4398_ER_Total_Singlets)
rm(Swarbrick_CID4461_ER_Total_Singlets)
rm(Swarbrick_CID4463_ER_Total_Singlets)
rm(Swarbrick_CID4471_ER_Total_Singlets)
rm(Swarbrick_CID4530N_ER_Total_Singlets)
rm(Swarbrick_CID4535_ER_Total_Singlets)
rm(Swarbrick_CID3586_HER2_Total_Singlets)
rm(Swarbrick_CID3838_HER2_Total_Singlets)
rm(Swarbrick_CID3921_HER2_Total_Singlets)
rm(Swarbrick_CID4066_HER2_Total_Singlets)
rm(Swarbrick_CID45171_HER2_Total_Singlets)
rm(Swarbrick_CID3946_TNBC_Total_Singlets)
rm(Swarbrick_CID3963_TNBC_Total_Singlets)
rm(Swarbrick_CID4465_TNBC_Total_Singlets)
rm(Swarbrick_CID4495_TNBC_Total_Singlets)
rm(Swarbrick_CID4513_TNBC_Total_Singlets)
rm(Swarbrick_CID4515_TNBC_Total_Singlets)
rm(Swarbrick_CID4523_TNBC_Total_Singlets)
rm(Swarbrick_CID44041_TNBC_Total_Singlets)
rm(Swarbrick_CID44971_TNBC_Total_Singlets)
rm(Swarbrick_CID44991_TNBC_Total_Singlets)
gc()

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples before merge

# Find Variable Features
all_Swarbrick_breast_tumors_Total_Singlets <- FindVariableFeatures(all_Swarbrick_breast_tumors_Total_Singlets, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(all_Swarbrick_breast_tumors_Total_Singlets)
all_Swarbrick_breast_tumors_Total_Singlets <- ScaleData(all_Swarbrick_breast_tumors_Total_Singlets, features = all.genes)


#####################
# REDUCE DIMENSIONS #
#####################
all_Swarbrick_breast_tumors_Total_Singlets <- RunPCA(all_Swarbrick_breast_tumors_Total_Singlets, features = VariableFeatures(object = all_Swarbrick_breast_tumors_Total_Singlets))

# Examine and visualize PCA results a few different ways
print(all_Swarbrick_breast_tumors_Total_Singlets[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(all_Swarbrick_breast_tumors_Total_Singlets, dims = 1:2, reduction = "pca")
DimPlot(all_Swarbrick_breast_tumors_Total_Singlets, reduction = "pca")

# Visualize PCs
ElbowPlot(all_Swarbrick_breast_tumors_Total_Singlets, ndims = 50)

all_Swarbrick_breast_tumors_Total_Singlets <- FindNeighbors(all_Swarbrick_breast_tumors_Total_Singlets, dims = 1:40)
all_Swarbrick_breast_tumors_Total_Singlets <- FindClusters(all_Swarbrick_breast_tumors_Total_Singlets, resolution = 0.09)

# Umap clustering
all_Swarbrick_breast_tumors_Total_Singlets <- RunUMAP(all_Swarbrick_breast_tumors_Total_Singlets, dims = 1:40)
DimPlot(all_Swarbrick_breast_tumors_Total_Singlets, reduction = "umap", raster = FALSE)
DimPlot(all_Swarbrick_breast_tumors_Total_Singlets, reduction = "umap", group.by = "orig.ident", raster = FALSE)

########
# SAVE #
########
saveRDS(all_Swarbrick_breast_tumors_Total_Singlets, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Merged/all_Swarbrick_breast_tumors_Total_Singlets.rds")

####################################
# Find Markers & Annotate Clusters #
####################################

# Perform Marker Analysis
all_Swarbrick_breast_tumors_Total_Singlets.cluster.markers <- FindAllMarkers(all_Swarbrick_breast_tumors_Total_Singlets, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(all_Swarbrick_breast_tumors_Total_Singlets.cluster.markers, file = "/R/R_Swarbrick/Swarbrick_Output/all_Swarbrick_breast_tumors_Total_Singlets.cluster.markers.csv")
rm(all_Swarbrick_breast_tumors_Total_Singlets.cluster.markers)

# Save Example Plots
DimPlot(all_Swarbrick_breast_tumors_Total_Singlets, reduction = "umap", raster = FALSE)
ggsave("/R/R_Swarbrick/Swarbrick_Output/all_Swarbrick_breast_tumors_Total_Singlets_Clusters.tiff", plot = last_plot(), device = "tiff",
       scale = 1, width = 16, height = 10,
       dpi = 200, limitsize = TRUE)

DimPlot(all_Swarbrick_breast_tumors_Total_Singlets, reduction = "umap", group.by = "orig.ident", raster = FALSE)
ggsave("/R/R_Swarbrick/Swarbrick_Output/all_Swarbrick_breast_tumors_Total_Singlets_Clusters_ident.tiff", plot = last_plot(), device = "tiff",
       scale = 1, width = 16, height = 10,
       dpi = 200, limitsize = TRUE)


########################################################################
##### Run SCSA in Python & then annotate clusters using the results ####
########################################################################
##################################################
######## Run in Terminal #########################
##################################################
#### Install SCSA to Home Directory, and create folders for "Swarbrick_Markers" and "Swarbrick_Annotations" within it
#### Place cluster markers file into "Swarbrick_Markers" and change "avg_log2FC" column to "avg_logFC"
#### Run the following within the terminal:
#### Change working directory
# cd ~/SCSA/

#### Run SCSA
# python3 SCSA.py -d whole_v2.db -i ~/SCSA/Swarbrick_Markers/all_Swarbrick_breast_tumors_Total_Singlets.cluster.markers.csv -s seurat -k All -E -g Human -p 0.01 -f 1.5 -o ~/SCSA/Swarbrick_Annotations/all_Swarbrick_breast_tumors_Total_Singlets.cluster.markers.txt -m txt

### Use SCSA results and known marker genes to facilitate annotations

############################################
######## Back to R #########################
############################################
# all_Swarbrick_breast_tumors_Total_Singlets <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Merged/all_Swarbrick_breast_tumors_Total_Singlets.rds")


# Count Total Cells 
Total <- Idents(all_Swarbrick_breast_tumors_Total_Singlets)
# Total = 95,062 total cells from 26 tumors



# Subpopulations will be annotated here next
# Draw Plots
DimPlot(all_Swarbrick_breast_tumors_Total_Singlets, reduction = "umap", raster = FALSE)
FeaturePlot(all_Swarbrick_breast_tumors_Total_Singlets, c("KRT5", "KRT14", "KRT8", "KRT18", "KRT19", 
                                                          "VIM", "PTPRC", "CD3D", "APOE", "MS4A1", "PECAM1", "COL1A1"),  raster = FALSE)
# Annotate Clusters
all_Swarbrick_breast_tumors_Total_Singlets_IDs <- c("T & NK Cells", "Monocytes & Macrophages", "Endothelial Cells", "Fibroblasts", "Pericytes",
                                                    "T & NK Cells", "Epithelial/Neoplastic Cells", "Epithelial/Neoplastic Cells", "Plasma Cells", "B Cells",
                                                    "Epithelial/Neoplastic Cells", "Epithelial/Neoplastic Cells", "Epithelial/Neoplastic Cells", "Epithelial/Neoplastic Cells", "Mesenchymal Cells",
                                                    "Epithelial/Neoplastic Cells", "Epithelial/Neoplastic Cells", "Epithelial/Neoplastic Cells", "Epithelial/Neoplastic Cells", "Epithelial/Neoplastic Cells",
                                                    "Epithelial/Neoplastic Cells", "Mesenchymal Cells", "Epithelial/Neoplastic Cells", "Epithelial/Neoplastic Cells", "Dendritic Cells")

names(all_Swarbrick_breast_tumors_Total_Singlets_IDs) <- levels(all_Swarbrick_breast_tumors_Total_Singlets)
all_Swarbrick_breast_tumors_Total_Singlets <- RenameIdents(all_Swarbrick_breast_tumors_Total_Singlets, all_Swarbrick_breast_tumors_Total_Singlets_IDs)

# Store Ident under Classification
all_Swarbrick_breast_tumors_Total_Singlets[["cell.ident"]] <- Idents(object = all_Swarbrick_breast_tumors_Total_Singlets)

# Draw plot
DimPlot(all_Swarbrick_breast_tumors_Total_Singlets, reduction = "umap", raster = FALSE)


# Draw plot with changed colors
DimPlot(all_Swarbrick_breast_tumors_Total_Singlets, reduction = "umap", raster = FALSE, cols = c('T & NK Cells' = '#0CB702', 'Monocytes & Macrophages' = '#E68613', 'Fibroblasts' = '#8494FF', 'Epithelial/Neoplastic Cells' = '#F8766D',
                                                                                                 'Plasma Cells' = '#00B8E7', 'B Cells' = '#00A9FF', 'Pericytes' = '#ED68ED', 'Endothelial Cells' = '#FF61CC',
                                                                                                 'Dendritic Cells' = '#ABA300',  'Mesenchymal Cells' = '#7CAE00'))
# Save plots
ggsave("/R/R_Swarbrick/Swarbrick_Output/all_Swarbrick_breast_tumors_Total_Singlets_Annotated_UMAP_Custom_Colors.tiff", plot = last_plot(), device = "tiff", 
       scale = 1, width = 16, height = 10,
       dpi = 200, limitsize = TRUE)


FeaturePlot(all_Swarbrick_breast_tumors_Total_Singlets, c("KRT14"),  raster = FALSE, cols = c("grey90", "firebrick"))
ggsave("/R/R_Swarbrick/Swarbrick_Output/all_Swarbrick_breast_tumors_Total_Singlets_KRT14_UMAP_Custom_Colors.tiff", plot = last_plot(), device = "tiff", 
       scale = 1, width = 16, height = 10,
       dpi = 200, limitsize = TRUE)

FeaturePlot(all_Swarbrick_breast_tumors_Total_Singlets, c("KRT18"),  raster = FALSE, cols = c("grey90", "firebrick"))
ggsave("/R/R_Swarbrick/Swarbrick_Output/all_Swarbrick_breast_tumors_Total_Singlets_KRT18_UMAP_Custom_Colors.tiff", plot = last_plot(), device = "tiff", 
       scale = 1, width = 16, height = 10,
       dpi = 200, limitsize = TRUE)

FeaturePlot(all_Swarbrick_breast_tumors_Total_Singlets, c("KRT19"),  raster = FALSE, cols = c("grey90", "firebrick"))
ggsave("/R/R_Swarbrick/Swarbrick_Output/all_Swarbrick_breast_tumors_Total_Singlets_KRT19_UMAP_Custom_Colors.tiff", plot = last_plot(), device = "tiff", 
       scale = 1, width = 16, height = 10,
       dpi = 200, limitsize = TRUE)

# Save RDS
saveRDS(all_Swarbrick_breast_tumors_Total_Singlets, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Annotated/all_Swarbrick_breast_tumors_Total_Singlets_Annotated.rds")



####################################
# Step 4: Extract Epithelial Cells # 
####################################
# all_Swarbrick_breast_tumors_Total_Singlets <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Annotated/all_Swarbrick_breast_tumors_Total_Singlets_Annotated.rds")

# Create ident order for plots
singlet_ident <- c("Epithelial/Neoplastic Cells", "B Cells", "Dendritic Cells", "Endothelial Cells", "Fibroblasts", "Mast Cells", "Mesenchymal Cells", "Monocytes & Macrophages", "Pericytes", "Plasma Cells", "T & NK Cells")
all_Swarbrick_breast_tumors_Total_Singlets@active.ident <- factor(x = all_Swarbrick_breast_tumors_Total_Singlets@active.ident, levels = singlet_ident)
rm(singlet_ident)

# Visualize Cutoffs & Save Example Plots
VlnPlot(all_Swarbrick_breast_tumors_Total_Singlets, c("KRT14"), pt.size = 0.01, raster = FALSE)
ggsave("/R/R_Swarbrick/Swarbrick_Output/all_Swarbrick_breast_tumors_Total_Singlets_VLN_KRT14.tiff", plot = last_plot(), device = "tiff",
       scale = 1, width = 16, height = 10,
       dpi = 200, limitsize = TRUE)

VlnPlot(all_Swarbrick_breast_tumors_Total_Singlets, c("KRT18"), pt.size = 0.01, raster = FALSE)
# Save Example Plots
ggsave("/R/R_Swarbrick/Swarbrick_Output/all_Swarbrick_breast_tumors_Total_Singlets_VLN_KRT18.tiff", plot = last_plot(), device = "tiff",
       scale = 1, width = 16, height = 10,
       dpi = 200, limitsize = TRUE)

VlnPlot(all_Swarbrick_breast_tumors_Total_Singlets, c("KRT19"), pt.size = 0.01, raster = FALSE)
# Save Example Plots
ggsave("/R/R_Swarbrick/Swarbrick_Output/all_Swarbrick_breast_tumors_Total_Singlets_VLN_KRT19.tiff", plot = last_plot(), device = "tiff",
       scale = 1, width = 16, height = 10,
       dpi = 200, limitsize = TRUE)

# Visualize Cutoffs & Save Example Plots - No Points
VlnPlot(all_Swarbrick_breast_tumors_Total_Singlets, c("KRT14"), pt.size = 0, raster = FALSE)
ggsave("/R/R_Swarbrick/Swarbrick_Output/all_Swarbrick_breast_tumors_Total_Singlets_VLN_KRT14_No_Points.tiff", plot = last_plot(), device = "tiff",
       scale = 1, width = 16, height = 10,
       dpi = 200, limitsize = TRUE)

VlnPlot(all_Swarbrick_breast_tumors_Total_Singlets, c("KRT18"), pt.size = 0, raster = FALSE)
# Save Example Plots
ggsave("/R/R_Swarbrick/Swarbrick_Output/all_Swarbrick_breast_tumors_Total_Singlets_VLN_KRT18_No_Points.tiff", plot = last_plot(), device = "tiff",
       scale = 1, width = 16, height = 10,
       dpi = 200, limitsize = TRUE)

VlnPlot(all_Swarbrick_breast_tumors_Total_Singlets, c("KRT19"), pt.size = 0, raster = FALSE)
# Save Example Plots
ggsave("/R/R_Swarbrick/Swarbrick_Output/all_Swarbrick_breast_tumors_Total_Singlets_VLN_KRT19_No_Points.tiff", plot = last_plot(), device = "tiff",
       scale = 1, width = 16, height = 10,
       dpi = 200, limitsize = TRUE)

# KRT14, KRT18, KRT19
VlnPlot(all_Swarbrick_breast_tumors_Total_Singlets, c("KRT14", "KRT18", "KRT19"), ncol = 3, pt.size = 0.1)
ggsave("/R/R_Swarbrick/Swarbrick_Output/all_Swarbrick_breast_tumors_Total_Singlets_VLN_All_KRT_VLN.tiff", plot = last_plot(), device = "tiff",
       scale = 1, width = 20, height = 16,
       dpi = 200, limitsize = TRUE)

# KRT14, KRT18, KRT19 - No Points 
VlnPlot(all_Swarbrick_breast_tumors_Total_Singlets, c("KRT14", "KRT18", "KRT19"), ncol = 3, pt.size = 0)
ggsave("/R/R_Swarbrick/Swarbrick_Output/all_Swarbrick_breast_tumors_Total_Singlets_VLN_All_KRT_VLN_No_Points.tiff", plot = last_plot(), device = "tiff",
       scale = 1, width = 20, height = 16,
       dpi = 200, limitsize = TRUE)


# KRT14, KRT18, KRT19 - Transparent Points 
# Create violin plots with geom_jitter for each feature
# Extra margin added to left-most side
features <- c("KRT14", "KRT18", "KRT19")
plots <- lapply(features, function(feature) {
  VlnPlot(all_Swarbrick_breast_tumors_Total_Singlets, features = feature, pt.size = 0, raster = FALSE) +
    geom_jitter(width = 0.35, alpha = 0.1, size = 0.1) + 
    theme(legend.position = "none", plot.margin = margin(5, 5, 5, 50))
})

# Combine the plots into one
wrap_plots(plots, ncol = 3)
ggsave("/R/R_Swarbrick/Swarbrick_Output/all_Swarbrick_breast_tumors_Total_Singlets_VLN_All_KRT_VLN_Transparent_Points.tiff", plot = last_plot(), device = "tiff",
       scale = 1, width = 20, height = 16,
       dpi = 200, limitsize = TRUE)

rm(features)
rm(plots)

###############################################
########### Subset Epithelial Cells ###########
###############################################

# Select cells greater than cutoffs for at least one of the three main mammary epithelial cell types
all_Swarbrick_breast_tumors_panCK <- subset(x = all_Swarbrick_breast_tumors_Total_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples before merge

# Find Variable Features
all_Swarbrick_breast_tumors_panCK <- FindVariableFeatures(all_Swarbrick_breast_tumors_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(all_Swarbrick_breast_tumors_panCK)
all_Swarbrick_breast_tumors_panCK <- ScaleData(all_Swarbrick_breast_tumors_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
all_Swarbrick_breast_tumors_panCK <- RunPCA(all_Swarbrick_breast_tumors_panCK, features = VariableFeatures(object = all_Swarbrick_breast_tumors_panCK))

# Examine and visualize PCA results a few different ways
print(all_Swarbrick_breast_tumors_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(all_Swarbrick_breast_tumors_panCK, dims = 1:2, reduction = "pca")
DimPlot(all_Swarbrick_breast_tumors_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(all_Swarbrick_breast_tumors_panCK, ndims = 50)

all_Swarbrick_breast_tumors_panCK <- FindNeighbors(all_Swarbrick_breast_tumors_panCK, dims = 1:40)
all_Swarbrick_breast_tumors_panCK <- FindClusters(all_Swarbrick_breast_tumors_panCK, resolution = 0.15)

# Umap clustering
all_Swarbrick_breast_tumors_panCK <- RunUMAP(all_Swarbrick_breast_tumors_panCK, dims = 1:40)
DimPlot(all_Swarbrick_breast_tumors_panCK, reduction = "umap", raster = FALSE)
DimPlot(all_Swarbrick_breast_tumors_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(all_Swarbrick_breast_tumors_panCK, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Merged/all_Swarbrick_breast_tumors_panCK.rds")

########################################################################
# Extract SeuratObject by original identity and save as epithelial RDS #
########################################################################

# HR-pos
Swarbrick_CID3941_ER_panCK <- subset(x = all_Swarbrick_breast_tumors_panCK, subset = orig.ident == "CID3941")
saveRDS(Swarbrick_CID3941_ER_panCK, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_panCK/Swarbrick_CID3941_ER_panCK.rds")

Swarbrick_CID3948_ER_panCK <- subset(x = all_Swarbrick_breast_tumors_panCK, subset = orig.ident == "CID3948")
saveRDS(Swarbrick_CID3948_ER_panCK, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_panCK/Swarbrick_CID3948_ER_panCK.rds")

Swarbrick_CID4040_ER_panCK <- subset(x = all_Swarbrick_breast_tumors_panCK, subset = orig.ident == "CID4040")
saveRDS(Swarbrick_CID4040_ER_panCK, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_panCK/Swarbrick_CID4040_ER_panCK.rds")

Swarbrick_CID4067_ER_panCK <- subset(x = all_Swarbrick_breast_tumors_panCK, subset = orig.ident == "CID4067")
saveRDS(Swarbrick_CID4067_ER_panCK, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_panCK/Swarbrick_CID4067_ER_panCK.rds")

Swarbrick_CID4290A_ER_panCK <- subset(x = all_Swarbrick_breast_tumors_panCK, subset = orig.ident == "CID4290A")
saveRDS(Swarbrick_CID4290A_ER_panCK, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_panCK/Swarbrick_CID4290A_ER_panCK.rds")

Swarbrick_CID4398_ER_panCK <- subset(x = all_Swarbrick_breast_tumors_panCK, subset = orig.ident == "CID4398")
saveRDS(Swarbrick_CID4398_ER_panCK, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_panCK/Swarbrick_CID4398_ER_panCK.rds")

Swarbrick_CID4461_ER_panCK <- subset(x = all_Swarbrick_breast_tumors_panCK, subset = orig.ident == "CID4461")
saveRDS(Swarbrick_CID4461_ER_panCK, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_panCK/Swarbrick_CID4461_ER_panCK.rds")

Swarbrick_CID4463_ER_panCK <- subset(x = all_Swarbrick_breast_tumors_panCK, subset = orig.ident == "CID4463")
saveRDS(Swarbrick_CID4463_ER_panCK, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_panCK/Swarbrick_CID4463_ER_panCK.rds")

Swarbrick_CID4471_ER_panCK <- subset(x = all_Swarbrick_breast_tumors_panCK, subset = orig.ident == "CID4471")
saveRDS(Swarbrick_CID4471_ER_panCK, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_panCK/Swarbrick_CID4471_ER_panCK.rds")

Swarbrick_CID4530N_ER_panCK <- subset(x = all_Swarbrick_breast_tumors_panCK, subset = orig.ident == "CID4530N")
saveRDS(Swarbrick_CID4530N_ER_panCK, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_panCK/Swarbrick_CID4530N_ER_panCK.rds")

Swarbrick_CID4535_ER_panCK <- subset(x = all_Swarbrick_breast_tumors_panCK, subset = orig.ident == "CID4535")
saveRDS(Swarbrick_CID4535_ER_panCK, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_panCK/Swarbrick_CID4535_ER_panCK.rds")

# HER2-pos
Swarbrick_CID3586_HER2_panCK <- subset(x = all_Swarbrick_breast_tumors_panCK, subset = orig.ident == "CID3586")
saveRDS(Swarbrick_CID3586_HER2_panCK, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_panCK/Swarbrick_CID3586_HER2_panCK.rds")

Swarbrick_CID3838_HER2_panCK <- subset(x = all_Swarbrick_breast_tumors_panCK, subset = orig.ident == "CID3838")
saveRDS(Swarbrick_CID3838_HER2_panCK, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_panCK/Swarbrick_CID3838_HER2_panCK.rds")

Swarbrick_CID3921_HER2_panCK <- subset(x = all_Swarbrick_breast_tumors_panCK, subset = orig.ident == "CID3921")
saveRDS(Swarbrick_CID3921_HER2_panCK, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_panCK/Swarbrick_CID3921_HER2_panCK.rds")

Swarbrick_CID4066_HER2_panCK <- subset(x = all_Swarbrick_breast_tumors_panCK, subset = orig.ident == "CID4066")
saveRDS(Swarbrick_CID4066_HER2_panCK, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_panCK/Swarbrick_CID4066_HER2_panCK.rds")

Swarbrick_CID45171_HER2_panCK <- subset(x = all_Swarbrick_breast_tumors_panCK, subset = orig.ident == "CID45171")
saveRDS(Swarbrick_CID45171_HER2_panCK, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_panCK/Swarbrick_CID45171_HER2_panCK.rds")

# TNBC
Swarbrick_CID3946_TNBC_panCK <- subset(x = all_Swarbrick_breast_tumors_panCK, subset = orig.ident == "CID3946")
saveRDS(Swarbrick_CID3946_TNBC_panCK, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_panCK/Swarbrick_CID3946_TNBC_panCK.rds")

Swarbrick_CID3963_TNBC_panCK <- subset(x = all_Swarbrick_breast_tumors_panCK, subset = orig.ident == "CID3963")
saveRDS(Swarbrick_CID3963_TNBC_panCK, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_panCK/Swarbrick_CID3963_TNBC_panCK.rds")

Swarbrick_CID4465_TNBC_panCK <- subset(x = all_Swarbrick_breast_tumors_panCK, subset = orig.ident == "CID4465")
saveRDS(Swarbrick_CID4465_TNBC_panCK, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_panCK/Swarbrick_CID4465_TNBC_panCK.rds")

Swarbrick_CID4495_TNBC_panCK <- subset(x = all_Swarbrick_breast_tumors_panCK, subset = orig.ident == "CID4495")
saveRDS(Swarbrick_CID4495_TNBC_panCK, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_panCK/Swarbrick_CID4495_TNBC_panCK.rds")

Swarbrick_CID4513_TNBC_panCK <- subset(x = all_Swarbrick_breast_tumors_panCK, subset = orig.ident == "CID4513")
saveRDS(Swarbrick_CID4513_TNBC_panCK, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_panCK/Swarbrick_CID4513_TNBC_panCK.rds")

Swarbrick_CID4515_TNBC_panCK <- subset(x = all_Swarbrick_breast_tumors_panCK, subset = orig.ident == "CID4515")
saveRDS(Swarbrick_CID4515_TNBC_panCK, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_panCK/Swarbrick_CID4515_TNBC_panCK.rds")

Swarbrick_CID4523_TNBC_panCK <- subset(x = all_Swarbrick_breast_tumors_panCK, subset = orig.ident == "CID4523")
saveRDS(Swarbrick_CID4523_TNBC_panCK, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_panCK/Swarbrick_CID4523_TNBC_panCK.rds")

Swarbrick_CID44041_TNBC_panCK <- subset(x = all_Swarbrick_breast_tumors_panCK, subset = orig.ident == "CID44041")
saveRDS(Swarbrick_CID44041_TNBC_panCK, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_panCK/Swarbrick_CID44041_TNBC_panCK.rds")

Swarbrick_CID44971_TNBC_panCK <- subset(x = all_Swarbrick_breast_tumors_panCK, subset = orig.ident == "CID44971")
saveRDS(Swarbrick_CID44971_TNBC_panCK, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_panCK/Swarbrick_CID44971_TNBC_panCK.rds")

Swarbrick_CID44991_TNBC_panCK <- subset(x = all_Swarbrick_breast_tumors_panCK, subset = orig.ident == "CID44991")
saveRDS(Swarbrick_CID44991_TNBC_panCK, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_panCK/Swarbrick_CID44991_TNBC_panCK.rds")


# Remove individual objects
rm(all_Swarbrick_breast_tumors_panCK)
rm(all_Swarbrick_breast_tumors_panCK)
rm(Swarbrick_CID3941_ER_panCK)
rm(Swarbrick_CID3948_ER_panCK)
rm(Swarbrick_CID4040_ER_panCK)
rm(Swarbrick_CID4067_ER_panCK) 
rm(Swarbrick_CID4290A_ER_panCK) 
rm(Swarbrick_CID4398_ER_panCK) 
rm(Swarbrick_CID4461_ER_panCK)
rm(Swarbrick_CID4463_ER_panCK)
rm(Swarbrick_CID4471_ER_panCK)
rm(Swarbrick_CID4530N_ER_panCK)
rm(Swarbrick_CID4535_ER_panCK)
rm(Swarbrick_CID3586_HER2_panCK)
rm(Swarbrick_CID3838_HER2_panCK)
rm(Swarbrick_CID3921_HER2_panCK)
rm(Swarbrick_CID4066_HER2_panCK)
rm(Swarbrick_CID45171_HER2_panCK)
rm(Swarbrick_CID3946_TNBC_panCK)
rm(Swarbrick_CID3963_TNBC_panCK)
rm(Swarbrick_CID4465_TNBC_panCK)
rm(Swarbrick_CID4495_TNBC_panCK)
rm(Swarbrick_CID4513_TNBC_panCK)
rm(Swarbrick_CID4515_TNBC_panCK)
rm(Swarbrick_CID4523_TNBC_panCK)
rm(Swarbrick_CID44041_TNBC_panCK)
rm(Swarbrick_CID44971_TNBC_panCK)
rm(Swarbrick_CID44991_TNBC_panCK)
gc()





