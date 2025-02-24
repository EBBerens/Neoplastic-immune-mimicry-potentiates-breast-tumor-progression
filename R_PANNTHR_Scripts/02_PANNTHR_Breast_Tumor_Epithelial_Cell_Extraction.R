#####################################################################################################################
#                              PANNTHR Dataset Breast Tumor Analysis Steps                                          #
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
# PANNTHR Samples
PANNTHR_Pt1_PreTreat_Singlets <- readRDS(file = "/R/R_PANNTHR/PANNTHR_RDS/RDS_Total_Singlets/PANNTHR_Pt1_PreTreat_Singlets.rds")
PANNTHR_Pt2_PreTreat_Singlets <- readRDS(file = "/R/R_PANNTHR/PANNTHR_RDS/RDS_Total_Singlets/PANNTHR_Pt2_PreTreat_Singlets.rds")
PANNTHR_Pt4_PreTreat_Singlets <- readRDS(file = "/R/R_PANNTHR/PANNTHR_RDS/RDS_Total_Singlets/PANNTHR_Pt4_PreTreat_Singlets.rds")
PANNTHR_Pt5_PreTreat_Singlets <- readRDS(file = "/R/R_PANNTHR/PANNTHR_RDS/RDS_Total_Singlets/PANNTHR_Pt5_PreTreat_Singlets.rds")
PANNTHR_Pt2_OnTreat_Singlets <- readRDS(file = "/R/R_PANNTHR/PANNTHR_RDS/RDS_Total_Singlets/PANNTHR_Pt2_OnTreat_Singlets.rds")
PANNTHR_Pt3_OnTreat_Singlets <- readRDS(file = "/R/R_PANNTHR/PANNTHR_RDS/RDS_Total_Singlets/PANNTHR_Pt3_OnTreat_Singlets.rds")
PANNTHR_Pt4_OnTreat_Singlets <- readRDS(file = "/R/R_PANNTHR/PANNTHR_RDS/RDS_Total_Singlets/PANNTHR_Pt4_OnTreat_Singlets.rds")
PANNTHR_Pt6_OnTreat_Singlets <- readRDS(file = "/R/R_PANNTHR/PANNTHR_RDS/RDS_Total_Singlets/PANNTHR_Pt6_OnTreat_Singlets.rds")
PANNTHR_Pt7_OnTreat_Singlets <- readRDS(file = "/R/R_PANNTHR/PANNTHR_RDS/RDS_Total_Singlets/PANNTHR_Pt7_OnTreat_Singlets.rds")
PANNTHR_Pt5_PostTreat_Singlets <- readRDS(file = "/R/R_PANNTHR/PANNTHR_RDS/RDS_Total_Singlets/PANNTHR_Pt5_PostTreat_Singlets.rds")



##############################################
# Merge Individual Files into Single Dataset #
##############################################
# Merge Datasets
all_PANNTHR_breast_tumors_Total_Singlets <- merge(PANNTHR_Pt1_PreTreat_Singlets, y = c(PANNTHR_Pt2_PreTreat_Singlets, PANNTHR_Pt4_PreTreat_Singlets, PANNTHR_Pt5_PreTreat_Singlets, 
                                                                                       PANNTHR_Pt2_OnTreat_Singlets, PANNTHR_Pt3_OnTreat_Singlets, PANNTHR_Pt4_OnTreat_Singlets,
                                                                                       PANNTHR_Pt6_OnTreat_Singlets, PANNTHR_Pt7_OnTreat_Singlets, PANNTHR_Pt5_PostTreat_Singlets))
# Join Layers
all_PANNTHR_breast_tumors_Total_Singlets <- JoinLayers(all_PANNTHR_breast_tumors_Total_Singlets)

# Remove individual objects
rm(PANNTHR_Pt1_PreTreat_Singlets)
rm(PANNTHR_Pt2_PreTreat_Singlets)
rm(PANNTHR_Pt4_PreTreat_Singlets)
rm(PANNTHR_Pt5_PreTreat_Singlets)
rm(PANNTHR_Pt2_OnTreat_Singlets)
rm(PANNTHR_Pt4_OnTreat_Singlets)
rm(PANNTHR_Pt6_OnTreat_Singlets)
rm(PANNTHR_Pt7_OnTreat_Singlets)
rm(PANNTHR_Pt5_PostTreat_Singlets)
gc()



####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
all_PANNTHR_breast_tumors_Total_Singlets <- FindVariableFeatures(all_PANNTHR_breast_tumors_Total_Singlets, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(all_PANNTHR_breast_tumors_Total_Singlets)
all_PANNTHR_breast_tumors_Total_Singlets <- ScaleData(all_PANNTHR_breast_tumors_Total_Singlets, features = all.genes)


#####################
# REDUCE DIMENSIONS #
#####################
all_PANNTHR_breast_tumors_Total_Singlets <- RunPCA(all_PANNTHR_breast_tumors_Total_Singlets, features = VariableFeatures(object = all_PANNTHR_breast_tumors_Total_Singlets))

# Examine and visualize PCA results a few different ways
print(all_PANNTHR_breast_tumors_Total_Singlets[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(all_PANNTHR_breast_tumors_Total_Singlets, dims = 1:2, reduction = "pca")
DimPlot(all_PANNTHR_breast_tumors_Total_Singlets, reduction = "pca")

# Visualize PCs
ElbowPlot(all_PANNTHR_breast_tumors_Total_Singlets, ndims = 50)

all_PANNTHR_breast_tumors_Total_Singlets <- FindNeighbors(all_PANNTHR_breast_tumors_Total_Singlets, dims = 1:50)
all_PANNTHR_breast_tumors_Total_Singlets <- FindClusters(all_PANNTHR_breast_tumors_Total_Singlets, resolution = 0.15)

# Umap clustering
all_PANNTHR_breast_tumors_Total_Singlets <- RunUMAP(all_PANNTHR_breast_tumors_Total_Singlets, dims = 1:50)
DimPlot(all_PANNTHR_breast_tumors_Total_Singlets, reduction = "umap", raster = FALSE, label = TRUE)
DimPlot(all_PANNTHR_breast_tumors_Total_Singlets, reduction = "umap", group.by = "orig.ident", raster = FALSE)


########
# SAVE #
########
saveRDS(all_PANNTHR_breast_tumors_Total_Singlets, file = "/R/R_PANNTHR/PANNTHR_RDS/RDS_Merged/all_PANNTHR_breast_tumors_Total_Singlets.rds")

####################################
# Find Markers & Annotate Clusters #
####################################

# Perform Marker Analysis
all_PANNTHR_breast_tumors_Total_Singlets.cluster.markers <- FindAllMarkers(all_PANNTHR_breast_tumors_Total_Singlets, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(all_PANNTHR_breast_tumors_Total_Singlets.cluster.markers, file = "/R/R_PANNTHR/PANNTHR_Output/all_PANNTHR_breast_tumors_Total_Singlets.cluster.markers.csv")
rm(all_PANNTHR_breast_tumors_Total_Singlets.cluster.markers)

# Save Example Plots
DimPlot(all_PANNTHR_breast_tumors_Total_Singlets, reduction = "umap", raster = FALSE)
ggsave("/R/R_PANNTHR/PANNTHR_Output/all_PANNTHR_breast_tumors_Total_Singlets_Clusters.tiff", plot = last_plot(), device = "tiff",
       scale = 1, width = 16, height = 10,
       dpi = 200, limitsize = TRUE)

DimPlot(all_PANNTHR_breast_tumors_Total_Singlets, reduction = "umap", raster = FALSE, label = TRUE)
ggsave("/R/R_PANNTHR/PANNTHR_Output/all_PANNTHR_breast_tumors_Total_Singlets_Clusters_Labels.tiff", plot = last_plot(), device = "tiff",
       scale = 1, width = 16, height = 10,
       dpi = 200, limitsize = TRUE)

DimPlot(all_PANNTHR_breast_tumors_Total_Singlets, reduction = "umap", group.by = "orig.ident", raster = FALSE)
ggsave("/R/R_PANNTHR/PANNTHR_Output/all_PANNTHR_breast_tumors_Total_Singlets_Clusters_ident.tiff", plot = last_plot(), device = "tiff",
       scale = 1, width = 16, height = 10,
       dpi = 200, limitsize = TRUE)


########################################################################
##### Run SCSA in Python & then annotate clusters using the results ####
########################################################################
##################################################
######## Run in Terminal #########################
##################################################
#### Install SCSA to Home Directory, and create folders for "PANNTHR_Markers" and "PANNTHR_Annotations" within it
#### Place cluster markers file into "PANNTHR_Markers" and change "avg_log2FC" column to "avg_logFC"
#### Run the following within the terminal:
#### Change working directory
# cd ~/SCSA/

#### Run SCSA
# python3 SCSA.py -d whole_v2.db -i ~/SCSA/PANNTHR_Markers/all_PANNTHR_breast_tumors_Total_Singlets.cluster.markers.csv -s seurat -k All -E -g Human -p 0.01 -f 1.5 -o ~/SCSA/PANNTHR_Annotations/all_PANNTHR_breast_tumors_Total_Singlets.cluster.markers.txt -m txt

### Use SCSA results and known marker genes to facilitate annotations

############################################
######## Back to R #########################
############################################
# all_PANNTHR_breast_tumors_Total_Singlets <- readRDS(file = "/R/R_PANNTHR/PANNTHR_RDS/RDS_Merged/all_PANNTHR_breast_tumors_Total_Singlets.rds")


# Count Total Cells 
Total <- Idents(all_PANNTHR_breast_tumors_Total_Singlets)
# Total = 112,883 total cells from 10 tumor samples


# Subpopulations will be annotated here next
# Draw Plots
DimPlot(all_PANNTHR_breast_tumors_Total_Singlets, reduction = "umap", raster = FALSE)
FeaturePlot(all_PANNTHR_breast_tumors_Total_Singlets, c("KRT5", "KRT14", "KRT8", "KRT18", "KRT19", 
                                                          "VIM", "PTPRC", "CD3D", "APOE", "MS4A1", "PECAM1", "COL1A1"),  raster = FALSE)
# Annotate Clusters
all_PANNTHR_breast_tumors_Total_Singlets_IDs <- c("T & NK Cells", "B Cells", "Epithelial/Neoplastic Cells", "T & NK Cells", "Plasma Cells",
                                                  "T & NK Cells", "Monocytes & Macrophages", "T & NK Cells", "Epithelial/Neoplastic Cells", "Fibroblasts",
                                                  "Monocytes & Macrophages", "Epithelial/Neoplastic Cells", "Monocytes & Macrophages", "Endothelial Cells", "B Cells",
                                                  "Mast Cells", "Epithelial/Neoplastic Cells", "Mast Cells", "Dendritic Cells", "Dendritic Cells")

names(all_PANNTHR_breast_tumors_Total_Singlets_IDs) <- levels(all_PANNTHR_breast_tumors_Total_Singlets)
all_PANNTHR_breast_tumors_Total_Singlets <- RenameIdents(all_PANNTHR_breast_tumors_Total_Singlets, all_PANNTHR_breast_tumors_Total_Singlets_IDs)

# Store Ident under Classification
all_PANNTHR_breast_tumors_Total_Singlets[["cell.ident"]] <- Idents(object = all_PANNTHR_breast_tumors_Total_Singlets)

# Draw plot
DimPlot(all_PANNTHR_breast_tumors_Total_Singlets, reduction = "umap", raster = FALSE)


# Draw plot with changed colors
DimPlot(all_PANNTHR_breast_tumors_Total_Singlets, reduction = "umap", raster = FALSE, cols = c('T & NK Cells' = '#0CB702', 'Monocytes & Macrophages' = '#E68613', 'Fibroblasts' = '#8494FF', 'Epithelial/Neoplastic Cells' = '#F8766D',
                                                                                                 'Plasma Cells' = '#00B8E7', 'B Cells' = '#00A9FF', 'Endothelial Cells' = '#FF61CC',
                                                                                                 'Dendritic Cells' = '#ABA300',  'Mast Cells' = '#7CAE00'))
# Save plots
ggsave("/R/R_PANNTHR/PANNTHR_Output/all_PANNTHR_breast_tumors_Total_Singlets_Annotated_UMAP_Custom_Colors.tiff", plot = last_plot(), device = "tiff", 
       scale = 1, width = 16, height = 10,
       dpi = 200, limitsize = TRUE)


FeaturePlot(all_PANNTHR_breast_tumors_Total_Singlets, c("KRT14"),  raster = FALSE, cols = c("grey90", "firebrick"))
ggsave("/R/R_PANNTHR/PANNTHR_Output/all_PANNTHR_breast_tumors_Total_Singlets_KRT14_UMAP_Custom_Colors.tiff", plot = last_plot(), device = "tiff", 
       scale = 1, width = 16, height = 10,
       dpi = 200, limitsize = TRUE)

FeaturePlot(all_PANNTHR_breast_tumors_Total_Singlets, c("KRT18"),  raster = FALSE, cols = c("grey90", "firebrick"))
ggsave("/R/R_PANNTHR/PANNTHR_Output/all_PANNTHR_breast_tumors_Total_Singlets_KRT18_UMAP_Custom_Colors.tiff", plot = last_plot(), device = "tiff", 
       scale = 1, width = 16, height = 10,
       dpi = 200, limitsize = TRUE)

FeaturePlot(all_PANNTHR_breast_tumors_Total_Singlets, c("KRT19"),  raster = FALSE, cols = c("grey90", "firebrick"))
ggsave("/R/R_PANNTHR/PANNTHR_Output/all_PANNTHR_breast_tumors_Total_Singlets_KRT19_UMAP_Custom_Colors.tiff", plot = last_plot(), device = "tiff", 
       scale = 1, width = 16, height = 10,
       dpi = 200, limitsize = TRUE)

# Save RDS
saveRDS(all_PANNTHR_breast_tumors_Total_Singlets, file = "/R/R_PANNTHR/PANNTHR_RDS/RDS_Annotated/all_PANNTHR_breast_tumors_Total_Singlets_Annotated.rds")



####################################
# Step 4: Extract Epithelial Cells # 
####################################
# all_PANNTHR_breast_tumors_Total_Singlets <- readRDS(file = "/R/R_PANNTHR/PANNTHR_RDS/RDS_Annotated/all_PANNTHR_breast_tumors_Total_Singlets_Annotated.rds")

# Create ident order for plots
singlet_ident <- c("Epithelial/Neoplastic Cells", "B Cells", "Dendritic Cells", "Endothelial Cells", "Fibroblasts", "Mast Cells", "Monocytes & Macrophages", "Plasma Cells", "T & NK Cells")
all_PANNTHR_breast_tumors_Total_Singlets@active.ident <- factor(x = all_PANNTHR_breast_tumors_Total_Singlets@active.ident, levels = singlet_ident)
rm(singlet_ident)

# Visualize Cutoffs & Save Example Plots
VlnPlot(all_PANNTHR_breast_tumors_Total_Singlets, c("KRT14"), pt.size = 0.01, raster = FALSE)
ggsave("/R/R_PANNTHR/PANNTHR_Output/all_PANNTHR_breast_tumors_Total_Singlets_VLN_KRT14.tiff", plot = last_plot(), device = "tiff",
       scale = 1, width = 16, height = 10,
       dpi = 200, limitsize = TRUE)

VlnPlot(all_PANNTHR_breast_tumors_Total_Singlets, c("KRT18"), pt.size = 0.01, raster = FALSE)
# Save Example Plots
ggsave("/R/R_PANNTHR/PANNTHR_Output/all_PANNTHR_breast_tumors_Total_Singlets_VLN_KRT18.tiff", plot = last_plot(), device = "tiff",
       scale = 1, width = 16, height = 10,
       dpi = 200, limitsize = TRUE)

VlnPlot(all_PANNTHR_breast_tumors_Total_Singlets, c("KRT19"), pt.size = 0.01, raster = FALSE)
# Save Example Plots
ggsave("/R/R_PANNTHR/PANNTHR_Output/all_PANNTHR_breast_tumors_Total_Singlets_VLN_KRT19.tiff", plot = last_plot(), device = "tiff",
       scale = 1, width = 16, height = 10,
       dpi = 200, limitsize = TRUE)

# Visualize Cutoffs & Save Example Plots - No Points
VlnPlot(all_PANNTHR_breast_tumors_Total_Singlets, c("KRT14"), pt.size = 0, raster = FALSE)
ggsave("/R/R_PANNTHR/PANNTHR_Output/all_PANNTHR_breast_tumors_Total_Singlets_VLN_KRT14_No_Points.tiff", plot = last_plot(), device = "tiff",
       scale = 1, width = 16, height = 10,
       dpi = 200, limitsize = TRUE)

VlnPlot(all_PANNTHR_breast_tumors_Total_Singlets, c("KRT18"), pt.size = 0, raster = FALSE)
# Save Example Plots
ggsave("/R/R_PANNTHR/PANNTHR_Output/all_PANNTHR_breast_tumors_Total_Singlets_VLN_KRT18_No_Points.tiff", plot = last_plot(), device = "tiff",
       scale = 1, width = 16, height = 10,
       dpi = 200, limitsize = TRUE)

VlnPlot(all_PANNTHR_breast_tumors_Total_Singlets, c("KRT19"), pt.size = 0, raster = FALSE)
# Save Example Plots
ggsave("/R/R_PANNTHR/PANNTHR_Output/all_PANNTHR_breast_tumors_Total_Singlets_VLN_KRT19_No_Points.tiff", plot = last_plot(), device = "tiff",
       scale = 1, width = 16, height = 10,
       dpi = 200, limitsize = TRUE)

# KRT14, KRT18, KRT19
VlnPlot(all_PANNTHR_breast_tumors_Total_Singlets, c("KRT14", "KRT18", "KRT19"), ncol = 3, pt.size = 0.1)
ggsave("/R/R_PANNTHR/PANNTHR_Output/all_PANNTHR_breast_tumors_Total_Singlets_VLN_All_KRT_VLN.tiff", plot = last_plot(), device = "tiff",
       scale = 1, width = 20, height = 16,
       dpi = 200, limitsize = TRUE)

# KRT14, KRT18, KRT19 - No Points 
VlnPlot(all_PANNTHR_breast_tumors_Total_Singlets, c("KRT14", "KRT18", "KRT19"), ncol = 3, pt.size = 0)
ggsave("/R/R_PANNTHR/PANNTHR_Output/all_PANNTHR_breast_tumors_Total_Singlets_VLN_All_KRT_VLN_No_Points.tiff", plot = last_plot(), device = "tiff",
       scale = 1, width = 20, height = 16,
       dpi = 200, limitsize = TRUE)


# KRT14, KRT18, KRT19 - Transparent Points 
# Create violin plots with geom_jitter for each feature
# Extra margin added to left-most side
features <- c("KRT14", "KRT18", "KRT19")
plots <- lapply(features, function(feature) {
  VlnPlot(all_PANNTHR_breast_tumors_Total_Singlets, features = feature, pt.size = 0, raster = FALSE) +
    geom_jitter(width = 0.35, alpha = 0.1, size = 0.1) + 
    theme(legend.position = "none", plot.margin = margin(5, 5, 5, 50))
})

# Combine the plots into one
wrap_plots(plots, ncol = 3)
ggsave("/R/R_PANNTHR/PANNTHR_Output/all_PANNTHR_breast_tumors_Total_Singlets_VLN_All_KRT_VLN_Transparent_Points.tiff", plot = last_plot(), device = "tiff",
       scale = 1, width = 20, height = 16,
       dpi = 200, limitsize = TRUE)

rm(features)
rm(plots)

###############################################
########### Subset Epithelial Cells ###########
###############################################

# Select cells greater than cutoffs for at least one of the three main mammary epithelial cell types
all_PANNTHR_breast_tumors_panCK <- subset(x = all_PANNTHR_breast_tumors_Total_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
all_PANNTHR_breast_tumors_panCK <- FindVariableFeatures(all_PANNTHR_breast_tumors_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(all_PANNTHR_breast_tumors_panCK)
all_PANNTHR_breast_tumors_panCK <- ScaleData(all_PANNTHR_breast_tumors_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
all_PANNTHR_breast_tumors_panCK <- RunPCA(all_PANNTHR_breast_tumors_panCK, features = VariableFeatures(object = all_PANNTHR_breast_tumors_panCK))

# Examine and visualize PCA results a few different ways
print(all_PANNTHR_breast_tumors_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(all_PANNTHR_breast_tumors_panCK, dims = 1:2, reduction = "pca")
DimPlot(all_PANNTHR_breast_tumors_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(all_PANNTHR_breast_tumors_panCK, ndims = 50)

all_PANNTHR_breast_tumors_panCK <- FindNeighbors(all_PANNTHR_breast_tumors_panCK, dims = 1:40)
all_PANNTHR_breast_tumors_panCK <- FindClusters(all_PANNTHR_breast_tumors_panCK, resolution = 0.15)

# Umap clustering
all_PANNTHR_breast_tumors_panCK <- RunUMAP(all_PANNTHR_breast_tumors_panCK, dims = 1:40)
DimPlot(all_PANNTHR_breast_tumors_panCK, reduction = "umap", raster = FALSE)
DimPlot(all_PANNTHR_breast_tumors_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(all_PANNTHR_breast_tumors_panCK, file = "/R/R_PANNTHR/PANNTHR_RDS/RDS_Merged/all_PANNTHR_breast_tumors_panCK.rds")




########################################################################
# Extract SeuratObject by original identity and save as epithelial RDS #
########################################################################

PANNTHR_Pt1_PreTreat_panCK <- subset(x = all_PANNTHR_breast_tumors_panCK, subset = orig.ident == "PANNTHR_Pt1_PreTreat")
saveRDS(PANNTHR_Pt1_PreTreat_panCK, file = "/R/R_PANNTHR/Visvader_RDS/RDS_panCK/PANNTHR_Pt1_PreTreat_panCK.rds")

PANNTHR_Pt2_PreTreat_panCK <- subset(x = all_PANNTHR_breast_tumors_panCK, subset = orig.ident == "PANNTHR_Pt2_PreTreat")
saveRDS(PANNTHR_Pt2_PreTreat_panCK, file = "/R/R_PANNTHR/Visvader_RDS/RDS_panCK/PANNTHR_Pt2_PreTreat_panCK.rds")

PANNTHR_Pt4_PreTreat_panCK <- subset(x = all_PANNTHR_breast_tumors_panCK, subset = orig.ident == "PANNTHR_Pt4_PreTreat")
saveRDS(PANNTHR_Pt4_PreTreat_panCK, file = "/R/R_PANNTHR/Visvader_RDS/RDS_panCK/PANNTHR_Pt4_PreTreat_panCK.rds")


PANNTHR_Pt5_PreTreat_panCK <- subset(x = all_PANNTHR_breast_tumors_panCK, subset = orig.ident == "PANNTHR_Pt5_PreTreat")
saveRDS(PANNTHR_Pt5_PreTreat_panCK, file = "/R/R_PANNTHR/Visvader_RDS/RDS_panCK/PANNTHR_Pt5_PreTreat_panCK.rds")

PANNTHR_Pt2_OnTreat_panCK <- subset(x = all_PANNTHR_breast_tumors_panCK, subset = orig.ident == "PANNTHR_Pt2_OnTreat")
saveRDS(PANNTHR_Pt2_OnTreat_panCK, file = "/R/R_PANNTHR/Visvader_RDS/RDS_panCK/PANNTHR_Pt2_OnTreat_panCK.rds")

PANNTHR_Pt3_OnTreat_panCK <- subset(x = all_PANNTHR_breast_tumors_panCK, subset = orig.ident == "PANNTHR_Pt3_OnTreat")
saveRDS(PANNTHR_Pt3_OnTreat_panCK, file = "/R/R_PANNTHR/Visvader_RDS/RDS_panCK/PANNTHR_Pt3_OnTreat_panCK.rds")

PANNTHR_Pt4_OnTreat_panCK <- subset(x = all_PANNTHR_breast_tumors_panCK, subset = orig.ident == "PANNTHR_Pt4_OnTreat")
saveRDS(PANNTHR_Pt4_OnTreat_panCK, file = "/R/R_PANNTHR/Visvader_RDS/RDS_panCK/PANNTHR_Pt4_OnTreat_panCK.rds")

PANNTHR_Pt6_OnTreat_panCK <- subset(x = all_PANNTHR_breast_tumors_panCK, subset = orig.ident == "PANNTHR_Pt6_OnTreat")
saveRDS(PANNTHR_Pt6_OnTreat_panCK, file = "/R/R_PANNTHR/Visvader_RDS/RDS_panCK/PANNTHR_Pt6_OnTreat_panCK.rds")

PANNTHR_Pt7_OnTreat_panCK <- subset(x = all_PANNTHR_breast_tumors_panCK, subset = orig.ident == "PANNTHR_Pt7_OnTreat")
saveRDS(PANNTHR_Pt7_OnTreat_panCK, file = "/R/R_PANNTHR/Visvader_RDS/RDS_panCK/PANNTHR_Pt7_OnTreat_panCK.rds")

PANNTHR_Pt5_PostTreat_panCK <- subset(x = all_PANNTHR_breast_tumors_panCK, subset = orig.ident == "PANNTHR_Pt5_PostTreat")
saveRDS(PANNTHR_Pt5_PostTreat_panCK, file = "/R/R_PANNTHR/Visvader_RDS/RDS_panCK/PANNTHR_Pt5_PostTreat_panCK.rds")

# Remove individual objects
rm(all_PANNTHR_breast_tumors_Total_Singlets)
rm(all_PANNTHR_breast_tumors_panCK)
rm(PANNTHR_Pt1_PreTreat_panCK)
rm(PANNTHR_Pt2_PreTreat_panCK)
rm(PANNTHR_Pt4_PreTreat_panCK)
rm(PANNTHR_Pt5_PreTreat_panCK)
rm(PANNTHR_Pt2_OnTreat_panCK)
rm(PANNTHR_Pt3_OnTreat_panCK)
rm(PANNTHR_Pt4_OnTreat_panCK)
rm(PANNTHR_Pt5_PostTreat_panCK)
rm(PANNTHR_Pt6_OnTreat_panCK)
rm(PANNTHR_Pt7_OnTreat_panCK)
gc()





