#####################################################################################################################
#                       Visvader (Pal et al) Dataset Breast Tumor Analysis Steps                                    #
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
Visvader_0001_ER_Total_Singlets <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Total_Singlets/Visvader_0001_ER_Total_Singlets.rds")
Visvader_0025_ER_Total_Singlets <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Total_Singlets/Visvader_0025_ER_Total_Singlets.rds")
Visvader_0029_7C_ER_Total_Singlets <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Total_Singlets/Visvader_0029_7C_ER_Total_Singlets.rds")
Visvader_0029_9C_ER_Total_Singlets <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Total_Singlets/Visvader_0029_9C_ER_Total_Singlets.rds")
Visvader_0032_ER_Total_Singlets <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Total_Singlets/Visvader_0032_ER_Total_Singlets.rds")
Visvader_0040_ER_Total_Singlets <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Total_Singlets/Visvader_0040_ER_Total_Singlets.rds")
Visvader_0042_ER_Total_Singlets <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Total_Singlets/Visvader_0042_ER_Total_Singlets.rds")
Visvader_0043_ER_Total_Singlets <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Total_Singlets/Visvader_0043_ER_Total_Singlets.rds")
Visvader_0056_ER_Total_Singlets <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Total_Singlets/Visvader_0056_ER_Total_Singlets.rds")
Visvader_0064_ER_Total_Singlets <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Total_Singlets/Visvader_0064_ER_Total_Singlets.rds")
Visvader_0068_ER_Total_Singlets <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Total_Singlets/Visvader_0068_ER_Total_Singlets.rds")
Visvader_0114_ER_Total_Singlets <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Total_Singlets/Visvader_0114_ER_Total_Singlets.rds")
Visvader_0125_ER_Total_Singlets <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Total_Singlets/Visvader_0125_ER_Total_Singlets.rds")
Visvader_0151_ER_Total_Singlets <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Total_Singlets/Visvader_0151_ER_Total_Singlets.rds")
Visvader_0163_ER_Total_Singlets <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Total_Singlets/Visvader_0163_ER_Total_Singlets.rds")
Visvader_0167_ER_Total_Singlets <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Total_Singlets/Visvader_0167_ER_Total_Singlets.rds")
Visvader_0173_ER_Total_Singlets <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Total_Singlets/Visvader_0173_ER_Total_Singlets.rds")
Visvader_0178_ER_Total_Singlets <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Total_Singlets/Visvader_0178_ER_Total_Singlets.rds")
Visvader_0360_ER_Total_Singlets <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Total_Singlets/Visvader_0360_ER_Total_Singlets.rds")
Visvader_0319_PR_Total_Singlets <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Total_Singlets/Visvader_0319_PR_Total_Singlets.rds")

# HER-pos Tumors
Visvader_0031_HER2_Total_Singlets <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Total_Singlets/Visvader_0031_HER2_Total_Singlets.rds")
Visvader_0069_HER2_Total_Singlets <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Total_Singlets/Visvader_0069_HER2_Total_Singlets.rds")
Visvader_0161_HER2_Total_Singlets <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Total_Singlets/Visvader_0161_HER2_Total_Singlets.rds")
Visvader_0176_HER2_Total_Singlets <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Total_Singlets/Visvader_0176_HER2_Total_Singlets.rds")
Visvader_0308_HER2_Total_Singlets <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Total_Singlets/Visvader_0308_HER2_Total_Singlets.rds")
Visvader_0337_HER2_Total_Singlets <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Total_Singlets/Visvader_0337_HER2_Total_Singlets.rds")

# TNBC Tumors
Visvader_106_TNBC_Total_Singlets <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Total_Singlets/Visvader_106_TNBC_Total_Singlets.rds")
Visvader_114_TNBC_Total_Singlets <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Total_Singlets/Visvader_114_TNBC_Total_Singlets.rds")
Visvader_126_TNBC_Total_Singlets <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Total_Singlets/Visvader_126_TNBC_Total_Singlets.rds")
Visvader_0131_TNBC_Total_Singlets <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Total_Singlets/Visvader_0131_TNBC_Total_Singlets.rds")
Visvader_135_TNBC_Total_Singlets <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Total_Singlets/Visvader_135_TNBC_Total_Singlets.rds")
Visvader_0177_TNBC_Total_Singlets <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Total_Singlets/Visvader_0177_TNBC_Total_Singlets.rds")
Visvader_0554_TNBC_Total_Singlets <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Total_Singlets/Visvader_0554_TNBC_Total_Singlets.rds")
Visvader_4031_TNBC_Total_Singlets <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Total_Singlets/Visvader_4031_TNBC_Total_Singlets.rds")


##############################################
# Merge Individual Files into Single Dataset #
##############################################
# Merge Datasets
all_Visvader_breast_tumors_Total_Singlets <- merge(Visvader_0001_ER_Total_Singlets, y = c(Visvader_0025_ER_Total_Singlets, Visvader_0029_7C_ER_Total_Singlets, Visvader_0029_9C_ER_Total_Singlets, Visvader_0032_ER_Total_Singlets,
                                                                                          Visvader_0040_ER_Total_Singlets, Visvader_0042_ER_Total_Singlets, Visvader_0043_ER_Total_Singlets, Visvader_0056_ER_Total_Singlets, Visvader_0064_ER_Total_Singlets,
                                                                                          Visvader_0068_ER_Total_Singlets, Visvader_0114_ER_Total_Singlets, Visvader_0125_ER_Total_Singlets, Visvader_0151_ER_Total_Singlets, Visvader_0163_ER_Total_Singlets,
                                                                                          Visvader_0167_ER_Total_Singlets, Visvader_0173_ER_Total_Singlets, Visvader_0178_ER_Total_Singlets, Visvader_0360_ER_Total_Singlets, Visvader_0319_PR_Total_Singlets,
                                                                                          Visvader_0031_HER2_Total_Singlets, Visvader_0069_HER2_Total_Singlets, Visvader_0161_HER2_Total_Singlets, Visvader_0176_HER2_Total_Singlets, Visvader_0308_HER2_Total_Singlets, 
                                                                                          Visvader_0337_HER2_Total_Singlets, Visvader_106_TNBC_Total_Singlets, Visvader_114_TNBC_Total_Singlets, Visvader_126_TNBC_Total_Singlets, Visvader_0131_TNBC_Total_Singlets, 
                                                                                          Visvader_135_TNBC_Total_Singlets, Visvader_0177_TNBC_Total_Singlets, Visvader_0554_TNBC_Total_Singlets, Visvader_4031_TNBC_Total_Singlets))
# Join Layers
all_Visvader_breast_tumors_Total_Singlets <- JoinLayers(all_Visvader_breast_tumors_Total_Singlets)

# Remove individual objects
rm(Visvader_0001_ER_Total_Singlets)
rm(Visvader_0025_ER_Total_Singlets)
rm(Visvader_0029_7C_ER_Total_Singlets)
rm(Visvader_0029_9C_ER_Total_Singlets)
rm(Visvader_0032_ER_Total_Singlets)
rm(Visvader_0040_ER_Total_Singlets)
rm(Visvader_0042_ER_Total_Singlets)
rm(Visvader_0043_ER_Total_Singlets)
rm(Visvader_0056_ER_Total_Singlets)
rm(Visvader_0064_ER_Total_Singlets)
rm(Visvader_0068_ER_Total_Singlets)
rm(Visvader_0114_ER_Total_Singlets)
rm(Visvader_0125_ER_Total_Singlets)
rm(Visvader_0151_ER_Total_Singlets)
rm(Visvader_0163_ER_Total_Singlets)
rm(Visvader_0167_ER_Total_Singlets)
rm(Visvader_0173_ER_Total_Singlets)
rm(Visvader_0178_ER_Total_Singlets)
rm(Visvader_0360_ER_Total_Singlets)
rm(Visvader_0319_PR_Total_Singlets)
rm(Visvader_0031_HER2_Total_Singlets)
rm(Visvader_0069_HER2_Total_Singlets)
rm(Visvader_0161_HER2_Total_Singlets)
rm(Visvader_0176_HER2_Total_Singlets)
rm(Visvader_0308_HER2_Total_Singlets)
rm(Visvader_0337_HER2_Total_Singlets)
rm(Visvader_106_TNBC_Total_Singlets)
rm(Visvader_114_TNBC_Total_Singlets)
rm(Visvader_126_TNBC_Total_Singlets)
rm(Visvader_0131_TNBC_Total_Singlets)
rm(Visvader_135_TNBC_Total_Singlets)
rm(Visvader_0177_TNBC_Total_Singlets)
rm(Visvader_0554_TNBC_Total_Singlets)
rm(Visvader_4031_TNBC_Total_Singlets)
gc()

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
all_Visvader_breast_tumors_Total_Singlets <- FindVariableFeatures(all_Visvader_breast_tumors_Total_Singlets, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(all_Visvader_breast_tumors_Total_Singlets)
all_Visvader_breast_tumors_Total_Singlets <- ScaleData(all_Visvader_breast_tumors_Total_Singlets, features = all.genes)


#####################
# REDUCE DIMENSIONS #
#####################
all_Visvader_breast_tumors_Total_Singlets <- RunPCA(all_Visvader_breast_tumors_Total_Singlets, features = VariableFeatures(object = all_Visvader_breast_tumors_Total_Singlets))

# Examine and visualize PCA results a few different ways
print(all_Visvader_breast_tumors_Total_Singlets[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(all_Visvader_breast_tumors_Total_Singlets, dims = 1:2, reduction = "pca")
DimPlot(all_Visvader_breast_tumors_Total_Singlets, reduction = "pca")

# Visualize PCs
ElbowPlot(all_Visvader_breast_tumors_Total_Singlets, ndims = 50)

all_Visvader_breast_tumors_Total_Singlets <- FindNeighbors(all_Visvader_breast_tumors_Total_Singlets, dims = 1:50)
all_Visvader_breast_tumors_Total_Singlets <- FindClusters(all_Visvader_breast_tumors_Total_Singlets, resolution = 0.12)

# Umap clustering
all_Visvader_breast_tumors_Total_Singlets <- RunUMAP(all_Visvader_breast_tumors_Total_Singlets, dims = 1:50)
DimPlot(all_Visvader_breast_tumors_Total_Singlets, reduction = "umap", raster = FALSE)
DimPlot(all_Visvader_breast_tumors_Total_Singlets, reduction = "umap", group.by = "orig.ident", raster = FALSE)

########
# SAVE #
########
saveRDS(all_Visvader_breast_tumors_Total_Singlets, file = "/R/R_Visvader/Visvader_RDS/RDS_Merged/all_Visvader_breast_tumors_Total_Singlets.rds")

#all_Visvader_breast_tumors_Total_Singlets <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Merged/all_Visvader_breast_tumors_Total_Singlets.rds")


# Count Total Cells 
Total <- Idents(all_Visvader_breast_tumors_Total_Singlets)
# Total = 156,613 total cells from 34 tumors

####################################
# Find Markers & Annotate Clusters #
####################################

# Perform Marker Analysis
all_Visvader_breast_tumors_Total_Singlets.cluster.markers <- FindAllMarkers(all_Visvader_breast_tumors_Total_Singlets, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(all_Visvader_breast_tumors_Total_Singlets.cluster.markers, file = "/R/R_Visvader/Visvader_Output/all_Visvader_breast_tumors_Total_Singlets.cluster.markers.csv")
rm(all_Visvader_breast_tumors_Total_Singlets.cluster.markers)

# Save Example Plots
DimPlot(all_Visvader_breast_tumors_Total_Singlets, reduction = "umap", raster = FALSE)
ggsave("/R/R_Visvader/Visvader_Output/all_Visvader_breast_tumors_Total_Singlets_Clusters.tiff", plot = last_plot(), device = "tiff",
       scale = 1, width = 16, height = 10,
       dpi = 200, limitsize = TRUE)

DimPlot(all_Visvader_breast_tumors_Total_Singlets, reduction = "umap", group.by = "orig.ident", raster = FALSE)
ggsave("/R/R_Visvader/Visvader_Output/all_Visvader_breast_tumors_Total_Singlets_Clusters_ident.tiff", plot = last_plot(), device = "tiff",
       scale = 1, width = 16, height = 10,
       dpi = 200, limitsize = TRUE)


########################################################################
##### Run SCSA in Python & then annotate clusters using the results ####
########################################################################
##################################################
######## Run in Terminal #########################
##################################################
#### Install SCSA to Home Directory, and create folders for "Visvader_Markers" and "Visvader_Annotations" within it
#### Place cluster markers file into "Visvader_Markers" and change "avg_log2FC" column to "avg_logFC"
#### Run the following within the terminal:
#### Change working directory
# cd ~/SCSA/

#### Run SCSA
# python3 SCSA.py -d whole_v2.db -i ~/SCSA/Visvader_Markers/all_Visvader_breast_tumors_Total_Singlets.cluster.markers.csv -s seurat -k All -E -g Human -p 0.01 -f 1.5 -o ~/SCSA/Visvader_Annotations/all_Visvader_breast_tumors_Total_Singlets.cluster.markers.txt -m txt

### Use SCSA results and known marker genes to facilitate annotations

############################################
######## Back to R #########################
############################################
# all_Visvader_breast_tumors_Total_Singlets <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Merged/all_Visvader_breast_tumors_Total_Singlets.rds")

# Subpopulations will be annotated here next
# Draw Plots
DimPlot(all_Visvader_breast_tumors_Total_Singlets, reduction = "umap", raster = FALSE)
FeaturePlot(all_Visvader_breast_tumors_Total_Singlets, c("KRT5", "KRT14", "KRT8", "KRT18", "KRT19", 
                                                         "VIM", "PTPRC", "CD3D", "APOE", "MS4A1", "PECAM1", "COL1A1"),  raster = FALSE)
# Annotate Clusters
all_Visvader_breast_tumors_Total_Singlets_IDs <- c("T & NK Cells", "Monocytes & Macrophages", "Epithelial/Neoplastic Cells", "Epithelial/Neoplastic Cells", "Fibroblasts",
                                                   "Epithelial/Neoplastic Cells", "B Cells", "Epithelial/Neoplastic Cells", "Epithelial/Neoplastic Cells", "Epithelial/Neoplastic Cells",
                                                   "Epithelial/Neoplastic Cells", "Epithelial/Neoplastic Cells", "Epithelial/Neoplastic Cells", "Epithelial/Neoplastic Cells", "Epithelial/Neoplastic Cells",
                                                   "Epithelial/Neoplastic Cells", "B Cells", "Epithelial/Neoplastic Cells", "Epithelial/Neoplastic Cells", "Pericytes",
                                                   "Endothelial Cells", "Epithelial/Neoplastic Cells", "Epithelial/Neoplastic Cells", "Epithelial/Neoplastic Cells", "Epithelial/Neoplastic Cells",
                                                   "Mast Cells", "Dendritic Cells")

names(all_Visvader_breast_tumors_Total_Singlets_IDs) <- levels(all_Visvader_breast_tumors_Total_Singlets)
all_Visvader_breast_tumors_Total_Singlets <- RenameIdents(all_Visvader_breast_tumors_Total_Singlets, all_Visvader_breast_tumors_Total_Singlets_IDs)

# Store Ident under Classification
all_Visvader_breast_tumors_Total_Singlets[["cell.ident"]] <- Idents(object = all_Visvader_breast_tumors_Total_Singlets)

# Draw plot
DimPlot(all_Visvader_breast_tumors_Total_Singlets, reduction = "umap", raster = FALSE)


# Draw plot with changed colors
DimPlot(all_Visvader_breast_tumors_Total_Singlets, reduction = "umap", raster = FALSE, cols = c('Epithelial/Neoplastic Cells' = '#F8766D',  'B Cells' = '#00A9FF', 'Dendritic Cells' = '#ABA300', 
                                                                                                'Endothelial Cells' = '#FF61CC','Fibroblasts' = '#8494FF', 'Mast Cells' = '#7CAE00', 
                                                                                                'Monocytes & Macrophages' = '#E68613', 'Pericytes' = '#ED68ED', 'T & NK Cells' = '#0CB702'))

# Save plots
ggsave("/R/R_Visvader/Visvader_Output/all_Visvader_breast_tumors_Total_Singlets_Annotated_UMAP_Custom_Colors.tiff", plot = last_plot(), device = "tiff", 
       scale = 1, width = 16, height = 10,
       dpi = 200, limitsize = TRUE)


FeaturePlot(all_Visvader_breast_tumors_Total_Singlets, c("KRT14"),  raster = FALSE, cols = c("grey90", "firebrick"))
ggsave("/R/R_Visvader/Visvader_Output/all_Visvader_breast_tumors_Total_Singlets_KRT14_UMAP_Custom_Colors.tiff", plot = last_plot(), device = "tiff", 
       scale = 1, width = 16, height = 10,
       dpi = 200, limitsize = TRUE)

FeaturePlot(all_Visvader_breast_tumors_Total_Singlets, c("KRT18"),  raster = FALSE, cols = c("grey90", "firebrick"))
ggsave("/R/R_Visvader/Visvader_Output/all_Visvader_breast_tumors_Total_Singlets_KRT18_UMAP_Custom_Colors.tiff", plot = last_plot(), device = "tiff", 
       scale = 1, width = 16, height = 10,
       dpi = 200, limitsize = TRUE)

FeaturePlot(all_Visvader_breast_tumors_Total_Singlets, c("KRT19"),  raster = FALSE, cols = c("grey90", "firebrick"))
ggsave("/R/R_Visvader/Visvader_Output/all_Visvader_breast_tumors_Total_Singlets_KRT19_UMAP_Custom_Colors.tiff", plot = last_plot(), device = "tiff", 
       scale = 1, width = 16, height = 10,
       dpi = 200, limitsize = TRUE)

# Save RDS
saveRDS(all_Visvader_breast_tumors_Total_Singlets, file = "/R/R_Visvader/Visvader_RDS/RDS_Annotated/all_Visvader_breast_tumors_Total_Singlets_Annotated.rds")



####################################
# Step 4: Extract Epithelial Cells # 
####################################
#all_Visvader_breast_tumors_Total_Singlets <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Annotated/all_Visvader_breast_tumors_Total_Singlets_Annotated.rds")

# Create ident order for plots
singlet_ident <- c("Epithelial/Neoplastic Cells", "B Cells", "Dendritic Cells", "Endothelial Cells", "Fibroblasts", "Mast Cells", "Monocytes & Macrophages", "Pericytes", "T & NK Cells")
all_Visvader_breast_tumors_Total_Singlets@active.ident <- factor(x = all_Visvader_breast_tumors_Total_Singlets@active.ident, levels = singlet_ident)
rm(singlet_ident)

# Visualize Cutoffs & Save Example Plots
VlnPlot(all_Visvader_breast_tumors_Total_Singlets, c("KRT14"), pt.size = 0.01, raster = FALSE)
ggsave("/R/R_Visvader/Visvader_Output/all_Visvader_breast_tumors_Total_Singlets_VLN_KRT14.tiff", plot = last_plot(), device = "tiff",
       scale = 1, width = 16, height = 10,
       dpi = 200, limitsize = TRUE)

VlnPlot(all_Visvader_breast_tumors_Total_Singlets, c("KRT18"), pt.size = 0.01, raster = FALSE)
# Save Example Plots
ggsave("/R/R_Visvader/Visvader_Output/all_Visvader_breast_tumors_Total_Singlets_VLN_KRT18.tiff", plot = last_plot(), device = "tiff",
       scale = 1, width = 16, height = 10,
       dpi = 200, limitsize = TRUE)

VlnPlot(all_Visvader_breast_tumors_Total_Singlets, c("KRT19"), pt.size = 0.01, raster = FALSE)
# Save Example Plots
ggsave("/R/R_Visvader/Visvader_Output/all_Visvader_breast_tumors_Total_Singlets_VLN_KRT19.tiff", plot = last_plot(), device = "tiff",
       scale = 1, width = 16, height = 10,
       dpi = 200, limitsize = TRUE)

# Visualize Cutoffs & Save Example Plots - No Points
VlnPlot(all_Visvader_breast_tumors_Total_Singlets, c("KRT14"), pt.size = 0, raster = FALSE)
ggsave("/R/R_Visvader/Visvader_Output/all_Visvader_breast_tumors_Total_Singlets_VLN_KRT14_No_Points.tiff", plot = last_plot(), device = "tiff",
       scale = 1, width = 16, height = 10,
       dpi = 200, limitsize = TRUE)

VlnPlot(all_Visvader_breast_tumors_Total_Singlets, c("KRT18"), pt.size = 0, raster = FALSE)
# Save Example Plots
ggsave("/R/R_Visvader/Visvader_Output/all_Visvader_breast_tumors_Total_Singlets_VLN_KRT18_No_Points.tiff", plot = last_plot(), device = "tiff",
       scale = 1, width = 16, height = 10,
       dpi = 200, limitsize = TRUE)

VlnPlot(all_Visvader_breast_tumors_Total_Singlets, c("KRT19"), pt.size = 0, raster = FALSE)
# Save Example Plots
ggsave("/R/R_Visvader/Visvader_Output/all_Visvader_breast_tumors_Total_Singlets_VLN_KRT19_No_Points.tiff", plot = last_plot(), device = "tiff",
       scale = 1, width = 16, height = 10,
       dpi = 200, limitsize = TRUE)

# KRT14, KRT18, KRT19
VlnPlot(all_Visvader_breast_tumors_Total_Singlets, c("KRT14", "KRT18", "KRT19"), ncol = 3, pt.size = 0.1)
ggsave("/R/R_Visvader/Visvader_Output/all_Visvader_breast_tumors_Total_Singlets_VLN_All_KRT_VLN.tiff", plot = last_plot(), device = "tiff",
       scale = 1, width = 20, height = 16,
       dpi = 200, limitsize = TRUE)

# KRT14, KRT18, KRT19 - No Points 
VlnPlot(all_Visvader_breast_tumors_Total_Singlets, c("KRT14", "KRT18", "KRT19"), ncol = 3, pt.size = 0)
ggsave("/R/R_Visvader/Visvader_Output/all_Visvader_breast_tumors_Total_Singlets_VLN_All_KRT_VLN_No_Points.tiff", plot = last_plot(), device = "tiff",
       scale = 1, width = 20, height = 16,
       dpi = 200, limitsize = TRUE)


# KRT14, KRT18, KRT19 - Transparent Points 
# Create violin plots with geom_jitter for each feature
# Extra margin added to left-most side
features <- c("KRT14", "KRT18", "KRT19")
plots <- lapply(features, function(feature) {
  VlnPlot(all_Visvader_breast_tumors_Total_Singlets, features = feature, pt.size = 0, raster = FALSE) +
    geom_jitter(width = 0.35, alpha = 0.1, size = 0.1) + 
    theme(legend.position = "none", plot.margin = margin(5, 5, 5, 50))
})

# Combine the plots into one
wrap_plots(plots, ncol = 3)
ggsave("/R/R_Visvader/Visvader_Output/all_Visvader_breast_tumors_Total_Singlets_VLN_All_KRT_VLN_Transparent_Points.tiff", plot = last_plot(), device = "tiff",
       scale = 1, width = 20, height = 16,
       dpi = 200, limitsize = TRUE)

rm(features)
rm(plots)

###############################################
########### Subset Epithelial Cells ###########
###############################################

# Select cells greater than cutoffs for at least one of the three main mammary epithelial cell types
all_Visvader_breast_tumors_panCK <- subset(x = all_Visvader_breast_tumors_Total_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
all_Visvader_breast_tumors_panCK <- FindVariableFeatures(all_Visvader_breast_tumors_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(all_Visvader_breast_tumors_panCK)
all_Visvader_breast_tumors_panCK <- ScaleData(all_Visvader_breast_tumors_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
all_Visvader_breast_tumors_panCK <- RunPCA(all_Visvader_breast_tumors_panCK, features = VariableFeatures(object = all_Visvader_breast_tumors_panCK))

# Examine and visualize PCA results a few different ways
print(all_Visvader_breast_tumors_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(all_Visvader_breast_tumors_panCK, dims = 1:2, reduction = "pca")
DimPlot(all_Visvader_breast_tumors_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(all_Visvader_breast_tumors_panCK, ndims = 50)

all_Visvader_breast_tumors_panCK <- FindNeighbors(all_Visvader_breast_tumors_panCK, dims = 1:40)
all_Visvader_breast_tumors_panCK <- FindClusters(all_Visvader_breast_tumors_panCK, resolution = 0.15)

# Umap clustering
all_Visvader_breast_tumors_panCK <- RunUMAP(all_Visvader_breast_tumors_panCK, dims = 1:40)
DimPlot(all_Visvader_breast_tumors_panCK, reduction = "umap", raster = FALSE)
DimPlot(all_Visvader_breast_tumors_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(all_Visvader_breast_tumors_panCK, file = "/R/R_Visvader/Visvader_RDS/RDS_Merged/all_Visvader_breast_tumors_panCK.rds")


########################################################################
# Extract SeuratObject by original identity and save as epithelial RDS #
########################################################################
# HR-pos
Visvader_0001_ER_panCK <- subset(x = all_Visvader_breast_tumors_panCK, subset = orig.ident == "Visvader_0001_ER_Total")
saveRDS(Visvader_0001_ER_panCK, file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0001_ER_panCK.rds")

Visvader_0025_ER_panCK <- subset(x = all_Visvader_breast_tumors_panCK, subset = orig.ident == "Visvader_0025_ER_Total")
saveRDS(Visvader_0025_ER_panCK, file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0025_ER_panCK.rds")

Visvader_0029_7C_ER_panCK <- subset(x = all_Visvader_breast_tumors_panCK, subset = orig.ident == "Visvader_0029_7C_ER_Total")
saveRDS(Visvader_0029_7C_ER_panCK, file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0029_7C_ER_panCK.rds")

Visvader_0029_9C_ER_panCK <- subset(x = all_Visvader_breast_tumors_panCK, subset = orig.ident == "Visvader_0029_9C_ER_Total")
saveRDS(Visvader_0029_9C_ER_panCK, file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0029_9C_ER_panCK.rds")

Visvader_0032_ER_panCK <- subset(x = all_Visvader_breast_tumors_panCK, subset = orig.ident == "Visvader_0032_ER_Total")
saveRDS(Visvader_0032_ER_panCK, file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0032_ER_panCK.rds")

Visvader_0040_ER_panCK <- subset(x = all_Visvader_breast_tumors_panCK, subset = orig.ident == "Visvader_0040_ER_Total")
saveRDS(Visvader_0040_ER_panCK, file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0040_ER_panCK.rds")

Visvader_0042_ER_panCK <- subset(x = all_Visvader_breast_tumors_panCK, subset = orig.ident == "Visvader_0042_ER_Total")
saveRDS(Visvader_0042_ER_panCK, file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0042_ER_panCK.rds")

Visvader_0043_ER_panCK <- subset(x = all_Visvader_breast_tumors_panCK, subset = orig.ident == "Visvader_0043_ER_Total")
saveRDS(Visvader_0043_ER_panCK, file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0043_ER_panCK.rds")

Visvader_0056_ER_panCK <- subset(x = all_Visvader_breast_tumors_panCK, subset = orig.ident == "Visvader_0056_ER_Total")
saveRDS(Visvader_0056_ER_panCK, file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0056_ER_panCK.rds")

Visvader_0064_ER_panCK <- subset(x = all_Visvader_breast_tumors_panCK, subset = orig.ident == "Visvader_0064_ER_Total")
saveRDS(Visvader_0064_ER_panCK, file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0064_ER_panCK.rds")

Visvader_0068_ER_panCK <- subset(x = all_Visvader_breast_tumors_panCK, subset = orig.ident == "Visvader_0068_ER_Total")
saveRDS(Visvader_0068_ER_panCK, file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0068_ER_panCK.rds")

Visvader_0114_ER_panCK <- subset(x = all_Visvader_breast_tumors_panCK, subset = orig.ident == "Visvader_0114_ER_Total")
saveRDS(Visvader_0114_ER_panCK, file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0114_ER_panCK.rds")

Visvader_0125_ER_panCK <- subset(x = all_Visvader_breast_tumors_panCK, subset = orig.ident == "Visvader_0125_ER_Total")
saveRDS(Visvader_0125_ER_panCK, file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0125_ER_panCK.rds")

Visvader_0151_ER_panCK <- subset(x = all_Visvader_breast_tumors_panCK, subset = orig.ident == "Visvader_0151_ER_Total")
saveRDS(Visvader_0151_ER_panCK, file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0151_ER_panCK.rds")

Visvader_0163_ER_panCK <- subset(x = all_Visvader_breast_tumors_panCK, subset = orig.ident == "Visvader_0163_ER_Total")
saveRDS(Visvader_0163_ER_panCK, file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0163_ER_panCK.rds")

Visvader_0167_ER_panCK <- subset(x = all_Visvader_breast_tumors_panCK, subset = orig.ident == "Visvader_0167_ER_Total")
saveRDS(Visvader_0167_ER_panCK, file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0167_ER_panCK.rds")

Visvader_0173_ER_panCK <- subset(x = all_Visvader_breast_tumors_panCK, subset = orig.ident == "Visvader_0173_ER_Total")
saveRDS(Visvader_0173_ER_panCK, file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0173_ER_panCK.rds")

Visvader_0178_ER_panCK <- subset(x = all_Visvader_breast_tumors_panCK, subset = orig.ident == "Visvader_0178_ER_Total")
saveRDS(Visvader_0178_ER_panCK, file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0178_ER_panCK.rds")

Visvader_0360_ER_panCK <- subset(x = all_Visvader_breast_tumors_panCK, subset = orig.ident == "Visvader_0360_ER_Total")
saveRDS(Visvader_0360_ER_panCK, file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0360_ER_panCK.rds")

Visvader_0319_PR_panCK <- subset(x = all_Visvader_breast_tumors_panCK, subset = orig.ident == "Visvader_0319_PR_Total")
saveRDS(Visvader_0319_PR_panCK, file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0319_PR_panCK.rds")

# HER2-pos
Visvader_0031_HER2_panCK <- subset(x = all_Visvader_breast_tumors_panCK, subset = orig.ident == "Visvader_0031_HER2_Total")
saveRDS(Visvader_0031_HER2_panCK, file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0031_HER2_panCK.rds")

Visvader_0069_HER2_panCK <- subset(x = all_Visvader_breast_tumors_panCK, subset = orig.ident == "Visvader_0069_HER2_Total")
saveRDS(Visvader_0069_HER2_panCK, file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0069_HER2_panCK.rds")

Visvader_0161_HER2_panCK <- subset(x = all_Visvader_breast_tumors_panCK, subset = orig.ident == "Visvader_0161_HER2_Total")
saveRDS(Visvader_0161_HER2_panCK, file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0161_HER2_panCK.rds")

Visvader_0176_HER2_panCK <- subset(x = all_Visvader_breast_tumors_panCK, subset = orig.ident == "Visvader_0176_HER2_Total")
saveRDS(Visvader_0176_HER2_panCK, file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0176_HER2_panCK.rds")

Visvader_0308_HER2_panCK <- subset(x = all_Visvader_breast_tumors_panCK, subset = orig.ident == "Visvader_0308_HER2_Total")
saveRDS(Visvader_0308_HER2_panCK, file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0308_HER2_panCK.rds")

Visvader_0337_HER2_panCK <- subset(x = all_Visvader_breast_tumors_panCK, subset = orig.ident == "Visvader_0337_HER2_Total")
saveRDS(Visvader_0337_HER2_panCK, file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0337_HER2_panCK.rds")

# TNBC
Visvader_106_TNBC_panCK <- subset(x = all_Visvader_breast_tumors_panCK, subset = orig.ident == "Visvader_106_TNBC_Total")
saveRDS(Visvader_106_TNBC_panCK, file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_106_TNBC_panCK.rds")

Visvader_114_TNBC_panCK <- subset(x = all_Visvader_breast_tumors_panCK, subset = orig.ident == "Visvader_114_TNBC_Total")
saveRDS(Visvader_114_TNBC_panCK, file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_114_TNBC_panCK.rds")

Visvader_126_TNBC_panCK <- subset(x = all_Visvader_breast_tumors_panCK, subset = orig.ident == "Visvader_126_TNBC_Total")
saveRDS(Visvader_126_TNBC_panCK, file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_126_TNBC_panCK.rds")

Visvader_0131_TNBC_panCK <- subset(x = all_Visvader_breast_tumors_panCK, subset = orig.ident == "Visvader_0131_TNBC_Total")
saveRDS(Visvader_0131_TNBC_panCK, file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0131_TNBC_panCK.rds")

Visvader_135_TNBC_panCK <- subset(x = all_Visvader_breast_tumors_panCK, subset = orig.ident == "Visvader_135_TNBC_Total")
saveRDS(Visvader_135_TNBC_panCK, file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_135_TNBC_panCK.rds")

Visvader_0177_TNBC_panCK <- subset(x = all_Visvader_breast_tumors_panCK, subset = orig.ident == "Visvader_0177_TNBC_Total")
saveRDS(Visvader_0177_TNBC_panCK, file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0177_TNBC_panCK.rds")

Visvader_0554_TNBC_panCK <- subset(x = all_Visvader_breast_tumors_panCK, subset = orig.ident == "Visvader_0554_TNBC_Total")
saveRDS(Visvader_0554_TNBC_panCK, file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0554_TNBC_panCK.rds")

Visvader_4031_TNBC_panCK <- subset(x = all_Visvader_breast_tumors_panCK, subset = orig.ident == "Visvader_4031_TNBC_Total")
saveRDS(Visvader_4031_TNBC_panCK, file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_4031_TNBC_panCK.rds")

# Remove individual objects
rm(all_Visvader_breast_tumors_Total_Singlets)
rm(all_Visvader_breast_tumors_panCK)
rm(Visvader_0001_ER_panCK)
rm(Visvader_0025_ER_panCK)
rm(Visvader_0029_7C_ER_panCK)
rm(Visvader_0029_9C_ER_panCK)
rm(Visvader_0032_ER_panCK)
rm(Visvader_0040_ER_panCK)
rm(Visvader_0042_ER_panCK)
rm(Visvader_0043_ER_panCK)
rm(Visvader_0056_ER_panCK)
rm(Visvader_0064_ER_panCK)
rm(Visvader_0068_ER_panCK)
rm(Visvader_0114_ER_panCK)
rm(Visvader_0125_ER_panCK)
rm(Visvader_0151_ER_panCK)
rm(Visvader_0163_ER_panCK)
rm(Visvader_0167_ER_panCK)
rm(Visvader_0173_ER_panCK)
rm(Visvader_0178_ER_panCK)
rm(Visvader_0360_ER_panCK)
rm(Visvader_0319_PR_panCK)
rm(Visvader_0031_HER2_panCK)
rm(Visvader_0069_HER2_panCK)
rm(Visvader_0161_HER2_panCK)
rm(Visvader_0176_HER2_panCK)
rm(Visvader_0308_HER2_panCK)
rm(Visvader_0337_HER2_panCK)
rm(Visvader_106_TNBC_panCK)
rm(Visvader_114_TNBC_panCK)
rm(Visvader_126_TNBC_panCK)
rm(Visvader_0131_TNBC_panCK)
rm(Visvader_135_TNBC_panCK)
rm(Visvader_0177_TNBC_panCK)
rm(Visvader_0554_TNBC_panCK)
rm(Visvader_4031_TNBC_panCK)
gc()



