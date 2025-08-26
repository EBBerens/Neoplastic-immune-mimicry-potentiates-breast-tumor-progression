#####################################################################################################################
#                   Visvader (Pal et al) Dataset Reduction Mammoplasty Analysis Steps                               #
#####################################################################################################################
# Step 3: Merge Reduction Mammoplasties into Combined Seurat Object                                                 #
# Step 4: Extract Epithelial Cells & Merge into Combined Seurat Object                                              #
#####################################################################################################################


##################################################################################################
# Step 3: Merge Reduction Mammoplasties into Combined Seurat Object                              # 
##################################################################################################
# Load Libraries
library(Seurat)
library(patchwork)
library(dplyr)
library(ggplot2)


#########################
# Load Individual Files #
#########################
Visvader_0019_NORM_Total_Singlets <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Total_Singlets/Visvader_0019_NORM_Total_Singlets.rds")
Visvader_0021_NORM_Total_Singlets <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Total_Singlets/Visvader_0021_NORM_Total_Singlets.rds")
Visvader_0064_NORM_Total_Singlets <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Total_Singlets/Visvader_0064_NORM_Total_Singlets.rds")
Visvader_0092_NORM_Total_Singlets <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Total_Singlets/Visvader_0092_NORM_Total_Singlets.rds")
Visvader_0093_NORM_Total_Singlets <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Total_Singlets/Visvader_0093_NORM_Total_Singlets.rds")
Visvader_0123_NORM_Total_Singlets <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Total_Singlets/Visvader_0123_NORM_Total_Singlets.rds")
Visvader_0169_NORM_Total_Singlets <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Total_Singlets/Visvader_0169_NORM_Total_Singlets.rds")
Visvader_0230_NORM_Total_Singlets <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Total_Singlets/Visvader_0230_NORM_Total_Singlets.rds")
Visvader_0233_NORM_Total_Singlets <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Total_Singlets/Visvader_0233_NORM_Total_Singlets.rds")
Visvader_0275_NORM_Total_Singlets <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Total_Singlets/Visvader_0275_NORM_Total_Singlets.rds")
Visvader_0288_NORM_Total_Singlets <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Total_Singlets/Visvader_0288_NORM_Total_Singlets.rds")
Visvader_0342_NORM_Total_Singlets <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Total_Singlets/Visvader_0342_NORM_Total_Singlets.rds")
Visvader_0372_NORM_Total_Singlets <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Total_Singlets/Visvader_0372_NORM_Total_Singlets.rds")


##############################################
# Merge Individual Files into Single Dataset #
##############################################
# Merge Datasets
Visvader_NORM_Total_Singlets <- merge(Visvader_0019_NORM_Total_Singlets, y = c(Visvader_0021_NORM_Total_Singlets, Visvader_0064_NORM_Total_Singlets, Visvader_0092_NORM_Total_Singlets, Visvader_0093_NORM_Total_Singlets,
                                                                               Visvader_0123_NORM_Total_Singlets, Visvader_0169_NORM_Total_Singlets, Visvader_0230_NORM_Total_Singlets, Visvader_0233_NORM_Total_Singlets,
                                                                               Visvader_0275_NORM_Total_Singlets, Visvader_0288_NORM_Total_Singlets, Visvader_0342_NORM_Total_Singlets, Visvader_0372_NORM_Total_Singlets))

# Join Layers
Visvader_NORM_Total_Singlets <- JoinLayers(Visvader_NORM_Total_Singlets)


# Remove individual objects
rm(Visvader_0019_NORM_Total_Singlets)
rm(Visvader_0021_NORM_Total_Singlets)
rm(Visvader_0064_NORM_Total_Singlets)
rm(Visvader_0092_NORM_Total_Singlets)
rm(Visvader_0093_NORM_Total_Singlets)
rm(Visvader_0123_NORM_Total_Singlets)
rm(Visvader_0169_NORM_Total_Singlets)
rm(Visvader_0230_NORM_Total_Singlets)
rm(Visvader_0233_NORM_Total_Singlets)
rm(Visvader_0275_NORM_Total_Singlets)
rm(Visvader_0288_NORM_Total_Singlets)
rm(Visvader_0342_NORM_Total_Singlets)
rm(Visvader_0372_NORM_Total_Singlets)


####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Visvader_NORM_Total_Singlets <- FindVariableFeatures(Visvader_NORM_Total_Singlets, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Visvader_NORM_Total_Singlets)
Visvader_NORM_Total_Singlets <- ScaleData(Visvader_NORM_Total_Singlets, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Visvader_NORM_Total_Singlets <- RunPCA(Visvader_NORM_Total_Singlets, features = VariableFeatures(object = Visvader_NORM_Total_Singlets))

# Examine and visualize PCA results a few different ways
print(Visvader_NORM_Total_Singlets[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Visvader_NORM_Total_Singlets, dims = 1:2, reduction = "pca")
DimPlot(Visvader_NORM_Total_Singlets, reduction = "pca")

# Visualize PCs
ElbowPlot(Visvader_NORM_Total_Singlets, ndims = 50)

# Choose dimensions
Visvader_NORM_Total_Singlets <- FindNeighbors(Visvader_NORM_Total_Singlets, dims = 1:30)
Visvader_NORM_Total_Singlets <- FindClusters(Visvader_NORM_Total_Singlets, resolution = 0.15)

# Umap clustering
Visvader_NORM_Total_Singlets <- RunUMAP(Visvader_NORM_Total_Singlets, dims = 1:30)
DimPlot(Visvader_NORM_Total_Singlets, reduction = "umap")

# Identify Samples
DimPlot(Visvader_NORM_Total_Singlets, reduction = "umap", group.by = "orig.ident")

# Save
saveRDS(Visvader_NORM_Total_Singlets, file = "/R/R_Visvader/Visvader_RDS/RDS_Merged/Visvader_NORM_Total_Singlets.rds")






###############################################################################
# Step 4: Extract Epithelial Cells & Merge into Combined Seurat Object        #
###############################################################################
# Read RDS
#Visvader_NORM_Total_Singlets <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Merged/Visvader_NORM_Total_Singlets.rds")
  
# Visualize Cutoffs & Save Example Plots
VlnPlot(Visvader_NORM_Total_Singlets, c("KRT14"), pt.size = 0.1, raster = FALSE)
ggsave("/R/R_Visvader/Visvader_Output/Visvader_NORM_Total_Singlets_VLN_KRT14.tiff", plot = last_plot(), device = "tiff",
       scale = 1, width = 16, height = 10,
       dpi = 200, limitsize = TRUE)

VlnPlot(Visvader_NORM_Total_Singlets, c("KRT18"), pt.size = 0.1, raster = FALSE)
# Save Example Plots
ggsave("/R/R_Visvader/Visvader_Output/Visvader_NORM_Total_Singlets_VLN_KRT18.tiff", plot = last_plot(), device = "tiff",
       scale = 1, width = 16, height = 10,
       dpi = 200, limitsize = TRUE)

VlnPlot(Visvader_NORM_Total_Singlets, c("KRT19"), pt.size = 0.1, raster = FALSE)
# Save Example Plots
ggsave("/R/R_Visvader/Visvader_Output/Visvader_NORM_Total_Singlets_VLN_KRT19.tiff", plot = last_plot(), device = "tiff",
       scale = 1, width = 16, height = 10,
       dpi = 200, limitsize = TRUE)

# Visualize Cutoffs & Save Example Plots - No Points
VlnPlot(Visvader_NORM_Total_Singlets, c("KRT14"), pt.size = 0, raster = FALSE)
ggsave("/R/R_Visvader/Visvader_Output/Visvader_NORM_Total_Singlets_VLN_KRT14_No_Points.tiff", plot = last_plot(), device = "tiff",
       scale = 1, width = 16, height = 10,
       dpi = 200, limitsize = TRUE)

VlnPlot(Visvader_NORM_Total_Singlets, c("KRT18"), pt.size = 0, raster = FALSE)
# Save Example Plots
ggsave("/R/R_Visvader/Visvader_Output/Visvader_NORM_Total_Singlets_VLN_KRT18_No_Points.tiff", plot = last_plot(), device = "tiff",
       scale = 1, width = 16, height = 10,
       dpi = 200, limitsize = TRUE)

VlnPlot(Visvader_NORM_Total_Singlets, c("KRT19"), pt.size = 0, raster = FALSE)
# Save Example Plots
ggsave("/R/R_Visvader/Visvader_Output/Visvader_NORM_Total_Singlets_VLN_KRT19_No_Points.tiff", plot = last_plot(), device = "tiff",
       scale = 1, width = 16, height = 10,
       dpi = 200, limitsize = TRUE)


###############################################
########### Subset Epithelial Cells ###########
###############################################

# Cells greater than cutoffs for at least one of the three main mammary epithelial cell types
all_Visvader_NORM_panCK <- subset(x = Visvader_NORM_Total_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)
# Remove object
rm(Visvader_NORM_Total_Singlets)


####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
all_Visvader_NORM_panCK <- FindVariableFeatures(all_Visvader_NORM_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(all_Visvader_NORM_panCK)
all_Visvader_NORM_panCK <- ScaleData(all_Visvader_NORM_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
all_Visvader_NORM_panCK <- RunPCA(all_Visvader_NORM_panCK, features = VariableFeatures(object = all_Visvader_NORM_panCK))

# Examine and visualize PCA results a few different ways
print(all_Visvader_NORM_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(all_Visvader_NORM_panCK, dims = 1:2, reduction = "pca")
DimPlot(all_Visvader_NORM_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(all_Visvader_NORM_panCK, ndims = 60)

all_Visvader_NORM_panCK <- FindNeighbors(all_Visvader_NORM_panCK, dims = 1:40)
all_Visvader_NORM_panCK <- FindClusters(all_Visvader_NORM_panCK, resolution = 0.15)

# Umap clustering
all_Visvader_NORM_panCK <- RunUMAP(all_Visvader_NORM_panCK, dims = 1:40)
DimPlot(all_Visvader_NORM_panCK, reduction = "umap", raster = FALSE)
DimPlot(all_Visvader_NORM_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(all_Visvader_NORM_panCK, file = "/R/R_Visvader/Visvader_RDS/RDS_Merged/all_Visvader_NORM_panCK.rds")

########################################################################
# Extract SeuratObject by original identity and save as epithelial RDS #
########################################################################
# Normal Samples
Visvader_0019_NORM_panCK <- subset(x = all_Visvader_NORM_panCK, subset = orig.ident == "Visvader_0019_NORM_Total")
saveRDS(Visvader_0019_NORM_panCK, file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0019_NORM_panCK.rds")

Visvader_0021_NORM_panCK <- subset(x = all_Visvader_NORM_panCK, subset = orig.ident == "Visvader_0021_NORM_Total")
saveRDS(Visvader_0021_NORM_panCK, file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0021_NORM_panCK.rds")

Visvader_0064_NORM_panCK <- subset(x = all_Visvader_NORM_panCK, subset = orig.ident == "Visvader_0064_NORM_Total")
saveRDS(Visvader_0064_NORM_panCK, file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0064_NORM_panCK.rds")

Visvader_0092_NORM_panCK <- subset(x = all_Visvader_NORM_panCK, subset = orig.ident == "Visvader_0092_NORM_Total")
saveRDS(Visvader_0092_NORM_panCK, file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0092_NORM_panCK.rds")

Visvader_0092_NORM_panCK <- subset(x = all_Visvader_NORM_panCK, subset = orig.ident == "Visvader_0092_NORM_Total")
saveRDS(Visvader_0092_NORM_panCK, file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0092_NORM_panCK.rds")

Visvader_0093_NORM_panCK <- subset(x = all_Visvader_NORM_panCK, subset = orig.ident == "Visvader_0093_NORM_Total")
saveRDS(Visvader_0093_NORM_panCK, file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0093_NORM_panCK.rds")

Visvader_0123_NORM_panCK <- subset(x = all_Visvader_NORM_panCK, subset = orig.ident == "Visvader_0123_NORM_Total")
saveRDS(Visvader_0123_NORM_panCK, file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0123_NORM_panCK.rds")

Visvader_0169_NORM_panCK <- subset(x = all_Visvader_NORM_panCK, subset = orig.ident == "Visvader_0169_NORM_Total")
saveRDS(Visvader_0169_NORM_panCK, file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0169_NORM_panCK.rds")

Visvader_0230_NORM_panCK <- subset(x = all_Visvader_NORM_panCK, subset = orig.ident == "Visvader_0230_NORM_Total")
saveRDS(Visvader_0230_NORM_panCK, file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0230_NORM_panCK.rds")

Visvader_0233_NORM_panCK <- subset(x = all_Visvader_NORM_panCK, subset = orig.ident == "Visvader_0233_NORM_Total")
saveRDS(Visvader_0233_NORM_panCK, file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0233_NORM_panCK.rds")

Visvader_0275_NORM_panCK <- subset(x = all_Visvader_NORM_panCK, subset = orig.ident == "Visvader_0275_NORM_Total")
saveRDS(Visvader_0275_NORM_panCK, file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0275_NORM_panCK.rds")

Visvader_0275_NORM_panCK <- subset(x = all_Visvader_NORM_panCK, subset = orig.ident == "Visvader_0275_NORM_Total")
saveRDS(Visvader_0275_NORM_panCK, file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0275_NORM_panCK.rds")

Visvader_0288_NORM_panCK <- subset(x = all_Visvader_NORM_panCK, subset = orig.ident == "Visvader_0288_NORM_Total")
saveRDS(Visvader_0288_NORM_panCK, file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0288_NORM_panCK.rds")

Visvader_0342_NORM_panCK <- subset(x = all_Visvader_NORM_panCK, subset = orig.ident == "Visvader_0342_NORM_Total")
saveRDS(Visvader_0342_NORM_panCK, file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0342_NORM_panCK.rds")

Visvader_0372_NORM_panCK <- subset(x = all_Visvader_NORM_panCK, subset = orig.ident == "Visvader_0372_NORM_Total")
saveRDS(Visvader_0372_NORM_panCK, file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0372_NORM_panCK.rds")

# Remove individual objects
rm(all_Visvader_NORM_panCK)
rm(Visvader_0019_NORM_panCK)
rm(Visvader_0021_NORM_panCK)
rm(Visvader_0064_NORM_panCK)
rm(Visvader_0092_NORM_panCK)
rm(Visvader_0093_NORM_panCK)
rm(Visvader_0123_NORM_panCK)
rm(Visvader_0169_NORM_panCK)
rm(Visvader_0230_NORM_panCK)
rm(Visvader_0233_NORM_panCK)
rm(Visvader_0275_NORM_panCK)
rm(Visvader_0288_NORM_panCK)
rm(Visvader_0342_NORM_panCK)
rm(Visvader_0372_NORM_panCK)
gc()



######################################################
# Merge Epithelial Cells into Combined Seurat Object #
######################################################

#########################
# Load Individual Files #
#########################
Visvader_0019_NORM_panCK <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0019_NORM_panCK.rds")
Visvader_0021_NORM_panCK <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0021_NORM_panCK.rds")
Visvader_0064_NORM_panCK <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0064_NORM_panCK.rds")
Visvader_0092_NORM_panCK <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0092_NORM_panCK.rds")
Visvader_0093_NORM_panCK <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0093_NORM_panCK.rds")
Visvader_0123_NORM_panCK <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0123_NORM_panCK.rds")
Visvader_0169_NORM_panCK <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0169_NORM_panCK.rds")
Visvader_0230_NORM_panCK <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0230_NORM_panCK.rds")
Visvader_0233_NORM_panCK <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0233_NORM_panCK.rds")
Visvader_0275_NORM_panCK <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0275_NORM_panCK.rds")
Visvader_0288_NORM_panCK <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0288_NORM_panCK.rds")
Visvader_0342_NORM_panCK <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0342_NORM_panCK.rds")
Visvader_0372_NORM_panCK <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0372_NORM_panCK.rds")

##############################################
# Merge Individual Files into Single Dataset #
##############################################
# Merge Datasets
Visvader_NORM_panCK_merged <- merge(Visvader_0019_NORM_panCK, y = c(Visvader_0021_NORM_panCK, Visvader_0064_NORM_panCK, Visvader_0092_NORM_panCK, Visvader_0093_NORM_panCK,
                                                                            Visvader_0123_NORM_panCK, Visvader_0169_NORM_panCK, Visvader_0230_NORM_panCK, Visvader_0233_NORM_panCK,
                                                                            Visvader_0275_NORM_panCK, Visvader_0288_NORM_panCK, Visvader_0342_NORM_panCK, Visvader_0372_NORM_panCK))

# Join Layers
Visvader_NORM_panCK_merged <- JoinLayers(Visvader_NORM_panCK_merged)

# Remove individual objects
rm(Visvader_0019_NORM_panCK)
rm(Visvader_0021_NORM_panCK)
rm(Visvader_0064_NORM_panCK)
rm(Visvader_0092_NORM_panCK)
rm(Visvader_0093_NORM_panCK)
rm(Visvader_0123_NORM_panCK)
rm(Visvader_0169_NORM_panCK)
rm(Visvader_0230_NORM_panCK)
rm(Visvader_0233_NORM_panCK)
rm(Visvader_0275_NORM_panCK)
rm(Visvader_0288_NORM_panCK)
rm(Visvader_0342_NORM_panCK)
rm(Visvader_0372_NORM_panCK)


####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Visvader_NORM_panCK_merged <- FindVariableFeatures(Visvader_NORM_panCK_merged, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Visvader_NORM_panCK_merged)
Visvader_NORM_panCK_merged <- ScaleData(Visvader_NORM_panCK_merged, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Visvader_NORM_panCK_merged <- RunPCA(Visvader_NORM_panCK_merged, features = VariableFeatures(object = Visvader_NORM_panCK_merged))

# Examine and visualize PCA results a few different ways
print(Visvader_NORM_panCK_merged[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Visvader_NORM_panCK_merged, dims = 1:2, reduction = "pca")
DimPlot(Visvader_NORM_panCK_merged, reduction = "pca")

# Visualize PCs
ElbowPlot(Visvader_NORM_panCK_merged, ndims = 50)

# Choose dimensions
Visvader_NORM_panCK_merged <- FindNeighbors(Visvader_NORM_panCK_merged, dims = 1:40)
Visvader_NORM_panCK_merged <- FindClusters(Visvader_NORM_panCK_merged, resolution = 0.1)

# Umap clustering
Visvader_NORM_panCK_merged <- RunUMAP(Visvader_NORM_panCK_merged, dims = 1:40)
DimPlot(Visvader_NORM_panCK_merged, reduction = "umap")

# Identify Samples
DimPlot(Visvader_NORM_panCK_merged, reduction = "umap", group.by = "orig.ident")

# Save
saveRDS(Visvader_NORM_panCK_merged, file = "/R/R_Visvader/Visvader_RDS/RDS_Merged/Visvader_NORM_panCK_merged.rds")

# Remove Objects
rm(Visvader_NORM_panCK_merged)
gc()


