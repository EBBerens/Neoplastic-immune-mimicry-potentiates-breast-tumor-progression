#####################################################################################################################
#                Navin (Kumar et al) Dataset Reduction Mammoplasty Analysis Steps                                   #
########################################################################## ##########################################
# Step 1: Process Each Reduction Mammoplasty Count Matrix into Seurat Object                                        #
# Step 2: Run DoubletFinder and Subset Singlets                                                                     #
#####################################################################################################################


######################################################################
# Step 1: Process Each Reduction Mammoplasty h5 into Seurat Object   #
######################################################################

# Load Libraries 
library(Seurat)
library(patchwork)
library(dplyr)

# Data downloaded from GEO (http://ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE235326)

##########################################################
########################### Navin_hbca_c14 ###############
##########################################################

# Load the dataset
Navin_hbca_c14 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c14/GSM7500506_hbca_c14_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c14 <- CreateSeuratObject(counts = Navin_hbca_c14, project = "Navin_hbca_c14", min.cells = 3, min.features = 200)
Navin_hbca_c14

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c14[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c14, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c14, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c14, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c14, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c14 <- subset(Navin_hbca_c14, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c14 <- NormalizeData(Navin_hbca_c14, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c14 <- FindVariableFeatures(Navin_hbca_c14, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c14)
Navin_hbca_c14 <- ScaleData(Navin_hbca_c14, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c14 <- RunPCA(Navin_hbca_c14, features = VariableFeatures(object = Navin_hbca_c14))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c14[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c14, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c14, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c14 <- JackStraw(Navin_hbca_c14, num.replicate = 100)
Navin_hbca_c14 <- ScoreJackStraw(Navin_hbca_c14, dims = 1:20)
ElbowPlot(Navin_hbca_c14)

# Cluster cells
Navin_hbca_c14 <- FindNeighbors(Navin_hbca_c14, dims = 1:20)
Navin_hbca_c14 <- FindClusters(Navin_hbca_c14, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c14 <- RunUMAP(Navin_hbca_c14, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c14, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c14, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c14.rds")

# Remove Object
rm(Navin_hbca_c14)



##########################################################
########################### Navin_hbca_c15 ###############
##########################################################

# Load the dataset
Navin_hbca_c15 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c15/GSM7500507_hbca_c15_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c15 <- CreateSeuratObject(counts = Navin_hbca_c15, project = "Navin_hbca_c15", min.cells = 3, min.features = 200)
Navin_hbca_c15

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c15[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c15, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c15, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c15, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c15, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c15 <- subset(Navin_hbca_c15, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c15 <- NormalizeData(Navin_hbca_c15, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c15 <- FindVariableFeatures(Navin_hbca_c15, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c15)
Navin_hbca_c15 <- ScaleData(Navin_hbca_c15, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c15 <- RunPCA(Navin_hbca_c15, features = VariableFeatures(object = Navin_hbca_c15))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c15[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c15, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c15, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c15 <- JackStraw(Navin_hbca_c15, num.replicate = 100)
Navin_hbca_c15 <- ScoreJackStraw(Navin_hbca_c15, dims = 1:20)
ElbowPlot(Navin_hbca_c15)

# Cluster cells
Navin_hbca_c15 <- FindNeighbors(Navin_hbca_c15, dims = 1:20)
Navin_hbca_c15 <- FindClusters(Navin_hbca_c15, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c15 <- RunUMAP(Navin_hbca_c15, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c15, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c15, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c15.rds")

# Remove Object
rm(Navin_hbca_c15)



##########################################################
########################### Navin_hbca_c19 ###############
##########################################################

# Load the dataset
Navin_hbca_c19 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c19/GSM7500511_hbca_c19_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c19 <- CreateSeuratObject(counts = Navin_hbca_c19, project = "Navin_hbca_c19", min.cells = 3, min.features = 200)
Navin_hbca_c19

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c19[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c19, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c19, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c19, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c19, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c19 <- subset(Navin_hbca_c19, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c19 <- NormalizeData(Navin_hbca_c19, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c19 <- FindVariableFeatures(Navin_hbca_c19, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c19)
Navin_hbca_c19 <- ScaleData(Navin_hbca_c19, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c19 <- RunPCA(Navin_hbca_c19, features = VariableFeatures(object = Navin_hbca_c19))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c19[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c19, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c19, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c19 <- JackStraw(Navin_hbca_c19, num.replicate = 100)
Navin_hbca_c19 <- ScoreJackStraw(Navin_hbca_c19, dims = 1:20)
ElbowPlot(Navin_hbca_c19)

# Cluster cells
Navin_hbca_c19 <- FindNeighbors(Navin_hbca_c19, dims = 1:20)
Navin_hbca_c19 <- FindClusters(Navin_hbca_c19, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c19 <- RunUMAP(Navin_hbca_c19, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c19, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c19, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c19.rds")

# Remove Object
rm(Navin_hbca_c19)



##########################################################
########################### Navin_hbca_c20 ###############
##########################################################

# Load the dataset
Navin_hbca_c20 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c20/GSM7500505_hbca_c20_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c20 <- CreateSeuratObject(counts = Navin_hbca_c20, project = "Navin_hbca_c20", min.cells = 3, min.features = 200)
Navin_hbca_c20

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c20[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c20, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c20, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c20, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c20, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c20 <- subset(Navin_hbca_c20, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c20 <- NormalizeData(Navin_hbca_c20, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c20 <- FindVariableFeatures(Navin_hbca_c20, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c20)
Navin_hbca_c20 <- ScaleData(Navin_hbca_c20, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c20 <- RunPCA(Navin_hbca_c20, features = VariableFeatures(object = Navin_hbca_c20))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c20[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c20, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c20, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c20 <- JackStraw(Navin_hbca_c20, num.replicate = 100)
Navin_hbca_c20 <- ScoreJackStraw(Navin_hbca_c20, dims = 1:20)
ElbowPlot(Navin_hbca_c20)

# Cluster cells
Navin_hbca_c20 <- FindNeighbors(Navin_hbca_c20, dims = 1:20)
Navin_hbca_c20 <- FindClusters(Navin_hbca_c20, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c20 <- RunUMAP(Navin_hbca_c20, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c20, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c20, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c20.rds")

# Remove Object
rm(Navin_hbca_c20)



##########################################################
########################### Navin_hbca_c22 ###############
##########################################################

# Load the dataset
Navin_hbca_c22 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c22/GSM7500508_hbca_c22_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c22 <- CreateSeuratObject(counts = Navin_hbca_c22, project = "Navin_hbca_c22", min.cells = 3, min.features = 200)
Navin_hbca_c22

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c22[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c22, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c22, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c22, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c22, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c22 <- subset(Navin_hbca_c22, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c22 <- NormalizeData(Navin_hbca_c22, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c22 <- FindVariableFeatures(Navin_hbca_c22, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c22)
Navin_hbca_c22 <- ScaleData(Navin_hbca_c22, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c22 <- RunPCA(Navin_hbca_c22, features = VariableFeatures(object = Navin_hbca_c22))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c22[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c22, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c22, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c22 <- JackStraw(Navin_hbca_c22, num.replicate = 100)
Navin_hbca_c22 <- ScoreJackStraw(Navin_hbca_c22, dims = 1:20)
ElbowPlot(Navin_hbca_c22)

# Cluster cells
Navin_hbca_c22 <- FindNeighbors(Navin_hbca_c22, dims = 1:20)
Navin_hbca_c22 <- FindClusters(Navin_hbca_c22, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c22 <- RunUMAP(Navin_hbca_c22, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c22, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c22, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c22.rds")

# Remove Object
rm(Navin_hbca_c22)



##########################################################
########################### Navin_hbca_c23 ###############
##########################################################

# Load the dataset
Navin_hbca_c23 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c23/GSM7500509_hbca_c23_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c23 <- CreateSeuratObject(counts = Navin_hbca_c23, project = "Navin_hbca_c23", min.cells = 3, min.features = 200)
Navin_hbca_c23

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c23[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c23, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c23, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c23, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c23, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c23 <- subset(Navin_hbca_c23, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c23 <- NormalizeData(Navin_hbca_c23, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c23 <- FindVariableFeatures(Navin_hbca_c23, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c23)
Navin_hbca_c23 <- ScaleData(Navin_hbca_c23, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c23 <- RunPCA(Navin_hbca_c23, features = VariableFeatures(object = Navin_hbca_c23))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c23[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c23, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c23, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c23 <- JackStraw(Navin_hbca_c23, num.replicate = 100)
Navin_hbca_c23 <- ScoreJackStraw(Navin_hbca_c23, dims = 1:20)
ElbowPlot(Navin_hbca_c23)

# Cluster cells
Navin_hbca_c23 <- FindNeighbors(Navin_hbca_c23, dims = 1:20)
Navin_hbca_c23 <- FindClusters(Navin_hbca_c23, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c23 <- RunUMAP(Navin_hbca_c23, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c23, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c23, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c23.rds")

# Remove Object
rm(Navin_hbca_c23)



##########################################################
########################### Navin_hbca_c24 ###############
##########################################################

# Load the dataset
Navin_hbca_c24 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c24/GSM7500510_hbca_c24_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c24 <- CreateSeuratObject(counts = Navin_hbca_c24, project = "Navin_hbca_c24", min.cells = 3, min.features = 200)
Navin_hbca_c24

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c24[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c24, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c24, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c24, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c24, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c24 <- subset(Navin_hbca_c24, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c24 <- NormalizeData(Navin_hbca_c24, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c24 <- FindVariableFeatures(Navin_hbca_c24, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c24)
Navin_hbca_c24 <- ScaleData(Navin_hbca_c24, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c24 <- RunPCA(Navin_hbca_c24, features = VariableFeatures(object = Navin_hbca_c24))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c24[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c24, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c24, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c24 <- JackStraw(Navin_hbca_c24, num.replicate = 100)
Navin_hbca_c24 <- ScoreJackStraw(Navin_hbca_c24, dims = 1:20)
ElbowPlot(Navin_hbca_c24)

# Cluster cells
Navin_hbca_c24 <- FindNeighbors(Navin_hbca_c24, dims = 1:20)
Navin_hbca_c24 <- FindClusters(Navin_hbca_c24, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c24 <- RunUMAP(Navin_hbca_c24, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c24, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c24, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c24.rds")

# Remove Object
rm(Navin_hbca_c24)



##########################################################
########################### Navin_hbca_c25 ###############
##########################################################

# Load the dataset
Navin_hbca_c25 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c25/GSM7500484_hbca_c25_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c25 <- CreateSeuratObject(counts = Navin_hbca_c25, project = "Navin_hbca_c25", min.cells = 3, min.features = 200)
Navin_hbca_c25

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c25[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c25, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c25, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c25, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c25, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c25 <- subset(Navin_hbca_c25, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c25 <- NormalizeData(Navin_hbca_c25, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c25 <- FindVariableFeatures(Navin_hbca_c25, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c25)
Navin_hbca_c25 <- ScaleData(Navin_hbca_c25, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c25 <- RunPCA(Navin_hbca_c25, features = VariableFeatures(object = Navin_hbca_c25))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c25[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c25, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c25, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c25 <- JackStraw(Navin_hbca_c25, num.replicate = 100)
Navin_hbca_c25 <- ScoreJackStraw(Navin_hbca_c25, dims = 1:20)
ElbowPlot(Navin_hbca_c25)

# Cluster cells
Navin_hbca_c25 <- FindNeighbors(Navin_hbca_c25, dims = 1:20)
Navin_hbca_c25 <- FindClusters(Navin_hbca_c25, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c25 <- RunUMAP(Navin_hbca_c25, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c25, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c25, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c25.rds")

# Remove Object
rm(Navin_hbca_c25)



##########################################################
########################### Navin_hbca_c26 ###############
##########################################################

# Load the dataset
Navin_hbca_c26 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c26/GSM7500485_hbca_c26_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c26 <- CreateSeuratObject(counts = Navin_hbca_c26, project = "Navin_hbca_c26", min.cells = 3, min.features = 200)
Navin_hbca_c26

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c26[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c26, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c26, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c26, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c26, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c26 <- subset(Navin_hbca_c26, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c26 <- NormalizeData(Navin_hbca_c26, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c26 <- FindVariableFeatures(Navin_hbca_c26, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c26)
Navin_hbca_c26 <- ScaleData(Navin_hbca_c26, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c26 <- RunPCA(Navin_hbca_c26, features = VariableFeatures(object = Navin_hbca_c26))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c26[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c26, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c26, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c26 <- JackStraw(Navin_hbca_c26, num.replicate = 100)
Navin_hbca_c26 <- ScoreJackStraw(Navin_hbca_c26, dims = 1:20)
ElbowPlot(Navin_hbca_c26)

# Cluster cells
Navin_hbca_c26 <- FindNeighbors(Navin_hbca_c26, dims = 1:20)
Navin_hbca_c26 <- FindClusters(Navin_hbca_c26, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c26 <- RunUMAP(Navin_hbca_c26, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c26, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c26, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c26.rds")

# Remove Object
rm(Navin_hbca_c26)



##########################################################
########################### Navin_hbca_c31 ###############
##########################################################

# Load the dataset
Navin_hbca_c31 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c31/GSM7500493_hbca_c31_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c31 <- CreateSeuratObject(counts = Navin_hbca_c31, project = "Navin_hbca_c31", min.cells = 3, min.features = 200)
Navin_hbca_c31

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c31[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c31, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c31, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c31, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c31, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c31 <- subset(Navin_hbca_c31, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c31 <- NormalizeData(Navin_hbca_c31, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c31 <- FindVariableFeatures(Navin_hbca_c31, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c31)
Navin_hbca_c31 <- ScaleData(Navin_hbca_c31, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c31 <- RunPCA(Navin_hbca_c31, features = VariableFeatures(object = Navin_hbca_c31))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c31[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c31, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c31, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c31 <- JackStraw(Navin_hbca_c31, num.replicate = 100)
Navin_hbca_c31 <- ScoreJackStraw(Navin_hbca_c31, dims = 1:20)
ElbowPlot(Navin_hbca_c31)

# Cluster cells
Navin_hbca_c31 <- FindNeighbors(Navin_hbca_c31, dims = 1:20)
Navin_hbca_c31 <- FindClusters(Navin_hbca_c31, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c31 <- RunUMAP(Navin_hbca_c31, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c31, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c31, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c31.rds")

# Remove Object
rm(Navin_hbca_c31)



##########################################################
########################### Navin_hbca_c32 ###############
##########################################################

# Load the dataset
Navin_hbca_c32 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c32/GSM7500494_hbca_c32_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c32 <- CreateSeuratObject(counts = Navin_hbca_c32, project = "Navin_hbca_c32", min.cells = 3, min.features = 200)
Navin_hbca_c32

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c32[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c32, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c32, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c32, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c32, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c32 <- subset(Navin_hbca_c32, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c32 <- NormalizeData(Navin_hbca_c32, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c32 <- FindVariableFeatures(Navin_hbca_c32, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c32)
Navin_hbca_c32 <- ScaleData(Navin_hbca_c32, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c32 <- RunPCA(Navin_hbca_c32, features = VariableFeatures(object = Navin_hbca_c32))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c32[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c32, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c32, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c32 <- JackStraw(Navin_hbca_c32, num.replicate = 100)
Navin_hbca_c32 <- ScoreJackStraw(Navin_hbca_c32, dims = 1:20)
ElbowPlot(Navin_hbca_c32)

# Cluster cells
Navin_hbca_c32 <- FindNeighbors(Navin_hbca_c32, dims = 1:20)
Navin_hbca_c32 <- FindClusters(Navin_hbca_c32, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c32 <- RunUMAP(Navin_hbca_c32, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c32, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c32, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c32.rds")

# Remove Object
rm(Navin_hbca_c32)



##########################################################
########################### Navin_hbca_c50 ###############
##########################################################

# Load the dataset
Navin_hbca_c50 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c50/GSM7500359_hbca_c50_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c50 <- CreateSeuratObject(counts = Navin_hbca_c50, project = "Navin_hbca_c50", min.cells = 3, min.features = 200)
Navin_hbca_c50

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c50[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c50, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c50, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c50, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c50, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c50 <- subset(Navin_hbca_c50, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c50 <- NormalizeData(Navin_hbca_c50, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c50 <- FindVariableFeatures(Navin_hbca_c50, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c50)
Navin_hbca_c50 <- ScaleData(Navin_hbca_c50, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c50 <- RunPCA(Navin_hbca_c50, features = VariableFeatures(object = Navin_hbca_c50))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c50[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c50, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c50, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c50 <- JackStraw(Navin_hbca_c50, num.replicate = 100)
Navin_hbca_c50 <- ScoreJackStraw(Navin_hbca_c50, dims = 1:20)
ElbowPlot(Navin_hbca_c50)

# Cluster cells
Navin_hbca_c50 <- FindNeighbors(Navin_hbca_c50, dims = 1:20)
Navin_hbca_c50 <- FindClusters(Navin_hbca_c50, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c50 <- RunUMAP(Navin_hbca_c50, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c50, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c50, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c50.rds")

# Remove Object
rm(Navin_hbca_c50)



##########################################################
########################### Navin_hbca_c51 ###############
##########################################################

# Load the dataset
Navin_hbca_c51 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c51/GSM7500360_hbca_c51_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c51 <- CreateSeuratObject(counts = Navin_hbca_c51, project = "Navin_hbca_c51", min.cells = 3, min.features = 200)
Navin_hbca_c51

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c51[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c51, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c51, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c51, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c51, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c51 <- subset(Navin_hbca_c51, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c51 <- NormalizeData(Navin_hbca_c51, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c51 <- FindVariableFeatures(Navin_hbca_c51, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c51)
Navin_hbca_c51 <- ScaleData(Navin_hbca_c51, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c51 <- RunPCA(Navin_hbca_c51, features = VariableFeatures(object = Navin_hbca_c51))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c51[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c51, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c51, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c51 <- JackStraw(Navin_hbca_c51, num.replicate = 100)
Navin_hbca_c51 <- ScoreJackStraw(Navin_hbca_c51, dims = 1:20)
ElbowPlot(Navin_hbca_c51)

# Cluster cells
Navin_hbca_c51 <- FindNeighbors(Navin_hbca_c51, dims = 1:20)
Navin_hbca_c51 <- FindClusters(Navin_hbca_c51, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c51 <- RunUMAP(Navin_hbca_c51, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c51, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c51, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c51.rds")

# Remove Object
rm(Navin_hbca_c51)



##########################################################
########################### Navin_hbca_c52 ###############
##########################################################

# Load the dataset
Navin_hbca_c52 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c52/GSM7500361_hbca_c52_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c52 <- CreateSeuratObject(counts = Navin_hbca_c52, project = "Navin_hbca_c52", min.cells = 3, min.features = 200)
Navin_hbca_c52

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c52[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c52, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c52, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c52, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c52, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c52 <- subset(Navin_hbca_c52, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c52 <- NormalizeData(Navin_hbca_c52, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c52 <- FindVariableFeatures(Navin_hbca_c52, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c52)
Navin_hbca_c52 <- ScaleData(Navin_hbca_c52, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c52 <- RunPCA(Navin_hbca_c52, features = VariableFeatures(object = Navin_hbca_c52))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c52[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c52, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c52, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c52 <- JackStraw(Navin_hbca_c52, num.replicate = 100)
Navin_hbca_c52 <- ScoreJackStraw(Navin_hbca_c52, dims = 1:20)
ElbowPlot(Navin_hbca_c52)

# Cluster cells
Navin_hbca_c52 <- FindNeighbors(Navin_hbca_c52, dims = 1:20)
Navin_hbca_c52 <- FindClusters(Navin_hbca_c52, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c52 <- RunUMAP(Navin_hbca_c52, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c52, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c52, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c52.rds")

# Remove Object
rm(Navin_hbca_c52)



##########################################################
########################### Navin_hbca_c53 ###############
##########################################################

# Load the dataset
Navin_hbca_c53 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c53/GSM7500362_hbca_c53_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c53 <- CreateSeuratObject(counts = Navin_hbca_c53, project = "Navin_hbca_c53", min.cells = 3, min.features = 200)
Navin_hbca_c53

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c53[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c53, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c53, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c53, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c53, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c53 <- subset(Navin_hbca_c53, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c53 <- NormalizeData(Navin_hbca_c53, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c53 <- FindVariableFeatures(Navin_hbca_c53, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c53)
Navin_hbca_c53 <- ScaleData(Navin_hbca_c53, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c53 <- RunPCA(Navin_hbca_c53, features = VariableFeatures(object = Navin_hbca_c53))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c53[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c53, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c53, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c53 <- JackStraw(Navin_hbca_c53, num.replicate = 100)
Navin_hbca_c53 <- ScoreJackStraw(Navin_hbca_c53, dims = 1:20)
ElbowPlot(Navin_hbca_c53)

# Cluster cells
Navin_hbca_c53 <- FindNeighbors(Navin_hbca_c53, dims = 1:20)
Navin_hbca_c53 <- FindClusters(Navin_hbca_c53, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c53 <- RunUMAP(Navin_hbca_c53, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c53, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c53, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c53.rds")

# Remove Object
rm(Navin_hbca_c53)



##########################################################
########################### Navin_hbca_c54 ###############
##########################################################

# Load the dataset
Navin_hbca_c54 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c54/GSM7500363_hbca_c54_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c54 <- CreateSeuratObject(counts = Navin_hbca_c54, project = "Navin_hbca_c54", min.cells = 3, min.features = 200)
Navin_hbca_c54

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c54[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c54, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c54, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c54, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c54, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c54 <- subset(Navin_hbca_c54, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c54 <- NormalizeData(Navin_hbca_c54, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c54 <- FindVariableFeatures(Navin_hbca_c54, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c54)
Navin_hbca_c54 <- ScaleData(Navin_hbca_c54, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c54 <- RunPCA(Navin_hbca_c54, features = VariableFeatures(object = Navin_hbca_c54))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c54[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c54, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c54, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c54 <- JackStraw(Navin_hbca_c54, num.replicate = 100)
Navin_hbca_c54 <- ScoreJackStraw(Navin_hbca_c54, dims = 1:20)
ElbowPlot(Navin_hbca_c54)

# Cluster cells
Navin_hbca_c54 <- FindNeighbors(Navin_hbca_c54, dims = 1:20)
Navin_hbca_c54 <- FindClusters(Navin_hbca_c54, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c54 <- RunUMAP(Navin_hbca_c54, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c54, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c54, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c54.rds")

# Remove Object
rm(Navin_hbca_c54)



##########################################################
########################### Navin_hbca_c55 ###############
##########################################################

# Load the dataset
Navin_hbca_c55 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c55/GSM7500364_hbca_c55_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c55 <- CreateSeuratObject(counts = Navin_hbca_c55, project = "Navin_hbca_c55", min.cells = 3, min.features = 200)
Navin_hbca_c55

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c55[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c55, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c55, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c55, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c55, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c55 <- subset(Navin_hbca_c55, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c55 <- NormalizeData(Navin_hbca_c55, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c55 <- FindVariableFeatures(Navin_hbca_c55, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c55)
Navin_hbca_c55 <- ScaleData(Navin_hbca_c55, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c55 <- RunPCA(Navin_hbca_c55, features = VariableFeatures(object = Navin_hbca_c55))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c55[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c55, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c55, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c55 <- JackStraw(Navin_hbca_c55, num.replicate = 100)
Navin_hbca_c55 <- ScoreJackStraw(Navin_hbca_c55, dims = 1:20)
ElbowPlot(Navin_hbca_c55)

# Cluster cells
Navin_hbca_c55 <- FindNeighbors(Navin_hbca_c55, dims = 1:20)
Navin_hbca_c55 <- FindClusters(Navin_hbca_c55, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c55 <- RunUMAP(Navin_hbca_c55, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c55, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c55, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c55.rds")

# Remove Object
rm(Navin_hbca_c55)



##########################################################
########################### Navin_hbca_c56 ###############
##########################################################

# Load the dataset
Navin_hbca_c56 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c56/GSM7500365_hbca_c56_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c56 <- CreateSeuratObject(counts = Navin_hbca_c56, project = "Navin_hbca_c56", min.cells = 3, min.features = 200)
Navin_hbca_c56

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c56[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c56, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c56, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c56, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c56, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c56 <- subset(Navin_hbca_c56, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c56 <- NormalizeData(Navin_hbca_c56, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c56 <- FindVariableFeatures(Navin_hbca_c56, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c56)
Navin_hbca_c56 <- ScaleData(Navin_hbca_c56, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c56 <- RunPCA(Navin_hbca_c56, features = VariableFeatures(object = Navin_hbca_c56))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c56[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c56, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c56, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c56 <- JackStraw(Navin_hbca_c56, num.replicate = 100)
Navin_hbca_c56 <- ScoreJackStraw(Navin_hbca_c56, dims = 1:20)
ElbowPlot(Navin_hbca_c56)

# Cluster cells
Navin_hbca_c56 <- FindNeighbors(Navin_hbca_c56, dims = 1:20)
Navin_hbca_c56 <- FindClusters(Navin_hbca_c56, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c56 <- RunUMAP(Navin_hbca_c56, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c56, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c56, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c56.rds")

# Remove Object
rm(Navin_hbca_c56)



##########################################################
########################### Navin_hbca_c57 ###############
##########################################################

# Load the dataset
Navin_hbca_c57 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c57/GSM7500366_hbca_c57_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c57 <- CreateSeuratObject(counts = Navin_hbca_c57, project = "Navin_hbca_c57", min.cells = 3, min.features = 200)
Navin_hbca_c57

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c57[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c57, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c57, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c57, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c57, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c57 <- subset(Navin_hbca_c57, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c57 <- NormalizeData(Navin_hbca_c57, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c57 <- FindVariableFeatures(Navin_hbca_c57, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c57)
Navin_hbca_c57 <- ScaleData(Navin_hbca_c57, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c57 <- RunPCA(Navin_hbca_c57, features = VariableFeatures(object = Navin_hbca_c57))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c57[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c57, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c57, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c57 <- JackStraw(Navin_hbca_c57, num.replicate = 100)
Navin_hbca_c57 <- ScoreJackStraw(Navin_hbca_c57, dims = 1:20)
ElbowPlot(Navin_hbca_c57)

# Cluster cells
Navin_hbca_c57 <- FindNeighbors(Navin_hbca_c57, dims = 1:20)
Navin_hbca_c57 <- FindClusters(Navin_hbca_c57, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c57 <- RunUMAP(Navin_hbca_c57, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c57, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c57, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c57.rds")

# Remove Object
rm(Navin_hbca_c57)



##########################################################
########################### Navin_hbca_c58 ###############
##########################################################

# Load the dataset
Navin_hbca_c58 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c58/GSM7500367_hbca_c58_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c58 <- CreateSeuratObject(counts = Navin_hbca_c58, project = "Navin_hbca_c58", min.cells = 3, min.features = 200)
Navin_hbca_c58

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c58[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c58, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c58, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c58, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c58, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c58 <- subset(Navin_hbca_c58, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c58 <- NormalizeData(Navin_hbca_c58, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c58 <- FindVariableFeatures(Navin_hbca_c58, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c58)
Navin_hbca_c58 <- ScaleData(Navin_hbca_c58, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c58 <- RunPCA(Navin_hbca_c58, features = VariableFeatures(object = Navin_hbca_c58))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c58[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c58, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c58, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c58 <- JackStraw(Navin_hbca_c58, num.replicate = 100)
Navin_hbca_c58 <- ScoreJackStraw(Navin_hbca_c58, dims = 1:20)
ElbowPlot(Navin_hbca_c58)

# Cluster cells
Navin_hbca_c58 <- FindNeighbors(Navin_hbca_c58, dims = 1:20)
Navin_hbca_c58 <- FindClusters(Navin_hbca_c58, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c58 <- RunUMAP(Navin_hbca_c58, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c58, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c58, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c58.rds")

# Remove Object
rm(Navin_hbca_c58)



##########################################################
########################### Navin_hbca_c59 ###############
##########################################################

# Load the dataset
Navin_hbca_c59 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c59/GSM7500370_hbca_c59_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c59 <- CreateSeuratObject(counts = Navin_hbca_c59, project = "Navin_hbca_c59", min.cells = 3, min.features = 200)
Navin_hbca_c59

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c59[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c59, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c59, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c59, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c59, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c59 <- subset(Navin_hbca_c59, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c59 <- NormalizeData(Navin_hbca_c59, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c59 <- FindVariableFeatures(Navin_hbca_c59, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c59)
Navin_hbca_c59 <- ScaleData(Navin_hbca_c59, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c59 <- RunPCA(Navin_hbca_c59, features = VariableFeatures(object = Navin_hbca_c59))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c59[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c59, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c59, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c59 <- JackStraw(Navin_hbca_c59, num.replicate = 100)
Navin_hbca_c59 <- ScoreJackStraw(Navin_hbca_c59, dims = 1:20)
ElbowPlot(Navin_hbca_c59)

# Cluster cells
Navin_hbca_c59 <- FindNeighbors(Navin_hbca_c59, dims = 1:20)
Navin_hbca_c59 <- FindClusters(Navin_hbca_c59, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c59 <- RunUMAP(Navin_hbca_c59, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c59, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c59, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c59.rds")

# Remove Object
rm(Navin_hbca_c59)



##########################################################
########################### Navin_hbca_c60 ###############
##########################################################

# Load the dataset
Navin_hbca_c60 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c60/GSM7500368_hbca_c60_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c60 <- CreateSeuratObject(counts = Navin_hbca_c60, project = "Navin_hbca_c60", min.cells = 3, min.features = 200)
Navin_hbca_c60

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c60[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c60, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c60, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c60, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c60, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c60 <- subset(Navin_hbca_c60, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c60 <- NormalizeData(Navin_hbca_c60, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c60 <- FindVariableFeatures(Navin_hbca_c60, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c60)
Navin_hbca_c60 <- ScaleData(Navin_hbca_c60, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c60 <- RunPCA(Navin_hbca_c60, features = VariableFeatures(object = Navin_hbca_c60))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c60[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c60, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c60, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c60 <- JackStraw(Navin_hbca_c60, num.replicate = 100)
Navin_hbca_c60 <- ScoreJackStraw(Navin_hbca_c60, dims = 1:20)
ElbowPlot(Navin_hbca_c60)

# Cluster cells
Navin_hbca_c60 <- FindNeighbors(Navin_hbca_c60, dims = 1:20)
Navin_hbca_c60 <- FindClusters(Navin_hbca_c60, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c60 <- RunUMAP(Navin_hbca_c60, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c60, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c60, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c60.rds")

# Remove Object
rm(Navin_hbca_c60)



##########################################################
########################### Navin_hbca_c61 ###############
##########################################################

# Load the dataset
Navin_hbca_c61 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c61/GSM7500369_hbca_c61_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c61 <- CreateSeuratObject(counts = Navin_hbca_c61, project = "Navin_hbca_c61", min.cells = 3, min.features = 200)
Navin_hbca_c61

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c61[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c61, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c61, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c61, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c61, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c61 <- subset(Navin_hbca_c61, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c61 <- NormalizeData(Navin_hbca_c61, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c61 <- FindVariableFeatures(Navin_hbca_c61, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c61)
Navin_hbca_c61 <- ScaleData(Navin_hbca_c61, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c61 <- RunPCA(Navin_hbca_c61, features = VariableFeatures(object = Navin_hbca_c61))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c61[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c61, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c61, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c61 <- JackStraw(Navin_hbca_c61, num.replicate = 100)
Navin_hbca_c61 <- ScoreJackStraw(Navin_hbca_c61, dims = 1:20)
ElbowPlot(Navin_hbca_c61)

# Cluster cells
Navin_hbca_c61 <- FindNeighbors(Navin_hbca_c61, dims = 1:20)
Navin_hbca_c61 <- FindClusters(Navin_hbca_c61, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c61 <- RunUMAP(Navin_hbca_c61, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c61, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c61, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c61.rds")

# Remove Object
rm(Navin_hbca_c61)



##########################################################
########################### Navin_hbca_c62 ###############
##########################################################

# Load the dataset
Navin_hbca_c62 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c62/GSM7500371_hbca_c62_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c62 <- CreateSeuratObject(counts = Navin_hbca_c62, project = "Navin_hbca_c62", min.cells = 3, min.features = 200)
Navin_hbca_c62

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c62[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c62, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c62, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c62, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c62, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c62 <- subset(Navin_hbca_c62, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c62 <- NormalizeData(Navin_hbca_c62, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c62 <- FindVariableFeatures(Navin_hbca_c62, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c62)
Navin_hbca_c62 <- ScaleData(Navin_hbca_c62, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c62 <- RunPCA(Navin_hbca_c62, features = VariableFeatures(object = Navin_hbca_c62))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c62[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c62, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c62, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c62 <- JackStraw(Navin_hbca_c62, num.replicate = 100)
Navin_hbca_c62 <- ScoreJackStraw(Navin_hbca_c62, dims = 1:20)
ElbowPlot(Navin_hbca_c62)

# Cluster cells
Navin_hbca_c62 <- FindNeighbors(Navin_hbca_c62, dims = 1:20)
Navin_hbca_c62 <- FindClusters(Navin_hbca_c62, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c62 <- RunUMAP(Navin_hbca_c62, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c62, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c62, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c62.rds")

# Remove Object
rm(Navin_hbca_c62)



##########################################################
########################### Navin_hbca_c63 ###############
##########################################################

# Load the dataset
Navin_hbca_c63 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c63/GSM7500372_hbca_c63_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c63 <- CreateSeuratObject(counts = Navin_hbca_c63, project = "Navin_hbca_c63", min.cells = 3, min.features = 200)
Navin_hbca_c63

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c63[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c63, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c63, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c63, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c63, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c63 <- subset(Navin_hbca_c63, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c63 <- NormalizeData(Navin_hbca_c63, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c63 <- FindVariableFeatures(Navin_hbca_c63, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c63)
Navin_hbca_c63 <- ScaleData(Navin_hbca_c63, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c63 <- RunPCA(Navin_hbca_c63, features = VariableFeatures(object = Navin_hbca_c63))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c63[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c63, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c63, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c63 <- JackStraw(Navin_hbca_c63, num.replicate = 100)
Navin_hbca_c63 <- ScoreJackStraw(Navin_hbca_c63, dims = 1:20)
ElbowPlot(Navin_hbca_c63)

# Cluster cells
Navin_hbca_c63 <- FindNeighbors(Navin_hbca_c63, dims = 1:20)
Navin_hbca_c63 <- FindClusters(Navin_hbca_c63, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c63 <- RunUMAP(Navin_hbca_c63, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c63, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c63, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c63.rds")

# Remove Object
rm(Navin_hbca_c63)



##########################################################
########################### Navin_hbca_c64 ###############
##########################################################

# Load the dataset
Navin_hbca_c64 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c64/GSM7500373_hbca_c64_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c64 <- CreateSeuratObject(counts = Navin_hbca_c64, project = "Navin_hbca_c64", min.cells = 3, min.features = 200)
Navin_hbca_c64

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c64[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c64, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c64, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c64, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c64, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c64 <- subset(Navin_hbca_c64, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c64 <- NormalizeData(Navin_hbca_c64, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c64 <- FindVariableFeatures(Navin_hbca_c64, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c64)
Navin_hbca_c64 <- ScaleData(Navin_hbca_c64, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c64 <- RunPCA(Navin_hbca_c64, features = VariableFeatures(object = Navin_hbca_c64))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c64[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c64, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c64, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c64 <- JackStraw(Navin_hbca_c64, num.replicate = 100)
Navin_hbca_c64 <- ScoreJackStraw(Navin_hbca_c64, dims = 1:20)
ElbowPlot(Navin_hbca_c64)

# Cluster cells
Navin_hbca_c64 <- FindNeighbors(Navin_hbca_c64, dims = 1:20)
Navin_hbca_c64 <- FindClusters(Navin_hbca_c64, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c64 <- RunUMAP(Navin_hbca_c64, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c64, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c64, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c64.rds")

# Remove Object
rm(Navin_hbca_c64)



##########################################################
########################### Navin_hbca_c65 ###############
##########################################################

# Load the dataset
Navin_hbca_c65 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c65/GSM7500374_hbca_c65_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c65 <- CreateSeuratObject(counts = Navin_hbca_c65, project = "Navin_hbca_c65", min.cells = 3, min.features = 200)
Navin_hbca_c65

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c65[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c65, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c65, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c65, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c65, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c65 <- subset(Navin_hbca_c65, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c65 <- NormalizeData(Navin_hbca_c65, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c65 <- FindVariableFeatures(Navin_hbca_c65, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c65)
Navin_hbca_c65 <- ScaleData(Navin_hbca_c65, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c65 <- RunPCA(Navin_hbca_c65, features = VariableFeatures(object = Navin_hbca_c65))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c65[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c65, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c65, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c65 <- JackStraw(Navin_hbca_c65, num.replicate = 100)
Navin_hbca_c65 <- ScoreJackStraw(Navin_hbca_c65, dims = 1:20)
ElbowPlot(Navin_hbca_c65)

# Cluster cells
Navin_hbca_c65 <- FindNeighbors(Navin_hbca_c65, dims = 1:20)
Navin_hbca_c65 <- FindClusters(Navin_hbca_c65, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c65 <- RunUMAP(Navin_hbca_c65, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c65, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c65, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c65.rds")

# Remove Object
rm(Navin_hbca_c65)



##########################################################
########################### Navin_hbca_c66 ###############
##########################################################

# Load the dataset
Navin_hbca_c66 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c66/GSM7500376_hbca_c66_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c66 <- CreateSeuratObject(counts = Navin_hbca_c66, project = "Navin_hbca_c66", min.cells = 3, min.features = 200)
Navin_hbca_c66

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c66[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c66, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c66, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c66, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c66, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c66 <- subset(Navin_hbca_c66, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c66 <- NormalizeData(Navin_hbca_c66, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c66 <- FindVariableFeatures(Navin_hbca_c66, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c66)
Navin_hbca_c66 <- ScaleData(Navin_hbca_c66, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c66 <- RunPCA(Navin_hbca_c66, features = VariableFeatures(object = Navin_hbca_c66))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c66[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c66, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c66, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c66 <- JackStraw(Navin_hbca_c66, num.replicate = 100)
Navin_hbca_c66 <- ScoreJackStraw(Navin_hbca_c66, dims = 1:20)
ElbowPlot(Navin_hbca_c66)

# Cluster cells
Navin_hbca_c66 <- FindNeighbors(Navin_hbca_c66, dims = 1:20)
Navin_hbca_c66 <- FindClusters(Navin_hbca_c66, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c66 <- RunUMAP(Navin_hbca_c66, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c66, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c66, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c66.rds")

# Remove Object
rm(Navin_hbca_c66)



##########################################################
########################### Navin_hbca_c67 ###############
##########################################################

# Load the dataset
Navin_hbca_c67 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c67/GSM7500375_hbca_c67_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c67 <- CreateSeuratObject(counts = Navin_hbca_c67, project = "Navin_hbca_c67", min.cells = 3, min.features = 200)
Navin_hbca_c67

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c67[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c67, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c67, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c67, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c67, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c67 <- subset(Navin_hbca_c67, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c67 <- NormalizeData(Navin_hbca_c67, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c67 <- FindVariableFeatures(Navin_hbca_c67, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c67)
Navin_hbca_c67 <- ScaleData(Navin_hbca_c67, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c67 <- RunPCA(Navin_hbca_c67, features = VariableFeatures(object = Navin_hbca_c67))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c67[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c67, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c67, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c67 <- JackStraw(Navin_hbca_c67, num.replicate = 100)
Navin_hbca_c67 <- ScoreJackStraw(Navin_hbca_c67, dims = 1:20)
ElbowPlot(Navin_hbca_c67)

# Cluster cells
Navin_hbca_c67 <- FindNeighbors(Navin_hbca_c67, dims = 1:20)
Navin_hbca_c67 <- FindClusters(Navin_hbca_c67, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c67 <- RunUMAP(Navin_hbca_c67, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c67, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c67, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c67.rds")

# Remove Object
rm(Navin_hbca_c67)



##########################################################
########################### Navin_hbca_c68 ###############
##########################################################

# Load the dataset
Navin_hbca_c68 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c68/GSM7500378_hbca_c68_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c68 <- CreateSeuratObject(counts = Navin_hbca_c68, project = "Navin_hbca_c68", min.cells = 3, min.features = 200)
Navin_hbca_c68

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c68[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c68, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c68, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c68, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c68, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c68 <- subset(Navin_hbca_c68, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c68 <- NormalizeData(Navin_hbca_c68, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c68 <- FindVariableFeatures(Navin_hbca_c68, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c68)
Navin_hbca_c68 <- ScaleData(Navin_hbca_c68, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c68 <- RunPCA(Navin_hbca_c68, features = VariableFeatures(object = Navin_hbca_c68))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c68[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c68, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c68, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c68 <- JackStraw(Navin_hbca_c68, num.replicate = 100)
Navin_hbca_c68 <- ScoreJackStraw(Navin_hbca_c68, dims = 1:20)
ElbowPlot(Navin_hbca_c68)

# Cluster cells
Navin_hbca_c68 <- FindNeighbors(Navin_hbca_c68, dims = 1:20)
Navin_hbca_c68 <- FindClusters(Navin_hbca_c68, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c68 <- RunUMAP(Navin_hbca_c68, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c68, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c68, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c68.rds")

# Remove Object
rm(Navin_hbca_c68)



##########################################################
########################### Navin_hbca_c69 ###############
##########################################################

# Load the dataset
Navin_hbca_c69 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c69/GSM7500377_hbca_c69_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c69 <- CreateSeuratObject(counts = Navin_hbca_c69, project = "Navin_hbca_c69", min.cells = 3, min.features = 200)
Navin_hbca_c69

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c69[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c69, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c69, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c69, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c69, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c69 <- subset(Navin_hbca_c69, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c69 <- NormalizeData(Navin_hbca_c69, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c69 <- FindVariableFeatures(Navin_hbca_c69, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c69)
Navin_hbca_c69 <- ScaleData(Navin_hbca_c69, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c69 <- RunPCA(Navin_hbca_c69, features = VariableFeatures(object = Navin_hbca_c69))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c69[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c69, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c69, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c69 <- JackStraw(Navin_hbca_c69, num.replicate = 100)
Navin_hbca_c69 <- ScoreJackStraw(Navin_hbca_c69, dims = 1:20)
ElbowPlot(Navin_hbca_c69)

# Cluster cells
Navin_hbca_c69 <- FindNeighbors(Navin_hbca_c69, dims = 1:20)
Navin_hbca_c69 <- FindClusters(Navin_hbca_c69, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c69 <- RunUMAP(Navin_hbca_c69, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c69, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c69, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c69.rds")

# Remove Object
rm(Navin_hbca_c69)



##########################################################
########################### Navin_hbca_c70 ###############
##########################################################

# Load the dataset
Navin_hbca_c70 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c70/GSM7500379_hbca_c70_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c70 <- CreateSeuratObject(counts = Navin_hbca_c70, project = "Navin_hbca_c70", min.cells = 3, min.features = 200)
Navin_hbca_c70

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c70[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c70, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c70, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c70, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c70, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c70 <- subset(Navin_hbca_c70, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c70 <- NormalizeData(Navin_hbca_c70, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c70 <- FindVariableFeatures(Navin_hbca_c70, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c70)
Navin_hbca_c70 <- ScaleData(Navin_hbca_c70, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c70 <- RunPCA(Navin_hbca_c70, features = VariableFeatures(object = Navin_hbca_c70))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c70[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c70, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c70, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c70 <- JackStraw(Navin_hbca_c70, num.replicate = 100)
Navin_hbca_c70 <- ScoreJackStraw(Navin_hbca_c70, dims = 1:20)
ElbowPlot(Navin_hbca_c70)

# Cluster cells
Navin_hbca_c70 <- FindNeighbors(Navin_hbca_c70, dims = 1:20)
Navin_hbca_c70 <- FindClusters(Navin_hbca_c70, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c70 <- RunUMAP(Navin_hbca_c70, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c70, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c70, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c70.rds")

# Remove Object
rm(Navin_hbca_c70)



##########################################################
########################### Navin_hbca_c71 ###############
##########################################################

# Load the dataset
Navin_hbca_c71 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c71/GSM7500380_hbca_c71_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c71 <- CreateSeuratObject(counts = Navin_hbca_c71, project = "Navin_hbca_c71", min.cells = 3, min.features = 200)
Navin_hbca_c71

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c71[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c71, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c71, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c71, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c71, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c71 <- subset(Navin_hbca_c71, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c71 <- NormalizeData(Navin_hbca_c71, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c71 <- FindVariableFeatures(Navin_hbca_c71, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c71)
Navin_hbca_c71 <- ScaleData(Navin_hbca_c71, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c71 <- RunPCA(Navin_hbca_c71, features = VariableFeatures(object = Navin_hbca_c71))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c71[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c71, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c71, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c71 <- JackStraw(Navin_hbca_c71, num.replicate = 100)
Navin_hbca_c71 <- ScoreJackStraw(Navin_hbca_c71, dims = 1:20)
ElbowPlot(Navin_hbca_c71)

# Cluster cells
Navin_hbca_c71 <- FindNeighbors(Navin_hbca_c71, dims = 1:20)
Navin_hbca_c71 <- FindClusters(Navin_hbca_c71, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c71 <- RunUMAP(Navin_hbca_c71, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c71, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c71, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c71.rds")

# Remove Object
rm(Navin_hbca_c71)



##########################################################
########################### Navin_hbca_c72 ###############
##########################################################

# Load the dataset
Navin_hbca_c72 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c72/GSM7500381_hbca_c72_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c72 <- CreateSeuratObject(counts = Navin_hbca_c72, project = "Navin_hbca_c72", min.cells = 3, min.features = 200)
Navin_hbca_c72

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c72[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c72, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c72, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c72, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c72, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c72 <- subset(Navin_hbca_c72, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c72 <- NormalizeData(Navin_hbca_c72, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c72 <- FindVariableFeatures(Navin_hbca_c72, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c72)
Navin_hbca_c72 <- ScaleData(Navin_hbca_c72, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c72 <- RunPCA(Navin_hbca_c72, features = VariableFeatures(object = Navin_hbca_c72))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c72[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c72, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c72, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c72 <- JackStraw(Navin_hbca_c72, num.replicate = 100)
Navin_hbca_c72 <- ScoreJackStraw(Navin_hbca_c72, dims = 1:20)
ElbowPlot(Navin_hbca_c72)

# Cluster cells
Navin_hbca_c72 <- FindNeighbors(Navin_hbca_c72, dims = 1:20)
Navin_hbca_c72 <- FindClusters(Navin_hbca_c72, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c72 <- RunUMAP(Navin_hbca_c72, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c72, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c72, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c72.rds")

# Remove Object
rm(Navin_hbca_c72)



##########################################################
########################### Navin_hbca_c73 ###############
##########################################################

# Load the dataset
Navin_hbca_c73 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c73/GSM7500382_hbca_c73_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c73 <- CreateSeuratObject(counts = Navin_hbca_c73, project = "Navin_hbca_c73", min.cells = 3, min.features = 200)
Navin_hbca_c73

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c73[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c73, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c73, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c73, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c73, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c73 <- subset(Navin_hbca_c73, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c73 <- NormalizeData(Navin_hbca_c73, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c73 <- FindVariableFeatures(Navin_hbca_c73, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c73)
Navin_hbca_c73 <- ScaleData(Navin_hbca_c73, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c73 <- RunPCA(Navin_hbca_c73, features = VariableFeatures(object = Navin_hbca_c73))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c73[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c73, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c73, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c73 <- JackStraw(Navin_hbca_c73, num.replicate = 100)
Navin_hbca_c73 <- ScoreJackStraw(Navin_hbca_c73, dims = 1:20)
ElbowPlot(Navin_hbca_c73)

# Cluster cells
Navin_hbca_c73 <- FindNeighbors(Navin_hbca_c73, dims = 1:20)
Navin_hbca_c73 <- FindClusters(Navin_hbca_c73, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c73 <- RunUMAP(Navin_hbca_c73, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c73, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c73, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c73.rds")

# Remove Object
rm(Navin_hbca_c73)



##########################################################
########################### Navin_hbca_c74 ###############
##########################################################

# Load the dataset
Navin_hbca_c74 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c74/GSM7500383_hbca_c74_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c74 <- CreateSeuratObject(counts = Navin_hbca_c74, project = "Navin_hbca_c74", min.cells = 3, min.features = 200)
Navin_hbca_c74

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c74[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c74, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c74, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c74, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c74, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c74 <- subset(Navin_hbca_c74, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c74 <- NormalizeData(Navin_hbca_c74, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c74 <- FindVariableFeatures(Navin_hbca_c74, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c74)
Navin_hbca_c74 <- ScaleData(Navin_hbca_c74, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c74 <- RunPCA(Navin_hbca_c74, features = VariableFeatures(object = Navin_hbca_c74))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c74[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c74, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c74, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c74 <- JackStraw(Navin_hbca_c74, num.replicate = 100)
Navin_hbca_c74 <- ScoreJackStraw(Navin_hbca_c74, dims = 1:20)
ElbowPlot(Navin_hbca_c74)

# Cluster cells
Navin_hbca_c74 <- FindNeighbors(Navin_hbca_c74, dims = 1:20)
Navin_hbca_c74 <- FindClusters(Navin_hbca_c74, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c74 <- RunUMAP(Navin_hbca_c74, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c74, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c74, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c74.rds")

# Remove Object
rm(Navin_hbca_c74)



##########################################################
########################### Navin_hbca_c75 ###############
##########################################################

# Load the dataset
Navin_hbca_c75 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c75/GSM7500384_hbca_c75_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c75 <- CreateSeuratObject(counts = Navin_hbca_c75, project = "Navin_hbca_c75", min.cells = 3, min.features = 200)
Navin_hbca_c75

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c75[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c75, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c75, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c75, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c75, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c75 <- subset(Navin_hbca_c75, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c75 <- NormalizeData(Navin_hbca_c75, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c75 <- FindVariableFeatures(Navin_hbca_c75, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c75)
Navin_hbca_c75 <- ScaleData(Navin_hbca_c75, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c75 <- RunPCA(Navin_hbca_c75, features = VariableFeatures(object = Navin_hbca_c75))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c75[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c75, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c75, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c75 <- JackStraw(Navin_hbca_c75, num.replicate = 100)
Navin_hbca_c75 <- ScoreJackStraw(Navin_hbca_c75, dims = 1:20)
ElbowPlot(Navin_hbca_c75)

# Cluster cells
Navin_hbca_c75 <- FindNeighbors(Navin_hbca_c75, dims = 1:20)
Navin_hbca_c75 <- FindClusters(Navin_hbca_c75, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c75 <- RunUMAP(Navin_hbca_c75, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c75, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c75, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c75.rds")

# Remove Object
rm(Navin_hbca_c75)



##########################################################
########################### Navin_hbca_c76 ###############
##########################################################

# Load the dataset
Navin_hbca_c76 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c76/GSM7500385_hbca_c76_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c76 <- CreateSeuratObject(counts = Navin_hbca_c76, project = "Navin_hbca_c76", min.cells = 3, min.features = 200)
Navin_hbca_c76

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c76[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c76, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c76, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c76, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c76, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c76 <- subset(Navin_hbca_c76, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c76 <- NormalizeData(Navin_hbca_c76, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c76 <- FindVariableFeatures(Navin_hbca_c76, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c76)
Navin_hbca_c76 <- ScaleData(Navin_hbca_c76, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c76 <- RunPCA(Navin_hbca_c76, features = VariableFeatures(object = Navin_hbca_c76))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c76[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c76, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c76, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c76 <- JackStraw(Navin_hbca_c76, num.replicate = 100)
Navin_hbca_c76 <- ScoreJackStraw(Navin_hbca_c76, dims = 1:20)
ElbowPlot(Navin_hbca_c76)

# Cluster cells
Navin_hbca_c76 <- FindNeighbors(Navin_hbca_c76, dims = 1:20)
Navin_hbca_c76 <- FindClusters(Navin_hbca_c76, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c76 <- RunUMAP(Navin_hbca_c76, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c76, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c76, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c76.rds")

# Remove Object
rm(Navin_hbca_c76)



##########################################################
########################### Navin_hbca_c77 ###############
##########################################################

# Load the dataset
Navin_hbca_c77 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c77/GSM7500386_hbca_c77_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c77 <- CreateSeuratObject(counts = Navin_hbca_c77, project = "Navin_hbca_c77", min.cells = 3, min.features = 200)
Navin_hbca_c77

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c77[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c77, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c77, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c77, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c77, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c77 <- subset(Navin_hbca_c77, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c77 <- NormalizeData(Navin_hbca_c77, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c77 <- FindVariableFeatures(Navin_hbca_c77, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c77)
Navin_hbca_c77 <- ScaleData(Navin_hbca_c77, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c77 <- RunPCA(Navin_hbca_c77, features = VariableFeatures(object = Navin_hbca_c77))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c77[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c77, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c77, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c77 <- JackStraw(Navin_hbca_c77, num.replicate = 100)
Navin_hbca_c77 <- ScoreJackStraw(Navin_hbca_c77, dims = 1:20)
ElbowPlot(Navin_hbca_c77)

# Cluster cells
Navin_hbca_c77 <- FindNeighbors(Navin_hbca_c77, dims = 1:20)
Navin_hbca_c77 <- FindClusters(Navin_hbca_c77, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c77 <- RunUMAP(Navin_hbca_c77, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c77, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c77, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c77.rds")

# Remove Object
rm(Navin_hbca_c77)



##########################################################
########################### Navin_hbca_c78 ###############
##########################################################

# Load the dataset
Navin_hbca_c78 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c78/GSM7500387_hbca_c78_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c78 <- CreateSeuratObject(counts = Navin_hbca_c78, project = "Navin_hbca_c78", min.cells = 3, min.features = 200)
Navin_hbca_c78

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c78[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c78, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c78, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c78, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c78, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c78 <- subset(Navin_hbca_c78, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c78 <- NormalizeData(Navin_hbca_c78, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c78 <- FindVariableFeatures(Navin_hbca_c78, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c78)
Navin_hbca_c78 <- ScaleData(Navin_hbca_c78, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c78 <- RunPCA(Navin_hbca_c78, features = VariableFeatures(object = Navin_hbca_c78))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c78[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c78, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c78, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c78 <- JackStraw(Navin_hbca_c78, num.replicate = 100)
Navin_hbca_c78 <- ScoreJackStraw(Navin_hbca_c78, dims = 1:20)
ElbowPlot(Navin_hbca_c78)

# Cluster cells
Navin_hbca_c78 <- FindNeighbors(Navin_hbca_c78, dims = 1:20)
Navin_hbca_c78 <- FindClusters(Navin_hbca_c78, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c78 <- RunUMAP(Navin_hbca_c78, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c78, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c78, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c78.rds")

# Remove Object
rm(Navin_hbca_c78)



##########################################################
########################### Navin_hbca_c79 ###############
##########################################################

# Load the dataset
Navin_hbca_c79 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c79/GSM7500388_hbca_c79_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c79 <- CreateSeuratObject(counts = Navin_hbca_c79, project = "Navin_hbca_c79", min.cells = 3, min.features = 200)
Navin_hbca_c79

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c79[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c79, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c79, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c79, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c79, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c79 <- subset(Navin_hbca_c79, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c79 <- NormalizeData(Navin_hbca_c79, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c79 <- FindVariableFeatures(Navin_hbca_c79, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c79)
Navin_hbca_c79 <- ScaleData(Navin_hbca_c79, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c79 <- RunPCA(Navin_hbca_c79, features = VariableFeatures(object = Navin_hbca_c79))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c79[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c79, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c79, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c79 <- JackStraw(Navin_hbca_c79, num.replicate = 100)
Navin_hbca_c79 <- ScoreJackStraw(Navin_hbca_c79, dims = 1:20)
ElbowPlot(Navin_hbca_c79)

# Cluster cells
Navin_hbca_c79 <- FindNeighbors(Navin_hbca_c79, dims = 1:20)
Navin_hbca_c79 <- FindClusters(Navin_hbca_c79, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c79 <- RunUMAP(Navin_hbca_c79, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c79, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c79, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c79.rds")

# Remove Object
rm(Navin_hbca_c79)



##########################################################
########################### Navin_hbca_c80 ###############
##########################################################

# Load the dataset
Navin_hbca_c80 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c80/GSM7500389_hbca_c80_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c80 <- CreateSeuratObject(counts = Navin_hbca_c80, project = "Navin_hbca_c80", min.cells = 3, min.features = 200)
Navin_hbca_c80

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c80[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c80, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c80, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c80, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c80, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c80 <- subset(Navin_hbca_c80, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c80 <- NormalizeData(Navin_hbca_c80, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c80 <- FindVariableFeatures(Navin_hbca_c80, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c80)
Navin_hbca_c80 <- ScaleData(Navin_hbca_c80, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c80 <- RunPCA(Navin_hbca_c80, features = VariableFeatures(object = Navin_hbca_c80))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c80[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c80, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c80, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c80 <- JackStraw(Navin_hbca_c80, num.replicate = 100)
Navin_hbca_c80 <- ScoreJackStraw(Navin_hbca_c80, dims = 1:20)
ElbowPlot(Navin_hbca_c80)

# Cluster cells
Navin_hbca_c80 <- FindNeighbors(Navin_hbca_c80, dims = 1:20)
Navin_hbca_c80 <- FindClusters(Navin_hbca_c80, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c80 <- RunUMAP(Navin_hbca_c80, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c80, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c80, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c80.rds")

# Remove Object
rm(Navin_hbca_c80)



##########################################################
########################### Navin_hbca_c81 ###############
##########################################################

# Load the dataset
Navin_hbca_c81 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c81/GSM7500390_hbca_c81_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c81 <- CreateSeuratObject(counts = Navin_hbca_c81, project = "Navin_hbca_c81", min.cells = 3, min.features = 200)
Navin_hbca_c81

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c81[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c81, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c81, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c81, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c81, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c81 <- subset(Navin_hbca_c81, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c81 <- NormalizeData(Navin_hbca_c81, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c81 <- FindVariableFeatures(Navin_hbca_c81, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c81)
Navin_hbca_c81 <- ScaleData(Navin_hbca_c81, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c81 <- RunPCA(Navin_hbca_c81, features = VariableFeatures(object = Navin_hbca_c81))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c81[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c81, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c81, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c81 <- JackStraw(Navin_hbca_c81, num.replicate = 100)
Navin_hbca_c81 <- ScoreJackStraw(Navin_hbca_c81, dims = 1:20)
ElbowPlot(Navin_hbca_c81)

# Cluster cells
Navin_hbca_c81 <- FindNeighbors(Navin_hbca_c81, dims = 1:20)
Navin_hbca_c81 <- FindClusters(Navin_hbca_c81, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c81 <- RunUMAP(Navin_hbca_c81, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c81, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c81, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c81.rds")

# Remove Object
rm(Navin_hbca_c81)



##########################################################
########################### Navin_hbca_c82 ###############
##########################################################

# Load the dataset
Navin_hbca_c82 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c82/GSM7500391_hbca_c82_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c82 <- CreateSeuratObject(counts = Navin_hbca_c82, project = "Navin_hbca_c82", min.cells = 3, min.features = 200)
Navin_hbca_c82

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c82[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c82, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c82, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c82, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c82, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c82 <- subset(Navin_hbca_c82, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c82 <- NormalizeData(Navin_hbca_c82, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c82 <- FindVariableFeatures(Navin_hbca_c82, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c82)
Navin_hbca_c82 <- ScaleData(Navin_hbca_c82, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c82 <- RunPCA(Navin_hbca_c82, features = VariableFeatures(object = Navin_hbca_c82))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c82[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c82, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c82, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c82 <- JackStraw(Navin_hbca_c82, num.replicate = 100)
Navin_hbca_c82 <- ScoreJackStraw(Navin_hbca_c82, dims = 1:20)
ElbowPlot(Navin_hbca_c82)

# Cluster cells
Navin_hbca_c82 <- FindNeighbors(Navin_hbca_c82, dims = 1:20)
Navin_hbca_c82 <- FindClusters(Navin_hbca_c82, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c82 <- RunUMAP(Navin_hbca_c82, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c82, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c82, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c82.rds")

# Remove Object
rm(Navin_hbca_c82)



##########################################################
########################### Navin_hbca_c83 ###############
##########################################################

# Load the dataset
Navin_hbca_c83 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c83/GSM7500392_hbca_c83_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c83 <- CreateSeuratObject(counts = Navin_hbca_c83, project = "Navin_hbca_c83", min.cells = 3, min.features = 200)
Navin_hbca_c83

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c83[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c83, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c83, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c83, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c83, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c83 <- subset(Navin_hbca_c83, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c83 <- NormalizeData(Navin_hbca_c83, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c83 <- FindVariableFeatures(Navin_hbca_c83, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c83)
Navin_hbca_c83 <- ScaleData(Navin_hbca_c83, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c83 <- RunPCA(Navin_hbca_c83, features = VariableFeatures(object = Navin_hbca_c83))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c83[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c83, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c83, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c83 <- JackStraw(Navin_hbca_c83, num.replicate = 100)
Navin_hbca_c83 <- ScoreJackStraw(Navin_hbca_c83, dims = 1:20)
ElbowPlot(Navin_hbca_c83)

# Cluster cells
Navin_hbca_c83 <- FindNeighbors(Navin_hbca_c83, dims = 1:20)
Navin_hbca_c83 <- FindClusters(Navin_hbca_c83, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c83 <- RunUMAP(Navin_hbca_c83, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c83, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c83, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c83.rds")

# Remove Object
rm(Navin_hbca_c83)



##########################################################
########################### Navin_hbca_c84 ###############
##########################################################

# Load the dataset
Navin_hbca_c84 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c84/GSM7500393_hbca_c84_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c84 <- CreateSeuratObject(counts = Navin_hbca_c84, project = "Navin_hbca_c84", min.cells = 3, min.features = 200)
Navin_hbca_c84

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c84[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c84, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c84, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c84, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c84, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c84 <- subset(Navin_hbca_c84, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c84 <- NormalizeData(Navin_hbca_c84, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c84 <- FindVariableFeatures(Navin_hbca_c84, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c84)
Navin_hbca_c84 <- ScaleData(Navin_hbca_c84, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c84 <- RunPCA(Navin_hbca_c84, features = VariableFeatures(object = Navin_hbca_c84))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c84[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c84, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c84, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c84 <- JackStraw(Navin_hbca_c84, num.replicate = 100)
Navin_hbca_c84 <- ScoreJackStraw(Navin_hbca_c84, dims = 1:20)
ElbowPlot(Navin_hbca_c84)

# Cluster cells
Navin_hbca_c84 <- FindNeighbors(Navin_hbca_c84, dims = 1:20)
Navin_hbca_c84 <- FindClusters(Navin_hbca_c84, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c84 <- RunUMAP(Navin_hbca_c84, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c84, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c84, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c84.rds")

# Remove Object
rm(Navin_hbca_c84)



##########################################################
########################### Navin_hbca_c85 ###############
##########################################################

# Load the dataset
Navin_hbca_c85 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c85/GSM7500394_hbca_c85_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c85 <- CreateSeuratObject(counts = Navin_hbca_c85, project = "Navin_hbca_c85", min.cells = 3, min.features = 200)
Navin_hbca_c85

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c85[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c85, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c85, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c85, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c85, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c85 <- subset(Navin_hbca_c85, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c85 <- NormalizeData(Navin_hbca_c85, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c85 <- FindVariableFeatures(Navin_hbca_c85, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c85)
Navin_hbca_c85 <- ScaleData(Navin_hbca_c85, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c85 <- RunPCA(Navin_hbca_c85, features = VariableFeatures(object = Navin_hbca_c85))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c85[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c85, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c85, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c85 <- JackStraw(Navin_hbca_c85, num.replicate = 100)
Navin_hbca_c85 <- ScoreJackStraw(Navin_hbca_c85, dims = 1:20)
ElbowPlot(Navin_hbca_c85)

# Cluster cells
Navin_hbca_c85 <- FindNeighbors(Navin_hbca_c85, dims = 1:20)
Navin_hbca_c85 <- FindClusters(Navin_hbca_c85, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c85 <- RunUMAP(Navin_hbca_c85, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c85, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c85, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c85.rds")

# Remove Object
rm(Navin_hbca_c85)



##########################################################
########################### Navin_hbca_c86 ###############
##########################################################

# Load the dataset
Navin_hbca_c86 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c86/GSM7500395_hbca_c86_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c86 <- CreateSeuratObject(counts = Navin_hbca_c86, project = "Navin_hbca_c86", min.cells = 3, min.features = 200)
Navin_hbca_c86

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c86[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c86, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c86, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c86, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c86, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c86 <- subset(Navin_hbca_c86, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c86 <- NormalizeData(Navin_hbca_c86, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c86 <- FindVariableFeatures(Navin_hbca_c86, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c86)
Navin_hbca_c86 <- ScaleData(Navin_hbca_c86, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c86 <- RunPCA(Navin_hbca_c86, features = VariableFeatures(object = Navin_hbca_c86))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c86[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c86, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c86, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c86 <- JackStraw(Navin_hbca_c86, num.replicate = 100)
Navin_hbca_c86 <- ScoreJackStraw(Navin_hbca_c86, dims = 1:20)
ElbowPlot(Navin_hbca_c86)

# Cluster cells
Navin_hbca_c86 <- FindNeighbors(Navin_hbca_c86, dims = 1:20)
Navin_hbca_c86 <- FindClusters(Navin_hbca_c86, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c86 <- RunUMAP(Navin_hbca_c86, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c86, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c86, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c86.rds")

# Remove Object
rm(Navin_hbca_c86)



##########################################################
########################### Navin_hbca_c87 ###############
##########################################################

# Load the dataset
Navin_hbca_c87 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c87/GSM7500396_hbca_c87_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c87 <- CreateSeuratObject(counts = Navin_hbca_c87, project = "Navin_hbca_c87", min.cells = 3, min.features = 200)
Navin_hbca_c87

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c87[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c87, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c87, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c87, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c87, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c87 <- subset(Navin_hbca_c87, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c87 <- NormalizeData(Navin_hbca_c87, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c87 <- FindVariableFeatures(Navin_hbca_c87, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c87)
Navin_hbca_c87 <- ScaleData(Navin_hbca_c87, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c87 <- RunPCA(Navin_hbca_c87, features = VariableFeatures(object = Navin_hbca_c87))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c87[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c87, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c87, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c87 <- JackStraw(Navin_hbca_c87, num.replicate = 100)
Navin_hbca_c87 <- ScoreJackStraw(Navin_hbca_c87, dims = 1:20)
ElbowPlot(Navin_hbca_c87)

# Cluster cells
Navin_hbca_c87 <- FindNeighbors(Navin_hbca_c87, dims = 1:20)
Navin_hbca_c87 <- FindClusters(Navin_hbca_c87, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c87 <- RunUMAP(Navin_hbca_c87, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c87, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c87, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c87.rds")

# Remove Object
rm(Navin_hbca_c87)



##########################################################
########################### Navin_hbca_c88 ###############
##########################################################

# Load the dataset
Navin_hbca_c88 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c88/GSM7500397_hbca_c88_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c88 <- CreateSeuratObject(counts = Navin_hbca_c88, project = "Navin_hbca_c88", min.cells = 3, min.features = 200)
Navin_hbca_c88

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c88[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c88, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c88, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c88, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c88, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c88 <- subset(Navin_hbca_c88, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c88 <- NormalizeData(Navin_hbca_c88, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c88 <- FindVariableFeatures(Navin_hbca_c88, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c88)
Navin_hbca_c88 <- ScaleData(Navin_hbca_c88, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c88 <- RunPCA(Navin_hbca_c88, features = VariableFeatures(object = Navin_hbca_c88))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c88[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c88, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c88, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c88 <- JackStraw(Navin_hbca_c88, num.replicate = 100)
Navin_hbca_c88 <- ScoreJackStraw(Navin_hbca_c88, dims = 1:20)
ElbowPlot(Navin_hbca_c88)

# Cluster cells
Navin_hbca_c88 <- FindNeighbors(Navin_hbca_c88, dims = 1:20)
Navin_hbca_c88 <- FindClusters(Navin_hbca_c88, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c88 <- RunUMAP(Navin_hbca_c88, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c88, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c88, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c88.rds")

# Remove Object
rm(Navin_hbca_c88)



##########################################################
########################### Navin_hbca_c89 ###############
##########################################################

# Load the dataset
Navin_hbca_c89 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c89/GSM7500398_hbca_c89_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c89 <- CreateSeuratObject(counts = Navin_hbca_c89, project = "Navin_hbca_c89", min.cells = 3, min.features = 200)
Navin_hbca_c89

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c89[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c89, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c89, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c89, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c89, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c89 <- subset(Navin_hbca_c89, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c89 <- NormalizeData(Navin_hbca_c89, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c89 <- FindVariableFeatures(Navin_hbca_c89, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c89)
Navin_hbca_c89 <- ScaleData(Navin_hbca_c89, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c89 <- RunPCA(Navin_hbca_c89, features = VariableFeatures(object = Navin_hbca_c89))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c89[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c89, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c89, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c89 <- JackStraw(Navin_hbca_c89, num.replicate = 100)
Navin_hbca_c89 <- ScoreJackStraw(Navin_hbca_c89, dims = 1:20)
ElbowPlot(Navin_hbca_c89)

# Cluster cells
Navin_hbca_c89 <- FindNeighbors(Navin_hbca_c89, dims = 1:20)
Navin_hbca_c89 <- FindClusters(Navin_hbca_c89, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c89 <- RunUMAP(Navin_hbca_c89, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c89, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c89, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c89.rds")

# Remove Object
rm(Navin_hbca_c89)



##########################################################
########################### Navin_hbca_c90 ###############
##########################################################

# Load the dataset
Navin_hbca_c90 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c90/GSM7500399_hbca_c90_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c90 <- CreateSeuratObject(counts = Navin_hbca_c90, project = "Navin_hbca_c90", min.cells = 3, min.features = 200)
Navin_hbca_c90

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c90[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c90, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c90, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c90, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c90, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c90 <- subset(Navin_hbca_c90, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c90 <- NormalizeData(Navin_hbca_c90, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c90 <- FindVariableFeatures(Navin_hbca_c90, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c90)
Navin_hbca_c90 <- ScaleData(Navin_hbca_c90, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c90 <- RunPCA(Navin_hbca_c90, features = VariableFeatures(object = Navin_hbca_c90))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c90[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c90, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c90, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c90 <- JackStraw(Navin_hbca_c90, num.replicate = 100)
Navin_hbca_c90 <- ScoreJackStraw(Navin_hbca_c90, dims = 1:20)
ElbowPlot(Navin_hbca_c90)

# Cluster cells
Navin_hbca_c90 <- FindNeighbors(Navin_hbca_c90, dims = 1:20)
Navin_hbca_c90 <- FindClusters(Navin_hbca_c90, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c90 <- RunUMAP(Navin_hbca_c90, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c90, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c90, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c90.rds")

# Remove Object
rm(Navin_hbca_c90)



##########################################################
########################### Navin_hbca_c91 ###############
##########################################################

# Load the dataset
Navin_hbca_c91 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c91/GSM7500400_hbca_c91_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c91 <- CreateSeuratObject(counts = Navin_hbca_c91, project = "Navin_hbca_c91", min.cells = 3, min.features = 200)
Navin_hbca_c91

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c91[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c91, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c91, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c91, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c91, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c91 <- subset(Navin_hbca_c91, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c91 <- NormalizeData(Navin_hbca_c91, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c91 <- FindVariableFeatures(Navin_hbca_c91, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c91)
Navin_hbca_c91 <- ScaleData(Navin_hbca_c91, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c91 <- RunPCA(Navin_hbca_c91, features = VariableFeatures(object = Navin_hbca_c91))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c91[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c91, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c91, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c91 <- JackStraw(Navin_hbca_c91, num.replicate = 100)
Navin_hbca_c91 <- ScoreJackStraw(Navin_hbca_c91, dims = 1:20)
ElbowPlot(Navin_hbca_c91)

# Cluster cells
Navin_hbca_c91 <- FindNeighbors(Navin_hbca_c91, dims = 1:20)
Navin_hbca_c91 <- FindClusters(Navin_hbca_c91, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c91 <- RunUMAP(Navin_hbca_c91, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c91, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c91, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c91.rds")

# Remove Object
rm(Navin_hbca_c91)



##########################################################
########################### Navin_hbca_c92 ###############
##########################################################

# Load the dataset
Navin_hbca_c92 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c92/GSM7500401_hbca_c92_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c92 <- CreateSeuratObject(counts = Navin_hbca_c92, project = "Navin_hbca_c92", min.cells = 3, min.features = 200)
Navin_hbca_c92

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c92[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c92, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c92, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c92, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c92, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c92 <- subset(Navin_hbca_c92, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c92 <- NormalizeData(Navin_hbca_c92, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c92 <- FindVariableFeatures(Navin_hbca_c92, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c92)
Navin_hbca_c92 <- ScaleData(Navin_hbca_c92, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c92 <- RunPCA(Navin_hbca_c92, features = VariableFeatures(object = Navin_hbca_c92))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c92[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c92, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c92, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c92 <- JackStraw(Navin_hbca_c92, num.replicate = 100)
Navin_hbca_c92 <- ScoreJackStraw(Navin_hbca_c92, dims = 1:20)
ElbowPlot(Navin_hbca_c92)

# Cluster cells
Navin_hbca_c92 <- FindNeighbors(Navin_hbca_c92, dims = 1:20)
Navin_hbca_c92 <- FindClusters(Navin_hbca_c92, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c92 <- RunUMAP(Navin_hbca_c92, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c92, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c92, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c92.rds")

# Remove Object
rm(Navin_hbca_c92)



##########################################################
########################### Navin_hbca_c93 ###############
##########################################################

# Load the dataset
Navin_hbca_c93 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c93/GSM7500403_hbca_c93_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c93 <- CreateSeuratObject(counts = Navin_hbca_c93, project = "Navin_hbca_c93", min.cells = 3, min.features = 200)
Navin_hbca_c93

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c93[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c93, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c93, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c93, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c93, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c93 <- subset(Navin_hbca_c93, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c93 <- NormalizeData(Navin_hbca_c93, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c93 <- FindVariableFeatures(Navin_hbca_c93, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c93)
Navin_hbca_c93 <- ScaleData(Navin_hbca_c93, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c93 <- RunPCA(Navin_hbca_c93, features = VariableFeatures(object = Navin_hbca_c93))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c93[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c93, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c93, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c93 <- JackStraw(Navin_hbca_c93, num.replicate = 100)
Navin_hbca_c93 <- ScoreJackStraw(Navin_hbca_c93, dims = 1:20)
ElbowPlot(Navin_hbca_c93)

# Cluster cells
Navin_hbca_c93 <- FindNeighbors(Navin_hbca_c93, dims = 1:20)
Navin_hbca_c93 <- FindClusters(Navin_hbca_c93, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c93 <- RunUMAP(Navin_hbca_c93, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c93, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c93, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c93.rds")

# Remove Object
rm(Navin_hbca_c93)



##########################################################
########################### Navin_hbca_c94 ###############
##########################################################

# Load the dataset
Navin_hbca_c94 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c94/GSM7500404_hbca_c94_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c94 <- CreateSeuratObject(counts = Navin_hbca_c94, project = "Navin_hbca_c94", min.cells = 3, min.features = 200)
Navin_hbca_c94

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c94[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c94, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c94, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c94, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c94, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c94 <- subset(Navin_hbca_c94, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c94 <- NormalizeData(Navin_hbca_c94, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c94 <- FindVariableFeatures(Navin_hbca_c94, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c94)
Navin_hbca_c94 <- ScaleData(Navin_hbca_c94, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c94 <- RunPCA(Navin_hbca_c94, features = VariableFeatures(object = Navin_hbca_c94))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c94[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c94, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c94, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c94 <- JackStraw(Navin_hbca_c94, num.replicate = 100)
Navin_hbca_c94 <- ScoreJackStraw(Navin_hbca_c94, dims = 1:20)
ElbowPlot(Navin_hbca_c94)

# Cluster cells
Navin_hbca_c94 <- FindNeighbors(Navin_hbca_c94, dims = 1:20)
Navin_hbca_c94 <- FindClusters(Navin_hbca_c94, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c94 <- RunUMAP(Navin_hbca_c94, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c94, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c94, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c94.rds")

# Remove Object
rm(Navin_hbca_c94)



##########################################################
########################### Navin_hbca_c95 ###############
##########################################################

# Load the dataset
Navin_hbca_c95 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c95/GSM7500405_hbca_c95_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c95 <- CreateSeuratObject(counts = Navin_hbca_c95, project = "Navin_hbca_c95", min.cells = 3, min.features = 200)
Navin_hbca_c95

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c95[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c95, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c95, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c95, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c95, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c95 <- subset(Navin_hbca_c95, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c95 <- NormalizeData(Navin_hbca_c95, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c95 <- FindVariableFeatures(Navin_hbca_c95, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c95)
Navin_hbca_c95 <- ScaleData(Navin_hbca_c95, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c95 <- RunPCA(Navin_hbca_c95, features = VariableFeatures(object = Navin_hbca_c95))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c95[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c95, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c95, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c95 <- JackStraw(Navin_hbca_c95, num.replicate = 100)
Navin_hbca_c95 <- ScoreJackStraw(Navin_hbca_c95, dims = 1:20)
ElbowPlot(Navin_hbca_c95)

# Cluster cells
Navin_hbca_c95 <- FindNeighbors(Navin_hbca_c95, dims = 1:20)
Navin_hbca_c95 <- FindClusters(Navin_hbca_c95, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c95 <- RunUMAP(Navin_hbca_c95, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c95, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c95, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c95.rds")

# Remove Object
rm(Navin_hbca_c95)



##########################################################
########################### Navin_hbca_c96 ###############
##########################################################

# Load the dataset
Navin_hbca_c96 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c96/GSM7500406_hbca_c96_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c96 <- CreateSeuratObject(counts = Navin_hbca_c96, project = "Navin_hbca_c96", min.cells = 3, min.features = 200)
Navin_hbca_c96

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c96[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c96, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c96, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c96, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c96, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c96 <- subset(Navin_hbca_c96, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c96 <- NormalizeData(Navin_hbca_c96, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c96 <- FindVariableFeatures(Navin_hbca_c96, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c96)
Navin_hbca_c96 <- ScaleData(Navin_hbca_c96, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c96 <- RunPCA(Navin_hbca_c96, features = VariableFeatures(object = Navin_hbca_c96))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c96[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c96, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c96, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c96 <- JackStraw(Navin_hbca_c96, num.replicate = 100)
Navin_hbca_c96 <- ScoreJackStraw(Navin_hbca_c96, dims = 1:20)
ElbowPlot(Navin_hbca_c96)

# Cluster cells
Navin_hbca_c96 <- FindNeighbors(Navin_hbca_c96, dims = 1:20)
Navin_hbca_c96 <- FindClusters(Navin_hbca_c96, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c96 <- RunUMAP(Navin_hbca_c96, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c96, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c96, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c96.rds")

# Remove Object
rm(Navin_hbca_c96)



##########################################################
########################### Navin_hbca_c97 ###############
##########################################################

# Load the dataset
Navin_hbca_c97 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c97/GSM7500407_hbca_c97_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c97 <- CreateSeuratObject(counts = Navin_hbca_c97, project = "Navin_hbca_c97", min.cells = 3, min.features = 200)
Navin_hbca_c97

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c97[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c97, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c97, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c97, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c97, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c97 <- subset(Navin_hbca_c97, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c97 <- NormalizeData(Navin_hbca_c97, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c97 <- FindVariableFeatures(Navin_hbca_c97, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c97)
Navin_hbca_c97 <- ScaleData(Navin_hbca_c97, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c97 <- RunPCA(Navin_hbca_c97, features = VariableFeatures(object = Navin_hbca_c97))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c97[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c97, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c97, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c97 <- JackStraw(Navin_hbca_c97, num.replicate = 100)
Navin_hbca_c97 <- ScoreJackStraw(Navin_hbca_c97, dims = 1:20)
ElbowPlot(Navin_hbca_c97)

# Cluster cells
Navin_hbca_c97 <- FindNeighbors(Navin_hbca_c97, dims = 1:20)
Navin_hbca_c97 <- FindClusters(Navin_hbca_c97, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c97 <- RunUMAP(Navin_hbca_c97, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c97, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c97, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c97.rds")

# Remove Object
rm(Navin_hbca_c97)



##########################################################
########################### Navin_hbca_c98 ###############
##########################################################

# Load the dataset
Navin_hbca_c98 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c98/GSM7500408_hbca_c98_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c98 <- CreateSeuratObject(counts = Navin_hbca_c98, project = "Navin_hbca_c98", min.cells = 3, min.features = 200)
Navin_hbca_c98

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c98[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c98, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c98, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c98, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c98, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c98 <- subset(Navin_hbca_c98, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c98 <- NormalizeData(Navin_hbca_c98, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c98 <- FindVariableFeatures(Navin_hbca_c98, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c98)
Navin_hbca_c98 <- ScaleData(Navin_hbca_c98, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c98 <- RunPCA(Navin_hbca_c98, features = VariableFeatures(object = Navin_hbca_c98))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c98[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c98, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c98, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c98 <- JackStraw(Navin_hbca_c98, num.replicate = 100)
Navin_hbca_c98 <- ScoreJackStraw(Navin_hbca_c98, dims = 1:20)
ElbowPlot(Navin_hbca_c98)

# Cluster cells
Navin_hbca_c98 <- FindNeighbors(Navin_hbca_c98, dims = 1:20)
Navin_hbca_c98 <- FindClusters(Navin_hbca_c98, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c98 <- RunUMAP(Navin_hbca_c98, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c98, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c98, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c98.rds")

# Remove Object
rm(Navin_hbca_c98)



##########################################################
########################### Navin_hbca_c99 ###############
##########################################################

# Load the dataset
Navin_hbca_c99 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c99/GSM7500409_hbca_c99_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c99 <- CreateSeuratObject(counts = Navin_hbca_c99, project = "Navin_hbca_c99", min.cells = 3, min.features = 200)
Navin_hbca_c99

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c99[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c99, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c99, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c99, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c99, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c99 <- subset(Navin_hbca_c99, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c99 <- NormalizeData(Navin_hbca_c99, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c99 <- FindVariableFeatures(Navin_hbca_c99, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c99)
Navin_hbca_c99 <- ScaleData(Navin_hbca_c99, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c99 <- RunPCA(Navin_hbca_c99, features = VariableFeatures(object = Navin_hbca_c99))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c99[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c99, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c99, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c99 <- JackStraw(Navin_hbca_c99, num.replicate = 100)
Navin_hbca_c99 <- ScoreJackStraw(Navin_hbca_c99, dims = 1:20)
ElbowPlot(Navin_hbca_c99)

# Cluster cells
Navin_hbca_c99 <- FindNeighbors(Navin_hbca_c99, dims = 1:20)
Navin_hbca_c99 <- FindClusters(Navin_hbca_c99, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c99 <- RunUMAP(Navin_hbca_c99, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c99, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c99, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c99.rds")

# Remove Object
rm(Navin_hbca_c99)



##########################################################
########################### Navin_hbca_c100 ##############
##########################################################

# Load the dataset
Navin_hbca_c100 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c100/GSM7500402_hbca_c100_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c100 <- CreateSeuratObject(counts = Navin_hbca_c100, project = "Navin_hbca_c100", min.cells = 3, min.features = 200)
Navin_hbca_c100

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c100[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c100, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c100, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c100, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c100, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c100 <- subset(Navin_hbca_c100, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c100 <- NormalizeData(Navin_hbca_c100, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c100 <- FindVariableFeatures(Navin_hbca_c100, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c100)
Navin_hbca_c100 <- ScaleData(Navin_hbca_c100, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c100 <- RunPCA(Navin_hbca_c100, features = VariableFeatures(object = Navin_hbca_c100))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c100[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c100, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c100, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c100 <- JackStraw(Navin_hbca_c100, num.replicate = 100)
Navin_hbca_c100 <- ScoreJackStraw(Navin_hbca_c100, dims = 1:20)
ElbowPlot(Navin_hbca_c100)

# Cluster cells
Navin_hbca_c100 <- FindNeighbors(Navin_hbca_c100, dims = 1:20)
Navin_hbca_c100 <- FindClusters(Navin_hbca_c100, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c100 <- RunUMAP(Navin_hbca_c100, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c100, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c100, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c100.rds")

# Remove Object
rm(Navin_hbca_c100)



##########################################################
########################### Navin_hbca_c101 ##############
##########################################################

# Load the dataset
Navin_hbca_c101 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c101/GSM7500410_hbca_c101_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c101 <- CreateSeuratObject(counts = Navin_hbca_c101, project = "Navin_hbca_c101", min.cells = 3, min.features = 200)
Navin_hbca_c101

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c101[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c101, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c101, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c101, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c101, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c101 <- subset(Navin_hbca_c101, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c101 <- NormalizeData(Navin_hbca_c101, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c101 <- FindVariableFeatures(Navin_hbca_c101, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c101)
Navin_hbca_c101 <- ScaleData(Navin_hbca_c101, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c101 <- RunPCA(Navin_hbca_c101, features = VariableFeatures(object = Navin_hbca_c101))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c101[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c101, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c101, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c101 <- JackStraw(Navin_hbca_c101, num.replicate = 100)
Navin_hbca_c101 <- ScoreJackStraw(Navin_hbca_c101, dims = 1:20)
ElbowPlot(Navin_hbca_c101)

# Cluster cells
Navin_hbca_c101 <- FindNeighbors(Navin_hbca_c101, dims = 1:20)
Navin_hbca_c101 <- FindClusters(Navin_hbca_c101, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c101 <- RunUMAP(Navin_hbca_c101, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c101, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c101, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c101.rds")

# Remove Object
rm(Navin_hbca_c101)



##########################################################
########################### Navin_hbca_c102 ##############
##########################################################

# Load the dataset
Navin_hbca_c102 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c102/GSM7500411_hbca_c102_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c102 <- CreateSeuratObject(counts = Navin_hbca_c102, project = "Navin_hbca_c102", min.cells = 3, min.features = 200)
Navin_hbca_c102

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c102[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c102, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c102, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c102, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c102, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c102 <- subset(Navin_hbca_c102, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c102 <- NormalizeData(Navin_hbca_c102, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c102 <- FindVariableFeatures(Navin_hbca_c102, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c102)
Navin_hbca_c102 <- ScaleData(Navin_hbca_c102, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c102 <- RunPCA(Navin_hbca_c102, features = VariableFeatures(object = Navin_hbca_c102))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c102[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c102, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c102, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c102 <- JackStraw(Navin_hbca_c102, num.replicate = 100)
Navin_hbca_c102 <- ScoreJackStraw(Navin_hbca_c102, dims = 1:20)
ElbowPlot(Navin_hbca_c102)

# Cluster cells
Navin_hbca_c102 <- FindNeighbors(Navin_hbca_c102, dims = 1:20)
Navin_hbca_c102 <- FindClusters(Navin_hbca_c102, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c102 <- RunUMAP(Navin_hbca_c102, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c102, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c102, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c102.rds")

# Remove Object
rm(Navin_hbca_c102)



##########################################################
########################### Navin_hbca_c103 ##############
##########################################################

# Load the dataset
Navin_hbca_c103 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c103/GSM7500412_hbca_c103_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c103 <- CreateSeuratObject(counts = Navin_hbca_c103, project = "Navin_hbca_c103", min.cells = 3, min.features = 200)
Navin_hbca_c103

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c103[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c103, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c103, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c103, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c103, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c103 <- subset(Navin_hbca_c103, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c103 <- NormalizeData(Navin_hbca_c103, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c103 <- FindVariableFeatures(Navin_hbca_c103, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c103)
Navin_hbca_c103 <- ScaleData(Navin_hbca_c103, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c103 <- RunPCA(Navin_hbca_c103, features = VariableFeatures(object = Navin_hbca_c103))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c103[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c103, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c103, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c103 <- JackStraw(Navin_hbca_c103, num.replicate = 100)
Navin_hbca_c103 <- ScoreJackStraw(Navin_hbca_c103, dims = 1:20)
ElbowPlot(Navin_hbca_c103)

# Cluster cells
Navin_hbca_c103 <- FindNeighbors(Navin_hbca_c103, dims = 1:20)
Navin_hbca_c103 <- FindClusters(Navin_hbca_c103, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c103 <- RunUMAP(Navin_hbca_c103, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c103, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c103, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c103.rds")

# Remove Object
rm(Navin_hbca_c103)



##########################################################
########################### Navin_hbca_c104 ##############
##########################################################

# Load the dataset
Navin_hbca_c104 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c104/GSM7500413_hbca_c104_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c104 <- CreateSeuratObject(counts = Navin_hbca_c104, project = "Navin_hbca_c104", min.cells = 3, min.features = 200)
Navin_hbca_c104

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c104[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c104, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c104, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c104, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c104, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c104 <- subset(Navin_hbca_c104, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c104 <- NormalizeData(Navin_hbca_c104, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c104 <- FindVariableFeatures(Navin_hbca_c104, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c104)
Navin_hbca_c104 <- ScaleData(Navin_hbca_c104, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c104 <- RunPCA(Navin_hbca_c104, features = VariableFeatures(object = Navin_hbca_c104))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c104[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c104, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c104, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c104 <- JackStraw(Navin_hbca_c104, num.replicate = 100)
Navin_hbca_c104 <- ScoreJackStraw(Navin_hbca_c104, dims = 1:20)
ElbowPlot(Navin_hbca_c104)

# Cluster cells
Navin_hbca_c104 <- FindNeighbors(Navin_hbca_c104, dims = 1:20)
Navin_hbca_c104 <- FindClusters(Navin_hbca_c104, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c104 <- RunUMAP(Navin_hbca_c104, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c104, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c104, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c104.rds")

# Remove Object
rm(Navin_hbca_c104)



##########################################################
########################### Navin_hbca_c105 ##############
##########################################################

# Load the dataset
Navin_hbca_c105 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c105/GSM7500414_hbca_c105_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c105 <- CreateSeuratObject(counts = Navin_hbca_c105, project = "Navin_hbca_c105", min.cells = 3, min.features = 200)
Navin_hbca_c105

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c105[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c105, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c105, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c105, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c105, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c105 <- subset(Navin_hbca_c105, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c105 <- NormalizeData(Navin_hbca_c105, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c105 <- FindVariableFeatures(Navin_hbca_c105, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c105)
Navin_hbca_c105 <- ScaleData(Navin_hbca_c105, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c105 <- RunPCA(Navin_hbca_c105, features = VariableFeatures(object = Navin_hbca_c105))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c105[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c105, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c105, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c105 <- JackStraw(Navin_hbca_c105, num.replicate = 100)
Navin_hbca_c105 <- ScoreJackStraw(Navin_hbca_c105, dims = 1:20)
ElbowPlot(Navin_hbca_c105)

# Cluster cells
Navin_hbca_c105 <- FindNeighbors(Navin_hbca_c105, dims = 1:20)
Navin_hbca_c105 <- FindClusters(Navin_hbca_c105, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c105 <- RunUMAP(Navin_hbca_c105, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c105, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c105, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c105.rds")

# Remove Object
rm(Navin_hbca_c105)



##########################################################
########################### Navin_hbca_c106 ##############
##########################################################

# Load the dataset
Navin_hbca_c106 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c106/GSM7500415_hbca_c106_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c106 <- CreateSeuratObject(counts = Navin_hbca_c106, project = "Navin_hbca_c106", min.cells = 3, min.features = 200)
Navin_hbca_c106

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c106[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c106, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c106, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c106, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c106, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c106 <- subset(Navin_hbca_c106, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c106 <- NormalizeData(Navin_hbca_c106, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c106 <- FindVariableFeatures(Navin_hbca_c106, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c106)
Navin_hbca_c106 <- ScaleData(Navin_hbca_c106, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c106 <- RunPCA(Navin_hbca_c106, features = VariableFeatures(object = Navin_hbca_c106))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c106[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c106, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c106, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c106 <- JackStraw(Navin_hbca_c106, num.replicate = 100)
Navin_hbca_c106 <- ScoreJackStraw(Navin_hbca_c106, dims = 1:20)
ElbowPlot(Navin_hbca_c106)

# Cluster cells
Navin_hbca_c106 <- FindNeighbors(Navin_hbca_c106, dims = 1:20)
Navin_hbca_c106 <- FindClusters(Navin_hbca_c106, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c106 <- RunUMAP(Navin_hbca_c106, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c106, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c106, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c106.rds")

# Remove Object
rm(Navin_hbca_c106)



##########################################################
########################### Navin_hbca_c107 ##############
##########################################################

# Load the dataset
Navin_hbca_c107 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c107/GSM7500416_hbca_c107_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c107 <- CreateSeuratObject(counts = Navin_hbca_c107, project = "Navin_hbca_c107", min.cells = 3, min.features = 200)
Navin_hbca_c107

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c107[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c107, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c107, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c107, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c107, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c107 <- subset(Navin_hbca_c107, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c107 <- NormalizeData(Navin_hbca_c107, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c107 <- FindVariableFeatures(Navin_hbca_c107, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c107)
Navin_hbca_c107 <- ScaleData(Navin_hbca_c107, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c107 <- RunPCA(Navin_hbca_c107, features = VariableFeatures(object = Navin_hbca_c107))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c107[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c107, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c107, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c107 <- JackStraw(Navin_hbca_c107, num.replicate = 100)
Navin_hbca_c107 <- ScoreJackStraw(Navin_hbca_c107, dims = 1:20)
ElbowPlot(Navin_hbca_c107)

# Cluster cells
Navin_hbca_c107 <- FindNeighbors(Navin_hbca_c107, dims = 1:20)
Navin_hbca_c107 <- FindClusters(Navin_hbca_c107, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c107 <- RunUMAP(Navin_hbca_c107, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c107, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c107, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c107.rds")

# Remove Object
rm(Navin_hbca_c107)



##########################################################
########################### Navin_hbca_c108 ##############
##########################################################

# Load the dataset
Navin_hbca_c108 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c108/GSM7500417_hbca_c108_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c108 <- CreateSeuratObject(counts = Navin_hbca_c108, project = "Navin_hbca_c108", min.cells = 3, min.features = 200)
Navin_hbca_c108

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c108[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c108, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c108, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c108, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c108, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c108 <- subset(Navin_hbca_c108, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c108 <- NormalizeData(Navin_hbca_c108, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c108 <- FindVariableFeatures(Navin_hbca_c108, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c108)
Navin_hbca_c108 <- ScaleData(Navin_hbca_c108, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c108 <- RunPCA(Navin_hbca_c108, features = VariableFeatures(object = Navin_hbca_c108))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c108[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c108, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c108, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c108 <- JackStraw(Navin_hbca_c108, num.replicate = 100)
Navin_hbca_c108 <- ScoreJackStraw(Navin_hbca_c108, dims = 1:20)
ElbowPlot(Navin_hbca_c108)

# Cluster cells
Navin_hbca_c108 <- FindNeighbors(Navin_hbca_c108, dims = 1:20)
Navin_hbca_c108 <- FindClusters(Navin_hbca_c108, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c108 <- RunUMAP(Navin_hbca_c108, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c108, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c108, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c108.rds")

# Remove Object
rm(Navin_hbca_c108)



##########################################################
########################### Navin_hbca_c109 ##############
##########################################################

# Load the dataset
Navin_hbca_c109 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c109/GSM7500418_hbca_c109_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c109 <- CreateSeuratObject(counts = Navin_hbca_c109, project = "Navin_hbca_c109", min.cells = 3, min.features = 200)
Navin_hbca_c109

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c109[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c109, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c109, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c109, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c109, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c109 <- subset(Navin_hbca_c109, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c109 <- NormalizeData(Navin_hbca_c109, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c109 <- FindVariableFeatures(Navin_hbca_c109, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c109)
Navin_hbca_c109 <- ScaleData(Navin_hbca_c109, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c109 <- RunPCA(Navin_hbca_c109, features = VariableFeatures(object = Navin_hbca_c109))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c109[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c109, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c109, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c109 <- JackStraw(Navin_hbca_c109, num.replicate = 100)
Navin_hbca_c109 <- ScoreJackStraw(Navin_hbca_c109, dims = 1:20)
ElbowPlot(Navin_hbca_c109)

# Cluster cells
Navin_hbca_c109 <- FindNeighbors(Navin_hbca_c109, dims = 1:20)
Navin_hbca_c109 <- FindClusters(Navin_hbca_c109, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c109 <- RunUMAP(Navin_hbca_c109, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c109, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c109, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c109.rds")

# Remove Object
rm(Navin_hbca_c109)



##########################################################
########################### Navin_hbca_c110 ##############
##########################################################

# Load the dataset
Navin_hbca_c110 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c110/GSM7500419_hbca_c110_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c110 <- CreateSeuratObject(counts = Navin_hbca_c110, project = "Navin_hbca_c110", min.cells = 3, min.features = 200)
Navin_hbca_c110

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c110[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c110, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c110, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c110, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c110, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c110 <- subset(Navin_hbca_c110, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c110 <- NormalizeData(Navin_hbca_c110, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c110 <- FindVariableFeatures(Navin_hbca_c110, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c110)
Navin_hbca_c110 <- ScaleData(Navin_hbca_c110, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c110 <- RunPCA(Navin_hbca_c110, features = VariableFeatures(object = Navin_hbca_c110))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c110[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c110, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c110, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c110 <- JackStraw(Navin_hbca_c110, num.replicate = 100)
Navin_hbca_c110 <- ScoreJackStraw(Navin_hbca_c110, dims = 1:20)
ElbowPlot(Navin_hbca_c110)

# Cluster cells
Navin_hbca_c110 <- FindNeighbors(Navin_hbca_c110, dims = 1:20)
Navin_hbca_c110 <- FindClusters(Navin_hbca_c110, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c110 <- RunUMAP(Navin_hbca_c110, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c110, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c110, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c110.rds")

# Remove Object
rm(Navin_hbca_c110)



##########################################################
########################### Navin_hbca_c111 ##############
##########################################################

# Load the dataset
Navin_hbca_c111 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c111/GSM7500420_hbca_c111_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c111 <- CreateSeuratObject(counts = Navin_hbca_c111, project = "Navin_hbca_c111", min.cells = 3, min.features = 200)
Navin_hbca_c111

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c111[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c111, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c111, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c111, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c111, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c111 <- subset(Navin_hbca_c111, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c111 <- NormalizeData(Navin_hbca_c111, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c111 <- FindVariableFeatures(Navin_hbca_c111, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c111)
Navin_hbca_c111 <- ScaleData(Navin_hbca_c111, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c111 <- RunPCA(Navin_hbca_c111, features = VariableFeatures(object = Navin_hbca_c111))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c111[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c111, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c111, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c111 <- JackStraw(Navin_hbca_c111, num.replicate = 100)
Navin_hbca_c111 <- ScoreJackStraw(Navin_hbca_c111, dims = 1:20)
ElbowPlot(Navin_hbca_c111)

# Cluster cells
Navin_hbca_c111 <- FindNeighbors(Navin_hbca_c111, dims = 1:20)
Navin_hbca_c111 <- FindClusters(Navin_hbca_c111, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c111 <- RunUMAP(Navin_hbca_c111, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c111, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c111, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c111.rds")

# Remove Object
rm(Navin_hbca_c111)



##########################################################
########################### Navin_hbca_c113 ##############
##########################################################

# Load the dataset
Navin_hbca_c113 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c113/GSM7500422_hbca_c113_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c113 <- CreateSeuratObject(counts = Navin_hbca_c113, project = "Navin_hbca_c113", min.cells = 3, min.features = 200)
Navin_hbca_c113

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c113[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c113, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c113, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c113, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c113, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c113 <- subset(Navin_hbca_c113, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c113 <- NormalizeData(Navin_hbca_c113, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c113 <- FindVariableFeatures(Navin_hbca_c113, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c113)
Navin_hbca_c113 <- ScaleData(Navin_hbca_c113, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c113 <- RunPCA(Navin_hbca_c113, features = VariableFeatures(object = Navin_hbca_c113))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c113[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c113, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c113, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c113 <- JackStraw(Navin_hbca_c113, num.replicate = 100)
Navin_hbca_c113 <- ScoreJackStraw(Navin_hbca_c113, dims = 1:20)
ElbowPlot(Navin_hbca_c113)

# Cluster cells
Navin_hbca_c113 <- FindNeighbors(Navin_hbca_c113, dims = 1:20)
Navin_hbca_c113 <- FindClusters(Navin_hbca_c113, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c113 <- RunUMAP(Navin_hbca_c113, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c113, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c113, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c113.rds")

# Remove Object
rm(Navin_hbca_c113)



##########################################################
########################### Navin_hbca_c114 ##############
##########################################################

# Load the dataset
Navin_hbca_c114 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c114/GSM7500423_hbca_c114_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c114 <- CreateSeuratObject(counts = Navin_hbca_c114, project = "Navin_hbca_c114", min.cells = 3, min.features = 200)
Navin_hbca_c114

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c114[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c114, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c114, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c114, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c114, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c114 <- subset(Navin_hbca_c114, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c114 <- NormalizeData(Navin_hbca_c114, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c114 <- FindVariableFeatures(Navin_hbca_c114, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c114)
Navin_hbca_c114 <- ScaleData(Navin_hbca_c114, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c114 <- RunPCA(Navin_hbca_c114, features = VariableFeatures(object = Navin_hbca_c114))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c114[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c114, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c114, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c114 <- JackStraw(Navin_hbca_c114, num.replicate = 100)
Navin_hbca_c114 <- ScoreJackStraw(Navin_hbca_c114, dims = 1:20)
ElbowPlot(Navin_hbca_c114)

# Cluster cells
Navin_hbca_c114 <- FindNeighbors(Navin_hbca_c114, dims = 1:20)
Navin_hbca_c114 <- FindClusters(Navin_hbca_c114, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c114 <- RunUMAP(Navin_hbca_c114, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c114, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c114, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c114.rds")

# Remove Object
rm(Navin_hbca_c114)



##########################################################
########################### Navin_hbca_c115 ##############
##########################################################

# Load the dataset
Navin_hbca_c115 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c115/GSM7500424_hbca_c115_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c115 <- CreateSeuratObject(counts = Navin_hbca_c115, project = "Navin_hbca_c115", min.cells = 3, min.features = 200)
Navin_hbca_c115

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c115[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c115, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c115, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c115, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c115, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c115 <- subset(Navin_hbca_c115, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c115 <- NormalizeData(Navin_hbca_c115, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c115 <- FindVariableFeatures(Navin_hbca_c115, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c115)
Navin_hbca_c115 <- ScaleData(Navin_hbca_c115, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c115 <- RunPCA(Navin_hbca_c115, features = VariableFeatures(object = Navin_hbca_c115))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c115[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c115, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c115, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c115 <- JackStraw(Navin_hbca_c115, num.replicate = 100)
Navin_hbca_c115 <- ScoreJackStraw(Navin_hbca_c115, dims = 1:20)
ElbowPlot(Navin_hbca_c115)

# Cluster cells
Navin_hbca_c115 <- FindNeighbors(Navin_hbca_c115, dims = 1:20)
Navin_hbca_c115 <- FindClusters(Navin_hbca_c115, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c115 <- RunUMAP(Navin_hbca_c115, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c115, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c115, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c115.rds")

# Remove Object
rm(Navin_hbca_c115)



##########################################################
########################### Navin_hbca_c118 ##############
##########################################################

# Load the dataset
Navin_hbca_c118 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c118/GSM7500427_hbca_c118_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c118 <- CreateSeuratObject(counts = Navin_hbca_c118, project = "Navin_hbca_c118", min.cells = 3, min.features = 200)
Navin_hbca_c118

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c118[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c118, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c118, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c118, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c118, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c118 <- subset(Navin_hbca_c118, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c118 <- NormalizeData(Navin_hbca_c118, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c118 <- FindVariableFeatures(Navin_hbca_c118, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c118)
Navin_hbca_c118 <- ScaleData(Navin_hbca_c118, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c118 <- RunPCA(Navin_hbca_c118, features = VariableFeatures(object = Navin_hbca_c118))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c118[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c118, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c118, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c118 <- JackStraw(Navin_hbca_c118, num.replicate = 100)
Navin_hbca_c118 <- ScoreJackStraw(Navin_hbca_c118, dims = 1:20)
ElbowPlot(Navin_hbca_c118)

# Cluster cells
Navin_hbca_c118 <- FindNeighbors(Navin_hbca_c118, dims = 1:20)
Navin_hbca_c118 <- FindClusters(Navin_hbca_c118, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c118 <- RunUMAP(Navin_hbca_c118, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c118, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c118, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c118.rds")

# Remove Object
rm(Navin_hbca_c118)



##########################################################
########################### Navin_hbca_c119 ##############
##########################################################

# Load the dataset
Navin_hbca_c119 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c119/GSM7500428_hbca_c119_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c119 <- CreateSeuratObject(counts = Navin_hbca_c119, project = "Navin_hbca_c119", min.cells = 3, min.features = 200)
Navin_hbca_c119

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c119[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c119, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c119, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c119, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c119, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c119 <- subset(Navin_hbca_c119, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c119 <- NormalizeData(Navin_hbca_c119, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c119 <- FindVariableFeatures(Navin_hbca_c119, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c119)
Navin_hbca_c119 <- ScaleData(Navin_hbca_c119, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c119 <- RunPCA(Navin_hbca_c119, features = VariableFeatures(object = Navin_hbca_c119))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c119[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c119, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c119, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c119 <- JackStraw(Navin_hbca_c119, num.replicate = 100)
Navin_hbca_c119 <- ScoreJackStraw(Navin_hbca_c119, dims = 1:20)
ElbowPlot(Navin_hbca_c119)

# Cluster cells
Navin_hbca_c119 <- FindNeighbors(Navin_hbca_c119, dims = 1:20)
Navin_hbca_c119 <- FindClusters(Navin_hbca_c119, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c119 <- RunUMAP(Navin_hbca_c119, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c119, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c119, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c119.rds")

# Remove Object
rm(Navin_hbca_c119)



##########################################################
########################### Navin_hbca_c121 ##############
##########################################################

# Load the dataset
Navin_hbca_c121 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c121/GSM7500430_hbca_c121_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c121 <- CreateSeuratObject(counts = Navin_hbca_c121, project = "Navin_hbca_c121", min.cells = 3, min.features = 200)
Navin_hbca_c121

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c121[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c121, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c121, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c121, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c121, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c121 <- subset(Navin_hbca_c121, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c121 <- NormalizeData(Navin_hbca_c121, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c121 <- FindVariableFeatures(Navin_hbca_c121, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c121)
Navin_hbca_c121 <- ScaleData(Navin_hbca_c121, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c121 <- RunPCA(Navin_hbca_c121, features = VariableFeatures(object = Navin_hbca_c121))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c121[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c121, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c121, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c121 <- JackStraw(Navin_hbca_c121, num.replicate = 100)
Navin_hbca_c121 <- ScoreJackStraw(Navin_hbca_c121, dims = 1:20)
ElbowPlot(Navin_hbca_c121)

# Cluster cells
Navin_hbca_c121 <- FindNeighbors(Navin_hbca_c121, dims = 1:20)
Navin_hbca_c121 <- FindClusters(Navin_hbca_c121, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c121 <- RunUMAP(Navin_hbca_c121, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c121, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c121, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c121.rds")

# Remove Object
rm(Navin_hbca_c121)



##########################################################
########################### Navin_hbca_c122 ##############
##########################################################

# Load the dataset
Navin_hbca_c122 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c122/GSM7500431_hbca_c122_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c122 <- CreateSeuratObject(counts = Navin_hbca_c122, project = "Navin_hbca_c122", min.cells = 3, min.features = 200)
Navin_hbca_c122

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c122[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c122, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c122, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c122, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c122, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c122 <- subset(Navin_hbca_c122, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c122 <- NormalizeData(Navin_hbca_c122, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c122 <- FindVariableFeatures(Navin_hbca_c122, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c122)
Navin_hbca_c122 <- ScaleData(Navin_hbca_c122, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c122 <- RunPCA(Navin_hbca_c122, features = VariableFeatures(object = Navin_hbca_c122))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c122[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c122, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c122, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c122 <- JackStraw(Navin_hbca_c122, num.replicate = 100)
Navin_hbca_c122 <- ScoreJackStraw(Navin_hbca_c122, dims = 1:20)
ElbowPlot(Navin_hbca_c122)

# Cluster cells
Navin_hbca_c122 <- FindNeighbors(Navin_hbca_c122, dims = 1:20)
Navin_hbca_c122 <- FindClusters(Navin_hbca_c122, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c122 <- RunUMAP(Navin_hbca_c122, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c122, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c122, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c122.rds")

# Remove Object
rm(Navin_hbca_c122)



##########################################################
########################### Navin_hbca_c123 ##############
##########################################################

# Load the dataset
Navin_hbca_c123 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c123/GSM7500432_hbca_c123_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c123 <- CreateSeuratObject(counts = Navin_hbca_c123, project = "Navin_hbca_c123", min.cells = 3, min.features = 200)
Navin_hbca_c123

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c123[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c123, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c123, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c123, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c123, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c123 <- subset(Navin_hbca_c123, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c123 <- NormalizeData(Navin_hbca_c123, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c123 <- FindVariableFeatures(Navin_hbca_c123, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c123)
Navin_hbca_c123 <- ScaleData(Navin_hbca_c123, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c123 <- RunPCA(Navin_hbca_c123, features = VariableFeatures(object = Navin_hbca_c123))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c123[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c123, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c123, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c123 <- JackStraw(Navin_hbca_c123, num.replicate = 100)
Navin_hbca_c123 <- ScoreJackStraw(Navin_hbca_c123, dims = 1:20)
ElbowPlot(Navin_hbca_c123)

# Cluster cells
Navin_hbca_c123 <- FindNeighbors(Navin_hbca_c123, dims = 1:20)
Navin_hbca_c123 <- FindClusters(Navin_hbca_c123, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c123 <- RunUMAP(Navin_hbca_c123, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c123, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c123, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c123.rds")

# Remove Object
rm(Navin_hbca_c123)



##########################################################
########################### Navin_hbca_c124 ##############
##########################################################

# Load the dataset
Navin_hbca_c124 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c124/GSM7500433_hbca_c124_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c124 <- CreateSeuratObject(counts = Navin_hbca_c124, project = "Navin_hbca_c124", min.cells = 3, min.features = 200)
Navin_hbca_c124

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c124[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c124, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c124, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c124, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c124, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c124 <- subset(Navin_hbca_c124, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c124 <- NormalizeData(Navin_hbca_c124, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c124 <- FindVariableFeatures(Navin_hbca_c124, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c124)
Navin_hbca_c124 <- ScaleData(Navin_hbca_c124, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c124 <- RunPCA(Navin_hbca_c124, features = VariableFeatures(object = Navin_hbca_c124))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c124[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c124, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c124, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c124 <- JackStraw(Navin_hbca_c124, num.replicate = 100)
Navin_hbca_c124 <- ScoreJackStraw(Navin_hbca_c124, dims = 1:20)
ElbowPlot(Navin_hbca_c124)

# Cluster cells
Navin_hbca_c124 <- FindNeighbors(Navin_hbca_c124, dims = 1:20)
Navin_hbca_c124 <- FindClusters(Navin_hbca_c124, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c124 <- RunUMAP(Navin_hbca_c124, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c124, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c124, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c124.rds")

# Remove Object
rm(Navin_hbca_c124)



##########################################################
########################### Navin_hbca_c125 ##############
##########################################################

# Load the dataset
Navin_hbca_c125 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c125/GSM7500434_hbca_c125_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c125 <- CreateSeuratObject(counts = Navin_hbca_c125, project = "Navin_hbca_c125", min.cells = 3, min.features = 200)
Navin_hbca_c125

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c125[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c125, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c125, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c125, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c125, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c125 <- subset(Navin_hbca_c125, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c125 <- NormalizeData(Navin_hbca_c125, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c125 <- FindVariableFeatures(Navin_hbca_c125, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c125)
Navin_hbca_c125 <- ScaleData(Navin_hbca_c125, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c125 <- RunPCA(Navin_hbca_c125, features = VariableFeatures(object = Navin_hbca_c125))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c125[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c125, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c125, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c125 <- JackStraw(Navin_hbca_c125, num.replicate = 100)
Navin_hbca_c125 <- ScoreJackStraw(Navin_hbca_c125, dims = 1:20)
ElbowPlot(Navin_hbca_c125)

# Cluster cells
Navin_hbca_c125 <- FindNeighbors(Navin_hbca_c125, dims = 1:20)
Navin_hbca_c125 <- FindClusters(Navin_hbca_c125, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c125 <- RunUMAP(Navin_hbca_c125, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c125, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c125, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c125.rds")

# Remove Object
rm(Navin_hbca_c125)



##########################################################
########################### Navin_hbca_c126 ##############
##########################################################

# Load the dataset
Navin_hbca_c126 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c126/GSM7500435_hbca_c126_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c126 <- CreateSeuratObject(counts = Navin_hbca_c126, project = "Navin_hbca_c126", min.cells = 3, min.features = 200)
Navin_hbca_c126

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c126[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c126, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c126, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c126, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c126, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c126 <- subset(Navin_hbca_c126, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c126 <- NormalizeData(Navin_hbca_c126, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c126 <- FindVariableFeatures(Navin_hbca_c126, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c126)
Navin_hbca_c126 <- ScaleData(Navin_hbca_c126, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c126 <- RunPCA(Navin_hbca_c126, features = VariableFeatures(object = Navin_hbca_c126))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c126[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c126, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c126, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c126 <- JackStraw(Navin_hbca_c126, num.replicate = 100)
Navin_hbca_c126 <- ScoreJackStraw(Navin_hbca_c126, dims = 1:20)
ElbowPlot(Navin_hbca_c126)

# Cluster cells
Navin_hbca_c126 <- FindNeighbors(Navin_hbca_c126, dims = 1:20)
Navin_hbca_c126 <- FindClusters(Navin_hbca_c126, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c126 <- RunUMAP(Navin_hbca_c126, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c126, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c126, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c126.rds")

# Remove Object
rm(Navin_hbca_c126)



##########################################################
########################### Navin_hbca_c127 ##############
##########################################################

# Load the dataset
Navin_hbca_c127 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c127/GSM7500436_hbca_c127_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c127 <- CreateSeuratObject(counts = Navin_hbca_c127, project = "Navin_hbca_c127", min.cells = 3, min.features = 200)
Navin_hbca_c127

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c127[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c127, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c127, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c127, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c127, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c127 <- subset(Navin_hbca_c127, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c127 <- NormalizeData(Navin_hbca_c127, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c127 <- FindVariableFeatures(Navin_hbca_c127, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c127)
Navin_hbca_c127 <- ScaleData(Navin_hbca_c127, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c127 <- RunPCA(Navin_hbca_c127, features = VariableFeatures(object = Navin_hbca_c127))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c127[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c127, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c127, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c127 <- JackStraw(Navin_hbca_c127, num.replicate = 100)
Navin_hbca_c127 <- ScoreJackStraw(Navin_hbca_c127, dims = 1:20)
ElbowPlot(Navin_hbca_c127)

# Cluster cells
Navin_hbca_c127 <- FindNeighbors(Navin_hbca_c127, dims = 1:20)
Navin_hbca_c127 <- FindClusters(Navin_hbca_c127, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c127 <- RunUMAP(Navin_hbca_c127, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c127, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c127, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c127.rds")

# Remove Object
rm(Navin_hbca_c127)



##########################################################
########################### Navin_hbca_c128 ##############
##########################################################

# Load the dataset
Navin_hbca_c128 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c128/GSM7500437_hbca_c128_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c128 <- CreateSeuratObject(counts = Navin_hbca_c128, project = "Navin_hbca_c128", min.cells = 3, min.features = 200)
Navin_hbca_c128

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c128[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c128, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c128, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c128, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c128, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c128 <- subset(Navin_hbca_c128, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c128 <- NormalizeData(Navin_hbca_c128, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c128 <- FindVariableFeatures(Navin_hbca_c128, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c128)
Navin_hbca_c128 <- ScaleData(Navin_hbca_c128, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c128 <- RunPCA(Navin_hbca_c128, features = VariableFeatures(object = Navin_hbca_c128))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c128[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c128, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c128, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c128 <- JackStraw(Navin_hbca_c128, num.replicate = 100)
Navin_hbca_c128 <- ScoreJackStraw(Navin_hbca_c128, dims = 1:20)
ElbowPlot(Navin_hbca_c128)

# Cluster cells
Navin_hbca_c128 <- FindNeighbors(Navin_hbca_c128, dims = 1:20)
Navin_hbca_c128 <- FindClusters(Navin_hbca_c128, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c128 <- RunUMAP(Navin_hbca_c128, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c128, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c128, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c128.rds")

# Remove Object
rm(Navin_hbca_c128)



##########################################################
########################### Navin_hbca_c129 ##############
##########################################################

# Load the dataset
Navin_hbca_c129 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c129/GSM7500438_hbca_c129_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c129 <- CreateSeuratObject(counts = Navin_hbca_c129, project = "Navin_hbca_c129", min.cells = 3, min.features = 200)
Navin_hbca_c129

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c129[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c129, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c129, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c129, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c129, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c129 <- subset(Navin_hbca_c129, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c129 <- NormalizeData(Navin_hbca_c129, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c129 <- FindVariableFeatures(Navin_hbca_c129, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c129)
Navin_hbca_c129 <- ScaleData(Navin_hbca_c129, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c129 <- RunPCA(Navin_hbca_c129, features = VariableFeatures(object = Navin_hbca_c129))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c129[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c129, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c129, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c129 <- JackStraw(Navin_hbca_c129, num.replicate = 100)
Navin_hbca_c129 <- ScoreJackStraw(Navin_hbca_c129, dims = 1:20)
ElbowPlot(Navin_hbca_c129)

# Cluster cells
Navin_hbca_c129 <- FindNeighbors(Navin_hbca_c129, dims = 1:20)
Navin_hbca_c129 <- FindClusters(Navin_hbca_c129, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c129 <- RunUMAP(Navin_hbca_c129, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c129, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c129, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c129.rds")

# Remove Object
rm(Navin_hbca_c129)



##########################################################
########################### Navin_hbca_c130 ##############
##########################################################

# Load the dataset
Navin_hbca_c130 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c130/GSM7500439_hbca_c130_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c130 <- CreateSeuratObject(counts = Navin_hbca_c130, project = "Navin_hbca_c130", min.cells = 3, min.features = 200)
Navin_hbca_c130

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c130[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c130, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c130, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c130, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c130, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c130 <- subset(Navin_hbca_c130, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c130 <- NormalizeData(Navin_hbca_c130, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c130 <- FindVariableFeatures(Navin_hbca_c130, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c130)
Navin_hbca_c130 <- ScaleData(Navin_hbca_c130, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c130 <- RunPCA(Navin_hbca_c130, features = VariableFeatures(object = Navin_hbca_c130))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c130[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c130, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c130, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c130 <- JackStraw(Navin_hbca_c130, num.replicate = 100)
Navin_hbca_c130 <- ScoreJackStraw(Navin_hbca_c130, dims = 1:20)
ElbowPlot(Navin_hbca_c130)

# Cluster cells
Navin_hbca_c130 <- FindNeighbors(Navin_hbca_c130, dims = 1:20)
Navin_hbca_c130 <- FindClusters(Navin_hbca_c130, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c130 <- RunUMAP(Navin_hbca_c130, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c130, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c130, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c130.rds")

# Remove Object
rm(Navin_hbca_c130)



##########################################################
########################### Navin_hbca_c131 ##############
##########################################################

# Load the dataset
Navin_hbca_c131 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c131/GSM7500440_hbca_c131_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c131 <- CreateSeuratObject(counts = Navin_hbca_c131, project = "Navin_hbca_c131", min.cells = 3, min.features = 200)
Navin_hbca_c131

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c131[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c131, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c131, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c131, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c131, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c131 <- subset(Navin_hbca_c131, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c131 <- NormalizeData(Navin_hbca_c131, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c131 <- FindVariableFeatures(Navin_hbca_c131, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c131)
Navin_hbca_c131 <- ScaleData(Navin_hbca_c131, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c131 <- RunPCA(Navin_hbca_c131, features = VariableFeatures(object = Navin_hbca_c131))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c131[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c131, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c131, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c131 <- JackStraw(Navin_hbca_c131, num.replicate = 100)
Navin_hbca_c131 <- ScoreJackStraw(Navin_hbca_c131, dims = 1:20)
ElbowPlot(Navin_hbca_c131)

# Cluster cells
Navin_hbca_c131 <- FindNeighbors(Navin_hbca_c131, dims = 1:20)
Navin_hbca_c131 <- FindClusters(Navin_hbca_c131, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c131 <- RunUMAP(Navin_hbca_c131, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c131, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c131, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c131.rds")

# Remove Object
rm(Navin_hbca_c131)



##########################################################
########################### Navin_hbca_c132 ##############
##########################################################

# Load the dataset
Navin_hbca_c132 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c132/GSM7500441_hbca_c132_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c132 <- CreateSeuratObject(counts = Navin_hbca_c132, project = "Navin_hbca_c132", min.cells = 3, min.features = 200)
Navin_hbca_c132

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c132[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c132, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c132, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c132, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c132, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c132 <- subset(Navin_hbca_c132, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c132 <- NormalizeData(Navin_hbca_c132, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c132 <- FindVariableFeatures(Navin_hbca_c132, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c132)
Navin_hbca_c132 <- ScaleData(Navin_hbca_c132, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c132 <- RunPCA(Navin_hbca_c132, features = VariableFeatures(object = Navin_hbca_c132))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c132[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c132, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c132, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c132 <- JackStraw(Navin_hbca_c132, num.replicate = 100)
Navin_hbca_c132 <- ScoreJackStraw(Navin_hbca_c132, dims = 1:20)
ElbowPlot(Navin_hbca_c132)

# Cluster cells
Navin_hbca_c132 <- FindNeighbors(Navin_hbca_c132, dims = 1:20)
Navin_hbca_c132 <- FindClusters(Navin_hbca_c132, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c132 <- RunUMAP(Navin_hbca_c132, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c132, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c132, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c132.rds")

# Remove Object
rm(Navin_hbca_c132)



##########################################################
########################### Navin_hbca_c133 ##############
##########################################################

# Load the dataset
Navin_hbca_c133 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c133/GSM7500442_hbca_c133_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c133 <- CreateSeuratObject(counts = Navin_hbca_c133, project = "Navin_hbca_c133", min.cells = 3, min.features = 200)
Navin_hbca_c133

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c133[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c133, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c133, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c133, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c133, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c133 <- subset(Navin_hbca_c133, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c133 <- NormalizeData(Navin_hbca_c133, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c133 <- FindVariableFeatures(Navin_hbca_c133, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c133)
Navin_hbca_c133 <- ScaleData(Navin_hbca_c133, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c133 <- RunPCA(Navin_hbca_c133, features = VariableFeatures(object = Navin_hbca_c133))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c133[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c133, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c133, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c133 <- JackStraw(Navin_hbca_c133, num.replicate = 100)
Navin_hbca_c133 <- ScoreJackStraw(Navin_hbca_c133, dims = 1:20)
ElbowPlot(Navin_hbca_c133)

# Cluster cells
Navin_hbca_c133 <- FindNeighbors(Navin_hbca_c133, dims = 1:20)
Navin_hbca_c133 <- FindClusters(Navin_hbca_c133, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c133 <- RunUMAP(Navin_hbca_c133, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c133, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c133, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c133.rds")

# Remove Object
rm(Navin_hbca_c133)



##########################################################
########################### Navin_hbca_c134 ##############
##########################################################

# Load the dataset
Navin_hbca_c134 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c134/GSM7500443_hbca_c134_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c134 <- CreateSeuratObject(counts = Navin_hbca_c134, project = "Navin_hbca_c134", min.cells = 3, min.features = 200)
Navin_hbca_c134

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c134[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c134, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c134, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c134, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c134, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c134 <- subset(Navin_hbca_c134, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c134 <- NormalizeData(Navin_hbca_c134, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c134 <- FindVariableFeatures(Navin_hbca_c134, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c134)
Navin_hbca_c134 <- ScaleData(Navin_hbca_c134, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c134 <- RunPCA(Navin_hbca_c134, features = VariableFeatures(object = Navin_hbca_c134))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c134[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c134, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c134, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c134 <- JackStraw(Navin_hbca_c134, num.replicate = 100)
Navin_hbca_c134 <- ScoreJackStraw(Navin_hbca_c134, dims = 1:20)
ElbowPlot(Navin_hbca_c134)

# Cluster cells
Navin_hbca_c134 <- FindNeighbors(Navin_hbca_c134, dims = 1:20)
Navin_hbca_c134 <- FindClusters(Navin_hbca_c134, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c134 <- RunUMAP(Navin_hbca_c134, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c134, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c134, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c134.rds")

# Remove Object
rm(Navin_hbca_c134)



##########################################################
########################### Navin_hbca_c135 ##############
##########################################################

# Load the dataset
Navin_hbca_c135 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c135/GSM7500444_hbca_c135_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c135 <- CreateSeuratObject(counts = Navin_hbca_c135, project = "Navin_hbca_c135", min.cells = 3, min.features = 200)
Navin_hbca_c135

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c135[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c135, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c135, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c135, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c135, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c135 <- subset(Navin_hbca_c135, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c135 <- NormalizeData(Navin_hbca_c135, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c135 <- FindVariableFeatures(Navin_hbca_c135, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c135)
Navin_hbca_c135 <- ScaleData(Navin_hbca_c135, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c135 <- RunPCA(Navin_hbca_c135, features = VariableFeatures(object = Navin_hbca_c135))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c135[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c135, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c135, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c135 <- JackStraw(Navin_hbca_c135, num.replicate = 100)
Navin_hbca_c135 <- ScoreJackStraw(Navin_hbca_c135, dims = 1:20)
ElbowPlot(Navin_hbca_c135)

# Cluster cells
Navin_hbca_c135 <- FindNeighbors(Navin_hbca_c135, dims = 1:20)
Navin_hbca_c135 <- FindClusters(Navin_hbca_c135, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c135 <- RunUMAP(Navin_hbca_c135, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c135, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c135, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c135.rds")

# Remove Object
rm(Navin_hbca_c135)



##########################################################
########################### Navin_hbca_c136 ##############
##########################################################

# Load the dataset
Navin_hbca_c136 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c136/GSM7500445_hbca_c136_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c136 <- CreateSeuratObject(counts = Navin_hbca_c136, project = "Navin_hbca_c136", min.cells = 3, min.features = 200)
Navin_hbca_c136

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c136[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c136, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c136, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c136, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c136, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c136 <- subset(Navin_hbca_c136, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c136 <- NormalizeData(Navin_hbca_c136, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c136 <- FindVariableFeatures(Navin_hbca_c136, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c136)
Navin_hbca_c136 <- ScaleData(Navin_hbca_c136, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c136 <- RunPCA(Navin_hbca_c136, features = VariableFeatures(object = Navin_hbca_c136))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c136[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c136, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c136, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c136 <- JackStraw(Navin_hbca_c136, num.replicate = 100)
Navin_hbca_c136 <- ScoreJackStraw(Navin_hbca_c136, dims = 1:20)
ElbowPlot(Navin_hbca_c136)

# Cluster cells
Navin_hbca_c136 <- FindNeighbors(Navin_hbca_c136, dims = 1:20)
Navin_hbca_c136 <- FindClusters(Navin_hbca_c136, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c136 <- RunUMAP(Navin_hbca_c136, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c136, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c136, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c136.rds")

# Remove Object
rm(Navin_hbca_c136)



##########################################################
########################### Navin_hbca_c137 ##############
##########################################################

# Load the dataset
Navin_hbca_c137 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c137/GSM7500446_hbca_c137_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c137 <- CreateSeuratObject(counts = Navin_hbca_c137, project = "Navin_hbca_c137", min.cells = 3, min.features = 200)
Navin_hbca_c137

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c137[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c137, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c137, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c137, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c137, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c137 <- subset(Navin_hbca_c137, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c137 <- NormalizeData(Navin_hbca_c137, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c137 <- FindVariableFeatures(Navin_hbca_c137, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c137)
Navin_hbca_c137 <- ScaleData(Navin_hbca_c137, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c137 <- RunPCA(Navin_hbca_c137, features = VariableFeatures(object = Navin_hbca_c137))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c137[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c137, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c137, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c137 <- JackStraw(Navin_hbca_c137, num.replicate = 100)
Navin_hbca_c137 <- ScoreJackStraw(Navin_hbca_c137, dims = 1:20)
ElbowPlot(Navin_hbca_c137)

# Cluster cells
Navin_hbca_c137 <- FindNeighbors(Navin_hbca_c137, dims = 1:20)
Navin_hbca_c137 <- FindClusters(Navin_hbca_c137, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c137 <- RunUMAP(Navin_hbca_c137, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c137, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c137, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c137.rds")

# Remove Object
rm(Navin_hbca_c137)



##########################################################
########################### Navin_hbca_c138 ##############
##########################################################

# Load the dataset
Navin_hbca_c138 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c138/GSM7500447_hbca_c138_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c138 <- CreateSeuratObject(counts = Navin_hbca_c138, project = "Navin_hbca_c138", min.cells = 3, min.features = 200)
Navin_hbca_c138

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c138[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c138, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c138, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c138, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c138, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c138 <- subset(Navin_hbca_c138, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c138 <- NormalizeData(Navin_hbca_c138, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c138 <- FindVariableFeatures(Navin_hbca_c138, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c138)
Navin_hbca_c138 <- ScaleData(Navin_hbca_c138, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c138 <- RunPCA(Navin_hbca_c138, features = VariableFeatures(object = Navin_hbca_c138))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c138[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c138, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c138, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c138 <- JackStraw(Navin_hbca_c138, num.replicate = 100)
Navin_hbca_c138 <- ScoreJackStraw(Navin_hbca_c138, dims = 1:20)
ElbowPlot(Navin_hbca_c138)

# Cluster cells
Navin_hbca_c138 <- FindNeighbors(Navin_hbca_c138, dims = 1:20)
Navin_hbca_c138 <- FindClusters(Navin_hbca_c138, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c138 <- RunUMAP(Navin_hbca_c138, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c138, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c138, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c138.rds")


# Remove Object
rm(Navin_hbca_c138)



##########################################################
########################### Navin_hbca_c139 ##############
##########################################################

# Load the dataset
Navin_hbca_c139 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c139/GSM7500448_hbca_c139_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c139 <- CreateSeuratObject(counts = Navin_hbca_c139, project = "Navin_hbca_c139", min.cells = 3, min.features = 200)
Navin_hbca_c139

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c139[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c139, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c139, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c139, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c139, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c139 <- subset(Navin_hbca_c139, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c139 <- NormalizeData(Navin_hbca_c139, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c139 <- FindVariableFeatures(Navin_hbca_c139, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c139)
Navin_hbca_c139 <- ScaleData(Navin_hbca_c139, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c139 <- RunPCA(Navin_hbca_c139, features = VariableFeatures(object = Navin_hbca_c139))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c139[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c139, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c139, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c139 <- JackStraw(Navin_hbca_c139, num.replicate = 100)
Navin_hbca_c139 <- ScoreJackStraw(Navin_hbca_c139, dims = 1:20)
ElbowPlot(Navin_hbca_c139)

# Cluster cells
Navin_hbca_c139 <- FindNeighbors(Navin_hbca_c139, dims = 1:20)
Navin_hbca_c139 <- FindClusters(Navin_hbca_c139, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c139 <- RunUMAP(Navin_hbca_c139, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c139, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c139, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c139.rds")

# Remove Object
rm(Navin_hbca_c139)



##########################################################
########################### Navin_hbca_c140 ##############
##########################################################

# Load the dataset
Navin_hbca_c140 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c140/GSM7500449_hbca_c140_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c140 <- CreateSeuratObject(counts = Navin_hbca_c140, project = "Navin_hbca_c140", min.cells = 3, min.features = 200)
Navin_hbca_c140

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c140[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c140, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c140, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c140, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c140, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c140 <- subset(Navin_hbca_c140, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c140 <- NormalizeData(Navin_hbca_c140, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c140 <- FindVariableFeatures(Navin_hbca_c140, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c140)
Navin_hbca_c140 <- ScaleData(Navin_hbca_c140, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c140 <- RunPCA(Navin_hbca_c140, features = VariableFeatures(object = Navin_hbca_c140))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c140[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c140, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c140, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c140 <- JackStraw(Navin_hbca_c140, num.replicate = 100)
Navin_hbca_c140 <- ScoreJackStraw(Navin_hbca_c140, dims = 1:20)
ElbowPlot(Navin_hbca_c140)

# Cluster cells
Navin_hbca_c140 <- FindNeighbors(Navin_hbca_c140, dims = 1:20)
Navin_hbca_c140 <- FindClusters(Navin_hbca_c140, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c140 <- RunUMAP(Navin_hbca_c140, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c140, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c140, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c140.rds")

# Remove Object
rm(Navin_hbca_c140)



##########################################################
########################### Navin_hbca_c141 ##############
##########################################################

# Load the dataset
Navin_hbca_c141 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c141/GSM7500450_hbca_c141_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c141 <- CreateSeuratObject(counts = Navin_hbca_c141, project = "Navin_hbca_c141", min.cells = 3, min.features = 200)
Navin_hbca_c141

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c141[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c141, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c141, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c141, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c141, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c141 <- subset(Navin_hbca_c141, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c141 <- NormalizeData(Navin_hbca_c141, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c141 <- FindVariableFeatures(Navin_hbca_c141, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c141)
Navin_hbca_c141 <- ScaleData(Navin_hbca_c141, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c141 <- RunPCA(Navin_hbca_c141, features = VariableFeatures(object = Navin_hbca_c141))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c141[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c141, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c141, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c141 <- JackStraw(Navin_hbca_c141, num.replicate = 100)
Navin_hbca_c141 <- ScoreJackStraw(Navin_hbca_c141, dims = 1:20)
ElbowPlot(Navin_hbca_c141)

# Cluster cells
Navin_hbca_c141 <- FindNeighbors(Navin_hbca_c141, dims = 1:20)
Navin_hbca_c141 <- FindClusters(Navin_hbca_c141, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c141 <- RunUMAP(Navin_hbca_c141, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c141, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c141, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c141.rds")

# Remove Object
rm(Navin_hbca_c141)



##########################################################
########################### Navin_hbca_c142 ##############
##########################################################

# Load the dataset
Navin_hbca_c142 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c142/GSM7500451_hbca_c142_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c142 <- CreateSeuratObject(counts = Navin_hbca_c142, project = "Navin_hbca_c142", min.cells = 3, min.features = 200)
Navin_hbca_c142

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c142[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c142, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c142, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c142, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c142, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c142 <- subset(Navin_hbca_c142, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c142 <- NormalizeData(Navin_hbca_c142, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c142 <- FindVariableFeatures(Navin_hbca_c142, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c142)
Navin_hbca_c142 <- ScaleData(Navin_hbca_c142, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c142 <- RunPCA(Navin_hbca_c142, features = VariableFeatures(object = Navin_hbca_c142))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c142[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c142, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c142, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c142 <- JackStraw(Navin_hbca_c142, num.replicate = 100)
Navin_hbca_c142 <- ScoreJackStraw(Navin_hbca_c142, dims = 1:20)
ElbowPlot(Navin_hbca_c142)

# Cluster cells
Navin_hbca_c142 <- FindNeighbors(Navin_hbca_c142, dims = 1:20)
Navin_hbca_c142 <- FindClusters(Navin_hbca_c142, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c142 <- RunUMAP(Navin_hbca_c142, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c142, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c142, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c142.rds")

# Remove Object
rm(Navin_hbca_c142)



##########################################################
########################### Navin_hbca_c143 ##############
##########################################################

# Load the dataset
Navin_hbca_c143 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c143/GSM7500452_hbca_c143_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c143 <- CreateSeuratObject(counts = Navin_hbca_c143, project = "Navin_hbca_c143", min.cells = 3, min.features = 200)
Navin_hbca_c143

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c143[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c143, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c143, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c143, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c143, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c143 <- subset(Navin_hbca_c143, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c143 <- NormalizeData(Navin_hbca_c143, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c143 <- FindVariableFeatures(Navin_hbca_c143, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c143)
Navin_hbca_c143 <- ScaleData(Navin_hbca_c143, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c143 <- RunPCA(Navin_hbca_c143, features = VariableFeatures(object = Navin_hbca_c143))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c143[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c143, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c143, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c143 <- JackStraw(Navin_hbca_c143, num.replicate = 100)
Navin_hbca_c143 <- ScoreJackStraw(Navin_hbca_c143, dims = 1:20)
ElbowPlot(Navin_hbca_c143)

# Cluster cells
Navin_hbca_c143 <- FindNeighbors(Navin_hbca_c143, dims = 1:20)
Navin_hbca_c143 <- FindClusters(Navin_hbca_c143, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c143 <- RunUMAP(Navin_hbca_c143, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c143, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c143, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c143.rds")

# Remove Object
rm(Navin_hbca_c143)



##########################################################
########################### Navin_hbca_c144 ##############
##########################################################

# Load the dataset
Navin_hbca_c144 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c144/GSM7500453_hbca_c144_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c144 <- CreateSeuratObject(counts = Navin_hbca_c144, project = "Navin_hbca_c144", min.cells = 3, min.features = 200)
Navin_hbca_c144

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c144[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c144, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c144, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c144, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c144, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c144 <- subset(Navin_hbca_c144, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c144 <- NormalizeData(Navin_hbca_c144, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c144 <- FindVariableFeatures(Navin_hbca_c144, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c144)
Navin_hbca_c144 <- ScaleData(Navin_hbca_c144, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c144 <- RunPCA(Navin_hbca_c144, features = VariableFeatures(object = Navin_hbca_c144))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c144[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c144, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c144, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c144 <- JackStraw(Navin_hbca_c144, num.replicate = 100)
Navin_hbca_c144 <- ScoreJackStraw(Navin_hbca_c144, dims = 1:20)
ElbowPlot(Navin_hbca_c144)

# Cluster cells
Navin_hbca_c144 <- FindNeighbors(Navin_hbca_c144, dims = 1:20)
Navin_hbca_c144 <- FindClusters(Navin_hbca_c144, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c144 <- RunUMAP(Navin_hbca_c144, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c144, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c144, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c144.rds")

# Remove Object
rm(Navin_hbca_c144)



##########################################################
########################### Navin_hbca_c145 ##############
##########################################################

# Load the dataset
Navin_hbca_c145 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c145/GSM7500454_hbca_c145_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c145 <- CreateSeuratObject(counts = Navin_hbca_c145, project = "Navin_hbca_c145", min.cells = 3, min.features = 200)
Navin_hbca_c145

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c145[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c145, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c145, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c145, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c145, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c145 <- subset(Navin_hbca_c145, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c145 <- NormalizeData(Navin_hbca_c145, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c145 <- FindVariableFeatures(Navin_hbca_c145, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c145)
Navin_hbca_c145 <- ScaleData(Navin_hbca_c145, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c145 <- RunPCA(Navin_hbca_c145, features = VariableFeatures(object = Navin_hbca_c145))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c145[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c145, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c145, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c145 <- JackStraw(Navin_hbca_c145, num.replicate = 100)
Navin_hbca_c145 <- ScoreJackStraw(Navin_hbca_c145, dims = 1:20)
ElbowPlot(Navin_hbca_c145)

# Cluster cells
Navin_hbca_c145 <- FindNeighbors(Navin_hbca_c145, dims = 1:20)
Navin_hbca_c145 <- FindClusters(Navin_hbca_c145, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c145 <- RunUMAP(Navin_hbca_c145, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c145, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c145, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c145.rds")

# Remove Object
rm(Navin_hbca_c145)



##########################################################
########################### Navin_hbca_c146 ##############
##########################################################

# Load the dataset
Navin_hbca_c146 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c146/GSM7500455_hbca_c146_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c146 <- CreateSeuratObject(counts = Navin_hbca_c146, project = "Navin_hbca_c146", min.cells = 3, min.features = 200)
Navin_hbca_c146

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c146[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c146, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c146, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c146, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c146, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c146 <- subset(Navin_hbca_c146, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c146 <- NormalizeData(Navin_hbca_c146, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c146 <- FindVariableFeatures(Navin_hbca_c146, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c146)
Navin_hbca_c146 <- ScaleData(Navin_hbca_c146, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c146 <- RunPCA(Navin_hbca_c146, features = VariableFeatures(object = Navin_hbca_c146))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c146[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c146, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c146, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c146 <- JackStraw(Navin_hbca_c146, num.replicate = 100)
Navin_hbca_c146 <- ScoreJackStraw(Navin_hbca_c146, dims = 1:20)
ElbowPlot(Navin_hbca_c146)

# Cluster cells
Navin_hbca_c146 <- FindNeighbors(Navin_hbca_c146, dims = 1:20)
Navin_hbca_c146 <- FindClusters(Navin_hbca_c146, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c146 <- RunUMAP(Navin_hbca_c146, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c146, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c146, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c146.rds")

# Remove Object
rm(Navin_hbca_c146)



##########################################################
########################### Navin_hbca_c147 ##############
##########################################################

# Load the dataset
Navin_hbca_c147 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c147/GSM7500456_hbca_c147_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c147 <- CreateSeuratObject(counts = Navin_hbca_c147, project = "Navin_hbca_c147", min.cells = 3, min.features = 200)
Navin_hbca_c147

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c147[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c147, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c147, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c147, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c147, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c147 <- subset(Navin_hbca_c147, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c147 <- NormalizeData(Navin_hbca_c147, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c147 <- FindVariableFeatures(Navin_hbca_c147, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c147)
Navin_hbca_c147 <- ScaleData(Navin_hbca_c147, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c147 <- RunPCA(Navin_hbca_c147, features = VariableFeatures(object = Navin_hbca_c147))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c147[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c147, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c147, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c147 <- JackStraw(Navin_hbca_c147, num.replicate = 100)
Navin_hbca_c147 <- ScoreJackStraw(Navin_hbca_c147, dims = 1:20)
ElbowPlot(Navin_hbca_c147)

# Cluster cells
Navin_hbca_c147 <- FindNeighbors(Navin_hbca_c147, dims = 1:20)
Navin_hbca_c147 <- FindClusters(Navin_hbca_c147, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c147 <- RunUMAP(Navin_hbca_c147, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c147, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c147, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c147.rds")

# Remove Object
rm(Navin_hbca_c147)



##########################################################
########################### Navin_hbca_c148 ##############
##########################################################

# Load the dataset
Navin_hbca_c148 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c148/GSM7500457_hbca_c148_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c148 <- CreateSeuratObject(counts = Navin_hbca_c148, project = "Navin_hbca_c148", min.cells = 3, min.features = 200)
Navin_hbca_c148

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c148[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c148, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c148, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c148, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c148, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c148 <- subset(Navin_hbca_c148, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c148 <- NormalizeData(Navin_hbca_c148, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c148 <- FindVariableFeatures(Navin_hbca_c148, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c148)
Navin_hbca_c148 <- ScaleData(Navin_hbca_c148, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c148 <- RunPCA(Navin_hbca_c148, features = VariableFeatures(object = Navin_hbca_c148))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c148[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c148, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c148, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c148 <- JackStraw(Navin_hbca_c148, num.replicate = 100)
Navin_hbca_c148 <- ScoreJackStraw(Navin_hbca_c148, dims = 1:20)
ElbowPlot(Navin_hbca_c148)

# Cluster cells
Navin_hbca_c148 <- FindNeighbors(Navin_hbca_c148, dims = 1:20)
Navin_hbca_c148 <- FindClusters(Navin_hbca_c148, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c148 <- RunUMAP(Navin_hbca_c148, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c148, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c148, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c148.rds")

# Remove Object
rm(Navin_hbca_c148)



##########################################################
########################### Navin_hbca_c149 ##############
##########################################################

# Load the dataset
Navin_hbca_c149 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c149/GSM7500458_hbca_c149_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c149 <- CreateSeuratObject(counts = Navin_hbca_c149, project = "Navin_hbca_c149", min.cells = 3, min.features = 200)
Navin_hbca_c149

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c149[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c149, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c149, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c149, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c149, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c149 <- subset(Navin_hbca_c149, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c149 <- NormalizeData(Navin_hbca_c149, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c149 <- FindVariableFeatures(Navin_hbca_c149, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c149)
Navin_hbca_c149 <- ScaleData(Navin_hbca_c149, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c149 <- RunPCA(Navin_hbca_c149, features = VariableFeatures(object = Navin_hbca_c149))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c149[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c149, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c149, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c149 <- JackStraw(Navin_hbca_c149, num.replicate = 100)
Navin_hbca_c149 <- ScoreJackStraw(Navin_hbca_c149, dims = 1:20)
ElbowPlot(Navin_hbca_c149)

# Cluster cells
Navin_hbca_c149 <- FindNeighbors(Navin_hbca_c149, dims = 1:20)
Navin_hbca_c149 <- FindClusters(Navin_hbca_c149, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c149 <- RunUMAP(Navin_hbca_c149, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c149, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c149, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c149.rds")

# Remove Object
rm(Navin_hbca_c149)



##########################################################
########################### Navin_hbca_c150 ##############
##########################################################

# Load the dataset
Navin_hbca_c150 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c150/GSM7500459_hbca_c150_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c150 <- CreateSeuratObject(counts = Navin_hbca_c150, project = "Navin_hbca_c150", min.cells = 3, min.features = 200)
Navin_hbca_c150

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c150[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c150, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c150, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c150, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c150, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c150 <- subset(Navin_hbca_c150, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c150 <- NormalizeData(Navin_hbca_c150, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c150 <- FindVariableFeatures(Navin_hbca_c150, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c150)
Navin_hbca_c150 <- ScaleData(Navin_hbca_c150, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c150 <- RunPCA(Navin_hbca_c150, features = VariableFeatures(object = Navin_hbca_c150))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c150[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c150, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c150, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c150 <- JackStraw(Navin_hbca_c150, num.replicate = 100)
Navin_hbca_c150 <- ScoreJackStraw(Navin_hbca_c150, dims = 1:20)
ElbowPlot(Navin_hbca_c150)

# Cluster cells
Navin_hbca_c150 <- FindNeighbors(Navin_hbca_c150, dims = 1:20)
Navin_hbca_c150 <- FindClusters(Navin_hbca_c150, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c150 <- RunUMAP(Navin_hbca_c150, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c150, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c150, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c150.rds")

# Remove Object
rm(Navin_hbca_c150)



##########################################################
########################### Navin_hbca_c151 ##############
##########################################################

# Load the dataset
Navin_hbca_c151 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c151/GSM7500460_hbca_c151_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c151 <- CreateSeuratObject(counts = Navin_hbca_c151, project = "Navin_hbca_c151", min.cells = 3, min.features = 200)
Navin_hbca_c151

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c151[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c151, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c151, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c151, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c151, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c151 <- subset(Navin_hbca_c151, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c151 <- NormalizeData(Navin_hbca_c151, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c151 <- FindVariableFeatures(Navin_hbca_c151, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c151)
Navin_hbca_c151 <- ScaleData(Navin_hbca_c151, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c151 <- RunPCA(Navin_hbca_c151, features = VariableFeatures(object = Navin_hbca_c151))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c151[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c151, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c151, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c151 <- JackStraw(Navin_hbca_c151, num.replicate = 100)
Navin_hbca_c151 <- ScoreJackStraw(Navin_hbca_c151, dims = 1:20)
ElbowPlot(Navin_hbca_c151)

# Cluster cells
Navin_hbca_c151 <- FindNeighbors(Navin_hbca_c151, dims = 1:20)
Navin_hbca_c151 <- FindClusters(Navin_hbca_c151, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c151 <- RunUMAP(Navin_hbca_c151, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c151, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c151, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c151.rds")

# Remove Object
rm(Navin_hbca_c151)



##########################################################
########################### Navin_hbca_c152 ##############
##########################################################

# Load the dataset
Navin_hbca_c152 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c152/GSM7500461_hbca_c152_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c152 <- CreateSeuratObject(counts = Navin_hbca_c152, project = "Navin_hbca_c152", min.cells = 3, min.features = 200)
Navin_hbca_c152

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c152[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c152, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c152, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c152, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c152, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c152 <- subset(Navin_hbca_c152, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c152 <- NormalizeData(Navin_hbca_c152, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c152 <- FindVariableFeatures(Navin_hbca_c152, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c152)
Navin_hbca_c152 <- ScaleData(Navin_hbca_c152, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c152 <- RunPCA(Navin_hbca_c152, features = VariableFeatures(object = Navin_hbca_c152))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c152[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c152, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c152, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c152 <- JackStraw(Navin_hbca_c152, num.replicate = 100)
Navin_hbca_c152 <- ScoreJackStraw(Navin_hbca_c152, dims = 1:20)
ElbowPlot(Navin_hbca_c152)

# Cluster cells
Navin_hbca_c152 <- FindNeighbors(Navin_hbca_c152, dims = 1:20)
Navin_hbca_c152 <- FindClusters(Navin_hbca_c152, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c152 <- RunUMAP(Navin_hbca_c152, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c152, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c152, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c152.rds")

# Remove Object
rm(Navin_hbca_c152)



##########################################################
########################### Navin_hbca_c153 ##############
##########################################################

# Load the dataset
Navin_hbca_c153 <- Read10X_h5("/R/R_Navin/Navin_Input/NORM/hbca_c153/GSM7500462_hbca_c153_filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
Navin_hbca_c153 <- CreateSeuratObject(counts = Navin_hbca_c153, project = "Navin_hbca_c153", min.cells = 3, min.features = 200)
Navin_hbca_c153

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Navin_hbca_c153[["percent.mt"]] <- PercentageFeatureSet(Navin_hbca_c153, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Navin_hbca_c153, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Navin_hbca_c153, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Navin_hbca_c153, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Navin_hbca_c153 <- subset(Navin_hbca_c153, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Navin_hbca_c153 <- NormalizeData(Navin_hbca_c153, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Navin_hbca_c153 <- FindVariableFeatures(Navin_hbca_c153, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Navin_hbca_c153)
Navin_hbca_c153 <- ScaleData(Navin_hbca_c153, features = all.genes)

# Perform linear dimensional reduction
Navin_hbca_c153 <- RunPCA(Navin_hbca_c153, features = VariableFeatures(object = Navin_hbca_c153))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c153[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c153, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c153, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Navin_hbca_c153 <- JackStraw(Navin_hbca_c153, num.replicate = 100)
Navin_hbca_c153 <- ScoreJackStraw(Navin_hbca_c153, dims = 1:20)
ElbowPlot(Navin_hbca_c153)

# Cluster cells
Navin_hbca_c153 <- FindNeighbors(Navin_hbca_c153, dims = 1:20)
Navin_hbca_c153 <- FindClusters(Navin_hbca_c153, resolution = 0.1)

# Run non-linear dimensional reduction
Navin_hbca_c153 <- RunUMAP(Navin_hbca_c153, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Navin_hbca_c153, reduction = "umap")

# Save RDS
saveRDS(Navin_hbca_c153, file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c153.rds")

# Remove Object
rm(Navin_hbca_c153)




#########################################################################################
##################### Step 2: Run DoubletFinder and Subset Singlets #####################
#########################################################################################

# Load Libraries
library(Seurat)
library(DoubletFinder)

# Assuming 9.6% doublet rate per sample based on averaged number of sequenced cells
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-

# Patient Samples
#hbca_c14
#hbca_c15
#hbca_c19
#hbca_c20
#hbca_c22
#hbca_c23
#hbca_c24
#hbca_c25
#hbca_c26
#hbca_c31
#hbca_c32
#hbca_c50
#hbca_c51
#hbca_c52
#hbca_c53
#hbca_c54
#hbca_c55
#hbca_c56
#hbca_c57
#hbca_c58
#hbca_c59
#hbca_c60
#hbca_c61
#hbca_c62
#hbca_c63
#hbca_c64
#hbca_c65
#hbca_c66 
#hbca_c67
#hbca_c68
#hbca_c69
#hbca_c70
#hbca_c71
#hbca_c72
#hbca_c73
#hbca_c74
#hbca_c75
#hbca_c76
#hbca_c77
#hbca_c78
#hbca_c79
#hbca_c80
#hbca_c81
#hbca_c82
#hbca_c83
#hbca_c84
#hbca_c85
#hbca_c86
#hbca_c87
#hbca_c88
#hbca_c89
#hbca_c90
#hbca_c91
#hbca_c92
#hbca_c93
#hbca_c94
#hbca_c95
#hbca_c96
#hbca_c97
#hbca_c98
#hbca_c99
#hbca_c100
#hbca_c101
#hbca_c102
#hbca_c103
#hbca_c104
#hbca_c105
#hbca_c106
#hbca_c107
#hbca_c108
#hbca_c109
#hbca_c110
#hbca_c111
#hbca_c112
#hbca_c113
#hbca_c114
#hbca_c117
#hbca_c118
#hbca_c120
#hbca_c121
#hbca_c122
#hbca_c123
#hbca_c124
#hbca_c125
#hbca_c126
#hbca_c127
#hbca_c128
#hbca_c129
#hbca_c130
#hbca_c131
#hbca_c132
#hbca_c133
#hbca_c134
#hbca_c135
#hbca_c136
#hbca_c137
#hbca_c138
#hbca_c139
#hbca_c140
#hbca_c141
#hbca_c142
#hbca_c143
#hbca_c144
#hbca_c145
#hbca_c146
#hbca_c147
#hbca_c148
#hbca_c149
#hbca_c150
#hbca_c151
#hbca_c152



################################################################################################
####################                  Reduction Mammoplasties                ###################
################################################################################################

################################################################################################
########################                Navin_hbca_c14                  ########################
################################################################################################


# Load Data
Navin_hbca_c14 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c14.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c14 <- paramSweep(Navin_hbca_c14, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c14 <- summarizeSweep(sweep.res.Navin_hbca_c14, GT = FALSE)
bcmvn_Navin_hbca_c14 <- find.pK(sweep.stats.Navin_hbca_c14)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c14[which.max(bcmvn_Navin_hbca_c14$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c14@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c14@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c14 <- doubletFinder(Navin_hbca_c14, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c14_Quant <- (Navin_hbca_c14@meta.data$DF.classification == "Singlet")
Navin_hbca_c14_Quant_Singlets <- length(Navin_hbca_c14_Quant[Navin_hbca_c14_Quant== TRUE])
Navin_hbca_c14_Quant_Doublets <- length(Navin_hbca_c14_Quant[Navin_hbca_c14_Quant== FALSE])
Navin_hbca_c14_Quant_Doublets_Percent <- Navin_hbca_c14_Quant_Doublets / (Navin_hbca_c14_Quant_Doublets + Navin_hbca_c14_Quant_Singlets) * 100
Navin_hbca_c14_Quant <- as.data.frame(c(Navin_hbca_c14_Quant_Singlets, Navin_hbca_c14_Quant_Doublets, Navin_hbca_c14_Quant_Doublets_Percent))
colnames(Navin_hbca_c14_Quant) <- c("Navin_hbca_c14_Quant")
rownames(Navin_hbca_c14_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c14_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c14_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c14_Singlets <- subset(Navin_hbca_c14, cells=rownames(Navin_hbca_c14@meta.data)[which(Navin_hbca_c14@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c14_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c14_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c14)
rm(Navin_hbca_c14_Singlets)
rm(Navin_hbca_c14_Quant)
rm(Navin_hbca_c14_Quant_Singlets)
rm(Navin_hbca_c14_Quant_Doublets)
rm(Navin_hbca_c14_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c14)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c14)
rm(sweep.stats.Navin_hbca_c14)
gc()



################################################################################################
########################                Navin_hbca_c15                  ########################
################################################################################################


# Load Data
Navin_hbca_c15 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c15.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c15 <- paramSweep(Navin_hbca_c15, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c15 <- summarizeSweep(sweep.res.Navin_hbca_c15, GT = FALSE)
bcmvn_Navin_hbca_c15 <- find.pK(sweep.stats.Navin_hbca_c15)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c15[which.max(bcmvn_Navin_hbca_c15$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c15@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c15@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c15 <- doubletFinder(Navin_hbca_c15, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c15_Quant <- (Navin_hbca_c15@meta.data$DF.classification == "Singlet")
Navin_hbca_c15_Quant_Singlets <- length(Navin_hbca_c15_Quant[Navin_hbca_c15_Quant== TRUE])
Navin_hbca_c15_Quant_Doublets <- length(Navin_hbca_c15_Quant[Navin_hbca_c15_Quant== FALSE])
Navin_hbca_c15_Quant_Doublets_Percent <- Navin_hbca_c15_Quant_Doublets / (Navin_hbca_c15_Quant_Doublets + Navin_hbca_c15_Quant_Singlets) * 100
Navin_hbca_c15_Quant <- as.data.frame(c(Navin_hbca_c15_Quant_Singlets, Navin_hbca_c15_Quant_Doublets, Navin_hbca_c15_Quant_Doublets_Percent))
colnames(Navin_hbca_c15_Quant) <- c("Navin_hbca_c15_Quant")
rownames(Navin_hbca_c15_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c15_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c15_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c15_Singlets <- subset(Navin_hbca_c15, cells=rownames(Navin_hbca_c15@meta.data)[which(Navin_hbca_c15@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c15_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c15_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c15)
rm(Navin_hbca_c15_Singlets)
rm(Navin_hbca_c15_Quant)
rm(Navin_hbca_c15_Quant_Singlets)
rm(Navin_hbca_c15_Quant_Doublets)
rm(Navin_hbca_c15_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c15)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c15)
rm(sweep.stats.Navin_hbca_c15)
gc()




################################################################################################
########################                Navin_hbca_c19                  ########################
################################################################################################


# Load Data
Navin_hbca_c19 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c19.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c19 <- paramSweep(Navin_hbca_c19, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c19 <- summarizeSweep(sweep.res.Navin_hbca_c19, GT = FALSE)
bcmvn_Navin_hbca_c19 <- find.pK(sweep.stats.Navin_hbca_c19)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c19[which.max(bcmvn_Navin_hbca_c19$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c19@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c19@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c19 <- doubletFinder(Navin_hbca_c19, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c19_Quant <- (Navin_hbca_c19@meta.data$DF.classification == "Singlet")
Navin_hbca_c19_Quant_Singlets <- length(Navin_hbca_c19_Quant[Navin_hbca_c19_Quant== TRUE])
Navin_hbca_c19_Quant_Doublets <- length(Navin_hbca_c19_Quant[Navin_hbca_c19_Quant== FALSE])
Navin_hbca_c19_Quant_Doublets_Percent <- Navin_hbca_c19_Quant_Doublets / (Navin_hbca_c19_Quant_Doublets + Navin_hbca_c19_Quant_Singlets) * 100
Navin_hbca_c19_Quant <- as.data.frame(c(Navin_hbca_c19_Quant_Singlets, Navin_hbca_c19_Quant_Doublets, Navin_hbca_c19_Quant_Doublets_Percent))
colnames(Navin_hbca_c19_Quant) <- c("Navin_hbca_c19_Quant")
rownames(Navin_hbca_c19_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c19_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c19_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c19_Singlets <- subset(Navin_hbca_c19, cells=rownames(Navin_hbca_c19@meta.data)[which(Navin_hbca_c19@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c19_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c19_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c19)
rm(Navin_hbca_c19_Singlets)
rm(Navin_hbca_c19_Quant)
rm(Navin_hbca_c19_Quant_Singlets)
rm(Navin_hbca_c19_Quant_Doublets)
rm(Navin_hbca_c19_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c19)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c19)
rm(sweep.stats.Navin_hbca_c19)
gc()



################################################################################################
########################                Navin_hbca_c20                  ########################
################################################################################################


# Load Data
Navin_hbca_c20 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c20.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c20 <- paramSweep(Navin_hbca_c20, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c20 <- summarizeSweep(sweep.res.Navin_hbca_c20, GT = FALSE)
bcmvn_Navin_hbca_c20 <- find.pK(sweep.stats.Navin_hbca_c20)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c20[which.max(bcmvn_Navin_hbca_c20$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c20@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c20@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c20 <- doubletFinder(Navin_hbca_c20, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c20_Quant <- (Navin_hbca_c20@meta.data$DF.classification == "Singlet")
Navin_hbca_c20_Quant_Singlets <- length(Navin_hbca_c20_Quant[Navin_hbca_c20_Quant== TRUE])
Navin_hbca_c20_Quant_Doublets <- length(Navin_hbca_c20_Quant[Navin_hbca_c20_Quant== FALSE])
Navin_hbca_c20_Quant_Doublets_Percent <- Navin_hbca_c20_Quant_Doublets / (Navin_hbca_c20_Quant_Doublets + Navin_hbca_c20_Quant_Singlets) * 100
Navin_hbca_c20_Quant <- as.data.frame(c(Navin_hbca_c20_Quant_Singlets, Navin_hbca_c20_Quant_Doublets, Navin_hbca_c20_Quant_Doublets_Percent))
colnames(Navin_hbca_c20_Quant) <- c("Navin_hbca_c20_Quant")
rownames(Navin_hbca_c20_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c20_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c20_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c20_Singlets <- subset(Navin_hbca_c20, cells=rownames(Navin_hbca_c20@meta.data)[which(Navin_hbca_c20@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c20_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c20_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c20)
rm(Navin_hbca_c20_Singlets)
rm(Navin_hbca_c20_Quant)
rm(Navin_hbca_c20_Quant_Singlets)
rm(Navin_hbca_c20_Quant_Doublets)
rm(Navin_hbca_c20_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c20)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c20)
rm(sweep.stats.Navin_hbca_c20)
gc()



################################################################################################
########################                Navin_hbca_c22                  ########################
################################################################################################


# Load Data
Navin_hbca_c22 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c22.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c22 <- paramSweep(Navin_hbca_c22, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c22 <- summarizeSweep(sweep.res.Navin_hbca_c22, GT = FALSE)
bcmvn_Navin_hbca_c22 <- find.pK(sweep.stats.Navin_hbca_c22)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c22[which.max(bcmvn_Navin_hbca_c22$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c22@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c22@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c22 <- doubletFinder(Navin_hbca_c22, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c22_Quant <- (Navin_hbca_c22@meta.data$DF.classification == "Singlet")
Navin_hbca_c22_Quant_Singlets <- length(Navin_hbca_c22_Quant[Navin_hbca_c22_Quant== TRUE])
Navin_hbca_c22_Quant_Doublets <- length(Navin_hbca_c22_Quant[Navin_hbca_c22_Quant== FALSE])
Navin_hbca_c22_Quant_Doublets_Percent <- Navin_hbca_c22_Quant_Doublets / (Navin_hbca_c22_Quant_Doublets + Navin_hbca_c22_Quant_Singlets) * 100
Navin_hbca_c22_Quant <- as.data.frame(c(Navin_hbca_c22_Quant_Singlets, Navin_hbca_c22_Quant_Doublets, Navin_hbca_c22_Quant_Doublets_Percent))
colnames(Navin_hbca_c22_Quant) <- c("Navin_hbca_c22_Quant")
rownames(Navin_hbca_c22_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c22_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c22_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c22_Singlets <- subset(Navin_hbca_c22, cells=rownames(Navin_hbca_c22@meta.data)[which(Navin_hbca_c22@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c22_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c22_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c22)
rm(Navin_hbca_c22_Singlets)
rm(Navin_hbca_c22_Quant)
rm(Navin_hbca_c22_Quant_Singlets)
rm(Navin_hbca_c22_Quant_Doublets)
rm(Navin_hbca_c22_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c22)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c22)
rm(sweep.stats.Navin_hbca_c22)
gc()



################################################################################################
########################                Navin_hbca_c23                  ########################
################################################################################################


# Load Data
Navin_hbca_c23 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c23.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c23 <- paramSweep(Navin_hbca_c23, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c23 <- summarizeSweep(sweep.res.Navin_hbca_c23, GT = FALSE)
bcmvn_Navin_hbca_c23 <- find.pK(sweep.stats.Navin_hbca_c23)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c23[which.max(bcmvn_Navin_hbca_c23$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c23@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c23@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c23 <- doubletFinder(Navin_hbca_c23, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c23_Quant <- (Navin_hbca_c23@meta.data$DF.classification == "Singlet")
Navin_hbca_c23_Quant_Singlets <- length(Navin_hbca_c23_Quant[Navin_hbca_c23_Quant== TRUE])
Navin_hbca_c23_Quant_Doublets <- length(Navin_hbca_c23_Quant[Navin_hbca_c23_Quant== FALSE])
Navin_hbca_c23_Quant_Doublets_Percent <- Navin_hbca_c23_Quant_Doublets / (Navin_hbca_c23_Quant_Doublets + Navin_hbca_c23_Quant_Singlets) * 100
Navin_hbca_c23_Quant <- as.data.frame(c(Navin_hbca_c23_Quant_Singlets, Navin_hbca_c23_Quant_Doublets, Navin_hbca_c23_Quant_Doublets_Percent))
colnames(Navin_hbca_c23_Quant) <- c("Navin_hbca_c23_Quant")
rownames(Navin_hbca_c23_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c23_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c23_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c23_Singlets <- subset(Navin_hbca_c23, cells=rownames(Navin_hbca_c23@meta.data)[which(Navin_hbca_c23@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c23_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c23_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c23)
rm(Navin_hbca_c23_Singlets)
rm(Navin_hbca_c23_Quant)
rm(Navin_hbca_c23_Quant_Singlets)
rm(Navin_hbca_c23_Quant_Doublets)
rm(Navin_hbca_c23_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c23)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c23)
rm(sweep.stats.Navin_hbca_c23)
gc()



################################################################################################
########################                Navin_hbca_c24                  ########################
################################################################################################


# Load Data
Navin_hbca_c24 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c24.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c24 <- paramSweep(Navin_hbca_c24, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c24 <- summarizeSweep(sweep.res.Navin_hbca_c24, GT = FALSE)
bcmvn_Navin_hbca_c24 <- find.pK(sweep.stats.Navin_hbca_c24)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c24[which.max(bcmvn_Navin_hbca_c24$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c24@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c24@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c24 <- doubletFinder(Navin_hbca_c24, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c24_Quant <- (Navin_hbca_c24@meta.data$DF.classification == "Singlet")
Navin_hbca_c24_Quant_Singlets <- length(Navin_hbca_c24_Quant[Navin_hbca_c24_Quant== TRUE])
Navin_hbca_c24_Quant_Doublets <- length(Navin_hbca_c24_Quant[Navin_hbca_c24_Quant== FALSE])
Navin_hbca_c24_Quant_Doublets_Percent <- Navin_hbca_c24_Quant_Doublets / (Navin_hbca_c24_Quant_Doublets + Navin_hbca_c24_Quant_Singlets) * 100
Navin_hbca_c24_Quant <- as.data.frame(c(Navin_hbca_c24_Quant_Singlets, Navin_hbca_c24_Quant_Doublets, Navin_hbca_c24_Quant_Doublets_Percent))
colnames(Navin_hbca_c24_Quant) <- c("Navin_hbca_c24_Quant")
rownames(Navin_hbca_c24_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c24_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c24_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c24_Singlets <- subset(Navin_hbca_c24, cells=rownames(Navin_hbca_c24@meta.data)[which(Navin_hbca_c24@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c24_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c24_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c24)
rm(Navin_hbca_c24_Singlets)
rm(Navin_hbca_c24_Quant)
rm(Navin_hbca_c24_Quant_Singlets)
rm(Navin_hbca_c24_Quant_Doublets)
rm(Navin_hbca_c24_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c24)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c24)
rm(sweep.stats.Navin_hbca_c24)
gc()



################################################################################################
########################                Navin_hbca_c25                  ########################
################################################################################################


# Load Data
Navin_hbca_c25 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c25.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c25 <- paramSweep(Navin_hbca_c25, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c25 <- summarizeSweep(sweep.res.Navin_hbca_c25, GT = FALSE)
bcmvn_Navin_hbca_c25 <- find.pK(sweep.stats.Navin_hbca_c25)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c25[which.max(bcmvn_Navin_hbca_c25$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c25@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c25@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c25 <- doubletFinder(Navin_hbca_c25, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c25_Quant <- (Navin_hbca_c25@meta.data$DF.classification == "Singlet")
Navin_hbca_c25_Quant_Singlets <- length(Navin_hbca_c25_Quant[Navin_hbca_c25_Quant== TRUE])
Navin_hbca_c25_Quant_Doublets <- length(Navin_hbca_c25_Quant[Navin_hbca_c25_Quant== FALSE])
Navin_hbca_c25_Quant_Doublets_Percent <- Navin_hbca_c25_Quant_Doublets / (Navin_hbca_c25_Quant_Doublets + Navin_hbca_c25_Quant_Singlets) * 100
Navin_hbca_c25_Quant <- as.data.frame(c(Navin_hbca_c25_Quant_Singlets, Navin_hbca_c25_Quant_Doublets, Navin_hbca_c25_Quant_Doublets_Percent))
colnames(Navin_hbca_c25_Quant) <- c("Navin_hbca_c25_Quant")
rownames(Navin_hbca_c25_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c25_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c25_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c25_Singlets <- subset(Navin_hbca_c25, cells=rownames(Navin_hbca_c25@meta.data)[which(Navin_hbca_c25@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c25_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c25_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c25)
rm(Navin_hbca_c25_Singlets)
rm(Navin_hbca_c25_Quant)
rm(Navin_hbca_c25_Quant_Singlets)
rm(Navin_hbca_c25_Quant_Doublets)
rm(Navin_hbca_c25_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c25)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c25)
rm(sweep.stats.Navin_hbca_c25)
gc()



################################################################################################
########################                Navin_hbca_c26                  ########################
################################################################################################


# Load Data
Navin_hbca_c26 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c26.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c26 <- paramSweep(Navin_hbca_c26, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c26 <- summarizeSweep(sweep.res.Navin_hbca_c26, GT = FALSE)
bcmvn_Navin_hbca_c26 <- find.pK(sweep.stats.Navin_hbca_c26)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c26[which.max(bcmvn_Navin_hbca_c26$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c26@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c26@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c26 <- doubletFinder(Navin_hbca_c26, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c26_Quant <- (Navin_hbca_c26@meta.data$DF.classification == "Singlet")
Navin_hbca_c26_Quant_Singlets <- length(Navin_hbca_c26_Quant[Navin_hbca_c26_Quant== TRUE])
Navin_hbca_c26_Quant_Doublets <- length(Navin_hbca_c26_Quant[Navin_hbca_c26_Quant== FALSE])
Navin_hbca_c26_Quant_Doublets_Percent <- Navin_hbca_c26_Quant_Doublets / (Navin_hbca_c26_Quant_Doublets + Navin_hbca_c26_Quant_Singlets) * 100
Navin_hbca_c26_Quant <- as.data.frame(c(Navin_hbca_c26_Quant_Singlets, Navin_hbca_c26_Quant_Doublets, Navin_hbca_c26_Quant_Doublets_Percent))
colnames(Navin_hbca_c26_Quant) <- c("Navin_hbca_c26_Quant")
rownames(Navin_hbca_c26_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c26_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c26_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c26_Singlets <- subset(Navin_hbca_c26, cells=rownames(Navin_hbca_c26@meta.data)[which(Navin_hbca_c26@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c26_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c26_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c26)
rm(Navin_hbca_c26_Singlets)
rm(Navin_hbca_c26_Quant)
rm(Navin_hbca_c26_Quant_Singlets)
rm(Navin_hbca_c26_Quant_Doublets)
rm(Navin_hbca_c26_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c26)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c26)
rm(sweep.stats.Navin_hbca_c26)
gc()



################################################################################################
########################                Navin_hbca_c31                  ########################
################################################################################################


# Load Data
Navin_hbca_c31 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c31.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c31 <- paramSweep(Navin_hbca_c31, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c31 <- summarizeSweep(sweep.res.Navin_hbca_c31, GT = FALSE)
bcmvn_Navin_hbca_c31 <- find.pK(sweep.stats.Navin_hbca_c31)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c31[which.max(bcmvn_Navin_hbca_c31$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c31@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c31@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c31 <- doubletFinder(Navin_hbca_c31, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c31_Quant <- (Navin_hbca_c31@meta.data$DF.classification == "Singlet")
Navin_hbca_c31_Quant_Singlets <- length(Navin_hbca_c31_Quant[Navin_hbca_c31_Quant== TRUE])
Navin_hbca_c31_Quant_Doublets <- length(Navin_hbca_c31_Quant[Navin_hbca_c31_Quant== FALSE])
Navin_hbca_c31_Quant_Doublets_Percent <- Navin_hbca_c31_Quant_Doublets / (Navin_hbca_c31_Quant_Doublets + Navin_hbca_c31_Quant_Singlets) * 100
Navin_hbca_c31_Quant <- as.data.frame(c(Navin_hbca_c31_Quant_Singlets, Navin_hbca_c31_Quant_Doublets, Navin_hbca_c31_Quant_Doublets_Percent))
colnames(Navin_hbca_c31_Quant) <- c("Navin_hbca_c31_Quant")
rownames(Navin_hbca_c31_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c31_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c31_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c31_Singlets <- subset(Navin_hbca_c31, cells=rownames(Navin_hbca_c31@meta.data)[which(Navin_hbca_c31@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c31_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c31_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c31)
rm(Navin_hbca_c31_Singlets)
rm(Navin_hbca_c31_Quant)
rm(Navin_hbca_c31_Quant_Singlets)
rm(Navin_hbca_c31_Quant_Doublets)
rm(Navin_hbca_c31_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c31)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c31)
rm(sweep.stats.Navin_hbca_c31)
gc()



################################################################################################
########################                Navin_hbca_c32                  ########################
################################################################################################


# Load Data
Navin_hbca_c32 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c32.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c32 <- paramSweep(Navin_hbca_c32, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c32 <- summarizeSweep(sweep.res.Navin_hbca_c32, GT = FALSE)
bcmvn_Navin_hbca_c32 <- find.pK(sweep.stats.Navin_hbca_c32)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c32[which.max(bcmvn_Navin_hbca_c32$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c32@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c32@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c32 <- doubletFinder(Navin_hbca_c32, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c32_Quant <- (Navin_hbca_c32@meta.data$DF.classification == "Singlet")
Navin_hbca_c32_Quant_Singlets <- length(Navin_hbca_c32_Quant[Navin_hbca_c32_Quant== TRUE])
Navin_hbca_c32_Quant_Doublets <- length(Navin_hbca_c32_Quant[Navin_hbca_c32_Quant== FALSE])
Navin_hbca_c32_Quant_Doublets_Percent <- Navin_hbca_c32_Quant_Doublets / (Navin_hbca_c32_Quant_Doublets + Navin_hbca_c32_Quant_Singlets) * 100
Navin_hbca_c32_Quant <- as.data.frame(c(Navin_hbca_c32_Quant_Singlets, Navin_hbca_c32_Quant_Doublets, Navin_hbca_c32_Quant_Doublets_Percent))
colnames(Navin_hbca_c32_Quant) <- c("Navin_hbca_c32_Quant")
rownames(Navin_hbca_c32_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c32_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c32_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c32_Singlets <- subset(Navin_hbca_c32, cells=rownames(Navin_hbca_c32@meta.data)[which(Navin_hbca_c32@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c32_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c32_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c32)
rm(Navin_hbca_c32_Singlets)
rm(Navin_hbca_c32_Quant)
rm(Navin_hbca_c32_Quant_Singlets)
rm(Navin_hbca_c32_Quant_Doublets)
rm(Navin_hbca_c32_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c32)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c32)
rm(sweep.stats.Navin_hbca_c32)
gc()



################################################################################################
########################                Navin_hbca_c50                  ########################
################################################################################################


# Load Data
Navin_hbca_c50 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c50.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c50 <- paramSweep(Navin_hbca_c50, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c50 <- summarizeSweep(sweep.res.Navin_hbca_c50, GT = FALSE)
bcmvn_Navin_hbca_c50 <- find.pK(sweep.stats.Navin_hbca_c50)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c50[which.max(bcmvn_Navin_hbca_c50$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c50@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c50@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c50 <- doubletFinder(Navin_hbca_c50, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c50_Quant <- (Navin_hbca_c50@meta.data$DF.classification == "Singlet")
Navin_hbca_c50_Quant_Singlets <- length(Navin_hbca_c50_Quant[Navin_hbca_c50_Quant== TRUE])
Navin_hbca_c50_Quant_Doublets <- length(Navin_hbca_c50_Quant[Navin_hbca_c50_Quant== FALSE])
Navin_hbca_c50_Quant_Doublets_Percent <- Navin_hbca_c50_Quant_Doublets / (Navin_hbca_c50_Quant_Doublets + Navin_hbca_c50_Quant_Singlets) * 100
Navin_hbca_c50_Quant <- as.data.frame(c(Navin_hbca_c50_Quant_Singlets, Navin_hbca_c50_Quant_Doublets, Navin_hbca_c50_Quant_Doublets_Percent))
colnames(Navin_hbca_c50_Quant) <- c("Navin_hbca_c50_Quant")
rownames(Navin_hbca_c50_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c50_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c50_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c50_Singlets <- subset(Navin_hbca_c50, cells=rownames(Navin_hbca_c50@meta.data)[which(Navin_hbca_c50@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c50_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c50_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c50)
rm(Navin_hbca_c50_Singlets)
rm(Navin_hbca_c50_Quant)
rm(Navin_hbca_c50_Quant_Singlets)
rm(Navin_hbca_c50_Quant_Doublets)
rm(Navin_hbca_c50_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c50)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c50)
rm(sweep.stats.Navin_hbca_c50)
gc()



################################################################################################
########################                Navin_hbca_c51                  ########################
################################################################################################


# Load Data
Navin_hbca_c51 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c51.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c51 <- paramSweep(Navin_hbca_c51, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c51 <- summarizeSweep(sweep.res.Navin_hbca_c51, GT = FALSE)
bcmvn_Navin_hbca_c51 <- find.pK(sweep.stats.Navin_hbca_c51)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c51[which.max(bcmvn_Navin_hbca_c51$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c51@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c51@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c51 <- doubletFinder(Navin_hbca_c51, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c51_Quant <- (Navin_hbca_c51@meta.data$DF.classification == "Singlet")
Navin_hbca_c51_Quant_Singlets <- length(Navin_hbca_c51_Quant[Navin_hbca_c51_Quant== TRUE])
Navin_hbca_c51_Quant_Doublets <- length(Navin_hbca_c51_Quant[Navin_hbca_c51_Quant== FALSE])
Navin_hbca_c51_Quant_Doublets_Percent <- Navin_hbca_c51_Quant_Doublets / (Navin_hbca_c51_Quant_Doublets + Navin_hbca_c51_Quant_Singlets) * 100
Navin_hbca_c51_Quant <- as.data.frame(c(Navin_hbca_c51_Quant_Singlets, Navin_hbca_c51_Quant_Doublets, Navin_hbca_c51_Quant_Doublets_Percent))
colnames(Navin_hbca_c51_Quant) <- c("Navin_hbca_c51_Quant")
rownames(Navin_hbca_c51_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c51_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c51_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c51_Singlets <- subset(Navin_hbca_c51, cells=rownames(Navin_hbca_c51@meta.data)[which(Navin_hbca_c51@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c51_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c51_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c51)
rm(Navin_hbca_c51_Singlets)
rm(Navin_hbca_c51_Quant)
rm(Navin_hbca_c51_Quant_Singlets)
rm(Navin_hbca_c51_Quant_Doublets)
rm(Navin_hbca_c51_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c51)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c51)
rm(sweep.stats.Navin_hbca_c51)
gc()



################################################################################################
########################                Navin_hbca_c52                  ########################
################################################################################################


# Load Data
Navin_hbca_c52 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c52.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c52 <- paramSweep(Navin_hbca_c52, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c52 <- summarizeSweep(sweep.res.Navin_hbca_c52, GT = FALSE)
bcmvn_Navin_hbca_c52 <- find.pK(sweep.stats.Navin_hbca_c52)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c52[which.max(bcmvn_Navin_hbca_c52$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c52@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c52@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c52 <- doubletFinder(Navin_hbca_c52, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c52_Quant <- (Navin_hbca_c52@meta.data$DF.classification == "Singlet")
Navin_hbca_c52_Quant_Singlets <- length(Navin_hbca_c52_Quant[Navin_hbca_c52_Quant== TRUE])
Navin_hbca_c52_Quant_Doublets <- length(Navin_hbca_c52_Quant[Navin_hbca_c52_Quant== FALSE])
Navin_hbca_c52_Quant_Doublets_Percent <- Navin_hbca_c52_Quant_Doublets / (Navin_hbca_c52_Quant_Doublets + Navin_hbca_c52_Quant_Singlets) * 100
Navin_hbca_c52_Quant <- as.data.frame(c(Navin_hbca_c52_Quant_Singlets, Navin_hbca_c52_Quant_Doublets, Navin_hbca_c52_Quant_Doublets_Percent))
colnames(Navin_hbca_c52_Quant) <- c("Navin_hbca_c52_Quant")
rownames(Navin_hbca_c52_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c52_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c52_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c52_Singlets <- subset(Navin_hbca_c52, cells=rownames(Navin_hbca_c52@meta.data)[which(Navin_hbca_c52@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c52_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c52_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c52)
rm(Navin_hbca_c52_Singlets)
rm(Navin_hbca_c52_Quant)
rm(Navin_hbca_c52_Quant_Singlets)
rm(Navin_hbca_c52_Quant_Doublets)
rm(Navin_hbca_c52_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c52)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c52)
rm(sweep.stats.Navin_hbca_c52)
gc()



################################################################################################
########################                Navin_hbca_c53                  ########################
################################################################################################


# Load Data
Navin_hbca_c53 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c53.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c53 <- paramSweep(Navin_hbca_c53, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c53 <- summarizeSweep(sweep.res.Navin_hbca_c53, GT = FALSE)
bcmvn_Navin_hbca_c53 <- find.pK(sweep.stats.Navin_hbca_c53)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c53[which.max(bcmvn_Navin_hbca_c53$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c53@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c53@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c53 <- doubletFinder(Navin_hbca_c53, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c53_Quant <- (Navin_hbca_c53@meta.data$DF.classification == "Singlet")
Navin_hbca_c53_Quant_Singlets <- length(Navin_hbca_c53_Quant[Navin_hbca_c53_Quant== TRUE])
Navin_hbca_c53_Quant_Doublets <- length(Navin_hbca_c53_Quant[Navin_hbca_c53_Quant== FALSE])
Navin_hbca_c53_Quant_Doublets_Percent <- Navin_hbca_c53_Quant_Doublets / (Navin_hbca_c53_Quant_Doublets + Navin_hbca_c53_Quant_Singlets) * 100
Navin_hbca_c53_Quant <- as.data.frame(c(Navin_hbca_c53_Quant_Singlets, Navin_hbca_c53_Quant_Doublets, Navin_hbca_c53_Quant_Doublets_Percent))
colnames(Navin_hbca_c53_Quant) <- c("Navin_hbca_c53_Quant")
rownames(Navin_hbca_c53_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c53_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c53_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c53_Singlets <- subset(Navin_hbca_c53, cells=rownames(Navin_hbca_c53@meta.data)[which(Navin_hbca_c53@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c53_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c53_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c53)
rm(Navin_hbca_c53_Singlets)
rm(Navin_hbca_c53_Quant)
rm(Navin_hbca_c53_Quant_Singlets)
rm(Navin_hbca_c53_Quant_Doublets)
rm(Navin_hbca_c53_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c53)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c53)
rm(sweep.stats.Navin_hbca_c53)
gc()



################################################################################################
########################                Navin_hbca_c54                  ########################
################################################################################################


# Load Data
Navin_hbca_c54 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c54.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c54 <- paramSweep(Navin_hbca_c54, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c54 <- summarizeSweep(sweep.res.Navin_hbca_c54, GT = FALSE)
bcmvn_Navin_hbca_c54 <- find.pK(sweep.stats.Navin_hbca_c54)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c54[which.max(bcmvn_Navin_hbca_c54$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c54@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c54@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c54 <- doubletFinder(Navin_hbca_c54, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c54_Quant <- (Navin_hbca_c54@meta.data$DF.classification == "Singlet")
Navin_hbca_c54_Quant_Singlets <- length(Navin_hbca_c54_Quant[Navin_hbca_c54_Quant== TRUE])
Navin_hbca_c54_Quant_Doublets <- length(Navin_hbca_c54_Quant[Navin_hbca_c54_Quant== FALSE])
Navin_hbca_c54_Quant_Doublets_Percent <- Navin_hbca_c54_Quant_Doublets / (Navin_hbca_c54_Quant_Doublets + Navin_hbca_c54_Quant_Singlets) * 100
Navin_hbca_c54_Quant <- as.data.frame(c(Navin_hbca_c54_Quant_Singlets, Navin_hbca_c54_Quant_Doublets, Navin_hbca_c54_Quant_Doublets_Percent))
colnames(Navin_hbca_c54_Quant) <- c("Navin_hbca_c54_Quant")
rownames(Navin_hbca_c54_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c54_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c54_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c54_Singlets <- subset(Navin_hbca_c54, cells=rownames(Navin_hbca_c54@meta.data)[which(Navin_hbca_c54@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c54_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c54_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c54)
rm(Navin_hbca_c54_Singlets)
rm(Navin_hbca_c54_Quant)
rm(Navin_hbca_c54_Quant_Singlets)
rm(Navin_hbca_c54_Quant_Doublets)
rm(Navin_hbca_c54_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c54)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c54)
rm(sweep.stats.Navin_hbca_c54)
gc()



################################################################################################
########################                Navin_hbca_c55                  ########################
################################################################################################


# Load Data
Navin_hbca_c55 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c55.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c55 <- paramSweep(Navin_hbca_c55, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c55 <- summarizeSweep(sweep.res.Navin_hbca_c55, GT = FALSE)
bcmvn_Navin_hbca_c55 <- find.pK(sweep.stats.Navin_hbca_c55)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c55[which.max(bcmvn_Navin_hbca_c55$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c55@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c55@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c55 <- doubletFinder(Navin_hbca_c55, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c55_Quant <- (Navin_hbca_c55@meta.data$DF.classification == "Singlet")
Navin_hbca_c55_Quant_Singlets <- length(Navin_hbca_c55_Quant[Navin_hbca_c55_Quant== TRUE])
Navin_hbca_c55_Quant_Doublets <- length(Navin_hbca_c55_Quant[Navin_hbca_c55_Quant== FALSE])
Navin_hbca_c55_Quant_Doublets_Percent <- Navin_hbca_c55_Quant_Doublets / (Navin_hbca_c55_Quant_Doublets + Navin_hbca_c55_Quant_Singlets) * 100
Navin_hbca_c55_Quant <- as.data.frame(c(Navin_hbca_c55_Quant_Singlets, Navin_hbca_c55_Quant_Doublets, Navin_hbca_c55_Quant_Doublets_Percent))
colnames(Navin_hbca_c55_Quant) <- c("Navin_hbca_c55_Quant")
rownames(Navin_hbca_c55_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c55_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c55_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c55_Singlets <- subset(Navin_hbca_c55, cells=rownames(Navin_hbca_c55@meta.data)[which(Navin_hbca_c55@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c55_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c55_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c55)
rm(Navin_hbca_c55_Singlets)
rm(Navin_hbca_c55_Quant)
rm(Navin_hbca_c55_Quant_Singlets)
rm(Navin_hbca_c55_Quant_Doublets)
rm(Navin_hbca_c55_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c55)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c55)
rm(sweep.stats.Navin_hbca_c55)
gc()



################################################################################################
########################                Navin_hbca_c56                  ########################
################################################################################################


# Load Data
Navin_hbca_c56 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c56.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c56 <- paramSweep(Navin_hbca_c56, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c56 <- summarizeSweep(sweep.res.Navin_hbca_c56, GT = FALSE)
bcmvn_Navin_hbca_c56 <- find.pK(sweep.stats.Navin_hbca_c56)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c56[which.max(bcmvn_Navin_hbca_c56$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c56@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c56@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c56 <- doubletFinder(Navin_hbca_c56, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c56_Quant <- (Navin_hbca_c56@meta.data$DF.classification == "Singlet")
Navin_hbca_c56_Quant_Singlets <- length(Navin_hbca_c56_Quant[Navin_hbca_c56_Quant== TRUE])
Navin_hbca_c56_Quant_Doublets <- length(Navin_hbca_c56_Quant[Navin_hbca_c56_Quant== FALSE])
Navin_hbca_c56_Quant_Doublets_Percent <- Navin_hbca_c56_Quant_Doublets / (Navin_hbca_c56_Quant_Doublets + Navin_hbca_c56_Quant_Singlets) * 100
Navin_hbca_c56_Quant <- as.data.frame(c(Navin_hbca_c56_Quant_Singlets, Navin_hbca_c56_Quant_Doublets, Navin_hbca_c56_Quant_Doublets_Percent))
colnames(Navin_hbca_c56_Quant) <- c("Navin_hbca_c56_Quant")
rownames(Navin_hbca_c56_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c56_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c56_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c56_Singlets <- subset(Navin_hbca_c56, cells=rownames(Navin_hbca_c56@meta.data)[which(Navin_hbca_c56@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c56_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c56_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c56)
rm(Navin_hbca_c56_Singlets)
rm(Navin_hbca_c56_Quant)
rm(Navin_hbca_c56_Quant_Singlets)
rm(Navin_hbca_c56_Quant_Doublets)
rm(Navin_hbca_c56_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c56)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c56)
rm(sweep.stats.Navin_hbca_c56)
gc()



################################################################################################
########################                Navin_hbca_c57                  ########################
################################################################################################


# Load Data
Navin_hbca_c57 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c57.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c57 <- paramSweep(Navin_hbca_c57, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c57 <- summarizeSweep(sweep.res.Navin_hbca_c57, GT = FALSE)
bcmvn_Navin_hbca_c57 <- find.pK(sweep.stats.Navin_hbca_c57)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c57[which.max(bcmvn_Navin_hbca_c57$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c57@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c57@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c57 <- doubletFinder(Navin_hbca_c57, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c57_Quant <- (Navin_hbca_c57@meta.data$DF.classification == "Singlet")
Navin_hbca_c57_Quant_Singlets <- length(Navin_hbca_c57_Quant[Navin_hbca_c57_Quant== TRUE])
Navin_hbca_c57_Quant_Doublets <- length(Navin_hbca_c57_Quant[Navin_hbca_c57_Quant== FALSE])
Navin_hbca_c57_Quant_Doublets_Percent <- Navin_hbca_c57_Quant_Doublets / (Navin_hbca_c57_Quant_Doublets + Navin_hbca_c57_Quant_Singlets) * 100
Navin_hbca_c57_Quant <- as.data.frame(c(Navin_hbca_c57_Quant_Singlets, Navin_hbca_c57_Quant_Doublets, Navin_hbca_c57_Quant_Doublets_Percent))
colnames(Navin_hbca_c57_Quant) <- c("Navin_hbca_c57_Quant")
rownames(Navin_hbca_c57_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c57_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c57_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c57_Singlets <- subset(Navin_hbca_c57, cells=rownames(Navin_hbca_c57@meta.data)[which(Navin_hbca_c57@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c57_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c57_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c57)
rm(Navin_hbca_c57_Singlets)
rm(Navin_hbca_c57_Quant)
rm(Navin_hbca_c57_Quant_Singlets)
rm(Navin_hbca_c57_Quant_Doublets)
rm(Navin_hbca_c57_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c57)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c57)
rm(sweep.stats.Navin_hbca_c57)
gc()



################################################################################################
########################                Navin_hbca_c58                  ########################
################################################################################################


# Load Data
Navin_hbca_c58 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c58.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c58 <- paramSweep(Navin_hbca_c58, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c58 <- summarizeSweep(sweep.res.Navin_hbca_c58, GT = FALSE)
bcmvn_Navin_hbca_c58 <- find.pK(sweep.stats.Navin_hbca_c58)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c58[which.max(bcmvn_Navin_hbca_c58$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c58@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c58@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c58 <- doubletFinder(Navin_hbca_c58, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c58_Quant <- (Navin_hbca_c58@meta.data$DF.classification == "Singlet")
Navin_hbca_c58_Quant_Singlets <- length(Navin_hbca_c58_Quant[Navin_hbca_c58_Quant== TRUE])
Navin_hbca_c58_Quant_Doublets <- length(Navin_hbca_c58_Quant[Navin_hbca_c58_Quant== FALSE])
Navin_hbca_c58_Quant_Doublets_Percent <- Navin_hbca_c58_Quant_Doublets / (Navin_hbca_c58_Quant_Doublets + Navin_hbca_c58_Quant_Singlets) * 100
Navin_hbca_c58_Quant <- as.data.frame(c(Navin_hbca_c58_Quant_Singlets, Navin_hbca_c58_Quant_Doublets, Navin_hbca_c58_Quant_Doublets_Percent))
colnames(Navin_hbca_c58_Quant) <- c("Navin_hbca_c58_Quant")
rownames(Navin_hbca_c58_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c58_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c58_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c58_Singlets <- subset(Navin_hbca_c58, cells=rownames(Navin_hbca_c58@meta.data)[which(Navin_hbca_c58@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c58_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c58_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c58)
rm(Navin_hbca_c58_Singlets)
rm(Navin_hbca_c58_Quant)
rm(Navin_hbca_c58_Quant_Singlets)
rm(Navin_hbca_c58_Quant_Doublets)
rm(Navin_hbca_c58_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c58)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c58)
rm(sweep.stats.Navin_hbca_c58)
gc()



################################################################################################
########################                Navin_hbca_c59                  ########################
################################################################################################


# Load Data
Navin_hbca_c59 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c59.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c59 <- paramSweep(Navin_hbca_c59, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c59 <- summarizeSweep(sweep.res.Navin_hbca_c59, GT = FALSE)
bcmvn_Navin_hbca_c59 <- find.pK(sweep.stats.Navin_hbca_c59)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c59[which.max(bcmvn_Navin_hbca_c59$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c59@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c59@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c59 <- doubletFinder(Navin_hbca_c59, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c59_Quant <- (Navin_hbca_c59@meta.data$DF.classification == "Singlet")
Navin_hbca_c59_Quant_Singlets <- length(Navin_hbca_c59_Quant[Navin_hbca_c59_Quant== TRUE])
Navin_hbca_c59_Quant_Doublets <- length(Navin_hbca_c59_Quant[Navin_hbca_c59_Quant== FALSE])
Navin_hbca_c59_Quant_Doublets_Percent <- Navin_hbca_c59_Quant_Doublets / (Navin_hbca_c59_Quant_Doublets + Navin_hbca_c59_Quant_Singlets) * 100
Navin_hbca_c59_Quant <- as.data.frame(c(Navin_hbca_c59_Quant_Singlets, Navin_hbca_c59_Quant_Doublets, Navin_hbca_c59_Quant_Doublets_Percent))
colnames(Navin_hbca_c59_Quant) <- c("Navin_hbca_c59_Quant")
rownames(Navin_hbca_c59_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c59_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c59_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c59_Singlets <- subset(Navin_hbca_c59, cells=rownames(Navin_hbca_c59@meta.data)[which(Navin_hbca_c59@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c59_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c59_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c59)
rm(Navin_hbca_c59_Singlets)
rm(Navin_hbca_c59_Quant)
rm(Navin_hbca_c59_Quant_Singlets)
rm(Navin_hbca_c59_Quant_Doublets)
rm(Navin_hbca_c59_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c59)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c59)
rm(sweep.stats.Navin_hbca_c59)
gc()



################################################################################################
########################                Navin_hbca_c60                  ########################
################################################################################################


# Load Data
Navin_hbca_c60 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c60.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c60 <- paramSweep(Navin_hbca_c60, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c60 <- summarizeSweep(sweep.res.Navin_hbca_c60, GT = FALSE)
bcmvn_Navin_hbca_c60 <- find.pK(sweep.stats.Navin_hbca_c60)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c60[which.max(bcmvn_Navin_hbca_c60$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c60@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c60@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c60 <- doubletFinder(Navin_hbca_c60, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c60_Quant <- (Navin_hbca_c60@meta.data$DF.classification == "Singlet")
Navin_hbca_c60_Quant_Singlets <- length(Navin_hbca_c60_Quant[Navin_hbca_c60_Quant== TRUE])
Navin_hbca_c60_Quant_Doublets <- length(Navin_hbca_c60_Quant[Navin_hbca_c60_Quant== FALSE])
Navin_hbca_c60_Quant_Doublets_Percent <- Navin_hbca_c60_Quant_Doublets / (Navin_hbca_c60_Quant_Doublets + Navin_hbca_c60_Quant_Singlets) * 100
Navin_hbca_c60_Quant <- as.data.frame(c(Navin_hbca_c60_Quant_Singlets, Navin_hbca_c60_Quant_Doublets, Navin_hbca_c60_Quant_Doublets_Percent))
colnames(Navin_hbca_c60_Quant) <- c("Navin_hbca_c60_Quant")
rownames(Navin_hbca_c60_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c60_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c60_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c60_Singlets <- subset(Navin_hbca_c60, cells=rownames(Navin_hbca_c60@meta.data)[which(Navin_hbca_c60@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c60_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c60_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c60)
rm(Navin_hbca_c60_Singlets)
rm(Navin_hbca_c60_Quant)
rm(Navin_hbca_c60_Quant_Singlets)
rm(Navin_hbca_c60_Quant_Doublets)
rm(Navin_hbca_c60_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c60)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c60)
rm(sweep.stats.Navin_hbca_c60)
gc()



################################################################################################
########################                Navin_hbca_c61                  ########################
################################################################################################


# Load Data
Navin_hbca_c61 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c61.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c61 <- paramSweep(Navin_hbca_c61, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c61 <- summarizeSweep(sweep.res.Navin_hbca_c61, GT = FALSE)
bcmvn_Navin_hbca_c61 <- find.pK(sweep.stats.Navin_hbca_c61)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c61[which.max(bcmvn_Navin_hbca_c61$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c61@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c61@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c61 <- doubletFinder(Navin_hbca_c61, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c61_Quant <- (Navin_hbca_c61@meta.data$DF.classification == "Singlet")
Navin_hbca_c61_Quant_Singlets <- length(Navin_hbca_c61_Quant[Navin_hbca_c61_Quant== TRUE])
Navin_hbca_c61_Quant_Doublets <- length(Navin_hbca_c61_Quant[Navin_hbca_c61_Quant== FALSE])
Navin_hbca_c61_Quant_Doublets_Percent <- Navin_hbca_c61_Quant_Doublets / (Navin_hbca_c61_Quant_Doublets + Navin_hbca_c61_Quant_Singlets) * 100
Navin_hbca_c61_Quant <- as.data.frame(c(Navin_hbca_c61_Quant_Singlets, Navin_hbca_c61_Quant_Doublets, Navin_hbca_c61_Quant_Doublets_Percent))
colnames(Navin_hbca_c61_Quant) <- c("Navin_hbca_c61_Quant")
rownames(Navin_hbca_c61_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c61_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c61_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c61_Singlets <- subset(Navin_hbca_c61, cells=rownames(Navin_hbca_c61@meta.data)[which(Navin_hbca_c61@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c61_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c61_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c61)
rm(Navin_hbca_c61_Singlets)
rm(Navin_hbca_c61_Quant)
rm(Navin_hbca_c61_Quant_Singlets)
rm(Navin_hbca_c61_Quant_Doublets)
rm(Navin_hbca_c61_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c61)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c61)
rm(sweep.stats.Navin_hbca_c61)
gc()



################################################################################################
########################                Navin_hbca_c62                  ########################
################################################################################################


# Load Data
Navin_hbca_c62 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c62.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c62 <- paramSweep(Navin_hbca_c62, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c62 <- summarizeSweep(sweep.res.Navin_hbca_c62, GT = FALSE)
bcmvn_Navin_hbca_c62 <- find.pK(sweep.stats.Navin_hbca_c62)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c62[which.max(bcmvn_Navin_hbca_c62$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c62@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c62@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c62 <- doubletFinder(Navin_hbca_c62, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c62_Quant <- (Navin_hbca_c62@meta.data$DF.classification == "Singlet")
Navin_hbca_c62_Quant_Singlets <- length(Navin_hbca_c62_Quant[Navin_hbca_c62_Quant== TRUE])
Navin_hbca_c62_Quant_Doublets <- length(Navin_hbca_c62_Quant[Navin_hbca_c62_Quant== FALSE])
Navin_hbca_c62_Quant_Doublets_Percent <- Navin_hbca_c62_Quant_Doublets / (Navin_hbca_c62_Quant_Doublets + Navin_hbca_c62_Quant_Singlets) * 100
Navin_hbca_c62_Quant <- as.data.frame(c(Navin_hbca_c62_Quant_Singlets, Navin_hbca_c62_Quant_Doublets, Navin_hbca_c62_Quant_Doublets_Percent))
colnames(Navin_hbca_c62_Quant) <- c("Navin_hbca_c62_Quant")
rownames(Navin_hbca_c62_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c62_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c62_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c62_Singlets <- subset(Navin_hbca_c62, cells=rownames(Navin_hbca_c62@meta.data)[which(Navin_hbca_c62@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c62_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c62_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c62)
rm(Navin_hbca_c62_Singlets)
rm(Navin_hbca_c62_Quant)
rm(Navin_hbca_c62_Quant_Singlets)
rm(Navin_hbca_c62_Quant_Doublets)
rm(Navin_hbca_c62_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c62)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c62)
rm(sweep.stats.Navin_hbca_c62)
gc()



################################################################################################
########################                Navin_hbca_c63                  ########################
################################################################################################


# Load Data
Navin_hbca_c63 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c63.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c63 <- paramSweep(Navin_hbca_c63, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c63 <- summarizeSweep(sweep.res.Navin_hbca_c63, GT = FALSE)
bcmvn_Navin_hbca_c63 <- find.pK(sweep.stats.Navin_hbca_c63)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c63[which.max(bcmvn_Navin_hbca_c63$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c63@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c63@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c63 <- doubletFinder(Navin_hbca_c63, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c63_Quant <- (Navin_hbca_c63@meta.data$DF.classification == "Singlet")
Navin_hbca_c63_Quant_Singlets <- length(Navin_hbca_c63_Quant[Navin_hbca_c63_Quant== TRUE])
Navin_hbca_c63_Quant_Doublets <- length(Navin_hbca_c63_Quant[Navin_hbca_c63_Quant== FALSE])
Navin_hbca_c63_Quant_Doublets_Percent <- Navin_hbca_c63_Quant_Doublets / (Navin_hbca_c63_Quant_Doublets + Navin_hbca_c63_Quant_Singlets) * 100
Navin_hbca_c63_Quant <- as.data.frame(c(Navin_hbca_c63_Quant_Singlets, Navin_hbca_c63_Quant_Doublets, Navin_hbca_c63_Quant_Doublets_Percent))
colnames(Navin_hbca_c63_Quant) <- c("Navin_hbca_c63_Quant")
rownames(Navin_hbca_c63_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c63_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c63_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c63_Singlets <- subset(Navin_hbca_c63, cells=rownames(Navin_hbca_c63@meta.data)[which(Navin_hbca_c63@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c63_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c63_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c63)
rm(Navin_hbca_c63_Singlets)
rm(Navin_hbca_c63_Quant)
rm(Navin_hbca_c63_Quant_Singlets)
rm(Navin_hbca_c63_Quant_Doublets)
rm(Navin_hbca_c63_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c63)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c63)
rm(sweep.stats.Navin_hbca_c63)
gc()



################################################################################################
########################                Navin_hbca_c64                  ########################
################################################################################################


# Load Data
Navin_hbca_c64 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c64.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c64 <- paramSweep(Navin_hbca_c64, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c64 <- summarizeSweep(sweep.res.Navin_hbca_c64, GT = FALSE)
bcmvn_Navin_hbca_c64 <- find.pK(sweep.stats.Navin_hbca_c64)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c64[which.max(bcmvn_Navin_hbca_c64$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c64@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c64@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c64 <- doubletFinder(Navin_hbca_c64, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c64_Quant <- (Navin_hbca_c64@meta.data$DF.classification == "Singlet")
Navin_hbca_c64_Quant_Singlets <- length(Navin_hbca_c64_Quant[Navin_hbca_c64_Quant== TRUE])
Navin_hbca_c64_Quant_Doublets <- length(Navin_hbca_c64_Quant[Navin_hbca_c64_Quant== FALSE])
Navin_hbca_c64_Quant_Doublets_Percent <- Navin_hbca_c64_Quant_Doublets / (Navin_hbca_c64_Quant_Doublets + Navin_hbca_c64_Quant_Singlets) * 100
Navin_hbca_c64_Quant <- as.data.frame(c(Navin_hbca_c64_Quant_Singlets, Navin_hbca_c64_Quant_Doublets, Navin_hbca_c64_Quant_Doublets_Percent))
colnames(Navin_hbca_c64_Quant) <- c("Navin_hbca_c64_Quant")
rownames(Navin_hbca_c64_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c64_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c64_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c64_Singlets <- subset(Navin_hbca_c64, cells=rownames(Navin_hbca_c64@meta.data)[which(Navin_hbca_c64@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c64_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c64_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c64)
rm(Navin_hbca_c64_Singlets)
rm(Navin_hbca_c64_Quant)
rm(Navin_hbca_c64_Quant_Singlets)
rm(Navin_hbca_c64_Quant_Doublets)
rm(Navin_hbca_c64_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c64)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c64)
rm(sweep.stats.Navin_hbca_c64)
gc()



################################################################################################
########################                Navin_hbca_c65                  ########################
################################################################################################


# Load Data
Navin_hbca_c65 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c65.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c65 <- paramSweep(Navin_hbca_c65, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c65 <- summarizeSweep(sweep.res.Navin_hbca_c65, GT = FALSE)
bcmvn_Navin_hbca_c65 <- find.pK(sweep.stats.Navin_hbca_c65)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c65[which.max(bcmvn_Navin_hbca_c65$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c65@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c65@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c65 <- doubletFinder(Navin_hbca_c65, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c65_Quant <- (Navin_hbca_c65@meta.data$DF.classification == "Singlet")
Navin_hbca_c65_Quant_Singlets <- length(Navin_hbca_c65_Quant[Navin_hbca_c65_Quant== TRUE])
Navin_hbca_c65_Quant_Doublets <- length(Navin_hbca_c65_Quant[Navin_hbca_c65_Quant== FALSE])
Navin_hbca_c65_Quant_Doublets_Percent <- Navin_hbca_c65_Quant_Doublets / (Navin_hbca_c65_Quant_Doublets + Navin_hbca_c65_Quant_Singlets) * 100
Navin_hbca_c65_Quant <- as.data.frame(c(Navin_hbca_c65_Quant_Singlets, Navin_hbca_c65_Quant_Doublets, Navin_hbca_c65_Quant_Doublets_Percent))
colnames(Navin_hbca_c65_Quant) <- c("Navin_hbca_c65_Quant")
rownames(Navin_hbca_c65_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c65_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c65_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c65_Singlets <- subset(Navin_hbca_c65, cells=rownames(Navin_hbca_c65@meta.data)[which(Navin_hbca_c65@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c65_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c65_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c65)
rm(Navin_hbca_c65_Singlets)
rm(Navin_hbca_c65_Quant)
rm(Navin_hbca_c65_Quant_Singlets)
rm(Navin_hbca_c65_Quant_Doublets)
rm(Navin_hbca_c65_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c65)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c65)
rm(sweep.stats.Navin_hbca_c65)
gc()



################################################################################################
########################                Navin_hbca_c66                  ########################
################################################################################################


# Load Data
Navin_hbca_c66 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c66.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c66 <- paramSweep(Navin_hbca_c66, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c66 <- summarizeSweep(sweep.res.Navin_hbca_c66, GT = FALSE)
bcmvn_Navin_hbca_c66 <- find.pK(sweep.stats.Navin_hbca_c66)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c66[which.max(bcmvn_Navin_hbca_c66$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c66@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c66@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c66 <- doubletFinder(Navin_hbca_c66, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c66_Quant <- (Navin_hbca_c66@meta.data$DF.classification == "Singlet")
Navin_hbca_c66_Quant_Singlets <- length(Navin_hbca_c66_Quant[Navin_hbca_c66_Quant== TRUE])
Navin_hbca_c66_Quant_Doublets <- length(Navin_hbca_c66_Quant[Navin_hbca_c66_Quant== FALSE])
Navin_hbca_c66_Quant_Doublets_Percent <- Navin_hbca_c66_Quant_Doublets / (Navin_hbca_c66_Quant_Doublets + Navin_hbca_c66_Quant_Singlets) * 100
Navin_hbca_c66_Quant <- as.data.frame(c(Navin_hbca_c66_Quant_Singlets, Navin_hbca_c66_Quant_Doublets, Navin_hbca_c66_Quant_Doublets_Percent))
colnames(Navin_hbca_c66_Quant) <- c("Navin_hbca_c66_Quant")
rownames(Navin_hbca_c66_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c66_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c66_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c66_Singlets <- subset(Navin_hbca_c66, cells=rownames(Navin_hbca_c66@meta.data)[which(Navin_hbca_c66@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c66_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c66_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c66)
rm(Navin_hbca_c66_Singlets)
rm(Navin_hbca_c66_Quant)
rm(Navin_hbca_c66_Quant_Singlets)
rm(Navin_hbca_c66_Quant_Doublets)
rm(Navin_hbca_c66_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c66)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c66)
rm(sweep.stats.Navin_hbca_c66)
gc()



################################################################################################
########################                Navin_hbca_c67                  ########################
################################################################################################


# Load Data
Navin_hbca_c67 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c67.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c67 <- paramSweep(Navin_hbca_c67, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c67 <- summarizeSweep(sweep.res.Navin_hbca_c67, GT = FALSE)
bcmvn_Navin_hbca_c67 <- find.pK(sweep.stats.Navin_hbca_c67)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c67[which.max(bcmvn_Navin_hbca_c67$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c67@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c67@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c67 <- doubletFinder(Navin_hbca_c67, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c67_Quant <- (Navin_hbca_c67@meta.data$DF.classification == "Singlet")
Navin_hbca_c67_Quant_Singlets <- length(Navin_hbca_c67_Quant[Navin_hbca_c67_Quant== TRUE])
Navin_hbca_c67_Quant_Doublets <- length(Navin_hbca_c67_Quant[Navin_hbca_c67_Quant== FALSE])
Navin_hbca_c67_Quant_Doublets_Percent <- Navin_hbca_c67_Quant_Doublets / (Navin_hbca_c67_Quant_Doublets + Navin_hbca_c67_Quant_Singlets) * 100
Navin_hbca_c67_Quant <- as.data.frame(c(Navin_hbca_c67_Quant_Singlets, Navin_hbca_c67_Quant_Doublets, Navin_hbca_c67_Quant_Doublets_Percent))
colnames(Navin_hbca_c67_Quant) <- c("Navin_hbca_c67_Quant")
rownames(Navin_hbca_c67_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c67_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c67_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c67_Singlets <- subset(Navin_hbca_c67, cells=rownames(Navin_hbca_c67@meta.data)[which(Navin_hbca_c67@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c67_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c67_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c67)
rm(Navin_hbca_c67_Singlets)
rm(Navin_hbca_c67_Quant)
rm(Navin_hbca_c67_Quant_Singlets)
rm(Navin_hbca_c67_Quant_Doublets)
rm(Navin_hbca_c67_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c67)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c67)
rm(sweep.stats.Navin_hbca_c67)
gc()



################################################################################################
########################                Navin_hbca_c68                  ########################
################################################################################################


# Load Data
Navin_hbca_c68 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c68.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c68 <- paramSweep(Navin_hbca_c68, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c68 <- summarizeSweep(sweep.res.Navin_hbca_c68, GT = FALSE)
bcmvn_Navin_hbca_c68 <- find.pK(sweep.stats.Navin_hbca_c68)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c68[which.max(bcmvn_Navin_hbca_c68$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c68@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c68@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c68 <- doubletFinder(Navin_hbca_c68, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c68_Quant <- (Navin_hbca_c68@meta.data$DF.classification == "Singlet")
Navin_hbca_c68_Quant_Singlets <- length(Navin_hbca_c68_Quant[Navin_hbca_c68_Quant== TRUE])
Navin_hbca_c68_Quant_Doublets <- length(Navin_hbca_c68_Quant[Navin_hbca_c68_Quant== FALSE])
Navin_hbca_c68_Quant_Doublets_Percent <- Navin_hbca_c68_Quant_Doublets / (Navin_hbca_c68_Quant_Doublets + Navin_hbca_c68_Quant_Singlets) * 100
Navin_hbca_c68_Quant <- as.data.frame(c(Navin_hbca_c68_Quant_Singlets, Navin_hbca_c68_Quant_Doublets, Navin_hbca_c68_Quant_Doublets_Percent))
colnames(Navin_hbca_c68_Quant) <- c("Navin_hbca_c68_Quant")
rownames(Navin_hbca_c68_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c68_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c68_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c68_Singlets <- subset(Navin_hbca_c68, cells=rownames(Navin_hbca_c68@meta.data)[which(Navin_hbca_c68@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c68_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c68_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c68)
rm(Navin_hbca_c68_Singlets)
rm(Navin_hbca_c68_Quant)
rm(Navin_hbca_c68_Quant_Singlets)
rm(Navin_hbca_c68_Quant_Doublets)
rm(Navin_hbca_c68_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c68)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c68)
rm(sweep.stats.Navin_hbca_c68)
gc()



################################################################################################
########################                Navin_hbca_c69                  ########################
################################################################################################


# Load Data
Navin_hbca_c69 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c69.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c69 <- paramSweep(Navin_hbca_c69, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c69 <- summarizeSweep(sweep.res.Navin_hbca_c69, GT = FALSE)
bcmvn_Navin_hbca_c69 <- find.pK(sweep.stats.Navin_hbca_c69)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c69[which.max(bcmvn_Navin_hbca_c69$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c69@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c69@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c69 <- doubletFinder(Navin_hbca_c69, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c69_Quant <- (Navin_hbca_c69@meta.data$DF.classification == "Singlet")
Navin_hbca_c69_Quant_Singlets <- length(Navin_hbca_c69_Quant[Navin_hbca_c69_Quant== TRUE])
Navin_hbca_c69_Quant_Doublets <- length(Navin_hbca_c69_Quant[Navin_hbca_c69_Quant== FALSE])
Navin_hbca_c69_Quant_Doublets_Percent <- Navin_hbca_c69_Quant_Doublets / (Navin_hbca_c69_Quant_Doublets + Navin_hbca_c69_Quant_Singlets) * 100
Navin_hbca_c69_Quant <- as.data.frame(c(Navin_hbca_c69_Quant_Singlets, Navin_hbca_c69_Quant_Doublets, Navin_hbca_c69_Quant_Doublets_Percent))
colnames(Navin_hbca_c69_Quant) <- c("Navin_hbca_c69_Quant")
rownames(Navin_hbca_c69_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c69_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c69_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c69_Singlets <- subset(Navin_hbca_c69, cells=rownames(Navin_hbca_c69@meta.data)[which(Navin_hbca_c69@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c69_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c69_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c69)
rm(Navin_hbca_c69_Singlets)
rm(Navin_hbca_c69_Quant)
rm(Navin_hbca_c69_Quant_Singlets)
rm(Navin_hbca_c69_Quant_Doublets)
rm(Navin_hbca_c69_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c69)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c69)
rm(sweep.stats.Navin_hbca_c69)
gc()



################################################################################################
########################                Navin_hbca_c70                  ########################
################################################################################################


# Load Data
Navin_hbca_c70 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c70.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c70 <- paramSweep(Navin_hbca_c70, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c70 <- summarizeSweep(sweep.res.Navin_hbca_c70, GT = FALSE)
bcmvn_Navin_hbca_c70 <- find.pK(sweep.stats.Navin_hbca_c70)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c70[which.max(bcmvn_Navin_hbca_c70$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c70@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c70@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c70 <- doubletFinder(Navin_hbca_c70, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c70_Quant <- (Navin_hbca_c70@meta.data$DF.classification == "Singlet")
Navin_hbca_c70_Quant_Singlets <- length(Navin_hbca_c70_Quant[Navin_hbca_c70_Quant== TRUE])
Navin_hbca_c70_Quant_Doublets <- length(Navin_hbca_c70_Quant[Navin_hbca_c70_Quant== FALSE])
Navin_hbca_c70_Quant_Doublets_Percent <- Navin_hbca_c70_Quant_Doublets / (Navin_hbca_c70_Quant_Doublets + Navin_hbca_c70_Quant_Singlets) * 100
Navin_hbca_c70_Quant <- as.data.frame(c(Navin_hbca_c70_Quant_Singlets, Navin_hbca_c70_Quant_Doublets, Navin_hbca_c70_Quant_Doublets_Percent))
colnames(Navin_hbca_c70_Quant) <- c("Navin_hbca_c70_Quant")
rownames(Navin_hbca_c70_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c70_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c70_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c70_Singlets <- subset(Navin_hbca_c70, cells=rownames(Navin_hbca_c70@meta.data)[which(Navin_hbca_c70@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c70_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c70_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c70)
rm(Navin_hbca_c70_Singlets)
rm(Navin_hbca_c70_Quant)
rm(Navin_hbca_c70_Quant_Singlets)
rm(Navin_hbca_c70_Quant_Doublets)
rm(Navin_hbca_c70_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c70)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c70)
rm(sweep.stats.Navin_hbca_c70)
gc()



################################################################################################
########################                Navin_hbca_c71                  ########################
################################################################################################


# Load Data
Navin_hbca_c71 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c71.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c71 <- paramSweep(Navin_hbca_c71, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c71 <- summarizeSweep(sweep.res.Navin_hbca_c71, GT = FALSE)
bcmvn_Navin_hbca_c71 <- find.pK(sweep.stats.Navin_hbca_c71)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c71[which.max(bcmvn_Navin_hbca_c71$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c71@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c71@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c71 <- doubletFinder(Navin_hbca_c71, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c71_Quant <- (Navin_hbca_c71@meta.data$DF.classification == "Singlet")
Navin_hbca_c71_Quant_Singlets <- length(Navin_hbca_c71_Quant[Navin_hbca_c71_Quant== TRUE])
Navin_hbca_c71_Quant_Doublets <- length(Navin_hbca_c71_Quant[Navin_hbca_c71_Quant== FALSE])
Navin_hbca_c71_Quant_Doublets_Percent <- Navin_hbca_c71_Quant_Doublets / (Navin_hbca_c71_Quant_Doublets + Navin_hbca_c71_Quant_Singlets) * 100
Navin_hbca_c71_Quant <- as.data.frame(c(Navin_hbca_c71_Quant_Singlets, Navin_hbca_c71_Quant_Doublets, Navin_hbca_c71_Quant_Doublets_Percent))
colnames(Navin_hbca_c71_Quant) <- c("Navin_hbca_c71_Quant")
rownames(Navin_hbca_c71_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c71_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c71_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c71_Singlets <- subset(Navin_hbca_c71, cells=rownames(Navin_hbca_c71@meta.data)[which(Navin_hbca_c71@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c71_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c71_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c71)
rm(Navin_hbca_c71_Singlets)
rm(Navin_hbca_c71_Quant)
rm(Navin_hbca_c71_Quant_Singlets)
rm(Navin_hbca_c71_Quant_Doublets)
rm(Navin_hbca_c71_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c71)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c71)
rm(sweep.stats.Navin_hbca_c71)
gc()



################################################################################################
########################                Navin_hbca_c72                  ########################
################################################################################################


# Load Data
Navin_hbca_c72 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c72.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c72 <- paramSweep(Navin_hbca_c72, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c72 <- summarizeSweep(sweep.res.Navin_hbca_c72, GT = FALSE)
bcmvn_Navin_hbca_c72 <- find.pK(sweep.stats.Navin_hbca_c72)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c72[which.max(bcmvn_Navin_hbca_c72$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c72@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c72@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c72 <- doubletFinder(Navin_hbca_c72, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c72_Quant <- (Navin_hbca_c72@meta.data$DF.classification == "Singlet")
Navin_hbca_c72_Quant_Singlets <- length(Navin_hbca_c72_Quant[Navin_hbca_c72_Quant== TRUE])
Navin_hbca_c72_Quant_Doublets <- length(Navin_hbca_c72_Quant[Navin_hbca_c72_Quant== FALSE])
Navin_hbca_c72_Quant_Doublets_Percent <- Navin_hbca_c72_Quant_Doublets / (Navin_hbca_c72_Quant_Doublets + Navin_hbca_c72_Quant_Singlets) * 100
Navin_hbca_c72_Quant <- as.data.frame(c(Navin_hbca_c72_Quant_Singlets, Navin_hbca_c72_Quant_Doublets, Navin_hbca_c72_Quant_Doublets_Percent))
colnames(Navin_hbca_c72_Quant) <- c("Navin_hbca_c72_Quant")
rownames(Navin_hbca_c72_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c72_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c72_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c72_Singlets <- subset(Navin_hbca_c72, cells=rownames(Navin_hbca_c72@meta.data)[which(Navin_hbca_c72@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c72_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c72_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c72)
rm(Navin_hbca_c72_Singlets)
rm(Navin_hbca_c72_Quant)
rm(Navin_hbca_c72_Quant_Singlets)
rm(Navin_hbca_c72_Quant_Doublets)
rm(Navin_hbca_c72_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c72)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c72)
rm(sweep.stats.Navin_hbca_c72)
gc()



################################################################################################
########################                Navin_hbca_c73                  ########################
################################################################################################


# Load Data
Navin_hbca_c73 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c73.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c73 <- paramSweep(Navin_hbca_c73, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c73 <- summarizeSweep(sweep.res.Navin_hbca_c73, GT = FALSE)
bcmvn_Navin_hbca_c73 <- find.pK(sweep.stats.Navin_hbca_c73)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c73[which.max(bcmvn_Navin_hbca_c73$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c73@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c73@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c73 <- doubletFinder(Navin_hbca_c73, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c73_Quant <- (Navin_hbca_c73@meta.data$DF.classification == "Singlet")
Navin_hbca_c73_Quant_Singlets <- length(Navin_hbca_c73_Quant[Navin_hbca_c73_Quant== TRUE])
Navin_hbca_c73_Quant_Doublets <- length(Navin_hbca_c73_Quant[Navin_hbca_c73_Quant== FALSE])
Navin_hbca_c73_Quant_Doublets_Percent <- Navin_hbca_c73_Quant_Doublets / (Navin_hbca_c73_Quant_Doublets + Navin_hbca_c73_Quant_Singlets) * 100
Navin_hbca_c73_Quant <- as.data.frame(c(Navin_hbca_c73_Quant_Singlets, Navin_hbca_c73_Quant_Doublets, Navin_hbca_c73_Quant_Doublets_Percent))
colnames(Navin_hbca_c73_Quant) <- c("Navin_hbca_c73_Quant")
rownames(Navin_hbca_c73_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c73_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c73_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c73_Singlets <- subset(Navin_hbca_c73, cells=rownames(Navin_hbca_c73@meta.data)[which(Navin_hbca_c73@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c73_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c73_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c73)
rm(Navin_hbca_c73_Singlets)
rm(Navin_hbca_c73_Quant)
rm(Navin_hbca_c73_Quant_Singlets)
rm(Navin_hbca_c73_Quant_Doublets)
rm(Navin_hbca_c73_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c73)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c73)
rm(sweep.stats.Navin_hbca_c73)
gc()



################################################################################################
########################                Navin_hbca_c74                  ########################
################################################################################################


# Load Data
Navin_hbca_c74 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c74.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c74 <- paramSweep(Navin_hbca_c74, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c74 <- summarizeSweep(sweep.res.Navin_hbca_c74, GT = FALSE)
bcmvn_Navin_hbca_c74 <- find.pK(sweep.stats.Navin_hbca_c74)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c74[which.max(bcmvn_Navin_hbca_c74$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c74@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c74@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c74 <- doubletFinder(Navin_hbca_c74, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c74_Quant <- (Navin_hbca_c74@meta.data$DF.classification == "Singlet")
Navin_hbca_c74_Quant_Singlets <- length(Navin_hbca_c74_Quant[Navin_hbca_c74_Quant== TRUE])
Navin_hbca_c74_Quant_Doublets <- length(Navin_hbca_c74_Quant[Navin_hbca_c74_Quant== FALSE])
Navin_hbca_c74_Quant_Doublets_Percent <- Navin_hbca_c74_Quant_Doublets / (Navin_hbca_c74_Quant_Doublets + Navin_hbca_c74_Quant_Singlets) * 100
Navin_hbca_c74_Quant <- as.data.frame(c(Navin_hbca_c74_Quant_Singlets, Navin_hbca_c74_Quant_Doublets, Navin_hbca_c74_Quant_Doublets_Percent))
colnames(Navin_hbca_c74_Quant) <- c("Navin_hbca_c74_Quant")
rownames(Navin_hbca_c74_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c74_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c74_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c74_Singlets <- subset(Navin_hbca_c74, cells=rownames(Navin_hbca_c74@meta.data)[which(Navin_hbca_c74@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c74_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c74_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c74)
rm(Navin_hbca_c74_Singlets)
rm(Navin_hbca_c74_Quant)
rm(Navin_hbca_c74_Quant_Singlets)
rm(Navin_hbca_c74_Quant_Doublets)
rm(Navin_hbca_c74_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c74)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c74)
rm(sweep.stats.Navin_hbca_c74)
gc()



################################################################################################
########################                Navin_hbca_c75                  ########################
################################################################################################


# Load Data
Navin_hbca_c75 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c75.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c75 <- paramSweep(Navin_hbca_c75, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c75 <- summarizeSweep(sweep.res.Navin_hbca_c75, GT = FALSE)
bcmvn_Navin_hbca_c75 <- find.pK(sweep.stats.Navin_hbca_c75)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c75[which.max(bcmvn_Navin_hbca_c75$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c75@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c75@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c75 <- doubletFinder(Navin_hbca_c75, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c75_Quant <- (Navin_hbca_c75@meta.data$DF.classification == "Singlet")
Navin_hbca_c75_Quant_Singlets <- length(Navin_hbca_c75_Quant[Navin_hbca_c75_Quant== TRUE])
Navin_hbca_c75_Quant_Doublets <- length(Navin_hbca_c75_Quant[Navin_hbca_c75_Quant== FALSE])
Navin_hbca_c75_Quant_Doublets_Percent <- Navin_hbca_c75_Quant_Doublets / (Navin_hbca_c75_Quant_Doublets + Navin_hbca_c75_Quant_Singlets) * 100
Navin_hbca_c75_Quant <- as.data.frame(c(Navin_hbca_c75_Quant_Singlets, Navin_hbca_c75_Quant_Doublets, Navin_hbca_c75_Quant_Doublets_Percent))
colnames(Navin_hbca_c75_Quant) <- c("Navin_hbca_c75_Quant")
rownames(Navin_hbca_c75_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c75_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c75_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c75_Singlets <- subset(Navin_hbca_c75, cells=rownames(Navin_hbca_c75@meta.data)[which(Navin_hbca_c75@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c75_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c75_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c75)
rm(Navin_hbca_c75_Singlets)
rm(Navin_hbca_c75_Quant)
rm(Navin_hbca_c75_Quant_Singlets)
rm(Navin_hbca_c75_Quant_Doublets)
rm(Navin_hbca_c75_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c75)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c75)
rm(sweep.stats.Navin_hbca_c75)
gc()



################################################################################################
########################                Navin_hbca_c76                  ########################
################################################################################################


# Load Data
Navin_hbca_c76 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c76.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c76 <- paramSweep(Navin_hbca_c76, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c76 <- summarizeSweep(sweep.res.Navin_hbca_c76, GT = FALSE)
bcmvn_Navin_hbca_c76 <- find.pK(sweep.stats.Navin_hbca_c76)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c76[which.max(bcmvn_Navin_hbca_c76$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c76@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c76@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c76 <- doubletFinder(Navin_hbca_c76, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c76_Quant <- (Navin_hbca_c76@meta.data$DF.classification == "Singlet")
Navin_hbca_c76_Quant_Singlets <- length(Navin_hbca_c76_Quant[Navin_hbca_c76_Quant== TRUE])
Navin_hbca_c76_Quant_Doublets <- length(Navin_hbca_c76_Quant[Navin_hbca_c76_Quant== FALSE])
Navin_hbca_c76_Quant_Doublets_Percent <- Navin_hbca_c76_Quant_Doublets / (Navin_hbca_c76_Quant_Doublets + Navin_hbca_c76_Quant_Singlets) * 100
Navin_hbca_c76_Quant <- as.data.frame(c(Navin_hbca_c76_Quant_Singlets, Navin_hbca_c76_Quant_Doublets, Navin_hbca_c76_Quant_Doublets_Percent))
colnames(Navin_hbca_c76_Quant) <- c("Navin_hbca_c76_Quant")
rownames(Navin_hbca_c76_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c76_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c76_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c76_Singlets <- subset(Navin_hbca_c76, cells=rownames(Navin_hbca_c76@meta.data)[which(Navin_hbca_c76@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c76_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c76_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c76)
rm(Navin_hbca_c76_Singlets)
rm(Navin_hbca_c76_Quant)
rm(Navin_hbca_c76_Quant_Singlets)
rm(Navin_hbca_c76_Quant_Doublets)
rm(Navin_hbca_c76_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c76)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c76)
rm(sweep.stats.Navin_hbca_c76)
gc()



################################################################################################
########################                Navin_hbca_c77                  ########################
################################################################################################


# Load Data
Navin_hbca_c77 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c77.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c77 <- paramSweep(Navin_hbca_c77, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c77 <- summarizeSweep(sweep.res.Navin_hbca_c77, GT = FALSE)
bcmvn_Navin_hbca_c77 <- find.pK(sweep.stats.Navin_hbca_c77)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c77[which.max(bcmvn_Navin_hbca_c77$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c77@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c77@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c77 <- doubletFinder(Navin_hbca_c77, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c77_Quant <- (Navin_hbca_c77@meta.data$DF.classification == "Singlet")
Navin_hbca_c77_Quant_Singlets <- length(Navin_hbca_c77_Quant[Navin_hbca_c77_Quant== TRUE])
Navin_hbca_c77_Quant_Doublets <- length(Navin_hbca_c77_Quant[Navin_hbca_c77_Quant== FALSE])
Navin_hbca_c77_Quant_Doublets_Percent <- Navin_hbca_c77_Quant_Doublets / (Navin_hbca_c77_Quant_Doublets + Navin_hbca_c77_Quant_Singlets) * 100
Navin_hbca_c77_Quant <- as.data.frame(c(Navin_hbca_c77_Quant_Singlets, Navin_hbca_c77_Quant_Doublets, Navin_hbca_c77_Quant_Doublets_Percent))
colnames(Navin_hbca_c77_Quant) <- c("Navin_hbca_c77_Quant")
rownames(Navin_hbca_c77_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c77_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c77_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c77_Singlets <- subset(Navin_hbca_c77, cells=rownames(Navin_hbca_c77@meta.data)[which(Navin_hbca_c77@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c77_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c77_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c77)
rm(Navin_hbca_c77_Singlets)
rm(Navin_hbca_c77_Quant)
rm(Navin_hbca_c77_Quant_Singlets)
rm(Navin_hbca_c77_Quant_Doublets)
rm(Navin_hbca_c77_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c77)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c77)
rm(sweep.stats.Navin_hbca_c77)
gc()



################################################################################################
########################                Navin_hbca_c78                  ########################
################################################################################################


# Load Data
Navin_hbca_c78 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c78.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c78 <- paramSweep(Navin_hbca_c78, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c78 <- summarizeSweep(sweep.res.Navin_hbca_c78, GT = FALSE)
bcmvn_Navin_hbca_c78 <- find.pK(sweep.stats.Navin_hbca_c78)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c78[which.max(bcmvn_Navin_hbca_c78$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c78@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c78@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c78 <- doubletFinder(Navin_hbca_c78, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c78_Quant <- (Navin_hbca_c78@meta.data$DF.classification == "Singlet")
Navin_hbca_c78_Quant_Singlets <- length(Navin_hbca_c78_Quant[Navin_hbca_c78_Quant== TRUE])
Navin_hbca_c78_Quant_Doublets <- length(Navin_hbca_c78_Quant[Navin_hbca_c78_Quant== FALSE])
Navin_hbca_c78_Quant_Doublets_Percent <- Navin_hbca_c78_Quant_Doublets / (Navin_hbca_c78_Quant_Doublets + Navin_hbca_c78_Quant_Singlets) * 100
Navin_hbca_c78_Quant <- as.data.frame(c(Navin_hbca_c78_Quant_Singlets, Navin_hbca_c78_Quant_Doublets, Navin_hbca_c78_Quant_Doublets_Percent))
colnames(Navin_hbca_c78_Quant) <- c("Navin_hbca_c78_Quant")
rownames(Navin_hbca_c78_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c78_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c78_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c78_Singlets <- subset(Navin_hbca_c78, cells=rownames(Navin_hbca_c78@meta.data)[which(Navin_hbca_c78@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c78_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c78_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c78)
rm(Navin_hbca_c78_Singlets)
rm(Navin_hbca_c78_Quant)
rm(Navin_hbca_c78_Quant_Singlets)
rm(Navin_hbca_c78_Quant_Doublets)
rm(Navin_hbca_c78_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c78)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c78)
rm(sweep.stats.Navin_hbca_c78)
gc()



################################################################################################
########################                Navin_hbca_c79                  ########################
################################################################################################


# Load Data
Navin_hbca_c79 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c79.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c79 <- paramSweep(Navin_hbca_c79, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c79 <- summarizeSweep(sweep.res.Navin_hbca_c79, GT = FALSE)
bcmvn_Navin_hbca_c79 <- find.pK(sweep.stats.Navin_hbca_c79)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c79[which.max(bcmvn_Navin_hbca_c79$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c79@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c79@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c79 <- doubletFinder(Navin_hbca_c79, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c79_Quant <- (Navin_hbca_c79@meta.data$DF.classification == "Singlet")
Navin_hbca_c79_Quant_Singlets <- length(Navin_hbca_c79_Quant[Navin_hbca_c79_Quant== TRUE])
Navin_hbca_c79_Quant_Doublets <- length(Navin_hbca_c79_Quant[Navin_hbca_c79_Quant== FALSE])
Navin_hbca_c79_Quant_Doublets_Percent <- Navin_hbca_c79_Quant_Doublets / (Navin_hbca_c79_Quant_Doublets + Navin_hbca_c79_Quant_Singlets) * 100
Navin_hbca_c79_Quant <- as.data.frame(c(Navin_hbca_c79_Quant_Singlets, Navin_hbca_c79_Quant_Doublets, Navin_hbca_c79_Quant_Doublets_Percent))
colnames(Navin_hbca_c79_Quant) <- c("Navin_hbca_c79_Quant")
rownames(Navin_hbca_c79_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c79_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c79_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c79_Singlets <- subset(Navin_hbca_c79, cells=rownames(Navin_hbca_c79@meta.data)[which(Navin_hbca_c79@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c79_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c79_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c79)
rm(Navin_hbca_c79_Singlets)
rm(Navin_hbca_c79_Quant)
rm(Navin_hbca_c79_Quant_Singlets)
rm(Navin_hbca_c79_Quant_Doublets)
rm(Navin_hbca_c79_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c79)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c79)
rm(sweep.stats.Navin_hbca_c79)
gc()



################################################################################################
########################                Navin_hbca_c80                  ########################
################################################################################################


# Load Data
Navin_hbca_c80 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c80.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c80 <- paramSweep(Navin_hbca_c80, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c80 <- summarizeSweep(sweep.res.Navin_hbca_c80, GT = FALSE)
bcmvn_Navin_hbca_c80 <- find.pK(sweep.stats.Navin_hbca_c80)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c80[which.max(bcmvn_Navin_hbca_c80$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c80@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c80@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c80 <- doubletFinder(Navin_hbca_c80, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c80_Quant <- (Navin_hbca_c80@meta.data$DF.classification == "Singlet")
Navin_hbca_c80_Quant_Singlets <- length(Navin_hbca_c80_Quant[Navin_hbca_c80_Quant== TRUE])
Navin_hbca_c80_Quant_Doublets <- length(Navin_hbca_c80_Quant[Navin_hbca_c80_Quant== FALSE])
Navin_hbca_c80_Quant_Doublets_Percent <- Navin_hbca_c80_Quant_Doublets / (Navin_hbca_c80_Quant_Doublets + Navin_hbca_c80_Quant_Singlets) * 100
Navin_hbca_c80_Quant <- as.data.frame(c(Navin_hbca_c80_Quant_Singlets, Navin_hbca_c80_Quant_Doublets, Navin_hbca_c80_Quant_Doublets_Percent))
colnames(Navin_hbca_c80_Quant) <- c("Navin_hbca_c80_Quant")
rownames(Navin_hbca_c80_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c80_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c80_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c80_Singlets <- subset(Navin_hbca_c80, cells=rownames(Navin_hbca_c80@meta.data)[which(Navin_hbca_c80@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c80_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c80_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c80)
rm(Navin_hbca_c80_Singlets)
rm(Navin_hbca_c80_Quant)
rm(Navin_hbca_c80_Quant_Singlets)
rm(Navin_hbca_c80_Quant_Doublets)
rm(Navin_hbca_c80_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c80)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c80)
rm(sweep.stats.Navin_hbca_c80)
gc()



################################################################################################
########################                Navin_hbca_c81                  ########################
################################################################################################


# Load Data
Navin_hbca_c81 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c81.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c81 <- paramSweep(Navin_hbca_c81, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c81 <- summarizeSweep(sweep.res.Navin_hbca_c81, GT = FALSE)
bcmvn_Navin_hbca_c81 <- find.pK(sweep.stats.Navin_hbca_c81)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c81[which.max(bcmvn_Navin_hbca_c81$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c81@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c81@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c81 <- doubletFinder(Navin_hbca_c81, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c81_Quant <- (Navin_hbca_c81@meta.data$DF.classification == "Singlet")
Navin_hbca_c81_Quant_Singlets <- length(Navin_hbca_c81_Quant[Navin_hbca_c81_Quant== TRUE])
Navin_hbca_c81_Quant_Doublets <- length(Navin_hbca_c81_Quant[Navin_hbca_c81_Quant== FALSE])
Navin_hbca_c81_Quant_Doublets_Percent <- Navin_hbca_c81_Quant_Doublets / (Navin_hbca_c81_Quant_Doublets + Navin_hbca_c81_Quant_Singlets) * 100
Navin_hbca_c81_Quant <- as.data.frame(c(Navin_hbca_c81_Quant_Singlets, Navin_hbca_c81_Quant_Doublets, Navin_hbca_c81_Quant_Doublets_Percent))
colnames(Navin_hbca_c81_Quant) <- c("Navin_hbca_c81_Quant")
rownames(Navin_hbca_c81_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c81_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c81_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c81_Singlets <- subset(Navin_hbca_c81, cells=rownames(Navin_hbca_c81@meta.data)[which(Navin_hbca_c81@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c81_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c81_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c81)
rm(Navin_hbca_c81_Singlets)
rm(Navin_hbca_c81_Quant)
rm(Navin_hbca_c81_Quant_Singlets)
rm(Navin_hbca_c81_Quant_Doublets)
rm(Navin_hbca_c81_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c81)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c81)
rm(sweep.stats.Navin_hbca_c81)
gc()



################################################################################################
########################                Navin_hbca_c82                  ########################
################################################################################################


# Load Data
Navin_hbca_c82 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c82.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c82 <- paramSweep(Navin_hbca_c82, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c82 <- summarizeSweep(sweep.res.Navin_hbca_c82, GT = FALSE)
bcmvn_Navin_hbca_c82 <- find.pK(sweep.stats.Navin_hbca_c82)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c82[which.max(bcmvn_Navin_hbca_c82$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c82@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c82@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c82 <- doubletFinder(Navin_hbca_c82, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c82_Quant <- (Navin_hbca_c82@meta.data$DF.classification == "Singlet")
Navin_hbca_c82_Quant_Singlets <- length(Navin_hbca_c82_Quant[Navin_hbca_c82_Quant== TRUE])
Navin_hbca_c82_Quant_Doublets <- length(Navin_hbca_c82_Quant[Navin_hbca_c82_Quant== FALSE])
Navin_hbca_c82_Quant_Doublets_Percent <- Navin_hbca_c82_Quant_Doublets / (Navin_hbca_c82_Quant_Doublets + Navin_hbca_c82_Quant_Singlets) * 100
Navin_hbca_c82_Quant <- as.data.frame(c(Navin_hbca_c82_Quant_Singlets, Navin_hbca_c82_Quant_Doublets, Navin_hbca_c82_Quant_Doublets_Percent))
colnames(Navin_hbca_c82_Quant) <- c("Navin_hbca_c82_Quant")
rownames(Navin_hbca_c82_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c82_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c82_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c82_Singlets <- subset(Navin_hbca_c82, cells=rownames(Navin_hbca_c82@meta.data)[which(Navin_hbca_c82@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c82_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c82_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c82)
rm(Navin_hbca_c82_Singlets)
rm(Navin_hbca_c82_Quant)
rm(Navin_hbca_c82_Quant_Singlets)
rm(Navin_hbca_c82_Quant_Doublets)
rm(Navin_hbca_c82_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c82)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c82)
rm(sweep.stats.Navin_hbca_c82)
gc()



################################################################################################
########################                Navin_hbca_c83                  ########################
################################################################################################


# Load Data
Navin_hbca_c83 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c83.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c83 <- paramSweep(Navin_hbca_c83, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c83 <- summarizeSweep(sweep.res.Navin_hbca_c83, GT = FALSE)
bcmvn_Navin_hbca_c83 <- find.pK(sweep.stats.Navin_hbca_c83)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c83[which.max(bcmvn_Navin_hbca_c83$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c83@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c83@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c83 <- doubletFinder(Navin_hbca_c83, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c83_Quant <- (Navin_hbca_c83@meta.data$DF.classification == "Singlet")
Navin_hbca_c83_Quant_Singlets <- length(Navin_hbca_c83_Quant[Navin_hbca_c83_Quant== TRUE])
Navin_hbca_c83_Quant_Doublets <- length(Navin_hbca_c83_Quant[Navin_hbca_c83_Quant== FALSE])
Navin_hbca_c83_Quant_Doublets_Percent <- Navin_hbca_c83_Quant_Doublets / (Navin_hbca_c83_Quant_Doublets + Navin_hbca_c83_Quant_Singlets) * 100
Navin_hbca_c83_Quant <- as.data.frame(c(Navin_hbca_c83_Quant_Singlets, Navin_hbca_c83_Quant_Doublets, Navin_hbca_c83_Quant_Doublets_Percent))
colnames(Navin_hbca_c83_Quant) <- c("Navin_hbca_c83_Quant")
rownames(Navin_hbca_c83_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c83_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c83_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c83_Singlets <- subset(Navin_hbca_c83, cells=rownames(Navin_hbca_c83@meta.data)[which(Navin_hbca_c83@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c83_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c83_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c83)
rm(Navin_hbca_c83_Singlets)
rm(Navin_hbca_c83_Quant)
rm(Navin_hbca_c83_Quant_Singlets)
rm(Navin_hbca_c83_Quant_Doublets)
rm(Navin_hbca_c83_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c83)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c83)
rm(sweep.stats.Navin_hbca_c83)
gc()



################################################################################################
########################                Navin_hbca_c84                  ########################
################################################################################################


# Load Data
Navin_hbca_c84 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c84.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c84 <- paramSweep(Navin_hbca_c84, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c84 <- summarizeSweep(sweep.res.Navin_hbca_c84, GT = FALSE)
bcmvn_Navin_hbca_c84 <- find.pK(sweep.stats.Navin_hbca_c84)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c84[which.max(bcmvn_Navin_hbca_c84$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c84@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c84@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c84 <- doubletFinder(Navin_hbca_c84, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c84_Quant <- (Navin_hbca_c84@meta.data$DF.classification == "Singlet")
Navin_hbca_c84_Quant_Singlets <- length(Navin_hbca_c84_Quant[Navin_hbca_c84_Quant== TRUE])
Navin_hbca_c84_Quant_Doublets <- length(Navin_hbca_c84_Quant[Navin_hbca_c84_Quant== FALSE])
Navin_hbca_c84_Quant_Doublets_Percent <- Navin_hbca_c84_Quant_Doublets / (Navin_hbca_c84_Quant_Doublets + Navin_hbca_c84_Quant_Singlets) * 100
Navin_hbca_c84_Quant <- as.data.frame(c(Navin_hbca_c84_Quant_Singlets, Navin_hbca_c84_Quant_Doublets, Navin_hbca_c84_Quant_Doublets_Percent))
colnames(Navin_hbca_c84_Quant) <- c("Navin_hbca_c84_Quant")
rownames(Navin_hbca_c84_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c84_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c84_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c84_Singlets <- subset(Navin_hbca_c84, cells=rownames(Navin_hbca_c84@meta.data)[which(Navin_hbca_c84@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c84_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c84_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c84)
rm(Navin_hbca_c84_Singlets)
rm(Navin_hbca_c84_Quant)
rm(Navin_hbca_c84_Quant_Singlets)
rm(Navin_hbca_c84_Quant_Doublets)
rm(Navin_hbca_c84_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c84)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c84)
rm(sweep.stats.Navin_hbca_c84)
gc()



################################################################################################
########################                Navin_hbca_c85                  ########################
################################################################################################


# Load Data
Navin_hbca_c85 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c85.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c85 <- paramSweep(Navin_hbca_c85, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c85 <- summarizeSweep(sweep.res.Navin_hbca_c85, GT = FALSE)
bcmvn_Navin_hbca_c85 <- find.pK(sweep.stats.Navin_hbca_c85)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c85[which.max(bcmvn_Navin_hbca_c85$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c85@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c85@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c85 <- doubletFinder(Navin_hbca_c85, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c85_Quant <- (Navin_hbca_c85@meta.data$DF.classification == "Singlet")
Navin_hbca_c85_Quant_Singlets <- length(Navin_hbca_c85_Quant[Navin_hbca_c85_Quant== TRUE])
Navin_hbca_c85_Quant_Doublets <- length(Navin_hbca_c85_Quant[Navin_hbca_c85_Quant== FALSE])
Navin_hbca_c85_Quant_Doublets_Percent <- Navin_hbca_c85_Quant_Doublets / (Navin_hbca_c85_Quant_Doublets + Navin_hbca_c85_Quant_Singlets) * 100
Navin_hbca_c85_Quant <- as.data.frame(c(Navin_hbca_c85_Quant_Singlets, Navin_hbca_c85_Quant_Doublets, Navin_hbca_c85_Quant_Doublets_Percent))
colnames(Navin_hbca_c85_Quant) <- c("Navin_hbca_c85_Quant")
rownames(Navin_hbca_c85_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c85_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c85_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c85_Singlets <- subset(Navin_hbca_c85, cells=rownames(Navin_hbca_c85@meta.data)[which(Navin_hbca_c85@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c85_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c85_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c85)
rm(Navin_hbca_c85_Singlets)
rm(Navin_hbca_c85_Quant)
rm(Navin_hbca_c85_Quant_Singlets)
rm(Navin_hbca_c85_Quant_Doublets)
rm(Navin_hbca_c85_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c85)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c85)
rm(sweep.stats.Navin_hbca_c85)
gc()



################################################################################################
########################                Navin_hbca_c86                  ########################
################################################################################################


# Load Data
Navin_hbca_c86 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c86.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c86 <- paramSweep(Navin_hbca_c86, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c86 <- summarizeSweep(sweep.res.Navin_hbca_c86, GT = FALSE)
bcmvn_Navin_hbca_c86 <- find.pK(sweep.stats.Navin_hbca_c86)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c86[which.max(bcmvn_Navin_hbca_c86$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c86@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c86@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c86 <- doubletFinder(Navin_hbca_c86, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c86_Quant <- (Navin_hbca_c86@meta.data$DF.classification == "Singlet")
Navin_hbca_c86_Quant_Singlets <- length(Navin_hbca_c86_Quant[Navin_hbca_c86_Quant== TRUE])
Navin_hbca_c86_Quant_Doublets <- length(Navin_hbca_c86_Quant[Navin_hbca_c86_Quant== FALSE])
Navin_hbca_c86_Quant_Doublets_Percent <- Navin_hbca_c86_Quant_Doublets / (Navin_hbca_c86_Quant_Doublets + Navin_hbca_c86_Quant_Singlets) * 100
Navin_hbca_c86_Quant <- as.data.frame(c(Navin_hbca_c86_Quant_Singlets, Navin_hbca_c86_Quant_Doublets, Navin_hbca_c86_Quant_Doublets_Percent))
colnames(Navin_hbca_c86_Quant) <- c("Navin_hbca_c86_Quant")
rownames(Navin_hbca_c86_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c86_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c86_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c86_Singlets <- subset(Navin_hbca_c86, cells=rownames(Navin_hbca_c86@meta.data)[which(Navin_hbca_c86@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c86_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c86_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c86)
rm(Navin_hbca_c86_Singlets)
rm(Navin_hbca_c86_Quant)
rm(Navin_hbca_c86_Quant_Singlets)
rm(Navin_hbca_c86_Quant_Doublets)
rm(Navin_hbca_c86_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c86)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c86)
rm(sweep.stats.Navin_hbca_c86)
gc()



################################################################################################
########################                Navin_hbca_c87                  ########################
################################################################################################


# Load Data
Navin_hbca_c87 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c87.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c87 <- paramSweep(Navin_hbca_c87, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c87 <- summarizeSweep(sweep.res.Navin_hbca_c87, GT = FALSE)
bcmvn_Navin_hbca_c87 <- find.pK(sweep.stats.Navin_hbca_c87)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c87[which.max(bcmvn_Navin_hbca_c87$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c87@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c87@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c87 <- doubletFinder(Navin_hbca_c87, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c87_Quant <- (Navin_hbca_c87@meta.data$DF.classification == "Singlet")
Navin_hbca_c87_Quant_Singlets <- length(Navin_hbca_c87_Quant[Navin_hbca_c87_Quant== TRUE])
Navin_hbca_c87_Quant_Doublets <- length(Navin_hbca_c87_Quant[Navin_hbca_c87_Quant== FALSE])
Navin_hbca_c87_Quant_Doublets_Percent <- Navin_hbca_c87_Quant_Doublets / (Navin_hbca_c87_Quant_Doublets + Navin_hbca_c87_Quant_Singlets) * 100
Navin_hbca_c87_Quant <- as.data.frame(c(Navin_hbca_c87_Quant_Singlets, Navin_hbca_c87_Quant_Doublets, Navin_hbca_c87_Quant_Doublets_Percent))
colnames(Navin_hbca_c87_Quant) <- c("Navin_hbca_c87_Quant")
rownames(Navin_hbca_c87_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c87_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c87_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c87_Singlets <- subset(Navin_hbca_c87, cells=rownames(Navin_hbca_c87@meta.data)[which(Navin_hbca_c87@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c87_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c87_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c87)
rm(Navin_hbca_c87_Singlets)
rm(Navin_hbca_c87_Quant)
rm(Navin_hbca_c87_Quant_Singlets)
rm(Navin_hbca_c87_Quant_Doublets)
rm(Navin_hbca_c87_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c87)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c87)
rm(sweep.stats.Navin_hbca_c87)
gc()



################################################################################################
########################                Navin_hbca_c88                  ########################
################################################################################################


# Load Data
Navin_hbca_c88 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c88.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c88 <- paramSweep(Navin_hbca_c88, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c88 <- summarizeSweep(sweep.res.Navin_hbca_c88, GT = FALSE)
bcmvn_Navin_hbca_c88 <- find.pK(sweep.stats.Navin_hbca_c88)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c88[which.max(bcmvn_Navin_hbca_c88$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c88@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c88@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c88 <- doubletFinder(Navin_hbca_c88, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c88_Quant <- (Navin_hbca_c88@meta.data$DF.classification == "Singlet")
Navin_hbca_c88_Quant_Singlets <- length(Navin_hbca_c88_Quant[Navin_hbca_c88_Quant== TRUE])
Navin_hbca_c88_Quant_Doublets <- length(Navin_hbca_c88_Quant[Navin_hbca_c88_Quant== FALSE])
Navin_hbca_c88_Quant_Doublets_Percent <- Navin_hbca_c88_Quant_Doublets / (Navin_hbca_c88_Quant_Doublets + Navin_hbca_c88_Quant_Singlets) * 100
Navin_hbca_c88_Quant <- as.data.frame(c(Navin_hbca_c88_Quant_Singlets, Navin_hbca_c88_Quant_Doublets, Navin_hbca_c88_Quant_Doublets_Percent))
colnames(Navin_hbca_c88_Quant) <- c("Navin_hbca_c88_Quant")
rownames(Navin_hbca_c88_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c88_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c88_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c88_Singlets <- subset(Navin_hbca_c88, cells=rownames(Navin_hbca_c88@meta.data)[which(Navin_hbca_c88@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c88_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c88_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c88)
rm(Navin_hbca_c88_Singlets)
rm(Navin_hbca_c88_Quant)
rm(Navin_hbca_c88_Quant_Singlets)
rm(Navin_hbca_c88_Quant_Doublets)
rm(Navin_hbca_c88_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c88)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c88)
rm(sweep.stats.Navin_hbca_c88)
gc()



################################################################################################
########################                Navin_hbca_c89                  ########################
################################################################################################


# Load Data
Navin_hbca_c89 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c89.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c89 <- paramSweep(Navin_hbca_c89, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c89 <- summarizeSweep(sweep.res.Navin_hbca_c89, GT = FALSE)
bcmvn_Navin_hbca_c89 <- find.pK(sweep.stats.Navin_hbca_c89)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c89[which.max(bcmvn_Navin_hbca_c89$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c89@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c89@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c89 <- doubletFinder(Navin_hbca_c89, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c89_Quant <- (Navin_hbca_c89@meta.data$DF.classification == "Singlet")
Navin_hbca_c89_Quant_Singlets <- length(Navin_hbca_c89_Quant[Navin_hbca_c89_Quant== TRUE])
Navin_hbca_c89_Quant_Doublets <- length(Navin_hbca_c89_Quant[Navin_hbca_c89_Quant== FALSE])
Navin_hbca_c89_Quant_Doublets_Percent <- Navin_hbca_c89_Quant_Doublets / (Navin_hbca_c89_Quant_Doublets + Navin_hbca_c89_Quant_Singlets) * 100
Navin_hbca_c89_Quant <- as.data.frame(c(Navin_hbca_c89_Quant_Singlets, Navin_hbca_c89_Quant_Doublets, Navin_hbca_c89_Quant_Doublets_Percent))
colnames(Navin_hbca_c89_Quant) <- c("Navin_hbca_c89_Quant")
rownames(Navin_hbca_c89_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c89_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c89_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c89_Singlets <- subset(Navin_hbca_c89, cells=rownames(Navin_hbca_c89@meta.data)[which(Navin_hbca_c89@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c89_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c89_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c89)
rm(Navin_hbca_c89_Singlets)
rm(Navin_hbca_c89_Quant)
rm(Navin_hbca_c89_Quant_Singlets)
rm(Navin_hbca_c89_Quant_Doublets)
rm(Navin_hbca_c89_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c89)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c89)
rm(sweep.stats.Navin_hbca_c89)
gc()



################################################################################################
########################                Navin_hbca_c90                  ########################
################################################################################################


# Load Data
Navin_hbca_c90 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c90.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c90 <- paramSweep(Navin_hbca_c90, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c90 <- summarizeSweep(sweep.res.Navin_hbca_c90, GT = FALSE)
bcmvn_Navin_hbca_c90 <- find.pK(sweep.stats.Navin_hbca_c90)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c90[which.max(bcmvn_Navin_hbca_c90$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c90@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c90@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c90 <- doubletFinder(Navin_hbca_c90, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c90_Quant <- (Navin_hbca_c90@meta.data$DF.classification == "Singlet")
Navin_hbca_c90_Quant_Singlets <- length(Navin_hbca_c90_Quant[Navin_hbca_c90_Quant== TRUE])
Navin_hbca_c90_Quant_Doublets <- length(Navin_hbca_c90_Quant[Navin_hbca_c90_Quant== FALSE])
Navin_hbca_c90_Quant_Doublets_Percent <- Navin_hbca_c90_Quant_Doublets / (Navin_hbca_c90_Quant_Doublets + Navin_hbca_c90_Quant_Singlets) * 100
Navin_hbca_c90_Quant <- as.data.frame(c(Navin_hbca_c90_Quant_Singlets, Navin_hbca_c90_Quant_Doublets, Navin_hbca_c90_Quant_Doublets_Percent))
colnames(Navin_hbca_c90_Quant) <- c("Navin_hbca_c90_Quant")
rownames(Navin_hbca_c90_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c90_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c90_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c90_Singlets <- subset(Navin_hbca_c90, cells=rownames(Navin_hbca_c90@meta.data)[which(Navin_hbca_c90@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c90_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c90_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c90)
rm(Navin_hbca_c90_Singlets)
rm(Navin_hbca_c90_Quant)
rm(Navin_hbca_c90_Quant_Singlets)
rm(Navin_hbca_c90_Quant_Doublets)
rm(Navin_hbca_c90_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c90)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c90)
rm(sweep.stats.Navin_hbca_c90)
gc()



################################################################################################
########################                Navin_hbca_c91                  ########################
################################################################################################


# Load Data
Navin_hbca_c91 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c91.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c91 <- paramSweep(Navin_hbca_c91, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c91 <- summarizeSweep(sweep.res.Navin_hbca_c91, GT = FALSE)
bcmvn_Navin_hbca_c91 <- find.pK(sweep.stats.Navin_hbca_c91)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c91[which.max(bcmvn_Navin_hbca_c91$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c91@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c91@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c91 <- doubletFinder(Navin_hbca_c91, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c91_Quant <- (Navin_hbca_c91@meta.data$DF.classification == "Singlet")
Navin_hbca_c91_Quant_Singlets <- length(Navin_hbca_c91_Quant[Navin_hbca_c91_Quant== TRUE])
Navin_hbca_c91_Quant_Doublets <- length(Navin_hbca_c91_Quant[Navin_hbca_c91_Quant== FALSE])
Navin_hbca_c91_Quant_Doublets_Percent <- Navin_hbca_c91_Quant_Doublets / (Navin_hbca_c91_Quant_Doublets + Navin_hbca_c91_Quant_Singlets) * 100
Navin_hbca_c91_Quant <- as.data.frame(c(Navin_hbca_c91_Quant_Singlets, Navin_hbca_c91_Quant_Doublets, Navin_hbca_c91_Quant_Doublets_Percent))
colnames(Navin_hbca_c91_Quant) <- c("Navin_hbca_c91_Quant")
rownames(Navin_hbca_c91_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c91_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c91_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c91_Singlets <- subset(Navin_hbca_c91, cells=rownames(Navin_hbca_c91@meta.data)[which(Navin_hbca_c91@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c91_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c91_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c91)
rm(Navin_hbca_c91_Singlets)
rm(Navin_hbca_c91_Quant)
rm(Navin_hbca_c91_Quant_Singlets)
rm(Navin_hbca_c91_Quant_Doublets)
rm(Navin_hbca_c91_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c91)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c91)
rm(sweep.stats.Navin_hbca_c91)
gc()



################################################################################################
########################                Navin_hbca_c92                  ########################
################################################################################################


# Load Data
Navin_hbca_c92 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c92.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c92 <- paramSweep(Navin_hbca_c92, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c92 <- summarizeSweep(sweep.res.Navin_hbca_c92, GT = FALSE)
bcmvn_Navin_hbca_c92 <- find.pK(sweep.stats.Navin_hbca_c92)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c92[which.max(bcmvn_Navin_hbca_c92$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c92@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c92@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c92 <- doubletFinder(Navin_hbca_c92, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c92_Quant <- (Navin_hbca_c92@meta.data$DF.classification == "Singlet")
Navin_hbca_c92_Quant_Singlets <- length(Navin_hbca_c92_Quant[Navin_hbca_c92_Quant== TRUE])
Navin_hbca_c92_Quant_Doublets <- length(Navin_hbca_c92_Quant[Navin_hbca_c92_Quant== FALSE])
Navin_hbca_c92_Quant_Doublets_Percent <- Navin_hbca_c92_Quant_Doublets / (Navin_hbca_c92_Quant_Doublets + Navin_hbca_c92_Quant_Singlets) * 100
Navin_hbca_c92_Quant <- as.data.frame(c(Navin_hbca_c92_Quant_Singlets, Navin_hbca_c92_Quant_Doublets, Navin_hbca_c92_Quant_Doublets_Percent))
colnames(Navin_hbca_c92_Quant) <- c("Navin_hbca_c92_Quant")
rownames(Navin_hbca_c92_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c92_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c92_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c92_Singlets <- subset(Navin_hbca_c92, cells=rownames(Navin_hbca_c92@meta.data)[which(Navin_hbca_c92@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c92_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c92_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c92)
rm(Navin_hbca_c92_Singlets)
rm(Navin_hbca_c92_Quant)
rm(Navin_hbca_c92_Quant_Singlets)
rm(Navin_hbca_c92_Quant_Doublets)
rm(Navin_hbca_c92_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c92)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c92)
rm(sweep.stats.Navin_hbca_c92)
gc()



################################################################################################
########################                Navin_hbca_c93                  ########################
################################################################################################


# Load Data
Navin_hbca_c93 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c93.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c93 <- paramSweep(Navin_hbca_c93, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c93 <- summarizeSweep(sweep.res.Navin_hbca_c93, GT = FALSE)
bcmvn_Navin_hbca_c93 <- find.pK(sweep.stats.Navin_hbca_c93)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c93[which.max(bcmvn_Navin_hbca_c93$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c93@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c93@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c93 <- doubletFinder(Navin_hbca_c93, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c93_Quant <- (Navin_hbca_c93@meta.data$DF.classification == "Singlet")
Navin_hbca_c93_Quant_Singlets <- length(Navin_hbca_c93_Quant[Navin_hbca_c93_Quant== TRUE])
Navin_hbca_c93_Quant_Doublets <- length(Navin_hbca_c93_Quant[Navin_hbca_c93_Quant== FALSE])
Navin_hbca_c93_Quant_Doublets_Percent <- Navin_hbca_c93_Quant_Doublets / (Navin_hbca_c93_Quant_Doublets + Navin_hbca_c93_Quant_Singlets) * 100
Navin_hbca_c93_Quant <- as.data.frame(c(Navin_hbca_c93_Quant_Singlets, Navin_hbca_c93_Quant_Doublets, Navin_hbca_c93_Quant_Doublets_Percent))
colnames(Navin_hbca_c93_Quant) <- c("Navin_hbca_c93_Quant")
rownames(Navin_hbca_c93_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c93_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c93_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c93_Singlets <- subset(Navin_hbca_c93, cells=rownames(Navin_hbca_c93@meta.data)[which(Navin_hbca_c93@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c93_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c93_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c93)
rm(Navin_hbca_c93_Singlets)
rm(Navin_hbca_c93_Quant)
rm(Navin_hbca_c93_Quant_Singlets)
rm(Navin_hbca_c93_Quant_Doublets)
rm(Navin_hbca_c93_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c93)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c93)
rm(sweep.stats.Navin_hbca_c93)
gc()



################################################################################################
########################                Navin_hbca_c94                  ########################
################################################################################################


# Load Data
Navin_hbca_c94 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c94.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c94 <- paramSweep(Navin_hbca_c94, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c94 <- summarizeSweep(sweep.res.Navin_hbca_c94, GT = FALSE)
bcmvn_Navin_hbca_c94 <- find.pK(sweep.stats.Navin_hbca_c94)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c94[which.max(bcmvn_Navin_hbca_c94$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c94@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c94@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c94 <- doubletFinder(Navin_hbca_c94, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c94_Quant <- (Navin_hbca_c94@meta.data$DF.classification == "Singlet")
Navin_hbca_c94_Quant_Singlets <- length(Navin_hbca_c94_Quant[Navin_hbca_c94_Quant== TRUE])
Navin_hbca_c94_Quant_Doublets <- length(Navin_hbca_c94_Quant[Navin_hbca_c94_Quant== FALSE])
Navin_hbca_c94_Quant_Doublets_Percent <- Navin_hbca_c94_Quant_Doublets / (Navin_hbca_c94_Quant_Doublets + Navin_hbca_c94_Quant_Singlets) * 100
Navin_hbca_c94_Quant <- as.data.frame(c(Navin_hbca_c94_Quant_Singlets, Navin_hbca_c94_Quant_Doublets, Navin_hbca_c94_Quant_Doublets_Percent))
colnames(Navin_hbca_c94_Quant) <- c("Navin_hbca_c94_Quant")
rownames(Navin_hbca_c94_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c94_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c94_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c94_Singlets <- subset(Navin_hbca_c94, cells=rownames(Navin_hbca_c94@meta.data)[which(Navin_hbca_c94@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c94_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c94_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c94)
rm(Navin_hbca_c94_Singlets)
rm(Navin_hbca_c94_Quant)
rm(Navin_hbca_c94_Quant_Singlets)
rm(Navin_hbca_c94_Quant_Doublets)
rm(Navin_hbca_c94_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c94)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c94)
rm(sweep.stats.Navin_hbca_c94)
gc()



################################################################################################
########################                Navin_hbca_c95                  ########################
################################################################################################


# Load Data
Navin_hbca_c95 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c95.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c95 <- paramSweep(Navin_hbca_c95, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c95 <- summarizeSweep(sweep.res.Navin_hbca_c95, GT = FALSE)
bcmvn_Navin_hbca_c95 <- find.pK(sweep.stats.Navin_hbca_c95)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c95[which.max(bcmvn_Navin_hbca_c95$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c95@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c95@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c95 <- doubletFinder(Navin_hbca_c95, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c95_Quant <- (Navin_hbca_c95@meta.data$DF.classification == "Singlet")
Navin_hbca_c95_Quant_Singlets <- length(Navin_hbca_c95_Quant[Navin_hbca_c95_Quant== TRUE])
Navin_hbca_c95_Quant_Doublets <- length(Navin_hbca_c95_Quant[Navin_hbca_c95_Quant== FALSE])
Navin_hbca_c95_Quant_Doublets_Percent <- Navin_hbca_c95_Quant_Doublets / (Navin_hbca_c95_Quant_Doublets + Navin_hbca_c95_Quant_Singlets) * 100
Navin_hbca_c95_Quant <- as.data.frame(c(Navin_hbca_c95_Quant_Singlets, Navin_hbca_c95_Quant_Doublets, Navin_hbca_c95_Quant_Doublets_Percent))
colnames(Navin_hbca_c95_Quant) <- c("Navin_hbca_c95_Quant")
rownames(Navin_hbca_c95_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c95_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c95_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c95_Singlets <- subset(Navin_hbca_c95, cells=rownames(Navin_hbca_c95@meta.data)[which(Navin_hbca_c95@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c95_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c95_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c95)
rm(Navin_hbca_c95_Singlets)
rm(Navin_hbca_c95_Quant)
rm(Navin_hbca_c95_Quant_Singlets)
rm(Navin_hbca_c95_Quant_Doublets)
rm(Navin_hbca_c95_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c95)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c95)
rm(sweep.stats.Navin_hbca_c95)
gc()



################################################################################################
########################                Navin_hbca_c96                  ########################
################################################################################################


# Load Data
Navin_hbca_c96 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c96.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c96 <- paramSweep(Navin_hbca_c96, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c96 <- summarizeSweep(sweep.res.Navin_hbca_c96, GT = FALSE)
bcmvn_Navin_hbca_c96 <- find.pK(sweep.stats.Navin_hbca_c96)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c96[which.max(bcmvn_Navin_hbca_c96$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c96@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c96@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c96 <- doubletFinder(Navin_hbca_c96, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c96_Quant <- (Navin_hbca_c96@meta.data$DF.classification == "Singlet")
Navin_hbca_c96_Quant_Singlets <- length(Navin_hbca_c96_Quant[Navin_hbca_c96_Quant== TRUE])
Navin_hbca_c96_Quant_Doublets <- length(Navin_hbca_c96_Quant[Navin_hbca_c96_Quant== FALSE])
Navin_hbca_c96_Quant_Doublets_Percent <- Navin_hbca_c96_Quant_Doublets / (Navin_hbca_c96_Quant_Doublets + Navin_hbca_c96_Quant_Singlets) * 100
Navin_hbca_c96_Quant <- as.data.frame(c(Navin_hbca_c96_Quant_Singlets, Navin_hbca_c96_Quant_Doublets, Navin_hbca_c96_Quant_Doublets_Percent))
colnames(Navin_hbca_c96_Quant) <- c("Navin_hbca_c96_Quant")
rownames(Navin_hbca_c96_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c96_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c96_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c96_Singlets <- subset(Navin_hbca_c96, cells=rownames(Navin_hbca_c96@meta.data)[which(Navin_hbca_c96@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c96_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c96_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c96)
rm(Navin_hbca_c96_Singlets)
rm(Navin_hbca_c96_Quant)
rm(Navin_hbca_c96_Quant_Singlets)
rm(Navin_hbca_c96_Quant_Doublets)
rm(Navin_hbca_c96_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c96)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c96)
rm(sweep.stats.Navin_hbca_c96)
gc()



################################################################################################
########################                Navin_hbca_c97                  ########################
################################################################################################


# Load Data
Navin_hbca_c97 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c97.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c97 <- paramSweep(Navin_hbca_c97, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c97 <- summarizeSweep(sweep.res.Navin_hbca_c97, GT = FALSE)
bcmvn_Navin_hbca_c97 <- find.pK(sweep.stats.Navin_hbca_c97)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c97[which.max(bcmvn_Navin_hbca_c97$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c97@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c97@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c97 <- doubletFinder(Navin_hbca_c97, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c97_Quant <- (Navin_hbca_c97@meta.data$DF.classification == "Singlet")
Navin_hbca_c97_Quant_Singlets <- length(Navin_hbca_c97_Quant[Navin_hbca_c97_Quant== TRUE])
Navin_hbca_c97_Quant_Doublets <- length(Navin_hbca_c97_Quant[Navin_hbca_c97_Quant== FALSE])
Navin_hbca_c97_Quant_Doublets_Percent <- Navin_hbca_c97_Quant_Doublets / (Navin_hbca_c97_Quant_Doublets + Navin_hbca_c97_Quant_Singlets) * 100
Navin_hbca_c97_Quant <- as.data.frame(c(Navin_hbca_c97_Quant_Singlets, Navin_hbca_c97_Quant_Doublets, Navin_hbca_c97_Quant_Doublets_Percent))
colnames(Navin_hbca_c97_Quant) <- c("Navin_hbca_c97_Quant")
rownames(Navin_hbca_c97_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c97_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c97_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c97_Singlets <- subset(Navin_hbca_c97, cells=rownames(Navin_hbca_c97@meta.data)[which(Navin_hbca_c97@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c97_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c97_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c97)
rm(Navin_hbca_c97_Singlets)
rm(Navin_hbca_c97_Quant)
rm(Navin_hbca_c97_Quant_Singlets)
rm(Navin_hbca_c97_Quant_Doublets)
rm(Navin_hbca_c97_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c97)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c97)
rm(sweep.stats.Navin_hbca_c97)
gc()



################################################################################################
########################                Navin_hbca_c98                  ########################
################################################################################################


# Load Data
Navin_hbca_c98 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c98.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c98 <- paramSweep(Navin_hbca_c98, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c98 <- summarizeSweep(sweep.res.Navin_hbca_c98, GT = FALSE)
bcmvn_Navin_hbca_c98 <- find.pK(sweep.stats.Navin_hbca_c98)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c98[which.max(bcmvn_Navin_hbca_c98$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c98@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c98@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c98 <- doubletFinder(Navin_hbca_c98, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c98_Quant <- (Navin_hbca_c98@meta.data$DF.classification == "Singlet")
Navin_hbca_c98_Quant_Singlets <- length(Navin_hbca_c98_Quant[Navin_hbca_c98_Quant== TRUE])
Navin_hbca_c98_Quant_Doublets <- length(Navin_hbca_c98_Quant[Navin_hbca_c98_Quant== FALSE])
Navin_hbca_c98_Quant_Doublets_Percent <- Navin_hbca_c98_Quant_Doublets / (Navin_hbca_c98_Quant_Doublets + Navin_hbca_c98_Quant_Singlets) * 100
Navin_hbca_c98_Quant <- as.data.frame(c(Navin_hbca_c98_Quant_Singlets, Navin_hbca_c98_Quant_Doublets, Navin_hbca_c98_Quant_Doublets_Percent))
colnames(Navin_hbca_c98_Quant) <- c("Navin_hbca_c98_Quant")
rownames(Navin_hbca_c98_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c98_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c98_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c98_Singlets <- subset(Navin_hbca_c98, cells=rownames(Navin_hbca_c98@meta.data)[which(Navin_hbca_c98@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c98_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c98_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c98)
rm(Navin_hbca_c98_Singlets)
rm(Navin_hbca_c98_Quant)
rm(Navin_hbca_c98_Quant_Singlets)
rm(Navin_hbca_c98_Quant_Doublets)
rm(Navin_hbca_c98_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c98)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c98)
rm(sweep.stats.Navin_hbca_c98)
gc()



################################################################################################
########################                Navin_hbca_c99                  ########################
################################################################################################


# Load Data
Navin_hbca_c99 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c99.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c99 <- paramSweep(Navin_hbca_c99, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c99 <- summarizeSweep(sweep.res.Navin_hbca_c99, GT = FALSE)
bcmvn_Navin_hbca_c99 <- find.pK(sweep.stats.Navin_hbca_c99)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c99[which.max(bcmvn_Navin_hbca_c99$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c99@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c99@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c99 <- doubletFinder(Navin_hbca_c99, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c99_Quant <- (Navin_hbca_c99@meta.data$DF.classification == "Singlet")
Navin_hbca_c99_Quant_Singlets <- length(Navin_hbca_c99_Quant[Navin_hbca_c99_Quant== TRUE])
Navin_hbca_c99_Quant_Doublets <- length(Navin_hbca_c99_Quant[Navin_hbca_c99_Quant== FALSE])
Navin_hbca_c99_Quant_Doublets_Percent <- Navin_hbca_c99_Quant_Doublets / (Navin_hbca_c99_Quant_Doublets + Navin_hbca_c99_Quant_Singlets) * 100
Navin_hbca_c99_Quant <- as.data.frame(c(Navin_hbca_c99_Quant_Singlets, Navin_hbca_c99_Quant_Doublets, Navin_hbca_c99_Quant_Doublets_Percent))
colnames(Navin_hbca_c99_Quant) <- c("Navin_hbca_c99_Quant")
rownames(Navin_hbca_c99_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c99_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c99_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c99_Singlets <- subset(Navin_hbca_c99, cells=rownames(Navin_hbca_c99@meta.data)[which(Navin_hbca_c99@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c99_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c99_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c99)
rm(Navin_hbca_c99_Singlets)
rm(Navin_hbca_c99_Quant)
rm(Navin_hbca_c99_Quant_Singlets)
rm(Navin_hbca_c99_Quant_Doublets)
rm(Navin_hbca_c99_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c99)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c99)
rm(sweep.stats.Navin_hbca_c99)
gc()



################################################################################################
########################                Navin_hbca_c100                  #######################
################################################################################################


# Load Data
Navin_hbca_c100 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c100.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c100 <- paramSweep(Navin_hbca_c100, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c100 <- summarizeSweep(sweep.res.Navin_hbca_c100, GT = FALSE)
bcmvn_Navin_hbca_c100 <- find.pK(sweep.stats.Navin_hbca_c100)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c100[which.max(bcmvn_Navin_hbca_c100$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c100@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c100@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c100 <- doubletFinder(Navin_hbca_c100, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c100_Quant <- (Navin_hbca_c100@meta.data$DF.classification == "Singlet")
Navin_hbca_c100_Quant_Singlets <- length(Navin_hbca_c100_Quant[Navin_hbca_c100_Quant== TRUE])
Navin_hbca_c100_Quant_Doublets <- length(Navin_hbca_c100_Quant[Navin_hbca_c100_Quant== FALSE])
Navin_hbca_c100_Quant_Doublets_Percent <- Navin_hbca_c100_Quant_Doublets / (Navin_hbca_c100_Quant_Doublets + Navin_hbca_c100_Quant_Singlets) * 100
Navin_hbca_c100_Quant <- as.data.frame(c(Navin_hbca_c100_Quant_Singlets, Navin_hbca_c100_Quant_Doublets, Navin_hbca_c100_Quant_Doublets_Percent))
colnames(Navin_hbca_c100_Quant) <- c("Navin_hbca_c100_Quant")
rownames(Navin_hbca_c100_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c100_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c100_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c100_Singlets <- subset(Navin_hbca_c100, cells=rownames(Navin_hbca_c100@meta.data)[which(Navin_hbca_c100@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c100_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c100_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c100)
rm(Navin_hbca_c100_Singlets)
rm(Navin_hbca_c100_Quant)
rm(Navin_hbca_c100_Quant_Singlets)
rm(Navin_hbca_c100_Quant_Doublets)
rm(Navin_hbca_c100_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c100)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c100)
rm(sweep.stats.Navin_hbca_c100)
gc()



################################################################################################
########################                Navin_hbca_c101                  #######################
################################################################################################


# Load Data
Navin_hbca_c101 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c101.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c101 <- paramSweep(Navin_hbca_c101, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c101 <- summarizeSweep(sweep.res.Navin_hbca_c101, GT = FALSE)
bcmvn_Navin_hbca_c101 <- find.pK(sweep.stats.Navin_hbca_c101)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c101[which.max(bcmvn_Navin_hbca_c101$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c101@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c101@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c101 <- doubletFinder(Navin_hbca_c101, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c101_Quant <- (Navin_hbca_c101@meta.data$DF.classification == "Singlet")
Navin_hbca_c101_Quant_Singlets <- length(Navin_hbca_c101_Quant[Navin_hbca_c101_Quant== TRUE])
Navin_hbca_c101_Quant_Doublets <- length(Navin_hbca_c101_Quant[Navin_hbca_c101_Quant== FALSE])
Navin_hbca_c101_Quant_Doublets_Percent <- Navin_hbca_c101_Quant_Doublets / (Navin_hbca_c101_Quant_Doublets + Navin_hbca_c101_Quant_Singlets) * 100
Navin_hbca_c101_Quant <- as.data.frame(c(Navin_hbca_c101_Quant_Singlets, Navin_hbca_c101_Quant_Doublets, Navin_hbca_c101_Quant_Doublets_Percent))
colnames(Navin_hbca_c101_Quant) <- c("Navin_hbca_c101_Quant")
rownames(Navin_hbca_c101_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c101_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c101_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c101_Singlets <- subset(Navin_hbca_c101, cells=rownames(Navin_hbca_c101@meta.data)[which(Navin_hbca_c101@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c101_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c101_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c101)
rm(Navin_hbca_c101_Singlets)
rm(Navin_hbca_c101_Quant)
rm(Navin_hbca_c101_Quant_Singlets)
rm(Navin_hbca_c101_Quant_Doublets)
rm(Navin_hbca_c101_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c101)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c101)
rm(sweep.stats.Navin_hbca_c101)
gc()



################################################################################################
########################                Navin_hbca_c102                  #######################
################################################################################################


# Load Data
Navin_hbca_c102 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c102.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c102 <- paramSweep(Navin_hbca_c102, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c102 <- summarizeSweep(sweep.res.Navin_hbca_c102, GT = FALSE)
bcmvn_Navin_hbca_c102 <- find.pK(sweep.stats.Navin_hbca_c102)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c102[which.max(bcmvn_Navin_hbca_c102$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c102@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c102@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c102 <- doubletFinder(Navin_hbca_c102, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c102_Quant <- (Navin_hbca_c102@meta.data$DF.classification == "Singlet")
Navin_hbca_c102_Quant_Singlets <- length(Navin_hbca_c102_Quant[Navin_hbca_c102_Quant== TRUE])
Navin_hbca_c102_Quant_Doublets <- length(Navin_hbca_c102_Quant[Navin_hbca_c102_Quant== FALSE])
Navin_hbca_c102_Quant_Doublets_Percent <- Navin_hbca_c102_Quant_Doublets / (Navin_hbca_c102_Quant_Doublets + Navin_hbca_c102_Quant_Singlets) * 100
Navin_hbca_c102_Quant <- as.data.frame(c(Navin_hbca_c102_Quant_Singlets, Navin_hbca_c102_Quant_Doublets, Navin_hbca_c102_Quant_Doublets_Percent))
colnames(Navin_hbca_c102_Quant) <- c("Navin_hbca_c102_Quant")
rownames(Navin_hbca_c102_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c102_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c102_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c102_Singlets <- subset(Navin_hbca_c102, cells=rownames(Navin_hbca_c102@meta.data)[which(Navin_hbca_c102@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c102_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c102_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c102)
rm(Navin_hbca_c102_Singlets)
rm(Navin_hbca_c102_Quant)
rm(Navin_hbca_c102_Quant_Singlets)
rm(Navin_hbca_c102_Quant_Doublets)
rm(Navin_hbca_c102_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c102)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c102)
rm(sweep.stats.Navin_hbca_c102)
gc()



################################################################################################
########################                Navin_hbca_c103                  #######################
################################################################################################


# Load Data
Navin_hbca_c103 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c103.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c103 <- paramSweep(Navin_hbca_c103, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c103 <- summarizeSweep(sweep.res.Navin_hbca_c103, GT = FALSE)
bcmvn_Navin_hbca_c103 <- find.pK(sweep.stats.Navin_hbca_c103)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c103[which.max(bcmvn_Navin_hbca_c103$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c103@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c103@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c103 <- doubletFinder(Navin_hbca_c103, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c103_Quant <- (Navin_hbca_c103@meta.data$DF.classification == "Singlet")
Navin_hbca_c103_Quant_Singlets <- length(Navin_hbca_c103_Quant[Navin_hbca_c103_Quant== TRUE])
Navin_hbca_c103_Quant_Doublets <- length(Navin_hbca_c103_Quant[Navin_hbca_c103_Quant== FALSE])
Navin_hbca_c103_Quant_Doublets_Percent <- Navin_hbca_c103_Quant_Doublets / (Navin_hbca_c103_Quant_Doublets + Navin_hbca_c103_Quant_Singlets) * 100
Navin_hbca_c103_Quant <- as.data.frame(c(Navin_hbca_c103_Quant_Singlets, Navin_hbca_c103_Quant_Doublets, Navin_hbca_c103_Quant_Doublets_Percent))
colnames(Navin_hbca_c103_Quant) <- c("Navin_hbca_c103_Quant")
rownames(Navin_hbca_c103_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c103_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c103_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c103_Singlets <- subset(Navin_hbca_c103, cells=rownames(Navin_hbca_c103@meta.data)[which(Navin_hbca_c103@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c103_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c103_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c103)
rm(Navin_hbca_c103_Singlets)
rm(Navin_hbca_c103_Quant)
rm(Navin_hbca_c103_Quant_Singlets)
rm(Navin_hbca_c103_Quant_Doublets)
rm(Navin_hbca_c103_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c103)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c103)
rm(sweep.stats.Navin_hbca_c103)
gc()



################################################################################################
########################                Navin_hbca_c104                  #######################
################################################################################################


# Load Data
Navin_hbca_c104 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c104.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c104 <- paramSweep(Navin_hbca_c104, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c104 <- summarizeSweep(sweep.res.Navin_hbca_c104, GT = FALSE)
bcmvn_Navin_hbca_c104 <- find.pK(sweep.stats.Navin_hbca_c104)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c104[which.max(bcmvn_Navin_hbca_c104$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c104@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c104@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c104 <- doubletFinder(Navin_hbca_c104, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c104_Quant <- (Navin_hbca_c104@meta.data$DF.classification == "Singlet")
Navin_hbca_c104_Quant_Singlets <- length(Navin_hbca_c104_Quant[Navin_hbca_c104_Quant== TRUE])
Navin_hbca_c104_Quant_Doublets <- length(Navin_hbca_c104_Quant[Navin_hbca_c104_Quant== FALSE])
Navin_hbca_c104_Quant_Doublets_Percent <- Navin_hbca_c104_Quant_Doublets / (Navin_hbca_c104_Quant_Doublets + Navin_hbca_c104_Quant_Singlets) * 100
Navin_hbca_c104_Quant <- as.data.frame(c(Navin_hbca_c104_Quant_Singlets, Navin_hbca_c104_Quant_Doublets, Navin_hbca_c104_Quant_Doublets_Percent))
colnames(Navin_hbca_c104_Quant) <- c("Navin_hbca_c104_Quant")
rownames(Navin_hbca_c104_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c104_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c104_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c104_Singlets <- subset(Navin_hbca_c104, cells=rownames(Navin_hbca_c104@meta.data)[which(Navin_hbca_c104@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c104_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c104_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c104)
rm(Navin_hbca_c104_Singlets)
rm(Navin_hbca_c104_Quant)
rm(Navin_hbca_c104_Quant_Singlets)
rm(Navin_hbca_c104_Quant_Doublets)
rm(Navin_hbca_c104_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c104)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c104)
rm(sweep.stats.Navin_hbca_c104)
gc()



################################################################################################
########################                Navin_hbca_c105                  #######################
################################################################################################


# Load Data
Navin_hbca_c105 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c105.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c105 <- paramSweep(Navin_hbca_c105, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c105 <- summarizeSweep(sweep.res.Navin_hbca_c105, GT = FALSE)
bcmvn_Navin_hbca_c105 <- find.pK(sweep.stats.Navin_hbca_c105)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c105[which.max(bcmvn_Navin_hbca_c105$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c105@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c105@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c105 <- doubletFinder(Navin_hbca_c105, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c105_Quant <- (Navin_hbca_c105@meta.data$DF.classification == "Singlet")
Navin_hbca_c105_Quant_Singlets <- length(Navin_hbca_c105_Quant[Navin_hbca_c105_Quant== TRUE])
Navin_hbca_c105_Quant_Doublets <- length(Navin_hbca_c105_Quant[Navin_hbca_c105_Quant== FALSE])
Navin_hbca_c105_Quant_Doublets_Percent <- Navin_hbca_c105_Quant_Doublets / (Navin_hbca_c105_Quant_Doublets + Navin_hbca_c105_Quant_Singlets) * 100
Navin_hbca_c105_Quant <- as.data.frame(c(Navin_hbca_c105_Quant_Singlets, Navin_hbca_c105_Quant_Doublets, Navin_hbca_c105_Quant_Doublets_Percent))
colnames(Navin_hbca_c105_Quant) <- c("Navin_hbca_c105_Quant")
rownames(Navin_hbca_c105_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c105_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c105_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c105_Singlets <- subset(Navin_hbca_c105, cells=rownames(Navin_hbca_c105@meta.data)[which(Navin_hbca_c105@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c105_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c105_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c105)
rm(Navin_hbca_c105_Singlets)
rm(Navin_hbca_c105_Quant)
rm(Navin_hbca_c105_Quant_Singlets)
rm(Navin_hbca_c105_Quant_Doublets)
rm(Navin_hbca_c105_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c105)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c105)
rm(sweep.stats.Navin_hbca_c105)
gc()



################################################################################################
########################                Navin_hbca_c106                  #######################
################################################################################################


# Load Data
Navin_hbca_c106 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c106.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c106 <- paramSweep(Navin_hbca_c106, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c106 <- summarizeSweep(sweep.res.Navin_hbca_c106, GT = FALSE)
bcmvn_Navin_hbca_c106 <- find.pK(sweep.stats.Navin_hbca_c106)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c106[which.max(bcmvn_Navin_hbca_c106$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c106@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c106@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c106 <- doubletFinder(Navin_hbca_c106, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c106_Quant <- (Navin_hbca_c106@meta.data$DF.classification == "Singlet")
Navin_hbca_c106_Quant_Singlets <- length(Navin_hbca_c106_Quant[Navin_hbca_c106_Quant== TRUE])
Navin_hbca_c106_Quant_Doublets <- length(Navin_hbca_c106_Quant[Navin_hbca_c106_Quant== FALSE])
Navin_hbca_c106_Quant_Doublets_Percent <- Navin_hbca_c106_Quant_Doublets / (Navin_hbca_c106_Quant_Doublets + Navin_hbca_c106_Quant_Singlets) * 100
Navin_hbca_c106_Quant <- as.data.frame(c(Navin_hbca_c106_Quant_Singlets, Navin_hbca_c106_Quant_Doublets, Navin_hbca_c106_Quant_Doublets_Percent))
colnames(Navin_hbca_c106_Quant) <- c("Navin_hbca_c106_Quant")
rownames(Navin_hbca_c106_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c106_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c106_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c106_Singlets <- subset(Navin_hbca_c106, cells=rownames(Navin_hbca_c106@meta.data)[which(Navin_hbca_c106@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c106_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c106_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c106)
rm(Navin_hbca_c106_Singlets)
rm(Navin_hbca_c106_Quant)
rm(Navin_hbca_c106_Quant_Singlets)
rm(Navin_hbca_c106_Quant_Doublets)
rm(Navin_hbca_c106_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c106)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c106)
rm(sweep.stats.Navin_hbca_c106)
gc()



################################################################################################
########################                Navin_hbca_c107                  #######################
################################################################################################


# Load Data
Navin_hbca_c107 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c107.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c107 <- paramSweep(Navin_hbca_c107, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c107 <- summarizeSweep(sweep.res.Navin_hbca_c107, GT = FALSE)
bcmvn_Navin_hbca_c107 <- find.pK(sweep.stats.Navin_hbca_c107)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c107[which.max(bcmvn_Navin_hbca_c107$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c107@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c107@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c107 <- doubletFinder(Navin_hbca_c107, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c107_Quant <- (Navin_hbca_c107@meta.data$DF.classification == "Singlet")
Navin_hbca_c107_Quant_Singlets <- length(Navin_hbca_c107_Quant[Navin_hbca_c107_Quant== TRUE])
Navin_hbca_c107_Quant_Doublets <- length(Navin_hbca_c107_Quant[Navin_hbca_c107_Quant== FALSE])
Navin_hbca_c107_Quant_Doublets_Percent <- Navin_hbca_c107_Quant_Doublets / (Navin_hbca_c107_Quant_Doublets + Navin_hbca_c107_Quant_Singlets) * 100
Navin_hbca_c107_Quant <- as.data.frame(c(Navin_hbca_c107_Quant_Singlets, Navin_hbca_c107_Quant_Doublets, Navin_hbca_c107_Quant_Doublets_Percent))
colnames(Navin_hbca_c107_Quant) <- c("Navin_hbca_c107_Quant")
rownames(Navin_hbca_c107_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c107_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c107_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c107_Singlets <- subset(Navin_hbca_c107, cells=rownames(Navin_hbca_c107@meta.data)[which(Navin_hbca_c107@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c107_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c107_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c107)
rm(Navin_hbca_c107_Singlets)
rm(Navin_hbca_c107_Quant)
rm(Navin_hbca_c107_Quant_Singlets)
rm(Navin_hbca_c107_Quant_Doublets)
rm(Navin_hbca_c107_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c107)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c107)
rm(sweep.stats.Navin_hbca_c107)
gc()



################################################################################################
########################                Navin_hbca_c108                  #######################
################################################################################################


# Load Data
Navin_hbca_c108 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c108.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c108 <- paramSweep(Navin_hbca_c108, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c108 <- summarizeSweep(sweep.res.Navin_hbca_c108, GT = FALSE)
bcmvn_Navin_hbca_c108 <- find.pK(sweep.stats.Navin_hbca_c108)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c108[which.max(bcmvn_Navin_hbca_c108$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c108@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c108@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c108 <- doubletFinder(Navin_hbca_c108, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c108_Quant <- (Navin_hbca_c108@meta.data$DF.classification == "Singlet")
Navin_hbca_c108_Quant_Singlets <- length(Navin_hbca_c108_Quant[Navin_hbca_c108_Quant== TRUE])
Navin_hbca_c108_Quant_Doublets <- length(Navin_hbca_c108_Quant[Navin_hbca_c108_Quant== FALSE])
Navin_hbca_c108_Quant_Doublets_Percent <- Navin_hbca_c108_Quant_Doublets / (Navin_hbca_c108_Quant_Doublets + Navin_hbca_c108_Quant_Singlets) * 100
Navin_hbca_c108_Quant <- as.data.frame(c(Navin_hbca_c108_Quant_Singlets, Navin_hbca_c108_Quant_Doublets, Navin_hbca_c108_Quant_Doublets_Percent))
colnames(Navin_hbca_c108_Quant) <- c("Navin_hbca_c108_Quant")
rownames(Navin_hbca_c108_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c108_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c108_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c108_Singlets <- subset(Navin_hbca_c108, cells=rownames(Navin_hbca_c108@meta.data)[which(Navin_hbca_c108@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c108_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c108_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c108)
rm(Navin_hbca_c108_Singlets)
rm(Navin_hbca_c108_Quant)
rm(Navin_hbca_c108_Quant_Singlets)
rm(Navin_hbca_c108_Quant_Doublets)
rm(Navin_hbca_c108_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c108)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c108)
rm(sweep.stats.Navin_hbca_c108)
gc()



################################################################################################
########################                Navin_hbca_c109                  #######################
################################################################################################


# Load Data
Navin_hbca_c109 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c109.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c109 <- paramSweep(Navin_hbca_c109, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c109 <- summarizeSweep(sweep.res.Navin_hbca_c109, GT = FALSE)
bcmvn_Navin_hbca_c109 <- find.pK(sweep.stats.Navin_hbca_c109)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c109[which.max(bcmvn_Navin_hbca_c109$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c109@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c109@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c109 <- doubletFinder(Navin_hbca_c109, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c109_Quant <- (Navin_hbca_c109@meta.data$DF.classification == "Singlet")
Navin_hbca_c109_Quant_Singlets <- length(Navin_hbca_c109_Quant[Navin_hbca_c109_Quant== TRUE])
Navin_hbca_c109_Quant_Doublets <- length(Navin_hbca_c109_Quant[Navin_hbca_c109_Quant== FALSE])
Navin_hbca_c109_Quant_Doublets_Percent <- Navin_hbca_c109_Quant_Doublets / (Navin_hbca_c109_Quant_Doublets + Navin_hbca_c109_Quant_Singlets) * 100
Navin_hbca_c109_Quant <- as.data.frame(c(Navin_hbca_c109_Quant_Singlets, Navin_hbca_c109_Quant_Doublets, Navin_hbca_c109_Quant_Doublets_Percent))
colnames(Navin_hbca_c109_Quant) <- c("Navin_hbca_c109_Quant")
rownames(Navin_hbca_c109_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c109_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c109_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c109_Singlets <- subset(Navin_hbca_c109, cells=rownames(Navin_hbca_c109@meta.data)[which(Navin_hbca_c109@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c109_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c109_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c109)
rm(Navin_hbca_c109_Singlets)
rm(Navin_hbca_c109_Quant)
rm(Navin_hbca_c109_Quant_Singlets)
rm(Navin_hbca_c109_Quant_Doublets)
rm(Navin_hbca_c109_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c109)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c109)
rm(sweep.stats.Navin_hbca_c109)
gc()



################################################################################################
########################                Navin_hbca_c110                  #######################
################################################################################################


# Load Data
Navin_hbca_c110 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c110.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c110 <- paramSweep(Navin_hbca_c110, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c110 <- summarizeSweep(sweep.res.Navin_hbca_c110, GT = FALSE)
bcmvn_Navin_hbca_c110 <- find.pK(sweep.stats.Navin_hbca_c110)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c110[which.max(bcmvn_Navin_hbca_c110$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c110@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c110@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c110 <- doubletFinder(Navin_hbca_c110, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c110_Quant <- (Navin_hbca_c110@meta.data$DF.classification == "Singlet")
Navin_hbca_c110_Quant_Singlets <- length(Navin_hbca_c110_Quant[Navin_hbca_c110_Quant== TRUE])
Navin_hbca_c110_Quant_Doublets <- length(Navin_hbca_c110_Quant[Navin_hbca_c110_Quant== FALSE])
Navin_hbca_c110_Quant_Doublets_Percent <- Navin_hbca_c110_Quant_Doublets / (Navin_hbca_c110_Quant_Doublets + Navin_hbca_c110_Quant_Singlets) * 100
Navin_hbca_c110_Quant <- as.data.frame(c(Navin_hbca_c110_Quant_Singlets, Navin_hbca_c110_Quant_Doublets, Navin_hbca_c110_Quant_Doublets_Percent))
colnames(Navin_hbca_c110_Quant) <- c("Navin_hbca_c110_Quant")
rownames(Navin_hbca_c110_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c110_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c110_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c110_Singlets <- subset(Navin_hbca_c110, cells=rownames(Navin_hbca_c110@meta.data)[which(Navin_hbca_c110@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c110_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c110_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c110)
rm(Navin_hbca_c110_Singlets)
rm(Navin_hbca_c110_Quant)
rm(Navin_hbca_c110_Quant_Singlets)
rm(Navin_hbca_c110_Quant_Doublets)
rm(Navin_hbca_c110_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c110)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c110)
rm(sweep.stats.Navin_hbca_c110)
gc()



################################################################################################
########################                Navin_hbca_c111                  #######################
################################################################################################


# Load Data
Navin_hbca_c111 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c111.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c111 <- paramSweep(Navin_hbca_c111, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c111 <- summarizeSweep(sweep.res.Navin_hbca_c111, GT = FALSE)
bcmvn_Navin_hbca_c111 <- find.pK(sweep.stats.Navin_hbca_c111)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c111[which.max(bcmvn_Navin_hbca_c111$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c111@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c111@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c111 <- doubletFinder(Navin_hbca_c111, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c111_Quant <- (Navin_hbca_c111@meta.data$DF.classification == "Singlet")
Navin_hbca_c111_Quant_Singlets <- length(Navin_hbca_c111_Quant[Navin_hbca_c111_Quant== TRUE])
Navin_hbca_c111_Quant_Doublets <- length(Navin_hbca_c111_Quant[Navin_hbca_c111_Quant== FALSE])
Navin_hbca_c111_Quant_Doublets_Percent <- Navin_hbca_c111_Quant_Doublets / (Navin_hbca_c111_Quant_Doublets + Navin_hbca_c111_Quant_Singlets) * 100
Navin_hbca_c111_Quant <- as.data.frame(c(Navin_hbca_c111_Quant_Singlets, Navin_hbca_c111_Quant_Doublets, Navin_hbca_c111_Quant_Doublets_Percent))
colnames(Navin_hbca_c111_Quant) <- c("Navin_hbca_c111_Quant")
rownames(Navin_hbca_c111_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c111_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c111_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c111_Singlets <- subset(Navin_hbca_c111, cells=rownames(Navin_hbca_c111@meta.data)[which(Navin_hbca_c111@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c111_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c111_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c111)
rm(Navin_hbca_c111_Singlets)
rm(Navin_hbca_c111_Quant)
rm(Navin_hbca_c111_Quant_Singlets)
rm(Navin_hbca_c111_Quant_Doublets)
rm(Navin_hbca_c111_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c111)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c111)
rm(sweep.stats.Navin_hbca_c111)
gc()



################################################################################################
########################                Navin_hbca_c113                  #######################
################################################################################################


# Load Data
Navin_hbca_c113 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c113.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c113 <- paramSweep(Navin_hbca_c113, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c113 <- summarizeSweep(sweep.res.Navin_hbca_c113, GT = FALSE)
bcmvn_Navin_hbca_c113 <- find.pK(sweep.stats.Navin_hbca_c113)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c113[which.max(bcmvn_Navin_hbca_c113$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c113@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c113@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c113 <- doubletFinder(Navin_hbca_c113, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c113_Quant <- (Navin_hbca_c113@meta.data$DF.classification == "Singlet")
Navin_hbca_c113_Quant_Singlets <- length(Navin_hbca_c113_Quant[Navin_hbca_c113_Quant== TRUE])
Navin_hbca_c113_Quant_Doublets <- length(Navin_hbca_c113_Quant[Navin_hbca_c113_Quant== FALSE])
Navin_hbca_c113_Quant_Doublets_Percent <- Navin_hbca_c113_Quant_Doublets / (Navin_hbca_c113_Quant_Doublets + Navin_hbca_c113_Quant_Singlets) * 100
Navin_hbca_c113_Quant <- as.data.frame(c(Navin_hbca_c113_Quant_Singlets, Navin_hbca_c113_Quant_Doublets, Navin_hbca_c113_Quant_Doublets_Percent))
colnames(Navin_hbca_c113_Quant) <- c("Navin_hbca_c113_Quant")
rownames(Navin_hbca_c113_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c113_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c113_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c113_Singlets <- subset(Navin_hbca_c113, cells=rownames(Navin_hbca_c113@meta.data)[which(Navin_hbca_c113@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c113_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c113_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c113)
rm(Navin_hbca_c113_Singlets)
rm(Navin_hbca_c113_Quant)
rm(Navin_hbca_c113_Quant_Singlets)
rm(Navin_hbca_c113_Quant_Doublets)
rm(Navin_hbca_c113_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c113)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c113)
rm(sweep.stats.Navin_hbca_c113)
gc()



################################################################################################
########################                Navin_hbca_c114                  #######################
################################################################################################


# Load Data
Navin_hbca_c114 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c114.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c114 <- paramSweep(Navin_hbca_c114, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c114 <- summarizeSweep(sweep.res.Navin_hbca_c114, GT = FALSE)
bcmvn_Navin_hbca_c114 <- find.pK(sweep.stats.Navin_hbca_c114)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c114[which.max(bcmvn_Navin_hbca_c114$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c114@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c114@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c114 <- doubletFinder(Navin_hbca_c114, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c114_Quant <- (Navin_hbca_c114@meta.data$DF.classification == "Singlet")
Navin_hbca_c114_Quant_Singlets <- length(Navin_hbca_c114_Quant[Navin_hbca_c114_Quant== TRUE])
Navin_hbca_c114_Quant_Doublets <- length(Navin_hbca_c114_Quant[Navin_hbca_c114_Quant== FALSE])
Navin_hbca_c114_Quant_Doublets_Percent <- Navin_hbca_c114_Quant_Doublets / (Navin_hbca_c114_Quant_Doublets + Navin_hbca_c114_Quant_Singlets) * 100
Navin_hbca_c114_Quant <- as.data.frame(c(Navin_hbca_c114_Quant_Singlets, Navin_hbca_c114_Quant_Doublets, Navin_hbca_c114_Quant_Doublets_Percent))
colnames(Navin_hbca_c114_Quant) <- c("Navin_hbca_c114_Quant")
rownames(Navin_hbca_c114_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c114_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c114_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c114_Singlets <- subset(Navin_hbca_c114, cells=rownames(Navin_hbca_c114@meta.data)[which(Navin_hbca_c114@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c114_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c114_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c114)
rm(Navin_hbca_c114_Singlets)
rm(Navin_hbca_c114_Quant)
rm(Navin_hbca_c114_Quant_Singlets)
rm(Navin_hbca_c114_Quant_Doublets)
rm(Navin_hbca_c114_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c114)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c114)
rm(sweep.stats.Navin_hbca_c114)
gc()



################################################################################################
########################                Navin_hbca_c115                  #######################
################################################################################################


# Load Data
Navin_hbca_c115 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c115.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c115 <- paramSweep(Navin_hbca_c115, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c115 <- summarizeSweep(sweep.res.Navin_hbca_c115, GT = FALSE)
bcmvn_Navin_hbca_c115 <- find.pK(sweep.stats.Navin_hbca_c115)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c115[which.max(bcmvn_Navin_hbca_c115$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c115@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c115@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c115 <- doubletFinder(Navin_hbca_c115, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c115_Quant <- (Navin_hbca_c115@meta.data$DF.classification == "Singlet")
Navin_hbca_c115_Quant_Singlets <- length(Navin_hbca_c115_Quant[Navin_hbca_c115_Quant== TRUE])
Navin_hbca_c115_Quant_Doublets <- length(Navin_hbca_c115_Quant[Navin_hbca_c115_Quant== FALSE])
Navin_hbca_c115_Quant_Doublets_Percent <- Navin_hbca_c115_Quant_Doublets / (Navin_hbca_c115_Quant_Doublets + Navin_hbca_c115_Quant_Singlets) * 100
Navin_hbca_c115_Quant <- as.data.frame(c(Navin_hbca_c115_Quant_Singlets, Navin_hbca_c115_Quant_Doublets, Navin_hbca_c115_Quant_Doublets_Percent))
colnames(Navin_hbca_c115_Quant) <- c("Navin_hbca_c115_Quant")
rownames(Navin_hbca_c115_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c115_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c115_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c115_Singlets <- subset(Navin_hbca_c115, cells=rownames(Navin_hbca_c115@meta.data)[which(Navin_hbca_c115@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c115_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c115_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c115)
rm(Navin_hbca_c115_Singlets)
rm(Navin_hbca_c115_Quant)
rm(Navin_hbca_c115_Quant_Singlets)
rm(Navin_hbca_c115_Quant_Doublets)
rm(Navin_hbca_c115_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c115)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c115)
rm(sweep.stats.Navin_hbca_c115)
gc()



################################################################################################
########################                Navin_hbca_c118                  #######################
################################################################################################


# Load Data
Navin_hbca_c118 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c118.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c118 <- paramSweep(Navin_hbca_c118, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c118 <- summarizeSweep(sweep.res.Navin_hbca_c118, GT = FALSE)
bcmvn_Navin_hbca_c118 <- find.pK(sweep.stats.Navin_hbca_c118)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c118[which.max(bcmvn_Navin_hbca_c118$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c118@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c118@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c118 <- doubletFinder(Navin_hbca_c118, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c118_Quant <- (Navin_hbca_c118@meta.data$DF.classification == "Singlet")
Navin_hbca_c118_Quant_Singlets <- length(Navin_hbca_c118_Quant[Navin_hbca_c118_Quant== TRUE])
Navin_hbca_c118_Quant_Doublets <- length(Navin_hbca_c118_Quant[Navin_hbca_c118_Quant== FALSE])
Navin_hbca_c118_Quant_Doublets_Percent <- Navin_hbca_c118_Quant_Doublets / (Navin_hbca_c118_Quant_Doublets + Navin_hbca_c118_Quant_Singlets) * 100
Navin_hbca_c118_Quant <- as.data.frame(c(Navin_hbca_c118_Quant_Singlets, Navin_hbca_c118_Quant_Doublets, Navin_hbca_c118_Quant_Doublets_Percent))
colnames(Navin_hbca_c118_Quant) <- c("Navin_hbca_c118_Quant")
rownames(Navin_hbca_c118_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c118_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c118_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c118_Singlets <- subset(Navin_hbca_c118, cells=rownames(Navin_hbca_c118@meta.data)[which(Navin_hbca_c118@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c118_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c118_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c118)
rm(Navin_hbca_c118_Singlets)
rm(Navin_hbca_c118_Quant)
rm(Navin_hbca_c118_Quant_Singlets)
rm(Navin_hbca_c118_Quant_Doublets)
rm(Navin_hbca_c118_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c118)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c118)
rm(sweep.stats.Navin_hbca_c118)
gc()



################################################################################################
########################                Navin_hbca_c119                  #######################
################################################################################################


# Load Data
Navin_hbca_c119 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c119.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c119 <- paramSweep(Navin_hbca_c119, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c119 <- summarizeSweep(sweep.res.Navin_hbca_c119, GT = FALSE)
bcmvn_Navin_hbca_c119 <- find.pK(sweep.stats.Navin_hbca_c119)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c119[which.max(bcmvn_Navin_hbca_c119$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c119@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c119@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c119 <- doubletFinder(Navin_hbca_c119, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c119_Quant <- (Navin_hbca_c119@meta.data$DF.classification == "Singlet")
Navin_hbca_c119_Quant_Singlets <- length(Navin_hbca_c119_Quant[Navin_hbca_c119_Quant== TRUE])
Navin_hbca_c119_Quant_Doublets <- length(Navin_hbca_c119_Quant[Navin_hbca_c119_Quant== FALSE])
Navin_hbca_c119_Quant_Doublets_Percent <- Navin_hbca_c119_Quant_Doublets / (Navin_hbca_c119_Quant_Doublets + Navin_hbca_c119_Quant_Singlets) * 100
Navin_hbca_c119_Quant <- as.data.frame(c(Navin_hbca_c119_Quant_Singlets, Navin_hbca_c119_Quant_Doublets, Navin_hbca_c119_Quant_Doublets_Percent))
colnames(Navin_hbca_c119_Quant) <- c("Navin_hbca_c119_Quant")
rownames(Navin_hbca_c119_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c119_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c119_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c119_Singlets <- subset(Navin_hbca_c119, cells=rownames(Navin_hbca_c119@meta.data)[which(Navin_hbca_c119@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c119_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c119_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c119)
rm(Navin_hbca_c119_Singlets)
rm(Navin_hbca_c119_Quant)
rm(Navin_hbca_c119_Quant_Singlets)
rm(Navin_hbca_c119_Quant_Doublets)
rm(Navin_hbca_c119_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c119)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c119)
rm(sweep.stats.Navin_hbca_c119)
gc()



################################################################################################
########################                Navin_hbca_c121                  #######################
################################################################################################


# Load Data
Navin_hbca_c121 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c121.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c121 <- paramSweep(Navin_hbca_c121, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c121 <- summarizeSweep(sweep.res.Navin_hbca_c121, GT = FALSE)
bcmvn_Navin_hbca_c121 <- find.pK(sweep.stats.Navin_hbca_c121)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c121[which.max(bcmvn_Navin_hbca_c121$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c121@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c121@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c121 <- doubletFinder(Navin_hbca_c121, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c121_Quant <- (Navin_hbca_c121@meta.data$DF.classification == "Singlet")
Navin_hbca_c121_Quant_Singlets <- length(Navin_hbca_c121_Quant[Navin_hbca_c121_Quant== TRUE])
Navin_hbca_c121_Quant_Doublets <- length(Navin_hbca_c121_Quant[Navin_hbca_c121_Quant== FALSE])
Navin_hbca_c121_Quant_Doublets_Percent <- Navin_hbca_c121_Quant_Doublets / (Navin_hbca_c121_Quant_Doublets + Navin_hbca_c121_Quant_Singlets) * 100
Navin_hbca_c121_Quant <- as.data.frame(c(Navin_hbca_c121_Quant_Singlets, Navin_hbca_c121_Quant_Doublets, Navin_hbca_c121_Quant_Doublets_Percent))
colnames(Navin_hbca_c121_Quant) <- c("Navin_hbca_c121_Quant")
rownames(Navin_hbca_c121_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c121_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c121_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c121_Singlets <- subset(Navin_hbca_c121, cells=rownames(Navin_hbca_c121@meta.data)[which(Navin_hbca_c121@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c121_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c121_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c121)
rm(Navin_hbca_c121_Singlets)
rm(Navin_hbca_c121_Quant)
rm(Navin_hbca_c121_Quant_Singlets)
rm(Navin_hbca_c121_Quant_Doublets)
rm(Navin_hbca_c121_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c121)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c121)
rm(sweep.stats.Navin_hbca_c121)
gc()



################################################################################################
########################                Navin_hbca_c122                  #######################
################################################################################################


# Load Data
Navin_hbca_c122 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c122.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c122 <- paramSweep(Navin_hbca_c122, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c122 <- summarizeSweep(sweep.res.Navin_hbca_c122, GT = FALSE)
bcmvn_Navin_hbca_c122 <- find.pK(sweep.stats.Navin_hbca_c122)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c122[which.max(bcmvn_Navin_hbca_c122$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c122@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c122@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c122 <- doubletFinder(Navin_hbca_c122, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c122_Quant <- (Navin_hbca_c122@meta.data$DF.classification == "Singlet")
Navin_hbca_c122_Quant_Singlets <- length(Navin_hbca_c122_Quant[Navin_hbca_c122_Quant== TRUE])
Navin_hbca_c122_Quant_Doublets <- length(Navin_hbca_c122_Quant[Navin_hbca_c122_Quant== FALSE])
Navin_hbca_c122_Quant_Doublets_Percent <- Navin_hbca_c122_Quant_Doublets / (Navin_hbca_c122_Quant_Doublets + Navin_hbca_c122_Quant_Singlets) * 100
Navin_hbca_c122_Quant <- as.data.frame(c(Navin_hbca_c122_Quant_Singlets, Navin_hbca_c122_Quant_Doublets, Navin_hbca_c122_Quant_Doublets_Percent))
colnames(Navin_hbca_c122_Quant) <- c("Navin_hbca_c122_Quant")
rownames(Navin_hbca_c122_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c122_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c122_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c122_Singlets <- subset(Navin_hbca_c122, cells=rownames(Navin_hbca_c122@meta.data)[which(Navin_hbca_c122@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c122_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c122_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c122)
rm(Navin_hbca_c122_Singlets)
rm(Navin_hbca_c122_Quant)
rm(Navin_hbca_c122_Quant_Singlets)
rm(Navin_hbca_c122_Quant_Doublets)
rm(Navin_hbca_c122_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c122)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c122)
rm(sweep.stats.Navin_hbca_c122)
gc()



################################################################################################
########################                Navin_hbca_c123                  #######################
################################################################################################


# Load Data
Navin_hbca_c123 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c123.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c123 <- paramSweep(Navin_hbca_c123, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c123 <- summarizeSweep(sweep.res.Navin_hbca_c123, GT = FALSE)
bcmvn_Navin_hbca_c123 <- find.pK(sweep.stats.Navin_hbca_c123)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c123[which.max(bcmvn_Navin_hbca_c123$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c123@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c123@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c123 <- doubletFinder(Navin_hbca_c123, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c123_Quant <- (Navin_hbca_c123@meta.data$DF.classification == "Singlet")
Navin_hbca_c123_Quant_Singlets <- length(Navin_hbca_c123_Quant[Navin_hbca_c123_Quant== TRUE])
Navin_hbca_c123_Quant_Doublets <- length(Navin_hbca_c123_Quant[Navin_hbca_c123_Quant== FALSE])
Navin_hbca_c123_Quant_Doublets_Percent <- Navin_hbca_c123_Quant_Doublets / (Navin_hbca_c123_Quant_Doublets + Navin_hbca_c123_Quant_Singlets) * 100
Navin_hbca_c123_Quant <- as.data.frame(c(Navin_hbca_c123_Quant_Singlets, Navin_hbca_c123_Quant_Doublets, Navin_hbca_c123_Quant_Doublets_Percent))
colnames(Navin_hbca_c123_Quant) <- c("Navin_hbca_c123_Quant")
rownames(Navin_hbca_c123_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c123_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c123_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c123_Singlets <- subset(Navin_hbca_c123, cells=rownames(Navin_hbca_c123@meta.data)[which(Navin_hbca_c123@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c123_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c123_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c123)
rm(Navin_hbca_c123_Singlets)
rm(Navin_hbca_c123_Quant)
rm(Navin_hbca_c123_Quant_Singlets)
rm(Navin_hbca_c123_Quant_Doublets)
rm(Navin_hbca_c123_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c123)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c123)
rm(sweep.stats.Navin_hbca_c123)
gc()



################################################################################################
########################                Navin_hbca_c124                  #######################
################################################################################################


# Load Data
Navin_hbca_c124 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c124.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c124 <- paramSweep(Navin_hbca_c124, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c124 <- summarizeSweep(sweep.res.Navin_hbca_c124, GT = FALSE)
bcmvn_Navin_hbca_c124 <- find.pK(sweep.stats.Navin_hbca_c124)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c124[which.max(bcmvn_Navin_hbca_c124$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c124@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c124@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c124 <- doubletFinder(Navin_hbca_c124, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c124_Quant <- (Navin_hbca_c124@meta.data$DF.classification == "Singlet")
Navin_hbca_c124_Quant_Singlets <- length(Navin_hbca_c124_Quant[Navin_hbca_c124_Quant== TRUE])
Navin_hbca_c124_Quant_Doublets <- length(Navin_hbca_c124_Quant[Navin_hbca_c124_Quant== FALSE])
Navin_hbca_c124_Quant_Doublets_Percent <- Navin_hbca_c124_Quant_Doublets / (Navin_hbca_c124_Quant_Doublets + Navin_hbca_c124_Quant_Singlets) * 100
Navin_hbca_c124_Quant <- as.data.frame(c(Navin_hbca_c124_Quant_Singlets, Navin_hbca_c124_Quant_Doublets, Navin_hbca_c124_Quant_Doublets_Percent))
colnames(Navin_hbca_c124_Quant) <- c("Navin_hbca_c124_Quant")
rownames(Navin_hbca_c124_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c124_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c124_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c124_Singlets <- subset(Navin_hbca_c124, cells=rownames(Navin_hbca_c124@meta.data)[which(Navin_hbca_c124@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c124_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c124_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c124)
rm(Navin_hbca_c124_Singlets)
rm(Navin_hbca_c124_Quant)
rm(Navin_hbca_c124_Quant_Singlets)
rm(Navin_hbca_c124_Quant_Doublets)
rm(Navin_hbca_c124_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c124)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c124)
rm(sweep.stats.Navin_hbca_c124)
gc()



################################################################################################
########################                Navin_hbca_c125                  #######################
################################################################################################


# Load Data
Navin_hbca_c125 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c125.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c125 <- paramSweep(Navin_hbca_c125, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c125 <- summarizeSweep(sweep.res.Navin_hbca_c125, GT = FALSE)
bcmvn_Navin_hbca_c125 <- find.pK(sweep.stats.Navin_hbca_c125)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c125[which.max(bcmvn_Navin_hbca_c125$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c125@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c125@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c125 <- doubletFinder(Navin_hbca_c125, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c125_Quant <- (Navin_hbca_c125@meta.data$DF.classification == "Singlet")
Navin_hbca_c125_Quant_Singlets <- length(Navin_hbca_c125_Quant[Navin_hbca_c125_Quant== TRUE])
Navin_hbca_c125_Quant_Doublets <- length(Navin_hbca_c125_Quant[Navin_hbca_c125_Quant== FALSE])
Navin_hbca_c125_Quant_Doublets_Percent <- Navin_hbca_c125_Quant_Doublets / (Navin_hbca_c125_Quant_Doublets + Navin_hbca_c125_Quant_Singlets) * 100
Navin_hbca_c125_Quant <- as.data.frame(c(Navin_hbca_c125_Quant_Singlets, Navin_hbca_c125_Quant_Doublets, Navin_hbca_c125_Quant_Doublets_Percent))
colnames(Navin_hbca_c125_Quant) <- c("Navin_hbca_c125_Quant")
rownames(Navin_hbca_c125_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c125_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c125_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c125_Singlets <- subset(Navin_hbca_c125, cells=rownames(Navin_hbca_c125@meta.data)[which(Navin_hbca_c125@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c125_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c125_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c125)
rm(Navin_hbca_c125_Singlets)
rm(Navin_hbca_c125_Quant)
rm(Navin_hbca_c125_Quant_Singlets)
rm(Navin_hbca_c125_Quant_Doublets)
rm(Navin_hbca_c125_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c125)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c125)
rm(sweep.stats.Navin_hbca_c125)
gc()



################################################################################################
########################                Navin_hbca_c126                  #######################
################################################################################################


# Load Data
Navin_hbca_c126 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c126.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c126 <- paramSweep(Navin_hbca_c126, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c126 <- summarizeSweep(sweep.res.Navin_hbca_c126, GT = FALSE)
bcmvn_Navin_hbca_c126 <- find.pK(sweep.stats.Navin_hbca_c126)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c126[which.max(bcmvn_Navin_hbca_c126$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c126@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c126@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c126 <- doubletFinder(Navin_hbca_c126, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c126_Quant <- (Navin_hbca_c126@meta.data$DF.classification == "Singlet")
Navin_hbca_c126_Quant_Singlets <- length(Navin_hbca_c126_Quant[Navin_hbca_c126_Quant== TRUE])
Navin_hbca_c126_Quant_Doublets <- length(Navin_hbca_c126_Quant[Navin_hbca_c126_Quant== FALSE])
Navin_hbca_c126_Quant_Doublets_Percent <- Navin_hbca_c126_Quant_Doublets / (Navin_hbca_c126_Quant_Doublets + Navin_hbca_c126_Quant_Singlets) * 100
Navin_hbca_c126_Quant <- as.data.frame(c(Navin_hbca_c126_Quant_Singlets, Navin_hbca_c126_Quant_Doublets, Navin_hbca_c126_Quant_Doublets_Percent))
colnames(Navin_hbca_c126_Quant) <- c("Navin_hbca_c126_Quant")
rownames(Navin_hbca_c126_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c126_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c126_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c126_Singlets <- subset(Navin_hbca_c126, cells=rownames(Navin_hbca_c126@meta.data)[which(Navin_hbca_c126@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c126_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c126_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c126)
rm(Navin_hbca_c126_Singlets)
rm(Navin_hbca_c126_Quant)
rm(Navin_hbca_c126_Quant_Singlets)
rm(Navin_hbca_c126_Quant_Doublets)
rm(Navin_hbca_c126_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c126)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c126)
rm(sweep.stats.Navin_hbca_c126)
gc()



################################################################################################
########################                Navin_hbca_c127                  #######################
################################################################################################


# Load Data
Navin_hbca_c127 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c127.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c127 <- paramSweep(Navin_hbca_c127, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c127 <- summarizeSweep(sweep.res.Navin_hbca_c127, GT = FALSE)
bcmvn_Navin_hbca_c127 <- find.pK(sweep.stats.Navin_hbca_c127)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c127[which.max(bcmvn_Navin_hbca_c127$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c127@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c127@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c127 <- doubletFinder(Navin_hbca_c127, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c127_Quant <- (Navin_hbca_c127@meta.data$DF.classification == "Singlet")
Navin_hbca_c127_Quant_Singlets <- length(Navin_hbca_c127_Quant[Navin_hbca_c127_Quant== TRUE])
Navin_hbca_c127_Quant_Doublets <- length(Navin_hbca_c127_Quant[Navin_hbca_c127_Quant== FALSE])
Navin_hbca_c127_Quant_Doublets_Percent <- Navin_hbca_c127_Quant_Doublets / (Navin_hbca_c127_Quant_Doublets + Navin_hbca_c127_Quant_Singlets) * 100
Navin_hbca_c127_Quant <- as.data.frame(c(Navin_hbca_c127_Quant_Singlets, Navin_hbca_c127_Quant_Doublets, Navin_hbca_c127_Quant_Doublets_Percent))
colnames(Navin_hbca_c127_Quant) <- c("Navin_hbca_c127_Quant")
rownames(Navin_hbca_c127_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c127_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c127_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c127_Singlets <- subset(Navin_hbca_c127, cells=rownames(Navin_hbca_c127@meta.data)[which(Navin_hbca_c127@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c127_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c127_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c127)
rm(Navin_hbca_c127_Singlets)
rm(Navin_hbca_c127_Quant)
rm(Navin_hbca_c127_Quant_Singlets)
rm(Navin_hbca_c127_Quant_Doublets)
rm(Navin_hbca_c127_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c127)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c127)
rm(sweep.stats.Navin_hbca_c127)
gc()



################################################################################################
########################                Navin_hbca_c128                  #######################
################################################################################################


# Load Data
Navin_hbca_c128 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c128.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c128 <- paramSweep(Navin_hbca_c128, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c128 <- summarizeSweep(sweep.res.Navin_hbca_c128, GT = FALSE)
bcmvn_Navin_hbca_c128 <- find.pK(sweep.stats.Navin_hbca_c128)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c128[which.max(bcmvn_Navin_hbca_c128$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c128@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c128@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c128 <- doubletFinder(Navin_hbca_c128, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c128_Quant <- (Navin_hbca_c128@meta.data$DF.classification == "Singlet")
Navin_hbca_c128_Quant_Singlets <- length(Navin_hbca_c128_Quant[Navin_hbca_c128_Quant== TRUE])
Navin_hbca_c128_Quant_Doublets <- length(Navin_hbca_c128_Quant[Navin_hbca_c128_Quant== FALSE])
Navin_hbca_c128_Quant_Doublets_Percent <- Navin_hbca_c128_Quant_Doublets / (Navin_hbca_c128_Quant_Doublets + Navin_hbca_c128_Quant_Singlets) * 100
Navin_hbca_c128_Quant <- as.data.frame(c(Navin_hbca_c128_Quant_Singlets, Navin_hbca_c128_Quant_Doublets, Navin_hbca_c128_Quant_Doublets_Percent))
colnames(Navin_hbca_c128_Quant) <- c("Navin_hbca_c128_Quant")
rownames(Navin_hbca_c128_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c128_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c128_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c128_Singlets <- subset(Navin_hbca_c128, cells=rownames(Navin_hbca_c128@meta.data)[which(Navin_hbca_c128@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c128_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c128_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c128)
rm(Navin_hbca_c128_Singlets)
rm(Navin_hbca_c128_Quant)
rm(Navin_hbca_c128_Quant_Singlets)
rm(Navin_hbca_c128_Quant_Doublets)
rm(Navin_hbca_c128_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c128)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c128)
rm(sweep.stats.Navin_hbca_c128)
gc()



################################################################################################
########################                Navin_hbca_c129                  #######################
################################################################################################


# Load Data
Navin_hbca_c129 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c129.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c129 <- paramSweep(Navin_hbca_c129, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c129 <- summarizeSweep(sweep.res.Navin_hbca_c129, GT = FALSE)
bcmvn_Navin_hbca_c129 <- find.pK(sweep.stats.Navin_hbca_c129)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c129[which.max(bcmvn_Navin_hbca_c129$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c129@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c129@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c129 <- doubletFinder(Navin_hbca_c129, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c129_Quant <- (Navin_hbca_c129@meta.data$DF.classification == "Singlet")
Navin_hbca_c129_Quant_Singlets <- length(Navin_hbca_c129_Quant[Navin_hbca_c129_Quant== TRUE])
Navin_hbca_c129_Quant_Doublets <- length(Navin_hbca_c129_Quant[Navin_hbca_c129_Quant== FALSE])
Navin_hbca_c129_Quant_Doublets_Percent <- Navin_hbca_c129_Quant_Doublets / (Navin_hbca_c129_Quant_Doublets + Navin_hbca_c129_Quant_Singlets) * 100
Navin_hbca_c129_Quant <- as.data.frame(c(Navin_hbca_c129_Quant_Singlets, Navin_hbca_c129_Quant_Doublets, Navin_hbca_c129_Quant_Doublets_Percent))
colnames(Navin_hbca_c129_Quant) <- c("Navin_hbca_c129_Quant")
rownames(Navin_hbca_c129_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c129_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c129_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c129_Singlets <- subset(Navin_hbca_c129, cells=rownames(Navin_hbca_c129@meta.data)[which(Navin_hbca_c129@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c129_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c129_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c129)
rm(Navin_hbca_c129_Singlets)
rm(Navin_hbca_c129_Quant)
rm(Navin_hbca_c129_Quant_Singlets)
rm(Navin_hbca_c129_Quant_Doublets)
rm(Navin_hbca_c129_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c129)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c129)
rm(sweep.stats.Navin_hbca_c129)
gc()



################################################################################################
########################                Navin_hbca_c130                  #######################
################################################################################################


# Load Data
Navin_hbca_c130 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c130.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c130 <- paramSweep(Navin_hbca_c130, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c130 <- summarizeSweep(sweep.res.Navin_hbca_c130, GT = FALSE)
bcmvn_Navin_hbca_c130 <- find.pK(sweep.stats.Navin_hbca_c130)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c130[which.max(bcmvn_Navin_hbca_c130$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c130@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c130@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c130 <- doubletFinder(Navin_hbca_c130, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c130_Quant <- (Navin_hbca_c130@meta.data$DF.classification == "Singlet")
Navin_hbca_c130_Quant_Singlets <- length(Navin_hbca_c130_Quant[Navin_hbca_c130_Quant== TRUE])
Navin_hbca_c130_Quant_Doublets <- length(Navin_hbca_c130_Quant[Navin_hbca_c130_Quant== FALSE])
Navin_hbca_c130_Quant_Doublets_Percent <- Navin_hbca_c130_Quant_Doublets / (Navin_hbca_c130_Quant_Doublets + Navin_hbca_c130_Quant_Singlets) * 100
Navin_hbca_c130_Quant <- as.data.frame(c(Navin_hbca_c130_Quant_Singlets, Navin_hbca_c130_Quant_Doublets, Navin_hbca_c130_Quant_Doublets_Percent))
colnames(Navin_hbca_c130_Quant) <- c("Navin_hbca_c130_Quant")
rownames(Navin_hbca_c130_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c130_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c130_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c130_Singlets <- subset(Navin_hbca_c130, cells=rownames(Navin_hbca_c130@meta.data)[which(Navin_hbca_c130@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c130_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c130_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c130)
rm(Navin_hbca_c130_Singlets)
rm(Navin_hbca_c130_Quant)
rm(Navin_hbca_c130_Quant_Singlets)
rm(Navin_hbca_c130_Quant_Doublets)
rm(Navin_hbca_c130_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c130)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c130)
rm(sweep.stats.Navin_hbca_c130)
gc()



################################################################################################
########################                Navin_hbca_c131                  #######################
################################################################################################


# Load Data
Navin_hbca_c131 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c131.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c131 <- paramSweep(Navin_hbca_c131, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c131 <- summarizeSweep(sweep.res.Navin_hbca_c131, GT = FALSE)
bcmvn_Navin_hbca_c131 <- find.pK(sweep.stats.Navin_hbca_c131)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c131[which.max(bcmvn_Navin_hbca_c131$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c131@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c131@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c131 <- doubletFinder(Navin_hbca_c131, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c131_Quant <- (Navin_hbca_c131@meta.data$DF.classification == "Singlet")
Navin_hbca_c131_Quant_Singlets <- length(Navin_hbca_c131_Quant[Navin_hbca_c131_Quant== TRUE])
Navin_hbca_c131_Quant_Doublets <- length(Navin_hbca_c131_Quant[Navin_hbca_c131_Quant== FALSE])
Navin_hbca_c131_Quant_Doublets_Percent <- Navin_hbca_c131_Quant_Doublets / (Navin_hbca_c131_Quant_Doublets + Navin_hbca_c131_Quant_Singlets) * 100
Navin_hbca_c131_Quant <- as.data.frame(c(Navin_hbca_c131_Quant_Singlets, Navin_hbca_c131_Quant_Doublets, Navin_hbca_c131_Quant_Doublets_Percent))
colnames(Navin_hbca_c131_Quant) <- c("Navin_hbca_c131_Quant")
rownames(Navin_hbca_c131_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c131_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c131_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c131_Singlets <- subset(Navin_hbca_c131, cells=rownames(Navin_hbca_c131@meta.data)[which(Navin_hbca_c131@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c131_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c131_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c131)
rm(Navin_hbca_c131_Singlets)
rm(Navin_hbca_c131_Quant)
rm(Navin_hbca_c131_Quant_Singlets)
rm(Navin_hbca_c131_Quant_Doublets)
rm(Navin_hbca_c131_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c131)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c131)
rm(sweep.stats.Navin_hbca_c131)
gc()



################################################################################################
########################                Navin_hbca_c132                  #######################
################################################################################################


# Load Data
Navin_hbca_c132 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c132.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c132 <- paramSweep(Navin_hbca_c132, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c132 <- summarizeSweep(sweep.res.Navin_hbca_c132, GT = FALSE)
bcmvn_Navin_hbca_c132 <- find.pK(sweep.stats.Navin_hbca_c132)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c132[which.max(bcmvn_Navin_hbca_c132$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c132@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c132@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c132 <- doubletFinder(Navin_hbca_c132, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c132_Quant <- (Navin_hbca_c132@meta.data$DF.classification == "Singlet")
Navin_hbca_c132_Quant_Singlets <- length(Navin_hbca_c132_Quant[Navin_hbca_c132_Quant== TRUE])
Navin_hbca_c132_Quant_Doublets <- length(Navin_hbca_c132_Quant[Navin_hbca_c132_Quant== FALSE])
Navin_hbca_c132_Quant_Doublets_Percent <- Navin_hbca_c132_Quant_Doublets / (Navin_hbca_c132_Quant_Doublets + Navin_hbca_c132_Quant_Singlets) * 100
Navin_hbca_c132_Quant <- as.data.frame(c(Navin_hbca_c132_Quant_Singlets, Navin_hbca_c132_Quant_Doublets, Navin_hbca_c132_Quant_Doublets_Percent))
colnames(Navin_hbca_c132_Quant) <- c("Navin_hbca_c132_Quant")
rownames(Navin_hbca_c132_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c132_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c132_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c132_Singlets <- subset(Navin_hbca_c132, cells=rownames(Navin_hbca_c132@meta.data)[which(Navin_hbca_c132@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c132_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c132_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c132)
rm(Navin_hbca_c132_Singlets)
rm(Navin_hbca_c132_Quant)
rm(Navin_hbca_c132_Quant_Singlets)
rm(Navin_hbca_c132_Quant_Doublets)
rm(Navin_hbca_c132_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c132)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c132)
rm(sweep.stats.Navin_hbca_c132)
gc()



################################################################################################
########################                Navin_hbca_c133                  #######################
################################################################################################


# Load Data
Navin_hbca_c133 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c133.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c133 <- paramSweep(Navin_hbca_c133, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c133 <- summarizeSweep(sweep.res.Navin_hbca_c133, GT = FALSE)
bcmvn_Navin_hbca_c133 <- find.pK(sweep.stats.Navin_hbca_c133)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c133[which.max(bcmvn_Navin_hbca_c133$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c133@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c133@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c133 <- doubletFinder(Navin_hbca_c133, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c133_Quant <- (Navin_hbca_c133@meta.data$DF.classification == "Singlet")
Navin_hbca_c133_Quant_Singlets <- length(Navin_hbca_c133_Quant[Navin_hbca_c133_Quant== TRUE])
Navin_hbca_c133_Quant_Doublets <- length(Navin_hbca_c133_Quant[Navin_hbca_c133_Quant== FALSE])
Navin_hbca_c133_Quant_Doublets_Percent <- Navin_hbca_c133_Quant_Doublets / (Navin_hbca_c133_Quant_Doublets + Navin_hbca_c133_Quant_Singlets) * 100
Navin_hbca_c133_Quant <- as.data.frame(c(Navin_hbca_c133_Quant_Singlets, Navin_hbca_c133_Quant_Doublets, Navin_hbca_c133_Quant_Doublets_Percent))
colnames(Navin_hbca_c133_Quant) <- c("Navin_hbca_c133_Quant")
rownames(Navin_hbca_c133_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c133_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c133_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c133_Singlets <- subset(Navin_hbca_c133, cells=rownames(Navin_hbca_c133@meta.data)[which(Navin_hbca_c133@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c133_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c133_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c133)
rm(Navin_hbca_c133_Singlets)
rm(Navin_hbca_c133_Quant)
rm(Navin_hbca_c133_Quant_Singlets)
rm(Navin_hbca_c133_Quant_Doublets)
rm(Navin_hbca_c133_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c133)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c133)
rm(sweep.stats.Navin_hbca_c133)
gc()



################################################################################################
########################                Navin_hbca_c134                  #######################
################################################################################################


# Load Data
Navin_hbca_c134 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c134.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c134 <- paramSweep(Navin_hbca_c134, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c134 <- summarizeSweep(sweep.res.Navin_hbca_c134, GT = FALSE)
bcmvn_Navin_hbca_c134 <- find.pK(sweep.stats.Navin_hbca_c134)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c134[which.max(bcmvn_Navin_hbca_c134$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c134@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c134@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c134 <- doubletFinder(Navin_hbca_c134, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c134_Quant <- (Navin_hbca_c134@meta.data$DF.classification == "Singlet")
Navin_hbca_c134_Quant_Singlets <- length(Navin_hbca_c134_Quant[Navin_hbca_c134_Quant== TRUE])
Navin_hbca_c134_Quant_Doublets <- length(Navin_hbca_c134_Quant[Navin_hbca_c134_Quant== FALSE])
Navin_hbca_c134_Quant_Doublets_Percent <- Navin_hbca_c134_Quant_Doublets / (Navin_hbca_c134_Quant_Doublets + Navin_hbca_c134_Quant_Singlets) * 100
Navin_hbca_c134_Quant <- as.data.frame(c(Navin_hbca_c134_Quant_Singlets, Navin_hbca_c134_Quant_Doublets, Navin_hbca_c134_Quant_Doublets_Percent))
colnames(Navin_hbca_c134_Quant) <- c("Navin_hbca_c134_Quant")
rownames(Navin_hbca_c134_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c134_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c134_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c134_Singlets <- subset(Navin_hbca_c134, cells=rownames(Navin_hbca_c134@meta.data)[which(Navin_hbca_c134@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c134_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c134_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c134)
rm(Navin_hbca_c134_Singlets)
rm(Navin_hbca_c134_Quant)
rm(Navin_hbca_c134_Quant_Singlets)
rm(Navin_hbca_c134_Quant_Doublets)
rm(Navin_hbca_c134_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c134)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c134)
rm(sweep.stats.Navin_hbca_c134)
gc()



################################################################################################
########################                Navin_hbca_c135                  #######################
################################################################################################


# Load Data
Navin_hbca_c135 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c135.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c135 <- paramSweep(Navin_hbca_c135, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c135 <- summarizeSweep(sweep.res.Navin_hbca_c135, GT = FALSE)
bcmvn_Navin_hbca_c135 <- find.pK(sweep.stats.Navin_hbca_c135)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c135[which.max(bcmvn_Navin_hbca_c135$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c135@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c135@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c135 <- doubletFinder(Navin_hbca_c135, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c135_Quant <- (Navin_hbca_c135@meta.data$DF.classification == "Singlet")
Navin_hbca_c135_Quant_Singlets <- length(Navin_hbca_c135_Quant[Navin_hbca_c135_Quant== TRUE])
Navin_hbca_c135_Quant_Doublets <- length(Navin_hbca_c135_Quant[Navin_hbca_c135_Quant== FALSE])
Navin_hbca_c135_Quant_Doublets_Percent <- Navin_hbca_c135_Quant_Doublets / (Navin_hbca_c135_Quant_Doublets + Navin_hbca_c135_Quant_Singlets) * 100
Navin_hbca_c135_Quant <- as.data.frame(c(Navin_hbca_c135_Quant_Singlets, Navin_hbca_c135_Quant_Doublets, Navin_hbca_c135_Quant_Doublets_Percent))
colnames(Navin_hbca_c135_Quant) <- c("Navin_hbca_c135_Quant")
rownames(Navin_hbca_c135_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c135_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c135_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c135_Singlets <- subset(Navin_hbca_c135, cells=rownames(Navin_hbca_c135@meta.data)[which(Navin_hbca_c135@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c135_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c135_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c135)
rm(Navin_hbca_c135_Singlets)
rm(Navin_hbca_c135_Quant)
rm(Navin_hbca_c135_Quant_Singlets)
rm(Navin_hbca_c135_Quant_Doublets)
rm(Navin_hbca_c135_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c135)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c135)
rm(sweep.stats.Navin_hbca_c135)
gc()



################################################################################################
########################                Navin_hbca_c136                  #######################
################################################################################################


# Load Data
Navin_hbca_c136 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c136.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c136 <- paramSweep(Navin_hbca_c136, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c136 <- summarizeSweep(sweep.res.Navin_hbca_c136, GT = FALSE)
bcmvn_Navin_hbca_c136 <- find.pK(sweep.stats.Navin_hbca_c136)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c136[which.max(bcmvn_Navin_hbca_c136$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c136@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c136@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c136 <- doubletFinder(Navin_hbca_c136, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c136_Quant <- (Navin_hbca_c136@meta.data$DF.classification == "Singlet")
Navin_hbca_c136_Quant_Singlets <- length(Navin_hbca_c136_Quant[Navin_hbca_c136_Quant== TRUE])
Navin_hbca_c136_Quant_Doublets <- length(Navin_hbca_c136_Quant[Navin_hbca_c136_Quant== FALSE])
Navin_hbca_c136_Quant_Doublets_Percent <- Navin_hbca_c136_Quant_Doublets / (Navin_hbca_c136_Quant_Doublets + Navin_hbca_c136_Quant_Singlets) * 100
Navin_hbca_c136_Quant <- as.data.frame(c(Navin_hbca_c136_Quant_Singlets, Navin_hbca_c136_Quant_Doublets, Navin_hbca_c136_Quant_Doublets_Percent))
colnames(Navin_hbca_c136_Quant) <- c("Navin_hbca_c136_Quant")
rownames(Navin_hbca_c136_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c136_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c136_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c136_Singlets <- subset(Navin_hbca_c136, cells=rownames(Navin_hbca_c136@meta.data)[which(Navin_hbca_c136@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c136_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c136_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c136)
rm(Navin_hbca_c136_Singlets)
rm(Navin_hbca_c136_Quant)
rm(Navin_hbca_c136_Quant_Singlets)
rm(Navin_hbca_c136_Quant_Doublets)
rm(Navin_hbca_c136_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c136)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c136)
rm(sweep.stats.Navin_hbca_c136)
gc()



################################################################################################
########################                Navin_hbca_c137                  #######################
################################################################################################


# Load Data
Navin_hbca_c137 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c137.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c137 <- paramSweep(Navin_hbca_c137, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c137 <- summarizeSweep(sweep.res.Navin_hbca_c137, GT = FALSE)
bcmvn_Navin_hbca_c137 <- find.pK(sweep.stats.Navin_hbca_c137)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c137[which.max(bcmvn_Navin_hbca_c137$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c137@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c137@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c137 <- doubletFinder(Navin_hbca_c137, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c137_Quant <- (Navin_hbca_c137@meta.data$DF.classification == "Singlet")
Navin_hbca_c137_Quant_Singlets <- length(Navin_hbca_c137_Quant[Navin_hbca_c137_Quant== TRUE])
Navin_hbca_c137_Quant_Doublets <- length(Navin_hbca_c137_Quant[Navin_hbca_c137_Quant== FALSE])
Navin_hbca_c137_Quant_Doublets_Percent <- Navin_hbca_c137_Quant_Doublets / (Navin_hbca_c137_Quant_Doublets + Navin_hbca_c137_Quant_Singlets) * 100
Navin_hbca_c137_Quant <- as.data.frame(c(Navin_hbca_c137_Quant_Singlets, Navin_hbca_c137_Quant_Doublets, Navin_hbca_c137_Quant_Doublets_Percent))
colnames(Navin_hbca_c137_Quant) <- c("Navin_hbca_c137_Quant")
rownames(Navin_hbca_c137_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c137_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c137_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c137_Singlets <- subset(Navin_hbca_c137, cells=rownames(Navin_hbca_c137@meta.data)[which(Navin_hbca_c137@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c137_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c137_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c137)
rm(Navin_hbca_c137_Singlets)
rm(Navin_hbca_c137_Quant)
rm(Navin_hbca_c137_Quant_Singlets)
rm(Navin_hbca_c137_Quant_Doublets)
rm(Navin_hbca_c137_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c137)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c137)
rm(sweep.stats.Navin_hbca_c137)
gc()



################################################################################################
########################                Navin_hbca_c138                  #######################
################################################################################################


# Load Data
Navin_hbca_c138 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c138.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c138 <- paramSweep(Navin_hbca_c138, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c138 <- summarizeSweep(sweep.res.Navin_hbca_c138, GT = FALSE)
bcmvn_Navin_hbca_c138 <- find.pK(sweep.stats.Navin_hbca_c138)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c138[which.max(bcmvn_Navin_hbca_c138$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c138@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c138@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c138 <- doubletFinder(Navin_hbca_c138, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c138_Quant <- (Navin_hbca_c138@meta.data$DF.classification == "Singlet")
Navin_hbca_c138_Quant_Singlets <- length(Navin_hbca_c138_Quant[Navin_hbca_c138_Quant== TRUE])
Navin_hbca_c138_Quant_Doublets <- length(Navin_hbca_c138_Quant[Navin_hbca_c138_Quant== FALSE])
Navin_hbca_c138_Quant_Doublets_Percent <- Navin_hbca_c138_Quant_Doublets / (Navin_hbca_c138_Quant_Doublets + Navin_hbca_c138_Quant_Singlets) * 100
Navin_hbca_c138_Quant <- as.data.frame(c(Navin_hbca_c138_Quant_Singlets, Navin_hbca_c138_Quant_Doublets, Navin_hbca_c138_Quant_Doublets_Percent))
colnames(Navin_hbca_c138_Quant) <- c("Navin_hbca_c138_Quant")
rownames(Navin_hbca_c138_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c138_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c138_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c138_Singlets <- subset(Navin_hbca_c138, cells=rownames(Navin_hbca_c138@meta.data)[which(Navin_hbca_c138@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c138_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c138_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c138)
rm(Navin_hbca_c138_Singlets)
rm(Navin_hbca_c138_Quant)
rm(Navin_hbca_c138_Quant_Singlets)
rm(Navin_hbca_c138_Quant_Doublets)
rm(Navin_hbca_c138_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c138)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c138)
rm(sweep.stats.Navin_hbca_c138)
gc()



################################################################################################
########################                Navin_hbca_c139                  #######################
################################################################################################


# Load Data
Navin_hbca_c139 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c139.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c139 <- paramSweep(Navin_hbca_c139, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c139 <- summarizeSweep(sweep.res.Navin_hbca_c139, GT = FALSE)
bcmvn_Navin_hbca_c139 <- find.pK(sweep.stats.Navin_hbca_c139)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c139[which.max(bcmvn_Navin_hbca_c139$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c139@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c139@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c139 <- doubletFinder(Navin_hbca_c139, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c139_Quant <- (Navin_hbca_c139@meta.data$DF.classification == "Singlet")
Navin_hbca_c139_Quant_Singlets <- length(Navin_hbca_c139_Quant[Navin_hbca_c139_Quant== TRUE])
Navin_hbca_c139_Quant_Doublets <- length(Navin_hbca_c139_Quant[Navin_hbca_c139_Quant== FALSE])
Navin_hbca_c139_Quant_Doublets_Percent <- Navin_hbca_c139_Quant_Doublets / (Navin_hbca_c139_Quant_Doublets + Navin_hbca_c139_Quant_Singlets) * 100
Navin_hbca_c139_Quant <- as.data.frame(c(Navin_hbca_c139_Quant_Singlets, Navin_hbca_c139_Quant_Doublets, Navin_hbca_c139_Quant_Doublets_Percent))
colnames(Navin_hbca_c139_Quant) <- c("Navin_hbca_c139_Quant")
rownames(Navin_hbca_c139_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c139_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c139_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c139_Singlets <- subset(Navin_hbca_c139, cells=rownames(Navin_hbca_c139@meta.data)[which(Navin_hbca_c139@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c139_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c139_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c139)
rm(Navin_hbca_c139_Singlets)
rm(Navin_hbca_c139_Quant)
rm(Navin_hbca_c139_Quant_Singlets)
rm(Navin_hbca_c139_Quant_Doublets)
rm(Navin_hbca_c139_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c139)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c139)
rm(sweep.stats.Navin_hbca_c139)
gc()



################################################################################################
########################                Navin_hbca_c140                  #######################
################################################################################################


# Load Data
Navin_hbca_c140 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c140.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c140 <- paramSweep(Navin_hbca_c140, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c140 <- summarizeSweep(sweep.res.Navin_hbca_c140, GT = FALSE)
bcmvn_Navin_hbca_c140 <- find.pK(sweep.stats.Navin_hbca_c140)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c140[which.max(bcmvn_Navin_hbca_c140$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c140@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c140@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c140 <- doubletFinder(Navin_hbca_c140, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c140_Quant <- (Navin_hbca_c140@meta.data$DF.classification == "Singlet")
Navin_hbca_c140_Quant_Singlets <- length(Navin_hbca_c140_Quant[Navin_hbca_c140_Quant== TRUE])
Navin_hbca_c140_Quant_Doublets <- length(Navin_hbca_c140_Quant[Navin_hbca_c140_Quant== FALSE])
Navin_hbca_c140_Quant_Doublets_Percent <- Navin_hbca_c140_Quant_Doublets / (Navin_hbca_c140_Quant_Doublets + Navin_hbca_c140_Quant_Singlets) * 100
Navin_hbca_c140_Quant <- as.data.frame(c(Navin_hbca_c140_Quant_Singlets, Navin_hbca_c140_Quant_Doublets, Navin_hbca_c140_Quant_Doublets_Percent))
colnames(Navin_hbca_c140_Quant) <- c("Navin_hbca_c140_Quant")
rownames(Navin_hbca_c140_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c140_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c140_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c140_Singlets <- subset(Navin_hbca_c140, cells=rownames(Navin_hbca_c140@meta.data)[which(Navin_hbca_c140@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c140_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c140_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c140)
rm(Navin_hbca_c140_Singlets)
rm(Navin_hbca_c140_Quant)
rm(Navin_hbca_c140_Quant_Singlets)
rm(Navin_hbca_c140_Quant_Doublets)
rm(Navin_hbca_c140_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c140)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c140)
rm(sweep.stats.Navin_hbca_c140)
gc()



################################################################################################
########################                Navin_hbca_c141                  #######################
################################################################################################


# Load Data
Navin_hbca_c141 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c141.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c141 <- paramSweep(Navin_hbca_c141, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c141 <- summarizeSweep(sweep.res.Navin_hbca_c141, GT = FALSE)
bcmvn_Navin_hbca_c141 <- find.pK(sweep.stats.Navin_hbca_c141)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c141[which.max(bcmvn_Navin_hbca_c141$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c141@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c141@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c141 <- doubletFinder(Navin_hbca_c141, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c141_Quant <- (Navin_hbca_c141@meta.data$DF.classification == "Singlet")
Navin_hbca_c141_Quant_Singlets <- length(Navin_hbca_c141_Quant[Navin_hbca_c141_Quant== TRUE])
Navin_hbca_c141_Quant_Doublets <- length(Navin_hbca_c141_Quant[Navin_hbca_c141_Quant== FALSE])
Navin_hbca_c141_Quant_Doublets_Percent <- Navin_hbca_c141_Quant_Doublets / (Navin_hbca_c141_Quant_Doublets + Navin_hbca_c141_Quant_Singlets) * 100
Navin_hbca_c141_Quant <- as.data.frame(c(Navin_hbca_c141_Quant_Singlets, Navin_hbca_c141_Quant_Doublets, Navin_hbca_c141_Quant_Doublets_Percent))
colnames(Navin_hbca_c141_Quant) <- c("Navin_hbca_c141_Quant")
rownames(Navin_hbca_c141_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c141_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c141_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c141_Singlets <- subset(Navin_hbca_c141, cells=rownames(Navin_hbca_c141@meta.data)[which(Navin_hbca_c141@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c141_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c141_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c141)
rm(Navin_hbca_c141_Singlets)
rm(Navin_hbca_c141_Quant)
rm(Navin_hbca_c141_Quant_Singlets)
rm(Navin_hbca_c141_Quant_Doublets)
rm(Navin_hbca_c141_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c141)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c141)
rm(sweep.stats.Navin_hbca_c141)
gc()



################################################################################################
########################                Navin_hbca_c142                  #######################
################################################################################################


# Load Data
Navin_hbca_c142 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c142.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c142 <- paramSweep(Navin_hbca_c142, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c142 <- summarizeSweep(sweep.res.Navin_hbca_c142, GT = FALSE)
bcmvn_Navin_hbca_c142 <- find.pK(sweep.stats.Navin_hbca_c142)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c142[which.max(bcmvn_Navin_hbca_c142$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c142@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c142@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c142 <- doubletFinder(Navin_hbca_c142, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c142_Quant <- (Navin_hbca_c142@meta.data$DF.classification == "Singlet")
Navin_hbca_c142_Quant_Singlets <- length(Navin_hbca_c142_Quant[Navin_hbca_c142_Quant== TRUE])
Navin_hbca_c142_Quant_Doublets <- length(Navin_hbca_c142_Quant[Navin_hbca_c142_Quant== FALSE])
Navin_hbca_c142_Quant_Doublets_Percent <- Navin_hbca_c142_Quant_Doublets / (Navin_hbca_c142_Quant_Doublets + Navin_hbca_c142_Quant_Singlets) * 100
Navin_hbca_c142_Quant <- as.data.frame(c(Navin_hbca_c142_Quant_Singlets, Navin_hbca_c142_Quant_Doublets, Navin_hbca_c142_Quant_Doublets_Percent))
colnames(Navin_hbca_c142_Quant) <- c("Navin_hbca_c142_Quant")
rownames(Navin_hbca_c142_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c142_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c142_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c142_Singlets <- subset(Navin_hbca_c142, cells=rownames(Navin_hbca_c142@meta.data)[which(Navin_hbca_c142@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c142_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c142_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c142)
rm(Navin_hbca_c142_Singlets)
rm(Navin_hbca_c142_Quant)
rm(Navin_hbca_c142_Quant_Singlets)
rm(Navin_hbca_c142_Quant_Doublets)
rm(Navin_hbca_c142_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c142)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c142)
rm(sweep.stats.Navin_hbca_c142)
gc()



################################################################################################
########################                Navin_hbca_c143                  #######################
################################################################################################


# Load Data
Navin_hbca_c143 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c143.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c143 <- paramSweep(Navin_hbca_c143, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c143 <- summarizeSweep(sweep.res.Navin_hbca_c143, GT = FALSE)
bcmvn_Navin_hbca_c143 <- find.pK(sweep.stats.Navin_hbca_c143)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c143[which.max(bcmvn_Navin_hbca_c143$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c143@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c143@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c143 <- doubletFinder(Navin_hbca_c143, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c143_Quant <- (Navin_hbca_c143@meta.data$DF.classification == "Singlet")
Navin_hbca_c143_Quant_Singlets <- length(Navin_hbca_c143_Quant[Navin_hbca_c143_Quant== TRUE])
Navin_hbca_c143_Quant_Doublets <- length(Navin_hbca_c143_Quant[Navin_hbca_c143_Quant== FALSE])
Navin_hbca_c143_Quant_Doublets_Percent <- Navin_hbca_c143_Quant_Doublets / (Navin_hbca_c143_Quant_Doublets + Navin_hbca_c143_Quant_Singlets) * 100
Navin_hbca_c143_Quant <- as.data.frame(c(Navin_hbca_c143_Quant_Singlets, Navin_hbca_c143_Quant_Doublets, Navin_hbca_c143_Quant_Doublets_Percent))
colnames(Navin_hbca_c143_Quant) <- c("Navin_hbca_c143_Quant")
rownames(Navin_hbca_c143_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c143_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c143_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c143_Singlets <- subset(Navin_hbca_c143, cells=rownames(Navin_hbca_c143@meta.data)[which(Navin_hbca_c143@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c143_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c143_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c143)
rm(Navin_hbca_c143_Singlets)
rm(Navin_hbca_c143_Quant)
rm(Navin_hbca_c143_Quant_Singlets)
rm(Navin_hbca_c143_Quant_Doublets)
rm(Navin_hbca_c143_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c143)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c143)
rm(sweep.stats.Navin_hbca_c143)
gc()



################################################################################################
########################                Navin_hbca_c144                  #######################
################################################################################################


# Load Data
Navin_hbca_c144 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c144.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c144 <- paramSweep(Navin_hbca_c144, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c144 <- summarizeSweep(sweep.res.Navin_hbca_c144, GT = FALSE)
bcmvn_Navin_hbca_c144 <- find.pK(sweep.stats.Navin_hbca_c144)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c144[which.max(bcmvn_Navin_hbca_c144$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c144@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c144@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c144 <- doubletFinder(Navin_hbca_c144, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c144_Quant <- (Navin_hbca_c144@meta.data$DF.classification == "Singlet")
Navin_hbca_c144_Quant_Singlets <- length(Navin_hbca_c144_Quant[Navin_hbca_c144_Quant== TRUE])
Navin_hbca_c144_Quant_Doublets <- length(Navin_hbca_c144_Quant[Navin_hbca_c144_Quant== FALSE])
Navin_hbca_c144_Quant_Doublets_Percent <- Navin_hbca_c144_Quant_Doublets / (Navin_hbca_c144_Quant_Doublets + Navin_hbca_c144_Quant_Singlets) * 100
Navin_hbca_c144_Quant <- as.data.frame(c(Navin_hbca_c144_Quant_Singlets, Navin_hbca_c144_Quant_Doublets, Navin_hbca_c144_Quant_Doublets_Percent))
colnames(Navin_hbca_c144_Quant) <- c("Navin_hbca_c144_Quant")
rownames(Navin_hbca_c144_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c144_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c144_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c144_Singlets <- subset(Navin_hbca_c144, cells=rownames(Navin_hbca_c144@meta.data)[which(Navin_hbca_c144@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c144_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c144_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c144)
rm(Navin_hbca_c144_Singlets)
rm(Navin_hbca_c144_Quant)
rm(Navin_hbca_c144_Quant_Singlets)
rm(Navin_hbca_c144_Quant_Doublets)
rm(Navin_hbca_c144_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c144)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c144)
rm(sweep.stats.Navin_hbca_c144)
gc()



################################################################################################
########################                Navin_hbca_c145                  #######################
################################################################################################


# Load Data
Navin_hbca_c145 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c145.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c145 <- paramSweep(Navin_hbca_c145, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c145 <- summarizeSweep(sweep.res.Navin_hbca_c145, GT = FALSE)
bcmvn_Navin_hbca_c145 <- find.pK(sweep.stats.Navin_hbca_c145)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c145[which.max(bcmvn_Navin_hbca_c145$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c145@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c145@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c145 <- doubletFinder(Navin_hbca_c145, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c145_Quant <- (Navin_hbca_c145@meta.data$DF.classification == "Singlet")
Navin_hbca_c145_Quant_Singlets <- length(Navin_hbca_c145_Quant[Navin_hbca_c145_Quant== TRUE])
Navin_hbca_c145_Quant_Doublets <- length(Navin_hbca_c145_Quant[Navin_hbca_c145_Quant== FALSE])
Navin_hbca_c145_Quant_Doublets_Percent <- Navin_hbca_c145_Quant_Doublets / (Navin_hbca_c145_Quant_Doublets + Navin_hbca_c145_Quant_Singlets) * 100
Navin_hbca_c145_Quant <- as.data.frame(c(Navin_hbca_c145_Quant_Singlets, Navin_hbca_c145_Quant_Doublets, Navin_hbca_c145_Quant_Doublets_Percent))
colnames(Navin_hbca_c145_Quant) <- c("Navin_hbca_c145_Quant")
rownames(Navin_hbca_c145_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c145_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c145_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c145_Singlets <- subset(Navin_hbca_c145, cells=rownames(Navin_hbca_c145@meta.data)[which(Navin_hbca_c145@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c145_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c145_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c145)
rm(Navin_hbca_c145_Singlets)
rm(Navin_hbca_c145_Quant)
rm(Navin_hbca_c145_Quant_Singlets)
rm(Navin_hbca_c145_Quant_Doublets)
rm(Navin_hbca_c145_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c145)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c145)
rm(sweep.stats.Navin_hbca_c145)
gc()



################################################################################################
########################                Navin_hbca_c146                  #######################
################################################################################################


# Load Data
Navin_hbca_c146 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c146.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c146 <- paramSweep(Navin_hbca_c146, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c146 <- summarizeSweep(sweep.res.Navin_hbca_c146, GT = FALSE)
bcmvn_Navin_hbca_c146 <- find.pK(sweep.stats.Navin_hbca_c146)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c146[which.max(bcmvn_Navin_hbca_c146$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c146@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c146@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c146 <- doubletFinder(Navin_hbca_c146, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c146_Quant <- (Navin_hbca_c146@meta.data$DF.classification == "Singlet")
Navin_hbca_c146_Quant_Singlets <- length(Navin_hbca_c146_Quant[Navin_hbca_c146_Quant== TRUE])
Navin_hbca_c146_Quant_Doublets <- length(Navin_hbca_c146_Quant[Navin_hbca_c146_Quant== FALSE])
Navin_hbca_c146_Quant_Doublets_Percent <- Navin_hbca_c146_Quant_Doublets / (Navin_hbca_c146_Quant_Doublets + Navin_hbca_c146_Quant_Singlets) * 100
Navin_hbca_c146_Quant <- as.data.frame(c(Navin_hbca_c146_Quant_Singlets, Navin_hbca_c146_Quant_Doublets, Navin_hbca_c146_Quant_Doublets_Percent))
colnames(Navin_hbca_c146_Quant) <- c("Navin_hbca_c146_Quant")
rownames(Navin_hbca_c146_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c146_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c146_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c146_Singlets <- subset(Navin_hbca_c146, cells=rownames(Navin_hbca_c146@meta.data)[which(Navin_hbca_c146@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c146_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c146_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c146)
rm(Navin_hbca_c146_Singlets)
rm(Navin_hbca_c146_Quant)
rm(Navin_hbca_c146_Quant_Singlets)
rm(Navin_hbca_c146_Quant_Doublets)
rm(Navin_hbca_c146_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c146)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c146)
rm(sweep.stats.Navin_hbca_c146)
gc()



################################################################################################
########################                Navin_hbca_c147                  #######################
################################################################################################


# Load Data
Navin_hbca_c147 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c147.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c147 <- paramSweep(Navin_hbca_c147, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c147 <- summarizeSweep(sweep.res.Navin_hbca_c147, GT = FALSE)
bcmvn_Navin_hbca_c147 <- find.pK(sweep.stats.Navin_hbca_c147)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c147[which.max(bcmvn_Navin_hbca_c147$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c147@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c147@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c147 <- doubletFinder(Navin_hbca_c147, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c147_Quant <- (Navin_hbca_c147@meta.data$DF.classification == "Singlet")
Navin_hbca_c147_Quant_Singlets <- length(Navin_hbca_c147_Quant[Navin_hbca_c147_Quant== TRUE])
Navin_hbca_c147_Quant_Doublets <- length(Navin_hbca_c147_Quant[Navin_hbca_c147_Quant== FALSE])
Navin_hbca_c147_Quant_Doublets_Percent <- Navin_hbca_c147_Quant_Doublets / (Navin_hbca_c147_Quant_Doublets + Navin_hbca_c147_Quant_Singlets) * 100
Navin_hbca_c147_Quant <- as.data.frame(c(Navin_hbca_c147_Quant_Singlets, Navin_hbca_c147_Quant_Doublets, Navin_hbca_c147_Quant_Doublets_Percent))
colnames(Navin_hbca_c147_Quant) <- c("Navin_hbca_c147_Quant")
rownames(Navin_hbca_c147_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c147_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c147_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c147_Singlets <- subset(Navin_hbca_c147, cells=rownames(Navin_hbca_c147@meta.data)[which(Navin_hbca_c147@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c147_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c147_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c147)
rm(Navin_hbca_c147_Singlets)
rm(Navin_hbca_c147_Quant)
rm(Navin_hbca_c147_Quant_Singlets)
rm(Navin_hbca_c147_Quant_Doublets)
rm(Navin_hbca_c147_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c147)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c147)
rm(sweep.stats.Navin_hbca_c147)
gc()



################################################################################################
########################                Navin_hbca_c148                  #######################
################################################################################################


# Load Data
Navin_hbca_c148 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c148.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c148 <- paramSweep(Navin_hbca_c148, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c148 <- summarizeSweep(sweep.res.Navin_hbca_c148, GT = FALSE)
bcmvn_Navin_hbca_c148 <- find.pK(sweep.stats.Navin_hbca_c148)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c148[which.max(bcmvn_Navin_hbca_c148$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c148@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c148@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c148 <- doubletFinder(Navin_hbca_c148, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c148_Quant <- (Navin_hbca_c148@meta.data$DF.classification == "Singlet")
Navin_hbca_c148_Quant_Singlets <- length(Navin_hbca_c148_Quant[Navin_hbca_c148_Quant== TRUE])
Navin_hbca_c148_Quant_Doublets <- length(Navin_hbca_c148_Quant[Navin_hbca_c148_Quant== FALSE])
Navin_hbca_c148_Quant_Doublets_Percent <- Navin_hbca_c148_Quant_Doublets / (Navin_hbca_c148_Quant_Doublets + Navin_hbca_c148_Quant_Singlets) * 100
Navin_hbca_c148_Quant <- as.data.frame(c(Navin_hbca_c148_Quant_Singlets, Navin_hbca_c148_Quant_Doublets, Navin_hbca_c148_Quant_Doublets_Percent))
colnames(Navin_hbca_c148_Quant) <- c("Navin_hbca_c148_Quant")
rownames(Navin_hbca_c148_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c148_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c148_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c148_Singlets <- subset(Navin_hbca_c148, cells=rownames(Navin_hbca_c148@meta.data)[which(Navin_hbca_c148@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c148_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c148_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c148)
rm(Navin_hbca_c148_Singlets)
rm(Navin_hbca_c148_Quant)
rm(Navin_hbca_c148_Quant_Singlets)
rm(Navin_hbca_c148_Quant_Doublets)
rm(Navin_hbca_c148_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c148)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c148)
rm(sweep.stats.Navin_hbca_c148)
gc()



################################################################################################
########################                Navin_hbca_c149                  #######################
################################################################################################


# Load Data
Navin_hbca_c149 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c149.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c149 <- paramSweep(Navin_hbca_c149, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c149 <- summarizeSweep(sweep.res.Navin_hbca_c149, GT = FALSE)
bcmvn_Navin_hbca_c149 <- find.pK(sweep.stats.Navin_hbca_c149)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c149[which.max(bcmvn_Navin_hbca_c149$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c149@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c149@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c149 <- doubletFinder(Navin_hbca_c149, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c149_Quant <- (Navin_hbca_c149@meta.data$DF.classification == "Singlet")
Navin_hbca_c149_Quant_Singlets <- length(Navin_hbca_c149_Quant[Navin_hbca_c149_Quant== TRUE])
Navin_hbca_c149_Quant_Doublets <- length(Navin_hbca_c149_Quant[Navin_hbca_c149_Quant== FALSE])
Navin_hbca_c149_Quant_Doublets_Percent <- Navin_hbca_c149_Quant_Doublets / (Navin_hbca_c149_Quant_Doublets + Navin_hbca_c149_Quant_Singlets) * 100
Navin_hbca_c149_Quant <- as.data.frame(c(Navin_hbca_c149_Quant_Singlets, Navin_hbca_c149_Quant_Doublets, Navin_hbca_c149_Quant_Doublets_Percent))
colnames(Navin_hbca_c149_Quant) <- c("Navin_hbca_c149_Quant")
rownames(Navin_hbca_c149_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c149_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c149_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c149_Singlets <- subset(Navin_hbca_c149, cells=rownames(Navin_hbca_c149@meta.data)[which(Navin_hbca_c149@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c149_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c149_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c149)
rm(Navin_hbca_c149_Singlets)
rm(Navin_hbca_c149_Quant)
rm(Navin_hbca_c149_Quant_Singlets)
rm(Navin_hbca_c149_Quant_Doublets)
rm(Navin_hbca_c149_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c149)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c149)
rm(sweep.stats.Navin_hbca_c149)
gc()



################################################################################################
########################                Navin_hbca_c150                  #######################
################################################################################################


# Load Data
Navin_hbca_c150 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c150.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c150 <- paramSweep(Navin_hbca_c150, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c150 <- summarizeSweep(sweep.res.Navin_hbca_c150, GT = FALSE)
bcmvn_Navin_hbca_c150 <- find.pK(sweep.stats.Navin_hbca_c150)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c150[which.max(bcmvn_Navin_hbca_c150$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c150@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c150@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c150 <- doubletFinder(Navin_hbca_c150, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c150_Quant <- (Navin_hbca_c150@meta.data$DF.classification == "Singlet")
Navin_hbca_c150_Quant_Singlets <- length(Navin_hbca_c150_Quant[Navin_hbca_c150_Quant== TRUE])
Navin_hbca_c150_Quant_Doublets <- length(Navin_hbca_c150_Quant[Navin_hbca_c150_Quant== FALSE])
Navin_hbca_c150_Quant_Doublets_Percent <- Navin_hbca_c150_Quant_Doublets / (Navin_hbca_c150_Quant_Doublets + Navin_hbca_c150_Quant_Singlets) * 100
Navin_hbca_c150_Quant <- as.data.frame(c(Navin_hbca_c150_Quant_Singlets, Navin_hbca_c150_Quant_Doublets, Navin_hbca_c150_Quant_Doublets_Percent))
colnames(Navin_hbca_c150_Quant) <- c("Navin_hbca_c150_Quant")
rownames(Navin_hbca_c150_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c150_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c150_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c150_Singlets <- subset(Navin_hbca_c150, cells=rownames(Navin_hbca_c150@meta.data)[which(Navin_hbca_c150@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c150_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c150_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c150)
rm(Navin_hbca_c150_Singlets)
rm(Navin_hbca_c150_Quant)
rm(Navin_hbca_c150_Quant_Singlets)
rm(Navin_hbca_c150_Quant_Doublets)
rm(Navin_hbca_c150_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c150)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c150)
rm(sweep.stats.Navin_hbca_c150)
gc()



################################################################################################
########################                Navin_hbca_c151                  #######################
################################################################################################


# Load Data
Navin_hbca_c151 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c151.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c151 <- paramSweep(Navin_hbca_c151, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c151 <- summarizeSweep(sweep.res.Navin_hbca_c151, GT = FALSE)
bcmvn_Navin_hbca_c151 <- find.pK(sweep.stats.Navin_hbca_c151)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c151[which.max(bcmvn_Navin_hbca_c151$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c151@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c151@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c151 <- doubletFinder(Navin_hbca_c151, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c151_Quant <- (Navin_hbca_c151@meta.data$DF.classification == "Singlet")
Navin_hbca_c151_Quant_Singlets <- length(Navin_hbca_c151_Quant[Navin_hbca_c151_Quant== TRUE])
Navin_hbca_c151_Quant_Doublets <- length(Navin_hbca_c151_Quant[Navin_hbca_c151_Quant== FALSE])
Navin_hbca_c151_Quant_Doublets_Percent <- Navin_hbca_c151_Quant_Doublets / (Navin_hbca_c151_Quant_Doublets + Navin_hbca_c151_Quant_Singlets) * 100
Navin_hbca_c151_Quant <- as.data.frame(c(Navin_hbca_c151_Quant_Singlets, Navin_hbca_c151_Quant_Doublets, Navin_hbca_c151_Quant_Doublets_Percent))
colnames(Navin_hbca_c151_Quant) <- c("Navin_hbca_c151_Quant")
rownames(Navin_hbca_c151_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c151_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c151_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c151_Singlets <- subset(Navin_hbca_c151, cells=rownames(Navin_hbca_c151@meta.data)[which(Navin_hbca_c151@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c151_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c151_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c151)
rm(Navin_hbca_c151_Singlets)
rm(Navin_hbca_c151_Quant)
rm(Navin_hbca_c151_Quant_Singlets)
rm(Navin_hbca_c151_Quant_Doublets)
rm(Navin_hbca_c151_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c151)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c151)
rm(sweep.stats.Navin_hbca_c151)
gc()



################################################################################################
########################                Navin_hbca_c152                  #######################
################################################################################################


# Load Data
Navin_hbca_c152 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c152.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c152 <- paramSweep(Navin_hbca_c152, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c152 <- summarizeSweep(sweep.res.Navin_hbca_c152, GT = FALSE)
bcmvn_Navin_hbca_c152 <- find.pK(sweep.stats.Navin_hbca_c152)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c152[which.max(bcmvn_Navin_hbca_c152$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c152@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c152@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c152 <- doubletFinder(Navin_hbca_c152, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c152_Quant <- (Navin_hbca_c152@meta.data$DF.classification == "Singlet")
Navin_hbca_c152_Quant_Singlets <- length(Navin_hbca_c152_Quant[Navin_hbca_c152_Quant== TRUE])
Navin_hbca_c152_Quant_Doublets <- length(Navin_hbca_c152_Quant[Navin_hbca_c152_Quant== FALSE])
Navin_hbca_c152_Quant_Doublets_Percent <- Navin_hbca_c152_Quant_Doublets / (Navin_hbca_c152_Quant_Doublets + Navin_hbca_c152_Quant_Singlets) * 100
Navin_hbca_c152_Quant <- as.data.frame(c(Navin_hbca_c152_Quant_Singlets, Navin_hbca_c152_Quant_Doublets, Navin_hbca_c152_Quant_Doublets_Percent))
colnames(Navin_hbca_c152_Quant) <- c("Navin_hbca_c152_Quant")
rownames(Navin_hbca_c152_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c152_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c152_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c152_Singlets <- subset(Navin_hbca_c152, cells=rownames(Navin_hbca_c152@meta.data)[which(Navin_hbca_c152@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c152_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c152_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c152)
rm(Navin_hbca_c152_Singlets)
rm(Navin_hbca_c152_Quant)
rm(Navin_hbca_c152_Quant_Singlets)
rm(Navin_hbca_c152_Quant_Doublets)
rm(Navin_hbca_c152_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c152)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c152)
rm(sweep.stats.Navin_hbca_c152)
gc()



################################################################################################
########################                Navin_hbca_c153                  #######################
################################################################################################


# Load Data
Navin_hbca_c153 <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total/Navin_hbca_c153.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Navin_hbca_c153 <- paramSweep(Navin_hbca_c153, PCs = 1:10, sct = FALSE)
sweep.stats.Navin_hbca_c153 <- summarizeSweep(sweep.res.Navin_hbca_c153, GT = FALSE)
bcmvn_Navin_hbca_c153 <- find.pK(sweep.stats.Navin_hbca_c153)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Navin_hbca_c153[which.max(bcmvn_Navin_hbca_c153$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Navin_hbca_c153@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 9.6% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Navin dataset = 10k per sample, will use doublet rate for 12k cells to err on the side of caution
nExp_poi <- round(0.096*nrow(Navin_hbca_c153@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Navin_hbca_c153 <- doubletFinder(Navin_hbca_c153, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Navin_hbca_c153_Quant <- (Navin_hbca_c153@meta.data$DF.classification == "Singlet")
Navin_hbca_c153_Quant_Singlets <- length(Navin_hbca_c153_Quant[Navin_hbca_c153_Quant== TRUE])
Navin_hbca_c153_Quant_Doublets <- length(Navin_hbca_c153_Quant[Navin_hbca_c153_Quant== FALSE])
Navin_hbca_c153_Quant_Doublets_Percent <- Navin_hbca_c153_Quant_Doublets / (Navin_hbca_c153_Quant_Doublets + Navin_hbca_c153_Quant_Singlets) * 100
Navin_hbca_c153_Quant <- as.data.frame(c(Navin_hbca_c153_Quant_Singlets, Navin_hbca_c153_Quant_Doublets, Navin_hbca_c153_Quant_Doublets_Percent))
colnames(Navin_hbca_c153_Quant) <- c("Navin_hbca_c153_Quant")
rownames(Navin_hbca_c153_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Navin_hbca_c153_Quant, file = "/R/R_Navin/Navin_Doublet_Tables/Navin_hbca_c153_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Navin_hbca_c153_Singlets <- subset(Navin_hbca_c153, cells=rownames(Navin_hbca_c153@meta.data)[which(Navin_hbca_c153@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Navin_hbca_c153_Singlets, file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c153_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Navin_hbca_c153)
rm(Navin_hbca_c153_Singlets)
rm(Navin_hbca_c153_Quant)
rm(Navin_hbca_c153_Quant_Singlets)
rm(Navin_hbca_c153_Quant_Doublets)
rm(Navin_hbca_c153_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Navin_hbca_c153)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Navin_hbca_c153)
rm(sweep.stats.Navin_hbca_c153)
gc()





