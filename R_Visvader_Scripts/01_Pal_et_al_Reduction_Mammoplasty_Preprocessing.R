#####################################################################################################################
#                Visvader (Pal et al) Dataset Reduction Mammoplasty Analysis Steps                                  #
#####################################################################################################################
# Step 1: Process Each Reduction Mammoplasty Count Matrix into Seurat Object                                        #
# Step 2: Run DoubletFinder and Subset Singlets                                                                     #
#####################################################################################################################


######################################################################
# Step 1: Process Each Breast Cancer Count Matrix into Seurat Object #
######################################################################

# Load Libraries 
library(Seurat)
library(patchwork)
library(dplyr)

# Data downloaded from GEO (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE161529)

#################################
# Reduction Mammoplasty Samples #
#################################

##################################################################
########################### Visvader_0019_NORM_Total #############
##################################################################

# Load the dataset
Visvader_0019_NORM_Total <- Read10X(data.dir = "/R/R_Visvader/Visvader_Input/NORM/Visvader_0019_NORM_Total/")

# Initialize the Seurat object with the raw (non-normalized data).
Visvader_0019_NORM_Total <- CreateSeuratObject(counts = Visvader_0019_NORM_Total, project = "Visvader_0019_NORM_Total", min.cells = 3, min.features = 200)
Visvader_0019_NORM_Total

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Visvader_0019_NORM_Total[["percent.mt"]] <- PercentageFeatureSet(Visvader_0019_NORM_Total, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Visvader_0019_NORM_Total, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Visvader_0019_NORM_Total, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Visvader_0019_NORM_Total, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Visvader_0019_NORM_Total <- subset(Visvader_0019_NORM_Total, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Visvader_0019_NORM_Total <- NormalizeData(Visvader_0019_NORM_Total, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Visvader_0019_NORM_Total <- FindVariableFeatures(Visvader_0019_NORM_Total, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Visvader_0019_NORM_Total)
Visvader_0019_NORM_Total <- ScaleData(Visvader_0019_NORM_Total, features = all.genes)

# Perform linear dimensional reduction
Visvader_0019_NORM_Total <- RunPCA(Visvader_0019_NORM_Total, features = VariableFeatures(object = Visvader_0019_NORM_Total))

# Examine and visualize PCA results a few different ways
print(Visvader_0019_NORM_Total[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Visvader_0019_NORM_Total, dims = 1:2, reduction = "pca")
DimPlot(Visvader_0019_NORM_Total, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Visvader_0019_NORM_Total <- JackStraw(Visvader_0019_NORM_Total, num.replicate = 100)
Visvader_0019_NORM_Total <- ScoreJackStraw(Visvader_0019_NORM_Total, dims = 1:20)
ElbowPlot(Visvader_0019_NORM_Total)

# Cluster cells
Visvader_0019_NORM_Total <- FindNeighbors(Visvader_0019_NORM_Total, dims = 1:20)
Visvader_0019_NORM_Total <- FindClusters(Visvader_0019_NORM_Total, resolution = 0.5)

# Run non-linear dimensional reduction
Visvader_0019_NORM_Total <- RunUMAP(Visvader_0019_NORM_Total, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Visvader_0019_NORM_Total, reduction = "umap")

# Save RDS
saveRDS(Visvader_0019_NORM_Total, file = "/R/R_Visvader/Visvader_RDS/RDS_Total/Visvader_0019_NORM_Total.rds")




##################################################################
########################### Visvader_0021_NORM_Total #############
##################################################################

# Load the dataset
Visvader_0021_NORM_Total <- Read10X(data.dir = "/R/R_Visvader/Visvader_Input/NORM/Visvader_0021_NORM_Total/")

# Initialize the Seurat object with the raw (non-normalized data).
Visvader_0021_NORM_Total <- CreateSeuratObject(counts = Visvader_0021_NORM_Total, project = "Visvader_0021_NORM_Total", min.cells = 3, min.features = 200)
Visvader_0021_NORM_Total

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Visvader_0021_NORM_Total[["percent.mt"]] <- PercentageFeatureSet(Visvader_0021_NORM_Total, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Visvader_0021_NORM_Total, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Visvader_0021_NORM_Total, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Visvader_0021_NORM_Total, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Visvader_0021_NORM_Total <- subset(Visvader_0021_NORM_Total, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Visvader_0021_NORM_Total <- NormalizeData(Visvader_0021_NORM_Total, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Visvader_0021_NORM_Total <- FindVariableFeatures(Visvader_0021_NORM_Total, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Visvader_0021_NORM_Total)
Visvader_0021_NORM_Total <- ScaleData(Visvader_0021_NORM_Total, features = all.genes)

# Perform linear dimensional reduction
Visvader_0021_NORM_Total <- RunPCA(Visvader_0021_NORM_Total, features = VariableFeatures(object = Visvader_0021_NORM_Total))

# Examine and visualize PCA results a few different ways
print(Visvader_0021_NORM_Total[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Visvader_0021_NORM_Total, dims = 1:2, reduction = "pca")
DimPlot(Visvader_0021_NORM_Total, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Visvader_0021_NORM_Total <- JackStraw(Visvader_0021_NORM_Total, num.replicate = 100)
Visvader_0021_NORM_Total <- ScoreJackStraw(Visvader_0021_NORM_Total, dims = 1:20)
ElbowPlot(Visvader_0021_NORM_Total)

# Cluster cells
Visvader_0021_NORM_Total <- FindNeighbors(Visvader_0021_NORM_Total, dims = 1:20)
Visvader_0021_NORM_Total <- FindClusters(Visvader_0021_NORM_Total, resolution = 0.5)

# Run non-linear dimensional reduction
Visvader_0021_NORM_Total <- RunUMAP(Visvader_0021_NORM_Total, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Visvader_0021_NORM_Total, reduction = "umap")

# Save RDS
saveRDS(Visvader_0021_NORM_Total, file = "/R/R_Visvader/Visvader_RDS/RDS_Total/Visvader_0021_NORM_Total.rds")




##################################################################
########################### Visvader_0064_NORM_Total #############
##################################################################

# Load the dataset
Visvader_0064_NORM_Total <- Read10X(data.dir = "/R/R_Visvader/Visvader_Input/NORM/Visvader_0064_NORM_Total/")

# Initialize the Seurat object with the raw (non-normalized data).
Visvader_0064_NORM_Total <- CreateSeuratObject(counts = Visvader_0064_NORM_Total, project = "Visvader_0064_NORM_Total", min.cells = 3, min.features = 200)
Visvader_0064_NORM_Total

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Visvader_0064_NORM_Total[["percent.mt"]] <- PercentageFeatureSet(Visvader_0064_NORM_Total, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Visvader_0064_NORM_Total, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Visvader_0064_NORM_Total, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Visvader_0064_NORM_Total, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Visvader_0064_NORM_Total <- subset(Visvader_0064_NORM_Total, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Visvader_0064_NORM_Total <- NormalizeData(Visvader_0064_NORM_Total, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Visvader_0064_NORM_Total <- FindVariableFeatures(Visvader_0064_NORM_Total, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Visvader_0064_NORM_Total)
Visvader_0064_NORM_Total <- ScaleData(Visvader_0064_NORM_Total, features = all.genes)

# Perform linear dimensional reduction
Visvader_0064_NORM_Total <- RunPCA(Visvader_0064_NORM_Total, features = VariableFeatures(object = Visvader_0064_NORM_Total))

# Examine and visualize PCA results a few different ways
print(Visvader_0064_NORM_Total[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Visvader_0064_NORM_Total, dims = 1:2, reduction = "pca")
DimPlot(Visvader_0064_NORM_Total, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Visvader_0064_NORM_Total <- JackStraw(Visvader_0064_NORM_Total, num.replicate = 100)
Visvader_0064_NORM_Total <- ScoreJackStraw(Visvader_0064_NORM_Total, dims = 1:20)
ElbowPlot(Visvader_0064_NORM_Total)

# Cluster cells
Visvader_0064_NORM_Total <- FindNeighbors(Visvader_0064_NORM_Total, dims = 1:20)
Visvader_0064_NORM_Total <- FindClusters(Visvader_0064_NORM_Total, resolution = 0.5)

# Run non-linear dimensional reduction
Visvader_0064_NORM_Total <- RunUMAP(Visvader_0064_NORM_Total, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Visvader_0064_NORM_Total, reduction = "umap")

# Save RDS
saveRDS(Visvader_0064_NORM_Total, file = "/R/R_Visvader/Visvader_RDS/RDS_Total/Visvader_0064_NORM_Total.rds")






##################################################################
########################### Visvader_0092_NORM_Total #############
##################################################################

# Load the dataset
Visvader_0092_NORM_Total <- Read10X(data.dir = "/R/R_Visvader/Visvader_Input/NORM/Visvader_0092_NORM_Total/")

# Initialize the Seurat object with the raw (non-normalized data).
Visvader_0092_NORM_Total <- CreateSeuratObject(counts = Visvader_0092_NORM_Total, project = "Visvader_0092_NORM_Total", min.cells = 3, min.features = 200)
Visvader_0092_NORM_Total

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Visvader_0092_NORM_Total[["percent.mt"]] <- PercentageFeatureSet(Visvader_0092_NORM_Total, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Visvader_0092_NORM_Total, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Visvader_0092_NORM_Total, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Visvader_0092_NORM_Total, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Visvader_0092_NORM_Total <- subset(Visvader_0092_NORM_Total, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Visvader_0092_NORM_Total <- NormalizeData(Visvader_0092_NORM_Total, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Visvader_0092_NORM_Total <- FindVariableFeatures(Visvader_0092_NORM_Total, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Visvader_0092_NORM_Total)
Visvader_0092_NORM_Total <- ScaleData(Visvader_0092_NORM_Total, features = all.genes)

# Perform linear dimensional reduction
Visvader_0092_NORM_Total <- RunPCA(Visvader_0092_NORM_Total, features = VariableFeatures(object = Visvader_0092_NORM_Total))

# Examine and visualize PCA results a few different ways
print(Visvader_0092_NORM_Total[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Visvader_0092_NORM_Total, dims = 1:2, reduction = "pca")
DimPlot(Visvader_0092_NORM_Total, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Visvader_0092_NORM_Total <- JackStraw(Visvader_0092_NORM_Total, num.replicate = 100)
Visvader_0092_NORM_Total <- ScoreJackStraw(Visvader_0092_NORM_Total, dims = 1:20)
ElbowPlot(Visvader_0092_NORM_Total)

# Cluster cells
Visvader_0092_NORM_Total <- FindNeighbors(Visvader_0092_NORM_Total, dims = 1:20)
Visvader_0092_NORM_Total <- FindClusters(Visvader_0092_NORM_Total, resolution = 0.5)

# Run non-linear dimensional reduction
Visvader_0092_NORM_Total <- RunUMAP(Visvader_0092_NORM_Total, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Visvader_0092_NORM_Total, reduction = "umap")

# Save RDS
saveRDS(Visvader_0092_NORM_Total, file = "/R/R_Visvader/Visvader_RDS/RDS_Total/Visvader_0092_NORM_Total.rds")







##################################################################
########################### Visvader_0093_NORM_Total #############
##################################################################

# Load the dataset
Visvader_0093_NORM_Total <- Read10X(data.dir = "/R/R_Visvader/Visvader_Input/NORM/Visvader_0093_NORM_Total/")

# Initialize the Seurat object with the raw (non-normalized data).
Visvader_0093_NORM_Total <- CreateSeuratObject(counts = Visvader_0093_NORM_Total, project = "Visvader_0093_NORM_Total", min.cells = 3, min.features = 200)
Visvader_0093_NORM_Total

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Visvader_0093_NORM_Total[["percent.mt"]] <- PercentageFeatureSet(Visvader_0093_NORM_Total, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Visvader_0093_NORM_Total, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Visvader_0093_NORM_Total, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Visvader_0093_NORM_Total, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Visvader_0093_NORM_Total <- subset(Visvader_0093_NORM_Total, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Visvader_0093_NORM_Total <- NormalizeData(Visvader_0093_NORM_Total, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Visvader_0093_NORM_Total <- FindVariableFeatures(Visvader_0093_NORM_Total, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Visvader_0093_NORM_Total)
Visvader_0093_NORM_Total <- ScaleData(Visvader_0093_NORM_Total, features = all.genes)

# Perform linear dimensional reduction
Visvader_0093_NORM_Total <- RunPCA(Visvader_0093_NORM_Total, features = VariableFeatures(object = Visvader_0093_NORM_Total))

# Examine and visualize PCA results a few different ways
print(Visvader_0093_NORM_Total[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Visvader_0093_NORM_Total, dims = 1:2, reduction = "pca")
DimPlot(Visvader_0093_NORM_Total, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Visvader_0093_NORM_Total <- JackStraw(Visvader_0093_NORM_Total, num.replicate = 100)
Visvader_0093_NORM_Total <- ScoreJackStraw(Visvader_0093_NORM_Total, dims = 1:20)
ElbowPlot(Visvader_0093_NORM_Total)

# Cluster cells
Visvader_0093_NORM_Total <- FindNeighbors(Visvader_0093_NORM_Total, dims = 1:20)
Visvader_0093_NORM_Total <- FindClusters(Visvader_0093_NORM_Total, resolution = 0.5)

# Run non-linear dimensional reduction
Visvader_0093_NORM_Total <- RunUMAP(Visvader_0093_NORM_Total, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Visvader_0093_NORM_Total, reduction = "umap")

# Save RDS
saveRDS(Visvader_0093_NORM_Total, file = "/R/R_Visvader/Visvader_RDS/RDS_Total/Visvader_0093_NORM_Total.rds")







##################################################################
########################### Visvader_0123_NORM_Total #############
##################################################################

# Load the dataset
Visvader_0123_NORM_Total <- Read10X(data.dir = "/R/R_Visvader/Visvader_Input/NORM/Visvader_0123_NORM_Total/")

# Initialize the Seurat object with the raw (non-normalized data).
Visvader_0123_NORM_Total <- CreateSeuratObject(counts = Visvader_0123_NORM_Total, project = "Visvader_0123_NORM_Total", min.cells = 3, min.features = 200)
Visvader_0123_NORM_Total

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Visvader_0123_NORM_Total[["percent.mt"]] <- PercentageFeatureSet(Visvader_0123_NORM_Total, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Visvader_0123_NORM_Total, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Visvader_0123_NORM_Total, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Visvader_0123_NORM_Total, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Visvader_0123_NORM_Total <- subset(Visvader_0123_NORM_Total, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Visvader_0123_NORM_Total <- NormalizeData(Visvader_0123_NORM_Total, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Visvader_0123_NORM_Total <- FindVariableFeatures(Visvader_0123_NORM_Total, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Visvader_0123_NORM_Total)
Visvader_0123_NORM_Total <- ScaleData(Visvader_0123_NORM_Total, features = all.genes)

# Perform linear dimensional reduction
Visvader_0123_NORM_Total <- RunPCA(Visvader_0123_NORM_Total, features = VariableFeatures(object = Visvader_0123_NORM_Total))

# Examine and visualize PCA results a few different ways
print(Visvader_0123_NORM_Total[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Visvader_0123_NORM_Total, dims = 1:2, reduction = "pca")
DimPlot(Visvader_0123_NORM_Total, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Visvader_0123_NORM_Total <- JackStraw(Visvader_0123_NORM_Total, num.replicate = 100)
Visvader_0123_NORM_Total <- ScoreJackStraw(Visvader_0123_NORM_Total, dims = 1:20)
ElbowPlot(Visvader_0123_NORM_Total)

# Cluster cells
Visvader_0123_NORM_Total <- FindNeighbors(Visvader_0123_NORM_Total, dims = 1:20)
Visvader_0123_NORM_Total <- FindClusters(Visvader_0123_NORM_Total, resolution = 0.5)

# Run non-linear dimensional reduction
Visvader_0123_NORM_Total <- RunUMAP(Visvader_0123_NORM_Total, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Visvader_0123_NORM_Total, reduction = "umap")

# Save RDS
saveRDS(Visvader_0123_NORM_Total, file = "/R/R_Visvader/Visvader_RDS/RDS_Total/Visvader_0123_NORM_Total.rds")







##################################################################
########################### Visvader_0169_NORM_Total #############
##################################################################

# Load the dataset
Visvader_0169_NORM_Total <- Read10X(data.dir = "/R/R_Visvader/Visvader_Input/NORM/Visvader_0169_NORM_Total/")

# Initialize the Seurat object with the raw (non-normalized data).
Visvader_0169_NORM_Total <- CreateSeuratObject(counts = Visvader_0169_NORM_Total, project = "Visvader_0169_NORM_Total", min.cells = 3, min.features = 200)
Visvader_0169_NORM_Total

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Visvader_0169_NORM_Total[["percent.mt"]] <- PercentageFeatureSet(Visvader_0169_NORM_Total, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Visvader_0169_NORM_Total, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Visvader_0169_NORM_Total, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Visvader_0169_NORM_Total, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Visvader_0169_NORM_Total <- subset(Visvader_0169_NORM_Total, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Visvader_0169_NORM_Total <- NormalizeData(Visvader_0169_NORM_Total, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Visvader_0169_NORM_Total <- FindVariableFeatures(Visvader_0169_NORM_Total, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Visvader_0169_NORM_Total)
Visvader_0169_NORM_Total <- ScaleData(Visvader_0169_NORM_Total, features = all.genes)

# Perform linear dimensional reduction
Visvader_0169_NORM_Total <- RunPCA(Visvader_0169_NORM_Total, features = VariableFeatures(object = Visvader_0169_NORM_Total))

# Examine and visualize PCA results a few different ways
print(Visvader_0169_NORM_Total[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Visvader_0169_NORM_Total, dims = 1:2, reduction = "pca")
DimPlot(Visvader_0169_NORM_Total, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Visvader_0169_NORM_Total <- JackStraw(Visvader_0169_NORM_Total, num.replicate = 100)
Visvader_0169_NORM_Total <- ScoreJackStraw(Visvader_0169_NORM_Total, dims = 1:20)
ElbowPlot(Visvader_0169_NORM_Total)

# Cluster cells
Visvader_0169_NORM_Total <- FindNeighbors(Visvader_0169_NORM_Total, dims = 1:20)
Visvader_0169_NORM_Total <- FindClusters(Visvader_0169_NORM_Total, resolution = 0.5)

# Run non-linear dimensional reduction
Visvader_0169_NORM_Total <- RunUMAP(Visvader_0169_NORM_Total, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Visvader_0169_NORM_Total, reduction = "umap")

# Save RDS
saveRDS(Visvader_0169_NORM_Total, file = "/R/R_Visvader/Visvader_RDS/RDS_Total/Visvader_0169_NORM_Total.rds")







##################################################################
########################### Visvader_0230_NORM_Total #############
##################################################################

# Load the dataset
Visvader_0230_NORM_Total <- Read10X(data.dir = "/R/R_Visvader/Visvader_Input/NORM/Visvader_0230_NORM_Total/")

# Initialize the Seurat object with the raw (non-normalized data).
Visvader_0230_NORM_Total <- CreateSeuratObject(counts = Visvader_0230_NORM_Total, project = "Visvader_0230_NORM_Total", min.cells = 3, min.features = 200)
Visvader_0230_NORM_Total

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Visvader_0230_NORM_Total[["percent.mt"]] <- PercentageFeatureSet(Visvader_0230_NORM_Total, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Visvader_0230_NORM_Total, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Visvader_0230_NORM_Total, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Visvader_0230_NORM_Total, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Visvader_0230_NORM_Total <- subset(Visvader_0230_NORM_Total, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Visvader_0230_NORM_Total <- NormalizeData(Visvader_0230_NORM_Total, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Visvader_0230_NORM_Total <- FindVariableFeatures(Visvader_0230_NORM_Total, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Visvader_0230_NORM_Total)
Visvader_0230_NORM_Total <- ScaleData(Visvader_0230_NORM_Total, features = all.genes)

# Perform linear dimensional reduction
Visvader_0230_NORM_Total <- RunPCA(Visvader_0230_NORM_Total, features = VariableFeatures(object = Visvader_0230_NORM_Total))

# Examine and visualize PCA results a few different ways
print(Visvader_0230_NORM_Total[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Visvader_0230_NORM_Total, dims = 1:2, reduction = "pca")
DimPlot(Visvader_0230_NORM_Total, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Visvader_0230_NORM_Total <- JackStraw(Visvader_0230_NORM_Total, num.replicate = 100)
Visvader_0230_NORM_Total <- ScoreJackStraw(Visvader_0230_NORM_Total, dims = 1:20)
ElbowPlot(Visvader_0230_NORM_Total)

# Cluster cells
Visvader_0230_NORM_Total <- FindNeighbors(Visvader_0230_NORM_Total, dims = 1:20)
Visvader_0230_NORM_Total <- FindClusters(Visvader_0230_NORM_Total, resolution = 0.5)

# Run non-linear dimensional reduction
Visvader_0230_NORM_Total <- RunUMAP(Visvader_0230_NORM_Total, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Visvader_0230_NORM_Total, reduction = "umap")

# Save RDS
saveRDS(Visvader_0230_NORM_Total, file = "/R/R_Visvader/Visvader_RDS/RDS_Total/Visvader_0230_NORM_Total.rds")







##################################################################
########################### Visvader_0233_NORM_Total #############
##################################################################

# Load the dataset
Visvader_0233_NORM_Total <- Read10X(data.dir = "/R/R_Visvader/Visvader_Input/NORM/Visvader_0233_NORM_Total/")

# Initialize the Seurat object with the raw (non-normalized data).
Visvader_0233_NORM_Total <- CreateSeuratObject(counts = Visvader_0233_NORM_Total, project = "Visvader_0233_NORM_Total", min.cells = 3, min.features = 200)
Visvader_0233_NORM_Total

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Visvader_0233_NORM_Total[["percent.mt"]] <- PercentageFeatureSet(Visvader_0233_NORM_Total, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Visvader_0233_NORM_Total, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Visvader_0233_NORM_Total, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Visvader_0233_NORM_Total, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Visvader_0233_NORM_Total <- subset(Visvader_0233_NORM_Total, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Visvader_0233_NORM_Total <- NormalizeData(Visvader_0233_NORM_Total, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Visvader_0233_NORM_Total <- FindVariableFeatures(Visvader_0233_NORM_Total, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Visvader_0233_NORM_Total)
Visvader_0233_NORM_Total <- ScaleData(Visvader_0233_NORM_Total, features = all.genes)

# Perform linear dimensional reduction
Visvader_0233_NORM_Total <- RunPCA(Visvader_0233_NORM_Total, features = VariableFeatures(object = Visvader_0233_NORM_Total))

# Examine and visualize PCA results a few different ways
print(Visvader_0233_NORM_Total[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Visvader_0233_NORM_Total, dims = 1:2, reduction = "pca")
DimPlot(Visvader_0233_NORM_Total, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Visvader_0233_NORM_Total <- JackStraw(Visvader_0233_NORM_Total, num.replicate = 100)
Visvader_0233_NORM_Total <- ScoreJackStraw(Visvader_0233_NORM_Total, dims = 1:20)
ElbowPlot(Visvader_0233_NORM_Total)

# Cluster cells
Visvader_0233_NORM_Total <- FindNeighbors(Visvader_0233_NORM_Total, dims = 1:20)
Visvader_0233_NORM_Total <- FindClusters(Visvader_0233_NORM_Total, resolution = 0.5)

# Run non-linear dimensional reduction
Visvader_0233_NORM_Total <- RunUMAP(Visvader_0233_NORM_Total, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Visvader_0233_NORM_Total, reduction = "umap")

# Save RDS
saveRDS(Visvader_0233_NORM_Total, file = "/R/R_Visvader/Visvader_RDS/RDS_Total/Visvader_0233_NORM_Total.rds")







##################################################################
########################### Visvader_0275_NORM_Total #############
##################################################################

# Load the dataset
Visvader_0275_NORM_Total <- Read10X(data.dir = "/R/R_Visvader/Visvader_Input/NORM/Visvader_0275_NORM_Total/")

# Initialize the Seurat object with the raw (non-normalized data).
Visvader_0275_NORM_Total <- CreateSeuratObject(counts = Visvader_0275_NORM_Total, project = "Visvader_0275_NORM_Total", min.cells = 3, min.features = 200)
Visvader_0275_NORM_Total

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Visvader_0275_NORM_Total[["percent.mt"]] <- PercentageFeatureSet(Visvader_0275_NORM_Total, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Visvader_0275_NORM_Total, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Visvader_0275_NORM_Total, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Visvader_0275_NORM_Total, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Visvader_0275_NORM_Total <- subset(Visvader_0275_NORM_Total, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Visvader_0275_NORM_Total <- NormalizeData(Visvader_0275_NORM_Total, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Visvader_0275_NORM_Total <- FindVariableFeatures(Visvader_0275_NORM_Total, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Visvader_0275_NORM_Total)
Visvader_0275_NORM_Total <- ScaleData(Visvader_0275_NORM_Total, features = all.genes)

# Perform linear dimensional reduction
Visvader_0275_NORM_Total <- RunPCA(Visvader_0275_NORM_Total, features = VariableFeatures(object = Visvader_0275_NORM_Total))

# Examine and visualize PCA results a few different ways
print(Visvader_0275_NORM_Total[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Visvader_0275_NORM_Total, dims = 1:2, reduction = "pca")
DimPlot(Visvader_0275_NORM_Total, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Visvader_0275_NORM_Total <- JackStraw(Visvader_0275_NORM_Total, num.replicate = 100)
Visvader_0275_NORM_Total <- ScoreJackStraw(Visvader_0275_NORM_Total, dims = 1:20)
ElbowPlot(Visvader_0275_NORM_Total)

# Cluster cells
Visvader_0275_NORM_Total <- FindNeighbors(Visvader_0275_NORM_Total, dims = 1:20)
Visvader_0275_NORM_Total <- FindClusters(Visvader_0275_NORM_Total, resolution = 0.5)

# Run non-linear dimensional reduction
Visvader_0275_NORM_Total <- RunUMAP(Visvader_0275_NORM_Total, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Visvader_0275_NORM_Total, reduction = "umap")

# Save RDS
saveRDS(Visvader_0275_NORM_Total, file = "/R/R_Visvader/Visvader_RDS/RDS_Total/Visvader_0275_NORM_Total.rds")







##################################################################
########################### Visvader_0288_NORM_Total #############
##################################################################

# Load the dataset
Visvader_0288_NORM_Total <- Read10X(data.dir = "/R/R_Visvader/Visvader_Input/NORM/Visvader_0288_NORM_Total/")

# Initialize the Seurat object with the raw (non-normalized data).
Visvader_0288_NORM_Total <- CreateSeuratObject(counts = Visvader_0288_NORM_Total, project = "Visvader_0288_NORM_Total", min.cells = 3, min.features = 200)
Visvader_0288_NORM_Total

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Visvader_0288_NORM_Total[["percent.mt"]] <- PercentageFeatureSet(Visvader_0288_NORM_Total, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Visvader_0288_NORM_Total, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Visvader_0288_NORM_Total, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Visvader_0288_NORM_Total, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Visvader_0288_NORM_Total <- subset(Visvader_0288_NORM_Total, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Visvader_0288_NORM_Total <- NormalizeData(Visvader_0288_NORM_Total, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Visvader_0288_NORM_Total <- FindVariableFeatures(Visvader_0288_NORM_Total, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Visvader_0288_NORM_Total)
Visvader_0288_NORM_Total <- ScaleData(Visvader_0288_NORM_Total, features = all.genes)

# Perform linear dimensional reduction
Visvader_0288_NORM_Total <- RunPCA(Visvader_0288_NORM_Total, features = VariableFeatures(object = Visvader_0288_NORM_Total))

# Examine and visualize PCA results a few different ways
print(Visvader_0288_NORM_Total[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Visvader_0288_NORM_Total, dims = 1:2, reduction = "pca")
DimPlot(Visvader_0288_NORM_Total, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Visvader_0288_NORM_Total <- JackStraw(Visvader_0288_NORM_Total, num.replicate = 100)
Visvader_0288_NORM_Total <- ScoreJackStraw(Visvader_0288_NORM_Total, dims = 1:20)
ElbowPlot(Visvader_0288_NORM_Total)

# Cluster cells
Visvader_0288_NORM_Total <- FindNeighbors(Visvader_0288_NORM_Total, dims = 1:20)
Visvader_0288_NORM_Total <- FindClusters(Visvader_0288_NORM_Total, resolution = 0.5)

# Run non-linear dimensional reduction
Visvader_0288_NORM_Total <- RunUMAP(Visvader_0288_NORM_Total, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Visvader_0288_NORM_Total, reduction = "umap")

# Save RDS
saveRDS(Visvader_0288_NORM_Total, file = "/R/R_Visvader/Visvader_RDS/RDS_Total/Visvader_0288_NORM_Total.rds")







##################################################################
########################### Visvader_0342_NORM_Total #############
##################################################################

# Load the dataset
Visvader_0342_NORM_Total <- Read10X(data.dir = "/R/R_Visvader/Visvader_Input/NORM/Visvader_0342_NORM_Total/")

# Initialize the Seurat object with the raw (non-normalized data).
Visvader_0342_NORM_Total <- CreateSeuratObject(counts = Visvader_0342_NORM_Total, project = "Visvader_0342_NORM_Total", min.cells = 3, min.features = 200)
Visvader_0342_NORM_Total

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Visvader_0342_NORM_Total[["percent.mt"]] <- PercentageFeatureSet(Visvader_0342_NORM_Total, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Visvader_0342_NORM_Total, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Visvader_0342_NORM_Total, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Visvader_0342_NORM_Total, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Visvader_0342_NORM_Total <- subset(Visvader_0342_NORM_Total, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Visvader_0342_NORM_Total <- NormalizeData(Visvader_0342_NORM_Total, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Visvader_0342_NORM_Total <- FindVariableFeatures(Visvader_0342_NORM_Total, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Visvader_0342_NORM_Total)
Visvader_0342_NORM_Total <- ScaleData(Visvader_0342_NORM_Total, features = all.genes)

# Perform linear dimensional reduction
Visvader_0342_NORM_Total <- RunPCA(Visvader_0342_NORM_Total, features = VariableFeatures(object = Visvader_0342_NORM_Total))

# Examine and visualize PCA results a few different ways
print(Visvader_0342_NORM_Total[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Visvader_0342_NORM_Total, dims = 1:2, reduction = "pca")
DimPlot(Visvader_0342_NORM_Total, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Visvader_0342_NORM_Total <- JackStraw(Visvader_0342_NORM_Total, num.replicate = 100)
Visvader_0342_NORM_Total <- ScoreJackStraw(Visvader_0342_NORM_Total, dims = 1:20)
ElbowPlot(Visvader_0342_NORM_Total)

# Cluster cells
Visvader_0342_NORM_Total <- FindNeighbors(Visvader_0342_NORM_Total, dims = 1:20)
Visvader_0342_NORM_Total <- FindClusters(Visvader_0342_NORM_Total, resolution = 0.5)

# Run non-linear dimensional reduction
Visvader_0342_NORM_Total <- RunUMAP(Visvader_0342_NORM_Total, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Visvader_0342_NORM_Total, reduction = "umap")

# Save RDS
saveRDS(Visvader_0342_NORM_Total, file = "/R/R_Visvader/Visvader_RDS/RDS_Total/Visvader_0342_NORM_Total.rds")







##################################################################
########################### Visvader_0372_NORM_Total #############
##################################################################

# Load the dataset
Visvader_0372_NORM_Total <- Read10X(data.dir = "/R/R_Visvader/Visvader_Input/NORM/Visvader_0372_NORM_Total/")

# Initialize the Seurat object with the raw (non-normalized data).
Visvader_0372_NORM_Total <- CreateSeuratObject(counts = Visvader_0372_NORM_Total, project = "Visvader_0372_NORM_Total", min.cells = 3, min.features = 200)
Visvader_0372_NORM_Total

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Visvader_0372_NORM_Total[["percent.mt"]] <- PercentageFeatureSet(Visvader_0372_NORM_Total, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Visvader_0372_NORM_Total, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Visvader_0372_NORM_Total, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Visvader_0372_NORM_Total, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Visvader_0372_NORM_Total <- subset(Visvader_0372_NORM_Total, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Visvader_0372_NORM_Total <- NormalizeData(Visvader_0372_NORM_Total, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Visvader_0372_NORM_Total <- FindVariableFeatures(Visvader_0372_NORM_Total, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Visvader_0372_NORM_Total)
Visvader_0372_NORM_Total <- ScaleData(Visvader_0372_NORM_Total, features = all.genes)

# Perform linear dimensional reduction
Visvader_0372_NORM_Total <- RunPCA(Visvader_0372_NORM_Total, features = VariableFeatures(object = Visvader_0372_NORM_Total))

# Examine and visualize PCA results a few different ways
print(Visvader_0372_NORM_Total[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Visvader_0372_NORM_Total, dims = 1:2, reduction = "pca")
DimPlot(Visvader_0372_NORM_Total, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Visvader_0372_NORM_Total <- JackStraw(Visvader_0372_NORM_Total, num.replicate = 100)
Visvader_0372_NORM_Total <- ScoreJackStraw(Visvader_0372_NORM_Total, dims = 1:20)
ElbowPlot(Visvader_0372_NORM_Total)

# Cluster cells
Visvader_0372_NORM_Total <- FindNeighbors(Visvader_0372_NORM_Total, dims = 1:20)
Visvader_0372_NORM_Total <- FindClusters(Visvader_0372_NORM_Total, resolution = 0.5)

# Run non-linear dimensional reduction
Visvader_0372_NORM_Total <- RunUMAP(Visvader_0372_NORM_Total, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Visvader_0372_NORM_Total, reduction = "umap")

# Save RDS
saveRDS(Visvader_0372_NORM_Total, file = "/R/R_Visvader/Visvader_RDS/RDS_Total/Visvader_0372_NORM_Total.rds")



#####################################################################
########### Step 2: Run DoubletFinder and Subset Singlets ########### 
#####################################################################

# Load Libraries
library(Seurat)
library(DoubletFinder)

# Assuming 5% doublet rate per sample based on average number of sequenced cells
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-

# Patient Samples
############################  Reduction Mammoplasties
# Visvader_0019_NORM_Total
# Visvader_0021_NORM_Total
# Visvader_0064_NORM_Total
# Visvader_0092_NORM_Total
# Visvader_0093_NORM_Total
# Visvader_0123_NORM_Total
# Visvader_0169_NORM_Total
# Visvader_0230_NORM_Total
# Visvader_0233_NORM_Total
# Visvader_0275_NORM_Total
# Visvader_0288_NORM_Total
# Visvader_0342_NORM_Total
# Visvader_0372_NORM_Total



#######################################################################################################
########################                Reduction Mammoplasty                  ########################
#######################################################################################################

################################################################################################
########################              Visvader_0019_NORM_Total             #####################
################################################################################################

# Load Data
Visvader_0019_NORM_Total <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Total/Visvader_0019_NORM_Total.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Visvader_0019_NORM_Total <- paramSweep(Visvader_0019_NORM_Total, PCs = 1:10, sct = FALSE)
sweep.stats.Visvader_0019_NORM_Total <- summarizeSweep(sweep.res.Visvader_0019_NORM_Total, GT = FALSE)
bcmvn_Visvader_0019_NORM_Total <- find.pK(sweep.stats.Visvader_0019_NORM_Total)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Visvader_0019_NORM_Total[which.max(bcmvn_Visvader_0019_NORM_Total$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Visvader_0019_NORM_Total@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 5% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Visvader dataset = 5k per sample, will use doublet rate for 6,000 cells to err on the side of caution
nExp_poi <- round(0.05*nrow(Visvader_0019_NORM_Total@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Visvader_0019_NORM_Total <- doubletFinder(Visvader_0019_NORM_Total, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Visvader_0019_NORM_Total_Quant <- (Visvader_0019_NORM_Total@meta.data$DF.classification == "Singlet")
Visvader_0019_NORM_Total_Quant_Singlets <- length(Visvader_0019_NORM_Total_Quant[Visvader_0019_NORM_Total_Quant== TRUE])
Visvader_0019_NORM_Total_Quant_Doublets <- length(Visvader_0019_NORM_Total_Quant[Visvader_0019_NORM_Total_Quant== FALSE])
Visvader_0019_NORM_Total_Quant_Doublets_Percent <- Visvader_0019_NORM_Total_Quant_Doublets /  (Visvader_0019_NORM_Total_Quant_Doublets + Visvader_0019_NORM_Total_Quant_Singlets) * 100
Visvader_0019_NORM_Total_Quant <- as.data.frame(c(Visvader_0019_NORM_Total_Quant_Singlets, Visvader_0019_NORM_Total_Quant_Doublets, Visvader_0019_NORM_Total_Quant_Doublets_Percent))
colnames(Visvader_0019_NORM_Total_Quant) <- c("Visvader_0019_NORM_Total_Quant")
rownames(Visvader_0019_NORM_Total_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Visvader_0019_NORM_Total_Quant, file = "/R/R_Visvader/Visvader_Doublet_Tables/Visvader_0019_NORM_Total_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Visvader_0019_NORM_Total_Singlets <- subset(Visvader_0019_NORM_Total, cells=rownames(Visvader_0019_NORM_Total@meta.data)[which(Visvader_0019_NORM_Total@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Visvader_0019_NORM_Total_Singlets, file = "/R/R_Visvader/Visvader_RDS/RDS_Total_Singlets/Visvader_0019_NORM_Total_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Visvader_0019_NORM_Total)
rm(Visvader_0019_NORM_Total_Singlets)
rm(Visvader_0019_NORM_Total_Quant)
rm(Visvader_0019_NORM_Total_Quant_Singlets)
rm(Visvader_0019_NORM_Total_Quant_Doublets)
rm(Visvader_0019_NORM_Total_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Visvader_0019_NORM_Total)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Visvader_0019_NORM_Total)
rm(sweep.stats.Visvader_0019_NORM_Total)
gc()



################################################################################################
########################              Visvader_0021_NORM_Total             #####################
################################################################################################

# Load Data
Visvader_0021_NORM_Total <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Total/Visvader_0021_NORM_Total.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Visvader_0021_NORM_Total <- paramSweep(Visvader_0021_NORM_Total, PCs = 1:10, sct = FALSE)
sweep.stats.Visvader_0021_NORM_Total <- summarizeSweep(sweep.res.Visvader_0021_NORM_Total, GT = FALSE)
bcmvn_Visvader_0021_NORM_Total <- find.pK(sweep.stats.Visvader_0021_NORM_Total)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Visvader_0021_NORM_Total[which.max(bcmvn_Visvader_0021_NORM_Total$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Visvader_0021_NORM_Total@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 5% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Visvader dataset = 5k per sample, will use doublet rate for 6,000 cells to err on the side of caution
nExp_poi <- round(0.05*nrow(Visvader_0021_NORM_Total@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Visvader_0021_NORM_Total <- doubletFinder(Visvader_0021_NORM_Total, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Visvader_0021_NORM_Total_Quant <- (Visvader_0021_NORM_Total@meta.data$DF.classification == "Singlet")
Visvader_0021_NORM_Total_Quant_Singlets <- length(Visvader_0021_NORM_Total_Quant[Visvader_0021_NORM_Total_Quant== TRUE])
Visvader_0021_NORM_Total_Quant_Doublets <- length(Visvader_0021_NORM_Total_Quant[Visvader_0021_NORM_Total_Quant== FALSE])
Visvader_0021_NORM_Total_Quant_Doublets_Percent <- Visvader_0021_NORM_Total_Quant_Doublets /  (Visvader_0021_NORM_Total_Quant_Doublets + Visvader_0021_NORM_Total_Quant_Singlets) * 100
Visvader_0021_NORM_Total_Quant <- as.data.frame(c(Visvader_0021_NORM_Total_Quant_Singlets, Visvader_0021_NORM_Total_Quant_Doublets, Visvader_0021_NORM_Total_Quant_Doublets_Percent))
colnames(Visvader_0021_NORM_Total_Quant) <- c("Visvader_0021_NORM_Total_Quant")
rownames(Visvader_0021_NORM_Total_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Visvader_0021_NORM_Total_Quant, file = "/R/R_Visvader/Visvader_Doublet_Tables/Visvader_0021_NORM_Total_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Visvader_0021_NORM_Total_Singlets <- subset(Visvader_0021_NORM_Total, cells=rownames(Visvader_0021_NORM_Total@meta.data)[which(Visvader_0021_NORM_Total@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Visvader_0021_NORM_Total_Singlets, file = "/R/R_Visvader/Visvader_RDS/RDS_Total_Singlets/Visvader_0021_NORM_Total_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Visvader_0021_NORM_Total)
rm(Visvader_0021_NORM_Total_Singlets)
rm(Visvader_0021_NORM_Total_Quant)
rm(Visvader_0021_NORM_Total_Quant_Singlets)
rm(Visvader_0021_NORM_Total_Quant_Doublets)
rm(Visvader_0021_NORM_Total_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Visvader_0021_NORM_Total)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Visvader_0021_NORM_Total)
rm(sweep.stats.Visvader_0021_NORM_Total)
gc()




################################################################################################
########################              Visvader_0064_NORM_Total             #####################
################################################################################################

# Load Data
Visvader_0064_NORM_Total <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Total/Visvader_0064_NORM_Total.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Visvader_0064_NORM_Total <- paramSweep(Visvader_0064_NORM_Total, PCs = 1:10, sct = FALSE)
sweep.stats.Visvader_0064_NORM_Total <- summarizeSweep(sweep.res.Visvader_0064_NORM_Total, GT = FALSE)
bcmvn_Visvader_0064_NORM_Total <- find.pK(sweep.stats.Visvader_0064_NORM_Total)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Visvader_0064_NORM_Total[which.max(bcmvn_Visvader_0064_NORM_Total$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Visvader_0064_NORM_Total@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 5% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Visvader dataset = 5k per sample, will use doublet rate for 6,000 cells to err on the side of caution
nExp_poi <- round(0.05*nrow(Visvader_0064_NORM_Total@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Visvader_0064_NORM_Total <- doubletFinder(Visvader_0064_NORM_Total, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Visvader_0064_NORM_Total_Quant <- (Visvader_0064_NORM_Total@meta.data$DF.classification == "Singlet")
Visvader_0064_NORM_Total_Quant_Singlets <- length(Visvader_0064_NORM_Total_Quant[Visvader_0064_NORM_Total_Quant== TRUE])
Visvader_0064_NORM_Total_Quant_Doublets <- length(Visvader_0064_NORM_Total_Quant[Visvader_0064_NORM_Total_Quant== FALSE])
Visvader_0064_NORM_Total_Quant_Doublets_Percent <- Visvader_0064_NORM_Total_Quant_Doublets /  (Visvader_0064_NORM_Total_Quant_Doublets + Visvader_0064_NORM_Total_Quant_Singlets) * 100
Visvader_0064_NORM_Total_Quant <- as.data.frame(c(Visvader_0064_NORM_Total_Quant_Singlets, Visvader_0064_NORM_Total_Quant_Doublets, Visvader_0064_NORM_Total_Quant_Doublets_Percent))
colnames(Visvader_0064_NORM_Total_Quant) <- c("Visvader_0064_NORM_Total_Quant")
rownames(Visvader_0064_NORM_Total_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Visvader_0064_NORM_Total_Quant, file = "/R/R_Visvader/Visvader_Doublet_Tables/Visvader_0064_NORM_Total_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Visvader_0064_NORM_Total_Singlets <- subset(Visvader_0064_NORM_Total, cells=rownames(Visvader_0064_NORM_Total@meta.data)[which(Visvader_0064_NORM_Total@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Visvader_0064_NORM_Total_Singlets, file = "/R/R_Visvader/Visvader_RDS/RDS_Total_Singlets/Visvader_0064_NORM_Total_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Visvader_0064_NORM_Total)
rm(Visvader_0064_NORM_Total_Singlets)
rm(Visvader_0064_NORM_Total_Quant)
rm(Visvader_0064_NORM_Total_Quant_Singlets)
rm(Visvader_0064_NORM_Total_Quant_Doublets)
rm(Visvader_0064_NORM_Total_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Visvader_0064_NORM_Total)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Visvader_0064_NORM_Total)
rm(sweep.stats.Visvader_0064_NORM_Total)
gc()




################################################################################################
########################              Visvader_0092_NORM_Total             #####################
################################################################################################

# Load Data
Visvader_0092_NORM_Total <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Total/Visvader_0092_NORM_Total.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Visvader_0092_NORM_Total <- paramSweep(Visvader_0092_NORM_Total, PCs = 1:10, sct = FALSE)
sweep.stats.Visvader_0092_NORM_Total <- summarizeSweep(sweep.res.Visvader_0092_NORM_Total, GT = FALSE)
bcmvn_Visvader_0092_NORM_Total <- find.pK(sweep.stats.Visvader_0092_NORM_Total)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Visvader_0092_NORM_Total[which.max(bcmvn_Visvader_0092_NORM_Total$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Visvader_0092_NORM_Total@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 5% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Visvader dataset = 5k per sample, will use doublet rate for 6,000 cells to err on the side of caution
nExp_poi <- round(0.05*nrow(Visvader_0092_NORM_Total@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Visvader_0092_NORM_Total <- doubletFinder(Visvader_0092_NORM_Total, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Visvader_0092_NORM_Total_Quant <- (Visvader_0092_NORM_Total@meta.data$DF.classification == "Singlet")
Visvader_0092_NORM_Total_Quant_Singlets <- length(Visvader_0092_NORM_Total_Quant[Visvader_0092_NORM_Total_Quant== TRUE])
Visvader_0092_NORM_Total_Quant_Doublets <- length(Visvader_0092_NORM_Total_Quant[Visvader_0092_NORM_Total_Quant== FALSE])
Visvader_0092_NORM_Total_Quant_Doublets_Percent <- Visvader_0092_NORM_Total_Quant_Doublets /  (Visvader_0092_NORM_Total_Quant_Doublets + Visvader_0092_NORM_Total_Quant_Singlets) * 100
Visvader_0092_NORM_Total_Quant <- as.data.frame(c(Visvader_0092_NORM_Total_Quant_Singlets, Visvader_0092_NORM_Total_Quant_Doublets, Visvader_0092_NORM_Total_Quant_Doublets_Percent))
colnames(Visvader_0092_NORM_Total_Quant) <- c("Visvader_0092_NORM_Total_Quant")
rownames(Visvader_0092_NORM_Total_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Visvader_0092_NORM_Total_Quant, file = "/R/R_Visvader/Visvader_Doublet_Tables/Visvader_0092_NORM_Total_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Visvader_0092_NORM_Total_Singlets <- subset(Visvader_0092_NORM_Total, cells=rownames(Visvader_0092_NORM_Total@meta.data)[which(Visvader_0092_NORM_Total@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Visvader_0092_NORM_Total_Singlets, file = "/R/R_Visvader/Visvader_RDS/RDS_Total_Singlets/Visvader_0092_NORM_Total_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Visvader_0092_NORM_Total)
rm(Visvader_0092_NORM_Total_Singlets)
rm(Visvader_0092_NORM_Total_Quant)
rm(Visvader_0092_NORM_Total_Quant_Singlets)
rm(Visvader_0092_NORM_Total_Quant_Doublets)
rm(Visvader_0092_NORM_Total_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Visvader_0092_NORM_Total)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Visvader_0092_NORM_Total)
rm(sweep.stats.Visvader_0092_NORM_Total)
gc()




################################################################################################
########################              Visvader_0093_NORM_Total             #####################
################################################################################################

# Load Data
Visvader_0093_NORM_Total <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Total/Visvader_0093_NORM_Total.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Visvader_0093_NORM_Total <- paramSweep(Visvader_0093_NORM_Total, PCs = 1:10, sct = FALSE)
sweep.stats.Visvader_0093_NORM_Total <- summarizeSweep(sweep.res.Visvader_0093_NORM_Total, GT = FALSE)
bcmvn_Visvader_0093_NORM_Total <- find.pK(sweep.stats.Visvader_0093_NORM_Total)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Visvader_0093_NORM_Total[which.max(bcmvn_Visvader_0093_NORM_Total$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Visvader_0093_NORM_Total@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 5% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Visvader dataset = 5k per sample, will use doublet rate for 6,000 cells to err on the side of caution
nExp_poi <- round(0.05*nrow(Visvader_0093_NORM_Total@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Visvader_0093_NORM_Total <- doubletFinder(Visvader_0093_NORM_Total, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Visvader_0093_NORM_Total_Quant <- (Visvader_0093_NORM_Total@meta.data$DF.classification == "Singlet")
Visvader_0093_NORM_Total_Quant_Singlets <- length(Visvader_0093_NORM_Total_Quant[Visvader_0093_NORM_Total_Quant== TRUE])
Visvader_0093_NORM_Total_Quant_Doublets <- length(Visvader_0093_NORM_Total_Quant[Visvader_0093_NORM_Total_Quant== FALSE])
Visvader_0093_NORM_Total_Quant_Doublets_Percent <- Visvader_0093_NORM_Total_Quant_Doublets /  (Visvader_0093_NORM_Total_Quant_Doublets + Visvader_0093_NORM_Total_Quant_Singlets) * 100
Visvader_0093_NORM_Total_Quant <- as.data.frame(c(Visvader_0093_NORM_Total_Quant_Singlets, Visvader_0093_NORM_Total_Quant_Doublets, Visvader_0093_NORM_Total_Quant_Doublets_Percent))
colnames(Visvader_0093_NORM_Total_Quant) <- c("Visvader_0093_NORM_Total_Quant")
rownames(Visvader_0093_NORM_Total_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Visvader_0093_NORM_Total_Quant, file = "/R/R_Visvader/Visvader_Doublet_Tables/Visvader_0093_NORM_Total_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Visvader_0093_NORM_Total_Singlets <- subset(Visvader_0093_NORM_Total, cells=rownames(Visvader_0093_NORM_Total@meta.data)[which(Visvader_0093_NORM_Total@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Visvader_0093_NORM_Total_Singlets, file = "/R/R_Visvader/Visvader_RDS/RDS_Total_Singlets/Visvader_0093_NORM_Total_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Visvader_0093_NORM_Total)
rm(Visvader_0093_NORM_Total_Singlets)
rm(Visvader_0093_NORM_Total_Quant)
rm(Visvader_0093_NORM_Total_Quant_Singlets)
rm(Visvader_0093_NORM_Total_Quant_Doublets)
rm(Visvader_0093_NORM_Total_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Visvader_0093_NORM_Total)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Visvader_0093_NORM_Total)
rm(sweep.stats.Visvader_0093_NORM_Total)
gc()




################################################################################################
########################              Visvader_0123_NORM_Total             #####################
################################################################################################

# Load Data
Visvader_0123_NORM_Total <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Total/Visvader_0123_NORM_Total.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Visvader_0123_NORM_Total <- paramSweep(Visvader_0123_NORM_Total, PCs = 1:10, sct = FALSE)
sweep.stats.Visvader_0123_NORM_Total <- summarizeSweep(sweep.res.Visvader_0123_NORM_Total, GT = FALSE)
bcmvn_Visvader_0123_NORM_Total <- find.pK(sweep.stats.Visvader_0123_NORM_Total)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Visvader_0123_NORM_Total[which.max(bcmvn_Visvader_0123_NORM_Total$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Visvader_0123_NORM_Total@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 5% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Visvader dataset = 5k per sample, will use doublet rate for 6,000 cells to err on the side of caution
nExp_poi <- round(0.05*nrow(Visvader_0123_NORM_Total@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Visvader_0123_NORM_Total <- doubletFinder(Visvader_0123_NORM_Total, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Visvader_0123_NORM_Total_Quant <- (Visvader_0123_NORM_Total@meta.data$DF.classification == "Singlet")
Visvader_0123_NORM_Total_Quant_Singlets <- length(Visvader_0123_NORM_Total_Quant[Visvader_0123_NORM_Total_Quant== TRUE])
Visvader_0123_NORM_Total_Quant_Doublets <- length(Visvader_0123_NORM_Total_Quant[Visvader_0123_NORM_Total_Quant== FALSE])
Visvader_0123_NORM_Total_Quant_Doublets_Percent <- Visvader_0123_NORM_Total_Quant_Doublets /  (Visvader_0123_NORM_Total_Quant_Doublets + Visvader_0123_NORM_Total_Quant_Singlets) * 100
Visvader_0123_NORM_Total_Quant <- as.data.frame(c(Visvader_0123_NORM_Total_Quant_Singlets, Visvader_0123_NORM_Total_Quant_Doublets, Visvader_0123_NORM_Total_Quant_Doublets_Percent))
colnames(Visvader_0123_NORM_Total_Quant) <- c("Visvader_0123_NORM_Total_Quant")
rownames(Visvader_0123_NORM_Total_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Visvader_0123_NORM_Total_Quant, file = "/R/R_Visvader/Visvader_Doublet_Tables/Visvader_0123_NORM_Total_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Visvader_0123_NORM_Total_Singlets <- subset(Visvader_0123_NORM_Total, cells=rownames(Visvader_0123_NORM_Total@meta.data)[which(Visvader_0123_NORM_Total@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Visvader_0123_NORM_Total_Singlets, file = "/R/R_Visvader/Visvader_RDS/RDS_Total_Singlets/Visvader_0123_NORM_Total_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Visvader_0123_NORM_Total)
rm(Visvader_0123_NORM_Total_Singlets)
rm(Visvader_0123_NORM_Total_Quant)
rm(Visvader_0123_NORM_Total_Quant_Singlets)
rm(Visvader_0123_NORM_Total_Quant_Doublets)
rm(Visvader_0123_NORM_Total_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Visvader_0123_NORM_Total)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Visvader_0123_NORM_Total)
rm(sweep.stats.Visvader_0123_NORM_Total)
gc()




################################################################################################
########################              Visvader_0169_NORM_Total             #####################
################################################################################################

# Load Data
Visvader_0169_NORM_Total <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Total/Visvader_0169_NORM_Total.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Visvader_0169_NORM_Total <- paramSweep(Visvader_0169_NORM_Total, PCs = 1:10, sct = FALSE)
sweep.stats.Visvader_0169_NORM_Total <- summarizeSweep(sweep.res.Visvader_0169_NORM_Total, GT = FALSE)
bcmvn_Visvader_0169_NORM_Total <- find.pK(sweep.stats.Visvader_0169_NORM_Total)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Visvader_0169_NORM_Total[which.max(bcmvn_Visvader_0169_NORM_Total$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Visvader_0169_NORM_Total@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 5% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Visvader dataset = 5k per sample, will use doublet rate for 6,000 cells to err on the side of caution
nExp_poi <- round(0.05*nrow(Visvader_0169_NORM_Total@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Visvader_0169_NORM_Total <- doubletFinder(Visvader_0169_NORM_Total, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Visvader_0169_NORM_Total_Quant <- (Visvader_0169_NORM_Total@meta.data$DF.classification == "Singlet")
Visvader_0169_NORM_Total_Quant_Singlets <- length(Visvader_0169_NORM_Total_Quant[Visvader_0169_NORM_Total_Quant== TRUE])
Visvader_0169_NORM_Total_Quant_Doublets <- length(Visvader_0169_NORM_Total_Quant[Visvader_0169_NORM_Total_Quant== FALSE])
Visvader_0169_NORM_Total_Quant_Doublets_Percent <- Visvader_0169_NORM_Total_Quant_Doublets /  (Visvader_0169_NORM_Total_Quant_Doublets + Visvader_0169_NORM_Total_Quant_Singlets) * 100
Visvader_0169_NORM_Total_Quant <- as.data.frame(c(Visvader_0169_NORM_Total_Quant_Singlets, Visvader_0169_NORM_Total_Quant_Doublets, Visvader_0169_NORM_Total_Quant_Doublets_Percent))
colnames(Visvader_0169_NORM_Total_Quant) <- c("Visvader_0169_NORM_Total_Quant")
rownames(Visvader_0169_NORM_Total_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Visvader_0169_NORM_Total_Quant, file = "/R/R_Visvader/Visvader_Doublet_Tables/Visvader_0169_NORM_Total_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Visvader_0169_NORM_Total_Singlets <- subset(Visvader_0169_NORM_Total, cells=rownames(Visvader_0169_NORM_Total@meta.data)[which(Visvader_0169_NORM_Total@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Visvader_0169_NORM_Total_Singlets, file = "/R/R_Visvader/Visvader_RDS/RDS_Total_Singlets/Visvader_0169_NORM_Total_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Visvader_0169_NORM_Total)
rm(Visvader_0169_NORM_Total_Singlets)
rm(Visvader_0169_NORM_Total_Quant)
rm(Visvader_0169_NORM_Total_Quant_Singlets)
rm(Visvader_0169_NORM_Total_Quant_Doublets)
rm(Visvader_0169_NORM_Total_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Visvader_0169_NORM_Total)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Visvader_0169_NORM_Total)
rm(sweep.stats.Visvader_0169_NORM_Total)
gc()




################################################################################################
########################              Visvader_0230_NORM_Total             #####################
################################################################################################

# Load Data
Visvader_0230_NORM_Total <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Total/Visvader_0230_NORM_Total.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Visvader_0230_NORM_Total <- paramSweep(Visvader_0230_NORM_Total, PCs = 1:10, sct = FALSE)
sweep.stats.Visvader_0230_NORM_Total <- summarizeSweep(sweep.res.Visvader_0230_NORM_Total, GT = FALSE)
bcmvn_Visvader_0230_NORM_Total <- find.pK(sweep.stats.Visvader_0230_NORM_Total)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Visvader_0230_NORM_Total[which.max(bcmvn_Visvader_0230_NORM_Total$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Visvader_0230_NORM_Total@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 5% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Visvader dataset = 5k per sample, will use doublet rate for 6,000 cells to err on the side of caution
nExp_poi <- round(0.05*nrow(Visvader_0230_NORM_Total@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Visvader_0230_NORM_Total <- doubletFinder(Visvader_0230_NORM_Total, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Visvader_0230_NORM_Total_Quant <- (Visvader_0230_NORM_Total@meta.data$DF.classification == "Singlet")
Visvader_0230_NORM_Total_Quant_Singlets <- length(Visvader_0230_NORM_Total_Quant[Visvader_0230_NORM_Total_Quant== TRUE])
Visvader_0230_NORM_Total_Quant_Doublets <- length(Visvader_0230_NORM_Total_Quant[Visvader_0230_NORM_Total_Quant== FALSE])
Visvader_0230_NORM_Total_Quant_Doublets_Percent <- Visvader_0230_NORM_Total_Quant_Doublets /  (Visvader_0230_NORM_Total_Quant_Doublets + Visvader_0230_NORM_Total_Quant_Singlets) * 100
Visvader_0230_NORM_Total_Quant <- as.data.frame(c(Visvader_0230_NORM_Total_Quant_Singlets, Visvader_0230_NORM_Total_Quant_Doublets, Visvader_0230_NORM_Total_Quant_Doublets_Percent))
colnames(Visvader_0230_NORM_Total_Quant) <- c("Visvader_0230_NORM_Total_Quant")
rownames(Visvader_0230_NORM_Total_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Visvader_0230_NORM_Total_Quant, file = "/R/R_Visvader/Visvader_Doublet_Tables/Visvader_0230_NORM_Total_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Visvader_0230_NORM_Total_Singlets <- subset(Visvader_0230_NORM_Total, cells=rownames(Visvader_0230_NORM_Total@meta.data)[which(Visvader_0230_NORM_Total@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Visvader_0230_NORM_Total_Singlets, file = "/R/R_Visvader/Visvader_RDS/RDS_Total_Singlets/Visvader_0230_NORM_Total_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Visvader_0230_NORM_Total)
rm(Visvader_0230_NORM_Total_Singlets)
rm(Visvader_0230_NORM_Total_Quant)
rm(Visvader_0230_NORM_Total_Quant_Singlets)
rm(Visvader_0230_NORM_Total_Quant_Doublets)
rm(Visvader_0230_NORM_Total_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Visvader_0230_NORM_Total)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Visvader_0230_NORM_Total)
rm(sweep.stats.Visvader_0230_NORM_Total)
gc()




################################################################################################
########################              Visvader_0233_NORM_Total             #####################
################################################################################################

# Load Data
Visvader_0233_NORM_Total <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Total/Visvader_0233_NORM_Total.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Visvader_0233_NORM_Total <- paramSweep(Visvader_0233_NORM_Total, PCs = 1:10, sct = FALSE)
sweep.stats.Visvader_0233_NORM_Total <- summarizeSweep(sweep.res.Visvader_0233_NORM_Total, GT = FALSE)
bcmvn_Visvader_0233_NORM_Total <- find.pK(sweep.stats.Visvader_0233_NORM_Total)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Visvader_0233_NORM_Total[which.max(bcmvn_Visvader_0233_NORM_Total$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Visvader_0233_NORM_Total@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 5% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Visvader dataset = 5k per sample, will use doublet rate for 6,000 cells to err on the side of caution
nExp_poi <- round(0.05*nrow(Visvader_0233_NORM_Total@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Visvader_0233_NORM_Total <- doubletFinder(Visvader_0233_NORM_Total, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Visvader_0233_NORM_Total_Quant <- (Visvader_0233_NORM_Total@meta.data$DF.classification == "Singlet")
Visvader_0233_NORM_Total_Quant_Singlets <- length(Visvader_0233_NORM_Total_Quant[Visvader_0233_NORM_Total_Quant== TRUE])
Visvader_0233_NORM_Total_Quant_Doublets <- length(Visvader_0233_NORM_Total_Quant[Visvader_0233_NORM_Total_Quant== FALSE])
Visvader_0233_NORM_Total_Quant_Doublets_Percent <- Visvader_0233_NORM_Total_Quant_Doublets /  (Visvader_0233_NORM_Total_Quant_Doublets + Visvader_0233_NORM_Total_Quant_Singlets) * 100
Visvader_0233_NORM_Total_Quant <- as.data.frame(c(Visvader_0233_NORM_Total_Quant_Singlets, Visvader_0233_NORM_Total_Quant_Doublets, Visvader_0233_NORM_Total_Quant_Doublets_Percent))
colnames(Visvader_0233_NORM_Total_Quant) <- c("Visvader_0233_NORM_Total_Quant")
rownames(Visvader_0233_NORM_Total_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Visvader_0233_NORM_Total_Quant, file = "/R/R_Visvader/Visvader_Doublet_Tables/Visvader_0233_NORM_Total_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Visvader_0233_NORM_Total_Singlets <- subset(Visvader_0233_NORM_Total, cells=rownames(Visvader_0233_NORM_Total@meta.data)[which(Visvader_0233_NORM_Total@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Visvader_0233_NORM_Total_Singlets, file = "/R/R_Visvader/Visvader_RDS/RDS_Total_Singlets/Visvader_0233_NORM_Total_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Visvader_0233_NORM_Total)
rm(Visvader_0233_NORM_Total_Singlets)
rm(Visvader_0233_NORM_Total_Quant)
rm(Visvader_0233_NORM_Total_Quant_Singlets)
rm(Visvader_0233_NORM_Total_Quant_Doublets)
rm(Visvader_0233_NORM_Total_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Visvader_0233_NORM_Total)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Visvader_0233_NORM_Total)
rm(sweep.stats.Visvader_0233_NORM_Total)
gc()




################################################################################################
########################              Visvader_0275_NORM_Total             #####################
################################################################################################

# Load Data
Visvader_0275_NORM_Total <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Total/Visvader_0275_NORM_Total.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Visvader_0275_NORM_Total <- paramSweep(Visvader_0275_NORM_Total, PCs = 1:10, sct = FALSE)
sweep.stats.Visvader_0275_NORM_Total <- summarizeSweep(sweep.res.Visvader_0275_NORM_Total, GT = FALSE)
bcmvn_Visvader_0275_NORM_Total <- find.pK(sweep.stats.Visvader_0275_NORM_Total)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Visvader_0275_NORM_Total[which.max(bcmvn_Visvader_0275_NORM_Total$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Visvader_0275_NORM_Total@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 5% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Visvader dataset = 5k per sample, will use doublet rate for 6,000 cells to err on the side of caution
nExp_poi <- round(0.05*nrow(Visvader_0275_NORM_Total@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Visvader_0275_NORM_Total <- doubletFinder(Visvader_0275_NORM_Total, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Visvader_0275_NORM_Total_Quant <- (Visvader_0275_NORM_Total@meta.data$DF.classification == "Singlet")
Visvader_0275_NORM_Total_Quant_Singlets <- length(Visvader_0275_NORM_Total_Quant[Visvader_0275_NORM_Total_Quant== TRUE])
Visvader_0275_NORM_Total_Quant_Doublets <- length(Visvader_0275_NORM_Total_Quant[Visvader_0275_NORM_Total_Quant== FALSE])
Visvader_0275_NORM_Total_Quant_Doublets_Percent <- Visvader_0275_NORM_Total_Quant_Doublets /  (Visvader_0275_NORM_Total_Quant_Doublets + Visvader_0275_NORM_Total_Quant_Singlets) * 100
Visvader_0275_NORM_Total_Quant <- as.data.frame(c(Visvader_0275_NORM_Total_Quant_Singlets, Visvader_0275_NORM_Total_Quant_Doublets, Visvader_0275_NORM_Total_Quant_Doublets_Percent))
colnames(Visvader_0275_NORM_Total_Quant) <- c("Visvader_0275_NORM_Total_Quant")
rownames(Visvader_0275_NORM_Total_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Visvader_0275_NORM_Total_Quant, file = "/R/R_Visvader/Visvader_Doublet_Tables/Visvader_0275_NORM_Total_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Visvader_0275_NORM_Total_Singlets <- subset(Visvader_0275_NORM_Total, cells=rownames(Visvader_0275_NORM_Total@meta.data)[which(Visvader_0275_NORM_Total@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Visvader_0275_NORM_Total_Singlets, file = "/R/R_Visvader/Visvader_RDS/RDS_Total_Singlets/Visvader_0275_NORM_Total_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Visvader_0275_NORM_Total)
rm(Visvader_0275_NORM_Total_Singlets)
rm(Visvader_0275_NORM_Total_Quant)
rm(Visvader_0275_NORM_Total_Quant_Singlets)
rm(Visvader_0275_NORM_Total_Quant_Doublets)
rm(Visvader_0275_NORM_Total_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Visvader_0275_NORM_Total)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Visvader_0275_NORM_Total)
rm(sweep.stats.Visvader_0275_NORM_Total)
gc()




################################################################################################
########################              Visvader_0288_NORM_Total             #####################
################################################################################################

# Load Data
Visvader_0288_NORM_Total <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Total/Visvader_0288_NORM_Total.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Visvader_0288_NORM_Total <- paramSweep(Visvader_0288_NORM_Total, PCs = 1:10, sct = FALSE)
sweep.stats.Visvader_0288_NORM_Total <- summarizeSweep(sweep.res.Visvader_0288_NORM_Total, GT = FALSE)
bcmvn_Visvader_0288_NORM_Total <- find.pK(sweep.stats.Visvader_0288_NORM_Total)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Visvader_0288_NORM_Total[which.max(bcmvn_Visvader_0288_NORM_Total$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Visvader_0288_NORM_Total@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 5% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Visvader dataset = 5k per sample, will use doublet rate for 6,000 cells to err on the side of caution
nExp_poi <- round(0.05*nrow(Visvader_0288_NORM_Total@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Visvader_0288_NORM_Total <- doubletFinder(Visvader_0288_NORM_Total, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Visvader_0288_NORM_Total_Quant <- (Visvader_0288_NORM_Total@meta.data$DF.classification == "Singlet")
Visvader_0288_NORM_Total_Quant_Singlets <- length(Visvader_0288_NORM_Total_Quant[Visvader_0288_NORM_Total_Quant== TRUE])
Visvader_0288_NORM_Total_Quant_Doublets <- length(Visvader_0288_NORM_Total_Quant[Visvader_0288_NORM_Total_Quant== FALSE])
Visvader_0288_NORM_Total_Quant_Doublets_Percent <- Visvader_0288_NORM_Total_Quant_Doublets /  (Visvader_0288_NORM_Total_Quant_Doublets + Visvader_0288_NORM_Total_Quant_Singlets) * 100
Visvader_0288_NORM_Total_Quant <- as.data.frame(c(Visvader_0288_NORM_Total_Quant_Singlets, Visvader_0288_NORM_Total_Quant_Doublets, Visvader_0288_NORM_Total_Quant_Doublets_Percent))
colnames(Visvader_0288_NORM_Total_Quant) <- c("Visvader_0288_NORM_Total_Quant")
rownames(Visvader_0288_NORM_Total_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Visvader_0288_NORM_Total_Quant, file = "/R/R_Visvader/Visvader_Doublet_Tables/Visvader_0288_NORM_Total_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Visvader_0288_NORM_Total_Singlets <- subset(Visvader_0288_NORM_Total, cells=rownames(Visvader_0288_NORM_Total@meta.data)[which(Visvader_0288_NORM_Total@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Visvader_0288_NORM_Total_Singlets, file = "/R/R_Visvader/Visvader_RDS/RDS_Total_Singlets/Visvader_0288_NORM_Total_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Visvader_0288_NORM_Total)
rm(Visvader_0288_NORM_Total_Singlets)
rm(Visvader_0288_NORM_Total_Quant)
rm(Visvader_0288_NORM_Total_Quant_Singlets)
rm(Visvader_0288_NORM_Total_Quant_Doublets)
rm(Visvader_0288_NORM_Total_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Visvader_0288_NORM_Total)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Visvader_0288_NORM_Total)
rm(sweep.stats.Visvader_0288_NORM_Total)
gc()




################################################################################################
########################              Visvader_0342_NORM_Total             #####################
################################################################################################

# Load Data
Visvader_0342_NORM_Total <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Total/Visvader_0342_NORM_Total.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Visvader_0342_NORM_Total <- paramSweep(Visvader_0342_NORM_Total, PCs = 1:10, sct = FALSE)
sweep.stats.Visvader_0342_NORM_Total <- summarizeSweep(sweep.res.Visvader_0342_NORM_Total, GT = FALSE)
bcmvn_Visvader_0342_NORM_Total <- find.pK(sweep.stats.Visvader_0342_NORM_Total)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Visvader_0342_NORM_Total[which.max(bcmvn_Visvader_0342_NORM_Total$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Visvader_0342_NORM_Total@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 5% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Visvader dataset = 5k per sample, will use doublet rate for 6,000 cells to err on the side of caution
nExp_poi <- round(0.05*nrow(Visvader_0342_NORM_Total@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Visvader_0342_NORM_Total <- doubletFinder(Visvader_0342_NORM_Total, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Visvader_0342_NORM_Total_Quant <- (Visvader_0342_NORM_Total@meta.data$DF.classification == "Singlet")
Visvader_0342_NORM_Total_Quant_Singlets <- length(Visvader_0342_NORM_Total_Quant[Visvader_0342_NORM_Total_Quant== TRUE])
Visvader_0342_NORM_Total_Quant_Doublets <- length(Visvader_0342_NORM_Total_Quant[Visvader_0342_NORM_Total_Quant== FALSE])
Visvader_0342_NORM_Total_Quant_Doublets_Percent <- Visvader_0342_NORM_Total_Quant_Doublets /  (Visvader_0342_NORM_Total_Quant_Doublets + Visvader_0342_NORM_Total_Quant_Singlets) * 100
Visvader_0342_NORM_Total_Quant <- as.data.frame(c(Visvader_0342_NORM_Total_Quant_Singlets, Visvader_0342_NORM_Total_Quant_Doublets, Visvader_0342_NORM_Total_Quant_Doublets_Percent))
colnames(Visvader_0342_NORM_Total_Quant) <- c("Visvader_0342_NORM_Total_Quant")
rownames(Visvader_0342_NORM_Total_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Visvader_0342_NORM_Total_Quant, file = "/R/R_Visvader/Visvader_Doublet_Tables/Visvader_0342_NORM_Total_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Visvader_0342_NORM_Total_Singlets <- subset(Visvader_0342_NORM_Total, cells=rownames(Visvader_0342_NORM_Total@meta.data)[which(Visvader_0342_NORM_Total@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Visvader_0342_NORM_Total_Singlets, file = "/R/R_Visvader/Visvader_RDS/RDS_Total_Singlets/Visvader_0342_NORM_Total_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Visvader_0342_NORM_Total)
rm(Visvader_0342_NORM_Total_Singlets)
rm(Visvader_0342_NORM_Total_Quant)
rm(Visvader_0342_NORM_Total_Quant_Singlets)
rm(Visvader_0342_NORM_Total_Quant_Doublets)
rm(Visvader_0342_NORM_Total_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Visvader_0342_NORM_Total)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Visvader_0342_NORM_Total)
rm(sweep.stats.Visvader_0342_NORM_Total)
gc()



################################################################################################
########################              Visvader_0372_NORM_Total             #####################
################################################################################################

# Load Data
Visvader_0372_NORM_Total <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Total/Visvader_0372_NORM_Total.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Visvader_0372_NORM_Total <- paramSweep(Visvader_0372_NORM_Total, PCs = 1:10, sct = FALSE)
sweep.stats.Visvader_0372_NORM_Total <- summarizeSweep(sweep.res.Visvader_0372_NORM_Total, GT = FALSE)
bcmvn_Visvader_0372_NORM_Total <- find.pK(sweep.stats.Visvader_0372_NORM_Total)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Visvader_0372_NORM_Total[which.max(bcmvn_Visvader_0372_NORM_Total$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Visvader_0372_NORM_Total@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 5% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Visvader dataset = 5k per sample, will use doublet rate for 6,000 cells to err on the side of caution
nExp_poi <- round(0.05*nrow(Visvader_0372_NORM_Total@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Visvader_0372_NORM_Total <- doubletFinder(Visvader_0372_NORM_Total, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Visvader_0372_NORM_Total_Quant <- (Visvader_0372_NORM_Total@meta.data$DF.classification == "Singlet")
Visvader_0372_NORM_Total_Quant_Singlets <- length(Visvader_0372_NORM_Total_Quant[Visvader_0372_NORM_Total_Quant== TRUE])
Visvader_0372_NORM_Total_Quant_Doublets <- length(Visvader_0372_NORM_Total_Quant[Visvader_0372_NORM_Total_Quant== FALSE])
Visvader_0372_NORM_Total_Quant_Doublets_Percent <- Visvader_0372_NORM_Total_Quant_Doublets /  (Visvader_0372_NORM_Total_Quant_Doublets + Visvader_0372_NORM_Total_Quant_Singlets) * 100
Visvader_0372_NORM_Total_Quant <- as.data.frame(c(Visvader_0372_NORM_Total_Quant_Singlets, Visvader_0372_NORM_Total_Quant_Doublets, Visvader_0372_NORM_Total_Quant_Doublets_Percent))
colnames(Visvader_0372_NORM_Total_Quant) <- c("Visvader_0372_NORM_Total_Quant")
rownames(Visvader_0372_NORM_Total_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Visvader_0372_NORM_Total_Quant, file = "/R/R_Visvader/Visvader_Doublet_Tables/Visvader_0372_NORM_Total_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Visvader_0372_NORM_Total_Singlets <- subset(Visvader_0372_NORM_Total, cells=rownames(Visvader_0372_NORM_Total@meta.data)[which(Visvader_0372_NORM_Total@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Visvader_0372_NORM_Total_Singlets, file = "/R/R_Visvader/Visvader_RDS/RDS_Total_Singlets/Visvader_0372_NORM_Total_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Visvader_0372_NORM_Total)
rm(Visvader_0372_NORM_Total_Singlets)
rm(Visvader_0372_NORM_Total_Quant)
rm(Visvader_0372_NORM_Total_Quant_Singlets)
rm(Visvader_0372_NORM_Total_Quant_Doublets)
rm(Visvader_0372_NORM_Total_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Visvader_0372_NORM_Total)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Visvader_0372_NORM_Total)
rm(sweep.stats.Visvader_0372_NORM_Total)
gc()


