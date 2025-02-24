#####################################################################################################################
#                       Visvader (Pal et al) Dataset Breast Tumor Analysis Steps                                    #
#####################################################################################################################
# Step 1: Process Each Breast Tumor Count Matrix into Seurat Object                                                 #
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


#########################################################
############# This is ER/PR Tumor Batch # ###############
#########################################################

##################################################################
########################### Visvader_0001_ER_Total ###############
##################################################################

# Load the dataset
Visvader_0001_ER_Total <- Read10X(data.dir = "/R/R_Visvader/Visvader_Input/ER/Visvader_0001_ER_Total/")

# Initialize the Seurat object with the raw (non-normalized data).
Visvader_0001_ER_Total <- CreateSeuratObject(counts = Visvader_0001_ER_Total, project = "Visvader_0001_ER_Total", min.cells = 3, min.features = 200)
Visvader_0001_ER_Total

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Visvader_0001_ER_Total[["percent.mt"]] <- PercentageFeatureSet(Visvader_0001_ER_Total, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Visvader_0001_ER_Total, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Visvader_0001_ER_Total, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Visvader_0001_ER_Total, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Visvader_0001_ER_Total <- subset(Visvader_0001_ER_Total, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Visvader_0001_ER_Total <- NormalizeData(Visvader_0001_ER_Total, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Visvader_0001_ER_Total <- FindVariableFeatures(Visvader_0001_ER_Total, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Visvader_0001_ER_Total)
Visvader_0001_ER_Total <- ScaleData(Visvader_0001_ER_Total, features = all.genes)

# Perform linear dimensional reduction
Visvader_0001_ER_Total <- RunPCA(Visvader_0001_ER_Total, features = VariableFeatures(object = Visvader_0001_ER_Total))

# Examine and visualize PCA results a few different ways
print(Visvader_0001_ER_Total[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Visvader_0001_ER_Total, dims = 1:2, reduction = "pca")
DimPlot(Visvader_0001_ER_Total, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Visvader_0001_ER_Total <- JackStraw(Visvader_0001_ER_Total, num.replicate = 100)
Visvader_0001_ER_Total <- ScoreJackStraw(Visvader_0001_ER_Total, dims = 1:20)
ElbowPlot(Visvader_0001_ER_Total)

# Cluster cells
Visvader_0001_ER_Total <- FindNeighbors(Visvader_0001_ER_Total, dims = 1:20)
Visvader_0001_ER_Total <- FindClusters(Visvader_0001_ER_Total, resolution = 0.5)

# Run non-linear dimensional reduction
Visvader_0001_ER_Total <- RunUMAP(Visvader_0001_ER_Total, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Visvader_0001_ER_Total, reduction = "umap")

# Save RDS
saveRDS(Visvader_0001_ER_Total, file = "/R/R_Visvader/Visvader_RDS/RDS_Total/Visvader_0001_ER_Total.rds")



##################################################################
########################### Visvader_0025_ER_Total ###############
##################################################################

# Load the dataset
Visvader_0025_ER_Total <- Read10X(data.dir = "/R/R_Visvader/Visvader_Input/ER/Visvader_0025_ER_Total/")

# Initialize the Seurat object with the raw (non-normalized data).
Visvader_0025_ER_Total <- CreateSeuratObject(counts = Visvader_0025_ER_Total, project = "Visvader_0025_ER_Total", min.cells = 3, min.features = 200)
Visvader_0025_ER_Total

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Visvader_0025_ER_Total[["percent.mt"]] <- PercentageFeatureSet(Visvader_0025_ER_Total, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Visvader_0025_ER_Total, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Visvader_0025_ER_Total, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Visvader_0025_ER_Total, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Visvader_0025_ER_Total <- subset(Visvader_0025_ER_Total, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Visvader_0025_ER_Total <- NormalizeData(Visvader_0025_ER_Total, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Visvader_0025_ER_Total <- FindVariableFeatures(Visvader_0025_ER_Total, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Visvader_0025_ER_Total)
Visvader_0025_ER_Total <- ScaleData(Visvader_0025_ER_Total, features = all.genes)

# Perform linear dimensional reduction
Visvader_0025_ER_Total <- RunPCA(Visvader_0025_ER_Total, features = VariableFeatures(object = Visvader_0025_ER_Total))

# Examine and visualize PCA results a few different ways
print(Visvader_0025_ER_Total[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Visvader_0025_ER_Total, dims = 1:2, reduction = "pca")
DimPlot(Visvader_0025_ER_Total, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Visvader_0025_ER_Total <- JackStraw(Visvader_0025_ER_Total, num.replicate = 100)
Visvader_0025_ER_Total <- ScoreJackStraw(Visvader_0025_ER_Total, dims = 1:20)
ElbowPlot(Visvader_0025_ER_Total)

# Cluster cells
Visvader_0025_ER_Total <- FindNeighbors(Visvader_0025_ER_Total, dims = 1:20)
Visvader_0025_ER_Total <- FindClusters(Visvader_0025_ER_Total, resolution = 0.5)

# Run non-linear dimensional reduction
Visvader_0025_ER_Total <- RunUMAP(Visvader_0025_ER_Total, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Visvader_0025_ER_Total, reduction = "umap")

# Save RDS
saveRDS(Visvader_0025_ER_Total, file = "/R/R_Visvader/Visvader_RDS/RDS_Total/Visvader_0025_ER_Total.rds")



##################################################################
########################### Visvader_0029_7C_ER_Total ############
##################################################################

# Load the dataset
Visvader_0029_7C_ER_Total <- Read10X(data.dir = "/R/R_Visvader/Visvader_Input/ER/Visvader_0029-7C_ER_Total/")

# Initialize the Seurat object with the raw (non-normalized data).
Visvader_0029_7C_ER_Total <- CreateSeuratObject(counts = Visvader_0029_7C_ER_Total, project = "Visvader_0029_7C_ER_Total", min.cells = 3, min.features = 200)
Visvader_0029_7C_ER_Total

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Visvader_0029_7C_ER_Total[["percent.mt"]] <- PercentageFeatureSet(Visvader_0029_7C_ER_Total, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Visvader_0029_7C_ER_Total, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Visvader_0029_7C_ER_Total, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Visvader_0029_7C_ER_Total, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Visvader_0029_7C_ER_Total <- subset(Visvader_0029_7C_ER_Total, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Visvader_0029_7C_ER_Total <- NormalizeData(Visvader_0029_7C_ER_Total, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Visvader_0029_7C_ER_Total <- FindVariableFeatures(Visvader_0029_7C_ER_Total, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Visvader_0029_7C_ER_Total)
Visvader_0029_7C_ER_Total <- ScaleData(Visvader_0029_7C_ER_Total, features = all.genes)

# Perform linear dimensional reduction
Visvader_0029_7C_ER_Total <- RunPCA(Visvader_0029_7C_ER_Total, features = VariableFeatures(object = Visvader_0029_7C_ER_Total))

# Examine and visualize PCA results a few different ways
print(Visvader_0029_7C_ER_Total[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Visvader_0029_7C_ER_Total, dims = 1:2, reduction = "pca")
DimPlot(Visvader_0029_7C_ER_Total, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Visvader_0029_7C_ER_Total <- JackStraw(Visvader_0029_7C_ER_Total, num.replicate = 100)
Visvader_0029_7C_ER_Total <- ScoreJackStraw(Visvader_0029_7C_ER_Total, dims = 1:20)
ElbowPlot(Visvader_0029_7C_ER_Total)

# Cluster cells
Visvader_0029_7C_ER_Total <- FindNeighbors(Visvader_0029_7C_ER_Total, dims = 1:20)
Visvader_0029_7C_ER_Total <- FindClusters(Visvader_0029_7C_ER_Total, resolution = 0.5)

# Run non-linear dimensional reduction
Visvader_0029_7C_ER_Total <- RunUMAP(Visvader_0029_7C_ER_Total, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Visvader_0029_7C_ER_Total, reduction = "umap")

# Save RDS
saveRDS(Visvader_0029_7C_ER_Total, file = "/R/R_Visvader/Visvader_RDS/RDS_Total/Visvader_0029_7C_ER_Total.rds")




##################################################################
########################### Visvader_0029_9C_ER_Total ############
##################################################################

# Load the dataset
Visvader_0029_9C_ER_Total <- Read10X(data.dir = "/R/R_Visvader/Visvader_Input/ER/Visvader_0029-9C_ER_Total/")

# Initialize the Seurat object with the raw (non-normalized data).
Visvader_0029_9C_ER_Total <- CreateSeuratObject(counts = Visvader_0029_9C_ER_Total, project = "Visvader_0029_9C_ER_Total", min.cells = 3, min.features = 200)
Visvader_0029_9C_ER_Total

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Visvader_0029_9C_ER_Total[["percent.mt"]] <- PercentageFeatureSet(Visvader_0029_9C_ER_Total, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Visvader_0029_9C_ER_Total, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Visvader_0029_9C_ER_Total, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Visvader_0029_9C_ER_Total, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Visvader_0029_9C_ER_Total <- subset(Visvader_0029_9C_ER_Total, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Visvader_0029_9C_ER_Total <- NormalizeData(Visvader_0029_9C_ER_Total, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Visvader_0029_9C_ER_Total <- FindVariableFeatures(Visvader_0029_9C_ER_Total, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Visvader_0029_9C_ER_Total)
Visvader_0029_9C_ER_Total <- ScaleData(Visvader_0029_9C_ER_Total, features = all.genes)

# Perform linear dimensional reduction
Visvader_0029_9C_ER_Total <- RunPCA(Visvader_0029_9C_ER_Total, features = VariableFeatures(object = Visvader_0029_9C_ER_Total))

# Examine and visualize PCA results a few different ways
print(Visvader_0029_9C_ER_Total[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Visvader_0029_9C_ER_Total, dims = 1:2, reduction = "pca")
DimPlot(Visvader_0029_9C_ER_Total, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Visvader_0029_9C_ER_Total <- JackStraw(Visvader_0029_9C_ER_Total, num.replicate = 100)
Visvader_0029_9C_ER_Total <- ScoreJackStraw(Visvader_0029_9C_ER_Total, dims = 1:20)
ElbowPlot(Visvader_0029_9C_ER_Total)

# Cluster cells
Visvader_0029_9C_ER_Total <- FindNeighbors(Visvader_0029_9C_ER_Total, dims = 1:20)
Visvader_0029_9C_ER_Total <- FindClusters(Visvader_0029_9C_ER_Total, resolution = 0.5)

# Run non-linear dimensional reduction
Visvader_0029_9C_ER_Total <- RunUMAP(Visvader_0029_9C_ER_Total, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Visvader_0029_9C_ER_Total, reduction = "umap")

# Save RDS
saveRDS(Visvader_0029_9C_ER_Total, file = "/R/R_Visvader/Visvader_RDS/RDS_Total/Visvader_0029_9C_ER_Total.rds")




##################################################################
########################### Visvader_0032_ER_Total ###############
##################################################################

# Load the dataset
Visvader_0032_ER_Total <- Read10X(data.dir = "/R/R_Visvader/Visvader_Input/ER/Visvader_0032_ER_Total/")

# Initialize the Seurat object with the raw (non-normalized data).
Visvader_0032_ER_Total <- CreateSeuratObject(counts = Visvader_0032_ER_Total, project = "Visvader_0032_ER_Total", min.cells = 3, min.features = 200)
Visvader_0032_ER_Total

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Visvader_0032_ER_Total[["percent.mt"]] <- PercentageFeatureSet(Visvader_0032_ER_Total, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Visvader_0032_ER_Total, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Visvader_0032_ER_Total, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Visvader_0032_ER_Total, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Visvader_0032_ER_Total <- subset(Visvader_0032_ER_Total, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Visvader_0032_ER_Total <- NormalizeData(Visvader_0032_ER_Total, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Visvader_0032_ER_Total <- FindVariableFeatures(Visvader_0032_ER_Total, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Visvader_0032_ER_Total)
Visvader_0032_ER_Total <- ScaleData(Visvader_0032_ER_Total, features = all.genes)

# Perform linear dimensional reduction
Visvader_0032_ER_Total <- RunPCA(Visvader_0032_ER_Total, features = VariableFeatures(object = Visvader_0032_ER_Total))

# Examine and visualize PCA results a few different ways
print(Visvader_0032_ER_Total[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Visvader_0032_ER_Total, dims = 1:2, reduction = "pca")
DimPlot(Visvader_0032_ER_Total, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Visvader_0032_ER_Total <- JackStraw(Visvader_0032_ER_Total, num.replicate = 100)
Visvader_0032_ER_Total <- ScoreJackStraw(Visvader_0032_ER_Total, dims = 1:20)
ElbowPlot(Visvader_0032_ER_Total)

# Cluster cells
Visvader_0032_ER_Total <- FindNeighbors(Visvader_0032_ER_Total, dims = 1:20)
Visvader_0032_ER_Total <- FindClusters(Visvader_0032_ER_Total, resolution = 0.5)

# Run non-linear dimensional reduction
Visvader_0032_ER_Total <- RunUMAP(Visvader_0032_ER_Total, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Visvader_0032_ER_Total, reduction = "umap")

# Save RDS
saveRDS(Visvader_0032_ER_Total, file = "/R/R_Visvader/Visvader_RDS/RDS_Total/Visvader_0032_ER_Total.rds")



##################################################################
########################### Visvader_0040_ER_Total ###############
##################################################################

# Load the dataset
Visvader_0040_ER_Total <- Read10X(data.dir = "/R/R_Visvader/Visvader_Input/ER/Visvader_0040_ER_Total/")

# Initialize the Seurat object with the raw (non-normalized data).
Visvader_0040_ER_Total <- CreateSeuratObject(counts = Visvader_0040_ER_Total, project = "Visvader_0040_ER_Total", min.cells = 3, min.features = 200)
Visvader_0040_ER_Total

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Visvader_0040_ER_Total[["percent.mt"]] <- PercentageFeatureSet(Visvader_0040_ER_Total, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Visvader_0040_ER_Total, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Visvader_0040_ER_Total, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Visvader_0040_ER_Total, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Visvader_0040_ER_Total <- subset(Visvader_0040_ER_Total, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Visvader_0040_ER_Total <- NormalizeData(Visvader_0040_ER_Total, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Visvader_0040_ER_Total <- FindVariableFeatures(Visvader_0040_ER_Total, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Visvader_0040_ER_Total)
Visvader_0040_ER_Total <- ScaleData(Visvader_0040_ER_Total, features = all.genes)

# Perform linear dimensional reduction
Visvader_0040_ER_Total <- RunPCA(Visvader_0040_ER_Total, features = VariableFeatures(object = Visvader_0040_ER_Total))

# Examine and visualize PCA results a few different ways
print(Visvader_0040_ER_Total[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Visvader_0040_ER_Total, dims = 1:2, reduction = "pca")
DimPlot(Visvader_0040_ER_Total, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Visvader_0040_ER_Total <- JackStraw(Visvader_0040_ER_Total, num.replicate = 100)
Visvader_0040_ER_Total <- ScoreJackStraw(Visvader_0040_ER_Total, dims = 1:20)
ElbowPlot(Visvader_0040_ER_Total)

# Cluster cells
Visvader_0040_ER_Total <- FindNeighbors(Visvader_0040_ER_Total, dims = 1:20)
Visvader_0040_ER_Total <- FindClusters(Visvader_0040_ER_Total, resolution = 0.5)

# Run non-linear dimensional reduction
Visvader_0040_ER_Total <- RunUMAP(Visvader_0040_ER_Total, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Visvader_0040_ER_Total, reduction = "umap")

# Save RDS
saveRDS(Visvader_0040_ER_Total, file = "/R/R_Visvader/Visvader_RDS/RDS_Total/Visvader_0040_ER_Total.rds")






##################################################################
########################### Visvader_0042_ER_Total ###############
##################################################################

# Load the dataset
Visvader_0042_ER_Total <- Read10X(data.dir = "/R/R_Visvader/Visvader_Input/ER/Visvader_0042_ER_Total/")

# Initialize the Seurat object with the raw (non-normalized data).
Visvader_0042_ER_Total <- CreateSeuratObject(counts = Visvader_0042_ER_Total, project = "Visvader_0042_ER_Total", min.cells = 3, min.features = 200)
Visvader_0042_ER_Total

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Visvader_0042_ER_Total[["percent.mt"]] <- PercentageFeatureSet(Visvader_0042_ER_Total, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Visvader_0042_ER_Total, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Visvader_0042_ER_Total, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Visvader_0042_ER_Total, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Visvader_0042_ER_Total <- subset(Visvader_0042_ER_Total, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Visvader_0042_ER_Total <- NormalizeData(Visvader_0042_ER_Total, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Visvader_0042_ER_Total <- FindVariableFeatures(Visvader_0042_ER_Total, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Visvader_0042_ER_Total)
Visvader_0042_ER_Total <- ScaleData(Visvader_0042_ER_Total, features = all.genes)

# Perform linear dimensional reduction
Visvader_0042_ER_Total <- RunPCA(Visvader_0042_ER_Total, features = VariableFeatures(object = Visvader_0042_ER_Total))

# Examine and visualize PCA results a few different ways
print(Visvader_0042_ER_Total[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Visvader_0042_ER_Total, dims = 1:2, reduction = "pca")
DimPlot(Visvader_0042_ER_Total, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Visvader_0042_ER_Total <- JackStraw(Visvader_0042_ER_Total, num.replicate = 100)
Visvader_0042_ER_Total <- ScoreJackStraw(Visvader_0042_ER_Total, dims = 1:20)
ElbowPlot(Visvader_0042_ER_Total)

# Cluster cells
Visvader_0042_ER_Total <- FindNeighbors(Visvader_0042_ER_Total, dims = 1:20)
Visvader_0042_ER_Total <- FindClusters(Visvader_0042_ER_Total, resolution = 0.5)

# Run non-linear dimensional reduction
Visvader_0042_ER_Total <- RunUMAP(Visvader_0042_ER_Total, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Visvader_0042_ER_Total, reduction = "umap")

# Save RDS
saveRDS(Visvader_0042_ER_Total, file = "/R/R_Visvader/Visvader_RDS/RDS_Total/Visvader_0042_ER_Total.rds")


##################################################################
########################### Visvader_0043_ER_Total ###############
##################################################################

# Load the dataset
Visvader_0043_ER_Total <- Read10X(data.dir = "/R/R_Visvader/Visvader_Input/ER/Visvader_0043_ER_Total/")

# Initialize the Seurat object with the raw (non-normalized data).
Visvader_0043_ER_Total <- CreateSeuratObject(counts = Visvader_0043_ER_Total, project = "Visvader_0043_ER_Total", min.cells = 3, min.features = 200)
Visvader_0043_ER_Total

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Visvader_0043_ER_Total[["percent.mt"]] <- PercentageFeatureSet(Visvader_0043_ER_Total, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Visvader_0043_ER_Total, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Visvader_0043_ER_Total, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Visvader_0043_ER_Total, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Visvader_0043_ER_Total <- subset(Visvader_0043_ER_Total, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Visvader_0043_ER_Total <- NormalizeData(Visvader_0043_ER_Total, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Visvader_0043_ER_Total <- FindVariableFeatures(Visvader_0043_ER_Total, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Visvader_0043_ER_Total)
Visvader_0043_ER_Total <- ScaleData(Visvader_0043_ER_Total, features = all.genes)

# Perform linear dimensional reduction
Visvader_0043_ER_Total <- RunPCA(Visvader_0043_ER_Total, features = VariableFeatures(object = Visvader_0043_ER_Total))

# Examine and visualize PCA results a few different ways
print(Visvader_0043_ER_Total[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Visvader_0043_ER_Total, dims = 1:2, reduction = "pca")
DimPlot(Visvader_0043_ER_Total, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Visvader_0043_ER_Total <- JackStraw(Visvader_0043_ER_Total, num.replicate = 100)
Visvader_0043_ER_Total <- ScoreJackStraw(Visvader_0043_ER_Total, dims = 1:20)
ElbowPlot(Visvader_0043_ER_Total)

# Cluster cells
Visvader_0043_ER_Total <- FindNeighbors(Visvader_0043_ER_Total, dims = 1:20)
Visvader_0043_ER_Total <- FindClusters(Visvader_0043_ER_Total, resolution = 0.5)

# Run non-linear dimensional reduction
Visvader_0043_ER_Total <- RunUMAP(Visvader_0043_ER_Total, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Visvader_0043_ER_Total, reduction = "umap")

# Save RDS
saveRDS(Visvader_0043_ER_Total, file = "/R/R_Visvader/Visvader_RDS/RDS_Total/Visvader_0043_ER_Total.rds")




##################################################################
########################### Visvader_0056_ER_Total ###############
##################################################################

# Load the dataset
Visvader_0056_ER_Total <- Read10X(data.dir = "/R/R_Visvader/Visvader_Input/ER/Visvader_0056_ER_Total/")

# Initialize the Seurat object with the raw (non-normalized data).
Visvader_0056_ER_Total <- CreateSeuratObject(counts = Visvader_0056_ER_Total, project = "Visvader_0056_ER_Total", min.cells = 3, min.features = 200)
Visvader_0056_ER_Total

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Visvader_0056_ER_Total[["percent.mt"]] <- PercentageFeatureSet(Visvader_0056_ER_Total, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Visvader_0056_ER_Total, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Visvader_0056_ER_Total, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Visvader_0056_ER_Total, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Visvader_0056_ER_Total <- subset(Visvader_0056_ER_Total, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Visvader_0056_ER_Total <- NormalizeData(Visvader_0056_ER_Total, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Visvader_0056_ER_Total <- FindVariableFeatures(Visvader_0056_ER_Total, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Visvader_0056_ER_Total)
Visvader_0056_ER_Total <- ScaleData(Visvader_0056_ER_Total, features = all.genes)

# Perform linear dimensional reduction
Visvader_0056_ER_Total <- RunPCA(Visvader_0056_ER_Total, features = VariableFeatures(object = Visvader_0056_ER_Total))

# Examine and visualize PCA results a few different ways
print(Visvader_0056_ER_Total[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Visvader_0056_ER_Total, dims = 1:2, reduction = "pca")
DimPlot(Visvader_0056_ER_Total, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Visvader_0056_ER_Total <- JackStraw(Visvader_0056_ER_Total, num.replicate = 100)
Visvader_0056_ER_Total <- ScoreJackStraw(Visvader_0056_ER_Total, dims = 1:20)
ElbowPlot(Visvader_0056_ER_Total)

# Cluster cells
Visvader_0056_ER_Total <- FindNeighbors(Visvader_0056_ER_Total, dims = 1:20)
Visvader_0056_ER_Total <- FindClusters(Visvader_0056_ER_Total, resolution = 0.5)

# Run non-linear dimensional reduction
Visvader_0056_ER_Total <- RunUMAP(Visvader_0056_ER_Total, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Visvader_0056_ER_Total, reduction = "umap")

# Save RDS
saveRDS(Visvader_0056_ER_Total, file = "/R/R_Visvader/Visvader_RDS/RDS_Total/Visvader_0056_ER_Total.rds")



##################################################################
########################### Visvader_0064_ER_Total ###############
##################################################################

# Load the dataset
Visvader_0064_ER_Total <- Read10X(data.dir = "/R/R_Visvader/Visvader_Input/ER/Visvader_0064_ER_Total/")

# Initialize the Seurat object with the raw (non-normalized data).
Visvader_0064_ER_Total <- CreateSeuratObject(counts = Visvader_0064_ER_Total, project = "Visvader_0064_ER_Total", min.cells = 3, min.features = 200)
Visvader_0064_ER_Total

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Visvader_0064_ER_Total[["percent.mt"]] <- PercentageFeatureSet(Visvader_0064_ER_Total, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Visvader_0064_ER_Total, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Visvader_0064_ER_Total, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Visvader_0064_ER_Total, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Visvader_0064_ER_Total <- subset(Visvader_0064_ER_Total, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Visvader_0064_ER_Total <- NormalizeData(Visvader_0064_ER_Total, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Visvader_0064_ER_Total <- FindVariableFeatures(Visvader_0064_ER_Total, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Visvader_0064_ER_Total)
Visvader_0064_ER_Total <- ScaleData(Visvader_0064_ER_Total, features = all.genes)

# Perform linear dimensional reduction
Visvader_0064_ER_Total <- RunPCA(Visvader_0064_ER_Total, features = VariableFeatures(object = Visvader_0064_ER_Total))

# Examine and visualize PCA results a few different ways
print(Visvader_0064_ER_Total[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Visvader_0064_ER_Total, dims = 1:2, reduction = "pca")
DimPlot(Visvader_0064_ER_Total, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Visvader_0064_ER_Total <- JackStraw(Visvader_0064_ER_Total, num.replicate = 100)
Visvader_0064_ER_Total <- ScoreJackStraw(Visvader_0064_ER_Total, dims = 1:20)
ElbowPlot(Visvader_0064_ER_Total)

# Cluster cells
Visvader_0064_ER_Total <- FindNeighbors(Visvader_0064_ER_Total, dims = 1:20)
Visvader_0064_ER_Total <- FindClusters(Visvader_0064_ER_Total, resolution = 0.5)

# Run non-linear dimensional reduction
Visvader_0064_ER_Total <- RunUMAP(Visvader_0064_ER_Total, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Visvader_0064_ER_Total, reduction = "umap")

# Save RDS
saveRDS(Visvader_0064_ER_Total, file = "/R/R_Visvader/Visvader_RDS/RDS_Total/Visvader_0064_ER_Total.rds")



##################################################################
########################### Visvader_0068_ER_Total ###############
##################################################################

# Load the dataset
Visvader_0068_ER_Total <- Read10X(data.dir = "/R/R_Visvader/Visvader_Input/ER/Visvader_0068_ER_Total/")

# Initialize the Seurat object with the raw (non-normalized data).
Visvader_0068_ER_Total <- CreateSeuratObject(counts = Visvader_0068_ER_Total, project = "Visvader_0068_ER_Total", min.cells = 3, min.features = 200)
Visvader_0068_ER_Total

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Visvader_0068_ER_Total[["percent.mt"]] <- PercentageFeatureSet(Visvader_0068_ER_Total, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Visvader_0068_ER_Total, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Visvader_0068_ER_Total, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Visvader_0068_ER_Total, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Visvader_0068_ER_Total <- subset(Visvader_0068_ER_Total, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Visvader_0068_ER_Total <- NormalizeData(Visvader_0068_ER_Total, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Visvader_0068_ER_Total <- FindVariableFeatures(Visvader_0068_ER_Total, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Visvader_0068_ER_Total)
Visvader_0068_ER_Total <- ScaleData(Visvader_0068_ER_Total, features = all.genes)

# Perform linear dimensional reduction
Visvader_0068_ER_Total <- RunPCA(Visvader_0068_ER_Total, features = VariableFeatures(object = Visvader_0068_ER_Total))

# Examine and visualize PCA results a few different ways
print(Visvader_0068_ER_Total[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Visvader_0068_ER_Total, dims = 1:2, reduction = "pca")
DimPlot(Visvader_0068_ER_Total, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Visvader_0068_ER_Total <- JackStraw(Visvader_0068_ER_Total, num.replicate = 100)
Visvader_0068_ER_Total <- ScoreJackStraw(Visvader_0068_ER_Total, dims = 1:20)
ElbowPlot(Visvader_0068_ER_Total)

# Cluster cells
Visvader_0068_ER_Total <- FindNeighbors(Visvader_0068_ER_Total, dims = 1:20)
Visvader_0068_ER_Total <- FindClusters(Visvader_0068_ER_Total, resolution = 0.5)

# Run non-linear dimensional reduction
Visvader_0068_ER_Total <- RunUMAP(Visvader_0068_ER_Total, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Visvader_0068_ER_Total, reduction = "umap")

# Save RDS
saveRDS(Visvader_0068_ER_Total, file = "/R/R_Visvader/Visvader_RDS/RDS_Total/Visvader_0068_ER_Total.rds")



##################################################################
########################### Visvader_0114_ER_Total ###############
##################################################################

# Load the dataset
Visvader_0114_ER_Total <- Read10X(data.dir = "/R/R_Visvader/Visvader_Input/ER/Visvader_0114_ER_Total/")

# Initialize the Seurat object with the raw (non-normalized data).
Visvader_0114_ER_Total <- CreateSeuratObject(counts = Visvader_0114_ER_Total, project = "Visvader_0114_ER_Total", min.cells = 3, min.features = 200)
Visvader_0114_ER_Total

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Visvader_0114_ER_Total[["percent.mt"]] <- PercentageFeatureSet(Visvader_0114_ER_Total, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Visvader_0114_ER_Total, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Visvader_0114_ER_Total, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Visvader_0114_ER_Total, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Visvader_0114_ER_Total <- subset(Visvader_0114_ER_Total, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Visvader_0114_ER_Total <- NormalizeData(Visvader_0114_ER_Total, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Visvader_0114_ER_Total <- FindVariableFeatures(Visvader_0114_ER_Total, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Visvader_0114_ER_Total)
Visvader_0114_ER_Total <- ScaleData(Visvader_0114_ER_Total, features = all.genes)

# Perform linear dimensional reduction
Visvader_0114_ER_Total <- RunPCA(Visvader_0114_ER_Total, features = VariableFeatures(object = Visvader_0114_ER_Total))

# Examine and visualize PCA results a few different ways
print(Visvader_0114_ER_Total[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Visvader_0114_ER_Total, dims = 1:2, reduction = "pca")
DimPlot(Visvader_0114_ER_Total, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Visvader_0114_ER_Total <- JackStraw(Visvader_0114_ER_Total, num.replicate = 100)
Visvader_0114_ER_Total <- ScoreJackStraw(Visvader_0114_ER_Total, dims = 1:20)
ElbowPlot(Visvader_0114_ER_Total)

# Cluster cells
Visvader_0114_ER_Total <- FindNeighbors(Visvader_0114_ER_Total, dims = 1:20)
Visvader_0114_ER_Total <- FindClusters(Visvader_0114_ER_Total, resolution = 0.5)

# Run non-linear dimensional reduction
Visvader_0114_ER_Total <- RunUMAP(Visvader_0114_ER_Total, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Visvader_0114_ER_Total, reduction = "umap")

# Save RDS
saveRDS(Visvader_0114_ER_Total, file = "/R/R_Visvader/Visvader_RDS/RDS_Total/Visvader_0114_ER_Total.rds")



##################################################################
########################### Visvader_0125_ER_Total ###############
##################################################################

# Load the dataset
Visvader_0125_ER_Total <- Read10X(data.dir = "/R/R_Visvader/Visvader_Input/ER/Visvader_0125_ER_Total/")

# Initialize the Seurat object with the raw (non-normalized data).
Visvader_0125_ER_Total <- CreateSeuratObject(counts = Visvader_0125_ER_Total, project = "Visvader_0125_ER_Total", min.cells = 3, min.features = 200)
Visvader_0125_ER_Total

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Visvader_0125_ER_Total[["percent.mt"]] <- PercentageFeatureSet(Visvader_0125_ER_Total, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Visvader_0125_ER_Total, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Visvader_0125_ER_Total, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Visvader_0125_ER_Total, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Visvader_0125_ER_Total <- subset(Visvader_0125_ER_Total, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Visvader_0125_ER_Total <- NormalizeData(Visvader_0125_ER_Total, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Visvader_0125_ER_Total <- FindVariableFeatures(Visvader_0125_ER_Total, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Visvader_0125_ER_Total)
Visvader_0125_ER_Total <- ScaleData(Visvader_0125_ER_Total, features = all.genes)

# Perform linear dimensional reduction
Visvader_0125_ER_Total <- RunPCA(Visvader_0125_ER_Total, features = VariableFeatures(object = Visvader_0125_ER_Total))

# Examine and visualize PCA results a few different ways
print(Visvader_0125_ER_Total[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Visvader_0125_ER_Total, dims = 1:2, reduction = "pca")
DimPlot(Visvader_0125_ER_Total, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Visvader_0125_ER_Total <- JackStraw(Visvader_0125_ER_Total, num.replicate = 100)
Visvader_0125_ER_Total <- ScoreJackStraw(Visvader_0125_ER_Total, dims = 1:20)
ElbowPlot(Visvader_0125_ER_Total)

# Cluster cells
Visvader_0125_ER_Total <- FindNeighbors(Visvader_0125_ER_Total, dims = 1:20)
Visvader_0125_ER_Total <- FindClusters(Visvader_0125_ER_Total, resolution = 0.5)

# Run non-linear dimensional reduction
Visvader_0125_ER_Total <- RunUMAP(Visvader_0125_ER_Total, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Visvader_0125_ER_Total, reduction = "umap")

# Save RDS
saveRDS(Visvader_0125_ER_Total, file = "/R/R_Visvader/Visvader_RDS/RDS_Total/Visvader_0125_ER_Total.rds")

##################################################################
########################### Visvader_0151_ER_Total ###############
##################################################################

# Load the dataset
Visvader_0151_ER_Total <- Read10X(data.dir = "/R/R_Visvader/Visvader_Input/ER/Visvader_0151_ER_Total/")

# Initialize the Seurat object with the raw (non-normalized data).
Visvader_0151_ER_Total <- CreateSeuratObject(counts = Visvader_0151_ER_Total, project = "Visvader_0151_ER_Total", min.cells = 3, min.features = 200)
Visvader_0151_ER_Total

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Visvader_0151_ER_Total[["percent.mt"]] <- PercentageFeatureSet(Visvader_0151_ER_Total, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Visvader_0151_ER_Total, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Visvader_0151_ER_Total, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Visvader_0151_ER_Total, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Visvader_0151_ER_Total <- subset(Visvader_0151_ER_Total, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Visvader_0151_ER_Total <- NormalizeData(Visvader_0151_ER_Total, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Visvader_0151_ER_Total <- FindVariableFeatures(Visvader_0151_ER_Total, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Visvader_0151_ER_Total)
Visvader_0151_ER_Total <- ScaleData(Visvader_0151_ER_Total, features = all.genes)

# Perform linear dimensional reduction
Visvader_0151_ER_Total <- RunPCA(Visvader_0151_ER_Total, features = VariableFeatures(object = Visvader_0151_ER_Total))

# Examine and visualize PCA results a few different ways
print(Visvader_0151_ER_Total[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Visvader_0151_ER_Total, dims = 1:2, reduction = "pca")
DimPlot(Visvader_0151_ER_Total, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Visvader_0151_ER_Total <- JackStraw(Visvader_0151_ER_Total, num.replicate = 100)
Visvader_0151_ER_Total <- ScoreJackStraw(Visvader_0151_ER_Total, dims = 1:20)
ElbowPlot(Visvader_0151_ER_Total)

# Cluster cells
Visvader_0151_ER_Total <- FindNeighbors(Visvader_0151_ER_Total, dims = 1:20)
Visvader_0151_ER_Total <- FindClusters(Visvader_0151_ER_Total, resolution = 0.5)

# Run non-linear dimensional reduction
Visvader_0151_ER_Total <- RunUMAP(Visvader_0151_ER_Total, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Visvader_0151_ER_Total, reduction = "umap")

# Save RDS
saveRDS(Visvader_0151_ER_Total, file = "/R/R_Visvader/Visvader_RDS/RDS_Total/Visvader_0151_ER_Total.rds")


##################################################################
########################### Visvader_0163_ER_Total ###############
##################################################################

# Load the dataset
Visvader_0163_ER_Total <- Read10X(data.dir = "/R/R_Visvader/Visvader_Input/ER/Visvader_0163_ER_Total/")

# Initialize the Seurat object with the raw (non-normalized data).
Visvader_0163_ER_Total <- CreateSeuratObject(counts = Visvader_0163_ER_Total, project = "Visvader_0163_ER_Total", min.cells = 3, min.features = 200)
Visvader_0163_ER_Total

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Visvader_0163_ER_Total[["percent.mt"]] <- PercentageFeatureSet(Visvader_0163_ER_Total, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Visvader_0163_ER_Total, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Visvader_0163_ER_Total, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Visvader_0163_ER_Total, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Visvader_0163_ER_Total <- subset(Visvader_0163_ER_Total, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Visvader_0163_ER_Total <- NormalizeData(Visvader_0163_ER_Total, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Visvader_0163_ER_Total <- FindVariableFeatures(Visvader_0163_ER_Total, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Visvader_0163_ER_Total)
Visvader_0163_ER_Total <- ScaleData(Visvader_0163_ER_Total, features = all.genes)

# Perform linear dimensional reduction
Visvader_0163_ER_Total <- RunPCA(Visvader_0163_ER_Total, features = VariableFeatures(object = Visvader_0163_ER_Total))

# Examine and visualize PCA results a few different ways
print(Visvader_0163_ER_Total[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Visvader_0163_ER_Total, dims = 1:2, reduction = "pca")
DimPlot(Visvader_0163_ER_Total, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Visvader_0163_ER_Total <- JackStraw(Visvader_0163_ER_Total, num.replicate = 100)
Visvader_0163_ER_Total <- ScoreJackStraw(Visvader_0163_ER_Total, dims = 1:20)
ElbowPlot(Visvader_0163_ER_Total)

# Cluster cells
Visvader_0163_ER_Total <- FindNeighbors(Visvader_0163_ER_Total, dims = 1:20)
Visvader_0163_ER_Total <- FindClusters(Visvader_0163_ER_Total, resolution = 0.5)

# Run non-linear dimensional reduction
Visvader_0163_ER_Total <- RunUMAP(Visvader_0163_ER_Total, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Visvader_0163_ER_Total, reduction = "umap")

# Save RDS
saveRDS(Visvader_0163_ER_Total, file = "/R/R_Visvader/Visvader_RDS/RDS_Total/Visvader_0163_ER_Total.rds")


##################################################################
########################### Visvader_0167_ER_Total ###############
##################################################################

# Load the dataset
Visvader_0167_ER_Total <- Read10X(data.dir = "/R/R_Visvader/Visvader_Input/ER/Visvader_0167_ER_Total/")

# Initialize the Seurat object with the raw (non-normalized data).
Visvader_0167_ER_Total <- CreateSeuratObject(counts = Visvader_0167_ER_Total, project = "Visvader_0167_ER_Total", min.cells = 3, min.features = 200)
Visvader_0167_ER_Total

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Visvader_0167_ER_Total[["percent.mt"]] <- PercentageFeatureSet(Visvader_0167_ER_Total, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Visvader_0167_ER_Total, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Visvader_0167_ER_Total, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Visvader_0167_ER_Total, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Visvader_0167_ER_Total <- subset(Visvader_0167_ER_Total, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Visvader_0167_ER_Total <- NormalizeData(Visvader_0167_ER_Total, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Visvader_0167_ER_Total <- FindVariableFeatures(Visvader_0167_ER_Total, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Visvader_0167_ER_Total)
Visvader_0167_ER_Total <- ScaleData(Visvader_0167_ER_Total, features = all.genes)

# Perform linear dimensional reduction
Visvader_0167_ER_Total <- RunPCA(Visvader_0167_ER_Total, features = VariableFeatures(object = Visvader_0167_ER_Total))

# Examine and visualize PCA results a few different ways
print(Visvader_0167_ER_Total[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Visvader_0167_ER_Total, dims = 1:2, reduction = "pca")
DimPlot(Visvader_0167_ER_Total, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Visvader_0167_ER_Total <- JackStraw(Visvader_0167_ER_Total, num.replicate = 100)
Visvader_0167_ER_Total <- ScoreJackStraw(Visvader_0167_ER_Total, dims = 1:20)
ElbowPlot(Visvader_0167_ER_Total)

# Cluster cells
Visvader_0167_ER_Total <- FindNeighbors(Visvader_0167_ER_Total, dims = 1:20)
Visvader_0167_ER_Total <- FindClusters(Visvader_0167_ER_Total, resolution = 0.5)

# Run non-linear dimensional reduction
Visvader_0167_ER_Total <- RunUMAP(Visvader_0167_ER_Total, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Visvader_0167_ER_Total, reduction = "umap")

# Save RDS
saveRDS(Visvader_0167_ER_Total, file = "/R/R_Visvader/Visvader_RDS/RDS_Total/Visvader_0167_ER_Total.rds")


##################################################################
########################### Visvader_0173_ER_Total ###############
##################################################################

# Load the dataset
Visvader_0173_ER_Total <- Read10X(data.dir = "/R/R_Visvader/Visvader_Input/ER/Visvader_0173_ER_Total/")

# Initialize the Seurat object with the raw (non-normalized data).
Visvader_0173_ER_Total <- CreateSeuratObject(counts = Visvader_0173_ER_Total, project = "Visvader_0173_ER_Total", min.cells = 3, min.features = 200)
Visvader_0173_ER_Total

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Visvader_0173_ER_Total[["percent.mt"]] <- PercentageFeatureSet(Visvader_0173_ER_Total, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Visvader_0173_ER_Total, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Visvader_0173_ER_Total, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Visvader_0173_ER_Total, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Visvader_0173_ER_Total <- subset(Visvader_0173_ER_Total, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Visvader_0173_ER_Total <- NormalizeData(Visvader_0173_ER_Total, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Visvader_0173_ER_Total <- FindVariableFeatures(Visvader_0173_ER_Total, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Visvader_0173_ER_Total)
Visvader_0173_ER_Total <- ScaleData(Visvader_0173_ER_Total, features = all.genes)

# Perform linear dimensional reduction
Visvader_0173_ER_Total <- RunPCA(Visvader_0173_ER_Total, features = VariableFeatures(object = Visvader_0173_ER_Total))

# Examine and visualize PCA results a few different ways
print(Visvader_0173_ER_Total[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Visvader_0173_ER_Total, dims = 1:2, reduction = "pca")
DimPlot(Visvader_0173_ER_Total, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Visvader_0173_ER_Total <- JackStraw(Visvader_0173_ER_Total, num.replicate = 100)
Visvader_0173_ER_Total <- ScoreJackStraw(Visvader_0173_ER_Total, dims = 1:20)
ElbowPlot(Visvader_0173_ER_Total)

# Cluster cells
Visvader_0173_ER_Total <- FindNeighbors(Visvader_0173_ER_Total, dims = 1:20)
Visvader_0173_ER_Total <- FindClusters(Visvader_0173_ER_Total, resolution = 0.5)

# Run non-linear dimensional reduction
Visvader_0173_ER_Total <- RunUMAP(Visvader_0173_ER_Total, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Visvader_0173_ER_Total, reduction = "umap")

# Save RDS
saveRDS(Visvader_0173_ER_Total, file = "/R/R_Visvader/Visvader_RDS/RDS_Total/Visvader_0173_ER_Total.rds")


##################################################################
########################### Visvader_0178_ER_Total ###############
##################################################################

# Load the dataset
Visvader_0178_ER_Total <- Read10X(data.dir = "/R/R_Visvader/Visvader_Input/ER/Visvader_0178_ER_Total/")

# Initialize the Seurat object with the raw (non-normalized data).
Visvader_0178_ER_Total <- CreateSeuratObject(counts = Visvader_0178_ER_Total, project = "Visvader_0178_ER_Total", min.cells = 3, min.features = 200)
Visvader_0178_ER_Total

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Visvader_0178_ER_Total[["percent.mt"]] <- PercentageFeatureSet(Visvader_0178_ER_Total, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Visvader_0178_ER_Total, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Visvader_0178_ER_Total, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Visvader_0178_ER_Total, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Visvader_0178_ER_Total <- subset(Visvader_0178_ER_Total, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Visvader_0178_ER_Total <- NormalizeData(Visvader_0178_ER_Total, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Visvader_0178_ER_Total <- FindVariableFeatures(Visvader_0178_ER_Total, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Visvader_0178_ER_Total)
Visvader_0178_ER_Total <- ScaleData(Visvader_0178_ER_Total, features = all.genes)

# Perform linear dimensional reduction
Visvader_0178_ER_Total <- RunPCA(Visvader_0178_ER_Total, features = VariableFeatures(object = Visvader_0178_ER_Total))

# Examine and visualize PCA results a few different ways
print(Visvader_0178_ER_Total[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Visvader_0178_ER_Total, dims = 1:2, reduction = "pca")
DimPlot(Visvader_0178_ER_Total, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Visvader_0178_ER_Total <- JackStraw(Visvader_0178_ER_Total, num.replicate = 100)
Visvader_0178_ER_Total <- ScoreJackStraw(Visvader_0178_ER_Total, dims = 1:20)
ElbowPlot(Visvader_0178_ER_Total)

# Cluster cells
Visvader_0178_ER_Total <- FindNeighbors(Visvader_0178_ER_Total, dims = 1:20)
Visvader_0178_ER_Total <- FindClusters(Visvader_0178_ER_Total, resolution = 0.5)

# Run non-linear dimensional reduction
Visvader_0178_ER_Total <- RunUMAP(Visvader_0178_ER_Total, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Visvader_0178_ER_Total, reduction = "umap")

# Save RDS
saveRDS(Visvader_0178_ER_Total, file = "/R/R_Visvader/Visvader_RDS/RDS_Total/Visvader_0178_ER_Total.rds")



##################################################################
########################### Visvader_0360_ER_Total ###############
##################################################################

# Load the dataset
Visvader_0360_ER_Total <- Read10X(data.dir = "/R/R_Visvader/Visvader_Input/ER/Visvader_0360_ER_Total/")

# Initialize the Seurat object with the raw (non-normalized data).
Visvader_0360_ER_Total <- CreateSeuratObject(counts = Visvader_0360_ER_Total, project = "Visvader_0360_ER_Total", min.cells = 3, min.features = 200)
Visvader_0360_ER_Total

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Visvader_0360_ER_Total[["percent.mt"]] <- PercentageFeatureSet(Visvader_0360_ER_Total, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Visvader_0360_ER_Total, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Visvader_0360_ER_Total, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Visvader_0360_ER_Total, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Visvader_0360_ER_Total <- subset(Visvader_0360_ER_Total, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Visvader_0360_ER_Total <- NormalizeData(Visvader_0360_ER_Total, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Visvader_0360_ER_Total <- FindVariableFeatures(Visvader_0360_ER_Total, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Visvader_0360_ER_Total)
Visvader_0360_ER_Total <- ScaleData(Visvader_0360_ER_Total, features = all.genes)

# Perform linear dimensional reduction
Visvader_0360_ER_Total <- RunPCA(Visvader_0360_ER_Total, features = VariableFeatures(object = Visvader_0360_ER_Total))

# Examine and visualize PCA results a few different ways
print(Visvader_0360_ER_Total[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Visvader_0360_ER_Total, dims = 1:2, reduction = "pca")
DimPlot(Visvader_0360_ER_Total, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Visvader_0360_ER_Total <- JackStraw(Visvader_0360_ER_Total, num.replicate = 100)
Visvader_0360_ER_Total <- ScoreJackStraw(Visvader_0360_ER_Total, dims = 1:20)
ElbowPlot(Visvader_0360_ER_Total)

# Cluster cells
Visvader_0360_ER_Total <- FindNeighbors(Visvader_0360_ER_Total, dims = 1:20)
Visvader_0360_ER_Total <- FindClusters(Visvader_0360_ER_Total, resolution = 0.5)

# Run non-linear dimensional reduction
Visvader_0360_ER_Total <- RunUMAP(Visvader_0360_ER_Total, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Visvader_0360_ER_Total, reduction = "umap")

# Save RDS
saveRDS(Visvader_0360_ER_Total, file = "/R/R_Visvader/Visvader_RDS/RDS_Total/Visvader_0360_ER_Total.rds")


##################################################################
########################### Visvader_0319_PR_Total ###############
##################################################################

# Load the dataset
Visvader_0319_PR_Total <- Read10X(data.dir = "/R/R_Visvader/Visvader_Input/PR/Visvader_0319_PR_Total/")

# Initialize the Seurat object with the raw (non-normalized data).
Visvader_0319_PR_Total <- CreateSeuratObject(counts = Visvader_0319_PR_Total, project = "Visvader_0319_PR_Total", min.cells = 3, min.features = 200)
Visvader_0319_PR_Total

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Visvader_0319_PR_Total[["percent.mt"]] <- PercentageFeatureSet(Visvader_0319_PR_Total, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Visvader_0319_PR_Total, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Visvader_0319_PR_Total, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Visvader_0319_PR_Total, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Visvader_0319_PR_Total <- subset(Visvader_0319_PR_Total, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Visvader_0319_PR_Total <- NormalizeData(Visvader_0319_PR_Total, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Visvader_0319_PR_Total <- FindVariableFeatures(Visvader_0319_PR_Total, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Visvader_0319_PR_Total)
Visvader_0319_PR_Total <- ScaleData(Visvader_0319_PR_Total, features = all.genes)

# Perform linear dimensional reduction
Visvader_0319_PR_Total <- RunPCA(Visvader_0319_PR_Total, features = VariableFeatures(object = Visvader_0319_PR_Total))

# Examine and visualize PCA results a few different ways
print(Visvader_0319_PR_Total[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Visvader_0319_PR_Total, dims = 1:2, reduction = "pca")
DimPlot(Visvader_0319_PR_Total, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Visvader_0319_PR_Total <- JackStraw(Visvader_0319_PR_Total, num.replicate = 100)
Visvader_0319_PR_Total <- ScoreJackStraw(Visvader_0319_PR_Total, dims = 1:20)
ElbowPlot(Visvader_0319_PR_Total)

# Cluster cells
Visvader_0319_PR_Total <- FindNeighbors(Visvader_0319_PR_Total, dims = 1:20)
Visvader_0319_PR_Total <- FindClusters(Visvader_0319_PR_Total, resolution = 0.5)

# Run non-linear dimensional reduction
Visvader_0319_PR_Total <- RunUMAP(Visvader_0319_PR_Total, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Visvader_0319_PR_Total, reduction = "umap")

# Save RDS
saveRDS(Visvader_0319_PR_Total, file = "/R/R_Visvader/Visvader_RDS/RDS_Total/Visvader_0319_PR_Total.rds")


#########################################################
############# This is HER2 Tumor Batch # ################
#########################################################

##################################################################
########################### Visvader_0031_HER2_Total #############
##################################################################

# Load the dataset
Visvader_0031_HER2_Total <- Read10X(data.dir = "/R/R_Visvader/Visvader_Input/HER2/Visvader_0031_HER2_Total/")

# Initialize the Seurat object with the raw (non-normalized data).
Visvader_0031_HER2_Total <- CreateSeuratObject(counts = Visvader_0031_HER2_Total, project = "Visvader_0031_HER2_Total", min.cells = 3, min.features = 200)
Visvader_0031_HER2_Total

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Visvader_0031_HER2_Total[["percent.mt"]] <- PercentageFeatureSet(Visvader_0031_HER2_Total, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Visvader_0031_HER2_Total, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Visvader_0031_HER2_Total, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Visvader_0031_HER2_Total, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Visvader_0031_HER2_Total <- subset(Visvader_0031_HER2_Total, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Visvader_0031_HER2_Total <- NormalizeData(Visvader_0031_HER2_Total, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Visvader_0031_HER2_Total <- FindVariableFeatures(Visvader_0031_HER2_Total, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Visvader_0031_HER2_Total)
Visvader_0031_HER2_Total <- ScaleData(Visvader_0031_HER2_Total, features = all.genes)

# Perform linear dimensional reduction
Visvader_0031_HER2_Total <- RunPCA(Visvader_0031_HER2_Total, features = VariableFeatures(object = Visvader_0031_HER2_Total))

# Examine and visualize PCA results a few different ways
print(Visvader_0031_HER2_Total[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Visvader_0031_HER2_Total, dims = 1:2, reduction = "pca")
DimPlot(Visvader_0031_HER2_Total, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Visvader_0031_HER2_Total <- JackStraw(Visvader_0031_HER2_Total, num.replicate = 100)
Visvader_0031_HER2_Total <- ScoreJackStraw(Visvader_0031_HER2_Total, dims = 1:20)
ElbowPlot(Visvader_0031_HER2_Total)

# Cluster cells
Visvader_0031_HER2_Total <- FindNeighbors(Visvader_0031_HER2_Total, dims = 1:20)
Visvader_0031_HER2_Total <- FindClusters(Visvader_0031_HER2_Total, resolution = 0.5)

# Run non-linear dimensional reduction
Visvader_0031_HER2_Total <- RunUMAP(Visvader_0031_HER2_Total, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Visvader_0031_HER2_Total, reduction = "umap")

# Save RDS
saveRDS(Visvader_0031_HER2_Total, file = "/R/R_Visvader/Visvader_RDS/RDS_Total/Visvader_0031_HER2_Total.rds")




##################################################################
########################### Visvader_0069_HER2_Total #############
##################################################################

# Load the dataset
Visvader_0069_HER2_Total <- Read10X(data.dir = "/R/R_Visvader/Visvader_Input/HER2/Visvader_0069_HER2_Total/")

# Initialize the Seurat object with the raw (non-normalized data).
Visvader_0069_HER2_Total <- CreateSeuratObject(counts = Visvader_0069_HER2_Total, project = "Visvader_0069_HER2_Total", min.cells = 3, min.features = 200)
Visvader_0069_HER2_Total

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Visvader_0069_HER2_Total[["percent.mt"]] <- PercentageFeatureSet(Visvader_0069_HER2_Total, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Visvader_0069_HER2_Total, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Visvader_0069_HER2_Total, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Visvader_0069_HER2_Total, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Visvader_0069_HER2_Total <- subset(Visvader_0069_HER2_Total, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Visvader_0069_HER2_Total <- NormalizeData(Visvader_0069_HER2_Total, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Visvader_0069_HER2_Total <- FindVariableFeatures(Visvader_0069_HER2_Total, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Visvader_0069_HER2_Total)
Visvader_0069_HER2_Total <- ScaleData(Visvader_0069_HER2_Total, features = all.genes)

# Perform linear dimensional reduction
Visvader_0069_HER2_Total <- RunPCA(Visvader_0069_HER2_Total, features = VariableFeatures(object = Visvader_0069_HER2_Total))

# Examine and visualize PCA results a few different ways
print(Visvader_0069_HER2_Total[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Visvader_0069_HER2_Total, dims = 1:2, reduction = "pca")
DimPlot(Visvader_0069_HER2_Total, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Visvader_0069_HER2_Total <- JackStraw(Visvader_0069_HER2_Total, num.replicate = 100)
Visvader_0069_HER2_Total <- ScoreJackStraw(Visvader_0069_HER2_Total, dims = 1:20)
ElbowPlot(Visvader_0069_HER2_Total)

# Cluster cells
Visvader_0069_HER2_Total <- FindNeighbors(Visvader_0069_HER2_Total, dims = 1:20)
Visvader_0069_HER2_Total <- FindClusters(Visvader_0069_HER2_Total, resolution = 0.5)

# Run non-linear dimensional reduction
Visvader_0069_HER2_Total <- RunUMAP(Visvader_0069_HER2_Total, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Visvader_0069_HER2_Total, reduction = "umap")

# Save RDS
saveRDS(Visvader_0069_HER2_Total, file = "/R/R_Visvader/Visvader_RDS/RDS_Total/Visvader_0069_HER2_Total.rds")




##################################################################
########################### Visvader_0161_HER2_Total #############
##################################################################

# Load the dataset
Visvader_0161_HER2_Total <- Read10X(data.dir = "/R/R_Visvader/Visvader_Input/HER2/Visvader_0161_HER2_Total/")

# Initialize the Seurat object with the raw (non-normalized data).
Visvader_0161_HER2_Total <- CreateSeuratObject(counts = Visvader_0161_HER2_Total, project = "Visvader_0161_HER2_Total", min.cells = 3, min.features = 200)
Visvader_0161_HER2_Total

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Visvader_0161_HER2_Total[["percent.mt"]] <- PercentageFeatureSet(Visvader_0161_HER2_Total, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Visvader_0161_HER2_Total, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Visvader_0161_HER2_Total, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Visvader_0161_HER2_Total, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Visvader_0161_HER2_Total <- subset(Visvader_0161_HER2_Total, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Visvader_0161_HER2_Total <- NormalizeData(Visvader_0161_HER2_Total, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Visvader_0161_HER2_Total <- FindVariableFeatures(Visvader_0161_HER2_Total, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Visvader_0161_HER2_Total)
Visvader_0161_HER2_Total <- ScaleData(Visvader_0161_HER2_Total, features = all.genes)

# Perform linear dimensional reduction
Visvader_0161_HER2_Total <- RunPCA(Visvader_0161_HER2_Total, features = VariableFeatures(object = Visvader_0161_HER2_Total))

# Examine and visualize PCA results a few different ways
print(Visvader_0161_HER2_Total[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Visvader_0161_HER2_Total, dims = 1:2, reduction = "pca")
DimPlot(Visvader_0161_HER2_Total, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Visvader_0161_HER2_Total <- JackStraw(Visvader_0161_HER2_Total, num.replicate = 100)
Visvader_0161_HER2_Total <- ScoreJackStraw(Visvader_0161_HER2_Total, dims = 1:20)
ElbowPlot(Visvader_0161_HER2_Total)

# Cluster cells
Visvader_0161_HER2_Total <- FindNeighbors(Visvader_0161_HER2_Total, dims = 1:20)
Visvader_0161_HER2_Total <- FindClusters(Visvader_0161_HER2_Total, resolution = 0.5)

# Run non-linear dimensional reduction
Visvader_0161_HER2_Total <- RunUMAP(Visvader_0161_HER2_Total, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Visvader_0161_HER2_Total, reduction = "umap")

# Save RDS
saveRDS(Visvader_0161_HER2_Total, file = "/R/R_Visvader/Visvader_RDS/RDS_Total/Visvader_0161_HER2_Total.rds")




##################################################################
########################### Visvader_0176_HER2_Total #############
##################################################################

# Load the dataset
Visvader_0176_HER2_Total <- Read10X(data.dir = "/R/R_Visvader/Visvader_Input/HER2/Visvader_0176_HER2_Total/")

# Initialize the Seurat object with the raw (non-normalized data).
Visvader_0176_HER2_Total <- CreateSeuratObject(counts = Visvader_0176_HER2_Total, project = "Visvader_0176_HER2_Total", min.cells = 3, min.features = 200)
Visvader_0176_HER2_Total

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Visvader_0176_HER2_Total[["percent.mt"]] <- PercentageFeatureSet(Visvader_0176_HER2_Total, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Visvader_0176_HER2_Total, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Visvader_0176_HER2_Total, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Visvader_0176_HER2_Total, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Visvader_0176_HER2_Total <- subset(Visvader_0176_HER2_Total, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Visvader_0176_HER2_Total <- NormalizeData(Visvader_0176_HER2_Total, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Visvader_0176_HER2_Total <- FindVariableFeatures(Visvader_0176_HER2_Total, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Visvader_0176_HER2_Total)
Visvader_0176_HER2_Total <- ScaleData(Visvader_0176_HER2_Total, features = all.genes)

# Perform linear dimensional reduction
Visvader_0176_HER2_Total <- RunPCA(Visvader_0176_HER2_Total, features = VariableFeatures(object = Visvader_0176_HER2_Total))

# Examine and visualize PCA results a few different ways
print(Visvader_0176_HER2_Total[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Visvader_0176_HER2_Total, dims = 1:2, reduction = "pca")
DimPlot(Visvader_0176_HER2_Total, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Visvader_0176_HER2_Total <- JackStraw(Visvader_0176_HER2_Total, num.replicate = 100)
Visvader_0176_HER2_Total <- ScoreJackStraw(Visvader_0176_HER2_Total, dims = 1:20)
ElbowPlot(Visvader_0176_HER2_Total)

# Cluster cells
Visvader_0176_HER2_Total <- FindNeighbors(Visvader_0176_HER2_Total, dims = 1:20)
Visvader_0176_HER2_Total <- FindClusters(Visvader_0176_HER2_Total, resolution = 0.5)

# Run non-linear dimensional reduction
Visvader_0176_HER2_Total <- RunUMAP(Visvader_0176_HER2_Total, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Visvader_0176_HER2_Total, reduction = "umap")

# Save RDS
saveRDS(Visvader_0176_HER2_Total, file = "/R/R_Visvader/Visvader_RDS/RDS_Total/Visvader_0176_HER2_Total.rds")




##################################################################
########################### Visvader_0308_HER2_Total #############
##################################################################

# Load the dataset
Visvader_0308_HER2_Total <- Read10X(data.dir = "/R/R_Visvader/Visvader_Input/HER2/Visvader_0308_HER2_Total/")

# Initialize the Seurat object with the raw (non-normalized data).
Visvader_0308_HER2_Total <- CreateSeuratObject(counts = Visvader_0308_HER2_Total, project = "Visvader_0308_HER2_Total", min.cells = 3, min.features = 200)
Visvader_0308_HER2_Total

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Visvader_0308_HER2_Total[["percent.mt"]] <- PercentageFeatureSet(Visvader_0308_HER2_Total, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Visvader_0308_HER2_Total, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Visvader_0308_HER2_Total, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Visvader_0308_HER2_Total, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Visvader_0308_HER2_Total <- subset(Visvader_0308_HER2_Total, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Visvader_0308_HER2_Total <- NormalizeData(Visvader_0308_HER2_Total, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Visvader_0308_HER2_Total <- FindVariableFeatures(Visvader_0308_HER2_Total, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Visvader_0308_HER2_Total)
Visvader_0308_HER2_Total <- ScaleData(Visvader_0308_HER2_Total, features = all.genes)

# Perform linear dimensional reduction
Visvader_0308_HER2_Total <- RunPCA(Visvader_0308_HER2_Total, features = VariableFeatures(object = Visvader_0308_HER2_Total))

# Examine and visualize PCA results a few different ways
print(Visvader_0308_HER2_Total[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Visvader_0308_HER2_Total, dims = 1:2, reduction = "pca")
DimPlot(Visvader_0308_HER2_Total, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Visvader_0308_HER2_Total <- JackStraw(Visvader_0308_HER2_Total, num.replicate = 100)
Visvader_0308_HER2_Total <- ScoreJackStraw(Visvader_0308_HER2_Total, dims = 1:20)
ElbowPlot(Visvader_0308_HER2_Total)

# Cluster cells
Visvader_0308_HER2_Total <- FindNeighbors(Visvader_0308_HER2_Total, dims = 1:20)
Visvader_0308_HER2_Total <- FindClusters(Visvader_0308_HER2_Total, resolution = 0.5)

# Run non-linear dimensional reduction
Visvader_0308_HER2_Total <- RunUMAP(Visvader_0308_HER2_Total, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Visvader_0308_HER2_Total, reduction = "umap")

# Save RDS
saveRDS(Visvader_0308_HER2_Total, file = "/R/R_Visvader/Visvader_RDS/RDS_Total/Visvader_0308_HER2_Total.rds")




##################################################################
########################### Visvader_0337_HER2_Total #############
##################################################################

# Load the dataset
Visvader_0337_HER2_Total <- Read10X(data.dir = "/R/R_Visvader/Visvader_Input/HER2/Visvader_0337_HER2_Total/")

# Initialize the Seurat object with the raw (non-normalized data).
Visvader_0337_HER2_Total <- CreateSeuratObject(counts = Visvader_0337_HER2_Total, project = "Visvader_0337_HER2_Total", min.cells = 3, min.features = 200)
Visvader_0337_HER2_Total

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Visvader_0337_HER2_Total[["percent.mt"]] <- PercentageFeatureSet(Visvader_0337_HER2_Total, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Visvader_0337_HER2_Total, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Visvader_0337_HER2_Total, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Visvader_0337_HER2_Total, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Visvader_0337_HER2_Total <- subset(Visvader_0337_HER2_Total, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Visvader_0337_HER2_Total <- NormalizeData(Visvader_0337_HER2_Total, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Visvader_0337_HER2_Total <- FindVariableFeatures(Visvader_0337_HER2_Total, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Visvader_0337_HER2_Total)
Visvader_0337_HER2_Total <- ScaleData(Visvader_0337_HER2_Total, features = all.genes)

# Perform linear dimensional reduction
Visvader_0337_HER2_Total <- RunPCA(Visvader_0337_HER2_Total, features = VariableFeatures(object = Visvader_0337_HER2_Total))

# Examine and visualize PCA results a few different ways
print(Visvader_0337_HER2_Total[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Visvader_0337_HER2_Total, dims = 1:2, reduction = "pca")
DimPlot(Visvader_0337_HER2_Total, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Visvader_0337_HER2_Total <- JackStraw(Visvader_0337_HER2_Total, num.replicate = 100)
Visvader_0337_HER2_Total <- ScoreJackStraw(Visvader_0337_HER2_Total, dims = 1:20)
ElbowPlot(Visvader_0337_HER2_Total)

# Cluster cells
Visvader_0337_HER2_Total <- FindNeighbors(Visvader_0337_HER2_Total, dims = 1:20)
Visvader_0337_HER2_Total <- FindClusters(Visvader_0337_HER2_Total, resolution = 0.5)

# Run non-linear dimensional reduction
Visvader_0337_HER2_Total <- RunUMAP(Visvader_0337_HER2_Total, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Visvader_0337_HER2_Total, reduction = "umap")

# Save RDS
saveRDS(Visvader_0337_HER2_Total, file = "/R/R_Visvader/Visvader_RDS/RDS_Total/Visvader_0337_HER2_Total.rds")


#########################################################
############# This is TNBC Tumor Batch # ################
#########################################################

##################################################################
########################### Visvader_106_TNBC_Total ##############
##################################################################

# Load the dataset
Visvader_106_TNBC_Total <- Read10X(data.dir = "/R/R_Visvader/Visvader_Input/TNBC/Visvader_106_TNBC_Total/")

# Initialize the Seurat object with the raw (non-normalized data).
Visvader_106_TNBC_Total <- CreateSeuratObject(counts = Visvader_106_TNBC_Total, project = "Visvader_106_TNBC_Total", min.cells = 3, min.features = 200)
Visvader_106_TNBC_Total

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Visvader_106_TNBC_Total[["percent.mt"]] <- PercentageFeatureSet(Visvader_106_TNBC_Total, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Visvader_106_TNBC_Total, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Visvader_106_TNBC_Total, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Visvader_106_TNBC_Total, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Visvader_106_TNBC_Total <- subset(Visvader_106_TNBC_Total, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Visvader_106_TNBC_Total <- NormalizeData(Visvader_106_TNBC_Total, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Visvader_106_TNBC_Total <- FindVariableFeatures(Visvader_106_TNBC_Total, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Visvader_106_TNBC_Total)
Visvader_106_TNBC_Total <- ScaleData(Visvader_106_TNBC_Total, features = all.genes)

# Perform linear dimensional reduction
Visvader_106_TNBC_Total <- RunPCA(Visvader_106_TNBC_Total, features = VariableFeatures(object = Visvader_106_TNBC_Total))

# Examine and visualize PCA results a few different ways
print(Visvader_106_TNBC_Total[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Visvader_106_TNBC_Total, dims = 1:2, reduction = "pca")
DimPlot(Visvader_106_TNBC_Total, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Visvader_106_TNBC_Total <- JackStraw(Visvader_106_TNBC_Total, num.replicate = 100)
Visvader_106_TNBC_Total <- ScoreJackStraw(Visvader_106_TNBC_Total, dims = 1:20)
ElbowPlot(Visvader_106_TNBC_Total)

# Cluster cells
Visvader_106_TNBC_Total <- FindNeighbors(Visvader_106_TNBC_Total, dims = 1:20)
Visvader_106_TNBC_Total <- FindClusters(Visvader_106_TNBC_Total, resolution = 0.5)

# Run non-linear dimensional reduction
Visvader_106_TNBC_Total <- RunUMAP(Visvader_106_TNBC_Total, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Visvader_106_TNBC_Total, reduction = "umap")

# Save RDS
saveRDS(Visvader_106_TNBC_Total, file = "/R/R_Visvader/Visvader_RDS/RDS_Total/Visvader_106_TNBC_Total.rds")






##################################################################
########################### Visvader_114_TNBC_Total ##############
##################################################################

# Load the dataset
Visvader_114_TNBC_Total <- Read10X(data.dir = "/R/R_Visvader/Visvader_Input/TNBC/Visvader_114_TNBC_Total/")

# Initialize the Seurat object with the raw (non-normalized data).
Visvader_114_TNBC_Total <- CreateSeuratObject(counts = Visvader_114_TNBC_Total, project = "Visvader_114_TNBC_Total", min.cells = 3, min.features = 200)
Visvader_114_TNBC_Total

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Visvader_114_TNBC_Total[["percent.mt"]] <- PercentageFeatureSet(Visvader_114_TNBC_Total, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Visvader_114_TNBC_Total, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Visvader_114_TNBC_Total, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Visvader_114_TNBC_Total, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Visvader_114_TNBC_Total <- subset(Visvader_114_TNBC_Total, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Visvader_114_TNBC_Total <- NormalizeData(Visvader_114_TNBC_Total, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Visvader_114_TNBC_Total <- FindVariableFeatures(Visvader_114_TNBC_Total, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Visvader_114_TNBC_Total)
Visvader_114_TNBC_Total <- ScaleData(Visvader_114_TNBC_Total, features = all.genes)

# Perform linear dimensional reduction
Visvader_114_TNBC_Total <- RunPCA(Visvader_114_TNBC_Total, features = VariableFeatures(object = Visvader_114_TNBC_Total))

# Examine and visualize PCA results a few different ways
print(Visvader_114_TNBC_Total[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Visvader_114_TNBC_Total, dims = 1:2, reduction = "pca")
DimPlot(Visvader_114_TNBC_Total, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Visvader_114_TNBC_Total <- JackStraw(Visvader_114_TNBC_Total, num.replicate = 100)
Visvader_114_TNBC_Total <- ScoreJackStraw(Visvader_114_TNBC_Total, dims = 1:20)
ElbowPlot(Visvader_114_TNBC_Total)

# Cluster cells
Visvader_114_TNBC_Total <- FindNeighbors(Visvader_114_TNBC_Total, dims = 1:20)
Visvader_114_TNBC_Total <- FindClusters(Visvader_114_TNBC_Total, resolution = 0.5)

# Run non-linear dimensional reduction
Visvader_114_TNBC_Total <- RunUMAP(Visvader_114_TNBC_Total, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Visvader_114_TNBC_Total, reduction = "umap")

# Save RDS
saveRDS(Visvader_114_TNBC_Total, file = "/R/R_Visvader/Visvader_RDS/RDS_Total/Visvader_114_TNBC_Total.rds")






##################################################################
########################### Visvader_126_TNBC_Total ##############
##################################################################

# Load the dataset
Visvader_126_TNBC_Total <- Read10X(data.dir = "/R/R_Visvader/Visvader_Input/TNBC/Visvader_126_TNBC_Total/")

# Initialize the Seurat object with the raw (non-normalized data).
Visvader_126_TNBC_Total <- CreateSeuratObject(counts = Visvader_126_TNBC_Total, project = "Visvader_126_TNBC_Total", min.cells = 3, min.features = 200)
Visvader_126_TNBC_Total

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Visvader_126_TNBC_Total[["percent.mt"]] <- PercentageFeatureSet(Visvader_126_TNBC_Total, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Visvader_126_TNBC_Total, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Visvader_126_TNBC_Total, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Visvader_126_TNBC_Total, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Visvader_126_TNBC_Total <- subset(Visvader_126_TNBC_Total, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)

# Normalize data
Visvader_126_TNBC_Total <- NormalizeData(Visvader_126_TNBC_Total, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Visvader_126_TNBC_Total <- FindVariableFeatures(Visvader_126_TNBC_Total, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Visvader_126_TNBC_Total)
Visvader_126_TNBC_Total <- ScaleData(Visvader_126_TNBC_Total, features = all.genes)

# Perform linear dimensional reduction
Visvader_126_TNBC_Total <- RunPCA(Visvader_126_TNBC_Total, features = VariableFeatures(object = Visvader_126_TNBC_Total))

# Examine and visualize PCA results a few different ways
print(Visvader_126_TNBC_Total[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Visvader_126_TNBC_Total, dims = 1:2, reduction = "pca")
DimPlot(Visvader_126_TNBC_Total, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Visvader_126_TNBC_Total <- JackStraw(Visvader_126_TNBC_Total, num.replicate = 100)
Visvader_126_TNBC_Total <- ScoreJackStraw(Visvader_126_TNBC_Total, dims = 1:20)
ElbowPlot(Visvader_126_TNBC_Total)

# Cluster cells
Visvader_126_TNBC_Total <- FindNeighbors(Visvader_126_TNBC_Total, dims = 1:20)
Visvader_126_TNBC_Total <- FindClusters(Visvader_126_TNBC_Total, resolution = 0.5)

# Run non-linear dimensional reduction
Visvader_126_TNBC_Total <- RunUMAP(Visvader_126_TNBC_Total, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Visvader_126_TNBC_Total, reduction = "umap")

# Save RDS
saveRDS(Visvader_126_TNBC_Total, file = "/R/R_Visvader/Visvader_RDS/RDS_Total/Visvader_126_TNBC_Total.rds")





##################################################################
########################### Visvader_0131_TNBC_Total #############
##################################################################

# Load the dataset
Visvader_0131_TNBC_Total <- Read10X(data.dir = "/R/R_Visvader/Visvader_Input/TNBC/Visvader_0131_TNBC_Total/")

# Initialize the Seurat object with the raw (non-normalized data).
Visvader_0131_TNBC_Total <- CreateSeuratObject(counts = Visvader_0131_TNBC_Total, project = "Visvader_0131_TNBC_Total", min.cells = 3, min.features = 200)
Visvader_0131_TNBC_Total

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Visvader_0131_TNBC_Total[["percent.mt"]] <- PercentageFeatureSet(Visvader_0131_TNBC_Total, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Visvader_0131_TNBC_Total, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Visvader_0131_TNBC_Total, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Visvader_0131_TNBC_Total, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Visvader_0131_TNBC_Total <- subset(Visvader_0131_TNBC_Total, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Visvader_0131_TNBC_Total <- NormalizeData(Visvader_0131_TNBC_Total, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Visvader_0131_TNBC_Total <- FindVariableFeatures(Visvader_0131_TNBC_Total, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Visvader_0131_TNBC_Total)
Visvader_0131_TNBC_Total <- ScaleData(Visvader_0131_TNBC_Total, features = all.genes)

# Perform linear dimensional reduction
Visvader_0131_TNBC_Total <- RunPCA(Visvader_0131_TNBC_Total, features = VariableFeatures(object = Visvader_0131_TNBC_Total))

# Examine and visualize PCA results a few different ways
print(Visvader_0131_TNBC_Total[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Visvader_0131_TNBC_Total, dims = 1:2, reduction = "pca")
DimPlot(Visvader_0131_TNBC_Total, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Visvader_0131_TNBC_Total <- JackStraw(Visvader_0131_TNBC_Total, num.replicate = 100)
Visvader_0131_TNBC_Total <- ScoreJackStraw(Visvader_0131_TNBC_Total, dims = 1:20)
ElbowPlot(Visvader_0131_TNBC_Total)

# Cluster cells
Visvader_0131_TNBC_Total <- FindNeighbors(Visvader_0131_TNBC_Total, dims = 1:20)
Visvader_0131_TNBC_Total <- FindClusters(Visvader_0131_TNBC_Total, resolution = 0.5)

# Run non-linear dimensional reduction
Visvader_0131_TNBC_Total <- RunUMAP(Visvader_0131_TNBC_Total, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Visvader_0131_TNBC_Total, reduction = "umap")

# Save RDS
saveRDS(Visvader_0131_TNBC_Total, file = "/R/R_Visvader/Visvader_RDS/RDS_Total/Visvader_0131_TNBC_Total.rds")


##################################################################
########################### Visvader_135_TNBC_Total ##############
##################################################################

# Load the dataset
Visvader_135_TNBC_Total <- Read10X(data.dir = "/R/R_Visvader/Visvader_Input/TNBC/Visvader_135_TNBC_Total/")

# Initialize the Seurat object with the raw (non-normalized data).
Visvader_135_TNBC_Total <- CreateSeuratObject(counts = Visvader_135_TNBC_Total, project = "Visvader_135_TNBC_Total", min.cells = 3, min.features = 200)
Visvader_135_TNBC_Total

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Visvader_135_TNBC_Total[["percent.mt"]] <- PercentageFeatureSet(Visvader_135_TNBC_Total, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Visvader_135_TNBC_Total, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Visvader_135_TNBC_Total, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Visvader_135_TNBC_Total, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Visvader_135_TNBC_Total <- subset(Visvader_135_TNBC_Total, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Visvader_135_TNBC_Total <- NormalizeData(Visvader_135_TNBC_Total, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Visvader_135_TNBC_Total <- FindVariableFeatures(Visvader_135_TNBC_Total, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Visvader_135_TNBC_Total)
Visvader_135_TNBC_Total <- ScaleData(Visvader_135_TNBC_Total, features = all.genes)

# Perform linear dimensional reduction
Visvader_135_TNBC_Total <- RunPCA(Visvader_135_TNBC_Total, features = VariableFeatures(object = Visvader_135_TNBC_Total))

# Examine and visualize PCA results a few different ways
print(Visvader_135_TNBC_Total[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Visvader_135_TNBC_Total, dims = 1:2, reduction = "pca")
DimPlot(Visvader_135_TNBC_Total, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Visvader_135_TNBC_Total <- JackStraw(Visvader_135_TNBC_Total, num.replicate = 100)
Visvader_135_TNBC_Total <- ScoreJackStraw(Visvader_135_TNBC_Total, dims = 1:20)
ElbowPlot(Visvader_135_TNBC_Total)

# Cluster cells
Visvader_135_TNBC_Total <- FindNeighbors(Visvader_135_TNBC_Total, dims = 1:20)
Visvader_135_TNBC_Total <- FindClusters(Visvader_135_TNBC_Total, resolution = 0.5)

# Run non-linear dimensional reduction
Visvader_135_TNBC_Total <- RunUMAP(Visvader_135_TNBC_Total, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Visvader_135_TNBC_Total, reduction = "umap")

# Save RDS
saveRDS(Visvader_135_TNBC_Total, file = "/R/R_Visvader/Visvader_RDS/RDS_Total/Visvader_135_TNBC_Total.rds")




##################################################################
########################### Visvader_0177_TNBC_Total #############
##################################################################

# Load the dataset
Visvader_0177_TNBC_Total <- Read10X(data.dir = "/R/R_Visvader/Visvader_Input/TNBC/Visvader_0177_TNBC_Total/")

# Initialize the Seurat object with the raw (non-normalized data).
Visvader_0177_TNBC_Total <- CreateSeuratObject(counts = Visvader_0177_TNBC_Total, project = "Visvader_0177_TNBC_Total", min.cells = 3, min.features = 200)
Visvader_0177_TNBC_Total

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Visvader_0177_TNBC_Total[["percent.mt"]] <- PercentageFeatureSet(Visvader_0177_TNBC_Total, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Visvader_0177_TNBC_Total, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Visvader_0177_TNBC_Total, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Visvader_0177_TNBC_Total, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Visvader_0177_TNBC_Total <- subset(Visvader_0177_TNBC_Total, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
Visvader_0177_TNBC_Total <- NormalizeData(Visvader_0177_TNBC_Total, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Visvader_0177_TNBC_Total <- FindVariableFeatures(Visvader_0177_TNBC_Total, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Visvader_0177_TNBC_Total)
Visvader_0177_TNBC_Total <- ScaleData(Visvader_0177_TNBC_Total, features = all.genes)

# Perform linear dimensional reduction
Visvader_0177_TNBC_Total <- RunPCA(Visvader_0177_TNBC_Total, features = VariableFeatures(object = Visvader_0177_TNBC_Total))

# Examine and visualize PCA results a few different ways
print(Visvader_0177_TNBC_Total[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Visvader_0177_TNBC_Total, dims = 1:2, reduction = "pca")
DimPlot(Visvader_0177_TNBC_Total, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Visvader_0177_TNBC_Total <- JackStraw(Visvader_0177_TNBC_Total, num.replicate = 100)
Visvader_0177_TNBC_Total <- ScoreJackStraw(Visvader_0177_TNBC_Total, dims = 1:20)
ElbowPlot(Visvader_0177_TNBC_Total)

# Cluster cells
Visvader_0177_TNBC_Total <- FindNeighbors(Visvader_0177_TNBC_Total, dims = 1:20)
Visvader_0177_TNBC_Total <- FindClusters(Visvader_0177_TNBC_Total, resolution = 0.5)

# Run non-linear dimensional reduction
Visvader_0177_TNBC_Total <- RunUMAP(Visvader_0177_TNBC_Total, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Visvader_0177_TNBC_Total, reduction = "umap")

# Save RDS
saveRDS(Visvader_0177_TNBC_Total, file = "/R/R_Visvader/Visvader_RDS/RDS_Total/Visvader_0177_TNBC_Total.rds")




##################################################################
########################### Visvader_0554_TNBC_Total #############
##################################################################

# Load the dataset
Visvader_0554_TNBC_Total <- Read10X(data.dir = "/R/R_Visvader/Visvader_Input/TNBC/Visvader_0554_TNBC_Total/")

# Initialize the Seurat object with the raw (non-normalized data).
Visvader_0554_TNBC_Total <- CreateSeuratObject(counts = Visvader_0554_TNBC_Total, project = "Visvader_0554_TNBC_Total", min.cells = 3, min.features = 200)
Visvader_0554_TNBC_Total

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Visvader_0554_TNBC_Total[["percent.mt"]] <- PercentageFeatureSet(Visvader_0554_TNBC_Total, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Visvader_0554_TNBC_Total, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Visvader_0554_TNBC_Total, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Visvader_0554_TNBC_Total, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Visvader_0554_TNBC_Total <- subset(Visvader_0554_TNBC_Total, subset = nFeature_RNA > 200 & nFeature_RNA < 75000 & percent.mt < 10)

# Normalize data
Visvader_0554_TNBC_Total <- NormalizeData(Visvader_0554_TNBC_Total, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Visvader_0554_TNBC_Total <- FindVariableFeatures(Visvader_0554_TNBC_Total, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Visvader_0554_TNBC_Total)
Visvader_0554_TNBC_Total <- ScaleData(Visvader_0554_TNBC_Total, features = all.genes)

# Perform linear dimensional reduction
Visvader_0554_TNBC_Total <- RunPCA(Visvader_0554_TNBC_Total, features = VariableFeatures(object = Visvader_0554_TNBC_Total))

# Examine and visualize PCA results a few different ways
print(Visvader_0554_TNBC_Total[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Visvader_0554_TNBC_Total, dims = 1:2, reduction = "pca")
DimPlot(Visvader_0554_TNBC_Total, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Visvader_0554_TNBC_Total <- JackStraw(Visvader_0554_TNBC_Total, num.replicate = 100)
Visvader_0554_TNBC_Total <- ScoreJackStraw(Visvader_0554_TNBC_Total, dims = 1:20)
ElbowPlot(Visvader_0554_TNBC_Total)

# Cluster cells
Visvader_0554_TNBC_Total <- FindNeighbors(Visvader_0554_TNBC_Total, dims = 1:20)
Visvader_0554_TNBC_Total <- FindClusters(Visvader_0554_TNBC_Total, resolution = 0.5)

# Run non-linear dimensional reduction
Visvader_0554_TNBC_Total <- RunUMAP(Visvader_0554_TNBC_Total, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Visvader_0554_TNBC_Total, reduction = "umap")

# Save RDS
saveRDS(Visvader_0554_TNBC_Total, file = "/R/R_Visvader/Visvader_RDS/RDS_Total/Visvader_0554_TNBC_Total.rds")




##################################################################
########################### Visvader_4031_TNBC_Total #############
##################################################################

# Load the dataset
Visvader_4031_TNBC_Total <- Read10X(data.dir = "/R/R_Visvader/Visvader_Input/TNBC/Visvader_4031_TNBC_Total/")

# Initialize the Seurat object with the raw (non-normalized data).
Visvader_4031_TNBC_Total <- CreateSeuratObject(counts = Visvader_4031_TNBC_Total, project = "Visvader_4031_TNBC_Total", min.cells = 3, min.features = 200)
Visvader_4031_TNBC_Total

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Visvader_4031_TNBC_Total[["percent.mt"]] <- PercentageFeatureSet(Visvader_4031_TNBC_Total, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Visvader_4031_TNBC_Total, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Visvader_4031_TNBC_Total, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Visvader_4031_TNBC_Total, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Visvader_4031_TNBC_Total <- subset(Visvader_4031_TNBC_Total, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)

# Normalize data
Visvader_4031_TNBC_Total <- NormalizeData(Visvader_4031_TNBC_Total, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Visvader_4031_TNBC_Total <- FindVariableFeatures(Visvader_4031_TNBC_Total, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Visvader_4031_TNBC_Total)
Visvader_4031_TNBC_Total <- ScaleData(Visvader_4031_TNBC_Total, features = all.genes)

# Perform linear dimensional reduction
Visvader_4031_TNBC_Total <- RunPCA(Visvader_4031_TNBC_Total, features = VariableFeatures(object = Visvader_4031_TNBC_Total))

# Examine and visualize PCA results a few different ways
print(Visvader_4031_TNBC_Total[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Visvader_4031_TNBC_Total, dims = 1:2, reduction = "pca")
DimPlot(Visvader_4031_TNBC_Total, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Visvader_4031_TNBC_Total <- JackStraw(Visvader_4031_TNBC_Total, num.replicate = 100)
Visvader_4031_TNBC_Total <- ScoreJackStraw(Visvader_4031_TNBC_Total, dims = 1:20)
ElbowPlot(Visvader_4031_TNBC_Total)

# Cluster cells
Visvader_4031_TNBC_Total <- FindNeighbors(Visvader_4031_TNBC_Total, dims = 1:20)
Visvader_4031_TNBC_Total <- FindClusters(Visvader_4031_TNBC_Total, resolution = 0.5)

# Run non-linear dimensional reduction
Visvader_4031_TNBC_Total <- RunUMAP(Visvader_4031_TNBC_Total, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Visvader_4031_TNBC_Total, reduction = "umap")

# Save RDS
saveRDS(Visvader_4031_TNBC_Total, file = "/R/R_Visvader/Visvader_RDS/RDS_Total/Visvader_4031_TNBC_Total.rds")



#####################################################################
########### Step 2: Run DoubletFinder and Subset Singlets ########### 
#####################################################################

# Load Libraries
library(Seurat)
library(DoubletFinder)

# Assuming 5% doublet rate per sample based on average number of sequenced cells
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-

# Patient Samples
############################ ER
# Visvader_0001_ER_Total
# Visvader_0025_ER_Total
# Visvader_0029_7C_ER_Total
# Visvader_0029_9C_ER_Total
# Visvader_0032_ER_Total
# Visvader_0040_ER_Total
# Visvader_0042_ER_Total
# Visvader_0043_ER_Total
# Visvader_0056_ER_Total
# Visvader_0064_ER_Total
# Visvader_0068_ER_Total
# Visvader_0114_ER_Total
# Visvader_0125_ER_Total
# Visvader_0151_ER_Total
# Visvader_0163_ER_Total
# Visvader_0167_ER_Total
# Visvader_0173_ER_Total
# Visvader_0178_ER_Total
# Visvader_0360_ER_Total
############################  PR
# Visvader_0319_PR_Total
############################  HER2
# Visvader_0031_HER2_Total
# Visvader_0069_HER2_Total
# Visvader_0161_HER2_Total
# Visvader_0176_HER2_Total
# Visvader_0308_HER2_Total
# Visvader_0337_HER2_Total
############################  TNBC
# Visvader_106_TNBC_Total
# Visvader_114_TNBC_Total
# Visvader_126_TNBC_Total
# Visvader_0131_TNBC_Total
# Visvader_135_TNBC_Total
# Visvader_0177_TNBC_Total
# Visvader_0554_TNBC_Total
# Visvader_4031_TNBC_Total



################################################################################################
########################                   ER-pos Tumors                ########################
################################################################################################

################################################################################################
########################              Visvader_0001_ER_Total             #######################
################################################################################################

# Load Data
Visvader_0001_ER_Total <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Total/Visvader_0001_ER_Total.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Visvader_0001_ER_Total <- paramSweep(Visvader_0001_ER_Total, PCs = 1:10, sct = FALSE)
sweep.stats.Visvader_0001_ER_Total <- summarizeSweep(sweep.res.Visvader_0001_ER_Total, GT = FALSE)
bcmvn_Visvader_0001_ER_Total <- find.pK(sweep.stats.Visvader_0001_ER_Total)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Visvader_0001_ER_Total[which.max(bcmvn_Visvader_0001_ER_Total$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Visvader_0001_ER_Total@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 5% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Visvader dataset = 5k per sample, will use doublet rate for 6,000 cells to err on the side of caution
nExp_poi <- round(0.05*nrow(Visvader_0001_ER_Total@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Visvader_0001_ER_Total <- doubletFinder(Visvader_0001_ER_Total, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Visvader_0001_ER_Total_Quant <- (Visvader_0001_ER_Total@meta.data$DF.classification == "Singlet")
Visvader_0001_ER_Total_Quant_Singlets <- length(Visvader_0001_ER_Total_Quant[Visvader_0001_ER_Total_Quant== TRUE])
Visvader_0001_ER_Total_Quant_Doublets <- length(Visvader_0001_ER_Total_Quant[Visvader_0001_ER_Total_Quant== FALSE])
Visvader_0001_ER_Total_Quant_Doublets_Percent <- Visvader_0001_ER_Total_Quant_Doublets / (Visvader_0001_ER_Total_Quant_Doublets + Visvader_0001_ER_Total_Quant_Singlets) * 100
Visvader_0001_ER_Total_Quant <- as.data.frame(c(Visvader_0001_ER_Total_Quant_Singlets, Visvader_0001_ER_Total_Quant_Doublets, Visvader_0001_ER_Total_Quant_Doublets_Percent))
colnames(Visvader_0001_ER_Total_Quant) <- c("Visvader_0001_ER_Total_Quant")
rownames(Visvader_0001_ER_Total_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Visvader_0001_ER_Total_Quant, file = "/R/R_Visvader/Visvader_Doublet_Tables/Visvader_0001_ER_Total_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Visvader_0001_ER_Total_Singlets <- subset(Visvader_0001_ER_Total, cells=rownames(Visvader_0001_ER_Total@meta.data)[which(Visvader_0001_ER_Total@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Visvader_0001_ER_Total_Singlets, file = "/R/R_Visvader/Visvader_RDS/RDS_Total_Singlets/Visvader_0001_ER_Total_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Visvader_0001_ER_Total)
rm(Visvader_0001_ER_Total_Singlets)
rm(Visvader_0001_ER_Total_Quant)
rm(Visvader_0001_ER_Total_Quant_Singlets)
rm(Visvader_0001_ER_Total_Quant_Doublets)
rm(Visvader_0001_ER_Total_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Visvader_0001_ER_Total)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Visvader_0001_ER_Total)
rm(sweep.stats.Visvader_0001_ER_Total)
gc()




################################################################################################
########################              Visvader_0025_ER_Total             #######################
################################################################################################

# Load Data
Visvader_0025_ER_Total <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Total/Visvader_0025_ER_Total.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Visvader_0025_ER_Total <- paramSweep(Visvader_0025_ER_Total, PCs = 1:10, sct = FALSE)
sweep.stats.Visvader_0025_ER_Total <- summarizeSweep(sweep.res.Visvader_0025_ER_Total, GT = FALSE)
bcmvn_Visvader_0025_ER_Total <- find.pK(sweep.stats.Visvader_0025_ER_Total)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Visvader_0025_ER_Total[which.max(bcmvn_Visvader_0025_ER_Total$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Visvader_0025_ER_Total@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 5% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Visvader dataset = 5k per sample, will use doublet rate for 6,000 cells to err on the side of caution
nExp_poi <- round(0.05*nrow(Visvader_0025_ER_Total@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Visvader_0025_ER_Total <- doubletFinder(Visvader_0025_ER_Total, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Visvader_0025_ER_Total_Quant <- (Visvader_0025_ER_Total@meta.data$DF.classification == "Singlet")
Visvader_0025_ER_Total_Quant_Singlets <- length(Visvader_0025_ER_Total_Quant[Visvader_0025_ER_Total_Quant== TRUE])
Visvader_0025_ER_Total_Quant_Doublets <- length(Visvader_0025_ER_Total_Quant[Visvader_0025_ER_Total_Quant== FALSE])
Visvader_0025_ER_Total_Quant_Doublets_Percent <- Visvader_0025_ER_Total_Quant_Doublets /  (Visvader_0025_ER_Total_Quant_Doublets + Visvader_0025_ER_Total_Quant_Singlets) * 100
Visvader_0025_ER_Total_Quant <- as.data.frame(c(Visvader_0025_ER_Total_Quant_Singlets, Visvader_0025_ER_Total_Quant_Doublets, Visvader_0025_ER_Total_Quant_Doublets_Percent))
colnames(Visvader_0025_ER_Total_Quant) <- c("Visvader_0025_ER_Total_Quant")
rownames(Visvader_0025_ER_Total_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Visvader_0025_ER_Total_Quant, file = "/R/R_Visvader/Visvader_Doublet_Tables/Visvader_0025_ER_Total_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Visvader_0025_ER_Total_Singlets <- subset(Visvader_0025_ER_Total, cells=rownames(Visvader_0025_ER_Total@meta.data)[which(Visvader_0025_ER_Total@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Visvader_0025_ER_Total_Singlets, file = "/R/R_Visvader/Visvader_RDS/RDS_Total_Singlets/Visvader_0025_ER_Total_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Visvader_0025_ER_Total)
rm(Visvader_0025_ER_Total_Singlets)
rm(Visvader_0025_ER_Total_Quant)
rm(Visvader_0025_ER_Total_Quant_Singlets)
rm(Visvader_0025_ER_Total_Quant_Doublets)
rm(Visvader_0025_ER_Total_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Visvader_0025_ER_Total)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Visvader_0025_ER_Total)
rm(sweep.stats.Visvader_0025_ER_Total)
gc()




################################################################################################
########################              Visvader_0029_7C_ER_Total             ####################
################################################################################################

# Load Data
Visvader_0029_7C_ER_Total <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Total/Visvader_0029_7C_ER_Total.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Visvader_0029_7C_ER_Total <- paramSweep(Visvader_0029_7C_ER_Total, PCs = 1:10, sct = FALSE)
sweep.stats.Visvader_0029_7C_ER_Total <- summarizeSweep(sweep.res.Visvader_0029_7C_ER_Total, GT = FALSE)
bcmvn_Visvader_0029_7C_ER_Total <- find.pK(sweep.stats.Visvader_0029_7C_ER_Total)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Visvader_0029_7C_ER_Total[which.max(bcmvn_Visvader_0029_7C_ER_Total$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Visvader_0029_7C_ER_Total@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 5% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Visvader dataset = 5k per sample, will use doublet rate for 6,000 cells to err on the side of caution
nExp_poi <- round(0.05*nrow(Visvader_0029_7C_ER_Total@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Visvader_0029_7C_ER_Total <- doubletFinder(Visvader_0029_7C_ER_Total, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Visvader_0029_7C_ER_Total_Quant <- (Visvader_0029_7C_ER_Total@meta.data$DF.classification == "Singlet")
Visvader_0029_7C_ER_Total_Quant_Singlets <- length(Visvader_0029_7C_ER_Total_Quant[Visvader_0029_7C_ER_Total_Quant== TRUE])
Visvader_0029_7C_ER_Total_Quant_Doublets <- length(Visvader_0029_7C_ER_Total_Quant[Visvader_0029_7C_ER_Total_Quant== FALSE])
Visvader_0029_7C_ER_Total_Quant_Doublets_Percent <- Visvader_0029_7C_ER_Total_Quant_Doublets /  (Visvader_0029_7C_ER_Total_Quant_Doublets + Visvader_0029_7C_ER_Total_Quant_Singlets) * 100
Visvader_0029_7C_ER_Total_Quant <- as.data.frame(c(Visvader_0029_7C_ER_Total_Quant_Singlets, Visvader_0029_7C_ER_Total_Quant_Doublets, Visvader_0029_7C_ER_Total_Quant_Doublets_Percent))
colnames(Visvader_0029_7C_ER_Total_Quant) <- c("Visvader_0029_7C_ER_Total_Quant")
rownames(Visvader_0029_7C_ER_Total_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Visvader_0029_7C_ER_Total_Quant, file = "/R/R_Visvader/Visvader_Doublet_Tables/Visvader_0029_7C_ER_Total_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Visvader_0029_7C_ER_Total_Singlets <- subset(Visvader_0029_7C_ER_Total, cells=rownames(Visvader_0029_7C_ER_Total@meta.data)[which(Visvader_0029_7C_ER_Total@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Visvader_0029_7C_ER_Total_Singlets, file = "/R/R_Visvader/Visvader_RDS/RDS_Total_Singlets/Visvader_0029_7C_ER_Total_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Visvader_0029_7C_ER_Total)
rm(Visvader_0029_7C_ER_Total_Singlets)
rm(Visvader_0029_7C_ER_Total_Quant)
rm(Visvader_0029_7C_ER_Total_Quant_Singlets)
rm(Visvader_0029_7C_ER_Total_Quant_Doublets)
rm(Visvader_0029_7C_ER_Total_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Visvader_0029_7C_ER_Total)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Visvader_0029_7C_ER_Total)
rm(sweep.stats.Visvader_0029_7C_ER_Total)
gc()




################################################################################################
########################              Visvader_0029_9C_ER_Total             ####################
################################################################################################

# Load Data
Visvader_0029_9C_ER_Total <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Total/Visvader_0029_9C_ER_Total.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Visvader_0029_9C_ER_Total <- paramSweep(Visvader_0029_9C_ER_Total, PCs = 1:10, sct = FALSE)
sweep.stats.Visvader_0029_9C_ER_Total <- summarizeSweep(sweep.res.Visvader_0029_9C_ER_Total, GT = FALSE)
bcmvn_Visvader_0029_9C_ER_Total <- find.pK(sweep.stats.Visvader_0029_9C_ER_Total)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Visvader_0029_9C_ER_Total[which.max(bcmvn_Visvader_0029_9C_ER_Total$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Visvader_0029_9C_ER_Total@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 5% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Visvader dataset = 5k per sample, will use doublet rate for 6,000 cells to err on the side of caution
nExp_poi <- round(0.05*nrow(Visvader_0029_9C_ER_Total@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Visvader_0029_9C_ER_Total <- doubletFinder(Visvader_0029_9C_ER_Total, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Visvader_0029_9C_ER_Total_Quant <- (Visvader_0029_9C_ER_Total@meta.data$DF.classification == "Singlet")
Visvader_0029_9C_ER_Total_Quant_Singlets <- length(Visvader_0029_9C_ER_Total_Quant[Visvader_0029_9C_ER_Total_Quant== TRUE])
Visvader_0029_9C_ER_Total_Quant_Doublets <- length(Visvader_0029_9C_ER_Total_Quant[Visvader_0029_9C_ER_Total_Quant== FALSE])
Visvader_0029_9C_ER_Total_Quant_Doublets_Percent <- Visvader_0029_9C_ER_Total_Quant_Doublets /  (Visvader_0029_9C_ER_Total_Quant_Doublets + Visvader_0029_9C_ER_Total_Quant_Singlets) * 100
Visvader_0029_9C_ER_Total_Quant <- as.data.frame(c(Visvader_0029_9C_ER_Total_Quant_Singlets, Visvader_0029_9C_ER_Total_Quant_Doublets, Visvader_0029_9C_ER_Total_Quant_Doublets_Percent))
colnames(Visvader_0029_9C_ER_Total_Quant) <- c("Visvader_0029_9C_ER_Total_Quant")
rownames(Visvader_0029_9C_ER_Total_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Visvader_0029_9C_ER_Total_Quant, file = "/R/R_Visvader/Visvader_Doublet_Tables/Visvader_0029_9C_ER_Total_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Visvader_0029_9C_ER_Total_Singlets <- subset(Visvader_0029_9C_ER_Total, cells=rownames(Visvader_0029_9C_ER_Total@meta.data)[which(Visvader_0029_9C_ER_Total@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Visvader_0029_9C_ER_Total_Singlets, file = "/R/R_Visvader/Visvader_RDS/RDS_Total_Singlets/Visvader_0029_9C_ER_Total_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Visvader_0029_9C_ER_Total)
rm(Visvader_0029_9C_ER_Total_Singlets)
rm(Visvader_0029_9C_ER_Total_Quant)
rm(Visvader_0029_9C_ER_Total_Quant_Singlets)
rm(Visvader_0029_9C_ER_Total_Quant_Doublets)
rm(Visvader_0029_9C_ER_Total_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Visvader_0029_9C_ER_Total)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Visvader_0029_9C_ER_Total)
rm(sweep.stats.Visvader_0029_9C_ER_Total)
gc()




################################################################################################
########################              Visvader_0032_ER_Total             #######################
################################################################################################

# Load Data
Visvader_0032_ER_Total <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Total/Visvader_0032_ER_Total.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Visvader_0032_ER_Total <- paramSweep(Visvader_0032_ER_Total, PCs = 1:10, sct = FALSE)
sweep.stats.Visvader_0032_ER_Total <- summarizeSweep(sweep.res.Visvader_0032_ER_Total, GT = FALSE)
bcmvn_Visvader_0032_ER_Total <- find.pK(sweep.stats.Visvader_0032_ER_Total)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Visvader_0032_ER_Total[which.max(bcmvn_Visvader_0032_ER_Total$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Visvader_0032_ER_Total@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 5% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Visvader dataset = 5k per sample, will use doublet rate for 6,000 cells to err on the side of caution
nExp_poi <- round(0.05*nrow(Visvader_0032_ER_Total@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Visvader_0032_ER_Total <- doubletFinder(Visvader_0032_ER_Total, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Visvader_0032_ER_Total_Quant <- (Visvader_0032_ER_Total@meta.data$DF.classification == "Singlet")
Visvader_0032_ER_Total_Quant_Singlets <- length(Visvader_0032_ER_Total_Quant[Visvader_0032_ER_Total_Quant== TRUE])
Visvader_0032_ER_Total_Quant_Doublets <- length(Visvader_0032_ER_Total_Quant[Visvader_0032_ER_Total_Quant== FALSE])
Visvader_0032_ER_Total_Quant_Doublets_Percent <- Visvader_0032_ER_Total_Quant_Doublets /  (Visvader_0032_ER_Total_Quant_Doublets + Visvader_0032_ER_Total_Quant_Singlets) * 100
Visvader_0032_ER_Total_Quant <- as.data.frame(c(Visvader_0032_ER_Total_Quant_Singlets, Visvader_0032_ER_Total_Quant_Doublets, Visvader_0032_ER_Total_Quant_Doublets_Percent))
colnames(Visvader_0032_ER_Total_Quant) <- c("Visvader_0032_ER_Total_Quant")
rownames(Visvader_0032_ER_Total_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Visvader_0032_ER_Total_Quant, file = "/R/R_Visvader/Visvader_Doublet_Tables/Visvader_0032_ER_Total_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Visvader_0032_ER_Total_Singlets <- subset(Visvader_0032_ER_Total, cells=rownames(Visvader_0032_ER_Total@meta.data)[which(Visvader_0032_ER_Total@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Visvader_0032_ER_Total_Singlets, file = "/R/R_Visvader/Visvader_RDS/RDS_Total_Singlets/Visvader_0032_ER_Total_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Visvader_0032_ER_Total)
rm(Visvader_0032_ER_Total_Singlets)
rm(Visvader_0032_ER_Total_Quant)
rm(Visvader_0032_ER_Total_Quant_Singlets)
rm(Visvader_0032_ER_Total_Quant_Doublets)
rm(Visvader_0032_ER_Total_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Visvader_0032_ER_Total)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Visvader_0032_ER_Total)
rm(sweep.stats.Visvader_0032_ER_Total)
gc()




################################################################################################
########################              Visvader_0040_ER_Total             #######################
################################################################################################

# Load Data
Visvader_0040_ER_Total <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Total/Visvader_0040_ER_Total.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Visvader_0040_ER_Total <- paramSweep(Visvader_0040_ER_Total, PCs = 1:10, sct = FALSE)
sweep.stats.Visvader_0040_ER_Total <- summarizeSweep(sweep.res.Visvader_0040_ER_Total, GT = FALSE)
bcmvn_Visvader_0040_ER_Total <- find.pK(sweep.stats.Visvader_0040_ER_Total)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Visvader_0040_ER_Total[which.max(bcmvn_Visvader_0040_ER_Total$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Visvader_0040_ER_Total@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 5% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Visvader dataset = 5k per sample, will use doublet rate for 6,000 cells to err on the side of caution
nExp_poi <- round(0.05*nrow(Visvader_0040_ER_Total@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Visvader_0040_ER_Total <- doubletFinder(Visvader_0040_ER_Total, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Visvader_0040_ER_Total_Quant <- (Visvader_0040_ER_Total@meta.data$DF.classification == "Singlet")
Visvader_0040_ER_Total_Quant_Singlets <- length(Visvader_0040_ER_Total_Quant[Visvader_0040_ER_Total_Quant== TRUE])
Visvader_0040_ER_Total_Quant_Doublets <- length(Visvader_0040_ER_Total_Quant[Visvader_0040_ER_Total_Quant== FALSE])
Visvader_0040_ER_Total_Quant_Doublets_Percent <- Visvader_0040_ER_Total_Quant_Doublets /  (Visvader_0040_ER_Total_Quant_Doublets + Visvader_0040_ER_Total_Quant_Singlets) * 100
Visvader_0040_ER_Total_Quant <- as.data.frame(c(Visvader_0040_ER_Total_Quant_Singlets, Visvader_0040_ER_Total_Quant_Doublets, Visvader_0040_ER_Total_Quant_Doublets_Percent))
colnames(Visvader_0040_ER_Total_Quant) <- c("Visvader_0040_ER_Total_Quant")
rownames(Visvader_0040_ER_Total_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Visvader_0040_ER_Total_Quant, file = "/R/R_Visvader/Visvader_Doublet_Tables/Visvader_0040_ER_Total_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Visvader_0040_ER_Total_Singlets <- subset(Visvader_0040_ER_Total, cells=rownames(Visvader_0040_ER_Total@meta.data)[which(Visvader_0040_ER_Total@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Visvader_0040_ER_Total_Singlets, file = "/R/R_Visvader/Visvader_RDS/RDS_Total_Singlets/Visvader_0040_ER_Total_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Visvader_0040_ER_Total)
rm(Visvader_0040_ER_Total_Singlets)
rm(Visvader_0040_ER_Total_Quant)
rm(Visvader_0040_ER_Total_Quant_Singlets)
rm(Visvader_0040_ER_Total_Quant_Doublets)
rm(Visvader_0040_ER_Total_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Visvader_0040_ER_Total)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Visvader_0040_ER_Total)
rm(sweep.stats.Visvader_0040_ER_Total)
gc()




################################################################################################
########################              Visvader_0042_ER_Total             #######################
################################################################################################

# Load Data
Visvader_0042_ER_Total <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Total/Visvader_0042_ER_Total.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Visvader_0042_ER_Total <- paramSweep(Visvader_0042_ER_Total, PCs = 1:10, sct = FALSE)
sweep.stats.Visvader_0042_ER_Total <- summarizeSweep(sweep.res.Visvader_0042_ER_Total, GT = FALSE)
bcmvn_Visvader_0042_ER_Total <- find.pK(sweep.stats.Visvader_0042_ER_Total)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Visvader_0042_ER_Total[which.max(bcmvn_Visvader_0042_ER_Total$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Visvader_0042_ER_Total@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 5% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Visvader dataset = 5k per sample, will use doublet rate for 6,000 cells to err on the side of caution
nExp_poi <- round(0.05*nrow(Visvader_0042_ER_Total@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Visvader_0042_ER_Total <- doubletFinder(Visvader_0042_ER_Total, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Visvader_0042_ER_Total_Quant <- (Visvader_0042_ER_Total@meta.data$DF.classification == "Singlet")
Visvader_0042_ER_Total_Quant_Singlets <- length(Visvader_0042_ER_Total_Quant[Visvader_0042_ER_Total_Quant== TRUE])
Visvader_0042_ER_Total_Quant_Doublets <- length(Visvader_0042_ER_Total_Quant[Visvader_0042_ER_Total_Quant== FALSE])
Visvader_0042_ER_Total_Quant_Doublets_Percent <- Visvader_0042_ER_Total_Quant_Doublets /  (Visvader_0042_ER_Total_Quant_Doublets + Visvader_0042_ER_Total_Quant_Singlets) * 100
Visvader_0042_ER_Total_Quant <- as.data.frame(c(Visvader_0042_ER_Total_Quant_Singlets, Visvader_0042_ER_Total_Quant_Doublets, Visvader_0042_ER_Total_Quant_Doublets_Percent))
colnames(Visvader_0042_ER_Total_Quant) <- c("Visvader_0042_ER_Total_Quant")
rownames(Visvader_0042_ER_Total_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Visvader_0042_ER_Total_Quant, file = "/R/R_Visvader/Visvader_Doublet_Tables/Visvader_0042_ER_Total_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Visvader_0042_ER_Total_Singlets <- subset(Visvader_0042_ER_Total, cells=rownames(Visvader_0042_ER_Total@meta.data)[which(Visvader_0042_ER_Total@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Visvader_0042_ER_Total_Singlets, file = "/R/R_Visvader/Visvader_RDS/RDS_Total_Singlets/Visvader_0042_ER_Total_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Visvader_0042_ER_Total)
rm(Visvader_0042_ER_Total_Singlets)
rm(Visvader_0042_ER_Total_Quant)
rm(Visvader_0042_ER_Total_Quant_Singlets)
rm(Visvader_0042_ER_Total_Quant_Doublets)
rm(Visvader_0042_ER_Total_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Visvader_0042_ER_Total)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Visvader_0042_ER_Total)
rm(sweep.stats.Visvader_0042_ER_Total)
gc()




################################################################################################
########################              Visvader_0043_ER_Total             #######################
################################################################################################

# Load Data
Visvader_0043_ER_Total <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Total/Visvader_0043_ER_Total.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Visvader_0043_ER_Total <- paramSweep(Visvader_0043_ER_Total, PCs = 1:10, sct = FALSE)
sweep.stats.Visvader_0043_ER_Total <- summarizeSweep(sweep.res.Visvader_0043_ER_Total, GT = FALSE)
bcmvn_Visvader_0043_ER_Total <- find.pK(sweep.stats.Visvader_0043_ER_Total)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Visvader_0043_ER_Total[which.max(bcmvn_Visvader_0043_ER_Total$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Visvader_0043_ER_Total@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 5% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Visvader dataset = 5k per sample, will use doublet rate for 6,000 cells to err on the side of caution
nExp_poi <- round(0.05*nrow(Visvader_0043_ER_Total@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Visvader_0043_ER_Total <- doubletFinder(Visvader_0043_ER_Total, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Visvader_0043_ER_Total_Quant <- (Visvader_0043_ER_Total@meta.data$DF.classification == "Singlet")
Visvader_0043_ER_Total_Quant_Singlets <- length(Visvader_0043_ER_Total_Quant[Visvader_0043_ER_Total_Quant== TRUE])
Visvader_0043_ER_Total_Quant_Doublets <- length(Visvader_0043_ER_Total_Quant[Visvader_0043_ER_Total_Quant== FALSE])
Visvader_0043_ER_Total_Quant_Doublets_Percent <- Visvader_0043_ER_Total_Quant_Doublets /  (Visvader_0043_ER_Total_Quant_Doublets + Visvader_0043_ER_Total_Quant_Singlets) * 100
Visvader_0043_ER_Total_Quant <- as.data.frame(c(Visvader_0043_ER_Total_Quant_Singlets, Visvader_0043_ER_Total_Quant_Doublets, Visvader_0043_ER_Total_Quant_Doublets_Percent))
colnames(Visvader_0043_ER_Total_Quant) <- c("Visvader_0043_ER_Total_Quant")
rownames(Visvader_0043_ER_Total_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Visvader_0043_ER_Total_Quant, file = "/R/R_Visvader/Visvader_Doublet_Tables/Visvader_0043_ER_Total_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Visvader_0043_ER_Total_Singlets <- subset(Visvader_0043_ER_Total, cells=rownames(Visvader_0043_ER_Total@meta.data)[which(Visvader_0043_ER_Total@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Visvader_0043_ER_Total_Singlets, file = "/R/R_Visvader/Visvader_RDS/RDS_Total_Singlets/Visvader_0043_ER_Total_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Visvader_0043_ER_Total)
rm(Visvader_0043_ER_Total_Singlets)
rm(Visvader_0043_ER_Total_Quant)
rm(Visvader_0043_ER_Total_Quant_Singlets)
rm(Visvader_0043_ER_Total_Quant_Doublets)
rm(Visvader_0043_ER_Total_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Visvader_0043_ER_Total)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Visvader_0043_ER_Total)
rm(sweep.stats.Visvader_0043_ER_Total)
gc()




################################################################################################
########################              Visvader_0056_ER_Total             #######################
################################################################################################

# Load Data
Visvader_0056_ER_Total <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Total/Visvader_0056_ER_Total.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Visvader_0056_ER_Total <- paramSweep(Visvader_0056_ER_Total, PCs = 1:10, sct = FALSE)
sweep.stats.Visvader_0056_ER_Total <- summarizeSweep(sweep.res.Visvader_0056_ER_Total, GT = FALSE)
bcmvn_Visvader_0056_ER_Total <- find.pK(sweep.stats.Visvader_0056_ER_Total)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Visvader_0056_ER_Total[which.max(bcmvn_Visvader_0056_ER_Total$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Visvader_0056_ER_Total@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 5% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Visvader dataset = 5k per sample, will use doublet rate for 6,000 cells to err on the side of caution
nExp_poi <- round(0.05*nrow(Visvader_0056_ER_Total@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Visvader_0056_ER_Total <- doubletFinder(Visvader_0056_ER_Total, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Visvader_0056_ER_Total_Quant <- (Visvader_0056_ER_Total@meta.data$DF.classification == "Singlet")
Visvader_0056_ER_Total_Quant_Singlets <- length(Visvader_0056_ER_Total_Quant[Visvader_0056_ER_Total_Quant== TRUE])
Visvader_0056_ER_Total_Quant_Doublets <- length(Visvader_0056_ER_Total_Quant[Visvader_0056_ER_Total_Quant== FALSE])
Visvader_0056_ER_Total_Quant_Doublets_Percent <- Visvader_0056_ER_Total_Quant_Doublets /  (Visvader_0056_ER_Total_Quant_Doublets + Visvader_0056_ER_Total_Quant_Singlets) * 100
Visvader_0056_ER_Total_Quant <- as.data.frame(c(Visvader_0056_ER_Total_Quant_Singlets, Visvader_0056_ER_Total_Quant_Doublets, Visvader_0056_ER_Total_Quant_Doublets_Percent))
colnames(Visvader_0056_ER_Total_Quant) <- c("Visvader_0056_ER_Total_Quant")
rownames(Visvader_0056_ER_Total_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Visvader_0056_ER_Total_Quant, file = "/R/R_Visvader/Visvader_Doublet_Tables/Visvader_0056_ER_Total_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Visvader_0056_ER_Total_Singlets <- subset(Visvader_0056_ER_Total, cells=rownames(Visvader_0056_ER_Total@meta.data)[which(Visvader_0056_ER_Total@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Visvader_0056_ER_Total_Singlets, file = "/R/R_Visvader/Visvader_RDS/RDS_Total_Singlets/Visvader_0056_ER_Total_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Visvader_0056_ER_Total)
rm(Visvader_0056_ER_Total_Singlets)
rm(Visvader_0056_ER_Total_Quant)
rm(Visvader_0056_ER_Total_Quant_Singlets)
rm(Visvader_0056_ER_Total_Quant_Doublets)
rm(Visvader_0056_ER_Total_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Visvader_0056_ER_Total)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Visvader_0056_ER_Total)
rm(sweep.stats.Visvader_0056_ER_Total)
gc()




################################################################################################
########################              Visvader_0064_ER_Total             #######################
################################################################################################

# Load Data
Visvader_0064_ER_Total <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Total/Visvader_0064_ER_Total.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Visvader_0064_ER_Total <- paramSweep(Visvader_0064_ER_Total, PCs = 1:10, sct = FALSE)
sweep.stats.Visvader_0064_ER_Total <- summarizeSweep(sweep.res.Visvader_0064_ER_Total, GT = FALSE)
bcmvn_Visvader_0064_ER_Total <- find.pK(sweep.stats.Visvader_0064_ER_Total)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Visvader_0064_ER_Total[which.max(bcmvn_Visvader_0064_ER_Total$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Visvader_0064_ER_Total@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 5% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Visvader dataset = 5k per sample, will use doublet rate for 6,000 cells to err on the side of caution
nExp_poi <- round(0.05*nrow(Visvader_0064_ER_Total@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Visvader_0064_ER_Total <- doubletFinder(Visvader_0064_ER_Total, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Visvader_0064_ER_Total_Quant <- (Visvader_0064_ER_Total@meta.data$DF.classification == "Singlet")
Visvader_0064_ER_Total_Quant_Singlets <- length(Visvader_0064_ER_Total_Quant[Visvader_0064_ER_Total_Quant== TRUE])
Visvader_0064_ER_Total_Quant_Doublets <- length(Visvader_0064_ER_Total_Quant[Visvader_0064_ER_Total_Quant== FALSE])
Visvader_0064_ER_Total_Quant_Doublets_Percent <- Visvader_0064_ER_Total_Quant_Doublets /  (Visvader_0064_ER_Total_Quant_Doublets + Visvader_0064_ER_Total_Quant_Singlets) * 100
Visvader_0064_ER_Total_Quant <- as.data.frame(c(Visvader_0064_ER_Total_Quant_Singlets, Visvader_0064_ER_Total_Quant_Doublets, Visvader_0064_ER_Total_Quant_Doublets_Percent))
colnames(Visvader_0064_ER_Total_Quant) <- c("Visvader_0064_ER_Total_Quant")
rownames(Visvader_0064_ER_Total_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Visvader_0064_ER_Total_Quant, file = "/R/R_Visvader/Visvader_Doublet_Tables/Visvader_0064_ER_Total_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Visvader_0064_ER_Total_Singlets <- subset(Visvader_0064_ER_Total, cells=rownames(Visvader_0064_ER_Total@meta.data)[which(Visvader_0064_ER_Total@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Visvader_0064_ER_Total_Singlets, file = "/R/R_Visvader/Visvader_RDS/RDS_Total_Singlets/Visvader_0064_ER_Total_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Visvader_0064_ER_Total)
rm(Visvader_0064_ER_Total_Singlets)
rm(Visvader_0064_ER_Total_Quant)
rm(Visvader_0064_ER_Total_Quant_Singlets)
rm(Visvader_0064_ER_Total_Quant_Doublets)
rm(Visvader_0064_ER_Total_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Visvader_0064_ER_Total)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Visvader_0064_ER_Total)
rm(sweep.stats.Visvader_0064_ER_Total)
gc()




################################################################################################
########################              Visvader_0068_ER_Total             #######################
################################################################################################

# Load Data
Visvader_0068_ER_Total <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Total/Visvader_0068_ER_Total.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Visvader_0068_ER_Total <- paramSweep(Visvader_0068_ER_Total, PCs = 1:10, sct = FALSE)
sweep.stats.Visvader_0068_ER_Total <- summarizeSweep(sweep.res.Visvader_0068_ER_Total, GT = FALSE)
bcmvn_Visvader_0068_ER_Total <- find.pK(sweep.stats.Visvader_0068_ER_Total)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Visvader_0068_ER_Total[which.max(bcmvn_Visvader_0068_ER_Total$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Visvader_0068_ER_Total@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 5% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Visvader dataset = 5k per sample, will use doublet rate for 6,000 cells to err on the side of caution
nExp_poi <- round(0.05*nrow(Visvader_0068_ER_Total@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Visvader_0068_ER_Total <- doubletFinder(Visvader_0068_ER_Total, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Visvader_0068_ER_Total_Quant <- (Visvader_0068_ER_Total@meta.data$DF.classification == "Singlet")
Visvader_0068_ER_Total_Quant_Singlets <- length(Visvader_0068_ER_Total_Quant[Visvader_0068_ER_Total_Quant== TRUE])
Visvader_0068_ER_Total_Quant_Doublets <- length(Visvader_0068_ER_Total_Quant[Visvader_0068_ER_Total_Quant== FALSE])
Visvader_0068_ER_Total_Quant_Doublets_Percent <- Visvader_0068_ER_Total_Quant_Doublets /  (Visvader_0068_ER_Total_Quant_Doublets + Visvader_0068_ER_Total_Quant_Singlets) * 100
Visvader_0068_ER_Total_Quant <- as.data.frame(c(Visvader_0068_ER_Total_Quant_Singlets, Visvader_0068_ER_Total_Quant_Doublets, Visvader_0068_ER_Total_Quant_Doublets_Percent))
colnames(Visvader_0068_ER_Total_Quant) <- c("Visvader_0068_ER_Total_Quant")
rownames(Visvader_0068_ER_Total_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Visvader_0068_ER_Total_Quant, file = "/R/R_Visvader/Visvader_Doublet_Tables/Visvader_0068_ER_Total_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Visvader_0068_ER_Total_Singlets <- subset(Visvader_0068_ER_Total, cells=rownames(Visvader_0068_ER_Total@meta.data)[which(Visvader_0068_ER_Total@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Visvader_0068_ER_Total_Singlets, file = "/R/R_Visvader/Visvader_RDS/RDS_Total_Singlets/Visvader_0068_ER_Total_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Visvader_0068_ER_Total)
rm(Visvader_0068_ER_Total_Singlets)
rm(Visvader_0068_ER_Total_Quant)
rm(Visvader_0068_ER_Total_Quant_Singlets)
rm(Visvader_0068_ER_Total_Quant_Doublets)
rm(Visvader_0068_ER_Total_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Visvader_0068_ER_Total)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Visvader_0068_ER_Total)
rm(sweep.stats.Visvader_0068_ER_Total)
gc()




################################################################################################
########################              Visvader_0114_ER_Total             #######################
################################################################################################

# Load Data
Visvader_0114_ER_Total <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Total/Visvader_0114_ER_Total.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Visvader_0114_ER_Total <- paramSweep(Visvader_0114_ER_Total, PCs = 1:10, sct = FALSE)
sweep.stats.Visvader_0114_ER_Total <- summarizeSweep(sweep.res.Visvader_0114_ER_Total, GT = FALSE)
bcmvn_Visvader_0114_ER_Total <- find.pK(sweep.stats.Visvader_0114_ER_Total)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Visvader_0114_ER_Total[which.max(bcmvn_Visvader_0114_ER_Total$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Visvader_0114_ER_Total@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 5% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Visvader dataset = 5k per sample, will use doublet rate for 6,000 cells to err on the side of caution
nExp_poi <- round(0.05*nrow(Visvader_0114_ER_Total@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Visvader_0114_ER_Total <- doubletFinder(Visvader_0114_ER_Total, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Visvader_0114_ER_Total_Quant <- (Visvader_0114_ER_Total@meta.data$DF.classification == "Singlet")
Visvader_0114_ER_Total_Quant_Singlets <- length(Visvader_0114_ER_Total_Quant[Visvader_0114_ER_Total_Quant== TRUE])
Visvader_0114_ER_Total_Quant_Doublets <- length(Visvader_0114_ER_Total_Quant[Visvader_0114_ER_Total_Quant== FALSE])
Visvader_0114_ER_Total_Quant_Doublets_Percent <- Visvader_0114_ER_Total_Quant_Doublets /  (Visvader_0114_ER_Total_Quant_Doublets + Visvader_0114_ER_Total_Quant_Singlets) * 100
Visvader_0114_ER_Total_Quant <- as.data.frame(c(Visvader_0114_ER_Total_Quant_Singlets, Visvader_0114_ER_Total_Quant_Doublets, Visvader_0114_ER_Total_Quant_Doublets_Percent))
colnames(Visvader_0114_ER_Total_Quant) <- c("Visvader_0114_ER_Total_Quant")
rownames(Visvader_0114_ER_Total_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Visvader_0114_ER_Total_Quant, file = "/R/R_Visvader/Visvader_Doublet_Tables/Visvader_0114_ER_Total_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Visvader_0114_ER_Total_Singlets <- subset(Visvader_0114_ER_Total, cells=rownames(Visvader_0114_ER_Total@meta.data)[which(Visvader_0114_ER_Total@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Visvader_0114_ER_Total_Singlets, file = "/R/R_Visvader/Visvader_RDS/RDS_Total_Singlets/Visvader_0114_ER_Total_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Visvader_0114_ER_Total)
rm(Visvader_0114_ER_Total_Singlets)
rm(Visvader_0114_ER_Total_Quant)
rm(Visvader_0114_ER_Total_Quant_Singlets)
rm(Visvader_0114_ER_Total_Quant_Doublets)
rm(Visvader_0114_ER_Total_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Visvader_0114_ER_Total)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Visvader_0114_ER_Total)
rm(sweep.stats.Visvader_0114_ER_Total)
gc()




################################################################################################
########################              Visvader_0125_ER_Total             #######################
################################################################################################

# Load Data
Visvader_0125_ER_Total <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Total/Visvader_0125_ER_Total.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Visvader_0125_ER_Total <- paramSweep(Visvader_0125_ER_Total, PCs = 1:10, sct = FALSE)
sweep.stats.Visvader_0125_ER_Total <- summarizeSweep(sweep.res.Visvader_0125_ER_Total, GT = FALSE)
bcmvn_Visvader_0125_ER_Total <- find.pK(sweep.stats.Visvader_0125_ER_Total)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Visvader_0125_ER_Total[which.max(bcmvn_Visvader_0125_ER_Total$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Visvader_0125_ER_Total@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 5% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Visvader dataset = 5k per sample, will use doublet rate for 6,000 cells to err on the side of caution
nExp_poi <- round(0.05*nrow(Visvader_0125_ER_Total@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Visvader_0125_ER_Total <- doubletFinder(Visvader_0125_ER_Total, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Visvader_0125_ER_Total_Quant <- (Visvader_0125_ER_Total@meta.data$DF.classification == "Singlet")
Visvader_0125_ER_Total_Quant_Singlets <- length(Visvader_0125_ER_Total_Quant[Visvader_0125_ER_Total_Quant== TRUE])
Visvader_0125_ER_Total_Quant_Doublets <- length(Visvader_0125_ER_Total_Quant[Visvader_0125_ER_Total_Quant== FALSE])
Visvader_0125_ER_Total_Quant_Doublets_Percent <- Visvader_0125_ER_Total_Quant_Doublets /  (Visvader_0125_ER_Total_Quant_Doublets + Visvader_0125_ER_Total_Quant_Singlets) * 100
Visvader_0125_ER_Total_Quant <- as.data.frame(c(Visvader_0125_ER_Total_Quant_Singlets, Visvader_0125_ER_Total_Quant_Doublets, Visvader_0125_ER_Total_Quant_Doublets_Percent))
colnames(Visvader_0125_ER_Total_Quant) <- c("Visvader_0125_ER_Total_Quant")
rownames(Visvader_0125_ER_Total_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Visvader_0125_ER_Total_Quant, file = "/R/R_Visvader/Visvader_Doublet_Tables/Visvader_0125_ER_Total_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Visvader_0125_ER_Total_Singlets <- subset(Visvader_0125_ER_Total, cells=rownames(Visvader_0125_ER_Total@meta.data)[which(Visvader_0125_ER_Total@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Visvader_0125_ER_Total_Singlets, file = "/R/R_Visvader/Visvader_RDS/RDS_Total_Singlets/Visvader_0125_ER_Total_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Visvader_0125_ER_Total)
rm(Visvader_0125_ER_Total_Singlets)
rm(Visvader_0125_ER_Total_Quant)
rm(Visvader_0125_ER_Total_Quant_Singlets)
rm(Visvader_0125_ER_Total_Quant_Doublets)
rm(Visvader_0125_ER_Total_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Visvader_0125_ER_Total)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Visvader_0125_ER_Total)
rm(sweep.stats.Visvader_0125_ER_Total)
gc()




################################################################################################
########################              Visvader_0151_ER_Total             #######################
################################################################################################

# Load Data
Visvader_0151_ER_Total <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Total/Visvader_0151_ER_Total.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Visvader_0151_ER_Total <- paramSweep(Visvader_0151_ER_Total, PCs = 1:10, sct = FALSE)
sweep.stats.Visvader_0151_ER_Total <- summarizeSweep(sweep.res.Visvader_0151_ER_Total, GT = FALSE)
bcmvn_Visvader_0151_ER_Total <- find.pK(sweep.stats.Visvader_0151_ER_Total)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Visvader_0151_ER_Total[which.max(bcmvn_Visvader_0151_ER_Total$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Visvader_0151_ER_Total@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 5% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Visvader dataset = 5k per sample, will use doublet rate for 6,000 cells to err on the side of caution
nExp_poi <- round(0.05*nrow(Visvader_0151_ER_Total@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Visvader_0151_ER_Total <- doubletFinder(Visvader_0151_ER_Total, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Visvader_0151_ER_Total_Quant <- (Visvader_0151_ER_Total@meta.data$DF.classification == "Singlet")
Visvader_0151_ER_Total_Quant_Singlets <- length(Visvader_0151_ER_Total_Quant[Visvader_0151_ER_Total_Quant== TRUE])
Visvader_0151_ER_Total_Quant_Doublets <- length(Visvader_0151_ER_Total_Quant[Visvader_0151_ER_Total_Quant== FALSE])
Visvader_0151_ER_Total_Quant_Doublets_Percent <- Visvader_0151_ER_Total_Quant_Doublets /  (Visvader_0151_ER_Total_Quant_Doublets + Visvader_0151_ER_Total_Quant_Singlets) * 100
Visvader_0151_ER_Total_Quant <- as.data.frame(c(Visvader_0151_ER_Total_Quant_Singlets, Visvader_0151_ER_Total_Quant_Doublets, Visvader_0151_ER_Total_Quant_Doublets_Percent))
colnames(Visvader_0151_ER_Total_Quant) <- c("Visvader_0151_ER_Total_Quant")
rownames(Visvader_0151_ER_Total_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Visvader_0151_ER_Total_Quant, file = "/R/R_Visvader/Visvader_Doublet_Tables/Visvader_0151_ER_Total_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Visvader_0151_ER_Total_Singlets <- subset(Visvader_0151_ER_Total, cells=rownames(Visvader_0151_ER_Total@meta.data)[which(Visvader_0151_ER_Total@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Visvader_0151_ER_Total_Singlets, file = "/R/R_Visvader/Visvader_RDS/RDS_Total_Singlets/Visvader_0151_ER_Total_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Visvader_0151_ER_Total)
rm(Visvader_0151_ER_Total_Singlets)
rm(Visvader_0151_ER_Total_Quant)
rm(Visvader_0151_ER_Total_Quant_Singlets)
rm(Visvader_0151_ER_Total_Quant_Doublets)
rm(Visvader_0151_ER_Total_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Visvader_0151_ER_Total)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Visvader_0151_ER_Total)
rm(sweep.stats.Visvader_0151_ER_Total)
gc()




################################################################################################
########################              Visvader_0163_ER_Total             #######################
################################################################################################

# Load Data
Visvader_0163_ER_Total <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Total/Visvader_0163_ER_Total.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Visvader_0163_ER_Total <- paramSweep(Visvader_0163_ER_Total, PCs = 1:10, sct = FALSE)
sweep.stats.Visvader_0163_ER_Total <- summarizeSweep(sweep.res.Visvader_0163_ER_Total, GT = FALSE)
bcmvn_Visvader_0163_ER_Total <- find.pK(sweep.stats.Visvader_0163_ER_Total)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Visvader_0163_ER_Total[which.max(bcmvn_Visvader_0163_ER_Total$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Visvader_0163_ER_Total@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 5% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Visvader dataset = 5k per sample, will use doublet rate for 6,000 cells to err on the side of caution
nExp_poi <- round(0.05*nrow(Visvader_0163_ER_Total@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Visvader_0163_ER_Total <- doubletFinder(Visvader_0163_ER_Total, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Visvader_0163_ER_Total_Quant <- (Visvader_0163_ER_Total@meta.data$DF.classification == "Singlet")
Visvader_0163_ER_Total_Quant_Singlets <- length(Visvader_0163_ER_Total_Quant[Visvader_0163_ER_Total_Quant== TRUE])
Visvader_0163_ER_Total_Quant_Doublets <- length(Visvader_0163_ER_Total_Quant[Visvader_0163_ER_Total_Quant== FALSE])
Visvader_0163_ER_Total_Quant_Doublets_Percent <- Visvader_0163_ER_Total_Quant_Doublets /  (Visvader_0163_ER_Total_Quant_Doublets + Visvader_0163_ER_Total_Quant_Singlets) * 100
Visvader_0163_ER_Total_Quant <- as.data.frame(c(Visvader_0163_ER_Total_Quant_Singlets, Visvader_0163_ER_Total_Quant_Doublets, Visvader_0163_ER_Total_Quant_Doublets_Percent))
colnames(Visvader_0163_ER_Total_Quant) <- c("Visvader_0163_ER_Total_Quant")
rownames(Visvader_0163_ER_Total_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Visvader_0163_ER_Total_Quant, file = "/R/R_Visvader/Visvader_Doublet_Tables/Visvader_0163_ER_Total_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Visvader_0163_ER_Total_Singlets <- subset(Visvader_0163_ER_Total, cells=rownames(Visvader_0163_ER_Total@meta.data)[which(Visvader_0163_ER_Total@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Visvader_0163_ER_Total_Singlets, file = "/R/R_Visvader/Visvader_RDS/RDS_Total_Singlets/Visvader_0163_ER_Total_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Visvader_0163_ER_Total)
rm(Visvader_0163_ER_Total_Singlets)
rm(Visvader_0163_ER_Total_Quant)
rm(Visvader_0163_ER_Total_Quant_Singlets)
rm(Visvader_0163_ER_Total_Quant_Doublets)
rm(Visvader_0163_ER_Total_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Visvader_0163_ER_Total)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Visvader_0163_ER_Total)
rm(sweep.stats.Visvader_0163_ER_Total)
gc()




################################################################################################
########################              Visvader_0167_ER_Total             #######################
################################################################################################

# Load Data
Visvader_0167_ER_Total <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Total/Visvader_0167_ER_Total.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Visvader_0167_ER_Total <- paramSweep(Visvader_0167_ER_Total, PCs = 1:10, sct = FALSE)
sweep.stats.Visvader_0167_ER_Total <- summarizeSweep(sweep.res.Visvader_0167_ER_Total, GT = FALSE)
bcmvn_Visvader_0167_ER_Total <- find.pK(sweep.stats.Visvader_0167_ER_Total)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Visvader_0167_ER_Total[which.max(bcmvn_Visvader_0167_ER_Total$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Visvader_0167_ER_Total@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 5% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Visvader dataset = 5k per sample, will use doublet rate for 6,000 cells to err on the side of caution
nExp_poi <- round(0.05*nrow(Visvader_0167_ER_Total@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Visvader_0167_ER_Total <- doubletFinder(Visvader_0167_ER_Total, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Visvader_0167_ER_Total_Quant <- (Visvader_0167_ER_Total@meta.data$DF.classification == "Singlet")
Visvader_0167_ER_Total_Quant_Singlets <- length(Visvader_0167_ER_Total_Quant[Visvader_0167_ER_Total_Quant== TRUE])
Visvader_0167_ER_Total_Quant_Doublets <- length(Visvader_0167_ER_Total_Quant[Visvader_0167_ER_Total_Quant== FALSE])
Visvader_0167_ER_Total_Quant_Doublets_Percent <- Visvader_0167_ER_Total_Quant_Doublets /  (Visvader_0167_ER_Total_Quant_Doublets + Visvader_0167_ER_Total_Quant_Singlets) * 100
Visvader_0167_ER_Total_Quant <- as.data.frame(c(Visvader_0167_ER_Total_Quant_Singlets, Visvader_0167_ER_Total_Quant_Doublets, Visvader_0167_ER_Total_Quant_Doublets_Percent))
colnames(Visvader_0167_ER_Total_Quant) <- c("Visvader_0167_ER_Total_Quant")
rownames(Visvader_0167_ER_Total_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Visvader_0167_ER_Total_Quant, file = "/R/R_Visvader/Visvader_Doublet_Tables/Visvader_0167_ER_Total_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Visvader_0167_ER_Total_Singlets <- subset(Visvader_0167_ER_Total, cells=rownames(Visvader_0167_ER_Total@meta.data)[which(Visvader_0167_ER_Total@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Visvader_0167_ER_Total_Singlets, file = "/R/R_Visvader/Visvader_RDS/RDS_Total_Singlets/Visvader_0167_ER_Total_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Visvader_0167_ER_Total)
rm(Visvader_0167_ER_Total_Singlets)
rm(Visvader_0167_ER_Total_Quant)
rm(Visvader_0167_ER_Total_Quant_Singlets)
rm(Visvader_0167_ER_Total_Quant_Doublets)
rm(Visvader_0167_ER_Total_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Visvader_0167_ER_Total)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Visvader_0167_ER_Total)
rm(sweep.stats.Visvader_0167_ER_Total)
gc()




################################################################################################
########################              Visvader_0173_ER_Total             #######################
################################################################################################

# Load Data
Visvader_0173_ER_Total <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Total/Visvader_0173_ER_Total.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Visvader_0173_ER_Total <- paramSweep(Visvader_0173_ER_Total, PCs = 1:10, sct = FALSE)
sweep.stats.Visvader_0173_ER_Total <- summarizeSweep(sweep.res.Visvader_0173_ER_Total, GT = FALSE)
bcmvn_Visvader_0173_ER_Total <- find.pK(sweep.stats.Visvader_0173_ER_Total)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Visvader_0173_ER_Total[which.max(bcmvn_Visvader_0173_ER_Total$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Visvader_0173_ER_Total@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 5% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Visvader dataset = 5k per sample, will use doublet rate for 6,000 cells to err on the side of caution
nExp_poi <- round(0.05*nrow(Visvader_0173_ER_Total@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Visvader_0173_ER_Total <- doubletFinder(Visvader_0173_ER_Total, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Visvader_0173_ER_Total_Quant <- (Visvader_0173_ER_Total@meta.data$DF.classification == "Singlet")
Visvader_0173_ER_Total_Quant_Singlets <- length(Visvader_0173_ER_Total_Quant[Visvader_0173_ER_Total_Quant== TRUE])
Visvader_0173_ER_Total_Quant_Doublets <- length(Visvader_0173_ER_Total_Quant[Visvader_0173_ER_Total_Quant== FALSE])
Visvader_0173_ER_Total_Quant_Doublets_Percent <- Visvader_0173_ER_Total_Quant_Doublets /  (Visvader_0173_ER_Total_Quant_Doublets + Visvader_0173_ER_Total_Quant_Singlets) * 100
Visvader_0173_ER_Total_Quant <- as.data.frame(c(Visvader_0173_ER_Total_Quant_Singlets, Visvader_0173_ER_Total_Quant_Doublets, Visvader_0173_ER_Total_Quant_Doublets_Percent))
colnames(Visvader_0173_ER_Total_Quant) <- c("Visvader_0173_ER_Total_Quant")
rownames(Visvader_0173_ER_Total_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Visvader_0173_ER_Total_Quant, file = "/R/R_Visvader/Visvader_Doublet_Tables/Visvader_0173_ER_Total_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Visvader_0173_ER_Total_Singlets <- subset(Visvader_0173_ER_Total, cells=rownames(Visvader_0173_ER_Total@meta.data)[which(Visvader_0173_ER_Total@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Visvader_0173_ER_Total_Singlets, file = "/R/R_Visvader/Visvader_RDS/RDS_Total_Singlets/Visvader_0173_ER_Total_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Visvader_0173_ER_Total)
rm(Visvader_0173_ER_Total_Singlets)
rm(Visvader_0173_ER_Total_Quant)
rm(Visvader_0173_ER_Total_Quant_Singlets)
rm(Visvader_0173_ER_Total_Quant_Doublets)
rm(Visvader_0173_ER_Total_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Visvader_0173_ER_Total)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Visvader_0173_ER_Total)
rm(sweep.stats.Visvader_0173_ER_Total)
gc()




################################################################################################
########################              Visvader_0178_ER_Total             #######################
################################################################################################

# Load Data
Visvader_0178_ER_Total <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Total/Visvader_0178_ER_Total.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Visvader_0178_ER_Total <- paramSweep(Visvader_0178_ER_Total, PCs = 1:10, sct = FALSE)
sweep.stats.Visvader_0178_ER_Total <- summarizeSweep(sweep.res.Visvader_0178_ER_Total, GT = FALSE)
bcmvn_Visvader_0178_ER_Total <- find.pK(sweep.stats.Visvader_0178_ER_Total)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Visvader_0178_ER_Total[which.max(bcmvn_Visvader_0178_ER_Total$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Visvader_0178_ER_Total@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 5% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Visvader dataset = 5k per sample, will use doublet rate for 6,000 cells to err on the side of caution
nExp_poi <- round(0.05*nrow(Visvader_0178_ER_Total@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Visvader_0178_ER_Total <- doubletFinder(Visvader_0178_ER_Total, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Visvader_0178_ER_Total_Quant <- (Visvader_0178_ER_Total@meta.data$DF.classification == "Singlet")
Visvader_0178_ER_Total_Quant_Singlets <- length(Visvader_0178_ER_Total_Quant[Visvader_0178_ER_Total_Quant== TRUE])
Visvader_0178_ER_Total_Quant_Doublets <- length(Visvader_0178_ER_Total_Quant[Visvader_0178_ER_Total_Quant== FALSE])
Visvader_0178_ER_Total_Quant_Doublets_Percent <- Visvader_0178_ER_Total_Quant_Doublets /  (Visvader_0178_ER_Total_Quant_Doublets + Visvader_0178_ER_Total_Quant_Singlets) * 100
Visvader_0178_ER_Total_Quant <- as.data.frame(c(Visvader_0178_ER_Total_Quant_Singlets, Visvader_0178_ER_Total_Quant_Doublets, Visvader_0178_ER_Total_Quant_Doublets_Percent))
colnames(Visvader_0178_ER_Total_Quant) <- c("Visvader_0178_ER_Total_Quant")
rownames(Visvader_0178_ER_Total_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Visvader_0178_ER_Total_Quant, file = "/R/R_Visvader/Visvader_Doublet_Tables/Visvader_0178_ER_Total_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Visvader_0178_ER_Total_Singlets <- subset(Visvader_0178_ER_Total, cells=rownames(Visvader_0178_ER_Total@meta.data)[which(Visvader_0178_ER_Total@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Visvader_0178_ER_Total_Singlets, file = "/R/R_Visvader/Visvader_RDS/RDS_Total_Singlets/Visvader_0178_ER_Total_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Visvader_0178_ER_Total)
rm(Visvader_0178_ER_Total_Singlets)
rm(Visvader_0178_ER_Total_Quant)
rm(Visvader_0178_ER_Total_Quant_Singlets)
rm(Visvader_0178_ER_Total_Quant_Doublets)
rm(Visvader_0178_ER_Total_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Visvader_0178_ER_Total)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Visvader_0178_ER_Total)
rm(sweep.stats.Visvader_0178_ER_Total)
gc()




################################################################################################
########################              Visvader_0360_ER_Total             #######################
################################################################################################

# Load Data
Visvader_0360_ER_Total <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Total/Visvader_0360_ER_Total.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Visvader_0360_ER_Total <- paramSweep(Visvader_0360_ER_Total, PCs = 1:10, sct = FALSE)
sweep.stats.Visvader_0360_ER_Total <- summarizeSweep(sweep.res.Visvader_0360_ER_Total, GT = FALSE)
bcmvn_Visvader_0360_ER_Total <- find.pK(sweep.stats.Visvader_0360_ER_Total)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Visvader_0360_ER_Total[which.max(bcmvn_Visvader_0360_ER_Total$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Visvader_0360_ER_Total@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 5% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Visvader dataset = 5k per sample, will use doublet rate for 6,000 cells to err on the side of caution
nExp_poi <- round(0.05*nrow(Visvader_0360_ER_Total@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Visvader_0360_ER_Total <- doubletFinder(Visvader_0360_ER_Total, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Visvader_0360_ER_Total_Quant <- (Visvader_0360_ER_Total@meta.data$DF.classification == "Singlet")
Visvader_0360_ER_Total_Quant_Singlets <- length(Visvader_0360_ER_Total_Quant[Visvader_0360_ER_Total_Quant== TRUE])
Visvader_0360_ER_Total_Quant_Doublets <- length(Visvader_0360_ER_Total_Quant[Visvader_0360_ER_Total_Quant== FALSE])
Visvader_0360_ER_Total_Quant_Doublets_Percent <- Visvader_0360_ER_Total_Quant_Doublets /  (Visvader_0360_ER_Total_Quant_Doublets + Visvader_0360_ER_Total_Quant_Singlets) * 100
Visvader_0360_ER_Total_Quant <- as.data.frame(c(Visvader_0360_ER_Total_Quant_Singlets, Visvader_0360_ER_Total_Quant_Doublets, Visvader_0360_ER_Total_Quant_Doublets_Percent))
colnames(Visvader_0360_ER_Total_Quant) <- c("Visvader_0360_ER_Total_Quant")
rownames(Visvader_0360_ER_Total_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Visvader_0360_ER_Total_Quant, file = "/R/R_Visvader/Visvader_Doublet_Tables/Visvader_0360_ER_Total_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Visvader_0360_ER_Total_Singlets <- subset(Visvader_0360_ER_Total, cells=rownames(Visvader_0360_ER_Total@meta.data)[which(Visvader_0360_ER_Total@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Visvader_0360_ER_Total_Singlets, file = "/R/R_Visvader/Visvader_RDS/RDS_Total_Singlets/Visvader_0360_ER_Total_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Visvader_0360_ER_Total)
rm(Visvader_0360_ER_Total_Singlets)
rm(Visvader_0360_ER_Total_Quant)
rm(Visvader_0360_ER_Total_Quant_Singlets)
rm(Visvader_0360_ER_Total_Quant_Doublets)
rm(Visvader_0360_ER_Total_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Visvader_0360_ER_Total)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Visvader_0360_ER_Total)
rm(sweep.stats.Visvader_0360_ER_Total)
gc()



################################################################################################
########################              Visvader_0319_PR_Total             #######################
################################################################################################

# Load Data
Visvader_0319_PR_Total <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Total/Visvader_0319_PR_Total.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Visvader_0319_PR_Total <- paramSweep(Visvader_0319_PR_Total, PCs = 1:10, sct = FALSE)
sweep.stats.Visvader_0319_PR_Total <- summarizeSweep(sweep.res.Visvader_0319_PR_Total, GT = FALSE)
bcmvn_Visvader_0319_PR_Total <- find.pK(sweep.stats.Visvader_0319_PR_Total)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Visvader_0319_PR_Total[which.max(bcmvn_Visvader_0319_PR_Total$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Visvader_0319_PR_Total@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 5% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Visvader dataset = 5k per sample, will use doublet rate for 6,000 cells to err on the side of caution
nExp_poi <- round(0.05*nrow(Visvader_0319_PR_Total@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Visvader_0319_PR_Total <- doubletFinder(Visvader_0319_PR_Total, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Visvader_0319_PR_Total_Quant <- (Visvader_0319_PR_Total@meta.data$DF.classification == "Singlet")
Visvader_0319_PR_Total_Quant_Singlets <- length(Visvader_0319_PR_Total_Quant[Visvader_0319_PR_Total_Quant== TRUE])
Visvader_0319_PR_Total_Quant_Doublets <- length(Visvader_0319_PR_Total_Quant[Visvader_0319_PR_Total_Quant== FALSE])
Visvader_0319_PR_Total_Quant_Doublets_Percent <- Visvader_0319_PR_Total_Quant_Doublets /  (Visvader_0319_PR_Total_Quant_Doublets + Visvader_0319_PR_Total_Quant_Singlets) * 100
Visvader_0319_PR_Total_Quant <- as.data.frame(c(Visvader_0319_PR_Total_Quant_Singlets, Visvader_0319_PR_Total_Quant_Doublets, Visvader_0319_PR_Total_Quant_Doublets_Percent))
colnames(Visvader_0319_PR_Total_Quant) <- c("Visvader_0319_PR_Total_Quant")
rownames(Visvader_0319_PR_Total_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Visvader_0319_PR_Total_Quant, file = "/R/R_Visvader/Visvader_Doublet_Tables/Visvader_0319_PR_Total_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Visvader_0319_PR_Total_Singlets <- subset(Visvader_0319_PR_Total, cells=rownames(Visvader_0319_PR_Total@meta.data)[which(Visvader_0319_PR_Total@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Visvader_0319_PR_Total_Singlets, file = "/R/R_Visvader/Visvader_RDS/RDS_Total_Singlets/Visvader_0319_PR_Total_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Visvader_0319_PR_Total)
rm(Visvader_0319_PR_Total_Singlets)
rm(Visvader_0319_PR_Total_Quant)
rm(Visvader_0319_PR_Total_Quant_Singlets)
rm(Visvader_0319_PR_Total_Quant_Doublets)
rm(Visvader_0319_PR_Total_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Visvader_0319_PR_Total)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Visvader_0319_PR_Total)
rm(sweep.stats.Visvader_0319_PR_Total)
gc()


################################################################################################
########################                HER2-pos Tumors                #########################
################################################################################################

################################################################################################
########################              Visvader_0031_HER2_Total             #####################
################################################################################################

# Load Data
Visvader_0031_HER2_Total <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Total/Visvader_0031_HER2_Total.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Visvader_0031_HER2_Total <- paramSweep(Visvader_0031_HER2_Total, PCs = 1:10, sct = FALSE)
sweep.stats.Visvader_0031_HER2_Total <- summarizeSweep(sweep.res.Visvader_0031_HER2_Total, GT = FALSE)
bcmvn_Visvader_0031_HER2_Total <- find.pK(sweep.stats.Visvader_0031_HER2_Total)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Visvader_0031_HER2_Total[which.max(bcmvn_Visvader_0031_HER2_Total$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Visvader_0031_HER2_Total@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 5% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Visvader dataset = 5k per sample, will use doublet rate for 6,000 cells to err on the side of caution
nExp_poi <- round(0.05*nrow(Visvader_0031_HER2_Total@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Visvader_0031_HER2_Total <- doubletFinder(Visvader_0031_HER2_Total, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Visvader_0031_HER2_Total_Quant <- (Visvader_0031_HER2_Total@meta.data$DF.classification == "Singlet")
Visvader_0031_HER2_Total_Quant_Singlets <- length(Visvader_0031_HER2_Total_Quant[Visvader_0031_HER2_Total_Quant== TRUE])
Visvader_0031_HER2_Total_Quant_Doublets <- length(Visvader_0031_HER2_Total_Quant[Visvader_0031_HER2_Total_Quant== FALSE])
Visvader_0031_HER2_Total_Quant_Doublets_Percent <- Visvader_0031_HER2_Total_Quant_Doublets /  (Visvader_0031_HER2_Total_Quant_Doublets + Visvader_0031_HER2_Total_Quant_Singlets) * 100
Visvader_0031_HER2_Total_Quant <- as.data.frame(c(Visvader_0031_HER2_Total_Quant_Singlets, Visvader_0031_HER2_Total_Quant_Doublets, Visvader_0031_HER2_Total_Quant_Doublets_Percent))
colnames(Visvader_0031_HER2_Total_Quant) <- c("Visvader_0031_HER2_Total_Quant")
rownames(Visvader_0031_HER2_Total_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Visvader_0031_HER2_Total_Quant, file = "/R/R_Visvader/Visvader_Doublet_Tables/Visvader_0031_HER2_Total_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Visvader_0031_HER2_Total_Singlets <- subset(Visvader_0031_HER2_Total, cells=rownames(Visvader_0031_HER2_Total@meta.data)[which(Visvader_0031_HER2_Total@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Visvader_0031_HER2_Total_Singlets, file = "/R/R_Visvader/Visvader_RDS/RDS_Total_Singlets/Visvader_0031_HER2_Total_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Visvader_0031_HER2_Total)
rm(Visvader_0031_HER2_Total_Singlets)
rm(Visvader_0031_HER2_Total_Quant)
rm(Visvader_0031_HER2_Total_Quant_Singlets)
rm(Visvader_0031_HER2_Total_Quant_Doublets)
rm(Visvader_0031_HER2_Total_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Visvader_0031_HER2_Total)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Visvader_0031_HER2_Total)
rm(sweep.stats.Visvader_0031_HER2_Total)
gc()



################################################################################################
########################              Visvader_0069_HER2_Total             #####################
################################################################################################

# Load Data
Visvader_0069_HER2_Total <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Total/Visvader_0069_HER2_Total.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Visvader_0069_HER2_Total <- paramSweep(Visvader_0069_HER2_Total, PCs = 1:10, sct = FALSE)
sweep.stats.Visvader_0069_HER2_Total <- summarizeSweep(sweep.res.Visvader_0069_HER2_Total, GT = FALSE)
bcmvn_Visvader_0069_HER2_Total <- find.pK(sweep.stats.Visvader_0069_HER2_Total)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Visvader_0069_HER2_Total[which.max(bcmvn_Visvader_0069_HER2_Total$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Visvader_0069_HER2_Total@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 5% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Visvader dataset = 5k per sample, will use doublet rate for 6,000 cells to err on the side of caution
nExp_poi <- round(0.05*nrow(Visvader_0069_HER2_Total@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Visvader_0069_HER2_Total <- doubletFinder(Visvader_0069_HER2_Total, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Visvader_0069_HER2_Total_Quant <- (Visvader_0069_HER2_Total@meta.data$DF.classification == "Singlet")
Visvader_0069_HER2_Total_Quant_Singlets <- length(Visvader_0069_HER2_Total_Quant[Visvader_0069_HER2_Total_Quant== TRUE])
Visvader_0069_HER2_Total_Quant_Doublets <- length(Visvader_0069_HER2_Total_Quant[Visvader_0069_HER2_Total_Quant== FALSE])
Visvader_0069_HER2_Total_Quant_Doublets_Percent <- Visvader_0069_HER2_Total_Quant_Doublets /  (Visvader_0069_HER2_Total_Quant_Doublets + Visvader_0069_HER2_Total_Quant_Singlets) * 100
Visvader_0069_HER2_Total_Quant <- as.data.frame(c(Visvader_0069_HER2_Total_Quant_Singlets, Visvader_0069_HER2_Total_Quant_Doublets, Visvader_0069_HER2_Total_Quant_Doublets_Percent))
colnames(Visvader_0069_HER2_Total_Quant) <- c("Visvader_0069_HER2_Total_Quant")
rownames(Visvader_0069_HER2_Total_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Visvader_0069_HER2_Total_Quant, file = "/R/R_Visvader/Visvader_Doublet_Tables/Visvader_0069_HER2_Total_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Visvader_0069_HER2_Total_Singlets <- subset(Visvader_0069_HER2_Total, cells=rownames(Visvader_0069_HER2_Total@meta.data)[which(Visvader_0069_HER2_Total@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Visvader_0069_HER2_Total_Singlets, file = "/R/R_Visvader/Visvader_RDS/RDS_Total_Singlets/Visvader_0069_HER2_Total_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Visvader_0069_HER2_Total)
rm(Visvader_0069_HER2_Total_Singlets)
rm(Visvader_0069_HER2_Total_Quant)
rm(Visvader_0069_HER2_Total_Quant_Singlets)
rm(Visvader_0069_HER2_Total_Quant_Doublets)
rm(Visvader_0069_HER2_Total_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Visvader_0069_HER2_Total)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Visvader_0069_HER2_Total)
rm(sweep.stats.Visvader_0069_HER2_Total)
gc()



################################################################################################
########################              Visvader_0161_HER2_Total             #####################
################################################################################################

# Load Data
Visvader_0161_HER2_Total <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Total/Visvader_0161_HER2_Total.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Visvader_0161_HER2_Total <- paramSweep(Visvader_0161_HER2_Total, PCs = 1:10, sct = FALSE)
sweep.stats.Visvader_0161_HER2_Total <- summarizeSweep(sweep.res.Visvader_0161_HER2_Total, GT = FALSE)
bcmvn_Visvader_0161_HER2_Total <- find.pK(sweep.stats.Visvader_0161_HER2_Total)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Visvader_0161_HER2_Total[which.max(bcmvn_Visvader_0161_HER2_Total$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Visvader_0161_HER2_Total@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 5% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Visvader dataset = 5k per sample, will use doublet rate for 6,000 cells to err on the side of caution
nExp_poi <- round(0.05*nrow(Visvader_0161_HER2_Total@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Visvader_0161_HER2_Total <- doubletFinder(Visvader_0161_HER2_Total, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Visvader_0161_HER2_Total_Quant <- (Visvader_0161_HER2_Total@meta.data$DF.classification == "Singlet")
Visvader_0161_HER2_Total_Quant_Singlets <- length(Visvader_0161_HER2_Total_Quant[Visvader_0161_HER2_Total_Quant== TRUE])
Visvader_0161_HER2_Total_Quant_Doublets <- length(Visvader_0161_HER2_Total_Quant[Visvader_0161_HER2_Total_Quant== FALSE])
Visvader_0161_HER2_Total_Quant_Doublets_Percent <- Visvader_0161_HER2_Total_Quant_Doublets /  (Visvader_0161_HER2_Total_Quant_Doublets + Visvader_0161_HER2_Total_Quant_Singlets) * 100
Visvader_0161_HER2_Total_Quant <- as.data.frame(c(Visvader_0161_HER2_Total_Quant_Singlets, Visvader_0161_HER2_Total_Quant_Doublets, Visvader_0161_HER2_Total_Quant_Doublets_Percent))
colnames(Visvader_0161_HER2_Total_Quant) <- c("Visvader_0161_HER2_Total_Quant")
rownames(Visvader_0161_HER2_Total_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Visvader_0161_HER2_Total_Quant, file = "/R/R_Visvader/Visvader_Doublet_Tables/Visvader_0161_HER2_Total_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Visvader_0161_HER2_Total_Singlets <- subset(Visvader_0161_HER2_Total, cells=rownames(Visvader_0161_HER2_Total@meta.data)[which(Visvader_0161_HER2_Total@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Visvader_0161_HER2_Total_Singlets, file = "/R/R_Visvader/Visvader_RDS/RDS_Total_Singlets/Visvader_0161_HER2_Total_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Visvader_0161_HER2_Total)
rm(Visvader_0161_HER2_Total_Singlets)
rm(Visvader_0161_HER2_Total_Quant)
rm(Visvader_0161_HER2_Total_Quant_Singlets)
rm(Visvader_0161_HER2_Total_Quant_Doublets)
rm(Visvader_0161_HER2_Total_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Visvader_0161_HER2_Total)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Visvader_0161_HER2_Total)
rm(sweep.stats.Visvader_0161_HER2_Total)
gc()



################################################################################################
########################              Visvader_0176_HER2_Total             #####################
################################################################################################

# Load Data
Visvader_0176_HER2_Total <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Total/Visvader_0176_HER2_Total.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Visvader_0176_HER2_Total <- paramSweep(Visvader_0176_HER2_Total, PCs = 1:10, sct = FALSE)
sweep.stats.Visvader_0176_HER2_Total <- summarizeSweep(sweep.res.Visvader_0176_HER2_Total, GT = FALSE)
bcmvn_Visvader_0176_HER2_Total <- find.pK(sweep.stats.Visvader_0176_HER2_Total)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Visvader_0176_HER2_Total[which.max(bcmvn_Visvader_0176_HER2_Total$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Visvader_0176_HER2_Total@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 5% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Visvader dataset = 5k per sample, will use doublet rate for 6,000 cells to err on the side of caution
nExp_poi <- round(0.05*nrow(Visvader_0176_HER2_Total@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Visvader_0176_HER2_Total <- doubletFinder(Visvader_0176_HER2_Total, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Visvader_0176_HER2_Total_Quant <- (Visvader_0176_HER2_Total@meta.data$DF.classification == "Singlet")
Visvader_0176_HER2_Total_Quant_Singlets <- length(Visvader_0176_HER2_Total_Quant[Visvader_0176_HER2_Total_Quant== TRUE])
Visvader_0176_HER2_Total_Quant_Doublets <- length(Visvader_0176_HER2_Total_Quant[Visvader_0176_HER2_Total_Quant== FALSE])
Visvader_0176_HER2_Total_Quant_Doublets_Percent <- Visvader_0176_HER2_Total_Quant_Doublets /  (Visvader_0176_HER2_Total_Quant_Doublets + Visvader_0176_HER2_Total_Quant_Singlets) * 100
Visvader_0176_HER2_Total_Quant <- as.data.frame(c(Visvader_0176_HER2_Total_Quant_Singlets, Visvader_0176_HER2_Total_Quant_Doublets, Visvader_0176_HER2_Total_Quant_Doublets_Percent))
colnames(Visvader_0176_HER2_Total_Quant) <- c("Visvader_0176_HER2_Total_Quant")
rownames(Visvader_0176_HER2_Total_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Visvader_0176_HER2_Total_Quant, file = "/R/R_Visvader/Visvader_Doublet_Tables/Visvader_0176_HER2_Total_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Visvader_0176_HER2_Total_Singlets <- subset(Visvader_0176_HER2_Total, cells=rownames(Visvader_0176_HER2_Total@meta.data)[which(Visvader_0176_HER2_Total@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Visvader_0176_HER2_Total_Singlets, file = "/R/R_Visvader/Visvader_RDS/RDS_Total_Singlets/Visvader_0176_HER2_Total_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Visvader_0176_HER2_Total)
rm(Visvader_0176_HER2_Total_Singlets)
rm(Visvader_0176_HER2_Total_Quant)
rm(Visvader_0176_HER2_Total_Quant_Singlets)
rm(Visvader_0176_HER2_Total_Quant_Doublets)
rm(Visvader_0176_HER2_Total_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Visvader_0176_HER2_Total)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Visvader_0176_HER2_Total)
rm(sweep.stats.Visvader_0176_HER2_Total)
gc()



################################################################################################
########################              Visvader_0308_HER2_Total             #####################
################################################################################################

# Load Data
Visvader_0308_HER2_Total <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Total/Visvader_0308_HER2_Total.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Visvader_0308_HER2_Total <- paramSweep(Visvader_0308_HER2_Total, PCs = 1:10, sct = FALSE)
sweep.stats.Visvader_0308_HER2_Total <- summarizeSweep(sweep.res.Visvader_0308_HER2_Total, GT = FALSE)
bcmvn_Visvader_0308_HER2_Total <- find.pK(sweep.stats.Visvader_0308_HER2_Total)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Visvader_0308_HER2_Total[which.max(bcmvn_Visvader_0308_HER2_Total$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Visvader_0308_HER2_Total@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 5% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Visvader dataset = 5k per sample, will use doublet rate for 6,000 cells to err on the side of caution
nExp_poi <- round(0.05*nrow(Visvader_0308_HER2_Total@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Visvader_0308_HER2_Total <- doubletFinder(Visvader_0308_HER2_Total, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Visvader_0308_HER2_Total_Quant <- (Visvader_0308_HER2_Total@meta.data$DF.classification == "Singlet")
Visvader_0308_HER2_Total_Quant_Singlets <- length(Visvader_0308_HER2_Total_Quant[Visvader_0308_HER2_Total_Quant== TRUE])
Visvader_0308_HER2_Total_Quant_Doublets <- length(Visvader_0308_HER2_Total_Quant[Visvader_0308_HER2_Total_Quant== FALSE])
Visvader_0308_HER2_Total_Quant_Doublets_Percent <- Visvader_0308_HER2_Total_Quant_Doublets /  (Visvader_0308_HER2_Total_Quant_Doublets + Visvader_0308_HER2_Total_Quant_Singlets) * 100
Visvader_0308_HER2_Total_Quant <- as.data.frame(c(Visvader_0308_HER2_Total_Quant_Singlets, Visvader_0308_HER2_Total_Quant_Doublets, Visvader_0308_HER2_Total_Quant_Doublets_Percent))
colnames(Visvader_0308_HER2_Total_Quant) <- c("Visvader_0308_HER2_Total_Quant")
rownames(Visvader_0308_HER2_Total_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Visvader_0308_HER2_Total_Quant, file = "/R/R_Visvader/Visvader_Doublet_Tables/Visvader_0308_HER2_Total_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Visvader_0308_HER2_Total_Singlets <- subset(Visvader_0308_HER2_Total, cells=rownames(Visvader_0308_HER2_Total@meta.data)[which(Visvader_0308_HER2_Total@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Visvader_0308_HER2_Total_Singlets, file = "/R/R_Visvader/Visvader_RDS/RDS_Total_Singlets/Visvader_0308_HER2_Total_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Visvader_0308_HER2_Total)
rm(Visvader_0308_HER2_Total_Singlets)
rm(Visvader_0308_HER2_Total_Quant)
rm(Visvader_0308_HER2_Total_Quant_Singlets)
rm(Visvader_0308_HER2_Total_Quant_Doublets)
rm(Visvader_0308_HER2_Total_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Visvader_0308_HER2_Total)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Visvader_0308_HER2_Total)
rm(sweep.stats.Visvader_0308_HER2_Total)
gc()



################################################################################################
########################              Visvader_0337_HER2_Total             #####################
################################################################################################

# Load Data
Visvader_0337_HER2_Total <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Total/Visvader_0337_HER2_Total.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Visvader_0337_HER2_Total <- paramSweep(Visvader_0337_HER2_Total, PCs = 1:10, sct = FALSE)
sweep.stats.Visvader_0337_HER2_Total <- summarizeSweep(sweep.res.Visvader_0337_HER2_Total, GT = FALSE)
bcmvn_Visvader_0337_HER2_Total <- find.pK(sweep.stats.Visvader_0337_HER2_Total)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Visvader_0337_HER2_Total[which.max(bcmvn_Visvader_0337_HER2_Total$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Visvader_0337_HER2_Total@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 5% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Visvader dataset = 5k per sample, will use doublet rate for 6,000 cells to err on the side of caution
nExp_poi <- round(0.05*nrow(Visvader_0337_HER2_Total@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Visvader_0337_HER2_Total <- doubletFinder(Visvader_0337_HER2_Total, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Visvader_0337_HER2_Total_Quant <- (Visvader_0337_HER2_Total@meta.data$DF.classification == "Singlet")
Visvader_0337_HER2_Total_Quant_Singlets <- length(Visvader_0337_HER2_Total_Quant[Visvader_0337_HER2_Total_Quant== TRUE])
Visvader_0337_HER2_Total_Quant_Doublets <- length(Visvader_0337_HER2_Total_Quant[Visvader_0337_HER2_Total_Quant== FALSE])
Visvader_0337_HER2_Total_Quant_Doublets_Percent <- Visvader_0337_HER2_Total_Quant_Doublets /  (Visvader_0337_HER2_Total_Quant_Doublets + Visvader_0337_HER2_Total_Quant_Singlets) * 100
Visvader_0337_HER2_Total_Quant <- as.data.frame(c(Visvader_0337_HER2_Total_Quant_Singlets, Visvader_0337_HER2_Total_Quant_Doublets, Visvader_0337_HER2_Total_Quant_Doublets_Percent))
colnames(Visvader_0337_HER2_Total_Quant) <- c("Visvader_0337_HER2_Total_Quant")
rownames(Visvader_0337_HER2_Total_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Visvader_0337_HER2_Total_Quant, file = "/R/R_Visvader/Visvader_Doublet_Tables/Visvader_0337_HER2_Total_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Visvader_0337_HER2_Total_Singlets <- subset(Visvader_0337_HER2_Total, cells=rownames(Visvader_0337_HER2_Total@meta.data)[which(Visvader_0337_HER2_Total@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Visvader_0337_HER2_Total_Singlets, file = "/R/R_Visvader/Visvader_RDS/RDS_Total_Singlets/Visvader_0337_HER2_Total_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Visvader_0337_HER2_Total)
rm(Visvader_0337_HER2_Total_Singlets)
rm(Visvader_0337_HER2_Total_Quant)
rm(Visvader_0337_HER2_Total_Quant_Singlets)
rm(Visvader_0337_HER2_Total_Quant_Doublets)
rm(Visvader_0337_HER2_Total_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Visvader_0337_HER2_Total)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Visvader_0337_HER2_Total)
rm(sweep.stats.Visvader_0337_HER2_Total)
gc()





################################################################################################
########################                   TNBC Tumors                  ########################
################################################################################################

################################################################################################
########################              Visvader_106_TNBC_Total             ######################
################################################################################################

# Load Data
Visvader_106_TNBC_Total <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Total/Visvader_106_TNBC_Total.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Visvader_106_TNBC_Total <- paramSweep(Visvader_106_TNBC_Total, PCs = 1:10, sct = FALSE)
sweep.stats.Visvader_106_TNBC_Total <- summarizeSweep(sweep.res.Visvader_106_TNBC_Total, GT = FALSE)
bcmvn_Visvader_106_TNBC_Total <- find.pK(sweep.stats.Visvader_106_TNBC_Total)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Visvader_106_TNBC_Total[which.max(bcmvn_Visvader_106_TNBC_Total$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Visvader_106_TNBC_Total@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 5% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Visvader dataset = 5k per sample, will use doublet rate for 6,000 cells to err on the side of caution
nExp_poi <- round(0.05*nrow(Visvader_106_TNBC_Total@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Visvader_106_TNBC_Total <- doubletFinder(Visvader_106_TNBC_Total, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Visvader_106_TNBC_Total_Quant <- (Visvader_106_TNBC_Total@meta.data$DF.classification == "Singlet")
Visvader_106_TNBC_Total_Quant_Singlets <- length(Visvader_106_TNBC_Total_Quant[Visvader_106_TNBC_Total_Quant== TRUE])
Visvader_106_TNBC_Total_Quant_Doublets <- length(Visvader_106_TNBC_Total_Quant[Visvader_106_TNBC_Total_Quant== FALSE])
Visvader_106_TNBC_Total_Quant_Doublets_Percent <- Visvader_106_TNBC_Total_Quant_Doublets /  (Visvader_106_TNBC_Total_Quant_Doublets + Visvader_106_TNBC_Total_Quant_Singlets) * 100
Visvader_106_TNBC_Total_Quant <- as.data.frame(c(Visvader_106_TNBC_Total_Quant_Singlets, Visvader_106_TNBC_Total_Quant_Doublets, Visvader_106_TNBC_Total_Quant_Doublets_Percent))
colnames(Visvader_106_TNBC_Total_Quant) <- c("Visvader_106_TNBC_Total_Quant")
rownames(Visvader_106_TNBC_Total_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Visvader_106_TNBC_Total_Quant, file = "/R/R_Visvader/Visvader_Doublet_Tables/Visvader_106_TNBC_Total_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Visvader_106_TNBC_Total_Singlets <- subset(Visvader_106_TNBC_Total, cells=rownames(Visvader_106_TNBC_Total@meta.data)[which(Visvader_106_TNBC_Total@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Visvader_106_TNBC_Total_Singlets, file = "/R/R_Visvader/Visvader_RDS/RDS_Total_Singlets/Visvader_106_TNBC_Total_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Visvader_106_TNBC_Total)
rm(Visvader_106_TNBC_Total_Singlets)
rm(Visvader_106_TNBC_Total_Quant)
rm(Visvader_106_TNBC_Total_Quant_Singlets)
rm(Visvader_106_TNBC_Total_Quant_Doublets)
rm(Visvader_106_TNBC_Total_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Visvader_106_TNBC_Total)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Visvader_106_TNBC_Total)
rm(sweep.stats.Visvader_106_TNBC_Total)
gc()



################################################################################################
########################              Visvader_114_TNBC_Total             ######################
################################################################################################

# Load Data
Visvader_114_TNBC_Total <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Total/Visvader_114_TNBC_Total.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Visvader_114_TNBC_Total <- paramSweep(Visvader_114_TNBC_Total, PCs = 1:10, sct = FALSE)
sweep.stats.Visvader_114_TNBC_Total <- summarizeSweep(sweep.res.Visvader_114_TNBC_Total, GT = FALSE)
bcmvn_Visvader_114_TNBC_Total <- find.pK(sweep.stats.Visvader_114_TNBC_Total)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Visvader_114_TNBC_Total[which.max(bcmvn_Visvader_114_TNBC_Total$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Visvader_114_TNBC_Total@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 5% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Visvader dataset = 5k per sample, will use doublet rate for 6,000 cells to err on the side of caution
nExp_poi <- round(0.05*nrow(Visvader_114_TNBC_Total@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Visvader_114_TNBC_Total <- doubletFinder(Visvader_114_TNBC_Total, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Visvader_114_TNBC_Total_Quant <- (Visvader_114_TNBC_Total@meta.data$DF.classification == "Singlet")
Visvader_114_TNBC_Total_Quant_Singlets <- length(Visvader_114_TNBC_Total_Quant[Visvader_114_TNBC_Total_Quant== TRUE])
Visvader_114_TNBC_Total_Quant_Doublets <- length(Visvader_114_TNBC_Total_Quant[Visvader_114_TNBC_Total_Quant== FALSE])
Visvader_114_TNBC_Total_Quant_Doublets_Percent <- Visvader_114_TNBC_Total_Quant_Doublets /  (Visvader_114_TNBC_Total_Quant_Doublets + Visvader_114_TNBC_Total_Quant_Singlets) * 100
Visvader_114_TNBC_Total_Quant <- as.data.frame(c(Visvader_114_TNBC_Total_Quant_Singlets, Visvader_114_TNBC_Total_Quant_Doublets, Visvader_114_TNBC_Total_Quant_Doublets_Percent))
colnames(Visvader_114_TNBC_Total_Quant) <- c("Visvader_114_TNBC_Total_Quant")
rownames(Visvader_114_TNBC_Total_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Visvader_114_TNBC_Total_Quant, file = "/R/R_Visvader/Visvader_Doublet_Tables/Visvader_114_TNBC_Total_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Visvader_114_TNBC_Total_Singlets <- subset(Visvader_114_TNBC_Total, cells=rownames(Visvader_114_TNBC_Total@meta.data)[which(Visvader_114_TNBC_Total@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Visvader_114_TNBC_Total_Singlets, file = "/R/R_Visvader/Visvader_RDS/RDS_Total_Singlets/Visvader_114_TNBC_Total_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Visvader_114_TNBC_Total)
rm(Visvader_114_TNBC_Total_Singlets)
rm(Visvader_114_TNBC_Total_Quant)
rm(Visvader_114_TNBC_Total_Quant_Singlets)
rm(Visvader_114_TNBC_Total_Quant_Doublets)
rm(Visvader_114_TNBC_Total_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Visvader_114_TNBC_Total)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Visvader_114_TNBC_Total)
rm(sweep.stats.Visvader_114_TNBC_Total)
gc()



################################################################################################
########################              Visvader_126_TNBC_Total             ######################
################################################################################################

# Load Data
Visvader_126_TNBC_Total <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Total/Visvader_126_TNBC_Total.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Visvader_126_TNBC_Total <- paramSweep(Visvader_126_TNBC_Total, PCs = 1:10, sct = FALSE)
sweep.stats.Visvader_126_TNBC_Total <- summarizeSweep(sweep.res.Visvader_126_TNBC_Total, GT = FALSE)
bcmvn_Visvader_126_TNBC_Total <- find.pK(sweep.stats.Visvader_126_TNBC_Total)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Visvader_126_TNBC_Total[which.max(bcmvn_Visvader_126_TNBC_Total$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Visvader_126_TNBC_Total@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 5% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Visvader dataset = 5k per sample, will use doublet rate for 6,000 cells to err on the side of caution
nExp_poi <- round(0.05*nrow(Visvader_126_TNBC_Total@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Visvader_126_TNBC_Total <- doubletFinder(Visvader_126_TNBC_Total, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Visvader_126_TNBC_Total_Quant <- (Visvader_126_TNBC_Total@meta.data$DF.classification == "Singlet")
Visvader_126_TNBC_Total_Quant_Singlets <- length(Visvader_126_TNBC_Total_Quant[Visvader_126_TNBC_Total_Quant== TRUE])
Visvader_126_TNBC_Total_Quant_Doublets <- length(Visvader_126_TNBC_Total_Quant[Visvader_126_TNBC_Total_Quant== FALSE])
Visvader_126_TNBC_Total_Quant_Doublets_Percent <- Visvader_126_TNBC_Total_Quant_Doublets /  (Visvader_126_TNBC_Total_Quant_Doublets + Visvader_126_TNBC_Total_Quant_Singlets) * 100
Visvader_126_TNBC_Total_Quant <- as.data.frame(c(Visvader_126_TNBC_Total_Quant_Singlets, Visvader_126_TNBC_Total_Quant_Doublets, Visvader_126_TNBC_Total_Quant_Doublets_Percent))
colnames(Visvader_126_TNBC_Total_Quant) <- c("Visvader_126_TNBC_Total_Quant")
rownames(Visvader_126_TNBC_Total_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Visvader_126_TNBC_Total_Quant, file = "/R/R_Visvader/Visvader_Doublet_Tables/Visvader_126_TNBC_Total_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Visvader_126_TNBC_Total_Singlets <- subset(Visvader_126_TNBC_Total, cells=rownames(Visvader_126_TNBC_Total@meta.data)[which(Visvader_126_TNBC_Total@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Visvader_126_TNBC_Total_Singlets, file = "/R/R_Visvader/Visvader_RDS/RDS_Total_Singlets/Visvader_126_TNBC_Total_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Visvader_126_TNBC_Total)
rm(Visvader_126_TNBC_Total_Singlets)
rm(Visvader_126_TNBC_Total_Quant)
rm(Visvader_126_TNBC_Total_Quant_Singlets)
rm(Visvader_126_TNBC_Total_Quant_Doublets)
rm(Visvader_126_TNBC_Total_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Visvader_126_TNBC_Total)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Visvader_126_TNBC_Total)
rm(sweep.stats.Visvader_126_TNBC_Total)
gc()



################################################################################################
########################              Visvader_0131_TNBC_Total             #####################
################################################################################################

# Load Data
Visvader_0131_TNBC_Total <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Total/Visvader_0131_TNBC_Total.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Visvader_0131_TNBC_Total <- paramSweep(Visvader_0131_TNBC_Total, PCs = 1:10, sct = FALSE)
sweep.stats.Visvader_0131_TNBC_Total <- summarizeSweep(sweep.res.Visvader_0131_TNBC_Total, GT = FALSE)
bcmvn_Visvader_0131_TNBC_Total <- find.pK(sweep.stats.Visvader_0131_TNBC_Total)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Visvader_0131_TNBC_Total[which.max(bcmvn_Visvader_0131_TNBC_Total$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Visvader_0131_TNBC_Total@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 5% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Visvader dataset = 5k per sample, will use doublet rate for 6,000 cells to err on the side of caution
nExp_poi <- round(0.05*nrow(Visvader_0131_TNBC_Total@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Visvader_0131_TNBC_Total <- doubletFinder(Visvader_0131_TNBC_Total, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Visvader_0131_TNBC_Total_Quant <- (Visvader_0131_TNBC_Total@meta.data$DF.classification == "Singlet")
Visvader_0131_TNBC_Total_Quant_Singlets <- length(Visvader_0131_TNBC_Total_Quant[Visvader_0131_TNBC_Total_Quant== TRUE])
Visvader_0131_TNBC_Total_Quant_Doublets <- length(Visvader_0131_TNBC_Total_Quant[Visvader_0131_TNBC_Total_Quant== FALSE])
Visvader_0131_TNBC_Total_Quant_Doublets_Percent <- Visvader_0131_TNBC_Total_Quant_Doublets /  (Visvader_0131_TNBC_Total_Quant_Doublets + Visvader_0131_TNBC_Total_Quant_Singlets) * 100
Visvader_0131_TNBC_Total_Quant <- as.data.frame(c(Visvader_0131_TNBC_Total_Quant_Singlets, Visvader_0131_TNBC_Total_Quant_Doublets, Visvader_0131_TNBC_Total_Quant_Doublets_Percent))
colnames(Visvader_0131_TNBC_Total_Quant) <- c("Visvader_0131_TNBC_Total_Quant")
rownames(Visvader_0131_TNBC_Total_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Visvader_0131_TNBC_Total_Quant, file = "/R/R_Visvader/Visvader_Doublet_Tables/Visvader_0131_TNBC_Total_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Visvader_0131_TNBC_Total_Singlets <- subset(Visvader_0131_TNBC_Total, cells=rownames(Visvader_0131_TNBC_Total@meta.data)[which(Visvader_0131_TNBC_Total@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Visvader_0131_TNBC_Total_Singlets, file = "/R/R_Visvader/Visvader_RDS/RDS_Total_Singlets/Visvader_0131_TNBC_Total_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Visvader_0131_TNBC_Total)
rm(Visvader_0131_TNBC_Total_Singlets)
rm(Visvader_0131_TNBC_Total_Quant)
rm(Visvader_0131_TNBC_Total_Quant_Singlets)
rm(Visvader_0131_TNBC_Total_Quant_Doublets)
rm(Visvader_0131_TNBC_Total_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Visvader_0131_TNBC_Total)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Visvader_0131_TNBC_Total)
rm(sweep.stats.Visvader_0131_TNBC_Total)
gc()



################################################################################################
########################              Visvader_135_TNBC_Total             ######################
################################################################################################

# Load Data
Visvader_135_TNBC_Total <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Total/Visvader_135_TNBC_Total.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Visvader_135_TNBC_Total <- paramSweep(Visvader_135_TNBC_Total, PCs = 1:10, sct = FALSE)
sweep.stats.Visvader_135_TNBC_Total <- summarizeSweep(sweep.res.Visvader_135_TNBC_Total, GT = FALSE)
bcmvn_Visvader_135_TNBC_Total <- find.pK(sweep.stats.Visvader_135_TNBC_Total)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Visvader_135_TNBC_Total[which.max(bcmvn_Visvader_135_TNBC_Total$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Visvader_135_TNBC_Total@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 5% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Visvader dataset = 5k per sample, will use doublet rate for 6,000 cells to err on the side of caution
nExp_poi <- round(0.05*nrow(Visvader_135_TNBC_Total@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Visvader_135_TNBC_Total <- doubletFinder(Visvader_135_TNBC_Total, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Visvader_135_TNBC_Total_Quant <- (Visvader_135_TNBC_Total@meta.data$DF.classification == "Singlet")
Visvader_135_TNBC_Total_Quant_Singlets <- length(Visvader_135_TNBC_Total_Quant[Visvader_135_TNBC_Total_Quant== TRUE])
Visvader_135_TNBC_Total_Quant_Doublets <- length(Visvader_135_TNBC_Total_Quant[Visvader_135_TNBC_Total_Quant== FALSE])
Visvader_135_TNBC_Total_Quant_Doublets_Percent <- Visvader_135_TNBC_Total_Quant_Doublets /  (Visvader_135_TNBC_Total_Quant_Doublets + Visvader_135_TNBC_Total_Quant_Singlets) * 100
Visvader_135_TNBC_Total_Quant <- as.data.frame(c(Visvader_135_TNBC_Total_Quant_Singlets, Visvader_135_TNBC_Total_Quant_Doublets, Visvader_135_TNBC_Total_Quant_Doublets_Percent))
colnames(Visvader_135_TNBC_Total_Quant) <- c("Visvader_135_TNBC_Total_Quant")
rownames(Visvader_135_TNBC_Total_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Visvader_135_TNBC_Total_Quant, file = "/R/R_Visvader/Visvader_Doublet_Tables/Visvader_135_TNBC_Total_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Visvader_135_TNBC_Total_Singlets <- subset(Visvader_135_TNBC_Total, cells=rownames(Visvader_135_TNBC_Total@meta.data)[which(Visvader_135_TNBC_Total@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Visvader_135_TNBC_Total_Singlets, file = "/R/R_Visvader/Visvader_RDS/RDS_Total_Singlets/Visvader_135_TNBC_Total_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Visvader_135_TNBC_Total)
rm(Visvader_135_TNBC_Total_Singlets)
rm(Visvader_135_TNBC_Total_Quant)
rm(Visvader_135_TNBC_Total_Quant_Singlets)
rm(Visvader_135_TNBC_Total_Quant_Doublets)
rm(Visvader_135_TNBC_Total_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Visvader_135_TNBC_Total)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Visvader_135_TNBC_Total)
rm(sweep.stats.Visvader_135_TNBC_Total)
gc()



################################################################################################
########################              Visvader_0177_TNBC_Total             #####################
################################################################################################

# Load Data
Visvader_0177_TNBC_Total <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Total/Visvader_0177_TNBC_Total.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Visvader_0177_TNBC_Total <- paramSweep(Visvader_0177_TNBC_Total, PCs = 1:10, sct = FALSE)
sweep.stats.Visvader_0177_TNBC_Total <- summarizeSweep(sweep.res.Visvader_0177_TNBC_Total, GT = FALSE)
bcmvn_Visvader_0177_TNBC_Total <- find.pK(sweep.stats.Visvader_0177_TNBC_Total)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Visvader_0177_TNBC_Total[which.max(bcmvn_Visvader_0177_TNBC_Total$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Visvader_0177_TNBC_Total@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 5% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Visvader dataset = 5k per sample, will use doublet rate for 6,000 cells to err on the side of caution
nExp_poi <- round(0.05*nrow(Visvader_0177_TNBC_Total@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Visvader_0177_TNBC_Total <- doubletFinder(Visvader_0177_TNBC_Total, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Visvader_0177_TNBC_Total_Quant <- (Visvader_0177_TNBC_Total@meta.data$DF.classification == "Singlet")
Visvader_0177_TNBC_Total_Quant_Singlets <- length(Visvader_0177_TNBC_Total_Quant[Visvader_0177_TNBC_Total_Quant== TRUE])
Visvader_0177_TNBC_Total_Quant_Doublets <- length(Visvader_0177_TNBC_Total_Quant[Visvader_0177_TNBC_Total_Quant== FALSE])
Visvader_0177_TNBC_Total_Quant_Doublets_Percent <- Visvader_0177_TNBC_Total_Quant_Doublets /  (Visvader_0177_TNBC_Total_Quant_Doublets + Visvader_0177_TNBC_Total_Quant_Singlets) * 100
Visvader_0177_TNBC_Total_Quant <- as.data.frame(c(Visvader_0177_TNBC_Total_Quant_Singlets, Visvader_0177_TNBC_Total_Quant_Doublets, Visvader_0177_TNBC_Total_Quant_Doublets_Percent))
colnames(Visvader_0177_TNBC_Total_Quant) <- c("Visvader_0177_TNBC_Total_Quant")
rownames(Visvader_0177_TNBC_Total_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Visvader_0177_TNBC_Total_Quant, file = "/R/R_Visvader/Visvader_Doublet_Tables/Visvader_0177_TNBC_Total_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Visvader_0177_TNBC_Total_Singlets <- subset(Visvader_0177_TNBC_Total, cells=rownames(Visvader_0177_TNBC_Total@meta.data)[which(Visvader_0177_TNBC_Total@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Visvader_0177_TNBC_Total_Singlets, file = "/R/R_Visvader/Visvader_RDS/RDS_Total_Singlets/Visvader_0177_TNBC_Total_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Visvader_0177_TNBC_Total)
rm(Visvader_0177_TNBC_Total_Singlets)
rm(Visvader_0177_TNBC_Total_Quant)
rm(Visvader_0177_TNBC_Total_Quant_Singlets)
rm(Visvader_0177_TNBC_Total_Quant_Doublets)
rm(Visvader_0177_TNBC_Total_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Visvader_0177_TNBC_Total)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Visvader_0177_TNBC_Total)
rm(sweep.stats.Visvader_0177_TNBC_Total)
gc()



################################################################################################
########################              Visvader_0554_TNBC_Total             #####################
################################################################################################

# Load Data
Visvader_0554_TNBC_Total <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Total/Visvader_0554_TNBC_Total.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Visvader_0554_TNBC_Total <- paramSweep(Visvader_0554_TNBC_Total, PCs = 1:10, sct = FALSE)
sweep.stats.Visvader_0554_TNBC_Total <- summarizeSweep(sweep.res.Visvader_0554_TNBC_Total, GT = FALSE)
bcmvn_Visvader_0554_TNBC_Total <- find.pK(sweep.stats.Visvader_0554_TNBC_Total)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Visvader_0554_TNBC_Total[which.max(bcmvn_Visvader_0554_TNBC_Total$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Visvader_0554_TNBC_Total@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 5% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Visvader dataset = 5k per sample, will use doublet rate for 6,000 cells to err on the side of caution
nExp_poi <- round(0.05*nrow(Visvader_0554_TNBC_Total@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Visvader_0554_TNBC_Total <- doubletFinder(Visvader_0554_TNBC_Total, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Visvader_0554_TNBC_Total_Quant <- (Visvader_0554_TNBC_Total@meta.data$DF.classification == "Singlet")
Visvader_0554_TNBC_Total_Quant_Singlets <- length(Visvader_0554_TNBC_Total_Quant[Visvader_0554_TNBC_Total_Quant== TRUE])
Visvader_0554_TNBC_Total_Quant_Doublets <- length(Visvader_0554_TNBC_Total_Quant[Visvader_0554_TNBC_Total_Quant== FALSE])
Visvader_0554_TNBC_Total_Quant_Doublets_Percent <- Visvader_0554_TNBC_Total_Quant_Doublets /  (Visvader_0554_TNBC_Total_Quant_Doublets + Visvader_0554_TNBC_Total_Quant_Singlets) * 100
Visvader_0554_TNBC_Total_Quant <- as.data.frame(c(Visvader_0554_TNBC_Total_Quant_Singlets, Visvader_0554_TNBC_Total_Quant_Doublets, Visvader_0554_TNBC_Total_Quant_Doublets_Percent))
colnames(Visvader_0554_TNBC_Total_Quant) <- c("Visvader_0554_TNBC_Total_Quant")
rownames(Visvader_0554_TNBC_Total_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Visvader_0554_TNBC_Total_Quant, file = "/R/R_Visvader/Visvader_Doublet_Tables/Visvader_0554_TNBC_Total_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Visvader_0554_TNBC_Total_Singlets <- subset(Visvader_0554_TNBC_Total, cells=rownames(Visvader_0554_TNBC_Total@meta.data)[which(Visvader_0554_TNBC_Total@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Visvader_0554_TNBC_Total_Singlets, file = "/R/R_Visvader/Visvader_RDS/RDS_Total_Singlets/Visvader_0554_TNBC_Total_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Visvader_0554_TNBC_Total)
rm(Visvader_0554_TNBC_Total_Singlets)
rm(Visvader_0554_TNBC_Total_Quant)
rm(Visvader_0554_TNBC_Total_Quant_Singlets)
rm(Visvader_0554_TNBC_Total_Quant_Doublets)
rm(Visvader_0554_TNBC_Total_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Visvader_0554_TNBC_Total)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Visvader_0554_TNBC_Total)
rm(sweep.stats.Visvader_0554_TNBC_Total)
gc()



################################################################################################
########################              Visvader_4031_TNBC_Total             #####################
################################################################################################

# Load Data
Visvader_4031_TNBC_Total <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Total/Visvader_4031_TNBC_Total.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Visvader_4031_TNBC_Total <- paramSweep(Visvader_4031_TNBC_Total, PCs = 1:10, sct = FALSE)
sweep.stats.Visvader_4031_TNBC_Total <- summarizeSweep(sweep.res.Visvader_4031_TNBC_Total, GT = FALSE)
bcmvn_Visvader_4031_TNBC_Total <- find.pK(sweep.stats.Visvader_4031_TNBC_Total)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Visvader_4031_TNBC_Total[which.max(bcmvn_Visvader_4031_TNBC_Total$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Visvader_4031_TNBC_Total@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 5% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Visvader dataset = 5k per sample, will use doublet rate for 6,000 cells to err on the side of caution
nExp_poi <- round(0.05*nrow(Visvader_4031_TNBC_Total@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Visvader_4031_TNBC_Total <- doubletFinder(Visvader_4031_TNBC_Total, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Visvader_4031_TNBC_Total_Quant <- (Visvader_4031_TNBC_Total@meta.data$DF.classification == "Singlet")
Visvader_4031_TNBC_Total_Quant_Singlets <- length(Visvader_4031_TNBC_Total_Quant[Visvader_4031_TNBC_Total_Quant== TRUE])
Visvader_4031_TNBC_Total_Quant_Doublets <- length(Visvader_4031_TNBC_Total_Quant[Visvader_4031_TNBC_Total_Quant== FALSE])
Visvader_4031_TNBC_Total_Quant_Doublets_Percent <- Visvader_4031_TNBC_Total_Quant_Doublets /  (Visvader_4031_TNBC_Total_Quant_Doublets + Visvader_4031_TNBC_Total_Quant_Singlets) * 100
Visvader_4031_TNBC_Total_Quant <- as.data.frame(c(Visvader_4031_TNBC_Total_Quant_Singlets, Visvader_4031_TNBC_Total_Quant_Doublets, Visvader_4031_TNBC_Total_Quant_Doublets_Percent))
colnames(Visvader_4031_TNBC_Total_Quant) <- c("Visvader_4031_TNBC_Total_Quant")
rownames(Visvader_4031_TNBC_Total_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Visvader_4031_TNBC_Total_Quant, file = "/R/R_Visvader/Visvader_Doublet_Tables/Visvader_4031_TNBC_Total_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Visvader_4031_TNBC_Total_Singlets <- subset(Visvader_4031_TNBC_Total, cells=rownames(Visvader_4031_TNBC_Total@meta.data)[which(Visvader_4031_TNBC_Total@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Visvader_4031_TNBC_Total_Singlets, file = "/R/R_Visvader/Visvader_RDS/RDS_Total_Singlets/Visvader_4031_TNBC_Total_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Visvader_4031_TNBC_Total)
rm(Visvader_4031_TNBC_Total_Singlets)
rm(Visvader_4031_TNBC_Total_Quant)
rm(Visvader_4031_TNBC_Total_Quant_Singlets)
rm(Visvader_4031_TNBC_Total_Quant_Doublets)
rm(Visvader_4031_TNBC_Total_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Visvader_4031_TNBC_Total)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Visvader_4031_TNBC_Total)
rm(sweep.stats.Visvader_4031_TNBC_Total)
gc()



