#####################################################################################################################
#                        Swarbrick (Wu et al) Dataset Breast Tumor Analysis Steps                                   #
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
library(Matrix)

# Employed workaround to load sparse matrix into Seurat
# See: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/matrices
# Files from the GEO entry (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176078) were renamed and zipped via gzip in terminal

#########################################################
############# This is ER Tumor Batch # ##################
#########################################################


#########################################################
#########      Swarbrick_CID3941_ER_Total       ######### 
#########################################################
# Load in matrix
matrix_dir = "/R/R_Swarbrick/Swarbrick_Input/ER/Swarbrick_CID3941/"
barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "features.tsv.gz")
matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")
mat <- readMM(file = matrix.path)
feature.names = read.delim(features.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
colnames(mat) = barcode.names$V1
rownames(mat) = feature.names$V1

# Read Metadata
metadata <- read.csv(file =  "/R/R_Swarbrick/Swarbrick_Input/ER/Swarbrick_CID3941/metadata.csv", header = TRUE, row.names = 1)

# Create Seurat Object
Swarbrick_CID3941_ER_Total <- CreateSeuratObject(mat, project="Swarbrick_CID3941_ER_Total", min.cells = 3, min.features = 200, meta.data = metadata)
Swarbrick_CID3941_ER_Total

# Remove matrix files
rm(mat)
rm(feature.names)
rm(barcode.names)
rm(metadata)

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Swarbrick_CID3941_ER_Total[["percent.mt"]] <- PercentageFeatureSet(Swarbrick_CID3941_ER_Total, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Swarbrick_CID3941_ER_Total, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Swarbrick_CID3941_ER_Total, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Swarbrick_CID3941_ER_Total, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Swarbrick_CID3941_ER_Total <- subset(Swarbrick_CID3941_ER_Total, subset = nFeature_RNA > 200 & percent.mt < 20)

# Normalize data
Swarbrick_CID3941_ER_Total <- NormalizeData(Swarbrick_CID3941_ER_Total, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Swarbrick_CID3941_ER_Total <- FindVariableFeatures(Swarbrick_CID3941_ER_Total, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Swarbrick_CID3941_ER_Total)
Swarbrick_CID3941_ER_Total <- ScaleData(Swarbrick_CID3941_ER_Total, features = all.genes)

# Perform linear dimensional reduction
Swarbrick_CID3941_ER_Total <- RunPCA(Swarbrick_CID3941_ER_Total, features = VariableFeatures(object = Swarbrick_CID3941_ER_Total))

# Examine and visualize PCA results a few different ways
print(Swarbrick_CID3941_ER_Total[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Swarbrick_CID3941_ER_Total, dims = 1:2, reduction = "pca")
DimPlot(Swarbrick_CID3941_ER_Total, reduction = "pca")

# Visualize PCs
ElbowPlot(Swarbrick_CID3941_ER_Total)

# Cluster cells
Swarbrick_CID3941_ER_Total <- FindNeighbors(Swarbrick_CID3941_ER_Total, dims = 1:20)
Swarbrick_CID3941_ER_Total <- FindClusters(Swarbrick_CID3941_ER_Total, resolution = 0.5)

# Run non-linear dimensional reduction
Swarbrick_CID3941_ER_Total <- RunUMAP(Swarbrick_CID3941_ER_Total, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Swarbrick_CID3941_ER_Total, reduction = "umap")

# Save RDS
saveRDS(Swarbrick_CID3941_ER_Total, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Total/Swarbrick_CID3941_ER_Total.rds")
rm(Swarbrick_CID3941_ER_Total)


#########################################################
#########     Swarbrick_CID3948_ER_Total        ######### 
#########################################################
# Load in matrix
matrix_dir = "/R/R_Swarbrick/Swarbrick_Input/ER/Swarbrick_CID3948/"
barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "features.tsv.gz")
matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")
mat <- readMM(file = matrix.path)
feature.names = read.delim(features.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
colnames(mat) = barcode.names$V1
rownames(mat) = feature.names$V1

# Read Metadata
metadata <- read.csv(file =  "/R/R_Swarbrick/Swarbrick_Input/ER/Swarbrick_CID3948/metadata.csv", header = TRUE, row.names = 1)

# Create Seurat Object
Swarbrick_CID3948_ER_Total <- CreateSeuratObject(mat, project="Swarbrick_CID3948_ER_Total", min.cells = 3, min.features = 200, meta.data = metadata)
Swarbrick_CID3948_ER_Total

# Remove matrix files
rm(mat)
rm(feature.names)
rm(barcode.names)
rm(metadata)

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Swarbrick_CID3948_ER_Total[["percent.mt"]] <- PercentageFeatureSet(Swarbrick_CID3948_ER_Total, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Swarbrick_CID3948_ER_Total, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Swarbrick_CID3948_ER_Total, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Swarbrick_CID3948_ER_Total, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Swarbrick_CID3948_ER_Total <- subset(Swarbrick_CID3948_ER_Total, subset = nFeature_RNA > 200 & percent.mt < 20)

# Normalize data
Swarbrick_CID3948_ER_Total <- NormalizeData(Swarbrick_CID3948_ER_Total, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Swarbrick_CID3948_ER_Total <- FindVariableFeatures(Swarbrick_CID3948_ER_Total, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Swarbrick_CID3948_ER_Total)
Swarbrick_CID3948_ER_Total <- ScaleData(Swarbrick_CID3948_ER_Total, features = all.genes)

# Perform linear dimensional reduction
Swarbrick_CID3948_ER_Total <- RunPCA(Swarbrick_CID3948_ER_Total, features = VariableFeatures(object = Swarbrick_CID3948_ER_Total))

# Examine and visualize PCA results a few different ways
print(Swarbrick_CID3948_ER_Total[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Swarbrick_CID3948_ER_Total, dims = 1:2, reduction = "pca")
DimPlot(Swarbrick_CID3948_ER_Total, reduction = "pca")

# Visualize PCs
ElbowPlot(Swarbrick_CID3948_ER_Total)

# Cluster cells
Swarbrick_CID3948_ER_Total <- FindNeighbors(Swarbrick_CID3948_ER_Total, dims = 1:20)
Swarbrick_CID3948_ER_Total <- FindClusters(Swarbrick_CID3948_ER_Total, resolution = 0.5)

# Run non-linear dimensional reduction
Swarbrick_CID3948_ER_Total <- RunUMAP(Swarbrick_CID3948_ER_Total, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Swarbrick_CID3948_ER_Total, reduction = "umap")

# Save RDS
saveRDS(Swarbrick_CID3948_ER_Total, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Total/Swarbrick_CID3948_ER_Total.rds")
rm(Swarbrick_CID3948_ER_Total)



#########################################################
#########      Swarbrick_CID4040_ER_Total       ######### 
#########################################################
# Load in matrix
matrix_dir = "/R/R_Swarbrick/Swarbrick_Input/ER/Swarbrick_CID4040/"
barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "features.tsv.gz")
matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")
mat <- readMM(file = matrix.path)
feature.names = read.delim(features.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
colnames(mat) = barcode.names$V1
rownames(mat) = feature.names$V1

# Read Metadata
metadata <- read.csv(file =  "/R/R_Swarbrick/Swarbrick_Input/ER/Swarbrick_CID4040/metadata.csv", header = TRUE, row.names = 1)

# Create Seurat Object
Swarbrick_CID4040_ER_Total <- CreateSeuratObject(mat, project="Swarbrick_CID4040_ER_Total", min.cells = 3, min.features = 200, meta.data = metadata)
Swarbrick_CID4040_ER_Total

# Remove matrix files
rm(mat)
rm(feature.names)
rm(barcode.names)
rm(metadata)

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Swarbrick_CID4040_ER_Total[["percent.mt"]] <- PercentageFeatureSet(Swarbrick_CID4040_ER_Total, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Swarbrick_CID4040_ER_Total, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Swarbrick_CID4040_ER_Total, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Swarbrick_CID4040_ER_Total, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Swarbrick_CID4040_ER_Total <- subset(Swarbrick_CID4040_ER_Total, subset = nFeature_RNA > 200 & percent.mt < 20)

# Normalize data
Swarbrick_CID4040_ER_Total <- NormalizeData(Swarbrick_CID4040_ER_Total, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Swarbrick_CID4040_ER_Total <- FindVariableFeatures(Swarbrick_CID4040_ER_Total, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Swarbrick_CID4040_ER_Total)
Swarbrick_CID4040_ER_Total <- ScaleData(Swarbrick_CID4040_ER_Total, features = all.genes)

# Perform linear dimensional reduction
Swarbrick_CID4040_ER_Total <- RunPCA(Swarbrick_CID4040_ER_Total, features = VariableFeatures(object = Swarbrick_CID4040_ER_Total))

# Examine and visualize PCA results a few different ways
print(Swarbrick_CID4040_ER_Total[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Swarbrick_CID4040_ER_Total, dims = 1:2, reduction = "pca")
DimPlot(Swarbrick_CID4040_ER_Total, reduction = "pca")

# Visualize PCs
ElbowPlot(Swarbrick_CID4040_ER_Total)

# Cluster cells
Swarbrick_CID4040_ER_Total <- FindNeighbors(Swarbrick_CID4040_ER_Total, dims = 1:20)
Swarbrick_CID4040_ER_Total <- FindClusters(Swarbrick_CID4040_ER_Total, resolution = 0.5)

# Run non-linear dimensional reduction
Swarbrick_CID4040_ER_Total <- RunUMAP(Swarbrick_CID4040_ER_Total, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Swarbrick_CID4040_ER_Total, reduction = "umap")

# Save RDS
saveRDS(Swarbrick_CID4040_ER_Total, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Total/Swarbrick_CID4040_ER_Total.rds")
rm(Swarbrick_CID4040_ER_Total)



#########################################################
#########       Swarbrick_CID4067_ER_Total      ######### 
#########################################################
# Load in matrix
matrix_dir = "/R/R_Swarbrick/Swarbrick_Input/ER/Swarbrick_CID4067/"
barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "features.tsv.gz")
matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")
mat <- readMM(file = matrix.path)
feature.names = read.delim(features.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
colnames(mat) = barcode.names$V1
rownames(mat) = feature.names$V1

# Read Metadata
metadata <- read.csv(file =  "/R/R_Swarbrick/Swarbrick_Input/ER/Swarbrick_CID4067/metadata.csv", header = TRUE, row.names = 1)

# Create Seurat Object
Swarbrick_CID4067_ER_Total <- CreateSeuratObject(mat, project="Swarbrick_CID4067_ER_Total", min.cells = 3, min.features = 200, meta.data = metadata)
Swarbrick_CID4067_ER_Total

# Remove matrix files
rm(mat)
rm(feature.names)
rm(barcode.names)
rm(metadata)

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Swarbrick_CID4067_ER_Total[["percent.mt"]] <- PercentageFeatureSet(Swarbrick_CID4067_ER_Total, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Swarbrick_CID4067_ER_Total, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Swarbrick_CID4067_ER_Total, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Swarbrick_CID4067_ER_Total, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Swarbrick_CID4067_ER_Total <- subset(Swarbrick_CID4067_ER_Total, subset = nFeature_RNA > 200 & percent.mt < 20)

# Normalize data
Swarbrick_CID4067_ER_Total <- NormalizeData(Swarbrick_CID4067_ER_Total, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Swarbrick_CID4067_ER_Total <- FindVariableFeatures(Swarbrick_CID4067_ER_Total, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Swarbrick_CID4067_ER_Total)
Swarbrick_CID4067_ER_Total <- ScaleData(Swarbrick_CID4067_ER_Total, features = all.genes)

# Perform linear dimensional reduction
Swarbrick_CID4067_ER_Total <- RunPCA(Swarbrick_CID4067_ER_Total, features = VariableFeatures(object = Swarbrick_CID4067_ER_Total))

# Examine and visualize PCA results a few different ways
print(Swarbrick_CID4067_ER_Total[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Swarbrick_CID4067_ER_Total, dims = 1:2, reduction = "pca")
DimPlot(Swarbrick_CID4067_ER_Total, reduction = "pca")

# Visualize PCs
ElbowPlot(Swarbrick_CID4067_ER_Total)

# Cluster cells
Swarbrick_CID4067_ER_Total <- FindNeighbors(Swarbrick_CID4067_ER_Total, dims = 1:20)
Swarbrick_CID4067_ER_Total <- FindClusters(Swarbrick_CID4067_ER_Total, resolution = 0.5)

# Run non-linear dimensional reduction
Swarbrick_CID4067_ER_Total <- RunUMAP(Swarbrick_CID4067_ER_Total, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Swarbrick_CID4067_ER_Total, reduction = "umap")

# Save RDS
saveRDS(Swarbrick_CID4067_ER_Total, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Total/Swarbrick_CID4067_ER_Total.rds")
rm(Swarbrick_CID4067_ER_Total)



#########################################################
#########     Swarbrick_CID4290A_ER_Total       ######### 
#########################################################
# Load in matrix
matrix_dir = "/R/R_Swarbrick/Swarbrick_Input/ER/Swarbrick_CID4290A/"
barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "features.tsv.gz")
matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")
mat <- readMM(file = matrix.path)
feature.names = read.delim(features.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
colnames(mat) = barcode.names$V1
rownames(mat) = feature.names$V1

# Read Metadata
metadata <- read.csv(file =  "/R/R_Swarbrick/Swarbrick_Input/ER/Swarbrick_CID4290A/metadata.csv", header = TRUE, row.names = 1)

# Create Seurat Object
Swarbrick_CID4290A_ER_Total <- CreateSeuratObject(mat, project="Swarbrick_CID4290A_ER_Total", min.cells = 3, min.features = 200, meta.data = metadata)
Swarbrick_CID4290A_ER_Total

# Remove matrix files
rm(mat)
rm(feature.names)
rm(barcode.names)
rm(metadata)

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Swarbrick_CID4290A_ER_Total[["percent.mt"]] <- PercentageFeatureSet(Swarbrick_CID4290A_ER_Total, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Swarbrick_CID4290A_ER_Total, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Swarbrick_CID4290A_ER_Total, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Swarbrick_CID4290A_ER_Total, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Swarbrick_CID4290A_ER_Total <- subset(Swarbrick_CID4290A_ER_Total, subset = nFeature_RNA > 200 & percent.mt < 20)

# Normalize data
Swarbrick_CID4290A_ER_Total <- NormalizeData(Swarbrick_CID4290A_ER_Total, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Swarbrick_CID4290A_ER_Total <- FindVariableFeatures(Swarbrick_CID4290A_ER_Total, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Swarbrick_CID4290A_ER_Total)
Swarbrick_CID4290A_ER_Total <- ScaleData(Swarbrick_CID4290A_ER_Total, features = all.genes)

# Perform linear dimensional reduction
Swarbrick_CID4290A_ER_Total <- RunPCA(Swarbrick_CID4290A_ER_Total, features = VariableFeatures(object = Swarbrick_CID4290A_ER_Total))

# Examine and visualize PCA results a few different ways
print(Swarbrick_CID4290A_ER_Total[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Swarbrick_CID4290A_ER_Total, dims = 1:2, reduction = "pca")
DimPlot(Swarbrick_CID4290A_ER_Total, reduction = "pca")

# Visualize PCs
ElbowPlot(Swarbrick_CID4290A_ER_Total)

# Cluster cells
Swarbrick_CID4290A_ER_Total <- FindNeighbors(Swarbrick_CID4290A_ER_Total, dims = 1:20)
Swarbrick_CID4290A_ER_Total <- FindClusters(Swarbrick_CID4290A_ER_Total, resolution = 0.5)

# Run non-linear dimensional reduction
Swarbrick_CID4290A_ER_Total <- RunUMAP(Swarbrick_CID4290A_ER_Total, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Swarbrick_CID4290A_ER_Total, reduction = "umap")

# Save RDS
saveRDS(Swarbrick_CID4290A_ER_Total, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Total/Swarbrick_CID4290A_ER_Total.rds")
rm(Swarbrick_CID4290A_ER_Total)



#########################################################
#########     Swarbrick_CID4398_ER_Total        ######### 
#########################################################
# Load in matrix
matrix_dir = "/R/R_Swarbrick/Swarbrick_Input/ER/Swarbrick_CID4398/"
barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "features.tsv.gz")
matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")
mat <- readMM(file = matrix.path)
feature.names = read.delim(features.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
colnames(mat) = barcode.names$V1
rownames(mat) = feature.names$V1

# Read Metadata
metadata <- read.csv(file =  "/R/R_Swarbrick/Swarbrick_Input/ER/Swarbrick_CID4398/metadata.csv", header = TRUE, row.names = 1)

# Create Seurat Object
Swarbrick_CID4398_ER_Total <- CreateSeuratObject(mat, project="Swarbrick_CID4398_ER_Total", min.cells = 3, min.features = 200, meta.data = metadata)
Swarbrick_CID4398_ER_Total

# Remove matrix files
rm(mat)
rm(feature.names)
rm(barcode.names)
rm(metadata)

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Swarbrick_CID4398_ER_Total[["percent.mt"]] <- PercentageFeatureSet(Swarbrick_CID4398_ER_Total, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Swarbrick_CID4398_ER_Total, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Swarbrick_CID4398_ER_Total, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Swarbrick_CID4398_ER_Total, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Swarbrick_CID4398_ER_Total <- subset(Swarbrick_CID4398_ER_Total, subset = nFeature_RNA > 200 & percent.mt < 20)

# Normalize data
Swarbrick_CID4398_ER_Total <- NormalizeData(Swarbrick_CID4398_ER_Total, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Swarbrick_CID4398_ER_Total <- FindVariableFeatures(Swarbrick_CID4398_ER_Total, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Swarbrick_CID4398_ER_Total)
Swarbrick_CID4398_ER_Total <- ScaleData(Swarbrick_CID4398_ER_Total, features = all.genes)

# Perform linear dimensional reduction
Swarbrick_CID4398_ER_Total <- RunPCA(Swarbrick_CID4398_ER_Total, features = VariableFeatures(object = Swarbrick_CID4398_ER_Total))

# Examine and visualize PCA results a few different ways
print(Swarbrick_CID4398_ER_Total[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Swarbrick_CID4398_ER_Total, dims = 1:2, reduction = "pca")
DimPlot(Swarbrick_CID4398_ER_Total, reduction = "pca")

# Visualize PCs
ElbowPlot(Swarbrick_CID4398_ER_Total)

# Cluster cells
Swarbrick_CID4398_ER_Total <- FindNeighbors(Swarbrick_CID4398_ER_Total, dims = 1:20)
Swarbrick_CID4398_ER_Total <- FindClusters(Swarbrick_CID4398_ER_Total, resolution = 0.5)

# Run non-linear dimensional reduction
Swarbrick_CID4398_ER_Total <- RunUMAP(Swarbrick_CID4398_ER_Total, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Swarbrick_CID4398_ER_Total, reduction = "umap")

# Save RDS
saveRDS(Swarbrick_CID4398_ER_Total, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Total/Swarbrick_CID4398_ER_Total.rds")
rm(Swarbrick_CID4398_ER_Total)



#########################################################
#########     Swarbrick_CID4461_ER_Total        ######### 
#########################################################
# Load in matrix
matrix_dir = "/R/R_Swarbrick/Swarbrick_Input/ER/Swarbrick_CID4461/"
barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "features.tsv.gz")
matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")
mat <- readMM(file = matrix.path)
feature.names = read.delim(features.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
colnames(mat) = barcode.names$V1
rownames(mat) = feature.names$V1

# Read Metadata
metadata <- read.csv(file =  "/R/R_Swarbrick/Swarbrick_Input/ER/Swarbrick_CID4461/metadata.csv", header = TRUE, row.names = 1)

# Create Seurat Object
Swarbrick_CID4461_ER_Total <- CreateSeuratObject(mat, project="Swarbrick_CID4461_ER_Total", min.cells = 3, min.features = 200, meta.data = metadata)
Swarbrick_CID4461_ER_Total

# Remove matrix files
rm(mat)
rm(feature.names)
rm(barcode.names)
rm(metadata)

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Swarbrick_CID4461_ER_Total[["percent.mt"]] <- PercentageFeatureSet(Swarbrick_CID4461_ER_Total, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Swarbrick_CID4461_ER_Total, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Swarbrick_CID4461_ER_Total, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Swarbrick_CID4461_ER_Total, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Swarbrick_CID4461_ER_Total <- subset(Swarbrick_CID4461_ER_Total, subset = nFeature_RNA > 200 & percent.mt < 20)

# Normalize data
Swarbrick_CID4461_ER_Total <- NormalizeData(Swarbrick_CID4461_ER_Total, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Swarbrick_CID4461_ER_Total <- FindVariableFeatures(Swarbrick_CID4461_ER_Total, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Swarbrick_CID4461_ER_Total)
Swarbrick_CID4461_ER_Total <- ScaleData(Swarbrick_CID4461_ER_Total, features = all.genes)

# Perform linear dimensional reduction
Swarbrick_CID4461_ER_Total <- RunPCA(Swarbrick_CID4461_ER_Total, features = VariableFeatures(object = Swarbrick_CID4461_ER_Total))

# Examine and visualize PCA results a few different ways
print(Swarbrick_CID4461_ER_Total[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Swarbrick_CID4461_ER_Total, dims = 1:2, reduction = "pca")
DimPlot(Swarbrick_CID4461_ER_Total, reduction = "pca")

# Visualize PCs
ElbowPlot(Swarbrick_CID4461_ER_Total)

# Cluster cells
Swarbrick_CID4461_ER_Total <- FindNeighbors(Swarbrick_CID4461_ER_Total, dims = 1:20)
Swarbrick_CID4461_ER_Total <- FindClusters(Swarbrick_CID4461_ER_Total, resolution = 0.5)

# Run non-linear dimensional reduction
Swarbrick_CID4461_ER_Total <- RunUMAP(Swarbrick_CID4461_ER_Total, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Swarbrick_CID4461_ER_Total, reduction = "umap")

# Save RDS
saveRDS(Swarbrick_CID4461_ER_Total, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Total/Swarbrick_CID4461_ER_Total.rds")
rm(Swarbrick_CID4461_ER_Total)



#########################################################
#########      Swarbrick_CID4463_ER_Total       ######### 
#########################################################
# Load in matrix
matrix_dir = "/R/R_Swarbrick/Swarbrick_Input/ER/Swarbrick_CID4463/"
barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "features.tsv.gz")
matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")
mat <- readMM(file = matrix.path)
feature.names = read.delim(features.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
colnames(mat) = barcode.names$V1
rownames(mat) = feature.names$V1

# Read Metadata
metadata <- read.csv(file =  "/R/R_Swarbrick/Swarbrick_Input/ER/Swarbrick_CID4463/metadata.csv", header = TRUE, row.names = 1)

# Create Seurat Object
Swarbrick_CID4463_ER_Total <- CreateSeuratObject(mat, project="Swarbrick_CID4463_ER_Total", min.cells = 3, min.features = 200, meta.data = metadata)
Swarbrick_CID4463_ER_Total

# Remove matrix files
rm(mat)
rm(feature.names)
rm(barcode.names)
rm(metadata)

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Swarbrick_CID4463_ER_Total[["percent.mt"]] <- PercentageFeatureSet(Swarbrick_CID4463_ER_Total, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Swarbrick_CID4463_ER_Total, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Swarbrick_CID4463_ER_Total, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Swarbrick_CID4463_ER_Total, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Swarbrick_CID4463_ER_Total <- subset(Swarbrick_CID4463_ER_Total, subset = nFeature_RNA > 200 & percent.mt < 20)

# Normalize data
Swarbrick_CID4463_ER_Total <- NormalizeData(Swarbrick_CID4463_ER_Total, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Swarbrick_CID4463_ER_Total <- FindVariableFeatures(Swarbrick_CID4463_ER_Total, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Swarbrick_CID4463_ER_Total)
Swarbrick_CID4463_ER_Total <- ScaleData(Swarbrick_CID4463_ER_Total, features = all.genes)

# Perform linear dimensional reduction
Swarbrick_CID4463_ER_Total <- RunPCA(Swarbrick_CID4463_ER_Total, features = VariableFeatures(object = Swarbrick_CID4463_ER_Total))

# Examine and visualize PCA results a few different ways
print(Swarbrick_CID4463_ER_Total[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Swarbrick_CID4463_ER_Total, dims = 1:2, reduction = "pca")
DimPlot(Swarbrick_CID4463_ER_Total, reduction = "pca")

# Visualize PCs
ElbowPlot(Swarbrick_CID4463_ER_Total)

# Cluster cells
Swarbrick_CID4463_ER_Total <- FindNeighbors(Swarbrick_CID4463_ER_Total, dims = 1:20)
Swarbrick_CID4463_ER_Total <- FindClusters(Swarbrick_CID4463_ER_Total, resolution = 0.5)

# Run non-linear dimensional reduction
Swarbrick_CID4463_ER_Total <- RunUMAP(Swarbrick_CID4463_ER_Total, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Swarbrick_CID4463_ER_Total, reduction = "umap")

# Save RDS
saveRDS(Swarbrick_CID4463_ER_Total, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Total/Swarbrick_CID4463_ER_Total.rds")
rm(Swarbrick_CID4463_ER_Total)



#########################################################
#########     Swarbrick_CID4471_ER_Total        ######### 
#########################################################
# Load in matrix
matrix_dir = "/R/R_Swarbrick/Swarbrick_Input/ER/Swarbrick_CID4471/"
barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "features.tsv.gz")
matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")
mat <- readMM(file = matrix.path)
feature.names = read.delim(features.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
colnames(mat) = barcode.names$V1
rownames(mat) = feature.names$V1

# Read Metadata
metadata <- read.csv(file =  "/R/R_Swarbrick/Swarbrick_Input/ER/Swarbrick_CID4471/metadata.csv", header = TRUE, row.names = 1)

# Create Seurat Object
Swarbrick_CID4471_ER_Total <- CreateSeuratObject(mat, project="Swarbrick_CID4471_ER_Total", min.cells = 3, min.features = 200, meta.data = metadata)
Swarbrick_CID4471_ER_Total

# Remove matrix files
rm(mat)
rm(feature.names)
rm(barcode.names)
rm(metadata)

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Swarbrick_CID4471_ER_Total[["percent.mt"]] <- PercentageFeatureSet(Swarbrick_CID4471_ER_Total, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Swarbrick_CID4471_ER_Total, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Swarbrick_CID4471_ER_Total, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Swarbrick_CID4471_ER_Total, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Swarbrick_CID4471_ER_Total <- subset(Swarbrick_CID4471_ER_Total, subset = nFeature_RNA > 200 & percent.mt < 20)

# Normalize data
Swarbrick_CID4471_ER_Total <- NormalizeData(Swarbrick_CID4471_ER_Total, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Swarbrick_CID4471_ER_Total <- FindVariableFeatures(Swarbrick_CID4471_ER_Total, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Swarbrick_CID4471_ER_Total)
Swarbrick_CID4471_ER_Total <- ScaleData(Swarbrick_CID4471_ER_Total, features = all.genes)

# Perform linear dimensional reduction
Swarbrick_CID4471_ER_Total <- RunPCA(Swarbrick_CID4471_ER_Total, features = VariableFeatures(object = Swarbrick_CID4471_ER_Total))

# Examine and visualize PCA results a few different ways
print(Swarbrick_CID4471_ER_Total[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Swarbrick_CID4471_ER_Total, dims = 1:2, reduction = "pca")
DimPlot(Swarbrick_CID4471_ER_Total, reduction = "pca")

# Visualize PCs
ElbowPlot(Swarbrick_CID4471_ER_Total)

# Cluster cells
Swarbrick_CID4471_ER_Total <- FindNeighbors(Swarbrick_CID4471_ER_Total, dims = 1:20)
Swarbrick_CID4471_ER_Total <- FindClusters(Swarbrick_CID4471_ER_Total, resolution = 0.5)

# Run non-linear dimensional reduction
Swarbrick_CID4471_ER_Total <- RunUMAP(Swarbrick_CID4471_ER_Total, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Swarbrick_CID4471_ER_Total, reduction = "umap")

# Save RDS
saveRDS(Swarbrick_CID4471_ER_Total, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Total/Swarbrick_CID4471_ER_Total.rds")
rm(Swarbrick_CID4471_ER_Total)



#########################################################
#########     Swarbrick_CID4530N_ER_Total       ######### 
#########################################################
# Load in matrix
matrix_dir = "/R/R_Swarbrick/Swarbrick_Input/ER/Swarbrick_CID4530N/"
barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "features.tsv.gz")
matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")
mat <- readMM(file = matrix.path)
feature.names = read.delim(features.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
colnames(mat) = barcode.names$V1
rownames(mat) = feature.names$V1

# Read Metadata
metadata <- read.csv(file =  "/R/R_Swarbrick/Swarbrick_Input/ER/Swarbrick_CID4530N/metadata.csv", header = TRUE, row.names = 1)

# Create Seurat Object
Swarbrick_CID4530N_ER_Total <- CreateSeuratObject(mat, project="Swarbrick_CID4530N_ER_Total", min.cells = 3, min.features = 200, meta.data = metadata)
Swarbrick_CID4530N_ER_Total

# Remove matrix files
rm(mat)
rm(feature.names)
rm(barcode.names)
rm(metadata)

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Swarbrick_CID4530N_ER_Total[["percent.mt"]] <- PercentageFeatureSet(Swarbrick_CID4530N_ER_Total, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Swarbrick_CID4530N_ER_Total, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Swarbrick_CID4530N_ER_Total, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Swarbrick_CID4530N_ER_Total, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Swarbrick_CID4530N_ER_Total <- subset(Swarbrick_CID4530N_ER_Total, subset = nFeature_RNA > 200 & percent.mt < 20)

# Normalize data
Swarbrick_CID4530N_ER_Total <- NormalizeData(Swarbrick_CID4530N_ER_Total, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Swarbrick_CID4530N_ER_Total <- FindVariableFeatures(Swarbrick_CID4530N_ER_Total, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Swarbrick_CID4530N_ER_Total)
Swarbrick_CID4530N_ER_Total <- ScaleData(Swarbrick_CID4530N_ER_Total, features = all.genes)

# Perform linear dimensional reduction
Swarbrick_CID4530N_ER_Total <- RunPCA(Swarbrick_CID4530N_ER_Total, features = VariableFeatures(object = Swarbrick_CID4530N_ER_Total))

# Examine and visualize PCA results a few different ways
print(Swarbrick_CID4530N_ER_Total[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Swarbrick_CID4530N_ER_Total, dims = 1:2, reduction = "pca")
DimPlot(Swarbrick_CID4530N_ER_Total, reduction = "pca")

# Visualize PCs
ElbowPlot(Swarbrick_CID4530N_ER_Total)

# Cluster cells
Swarbrick_CID4530N_ER_Total <- FindNeighbors(Swarbrick_CID4530N_ER_Total, dims = 1:20)
Swarbrick_CID4530N_ER_Total <- FindClusters(Swarbrick_CID4530N_ER_Total, resolution = 0.5)

# Run non-linear dimensional reduction
Swarbrick_CID4530N_ER_Total <- RunUMAP(Swarbrick_CID4530N_ER_Total, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Swarbrick_CID4530N_ER_Total, reduction = "umap")

# Save RDS
saveRDS(Swarbrick_CID4530N_ER_Total, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Total/Swarbrick_CID4530N_ER_Total.rds")
rm(Swarbrick_CID4530N_ER_Total)



#########################################################
#########     Swarbrick_CID4535_ER_Total        ######### 
#########################################################
# Load in matrix
matrix_dir = "/R/R_Swarbrick/Swarbrick_Input/ER/Swarbrick_CID4535/"
barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "features.tsv.gz")
matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")
mat <- readMM(file = matrix.path)
feature.names = read.delim(features.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
colnames(mat) = barcode.names$V1
rownames(mat) = feature.names$V1

# Read Metadata
metadata <- read.csv(file =  "/R/R_Swarbrick/Swarbrick_Input/ER/Swarbrick_CID4535/metadata.csv", header = TRUE, row.names = 1)

# Create Seurat Object
Swarbrick_CID4535_ER_Total <- CreateSeuratObject(mat, project="Swarbrick_CID4535_ER_Total", min.cells = 3, min.features = 200, meta.data = metadata)
Swarbrick_CID4535_ER_Total

# Remove matrix files
rm(mat)
rm(feature.names)
rm(barcode.names)
rm(metadata)

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Swarbrick_CID4535_ER_Total[["percent.mt"]] <- PercentageFeatureSet(Swarbrick_CID4535_ER_Total, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Swarbrick_CID4535_ER_Total, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Swarbrick_CID4535_ER_Total, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Swarbrick_CID4535_ER_Total, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Swarbrick_CID4535_ER_Total <- subset(Swarbrick_CID4535_ER_Total, subset = nFeature_RNA > 200 & percent.mt < 20)

# Normalize data
Swarbrick_CID4535_ER_Total <- NormalizeData(Swarbrick_CID4535_ER_Total, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Swarbrick_CID4535_ER_Total <- FindVariableFeatures(Swarbrick_CID4535_ER_Total, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Swarbrick_CID4535_ER_Total)
Swarbrick_CID4535_ER_Total <- ScaleData(Swarbrick_CID4535_ER_Total, features = all.genes)

# Perform linear dimensional reduction
Swarbrick_CID4535_ER_Total <- RunPCA(Swarbrick_CID4535_ER_Total, features = VariableFeatures(object = Swarbrick_CID4535_ER_Total))

# Examine and visualize PCA results a few different ways
print(Swarbrick_CID4535_ER_Total[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Swarbrick_CID4535_ER_Total, dims = 1:2, reduction = "pca")
DimPlot(Swarbrick_CID4535_ER_Total, reduction = "pca")

# Visualize PCs
ElbowPlot(Swarbrick_CID4535_ER_Total)

# Cluster cells
Swarbrick_CID4535_ER_Total <- FindNeighbors(Swarbrick_CID4535_ER_Total, dims = 1:20)
Swarbrick_CID4535_ER_Total <- FindClusters(Swarbrick_CID4535_ER_Total, resolution = 0.5)

# Run non-linear dimensional reduction
Swarbrick_CID4535_ER_Total <- RunUMAP(Swarbrick_CID4535_ER_Total, dims = 1:20)

# Visualize clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Swarbrick_CID4535_ER_Total, reduction = "umap")

# Save RDS
saveRDS(Swarbrick_CID4535_ER_Total, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Total/Swarbrick_CID4535_ER_Total.rds")
rm(Swarbrick_CID4535_ER_Total)





#########################################################
############# This is HER2 Tumor Batch # ################
#########################################################


#########################################################
#########     Swarbrick_CID3586_HER2_Total      ######### 
#########################################################
# Load in matrix
matrix_dir = "/R/R_Swarbrick/Swarbrick_Input/HER2/Swarbrick_CID3586/"
barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "features.tsv.gz")
matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")
mat <- readMM(file = matrix.path)
feature.names = read.delim(features.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
colnames(mat) = barcode.names$V1
rownames(mat) = feature.names$V1

# Read Metadata
metadata <- read.csv(file =  "/R/R_Swarbrick/Swarbrick_Input/HER2/Swarbrick_CID3586/metadata.csv", header = TRUE, row.names = 1)

# Create Seurat Object
Swarbrick_CID3586_HER2_Total <- CreateSeuratObject(mat, project="Swarbrick_CID3586_HER2_Total", min.cells = 3, min.features = 200, meta.data = metadata)
Swarbrick_CID3586_HER2_Total

# Remove matrix files
rm(mat)
rm(feature.names)
rm(barcode.names)
rm(metadata)

# The [[ opperator can add columns to object metadata. This is a great place to stash QC stats
Swarbrick_CID3586_HER2_Total[["percent.mt"]] <- PercentageFeatureSet(Swarbrick_CID3586_HER2_Total, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Swarbrick_CID3586_HER2_Total, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Swarbrick_CID3586_HER2_Total, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Swarbrick_CID3586_HER2_Total, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Swarbrick_CID3586_HER2_Total <- subset(Swarbrick_CID3586_HER2_Total, subset = nFeature_RNA > 200 & percent.mt < 20)

# Normalize data
Swarbrick_CID3586_HER2_Total <- NormalizeData(Swarbrick_CID3586_HER2_Total, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Swarbrick_CID3586_HER2_Total <- FindVariableFeatures(Swarbrick_CID3586_HER2_Total, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Swarbrick_CID3586_HER2_Total)
Swarbrick_CID3586_HER2_Total <- ScaleData(Swarbrick_CID3586_HER2_Total, features = all.genes)

# Perform linear dimensional reduction
Swarbrick_CID3586_HER2_Total <- RunPCA(Swarbrick_CID3586_HER2_Total, features = VariableFeatures(object = Swarbrick_CID3586_HER2_Total))

# Examine and visualize PCA results a few different ways
print(Swarbrick_CID3586_HER2_Total[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Swarbrick_CID3586_HER2_Total, dims = 1:2, reduction = "pca")
DimPlot(Swarbrick_CID3586_HER2_Total, reduction = "pca")

# Visualize PCs
ElbowPlot(Swarbrick_CID3586_HER2_Total)

# Cluster cells
Swarbrick_CID3586_HER2_Total <- FindNeighbors(Swarbrick_CID3586_HER2_Total, dims = 1:20)
Swarbrick_CID3586_HER2_Total <- FindClusters(Swarbrick_CID3586_HER2_Total, resolution = 0.5)

# Run non-linear dimensional reduction
Swarbrick_CID3586_HER2_Total <- RunUMAP(Swarbrick_CID3586_HER2_Total, dims = 1:20)

# Visualize Clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual Clusters
DimPlot(Swarbrick_CID3586_HER2_Total, reduction = "umap")

# Save RDS
saveRDS(Swarbrick_CID3586_HER2_Total, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Total/Swarbrick_CID3586_HER2_Total.rds")
rm(Swarbrick_CID3586_HER2_Total)






#########################################################
#########     Swarbrick_CID3838_HER2_Total      ######### 
#########################################################
# Load in matrix
matrix_dir = "/R/R_Swarbrick/Swarbrick_Input/HER2/Swarbrick_CID3838/"
barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "features.tsv.gz")
matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")
mat <- readMM(file = matrix.path)
feature.names = read.delim(features.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
colnames(mat) = barcode.names$V1
rownames(mat) = feature.names$V1

# Read Metadata
metadata <- read.csv(file =  "/R/R_Swarbrick/Swarbrick_Input/HER2/Swarbrick_CID3838/metadata.csv", header = TRUE, row.names = 1)

# Create Seurat Object
Swarbrick_CID3838_HER2_Total <- CreateSeuratObject(mat, project="Swarbrick_CID3838_HER2_Total", min.cells = 3, min.features = 200, meta.data = metadata)
Swarbrick_CID3838_HER2_Total

# Remove matrix files
rm(mat)
rm(feature.names)
rm(barcode.names)
rm(metadata)

# The [[ opperator can add columns to object metadata. This is a great place to stash QC stats
Swarbrick_CID3838_HER2_Total[["percent.mt"]] <- PercentageFeatureSet(Swarbrick_CID3838_HER2_Total, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Swarbrick_CID3838_HER2_Total, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Swarbrick_CID3838_HER2_Total, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Swarbrick_CID3838_HER2_Total, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Swarbrick_CID3838_HER2_Total <- subset(Swarbrick_CID3838_HER2_Total, subset = nFeature_RNA > 200 & percent.mt < 20)

# Normalize data
Swarbrick_CID3838_HER2_Total <- NormalizeData(Swarbrick_CID3838_HER2_Total, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Swarbrick_CID3838_HER2_Total <- FindVariableFeatures(Swarbrick_CID3838_HER2_Total, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Swarbrick_CID3838_HER2_Total)
Swarbrick_CID3838_HER2_Total <- ScaleData(Swarbrick_CID3838_HER2_Total, features = all.genes)

# Perform linear dimensional reduction
Swarbrick_CID3838_HER2_Total <- RunPCA(Swarbrick_CID3838_HER2_Total, features = VariableFeatures(object = Swarbrick_CID3838_HER2_Total))

# Examine and visualize PCA results a few different ways
print(Swarbrick_CID3838_HER2_Total[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Swarbrick_CID3838_HER2_Total, dims = 1:2, reduction = "pca")
DimPlot(Swarbrick_CID3838_HER2_Total, reduction = "pca")

# Visualize PCs
ElbowPlot(Swarbrick_CID3838_HER2_Total)

# Cluster cells
Swarbrick_CID3838_HER2_Total <- FindNeighbors(Swarbrick_CID3838_HER2_Total, dims = 1:20)
Swarbrick_CID3838_HER2_Total <- FindClusters(Swarbrick_CID3838_HER2_Total, resolution = 0.5)

# Run non-linear dimensional reduction
Swarbrick_CID3838_HER2_Total <- RunUMAP(Swarbrick_CID3838_HER2_Total, dims = 1:20)

# Visualize Clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual Clusters
DimPlot(Swarbrick_CID3838_HER2_Total, reduction = "umap")

# Save RDS
saveRDS(Swarbrick_CID3838_HER2_Total, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Total/Swarbrick_CID3838_HER2_Total.rds")
rm(Swarbrick_CID3838_HER2_Total)





#########################################################
#########     Swarbrick_CID3921_HER2_Total      ######### 
#########################################################
# Load in matrix
matrix_dir = "/R/R_Swarbrick/Swarbrick_Input/HER2/Swarbrick_CID3921/"
barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "features.tsv.gz")
matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")
mat <- readMM(file = matrix.path)
feature.names = read.delim(features.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
colnames(mat) = barcode.names$V1
rownames(mat) = feature.names$V1

# Read Metadata
metadata <- read.csv(file =  "/R/R_Swarbrick/Swarbrick_Input/HER2/Swarbrick_CID3921/metadata.csv", header = TRUE, row.names = 1)

# Create Seurat Object
Swarbrick_CID3921_HER2_Total <- CreateSeuratObject(mat, project="Swarbrick_CID3921_HER2_Total", min.cells = 3, min.features = 200, meta.data = metadata)
Swarbrick_CID3921_HER2_Total

# Remove matrix files
rm(mat)
rm(feature.names)
rm(barcode.names)
rm(metadata)

# The [[ opperator can add columns to object metadata. This is a great place to stash QC stats
Swarbrick_CID3921_HER2_Total[["percent.mt"]] <- PercentageFeatureSet(Swarbrick_CID3921_HER2_Total, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Swarbrick_CID3921_HER2_Total, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Swarbrick_CID3921_HER2_Total, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Swarbrick_CID3921_HER2_Total, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Swarbrick_CID3921_HER2_Total <- subset(Swarbrick_CID3921_HER2_Total, subset = nFeature_RNA > 200 & percent.mt < 20)

# Normalize data
Swarbrick_CID3921_HER2_Total <- NormalizeData(Swarbrick_CID3921_HER2_Total, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Swarbrick_CID3921_HER2_Total <- FindVariableFeatures(Swarbrick_CID3921_HER2_Total, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Swarbrick_CID3921_HER2_Total)
Swarbrick_CID3921_HER2_Total <- ScaleData(Swarbrick_CID3921_HER2_Total, features = all.genes)

# Perform linear dimensional reduction
Swarbrick_CID3921_HER2_Total <- RunPCA(Swarbrick_CID3921_HER2_Total, features = VariableFeatures(object = Swarbrick_CID3921_HER2_Total))

# Examine and visualize PCA results a few different ways
print(Swarbrick_CID3921_HER2_Total[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Swarbrick_CID3921_HER2_Total, dims = 1:2, reduction = "pca")
DimPlot(Swarbrick_CID3921_HER2_Total, reduction = "pca")

# Visualize PCs
ElbowPlot(Swarbrick_CID3921_HER2_Total)

# Cluster cells
Swarbrick_CID3921_HER2_Total <- FindNeighbors(Swarbrick_CID3921_HER2_Total, dims = 1:20)
Swarbrick_CID3921_HER2_Total <- FindClusters(Swarbrick_CID3921_HER2_Total, resolution = 0.5)

# Run non-linear dimensional reduction
Swarbrick_CID3921_HER2_Total <- RunUMAP(Swarbrick_CID3921_HER2_Total, dims = 1:20)

# Visualize Clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual Clusters
DimPlot(Swarbrick_CID3921_HER2_Total, reduction = "umap")

# Save RDS
saveRDS(Swarbrick_CID3921_HER2_Total, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Total/Swarbrick_CID3921_HER2_Total.rds")
rm(Swarbrick_CID3921_HER2_Total)





#########################################################
#########     Swarbrick_CID4066_HER2_Total      ######### 
#########################################################
# Load in matrix
matrix_dir = "/R/R_Swarbrick/Swarbrick_Input/HER2/Swarbrick_CID4066/"
barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "features.tsv.gz")
matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")
mat <- readMM(file = matrix.path)
feature.names = read.delim(features.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
colnames(mat) = barcode.names$V1
rownames(mat) = feature.names$V1

# Read Metadata
metadata <- read.csv(file =  "/R/R_Swarbrick/Swarbrick_Input/HER2/Swarbrick_CID4066/metadata.csv", header = TRUE, row.names = 1)

# Create Seurat Object
Swarbrick_CID4066_HER2_Total <- CreateSeuratObject(mat, project="Swarbrick_CID4066_HER2_Total", min.cells = 3, min.features = 200, meta.data = metadata)
Swarbrick_CID4066_HER2_Total

# Remove matrix files
rm(mat)
rm(feature.names)
rm(barcode.names)
rm(metadata)

# The [[ opperator can add columns to object metadata. This is a great place to stash QC stats
Swarbrick_CID4066_HER2_Total[["percent.mt"]] <- PercentageFeatureSet(Swarbrick_CID4066_HER2_Total, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Swarbrick_CID4066_HER2_Total, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Swarbrick_CID4066_HER2_Total, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Swarbrick_CID4066_HER2_Total, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Swarbrick_CID4066_HER2_Total <- subset(Swarbrick_CID4066_HER2_Total, subset = nFeature_RNA > 200 & percent.mt < 20)

# Normalize data
Swarbrick_CID4066_HER2_Total <- NormalizeData(Swarbrick_CID4066_HER2_Total, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Swarbrick_CID4066_HER2_Total <- FindVariableFeatures(Swarbrick_CID4066_HER2_Total, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Swarbrick_CID4066_HER2_Total)
Swarbrick_CID4066_HER2_Total <- ScaleData(Swarbrick_CID4066_HER2_Total, features = all.genes)

# Perform linear dimensional reduction
Swarbrick_CID4066_HER2_Total <- RunPCA(Swarbrick_CID4066_HER2_Total, features = VariableFeatures(object = Swarbrick_CID4066_HER2_Total))

# Examine and visualize PCA results a few different ways
print(Swarbrick_CID4066_HER2_Total[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Swarbrick_CID4066_HER2_Total, dims = 1:2, reduction = "pca")
DimPlot(Swarbrick_CID4066_HER2_Total, reduction = "pca")

# Visualize PCs
ElbowPlot(Swarbrick_CID4066_HER2_Total)

# Cluster cells
Swarbrick_CID4066_HER2_Total <- FindNeighbors(Swarbrick_CID4066_HER2_Total, dims = 1:20)
Swarbrick_CID4066_HER2_Total <- FindClusters(Swarbrick_CID4066_HER2_Total, resolution = 0.5)

# Run non-linear dimensional reduction
Swarbrick_CID4066_HER2_Total <- RunUMAP(Swarbrick_CID4066_HER2_Total, dims = 1:20)

# Visualize Clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual Clusters
DimPlot(Swarbrick_CID4066_HER2_Total, reduction = "umap")

# Save RDS
saveRDS(Swarbrick_CID4066_HER2_Total, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Total/Swarbrick_CID4066_HER2_Total.rds")
rm(Swarbrick_CID4066_HER2_Total)





#########################################################
#########     Swarbrick_CID45171_HER2_Total     ######### 
#########################################################
# Load in matrix
matrix_dir = "/R/R_Swarbrick/Swarbrick_Input/HER2/Swarbrick_CID45171/"
barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "features.tsv.gz")
matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")
mat <- readMM(file = matrix.path)
feature.names = read.delim(features.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
colnames(mat) = barcode.names$V1
rownames(mat) = feature.names$V1

# Read Metadata
metadata <- read.csv(file =  "/R/R_Swarbrick/Swarbrick_Input/HER2/Swarbrick_CID45171/metadata.csv", header = TRUE, row.names = 1)

# Create Seurat Object
Swarbrick_CID45171_HER2_Total <- CreateSeuratObject(mat, project="Swarbrick_CID45171_HER2_Total", min.cells = 3, min.features = 200, meta.data = metadata)
Swarbrick_CID45171_HER2_Total

# Remove matrix files
rm(mat)
rm(feature.names)
rm(barcode.names)
rm(metadata)

# The [[ opperator can add columns to object metadata. This is a great place to stash QC stats
Swarbrick_CID45171_HER2_Total[["percent.mt"]] <- PercentageFeatureSet(Swarbrick_CID45171_HER2_Total, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Swarbrick_CID45171_HER2_Total, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Swarbrick_CID45171_HER2_Total, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Swarbrick_CID45171_HER2_Total, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Swarbrick_CID45171_HER2_Total <- subset(Swarbrick_CID45171_HER2_Total, subset = nFeature_RNA > 200 & percent.mt < 20)

# Normalize data
Swarbrick_CID45171_HER2_Total <- NormalizeData(Swarbrick_CID45171_HER2_Total, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Swarbrick_CID45171_HER2_Total <- FindVariableFeatures(Swarbrick_CID45171_HER2_Total, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Swarbrick_CID45171_HER2_Total)
Swarbrick_CID45171_HER2_Total <- ScaleData(Swarbrick_CID45171_HER2_Total, features = all.genes)

# Perform linear dimensional reduction
Swarbrick_CID45171_HER2_Total <- RunPCA(Swarbrick_CID45171_HER2_Total, features = VariableFeatures(object = Swarbrick_CID45171_HER2_Total))

# Examine and visualize PCA results a few different ways
print(Swarbrick_CID45171_HER2_Total[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Swarbrick_CID45171_HER2_Total, dims = 1:2, reduction = "pca")
DimPlot(Swarbrick_CID45171_HER2_Total, reduction = "pca")

# Visualize PCs
ElbowPlot(Swarbrick_CID45171_HER2_Total)

# Cluster cells
Swarbrick_CID45171_HER2_Total <- FindNeighbors(Swarbrick_CID45171_HER2_Total, dims = 1:20)
Swarbrick_CID45171_HER2_Total <- FindClusters(Swarbrick_CID45171_HER2_Total, resolution = 0.5)

# Run non-linear dimensional reduction
Swarbrick_CID45171_HER2_Total <- RunUMAP(Swarbrick_CID45171_HER2_Total, dims = 1:20)

# Visualize Clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual Clusters
DimPlot(Swarbrick_CID45171_HER2_Total, reduction = "umap")

# Save RDS
saveRDS(Swarbrick_CID45171_HER2_Total, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Total/Swarbrick_CID45171_HER2_Total.rds")
rm(Swarbrick_CID45171_HER2_Total)







#########################################################
############# This is TNBC Tumor Batch # ################
#########################################################


#########################################################
#########    Swarbrick_CID3946_TNBC_Total       ######### 
#########################################################
# Load in matrix
matrix_dir = "/R/R_Swarbrick/Swarbrick_Input/TNBC/Swarbrick_CID3946/"
barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "features.tsv.gz")
matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")
mat <- readMM(file = matrix.path)
feature.names = read.delim(features.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
colnames(mat) = barcode.names$V1
rownames(mat) = feature.names$V1

# Read Metadata
metadata <- read.csv(file =  "/R/R_Swarbrick/Swarbrick_Input/TNBC/Swarbrick_CID3946/metadata.csv", header = TRUE, row.names = 1)

# Create Seurat Object
Swarbrick_CID3946_TNBC_Total <- CreateSeuratObject(mat, project="Swarbrick_CID3946_TNBC_Total", min.cells = 3, min.features = 200, meta.data = metadata)
Swarbrick_CID3946_TNBC_Total

# Remove matrix files
rm(mat)
rm(feature.names)
rm(barcode.names)
rm(metadata)

# The [[ opperator can add columns to object metadata. This is a great place to stash QC stats
Swarbrick_CID3946_TNBC_Total[["percent.mt"]] <- PercentageFeatureSet(Swarbrick_CID3946_TNBC_Total, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Swarbrick_CID3946_TNBC_Total, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Swarbrick_CID3946_TNBC_Total, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Swarbrick_CID3946_TNBC_Total, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Swarbrick_CID3946_TNBC_Total <- subset(Swarbrick_CID3946_TNBC_Total, subset = nFeature_RNA > 200 & percent.mt < 20)

# Normalize data
Swarbrick_CID3946_TNBC_Total <- NormalizeData(Swarbrick_CID3946_TNBC_Total, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Swarbrick_CID3946_TNBC_Total <- FindVariableFeatures(Swarbrick_CID3946_TNBC_Total, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Swarbrick_CID3946_TNBC_Total)
Swarbrick_CID3946_TNBC_Total <- ScaleData(Swarbrick_CID3946_TNBC_Total, features = all.genes)

# Perform linear dimensional reduction
Swarbrick_CID3946_TNBC_Total <- RunPCA(Swarbrick_CID3946_TNBC_Total, features = VariableFeatures(object = Swarbrick_CID3946_TNBC_Total))

# Examine and visualize PCA results a few different ways
print(Swarbrick_CID3946_TNBC_Total[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Swarbrick_CID3946_TNBC_Total, dims = 1:2, reduction = "pca")
DimPlot(Swarbrick_CID3946_TNBC_Total, reduction = "pca")

# Visualize PCs
ElbowPlot(Swarbrick_CID3946_TNBC_Total)

# Cluster cells
Swarbrick_CID3946_TNBC_Total <- FindNeighbors(Swarbrick_CID3946_TNBC_Total, dims = 1:20)
Swarbrick_CID3946_TNBC_Total <- FindClusters(Swarbrick_CID3946_TNBC_Total, resolution = 0.5)

# Run non-linear dimensional reduction
Swarbrick_CID3946_TNBC_Total <- RunUMAP(Swarbrick_CID3946_TNBC_Total, dims = 1:20)

# Visualize Clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual Clusters
DimPlot(Swarbrick_CID3946_TNBC_Total, reduction = "umap")

# Save RDS
saveRDS(Swarbrick_CID3946_TNBC_Total, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Total/Swarbrick_CID3946_TNBC_Total.rds")
rm(Swarbrick_CID3946_TNBC_Total)





#########################################################
#########     Swarbrick_CID3963_TNBC_Total      ######### 
#########################################################
# Load in matrix
matrix_dir = "/R/R_Swarbrick/Swarbrick_Input/TNBC/Swarbrick_CID3963/"
barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "features.tsv.gz")
matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")
mat <- readMM(file = matrix.path)
feature.names = read.delim(features.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
colnames(mat) = barcode.names$V1
rownames(mat) = feature.names$V1

# Read Metadata
metadata <- read.csv(file =  "/R/R_Swarbrick/Swarbrick_Input/TNBC/Swarbrick_CID3963/metadata.csv", header = TRUE, row.names = 1)

# Create Seurat Object
Swarbrick_CID3963_TNBC_Total <- CreateSeuratObject(mat, project="Swarbrick_CID3963_TNBC_Total", min.cells = 3, min.features = 200, meta.data = metadata)
Swarbrick_CID3963_TNBC_Total

# Remove matrix files
rm(mat)
rm(feature.names)
rm(barcode.names)
rm(metadata)

# The [[ opperator can add columns to object metadata. This is a great place to stash QC stats
Swarbrick_CID3963_TNBC_Total[["percent.mt"]] <- PercentageFeatureSet(Swarbrick_CID3963_TNBC_Total, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Swarbrick_CID3963_TNBC_Total, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Swarbrick_CID3963_TNBC_Total, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Swarbrick_CID3963_TNBC_Total, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Swarbrick_CID3963_TNBC_Total <- subset(Swarbrick_CID3963_TNBC_Total, subset = nFeature_RNA > 200 & percent.mt < 20)

# Normalize data
Swarbrick_CID3963_TNBC_Total <- NormalizeData(Swarbrick_CID3963_TNBC_Total, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Swarbrick_CID3963_TNBC_Total <- FindVariableFeatures(Swarbrick_CID3963_TNBC_Total, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Swarbrick_CID3963_TNBC_Total)
Swarbrick_CID3963_TNBC_Total <- ScaleData(Swarbrick_CID3963_TNBC_Total, features = all.genes)

# Perform linear dimensional reduction
Swarbrick_CID3963_TNBC_Total <- RunPCA(Swarbrick_CID3963_TNBC_Total, features = VariableFeatures(object = Swarbrick_CID3963_TNBC_Total))

# Examine and visualize PCA results a few different ways
print(Swarbrick_CID3963_TNBC_Total[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Swarbrick_CID3963_TNBC_Total, dims = 1:2, reduction = "pca")
DimPlot(Swarbrick_CID3963_TNBC_Total, reduction = "pca")

# Visualize PCs
ElbowPlot(Swarbrick_CID3963_TNBC_Total)

# Cluster cells
Swarbrick_CID3963_TNBC_Total <- FindNeighbors(Swarbrick_CID3963_TNBC_Total, dims = 1:20)
Swarbrick_CID3963_TNBC_Total <- FindClusters(Swarbrick_CID3963_TNBC_Total, resolution = 0.5)

# Run non-linear dimensional reduction
Swarbrick_CID3963_TNBC_Total <- RunUMAP(Swarbrick_CID3963_TNBC_Total, dims = 1:20)

# Visualize Clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual Clusters
DimPlot(Swarbrick_CID3963_TNBC_Total, reduction = "umap")

# Save RDS
saveRDS(Swarbrick_CID3963_TNBC_Total, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Total/Swarbrick_CID3963_TNBC_Total.rds")
rm(Swarbrick_CID3963_TNBC_Total)





#########################################################
#########     Swarbrick_CID4465_TNBC_Total      ######### 
#########################################################
# Load in matrix
matrix_dir = "/R/R_Swarbrick/Swarbrick_Input/TNBC/Swarbrick_CID4465/"
barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "features.tsv.gz")
matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")
mat <- readMM(file = matrix.path)
feature.names = read.delim(features.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
colnames(mat) = barcode.names$V1
rownames(mat) = feature.names$V1

# Read Metadata
metadata <- read.csv(file =  "/R/R_Swarbrick/Swarbrick_Input/TNBC/Swarbrick_CID4465/metadata.csv", header = TRUE, row.names = 1)

# Create Seurat Object
Swarbrick_CID4465_TNBC_Total <- CreateSeuratObject(mat, project="Swarbrick_CID4465_TNBC_Total", min.cells = 3, min.features = 200, meta.data = metadata)
Swarbrick_CID4465_TNBC_Total

# Remove matrix files
rm(mat)
rm(feature.names)
rm(barcode.names)
rm(metadata)

# The [[ opperator can add columns to object metadata. This is a great place to stash QC stats
Swarbrick_CID4465_TNBC_Total[["percent.mt"]] <- PercentageFeatureSet(Swarbrick_CID4465_TNBC_Total, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Swarbrick_CID4465_TNBC_Total, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Swarbrick_CID4465_TNBC_Total, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Swarbrick_CID4465_TNBC_Total, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Swarbrick_CID4465_TNBC_Total <- subset(Swarbrick_CID4465_TNBC_Total, subset = nFeature_RNA > 200 & percent.mt < 20)

# Normalize data
Swarbrick_CID4465_TNBC_Total <- NormalizeData(Swarbrick_CID4465_TNBC_Total, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Swarbrick_CID4465_TNBC_Total <- FindVariableFeatures(Swarbrick_CID4465_TNBC_Total, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Swarbrick_CID4465_TNBC_Total)
Swarbrick_CID4465_TNBC_Total <- ScaleData(Swarbrick_CID4465_TNBC_Total, features = all.genes)

# Perform linear dimensional reduction
Swarbrick_CID4465_TNBC_Total <- RunPCA(Swarbrick_CID4465_TNBC_Total, features = VariableFeatures(object = Swarbrick_CID4465_TNBC_Total))

# Examine and visualize PCA results a few different ways
print(Swarbrick_CID4465_TNBC_Total[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Swarbrick_CID4465_TNBC_Total, dims = 1:2, reduction = "pca")
DimPlot(Swarbrick_CID4465_TNBC_Total, reduction = "pca")

# Visualize PCs
ElbowPlot(Swarbrick_CID4465_TNBC_Total)

# Cluster cells
Swarbrick_CID4465_TNBC_Total <- FindNeighbors(Swarbrick_CID4465_TNBC_Total, dims = 1:20)
Swarbrick_CID4465_TNBC_Total <- FindClusters(Swarbrick_CID4465_TNBC_Total, resolution = 0.5)

# Run non-linear dimensional reduction
Swarbrick_CID4465_TNBC_Total <- RunUMAP(Swarbrick_CID4465_TNBC_Total, dims = 1:20)

# Visualize Clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual Clusters
DimPlot(Swarbrick_CID4465_TNBC_Total, reduction = "umap")

# Save RDS
saveRDS(Swarbrick_CID4465_TNBC_Total, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Total/Swarbrick_CID4465_TNBC_Total.rds")
rm(Swarbrick_CID4465_TNBC_Total)





#########################################################
#########     Swarbrick_CID4495_TNBC_Total      ######### 
#########################################################
# Load in matrix
matrix_dir = "/R/R_Swarbrick/Swarbrick_Input/TNBC/Swarbrick_CID4495/"
barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "features.tsv.gz")
matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")
mat <- readMM(file = matrix.path)
feature.names = read.delim(features.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
colnames(mat) = barcode.names$V1
rownames(mat) = feature.names$V1

# Read Metadata
metadata <- read.csv(file =  "/R/R_Swarbrick/Swarbrick_Input/TNBC/Swarbrick_CID4495/metadata.csv", header = TRUE, row.names = 1)

# Create Seurat Object
Swarbrick_CID4495_TNBC_Total <- CreateSeuratObject(mat, project="Swarbrick_CID4495_TNBC_Total", min.cells = 3, min.features = 200, meta.data = metadata)
Swarbrick_CID4495_TNBC_Total

# Remove matrix files
rm(mat)
rm(feature.names)
rm(barcode.names)
rm(metadata)

# The [[ opperator can add columns to object metadata. This is a great place to stash QC stats
Swarbrick_CID4495_TNBC_Total[["percent.mt"]] <- PercentageFeatureSet(Swarbrick_CID4495_TNBC_Total, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Swarbrick_CID4495_TNBC_Total, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Swarbrick_CID4495_TNBC_Total, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Swarbrick_CID4495_TNBC_Total, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Swarbrick_CID4495_TNBC_Total <- subset(Swarbrick_CID4495_TNBC_Total, subset = nFeature_RNA > 200 & percent.mt < 20)

# Normalize data
Swarbrick_CID4495_TNBC_Total <- NormalizeData(Swarbrick_CID4495_TNBC_Total, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Swarbrick_CID4495_TNBC_Total <- FindVariableFeatures(Swarbrick_CID4495_TNBC_Total, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Swarbrick_CID4495_TNBC_Total)
Swarbrick_CID4495_TNBC_Total <- ScaleData(Swarbrick_CID4495_TNBC_Total, features = all.genes)

# Perform linear dimensional reduction
Swarbrick_CID4495_TNBC_Total <- RunPCA(Swarbrick_CID4495_TNBC_Total, features = VariableFeatures(object = Swarbrick_CID4495_TNBC_Total))

# Examine and visualize PCA results a few different ways
print(Swarbrick_CID4495_TNBC_Total[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Swarbrick_CID4495_TNBC_Total, dims = 1:2, reduction = "pca")
DimPlot(Swarbrick_CID4495_TNBC_Total, reduction = "pca")

# Visualize PCs
ElbowPlot(Swarbrick_CID4495_TNBC_Total)

# Cluster cells
Swarbrick_CID4495_TNBC_Total <- FindNeighbors(Swarbrick_CID4495_TNBC_Total, dims = 1:20)
Swarbrick_CID4495_TNBC_Total <- FindClusters(Swarbrick_CID4495_TNBC_Total, resolution = 0.5)

# Run non-linear dimensional reduction
Swarbrick_CID4495_TNBC_Total <- RunUMAP(Swarbrick_CID4495_TNBC_Total, dims = 1:20)

# Visualize Clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual Clusters
DimPlot(Swarbrick_CID4495_TNBC_Total, reduction = "umap")

# Save RDS
saveRDS(Swarbrick_CID4495_TNBC_Total, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Total/Swarbrick_CID4495_TNBC_Total.rds")
rm(Swarbrick_CID4495_TNBC_Total)





#########################################################
#########     Swarbrick_CID4513_TNBC_Total      ######### 
#########################################################
# Load in matrix
matrix_dir = "/R/R_Swarbrick/Swarbrick_Input/TNBC/Swarbrick_CID4513/"
barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "features.tsv.gz")
matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")
mat <- readMM(file = matrix.path)
feature.names = read.delim(features.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
colnames(mat) = barcode.names$V1
rownames(mat) = feature.names$V1

# Read Metadata
metadata <- read.csv(file =  "/R/R_Swarbrick/Swarbrick_Input/TNBC/Swarbrick_CID4513/metadata.csv", header = TRUE, row.names = 1)

# Create Seurat Object
Swarbrick_CID4513_TNBC_Total <- CreateSeuratObject(mat, project="Swarbrick_CID4513_TNBC_Total", min.cells = 3, min.features = 200, meta.data = metadata)
Swarbrick_CID4513_TNBC_Total

# Remove matrix files
rm(mat)
rm(feature.names)
rm(barcode.names)
rm(metadata)

# The [[ opperator can add columns to object metadata. This is a great place to stash QC stats
Swarbrick_CID4513_TNBC_Total[["percent.mt"]] <- PercentageFeatureSet(Swarbrick_CID4513_TNBC_Total, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Swarbrick_CID4513_TNBC_Total, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Swarbrick_CID4513_TNBC_Total, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Swarbrick_CID4513_TNBC_Total, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Swarbrick_CID4513_TNBC_Total <- subset(Swarbrick_CID4513_TNBC_Total, subset = nFeature_RNA > 200 & percent.mt < 20)

# Normalize data
Swarbrick_CID4513_TNBC_Total <- NormalizeData(Swarbrick_CID4513_TNBC_Total, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Swarbrick_CID4513_TNBC_Total <- FindVariableFeatures(Swarbrick_CID4513_TNBC_Total, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Swarbrick_CID4513_TNBC_Total)
Swarbrick_CID4513_TNBC_Total <- ScaleData(Swarbrick_CID4513_TNBC_Total, features = all.genes)

# Perform linear dimensional reduction
Swarbrick_CID4513_TNBC_Total <- RunPCA(Swarbrick_CID4513_TNBC_Total, features = VariableFeatures(object = Swarbrick_CID4513_TNBC_Total))

# Examine and visualize PCA results a few different ways
print(Swarbrick_CID4513_TNBC_Total[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Swarbrick_CID4513_TNBC_Total, dims = 1:2, reduction = "pca")
DimPlot(Swarbrick_CID4513_TNBC_Total, reduction = "pca")

# Visualize PCs
ElbowPlot(Swarbrick_CID4513_TNBC_Total)

# Cluster cells
Swarbrick_CID4513_TNBC_Total <- FindNeighbors(Swarbrick_CID4513_TNBC_Total, dims = 1:20)
Swarbrick_CID4513_TNBC_Total <- FindClusters(Swarbrick_CID4513_TNBC_Total, resolution = 0.5)

# Run non-linear dimensional reduction
Swarbrick_CID4513_TNBC_Total <- RunUMAP(Swarbrick_CID4513_TNBC_Total, dims = 1:20)

# Visualize Clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual Clusters
DimPlot(Swarbrick_CID4513_TNBC_Total, reduction = "umap")

# Save RDS
saveRDS(Swarbrick_CID4513_TNBC_Total, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Total/Swarbrick_CID4513_TNBC_Total.rds")
rm(Swarbrick_CID4513_TNBC_Total)





#########################################################
#########     Swarbrick_CID4515_TNBC_Total      ######### 
#########################################################
# Load in matrix
matrix_dir = "/R/R_Swarbrick/Swarbrick_Input/TNBC/Swarbrick_CID4515/"
barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "features.tsv.gz")
matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")
mat <- readMM(file = matrix.path)
feature.names = read.delim(features.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
colnames(mat) = barcode.names$V1
rownames(mat) = feature.names$V1

# Read Metadata
metadata <- read.csv(file =  "/R/R_Swarbrick/Swarbrick_Input/TNBC/Swarbrick_CID4515/metadata.csv", header = TRUE, row.names = 1)

# Create Seurat Object
Swarbrick_CID4515_TNBC_Total <- CreateSeuratObject(mat, project="Swarbrick_CID4515_TNBC_Total", min.cells = 3, min.features = 200, meta.data = metadata)
Swarbrick_CID4515_TNBC_Total

# Remove matrix files
rm(mat)
rm(feature.names)
rm(barcode.names)
rm(metadata)

# The [[ opperator can add columns to object metadata. This is a great place to stash QC stats
Swarbrick_CID4515_TNBC_Total[["percent.mt"]] <- PercentageFeatureSet(Swarbrick_CID4515_TNBC_Total, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Swarbrick_CID4515_TNBC_Total, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Swarbrick_CID4515_TNBC_Total, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Swarbrick_CID4515_TNBC_Total, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Swarbrick_CID4515_TNBC_Total <- subset(Swarbrick_CID4515_TNBC_Total, subset = nFeature_RNA > 200 & percent.mt < 20)

# Normalize data
Swarbrick_CID4515_TNBC_Total <- NormalizeData(Swarbrick_CID4515_TNBC_Total, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Swarbrick_CID4515_TNBC_Total <- FindVariableFeatures(Swarbrick_CID4515_TNBC_Total, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Swarbrick_CID4515_TNBC_Total)
Swarbrick_CID4515_TNBC_Total <- ScaleData(Swarbrick_CID4515_TNBC_Total, features = all.genes)

# Perform linear dimensional reduction
Swarbrick_CID4515_TNBC_Total <- RunPCA(Swarbrick_CID4515_TNBC_Total, features = VariableFeatures(object = Swarbrick_CID4515_TNBC_Total))

# Examine and visualize PCA results a few different ways
print(Swarbrick_CID4515_TNBC_Total[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Swarbrick_CID4515_TNBC_Total, dims = 1:2, reduction = "pca")
DimPlot(Swarbrick_CID4515_TNBC_Total, reduction = "pca")

# Visualize PCs
ElbowPlot(Swarbrick_CID4515_TNBC_Total)

# Cluster cells
Swarbrick_CID4515_TNBC_Total <- FindNeighbors(Swarbrick_CID4515_TNBC_Total, dims = 1:20)
Swarbrick_CID4515_TNBC_Total <- FindClusters(Swarbrick_CID4515_TNBC_Total, resolution = 0.5)

# Run non-linear dimensional reduction
Swarbrick_CID4515_TNBC_Total <- RunUMAP(Swarbrick_CID4515_TNBC_Total, dims = 1:20)

# Visualize Clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual Clusters
DimPlot(Swarbrick_CID4515_TNBC_Total, reduction = "umap")

# Save RDS
saveRDS(Swarbrick_CID4515_TNBC_Total, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Total/Swarbrick_CID4515_TNBC_Total.rds")
rm(Swarbrick_CID4515_TNBC_Total)





#########################################################
#########     Swarbrick_CID4523_TNBC_Total      ######### 
#########################################################
# Load in matrix
matrix_dir = "/R/R_Swarbrick/Swarbrick_Input/TNBC/Swarbrick_CID4523/"
barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "features.tsv.gz")
matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")
mat <- readMM(file = matrix.path)
feature.names = read.delim(features.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
colnames(mat) = barcode.names$V1
rownames(mat) = feature.names$V1

# Read Metadata
metadata <- read.csv(file =  "/R/R_Swarbrick/Swarbrick_Input/TNBC/Swarbrick_CID4523/metadata.csv", header = TRUE, row.names = 1)

# Create Seurat Object
Swarbrick_CID4523_TNBC_Total <- CreateSeuratObject(mat, project="Swarbrick_CID4523_TNBC_Total", min.cells = 3, min.features = 200, meta.data = metadata)
Swarbrick_CID4523_TNBC_Total

# Remove matrix files
rm(mat)
rm(feature.names)
rm(barcode.names)
rm(metadata)

# The [[ opperator can add columns to object metadata. This is a great place to stash QC stats
Swarbrick_CID4523_TNBC_Total[["percent.mt"]] <- PercentageFeatureSet(Swarbrick_CID4523_TNBC_Total, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Swarbrick_CID4523_TNBC_Total, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Swarbrick_CID4523_TNBC_Total, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Swarbrick_CID4523_TNBC_Total, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Swarbrick_CID4523_TNBC_Total <- subset(Swarbrick_CID4523_TNBC_Total, subset = nFeature_RNA > 200 & percent.mt < 20)

# Normalize data
Swarbrick_CID4523_TNBC_Total <- NormalizeData(Swarbrick_CID4523_TNBC_Total, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Swarbrick_CID4523_TNBC_Total <- FindVariableFeatures(Swarbrick_CID4523_TNBC_Total, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Swarbrick_CID4523_TNBC_Total)
Swarbrick_CID4523_TNBC_Total <- ScaleData(Swarbrick_CID4523_TNBC_Total, features = all.genes)

# Perform linear dimensional reduction
Swarbrick_CID4523_TNBC_Total <- RunPCA(Swarbrick_CID4523_TNBC_Total, features = VariableFeatures(object = Swarbrick_CID4523_TNBC_Total))

# Examine and visualize PCA results a few different ways
print(Swarbrick_CID4523_TNBC_Total[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Swarbrick_CID4523_TNBC_Total, dims = 1:2, reduction = "pca")
DimPlot(Swarbrick_CID4523_TNBC_Total, reduction = "pca")

# Visualize PCs
ElbowPlot(Swarbrick_CID4523_TNBC_Total)

# Cluster cells
Swarbrick_CID4523_TNBC_Total <- FindNeighbors(Swarbrick_CID4523_TNBC_Total, dims = 1:20)
Swarbrick_CID4523_TNBC_Total <- FindClusters(Swarbrick_CID4523_TNBC_Total, resolution = 0.5)

# Run non-linear dimensional reduction
Swarbrick_CID4523_TNBC_Total <- RunUMAP(Swarbrick_CID4523_TNBC_Total, dims = 1:20)

# Visualize Clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual Clusters
DimPlot(Swarbrick_CID4523_TNBC_Total, reduction = "umap")

# Save RDS
saveRDS(Swarbrick_CID4523_TNBC_Total, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Total/Swarbrick_CID4523_TNBC_Total.rds")
rm(Swarbrick_CID4523_TNBC_Total)




#########################################################
#########     Swarbrick_CID44041_TNBC_Total     ######### 
#########################################################
# Load in matrix
matrix_dir = "/R/R_Swarbrick/Swarbrick_Input/TNBC/Swarbrick_CID44041/"
barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "features.tsv.gz")
matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")
mat <- readMM(file = matrix.path)
feature.names = read.delim(features.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
colnames(mat) = barcode.names$V1
rownames(mat) = feature.names$V1

# Read Metadata
metadata <- read.csv(file =  "/R/R_Swarbrick/Swarbrick_Input/TNBC/Swarbrick_CID44041/metadata.csv", header = TRUE, row.names = 1)

# Create Seurat Object
Swarbrick_CID44041_TNBC_Total <- CreateSeuratObject(mat, project="Swarbrick_CID44041_TNBC_Total", min.cells = 3, min.features = 200, meta.data = metadata)
Swarbrick_CID44041_TNBC_Total

# Remove matrix files
rm(mat)
rm(feature.names)
rm(barcode.names)
rm(metadata)

# The [[ opperator can add columns to object metadata. This is a great place to stash QC stats
Swarbrick_CID44041_TNBC_Total[["percent.mt"]] <- PercentageFeatureSet(Swarbrick_CID44041_TNBC_Total, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Swarbrick_CID44041_TNBC_Total, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Swarbrick_CID44041_TNBC_Total, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Swarbrick_CID44041_TNBC_Total, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Swarbrick_CID44041_TNBC_Total <- subset(Swarbrick_CID44041_TNBC_Total, subset = nFeature_RNA > 200 & percent.mt < 20)

# Normalize data
Swarbrick_CID44041_TNBC_Total <- NormalizeData(Swarbrick_CID44041_TNBC_Total, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Swarbrick_CID44041_TNBC_Total <- FindVariableFeatures(Swarbrick_CID44041_TNBC_Total, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Swarbrick_CID44041_TNBC_Total)
Swarbrick_CID44041_TNBC_Total <- ScaleData(Swarbrick_CID44041_TNBC_Total, features = all.genes)

# Perform linear dimensional reduction
Swarbrick_CID44041_TNBC_Total <- RunPCA(Swarbrick_CID44041_TNBC_Total, features = VariableFeatures(object = Swarbrick_CID44041_TNBC_Total))

# Examine and visualize PCA results a few different ways
print(Swarbrick_CID44041_TNBC_Total[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Swarbrick_CID44041_TNBC_Total, dims = 1:2, reduction = "pca")
DimPlot(Swarbrick_CID44041_TNBC_Total, reduction = "pca")

# Visualize PCs
ElbowPlot(Swarbrick_CID44041_TNBC_Total)

# Cluster cells
Swarbrick_CID44041_TNBC_Total <- FindNeighbors(Swarbrick_CID44041_TNBC_Total, dims = 1:20)
Swarbrick_CID44041_TNBC_Total <- FindClusters(Swarbrick_CID44041_TNBC_Total, resolution = 0.5)

# Run non-linear dimensional reduction
Swarbrick_CID44041_TNBC_Total <- RunUMAP(Swarbrick_CID44041_TNBC_Total, dims = 1:20)

# Visualize Clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual Clusters
DimPlot(Swarbrick_CID44041_TNBC_Total, reduction = "umap")

# Save RDS
saveRDS(Swarbrick_CID44041_TNBC_Total, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Total/Swarbrick_CID44041_TNBC_Total.rds")
rm(Swarbrick_CID44041_TNBC_Total)




#########################################################
#########     Swarbrick_CID44971_TNBC_Total     ######### 
#########################################################
# Load in matrix
matrix_dir = "/R/R_Swarbrick/Swarbrick_Input/TNBC/Swarbrick_CID44971/"
barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "features.tsv.gz")
matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")
mat <- readMM(file = matrix.path)
feature.names = read.delim(features.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
colnames(mat) = barcode.names$V1
rownames(mat) = feature.names$V1

# Read Metadata
metadata <- read.csv(file =  "/R/R_Swarbrick/Swarbrick_Input/TNBC/Swarbrick_CID44971/metadata.csv", header = TRUE, row.names = 1)

# Create Seurat Object
Swarbrick_CID44971_TNBC_Total <- CreateSeuratObject(mat, project="Swarbrick_CID44971_TNBC_Total", min.cells = 3, min.features = 200, meta.data = metadata)
Swarbrick_CID44971_TNBC_Total

# Remove matrix files
rm(mat)
rm(feature.names)
rm(barcode.names)
rm(metadata)

# The [[ opperator can add columns to object metadata. This is a great place to stash QC stats
Swarbrick_CID44971_TNBC_Total[["percent.mt"]] <- PercentageFeatureSet(Swarbrick_CID44971_TNBC_Total, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Swarbrick_CID44971_TNBC_Total, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Swarbrick_CID44971_TNBC_Total, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Swarbrick_CID44971_TNBC_Total, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Swarbrick_CID44971_TNBC_Total <- subset(Swarbrick_CID44971_TNBC_Total, subset = nFeature_RNA > 200 & percent.mt < 20)

# Normalize data
Swarbrick_CID44971_TNBC_Total <- NormalizeData(Swarbrick_CID44971_TNBC_Total, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Swarbrick_CID44971_TNBC_Total <- FindVariableFeatures(Swarbrick_CID44971_TNBC_Total, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Swarbrick_CID44971_TNBC_Total)
Swarbrick_CID44971_TNBC_Total <- ScaleData(Swarbrick_CID44971_TNBC_Total, features = all.genes)

# Perform linear dimensional reduction
Swarbrick_CID44971_TNBC_Total <- RunPCA(Swarbrick_CID44971_TNBC_Total, features = VariableFeatures(object = Swarbrick_CID44971_TNBC_Total))

# Examine and visualize PCA results a few different ways
print(Swarbrick_CID44971_TNBC_Total[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Swarbrick_CID44971_TNBC_Total, dims = 1:2, reduction = "pca")
DimPlot(Swarbrick_CID44971_TNBC_Total, reduction = "pca")

# Visualize PCs
ElbowPlot(Swarbrick_CID44971_TNBC_Total)

# Cluster cells
Swarbrick_CID44971_TNBC_Total <- FindNeighbors(Swarbrick_CID44971_TNBC_Total, dims = 1:20)
Swarbrick_CID44971_TNBC_Total <- FindClusters(Swarbrick_CID44971_TNBC_Total, resolution = 0.5)

# Run non-linear dimensional reduction
Swarbrick_CID44971_TNBC_Total <- RunUMAP(Swarbrick_CID44971_TNBC_Total, dims = 1:20)

# Visualize Clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual Clusters
DimPlot(Swarbrick_CID44971_TNBC_Total, reduction = "umap")

# Save RDS
saveRDS(Swarbrick_CID44971_TNBC_Total, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Total/Swarbrick_CID44971_TNBC_Total.rds")
rm(Swarbrick_CID44971_TNBC_Total)





#########################################################
#########     Swarbrick_CID44991_TNBC_Total     ######### 
#########################################################
# Load in matrix
matrix_dir = "/R/R_Swarbrick/Swarbrick_Input/TNBC/Swarbrick_CID44991/"
barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "features.tsv.gz")
matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")
mat <- readMM(file = matrix.path)
feature.names = read.delim(features.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
colnames(mat) = barcode.names$V1
rownames(mat) = feature.names$V1

# Read Metadata
metadata <- read.csv(file =  "/R/R_Swarbrick/Swarbrick_Input/TNBC/Swarbrick_CID44991/metadata.csv", header = TRUE, row.names = 1)

# Create Seurat Object
Swarbrick_CID44991_TNBC_Total <- CreateSeuratObject(mat, project="Swarbrick_CID44991_TNBC_Total", min.cells = 3, min.features = 200, meta.data = metadata)
Swarbrick_CID44991_TNBC_Total

# Remove matrix files
rm(mat)
rm(feature.names)
rm(barcode.names)
rm(metadata)

# The [[ opperator can add columns to object metadata. This is a great place to stash QC stats
Swarbrick_CID44991_TNBC_Total[["percent.mt"]] <- PercentageFeatureSet(Swarbrick_CID44991_TNBC_Total, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Swarbrick_CID44991_TNBC_Total, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Swarbrick_CID44991_TNBC_Total, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Swarbrick_CID44991_TNBC_Total, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
Swarbrick_CID44991_TNBC_Total <- subset(Swarbrick_CID44991_TNBC_Total, subset = nFeature_RNA > 200 & percent.mt < 20)

# Normalize data
Swarbrick_CID44991_TNBC_Total <- NormalizeData(Swarbrick_CID44991_TNBC_Total, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
Swarbrick_CID44991_TNBC_Total <- FindVariableFeatures(Swarbrick_CID44991_TNBC_Total, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(Swarbrick_CID44991_TNBC_Total)
Swarbrick_CID44991_TNBC_Total <- ScaleData(Swarbrick_CID44991_TNBC_Total, features = all.genes)

# Perform linear dimensional reduction
Swarbrick_CID44991_TNBC_Total <- RunPCA(Swarbrick_CID44991_TNBC_Total, features = VariableFeatures(object = Swarbrick_CID44991_TNBC_Total))

# Examine and visualize PCA results a few different ways
print(Swarbrick_CID44991_TNBC_Total[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Swarbrick_CID44991_TNBC_Total, dims = 1:2, reduction = "pca")
DimPlot(Swarbrick_CID44991_TNBC_Total, reduction = "pca")

# Visualize PCs
ElbowPlot(Swarbrick_CID44991_TNBC_Total)

# Cluster cells
Swarbrick_CID44991_TNBC_Total <- FindNeighbors(Swarbrick_CID44991_TNBC_Total, dims = 1:20)
Swarbrick_CID44991_TNBC_Total <- FindClusters(Swarbrick_CID44991_TNBC_Total, resolution = 0.5)

# Run non-linear dimensional reduction
Swarbrick_CID44991_TNBC_Total <- RunUMAP(Swarbrick_CID44991_TNBC_Total, dims = 1:20)

# Visualize Clusters
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual Clusters
DimPlot(Swarbrick_CID44991_TNBC_Total, reduction = "umap")

# Save RDS
saveRDS(Swarbrick_CID44991_TNBC_Total, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Total/Swarbrick_CID44991_TNBC_Total.rds")
rm(Swarbrick_CID44991_TNBC_Total)


#####################################################################
########### Step 2: Run DoubletFinder and Subset Singlets ########### 
#####################################################################

# Load Libraries
library(Seurat)
library(DoubletFinder)

# Assuming 5% doublet rate per sample based on average number of sequenced cells
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-

# Patient Samples
################################# ER
# Swarbrick_CID3941_ER_Total
# Swarbrick_CID3948_ER_Total
# Swarbrick_CID4040_ER_Total
# Swarbrick_CID4067_ER_Total
# Swarbrick_CID4290A_ER_Total
# Swarbrick_CID4398_ER_Total
# Swarbrick_CID4461_ER_Total
# Swarbrick_CID4463_ER_Total
# Swarbrick_CID4471_ER_Total
# Swarbrick_CID4530N_ER_Total
# Swarbrick_CID4535_ER_Total
#################################  HER2
# Swarbrick_CID3586_HER2_Total
# Swarbrick_CID3838_HER2_Total
# Swarbrick_CID3921_HER2_Total
# Swarbrick_CID4066_HER2_Total
# Swarbrick_CID45171_HER2_Total
#################################  TNBC
# Swarbrick_CID3946_TNBC_Total
# Swarbrick_CID3963_TNBC_Total
# Swarbrick_CID4465_TNBC_Total
# Swarbrick_CID4495_TNBC_Total
# Swarbrick_CID4513_TNBC_Total
# Swarbrick_CID4515_TNBC_Total
# Swarbrick_CID4523_TNBC_Total
# Swarbrick_CID44041_TNBC_Total
# Swarbrick_CID44971_TNBC_Total
# Swarbrick_CID44991_TNBC_Total

################################################################################################
########################                   ER-pos Tumors                ########################
################################################################################################

################################################################################################
########################              Swarbrick_CID3941_ER_Total             ###################
################################################################################################

# Load Data
Swarbrick_CID3941_ER_Total <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Total/Swarbrick_CID3941_ER_Total.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Swarbrick_CID3941_ER_Total <- paramSweep(Swarbrick_CID3941_ER_Total, PCs = 1:10, sct = TRUE)
sweep.stats.Swarbrick_CID3941_ER_Total <- summarizeSweep(sweep.res.Swarbrick_CID3941_ER_Total, GT = FALSE)
bcmvn_Swarbrick_CID3941_ER_Total <- find.pK(sweep.stats.Swarbrick_CID3941_ER_Total)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Swarbrick_CID3941_ER_Total[which.max(bcmvn_Swarbrick_CID3941_ER_Total$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Swarbrick_CID3941_ER_Total@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 5% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Swarbrick dataset = 5k per sample, will use doublet rate for 6,000 cells to err on the side of caution
nExp_poi <- round(0.05*nrow(Swarbrick_CID3941_ER_Total@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Swarbrick_CID3941_ER_Total <- doubletFinder(Swarbrick_CID3941_ER_Total, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Swarbrick_CID3941_ER_Total_Quant <- (Swarbrick_CID3941_ER_Total@meta.data$DF.classification == "Singlet")
Swarbrick_CID3941_ER_Total_Quant_Singlets <- length(Swarbrick_CID3941_ER_Total_Quant[Swarbrick_CID3941_ER_Total_Quant== TRUE])
Swarbrick_CID3941_ER_Total_Quant_Doublets <- length(Swarbrick_CID3941_ER_Total_Quant[Swarbrick_CID3941_ER_Total_Quant== FALSE])
Swarbrick_CID3941_ER_Total_Quant_Doublets_Percent <- Swarbrick_CID3941_ER_Total_Quant_Doublets / (Swarbrick_CID3941_ER_Total_Quant_Doublets + Swarbrick_CID3941_ER_Total_Quant_Singlets)* 100
Swarbrick_CID3941_ER_Total_Quant <- as.data.frame(c(Swarbrick_CID3941_ER_Total_Quant_Singlets, Swarbrick_CID3941_ER_Total_Quant_Doublets, Swarbrick_CID3941_ER_Total_Quant_Doublets_Percent))
colnames(Swarbrick_CID3941_ER_Total_Quant) <- c("Swarbrick_CID3941_ER_Total_Quant")
rownames(Swarbrick_CID3941_ER_Total_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Swarbrick_CID3941_ER_Total_Quant, file = "/R/R_Swarbrick/Swarbrick_Doublet_Tables/Swarbrick_CID3941_ER_Total_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Swarbrick_CID3941_ER_Total_Singlets <- subset(Swarbrick_CID3941_ER_Total, cells=rownames(Swarbrick_CID3941_ER_Total@meta.data)[which(Swarbrick_CID3941_ER_Total@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Swarbrick_CID3941_ER_Total_Singlets, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Total_Singlets/Swarbrick_CID3941_ER_Total_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Swarbrick_CID3941_ER_Total)
rm(Swarbrick_CID3941_ER_Total_Singlets)
rm(Swarbrick_CID3941_ER_Total_Quant)
rm(Swarbrick_CID3941_ER_Total_Quant_Singlets)
rm(Swarbrick_CID3941_ER_Total_Quant_Doublets)
rm(Swarbrick_CID3941_ER_Total_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Swarbrick_CID3941_ER_Total)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Swarbrick_CID3941_ER_Total)
rm(sweep.stats.Swarbrick_CID3941_ER_Total)
gc()




################################################################################################
########################              Swarbrick_CID3948_ER_Total             ###################
################################################################################################

# Load Data
Swarbrick_CID3948_ER_Total <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Total/Swarbrick_CID3948_ER_Total.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Swarbrick_CID3948_ER_Total <- paramSweep(Swarbrick_CID3948_ER_Total, PCs = 1:10, sct = TRUE)
sweep.stats.Swarbrick_CID3948_ER_Total <- summarizeSweep(sweep.res.Swarbrick_CID3948_ER_Total, GT = FALSE)
bcmvn_Swarbrick_CID3948_ER_Total <- find.pK(sweep.stats.Swarbrick_CID3948_ER_Total)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Swarbrick_CID3948_ER_Total[which.max(bcmvn_Swarbrick_CID3948_ER_Total$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Swarbrick_CID3948_ER_Total@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 5% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Swarbrick dataset = 5k per sample, will use doublet rate for 6,000 cells to err on the side of caution
nExp_poi <- round(0.05*nrow(Swarbrick_CID3948_ER_Total@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Swarbrick_CID3948_ER_Total <- doubletFinder(Swarbrick_CID3948_ER_Total, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Swarbrick_CID3948_ER_Total_Quant <- (Swarbrick_CID3948_ER_Total@meta.data$DF.classification == "Singlet")
Swarbrick_CID3948_ER_Total_Quant_Singlets <- length(Swarbrick_CID3948_ER_Total_Quant[Swarbrick_CID3948_ER_Total_Quant== TRUE])
Swarbrick_CID3948_ER_Total_Quant_Doublets <- length(Swarbrick_CID3948_ER_Total_Quant[Swarbrick_CID3948_ER_Total_Quant== FALSE])
Swarbrick_CID3948_ER_Total_Quant_Doublets_Percent <- Swarbrick_CID3948_ER_Total_Quant_Doublets / (Swarbrick_CID3948_ER_Total_Quant_Doublets + Swarbrick_CID3948_ER_Total_Quant_Singlets) * 100
Swarbrick_CID3948_ER_Total_Quant <- as.data.frame(c(Swarbrick_CID3948_ER_Total_Quant_Singlets, Swarbrick_CID3948_ER_Total_Quant_Doublets, Swarbrick_CID3948_ER_Total_Quant_Doublets_Percent))
colnames(Swarbrick_CID3948_ER_Total_Quant) <- c("Swarbrick_CID3948_ER_Total_Quant")
rownames(Swarbrick_CID3948_ER_Total_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Swarbrick_CID3948_ER_Total_Quant, file = "/R/R_Swarbrick/Swarbrick_Doublet_Tables/Swarbrick_CID3948_ER_Total_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Swarbrick_CID3948_ER_Total_Singlets <- subset(Swarbrick_CID3948_ER_Total, cells=rownames(Swarbrick_CID3948_ER_Total@meta.data)[which(Swarbrick_CID3948_ER_Total@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Swarbrick_CID3948_ER_Total_Singlets, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Total_Singlets/Swarbrick_CID3948_ER_Total_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Swarbrick_CID3948_ER_Total)
rm(Swarbrick_CID3948_ER_Total_Singlets)
rm(Swarbrick_CID3948_ER_Total_Quant)
rm(Swarbrick_CID3948_ER_Total_Quant_Singlets)
rm(Swarbrick_CID3948_ER_Total_Quant_Doublets)
rm(Swarbrick_CID3948_ER_Total_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Swarbrick_CID3948_ER_Total)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Swarbrick_CID3948_ER_Total)
rm(sweep.stats.Swarbrick_CID3948_ER_Total)
gc()




################################################################################################
########################              Swarbrick_CID4040_ER_Total             ###################
################################################################################################

# Load Data
Swarbrick_CID4040_ER_Total <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Total/Swarbrick_CID4040_ER_Total.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Swarbrick_CID4040_ER_Total <- paramSweep(Swarbrick_CID4040_ER_Total, PCs = 1:10, sct = TRUE)
sweep.stats.Swarbrick_CID4040_ER_Total <- summarizeSweep(sweep.res.Swarbrick_CID4040_ER_Total, GT = FALSE)
bcmvn_Swarbrick_CID4040_ER_Total <- find.pK(sweep.stats.Swarbrick_CID4040_ER_Total)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Swarbrick_CID4040_ER_Total[which.max(bcmvn_Swarbrick_CID4040_ER_Total$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Swarbrick_CID4040_ER_Total@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 5% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Swarbrick dataset = 5k per sample, will use doublet rate for 6,000 cells to err on the side of caution
nExp_poi <- round(0.05*nrow(Swarbrick_CID4040_ER_Total@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Swarbrick_CID4040_ER_Total <- doubletFinder(Swarbrick_CID4040_ER_Total, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Swarbrick_CID4040_ER_Total_Quant <- (Swarbrick_CID4040_ER_Total@meta.data$DF.classification == "Singlet")
Swarbrick_CID4040_ER_Total_Quant_Singlets <- length(Swarbrick_CID4040_ER_Total_Quant[Swarbrick_CID4040_ER_Total_Quant== TRUE])
Swarbrick_CID4040_ER_Total_Quant_Doublets <- length(Swarbrick_CID4040_ER_Total_Quant[Swarbrick_CID4040_ER_Total_Quant== FALSE])
Swarbrick_CID4040_ER_Total_Quant_Doublets_Percent <- Swarbrick_CID4040_ER_Total_Quant_Doublets / (Swarbrick_CID4040_ER_Total_Quant_Doublets + Swarbrick_CID4040_ER_Total_Quant_Singlets)* 100
Swarbrick_CID4040_ER_Total_Quant <- as.data.frame(c(Swarbrick_CID4040_ER_Total_Quant_Singlets, Swarbrick_CID4040_ER_Total_Quant_Doublets, Swarbrick_CID4040_ER_Total_Quant_Doublets_Percent))
colnames(Swarbrick_CID4040_ER_Total_Quant) <- c("Swarbrick_CID4040_ER_Total_Quant")
rownames(Swarbrick_CID4040_ER_Total_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Swarbrick_CID4040_ER_Total_Quant, file = "/R/R_Swarbrick/Swarbrick_Doublet_Tables/Swarbrick_CID4040_ER_Total_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Swarbrick_CID4040_ER_Total_Singlets <- subset(Swarbrick_CID4040_ER_Total, cells=rownames(Swarbrick_CID4040_ER_Total@meta.data)[which(Swarbrick_CID4040_ER_Total@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Swarbrick_CID4040_ER_Total_Singlets, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Total_Singlets/Swarbrick_CID4040_ER_Total_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Swarbrick_CID4040_ER_Total)
rm(Swarbrick_CID4040_ER_Total_Singlets)
rm(Swarbrick_CID4040_ER_Total_Quant)
rm(Swarbrick_CID4040_ER_Total_Quant_Singlets)
rm(Swarbrick_CID4040_ER_Total_Quant_Doublets)
rm(Swarbrick_CID4040_ER_Total_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Swarbrick_CID4040_ER_Total)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Swarbrick_CID4040_ER_Total)
rm(sweep.stats.Swarbrick_CID4040_ER_Total)
gc()




################################################################################################
########################              Swarbrick_CID4067_ER_Total             ###################
################################################################################################

# Load Data
Swarbrick_CID4067_ER_Total <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Total/Swarbrick_CID4067_ER_Total.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Swarbrick_CID4067_ER_Total <- paramSweep(Swarbrick_CID4067_ER_Total, PCs = 1:10, sct = TRUE)
sweep.stats.Swarbrick_CID4067_ER_Total <- summarizeSweep(sweep.res.Swarbrick_CID4067_ER_Total, GT = FALSE)
bcmvn_Swarbrick_CID4067_ER_Total <- find.pK(sweep.stats.Swarbrick_CID4067_ER_Total)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Swarbrick_CID4067_ER_Total[which.max(bcmvn_Swarbrick_CID4067_ER_Total$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Swarbrick_CID4067_ER_Total@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 5% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Swarbrick dataset = 5k per sample, will use doublet rate for 6,000 cells to err on the side of caution
nExp_poi <- round(0.05*nrow(Swarbrick_CID4067_ER_Total@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Swarbrick_CID4067_ER_Total <- doubletFinder(Swarbrick_CID4067_ER_Total, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Swarbrick_CID4067_ER_Total_Quant <- (Swarbrick_CID4067_ER_Total@meta.data$DF.classification == "Singlet")
Swarbrick_CID4067_ER_Total_Quant_Singlets <- length(Swarbrick_CID4067_ER_Total_Quant[Swarbrick_CID4067_ER_Total_Quant== TRUE])
Swarbrick_CID4067_ER_Total_Quant_Doublets <- length(Swarbrick_CID4067_ER_Total_Quant[Swarbrick_CID4067_ER_Total_Quant== FALSE])
Swarbrick_CID4067_ER_Total_Quant_Doublets_Percent <- Swarbrick_CID4067_ER_Total_Quant_Doublets / (Swarbrick_CID4067_ER_Total_Quant_Doublets + Swarbrick_CID4067_ER_Total_Quant_Singlets) * 100
Swarbrick_CID4067_ER_Total_Quant <- as.data.frame(c(Swarbrick_CID4067_ER_Total_Quant_Singlets, Swarbrick_CID4067_ER_Total_Quant_Doublets, Swarbrick_CID4067_ER_Total_Quant_Doublets_Percent))
colnames(Swarbrick_CID4067_ER_Total_Quant) <- c("Swarbrick_CID4067_ER_Total_Quant")
rownames(Swarbrick_CID4067_ER_Total_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Swarbrick_CID4067_ER_Total_Quant, file = "/R/R_Swarbrick/Swarbrick_Doublet_Tables/Swarbrick_CID4067_ER_Total_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Swarbrick_CID4067_ER_Total_Singlets <- subset(Swarbrick_CID4067_ER_Total, cells=rownames(Swarbrick_CID4067_ER_Total@meta.data)[which(Swarbrick_CID4067_ER_Total@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Swarbrick_CID4067_ER_Total_Singlets, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Total_Singlets/Swarbrick_CID4067_ER_Total_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Swarbrick_CID4067_ER_Total)
rm(Swarbrick_CID4067_ER_Total_Singlets)
rm(Swarbrick_CID4067_ER_Total_Quant)
rm(Swarbrick_CID4067_ER_Total_Quant_Singlets)
rm(Swarbrick_CID4067_ER_Total_Quant_Doublets)
rm(Swarbrick_CID4067_ER_Total_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Swarbrick_CID4067_ER_Total)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Swarbrick_CID4067_ER_Total)
rm(sweep.stats.Swarbrick_CID4067_ER_Total)
gc()




################################################################################################
########################              Swarbrick_CID4290A_ER_Total             ##################
################################################################################################

# Load Data
Swarbrick_CID4290A_ER_Total <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Total/Swarbrick_CID4290A_ER_Total.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Swarbrick_CID4290A_ER_Total <- paramSweep(Swarbrick_CID4290A_ER_Total, PCs = 1:10, sct = TRUE)
sweep.stats.Swarbrick_CID4290A_ER_Total <- summarizeSweep(sweep.res.Swarbrick_CID4290A_ER_Total, GT = FALSE)
bcmvn_Swarbrick_CID4290A_ER_Total <- find.pK(sweep.stats.Swarbrick_CID4290A_ER_Total)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Swarbrick_CID4290A_ER_Total[which.max(bcmvn_Swarbrick_CID4290A_ER_Total$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Swarbrick_CID4290A_ER_Total@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 5% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Swarbrick dataset = 5k per sample, will use doublet rate for 6,000 cells to err on the side of caution
nExp_poi <- round(0.05*nrow(Swarbrick_CID4290A_ER_Total@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Swarbrick_CID4290A_ER_Total <- doubletFinder(Swarbrick_CID4290A_ER_Total, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Swarbrick_CID4290A_ER_Total_Quant <- (Swarbrick_CID4290A_ER_Total@meta.data$DF.classification == "Singlet")
Swarbrick_CID4290A_ER_Total_Quant_Singlets <- length(Swarbrick_CID4290A_ER_Total_Quant[Swarbrick_CID4290A_ER_Total_Quant== TRUE])
Swarbrick_CID4290A_ER_Total_Quant_Doublets <- length(Swarbrick_CID4290A_ER_Total_Quant[Swarbrick_CID4290A_ER_Total_Quant== FALSE])
Swarbrick_CID4290A_ER_Total_Quant_Doublets_Percent <- Swarbrick_CID4290A_ER_Total_Quant_Doublets / (Swarbrick_CID4290A_ER_Total_Quant_Doublets + Swarbrick_CID4290A_ER_Total_Quant_Singlets) * 100
Swarbrick_CID4290A_ER_Total_Quant <- as.data.frame(c(Swarbrick_CID4290A_ER_Total_Quant_Singlets, Swarbrick_CID4290A_ER_Total_Quant_Doublets, Swarbrick_CID4290A_ER_Total_Quant_Doublets_Percent))
colnames(Swarbrick_CID4290A_ER_Total_Quant) <- c("Swarbrick_CID4290A_ER_Total_Quant")
rownames(Swarbrick_CID4290A_ER_Total_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Swarbrick_CID4290A_ER_Total_Quant, file = "/R/R_Swarbrick/Swarbrick_Doublet_Tables/Swarbrick_CID4290A_ER_Total_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Swarbrick_CID4290A_ER_Total_Singlets <- subset(Swarbrick_CID4290A_ER_Total, cells=rownames(Swarbrick_CID4290A_ER_Total@meta.data)[which(Swarbrick_CID4290A_ER_Total@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Swarbrick_CID4290A_ER_Total_Singlets, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Total_Singlets/Swarbrick_CID4290A_ER_Total_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Swarbrick_CID4290A_ER_Total)
rm(Swarbrick_CID4290A_ER_Total_Singlets)
rm(Swarbrick_CID4290A_ER_Total_Quant)
rm(Swarbrick_CID4290A_ER_Total_Quant_Singlets)
rm(Swarbrick_CID4290A_ER_Total_Quant_Doublets)
rm(Swarbrick_CID4290A_ER_Total_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Swarbrick_CID4290A_ER_Total)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Swarbrick_CID4290A_ER_Total)
rm(sweep.stats.Swarbrick_CID4290A_ER_Total)
gc()




################################################################################################
########################              Swarbrick_CID4398_ER_Total             ###################
################################################################################################

# Load Data
Swarbrick_CID4398_ER_Total <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Total/Swarbrick_CID4398_ER_Total.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Swarbrick_CID4398_ER_Total <- paramSweep(Swarbrick_CID4398_ER_Total, PCs = 1:10, sct = TRUE)
sweep.stats.Swarbrick_CID4398_ER_Total <- summarizeSweep(sweep.res.Swarbrick_CID4398_ER_Total, GT = FALSE)
bcmvn_Swarbrick_CID4398_ER_Total <- find.pK(sweep.stats.Swarbrick_CID4398_ER_Total)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Swarbrick_CID4398_ER_Total[which.max(bcmvn_Swarbrick_CID4398_ER_Total$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Swarbrick_CID4398_ER_Total@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 5% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Swarbrick dataset = 5k per sample, will use doublet rate for 6,000 cells to err on the side of caution
nExp_poi <- round(0.05*nrow(Swarbrick_CID4398_ER_Total@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Swarbrick_CID4398_ER_Total <- doubletFinder(Swarbrick_CID4398_ER_Total, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Swarbrick_CID4398_ER_Total_Quant <- (Swarbrick_CID4398_ER_Total@meta.data$DF.classification == "Singlet")
Swarbrick_CID4398_ER_Total_Quant_Singlets <- length(Swarbrick_CID4398_ER_Total_Quant[Swarbrick_CID4398_ER_Total_Quant== TRUE])
Swarbrick_CID4398_ER_Total_Quant_Doublets <- length(Swarbrick_CID4398_ER_Total_Quant[Swarbrick_CID4398_ER_Total_Quant== FALSE])
Swarbrick_CID4398_ER_Total_Quant_Doublets_Percent <- Swarbrick_CID4398_ER_Total_Quant_Doublets / (Swarbrick_CID4398_ER_Total_Quant_Doublets + Swarbrick_CID4398_ER_Total_Quant_Singlets) * 100
Swarbrick_CID4398_ER_Total_Quant <- as.data.frame(c(Swarbrick_CID4398_ER_Total_Quant_Singlets, Swarbrick_CID4398_ER_Total_Quant_Doublets, Swarbrick_CID4398_ER_Total_Quant_Doublets_Percent))
colnames(Swarbrick_CID4398_ER_Total_Quant) <- c("Swarbrick_CID4398_ER_Total_Quant")
rownames(Swarbrick_CID4398_ER_Total_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Swarbrick_CID4398_ER_Total_Quant, file = "/R/R_Swarbrick/Swarbrick_Doublet_Tables/Swarbrick_CID4398_ER_Total_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Swarbrick_CID4398_ER_Total_Singlets <- subset(Swarbrick_CID4398_ER_Total, cells=rownames(Swarbrick_CID4398_ER_Total@meta.data)[which(Swarbrick_CID4398_ER_Total@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Swarbrick_CID4398_ER_Total_Singlets, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Total_Singlets/Swarbrick_CID4398_ER_Total_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Swarbrick_CID4398_ER_Total)
rm(Swarbrick_CID4398_ER_Total_Singlets)
rm(Swarbrick_CID4398_ER_Total_Quant)
rm(Swarbrick_CID4398_ER_Total_Quant_Singlets)
rm(Swarbrick_CID4398_ER_Total_Quant_Doublets)
rm(Swarbrick_CID4398_ER_Total_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Swarbrick_CID4398_ER_Total)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Swarbrick_CID4398_ER_Total)
rm(sweep.stats.Swarbrick_CID4398_ER_Total)
gc()




################################################################################################
########################              Swarbrick_CID4461_ER_Total             ###################
################################################################################################

# Load Data
Swarbrick_CID4461_ER_Total <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Total/Swarbrick_CID4461_ER_Total.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Swarbrick_CID4461_ER_Total <- paramSweep(Swarbrick_CID4461_ER_Total, PCs = 1:10, sct = TRUE)
sweep.stats.Swarbrick_CID4461_ER_Total <- summarizeSweep(sweep.res.Swarbrick_CID4461_ER_Total, GT = FALSE)
bcmvn_Swarbrick_CID4461_ER_Total <- find.pK(sweep.stats.Swarbrick_CID4461_ER_Total)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Swarbrick_CID4461_ER_Total[which.max(bcmvn_Swarbrick_CID4461_ER_Total$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Swarbrick_CID4461_ER_Total@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 5% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Swarbrick dataset = 5k per sample, will use doublet rate for 6,000 cells to err on the side of caution
nExp_poi <- round(0.05*nrow(Swarbrick_CID4461_ER_Total@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Swarbrick_CID4461_ER_Total <- doubletFinder(Swarbrick_CID4461_ER_Total, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Swarbrick_CID4461_ER_Total_Quant <- (Swarbrick_CID4461_ER_Total@meta.data$DF.classification == "Singlet")
Swarbrick_CID4461_ER_Total_Quant_Singlets <- length(Swarbrick_CID4461_ER_Total_Quant[Swarbrick_CID4461_ER_Total_Quant== TRUE])
Swarbrick_CID4461_ER_Total_Quant_Doublets <- length(Swarbrick_CID4461_ER_Total_Quant[Swarbrick_CID4461_ER_Total_Quant== FALSE])
Swarbrick_CID4461_ER_Total_Quant_Doublets_Percent <- Swarbrick_CID4461_ER_Total_Quant_Doublets / (Swarbrick_CID4461_ER_Total_Quant_Doublets + Swarbrick_CID4461_ER_Total_Quant_Singlets) * 100
Swarbrick_CID4461_ER_Total_Quant <- as.data.frame(c(Swarbrick_CID4461_ER_Total_Quant_Singlets, Swarbrick_CID4461_ER_Total_Quant_Doublets, Swarbrick_CID4461_ER_Total_Quant_Doublets_Percent))
colnames(Swarbrick_CID4461_ER_Total_Quant) <- c("Swarbrick_CID4461_ER_Total_Quant")
rownames(Swarbrick_CID4461_ER_Total_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Swarbrick_CID4461_ER_Total_Quant, file = "/R/R_Swarbrick/Swarbrick_Doublet_Tables/Swarbrick_CID4461_ER_Total_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Swarbrick_CID4461_ER_Total_Singlets <- subset(Swarbrick_CID4461_ER_Total, cells=rownames(Swarbrick_CID4461_ER_Total@meta.data)[which(Swarbrick_CID4461_ER_Total@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Swarbrick_CID4461_ER_Total_Singlets, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Total_Singlets/Swarbrick_CID4461_ER_Total_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Swarbrick_CID4461_ER_Total)
rm(Swarbrick_CID4461_ER_Total_Singlets)
rm(Swarbrick_CID4461_ER_Total_Quant)
rm(Swarbrick_CID4461_ER_Total_Quant_Singlets)
rm(Swarbrick_CID4461_ER_Total_Quant_Doublets)
rm(Swarbrick_CID4461_ER_Total_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Swarbrick_CID4461_ER_Total)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Swarbrick_CID4461_ER_Total)
rm(sweep.stats.Swarbrick_CID4461_ER_Total)
gc()




################################################################################################
########################              Swarbrick_CID4463_ER_Total             ###################
################################################################################################

# Load Data
Swarbrick_CID4463_ER_Total <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Total/Swarbrick_CID4463_ER_Total.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Swarbrick_CID4463_ER_Total <- paramSweep(Swarbrick_CID4463_ER_Total, PCs = 1:10, sct = TRUE)
sweep.stats.Swarbrick_CID4463_ER_Total <- summarizeSweep(sweep.res.Swarbrick_CID4463_ER_Total, GT = FALSE)
bcmvn_Swarbrick_CID4463_ER_Total <- find.pK(sweep.stats.Swarbrick_CID4463_ER_Total)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Swarbrick_CID4463_ER_Total[which.max(bcmvn_Swarbrick_CID4463_ER_Total$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Swarbrick_CID4463_ER_Total@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 5% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Swarbrick dataset = 5k per sample, will use doublet rate for 6,000 cells to err on the side of caution
nExp_poi <- round(0.05*nrow(Swarbrick_CID4463_ER_Total@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Swarbrick_CID4463_ER_Total <- doubletFinder(Swarbrick_CID4463_ER_Total, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Swarbrick_CID4463_ER_Total_Quant <- (Swarbrick_CID4463_ER_Total@meta.data$DF.classification == "Singlet")
Swarbrick_CID4463_ER_Total_Quant_Singlets <- length(Swarbrick_CID4463_ER_Total_Quant[Swarbrick_CID4463_ER_Total_Quant== TRUE])
Swarbrick_CID4463_ER_Total_Quant_Doublets <- length(Swarbrick_CID4463_ER_Total_Quant[Swarbrick_CID4463_ER_Total_Quant== FALSE])
Swarbrick_CID4463_ER_Total_Quant_Doublets_Percent <- Swarbrick_CID4463_ER_Total_Quant_Doublets / (Swarbrick_CID4463_ER_Total_Quant_Doublets + Swarbrick_CID4463_ER_Total_Quant_Singlets) * 100
Swarbrick_CID4463_ER_Total_Quant <- as.data.frame(c(Swarbrick_CID4463_ER_Total_Quant_Singlets, Swarbrick_CID4463_ER_Total_Quant_Doublets, Swarbrick_CID4463_ER_Total_Quant_Doublets_Percent))
colnames(Swarbrick_CID4463_ER_Total_Quant) <- c("Swarbrick_CID4463_ER_Total_Quant")
rownames(Swarbrick_CID4463_ER_Total_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Swarbrick_CID4463_ER_Total_Quant, file = "/R/R_Swarbrick/Swarbrick_Doublet_Tables/Swarbrick_CID4463_ER_Total_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Swarbrick_CID4463_ER_Total_Singlets <- subset(Swarbrick_CID4463_ER_Total, cells=rownames(Swarbrick_CID4463_ER_Total@meta.data)[which(Swarbrick_CID4463_ER_Total@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Swarbrick_CID4463_ER_Total_Singlets, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Total_Singlets/Swarbrick_CID4463_ER_Total_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Swarbrick_CID4463_ER_Total)
rm(Swarbrick_CID4463_ER_Total_Singlets)
rm(Swarbrick_CID4463_ER_Total_Quant)
rm(Swarbrick_CID4463_ER_Total_Quant_Singlets)
rm(Swarbrick_CID4463_ER_Total_Quant_Doublets)
rm(Swarbrick_CID4463_ER_Total_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Swarbrick_CID4463_ER_Total)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Swarbrick_CID4463_ER_Total)
rm(sweep.stats.Swarbrick_CID4463_ER_Total)
gc()




################################################################################################
########################              Swarbrick_CID4471_ER_Total             ###################
################################################################################################

# Load Data
Swarbrick_CID4471_ER_Total <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Total/Swarbrick_CID4471_ER_Total.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Swarbrick_CID4471_ER_Total <- paramSweep(Swarbrick_CID4471_ER_Total, PCs = 1:10, sct = TRUE)
sweep.stats.Swarbrick_CID4471_ER_Total <- summarizeSweep(sweep.res.Swarbrick_CID4471_ER_Total, GT = FALSE)
bcmvn_Swarbrick_CID4471_ER_Total <- find.pK(sweep.stats.Swarbrick_CID4471_ER_Total)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Swarbrick_CID4471_ER_Total[which.max(bcmvn_Swarbrick_CID4471_ER_Total$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Swarbrick_CID4471_ER_Total@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 5% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Swarbrick dataset = 5k per sample, will use doublet rate for 6,000 cells to err on the side of caution
nExp_poi <- round(0.05*nrow(Swarbrick_CID4471_ER_Total@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Swarbrick_CID4471_ER_Total <- doubletFinder(Swarbrick_CID4471_ER_Total, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Swarbrick_CID4471_ER_Total_Quant <- (Swarbrick_CID4471_ER_Total@meta.data$DF.classification == "Singlet")
Swarbrick_CID4471_ER_Total_Quant_Singlets <- length(Swarbrick_CID4471_ER_Total_Quant[Swarbrick_CID4471_ER_Total_Quant== TRUE])
Swarbrick_CID4471_ER_Total_Quant_Doublets <- length(Swarbrick_CID4471_ER_Total_Quant[Swarbrick_CID4471_ER_Total_Quant== FALSE])
Swarbrick_CID4471_ER_Total_Quant_Doublets_Percent <- Swarbrick_CID4471_ER_Total_Quant_Doublets / (Swarbrick_CID4471_ER_Total_Quant_Doublets + Swarbrick_CID4471_ER_Total_Quant_Singlets) * 100
Swarbrick_CID4471_ER_Total_Quant <- as.data.frame(c(Swarbrick_CID4471_ER_Total_Quant_Singlets, Swarbrick_CID4471_ER_Total_Quant_Doublets, Swarbrick_CID4471_ER_Total_Quant_Doublets_Percent))
colnames(Swarbrick_CID4471_ER_Total_Quant) <- c("Swarbrick_CID4471_ER_Total_Quant")
rownames(Swarbrick_CID4471_ER_Total_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Swarbrick_CID4471_ER_Total_Quant, file = "/R/R_Swarbrick/Swarbrick_Doublet_Tables/Swarbrick_CID4471_ER_Total_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Swarbrick_CID4471_ER_Total_Singlets <- subset(Swarbrick_CID4471_ER_Total, cells=rownames(Swarbrick_CID4471_ER_Total@meta.data)[which(Swarbrick_CID4471_ER_Total@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Swarbrick_CID4471_ER_Total_Singlets, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Total_Singlets/Swarbrick_CID4471_ER_Total_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Swarbrick_CID4471_ER_Total)
rm(Swarbrick_CID4471_ER_Total_Singlets)
rm(Swarbrick_CID4471_ER_Total_Quant)
rm(Swarbrick_CID4471_ER_Total_Quant_Singlets)
rm(Swarbrick_CID4471_ER_Total_Quant_Doublets)
rm(Swarbrick_CID4471_ER_Total_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Swarbrick_CID4471_ER_Total)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Swarbrick_CID4471_ER_Total)
rm(sweep.stats.Swarbrick_CID4471_ER_Total)
gc()




################################################################################################
########################              Swarbrick_CID4530N_ER_Total             ##################
################################################################################################

# Load Data
Swarbrick_CID4530N_ER_Total <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Total/Swarbrick_CID4530N_ER_Total.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Swarbrick_CID4530N_ER_Total <- paramSweep(Swarbrick_CID4530N_ER_Total, PCs = 1:10, sct = TRUE)
sweep.stats.Swarbrick_CID4530N_ER_Total <- summarizeSweep(sweep.res.Swarbrick_CID4530N_ER_Total, GT = FALSE)
bcmvn_Swarbrick_CID4530N_ER_Total <- find.pK(sweep.stats.Swarbrick_CID4530N_ER_Total)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Swarbrick_CID4530N_ER_Total[which.max(bcmvn_Swarbrick_CID4530N_ER_Total$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Swarbrick_CID4530N_ER_Total@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 5% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Swarbrick dataset = 5k per sample, will use doublet rate for 6,000 cells to err on the side of caution
nExp_poi <- round(0.05*nrow(Swarbrick_CID4530N_ER_Total@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Swarbrick_CID4530N_ER_Total <- doubletFinder(Swarbrick_CID4530N_ER_Total, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Swarbrick_CID4530N_ER_Total_Quant <- (Swarbrick_CID4530N_ER_Total@meta.data$DF.classification == "Singlet")
Swarbrick_CID4530N_ER_Total_Quant_Singlets <- length(Swarbrick_CID4530N_ER_Total_Quant[Swarbrick_CID4530N_ER_Total_Quant== TRUE])
Swarbrick_CID4530N_ER_Total_Quant_Doublets <- length(Swarbrick_CID4530N_ER_Total_Quant[Swarbrick_CID4530N_ER_Total_Quant== FALSE])
Swarbrick_CID4530N_ER_Total_Quant_Doublets_Percent <- Swarbrick_CID4530N_ER_Total_Quant_Doublets / (Swarbrick_CID4530N_ER_Total_Quant_Doublets + Swarbrick_CID4530N_ER_Total_Quant_Singlets) * 100
Swarbrick_CID4530N_ER_Total_Quant <- as.data.frame(c(Swarbrick_CID4530N_ER_Total_Quant_Singlets, Swarbrick_CID4530N_ER_Total_Quant_Doublets, Swarbrick_CID4530N_ER_Total_Quant_Doublets_Percent))
colnames(Swarbrick_CID4530N_ER_Total_Quant) <- c("Swarbrick_CID4530N_ER_Total_Quant")
rownames(Swarbrick_CID4530N_ER_Total_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Swarbrick_CID4530N_ER_Total_Quant, file = "/R/R_Swarbrick/Swarbrick_Doublet_Tables/Swarbrick_CID4530N_ER_Total_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Swarbrick_CID4530N_ER_Total_Singlets <- subset(Swarbrick_CID4530N_ER_Total, cells=rownames(Swarbrick_CID4530N_ER_Total@meta.data)[which(Swarbrick_CID4530N_ER_Total@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Swarbrick_CID4530N_ER_Total_Singlets, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Total_Singlets/Swarbrick_CID4530N_ER_Total_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Swarbrick_CID4530N_ER_Total)
rm(Swarbrick_CID4530N_ER_Total_Singlets)
rm(Swarbrick_CID4530N_ER_Total_Quant)
rm(Swarbrick_CID4530N_ER_Total_Quant_Singlets)
rm(Swarbrick_CID4530N_ER_Total_Quant_Doublets)
rm(Swarbrick_CID4530N_ER_Total_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Swarbrick_CID4530N_ER_Total)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Swarbrick_CID4530N_ER_Total)
rm(sweep.stats.Swarbrick_CID4530N_ER_Total)
gc()




################################################################################################
########################              Swarbrick_CID4535_ER_Total             ###################
################################################################################################

# Load Data
Swarbrick_CID4535_ER_Total <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Total/Swarbrick_CID4535_ER_Total.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Swarbrick_CID4535_ER_Total <- paramSweep(Swarbrick_CID4535_ER_Total, PCs = 1:10, sct = TRUE)
sweep.stats.Swarbrick_CID4535_ER_Total <- summarizeSweep(sweep.res.Swarbrick_CID4535_ER_Total, GT = FALSE)
bcmvn_Swarbrick_CID4535_ER_Total <- find.pK(sweep.stats.Swarbrick_CID4535_ER_Total)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Swarbrick_CID4535_ER_Total[which.max(bcmvn_Swarbrick_CID4535_ER_Total$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Swarbrick_CID4535_ER_Total@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 5% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Swarbrick dataset = 5k per sample, will use doublet rate for 6,000 cells to err on the side of caution
nExp_poi <- round(0.05*nrow(Swarbrick_CID4535_ER_Total@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Swarbrick_CID4535_ER_Total <- doubletFinder(Swarbrick_CID4535_ER_Total, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Swarbrick_CID4535_ER_Total_Quant <- (Swarbrick_CID4535_ER_Total@meta.data$DF.classification == "Singlet")
Swarbrick_CID4535_ER_Total_Quant_Singlets <- length(Swarbrick_CID4535_ER_Total_Quant[Swarbrick_CID4535_ER_Total_Quant== TRUE])
Swarbrick_CID4535_ER_Total_Quant_Doublets <- length(Swarbrick_CID4535_ER_Total_Quant[Swarbrick_CID4535_ER_Total_Quant== FALSE])
Swarbrick_CID4535_ER_Total_Quant_Doublets_Percent <- Swarbrick_CID4535_ER_Total_Quant_Doublets / (Swarbrick_CID4535_ER_Total_Quant_Doublets + Swarbrick_CID4535_ER_Total_Quant_Singlets) * 100
Swarbrick_CID4535_ER_Total_Quant <- as.data.frame(c(Swarbrick_CID4535_ER_Total_Quant_Singlets, Swarbrick_CID4535_ER_Total_Quant_Doublets, Swarbrick_CID4535_ER_Total_Quant_Doublets_Percent))
colnames(Swarbrick_CID4535_ER_Total_Quant) <- c("Swarbrick_CID4535_ER_Total_Quant")
rownames(Swarbrick_CID4535_ER_Total_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Swarbrick_CID4535_ER_Total_Quant, file = "/R/R_Swarbrick/Swarbrick_Doublet_Tables/Swarbrick_CID4535_ER_Total_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Swarbrick_CID4535_ER_Total_Singlets <- subset(Swarbrick_CID4535_ER_Total, cells=rownames(Swarbrick_CID4535_ER_Total@meta.data)[which(Swarbrick_CID4535_ER_Total@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Swarbrick_CID4535_ER_Total_Singlets, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Total_Singlets/Swarbrick_CID4535_ER_Total_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Swarbrick_CID4535_ER_Total)
rm(Swarbrick_CID4535_ER_Total_Singlets)
rm(Swarbrick_CID4535_ER_Total_Quant)
rm(Swarbrick_CID4535_ER_Total_Quant_Singlets)
rm(Swarbrick_CID4535_ER_Total_Quant_Doublets)
rm(Swarbrick_CID4535_ER_Total_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Swarbrick_CID4535_ER_Total)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Swarbrick_CID4535_ER_Total)
rm(sweep.stats.Swarbrick_CID4535_ER_Total)
gc()


################################################################################################
########################                HER2-pos Tumors                #########################
################################################################################################

################################################################################################
########################              Swarbrick_CID3586_HER2_Total             #################
################################################################################################

# Load Data
Swarbrick_CID3586_HER2_Total <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Total/Swarbrick_CID3586_HER2_Total.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Swarbrick_CID3586_HER2_Total <- paramSweep(Swarbrick_CID3586_HER2_Total, PCs = 1:10, sct = TRUE)
sweep.stats.Swarbrick_CID3586_HER2_Total <- summarizeSweep(sweep.res.Swarbrick_CID3586_HER2_Total, GT = FALSE)
bcmvn_Swarbrick_CID3586_HER2_Total <- find.pK(sweep.stats.Swarbrick_CID3586_HER2_Total)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Swarbrick_CID3586_HER2_Total[which.max(bcmvn_Swarbrick_CID3586_HER2_Total$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Swarbrick_CID3586_HER2_Total@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 5% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Swarbrick dataset = 5k per sample, will use doublet rate for 6,000 cells to err on the side of caution
nExp_poi <- round(0.05*nrow(Swarbrick_CID3586_HER2_Total@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Swarbrick_CID3586_HER2_Total <- doubletFinder(Swarbrick_CID3586_HER2_Total, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Swarbrick_CID3586_HER2_Total_Quant <- (Swarbrick_CID3586_HER2_Total@meta.data$DF.classification == "Singlet")
Swarbrick_CID3586_HER2_Total_Quant_Singlets <- length(Swarbrick_CID3586_HER2_Total_Quant[Swarbrick_CID3586_HER2_Total_Quant== TRUE])
Swarbrick_CID3586_HER2_Total_Quant_Doublets <- length(Swarbrick_CID3586_HER2_Total_Quant[Swarbrick_CID3586_HER2_Total_Quant== FALSE])
Swarbrick_CID3586_HER2_Total_Quant_Doublets_Percent <- Swarbrick_CID3586_HER2_Total_Quant_Doublets / (Swarbrick_CID3586_HER2_Total_Quant_Doublets + Swarbrick_CID3586_HER2_Total_Quant_Singlets) * 100
Swarbrick_CID3586_HER2_Total_Quant <- as.data.frame(c(Swarbrick_CID3586_HER2_Total_Quant_Singlets, Swarbrick_CID3586_HER2_Total_Quant_Doublets, Swarbrick_CID3586_HER2_Total_Quant_Doublets_Percent))
colnames(Swarbrick_CID3586_HER2_Total_Quant) <- c("Swarbrick_CID3586_HER2_Total_Quant")
rownames(Swarbrick_CID3586_HER2_Total_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Swarbrick_CID3586_HER2_Total_Quant, file = "/R/R_Swarbrick/Swarbrick_Doublet_Tables/Swarbrick_CID3586_HER2_Total_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Swarbrick_CID3586_HER2_Total_Singlets <- subset(Swarbrick_CID3586_HER2_Total, cells=rownames(Swarbrick_CID3586_HER2_Total@meta.data)[which(Swarbrick_CID3586_HER2_Total@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Swarbrick_CID3586_HER2_Total_Singlets, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Total_Singlets/Swarbrick_CID3586_HER2_Total_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Swarbrick_CID3586_HER2_Total)
rm(Swarbrick_CID3586_HER2_Total_Singlets)
rm(Swarbrick_CID3586_HER2_Total_Quant)
rm(Swarbrick_CID3586_HER2_Total_Quant_Singlets)
rm(Swarbrick_CID3586_HER2_Total_Quant_Doublets)
rm(Swarbrick_CID3586_HER2_Total_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Swarbrick_CID3586_HER2_Total)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Swarbrick_CID3586_HER2_Total)
rm(sweep.stats.Swarbrick_CID3586_HER2_Total)
gc()



################################################################################################
########################              Swarbrick_CID3838_HER2_Total             #################
################################################################################################

# Load Data
Swarbrick_CID3838_HER2_Total <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Total/Swarbrick_CID3838_HER2_Total.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Swarbrick_CID3838_HER2_Total <- paramSweep(Swarbrick_CID3838_HER2_Total, PCs = 1:10, sct = TRUE)
sweep.stats.Swarbrick_CID3838_HER2_Total <- summarizeSweep(sweep.res.Swarbrick_CID3838_HER2_Total, GT = FALSE)
bcmvn_Swarbrick_CID3838_HER2_Total <- find.pK(sweep.stats.Swarbrick_CID3838_HER2_Total)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Swarbrick_CID3838_HER2_Total[which.max(bcmvn_Swarbrick_CID3838_HER2_Total$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Swarbrick_CID3838_HER2_Total@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 5% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Swarbrick dataset = 5k per sample, will use doublet rate for 6,000 cells to err on the side of caution
nExp_poi <- round(0.05*nrow(Swarbrick_CID3838_HER2_Total@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Swarbrick_CID3838_HER2_Total <- doubletFinder(Swarbrick_CID3838_HER2_Total, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Swarbrick_CID3838_HER2_Total_Quant <- (Swarbrick_CID3838_HER2_Total@meta.data$DF.classification == "Singlet")
Swarbrick_CID3838_HER2_Total_Quant_Singlets <- length(Swarbrick_CID3838_HER2_Total_Quant[Swarbrick_CID3838_HER2_Total_Quant== TRUE])
Swarbrick_CID3838_HER2_Total_Quant_Doublets <- length(Swarbrick_CID3838_HER2_Total_Quant[Swarbrick_CID3838_HER2_Total_Quant== FALSE])
Swarbrick_CID3838_HER2_Total_Quant_Doublets_Percent <- Swarbrick_CID3838_HER2_Total_Quant_Doublets / (Swarbrick_CID3838_HER2_Total_Quant_Doublets + Swarbrick_CID3838_HER2_Total_Quant_Singlets) * 100
Swarbrick_CID3838_HER2_Total_Quant <- as.data.frame(c(Swarbrick_CID3838_HER2_Total_Quant_Singlets, Swarbrick_CID3838_HER2_Total_Quant_Doublets, Swarbrick_CID3838_HER2_Total_Quant_Doublets_Percent))
colnames(Swarbrick_CID3838_HER2_Total_Quant) <- c("Swarbrick_CID3838_HER2_Total_Quant")
rownames(Swarbrick_CID3838_HER2_Total_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Swarbrick_CID3838_HER2_Total_Quant, file = "/R/R_Swarbrick/Swarbrick_Doublet_Tables/Swarbrick_CID3838_HER2_Total_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Swarbrick_CID3838_HER2_Total_Singlets <- subset(Swarbrick_CID3838_HER2_Total, cells=rownames(Swarbrick_CID3838_HER2_Total@meta.data)[which(Swarbrick_CID3838_HER2_Total@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Swarbrick_CID3838_HER2_Total_Singlets, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Total_Singlets/Swarbrick_CID3838_HER2_Total_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Swarbrick_CID3838_HER2_Total)
rm(Swarbrick_CID3838_HER2_Total_Singlets)
rm(Swarbrick_CID3838_HER2_Total_Quant)
rm(Swarbrick_CID3838_HER2_Total_Quant_Singlets)
rm(Swarbrick_CID3838_HER2_Total_Quant_Doublets)
rm(Swarbrick_CID3838_HER2_Total_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Swarbrick_CID3838_HER2_Total)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Swarbrick_CID3838_HER2_Total)
rm(sweep.stats.Swarbrick_CID3838_HER2_Total)
gc()



################################################################################################
########################              Swarbrick_CID3921_HER2_Total             #################
################################################################################################

# Load Data
Swarbrick_CID3921_HER2_Total <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Total/Swarbrick_CID3921_HER2_Total.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Swarbrick_CID3921_HER2_Total <- paramSweep(Swarbrick_CID3921_HER2_Total, PCs = 1:10, sct = TRUE)
sweep.stats.Swarbrick_CID3921_HER2_Total <- summarizeSweep(sweep.res.Swarbrick_CID3921_HER2_Total, GT = FALSE)
bcmvn_Swarbrick_CID3921_HER2_Total <- find.pK(sweep.stats.Swarbrick_CID3921_HER2_Total)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Swarbrick_CID3921_HER2_Total[which.max(bcmvn_Swarbrick_CID3921_HER2_Total$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Swarbrick_CID3921_HER2_Total@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 5% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Swarbrick dataset = 5k per sample, will use doublet rate for 6,000 cells to err on the side of caution
nExp_poi <- round(0.05*nrow(Swarbrick_CID3921_HER2_Total@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Swarbrick_CID3921_HER2_Total <- doubletFinder(Swarbrick_CID3921_HER2_Total, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Swarbrick_CID3921_HER2_Total_Quant <- (Swarbrick_CID3921_HER2_Total@meta.data$DF.classification == "Singlet")
Swarbrick_CID3921_HER2_Total_Quant_Singlets <- length(Swarbrick_CID3921_HER2_Total_Quant[Swarbrick_CID3921_HER2_Total_Quant== TRUE])
Swarbrick_CID3921_HER2_Total_Quant_Doublets <- length(Swarbrick_CID3921_HER2_Total_Quant[Swarbrick_CID3921_HER2_Total_Quant== FALSE])
Swarbrick_CID3921_HER2_Total_Quant_Doublets_Percent <- Swarbrick_CID3921_HER2_Total_Quant_Doublets / (Swarbrick_CID3921_HER2_Total_Quant_Doublets + Swarbrick_CID3921_HER2_Total_Quant_Singlets) * 100
Swarbrick_CID3921_HER2_Total_Quant <- as.data.frame(c(Swarbrick_CID3921_HER2_Total_Quant_Singlets, Swarbrick_CID3921_HER2_Total_Quant_Doublets, Swarbrick_CID3921_HER2_Total_Quant_Doublets_Percent))
colnames(Swarbrick_CID3921_HER2_Total_Quant) <- c("Swarbrick_CID3921_HER2_Total_Quant")
rownames(Swarbrick_CID3921_HER2_Total_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Swarbrick_CID3921_HER2_Total_Quant, file = "/R/R_Swarbrick/Swarbrick_Doublet_Tables/Swarbrick_CID3921_HER2_Total_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Swarbrick_CID3921_HER2_Total_Singlets <- subset(Swarbrick_CID3921_HER2_Total, cells=rownames(Swarbrick_CID3921_HER2_Total@meta.data)[which(Swarbrick_CID3921_HER2_Total@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Swarbrick_CID3921_HER2_Total_Singlets, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Total_Singlets/Swarbrick_CID3921_HER2_Total_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Swarbrick_CID3921_HER2_Total)
rm(Swarbrick_CID3921_HER2_Total_Singlets)
rm(Swarbrick_CID3921_HER2_Total_Quant)
rm(Swarbrick_CID3921_HER2_Total_Quant_Singlets)
rm(Swarbrick_CID3921_HER2_Total_Quant_Doublets)
rm(Swarbrick_CID3921_HER2_Total_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Swarbrick_CID3921_HER2_Total)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Swarbrick_CID3921_HER2_Total)
rm(sweep.stats.Swarbrick_CID3921_HER2_Total)
gc()



################################################################################################
########################              Swarbrick_CID4066_HER2_Total             #################
################################################################################################

# Load Data
Swarbrick_CID4066_HER2_Total <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Total/Swarbrick_CID4066_HER2_Total.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Swarbrick_CID4066_HER2_Total <- paramSweep(Swarbrick_CID4066_HER2_Total, PCs = 1:10, sct = TRUE)
sweep.stats.Swarbrick_CID4066_HER2_Total <- summarizeSweep(sweep.res.Swarbrick_CID4066_HER2_Total, GT = FALSE)
bcmvn_Swarbrick_CID4066_HER2_Total <- find.pK(sweep.stats.Swarbrick_CID4066_HER2_Total)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Swarbrick_CID4066_HER2_Total[which.max(bcmvn_Swarbrick_CID4066_HER2_Total$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Swarbrick_CID4066_HER2_Total@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 5% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Swarbrick dataset = 5k per sample, will use doublet rate for 6,000 cells to err on the side of caution
nExp_poi <- round(0.05*nrow(Swarbrick_CID4066_HER2_Total@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Swarbrick_CID4066_HER2_Total <- doubletFinder(Swarbrick_CID4066_HER2_Total, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Swarbrick_CID4066_HER2_Total_Quant <- (Swarbrick_CID4066_HER2_Total@meta.data$DF.classification == "Singlet")
Swarbrick_CID4066_HER2_Total_Quant_Singlets <- length(Swarbrick_CID4066_HER2_Total_Quant[Swarbrick_CID4066_HER2_Total_Quant== TRUE])
Swarbrick_CID4066_HER2_Total_Quant_Doublets <- length(Swarbrick_CID4066_HER2_Total_Quant[Swarbrick_CID4066_HER2_Total_Quant== FALSE])
Swarbrick_CID4066_HER2_Total_Quant_Doublets_Percent <- Swarbrick_CID4066_HER2_Total_Quant_Doublets / (Swarbrick_CID4066_HER2_Total_Quant_Doublets + Swarbrick_CID4066_HER2_Total_Quant_Singlets) * 100
Swarbrick_CID4066_HER2_Total_Quant <- as.data.frame(c(Swarbrick_CID4066_HER2_Total_Quant_Singlets, Swarbrick_CID4066_HER2_Total_Quant_Doublets, Swarbrick_CID4066_HER2_Total_Quant_Doublets_Percent))
colnames(Swarbrick_CID4066_HER2_Total_Quant) <- c("Swarbrick_CID4066_HER2_Total_Quant")
rownames(Swarbrick_CID4066_HER2_Total_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Swarbrick_CID4066_HER2_Total_Quant, file = "/R/R_Swarbrick/Swarbrick_Doublet_Tables/Swarbrick_CID4066_HER2_Total_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Swarbrick_CID4066_HER2_Total_Singlets <- subset(Swarbrick_CID4066_HER2_Total, cells=rownames(Swarbrick_CID4066_HER2_Total@meta.data)[which(Swarbrick_CID4066_HER2_Total@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Swarbrick_CID4066_HER2_Total_Singlets, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Total_Singlets/Swarbrick_CID4066_HER2_Total_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Swarbrick_CID4066_HER2_Total)
rm(Swarbrick_CID4066_HER2_Total_Singlets)
rm(Swarbrick_CID4066_HER2_Total_Quant)
rm(Swarbrick_CID4066_HER2_Total_Quant_Singlets)
rm(Swarbrick_CID4066_HER2_Total_Quant_Doublets)
rm(Swarbrick_CID4066_HER2_Total_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Swarbrick_CID4066_HER2_Total)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Swarbrick_CID4066_HER2_Total)
rm(sweep.stats.Swarbrick_CID4066_HER2_Total)
gc()



################################################################################################
########################              Swarbrick_CID45171_HER2_Total             ################
################################################################################################

# Load Data
Swarbrick_CID45171_HER2_Total <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Total/Swarbrick_CID45171_HER2_Total.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Swarbrick_CID45171_HER2_Total <- paramSweep(Swarbrick_CID45171_HER2_Total, PCs = 1:10, sct = TRUE)
sweep.stats.Swarbrick_CID45171_HER2_Total <- summarizeSweep(sweep.res.Swarbrick_CID45171_HER2_Total, GT = FALSE)
bcmvn_Swarbrick_CID45171_HER2_Total <- find.pK(sweep.stats.Swarbrick_CID45171_HER2_Total)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Swarbrick_CID45171_HER2_Total[which.max(bcmvn_Swarbrick_CID45171_HER2_Total$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Swarbrick_CID45171_HER2_Total@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 5% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Swarbrick dataset = 5k per sample, will use doublet rate for 6,000 cells to err on the side of caution
nExp_poi <- round(0.05*nrow(Swarbrick_CID45171_HER2_Total@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Swarbrick_CID45171_HER2_Total <- doubletFinder(Swarbrick_CID45171_HER2_Total, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Swarbrick_CID45171_HER2_Total_Quant <- (Swarbrick_CID45171_HER2_Total@meta.data$DF.classification == "Singlet")
Swarbrick_CID45171_HER2_Total_Quant_Singlets <- length(Swarbrick_CID45171_HER2_Total_Quant[Swarbrick_CID45171_HER2_Total_Quant== TRUE])
Swarbrick_CID45171_HER2_Total_Quant_Doublets <- length(Swarbrick_CID45171_HER2_Total_Quant[Swarbrick_CID45171_HER2_Total_Quant== FALSE])
Swarbrick_CID45171_HER2_Total_Quant_Doublets_Percent <- Swarbrick_CID45171_HER2_Total_Quant_Doublets / (Swarbrick_CID45171_HER2_Total_Quant_Doublets + Swarbrick_CID45171_HER2_Total_Quant_Singlets) * 100
Swarbrick_CID45171_HER2_Total_Quant <- as.data.frame(c(Swarbrick_CID45171_HER2_Total_Quant_Singlets, Swarbrick_CID45171_HER2_Total_Quant_Doublets, Swarbrick_CID45171_HER2_Total_Quant_Doublets_Percent))
colnames(Swarbrick_CID45171_HER2_Total_Quant) <- c("Swarbrick_CID45171_HER2_Total_Quant")
rownames(Swarbrick_CID45171_HER2_Total_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Swarbrick_CID45171_HER2_Total_Quant, file = "/R/R_Swarbrick/Swarbrick_Doublet_Tables/Swarbrick_CID45171_HER2_Total_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Swarbrick_CID45171_HER2_Total_Singlets <- subset(Swarbrick_CID45171_HER2_Total, cells=rownames(Swarbrick_CID45171_HER2_Total@meta.data)[which(Swarbrick_CID45171_HER2_Total@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Swarbrick_CID45171_HER2_Total_Singlets, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Total_Singlets/Swarbrick_CID45171_HER2_Total_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Swarbrick_CID45171_HER2_Total)
rm(Swarbrick_CID45171_HER2_Total_Singlets)
rm(Swarbrick_CID45171_HER2_Total_Quant)
rm(Swarbrick_CID45171_HER2_Total_Quant_Singlets)
rm(Swarbrick_CID45171_HER2_Total_Quant_Doublets)
rm(Swarbrick_CID45171_HER2_Total_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Swarbrick_CID45171_HER2_Total)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Swarbrick_CID45171_HER2_Total)
rm(sweep.stats.Swarbrick_CID45171_HER2_Total)
gc()



################################################################################################
########################                   TNBC Tumors                  ########################
################################################################################################

################################################################################################
########################              Swarbrick_CID3946_TNBC_Total             #################
################################################################################################

# Load Data
Swarbrick_CID3946_TNBC_Total <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Total/Swarbrick_CID3946_TNBC_Total.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Swarbrick_CID3946_TNBC_Total <- paramSweep(Swarbrick_CID3946_TNBC_Total, PCs = 1:10, sct = TRUE)
sweep.stats.Swarbrick_CID3946_TNBC_Total <- summarizeSweep(sweep.res.Swarbrick_CID3946_TNBC_Total, GT = FALSE)
bcmvn_Swarbrick_CID3946_TNBC_Total <- find.pK(sweep.stats.Swarbrick_CID3946_TNBC_Total)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Swarbrick_CID3946_TNBC_Total[which.max(bcmvn_Swarbrick_CID3946_TNBC_Total$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Swarbrick_CID3946_TNBC_Total@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 5% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Swarbrick dataset = 5k per sample, will use doublet rate for 6,000 cells to err on the side of caution
nExp_poi <- round(0.05*nrow(Swarbrick_CID3946_TNBC_Total@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Swarbrick_CID3946_TNBC_Total <- doubletFinder(Swarbrick_CID3946_TNBC_Total, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Swarbrick_CID3946_TNBC_Total_Quant <- (Swarbrick_CID3946_TNBC_Total@meta.data$DF.classification == "Singlet")
Swarbrick_CID3946_TNBC_Total_Quant_Singlets <- length(Swarbrick_CID3946_TNBC_Total_Quant[Swarbrick_CID3946_TNBC_Total_Quant== TRUE])
Swarbrick_CID3946_TNBC_Total_Quant_Doublets <- length(Swarbrick_CID3946_TNBC_Total_Quant[Swarbrick_CID3946_TNBC_Total_Quant== FALSE])
Swarbrick_CID3946_TNBC_Total_Quant_Doublets_Percent <- Swarbrick_CID3946_TNBC_Total_Quant_Doublets / (Swarbrick_CID3946_TNBC_Total_Quant_Doublets + Swarbrick_CID3946_TNBC_Total_Quant_Singlets) * 100
Swarbrick_CID3946_TNBC_Total_Quant <- as.data.frame(c(Swarbrick_CID3946_TNBC_Total_Quant_Singlets, Swarbrick_CID3946_TNBC_Total_Quant_Doublets, Swarbrick_CID3946_TNBC_Total_Quant_Doublets_Percent))
colnames(Swarbrick_CID3946_TNBC_Total_Quant) <- c("Swarbrick_CID3946_TNBC_Total_Quant")
rownames(Swarbrick_CID3946_TNBC_Total_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Swarbrick_CID3946_TNBC_Total_Quant, file = "/R/R_Swarbrick/Swarbrick_Doublet_Tables/Swarbrick_CID3946_TNBC_Total_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Swarbrick_CID3946_TNBC_Total_Singlets <- subset(Swarbrick_CID3946_TNBC_Total, cells=rownames(Swarbrick_CID3946_TNBC_Total@meta.data)[which(Swarbrick_CID3946_TNBC_Total@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Swarbrick_CID3946_TNBC_Total_Singlets, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Total_Singlets/Swarbrick_CID3946_TNBC_Total_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Swarbrick_CID3946_TNBC_Total)
rm(Swarbrick_CID3946_TNBC_Total_Singlets)
rm(Swarbrick_CID3946_TNBC_Total_Quant)
rm(Swarbrick_CID3946_TNBC_Total_Quant_Singlets)
rm(Swarbrick_CID3946_TNBC_Total_Quant_Doublets)
rm(Swarbrick_CID3946_TNBC_Total_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Swarbrick_CID3946_TNBC_Total)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Swarbrick_CID3946_TNBC_Total)
rm(sweep.stats.Swarbrick_CID3946_TNBC_Total)
gc()



################################################################################################
########################              Swarbrick_CID3963_TNBC_Total             #################
################################################################################################

# Load Data
Swarbrick_CID3963_TNBC_Total <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Total/Swarbrick_CID3963_TNBC_Total.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Swarbrick_CID3963_TNBC_Total <- paramSweep(Swarbrick_CID3963_TNBC_Total, PCs = 1:10, sct = TRUE)
sweep.stats.Swarbrick_CID3963_TNBC_Total <- summarizeSweep(sweep.res.Swarbrick_CID3963_TNBC_Total, GT = FALSE)
bcmvn_Swarbrick_CID3963_TNBC_Total <- find.pK(sweep.stats.Swarbrick_CID3963_TNBC_Total)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Swarbrick_CID3963_TNBC_Total[which.max(bcmvn_Swarbrick_CID3963_TNBC_Total$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Swarbrick_CID3963_TNBC_Total@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 5% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Swarbrick dataset = 5k per sample, will use doublet rate for 6,000 cells to err on the side of caution
nExp_poi <- round(0.05*nrow(Swarbrick_CID3963_TNBC_Total@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Swarbrick_CID3963_TNBC_Total <- doubletFinder(Swarbrick_CID3963_TNBC_Total, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Swarbrick_CID3963_TNBC_Total_Quant <- (Swarbrick_CID3963_TNBC_Total@meta.data$DF.classification == "Singlet")
Swarbrick_CID3963_TNBC_Total_Quant_Singlets <- length(Swarbrick_CID3963_TNBC_Total_Quant[Swarbrick_CID3963_TNBC_Total_Quant== TRUE])
Swarbrick_CID3963_TNBC_Total_Quant_Doublets <- length(Swarbrick_CID3963_TNBC_Total_Quant[Swarbrick_CID3963_TNBC_Total_Quant== FALSE])
Swarbrick_CID3963_TNBC_Total_Quant_Doublets_Percent <- Swarbrick_CID3963_TNBC_Total_Quant_Doublets / (Swarbrick_CID3963_TNBC_Total_Quant_Doublets + Swarbrick_CID3963_TNBC_Total_Quant_Singlets) * 100
Swarbrick_CID3963_TNBC_Total_Quant <- as.data.frame(c(Swarbrick_CID3963_TNBC_Total_Quant_Singlets, Swarbrick_CID3963_TNBC_Total_Quant_Doublets, Swarbrick_CID3963_TNBC_Total_Quant_Doublets_Percent))
colnames(Swarbrick_CID3963_TNBC_Total_Quant) <- c("Swarbrick_CID3963_TNBC_Total_Quant")
rownames(Swarbrick_CID3963_TNBC_Total_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Swarbrick_CID3963_TNBC_Total_Quant, file = "/R/R_Swarbrick/Swarbrick_Doublet_Tables/Swarbrick_CID3963_TNBC_Total_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Swarbrick_CID3963_TNBC_Total_Singlets <- subset(Swarbrick_CID3963_TNBC_Total, cells=rownames(Swarbrick_CID3963_TNBC_Total@meta.data)[which(Swarbrick_CID3963_TNBC_Total@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Swarbrick_CID3963_TNBC_Total_Singlets, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Total_Singlets/Swarbrick_CID3963_TNBC_Total_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Swarbrick_CID3963_TNBC_Total)
rm(Swarbrick_CID3963_TNBC_Total_Singlets)
rm(Swarbrick_CID3963_TNBC_Total_Quant)
rm(Swarbrick_CID3963_TNBC_Total_Quant_Singlets)
rm(Swarbrick_CID3963_TNBC_Total_Quant_Doublets)
rm(Swarbrick_CID3963_TNBC_Total_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Swarbrick_CID3963_TNBC_Total)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Swarbrick_CID3963_TNBC_Total)
rm(sweep.stats.Swarbrick_CID3963_TNBC_Total)
gc()



################################################################################################
########################              Swarbrick_CID4465_TNBC_Total             #################
################################################################################################

# Load Data
Swarbrick_CID4465_TNBC_Total <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Total/Swarbrick_CID4465_TNBC_Total.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Swarbrick_CID4465_TNBC_Total <- paramSweep(Swarbrick_CID4465_TNBC_Total, PCs = 1:10, sct = TRUE)
sweep.stats.Swarbrick_CID4465_TNBC_Total <- summarizeSweep(sweep.res.Swarbrick_CID4465_TNBC_Total, GT = FALSE)
bcmvn_Swarbrick_CID4465_TNBC_Total <- find.pK(sweep.stats.Swarbrick_CID4465_TNBC_Total)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Swarbrick_CID4465_TNBC_Total[which.max(bcmvn_Swarbrick_CID4465_TNBC_Total$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Swarbrick_CID4465_TNBC_Total@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 5% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Swarbrick dataset = 5k per sample, will use doublet rate for 6,000 cells to err on the side of caution
nExp_poi <- round(0.05*nrow(Swarbrick_CID4465_TNBC_Total@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Swarbrick_CID4465_TNBC_Total <- doubletFinder(Swarbrick_CID4465_TNBC_Total, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Swarbrick_CID4465_TNBC_Total_Quant <- (Swarbrick_CID4465_TNBC_Total@meta.data$DF.classification == "Singlet")
Swarbrick_CID4465_TNBC_Total_Quant_Singlets <- length(Swarbrick_CID4465_TNBC_Total_Quant[Swarbrick_CID4465_TNBC_Total_Quant== TRUE])
Swarbrick_CID4465_TNBC_Total_Quant_Doublets <- length(Swarbrick_CID4465_TNBC_Total_Quant[Swarbrick_CID4465_TNBC_Total_Quant== FALSE])
Swarbrick_CID4465_TNBC_Total_Quant_Doublets_Percent <- Swarbrick_CID4465_TNBC_Total_Quant_Doublets / (Swarbrick_CID4465_TNBC_Total_Quant_Doublets + Swarbrick_CID4465_TNBC_Total_Quant_Singlets) * 100
Swarbrick_CID4465_TNBC_Total_Quant <- as.data.frame(c(Swarbrick_CID4465_TNBC_Total_Quant_Singlets, Swarbrick_CID4465_TNBC_Total_Quant_Doublets, Swarbrick_CID4465_TNBC_Total_Quant_Doublets_Percent))
colnames(Swarbrick_CID4465_TNBC_Total_Quant) <- c("Swarbrick_CID4465_TNBC_Total_Quant")
rownames(Swarbrick_CID4465_TNBC_Total_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Swarbrick_CID4465_TNBC_Total_Quant, file = "/R/R_Swarbrick/Swarbrick_Doublet_Tables/Swarbrick_CID4465_TNBC_Total_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Swarbrick_CID4465_TNBC_Total_Singlets <- subset(Swarbrick_CID4465_TNBC_Total, cells=rownames(Swarbrick_CID4465_TNBC_Total@meta.data)[which(Swarbrick_CID4465_TNBC_Total@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Swarbrick_CID4465_TNBC_Total_Singlets, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Total_Singlets/Swarbrick_CID4465_TNBC_Total_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Swarbrick_CID4465_TNBC_Total)
rm(Swarbrick_CID4465_TNBC_Total_Singlets)
rm(Swarbrick_CID4465_TNBC_Total_Quant)
rm(Swarbrick_CID4465_TNBC_Total_Quant_Singlets)
rm(Swarbrick_CID4465_TNBC_Total_Quant_Doublets)
rm(Swarbrick_CID4465_TNBC_Total_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Swarbrick_CID4465_TNBC_Total)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Swarbrick_CID4465_TNBC_Total)
rm(sweep.stats.Swarbrick_CID4465_TNBC_Total)
gc()



################################################################################################
########################              Swarbrick_CID4495_TNBC_Total             #################
################################################################################################

# Load Data
Swarbrick_CID4495_TNBC_Total <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Total/Swarbrick_CID4495_TNBC_Total.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Swarbrick_CID4495_TNBC_Total <- paramSweep(Swarbrick_CID4495_TNBC_Total, PCs = 1:10, sct = TRUE)
sweep.stats.Swarbrick_CID4495_TNBC_Total <- summarizeSweep(sweep.res.Swarbrick_CID4495_TNBC_Total, GT = FALSE)
bcmvn_Swarbrick_CID4495_TNBC_Total <- find.pK(sweep.stats.Swarbrick_CID4495_TNBC_Total)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Swarbrick_CID4495_TNBC_Total[which.max(bcmvn_Swarbrick_CID4495_TNBC_Total$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Swarbrick_CID4495_TNBC_Total@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 5% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Swarbrick dataset = 5k per sample, will use doublet rate for 6,000 cells to err on the side of caution
nExp_poi <- round(0.05*nrow(Swarbrick_CID4495_TNBC_Total@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Swarbrick_CID4495_TNBC_Total <- doubletFinder(Swarbrick_CID4495_TNBC_Total, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Swarbrick_CID4495_TNBC_Total_Quant <- (Swarbrick_CID4495_TNBC_Total@meta.data$DF.classification == "Singlet")
Swarbrick_CID4495_TNBC_Total_Quant_Singlets <- length(Swarbrick_CID4495_TNBC_Total_Quant[Swarbrick_CID4495_TNBC_Total_Quant== TRUE])
Swarbrick_CID4495_TNBC_Total_Quant_Doublets <- length(Swarbrick_CID4495_TNBC_Total_Quant[Swarbrick_CID4495_TNBC_Total_Quant== FALSE])
Swarbrick_CID4495_TNBC_Total_Quant_Doublets_Percent <- Swarbrick_CID4495_TNBC_Total_Quant_Doublets / (Swarbrick_CID4495_TNBC_Total_Quant_Doublets + Swarbrick_CID4495_TNBC_Total_Quant_Singlets) * 100
Swarbrick_CID4495_TNBC_Total_Quant <- as.data.frame(c(Swarbrick_CID4495_TNBC_Total_Quant_Singlets, Swarbrick_CID4495_TNBC_Total_Quant_Doublets, Swarbrick_CID4495_TNBC_Total_Quant_Doublets_Percent))
colnames(Swarbrick_CID4495_TNBC_Total_Quant) <- c("Swarbrick_CID4495_TNBC_Total_Quant")
rownames(Swarbrick_CID4495_TNBC_Total_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Swarbrick_CID4495_TNBC_Total_Quant, file = "/R/R_Swarbrick/Swarbrick_Doublet_Tables/Swarbrick_CID4495_TNBC_Total_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Swarbrick_CID4495_TNBC_Total_Singlets <- subset(Swarbrick_CID4495_TNBC_Total, cells=rownames(Swarbrick_CID4495_TNBC_Total@meta.data)[which(Swarbrick_CID4495_TNBC_Total@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Swarbrick_CID4495_TNBC_Total_Singlets, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Total_Singlets/Swarbrick_CID4495_TNBC_Total_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Swarbrick_CID4495_TNBC_Total)
rm(Swarbrick_CID4495_TNBC_Total_Singlets)
rm(Swarbrick_CID4495_TNBC_Total_Quant)
rm(Swarbrick_CID4495_TNBC_Total_Quant_Singlets)
rm(Swarbrick_CID4495_TNBC_Total_Quant_Doublets)
rm(Swarbrick_CID4495_TNBC_Total_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Swarbrick_CID4495_TNBC_Total)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Swarbrick_CID4495_TNBC_Total)
rm(sweep.stats.Swarbrick_CID4495_TNBC_Total)
gc()



################################################################################################
########################              Swarbrick_CID4513_TNBC_Total             #################
################################################################################################

# Load Data
Swarbrick_CID4513_TNBC_Total <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Total/Swarbrick_CID4513_TNBC_Total.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Swarbrick_CID4513_TNBC_Total <- paramSweep(Swarbrick_CID4513_TNBC_Total, PCs = 1:10, sct = TRUE)
sweep.stats.Swarbrick_CID4513_TNBC_Total <- summarizeSweep(sweep.res.Swarbrick_CID4513_TNBC_Total, GT = FALSE)
bcmvn_Swarbrick_CID4513_TNBC_Total <- find.pK(sweep.stats.Swarbrick_CID4513_TNBC_Total)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Swarbrick_CID4513_TNBC_Total[which.max(bcmvn_Swarbrick_CID4513_TNBC_Total$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Swarbrick_CID4513_TNBC_Total@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 5% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Swarbrick dataset = 5k per sample, will use doublet rate for 6,000 cells to err on the side of caution
nExp_poi <- round(0.05*nrow(Swarbrick_CID4513_TNBC_Total@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Swarbrick_CID4513_TNBC_Total <- doubletFinder(Swarbrick_CID4513_TNBC_Total, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Swarbrick_CID4513_TNBC_Total_Quant <- (Swarbrick_CID4513_TNBC_Total@meta.data$DF.classification == "Singlet")
Swarbrick_CID4513_TNBC_Total_Quant_Singlets <- length(Swarbrick_CID4513_TNBC_Total_Quant[Swarbrick_CID4513_TNBC_Total_Quant== TRUE])
Swarbrick_CID4513_TNBC_Total_Quant_Doublets <- length(Swarbrick_CID4513_TNBC_Total_Quant[Swarbrick_CID4513_TNBC_Total_Quant== FALSE])
Swarbrick_CID4513_TNBC_Total_Quant_Doublets_Percent <- Swarbrick_CID4513_TNBC_Total_Quant_Doublets / (Swarbrick_CID4513_TNBC_Total_Quant_Doublets + Swarbrick_CID4513_TNBC_Total_Quant_Singlets) * 100
Swarbrick_CID4513_TNBC_Total_Quant <- as.data.frame(c(Swarbrick_CID4513_TNBC_Total_Quant_Singlets, Swarbrick_CID4513_TNBC_Total_Quant_Doublets, Swarbrick_CID4513_TNBC_Total_Quant_Doublets_Percent))
colnames(Swarbrick_CID4513_TNBC_Total_Quant) <- c("Swarbrick_CID4513_TNBC_Total_Quant")
rownames(Swarbrick_CID4513_TNBC_Total_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Swarbrick_CID4513_TNBC_Total_Quant, file = "/R/R_Swarbrick/Swarbrick_Doublet_Tables/Swarbrick_CID4513_TNBC_Total_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Swarbrick_CID4513_TNBC_Total_Singlets <- subset(Swarbrick_CID4513_TNBC_Total, cells=rownames(Swarbrick_CID4513_TNBC_Total@meta.data)[which(Swarbrick_CID4513_TNBC_Total@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Swarbrick_CID4513_TNBC_Total_Singlets, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Total_Singlets/Swarbrick_CID4513_TNBC_Total_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Swarbrick_CID4513_TNBC_Total)
rm(Swarbrick_CID4513_TNBC_Total_Singlets)
rm(Swarbrick_CID4513_TNBC_Total_Quant)
rm(Swarbrick_CID4513_TNBC_Total_Quant_Singlets)
rm(Swarbrick_CID4513_TNBC_Total_Quant_Doublets)
rm(Swarbrick_CID4513_TNBC_Total_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Swarbrick_CID4513_TNBC_Total)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Swarbrick_CID4513_TNBC_Total)
rm(sweep.stats.Swarbrick_CID4513_TNBC_Total)
gc()



################################################################################################
########################              Swarbrick_CID4515_TNBC_Total             #################
################################################################################################

# Load Data
Swarbrick_CID4515_TNBC_Total <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Total/Swarbrick_CID4515_TNBC_Total.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Swarbrick_CID4515_TNBC_Total <- paramSweep(Swarbrick_CID4515_TNBC_Total, PCs = 1:10, sct = TRUE)
sweep.stats.Swarbrick_CID4515_TNBC_Total <- summarizeSweep(sweep.res.Swarbrick_CID4515_TNBC_Total, GT = FALSE)
bcmvn_Swarbrick_CID4515_TNBC_Total <- find.pK(sweep.stats.Swarbrick_CID4515_TNBC_Total)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Swarbrick_CID4515_TNBC_Total[which.max(bcmvn_Swarbrick_CID4515_TNBC_Total$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Swarbrick_CID4515_TNBC_Total@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 5% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Swarbrick dataset = 5k per sample, will use doublet rate for 6,000 cells to err on the side of caution
nExp_poi <- round(0.05*nrow(Swarbrick_CID4515_TNBC_Total@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Swarbrick_CID4515_TNBC_Total <- doubletFinder(Swarbrick_CID4515_TNBC_Total, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Swarbrick_CID4515_TNBC_Total_Quant <- (Swarbrick_CID4515_TNBC_Total@meta.data$DF.classification == "Singlet")
Swarbrick_CID4515_TNBC_Total_Quant_Singlets <- length(Swarbrick_CID4515_TNBC_Total_Quant[Swarbrick_CID4515_TNBC_Total_Quant== TRUE])
Swarbrick_CID4515_TNBC_Total_Quant_Doublets <- length(Swarbrick_CID4515_TNBC_Total_Quant[Swarbrick_CID4515_TNBC_Total_Quant== FALSE])
Swarbrick_CID4515_TNBC_Total_Quant_Doublets_Percent <- Swarbrick_CID4515_TNBC_Total_Quant_Doublets / (Swarbrick_CID4515_TNBC_Total_Quant_Doublets + Swarbrick_CID4515_TNBC_Total_Quant_Singlets) * 100
Swarbrick_CID4515_TNBC_Total_Quant <- as.data.frame(c(Swarbrick_CID4515_TNBC_Total_Quant_Singlets, Swarbrick_CID4515_TNBC_Total_Quant_Doublets, Swarbrick_CID4515_TNBC_Total_Quant_Doublets_Percent))
colnames(Swarbrick_CID4515_TNBC_Total_Quant) <- c("Swarbrick_CID4515_TNBC_Total_Quant")
rownames(Swarbrick_CID4515_TNBC_Total_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Swarbrick_CID4515_TNBC_Total_Quant, file = "/R/R_Swarbrick/Swarbrick_Doublet_Tables/Swarbrick_CID4515_TNBC_Total_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Swarbrick_CID4515_TNBC_Total_Singlets <- subset(Swarbrick_CID4515_TNBC_Total, cells=rownames(Swarbrick_CID4515_TNBC_Total@meta.data)[which(Swarbrick_CID4515_TNBC_Total@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Swarbrick_CID4515_TNBC_Total_Singlets, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Total_Singlets/Swarbrick_CID4515_TNBC_Total_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Swarbrick_CID4515_TNBC_Total)
rm(Swarbrick_CID4515_TNBC_Total_Singlets)
rm(Swarbrick_CID4515_TNBC_Total_Quant)
rm(Swarbrick_CID4515_TNBC_Total_Quant_Singlets)
rm(Swarbrick_CID4515_TNBC_Total_Quant_Doublets)
rm(Swarbrick_CID4515_TNBC_Total_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Swarbrick_CID4515_TNBC_Total)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Swarbrick_CID4515_TNBC_Total)
rm(sweep.stats.Swarbrick_CID4515_TNBC_Total)
gc()



################################################################################################
########################              Swarbrick_CID4523_TNBC_Total             #################
################################################################################################

# Load Data
Swarbrick_CID4523_TNBC_Total <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Total/Swarbrick_CID4523_TNBC_Total.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Swarbrick_CID4523_TNBC_Total <- paramSweep(Swarbrick_CID4523_TNBC_Total, PCs = 1:10, sct = TRUE)
sweep.stats.Swarbrick_CID4523_TNBC_Total <- summarizeSweep(sweep.res.Swarbrick_CID4523_TNBC_Total, GT = FALSE)
bcmvn_Swarbrick_CID4523_TNBC_Total <- find.pK(sweep.stats.Swarbrick_CID4523_TNBC_Total)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Swarbrick_CID4523_TNBC_Total[which.max(bcmvn_Swarbrick_CID4523_TNBC_Total$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Swarbrick_CID4523_TNBC_Total@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 5% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Swarbrick dataset = 5k per sample, will use doublet rate for 6,000 cells to err on the side of caution
nExp_poi <- round(0.05*nrow(Swarbrick_CID4523_TNBC_Total@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Swarbrick_CID4523_TNBC_Total <- doubletFinder(Swarbrick_CID4523_TNBC_Total, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Swarbrick_CID4523_TNBC_Total_Quant <- (Swarbrick_CID4523_TNBC_Total@meta.data$DF.classification == "Singlet")
Swarbrick_CID4523_TNBC_Total_Quant_Singlets <- length(Swarbrick_CID4523_TNBC_Total_Quant[Swarbrick_CID4523_TNBC_Total_Quant== TRUE])
Swarbrick_CID4523_TNBC_Total_Quant_Doublets <- length(Swarbrick_CID4523_TNBC_Total_Quant[Swarbrick_CID4523_TNBC_Total_Quant== FALSE])
Swarbrick_CID4523_TNBC_Total_Quant_Doublets_Percent <- Swarbrick_CID4523_TNBC_Total_Quant_Doublets / (Swarbrick_CID4523_TNBC_Total_Quant_Doublets + Swarbrick_CID4523_TNBC_Total_Quant_Singlets) * 100
Swarbrick_CID4523_TNBC_Total_Quant <- as.data.frame(c(Swarbrick_CID4523_TNBC_Total_Quant_Singlets, Swarbrick_CID4523_TNBC_Total_Quant_Doublets, Swarbrick_CID4523_TNBC_Total_Quant_Doublets_Percent))
colnames(Swarbrick_CID4523_TNBC_Total_Quant) <- c("Swarbrick_CID4523_TNBC_Total_Quant")
rownames(Swarbrick_CID4523_TNBC_Total_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Swarbrick_CID4523_TNBC_Total_Quant, file = "/R/R_Swarbrick/Swarbrick_Doublet_Tables/Swarbrick_CID4523_TNBC_Total_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Swarbrick_CID4523_TNBC_Total_Singlets <- subset(Swarbrick_CID4523_TNBC_Total, cells=rownames(Swarbrick_CID4523_TNBC_Total@meta.data)[which(Swarbrick_CID4523_TNBC_Total@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Swarbrick_CID4523_TNBC_Total_Singlets, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Total_Singlets/Swarbrick_CID4523_TNBC_Total_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Swarbrick_CID4523_TNBC_Total)
rm(Swarbrick_CID4523_TNBC_Total_Singlets)
rm(Swarbrick_CID4523_TNBC_Total_Quant)
rm(Swarbrick_CID4523_TNBC_Total_Quant_Singlets)
rm(Swarbrick_CID4523_TNBC_Total_Quant_Doublets)
rm(Swarbrick_CID4523_TNBC_Total_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Swarbrick_CID4523_TNBC_Total)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Swarbrick_CID4523_TNBC_Total)
rm(sweep.stats.Swarbrick_CID4523_TNBC_Total)
gc()



################################################################################################
########################              Swarbrick_CID44041_TNBC_Total             ################
################################################################################################

# Load Data
Swarbrick_CID44041_TNBC_Total <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Total/Swarbrick_CID44041_TNBC_Total.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Swarbrick_CID44041_TNBC_Total <- paramSweep(Swarbrick_CID44041_TNBC_Total, PCs = 1:10, sct = TRUE)
sweep.stats.Swarbrick_CID44041_TNBC_Total <- summarizeSweep(sweep.res.Swarbrick_CID44041_TNBC_Total, GT = FALSE)
bcmvn_Swarbrick_CID44041_TNBC_Total <- find.pK(sweep.stats.Swarbrick_CID44041_TNBC_Total)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Swarbrick_CID44041_TNBC_Total[which.max(bcmvn_Swarbrick_CID44041_TNBC_Total$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Swarbrick_CID44041_TNBC_Total@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 5% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Swarbrick dataset = 5k per sample, will use doublet rate for 6,000 cells to err on the side of caution
nExp_poi <- round(0.05*nrow(Swarbrick_CID44041_TNBC_Total@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Swarbrick_CID44041_TNBC_Total <- doubletFinder(Swarbrick_CID44041_TNBC_Total, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Swarbrick_CID44041_TNBC_Total_Quant <- (Swarbrick_CID44041_TNBC_Total@meta.data$DF.classification == "Singlet")
Swarbrick_CID44041_TNBC_Total_Quant_Singlets <- length(Swarbrick_CID44041_TNBC_Total_Quant[Swarbrick_CID44041_TNBC_Total_Quant== TRUE])
Swarbrick_CID44041_TNBC_Total_Quant_Doublets <- length(Swarbrick_CID44041_TNBC_Total_Quant[Swarbrick_CID44041_TNBC_Total_Quant== FALSE])
Swarbrick_CID44041_TNBC_Total_Quant_Doublets_Percent <- Swarbrick_CID44041_TNBC_Total_Quant_Doublets / (Swarbrick_CID44041_TNBC_Total_Quant_Doublets + Swarbrick_CID44041_TNBC_Total_Quant_Singlets) * 100
Swarbrick_CID44041_TNBC_Total_Quant <- as.data.frame(c(Swarbrick_CID44041_TNBC_Total_Quant_Singlets, Swarbrick_CID44041_TNBC_Total_Quant_Doublets, Swarbrick_CID44041_TNBC_Total_Quant_Doublets_Percent))
colnames(Swarbrick_CID44041_TNBC_Total_Quant) <- c("Swarbrick_CID44041_TNBC_Total_Quant")
rownames(Swarbrick_CID44041_TNBC_Total_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Swarbrick_CID44041_TNBC_Total_Quant, file = "/R/R_Swarbrick/Swarbrick_Doublet_Tables/Swarbrick_CID44041_TNBC_Total_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Swarbrick_CID44041_TNBC_Total_Singlets <- subset(Swarbrick_CID44041_TNBC_Total, cells=rownames(Swarbrick_CID44041_TNBC_Total@meta.data)[which(Swarbrick_CID44041_TNBC_Total@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Swarbrick_CID44041_TNBC_Total_Singlets, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Total_Singlets/Swarbrick_CID44041_TNBC_Total_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Swarbrick_CID44041_TNBC_Total)
rm(Swarbrick_CID44041_TNBC_Total_Singlets)
rm(Swarbrick_CID44041_TNBC_Total_Quant)
rm(Swarbrick_CID44041_TNBC_Total_Quant_Singlets)
rm(Swarbrick_CID44041_TNBC_Total_Quant_Doublets)
rm(Swarbrick_CID44041_TNBC_Total_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Swarbrick_CID44041_TNBC_Total)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Swarbrick_CID44041_TNBC_Total)
rm(sweep.stats.Swarbrick_CID44041_TNBC_Total)
gc()





################################################################################################
########################              Swarbrick_CID44971_TNBC_Total             ################
################################################################################################

# Load Data
Swarbrick_CID44971_TNBC_Total <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Total/Swarbrick_CID44971_TNBC_Total.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Swarbrick_CID44971_TNBC_Total <- paramSweep(Swarbrick_CID44971_TNBC_Total, PCs = 1:10, sct = TRUE)
sweep.stats.Swarbrick_CID44971_TNBC_Total <- summarizeSweep(sweep.res.Swarbrick_CID44971_TNBC_Total, GT = FALSE)
bcmvn_Swarbrick_CID44971_TNBC_Total <- find.pK(sweep.stats.Swarbrick_CID44971_TNBC_Total)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Swarbrick_CID44971_TNBC_Total[which.max(bcmvn_Swarbrick_CID44971_TNBC_Total$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Swarbrick_CID44971_TNBC_Total@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 5% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Swarbrick dataset = 5k per sample, will use doublet rate for 6,000 cells to err on the side of caution
nExp_poi <- round(0.05*nrow(Swarbrick_CID44971_TNBC_Total@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Swarbrick_CID44971_TNBC_Total <- doubletFinder(Swarbrick_CID44971_TNBC_Total, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Swarbrick_CID44971_TNBC_Total_Quant <- (Swarbrick_CID44971_TNBC_Total@meta.data$DF.classification == "Singlet")
Swarbrick_CID44971_TNBC_Total_Quant_Singlets <- length(Swarbrick_CID44971_TNBC_Total_Quant[Swarbrick_CID44971_TNBC_Total_Quant== TRUE])
Swarbrick_CID44971_TNBC_Total_Quant_Doublets <- length(Swarbrick_CID44971_TNBC_Total_Quant[Swarbrick_CID44971_TNBC_Total_Quant== FALSE])
Swarbrick_CID44971_TNBC_Total_Quant_Doublets_Percent <- Swarbrick_CID44971_TNBC_Total_Quant_Doublets / (Swarbrick_CID44971_TNBC_Total_Quant_Doublets + Swarbrick_CID44971_TNBC_Total_Quant_Singlets) * 100
Swarbrick_CID44971_TNBC_Total_Quant <- as.data.frame(c(Swarbrick_CID44971_TNBC_Total_Quant_Singlets, Swarbrick_CID44971_TNBC_Total_Quant_Doublets, Swarbrick_CID44971_TNBC_Total_Quant_Doublets_Percent))
colnames(Swarbrick_CID44971_TNBC_Total_Quant) <- c("Swarbrick_CID44971_TNBC_Total_Quant")
rownames(Swarbrick_CID44971_TNBC_Total_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Swarbrick_CID44971_TNBC_Total_Quant, file = "/R/R_Swarbrick/Swarbrick_Doublet_Tables/Swarbrick_CID44971_TNBC_Total_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Swarbrick_CID44971_TNBC_Total_Singlets <- subset(Swarbrick_CID44971_TNBC_Total, cells=rownames(Swarbrick_CID44971_TNBC_Total@meta.data)[which(Swarbrick_CID44971_TNBC_Total@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Swarbrick_CID44971_TNBC_Total_Singlets, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Total_Singlets/Swarbrick_CID44971_TNBC_Total_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Swarbrick_CID44971_TNBC_Total)
rm(Swarbrick_CID44971_TNBC_Total_Singlets)
rm(Swarbrick_CID44971_TNBC_Total_Quant)
rm(Swarbrick_CID44971_TNBC_Total_Quant_Singlets)
rm(Swarbrick_CID44971_TNBC_Total_Quant_Doublets)
rm(Swarbrick_CID44971_TNBC_Total_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Swarbrick_CID44971_TNBC_Total)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Swarbrick_CID44971_TNBC_Total)
rm(sweep.stats.Swarbrick_CID44971_TNBC_Total)
gc()




################################################################################################
########################              Swarbrick_CID44991_TNBC_Total             ################
################################################################################################

# Load Data
Swarbrick_CID44991_TNBC_Total <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Total/Swarbrick_CID44991_TNBC_Total.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.Swarbrick_CID44991_TNBC_Total <- paramSweep(Swarbrick_CID44991_TNBC_Total, PCs = 1:10, sct = TRUE)
sweep.stats.Swarbrick_CID44991_TNBC_Total <- summarizeSweep(sweep.res.Swarbrick_CID44991_TNBC_Total, GT = FALSE)
bcmvn_Swarbrick_CID44991_TNBC_Total <- find.pK(sweep.stats.Swarbrick_CID44991_TNBC_Total)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_Swarbrick_CID44991_TNBC_Total[which.max(bcmvn_Swarbrick_CID44991_TNBC_Total$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Swarbrick_CID44991_TNBC_Total@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 5% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the Swarbrick dataset = 5k per sample, will use doublet rate for 6,000 cells to err on the side of caution
nExp_poi <- round(0.05*nrow(Swarbrick_CID44991_TNBC_Total@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
Swarbrick_CID44991_TNBC_Total <- doubletFinder(Swarbrick_CID44991_TNBC_Total, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
Swarbrick_CID44991_TNBC_Total_Quant <- (Swarbrick_CID44991_TNBC_Total@meta.data$DF.classification == "Singlet")
Swarbrick_CID44991_TNBC_Total_Quant_Singlets <- length(Swarbrick_CID44991_TNBC_Total_Quant[Swarbrick_CID44991_TNBC_Total_Quant== TRUE])
Swarbrick_CID44991_TNBC_Total_Quant_Doublets <- length(Swarbrick_CID44991_TNBC_Total_Quant[Swarbrick_CID44991_TNBC_Total_Quant== FALSE])
Swarbrick_CID44991_TNBC_Total_Quant_Doublets_Percent <- Swarbrick_CID44991_TNBC_Total_Quant_Doublets / (Swarbrick_CID44991_TNBC_Total_Quant_Doublets + Swarbrick_CID44991_TNBC_Total_Quant_Singlets) * 100
Swarbrick_CID44991_TNBC_Total_Quant <- as.data.frame(c(Swarbrick_CID44991_TNBC_Total_Quant_Singlets, Swarbrick_CID44991_TNBC_Total_Quant_Doublets, Swarbrick_CID44991_TNBC_Total_Quant_Doublets_Percent))
colnames(Swarbrick_CID44991_TNBC_Total_Quant) <- c("Swarbrick_CID44991_TNBC_Total_Quant")
rownames(Swarbrick_CID44991_TNBC_Total_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(Swarbrick_CID44991_TNBC_Total_Quant, file = "/R/R_Swarbrick/Swarbrick_Doublet_Tables/Swarbrick_CID44991_TNBC_Total_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
Swarbrick_CID44991_TNBC_Total_Singlets <- subset(Swarbrick_CID44991_TNBC_Total, cells=rownames(Swarbrick_CID44991_TNBC_Total@meta.data)[which(Swarbrick_CID44991_TNBC_Total@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(Swarbrick_CID44991_TNBC_Total_Singlets, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Total_Singlets/Swarbrick_CID44991_TNBC_Total_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(Swarbrick_CID44991_TNBC_Total)
rm(Swarbrick_CID44991_TNBC_Total_Singlets)
rm(Swarbrick_CID44991_TNBC_Total_Quant)
rm(Swarbrick_CID44991_TNBC_Total_Quant_Singlets)
rm(Swarbrick_CID44991_TNBC_Total_Quant_Doublets)
rm(Swarbrick_CID44991_TNBC_Total_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_Swarbrick_CID44991_TNBC_Total)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.Swarbrick_CID44991_TNBC_Total)
rm(sweep.stats.Swarbrick_CID44991_TNBC_Total)
gc()





