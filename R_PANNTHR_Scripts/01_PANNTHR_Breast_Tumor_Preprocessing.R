#####################################################################################################################
#                              PANNTHR Dataset Breast Tumor Analysis Steps                                          #
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

#########################################################
############# This Pre-Treatment Batch  # ###############
#########################################################

###############################################################
########################### PANNTHR_Pt1_PreTreat ##############
###############################################################

# Load the dataset
PANNTHR_Pt1_PreTreat <- Read10X(data.dir = "/R/R_PANNTHR/PANNTHR_Input/Pre_treatment/PANNTHR_Pt1_PreTreat/filtered_feature_bc_matrix/")

# Initialize the Seurat object with the raw (non-normalized data).
PANNTHR_Pt1_PreTreat <- CreateSeuratObject(counts = PANNTHR_Pt1_PreTreat, project = "PANNTHR_Pt1_PreTreat", min.cells = 3, min.features = 200)
PANNTHR_Pt1_PreTreat

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
PANNTHR_Pt1_PreTreat[["percent.mt"]] <- PercentageFeatureSet(PANNTHR_Pt1_PreTreat, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(PANNTHR_Pt1_PreTreat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(PANNTHR_Pt1_PreTreat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(PANNTHR_Pt1_PreTreat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
PANNTHR_Pt1_PreTreat <- subset(PANNTHR_Pt1_PreTreat, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & percent.mt < 10)

# Normalize data
PANNTHR_Pt1_PreTreat <- NormalizeData(PANNTHR_Pt1_PreTreat, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
PANNTHR_Pt1_PreTreat <- FindVariableFeatures(PANNTHR_Pt1_PreTreat, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(PANNTHR_Pt1_PreTreat)
PANNTHR_Pt1_PreTreat <- ScaleData(PANNTHR_Pt1_PreTreat, features = all.genes)

# Perform linear dimensional reduction
PANNTHR_Pt1_PreTreat <- RunPCA(PANNTHR_Pt1_PreTreat, features = VariableFeatures(object = PANNTHR_Pt1_PreTreat))

# Examine and visualize PCA results a few different ways
print(PANNTHR_Pt1_PreTreat[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(PANNTHR_Pt1_PreTreat, dims = 1:2, reduction = "pca")
DimPlot(PANNTHR_Pt1_PreTreat, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
PANNTHR_Pt1_PreTreat <- JackStraw(PANNTHR_Pt1_PreTreat, num.replicate = 100)
PANNTHR_Pt1_PreTreat <- ScoreJackStraw(PANNTHR_Pt1_PreTreat, dims = 1:20)
ElbowPlot(PANNTHR_Pt1_PreTreat)

# Cluster cells
PANNTHR_Pt1_PreTreat <- FindNeighbors(PANNTHR_Pt1_PreTreat, dims = 1:20)
PANNTHR_Pt1_PreTreat <- FindClusters(PANNTHR_Pt1_PreTreat, resolution = 0.5)

# Run non-linear dimensional reduction
PANNTHR_Pt1_PreTreat <- RunUMAP(PANNTHR_Pt1_PreTreat, dims = 1:20)

# Visualize clusters
DimPlot(PANNTHR_Pt1_PreTreat, reduction = "umap")

# Save RDS
saveRDS(PANNTHR_Pt1_PreTreat, file = "/R/R_PANNTHR/PANNTHR_RDS/RDS_Total/PANNTHR_Pt1_PreTreat.rds")



###############################################################
########################### PANNTHR_Pt2_PreTreat ##############
###############################################################

# Load the dataset
PANNTHR_Pt2_PreTreat <- Read10X(data.dir = "/R/R_PANNTHR/PANNTHR_Input/Pre_treatment/PANNTHR_Pt2_PreTreat/filtered_feature_bc_matrix/")

# Initialize the Seurat object with the raw (non-normalized data).
PANNTHR_Pt2_PreTreat <- CreateSeuratObject(counts = PANNTHR_Pt2_PreTreat, project = "PANNTHR_Pt2_PreTreat", min.cells = 3, min.features = 200)
PANNTHR_Pt2_PreTreat

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
PANNTHR_Pt2_PreTreat[["percent.mt"]] <- PercentageFeatureSet(PANNTHR_Pt2_PreTreat, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(PANNTHR_Pt2_PreTreat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(PANNTHR_Pt2_PreTreat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(PANNTHR_Pt2_PreTreat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
PANNTHR_Pt2_PreTreat <- subset(PANNTHR_Pt2_PreTreat, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 10)

# Normalize data
PANNTHR_Pt2_PreTreat <- NormalizeData(PANNTHR_Pt2_PreTreat, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
PANNTHR_Pt2_PreTreat <- FindVariableFeatures(PANNTHR_Pt2_PreTreat, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(PANNTHR_Pt2_PreTreat)
PANNTHR_Pt2_PreTreat <- ScaleData(PANNTHR_Pt2_PreTreat, features = all.genes)

# Perform linear dimensional reduction
PANNTHR_Pt2_PreTreat <- RunPCA(PANNTHR_Pt2_PreTreat, features = VariableFeatures(object = PANNTHR_Pt2_PreTreat))

# Examine and visualize PCA results a few different ways
print(PANNTHR_Pt2_PreTreat[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(PANNTHR_Pt2_PreTreat, dims = 1:2, reduction = "pca")
DimPlot(PANNTHR_Pt2_PreTreat, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
PANNTHR_Pt2_PreTreat <- JackStraw(PANNTHR_Pt2_PreTreat, num.replicate = 100)
PANNTHR_Pt2_PreTreat <- ScoreJackStraw(PANNTHR_Pt2_PreTreat, dims = 1:20)
ElbowPlot(PANNTHR_Pt2_PreTreat)

# Cluster cells
PANNTHR_Pt2_PreTreat <- FindNeighbors(PANNTHR_Pt2_PreTreat, dims = 1:20)
PANNTHR_Pt2_PreTreat <- FindClusters(PANNTHR_Pt2_PreTreat, resolution = 0.5)

# Run non-linear dimensional reduction
PANNTHR_Pt2_PreTreat <- RunUMAP(PANNTHR_Pt2_PreTreat, dims = 1:20)

# Visualize clusters
DimPlot(PANNTHR_Pt2_PreTreat, reduction = "umap")

# Save RDS
saveRDS(PANNTHR_Pt2_PreTreat, file = "/R/R_PANNTHR/PANNTHR_RDS/RDS_Total/PANNTHR_Pt2_PreTreat.rds")



###############################################################
########################### PANNTHR_Pt4_PreTreat ##############
###############################################################

# Load the dataset
PANNTHR_Pt4_PreTreat <- Read10X(data.dir = "/R/R_PANNTHR/PANNTHR_Input/Pre_treatment/PANNTHR_Pt4_PreTreat/filtered_feature_bc_matrix/")

# Initialize the Seurat object with the raw (non-normalized data).
PANNTHR_Pt4_PreTreat <- CreateSeuratObject(counts = PANNTHR_Pt4_PreTreat, project = "PANNTHR_Pt4_PreTreat", min.cells = 3, min.features = 200)
PANNTHR_Pt4_PreTreat

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
PANNTHR_Pt4_PreTreat[["percent.mt"]] <- PercentageFeatureSet(PANNTHR_Pt4_PreTreat, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(PANNTHR_Pt4_PreTreat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(PANNTHR_Pt4_PreTreat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(PANNTHR_Pt4_PreTreat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
PANNTHR_Pt4_PreTreat <- subset(PANNTHR_Pt4_PreTreat, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)

# Normalize data
PANNTHR_Pt4_PreTreat <- NormalizeData(PANNTHR_Pt4_PreTreat, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
PANNTHR_Pt4_PreTreat <- FindVariableFeatures(PANNTHR_Pt4_PreTreat, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(PANNTHR_Pt4_PreTreat)
PANNTHR_Pt4_PreTreat <- ScaleData(PANNTHR_Pt4_PreTreat, features = all.genes)

# Perform linear dimensional reduction
PANNTHR_Pt4_PreTreat <- RunPCA(PANNTHR_Pt4_PreTreat, features = VariableFeatures(object = PANNTHR_Pt4_PreTreat))

# Examine and visualize PCA results a few different ways
print(PANNTHR_Pt4_PreTreat[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(PANNTHR_Pt4_PreTreat, dims = 1:2, reduction = "pca")
DimPlot(PANNTHR_Pt4_PreTreat, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
PANNTHR_Pt4_PreTreat <- JackStraw(PANNTHR_Pt4_PreTreat, num.replicate = 100)
PANNTHR_Pt4_PreTreat <- ScoreJackStraw(PANNTHR_Pt4_PreTreat, dims = 1:20)
ElbowPlot(PANNTHR_Pt4_PreTreat)

# Cluster cells
PANNTHR_Pt4_PreTreat <- FindNeighbors(PANNTHR_Pt4_PreTreat, dims = 1:20)
PANNTHR_Pt4_PreTreat <- FindClusters(PANNTHR_Pt4_PreTreat, resolution = 0.5)

# Run non-linear dimensional reduction
PANNTHR_Pt4_PreTreat <- RunUMAP(PANNTHR_Pt4_PreTreat, dims = 1:20)

# Visualize clusters
DimPlot(PANNTHR_Pt4_PreTreat, reduction = "umap")

# Save RDS
saveRDS(PANNTHR_Pt4_PreTreat, file = "/R/R_PANNTHR/PANNTHR_RDS/RDS_Total/PANNTHR_Pt4_PreTreat.rds")



###############################################################
########################### PANNTHR_Pt5_PreTreat ##############
###############################################################

# Load the dataset
PANNTHR_Pt5_PreTreat <- Read10X(data.dir = "/R/R_PANNTHR/PANNTHR_Input/Pre_treatment/PANNTHR_Pt5_PreTreat/filtered_feature_bc_matrix/")

# Initialize the Seurat object with the raw (non-normalized data).
PANNTHR_Pt5_PreTreat <- CreateSeuratObject(counts = PANNTHR_Pt5_PreTreat, project = "PANNTHR_Pt5_PreTreat", min.cells = 3, min.features = 200)
PANNTHR_Pt5_PreTreat

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
PANNTHR_Pt5_PreTreat[["percent.mt"]] <- PercentageFeatureSet(PANNTHR_Pt5_PreTreat, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(PANNTHR_Pt5_PreTreat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(PANNTHR_Pt5_PreTreat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(PANNTHR_Pt5_PreTreat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
PANNTHR_Pt5_PreTreat <- subset(PANNTHR_Pt5_PreTreat, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & percent.mt < 10)

# Normalize data
PANNTHR_Pt5_PreTreat <- NormalizeData(PANNTHR_Pt5_PreTreat, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
PANNTHR_Pt5_PreTreat <- FindVariableFeatures(PANNTHR_Pt5_PreTreat, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(PANNTHR_Pt5_PreTreat)
PANNTHR_Pt5_PreTreat <- ScaleData(PANNTHR_Pt5_PreTreat, features = all.genes)

# Perform linear dimensional reduction
PANNTHR_Pt5_PreTreat <- RunPCA(PANNTHR_Pt5_PreTreat, features = VariableFeatures(object = PANNTHR_Pt5_PreTreat))

# Examine and visualize PCA results a few different ways
print(PANNTHR_Pt5_PreTreat[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(PANNTHR_Pt5_PreTreat, dims = 1:2, reduction = "pca")
DimPlot(PANNTHR_Pt5_PreTreat, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
PANNTHR_Pt5_PreTreat <- JackStraw(PANNTHR_Pt5_PreTreat, num.replicate = 100)
PANNTHR_Pt5_PreTreat <- ScoreJackStraw(PANNTHR_Pt5_PreTreat, dims = 1:20)
ElbowPlot(PANNTHR_Pt5_PreTreat)

# Cluster cells
PANNTHR_Pt5_PreTreat <- FindNeighbors(PANNTHR_Pt5_PreTreat, dims = 1:20)
PANNTHR_Pt5_PreTreat <- FindClusters(PANNTHR_Pt5_PreTreat, resolution = 0.5)

# Run non-linear dimensional reduction
PANNTHR_Pt5_PreTreat <- RunUMAP(PANNTHR_Pt5_PreTreat, dims = 1:20)

# Visualize clusters
DimPlot(PANNTHR_Pt5_PreTreat, reduction = "umap")

# Save RDS
saveRDS(PANNTHR_Pt5_PreTreat, file = "/R/R_PANNTHR/PANNTHR_RDS/RDS_Total/PANNTHR_Pt5_PreTreat.rds")



#########################################################
#############  This On-Treatment Batch   ################
#########################################################

##############################################################
########################### PANNTHR_Pt2_OnTreat ##############
##############################################################

# Load the dataset
PANNTHR_Pt2_OnTreat <- Read10X(data.dir = "/R/R_PANNTHR/PANNTHR_Input/On_treatment/PANNTHR_Pt2_OnTreat/filtered_feature_bc_matrix/")

# Initialize the Seurat object with the raw (non-normalized data).
PANNTHR_Pt2_OnTreat <- CreateSeuratObject(counts = PANNTHR_Pt2_OnTreat, project = "PANNTHR_Pt2_OnTreat", min.cells = 3, min.features = 200)
PANNTHR_Pt2_OnTreat

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
PANNTHR_Pt2_OnTreat[["percent.mt"]] <- PercentageFeatureSet(PANNTHR_Pt2_OnTreat, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(PANNTHR_Pt2_OnTreat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(PANNTHR_Pt2_OnTreat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(PANNTHR_Pt2_OnTreat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
PANNTHR_Pt2_OnTreat <- subset(PANNTHR_Pt2_OnTreat, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & percent.mt < 10)

# Normalize data
PANNTHR_Pt2_OnTreat <- NormalizeData(PANNTHR_Pt2_OnTreat, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
PANNTHR_Pt2_OnTreat <- FindVariableFeatures(PANNTHR_Pt2_OnTreat, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(PANNTHR_Pt2_OnTreat)
PANNTHR_Pt2_OnTreat <- ScaleData(PANNTHR_Pt2_OnTreat, features = all.genes)

# Perform linear dimensional reduction
PANNTHR_Pt2_OnTreat <- RunPCA(PANNTHR_Pt2_OnTreat, features = VariableFeatures(object = PANNTHR_Pt2_OnTreat))

# Examine and visualize PCA results a few different ways
print(PANNTHR_Pt2_OnTreat[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(PANNTHR_Pt2_OnTreat, dims = 1:2, reduction = "pca")
DimPlot(PANNTHR_Pt2_OnTreat, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
PANNTHR_Pt2_OnTreat <- JackStraw(PANNTHR_Pt2_OnTreat, num.replicate = 100)
PANNTHR_Pt2_OnTreat <- ScoreJackStraw(PANNTHR_Pt2_OnTreat, dims = 1:20)
ElbowPlot(PANNTHR_Pt2_OnTreat)

# Cluster cells
PANNTHR_Pt2_OnTreat <- FindNeighbors(PANNTHR_Pt2_OnTreat, dims = 1:20)
PANNTHR_Pt2_OnTreat <- FindClusters(PANNTHR_Pt2_OnTreat, resolution = 0.5)

# Run non-linear dimensional reduction
PANNTHR_Pt2_OnTreat <- RunUMAP(PANNTHR_Pt2_OnTreat, dims = 1:20)

# Visualize clusters
DimPlot(PANNTHR_Pt2_OnTreat, reduction = "umap")

# Save RDS
saveRDS(PANNTHR_Pt2_OnTreat, file = "/R/R_PANNTHR/PANNTHR_RDS/RDS_Total/PANNTHR_Pt2_OnTreat.rds")



##############################################################
########################### PANNTHR_Pt3_OnTreat ##############
##############################################################

# Load the dataset
PANNTHR_Pt3_OnTreat <- Read10X(data.dir = "/R/R_PANNTHR/PANNTHR_Input/On_treatment/PANNTHR_Pt3_OnTreat/filtered_feature_bc_matrix/")

# Initialize the Seurat object with the raw (non-normalized data).
PANNTHR_Pt3_OnTreat <- CreateSeuratObject(counts = PANNTHR_Pt3_OnTreat, project = "PANNTHR_Pt3_OnTreat", min.cells = 3, min.features = 200)
PANNTHR_Pt3_OnTreat

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
PANNTHR_Pt3_OnTreat[["percent.mt"]] <- PercentageFeatureSet(PANNTHR_Pt3_OnTreat, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(PANNTHR_Pt3_OnTreat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(PANNTHR_Pt3_OnTreat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(PANNTHR_Pt3_OnTreat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
PANNTHR_Pt3_OnTreat <- subset(PANNTHR_Pt3_OnTreat, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & percent.mt < 10)

# Normalize data
PANNTHR_Pt3_OnTreat <- NormalizeData(PANNTHR_Pt3_OnTreat, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
PANNTHR_Pt3_OnTreat <- FindVariableFeatures(PANNTHR_Pt3_OnTreat, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(PANNTHR_Pt3_OnTreat)
PANNTHR_Pt3_OnTreat <- ScaleData(PANNTHR_Pt3_OnTreat, features = all.genes)

# Perform linear dimensional reduction
PANNTHR_Pt3_OnTreat <- RunPCA(PANNTHR_Pt3_OnTreat, features = VariableFeatures(object = PANNTHR_Pt3_OnTreat))

# Examine and visualize PCA results a few different ways
print(PANNTHR_Pt3_OnTreat[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(PANNTHR_Pt3_OnTreat, dims = 1:2, reduction = "pca")
DimPlot(PANNTHR_Pt3_OnTreat, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
PANNTHR_Pt3_OnTreat <- JackStraw(PANNTHR_Pt3_OnTreat, num.replicate = 100)
PANNTHR_Pt3_OnTreat <- ScoreJackStraw(PANNTHR_Pt3_OnTreat, dims = 1:20)
ElbowPlot(PANNTHR_Pt3_OnTreat)

# Cluster cells
PANNTHR_Pt3_OnTreat <- FindNeighbors(PANNTHR_Pt3_OnTreat, dims = 1:20)
PANNTHR_Pt3_OnTreat <- FindClusters(PANNTHR_Pt3_OnTreat, resolution = 0.5)

# Run non-linear dimensional reduction
PANNTHR_Pt3_OnTreat <- RunUMAP(PANNTHR_Pt3_OnTreat, dims = 1:20)

# Visualize clusters
DimPlot(PANNTHR_Pt3_OnTreat, reduction = "umap")

# Save RDS
saveRDS(PANNTHR_Pt3_OnTreat, file = "/R/R_PANNTHR/PANNTHR_RDS/RDS_Total/PANNTHR_Pt3_OnTreat.rds")



##############################################################
########################### PANNTHR_Pt4_OnTreat ##############
##############################################################

# Load the dataset
PANNTHR_Pt4_OnTreat <- Read10X(data.dir = "/R/R_PANNTHR/PANNTHR_Input/On_treatment/PANNTHR_Pt4_OnTreat/filtered_feature_bc_matrix/")

# Initialize the Seurat object with the raw (non-normalized data).
PANNTHR_Pt4_OnTreat <- CreateSeuratObject(counts = PANNTHR_Pt4_OnTreat, project = "PANNTHR_Pt4_OnTreat", min.cells = 3, min.features = 200)
PANNTHR_Pt4_OnTreat

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
PANNTHR_Pt4_OnTreat[["percent.mt"]] <- PercentageFeatureSet(PANNTHR_Pt4_OnTreat, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(PANNTHR_Pt4_OnTreat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(PANNTHR_Pt4_OnTreat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(PANNTHR_Pt4_OnTreat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
PANNTHR_Pt4_OnTreat <- subset(PANNTHR_Pt4_OnTreat, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & percent.mt < 10)

# Normalize data
PANNTHR_Pt4_OnTreat <- NormalizeData(PANNTHR_Pt4_OnTreat, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
PANNTHR_Pt4_OnTreat <- FindVariableFeatures(PANNTHR_Pt4_OnTreat, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(PANNTHR_Pt4_OnTreat)
PANNTHR_Pt4_OnTreat <- ScaleData(PANNTHR_Pt4_OnTreat, features = all.genes)

# Perform linear dimensional reduction
PANNTHR_Pt4_OnTreat <- RunPCA(PANNTHR_Pt4_OnTreat, features = VariableFeatures(object = PANNTHR_Pt4_OnTreat))

# Examine and visualize PCA results a few different ways
print(PANNTHR_Pt4_OnTreat[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(PANNTHR_Pt4_OnTreat, dims = 1:2, reduction = "pca")
DimPlot(PANNTHR_Pt4_OnTreat, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
PANNTHR_Pt4_OnTreat <- JackStraw(PANNTHR_Pt4_OnTreat, num.replicate = 100)
PANNTHR_Pt4_OnTreat <- ScoreJackStraw(PANNTHR_Pt4_OnTreat, dims = 1:20)
ElbowPlot(PANNTHR_Pt4_OnTreat)

# Cluster cells
PANNTHR_Pt4_OnTreat <- FindNeighbors(PANNTHR_Pt4_OnTreat, dims = 1:20)
PANNTHR_Pt4_OnTreat <- FindClusters(PANNTHR_Pt4_OnTreat, resolution = 0.5)

# Run non-linear dimensional reduction
PANNTHR_Pt4_OnTreat <- RunUMAP(PANNTHR_Pt4_OnTreat, dims = 1:20)

# Visualize clusters
DimPlot(PANNTHR_Pt4_OnTreat, reduction = "umap")

# Save RDS
saveRDS(PANNTHR_Pt4_OnTreat, file = "/R/R_PANNTHR/PANNTHR_RDS/RDS_Total/PANNTHR_Pt4_OnTreat.rds")



##############################################################
########################### PANNTHR_Pt6_OnTreat ##############
##############################################################

# Load the dataset
PANNTHR_Pt6_OnTreat <- Read10X(data.dir = "/R/R_PANNTHR/PANNTHR_Input/On_treatment/PANNTHR_Pt6_OnTreat/filtered_feature_bc_matrix/")

# Initialize the Seurat object with the raw (non-normalized data).
PANNTHR_Pt6_OnTreat <- CreateSeuratObject(counts = PANNTHR_Pt6_OnTreat, project = "PANNTHR_Pt6_OnTreat", min.cells = 3, min.features = 200)
PANNTHR_Pt6_OnTreat

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
PANNTHR_Pt6_OnTreat[["percent.mt"]] <- PercentageFeatureSet(PANNTHR_Pt6_OnTreat, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(PANNTHR_Pt6_OnTreat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(PANNTHR_Pt6_OnTreat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(PANNTHR_Pt6_OnTreat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
PANNTHR_Pt6_OnTreat <- subset(PANNTHR_Pt6_OnTreat, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
PANNTHR_Pt6_OnTreat <- NormalizeData(PANNTHR_Pt6_OnTreat, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
PANNTHR_Pt6_OnTreat <- FindVariableFeatures(PANNTHR_Pt6_OnTreat, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(PANNTHR_Pt6_OnTreat)
PANNTHR_Pt6_OnTreat <- ScaleData(PANNTHR_Pt6_OnTreat, features = all.genes)

# Perform linear dimensional reduction
PANNTHR_Pt6_OnTreat <- RunPCA(PANNTHR_Pt6_OnTreat, features = VariableFeatures(object = PANNTHR_Pt6_OnTreat))

# Examine and visualize PCA results a few different ways
print(PANNTHR_Pt6_OnTreat[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(PANNTHR_Pt6_OnTreat, dims = 1:2, reduction = "pca")
DimPlot(PANNTHR_Pt6_OnTreat, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
PANNTHR_Pt6_OnTreat <- JackStraw(PANNTHR_Pt6_OnTreat, num.replicate = 100)
PANNTHR_Pt6_OnTreat <- ScoreJackStraw(PANNTHR_Pt6_OnTreat, dims = 1:20)
ElbowPlot(PANNTHR_Pt6_OnTreat)

# Cluster cells
PANNTHR_Pt6_OnTreat <- FindNeighbors(PANNTHR_Pt6_OnTreat, dims = 1:20)
PANNTHR_Pt6_OnTreat <- FindClusters(PANNTHR_Pt6_OnTreat, resolution = 0.5)

# Run non-linear dimensional reduction
PANNTHR_Pt6_OnTreat <- RunUMAP(PANNTHR_Pt6_OnTreat, dims = 1:20)

# Visualize clusters
DimPlot(PANNTHR_Pt6_OnTreat, reduction = "umap")

# Save RDS
saveRDS(PANNTHR_Pt6_OnTreat, file = "/R/R_PANNTHR/PANNTHR_RDS/RDS_Total/PANNTHR_Pt6_OnTreat.rds")



##############################################################
########################### PANNTHR_Pt7_OnTreat ##############
##############################################################

# Load the dataset
PANNTHR_Pt7_OnTreat <- Read10X(data.dir = "/R/R_PANNTHR/PANNTHR_Input/On_treatment/PANNTHR_Pt7_OnTreat/filtered_feature_bc_matrix/")

# Initialize the Seurat object with the raw (non-normalized data).
PANNTHR_Pt7_OnTreat <- CreateSeuratObject(counts = PANNTHR_Pt7_OnTreat, project = "PANNTHR_Pt7_OnTreat", min.cells = 3, min.features = 200)
PANNTHR_Pt7_OnTreat

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
PANNTHR_Pt7_OnTreat[["percent.mt"]] <- PercentageFeatureSet(PANNTHR_Pt7_OnTreat, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(PANNTHR_Pt7_OnTreat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(PANNTHR_Pt7_OnTreat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(PANNTHR_Pt7_OnTreat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
PANNTHR_Pt7_OnTreat <- subset(PANNTHR_Pt7_OnTreat, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
PANNTHR_Pt7_OnTreat <- NormalizeData(PANNTHR_Pt7_OnTreat, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
PANNTHR_Pt7_OnTreat <- FindVariableFeatures(PANNTHR_Pt7_OnTreat, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(PANNTHR_Pt7_OnTreat)
PANNTHR_Pt7_OnTreat <- ScaleData(PANNTHR_Pt7_OnTreat, features = all.genes)

# Perform linear dimensional reduction
PANNTHR_Pt7_OnTreat <- RunPCA(PANNTHR_Pt7_OnTreat, features = VariableFeatures(object = PANNTHR_Pt7_OnTreat))

# Examine and visualize PCA results a few different ways
print(PANNTHR_Pt7_OnTreat[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(PANNTHR_Pt7_OnTreat, dims = 1:2, reduction = "pca")
DimPlot(PANNTHR_Pt7_OnTreat, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
PANNTHR_Pt7_OnTreat <- JackStraw(PANNTHR_Pt7_OnTreat, num.replicate = 100)
PANNTHR_Pt7_OnTreat <- ScoreJackStraw(PANNTHR_Pt7_OnTreat, dims = 1:20)
ElbowPlot(PANNTHR_Pt7_OnTreat)

# Cluster cells
PANNTHR_Pt7_OnTreat <- FindNeighbors(PANNTHR_Pt7_OnTreat, dims = 1:20)
PANNTHR_Pt7_OnTreat <- FindClusters(PANNTHR_Pt7_OnTreat, resolution = 0.5)

# Run non-linear dimensional reduction
PANNTHR_Pt7_OnTreat <- RunUMAP(PANNTHR_Pt7_OnTreat, dims = 1:20)

# Visualize clusters
DimPlot(PANNTHR_Pt7_OnTreat, reduction = "umap")

# Save RDS
saveRDS(PANNTHR_Pt7_OnTreat, file = "/R/R_PANNTHR/PANNTHR_RDS/RDS_Total/PANNTHR_Pt7_OnTreat.rds")



#########################################################
############# This Post-Treatment Batch  ################
#########################################################

################################################################
########################### PANNTHR_Pt5_PostTreat ##############
################################################################

# Load the dataset
PANNTHR_Pt5_PostTreat <- Read10X(data.dir = "/R/R_PANNTHR/PANNTHR_Input/Post_treatment/PANNTHR_Pt5_PostTreat/filtered_feature_bc_matrix/")

# Initialize the Seurat object with the raw (non-normalized data).
PANNTHR_Pt5_PostTreat <- CreateSeuratObject(counts = PANNTHR_Pt5_PostTreat, project = "PANNTHR_Pt5_PostTreat", min.cells = 3, min.features = 200)
PANNTHR_Pt5_PostTreat

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
PANNTHR_Pt5_PostTreat[["percent.mt"]] <- PercentageFeatureSet(PANNTHR_Pt5_PostTreat, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(PANNTHR_Pt5_PostTreat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(PANNTHR_Pt5_PostTreat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(PANNTHR_Pt5_PostTreat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
PANNTHR_Pt5_PostTreat <- subset(PANNTHR_Pt5_PostTreat, subset = nFeature_RNA > 200 & nFeature_RNA < 9000 & percent.mt < 10)

# Normalize data
PANNTHR_Pt5_PostTreat <- NormalizeData(PANNTHR_Pt5_PostTreat, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
PANNTHR_Pt5_PostTreat <- FindVariableFeatures(PANNTHR_Pt5_PostTreat, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(PANNTHR_Pt5_PostTreat)
PANNTHR_Pt5_PostTreat <- ScaleData(PANNTHR_Pt5_PostTreat, features = all.genes)

# Perform linear dimensional reduction
PANNTHR_Pt5_PostTreat <- RunPCA(PANNTHR_Pt5_PostTreat, features = VariableFeatures(object = PANNTHR_Pt5_PostTreat))

# Examine and visualize PCA results a few different ways
print(PANNTHR_Pt5_PostTreat[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(PANNTHR_Pt5_PostTreat, dims = 1:2, reduction = "pca")
DimPlot(PANNTHR_Pt5_PostTreat, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
PANNTHR_Pt5_PostTreat <- JackStraw(PANNTHR_Pt5_PostTreat, num.replicate = 100)
PANNTHR_Pt5_PostTreat <- ScoreJackStraw(PANNTHR_Pt5_PostTreat, dims = 1:20)
ElbowPlot(PANNTHR_Pt5_PostTreat)

# Cluster cells
PANNTHR_Pt5_PostTreat <- FindNeighbors(PANNTHR_Pt5_PostTreat, dims = 1:20)
PANNTHR_Pt5_PostTreat <- FindClusters(PANNTHR_Pt5_PostTreat, resolution = 0.5)

# Run non-linear dimensional reduction
PANNTHR_Pt5_PostTreat <- RunUMAP(PANNTHR_Pt5_PostTreat, dims = 1:20)

# Visualize clusters
DimPlot(PANNTHR_Pt5_PostTreat, reduction = "umap")

# Save RDS
saveRDS(PANNTHR_Pt5_PostTreat, file = "/R/R_PANNTHR/PANNTHR_RDS/RDS_Total/PANNTHR_Pt5_PostTreat.rds")


PANNTHR_Pt5_PostTreat <- readRDS(file = "/R/R_PANNTHR/PANNTHR_RDS/RDS_Total/PANNTHR_Pt5_PostTreat.rds")


#####################################################################
########### Step 2: Run DoubletFinder and Subset Singlets ########### 
#####################################################################

# Load Libraries
library(Seurat)
library(DoubletFinder)

# Tailoring doublet rate to number of cells per sample, due to highly variable cell counts
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-

# Key: cells (expected doublet rate, %)
########## Pre-Treatment
# PANNTHR_Pt4_PreTreat: 21575 (17.3)
# PANNTHR_Pt2_PreTreat: 29896 (23.9)
# PANNTHR_Pt1_PreTreat: 9658 (7.7)
# PANNTHR_Pt5_PreTreat: 7149 (5.7)
########## On-Treatment
# PANNTHR_Pt6_OnTreat: 1083 (0.866)
# PANNTHR_Pt4_OnTreat: 23784 (19.0)
# PANNTHR_Pt3_OnTreat: 8593 (6.9)
# PANNTHR_Pt2_OnTreat: 38949 (31.2)
# PANNTHR_Pt7_OnTreat: 195 (0.156)
########## Post-Treatment
# PANNTHR_Pt5_PostTreat: 1265 (1.01)


################################################################################################
########################             Pre-Treatment Tumors                #######################
################################################################################################

################################################################################################
########################              PANNTHR_Pt1_PreTreat              ########################
################################################################################################

# Load Data
PANNTHR_Pt1_PreTreat <- readRDS(file = "/R/R_PANNTHR/PANNTHR_RDS/RDS_Total/PANNTHR_Pt1_PreTreat.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.PANNTHR_Pt1_PreTreat <- paramSweep(PANNTHR_Pt1_PreTreat, PCs = 1:10, sct = FALSE)
sweep.stats.PANNTHR_Pt1_PreTreat <- summarizeSweep(sweep.res.PANNTHR_Pt1_PreTreat, GT = FALSE)
bcmvn_PANNTHR_Pt1_PreTreat <- find.pK(sweep.stats.PANNTHR_Pt1_PreTreat)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distributiPre
bcmvn.max <- bcmvn_PANNTHR_Pt1_PreTreat[which.max(bcmvn_PANNTHR_Pt1_PreTreat$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet ProportiPre Estimate -------------------------------------------------------------------------------------
annotatiPres <- PANNTHR_Pt1_PreTreat@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotatiPres)
## Assuming 3% doublet formatiPre rate based Pre table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Tailoring doublet rate to cell number using linear regression of doublet table
nExp_poi <- round(0.08*nrow(PANNTHR_Pt1_PreTreat@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
PANNTHR_Pt1_PreTreat <- doubletFinder(PANNTHR_Pt1_PreTreat, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
PANNTHR_Pt1_PreTreat_Quant <- (PANNTHR_Pt1_PreTreat@meta.data$DF.classification == "Singlet")
PANNTHR_Pt1_PreTreat_Quant_Singlets <- length(PANNTHR_Pt1_PreTreat_Quant[PANNTHR_Pt1_PreTreat_Quant== TRUE])
PANNTHR_Pt1_PreTreat_Quant_Doublets <- length(PANNTHR_Pt1_PreTreat_Quant[PANNTHR_Pt1_PreTreat_Quant== FALSE])
PANNTHR_Pt1_PreTreat_Quant_Doublets_Percent <- PANNTHR_Pt1_PreTreat_Quant_Doublets / (PANNTHR_Pt1_PreTreat_Quant_Doublets + PANNTHR_Pt1_PreTreat_Quant_Singlets) * 100
PANNTHR_Pt1_PreTreat_Quant <- as.data.frame(c(PANNTHR_Pt1_PreTreat_Quant_Singlets, PANNTHR_Pt1_PreTreat_Quant_Doublets, PANNTHR_Pt1_PreTreat_Quant_Doublets_Percent))
colnames(PANNTHR_Pt1_PreTreat_Quant) <- c("PANNTHR_Pt1_PreTreat_Quant")
rownames(PANNTHR_Pt1_PreTreat_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(PANNTHR_Pt1_PreTreat_Quant, file = "/R/R_PANNTHR/PANNTHR_Doublet_Tables/PANNTHR_Pt1_PreTreat_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
PANNTHR_Pt1_PreTreat_Singlets <- subset(PANNTHR_Pt1_PreTreat, cells=rownames(PANNTHR_Pt1_PreTreat@meta.data)[which(PANNTHR_Pt1_PreTreat@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(PANNTHR_Pt1_PreTreat_Singlets, file = "/R/R_PANNTHR/PANNTHR_RDS/RDS_Total_Singlets/PANNTHR_Pt1_PreTreat_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(PANNTHR_Pt1_PreTreat)
rm(PANNTHR_Pt1_PreTreat_Singlets)
rm(PANNTHR_Pt1_PreTreat_Quant)
rm(PANNTHR_Pt1_PreTreat_Quant_Singlets)
rm(PANNTHR_Pt1_PreTreat_Quant_Doublets)
rm(PANNTHR_Pt1_PreTreat_Quant_Doublets_Percent)
rm(annotatiPres)
rm(bcmvn_PANNTHR_Pt1_PreTreat)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.PANNTHR_Pt1_PreTreat)
rm(sweep.stats.PANNTHR_Pt1_PreTreat)
gc()



################################################################################################
########################              PANNTHR_Pt2_PreTreat              ########################
################################################################################################

# Load Data
PANNTHR_Pt2_PreTreat <- readRDS(file = "/R/R_PANNTHR/PANNTHR_RDS/RDS_Total/PANNTHR_Pt2_PreTreat.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.PANNTHR_Pt2_PreTreat <- paramSweep(PANNTHR_Pt2_PreTreat, PCs = 1:10, sct = FALSE)
sweep.stats.PANNTHR_Pt2_PreTreat <- summarizeSweep(sweep.res.PANNTHR_Pt2_PreTreat, GT = FALSE)
bcmvn_PANNTHR_Pt2_PreTreat <- find.pK(sweep.stats.PANNTHR_Pt2_PreTreat)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distributiPre
bcmvn.max <- bcmvn_PANNTHR_Pt2_PreTreat[which.max(bcmvn_PANNTHR_Pt2_PreTreat$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet ProportiPre Estimate -------------------------------------------------------------------------------------
annotatiPres <- PANNTHR_Pt2_PreTreat@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotatiPres)
## Assuming 3% doublet formatiPre rate based Pre table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Tailoring doublet rate to cell number using linear regression of doublet table
nExp_poi <- round(0.24*nrow(PANNTHR_Pt2_PreTreat@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
PANNTHR_Pt2_PreTreat <- doubletFinder(PANNTHR_Pt2_PreTreat, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
PANNTHR_Pt2_PreTreat_Quant <- (PANNTHR_Pt2_PreTreat@meta.data$DF.classification == "Singlet")
PANNTHR_Pt2_PreTreat_Quant_Singlets <- length(PANNTHR_Pt2_PreTreat_Quant[PANNTHR_Pt2_PreTreat_Quant== TRUE])
PANNTHR_Pt2_PreTreat_Quant_Doublets <- length(PANNTHR_Pt2_PreTreat_Quant[PANNTHR_Pt2_PreTreat_Quant== FALSE])
PANNTHR_Pt2_PreTreat_Quant_Doublets_Percent <- PANNTHR_Pt2_PreTreat_Quant_Doublets / (PANNTHR_Pt2_PreTreat_Quant_Doublets + PANNTHR_Pt2_PreTreat_Quant_Singlets) * 100
PANNTHR_Pt2_PreTreat_Quant <- as.data.frame(c(PANNTHR_Pt2_PreTreat_Quant_Singlets, PANNTHR_Pt2_PreTreat_Quant_Doublets, PANNTHR_Pt2_PreTreat_Quant_Doublets_Percent))
colnames(PANNTHR_Pt2_PreTreat_Quant) <- c("PANNTHR_Pt2_PreTreat_Quant")
rownames(PANNTHR_Pt2_PreTreat_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(PANNTHR_Pt2_PreTreat_Quant, file = "/R/R_PANNTHR/PANNTHR_Doublet_Tables/PANNTHR_Pt2_PreTreat_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
PANNTHR_Pt2_PreTreat_Singlets <- subset(PANNTHR_Pt2_PreTreat, cells=rownames(PANNTHR_Pt2_PreTreat@meta.data)[which(PANNTHR_Pt2_PreTreat@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(PANNTHR_Pt2_PreTreat_Singlets, file = "/R/R_PANNTHR/PANNTHR_RDS/RDS_Total_Singlets/PANNTHR_Pt2_PreTreat_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(PANNTHR_Pt2_PreTreat)
rm(PANNTHR_Pt2_PreTreat_Singlets)
rm(PANNTHR_Pt2_PreTreat_Quant)
rm(PANNTHR_Pt2_PreTreat_Quant_Singlets)
rm(PANNTHR_Pt2_PreTreat_Quant_Doublets)
rm(PANNTHR_Pt2_PreTreat_Quant_Doublets_Percent)
rm(annotatiPres)
rm(bcmvn_PANNTHR_Pt2_PreTreat)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.PANNTHR_Pt2_PreTreat)
rm(sweep.stats.PANNTHR_Pt2_PreTreat)
gc()



################################################################################################
########################              PANNTHR_Pt4_PreTreat              ########################
################################################################################################

# Load Data
PANNTHR_Pt4_PreTreat <- readRDS(file = "/R/R_PANNTHR/PANNTHR_RDS/RDS_Total/PANNTHR_Pt4_PreTreat.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.PANNTHR_Pt4_PreTreat <- paramSweep(PANNTHR_Pt4_PreTreat, PCs = 1:10, sct = FALSE)
sweep.stats.PANNTHR_Pt4_PreTreat <- summarizeSweep(sweep.res.PANNTHR_Pt4_PreTreat, GT = FALSE)
bcmvn_PANNTHR_Pt4_PreTreat <- find.pK(sweep.stats.PANNTHR_Pt4_PreTreat)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distributiPre
bcmvn.max <- bcmvn_PANNTHR_Pt4_PreTreat[which.max(bcmvn_PANNTHR_Pt4_PreTreat$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet ProportiPre Estimate -------------------------------------------------------------------------------------
annotatiPres <- PANNTHR_Pt4_PreTreat@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotatiPres)
## Assuming 3% doublet formatiPre rate based Pre table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Tailoring doublet rate to cell number using linear regression of doublet table
nExp_poi <- round(0.17*nrow(PANNTHR_Pt4_PreTreat@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
PANNTHR_Pt4_PreTreat <- doubletFinder(PANNTHR_Pt4_PreTreat, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
PANNTHR_Pt4_PreTreat_Quant <- (PANNTHR_Pt4_PreTreat@meta.data$DF.classification == "Singlet")
PANNTHR_Pt4_PreTreat_Quant_Singlets <- length(PANNTHR_Pt4_PreTreat_Quant[PANNTHR_Pt4_PreTreat_Quant== TRUE])
PANNTHR_Pt4_PreTreat_Quant_Doublets <- length(PANNTHR_Pt4_PreTreat_Quant[PANNTHR_Pt4_PreTreat_Quant== FALSE])
PANNTHR_Pt4_PreTreat_Quant_Doublets_Percent <- PANNTHR_Pt4_PreTreat_Quant_Doublets / (PANNTHR_Pt4_PreTreat_Quant_Doublets + PANNTHR_Pt4_PreTreat_Quant_Singlets) * 100
PANNTHR_Pt4_PreTreat_Quant <- as.data.frame(c(PANNTHR_Pt4_PreTreat_Quant_Singlets, PANNTHR_Pt4_PreTreat_Quant_Doublets, PANNTHR_Pt4_PreTreat_Quant_Doublets_Percent))
colnames(PANNTHR_Pt4_PreTreat_Quant) <- c("PANNTHR_Pt4_PreTreat_Quant")
rownames(PANNTHR_Pt4_PreTreat_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(PANNTHR_Pt4_PreTreat_Quant, file = "/R/R_PANNTHR/PANNTHR_Doublet_Tables/PANNTHR_Pt4_PreTreat_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
PANNTHR_Pt4_PreTreat_Singlets <- subset(PANNTHR_Pt4_PreTreat, cells=rownames(PANNTHR_Pt4_PreTreat@meta.data)[which(PANNTHR_Pt4_PreTreat@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(PANNTHR_Pt4_PreTreat_Singlets, file = "/R/R_PANNTHR/PANNTHR_RDS/RDS_Total_Singlets/PANNTHR_Pt4_PreTreat_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(PANNTHR_Pt4_PreTreat)
rm(PANNTHR_Pt4_PreTreat_Singlets)
rm(PANNTHR_Pt4_PreTreat_Quant)
rm(PANNTHR_Pt4_PreTreat_Quant_Singlets)
rm(PANNTHR_Pt4_PreTreat_Quant_Doublets)
rm(PANNTHR_Pt4_PreTreat_Quant_Doublets_Percent)
rm(annotatiPres)
rm(bcmvn_PANNTHR_Pt4_PreTreat)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.PANNTHR_Pt4_PreTreat)
rm(sweep.stats.PANNTHR_Pt4_PreTreat)
gc()



################################################################################################
########################              PANNTHR_Pt5_PreTreat              ########################
################################################################################################

# Load Data
PANNTHR_Pt5_PreTreat <- readRDS(file = "/R/R_PANNTHR/PANNTHR_RDS/RDS_Total/PANNTHR_Pt5_PreTreat.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.PANNTHR_Pt5_PreTreat <- paramSweep(PANNTHR_Pt5_PreTreat, PCs = 1:10, sct = FALSE)
sweep.stats.PANNTHR_Pt5_PreTreat <- summarizeSweep(sweep.res.PANNTHR_Pt5_PreTreat, GT = FALSE)
bcmvn_PANNTHR_Pt5_PreTreat <- find.pK(sweep.stats.PANNTHR_Pt5_PreTreat)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distributiPre
bcmvn.max <- bcmvn_PANNTHR_Pt5_PreTreat[which.max(bcmvn_PANNTHR_Pt5_PreTreat$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet ProportiPre Estimate -------------------------------------------------------------------------------------
annotatiPres <- PANNTHR_Pt5_PreTreat@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotatiPres)
## Assuming 3% doublet formatiPre rate based Pre table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Tailoring doublet rate to cell number using linear regression of doublet table
nExp_poi <- round(0.06*nrow(PANNTHR_Pt5_PreTreat@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
PANNTHR_Pt5_PreTreat <- doubletFinder(PANNTHR_Pt5_PreTreat, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
PANNTHR_Pt5_PreTreat_Quant <- (PANNTHR_Pt5_PreTreat@meta.data$DF.classification == "Singlet")
PANNTHR_Pt5_PreTreat_Quant_Singlets <- length(PANNTHR_Pt5_PreTreat_Quant[PANNTHR_Pt5_PreTreat_Quant== TRUE])
PANNTHR_Pt5_PreTreat_Quant_Doublets <- length(PANNTHR_Pt5_PreTreat_Quant[PANNTHR_Pt5_PreTreat_Quant== FALSE])
PANNTHR_Pt5_PreTreat_Quant_Doublets_Percent <- PANNTHR_Pt5_PreTreat_Quant_Doublets / (PANNTHR_Pt5_PreTreat_Quant_Doublets + PANNTHR_Pt5_PreTreat_Quant_Singlets) * 100
PANNTHR_Pt5_PreTreat_Quant <- as.data.frame(c(PANNTHR_Pt5_PreTreat_Quant_Singlets, PANNTHR_Pt5_PreTreat_Quant_Doublets, PANNTHR_Pt5_PreTreat_Quant_Doublets_Percent))
colnames(PANNTHR_Pt5_PreTreat_Quant) <- c("PANNTHR_Pt5_PreTreat_Quant")
rownames(PANNTHR_Pt5_PreTreat_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(PANNTHR_Pt5_PreTreat_Quant, file = "/R/R_PANNTHR/PANNTHR_Doublet_Tables/PANNTHR_Pt5_PreTreat_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
PANNTHR_Pt5_PreTreat_Singlets <- subset(PANNTHR_Pt5_PreTreat, cells=rownames(PANNTHR_Pt5_PreTreat@meta.data)[which(PANNTHR_Pt5_PreTreat@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(PANNTHR_Pt5_PreTreat_Singlets, file = "/R/R_PANNTHR/PANNTHR_RDS/RDS_Total_Singlets/PANNTHR_Pt5_PreTreat_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(PANNTHR_Pt5_PreTreat)
rm(PANNTHR_Pt5_PreTreat_Singlets)
rm(PANNTHR_Pt5_PreTreat_Quant)
rm(PANNTHR_Pt5_PreTreat_Quant_Singlets)
rm(PANNTHR_Pt5_PreTreat_Quant_Doublets)
rm(PANNTHR_Pt5_PreTreat_Quant_Doublets_Percent)
rm(annotatiPres)
rm(bcmvn_PANNTHR_Pt5_PreTreat)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.PANNTHR_Pt5_PreTreat)
rm(sweep.stats.PANNTHR_Pt5_PreTreat)
gc()



################################################################################################
########################             On-Treatment Tumors                ########################
################################################################################################

################################################################################################
########################              PANNTHR_Pt2_OnTreat              ########################
################################################################################################

# Load Data
PANNTHR_Pt2_OnTreat <- readRDS(file = "/R/R_PANNTHR/PANNTHR_RDS/RDS_Total/PANNTHR_Pt2_OnTreat.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.PANNTHR_Pt2_OnTreat <- paramSweep(PANNTHR_Pt2_OnTreat, PCs = 1:10, sct = FALSE)
sweep.stats.PANNTHR_Pt2_OnTreat <- summarizeSweep(sweep.res.PANNTHR_Pt2_OnTreat, GT = FALSE)
bcmvn_PANNTHR_Pt2_OnTreat <- find.pK(sweep.stats.PANNTHR_Pt2_OnTreat)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_PANNTHR_Pt2_OnTreat[which.max(bcmvn_PANNTHR_Pt2_OnTreat$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- PANNTHR_Pt2_OnTreat@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 3% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the PANNTHR dataset = 5k per sample, will use doublet rate for 6,000 cells to err on the side of caution
# Tailoring doublet rate to cell number using linear regression of doublet table
nExp_poi <- round(0.31*nrow(PANNTHR_Pt2_OnTreat@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
PANNTHR_Pt2_OnTreat <- doubletFinder(PANNTHR_Pt2_OnTreat, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
PANNTHR_Pt2_OnTreat_Quant <- (PANNTHR_Pt2_OnTreat@meta.data$DF.classification == "Singlet")
PANNTHR_Pt2_OnTreat_Quant_Singlets <- length(PANNTHR_Pt2_OnTreat_Quant[PANNTHR_Pt2_OnTreat_Quant== TRUE])
PANNTHR_Pt2_OnTreat_Quant_Doublets <- length(PANNTHR_Pt2_OnTreat_Quant[PANNTHR_Pt2_OnTreat_Quant== FALSE])
PANNTHR_Pt2_OnTreat_Quant_Doublets_Percent <- PANNTHR_Pt2_OnTreat_Quant_Doublets / (PANNTHR_Pt2_OnTreat_Quant_Doublets + PANNTHR_Pt2_OnTreat_Quant_Singlets) * 100
PANNTHR_Pt2_OnTreat_Quant <- as.data.frame(c(PANNTHR_Pt2_OnTreat_Quant_Singlets, PANNTHR_Pt2_OnTreat_Quant_Doublets, PANNTHR_Pt2_OnTreat_Quant_Doublets_Percent))
colnames(PANNTHR_Pt2_OnTreat_Quant) <- c("PANNTHR_Pt2_OnTreat_Quant")
rownames(PANNTHR_Pt2_OnTreat_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(PANNTHR_Pt2_OnTreat_Quant, file = "/R/R_PANNTHR/PANNTHR_Doublet_Tables/PANNTHR_Pt2_OnTreat_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
PANNTHR_Pt2_OnTreat_Singlets <- subset(PANNTHR_Pt2_OnTreat, cells=rownames(PANNTHR_Pt2_OnTreat@meta.data)[which(PANNTHR_Pt2_OnTreat@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(PANNTHR_Pt2_OnTreat_Singlets, file = "/R/R_PANNTHR/PANNTHR_RDS/RDS_Total_Singlets/PANNTHR_Pt2_OnTreat_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(PANNTHR_Pt2_OnTreat)
rm(PANNTHR_Pt2_OnTreat_Singlets)
rm(PANNTHR_Pt2_OnTreat_Quant)
rm(PANNTHR_Pt2_OnTreat_Quant_Singlets)
rm(PANNTHR_Pt2_OnTreat_Quant_Doublets)
rm(PANNTHR_Pt2_OnTreat_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_PANNTHR_Pt2_OnTreat)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.PANNTHR_Pt2_OnTreat)
rm(sweep.stats.PANNTHR_Pt2_OnTreat)
gc()



################################################################################################
########################              PANNTHR_Pt3_OnTreat              ########################
################################################################################################

# Load Data
PANNTHR_Pt3_OnTreat <- readRDS(file = "/R/R_PANNTHR/PANNTHR_RDS/RDS_Total/PANNTHR_Pt3_OnTreat.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.PANNTHR_Pt3_OnTreat <- paramSweep(PANNTHR_Pt3_OnTreat, PCs = 1:10, sct = FALSE)
sweep.stats.PANNTHR_Pt3_OnTreat <- summarizeSweep(sweep.res.PANNTHR_Pt3_OnTreat, GT = FALSE)
bcmvn_PANNTHR_Pt3_OnTreat <- find.pK(sweep.stats.PANNTHR_Pt3_OnTreat)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_PANNTHR_Pt3_OnTreat[which.max(bcmvn_PANNTHR_Pt3_OnTreat$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- PANNTHR_Pt3_OnTreat@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 3% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the PANNTHR dataset = 5k per sample, will use doublet rate for 6,000 cells to err on the side of caution
# Tailoring doublet rate to cell number using linear regression of doublet table
nExp_poi <- round(0.07*nrow(PANNTHR_Pt3_OnTreat@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
PANNTHR_Pt3_OnTreat <- doubletFinder(PANNTHR_Pt3_OnTreat, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
PANNTHR_Pt3_OnTreat_Quant <- (PANNTHR_Pt3_OnTreat@meta.data$DF.classification == "Singlet")
PANNTHR_Pt3_OnTreat_Quant_Singlets <- length(PANNTHR_Pt3_OnTreat_Quant[PANNTHR_Pt3_OnTreat_Quant== TRUE])
PANNTHR_Pt3_OnTreat_Quant_Doublets <- length(PANNTHR_Pt3_OnTreat_Quant[PANNTHR_Pt3_OnTreat_Quant== FALSE])
PANNTHR_Pt3_OnTreat_Quant_Doublets_Percent <- PANNTHR_Pt3_OnTreat_Quant_Doublets / (PANNTHR_Pt3_OnTreat_Quant_Doublets + PANNTHR_Pt3_OnTreat_Quant_Singlets) * 100
PANNTHR_Pt3_OnTreat_Quant <- as.data.frame(c(PANNTHR_Pt3_OnTreat_Quant_Singlets, PANNTHR_Pt3_OnTreat_Quant_Doublets, PANNTHR_Pt3_OnTreat_Quant_Doublets_Percent))
colnames(PANNTHR_Pt3_OnTreat_Quant) <- c("PANNTHR_Pt3_OnTreat_Quant")
rownames(PANNTHR_Pt3_OnTreat_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(PANNTHR_Pt3_OnTreat_Quant, file = "/R/R_PANNTHR/PANNTHR_Doublet_Tables/PANNTHR_Pt3_OnTreat_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
PANNTHR_Pt3_OnTreat_Singlets <- subset(PANNTHR_Pt3_OnTreat, cells=rownames(PANNTHR_Pt3_OnTreat@meta.data)[which(PANNTHR_Pt3_OnTreat@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(PANNTHR_Pt3_OnTreat_Singlets, file = "/R/R_PANNTHR/PANNTHR_RDS/RDS_Total_Singlets/PANNTHR_Pt3_OnTreat_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(PANNTHR_Pt3_OnTreat)
rm(PANNTHR_Pt3_OnTreat_Singlets)
rm(PANNTHR_Pt3_OnTreat_Quant)
rm(PANNTHR_Pt3_OnTreat_Quant_Singlets)
rm(PANNTHR_Pt3_OnTreat_Quant_Doublets)
rm(PANNTHR_Pt3_OnTreat_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_PANNTHR_Pt3_OnTreat)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.PANNTHR_Pt3_OnTreat)
rm(sweep.stats.PANNTHR_Pt3_OnTreat)
gc()



################################################################################################
########################              PANNTHR_Pt4_OnTreat              ########################
################################################################################################

# Load Data
PANNTHR_Pt4_OnTreat <- readRDS(file = "/R/R_PANNTHR/PANNTHR_RDS/RDS_Total/PANNTHR_Pt4_OnTreat.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.PANNTHR_Pt4_OnTreat <- paramSweep(PANNTHR_Pt4_OnTreat, PCs = 1:10, sct = FALSE)
sweep.stats.PANNTHR_Pt4_OnTreat <- summarizeSweep(sweep.res.PANNTHR_Pt4_OnTreat, GT = FALSE)
bcmvn_PANNTHR_Pt4_OnTreat <- find.pK(sweep.stats.PANNTHR_Pt4_OnTreat)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_PANNTHR_Pt4_OnTreat[which.max(bcmvn_PANNTHR_Pt4_OnTreat$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- PANNTHR_Pt4_OnTreat@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 3% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the PANNTHR dataset = 5k per sample, will use doublet rate for 6,000 cells to err on the side of caution
# Tailoring doublet rate to cell number using linear regression of doublet table
nExp_poi <- round(0.19*nrow(PANNTHR_Pt4_OnTreat@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
PANNTHR_Pt4_OnTreat <- doubletFinder(PANNTHR_Pt4_OnTreat, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
PANNTHR_Pt4_OnTreat_Quant <- (PANNTHR_Pt4_OnTreat@meta.data$DF.classification == "Singlet")
PANNTHR_Pt4_OnTreat_Quant_Singlets <- length(PANNTHR_Pt4_OnTreat_Quant[PANNTHR_Pt4_OnTreat_Quant== TRUE])
PANNTHR_Pt4_OnTreat_Quant_Doublets <- length(PANNTHR_Pt4_OnTreat_Quant[PANNTHR_Pt4_OnTreat_Quant== FALSE])
PANNTHR_Pt4_OnTreat_Quant_Doublets_Percent <- PANNTHR_Pt4_OnTreat_Quant_Doublets / (PANNTHR_Pt4_OnTreat_Quant_Doublets + PANNTHR_Pt4_OnTreat_Quant_Singlets) * 100
PANNTHR_Pt4_OnTreat_Quant <- as.data.frame(c(PANNTHR_Pt4_OnTreat_Quant_Singlets, PANNTHR_Pt4_OnTreat_Quant_Doublets, PANNTHR_Pt4_OnTreat_Quant_Doublets_Percent))
colnames(PANNTHR_Pt4_OnTreat_Quant) <- c("PANNTHR_Pt4_OnTreat_Quant")
rownames(PANNTHR_Pt4_OnTreat_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(PANNTHR_Pt4_OnTreat_Quant, file = "/R/R_PANNTHR/PANNTHR_Doublet_Tables/PANNTHR_Pt4_OnTreat_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
PANNTHR_Pt4_OnTreat_Singlets <- subset(PANNTHR_Pt4_OnTreat, cells=rownames(PANNTHR_Pt4_OnTreat@meta.data)[which(PANNTHR_Pt4_OnTreat@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(PANNTHR_Pt4_OnTreat_Singlets, file = "/R/R_PANNTHR/PANNTHR_RDS/RDS_Total_Singlets/PANNTHR_Pt4_OnTreat_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(PANNTHR_Pt4_OnTreat)
rm(PANNTHR_Pt4_OnTreat_Singlets)
rm(PANNTHR_Pt4_OnTreat_Quant)
rm(PANNTHR_Pt4_OnTreat_Quant_Singlets)
rm(PANNTHR_Pt4_OnTreat_Quant_Doublets)
rm(PANNTHR_Pt4_OnTreat_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_PANNTHR_Pt4_OnTreat)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.PANNTHR_Pt4_OnTreat)
rm(sweep.stats.PANNTHR_Pt4_OnTreat)
gc()



################################################################################################
########################              PANNTHR_Pt6_OnTreat              ########################
################################################################################################

# Load Data
PANNTHR_Pt6_OnTreat <- readRDS(file = "/R/R_PANNTHR/PANNTHR_RDS/RDS_Total/PANNTHR_Pt6_OnTreat.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.PANNTHR_Pt6_OnTreat <- paramSweep(PANNTHR_Pt6_OnTreat, PCs = 1:10, sct = FALSE)
sweep.stats.PANNTHR_Pt6_OnTreat <- summarizeSweep(sweep.res.PANNTHR_Pt6_OnTreat, GT = FALSE)
bcmvn_PANNTHR_Pt6_OnTreat <- find.pK(sweep.stats.PANNTHR_Pt6_OnTreat)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_PANNTHR_Pt6_OnTreat[which.max(bcmvn_PANNTHR_Pt6_OnTreat$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- PANNTHR_Pt6_OnTreat@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 3% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the PANNTHR dataset = 5k per sample, will use doublet rate for 6,000 cells to err on the side of caution
# Tailoring doublet rate to cell number using linear regression of doublet table
nExp_poi <- round(0.009*nrow(PANNTHR_Pt6_OnTreat@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
PANNTHR_Pt6_OnTreat <- doubletFinder(PANNTHR_Pt6_OnTreat, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
PANNTHR_Pt6_OnTreat_Quant <- (PANNTHR_Pt6_OnTreat@meta.data$DF.classification == "Singlet")
PANNTHR_Pt6_OnTreat_Quant_Singlets <- length(PANNTHR_Pt6_OnTreat_Quant[PANNTHR_Pt6_OnTreat_Quant== TRUE])
PANNTHR_Pt6_OnTreat_Quant_Doublets <- length(PANNTHR_Pt6_OnTreat_Quant[PANNTHR_Pt6_OnTreat_Quant== FALSE])
PANNTHR_Pt6_OnTreat_Quant_Doublets_Percent <- PANNTHR_Pt6_OnTreat_Quant_Doublets / (PANNTHR_Pt6_OnTreat_Quant_Doublets + PANNTHR_Pt6_OnTreat_Quant_Singlets) * 100
PANNTHR_Pt6_OnTreat_Quant <- as.data.frame(c(PANNTHR_Pt6_OnTreat_Quant_Singlets, PANNTHR_Pt6_OnTreat_Quant_Doublets, PANNTHR_Pt6_OnTreat_Quant_Doublets_Percent))
colnames(PANNTHR_Pt6_OnTreat_Quant) <- c("PANNTHR_Pt6_OnTreat_Quant")
rownames(PANNTHR_Pt6_OnTreat_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(PANNTHR_Pt6_OnTreat_Quant, file = "/R/R_PANNTHR/PANNTHR_Doublet_Tables/PANNTHR_Pt6_OnTreat_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
PANNTHR_Pt6_OnTreat_Singlets <- subset(PANNTHR_Pt6_OnTreat, cells=rownames(PANNTHR_Pt6_OnTreat@meta.data)[which(PANNTHR_Pt6_OnTreat@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(PANNTHR_Pt6_OnTreat_Singlets, file = "/R/R_PANNTHR/PANNTHR_RDS/RDS_Total_Singlets/PANNTHR_Pt6_OnTreat_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(PANNTHR_Pt6_OnTreat)
rm(PANNTHR_Pt6_OnTreat_Singlets)
rm(PANNTHR_Pt6_OnTreat_Quant)
rm(PANNTHR_Pt6_OnTreat_Quant_Singlets)
rm(PANNTHR_Pt6_OnTreat_Quant_Doublets)
rm(PANNTHR_Pt6_OnTreat_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_PANNTHR_Pt6_OnTreat)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.PANNTHR_Pt6_OnTreat)
rm(sweep.stats.PANNTHR_Pt6_OnTreat)
gc()



################################################################################################
########################              PANNTHR_Pt7_OnTreat              ########################
################################################################################################

# Load Data
PANNTHR_Pt7_OnTreat <- readRDS(file = "/R/R_PANNTHR/PANNTHR_RDS/RDS_Total/PANNTHR_Pt7_OnTreat.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.PANNTHR_Pt7_OnTreat <- paramSweep(PANNTHR_Pt7_OnTreat, PCs = 1:10, sct = FALSE)
sweep.stats.PANNTHR_Pt7_OnTreat <- summarizeSweep(sweep.res.PANNTHR_Pt7_OnTreat, GT = FALSE)
bcmvn_PANNTHR_Pt7_OnTreat <- find.pK(sweep.stats.PANNTHR_Pt7_OnTreat)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn_PANNTHR_Pt7_OnTreat[which.max(bcmvn_PANNTHR_Pt7_OnTreat$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- PANNTHR_Pt7_OnTreat@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
## Assuming 3% doublet formation rate based on table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Average cells per sample in the PANNTHR dataset = 5k per sample, will use doublet rate for 6,000 cells to err on the side of caution
# Tailoring doublet rate to cell number using linear regression of doublet table
nExp_poi <- round(0.002*nrow(PANNTHR_Pt7_OnTreat@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
PANNTHR_Pt7_OnTreat <- doubletFinder(PANNTHR_Pt7_OnTreat, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
PANNTHR_Pt7_OnTreat_Quant <- (PANNTHR_Pt7_OnTreat@meta.data$DF.classification == "Singlet")
PANNTHR_Pt7_OnTreat_Quant_Singlets <- length(PANNTHR_Pt7_OnTreat_Quant[PANNTHR_Pt7_OnTreat_Quant== TRUE])
PANNTHR_Pt7_OnTreat_Quant_Doublets <- length(PANNTHR_Pt7_OnTreat_Quant[PANNTHR_Pt7_OnTreat_Quant== FALSE])
PANNTHR_Pt7_OnTreat_Quant_Doublets_Percent <- PANNTHR_Pt7_OnTreat_Quant_Doublets / (PANNTHR_Pt7_OnTreat_Quant_Doublets + PANNTHR_Pt7_OnTreat_Quant_Singlets) * 100
PANNTHR_Pt7_OnTreat_Quant <- as.data.frame(c(PANNTHR_Pt7_OnTreat_Quant_Singlets, PANNTHR_Pt7_OnTreat_Quant_Doublets, PANNTHR_Pt7_OnTreat_Quant_Doublets_Percent))
colnames(PANNTHR_Pt7_OnTreat_Quant) <- c("PANNTHR_Pt7_OnTreat_Quant")
rownames(PANNTHR_Pt7_OnTreat_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(PANNTHR_Pt7_OnTreat_Quant, file = "/R/R_PANNTHR/PANNTHR_Doublet_Tables/PANNTHR_Pt7_OnTreat_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
PANNTHR_Pt7_OnTreat_Singlets <- subset(PANNTHR_Pt7_OnTreat, cells=rownames(PANNTHR_Pt7_OnTreat@meta.data)[which(PANNTHR_Pt7_OnTreat@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(PANNTHR_Pt7_OnTreat_Singlets, file = "/R/R_PANNTHR/PANNTHR_RDS/RDS_Total_Singlets/PANNTHR_Pt7_OnTreat_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(PANNTHR_Pt7_OnTreat)
rm(PANNTHR_Pt7_OnTreat_Singlets)
rm(PANNTHR_Pt7_OnTreat_Quant)
rm(PANNTHR_Pt7_OnTreat_Quant_Singlets)
rm(PANNTHR_Pt7_OnTreat_Quant_Doublets)
rm(PANNTHR_Pt7_OnTreat_Quant_Doublets_Percent)
rm(annotations)
rm(bcmvn_PANNTHR_Pt7_OnTreat)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.PANNTHR_Pt7_OnTreat)
rm(sweep.stats.PANNTHR_Pt7_OnTreat)
gc()



################################################################################################
########################             Post-Treatment Tumors                ######################
################################################################################################

################################################################################################
########################              PANNTHR_Pt5_PostTreat              ########################
################################################################################################

# Load Data
PANNTHR_Pt5_PostTreat <- readRDS(file = "/R/R_PANNTHR/PANNTHR_RDS/RDS_Total/PANNTHR_Pt5_PostTreat.rds")

# DoubletFinder
# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.PANNTHR_Pt5_PostTreat <- paramSweep(PANNTHR_Pt5_PostTreat, PCs = 1:10, sct = FALSE)
sweep.stats.PANNTHR_Pt5_PostTreat <- summarizeSweep(sweep.res.PANNTHR_Pt5_PostTreat, GT = FALSE)
bcmvn_PANNTHR_Pt5_PostTreat <- find.pK(sweep.stats.PANNTHR_Pt5_PostTreat)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distributiPre
bcmvn.max <- bcmvn_PANNTHR_Pt5_PostTreat[which.max(bcmvn_PANNTHR_Pt5_PostTreat$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic Doublet ProportiPre Estimate -------------------------------------------------------------------------------------
annotatiPres <- PANNTHR_Pt5_PostTreat@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotatiPres)
## Assuming 3% doublet formatiPre rate based Pre table from 10X website
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
# Table indicates multiplet rate is ~2.4% for ~3k cells recovered, ~4% for 5k cells recovered; 
# Tailoring doublet rate to cell number using linear regression of doublet table
nExp_poi <- round(0.01*nrow(PANNTHR_Pt5_PostTreat@meta.data))   
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder ----------------------------------------------------------------------------------------------------------
PANNTHR_Pt5_PostTreat <- doubletFinder(PANNTHR_Pt5_PostTreat, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Quantify Doublets ----------------------------------------------------------------------------------------------------------
PANNTHR_Pt5_PostTreat_Quant <- (PANNTHR_Pt5_PostTreat@meta.data$DF.classification == "Singlet")
PANNTHR_Pt5_PostTreat_Quant_Singlets <- length(PANNTHR_Pt5_PostTreat_Quant[PANNTHR_Pt5_PostTreat_Quant== TRUE])
PANNTHR_Pt5_PostTreat_Quant_Doublets <- length(PANNTHR_Pt5_PostTreat_Quant[PANNTHR_Pt5_PostTreat_Quant== FALSE])
PANNTHR_Pt5_PostTreat_Quant_Doublets_Percent <- PANNTHR_Pt5_PostTreat_Quant_Doublets / (PANNTHR_Pt5_PostTreat_Quant_Doublets + PANNTHR_Pt5_PostTreat_Quant_Singlets) * 100
PANNTHR_Pt5_PostTreat_Quant <- as.data.frame(c(PANNTHR_Pt5_PostTreat_Quant_Singlets, PANNTHR_Pt5_PostTreat_Quant_Doublets, PANNTHR_Pt5_PostTreat_Quant_Doublets_Percent))
colnames(PANNTHR_Pt5_PostTreat_Quant) <- c("PANNTHR_Pt5_PostTreat_Quant")
rownames(PANNTHR_Pt5_PostTreat_Quant) <- c("Singlets","Doublets", "Percent_Doublets")
write.csv(PANNTHR_Pt5_PostTreat_Quant, file = "/R/R_PANNTHR/PANNTHR_Doublet_Tables/PANNTHR_Pt5_PostTreat_Quant.csv")

# Subset Singlets -------------------------------------------------------------------------------------------------------------
PANNTHR_Pt5_PostTreat_Singlets <- subset(PANNTHR_Pt5_PostTreat, cells=rownames(PANNTHR_Pt5_PostTreat@meta.data)[which(PANNTHR_Pt5_PostTreat@meta.data$DF.classification == "Singlet")])
# Save Singlet RDS
saveRDS(PANNTHR_Pt5_PostTreat_Singlets, file = "/R/R_PANNTHR/PANNTHR_RDS/RDS_Total_Singlets/PANNTHR_Pt5_PostTreat_Singlets.rds")

# Remove Objects to save memory ------------------------------------------------------------------------------------------------- 
rm(PANNTHR_Pt5_PostTreat)
rm(PANNTHR_Pt5_PostTreat_Singlets)
rm(PANNTHR_Pt5_PostTreat_Quant)
rm(PANNTHR_Pt5_PostTreat_Quant_Singlets)
rm(PANNTHR_Pt5_PostTreat_Quant_Doublets)
rm(PANNTHR_Pt5_PostTreat_Quant_Doublets_Percent)
rm(annotatiPres)
rm(bcmvn_PANNTHR_Pt5_PostTreat)
rm(bcmvn.max)
rm(homotypic.prop)
rm(nExp_poi)
rm(nExp_poi.adj)
rm(optimal.pk)
rm(sweep.res.PANNTHR_Pt5_PostTreat)
rm(sweep.stats.PANNTHR_Pt5_PostTreat)
gc()






