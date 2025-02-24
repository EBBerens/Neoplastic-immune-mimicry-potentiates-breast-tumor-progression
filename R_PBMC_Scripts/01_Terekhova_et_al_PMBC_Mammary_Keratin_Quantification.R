#####################################################################################################################
#                                      PBMC Dataset Analysis Steps                                                  #
#####################################################################################################################
# Step 1: Process Each Leukocyte Subset into Seurat Object                                                          #
# Step 2: Downsample Leukocyte Subsets and then Merge into Combined Seurat Object                                   #
# Step 3: Extract Expression Matrices for Mammary Keratin Genes per Leukocyte Subset                                #
#####################################################################################################################

# The PBMC dataset is very large
# Thus, Seurat example Umap plots are created from a downsampled object of 500,000 cells
# Next, expression matrices are extracted from individual leukoyte subsets to quantify cells with keratin expression; quantification done on the whole

############################################################
# Step 1: Process Each Leukocyte Subset into Seurat Object #
############################################################

#################
# Unzip Subsets #
#################
# Load Libraries 
library(utils)

# Unzip
untar("/R/R_PBMC/PBMC_Input/Artyomov_PBMC/Unzip/b_cells.tar.gz", compressed = 'gzip', exdir = "/R/R_PBMC/PBMC_Input/Artyomov_PBMC/Unzip/Extracted")
untar("/R/R_PBMC/PBMC_Input/Artyomov_PBMC/Unzip/cd4_t_cells.tar.tgz", compressed = 'gzip', exdir = "/R/R_PBMC/PBMC_Input/Artyomov_PBMC/Unzip/Extracted")
untar("/R/R_PBMC/PBMC_Input/Artyomov_PBMC/Unzip/cd4_t_helper_memory_cells.tar.tgz", compressed = 'gzip', exdir = "/R/R_PBMC/PBMC_Input/Artyomov_PBMC/Unzip/Extracted")
untar("/R/R_PBMC/PBMC_Input/Artyomov_PBMC/Unzip/conventional_cd8_t_cells.tar.tgz", compressed = 'gzip', exdir = "/R/R_PBMC/PBMC_Input/Artyomov_PBMC/Unzip/Extracted")
untar("/R/R_PBMC/PBMC_Input/Artyomov_PBMC/Unzip/gd_t_cells.tar.gz", compressed = 'gzip', exdir = "/R/R_PBMC/PBMC_Input/Artyomov_PBMC/Unzip/Extracted")
untar("/R/R_PBMC/PBMC_Input/Artyomov_PBMC/Unzip/mait_cells.tar.tgz", compressed = 'gzip', exdir = "/R/R_PBMC/PBMC_Input/Artyomov_PBMC/Unzip/Extracted")
untar("/R/R_PBMC/PBMC_Input/Artyomov_PBMC/Unzip/myeloid_cells.tar.gz", compressed = 'gzip', exdir = "/R/R_PBMC/PBMC_Input/Artyomov_PBMC/Unzip/Extracted")
untar("/R/R_PBMC/PBMC_Input/Artyomov_PBMC/Unzip/nk_cells.tar.gz", compressed = 'gzip', exdir = "/R/R_PBMC/PBMC_Input/Artyomov_PBMC/Unzip/Extracted")
untar("/R/R_PBMC/PBMC_Input/Artyomov_PBMC/Unzip/progenitor_cells.tar.tgz", compressed = 'gzip', exdir = "/R/R_PBMC/PBMC_Input/Artyomov_PBMC/Unzip/Extracted")

# Locate "rna.rds" files and place into /R/R_PBMC/PBMC_Input/Artyomov_PBMC/Subset_RDS/ 

######################################################
# Process Each Subset Separately into Seurat Objects #
######################################################
# Load Libraries 
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(BPCells)

#########################
# Load Individual Files #
#########################
b_cells <- readRDS(file = "/R/R_PBMC/PBMC_Input/Artyomov_PBMC/Subset_RDS/b_cells_rna.rds")
cd4_helper_memory_cells <- readRDS(file = "/R/R_PBMC/PBMC_Input/Artyomov_PBMC/Subset_RDS/cd4_helper_memory_rna.rds")
cd4_cells <- readRDS(file = "/R/R_PBMC/PBMC_Input/Artyomov_PBMC/Subset_RDS/cd4_rna.rds")
conventional_cd8_cells <- readRDS(file = "/R/R_PBMC/PBMC_Input/Artyomov_PBMC/Subset_RDS/conventional_cd8_rna.rds")
gd_t_cells <- readRDS(file = "/R/R_PBMC/PBMC_Input/Artyomov_PBMC/Subset_RDS/gd_t_cells_rna.rds")
mait_cells <- readRDS(file = "/R/R_PBMC/PBMC_Input/Artyomov_PBMC/Subset_RDS/mait_cells_rna.rds")
myeloid_cells <- readRDS(file = "/R/R_PBMC/PBMC_Input/Artyomov_PBMC/Subset_RDS/myeloid_cells_rna.rds")
nk_cells <- readRDS(file = "/R/R_PBMC/PBMC_Input/Artyomov_PBMC/Subset_RDS/nk_cells_rna.rds")
progenitor_cells <- readRDS(file = "/R/R_PBMC/PBMC_Input/Artyomov_PBMC/Subset_RDS/progenitor_cells_rna.rds")

# Initialize Seurat Objects for each
b_cells <- CreateSeuratObject(counts = b_cells, project = "b_cells", min.cells = 3, min.features = 200)
cd4_helper_memory_cells <- CreateSeuratObject(counts = cd4_helper_memory_cells, project = "cd4_helper_memory_cells", min.cells = 3, min.features = 200)
cd4_cells <- CreateSeuratObject(counts = cd4_cells, project = "cd4_cells", min.cells = 3, min.features = 200)
conventional_cd8_cells <- CreateSeuratObject(counts = conventional_cd8_cells, project = "conventional_cd8_cells", min.cells = 3, min.features = 200)
gd_t_cells <- CreateSeuratObject(counts = gd_t_cells, project = "gd_t_cells", min.cells = 3, min.features = 200)
mait_cells <- CreateSeuratObject(counts = mait_cells, project = "mait_cells", min.cells = 3, min.features = 200)
myeloid_cells <- CreateSeuratObject(counts = myeloid_cells, project = "myeloid_cells", min.cells = 3, min.features = 200)
nk_cells <- CreateSeuratObject(counts = nk_cells, project = "nk_cells", min.cells = 3, min.features = 200)
progenitor_cells <- CreateSeuratObject(counts = progenitor_cells, project = "progenitor_cells", min.cells = 3, min.features = 200)



###########
# b_cells #
###########

# Load
b_cells <- readRDS(file = "/R/R_PBMC/PBMC_Input/Artyomov_PBMC/Subset_RDS/b_cells_rna.rds")

# Create Seurat Object
b_cells <- CreateSeuratObject(counts = b_cells, project = "b_cells", min.cells = 3, min.features = 200)

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
b_cells[["percent.mt"]] <- PercentageFeatureSet(b_cells, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(b_cells, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(b_cells, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(b_cells, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
# Skipped because data already filtered
#b_cells <- subset(b_cells, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
b_cells <- NormalizeData(b_cells, normalization.method = "LogNormalize", scale.factor = 10000)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Find Variable Features
b_cells <- FindVariableFeatures(b_cells, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(b_cells)
b_cells <- ScaleData(b_cells, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
b_cells <- RunPCA(b_cells, features = VariableFeatures(object = b_cells))

# Examine and visualize PCA results a few different ways
print(b_cells[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(b_cells, dims = 1:2, reduction = "pca")
DimPlot(b_cells, reduction = "pca")

# Visualize PCs
ElbowPlot(b_cells)

# Choose dimensions
b_cells <- FindNeighbors(b_cells, dims = 1:15)
b_cells <- FindClusters(b_cells, resolution = 0.12)

# Umap clustering
b_cells <- RunUMAP(b_cells, dims = 1:15)
DimPlot(b_cells, reduction = "umap")

# Identify Samples
DimPlot(b_cells, reduction = "umap", group.by = "orig.ident")

# Save
saveRDS(b_cells, file = "/R/R_PBMC/PBMC_RDS/RDS_Total/b_cells.rds")


# Clear global environment
rm(b_cells)
gc()



###########################
# cd4_helper_memory_cells #
###########################

# Load
cd4_helper_memory_cells <- readRDS(file = "/R/R_PBMC/PBMC_Input/Artyomov_PBMC/Subset_RDS/cd4_helper_memory_rna.rds")

# Create Seurat Object
cd4_helper_memory_cells <- CreateSeuratObject(counts = cd4_helper_memory_cells, project = "cd4_helper_memory_cells", min.cells = 3, min.features = 200)



# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
cd4_helper_memory_cells[["percent.mt"]] <- PercentageFeatureSet(cd4_helper_memory_cells, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(cd4_helper_memory_cells, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(cd4_helper_memory_cells, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(cd4_helper_memory_cells, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
# Skipped because data already filtered
#cd4_helper_memory_cells <- subset(cd4_helper_memory_cells, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
cd4_helper_memory_cells <- NormalizeData(cd4_helper_memory_cells, normalization.method = "LogNormalize", scale.factor = 10000)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Find Variable Features
cd4_helper_memory_cells <- FindVariableFeatures(cd4_helper_memory_cells, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(cd4_helper_memory_cells)
cd4_helper_memory_cells <- ScaleData(cd4_helper_memory_cells, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
cd4_helper_memory_cells <- RunPCA(cd4_helper_memory_cells, features = VariableFeatures(object = cd4_helper_memory_cells))

# Examine and visualize PCA results a few different ways
print(cd4_helper_memory_cells[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(cd4_helper_memory_cells, dims = 1:2, reduction = "pca")
DimPlot(cd4_helper_memory_cells, reduction = "pca")

# Visualize PCs
ElbowPlot(cd4_helper_memory_cells)

# Choose dimensions
cd4_helper_memory_cells <- FindNeighbors(cd4_helper_memory_cells, dims = 1:15)
cd4_helper_memory_cells <- FindClusters(cd4_helper_memory_cells, resolution = 0.12)

# Umap clustering
cd4_helper_memory_cells <- RunUMAP(cd4_helper_memory_cells, dims = 1:15)
DimPlot(cd4_helper_memory_cells, reduction = "umap")

# Identify Samples
DimPlot(cd4_helper_memory_cells, reduction = "umap", group.by = "orig.ident")

# Save
saveRDS(cd4_helper_memory_cells, file = "/R/R_PBMC/PBMC_RDS/RDS_Total/cd4_helper_memory_cells.rds")

# Clear global environment
rm(cd4_helper_memory_cells)
gc()


#############
# cd4_cells #
#############

# Load
cd4_cells <- readRDS(file = "/R/R_PBMC/PBMC_Input/Artyomov_PBMC/Subset_RDS/cd4_rna.rds")


# Create Seurat Object
cd4_cells <- CreateSeuratObject(counts = cd4_cells, project = "cd4_cells", min.cells = 3, min.features = 200)



# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
cd4_cells[["percent.mt"]] <- PercentageFeatureSet(cd4_cells, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(cd4_cells, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(cd4_cells, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(cd4_cells, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
# Skipped because data already filtered
#cd4_cells <- subset(cd4_cells, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
cd4_cells <- NormalizeData(cd4_cells, normalization.method = "LogNormalize", scale.factor = 10000)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Find Variable Features
cd4_cells <- FindVariableFeatures(cd4_cells, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(cd4_cells)
cd4_cells <- ScaleData(cd4_cells, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
cd4_cells <- RunPCA(cd4_cells, features = VariableFeatures(object = cd4_cells))

# Examine and visualize PCA results a few different ways
print(cd4_cells[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(cd4_cells, dims = 1:2, reduction = "pca")
DimPlot(cd4_cells, reduction = "pca")

# Visualize PCs
ElbowPlot(cd4_cells)

# Choose dimensions
cd4_cells <- FindNeighbors(cd4_cells, dims = 1:15)
cd4_cells <- FindClusters(cd4_cells, resolution = 0.12)

# Umap clustering
cd4_cells <- RunUMAP(cd4_cells, dims = 1:15)
DimPlot(cd4_cells, reduction = "umap")

# Identify Samples
DimPlot(cd4_cells, reduction = "umap", group.by = "orig.ident")

# Save
saveRDS(cd4_cells, file = "/R/R_PBMC/PBMC_RDS/RDS_Total/cd4_cells.rds")

# Clear global environment
rm(cd4_cells)
gc()



##########################
# conventional_cd8_cells #
##########################

# Load
conventional_cd8_cells <- readRDS(file = "/R/R_PBMC/PBMC_Input/Artyomov_PBMC/Subset_RDS/conventional_cd8_rna.rds")

# Create Seurat Object
conventional_cd8_cells <- CreateSeuratObject(counts = conventional_cd8_cells, project = "conventional_cd8_cells", min.cells = 3, min.features = 200)



# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
conventional_cd8_cells[["percent.mt"]] <- PercentageFeatureSet(conventional_cd8_cells, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(conventional_cd8_cells, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(conventional_cd8_cells, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(conventional_cd8_cells, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
# Skipped because data already filtered
#conventional_cd8_cells <- subset(conventional_cd8_cells, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
conventional_cd8_cells <- NormalizeData(conventional_cd8_cells, normalization.method = "LogNormalize", scale.factor = 10000)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Find Variable Features
conventional_cd8_cells <- FindVariableFeatures(conventional_cd8_cells, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(conventional_cd8_cells)
conventional_cd8_cells <- ScaleData(conventional_cd8_cells, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
conventional_cd8_cells <- RunPCA(conventional_cd8_cells, features = VariableFeatures(object = conventional_cd8_cells))

# Examine and visualize PCA results a few different ways
print(conventional_cd8_cells[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(conventional_cd8_cells, dims = 1:2, reduction = "pca")
DimPlot(conventional_cd8_cells, reduction = "pca")

# Visualize PCs
ElbowPlot(conventional_cd8_cells)

# Choose dimensions
conventional_cd8_cells <- FindNeighbors(conventional_cd8_cells, dims = 1:15)
conventional_cd8_cells <- FindClusters(conventional_cd8_cells, resolution = 0.12)

# Umap clustering
conventional_cd8_cells <- RunUMAP(conventional_cd8_cells, dims = 1:15)
DimPlot(conventional_cd8_cells, reduction = "umap")

# Identify Samples
DimPlot(conventional_cd8_cells, reduction = "umap", group.by = "orig.ident")

# Save
saveRDS(conventional_cd8_cells, file = "/R/R_PBMC/PBMC_RDS/RDS_Total/conventional_cd8_cells.rds")

# Clear global environment
rm(conventional_cd8_cells)
gc()



##############
# gd_t_cells #
##############

# Load
gd_t_cells <- readRDS(file = "/R/R_PBMC/PBMC_Input/Artyomov_PBMC/Subset_RDS/gd_t_cells_rna.rds")

# Create Seurat Object
gd_t_cells <- CreateSeuratObject(counts = gd_t_cells, project = "gd_t_cells", min.cells = 3, min.features = 200)



# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
gd_t_cells[["percent.mt"]] <- PercentageFeatureSet(gd_t_cells, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(gd_t_cells, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(gd_t_cells, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(gd_t_cells, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
# Skipped because data already filtered
#gd_t_cells <- subset(gd_t_cells, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
gd_t_cells <- NormalizeData(gd_t_cells, normalization.method = "LogNormalize", scale.factor = 10000)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Find Variable Features
gd_t_cells <- FindVariableFeatures(gd_t_cells, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(gd_t_cells)
gd_t_cells <- ScaleData(gd_t_cells, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
gd_t_cells <- RunPCA(gd_t_cells, features = VariableFeatures(object = gd_t_cells))

# Examine and visualize PCA results a few different ways
print(gd_t_cells[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(gd_t_cells, dims = 1:2, reduction = "pca")
DimPlot(gd_t_cells, reduction = "pca")

# Visualize PCs
ElbowPlot(gd_t_cells)

# Choose dimensions
gd_t_cells <- FindNeighbors(gd_t_cells, dims = 1:15)
gd_t_cells <- FindClusters(gd_t_cells, resolution = 0.12)

# Umap clustering
gd_t_cells <- RunUMAP(gd_t_cells, dims = 1:15)
DimPlot(gd_t_cells, reduction = "umap")

# Identify Samples
DimPlot(gd_t_cells, reduction = "umap", group.by = "orig.ident")

# Save
saveRDS(gd_t_cells, file = "/R/R_PBMC/PBMC_RDS/RDS_Total/gd_t_cells.rds")

# Clear global environment
rm(gd_t_cells)
gc()



##############
# mait_cells #
##############

# Load
mait_cells <- readRDS(file = "/R/R_PBMC/PBMC_Input/Artyomov_PBMC/Subset_RDS/mait_cells_rna.rds")

# Create Seurat Object
mait_cells <- CreateSeuratObject(counts = mait_cells, project = "mait_cells", min.cells = 3, min.features = 200)



# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
mait_cells[["percent.mt"]] <- PercentageFeatureSet(mait_cells, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(mait_cells, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(mait_cells, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(mait_cells, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
# Skipped because data already filtered
#mait_cells <- subset(mait_cells, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
mait_cells <- NormalizeData(mait_cells, normalization.method = "LogNormalize", scale.factor = 10000)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Find Variable Features
mait_cells <- FindVariableFeatures(mait_cells, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(mait_cells)
mait_cells <- ScaleData(mait_cells, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
mait_cells <- RunPCA(mait_cells, features = VariableFeatures(object = mait_cells))

# Examine and visualize PCA results a few different ways
print(mait_cells[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(mait_cells, dims = 1:2, reduction = "pca")
DimPlot(mait_cells, reduction = "pca")

# Visualize PCs
ElbowPlot(mait_cells)

# Choose dimensions
mait_cells <- FindNeighbors(mait_cells, dims = 1:15)
mait_cells <- FindClusters(mait_cells, resolution = 0.12)

# Umap clustering
mait_cells <- RunUMAP(mait_cells, dims = 1:15)
DimPlot(mait_cells, reduction = "umap")

# Identify Samples
DimPlot(mait_cells, reduction = "umap", group.by = "orig.ident")

# Save
saveRDS(mait_cells, file = "/R/R_PBMC/PBMC_RDS/RDS_Total/mait_cells.rds")

# Clear global environment
rm(mait_cells)
gc()



#################
# myeloid_cells #
#################

# Load
myeloid_cells <- readRDS(file = "/R/R_PBMC/PBMC_Input/Artyomov_PBMC/Subset_RDS/myeloid_cells_rna.rds")

# Create Seurat Object
myeloid_cells <- CreateSeuratObject(counts = myeloid_cells, project = "myeloid_cells", min.cells = 3, min.features = 200)



# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
myeloid_cells[["percent.mt"]] <- PercentageFeatureSet(myeloid_cells, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(myeloid_cells, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(myeloid_cells, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(myeloid_cells, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
# Skipped because data already filtered
#myeloid_cells <- subset(myeloid_cells, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
myeloid_cells <- NormalizeData(myeloid_cells, normalization.method = "LogNormalize", scale.factor = 10000)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Find Variable Features
myeloid_cells <- FindVariableFeatures(myeloid_cells, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(myeloid_cells)
myeloid_cells <- ScaleData(myeloid_cells, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
myeloid_cells <- RunPCA(myeloid_cells, features = VariableFeatures(object = myeloid_cells))

# Examine and visualize PCA results a few different ways
print(myeloid_cells[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(myeloid_cells, dims = 1:2, reduction = "pca")
DimPlot(myeloid_cells, reduction = "pca")

# Visualize PCs
ElbowPlot(myeloid_cells)

# Choose dimensions
myeloid_cells <- FindNeighbors(myeloid_cells, dims = 1:15)
myeloid_cells <- FindClusters(myeloid_cells, resolution = 0.12)

# Umap clustering
myeloid_cells <- RunUMAP(myeloid_cells, dims = 1:15)
DimPlot(myeloid_cells, reduction = "umap")

# Identify Samples
DimPlot(myeloid_cells, reduction = "umap", group.by = "orig.ident")

# Save
saveRDS(myeloid_cells, file = "/R/R_PBMC/PBMC_RDS/RDS_Total/myeloid_cells.rds")

# Clear global environment
rm(myeloid_cells)
gc()


############
# nk_cells #
############

# Load
nk_cells <- readRDS(file = "/R/R_PBMC/PBMC_Input/Artyomov_PBMC/Subset_RDS/nk_cells_rna.rds")

# Create Seurat Object
nk_cells <- CreateSeuratObject(counts = nk_cells, project = "nk_cells", min.cells = 3, min.features = 200)



# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
nk_cells[["percent.mt"]] <- PercentageFeatureSet(nk_cells, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(nk_cells, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(nk_cells, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(nk_cells, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
# Skipped because data already filtered
#nk_cells <- subset(nk_cells, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
nk_cells <- NormalizeData(nk_cells, normalization.method = "LogNormalize", scale.factor = 10000)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Find Variable Features
nk_cells <- FindVariableFeatures(nk_cells, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(nk_cells)
nk_cells <- ScaleData(nk_cells, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
nk_cells <- RunPCA(nk_cells, features = VariableFeatures(object = nk_cells))

# Examine and visualize PCA results a few different ways
print(nk_cells[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(nk_cells, dims = 1:2, reduction = "pca")
DimPlot(nk_cells, reduction = "pca")

# Visualize PCs
ElbowPlot(nk_cells)

# Choose dimensions
nk_cells <- FindNeighbors(nk_cells, dims = 1:15)
nk_cells <- FindClusters(nk_cells, resolution = 0.12)

# Umap clustering
nk_cells <- RunUMAP(nk_cells, dims = 1:15)
DimPlot(nk_cells, reduction = "umap")

# Identify Samples
DimPlot(nk_cells, reduction = "umap", group.by = "orig.ident")

# Save
saveRDS(nk_cells, file = "/R/R_PBMC/PBMC_RDS/RDS_Total/nk_cells.rds")

# Clear global environment
rm(nk_cells)
gc()



####################
# progenitor_cells #
####################

# Load
progenitor_cells <- readRDS(file = "/R/R_PBMC/PBMC_Input/Artyomov_PBMC/Subset_RDS/progenitor_cells_rna.rds")

# Create Seurat Object
progenitor_cells <- CreateSeuratObject(counts = progenitor_cells, project = "progenitor_cells", min.cells = 3, min.features = 200)



# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
progenitor_cells[["percent.mt"]] <- PercentageFeatureSet(progenitor_cells, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(progenitor_cells, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(progenitor_cells, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(progenitor_cells, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data; mitochondrial counts higher than normal; losing many cells
# Skipped because data already filtered
#progenitor_cells <- subset(progenitor_cells, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalize data
progenitor_cells <- NormalizeData(progenitor_cells, normalization.method = "LogNormalize", scale.factor = 10000)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Find Variable Features
progenitor_cells <- FindVariableFeatures(progenitor_cells, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(progenitor_cells)
progenitor_cells <- ScaleData(progenitor_cells, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
progenitor_cells <- RunPCA(progenitor_cells, features = VariableFeatures(object = progenitor_cells))

# Examine and visualize PCA results a few different ways
print(progenitor_cells[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(progenitor_cells, dims = 1:2, reduction = "pca")
DimPlot(progenitor_cells, reduction = "pca")

# Visualize PCs
ElbowPlot(progenitor_cells)

# Choose dimensions
progenitor_cells <- FindNeighbors(progenitor_cells, dims = 1:15)
progenitor_cells <- FindClusters(progenitor_cells, resolution = 0.12)

# Umap clustering
progenitor_cells <- RunUMAP(progenitor_cells, dims = 1:15)
DimPlot(progenitor_cells, reduction = "umap")

# Identify Samples
DimPlot(progenitor_cells, reduction = "umap", group.by = "orig.ident")

# Save
saveRDS(progenitor_cells, file = "/R/R_PBMC/PBMC_RDS/RDS_Total/progenitor_cells.rds")

# Clear global environment
rm(progenitor_cells)
gc()


progenitor_cells <- readRDS(file = "/R/R_PBMC/PBMC_RDS/RDS_Total/progenitor_cells.rds") 




###################################################################################
# Step 2: Downsample Leukocyte Subsets and then Merge into Combined Seurat Object #
###################################################################################

############################
# Load Immune Cell Subsets #
############################
b_cells <- readRDS(file = "/R/R_PBMC/PBMC_RDS/RDS_Total/b_cells.rds") 
cd4_helper_memory_cells <- readRDS(file = "/R/R_PBMC/PBMC_RDS/RDS_Total/cd4_helper_memory_cells.rds") 
cd4_cells <- readRDS(file = "/R/R_PBMC/PBMC_RDS/RDS_Total/cd4_cells.rds") 
conventional_cd8_cells <- readRDS(file = "/R/R_PBMC/PBMC_RDS/RDS_Total/conventional_cd8_cells.rds") 
gd_t_cells <- readRDS(file = "/R/R_PBMC/PBMC_RDS/RDS_Total/gd_t_cells.rds") 
mait_cells <- readRDS(file = "/R/R_PBMC/PBMC_RDS/RDS_Total/mait_cells.rds") 
myeloid_cells <- readRDS(file = "/R/R_PBMC/PBMC_RDS/RDS_Total/myeloid_cells.rds") 
nk_cells <- readRDS(file = "/R/R_PBMC/PBMC_RDS/RDS_Total/nk_cells.rds") 
progenitor_cells <- readRDS(file = "/R/R_PBMC/PBMC_RDS/RDS_Total/progenitor_cells.rds") 


# Count Total Cells for each subset
b_cells_Total <- Idents(b_cells)
cd4_helper_memory_cells_Total <- Idents(cd4_helper_memory_cells)
cd4_cells_Total <- Idents(cd4_cells)
conventional_cd8_cells_Total <- Idents(conventional_cd8_cells)
gd_t_cells_Total <- Idents(gd_t_cells)
mait_cells_Total <- Idents(mait_cells)
myeloid_cells_Total <- Idents(myeloid_cells)
nk_cells_Total <- Idents(nk_cells)
progenitor_cells_Total <- Idents(progenitor_cells)


# Counts
# b_cells_Total = 71,614
# cd4_helper_memory_cells_Total = 352,931
# cd4_cells_Total = 901,152
# conventional_cd8_cells_Total = 313,343
# gd_t_cells_Total = 60,325
# mait_cells_Total = 24,245
# myeloid_cells_Total = 336,935
# nk_cells_Total = 205,469
# progenitor_cells_Total = 1,794
# GRAND TOTAL = 2,267,772

# Remove counts
rm(b_cells_Total)
rm(cd4_helper_memory_cells_Total)
rm(cd4_cells_Total)
rm(conventional_cd8_cells_Total)
rm(gd_t_cells_Total)
rm(mait_cells_Total)
rm(myeloid_cells_Total)
rm(nk_cells_Total)
rm(progenitor_cells_Total)


# Downsample a proportion of each object, to make a representative subsample using these numbers:
# Subsample to achieve 500,000 cells (~22%)
# b_cells = 15938
# cd4_helper_memory_cells = 77744
# cd4_cells = 198353
# conventional_cd8_cells = 69035
# gd_t_cells = 13371
# mait_cells = 5433
# myeloid_cells = 74325
# nk_cells = 45403
# progenitor_cells = 398
# GRAND TOTAL = 200,000

# Set seed for reproducibility
set.seed(123)

# Subsetting by Downsampling
b_cells <- b_cells[, sample(colnames(b_cells), size = 15938, replace=FALSE)]
cd4_helper_memory_cells <- cd4_helper_memory_cells[, sample(colnames(cd4_helper_memory_cells), size = 77744, replace=FALSE)]
cd4_cells <- cd4_cells[, sample(colnames(cd4_cells), size = 198353, replace=FALSE)]
conventional_cd8_cells <- conventional_cd8_cells[, sample(colnames(conventional_cd8_cells), size = 69035, replace=FALSE)]
gd_t_cells <- gd_t_cells[, sample(colnames(gd_t_cells), size = 13371, replace=FALSE)]
mait_cells <- mait_cells[, sample(colnames(mait_cells), size = 5433, replace=FALSE)]
myeloid_cells <- myeloid_cells[, sample(colnames(myeloid_cells), size = 74325, replace=FALSE)]
nk_cells <- nk_cells[, sample(colnames(nk_cells), size = 45403, replace=FALSE)]
progenitor_cells <- progenitor_cells[, sample(colnames(progenitor_cells), size = 398, replace=FALSE)]

# Merge Datasets
Artyomov_PBMC_500k_SKETCHED <- merge(x = b_cells, y = c(cd4_helper_memory_cells, cd4_cells, conventional_cd8_cells,
                                          gd_t_cells, mait_cells, myeloid_cells, nk_cells, progenitor_cells))

# Remove Objects
rm(b_cells)
rm(cd4_helper_memory_cells)
rm(cd4_cells)
rm(conventional_cd8_cells)
rm(gd_t_cells)
rm(mait_cells)
rm(myeloid_cells)
rm(nk_cells)
rm(progenitor_cells)
gc()


# Join Layers
Artyomov_PBMC_500k_Sampled <- JoinLayers(Artyomov_PBMC_500k_Sampled)

# Normalize data
Artyomov_PBMC_500k_Sampled <- NormalizeData(Artyomov_PBMC_500k_Sampled, normalization.method = "LogNormalize", scale.factor = 10000)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Find Variable Features
Artyomov_PBMC_500k_Sampled <- FindVariableFeatures(Artyomov_PBMC_500k_Sampled, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Artyomov_PBMC_500k_Sampled)
Artyomov_PBMC_500k_Sampled <- ScaleData(Artyomov_PBMC_500k_Sampled, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Artyomov_PBMC_500k_Sampled <- RunPCA(Artyomov_PBMC_500k_Sampled, features = VariableFeatures(object = Artyomov_PBMC_500k_Sampled))

# Examine and visualize PCA results a few different ways
print(Artyomov_PBMC_500k_Sampled[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Artyomov_PBMC_500k_Sampled, dims = 1:2, reduction = "pca")
DimPlot(Artyomov_PBMC_500k_Sampled, reduction = "pca")

# Visualize PCs
ElbowPlot(Artyomov_PBMC_500k_Sampled)

# Choose dimensions
Artyomov_PBMC_500k_Sampled <- FindNeighbors(Artyomov_PBMC_500k_Sampled, dims = 1:15)
Artyomov_PBMC_500k_Sampled <- FindClusters(Artyomov_PBMC_500k_Sampled, resolution = 0.12)

# Umap clustering
Artyomov_PBMC_500k_Sampled <- RunUMAP(Artyomov_PBMC_500k_Sampled, dims = 1:15)
DimPlot(Artyomov_PBMC_500k_Sampled, reduction = "umap", raster = FALSE)

# Identify Samples
DimPlot(Artyomov_PBMC_500k_Sampled, reduction = "umap", group.by = "orig.ident")

# Save
saveRDS(Artyomov_PBMC_500k_Sampled, file = "/R/R_PBMC/PBMC_RDS/RDS_Total/Artyomov_PBMC_500k_Sampled.rds")

# Clear global environment
rm(Artyomov_PBMC_500k_Sampled)
gc()



##########################
# Generate Example Plots #
##########################

#Artyomov_PBMC_500k_Sampled <- readRDS(file = "/R/R_PBMC/PBMC_RDS/RDS_Total/Artyomov_PBMC_500k_Sampled.rds")


library(ggplot2)
FeaturePlot(Artyomov_PBMC_500k_Sampled, c("KRT18"),  raster = FALSE, cols = c("grey90", "firebrick"))
ggsave("/R/R_PBMC/PBMC_Output/Artyomov_PBMC_500k_Sampled_KRT18.tiff", plot = last_plot(), device = "tiff", 
       scale = 1, width = 16, height = 10,
       dpi = 200, limitsize = TRUE)

library(ggplot2)
FeaturePlot(Artyomov_PBMC_500k_Sampled, c("KRT19"),  raster = FALSE, cols = c("grey90", "firebrick"))
ggsave("/R/R_PBMC/PBMC_Output/Artyomov_PBMC_500k_Sampled_KRT19.tiff", plot = last_plot(), device = "tiff", 
       scale = 1, width = 16, height = 10,
       dpi = 200, limitsize = TRUE)

DimPlot(Artyomov_PBMC_500k_Sampled, raster = FALSE)
ggsave("/R/R_PBMC/PBMC_Output/Artyomov_PBMC_500k_Sampled_Clusters.tiff", plot = last_plot(), device = "tiff", 
       scale = 1, width = 16, height = 10,
       dpi = 200, limitsize = TRUE)

DimPlot(Artyomov_PBMC_500k_Sampled, raster = FALSE, group.by = "orig.ident")
ggsave("/R/R_PBMC/PBMC_Output/Artyomov_PBMC_500k_Sampled_Clusters_Ident.tiff", plot = last_plot(), device = "tiff", 
       scale = 1, width = 16, height = 10,
       dpi = 200, limitsize = TRUE)



######################################################################################
# Step 3: Extract Expression Matrices for Mammary Keratin Genes per Leukocyte Subset #
######################################################################################


# Extract Expression Matrix
b_cells_matrix <- b_cells[["RNA"]]$counts
b_cells_matrix <- as.matrix(b_cells_matrix, 'sparseMatrix')
b_cells_matrix <- subset(b_cells_matrix, rownames(b_cells_matrix) %in% c("KRT14", "KRT18", "KRT19"))
b_cells_matrix <-t(b_cells_matrix)
write.csv(b_cells_matrix, file = "/R/R_PBMC/PBMC_Expression_Matrices/PBMC_Expression_Matrices_Output/b_cells_KRT_Expression.csv")
rm(b_cells_matrix)
gc()

# Extract Expression Matrix
cd4_helper_memory_cells_matrix <- cd4_helper_memory_cells[["RNA"]]$counts
cd4_helper_memory_cells_matrix <- as.matrix(cd4_helper_memory_cells_matrix, 'sparseMatrix')
cd4_helper_memory_cells_matrix <- subset(cd4_helper_memory_cells_matrix, rownames(cd4_helper_memory_cells_matrix) %in% c("KRT14", "KRT18", "KRT19"))
cd4_helper_memory_cells_matrix <-t(cd4_helper_memory_cells_matrix)
write.csv(cd4_helper_memory_cells_matrix, file = "/R/R_PBMC/PBMC_Expression_Matrices/PBMC_Expression_Matrices_Output/cd4_helper_memory_cells_KRT_Expression.csv")
rm(cd4_helper_memory_cells_matrix)
gc()

# Extract Expression Matrix
cd4_cells_matrix <- cd4_cells[["RNA"]]$counts
cd4_cells_matrix <- as.matrix(cd4_cells_matrix, 'sparseMatrix')
cd4_cells_matrix <- subset(cd4_cells_matrix, rownames(cd4_cells_matrix) %in% c("KRT14", "KRT18", "KRT19"))
cd4_cells_matrix <-t(cd4_cells_matrix)
write.csv(cd4_cells_matrix, file = "/R/R_PBMC/PBMC_Expression_Matrices/PBMC_Expression_Matrices_Output/cd4_cells_KRT_Expression.csv")
rm(cd4_cells_matrix)
gc()

# Extract Expression Matrix
conventional_cd8_cells_matrix <- conventional_cd8_cells[["RNA"]]$counts
conventional_cd8_cells_matrix <- as.matrix(conventional_cd8_cells_matrix, 'sparseMatrix')
conventional_cd8_cells_matrix <- subset(conventional_cd8_cells_matrix, rownames(conventional_cd8_cells_matrix) %in% c("KRT14", "KRT18", "KRT19"))
conventional_cd8_cells_matrix <-t(conventional_cd8_cells_matrix)
write.csv(conventional_cd8_cells_matrix, file = "/R/R_PBMC/PBMC_Expression_Matrices/PBMC_Expression_Matrices_Output/conventional_cd8_cells_KRT_Expression.csv")
rm(conventional_cd8_cells_matrix)
gc()

# Extract Expression Matrix
gd_t_cells_matrix <- gd_t_cells[["RNA"]]$counts
gd_t_cells_matrix <- as.matrix(gd_t_cells_matrix, 'sparseMatrix')
gd_t_cells_matrix <- subset(gd_t_cells_matrix, rownames(gd_t_cells_matrix) %in% c("KRT14", "KRT18", "KRT19"))
gd_t_cells_matrix <-t(gd_t_cells_matrix)
write.csv(gd_t_cells_matrix, file = "/R/R_PBMC/PBMC_Expression_Matrices/PBMC_Expression_Matrices_Output/gd_t_cells_KRT_Expression.csv")
rm(gd_t_cells_matrix)
gc()

# Extract Expression Matrix
mait_cells_matrix <- mait_cells[["RNA"]]$counts
mait_cells_matrix <- as.matrix(mait_cells_matrix, 'sparseMatrix')
mait_cells_matrix <- subset(mait_cells_matrix, rownames(mait_cells_matrix) %in% c("KRT14", "KRT18", "KRT19"))
mait_cells_matrix <-t(mait_cells_matrix)
write.csv(mait_cells_matrix, file = "/R/R_PBMC/PBMC_Expression_Matrices/PBMC_Expression_Matrices_Output/mait_cells_KRT_Expression.csv")
rm(mait_cells_matrix)
gc()

# Extract Expression Matrix
myeloid_cells_matrix <- myeloid_cells[["RNA"]]$counts
myeloid_cells_matrix <- as.matrix(myeloid_cells_matrix, 'sparseMatrix')
myeloid_cells_matrix <- subset(myeloid_cells_matrix, rownames(myeloid_cells_matrix) %in% c("KRT14", "KRT18", "KRT19"))
myeloid_cells_matrix <-t(myeloid_cells_matrix)
write.csv(myeloid_cells_matrix, file = "/R/R_PBMC/PBMC_Expression_Matrices/PBMC_Expression_Matrices_Output/myeloid_cells_KRT_Expression.csv")
rm(myeloid_cells_matrix)
gc()

# Extract Expression Matrix
nk_cells_matrix <- nk_cells[["RNA"]]$counts
nk_cells_matrix <- as.matrix(nk_cells_matrix, 'sparseMatrix')
nk_cells_matrix <- subset(nk_cells_matrix, rownames(nk_cells_matrix) %in% c("KRT14", "KRT18", "KRT19"))
nk_cells_matrix <-t(nk_cells_matrix)
write.csv(nk_cells_matrix, file = "/R/R_PBMC/PBMC_Expression_Matrices/PBMC_Expression_Matrices_Output/nk_cells_KRT_Expression.csv")
rm(nk_cells_matrix)
gc()

# Extract Expression Matrix
progenitor_cells_matrix <- progenitor_cells[["RNA"]]$counts
progenitor_cells_matrix <- as.matrix(progenitor_cells_matrix, 'sparseMatrix')
progenitor_cells_matrix <- subset(progenitor_cells_matrix, rownames(progenitor_cells_matrix) %in% c("KRT14", "KRT18", "KRT19"))
progenitor_cells_matrix <-t(progenitor_cells_matrix)
write.csv(progenitor_cells_matrix, file = "/R/R_PBMC/PBMC_Expression_Matrices/PBMC_Expression_Matrices_Output/progenitor_cells_KRT_Expression.csv")
rm(progenitor_cells_matrix)
gc()


