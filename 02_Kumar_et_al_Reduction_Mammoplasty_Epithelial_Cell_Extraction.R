#####################################################################################################################
#                   Navin (Kumar et al) Dataset Reduction Mammoplasty Analysis Steps                                #
#####################################################################################################################
# Step 3: Extract Epithelial Cells from Mammoplasty Samples                                                         #
# Step 4: Merge Epithelial Cells into Combined Seurat Object                                                        #
#####################################################################################################################



#################################################################
# Step 3: Extract Epithelial Cells from Mammoplasty Samples     # 
#################################################################
# Load Libraries
library(Seurat)
library(patchwork)
library(dplyr)

############################ 
# Navin_hbca_c14_Singlets  #
############################ 
# Load Object
Navin_hbca_c14_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c14_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c14_panCK <- subset(x = Navin_hbca_c14_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c14_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c14_panCK <- FindVariableFeatures(Navin_hbca_c14_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c14_panCK)
Navin_hbca_c14_panCK <- ScaleData(Navin_hbca_c14_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c14_panCK <- RunPCA(Navin_hbca_c14_panCK, features = VariableFeatures(object = Navin_hbca_c14_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c14_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c14_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c14_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c14_panCK, ndims = 60)

Navin_hbca_c14_panCK <- FindNeighbors(Navin_hbca_c14_panCK, dims = 1:40)
Navin_hbca_c14_panCK <- FindClusters(Navin_hbca_c14_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c14_panCK <- RunUMAP(Navin_hbca_c14_panCK, dims = 1:40)
DimPlot(Navin_hbca_c14_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c14_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c14_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c14_panCK.rds")

# Remove Object
rm(Navin_hbca_c14_panCK)
gc()



############################ 
# Navin_hbca_c15_Singlets  #
############################ 
# Load Object
Navin_hbca_c15_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c15_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c15_panCK <- subset(x = Navin_hbca_c15_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c15_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c15_panCK <- FindVariableFeatures(Navin_hbca_c15_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c15_panCK)
Navin_hbca_c15_panCK <- ScaleData(Navin_hbca_c15_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c15_panCK <- RunPCA(Navin_hbca_c15_panCK, features = VariableFeatures(object = Navin_hbca_c15_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c15_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c15_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c15_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c15_panCK, ndims = 60)

Navin_hbca_c15_panCK <- FindNeighbors(Navin_hbca_c15_panCK, dims = 1:40)
Navin_hbca_c15_panCK <- FindClusters(Navin_hbca_c15_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c15_panCK <- RunUMAP(Navin_hbca_c15_panCK, dims = 1:40)
DimPlot(Navin_hbca_c15_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c15_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c15_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c15_panCK.rds")

# Remove Object
rm(Navin_hbca_c15_panCK)
gc()



############################ 
# Navin_hbca_c19_Singlets  #
############################ 
# Load Object
Navin_hbca_c19_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c19_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c19_panCK <- subset(x = Navin_hbca_c19_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c19_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c19_panCK <- FindVariableFeatures(Navin_hbca_c19_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c19_panCK)
Navin_hbca_c19_panCK <- ScaleData(Navin_hbca_c19_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c19_panCK <- RunPCA(Navin_hbca_c19_panCK, features = VariableFeatures(object = Navin_hbca_c19_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c19_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c19_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c19_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c19_panCK, ndims = 60)

Navin_hbca_c19_panCK <- FindNeighbors(Navin_hbca_c19_panCK, dims = 1:40)
Navin_hbca_c19_panCK <- FindClusters(Navin_hbca_c19_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c19_panCK <- RunUMAP(Navin_hbca_c19_panCK, dims = 1:40)
DimPlot(Navin_hbca_c19_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c19_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c19_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c19_panCK.rds")

# Remove Object
rm(Navin_hbca_c19_panCK)
gc()



############################ 
# Navin_hbca_c20_Singlets  #
############################ 
# Load Object
Navin_hbca_c20_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c20_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c20_panCK <- subset(x = Navin_hbca_c20_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c20_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c20_panCK <- FindVariableFeatures(Navin_hbca_c20_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c20_panCK)
Navin_hbca_c20_panCK <- ScaleData(Navin_hbca_c20_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c20_panCK <- RunPCA(Navin_hbca_c20_panCK, features = VariableFeatures(object = Navin_hbca_c20_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c20_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c20_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c20_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c20_panCK, ndims = 60)

Navin_hbca_c20_panCK <- FindNeighbors(Navin_hbca_c20_panCK, dims = 1:40)
Navin_hbca_c20_panCK <- FindClusters(Navin_hbca_c20_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c20_panCK <- RunUMAP(Navin_hbca_c20_panCK, dims = 1:40)
DimPlot(Navin_hbca_c20_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c20_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c20_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c20_panCK.rds")

# Remove Object
rm(Navin_hbca_c20_panCK)
gc()



############################ 
# Navin_hbca_c22_Singlets  #
############################ 
# Load Object
Navin_hbca_c22_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c22_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c22_panCK <- subset(x = Navin_hbca_c22_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c22_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c22_panCK <- FindVariableFeatures(Navin_hbca_c22_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c22_panCK)
Navin_hbca_c22_panCK <- ScaleData(Navin_hbca_c22_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c22_panCK <- RunPCA(Navin_hbca_c22_panCK, features = VariableFeatures(object = Navin_hbca_c22_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c22_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c22_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c22_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c22_panCK, ndims = 60)

Navin_hbca_c22_panCK <- FindNeighbors(Navin_hbca_c22_panCK, dims = 1:40)
Navin_hbca_c22_panCK <- FindClusters(Navin_hbca_c22_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c22_panCK <- RunUMAP(Navin_hbca_c22_panCK, dims = 1:40)
DimPlot(Navin_hbca_c22_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c22_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c22_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c22_panCK.rds")

# Remove Object
rm(Navin_hbca_c22_panCK)
gc()



############################ 
# Navin_hbca_c23_Singlets  #
############################ 
# Load Object
Navin_hbca_c23_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c23_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c23_panCK <- subset(x = Navin_hbca_c23_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c23_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c23_panCK <- FindVariableFeatures(Navin_hbca_c23_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c23_panCK)
Navin_hbca_c23_panCK <- ScaleData(Navin_hbca_c23_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c23_panCK <- RunPCA(Navin_hbca_c23_panCK, features = VariableFeatures(object = Navin_hbca_c23_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c23_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c23_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c23_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c23_panCK, ndims = 60)

Navin_hbca_c23_panCK <- FindNeighbors(Navin_hbca_c23_panCK, dims = 1:40)
Navin_hbca_c23_panCK <- FindClusters(Navin_hbca_c23_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c23_panCK <- RunUMAP(Navin_hbca_c23_panCK, dims = 1:40)
DimPlot(Navin_hbca_c23_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c23_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c23_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c23_panCK.rds")

# Remove Object
rm(Navin_hbca_c23_panCK)
gc()



############################ 
# Navin_hbca_c24_Singlets  #
############################ 
# Load Object
Navin_hbca_c24_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c24_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c24_panCK <- subset(x = Navin_hbca_c24_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c24_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c24_panCK <- FindVariableFeatures(Navin_hbca_c24_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c24_panCK)
Navin_hbca_c24_panCK <- ScaleData(Navin_hbca_c24_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c24_panCK <- RunPCA(Navin_hbca_c24_panCK, features = VariableFeatures(object = Navin_hbca_c24_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c24_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c24_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c24_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c24_panCK, ndims = 60)

Navin_hbca_c24_panCK <- FindNeighbors(Navin_hbca_c24_panCK, dims = 1:40)
Navin_hbca_c24_panCK <- FindClusters(Navin_hbca_c24_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c24_panCK <- RunUMAP(Navin_hbca_c24_panCK, dims = 1:40)
DimPlot(Navin_hbca_c24_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c24_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c24_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c24_panCK.rds")

# Remove Object
rm(Navin_hbca_c24_panCK)
gc()



############################ 
# Navin_hbca_c25_Singlets  #
############################ 

# Note: Not enough cells in Navin_hbca_c25_panCK to redo clustering; matrices for this sample were exported directly after subsetting

# Load Object
Navin_hbca_c25_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c25_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c25_panCK <- subset(x = Navin_hbca_c25_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c25_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
#Navin_hbca_c25_panCK <- FindVariableFeatures(Navin_hbca_c25_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
#all.genes <- rownames(Navin_hbca_c25_panCK)
#Navin_hbca_c25_panCK <- ScaleData(Navin_hbca_c25_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
#Navin_hbca_c25_panCK <- RunPCA(Navin_hbca_c25_panCK, features = VariableFeatures(object = Navin_hbca_c25_panCK))

# Examine and visualize PCA results a few different ways
#print(Navin_hbca_c25_panCK[["pca"]], dims = 1:5, nfeatures = 5)
#VizDimLoadings(Navin_hbca_c25_panCK, dims = 1:2, reduction = "pca")
#DimPlot(Navin_hbca_c25_panCK, reduction = "pca")

# Visualize PCs
#ElbowPlot(Navin_hbca_c25_panCK, ndims = 60)

#Navin_hbca_c25_panCK <- FindNeighbors(Navin_hbca_c25_panCK, dims = 1:40)
#Navin_hbca_c25_panCK <- FindClusters(Navin_hbca_c25_panCK, resolution = 0.15)

# Umap clustering
#Navin_hbca_c25_panCK <- RunUMAP(Navin_hbca_c25_panCK, dims = 1:40)
#DimPlot(Navin_hbca_c25_panCK, reduction = "umap", raster = FALSE)
#DimPlot(Navin_hbca_c25_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c25_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c25_panCK.rds")

# Remove Object
rm(Navin_hbca_c25_panCK)
gc()



############################ 
# Navin_hbca_c26_Singlets  #
############################ 
# Load Object
Navin_hbca_c26_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c26_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c26_panCK <- subset(x = Navin_hbca_c26_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c26_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c26_panCK <- FindVariableFeatures(Navin_hbca_c26_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c26_panCK)
Navin_hbca_c26_panCK <- ScaleData(Navin_hbca_c26_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c26_panCK <- RunPCA(Navin_hbca_c26_panCK, features = VariableFeatures(object = Navin_hbca_c26_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c26_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c26_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c26_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c26_panCK, ndims = 60)

Navin_hbca_c26_panCK <- FindNeighbors(Navin_hbca_c26_panCK, dims = 1:40)
Navin_hbca_c26_panCK <- FindClusters(Navin_hbca_c26_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c26_panCK <- RunUMAP(Navin_hbca_c26_panCK, dims = 1:40)
DimPlot(Navin_hbca_c26_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c26_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c26_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c26_panCK.rds")

# Remove Object
rm(Navin_hbca_c26_panCK)
gc()



############################ 
# Navin_hbca_c31_Singlets  #
############################ 
# Load Object
Navin_hbca_c31_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c31_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c31_panCK <- subset(x = Navin_hbca_c31_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c31_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c31_panCK <- FindVariableFeatures(Navin_hbca_c31_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c31_panCK)
Navin_hbca_c31_panCK <- ScaleData(Navin_hbca_c31_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c31_panCK <- RunPCA(Navin_hbca_c31_panCK, features = VariableFeatures(object = Navin_hbca_c31_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c31_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c31_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c31_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c31_panCK, ndims = 60)

Navin_hbca_c31_panCK <- FindNeighbors(Navin_hbca_c31_panCK, dims = 1:40)
Navin_hbca_c31_panCK <- FindClusters(Navin_hbca_c31_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c31_panCK <- RunUMAP(Navin_hbca_c31_panCK, dims = 1:40)
DimPlot(Navin_hbca_c31_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c31_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c31_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c31_panCK.rds")

# Remove Object
rm(Navin_hbca_c31_panCK)
gc()





############################ 
# Navin_hbca_c32_Singlets  #
############################ 
# Load Object
Navin_hbca_c32_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c32_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c32_panCK <- subset(x = Navin_hbca_c32_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c32_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c32_panCK <- FindVariableFeatures(Navin_hbca_c32_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c32_panCK)
Navin_hbca_c32_panCK <- ScaleData(Navin_hbca_c32_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c32_panCK <- RunPCA(Navin_hbca_c32_panCK, features = VariableFeatures(object = Navin_hbca_c32_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c32_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c32_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c32_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c32_panCK, ndims = 60)

Navin_hbca_c32_panCK <- FindNeighbors(Navin_hbca_c32_panCK, dims = 1:40)
Navin_hbca_c32_panCK <- FindClusters(Navin_hbca_c32_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c32_panCK <- RunUMAP(Navin_hbca_c32_panCK, dims = 1:40)
DimPlot(Navin_hbca_c32_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c32_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c32_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c32_panCK.rds")

# Remove Object
rm(Navin_hbca_c32_panCK)
gc()



############################ 
# Navin_hbca_c50_Singlets  #
############################ 
# Load Object
Navin_hbca_c50_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c50_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c50_panCK <- subset(x = Navin_hbca_c50_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c50_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c50_panCK <- FindVariableFeatures(Navin_hbca_c50_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c50_panCK)
Navin_hbca_c50_panCK <- ScaleData(Navin_hbca_c50_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c50_panCK <- RunPCA(Navin_hbca_c50_panCK, features = VariableFeatures(object = Navin_hbca_c50_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c50_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c50_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c50_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c50_panCK, ndims = 60)

Navin_hbca_c50_panCK <- FindNeighbors(Navin_hbca_c50_panCK, dims = 1:40)
Navin_hbca_c50_panCK <- FindClusters(Navin_hbca_c50_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c50_panCK <- RunUMAP(Navin_hbca_c50_panCK, dims = 1:40)
DimPlot(Navin_hbca_c50_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c50_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c50_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c50_panCK.rds")

# Remove Object
rm(Navin_hbca_c50_panCK)
gc()



############################ 
# Navin_hbca_c51_Singlets  #
############################ 
# Load Object
Navin_hbca_c51_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c51_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c51_panCK <- subset(x = Navin_hbca_c51_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c51_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c51_panCK <- FindVariableFeatures(Navin_hbca_c51_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c51_panCK)
Navin_hbca_c51_panCK <- ScaleData(Navin_hbca_c51_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c51_panCK <- RunPCA(Navin_hbca_c51_panCK, features = VariableFeatures(object = Navin_hbca_c51_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c51_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c51_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c51_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c51_panCK, ndims = 60)

Navin_hbca_c51_panCK <- FindNeighbors(Navin_hbca_c51_panCK, dims = 1:40)
Navin_hbca_c51_panCK <- FindClusters(Navin_hbca_c51_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c51_panCK <- RunUMAP(Navin_hbca_c51_panCK, dims = 1:40)
DimPlot(Navin_hbca_c51_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c51_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c51_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c51_panCK.rds")

# Remove Object
rm(Navin_hbca_c51_panCK)
gc()



############################ 
# Navin_hbca_c52_Singlets  #
############################ 
# Load Object
Navin_hbca_c52_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c52_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c52_panCK <- subset(x = Navin_hbca_c52_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c52_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c52_panCK <- FindVariableFeatures(Navin_hbca_c52_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c52_panCK)
Navin_hbca_c52_panCK <- ScaleData(Navin_hbca_c52_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c52_panCK <- RunPCA(Navin_hbca_c52_panCK, features = VariableFeatures(object = Navin_hbca_c52_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c52_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c52_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c52_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c52_panCK, ndims = 60)

Navin_hbca_c52_panCK <- FindNeighbors(Navin_hbca_c52_panCK, dims = 1:40)
Navin_hbca_c52_panCK <- FindClusters(Navin_hbca_c52_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c52_panCK <- RunUMAP(Navin_hbca_c52_panCK, dims = 1:40)
DimPlot(Navin_hbca_c52_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c52_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c52_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c52_panCK.rds")

# Remove Object
rm(Navin_hbca_c52_panCK)
gc()



############################ 
# Navin_hbca_c53_Singlets  #
############################ 
# Load Object
Navin_hbca_c53_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c53_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c53_panCK <- subset(x = Navin_hbca_c53_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c53_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c53_panCK <- FindVariableFeatures(Navin_hbca_c53_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c53_panCK)
Navin_hbca_c53_panCK <- ScaleData(Navin_hbca_c53_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c53_panCK <- RunPCA(Navin_hbca_c53_panCK, features = VariableFeatures(object = Navin_hbca_c53_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c53_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c53_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c53_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c53_panCK, ndims = 60)

Navin_hbca_c53_panCK <- FindNeighbors(Navin_hbca_c53_panCK, dims = 1:40)
Navin_hbca_c53_panCK <- FindClusters(Navin_hbca_c53_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c53_panCK <- RunUMAP(Navin_hbca_c53_panCK, dims = 1:40)
DimPlot(Navin_hbca_c53_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c53_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c53_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c53_panCK.rds")

# Remove Object
rm(Navin_hbca_c53_panCK)
gc()



############################ 
# Navin_hbca_c54_Singlets  #
############################ 
# Load Object
Navin_hbca_c54_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c54_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c54_panCK <- subset(x = Navin_hbca_c54_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c54_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c54_panCK <- FindVariableFeatures(Navin_hbca_c54_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c54_panCK)
Navin_hbca_c54_panCK <- ScaleData(Navin_hbca_c54_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c54_panCK <- RunPCA(Navin_hbca_c54_panCK, features = VariableFeatures(object = Navin_hbca_c54_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c54_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c54_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c54_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c54_panCK, ndims = 60)

Navin_hbca_c54_panCK <- FindNeighbors(Navin_hbca_c54_panCK, dims = 1:40)
Navin_hbca_c54_panCK <- FindClusters(Navin_hbca_c54_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c54_panCK <- RunUMAP(Navin_hbca_c54_panCK, dims = 1:40)
DimPlot(Navin_hbca_c54_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c54_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c54_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c54_panCK.rds")

# Remove Object
rm(Navin_hbca_c54_panCK)
gc()



############################ 
# Navin_hbca_c55_Singlets  #
############################ 
# Load Object
Navin_hbca_c55_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c55_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c55_panCK <- subset(x = Navin_hbca_c55_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c55_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c55_panCK <- FindVariableFeatures(Navin_hbca_c55_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c55_panCK)
Navin_hbca_c55_panCK <- ScaleData(Navin_hbca_c55_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c55_panCK <- RunPCA(Navin_hbca_c55_panCK, features = VariableFeatures(object = Navin_hbca_c55_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c55_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c55_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c55_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c55_panCK, ndims = 60)

Navin_hbca_c55_panCK <- FindNeighbors(Navin_hbca_c55_panCK, dims = 1:40)
Navin_hbca_c55_panCK <- FindClusters(Navin_hbca_c55_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c55_panCK <- RunUMAP(Navin_hbca_c55_panCK, dims = 1:40)
DimPlot(Navin_hbca_c55_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c55_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c55_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c55_panCK.rds")

# Remove Object
rm(Navin_hbca_c55_panCK)
gc()



############################ 
# Navin_hbca_c56_Singlets  #
############################ 
# Load Object
Navin_hbca_c56_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c56_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c56_panCK <- subset(x = Navin_hbca_c56_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c56_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c56_panCK <- FindVariableFeatures(Navin_hbca_c56_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c56_panCK)
Navin_hbca_c56_panCK <- ScaleData(Navin_hbca_c56_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c56_panCK <- RunPCA(Navin_hbca_c56_panCK, features = VariableFeatures(object = Navin_hbca_c56_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c56_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c56_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c56_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c56_panCK, ndims = 60)

Navin_hbca_c56_panCK <- FindNeighbors(Navin_hbca_c56_panCK, dims = 1:40)
Navin_hbca_c56_panCK <- FindClusters(Navin_hbca_c56_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c56_panCK <- RunUMAP(Navin_hbca_c56_panCK, dims = 1:40)
DimPlot(Navin_hbca_c56_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c56_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c56_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c56_panCK.rds")

# Remove Object
rm(Navin_hbca_c56_panCK)
gc()



############################ 
# Navin_hbca_c57_Singlets  #
############################ 
# Load Object
Navin_hbca_c57_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c57_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c57_panCK <- subset(x = Navin_hbca_c57_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c57_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c57_panCK <- FindVariableFeatures(Navin_hbca_c57_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c57_panCK)
Navin_hbca_c57_panCK <- ScaleData(Navin_hbca_c57_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c57_panCK <- RunPCA(Navin_hbca_c57_panCK, features = VariableFeatures(object = Navin_hbca_c57_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c57_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c57_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c57_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c57_panCK, ndims = 60)

Navin_hbca_c57_panCK <- FindNeighbors(Navin_hbca_c57_panCK, dims = 1:40)
Navin_hbca_c57_panCK <- FindClusters(Navin_hbca_c57_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c57_panCK <- RunUMAP(Navin_hbca_c57_panCK, dims = 1:40)
DimPlot(Navin_hbca_c57_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c57_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c57_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c57_panCK.rds")

# Remove Object
rm(Navin_hbca_c57_panCK)
gc()



############################ 
# Navin_hbca_c58_Singlets  #
############################ 
# Load Object
Navin_hbca_c58_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c58_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c58_panCK <- subset(x = Navin_hbca_c58_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c58_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c58_panCK <- FindVariableFeatures(Navin_hbca_c58_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c58_panCK)
Navin_hbca_c58_panCK <- ScaleData(Navin_hbca_c58_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c58_panCK <- RunPCA(Navin_hbca_c58_panCK, features = VariableFeatures(object = Navin_hbca_c58_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c58_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c58_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c58_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c58_panCK, ndims = 60)

Navin_hbca_c58_panCK <- FindNeighbors(Navin_hbca_c58_panCK, dims = 1:40)
Navin_hbca_c58_panCK <- FindClusters(Navin_hbca_c58_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c58_panCK <- RunUMAP(Navin_hbca_c58_panCK, dims = 1:40)
DimPlot(Navin_hbca_c58_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c58_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c58_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c58_panCK.rds")

# Remove Object
rm(Navin_hbca_c58_panCK)
gc()





############################ 
# Navin_hbca_c59_Singlets  #
############################ 
# Load Object
Navin_hbca_c59_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c59_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c59_panCK <- subset(x = Navin_hbca_c59_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c59_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c59_panCK <- FindVariableFeatures(Navin_hbca_c59_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c59_panCK)
Navin_hbca_c59_panCK <- ScaleData(Navin_hbca_c59_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c59_panCK <- RunPCA(Navin_hbca_c59_panCK, features = VariableFeatures(object = Navin_hbca_c59_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c59_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c59_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c59_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c59_panCK, ndims = 60)

Navin_hbca_c59_panCK <- FindNeighbors(Navin_hbca_c59_panCK, dims = 1:40)
Navin_hbca_c59_panCK <- FindClusters(Navin_hbca_c59_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c59_panCK <- RunUMAP(Navin_hbca_c59_panCK, dims = 1:40)
DimPlot(Navin_hbca_c59_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c59_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c59_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c59_panCK.rds")

# Remove Object
rm(Navin_hbca_c59_panCK)
gc()



############################ 
# Navin_hbca_c60_Singlets  #
############################ 
# Load Object
Navin_hbca_c60_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c60_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c60_panCK <- subset(x = Navin_hbca_c60_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c60_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c60_panCK <- FindVariableFeatures(Navin_hbca_c60_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c60_panCK)
Navin_hbca_c60_panCK <- ScaleData(Navin_hbca_c60_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c60_panCK <- RunPCA(Navin_hbca_c60_panCK, features = VariableFeatures(object = Navin_hbca_c60_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c60_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c60_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c60_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c60_panCK, ndims = 60)

Navin_hbca_c60_panCK <- FindNeighbors(Navin_hbca_c60_panCK, dims = 1:40)
Navin_hbca_c60_panCK <- FindClusters(Navin_hbca_c60_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c60_panCK <- RunUMAP(Navin_hbca_c60_panCK, dims = 1:40)
DimPlot(Navin_hbca_c60_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c60_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c60_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c60_panCK.rds")

# Remove Object
rm(Navin_hbca_c60_panCK)
gc()



############################ 
# Navin_hbca_c61_Singlets  #
############################ 
# Load Object
Navin_hbca_c61_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c61_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c61_panCK <- subset(x = Navin_hbca_c61_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c61_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c61_panCK <- FindVariableFeatures(Navin_hbca_c61_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c61_panCK)
Navin_hbca_c61_panCK <- ScaleData(Navin_hbca_c61_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c61_panCK <- RunPCA(Navin_hbca_c61_panCK, features = VariableFeatures(object = Navin_hbca_c61_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c61_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c61_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c61_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c61_panCK, ndims = 60)

Navin_hbca_c61_panCK <- FindNeighbors(Navin_hbca_c61_panCK, dims = 1:40)
Navin_hbca_c61_panCK <- FindClusters(Navin_hbca_c61_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c61_panCK <- RunUMAP(Navin_hbca_c61_panCK, dims = 1:40)
DimPlot(Navin_hbca_c61_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c61_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c61_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c61_panCK.rds")

# Remove Object
rm(Navin_hbca_c61_panCK)
gc()



############################ 
# Navin_hbca_c62_Singlets  #
############################ 
# Load Object
Navin_hbca_c62_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c62_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c62_panCK <- subset(x = Navin_hbca_c62_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c62_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c62_panCK <- FindVariableFeatures(Navin_hbca_c62_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c62_panCK)
Navin_hbca_c62_panCK <- ScaleData(Navin_hbca_c62_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c62_panCK <- RunPCA(Navin_hbca_c62_panCK, features = VariableFeatures(object = Navin_hbca_c62_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c62_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c62_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c62_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c62_panCK, ndims = 60)

Navin_hbca_c62_panCK <- FindNeighbors(Navin_hbca_c62_panCK, dims = 1:40)
Navin_hbca_c62_panCK <- FindClusters(Navin_hbca_c62_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c62_panCK <- RunUMAP(Navin_hbca_c62_panCK, dims = 1:40)
DimPlot(Navin_hbca_c62_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c62_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c62_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c62_panCK.rds")

# Remove Object
rm(Navin_hbca_c62_panCK)
gc()



############################ 
# Navin_hbca_c63_Singlets  #
############################ 
# Load Object
Navin_hbca_c63_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c63_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c63_panCK <- subset(x = Navin_hbca_c63_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c63_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c63_panCK <- FindVariableFeatures(Navin_hbca_c63_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c63_panCK)
Navin_hbca_c63_panCK <- ScaleData(Navin_hbca_c63_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c63_panCK <- RunPCA(Navin_hbca_c63_panCK, features = VariableFeatures(object = Navin_hbca_c63_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c63_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c63_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c63_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c63_panCK, ndims = 60)

Navin_hbca_c63_panCK <- FindNeighbors(Navin_hbca_c63_panCK, dims = 1:40)
Navin_hbca_c63_panCK <- FindClusters(Navin_hbca_c63_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c63_panCK <- RunUMAP(Navin_hbca_c63_panCK, dims = 1:40)
DimPlot(Navin_hbca_c63_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c63_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c63_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c63_panCK.rds")

# Remove Object
rm(Navin_hbca_c63_panCK)
gc()



############################ 
# Navin_hbca_c64_Singlets  #
############################ 
# Load Object
Navin_hbca_c64_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c64_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c64_panCK <- subset(x = Navin_hbca_c64_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c64_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c64_panCK <- FindVariableFeatures(Navin_hbca_c64_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c64_panCK)
Navin_hbca_c64_panCK <- ScaleData(Navin_hbca_c64_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c64_panCK <- RunPCA(Navin_hbca_c64_panCK, features = VariableFeatures(object = Navin_hbca_c64_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c64_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c64_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c64_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c64_panCK, ndims = 60)

Navin_hbca_c64_panCK <- FindNeighbors(Navin_hbca_c64_panCK, dims = 1:40)
Navin_hbca_c64_panCK <- FindClusters(Navin_hbca_c64_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c64_panCK <- RunUMAP(Navin_hbca_c64_panCK, dims = 1:40)
DimPlot(Navin_hbca_c64_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c64_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c64_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c64_panCK.rds")

# Remove Object
rm(Navin_hbca_c64_panCK)
gc()



############################ 
# Navin_hbca_c65_Singlets  #
############################ 
# Load Object
Navin_hbca_c65_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c65_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c65_panCK <- subset(x = Navin_hbca_c65_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c65_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c65_panCK <- FindVariableFeatures(Navin_hbca_c65_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c65_panCK)
Navin_hbca_c65_panCK <- ScaleData(Navin_hbca_c65_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c65_panCK <- RunPCA(Navin_hbca_c65_panCK, features = VariableFeatures(object = Navin_hbca_c65_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c65_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c65_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c65_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c65_panCK, ndims = 60)

Navin_hbca_c65_panCK <- FindNeighbors(Navin_hbca_c65_panCK, dims = 1:40)
Navin_hbca_c65_panCK <- FindClusters(Navin_hbca_c65_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c65_panCK <- RunUMAP(Navin_hbca_c65_panCK, dims = 1:40)
DimPlot(Navin_hbca_c65_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c65_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c65_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c65_panCK.rds")

# Remove Object
rm(Navin_hbca_c65_panCK)
gc()



############################ 
# Navin_hbca_c66_Singlets  #
############################ 
# Load Object
Navin_hbca_c66_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c66_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c66_panCK <- subset(x = Navin_hbca_c66_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c66_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c66_panCK <- FindVariableFeatures(Navin_hbca_c66_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c66_panCK)
Navin_hbca_c66_panCK <- ScaleData(Navin_hbca_c66_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c66_panCK <- RunPCA(Navin_hbca_c66_panCK, features = VariableFeatures(object = Navin_hbca_c66_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c66_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c66_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c66_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c66_panCK, ndims = 60)

Navin_hbca_c66_panCK <- FindNeighbors(Navin_hbca_c66_panCK, dims = 1:40)
Navin_hbca_c66_panCK <- FindClusters(Navin_hbca_c66_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c66_panCK <- RunUMAP(Navin_hbca_c66_panCK, dims = 1:40)
DimPlot(Navin_hbca_c66_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c66_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c66_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c66_panCK.rds")

# Remove Object
rm(Navin_hbca_c66_panCK)
gc()



############################ 
# Navin_hbca_c67_Singlets  #
############################ 
# Load Object
Navin_hbca_c67_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c67_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c67_panCK <- subset(x = Navin_hbca_c67_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c67_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c67_panCK <- FindVariableFeatures(Navin_hbca_c67_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c67_panCK)
Navin_hbca_c67_panCK <- ScaleData(Navin_hbca_c67_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c67_panCK <- RunPCA(Navin_hbca_c67_panCK, features = VariableFeatures(object = Navin_hbca_c67_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c67_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c67_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c67_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c67_panCK, ndims = 60)

Navin_hbca_c67_panCK <- FindNeighbors(Navin_hbca_c67_panCK, dims = 1:40)
Navin_hbca_c67_panCK <- FindClusters(Navin_hbca_c67_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c67_panCK <- RunUMAP(Navin_hbca_c67_panCK, dims = 1:40)
DimPlot(Navin_hbca_c67_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c67_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c67_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c67_panCK.rds")

# Remove Object
rm(Navin_hbca_c67_panCK)
gc()



############################ 
# Navin_hbca_c68_Singlets  #
############################ 
# Load Object
Navin_hbca_c68_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c68_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c68_panCK <- subset(x = Navin_hbca_c68_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c68_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c68_panCK <- FindVariableFeatures(Navin_hbca_c68_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c68_panCK)
Navin_hbca_c68_panCK <- ScaleData(Navin_hbca_c68_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c68_panCK <- RunPCA(Navin_hbca_c68_panCK, features = VariableFeatures(object = Navin_hbca_c68_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c68_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c68_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c68_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c68_panCK, ndims = 60)

Navin_hbca_c68_panCK <- FindNeighbors(Navin_hbca_c68_panCK, dims = 1:40)
Navin_hbca_c68_panCK <- FindClusters(Navin_hbca_c68_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c68_panCK <- RunUMAP(Navin_hbca_c68_panCK, dims = 1:40)
DimPlot(Navin_hbca_c68_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c68_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c68_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c68_panCK.rds")

# Remove Object
rm(Navin_hbca_c68_panCK)
gc()





############################ 
# Navin_hbca_c69_Singlets  #
############################ 
# Load Object
Navin_hbca_c69_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c69_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c69_panCK <- subset(x = Navin_hbca_c69_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c69_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c69_panCK <- FindVariableFeatures(Navin_hbca_c69_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c69_panCK)
Navin_hbca_c69_panCK <- ScaleData(Navin_hbca_c69_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c69_panCK <- RunPCA(Navin_hbca_c69_panCK, features = VariableFeatures(object = Navin_hbca_c69_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c69_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c69_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c69_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c69_panCK, ndims = 60)

Navin_hbca_c69_panCK <- FindNeighbors(Navin_hbca_c69_panCK, dims = 1:40)
Navin_hbca_c69_panCK <- FindClusters(Navin_hbca_c69_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c69_panCK <- RunUMAP(Navin_hbca_c69_panCK, dims = 1:40)
DimPlot(Navin_hbca_c69_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c69_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c69_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c69_panCK.rds")

# Remove Object
rm(Navin_hbca_c69_panCK)
gc()



############################ 
# Navin_hbca_c70_Singlets  #
############################ 
# Load Object
Navin_hbca_c70_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c70_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c70_panCK <- subset(x = Navin_hbca_c70_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c70_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c70_panCK <- FindVariableFeatures(Navin_hbca_c70_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c70_panCK)
Navin_hbca_c70_panCK <- ScaleData(Navin_hbca_c70_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c70_panCK <- RunPCA(Navin_hbca_c70_panCK, features = VariableFeatures(object = Navin_hbca_c70_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c70_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c70_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c70_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c70_panCK, ndims = 60)

Navin_hbca_c70_panCK <- FindNeighbors(Navin_hbca_c70_panCK, dims = 1:40)
Navin_hbca_c70_panCK <- FindClusters(Navin_hbca_c70_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c70_panCK <- RunUMAP(Navin_hbca_c70_panCK, dims = 1:40)
DimPlot(Navin_hbca_c70_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c70_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c70_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c70_panCK.rds")

# Remove Object
rm(Navin_hbca_c70_panCK)
gc()



############################ 
# Navin_hbca_c71_Singlets  #
############################ 
# Load Object
Navin_hbca_c71_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c71_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c71_panCK <- subset(x = Navin_hbca_c71_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c71_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c71_panCK <- FindVariableFeatures(Navin_hbca_c71_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c71_panCK)
Navin_hbca_c71_panCK <- ScaleData(Navin_hbca_c71_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c71_panCK <- RunPCA(Navin_hbca_c71_panCK, features = VariableFeatures(object = Navin_hbca_c71_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c71_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c71_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c71_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c71_panCK, ndims = 60)

Navin_hbca_c71_panCK <- FindNeighbors(Navin_hbca_c71_panCK, dims = 1:40)
Navin_hbca_c71_panCK <- FindClusters(Navin_hbca_c71_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c71_panCK <- RunUMAP(Navin_hbca_c71_panCK, dims = 1:40)
DimPlot(Navin_hbca_c71_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c71_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c71_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c71_panCK.rds")

# Remove Object
rm(Navin_hbca_c71_panCK)
gc()



############################ 
# Navin_hbca_c72_Singlets  #
############################ 
# Load Object
Navin_hbca_c72_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c72_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c72_panCK <- subset(x = Navin_hbca_c72_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c72_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c72_panCK <- FindVariableFeatures(Navin_hbca_c72_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c72_panCK)
Navin_hbca_c72_panCK <- ScaleData(Navin_hbca_c72_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c72_panCK <- RunPCA(Navin_hbca_c72_panCK, features = VariableFeatures(object = Navin_hbca_c72_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c72_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c72_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c72_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c72_panCK, ndims = 60)

Navin_hbca_c72_panCK <- FindNeighbors(Navin_hbca_c72_panCK, dims = 1:40)
Navin_hbca_c72_panCK <- FindClusters(Navin_hbca_c72_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c72_panCK <- RunUMAP(Navin_hbca_c72_panCK, dims = 1:40)
DimPlot(Navin_hbca_c72_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c72_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c72_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c72_panCK.rds")

# Remove Object
rm(Navin_hbca_c72_panCK)
gc()



############################ 
# Navin_hbca_c73_Singlets  #
############################ 
# Load Object
Navin_hbca_c73_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c73_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c73_panCK <- subset(x = Navin_hbca_c73_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c73_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c73_panCK <- FindVariableFeatures(Navin_hbca_c73_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c73_panCK)
Navin_hbca_c73_panCK <- ScaleData(Navin_hbca_c73_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c73_panCK <- RunPCA(Navin_hbca_c73_panCK, features = VariableFeatures(object = Navin_hbca_c73_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c73_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c73_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c73_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c73_panCK, ndims = 60)

Navin_hbca_c73_panCK <- FindNeighbors(Navin_hbca_c73_panCK, dims = 1:40)
Navin_hbca_c73_panCK <- FindClusters(Navin_hbca_c73_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c73_panCK <- RunUMAP(Navin_hbca_c73_panCK, dims = 1:40)
DimPlot(Navin_hbca_c73_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c73_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c73_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c73_panCK.rds")

# Remove Object
rm(Navin_hbca_c73_panCK)
gc()



############################ 
# Navin_hbca_c74_Singlets  #
############################ 
# Load Object
Navin_hbca_c74_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c74_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c74_panCK <- subset(x = Navin_hbca_c74_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c74_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c74_panCK <- FindVariableFeatures(Navin_hbca_c74_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c74_panCK)
Navin_hbca_c74_panCK <- ScaleData(Navin_hbca_c74_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c74_panCK <- RunPCA(Navin_hbca_c74_panCK, features = VariableFeatures(object = Navin_hbca_c74_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c74_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c74_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c74_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c74_panCK, ndims = 60)

Navin_hbca_c74_panCK <- FindNeighbors(Navin_hbca_c74_panCK, dims = 1:40)
Navin_hbca_c74_panCK <- FindClusters(Navin_hbca_c74_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c74_panCK <- RunUMAP(Navin_hbca_c74_panCK, dims = 1:40)
DimPlot(Navin_hbca_c74_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c74_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c74_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c74_panCK.rds")

# Remove Object
rm(Navin_hbca_c74_panCK)
gc()



############################ 
# Navin_hbca_c75_Singlets  #
############################ 
# Load Object
Navin_hbca_c75_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c75_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c75_panCK <- subset(x = Navin_hbca_c75_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c75_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c75_panCK <- FindVariableFeatures(Navin_hbca_c75_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c75_panCK)
Navin_hbca_c75_panCK <- ScaleData(Navin_hbca_c75_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c75_panCK <- RunPCA(Navin_hbca_c75_panCK, features = VariableFeatures(object = Navin_hbca_c75_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c75_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c75_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c75_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c75_panCK, ndims = 60)

Navin_hbca_c75_panCK <- FindNeighbors(Navin_hbca_c75_panCK, dims = 1:40)
Navin_hbca_c75_panCK <- FindClusters(Navin_hbca_c75_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c75_panCK <- RunUMAP(Navin_hbca_c75_panCK, dims = 1:40)
DimPlot(Navin_hbca_c75_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c75_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c75_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c75_panCK.rds")

# Remove Object
rm(Navin_hbca_c75_panCK)
gc()



############################ 
# Navin_hbca_c76_Singlets  #
############################ 
# Load Object
Navin_hbca_c76_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c76_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c76_panCK <- subset(x = Navin_hbca_c76_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c76_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c76_panCK <- FindVariableFeatures(Navin_hbca_c76_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c76_panCK)
Navin_hbca_c76_panCK <- ScaleData(Navin_hbca_c76_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c76_panCK <- RunPCA(Navin_hbca_c76_panCK, features = VariableFeatures(object = Navin_hbca_c76_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c76_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c76_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c76_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c76_panCK, ndims = 60)

Navin_hbca_c76_panCK <- FindNeighbors(Navin_hbca_c76_panCK, dims = 1:40)
Navin_hbca_c76_panCK <- FindClusters(Navin_hbca_c76_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c76_panCK <- RunUMAP(Navin_hbca_c76_panCK, dims = 1:40)
DimPlot(Navin_hbca_c76_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c76_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c76_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c76_panCK.rds")

# Remove Object
rm(Navin_hbca_c76_panCK)
gc()



############################ 
# Navin_hbca_c77_Singlets  #
############################ 
# Load Object
Navin_hbca_c77_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c77_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c77_panCK <- subset(x = Navin_hbca_c77_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c77_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c77_panCK <- FindVariableFeatures(Navin_hbca_c77_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c77_panCK)
Navin_hbca_c77_panCK <- ScaleData(Navin_hbca_c77_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c77_panCK <- RunPCA(Navin_hbca_c77_panCK, features = VariableFeatures(object = Navin_hbca_c77_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c77_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c77_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c77_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c77_panCK, ndims = 60)

Navin_hbca_c77_panCK <- FindNeighbors(Navin_hbca_c77_panCK, dims = 1:40)
Navin_hbca_c77_panCK <- FindClusters(Navin_hbca_c77_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c77_panCK <- RunUMAP(Navin_hbca_c77_panCK, dims = 1:40)
DimPlot(Navin_hbca_c77_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c77_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c77_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c77_panCK.rds")

# Remove Object
rm(Navin_hbca_c77_panCK)
gc()



############################ 
# Navin_hbca_c78_Singlets  #
############################ 
# Load Object
Navin_hbca_c78_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c78_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c78_panCK <- subset(x = Navin_hbca_c78_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c78_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c78_panCK <- FindVariableFeatures(Navin_hbca_c78_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c78_panCK)
Navin_hbca_c78_panCK <- ScaleData(Navin_hbca_c78_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c78_panCK <- RunPCA(Navin_hbca_c78_panCK, features = VariableFeatures(object = Navin_hbca_c78_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c78_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c78_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c78_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c78_panCK, ndims = 60)

Navin_hbca_c78_panCK <- FindNeighbors(Navin_hbca_c78_panCK, dims = 1:40)
Navin_hbca_c78_panCK <- FindClusters(Navin_hbca_c78_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c78_panCK <- RunUMAP(Navin_hbca_c78_panCK, dims = 1:40)
DimPlot(Navin_hbca_c78_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c78_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c78_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c78_panCK.rds")

# Remove Object
rm(Navin_hbca_c78_panCK)
gc()





############################ 
# Navin_hbca_c79_Singlets  #
############################ 
# Load Object
Navin_hbca_c79_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c79_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c79_panCK <- subset(x = Navin_hbca_c79_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c79_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c79_panCK <- FindVariableFeatures(Navin_hbca_c79_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c79_panCK)
Navin_hbca_c79_panCK <- ScaleData(Navin_hbca_c79_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c79_panCK <- RunPCA(Navin_hbca_c79_panCK, features = VariableFeatures(object = Navin_hbca_c79_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c79_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c79_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c79_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c79_panCK, ndims = 60)

Navin_hbca_c79_panCK <- FindNeighbors(Navin_hbca_c79_panCK, dims = 1:40)
Navin_hbca_c79_panCK <- FindClusters(Navin_hbca_c79_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c79_panCK <- RunUMAP(Navin_hbca_c79_panCK, dims = 1:40)
DimPlot(Navin_hbca_c79_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c79_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c79_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c79_panCK.rds")

# Remove Object
rm(Navin_hbca_c79_panCK)
gc()



############################ 
# Navin_hbca_c80_Singlets  #
############################ 
# Load Object
Navin_hbca_c80_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c80_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c80_panCK <- subset(x = Navin_hbca_c80_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c80_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c80_panCK <- FindVariableFeatures(Navin_hbca_c80_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c80_panCK)
Navin_hbca_c80_panCK <- ScaleData(Navin_hbca_c80_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c80_panCK <- RunPCA(Navin_hbca_c80_panCK, features = VariableFeatures(object = Navin_hbca_c80_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c80_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c80_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c80_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c80_panCK, ndims = 60)

Navin_hbca_c80_panCK <- FindNeighbors(Navin_hbca_c80_panCK, dims = 1:40)
Navin_hbca_c80_panCK <- FindClusters(Navin_hbca_c80_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c80_panCK <- RunUMAP(Navin_hbca_c80_panCK, dims = 1:40)
DimPlot(Navin_hbca_c80_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c80_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c80_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c80_panCK.rds")

# Remove Object
rm(Navin_hbca_c80_panCK)
gc()



############################ 
# Navin_hbca_c81_Singlets  #
############################ 
# Load Object
Navin_hbca_c81_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c81_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c81_panCK <- subset(x = Navin_hbca_c81_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c81_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c81_panCK <- FindVariableFeatures(Navin_hbca_c81_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c81_panCK)
Navin_hbca_c81_panCK <- ScaleData(Navin_hbca_c81_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c81_panCK <- RunPCA(Navin_hbca_c81_panCK, features = VariableFeatures(object = Navin_hbca_c81_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c81_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c81_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c81_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c81_panCK, ndims = 60)

Navin_hbca_c81_panCK <- FindNeighbors(Navin_hbca_c81_panCK, dims = 1:40)
Navin_hbca_c81_panCK <- FindClusters(Navin_hbca_c81_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c81_panCK <- RunUMAP(Navin_hbca_c81_panCK, dims = 1:40)
DimPlot(Navin_hbca_c81_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c81_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c81_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c81_panCK.rds")

# Remove Object
rm(Navin_hbca_c81_panCK)
gc()



############################ 
# Navin_hbca_c82_Singlets  #
############################ 
# Load Object
Navin_hbca_c82_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c82_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c82_panCK <- subset(x = Navin_hbca_c82_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c82_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c82_panCK <- FindVariableFeatures(Navin_hbca_c82_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c82_panCK)
Navin_hbca_c82_panCK <- ScaleData(Navin_hbca_c82_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c82_panCK <- RunPCA(Navin_hbca_c82_panCK, features = VariableFeatures(object = Navin_hbca_c82_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c82_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c82_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c82_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c82_panCK, ndims = 60)

Navin_hbca_c82_panCK <- FindNeighbors(Navin_hbca_c82_panCK, dims = 1:40)
Navin_hbca_c82_panCK <- FindClusters(Navin_hbca_c82_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c82_panCK <- RunUMAP(Navin_hbca_c82_panCK, dims = 1:40)
DimPlot(Navin_hbca_c82_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c82_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c82_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c82_panCK.rds")

# Remove Object
rm(Navin_hbca_c82_panCK)
gc()



############################ 
# Navin_hbca_c83_Singlets  #
############################ 
# Load Object
Navin_hbca_c83_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c83_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c83_panCK <- subset(x = Navin_hbca_c83_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c83_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c83_panCK <- FindVariableFeatures(Navin_hbca_c83_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c83_panCK)
Navin_hbca_c83_panCK <- ScaleData(Navin_hbca_c83_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c83_panCK <- RunPCA(Navin_hbca_c83_panCK, features = VariableFeatures(object = Navin_hbca_c83_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c83_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c83_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c83_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c83_panCK, ndims = 60)

Navin_hbca_c83_panCK <- FindNeighbors(Navin_hbca_c83_panCK, dims = 1:40)
Navin_hbca_c83_panCK <- FindClusters(Navin_hbca_c83_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c83_panCK <- RunUMAP(Navin_hbca_c83_panCK, dims = 1:40)
DimPlot(Navin_hbca_c83_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c83_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c83_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c83_panCK.rds")

# Remove Object
rm(Navin_hbca_c83_panCK)
gc()



############################ 
# Navin_hbca_c84_Singlets  #
############################ 
# Load Object
Navin_hbca_c84_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c84_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c84_panCK <- subset(x = Navin_hbca_c84_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c84_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c84_panCK <- FindVariableFeatures(Navin_hbca_c84_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c84_panCK)
Navin_hbca_c84_panCK <- ScaleData(Navin_hbca_c84_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c84_panCK <- RunPCA(Navin_hbca_c84_panCK, features = VariableFeatures(object = Navin_hbca_c84_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c84_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c84_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c84_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c84_panCK, ndims = 60)

Navin_hbca_c84_panCK <- FindNeighbors(Navin_hbca_c84_panCK, dims = 1:40)
Navin_hbca_c84_panCK <- FindClusters(Navin_hbca_c84_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c84_panCK <- RunUMAP(Navin_hbca_c84_panCK, dims = 1:40)
DimPlot(Navin_hbca_c84_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c84_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c84_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c84_panCK.rds")

# Remove Object
rm(Navin_hbca_c84_panCK)
gc()



############################ 
# Navin_hbca_c85_Singlets  #
############################ 
# Load Object
Navin_hbca_c85_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c85_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c85_panCK <- subset(x = Navin_hbca_c85_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c85_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c85_panCK <- FindVariableFeatures(Navin_hbca_c85_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c85_panCK)
Navin_hbca_c85_panCK <- ScaleData(Navin_hbca_c85_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c85_panCK <- RunPCA(Navin_hbca_c85_panCK, features = VariableFeatures(object = Navin_hbca_c85_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c85_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c85_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c85_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c85_panCK, ndims = 60)

Navin_hbca_c85_panCK <- FindNeighbors(Navin_hbca_c85_panCK, dims = 1:40)
Navin_hbca_c85_panCK <- FindClusters(Navin_hbca_c85_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c85_panCK <- RunUMAP(Navin_hbca_c85_panCK, dims = 1:40)
DimPlot(Navin_hbca_c85_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c85_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c85_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c85_panCK.rds")

# Remove Object
rm(Navin_hbca_c85_panCK)
gc()



############################ 
# Navin_hbca_c86_Singlets  #
############################ 
# Load Object
Navin_hbca_c86_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c86_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c86_panCK <- subset(x = Navin_hbca_c86_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c86_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c86_panCK <- FindVariableFeatures(Navin_hbca_c86_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c86_panCK)
Navin_hbca_c86_panCK <- ScaleData(Navin_hbca_c86_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c86_panCK <- RunPCA(Navin_hbca_c86_panCK, features = VariableFeatures(object = Navin_hbca_c86_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c86_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c86_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c86_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c86_panCK, ndims = 60)

Navin_hbca_c86_panCK <- FindNeighbors(Navin_hbca_c86_panCK, dims = 1:40)
Navin_hbca_c86_panCK <- FindClusters(Navin_hbca_c86_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c86_panCK <- RunUMAP(Navin_hbca_c86_panCK, dims = 1:40)
DimPlot(Navin_hbca_c86_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c86_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c86_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c86_panCK.rds")

# Remove Object
rm(Navin_hbca_c86_panCK)
gc()



############################ 
# Navin_hbca_c87_Singlets  #
############################ 
# Load Object
Navin_hbca_c87_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c87_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c87_panCK <- subset(x = Navin_hbca_c87_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c87_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c87_panCK <- FindVariableFeatures(Navin_hbca_c87_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c87_panCK)
Navin_hbca_c87_panCK <- ScaleData(Navin_hbca_c87_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c87_panCK <- RunPCA(Navin_hbca_c87_panCK, features = VariableFeatures(object = Navin_hbca_c87_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c87_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c87_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c87_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c87_panCK, ndims = 60)

Navin_hbca_c87_panCK <- FindNeighbors(Navin_hbca_c87_panCK, dims = 1:40)
Navin_hbca_c87_panCK <- FindClusters(Navin_hbca_c87_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c87_panCK <- RunUMAP(Navin_hbca_c87_panCK, dims = 1:40)
DimPlot(Navin_hbca_c87_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c87_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c87_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c87_panCK.rds")

# Remove Object
rm(Navin_hbca_c87_panCK)
gc()



############################ 
# Navin_hbca_c88_Singlets  #
############################ 
# Load Object
Navin_hbca_c88_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c88_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c88_panCK <- subset(x = Navin_hbca_c88_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c88_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c88_panCK <- FindVariableFeatures(Navin_hbca_c88_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c88_panCK)
Navin_hbca_c88_panCK <- ScaleData(Navin_hbca_c88_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c88_panCK <- RunPCA(Navin_hbca_c88_panCK, features = VariableFeatures(object = Navin_hbca_c88_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c88_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c88_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c88_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c88_panCK, ndims = 60)

Navin_hbca_c88_panCK <- FindNeighbors(Navin_hbca_c88_panCK, dims = 1:40)
Navin_hbca_c88_panCK <- FindClusters(Navin_hbca_c88_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c88_panCK <- RunUMAP(Navin_hbca_c88_panCK, dims = 1:40)
DimPlot(Navin_hbca_c88_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c88_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c88_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c88_panCK.rds")

# Remove Object
rm(Navin_hbca_c88_panCK)
gc()





############################ 
# Navin_hbca_c89_Singlets  #
############################ 
# Load Object
Navin_hbca_c89_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c89_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c89_panCK <- subset(x = Navin_hbca_c89_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c89_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c89_panCK <- FindVariableFeatures(Navin_hbca_c89_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c89_panCK)
Navin_hbca_c89_panCK <- ScaleData(Navin_hbca_c89_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c89_panCK <- RunPCA(Navin_hbca_c89_panCK, features = VariableFeatures(object = Navin_hbca_c89_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c89_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c89_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c89_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c89_panCK, ndims = 60)

Navin_hbca_c89_panCK <- FindNeighbors(Navin_hbca_c89_panCK, dims = 1:40)
Navin_hbca_c89_panCK <- FindClusters(Navin_hbca_c89_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c89_panCK <- RunUMAP(Navin_hbca_c89_panCK, dims = 1:40)
DimPlot(Navin_hbca_c89_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c89_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c89_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c89_panCK.rds")

# Remove Object
rm(Navin_hbca_c89_panCK)
gc()



############################ 
# Navin_hbca_c90_Singlets  #
############################ 
# Load Object
Navin_hbca_c90_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c90_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c90_panCK <- subset(x = Navin_hbca_c90_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c90_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c90_panCK <- FindVariableFeatures(Navin_hbca_c90_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c90_panCK)
Navin_hbca_c90_panCK <- ScaleData(Navin_hbca_c90_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c90_panCK <- RunPCA(Navin_hbca_c90_panCK, features = VariableFeatures(object = Navin_hbca_c90_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c90_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c90_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c90_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c90_panCK, ndims = 60)

Navin_hbca_c90_panCK <- FindNeighbors(Navin_hbca_c90_panCK, dims = 1:40)
Navin_hbca_c90_panCK <- FindClusters(Navin_hbca_c90_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c90_panCK <- RunUMAP(Navin_hbca_c90_panCK, dims = 1:40)
DimPlot(Navin_hbca_c90_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c90_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c90_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c90_panCK.rds")

# Remove Object
rm(Navin_hbca_c90_panCK)
gc()



############################ 
# Navin_hbca_c91_Singlets  #
############################ 
# Load Object
Navin_hbca_c91_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c91_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c91_panCK <- subset(x = Navin_hbca_c91_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c91_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c91_panCK <- FindVariableFeatures(Navin_hbca_c91_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c91_panCK)
Navin_hbca_c91_panCK <- ScaleData(Navin_hbca_c91_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c91_panCK <- RunPCA(Navin_hbca_c91_panCK, features = VariableFeatures(object = Navin_hbca_c91_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c91_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c91_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c91_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c91_panCK, ndims = 60)

Navin_hbca_c91_panCK <- FindNeighbors(Navin_hbca_c91_panCK, dims = 1:40)
Navin_hbca_c91_panCK <- FindClusters(Navin_hbca_c91_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c91_panCK <- RunUMAP(Navin_hbca_c91_panCK, dims = 1:40)
DimPlot(Navin_hbca_c91_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c91_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c91_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c91_panCK.rds")

# Remove Object
rm(Navin_hbca_c91_panCK)
gc()



############################ 
# Navin_hbca_c92_Singlets  #
############################ 
# Load Object
Navin_hbca_c92_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c92_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c92_panCK <- subset(x = Navin_hbca_c92_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c92_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c92_panCK <- FindVariableFeatures(Navin_hbca_c92_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c92_panCK)
Navin_hbca_c92_panCK <- ScaleData(Navin_hbca_c92_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c92_panCK <- RunPCA(Navin_hbca_c92_panCK, features = VariableFeatures(object = Navin_hbca_c92_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c92_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c92_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c92_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c92_panCK, ndims = 60)

Navin_hbca_c92_panCK <- FindNeighbors(Navin_hbca_c92_panCK, dims = 1:40)
Navin_hbca_c92_panCK <- FindClusters(Navin_hbca_c92_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c92_panCK <- RunUMAP(Navin_hbca_c92_panCK, dims = 1:40)
DimPlot(Navin_hbca_c92_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c92_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c92_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c92_panCK.rds")

# Remove Object
rm(Navin_hbca_c92_panCK)
gc()



############################ 
# Navin_hbca_c93_Singlets  #
############################ 
# Load Object
Navin_hbca_c93_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c93_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c93_panCK <- subset(x = Navin_hbca_c93_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c93_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c93_panCK <- FindVariableFeatures(Navin_hbca_c93_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c93_panCK)
Navin_hbca_c93_panCK <- ScaleData(Navin_hbca_c93_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c93_panCK <- RunPCA(Navin_hbca_c93_panCK, features = VariableFeatures(object = Navin_hbca_c93_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c93_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c93_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c93_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c93_panCK, ndims = 60)

Navin_hbca_c93_panCK <- FindNeighbors(Navin_hbca_c93_panCK, dims = 1:40)
Navin_hbca_c93_panCK <- FindClusters(Navin_hbca_c93_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c93_panCK <- RunUMAP(Navin_hbca_c93_panCK, dims = 1:40)
DimPlot(Navin_hbca_c93_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c93_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c93_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c93_panCK.rds")

# Remove Object
rm(Navin_hbca_c93_panCK)
gc()



############################ 
# Navin_hbca_c94_Singlets  #
############################ 
# Load Object
Navin_hbca_c94_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c94_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c94_panCK <- subset(x = Navin_hbca_c94_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c94_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c94_panCK <- FindVariableFeatures(Navin_hbca_c94_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c94_panCK)
Navin_hbca_c94_panCK <- ScaleData(Navin_hbca_c94_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c94_panCK <- RunPCA(Navin_hbca_c94_panCK, features = VariableFeatures(object = Navin_hbca_c94_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c94_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c94_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c94_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c94_panCK, ndims = 60)

Navin_hbca_c94_panCK <- FindNeighbors(Navin_hbca_c94_panCK, dims = 1:40)
Navin_hbca_c94_panCK <- FindClusters(Navin_hbca_c94_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c94_panCK <- RunUMAP(Navin_hbca_c94_panCK, dims = 1:40)
DimPlot(Navin_hbca_c94_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c94_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c94_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c94_panCK.rds")

# Remove Object
rm(Navin_hbca_c94_panCK)
gc()



############################ 
# Navin_hbca_c95_Singlets  #
############################ 
# Load Object
Navin_hbca_c95_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c95_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c95_panCK <- subset(x = Navin_hbca_c95_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c95_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c95_panCK <- FindVariableFeatures(Navin_hbca_c95_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c95_panCK)
Navin_hbca_c95_panCK <- ScaleData(Navin_hbca_c95_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c95_panCK <- RunPCA(Navin_hbca_c95_panCK, features = VariableFeatures(object = Navin_hbca_c95_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c95_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c95_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c95_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c95_panCK, ndims = 60)

Navin_hbca_c95_panCK <- FindNeighbors(Navin_hbca_c95_panCK, dims = 1:40)
Navin_hbca_c95_panCK <- FindClusters(Navin_hbca_c95_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c95_panCK <- RunUMAP(Navin_hbca_c95_panCK, dims = 1:40)
DimPlot(Navin_hbca_c95_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c95_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c95_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c95_panCK.rds")

# Remove Object
rm(Navin_hbca_c95_panCK)
gc()



############################ 
# Navin_hbca_c96_Singlets  #
############################ 
# Load Object
Navin_hbca_c96_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c96_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c96_panCK <- subset(x = Navin_hbca_c96_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c96_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c96_panCK <- FindVariableFeatures(Navin_hbca_c96_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c96_panCK)
Navin_hbca_c96_panCK <- ScaleData(Navin_hbca_c96_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c96_panCK <- RunPCA(Navin_hbca_c96_panCK, features = VariableFeatures(object = Navin_hbca_c96_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c96_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c96_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c96_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c96_panCK, ndims = 60)

Navin_hbca_c96_panCK <- FindNeighbors(Navin_hbca_c96_panCK, dims = 1:40)
Navin_hbca_c96_panCK <- FindClusters(Navin_hbca_c96_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c96_panCK <- RunUMAP(Navin_hbca_c96_panCK, dims = 1:40)
DimPlot(Navin_hbca_c96_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c96_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c96_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c96_panCK.rds")

# Remove Object
rm(Navin_hbca_c96_panCK)
gc()



############################ 
# Navin_hbca_c97_Singlets  #
############################ 
# Load Object
Navin_hbca_c97_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c97_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c97_panCK <- subset(x = Navin_hbca_c97_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c97_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c97_panCK <- FindVariableFeatures(Navin_hbca_c97_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c97_panCK)
Navin_hbca_c97_panCK <- ScaleData(Navin_hbca_c97_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c97_panCK <- RunPCA(Navin_hbca_c97_panCK, features = VariableFeatures(object = Navin_hbca_c97_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c97_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c97_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c97_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c97_panCK, ndims = 60)

Navin_hbca_c97_panCK <- FindNeighbors(Navin_hbca_c97_panCK, dims = 1:40)
Navin_hbca_c97_panCK <- FindClusters(Navin_hbca_c97_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c97_panCK <- RunUMAP(Navin_hbca_c97_panCK, dims = 1:40)
DimPlot(Navin_hbca_c97_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c97_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c97_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c97_panCK.rds")

# Remove Object
rm(Navin_hbca_c97_panCK)
gc()



############################ 
# Navin_hbca_c98_Singlets  #
############################ 
# Load Object
Navin_hbca_c98_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c98_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c98_panCK <- subset(x = Navin_hbca_c98_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c98_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c98_panCK <- FindVariableFeatures(Navin_hbca_c98_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c98_panCK)
Navin_hbca_c98_panCK <- ScaleData(Navin_hbca_c98_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c98_panCK <- RunPCA(Navin_hbca_c98_panCK, features = VariableFeatures(object = Navin_hbca_c98_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c98_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c98_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c98_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c98_panCK, ndims = 60)

Navin_hbca_c98_panCK <- FindNeighbors(Navin_hbca_c98_panCK, dims = 1:40)
Navin_hbca_c98_panCK <- FindClusters(Navin_hbca_c98_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c98_panCK <- RunUMAP(Navin_hbca_c98_panCK, dims = 1:40)
DimPlot(Navin_hbca_c98_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c98_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c98_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c98_panCK.rds")

# Remove Object
rm(Navin_hbca_c98_panCK)
gc()





############################ 
# Navin_hbca_c99_Singlets  #
############################ 
# Load Object
Navin_hbca_c99_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c99_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c99_panCK <- subset(x = Navin_hbca_c99_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c99_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c99_panCK <- FindVariableFeatures(Navin_hbca_c99_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c99_panCK)
Navin_hbca_c99_panCK <- ScaleData(Navin_hbca_c99_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c99_panCK <- RunPCA(Navin_hbca_c99_panCK, features = VariableFeatures(object = Navin_hbca_c99_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c99_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c99_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c99_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c99_panCK, ndims = 60)

Navin_hbca_c99_panCK <- FindNeighbors(Navin_hbca_c99_panCK, dims = 1:40)
Navin_hbca_c99_panCK <- FindClusters(Navin_hbca_c99_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c99_panCK <- RunUMAP(Navin_hbca_c99_panCK, dims = 1:40)
DimPlot(Navin_hbca_c99_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c99_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c99_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c99_panCK.rds")

# Remove Object
rm(Navin_hbca_c99_panCK)
gc()



############################ 
# Navin_hbca_c100_Singlets  #
############################ 
# Load Object
Navin_hbca_c100_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c100_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c100_panCK <- subset(x = Navin_hbca_c100_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c100_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c100_panCK <- FindVariableFeatures(Navin_hbca_c100_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c100_panCK)
Navin_hbca_c100_panCK <- ScaleData(Navin_hbca_c100_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c100_panCK <- RunPCA(Navin_hbca_c100_panCK, features = VariableFeatures(object = Navin_hbca_c100_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c100_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c100_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c100_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c100_panCK, ndims = 60)

Navin_hbca_c100_panCK <- FindNeighbors(Navin_hbca_c100_panCK, dims = 1:40)
Navin_hbca_c100_panCK <- FindClusters(Navin_hbca_c100_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c100_panCK <- RunUMAP(Navin_hbca_c100_panCK, dims = 1:40)
DimPlot(Navin_hbca_c100_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c100_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c100_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c100_panCK.rds")

# Remove Object
rm(Navin_hbca_c100_panCK)
gc()



############################ 
# Navin_hbca_c101_Singlets  #
############################ 
# Load Object
Navin_hbca_c101_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c101_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c101_panCK <- subset(x = Navin_hbca_c101_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c101_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c101_panCK <- FindVariableFeatures(Navin_hbca_c101_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c101_panCK)
Navin_hbca_c101_panCK <- ScaleData(Navin_hbca_c101_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c101_panCK <- RunPCA(Navin_hbca_c101_panCK, features = VariableFeatures(object = Navin_hbca_c101_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c101_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c101_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c101_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c101_panCK, ndims = 60)

Navin_hbca_c101_panCK <- FindNeighbors(Navin_hbca_c101_panCK, dims = 1:40)
Navin_hbca_c101_panCK <- FindClusters(Navin_hbca_c101_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c101_panCK <- RunUMAP(Navin_hbca_c101_panCK, dims = 1:40)
DimPlot(Navin_hbca_c101_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c101_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c101_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c101_panCK.rds")

# Remove Object
rm(Navin_hbca_c101_panCK)
gc()



############################ 
# Navin_hbca_c102_Singlets  #
############################ 
# Load Object
Navin_hbca_c102_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c102_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c102_panCK <- subset(x = Navin_hbca_c102_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c102_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c102_panCK <- FindVariableFeatures(Navin_hbca_c102_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c102_panCK)
Navin_hbca_c102_panCK <- ScaleData(Navin_hbca_c102_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c102_panCK <- RunPCA(Navin_hbca_c102_panCK, features = VariableFeatures(object = Navin_hbca_c102_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c102_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c102_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c102_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c102_panCK, ndims = 60)

Navin_hbca_c102_panCK <- FindNeighbors(Navin_hbca_c102_panCK, dims = 1:40)
Navin_hbca_c102_panCK <- FindClusters(Navin_hbca_c102_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c102_panCK <- RunUMAP(Navin_hbca_c102_panCK, dims = 1:40)
DimPlot(Navin_hbca_c102_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c102_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c102_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c102_panCK.rds")

# Remove Object
rm(Navin_hbca_c102_panCK)
gc()



############################ 
# Navin_hbca_c103_Singlets  #
############################ 
# Load Object
Navin_hbca_c103_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c103_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c103_panCK <- subset(x = Navin_hbca_c103_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c103_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c103_panCK <- FindVariableFeatures(Navin_hbca_c103_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c103_panCK)
Navin_hbca_c103_panCK <- ScaleData(Navin_hbca_c103_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c103_panCK <- RunPCA(Navin_hbca_c103_panCK, features = VariableFeatures(object = Navin_hbca_c103_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c103_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c103_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c103_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c103_panCK, ndims = 60)

Navin_hbca_c103_panCK <- FindNeighbors(Navin_hbca_c103_panCK, dims = 1:40)
Navin_hbca_c103_panCK <- FindClusters(Navin_hbca_c103_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c103_panCK <- RunUMAP(Navin_hbca_c103_panCK, dims = 1:40)
DimPlot(Navin_hbca_c103_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c103_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c103_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c103_panCK.rds")

# Remove Object
rm(Navin_hbca_c103_panCK)
gc()



############################ 
# Navin_hbca_c104_Singlets  #
############################ 
# Load Object
Navin_hbca_c104_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c104_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c104_panCK <- subset(x = Navin_hbca_c104_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c104_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c104_panCK <- FindVariableFeatures(Navin_hbca_c104_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c104_panCK)
Navin_hbca_c104_panCK <- ScaleData(Navin_hbca_c104_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c104_panCK <- RunPCA(Navin_hbca_c104_panCK, features = VariableFeatures(object = Navin_hbca_c104_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c104_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c104_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c104_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c104_panCK, ndims = 60)

Navin_hbca_c104_panCK <- FindNeighbors(Navin_hbca_c104_panCK, dims = 1:40)
Navin_hbca_c104_panCK <- FindClusters(Navin_hbca_c104_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c104_panCK <- RunUMAP(Navin_hbca_c104_panCK, dims = 1:40)
DimPlot(Navin_hbca_c104_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c104_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c104_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c104_panCK.rds")

# Remove Object
rm(Navin_hbca_c104_panCK)
gc()



############################ 
# Navin_hbca_c105_Singlets  #
############################ 
# Load Object
Navin_hbca_c105_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c105_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c105_panCK <- subset(x = Navin_hbca_c105_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c105_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c105_panCK <- FindVariableFeatures(Navin_hbca_c105_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c105_panCK)
Navin_hbca_c105_panCK <- ScaleData(Navin_hbca_c105_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c105_panCK <- RunPCA(Navin_hbca_c105_panCK, features = VariableFeatures(object = Navin_hbca_c105_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c105_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c105_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c105_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c105_panCK, ndims = 60)

Navin_hbca_c105_panCK <- FindNeighbors(Navin_hbca_c105_panCK, dims = 1:40)
Navin_hbca_c105_panCK <- FindClusters(Navin_hbca_c105_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c105_panCK <- RunUMAP(Navin_hbca_c105_panCK, dims = 1:40)
DimPlot(Navin_hbca_c105_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c105_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c105_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c105_panCK.rds")

# Remove Object
rm(Navin_hbca_c105_panCK)
gc()



############################ 
# Navin_hbca_c106_Singlets  #
############################ 
# Load Object
Navin_hbca_c106_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c106_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c106_panCK <- subset(x = Navin_hbca_c106_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c106_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c106_panCK <- FindVariableFeatures(Navin_hbca_c106_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c106_panCK)
Navin_hbca_c106_panCK <- ScaleData(Navin_hbca_c106_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c106_panCK <- RunPCA(Navin_hbca_c106_panCK, features = VariableFeatures(object = Navin_hbca_c106_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c106_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c106_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c106_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c106_panCK, ndims = 60)

Navin_hbca_c106_panCK <- FindNeighbors(Navin_hbca_c106_panCK, dims = 1:40)
Navin_hbca_c106_panCK <- FindClusters(Navin_hbca_c106_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c106_panCK <- RunUMAP(Navin_hbca_c106_panCK, dims = 1:40)
DimPlot(Navin_hbca_c106_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c106_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c106_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c106_panCK.rds")

# Remove Object
rm(Navin_hbca_c106_panCK)
gc()



############################ 
# Navin_hbca_c107_Singlets  #
############################ 
# Load Object
Navin_hbca_c107_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c107_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c107_panCK <- subset(x = Navin_hbca_c107_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c107_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c107_panCK <- FindVariableFeatures(Navin_hbca_c107_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c107_panCK)
Navin_hbca_c107_panCK <- ScaleData(Navin_hbca_c107_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c107_panCK <- RunPCA(Navin_hbca_c107_panCK, features = VariableFeatures(object = Navin_hbca_c107_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c107_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c107_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c107_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c107_panCK, ndims = 60)

Navin_hbca_c107_panCK <- FindNeighbors(Navin_hbca_c107_panCK, dims = 1:40)
Navin_hbca_c107_panCK <- FindClusters(Navin_hbca_c107_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c107_panCK <- RunUMAP(Navin_hbca_c107_panCK, dims = 1:40)
DimPlot(Navin_hbca_c107_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c107_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c107_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c107_panCK.rds")

# Remove Object
rm(Navin_hbca_c107_panCK)
gc()



############################ 
# Navin_hbca_c108_Singlets  #
############################ 
# Load Object
Navin_hbca_c108_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c108_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c108_panCK <- subset(x = Navin_hbca_c108_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c108_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c108_panCK <- FindVariableFeatures(Navin_hbca_c108_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c108_panCK)
Navin_hbca_c108_panCK <- ScaleData(Navin_hbca_c108_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c108_panCK <- RunPCA(Navin_hbca_c108_panCK, features = VariableFeatures(object = Navin_hbca_c108_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c108_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c108_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c108_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c108_panCK, ndims = 60)

Navin_hbca_c108_panCK <- FindNeighbors(Navin_hbca_c108_panCK, dims = 1:40)
Navin_hbca_c108_panCK <- FindClusters(Navin_hbca_c108_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c108_panCK <- RunUMAP(Navin_hbca_c108_panCK, dims = 1:40)
DimPlot(Navin_hbca_c108_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c108_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c108_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c108_panCK.rds")

# Remove Object
rm(Navin_hbca_c108_panCK)
gc()





############################ 
# Navin_hbca_c109_Singlets  #
############################ 
# Load Object
Navin_hbca_c109_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c109_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c109_panCK <- subset(x = Navin_hbca_c109_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c109_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c109_panCK <- FindVariableFeatures(Navin_hbca_c109_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c109_panCK)
Navin_hbca_c109_panCK <- ScaleData(Navin_hbca_c109_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c109_panCK <- RunPCA(Navin_hbca_c109_panCK, features = VariableFeatures(object = Navin_hbca_c109_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c109_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c109_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c109_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c109_panCK, ndims = 60)

Navin_hbca_c109_panCK <- FindNeighbors(Navin_hbca_c109_panCK, dims = 1:40)
Navin_hbca_c109_panCK <- FindClusters(Navin_hbca_c109_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c109_panCK <- RunUMAP(Navin_hbca_c109_panCK, dims = 1:40)
DimPlot(Navin_hbca_c109_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c109_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c109_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c109_panCK.rds")

# Remove Object
rm(Navin_hbca_c109_panCK)
gc()



############################ 
# Navin_hbca_c110_Singlets  #
############################ 
# Load Object
Navin_hbca_c110_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c110_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c110_panCK <- subset(x = Navin_hbca_c110_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c110_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c110_panCK <- FindVariableFeatures(Navin_hbca_c110_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c110_panCK)
Navin_hbca_c110_panCK <- ScaleData(Navin_hbca_c110_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c110_panCK <- RunPCA(Navin_hbca_c110_panCK, features = VariableFeatures(object = Navin_hbca_c110_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c110_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c110_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c110_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c110_panCK, ndims = 60)

Navin_hbca_c110_panCK <- FindNeighbors(Navin_hbca_c110_panCK, dims = 1:40)
Navin_hbca_c110_panCK <- FindClusters(Navin_hbca_c110_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c110_panCK <- RunUMAP(Navin_hbca_c110_panCK, dims = 1:40)
DimPlot(Navin_hbca_c110_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c110_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c110_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c110_panCK.rds")

# Remove Object
rm(Navin_hbca_c110_panCK)
gc()



############################ 
# Navin_hbca_c111_Singlets  #
############################ 
# Load Object
Navin_hbca_c111_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c111_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c111_panCK <- subset(x = Navin_hbca_c111_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c111_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c111_panCK <- FindVariableFeatures(Navin_hbca_c111_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c111_panCK)
Navin_hbca_c111_panCK <- ScaleData(Navin_hbca_c111_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c111_panCK <- RunPCA(Navin_hbca_c111_panCK, features = VariableFeatures(object = Navin_hbca_c111_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c111_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c111_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c111_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c111_panCK, ndims = 60)

Navin_hbca_c111_panCK <- FindNeighbors(Navin_hbca_c111_panCK, dims = 1:40)
Navin_hbca_c111_panCK <- FindClusters(Navin_hbca_c111_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c111_panCK <- RunUMAP(Navin_hbca_c111_panCK, dims = 1:40)
DimPlot(Navin_hbca_c111_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c111_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c111_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c111_panCK.rds")

# Remove Object
rm(Navin_hbca_c111_panCK)
gc()



############################ 
# Navin_hbca_c113_Singlets  #
############################ 
# Load Object
Navin_hbca_c113_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c113_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c113_panCK <- subset(x = Navin_hbca_c113_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c113_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c113_panCK <- FindVariableFeatures(Navin_hbca_c113_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c113_panCK)
Navin_hbca_c113_panCK <- ScaleData(Navin_hbca_c113_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c113_panCK <- RunPCA(Navin_hbca_c113_panCK, features = VariableFeatures(object = Navin_hbca_c113_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c113_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c113_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c113_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c113_panCK, ndims = 60)

Navin_hbca_c113_panCK <- FindNeighbors(Navin_hbca_c113_panCK, dims = 1:40)
Navin_hbca_c113_panCK <- FindClusters(Navin_hbca_c113_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c113_panCK <- RunUMAP(Navin_hbca_c113_panCK, dims = 1:40)
DimPlot(Navin_hbca_c113_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c113_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c113_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c113_panCK.rds")

# Remove Object
rm(Navin_hbca_c113_panCK)
gc()



############################ 
# Navin_hbca_c114_Singlets  #
############################ 
# Load Object
Navin_hbca_c114_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c114_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c114_panCK <- subset(x = Navin_hbca_c114_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c114_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c114_panCK <- FindVariableFeatures(Navin_hbca_c114_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c114_panCK)
Navin_hbca_c114_panCK <- ScaleData(Navin_hbca_c114_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c114_panCK <- RunPCA(Navin_hbca_c114_panCK, features = VariableFeatures(object = Navin_hbca_c114_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c114_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c114_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c114_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c114_panCK, ndims = 60)

Navin_hbca_c114_panCK <- FindNeighbors(Navin_hbca_c114_panCK, dims = 1:40)
Navin_hbca_c114_panCK <- FindClusters(Navin_hbca_c114_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c114_panCK <- RunUMAP(Navin_hbca_c114_panCK, dims = 1:40)
DimPlot(Navin_hbca_c114_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c114_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c114_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c114_panCK.rds")

# Remove Object
rm(Navin_hbca_c114_panCK)
gc()



############################ 
# Navin_hbca_c115_Singlets  #
############################ 
# Load Object
Navin_hbca_c115_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c115_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c115_panCK <- subset(x = Navin_hbca_c115_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c115_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c115_panCK <- FindVariableFeatures(Navin_hbca_c115_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c115_panCK)
Navin_hbca_c115_panCK <- ScaleData(Navin_hbca_c115_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c115_panCK <- RunPCA(Navin_hbca_c115_panCK, features = VariableFeatures(object = Navin_hbca_c115_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c115_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c115_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c115_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c115_panCK, ndims = 60)

Navin_hbca_c115_panCK <- FindNeighbors(Navin_hbca_c115_panCK, dims = 1:40)
Navin_hbca_c115_panCK <- FindClusters(Navin_hbca_c115_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c115_panCK <- RunUMAP(Navin_hbca_c115_panCK, dims = 1:40)
DimPlot(Navin_hbca_c115_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c115_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c115_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c115_panCK.rds")

# Remove Object
rm(Navin_hbca_c115_panCK)
gc()



############################ 
# Navin_hbca_c118_Singlets  #
############################ 
# Load Object
Navin_hbca_c118_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c118_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c118_panCK <- subset(x = Navin_hbca_c118_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c118_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c118_panCK <- FindVariableFeatures(Navin_hbca_c118_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c118_panCK)
Navin_hbca_c118_panCK <- ScaleData(Navin_hbca_c118_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c118_panCK <- RunPCA(Navin_hbca_c118_panCK, features = VariableFeatures(object = Navin_hbca_c118_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c118_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c118_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c118_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c118_panCK, ndims = 60)

Navin_hbca_c118_panCK <- FindNeighbors(Navin_hbca_c118_panCK, dims = 1:40)
Navin_hbca_c118_panCK <- FindClusters(Navin_hbca_c118_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c118_panCK <- RunUMAP(Navin_hbca_c118_panCK, dims = 1:40)
DimPlot(Navin_hbca_c118_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c118_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c118_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c118_panCK.rds")

# Remove Object
rm(Navin_hbca_c118_panCK)
gc()



############################ 
# Navin_hbca_c119_Singlets  #
############################ 
# Load Object
Navin_hbca_c119_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c119_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c119_panCK <- subset(x = Navin_hbca_c119_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c119_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c119_panCK <- FindVariableFeatures(Navin_hbca_c119_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c119_panCK)
Navin_hbca_c119_panCK <- ScaleData(Navin_hbca_c119_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c119_panCK <- RunPCA(Navin_hbca_c119_panCK, features = VariableFeatures(object = Navin_hbca_c119_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c119_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c119_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c119_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c119_panCK, ndims = 60)

Navin_hbca_c119_panCK <- FindNeighbors(Navin_hbca_c119_panCK, dims = 1:40)
Navin_hbca_c119_panCK <- FindClusters(Navin_hbca_c119_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c119_panCK <- RunUMAP(Navin_hbca_c119_panCK, dims = 1:40)
DimPlot(Navin_hbca_c119_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c119_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c119_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c119_panCK.rds")

# Remove Object
rm(Navin_hbca_c119_panCK)
gc()



############################ 
# Navin_hbca_c121_Singlets  #
############################ 
# Load Object
Navin_hbca_c121_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c121_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c121_panCK <- subset(x = Navin_hbca_c121_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c121_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c121_panCK <- FindVariableFeatures(Navin_hbca_c121_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c121_panCK)
Navin_hbca_c121_panCK <- ScaleData(Navin_hbca_c121_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c121_panCK <- RunPCA(Navin_hbca_c121_panCK, features = VariableFeatures(object = Navin_hbca_c121_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c121_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c121_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c121_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c121_panCK, ndims = 60)

Navin_hbca_c121_panCK <- FindNeighbors(Navin_hbca_c121_panCK, dims = 1:40)
Navin_hbca_c121_panCK <- FindClusters(Navin_hbca_c121_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c121_panCK <- RunUMAP(Navin_hbca_c121_panCK, dims = 1:40)
DimPlot(Navin_hbca_c121_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c121_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c121_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c121_panCK.rds")

# Remove Object
rm(Navin_hbca_c121_panCK)
gc()



############################ 
# Navin_hbca_c122_Singlets  #
############################ 
# Load Object
Navin_hbca_c122_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c122_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c122_panCK <- subset(x = Navin_hbca_c122_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c122_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c122_panCK <- FindVariableFeatures(Navin_hbca_c122_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c122_panCK)
Navin_hbca_c122_panCK <- ScaleData(Navin_hbca_c122_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c122_panCK <- RunPCA(Navin_hbca_c122_panCK, features = VariableFeatures(object = Navin_hbca_c122_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c122_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c122_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c122_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c122_panCK, ndims = 60)

Navin_hbca_c122_panCK <- FindNeighbors(Navin_hbca_c122_panCK, dims = 1:40)
Navin_hbca_c122_panCK <- FindClusters(Navin_hbca_c122_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c122_panCK <- RunUMAP(Navin_hbca_c122_panCK, dims = 1:40)
DimPlot(Navin_hbca_c122_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c122_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c122_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c122_panCK.rds")

# Remove Object
rm(Navin_hbca_c122_panCK)
gc()





############################ 
# Navin_hbca_c123_Singlets  #
############################ 
# Load Object
Navin_hbca_c123_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c123_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c123_panCK <- subset(x = Navin_hbca_c123_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c123_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c123_panCK <- FindVariableFeatures(Navin_hbca_c123_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c123_panCK)
Navin_hbca_c123_panCK <- ScaleData(Navin_hbca_c123_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c123_panCK <- RunPCA(Navin_hbca_c123_panCK, features = VariableFeatures(object = Navin_hbca_c123_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c123_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c123_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c123_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c123_panCK, ndims = 60)

Navin_hbca_c123_panCK <- FindNeighbors(Navin_hbca_c123_panCK, dims = 1:40)
Navin_hbca_c123_panCK <- FindClusters(Navin_hbca_c123_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c123_panCK <- RunUMAP(Navin_hbca_c123_panCK, dims = 1:40)
DimPlot(Navin_hbca_c123_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c123_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c123_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c123_panCK.rds")

# Remove Object
rm(Navin_hbca_c123_panCK)
gc()



############################ 
# Navin_hbca_c124_Singlets  #
############################ 
# Load Object
Navin_hbca_c124_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c124_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c124_panCK <- subset(x = Navin_hbca_c124_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c124_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c124_panCK <- FindVariableFeatures(Navin_hbca_c124_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c124_panCK)
Navin_hbca_c124_panCK <- ScaleData(Navin_hbca_c124_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c124_panCK <- RunPCA(Navin_hbca_c124_panCK, features = VariableFeatures(object = Navin_hbca_c124_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c124_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c124_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c124_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c124_panCK, ndims = 60)

Navin_hbca_c124_panCK <- FindNeighbors(Navin_hbca_c124_panCK, dims = 1:40)
Navin_hbca_c124_panCK <- FindClusters(Navin_hbca_c124_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c124_panCK <- RunUMAP(Navin_hbca_c124_panCK, dims = 1:40)
DimPlot(Navin_hbca_c124_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c124_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c124_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c124_panCK.rds")

# Remove Object
rm(Navin_hbca_c124_panCK)
gc()



############################ 
# Navin_hbca_c125_Singlets  #
############################ 
# Load Object
Navin_hbca_c125_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c125_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c125_panCK <- subset(x = Navin_hbca_c125_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c125_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c125_panCK <- FindVariableFeatures(Navin_hbca_c125_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c125_panCK)
Navin_hbca_c125_panCK <- ScaleData(Navin_hbca_c125_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c125_panCK <- RunPCA(Navin_hbca_c125_panCK, features = VariableFeatures(object = Navin_hbca_c125_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c125_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c125_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c125_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c125_panCK, ndims = 60)

Navin_hbca_c125_panCK <- FindNeighbors(Navin_hbca_c125_panCK, dims = 1:40)
Navin_hbca_c125_panCK <- FindClusters(Navin_hbca_c125_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c125_panCK <- RunUMAP(Navin_hbca_c125_panCK, dims = 1:40)
DimPlot(Navin_hbca_c125_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c125_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c125_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c125_panCK.rds")

# Remove Object
rm(Navin_hbca_c125_panCK)
gc()



############################ 
# Navin_hbca_c126_Singlets  #
############################ 
# Load Object
Navin_hbca_c126_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c126_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c126_panCK <- subset(x = Navin_hbca_c126_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c126_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c126_panCK <- FindVariableFeatures(Navin_hbca_c126_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c126_panCK)
Navin_hbca_c126_panCK <- ScaleData(Navin_hbca_c126_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c126_panCK <- RunPCA(Navin_hbca_c126_panCK, features = VariableFeatures(object = Navin_hbca_c126_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c126_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c126_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c126_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c126_panCK, ndims = 60)

Navin_hbca_c126_panCK <- FindNeighbors(Navin_hbca_c126_panCK, dims = 1:40)
Navin_hbca_c126_panCK <- FindClusters(Navin_hbca_c126_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c126_panCK <- RunUMAP(Navin_hbca_c126_panCK, dims = 1:40)
DimPlot(Navin_hbca_c126_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c126_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c126_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c126_panCK.rds")

# Remove Object
rm(Navin_hbca_c126_panCK)
gc()



############################ 
# Navin_hbca_c127_Singlets  #
############################ 
# Load Object
Navin_hbca_c127_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c127_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c127_panCK <- subset(x = Navin_hbca_c127_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c127_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c127_panCK <- FindVariableFeatures(Navin_hbca_c127_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c127_panCK)
Navin_hbca_c127_panCK <- ScaleData(Navin_hbca_c127_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c127_panCK <- RunPCA(Navin_hbca_c127_panCK, features = VariableFeatures(object = Navin_hbca_c127_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c127_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c127_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c127_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c127_panCK, ndims = 60)

Navin_hbca_c127_panCK <- FindNeighbors(Navin_hbca_c127_panCK, dims = 1:40)
Navin_hbca_c127_panCK <- FindClusters(Navin_hbca_c127_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c127_panCK <- RunUMAP(Navin_hbca_c127_panCK, dims = 1:40)
DimPlot(Navin_hbca_c127_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c127_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c127_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c127_panCK.rds")

# Remove Object
rm(Navin_hbca_c127_panCK)
gc()



############################ 
# Navin_hbca_c128_Singlets  #
############################ 
# Load Object
Navin_hbca_c128_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c128_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c128_panCK <- subset(x = Navin_hbca_c128_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c128_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c128_panCK <- FindVariableFeatures(Navin_hbca_c128_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c128_panCK)
Navin_hbca_c128_panCK <- ScaleData(Navin_hbca_c128_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c128_panCK <- RunPCA(Navin_hbca_c128_panCK, features = VariableFeatures(object = Navin_hbca_c128_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c128_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c128_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c128_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c128_panCK, ndims = 60)

Navin_hbca_c128_panCK <- FindNeighbors(Navin_hbca_c128_panCK, dims = 1:40)
Navin_hbca_c128_panCK <- FindClusters(Navin_hbca_c128_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c128_panCK <- RunUMAP(Navin_hbca_c128_panCK, dims = 1:40)
DimPlot(Navin_hbca_c128_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c128_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c128_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c128_panCK.rds")

# Remove Object
rm(Navin_hbca_c128_panCK)
gc()



############################ 
# Navin_hbca_c129_Singlets  #
############################ 
# Load Object
Navin_hbca_c129_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c129_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c129_panCK <- subset(x = Navin_hbca_c129_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c129_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c129_panCK <- FindVariableFeatures(Navin_hbca_c129_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c129_panCK)
Navin_hbca_c129_panCK <- ScaleData(Navin_hbca_c129_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c129_panCK <- RunPCA(Navin_hbca_c129_panCK, features = VariableFeatures(object = Navin_hbca_c129_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c129_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c129_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c129_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c129_panCK, ndims = 60)

Navin_hbca_c129_panCK <- FindNeighbors(Navin_hbca_c129_panCK, dims = 1:40)
Navin_hbca_c129_panCK <- FindClusters(Navin_hbca_c129_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c129_panCK <- RunUMAP(Navin_hbca_c129_panCK, dims = 1:40)
DimPlot(Navin_hbca_c129_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c129_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c129_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c129_panCK.rds")

# Remove Object
rm(Navin_hbca_c129_panCK)
gc()



############################ 
# Navin_hbca_c130_Singlets  #
############################ 
# Load Object
Navin_hbca_c130_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c130_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c130_panCK <- subset(x = Navin_hbca_c130_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c130_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c130_panCK <- FindVariableFeatures(Navin_hbca_c130_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c130_panCK)
Navin_hbca_c130_panCK <- ScaleData(Navin_hbca_c130_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c130_panCK <- RunPCA(Navin_hbca_c130_panCK, features = VariableFeatures(object = Navin_hbca_c130_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c130_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c130_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c130_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c130_panCK, ndims = 60)

Navin_hbca_c130_panCK <- FindNeighbors(Navin_hbca_c130_panCK, dims = 1:40)
Navin_hbca_c130_panCK <- FindClusters(Navin_hbca_c130_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c130_panCK <- RunUMAP(Navin_hbca_c130_panCK, dims = 1:40)
DimPlot(Navin_hbca_c130_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c130_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c130_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c130_panCK.rds")

# Remove Object
rm(Navin_hbca_c130_panCK)
gc()



############################ 
# Navin_hbca_c131_Singlets  #
############################ 
# Load Object
Navin_hbca_c131_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c131_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c131_panCK <- subset(x = Navin_hbca_c131_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c131_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c131_panCK <- FindVariableFeatures(Navin_hbca_c131_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c131_panCK)
Navin_hbca_c131_panCK <- ScaleData(Navin_hbca_c131_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c131_panCK <- RunPCA(Navin_hbca_c131_panCK, features = VariableFeatures(object = Navin_hbca_c131_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c131_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c131_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c131_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c131_panCK, ndims = 60)

Navin_hbca_c131_panCK <- FindNeighbors(Navin_hbca_c131_panCK, dims = 1:40)
Navin_hbca_c131_panCK <- FindClusters(Navin_hbca_c131_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c131_panCK <- RunUMAP(Navin_hbca_c131_panCK, dims = 1:40)
DimPlot(Navin_hbca_c131_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c131_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c131_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c131_panCK.rds")

# Remove Object
rm(Navin_hbca_c131_panCK)
gc()



############################ 
# Navin_hbca_c132_Singlets  #
############################ 
# Load Object
Navin_hbca_c132_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c132_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c132_panCK <- subset(x = Navin_hbca_c132_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c132_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c132_panCK <- FindVariableFeatures(Navin_hbca_c132_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c132_panCK)
Navin_hbca_c132_panCK <- ScaleData(Navin_hbca_c132_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c132_panCK <- RunPCA(Navin_hbca_c132_panCK, features = VariableFeatures(object = Navin_hbca_c132_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c132_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c132_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c132_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c132_panCK, ndims = 60)

Navin_hbca_c132_panCK <- FindNeighbors(Navin_hbca_c132_panCK, dims = 1:40)
Navin_hbca_c132_panCK <- FindClusters(Navin_hbca_c132_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c132_panCK <- RunUMAP(Navin_hbca_c132_panCK, dims = 1:40)
DimPlot(Navin_hbca_c132_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c132_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c132_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c132_panCK.rds")

# Remove Object
rm(Navin_hbca_c132_panCK)
gc()





############################ 
# Navin_hbca_c133_Singlets  #
############################ 
# Load Object
Navin_hbca_c133_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c133_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c133_panCK <- subset(x = Navin_hbca_c133_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c133_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c133_panCK <- FindVariableFeatures(Navin_hbca_c133_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c133_panCK)
Navin_hbca_c133_panCK <- ScaleData(Navin_hbca_c133_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c133_panCK <- RunPCA(Navin_hbca_c133_panCK, features = VariableFeatures(object = Navin_hbca_c133_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c133_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c133_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c133_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c133_panCK, ndims = 60)

Navin_hbca_c133_panCK <- FindNeighbors(Navin_hbca_c133_panCK, dims = 1:40)
Navin_hbca_c133_panCK <- FindClusters(Navin_hbca_c133_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c133_panCK <- RunUMAP(Navin_hbca_c133_panCK, dims = 1:40)
DimPlot(Navin_hbca_c133_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c133_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c133_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c133_panCK.rds")

# Remove Object
rm(Navin_hbca_c133_panCK)
gc()



############################ 
# Navin_hbca_c134_Singlets  #
############################ 
# Load Object
Navin_hbca_c134_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c134_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c134_panCK <- subset(x = Navin_hbca_c134_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c134_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c134_panCK <- FindVariableFeatures(Navin_hbca_c134_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c134_panCK)
Navin_hbca_c134_panCK <- ScaleData(Navin_hbca_c134_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c134_panCK <- RunPCA(Navin_hbca_c134_panCK, features = VariableFeatures(object = Navin_hbca_c134_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c134_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c134_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c134_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c134_panCK, ndims = 60)

Navin_hbca_c134_panCK <- FindNeighbors(Navin_hbca_c134_panCK, dims = 1:40)
Navin_hbca_c134_panCK <- FindClusters(Navin_hbca_c134_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c134_panCK <- RunUMAP(Navin_hbca_c134_panCK, dims = 1:40)
DimPlot(Navin_hbca_c134_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c134_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c134_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c134_panCK.rds")

# Remove Object
rm(Navin_hbca_c134_panCK)
gc()



############################ 
# Navin_hbca_c135_Singlets  #
############################ 
# Load Object
Navin_hbca_c135_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c135_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c135_panCK <- subset(x = Navin_hbca_c135_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c135_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c135_panCK <- FindVariableFeatures(Navin_hbca_c135_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c135_panCK)
Navin_hbca_c135_panCK <- ScaleData(Navin_hbca_c135_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c135_panCK <- RunPCA(Navin_hbca_c135_panCK, features = VariableFeatures(object = Navin_hbca_c135_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c135_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c135_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c135_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c135_panCK, ndims = 60)

Navin_hbca_c135_panCK <- FindNeighbors(Navin_hbca_c135_panCK, dims = 1:40)
Navin_hbca_c135_panCK <- FindClusters(Navin_hbca_c135_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c135_panCK <- RunUMAP(Navin_hbca_c135_panCK, dims = 1:40)
DimPlot(Navin_hbca_c135_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c135_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c135_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c135_panCK.rds")

# Remove Object
rm(Navin_hbca_c135_panCK)
gc()



############################ 
# Navin_hbca_c136_Singlets  #
############################ 
# Load Object
Navin_hbca_c136_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c136_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c136_panCK <- subset(x = Navin_hbca_c136_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c136_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c136_panCK <- FindVariableFeatures(Navin_hbca_c136_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c136_panCK)
Navin_hbca_c136_panCK <- ScaleData(Navin_hbca_c136_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c136_panCK <- RunPCA(Navin_hbca_c136_panCK, features = VariableFeatures(object = Navin_hbca_c136_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c136_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c136_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c136_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c136_panCK, ndims = 60)

Navin_hbca_c136_panCK <- FindNeighbors(Navin_hbca_c136_panCK, dims = 1:40)
Navin_hbca_c136_panCK <- FindClusters(Navin_hbca_c136_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c136_panCK <- RunUMAP(Navin_hbca_c136_panCK, dims = 1:40)
DimPlot(Navin_hbca_c136_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c136_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c136_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c136_panCK.rds")

# Remove Object
rm(Navin_hbca_c136_panCK)
gc()



############################ 
# Navin_hbca_c137_Singlets  #
############################ 
# Load Object
Navin_hbca_c137_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c137_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c137_panCK <- subset(x = Navin_hbca_c137_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c137_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c137_panCK <- FindVariableFeatures(Navin_hbca_c137_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c137_panCK)
Navin_hbca_c137_panCK <- ScaleData(Navin_hbca_c137_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c137_panCK <- RunPCA(Navin_hbca_c137_panCK, features = VariableFeatures(object = Navin_hbca_c137_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c137_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c137_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c137_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c137_panCK, ndims = 60)

Navin_hbca_c137_panCK <- FindNeighbors(Navin_hbca_c137_panCK, dims = 1:40)
Navin_hbca_c137_panCK <- FindClusters(Navin_hbca_c137_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c137_panCK <- RunUMAP(Navin_hbca_c137_panCK, dims = 1:40)
DimPlot(Navin_hbca_c137_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c137_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c137_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c137_panCK.rds")

# Remove Object
rm(Navin_hbca_c137_panCK)
gc()



############################ 
# Navin_hbca_c138_Singlets  #
############################ 
# Load Object
Navin_hbca_c138_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c138_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c138_panCK <- subset(x = Navin_hbca_c138_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c138_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c138_panCK <- FindVariableFeatures(Navin_hbca_c138_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c138_panCK)
Navin_hbca_c138_panCK <- ScaleData(Navin_hbca_c138_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c138_panCK <- RunPCA(Navin_hbca_c138_panCK, features = VariableFeatures(object = Navin_hbca_c138_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c138_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c138_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c138_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c138_panCK, ndims = 60)

Navin_hbca_c138_panCK <- FindNeighbors(Navin_hbca_c138_panCK, dims = 1:40)
Navin_hbca_c138_panCK <- FindClusters(Navin_hbca_c138_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c138_panCK <- RunUMAP(Navin_hbca_c138_panCK, dims = 1:40)
DimPlot(Navin_hbca_c138_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c138_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c138_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c138_panCK.rds")

# Remove Object
rm(Navin_hbca_c138_panCK)
gc()



############################ 
# Navin_hbca_c139_Singlets  #
############################ 
# Load Object
Navin_hbca_c139_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c139_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c139_panCK <- subset(x = Navin_hbca_c139_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c139_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c139_panCK <- FindVariableFeatures(Navin_hbca_c139_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c139_panCK)
Navin_hbca_c139_panCK <- ScaleData(Navin_hbca_c139_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c139_panCK <- RunPCA(Navin_hbca_c139_panCK, features = VariableFeatures(object = Navin_hbca_c139_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c139_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c139_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c139_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c139_panCK, ndims = 60)

Navin_hbca_c139_panCK <- FindNeighbors(Navin_hbca_c139_panCK, dims = 1:40)
Navin_hbca_c139_panCK <- FindClusters(Navin_hbca_c139_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c139_panCK <- RunUMAP(Navin_hbca_c139_panCK, dims = 1:40)
DimPlot(Navin_hbca_c139_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c139_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c139_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c139_panCK.rds")

# Remove Object
rm(Navin_hbca_c139_panCK)
gc()



############################ 
# Navin_hbca_c140_Singlets  #
############################ 
# Load Object
Navin_hbca_c140_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c140_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c140_panCK <- subset(x = Navin_hbca_c140_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c140_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c140_panCK <- FindVariableFeatures(Navin_hbca_c140_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c140_panCK)
Navin_hbca_c140_panCK <- ScaleData(Navin_hbca_c140_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c140_panCK <- RunPCA(Navin_hbca_c140_panCK, features = VariableFeatures(object = Navin_hbca_c140_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c140_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c140_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c140_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c140_panCK, ndims = 60)

Navin_hbca_c140_panCK <- FindNeighbors(Navin_hbca_c140_panCK, dims = 1:40)
Navin_hbca_c140_panCK <- FindClusters(Navin_hbca_c140_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c140_panCK <- RunUMAP(Navin_hbca_c140_panCK, dims = 1:40)
DimPlot(Navin_hbca_c140_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c140_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c140_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c140_panCK.rds")

# Remove Object
rm(Navin_hbca_c140_panCK)
gc()



############################ 
# Navin_hbca_c141_Singlets  #
############################ 
# Load Object
Navin_hbca_c141_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c141_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c141_panCK <- subset(x = Navin_hbca_c141_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c141_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c141_panCK <- FindVariableFeatures(Navin_hbca_c141_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c141_panCK)
Navin_hbca_c141_panCK <- ScaleData(Navin_hbca_c141_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c141_panCK <- RunPCA(Navin_hbca_c141_panCK, features = VariableFeatures(object = Navin_hbca_c141_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c141_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c141_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c141_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c141_panCK, ndims = 60)

Navin_hbca_c141_panCK <- FindNeighbors(Navin_hbca_c141_panCK, dims = 1:40)
Navin_hbca_c141_panCK <- FindClusters(Navin_hbca_c141_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c141_panCK <- RunUMAP(Navin_hbca_c141_panCK, dims = 1:40)
DimPlot(Navin_hbca_c141_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c141_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c141_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c141_panCK.rds")

# Remove Object
rm(Navin_hbca_c141_panCK)
gc()



############################ 
# Navin_hbca_c142_Singlets  #
############################ 
# Load Object
Navin_hbca_c142_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c142_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c142_panCK <- subset(x = Navin_hbca_c142_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c142_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c142_panCK <- FindVariableFeatures(Navin_hbca_c142_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c142_panCK)
Navin_hbca_c142_panCK <- ScaleData(Navin_hbca_c142_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c142_panCK <- RunPCA(Navin_hbca_c142_panCK, features = VariableFeatures(object = Navin_hbca_c142_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c142_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c142_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c142_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c142_panCK, ndims = 60)

Navin_hbca_c142_panCK <- FindNeighbors(Navin_hbca_c142_panCK, dims = 1:40)
Navin_hbca_c142_panCK <- FindClusters(Navin_hbca_c142_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c142_panCK <- RunUMAP(Navin_hbca_c142_panCK, dims = 1:40)
DimPlot(Navin_hbca_c142_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c142_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c142_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c142_panCK.rds")

# Remove Object
rm(Navin_hbca_c142_panCK)
gc()



############################ 
# Navin_hbca_c143_Singlets  #
############################ 
# Load Object
Navin_hbca_c143_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c143_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c143_panCK <- subset(x = Navin_hbca_c143_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c143_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c143_panCK <- FindVariableFeatures(Navin_hbca_c143_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c143_panCK)
Navin_hbca_c143_panCK <- ScaleData(Navin_hbca_c143_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c143_panCK <- RunPCA(Navin_hbca_c143_panCK, features = VariableFeatures(object = Navin_hbca_c143_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c143_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c143_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c143_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c143_panCK, ndims = 60)

Navin_hbca_c143_panCK <- FindNeighbors(Navin_hbca_c143_panCK, dims = 1:40)
Navin_hbca_c143_panCK <- FindClusters(Navin_hbca_c143_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c143_panCK <- RunUMAP(Navin_hbca_c143_panCK, dims = 1:40)
DimPlot(Navin_hbca_c143_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c143_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c143_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c143_panCK.rds")

# Remove Object
rm(Navin_hbca_c143_panCK)
gc()




############################ 
# Navin_hbca_c144_Singlets  #
############################ 
# Load Object
Navin_hbca_c144_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c144_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c144_panCK <- subset(x = Navin_hbca_c144_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c144_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c144_panCK <- FindVariableFeatures(Navin_hbca_c144_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c144_panCK)
Navin_hbca_c144_panCK <- ScaleData(Navin_hbca_c144_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c144_panCK <- RunPCA(Navin_hbca_c144_panCK, features = VariableFeatures(object = Navin_hbca_c144_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c144_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c144_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c144_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c144_panCK, ndims = 60)

Navin_hbca_c144_panCK <- FindNeighbors(Navin_hbca_c144_panCK, dims = 1:40)
Navin_hbca_c144_panCK <- FindClusters(Navin_hbca_c144_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c144_panCK <- RunUMAP(Navin_hbca_c144_panCK, dims = 1:40)
DimPlot(Navin_hbca_c144_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c144_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c144_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c144_panCK.rds")

# Remove Object
rm(Navin_hbca_c144_panCK)
gc()




############################ 
# Navin_hbca_c145_Singlets  #
############################ 
# Load Object
Navin_hbca_c145_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c145_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c145_panCK <- subset(x = Navin_hbca_c145_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c145_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c145_panCK <- FindVariableFeatures(Navin_hbca_c145_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c145_panCK)
Navin_hbca_c145_panCK <- ScaleData(Navin_hbca_c145_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c145_panCK <- RunPCA(Navin_hbca_c145_panCK, features = VariableFeatures(object = Navin_hbca_c145_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c145_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c145_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c145_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c145_panCK, ndims = 60)

Navin_hbca_c145_panCK <- FindNeighbors(Navin_hbca_c145_panCK, dims = 1:40)
Navin_hbca_c145_panCK <- FindClusters(Navin_hbca_c145_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c145_panCK <- RunUMAP(Navin_hbca_c145_panCK, dims = 1:40)
DimPlot(Navin_hbca_c145_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c145_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c145_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c145_panCK.rds")

# Remove Object
rm(Navin_hbca_c145_panCK)
gc()




############################ 
# Navin_hbca_c146_Singlets  #
############################ 
# Load Object
Navin_hbca_c146_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c146_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c146_panCK <- subset(x = Navin_hbca_c146_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c146_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c146_panCK <- FindVariableFeatures(Navin_hbca_c146_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c146_panCK)
Navin_hbca_c146_panCK <- ScaleData(Navin_hbca_c146_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c146_panCK <- RunPCA(Navin_hbca_c146_panCK, features = VariableFeatures(object = Navin_hbca_c146_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c146_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c146_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c146_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c146_panCK, ndims = 60)

Navin_hbca_c146_panCK <- FindNeighbors(Navin_hbca_c146_panCK, dims = 1:40)
Navin_hbca_c146_panCK <- FindClusters(Navin_hbca_c146_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c146_panCK <- RunUMAP(Navin_hbca_c146_panCK, dims = 1:40)
DimPlot(Navin_hbca_c146_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c146_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c146_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c146_panCK.rds")

# Remove Object
rm(Navin_hbca_c146_panCK)
gc()




############################ 
# Navin_hbca_c147_Singlets  #
############################ 
# Load Object
Navin_hbca_c147_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c147_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c147_panCK <- subset(x = Navin_hbca_c147_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c147_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c147_panCK <- FindVariableFeatures(Navin_hbca_c147_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c147_panCK)
Navin_hbca_c147_panCK <- ScaleData(Navin_hbca_c147_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c147_panCK <- RunPCA(Navin_hbca_c147_panCK, features = VariableFeatures(object = Navin_hbca_c147_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c147_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c147_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c147_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c147_panCK, ndims = 60)

Navin_hbca_c147_panCK <- FindNeighbors(Navin_hbca_c147_panCK, dims = 1:40)
Navin_hbca_c147_panCK <- FindClusters(Navin_hbca_c147_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c147_panCK <- RunUMAP(Navin_hbca_c147_panCK, dims = 1:40)
DimPlot(Navin_hbca_c147_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c147_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c147_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c147_panCK.rds")

# Remove Object
rm(Navin_hbca_c147_panCK)
gc()




############################ 
# Navin_hbca_c148_Singlets  #
############################ 
# Load Object
Navin_hbca_c148_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c148_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c148_panCK <- subset(x = Navin_hbca_c148_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c148_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c148_panCK <- FindVariableFeatures(Navin_hbca_c148_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c148_panCK)
Navin_hbca_c148_panCK <- ScaleData(Navin_hbca_c148_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c148_panCK <- RunPCA(Navin_hbca_c148_panCK, features = VariableFeatures(object = Navin_hbca_c148_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c148_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c148_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c148_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c148_panCK, ndims = 60)

Navin_hbca_c148_panCK <- FindNeighbors(Navin_hbca_c148_panCK, dims = 1:40)
Navin_hbca_c148_panCK <- FindClusters(Navin_hbca_c148_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c148_panCK <- RunUMAP(Navin_hbca_c148_panCK, dims = 1:40)
DimPlot(Navin_hbca_c148_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c148_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c148_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c148_panCK.rds")

# Remove Object
rm(Navin_hbca_c148_panCK)
gc()




############################ 
# Navin_hbca_c149_Singlets  #
############################ 
# Load Object
Navin_hbca_c149_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c149_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c149_panCK <- subset(x = Navin_hbca_c149_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c149_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c149_panCK <- FindVariableFeatures(Navin_hbca_c149_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c149_panCK)
Navin_hbca_c149_panCK <- ScaleData(Navin_hbca_c149_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c149_panCK <- RunPCA(Navin_hbca_c149_panCK, features = VariableFeatures(object = Navin_hbca_c149_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c149_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c149_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c149_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c149_panCK, ndims = 60)

Navin_hbca_c149_panCK <- FindNeighbors(Navin_hbca_c149_panCK, dims = 1:40)
Navin_hbca_c149_panCK <- FindClusters(Navin_hbca_c149_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c149_panCK <- RunUMAP(Navin_hbca_c149_panCK, dims = 1:40)
DimPlot(Navin_hbca_c149_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c149_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c149_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c149_panCK.rds")

# Remove Object
rm(Navin_hbca_c149_panCK)
gc()




############################ 
# Navin_hbca_c150_Singlets  #
############################ 
# Load Object
Navin_hbca_c150_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c150_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c150_panCK <- subset(x = Navin_hbca_c150_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c150_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c150_panCK <- FindVariableFeatures(Navin_hbca_c150_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c150_panCK)
Navin_hbca_c150_panCK <- ScaleData(Navin_hbca_c150_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c150_panCK <- RunPCA(Navin_hbca_c150_panCK, features = VariableFeatures(object = Navin_hbca_c150_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c150_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c150_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c150_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c150_panCK, ndims = 60)

Navin_hbca_c150_panCK <- FindNeighbors(Navin_hbca_c150_panCK, dims = 1:40)
Navin_hbca_c150_panCK <- FindClusters(Navin_hbca_c150_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c150_panCK <- RunUMAP(Navin_hbca_c150_panCK, dims = 1:40)
DimPlot(Navin_hbca_c150_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c150_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c150_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c150_panCK.rds")

# Remove Object
rm(Navin_hbca_c150_panCK)
gc()




############################ 
# Navin_hbca_c151_Singlets  #
############################ 
# Load Object
Navin_hbca_c151_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c151_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c151_panCK <- subset(x = Navin_hbca_c151_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c151_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c151_panCK <- FindVariableFeatures(Navin_hbca_c151_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c151_panCK)
Navin_hbca_c151_panCK <- ScaleData(Navin_hbca_c151_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c151_panCK <- RunPCA(Navin_hbca_c151_panCK, features = VariableFeatures(object = Navin_hbca_c151_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c151_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c151_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c151_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c151_panCK, ndims = 60)

Navin_hbca_c151_panCK <- FindNeighbors(Navin_hbca_c151_panCK, dims = 1:40)
Navin_hbca_c151_panCK <- FindClusters(Navin_hbca_c151_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c151_panCK <- RunUMAP(Navin_hbca_c151_panCK, dims = 1:40)
DimPlot(Navin_hbca_c151_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c151_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c151_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c151_panCK.rds")

# Remove Object
rm(Navin_hbca_c151_panCK)
gc()




############################ 
# Navin_hbca_c152_Singlets  #
############################ 
# Load Object
Navin_hbca_c152_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c152_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c152_panCK <- subset(x = Navin_hbca_c152_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c152_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c152_panCK <- FindVariableFeatures(Navin_hbca_c152_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c152_panCK)
Navin_hbca_c152_panCK <- ScaleData(Navin_hbca_c152_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c152_panCK <- RunPCA(Navin_hbca_c152_panCK, features = VariableFeatures(object = Navin_hbca_c152_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c152_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c152_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c152_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c152_panCK, ndims = 60)

Navin_hbca_c152_panCK <- FindNeighbors(Navin_hbca_c152_panCK, dims = 1:40)
Navin_hbca_c152_panCK <- FindClusters(Navin_hbca_c152_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c152_panCK <- RunUMAP(Navin_hbca_c152_panCK, dims = 1:40)
DimPlot(Navin_hbca_c152_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c152_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c152_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c152_panCK.rds")

# Remove Object
rm(Navin_hbca_c152_panCK)
gc()




############################ 
# Navin_hbca_c153_Singlets  #
############################ 
# Load Object
Navin_hbca_c153_Singlets <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Total_Singlets/Navin_hbca_c153_Singlets.rds")

# Collect cells greater than cutoffs for at least one of the three main mammary epithelial cell types
Navin_hbca_c153_panCK <- subset(x = Navin_hbca_c153_Singlets, subset = KRT14 > 1 | KRT18 > 1 | KRT19 > 1)

# Remove object
rm(Navin_hbca_c153_Singlets)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_hbca_c153_panCK <- FindVariableFeatures(Navin_hbca_c153_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_hbca_c153_panCK)
Navin_hbca_c153_panCK <- ScaleData(Navin_hbca_c153_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_hbca_c153_panCK <- RunPCA(Navin_hbca_c153_panCK, features = VariableFeatures(object = Navin_hbca_c153_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_hbca_c153_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_hbca_c153_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_hbca_c153_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_hbca_c153_panCK, ndims = 60)

Navin_hbca_c153_panCK <- FindNeighbors(Navin_hbca_c153_panCK, dims = 1:40)
Navin_hbca_c153_panCK <- FindClusters(Navin_hbca_c153_panCK, resolution = 0.15)

# Umap clustering
Navin_hbca_c153_panCK <- RunUMAP(Navin_hbca_c153_panCK, dims = 1:40)
DimPlot(Navin_hbca_c153_panCK, reduction = "umap", raster = FALSE)
DimPlot(Navin_hbca_c153_panCK, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(Navin_hbca_c153_panCK, file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c153_panCK.rds")

# Remove Object
rm(Navin_hbca_c153_panCK)
gc()



##############################################################
# Step 4: Merge Epithelial Cells into Combined Seurat Object #
##############################################################

# Merge in Batches of 20 samples
# Done in 6 initial batches and then 1 final batch

# Load Libraries 
library(Seurat)

########################
# Initial Batch 1 of 6 #
########################

#########################
# Load Individual Files #
#########################
Navin_hbca_c14_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c14_panCK.rds")
Navin_hbca_c15_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c15_panCK.rds")
Navin_hbca_c19_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c19_panCK.rds")
Navin_hbca_c20_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c20_panCK.rds")
Navin_hbca_c22_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c22_panCK.rds")
Navin_hbca_c23_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c23_panCK.rds")
Navin_hbca_c24_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c24_panCK.rds")
Navin_hbca_c25_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c25_panCK.rds")
Navin_hbca_c26_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c26_panCK.rds")
Navin_hbca_c31_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c31_panCK.rds")
Navin_hbca_c32_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c32_panCK.rds")
Navin_hbca_c50_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c50_panCK.rds")
Navin_hbca_c51_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c51_panCK.rds")
Navin_hbca_c52_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c52_panCK.rds")
Navin_hbca_c53_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c53_panCK.rds")
Navin_hbca_c54_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c54_panCK.rds")
Navin_hbca_c55_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c55_panCK.rds")
Navin_hbca_c56_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c56_panCK.rds")
Navin_hbca_c57_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c57_panCK.rds")
Navin_hbca_c58_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c58_panCK.rds")


##############################################
# Merge Individual Files into Single Dataset #
##############################################
# Merge Datasets
Navin_NORM_Batch1_panCK <- merge(x = Navin_hbca_c14_panCK, y = c(Navin_hbca_c15_panCK, Navin_hbca_c19_panCK, Navin_hbca_c20_panCK, 
                                                                 Navin_hbca_c22_panCK, Navin_hbca_c23_panCK, Navin_hbca_c24_panCK,
                                                                 Navin_hbca_c25_panCK, Navin_hbca_c26_panCK, Navin_hbca_c31_panCK,
                                                                 Navin_hbca_c32_panCK, Navin_hbca_c50_panCK, Navin_hbca_c51_panCK,
                                                                 Navin_hbca_c52_panCK, Navin_hbca_c53_panCK, Navin_hbca_c54_panCK,
                                                                 Navin_hbca_c55_panCK, Navin_hbca_c56_panCK, Navin_hbca_c57_panCK,
                                                                 Navin_hbca_c58_panCK))

# Join Layers
Navin_NORM_Batch1_panCK <- JoinLayers(Navin_NORM_Batch1_panCK)

# Remove Objects
rm(Navin_hbca_c14_panCK)
rm(Navin_hbca_c15_panCK)
rm(Navin_hbca_c19_panCK)
rm(Navin_hbca_c20_panCK)
rm(Navin_hbca_c22_panCK)
rm(Navin_hbca_c23_panCK)
rm(Navin_hbca_c24_panCK)
rm(Navin_hbca_c25_panCK)
rm(Navin_hbca_c26_panCK)
rm(Navin_hbca_c31_panCK)
rm(Navin_hbca_c32_panCK)
rm(Navin_hbca_c50_panCK)
rm(Navin_hbca_c51_panCK)
rm(Navin_hbca_c52_panCK)
rm(Navin_hbca_c53_panCK)
rm(Navin_hbca_c54_panCK)
rm(Navin_hbca_c55_panCK)
rm(Navin_hbca_c56_panCK)
rm(Navin_hbca_c57_panCK)
rm(Navin_hbca_c58_panCK)
gc()


####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_NORM_Batch1_panCK <- FindVariableFeatures(Navin_NORM_Batch1_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_NORM_Batch1_panCK)
Navin_NORM_Batch1_panCK <- ScaleData(Navin_NORM_Batch1_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_NORM_Batch1_panCK <- RunPCA(Navin_NORM_Batch1_panCK, features = VariableFeatures(object = Navin_NORM_Batch1_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_NORM_Batch1_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_NORM_Batch1_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_NORM_Batch1_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_NORM_Batch1_panCK)

# Choose dimensions
Navin_NORM_Batch1_panCK <- FindNeighbors(Navin_NORM_Batch1_panCK, dims = 1:40)
Navin_NORM_Batch1_panCK <- FindClusters(Navin_NORM_Batch1_panCK, resolution = 0.12)

# Umap clustering
Navin_NORM_Batch1_panCK <- RunUMAP(Navin_NORM_Batch1_panCK, dims = 1:40)
DimPlot(Navin_NORM_Batch1_panCK, reduction = "umap")

# Identify Samples
DimPlot(Navin_NORM_Batch1_panCK, reduction = "umap", group.by = "orig.ident")

# Save
saveRDS(Navin_NORM_Batch1_panCK, file = "/R/R_Navin/Navin_RDS/RDS_Merged/Navin_NORM_Batch1_panCK.rds")

# Clear global environment
rm(Navin_NORM_Batch1_panCK)
gc()





########################
# Initial Batch 2 of 6 #
########################

#########################
# Load Individual Files #
#########################
Navin_hbca_c59_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c59_panCK.rds")
Navin_hbca_c60_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c60_panCK.rds")
Navin_hbca_c61_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c61_panCK.rds")
Navin_hbca_c62_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c62_panCK.rds")
Navin_hbca_c63_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c63_panCK.rds")
Navin_hbca_c64_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c64_panCK.rds")
Navin_hbca_c65_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c65_panCK.rds")
Navin_hbca_c66_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c66_panCK.rds")
Navin_hbca_c67_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c67_panCK.rds")
Navin_hbca_c68_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c68_panCK.rds")
Navin_hbca_c69_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c69_panCK.rds")
Navin_hbca_c70_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c70_panCK.rds")
Navin_hbca_c71_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c71_panCK.rds")
Navin_hbca_c72_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c72_panCK.rds")
Navin_hbca_c73_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c73_panCK.rds")
Navin_hbca_c74_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c74_panCK.rds")
Navin_hbca_c75_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c75_panCK.rds")
Navin_hbca_c76_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c76_panCK.rds")
Navin_hbca_c77_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c77_panCK.rds")
Navin_hbca_c78_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c78_panCK.rds")


##############################################
# Merge Individual Files into Single Dataset #
##############################################
# Merge Datasets
Navin_NORM_Batch2_panCK <- merge(x = Navin_hbca_c59_panCK, y = c(Navin_hbca_c60_panCK, Navin_hbca_c61_panCK, Navin_hbca_c62_panCK, 
                                                                 Navin_hbca_c63_panCK, Navin_hbca_c64_panCK, Navin_hbca_c65_panCK,
                                                                 Navin_hbca_c66_panCK, Navin_hbca_c67_panCK, Navin_hbca_c68_panCK,
                                                                 Navin_hbca_c69_panCK, Navin_hbca_c70_panCK, Navin_hbca_c71_panCK,
                                                                 Navin_hbca_c72_panCK, Navin_hbca_c73_panCK, Navin_hbca_c74_panCK,
                                                                 Navin_hbca_c75_panCK, Navin_hbca_c76_panCK, Navin_hbca_c77_panCK,
                                                                 Navin_hbca_c78_panCK))

# Join Layers
Navin_NORM_Batch2_panCK <- JoinLayers(Navin_NORM_Batch2_panCK)

# Remove Objects
rm(Navin_hbca_c59_panCK)
rm(Navin_hbca_c60_panCK)
rm(Navin_hbca_c61_panCK)
rm(Navin_hbca_c62_panCK)
rm(Navin_hbca_c63_panCK)
rm(Navin_hbca_c64_panCK)
rm(Navin_hbca_c65_panCK)
rm(Navin_hbca_c66_panCK)
rm(Navin_hbca_c67_panCK)
rm(Navin_hbca_c68_panCK)
rm(Navin_hbca_c69_panCK)
rm(Navin_hbca_c70_panCK)
rm(Navin_hbca_c71_panCK)
rm(Navin_hbca_c72_panCK)
rm(Navin_hbca_c73_panCK)
rm(Navin_hbca_c74_panCK)
rm(Navin_hbca_c75_panCK)
rm(Navin_hbca_c76_panCK)
rm(Navin_hbca_c77_panCK)
rm(Navin_hbca_c78_panCK)
gc()


####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_NORM_Batch2_panCK <- FindVariableFeatures(Navin_NORM_Batch2_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_NORM_Batch2_panCK)
Navin_NORM_Batch2_panCK <- ScaleData(Navin_NORM_Batch2_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_NORM_Batch2_panCK <- RunPCA(Navin_NORM_Batch2_panCK, features = VariableFeatures(object = Navin_NORM_Batch2_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_NORM_Batch2_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_NORM_Batch2_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_NORM_Batch2_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_NORM_Batch2_panCK)

# Choose dimensions
Navin_NORM_Batch2_panCK <- FindNeighbors(Navin_NORM_Batch2_panCK, dims = 1:40)
Navin_NORM_Batch2_panCK <- FindClusters(Navin_NORM_Batch2_panCK, resolution = 0.12)

# Umap clustering
Navin_NORM_Batch2_panCK <- RunUMAP(Navin_NORM_Batch2_panCK, dims = 1:40)
DimPlot(Navin_NORM_Batch2_panCK, reduction = "umap")

# Identify Samples
DimPlot(Navin_NORM_Batch2_panCK, reduction = "umap", group.by = "orig.ident")

# Save
saveRDS(Navin_NORM_Batch2_panCK, file = "/R/R_Navin/Navin_RDS/RDS_Merged/Navin_NORM_Batch2_panCK.rds")

# Clear global environment
rm(Navin_NORM_Batch2_panCK)
gc()







########################
# Initial Batch 3 of 6 #
########################

#########################
# Load Individual Files #
#########################
Navin_hbca_c79_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c79_panCK.rds")
Navin_hbca_c80_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c80_panCK.rds")
Navin_hbca_c81_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c81_panCK.rds")
Navin_hbca_c82_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c82_panCK.rds")
Navin_hbca_c83_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c83_panCK.rds")
Navin_hbca_c84_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c84_panCK.rds")
Navin_hbca_c85_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c85_panCK.rds")
Navin_hbca_c86_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c86_panCK.rds")
Navin_hbca_c87_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c87_panCK.rds")
Navin_hbca_c88_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c88_panCK.rds")
Navin_hbca_c89_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c89_panCK.rds")
Navin_hbca_c90_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c90_panCK.rds")
Navin_hbca_c91_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c91_panCK.rds")
Navin_hbca_c92_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c92_panCK.rds")
Navin_hbca_c93_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c93_panCK.rds")
Navin_hbca_c94_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c94_panCK.rds")
Navin_hbca_c95_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c95_panCK.rds")
Navin_hbca_c96_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c96_panCK.rds")
Navin_hbca_c97_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c97_panCK.rds")
Navin_hbca_c98_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c98_panCK.rds")


##############################################
# Merge Individual Files into Single Dataset #
##############################################
# Merge Datasets
Navin_NORM_Batch3_panCK <- merge(x = Navin_hbca_c79_panCK, y = c(Navin_hbca_c80_panCK, Navin_hbca_c81_panCK, Navin_hbca_c82_panCK, 
                                                                 Navin_hbca_c83_panCK, Navin_hbca_c84_panCK, Navin_hbca_c85_panCK,
                                                                 Navin_hbca_c86_panCK, Navin_hbca_c87_panCK, Navin_hbca_c88_panCK,
                                                                 Navin_hbca_c89_panCK, Navin_hbca_c90_panCK, Navin_hbca_c91_panCK,
                                                                 Navin_hbca_c92_panCK, Navin_hbca_c93_panCK, Navin_hbca_c94_panCK,
                                                                 Navin_hbca_c95_panCK, Navin_hbca_c96_panCK, Navin_hbca_c97_panCK,
                                                                 Navin_hbca_c98_panCK))

# Join Layers
Navin_NORM_Batch3_panCK <- JoinLayers(Navin_NORM_Batch3_panCK)

# Remove Objects
rm(Navin_hbca_c79_panCK)
rm(Navin_hbca_c80_panCK)
rm(Navin_hbca_c81_panCK)
rm(Navin_hbca_c82_panCK)
rm(Navin_hbca_c83_panCK)
rm(Navin_hbca_c84_panCK)
rm(Navin_hbca_c85_panCK)
rm(Navin_hbca_c86_panCK)
rm(Navin_hbca_c87_panCK)
rm(Navin_hbca_c88_panCK)
rm(Navin_hbca_c89_panCK)
rm(Navin_hbca_c90_panCK)
rm(Navin_hbca_c91_panCK)
rm(Navin_hbca_c92_panCK)
rm(Navin_hbca_c93_panCK)
rm(Navin_hbca_c94_panCK)
rm(Navin_hbca_c95_panCK)
rm(Navin_hbca_c96_panCK)
rm(Navin_hbca_c97_panCK)
rm(Navin_hbca_c98_panCK)
gc()


####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_NORM_Batch3_panCK <- FindVariableFeatures(Navin_NORM_Batch3_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_NORM_Batch3_panCK)
Navin_NORM_Batch3_panCK <- ScaleData(Navin_NORM_Batch3_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_NORM_Batch3_panCK <- RunPCA(Navin_NORM_Batch3_panCK, features = VariableFeatures(object = Navin_NORM_Batch3_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_NORM_Batch3_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_NORM_Batch3_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_NORM_Batch3_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_NORM_Batch3_panCK)

# Choose dimensions
Navin_NORM_Batch3_panCK <- FindNeighbors(Navin_NORM_Batch3_panCK, dims = 1:40)
Navin_NORM_Batch3_panCK <- FindClusters(Navin_NORM_Batch3_panCK, resolution = 0.12)

# Umap clustering
Navin_NORM_Batch3_panCK <- RunUMAP(Navin_NORM_Batch3_panCK, dims = 1:40)
DimPlot(Navin_NORM_Batch3_panCK, reduction = "umap")

# Identify Samples
DimPlot(Navin_NORM_Batch3_panCK, reduction = "umap", group.by = "orig.ident")

# Save
saveRDS(Navin_NORM_Batch3_panCK, file = "/R/R_Navin/Navin_RDS/RDS_Merged/Navin_NORM_Batch3_panCK.rds")

# Clear global environment
rm(Navin_NORM_Batch3_panCK)
gc()







########################
# Initial Batch 4 of 6 #
########################

#########################
# Load Individual Files #
#########################
Navin_hbca_c99_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c99_panCK.rds")
Navin_hbca_c100_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c100_panCK.rds")
Navin_hbca_c101_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c101_panCK.rds")
Navin_hbca_c102_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c102_panCK.rds")
Navin_hbca_c103_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c103_panCK.rds")
Navin_hbca_c104_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c104_panCK.rds")
Navin_hbca_c105_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c105_panCK.rds")
Navin_hbca_c106_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c106_panCK.rds")
Navin_hbca_c107_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c107_panCK.rds")
Navin_hbca_c108_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c108_panCK.rds")
Navin_hbca_c109_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c109_panCK.rds")
Navin_hbca_c110_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c110_panCK.rds")
Navin_hbca_c111_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c111_panCK.rds")
Navin_hbca_c113_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c113_panCK.rds")
Navin_hbca_c114_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c114_panCK.rds")
Navin_hbca_c115_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c115_panCK.rds")
Navin_hbca_c118_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c118_panCK.rds")
Navin_hbca_c119_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c119_panCK.rds")
Navin_hbca_c121_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c121_panCK.rds")
Navin_hbca_c122_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c122_panCK.rds")


##############################################
# Merge Individual Files into Single Dataset #
##############################################
# Merge Datasets
Navin_NORM_Batch4_panCK <- merge(x = Navin_hbca_c99_panCK, y = c(Navin_hbca_c100_panCK, Navin_hbca_c101_panCK, Navin_hbca_c102_panCK, 
                                                                 Navin_hbca_c103_panCK, Navin_hbca_c104_panCK, Navin_hbca_c105_panCK,
                                                                 Navin_hbca_c106_panCK, Navin_hbca_c107_panCK, Navin_hbca_c108_panCK,
                                                                 Navin_hbca_c109_panCK, Navin_hbca_c110_panCK, Navin_hbca_c111_panCK,
                                                                 Navin_hbca_c113_panCK, Navin_hbca_c114_panCK, Navin_hbca_c115_panCK,
                                                                 Navin_hbca_c118_panCK, Navin_hbca_c119_panCK, Navin_hbca_c121_panCK,
                                                                 Navin_hbca_c122_panCK))

# Join Layers
Navin_NORM_Batch4_panCK <- JoinLayers(Navin_NORM_Batch4_panCK)

# Remove Objects
rm(Navin_hbca_c99_panCK)
rm(Navin_hbca_c100_panCK)
rm(Navin_hbca_c101_panCK)
rm(Navin_hbca_c102_panCK)
rm(Navin_hbca_c103_panCK)
rm(Navin_hbca_c104_panCK)
rm(Navin_hbca_c105_panCK)
rm(Navin_hbca_c106_panCK)
rm(Navin_hbca_c107_panCK)
rm(Navin_hbca_c108_panCK)
rm(Navin_hbca_c109_panCK)
rm(Navin_hbca_c110_panCK)
rm(Navin_hbca_c111_panCK)
rm(Navin_hbca_c113_panCK)
rm(Navin_hbca_c114_panCK)
rm(Navin_hbca_c115_panCK)
rm(Navin_hbca_c118_panCK)
rm(Navin_hbca_c119_panCK)
rm(Navin_hbca_c121_panCK)
rm(Navin_hbca_c122_panCK)
gc()


####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_NORM_Batch4_panCK <- FindVariableFeatures(Navin_NORM_Batch4_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_NORM_Batch4_panCK)
Navin_NORM_Batch4_panCK <- ScaleData(Navin_NORM_Batch4_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_NORM_Batch4_panCK <- RunPCA(Navin_NORM_Batch4_panCK, features = VariableFeatures(object = Navin_NORM_Batch4_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_NORM_Batch4_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_NORM_Batch4_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_NORM_Batch4_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_NORM_Batch4_panCK)

# Choose dimensions
Navin_NORM_Batch4_panCK <- FindNeighbors(Navin_NORM_Batch4_panCK, dims = 1:40)
Navin_NORM_Batch4_panCK <- FindClusters(Navin_NORM_Batch4_panCK, resolution = 0.12)

# Umap clustering
Navin_NORM_Batch4_panCK <- RunUMAP(Navin_NORM_Batch4_panCK, dims = 1:40)
DimPlot(Navin_NORM_Batch4_panCK, reduction = "umap")

# Identify Samples
DimPlot(Navin_NORM_Batch4_panCK, reduction = "umap", group.by = "orig.ident")

# Save
saveRDS(Navin_NORM_Batch4_panCK, file = "/R/R_Navin/Navin_RDS/RDS_Merged/Navin_NORM_Batch4_panCK.rds")

# Clear global environment
rm(Navin_NORM_Batch4_panCK)
gc()






########################
# Initial Batch 5 of 6 #
########################

#########################
# Load Individual Files #
#########################
Navin_hbca_c123_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c123_panCK.rds")
Navin_hbca_c124_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c124_panCK.rds")
Navin_hbca_c125_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c125_panCK.rds")
Navin_hbca_c126_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c126_panCK.rds")
Navin_hbca_c127_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c127_panCK.rds")
Navin_hbca_c128_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c128_panCK.rds")
Navin_hbca_c129_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c129_panCK.rds")
Navin_hbca_c130_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c130_panCK.rds")
Navin_hbca_c131_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c131_panCK.rds")
Navin_hbca_c132_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c132_panCK.rds")
Navin_hbca_c133_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c133_panCK.rds")
Navin_hbca_c134_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c134_panCK.rds")
Navin_hbca_c135_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c135_panCK.rds")
Navin_hbca_c136_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c136_panCK.rds")
Navin_hbca_c137_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c137_panCK.rds")
Navin_hbca_c138_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c138_panCK.rds")
Navin_hbca_c139_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c139_panCK.rds")
Navin_hbca_c140_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c140_panCK.rds")
Navin_hbca_c141_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c141_panCK.rds")
Navin_hbca_c142_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c142_panCK.rds")


##############################################
# Merge Individual Files into Single Dataset #
##############################################
# Merge Datasets
Navin_NORM_Batch5_panCK <- merge(x = Navin_hbca_c123_panCK, y = c(Navin_hbca_c124_panCK, Navin_hbca_c125_panCK, Navin_hbca_c126_panCK, 
                                                                  Navin_hbca_c127_panCK, Navin_hbca_c128_panCK, Navin_hbca_c129_panCK, 
                                                                  Navin_hbca_c130_panCK, Navin_hbca_c131_panCK, Navin_hbca_c132_panCK, 
                                                                  Navin_hbca_c133_panCK, Navin_hbca_c134_panCK, Navin_hbca_c135_panCK, 
                                                                  Navin_hbca_c136_panCK, Navin_hbca_c137_panCK,  Navin_hbca_c138_panCK,
                                                                  Navin_hbca_c139_panCK, Navin_hbca_c140_panCK,  Navin_hbca_c141_panCK,
                                                                  Navin_hbca_c142_panCK))


# Join Layers
Navin_NORM_Batch5_panCK <- JoinLayers(Navin_NORM_Batch5_panCK)

# Remove Objects
rm(Navin_hbca_c123_panCK)
rm(Navin_hbca_c124_panCK)
rm(Navin_hbca_c125_panCK)
rm(Navin_hbca_c126_panCK)
rm(Navin_hbca_c127_panCK)
rm(Navin_hbca_c128_panCK)
rm(Navin_hbca_c129_panCK)
rm(Navin_hbca_c130_panCK)
rm(Navin_hbca_c131_panCK)
rm(Navin_hbca_c132_panCK)
rm(Navin_hbca_c133_panCK)
rm(Navin_hbca_c134_panCK)
rm(Navin_hbca_c135_panCK)
rm(Navin_hbca_c136_panCK)
rm(Navin_hbca_c137_panCK)
rm(Navin_hbca_c138_panCK)
rm(Navin_hbca_c139_panCK)
rm(Navin_hbca_c140_panCK)
rm(Navin_hbca_c141_panCK)
rm(Navin_hbca_c142_panCK)
gc()


####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_NORM_Batch5_panCK <- FindVariableFeatures(Navin_NORM_Batch5_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_NORM_Batch5_panCK)
Navin_NORM_Batch5_panCK <- ScaleData(Navin_NORM_Batch5_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_NORM_Batch5_panCK <- RunPCA(Navin_NORM_Batch5_panCK, features = VariableFeatures(object = Navin_NORM_Batch5_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_NORM_Batch5_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_NORM_Batch5_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_NORM_Batch5_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_NORM_Batch5_panCK)

# Choose dimensions
Navin_NORM_Batch5_panCK <- FindNeighbors(Navin_NORM_Batch5_panCK, dims = 1:40)
Navin_NORM_Batch5_panCK <- FindClusters(Navin_NORM_Batch5_panCK, resolution = 0.12)

# Umap clustering
Navin_NORM_Batch5_panCK <- RunUMAP(Navin_NORM_Batch5_panCK, dims = 1:40)
DimPlot(Navin_NORM_Batch5_panCK, reduction = "umap")

# Identify Samples
DimPlot(Navin_NORM_Batch5_panCK, reduction = "umap", group.by = "orig.ident")

# Save
saveRDS(Navin_NORM_Batch5_panCK, file = "/R/R_Navin/Navin_RDS/RDS_Merged/Navin_NORM_Batch5_panCK.rds")

# Clear global environment
rm(Navin_NORM_Batch5_panCK)
gc()







########################
# Initial Batch 6 of 6 #
########################

#########################
# Load Individual Files #
#########################
Navin_hbca_c143_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c143_panCK.rds")
Navin_hbca_c144_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c144_panCK.rds")
Navin_hbca_c145_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c145_panCK.rds")
Navin_hbca_c146_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c146_panCK.rds")
Navin_hbca_c147_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c147_panCK.rds")
Navin_hbca_c148_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c148_panCK.rds")
Navin_hbca_c149_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c149_panCK.rds")
Navin_hbca_c150_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c150_panCK.rds")
Navin_hbca_c151_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c151_panCK.rds")
Navin_hbca_c152_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c152_panCK.rds")
Navin_hbca_c153_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c153_panCK.rds")


##############################################
# Merge Individual Files into Single Dataset #
##############################################
# Merge Datasets
Navin_NORM_Batch6_panCK <- merge(x = Navin_hbca_c143_panCK, y = c(Navin_hbca_c144_panCK, Navin_hbca_c145_panCK, Navin_hbca_c146_panCK, 
                                                                  Navin_hbca_c147_panCK, Navin_hbca_c148_panCK, Navin_hbca_c149_panCK, 
                                                                  Navin_hbca_c150_panCK, Navin_hbca_c151_panCK, Navin_hbca_c152_panCK, 
                                                                  Navin_hbca_c153_panCK))


# Join Layers
Navin_NORM_Batch6_panCK <- JoinLayers(Navin_NORM_Batch6_panCK)

# Remove Objects
rm(Navin_hbca_c143_panCK)
rm(Navin_hbca_c144_panCK)
rm(Navin_hbca_c145_panCK)
rm(Navin_hbca_c146_panCK)
rm(Navin_hbca_c147_panCK)
rm(Navin_hbca_c148_panCK)
rm(Navin_hbca_c149_panCK)
rm(Navin_hbca_c150_panCK)
rm(Navin_hbca_c151_panCK)
rm(Navin_hbca_c152_panCK)
rm(Navin_hbca_c153_panCK)
gc()


####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_NORM_Batch6_panCK <- FindVariableFeatures(Navin_NORM_Batch6_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_NORM_Batch6_panCK)
Navin_NORM_Batch6_panCK <- ScaleData(Navin_NORM_Batch6_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_NORM_Batch6_panCK <- RunPCA(Navin_NORM_Batch6_panCK, features = VariableFeatures(object = Navin_NORM_Batch6_panCK))

# Examine and visualize PCA results a few different ways
print(Navin_NORM_Batch6_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_NORM_Batch6_panCK, dims = 1:2, reduction = "pca")
DimPlot(Navin_NORM_Batch6_panCK, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_NORM_Batch6_panCK)

# Choose dimensions
Navin_NORM_Batch6_panCK <- FindNeighbors(Navin_NORM_Batch6_panCK, dims = 1:40)
Navin_NORM_Batch6_panCK <- FindClusters(Navin_NORM_Batch6_panCK, resolution = 0.12)

# Umap clustering
Navin_NORM_Batch6_panCK <- RunUMAP(Navin_NORM_Batch6_panCK, dims = 1:40)
DimPlot(Navin_NORM_Batch6_panCK, reduction = "umap")

# Identify Samples
DimPlot(Navin_NORM_Batch6_panCK, reduction = "umap", group.by = "orig.ident")

# Save
saveRDS(Navin_NORM_Batch6_panCK, file = "/R/R_Navin/Navin_RDS/RDS_Merged/Navin_NORM_Batch6_panCK.rds")

# Clear global environment
rm(Navin_NORM_Batch6_panCK)
gc()




###############
# Final Merge #
###############

#########################
# Load Individual Files #
#########################
Navin_NORM_Batch1_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Merged/Navin_NORM_Batch1_panCK.rds")
Navin_NORM_Batch2_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Merged/Navin_NORM_Batch2_panCK.rds")
Navin_NORM_Batch3_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Merged/Navin_NORM_Batch3_panCK.rds")
Navin_NORM_Batch4_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Merged/Navin_NORM_Batch4_panCK.rds")
Navin_NORM_Batch5_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Merged/Navin_NORM_Batch5_panCK.rds")
Navin_NORM_Batch6_panCK <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Merged/Navin_NORM_Batch6_panCK.rds")



##############################################
# Merge Individual Files into Single Dataset #
##############################################
# Merge Datasets
Navin_NORM_panCK_merged <- merge(x = Navin_NORM_Batch1_panCK, y = c(Navin_NORM_Batch2_panCK, Navin_NORM_Batch3_panCK, Navin_NORM_Batch4_panCK, 
                                                                    Navin_NORM_Batch5_panCK, Navin_NORM_Batch6_panCK))


# Join Layers
Navin_NORM_panCK_merged <- JoinLayers(Navin_NORM_panCK_merged)

# Remove Objects
rm(Navin_NORM_Batch1_panCK)
rm(Navin_NORM_Batch2_panCK)
rm(Navin_NORM_Batch3_panCK)
rm(Navin_NORM_Batch4_panCK)
rm(Navin_NORM_Batch5_panCK)
rm(Navin_NORM_Batch6_panCK)
gc()


####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_NORM_panCK_merged <- FindVariableFeatures(Navin_NORM_panCK_merged, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_NORM_panCK_merged)
Navin_NORM_panCK_merged <- ScaleData(Navin_NORM_panCK_merged, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_NORM_panCK_merged <- RunPCA(Navin_NORM_panCK_merged, features = VariableFeatures(object = Navin_NORM_panCK_merged))

# Examine and visualize PCA results a few different ways
print(Navin_NORM_panCK_merged[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_NORM_panCK_merged, dims = 1:2, reduction = "pca")
DimPlot(Navin_NORM_panCK_merged, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_NORM_panCK_merged)

# Choose dimensions
Navin_NORM_panCK_merged <- FindNeighbors(Navin_NORM_panCK_merged, dims = 1:40)
Navin_NORM_panCK_merged <- FindClusters(Navin_NORM_panCK_merged, resolution = 0.12)

# Umap clustering
Navin_NORM_panCK_merged <- RunUMAP(Navin_NORM_panCK_merged, dims = 1:40)
DimPlot(Navin_NORM_panCK_merged, reduction = "umap")

# Identify Samples
DimPlot(Navin_NORM_panCK_merged, reduction = "umap", group.by = "orig.ident")

# Save
saveRDS(Navin_NORM_panCK_merged, file = "/R/R_Navin/Navin_RDS/RDS_Merged/Navin_NORM_panCK_merged.rds")

# Clear global environment
rm(Navin_NORM_panCK_merged)
gc()
