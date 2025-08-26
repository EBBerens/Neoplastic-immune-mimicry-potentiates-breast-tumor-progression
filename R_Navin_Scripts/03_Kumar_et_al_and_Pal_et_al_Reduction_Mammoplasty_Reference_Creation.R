#####################################################################################################################
#        Navin (Kumar et al) and Visvader (Pal et al) Datasets Reduction Mammoplasty Analysis Steps                 #
#####################################################################################################################
# Step 5: Merge Navin and Visvader Reduction Mammoplasty Seurat Objects                                             #
# Step 6: Create Reduction Mammoplasty Reference for inferCNV                                                       #
#####################################################################################################################

# Load Libraries
library(Seurat)
library(patchwork)
library(dplyr)

#########################################################################
# Step 5: Merge Navin and Visvader Reduction Mammoplasty Seurat Objects #
#########################################################################


#########################
# Load Individual Files #
#########################
Navin_NORM_panCK_merged <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Merged/Navin_NORM_panCK_merged.rds")
Visvader_NORM_panCK_merged <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Merged/Visvader_NORM_panCK_merged.rds")


##############################################
# Merge Individual Files into Single Dataset #
##############################################
# Merge Datasets
Navin_Visvader_NORM_panCK_combined <- merge(x = Navin_NORM_panCK_merged, y = c(Visvader_NORM_panCK_merged))


# Join Layers
Navin_Visvader_NORM_panCK_combined <- JoinLayers(Navin_Visvader_NORM_panCK_combined)

# Remove Objects
rm(Navin_NORM_panCK_merged)
rm(Visvader_NORM_panCK_merged)
gc()


####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_Visvader_NORM_panCK_combined <- FindVariableFeatures(Navin_Visvader_NORM_panCK_combined, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_Visvader_NORM_panCK_combined)
Navin_Visvader_NORM_panCK_combined <- ScaleData(Navin_Visvader_NORM_panCK_combined, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_Visvader_NORM_panCK_combined <- RunPCA(Navin_Visvader_NORM_panCK_combined, features = VariableFeatures(object = Navin_Visvader_NORM_panCK_combined))

# Examine and visualize PCA results a few different ways
print(Navin_Visvader_NORM_panCK_combined[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_Visvader_NORM_panCK_combined, dims = 1:2, reduction = "pca")
DimPlot(Navin_Visvader_NORM_panCK_combined, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_Visvader_NORM_panCK_combined)

# Choose dimensions
Navin_Visvader_NORM_panCK_combined <- FindNeighbors(Navin_Visvader_NORM_panCK_combined, dims = 1:40)
Navin_Visvader_NORM_panCK_combined <- FindClusters(Navin_Visvader_NORM_panCK_combined, resolution = 0.08)

# Umap clustering
Navin_Visvader_NORM_panCK_combined <- RunUMAP(Navin_Visvader_NORM_panCK_combined, dims = 1:40)
DimPlot(Navin_Visvader_NORM_panCK_combined, reduction = "umap", label = TRUE, raster = FALSE)

# Identify Samples
DimPlot(Navin_Visvader_NORM_panCK_combined, reduction = "umap", group.by = "orig.ident")

# Save
saveRDS(Navin_Visvader_NORM_panCK_combined, file = "/R/R_Navin/Navin_RDS/RDS_Merged/Navin_Visvader_NORM_panCK_combined.rds")

# Clear global environment
#rm(Navin_Visvader_NORM_panCK_combined)
gc()





################################################################
# Step 6: Create Reduction Mammoplasty Reference for inferCNV  #                                                       
################################################################

# Note: 8 clusters are detected after reference creation; 8 reference groups will thus be used during inferCNV 

# Load Navin & Visvader normal mammary epithelial SeuratObject
#Navin_Visvader_NORM_panCK_combined <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Merged/Navin_Visvader_NORM_panCK_combined.rds")

# Set seed for reproducibility
set.seed(123)

# Downsample cells from each epithelial sample using orig.ident to create reference
# Will target 100 cells from each sample for samplling, yielding a combined reference of ~12,000 cells
downsampled_cells <- Navin_Visvader_NORM_panCK_combined@meta.data %>%
  rownames_to_column(var = "cell") %>%
  group_by(orig.ident) %>%
  sample_n(size = min(100, n()), replace = FALSE) %>%
  pull(cell)

# Create subsetted SeuratObject with the downsampled cells
Navin_Visvader_NORM_panCK_Reference <- subset(Navin_Visvader_NORM_panCK_combined, cells = downsampled_cells)

# Removed starting object and downsampled intermediate
rm(Navin_Visvader_NORM_panCK_combined)
rm(downsampled_cells)

####################################
#  Reprocess the downsampled data  #
####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Navin_Visvader_NORM_panCK_Reference <- FindVariableFeatures(Navin_Visvader_NORM_panCK_Reference, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Navin_Visvader_NORM_panCK_Reference)
Navin_Visvader_NORM_panCK_Reference <- ScaleData(Navin_Visvader_NORM_panCK_Reference, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Navin_Visvader_NORM_panCK_Reference <- RunPCA(Navin_Visvader_NORM_panCK_Reference, features = VariableFeatures(object = Navin_Visvader_NORM_panCK_Reference))

# Examine and visualize PCA results a few different ways
print(Navin_Visvader_NORM_panCK_Reference[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Navin_Visvader_NORM_panCK_Reference, dims = 1:2, reduction = "pca")
DimPlot(Navin_Visvader_NORM_panCK_Reference, reduction = "pca")

# Visualize PCs
ElbowPlot(Navin_Visvader_NORM_panCK_Reference, ndims = 50)

# Choose dimensions
Navin_Visvader_NORM_panCK_Reference <- FindNeighbors(Navin_Visvader_NORM_panCK_Reference, dims = 1:30)
Navin_Visvader_NORM_panCK_Reference <- FindClusters(Navin_Visvader_NORM_panCK_Reference, resolution = 0.11)

# Umap clustering
Navin_Visvader_NORM_panCK_Reference <- RunUMAP(Navin_Visvader_NORM_panCK_Reference, dims = 1:30)
DimPlot(Navin_Visvader_NORM_panCK_Reference, reduction = "umap")

# Identify Samples
DimPlot(Navin_Visvader_NORM_panCK_Reference, reduction = "umap", group.by = "orig.ident")

# Save
saveRDS(Navin_Visvader_NORM_panCK_Reference, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Navin_Visvader_NORM_panCK_Reference.rds")

# Remove Objects
rm(Navin_Visvader_NORM_panCK_Reference)
gc()
