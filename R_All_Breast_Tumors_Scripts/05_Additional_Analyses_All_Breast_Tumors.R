#######################################################################################################################
#                          Additional Analyses with All Breast Tumor Samples Merged                                   #
#######################################################################################################################
# Step 19: Correlate Immune Mimicry Receptor Gene Expression Across Human Breast Tumors                               #
# Step 20: Recombine Annotated Stroma with Annotated Neoplastic Compartment for Each Dataset                          #
# Step 21: Merged All Annotated Cells from All Breast Tumors                                                          #  
# Step 22: Identify Differentially Expressed Genes in Immune-mimicked vs Neoplastic or Immune-mimicked vs Leukocytes  #
# Step 23: Determine Biological Processes Related to Cell Type DEGs                                                   #
# Step 24: Quantify Breast Tumor Cells Expressing IM Markers Inside vs Outside IM Clusters                            #
# Step 25: Evaluate Breast Cancer Stem Gene Expression per Immune Mimicry Marker & Generate Heatmap                   #
#######################################################################################################################

#########################################################################################
# Step 19: Correlate Immune Mimicry Receptor Gene Expression Across Human Breast Tumors #
#########################################################################################

# Load Libraries
library(Seurat)
library(dplyr)
library(corrplot)
library(ggplot2)

#################
# Preprocessing #
#################
##############################################################
# Extract Immune Mimicry Expression Matrices for Merged Data #
##############################################################

########################################
# all_Visvader_breast_tumors_panCK_cnv #
########################################
# Load data file
all_Visvader_breast_tumors_panCK_cnv <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Merged/all_Visvader_breast_tumors_panCK_cnv.rds")
# Extract Expression Matrix
all_Visvader_breast_tumors_panCK_cnv <- all_Visvader_breast_tumors_panCK_cnv[["RNA"]]$data
all_Visvader_breast_tumors_panCK_cnv <- as.matrix(all_Visvader_breast_tumors_panCK_cnv, 'sparseMatrix')
all_Visvader_breast_tumors_panCK_cnv <- subset(all_Visvader_breast_tumors_panCK_cnv, rownames(all_Visvader_breast_tumors_panCK_cnv) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
all_Visvader_breast_tumors_panCK_cnv <-t(all_Visvader_breast_tumors_panCK_cnv)
write.csv(all_Visvader_breast_tumors_panCK_cnv, file = "/R/R_Visvader/Visvader_Expression_Matrices/Visvader_Expression_Matrices_Output/all_Visvader_breast_tumors_panCK_cnv_Immune_Mimicry_Cell_Expression_Matrix_Normalized_Data_for_Correlations.csv")
gc()


# Load data file
all_Swarbrick_breast_tumors_panCK_cnv <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Merged/all_Swarbrick_breast_tumors_panCK_cnv.rds")
# Extract Expression Matrix
all_Swarbrick_breast_tumors_panCK_cnv <- all_Swarbrick_breast_tumors_panCK_cnv[["RNA"]]$data
all_Swarbrick_breast_tumors_panCK_cnv <- as.matrix(all_Swarbrick_breast_tumors_panCK_cnv, 'sparseMatrix')
all_Swarbrick_breast_tumors_panCK_cnv <- subset(all_Swarbrick_breast_tumors_panCK_cnv, rownames(all_Swarbrick_breast_tumors_panCK_cnv) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
all_Swarbrick_breast_tumors_panCK_cnv <-t(all_Swarbrick_breast_tumors_panCK_cnv)
write.csv(all_Swarbrick_breast_tumors_panCK_cnv, file = "/R/R_Swarbrick/Swarbrick_Expression_Matrices/Swarbrick_Expression_Matrices_Output/all_Swarbrick_breast_tumors_panCK_cnv_Immune_Mimicry_Cell_Expression_Matrix_Normalized_Data_for_Correlations.csv")
gc()


# Load data file
all_PANNTHR_breast_tumors_panCK_cnv <- readRDS(file = "/R/R_PANNTHR/PANNTHR_RDS/RDS_Merged/all_PANNTHR_breast_tumors_panCK_cnv.rds")
# Extract Expression Matrix
all_PANNTHR_breast_tumors_panCK_cnv <- all_PANNTHR_breast_tumors_panCK_cnv[["RNA"]]$data
all_PANNTHR_breast_tumors_panCK_cnv <- as.matrix(all_PANNTHR_breast_tumors_panCK_cnv, 'sparseMatrix')
all_PANNTHR_breast_tumors_panCK_cnv <- subset(all_PANNTHR_breast_tumors_panCK_cnv, rownames(all_PANNTHR_breast_tumors_panCK_cnv) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
all_PANNTHR_breast_tumors_panCK_cnv <-t(all_PANNTHR_breast_tumors_panCK_cnv)
write.csv(all_PANNTHR_breast_tumors_panCK_cnv, file = "/R/R_PANNTHR/PANNTHR_Expression_Matrices/PANNTHR_Expression_Matrices_Output/all_PANNTHR_breast_tumors_panCK_cnv_Immune_Mimicry_Cell_Expression_Matrix_Normalized_Data_for_Correlations.csv")
gc()

###################################
# Concatenate Expression Matrices #
###################################
immune_mimicry_correlations <- rbind(all_Visvader_breast_tumors_panCK_cnv, all_Swarbrick_breast_tumors_panCK_cnv, all_PANNTHR_breast_tumors_panCK_cnv)
rm(all_Visvader_breast_tumors_panCK_cnv)
rm(all_Swarbrick_breast_tumors_panCK_cnv)
rm(all_PANNTHR_breast_tumors_panCK_cnv)

############################################
# Perform Correlation Analysis in native R #
############################################
# http://www.sthda.com/english/wiki/correlation-matrix-a-quick-start-guide-to-analyze-format-and-visualize-a-correlation-matrix-using-r-software
cor.immune_mimicry_correlations<- cor(immune_mimicry_correlations, method=c("spearman"))

# Save Results
write.csv(cor.immune_mimicry_correlations, file = "/R/R_Visvader/Visvader_Expression_Matrices/Immune_mimicry_Marker_Correlation_Matrix_Visvader_Swarbrick_PANNTHR.csv")

##############################################################################################################

#################################
# Make Plots with Custom Colors #
#################################
# Set colors (e.g. reverse color palette)
col<- colorRampPalette(c("#2166ac", "#67a9cf", "#f7f7f7", "#ef8a62", "#b2182b"))(200)


# Compute p-values
# Will set significance threshold of 0.001 below
testRes = cor.mtest(cor.immune_mimicry_correlations)

# Make Plot without Coefficients
corrplot(cor.immune_mimicry_correlations, p.mat = testRes$p, sig.level = 0.001, order = 'hclust', hclust.method = "average", insig = 'pch', pch.cex = 0, method = "shade", shade.col = NA, col = col)
dev.copy(tiff, '/R/R_Visvader/Visvader_Output/Visvader_Swarbrick_IM_Correlations_custom_color.tiff', width=5000, height=5000, res=320)
dev.off()

# Make Plot with Coefficients
corrplot(cor.immune_mimicry_correlations, p.mat = testRes$p, sig.level = 0.001, order = 'hclust', hclust.method = "average", insig = 'pch', pch.cex = 0, addCoef.col = 'black', method = "shade", shade.col = NA,  col = col)
dev.copy(tiff, '/R/R_Visvader/Visvader_Output/Visvader_Swarbrick_IM_Correlations_with_coef_custom_color.tiff', width=5000, height=5000, res=320)
dev.off()

# Make Plot without Coefficients
corrplot(cor.immune_mimicry_correlations, p.mat = testRes$p, sig.level = 0.001, order = 'hclust', hclust.method = "average", insig = 'blank', pch.cex = 0, method = "shade", shade.col = NA,  col = col)
dev.copy(tiff, '/R/R_Visvader/Visvader_Output/Visvader_Swarbrick_IM_Correlations_only_sig_custom_color.tiff', width=5000, height=5000, res=320)
dev.off()

# Make Custom plot
corrplot(cor.immune_mimicry_correlations, p.mat = testRes$p, sig.level = 0.001, order = 'hclust', hclust.method = "average", insig = 'blank', addCoef.col ='black', number.cex = 0.8, pch.cex = 0, method = "shade", shade.col = NA,  col = col)
dev.copy(tiff, '/R/R_Visvader/Visvader_Output/Visvader_Swarbrick_IM_Correlations_only_sig_Custom_custom_color.tiff', width=5000, height=5000, res=320)
dev.off()

# Remove objects
rm(cor.immune_mimicry_correlations)
rm(testRes)
gc()



##############################################################################################
# Step 20: Recombine Annotated Stroma with Annotated Neoplastic Compartment for Each Dataset #
##############################################################################################

##################################
######## Visvader Dataset ########
##################################
# Load Annotated SeuratObjects
all_Visvader_breast_tumors_stroma_annotated <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Annotated/all_Visvader_breast_tumors_stroma_annotated.rds")
all_Visvader_breast_tumors_panCK_cnv_annotated <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Annotated/all_Visvader_breast_tumors_panCK_cnv_annotated.rds")


###############################################################################
# Initial annotations should be retained under cell.ident; update annotations #
###############################################################################

# Set Next Identity
all_Visvader_breast_tumors_stroma_annotated <- SetIdent(all_Visvader_breast_tumors_stroma_annotated, value = "cell.ident")

# Check Idents
Idents(object = all_Visvader_breast_tumors_stroma_annotated)
table(Idents(all_Visvader_breast_tumors_stroma_annotated))

# Annotate Clusters = main.ident
Stromal_IDs_main <- c("Leukocytes", "Leukocytes", "Mesenchymal Cells", "Leukocytes", "Mesenchymal Cells", "Mesenchymal Cells", "Leukocytes", "Leukocytes")

# Apply New Labels to Clusters
names(Stromal_IDs_main) <- levels(all_Visvader_breast_tumors_stroma_annotated)
all_Visvader_breast_tumors_stroma_annotated <- RenameIdents(all_Visvader_breast_tumors_stroma_annotated, Stromal_IDs_main)

# Remove old Seurat Object and ID list
rm(Stromal_IDs_main)

# Store Annotations Under Classification
all_Visvader_breast_tumors_stroma_annotated[["main.ident"]] <- Idents(object = all_Visvader_breast_tumors_stroma_annotated)


##############################################
# Merge Individual Files into Single Dataset #
##############################################
# Merge Datasets
Recombined_Visvader_Stroma_Neoplastic_annotated <- merge(all_Visvader_breast_tumors_stroma_annotated, y = c(all_Visvader_breast_tumors_panCK_cnv_annotated))

# Remove starting objects
rm(all_Visvader_breast_tumors_stroma_annotated)
rm(all_Visvader_breast_tumors_panCK_cnv_annotated)
gc()

# Join Layers
Recombined_Visvader_Stroma_Neoplastic_annotated <- JoinLayers(Recombined_Visvader_Stroma_Neoplastic_annotated)


####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Recombined_Visvader_Stroma_Neoplastic_annotated <- FindVariableFeatures(Recombined_Visvader_Stroma_Neoplastic_annotated, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Recombined_Visvader_Stroma_Neoplastic_annotated)
Recombined_Visvader_Stroma_Neoplastic_annotated <- ScaleData(Recombined_Visvader_Stroma_Neoplastic_annotated, features = all.genes)


#####################
# REDUCE DIMENSIONS #
#####################
Recombined_Visvader_Stroma_Neoplastic_annotated <- RunPCA(Recombined_Visvader_Stroma_Neoplastic_annotated, features = VariableFeatures(object = Recombined_Visvader_Stroma_Neoplastic_annotated))

# Examine and visualize PCA results a few different ways
print(Recombined_Visvader_Stroma_Neoplastic_annotated[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Recombined_Visvader_Stroma_Neoplastic_annotated, dims = 1:2, reduction = "pca")
DimPlot(Recombined_Visvader_Stroma_Neoplastic_annotated, reduction = "pca")

# Visualize PCs
ElbowPlot(Recombined_Visvader_Stroma_Neoplastic_annotated, ndims = 50)

Recombined_Visvader_Stroma_Neoplastic_annotated <- FindNeighbors(Recombined_Visvader_Stroma_Neoplastic_annotated, dims = 1:50)
Recombined_Visvader_Stroma_Neoplastic_annotated <- FindClusters(Recombined_Visvader_Stroma_Neoplastic_annotated, resolution = 0.12)

# Umap clustering
Recombined_Visvader_Stroma_Neoplastic_annotated <- RunUMAP(Recombined_Visvader_Stroma_Neoplastic_annotated, dims = 1:50)
DimPlot(Recombined_Visvader_Stroma_Neoplastic_annotated, reduction = "umap", raster = FALSE, group.by = "detailed.ident")
DimPlot(Recombined_Visvader_Stroma_Neoplastic_annotated, reduction = "umap", group.by = "orig.ident", raster = FALSE)

########
# SAVE #
########
saveRDS(Recombined_Visvader_Stroma_Neoplastic_annotated, file = "/R/R_Visvader/Visvader_RDS/RDS_Annotated/Recombined_Visvader_Stroma_Neoplastic_annotated.rds")

# Remove Object
rm(Recombined_Visvader_Stroma_Neoplastic_annotated)
gc()



###################################
######## Swarbrick Dataset ########
###################################
# Load Annotated SeuratObjects
all_Swarbrick_breast_tumors_stroma_annotated <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Annotated/all_Swarbrick_breast_tumors_stroma_annotated.rds")
all_Swarbrick_breast_tumors_panCK_cnv_annotated <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Annotated/all_Swarbrick_breast_tumors_panCK_cnv_annotated.rds")


###############################################################################
# Initial annotations should be retained under cell.ident; update annotations #
###############################################################################

# Set Next Identity
all_Swarbrick_breast_tumors_stroma_annotated <- SetIdent(all_Swarbrick_breast_tumors_stroma_annotated, value = "cell.ident")

# Check Idents
Idents(object = all_Swarbrick_breast_tumors_stroma_annotated)
table(Idents(all_Swarbrick_breast_tumors_stroma_annotated))

# Annotate Clusters = main.ident
Stromal_IDs_main <- c("Leukocytes", "Leukocytes", "Mesenchymal Cells", "Mesenchymal Cells", "Mesenchymal Cells", "Leukocytes", "Leukocytes", "Mesenchymal Cells", "Leukocytes")

# Apply New Labels to Clusters
names(Stromal_IDs_main) <- levels(all_Swarbrick_breast_tumors_stroma_annotated)
all_Swarbrick_breast_tumors_stroma_annotated <- RenameIdents(all_Swarbrick_breast_tumors_stroma_annotated, Stromal_IDs_main)

# Remove old Seurat Object and ID list
rm(Stromal_IDs_main)

# Store Annotations Under Classification
all_Swarbrick_breast_tumors_stroma_annotated[["main.ident"]] <- Idents(object = all_Swarbrick_breast_tumors_stroma_annotated)


##############################################
# Merge Individual Files into Single Dataset #
##############################################
# Merge Datasets
Recombined_Swarbrick_Stroma_Neoplastic_annotated <- merge(all_Swarbrick_breast_tumors_stroma_annotated, y = c(all_Swarbrick_breast_tumors_panCK_cnv_annotated))

# Remove starting objects
rm(all_Swarbrick_breast_tumors_stroma_annotated)
rm(all_Swarbrick_breast_tumors_panCK_cnv_annotated)
gc()

# Join Layers
Recombined_Swarbrick_Stroma_Neoplastic_annotated <- JoinLayers(Recombined_Swarbrick_Stroma_Neoplastic_annotated)


####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Recombined_Swarbrick_Stroma_Neoplastic_annotated <- FindVariableFeatures(Recombined_Swarbrick_Stroma_Neoplastic_annotated, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Recombined_Swarbrick_Stroma_Neoplastic_annotated)
Recombined_Swarbrick_Stroma_Neoplastic_annotated <- ScaleData(Recombined_Swarbrick_Stroma_Neoplastic_annotated, features = all.genes)


#####################
# REDUCE DIMENSIONS #
#####################
Recombined_Swarbrick_Stroma_Neoplastic_annotated <- RunPCA(Recombined_Swarbrick_Stroma_Neoplastic_annotated, features = VariableFeatures(object = Recombined_Swarbrick_Stroma_Neoplastic_annotated))

# Examine and visualize PCA results a few different ways
print(Recombined_Swarbrick_Stroma_Neoplastic_annotated[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Recombined_Swarbrick_Stroma_Neoplastic_annotated, dims = 1:2, reduction = "pca")
DimPlot(Recombined_Swarbrick_Stroma_Neoplastic_annotated, reduction = "pca")

# Visualize PCs
ElbowPlot(Recombined_Swarbrick_Stroma_Neoplastic_annotated, ndims = 50)

Recombined_Swarbrick_Stroma_Neoplastic_annotated <- FindNeighbors(Recombined_Swarbrick_Stroma_Neoplastic_annotated, dims = 1:50)
Recombined_Swarbrick_Stroma_Neoplastic_annotated <- FindClusters(Recombined_Swarbrick_Stroma_Neoplastic_annotated, resolution = 0.12)

# Umap clustering
Recombined_Swarbrick_Stroma_Neoplastic_annotated <- RunUMAP(Recombined_Swarbrick_Stroma_Neoplastic_annotated, dims = 1:50)
DimPlot(Recombined_Swarbrick_Stroma_Neoplastic_annotated, reduction = "umap", raster = FALSE, group.by = "detailed.ident")
DimPlot(Recombined_Swarbrick_Stroma_Neoplastic_annotated, reduction = "umap", group.by = "orig.ident", raster = FALSE)

########
# SAVE #
########
saveRDS(Recombined_Swarbrick_Stroma_Neoplastic_annotated, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Annotated/Recombined_Swarbrick_Stroma_Neoplastic_annotated.rds")

# Remove Object
rm(Recombined_Swarbrick_Stroma_Neoplastic_annotated)
gc()

##############################################################################
###### Extract Stromal Cells from Total Singlets #############################
##############################################################################
# Read RDS
all_PANNTHR_breast_tumors_Total_Singlets <- readRDS(file = "/R/R_PANNTHR/PANNTHR_RDS/RDS_Annotated/all_PANNTHR_breast_tumors_Total_Singlets_Annotated.rds")

#################################
# Extract Stromal Cells #########
#################################

# Extract cells not analyzed as neoplastic; this requires two rounds of subsetting
# First, remove epithelial/neoplastic clusters
all_PANNTHR_breast_tumors_stroma_annotated <- subset(x = all_PANNTHR_breast_tumors_Total_Singlets, idents = c("Epithelial/Neoplastic Cells"), invert =  TRUE)
# Next, remove KRT-pos cells that were gated and considered epithelial
all_PANNTHR_breast_tumors_stroma_annotated <- subset(x = all_PANNTHR_breast_tumors_stroma_annotated, subset = KRT14 < 1 & KRT18 < 1 & KRT19 < 1)

# Remove Initial Object
rm(all_PANNTHR_breast_tumors_Total_Singlets)
gc()

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
all_PANNTHR_breast_tumors_stroma_annotated <- FindVariableFeatures(all_PANNTHR_breast_tumors_stroma_annotated, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(all_PANNTHR_breast_tumors_stroma_annotated)
all_PANNTHR_breast_tumors_stroma_annotated <- ScaleData(all_PANNTHR_breast_tumors_stroma_annotated, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
all_PANNTHR_breast_tumors_stroma_annotated <- RunPCA(all_PANNTHR_breast_tumors_stroma_annotated, features = VariableFeatures(object = all_PANNTHR_breast_tumors_stroma_annotated))

# Examine and visualize PCA results a few different ways
print(all_PANNTHR_breast_tumors_stroma_annotated[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(all_PANNTHR_breast_tumors_stroma_annotated, dims = 1:2, reduction = "pca")
DimPlot(all_PANNTHR_breast_tumors_stroma_annotated, reduction = "pca")

# Visualize PCs
ElbowPlot(all_PANNTHR_breast_tumors_stroma_annotated, ndims = 50)

all_PANNTHR_breast_tumors_stroma_annotated <- FindNeighbors(all_PANNTHR_breast_tumors_stroma_annotated, dims = 1:40)
all_PANNTHR_breast_tumors_stroma_annotated <- FindClusters(all_PANNTHR_breast_tumors_stroma_annotated, resolution = 0.08)

# Umap clustering
all_PANNTHR_breast_tumors_stroma_annotated <- RunUMAP(all_PANNTHR_breast_tumors_stroma_annotated, dims = 1:40)
DimPlot(all_PANNTHR_breast_tumors_stroma_annotated, reduction = "umap", raster = FALSE)
DimPlot(all_PANNTHR_breast_tumors_stroma_annotated, reduction = "umap", group.by = "orig.ident", raster = FALSE)

# Save RDS
saveRDS(all_PANNTHR_breast_tumors_stroma_annotated, file = "/R/R_PANNTHR/PANNTHR_RDS/RDS_Annotated/all_PANNTHR_breast_tumors_stroma_annotated.rds")



###################################
######### PANNTHR Dataset #########
###################################
# Load Annotated SeuratObjects
all_PANNTHR_breast_tumors_stroma_annotated <- readRDS(file = "/R/R_PANNTHR/PANNTHR_RDS/RDS_Annotated/all_PANNTHR_breast_tumors_stroma_annotated.rds")
all_PANNTHR_breast_tumors_panCK_cnv_annotated <- readRDS(file = "/R/R_PANNTHR/PANNTHR_RDS/RDS_Annotated/all_PANNTHR_breast_tumors_panCK_cnv_annotated.rds")


###############################################################################
# Initial annotations should be retained under cell.ident; update annotations #
###############################################################################

# Set Next Identity
all_PANNTHR_breast_tumors_stroma_annotated <- SetIdent(all_PANNTHR_breast_tumors_stroma_annotated, value = "cell.ident")

# Check Idents
Idents(object = all_PANNTHR_breast_tumors_stroma_annotated)
table(Idents(all_PANNTHR_breast_tumors_stroma_annotated))

# Annotate Clusters = main.ident
Stromal_IDs_main <- c("Leukocytes", "Leukocytes", "Leukocytes", "Leukocytes", "Mesenchymal Cells", "Mesenchymal Cells", "Leukocytes", "Leukocytes")

# Apply New Labels to Clusters
names(Stromal_IDs_main) <- levels(all_PANNTHR_breast_tumors_stroma_annotated)
all_PANNTHR_breast_tumors_stroma_annotated <- RenameIdents(all_PANNTHR_breast_tumors_stroma_annotated, Stromal_IDs_main)

# Remove old Seurat Object and ID list
rm(Stromal_IDs_main)

# Store Annotations Under Classification
all_PANNTHR_breast_tumors_stroma_annotated[["main.ident"]] <- Idents(object = all_PANNTHR_breast_tumors_stroma_annotated)


##############################################
# Merge Individual Files into Single Dataset #
##############################################
# Merge Datasets
Recombined_PANNTHR_Stroma_Neoplastic_annotated <- merge(all_PANNTHR_breast_tumors_stroma_annotated, y = c(all_PANNTHR_breast_tumors_panCK_cnv_annotated))

# Remove starting objects
rm(all_PANNTHR_breast_tumors_stroma_annotated)
rm(all_PANNTHR_breast_tumors_panCK_cnv_annotated)
gc()

# Join Layers
Recombined_PANNTHR_Stroma_Neoplastic_annotated <- JoinLayers(Recombined_PANNTHR_Stroma_Neoplastic_annotated)


####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Recombined_PANNTHR_Stroma_Neoplastic_annotated <- FindVariableFeatures(Recombined_PANNTHR_Stroma_Neoplastic_annotated, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Recombined_PANNTHR_Stroma_Neoplastic_annotated)
Recombined_PANNTHR_Stroma_Neoplastic_annotated <- ScaleData(Recombined_PANNTHR_Stroma_Neoplastic_annotated, features = all.genes)


#####################
# REDUCE DIMENSIONS #
#####################
Recombined_PANNTHR_Stroma_Neoplastic_annotated <- RunPCA(Recombined_PANNTHR_Stroma_Neoplastic_annotated, features = VariableFeatures(object = Recombined_PANNTHR_Stroma_Neoplastic_annotated))

# Examine and visualize PCA results a few different ways
print(Recombined_PANNTHR_Stroma_Neoplastic_annotated[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Recombined_PANNTHR_Stroma_Neoplastic_annotated, dims = 1:2, reduction = "pca")
DimPlot(Recombined_PANNTHR_Stroma_Neoplastic_annotated, reduction = "pca")

# Visualize PCs
ElbowPlot(Recombined_PANNTHR_Stroma_Neoplastic_annotated, ndims = 50)

Recombined_PANNTHR_Stroma_Neoplastic_annotated <- FindNeighbors(Recombined_PANNTHR_Stroma_Neoplastic_annotated, dims = 1:50)
Recombined_PANNTHR_Stroma_Neoplastic_annotated <- FindClusters(Recombined_PANNTHR_Stroma_Neoplastic_annotated, resolution = 0.15)

# Umap clustering
Recombined_PANNTHR_Stroma_Neoplastic_annotated <- RunUMAP(Recombined_PANNTHR_Stroma_Neoplastic_annotated, dims = 1:50)
DimPlot(Recombined_PANNTHR_Stroma_Neoplastic_annotated, reduction = "umap", raster = FALSE, group.by = "detailed.ident")
DimPlot(Recombined_PANNTHR_Stroma_Neoplastic_annotated, reduction = "umap", group.by = "orig.ident", raster = FALSE)

########
# SAVE #
########
saveRDS(Recombined_PANNTHR_Stroma_Neoplastic_annotated, file = "/R/R_PANNTHR/PANNTHR_RDS/RDS_Annotated/Recombined_PANNTHR_Stroma_Neoplastic_annotated.rds")

# Remove Object
rm(Recombined_PANNTHR_Stroma_Neoplastic_annotated)
gc()



##############################################################
# Step 21: Merged All Annotated Cells from All Breast Tumors #  
##############################################################
# Load Libraries
library(Seurat)
library(dplyr)


# Load SeuratObjects with all cells and corresponding annotations
Recombined_Visvader_Stroma_Neoplastic_annotated <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Annotated/Recombined_Visvader_Stroma_Neoplastic_annotated.rds")
Recombined_Swarbrick_Stroma_Neoplastic_annotated <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Annotated/Recombined_Swarbrick_Stroma_Neoplastic_annotated.rds")
Recombined_PANNTHR_Stroma_Neoplastic_annotated <- readRDS(file = "/R/R_PANNTHR/PANNTHR_RDS/RDS_Annotated/Recombined_PANNTHR_Stroma_Neoplastic_annotated.rds")


##############################################
# Merge Individual Files into Single Dataset #
##############################################
# Merge Datasets
Recombined_All_Breast_Tumors_Stroma_Neoplastic_annotated <- merge(Recombined_Visvader_Stroma_Neoplastic_annotated, y = c(Recombined_Swarbrick_Stroma_Neoplastic_annotated, Recombined_PANNTHR_Stroma_Neoplastic_annotated))


# The PANNTHR dataset contains Ensembl IDs and LCRNAs, and so only genes common to all datasets are kept
# Identify genes that are present in all three objects
Shared_genes <- Reduce(intersect, list(rownames(Recombined_Visvader_Stroma_Neoplastic_annotated), rownames(Recombined_Swarbrick_Stroma_Neoplastic_annotated), rownames(Recombined_PANNTHR_Stroma_Neoplastic_annotated)))

# Subset the merged object to keep only common genes
Recombined_All_Breast_Tumors_Stroma_Neoplastic_annotated <- Recombined_All_Breast_Tumors_Stroma_Neoplastic_annotated[Shared_genes, ]

# Check the resulting object
print(Recombined_All_Breast_Tumors_Stroma_Neoplastic_annotated)


# Remove starting objects
rm(Recombined_Visvader_Stroma_Neoplastic_annotated)
rm(Recombined_Swarbrick_Stroma_Neoplastic_annotated)
rm(Recombined_PANNTHR_Stroma_Neoplastic_annotated)
rm(Shared_genes)
gc()

# Join Layers
Recombined_All_Breast_Tumors_Stroma_Neoplastic_annotated <- JoinLayers(Recombined_All_Breast_Tumors_Stroma_Neoplastic_annotated)


####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Recombined_All_Breast_Tumors_Stroma_Neoplastic_annotated <- FindVariableFeatures(Recombined_All_Breast_Tumors_Stroma_Neoplastic_annotated, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Recombined_All_Breast_Tumors_Stroma_Neoplastic_annotated)
Recombined_All_Breast_Tumors_Stroma_Neoplastic_annotated <- ScaleData(Recombined_All_Breast_Tumors_Stroma_Neoplastic_annotated, features = all.genes)


#####################
# REDUCE DIMENSIONS #
#####################
Recombined_All_Breast_Tumors_Stroma_Neoplastic_annotated <- RunPCA(Recombined_All_Breast_Tumors_Stroma_Neoplastic_annotated, features = VariableFeatures(object = Recombined_All_Breast_Tumors_Stroma_Neoplastic_annotated))

# Examine and visualize PCA results a few different ways
print(Recombined_All_Breast_Tumors_Stroma_Neoplastic_annotated[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Recombined_All_Breast_Tumors_Stroma_Neoplastic_annotated, dims = 1:2, reduction = "pca")
DimPlot(Recombined_All_Breast_Tumors_Stroma_Neoplastic_annotated, reduction = "pca")

# Visualize PCs
ElbowPlot(Recombined_All_Breast_Tumors_Stroma_Neoplastic_annotated, ndims = 50)

Recombined_All_Breast_Tumors_Stroma_Neoplastic_annotated <- FindNeighbors(Recombined_All_Breast_Tumors_Stroma_Neoplastic_annotated, dims = 1:50)
Recombined_All_Breast_Tumors_Stroma_Neoplastic_annotated <- FindClusters(Recombined_All_Breast_Tumors_Stroma_Neoplastic_annotated, resolution = 0.12)

# Umap clustering
Recombined_All_Breast_Tumors_Stroma_Neoplastic_annotated <- RunUMAP(Recombined_All_Breast_Tumors_Stroma_Neoplastic_annotated, dims = 1:50)
DimPlot(Recombined_All_Breast_Tumors_Stroma_Neoplastic_annotated, reduction = "umap", raster = FALSE, group.by = "detailed.ident")
DimPlot(Recombined_All_Breast_Tumors_Stroma_Neoplastic_annotated, reduction = "umap", group.by = "orig.ident", raster = FALSE)

########
# SAVE #
########
saveRDS(Recombined_All_Breast_Tumors_Stroma_Neoplastic_annotated, file = "/R/R_All_Breast_Tumors/All_Breast_Tumors_RDS/Recombined_All_Breast_Tumors_Stroma_Neoplastic_annotated.rds")
gc()



#####################################################################################################################
# Step 22: Identify Differentially Expressed Genes in Immune-mimicked vs Neoplastic or Immune-mimicked vs Leukocyte #
#####################################################################################################################
# Load Combined RDS with IM and Stromal Annotations
#Recombined_All_Breast_Tumors_Stroma_Neoplastic_annotated <- readRDS(file = "/R/R_All_Breast_Tumors/All_Breast_Tumors_RDS/Recombined_All_Breast_Tumors_Stroma_Neoplastic_annotated.rds")


# Set Identity
Recombined_All_Breast_Tumors_Stroma_Neoplastic_annotated <- SetIdent(Recombined_All_Breast_Tumors_Stroma_Neoplastic_annotated, value = "main.ident")

# Verify Idents
Idents(object = Recombined_All_Breast_Tumors_Stroma_Neoplastic_annotated)
table(Idents(Recombined_All_Breast_Tumors_Stroma_Neoplastic_annotated))


# Perform differential gene expression analysis for immune mimicry vs Neoplastic Cells
Immune_Mimicry_vs_Neoplastic.de.markers <- FindMarkers(Recombined_All_Breast_Tumors_Stroma_Neoplastic_annotated, ident.1 = "Immune-like", ident.2 = "Neoplastic Cells", only.pos = FALSE)
write.csv(Immune_Mimicry_vs_Neoplastic.de.markers, file = "/R/R_All_Breast_Tumors/All_Breast_Tumors_Output/Immune_Mimicry_vs_Neoplastic.de.markers.csv")
# Remove Markers
rm(Immune_Mimicry_vs_Neoplastic.de.markers)


# Perform differential gene expression analysis for immune mimicry vs Leukocytes
Immune_Mimicry_vs_Leukocytes.de.markers <- FindMarkers(Recombined_All_Breast_Tumors_Stroma_Neoplastic_annotated, ident.1 = "Immune-like", ident.2 = "Leukocytes", only.pos = FALSE)
write.csv(Immune_Mimicry_vs_Leukocytes.de.markers, file = "/R/R_All_Breast_Tumors/All_Breast_Tumors_Output/Immune_Mimicry_vs_Leukocytes.de.markers.csv")
# Remove Markers
rm(Immune_Mimicry_vs_Leukocytes.de.markers)


# Perform differential gene expression analysis for Leukocytes vs Neoplastic Cells
Leukocytes_vs_Neoplastic.de.markers <- FindMarkers(Recombined_All_Breast_Tumors_Stroma_Neoplastic_annotated, ident.1 = "Leukocytes", ident.2 = "Neoplastic Cells", only.pos = FALSE)
write.csv(Leukocytes_vs_Neoplastic.de.markers, file = "/R/R_All_Breast_Tumors/All_Breast_Tumors_Output/Leukocytes_vs_Neoplastic.de.markers.csv")
# Remove Markers
rm(Leukocytes_vs_Neoplastic.de.markers)


# Perform differential gene expression analysis for immune mimicry vs Leukocytes
Immune_Mimicry_vs_Leukocytes.de.markers <- FindMarkers(Recombined_All_Breast_Tumors_Stroma_Neoplastic_annotated, ident.1 = "Immune-like", ident.2 = "Leukocytes", only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(Immune_Mimicry_vs_Leukocytes.de.markers, file = "/R/R_All_Breast_Tumors/All_Breast_Tumors_Output/Immune_Mimicry_vs_Leukocytes.de.markers_0.25_filter.csv")
# Remove Markers
rm(Immune_Mimicry_vs_Leukocytes.de.markers)



######################################################################
# Step 23: Determine Biological Processes Related to Cell Type DEGs  #
######################################################################
# Load Libraries
library(tidyverse)
library(data.table)
library(fgsea)


######################################
# Immune_Mimicry_vs_Neoplastic Cells #
######################################

# Load Differential Expression Results
Immune_Mimicry_vs_Neoplastic.de.markers <- read.csv(file = "/R/R_All_Breast_Tumors/All_Breast_Tumors_Output/Immune_Mimicry_vs_Neoplastic.de.markers.csv")

# Filter DEGs to only significant genes with p < 0.01
Immune_Mimicry_vs_Neoplastic.de.markers <- Immune_Mimicry_vs_Neoplastic.de.markers %>% filter(p_val_adj<0.01)

# Make first column gene name
# Note column 1 is the gene name
colnames(Immune_Mimicry_vs_Neoplastic.de.markers)[1] <- "Gene"

# Collect columns: gene names, avg_log2FC, and cluster; sort by descending 
Immune_Mimicry_vs_Neoplastic.de.markers <- Immune_Mimicry_vs_Neoplastic.de.markers %>% dplyr::select(Gene, avg_log2FC)
Immune_Mimicry_vs_Neoplastic.de.markers <- Immune_Mimicry_vs_Neoplastic.de.markers[order(Immune_Mimicry_vs_Neoplastic.de.markers$avg_log2FC, decreasing = TRUE),]  

# Save Output
write.csv(Immune_Mimicry_vs_Neoplastic.de.markers, file = "/R/R_All_Breast_Tumors/All_Breast_Tumors_fgsea/All_Breast_Tumors_fgsea_Input/Immune_Mimicry_vs_Neoplastic.de.markers_filtered.csv")


# Collapse into gene name and LogFC
Immune_Mimicry_vs_Neoplastic.de.markers <- Immune_Mimicry_vs_Neoplastic.de.markers %>% dplyr::select(Gene, avg_log2FC)


# Create a vector of gene ranks
Immune_Mimicry_vs_Neoplastic.de.markers.ranks <- deframe(Immune_Mimicry_vs_Neoplastic.de.markers)


# Load the pathways into a named list
pathways.GO.BP <- gmtPathways("/R/R_All_Breast_Tumors/All_Breast_Tumors_fgsea/All_Breast_Tumors_fgsea_MSigDB/c5.go.bp.v2022.1.Hs.symbols.gmt")


########################################################
# Run fgsea on Immune_Mimicry_vs_Neoplastic.de.markers #
########################################################
fgsea.Immune_Mimicry_vs_Neoplastic.de.markers.res <- fgsea(pathways = pathways.GO.BP, 
                                                           stats    = Immune_Mimicry_vs_Neoplastic.de.markers.ranks,
                                                           eps      = 0.0,
                                                           minSize  = 15,
                                                           maxSize  = 500, nPermSimple = 10000)

# Summarize Immune_Mimicry_vs_Neoplastic.de.markers.results by collapsing into main pathways, to omit redundancies
Immune_Mimicry_vs_Neoplastic.de.markers.collapsedPathways <- collapsePathways(fgsea.Immune_Mimicry_vs_Neoplastic.de.markers.res[order(pval)][padj < 0.05], 
                                                                              pathways.GO.BP, Immune_Mimicry_vs_Neoplastic.de.markers.ranks)
Immune_Mimicry_vs_Neoplastic.de.markers.mainPathways <- fgsea.Immune_Mimicry_vs_Neoplastic.de.markers.res[pathway %in% Immune_Mimicry_vs_Neoplastic.de.markers.collapsedPathways$mainPathways][
  order(-NES), pathway]
plotGseaTable(pathways.GO.BP[Immune_Mimicry_vs_Neoplastic.de.markers.mainPathways], Immune_Mimicry_vs_Neoplastic.de.markers.ranks, fgsea.Immune_Mimicry_vs_Neoplastic.de.markers.res, 
              gseaParam = 0.5)

# Filter Immune_Mimicry_vs_Neoplastic.de.markers.results so they only include the main pathways
fgsea.Immune_Mimicry_vs_Neoplastic.de.markers.res_collapsed <- filter(fgsea.Immune_Mimicry_vs_Neoplastic.de.markers.res, pathway %in% Immune_Mimicry_vs_Neoplastic.de.markers.mainPathways)

# Save Immune_Mimicry_vs_Neoplastic.de.markers.results
fwrite(fgsea.Immune_Mimicry_vs_Neoplastic.de.markers.res, file="/R/R_All_Breast_Tumors/All_Breast_Tumors_fgsea/All_Breast_Tumors_fgsea_Output/fgsea.Immune_Mimicry_vs_Neoplastic.de.markers.res.csv")
fwrite(fgsea.Immune_Mimicry_vs_Neoplastic.de.markers.res_collapsed, file="/R/R_All_Breast_Tumors/All_Breast_Tumors_fgsea/All_Breast_Tumors_fgsea_Output/fgsea.Immune_Mimicry_vs_Neoplastic.de.markers.res_collapsed.csv")

# Retain pathways with adjusted p-values below <0.05 
fgsea.Immune_Mimicry_vs_Neoplastic.de.markers.res_sig <- fgsea.Immune_Mimicry_vs_Neoplastic.de.markers.res %>% filter(padj<0.05)

# Save Immune_Mimicry_vs_Neoplastic.de.markers.results
fwrite(fgsea.Immune_Mimicry_vs_Neoplastic.de.markers.res_sig, file="/R/R_All_Breast_Tumors/All_Breast_Tumors_fgsea/All_Breast_Tumors_fgsea_Output/fgsea.Immune_Mimicry_vs_Neoplastic.de.markers.res_sig.csv")


# Remove objects
rm(fgsea.Immune_Mimicry_vs_Neoplastic.de.markers.res)
rm(fgsea.Immune_Mimicry_vs_Neoplastic.de.markers.res_collapsed)
rm(Immune_Mimicry_vs_Neoplastic.de.markers.collapsedPathways)
rm(Immune_Mimicry_vs_Neoplastic.de.markers.mainPathways)
rm(Immune_Mimicry_vs_Neoplastic.de.markers.ranks)
rm(fgsea.Immune_Mimicry_vs_Neoplastic.de.markers.res_sig)



################################
# Immune_Mimicry_vs_Leukocytes #
################################

# Load Differential Expression Results
Immune_Mimicry_vs_Leukocytes.de.markers <- read.csv(file = "/R/R_All_Breast_Tumors/All_Breast_Tumors_Output/Immune_Mimicry_vs_Leukocytes.de.markers.csv")

# Filter DEGs to only significant genes with p < 0.01
Immune_Mimicry_vs_Leukocytes.de.markers <- Immune_Mimicry_vs_Leukocytes.de.markers %>% filter(p_val_adj<0.01)

# Make first column gene name
# Note column 1 is the gene name
colnames(Immune_Mimicry_vs_Leukocytes.de.markers)[1] <- "Gene"

# Collect columns: gene names, avg_log2FC, and cluster; sort by descending 
Immune_Mimicry_vs_Leukocytes.de.markers <- Immune_Mimicry_vs_Leukocytes.de.markers %>% dplyr::select(Gene, avg_log2FC)
Immune_Mimicry_vs_Leukocytes.de.markers <- Immune_Mimicry_vs_Leukocytes.de.markers[order(Immune_Mimicry_vs_Leukocytes.de.markers$avg_log2FC, decreasing = TRUE),]  

# Save Output
write.csv(Immune_Mimicry_vs_Leukocytes.de.markers, file = "/R/R_All_Breast_Tumors/All_Breast_Tumors_fgsea/All_Breast_Tumors_fgsea_Input/Immune_Mimicry_vs_Leukocytes.de.markers_filtered.csv")


# Collapse into gene name and LogFC
Immune_Mimicry_vs_Leukocytes.de.markers <- Immune_Mimicry_vs_Leukocytes.de.markers %>% dplyr::select(Gene, avg_log2FC)


# Create a vector of gene ranks
Immune_Mimicry_vs_Leukocytes.de.markers.ranks <- deframe(Immune_Mimicry_vs_Leukocytes.de.markers)


# Load the pathways into a named list
pathways.GO.BP <- gmtPathways("/R/R_All_Breast_Tumors/All_Breast_Tumors_fgsea/All_Breast_Tumors_fgsea_MSigDB/c5.go.bp.v2022.1.Hs.symbols.gmt")


########################################################
# Run fgsea on Immune_Mimicry_vs_Leukocytes.de.markers #
########################################################
fgsea.Immune_Mimicry_vs_Leukocytes.de.markers.res <- fgsea(pathways = pathways.GO.BP, 
                                                           stats    = Immune_Mimicry_vs_Leukocytes.de.markers.ranks,
                                                           eps      = 0.0,
                                                           minSize  = 15,
                                                           maxSize  = 500, nPermSimple = 10000)

# Summarize Immune_Mimicry_vs_Leukocytes.de.markers.results by collapsing into main pathways, to omit redundancies
Immune_Mimicry_vs_Leukocytes.de.markers.collapsedPathways <- collapsePathways(fgsea.Immune_Mimicry_vs_Leukocytes.de.markers.res[order(pval)][padj < 0.05], 
                                                                              pathways.GO.BP, Immune_Mimicry_vs_Leukocytes.de.markers.ranks)
Immune_Mimicry_vs_Leukocytes.de.markers.mainPathways <- fgsea.Immune_Mimicry_vs_Leukocytes.de.markers.res[pathway %in% Immune_Mimicry_vs_Leukocytes.de.markers.collapsedPathways$mainPathways][
  order(-NES), pathway]
plotGseaTable(pathways.GO.BP[Immune_Mimicry_vs_Leukocytes.de.markers.mainPathways], Immune_Mimicry_vs_Leukocytes.de.markers.ranks, fgsea.Immune_Mimicry_vs_Leukocytes.de.markers.res, 
              gseaParam = 0.5)

# Filter Immune_Mimicry_vs_Leukocytes.de.markers.results so they only include the main pathways
fgsea.Immune_Mimicry_vs_Leukocytes.de.markers.res_collapsed <- filter(fgsea.Immune_Mimicry_vs_Leukocytes.de.markers.res, pathway %in% Immune_Mimicry_vs_Leukocytes.de.markers.mainPathways)

# Save Immune_Mimicry_vs_Leukocytes.de.markers.results
fwrite(fgsea.Immune_Mimicry_vs_Leukocytes.de.markers.res, file="/R/R_All_Breast_Tumors/All_Breast_Tumors_fgsea/All_Breast_Tumors_fgsea_Output/fgsea.Immune_Mimicry_vs_Leukocytes.de.markers.res.csv")
fwrite(fgsea.Immune_Mimicry_vs_Leukocytes.de.markers.res_collapsed, file="/R/R_All_Breast_Tumors/All_Breast_Tumors_fgsea/All_Breast_Tumors_fgsea_Output/fgsea.Immune_Mimicry_vs_Leukocytes.de.markers.res_collapsed.csv")

# Retain pathways with adjusted p-values below <0.05 
fgsea.Immune_Mimicry_vs_Leukocytes.de.markers.res_sig <- fgsea.Immune_Mimicry_vs_Leukocytes.de.markers.res %>% filter(padj<0.05)

# Save Immune_Mimicry_vs_Leukocytes.de.markers.results
fwrite(fgsea.Immune_Mimicry_vs_Leukocytes.de.markers.res_sig, file="/R/R_All_Breast_Tumors/All_Breast_Tumors_fgsea/All_Breast_Tumors_fgsea_Output/fgsea.Immune_Mimicry_vs_Leukocytes.de.markers.res_sig.csv")


# Remove objects
rm(fgsea.Immune_Mimicry_vs_Leukocytes.de.markers.res)
rm(fgsea.Immune_Mimicry_vs_Leukocytes.de.markers.res_collapsed)
rm(Immune_Mimicry_vs_Leukocytes.de.markers.collapsedPathways)
rm(Immune_Mimicry_vs_Leukocytes.de.markers.mainPathways)
rm(Immune_Mimicry_vs_Leukocytes.de.markers.ranks)
rm(fgsea.Immune_Mimicry_vs_Leukocytes.de.markers.res_sig)



###############################################################################################
# Step 24: Quantify Breast Tumor Cells Expressing IM Markers Inside vs Outside IM Clusters    #
###############################################################################################
# Load Combined RDS with IM and Stromal Annotations
Recombined_All_Breast_Tumors_Stroma_Neoplastic_annotated <- readRDS(file = "/R/R_All_Breast_Tumors/All_Breast_Tumors_RDS/Recombined_All_Breast_Tumors_Stroma_Neoplastic_annotated.rds")


# Set Identity
Recombined_All_Breast_Tumors_Stroma_Neoplastic_annotated <- SetIdent(Recombined_All_Breast_Tumors_Stroma_Neoplastic_annotated, value = "main.ident")

# Verify Idents
Idents(object = Recombined_All_Breast_Tumors_Stroma_Neoplastic_annotated)
table(Idents(Recombined_All_Breast_Tumors_Stroma_Neoplastic_annotated))

# Collect Immune-like Cells and Neoplastic Cells
Immune_like <- subset(x = Recombined_All_Breast_Tumors_Stroma_Neoplastic_annotated, idents = c("Immune-like"))
Neoplastic_cells <- subset(x = Recombined_All_Breast_Tumors_Stroma_Neoplastic_annotated, idents = c("Neoplastic Cells"))

# Remove Starting Object
#rm(Recombined_All_Breast_Tumors_Stroma_Neoplastic_annotated)
gc()


##############################################
# Enumerate IMBC Genes for Immune_like Cells #
##############################################
# Extract the gene expression matrix
expression_matrix <- GetAssayData(Immune_like, slot = "data")
# Extract the metadata to identify samples
metadata <- Immune_like@meta.data


# Initialize a list to store counts for each gene
gene_counts_list <- list()

# Define a function to process a gene of interest
process_gene <- function(gene_symbol, expression_matrix, metadata) {
  if (gene_symbol %in% rownames(expression_matrix)) {
    # Identify cells expressing the gene (expression > 0)
    gene_positive_cells <- expression_matrix[gene_symbol, ] > 0
    
    # Count the total number of positive cells across all samples
    total_positive_cells <- sum(gene_positive_cells)
    
    # Store the count in the list
    gene_counts_list[[gene_symbol]] <<- total_positive_cells
    
    # Print the count
    print(paste(gene_symbol, ":", total_positive_cells, "cells express the gene."))
  } else {
    print(paste(gene_symbol, "gene is not present in the dataset."))
  }
}

# Example usage for multiple genes
genes_of_interest <- c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44")

for (gene in genes_of_interest) {
  process_gene(gene, expression_matrix, metadata)
}

# Combine all gene counts into a single data frame
combined_counts <- data.frame(gene_symbol = names(gene_counts_list), count = unlist(gene_counts_list))

# Sort the combined counts alphabetically by gene symbol
combined_counts <- combined_counts[order(combined_counts$gene_symbol), ]
print(combined_counts)

# Save the combined counts to a single CSV file
write.csv(combined_counts, file = "/R/R_All_Breast_Tumors/All_Breast_Tumors_Output/All_Breast_Tumors_Immune_like_IMBC_Genes.csv", row.names = FALSE)

# Remove Objects
rm(gene)
rm(gene_counts_list)
rm(genes_of_interest)
rm(expression_matrix)
rm(metadata)
rm(process_gene)
rm(combined_counts)
rm(Immune_like)


#######################################################################
# Enumerate IMBC Genes for Neoplastic_cells                           #
#######################################################################
# Extract the gene expression matrix
expression_matrix <- GetAssayData(Neoplastic_cells, slot = "data")
# Extract the metadata to identify samples
metadata <- Neoplastic_cells@meta.data


# Initialize a list to store counts for each gene
gene_counts_list <- list()

# Define a function to process a gene of interest
process_gene <- function(gene_symbol, expression_matrix, metadata) {
  if (gene_symbol %in% rownames(expression_matrix)) {
    # Identify cells expressing the gene (expression > 0)
    gene_positive_cells <- expression_matrix[gene_symbol, ] > 0
    
    # Count the total number of positive cells across all samples
    total_positive_cells <- sum(gene_positive_cells)
    
    # Store the count in the list
    gene_counts_list[[gene_symbol]] <<- total_positive_cells
    
    # Print the count
    print(paste(gene_symbol, ":", total_positive_cells, "cells express the gene."))
  } else {
    print(paste(gene_symbol, "gene is not present in the dataset."))
  }
}

# Example usage for multiple genes
genes_of_interest <- c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44")

for (gene in genes_of_interest) {
  process_gene(gene, expression_matrix, metadata)
}

# Combine all gene counts into a single data frame
combined_counts <- data.frame(gene_symbol = names(gene_counts_list), count = unlist(gene_counts_list))

# Sort the combined counts alphabetically by gene symbol
combined_counts <- combined_counts[order(combined_counts$gene_symbol), ]
print(combined_counts)

# Save the combined counts to a single CSV file
write.csv(combined_counts, file = "/R/R_All_Breast_Tumors/All_Breast_Tumors_Output/All_Breast_Tumors_Neoplastic_cells_IMBC_Genes.csv", row.names = FALSE)

# Remove Objects
rm(gene)
rm(gene_counts_list)
rm(genes_of_interest)
rm(expression_matrix)
rm(metadata)
rm(process_gene)
rm(combined_counts)
rm(Neoplastic_cells)



#######################################################################################################################
# Step 25: Evaluate Breast Cancer Stem Gene Expression per Immune Mimicry Marker & Generate Heatmap                   #
#######################################################################################################################
# Load Combined RDS with IM and Stromal Annotations
#Recombined_All_Breast_Tumors_Stroma_Neoplastic_annotated <- readRDS(file = "/R/R_All_Breast_Tumors/All_Breast_Tumors_RDS/Recombined_All_Breast_Tumors_Stroma_Neoplastic_annotated.rds")


# Set Identity
Recombined_All_Breast_Tumors_Stroma_Neoplastic_annotated <- SetIdent(Recombined_All_Breast_Tumors_Stroma_Neoplastic_annotated, value = "main.ident")

# Verify Idents
Idents(object = Recombined_All_Breast_Tumors_Stroma_Neoplastic_annotated)
table(Idents(Recombined_All_Breast_Tumors_Stroma_Neoplastic_annotated))

DimPlot(Recombined_All_Breast_Tumors_Stroma_Neoplastic_annotated)

# Extract Cell Type Subsets
Immune_like_and_Neoplastic <- subset(x = Recombined_All_Breast_Tumors_Stroma_Neoplastic_annotated, idents = c("Immune-like", "Neoplastic Cells"))
# Remove all cells expressing immune surface receptors
Neoplastic_cells_without_IM_Markers <- subset(x = Neoplastic_cells, subset = FCGR3A == 0 | CD2 == 0 | CD3E == 0 | CD7 == 0 | CD3G == 0 | CD3D == 0 | IL7R == 0 | FCGR2A == 0 | CD68 == 0 | CD52 == 0 | CD69 == 0 | CD14 == 0 | PTPRC == 0 | TNFSF13B == 0 | CD37 == 0 | ITGB2 == 0 | CXCR4 == 0 | CD74 == 0 | CLEC7A == 0 | CD53 == 0 | CD83 == 0 | PLAUR == 0 | CD44 == 0)

# Remove Starting Object
#rm(Recombined_All_Breast_Tumors_Stroma_Neoplastic_annotated)


########################
# Subset FCGR3A+ Cells #
########################
# Subsetting all cells for FCGR3A
FCGR3A.pos.subset <- subset(x = Immune_like_and_Neoplastic, subset = FCGR3A > 0)
# Export FCGR3A Expression Matrix for Stem-like Genes
FCGR3A.pos.subset <- FCGR3A.pos.subset[["RNA"]]$data
FCGR3A.pos.subset <- as.matrix(FCGR3A.pos.subset, 'sparseMatrix')
FCGR3A.pos.subset <- subset(FCGR3A.pos.subset, rownames(FCGR3A.pos.subset) %in% c("NANOG", "SOX2", "CD44", "CXCR4", "CD70",  "ITGB3", "SNAI1", "PROCR", "ITGB1", "ZEB1", "CD34", "THY1", "GLI2", "NOTCH1", "GLI1", "BMI1", "EPCAM", "PROM1", "ITGA6", "POU5F1", "LGR5", "KIT", "ALDH1A1"))
FCGR3A.pos.subset <- as.data.frame(rowMeans(FCGR3A.pos.subset))
# Further organize and save data frame
FCGR3A.pos.subset <-t(FCGR3A.pos.subset)
write.csv(FCGR3A.pos.subset, file = "/R/R_All_Breast_Tumors/All_Breast_Tumors_Output/FCGR3A.All.pos.subset_stem_markers.csv")
rm(FCGR3A.pos.subset)
gc()

#####################
# Subset CD2+ Cells #
#####################
# Subsetting all cells for CD2
CD2.pos.subset <- subset(x = Immune_like_and_Neoplastic, subset = CD2 > 0)
# Export CD2 Expression Matrix for Stem-like Genes
CD2.pos.subset <- CD2.pos.subset[["RNA"]]$data
CD2.pos.subset <- as.matrix(CD2.pos.subset, 'sparseMatrix')
CD2.pos.subset <- subset(CD2.pos.subset, rownames(CD2.pos.subset) %in% c("NANOG", "SOX2", "CD44", "CXCR4", "CD70",  "ITGB3", "SNAI1", "PROCR", "ITGB1", "ZEB1", "CD34", "THY1", "GLI2", "NOTCH1", "GLI1", "BMI1", "EPCAM", "PROM1", "ITGA6", "POU5F1", "LGR5", "KIT", "ALDH1A1"))
CD2.pos.subset <- as.data.frame(rowMeans(CD2.pos.subset))
# Further organize and save data frame
CD2.pos.subset <-t(CD2.pos.subset)
write.csv(CD2.pos.subset, file = "/R/R_All_Breast_Tumors/All_Breast_Tumors_Output/CD2.All.pos.subset_stem_markers.csv")
rm(CD2.pos.subset)
gc()

######################
# Subset CD3E+ Cells #
######################
# Subsetting all cells for CD3E
CD3E.pos.subset <- subset(x = Immune_like_and_Neoplastic, subset = CD3E > 0)
# Export CD3E Expression Matrix for Stem-like Genes
CD3E.pos.subset <- CD3E.pos.subset[["RNA"]]$data
CD3E.pos.subset <- as.matrix(CD3E.pos.subset, 'sparseMatrix')
CD3E.pos.subset <- subset(CD3E.pos.subset, rownames(CD3E.pos.subset) %in% c("NANOG", "SOX2", "CD44", "CXCR4", "CD70",  "ITGB3", "SNAI1", "PROCR", "ITGB1", "ZEB1", "CD34", "THY1", "GLI2", "NOTCH1", "GLI1", "BMI1", "EPCAM", "PROM1", "ITGA6", "POU5F1", "LGR5", "KIT", "ALDH1A1"))
CD3E.pos.subset <- as.data.frame(rowMeans(CD3E.pos.subset))
# Further organize and save data frame
CD3E.pos.subset <-t(CD3E.pos.subset)
write.csv(CD3E.pos.subset, file = "/R/R_All_Breast_Tumors/All_Breast_Tumors_Output/CD3E.All.pos.subset_stem_markers.csv")
rm(CD3E.pos.subset)
gc()

######################
# Subset CD7+ Cells #
######################
# Subsetting all cells for CD7
CD7.pos.subset <- subset(x = Immune_like_and_Neoplastic, subset = CD7 > 0)
# Export CD7 Expression Matrix for Stem-like Genes
CD7.pos.subset <- CD7.pos.subset[["RNA"]]$data
CD7.pos.subset <- as.matrix(CD7.pos.subset, 'sparseMatrix')
CD7.pos.subset <- subset(CD7.pos.subset, rownames(CD7.pos.subset) %in% c("NANOG", "SOX2", "CD44", "CXCR4", "CD70",  "ITGB3", "SNAI1", "PROCR", "ITGB1", "ZEB1", "CD34", "THY1", "GLI2", "NOTCH1", "GLI1", "BMI1", "EPCAM", "PROM1", "ITGA6", "POU5F1", "LGR5", "KIT", "ALDH1A1"))
CD7.pos.subset <- as.data.frame(rowMeans(CD7.pos.subset))
# Further organize and save data frame
CD7.pos.subset <-t(CD7.pos.subset)
write.csv(CD7.pos.subset, file = "/R/R_All_Breast_Tumors/All_Breast_Tumors_Output/CD7.All.pos.subset_stem_markers.csv")
rm(CD7.pos.subset)
gc()

######################
# Subset CD3G+ Cells #
######################
# Subsetting all cells for CD3G
CD3G.pos.subset <- subset(x = Immune_like_and_Neoplastic, subset = CD3G > 0)
# Export CD3G Expression Matrix for Stem-like Genes
CD3G.pos.subset <- CD3G.pos.subset[["RNA"]]$data
CD3G.pos.subset <- as.matrix(CD3G.pos.subset, 'sparseMatrix')
CD3G.pos.subset <- subset(CD3G.pos.subset, rownames(CD3G.pos.subset) %in% c("NANOG", "SOX2", "CD44", "CXCR4", "CD70",  "ITGB3", "SNAI1", "PROCR", "ITGB1", "ZEB1", "CD34", "THY1", "GLI2", "NOTCH1", "GLI1", "BMI1", "EPCAM", "PROM1", "ITGA6", "POU5F1", "LGR5", "KIT", "ALDH1A1"))
CD3G.pos.subset <- as.data.frame(rowMeans(CD3G.pos.subset))
# Further organize and save data frame
CD3G.pos.subset <-t(CD3G.pos.subset)
write.csv(CD3G.pos.subset, file = "/R/R_All_Breast_Tumors/All_Breast_Tumors_Output/CD3G.All.pos.subset_stem_markers.csv")
rm(CD3G.pos.subset)
gc()

######################
# Subset CD3D+ Cells #
######################
# Subsetting all cells for CD3D
CD3D.pos.subset <- subset(x = Immune_like_and_Neoplastic, subset = CD3D > 0)
# Export CD3D Expression Matrix for Stem-like Genes
CD3D.pos.subset <- CD3D.pos.subset[["RNA"]]$data
CD3D.pos.subset <- as.matrix(CD3D.pos.subset, 'sparseMatrix')
CD3D.pos.subset <- subset(CD3D.pos.subset, rownames(CD3D.pos.subset) %in% c("NANOG", "SOX2", "CD44", "CXCR4", "CD70",  "ITGB3", "SNAI1", "PROCR", "ITGB1", "ZEB1", "CD34", "THY1", "GLI2", "NOTCH1", "GLI1", "BMI1", "EPCAM", "PROM1", "ITGA6", "POU5F1", "LGR5", "KIT", "ALDH1A1"))
CD3D.pos.subset <- as.data.frame(rowMeans(CD3D.pos.subset))
# Further organize and save data frame
CD3D.pos.subset <-t(CD3D.pos.subset)
write.csv(CD3D.pos.subset, file = "/R/R_All_Breast_Tumors/All_Breast_Tumors_Output/CD3D.All.pos.subset_stem_markers.csv")
rm(CD3D.pos.subset)
gc()

######################
# Subset IL7R+ Cells #
######################
# Subsetting all cells for IL7R
IL7R.pos.subset <- subset(x = Immune_like_and_Neoplastic, subset = IL7R > 0)
# Export IL7R Expression Matrix for Stem-like Genes
IL7R.pos.subset <- IL7R.pos.subset[["RNA"]]$data
IL7R.pos.subset <- as.matrix(IL7R.pos.subset, 'sparseMatrix')
IL7R.pos.subset <- subset(IL7R.pos.subset, rownames(IL7R.pos.subset) %in% c("NANOG", "SOX2", "CD44", "CXCR4", "CD70",  "ITGB3", "SNAI1", "PROCR", "ITGB1", "ZEB1", "CD34", "THY1", "GLI2", "NOTCH1", "GLI1", "BMI1", "EPCAM", "PROM1", "ITGA6", "POU5F1", "LGR5", "KIT", "ALDH1A1"))
IL7R.pos.subset <- as.data.frame(rowMeans(IL7R.pos.subset))
# Further organize and save data frame
IL7R.pos.subset <-t(IL7R.pos.subset)
write.csv(IL7R.pos.subset, file = "/R/R_All_Breast_Tumors/All_Breast_Tumors_Output/IL7R.All.pos.subset_stem_markers.csv")
rm(IL7R.pos.subset)
gc()

########################
# Subset FCGR2A+ Cells #
########################
# Subsetting all cells for FCGR2A
FCGR2A.pos.subset <- subset(x = Immune_like_and_Neoplastic, subset = FCGR2A > 0)
# Export FCGR2A Expression Matrix for Stem-like Genes
FCGR2A.pos.subset <- FCGR2A.pos.subset[["RNA"]]$data
FCGR2A.pos.subset <- as.matrix(FCGR2A.pos.subset, 'sparseMatrix')
FCGR2A.pos.subset <- subset(FCGR2A.pos.subset, rownames(FCGR2A.pos.subset) %in% c("NANOG", "SOX2", "CD44", "CXCR4", "CD70",  "ITGB3", "SNAI1", "PROCR", "ITGB1", "ZEB1", "CD34", "THY1", "GLI2", "NOTCH1", "GLI1", "BMI1", "EPCAM", "PROM1", "ITGA6", "POU5F1", "LGR5", "KIT", "ALDH1A1"))
FCGR2A.pos.subset <- as.data.frame(rowMeans(FCGR2A.pos.subset))
# Further organize and save data frame
FCGR2A.pos.subset <-t(FCGR2A.pos.subset)
write.csv(FCGR2A.pos.subset, file = "/R/R_All_Breast_Tumors/All_Breast_Tumors_Output/FCGR2A.All.pos.subset_stem_markers.csv")
rm(FCGR2A.pos.subset)
gc()

######################
# Subset CD68+ Cells #
######################
# Subsetting all cells for CD68
CD68.pos.subset <- subset(x = Immune_like_and_Neoplastic, subset = CD68 > 0)
# Export CD68 Expression Matrix for Stem-like Genes
CD68.pos.subset <- CD68.pos.subset[["RNA"]]$data
CD68.pos.subset <- as.matrix(CD68.pos.subset, 'sparseMatrix')
CD68.pos.subset <- subset(CD68.pos.subset, rownames(CD68.pos.subset) %in% c("NANOG", "SOX2", "CD44", "CXCR4", "CD70",  "ITGB3", "SNAI1", "PROCR", "ITGB1", "ZEB1", "CD34", "THY1", "GLI2", "NOTCH1", "GLI1", "BMI1", "EPCAM", "PROM1", "ITGA6", "POU5F1", "LGR5", "KIT", "ALDH1A1"))
CD68.pos.subset <- as.data.frame(rowMeans(CD68.pos.subset))
# Further organize and save data frame
CD68.pos.subset <-t(CD68.pos.subset)
write.csv(CD68.pos.subset, file = "/R/R_All_Breast_Tumors/All_Breast_Tumors_Output/CD68.All.pos.subset_stem_markers.csv")
rm(CD68.pos.subset)
gc()

######################
# Subset CD52+ Cells #
######################
# Subsetting all cells for CD52
CD52.pos.subset <- subset(x = Immune_like_and_Neoplastic, subset = CD52 > 0)
# Export CD52 Expression Matrix for Stem-like Genes
CD52.pos.subset <- CD52.pos.subset[["RNA"]]$data
CD52.pos.subset <- as.matrix(CD52.pos.subset, 'sparseMatrix')
CD52.pos.subset <- subset(CD52.pos.subset, rownames(CD52.pos.subset) %in% c("NANOG", "SOX2", "CD44", "CXCR4", "CD70",  "ITGB3", "SNAI1", "PROCR", "ITGB1", "ZEB1", "CD34", "THY1", "GLI2", "NOTCH1", "GLI1", "BMI1", "EPCAM", "PROM1", "ITGA6", "POU5F1", "LGR5", "KIT", "ALDH1A1"))
CD52.pos.subset <- as.data.frame(rowMeans(CD52.pos.subset))
# Further organize and save data frame
CD52.pos.subset <-t(CD52.pos.subset)
write.csv(CD52.pos.subset, file = "/R/R_All_Breast_Tumors/All_Breast_Tumors_Output/CD52.All.pos.subset_stem_markers.csv")
rm(CD52.pos.subset)
gc()

######################
# Subset CD69+ Cells #
######################
# Subsetting all cells for CD69
CD69.pos.subset <- subset(x = Immune_like_and_Neoplastic, subset = CD69 > 0)
# Export CD69 Expression Matrix for Stem-like Genes
CD69.pos.subset <- CD69.pos.subset[["RNA"]]$data
CD69.pos.subset <- as.matrix(CD69.pos.subset, 'sparseMatrix')
CD69.pos.subset <- subset(CD69.pos.subset, rownames(CD69.pos.subset) %in% c("NANOG", "SOX2", "CD44", "CXCR4", "CD70",  "ITGB3", "SNAI1", "PROCR", "ITGB1", "ZEB1", "CD34", "THY1", "GLI2", "NOTCH1", "GLI1", "BMI1", "EPCAM", "PROM1", "ITGA6", "POU5F1", "LGR5", "KIT", "ALDH1A1"))
CD69.pos.subset <- as.data.frame(rowMeans(CD69.pos.subset))
# Further organize and save data frame
CD69.pos.subset <-t(CD69.pos.subset)
write.csv(CD69.pos.subset, file = "/R/R_All_Breast_Tumors/All_Breast_Tumors_Output/CD69.All.pos.subset_stem_markers.csv")
rm(CD69.pos.subset)
gc()

######################
# Subset CD14+ Cells #
######################
# Subsetting all cells for CD14
CD14.pos.subset <- subset(x = Immune_like_and_Neoplastic, subset = CD14 > 0)
# Export CD14 Expression Matrix for Stem-like Genes
CD14.pos.subset <- CD14.pos.subset[["RNA"]]$data
CD14.pos.subset <- as.matrix(CD14.pos.subset, 'sparseMatrix')
CD14.pos.subset <- subset(CD14.pos.subset, rownames(CD14.pos.subset) %in% c("NANOG", "SOX2", "CD44", "CXCR4", "CD70",  "ITGB3", "SNAI1", "PROCR", "ITGB1", "ZEB1", "CD34", "THY1", "GLI2", "NOTCH1", "GLI1", "BMI1", "EPCAM", "PROM1", "ITGA6", "POU5F1", "LGR5", "KIT", "ALDH1A1"))
CD14.pos.subset <- as.data.frame(rowMeans(CD14.pos.subset))
# Further organize and save data frame
CD14.pos.subset <-t(CD14.pos.subset)
write.csv(CD14.pos.subset, file = "/R/R_All_Breast_Tumors/All_Breast_Tumors_Output/CD14.All.pos.subset_stem_markers.csv")
rm(CD14.pos.subset)
gc()

#######################
# Subset PTPRC+ Cells #
#######################
# Subsetting all cells for PTPRC
PTPRC.pos.subset <- subset(x = Immune_like_and_Neoplastic, subset = PTPRC > 0)
# Export PTPRC Expression Matrix for Stem-like Genes
PTPRC.pos.subset <- PTPRC.pos.subset[["RNA"]]$data
PTPRC.pos.subset <- as.matrix(PTPRC.pos.subset, 'sparseMatrix')
PTPRC.pos.subset <- subset(PTPRC.pos.subset, rownames(PTPRC.pos.subset) %in% c("NANOG", "SOX2", "CD44", "CXCR4", "CD70",  "ITGB3", "SNAI1", "PROCR", "ITGB1", "ZEB1", "CD34", "THY1", "GLI2", "NOTCH1", "GLI1", "BMI1", "EPCAM", "PROM1", "ITGA6", "POU5F1", "LGR5", "KIT", "ALDH1A1"))
PTPRC.pos.subset <- as.data.frame(rowMeans(PTPRC.pos.subset))
# Further organize and save data frame
PTPRC.pos.subset <-t(PTPRC.pos.subset)
write.csv(PTPRC.pos.subset, file = "/R/R_All_Breast_Tumors/All_Breast_Tumors_Output/PTPRC.All.pos.subset_stem_markers.csv")
rm(PTPRC.pos.subset)
gc()

##########################
# Subset TNFSF13B+ Cells #
##########################
# Subsetting all cells for TNFSF13B
TNFSF13B.pos.subset <- subset(x = Immune_like_and_Neoplastic, subset = TNFSF13B > 0)
# Export TNFSF13B Expression Matrix for Stem-like Genes
TNFSF13B.pos.subset <- TNFSF13B.pos.subset[["RNA"]]$data
TNFSF13B.pos.subset <- as.matrix(TNFSF13B.pos.subset, 'sparseMatrix')
TNFSF13B.pos.subset <- subset(TNFSF13B.pos.subset, rownames(TNFSF13B.pos.subset) %in% c("NANOG", "SOX2", "CD44", "CXCR4", "CD70",  "ITGB3", "SNAI1", "PROCR", "ITGB1", "ZEB1", "CD34", "THY1", "GLI2", "NOTCH1", "GLI1", "BMI1", "EPCAM", "PROM1", "ITGA6", "POU5F1", "LGR5", "KIT", "ALDH1A1"))
TNFSF13B.pos.subset <- as.data.frame(rowMeans(TNFSF13B.pos.subset))
# Further organize and save data frame
TNFSF13B.pos.subset <-t(TNFSF13B.pos.subset)
write.csv(TNFSF13B.pos.subset, file = "/R/R_All_Breast_Tumors/All_Breast_Tumors_Output/TNFSF13B.All.pos.subset_stem_markers.csv")
rm(TNFSF13B.pos.subset)
gc()

######################
# Subset CD37+ Cells #
######################
# Subsetting all cells for CD37
CD37.pos.subset <- subset(x = Immune_like_and_Neoplastic, subset = CD37 > 0)
# Export CD37 Expression Matrix for Stem-like Genes
CD37.pos.subset <- CD37.pos.subset[["RNA"]]$data
CD37.pos.subset <- as.matrix(CD37.pos.subset, 'sparseMatrix')
CD37.pos.subset <- subset(CD37.pos.subset, rownames(CD37.pos.subset) %in% c("NANOG", "SOX2", "CD44", "CXCR4", "CD70",  "ITGB3", "SNAI1", "PROCR", "ITGB1", "ZEB1", "CD34", "THY1", "GLI2", "NOTCH1", "GLI1", "BMI1", "EPCAM", "PROM1", "ITGA6", "POU5F1", "LGR5", "KIT", "ALDH1A1"))
CD37.pos.subset <- as.data.frame(rowMeans(CD37.pos.subset))
# Further organize and save data frame
CD37.pos.subset <-t(CD37.pos.subset)
write.csv(CD37.pos.subset, file = "/R/R_All_Breast_Tumors/All_Breast_Tumors_Output/CD37.All.pos.subset_stem_markers.csv")
rm(CD37.pos.subset)
gc()

#######################
# Subset ITGB2+ Cells #
#######################
# Subsetting all cells for ITGB2
ITGB2.pos.subset <- subset(x = Immune_like_and_Neoplastic, subset = ITGB2 > 0)
# Export ITGB2 Expression Matrix for Stem-like Genes
ITGB2.pos.subset <- ITGB2.pos.subset[["RNA"]]$data
ITGB2.pos.subset <- as.matrix(ITGB2.pos.subset, 'sparseMatrix')
ITGB2.pos.subset <- subset(ITGB2.pos.subset, rownames(ITGB2.pos.subset) %in% c("NANOG", "SOX2", "CD44", "CXCR4", "CD70",  "ITGB3", "SNAI1", "PROCR", "ITGB1", "ZEB1", "CD34", "THY1", "GLI2", "NOTCH1", "GLI1", "BMI1", "EPCAM", "PROM1", "ITGA6", "POU5F1", "LGR5", "KIT", "ALDH1A1"))
ITGB2.pos.subset <- as.data.frame(rowMeans(ITGB2.pos.subset))
# Further organize and save data frame
ITGB2.pos.subset <-t(ITGB2.pos.subset)
write.csv(ITGB2.pos.subset, file = "/R/R_All_Breast_Tumors/All_Breast_Tumors_Output/ITGB2.All.pos.subset_stem_markers.csv")
rm(ITGB2.pos.subset)
gc()

#######################
# Subset CXCR4+ Cells #
#######################
# Subsetting all cells for CXCR4
CXCR4.pos.subset <- subset(x = Immune_like_and_Neoplastic, subset = CXCR4 > 0)
# Export CXCR4 Expression Matrix for Stem-like Genes
CXCR4.pos.subset <- CXCR4.pos.subset[["RNA"]]$data
CXCR4.pos.subset <- as.matrix(CXCR4.pos.subset, 'sparseMatrix')
CXCR4.pos.subset <- subset(CXCR4.pos.subset, rownames(CXCR4.pos.subset) %in% c("NANOG", "SOX2", "CD44", "CXCR4", "CD70",  "ITGB3", "SNAI1", "PROCR", "ITGB1", "ZEB1", "CD34", "THY1", "GLI2", "NOTCH1", "GLI1", "BMI1", "EPCAM", "PROM1", "ITGA6", "POU5F1", "LGR5", "KIT", "ALDH1A1"))
CXCR4.pos.subset <- as.data.frame(rowMeans(CXCR4.pos.subset))
# Further organize and save data frame
CXCR4.pos.subset <-t(CXCR4.pos.subset)
write.csv(CXCR4.pos.subset, file = "/R/R_All_Breast_Tumors/All_Breast_Tumors_Output/CXCR4.All.pos.subset_stem_markers.csv")
rm(CXCR4.pos.subset)
gc()

######################
# Subset CD74+ Cells #
######################
# Subsetting all cells for CD74
CD74.pos.subset <- subset(x = Immune_like_and_Neoplastic, subset = CD74 > 0)
# Export CD74 Expression Matrix for Stem-like Genes
CD74.pos.subset <- CD74.pos.subset[["RNA"]]$data
CD74.pos.subset <- as.matrix(CD74.pos.subset, 'sparseMatrix')
CD74.pos.subset <- subset(CD74.pos.subset, rownames(CD74.pos.subset) %in% c("NANOG", "SOX2", "CD44", "CXCR4", "CD70",  "ITGB3", "SNAI1", "PROCR", "ITGB1", "ZEB1", "CD34", "THY1", "GLI2", "NOTCH1", "GLI1", "BMI1", "EPCAM", "PROM1", "ITGA6", "POU5F1", "LGR5", "KIT", "ALDH1A1"))
CD74.pos.subset <- as.data.frame(rowMeans(CD74.pos.subset))
# Further organize and save data frame
CD74.pos.subset <-t(CD74.pos.subset)
write.csv(CD74.pos.subset, file = "/R/R_All_Breast_Tumors/All_Breast_Tumors_Output/CD74.All.pos.subset_stem_markers.csv")
rm(CD74.pos.subset)
gc()

########################
# Subset CLEC7A+ Cells #
########################
# Subsetting all cells for CLEC7A
CLEC7A.pos.subset <- subset(x = Immune_like_and_Neoplastic, subset = CLEC7A > 0)
# Export CLEC7A Expression Matrix for Stem-like Genes
CLEC7A.pos.subset <- CLEC7A.pos.subset[["RNA"]]$data
CLEC7A.pos.subset <- as.matrix(CLEC7A.pos.subset, 'sparseMatrix')
CLEC7A.pos.subset <- subset(CLEC7A.pos.subset, rownames(CLEC7A.pos.subset) %in% c("NANOG", "SOX2", "CD44", "CXCR4", "CD70",  "ITGB3", "SNAI1", "PROCR", "ITGB1", "ZEB1", "CD34", "THY1", "GLI2", "NOTCH1", "GLI1", "BMI1", "EPCAM", "PROM1", "ITGA6", "POU5F1", "LGR5", "KIT", "ALDH1A1"))
CLEC7A.pos.subset <- as.data.frame(rowMeans(CLEC7A.pos.subset))
# Further organize and save data frame
CLEC7A.pos.subset <-t(CLEC7A.pos.subset)
write.csv(CLEC7A.pos.subset, file = "/R/R_All_Breast_Tumors/All_Breast_Tumors_Output/CLEC7A.All.pos.subset_stem_markers.csv")
rm(CLEC7A.pos.subset)
gc()

######################
# Subset CD53+ Cells #
######################
# Subsetting all cells for CD53
CD53.pos.subset <- subset(x = Immune_like_and_Neoplastic, subset = CD53 > 0)
# Export CD53 Expression Matrix for Stem-like Genes
CD53.pos.subset <- CD53.pos.subset[["RNA"]]$data
CD53.pos.subset <- as.matrix(CD53.pos.subset, 'sparseMatrix')
CD53.pos.subset <- subset(CD53.pos.subset, rownames(CD53.pos.subset) %in% c("NANOG", "SOX2", "CD44", "CXCR4", "CD70",  "ITGB3", "SNAI1", "PROCR", "ITGB1", "ZEB1", "CD34", "THY1", "GLI2", "NOTCH1", "GLI1", "BMI1", "EPCAM", "PROM1", "ITGA6", "POU5F1", "LGR5", "KIT", "ALDH1A1"))
CD53.pos.subset <- as.data.frame(rowMeans(CD53.pos.subset))
# Further organize and save data frame
CD53.pos.subset <-t(CD53.pos.subset)
write.csv(CD53.pos.subset, file = "/R/R_All_Breast_Tumors/All_Breast_Tumors_Output/CD53.All.pos.subset_stem_markers.csv")
rm(CD53.pos.subset)
gc()

######################
# Subset CD83+ Cells #
######################
# Subsetting all cells for CD83
CD83.pos.subset <- subset(x = Immune_like_and_Neoplastic, subset = CD83 > 0)
# Export CD83 Expression Matrix for Stem-like Genes
CD83.pos.subset <- CD83.pos.subset[["RNA"]]$data
CD83.pos.subset <- as.matrix(CD83.pos.subset, 'sparseMatrix')
CD83.pos.subset <- subset(CD83.pos.subset, rownames(CD83.pos.subset) %in% c("NANOG", "SOX2", "CD44", "CXCR4", "CD70",  "ITGB3", "SNAI1", "PROCR", "ITGB1", "ZEB1", "CD34", "THY1", "GLI2", "NOTCH1", "GLI1", "BMI1", "EPCAM", "PROM1", "ITGA6", "POU5F1", "LGR5", "KIT", "ALDH1A1"))
CD83.pos.subset <- as.data.frame(rowMeans(CD83.pos.subset))
# Further organize and save data frame
CD83.pos.subset <-t(CD83.pos.subset)
write.csv(CD83.pos.subset, file = "/R/R_All_Breast_Tumors/All_Breast_Tumors_Output/CD83.All.pos.subset_stem_markers.csv")
rm(CD83.pos.subset)
gc()

#######################
# Subset PLAUR+ Cells #
#######################
# Subsetting all cells for PLAUR
PLAUR.pos.subset <- subset(x = Immune_like_and_Neoplastic, subset = PLAUR > 0)
# Export PLAUR Expression Matrix for Stem-like Genes
PLAUR.pos.subset <- PLAUR.pos.subset[["RNA"]]$data
PLAUR.pos.subset <- as.matrix(PLAUR.pos.subset, 'sparseMatrix')
PLAUR.pos.subset <- subset(PLAUR.pos.subset, rownames(PLAUR.pos.subset) %in% c("NANOG", "SOX2", "CD44", "CXCR4", "CD70",  "ITGB3", "SNAI1", "PROCR", "ITGB1", "ZEB1", "CD34", "THY1", "GLI2", "NOTCH1", "GLI1", "BMI1", "EPCAM", "PROM1", "ITGA6", "POU5F1", "LGR5", "KIT", "ALDH1A1"))
PLAUR.pos.subset <- as.data.frame(rowMeans(PLAUR.pos.subset))
# Further organize and save data frame
PLAUR.pos.subset <-t(PLAUR.pos.subset)
write.csv(PLAUR.pos.subset, file = "/R/R_All_Breast_Tumors/All_Breast_Tumors_Output/PLAUR.All.pos.subset_stem_markers.csv")
rm(PLAUR.pos.subset)
gc()

######################
# Subset CD44+ Cells #
######################
# Subsetting all cells for CD44
CD44.pos.subset <- subset(x = Immune_like_and_Neoplastic, subset = CD44 > 0)
# Export CD44 Expression Matrix for Stem-like Genes
CD44.pos.subset <- CD44.pos.subset[["RNA"]]$data
CD44.pos.subset <- as.matrix(CD44.pos.subset, 'sparseMatrix')
CD44.pos.subset <- subset(CD44.pos.subset, rownames(CD44.pos.subset) %in% c("NANOG", "SOX2", "CD44", "CXCR4", "CD70",  "ITGB3", "SNAI1", "PROCR", "ITGB1", "ZEB1", "CD34", "THY1", "GLI2", "NOTCH1", "GLI1", "BMI1", "EPCAM", "PROM1", "ITGA6", "POU5F1", "LGR5", "KIT", "ALDH1A1"))
CD44.pos.subset <- as.data.frame(rowMeans(CD44.pos.subset))
# Further organize and save data frame
CD44.pos.subset <-t(CD44.pos.subset)
write.csv(CD44.pos.subset, file = "/R/R_All_Breast_Tumors/All_Breast_Tumors_Output/CD44.All.pos.subset_stem_markers.csv")
rm(CD44.pos.subset)
gc()


##############################################
# Subset Neoplastic_cells_without_IM_Markers #
##############################################
# Export FCGR3A Expression Matrix for Stem-like Genes
Neoplastic_cells_without_IM_Markers <- Neoplastic_cells_without_IM_Markers[["RNA"]]$data
Neoplastic_cells_without_IM_Markers <- as.matrix(Neoplastic_cells_without_IM_Markers, 'sparseMatrix')
Neoplastic_cells_without_IM_Markers <- subset(Neoplastic_cells_without_IM_Markers, rownames(Neoplastic_cells_without_IM_Markers) %in% c("NANOG", "SOX2", "CD44", "CXCR4", "CD70",  "ITGB3", "SNAI1", "PROCR", "ITGB1", "ZEB1", "CD34", "THY1", "GLI2", "NOTCH1", "GLI1", "BMI1", "EPCAM", "PROM1", "ITGA6", "POU5F1", "LGR5", "KIT", "ALDH1A1"))
Neoplastic_cells_without_IM_Markers <- as.data.frame(rowMeans(Neoplastic_cells_without_IM_Markers))
# Further organize and save data frame
Neoplastic_cells_without_IM_Markers <-t(Neoplastic_cells_without_IM_Markers)
write.csv(Neoplastic_cells_without_IM_Markers, file = "/R/R_All_Breast_Tumors/All_Breast_Tumors_Output/Neoplastic_cells_without_IM_Markers_stem_markers.csv")
rm(Neoplastic_cells_without_IM_Markers)
gc()


###################################
# Concatenate Results for Heatmap #
###################################
# Load all files
file_paths <- c(
  "/R/R_All_Breast_Tumors/All_Breast_Tumors_Output/FCGR3A.All.pos.subset_stem_markers.csv",
  "/R/R_All_Breast_Tumors/All_Breast_Tumors_Output/CD2.All.pos.subset_stem_markers.csv",
  "/R/R_All_Breast_Tumors/All_Breast_Tumors_Output/CD3E.All.pos.subset_stem_markers.csv",
  "/R/R_All_Breast_Tumors/All_Breast_Tumors_Output/CD7.All.pos.subset_stem_markers.csv",
  "/R/R_All_Breast_Tumors/All_Breast_Tumors_Output/CD3G.All.pos.subset_stem_markers.csv",
  "/R/R_All_Breast_Tumors/All_Breast_Tumors_Output/CD3D.All.pos.subset_stem_markers.csv",
  "/R/R_All_Breast_Tumors/All_Breast_Tumors_Output/IL7R.All.pos.subset_stem_markers.csv",
  "/R/R_All_Breast_Tumors/All_Breast_Tumors_Output/FCGR2A.All.pos.subset_stem_markers.csv",
  "/R/R_All_Breast_Tumors/All_Breast_Tumors_Output/CD68.All.pos.subset_stem_markers.csv",
  "/R/R_All_Breast_Tumors/All_Breast_Tumors_Output/CD52.All.pos.subset_stem_markers.csv",
  "/R/R_All_Breast_Tumors/All_Breast_Tumors_Output/CD69.All.pos.subset_stem_markers.csv",
  "/R/R_All_Breast_Tumors/All_Breast_Tumors_Output/CD14.All.pos.subset_stem_markers.csv",
  "/R/R_All_Breast_Tumors/All_Breast_Tumors_Output/PTPRC.All.pos.subset_stem_markers.csv",
  "/R/R_All_Breast_Tumors/All_Breast_Tumors_Output/TNFSF13B.All.pos.subset_stem_markers.csv",
  "/R/R_All_Breast_Tumors/All_Breast_Tumors_Output/CD37.All.pos.subset_stem_markers.csv",
  "/R/R_All_Breast_Tumors/All_Breast_Tumors_Output/ITGB2.All.pos.subset_stem_markers.csv",
  "/R/R_All_Breast_Tumors/All_Breast_Tumors_Output/CXCR4.All.pos.subset_stem_markers.csv",
  "/R/R_All_Breast_Tumors/All_Breast_Tumors_Output/CD74.All.pos.subset_stem_markers.csv",
  "/R/R_All_Breast_Tumors/All_Breast_Tumors_Output/CLEC7A.All.pos.subset_stem_markers.csv",
  "/R/R_All_Breast_Tumors/All_Breast_Tumors_Output/CD53.All.pos.subset_stem_markers.csv",
  "/R/R_All_Breast_Tumors/All_Breast_Tumors_Output/CD83.All.pos.subset_stem_markers.csv",
  "/R/R_All_Breast_Tumors/All_Breast_Tumors_Output/PLAUR.All.pos.subset_stem_markers.csv",
  "/R/R_All_Breast_Tumors/All_Breast_Tumors_Output/CD44.All.pos.subset_stem_markers.csv",
  "/R/R_All_Breast_Tumors/All_Breast_Tumors_Output/Neoplastic_cells_without_IM_Markers_stem_markers.csv"
)


# Initialize an empty list to store data frames
IM_List <- list()

# Loop through each file path
for (file_path in file_paths) {
  # Read the CSV file
  IM_Compiled <- read.csv(file_path)
  
  # Extract the object name from the file name and create a new column
  object_name <- sub("(.subset_stem_markers|_without_IM_Markers_stem_markers)\\.csv$", "", basename(file_path))
  IM_Compiled <- cbind(object_name = object_name, IM_Compiled)
  
  # Append the data frame to the list
  IM_List[[length(IM_List) + 1]] <- IM_Compiled
}

# Combine all data frames into one
IM_Subsets_for_Heatmap_Stem_Markers <- do.call(rbind, IM_List)

# Print the combined data frame
print(IM_Subsets_for_Heatmap_Stem_Markers)

# Remove objects
rm(IM_List)
rm(IM_Compiled)
rm(file_path)
rm(file_paths)
rm(object_name)

# Delete Column "X"
IM_Subsets_for_Heatmap_Stem_Markers$X <- NULL

# Transpose data frame
IM_Subsets_for_Heatmap_Stem_Markers <- t(IM_Subsets_for_Heatmap_Stem_Markers)

# Set the first column as the new header
colnames(IM_Subsets_for_Heatmap_Stem_Markers) <- IM_Subsets_for_Heatmap_Stem_Markers[1, ]
# Remove the first row
IM_Subsets_for_Heatmap_Stem_Markers <- IM_Subsets_for_Heatmap_Stem_Markers[-1, ] 

# Save table
write.csv(IM_Subsets_for_Heatmap_Stem_Markers, "/R/R_All_Breast_Tumors/All_Breast_Tumors_Output/All.Neoplastic.Subset_for_Heatmap_Stem_Markers.csv")

#############################
##### Generate Heatmap  ##### 
#############################
# load packages
library("pheatmap")
library("RColorBrewer")
library("viridis")
library("dplyr")


# Read in data
IM_Subsets_for_Heatmap_Stem_Markers <- read.csv("/R/R_All_Breast_Tumors/All_Breast_Tumors_Output/All.Neoplastic.Subset_for_Heatmap_Stem_Markers.csv", header = TRUE, row.names = 1, sep = ",")
# Add column of all zeroes, done to preserve scaling; will be removed after
IM_Subsets_for_Heatmap_Stem_Markers$Zeroes <- 0

# Convert data to numeric format
IM_Subsets_for_Heatmap_Stem_Markers <- data.matrix(IM_Subsets_for_Heatmap_Stem_Markers)

# Transpose before scale
IM_Subsets_for_Heatmap_Stem_Markers <- t(IM_Subsets_for_Heatmap_Stem_Markers)
# Scale
IM_Subsets_for_Heatmap_Stem_Markers <- scale(IM_Subsets_for_Heatmap_Stem_Markers)
# Transpose after scale
IM_Subsets_for_Heatmap_Stem_Markers <- t(IM_Subsets_for_Heatmap_Stem_Markers)

# Remove zero column
IM_Subsets_for_Heatmap_Stem_Markers <- IM_Subsets_for_Heatmap_Stem_Markers[, -25]

# Change specific column names
colnames(IM_Subsets_for_Heatmap_Stem_Markers)[colnames(IM_Subsets_for_Heatmap_Stem_Markers) == "FCGR3A.All.pos"] <- "FCGR3A/CD16-pos"
colnames(IM_Subsets_for_Heatmap_Stem_Markers)[colnames(IM_Subsets_for_Heatmap_Stem_Markers) == "CD2.All.pos"] <- "CD2-pos"
colnames(IM_Subsets_for_Heatmap_Stem_Markers)[colnames(IM_Subsets_for_Heatmap_Stem_Markers) == "CD3E.All.pos"] <- "CD3E-pos"
colnames(IM_Subsets_for_Heatmap_Stem_Markers)[colnames(IM_Subsets_for_Heatmap_Stem_Markers) == "CD7.All.pos"] <- "CD7-pos"
colnames(IM_Subsets_for_Heatmap_Stem_Markers)[colnames(IM_Subsets_for_Heatmap_Stem_Markers) == "CD3G.All.pos"] <- "CD3G-pos"
colnames(IM_Subsets_for_Heatmap_Stem_Markers)[colnames(IM_Subsets_for_Heatmap_Stem_Markers) == "CD3D.All.pos"] <- "CD3D-pos"
colnames(IM_Subsets_for_Heatmap_Stem_Markers)[colnames(IM_Subsets_for_Heatmap_Stem_Markers) == "IL7R.All.pos"] <- "IL7R/CD127-pos"
colnames(IM_Subsets_for_Heatmap_Stem_Markers)[colnames(IM_Subsets_for_Heatmap_Stem_Markers) == "FCGR2A.All.pos"] <- "FCGR2A/CD32-pos"
colnames(IM_Subsets_for_Heatmap_Stem_Markers)[colnames(IM_Subsets_for_Heatmap_Stem_Markers) == "CD68.All.pos"] <- "CD68-pos"
colnames(IM_Subsets_for_Heatmap_Stem_Markers)[colnames(IM_Subsets_for_Heatmap_Stem_Markers) == "CD52.All.pos"] <- "CD52-pos"
colnames(IM_Subsets_for_Heatmap_Stem_Markers)[colnames(IM_Subsets_for_Heatmap_Stem_Markers) == "CD69.All.pos"] <- "CD69-pos"
colnames(IM_Subsets_for_Heatmap_Stem_Markers)[colnames(IM_Subsets_for_Heatmap_Stem_Markers) == "CD14.All.pos"] <- "CD14-pos"
colnames(IM_Subsets_for_Heatmap_Stem_Markers)[colnames(IM_Subsets_for_Heatmap_Stem_Markers) == "PTPRC.All.pos"] <- "CD45-pos"
colnames(IM_Subsets_for_Heatmap_Stem_Markers)[colnames(IM_Subsets_for_Heatmap_Stem_Markers) == "TNFSF13B.All.pos"] <- "TNFSF13B/CD257-pos"
colnames(IM_Subsets_for_Heatmap_Stem_Markers)[colnames(IM_Subsets_for_Heatmap_Stem_Markers) == "CD37.All.pos"] <- "CD37-pos"
colnames(IM_Subsets_for_Heatmap_Stem_Markers)[colnames(IM_Subsets_for_Heatmap_Stem_Markers) == "ITGB2.All.pos"] <- "ITGB2/CD18-pos"
colnames(IM_Subsets_for_Heatmap_Stem_Markers)[colnames(IM_Subsets_for_Heatmap_Stem_Markers) == "CXCR4.All.pos"] <- "CXCR4/CD184-pos"
colnames(IM_Subsets_for_Heatmap_Stem_Markers)[colnames(IM_Subsets_for_Heatmap_Stem_Markers) == "CD74.All.pos"] <- "CD74-pos"
colnames(IM_Subsets_for_Heatmap_Stem_Markers)[colnames(IM_Subsets_for_Heatmap_Stem_Markers) == "CLEC7A.All.pos"] <- "CLEC7A/CD369-pos"
colnames(IM_Subsets_for_Heatmap_Stem_Markers)[colnames(IM_Subsets_for_Heatmap_Stem_Markers) == "CD53.All.pos"] <- "CD53-pos"
colnames(IM_Subsets_for_Heatmap_Stem_Markers)[colnames(IM_Subsets_for_Heatmap_Stem_Markers) == "CD83.All.pos"] <- "CD83-pos"
colnames(IM_Subsets_for_Heatmap_Stem_Markers)[colnames(IM_Subsets_for_Heatmap_Stem_Markers) == "PLAUR.All.pos"] <- "PLAUR/CD87-pos"
colnames(IM_Subsets_for_Heatmap_Stem_Markers)[colnames(IM_Subsets_for_Heatmap_Stem_Markers) == "CD44.All.pos"] <- "CD44-pos"
colnames(IM_Subsets_for_Heatmap_Stem_Markers)[colnames(IM_Subsets_for_Heatmap_Stem_Markers) == "Neoplastic_cells"] <- "Other Neoplastic Cells"


# Change specific row names
rownames(IM_Subsets_for_Heatmap_Stem_Markers)[rownames(IM_Subsets_for_Heatmap_Stem_Markers) == "ITGB3"] <- "CD61"
rownames(IM_Subsets_for_Heatmap_Stem_Markers)[rownames(IM_Subsets_for_Heatmap_Stem_Markers) == "THY1"] <- "THY1/CD90"
rownames(IM_Subsets_for_Heatmap_Stem_Markers)[rownames(IM_Subsets_for_Heatmap_Stem_Markers) == "ITGB1"] <- "CD29"
rownames(IM_Subsets_for_Heatmap_Stem_Markers)[rownames(IM_Subsets_for_Heatmap_Stem_Markers) == "ITGA6"] <- "CD49f"
rownames(IM_Subsets_for_Heatmap_Stem_Markers)[rownames(IM_Subsets_for_Heatmap_Stem_Markers) == "POU5F1"] <- "OCT4"
rownames(IM_Subsets_for_Heatmap_Stem_Markers)[rownames(IM_Subsets_for_Heatmap_Stem_Markers) == "PROM1"] <- "CD133"




# Manually specify the order of the columns: Essentially redo the clustering while locking Neoplastic cells in the first row
IM_Subsets_for_Heatmap_Stem_Markers <- IM_Subsets_for_Heatmap_Stem_Markers[, c("Other Neoplastic Cells", 
                                                                               "CLEC7A/CD369-pos",
                                                                               "CD44-pos",
                                                                               "CD14-pos",
                                                                               "CXCR4/CD184-pos",
                                                                               "CD83-pos",
                                                                               "PLAUR/CD87-pos",
                                                                               "ITGB2/CD18-pos",
                                                                               "TNFSF13B/CD257-pos",
                                                                               "FCGR2A/CD32-pos",
                                                                               "CD68-pos",
                                                                               "CD74-pos",
                                                                               "FCGR3A/CD16-pos",
                                                                               "IL7R/CD127-pos",
                                                                               "CD3G-pos",
                                                                               "CD2-pos",
                                                                               "CD3E-pos",
                                                                               "CD7-pos",
                                                                               "CD3D-pos",
                                                                               "CD45-pos",
                                                                               "CD69-pos",
                                                                               "CD52-pos",
                                                                               "CD37-pos",
                                                                               "CD53-pos")]

# Set breaks & colors
the_breaks <- seq(-1, 1, by=0.25)
the_palette <- colorRampPalette(c("steelblue", "yellow3", "red4"))(8)

# Draw heatmap
IM_heatmap <- pheatmap(IM_Subsets_for_Heatmap_Stem_Markers, annotation_col = NA, show_colnames = TRUE, fontsize_row = 12, fontsize_col = 12,
                       breaks = the_breaks, color = the_palette, clustering_method = "median", cluster_rows = TRUE, 
                       cluster_cols = FALSE, treeheight_col = 0, scale = "none", cellwidth = 20, cellheight = 20,
                       main = "Stem/Progenitor-like Markers in Neoplastic Cells Expressing Immune Surface Receptors")

# Save
save_pheatmap_tif <- function(x, IM_heatmap, width=3000, height=3000, res = 300) {
  tiff(IM_heatmap, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

save_pheatmap_tif(IM_heatmap, "/R/R_All_Breast_Tumors/All_Breast_Tumors_Output/IM_Subsets_Stem_Markers_Heatmap.All.pos.subset.tif")












