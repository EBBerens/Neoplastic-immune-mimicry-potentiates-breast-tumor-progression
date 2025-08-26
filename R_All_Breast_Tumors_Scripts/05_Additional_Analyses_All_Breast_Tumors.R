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
# Step 26: Further evaluating inferred CNVs in CD45-pos vs CD45-neg Neoplastic Cells                                  #
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



#######################################################################################################################
# Step 26: Further evaluating inferred CNVs in CD45-pos vs CD45-neg Neoplastic Cells                                  #
#######################################################################################################################
# Part 1: InferCNV Group Analysis: CD45-neg vs CD45-pos                                                               #
# Part 2: CD45-pos Infer CNV Concordance & Coverage Analysis                                                          #
# Part 3: Patient Coverage Heatmap of CD45-pos CNV Regions Across Chromosomes                                         #
# Part 4: Within vs Between Patient CD45-pos CNV Region Correlations                                                  #
#######################################################################################################################



#################################################################################
# Part 1: InferCNV Group Analysis: CD45-neg vs CD45-pos (Streamlined)           #
#################################################################################

# Load libraries
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape2)

# Define output directory
output_dir <- "/R/R_All_Breast_Tumors/inferCNV_Analysis/CNV_Group_Analysis/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Dataset configuration
datasets <- list(
  visvader = list(
    name = "visvader",
    cd45_pos_path = "/R/R_Visvader/Visvader_RDS/RDS_Annotated/CD45.pos.Visvader_panCK_cnv.rds",
    cd45_neg_path = "/R/R_Visvader/Visvader_RDS/RDS_Annotated/CD45.neg.Visvader_panCK_cnv.rds",
    infercnv_base_path = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/"
  ),
  swarbrick = list(
    name = "swarbrick",
    cd45_pos_path = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Annotated/CD45.pos.Swarbrick_panCK_cnv.rds",
    cd45_neg_path = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Annotated/CD45.neg.Swarbrick_panCK_cnv.rds",
    infercnv_base_path = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/"
  ),
  pannthr = list(
    name = "pannthr",
    cd45_pos_path = "/R/R_PANNTHR/PANNTHR_RDS/RDS_Annotated/CD45.pos.PANNTHR_panCK_cnv.rds",
    cd45_neg_path = "/R/R_PANNTHR/PANNTHR_RDS/RDS_Annotated/CD45.neg.PANNTHR_panCK_cnv.rds",
    infercnv_base_path = "/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Output/"
  )
)

# Analysis parameters
MIN_CD45_POS <- 10
GROUPING_FILE <- "infercnv.17_HMM_predHMMi6.leiden.hmm_mode-subclusters.observation_groupings.txt"
REGIONS_FILE <- "HMM_CNV_predictions.HMMi6.leiden.hmm_mode-subclusters.Pnorm_0.5.pred_cnv_regions.dat"

# Create cell metadata
create_cell_metadata <- function(seurat_obj, cell_type, dataset_name) {
  data.frame(
    original_barcode = rownames(seurat_obj@meta.data),
    patient_id = seurat_obj@meta.data$orig.ident,
    cell_type = cell_type,
    dataset = dataset_name,
    stringsAsFactors = FALSE
  )
}

# Create patient folder mapping
create_patient_folder_mapping <- function(patient_ids, dataset_name, infercnv_base_path) {
  if (dataset_name == "visvader") {
    sapply(patient_ids, function(p) {
      base_id <- sub("_Total$", "", p)
      paste0(base_id, "_panCK")
    })
  } else if (dataset_name == "swarbrick") {
    result <- character(length(patient_ids))
    names(result) <- patient_ids
    all_dirs <- list.dirs(infercnv_base_path, full.names = FALSE, recursive = FALSE)
    
    for (i in seq_along(patient_ids)) {
      patient_id <- patient_ids[i]
      base_id <- sub("_Total$", "", patient_id)
      matching_dirs <- grep(paste0(base_id, ".*_panCK$"), all_dirs, value = TRUE)
      
      if (length(matching_dirs) >= 1) {
        result[i] <- matching_dirs[1]
      } else {
        result[i] <- paste0("NOTFOUND_", base_id, "_panCK")
      }
    }
    return(result)
  } else if (dataset_name == "pannthr") {
    sapply(patient_ids, function(p) {
      base_id <- sub("_Total$", "", p)
      paste0(base_id, "_panCK")
    })
  }
}

# Extract cell group assignments
extract_cell_group_assignments <- function(patient_dir, patient_id) {
  grouping_file <- file.path(patient_dir, GROUPING_FILE)
  
  if (!file.exists(grouping_file)) {
    return(NULL)
  }
  
  tryCatch({
    file_lines <- readLines(grouping_file, n = 10)
    first_line <- file_lines[1]
    has_headers <- grepl('"', first_line) && grepl('Color|Group', first_line)
    
    if (has_headers) {
      groupings <- read.table(grouping_file, header = FALSE, sep = "\t", 
                              stringsAsFactors = FALSE, quote = '"', skip = 1)
    } else {
      groupings <- read.table(grouping_file, header = FALSE, sep = "\t", 
                              stringsAsFactors = FALSE, quote = "")
    }
    
    if (ncol(groupings) == 1) {
      split_data <- strsplit(groupings$V1, "\\s+")
      if (length(split_data[[1]]) >= 2) {
        cell_barcodes <- sapply(split_data, function(x) x[1])
        cell_groups <- sapply(split_data, function(x) x[2])
        groupings <- data.frame(cell_barcode = cell_barcodes, cell_group_name = cell_groups, stringsAsFactors = FALSE)
      } else {
        return(NULL)
      }
    } else if (ncol(groupings) >= 2) {
      colnames(groupings)[1:2] <- c("cell_barcode", "cell_group_name")
      groupings <- groupings[, 1:2]
    } else {
      return(NULL)
    }
    
    groupings$cell_barcode <- gsub('"', '', groupings$cell_barcode)
    groupings$cell_group_name <- gsub('"', '', groupings$cell_group_name)
    groupings$cell_barcode <- gsub("-[0-9]+$", "", groupings$cell_barcode)
    
    valid_rows <- !is.na(groupings$cell_barcode) & !is.na(groupings$cell_group_name) & 
      groupings$cell_barcode != "" & groupings$cell_group_name != ""
    
    groupings <- groupings[valid_rows, ]
    if (nrow(groupings) == 0) return(NULL)
    
    return(groupings)
  }, error = function(e) {
    return(NULL)
  })
}

# Load CNV regions data
load_cnv_regions_data <- function(patient_dir) {
  regions_file <- file.path(patient_dir, REGIONS_FILE)
  
  if (!file.exists(regions_file)) {
    return(NULL)
  }
  
  tryCatch({
    regions_data <- read.table(regions_file, header = TRUE, sep = "\t", 
                               stringsAsFactors = FALSE, quote = "")
    
    required_cols <- c("cell_group_name", "cnv_name", "state", "chr", "start", "end")
    if (!all(required_cols %in% colnames(regions_data))) {
      return(NULL)
    }
    
    valid_chr_pattern <- "^chr([1-9]|1[0-9]|2[0-2]|X|Y)$"
    regions_data <- regions_data[grepl(valid_chr_pattern, regions_data$chr), ]
    
    if (nrow(regions_data) == 0) return(NULL)
    
    regions_data$cnv_call <- ifelse(regions_data$state <= 2, "deletion",
                                    ifelse(regions_data$state >= 4, "amplification", "neutral"))
    
    regions_data <- regions_data[regions_data$cnv_call != "neutral", ]
    if (nrow(regions_data) == 0) return(NULL)
    
    return(regions_data)
  }, error = function(e) {
    return(NULL)
  })
}

# Normalize group names
normalize_cnv_group_names <- function(group_names) {
  normalized <- character(length(group_names))
  
  for (i in seq_along(group_names)) {
    group_name <- group_names[i]
    
    if (grepl("Tumor\\.Tumor_s", group_name)) {
      parts <- strsplit(group_name, "\\.")[[1]]
      tumor_part <- parts[grepl("^Tumor_s", parts)]
      if (length(tumor_part) > 0) {
        normalized[i] <- tumor_part[1]
      } else {
        normalized[i] <- group_name
      }
    } else {
      normalized[i] <- group_name
    }
  }
  
  return(normalized)
}

# Extract CNV groups for cells (streamlined version)
extract_cnv_groups_for_cells <- function(cell_groupings, cell_barcodes, cell_type, patient_id) {
  if (is.null(cell_groupings) || length(cell_barcodes) == 0) {
    return(NULL)
  }
  
  seurat_base_barcodes <- gsub("-\\d+$", "", cell_barcodes)
  infercnv_barcodes <- cell_groupings$cell_barcode
  matched_indices <- which(infercnv_barcodes %in% seurat_base_barcodes)
  
  if (length(matched_indices) == 0) return(NULL)
  
  matched_groups <- cell_groupings$cell_group_name[matched_indices]
  if (length(matched_groups) == 0) return(NULL)
  
  group_counts <- table(matched_groups)
  
  result <- data.frame(
    patient_id = patient_id,
    cell_type = cell_type,
    cnv_group = names(group_counts),
    cell_count = as.numeric(group_counts),
    total_cells_matched = length(matched_indices),
    stringsAsFactors = FALSE
  )
  
  return(result)
}

# Process single patient (streamlined)
process_patient_cnv_groups <- function(patient_id, folder_name, patient_metadata, dataset_config) {
  patient_dir <- file.path(dataset_config$infercnv_base_path, folder_name)
  
  if (!dir.exists(patient_dir)) {
    return(NULL)
  }
  
  patient_cells <- patient_metadata %>% filter(patient_id == !!patient_id)
  pos_cells <- patient_cells$original_barcode[patient_cells$cell_type == "CD45-pos"]
  neg_cells <- patient_cells$original_barcode[patient_cells$cell_type == "CD45-neg"]
  
  if (length(pos_cells) < MIN_CD45_POS) {
    return(NULL)
  }
  
  cell_groupings <- extract_cell_group_assignments(patient_dir, patient_id)
  if (is.null(cell_groupings)) {
    return(NULL)
  }
  
  pos_groups <- extract_cnv_groups_for_cells(cell_groupings, pos_cells, "CD45-pos", patient_id)
  neg_groups <- extract_cnv_groups_for_cells(cell_groupings, neg_cells, "CD45-neg", patient_id)
  
  if (is.null(pos_groups) && is.null(neg_groups)) {
    return(NULL)
  }
  
  all_groups <- rbind(pos_groups, neg_groups)
  if (is.null(all_groups) || nrow(all_groups) == 0) {
    return(NULL)
  }
  
  all_groups$dataset <- dataset_config$name
  
  return(list(group_data = all_groups))
}

# Analyze patient group overlap (streamlined)
analyze_patient_group_overlap <- function(patient_groups) {
  if (is.null(patient_groups) || nrow(patient_groups) == 0) {
    return(NULL)
  }
  
  patient_id <- unique(patient_groups$patient_id)[1]
  dataset <- unique(patient_groups$dataset)[1]
  
  pos_groups <- patient_groups[patient_groups$cell_type == "CD45-pos", ]
  neg_groups <- patient_groups[patient_groups$cell_type == "CD45-neg", ]
  
  n_pos_groups <- nrow(pos_groups)
  n_neg_groups <- nrow(neg_groups)
  
  if (n_pos_groups > 0 && n_neg_groups > 0) {
    pos_group_names <- pos_groups$cnv_group
    neg_group_names <- neg_groups$cnv_group
    overlapping_groups <- intersect(pos_group_names, neg_group_names)
    n_overlapping <- length(overlapping_groups)
    total_pos_cells <- sum(pos_groups$cell_count)
    total_neg_cells <- sum(neg_groups$cell_count)
  } else {
    overlapping_groups <- character(0)
    n_overlapping <- 0
    total_pos_cells <- if (n_pos_groups > 0) sum(pos_groups$cell_count) else 0
    total_neg_cells <- if (n_neg_groups > 0) sum(neg_groups$cell_count) else 0
  }
  
  result <- data.frame(
    patient_id = patient_id,
    dataset = dataset,
    n_pos_groups = n_pos_groups,
    n_neg_groups = n_neg_groups,
    n_overlapping_groups = n_overlapping,
    total_pos_cells = total_pos_cells,
    total_neg_cells = total_neg_cells,
    overlap_fraction_of_pos = if (n_pos_groups > 0) n_overlapping / n_pos_groups else 0,
    overlap_fraction_of_neg = if (n_neg_groups > 0) n_overlapping / n_neg_groups else 0,
    stringsAsFactors = FALSE
  )
  
  return(result)
}

# Analyze dataset (streamlined)
analyze_dataset_cnv_groups <- function(dataset_config) {
  cd45_pos <- readRDS(dataset_config$cd45_pos_path)
  cd45_neg <- readRDS(dataset_config$cd45_neg_path)
  
  cell_metadata_pos <- create_cell_metadata(cd45_pos, "CD45-pos", dataset_config$name)
  cell_metadata_neg <- create_cell_metadata(cd45_neg, "CD45-neg", dataset_config$name)
  all_cell_metadata <- rbind(cell_metadata_pos, cell_metadata_neg)
  
  patient_ids_pos <- unique(cd45_pos@meta.data$orig.ident)
  patient_ids_neg <- unique(cd45_neg@meta.data$orig.ident)
  patient_ids <- intersect(patient_ids_pos, patient_ids_neg)
  
  rm(cd45_pos, cd45_neg)
  gc()
  
  patient_to_folder <- create_patient_folder_mapping(patient_ids, dataset_config$name, dataset_config$infercnv_base_path)
  
  suitable_patients <- character()
  for (patient in patient_ids) {
    pos_count <- sum(all_cell_metadata$patient_id == patient & all_cell_metadata$cell_type == "CD45-pos")
    if (pos_count >= MIN_CD45_POS) {
      suitable_patients <- c(suitable_patients, patient)
    }
  }
  
  if (length(suitable_patients) == 0) {
    return(NULL)
  }
  
  all_group_data <- list()
  all_overlap_data <- list()
  
  for (patient_id in suitable_patients) {
    folder_name <- patient_to_folder[patient_id]
    full_path <- file.path(dataset_config$infercnv_base_path, folder_name)
    if (!dir.exists(full_path)) next
    
    tryCatch({
      patient_result <- process_patient_cnv_groups(patient_id, folder_name, all_cell_metadata, dataset_config)
      
      if (!is.null(patient_result) && !is.null(patient_result$group_data)) {
        all_group_data[[patient_id]] <- patient_result$group_data
        
        overlap_result <- analyze_patient_group_overlap(patient_result$group_data)
        if (!is.null(overlap_result)) {
          all_overlap_data[[patient_id]] <- overlap_result
        }
      }
    }, error = function(e) {
      NULL
    })
  }
  
  if (length(all_group_data) > 0) {
    combined_groups <- do.call(rbind, all_group_data)
    combined_overlap <- if (length(all_overlap_data) > 0) do.call(rbind, all_overlap_data) else NULL
    
    return(list(
      group_data = combined_groups,
      overlap_data = combined_overlap
    ))
  }
  
  return(NULL)
}

# Create streamlined visualizations (3 plots only)
create_cnv_group_visualizations <- function(overlap_data, group_data) {
  if (nrow(overlap_data) == 0) {
    return(NULL)
  }
  
  tiff_file <- file.path(output_dir, "CNV_Group_Analysis_Plots.tiff")
  
  tryCatch({
    tiff(tiff_file, width = 16, height = 12, units = "in", res = 300)
    par(mfrow = c(2, 2), mar = c(5, 4, 3, 2))
    
    # Plot 1: Cell Counts per Patient
    cell_count_data <- overlap_data %>%
      select(patient_id, total_pos_cells, total_neg_cells) %>%
      gather(key = "cell_type", value = "cell_count", -patient_id) %>%
      mutate(cell_type = ifelse(cell_type == "total_pos_cells", "CD45-pos", "CD45-neg"))
    
    valid_cell_count_data <- cell_count_data %>%
      filter(!is.na(cell_count) & cell_count > 0)
    
    if (nrow(valid_cell_count_data) > 0) {
      boxplot(cell_count ~ cell_type, data = valid_cell_count_data,
              main = "Cell Counts per Patient",
              ylab = "Number of Cells",
              col = c("lightcoral", "lightblue"),
              log = "y")
    }
    
    # Plot 2: CNV Groups per Patient
    group_counts <- overlap_data %>%
      select(patient_id, n_pos_groups, n_neg_groups) %>%
      gather(key = "cell_type", value = "n_groups", -patient_id) %>%
      mutate(cell_type = ifelse(cell_type == "n_pos_groups", "CD45-pos", "CD45-neg"))
    
    boxplot(n_groups ~ cell_type, data = group_counts,
            main = "CNV Groups per Patient",
            ylab = "Number of CNV Groups",
            col = c("lightcoral", "lightblue"))
    
    # Plot 3: CNV Group Overlap Fractions
    boxplot(list(
      "Fraction of CD45- groups" = overlap_data$overlap_fraction_of_neg,
      "Fraction of CD45+ groups" = overlap_data$overlap_fraction_of_pos
    ),
    main = "CNV Group Overlap Fractions",
    ylab = "Overlap Fraction",
    col = c("orange", "yellow"))
    
    dev.off()
    return(tiff_file)
  }, error = function(e) {
    if (dev.cur() != 1) dev.off()
    return(NULL)
  })
}

# Main analysis function (streamlined)
main_cnv_group_analysis <- function() {
  all_group_data <- data.frame()
  all_overlap_data <- data.frame()
  
  for (dataset_name in names(datasets)) {
    dataset_config <- datasets[[dataset_name]]
    
    if (file.exists(dataset_config$cd45_pos_path) && 
        file.exists(dataset_config$cd45_neg_path) &&
        dir.exists(dataset_config$infercnv_base_path)) {
      
      result <- analyze_dataset_cnv_groups(dataset_config)
      if (!is.null(result)) {
        all_group_data <- rbind(all_group_data, result$group_data)
        if (!is.null(result$overlap_data)) {
          all_overlap_data <- rbind(all_overlap_data, result$overlap_data)
        }
      }
    }
  }
  
  return(list(
    group_data = all_group_data, 
    overlap_data = all_overlap_data
  ))
}

# Save streamlined analysis data
save_analysis_data <- function(results) {
  if (nrow(results$group_data) > 0) {
    write.csv(results$group_data, 
              file.path(output_dir, "cnv_group_data.csv"), 
              row.names = FALSE)
  }
  
  if (nrow(results$overlap_data) > 0) {
    write.csv(results$overlap_data, 
              file.path(output_dir, "patient_overlap_data.csv"), 
              row.names = FALSE)
  }
  
  # Create plot data for external recreation (streamlined)
  if (nrow(results$overlap_data) > 0) {
    # Cell counts data
    cell_count_data <- results$overlap_data %>%
      select(patient_id, dataset, total_pos_cells, total_neg_cells) %>%
      gather(key = "cell_type", value = "cell_count", -patient_id, -dataset) %>%
      mutate(cell_type = ifelse(cell_type == "total_pos_cells", "CD45-pos", "CD45-neg"))
    
    write.csv(cell_count_data, 
              file.path(output_dir, "plot_data_cell_counts.csv"), 
              row.names = FALSE)
    
    # Group counts data
    group_counts_data <- results$overlap_data %>%
      select(patient_id, dataset, n_pos_groups, n_neg_groups) %>%
      gather(key = "cell_type", value = "n_groups", -patient_id, -dataset) %>%
      mutate(cell_type = ifelse(cell_type == "n_pos_groups", "CD45-pos", "CD45-neg"))
    
    write.csv(group_counts_data, 
              file.path(output_dir, "plot_data_group_counts.csv"), 
              row.names = FALSE)
    
    # Overlap fractions data
    overlap_fractions_data <- results$overlap_data %>%
      select(patient_id, dataset, overlap_fraction_of_pos, overlap_fraction_of_neg) %>%
      gather(key = "fraction_type", value = "overlap_fraction", -patient_id, -dataset) %>%
      mutate(fraction_type = ifelse(fraction_type == "overlap_fraction_of_pos", 
                                    "Fraction of CD45+ groups", "Fraction of CD45- groups"))
    
    write.csv(overlap_fractions_data, 
              file.path(output_dir, "plot_data_overlap_fractions.csv"), 
              row.names = FALSE)
  }
}

# Execute analysis
results <- main_cnv_group_analysis()

if (nrow(results$group_data) > 0) {
  create_cnv_group_visualizations(results$overlap_data, results$group_data)
  save_analysis_data(results)
}







#################################################################################
# Part 2: CD45-pos Infer CNV Concordance & Coverage Analysis                    #
#################################################################################

# Load Libraries
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)

# Define output directory
output_dir <- "/R/R_All_Breast_Tumors/inferCNV_Analysis/CD45_CNV_Concordance_Analysis/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Dataset configuration
datasets <- list(
  visvader = list(
    name = "visvader",
    cd45_pos_path = "/R/R_Visvader/Visvader_RDS/RDS_Annotated/CD45.pos.Visvader_panCK_cnv.rds",
    cd45_neg_path = "/R/R_Visvader/Visvader_RDS/RDS_Annotated/CD45.neg.Visvader_panCK_cnv.rds",
    infercnv_base_path = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output"
  ),
  swarbrick = list(
    name = "swarbrick",
    cd45_pos_path = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Annotated/CD45.pos.Swarbrick_panCK_cnv.rds",
    cd45_neg_path = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Annotated/CD45.neg.Swarbrick_panCK_cnv.rds",
    infercnv_base_path = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output"
  ),
  pannthr = list(
    name = "pannthr",
    cd45_pos_path = "/R/R_PANNTHR/PANNTHR_RDS/RDS_Annotated/CD45.pos.PANNTHR_panCK_cnv.rds",
    cd45_neg_path = "/R/R_PANNTHR/PANNTHR_RDS/RDS_Annotated/CD45.neg.PANNTHR_panCK_cnv.rds",
    infercnv_base_path = "/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Output"
  )
)

# Analysis parameters
MIN_CD45_POS <- 10
REGIONS_FILE <- "HMM_CNV_predictions.HMMi6.leiden.hmm_mode-subclusters.Pnorm_0.5.pred_cnv_regions.dat"
GENES_FILE <- "HMM_CNV_predictions.HMMi6.leiden.hmm_mode-subclusters.Pnorm_0.5.pred_cnv_genes.dat"
GROUPING_FILE <- "infercnv.17_HMM_predHMMi6.leiden.hmm_mode-subclusters.observation_groupings.txt"

# Function to create cell metadata
create_cell_metadata <- function(seurat_obj, cell_type, dataset_name) {
  original_barcodes <- rownames(seurat_obj@meta.data)
  patient_ids <- seurat_obj@meta.data$orig.ident
  
  metadata <- data.frame(
    original_barcode = original_barcodes,
    patient_id = patient_ids,
    cell_type = cell_type,
    dataset = dataset_name,
    stringsAsFactors = FALSE
  )
  
  return(metadata)
}

# Function to map HMM states to CNV calls
map_hmm_states_to_cnv <- function(state_vector) {
  result <- character(length(state_vector))
  numeric_states <- as.numeric(state_vector)
  
  result[numeric_states <= 2] <- "deletion"
  result[numeric_states == 3] <- "neutral"
  result[numeric_states >= 4] <- "amplification"
  result[is.na(numeric_states)] <- "neutral"
  
  return(result)
}

# Function to normalize group names for matching
normalize_group_names <- function(group_names) {
  normalized <- character(length(group_names))
  
  for (i in seq_along(group_names)) {
    group_name <- group_names[i]
    
    if (grepl("Tumor\\.Tumor_s", group_name)) {
      parts <- strsplit(group_name, "\\.")[[1]]
      tumor_part <- parts[grepl("^Tumor_s", parts)]
      if (length(tumor_part) > 0) {
        normalized[i] <- tumor_part[1]
      } else {
        normalized[i] <- group_name
      }
    } else {
      normalized[i] <- group_name
    }
  }
  
  return(normalized)
}

# Function to create patient folder mapping
create_patient_folder_mapping <- function(patient_ids, dataset_name, infercnv_base_path) {
  if (dataset_name == "visvader") {
    sapply(patient_ids, function(p) {
      base_id <- sub("_Total$", "", p)
      paste0(base_id, "_panCK")
    })
  } else if (dataset_name == "swarbrick") {
    result <- character(length(patient_ids))
    names(result) <- patient_ids
    
    all_dirs <- list.dirs(infercnv_base_path, full.names = FALSE, recursive = FALSE)
    
    for (i in seq_along(patient_ids)) {
      patient_id <- patient_ids[i]
      base_id <- sub("_Total$", "", patient_id)
      
      matching_dirs <- grep(paste0(base_id, ".*_panCK$"), all_dirs, value = TRUE)
      
      if (length(matching_dirs) >= 1) {
        result[i] <- matching_dirs[1]
      } else {
        result[i] <- paste0("NOTFOUND_", base_id, "_panCK")
      }
    }
    
    return(result)
  } else if (dataset_name == "pannthr") {
    sapply(patient_ids, function(p) {
      base_id <- sub("_Total$", "", p)
      paste0(base_id, "_panCK")
    })
  }
}

# Function to extract cell group assignments
extract_cell_group_assignments <- function(patient_dir, patient_id) {
  grouping_file_path <- file.path(patient_dir, GROUPING_FILE)
  
  if (!file.exists(grouping_file_path)) {
    return(NULL)
  }
  
  tryCatch({
    raw_lines <- readLines(grouping_file_path)
    
    if (length(raw_lines) < 2) {
      return(NULL)
    }
    
    data_lines <- raw_lines[-1]
    parsed_data <- list()
    
    for (i in 1:length(data_lines)) {
      line <- data_lines[i]
      
      if (grepl('^"', line)) {
        parts <- strsplit(line, " ")[[1]]
        parts <- gsub('"', '', parts)
        parts <- parts[parts != ""]
        
        if (length(parts) >= 2) {
          barcode <- parts[1]
          group <- parts[2]
          
          if (barcode != "" && group != "") {
            parsed_data[[length(parsed_data) + 1]] <- data.frame(
              cell_barcode = barcode,
              cell_group_name = group,
              stringsAsFactors = FALSE
            )
          }
        }
      }
    }
    
    if (length(parsed_data) == 0) {
      return(NULL)
    }
    
    groupings <- do.call(rbind, parsed_data)
    groupings <- groupings[!is.na(groupings$cell_barcode) & groupings$cell_barcode != "", ]
    groupings <- groupings[!is.na(groupings$cell_group_name) & groupings$cell_group_name != "", ]
    
    return(groupings)
    
  }, error = function(e) {
    return(NULL)
  })
}

# Function to load CNV data
load_cnv_data <- function(patient_dir, data_type) {
  file_name <- if (data_type == "regions") REGIONS_FILE else GENES_FILE
  file_path <- file.path(patient_dir, file_name)
  
  if (!file.exists(file_path)) {
    return(NULL)
  }
  
  tryCatch({
    data <- read.table(file_path, header = TRUE, sep = "\t", 
                       stringsAsFactors = FALSE, quote = "", check.names = FALSE)
    return(data)
  }, error = function(e) {
    return(NULL)
  })
}

# Function to extract CNV data for specific cells
extract_cell_cnv_data <- function(cnv_data, cell_groupings, cell_barcodes, cell_type, patient_id, data_type) {
  if (is.null(cnv_data) || is.null(cell_groupings)) {
    return(NULL)
  }
  
  matched_cell_groups <- cell_groupings$cell_group_name[cell_groupings$cell_barcode %in% cell_barcodes]
  
  if (length(matched_cell_groups) == 0) {
    return(NULL)
  }
  
  unique_matched_groups <- unique(matched_cell_groups)
  cnv_data$normalized_group_name <- normalize_group_names(cnv_data$cell_group_name)
  
  cell_cnv_data <- cnv_data %>%
    filter(normalized_group_name %in% unique_matched_groups) %>%
    mutate(
      cnv_call = map_hmm_states_to_cnv(state),
      cell_type = cell_type,
      patient_id = patient_id,
      data_type = data_type
    ) %>%
    filter(cnv_call != "neutral")
  
  if (data_type == "regions") {
    cell_cnv_data$region_size <- cell_cnv_data$end - cell_cnv_data$start + 1
    cell_cnv_data$region_id <- paste0(cell_cnv_data$chr, ":", cell_cnv_data$start, "-", cell_cnv_data$end)
  } else if (data_type == "genes") {
    cell_cnv_data$gene_id <- paste0(cell_cnv_data$chr, ":", cell_cnv_data$gene)
  }
  
  return(cell_cnv_data)
}

# Function to calculate concordance and coverage metrics
calculate_cnv_metrics <- function(pos_cnv_data, neg_cnv_data, patient_id, data_type) {
  if (is.null(pos_cnv_data) || is.null(neg_cnv_data) || 
      nrow(pos_cnv_data) == 0 || nrow(neg_cnv_data) == 0) {
    return(NULL)
  }
  
  # Determine the identifier column based on data type
  if (data_type == "regions") {
    pos_items <- unique(pos_cnv_data$region_id)
    neg_items <- unique(neg_cnv_data$region_id)
    id_col <- "region_id"
  } else if (data_type == "genes") {
    pos_items <- unique(pos_cnv_data$gene_id)
    neg_items <- unique(neg_cnv_data$gene_id)
    id_col <- "gene_id"
  } else {
    return(NULL)
  }
  
  # Calculate basic set operations
  overlap_items <- intersect(pos_items, neg_items)
  pos_specific_items <- setdiff(pos_items, neg_items)
  neg_specific_items <- setdiff(neg_items, pos_items)
  
  # Concordance for overlapping items
  concordant_count <- 0
  total_overlapping <- length(overlap_items)
  
  directional_concordance <- if (total_overlapping > 0) {
    for (item in overlap_items) {
      pos_calls <- pos_cnv_data$cnv_call[pos_cnv_data[[id_col]] == item]
      neg_calls <- neg_cnv_data$cnv_call[neg_cnv_data[[id_col]] == item]
      
      pos_call <- names(sort(table(pos_calls), decreasing = TRUE))[1]
      neg_call <- names(sort(table(neg_calls), decreasing = TRUE))[1]
      
      if (pos_call == neg_call) {
        concordant_count <- concordant_count + 1
      }
    }
    
    concordant_count / total_overlapping
  } else {
    NA
  }
  
  # Coverage metrics
  pos_coverage_by_neg <- if (length(pos_items) > 0) {
    length(overlap_items) / length(pos_items)
  } else {
    0
  }
  
  neg_coverage_by_pos <- if (length(neg_items) > 0) {
    length(overlap_items) / length(neg_items)
  } else {
    0
  }
  
  # CNV type-specific concordance
  deletion_concordance <- calculate_cnv_type_concordance(
    pos_cnv_data, neg_cnv_data, overlap_items, "deletion", id_col
  )
  
  amplification_concordance <- calculate_cnv_type_concordance(
    pos_cnv_data, neg_cnv_data, overlap_items, "amplification", id_col
  )
  
  results <- list(
    directional_concordance = directional_concordance,
    pos_coverage_by_neg = pos_coverage_by_neg,
    neg_coverage_by_pos = neg_coverage_by_pos,
    deletion_concordance = deletion_concordance,
    amplification_concordance = amplification_concordance,
    total_overlapping_items = total_overlapping,
    pos_total_items = length(pos_items),
    neg_total_items = length(neg_items),
    pos_specific_items = length(pos_specific_items),
    neg_specific_items = length(neg_specific_items),
    concordant_items = concordant_count
  )
  
  return(results)
}

# Helper function to calculate CNV type-specific concordance
calculate_cnv_type_concordance <- function(pos_cnv_data, neg_cnv_data, overlap_items, cnv_type, id_col) {
  if (length(overlap_items) == 0) {
    return(NA)
  }
  
  concordant_count <- 0
  relevant_overlaps <- 0
  
  for (item in overlap_items) {
    pos_calls <- pos_cnv_data$cnv_call[pos_cnv_data[[id_col]] == item]
    neg_calls <- neg_cnv_data$cnv_call[neg_cnv_data[[id_col]] == item]
    
    pos_has_type <- any(pos_calls == cnv_type)
    neg_has_type <- any(neg_calls == cnv_type)
    
    if (pos_has_type || neg_has_type) {
      relevant_overlaps <- relevant_overlaps + 1
      
      if (pos_has_type && neg_has_type) {
        concordant_count <- concordant_count + 1
      }
    }
  }
  
  if (relevant_overlaps > 0) {
    return(concordant_count / relevant_overlaps)
  } else {
    return(NA)
  }
}

# Function to compare CNVs between populations
compare_cnv_populations <- function(pos_cnv_data, neg_cnv_data, patient_id, data_type) {
  if (is.null(pos_cnv_data) || is.null(neg_cnv_data) || 
      nrow(pos_cnv_data) == 0 || nrow(neg_cnv_data) == 0) {
    return(NULL)
  }
  
  metrics <- calculate_cnv_metrics(pos_cnv_data, neg_cnv_data, patient_id, data_type)
  
  if (is.null(metrics)) {
    return(NULL)
  }
  
  summary_metrics <- list(
    patient_id = patient_id,
    data_type = data_type,
    directional_concordance = metrics$directional_concordance,
    pos_coverage_by_neg = metrics$pos_coverage_by_neg,
    neg_coverage_by_pos = metrics$neg_coverage_by_pos,
    deletion_concordance = metrics$deletion_concordance,
    amplification_concordance = metrics$amplification_concordance,
    pos_total_cnvs = nrow(pos_cnv_data),
    neg_total_cnvs = nrow(neg_cnv_data),
    pos_deletions = sum(pos_cnv_data$cnv_call == "deletion"),
    pos_amplifications = sum(pos_cnv_data$cnv_call == "amplification"),
    neg_deletions = sum(neg_cnv_data$cnv_call == "deletion"),
    neg_amplifications = sum(neg_cnv_data$cnv_call == "amplification"),
    total_overlapping_items = metrics$total_overlapping_items,
    pos_total_items = metrics$pos_total_items,
    neg_total_items = metrics$neg_total_items,
    pos_specific_items = metrics$pos_specific_items,
    neg_specific_items = metrics$neg_specific_items,
    concordant_items = metrics$concordant_items
  )
  
  if (data_type == "regions") {
    summary_metrics$pos_mean_region_size <- if(nrow(pos_cnv_data) > 0) mean(pos_cnv_data$region_size) else 0
    summary_metrics$neg_mean_region_size <- if(nrow(neg_cnv_data) > 0) mean(neg_cnv_data$region_size) else 0
    summary_metrics$pos_total_genomic_size <- if(nrow(pos_cnv_data) > 0) sum(pos_cnv_data$region_size) else 0
    summary_metrics$neg_total_genomic_size <- if(nrow(neg_cnv_data) > 0) sum(neg_cnv_data$region_size) else 0
  }
  
  summary_metrics <- lapply(summary_metrics, function(x) {
    if (is.numeric(x) && (is.nan(x) || is.infinite(x))) return(0)
    return(x)
  })
  
  return(list(summary_metrics = summary_metrics))
}

# Function to filter patients by cell count
filter_patients_by_cell_count <- function(patient_ids, all_cell_metadata) {
  suitable_patients <- character()
  
  for (patient_id in patient_ids) {
    pos_count <- sum(all_cell_metadata$patient_id == patient_id & 
                       all_cell_metadata$cell_type == "CD45-pos")
    
    if (pos_count >= MIN_CD45_POS) {
      suitable_patients <- c(suitable_patients, patient_id)
    }
  }
  
  return(suitable_patients)
}

# Function to process patient
process_patient <- function(patient_id, folder_name, patient_metadata, dataset_config) {
  patient_dir <- file.path(dataset_config$infercnv_base_path, folder_name)
  
  if (!dir.exists(patient_dir)) {
    return(NULL)
  }
  
  patient_cells <- patient_metadata %>% filter(patient_id == !!patient_id)
  pos_cells <- patient_cells$original_barcode[patient_cells$cell_type == "CD45-pos"]
  neg_cells <- patient_cells$original_barcode[patient_cells$cell_type == "CD45-neg"]
  
  if (length(pos_cells) < MIN_CD45_POS) {
    return(NULL)
  }
  
  cell_groupings <- extract_cell_group_assignments(patient_dir, patient_id)
  if (is.null(cell_groupings)) {
    return(NULL)
  }
  
  infercnv_barcodes <- cell_groupings$cell_barcode
  pos_matches <- intersect(pos_cells, infercnv_barcodes)
  neg_matches <- intersect(neg_cells, infercnv_barcodes)
  
  if (length(pos_matches) < MIN_CD45_POS) {
    return(NULL)
  }
  
  results <- list()
  
  for (data_type in c("regions", "genes")) {
    cnv_data <- load_cnv_data(patient_dir, data_type)
    if (is.null(cnv_data)) {
      next
    }
    
    pos_cnv_data <- extract_cell_cnv_data(cnv_data, cell_groupings, pos_matches, "CD45_pos", patient_id, data_type)
    neg_cnv_data <- extract_cell_cnv_data(cnv_data, cell_groupings, neg_matches, "CD45_neg", patient_id, data_type)
    
    if (is.null(pos_cnv_data) || is.null(neg_cnv_data) || 
        nrow(pos_cnv_data) == 0 || nrow(neg_cnv_data) == 0) {
      next
    }
    
    comparison_results <- compare_cnv_populations(pos_cnv_data, neg_cnv_data, patient_id, data_type)
    
    if (!is.null(comparison_results)) {
      comparison_results$summary_metrics$dataset <- dataset_config$name
      comparison_results$summary_metrics$pos_cells_matched <- length(pos_matches)
      comparison_results$summary_metrics$neg_cells_matched <- length(neg_matches)
      
      results[[data_type]] <- comparison_results
    }
  }
  
  return(if(length(results) > 0) results else NULL)
}

# Function to analyze dataset
analyze_dataset <- function(dataset_config) {
  dataset_name <- dataset_config$name
  
  cd45_pos <- readRDS(dataset_config$cd45_pos_path)
  cd45_neg <- readRDS(dataset_config$cd45_neg_path)
  
  cell_metadata_pos <- create_cell_metadata(cd45_pos, "CD45-pos", dataset_name)
  cell_metadata_neg <- create_cell_metadata(cd45_neg, "CD45-neg", dataset_name)
  all_cell_metadata <- rbind(cell_metadata_pos, cell_metadata_neg)
  
  patient_ids_pos <- unique(cd45_pos@meta.data$orig.ident)
  patient_ids_neg <- unique(cd45_neg@meta.data$orig.ident)
  patient_ids <- intersect(patient_ids_pos, patient_ids_neg)
  
  rm(cd45_pos, cd45_neg)
  gc()
  
  patient_to_folder <- create_patient_folder_mapping(patient_ids, dataset_name, dataset_config$infercnv_base_path)
  
  suitable_patients <- filter_patients_by_cell_count(patient_ids, all_cell_metadata)
  
  if (length(suitable_patients) == 0) {
    return(NULL)
  }
  
  all_results <- list()
  successful_patients <- character()
  
  for (i in 1:length(suitable_patients)) {
    patient_id <- suitable_patients[i]
    folder_name <- patient_to_folder[patient_id]
    
    full_path <- file.path(dataset_config$infercnv_base_path, folder_name)
    if (!dir.exists(full_path)) {
      next
    }
    
    tryCatch({
      result <- process_patient(patient_id, folder_name, all_cell_metadata, dataset_config)
      
      if (!is.null(result)) {
        all_results[[patient_id]] <- result
        successful_patients <- c(successful_patients, patient_id)
      }
      
    }, error = function(e) {
      # Skip failed patients silently
    })
  }
  
  if (length(successful_patients) > 0) {
    regions_summary <- data.frame()
    genes_summary <- data.frame()
    
    for (patient_id in successful_patients) {
      result <- all_results[[patient_id]]
      
      for (data_type in c("regions", "genes")) {
        if (data_type %in% names(result)) {
          metrics <- result[[data_type]]$summary_metrics
          
          summary_row <- data.frame(
            patient_id = patient_id,
            dataset = dataset_name,
            data_type = data_type,
            pos_cells_matched = metrics$pos_cells_matched,
            neg_cells_matched = metrics$neg_cells_matched,
            directional_concordance = ifelse("directional_concordance" %in% names(metrics), metrics$directional_concordance, NA),
            pos_coverage_by_neg = ifelse("pos_coverage_by_neg" %in% names(metrics), metrics$pos_coverage_by_neg, 0),
            neg_coverage_by_pos = ifelse("neg_coverage_by_pos" %in% names(metrics), metrics$neg_coverage_by_pos, 0),
            deletion_concordance = ifelse("deletion_concordance" %in% names(metrics), metrics$deletion_concordance, NA),
            amplification_concordance = ifelse("amplification_concordance" %in% names(metrics), metrics$amplification_concordance, NA),
            total_overlapping_items = ifelse("total_overlapping_items" %in% names(metrics), metrics$total_overlapping_items, 0),
            pos_total_items = ifelse("pos_total_items" %in% names(metrics), metrics$pos_total_items, 0),
            neg_total_items = ifelse("neg_total_items" %in% names(metrics), metrics$neg_total_items, 0),
            pos_specific_items = ifelse("pos_specific_items" %in% names(metrics), metrics$pos_specific_items, 0),
            neg_specific_items = ifelse("neg_specific_items" %in% names(metrics), metrics$neg_specific_items, 0),
            concordant_items = ifelse("concordant_items" %in% names(metrics), metrics$concordant_items, 0),
            pos_total_cnvs = metrics$pos_total_cnvs,
            neg_total_cnvs = metrics$neg_total_cnvs,
            pos_deletions = metrics$pos_deletions,
            pos_amplifications = metrics$pos_amplifications,
            neg_deletions = metrics$neg_deletions,
            neg_amplifications = metrics$neg_amplifications,
            stringsAsFactors = FALSE
          )
          
          if (data_type == "regions") {
            summary_row$pos_mean_region_size <- ifelse("pos_mean_region_size" %in% names(metrics), metrics$pos_mean_region_size, 0)
            summary_row$neg_mean_region_size <- ifelse("neg_mean_region_size" %in% names(metrics), metrics$neg_mean_region_size, 0)
            summary_row$pos_total_genomic_size <- ifelse("pos_total_genomic_size" %in% names(metrics), metrics$pos_total_genomic_size, 0)
            summary_row$neg_total_genomic_size <- ifelse("neg_total_genomic_size" %in% names(metrics), metrics$neg_total_genomic_size, 0)
            
            regions_summary <- rbind(regions_summary, summary_row)
          } else {
            genes_summary <- rbind(genes_summary, summary_row)
          }
        }
      }
    }
    
    if (nrow(regions_summary) > 0) {
      regions_summary[is.na(regions_summary)] <- 0
      regions_summary[sapply(regions_summary, is.infinite)] <- 0
    }
    
    if (nrow(genes_summary) > 0) {
      genes_summary[is.na(genes_summary)] <- 0
      genes_summary[sapply(genes_summary, is.infinite)] <- 0
    }
    
    return(list(regions = regions_summary, genes = genes_summary))
  } else {
    return(NULL)
  }
}

# Run analysis for all datasets
all_results <- list()

for (dataset_name in names(datasets)) {
  dataset_config <- datasets[[dataset_name]]
  
  if (file.exists(dataset_config$cd45_pos_path) && 
      file.exists(dataset_config$cd45_neg_path) &&
      dir.exists(dataset_config$infercnv_base_path)) {
    
    result <- analyze_dataset(dataset_config)
    all_results[[dataset_name]] <- result
  }
}

# Combine results across datasets
combined_regions_summary <- data.frame()
combined_genes_summary <- data.frame()

for (dataset_name in names(all_results)) {
  if (!is.null(all_results[[dataset_name]])) {
    if (!is.null(all_results[[dataset_name]]$regions)) {
      combined_regions_summary <- rbind(combined_regions_summary, all_results[[dataset_name]]$regions)
    }
    if (!is.null(all_results[[dataset_name]]$genes)) {
      combined_genes_summary <- rbind(combined_genes_summary, all_results[[dataset_name]]$genes)
    }
  }
}

# Save combined results
if (nrow(combined_regions_summary) > 0) {
  write.csv(combined_regions_summary, 
            file.path(output_dir, "regions_summary.csv"), 
            row.names = FALSE)
}

if (nrow(combined_genes_summary) > 0) {
  write.csv(combined_genes_summary, 
            file.path(output_dir, "genes_summary.csv"), 
            row.names = FALSE)
}

# Create metric summaries
if (nrow(combined_regions_summary) > 0 || nrow(combined_genes_summary) > 0) {
  
  if (nrow(combined_regions_summary) > 0) {
    regions_metrics_summary <- combined_regions_summary %>%
      group_by(dataset, data_type) %>%
      summarise(
        n_patients = n(),
        mean_directional_concordance = mean(directional_concordance, na.rm = TRUE),
        median_directional_concordance = median(directional_concordance, na.rm = TRUE),
        mean_pos_coverage_by_neg = mean(pos_coverage_by_neg, na.rm = TRUE),
        median_pos_coverage_by_neg = median(pos_coverage_by_neg, na.rm = TRUE),
        mean_deletion_concordance = mean(deletion_concordance, na.rm = TRUE),
        mean_amplification_concordance = mean(amplification_concordance, na.rm = TRUE),
        mean_overlapping_items = mean(total_overlapping_items, na.rm = TRUE),
        mean_pos_specific_items = mean(pos_specific_items, na.rm = TRUE),
        mean_neg_specific_items = mean(neg_specific_items, na.rm = TRUE),
        high_concordance_patients = sum(directional_concordance >= 0.7, na.rm = TRUE),
        medium_concordance_patients = sum(directional_concordance >= 0.5 & directional_concordance < 0.7, na.rm = TRUE),
        low_concordance_patients = sum(directional_concordance < 0.5, na.rm = TRUE),
        high_pos_coverage_patients = sum(pos_coverage_by_neg >= 0.7, na.rm = TRUE),
        medium_pos_coverage_patients = sum(pos_coverage_by_neg >= 0.5 & pos_coverage_by_neg < 0.7, na.rm = TRUE),
        low_pos_coverage_patients = sum(pos_coverage_by_neg < 0.5, na.rm = TRUE),
        mean_pos_genomic_size = mean(pos_total_genomic_size, na.rm = TRUE),
        mean_neg_genomic_size = mean(neg_total_genomic_size, na.rm = TRUE),
        mean_pos_region_size = mean(pos_mean_region_size, na.rm = TRUE),
        mean_neg_region_size = mean(neg_mean_region_size, na.rm = TRUE),
        .groups = 'drop'
      )
    
    write.csv(regions_metrics_summary, 
              file.path(output_dir, "regions_metrics_summary.csv"), 
              row.names = FALSE)
  }
  
  if (nrow(combined_genes_summary) > 0) {
    genes_metrics_summary <- combined_genes_summary %>%
      group_by(dataset, data_type) %>%
      summarise(
        n_patients = n(),
        mean_directional_concordance = mean(directional_concordance, na.rm = TRUE),
        median_directional_concordance = median(directional_concordance, na.rm = TRUE),
        mean_pos_coverage_by_neg = mean(pos_coverage_by_neg, na.rm = TRUE),
        median_pos_coverage_by_neg = median(pos_coverage_by_neg, na.rm = TRUE),
        mean_deletion_concordance = mean(deletion_concordance, na.rm = TRUE),
        mean_amplification_concordance = mean(amplification_concordance, na.rm = TRUE),
        mean_overlapping_items = mean(total_overlapping_items, na.rm = TRUE),
        mean_pos_specific_items = mean(pos_specific_items, na.rm = TRUE),
        mean_neg_specific_items = mean(neg_specific_items, na.rm = TRUE),
        high_concordance_patients = sum(directional_concordance >= 0.7, na.rm = TRUE),
        medium_concordance_patients = sum(directional_concordance >= 0.5 & directional_concordance < 0.7, na.rm = TRUE),
        low_concordance_patients = sum(directional_concordance < 0.5, na.rm = TRUE),
        high_pos_coverage_patients = sum(pos_coverage_by_neg >= 0.7, na.rm = TRUE),
        medium_pos_coverage_patients = sum(pos_coverage_by_neg >= 0.5 & pos_coverage_by_neg < 0.7, na.rm = TRUE),
        low_pos_coverage_patients = sum(pos_coverage_by_neg < 0.5, na.rm = TRUE),
        .groups = 'drop'
      )
    
    write.csv(genes_metrics_summary, 
              file.path(output_dir, "genes_metrics_summary.csv"), 
              row.names = FALSE)
  }
  
  # Create combined summary with common columns
  if (nrow(combined_regions_summary) > 0 && nrow(combined_genes_summary) > 0) {
    regions_cols <- colnames(regions_metrics_summary)
    genes_cols <- colnames(genes_metrics_summary)
    common_cols <- intersect(regions_cols, genes_cols)
    
    if (length(common_cols) > 0) {
      combined_metrics_summary <- rbind(
        regions_metrics_summary[, common_cols, drop = FALSE],
        genes_metrics_summary[, common_cols, drop = FALSE]
      )
      
      write.csv(combined_metrics_summary, 
                file.path(output_dir, "combined_metrics_summary.csv"), 
                row.names = FALSE)
    }
  }
}

# Create visualizations
if (nrow(combined_regions_summary) > 2 || nrow(combined_genes_summary) > 2) {
  tiff(file.path(output_dir, "concordance_analysis_plots.tiff"), width = 16, height = 12, units = "in", res = 300)
  
  par(mfrow = c(3, 2), mar = c(4, 4, 3, 2))
  
  # Region-level directional concordance distribution
  if (nrow(combined_regions_summary) > 2) {
    regions_concordance <- combined_regions_summary$directional_concordance[!is.na(combined_regions_summary$directional_concordance)]
    if (length(regions_concordance) > 0) {
      hist(regions_concordance,
           main = paste("Region-level Directional Concordance (n =", length(regions_concordance), ")"),
           xlab = "Directional Concordance", 
           col = "lightblue", breaks = 15,
           xlim = c(0, 1))
      abline(v = mean(regions_concordance), col = "red", lwd = 2, lty = 2)
      abline(v = 0.7, col = "green", lwd = 2, lty = 3)
      abline(v = 0.5, col = "orange", lwd = 2, lty = 3)
    }
  }
  
  # Gene-level directional concordance distribution
  if (nrow(combined_genes_summary) > 2) {
    genes_concordance <- combined_genes_summary$directional_concordance[!is.na(combined_genes_summary$directional_concordance)]
    if (length(genes_concordance) > 0) {
      hist(genes_concordance,
           main = paste("Gene-level Directional Concordance (n =", length(genes_concordance), ")"),
           xlab = "Directional Concordance", 
           col = "lightgreen", breaks = 15,
           xlim = c(0, 1))
      abline(v = mean(genes_concordance), col = "red", lwd = 2, lty = 2)
      abline(v = 0.7, col = "green", lwd = 2, lty = 3)
      abline(v = 0.5, col = "orange", lwd = 2, lty = 3)
    }
  }
  
  # Region-level CD45+ coverage distribution
  if (nrow(combined_regions_summary) > 2) {
    regions_pos_coverage <- combined_regions_summary$pos_coverage_by_neg[!is.na(combined_regions_summary$pos_coverage_by_neg)]
    if (length(regions_pos_coverage) > 0) {
      hist(regions_pos_coverage,
           main = paste("Region-level CD45+ Coverage (n =", length(regions_pos_coverage), ")"),
           xlab = "CD45+ Coverage by CD45-", 
           col = "lightpink", breaks = 15,
           xlim = c(0, 1))
      abline(v = mean(regions_pos_coverage), col = "red", lwd = 2, lty = 2)
      abline(v = 0.7, col = "green", lwd = 2, lty = 3)
      abline(v = 0.5, col = "orange", lwd = 2, lty = 3)
    }
  }
  
  # Gene-level CD45+ coverage distribution
  if (nrow(combined_genes_summary) > 2) {
    genes_pos_coverage <- combined_genes_summary$pos_coverage_by_neg[!is.na(combined_genes_summary$pos_coverage_by_neg)]
    if (length(genes_pos_coverage) > 0) {
      hist(genes_pos_coverage,
           main = paste("Gene-level CD45+ Coverage (n =", length(genes_pos_coverage), ")"),
           xlab = "CD45+ Coverage by CD45-", 
           col = "lightsteelblue", breaks = 15,
           xlim = c(0, 1))
      abline(v = mean(genes_pos_coverage), col = "red", lwd = 2, lty = 2)
      abline(v = 0.7, col = "green", lwd = 2, lty = 3)
      abline(v = 0.5, col = "orange", lwd = 2, lty = 3)
    }
  }
  
  # CNV type-specific concordance comparison
  if (nrow(combined_regions_summary) > 2) {
    del_concordance <- combined_regions_summary$deletion_concordance[!is.na(combined_regions_summary$deletion_concordance)]
    amp_concordance <- combined_regions_summary$amplification_concordance[!is.na(combined_regions_summary$amplification_concordance)]
    
    if (length(del_concordance) > 0 && length(amp_concordance) > 0) {
      boxplot(list("Deletions" = del_concordance, "Amplifications" = amp_concordance),
              main = "CNV Type-Specific Concordance (Regions)",
              ylab = "Directional Concordance",
              col = c("red", "blue"))
    }
  }
  
  dev.off()
}







#################################################################################
# Part 3: Patient Coverage Heatmap of CD45-pos CNV Regions Across Chromosomes   #
#################################################################################

# Load required libraries
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

# Define output directory
output_dir <- "/R/R_All_Breast_Tumors/inferCNV_Analysis/CD45_Coverage_Chromosome_Heatmap/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Dataset configuration
datasets <- list(
  visvader = list(
    name = "visvader",
    cd45_pos_path = "/R/R_Visvader/Visvader_RDS/RDS_Annotated/CD45.pos.Visvader_panCK_cnv.rds",
    cd45_neg_path = "/R/R_Visvader/Visvader_RDS/RDS_Annotated/CD45.neg.Visvader_panCK_cnv.rds",
    infercnv_base_path = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output"
  ),
  swarbrick = list(
    name = "swarbrick",
    cd45_pos_path = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Annotated/CD45.pos.Swarbrick_panCK_cnv.rds",
    cd45_neg_path = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Annotated/CD45.neg.Swarbrick_panCK_cnv.rds",
    infercnv_base_path = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output"
  ),
  pannthr = list(
    name = "pannthr",
    cd45_pos_path = "/R/R_PANNTHR/PANNTHR_RDS/RDS_Annotated/CD45.pos.PANNTHR_panCK_cnv.rds",
    cd45_neg_path = "/R/R_PANNTHR/PANNTHR_RDS/RDS_Annotated/CD45.neg.PANNTHR_panCK_cnv.rds",
    infercnv_base_path = "/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Output"
  )
)

# Analysis parameters
MIN_CD45_POS <- 10
REGIONS_FILE <- "HMM_CNV_predictions.HMMi6.leiden.hmm_mode-subclusters.Pnorm_0.5.pred_cnv_regions.dat"
GROUPING_FILE <- "infercnv.17_HMM_predHMMi6.leiden.hmm_mode-subclusters.observation_groupings.txt"

# Convert original patient IDs to desired format
convert_patient_id_to_desired <- function(orig_id, dataset_name) {
  if (dataset_name == "visvader") {
    desired_id <- gsub("_Total$", "", orig_id)
    desired_id <- gsub("^Visvader_", "Pal_", desired_id)
    if (orig_id == "Visvader_0114_ER_Total") {
      desired_id <- "Pal_0114_ER_ER"
    }
    return(desired_id)
    
  } else if (dataset_name == "swarbrick") {
    swarbrick_mapping <- list(
      "CID3941" = "Wu_CID3941_ER", "CID3948" = "Wu_CID3948_ER", 
      "CID4040" = "Wu_CID4040_ER", "CID4067" = "Wu_CID4067_ER",
      "CID4290A" = "Wu_CID4290A_ER", "CID4398" = "Wu_CID4398_ER",
      "CID4461" = "Wu_CID4461_ER", "CID4463" = "Wu_CID4463_ER",
      "CID4471" = "Wu_CID4471_ER", "CID4530N" = "Wu_CID4530N_ER",
      "CID4535" = "Wu_CID4535_ER", "CID3586" = "Wu_CID3586_HER2",
      "CID3838" = "Wu_CID3838_HER2", "CID3921" = "Wu_CID3921_HER2", 
      "CID4066" = "Wu_CID4066_HER2", "CID45171" = "Wu_CID45171_HER2",
      "CID3946" = "Wu_CID3946_TNBC", "CID3963" = "Wu_CID3963_TNBC",
      "CID44041" = "Wu_CID44041_TNBC", "CID4465" = "Wu_CID4465_TNBC",
      "CID4495" = "Wu_CID4495_TNBC", "CID44971" = "Wu_CID44971_TNBC",
      "CID44991" = "Wu_CID44991_TNBC", "CID4513" = "Wu_CID4513_TNBC",
      "CID4515" = "Wu_CID4515_TNBC", "CID4523" = "Wu_CID4523_TNBC"
    )
    
    if (orig_id %in% names(swarbrick_mapping)) {
      return(swarbrick_mapping[[orig_id]])
    } else {
      return(paste0("Wu_", orig_id, "_ER"))
    }
    
  } else if (dataset_name == "pannthr") {
    return(paste0(orig_id, "_ER"))
  } else {
    return(orig_id)
  }
}

# Create cell metadata
create_cell_metadata <- function(seurat_obj, cell_type, dataset_name) {
  original_barcodes <- rownames(seurat_obj@meta.data)
  patient_ids <- seurat_obj@meta.data$orig.ident
  
  metadata <- data.frame(
    original_barcode = original_barcodes,
    patient_id = patient_ids,
    cell_type = cell_type,
    dataset = dataset_name,
    stringsAsFactors = FALSE
  )
  
  return(metadata)
}

# Map HMM states to CNV calls
map_hmm_states_to_cnv <- function(state_vector) {
  result <- character(length(state_vector))
  numeric_states <- as.numeric(state_vector)
  
  result[numeric_states <= 2] <- "deletion"
  result[numeric_states == 3] <- "neutral"
  result[numeric_states >= 4] <- "amplification"
  result[is.na(numeric_states)] <- "neutral"
  
  return(result)
}

# Normalize group names for matching
normalize_group_names <- function(group_names) {
  normalized <- character(length(group_names))
  
  for (i in seq_along(group_names)) {
    group_name <- group_names[i]
    
    if (grepl("Tumor\\.Tumor_s", group_name)) {
      parts <- strsplit(group_name, "\\.")[[1]]
      tumor_part <- parts[grepl("^Tumor_s", parts)]
      if (length(tumor_part) > 0) {
        normalized[i] <- tumor_part[1]
      } else {
        normalized[i] <- group_name
      }
    } else {
      normalized[i] <- group_name
    }
  }
  
  return(normalized)
}

# Create patient folder mapping
create_patient_folder_mapping <- function(patient_ids, dataset_name, infercnv_base_path) {
  if (dataset_name == "visvader") {
    sapply(patient_ids, function(p) {
      base_id <- sub("_Total$", "", p)
      paste0(base_id, "_panCK")
    })
  } else if (dataset_name == "swarbrick") {
    result <- character(length(patient_ids))
    names(result) <- patient_ids
    
    all_dirs <- list.dirs(infercnv_base_path, full.names = FALSE, recursive = FALSE)
    
    for (i in seq_along(patient_ids)) {
      patient_id <- patient_ids[i]
      base_id <- sub("_Total$", "", patient_id)
      
      matching_dirs <- grep(paste0(base_id, ".*_panCK$"), all_dirs, value = TRUE)
      
      if (length(matching_dirs) >= 1) {
        result[i] <- matching_dirs[1]
      } else {
        result[i] <- paste0("NOTFOUND_", base_id, "_panCK")
      }
    }
    
    return(result)
  } else if (dataset_name == "pannthr") {
    sapply(patient_ids, function(p) {
      base_id <- sub("_Total$", "", p)
      paste0(base_id, "_panCK")
    })
  }
}

# Extract cell group assignments
extract_cell_group_assignments <- function(patient_dir, patient_id) {
  grouping_file_path <- file.path(patient_dir, GROUPING_FILE)
  
  if (!file.exists(grouping_file_path)) {
    return(NULL)
  }
  
  tryCatch({
    raw_lines <- readLines(grouping_file_path)
    
    if (length(raw_lines) < 2) {
      return(NULL)
    }
    
    data_lines <- raw_lines[-1]
    parsed_data <- list()
    
    for (i in 1:length(data_lines)) {
      line <- data_lines[i]
      
      if (grepl('^"', line)) {
        parts <- strsplit(line, " ")[[1]]
        parts <- gsub('"', '', parts)
        parts <- parts[parts != ""]
        
        if (length(parts) >= 2) {
          barcode <- parts[1]
          group <- parts[2]
          
          if (barcode != "" && group != "") {
            parsed_data[[length(parsed_data) + 1]] <- data.frame(
              cell_barcode = barcode,
              cell_group_name = group,
              stringsAsFactors = FALSE
            )
          }
        }
      }
    }
    
    if (length(parsed_data) == 0) {
      return(NULL)
    }
    
    groupings <- do.call(rbind, parsed_data)
    groupings <- groupings[!is.na(groupings$cell_barcode) & groupings$cell_barcode != "", ]
    groupings <- groupings[!is.na(groupings$cell_group_name) & groupings$cell_group_name != "", ]
    
    return(groupings)
    
  }, error = function(e) {
    return(NULL)
  })
}

# Load CNV regions data
load_cnv_regions_data <- function(patient_dir) {
  file_path <- file.path(patient_dir, REGIONS_FILE)
  
  if (!file.exists(file_path)) {
    return(NULL)
  }
  
  tryCatch({
    data <- read.table(file_path, header = TRUE, sep = "\t", 
                       stringsAsFactors = FALSE, quote = "", check.names = FALSE)
    return(data)
  }, error = function(e) {
    return(NULL)
  })
}

# Extract CNV data for specific cells
extract_cell_cnv_data <- function(cnv_data, cell_groupings, cell_barcodes, cell_type, patient_id) {
  if (is.null(cnv_data) || is.null(cell_groupings)) {
    return(NULL)
  }
  
  matched_cell_groups <- cell_groupings$cell_group_name[cell_groupings$cell_barcode %in% cell_barcodes]
  
  if (length(matched_cell_groups) == 0) {
    return(NULL)
  }
  
  unique_matched_groups <- unique(matched_cell_groups)
  cnv_data$normalized_group_name <- normalize_group_names(cnv_data$cell_group_name)
  
  cell_cnv_data <- cnv_data %>%
    filter(normalized_group_name %in% unique_matched_groups) %>%
    mutate(
      cnv_call = map_hmm_states_to_cnv(state),
      cell_type = cell_type,
      patient_id = patient_id
    ) %>%
    filter(cnv_call != "neutral") %>%
    mutate(
      region_size = end - start + 1,
      region_id = paste0(chr, ":", start, "-", end)
    )
  
  return(cell_cnv_data)
}

# Calculate chromosome-specific CD45+ coverage
calculate_chromosome_cd45_coverage <- function(pos_cnv_data, neg_cnv_data, patient_id, desired_patient_id) {
  if (is.null(pos_cnv_data) || is.null(neg_cnv_data) || 
      nrow(pos_cnv_data) == 0 || nrow(neg_cnv_data) == 0) {
    return(NULL)
  }
  
  pos_cnv_data$chr <- as.character(pos_cnv_data$chr)
  neg_cnv_data$chr <- as.character(neg_cnv_data$chr)
  
  valid_chr_pattern <- "^chr([1-9]|1[0-9]|2[0-2]|X|Y)$"
  pos_cnv_data <- pos_cnv_data[grepl(valid_chr_pattern, pos_cnv_data$chr), ]
  neg_cnv_data <- neg_cnv_data[grepl(valid_chr_pattern, neg_cnv_data$chr), ]
  
  if (nrow(pos_cnv_data) == 0 || nrow(neg_cnv_data) == 0) {
    return(NULL)
  }
  
  pos_cnv_data$region_id <- make.unique(pos_cnv_data$region_id)
  neg_cnv_data$region_id <- make.unique(neg_cnv_data$region_id)
  
  pos_chrs <- unique(pos_cnv_data$chr)
  
  chr_coverage <- data.frame()
  
  for (chr in pos_chrs) {
    pos_chr_data <- pos_cnv_data[pos_cnv_data$chr == chr, ]
    neg_chr_data <- neg_cnv_data[neg_cnv_data$chr == chr, ]
    
    pos_regions <- unique(pos_chr_data$region_id)
    neg_regions <- if (nrow(neg_chr_data) > 0) unique(neg_chr_data$region_id) else character(0)
    
    overlap_regions <- intersect(pos_regions, neg_regions)
    coverage <- if (length(pos_regions) > 0) length(overlap_regions) / length(pos_regions) else 0
    
    if (is.na(coverage) || is.infinite(coverage) || is.nan(coverage)) {
      coverage <- 0
    }
    
    chr_coverage <- rbind(chr_coverage, data.frame(
      patient_id = patient_id,
      desired_patient_id = desired_patient_id,
      chromosome = chr,
      cd45_pos_coverage = coverage,
      pos_total_regions = length(pos_regions),
      neg_total_regions = length(neg_regions),
      overlap_regions = length(overlap_regions),
      stringsAsFactors = FALSE
    ))
  }
  
  return(chr_coverage)
}

# Process a single patient
process_patient_coverage <- function(patient_id, folder_name, patient_metadata, dataset_config) {
  patient_dir <- file.path(dataset_config$infercnv_base_path, folder_name)
  
  if (!dir.exists(patient_dir)) {
    return(NULL)
  }
  
  patient_cells <- patient_metadata %>% filter(patient_id == !!patient_id)
  pos_cells <- patient_cells$original_barcode[patient_cells$cell_type == "CD45-pos"]
  neg_cells <- patient_cells$original_barcode[patient_cells$cell_type == "CD45-neg"]
  
  if (length(pos_cells) < MIN_CD45_POS) {
    return(NULL)
  }
  
  cell_groupings <- extract_cell_group_assignments(patient_dir, patient_id)
  if (is.null(cell_groupings)) {
    return(NULL)
  }
  
  infercnv_barcodes <- cell_groupings$cell_barcode
  pos_matches <- intersect(pos_cells, infercnv_barcodes)
  neg_matches <- intersect(neg_cells, infercnv_barcodes)
  
  if (length(pos_matches) < MIN_CD45_POS) {
    return(NULL)
  }
  
  cnv_data <- load_cnv_regions_data(patient_dir)
  if (is.null(cnv_data)) {
    return(NULL)
  }
  
  pos_cnv_data <- extract_cell_cnv_data(cnv_data, cell_groupings, pos_matches, "CD45_pos", patient_id)
  neg_cnv_data <- extract_cell_cnv_data(cnv_data, cell_groupings, neg_matches, "CD45_neg", patient_id)
  
  if (is.null(pos_cnv_data) || is.null(neg_cnv_data) || 
      nrow(pos_cnv_data) == 0 || nrow(neg_cnv_data) == 0) {
    return(NULL)
  }
  
  desired_patient_id <- convert_patient_id_to_desired(patient_id, dataset_config$name)
  
  chr_coverage <- calculate_chromosome_cd45_coverage(pos_cnv_data, neg_cnv_data, patient_id, desired_patient_id)
  
  if (!is.null(chr_coverage)) {
    chr_coverage$dataset <- dataset_config$name
    chr_coverage$pos_cells_matched <- length(pos_matches)
    chr_coverage$neg_cells_matched <- length(neg_matches)
  }
  
  return(chr_coverage)
}

# Analyze dataset for chromosome coverage
analyze_dataset_chromosome_coverage <- function(dataset_config) {
  cd45_pos <- readRDS(dataset_config$cd45_pos_path)
  cd45_neg <- readRDS(dataset_config$cd45_neg_path)
  
  cell_metadata_pos <- create_cell_metadata(cd45_pos, "CD45-pos", dataset_config$name)
  cell_metadata_neg <- create_cell_metadata(cd45_neg, "CD45-neg", dataset_config$name)
  all_cell_metadata <- rbind(cell_metadata_pos, cell_metadata_neg)
  
  patient_ids_pos <- unique(cd45_pos@meta.data$orig.ident)
  patient_ids_neg <- unique(cd45_neg@meta.data$orig.ident)
  patient_ids <- intersect(patient_ids_pos, patient_ids_neg)
  
  rm(cd45_pos, cd45_neg)
  gc()
  
  patient_to_folder <- create_patient_folder_mapping(patient_ids, dataset_config$name, dataset_config$infercnv_base_path)
  
  suitable_patients <- character()
  for (patient in patient_ids) {
    pos_count <- sum(all_cell_metadata$patient_id == patient & all_cell_metadata$cell_type == "CD45-pos")
    if (pos_count >= MIN_CD45_POS) {
      suitable_patients <- c(suitable_patients, patient)
    }
  }
  
  if (length(suitable_patients) == 0) {
    return(NULL)
  }
  
  all_coverage_results <- data.frame()
  
  for (patient_id in suitable_patients) {
    folder_name <- patient_to_folder[patient_id]
    
    full_path <- file.path(dataset_config$infercnv_base_path, folder_name)
    if (!dir.exists(full_path)) {
      next
    }
    
    tryCatch({
      result <- process_patient_coverage(patient_id, folder_name, all_cell_metadata, dataset_config)
      
      if (!is.null(result)) {
        all_coverage_results <- rbind(all_coverage_results, result)
      }
      
    }, error = function(e) {
      next
    })
  }
  
  return(all_coverage_results)
}

# Create chromosome coverage heatmap
create_chromosome_coverage_heatmap <- function(all_coverage_data) {
  if (nrow(all_coverage_data) == 0) {
    return(NULL)
  }
  
  all_coverage_data <- all_coverage_data[
    is.finite(all_coverage_data$cd45_pos_coverage) & 
      !is.na(all_coverage_data$cd45_pos_coverage), ]
  
  if (nrow(all_coverage_data) == 0) {
    return(NULL)
  }
  
  all_coverage_data$cd45_pos_coverage <- pmax(0, pmin(1, all_coverage_data$cd45_pos_coverage))
  
  chromosome_order <- c(paste0("chr", 1:22), "chrX", "chrY")
  available_chromosomes <- intersect(chromosome_order, unique(all_coverage_data$chromosome))
  all_coverage_data$chromosome <- factor(all_coverage_data$chromosome, levels = available_chromosomes)
  
  # Order patients alphabetically within dataset groups (Pal, Wu, PANNTHR)
  dataset_order <- c("visvader", "swarbrick", "pannthr")
  all_coverage_data$dataset <- factor(all_coverage_data$dataset, levels = dataset_order)
  
  # Extract suffix for subtype grouping
  all_coverage_data$subtype <- case_when(
    grepl("_ER$|_ER_ER$", all_coverage_data$desired_patient_id) ~ "ER",
    grepl("_PR$", all_coverage_data$desired_patient_id) ~ "PR", 
    grepl("_HER2$", all_coverage_data$desired_patient_id) ~ "HER2",
    grepl("_TNBC$", all_coverage_data$desired_patient_id) ~ "TNBC",
    TRUE ~ "ER"
  )
  
  subtype_order <- c("ER", "PR", "HER2", "TNBC")
  all_coverage_data$subtype <- factor(all_coverage_data$subtype, levels = subtype_order)
  
  patient_order_data <- all_coverage_data %>%
    select(desired_patient_id, dataset, subtype) %>%
    distinct() %>%
    arrange(desc(dataset), desc(subtype), desc(desired_patient_id))
  
  all_coverage_data$desired_patient_id <- factor(all_coverage_data$desired_patient_id, 
                                                 levels = patient_order_data$desired_patient_id)
  
  p1 <- ggplot(all_coverage_data, aes(x = chromosome, y = desired_patient_id, fill = cd45_pos_coverage)) +
    geom_tile(color = "white", size = 0.3) +
    scale_fill_gradient(
      low = "blue", 
      high = "red",
      na.value = "grey90",
      name = "CD45+\nCoverage",
      labels = percent_format(accuracy = 1),
      guide = guide_colorbar(title.position = "top", title.hjust = 0.5)
    ) +
    labs(
      title = "CD45+ Coverage by Chromosome",
      subtitle = "Fraction of CD45+ CNV regions also found in CD45- cells",
      x = "Chromosome",
      y = "Patient"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 8),
      axis.title = element_text(size = 12),
      plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9),
      panel.grid = element_blank(),
      panel.border = element_blank()
    )
  
  ggsave(file.path(output_dir, "CD45_Coverage_Chromosome_Heatmap.tiff"), 
         plot = p1, width = 16, height = 12, dpi = 300)
  
  return(p1)
}

# Main analysis function
main_chromosome_coverage_analysis <- function() {
  all_coverage_data <- data.frame()
  
  for (dataset_name in names(datasets)) {
    dataset_config <- datasets[[dataset_name]]
    
    if (file.exists(dataset_config$cd45_pos_path) && 
        file.exists(dataset_config$cd45_neg_path) &&
        dir.exists(dataset_config$infercnv_base_path)) {
      
      result <- analyze_dataset_chromosome_coverage(dataset_config)
      
      if (!is.null(result)) {
        all_coverage_data <- rbind(all_coverage_data, result)
      }
    }
  }
  
  if (nrow(all_coverage_data) == 0) {
    return(NULL)
  }
  
  write.csv(all_coverage_data, file.path(output_dir, "complete_chromosome_coverage_data.csv"), row.names = FALSE)
  
  heatmap_plot <- create_chromosome_coverage_heatmap(all_coverage_data)
  
  return(all_coverage_data)
}

# Execute analysis
coverage_results <- main_chromosome_coverage_analysis()







#################################################################################
# Part 4: Within vs Between Patient CD45-pos CNV Region Correlations            #
#################################################################################

# Load required libraries
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)

# Define output directory
output_dir <- "/R/R_All_Breast_Tumors/inferCNV_Analysis/Patient_CNV_Region_Correlations/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Dataset configuration
datasets <- list(
  visvader = list(
    name = "visvader",
    cd45_pos_path = "/R/R_Visvader/Visvader_RDS/RDS_Annotated/CD45.pos.Visvader_panCK_cnv.rds",
    cd45_neg_path = "/R/R_Visvader/Visvader_RDS/RDS_Annotated/CD45.neg.Visvader_panCK_cnv.rds",
    infercnv_base_path = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/"
  ),
  swarbrick = list(
    name = "swarbrick",
    cd45_pos_path = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Annotated/CD45.pos.Swarbrick_panCK_cnv.rds",
    cd45_neg_path = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Annotated/CD45.neg.Swarbrick_panCK_cnv.rds",
    infercnv_base_path = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/"
  ),
  pannthr = list(
    name = "pannthr",
    cd45_pos_path = "/R/R_PANNTHR/PANNTHR_RDS/RDS_Annotated/CD45.pos.PANNTHR_panCK_cnv.rds",
    cd45_neg_path = "/R/R_PANNTHR/PANNTHR_RDS/RDS_Annotated/CD45.neg.PANNTHR_panCK_cnv.rds",
    infercnv_base_path = "/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Output/"
  )
)

# Analysis parameters
MIN_CD45_POS <- 10
REGIONS_FILE <- "HMM_CNV_predictions.HMMi6.leiden.hmm_mode-subclusters.Pnorm_0.5.pred_cnv_regions.dat"
GROUPING_FILE <- "infercnv.17_HMM_predHMMi6.leiden.hmm_mode-subclusters.observation_groupings.txt"

# Core analysis functions
create_cell_metadata <- function(seurat_obj, cell_type, dataset_name) {
  original_barcodes <- rownames(seurat_obj@meta.data)
  patient_ids <- seurat_obj@meta.data$orig.ident
  
  metadata <- data.frame(
    original_barcode = original_barcodes,
    patient_id = patient_ids,
    cell_type = cell_type,
    dataset = dataset_name,
    stringsAsFactors = FALSE
  )
  
  return(metadata)
}

map_hmm_states_to_cnv <- function(state_vector) {
  result <- character(length(state_vector))
  numeric_states <- as.numeric(state_vector)
  
  result[numeric_states <= 2] <- "deletion"
  result[numeric_states == 3] <- "neutral"
  result[numeric_states >= 4] <- "amplification"
  result[is.na(numeric_states)] <- "neutral"
  
  return(result)
}

normalize_group_names <- function(group_names) {
  normalized <- character(length(group_names))
  
  for (i in seq_along(group_names)) {
    group_name <- group_names[i]
    
    if (grepl("Tumor\\.Tumor_s", group_name)) {
      parts <- strsplit(group_name, "\\.")[[1]]
      tumor_part <- parts[grepl("^Tumor_s", parts)]
      if (length(tumor_part) > 0) {
        normalized[i] <- tumor_part[1]
      } else {
        normalized[i] <- group_name
      }
    } else {
      normalized[i] <- group_name
    }
  }
  
  return(normalized)
}

create_patient_folder_mapping <- function(patient_ids, dataset_name, infercnv_base_path) {
  if (dataset_name == "visvader") {
    sapply(patient_ids, function(p) {
      base_id <- sub("_Total$", "", p)
      paste0(base_id, "_panCK")
    })
  } else if (dataset_name == "swarbrick") {
    result <- character(length(patient_ids))
    names(result) <- patient_ids
    
    all_dirs <- list.dirs(infercnv_base_path, full.names = FALSE, recursive = FALSE)
    
    for (i in seq_along(patient_ids)) {
      patient_id <- patient_ids[i]
      base_id <- sub("_Total$", "", patient_id)
      
      matching_dirs <- grep(paste0(base_id, ".*_panCK$"), all_dirs, value = TRUE)
      
      if (length(matching_dirs) >= 1) {
        result[i] <- matching_dirs[1]
      } else {
        result[i] <- paste0("NOTFOUND_", base_id, "_panCK")
      }
    }
    
    return(result)
  } else if (dataset_name == "pannthr") {
    sapply(patient_ids, function(p) {
      base_id <- sub("_Total$", "", p)
      paste0(base_id, "_panCK")
    })
  }
}

extract_cell_group_assignments <- function(patient_dir, patient_id) {
  grouping_file_path <- file.path(patient_dir, GROUPING_FILE)
  
  if (!file.exists(grouping_file_path)) {
    return(NULL)
  }
  
  tryCatch({
    raw_lines <- readLines(grouping_file_path)
    
    if (length(raw_lines) < 2) {
      return(NULL)
    }
    
    data_lines <- raw_lines[-1]
    parsed_data <- list()
    
    for (i in 1:length(data_lines)) {
      line <- data_lines[i]
      
      if (grepl('^"', line)) {
        parts <- strsplit(line, " ")[[1]]
        parts <- gsub('"', '', parts)
        parts <- parts[parts != ""]
        
        if (length(parts) >= 2) {
          barcode <- parts[1]
          group <- parts[2]
          
          if (barcode != "" && group != "") {
            parsed_data[[length(parsed_data) + 1]] <- data.frame(
              cell_barcode = barcode,
              cell_group_name = group,
              stringsAsFactors = FALSE
            )
          }
        }
      }
    }
    
    if (length(parsed_data) == 0) {
      return(NULL)
    }
    
    groupings <- do.call(rbind, parsed_data)
    groupings <- groupings[!is.na(groupings$cell_barcode) & groupings$cell_barcode != "", ]
    groupings <- groupings[!is.na(groupings$cell_group_name) & groupings$cell_group_name != "", ]
    
    return(groupings)
    
  }, error = function(e) {
    return(NULL)
  })
}

load_cnv_regions_data <- function(patient_dir) {
  file_path <- file.path(patient_dir, REGIONS_FILE)
  
  if (!file.exists(file_path)) {
    return(NULL)
  }
  
  tryCatch({
    data <- read.table(file_path, header = TRUE, sep = "\t", 
                       stringsAsFactors = FALSE, quote = "", check.names = FALSE)
    return(data)
  }, error = function(e) {
    return(NULL)
  })
}

extract_cell_cnv_data <- function(cnv_data, cell_groupings, cell_barcodes, cell_type, patient_id) {
  if (is.null(cnv_data) || is.null(cell_groupings)) {
    return(NULL)
  }
  
  matched_cell_groups <- cell_groupings$cell_group_name[cell_groupings$cell_barcode %in% cell_barcodes]
  
  if (length(matched_cell_groups) == 0) {
    return(NULL)
  }
  
  unique_matched_groups <- unique(matched_cell_groups)
  cnv_data$normalized_group_name <- normalize_group_names(cnv_data$cell_group_name)
  
  cell_cnv_data <- cnv_data %>%
    filter(normalized_group_name %in% unique_matched_groups) %>%
    mutate(
      cnv_call = map_hmm_states_to_cnv(state),
      cell_type = cell_type,
      patient_id = patient_id
    ) %>%
    filter(cnv_call != "neutral") %>%
    mutate(
      region_size = end - start + 1,
      region_id = paste0(chr, ":", start, "-", end)
    )
  
  return(cell_cnv_data)
}

create_chromosome_cnv_profile <- function(cnv_data, patient_id, cell_type) {
  if (is.null(cnv_data) || nrow(cnv_data) == 0) {
    return(NULL)
  }
  
  chr_profile <- cnv_data %>%
    group_by(chr) %>%
    summarise(
      total_cnv_regions = n(),
      deletion_regions = sum(cnv_call == "deletion"),
      amplification_regions = sum(cnv_call == "amplification"),
      mean_region_size = mean(region_size, na.rm = TRUE),
      total_genomic_size = sum(region_size, na.rm = TRUE),
      deletion_prop = sum(cnv_call == "deletion") / n(),
      amplification_prop = sum(cnv_call == "amplification") / n(),
      .groups = 'drop'
    ) %>%
    mutate(
      patient_id = patient_id,
      cell_type = cell_type
    )
  
  return(chr_profile)
}

calculate_patient_cnv_correlation <- function(pos_profile, neg_profile, patient_id) {
  if (is.null(pos_profile) || is.null(neg_profile)) {
    return(data.frame(
      patient_id = patient_id,
      correlation = NA,
      n_common_chromosomes = 0,
      correlation_method = "missing_data",
      stringsAsFactors = FALSE
    ))
  }
  
  common_chromosomes <- intersect(pos_profile$chr, neg_profile$chr)
  
  if (length(common_chromosomes) < 2) {
    return(data.frame(
      patient_id = patient_id,
      correlation = NA,
      n_common_chromosomes = length(common_chromosomes),
      correlation_method = "insufficient_chromosomes",
      stringsAsFactors = FALSE
    ))
  }
  
  pos_common <- pos_profile %>% filter(chr %in% common_chromosomes) %>% arrange(chr)
  neg_common <- neg_profile %>% filter(chr %in% common_chromosomes) %>% arrange(chr)
  
  correlation_results <- list()
  
  if (all(pos_common$chr == neg_common$chr)) {
    pos_values <- pos_common$total_cnv_regions
    neg_values <- neg_common$total_cnv_regions
    
    pos_var <- var(pos_values)
    neg_var <- var(neg_values)
    
    if (pos_var > 0 && neg_var > 0) {
      cor_val <- tryCatch({
        cor(pos_values, neg_values, use = "complete.obs")
      }, error = function(e) NA)
      
      if (!is.na(cor_val) && is.finite(cor_val)) {
        correlation_results[["total_regions"]] <- list(
          correlation = cor_val,
          method = "total_cnv_regions"
        )
      }
    } else if (pos_var == 0 && neg_var == 0) {
      correlation_results[["total_regions"]] <- list(
        correlation = 1.0,
        method = "perfect_agreement"
      )
    }
  }
  
  if (length(correlation_results) == 0) {
    pos_values <- pos_common$total_genomic_size
    neg_values <- neg_common$total_genomic_size
    
    pos_var <- var(pos_values)
    neg_var <- var(neg_values)
    
    if (pos_var > 0 && neg_var > 0) {
      cor_val <- tryCatch({
        cor(pos_values, neg_values, use = "complete.obs")
      }, error = function(e) NA)
      
      if (!is.na(cor_val) && is.finite(cor_val)) {
        correlation_results[["genomic_size"]] <- list(
          correlation = cor_val,
          method = "total_genomic_size"
        )
      }
    } else if (pos_var == 0 && neg_var == 0) {
      correlation_results[["genomic_size"]] <- list(
        correlation = 1.0,
        method = "perfect_agreement"
      )
    }
  }
  
  if (length(correlation_results) == 0) {
    pos_values <- pos_common$deletion_prop
    neg_values <- neg_common$deletion_prop
    
    pos_var <- var(pos_values)
    neg_var <- var(neg_values)
    
    if (pos_var > 0 && neg_var > 0) {
      cor_val <- tryCatch({
        cor(pos_values, neg_values, use = "complete.obs")
      }, error = function(e) NA)
      
      if (!is.na(cor_val) && is.finite(cor_val)) {
        correlation_results[["deletion_prop"]] <- list(
          correlation = cor_val,
          method = "deletion_proportion"
        )
      }
    } else if (pos_var == 0 && neg_var == 0) {
      correlation_results[["deletion_prop"]] <- list(
        correlation = 1.0,
        method = "perfect_agreement"
      )
    }
  }
  
  if (length(correlation_results) > 0) {
    best_result <- correlation_results[[1]]
    final_correlation <- best_result$correlation
    final_method <- best_result$method
  } else {
    final_correlation <- NA
    final_method <- "calculation_failed"
  }
  
  return(data.frame(
    patient_id = patient_id,
    correlation = final_correlation,
    n_common_chromosomes = length(common_chromosomes),
    correlation_method = final_method,
    stringsAsFactors = FALSE
  ))
}

process_patient_correlation <- function(patient_id, folder_name, patient_metadata, dataset_config) {
  patient_dir <- file.path(dataset_config$infercnv_base_path, folder_name)
  
  if (!dir.exists(patient_dir)) {
    return(NULL)
  }
  
  patient_cells <- patient_metadata %>% filter(patient_id == !!patient_id)
  pos_cells <- patient_cells$original_barcode[patient_cells$cell_type == "CD45-pos"]
  neg_cells <- patient_cells$original_barcode[patient_cells$cell_type == "CD45-neg"]
  
  if (length(pos_cells) < MIN_CD45_POS) {
    return(NULL)
  }
  
  cell_groupings <- extract_cell_group_assignments(patient_dir, patient_id)
  if (is.null(cell_groupings)) {
    return(NULL)
  }
  
  infercnv_barcodes <- cell_groupings$cell_barcode
  pos_matches <- intersect(pos_cells, infercnv_barcodes)
  neg_matches <- intersect(neg_cells, infercnv_barcodes)
  
  if (length(pos_matches) < MIN_CD45_POS) {
    return(NULL)
  }
  
  cnv_data <- load_cnv_regions_data(patient_dir)
  if (is.null(cnv_data)) {
    return(NULL)
  }
  
  pos_cnv_data <- extract_cell_cnv_data(cnv_data, cell_groupings, pos_matches, "CD45_pos", patient_id)
  neg_cnv_data <- extract_cell_cnv_data(cnv_data, cell_groupings, neg_matches, "CD45_neg", patient_id)
  
  if (is.null(pos_cnv_data) || is.null(neg_cnv_data) || 
      nrow(pos_cnv_data) == 0 || nrow(neg_cnv_data) == 0) {
    return(NULL)
  }
  
  pos_chr_profile <- create_chromosome_cnv_profile(pos_cnv_data, patient_id, "CD45_pos")
  neg_chr_profile <- create_chromosome_cnv_profile(neg_cnv_data, patient_id, "CD45_neg")
  
  if (is.null(pos_chr_profile) || is.null(neg_chr_profile)) {
    return(NULL)
  }
  
  within_patient_cor <- calculate_patient_cnv_correlation(pos_chr_profile, neg_chr_profile, patient_id)
  
  return(list(
    patient_id = patient_id,
    dataset = dataset_config$name,
    pos_profile = pos_chr_profile,
    neg_profile = neg_chr_profile,
    within_patient_correlation = within_patient_cor,
    pos_cells_matched = length(pos_matches),
    neg_cells_matched = length(neg_matches)
  ))
}

create_patient_profiles <- function() {
  all_patient_profiles <- list()
  
  for (dataset_name in names(datasets)) {
    dataset_config <- datasets[[dataset_name]]
    
    if (!file.exists(dataset_config$cd45_pos_path) || 
        !file.exists(dataset_config$cd45_neg_path) ||
        !dir.exists(dataset_config$infercnv_base_path)) {
      next
    }
    
    cd45_pos <- readRDS(dataset_config$cd45_pos_path)
    cd45_neg <- readRDS(dataset_config$cd45_neg_path)
    
    cell_metadata_pos <- create_cell_metadata(cd45_pos, "CD45-pos", dataset_name)
    cell_metadata_neg <- create_cell_metadata(cd45_neg, "CD45-neg", dataset_name)
    all_cell_metadata <- rbind(cell_metadata_pos, cell_metadata_neg)
    
    patient_ids_pos <- unique(cd45_pos@meta.data$orig.ident)
    patient_ids_neg <- unique(cd45_neg@meta.data$orig.ident)
    patient_ids <- intersect(patient_ids_pos, patient_ids_neg)
    
    rm(cd45_pos, cd45_neg)
    gc()
    
    patient_to_folder <- create_patient_folder_mapping(patient_ids, dataset_name, dataset_config$infercnv_base_path)
    
    suitable_patients <- character()
    for (patient in patient_ids) {
      pos_count <- sum(all_cell_metadata$patient_id == patient & all_cell_metadata$cell_type == "CD45-pos")
      if (pos_count >= MIN_CD45_POS) {
        suitable_patients <- c(suitable_patients, patient)
      }
    }
    
    for (patient_id in suitable_patients) {
      folder_name <- patient_to_folder[patient_id]
      
      result <- tryCatch({
        process_patient_correlation(patient_id, folder_name, all_cell_metadata, dataset_config)
      }, error = function(e) NULL)
      
      if (!is.null(result)) {
        all_patient_profiles[[paste0(dataset_name, "_", patient_id)]] <- result
      }
    }
  }
  
  return(all_patient_profiles)
}

calculate_neoplastic_CD45_patient_correlations <- function(all_patient_profiles) {
  all_correlations <- data.frame()
  
  within_patient_correlations <- data.frame()
  
  for (profile_name in names(all_patient_profiles)) {
    profile <- all_patient_profiles[[profile_name]]
    
    if (!is.null(profile$within_patient_correlation)) {
      within_cor <- profile$within_patient_correlation
      within_cor$dataset <- profile$dataset
      within_cor$comparison_type <- "Within Patient\n(CD45+ vs CD45-)"
      
      within_patient_correlations <- rbind(within_patient_correlations, within_cor)
    }
  }
  
  between_pos_pos_correlations <- data.frame()
  
  patient_names <- names(all_patient_profiles)
  for (i in 1:length(patient_names)) {
    profile_i <- all_patient_profiles[[patient_names[i]]]
    
    correlations_for_patient_i <- c()
    
    for (j in 1:length(patient_names)) {
      if (i != j) {
        profile_j <- all_patient_profiles[[patient_names[j]]]
        
        cor_result <- calculate_patient_cnv_correlation(profile_i$pos_profile, profile_j$pos_profile, 
                                                        profile_i$patient_id)
        
        if (!is.null(cor_result) && !is.na(cor_result$correlation)) {
          correlations_for_patient_i <- c(correlations_for_patient_i, cor_result$correlation)
        }
      }
    }
    
    if (length(correlations_for_patient_i) > 0) {
      mean_cor <- mean(correlations_for_patient_i, na.rm = TRUE)
      
      between_pos_pos_correlations <- rbind(between_pos_pos_correlations, data.frame(
        patient_id = profile_i$patient_id,
        correlation = mean_cor,
        n_common_chromosomes = NA,
        correlation_method = "between_patient_mean",
        dataset = profile_i$dataset,
        comparison_type = "Between Patients\n(CD45+ vs CD45+)",
        stringsAsFactors = FALSE
      ))
    }
  }
  
  between_pos_neg_correlations <- data.frame()
  
  for (i in 1:length(patient_names)) {
    profile_i <- all_patient_profiles[[patient_names[i]]]
    
    correlations_for_patient_i <- c()
    
    for (j in 1:length(patient_names)) {
      if (i != j) {
        profile_j <- all_patient_profiles[[patient_names[j]]]
        
        cor_result <- calculate_patient_cnv_correlation(profile_i$pos_profile, profile_j$neg_profile, 
                                                        profile_i$patient_id)
        
        if (!is.null(cor_result) && !is.na(cor_result$correlation)) {
          correlations_for_patient_i <- c(correlations_for_patient_i, cor_result$correlation)
        }
      }
    }
    
    if (length(correlations_for_patient_i) > 0) {
      mean_cor <- mean(correlations_for_patient_i, na.rm = TRUE)
      
      between_pos_neg_correlations <- rbind(between_pos_neg_correlations, data.frame(
        patient_id = profile_i$patient_id,
        correlation = mean_cor,
        n_common_chromosomes = NA,
        correlation_method = "between_patient_mean",
        dataset = profile_i$dataset,
        comparison_type = "Between Patients\n(CD45+ vs CD45-)",
        stringsAsFactors = FALSE
      ))
    }
  }
  
  all_correlations <- rbind(within_patient_correlations, between_pos_pos_correlations, between_pos_neg_correlations)
  
  return(all_correlations)
}

create_neoplastic_CD45_patient_correlations_plot <- function(all_correlations) {
  if (nrow(all_correlations) == 0) {
    return(NULL)
  }
  
  write.csv(all_correlations, file.path(output_dir, "neoplastic_CD45_patient_correlations_data.csv"), row.names = FALSE)
  
  tiff(file.path(output_dir, "neoplastic_CD45_patient_correlations_plot.tiff"), width = 14, height = 6, units = "in", res = 300)
  
  layout(matrix(c(1, 1, 2), 1, 3, byrow = TRUE), widths = c(2, 1))
  
  par(mar = c(6, 5, 4, 2))
  
  valid_correlations <- all_correlations[!is.na(all_correlations$correlation), ]
  
  comparison_order <- c("Within Patient\n(CD45+ vs CD45-)", 
                        "Between Patients\n(CD45+ vs CD45+)", 
                        "Between Patients\n(CD45+ vs CD45-)")
  
  comparison_list <- split(valid_correlations$correlation, valid_correlations$comparison_type)
  comparison_list <- comparison_list[comparison_order[comparison_order %in% names(comparison_list)]]
  
  colors <- c("#3498db", "#e74c3c", "#f39c12")
  names(colors) <- comparison_order
  
  if (length(comparison_list) > 0) {
    bp <- boxplot(comparison_list,
                  main = "CNV Correlations: Within-Patient vs Between-Patient Comparisons",
                  ylab = "Correlation Coefficient",
                  xlab = "",
                  col = colors[names(comparison_list)],
                  border = "black",
                  outline = TRUE,
                  las = 1,
                  cex.main = 1.4,
                  cex.lab = 1.2,
                  cex.axis = 1.0)
    
    for (i in 1:length(comparison_list)) {
      if (length(comparison_list[[i]]) > 0) {
        points(jitter(rep(i, length(comparison_list[[i]])), factor = 0.3),
               comparison_list[[i]],
               pch = 16,
               col = alpha(colors[names(comparison_list)[i]], 0.6),
               cex = 0.8)
      }
    }
    
    for (i in 1:length(comparison_list)) {
      text(i, par("usr")[3] - 0.1 * diff(par("usr")[3:4]), 
           paste("n =", length(comparison_list[[i]])), 
           xpd = TRUE, cex = 0.9, adj = 0.5)
    }
    
    if (length(comparison_list) >= 2) {
      if ("Within Patient\n(CD45+ vs CD45-)" %in% names(comparison_list) &&
          "Between Patients\n(CD45+ vs CD45+)" %in% names(comparison_list)) {
        
        within_data <- comparison_list[["Within Patient\n(CD45+ vs CD45-)"]]
        between_pos_data <- comparison_list[["Between Patients\n(CD45+ vs CD45+)"]]
        
        if (length(within_data) > 0 && length(between_pos_data) > 0) {
          test_result <- wilcox.test(within_data, between_pos_data, alternative = "two.sided")
          
          y_max <- max(unlist(comparison_list), na.rm = TRUE)
          segments(1, y_max + 0.05, 2, y_max + 0.05, xpd = TRUE)
          segments(1, y_max + 0.03, 1, y_max + 0.05, xpd = TRUE)
          segments(2, y_max + 0.03, 2, y_max + 0.05, xpd = TRUE)
          
          p_text <- if (test_result$p.value < 0.001) "***" else 
            if (test_result$p.value < 0.01) "**" else 
              if (test_result$p.value < 0.05) "*" else "ns"
          
          text(1.5, y_max + 0.08, p_text, xpd = TRUE, cex = 1.2, adj = 0.5)
          text(1.5, y_max + 0.12, paste("p =", format.pval(test_result$p.value, digits = 3)), 
               xpd = TRUE, cex = 0.8, adj = 0.5)
        }
      }
    }
    
    grid(nx = NA, ny = NULL, col = "lightgray", lty = "dotted")
  }
  
  par(mar = c(5, 4, 3, 2))
  
  dataset_summary <- valid_correlations %>%
    group_by(dataset, comparison_type) %>%
    summarise(
      mean_correlation = mean(correlation, na.rm = TRUE),
      median_correlation = median(correlation, na.rm = TRUE),
      n_patients = n(),
      .groups = 'drop'
    )
  
  if (nrow(dataset_summary) > 0) {
    dataset_matrix <- dataset_summary %>%
      select(dataset, comparison_type, mean_correlation) %>%
      pivot_wider(names_from = comparison_type, values_from = mean_correlation, values_fill = NA) %>%
      column_to_rownames("dataset") %>%
      as.matrix()
    
    col_order <- intersect(comparison_order, colnames(dataset_matrix))
    if (length(col_order) > 0) {
      dataset_matrix <- dataset_matrix[, col_order, drop = FALSE]
      
      bp2 <- barplot(t(dataset_matrix),
                     main = "Mean Correlations by Dataset",
                     ylab = "Mean Correlation",
                     beside = TRUE,
                     col = colors[colnames(dataset_matrix)],
                     border = "black",
                     las = 1,
                     cex.main = 1.2,
                     ylim = c(min(dataset_matrix, na.rm = TRUE) - 0.1,
                              max(dataset_matrix, na.rm = TRUE) + 0.1))
      
      legend("topright",
             legend = colnames(dataset_matrix),
             fill = colors[colnames(dataset_matrix)],
             cex = 0.8,
             bty = "n")
    }
  }
  
  dev.off()
  
  return(all_correlations)
}

create_summary_statistics <- function(all_correlations) {
  summary_stats <- all_correlations %>%
    filter(!is.na(correlation)) %>%
    group_by(comparison_type) %>%
    summarise(
      n_patients = n(),
      mean_correlation = mean(correlation, na.rm = TRUE),
      median_correlation = median(correlation, na.rm = TRUE),
      sd_correlation = sd(correlation, na.rm = TRUE),
      min_correlation = min(correlation, na.rm = TRUE),
      max_correlation = max(correlation, na.rm = TRUE),
      q25_correlation = quantile(correlation, 0.25, na.rm = TRUE),
      q75_correlation = quantile(correlation, 0.75, na.rm = TRUE),
      .groups = 'drop'
    )
  
  write.csv(summary_stats, file.path(output_dir, "summary_statistics.csv"), row.names = FALSE)
  
  dataset_summary <- all_correlations %>%
    filter(!is.na(correlation)) %>%
    group_by(dataset, comparison_type) %>%
    summarise(
      n_patients = n(),
      mean_correlation = mean(correlation, na.rm = TRUE),
      median_correlation = median(correlation, na.rm = TRUE),
      sd_correlation = sd(correlation, na.rm = TRUE),
      .groups = 'drop'
    )
  
  write.csv(dataset_summary, file.path(output_dir, "dataset_summary.csv"), row.names = FALSE)
  
  method_stats <- all_correlations %>%
    filter(!is.na(correlation)) %>%
    group_by(correlation_method) %>%
    summarise(
      n_patients = n(),
      mean_correlation = mean(correlation, na.rm = TRUE),
      median_correlation = median(correlation, na.rm = TRUE),
      correlation_range = paste(round(min(correlation, na.rm = TRUE), 3), 
                                "to", 
                                round(max(correlation, na.rm = TRUE), 3)),
      .groups = 'drop'
    )
  
  write.csv(method_stats, file.path(output_dir, "correlation_method_statistics.csv"), row.names = FALSE)
  
  stat_tests <- data.frame()
  
  valid_data <- all_correlations[!is.na(all_correlations$correlation), ]
  within_data <- valid_data$correlation[valid_data$comparison_type == "Within Patient\n(CD45+ vs CD45-)"]
  between_pos_data <- valid_data$correlation[valid_data$comparison_type == "Between Patients\n(CD45+ vs CD45+)"]
  between_mixed_data <- valid_data$correlation[valid_data$comparison_type == "Between Patients\n(CD45+ vs CD45-)"]
  
  if (length(within_data) > 0 && length(between_pos_data) > 0) {
    test1 <- wilcox.test(within_data, between_pos_data, alternative = "two.sided")
    
    stat_tests <- rbind(stat_tests, data.frame(
      comparison = "Within vs Between (CD45+ vs CD45+)",
      test_type = "Wilcoxon rank-sum",
      p_value = test1$p.value,
      statistic = as.numeric(test1$statistic),
      n1 = length(within_data),
      n2 = length(between_pos_data),
      mean_diff = mean(within_data, na.rm = TRUE) - mean(between_pos_data, na.rm = TRUE),
      effect_interpretation = ifelse(test1$p.value < 0.05, 
                                     ifelse(mean(within_data) > mean(between_pos_data), 
                                            "Within significantly higher", 
                                            "Between significantly higher"),
                                     "No significant difference"),
      stringsAsFactors = FALSE
    ))
  }
  
  if (length(within_data) > 0 && length(between_mixed_data) > 0) {
    test2 <- wilcox.test(within_data, between_mixed_data, alternative = "two.sided")
    
    stat_tests <- rbind(stat_tests, data.frame(
      comparison = "Within vs Between (CD45+ vs CD45-)",
      test_type = "Wilcoxon rank-sum",
      p_value = test2$p.value,
      statistic = as.numeric(test2$statistic),
      n1 = length(within_data),
      n2 = length(between_mixed_data),
      mean_diff = mean(within_data, na.rm = TRUE) - mean(between_mixed_data, na.rm = TRUE),
      effect_interpretation = ifelse(test2$p.value < 0.05, 
                                     ifelse(mean(within_data) > mean(between_mixed_data), 
                                            "Within significantly higher", 
                                            "Between significantly higher"),
                                     "No significant difference"),
      stringsAsFactors = FALSE
    ))
  }
  
  if (length(between_pos_data) > 0 && length(between_mixed_data) > 0) {
    test3 <- wilcox.test(between_pos_data, between_mixed_data, alternative = "two.sided")
    
    stat_tests <- rbind(stat_tests, data.frame(
      comparison = "Between (CD45+ vs CD45+) vs Between (CD45+ vs CD45-)",
      test_type = "Wilcoxon rank-sum",
      p_value = test3$p.value,
      statistic = as.numeric(test3$statistic),
      n1 = length(between_pos_data),
      n2 = length(between_mixed_data),
      mean_diff = mean(between_pos_data, na.rm = TRUE) - mean(between_mixed_data, na.rm = TRUE),
      effect_interpretation = ifelse(test3$p.value < 0.05, 
                                     ifelse(mean(between_pos_data) > mean(between_mixed_data), 
                                            "Same cell type higher", 
                                            "Mixed cell type higher"),
                                     "No significant difference"),
      stringsAsFactors = FALSE
    ))
  }
  
  write.csv(stat_tests, file.path(output_dir, "statistical_tests.csv"), row.names = FALSE)
  
  return(list(
    summary_stats = summary_stats, 
    dataset_summary = dataset_summary, 
    method_stats = method_stats, 
    stat_tests = stat_tests
  ))
}

main_cnv_correlation_analysis <- function() {
  all_patient_profiles <- create_patient_profiles()
  
  if (length(all_patient_profiles) == 0) {
    return(NULL)
  }
  
  all_correlations <- calculate_neoplastic_CD45_patient_correlations(all_patient_profiles)
  
  if (nrow(all_correlations) == 0) {
    return(NULL)
  }
  
  plot_data <- create_neoplastic_CD45_patient_correlations_plot(all_correlations)
  stats_results <- create_summary_statistics(all_correlations)
  
  return(list(
    patient_profiles = all_patient_profiles,
    correlations = all_correlations,
    statistics = stats_results
  ))
}

main_results <- main_cnv_correlation_analysis()






