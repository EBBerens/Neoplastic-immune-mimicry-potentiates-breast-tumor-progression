#####################################################################################################################
#                                PANNTHR Dataset Breast Tumor Analysis Steps                                        #
#####################################################################################################################
# Step 5: Run InferCNV on Mammary Epithelial Compartment using Normal Mammary Reference and Subset Neoplastic Cells #
# Step 6: Merge Neoplastic Cells from Different Tumors into Combined Object                                         #
#####################################################################################################################


###########################################################################################################
# Step 5: Run InferCNV on Breast Tumors using Reduction Mammoplasty Reference and Subset Neoplastic Cells #
###########################################################################################################

# Load Libraries
library(Seurat)
library(dplyr)
library(patchwork)
library(infercnv)
library(igraph)
library(leiden)
library(ggplot2)

#########################################################################################################
######################################### Inferred CNV Workflow #########################################
#########################################################################################################
# Pan-keratin reduction mammoplasty epithelial reference from Navin and Visvader datasets is used here  #
# Seurat object from tumor is merged with the aforementioned reference                                  #
# CNVs are inferred in epithelial cells within the tumor specimen                                       #
# Normal mammary gland reference is removed from the Seurat object                                      #
# Tumor epithelia is considered neoplastic if it harbors at least one chromosome with an inferred CNV   # 
# Using hg19 file: (https://data.broadinstitute.org/Trinity/CTAT/cnv/)                                  #
#########################################################################################################

###################################################################################################
#############################     InferCNV: PANNTHR Tumor Specimens      ##########################
###################################################################################################


##############################################
#### InferCNV: PANNTHR_Pt1_PreTreat_panCK ####
##############################################

########## Merge with normal mammary glad, extract epithelium, and process via Seurat
# Load RDS
PANNTHR_Pt1_PreTreat_panCK <- readRDS(file = "/R/R_PANNTHR/PANNTHR_RDS/RDS_panCK/PANNTHR_Pt1_PreTreat_panCK.rds")
Visvader_merged_NORM_sampled_panCK <- readRDS(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Visvader_merged_NORM_sampled_panCK.rds")

# Set identities 
colnames(x = PANNTHR_Pt1_PreTreat_panCK[[]])
Idents(object = PANNTHR_Pt1_PreTreat_panCK) <- 'Tumor'
PANNTHR_Pt1_PreTreat_panCK[["cnv.ident"]] <- Idents(object = PANNTHR_Pt1_PreTreat_panCK)
levels(x = PANNTHR_Pt1_PreTreat_panCK)

colnames(x = Visvader_merged_NORM_sampled_panCK[[]])
Idents(object = Visvader_merged_NORM_sampled_panCK) <- 'Normal'
Visvader_merged_NORM_sampled_panCK[["cnv.ident"]] <- Idents(object = Visvader_merged_NORM_sampled_panCK)
levels(x = Visvader_merged_NORM_sampled_panCK)

# Merge subsets into recombined RDS
PANNTHR_Pt1_PreTreat_Visvader_merged_NORM_sampled_panCK <- merge(x = PANNTHR_Pt1_PreTreat_panCK, y = Visvader_merged_NORM_sampled_panCK)
# Join Layers
PANNTHR_Pt1_PreTreat_Visvader_merged_NORM_sampled_panCK <- JoinLayers(PANNTHR_Pt1_PreTreat_Visvader_merged_NORM_sampled_panCK)

# Remove subsetted objects
rm(PANNTHR_Pt1_PreTreat_panCK)
rm(Visvader_merged_NORM_sampled_panCK)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
PANNTHR_Pt1_PreTreat_Visvader_merged_NORM_sampled_panCK <- FindVariableFeatures(PANNTHR_Pt1_PreTreat_Visvader_merged_NORM_sampled_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(PANNTHR_Pt1_PreTreat_Visvader_merged_NORM_sampled_panCK)
PANNTHR_Pt1_PreTreat_Visvader_merged_NORM_sampled_panCK <- ScaleData(PANNTHR_Pt1_PreTreat_Visvader_merged_NORM_sampled_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
PANNTHR_Pt1_PreTreat_Visvader_merged_NORM_sampled_panCK <- RunPCA(PANNTHR_Pt1_PreTreat_Visvader_merged_NORM_sampled_panCK, features = VariableFeatures(object = PANNTHR_Pt1_PreTreat_Visvader_merged_NORM_sampled_panCK))

# Examine and visualize PCA results a few different ways
print(PANNTHR_Pt1_PreTreat_Visvader_merged_NORM_sampled_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(PANNTHR_Pt1_PreTreat_Visvader_merged_NORM_sampled_panCK, dims = 1:2, reduction = "pca")
DimPlot(PANNTHR_Pt1_PreTreat_Visvader_merged_NORM_sampled_panCK, reduction = "pca")

PANNTHR_Pt1_PreTreat_Visvader_merged_NORM_sampled_panCK <- FindNeighbors(PANNTHR_Pt1_PreTreat_Visvader_merged_NORM_sampled_panCK, dims = 1:20)
PANNTHR_Pt1_PreTreat_Visvader_merged_NORM_sampled_panCK <- FindClusters(PANNTHR_Pt1_PreTreat_Visvader_merged_NORM_sampled_panCK, resolution = 0.1)

# Umap clustering
PANNTHR_Pt1_PreTreat_Visvader_merged_NORM_sampled_panCK <- RunUMAP(PANNTHR_Pt1_PreTreat_Visvader_merged_NORM_sampled_panCK, dims = 1:20)
DimPlot(PANNTHR_Pt1_PreTreat_Visvader_merged_NORM_sampled_panCK, reduction = "umap")
DimPlot(PANNTHR_Pt1_PreTreat_Visvader_merged_NORM_sampled_panCK, reduction = "umap", group.by = "orig.ident")

# Save RDS
saveRDS(PANNTHR_Pt1_PreTreat_Visvader_merged_NORM_sampled_panCK, file = "/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Input/PANNTHR_inferCNV_Input_RDS/PANNTHR_Pt1_PreTreat_Visvader_merged_NORM_sampled_panCK.rds")

########## Prepare Cell Annotations
# Extract & Save Meta Data
PANNTHR_Pt1_PreTreat_Visvader_merged_NORM_sampled_panCK.meta.data <- as.data.frame(as.matrix(PANNTHR_Pt1_PreTreat_Visvader_merged_NORM_sampled_panCK@meta.data))
# Extract pertinent columns
# Make row names first column
PANNTHR_Pt1_PreTreat_Visvader_merged_NORM_sampled_panCK.meta.data <- tibble::rownames_to_column(PANNTHR_Pt1_PreTreat_Visvader_merged_NORM_sampled_panCK.meta.data, "Cells")
# In Seurat, "X" becomes the cell names while "cnv.ident" is specified above
PANNTHR_Pt1_PreTreat_Visvader_merged_NORM_sampled_panCK.annotations <- PANNTHR_Pt1_PreTreat_Visvader_merged_NORM_sampled_panCK.meta.data[, c("Cells", "cnv.ident")]
# Delete header
names(PANNTHR_Pt1_PreTreat_Visvader_merged_NORM_sampled_panCK.annotations) <- NULL
# Save as table
write.table(PANNTHR_Pt1_PreTreat_Visvader_merged_NORM_sampled_panCK.annotations,"/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Input/PANNTHR_Pt1_PreTreat_Visvader_merged_NORM_sampled_panCK.annotations.txt",sep="\t",row.names=FALSE)
# Remove metadata and annotations
rm(PANNTHR_Pt1_PreTreat_Visvader_merged_NORM_sampled_panCK.meta.data)
rm(PANNTHR_Pt1_PreTreat_Visvader_merged_NORM_sampled_panCK.annotations)

########## Generate Count Matrix
# Extract counts from the Seurat object
PANNTHR_Pt1_PreTreat_Visvader_merged_NORM_sampled_panCK.count <- PANNTHR_Pt1_PreTreat_Visvader_merged_NORM_sampled_panCK[["RNA"]]$counts
# Convert counts to a sparse matrix
PANNTHR_Pt1_PreTreat_Visvader_merged_NORM_sampled_panCK.count.matrix <- as(PANNTHR_Pt1_PreTreat_Visvader_merged_NORM_sampled_panCK.count, 'sparseMatrix')
rm(PANNTHR_Pt1_PreTreat_Visvader_merged_NORM_sampled_panCK.count)

# Convert to requisite data format
PANNTHR_Pt1_PreTreat_Visvader_merged_NORM_sampled_panCK.count.matrix <- as.data.frame(PANNTHR_Pt1_PreTreat_Visvader_merged_NORM_sampled_panCK.count.matrix)
# Save table
write.csv(PANNTHR_Pt1_PreTreat_Visvader_merged_NORM_sampled_panCK.count.matrix, file = "/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Input/PANNTHR_inferCNV_Input_Count_Matrix/PANNTHR_Pt1_PreTreat_Visvader_merged_NORM_sampled_panCK.count.matrix.csv")
# Remove Objects
rm(PANNTHR_Pt1_PreTreat_Visvader_merged_NORM_sampled_panCK)
rm(PANNTHR_Pt1_PreTreat_Visvader_merged_NORM_sampled_panCK.count.matrix)

########### Run inferCNV
# Note: check.names = FALSE prevents cell names from having dashes replaced with dots; essential to match annotations file
PANNTHR_Pt1_PreTreat_inferCNV <- CreateInfercnvObject(raw_counts_matrix=read.csv(file = "/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Input/PANNTHR_inferCNV_Input_Count_Matrix/PANNTHR_Pt1_PreTreat_Visvader_merged_NORM_sampled_panCK.count.matrix.csv", header = TRUE, row.names = 1,check.names = FALSE),
                                                      annotations_file=read.table(file = "/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Input/PANNTHR_Pt1_PreTreat_Visvader_merged_NORM_sampled_panCK.annotations.txt", header = FALSE, row.names = 1),
                                                      delim="\t",
                                                      gene_order_file="/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Input/hg38_gencode_v27.txt",
                                                      ref_group_names=c("Normal"))


# perform infercnv operations to reveal cnv signal
PANNTHR_Pt1_PreTreat_inferCNV = infercnv::run(PANNTHR_Pt1_PreTreat_inferCNV,
                                              cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                              out_dir="/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Output/PANNTHR_Pt1_PreTreat_panCK/",  # dir is auto-created for storing outputs
                                              cluster_by_groups=T,   # cluster
                                              denoise=T,
                                              HMM=T,
                                              num_ref_groups = 8, analysis_mode = "subclusters")

###########  Integrate data with Seurat object
# Read RDS and add inferCNV results to metadata
PANNTHR_Pt1_PreTreat_Visvader_merged_NORM_sampled_panCK <- readRDS(file = "/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Input/PANNTHR_inferCNV_Input_RDS/PANNTHR_Pt1_PreTreat_Visvader_merged_NORM_sampled_panCK.rds") 
PANNTHR_Pt1_PreTreat_panCK_cnv <- infercnv::add_to_seurat(infercnv_output_path="/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Output/PANNTHR_Pt1_PreTreat_panCK/",
                                                          seurat_obj=PANNTHR_Pt1_PreTreat_Visvader_merged_NORM_sampled_panCK, # optional
                                                          top_n=10)

# Collect Tumor cells from Seurat Object
PANNTHR_Pt1_PreTreat_panCK_cnv <- subset(x = PANNTHR_Pt1_PreTreat_panCK_cnv, subset = cnv.ident == "Tumor")
# Quantify chromosomes with alterations per cell
PANNTHR_Pt1_PreTreat_panCK_cnv_Pos_Chromosomes <- as.data.frame(rowSums(PANNTHR_Pt1_PreTreat_panCK_cnv@meta.data[,startsWith(names(PANNTHR_Pt1_PreTreat_panCK_cnv@meta.data),"has_cnv_")] > 0 ))
PANNTHR_Pt1_PreTreat_panCK_cnv_Neg_Chromosomes <- as.data.frame(rowSums(PANNTHR_Pt1_PreTreat_panCK_cnv@meta.data[,startsWith(names(PANNTHR_Pt1_PreTreat_panCK_cnv@meta.data),"has_cnv_")] == 0))
PANNTHR_Pt1_PreTreat_panCK_cnv_Chromosome_Quant <- as.data.frame(c(PANNTHR_Pt1_PreTreat_panCK_cnv_Pos_Chromosomes, PANNTHR_Pt1_PreTreat_panCK_cnv_Neg_Chromosomes))
colnames(PANNTHR_Pt1_PreTreat_panCK_cnv_Chromosome_Quant) <- c("Chromosomes CNV Pos", "Chromosomes CNV Neg")
write.csv(PANNTHR_Pt1_PreTreat_panCK_cnv_Chromosome_Quant, file = "/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Output/PANNTHR_Pt1_PreTreat_panCK/PANNTHR_Pt1_PreTreat_panCK_cnv_Chromosome_Quant.csv")

# Save CNV Metadata
PANNTHR_Pt1_PreTreat_panCK_cnv.meta.data <- as.data.frame(as.matrix(PANNTHR_Pt1_PreTreat_panCK_cnv@meta.data))
write.csv(PANNTHR_Pt1_PreTreat_panCK_cnv.meta.data, file = "/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Output/PANNTHR_Pt1_PreTreat_panCK/PANNTHR_Pt1_PreTreat_panCK_cnv.meta.data.csv")

# Subset Cells with inferred CNVs -------------------------------------------------------------------------------------------------------------
PANNTHR_Pt1_PreTreat_panCK_cnv <- subset(PANNTHR_Pt1_PreTreat_panCK_cnv, cells=rownames(PANNTHR_Pt1_PreTreat_panCK_cnv@meta.data)[rowSums(PANNTHR_Pt1_PreTreat_panCK_cnv@meta.data[,startsWith(names(PANNTHR_Pt1_PreTreat_panCK_cnv@meta.data),"has_cnv_")]) > 0 ])

# Save RDS with inferredCNVs
saveRDS(PANNTHR_Pt1_PreTreat_panCK_cnv, file = "/R/R_PANNTHR/PANNTHR_RDS/RDS_inferCNV/PANNTHR_Pt1_PreTreat_panCK_cnv.rds")


# Free memory
rm(PANNTHR_Pt1_PreTreat_inferCNV)
rm(PANNTHR_Pt1_PreTreat_Visvader_merged_NORM_sampled_panCK)
rm(PANNTHR_Pt1_PreTreat_panCK_cnv)
rm(PANNTHR_Pt1_PreTreat_panCK_cnv_Pos_Chromosomes)
rm(PANNTHR_Pt1_PreTreat_panCK_cnv_Neg_Chromosomes)
rm(PANNTHR_Pt1_PreTreat_panCK_cnv_Chromosome_Quant)
rm(PANNTHR_Pt1_PreTreat_panCK_cnv.meta.data)
gc()



##############################################
#### InferCNV: PANNTHR_Pt2_PreTreat_panCK ####
##############################################

########## Merge with normal mammary glad, extract epithelium, and process via Seurat
# Load RDS
PANNTHR_Pt2_PreTreat_panCK <- readRDS(file = "/R/R_PANNTHR/PANNTHR_RDS/RDS_panCK/PANNTHR_Pt2_PreTreat_panCK.rds")
Visvader_merged_NORM_sampled_panCK <- readRDS(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Visvader_merged_NORM_sampled_panCK.rds")

# Set identities 
colnames(x = PANNTHR_Pt2_PreTreat_panCK[[]])
Idents(object = PANNTHR_Pt2_PreTreat_panCK) <- 'Tumor'
PANNTHR_Pt2_PreTreat_panCK[["cnv.ident"]] <- Idents(object = PANNTHR_Pt2_PreTreat_panCK)
levels(x = PANNTHR_Pt2_PreTreat_panCK)

colnames(x = Visvader_merged_NORM_sampled_panCK[[]])
Idents(object = Visvader_merged_NORM_sampled_panCK) <- 'Normal'
Visvader_merged_NORM_sampled_panCK[["cnv.ident"]] <- Idents(object = Visvader_merged_NORM_sampled_panCK)
levels(x = Visvader_merged_NORM_sampled_panCK)

# Merge subsets into recombined RDS
PANNTHR_Pt2_PreTreat_Visvader_merged_NORM_sampled_panCK <- merge(x = PANNTHR_Pt2_PreTreat_panCK, y = Visvader_merged_NORM_sampled_panCK)
# Join Layers
PANNTHR_Pt2_PreTreat_Visvader_merged_NORM_sampled_panCK <- JoinLayers(PANNTHR_Pt2_PreTreat_Visvader_merged_NORM_sampled_panCK)

# Remove subsetted objects
rm(PANNTHR_Pt2_PreTreat_panCK)
rm(Visvader_merged_NORM_sampled_panCK)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
PANNTHR_Pt2_PreTreat_Visvader_merged_NORM_sampled_panCK <- FindVariableFeatures(PANNTHR_Pt2_PreTreat_Visvader_merged_NORM_sampled_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(PANNTHR_Pt2_PreTreat_Visvader_merged_NORM_sampled_panCK)
PANNTHR_Pt2_PreTreat_Visvader_merged_NORM_sampled_panCK <- ScaleData(PANNTHR_Pt2_PreTreat_Visvader_merged_NORM_sampled_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
PANNTHR_Pt2_PreTreat_Visvader_merged_NORM_sampled_panCK <- RunPCA(PANNTHR_Pt2_PreTreat_Visvader_merged_NORM_sampled_panCK, features = VariableFeatures(object = PANNTHR_Pt2_PreTreat_Visvader_merged_NORM_sampled_panCK))

# Examine and visualize PCA results a few different ways
print(PANNTHR_Pt2_PreTreat_Visvader_merged_NORM_sampled_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(PANNTHR_Pt2_PreTreat_Visvader_merged_NORM_sampled_panCK, dims = 1:2, reduction = "pca")
DimPlot(PANNTHR_Pt2_PreTreat_Visvader_merged_NORM_sampled_panCK, reduction = "pca")

PANNTHR_Pt2_PreTreat_Visvader_merged_NORM_sampled_panCK <- FindNeighbors(PANNTHR_Pt2_PreTreat_Visvader_merged_NORM_sampled_panCK, dims = 1:20)
PANNTHR_Pt2_PreTreat_Visvader_merged_NORM_sampled_panCK <- FindClusters(PANNTHR_Pt2_PreTreat_Visvader_merged_NORM_sampled_panCK, resolution = 0.1)

# Umap clustering
PANNTHR_Pt2_PreTreat_Visvader_merged_NORM_sampled_panCK <- RunUMAP(PANNTHR_Pt2_PreTreat_Visvader_merged_NORM_sampled_panCK, dims = 1:20)
DimPlot(PANNTHR_Pt2_PreTreat_Visvader_merged_NORM_sampled_panCK, reduction = "umap")
DimPlot(PANNTHR_Pt2_PreTreat_Visvader_merged_NORM_sampled_panCK, reduction = "umap", group.by = "orig.ident")

# Save RDS
saveRDS(PANNTHR_Pt2_PreTreat_Visvader_merged_NORM_sampled_panCK, file = "/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Input/PANNTHR_inferCNV_Input_RDS/PANNTHR_Pt2_PreTreat_Visvader_merged_NORM_sampled_panCK.rds")

########## Prepare Cell Annotations
# Extract & Save Meta Data
PANNTHR_Pt2_PreTreat_Visvader_merged_NORM_sampled_panCK.meta.data <- as.data.frame(as.matrix(PANNTHR_Pt2_PreTreat_Visvader_merged_NORM_sampled_panCK@meta.data))
# Extract pertinent columns
# Make row names first column
PANNTHR_Pt2_PreTreat_Visvader_merged_NORM_sampled_panCK.meta.data <- tibble::rownames_to_column(PANNTHR_Pt2_PreTreat_Visvader_merged_NORM_sampled_panCK.meta.data, "Cells")
# In Seurat, "X" becomes the cell names while "cnv.ident" is specified above
PANNTHR_Pt2_PreTreat_Visvader_merged_NORM_sampled_panCK.annotations <- PANNTHR_Pt2_PreTreat_Visvader_merged_NORM_sampled_panCK.meta.data[, c("Cells", "cnv.ident")]
# Delete header
names(PANNTHR_Pt2_PreTreat_Visvader_merged_NORM_sampled_panCK.annotations) <- NULL
# Save as table
write.table(PANNTHR_Pt2_PreTreat_Visvader_merged_NORM_sampled_panCK.annotations,"/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Input/PANNTHR_Pt2_PreTreat_Visvader_merged_NORM_sampled_panCK.annotations.txt",sep="\t",row.names=FALSE)
# Remove metadata and annotations
rm(PANNTHR_Pt2_PreTreat_Visvader_merged_NORM_sampled_panCK.meta.data)
rm(PANNTHR_Pt2_PreTreat_Visvader_merged_NORM_sampled_panCK.annotations)

########## Generate Count Matrix
# Extract counts from the Seurat object
PANNTHR_Pt2_PreTreat_Visvader_merged_NORM_sampled_panCK.count <- PANNTHR_Pt2_PreTreat_Visvader_merged_NORM_sampled_panCK[["RNA"]]$counts
# Convert counts to a sparse matrix
PANNTHR_Pt2_PreTreat_Visvader_merged_NORM_sampled_panCK.count.matrix <- as(PANNTHR_Pt2_PreTreat_Visvader_merged_NORM_sampled_panCK.count, 'sparseMatrix')
rm(PANNTHR_Pt2_PreTreat_Visvader_merged_NORM_sampled_panCK.count)

# Convert to requisite data format
PANNTHR_Pt2_PreTreat_Visvader_merged_NORM_sampled_panCK.count.matrix <- as.data.frame(PANNTHR_Pt2_PreTreat_Visvader_merged_NORM_sampled_panCK.count.matrix)
# Save table
write.csv(PANNTHR_Pt2_PreTreat_Visvader_merged_NORM_sampled_panCK.count.matrix, file = "/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Input/PANNTHR_inferCNV_Input_Count_Matrix/PANNTHR_Pt2_PreTreat_Visvader_merged_NORM_sampled_panCK.count.matrix.csv")
# Remove Objects
rm(PANNTHR_Pt2_PreTreat_Visvader_merged_NORM_sampled_panCK)
rm(PANNTHR_Pt2_PreTreat_Visvader_merged_NORM_sampled_panCK.count.matrix)

########### Run inferCNV
# Note: check.names = FALSE prevents cell names from having dashes replaced with dots; essential to match annotations file
PANNTHR_Pt2_PreTreat_inferCNV <- CreateInfercnvObject(raw_counts_matrix=read.csv(file = "/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Input/PANNTHR_inferCNV_Input_Count_Matrix/PANNTHR_Pt2_PreTreat_Visvader_merged_NORM_sampled_panCK.count.matrix.csv", header = TRUE, row.names = 1,check.names = FALSE),
                                                      annotations_file=read.table(file = "/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Input/PANNTHR_Pt2_PreTreat_Visvader_merged_NORM_sampled_panCK.annotations.txt", header = FALSE, row.names = 1),
                                                      delim="\t",
                                                      gene_order_file="/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Input/hg38_gencode_v27.txt",
                                                      ref_group_names=c("Normal"))


# perform infercnv operations to reveal cnv signal
PANNTHR_Pt2_PreTreat_inferCNV = infercnv::run(PANNTHR_Pt2_PreTreat_inferCNV,
                                              cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                              out_dir="/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Output/PANNTHR_Pt2_PreTreat_panCK/",  # dir is auto-created for storing outputs
                                              cluster_by_groups=T,   # cluster
                                              denoise=T,
                                              HMM=T,
                                              num_ref_groups = 8, analysis_mode = "subclusters")

###########  Integrate data with Seurat object
# Read RDS and add inferCNV results to metadata
PANNTHR_Pt2_PreTreat_Visvader_merged_NORM_sampled_panCK <- readRDS(file = "/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Input/PANNTHR_inferCNV_Input_RDS/PANNTHR_Pt2_PreTreat_Visvader_merged_NORM_sampled_panCK.rds") 
PANNTHR_Pt2_PreTreat_panCK_cnv <- infercnv::add_to_seurat(infercnv_output_path="/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Output/PANNTHR_Pt2_PreTreat_panCK/",
                                                          seurat_obj=PANNTHR_Pt2_PreTreat_Visvader_merged_NORM_sampled_panCK, # optional
                                                          top_n=10)

# Collect Tumor cells from Seurat Object
PANNTHR_Pt2_PreTreat_panCK_cnv <- subset(x = PANNTHR_Pt2_PreTreat_panCK_cnv, subset = cnv.ident == "Tumor")
# Quantify chromosomes with alterations per cell
PANNTHR_Pt2_PreTreat_panCK_cnv_Pos_Chromosomes <- as.data.frame(rowSums(PANNTHR_Pt2_PreTreat_panCK_cnv@meta.data[,startsWith(names(PANNTHR_Pt2_PreTreat_panCK_cnv@meta.data),"has_cnv_")] > 0 ))
PANNTHR_Pt2_PreTreat_panCK_cnv_Neg_Chromosomes <- as.data.frame(rowSums(PANNTHR_Pt2_PreTreat_panCK_cnv@meta.data[,startsWith(names(PANNTHR_Pt2_PreTreat_panCK_cnv@meta.data),"has_cnv_")] == 0))
PANNTHR_Pt2_PreTreat_panCK_cnv_Chromosome_Quant <- as.data.frame(c(PANNTHR_Pt2_PreTreat_panCK_cnv_Pos_Chromosomes, PANNTHR_Pt2_PreTreat_panCK_cnv_Neg_Chromosomes))
colnames(PANNTHR_Pt2_PreTreat_panCK_cnv_Chromosome_Quant) <- c("Chromosomes CNV Pos", "Chromosomes CNV Neg")
write.csv(PANNTHR_Pt2_PreTreat_panCK_cnv_Chromosome_Quant, file = "/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Output/PANNTHR_Pt2_PreTreat_panCK/PANNTHR_Pt2_PreTreat_panCK_cnv_Chromosome_Quant.csv")

# Save CNV Metadata
PANNTHR_Pt2_PreTreat_panCK_cnv.meta.data <- as.data.frame(as.matrix(PANNTHR_Pt2_PreTreat_panCK_cnv@meta.data))
write.csv(PANNTHR_Pt2_PreTreat_panCK_cnv.meta.data, file = "/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Output/PANNTHR_Pt2_PreTreat_panCK/PANNTHR_Pt2_PreTreat_panCK_cnv.meta.data.csv")

# Subset Cells with inferred CNVs -------------------------------------------------------------------------------------------------------------
PANNTHR_Pt2_PreTreat_panCK_cnv <- subset(PANNTHR_Pt2_PreTreat_panCK_cnv, cells=rownames(PANNTHR_Pt2_PreTreat_panCK_cnv@meta.data)[rowSums(PANNTHR_Pt2_PreTreat_panCK_cnv@meta.data[,startsWith(names(PANNTHR_Pt2_PreTreat_panCK_cnv@meta.data),"has_cnv_")]) > 0 ])

# Save RDS with inferredCNVs
saveRDS(PANNTHR_Pt2_PreTreat_panCK_cnv, file = "/R/R_PANNTHR/PANNTHR_RDS/RDS_inferCNV/PANNTHR_Pt2_PreTreat_panCK_cnv.rds")


# Free memory
rm(PANNTHR_Pt2_PreTreat_inferCNV)
rm(PANNTHR_Pt2_PreTreat_Visvader_merged_NORM_sampled_panCK)
rm(PANNTHR_Pt2_PreTreat_panCK_cnv)
rm(PANNTHR_Pt2_PreTreat_panCK_cnv_Pos_Chromosomes)
rm(PANNTHR_Pt2_PreTreat_panCK_cnv_Neg_Chromosomes)
rm(PANNTHR_Pt2_PreTreat_panCK_cnv_Chromosome_Quant)
rm(PANNTHR_Pt2_PreTreat_panCK_cnv.meta.data)
gc()



##############################################
#### InferCNV: PANNTHR_Pt4_PreTreat_panCK ####
##############################################

########## Merge with normal mammary glad, extract epithelium, and process via Seurat
# Load RDS
PANNTHR_Pt4_PreTreat_panCK <- readRDS(file = "/R/R_PANNTHR/PANNTHR_RDS/RDS_panCK/PANNTHR_Pt4_PreTreat_panCK.rds")
Visvader_merged_NORM_sampled_panCK <- readRDS(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Visvader_merged_NORM_sampled_panCK.rds")

# Set identities 
colnames(x = PANNTHR_Pt4_PreTreat_panCK[[]])
Idents(object = PANNTHR_Pt4_PreTreat_panCK) <- 'Tumor'
PANNTHR_Pt4_PreTreat_panCK[["cnv.ident"]] <- Idents(object = PANNTHR_Pt4_PreTreat_panCK)
levels(x = PANNTHR_Pt4_PreTreat_panCK)

colnames(x = Visvader_merged_NORM_sampled_panCK[[]])
Idents(object = Visvader_merged_NORM_sampled_panCK) <- 'Normal'
Visvader_merged_NORM_sampled_panCK[["cnv.ident"]] <- Idents(object = Visvader_merged_NORM_sampled_panCK)
levels(x = Visvader_merged_NORM_sampled_panCK)

# Merge subsets into recombined RDS
PANNTHR_Pt4_PreTreat_Visvader_merged_NORM_sampled_panCK <- merge(x = PANNTHR_Pt4_PreTreat_panCK, y = Visvader_merged_NORM_sampled_panCK)
# Join Layers
PANNTHR_Pt4_PreTreat_Visvader_merged_NORM_sampled_panCK <- JoinLayers(PANNTHR_Pt4_PreTreat_Visvader_merged_NORM_sampled_panCK)

# Remove subsetted objects
rm(PANNTHR_Pt4_PreTreat_panCK)
rm(Visvader_merged_NORM_sampled_panCK)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
PANNTHR_Pt4_PreTreat_Visvader_merged_NORM_sampled_panCK <- FindVariableFeatures(PANNTHR_Pt4_PreTreat_Visvader_merged_NORM_sampled_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(PANNTHR_Pt4_PreTreat_Visvader_merged_NORM_sampled_panCK)
PANNTHR_Pt4_PreTreat_Visvader_merged_NORM_sampled_panCK <- ScaleData(PANNTHR_Pt4_PreTreat_Visvader_merged_NORM_sampled_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
PANNTHR_Pt4_PreTreat_Visvader_merged_NORM_sampled_panCK <- RunPCA(PANNTHR_Pt4_PreTreat_Visvader_merged_NORM_sampled_panCK, features = VariableFeatures(object = PANNTHR_Pt4_PreTreat_Visvader_merged_NORM_sampled_panCK))

# Examine and visualize PCA results a few different ways
print(PANNTHR_Pt4_PreTreat_Visvader_merged_NORM_sampled_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(PANNTHR_Pt4_PreTreat_Visvader_merged_NORM_sampled_panCK, dims = 1:2, reduction = "pca")
DimPlot(PANNTHR_Pt4_PreTreat_Visvader_merged_NORM_sampled_panCK, reduction = "pca")

PANNTHR_Pt4_PreTreat_Visvader_merged_NORM_sampled_panCK <- FindNeighbors(PANNTHR_Pt4_PreTreat_Visvader_merged_NORM_sampled_panCK, dims = 1:20)
PANNTHR_Pt4_PreTreat_Visvader_merged_NORM_sampled_panCK <- FindClusters(PANNTHR_Pt4_PreTreat_Visvader_merged_NORM_sampled_panCK, resolution = 0.1)

# Umap clustering
PANNTHR_Pt4_PreTreat_Visvader_merged_NORM_sampled_panCK <- RunUMAP(PANNTHR_Pt4_PreTreat_Visvader_merged_NORM_sampled_panCK, dims = 1:20)
DimPlot(PANNTHR_Pt4_PreTreat_Visvader_merged_NORM_sampled_panCK, reduction = "umap")
DimPlot(PANNTHR_Pt4_PreTreat_Visvader_merged_NORM_sampled_panCK, reduction = "umap", group.by = "orig.ident")

# Save RDS
saveRDS(PANNTHR_Pt4_PreTreat_Visvader_merged_NORM_sampled_panCK, file = "/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Input/PANNTHR_inferCNV_Input_RDS/PANNTHR_Pt4_PreTreat_Visvader_merged_NORM_sampled_panCK.rds")

########## Prepare Cell Annotations
# Extract & Save Meta Data
PANNTHR_Pt4_PreTreat_Visvader_merged_NORM_sampled_panCK.meta.data <- as.data.frame(as.matrix(PANNTHR_Pt4_PreTreat_Visvader_merged_NORM_sampled_panCK@meta.data))
# Extract pertinent columns
# Make row names first column
PANNTHR_Pt4_PreTreat_Visvader_merged_NORM_sampled_panCK.meta.data <- tibble::rownames_to_column(PANNTHR_Pt4_PreTreat_Visvader_merged_NORM_sampled_panCK.meta.data, "Cells")
# In Seurat, "X" becomes the cell names while "cnv.ident" is specified above
PANNTHR_Pt4_PreTreat_Visvader_merged_NORM_sampled_panCK.annotations <- PANNTHR_Pt4_PreTreat_Visvader_merged_NORM_sampled_panCK.meta.data[, c("Cells", "cnv.ident")]
# Delete header
names(PANNTHR_Pt4_PreTreat_Visvader_merged_NORM_sampled_panCK.annotations) <- NULL
# Save as table
write.table(PANNTHR_Pt4_PreTreat_Visvader_merged_NORM_sampled_panCK.annotations,"/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Input/PANNTHR_Pt4_PreTreat_Visvader_merged_NORM_sampled_panCK.annotations.txt",sep="\t",row.names=FALSE)
# Remove metadata and annotations
rm(PANNTHR_Pt4_PreTreat_Visvader_merged_NORM_sampled_panCK.meta.data)
rm(PANNTHR_Pt4_PreTreat_Visvader_merged_NORM_sampled_panCK.annotations)

########## Generate Count Matrix
# Extract counts from the Seurat object
PANNTHR_Pt4_PreTreat_Visvader_merged_NORM_sampled_panCK.count <- PANNTHR_Pt4_PreTreat_Visvader_merged_NORM_sampled_panCK[["RNA"]]$counts
# Convert counts to a sparse matrix
PANNTHR_Pt4_PreTreat_Visvader_merged_NORM_sampled_panCK.count.matrix <- as(PANNTHR_Pt4_PreTreat_Visvader_merged_NORM_sampled_panCK.count, 'sparseMatrix')
rm(PANNTHR_Pt4_PreTreat_Visvader_merged_NORM_sampled_panCK.count)

# Convert to requisite data format
PANNTHR_Pt4_PreTreat_Visvader_merged_NORM_sampled_panCK.count.matrix <- as.data.frame(PANNTHR_Pt4_PreTreat_Visvader_merged_NORM_sampled_panCK.count.matrix)
# Save table
write.csv(PANNTHR_Pt4_PreTreat_Visvader_merged_NORM_sampled_panCK.count.matrix, file = "/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Input/PANNTHR_inferCNV_Input_Count_Matrix/PANNTHR_Pt4_PreTreat_Visvader_merged_NORM_sampled_panCK.count.matrix.csv")
# Remove Objects
rm(PANNTHR_Pt4_PreTreat_Visvader_merged_NORM_sampled_panCK)
rm(PANNTHR_Pt4_PreTreat_Visvader_merged_NORM_sampled_panCK.count.matrix)

########### Run inferCNV
# Note: check.names = FALSE prevents cell names from having dashes replaced with dots; essential to match annotations file
PANNTHR_Pt4_PreTreat_inferCNV <- CreateInfercnvObject(raw_counts_matrix=read.csv(file = "/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Input/PANNTHR_inferCNV_Input_Count_Matrix/PANNTHR_Pt4_PreTreat_Visvader_merged_NORM_sampled_panCK.count.matrix.csv", header = TRUE, row.names = 1,check.names = FALSE),
                                               annotations_file=read.table(file = "/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Input/PANNTHR_Pt4_PreTreat_Visvader_merged_NORM_sampled_panCK.annotations.txt", header = FALSE, row.names = 1),
                                               delim="\t",
                                               gene_order_file="/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Input/hg38_gencode_v27.txt",
                                               ref_group_names=c("Normal"))


# perform infercnv operations to reveal cnv signal
PANNTHR_Pt4_PreTreat_inferCNV = infercnv::run(PANNTHR_Pt4_PreTreat_inferCNV,
                                       cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                       out_dir="/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Output/PANNTHR_Pt4_PreTreat_panCK/",  # dir is auto-created for storing outputs
                                       cluster_by_groups=T,   # cluster
                                       denoise=T,
                                       HMM=T,
                                       num_ref_groups = 8, analysis_mode = "subclusters")

###########  Integrate data with Seurat object
# Read RDS and add inferCNV results to metadata
PANNTHR_Pt4_PreTreat_Visvader_merged_NORM_sampled_panCK <- readRDS(file = "/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Input/PANNTHR_inferCNV_Input_RDS/PANNTHR_Pt4_PreTreat_Visvader_merged_NORM_sampled_panCK.rds") 
PANNTHR_Pt4_PreTreat_panCK_cnv <- infercnv::add_to_seurat(infercnv_output_path="/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Output/PANNTHR_Pt4_PreTreat_panCK/",
                                                   seurat_obj=PANNTHR_Pt4_PreTreat_Visvader_merged_NORM_sampled_panCK, # optional
                                                   top_n=10)

# Collect Tumor cells from Seurat Object
PANNTHR_Pt4_PreTreat_panCK_cnv <- subset(x = PANNTHR_Pt4_PreTreat_panCK_cnv, subset = cnv.ident == "Tumor")
# Quantify chromosomes with alterations per cell
PANNTHR_Pt4_PreTreat_panCK_cnv_Pos_Chromosomes <- as.data.frame(rowSums(PANNTHR_Pt4_PreTreat_panCK_cnv@meta.data[,startsWith(names(PANNTHR_Pt4_PreTreat_panCK_cnv@meta.data),"has_cnv_")] > 0 ))
PANNTHR_Pt4_PreTreat_panCK_cnv_Neg_Chromosomes <- as.data.frame(rowSums(PANNTHR_Pt4_PreTreat_panCK_cnv@meta.data[,startsWith(names(PANNTHR_Pt4_PreTreat_panCK_cnv@meta.data),"has_cnv_")] == 0))
PANNTHR_Pt4_PreTreat_panCK_cnv_Chromosome_Quant <- as.data.frame(c(PANNTHR_Pt4_PreTreat_panCK_cnv_Pos_Chromosomes, PANNTHR_Pt4_PreTreat_panCK_cnv_Neg_Chromosomes))
colnames(PANNTHR_Pt4_PreTreat_panCK_cnv_Chromosome_Quant) <- c("Chromosomes CNV Pos", "Chromosomes CNV Neg")
write.csv(PANNTHR_Pt4_PreTreat_panCK_cnv_Chromosome_Quant, file = "/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Output/PANNTHR_Pt4_PreTreat_panCK/PANNTHR_Pt4_PreTreat_panCK_cnv_Chromosome_Quant.csv")

# Save CNV Metadata
PANNTHR_Pt4_PreTreat_panCK_cnv.meta.data <- as.data.frame(as.matrix(PANNTHR_Pt4_PreTreat_panCK_cnv@meta.data))
write.csv(PANNTHR_Pt4_PreTreat_panCK_cnv.meta.data, file = "/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Output/PANNTHR_Pt4_PreTreat_panCK/PANNTHR_Pt4_PreTreat_panCK_cnv.meta.data.csv")

# Subset Cells with inferred CNVs -------------------------------------------------------------------------------------------------------------
PANNTHR_Pt4_PreTreat_panCK_cnv <- subset(PANNTHR_Pt4_PreTreat_panCK_cnv, cells=rownames(PANNTHR_Pt4_PreTreat_panCK_cnv@meta.data)[rowSums(PANNTHR_Pt4_PreTreat_panCK_cnv@meta.data[,startsWith(names(PANNTHR_Pt4_PreTreat_panCK_cnv@meta.data),"has_cnv_")]) > 0 ])

# Save RDS with inferredCNVs
saveRDS(PANNTHR_Pt4_PreTreat_panCK_cnv, file = "/R/R_PANNTHR/PANNTHR_RDS/RDS_inferCNV/PANNTHR_Pt4_PreTreat_panCK_cnv.rds")


# Free memory
rm(PANNTHR_Pt4_PreTreat_inferCNV)
rm(PANNTHR_Pt4_PreTreat_Visvader_merged_NORM_sampled_panCK)
rm(PANNTHR_Pt4_PreTreat_panCK_cnv)
rm(PANNTHR_Pt4_PreTreat_panCK_cnv_Pos_Chromosomes)
rm(PANNTHR_Pt4_PreTreat_panCK_cnv_Neg_Chromosomes)
rm(PANNTHR_Pt4_PreTreat_panCK_cnv_Chromosome_Quant)
rm(PANNTHR_Pt4_PreTreat_panCK_cnv.meta.data)
gc()



##############################################
#### InferCNV: PANNTHR_Pt5_PreTreat_panCK ####
##############################################

########## Merge with normal mammary glad, extract epithelium, and process via Seurat
# Load RDS
PANNTHR_Pt5_PreTreat_panCK <- readRDS(file = "/R/R_PANNTHR/PANNTHR_RDS/RDS_panCK/PANNTHR_Pt5_PreTreat_panCK.rds")
Visvader_merged_NORM_sampled_panCK <- readRDS(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Visvader_merged_NORM_sampled_panCK.rds")

# Set identities 
colnames(x = PANNTHR_Pt5_PreTreat_panCK[[]])
Idents(object = PANNTHR_Pt5_PreTreat_panCK) <- 'Tumor'
PANNTHR_Pt5_PreTreat_panCK[["cnv.ident"]] <- Idents(object = PANNTHR_Pt5_PreTreat_panCK)
levels(x = PANNTHR_Pt5_PreTreat_panCK)

colnames(x = Visvader_merged_NORM_sampled_panCK[[]])
Idents(object = Visvader_merged_NORM_sampled_panCK) <- 'Normal'
Visvader_merged_NORM_sampled_panCK[["cnv.ident"]] <- Idents(object = Visvader_merged_NORM_sampled_panCK)
levels(x = Visvader_merged_NORM_sampled_panCK)

# Merge subsets into recombined RDS
PANNTHR_Pt5_PreTreat_Visvader_merged_NORM_sampled_panCK <- merge(x = PANNTHR_Pt5_PreTreat_panCK, y = Visvader_merged_NORM_sampled_panCK)
# Join Layers
PANNTHR_Pt5_PreTreat_Visvader_merged_NORM_sampled_panCK <- JoinLayers(PANNTHR_Pt5_PreTreat_Visvader_merged_NORM_sampled_panCK)

# Remove subsetted objects
rm(PANNTHR_Pt5_PreTreat_panCK)
rm(Visvader_merged_NORM_sampled_panCK)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
PANNTHR_Pt5_PreTreat_Visvader_merged_NORM_sampled_panCK <- FindVariableFeatures(PANNTHR_Pt5_PreTreat_Visvader_merged_NORM_sampled_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(PANNTHR_Pt5_PreTreat_Visvader_merged_NORM_sampled_panCK)
PANNTHR_Pt5_PreTreat_Visvader_merged_NORM_sampled_panCK <- ScaleData(PANNTHR_Pt5_PreTreat_Visvader_merged_NORM_sampled_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
PANNTHR_Pt5_PreTreat_Visvader_merged_NORM_sampled_panCK <- RunPCA(PANNTHR_Pt5_PreTreat_Visvader_merged_NORM_sampled_panCK, features = VariableFeatures(object = PANNTHR_Pt5_PreTreat_Visvader_merged_NORM_sampled_panCK))

# Examine and visualize PCA results a few different ways
print(PANNTHR_Pt5_PreTreat_Visvader_merged_NORM_sampled_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(PANNTHR_Pt5_PreTreat_Visvader_merged_NORM_sampled_panCK, dims = 1:2, reduction = "pca")
DimPlot(PANNTHR_Pt5_PreTreat_Visvader_merged_NORM_sampled_panCK, reduction = "pca")

PANNTHR_Pt5_PreTreat_Visvader_merged_NORM_sampled_panCK <- FindNeighbors(PANNTHR_Pt5_PreTreat_Visvader_merged_NORM_sampled_panCK, dims = 1:20)
PANNTHR_Pt5_PreTreat_Visvader_merged_NORM_sampled_panCK <- FindClusters(PANNTHR_Pt5_PreTreat_Visvader_merged_NORM_sampled_panCK, resolution = 0.1)

# Umap clustering
PANNTHR_Pt5_PreTreat_Visvader_merged_NORM_sampled_panCK <- RunUMAP(PANNTHR_Pt5_PreTreat_Visvader_merged_NORM_sampled_panCK, dims = 1:20)
DimPlot(PANNTHR_Pt5_PreTreat_Visvader_merged_NORM_sampled_panCK, reduction = "umap")
DimPlot(PANNTHR_Pt5_PreTreat_Visvader_merged_NORM_sampled_panCK, reduction = "umap", group.by = "orig.ident")

# Save RDS
saveRDS(PANNTHR_Pt5_PreTreat_Visvader_merged_NORM_sampled_panCK, file = "/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Input/PANNTHR_inferCNV_Input_RDS/PANNTHR_Pt5_PreTreat_Visvader_merged_NORM_sampled_panCK.rds")

########## Prepare Cell Annotations
# Extract & Save Meta Data
PANNTHR_Pt5_PreTreat_Visvader_merged_NORM_sampled_panCK.meta.data <- as.data.frame(as.matrix(PANNTHR_Pt5_PreTreat_Visvader_merged_NORM_sampled_panCK@meta.data))
# Extract pertinent columns
# Make row names first column
PANNTHR_Pt5_PreTreat_Visvader_merged_NORM_sampled_panCK.meta.data <- tibble::rownames_to_column(PANNTHR_Pt5_PreTreat_Visvader_merged_NORM_sampled_panCK.meta.data, "Cells")
# In Seurat, "X" becomes the cell names while "cnv.ident" is specified above
PANNTHR_Pt5_PreTreat_Visvader_merged_NORM_sampled_panCK.annotations <- PANNTHR_Pt5_PreTreat_Visvader_merged_NORM_sampled_panCK.meta.data[, c("Cells", "cnv.ident")]
# Delete header
names(PANNTHR_Pt5_PreTreat_Visvader_merged_NORM_sampled_panCK.annotations) <- NULL
# Save as table
write.table(PANNTHR_Pt5_PreTreat_Visvader_merged_NORM_sampled_panCK.annotations,"/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Input/PANNTHR_Pt5_PreTreat_Visvader_merged_NORM_sampled_panCK.annotations.txt",sep="\t",row.names=FALSE)
# Remove metadata and annotations
rm(PANNTHR_Pt5_PreTreat_Visvader_merged_NORM_sampled_panCK.meta.data)
rm(PANNTHR_Pt5_PreTreat_Visvader_merged_NORM_sampled_panCK.annotations)

########## Generate Count Matrix
# Extract counts from the Seurat object
PANNTHR_Pt5_PreTreat_Visvader_merged_NORM_sampled_panCK.count <- PANNTHR_Pt5_PreTreat_Visvader_merged_NORM_sampled_panCK[["RNA"]]$counts
# Convert counts to a sparse matrix
PANNTHR_Pt5_PreTreat_Visvader_merged_NORM_sampled_panCK.count.matrix <- as(PANNTHR_Pt5_PreTreat_Visvader_merged_NORM_sampled_panCK.count, 'sparseMatrix')
rm(PANNTHR_Pt5_PreTreat_Visvader_merged_NORM_sampled_panCK.count)

# Convert to requisite data format
PANNTHR_Pt5_PreTreat_Visvader_merged_NORM_sampled_panCK.count.matrix <- as.data.frame(PANNTHR_Pt5_PreTreat_Visvader_merged_NORM_sampled_panCK.count.matrix)
# Save table
write.csv(PANNTHR_Pt5_PreTreat_Visvader_merged_NORM_sampled_panCK.count.matrix, file = "/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Input/PANNTHR_inferCNV_Input_Count_Matrix/PANNTHR_Pt5_PreTreat_Visvader_merged_NORM_sampled_panCK.count.matrix.csv")
# Remove Objects
rm(PANNTHR_Pt5_PreTreat_Visvader_merged_NORM_sampled_panCK)
rm(PANNTHR_Pt5_PreTreat_Visvader_merged_NORM_sampled_panCK.count.matrix)

########### Run inferCNV
# Note: check.names = FALSE prevents cell names from having dashes replaced with dots; essential to match annotations file
PANNTHR_Pt5_PreTreat_inferCNV <- CreateInfercnvObject(raw_counts_matrix=read.csv(file = "/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Input/PANNTHR_inferCNV_Input_Count_Matrix/PANNTHR_Pt5_PreTreat_Visvader_merged_NORM_sampled_panCK.count.matrix.csv", header = TRUE, row.names = 1,check.names = FALSE),
                                                  annotations_file=read.table(file = "/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Input/PANNTHR_Pt5_PreTreat_Visvader_merged_NORM_sampled_panCK.annotations.txt", header = FALSE, row.names = 1),
                                                  delim="\t",
                                                  gene_order_file="/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Input/hg38_gencode_v27.txt",
                                                  ref_group_names=c("Normal"))


# perform infercnv operations to reveal cnv signal
PANNTHR_Pt5_PreTreat_inferCNV = infercnv::run(PANNTHR_Pt5_PreTreat_inferCNV,
                                          cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                          out_dir="/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Output/PANNTHR_Pt5_PreTreat_panCK/",  # dir is auto-created for storing outputs
                                          cluster_by_groups=T,   # cluster
                                          denoise=T,
                                          HMM=T,
                                          num_ref_groups = 8, analysis_mode = "subclusters")

###########  Integrate data with Seurat object
# Read RDS and add inferCNV results to metadata
PANNTHR_Pt5_PreTreat_Visvader_merged_NORM_sampled_panCK <- readRDS(file = "/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Input/PANNTHR_inferCNV_Input_RDS/PANNTHR_Pt5_PreTreat_Visvader_merged_NORM_sampled_panCK.rds") 
PANNTHR_Pt5_PreTreat_panCK_cnv <- infercnv::add_to_seurat(infercnv_output_path="/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Output/PANNTHR_Pt5_PreTreat_panCK/",
                                                      seurat_obj=PANNTHR_Pt5_PreTreat_Visvader_merged_NORM_sampled_panCK, # optional
                                                      top_n=10)

# Collect Tumor cells from Seurat Object
PANNTHR_Pt5_PreTreat_panCK_cnv <- subset(x = PANNTHR_Pt5_PreTreat_panCK_cnv, subset = cnv.ident == "Tumor")
# Quantify chromosomes with alterations per cell
PANNTHR_Pt5_PreTreat_panCK_cnv_Pos_Chromosomes <- as.data.frame(rowSums(PANNTHR_Pt5_PreTreat_panCK_cnv@meta.data[,startsWith(names(PANNTHR_Pt5_PreTreat_panCK_cnv@meta.data),"has_cnv_")] > 0 ))
PANNTHR_Pt5_PreTreat_panCK_cnv_Neg_Chromosomes <- as.data.frame(rowSums(PANNTHR_Pt5_PreTreat_panCK_cnv@meta.data[,startsWith(names(PANNTHR_Pt5_PreTreat_panCK_cnv@meta.data),"has_cnv_")] == 0))
PANNTHR_Pt5_PreTreat_panCK_cnv_Chromosome_Quant <- as.data.frame(c(PANNTHR_Pt5_PreTreat_panCK_cnv_Pos_Chromosomes, PANNTHR_Pt5_PreTreat_panCK_cnv_Neg_Chromosomes))
colnames(PANNTHR_Pt5_PreTreat_panCK_cnv_Chromosome_Quant) <- c("Chromosomes CNV Pos", "Chromosomes CNV Neg")
write.csv(PANNTHR_Pt5_PreTreat_panCK_cnv_Chromosome_Quant, file = "/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Output/PANNTHR_Pt5_PreTreat_panCK/PANNTHR_Pt5_PreTreat_panCK_cnv_Chromosome_Quant.csv")

# Save CNV Metadata
PANNTHR_Pt5_PreTreat_panCK_cnv.meta.data <- as.data.frame(as.matrix(PANNTHR_Pt5_PreTreat_panCK_cnv@meta.data))
write.csv(PANNTHR_Pt5_PreTreat_panCK_cnv.meta.data, file = "/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Output/PANNTHR_Pt5_PreTreat_panCK/PANNTHR_Pt5_PreTreat_panCK_cnv.meta.data.csv")

# Subset Cells with inferred CNVs -------------------------------------------------------------------------------------------------------------
PANNTHR_Pt5_PreTreat_panCK_cnv <- subset(PANNTHR_Pt5_PreTreat_panCK_cnv, cells=rownames(PANNTHR_Pt5_PreTreat_panCK_cnv@meta.data)[rowSums(PANNTHR_Pt5_PreTreat_panCK_cnv@meta.data[,startsWith(names(PANNTHR_Pt5_PreTreat_panCK_cnv@meta.data),"has_cnv_")]) > 0 ])

# Save RDS with inferredCNVs
saveRDS(PANNTHR_Pt5_PreTreat_panCK_cnv, file = "/R/R_PANNTHR/PANNTHR_RDS/RDS_inferCNV/PANNTHR_Pt5_PreTreat_panCK_cnv.rds")


# Free memory
rm(PANNTHR_Pt5_PreTreat_inferCNV)
rm(PANNTHR_Pt5_PreTreat_Visvader_merged_NORM_sampled_panCK)
rm(PANNTHR_Pt5_PreTreat_panCK_cnv)
rm(PANNTHR_Pt5_PreTreat_panCK_cnv_Pos_Chromosomes)
rm(PANNTHR_Pt5_PreTreat_panCK_cnv_Neg_Chromosomes)
rm(PANNTHR_Pt5_PreTreat_panCK_cnv_Chromosome_Quant)
rm(PANNTHR_Pt5_PreTreat_panCK_cnv.meta.data)
gc()



#############################################
#### InferCNV: PANNTHR_Pt2_OnTreat_panCK ####
#############################################

########## Merge with normal mammary glad, extract epithelium, and process via Seurat
# Load RDS
PANNTHR_Pt2_OnTreat_panCK <- readRDS(file = "/R/R_PANNTHR/PANNTHR_RDS/RDS_panCK/PANNTHR_Pt2_OnTreat_panCK.rds")
Visvader_merged_NORM_sampled_panCK <- readRDS(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Visvader_merged_NORM_sampled_panCK.rds")

# Set identities 
colnames(x = PANNTHR_Pt2_OnTreat_panCK[[]])
Idents(object = PANNTHR_Pt2_OnTreat_panCK) <- 'Tumor'
PANNTHR_Pt2_OnTreat_panCK[["cnv.ident"]] <- Idents(object = PANNTHR_Pt2_OnTreat_panCK)
levels(x = PANNTHR_Pt2_OnTreat_panCK)

colnames(x = Visvader_merged_NORM_sampled_panCK[[]])
Idents(object = Visvader_merged_NORM_sampled_panCK) <- 'Normal'
Visvader_merged_NORM_sampled_panCK[["cnv.ident"]] <- Idents(object = Visvader_merged_NORM_sampled_panCK)
levels(x = Visvader_merged_NORM_sampled_panCK)

# Merge subsets into recombined RDS
PANNTHR_Pt2_OnTreat_Visvader_merged_NORM_sampled_panCK <- merge(x = PANNTHR_Pt2_OnTreat_panCK, y = Visvader_merged_NORM_sampled_panCK)
# Join Layers
PANNTHR_Pt2_OnTreat_Visvader_merged_NORM_sampled_panCK <- JoinLayers(PANNTHR_Pt2_OnTreat_Visvader_merged_NORM_sampled_panCK)

# Remove subsetted objects
rm(PANNTHR_Pt2_OnTreat_panCK)
rm(Visvader_merged_NORM_sampled_panCK)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
PANNTHR_Pt2_OnTreat_Visvader_merged_NORM_sampled_panCK <- FindVariableFeatures(PANNTHR_Pt2_OnTreat_Visvader_merged_NORM_sampled_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(PANNTHR_Pt2_OnTreat_Visvader_merged_NORM_sampled_panCK)
PANNTHR_Pt2_OnTreat_Visvader_merged_NORM_sampled_panCK <- ScaleData(PANNTHR_Pt2_OnTreat_Visvader_merged_NORM_sampled_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
PANNTHR_Pt2_OnTreat_Visvader_merged_NORM_sampled_panCK <- RunPCA(PANNTHR_Pt2_OnTreat_Visvader_merged_NORM_sampled_panCK, features = VariableFeatures(object = PANNTHR_Pt2_OnTreat_Visvader_merged_NORM_sampled_panCK))

# Examine and visualize PCA results a few different ways
print(PANNTHR_Pt2_OnTreat_Visvader_merged_NORM_sampled_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(PANNTHR_Pt2_OnTreat_Visvader_merged_NORM_sampled_panCK, dims = 1:2, reduction = "pca")
DimPlot(PANNTHR_Pt2_OnTreat_Visvader_merged_NORM_sampled_panCK, reduction = "pca")

PANNTHR_Pt2_OnTreat_Visvader_merged_NORM_sampled_panCK <- FindNeighbors(PANNTHR_Pt2_OnTreat_Visvader_merged_NORM_sampled_panCK, dims = 1:20)
PANNTHR_Pt2_OnTreat_Visvader_merged_NORM_sampled_panCK <- FindClusters(PANNTHR_Pt2_OnTreat_Visvader_merged_NORM_sampled_panCK, resolution = 0.1)

# Umap clustering
PANNTHR_Pt2_OnTreat_Visvader_merged_NORM_sampled_panCK <- RunUMAP(PANNTHR_Pt2_OnTreat_Visvader_merged_NORM_sampled_panCK, dims = 1:20)
DimPlot(PANNTHR_Pt2_OnTreat_Visvader_merged_NORM_sampled_panCK, reduction = "umap")
DimPlot(PANNTHR_Pt2_OnTreat_Visvader_merged_NORM_sampled_panCK, reduction = "umap", group.by = "orig.ident")

# Save RDS
saveRDS(PANNTHR_Pt2_OnTreat_Visvader_merged_NORM_sampled_panCK, file = "/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Input/PANNTHR_inferCNV_Input_RDS/PANNTHR_Pt2_OnTreat_Visvader_merged_NORM_sampled_panCK.rds")

########## Prepare Cell Annotations
# Extract & Save Meta Data
PANNTHR_Pt2_OnTreat_Visvader_merged_NORM_sampled_panCK.meta.data <- as.data.frame(as.matrix(PANNTHR_Pt2_OnTreat_Visvader_merged_NORM_sampled_panCK@meta.data))
# Extract pertinent columns
# Make row names first column
PANNTHR_Pt2_OnTreat_Visvader_merged_NORM_sampled_panCK.meta.data <- tibble::rownames_to_column(PANNTHR_Pt2_OnTreat_Visvader_merged_NORM_sampled_panCK.meta.data, "Cells")
# In Seurat, "X" becomes the cell names while "cnv.ident" is specified above
PANNTHR_Pt2_OnTreat_Visvader_merged_NORM_sampled_panCK.annotations <- PANNTHR_Pt2_OnTreat_Visvader_merged_NORM_sampled_panCK.meta.data[, c("Cells", "cnv.ident")]
# Delete header
names(PANNTHR_Pt2_OnTreat_Visvader_merged_NORM_sampled_panCK.annotations) <- NULL
# Save as table
write.table(PANNTHR_Pt2_OnTreat_Visvader_merged_NORM_sampled_panCK.annotations,"/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Input/PANNTHR_Pt2_OnTreat_Visvader_merged_NORM_sampled_panCK.annotations.txt",sep="\t",row.names=FALSE)
# Remove metadata and annotations
rm(PANNTHR_Pt2_OnTreat_Visvader_merged_NORM_sampled_panCK.meta.data)
rm(PANNTHR_Pt2_OnTreat_Visvader_merged_NORM_sampled_panCK.annotations)

########## Generate Count Matrix
# Extract counts from the Seurat object
PANNTHR_Pt2_OnTreat_Visvader_merged_NORM_sampled_panCK.count <- PANNTHR_Pt2_OnTreat_Visvader_merged_NORM_sampled_panCK[["RNA"]]$counts
# Convert counts to a sparse matrix
PANNTHR_Pt2_OnTreat_Visvader_merged_NORM_sampled_panCK.count.matrix <- as(PANNTHR_Pt2_OnTreat_Visvader_merged_NORM_sampled_panCK.count, 'sparseMatrix')
rm(PANNTHR_Pt2_OnTreat_Visvader_merged_NORM_sampled_panCK.count)

# Convert to requisite data format
PANNTHR_Pt2_OnTreat_Visvader_merged_NORM_sampled_panCK.count.matrix <- as.data.frame(PANNTHR_Pt2_OnTreat_Visvader_merged_NORM_sampled_panCK.count.matrix)
# Save table
write.csv(PANNTHR_Pt2_OnTreat_Visvader_merged_NORM_sampled_panCK.count.matrix, file = "/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Input/PANNTHR_inferCNV_Input_Count_Matrix/PANNTHR_Pt2_OnTreat_Visvader_merged_NORM_sampled_panCK.count.matrix.csv")
# Remove Objects
rm(PANNTHR_Pt2_OnTreat_Visvader_merged_NORM_sampled_panCK)
rm(PANNTHR_Pt2_OnTreat_Visvader_merged_NORM_sampled_panCK.count.matrix)

########### Run inferCNV
# Note: check.names = FALSE prevents cell names from having dashes replaced with dots; essential to match annotations file
PANNTHR_Pt2_OnTreat_inferCNV <- CreateInfercnvObject(raw_counts_matrix=read.csv(file = "/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Input/PANNTHR_inferCNV_Input_Count_Matrix/PANNTHR_Pt2_OnTreat_Visvader_merged_NORM_sampled_panCK.count.matrix.csv", header = TRUE, row.names = 1,check.names = FALSE),
                                                     annotations_file=read.table(file = "/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Input/PANNTHR_Pt2_OnTreat_Visvader_merged_NORM_sampled_panCK.annotations.txt", header = FALSE, row.names = 1),
                                                     delim="\t",
                                                     gene_order_file="/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Input/hg38_gencode_v27.txt",
                                                     ref_group_names=c("Normal"))


# perform infercnv operations to reveal cnv signal
PANNTHR_Pt2_OnTreat_inferCNV = infercnv::run(PANNTHR_Pt2_OnTreat_inferCNV,
                                             cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                             out_dir="/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Output/PANNTHR_Pt2_OnTreat_panCK/",  # dir is auto-created for storing outputs
                                             cluster_by_groups=T,   # cluster
                                             denoise=T,
                                             HMM=T,
                                             num_ref_groups = 8, analysis_mode = "subclusters")

###########  Integrate data with Seurat object
# Read RDS and add inferCNV results to metadata
PANNTHR_Pt2_OnTreat_Visvader_merged_NORM_sampled_panCK <- readRDS(file = "/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Input/PANNTHR_inferCNV_Input_RDS/PANNTHR_Pt2_OnTreat_Visvader_merged_NORM_sampled_panCK.rds") 
PANNTHR_Pt2_OnTreat_panCK_cnv <- infercnv::add_to_seurat(infercnv_output_path="/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Output/PANNTHR_Pt2_OnTreat_panCK/",
                                                         seurat_obj=PANNTHR_Pt2_OnTreat_Visvader_merged_NORM_sampled_panCK, # optional
                                                         top_n=10)

# Collect Tumor cells from Seurat Object
PANNTHR_Pt2_OnTreat_panCK_cnv <- subset(x = PANNTHR_Pt2_OnTreat_panCK_cnv, subset = cnv.ident == "Tumor")
# Quantify chromosomes with alterations per cell
PANNTHR_Pt2_OnTreat_panCK_cnv_Pos_Chromosomes <- as.data.frame(rowSums(PANNTHR_Pt2_OnTreat_panCK_cnv@meta.data[,startsWith(names(PANNTHR_Pt2_OnTreat_panCK_cnv@meta.data),"has_cnv_")] > 0 ))
PANNTHR_Pt2_OnTreat_panCK_cnv_Neg_Chromosomes <- as.data.frame(rowSums(PANNTHR_Pt2_OnTreat_panCK_cnv@meta.data[,startsWith(names(PANNTHR_Pt2_OnTreat_panCK_cnv@meta.data),"has_cnv_")] == 0))
PANNTHR_Pt2_OnTreat_panCK_cnv_Chromosome_Quant <- as.data.frame(c(PANNTHR_Pt2_OnTreat_panCK_cnv_Pos_Chromosomes, PANNTHR_Pt2_OnTreat_panCK_cnv_Neg_Chromosomes))
colnames(PANNTHR_Pt2_OnTreat_panCK_cnv_Chromosome_Quant) <- c("Chromosomes CNV Pos", "Chromosomes CNV Neg")
write.csv(PANNTHR_Pt2_OnTreat_panCK_cnv_Chromosome_Quant, file = "/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Output/PANNTHR_Pt2_OnTreat_panCK/PANNTHR_Pt2_OnTreat_panCK_cnv_Chromosome_Quant.csv")

# Save CNV Metadata
PANNTHR_Pt2_OnTreat_panCK_cnv.meta.data <- as.data.frame(as.matrix(PANNTHR_Pt2_OnTreat_panCK_cnv@meta.data))
write.csv(PANNTHR_Pt2_OnTreat_panCK_cnv.meta.data, file = "/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Output/PANNTHR_Pt2_OnTreat_panCK/PANNTHR_Pt2_OnTreat_panCK_cnv.meta.data.csv")

# Subset Cells with inferred CNVs -------------------------------------------------------------------------------------------------------------
PANNTHR_Pt2_OnTreat_panCK_cnv <- subset(PANNTHR_Pt2_OnTreat_panCK_cnv, cells=rownames(PANNTHR_Pt2_OnTreat_panCK_cnv@meta.data)[rowSums(PANNTHR_Pt2_OnTreat_panCK_cnv@meta.data[,startsWith(names(PANNTHR_Pt2_OnTreat_panCK_cnv@meta.data),"has_cnv_")]) > 0 ])

# Save RDS with inferredCNVs
saveRDS(PANNTHR_Pt2_OnTreat_panCK_cnv, file = "/R/R_PANNTHR/PANNTHR_RDS/RDS_inferCNV/PANNTHR_Pt2_OnTreat_panCK_cnv.rds")


# Free memory
rm(PANNTHR_Pt2_OnTreat_inferCNV)
rm(PANNTHR_Pt2_OnTreat_Visvader_merged_NORM_sampled_panCK)
rm(PANNTHR_Pt2_OnTreat_panCK_cnv)
rm(PANNTHR_Pt2_OnTreat_panCK_cnv_Pos_Chromosomes)
rm(PANNTHR_Pt2_OnTreat_panCK_cnv_Neg_Chromosomes)
rm(PANNTHR_Pt2_OnTreat_panCK_cnv_Chromosome_Quant)
rm(PANNTHR_Pt2_OnTreat_panCK_cnv.meta.data)
gc()



#############################################
#### InferCNV: PANNTHR_Pt3_OnTreat_panCK ####
#############################################

########## Merge with normal mammary glad, extract epithelium, and process via Seurat
# Load RDS
PANNTHR_Pt3_OnTreat_panCK <- readRDS(file = "/R/R_PANNTHR/PANNTHR_RDS/RDS_panCK/PANNTHR_Pt3_OnTreat_panCK.rds")
Visvader_merged_NORM_sampled_panCK <- readRDS(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Visvader_merged_NORM_sampled_panCK.rds")

# Set identities 
colnames(x = PANNTHR_Pt3_OnTreat_panCK[[]])
Idents(object = PANNTHR_Pt3_OnTreat_panCK) <- 'Tumor'
PANNTHR_Pt3_OnTreat_panCK[["cnv.ident"]] <- Idents(object = PANNTHR_Pt3_OnTreat_panCK)
levels(x = PANNTHR_Pt3_OnTreat_panCK)

colnames(x = Visvader_merged_NORM_sampled_panCK[[]])
Idents(object = Visvader_merged_NORM_sampled_panCK) <- 'Normal'
Visvader_merged_NORM_sampled_panCK[["cnv.ident"]] <- Idents(object = Visvader_merged_NORM_sampled_panCK)
levels(x = Visvader_merged_NORM_sampled_panCK)

# Merge subsets into recombined RDS
PANNTHR_Pt3_OnTreat_Visvader_merged_NORM_sampled_panCK <- merge(x = PANNTHR_Pt3_OnTreat_panCK, y = Visvader_merged_NORM_sampled_panCK)
# Join Layers
PANNTHR_Pt3_OnTreat_Visvader_merged_NORM_sampled_panCK <- JoinLayers(PANNTHR_Pt3_OnTreat_Visvader_merged_NORM_sampled_panCK)

# Remove subsetted objects
rm(PANNTHR_Pt3_OnTreat_panCK)
rm(Visvader_merged_NORM_sampled_panCK)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
PANNTHR_Pt3_OnTreat_Visvader_merged_NORM_sampled_panCK <- FindVariableFeatures(PANNTHR_Pt3_OnTreat_Visvader_merged_NORM_sampled_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(PANNTHR_Pt3_OnTreat_Visvader_merged_NORM_sampled_panCK)
PANNTHR_Pt3_OnTreat_Visvader_merged_NORM_sampled_panCK <- ScaleData(PANNTHR_Pt3_OnTreat_Visvader_merged_NORM_sampled_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
PANNTHR_Pt3_OnTreat_Visvader_merged_NORM_sampled_panCK <- RunPCA(PANNTHR_Pt3_OnTreat_Visvader_merged_NORM_sampled_panCK, features = VariableFeatures(object = PANNTHR_Pt3_OnTreat_Visvader_merged_NORM_sampled_panCK))

# Examine and visualize PCA results a few different ways
print(PANNTHR_Pt3_OnTreat_Visvader_merged_NORM_sampled_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(PANNTHR_Pt3_OnTreat_Visvader_merged_NORM_sampled_panCK, dims = 1:2, reduction = "pca")
DimPlot(PANNTHR_Pt3_OnTreat_Visvader_merged_NORM_sampled_panCK, reduction = "pca")

PANNTHR_Pt3_OnTreat_Visvader_merged_NORM_sampled_panCK <- FindNeighbors(PANNTHR_Pt3_OnTreat_Visvader_merged_NORM_sampled_panCK, dims = 1:20)
PANNTHR_Pt3_OnTreat_Visvader_merged_NORM_sampled_panCK <- FindClusters(PANNTHR_Pt3_OnTreat_Visvader_merged_NORM_sampled_panCK, resolution = 0.1)

# Umap clustering
PANNTHR_Pt3_OnTreat_Visvader_merged_NORM_sampled_panCK <- RunUMAP(PANNTHR_Pt3_OnTreat_Visvader_merged_NORM_sampled_panCK, dims = 1:20)
DimPlot(PANNTHR_Pt3_OnTreat_Visvader_merged_NORM_sampled_panCK, reduction = "umap")
DimPlot(PANNTHR_Pt3_OnTreat_Visvader_merged_NORM_sampled_panCK, reduction = "umap", group.by = "orig.ident")

# Save RDS
saveRDS(PANNTHR_Pt3_OnTreat_Visvader_merged_NORM_sampled_panCK, file = "/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Input/PANNTHR_inferCNV_Input_RDS/PANNTHR_Pt3_OnTreat_Visvader_merged_NORM_sampled_panCK.rds")

########## Prepare Cell Annotations
# Extract & Save Meta Data
PANNTHR_Pt3_OnTreat_Visvader_merged_NORM_sampled_panCK.meta.data <- as.data.frame(as.matrix(PANNTHR_Pt3_OnTreat_Visvader_merged_NORM_sampled_panCK@meta.data))
# Extract pertinent columns
# Make row names first column
PANNTHR_Pt3_OnTreat_Visvader_merged_NORM_sampled_panCK.meta.data <- tibble::rownames_to_column(PANNTHR_Pt3_OnTreat_Visvader_merged_NORM_sampled_panCK.meta.data, "Cells")
# In Seurat, "X" becomes the cell names while "cnv.ident" is specified above
PANNTHR_Pt3_OnTreat_Visvader_merged_NORM_sampled_panCK.annotations <- PANNTHR_Pt3_OnTreat_Visvader_merged_NORM_sampled_panCK.meta.data[, c("Cells", "cnv.ident")]
# Delete header
names(PANNTHR_Pt3_OnTreat_Visvader_merged_NORM_sampled_panCK.annotations) <- NULL
# Save as table
write.table(PANNTHR_Pt3_OnTreat_Visvader_merged_NORM_sampled_panCK.annotations,"/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Input/PANNTHR_Pt3_OnTreat_Visvader_merged_NORM_sampled_panCK.annotations.txt",sep="\t",row.names=FALSE)
# Remove metadata and annotations
rm(PANNTHR_Pt3_OnTreat_Visvader_merged_NORM_sampled_panCK.meta.data)
rm(PANNTHR_Pt3_OnTreat_Visvader_merged_NORM_sampled_panCK.annotations)

########## Generate Count Matrix
# Extract counts from the Seurat object
PANNTHR_Pt3_OnTreat_Visvader_merged_NORM_sampled_panCK.count <- PANNTHR_Pt3_OnTreat_Visvader_merged_NORM_sampled_panCK[["RNA"]]$counts
# Convert counts to a sparse matrix
PANNTHR_Pt3_OnTreat_Visvader_merged_NORM_sampled_panCK.count.matrix <- as(PANNTHR_Pt3_OnTreat_Visvader_merged_NORM_sampled_panCK.count, 'sparseMatrix')
rm(PANNTHR_Pt3_OnTreat_Visvader_merged_NORM_sampled_panCK.count)

# Convert to requisite data format
PANNTHR_Pt3_OnTreat_Visvader_merged_NORM_sampled_panCK.count.matrix <- as.data.frame(PANNTHR_Pt3_OnTreat_Visvader_merged_NORM_sampled_panCK.count.matrix)
# Save table
write.csv(PANNTHR_Pt3_OnTreat_Visvader_merged_NORM_sampled_panCK.count.matrix, file = "/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Input/PANNTHR_inferCNV_Input_Count_Matrix/PANNTHR_Pt3_OnTreat_Visvader_merged_NORM_sampled_panCK.count.matrix.csv")
# Remove Objects
rm(PANNTHR_Pt3_OnTreat_Visvader_merged_NORM_sampled_panCK)
rm(PANNTHR_Pt3_OnTreat_Visvader_merged_NORM_sampled_panCK.count.matrix)

########### Run inferCNV
# Note: check.names = FALSE prevents cell names from having dashes replaced with dots; essential to match annotations file
PANNTHR_Pt3_OnTreat_inferCNV <- CreateInfercnvObject(raw_counts_matrix=read.csv(file = "/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Input/PANNTHR_inferCNV_Input_Count_Matrix/PANNTHR_Pt3_OnTreat_Visvader_merged_NORM_sampled_panCK.count.matrix.csv", header = TRUE, row.names = 1,check.names = FALSE),
                                                     annotations_file=read.table(file = "/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Input/PANNTHR_Pt3_OnTreat_Visvader_merged_NORM_sampled_panCK.annotations.txt", header = FALSE, row.names = 1),
                                                     delim="\t",
                                                     gene_order_file="/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Input/hg38_gencode_v27.txt",
                                                     ref_group_names=c("Normal"))


# perform infercnv operations to reveal cnv signal
PANNTHR_Pt3_OnTreat_inferCNV = infercnv::run(PANNTHR_Pt3_OnTreat_inferCNV,
                                             cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                             out_dir="/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Output/PANNTHR_Pt3_OnTreat_panCK/",  # dir is auto-created for storing outputs
                                             cluster_by_groups=T,   # cluster
                                             denoise=T,
                                             HMM=T,
                                             num_ref_groups = 8, analysis_mode = "subclusters")

###########  Integrate data with Seurat object
# Read RDS and add inferCNV results to metadata
PANNTHR_Pt3_OnTreat_Visvader_merged_NORM_sampled_panCK <- readRDS(file = "/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Input/PANNTHR_inferCNV_Input_RDS/PANNTHR_Pt3_OnTreat_Visvader_merged_NORM_sampled_panCK.rds") 
PANNTHR_Pt3_OnTreat_panCK_cnv <- infercnv::add_to_seurat(infercnv_output_path="/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Output/PANNTHR_Pt3_OnTreat_panCK/",
                                                         seurat_obj=PANNTHR_Pt3_OnTreat_Visvader_merged_NORM_sampled_panCK, # optional
                                                         top_n=10)

# Collect Tumor cells from Seurat Object
PANNTHR_Pt3_OnTreat_panCK_cnv <- subset(x = PANNTHR_Pt3_OnTreat_panCK_cnv, subset = cnv.ident == "Tumor")
# Quantify chromosomes with alterations per cell
PANNTHR_Pt3_OnTreat_panCK_cnv_Pos_Chromosomes <- as.data.frame(rowSums(PANNTHR_Pt3_OnTreat_panCK_cnv@meta.data[,startsWith(names(PANNTHR_Pt3_OnTreat_panCK_cnv@meta.data),"has_cnv_")] > 0 ))
PANNTHR_Pt3_OnTreat_panCK_cnv_Neg_Chromosomes <- as.data.frame(rowSums(PANNTHR_Pt3_OnTreat_panCK_cnv@meta.data[,startsWith(names(PANNTHR_Pt3_OnTreat_panCK_cnv@meta.data),"has_cnv_")] == 0))
PANNTHR_Pt3_OnTreat_panCK_cnv_Chromosome_Quant <- as.data.frame(c(PANNTHR_Pt3_OnTreat_panCK_cnv_Pos_Chromosomes, PANNTHR_Pt3_OnTreat_panCK_cnv_Neg_Chromosomes))
colnames(PANNTHR_Pt3_OnTreat_panCK_cnv_Chromosome_Quant) <- c("Chromosomes CNV Pos", "Chromosomes CNV Neg")
write.csv(PANNTHR_Pt3_OnTreat_panCK_cnv_Chromosome_Quant, file = "/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Output/PANNTHR_Pt3_OnTreat_panCK/PANNTHR_Pt3_OnTreat_panCK_cnv_Chromosome_Quant.csv")

# Save CNV Metadata
PANNTHR_Pt3_OnTreat_panCK_cnv.meta.data <- as.data.frame(as.matrix(PANNTHR_Pt3_OnTreat_panCK_cnv@meta.data))
write.csv(PANNTHR_Pt3_OnTreat_panCK_cnv.meta.data, file = "/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Output/PANNTHR_Pt3_OnTreat_panCK/PANNTHR_Pt3_OnTreat_panCK_cnv.meta.data.csv")

# Subset Cells with inferred CNVs -------------------------------------------------------------------------------------------------------------
PANNTHR_Pt3_OnTreat_panCK_cnv <- subset(PANNTHR_Pt3_OnTreat_panCK_cnv, cells=rownames(PANNTHR_Pt3_OnTreat_panCK_cnv@meta.data)[rowSums(PANNTHR_Pt3_OnTreat_panCK_cnv@meta.data[,startsWith(names(PANNTHR_Pt3_OnTreat_panCK_cnv@meta.data),"has_cnv_")]) > 0 ])

# Save RDS with inferredCNVs
saveRDS(PANNTHR_Pt3_OnTreat_panCK_cnv, file = "/R/R_PANNTHR/PANNTHR_RDS/RDS_inferCNV/PANNTHR_Pt3_OnTreat_panCK_cnv.rds")


# Free memory
rm(PANNTHR_Pt3_OnTreat_inferCNV)
rm(PANNTHR_Pt3_OnTreat_Visvader_merged_NORM_sampled_panCK)
rm(PANNTHR_Pt3_OnTreat_panCK_cnv)
rm(PANNTHR_Pt3_OnTreat_panCK_cnv_Pos_Chromosomes)
rm(PANNTHR_Pt3_OnTreat_panCK_cnv_Neg_Chromosomes)
rm(PANNTHR_Pt3_OnTreat_panCK_cnv_Chromosome_Quant)
rm(PANNTHR_Pt3_OnTreat_panCK_cnv.meta.data)
gc()



#############################################
#### InferCNV: PANNTHR_Pt4_OnTreat_panCK ####
#############################################

########## Merge with normal mammary glad, extract epithelium, and process via Seurat
# Load RDS
PANNTHR_Pt4_OnTreat_panCK <- readRDS(file = "/R/R_PANNTHR/PANNTHR_RDS/RDS_panCK/PANNTHR_Pt4_OnTreat_panCK.rds")
Visvader_merged_NORM_sampled_panCK <- readRDS(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Visvader_merged_NORM_sampled_panCK.rds")

# Set identities 
colnames(x = PANNTHR_Pt4_OnTreat_panCK[[]])
Idents(object = PANNTHR_Pt4_OnTreat_panCK) <- 'Tumor'
PANNTHR_Pt4_OnTreat_panCK[["cnv.ident"]] <- Idents(object = PANNTHR_Pt4_OnTreat_panCK)
levels(x = PANNTHR_Pt4_OnTreat_panCK)

colnames(x = Visvader_merged_NORM_sampled_panCK[[]])
Idents(object = Visvader_merged_NORM_sampled_panCK) <- 'Normal'
Visvader_merged_NORM_sampled_panCK[["cnv.ident"]] <- Idents(object = Visvader_merged_NORM_sampled_panCK)
levels(x = Visvader_merged_NORM_sampled_panCK)

# Merge subsets into recombined RDS
PANNTHR_Pt4_OnTreat_Visvader_merged_NORM_sampled_panCK <- merge(x = PANNTHR_Pt4_OnTreat_panCK, y = Visvader_merged_NORM_sampled_panCK)
# Join Layers
PANNTHR_Pt4_OnTreat_Visvader_merged_NORM_sampled_panCK <- JoinLayers(PANNTHR_Pt4_OnTreat_Visvader_merged_NORM_sampled_panCK)

# Remove subsetted objects
rm(PANNTHR_Pt4_OnTreat_panCK)
rm(Visvader_merged_NORM_sampled_panCK)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
PANNTHR_Pt4_OnTreat_Visvader_merged_NORM_sampled_panCK <- FindVariableFeatures(PANNTHR_Pt4_OnTreat_Visvader_merged_NORM_sampled_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(PANNTHR_Pt4_OnTreat_Visvader_merged_NORM_sampled_panCK)
PANNTHR_Pt4_OnTreat_Visvader_merged_NORM_sampled_panCK <- ScaleData(PANNTHR_Pt4_OnTreat_Visvader_merged_NORM_sampled_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
PANNTHR_Pt4_OnTreat_Visvader_merged_NORM_sampled_panCK <- RunPCA(PANNTHR_Pt4_OnTreat_Visvader_merged_NORM_sampled_panCK, features = VariableFeatures(object = PANNTHR_Pt4_OnTreat_Visvader_merged_NORM_sampled_panCK))

# Examine and visualize PCA results a few different ways
print(PANNTHR_Pt4_OnTreat_Visvader_merged_NORM_sampled_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(PANNTHR_Pt4_OnTreat_Visvader_merged_NORM_sampled_panCK, dims = 1:2, reduction = "pca")
DimPlot(PANNTHR_Pt4_OnTreat_Visvader_merged_NORM_sampled_panCK, reduction = "pca")

PANNTHR_Pt4_OnTreat_Visvader_merged_NORM_sampled_panCK <- FindNeighbors(PANNTHR_Pt4_OnTreat_Visvader_merged_NORM_sampled_panCK, dims = 1:20)
PANNTHR_Pt4_OnTreat_Visvader_merged_NORM_sampled_panCK <- FindClusters(PANNTHR_Pt4_OnTreat_Visvader_merged_NORM_sampled_panCK, resolution = 0.1)

# Umap clustering
PANNTHR_Pt4_OnTreat_Visvader_merged_NORM_sampled_panCK <- RunUMAP(PANNTHR_Pt4_OnTreat_Visvader_merged_NORM_sampled_panCK, dims = 1:20)
DimPlot(PANNTHR_Pt4_OnTreat_Visvader_merged_NORM_sampled_panCK, reduction = "umap")
DimPlot(PANNTHR_Pt4_OnTreat_Visvader_merged_NORM_sampled_panCK, reduction = "umap", group.by = "orig.ident")

# Save RDS
saveRDS(PANNTHR_Pt4_OnTreat_Visvader_merged_NORM_sampled_panCK, file = "/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Input/PANNTHR_inferCNV_Input_RDS/PANNTHR_Pt4_OnTreat_Visvader_merged_NORM_sampled_panCK.rds")

########## Prepare Cell Annotations
# Extract & Save Meta Data
PANNTHR_Pt4_OnTreat_Visvader_merged_NORM_sampled_panCK.meta.data <- as.data.frame(as.matrix(PANNTHR_Pt4_OnTreat_Visvader_merged_NORM_sampled_panCK@meta.data))
# Extract pertinent columns
# Make row names first column
PANNTHR_Pt4_OnTreat_Visvader_merged_NORM_sampled_panCK.meta.data <- tibble::rownames_to_column(PANNTHR_Pt4_OnTreat_Visvader_merged_NORM_sampled_panCK.meta.data, "Cells")
# In Seurat, "X" becomes the cell names while "cnv.ident" is specified above
PANNTHR_Pt4_OnTreat_Visvader_merged_NORM_sampled_panCK.annotations <- PANNTHR_Pt4_OnTreat_Visvader_merged_NORM_sampled_panCK.meta.data[, c("Cells", "cnv.ident")]
# Delete header
names(PANNTHR_Pt4_OnTreat_Visvader_merged_NORM_sampled_panCK.annotations) <- NULL
# Save as table
write.table(PANNTHR_Pt4_OnTreat_Visvader_merged_NORM_sampled_panCK.annotations,"/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Input/PANNTHR_Pt4_OnTreat_Visvader_merged_NORM_sampled_panCK.annotations.txt",sep="\t",row.names=FALSE)
# Remove metadata and annotations
rm(PANNTHR_Pt4_OnTreat_Visvader_merged_NORM_sampled_panCK.meta.data)
rm(PANNTHR_Pt4_OnTreat_Visvader_merged_NORM_sampled_panCK.annotations)

########## Generate Count Matrix
# Extract counts from the Seurat object
PANNTHR_Pt4_OnTreat_Visvader_merged_NORM_sampled_panCK.count <- PANNTHR_Pt4_OnTreat_Visvader_merged_NORM_sampled_panCK[["RNA"]]$counts
# Convert counts to a sparse matrix
PANNTHR_Pt4_OnTreat_Visvader_merged_NORM_sampled_panCK.count.matrix <- as(PANNTHR_Pt4_OnTreat_Visvader_merged_NORM_sampled_panCK.count, 'sparseMatrix')
rm(PANNTHR_Pt4_OnTreat_Visvader_merged_NORM_sampled_panCK.count)

# Convert to requisite data format
PANNTHR_Pt4_OnTreat_Visvader_merged_NORM_sampled_panCK.count.matrix <- as.data.frame(PANNTHR_Pt4_OnTreat_Visvader_merged_NORM_sampled_panCK.count.matrix)
# Save table
write.csv(PANNTHR_Pt4_OnTreat_Visvader_merged_NORM_sampled_panCK.count.matrix, file = "/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Input/PANNTHR_inferCNV_Input_Count_Matrix/PANNTHR_Pt4_OnTreat_Visvader_merged_NORM_sampled_panCK.count.matrix.csv")
# Remove Objects
rm(PANNTHR_Pt4_OnTreat_Visvader_merged_NORM_sampled_panCK)
rm(PANNTHR_Pt4_OnTreat_Visvader_merged_NORM_sampled_panCK.count.matrix)

########### Run inferCNV
# Note: check.names = FALSE prevents cell names from having dashes replaced with dots; essential to match annotations file
PANNTHR_Pt4_OnTreat_inferCNV <- CreateInfercnvObject(raw_counts_matrix=read.csv(file = "/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Input/PANNTHR_inferCNV_Input_Count_Matrix/PANNTHR_Pt4_OnTreat_Visvader_merged_NORM_sampled_panCK.count.matrix.csv", header = TRUE, row.names = 1,check.names = FALSE),
                                                     annotations_file=read.table(file = "/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Input/PANNTHR_Pt4_OnTreat_Visvader_merged_NORM_sampled_panCK.annotations.txt", header = FALSE, row.names = 1),
                                                     delim="\t",
                                                     gene_order_file="/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Input/hg38_gencode_v27.txt",
                                                     ref_group_names=c("Normal"))


# perform infercnv operations to reveal cnv signal
PANNTHR_Pt4_OnTreat_inferCNV = infercnv::run(PANNTHR_Pt4_OnTreat_inferCNV,
                                             cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                             out_dir="/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Output/PANNTHR_Pt4_OnTreat_panCK/",  # dir is auto-created for storing outputs
                                             cluster_by_groups=T,   # cluster
                                             denoise=T,
                                             HMM=T,
                                             num_ref_groups = 8, analysis_mode = "subclusters")

###########  Integrate data with Seurat object
# Read RDS and add inferCNV results to metadata
PANNTHR_Pt4_OnTreat_Visvader_merged_NORM_sampled_panCK <- readRDS(file = "/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Input/PANNTHR_inferCNV_Input_RDS/PANNTHR_Pt4_OnTreat_Visvader_merged_NORM_sampled_panCK.rds") 
PANNTHR_Pt4_OnTreat_panCK_cnv <- infercnv::add_to_seurat(infercnv_output_path="/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Output/PANNTHR_Pt4_OnTreat_panCK/",
                                                         seurat_obj=PANNTHR_Pt4_OnTreat_Visvader_merged_NORM_sampled_panCK, # optional
                                                         top_n=10)

# Collect Tumor cells from Seurat Object
PANNTHR_Pt4_OnTreat_panCK_cnv <- subset(x = PANNTHR_Pt4_OnTreat_panCK_cnv, subset = cnv.ident == "Tumor")
# Quantify chromosomes with alterations per cell
PANNTHR_Pt4_OnTreat_panCK_cnv_Pos_Chromosomes <- as.data.frame(rowSums(PANNTHR_Pt4_OnTreat_panCK_cnv@meta.data[,startsWith(names(PANNTHR_Pt4_OnTreat_panCK_cnv@meta.data),"has_cnv_")] > 0 ))
PANNTHR_Pt4_OnTreat_panCK_cnv_Neg_Chromosomes <- as.data.frame(rowSums(PANNTHR_Pt4_OnTreat_panCK_cnv@meta.data[,startsWith(names(PANNTHR_Pt4_OnTreat_panCK_cnv@meta.data),"has_cnv_")] == 0))
PANNTHR_Pt4_OnTreat_panCK_cnv_Chromosome_Quant <- as.data.frame(c(PANNTHR_Pt4_OnTreat_panCK_cnv_Pos_Chromosomes, PANNTHR_Pt4_OnTreat_panCK_cnv_Neg_Chromosomes))
colnames(PANNTHR_Pt4_OnTreat_panCK_cnv_Chromosome_Quant) <- c("Chromosomes CNV Pos", "Chromosomes CNV Neg")
write.csv(PANNTHR_Pt4_OnTreat_panCK_cnv_Chromosome_Quant, file = "/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Output/PANNTHR_Pt4_OnTreat_panCK/PANNTHR_Pt4_OnTreat_panCK_cnv_Chromosome_Quant.csv")

# Save CNV Metadata
PANNTHR_Pt4_OnTreat_panCK_cnv.meta.data <- as.data.frame(as.matrix(PANNTHR_Pt4_OnTreat_panCK_cnv@meta.data))
write.csv(PANNTHR_Pt4_OnTreat_panCK_cnv.meta.data, file = "/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Output/PANNTHR_Pt4_OnTreat_panCK/PANNTHR_Pt4_OnTreat_panCK_cnv.meta.data.csv")

# Subset Cells with inferred CNVs -------------------------------------------------------------------------------------------------------------
PANNTHR_Pt4_OnTreat_panCK_cnv <- subset(PANNTHR_Pt4_OnTreat_panCK_cnv, cells=rownames(PANNTHR_Pt4_OnTreat_panCK_cnv@meta.data)[rowSums(PANNTHR_Pt4_OnTreat_panCK_cnv@meta.data[,startsWith(names(PANNTHR_Pt4_OnTreat_panCK_cnv@meta.data),"has_cnv_")]) > 0 ])

# Save RDS with inferredCNVs
saveRDS(PANNTHR_Pt4_OnTreat_panCK_cnv, file = "/R/R_PANNTHR/PANNTHR_RDS/RDS_inferCNV/PANNTHR_Pt4_OnTreat_panCK_cnv.rds")


# Free memory
rm(PANNTHR_Pt4_OnTreat_inferCNV)
rm(PANNTHR_Pt4_OnTreat_Visvader_merged_NORM_sampled_panCK)
rm(PANNTHR_Pt4_OnTreat_panCK_cnv)
rm(PANNTHR_Pt4_OnTreat_panCK_cnv_Pos_Chromosomes)
rm(PANNTHR_Pt4_OnTreat_panCK_cnv_Neg_Chromosomes)
rm(PANNTHR_Pt4_OnTreat_panCK_cnv_Chromosome_Quant)
rm(PANNTHR_Pt4_OnTreat_panCK_cnv.meta.data)
gc()



#############################################
#### InferCNV: PANNTHR_Pt6_OnTreat_panCK ####
#############################################

########## Merge with normal mammary glad, extract epithelium, and process via Seurat
# Load RDS
PANNTHR_Pt6_OnTreat_panCK <- readRDS(file = "/R/R_PANNTHR/PANNTHR_RDS/RDS_panCK/PANNTHR_Pt6_OnTreat_panCK.rds")
Visvader_merged_NORM_sampled_panCK <- readRDS(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Visvader_merged_NORM_sampled_panCK.rds")

# Set identities 
colnames(x = PANNTHR_Pt6_OnTreat_panCK[[]])
Idents(object = PANNTHR_Pt6_OnTreat_panCK) <- 'Tumor'
PANNTHR_Pt6_OnTreat_panCK[["cnv.ident"]] <- Idents(object = PANNTHR_Pt6_OnTreat_panCK)
levels(x = PANNTHR_Pt6_OnTreat_panCK)

colnames(x = Visvader_merged_NORM_sampled_panCK[[]])
Idents(object = Visvader_merged_NORM_sampled_panCK) <- 'Normal'
Visvader_merged_NORM_sampled_panCK[["cnv.ident"]] <- Idents(object = Visvader_merged_NORM_sampled_panCK)
levels(x = Visvader_merged_NORM_sampled_panCK)

# Merge subsets into recombined RDS
PANNTHR_Pt6_OnTreat_Visvader_merged_NORM_sampled_panCK <- merge(x = PANNTHR_Pt6_OnTreat_panCK, y = Visvader_merged_NORM_sampled_panCK)
# Join Layers
PANNTHR_Pt6_OnTreat_Visvader_merged_NORM_sampled_panCK <- JoinLayers(PANNTHR_Pt6_OnTreat_Visvader_merged_NORM_sampled_panCK)

# Remove subsetted objects
rm(PANNTHR_Pt6_OnTreat_panCK)
rm(Visvader_merged_NORM_sampled_panCK)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
PANNTHR_Pt6_OnTreat_Visvader_merged_NORM_sampled_panCK <- FindVariableFeatures(PANNTHR_Pt6_OnTreat_Visvader_merged_NORM_sampled_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(PANNTHR_Pt6_OnTreat_Visvader_merged_NORM_sampled_panCK)
PANNTHR_Pt6_OnTreat_Visvader_merged_NORM_sampled_panCK <- ScaleData(PANNTHR_Pt6_OnTreat_Visvader_merged_NORM_sampled_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
PANNTHR_Pt6_OnTreat_Visvader_merged_NORM_sampled_panCK <- RunPCA(PANNTHR_Pt6_OnTreat_Visvader_merged_NORM_sampled_panCK, features = VariableFeatures(object = PANNTHR_Pt6_OnTreat_Visvader_merged_NORM_sampled_panCK))

# Examine and visualize PCA results a few different ways
print(PANNTHR_Pt6_OnTreat_Visvader_merged_NORM_sampled_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(PANNTHR_Pt6_OnTreat_Visvader_merged_NORM_sampled_panCK, dims = 1:2, reduction = "pca")
DimPlot(PANNTHR_Pt6_OnTreat_Visvader_merged_NORM_sampled_panCK, reduction = "pca")

PANNTHR_Pt6_OnTreat_Visvader_merged_NORM_sampled_panCK <- FindNeighbors(PANNTHR_Pt6_OnTreat_Visvader_merged_NORM_sampled_panCK, dims = 1:20)
PANNTHR_Pt6_OnTreat_Visvader_merged_NORM_sampled_panCK <- FindClusters(PANNTHR_Pt6_OnTreat_Visvader_merged_NORM_sampled_panCK, resolution = 0.1)

# Umap clustering
PANNTHR_Pt6_OnTreat_Visvader_merged_NORM_sampled_panCK <- RunUMAP(PANNTHR_Pt6_OnTreat_Visvader_merged_NORM_sampled_panCK, dims = 1:20)
DimPlot(PANNTHR_Pt6_OnTreat_Visvader_merged_NORM_sampled_panCK, reduction = "umap")
DimPlot(PANNTHR_Pt6_OnTreat_Visvader_merged_NORM_sampled_panCK, reduction = "umap", group.by = "orig.ident")

# Save RDS
saveRDS(PANNTHR_Pt6_OnTreat_Visvader_merged_NORM_sampled_panCK, file = "/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Input/PANNTHR_inferCNV_Input_RDS/PANNTHR_Pt6_OnTreat_Visvader_merged_NORM_sampled_panCK.rds")

########## Prepare Cell Annotations
# Extract & Save Meta Data
PANNTHR_Pt6_OnTreat_Visvader_merged_NORM_sampled_panCK.meta.data <- as.data.frame(as.matrix(PANNTHR_Pt6_OnTreat_Visvader_merged_NORM_sampled_panCK@meta.data))
# Extract pertinent columns
# Make row names first column
PANNTHR_Pt6_OnTreat_Visvader_merged_NORM_sampled_panCK.meta.data <- tibble::rownames_to_column(PANNTHR_Pt6_OnTreat_Visvader_merged_NORM_sampled_panCK.meta.data, "Cells")
# In Seurat, "X" becomes the cell names while "cnv.ident" is specified above
PANNTHR_Pt6_OnTreat_Visvader_merged_NORM_sampled_panCK.annotations <- PANNTHR_Pt6_OnTreat_Visvader_merged_NORM_sampled_panCK.meta.data[, c("Cells", "cnv.ident")]
# Delete header
names(PANNTHR_Pt6_OnTreat_Visvader_merged_NORM_sampled_panCK.annotations) <- NULL
# Save as table
write.table(PANNTHR_Pt6_OnTreat_Visvader_merged_NORM_sampled_panCK.annotations,"/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Input/PANNTHR_Pt6_OnTreat_Visvader_merged_NORM_sampled_panCK.annotations.txt",sep="\t",row.names=FALSE)
# Remove metadata and annotations
rm(PANNTHR_Pt6_OnTreat_Visvader_merged_NORM_sampled_panCK.meta.data)
rm(PANNTHR_Pt6_OnTreat_Visvader_merged_NORM_sampled_panCK.annotations)

########## Generate Count Matrix
# Extract counts from the Seurat object
PANNTHR_Pt6_OnTreat_Visvader_merged_NORM_sampled_panCK.count <- PANNTHR_Pt6_OnTreat_Visvader_merged_NORM_sampled_panCK[["RNA"]]$counts
# Convert counts to a sparse matrix
PANNTHR_Pt6_OnTreat_Visvader_merged_NORM_sampled_panCK.count.matrix <- as(PANNTHR_Pt6_OnTreat_Visvader_merged_NORM_sampled_panCK.count, 'sparseMatrix')
rm(PANNTHR_Pt6_OnTreat_Visvader_merged_NORM_sampled_panCK.count)

# Convert to requisite data format
PANNTHR_Pt6_OnTreat_Visvader_merged_NORM_sampled_panCK.count.matrix <- as.data.frame(PANNTHR_Pt6_OnTreat_Visvader_merged_NORM_sampled_panCK.count.matrix)
# Save table
write.csv(PANNTHR_Pt6_OnTreat_Visvader_merged_NORM_sampled_panCK.count.matrix, file = "/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Input/PANNTHR_inferCNV_Input_Count_Matrix/PANNTHR_Pt6_OnTreat_Visvader_merged_NORM_sampled_panCK.count.matrix.csv")
# Remove Objects
rm(PANNTHR_Pt6_OnTreat_Visvader_merged_NORM_sampled_panCK)
rm(PANNTHR_Pt6_OnTreat_Visvader_merged_NORM_sampled_panCK.count.matrix)

########### Run inferCNV
# Note: check.names = FALSE prevents cell names from having dashes replaced with dots; essential to match annotations file
PANNTHR_Pt6_OnTreat_inferCNV <- CreateInfercnvObject(raw_counts_matrix=read.csv(file = "/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Input/PANNTHR_inferCNV_Input_Count_Matrix/PANNTHR_Pt6_OnTreat_Visvader_merged_NORM_sampled_panCK.count.matrix.csv", header = TRUE, row.names = 1,check.names = FALSE),
                                               annotations_file=read.table(file = "/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Input/PANNTHR_Pt6_OnTreat_Visvader_merged_NORM_sampled_panCK.annotations.txt", header = FALSE, row.names = 1),
                                               delim="\t",
                                               gene_order_file="/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Input/hg38_gencode_v27.txt",
                                               ref_group_names=c("Normal"))


# perform infercnv operations to reveal cnv signal
PANNTHR_Pt6_OnTreat_inferCNV = infercnv::run(PANNTHR_Pt6_OnTreat_inferCNV,
                                       cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                       out_dir="/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Output/PANNTHR_Pt6_OnTreat_panCK/",  # dir is auto-created for storing outputs
                                       cluster_by_groups=T,   # cluster
                                       denoise=T,
                                       HMM=T,
                                       num_ref_groups = 8, analysis_mode = "subclusters")

###########  Integrate data with Seurat object
# Read RDS and add inferCNV results to metadata
PANNTHR_Pt6_OnTreat_Visvader_merged_NORM_sampled_panCK <- readRDS(file = "/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Input/PANNTHR_inferCNV_Input_RDS/PANNTHR_Pt6_OnTreat_Visvader_merged_NORM_sampled_panCK.rds") 
PANNTHR_Pt6_OnTreat_panCK_cnv <- infercnv::add_to_seurat(infercnv_output_path="/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Output/PANNTHR_Pt6_OnTreat_panCK/",
                                                   seurat_obj=PANNTHR_Pt6_OnTreat_Visvader_merged_NORM_sampled_panCK, # optional
                                                   top_n=10)

# Collect Tumor cells from Seurat Object
PANNTHR_Pt6_OnTreat_panCK_cnv <- subset(x = PANNTHR_Pt6_OnTreat_panCK_cnv, subset = cnv.ident == "Tumor")
# Quantify chromosomes with alterations per cell
PANNTHR_Pt6_OnTreat_panCK_cnv_Pos_Chromosomes <- as.data.frame(rowSums(PANNTHR_Pt6_OnTreat_panCK_cnv@meta.data[,startsWith(names(PANNTHR_Pt6_OnTreat_panCK_cnv@meta.data),"has_cnv_")] > 0 ))
PANNTHR_Pt6_OnTreat_panCK_cnv_Neg_Chromosomes <- as.data.frame(rowSums(PANNTHR_Pt6_OnTreat_panCK_cnv@meta.data[,startsWith(names(PANNTHR_Pt6_OnTreat_panCK_cnv@meta.data),"has_cnv_")] == 0))
PANNTHR_Pt6_OnTreat_panCK_cnv_Chromosome_Quant <- as.data.frame(c(PANNTHR_Pt6_OnTreat_panCK_cnv_Pos_Chromosomes, PANNTHR_Pt6_OnTreat_panCK_cnv_Neg_Chromosomes))
colnames(PANNTHR_Pt6_OnTreat_panCK_cnv_Chromosome_Quant) <- c("Chromosomes CNV Pos", "Chromosomes CNV Neg")
write.csv(PANNTHR_Pt6_OnTreat_panCK_cnv_Chromosome_Quant, file = "/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Output/PANNTHR_Pt6_OnTreat_panCK/PANNTHR_Pt6_OnTreat_panCK_cnv_Chromosome_Quant.csv")

# Save CNV Metadata
PANNTHR_Pt6_OnTreat_panCK_cnv.meta.data <- as.data.frame(as.matrix(PANNTHR_Pt6_OnTreat_panCK_cnv@meta.data))
write.csv(PANNTHR_Pt6_OnTreat_panCK_cnv.meta.data, file = "/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Output/PANNTHR_Pt6_OnTreat_panCK/PANNTHR_Pt6_OnTreat_panCK_cnv.meta.data.csv")

# Subset Cells with inferred CNVs -------------------------------------------------------------------------------------------------------------
PANNTHR_Pt6_OnTreat_panCK_cnv <- subset(PANNTHR_Pt6_OnTreat_panCK_cnv, cells=rownames(PANNTHR_Pt6_OnTreat_panCK_cnv@meta.data)[rowSums(PANNTHR_Pt6_OnTreat_panCK_cnv@meta.data[,startsWith(names(PANNTHR_Pt6_OnTreat_panCK_cnv@meta.data),"has_cnv_")]) > 0 ])

# Save RDS with inferredCNVs
saveRDS(PANNTHR_Pt6_OnTreat_panCK_cnv, file = "/R/R_PANNTHR/PANNTHR_RDS/RDS_inferCNV/PANNTHR_Pt6_OnTreat_panCK_cnv.rds")


# Free memory
rm(PANNTHR_Pt6_OnTreat_inferCNV)
rm(PANNTHR_Pt6_OnTreat_Visvader_merged_NORM_sampled_panCK)
rm(PANNTHR_Pt6_OnTreat_panCK_cnv)
rm(PANNTHR_Pt6_OnTreat_panCK_cnv_Pos_Chromosomes)
rm(PANNTHR_Pt6_OnTreat_panCK_cnv_Neg_Chromosomes)
rm(PANNTHR_Pt6_OnTreat_panCK_cnv_Chromosome_Quant)
rm(PANNTHR_Pt6_OnTreat_panCK_cnv.meta.data)
gc()



#############################################
#### InferCNV: PANNTHR_Pt7_OnTreat_panCK ####
#############################################

########## Merge with normal mammary glad, extract epithelium, and process via Seurat
# Load RDS
PANNTHR_Pt7_OnTreat_panCK <- readRDS(file = "/R/R_PANNTHR/PANNTHR_RDS/RDS_panCK/PANNTHR_Pt7_OnTreat_panCK.rds")
Visvader_merged_NORM_sampled_panCK <- readRDS(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Visvader_merged_NORM_sampled_panCK.rds")

# Set identities 
colnames(x = PANNTHR_Pt7_OnTreat_panCK[[]])
Idents(object = PANNTHR_Pt7_OnTreat_panCK) <- 'Tumor'
PANNTHR_Pt7_OnTreat_panCK[["cnv.ident"]] <- Idents(object = PANNTHR_Pt7_OnTreat_panCK)
levels(x = PANNTHR_Pt7_OnTreat_panCK)

colnames(x = Visvader_merged_NORM_sampled_panCK[[]])
Idents(object = Visvader_merged_NORM_sampled_panCK) <- 'Normal'
Visvader_merged_NORM_sampled_panCK[["cnv.ident"]] <- Idents(object = Visvader_merged_NORM_sampled_panCK)
levels(x = Visvader_merged_NORM_sampled_panCK)

# Merge subsets into recombined RDS
PANNTHR_Pt7_OnTreat_Visvader_merged_NORM_sampled_panCK <- merge(x = PANNTHR_Pt7_OnTreat_panCK, y = Visvader_merged_NORM_sampled_panCK)
# Join Layers
PANNTHR_Pt7_OnTreat_Visvader_merged_NORM_sampled_panCK <- JoinLayers(PANNTHR_Pt7_OnTreat_Visvader_merged_NORM_sampled_panCK)

# Remove subsetted objects
rm(PANNTHR_Pt7_OnTreat_panCK)
rm(Visvader_merged_NORM_sampled_panCK)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
PANNTHR_Pt7_OnTreat_Visvader_merged_NORM_sampled_panCK <- FindVariableFeatures(PANNTHR_Pt7_OnTreat_Visvader_merged_NORM_sampled_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(PANNTHR_Pt7_OnTreat_Visvader_merged_NORM_sampled_panCK)
PANNTHR_Pt7_OnTreat_Visvader_merged_NORM_sampled_panCK <- ScaleData(PANNTHR_Pt7_OnTreat_Visvader_merged_NORM_sampled_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
PANNTHR_Pt7_OnTreat_Visvader_merged_NORM_sampled_panCK <- RunPCA(PANNTHR_Pt7_OnTreat_Visvader_merged_NORM_sampled_panCK, features = VariableFeatures(object = PANNTHR_Pt7_OnTreat_Visvader_merged_NORM_sampled_panCK))

# Examine and visualize PCA results a few different ways
print(PANNTHR_Pt7_OnTreat_Visvader_merged_NORM_sampled_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(PANNTHR_Pt7_OnTreat_Visvader_merged_NORM_sampled_panCK, dims = 1:2, reduction = "pca")
DimPlot(PANNTHR_Pt7_OnTreat_Visvader_merged_NORM_sampled_panCK, reduction = "pca")

PANNTHR_Pt7_OnTreat_Visvader_merged_NORM_sampled_panCK <- FindNeighbors(PANNTHR_Pt7_OnTreat_Visvader_merged_NORM_sampled_panCK, dims = 1:20)
PANNTHR_Pt7_OnTreat_Visvader_merged_NORM_sampled_panCK <- FindClusters(PANNTHR_Pt7_OnTreat_Visvader_merged_NORM_sampled_panCK, resolution = 0.1)

# Umap clustering
PANNTHR_Pt7_OnTreat_Visvader_merged_NORM_sampled_panCK <- RunUMAP(PANNTHR_Pt7_OnTreat_Visvader_merged_NORM_sampled_panCK, dims = 1:20)
DimPlot(PANNTHR_Pt7_OnTreat_Visvader_merged_NORM_sampled_panCK, reduction = "umap")
DimPlot(PANNTHR_Pt7_OnTreat_Visvader_merged_NORM_sampled_panCK, reduction = "umap", group.by = "orig.ident")

# Save RDS
saveRDS(PANNTHR_Pt7_OnTreat_Visvader_merged_NORM_sampled_panCK, file = "/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Input/PANNTHR_inferCNV_Input_RDS/PANNTHR_Pt7_OnTreat_Visvader_merged_NORM_sampled_panCK.rds")

########## Prepare Cell Annotations
# Extract & Save Meta Data
PANNTHR_Pt7_OnTreat_Visvader_merged_NORM_sampled_panCK.meta.data <- as.data.frame(as.matrix(PANNTHR_Pt7_OnTreat_Visvader_merged_NORM_sampled_panCK@meta.data))
# Extract pertinent columns
# Make row names first column
PANNTHR_Pt7_OnTreat_Visvader_merged_NORM_sampled_panCK.meta.data <- tibble::rownames_to_column(PANNTHR_Pt7_OnTreat_Visvader_merged_NORM_sampled_panCK.meta.data, "Cells")
# In Seurat, "X" becomes the cell names while "cnv.ident" is specified above
PANNTHR_Pt7_OnTreat_Visvader_merged_NORM_sampled_panCK.annotations <- PANNTHR_Pt7_OnTreat_Visvader_merged_NORM_sampled_panCK.meta.data[, c("Cells", "cnv.ident")]
# Delete header
names(PANNTHR_Pt7_OnTreat_Visvader_merged_NORM_sampled_panCK.annotations) <- NULL
# Save as table
write.table(PANNTHR_Pt7_OnTreat_Visvader_merged_NORM_sampled_panCK.annotations,"/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Input/PANNTHR_Pt7_OnTreat_Visvader_merged_NORM_sampled_panCK.annotations.txt",sep="\t",row.names=FALSE)
# Remove metadata and annotations
rm(PANNTHR_Pt7_OnTreat_Visvader_merged_NORM_sampled_panCK.meta.data)
rm(PANNTHR_Pt7_OnTreat_Visvader_merged_NORM_sampled_panCK.annotations)

########## Generate Count Matrix
# Extract counts from the Seurat object
PANNTHR_Pt7_OnTreat_Visvader_merged_NORM_sampled_panCK.count <- PANNTHR_Pt7_OnTreat_Visvader_merged_NORM_sampled_panCK[["RNA"]]$counts
# Convert counts to a sparse matrix
PANNTHR_Pt7_OnTreat_Visvader_merged_NORM_sampled_panCK.count.matrix <- as(PANNTHR_Pt7_OnTreat_Visvader_merged_NORM_sampled_panCK.count, 'sparseMatrix')
rm(PANNTHR_Pt7_OnTreat_Visvader_merged_NORM_sampled_panCK.count)

# Convert to requisite data format
PANNTHR_Pt7_OnTreat_Visvader_merged_NORM_sampled_panCK.count.matrix <- as.data.frame(PANNTHR_Pt7_OnTreat_Visvader_merged_NORM_sampled_panCK.count.matrix)
# Save table
write.csv(PANNTHR_Pt7_OnTreat_Visvader_merged_NORM_sampled_panCK.count.matrix, file = "/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Input/PANNTHR_inferCNV_Input_Count_Matrix/PANNTHR_Pt7_OnTreat_Visvader_merged_NORM_sampled_panCK.count.matrix.csv")
# Remove Objects
rm(PANNTHR_Pt7_OnTreat_Visvader_merged_NORM_sampled_panCK)
rm(PANNTHR_Pt7_OnTreat_Visvader_merged_NORM_sampled_panCK.count.matrix)

########### Run inferCNV
# Note: check.names = FALSE prevents cell names from having dashes replaced with dots; essential to match annotations file
PANNTHR_Pt7_OnTreat_inferCNV <- CreateInfercnvObject(raw_counts_matrix=read.csv(file = "/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Input/PANNTHR_inferCNV_Input_Count_Matrix/PANNTHR_Pt7_OnTreat_Visvader_merged_NORM_sampled_panCK.count.matrix.csv", header = TRUE, row.names = 1,check.names = FALSE),
                                               annotations_file=read.table(file = "/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Input/PANNTHR_Pt7_OnTreat_Visvader_merged_NORM_sampled_panCK.annotations.txt", header = FALSE, row.names = 1),
                                               delim="\t",
                                               gene_order_file="/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Input/hg38_gencode_v27.txt",
                                               ref_group_names=c("Normal"))


# perform infercnv operations to reveal cnv signal
PANNTHR_Pt7_OnTreat_inferCNV = infercnv::run(PANNTHR_Pt7_OnTreat_inferCNV,
                                       cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                       out_dir="/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Output/PANNTHR_Pt7_OnTreat_panCK/",  # dir is auto-created for storing outputs
                                       cluster_by_groups=T,   # cluster
                                       denoise=T,
                                       HMM=T,
                                       num_ref_groups = 8, analysis_mode = "subclusters")

###########  Integrate data with Seurat object
# Read RDS and add inferCNV results to metadata
PANNTHR_Pt7_OnTreat_Visvader_merged_NORM_sampled_panCK <- readRDS(file = "/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Input/PANNTHR_inferCNV_Input_RDS/PANNTHR_Pt7_OnTreat_Visvader_merged_NORM_sampled_panCK.rds") 
PANNTHR_Pt7_OnTreat_panCK_cnv <- infercnv::add_to_seurat(infercnv_output_path="/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Output/PANNTHR_Pt7_OnTreat_panCK/",
                                                   seurat_obj=PANNTHR_Pt7_OnTreat_Visvader_merged_NORM_sampled_panCK, # optional
                                                   top_n=10)

# Collect Tumor cells from Seurat Object
PANNTHR_Pt7_OnTreat_panCK_cnv <- subset(x = PANNTHR_Pt7_OnTreat_panCK_cnv, subset = cnv.ident == "Tumor")
# Quantify chromosomes with alterations per cell
PANNTHR_Pt7_OnTreat_panCK_cnv_Pos_Chromosomes <- as.data.frame(rowSums(PANNTHR_Pt7_OnTreat_panCK_cnv@meta.data[,startsWith(names(PANNTHR_Pt7_OnTreat_panCK_cnv@meta.data),"has_cnv_")] > 0 ))
PANNTHR_Pt7_OnTreat_panCK_cnv_Neg_Chromosomes <- as.data.frame(rowSums(PANNTHR_Pt7_OnTreat_panCK_cnv@meta.data[,startsWith(names(PANNTHR_Pt7_OnTreat_panCK_cnv@meta.data),"has_cnv_")] == 0))
PANNTHR_Pt7_OnTreat_panCK_cnv_Chromosome_Quant <- as.data.frame(c(PANNTHR_Pt7_OnTreat_panCK_cnv_Pos_Chromosomes, PANNTHR_Pt7_OnTreat_panCK_cnv_Neg_Chromosomes))
colnames(PANNTHR_Pt7_OnTreat_panCK_cnv_Chromosome_Quant) <- c("Chromosomes CNV Pos", "Chromosomes CNV Neg")
write.csv(PANNTHR_Pt7_OnTreat_panCK_cnv_Chromosome_Quant, file = "/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Output/PANNTHR_Pt7_OnTreat_panCK/PANNTHR_Pt7_OnTreat_panCK_cnv_Chromosome_Quant.csv")

# Save CNV Metadata
PANNTHR_Pt7_OnTreat_panCK_cnv.meta.data <- as.data.frame(as.matrix(PANNTHR_Pt7_OnTreat_panCK_cnv@meta.data))
write.csv(PANNTHR_Pt7_OnTreat_panCK_cnv.meta.data, file = "/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Output/PANNTHR_Pt7_OnTreat_panCK/PANNTHR_Pt7_OnTreat_panCK_cnv.meta.data.csv")

# Subset Cells with inferred CNVs -------------------------------------------------------------------------------------------------------------
PANNTHR_Pt7_OnTreat_panCK_cnv <- subset(PANNTHR_Pt7_OnTreat_panCK_cnv, cells=rownames(PANNTHR_Pt7_OnTreat_panCK_cnv@meta.data)[rowSums(PANNTHR_Pt7_OnTreat_panCK_cnv@meta.data[,startsWith(names(PANNTHR_Pt7_OnTreat_panCK_cnv@meta.data),"has_cnv_")]) > 0 ])

# Save RDS with inferredCNVs
saveRDS(PANNTHR_Pt7_OnTreat_panCK_cnv, file = "/R/R_PANNTHR/PANNTHR_RDS/RDS_inferCNV/PANNTHR_Pt7_OnTreat_panCK_cnv.rds")


# Free memory
rm(PANNTHR_Pt7_OnTreat_inferCNV)
rm(PANNTHR_Pt7_OnTreat_Visvader_merged_NORM_sampled_panCK)
rm(PANNTHR_Pt7_OnTreat_panCK_cnv)
rm(PANNTHR_Pt7_OnTreat_panCK_cnv_Pos_Chromosomes)
rm(PANNTHR_Pt7_OnTreat_panCK_cnv_Neg_Chromosomes)
rm(PANNTHR_Pt7_OnTreat_panCK_cnv_Chromosome_Quant)
rm(PANNTHR_Pt7_OnTreat_panCK_cnv.meta.data)
gc()



###############################################
#### InferCNV: PANNTHR_Pt5_PostTreat_panCK ####
###############################################

########## Merge with normal mammary glad, extract epithelium, and process via Seurat
# Load RDS
PANNTHR_Pt5_PostTreat_panCK <- readRDS(file = "/R/R_PANNTHR/PANNTHR_RDS/RDS_panCK/PANNTHR_Pt5_PostTreat_panCK.rds")
Visvader_merged_NORM_sampled_panCK <- readRDS(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Visvader_merged_NORM_sampled_panCK.rds")

# Set identities 
colnames(x = PANNTHR_Pt5_PostTreat_panCK[[]])
Idents(object = PANNTHR_Pt5_PostTreat_panCK) <- 'Tumor'
PANNTHR_Pt5_PostTreat_panCK[["cnv.ident"]] <- Idents(object = PANNTHR_Pt5_PostTreat_panCK)
levels(x = PANNTHR_Pt5_PostTreat_panCK)

colnames(x = Visvader_merged_NORM_sampled_panCK[[]])
Idents(object = Visvader_merged_NORM_sampled_panCK) <- 'Normal'
Visvader_merged_NORM_sampled_panCK[["cnv.ident"]] <- Idents(object = Visvader_merged_NORM_sampled_panCK)
levels(x = Visvader_merged_NORM_sampled_panCK)

# Merge subsets into recombined RDS
PANNTHR_Pt5_PostTreat_Visvader_merged_NORM_sampled_panCK <- merge(x = PANNTHR_Pt5_PostTreat_panCK, y = Visvader_merged_NORM_sampled_panCK)
# Join Layers
PANNTHR_Pt5_PostTreat_Visvader_merged_NORM_sampled_panCK <- JoinLayers(PANNTHR_Pt5_PostTreat_Visvader_merged_NORM_sampled_panCK)

# Remove subsetted objects
rm(PANNTHR_Pt5_PostTreat_panCK)
rm(Visvader_merged_NORM_sampled_panCK)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
PANNTHR_Pt5_PostTreat_Visvader_merged_NORM_sampled_panCK <- FindVariableFeatures(PANNTHR_Pt5_PostTreat_Visvader_merged_NORM_sampled_panCK, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(PANNTHR_Pt5_PostTreat_Visvader_merged_NORM_sampled_panCK)
PANNTHR_Pt5_PostTreat_Visvader_merged_NORM_sampled_panCK <- ScaleData(PANNTHR_Pt5_PostTreat_Visvader_merged_NORM_sampled_panCK, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
PANNTHR_Pt5_PostTreat_Visvader_merged_NORM_sampled_panCK <- RunPCA(PANNTHR_Pt5_PostTreat_Visvader_merged_NORM_sampled_panCK, features = VariableFeatures(object = PANNTHR_Pt5_PostTreat_Visvader_merged_NORM_sampled_panCK))

# Examine and visualize PCA results a few different ways
print(PANNTHR_Pt5_PostTreat_Visvader_merged_NORM_sampled_panCK[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(PANNTHR_Pt5_PostTreat_Visvader_merged_NORM_sampled_panCK, dims = 1:2, reduction = "pca")
DimPlot(PANNTHR_Pt5_PostTreat_Visvader_merged_NORM_sampled_panCK, reduction = "pca")

PANNTHR_Pt5_PostTreat_Visvader_merged_NORM_sampled_panCK <- FindNeighbors(PANNTHR_Pt5_PostTreat_Visvader_merged_NORM_sampled_panCK, dims = 1:20)
PANNTHR_Pt5_PostTreat_Visvader_merged_NORM_sampled_panCK <- FindClusters(PANNTHR_Pt5_PostTreat_Visvader_merged_NORM_sampled_panCK, resolution = 0.1)

# Umap clustering
PANNTHR_Pt5_PostTreat_Visvader_merged_NORM_sampled_panCK <- RunUMAP(PANNTHR_Pt5_PostTreat_Visvader_merged_NORM_sampled_panCK, dims = 1:20)
DimPlot(PANNTHR_Pt5_PostTreat_Visvader_merged_NORM_sampled_panCK, reduction = "umap")
DimPlot(PANNTHR_Pt5_PostTreat_Visvader_merged_NORM_sampled_panCK, reduction = "umap", group.by = "orig.ident")

# Save RDS
saveRDS(PANNTHR_Pt5_PostTreat_Visvader_merged_NORM_sampled_panCK, file = "/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Input/PANNTHR_inferCNV_Input_RDS/PANNTHR_Pt5_PostTreat_Visvader_merged_NORM_sampled_panCK.rds")

########## Prepare Cell Annotations
# Extract & Save Meta Data
PANNTHR_Pt5_PostTreat_Visvader_merged_NORM_sampled_panCK.meta.data <- as.data.frame(as.matrix(PANNTHR_Pt5_PostTreat_Visvader_merged_NORM_sampled_panCK@meta.data))
# Extract pertinent columns
# Make row names first column
PANNTHR_Pt5_PostTreat_Visvader_merged_NORM_sampled_panCK.meta.data <- tibble::rownames_to_column(PANNTHR_Pt5_PostTreat_Visvader_merged_NORM_sampled_panCK.meta.data, "Cells")
# In Seurat, "X" becomes the cell names while "cnv.ident" is specified above
PANNTHR_Pt5_PostTreat_Visvader_merged_NORM_sampled_panCK.annotations <- PANNTHR_Pt5_PostTreat_Visvader_merged_NORM_sampled_panCK.meta.data[, c("Cells", "cnv.ident")]
# Delete header
names(PANNTHR_Pt5_PostTreat_Visvader_merged_NORM_sampled_panCK.annotations) <- NULL
# Save as table
write.table(PANNTHR_Pt5_PostTreat_Visvader_merged_NORM_sampled_panCK.annotations,"/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Input/PANNTHR_Pt5_PostTreat_Visvader_merged_NORM_sampled_panCK.annotations.txt",sep="\t",row.names=FALSE)
# Remove metadata and annotations
rm(PANNTHR_Pt5_PostTreat_Visvader_merged_NORM_sampled_panCK.meta.data)
rm(PANNTHR_Pt5_PostTreat_Visvader_merged_NORM_sampled_panCK.annotations)

########## Generate Count Matrix
# Extract counts from the Seurat object
PANNTHR_Pt5_PostTreat_Visvader_merged_NORM_sampled_panCK.count <- PANNTHR_Pt5_PostTreat_Visvader_merged_NORM_sampled_panCK[["RNA"]]$counts
# Convert counts to a sparse matrix
PANNTHR_Pt5_PostTreat_Visvader_merged_NORM_sampled_panCK.count.matrix <- as(PANNTHR_Pt5_PostTreat_Visvader_merged_NORM_sampled_panCK.count, 'sparseMatrix')
rm(PANNTHR_Pt5_PostTreat_Visvader_merged_NORM_sampled_panCK.count)

# Convert to requisite data format
PANNTHR_Pt5_PostTreat_Visvader_merged_NORM_sampled_panCK.count.matrix <- as.data.frame(PANNTHR_Pt5_PostTreat_Visvader_merged_NORM_sampled_panCK.count.matrix)
# Save table
write.csv(PANNTHR_Pt5_PostTreat_Visvader_merged_NORM_sampled_panCK.count.matrix, file = "/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Input/PANNTHR_inferCNV_Input_Count_Matrix/PANNTHR_Pt5_PostTreat_Visvader_merged_NORM_sampled_panCK.count.matrix.csv")
# Remove Objects
rm(PANNTHR_Pt5_PostTreat_Visvader_merged_NORM_sampled_panCK)
rm(PANNTHR_Pt5_PostTreat_Visvader_merged_NORM_sampled_panCK.count.matrix)

########### Run inferCNV
# Note: check.names = FALSE prevents cell names from having dashes replaced with dots; essential to match annotations file
PANNTHR_Pt5_PostTreat_inferCNV <- CreateInfercnvObject(raw_counts_matrix=read.csv(file = "/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Input/PANNTHR_inferCNV_Input_Count_Matrix/PANNTHR_Pt5_PostTreat_Visvader_merged_NORM_sampled_panCK.count.matrix.csv", header = TRUE, row.names = 1,check.names = FALSE),
                                               annotations_file=read.table(file = "/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Input/PANNTHR_Pt5_PostTreat_Visvader_merged_NORM_sampled_panCK.annotations.txt", header = FALSE, row.names = 1),
                                               delim="\t",
                                               gene_order_file="/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Input/hg38_gencode_v27.txt",
                                               ref_group_names=c("Normal"))


# perform infercnv operations to reveal cnv signal
PANNTHR_Pt5_PostTreat_inferCNV = infercnv::run(PANNTHR_Pt5_PostTreat_inferCNV,
                                       cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                       out_dir="/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Output/PANNTHR_Pt5_PostTreat_panCK/",  # dir is auto-created for storing outputs
                                       cluster_by_groups=T,   # cluster
                                       denoise=T,
                                       HMM=T,
                                       num_ref_groups = 8, analysis_mode = "subclusters")

###########  Integrate data with Seurat object
# Read RDS and add inferCNV results to metadata
PANNTHR_Pt5_PostTreat_Visvader_merged_NORM_sampled_panCK <- readRDS(file = "/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Input/PANNTHR_inferCNV_Input_RDS/PANNTHR_Pt5_PostTreat_Visvader_merged_NORM_sampled_panCK.rds") 
PANNTHR_Pt5_PostTreat_panCK_cnv <- infercnv::add_to_seurat(infercnv_output_path="/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Output/PANNTHR_Pt5_PostTreat_panCK/",
                                                   seurat_obj=PANNTHR_Pt5_PostTreat_Visvader_merged_NORM_sampled_panCK, # optional
                                                   top_n=10)

# Collect Tumor cells from Seurat Object
PANNTHR_Pt5_PostTreat_panCK_cnv <- subset(x = PANNTHR_Pt5_PostTreat_panCK_cnv, subset = cnv.ident == "Tumor")
# Quantify chromosomes with alterations per cell
PANNTHR_Pt5_PostTreat_panCK_cnv_Pos_Chromosomes <- as.data.frame(rowSums(PANNTHR_Pt5_PostTreat_panCK_cnv@meta.data[,startsWith(names(PANNTHR_Pt5_PostTreat_panCK_cnv@meta.data),"has_cnv_")] > 0 ))
PANNTHR_Pt5_PostTreat_panCK_cnv_Neg_Chromosomes <- as.data.frame(rowSums(PANNTHR_Pt5_PostTreat_panCK_cnv@meta.data[,startsWith(names(PANNTHR_Pt5_PostTreat_panCK_cnv@meta.data),"has_cnv_")] == 0))
PANNTHR_Pt5_PostTreat_panCK_cnv_Chromosome_Quant <- as.data.frame(c(PANNTHR_Pt5_PostTreat_panCK_cnv_Pos_Chromosomes, PANNTHR_Pt5_PostTreat_panCK_cnv_Neg_Chromosomes))
colnames(PANNTHR_Pt5_PostTreat_panCK_cnv_Chromosome_Quant) <- c("Chromosomes CNV Pos", "Chromosomes CNV Neg")
write.csv(PANNTHR_Pt5_PostTreat_panCK_cnv_Chromosome_Quant, file = "/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Output/PANNTHR_Pt5_PostTreat_panCK/PANNTHR_Pt5_PostTreat_panCK_cnv_Chromosome_Quant.csv")

# Save CNV Metadata
PANNTHR_Pt5_PostTreat_panCK_cnv.meta.data <- as.data.frame(as.matrix(PANNTHR_Pt5_PostTreat_panCK_cnv@meta.data))
write.csv(PANNTHR_Pt5_PostTreat_panCK_cnv.meta.data, file = "/R/R_PANNTHR/PANNTHR_inferCNV/PANNTHR_inferCNV_Output/PANNTHR_Pt5_PostTreat_panCK/PANNTHR_Pt5_PostTreat_panCK_cnv.meta.data.csv")

# Subset Cells with inferred CNVs -------------------------------------------------------------------------------------------------------------
PANNTHR_Pt5_PostTreat_panCK_cnv <- subset(PANNTHR_Pt5_PostTreat_panCK_cnv, cells=rownames(PANNTHR_Pt5_PostTreat_panCK_cnv@meta.data)[rowSums(PANNTHR_Pt5_PostTreat_panCK_cnv@meta.data[,startsWith(names(PANNTHR_Pt5_PostTreat_panCK_cnv@meta.data),"has_cnv_")]) > 0 ])

# Save RDS with inferredCNVs
saveRDS(PANNTHR_Pt5_PostTreat_panCK_cnv, file = "/R/R_PANNTHR/PANNTHR_RDS/RDS_inferCNV/PANNTHR_Pt5_PostTreat_panCK_cnv.rds")


# Free memory
rm(PANNTHR_Pt5_PostTreat_inferCNV)
rm(PANNTHR_Pt5_PostTreat_Visvader_merged_NORM_sampled_panCK)
rm(PANNTHR_Pt5_PostTreat_panCK_cnv)
rm(PANNTHR_Pt5_PostTreat_panCK_cnv_Pos_Chromosomes)
rm(PANNTHR_Pt5_PostTreat_panCK_cnv_Neg_Chromosomes)
rm(PANNTHR_Pt5_PostTreat_panCK_cnv_Chromosome_Quant)
rm(PANNTHR_Pt5_PostTreat_panCK_cnv.meta.data)
gc()


#############################################################################
# Step 6: Merge Neoplastic Cells from Different Tumors into Combined Object #  
#############################################################################
# Load Libraries
library(Seurat)
library(dplyr)
library(patchwork)

#########################################################
############# This is PANNTHR Tumor Batch # #############
#########################################################


#########################
# Load Individual Files #
#########################
PANNTHR_Pt1_PreTreat_panCK_cnv <- readRDS(file = "/R/R_PANNTHR/PANNTHR_RDS/RDS_inferCNV/PANNTHR_Pt1_PreTreat_panCK_cnv.rds")
PANNTHR_Pt2_PreTreat_panCK_cnv <- readRDS(file = "/R/R_PANNTHR/PANNTHR_RDS/RDS_inferCNV/PANNTHR_Pt2_PreTreat_panCK_cnv.rds")
PANNTHR_Pt4_PreTreat_panCK_cnv <- readRDS(file = "/R/R_PANNTHR/PANNTHR_RDS/RDS_inferCNV/PANNTHR_Pt4_PreTreat_panCK_cnv.rds")
PANNTHR_Pt5_PreTreat_panCK_cnv <- readRDS(file = "/R/R_PANNTHR/PANNTHR_RDS/RDS_inferCNV/PANNTHR_Pt5_PreTreat_panCK_cnv.rds")
PANNTHR_Pt2_OnTreat_panCK_cnv <- readRDS(file = "/R/R_PANNTHR/PANNTHR_RDS/RDS_inferCNV/PANNTHR_Pt2_OnTreat_panCK_cnv.rds")
PANNTHR_Pt3_OnTreat_panCK_cnv <- readRDS(file = "/R/R_PANNTHR/PANNTHR_RDS/RDS_inferCNV/PANNTHR_Pt3_OnTreat_panCK_cnv.rds")
PANNTHR_Pt4_OnTreat_panCK_cnv <- readRDS(file = "/R/R_PANNTHR/PANNTHR_RDS/RDS_inferCNV/PANNTHR_Pt4_OnTreat_panCK_cnv.rds")
PANNTHR_Pt6_OnTreat_panCK_cnv <- readRDS(file = "/R/R_PANNTHR/PANNTHR_RDS/RDS_inferCNV/PANNTHR_Pt6_OnTreat_panCK_cnv.rds")
PANNTHR_Pt7_OnTreat_panCK_cnv <- readRDS(file = "/R/R_PANNTHR/PANNTHR_RDS/RDS_inferCNV/PANNTHR_Pt7_OnTreat_panCK_cnv.rds")
PANNTHR_Pt5_PostTreat_panCK_cnv <- readRDS(file = "/R/R_PANNTHR/PANNTHR_RDS/RDS_inferCNV/PANNTHR_Pt5_PostTreat_panCK_cnv.rds")

##############################################
# Merge Individual Files into Single Dataset #
##############################################
# Merge Datasets
all_PANNTHR_breast_tumors_panCK_cnv <- merge(x = PANNTHR_Pt1_PreTreat_panCK_cnv, y = c(PANNTHR_Pt2_PreTreat_panCK_cnv, PANNTHR_Pt4_PreTreat_panCK_cnv, PANNTHR_Pt5_PreTreat_panCK_cnv, 
                                                                                       PANNTHR_Pt2_OnTreat_panCK_cnv, PANNTHR_Pt3_OnTreat_panCK_cnv, PANNTHR_Pt4_OnTreat_panCK_cnv,
                                                                                       PANNTHR_Pt6_OnTreat_panCK_cnv, PANNTHR_Pt7_OnTreat_panCK_cnv, PANNTHR_Pt5_PostTreat_panCK_cnv))

# Join Layers
all_PANNTHR_breast_tumors_panCK_cnv <- JoinLayers(all_PANNTHR_breast_tumors_panCK_cnv)

# Remove Objects
rm(PANNTHR_Pt1_PreTreat_panCK_cnv)
rm(PANNTHR_Pt2_PreTreat_panCK_cnv)
rm(PANNTHR_Pt4_PreTreat_panCK_cnv)
rm(PANNTHR_Pt5_PreTreat_panCK_cnv)
rm(PANNTHR_Pt2_OnTreat_panCK_cnv)
rm(PANNTHR_Pt3_OnTreat_panCK_cnv)
rm(PANNTHR_Pt4_OnTreat_panCK_cnv)
rm(PANNTHR_Pt6_OnTreat_panCK_cnv)
rm(PANNTHR_Pt7_OnTreat_panCK_cnv)
rm(PANNTHR_Pt5_PostTreat_panCK_cnv)
gc()



####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
all_PANNTHR_breast_tumors_panCK_cnv <- FindVariableFeatures(all_PANNTHR_breast_tumors_panCK_cnv, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(all_PANNTHR_breast_tumors_panCK_cnv)
all_PANNTHR_breast_tumors_panCK_cnv <- ScaleData(all_PANNTHR_breast_tumors_panCK_cnv, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
all_PANNTHR_breast_tumors_panCK_cnv <- RunPCA(all_PANNTHR_breast_tumors_panCK_cnv, features = VariableFeatures(object = all_PANNTHR_breast_tumors_panCK_cnv))

# Examine and visualize PCA results a few different ways
print(all_PANNTHR_breast_tumors_panCK_cnv[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(all_PANNTHR_breast_tumors_panCK_cnv, dims = 1:2, reduction = "pca")
DimPlot(all_PANNTHR_breast_tumors_panCK_cnv, reduction = "pca")

# Visualize PCs
ElbowPlot(all_PANNTHR_breast_tumors_panCK_cnv)

# Choose dimensions
all_PANNTHR_breast_tumors_panCK_cnv <- FindNeighbors(all_PANNTHR_breast_tumors_panCK_cnv, dims = 1:12)
all_PANNTHR_breast_tumors_panCK_cnv <- FindClusters(all_PANNTHR_breast_tumors_panCK_cnv, resolution = 0.1)

# Umap clustering
all_PANNTHR_breast_tumors_panCK_cnv <- RunUMAP(all_PANNTHR_breast_tumors_panCK_cnv, dims = 1:12)
DimPlot(all_PANNTHR_breast_tumors_panCK_cnv, reduction = "umap", label = TRUE)

# Identify Samples
DimPlot(all_PANNTHR_breast_tumors_panCK_cnv, reduction = "umap", group.by = "orig.ident")

# Save
saveRDS(all_PANNTHR_breast_tumors_panCK_cnv, file = "/R/R_PANNTHR/PANNTHR_RDS/RDS_Merged/all_PANNTHR_breast_tumors_panCK_cnv.rds")
# Clear global environment
#rm(all_PANNTHR_breast_tumors_panCK_cnv)
#gc()




