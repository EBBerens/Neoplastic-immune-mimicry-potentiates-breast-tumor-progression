#####################################################################################################################
#                       Visvader (Pal et al) Dataset Breast Tumor Analysis Steps                                    #
#####################################################################################################################
# Step 5: Run InferCNV on Breast Tumors using Reduction Mammoplasty Reference and Subset Neoplastic Cells           #
# Step 6: Merge Neoplastic Cells from Different Tumors into Combined Object for Analysis                            #
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
#############################     InferCNV: HR-pos Tumor Specimens      ###########################    
###################################################################################################

##########################################
#### InferCNV: Visvader_0001_ER_panCK ####
##########################################

########## Merge with normal mammary glad, extract epithelium, and process via Seurat
# Load RDS
Visvader_0001_ER_panCK <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0001_ER_panCK.rds")
Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Navin_Visvader_NORM_panCK_Reference.rds")

# Set identities 
colnames(x = Visvader_0001_ER_panCK[[]])
Idents(object = Visvader_0001_ER_panCK) <- 'Tumor'
Visvader_0001_ER_panCK[["cnv.ident"]] <- Idents(object = Visvader_0001_ER_panCK)
levels(x = Visvader_0001_ER_panCK)

colnames(x = Navin_Visvader_NORM_panCK_Reference[[]])
Idents(object = Navin_Visvader_NORM_panCK_Reference) <- 'Normal'
Navin_Visvader_NORM_panCK_Reference[["cnv.ident"]] <- Idents(object = Navin_Visvader_NORM_panCK_Reference)
levels(x = Navin_Visvader_NORM_panCK_Reference)

# Merge subsets into recombined RDS
Visvader_0001_ER_Navin_Visvader_NORM_panCK_Reference <- merge(x = Visvader_0001_ER_panCK, y = Navin_Visvader_NORM_panCK_Reference)
# Join Layers
Visvader_0001_ER_Navin_Visvader_NORM_panCK_Reference <- JoinLayers(Visvader_0001_ER_Navin_Visvader_NORM_panCK_Reference)

# Remove subsetted objects
rm(Visvader_0001_ER_panCK)
rm(Navin_Visvader_NORM_panCK_Reference)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Visvader_0001_ER_Navin_Visvader_NORM_panCK_Reference <- FindVariableFeatures(Visvader_0001_ER_Navin_Visvader_NORM_panCK_Reference, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Visvader_0001_ER_Navin_Visvader_NORM_panCK_Reference)
Visvader_0001_ER_Navin_Visvader_NORM_panCK_Reference <- ScaleData(Visvader_0001_ER_Navin_Visvader_NORM_panCK_Reference, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Visvader_0001_ER_Navin_Visvader_NORM_panCK_Reference <- RunPCA(Visvader_0001_ER_Navin_Visvader_NORM_panCK_Reference, features = VariableFeatures(object = Visvader_0001_ER_Navin_Visvader_NORM_panCK_Reference))

# Examine and visualize PCA results a few different ways
print(Visvader_0001_ER_Navin_Visvader_NORM_panCK_Reference[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Visvader_0001_ER_Navin_Visvader_NORM_panCK_Reference, dims = 1:2, reduction = "pca")
DimPlot(Visvader_0001_ER_Navin_Visvader_NORM_panCK_Reference, reduction = "pca")

Visvader_0001_ER_Navin_Visvader_NORM_panCK_Reference <- FindNeighbors(Visvader_0001_ER_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
Visvader_0001_ER_Navin_Visvader_NORM_panCK_Reference <- FindClusters(Visvader_0001_ER_Navin_Visvader_NORM_panCK_Reference, resolution = 0.1)

# Umap clustering
Visvader_0001_ER_Navin_Visvader_NORM_panCK_Reference <- RunUMAP(Visvader_0001_ER_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
DimPlot(Visvader_0001_ER_Navin_Visvader_NORM_panCK_Reference, reduction = "umap")
DimPlot(Visvader_0001_ER_Navin_Visvader_NORM_panCK_Reference, reduction = "umap", group.by = "orig.ident")

# Save RDS
saveRDS(Visvader_0001_ER_Navin_Visvader_NORM_panCK_Reference, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Visvader_0001_ER_Navin_Visvader_NORM_panCK_Reference.rds")

########## Prepare Cell Annotations
# Extract & Save Meta Data
Visvader_0001_ER_Navin_Visvader_NORM_panCK_Reference.meta.data <- as.data.frame(as.matrix(Visvader_0001_ER_Navin_Visvader_NORM_panCK_Reference@meta.data))
# Extract pertinent columns
# Make row names first column
Visvader_0001_ER_Navin_Visvader_NORM_panCK_Reference.meta.data <- tibble::rownames_to_column(Visvader_0001_ER_Navin_Visvader_NORM_panCK_Reference.meta.data, "Cells")
# In Seurat, "X" becomes the cell names while "cnv.ident" is specified above
Visvader_0001_ER_Navin_Visvader_NORM_panCK_Reference.annotations <- Visvader_0001_ER_Navin_Visvader_NORM_panCK_Reference.meta.data[, c("Cells", "cnv.ident")]
# Delete header
names(Visvader_0001_ER_Navin_Visvader_NORM_panCK_Reference.annotations) <- NULL
# Save as table
write.table(Visvader_0001_ER_Navin_Visvader_NORM_panCK_Reference.annotations,"/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_0001_ER_Navin_Visvader_NORM_panCK_Reference.annotations.txt",sep="\t",row.names=FALSE)
# Remove metadata and annotations
rm(Visvader_0001_ER_Navin_Visvader_NORM_panCK_Reference.meta.data)
rm(Visvader_0001_ER_Navin_Visvader_NORM_panCK_Reference.annotations)

########## Generate Count Matrix
# Extract counts from the Seurat object
Visvader_0001_ER_Navin_Visvader_NORM_panCK_Reference.count <- Visvader_0001_ER_Navin_Visvader_NORM_panCK_Reference[["RNA"]]$counts
# Convert counts to a sparse matrix
Visvader_0001_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as(Visvader_0001_ER_Navin_Visvader_NORM_panCK_Reference.count, 'sparseMatrix')
rm(Visvader_0001_ER_Navin_Visvader_NORM_panCK_Reference.count)

# Convert to requisite data format
Visvader_0001_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as.data.frame(Visvader_0001_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix)
# Save table
write.csv(Visvader_0001_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_Count_Matrix/Visvader_0001_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv")
# Remove Objects
rm(Visvader_0001_ER_Navin_Visvader_NORM_panCK_Reference)
rm(Visvader_0001_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix)

########### Run inferCNV
# Note: check.names = FALSE prevents cell names from having dashes replaced with dots; essential to match annotations file
Visvader_0001_ER_inferCNV <- CreateInfercnvObject(raw_counts_matrix=read.csv(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_Count_Matrix/Visvader_0001_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv", header = TRUE, row.names = 1,check.names = FALSE),
                                                  annotations_file=read.table(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_0001_ER_Navin_Visvader_NORM_panCK_Reference.annotations.txt", header = FALSE, row.names = 1),
                                                  delim="\t",
                                                  gene_order_file="/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/hg38_gencode_v27.txt",
                                                  ref_group_names=c("Normal"))


# perform infercnv operations to reveal cnv signal
Visvader_0001_ER_inferCNV = infercnv::run(Visvader_0001_ER_inferCNV,
                                          cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                          out_dir="/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0001_ER_panCK/",  # dir is auto-created for storing outputs
                                          cluster_by_groups=T,   # cluster
                                          denoise=T,
                                          HMM=T,
                                          num_ref_groups = 8, analysis_mode = "subclusters")

###########  Integrate data with Seurat object
# Read RDS and add inferCNV results to metadata
Visvader_0001_ER_Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Visvader_0001_ER_Navin_Visvader_NORM_panCK_Reference.rds") 
Visvader_0001_ER_panCK_cnv <- infercnv::add_to_seurat(infercnv_output_path="/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0001_ER_panCK/",
                                                      seurat_obj=Visvader_0001_ER_Navin_Visvader_NORM_panCK_Reference, # optional
                                                      top_n=10)

# Collect Tumor cells from Seurat Object
Visvader_0001_ER_panCK_cnv <- subset(x = Visvader_0001_ER_panCK_cnv, subset = cnv.ident == "Tumor")
# Quantify chromosomes with alterations per cell
Visvader_0001_ER_panCK_cnv_Pos_Chromosomes <- as.data.frame(rowSums(Visvader_0001_ER_panCK_cnv@meta.data[,startsWith(names(Visvader_0001_ER_panCK_cnv@meta.data),"has_cnv_")] > 0 ))
Visvader_0001_ER_panCK_cnv_Neg_Chromosomes <- as.data.frame(rowSums(Visvader_0001_ER_panCK_cnv@meta.data[,startsWith(names(Visvader_0001_ER_panCK_cnv@meta.data),"has_cnv_")] == 0))
Visvader_0001_ER_panCK_cnv_Chromosome_Quant <- as.data.frame(c(Visvader_0001_ER_panCK_cnv_Pos_Chromosomes, Visvader_0001_ER_panCK_cnv_Neg_Chromosomes))
colnames(Visvader_0001_ER_panCK_cnv_Chromosome_Quant) <- c("Chromosomes CNV Pos", "Chromosomes CNV Neg")
write.csv(Visvader_0001_ER_panCK_cnv_Chromosome_Quant, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0001_ER_panCK/Visvader_0001_ER_panCK_cnv_Chromosome_Quant.csv")

# Save CNV Metadata
Visvader_0001_ER_panCK_cnv.meta.data <- as.data.frame(as.matrix(Visvader_0001_ER_panCK_cnv@meta.data))
write.csv(Visvader_0001_ER_panCK_cnv.meta.data, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0001_ER_panCK/Visvader_0001_ER_panCK_cnv.meta.data.csv")

# Subset Cells with inferred CNVs -------------------------------------------------------------------------------------------------------------
Visvader_0001_ER_panCK_cnv <- subset(Visvader_0001_ER_panCK_cnv, cells=rownames(Visvader_0001_ER_panCK_cnv@meta.data)[rowSums(Visvader_0001_ER_panCK_cnv@meta.data[,startsWith(names(Visvader_0001_ER_panCK_cnv@meta.data),"has_cnv_")]) > 0 ])

# Save RDS with inferredCNVs
saveRDS(Visvader_0001_ER_panCK_cnv, file = "/R/R_Visvader/Visvader_RDS/RDS_inferCNV/Visvader_0001_ER_panCK_cnv.rds")


# Free memory
rm(Visvader_0001_ER_inferCNV)
rm(Visvader_0001_ER_Navin_Visvader_NORM_panCK_Reference)
rm(Visvader_0001_ER_panCK_cnv)
rm(Visvader_0001_ER_panCK_cnv_Pos_Chromosomes)
rm(Visvader_0001_ER_panCK_cnv_Neg_Chromosomes)
rm(Visvader_0001_ER_panCK_cnv_Chromosome_Quant)
rm(Visvader_0001_ER_panCK_cnv.meta.data)
gc()



##########################################
#### InferCNV: Visvader_0025_ER_panCK ####
##########################################

########## Merge with normal mammary glad, extract epithelium, and process via Seurat
# Load RDS
Visvader_0025_ER_panCK <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0025_ER_panCK.rds")
Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Navin_Visvader_NORM_panCK_Reference.rds")

# Set identities 
colnames(x = Visvader_0025_ER_panCK[[]])
Idents(object = Visvader_0025_ER_panCK) <- 'Tumor'
Visvader_0025_ER_panCK[["cnv.ident"]] <- Idents(object = Visvader_0025_ER_panCK)
levels(x = Visvader_0025_ER_panCK)

colnames(x = Navin_Visvader_NORM_panCK_Reference[[]])
Idents(object = Navin_Visvader_NORM_panCK_Reference) <- 'Normal'
Navin_Visvader_NORM_panCK_Reference[["cnv.ident"]] <- Idents(object = Navin_Visvader_NORM_panCK_Reference)
levels(x = Navin_Visvader_NORM_panCK_Reference)

# Merge subsets into recombined RDS
Visvader_0025_ER_Navin_Visvader_NORM_panCK_Reference <- merge(x = Visvader_0025_ER_panCK, y = Navin_Visvader_NORM_panCK_Reference)
# Join Layers
Visvader_0025_ER_Navin_Visvader_NORM_panCK_Reference <- JoinLayers(Visvader_0025_ER_Navin_Visvader_NORM_panCK_Reference)

# Remove subsetted objects
rm(Visvader_0025_ER_panCK)
rm(Navin_Visvader_NORM_panCK_Reference)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Visvader_0025_ER_Navin_Visvader_NORM_panCK_Reference <- FindVariableFeatures(Visvader_0025_ER_Navin_Visvader_NORM_panCK_Reference, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Visvader_0025_ER_Navin_Visvader_NORM_panCK_Reference)
Visvader_0025_ER_Navin_Visvader_NORM_panCK_Reference <- ScaleData(Visvader_0025_ER_Navin_Visvader_NORM_panCK_Reference, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Visvader_0025_ER_Navin_Visvader_NORM_panCK_Reference <- RunPCA(Visvader_0025_ER_Navin_Visvader_NORM_panCK_Reference, features = VariableFeatures(object = Visvader_0025_ER_Navin_Visvader_NORM_panCK_Reference))

# Examine and visualize PCA results a few different ways
print(Visvader_0025_ER_Navin_Visvader_NORM_panCK_Reference[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Visvader_0025_ER_Navin_Visvader_NORM_panCK_Reference, dims = 1:2, reduction = "pca")
DimPlot(Visvader_0025_ER_Navin_Visvader_NORM_panCK_Reference, reduction = "pca")

Visvader_0025_ER_Navin_Visvader_NORM_panCK_Reference <- FindNeighbors(Visvader_0025_ER_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
Visvader_0025_ER_Navin_Visvader_NORM_panCK_Reference <- FindClusters(Visvader_0025_ER_Navin_Visvader_NORM_panCK_Reference, resolution = 0.1)

# Umap clustering
Visvader_0025_ER_Navin_Visvader_NORM_panCK_Reference <- RunUMAP(Visvader_0025_ER_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
DimPlot(Visvader_0025_ER_Navin_Visvader_NORM_panCK_Reference, reduction = "umap")
DimPlot(Visvader_0025_ER_Navin_Visvader_NORM_panCK_Reference, reduction = "umap", group.by = "orig.ident")

# Save RDS
saveRDS(Visvader_0025_ER_Navin_Visvader_NORM_panCK_Reference, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Visvader_0025_ER_Navin_Visvader_NORM_panCK_Reference.rds")

########## Prepare Cell Annotations
# Extract & Save Meta Data
Visvader_0025_ER_Navin_Visvader_NORM_panCK_Reference.meta.data <- as.data.frame(as.matrix(Visvader_0025_ER_Navin_Visvader_NORM_panCK_Reference@meta.data))
# Extract pertinent columns
# Make row names first column
Visvader_0025_ER_Navin_Visvader_NORM_panCK_Reference.meta.data <- tibble::rownames_to_column(Visvader_0025_ER_Navin_Visvader_NORM_panCK_Reference.meta.data, "Cells")
# In Seurat, "X" becomes the cell names while "cnv.ident" is specified above
Visvader_0025_ER_Navin_Visvader_NORM_panCK_Reference.annotations <- Visvader_0025_ER_Navin_Visvader_NORM_panCK_Reference.meta.data[, c("Cells", "cnv.ident")]
# Delete header
names(Visvader_0025_ER_Navin_Visvader_NORM_panCK_Reference.annotations) <- NULL
# Save as table
write.table(Visvader_0025_ER_Navin_Visvader_NORM_panCK_Reference.annotations,"/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_0025_ER_Navin_Visvader_NORM_panCK_Reference.annotations.txt",sep="\t",row.names=FALSE)
# Remove metadata and annotations
rm(Visvader_0025_ER_Navin_Visvader_NORM_panCK_Reference.meta.data)
rm(Visvader_0025_ER_Navin_Visvader_NORM_panCK_Reference.annotations)

########## Generate Count Matrix
# Extract counts from the Seurat object
Visvader_0025_ER_Navin_Visvader_NORM_panCK_Reference.count <- Visvader_0025_ER_Navin_Visvader_NORM_panCK_Reference[["RNA"]]$counts
# Convert counts to a sparse matrix
Visvader_0025_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as(Visvader_0025_ER_Navin_Visvader_NORM_panCK_Reference.count, 'sparseMatrix')
rm(Visvader_0025_ER_Navin_Visvader_NORM_panCK_Reference.count)

# Convert to requisite data format
Visvader_0025_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as.data.frame(Visvader_0025_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix)
# Save table
write.csv(Visvader_0025_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_Count_Matrix/Visvader_0025_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv")
# Remove Objects
rm(Visvader_0025_ER_Navin_Visvader_NORM_panCK_Reference)
rm(Visvader_0025_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix)

########### Run inferCNV
# Note: check.names = FALSE prevents cell names from having dashes replaced with dots; essential to match annotations file
Visvader_0025_ER_inferCNV <- CreateInfercnvObject(raw_counts_matrix=read.csv(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_Count_Matrix/Visvader_0025_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv", header = TRUE, row.names = 1,check.names = FALSE),
                                                  annotations_file=read.table(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_0025_ER_Navin_Visvader_NORM_panCK_Reference.annotations.txt", header = FALSE, row.names = 1),
                                                  delim="\t",
                                                  gene_order_file="/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/hg38_gencode_v27.txt",
                                                  ref_group_names=c("Normal"))


# perform infercnv operations to reveal cnv signal
Visvader_0025_ER_inferCNV = infercnv::run(Visvader_0025_ER_inferCNV,
                                          cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                          out_dir="/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0025_ER_panCK/",  # dir is auto-created for storing outputs
                                          cluster_by_groups=T,   # cluster
                                          denoise=T,
                                          HMM=T,
                                          num_ref_groups = 8, analysis_mode = "subclusters")

###########  Integrate data with Seurat object
# Read RDS and add inferCNV results to metadata
Visvader_0025_ER_Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Visvader_0025_ER_Navin_Visvader_NORM_panCK_Reference.rds") 
Visvader_0025_ER_panCK_cnv <- infercnv::add_to_seurat(infercnv_output_path="/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0025_ER_panCK/",
                                                      seurat_obj=Visvader_0025_ER_Navin_Visvader_NORM_panCK_Reference, # optional
                                                      top_n=10)

# Collect Tumor cells from Seurat Object
Visvader_0025_ER_panCK_cnv <- subset(x = Visvader_0025_ER_panCK_cnv, subset = cnv.ident == "Tumor")
# Quantify chromosomes with alterations per cell
Visvader_0025_ER_panCK_cnv_Pos_Chromosomes <- as.data.frame(rowSums(Visvader_0025_ER_panCK_cnv@meta.data[,startsWith(names(Visvader_0025_ER_panCK_cnv@meta.data),"has_cnv_")] > 0 ))
Visvader_0025_ER_panCK_cnv_Neg_Chromosomes <- as.data.frame(rowSums(Visvader_0025_ER_panCK_cnv@meta.data[,startsWith(names(Visvader_0025_ER_panCK_cnv@meta.data),"has_cnv_")] == 0))
Visvader_0025_ER_panCK_cnv_Chromosome_Quant <- as.data.frame(c(Visvader_0025_ER_panCK_cnv_Pos_Chromosomes, Visvader_0025_ER_panCK_cnv_Neg_Chromosomes))
colnames(Visvader_0025_ER_panCK_cnv_Chromosome_Quant) <- c("Chromosomes CNV Pos", "Chromosomes CNV Neg")
write.csv(Visvader_0025_ER_panCK_cnv_Chromosome_Quant, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0025_ER_panCK/Visvader_0025_ER_panCK_cnv_Chromosome_Quant.csv")

# Save CNV Metadata
Visvader_0025_ER_panCK_cnv.meta.data <- as.data.frame(as.matrix(Visvader_0025_ER_panCK_cnv@meta.data))
write.csv(Visvader_0025_ER_panCK_cnv.meta.data, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0025_ER_panCK/Visvader_0025_ER_panCK_cnv.meta.data.csv")

# Subset Cells with inferred CNVs -------------------------------------------------------------------------------------------------------------
Visvader_0025_ER_panCK_cnv <- subset(Visvader_0025_ER_panCK_cnv, cells=rownames(Visvader_0025_ER_panCK_cnv@meta.data)[rowSums(Visvader_0025_ER_panCK_cnv@meta.data[,startsWith(names(Visvader_0025_ER_panCK_cnv@meta.data),"has_cnv_")]) > 0 ])

# Save RDS with inferredCNVs
saveRDS(Visvader_0025_ER_panCK_cnv, file = "/R/R_Visvader/Visvader_RDS/RDS_inferCNV/Visvader_0025_ER_panCK_cnv.rds")


# Free memory
rm(Visvader_0025_ER_inferCNV)
rm(Visvader_0025_ER_Navin_Visvader_NORM_panCK_Reference)
rm(Visvader_0025_ER_panCK_cnv)
rm(Visvader_0025_ER_panCK_cnv_Pos_Chromosomes)
rm(Visvader_0025_ER_panCK_cnv_Neg_Chromosomes)
rm(Visvader_0025_ER_panCK_cnv_Chromosome_Quant)
rm(Visvader_0025_ER_panCK_cnv.meta.data)
gc()



#############################################
#### InferCNV: Visvader_0029_7C_ER_panCK ####
#############################################

########## Merge with normal mammary glad, extract epithelium, and process via Seurat
# Load RDS
Visvader_0029_7C_ER_panCK <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0029_7C_ER_panCK.rds")
Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Navin_Visvader_NORM_panCK_Reference.rds")

# Set identities 
colnames(x = Visvader_0029_7C_ER_panCK[[]])
Idents(object = Visvader_0029_7C_ER_panCK) <- 'Tumor'
Visvader_0029_7C_ER_panCK[["cnv.ident"]] <- Idents(object = Visvader_0029_7C_ER_panCK)
levels(x = Visvader_0029_7C_ER_panCK)

colnames(x = Navin_Visvader_NORM_panCK_Reference[[]])
Idents(object = Navin_Visvader_NORM_panCK_Reference) <- 'Normal'
Navin_Visvader_NORM_panCK_Reference[["cnv.ident"]] <- Idents(object = Navin_Visvader_NORM_panCK_Reference)
levels(x = Navin_Visvader_NORM_panCK_Reference)

# Merge subsets into recombined RDS
Visvader_0029_7C_ER_Navin_Visvader_NORM_panCK_Reference <- merge(x = Visvader_0029_7C_ER_panCK, y = Navin_Visvader_NORM_panCK_Reference)
# Join Layers
Visvader_0029_7C_ER_Navin_Visvader_NORM_panCK_Reference <- JoinLayers(Visvader_0029_7C_ER_Navin_Visvader_NORM_panCK_Reference)

# Remove subsetted objects
rm(Visvader_0029_7C_ER_panCK)
rm(Navin_Visvader_NORM_panCK_Reference)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Visvader_0029_7C_ER_Navin_Visvader_NORM_panCK_Reference <- FindVariableFeatures(Visvader_0029_7C_ER_Navin_Visvader_NORM_panCK_Reference, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Visvader_0029_7C_ER_Navin_Visvader_NORM_panCK_Reference)
Visvader_0029_7C_ER_Navin_Visvader_NORM_panCK_Reference <- ScaleData(Visvader_0029_7C_ER_Navin_Visvader_NORM_panCK_Reference, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Visvader_0029_7C_ER_Navin_Visvader_NORM_panCK_Reference <- RunPCA(Visvader_0029_7C_ER_Navin_Visvader_NORM_panCK_Reference, features = VariableFeatures(object = Visvader_0029_7C_ER_Navin_Visvader_NORM_panCK_Reference))

# Examine and visualize PCA results a few different ways
print(Visvader_0029_7C_ER_Navin_Visvader_NORM_panCK_Reference[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Visvader_0029_7C_ER_Navin_Visvader_NORM_panCK_Reference, dims = 1:2, reduction = "pca")
DimPlot(Visvader_0029_7C_ER_Navin_Visvader_NORM_panCK_Reference, reduction = "pca")

Visvader_0029_7C_ER_Navin_Visvader_NORM_panCK_Reference <- FindNeighbors(Visvader_0029_7C_ER_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
Visvader_0029_7C_ER_Navin_Visvader_NORM_panCK_Reference <- FindClusters(Visvader_0029_7C_ER_Navin_Visvader_NORM_panCK_Reference, resolution = 0.1)

# Umap clustering
Visvader_0029_7C_ER_Navin_Visvader_NORM_panCK_Reference <- RunUMAP(Visvader_0029_7C_ER_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
DimPlot(Visvader_0029_7C_ER_Navin_Visvader_NORM_panCK_Reference, reduction = "umap")
DimPlot(Visvader_0029_7C_ER_Navin_Visvader_NORM_panCK_Reference, reduction = "umap", group.by = "orig.ident")

# Save RDS
saveRDS(Visvader_0029_7C_ER_Navin_Visvader_NORM_panCK_Reference, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Visvader_0029_7C_ER_Navin_Visvader_NORM_panCK_Reference.rds")

########## Prepare Cell Annotations
# Extract & Save Meta Data
Visvader_0029_7C_ER_Navin_Visvader_NORM_panCK_Reference.meta.data <- as.data.frame(as.matrix(Visvader_0029_7C_ER_Navin_Visvader_NORM_panCK_Reference@meta.data))
# Extract pertinent columns
# Make row names first column
Visvader_0029_7C_ER_Navin_Visvader_NORM_panCK_Reference.meta.data <- tibble::rownames_to_column(Visvader_0029_7C_ER_Navin_Visvader_NORM_panCK_Reference.meta.data, "Cells")
# In Seurat, "X" becomes the cell names while "cnv.ident" is specified above
Visvader_0029_7C_ER_Navin_Visvader_NORM_panCK_Reference.annotations <- Visvader_0029_7C_ER_Navin_Visvader_NORM_panCK_Reference.meta.data[, c("Cells", "cnv.ident")]
# Delete header
names(Visvader_0029_7C_ER_Navin_Visvader_NORM_panCK_Reference.annotations) <- NULL
# Save as table
write.table(Visvader_0029_7C_ER_Navin_Visvader_NORM_panCK_Reference.annotations,"/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_0029_7C_ER_Navin_Visvader_NORM_panCK_Reference.annotations.txt",sep="\t",row.names=FALSE)
# Remove metadata and annotations
rm(Visvader_0029_7C_ER_Navin_Visvader_NORM_panCK_Reference.meta.data)
rm(Visvader_0029_7C_ER_Navin_Visvader_NORM_panCK_Reference.annotations)

########## Generate Count Matrix
# Extract counts from the Seurat object
Visvader_0029_7C_ER_Navin_Visvader_NORM_panCK_Reference.count <- Visvader_0029_7C_ER_Navin_Visvader_NORM_panCK_Reference[["RNA"]]$counts
# Convert counts to a sparse matrix
Visvader_0029_7C_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as(Visvader_0029_7C_ER_Navin_Visvader_NORM_panCK_Reference.count, 'sparseMatrix')
rm(Visvader_0029_7C_ER_Navin_Visvader_NORM_panCK_Reference.count)

# Convert to requisite data format
Visvader_0029_7C_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as.data.frame(Visvader_0029_7C_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix)
# Save table
write.csv(Visvader_0029_7C_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_Count_Matrix/Visvader_0029_7C_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv")
# Remove Objects
rm(Visvader_0029_7C_ER_Navin_Visvader_NORM_panCK_Reference)
rm(Visvader_0029_7C_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix)

########### Run inferCNV
# Note: check.names = FALSE prevents cell names from having dashes replaced with dots; essential to match annotations file
Visvader_0029_7C_ER_inferCNV <- CreateInfercnvObject(raw_counts_matrix=read.csv(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_Count_Matrix/Visvader_0029_7C_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv", header = TRUE, row.names = 1,check.names = FALSE),
                                                     annotations_file=read.table(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_0029_7C_ER_Navin_Visvader_NORM_panCK_Reference.annotations.txt", header = FALSE, row.names = 1),
                                                     delim="\t",
                                                     gene_order_file="/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/hg38_gencode_v27.txt",
                                                     ref_group_names=c("Normal"))


# perform infercnv operations to reveal cnv signal
Visvader_0029_7C_ER_inferCNV = infercnv::run(Visvader_0029_7C_ER_inferCNV,
                                             cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                             out_dir="/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0029_7C_ER_panCK/",  # dir is auto-created for storing outputs
                                             cluster_by_groups=T,   # cluster
                                             denoise=T,
                                             HMM=T,
                                             num_ref_groups = 8, analysis_mode = "subclusters")

###########  Integrate data with Seurat object
# Read RDS and add inferCNV results to metadata
Visvader_0029_7C_ER_Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Visvader_0029_7C_ER_Navin_Visvader_NORM_panCK_Reference.rds") 
Visvader_0029_7C_ER_panCK_cnv <- infercnv::add_to_seurat(infercnv_output_path="/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0029_7C_ER_panCK/",
                                                         seurat_obj=Visvader_0029_7C_ER_Navin_Visvader_NORM_panCK_Reference, # optional
                                                         top_n=10)

# Collect Tumor cells from Seurat Object
Visvader_0029_7C_ER_panCK_cnv <- subset(x = Visvader_0029_7C_ER_panCK_cnv, subset = cnv.ident == "Tumor")
# Quantify chromosomes with alterations per cell
Visvader_0029_7C_ER_panCK_cnv_Pos_Chromosomes <- as.data.frame(rowSums(Visvader_0029_7C_ER_panCK_cnv@meta.data[,startsWith(names(Visvader_0029_7C_ER_panCK_cnv@meta.data),"has_cnv_")] > 0 ))
Visvader_0029_7C_ER_panCK_cnv_Neg_Chromosomes <- as.data.frame(rowSums(Visvader_0029_7C_ER_panCK_cnv@meta.data[,startsWith(names(Visvader_0029_7C_ER_panCK_cnv@meta.data),"has_cnv_")] == 0))
Visvader_0029_7C_ER_panCK_cnv_Chromosome_Quant <- as.data.frame(c(Visvader_0029_7C_ER_panCK_cnv_Pos_Chromosomes, Visvader_0029_7C_ER_panCK_cnv_Neg_Chromosomes))
colnames(Visvader_0029_7C_ER_panCK_cnv_Chromosome_Quant) <- c("Chromosomes CNV Pos", "Chromosomes CNV Neg")
write.csv(Visvader_0029_7C_ER_panCK_cnv_Chromosome_Quant, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0029_7C_ER_panCK/Visvader_0029_7C_ER_panCK_cnv_Chromosome_Quant.csv")

# Save CNV Metadata
Visvader_0029_7C_ER_panCK_cnv.meta.data <- as.data.frame(as.matrix(Visvader_0029_7C_ER_panCK_cnv@meta.data))
write.csv(Visvader_0029_7C_ER_panCK_cnv.meta.data, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0029_7C_ER_panCK/Visvader_0029_7C_ER_panCK_cnv.meta.data.csv")

# Subset Cells with inferred CNVs -------------------------------------------------------------------------------------------------------------
Visvader_0029_7C_ER_panCK_cnv <- subset(Visvader_0029_7C_ER_panCK_cnv, cells=rownames(Visvader_0029_7C_ER_panCK_cnv@meta.data)[rowSums(Visvader_0029_7C_ER_panCK_cnv@meta.data[,startsWith(names(Visvader_0029_7C_ER_panCK_cnv@meta.data),"has_cnv_")]) > 0 ])

# Save RDS with inferredCNVs
saveRDS(Visvader_0029_7C_ER_panCK_cnv, file = "/R/R_Visvader/Visvader_RDS/RDS_inferCNV/Visvader_0029_7C_ER_panCK_cnv.rds")


# Free memory
rm(Visvader_0029_7C_ER_inferCNV)
rm(Visvader_0029_7C_ER_Navin_Visvader_NORM_panCK_Reference)
rm(Visvader_0029_7C_ER_panCK_cnv)
rm(Visvader_0029_7C_ER_panCK_cnv_Pos_Chromosomes)
rm(Visvader_0029_7C_ER_panCK_cnv_Neg_Chromosomes)
rm(Visvader_0029_7C_ER_panCK_cnv_Chromosome_Quant)
rm(Visvader_0029_7C_ER_panCK_cnv.meta.data)
gc()



#############################################
#### InferCNV: Visvader_0029_9C_ER_panCK ####
#############################################

########## Merge with normal mammary glad, extract epithelium, and process via Seurat
# Load RDS
Visvader_0029_9C_ER_panCK <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0029_9C_ER_panCK.rds")
Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Navin_Visvader_NORM_panCK_Reference.rds")

# Set identities 
colnames(x = Visvader_0029_9C_ER_panCK[[]])
Idents(object = Visvader_0029_9C_ER_panCK) <- 'Tumor'
Visvader_0029_9C_ER_panCK[["cnv.ident"]] <- Idents(object = Visvader_0029_9C_ER_panCK)
levels(x = Visvader_0029_9C_ER_panCK)

colnames(x = Navin_Visvader_NORM_panCK_Reference[[]])
Idents(object = Navin_Visvader_NORM_panCK_Reference) <- 'Normal'
Navin_Visvader_NORM_panCK_Reference[["cnv.ident"]] <- Idents(object = Navin_Visvader_NORM_panCK_Reference)
levels(x = Navin_Visvader_NORM_panCK_Reference)

# Merge subsets into recombined RDS
Visvader_0029_9C_ER_Navin_Visvader_NORM_panCK_Reference <- merge(x = Visvader_0029_9C_ER_panCK, y = Navin_Visvader_NORM_panCK_Reference)
# Join Layers
Visvader_0029_9C_ER_Navin_Visvader_NORM_panCK_Reference <- JoinLayers(Visvader_0029_9C_ER_Navin_Visvader_NORM_panCK_Reference)

# Remove subsetted objects
rm(Visvader_0029_9C_ER_panCK)
rm(Navin_Visvader_NORM_panCK_Reference)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Visvader_0029_9C_ER_Navin_Visvader_NORM_panCK_Reference <- FindVariableFeatures(Visvader_0029_9C_ER_Navin_Visvader_NORM_panCK_Reference, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Visvader_0029_9C_ER_Navin_Visvader_NORM_panCK_Reference)
Visvader_0029_9C_ER_Navin_Visvader_NORM_panCK_Reference <- ScaleData(Visvader_0029_9C_ER_Navin_Visvader_NORM_panCK_Reference, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Visvader_0029_9C_ER_Navin_Visvader_NORM_panCK_Reference <- RunPCA(Visvader_0029_9C_ER_Navin_Visvader_NORM_panCK_Reference, features = VariableFeatures(object = Visvader_0029_9C_ER_Navin_Visvader_NORM_panCK_Reference))

# Examine and visualize PCA results a few different ways
print(Visvader_0029_9C_ER_Navin_Visvader_NORM_panCK_Reference[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Visvader_0029_9C_ER_Navin_Visvader_NORM_panCK_Reference, dims = 1:2, reduction = "pca")
DimPlot(Visvader_0029_9C_ER_Navin_Visvader_NORM_panCK_Reference, reduction = "pca")

Visvader_0029_9C_ER_Navin_Visvader_NORM_panCK_Reference <- FindNeighbors(Visvader_0029_9C_ER_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
Visvader_0029_9C_ER_Navin_Visvader_NORM_panCK_Reference <- FindClusters(Visvader_0029_9C_ER_Navin_Visvader_NORM_panCK_Reference, resolution = 0.1)

# Umap clustering
Visvader_0029_9C_ER_Navin_Visvader_NORM_panCK_Reference <- RunUMAP(Visvader_0029_9C_ER_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
DimPlot(Visvader_0029_9C_ER_Navin_Visvader_NORM_panCK_Reference, reduction = "umap")
DimPlot(Visvader_0029_9C_ER_Navin_Visvader_NORM_panCK_Reference, reduction = "umap", group.by = "orig.ident")

# Save RDS
saveRDS(Visvader_0029_9C_ER_Navin_Visvader_NORM_panCK_Reference, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Visvader_0029_9C_ER_Navin_Visvader_NORM_panCK_Reference.rds")

########## Prepare Cell Annotations
# Extract & Save Meta Data
Visvader_0029_9C_ER_Navin_Visvader_NORM_panCK_Reference.meta.data <- as.data.frame(as.matrix(Visvader_0029_9C_ER_Navin_Visvader_NORM_panCK_Reference@meta.data))
# Extract pertinent columns
# Make row names first column
Visvader_0029_9C_ER_Navin_Visvader_NORM_panCK_Reference.meta.data <- tibble::rownames_to_column(Visvader_0029_9C_ER_Navin_Visvader_NORM_panCK_Reference.meta.data, "Cells")
# In Seurat, "X" becomes the cell names while "cnv.ident" is specified above
Visvader_0029_9C_ER_Navin_Visvader_NORM_panCK_Reference.annotations <- Visvader_0029_9C_ER_Navin_Visvader_NORM_panCK_Reference.meta.data[, c("Cells", "cnv.ident")]
# Delete header
names(Visvader_0029_9C_ER_Navin_Visvader_NORM_panCK_Reference.annotations) <- NULL
# Save as table
write.table(Visvader_0029_9C_ER_Navin_Visvader_NORM_panCK_Reference.annotations,"/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_0029_9C_ER_Navin_Visvader_NORM_panCK_Reference.annotations.txt",sep="\t",row.names=FALSE)
# Remove metadata and annotations
rm(Visvader_0029_9C_ER_Navin_Visvader_NORM_panCK_Reference.meta.data)
rm(Visvader_0029_9C_ER_Navin_Visvader_NORM_panCK_Reference.annotations)

########## Generate Count Matrix
# Extract counts from the Seurat object
Visvader_0029_9C_ER_Navin_Visvader_NORM_panCK_Reference.count <- Visvader_0029_9C_ER_Navin_Visvader_NORM_panCK_Reference[["RNA"]]$counts
# Convert counts to a sparse matrix
Visvader_0029_9C_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as(Visvader_0029_9C_ER_Navin_Visvader_NORM_panCK_Reference.count, 'sparseMatrix')
rm(Visvader_0029_9C_ER_Navin_Visvader_NORM_panCK_Reference.count)

# Convert to requisite data format
Visvader_0029_9C_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as.data.frame(Visvader_0029_9C_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix)
# Save table
write.csv(Visvader_0029_9C_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_Count_Matrix/Visvader_0029_9C_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv")
# Remove Objects
rm(Visvader_0029_9C_ER_Navin_Visvader_NORM_panCK_Reference)
rm(Visvader_0029_9C_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix)

########### Run inferCNV
# Note: check.names = FALSE prevents cell names from having dashes replaced with dots; essential to match annotations file
Visvader_0029_9C_ER_inferCNV <- CreateInfercnvObject(raw_counts_matrix=read.csv(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_Count_Matrix/Visvader_0029_9C_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv", header = TRUE, row.names = 1,check.names = FALSE),
                                                     annotations_file=read.table(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_0029_9C_ER_Navin_Visvader_NORM_panCK_Reference.annotations.txt", header = FALSE, row.names = 1),
                                                     delim="\t",
                                                     gene_order_file="/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/hg38_gencode_v27.txt",
                                                     ref_group_names=c("Normal"))


# perform infercnv operations to reveal cnv signal
Visvader_0029_9C_ER_inferCNV = infercnv::run(Visvader_0029_9C_ER_inferCNV,
                                             cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                             out_dir="/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0029_9C_ER_panCK/",  # dir is auto-created for storing outputs
                                             cluster_by_groups=T,   # cluster
                                             denoise=T,
                                             HMM=T,
                                             num_ref_groups = 8, analysis_mode = "subclusters")

###########  Integrate data with Seurat object
# Read RDS and add inferCNV results to metadata
Visvader_0029_9C_ER_Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Visvader_0029_9C_ER_Navin_Visvader_NORM_panCK_Reference.rds") 
Visvader_0029_9C_ER_panCK_cnv <- infercnv::add_to_seurat(infercnv_output_path="/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0029_9C_ER_panCK/",
                                                         seurat_obj=Visvader_0029_9C_ER_Navin_Visvader_NORM_panCK_Reference, # optional
                                                         top_n=10)

# Collect Tumor cells from Seurat Object
Visvader_0029_9C_ER_panCK_cnv <- subset(x = Visvader_0029_9C_ER_panCK_cnv, subset = cnv.ident == "Tumor")
# Quantify chromosomes with alterations per cell
Visvader_0029_9C_ER_panCK_cnv_Pos_Chromosomes <- as.data.frame(rowSums(Visvader_0029_9C_ER_panCK_cnv@meta.data[,startsWith(names(Visvader_0029_9C_ER_panCK_cnv@meta.data),"has_cnv_")] > 0 ))
Visvader_0029_9C_ER_panCK_cnv_Neg_Chromosomes <- as.data.frame(rowSums(Visvader_0029_9C_ER_panCK_cnv@meta.data[,startsWith(names(Visvader_0029_9C_ER_panCK_cnv@meta.data),"has_cnv_")] == 0))
Visvader_0029_9C_ER_panCK_cnv_Chromosome_Quant <- as.data.frame(c(Visvader_0029_9C_ER_panCK_cnv_Pos_Chromosomes, Visvader_0029_9C_ER_panCK_cnv_Neg_Chromosomes))
colnames(Visvader_0029_9C_ER_panCK_cnv_Chromosome_Quant) <- c("Chromosomes CNV Pos", "Chromosomes CNV Neg")
write.csv(Visvader_0029_9C_ER_panCK_cnv_Chromosome_Quant, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0029_9C_ER_panCK/Visvader_0029_9C_ER_panCK_cnv_Chromosome_Quant.csv")

# Save CNV Metadata
Visvader_0029_9C_ER_panCK_cnv.meta.data <- as.data.frame(as.matrix(Visvader_0029_9C_ER_panCK_cnv@meta.data))
write.csv(Visvader_0029_9C_ER_panCK_cnv.meta.data, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0029_9C_ER_panCK/Visvader_0029_9C_ER_panCK_cnv.meta.data.csv")

# Subset Cells with inferred CNVs -------------------------------------------------------------------------------------------------------------
Visvader_0029_9C_ER_panCK_cnv <- subset(Visvader_0029_9C_ER_panCK_cnv, cells=rownames(Visvader_0029_9C_ER_panCK_cnv@meta.data)[rowSums(Visvader_0029_9C_ER_panCK_cnv@meta.data[,startsWith(names(Visvader_0029_9C_ER_panCK_cnv@meta.data),"has_cnv_")]) > 0 ])

# Save RDS with inferredCNVs
saveRDS(Visvader_0029_9C_ER_panCK_cnv, file = "/R/R_Visvader/Visvader_RDS/RDS_inferCNV/Visvader_0029_9C_ER_panCK_cnv.rds")


# Free memory
rm(Visvader_0029_9C_ER_inferCNV)
rm(Visvader_0029_9C_ER_Navin_Visvader_NORM_panCK_Reference)
rm(Visvader_0029_9C_ER_panCK_cnv)
rm(Visvader_0029_9C_ER_panCK_cnv_Pos_Chromosomes)
rm(Visvader_0029_9C_ER_panCK_cnv_Neg_Chromosomes)
rm(Visvader_0029_9C_ER_panCK_cnv_Chromosome_Quant)
rm(Visvader_0029_9C_ER_panCK_cnv.meta.data)
gc()



##########################################
#### InferCNV: Visvader_0032_ER_panCK ####
##########################################

########## Merge with normal mammary glad, extract epithelium, and process via Seurat
# Load RDS
Visvader_0032_ER_panCK <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0032_ER_panCK.rds")
Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Navin_Visvader_NORM_panCK_Reference.rds")

# Set identities 
colnames(x = Visvader_0032_ER_panCK[[]])
Idents(object = Visvader_0032_ER_panCK) <- 'Tumor'
Visvader_0032_ER_panCK[["cnv.ident"]] <- Idents(object = Visvader_0032_ER_panCK)
levels(x = Visvader_0032_ER_panCK)

colnames(x = Navin_Visvader_NORM_panCK_Reference[[]])
Idents(object = Navin_Visvader_NORM_panCK_Reference) <- 'Normal'
Navin_Visvader_NORM_panCK_Reference[["cnv.ident"]] <- Idents(object = Navin_Visvader_NORM_panCK_Reference)
levels(x = Navin_Visvader_NORM_panCK_Reference)

# Merge subsets into recombined RDS
Visvader_0032_ER_Navin_Visvader_NORM_panCK_Reference <- merge(x = Visvader_0032_ER_panCK, y = Navin_Visvader_NORM_panCK_Reference)
# Join Layers
Visvader_0032_ER_Navin_Visvader_NORM_panCK_Reference <- JoinLayers(Visvader_0032_ER_Navin_Visvader_NORM_panCK_Reference)

# Remove subsetted objects
rm(Visvader_0032_ER_panCK)
rm(Navin_Visvader_NORM_panCK_Reference)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Visvader_0032_ER_Navin_Visvader_NORM_panCK_Reference <- FindVariableFeatures(Visvader_0032_ER_Navin_Visvader_NORM_panCK_Reference, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Visvader_0032_ER_Navin_Visvader_NORM_panCK_Reference)
Visvader_0032_ER_Navin_Visvader_NORM_panCK_Reference <- ScaleData(Visvader_0032_ER_Navin_Visvader_NORM_panCK_Reference, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Visvader_0032_ER_Navin_Visvader_NORM_panCK_Reference <- RunPCA(Visvader_0032_ER_Navin_Visvader_NORM_panCK_Reference, features = VariableFeatures(object = Visvader_0032_ER_Navin_Visvader_NORM_panCK_Reference))

# Examine and visualize PCA results a few different ways
print(Visvader_0032_ER_Navin_Visvader_NORM_panCK_Reference[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Visvader_0032_ER_Navin_Visvader_NORM_panCK_Reference, dims = 1:2, reduction = "pca")
DimPlot(Visvader_0032_ER_Navin_Visvader_NORM_panCK_Reference, reduction = "pca")

Visvader_0032_ER_Navin_Visvader_NORM_panCK_Reference <- FindNeighbors(Visvader_0032_ER_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
Visvader_0032_ER_Navin_Visvader_NORM_panCK_Reference <- FindClusters(Visvader_0032_ER_Navin_Visvader_NORM_panCK_Reference, resolution = 0.1)

# Umap clustering
Visvader_0032_ER_Navin_Visvader_NORM_panCK_Reference <- RunUMAP(Visvader_0032_ER_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
DimPlot(Visvader_0032_ER_Navin_Visvader_NORM_panCK_Reference, reduction = "umap")
DimPlot(Visvader_0032_ER_Navin_Visvader_NORM_panCK_Reference, reduction = "umap", group.by = "orig.ident")

# Save RDS
saveRDS(Visvader_0032_ER_Navin_Visvader_NORM_panCK_Reference, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Visvader_0032_ER_Navin_Visvader_NORM_panCK_Reference.rds")

########## Prepare Cell Annotations
# Extract & Save Meta Data
Visvader_0032_ER_Navin_Visvader_NORM_panCK_Reference.meta.data <- as.data.frame(as.matrix(Visvader_0032_ER_Navin_Visvader_NORM_panCK_Reference@meta.data))
# Extract pertinent columns
# Make row names first column
Visvader_0032_ER_Navin_Visvader_NORM_panCK_Reference.meta.data <- tibble::rownames_to_column(Visvader_0032_ER_Navin_Visvader_NORM_panCK_Reference.meta.data, "Cells")
# In Seurat, "X" becomes the cell names while "cnv.ident" is specified above
Visvader_0032_ER_Navin_Visvader_NORM_panCK_Reference.annotations <- Visvader_0032_ER_Navin_Visvader_NORM_panCK_Reference.meta.data[, c("Cells", "cnv.ident")]
# Delete header
names(Visvader_0032_ER_Navin_Visvader_NORM_panCK_Reference.annotations) <- NULL
# Save as table
write.table(Visvader_0032_ER_Navin_Visvader_NORM_panCK_Reference.annotations,"/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_0032_ER_Navin_Visvader_NORM_panCK_Reference.annotations.txt",sep="\t",row.names=FALSE)
# Remove metadata and annotations
rm(Visvader_0032_ER_Navin_Visvader_NORM_panCK_Reference.meta.data)
rm(Visvader_0032_ER_Navin_Visvader_NORM_panCK_Reference.annotations)

########## Generate Count Matrix
# Extract counts from the Seurat object
Visvader_0032_ER_Navin_Visvader_NORM_panCK_Reference.count <- Visvader_0032_ER_Navin_Visvader_NORM_panCK_Reference[["RNA"]]$counts
# Convert counts to a sparse matrix
Visvader_0032_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as(Visvader_0032_ER_Navin_Visvader_NORM_panCK_Reference.count, 'sparseMatrix')
rm(Visvader_0032_ER_Navin_Visvader_NORM_panCK_Reference.count)

# Convert to requisite data format
Visvader_0032_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as.data.frame(Visvader_0032_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix)
# Save table
write.csv(Visvader_0032_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_Count_Matrix/Visvader_0032_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv")
# Remove Objects
rm(Visvader_0032_ER_Navin_Visvader_NORM_panCK_Reference)
rm(Visvader_0032_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix)

########### Run inferCNV
# Note: check.names = FALSE prevents cell names from having dashes replaced with dots; essential to match annotations file
Visvader_0032_ER_inferCNV <- CreateInfercnvObject(raw_counts_matrix=read.csv(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_Count_Matrix/Visvader_0032_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv", header = TRUE, row.names = 1,check.names = FALSE),
                                                  annotations_file=read.table(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_0032_ER_Navin_Visvader_NORM_panCK_Reference.annotations.txt", header = FALSE, row.names = 1),
                                                  delim="\t",
                                                  gene_order_file="/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/hg38_gencode_v27.txt",
                                                  ref_group_names=c("Normal"))


# perform infercnv operations to reveal cnv signal
Visvader_0032_ER_inferCNV = infercnv::run(Visvader_0032_ER_inferCNV,
                                          cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                          out_dir="/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0032_ER_panCK/",  # dir is auto-created for storing outputs
                                          cluster_by_groups=T,   # cluster
                                          denoise=T,
                                          HMM=T,
                                          num_ref_groups = 8, analysis_mode = "subclusters")

###########  Integrate data with Seurat object
# Read RDS and add inferCNV results to metadata
Visvader_0032_ER_Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Visvader_0032_ER_Navin_Visvader_NORM_panCK_Reference.rds") 
Visvader_0032_ER_panCK_cnv <- infercnv::add_to_seurat(infercnv_output_path="/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0032_ER_panCK/",
                                                      seurat_obj=Visvader_0032_ER_Navin_Visvader_NORM_panCK_Reference, # optional
                                                      top_n=10)

# Collect Tumor cells from Seurat Object
Visvader_0032_ER_panCK_cnv <- subset(x = Visvader_0032_ER_panCK_cnv, subset = cnv.ident == "Tumor")
# Quantify chromosomes with alterations per cell
Visvader_0032_ER_panCK_cnv_Pos_Chromosomes <- as.data.frame(rowSums(Visvader_0032_ER_panCK_cnv@meta.data[,startsWith(names(Visvader_0032_ER_panCK_cnv@meta.data),"has_cnv_")] > 0 ))
Visvader_0032_ER_panCK_cnv_Neg_Chromosomes <- as.data.frame(rowSums(Visvader_0032_ER_panCK_cnv@meta.data[,startsWith(names(Visvader_0032_ER_panCK_cnv@meta.data),"has_cnv_")] == 0))
Visvader_0032_ER_panCK_cnv_Chromosome_Quant <- as.data.frame(c(Visvader_0032_ER_panCK_cnv_Pos_Chromosomes, Visvader_0032_ER_panCK_cnv_Neg_Chromosomes))
colnames(Visvader_0032_ER_panCK_cnv_Chromosome_Quant) <- c("Chromosomes CNV Pos", "Chromosomes CNV Neg")
write.csv(Visvader_0032_ER_panCK_cnv_Chromosome_Quant, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0032_ER_panCK/Visvader_0032_ER_panCK_cnv_Chromosome_Quant.csv")

# Save CNV Metadata
Visvader_0032_ER_panCK_cnv.meta.data <- as.data.frame(as.matrix(Visvader_0032_ER_panCK_cnv@meta.data))
write.csv(Visvader_0032_ER_panCK_cnv.meta.data, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0032_ER_panCK/Visvader_0032_ER_panCK_cnv.meta.data.csv")

# Subset Cells with inferred CNVs -------------------------------------------------------------------------------------------------------------
Visvader_0032_ER_panCK_cnv <- subset(Visvader_0032_ER_panCK_cnv, cells=rownames(Visvader_0032_ER_panCK_cnv@meta.data)[rowSums(Visvader_0032_ER_panCK_cnv@meta.data[,startsWith(names(Visvader_0032_ER_panCK_cnv@meta.data),"has_cnv_")]) > 0 ])

# Save RDS with inferredCNVs
saveRDS(Visvader_0032_ER_panCK_cnv, file = "/R/R_Visvader/Visvader_RDS/RDS_inferCNV/Visvader_0032_ER_panCK_cnv.rds")


# Free memory
rm(Visvader_0032_ER_inferCNV)
rm(Visvader_0032_ER_Navin_Visvader_NORM_panCK_Reference)
rm(Visvader_0032_ER_panCK_cnv)
rm(Visvader_0032_ER_panCK_cnv_Pos_Chromosomes)
rm(Visvader_0032_ER_panCK_cnv_Neg_Chromosomes)
rm(Visvader_0032_ER_panCK_cnv_Chromosome_Quant)
rm(Visvader_0032_ER_panCK_cnv.meta.data)
gc()



##########################################
#### InferCNV: Visvader_0040_ER_panCK ####
##########################################

########## Merge with normal mammary glad, extract epithelium, and process via Seurat
# Load RDS
Visvader_0040_ER_panCK <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0040_ER_panCK.rds")
Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Navin_Visvader_NORM_panCK_Reference.rds")

# Set identities 
colnames(x = Visvader_0040_ER_panCK[[]])
Idents(object = Visvader_0040_ER_panCK) <- 'Tumor'
Visvader_0040_ER_panCK[["cnv.ident"]] <- Idents(object = Visvader_0040_ER_panCK)
levels(x = Visvader_0040_ER_panCK)

colnames(x = Navin_Visvader_NORM_panCK_Reference[[]])
Idents(object = Navin_Visvader_NORM_panCK_Reference) <- 'Normal'
Navin_Visvader_NORM_panCK_Reference[["cnv.ident"]] <- Idents(object = Navin_Visvader_NORM_panCK_Reference)
levels(x = Navin_Visvader_NORM_panCK_Reference)

# Merge subsets into recombined RDS
Visvader_0040_ER_Navin_Visvader_NORM_panCK_Reference <- merge(x = Visvader_0040_ER_panCK, y = Navin_Visvader_NORM_panCK_Reference)
# Join Layers
Visvader_0040_ER_Navin_Visvader_NORM_panCK_Reference <- JoinLayers(Visvader_0040_ER_Navin_Visvader_NORM_panCK_Reference)

# Remove subsetted objects
rm(Visvader_0040_ER_panCK)
rm(Navin_Visvader_NORM_panCK_Reference)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Visvader_0040_ER_Navin_Visvader_NORM_panCK_Reference <- FindVariableFeatures(Visvader_0040_ER_Navin_Visvader_NORM_panCK_Reference, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Visvader_0040_ER_Navin_Visvader_NORM_panCK_Reference)
Visvader_0040_ER_Navin_Visvader_NORM_panCK_Reference <- ScaleData(Visvader_0040_ER_Navin_Visvader_NORM_panCK_Reference, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Visvader_0040_ER_Navin_Visvader_NORM_panCK_Reference <- RunPCA(Visvader_0040_ER_Navin_Visvader_NORM_panCK_Reference, features = VariableFeatures(object = Visvader_0040_ER_Navin_Visvader_NORM_panCK_Reference))

# Examine and visualize PCA results a few different ways
print(Visvader_0040_ER_Navin_Visvader_NORM_panCK_Reference[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Visvader_0040_ER_Navin_Visvader_NORM_panCK_Reference, dims = 1:2, reduction = "pca")
DimPlot(Visvader_0040_ER_Navin_Visvader_NORM_panCK_Reference, reduction = "pca")

Visvader_0040_ER_Navin_Visvader_NORM_panCK_Reference <- FindNeighbors(Visvader_0040_ER_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
Visvader_0040_ER_Navin_Visvader_NORM_panCK_Reference <- FindClusters(Visvader_0040_ER_Navin_Visvader_NORM_panCK_Reference, resolution = 0.1)

# Umap clustering
Visvader_0040_ER_Navin_Visvader_NORM_panCK_Reference <- RunUMAP(Visvader_0040_ER_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
DimPlot(Visvader_0040_ER_Navin_Visvader_NORM_panCK_Reference, reduction = "umap")
DimPlot(Visvader_0040_ER_Navin_Visvader_NORM_panCK_Reference, reduction = "umap", group.by = "orig.ident")

# Save RDS
saveRDS(Visvader_0040_ER_Navin_Visvader_NORM_panCK_Reference, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Visvader_0040_ER_Navin_Visvader_NORM_panCK_Reference.rds")

########## Prepare Cell Annotations
# Extract & Save Meta Data
Visvader_0040_ER_Navin_Visvader_NORM_panCK_Reference.meta.data <- as.data.frame(as.matrix(Visvader_0040_ER_Navin_Visvader_NORM_panCK_Reference@meta.data))
# Extract pertinent columns
# Make row names first column
Visvader_0040_ER_Navin_Visvader_NORM_panCK_Reference.meta.data <- tibble::rownames_to_column(Visvader_0040_ER_Navin_Visvader_NORM_panCK_Reference.meta.data, "Cells")
# In Seurat, "X" becomes the cell names while "cnv.ident" is specified above
Visvader_0040_ER_Navin_Visvader_NORM_panCK_Reference.annotations <- Visvader_0040_ER_Navin_Visvader_NORM_panCK_Reference.meta.data[, c("Cells", "cnv.ident")]
# Delete header
names(Visvader_0040_ER_Navin_Visvader_NORM_panCK_Reference.annotations) <- NULL
# Save as table
write.table(Visvader_0040_ER_Navin_Visvader_NORM_panCK_Reference.annotations,"/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_0040_ER_Navin_Visvader_NORM_panCK_Reference.annotations.txt",sep="\t",row.names=FALSE)
# Remove metadata and annotations
rm(Visvader_0040_ER_Navin_Visvader_NORM_panCK_Reference.meta.data)
rm(Visvader_0040_ER_Navin_Visvader_NORM_panCK_Reference.annotations)

########## Generate Count Matrix
# Extract counts from the Seurat object
Visvader_0040_ER_Navin_Visvader_NORM_panCK_Reference.count <- Visvader_0040_ER_Navin_Visvader_NORM_panCK_Reference[["RNA"]]$counts
# Convert counts to a sparse matrix
Visvader_0040_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as(Visvader_0040_ER_Navin_Visvader_NORM_panCK_Reference.count, 'sparseMatrix')
rm(Visvader_0040_ER_Navin_Visvader_NORM_panCK_Reference.count)

# Convert to requisite data format
Visvader_0040_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as.data.frame(Visvader_0040_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix)
# Save table
write.csv(Visvader_0040_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_Count_Matrix/Visvader_0040_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv")
# Remove Objects
rm(Visvader_0040_ER_Navin_Visvader_NORM_panCK_Reference)
rm(Visvader_0040_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix)

########### Run inferCNV
# Note: check.names = FALSE prevents cell names from having dashes replaced with dots; essential to match annotations file
Visvader_0040_ER_inferCNV <- CreateInfercnvObject(raw_counts_matrix=read.csv(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_Count_Matrix/Visvader_0040_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv", header = TRUE, row.names = 1,check.names = FALSE),
                                                  annotations_file=read.table(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_0040_ER_Navin_Visvader_NORM_panCK_Reference.annotations.txt", header = FALSE, row.names = 1),
                                                  delim="\t",
                                                  gene_order_file="/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/hg38_gencode_v27.txt",
                                                  ref_group_names=c("Normal"))


# perform infercnv operations to reveal cnv signal
Visvader_0040_ER_inferCNV = infercnv::run(Visvader_0040_ER_inferCNV,
                                          cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                          out_dir="/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0040_ER_panCK/",  # dir is auto-created for storing outputs
                                          cluster_by_groups=T,   # cluster
                                          denoise=T,
                                          HMM=T,
                                          num_ref_groups = 8, analysis_mode = "subclusters")

###########  Integrate data with Seurat object
# Read RDS and add inferCNV results to metadata
Visvader_0040_ER_Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Visvader_0040_ER_Navin_Visvader_NORM_panCK_Reference.rds") 
Visvader_0040_ER_panCK_cnv <- infercnv::add_to_seurat(infercnv_output_path="/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0040_ER_panCK/",
                                                      seurat_obj=Visvader_0040_ER_Navin_Visvader_NORM_panCK_Reference, # optional
                                                      top_n=10)

# Collect Tumor cells from Seurat Object
Visvader_0040_ER_panCK_cnv <- subset(x = Visvader_0040_ER_panCK_cnv, subset = cnv.ident == "Tumor")
# Quantify chromosomes with alterations per cell
Visvader_0040_ER_panCK_cnv_Pos_Chromosomes <- as.data.frame(rowSums(Visvader_0040_ER_panCK_cnv@meta.data[,startsWith(names(Visvader_0040_ER_panCK_cnv@meta.data),"has_cnv_")] > 0 ))
Visvader_0040_ER_panCK_cnv_Neg_Chromosomes <- as.data.frame(rowSums(Visvader_0040_ER_panCK_cnv@meta.data[,startsWith(names(Visvader_0040_ER_panCK_cnv@meta.data),"has_cnv_")] == 0))
Visvader_0040_ER_panCK_cnv_Chromosome_Quant <- as.data.frame(c(Visvader_0040_ER_panCK_cnv_Pos_Chromosomes, Visvader_0040_ER_panCK_cnv_Neg_Chromosomes))
colnames(Visvader_0040_ER_panCK_cnv_Chromosome_Quant) <- c("Chromosomes CNV Pos", "Chromosomes CNV Neg")
write.csv(Visvader_0040_ER_panCK_cnv_Chromosome_Quant, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0040_ER_panCK/Visvader_0040_ER_panCK_cnv_Chromosome_Quant.csv")

# Save CNV Metadata
Visvader_0040_ER_panCK_cnv.meta.data <- as.data.frame(as.matrix(Visvader_0040_ER_panCK_cnv@meta.data))
write.csv(Visvader_0040_ER_panCK_cnv.meta.data, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0040_ER_panCK/Visvader_0040_ER_panCK_cnv.meta.data.csv")

# Subset Cells with inferred CNVs -------------------------------------------------------------------------------------------------------------
Visvader_0040_ER_panCK_cnv <- subset(Visvader_0040_ER_panCK_cnv, cells=rownames(Visvader_0040_ER_panCK_cnv@meta.data)[rowSums(Visvader_0040_ER_panCK_cnv@meta.data[,startsWith(names(Visvader_0040_ER_panCK_cnv@meta.data),"has_cnv_")]) > 0 ])

# Save RDS with inferredCNVs
saveRDS(Visvader_0040_ER_panCK_cnv, file = "/R/R_Visvader/Visvader_RDS/RDS_inferCNV/Visvader_0040_ER_panCK_cnv.rds")


# Free memory
rm(Visvader_0040_ER_inferCNV)
rm(Visvader_0040_ER_Navin_Visvader_NORM_panCK_Reference)
rm(Visvader_0040_ER_panCK_cnv)
rm(Visvader_0040_ER_panCK_cnv_Pos_Chromosomes)
rm(Visvader_0040_ER_panCK_cnv_Neg_Chromosomes)
rm(Visvader_0040_ER_panCK_cnv_Chromosome_Quant)
rm(Visvader_0040_ER_panCK_cnv.meta.data)
gc()



##########################################
#### InferCNV: Visvader_0042_ER_panCK ####
##########################################

########## Merge with normal mammary glad, extract epithelium, and process via Seurat
# Load RDS
Visvader_0042_ER_panCK <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0042_ER_panCK.rds")
Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Navin_Visvader_NORM_panCK_Reference.rds")

# Set identities 
colnames(x = Visvader_0042_ER_panCK[[]])
Idents(object = Visvader_0042_ER_panCK) <- 'Tumor'
Visvader_0042_ER_panCK[["cnv.ident"]] <- Idents(object = Visvader_0042_ER_panCK)
levels(x = Visvader_0042_ER_panCK)

colnames(x = Navin_Visvader_NORM_panCK_Reference[[]])
Idents(object = Navin_Visvader_NORM_panCK_Reference) <- 'Normal'
Navin_Visvader_NORM_panCK_Reference[["cnv.ident"]] <- Idents(object = Navin_Visvader_NORM_panCK_Reference)
levels(x = Navin_Visvader_NORM_panCK_Reference)

# Merge subsets into recombined RDS
Visvader_0042_ER_Navin_Visvader_NORM_panCK_Reference <- merge(x = Visvader_0042_ER_panCK, y = Navin_Visvader_NORM_panCK_Reference)
# Join Layers
Visvader_0042_ER_Navin_Visvader_NORM_panCK_Reference <- JoinLayers(Visvader_0042_ER_Navin_Visvader_NORM_panCK_Reference)

# Remove subsetted objects
rm(Visvader_0042_ER_panCK)
rm(Navin_Visvader_NORM_panCK_Reference)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Visvader_0042_ER_Navin_Visvader_NORM_panCK_Reference <- FindVariableFeatures(Visvader_0042_ER_Navin_Visvader_NORM_panCK_Reference, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Visvader_0042_ER_Navin_Visvader_NORM_panCK_Reference)
Visvader_0042_ER_Navin_Visvader_NORM_panCK_Reference <- ScaleData(Visvader_0042_ER_Navin_Visvader_NORM_panCK_Reference, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Visvader_0042_ER_Navin_Visvader_NORM_panCK_Reference <- RunPCA(Visvader_0042_ER_Navin_Visvader_NORM_panCK_Reference, features = VariableFeatures(object = Visvader_0042_ER_Navin_Visvader_NORM_panCK_Reference))

# Examine and visualize PCA results a few different ways
print(Visvader_0042_ER_Navin_Visvader_NORM_panCK_Reference[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Visvader_0042_ER_Navin_Visvader_NORM_panCK_Reference, dims = 1:2, reduction = "pca")
DimPlot(Visvader_0042_ER_Navin_Visvader_NORM_panCK_Reference, reduction = "pca")

Visvader_0042_ER_Navin_Visvader_NORM_panCK_Reference <- FindNeighbors(Visvader_0042_ER_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
Visvader_0042_ER_Navin_Visvader_NORM_panCK_Reference <- FindClusters(Visvader_0042_ER_Navin_Visvader_NORM_panCK_Reference, resolution = 0.1)

# Umap clustering
Visvader_0042_ER_Navin_Visvader_NORM_panCK_Reference <- RunUMAP(Visvader_0042_ER_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
DimPlot(Visvader_0042_ER_Navin_Visvader_NORM_panCK_Reference, reduction = "umap")
DimPlot(Visvader_0042_ER_Navin_Visvader_NORM_panCK_Reference, reduction = "umap", group.by = "orig.ident")

# Save RDS
saveRDS(Visvader_0042_ER_Navin_Visvader_NORM_panCK_Reference, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Visvader_0042_ER_Navin_Visvader_NORM_panCK_Reference.rds")

########## Prepare Cell Annotations
# Extract & Save Meta Data
Visvader_0042_ER_Navin_Visvader_NORM_panCK_Reference.meta.data <- as.data.frame(as.matrix(Visvader_0042_ER_Navin_Visvader_NORM_panCK_Reference@meta.data))
# Extract pertinent columns
# Make row names first column
Visvader_0042_ER_Navin_Visvader_NORM_panCK_Reference.meta.data <- tibble::rownames_to_column(Visvader_0042_ER_Navin_Visvader_NORM_panCK_Reference.meta.data, "Cells")
# In Seurat, "X" becomes the cell names while "cnv.ident" is specified above
Visvader_0042_ER_Navin_Visvader_NORM_panCK_Reference.annotations <- Visvader_0042_ER_Navin_Visvader_NORM_panCK_Reference.meta.data[, c("Cells", "cnv.ident")]
# Delete header
names(Visvader_0042_ER_Navin_Visvader_NORM_panCK_Reference.annotations) <- NULL
# Save as table
write.table(Visvader_0042_ER_Navin_Visvader_NORM_panCK_Reference.annotations,"/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_0042_ER_Navin_Visvader_NORM_panCK_Reference.annotations.txt",sep="\t",row.names=FALSE)
# Remove metadata and annotations
rm(Visvader_0042_ER_Navin_Visvader_NORM_panCK_Reference.meta.data)
rm(Visvader_0042_ER_Navin_Visvader_NORM_panCK_Reference.annotations)

########## Generate Count Matrix
# Extract counts from the Seurat object
Visvader_0042_ER_Navin_Visvader_NORM_panCK_Reference.count <- Visvader_0042_ER_Navin_Visvader_NORM_panCK_Reference[["RNA"]]$counts
# Convert counts to a sparse matrix
Visvader_0042_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as(Visvader_0042_ER_Navin_Visvader_NORM_panCK_Reference.count, 'sparseMatrix')
rm(Visvader_0042_ER_Navin_Visvader_NORM_panCK_Reference.count)

# Convert to requisite data format
Visvader_0042_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as.data.frame(Visvader_0042_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix)
# Save table
write.csv(Visvader_0042_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_Count_Matrix/Visvader_0042_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv")
# Remove Objects
rm(Visvader_0042_ER_Navin_Visvader_NORM_panCK_Reference)
rm(Visvader_0042_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix)

########### Run inferCNV
# Note: check.names = FALSE prevents cell names from having dashes replaced with dots; essential to match annotations file
Visvader_0042_ER_inferCNV <- CreateInfercnvObject(raw_counts_matrix=read.csv(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_Count_Matrix/Visvader_0042_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv", header = TRUE, row.names = 1,check.names = FALSE),
                                                  annotations_file=read.table(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_0042_ER_Navin_Visvader_NORM_panCK_Reference.annotations.txt", header = FALSE, row.names = 1),
                                                  delim="\t",
                                                  gene_order_file="/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/hg38_gencode_v27.txt",
                                                  ref_group_names=c("Normal"))


# perform infercnv operations to reveal cnv signal
Visvader_0042_ER_inferCNV = infercnv::run(Visvader_0042_ER_inferCNV,
                                          cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                          out_dir="/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0042_ER_panCK/",  # dir is auto-created for storing outputs
                                          cluster_by_groups=T,   # cluster
                                          denoise=T,
                                          HMM=T,
                                          num_ref_groups = 8, analysis_mode = "subclusters")

###########  Integrate data with Seurat object
# Read RDS and add inferCNV results to metadata
Visvader_0042_ER_Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Visvader_0042_ER_Navin_Visvader_NORM_panCK_Reference.rds") 
Visvader_0042_ER_panCK_cnv <- infercnv::add_to_seurat(infercnv_output_path="/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0042_ER_panCK/",
                                                      seurat_obj=Visvader_0042_ER_Navin_Visvader_NORM_panCK_Reference, # optional
                                                      top_n=10)

# Collect Tumor cells from Seurat Object
Visvader_0042_ER_panCK_cnv <- subset(x = Visvader_0042_ER_panCK_cnv, subset = cnv.ident == "Tumor")
# Quantify chromosomes with alterations per cell
Visvader_0042_ER_panCK_cnv_Pos_Chromosomes <- as.data.frame(rowSums(Visvader_0042_ER_panCK_cnv@meta.data[,startsWith(names(Visvader_0042_ER_panCK_cnv@meta.data),"has_cnv_")] > 0 ))
Visvader_0042_ER_panCK_cnv_Neg_Chromosomes <- as.data.frame(rowSums(Visvader_0042_ER_panCK_cnv@meta.data[,startsWith(names(Visvader_0042_ER_panCK_cnv@meta.data),"has_cnv_")] == 0))
Visvader_0042_ER_panCK_cnv_Chromosome_Quant <- as.data.frame(c(Visvader_0042_ER_panCK_cnv_Pos_Chromosomes, Visvader_0042_ER_panCK_cnv_Neg_Chromosomes))
colnames(Visvader_0042_ER_panCK_cnv_Chromosome_Quant) <- c("Chromosomes CNV Pos", "Chromosomes CNV Neg")
write.csv(Visvader_0042_ER_panCK_cnv_Chromosome_Quant, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0042_ER_panCK/Visvader_0042_ER_panCK_cnv_Chromosome_Quant.csv")

# Save CNV Metadata
Visvader_0042_ER_panCK_cnv.meta.data <- as.data.frame(as.matrix(Visvader_0042_ER_panCK_cnv@meta.data))
write.csv(Visvader_0042_ER_panCK_cnv.meta.data, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0042_ER_panCK/Visvader_0042_ER_panCK_cnv.meta.data.csv")

# Subset Cells with inferred CNVs -------------------------------------------------------------------------------------------------------------
Visvader_0042_ER_panCK_cnv <- subset(Visvader_0042_ER_panCK_cnv, cells=rownames(Visvader_0042_ER_panCK_cnv@meta.data)[rowSums(Visvader_0042_ER_panCK_cnv@meta.data[,startsWith(names(Visvader_0042_ER_panCK_cnv@meta.data),"has_cnv_")]) > 0 ])

# Save RDS with inferredCNVs
saveRDS(Visvader_0042_ER_panCK_cnv, file = "/R/R_Visvader/Visvader_RDS/RDS_inferCNV/Visvader_0042_ER_panCK_cnv.rds")


# Free memory
rm(Visvader_0042_ER_inferCNV)
rm(Visvader_0042_ER_Navin_Visvader_NORM_panCK_Reference)
rm(Visvader_0042_ER_panCK_cnv)
rm(Visvader_0042_ER_panCK_cnv_Pos_Chromosomes)
rm(Visvader_0042_ER_panCK_cnv_Neg_Chromosomes)
rm(Visvader_0042_ER_panCK_cnv_Chromosome_Quant)
rm(Visvader_0042_ER_panCK_cnv.meta.data)
gc()



##########################################
#### InferCNV: Visvader_0043_ER_panCK ####
##########################################

########## Merge with normal mammary glad, extract epithelium, and process via Seurat
# Load RDS
Visvader_0043_ER_panCK <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0043_ER_panCK.rds")
Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Navin_Visvader_NORM_panCK_Reference.rds")

# Set identities 
colnames(x = Visvader_0043_ER_panCK[[]])
Idents(object = Visvader_0043_ER_panCK) <- 'Tumor'
Visvader_0043_ER_panCK[["cnv.ident"]] <- Idents(object = Visvader_0043_ER_panCK)
levels(x = Visvader_0043_ER_panCK)

colnames(x = Navin_Visvader_NORM_panCK_Reference[[]])
Idents(object = Navin_Visvader_NORM_panCK_Reference) <- 'Normal'
Navin_Visvader_NORM_panCK_Reference[["cnv.ident"]] <- Idents(object = Navin_Visvader_NORM_panCK_Reference)
levels(x = Navin_Visvader_NORM_panCK_Reference)

# Merge subsets into recombined RDS
Visvader_0043_ER_Navin_Visvader_NORM_panCK_Reference <- merge(x = Visvader_0043_ER_panCK, y = Navin_Visvader_NORM_panCK_Reference)
# Join Layers
Visvader_0043_ER_Navin_Visvader_NORM_panCK_Reference <- JoinLayers(Visvader_0043_ER_Navin_Visvader_NORM_panCK_Reference)

# Remove subsetted objects
rm(Visvader_0043_ER_panCK)
rm(Navin_Visvader_NORM_panCK_Reference)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Visvader_0043_ER_Navin_Visvader_NORM_panCK_Reference <- FindVariableFeatures(Visvader_0043_ER_Navin_Visvader_NORM_panCK_Reference, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Visvader_0043_ER_Navin_Visvader_NORM_panCK_Reference)
Visvader_0043_ER_Navin_Visvader_NORM_panCK_Reference <- ScaleData(Visvader_0043_ER_Navin_Visvader_NORM_panCK_Reference, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Visvader_0043_ER_Navin_Visvader_NORM_panCK_Reference <- RunPCA(Visvader_0043_ER_Navin_Visvader_NORM_panCK_Reference, features = VariableFeatures(object = Visvader_0043_ER_Navin_Visvader_NORM_panCK_Reference))

# Examine and visualize PCA results a few different ways
print(Visvader_0043_ER_Navin_Visvader_NORM_panCK_Reference[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Visvader_0043_ER_Navin_Visvader_NORM_panCK_Reference, dims = 1:2, reduction = "pca")
DimPlot(Visvader_0043_ER_Navin_Visvader_NORM_panCK_Reference, reduction = "pca")

Visvader_0043_ER_Navin_Visvader_NORM_panCK_Reference <- FindNeighbors(Visvader_0043_ER_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
Visvader_0043_ER_Navin_Visvader_NORM_panCK_Reference <- FindClusters(Visvader_0043_ER_Navin_Visvader_NORM_panCK_Reference, resolution = 0.1)

# Umap clustering
Visvader_0043_ER_Navin_Visvader_NORM_panCK_Reference <- RunUMAP(Visvader_0043_ER_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
DimPlot(Visvader_0043_ER_Navin_Visvader_NORM_panCK_Reference, reduction = "umap")
DimPlot(Visvader_0043_ER_Navin_Visvader_NORM_panCK_Reference, reduction = "umap", group.by = "orig.ident")

# Save RDS
saveRDS(Visvader_0043_ER_Navin_Visvader_NORM_panCK_Reference, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Visvader_0043_ER_Navin_Visvader_NORM_panCK_Reference.rds")

########## Prepare Cell Annotations
# Extract & Save Meta Data
Visvader_0043_ER_Navin_Visvader_NORM_panCK_Reference.meta.data <- as.data.frame(as.matrix(Visvader_0043_ER_Navin_Visvader_NORM_panCK_Reference@meta.data))
# Extract pertinent columns
# Make row names first column
Visvader_0043_ER_Navin_Visvader_NORM_panCK_Reference.meta.data <- tibble::rownames_to_column(Visvader_0043_ER_Navin_Visvader_NORM_panCK_Reference.meta.data, "Cells")
# In Seurat, "X" becomes the cell names while "cnv.ident" is specified above
Visvader_0043_ER_Navin_Visvader_NORM_panCK_Reference.annotations <- Visvader_0043_ER_Navin_Visvader_NORM_panCK_Reference.meta.data[, c("Cells", "cnv.ident")]
# Delete header
names(Visvader_0043_ER_Navin_Visvader_NORM_panCK_Reference.annotations) <- NULL
# Save as table
write.table(Visvader_0043_ER_Navin_Visvader_NORM_panCK_Reference.annotations,"/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_0043_ER_Navin_Visvader_NORM_panCK_Reference.annotations.txt",sep="\t",row.names=FALSE)
# Remove metadata and annotations
rm(Visvader_0043_ER_Navin_Visvader_NORM_panCK_Reference.meta.data)
rm(Visvader_0043_ER_Navin_Visvader_NORM_panCK_Reference.annotations)

########## Generate Count Matrix
# Extract counts from the Seurat object
Visvader_0043_ER_Navin_Visvader_NORM_panCK_Reference.count <- Visvader_0043_ER_Navin_Visvader_NORM_panCK_Reference[["RNA"]]$counts
# Convert counts to a sparse matrix
Visvader_0043_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as(Visvader_0043_ER_Navin_Visvader_NORM_panCK_Reference.count, 'sparseMatrix')
rm(Visvader_0043_ER_Navin_Visvader_NORM_panCK_Reference.count)

# Convert to requisite data format
Visvader_0043_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as.data.frame(Visvader_0043_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix)
# Save table
write.csv(Visvader_0043_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_Count_Matrix/Visvader_0043_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv")
# Remove Objects
rm(Visvader_0043_ER_Navin_Visvader_NORM_panCK_Reference)
rm(Visvader_0043_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix)

########### Run inferCNV
# Note: check.names = FALSE prevents cell names from having dashes replaced with dots; essential to match annotations file
Visvader_0043_ER_inferCNV <- CreateInfercnvObject(raw_counts_matrix=read.csv(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_Count_Matrix/Visvader_0043_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv", header = TRUE, row.names = 1,check.names = FALSE),
                                                  annotations_file=read.table(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_0043_ER_Navin_Visvader_NORM_panCK_Reference.annotations.txt", header = FALSE, row.names = 1),
                                                  delim="\t",
                                                  gene_order_file="/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/hg38_gencode_v27.txt",
                                                  ref_group_names=c("Normal"))


# perform infercnv operations to reveal cnv signal
Visvader_0043_ER_inferCNV = infercnv::run(Visvader_0043_ER_inferCNV,
                                          cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                          out_dir="/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0043_ER_panCK/",  # dir is auto-created for storing outputs
                                          cluster_by_groups=T,   # cluster
                                          denoise=T,
                                          HMM=T,
                                          num_ref_groups = 8, analysis_mode = "subclusters")

###########  Integrate data with Seurat object
# Read RDS and add inferCNV results to metadata
Visvader_0043_ER_Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Visvader_0043_ER_Navin_Visvader_NORM_panCK_Reference.rds") 
Visvader_0043_ER_panCK_cnv <- infercnv::add_to_seurat(infercnv_output_path="/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0043_ER_panCK/",
                                                      seurat_obj=Visvader_0043_ER_Navin_Visvader_NORM_panCK_Reference, # optional
                                                      top_n=10)

# Collect Tumor cells from Seurat Object
Visvader_0043_ER_panCK_cnv <- subset(x = Visvader_0043_ER_panCK_cnv, subset = cnv.ident == "Tumor")
# Quantify chromosomes with alterations per cell
Visvader_0043_ER_panCK_cnv_Pos_Chromosomes <- as.data.frame(rowSums(Visvader_0043_ER_panCK_cnv@meta.data[,startsWith(names(Visvader_0043_ER_panCK_cnv@meta.data),"has_cnv_")] > 0 ))
Visvader_0043_ER_panCK_cnv_Neg_Chromosomes <- as.data.frame(rowSums(Visvader_0043_ER_panCK_cnv@meta.data[,startsWith(names(Visvader_0043_ER_panCK_cnv@meta.data),"has_cnv_")] == 0))
Visvader_0043_ER_panCK_cnv_Chromosome_Quant <- as.data.frame(c(Visvader_0043_ER_panCK_cnv_Pos_Chromosomes, Visvader_0043_ER_panCK_cnv_Neg_Chromosomes))
colnames(Visvader_0043_ER_panCK_cnv_Chromosome_Quant) <- c("Chromosomes CNV Pos", "Chromosomes CNV Neg")
write.csv(Visvader_0043_ER_panCK_cnv_Chromosome_Quant, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0043_ER_panCK/Visvader_0043_ER_panCK_cnv_Chromosome_Quant.csv")

# Save CNV Metadata
Visvader_0043_ER_panCK_cnv.meta.data <- as.data.frame(as.matrix(Visvader_0043_ER_panCK_cnv@meta.data))
write.csv(Visvader_0043_ER_panCK_cnv.meta.data, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0043_ER_panCK/Visvader_0043_ER_panCK_cnv.meta.data.csv")

# Subset Cells with inferred CNVs -------------------------------------------------------------------------------------------------------------
Visvader_0043_ER_panCK_cnv <- subset(Visvader_0043_ER_panCK_cnv, cells=rownames(Visvader_0043_ER_panCK_cnv@meta.data)[rowSums(Visvader_0043_ER_panCK_cnv@meta.data[,startsWith(names(Visvader_0043_ER_panCK_cnv@meta.data),"has_cnv_")]) > 0 ])

# Save RDS with inferredCNVs
saveRDS(Visvader_0043_ER_panCK_cnv, file = "/R/R_Visvader/Visvader_RDS/RDS_inferCNV/Visvader_0043_ER_panCK_cnv.rds")


# Free memory
rm(Visvader_0043_ER_inferCNV)
rm(Visvader_0043_ER_Navin_Visvader_NORM_panCK_Reference)
rm(Visvader_0043_ER_panCK_cnv)
rm(Visvader_0043_ER_panCK_cnv_Pos_Chromosomes)
rm(Visvader_0043_ER_panCK_cnv_Neg_Chromosomes)
rm(Visvader_0043_ER_panCK_cnv_Chromosome_Quant)
rm(Visvader_0043_ER_panCK_cnv.meta.data)
gc()



##########################################
#### InferCNV: Visvader_0056_ER_panCK ####
##########################################

########## Merge with normal mammary glad, extract epithelium, and process via Seurat
# Load RDS
Visvader_0056_ER_panCK <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0056_ER_panCK.rds")
Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Navin_Visvader_NORM_panCK_Reference.rds")

# Set identities 
colnames(x = Visvader_0056_ER_panCK[[]])
Idents(object = Visvader_0056_ER_panCK) <- 'Tumor'
Visvader_0056_ER_panCK[["cnv.ident"]] <- Idents(object = Visvader_0056_ER_panCK)
levels(x = Visvader_0056_ER_panCK)

colnames(x = Navin_Visvader_NORM_panCK_Reference[[]])
Idents(object = Navin_Visvader_NORM_panCK_Reference) <- 'Normal'
Navin_Visvader_NORM_panCK_Reference[["cnv.ident"]] <- Idents(object = Navin_Visvader_NORM_panCK_Reference)
levels(x = Navin_Visvader_NORM_panCK_Reference)

# Merge subsets into recombined RDS
Visvader_0056_ER_Navin_Visvader_NORM_panCK_Reference <- merge(x = Visvader_0056_ER_panCK, y = Navin_Visvader_NORM_panCK_Reference)
# Join Layers
Visvader_0056_ER_Navin_Visvader_NORM_panCK_Reference <- JoinLayers(Visvader_0056_ER_Navin_Visvader_NORM_panCK_Reference)

# Remove subsetted objects
rm(Visvader_0056_ER_panCK)
rm(Navin_Visvader_NORM_panCK_Reference)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Visvader_0056_ER_Navin_Visvader_NORM_panCK_Reference <- FindVariableFeatures(Visvader_0056_ER_Navin_Visvader_NORM_panCK_Reference, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Visvader_0056_ER_Navin_Visvader_NORM_panCK_Reference)
Visvader_0056_ER_Navin_Visvader_NORM_panCK_Reference <- ScaleData(Visvader_0056_ER_Navin_Visvader_NORM_panCK_Reference, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Visvader_0056_ER_Navin_Visvader_NORM_panCK_Reference <- RunPCA(Visvader_0056_ER_Navin_Visvader_NORM_panCK_Reference, features = VariableFeatures(object = Visvader_0056_ER_Navin_Visvader_NORM_panCK_Reference))

# Examine and visualize PCA results a few different ways
print(Visvader_0056_ER_Navin_Visvader_NORM_panCK_Reference[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Visvader_0056_ER_Navin_Visvader_NORM_panCK_Reference, dims = 1:2, reduction = "pca")
DimPlot(Visvader_0056_ER_Navin_Visvader_NORM_panCK_Reference, reduction = "pca")

Visvader_0056_ER_Navin_Visvader_NORM_panCK_Reference <- FindNeighbors(Visvader_0056_ER_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
Visvader_0056_ER_Navin_Visvader_NORM_panCK_Reference <- FindClusters(Visvader_0056_ER_Navin_Visvader_NORM_panCK_Reference, resolution = 0.1)

# Umap clustering
Visvader_0056_ER_Navin_Visvader_NORM_panCK_Reference <- RunUMAP(Visvader_0056_ER_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
DimPlot(Visvader_0056_ER_Navin_Visvader_NORM_panCK_Reference, reduction = "umap")
DimPlot(Visvader_0056_ER_Navin_Visvader_NORM_panCK_Reference, reduction = "umap", group.by = "orig.ident")

# Save RDS
saveRDS(Visvader_0056_ER_Navin_Visvader_NORM_panCK_Reference, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Visvader_0056_ER_Navin_Visvader_NORM_panCK_Reference.rds")

########## Prepare Cell Annotations
# Extract & Save Meta Data
Visvader_0056_ER_Navin_Visvader_NORM_panCK_Reference.meta.data <- as.data.frame(as.matrix(Visvader_0056_ER_Navin_Visvader_NORM_panCK_Reference@meta.data))
# Extract pertinent columns
# Make row names first column
Visvader_0056_ER_Navin_Visvader_NORM_panCK_Reference.meta.data <- tibble::rownames_to_column(Visvader_0056_ER_Navin_Visvader_NORM_panCK_Reference.meta.data, "Cells")
# In Seurat, "X" becomes the cell names while "cnv.ident" is specified above
Visvader_0056_ER_Navin_Visvader_NORM_panCK_Reference.annotations <- Visvader_0056_ER_Navin_Visvader_NORM_panCK_Reference.meta.data[, c("Cells", "cnv.ident")]
# Delete header
names(Visvader_0056_ER_Navin_Visvader_NORM_panCK_Reference.annotations) <- NULL
# Save as table
write.table(Visvader_0056_ER_Navin_Visvader_NORM_panCK_Reference.annotations,"/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_0056_ER_Navin_Visvader_NORM_panCK_Reference.annotations.txt",sep="\t",row.names=FALSE)
# Remove metadata and annotations
rm(Visvader_0056_ER_Navin_Visvader_NORM_panCK_Reference.meta.data)
rm(Visvader_0056_ER_Navin_Visvader_NORM_panCK_Reference.annotations)

########## Generate Count Matrix
# Extract counts from the Seurat object
Visvader_0056_ER_Navin_Visvader_NORM_panCK_Reference.count <- Visvader_0056_ER_Navin_Visvader_NORM_panCK_Reference[["RNA"]]$counts
# Convert counts to a sparse matrix
Visvader_0056_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as(Visvader_0056_ER_Navin_Visvader_NORM_panCK_Reference.count, 'sparseMatrix')
rm(Visvader_0056_ER_Navin_Visvader_NORM_panCK_Reference.count)

# Convert to requisite data format
Visvader_0056_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as.data.frame(Visvader_0056_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix)
# Save table
write.csv(Visvader_0056_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_Count_Matrix/Visvader_0056_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv")
# Remove Objects
rm(Visvader_0056_ER_Navin_Visvader_NORM_panCK_Reference)
rm(Visvader_0056_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix)

########### Run inferCNV
# Note: check.names = FALSE prevents cell names from having dashes replaced with dots; essential to match annotations file
Visvader_0056_ER_inferCNV <- CreateInfercnvObject(raw_counts_matrix=read.csv(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_Count_Matrix/Visvader_0056_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv", header = TRUE, row.names = 1,check.names = FALSE),
                                                  annotations_file=read.table(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_0056_ER_Navin_Visvader_NORM_panCK_Reference.annotations.txt", header = FALSE, row.names = 1),
                                                  delim="\t",
                                                  gene_order_file="/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/hg38_gencode_v27.txt",
                                                  ref_group_names=c("Normal"))


# perform infercnv operations to reveal cnv signal
Visvader_0056_ER_inferCNV = infercnv::run(Visvader_0056_ER_inferCNV,
                                          cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                          out_dir="/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0056_ER_panCK/",  # dir is auto-created for storing outputs
                                          cluster_by_groups=T,   # cluster
                                          denoise=T,
                                          HMM=T,
                                          num_ref_groups = 8, analysis_mode = "subclusters")

###########  Integrate data with Seurat object
# Read RDS and add inferCNV results to metadata
Visvader_0056_ER_Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Visvader_0056_ER_Navin_Visvader_NORM_panCK_Reference.rds") 
Visvader_0056_ER_panCK_cnv <- infercnv::add_to_seurat(infercnv_output_path="/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0056_ER_panCK/",
                                                      seurat_obj=Visvader_0056_ER_Navin_Visvader_NORM_panCK_Reference, # optional
                                                      top_n=10)

# Collect Tumor cells from Seurat Object
Visvader_0056_ER_panCK_cnv <- subset(x = Visvader_0056_ER_panCK_cnv, subset = cnv.ident == "Tumor")
# Quantify chromosomes with alterations per cell
Visvader_0056_ER_panCK_cnv_Pos_Chromosomes <- as.data.frame(rowSums(Visvader_0056_ER_panCK_cnv@meta.data[,startsWith(names(Visvader_0056_ER_panCK_cnv@meta.data),"has_cnv_")] > 0 ))
Visvader_0056_ER_panCK_cnv_Neg_Chromosomes <- as.data.frame(rowSums(Visvader_0056_ER_panCK_cnv@meta.data[,startsWith(names(Visvader_0056_ER_panCK_cnv@meta.data),"has_cnv_")] == 0))
Visvader_0056_ER_panCK_cnv_Chromosome_Quant <- as.data.frame(c(Visvader_0056_ER_panCK_cnv_Pos_Chromosomes, Visvader_0056_ER_panCK_cnv_Neg_Chromosomes))
colnames(Visvader_0056_ER_panCK_cnv_Chromosome_Quant) <- c("Chromosomes CNV Pos", "Chromosomes CNV Neg")
write.csv(Visvader_0056_ER_panCK_cnv_Chromosome_Quant, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0056_ER_panCK/Visvader_0056_ER_panCK_cnv_Chromosome_Quant.csv")

# Save CNV Metadata
Visvader_0056_ER_panCK_cnv.meta.data <- as.data.frame(as.matrix(Visvader_0056_ER_panCK_cnv@meta.data))
write.csv(Visvader_0056_ER_panCK_cnv.meta.data, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0056_ER_panCK/Visvader_0056_ER_panCK_cnv.meta.data.csv")

# Subset Cells with inferred CNVs -------------------------------------------------------------------------------------------------------------
Visvader_0056_ER_panCK_cnv <- subset(Visvader_0056_ER_panCK_cnv, cells=rownames(Visvader_0056_ER_panCK_cnv@meta.data)[rowSums(Visvader_0056_ER_panCK_cnv@meta.data[,startsWith(names(Visvader_0056_ER_panCK_cnv@meta.data),"has_cnv_")]) > 0 ])

# Save RDS with inferredCNVs
saveRDS(Visvader_0056_ER_panCK_cnv, file = "/R/R_Visvader/Visvader_RDS/RDS_inferCNV/Visvader_0056_ER_panCK_cnv.rds")


# Free memory
rm(Visvader_0056_ER_inferCNV)
rm(Visvader_0056_ER_Navin_Visvader_NORM_panCK_Reference)
rm(Visvader_0056_ER_panCK_cnv)
rm(Visvader_0056_ER_panCK_cnv_Pos_Chromosomes)
rm(Visvader_0056_ER_panCK_cnv_Neg_Chromosomes)
rm(Visvader_0056_ER_panCK_cnv_Chromosome_Quant)
rm(Visvader_0056_ER_panCK_cnv.meta.data)
gc()



##########################################
#### InferCNV: Visvader_0064_ER_panCK ####
##########################################

########## Merge with normal mammary glad, extract epithelium, and process via Seurat
# Load RDS
Visvader_0064_ER_panCK <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0064_ER_panCK.rds")
Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Navin_Visvader_NORM_panCK_Reference.rds")

# Set identities 
colnames(x = Visvader_0064_ER_panCK[[]])
Idents(object = Visvader_0064_ER_panCK) <- 'Tumor'
Visvader_0064_ER_panCK[["cnv.ident"]] <- Idents(object = Visvader_0064_ER_panCK)
levels(x = Visvader_0064_ER_panCK)

colnames(x = Navin_Visvader_NORM_panCK_Reference[[]])
Idents(object = Navin_Visvader_NORM_panCK_Reference) <- 'Normal'
Navin_Visvader_NORM_panCK_Reference[["cnv.ident"]] <- Idents(object = Navin_Visvader_NORM_panCK_Reference)
levels(x = Navin_Visvader_NORM_panCK_Reference)

# Merge subsets into recombined RDS
Visvader_0064_ER_Navin_Visvader_NORM_panCK_Reference <- merge(x = Visvader_0064_ER_panCK, y = Navin_Visvader_NORM_panCK_Reference)
# Join Layers
Visvader_0064_ER_Navin_Visvader_NORM_panCK_Reference <- JoinLayers(Visvader_0064_ER_Navin_Visvader_NORM_panCK_Reference)

# Remove subsetted objects
rm(Visvader_0064_ER_panCK)
rm(Navin_Visvader_NORM_panCK_Reference)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Visvader_0064_ER_Navin_Visvader_NORM_panCK_Reference <- FindVariableFeatures(Visvader_0064_ER_Navin_Visvader_NORM_panCK_Reference, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Visvader_0064_ER_Navin_Visvader_NORM_panCK_Reference)
Visvader_0064_ER_Navin_Visvader_NORM_panCK_Reference <- ScaleData(Visvader_0064_ER_Navin_Visvader_NORM_panCK_Reference, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Visvader_0064_ER_Navin_Visvader_NORM_panCK_Reference <- RunPCA(Visvader_0064_ER_Navin_Visvader_NORM_panCK_Reference, features = VariableFeatures(object = Visvader_0064_ER_Navin_Visvader_NORM_panCK_Reference))

# Examine and visualize PCA results a few different ways
print(Visvader_0064_ER_Navin_Visvader_NORM_panCK_Reference[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Visvader_0064_ER_Navin_Visvader_NORM_panCK_Reference, dims = 1:2, reduction = "pca")
DimPlot(Visvader_0064_ER_Navin_Visvader_NORM_panCK_Reference, reduction = "pca")

Visvader_0064_ER_Navin_Visvader_NORM_panCK_Reference <- FindNeighbors(Visvader_0064_ER_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
Visvader_0064_ER_Navin_Visvader_NORM_panCK_Reference <- FindClusters(Visvader_0064_ER_Navin_Visvader_NORM_panCK_Reference, resolution = 0.1)

# Umap clustering
Visvader_0064_ER_Navin_Visvader_NORM_panCK_Reference <- RunUMAP(Visvader_0064_ER_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
DimPlot(Visvader_0064_ER_Navin_Visvader_NORM_panCK_Reference, reduction = "umap")
DimPlot(Visvader_0064_ER_Navin_Visvader_NORM_panCK_Reference, reduction = "umap", group.by = "orig.ident")

# Save RDS
saveRDS(Visvader_0064_ER_Navin_Visvader_NORM_panCK_Reference, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Visvader_0064_ER_Navin_Visvader_NORM_panCK_Reference.rds")

########## Prepare Cell Annotations
# Extract & Save Meta Data
Visvader_0064_ER_Navin_Visvader_NORM_panCK_Reference.meta.data <- as.data.frame(as.matrix(Visvader_0064_ER_Navin_Visvader_NORM_panCK_Reference@meta.data))
# Extract pertinent columns
# Make row names first column
Visvader_0064_ER_Navin_Visvader_NORM_panCK_Reference.meta.data <- tibble::rownames_to_column(Visvader_0064_ER_Navin_Visvader_NORM_panCK_Reference.meta.data, "Cells")
# In Seurat, "X" becomes the cell names while "cnv.ident" is specified above
Visvader_0064_ER_Navin_Visvader_NORM_panCK_Reference.annotations <- Visvader_0064_ER_Navin_Visvader_NORM_panCK_Reference.meta.data[, c("Cells", "cnv.ident")]
# Delete header
names(Visvader_0064_ER_Navin_Visvader_NORM_panCK_Reference.annotations) <- NULL
# Save as table
write.table(Visvader_0064_ER_Navin_Visvader_NORM_panCK_Reference.annotations,"/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_0064_ER_Navin_Visvader_NORM_panCK_Reference.annotations.txt",sep="\t",row.names=FALSE)
# Remove metadata and annotations
rm(Visvader_0064_ER_Navin_Visvader_NORM_panCK_Reference.meta.data)
rm(Visvader_0064_ER_Navin_Visvader_NORM_panCK_Reference.annotations)

########## Generate Count Matrix
# Extract counts from the Seurat object
Visvader_0064_ER_Navin_Visvader_NORM_panCK_Reference.count <- Visvader_0064_ER_Navin_Visvader_NORM_panCK_Reference[["RNA"]]$counts
# Convert counts to a sparse matrix
Visvader_0064_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as(Visvader_0064_ER_Navin_Visvader_NORM_panCK_Reference.count, 'sparseMatrix')
rm(Visvader_0064_ER_Navin_Visvader_NORM_panCK_Reference.count)

# Convert to requisite data format
Visvader_0064_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as.data.frame(Visvader_0064_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix)
# Save table
write.csv(Visvader_0064_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_Count_Matrix/Visvader_0064_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv")
# Remove Objects
rm(Visvader_0064_ER_Navin_Visvader_NORM_panCK_Reference)
rm(Visvader_0064_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix)

########### Run inferCNV
# Note: check.names = FALSE prevents cell names from having dashes replaced with dots; essential to match annotations file
Visvader_0064_ER_inferCNV <- CreateInfercnvObject(raw_counts_matrix=read.csv(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_Count_Matrix/Visvader_0064_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv", header = TRUE, row.names = 1,check.names = FALSE),
                                                  annotations_file=read.table(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_0064_ER_Navin_Visvader_NORM_panCK_Reference.annotations.txt", header = FALSE, row.names = 1),
                                                  delim="\t",
                                                  gene_order_file="/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/hg38_gencode_v27.txt",
                                                  ref_group_names=c("Normal"))


# perform infercnv operations to reveal cnv signal
Visvader_0064_ER_inferCNV = infercnv::run(Visvader_0064_ER_inferCNV,
                                          cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                          out_dir="/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0064_ER_panCK/",  # dir is auto-created for storing outputs
                                          cluster_by_groups=T,   # cluster
                                          denoise=T,
                                          HMM=T,
                                          num_ref_groups = 8, analysis_mode = "subclusters")

###########  Integrate data with Seurat object
# Read RDS and add inferCNV results to metadata
Visvader_0064_ER_Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Visvader_0064_ER_Navin_Visvader_NORM_panCK_Reference.rds") 
Visvader_0064_ER_panCK_cnv <- infercnv::add_to_seurat(infercnv_output_path="/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0064_ER_panCK/",
                                                      seurat_obj=Visvader_0064_ER_Navin_Visvader_NORM_panCK_Reference, # optional
                                                      top_n=10)

# Collect Tumor cells from Seurat Object
Visvader_0064_ER_panCK_cnv <- subset(x = Visvader_0064_ER_panCK_cnv, subset = cnv.ident == "Tumor")
# Quantify chromosomes with alterations per cell
Visvader_0064_ER_panCK_cnv_Pos_Chromosomes <- as.data.frame(rowSums(Visvader_0064_ER_panCK_cnv@meta.data[,startsWith(names(Visvader_0064_ER_panCK_cnv@meta.data),"has_cnv_")] > 0 ))
Visvader_0064_ER_panCK_cnv_Neg_Chromosomes <- as.data.frame(rowSums(Visvader_0064_ER_panCK_cnv@meta.data[,startsWith(names(Visvader_0064_ER_panCK_cnv@meta.data),"has_cnv_")] == 0))
Visvader_0064_ER_panCK_cnv_Chromosome_Quant <- as.data.frame(c(Visvader_0064_ER_panCK_cnv_Pos_Chromosomes, Visvader_0064_ER_panCK_cnv_Neg_Chromosomes))
colnames(Visvader_0064_ER_panCK_cnv_Chromosome_Quant) <- c("Chromosomes CNV Pos", "Chromosomes CNV Neg")
write.csv(Visvader_0064_ER_panCK_cnv_Chromosome_Quant, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0064_ER_panCK/Visvader_0064_ER_panCK_cnv_Chromosome_Quant.csv")

# Save CNV Metadata
Visvader_0064_ER_panCK_cnv.meta.data <- as.data.frame(as.matrix(Visvader_0064_ER_panCK_cnv@meta.data))
write.csv(Visvader_0064_ER_panCK_cnv.meta.data, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0064_ER_panCK/Visvader_0064_ER_panCK_cnv.meta.data.csv")

# Subset Cells with inferred CNVs -------------------------------------------------------------------------------------------------------------
Visvader_0064_ER_panCK_cnv <- subset(Visvader_0064_ER_panCK_cnv, cells=rownames(Visvader_0064_ER_panCK_cnv@meta.data)[rowSums(Visvader_0064_ER_panCK_cnv@meta.data[,startsWith(names(Visvader_0064_ER_panCK_cnv@meta.data),"has_cnv_")]) > 0 ])

# Save RDS with inferredCNVs
saveRDS(Visvader_0064_ER_panCK_cnv, file = "/R/R_Visvader/Visvader_RDS/RDS_inferCNV/Visvader_0064_ER_panCK_cnv.rds")


# Free memory
rm(Visvader_0064_ER_inferCNV)
rm(Visvader_0064_ER_Navin_Visvader_NORM_panCK_Reference)
rm(Visvader_0064_ER_panCK_cnv)
rm(Visvader_0064_ER_panCK_cnv_Pos_Chromosomes)
rm(Visvader_0064_ER_panCK_cnv_Neg_Chromosomes)
rm(Visvader_0064_ER_panCK_cnv_Chromosome_Quant)
rm(Visvader_0064_ER_panCK_cnv.meta.data)
gc()



##########################################
#### InferCNV: Visvader_0068_ER_panCK ####
##########################################

########## Merge with normal mammary glad, extract epithelium, and process via Seurat
# Load RDS
Visvader_0068_ER_panCK <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0068_ER_panCK.rds")
Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Navin_Visvader_NORM_panCK_Reference.rds")

# Set identities 
colnames(x = Visvader_0068_ER_panCK[[]])
Idents(object = Visvader_0068_ER_panCK) <- 'Tumor'
Visvader_0068_ER_panCK[["cnv.ident"]] <- Idents(object = Visvader_0068_ER_panCK)
levels(x = Visvader_0068_ER_panCK)

colnames(x = Navin_Visvader_NORM_panCK_Reference[[]])
Idents(object = Navin_Visvader_NORM_panCK_Reference) <- 'Normal'
Navin_Visvader_NORM_panCK_Reference[["cnv.ident"]] <- Idents(object = Navin_Visvader_NORM_panCK_Reference)
levels(x = Navin_Visvader_NORM_panCK_Reference)

# Merge subsets into recombined RDS
Visvader_0068_ER_Navin_Visvader_NORM_panCK_Reference <- merge(x = Visvader_0068_ER_panCK, y = Navin_Visvader_NORM_panCK_Reference)
# Join Layers
Visvader_0068_ER_Navin_Visvader_NORM_panCK_Reference <- JoinLayers(Visvader_0068_ER_Navin_Visvader_NORM_panCK_Reference)

# Remove subsetted objects
rm(Visvader_0068_ER_panCK)
rm(Navin_Visvader_NORM_panCK_Reference)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Visvader_0068_ER_Navin_Visvader_NORM_panCK_Reference <- FindVariableFeatures(Visvader_0068_ER_Navin_Visvader_NORM_panCK_Reference, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Visvader_0068_ER_Navin_Visvader_NORM_panCK_Reference)
Visvader_0068_ER_Navin_Visvader_NORM_panCK_Reference <- ScaleData(Visvader_0068_ER_Navin_Visvader_NORM_panCK_Reference, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Visvader_0068_ER_Navin_Visvader_NORM_panCK_Reference <- RunPCA(Visvader_0068_ER_Navin_Visvader_NORM_panCK_Reference, features = VariableFeatures(object = Visvader_0068_ER_Navin_Visvader_NORM_panCK_Reference))

# Examine and visualize PCA results a few different ways
print(Visvader_0068_ER_Navin_Visvader_NORM_panCK_Reference[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Visvader_0068_ER_Navin_Visvader_NORM_panCK_Reference, dims = 1:2, reduction = "pca")
DimPlot(Visvader_0068_ER_Navin_Visvader_NORM_panCK_Reference, reduction = "pca")

Visvader_0068_ER_Navin_Visvader_NORM_panCK_Reference <- FindNeighbors(Visvader_0068_ER_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
Visvader_0068_ER_Navin_Visvader_NORM_panCK_Reference <- FindClusters(Visvader_0068_ER_Navin_Visvader_NORM_panCK_Reference, resolution = 0.1)

# Umap clustering
Visvader_0068_ER_Navin_Visvader_NORM_panCK_Reference <- RunUMAP(Visvader_0068_ER_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
DimPlot(Visvader_0068_ER_Navin_Visvader_NORM_panCK_Reference, reduction = "umap")
DimPlot(Visvader_0068_ER_Navin_Visvader_NORM_panCK_Reference, reduction = "umap", group.by = "orig.ident")

# Save RDS
saveRDS(Visvader_0068_ER_Navin_Visvader_NORM_panCK_Reference, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Visvader_0068_ER_Navin_Visvader_NORM_panCK_Reference.rds")

########## Prepare Cell Annotations
# Extract & Save Meta Data
Visvader_0068_ER_Navin_Visvader_NORM_panCK_Reference.meta.data <- as.data.frame(as.matrix(Visvader_0068_ER_Navin_Visvader_NORM_panCK_Reference@meta.data))
# Extract pertinent columns
# Make row names first column
Visvader_0068_ER_Navin_Visvader_NORM_panCK_Reference.meta.data <- tibble::rownames_to_column(Visvader_0068_ER_Navin_Visvader_NORM_panCK_Reference.meta.data, "Cells")
# In Seurat, "X" becomes the cell names while "cnv.ident" is specified above
Visvader_0068_ER_Navin_Visvader_NORM_panCK_Reference.annotations <- Visvader_0068_ER_Navin_Visvader_NORM_panCK_Reference.meta.data[, c("Cells", "cnv.ident")]
# Delete header
names(Visvader_0068_ER_Navin_Visvader_NORM_panCK_Reference.annotations) <- NULL
# Save as table
write.table(Visvader_0068_ER_Navin_Visvader_NORM_panCK_Reference.annotations,"/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_0068_ER_Navin_Visvader_NORM_panCK_Reference.annotations.txt",sep="\t",row.names=FALSE)
# Remove metadata and annotations
rm(Visvader_0068_ER_Navin_Visvader_NORM_panCK_Reference.meta.data)
rm(Visvader_0068_ER_Navin_Visvader_NORM_panCK_Reference.annotations)

########## Generate Count Matrix
# Extract counts from the Seurat object
Visvader_0068_ER_Navin_Visvader_NORM_panCK_Reference.count <- Visvader_0068_ER_Navin_Visvader_NORM_panCK_Reference[["RNA"]]$counts
# Convert counts to a sparse matrix
Visvader_0068_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as(Visvader_0068_ER_Navin_Visvader_NORM_panCK_Reference.count, 'sparseMatrix')
rm(Visvader_0068_ER_Navin_Visvader_NORM_panCK_Reference.count)

# Convert to requisite data format
Visvader_0068_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as.data.frame(Visvader_0068_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix)
# Save table
write.csv(Visvader_0068_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_Count_Matrix/Visvader_0068_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv")
# Remove Objects
rm(Visvader_0068_ER_Navin_Visvader_NORM_panCK_Reference)
rm(Visvader_0068_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix)

########### Run inferCNV
# Note: check.names = FALSE prevents cell names from having dashes replaced with dots; essential to match annotations file
Visvader_0068_ER_inferCNV <- CreateInfercnvObject(raw_counts_matrix=read.csv(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_Count_Matrix/Visvader_0068_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv", header = TRUE, row.names = 1,check.names = FALSE),
                                                  annotations_file=read.table(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_0068_ER_Navin_Visvader_NORM_panCK_Reference.annotations.txt", header = FALSE, row.names = 1),
                                                  delim="\t",
                                                  gene_order_file="/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/hg38_gencode_v27.txt",
                                                  ref_group_names=c("Normal"))


# perform infercnv operations to reveal cnv signal
Visvader_0068_ER_inferCNV = infercnv::run(Visvader_0068_ER_inferCNV,
                                          cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                          out_dir="/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0068_ER_panCK/",  # dir is auto-created for storing outputs
                                          cluster_by_groups=T,   # cluster
                                          denoise=T,
                                          HMM=T,
                                          num_ref_groups = 8, analysis_mode = "subclusters")

###########  Integrate data with Seurat object
# Read RDS and add inferCNV results to metadata
Visvader_0068_ER_Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Visvader_0068_ER_Navin_Visvader_NORM_panCK_Reference.rds") 
Visvader_0068_ER_panCK_cnv <- infercnv::add_to_seurat(infercnv_output_path="/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0068_ER_panCK/",
                                                      seurat_obj=Visvader_0068_ER_Navin_Visvader_NORM_panCK_Reference, # optional
                                                      top_n=10)

# Collect Tumor cells from Seurat Object
Visvader_0068_ER_panCK_cnv <- subset(x = Visvader_0068_ER_panCK_cnv, subset = cnv.ident == "Tumor")
# Quantify chromosomes with alterations per cell
Visvader_0068_ER_panCK_cnv_Pos_Chromosomes <- as.data.frame(rowSums(Visvader_0068_ER_panCK_cnv@meta.data[,startsWith(names(Visvader_0068_ER_panCK_cnv@meta.data),"has_cnv_")] > 0 ))
Visvader_0068_ER_panCK_cnv_Neg_Chromosomes <- as.data.frame(rowSums(Visvader_0068_ER_panCK_cnv@meta.data[,startsWith(names(Visvader_0068_ER_panCK_cnv@meta.data),"has_cnv_")] == 0))
Visvader_0068_ER_panCK_cnv_Chromosome_Quant <- as.data.frame(c(Visvader_0068_ER_panCK_cnv_Pos_Chromosomes, Visvader_0068_ER_panCK_cnv_Neg_Chromosomes))
colnames(Visvader_0068_ER_panCK_cnv_Chromosome_Quant) <- c("Chromosomes CNV Pos", "Chromosomes CNV Neg")
write.csv(Visvader_0068_ER_panCK_cnv_Chromosome_Quant, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0068_ER_panCK/Visvader_0068_ER_panCK_cnv_Chromosome_Quant.csv")

# Save CNV Metadata
Visvader_0068_ER_panCK_cnv.meta.data <- as.data.frame(as.matrix(Visvader_0068_ER_panCK_cnv@meta.data))
write.csv(Visvader_0068_ER_panCK_cnv.meta.data, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0068_ER_panCK/Visvader_0068_ER_panCK_cnv.meta.data.csv")

# Subset Cells with inferred CNVs -------------------------------------------------------------------------------------------------------------
Visvader_0068_ER_panCK_cnv <- subset(Visvader_0068_ER_panCK_cnv, cells=rownames(Visvader_0068_ER_panCK_cnv@meta.data)[rowSums(Visvader_0068_ER_panCK_cnv@meta.data[,startsWith(names(Visvader_0068_ER_panCK_cnv@meta.data),"has_cnv_")]) > 0 ])

# Save RDS with inferredCNVs
saveRDS(Visvader_0068_ER_panCK_cnv, file = "/R/R_Visvader/Visvader_RDS/RDS_inferCNV/Visvader_0068_ER_panCK_cnv.rds")


# Free memory
rm(Visvader_0068_ER_inferCNV)
rm(Visvader_0068_ER_Navin_Visvader_NORM_panCK_Reference)
rm(Visvader_0068_ER_panCK_cnv)
rm(Visvader_0068_ER_panCK_cnv_Pos_Chromosomes)
rm(Visvader_0068_ER_panCK_cnv_Neg_Chromosomes)
rm(Visvader_0068_ER_panCK_cnv_Chromosome_Quant)
rm(Visvader_0068_ER_panCK_cnv.meta.data)
gc()



##########################################
#### InferCNV: Visvader_0114_ER_panCK ####
##########################################

########## Merge with normal mammary glad, extract epithelium, and process via Seurat
# Load RDS
Visvader_0114_ER_panCK <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0114_ER_panCK.rds")
Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Navin_Visvader_NORM_panCK_Reference.rds")

# Set identities 
colnames(x = Visvader_0114_ER_panCK[[]])
Idents(object = Visvader_0114_ER_panCK) <- 'Tumor'
Visvader_0114_ER_panCK[["cnv.ident"]] <- Idents(object = Visvader_0114_ER_panCK)
levels(x = Visvader_0114_ER_panCK)

colnames(x = Navin_Visvader_NORM_panCK_Reference[[]])
Idents(object = Navin_Visvader_NORM_panCK_Reference) <- 'Normal'
Navin_Visvader_NORM_panCK_Reference[["cnv.ident"]] <- Idents(object = Navin_Visvader_NORM_panCK_Reference)
levels(x = Navin_Visvader_NORM_panCK_Reference)

# Merge subsets into recombined RDS
Visvader_0114_ER_Navin_Visvader_NORM_panCK_Reference <- merge(x = Visvader_0114_ER_panCK, y = Navin_Visvader_NORM_panCK_Reference)
# Join Layers
Visvader_0114_ER_Navin_Visvader_NORM_panCK_Reference <- JoinLayers(Visvader_0114_ER_Navin_Visvader_NORM_panCK_Reference)

# Remove subsetted objects
rm(Visvader_0114_ER_panCK)
rm(Navin_Visvader_NORM_panCK_Reference)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Visvader_0114_ER_Navin_Visvader_NORM_panCK_Reference <- FindVariableFeatures(Visvader_0114_ER_Navin_Visvader_NORM_panCK_Reference, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Visvader_0114_ER_Navin_Visvader_NORM_panCK_Reference)
Visvader_0114_ER_Navin_Visvader_NORM_panCK_Reference <- ScaleData(Visvader_0114_ER_Navin_Visvader_NORM_panCK_Reference, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Visvader_0114_ER_Navin_Visvader_NORM_panCK_Reference <- RunPCA(Visvader_0114_ER_Navin_Visvader_NORM_panCK_Reference, features = VariableFeatures(object = Visvader_0114_ER_Navin_Visvader_NORM_panCK_Reference))

# Examine and visualize PCA results a few different ways
print(Visvader_0114_ER_Navin_Visvader_NORM_panCK_Reference[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Visvader_0114_ER_Navin_Visvader_NORM_panCK_Reference, dims = 1:2, reduction = "pca")
DimPlot(Visvader_0114_ER_Navin_Visvader_NORM_panCK_Reference, reduction = "pca")

Visvader_0114_ER_Navin_Visvader_NORM_panCK_Reference <- FindNeighbors(Visvader_0114_ER_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
Visvader_0114_ER_Navin_Visvader_NORM_panCK_Reference <- FindClusters(Visvader_0114_ER_Navin_Visvader_NORM_panCK_Reference, resolution = 0.1)

# Umap clustering
Visvader_0114_ER_Navin_Visvader_NORM_panCK_Reference <- RunUMAP(Visvader_0114_ER_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
DimPlot(Visvader_0114_ER_Navin_Visvader_NORM_panCK_Reference, reduction = "umap")
DimPlot(Visvader_0114_ER_Navin_Visvader_NORM_panCK_Reference, reduction = "umap", group.by = "orig.ident")

# Save RDS
saveRDS(Visvader_0114_ER_Navin_Visvader_NORM_panCK_Reference, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Visvader_0114_ER_Navin_Visvader_NORM_panCK_Reference.rds")

########## Prepare Cell Annotations
# Extract & Save Meta Data
Visvader_0114_ER_Navin_Visvader_NORM_panCK_Reference.meta.data <- as.data.frame(as.matrix(Visvader_0114_ER_Navin_Visvader_NORM_panCK_Reference@meta.data))
# Extract pertinent columns
# Make row names first column
Visvader_0114_ER_Navin_Visvader_NORM_panCK_Reference.meta.data <- tibble::rownames_to_column(Visvader_0114_ER_Navin_Visvader_NORM_panCK_Reference.meta.data, "Cells")
# In Seurat, "X" becomes the cell names while "cnv.ident" is specified above
Visvader_0114_ER_Navin_Visvader_NORM_panCK_Reference.annotations <- Visvader_0114_ER_Navin_Visvader_NORM_panCK_Reference.meta.data[, c("Cells", "cnv.ident")]
# Delete header
names(Visvader_0114_ER_Navin_Visvader_NORM_panCK_Reference.annotations) <- NULL
# Save as table
write.table(Visvader_0114_ER_Navin_Visvader_NORM_panCK_Reference.annotations,"/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_0114_ER_Navin_Visvader_NORM_panCK_Reference.annotations.txt",sep="\t",row.names=FALSE)
# Remove metadata and annotations
rm(Visvader_0114_ER_Navin_Visvader_NORM_panCK_Reference.meta.data)
rm(Visvader_0114_ER_Navin_Visvader_NORM_panCK_Reference.annotations)

########## Generate Count Matrix
# Extract counts from the Seurat object
Visvader_0114_ER_Navin_Visvader_NORM_panCK_Reference.count <- Visvader_0114_ER_Navin_Visvader_NORM_panCK_Reference[["RNA"]]$counts
# Convert counts to a sparse matrix
Visvader_0114_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as(Visvader_0114_ER_Navin_Visvader_NORM_panCK_Reference.count, 'sparseMatrix')
rm(Visvader_0114_ER_Navin_Visvader_NORM_panCK_Reference.count)

# Convert to requisite data format
Visvader_0114_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as.data.frame(Visvader_0114_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix)
# Save table
write.csv(Visvader_0114_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_Count_Matrix/Visvader_0114_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv")
# Remove Objects
rm(Visvader_0114_ER_Navin_Visvader_NORM_panCK_Reference)
rm(Visvader_0114_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix)

########### Run inferCNV
# Note: check.names = FALSE prevents cell names from having dashes replaced with dots; essential to match annotations file
Visvader_0114_ER_inferCNV <- CreateInfercnvObject(raw_counts_matrix=read.csv(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_Count_Matrix/Visvader_0114_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv", header = TRUE, row.names = 1,check.names = FALSE),
                                                  annotations_file=read.table(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_0114_ER_Navin_Visvader_NORM_panCK_Reference.annotations.txt", header = FALSE, row.names = 1),
                                                  delim="\t",
                                                  gene_order_file="/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/hg38_gencode_v27.txt",
                                                  ref_group_names=c("Normal"))


# perform infercnv operations to reveal cnv signal
Visvader_0114_ER_inferCNV = infercnv::run(Visvader_0114_ER_inferCNV,
                                          cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                          out_dir="/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0114_ER_panCK/",  # dir is auto-created for storing outputs
                                          cluster_by_groups=T,   # cluster
                                          denoise=T,
                                          HMM=T,
                                          num_ref_groups = 8, analysis_mode = "subclusters")

###########  Integrate data with Seurat object
# Read RDS and add inferCNV results to metadata
Visvader_0114_ER_Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Visvader_0114_ER_Navin_Visvader_NORM_panCK_Reference.rds") 
Visvader_0114_ER_panCK_cnv <- infercnv::add_to_seurat(infercnv_output_path="/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0114_ER_panCK/",
                                                      seurat_obj=Visvader_0114_ER_Navin_Visvader_NORM_panCK_Reference, # optional
                                                      top_n=10)

# Collect Tumor cells from Seurat Object
Visvader_0114_ER_panCK_cnv <- subset(x = Visvader_0114_ER_panCK_cnv, subset = cnv.ident == "Tumor")
# Quantify chromosomes with alterations per cell
Visvader_0114_ER_panCK_cnv_Pos_Chromosomes <- as.data.frame(rowSums(Visvader_0114_ER_panCK_cnv@meta.data[,startsWith(names(Visvader_0114_ER_panCK_cnv@meta.data),"has_cnv_")] > 0 ))
Visvader_0114_ER_panCK_cnv_Neg_Chromosomes <- as.data.frame(rowSums(Visvader_0114_ER_panCK_cnv@meta.data[,startsWith(names(Visvader_0114_ER_panCK_cnv@meta.data),"has_cnv_")] == 0))
Visvader_0114_ER_panCK_cnv_Chromosome_Quant <- as.data.frame(c(Visvader_0114_ER_panCK_cnv_Pos_Chromosomes, Visvader_0114_ER_panCK_cnv_Neg_Chromosomes))
colnames(Visvader_0114_ER_panCK_cnv_Chromosome_Quant) <- c("Chromosomes CNV Pos", "Chromosomes CNV Neg")
write.csv(Visvader_0114_ER_panCK_cnv_Chromosome_Quant, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0114_ER_panCK/Visvader_0114_ER_panCK_cnv_Chromosome_Quant.csv")

# Save CNV Metadata
Visvader_0114_ER_panCK_cnv.meta.data <- as.data.frame(as.matrix(Visvader_0114_ER_panCK_cnv@meta.data))
write.csv(Visvader_0114_ER_panCK_cnv.meta.data, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0114_ER_panCK/Visvader_0114_ER_panCK_cnv.meta.data.csv")

# Subset Cells with inferred CNVs -------------------------------------------------------------------------------------------------------------
Visvader_0114_ER_panCK_cnv <- subset(Visvader_0114_ER_panCK_cnv, cells=rownames(Visvader_0114_ER_panCK_cnv@meta.data)[rowSums(Visvader_0114_ER_panCK_cnv@meta.data[,startsWith(names(Visvader_0114_ER_panCK_cnv@meta.data),"has_cnv_")]) > 0 ])

# Save RDS with inferredCNVs
saveRDS(Visvader_0114_ER_panCK_cnv, file = "/R/R_Visvader/Visvader_RDS/RDS_inferCNV/Visvader_0114_ER_panCK_cnv.rds")


# Free memory
rm(Visvader_0114_ER_inferCNV)
rm(Visvader_0114_ER_Navin_Visvader_NORM_panCK_Reference)
rm(Visvader_0114_ER_panCK_cnv)
rm(Visvader_0114_ER_panCK_cnv_Pos_Chromosomes)
rm(Visvader_0114_ER_panCK_cnv_Neg_Chromosomes)
rm(Visvader_0114_ER_panCK_cnv_Chromosome_Quant)
rm(Visvader_0114_ER_panCK_cnv.meta.data)
gc()



##########################################
#### InferCNV: Visvader_0125_ER_panCK ####
##########################################

########## Merge with normal mammary glad, extract epithelium, and process via Seurat
# Load RDS
Visvader_0125_ER_panCK <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0125_ER_panCK.rds")
Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Navin_Visvader_NORM_panCK_Reference.rds")

# Set identities 
colnames(x = Visvader_0125_ER_panCK[[]])
Idents(object = Visvader_0125_ER_panCK) <- 'Tumor'
Visvader_0125_ER_panCK[["cnv.ident"]] <- Idents(object = Visvader_0125_ER_panCK)
levels(x = Visvader_0125_ER_panCK)

colnames(x = Navin_Visvader_NORM_panCK_Reference[[]])
Idents(object = Navin_Visvader_NORM_panCK_Reference) <- 'Normal'
Navin_Visvader_NORM_panCK_Reference[["cnv.ident"]] <- Idents(object = Navin_Visvader_NORM_panCK_Reference)
levels(x = Navin_Visvader_NORM_panCK_Reference)

# Merge subsets into recombined RDS
Visvader_0125_ER_Navin_Visvader_NORM_panCK_Reference <- merge(x = Visvader_0125_ER_panCK, y = Navin_Visvader_NORM_panCK_Reference)
# Join Layers
Visvader_0125_ER_Navin_Visvader_NORM_panCK_Reference <- JoinLayers(Visvader_0125_ER_Navin_Visvader_NORM_panCK_Reference)

# Remove subsetted objects
rm(Visvader_0125_ER_panCK)
rm(Navin_Visvader_NORM_panCK_Reference)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Visvader_0125_ER_Navin_Visvader_NORM_panCK_Reference <- FindVariableFeatures(Visvader_0125_ER_Navin_Visvader_NORM_panCK_Reference, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Visvader_0125_ER_Navin_Visvader_NORM_panCK_Reference)
Visvader_0125_ER_Navin_Visvader_NORM_panCK_Reference <- ScaleData(Visvader_0125_ER_Navin_Visvader_NORM_panCK_Reference, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Visvader_0125_ER_Navin_Visvader_NORM_panCK_Reference <- RunPCA(Visvader_0125_ER_Navin_Visvader_NORM_panCK_Reference, features = VariableFeatures(object = Visvader_0125_ER_Navin_Visvader_NORM_panCK_Reference))

# Examine and visualize PCA results a few different ways
print(Visvader_0125_ER_Navin_Visvader_NORM_panCK_Reference[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Visvader_0125_ER_Navin_Visvader_NORM_panCK_Reference, dims = 1:2, reduction = "pca")
DimPlot(Visvader_0125_ER_Navin_Visvader_NORM_panCK_Reference, reduction = "pca")

Visvader_0125_ER_Navin_Visvader_NORM_panCK_Reference <- FindNeighbors(Visvader_0125_ER_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
Visvader_0125_ER_Navin_Visvader_NORM_panCK_Reference <- FindClusters(Visvader_0125_ER_Navin_Visvader_NORM_panCK_Reference, resolution = 0.1)

# Umap clustering
Visvader_0125_ER_Navin_Visvader_NORM_panCK_Reference <- RunUMAP(Visvader_0125_ER_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
DimPlot(Visvader_0125_ER_Navin_Visvader_NORM_panCK_Reference, reduction = "umap")
DimPlot(Visvader_0125_ER_Navin_Visvader_NORM_panCK_Reference, reduction = "umap", group.by = "orig.ident")

# Save RDS
saveRDS(Visvader_0125_ER_Navin_Visvader_NORM_panCK_Reference, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Visvader_0125_ER_Navin_Visvader_NORM_panCK_Reference.rds")

########## Prepare Cell Annotations
# Extract & Save Meta Data
Visvader_0125_ER_Navin_Visvader_NORM_panCK_Reference.meta.data <- as.data.frame(as.matrix(Visvader_0125_ER_Navin_Visvader_NORM_panCK_Reference@meta.data))
# Extract pertinent columns
# Make row names first column
Visvader_0125_ER_Navin_Visvader_NORM_panCK_Reference.meta.data <- tibble::rownames_to_column(Visvader_0125_ER_Navin_Visvader_NORM_panCK_Reference.meta.data, "Cells")
# In Seurat, "X" becomes the cell names while "cnv.ident" is specified above
Visvader_0125_ER_Navin_Visvader_NORM_panCK_Reference.annotations <- Visvader_0125_ER_Navin_Visvader_NORM_panCK_Reference.meta.data[, c("Cells", "cnv.ident")]
# Delete header
names(Visvader_0125_ER_Navin_Visvader_NORM_panCK_Reference.annotations) <- NULL
# Save as table
write.table(Visvader_0125_ER_Navin_Visvader_NORM_panCK_Reference.annotations,"/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_0125_ER_Navin_Visvader_NORM_panCK_Reference.annotations.txt",sep="\t",row.names=FALSE)
# Remove metadata and annotations
rm(Visvader_0125_ER_Navin_Visvader_NORM_panCK_Reference.meta.data)
rm(Visvader_0125_ER_Navin_Visvader_NORM_panCK_Reference.annotations)

########## Generate Count Matrix
# Extract counts from the Seurat object
Visvader_0125_ER_Navin_Visvader_NORM_panCK_Reference.count <- Visvader_0125_ER_Navin_Visvader_NORM_panCK_Reference[["RNA"]]$counts
# Convert counts to a sparse matrix
Visvader_0125_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as(Visvader_0125_ER_Navin_Visvader_NORM_panCK_Reference.count, 'sparseMatrix')
rm(Visvader_0125_ER_Navin_Visvader_NORM_panCK_Reference.count)

# Convert to requisite data format
Visvader_0125_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as.data.frame(Visvader_0125_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix)
# Save table
write.csv(Visvader_0125_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_Count_Matrix/Visvader_0125_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv")
# Remove Objects
rm(Visvader_0125_ER_Navin_Visvader_NORM_panCK_Reference)
rm(Visvader_0125_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix)

########### Run inferCNV
# Note: check.names = FALSE prevents cell names from having dashes replaced with dots; essential to match annotations file
Visvader_0125_ER_inferCNV <- CreateInfercnvObject(raw_counts_matrix=read.csv(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_Count_Matrix/Visvader_0125_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv", header = TRUE, row.names = 1,check.names = FALSE),
                                                  annotations_file=read.table(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_0125_ER_Navin_Visvader_NORM_panCK_Reference.annotations.txt", header = FALSE, row.names = 1),
                                                  delim="\t",
                                                  gene_order_file="/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/hg38_gencode_v27.txt",
                                                  ref_group_names=c("Normal"))


# perform infercnv operations to reveal cnv signal
Visvader_0125_ER_inferCNV = infercnv::run(Visvader_0125_ER_inferCNV,
                                          cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                          out_dir="/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0125_ER_panCK/",  # dir is auto-created for storing outputs
                                          cluster_by_groups=T,   # cluster
                                          denoise=T,
                                          HMM=T,
                                          num_ref_groups = 8, analysis_mode = "subclusters")

###########  Integrate data with Seurat object
# Read RDS and add inferCNV results to metadata
Visvader_0125_ER_Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Visvader_0125_ER_Navin_Visvader_NORM_panCK_Reference.rds") 
Visvader_0125_ER_panCK_cnv <- infercnv::add_to_seurat(infercnv_output_path="/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0125_ER_panCK/",
                                                      seurat_obj=Visvader_0125_ER_Navin_Visvader_NORM_panCK_Reference, # optional
                                                      top_n=10)

# Collect Tumor cells from Seurat Object
Visvader_0125_ER_panCK_cnv <- subset(x = Visvader_0125_ER_panCK_cnv, subset = cnv.ident == "Tumor")
# Quantify chromosomes with alterations per cell
Visvader_0125_ER_panCK_cnv_Pos_Chromosomes <- as.data.frame(rowSums(Visvader_0125_ER_panCK_cnv@meta.data[,startsWith(names(Visvader_0125_ER_panCK_cnv@meta.data),"has_cnv_")] > 0 ))
Visvader_0125_ER_panCK_cnv_Neg_Chromosomes <- as.data.frame(rowSums(Visvader_0125_ER_panCK_cnv@meta.data[,startsWith(names(Visvader_0125_ER_panCK_cnv@meta.data),"has_cnv_")] == 0))
Visvader_0125_ER_panCK_cnv_Chromosome_Quant <- as.data.frame(c(Visvader_0125_ER_panCK_cnv_Pos_Chromosomes, Visvader_0125_ER_panCK_cnv_Neg_Chromosomes))
colnames(Visvader_0125_ER_panCK_cnv_Chromosome_Quant) <- c("Chromosomes CNV Pos", "Chromosomes CNV Neg")
write.csv(Visvader_0125_ER_panCK_cnv_Chromosome_Quant, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0125_ER_panCK/Visvader_0125_ER_panCK_cnv_Chromosome_Quant.csv")

# Save CNV Metadata
Visvader_0125_ER_panCK_cnv.meta.data <- as.data.frame(as.matrix(Visvader_0125_ER_panCK_cnv@meta.data))
write.csv(Visvader_0125_ER_panCK_cnv.meta.data, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0125_ER_panCK/Visvader_0125_ER_panCK_cnv.meta.data.csv")

# Subset Cells with inferred CNVs -------------------------------------------------------------------------------------------------------------
Visvader_0125_ER_panCK_cnv <- subset(Visvader_0125_ER_panCK_cnv, cells=rownames(Visvader_0125_ER_panCK_cnv@meta.data)[rowSums(Visvader_0125_ER_panCK_cnv@meta.data[,startsWith(names(Visvader_0125_ER_panCK_cnv@meta.data),"has_cnv_")]) > 0 ])

# Save RDS with inferredCNVs
saveRDS(Visvader_0125_ER_panCK_cnv, file = "/R/R_Visvader/Visvader_RDS/RDS_inferCNV/Visvader_0125_ER_panCK_cnv.rds")


# Free memory
rm(Visvader_0125_ER_inferCNV)
rm(Visvader_0125_ER_Navin_Visvader_NORM_panCK_Reference)
rm(Visvader_0125_ER_panCK_cnv)
rm(Visvader_0125_ER_panCK_cnv_Pos_Chromosomes)
rm(Visvader_0125_ER_panCK_cnv_Neg_Chromosomes)
rm(Visvader_0125_ER_panCK_cnv_Chromosome_Quant)
rm(Visvader_0125_ER_panCK_cnv.meta.data)
gc()



##########################################
#### InferCNV: Visvader_0151_ER_panCK ####
##########################################

########## Merge with normal mammary glad, extract epithelium, and process via Seurat
# Load RDS
Visvader_0151_ER_panCK <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0151_ER_panCK.rds")
Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Navin_Visvader_NORM_panCK_Reference.rds")

# Set identities 
colnames(x = Visvader_0151_ER_panCK[[]])
Idents(object = Visvader_0151_ER_panCK) <- 'Tumor'
Visvader_0151_ER_panCK[["cnv.ident"]] <- Idents(object = Visvader_0151_ER_panCK)
levels(x = Visvader_0151_ER_panCK)

colnames(x = Navin_Visvader_NORM_panCK_Reference[[]])
Idents(object = Navin_Visvader_NORM_panCK_Reference) <- 'Normal'
Navin_Visvader_NORM_panCK_Reference[["cnv.ident"]] <- Idents(object = Navin_Visvader_NORM_panCK_Reference)
levels(x = Navin_Visvader_NORM_panCK_Reference)

# Merge subsets into recombined RDS
Visvader_0151_ER_Navin_Visvader_NORM_panCK_Reference <- merge(x = Visvader_0151_ER_panCK, y = Navin_Visvader_NORM_panCK_Reference)
# Join Layers
Visvader_0151_ER_Navin_Visvader_NORM_panCK_Reference <- JoinLayers(Visvader_0151_ER_Navin_Visvader_NORM_panCK_Reference)

# Remove subsetted objects
rm(Visvader_0151_ER_panCK)
rm(Navin_Visvader_NORM_panCK_Reference)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Visvader_0151_ER_Navin_Visvader_NORM_panCK_Reference <- FindVariableFeatures(Visvader_0151_ER_Navin_Visvader_NORM_panCK_Reference, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Visvader_0151_ER_Navin_Visvader_NORM_panCK_Reference)
Visvader_0151_ER_Navin_Visvader_NORM_panCK_Reference <- ScaleData(Visvader_0151_ER_Navin_Visvader_NORM_panCK_Reference, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Visvader_0151_ER_Navin_Visvader_NORM_panCK_Reference <- RunPCA(Visvader_0151_ER_Navin_Visvader_NORM_panCK_Reference, features = VariableFeatures(object = Visvader_0151_ER_Navin_Visvader_NORM_panCK_Reference))

# Examine and visualize PCA results a few different ways
print(Visvader_0151_ER_Navin_Visvader_NORM_panCK_Reference[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Visvader_0151_ER_Navin_Visvader_NORM_panCK_Reference, dims = 1:2, reduction = "pca")
DimPlot(Visvader_0151_ER_Navin_Visvader_NORM_panCK_Reference, reduction = "pca")

Visvader_0151_ER_Navin_Visvader_NORM_panCK_Reference <- FindNeighbors(Visvader_0151_ER_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
Visvader_0151_ER_Navin_Visvader_NORM_panCK_Reference <- FindClusters(Visvader_0151_ER_Navin_Visvader_NORM_panCK_Reference, resolution = 0.1)

# Umap clustering
Visvader_0151_ER_Navin_Visvader_NORM_panCK_Reference <- RunUMAP(Visvader_0151_ER_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
DimPlot(Visvader_0151_ER_Navin_Visvader_NORM_panCK_Reference, reduction = "umap")
DimPlot(Visvader_0151_ER_Navin_Visvader_NORM_panCK_Reference, reduction = "umap", group.by = "orig.ident")

# Save RDS
saveRDS(Visvader_0151_ER_Navin_Visvader_NORM_panCK_Reference, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Visvader_0151_ER_Navin_Visvader_NORM_panCK_Reference.rds")

########## Prepare Cell Annotations
# Extract & Save Meta Data
Visvader_0151_ER_Navin_Visvader_NORM_panCK_Reference.meta.data <- as.data.frame(as.matrix(Visvader_0151_ER_Navin_Visvader_NORM_panCK_Reference@meta.data))
# Extract pertinent columns
# Make row names first column
Visvader_0151_ER_Navin_Visvader_NORM_panCK_Reference.meta.data <- tibble::rownames_to_column(Visvader_0151_ER_Navin_Visvader_NORM_panCK_Reference.meta.data, "Cells")
# In Seurat, "X" becomes the cell names while "cnv.ident" is specified above
Visvader_0151_ER_Navin_Visvader_NORM_panCK_Reference.annotations <- Visvader_0151_ER_Navin_Visvader_NORM_panCK_Reference.meta.data[, c("Cells", "cnv.ident")]
# Delete header
names(Visvader_0151_ER_Navin_Visvader_NORM_panCK_Reference.annotations) <- NULL
# Save as table
write.table(Visvader_0151_ER_Navin_Visvader_NORM_panCK_Reference.annotations,"/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_0151_ER_Navin_Visvader_NORM_panCK_Reference.annotations.txt",sep="\t",row.names=FALSE)
# Remove metadata and annotations
rm(Visvader_0151_ER_Navin_Visvader_NORM_panCK_Reference.meta.data)
rm(Visvader_0151_ER_Navin_Visvader_NORM_panCK_Reference.annotations)

########## Generate Count Matrix
# Extract counts from the Seurat object
Visvader_0151_ER_Navin_Visvader_NORM_panCK_Reference.count <- Visvader_0151_ER_Navin_Visvader_NORM_panCK_Reference[["RNA"]]$counts
# Convert counts to a sparse matrix
Visvader_0151_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as(Visvader_0151_ER_Navin_Visvader_NORM_panCK_Reference.count, 'sparseMatrix')
rm(Visvader_0151_ER_Navin_Visvader_NORM_panCK_Reference.count)

# Convert to requisite data format
Visvader_0151_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as.data.frame(Visvader_0151_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix)
# Save table
write.csv(Visvader_0151_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_Count_Matrix/Visvader_0151_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv")
# Remove Objects
rm(Visvader_0151_ER_Navin_Visvader_NORM_panCK_Reference)
rm(Visvader_0151_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix)

########### Run inferCNV
# Note: check.names = FALSE prevents cell names from having dashes replaced with dots; essential to match annotations file
Visvader_0151_ER_inferCNV <- CreateInfercnvObject(raw_counts_matrix=read.csv(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_Count_Matrix/Visvader_0151_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv", header = TRUE, row.names = 1,check.names = FALSE),
                                                  annotations_file=read.table(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_0151_ER_Navin_Visvader_NORM_panCK_Reference.annotations.txt", header = FALSE, row.names = 1),
                                                  delim="\t",
                                                  gene_order_file="/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/hg38_gencode_v27.txt",
                                                  ref_group_names=c("Normal"))


# perform infercnv operations to reveal cnv signal
Visvader_0151_ER_inferCNV = infercnv::run(Visvader_0151_ER_inferCNV,
                                          cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                          out_dir="/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0151_ER_panCK/",  # dir is auto-created for storing outputs
                                          cluster_by_groups=T,   # cluster
                                          denoise=T,
                                          HMM=T,
                                          num_ref_groups = 8, analysis_mode = "subclusters")

###########  Integrate data with Seurat object
# Read RDS and add inferCNV results to metadata
Visvader_0151_ER_Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Visvader_0151_ER_Navin_Visvader_NORM_panCK_Reference.rds") 
Visvader_0151_ER_panCK_cnv <- infercnv::add_to_seurat(infercnv_output_path="/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0151_ER_panCK/",
                                                      seurat_obj=Visvader_0151_ER_Navin_Visvader_NORM_panCK_Reference, # optional
                                                      top_n=10)

# Collect Tumor cells from Seurat Object
Visvader_0151_ER_panCK_cnv <- subset(x = Visvader_0151_ER_panCK_cnv, subset = cnv.ident == "Tumor")
# Quantify chromosomes with alterations per cell
Visvader_0151_ER_panCK_cnv_Pos_Chromosomes <- as.data.frame(rowSums(Visvader_0151_ER_panCK_cnv@meta.data[,startsWith(names(Visvader_0151_ER_panCK_cnv@meta.data),"has_cnv_")] > 0 ))
Visvader_0151_ER_panCK_cnv_Neg_Chromosomes <- as.data.frame(rowSums(Visvader_0151_ER_panCK_cnv@meta.data[,startsWith(names(Visvader_0151_ER_panCK_cnv@meta.data),"has_cnv_")] == 0))
Visvader_0151_ER_panCK_cnv_Chromosome_Quant <- as.data.frame(c(Visvader_0151_ER_panCK_cnv_Pos_Chromosomes, Visvader_0151_ER_panCK_cnv_Neg_Chromosomes))
colnames(Visvader_0151_ER_panCK_cnv_Chromosome_Quant) <- c("Chromosomes CNV Pos", "Chromosomes CNV Neg")
write.csv(Visvader_0151_ER_panCK_cnv_Chromosome_Quant, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0151_ER_panCK/Visvader_0151_ER_panCK_cnv_Chromosome_Quant.csv")

# Save CNV Metadata
Visvader_0151_ER_panCK_cnv.meta.data <- as.data.frame(as.matrix(Visvader_0151_ER_panCK_cnv@meta.data))
write.csv(Visvader_0151_ER_panCK_cnv.meta.data, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0151_ER_panCK/Visvader_0151_ER_panCK_cnv.meta.data.csv")

# Subset Cells with inferred CNVs -------------------------------------------------------------------------------------------------------------
Visvader_0151_ER_panCK_cnv <- subset(Visvader_0151_ER_panCK_cnv, cells=rownames(Visvader_0151_ER_panCK_cnv@meta.data)[rowSums(Visvader_0151_ER_panCK_cnv@meta.data[,startsWith(names(Visvader_0151_ER_panCK_cnv@meta.data),"has_cnv_")]) > 0 ])

# Save RDS with inferredCNVs
saveRDS(Visvader_0151_ER_panCK_cnv, file = "/R/R_Visvader/Visvader_RDS/RDS_inferCNV/Visvader_0151_ER_panCK_cnv.rds")


# Free memory
rm(Visvader_0151_ER_inferCNV)
rm(Visvader_0151_ER_Navin_Visvader_NORM_panCK_Reference)
rm(Visvader_0151_ER_panCK_cnv)
rm(Visvader_0151_ER_panCK_cnv_Pos_Chromosomes)
rm(Visvader_0151_ER_panCK_cnv_Neg_Chromosomes)
rm(Visvader_0151_ER_panCK_cnv_Chromosome_Quant)
rm(Visvader_0151_ER_panCK_cnv.meta.data)
gc()



##########################################
#### InferCNV: Visvader_0163_ER_panCK ####
##########################################

########## Merge with normal mammary glad, extract epithelium, and process via Seurat
# Load RDS
Visvader_0163_ER_panCK <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0163_ER_panCK.rds")
Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Navin_Visvader_NORM_panCK_Reference.rds")

# Set identities 
colnames(x = Visvader_0163_ER_panCK[[]])
Idents(object = Visvader_0163_ER_panCK) <- 'Tumor'
Visvader_0163_ER_panCK[["cnv.ident"]] <- Idents(object = Visvader_0163_ER_panCK)
levels(x = Visvader_0163_ER_panCK)

colnames(x = Navin_Visvader_NORM_panCK_Reference[[]])
Idents(object = Navin_Visvader_NORM_panCK_Reference) <- 'Normal'
Navin_Visvader_NORM_panCK_Reference[["cnv.ident"]] <- Idents(object = Navin_Visvader_NORM_panCK_Reference)
levels(x = Navin_Visvader_NORM_panCK_Reference)

# Merge subsets into recombined RDS
Visvader_0163_ER_Navin_Visvader_NORM_panCK_Reference <- merge(x = Visvader_0163_ER_panCK, y = Navin_Visvader_NORM_panCK_Reference)
# Join Layers
Visvader_0163_ER_Navin_Visvader_NORM_panCK_Reference <- JoinLayers(Visvader_0163_ER_Navin_Visvader_NORM_panCK_Reference)

# Remove subsetted objects
rm(Visvader_0163_ER_panCK)
rm(Navin_Visvader_NORM_panCK_Reference)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Visvader_0163_ER_Navin_Visvader_NORM_panCK_Reference <- FindVariableFeatures(Visvader_0163_ER_Navin_Visvader_NORM_panCK_Reference, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Visvader_0163_ER_Navin_Visvader_NORM_panCK_Reference)
Visvader_0163_ER_Navin_Visvader_NORM_panCK_Reference <- ScaleData(Visvader_0163_ER_Navin_Visvader_NORM_panCK_Reference, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Visvader_0163_ER_Navin_Visvader_NORM_panCK_Reference <- RunPCA(Visvader_0163_ER_Navin_Visvader_NORM_panCK_Reference, features = VariableFeatures(object = Visvader_0163_ER_Navin_Visvader_NORM_panCK_Reference))

# Examine and visualize PCA results a few different ways
print(Visvader_0163_ER_Navin_Visvader_NORM_panCK_Reference[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Visvader_0163_ER_Navin_Visvader_NORM_panCK_Reference, dims = 1:2, reduction = "pca")
DimPlot(Visvader_0163_ER_Navin_Visvader_NORM_panCK_Reference, reduction = "pca")

Visvader_0163_ER_Navin_Visvader_NORM_panCK_Reference <- FindNeighbors(Visvader_0163_ER_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
Visvader_0163_ER_Navin_Visvader_NORM_panCK_Reference <- FindClusters(Visvader_0163_ER_Navin_Visvader_NORM_panCK_Reference, resolution = 0.1)

# Umap clustering
Visvader_0163_ER_Navin_Visvader_NORM_panCK_Reference <- RunUMAP(Visvader_0163_ER_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
DimPlot(Visvader_0163_ER_Navin_Visvader_NORM_panCK_Reference, reduction = "umap")
DimPlot(Visvader_0163_ER_Navin_Visvader_NORM_panCK_Reference, reduction = "umap", group.by = "orig.ident")

# Save RDS
saveRDS(Visvader_0163_ER_Navin_Visvader_NORM_panCK_Reference, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Visvader_0163_ER_Navin_Visvader_NORM_panCK_Reference.rds")

########## Prepare Cell Annotations
# Extract & Save Meta Data
Visvader_0163_ER_Navin_Visvader_NORM_panCK_Reference.meta.data <- as.data.frame(as.matrix(Visvader_0163_ER_Navin_Visvader_NORM_panCK_Reference@meta.data))
# Extract pertinent columns
# Make row names first column
Visvader_0163_ER_Navin_Visvader_NORM_panCK_Reference.meta.data <- tibble::rownames_to_column(Visvader_0163_ER_Navin_Visvader_NORM_panCK_Reference.meta.data, "Cells")
# In Seurat, "X" becomes the cell names while "cnv.ident" is specified above
Visvader_0163_ER_Navin_Visvader_NORM_panCK_Reference.annotations <- Visvader_0163_ER_Navin_Visvader_NORM_panCK_Reference.meta.data[, c("Cells", "cnv.ident")]
# Delete header
names(Visvader_0163_ER_Navin_Visvader_NORM_panCK_Reference.annotations) <- NULL
# Save as table
write.table(Visvader_0163_ER_Navin_Visvader_NORM_panCK_Reference.annotations,"/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_0163_ER_Navin_Visvader_NORM_panCK_Reference.annotations.txt",sep="\t",row.names=FALSE)
# Remove metadata and annotations
rm(Visvader_0163_ER_Navin_Visvader_NORM_panCK_Reference.meta.data)
rm(Visvader_0163_ER_Navin_Visvader_NORM_panCK_Reference.annotations)

########## Generate Count Matrix
# Extract counts from the Seurat object
Visvader_0163_ER_Navin_Visvader_NORM_panCK_Reference.count <- Visvader_0163_ER_Navin_Visvader_NORM_panCK_Reference[["RNA"]]$counts
# Convert counts to a sparse matrix
Visvader_0163_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as(Visvader_0163_ER_Navin_Visvader_NORM_panCK_Reference.count, 'sparseMatrix')
rm(Visvader_0163_ER_Navin_Visvader_NORM_panCK_Reference.count)

# Convert to requisite data format
Visvader_0163_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as.data.frame(Visvader_0163_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix)
# Save table
write.csv(Visvader_0163_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_Count_Matrix/Visvader_0163_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv")
# Remove Objects
rm(Visvader_0163_ER_Navin_Visvader_NORM_panCK_Reference)
rm(Visvader_0163_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix)

########### Run inferCNV
# Note: check.names = FALSE prevents cell names from having dashes replaced with dots; essential to match annotations file
Visvader_0163_ER_inferCNV <- CreateInfercnvObject(raw_counts_matrix=read.csv(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_Count_Matrix/Visvader_0163_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv", header = TRUE, row.names = 1,check.names = FALSE),
                                                  annotations_file=read.table(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_0163_ER_Navin_Visvader_NORM_panCK_Reference.annotations.txt", header = FALSE, row.names = 1),
                                                  delim="\t",
                                                  gene_order_file="/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/hg38_gencode_v27.txt",
                                                  ref_group_names=c("Normal"))


# perform infercnv operations to reveal cnv signal
Visvader_0163_ER_inferCNV = infercnv::run(Visvader_0163_ER_inferCNV,
                                          cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                          out_dir="/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0163_ER_panCK/",  # dir is auto-created for storing outputs
                                          cluster_by_groups=T,   # cluster
                                          denoise=T,
                                          HMM=T,
                                          num_ref_groups = 8, analysis_mode = "subclusters")

###########  Integrate data with Seurat object
# Read RDS and add inferCNV results to metadata
Visvader_0163_ER_Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Visvader_0163_ER_Navin_Visvader_NORM_panCK_Reference.rds") 
Visvader_0163_ER_panCK_cnv <- infercnv::add_to_seurat(infercnv_output_path="/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0163_ER_panCK/",
                                                      seurat_obj=Visvader_0163_ER_Navin_Visvader_NORM_panCK_Reference, # optional
                                                      top_n=10)

# Collect Tumor cells from Seurat Object
Visvader_0163_ER_panCK_cnv <- subset(x = Visvader_0163_ER_panCK_cnv, subset = cnv.ident == "Tumor")
# Quantify chromosomes with alterations per cell
Visvader_0163_ER_panCK_cnv_Pos_Chromosomes <- as.data.frame(rowSums(Visvader_0163_ER_panCK_cnv@meta.data[,startsWith(names(Visvader_0163_ER_panCK_cnv@meta.data),"has_cnv_")] > 0 ))
Visvader_0163_ER_panCK_cnv_Neg_Chromosomes <- as.data.frame(rowSums(Visvader_0163_ER_panCK_cnv@meta.data[,startsWith(names(Visvader_0163_ER_panCK_cnv@meta.data),"has_cnv_")] == 0))
Visvader_0163_ER_panCK_cnv_Chromosome_Quant <- as.data.frame(c(Visvader_0163_ER_panCK_cnv_Pos_Chromosomes, Visvader_0163_ER_panCK_cnv_Neg_Chromosomes))
colnames(Visvader_0163_ER_panCK_cnv_Chromosome_Quant) <- c("Chromosomes CNV Pos", "Chromosomes CNV Neg")
write.csv(Visvader_0163_ER_panCK_cnv_Chromosome_Quant, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0163_ER_panCK/Visvader_0163_ER_panCK_cnv_Chromosome_Quant.csv")

# Save CNV Metadata
Visvader_0163_ER_panCK_cnv.meta.data <- as.data.frame(as.matrix(Visvader_0163_ER_panCK_cnv@meta.data))
write.csv(Visvader_0163_ER_panCK_cnv.meta.data, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0163_ER_panCK/Visvader_0163_ER_panCK_cnv.meta.data.csv")

# Subset Cells with inferred CNVs -------------------------------------------------------------------------------------------------------------
Visvader_0163_ER_panCK_cnv <- subset(Visvader_0163_ER_panCK_cnv, cells=rownames(Visvader_0163_ER_panCK_cnv@meta.data)[rowSums(Visvader_0163_ER_panCK_cnv@meta.data[,startsWith(names(Visvader_0163_ER_panCK_cnv@meta.data),"has_cnv_")]) > 0 ])

# Save RDS with inferredCNVs
saveRDS(Visvader_0163_ER_panCK_cnv, file = "/R/R_Visvader/Visvader_RDS/RDS_inferCNV/Visvader_0163_ER_panCK_cnv.rds")


# Free memory
rm(Visvader_0163_ER_inferCNV)
rm(Visvader_0163_ER_Navin_Visvader_NORM_panCK_Reference)
rm(Visvader_0163_ER_panCK_cnv)
rm(Visvader_0163_ER_panCK_cnv_Pos_Chromosomes)
rm(Visvader_0163_ER_panCK_cnv_Neg_Chromosomes)
rm(Visvader_0163_ER_panCK_cnv_Chromosome_Quant)
rm(Visvader_0163_ER_panCK_cnv.meta.data)
gc()



##########################################
#### InferCNV: Visvader_0167_ER_panCK ####
##########################################

########## Merge with normal mammary glad, extract epithelium, and process via Seurat
# Load RDS
Visvader_0167_ER_panCK <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0167_ER_panCK.rds")
Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Navin_Visvader_NORM_panCK_Reference.rds")

# Set identities 
colnames(x = Visvader_0167_ER_panCK[[]])
Idents(object = Visvader_0167_ER_panCK) <- 'Tumor'
Visvader_0167_ER_panCK[["cnv.ident"]] <- Idents(object = Visvader_0167_ER_panCK)
levels(x = Visvader_0167_ER_panCK)

colnames(x = Navin_Visvader_NORM_panCK_Reference[[]])
Idents(object = Navin_Visvader_NORM_panCK_Reference) <- 'Normal'
Navin_Visvader_NORM_panCK_Reference[["cnv.ident"]] <- Idents(object = Navin_Visvader_NORM_panCK_Reference)
levels(x = Navin_Visvader_NORM_panCK_Reference)

# Merge subsets into recombined RDS
Visvader_0167_ER_Navin_Visvader_NORM_panCK_Reference <- merge(x = Visvader_0167_ER_panCK, y = Navin_Visvader_NORM_panCK_Reference)
# Join Layers
Visvader_0167_ER_Navin_Visvader_NORM_panCK_Reference <- JoinLayers(Visvader_0167_ER_Navin_Visvader_NORM_panCK_Reference)

# Remove subsetted objects
rm(Visvader_0167_ER_panCK)
rm(Navin_Visvader_NORM_panCK_Reference)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Visvader_0167_ER_Navin_Visvader_NORM_panCK_Reference <- FindVariableFeatures(Visvader_0167_ER_Navin_Visvader_NORM_panCK_Reference, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Visvader_0167_ER_Navin_Visvader_NORM_panCK_Reference)
Visvader_0167_ER_Navin_Visvader_NORM_panCK_Reference <- ScaleData(Visvader_0167_ER_Navin_Visvader_NORM_panCK_Reference, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Visvader_0167_ER_Navin_Visvader_NORM_panCK_Reference <- RunPCA(Visvader_0167_ER_Navin_Visvader_NORM_panCK_Reference, features = VariableFeatures(object = Visvader_0167_ER_Navin_Visvader_NORM_panCK_Reference))

# Examine and visualize PCA results a few different ways
print(Visvader_0167_ER_Navin_Visvader_NORM_panCK_Reference[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Visvader_0167_ER_Navin_Visvader_NORM_panCK_Reference, dims = 1:2, reduction = "pca")
DimPlot(Visvader_0167_ER_Navin_Visvader_NORM_panCK_Reference, reduction = "pca")

Visvader_0167_ER_Navin_Visvader_NORM_panCK_Reference <- FindNeighbors(Visvader_0167_ER_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
Visvader_0167_ER_Navin_Visvader_NORM_panCK_Reference <- FindClusters(Visvader_0167_ER_Navin_Visvader_NORM_panCK_Reference, resolution = 0.1)

# Umap clustering
Visvader_0167_ER_Navin_Visvader_NORM_panCK_Reference <- RunUMAP(Visvader_0167_ER_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
DimPlot(Visvader_0167_ER_Navin_Visvader_NORM_panCK_Reference, reduction = "umap")
DimPlot(Visvader_0167_ER_Navin_Visvader_NORM_panCK_Reference, reduction = "umap", group.by = "orig.ident")

# Save RDS
saveRDS(Visvader_0167_ER_Navin_Visvader_NORM_panCK_Reference, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Visvader_0167_ER_Navin_Visvader_NORM_panCK_Reference.rds")

########## Prepare Cell Annotations
# Extract & Save Meta Data
Visvader_0167_ER_Navin_Visvader_NORM_panCK_Reference.meta.data <- as.data.frame(as.matrix(Visvader_0167_ER_Navin_Visvader_NORM_panCK_Reference@meta.data))
# Extract pertinent columns
# Make row names first column
Visvader_0167_ER_Navin_Visvader_NORM_panCK_Reference.meta.data <- tibble::rownames_to_column(Visvader_0167_ER_Navin_Visvader_NORM_panCK_Reference.meta.data, "Cells")
# In Seurat, "X" becomes the cell names while "cnv.ident" is specified above
Visvader_0167_ER_Navin_Visvader_NORM_panCK_Reference.annotations <- Visvader_0167_ER_Navin_Visvader_NORM_panCK_Reference.meta.data[, c("Cells", "cnv.ident")]
# Delete header
names(Visvader_0167_ER_Navin_Visvader_NORM_panCK_Reference.annotations) <- NULL
# Save as table
write.table(Visvader_0167_ER_Navin_Visvader_NORM_panCK_Reference.annotations,"/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_0167_ER_Navin_Visvader_NORM_panCK_Reference.annotations.txt",sep="\t",row.names=FALSE)
# Remove metadata and annotations
rm(Visvader_0167_ER_Navin_Visvader_NORM_panCK_Reference.meta.data)
rm(Visvader_0167_ER_Navin_Visvader_NORM_panCK_Reference.annotations)

########## Generate Count Matrix
# Extract counts from the Seurat object
Visvader_0167_ER_Navin_Visvader_NORM_panCK_Reference.count <- Visvader_0167_ER_Navin_Visvader_NORM_panCK_Reference[["RNA"]]$counts
# Convert counts to a sparse matrix
Visvader_0167_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as(Visvader_0167_ER_Navin_Visvader_NORM_panCK_Reference.count, 'sparseMatrix')
rm(Visvader_0167_ER_Navin_Visvader_NORM_panCK_Reference.count)

# Convert to requisite data format
Visvader_0167_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as.data.frame(Visvader_0167_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix)
# Save table
write.csv(Visvader_0167_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_Count_Matrix/Visvader_0167_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv")
# Remove Objects
rm(Visvader_0167_ER_Navin_Visvader_NORM_panCK_Reference)
rm(Visvader_0167_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix)

########### Run inferCNV
# Note: check.names = FALSE prevents cell names from having dashes replaced with dots; essential to match annotations file
Visvader_0167_ER_inferCNV <- CreateInfercnvObject(raw_counts_matrix=read.csv(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_Count_Matrix/Visvader_0167_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv", header = TRUE, row.names = 1,check.names = FALSE),
                                                  annotations_file=read.table(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_0167_ER_Navin_Visvader_NORM_panCK_Reference.annotations.txt", header = FALSE, row.names = 1),
                                                  delim="\t",
                                                  gene_order_file="/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/hg38_gencode_v27.txt",
                                                  ref_group_names=c("Normal"))


# perform infercnv operations to reveal cnv signal
Visvader_0167_ER_inferCNV = infercnv::run(Visvader_0167_ER_inferCNV,
                                          cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                          out_dir="/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0167_ER_panCK/",  # dir is auto-created for storing outputs
                                          cluster_by_groups=T,   # cluster
                                          denoise=T,
                                          HMM=T,
                                          num_ref_groups = 8, analysis_mode = "subclusters")

###########  Integrate data with Seurat object
# Read RDS and add inferCNV results to metadata
Visvader_0167_ER_Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Visvader_0167_ER_Navin_Visvader_NORM_panCK_Reference.rds") 
Visvader_0167_ER_panCK_cnv <- infercnv::add_to_seurat(infercnv_output_path="/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0167_ER_panCK/",
                                                      seurat_obj=Visvader_0167_ER_Navin_Visvader_NORM_panCK_Reference, # optional
                                                      top_n=10)

# Collect Tumor cells from Seurat Object
Visvader_0167_ER_panCK_cnv <- subset(x = Visvader_0167_ER_panCK_cnv, subset = cnv.ident == "Tumor")
# Quantify chromosomes with alterations per cell
Visvader_0167_ER_panCK_cnv_Pos_Chromosomes <- as.data.frame(rowSums(Visvader_0167_ER_panCK_cnv@meta.data[,startsWith(names(Visvader_0167_ER_panCK_cnv@meta.data),"has_cnv_")] > 0 ))
Visvader_0167_ER_panCK_cnv_Neg_Chromosomes <- as.data.frame(rowSums(Visvader_0167_ER_panCK_cnv@meta.data[,startsWith(names(Visvader_0167_ER_panCK_cnv@meta.data),"has_cnv_")] == 0))
Visvader_0167_ER_panCK_cnv_Chromosome_Quant <- as.data.frame(c(Visvader_0167_ER_panCK_cnv_Pos_Chromosomes, Visvader_0167_ER_panCK_cnv_Neg_Chromosomes))
colnames(Visvader_0167_ER_panCK_cnv_Chromosome_Quant) <- c("Chromosomes CNV Pos", "Chromosomes CNV Neg")
write.csv(Visvader_0167_ER_panCK_cnv_Chromosome_Quant, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0167_ER_panCK/Visvader_0167_ER_panCK_cnv_Chromosome_Quant.csv")

# Save CNV Metadata
Visvader_0167_ER_panCK_cnv.meta.data <- as.data.frame(as.matrix(Visvader_0167_ER_panCK_cnv@meta.data))
write.csv(Visvader_0167_ER_panCK_cnv.meta.data, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0167_ER_panCK/Visvader_0167_ER_panCK_cnv.meta.data.csv")

# Subset Cells with inferred CNVs -------------------------------------------------------------------------------------------------------------
Visvader_0167_ER_panCK_cnv <- subset(Visvader_0167_ER_panCK_cnv, cells=rownames(Visvader_0167_ER_panCK_cnv@meta.data)[rowSums(Visvader_0167_ER_panCK_cnv@meta.data[,startsWith(names(Visvader_0167_ER_panCK_cnv@meta.data),"has_cnv_")]) > 0 ])

# Save RDS with inferredCNVs
saveRDS(Visvader_0167_ER_panCK_cnv, file = "/R/R_Visvader/Visvader_RDS/RDS_inferCNV/Visvader_0167_ER_panCK_cnv.rds")


# Free memory
rm(Visvader_0167_ER_inferCNV)
rm(Visvader_0167_ER_Navin_Visvader_NORM_panCK_Reference)
rm(Visvader_0167_ER_panCK_cnv)
rm(Visvader_0167_ER_panCK_cnv_Pos_Chromosomes)
rm(Visvader_0167_ER_panCK_cnv_Neg_Chromosomes)
rm(Visvader_0167_ER_panCK_cnv_Chromosome_Quant)
rm(Visvader_0167_ER_panCK_cnv.meta.data)
gc()



##########################################
#### InferCNV: Visvader_0173_ER_panCK ####
##########################################

########## Merge with normal mammary glad, extract epithelium, and process via Seurat
# Load RDS
Visvader_0173_ER_panCK <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0173_ER_panCK.rds")
Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Navin_Visvader_NORM_panCK_Reference.rds")

# Set identities 
colnames(x = Visvader_0173_ER_panCK[[]])
Idents(object = Visvader_0173_ER_panCK) <- 'Tumor'
Visvader_0173_ER_panCK[["cnv.ident"]] <- Idents(object = Visvader_0173_ER_panCK)
levels(x = Visvader_0173_ER_panCK)

colnames(x = Navin_Visvader_NORM_panCK_Reference[[]])
Idents(object = Navin_Visvader_NORM_panCK_Reference) <- 'Normal'
Navin_Visvader_NORM_panCK_Reference[["cnv.ident"]] <- Idents(object = Navin_Visvader_NORM_panCK_Reference)
levels(x = Navin_Visvader_NORM_panCK_Reference)

# Merge subsets into recombined RDS
Visvader_0173_ER_Navin_Visvader_NORM_panCK_Reference <- merge(x = Visvader_0173_ER_panCK, y = Navin_Visvader_NORM_panCK_Reference)
# Join Layers
Visvader_0173_ER_Navin_Visvader_NORM_panCK_Reference <- JoinLayers(Visvader_0173_ER_Navin_Visvader_NORM_panCK_Reference)

# Remove subsetted objects
rm(Visvader_0173_ER_panCK)
rm(Navin_Visvader_NORM_panCK_Reference)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Visvader_0173_ER_Navin_Visvader_NORM_panCK_Reference <- FindVariableFeatures(Visvader_0173_ER_Navin_Visvader_NORM_panCK_Reference, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Visvader_0173_ER_Navin_Visvader_NORM_panCK_Reference)
Visvader_0173_ER_Navin_Visvader_NORM_panCK_Reference <- ScaleData(Visvader_0173_ER_Navin_Visvader_NORM_panCK_Reference, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Visvader_0173_ER_Navin_Visvader_NORM_panCK_Reference <- RunPCA(Visvader_0173_ER_Navin_Visvader_NORM_panCK_Reference, features = VariableFeatures(object = Visvader_0173_ER_Navin_Visvader_NORM_panCK_Reference))

# Examine and visualize PCA results a few different ways
print(Visvader_0173_ER_Navin_Visvader_NORM_panCK_Reference[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Visvader_0173_ER_Navin_Visvader_NORM_panCK_Reference, dims = 1:2, reduction = "pca")
DimPlot(Visvader_0173_ER_Navin_Visvader_NORM_panCK_Reference, reduction = "pca")

Visvader_0173_ER_Navin_Visvader_NORM_panCK_Reference <- FindNeighbors(Visvader_0173_ER_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
Visvader_0173_ER_Navin_Visvader_NORM_panCK_Reference <- FindClusters(Visvader_0173_ER_Navin_Visvader_NORM_panCK_Reference, resolution = 0.1)

# Umap clustering
Visvader_0173_ER_Navin_Visvader_NORM_panCK_Reference <- RunUMAP(Visvader_0173_ER_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
DimPlot(Visvader_0173_ER_Navin_Visvader_NORM_panCK_Reference, reduction = "umap")
DimPlot(Visvader_0173_ER_Navin_Visvader_NORM_panCK_Reference, reduction = "umap", group.by = "orig.ident")

# Save RDS
saveRDS(Visvader_0173_ER_Navin_Visvader_NORM_panCK_Reference, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Visvader_0173_ER_Navin_Visvader_NORM_panCK_Reference.rds")

########## Prepare Cell Annotations
# Extract & Save Meta Data
Visvader_0173_ER_Navin_Visvader_NORM_panCK_Reference.meta.data <- as.data.frame(as.matrix(Visvader_0173_ER_Navin_Visvader_NORM_panCK_Reference@meta.data))
# Extract pertinent columns
# Make row names first column
Visvader_0173_ER_Navin_Visvader_NORM_panCK_Reference.meta.data <- tibble::rownames_to_column(Visvader_0173_ER_Navin_Visvader_NORM_panCK_Reference.meta.data, "Cells")
# In Seurat, "X" becomes the cell names while "cnv.ident" is specified above
Visvader_0173_ER_Navin_Visvader_NORM_panCK_Reference.annotations <- Visvader_0173_ER_Navin_Visvader_NORM_panCK_Reference.meta.data[, c("Cells", "cnv.ident")]
# Delete header
names(Visvader_0173_ER_Navin_Visvader_NORM_panCK_Reference.annotations) <- NULL
# Save as table
write.table(Visvader_0173_ER_Navin_Visvader_NORM_panCK_Reference.annotations,"/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_0173_ER_Navin_Visvader_NORM_panCK_Reference.annotations.txt",sep="\t",row.names=FALSE)
# Remove metadata and annotations
rm(Visvader_0173_ER_Navin_Visvader_NORM_panCK_Reference.meta.data)
rm(Visvader_0173_ER_Navin_Visvader_NORM_panCK_Reference.annotations)

########## Generate Count Matrix
# Extract counts from the Seurat object
Visvader_0173_ER_Navin_Visvader_NORM_panCK_Reference.count <- Visvader_0173_ER_Navin_Visvader_NORM_panCK_Reference[["RNA"]]$counts
# Convert counts to a sparse matrix
Visvader_0173_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as(Visvader_0173_ER_Navin_Visvader_NORM_panCK_Reference.count, 'sparseMatrix')
rm(Visvader_0173_ER_Navin_Visvader_NORM_panCK_Reference.count)

# Convert to requisite data format
Visvader_0173_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as.data.frame(Visvader_0173_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix)
# Save table
write.csv(Visvader_0173_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_Count_Matrix/Visvader_0173_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv")
# Remove Objects
rm(Visvader_0173_ER_Navin_Visvader_NORM_panCK_Reference)
rm(Visvader_0173_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix)

########### Run inferCNV
# Note: check.names = FALSE prevents cell names from having dashes replaced with dots; essential to match annotations file
Visvader_0173_ER_inferCNV <- CreateInfercnvObject(raw_counts_matrix=read.csv(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_Count_Matrix/Visvader_0173_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv", header = TRUE, row.names = 1,check.names = FALSE),
                                                  annotations_file=read.table(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_0173_ER_Navin_Visvader_NORM_panCK_Reference.annotations.txt", header = FALSE, row.names = 1),
                                                  delim="\t",
                                                  gene_order_file="/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/hg38_gencode_v27.txt",
                                                  ref_group_names=c("Normal"))


# perform infercnv operations to reveal cnv signal
Visvader_0173_ER_inferCNV = infercnv::run(Visvader_0173_ER_inferCNV,
                                          cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                          out_dir="/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0173_ER_panCK/",  # dir is auto-created for storing outputs
                                          cluster_by_groups=T,   # cluster
                                          denoise=T,
                                          HMM=T,
                                          num_ref_groups = 8, analysis_mode = "subclusters")

###########  Integrate data with Seurat object
# Read RDS and add inferCNV results to metadata
Visvader_0173_ER_Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Visvader_0173_ER_Navin_Visvader_NORM_panCK_Reference.rds") 
Visvader_0173_ER_panCK_cnv <- infercnv::add_to_seurat(infercnv_output_path="/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0173_ER_panCK/",
                                                      seurat_obj=Visvader_0173_ER_Navin_Visvader_NORM_panCK_Reference, # optional
                                                      top_n=10)

# Collect Tumor cells from Seurat Object
Visvader_0173_ER_panCK_cnv <- subset(x = Visvader_0173_ER_panCK_cnv, subset = cnv.ident == "Tumor")
# Quantify chromosomes with alterations per cell
Visvader_0173_ER_panCK_cnv_Pos_Chromosomes <- as.data.frame(rowSums(Visvader_0173_ER_panCK_cnv@meta.data[,startsWith(names(Visvader_0173_ER_panCK_cnv@meta.data),"has_cnv_")] > 0 ))
Visvader_0173_ER_panCK_cnv_Neg_Chromosomes <- as.data.frame(rowSums(Visvader_0173_ER_panCK_cnv@meta.data[,startsWith(names(Visvader_0173_ER_panCK_cnv@meta.data),"has_cnv_")] == 0))
Visvader_0173_ER_panCK_cnv_Chromosome_Quant <- as.data.frame(c(Visvader_0173_ER_panCK_cnv_Pos_Chromosomes, Visvader_0173_ER_panCK_cnv_Neg_Chromosomes))
colnames(Visvader_0173_ER_panCK_cnv_Chromosome_Quant) <- c("Chromosomes CNV Pos", "Chromosomes CNV Neg")
write.csv(Visvader_0173_ER_panCK_cnv_Chromosome_Quant, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0173_ER_panCK/Visvader_0173_ER_panCK_cnv_Chromosome_Quant.csv")

# Save CNV Metadata
Visvader_0173_ER_panCK_cnv.meta.data <- as.data.frame(as.matrix(Visvader_0173_ER_panCK_cnv@meta.data))
write.csv(Visvader_0173_ER_panCK_cnv.meta.data, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0173_ER_panCK/Visvader_0173_ER_panCK_cnv.meta.data.csv")

# Subset Cells with inferred CNVs -------------------------------------------------------------------------------------------------------------
Visvader_0173_ER_panCK_cnv <- subset(Visvader_0173_ER_panCK_cnv, cells=rownames(Visvader_0173_ER_panCK_cnv@meta.data)[rowSums(Visvader_0173_ER_panCK_cnv@meta.data[,startsWith(names(Visvader_0173_ER_panCK_cnv@meta.data),"has_cnv_")]) > 0 ])

# Save RDS with inferredCNVs
saveRDS(Visvader_0173_ER_panCK_cnv, file = "/R/R_Visvader/Visvader_RDS/RDS_inferCNV/Visvader_0173_ER_panCK_cnv.rds")


# Free memory
rm(Visvader_0173_ER_inferCNV)
rm(Visvader_0173_ER_Navin_Visvader_NORM_panCK_Reference)
rm(Visvader_0173_ER_panCK_cnv)
rm(Visvader_0173_ER_panCK_cnv_Pos_Chromosomes)
rm(Visvader_0173_ER_panCK_cnv_Neg_Chromosomes)
rm(Visvader_0173_ER_panCK_cnv_Chromosome_Quant)
rm(Visvader_0173_ER_panCK_cnv.meta.data)
gc()



##########################################
#### InferCNV: Visvader_0178_ER_panCK ####
##########################################

########## Merge with normal mammary glad, extract epithelium, and process via Seurat
# Load RDS
Visvader_0178_ER_panCK <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0178_ER_panCK.rds")
Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Navin_Visvader_NORM_panCK_Reference.rds")

# Set identities 
colnames(x = Visvader_0178_ER_panCK[[]])
Idents(object = Visvader_0178_ER_panCK) <- 'Tumor'
Visvader_0178_ER_panCK[["cnv.ident"]] <- Idents(object = Visvader_0178_ER_panCK)
levels(x = Visvader_0178_ER_panCK)

colnames(x = Navin_Visvader_NORM_panCK_Reference[[]])
Idents(object = Navin_Visvader_NORM_panCK_Reference) <- 'Normal'
Navin_Visvader_NORM_panCK_Reference[["cnv.ident"]] <- Idents(object = Navin_Visvader_NORM_panCK_Reference)
levels(x = Navin_Visvader_NORM_panCK_Reference)

# Merge subsets into recombined RDS
Visvader_0178_ER_Navin_Visvader_NORM_panCK_Reference <- merge(x = Visvader_0178_ER_panCK, y = Navin_Visvader_NORM_panCK_Reference)
# Join Layers
Visvader_0178_ER_Navin_Visvader_NORM_panCK_Reference <- JoinLayers(Visvader_0178_ER_Navin_Visvader_NORM_panCK_Reference)

# Remove subsetted objects
rm(Visvader_0178_ER_panCK)
rm(Navin_Visvader_NORM_panCK_Reference)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Visvader_0178_ER_Navin_Visvader_NORM_panCK_Reference <- FindVariableFeatures(Visvader_0178_ER_Navin_Visvader_NORM_panCK_Reference, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Visvader_0178_ER_Navin_Visvader_NORM_panCK_Reference)
Visvader_0178_ER_Navin_Visvader_NORM_panCK_Reference <- ScaleData(Visvader_0178_ER_Navin_Visvader_NORM_panCK_Reference, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Visvader_0178_ER_Navin_Visvader_NORM_panCK_Reference <- RunPCA(Visvader_0178_ER_Navin_Visvader_NORM_panCK_Reference, features = VariableFeatures(object = Visvader_0178_ER_Navin_Visvader_NORM_panCK_Reference))

# Examine and visualize PCA results a few different ways
print(Visvader_0178_ER_Navin_Visvader_NORM_panCK_Reference[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Visvader_0178_ER_Navin_Visvader_NORM_panCK_Reference, dims = 1:2, reduction = "pca")
DimPlot(Visvader_0178_ER_Navin_Visvader_NORM_panCK_Reference, reduction = "pca")

Visvader_0178_ER_Navin_Visvader_NORM_panCK_Reference <- FindNeighbors(Visvader_0178_ER_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
Visvader_0178_ER_Navin_Visvader_NORM_panCK_Reference <- FindClusters(Visvader_0178_ER_Navin_Visvader_NORM_panCK_Reference, resolution = 0.1)

# Umap clustering
Visvader_0178_ER_Navin_Visvader_NORM_panCK_Reference <- RunUMAP(Visvader_0178_ER_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
DimPlot(Visvader_0178_ER_Navin_Visvader_NORM_panCK_Reference, reduction = "umap")
DimPlot(Visvader_0178_ER_Navin_Visvader_NORM_panCK_Reference, reduction = "umap", group.by = "orig.ident")

# Save RDS
saveRDS(Visvader_0178_ER_Navin_Visvader_NORM_panCK_Reference, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Visvader_0178_ER_Navin_Visvader_NORM_panCK_Reference.rds")

########## Prepare Cell Annotations
# Extract & Save Meta Data
Visvader_0178_ER_Navin_Visvader_NORM_panCK_Reference.meta.data <- as.data.frame(as.matrix(Visvader_0178_ER_Navin_Visvader_NORM_panCK_Reference@meta.data))
# Extract pertinent columns
# Make row names first column
Visvader_0178_ER_Navin_Visvader_NORM_panCK_Reference.meta.data <- tibble::rownames_to_column(Visvader_0178_ER_Navin_Visvader_NORM_panCK_Reference.meta.data, "Cells")
# In Seurat, "X" becomes the cell names while "cnv.ident" is specified above
Visvader_0178_ER_Navin_Visvader_NORM_panCK_Reference.annotations <- Visvader_0178_ER_Navin_Visvader_NORM_panCK_Reference.meta.data[, c("Cells", "cnv.ident")]
# Delete header
names(Visvader_0178_ER_Navin_Visvader_NORM_panCK_Reference.annotations) <- NULL
# Save as table
write.table(Visvader_0178_ER_Navin_Visvader_NORM_panCK_Reference.annotations,"/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_0178_ER_Navin_Visvader_NORM_panCK_Reference.annotations.txt",sep="\t",row.names=FALSE)
# Remove metadata and annotations
rm(Visvader_0178_ER_Navin_Visvader_NORM_panCK_Reference.meta.data)
rm(Visvader_0178_ER_Navin_Visvader_NORM_panCK_Reference.annotations)

########## Generate Count Matrix
# Extract counts from the Seurat object
Visvader_0178_ER_Navin_Visvader_NORM_panCK_Reference.count <- Visvader_0178_ER_Navin_Visvader_NORM_panCK_Reference[["RNA"]]$counts
# Convert counts to a sparse matrix
Visvader_0178_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as(Visvader_0178_ER_Navin_Visvader_NORM_panCK_Reference.count, 'sparseMatrix')
rm(Visvader_0178_ER_Navin_Visvader_NORM_panCK_Reference.count)

# Convert to requisite data format
Visvader_0178_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as.data.frame(Visvader_0178_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix)
# Save table
write.csv(Visvader_0178_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_Count_Matrix/Visvader_0178_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv")
# Remove Objects
rm(Visvader_0178_ER_Navin_Visvader_NORM_panCK_Reference)
rm(Visvader_0178_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix)

########### Run inferCNV
# Note: check.names = FALSE prevents cell names from having dashes replaced with dots; essential to match annotations file
Visvader_0178_ER_inferCNV <- CreateInfercnvObject(raw_counts_matrix=read.csv(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_Count_Matrix/Visvader_0178_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv", header = TRUE, row.names = 1,check.names = FALSE),
                                                  annotations_file=read.table(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_0178_ER_Navin_Visvader_NORM_panCK_Reference.annotations.txt", header = FALSE, row.names = 1),
                                                  delim="\t",
                                                  gene_order_file="/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/hg38_gencode_v27.txt",
                                                  ref_group_names=c("Normal"))


# perform infercnv operations to reveal cnv signal
Visvader_0178_ER_inferCNV = infercnv::run(Visvader_0178_ER_inferCNV,
                                          cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                          out_dir="/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0178_ER_panCK/",  # dir is auto-created for storing outputs
                                          cluster_by_groups=T,   # cluster
                                          denoise=T,
                                          HMM=T,
                                          num_ref_groups = 8, analysis_mode = "subclusters")

###########  Integrate data with Seurat object
# Read RDS and add inferCNV results to metadata
Visvader_0178_ER_Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Visvader_0178_ER_Navin_Visvader_NORM_panCK_Reference.rds") 
Visvader_0178_ER_panCK_cnv <- infercnv::add_to_seurat(infercnv_output_path="/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0178_ER_panCK/",
                                                      seurat_obj=Visvader_0178_ER_Navin_Visvader_NORM_panCK_Reference, # optional
                                                      top_n=10)

# Collect Tumor cells from Seurat Object
Visvader_0178_ER_panCK_cnv <- subset(x = Visvader_0178_ER_panCK_cnv, subset = cnv.ident == "Tumor")
# Quantify chromosomes with alterations per cell
Visvader_0178_ER_panCK_cnv_Pos_Chromosomes <- as.data.frame(rowSums(Visvader_0178_ER_panCK_cnv@meta.data[,startsWith(names(Visvader_0178_ER_panCK_cnv@meta.data),"has_cnv_")] > 0 ))
Visvader_0178_ER_panCK_cnv_Neg_Chromosomes <- as.data.frame(rowSums(Visvader_0178_ER_panCK_cnv@meta.data[,startsWith(names(Visvader_0178_ER_panCK_cnv@meta.data),"has_cnv_")] == 0))
Visvader_0178_ER_panCK_cnv_Chromosome_Quant <- as.data.frame(c(Visvader_0178_ER_panCK_cnv_Pos_Chromosomes, Visvader_0178_ER_panCK_cnv_Neg_Chromosomes))
colnames(Visvader_0178_ER_panCK_cnv_Chromosome_Quant) <- c("Chromosomes CNV Pos", "Chromosomes CNV Neg")
write.csv(Visvader_0178_ER_panCK_cnv_Chromosome_Quant, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0178_ER_panCK/Visvader_0178_ER_panCK_cnv_Chromosome_Quant.csv")

# Save CNV Metadata
Visvader_0178_ER_panCK_cnv.meta.data <- as.data.frame(as.matrix(Visvader_0178_ER_panCK_cnv@meta.data))
write.csv(Visvader_0178_ER_panCK_cnv.meta.data, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0178_ER_panCK/Visvader_0178_ER_panCK_cnv.meta.data.csv")

# Subset Cells with inferred CNVs -------------------------------------------------------------------------------------------------------------
Visvader_0178_ER_panCK_cnv <- subset(Visvader_0178_ER_panCK_cnv, cells=rownames(Visvader_0178_ER_panCK_cnv@meta.data)[rowSums(Visvader_0178_ER_panCK_cnv@meta.data[,startsWith(names(Visvader_0178_ER_panCK_cnv@meta.data),"has_cnv_")]) > 0 ])

# Save RDS with inferredCNVs
saveRDS(Visvader_0178_ER_panCK_cnv, file = "/R/R_Visvader/Visvader_RDS/RDS_inferCNV/Visvader_0178_ER_panCK_cnv.rds")


# Free memory
rm(Visvader_0178_ER_inferCNV)
rm(Visvader_0178_ER_Navin_Visvader_NORM_panCK_Reference)
rm(Visvader_0178_ER_panCK_cnv)
rm(Visvader_0178_ER_panCK_cnv_Pos_Chromosomes)
rm(Visvader_0178_ER_panCK_cnv_Neg_Chromosomes)
rm(Visvader_0178_ER_panCK_cnv_Chromosome_Quant)
rm(Visvader_0178_ER_panCK_cnv.meta.data)
gc()



##########################################
#### InferCNV: Visvader_0360_ER_panCK ####
##########################################

########## Merge with normal mammary glad, extract epithelium, and process via Seurat
# Load RDS
Visvader_0360_ER_panCK <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0360_ER_panCK.rds")
Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Navin_Visvader_NORM_panCK_Reference.rds")

# Set identities 
colnames(x = Visvader_0360_ER_panCK[[]])
Idents(object = Visvader_0360_ER_panCK) <- 'Tumor'
Visvader_0360_ER_panCK[["cnv.ident"]] <- Idents(object = Visvader_0360_ER_panCK)
levels(x = Visvader_0360_ER_panCK)

colnames(x = Navin_Visvader_NORM_panCK_Reference[[]])
Idents(object = Navin_Visvader_NORM_panCK_Reference) <- 'Normal'
Navin_Visvader_NORM_panCK_Reference[["cnv.ident"]] <- Idents(object = Navin_Visvader_NORM_panCK_Reference)
levels(x = Navin_Visvader_NORM_panCK_Reference)

# Merge subsets into recombined RDS
Visvader_0360_ER_Navin_Visvader_NORM_panCK_Reference <- merge(x = Visvader_0360_ER_panCK, y = Navin_Visvader_NORM_panCK_Reference)
# Join Layers
Visvader_0360_ER_Navin_Visvader_NORM_panCK_Reference <- JoinLayers(Visvader_0360_ER_Navin_Visvader_NORM_panCK_Reference)

# Remove subsetted objects
rm(Visvader_0360_ER_panCK)
rm(Navin_Visvader_NORM_panCK_Reference)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Visvader_0360_ER_Navin_Visvader_NORM_panCK_Reference <- FindVariableFeatures(Visvader_0360_ER_Navin_Visvader_NORM_panCK_Reference, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Visvader_0360_ER_Navin_Visvader_NORM_panCK_Reference)
Visvader_0360_ER_Navin_Visvader_NORM_panCK_Reference <- ScaleData(Visvader_0360_ER_Navin_Visvader_NORM_panCK_Reference, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Visvader_0360_ER_Navin_Visvader_NORM_panCK_Reference <- RunPCA(Visvader_0360_ER_Navin_Visvader_NORM_panCK_Reference, features = VariableFeatures(object = Visvader_0360_ER_Navin_Visvader_NORM_panCK_Reference))

# Examine and visualize PCA results a few different ways
print(Visvader_0360_ER_Navin_Visvader_NORM_panCK_Reference[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Visvader_0360_ER_Navin_Visvader_NORM_panCK_Reference, dims = 1:2, reduction = "pca")
DimPlot(Visvader_0360_ER_Navin_Visvader_NORM_panCK_Reference, reduction = "pca")

Visvader_0360_ER_Navin_Visvader_NORM_panCK_Reference <- FindNeighbors(Visvader_0360_ER_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
Visvader_0360_ER_Navin_Visvader_NORM_panCK_Reference <- FindClusters(Visvader_0360_ER_Navin_Visvader_NORM_panCK_Reference, resolution = 0.1)

# Umap clustering
Visvader_0360_ER_Navin_Visvader_NORM_panCK_Reference <- RunUMAP(Visvader_0360_ER_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
DimPlot(Visvader_0360_ER_Navin_Visvader_NORM_panCK_Reference, reduction = "umap")
DimPlot(Visvader_0360_ER_Navin_Visvader_NORM_panCK_Reference, reduction = "umap", group.by = "orig.ident")

# Save RDS
saveRDS(Visvader_0360_ER_Navin_Visvader_NORM_panCK_Reference, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Visvader_0360_ER_Navin_Visvader_NORM_panCK_Reference.rds")

########## Prepare Cell Annotations
# Extract & Save Meta Data
Visvader_0360_ER_Navin_Visvader_NORM_panCK_Reference.meta.data <- as.data.frame(as.matrix(Visvader_0360_ER_Navin_Visvader_NORM_panCK_Reference@meta.data))
# Extract pertinent columns
# Make row names first column
Visvader_0360_ER_Navin_Visvader_NORM_panCK_Reference.meta.data <- tibble::rownames_to_column(Visvader_0360_ER_Navin_Visvader_NORM_panCK_Reference.meta.data, "Cells")
# In Seurat, "X" becomes the cell names while "cnv.ident" is specified above
Visvader_0360_ER_Navin_Visvader_NORM_panCK_Reference.annotations <- Visvader_0360_ER_Navin_Visvader_NORM_panCK_Reference.meta.data[, c("Cells", "cnv.ident")]
# Delete header
names(Visvader_0360_ER_Navin_Visvader_NORM_panCK_Reference.annotations) <- NULL
# Save as table
write.table(Visvader_0360_ER_Navin_Visvader_NORM_panCK_Reference.annotations,"/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_0360_ER_Navin_Visvader_NORM_panCK_Reference.annotations.txt",sep="\t",row.names=FALSE)
# Remove metadata and annotations
rm(Visvader_0360_ER_Navin_Visvader_NORM_panCK_Reference.meta.data)
rm(Visvader_0360_ER_Navin_Visvader_NORM_panCK_Reference.annotations)

########## Generate Count Matrix
# Extract counts from the Seurat object
Visvader_0360_ER_Navin_Visvader_NORM_panCK_Reference.count <- Visvader_0360_ER_Navin_Visvader_NORM_panCK_Reference[["RNA"]]$counts
# Convert counts to a sparse matrix
Visvader_0360_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as(Visvader_0360_ER_Navin_Visvader_NORM_panCK_Reference.count, 'sparseMatrix')
rm(Visvader_0360_ER_Navin_Visvader_NORM_panCK_Reference.count)

# Convert to requisite data format
Visvader_0360_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as.data.frame(Visvader_0360_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix)
# Save table
write.csv(Visvader_0360_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_Count_Matrix/Visvader_0360_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv")
# Remove Objects
rm(Visvader_0360_ER_Navin_Visvader_NORM_panCK_Reference)
rm(Visvader_0360_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix)

########### Run inferCNV
# Note: check.names = FALSE prevents cell names from having dashes replaced with dots; essential to match annotations file
Visvader_0360_ER_inferCNV <- CreateInfercnvObject(raw_counts_matrix=read.csv(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_Count_Matrix/Visvader_0360_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv", header = TRUE, row.names = 1,check.names = FALSE),
                                                  annotations_file=read.table(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_0360_ER_Navin_Visvader_NORM_panCK_Reference.annotations.txt", header = FALSE, row.names = 1),
                                                  delim="\t",
                                                  gene_order_file="/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/hg38_gencode_v27.txt",
                                                  ref_group_names=c("Normal"))


# perform infercnv operations to reveal cnv signal
Visvader_0360_ER_inferCNV = infercnv::run(Visvader_0360_ER_inferCNV,
                                          cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                          out_dir="/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0360_ER_panCK/",  # dir is auto-created for storing outputs
                                          cluster_by_groups=T,   # cluster
                                          denoise=T,
                                          HMM=T,
                                          num_ref_groups = 8, analysis_mode = "subclusters")

###########  Integrate data with Seurat object
# Read RDS and add inferCNV results to metadata
Visvader_0360_ER_Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Visvader_0360_ER_Navin_Visvader_NORM_panCK_Reference.rds") 
Visvader_0360_ER_panCK_cnv <- infercnv::add_to_seurat(infercnv_output_path="/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0360_ER_panCK/",
                                                      seurat_obj=Visvader_0360_ER_Navin_Visvader_NORM_panCK_Reference, # optional
                                                      top_n=10)

# Collect Tumor cells from Seurat Object
Visvader_0360_ER_panCK_cnv <- subset(x = Visvader_0360_ER_panCK_cnv, subset = cnv.ident == "Tumor")
# Quantify chromosomes with alterations per cell
Visvader_0360_ER_panCK_cnv_Pos_Chromosomes <- as.data.frame(rowSums(Visvader_0360_ER_panCK_cnv@meta.data[,startsWith(names(Visvader_0360_ER_panCK_cnv@meta.data),"has_cnv_")] > 0 ))
Visvader_0360_ER_panCK_cnv_Neg_Chromosomes <- as.data.frame(rowSums(Visvader_0360_ER_panCK_cnv@meta.data[,startsWith(names(Visvader_0360_ER_panCK_cnv@meta.data),"has_cnv_")] == 0))
Visvader_0360_ER_panCK_cnv_Chromosome_Quant <- as.data.frame(c(Visvader_0360_ER_panCK_cnv_Pos_Chromosomes, Visvader_0360_ER_panCK_cnv_Neg_Chromosomes))
colnames(Visvader_0360_ER_panCK_cnv_Chromosome_Quant) <- c("Chromosomes CNV Pos", "Chromosomes CNV Neg")
write.csv(Visvader_0360_ER_panCK_cnv_Chromosome_Quant, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0360_ER_panCK/Visvader_0360_ER_panCK_cnv_Chromosome_Quant.csv")

# Save CNV Metadata
Visvader_0360_ER_panCK_cnv.meta.data <- as.data.frame(as.matrix(Visvader_0360_ER_panCK_cnv@meta.data))
write.csv(Visvader_0360_ER_panCK_cnv.meta.data, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0360_ER_panCK/Visvader_0360_ER_panCK_cnv.meta.data.csv")

# Subset Cells with inferred CNVs -------------------------------------------------------------------------------------------------------------
Visvader_0360_ER_panCK_cnv <- subset(Visvader_0360_ER_panCK_cnv, cells=rownames(Visvader_0360_ER_panCK_cnv@meta.data)[rowSums(Visvader_0360_ER_panCK_cnv@meta.data[,startsWith(names(Visvader_0360_ER_panCK_cnv@meta.data),"has_cnv_")]) > 0 ])

# Save RDS with inferredCNVs
saveRDS(Visvader_0360_ER_panCK_cnv, file = "/R/R_Visvader/Visvader_RDS/RDS_inferCNV/Visvader_0360_ER_panCK_cnv.rds")


# Free memory
rm(Visvader_0360_ER_inferCNV)
rm(Visvader_0360_ER_Navin_Visvader_NORM_panCK_Reference)
rm(Visvader_0360_ER_panCK_cnv)
rm(Visvader_0360_ER_panCK_cnv_Pos_Chromosomes)
rm(Visvader_0360_ER_panCK_cnv_Neg_Chromosomes)
rm(Visvader_0360_ER_panCK_cnv_Chromosome_Quant)
rm(Visvader_0360_ER_panCK_cnv.meta.data)
gc()



##########################################
#### InferCNV: Visvader_0319_PR_panCK ####
##########################################

########## Merge with normal mammary glad, extract epithelium, and process via Seurat
# Load RDS
Visvader_0319_PR_panCK <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0319_PR_panCK.rds")
Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Navin_Visvader_NORM_panCK_Reference.rds")

# Set identities 
colnames(x = Visvader_0319_PR_panCK[[]])
Idents(object = Visvader_0319_PR_panCK) <- 'Tumor'
Visvader_0319_PR_panCK[["cnv.ident"]] <- Idents(object = Visvader_0319_PR_panCK)
levels(x = Visvader_0319_PR_panCK)

colnames(x = Navin_Visvader_NORM_panCK_Reference[[]])
Idents(object = Navin_Visvader_NORM_panCK_Reference) <- 'Normal'
Navin_Visvader_NORM_panCK_Reference[["cnv.ident"]] <- Idents(object = Navin_Visvader_NORM_panCK_Reference)
levels(x = Navin_Visvader_NORM_panCK_Reference)

# Merge subsets into recombined RDS
Visvader_0319_PR_Navin_Visvader_NORM_panCK_Reference <- merge(x = Visvader_0319_PR_panCK, y = Navin_Visvader_NORM_panCK_Reference)
# Join Layers
Visvader_0319_PR_Navin_Visvader_NORM_panCK_Reference <- JoinLayers(Visvader_0319_PR_Navin_Visvader_NORM_panCK_Reference)

# Remove subsetted objects
rm(Visvader_0319_PR_panCK)
rm(Navin_Visvader_NORM_panCK_Reference)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Visvader_0319_PR_Navin_Visvader_NORM_panCK_Reference <- FindVariableFeatures(Visvader_0319_PR_Navin_Visvader_NORM_panCK_Reference, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Visvader_0319_PR_Navin_Visvader_NORM_panCK_Reference)
Visvader_0319_PR_Navin_Visvader_NORM_panCK_Reference <- ScaleData(Visvader_0319_PR_Navin_Visvader_NORM_panCK_Reference, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Visvader_0319_PR_Navin_Visvader_NORM_panCK_Reference <- RunPCA(Visvader_0319_PR_Navin_Visvader_NORM_panCK_Reference, features = VariableFeatures(object = Visvader_0319_PR_Navin_Visvader_NORM_panCK_Reference))

# Examine and visualize PCA results a few different ways
print(Visvader_0319_PR_Navin_Visvader_NORM_panCK_Reference[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Visvader_0319_PR_Navin_Visvader_NORM_panCK_Reference, dims = 1:2, reduction = "pca")
DimPlot(Visvader_0319_PR_Navin_Visvader_NORM_panCK_Reference, reduction = "pca")

Visvader_0319_PR_Navin_Visvader_NORM_panCK_Reference <- FindNeighbors(Visvader_0319_PR_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
Visvader_0319_PR_Navin_Visvader_NORM_panCK_Reference <- FindClusters(Visvader_0319_PR_Navin_Visvader_NORM_panCK_Reference, resolution = 0.1)

# Umap clustering
Visvader_0319_PR_Navin_Visvader_NORM_panCK_Reference <- RunUMAP(Visvader_0319_PR_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
DimPlot(Visvader_0319_PR_Navin_Visvader_NORM_panCK_Reference, reduction = "umap")
DimPlot(Visvader_0319_PR_Navin_Visvader_NORM_panCK_Reference, reduction = "umap", group.by = "orig.ident")

# Save RDS
saveRDS(Visvader_0319_PR_Navin_Visvader_NORM_panCK_Reference, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Visvader_0319_PR_Navin_Visvader_NORM_panCK_Reference.rds")

########## Prepare Cell Annotations
# Extract & Save Meta Data
Visvader_0319_PR_Navin_Visvader_NORM_panCK_Reference.meta.data <- as.data.frame(as.matrix(Visvader_0319_PR_Navin_Visvader_NORM_panCK_Reference@meta.data))
# Extract pertinent columns
# Make row names first column
Visvader_0319_PR_Navin_Visvader_NORM_panCK_Reference.meta.data <- tibble::rownames_to_column(Visvader_0319_PR_Navin_Visvader_NORM_panCK_Reference.meta.data, "Cells")
# In Seurat, "X" becomes the cell names while "cnv.ident" is specified above
Visvader_0319_PR_Navin_Visvader_NORM_panCK_Reference.annotations <- Visvader_0319_PR_Navin_Visvader_NORM_panCK_Reference.meta.data[, c("Cells", "cnv.ident")]
# Delete header
names(Visvader_0319_PR_Navin_Visvader_NORM_panCK_Reference.annotations) <- NULL
# Save as table
write.table(Visvader_0319_PR_Navin_Visvader_NORM_panCK_Reference.annotations,"/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_0319_PR_Navin_Visvader_NORM_panCK_Reference.annotations.txt",sep="\t",row.names=FALSE)
# Remove metadata and annotations
rm(Visvader_0319_PR_Navin_Visvader_NORM_panCK_Reference.meta.data)
rm(Visvader_0319_PR_Navin_Visvader_NORM_panCK_Reference.annotations)

########## Generate Count Matrix
# Extract counts from the Seurat object
Visvader_0319_PR_Navin_Visvader_NORM_panCK_Reference.count <- Visvader_0319_PR_Navin_Visvader_NORM_panCK_Reference[["RNA"]]$counts
# Convert counts to a sparse matrix
Visvader_0319_PR_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as(Visvader_0319_PR_Navin_Visvader_NORM_panCK_Reference.count, 'sparseMatrix')
rm(Visvader_0319_PR_Navin_Visvader_NORM_panCK_Reference.count)

# Convert to requisite data format
Visvader_0319_PR_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as.data.frame(Visvader_0319_PR_Navin_Visvader_NORM_panCK_Reference.count.matrix)
# Save table
write.csv(Visvader_0319_PR_Navin_Visvader_NORM_panCK_Reference.count.matrix, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_Count_Matrix/Visvader_0319_PR_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv")
# Remove Objects
rm(Visvader_0319_PR_Navin_Visvader_NORM_panCK_Reference)
rm(Visvader_0319_PR_Navin_Visvader_NORM_panCK_Reference.count.matrix)

########### Run inferCNV
# Note: check.names = FALSE prevents cell names from having dashes replaced with dots; essential to match annotations file
Visvader_0319_PR_inferCNV <- CreateInfercnvObject(raw_counts_matrix=read.csv(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_Count_Matrix/Visvader_0319_PR_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv", header = TRUE, row.names = 1,check.names = FALSE),
                                                  annotations_file=read.table(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_0319_PR_Navin_Visvader_NORM_panCK_Reference.annotations.txt", header = FALSE, row.names = 1),
                                                  delim="\t",
                                                  gene_order_file="/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/hg38_gencode_v27.txt",
                                                  ref_group_names=c("Normal"))


# perform infercnv operations to reveal cnv signal
Visvader_0319_PR_inferCNV = infercnv::run(Visvader_0319_PR_inferCNV,
                                          cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                          out_dir="/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0319_PR_panCK/",  # dir is auto-created for storing outputs
                                          cluster_by_groups=T,   # cluster
                                          denoise=T,
                                          HMM=T,
                                          num_ref_groups = 8, analysis_mode = "subclusters")

###########  Integrate data with Seurat object
# Read RDS and add inferCNV results to metadata
Visvader_0319_PR_Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Visvader_0319_PR_Navin_Visvader_NORM_panCK_Reference.rds") 
Visvader_0319_PR_panCK_cnv <- infercnv::add_to_seurat(infercnv_output_path="/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0319_PR_panCK/",
                                                      seurat_obj=Visvader_0319_PR_Navin_Visvader_NORM_panCK_Reference, # optional
                                                      top_n=10)

# Collect Tumor cells from Seurat Object
Visvader_0319_PR_panCK_cnv <- subset(x = Visvader_0319_PR_panCK_cnv, subset = cnv.ident == "Tumor")
# Quantify chromosomes with alterations per cell
Visvader_0319_PR_panCK_cnv_Pos_Chromosomes <- as.data.frame(rowSums(Visvader_0319_PR_panCK_cnv@meta.data[,startsWith(names(Visvader_0319_PR_panCK_cnv@meta.data),"has_cnv_")] > 0 ))
Visvader_0319_PR_panCK_cnv_Neg_Chromosomes <- as.data.frame(rowSums(Visvader_0319_PR_panCK_cnv@meta.data[,startsWith(names(Visvader_0319_PR_panCK_cnv@meta.data),"has_cnv_")] == 0))
Visvader_0319_PR_panCK_cnv_Chromosome_Quant <- as.data.frame(c(Visvader_0319_PR_panCK_cnv_Pos_Chromosomes, Visvader_0319_PR_panCK_cnv_Neg_Chromosomes))
colnames(Visvader_0319_PR_panCK_cnv_Chromosome_Quant) <- c("Chromosomes CNV Pos", "Chromosomes CNV Neg")
write.csv(Visvader_0319_PR_panCK_cnv_Chromosome_Quant, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0319_PR_panCK/Visvader_0319_PR_panCK_cnv_Chromosome_Quant.csv")

# Save CNV Metadata
Visvader_0319_PR_panCK_cnv.meta.data <- as.data.frame(as.matrix(Visvader_0319_PR_panCK_cnv@meta.data))
write.csv(Visvader_0319_PR_panCK_cnv.meta.data, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0319_PR_panCK/Visvader_0319_PR_panCK_cnv.meta.data.csv")

# Subset Cells with inferred CNVs -------------------------------------------------------------------------------------------------------------
Visvader_0319_PR_panCK_cnv <- subset(Visvader_0319_PR_panCK_cnv, cells=rownames(Visvader_0319_PR_panCK_cnv@meta.data)[rowSums(Visvader_0319_PR_panCK_cnv@meta.data[,startsWith(names(Visvader_0319_PR_panCK_cnv@meta.data),"has_cnv_")]) > 0 ])

# Save RDS with inferredCNVs
saveRDS(Visvader_0319_PR_panCK_cnv, file = "/R/R_Visvader/Visvader_RDS/RDS_inferCNV/Visvader_0319_PR_panCK_cnv.rds")


# Free memory
rm(Visvader_0319_PR_inferCNV)
rm(Visvader_0319_PR_Navin_Visvader_NORM_panCK_Reference)
rm(Visvader_0319_PR_panCK_cnv)
rm(Visvader_0319_PR_panCK_cnv_Pos_Chromosomes)
rm(Visvader_0319_PR_panCK_cnv_Neg_Chromosomes)
rm(Visvader_0319_PR_panCK_cnv_Chromosome_Quant)
rm(Visvader_0319_PR_panCK_cnv.meta.data)
gc()



###################################################################################################
#############################     InferCNV: HER2 Tumor Specimens      #############################     
###################################################################################################

############################################
#### InferCNV: Visvader_0031_HER2_panCK ####
############################################

########## Merge with normal mammary glad, extract epithelium, and process via Seurat
# Load RDS
Visvader_0031_HER2_panCK <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0031_HER2_panCK.rds")
Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Navin_Visvader_NORM_panCK_Reference.rds")

# Set identities 
colnames(x = Visvader_0031_HER2_panCK[[]])
Idents(object = Visvader_0031_HER2_panCK) <- 'Tumor'
Visvader_0031_HER2_panCK[["cnv.ident"]] <- Idents(object = Visvader_0031_HER2_panCK)
levels(x = Visvader_0031_HER2_panCK)

colnames(x = Navin_Visvader_NORM_panCK_Reference[[]])
Idents(object = Navin_Visvader_NORM_panCK_Reference) <- 'Normal'
Navin_Visvader_NORM_panCK_Reference[["cnv.ident"]] <- Idents(object = Navin_Visvader_NORM_panCK_Reference)
levels(x = Navin_Visvader_NORM_panCK_Reference)

# Merge subsets into recombined RDS
Visvader_0031_HER2_Navin_Visvader_NORM_panCK_Reference <- merge(x = Visvader_0031_HER2_panCK, y = Navin_Visvader_NORM_panCK_Reference)
# Join Layers
Visvader_0031_HER2_Navin_Visvader_NORM_panCK_Reference <- JoinLayers(Visvader_0031_HER2_Navin_Visvader_NORM_panCK_Reference)

# Remove subsetted objects
rm(Visvader_0031_HER2_panCK)
rm(Navin_Visvader_NORM_panCK_Reference)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Visvader_0031_HER2_Navin_Visvader_NORM_panCK_Reference <- FindVariableFeatures(Visvader_0031_HER2_Navin_Visvader_NORM_panCK_Reference, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Visvader_0031_HER2_Navin_Visvader_NORM_panCK_Reference)
Visvader_0031_HER2_Navin_Visvader_NORM_panCK_Reference <- ScaleData(Visvader_0031_HER2_Navin_Visvader_NORM_panCK_Reference, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Visvader_0031_HER2_Navin_Visvader_NORM_panCK_Reference <- RunPCA(Visvader_0031_HER2_Navin_Visvader_NORM_panCK_Reference, features = VariableFeatures(object = Visvader_0031_HER2_Navin_Visvader_NORM_panCK_Reference))

# Examine and visualize PCA results a few different ways
print(Visvader_0031_HER2_Navin_Visvader_NORM_panCK_Reference[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Visvader_0031_HER2_Navin_Visvader_NORM_panCK_Reference, dims = 1:2, reduction = "pca")
DimPlot(Visvader_0031_HER2_Navin_Visvader_NORM_panCK_Reference, reduction = "pca")

Visvader_0031_HER2_Navin_Visvader_NORM_panCK_Reference <- FindNeighbors(Visvader_0031_HER2_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
Visvader_0031_HER2_Navin_Visvader_NORM_panCK_Reference <- FindClusters(Visvader_0031_HER2_Navin_Visvader_NORM_panCK_Reference, resolution = 0.1)

# Umap clustering
Visvader_0031_HER2_Navin_Visvader_NORM_panCK_Reference <- RunUMAP(Visvader_0031_HER2_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
DimPlot(Visvader_0031_HER2_Navin_Visvader_NORM_panCK_Reference, reduction = "umap")
DimPlot(Visvader_0031_HER2_Navin_Visvader_NORM_panCK_Reference, reduction = "umap", group.by = "orig.ident")

# Save RDS
saveRDS(Visvader_0031_HER2_Navin_Visvader_NORM_panCK_Reference, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Visvader_0031_HER2_Navin_Visvader_NORM_panCK_Reference.rds")

########## Prepare Cell Annotations
# Extract & Save Meta Data
Visvader_0031_HER2_Navin_Visvader_NORM_panCK_Reference.meta.data <- as.data.frame(as.matrix(Visvader_0031_HER2_Navin_Visvader_NORM_panCK_Reference@meta.data))
# Extract pertinent columns
# Make row names first column
Visvader_0031_HER2_Navin_Visvader_NORM_panCK_Reference.meta.data <- tibble::rownames_to_column(Visvader_0031_HER2_Navin_Visvader_NORM_panCK_Reference.meta.data, "Cells")
# In Seurat, "X" becomes the cell names while "cnv.ident" is specified above
Visvader_0031_HER2_Navin_Visvader_NORM_panCK_Reference.annotations <- Visvader_0031_HER2_Navin_Visvader_NORM_panCK_Reference.meta.data[, c("Cells", "cnv.ident")]
# Delete header
names(Visvader_0031_HER2_Navin_Visvader_NORM_panCK_Reference.annotations) <- NULL
# Save as table
write.table(Visvader_0031_HER2_Navin_Visvader_NORM_panCK_Reference.annotations,"/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_0031_HER2_Navin_Visvader_NORM_panCK_Reference.annotations.txt",sep="\t",row.names=FALSE)
# Remove metadata and annotations
rm(Visvader_0031_HER2_Navin_Visvader_NORM_panCK_Reference.meta.data)
rm(Visvader_0031_HER2_Navin_Visvader_NORM_panCK_Reference.annotations)

########## Generate Count Matrix
# Extract counts from the Seurat object
Visvader_0031_HER2_Navin_Visvader_NORM_panCK_Reference.count <- Visvader_0031_HER2_Navin_Visvader_NORM_panCK_Reference[["RNA"]]$counts
# Convert counts to a sparse matrix
Visvader_0031_HER2_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as(Visvader_0031_HER2_Navin_Visvader_NORM_panCK_Reference.count, 'sparseMatrix')
rm(Visvader_0031_HER2_Navin_Visvader_NORM_panCK_Reference.count)

# Convert to requisite data format
Visvader_0031_HER2_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as.data.frame(Visvader_0031_HER2_Navin_Visvader_NORM_panCK_Reference.count.matrix)
# Save table
write.csv(Visvader_0031_HER2_Navin_Visvader_NORM_panCK_Reference.count.matrix, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_Count_Matrix/Visvader_0031_HER2_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv")
# Remove Objects
rm(Visvader_0031_HER2_Navin_Visvader_NORM_panCK_Reference)
rm(Visvader_0031_HER2_Navin_Visvader_NORM_panCK_Reference.count.matrix)

########### Run inferCNV
# Note: check.names = FALSE prevents cell names from having dashes replaced with dots; essential to match annotations file
Visvader_0031_HER2_inferCNV <- CreateInfercnvObject(raw_counts_matrix=read.csv(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_Count_Matrix/Visvader_0031_HER2_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv", header = TRUE, row.names = 1,check.names = FALSE),
                                                    annotations_file=read.table(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_0031_HER2_Navin_Visvader_NORM_panCK_Reference.annotations.txt", header = FALSE, row.names = 1),
                                                    delim="\t",
                                                    gene_order_file="/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/hg38_gencode_v27.txt",
                                                    ref_group_names=c("Normal"))


# perform infercnv operations to reveal cnv signal
Visvader_0031_HER2_inferCNV = infercnv::run(Visvader_0031_HER2_inferCNV,
                                            cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                            out_dir="/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0031_HER2_panCK/",  # dir is auto-created for storing outputs
                                            cluster_by_groups=T,   # cluster
                                            denoise=T,
                                            HMM=T,
                                            num_ref_groups = 8, analysis_mode = "subclusters")

###########  Integrate data with Seurat object
# Read RDS and add inferCNV results to metadata
Visvader_0031_HER2_Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Visvader_0031_HER2_Navin_Visvader_NORM_panCK_Reference.rds") 
Visvader_0031_HER2_panCK_cnv <- infercnv::add_to_seurat(infercnv_output_path="/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0031_HER2_panCK/",
                                                        seurat_obj=Visvader_0031_HER2_Navin_Visvader_NORM_panCK_Reference, # optional
                                                        top_n=10)

# Collect Tumor cells from Seurat Object
Visvader_0031_HER2_panCK_cnv <- subset(x = Visvader_0031_HER2_panCK_cnv, subset = cnv.ident == "Tumor")
# Quantify chromosomes with alterations per cell
Visvader_0031_HER2_panCK_cnv_Pos_Chromosomes <- as.data.frame(rowSums(Visvader_0031_HER2_panCK_cnv@meta.data[,startsWith(names(Visvader_0031_HER2_panCK_cnv@meta.data),"has_cnv_")] > 0 ))
Visvader_0031_HER2_panCK_cnv_Neg_Chromosomes <- as.data.frame(rowSums(Visvader_0031_HER2_panCK_cnv@meta.data[,startsWith(names(Visvader_0031_HER2_panCK_cnv@meta.data),"has_cnv_")] == 0))
Visvader_0031_HER2_panCK_cnv_Chromosome_Quant <- as.data.frame(c(Visvader_0031_HER2_panCK_cnv_Pos_Chromosomes, Visvader_0031_HER2_panCK_cnv_Neg_Chromosomes))
colnames(Visvader_0031_HER2_panCK_cnv_Chromosome_Quant) <- c("Chromosomes CNV Pos", "Chromosomes CNV Neg")
write.csv(Visvader_0031_HER2_panCK_cnv_Chromosome_Quant, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0031_HER2_panCK/Visvader_0031_HER2_panCK_cnv_Chromosome_Quant.csv")

# Save CNV Metadata
Visvader_0031_HER2_panCK_cnv.meta.data <- as.data.frame(as.matrix(Visvader_0031_HER2_panCK_cnv@meta.data))
write.csv(Visvader_0031_HER2_panCK_cnv.meta.data, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0031_HER2_panCK/Visvader_0031_HER2_panCK_cnv.meta.data.csv")

# Subset Cells with inferred CNVs -------------------------------------------------------------------------------------------------------------
Visvader_0031_HER2_panCK_cnv <- subset(Visvader_0031_HER2_panCK_cnv, cells=rownames(Visvader_0031_HER2_panCK_cnv@meta.data)[rowSums(Visvader_0031_HER2_panCK_cnv@meta.data[,startsWith(names(Visvader_0031_HER2_panCK_cnv@meta.data),"has_cnv_")]) > 0 ])

# Save RDS with inferredCNVs
saveRDS(Visvader_0031_HER2_panCK_cnv, file = "/R/R_Visvader/Visvader_RDS/RDS_inferCNV/Visvader_0031_HER2_panCK_cnv.rds")

# Free memory
rm(Visvader_0031_HER2_inferCNV)
rm(Visvader_0031_HER2_Navin_Visvader_NORM_panCK_Reference)
rm(Visvader_0031_HER2_panCK_cnv)
rm(Visvader_0031_HER2_panCK_cnv_Pos_Chromosomes)
rm(Visvader_0031_HER2_panCK_cnv_Neg_Chromosomes)
rm(Visvader_0031_HER2_panCK_cnv_Chromosome_Quant)
rm(Visvader_0031_HER2_panCK_cnv.meta.data)
gc()


############################################
#### InferCNV: Visvader_0069_HER2_panCK ####
############################################

########## Merge with normal mammary glad, extract epithelium, and process via Seurat
# Load RDS
Visvader_0069_HER2_panCK <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0069_HER2_panCK.rds")
Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Navin_Visvader_NORM_panCK_Reference.rds")

# Set identities 
colnames(x = Visvader_0069_HER2_panCK[[]])
Idents(object = Visvader_0069_HER2_panCK) <- 'Tumor'
Visvader_0069_HER2_panCK[["cnv.ident"]] <- Idents(object = Visvader_0069_HER2_panCK)
levels(x = Visvader_0069_HER2_panCK)

colnames(x = Navin_Visvader_NORM_panCK_Reference[[]])
Idents(object = Navin_Visvader_NORM_panCK_Reference) <- 'Normal'
Navin_Visvader_NORM_panCK_Reference[["cnv.ident"]] <- Idents(object = Navin_Visvader_NORM_panCK_Reference)
levels(x = Navin_Visvader_NORM_panCK_Reference)

# Merge subsets into recombined RDS
Visvader_0069_HER2_Navin_Visvader_NORM_panCK_Reference <- merge(x = Visvader_0069_HER2_panCK, y = Navin_Visvader_NORM_panCK_Reference)
# Join Layers
Visvader_0069_HER2_Navin_Visvader_NORM_panCK_Reference <- JoinLayers(Visvader_0069_HER2_Navin_Visvader_NORM_panCK_Reference)

# Remove subsetted objects
rm(Visvader_0069_HER2_panCK)
rm(Navin_Visvader_NORM_panCK_Reference)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Visvader_0069_HER2_Navin_Visvader_NORM_panCK_Reference <- FindVariableFeatures(Visvader_0069_HER2_Navin_Visvader_NORM_panCK_Reference, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Visvader_0069_HER2_Navin_Visvader_NORM_panCK_Reference)
Visvader_0069_HER2_Navin_Visvader_NORM_panCK_Reference <- ScaleData(Visvader_0069_HER2_Navin_Visvader_NORM_panCK_Reference, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Visvader_0069_HER2_Navin_Visvader_NORM_panCK_Reference <- RunPCA(Visvader_0069_HER2_Navin_Visvader_NORM_panCK_Reference, features = VariableFeatures(object = Visvader_0069_HER2_Navin_Visvader_NORM_panCK_Reference))

# Examine and visualize PCA results a few different ways
print(Visvader_0069_HER2_Navin_Visvader_NORM_panCK_Reference[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Visvader_0069_HER2_Navin_Visvader_NORM_panCK_Reference, dims = 1:2, reduction = "pca")
DimPlot(Visvader_0069_HER2_Navin_Visvader_NORM_panCK_Reference, reduction = "pca")

Visvader_0069_HER2_Navin_Visvader_NORM_panCK_Reference <- FindNeighbors(Visvader_0069_HER2_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
Visvader_0069_HER2_Navin_Visvader_NORM_panCK_Reference <- FindClusters(Visvader_0069_HER2_Navin_Visvader_NORM_panCK_Reference, resolution = 0.1)

# Umap clustering
Visvader_0069_HER2_Navin_Visvader_NORM_panCK_Reference <- RunUMAP(Visvader_0069_HER2_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
DimPlot(Visvader_0069_HER2_Navin_Visvader_NORM_panCK_Reference, reduction = "umap")
DimPlot(Visvader_0069_HER2_Navin_Visvader_NORM_panCK_Reference, reduction = "umap", group.by = "orig.ident")

# Save RDS
saveRDS(Visvader_0069_HER2_Navin_Visvader_NORM_panCK_Reference, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Visvader_0069_HER2_Navin_Visvader_NORM_panCK_Reference.rds")

########## Prepare Cell Annotations
# Extract & Save Meta Data
Visvader_0069_HER2_Navin_Visvader_NORM_panCK_Reference.meta.data <- as.data.frame(as.matrix(Visvader_0069_HER2_Navin_Visvader_NORM_panCK_Reference@meta.data))
# Extract pertinent columns
# Make row names first column
Visvader_0069_HER2_Navin_Visvader_NORM_panCK_Reference.meta.data <- tibble::rownames_to_column(Visvader_0069_HER2_Navin_Visvader_NORM_panCK_Reference.meta.data, "Cells")
# In Seurat, "X" becomes the cell names while "cnv.ident" is specified above
Visvader_0069_HER2_Navin_Visvader_NORM_panCK_Reference.annotations <- Visvader_0069_HER2_Navin_Visvader_NORM_panCK_Reference.meta.data[, c("Cells", "cnv.ident")]
# Delete header
names(Visvader_0069_HER2_Navin_Visvader_NORM_panCK_Reference.annotations) <- NULL
# Save as table
write.table(Visvader_0069_HER2_Navin_Visvader_NORM_panCK_Reference.annotations,"/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_0069_HER2_Navin_Visvader_NORM_panCK_Reference.annotations.txt",sep="\t",row.names=FALSE)
# Remove metadata and annotations
rm(Visvader_0069_HER2_Navin_Visvader_NORM_panCK_Reference.meta.data)
rm(Visvader_0069_HER2_Navin_Visvader_NORM_panCK_Reference.annotations)

########## Generate Count Matrix
# Extract counts from the Seurat object
Visvader_0069_HER2_Navin_Visvader_NORM_panCK_Reference.count <- Visvader_0069_HER2_Navin_Visvader_NORM_panCK_Reference[["RNA"]]$counts
# Convert counts to a sparse matrix
Visvader_0069_HER2_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as(Visvader_0069_HER2_Navin_Visvader_NORM_panCK_Reference.count, 'sparseMatrix')
rm(Visvader_0069_HER2_Navin_Visvader_NORM_panCK_Reference.count)

# Convert to requisite data format
Visvader_0069_HER2_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as.data.frame(Visvader_0069_HER2_Navin_Visvader_NORM_panCK_Reference.count.matrix)
# Save table
write.csv(Visvader_0069_HER2_Navin_Visvader_NORM_panCK_Reference.count.matrix, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_Count_Matrix/Visvader_0069_HER2_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv")
# Remove Objects
rm(Visvader_0069_HER2_Navin_Visvader_NORM_panCK_Reference)
rm(Visvader_0069_HER2_Navin_Visvader_NORM_panCK_Reference.count.matrix)

########### Run inferCNV
# Note: check.names = FALSE prevents cell names from having dashes replaced with dots; essential to match annotations file
Visvader_0069_HER2_inferCNV <- CreateInfercnvObject(raw_counts_matrix=read.csv(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_Count_Matrix/Visvader_0069_HER2_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv", header = TRUE, row.names = 1,check.names = FALSE),
                                                    annotations_file=read.table(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_0069_HER2_Navin_Visvader_NORM_panCK_Reference.annotations.txt", header = FALSE, row.names = 1),
                                                    delim="\t",
                                                    gene_order_file="/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/hg38_gencode_v27.txt",
                                                    ref_group_names=c("Normal"))


# perform infercnv operations to reveal cnv signal
Visvader_0069_HER2_inferCNV = infercnv::run(Visvader_0069_HER2_inferCNV,
                                            cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                            out_dir="/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0069_HER2_panCK/",  # dir is auto-created for storing outputs
                                            cluster_by_groups=T,   # cluster
                                            denoise=T,
                                            HMM=T,
                                            num_ref_groups = 8, analysis_mode = "subclusters")

###########  Integrate data with Seurat object
# Read RDS and add inferCNV results to metadata
Visvader_0069_HER2_Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Visvader_0069_HER2_Navin_Visvader_NORM_panCK_Reference.rds") 
Visvader_0069_HER2_panCK_cnv <- infercnv::add_to_seurat(infercnv_output_path="/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0069_HER2_panCK/",
                                                        seurat_obj=Visvader_0069_HER2_Navin_Visvader_NORM_panCK_Reference, # optional
                                                        top_n=10)

# Collect Tumor cells from Seurat Object
Visvader_0069_HER2_panCK_cnv <- subset(x = Visvader_0069_HER2_panCK_cnv, subset = cnv.ident == "Tumor")
# Quantify chromosomes with alterations per cell
Visvader_0069_HER2_panCK_cnv_Pos_Chromosomes <- as.data.frame(rowSums(Visvader_0069_HER2_panCK_cnv@meta.data[,startsWith(names(Visvader_0069_HER2_panCK_cnv@meta.data),"has_cnv_")] > 0 ))
Visvader_0069_HER2_panCK_cnv_Neg_Chromosomes <- as.data.frame(rowSums(Visvader_0069_HER2_panCK_cnv@meta.data[,startsWith(names(Visvader_0069_HER2_panCK_cnv@meta.data),"has_cnv_")] == 0))
Visvader_0069_HER2_panCK_cnv_Chromosome_Quant <- as.data.frame(c(Visvader_0069_HER2_panCK_cnv_Pos_Chromosomes, Visvader_0069_HER2_panCK_cnv_Neg_Chromosomes))
colnames(Visvader_0069_HER2_panCK_cnv_Chromosome_Quant) <- c("Chromosomes CNV Pos", "Chromosomes CNV Neg")
write.csv(Visvader_0069_HER2_panCK_cnv_Chromosome_Quant, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0069_HER2_panCK/Visvader_0069_HER2_panCK_cnv_Chromosome_Quant.csv")

# Save CNV Metadata
Visvader_0069_HER2_panCK_cnv.meta.data <- as.data.frame(as.matrix(Visvader_0069_HER2_panCK_cnv@meta.data))
write.csv(Visvader_0069_HER2_panCK_cnv.meta.data, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0069_HER2_panCK/Visvader_0069_HER2_panCK_cnv.meta.data.csv")

# Subset Cells with inferred CNVs -------------------------------------------------------------------------------------------------------------
Visvader_0069_HER2_panCK_cnv <- subset(Visvader_0069_HER2_panCK_cnv, cells=rownames(Visvader_0069_HER2_panCK_cnv@meta.data)[rowSums(Visvader_0069_HER2_panCK_cnv@meta.data[,startsWith(names(Visvader_0069_HER2_panCK_cnv@meta.data),"has_cnv_")]) > 0 ])

# Save RDS with inferredCNVs
saveRDS(Visvader_0069_HER2_panCK_cnv, file = "/R/R_Visvader/Visvader_RDS/RDS_inferCNV/Visvader_0069_HER2_panCK_cnv.rds")


# Free memory
rm(Visvader_0069_HER2_inferCNV)
rm(Visvader_0069_HER2_Navin_Visvader_NORM_panCK_Reference)
rm(Visvader_0069_HER2_panCK_cnv)
rm(Visvader_0069_HER2_panCK_cnv_Pos_Chromosomes)
rm(Visvader_0069_HER2_panCK_cnv_Neg_Chromosomes)
rm(Visvader_0069_HER2_panCK_cnv_Chromosome_Quant)
rm(Visvader_0069_HER2_panCK_cnv.meta.data)
gc()


############################################
#### InferCNV: Visvader_0161_HER2_panCK ####
############################################

########## Merge with normal mammary glad, extract epithelium, and process via Seurat
# Load RDS
Visvader_0161_HER2_panCK <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0161_HER2_panCK.rds")
Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Navin_Visvader_NORM_panCK_Reference.rds")

# Set identities 
colnames(x = Visvader_0161_HER2_panCK[[]])
Idents(object = Visvader_0161_HER2_panCK) <- 'Tumor'
Visvader_0161_HER2_panCK[["cnv.ident"]] <- Idents(object = Visvader_0161_HER2_panCK)
levels(x = Visvader_0161_HER2_panCK)

colnames(x = Navin_Visvader_NORM_panCK_Reference[[]])
Idents(object = Navin_Visvader_NORM_panCK_Reference) <- 'Normal'
Navin_Visvader_NORM_panCK_Reference[["cnv.ident"]] <- Idents(object = Navin_Visvader_NORM_panCK_Reference)
levels(x = Navin_Visvader_NORM_panCK_Reference)

# Merge subsets into recombined RDS
Visvader_0161_HER2_Navin_Visvader_NORM_panCK_Reference <- merge(x = Visvader_0161_HER2_panCK, y = Navin_Visvader_NORM_panCK_Reference)
# Join Layers
Visvader_0161_HER2_Navin_Visvader_NORM_panCK_Reference <- JoinLayers(Visvader_0161_HER2_Navin_Visvader_NORM_panCK_Reference)

# Remove subsetted objects
rm(Visvader_0161_HER2_panCK)
rm(Navin_Visvader_NORM_panCK_Reference)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Visvader_0161_HER2_Navin_Visvader_NORM_panCK_Reference <- FindVariableFeatures(Visvader_0161_HER2_Navin_Visvader_NORM_panCK_Reference, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Visvader_0161_HER2_Navin_Visvader_NORM_panCK_Reference)
Visvader_0161_HER2_Navin_Visvader_NORM_panCK_Reference <- ScaleData(Visvader_0161_HER2_Navin_Visvader_NORM_panCK_Reference, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Visvader_0161_HER2_Navin_Visvader_NORM_panCK_Reference <- RunPCA(Visvader_0161_HER2_Navin_Visvader_NORM_panCK_Reference, features = VariableFeatures(object = Visvader_0161_HER2_Navin_Visvader_NORM_panCK_Reference))

# Examine and visualize PCA results a few different ways
print(Visvader_0161_HER2_Navin_Visvader_NORM_panCK_Reference[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Visvader_0161_HER2_Navin_Visvader_NORM_panCK_Reference, dims = 1:2, reduction = "pca")
DimPlot(Visvader_0161_HER2_Navin_Visvader_NORM_panCK_Reference, reduction = "pca")

Visvader_0161_HER2_Navin_Visvader_NORM_panCK_Reference <- FindNeighbors(Visvader_0161_HER2_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
Visvader_0161_HER2_Navin_Visvader_NORM_panCK_Reference <- FindClusters(Visvader_0161_HER2_Navin_Visvader_NORM_panCK_Reference, resolution = 0.1)

# Umap clustering
Visvader_0161_HER2_Navin_Visvader_NORM_panCK_Reference <- RunUMAP(Visvader_0161_HER2_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
DimPlot(Visvader_0161_HER2_Navin_Visvader_NORM_panCK_Reference, reduction = "umap")
DimPlot(Visvader_0161_HER2_Navin_Visvader_NORM_panCK_Reference, reduction = "umap", group.by = "orig.ident")

# Save RDS
saveRDS(Visvader_0161_HER2_Navin_Visvader_NORM_panCK_Reference, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Visvader_0161_HER2_Navin_Visvader_NORM_panCK_Reference.rds")

########## Prepare Cell Annotations
# Extract & Save Meta Data
Visvader_0161_HER2_Navin_Visvader_NORM_panCK_Reference.meta.data <- as.data.frame(as.matrix(Visvader_0161_HER2_Navin_Visvader_NORM_panCK_Reference@meta.data))
# Extract pertinent columns
# Make row names first column
Visvader_0161_HER2_Navin_Visvader_NORM_panCK_Reference.meta.data <- tibble::rownames_to_column(Visvader_0161_HER2_Navin_Visvader_NORM_panCK_Reference.meta.data, "Cells")
# In Seurat, "X" becomes the cell names while "cnv.ident" is specified above
Visvader_0161_HER2_Navin_Visvader_NORM_panCK_Reference.annotations <- Visvader_0161_HER2_Navin_Visvader_NORM_panCK_Reference.meta.data[, c("Cells", "cnv.ident")]
# Delete header
names(Visvader_0161_HER2_Navin_Visvader_NORM_panCK_Reference.annotations) <- NULL
# Save as table
write.table(Visvader_0161_HER2_Navin_Visvader_NORM_panCK_Reference.annotations,"/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_0161_HER2_Navin_Visvader_NORM_panCK_Reference.annotations.txt",sep="\t",row.names=FALSE)
# Remove metadata and annotations
rm(Visvader_0161_HER2_Navin_Visvader_NORM_panCK_Reference.meta.data)
rm(Visvader_0161_HER2_Navin_Visvader_NORM_panCK_Reference.annotations)

########## Generate Count Matrix
# Extract counts from the Seurat object
Visvader_0161_HER2_Navin_Visvader_NORM_panCK_Reference.count <- Visvader_0161_HER2_Navin_Visvader_NORM_panCK_Reference[["RNA"]]$counts
# Convert counts to a sparse matrix
Visvader_0161_HER2_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as(Visvader_0161_HER2_Navin_Visvader_NORM_panCK_Reference.count, 'sparseMatrix')
rm(Visvader_0161_HER2_Navin_Visvader_NORM_panCK_Reference.count)

# Convert to requisite data format
Visvader_0161_HER2_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as.data.frame(Visvader_0161_HER2_Navin_Visvader_NORM_panCK_Reference.count.matrix)
# Save table
write.csv(Visvader_0161_HER2_Navin_Visvader_NORM_panCK_Reference.count.matrix, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_Count_Matrix/Visvader_0161_HER2_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv")
# Remove Objects
rm(Visvader_0161_HER2_Navin_Visvader_NORM_panCK_Reference)
rm(Visvader_0161_HER2_Navin_Visvader_NORM_panCK_Reference.count.matrix)

########### Run inferCNV
# Note: check.names = FALSE prevents cell names from having dashes replaced with dots; essential to match annotations file
Visvader_0161_HER2_inferCNV <- CreateInfercnvObject(raw_counts_matrix=read.csv(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_Count_Matrix/Visvader_0161_HER2_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv", header = TRUE, row.names = 1,check.names = FALSE),
                                                    annotations_file=read.table(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_0161_HER2_Navin_Visvader_NORM_panCK_Reference.annotations.txt", header = FALSE, row.names = 1),
                                                    delim="\t",
                                                    gene_order_file="/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/hg38_gencode_v27.txt",
                                                    ref_group_names=c("Normal"))


# perform infercnv operations to reveal cnv signal
Visvader_0161_HER2_inferCNV = infercnv::run(Visvader_0161_HER2_inferCNV,
                                            cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                            out_dir="/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0161_HER2_panCK/",  # dir is auto-created for storing outputs
                                            cluster_by_groups=T,   # cluster
                                            denoise=T,
                                            HMM=T,
                                            num_ref_groups = 8, analysis_mode = "subclusters")

###########  Integrate data with Seurat object
# Read RDS and add inferCNV results to metadata
Visvader_0161_HER2_Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Visvader_0161_HER2_Navin_Visvader_NORM_panCK_Reference.rds") 
Visvader_0161_HER2_panCK_cnv <- infercnv::add_to_seurat(infercnv_output_path="/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0161_HER2_panCK/",
                                                        seurat_obj=Visvader_0161_HER2_Navin_Visvader_NORM_panCK_Reference, # optional
                                                        top_n=10)

# Collect Tumor cells from Seurat Object
Visvader_0161_HER2_panCK_cnv <- subset(x = Visvader_0161_HER2_panCK_cnv, subset = cnv.ident == "Tumor")
# Quantify chromosomes with alterations per cell
Visvader_0161_HER2_panCK_cnv_Pos_Chromosomes <- as.data.frame(rowSums(Visvader_0161_HER2_panCK_cnv@meta.data[,startsWith(names(Visvader_0161_HER2_panCK_cnv@meta.data),"has_cnv_")] > 0 ))
Visvader_0161_HER2_panCK_cnv_Neg_Chromosomes <- as.data.frame(rowSums(Visvader_0161_HER2_panCK_cnv@meta.data[,startsWith(names(Visvader_0161_HER2_panCK_cnv@meta.data),"has_cnv_")] == 0))
Visvader_0161_HER2_panCK_cnv_Chromosome_Quant <- as.data.frame(c(Visvader_0161_HER2_panCK_cnv_Pos_Chromosomes, Visvader_0161_HER2_panCK_cnv_Neg_Chromosomes))
colnames(Visvader_0161_HER2_panCK_cnv_Chromosome_Quant) <- c("Chromosomes CNV Pos", "Chromosomes CNV Neg")
write.csv(Visvader_0161_HER2_panCK_cnv_Chromosome_Quant, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0161_HER2_panCK/Visvader_0161_HER2_panCK_cnv_Chromosome_Quant.csv")

# Save CNV Metadata
Visvader_0161_HER2_panCK_cnv.meta.data <- as.data.frame(as.matrix(Visvader_0161_HER2_panCK_cnv@meta.data))
write.csv(Visvader_0161_HER2_panCK_cnv.meta.data, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0161_HER2_panCK/Visvader_0161_HER2_panCK_cnv.meta.data.csv")

# Subset Cells with inferred CNVs -------------------------------------------------------------------------------------------------------------
Visvader_0161_HER2_panCK_cnv <- subset(Visvader_0161_HER2_panCK_cnv, cells=rownames(Visvader_0161_HER2_panCK_cnv@meta.data)[rowSums(Visvader_0161_HER2_panCK_cnv@meta.data[,startsWith(names(Visvader_0161_HER2_panCK_cnv@meta.data),"has_cnv_")]) > 0 ])

# Save RDS with inferredCNVs
saveRDS(Visvader_0161_HER2_panCK_cnv, file = "/R/R_Visvader/Visvader_RDS/RDS_inferCNV/Visvader_0161_HER2_panCK_cnv.rds")


# Free memory
rm(Visvader_0161_HER2_inferCNV)
rm(Visvader_0161_HER2_Navin_Visvader_NORM_panCK_Reference)
rm(Visvader_0161_HER2_panCK_cnv)
rm(Visvader_0161_HER2_panCK_cnv_Pos_Chromosomes)
rm(Visvader_0161_HER2_panCK_cnv_Neg_Chromosomes)
rm(Visvader_0161_HER2_panCK_cnv_Chromosome_Quant)
rm(Visvader_0161_HER2_panCK_cnv.meta.data)
gc()


############################################
#### InferCNV: Visvader_0176_HER2_panCK ####
############################################

########## Merge with normal mammary glad, extract epithelium, and process via Seurat
# Load RDS
Visvader_0176_HER2_panCK <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0176_HER2_panCK.rds")
Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Navin_Visvader_NORM_panCK_Reference.rds")

# Set identities 
colnames(x = Visvader_0176_HER2_panCK[[]])
Idents(object = Visvader_0176_HER2_panCK) <- 'Tumor'
Visvader_0176_HER2_panCK[["cnv.ident"]] <- Idents(object = Visvader_0176_HER2_panCK)
levels(x = Visvader_0176_HER2_panCK)

colnames(x = Navin_Visvader_NORM_panCK_Reference[[]])
Idents(object = Navin_Visvader_NORM_panCK_Reference) <- 'Normal'
Navin_Visvader_NORM_panCK_Reference[["cnv.ident"]] <- Idents(object = Navin_Visvader_NORM_panCK_Reference)
levels(x = Navin_Visvader_NORM_panCK_Reference)

# Merge subsets into recombined RDS
Visvader_0176_HER2_Navin_Visvader_NORM_panCK_Reference <- merge(x = Visvader_0176_HER2_panCK, y = Navin_Visvader_NORM_panCK_Reference)
# Join Layers
Visvader_0176_HER2_Navin_Visvader_NORM_panCK_Reference <- JoinLayers(Visvader_0176_HER2_Navin_Visvader_NORM_panCK_Reference)

# Remove subsetted objects
rm(Visvader_0176_HER2_panCK)
rm(Navin_Visvader_NORM_panCK_Reference)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Visvader_0176_HER2_Navin_Visvader_NORM_panCK_Reference <- FindVariableFeatures(Visvader_0176_HER2_Navin_Visvader_NORM_panCK_Reference, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Visvader_0176_HER2_Navin_Visvader_NORM_panCK_Reference)
Visvader_0176_HER2_Navin_Visvader_NORM_panCK_Reference <- ScaleData(Visvader_0176_HER2_Navin_Visvader_NORM_panCK_Reference, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Visvader_0176_HER2_Navin_Visvader_NORM_panCK_Reference <- RunPCA(Visvader_0176_HER2_Navin_Visvader_NORM_panCK_Reference, features = VariableFeatures(object = Visvader_0176_HER2_Navin_Visvader_NORM_panCK_Reference))

# Examine and visualize PCA results a few different ways
print(Visvader_0176_HER2_Navin_Visvader_NORM_panCK_Reference[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Visvader_0176_HER2_Navin_Visvader_NORM_panCK_Reference, dims = 1:2, reduction = "pca")
DimPlot(Visvader_0176_HER2_Navin_Visvader_NORM_panCK_Reference, reduction = "pca")

Visvader_0176_HER2_Navin_Visvader_NORM_panCK_Reference <- FindNeighbors(Visvader_0176_HER2_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
Visvader_0176_HER2_Navin_Visvader_NORM_panCK_Reference <- FindClusters(Visvader_0176_HER2_Navin_Visvader_NORM_panCK_Reference, resolution = 0.1)

# Umap clustering
Visvader_0176_HER2_Navin_Visvader_NORM_panCK_Reference <- RunUMAP(Visvader_0176_HER2_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
DimPlot(Visvader_0176_HER2_Navin_Visvader_NORM_panCK_Reference, reduction = "umap")
DimPlot(Visvader_0176_HER2_Navin_Visvader_NORM_panCK_Reference, reduction = "umap", group.by = "orig.ident")

# Save RDS
saveRDS(Visvader_0176_HER2_Navin_Visvader_NORM_panCK_Reference, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Visvader_0176_HER2_Navin_Visvader_NORM_panCK_Reference.rds")

########## Prepare Cell Annotations
# Extract & Save Meta Data
Visvader_0176_HER2_Navin_Visvader_NORM_panCK_Reference.meta.data <- as.data.frame(as.matrix(Visvader_0176_HER2_Navin_Visvader_NORM_panCK_Reference@meta.data))
# Extract pertinent columns
# Make row names first column
Visvader_0176_HER2_Navin_Visvader_NORM_panCK_Reference.meta.data <- tibble::rownames_to_column(Visvader_0176_HER2_Navin_Visvader_NORM_panCK_Reference.meta.data, "Cells")
# In Seurat, "X" becomes the cell names while "cnv.ident" is specified above
Visvader_0176_HER2_Navin_Visvader_NORM_panCK_Reference.annotations <- Visvader_0176_HER2_Navin_Visvader_NORM_panCK_Reference.meta.data[, c("Cells", "cnv.ident")]
# Delete header
names(Visvader_0176_HER2_Navin_Visvader_NORM_panCK_Reference.annotations) <- NULL
# Save as table
write.table(Visvader_0176_HER2_Navin_Visvader_NORM_panCK_Reference.annotations,"/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_0176_HER2_Navin_Visvader_NORM_panCK_Reference.annotations.txt",sep="\t",row.names=FALSE)
# Remove metadata and annotations
rm(Visvader_0176_HER2_Navin_Visvader_NORM_panCK_Reference.meta.data)
rm(Visvader_0176_HER2_Navin_Visvader_NORM_panCK_Reference.annotations)

########## Generate Count Matrix
# Extract counts from the Seurat object
Visvader_0176_HER2_Navin_Visvader_NORM_panCK_Reference.count <- Visvader_0176_HER2_Navin_Visvader_NORM_panCK_Reference[["RNA"]]$counts
# Convert counts to a sparse matrix
Visvader_0176_HER2_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as(Visvader_0176_HER2_Navin_Visvader_NORM_panCK_Reference.count, 'sparseMatrix')
rm(Visvader_0176_HER2_Navin_Visvader_NORM_panCK_Reference.count)

# Convert to requisite data format
Visvader_0176_HER2_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as.data.frame(Visvader_0176_HER2_Navin_Visvader_NORM_panCK_Reference.count.matrix)
# Save table
write.csv(Visvader_0176_HER2_Navin_Visvader_NORM_panCK_Reference.count.matrix, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_Count_Matrix/Visvader_0176_HER2_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv")
# Remove Objects
rm(Visvader_0176_HER2_Navin_Visvader_NORM_panCK_Reference)
rm(Visvader_0176_HER2_Navin_Visvader_NORM_panCK_Reference.count.matrix)

########### Run inferCNV
# Note: check.names = FALSE prevents cell names from having dashes replaced with dots; essential to match annotations file
Visvader_0176_HER2_inferCNV <- CreateInfercnvObject(raw_counts_matrix=read.csv(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_Count_Matrix/Visvader_0176_HER2_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv", header = TRUE, row.names = 1,check.names = FALSE),
                                                    annotations_file=read.table(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_0176_HER2_Navin_Visvader_NORM_panCK_Reference.annotations.txt", header = FALSE, row.names = 1),
                                                    delim="\t",
                                                    gene_order_file="/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/hg38_gencode_v27.txt",
                                                    ref_group_names=c("Normal"))


# perform infercnv operations to reveal cnv signal
Visvader_0176_HER2_inferCNV = infercnv::run(Visvader_0176_HER2_inferCNV,
                                            cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                            out_dir="/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0176_HER2_panCK/",  # dir is auto-created for storing outputs
                                            cluster_by_groups=T,   # cluster
                                            denoise=T,
                                            HMM=T,
                                            num_ref_groups = 8, analysis_mode = "subclusters")

###########  Integrate data with Seurat object
# Read RDS and add inferCNV results to metadata
Visvader_0176_HER2_Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Visvader_0176_HER2_Navin_Visvader_NORM_panCK_Reference.rds") 
Visvader_0176_HER2_panCK_cnv <- infercnv::add_to_seurat(infercnv_output_path="/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0176_HER2_panCK/",
                                                        seurat_obj=Visvader_0176_HER2_Navin_Visvader_NORM_panCK_Reference, # optional
                                                        top_n=10)

# Collect Tumor cells from Seurat Object
Visvader_0176_HER2_panCK_cnv <- subset(x = Visvader_0176_HER2_panCK_cnv, subset = cnv.ident == "Tumor")
# Quantify chromosomes with alterations per cell
Visvader_0176_HER2_panCK_cnv_Pos_Chromosomes <- as.data.frame(rowSums(Visvader_0176_HER2_panCK_cnv@meta.data[,startsWith(names(Visvader_0176_HER2_panCK_cnv@meta.data),"has_cnv_")] > 0 ))
Visvader_0176_HER2_panCK_cnv_Neg_Chromosomes <- as.data.frame(rowSums(Visvader_0176_HER2_panCK_cnv@meta.data[,startsWith(names(Visvader_0176_HER2_panCK_cnv@meta.data),"has_cnv_")] == 0))
Visvader_0176_HER2_panCK_cnv_Chromosome_Quant <- as.data.frame(c(Visvader_0176_HER2_panCK_cnv_Pos_Chromosomes, Visvader_0176_HER2_panCK_cnv_Neg_Chromosomes))
colnames(Visvader_0176_HER2_panCK_cnv_Chromosome_Quant) <- c("Chromosomes CNV Pos", "Chromosomes CNV Neg")
write.csv(Visvader_0176_HER2_panCK_cnv_Chromosome_Quant, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0176_HER2_panCK/Visvader_0176_HER2_panCK_cnv_Chromosome_Quant.csv")

# Save CNV Metadata
Visvader_0176_HER2_panCK_cnv.meta.data <- as.data.frame(as.matrix(Visvader_0176_HER2_panCK_cnv@meta.data))
write.csv(Visvader_0176_HER2_panCK_cnv.meta.data, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0176_HER2_panCK/Visvader_0176_HER2_panCK_cnv.meta.data.csv")

# Subset Cells with inferred CNVs -------------------------------------------------------------------------------------------------------------
Visvader_0176_HER2_panCK_cnv <- subset(Visvader_0176_HER2_panCK_cnv, cells=rownames(Visvader_0176_HER2_panCK_cnv@meta.data)[rowSums(Visvader_0176_HER2_panCK_cnv@meta.data[,startsWith(names(Visvader_0176_HER2_panCK_cnv@meta.data),"has_cnv_")]) > 0 ])

# Save RDS with inferredCNVs
saveRDS(Visvader_0176_HER2_panCK_cnv, file = "/R/R_Visvader/Visvader_RDS/RDS_inferCNV/Visvader_0176_HER2_panCK_cnv.rds")


# Free memory
rm(Visvader_0176_HER2_inferCNV)
rm(Visvader_0176_HER2_Navin_Visvader_NORM_panCK_Reference)
rm(Visvader_0176_HER2_panCK_cnv)
rm(Visvader_0176_HER2_panCK_cnv_Pos_Chromosomes)
rm(Visvader_0176_HER2_panCK_cnv_Neg_Chromosomes)
rm(Visvader_0176_HER2_panCK_cnv_Chromosome_Quant)
rm(Visvader_0176_HER2_panCK_cnv.meta.data)
gc()


############################################
#### InferCNV: Visvader_0308_HER2_panCK ####
############################################

########## Merge with normal mammary glad, extract epithelium, and process via Seurat
# Load RDS
Visvader_0308_HER2_panCK <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0308_HER2_panCK.rds")
Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Navin_Visvader_NORM_panCK_Reference.rds")

# Set identities 
colnames(x = Visvader_0308_HER2_panCK[[]])
Idents(object = Visvader_0308_HER2_panCK) <- 'Tumor'
Visvader_0308_HER2_panCK[["cnv.ident"]] <- Idents(object = Visvader_0308_HER2_panCK)
levels(x = Visvader_0308_HER2_panCK)

colnames(x = Navin_Visvader_NORM_panCK_Reference[[]])
Idents(object = Navin_Visvader_NORM_panCK_Reference) <- 'Normal'
Navin_Visvader_NORM_panCK_Reference[["cnv.ident"]] <- Idents(object = Navin_Visvader_NORM_panCK_Reference)
levels(x = Navin_Visvader_NORM_panCK_Reference)

# Merge subsets into recombined RDS
Visvader_0308_HER2_Navin_Visvader_NORM_panCK_Reference <- merge(x = Visvader_0308_HER2_panCK, y = Navin_Visvader_NORM_panCK_Reference)
# Join Layers
Visvader_0308_HER2_Navin_Visvader_NORM_panCK_Reference <- JoinLayers(Visvader_0308_HER2_Navin_Visvader_NORM_panCK_Reference)

# Remove subsetted objects
rm(Visvader_0308_HER2_panCK)
rm(Navin_Visvader_NORM_panCK_Reference)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Visvader_0308_HER2_Navin_Visvader_NORM_panCK_Reference <- FindVariableFeatures(Visvader_0308_HER2_Navin_Visvader_NORM_panCK_Reference, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Visvader_0308_HER2_Navin_Visvader_NORM_panCK_Reference)
Visvader_0308_HER2_Navin_Visvader_NORM_panCK_Reference <- ScaleData(Visvader_0308_HER2_Navin_Visvader_NORM_panCK_Reference, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Visvader_0308_HER2_Navin_Visvader_NORM_panCK_Reference <- RunPCA(Visvader_0308_HER2_Navin_Visvader_NORM_panCK_Reference, features = VariableFeatures(object = Visvader_0308_HER2_Navin_Visvader_NORM_panCK_Reference))

# Examine and visualize PCA results a few different ways
print(Visvader_0308_HER2_Navin_Visvader_NORM_panCK_Reference[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Visvader_0308_HER2_Navin_Visvader_NORM_panCK_Reference, dims = 1:2, reduction = "pca")
DimPlot(Visvader_0308_HER2_Navin_Visvader_NORM_panCK_Reference, reduction = "pca")

Visvader_0308_HER2_Navin_Visvader_NORM_panCK_Reference <- FindNeighbors(Visvader_0308_HER2_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
Visvader_0308_HER2_Navin_Visvader_NORM_panCK_Reference <- FindClusters(Visvader_0308_HER2_Navin_Visvader_NORM_panCK_Reference, resolution = 0.1)

# Umap clustering
Visvader_0308_HER2_Navin_Visvader_NORM_panCK_Reference <- RunUMAP(Visvader_0308_HER2_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
DimPlot(Visvader_0308_HER2_Navin_Visvader_NORM_panCK_Reference, reduction = "umap")
DimPlot(Visvader_0308_HER2_Navin_Visvader_NORM_panCK_Reference, reduction = "umap", group.by = "orig.ident")

# Save RDS
saveRDS(Visvader_0308_HER2_Navin_Visvader_NORM_panCK_Reference, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Visvader_0308_HER2_Navin_Visvader_NORM_panCK_Reference.rds")

########## Prepare Cell Annotations
# Extract & Save Meta Data
Visvader_0308_HER2_Navin_Visvader_NORM_panCK_Reference.meta.data <- as.data.frame(as.matrix(Visvader_0308_HER2_Navin_Visvader_NORM_panCK_Reference@meta.data))
# Extract pertinent columns
# Make row names first column
Visvader_0308_HER2_Navin_Visvader_NORM_panCK_Reference.meta.data <- tibble::rownames_to_column(Visvader_0308_HER2_Navin_Visvader_NORM_panCK_Reference.meta.data, "Cells")
# In Seurat, "X" becomes the cell names while "cnv.ident" is specified above
Visvader_0308_HER2_Navin_Visvader_NORM_panCK_Reference.annotations <- Visvader_0308_HER2_Navin_Visvader_NORM_panCK_Reference.meta.data[, c("Cells", "cnv.ident")]
# Delete header
names(Visvader_0308_HER2_Navin_Visvader_NORM_panCK_Reference.annotations) <- NULL
# Save as table
write.table(Visvader_0308_HER2_Navin_Visvader_NORM_panCK_Reference.annotations,"/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_0308_HER2_Navin_Visvader_NORM_panCK_Reference.annotations.txt",sep="\t",row.names=FALSE)
# Remove metadata and annotations
rm(Visvader_0308_HER2_Navin_Visvader_NORM_panCK_Reference.meta.data)
rm(Visvader_0308_HER2_Navin_Visvader_NORM_panCK_Reference.annotations)

########## Generate Count Matrix
# Extract counts from the Seurat object
Visvader_0308_HER2_Navin_Visvader_NORM_panCK_Reference.count <- Visvader_0308_HER2_Navin_Visvader_NORM_panCK_Reference[["RNA"]]$counts
# Convert counts to a sparse matrix
Visvader_0308_HER2_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as(Visvader_0308_HER2_Navin_Visvader_NORM_panCK_Reference.count, 'sparseMatrix')
rm(Visvader_0308_HER2_Navin_Visvader_NORM_panCK_Reference.count)

# Convert to requisite data format
Visvader_0308_HER2_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as.data.frame(Visvader_0308_HER2_Navin_Visvader_NORM_panCK_Reference.count.matrix)
# Save table
write.csv(Visvader_0308_HER2_Navin_Visvader_NORM_panCK_Reference.count.matrix, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_Count_Matrix/Visvader_0308_HER2_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv")
# Remove Objects
rm(Visvader_0308_HER2_Navin_Visvader_NORM_panCK_Reference)
rm(Visvader_0308_HER2_Navin_Visvader_NORM_panCK_Reference.count.matrix)

########### Run inferCNV
# Note: check.names = FALSE prevents cell names from having dashes replaced with dots; essential to match annotations file
Visvader_0308_HER2_inferCNV <- CreateInfercnvObject(raw_counts_matrix=read.csv(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_Count_Matrix/Visvader_0308_HER2_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv", header = TRUE, row.names = 1,check.names = FALSE),
                                                    annotations_file=read.table(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_0308_HER2_Navin_Visvader_NORM_panCK_Reference.annotations.txt", header = FALSE, row.names = 1),
                                                    delim="\t",
                                                    gene_order_file="/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/hg38_gencode_v27.txt",
                                                    ref_group_names=c("Normal"))


# perform infercnv operations to reveal cnv signal
Visvader_0308_HER2_inferCNV = infercnv::run(Visvader_0308_HER2_inferCNV,
                                            cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                            out_dir="/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0308_HER2_panCK/",  # dir is auto-created for storing outputs
                                            cluster_by_groups=T,   # cluster
                                            denoise=T,
                                            HMM=T,
                                            num_ref_groups = 8, analysis_mode = "subclusters")

###########  Integrate data with Seurat object
# Read RDS and add inferCNV results to metadata
Visvader_0308_HER2_Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Visvader_0308_HER2_Navin_Visvader_NORM_panCK_Reference.rds") 
Visvader_0308_HER2_panCK_cnv <- infercnv::add_to_seurat(infercnv_output_path="/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0308_HER2_panCK/",
                                                        seurat_obj=Visvader_0308_HER2_Navin_Visvader_NORM_panCK_Reference, # optional
                                                        top_n=10)

# Collect Tumor cells from Seurat Object
Visvader_0308_HER2_panCK_cnv <- subset(x = Visvader_0308_HER2_panCK_cnv, subset = cnv.ident == "Tumor")
# Quantify chromosomes with alterations per cell
Visvader_0308_HER2_panCK_cnv_Pos_Chromosomes <- as.data.frame(rowSums(Visvader_0308_HER2_panCK_cnv@meta.data[,startsWith(names(Visvader_0308_HER2_panCK_cnv@meta.data),"has_cnv_")] > 0 ))
Visvader_0308_HER2_panCK_cnv_Neg_Chromosomes <- as.data.frame(rowSums(Visvader_0308_HER2_panCK_cnv@meta.data[,startsWith(names(Visvader_0308_HER2_panCK_cnv@meta.data),"has_cnv_")] == 0))
Visvader_0308_HER2_panCK_cnv_Chromosome_Quant <- as.data.frame(c(Visvader_0308_HER2_panCK_cnv_Pos_Chromosomes, Visvader_0308_HER2_panCK_cnv_Neg_Chromosomes))
colnames(Visvader_0308_HER2_panCK_cnv_Chromosome_Quant) <- c("Chromosomes CNV Pos", "Chromosomes CNV Neg")
write.csv(Visvader_0308_HER2_panCK_cnv_Chromosome_Quant, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0308_HER2_panCK/Visvader_0308_HER2_panCK_cnv_Chromosome_Quant.csv")

# Save CNV Metadata
Visvader_0308_HER2_panCK_cnv.meta.data <- as.data.frame(as.matrix(Visvader_0308_HER2_panCK_cnv@meta.data))
write.csv(Visvader_0308_HER2_panCK_cnv.meta.data, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0308_HER2_panCK/Visvader_0308_HER2_panCK_cnv.meta.data.csv")

# Subset Cells with inferred CNVs -------------------------------------------------------------------------------------------------------------
Visvader_0308_HER2_panCK_cnv <- subset(Visvader_0308_HER2_panCK_cnv, cells=rownames(Visvader_0308_HER2_panCK_cnv@meta.data)[rowSums(Visvader_0308_HER2_panCK_cnv@meta.data[,startsWith(names(Visvader_0308_HER2_panCK_cnv@meta.data),"has_cnv_")]) > 0 ])

# Save RDS with inferredCNVs
saveRDS(Visvader_0308_HER2_panCK_cnv, file = "/R/R_Visvader/Visvader_RDS/RDS_inferCNV/Visvader_0308_HER2_panCK_cnv.rds")


# Free memory
rm(Visvader_0308_HER2_inferCNV)
rm(Visvader_0308_HER2_Navin_Visvader_NORM_panCK_Reference)
rm(Visvader_0308_HER2_panCK_cnv)
rm(Visvader_0308_HER2_panCK_cnv_Pos_Chromosomes)
rm(Visvader_0308_HER2_panCK_cnv_Neg_Chromosomes)
rm(Visvader_0308_HER2_panCK_cnv_Chromosome_Quant)
rm(Visvader_0308_HER2_panCK_cnv.meta.data)
gc()


############################################
#### InferCNV: Visvader_0337_HER2_panCK ####
############################################

########## Merge with normal mammary glad, extract epithelium, and process via Seurat
# Load RDS
Visvader_0337_HER2_panCK <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0337_HER2_panCK.rds")
Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Navin_Visvader_NORM_panCK_Reference.rds")

# Set identities 
colnames(x = Visvader_0337_HER2_panCK[[]])
Idents(object = Visvader_0337_HER2_panCK) <- 'Tumor'
Visvader_0337_HER2_panCK[["cnv.ident"]] <- Idents(object = Visvader_0337_HER2_panCK)
levels(x = Visvader_0337_HER2_panCK)

colnames(x = Navin_Visvader_NORM_panCK_Reference[[]])
Idents(object = Navin_Visvader_NORM_panCK_Reference) <- 'Normal'
Navin_Visvader_NORM_panCK_Reference[["cnv.ident"]] <- Idents(object = Navin_Visvader_NORM_panCK_Reference)
levels(x = Navin_Visvader_NORM_panCK_Reference)

# Merge subsets into recombined RDS
Visvader_0337_HER2_Navin_Visvader_NORM_panCK_Reference <- merge(x = Visvader_0337_HER2_panCK, y = Navin_Visvader_NORM_panCK_Reference)
# Join Layers
Visvader_0337_HER2_Navin_Visvader_NORM_panCK_Reference <- JoinLayers(Visvader_0337_HER2_Navin_Visvader_NORM_panCK_Reference)

# Remove subsetted objects
rm(Visvader_0337_HER2_panCK)
rm(Navin_Visvader_NORM_panCK_Reference)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Visvader_0337_HER2_Navin_Visvader_NORM_panCK_Reference <- FindVariableFeatures(Visvader_0337_HER2_Navin_Visvader_NORM_panCK_Reference, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Visvader_0337_HER2_Navin_Visvader_NORM_panCK_Reference)
Visvader_0337_HER2_Navin_Visvader_NORM_panCK_Reference <- ScaleData(Visvader_0337_HER2_Navin_Visvader_NORM_panCK_Reference, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Visvader_0337_HER2_Navin_Visvader_NORM_panCK_Reference <- RunPCA(Visvader_0337_HER2_Navin_Visvader_NORM_panCK_Reference, features = VariableFeatures(object = Visvader_0337_HER2_Navin_Visvader_NORM_panCK_Reference))

# Examine and visualize PCA results a few different ways
print(Visvader_0337_HER2_Navin_Visvader_NORM_panCK_Reference[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Visvader_0337_HER2_Navin_Visvader_NORM_panCK_Reference, dims = 1:2, reduction = "pca")
DimPlot(Visvader_0337_HER2_Navin_Visvader_NORM_panCK_Reference, reduction = "pca")

Visvader_0337_HER2_Navin_Visvader_NORM_panCK_Reference <- FindNeighbors(Visvader_0337_HER2_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
Visvader_0337_HER2_Navin_Visvader_NORM_panCK_Reference <- FindClusters(Visvader_0337_HER2_Navin_Visvader_NORM_panCK_Reference, resolution = 0.1)

# Umap clustering
Visvader_0337_HER2_Navin_Visvader_NORM_panCK_Reference <- RunUMAP(Visvader_0337_HER2_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
DimPlot(Visvader_0337_HER2_Navin_Visvader_NORM_panCK_Reference, reduction = "umap")
DimPlot(Visvader_0337_HER2_Navin_Visvader_NORM_panCK_Reference, reduction = "umap", group.by = "orig.ident")

# Save RDS
saveRDS(Visvader_0337_HER2_Navin_Visvader_NORM_panCK_Reference, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Visvader_0337_HER2_Navin_Visvader_NORM_panCK_Reference.rds")

########## Prepare Cell Annotations
# Extract & Save Meta Data
Visvader_0337_HER2_Navin_Visvader_NORM_panCK_Reference.meta.data <- as.data.frame(as.matrix(Visvader_0337_HER2_Navin_Visvader_NORM_panCK_Reference@meta.data))
# Extract pertinent columns
# Make row names first column
Visvader_0337_HER2_Navin_Visvader_NORM_panCK_Reference.meta.data <- tibble::rownames_to_column(Visvader_0337_HER2_Navin_Visvader_NORM_panCK_Reference.meta.data, "Cells")
# In Seurat, "X" becomes the cell names while "cnv.ident" is specified above
Visvader_0337_HER2_Navin_Visvader_NORM_panCK_Reference.annotations <- Visvader_0337_HER2_Navin_Visvader_NORM_panCK_Reference.meta.data[, c("Cells", "cnv.ident")]
# Delete header
names(Visvader_0337_HER2_Navin_Visvader_NORM_panCK_Reference.annotations) <- NULL
# Save as table
write.table(Visvader_0337_HER2_Navin_Visvader_NORM_panCK_Reference.annotations,"/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_0337_HER2_Navin_Visvader_NORM_panCK_Reference.annotations.txt",sep="\t",row.names=FALSE)
# Remove metadata and annotations
rm(Visvader_0337_HER2_Navin_Visvader_NORM_panCK_Reference.meta.data)
rm(Visvader_0337_HER2_Navin_Visvader_NORM_panCK_Reference.annotations)

########## Generate Count Matrix
# Extract counts from the Seurat object
Visvader_0337_HER2_Navin_Visvader_NORM_panCK_Reference.count <- Visvader_0337_HER2_Navin_Visvader_NORM_panCK_Reference[["RNA"]]$counts
# Convert counts to a sparse matrix
Visvader_0337_HER2_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as(Visvader_0337_HER2_Navin_Visvader_NORM_panCK_Reference.count, 'sparseMatrix')
rm(Visvader_0337_HER2_Navin_Visvader_NORM_panCK_Reference.count)

# Convert to requisite data format
Visvader_0337_HER2_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as.data.frame(Visvader_0337_HER2_Navin_Visvader_NORM_panCK_Reference.count.matrix)
# Save table
write.csv(Visvader_0337_HER2_Navin_Visvader_NORM_panCK_Reference.count.matrix, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_Count_Matrix/Visvader_0337_HER2_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv")
# Remove Objects
rm(Visvader_0337_HER2_Navin_Visvader_NORM_panCK_Reference)
rm(Visvader_0337_HER2_Navin_Visvader_NORM_panCK_Reference.count.matrix)

########### Run inferCNV
# Note: check.names = FALSE prevents cell names from having dashes replaced with dots; essential to match annotations file
Visvader_0337_HER2_inferCNV <- CreateInfercnvObject(raw_counts_matrix=read.csv(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_Count_Matrix/Visvader_0337_HER2_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv", header = TRUE, row.names = 1,check.names = FALSE),
                                                    annotations_file=read.table(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_0337_HER2_Navin_Visvader_NORM_panCK_Reference.annotations.txt", header = FALSE, row.names = 1),
                                                    delim="\t",
                                                    gene_order_file="/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/hg38_gencode_v27.txt",
                                                    ref_group_names=c("Normal"))


# perform infercnv operations to reveal cnv signal
Visvader_0337_HER2_inferCNV = infercnv::run(Visvader_0337_HER2_inferCNV,
                                            cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                            out_dir="/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0337_HER2_panCK/",  # dir is auto-created for storing outputs
                                            cluster_by_groups=T,   # cluster
                                            denoise=T,
                                            HMM=T,
                                            num_ref_groups = 8, analysis_mode = "subclusters")

###########  Integrate data with Seurat object
# Read RDS and add inferCNV results to metadata
Visvader_0337_HER2_Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Visvader_0337_HER2_Navin_Visvader_NORM_panCK_Reference.rds") 
Visvader_0337_HER2_panCK_cnv <- infercnv::add_to_seurat(infercnv_output_path="/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0337_HER2_panCK/",
                                                        seurat_obj=Visvader_0337_HER2_Navin_Visvader_NORM_panCK_Reference, # optional
                                                        top_n=10)

# Collect Tumor cells from Seurat Object
Visvader_0337_HER2_panCK_cnv <- subset(x = Visvader_0337_HER2_panCK_cnv, subset = cnv.ident == "Tumor")
# Quantify chromosomes with alterations per cell
Visvader_0337_HER2_panCK_cnv_Pos_Chromosomes <- as.data.frame(rowSums(Visvader_0337_HER2_panCK_cnv@meta.data[,startsWith(names(Visvader_0337_HER2_panCK_cnv@meta.data),"has_cnv_")] > 0 ))
Visvader_0337_HER2_panCK_cnv_Neg_Chromosomes <- as.data.frame(rowSums(Visvader_0337_HER2_panCK_cnv@meta.data[,startsWith(names(Visvader_0337_HER2_panCK_cnv@meta.data),"has_cnv_")] == 0))
Visvader_0337_HER2_panCK_cnv_Chromosome_Quant <- as.data.frame(c(Visvader_0337_HER2_panCK_cnv_Pos_Chromosomes, Visvader_0337_HER2_panCK_cnv_Neg_Chromosomes))
colnames(Visvader_0337_HER2_panCK_cnv_Chromosome_Quant) <- c("Chromosomes CNV Pos", "Chromosomes CNV Neg")
write.csv(Visvader_0337_HER2_panCK_cnv_Chromosome_Quant, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0337_HER2_panCK/Visvader_0337_HER2_panCK_cnv_Chromosome_Quant.csv")

# Save CNV Metadata
Visvader_0337_HER2_panCK_cnv.meta.data <- as.data.frame(as.matrix(Visvader_0337_HER2_panCK_cnv@meta.data))
write.csv(Visvader_0337_HER2_panCK_cnv.meta.data, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0337_HER2_panCK/Visvader_0337_HER2_panCK_cnv.meta.data.csv")

# Subset Cells with inferred CNVs -------------------------------------------------------------------------------------------------------------
Visvader_0337_HER2_panCK_cnv <- subset(Visvader_0337_HER2_panCK_cnv, cells=rownames(Visvader_0337_HER2_panCK_cnv@meta.data)[rowSums(Visvader_0337_HER2_panCK_cnv@meta.data[,startsWith(names(Visvader_0337_HER2_panCK_cnv@meta.data),"has_cnv_")]) > 0 ])

# Save RDS with inferredCNVs
saveRDS(Visvader_0337_HER2_panCK_cnv, file = "/R/R_Visvader/Visvader_RDS/RDS_inferCNV/Visvader_0337_HER2_panCK_cnv.rds")


# Free memory
rm(Visvader_0337_HER2_inferCNV)
rm(Visvader_0337_HER2_Navin_Visvader_NORM_panCK_Reference)
rm(Visvader_0337_HER2_panCK_cnv)
rm(Visvader_0337_HER2_panCK_cnv_Pos_Chromosomes)
rm(Visvader_0337_HER2_panCK_cnv_Neg_Chromosomes)
rm(Visvader_0337_HER2_panCK_cnv_Chromosome_Quant)
rm(Visvader_0337_HER2_panCK_cnv.meta.data)
gc()




###################################################################################################
#############################     InferCNV: TNBC Tumor Specimens      #############################     
###################################################################################################

###########################################
#### InferCNV: Visvader_106_TNBC_panCK ####
###########################################

########## Merge with normal mammary glad, extract epithelium, and process via Seurat
# Load RDS
Visvader_106_TNBC_panCK <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_106_TNBC_panCK.rds")
Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Navin_Visvader_NORM_panCK_Reference.rds")

# Set identities 
colnames(x = Visvader_106_TNBC_panCK[[]])
Idents(object = Visvader_106_TNBC_panCK) <- 'Tumor'
Visvader_106_TNBC_panCK[["cnv.ident"]] <- Idents(object = Visvader_106_TNBC_panCK)
levels(x = Visvader_106_TNBC_panCK)

colnames(x = Navin_Visvader_NORM_panCK_Reference[[]])
Idents(object = Navin_Visvader_NORM_panCK_Reference) <- 'Normal'
Navin_Visvader_NORM_panCK_Reference[["cnv.ident"]] <- Idents(object = Navin_Visvader_NORM_panCK_Reference)
levels(x = Navin_Visvader_NORM_panCK_Reference)

# Merge subsets into recombined RDS
Visvader_106_TNBC_Navin_Visvader_NORM_panCK_Reference <- merge(x = Visvader_106_TNBC_panCK, y = Navin_Visvader_NORM_panCK_Reference)
# Join Layers
Visvader_106_TNBC_Navin_Visvader_NORM_panCK_Reference <- JoinLayers(Visvader_106_TNBC_Navin_Visvader_NORM_panCK_Reference)

# Remove subsetted objects
rm(Visvader_106_TNBC_panCK)
rm(Navin_Visvader_NORM_panCK_Reference)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Visvader_106_TNBC_Navin_Visvader_NORM_panCK_Reference <- FindVariableFeatures(Visvader_106_TNBC_Navin_Visvader_NORM_panCK_Reference, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Visvader_106_TNBC_Navin_Visvader_NORM_panCK_Reference)
Visvader_106_TNBC_Navin_Visvader_NORM_panCK_Reference <- ScaleData(Visvader_106_TNBC_Navin_Visvader_NORM_panCK_Reference, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Visvader_106_TNBC_Navin_Visvader_NORM_panCK_Reference <- RunPCA(Visvader_106_TNBC_Navin_Visvader_NORM_panCK_Reference, features = VariableFeatures(object = Visvader_106_TNBC_Navin_Visvader_NORM_panCK_Reference))

# Examine and visualize PCA results a few different ways
print(Visvader_106_TNBC_Navin_Visvader_NORM_panCK_Reference[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Visvader_106_TNBC_Navin_Visvader_NORM_panCK_Reference, dims = 1:2, reduction = "pca")
DimPlot(Visvader_106_TNBC_Navin_Visvader_NORM_panCK_Reference, reduction = "pca")

Visvader_106_TNBC_Navin_Visvader_NORM_panCK_Reference <- FindNeighbors(Visvader_106_TNBC_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
Visvader_106_TNBC_Navin_Visvader_NORM_panCK_Reference <- FindClusters(Visvader_106_TNBC_Navin_Visvader_NORM_panCK_Reference, resolution = 0.1)

# Umap clustering
Visvader_106_TNBC_Navin_Visvader_NORM_panCK_Reference <- RunUMAP(Visvader_106_TNBC_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
DimPlot(Visvader_106_TNBC_Navin_Visvader_NORM_panCK_Reference, reduction = "umap")
DimPlot(Visvader_106_TNBC_Navin_Visvader_NORM_panCK_Reference, reduction = "umap", group.by = "orig.ident")

# Save RDS
saveRDS(Visvader_106_TNBC_Navin_Visvader_NORM_panCK_Reference, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Visvader_106_TNBC_Navin_Visvader_NORM_panCK_Reference.rds")

########## Prepare Cell Annotations
# Extract & Save Meta Data
Visvader_106_TNBC_Navin_Visvader_NORM_panCK_Reference.meta.data <- as.data.frame(as.matrix(Visvader_106_TNBC_Navin_Visvader_NORM_panCK_Reference@meta.data))
# Extract pertinent columns
# Make row names first column
Visvader_106_TNBC_Navin_Visvader_NORM_panCK_Reference.meta.data <- tibble::rownames_to_column(Visvader_106_TNBC_Navin_Visvader_NORM_panCK_Reference.meta.data, "Cells")
# In Seurat, "X" becomes the cell names while "cnv.ident" is specified above
Visvader_106_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations <- Visvader_106_TNBC_Navin_Visvader_NORM_panCK_Reference.meta.data[, c("Cells", "cnv.ident")]
# Delete header
names(Visvader_106_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations) <- NULL
# Save as table
write.table(Visvader_106_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations,"/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_106_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations.txt",sep="\t",row.names=FALSE)
# Remove metadata and annotations
rm(Visvader_106_TNBC_Navin_Visvader_NORM_panCK_Reference.meta.data)
rm(Visvader_106_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations)

########## Generate Count Matrix
# Extract counts from the Seurat object
Visvader_106_TNBC_Navin_Visvader_NORM_panCK_Reference.count <- Visvader_106_TNBC_Navin_Visvader_NORM_panCK_Reference[["RNA"]]$counts
# Convert counts to a sparse matrix
Visvader_106_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as(Visvader_106_TNBC_Navin_Visvader_NORM_panCK_Reference.count, 'sparseMatrix')
rm(Visvader_106_TNBC_Navin_Visvader_NORM_panCK_Reference.count)

# Convert to requisite data format
Visvader_106_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as.data.frame(Visvader_106_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix)
# Save table
write.csv(Visvader_106_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_Count_Matrix/Visvader_106_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv")
# Remove Objects
rm(Visvader_106_TNBC_Navin_Visvader_NORM_panCK_Reference)
rm(Visvader_106_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix)

########### Run inferCNV
# Note: check.names = FALSE prevents cell names from having dashes replaced with dots; essential to match annotations file
Visvader_106_TNBC_inferCNV <- CreateInfercnvObject(raw_counts_matrix=read.csv(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_Count_Matrix/Visvader_106_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv", header = TRUE, row.names = 1,check.names = FALSE),
                                                   annotations_file=read.table(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_106_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations.txt", header = FALSE, row.names = 1),
                                                   delim="\t",
                                                   gene_order_file="/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/hg38_gencode_v27.txt",
                                                   ref_group_names=c("Normal"))


# perform infercnv operations to reveal cnv signal
Visvader_106_TNBC_inferCNV = infercnv::run(Visvader_106_TNBC_inferCNV,
                                           cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                           out_dir="/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_106_TNBC_panCK/",  # dir is auto-created for storing outputs
                                           cluster_by_groups=T,   # cluster
                                           denoise=T,
                                           HMM=T,
                                           num_ref_groups = 8, analysis_mode = "subclusters")

###########  Integrate data with Seurat object
# Read RDS and add inferCNV results to metadata
Visvader_106_TNBC_Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Visvader_106_TNBC_Navin_Visvader_NORM_panCK_Reference.rds") 
Visvader_106_TNBC_panCK_cnv <- infercnv::add_to_seurat(infercnv_output_path="/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_106_TNBC_panCK/",
                                                       seurat_obj=Visvader_106_TNBC_Navin_Visvader_NORM_panCK_Reference, # optional
                                                       top_n=10)

# Collect Tumor cells from Seurat Object
Visvader_106_TNBC_panCK_cnv <- subset(x = Visvader_106_TNBC_panCK_cnv, subset = cnv.ident == "Tumor")
# Quantify chromosomes with alterations per cell
Visvader_106_TNBC_panCK_cnv_Pos_Chromosomes <- as.data.frame(rowSums(Visvader_106_TNBC_panCK_cnv@meta.data[,startsWith(names(Visvader_106_TNBC_panCK_cnv@meta.data),"has_cnv_")] > 0 ))
Visvader_106_TNBC_panCK_cnv_Neg_Chromosomes <- as.data.frame(rowSums(Visvader_106_TNBC_panCK_cnv@meta.data[,startsWith(names(Visvader_106_TNBC_panCK_cnv@meta.data),"has_cnv_")] == 0))
Visvader_106_TNBC_panCK_cnv_Chromosome_Quant <- as.data.frame(c(Visvader_106_TNBC_panCK_cnv_Pos_Chromosomes, Visvader_106_TNBC_panCK_cnv_Neg_Chromosomes))
colnames(Visvader_106_TNBC_panCK_cnv_Chromosome_Quant) <- c("Chromosomes CNV Pos", "Chromosomes CNV Neg")
write.csv(Visvader_106_TNBC_panCK_cnv_Chromosome_Quant, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_106_TNBC_panCK/Visvader_106_TNBC_panCK_cnv_Chromosome_Quant.csv")

# Save CNV Metadata
Visvader_106_TNBC_panCK_cnv.meta.data <- as.data.frame(as.matrix(Visvader_106_TNBC_panCK_cnv@meta.data))
write.csv(Visvader_106_TNBC_panCK_cnv.meta.data, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_106_TNBC_panCK/Visvader_106_TNBC_panCK_cnv.meta.data.csv")

# Subset Cells with inferred CNVs -------------------------------------------------------------------------------------------------------------
Visvader_106_TNBC_panCK_cnv <- subset(Visvader_106_TNBC_panCK_cnv, cells=rownames(Visvader_106_TNBC_panCK_cnv@meta.data)[rowSums(Visvader_106_TNBC_panCK_cnv@meta.data[,startsWith(names(Visvader_106_TNBC_panCK_cnv@meta.data),"has_cnv_")]) > 0 ])

# Save RDS with inferredCNVs
saveRDS(Visvader_106_TNBC_panCK_cnv, file = "/R/R_Visvader/Visvader_RDS/RDS_inferCNV/Visvader_106_TNBC_panCK_cnv.rds")


# Free memory
rm(Visvader_106_TNBC_inferCNV)
rm(Visvader_106_TNBC_Navin_Visvader_NORM_panCK_Reference)
rm(Visvader_106_TNBC_panCK_cnv)
rm(Visvader_106_TNBC_panCK_cnv_Pos_Chromosomes)
rm(Visvader_106_TNBC_panCK_cnv_Neg_Chromosomes)
rm(Visvader_106_TNBC_panCK_cnv_Chromosome_Quant)
rm(Visvader_106_TNBC_panCK_cnv.meta.data)
gc()



###########################################
#### InferCNV: Visvader_114_TNBC_panCK ####
###########################################

########## Merge with normal mammary glad, extract epithelium, and process via Seurat
# Load RDS
Visvader_114_TNBC_panCK <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_114_TNBC_panCK.rds")
Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Navin_Visvader_NORM_panCK_Reference.rds")

# Set identities 
colnames(x = Visvader_114_TNBC_panCK[[]])
Idents(object = Visvader_114_TNBC_panCK) <- 'Tumor'
Visvader_114_TNBC_panCK[["cnv.ident"]] <- Idents(object = Visvader_114_TNBC_panCK)
levels(x = Visvader_114_TNBC_panCK)

colnames(x = Navin_Visvader_NORM_panCK_Reference[[]])
Idents(object = Navin_Visvader_NORM_panCK_Reference) <- 'Normal'
Navin_Visvader_NORM_panCK_Reference[["cnv.ident"]] <- Idents(object = Navin_Visvader_NORM_panCK_Reference)
levels(x = Navin_Visvader_NORM_panCK_Reference)

# Merge subsets into recombined RDS
Visvader_114_TNBC_Navin_Visvader_NORM_panCK_Reference <- merge(x = Visvader_114_TNBC_panCK, y = Navin_Visvader_NORM_panCK_Reference)
# Join Layers
Visvader_114_TNBC_Navin_Visvader_NORM_panCK_Reference <- JoinLayers(Visvader_114_TNBC_Navin_Visvader_NORM_panCK_Reference)

# Remove subsetted objects
rm(Visvader_114_TNBC_panCK)
rm(Navin_Visvader_NORM_panCK_Reference)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Visvader_114_TNBC_Navin_Visvader_NORM_panCK_Reference <- FindVariableFeatures(Visvader_114_TNBC_Navin_Visvader_NORM_panCK_Reference, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Visvader_114_TNBC_Navin_Visvader_NORM_panCK_Reference)
Visvader_114_TNBC_Navin_Visvader_NORM_panCK_Reference <- ScaleData(Visvader_114_TNBC_Navin_Visvader_NORM_panCK_Reference, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Visvader_114_TNBC_Navin_Visvader_NORM_panCK_Reference <- RunPCA(Visvader_114_TNBC_Navin_Visvader_NORM_panCK_Reference, features = VariableFeatures(object = Visvader_114_TNBC_Navin_Visvader_NORM_panCK_Reference))

# Examine and visualize PCA results a few different ways
print(Visvader_114_TNBC_Navin_Visvader_NORM_panCK_Reference[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Visvader_114_TNBC_Navin_Visvader_NORM_panCK_Reference, dims = 1:2, reduction = "pca")
DimPlot(Visvader_114_TNBC_Navin_Visvader_NORM_panCK_Reference, reduction = "pca")

Visvader_114_TNBC_Navin_Visvader_NORM_panCK_Reference <- FindNeighbors(Visvader_114_TNBC_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
Visvader_114_TNBC_Navin_Visvader_NORM_panCK_Reference <- FindClusters(Visvader_114_TNBC_Navin_Visvader_NORM_panCK_Reference, resolution = 0.1)

# Umap clustering
Visvader_114_TNBC_Navin_Visvader_NORM_panCK_Reference <- RunUMAP(Visvader_114_TNBC_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
DimPlot(Visvader_114_TNBC_Navin_Visvader_NORM_panCK_Reference, reduction = "umap")
DimPlot(Visvader_114_TNBC_Navin_Visvader_NORM_panCK_Reference, reduction = "umap", group.by = "orig.ident")

# Save RDS
saveRDS(Visvader_114_TNBC_Navin_Visvader_NORM_panCK_Reference, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Visvader_114_TNBC_Navin_Visvader_NORM_panCK_Reference.rds")

########## Prepare Cell Annotations
# Extract & Save Meta Data
Visvader_114_TNBC_Navin_Visvader_NORM_panCK_Reference.meta.data <- as.data.frame(as.matrix(Visvader_114_TNBC_Navin_Visvader_NORM_panCK_Reference@meta.data))
# Extract pertinent columns
# Make row names first column
Visvader_114_TNBC_Navin_Visvader_NORM_panCK_Reference.meta.data <- tibble::rownames_to_column(Visvader_114_TNBC_Navin_Visvader_NORM_panCK_Reference.meta.data, "Cells")
# In Seurat, "X" becomes the cell names while "cnv.ident" is specified above
Visvader_114_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations <- Visvader_114_TNBC_Navin_Visvader_NORM_panCK_Reference.meta.data[, c("Cells", "cnv.ident")]
# Delete header
names(Visvader_114_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations) <- NULL
# Save as table
write.table(Visvader_114_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations,"/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_114_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations.txt",sep="\t",row.names=FALSE)
# Remove metadata and annotations
rm(Visvader_114_TNBC_Navin_Visvader_NORM_panCK_Reference.meta.data)
rm(Visvader_114_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations)

########## Generate Count Matrix
# Extract counts from the Seurat object
Visvader_114_TNBC_Navin_Visvader_NORM_panCK_Reference.count <- Visvader_114_TNBC_Navin_Visvader_NORM_panCK_Reference[["RNA"]]$counts
# Convert counts to a sparse matrix
Visvader_114_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as(Visvader_114_TNBC_Navin_Visvader_NORM_panCK_Reference.count, 'sparseMatrix')
rm(Visvader_114_TNBC_Navin_Visvader_NORM_panCK_Reference.count)

# Convert to requisite data format
Visvader_114_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as.data.frame(Visvader_114_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix)
# Save table
write.csv(Visvader_114_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_Count_Matrix/Visvader_114_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv")
# Remove Objects
rm(Visvader_114_TNBC_Navin_Visvader_NORM_panCK_Reference)
rm(Visvader_114_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix)

########### Run inferCNV
# Note: check.names = FALSE prevents cell names from having dashes replaced with dots; essential to match annotations file
Visvader_114_TNBC_inferCNV <- CreateInfercnvObject(raw_counts_matrix=read.csv(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_Count_Matrix/Visvader_114_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv", header = TRUE, row.names = 1,check.names = FALSE),
                                                   annotations_file=read.table(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_114_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations.txt", header = FALSE, row.names = 1),
                                                   delim="\t",
                                                   gene_order_file="/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/hg38_gencode_v27.txt",
                                                   ref_group_names=c("Normal"))


# perform infercnv operations to reveal cnv signal
Visvader_114_TNBC_inferCNV = infercnv::run(Visvader_114_TNBC_inferCNV,
                                           cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                           out_dir="/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_114_TNBC_panCK/",  # dir is auto-created for storing outputs
                                           cluster_by_groups=T,   # cluster
                                           denoise=T,
                                           HMM=T,
                                           num_ref_groups = 8, analysis_mode = "subclusters")

###########  Integrate data with Seurat object
# Read RDS and add inferCNV results to metadata
Visvader_114_TNBC_Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Visvader_114_TNBC_Navin_Visvader_NORM_panCK_Reference.rds") 
Visvader_114_TNBC_panCK_cnv <- infercnv::add_to_seurat(infercnv_output_path="/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_114_TNBC_panCK/",
                                                       seurat_obj=Visvader_114_TNBC_Navin_Visvader_NORM_panCK_Reference, # optional
                                                       top_n=10)

# Collect Tumor cells from Seurat Object
Visvader_114_TNBC_panCK_cnv <- subset(x = Visvader_114_TNBC_panCK_cnv, subset = cnv.ident == "Tumor")
# Quantify chromosomes with alterations per cell
Visvader_114_TNBC_panCK_cnv_Pos_Chromosomes <- as.data.frame(rowSums(Visvader_114_TNBC_panCK_cnv@meta.data[,startsWith(names(Visvader_114_TNBC_panCK_cnv@meta.data),"has_cnv_")] > 0 ))
Visvader_114_TNBC_panCK_cnv_Neg_Chromosomes <- as.data.frame(rowSums(Visvader_114_TNBC_panCK_cnv@meta.data[,startsWith(names(Visvader_114_TNBC_panCK_cnv@meta.data),"has_cnv_")] == 0))
Visvader_114_TNBC_panCK_cnv_Chromosome_Quant <- as.data.frame(c(Visvader_114_TNBC_panCK_cnv_Pos_Chromosomes, Visvader_114_TNBC_panCK_cnv_Neg_Chromosomes))
colnames(Visvader_114_TNBC_panCK_cnv_Chromosome_Quant) <- c("Chromosomes CNV Pos", "Chromosomes CNV Neg")
write.csv(Visvader_114_TNBC_panCK_cnv_Chromosome_Quant, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_114_TNBC_panCK/Visvader_114_TNBC_panCK_cnv_Chromosome_Quant.csv")

# Save CNV Metadata
Visvader_114_TNBC_panCK_cnv.meta.data <- as.data.frame(as.matrix(Visvader_114_TNBC_panCK_cnv@meta.data))
write.csv(Visvader_114_TNBC_panCK_cnv.meta.data, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_114_TNBC_panCK/Visvader_114_TNBC_panCK_cnv.meta.data.csv")

# Subset Cells with inferred CNVs -------------------------------------------------------------------------------------------------------------
Visvader_114_TNBC_panCK_cnv <- subset(Visvader_114_TNBC_panCK_cnv, cells=rownames(Visvader_114_TNBC_panCK_cnv@meta.data)[rowSums(Visvader_114_TNBC_panCK_cnv@meta.data[,startsWith(names(Visvader_114_TNBC_panCK_cnv@meta.data),"has_cnv_")]) > 0 ])

# Save RDS with inferredCNVs
saveRDS(Visvader_114_TNBC_panCK_cnv, file = "/R/R_Visvader/Visvader_RDS/RDS_inferCNV/Visvader_114_TNBC_panCK_cnv.rds")


# Free memory
rm(Visvader_114_TNBC_inferCNV)
rm(Visvader_114_TNBC_Navin_Visvader_NORM_panCK_Reference)
rm(Visvader_114_TNBC_panCK_cnv)
rm(Visvader_114_TNBC_panCK_cnv_Pos_Chromosomes)
rm(Visvader_114_TNBC_panCK_cnv_Neg_Chromosomes)
rm(Visvader_114_TNBC_panCK_cnv_Chromosome_Quant)
rm(Visvader_114_TNBC_panCK_cnv.meta.data)
gc()



###########################################
#### InferCNV: Visvader_126_TNBC_panCK ####
###########################################

########## Merge with normal mammary glad, extract epithelium, and process via Seurat
# Load RDS
Visvader_126_TNBC_panCK <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_126_TNBC_panCK.rds")
Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Navin_Visvader_NORM_panCK_Reference.rds")

# Set identities 
colnames(x = Visvader_126_TNBC_panCK[[]])
Idents(object = Visvader_126_TNBC_panCK) <- 'Tumor'
Visvader_126_TNBC_panCK[["cnv.ident"]] <- Idents(object = Visvader_126_TNBC_panCK)
levels(x = Visvader_126_TNBC_panCK)

colnames(x = Navin_Visvader_NORM_panCK_Reference[[]])
Idents(object = Navin_Visvader_NORM_panCK_Reference) <- 'Normal'
Navin_Visvader_NORM_panCK_Reference[["cnv.ident"]] <- Idents(object = Navin_Visvader_NORM_panCK_Reference)
levels(x = Navin_Visvader_NORM_panCK_Reference)

# Merge subsets into recombined RDS
Visvader_126_TNBC_Navin_Visvader_NORM_panCK_Reference <- merge(x = Visvader_126_TNBC_panCK, y = Navin_Visvader_NORM_panCK_Reference)
# Join Layers
Visvader_126_TNBC_Navin_Visvader_NORM_panCK_Reference <- JoinLayers(Visvader_126_TNBC_Navin_Visvader_NORM_panCK_Reference)

# Remove subsetted objects
rm(Visvader_126_TNBC_panCK)
rm(Navin_Visvader_NORM_panCK_Reference)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Visvader_126_TNBC_Navin_Visvader_NORM_panCK_Reference <- FindVariableFeatures(Visvader_126_TNBC_Navin_Visvader_NORM_panCK_Reference, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Visvader_126_TNBC_Navin_Visvader_NORM_panCK_Reference)
Visvader_126_TNBC_Navin_Visvader_NORM_panCK_Reference <- ScaleData(Visvader_126_TNBC_Navin_Visvader_NORM_panCK_Reference, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Visvader_126_TNBC_Navin_Visvader_NORM_panCK_Reference <- RunPCA(Visvader_126_TNBC_Navin_Visvader_NORM_panCK_Reference, features = VariableFeatures(object = Visvader_126_TNBC_Navin_Visvader_NORM_panCK_Reference))

# Examine and visualize PCA results a few different ways
print(Visvader_126_TNBC_Navin_Visvader_NORM_panCK_Reference[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Visvader_126_TNBC_Navin_Visvader_NORM_panCK_Reference, dims = 1:2, reduction = "pca")
DimPlot(Visvader_126_TNBC_Navin_Visvader_NORM_panCK_Reference, reduction = "pca")

Visvader_126_TNBC_Navin_Visvader_NORM_panCK_Reference <- FindNeighbors(Visvader_126_TNBC_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
Visvader_126_TNBC_Navin_Visvader_NORM_panCK_Reference <- FindClusters(Visvader_126_TNBC_Navin_Visvader_NORM_panCK_Reference, resolution = 0.1)

# Umap clustering
Visvader_126_TNBC_Navin_Visvader_NORM_panCK_Reference <- RunUMAP(Visvader_126_TNBC_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
DimPlot(Visvader_126_TNBC_Navin_Visvader_NORM_panCK_Reference, reduction = "umap")
DimPlot(Visvader_126_TNBC_Navin_Visvader_NORM_panCK_Reference, reduction = "umap", group.by = "orig.ident")

# Save RDS
saveRDS(Visvader_126_TNBC_Navin_Visvader_NORM_panCK_Reference, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Visvader_126_TNBC_Navin_Visvader_NORM_panCK_Reference.rds")

########## Prepare Cell Annotations
# Extract & Save Meta Data
Visvader_126_TNBC_Navin_Visvader_NORM_panCK_Reference.meta.data <- as.data.frame(as.matrix(Visvader_126_TNBC_Navin_Visvader_NORM_panCK_Reference@meta.data))
# Extract pertinent columns
# Make row names first column
Visvader_126_TNBC_Navin_Visvader_NORM_panCK_Reference.meta.data <- tibble::rownames_to_column(Visvader_126_TNBC_Navin_Visvader_NORM_panCK_Reference.meta.data, "Cells")
# In Seurat, "X" becomes the cell names while "cnv.ident" is specified above
Visvader_126_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations <- Visvader_126_TNBC_Navin_Visvader_NORM_panCK_Reference.meta.data[, c("Cells", "cnv.ident")]
# Delete header
names(Visvader_126_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations) <- NULL
# Save as table
write.table(Visvader_126_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations,"/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_126_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations.txt",sep="\t",row.names=FALSE)
# Remove metadata and annotations
rm(Visvader_126_TNBC_Navin_Visvader_NORM_panCK_Reference.meta.data)
rm(Visvader_126_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations)

########## Generate Count Matrix
# Extract counts from the Seurat object
Visvader_126_TNBC_Navin_Visvader_NORM_panCK_Reference.count <- Visvader_126_TNBC_Navin_Visvader_NORM_panCK_Reference[["RNA"]]$counts
# Convert counts to a sparse matrix
Visvader_126_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as(Visvader_126_TNBC_Navin_Visvader_NORM_panCK_Reference.count, 'sparseMatrix')
rm(Visvader_126_TNBC_Navin_Visvader_NORM_panCK_Reference.count)

# Convert to requisite data format
Visvader_126_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as.data.frame(Visvader_126_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix)
# Save table
write.csv(Visvader_126_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_Count_Matrix/Visvader_126_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv")
# Remove Objects
rm(Visvader_126_TNBC_Navin_Visvader_NORM_panCK_Reference)
rm(Visvader_126_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix)

########### Run inferCNV
# Note: check.names = FALSE prevents cell names from having dashes replaced with dots; essential to match annotations file
Visvader_126_TNBC_inferCNV <- CreateInfercnvObject(raw_counts_matrix=read.csv(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_Count_Matrix/Visvader_126_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv", header = TRUE, row.names = 1,check.names = FALSE),
                                                   annotations_file=read.table(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_126_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations.txt", header = FALSE, row.names = 1),
                                                   delim="\t",
                                                   gene_order_file="/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/hg38_gencode_v27.txt",
                                                   ref_group_names=c("Normal"))


# perform infercnv operations to reveal cnv signal
Visvader_126_TNBC_inferCNV = infercnv::run(Visvader_126_TNBC_inferCNV,
                                           cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                           out_dir="/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_126_TNBC_panCK/",  # dir is auto-created for storing outputs
                                           cluster_by_groups=T,   # cluster
                                           denoise=T,
                                           HMM=T,
                                           num_ref_groups = 8, analysis_mode = "subclusters")

###########  Integrate data with Seurat object
# Read RDS and add inferCNV results to metadata
Visvader_126_TNBC_Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Visvader_126_TNBC_Navin_Visvader_NORM_panCK_Reference.rds") 
Visvader_126_TNBC_panCK_cnv <- infercnv::add_to_seurat(infercnv_output_path="/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_126_TNBC_panCK/",
                                                       seurat_obj=Visvader_126_TNBC_Navin_Visvader_NORM_panCK_Reference, # optional
                                                       top_n=10)

# Collect Tumor cells from Seurat Object
Visvader_126_TNBC_panCK_cnv <- subset(x = Visvader_126_TNBC_panCK_cnv, subset = cnv.ident == "Tumor")
# Quantify chromosomes with alterations per cell
Visvader_126_TNBC_panCK_cnv_Pos_Chromosomes <- as.data.frame(rowSums(Visvader_126_TNBC_panCK_cnv@meta.data[,startsWith(names(Visvader_126_TNBC_panCK_cnv@meta.data),"has_cnv_")] > 0 ))
Visvader_126_TNBC_panCK_cnv_Neg_Chromosomes <- as.data.frame(rowSums(Visvader_126_TNBC_panCK_cnv@meta.data[,startsWith(names(Visvader_126_TNBC_panCK_cnv@meta.data),"has_cnv_")] == 0))
Visvader_126_TNBC_panCK_cnv_Chromosome_Quant <- as.data.frame(c(Visvader_126_TNBC_panCK_cnv_Pos_Chromosomes, Visvader_126_TNBC_panCK_cnv_Neg_Chromosomes))
colnames(Visvader_126_TNBC_panCK_cnv_Chromosome_Quant) <- c("Chromosomes CNV Pos", "Chromosomes CNV Neg")
write.csv(Visvader_126_TNBC_panCK_cnv_Chromosome_Quant, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_126_TNBC_panCK/Visvader_126_TNBC_panCK_cnv_Chromosome_Quant.csv")

# Save CNV Metadata
Visvader_126_TNBC_panCK_cnv.meta.data <- as.data.frame(as.matrix(Visvader_126_TNBC_panCK_cnv@meta.data))
write.csv(Visvader_126_TNBC_panCK_cnv.meta.data, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_126_TNBC_panCK/Visvader_126_TNBC_panCK_cnv.meta.data.csv")

# Subset Cells with inferred CNVs -------------------------------------------------------------------------------------------------------------
Visvader_126_TNBC_panCK_cnv <- subset(Visvader_126_TNBC_panCK_cnv, cells=rownames(Visvader_126_TNBC_panCK_cnv@meta.data)[rowSums(Visvader_126_TNBC_panCK_cnv@meta.data[,startsWith(names(Visvader_126_TNBC_panCK_cnv@meta.data),"has_cnv_")]) > 0 ])

# Save RDS with inferredCNVs
saveRDS(Visvader_126_TNBC_panCK_cnv, file = "/R/R_Visvader/Visvader_RDS/RDS_inferCNV/Visvader_126_TNBC_panCK_cnv.rds")


# Free memory
rm(Visvader_126_TNBC_inferCNV)
rm(Visvader_126_TNBC_Navin_Visvader_NORM_panCK_Reference)
rm(Visvader_126_TNBC_panCK_cnv)
rm(Visvader_126_TNBC_panCK_cnv_Pos_Chromosomes)
rm(Visvader_126_TNBC_panCK_cnv_Neg_Chromosomes)
rm(Visvader_126_TNBC_panCK_cnv_Chromosome_Quant)
rm(Visvader_126_TNBC_panCK_cnv.meta.data)
gc()



############################################
#### InferCNV: Visvader_0131_TNBC_panCK ####
############################################

########## Merge with normal mammary glad, extract epithelium, and process via Seurat
# Load RDS
Visvader_0131_TNBC_panCK <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0131_TNBC_panCK.rds")
Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Navin_Visvader_NORM_panCK_Reference.rds")

# Set identities 
colnames(x = Visvader_0131_TNBC_panCK[[]])
Idents(object = Visvader_0131_TNBC_panCK) <- 'Tumor'
Visvader_0131_TNBC_panCK[["cnv.ident"]] <- Idents(object = Visvader_0131_TNBC_panCK)
levels(x = Visvader_0131_TNBC_panCK)

colnames(x = Navin_Visvader_NORM_panCK_Reference[[]])
Idents(object = Navin_Visvader_NORM_panCK_Reference) <- 'Normal'
Navin_Visvader_NORM_panCK_Reference[["cnv.ident"]] <- Idents(object = Navin_Visvader_NORM_panCK_Reference)
levels(x = Navin_Visvader_NORM_panCK_Reference)

# Merge subsets into recombined RDS
Visvader_0131_TNBC_Navin_Visvader_NORM_panCK_Reference <- merge(x = Visvader_0131_TNBC_panCK, y = Navin_Visvader_NORM_panCK_Reference)
# Join Layers
Visvader_0131_TNBC_Navin_Visvader_NORM_panCK_Reference <- JoinLayers(Visvader_0131_TNBC_Navin_Visvader_NORM_panCK_Reference)

# Remove subsetted objects
rm(Visvader_0131_TNBC_panCK)
rm(Navin_Visvader_NORM_panCK_Reference)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Visvader_0131_TNBC_Navin_Visvader_NORM_panCK_Reference <- FindVariableFeatures(Visvader_0131_TNBC_Navin_Visvader_NORM_panCK_Reference, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Visvader_0131_TNBC_Navin_Visvader_NORM_panCK_Reference)
Visvader_0131_TNBC_Navin_Visvader_NORM_panCK_Reference <- ScaleData(Visvader_0131_TNBC_Navin_Visvader_NORM_panCK_Reference, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Visvader_0131_TNBC_Navin_Visvader_NORM_panCK_Reference <- RunPCA(Visvader_0131_TNBC_Navin_Visvader_NORM_panCK_Reference, features = VariableFeatures(object = Visvader_0131_TNBC_Navin_Visvader_NORM_panCK_Reference))

# Examine and visualize PCA results a few different ways
print(Visvader_0131_TNBC_Navin_Visvader_NORM_panCK_Reference[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Visvader_0131_TNBC_Navin_Visvader_NORM_panCK_Reference, dims = 1:2, reduction = "pca")
DimPlot(Visvader_0131_TNBC_Navin_Visvader_NORM_panCK_Reference, reduction = "pca")

Visvader_0131_TNBC_Navin_Visvader_NORM_panCK_Reference <- FindNeighbors(Visvader_0131_TNBC_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
Visvader_0131_TNBC_Navin_Visvader_NORM_panCK_Reference <- FindClusters(Visvader_0131_TNBC_Navin_Visvader_NORM_panCK_Reference, resolution = 0.1)

# Umap clustering
Visvader_0131_TNBC_Navin_Visvader_NORM_panCK_Reference <- RunUMAP(Visvader_0131_TNBC_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
DimPlot(Visvader_0131_TNBC_Navin_Visvader_NORM_panCK_Reference, reduction = "umap")
DimPlot(Visvader_0131_TNBC_Navin_Visvader_NORM_panCK_Reference, reduction = "umap", group.by = "orig.ident")

# Save RDS
saveRDS(Visvader_0131_TNBC_Navin_Visvader_NORM_panCK_Reference, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Visvader_0131_TNBC_Navin_Visvader_NORM_panCK_Reference.rds")

########## Prepare Cell Annotations
# Extract & Save Meta Data
Visvader_0131_TNBC_Navin_Visvader_NORM_panCK_Reference.meta.data <- as.data.frame(as.matrix(Visvader_0131_TNBC_Navin_Visvader_NORM_panCK_Reference@meta.data))
# Extract pertinent columns
# Make row names first column
Visvader_0131_TNBC_Navin_Visvader_NORM_panCK_Reference.meta.data <- tibble::rownames_to_column(Visvader_0131_TNBC_Navin_Visvader_NORM_panCK_Reference.meta.data, "Cells")
# In Seurat, "X" becomes the cell names while "cnv.ident" is specified above
Visvader_0131_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations <- Visvader_0131_TNBC_Navin_Visvader_NORM_panCK_Reference.meta.data[, c("Cells", "cnv.ident")]
# Delete header
names(Visvader_0131_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations) <- NULL
# Save as table
write.table(Visvader_0131_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations,"/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_0131_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations.txt",sep="\t",row.names=FALSE)
# Remove metadata and annotations
rm(Visvader_0131_TNBC_Navin_Visvader_NORM_panCK_Reference.meta.data)
rm(Visvader_0131_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations)

########## Generate Count Matrix
# Extract counts from the Seurat object
Visvader_0131_TNBC_Navin_Visvader_NORM_panCK_Reference.count <- Visvader_0131_TNBC_Navin_Visvader_NORM_panCK_Reference[["RNA"]]$counts
# Convert counts to a sparse matrix
Visvader_0131_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as(Visvader_0131_TNBC_Navin_Visvader_NORM_panCK_Reference.count, 'sparseMatrix')
rm(Visvader_0131_TNBC_Navin_Visvader_NORM_panCK_Reference.count)

# Convert to requisite data format
Visvader_0131_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as.data.frame(Visvader_0131_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix)
# Save table
write.csv(Visvader_0131_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_Count_Matrix/Visvader_0131_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv")
# Remove Objects
rm(Visvader_0131_TNBC_Navin_Visvader_NORM_panCK_Reference)
rm(Visvader_0131_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix)

########### Run inferCNV
# Note: check.names = FALSE prevents cell names from having dashes replaced with dots; essential to match annotations file
Visvader_0131_TNBC_inferCNV <- CreateInfercnvObject(raw_counts_matrix=read.csv(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_Count_Matrix/Visvader_0131_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv", header = TRUE, row.names = 1,check.names = FALSE),
                                                    annotations_file=read.table(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_0131_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations.txt", header = FALSE, row.names = 1),
                                                    delim="\t",
                                                    gene_order_file="/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/hg38_gencode_v27.txt",
                                                    ref_group_names=c("Normal"))


# perform infercnv operations to reveal cnv signal
Visvader_0131_TNBC_inferCNV = infercnv::run(Visvader_0131_TNBC_inferCNV,
                                            cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                            out_dir="/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0131_TNBC_panCK/",  # dir is auto-created for storing outputs
                                            cluster_by_groups=T,   # cluster
                                            denoise=T,
                                            HMM=T,
                                            num_ref_groups = 8, analysis_mode = "subclusters")

###########  Integrate data with Seurat object
# Read RDS and add inferCNV results to metadata
Visvader_0131_TNBC_Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Visvader_0131_TNBC_Navin_Visvader_NORM_panCK_Reference.rds") 
Visvader_0131_TNBC_panCK_cnv <- infercnv::add_to_seurat(infercnv_output_path="/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0131_TNBC_panCK/",
                                                        seurat_obj=Visvader_0131_TNBC_Navin_Visvader_NORM_panCK_Reference, # optional
                                                        top_n=10)

# Collect Tumor cells from Seurat Object
Visvader_0131_TNBC_panCK_cnv <- subset(x = Visvader_0131_TNBC_panCK_cnv, subset = cnv.ident == "Tumor")
# Quantify chromosomes with alterations per cell
Visvader_0131_TNBC_panCK_cnv_Pos_Chromosomes <- as.data.frame(rowSums(Visvader_0131_TNBC_panCK_cnv@meta.data[,startsWith(names(Visvader_0131_TNBC_panCK_cnv@meta.data),"has_cnv_")] > 0 ))
Visvader_0131_TNBC_panCK_cnv_Neg_Chromosomes <- as.data.frame(rowSums(Visvader_0131_TNBC_panCK_cnv@meta.data[,startsWith(names(Visvader_0131_TNBC_panCK_cnv@meta.data),"has_cnv_")] == 0))
Visvader_0131_TNBC_panCK_cnv_Chromosome_Quant <- as.data.frame(c(Visvader_0131_TNBC_panCK_cnv_Pos_Chromosomes, Visvader_0131_TNBC_panCK_cnv_Neg_Chromosomes))
colnames(Visvader_0131_TNBC_panCK_cnv_Chromosome_Quant) <- c("Chromosomes CNV Pos", "Chromosomes CNV Neg")
write.csv(Visvader_0131_TNBC_panCK_cnv_Chromosome_Quant, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0131_TNBC_panCK/Visvader_0131_TNBC_panCK_cnv_Chromosome_Quant.csv")

# Save CNV Metadata
Visvader_0131_TNBC_panCK_cnv.meta.data <- as.data.frame(as.matrix(Visvader_0131_TNBC_panCK_cnv@meta.data))
write.csv(Visvader_0131_TNBC_panCK_cnv.meta.data, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0131_TNBC_panCK/Visvader_0131_TNBC_panCK_cnv.meta.data.csv")

# Subset Cells with inferred CNVs -------------------------------------------------------------------------------------------------------------
Visvader_0131_TNBC_panCK_cnv <- subset(Visvader_0131_TNBC_panCK_cnv, cells=rownames(Visvader_0131_TNBC_panCK_cnv@meta.data)[rowSums(Visvader_0131_TNBC_panCK_cnv@meta.data[,startsWith(names(Visvader_0131_TNBC_panCK_cnv@meta.data),"has_cnv_")]) > 0 ])

# Save RDS with inferredCNVs
saveRDS(Visvader_0131_TNBC_panCK_cnv, file = "/R/R_Visvader/Visvader_RDS/RDS_inferCNV/Visvader_0131_TNBC_panCK_cnv.rds")


# Free memory
rm(Visvader_0131_TNBC_inferCNV)
rm(Visvader_0131_TNBC_Navin_Visvader_NORM_panCK_Reference)
rm(Visvader_0131_TNBC_panCK_cnv)
rm(Visvader_0131_TNBC_panCK_cnv_Pos_Chromosomes)
rm(Visvader_0131_TNBC_panCK_cnv_Neg_Chromosomes)
rm(Visvader_0131_TNBC_panCK_cnv_Chromosome_Quant)
rm(Visvader_0131_TNBC_panCK_cnv.meta.data)
gc()



###########################################
#### InferCNV: Visvader_135_TNBC_panCK ####
###########################################

########## Merge with normal mammary glad, extract epithelium, and process via Seurat
# Load RDS
Visvader_135_TNBC_panCK <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_135_TNBC_panCK.rds")
Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Navin_Visvader_NORM_panCK_Reference.rds")

# Set identities 
colnames(x = Visvader_135_TNBC_panCK[[]])
Idents(object = Visvader_135_TNBC_panCK) <- 'Tumor'
Visvader_135_TNBC_panCK[["cnv.ident"]] <- Idents(object = Visvader_135_TNBC_panCK)
levels(x = Visvader_135_TNBC_panCK)

colnames(x = Navin_Visvader_NORM_panCK_Reference[[]])
Idents(object = Navin_Visvader_NORM_panCK_Reference) <- 'Normal'
Navin_Visvader_NORM_panCK_Reference[["cnv.ident"]] <- Idents(object = Navin_Visvader_NORM_panCK_Reference)
levels(x = Navin_Visvader_NORM_panCK_Reference)

# Merge subsets into recombined RDS
Visvader_135_TNBC_Navin_Visvader_NORM_panCK_Reference <- merge(x = Visvader_135_TNBC_panCK, y = Navin_Visvader_NORM_panCK_Reference)
# Join Layers
Visvader_135_TNBC_Navin_Visvader_NORM_panCK_Reference <- JoinLayers(Visvader_135_TNBC_Navin_Visvader_NORM_panCK_Reference)

# Remove subsetted objects
rm(Visvader_135_TNBC_panCK)
rm(Navin_Visvader_NORM_panCK_Reference)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Visvader_135_TNBC_Navin_Visvader_NORM_panCK_Reference <- FindVariableFeatures(Visvader_135_TNBC_Navin_Visvader_NORM_panCK_Reference, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Visvader_135_TNBC_Navin_Visvader_NORM_panCK_Reference)
Visvader_135_TNBC_Navin_Visvader_NORM_panCK_Reference <- ScaleData(Visvader_135_TNBC_Navin_Visvader_NORM_panCK_Reference, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Visvader_135_TNBC_Navin_Visvader_NORM_panCK_Reference <- RunPCA(Visvader_135_TNBC_Navin_Visvader_NORM_panCK_Reference, features = VariableFeatures(object = Visvader_135_TNBC_Navin_Visvader_NORM_panCK_Reference))

# Examine and visualize PCA results a few different ways
print(Visvader_135_TNBC_Navin_Visvader_NORM_panCK_Reference[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Visvader_135_TNBC_Navin_Visvader_NORM_panCK_Reference, dims = 1:2, reduction = "pca")
DimPlot(Visvader_135_TNBC_Navin_Visvader_NORM_panCK_Reference, reduction = "pca")

Visvader_135_TNBC_Navin_Visvader_NORM_panCK_Reference <- FindNeighbors(Visvader_135_TNBC_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
Visvader_135_TNBC_Navin_Visvader_NORM_panCK_Reference <- FindClusters(Visvader_135_TNBC_Navin_Visvader_NORM_panCK_Reference, resolution = 0.1)

# Umap clustering
Visvader_135_TNBC_Navin_Visvader_NORM_panCK_Reference <- RunUMAP(Visvader_135_TNBC_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
DimPlot(Visvader_135_TNBC_Navin_Visvader_NORM_panCK_Reference, reduction = "umap")
DimPlot(Visvader_135_TNBC_Navin_Visvader_NORM_panCK_Reference, reduction = "umap", group.by = "orig.ident")

# Save RDS
saveRDS(Visvader_135_TNBC_Navin_Visvader_NORM_panCK_Reference, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Visvader_135_TNBC_Navin_Visvader_NORM_panCK_Reference.rds")

########## Prepare Cell Annotations
# Extract & Save Meta Data
Visvader_135_TNBC_Navin_Visvader_NORM_panCK_Reference.meta.data <- as.data.frame(as.matrix(Visvader_135_TNBC_Navin_Visvader_NORM_panCK_Reference@meta.data))
# Extract pertinent columns
# Make row names first column
Visvader_135_TNBC_Navin_Visvader_NORM_panCK_Reference.meta.data <- tibble::rownames_to_column(Visvader_135_TNBC_Navin_Visvader_NORM_panCK_Reference.meta.data, "Cells")
# In Seurat, "X" becomes the cell names while "cnv.ident" is specified above
Visvader_135_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations <- Visvader_135_TNBC_Navin_Visvader_NORM_panCK_Reference.meta.data[, c("Cells", "cnv.ident")]
# Delete header
names(Visvader_135_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations) <- NULL
# Save as table
write.table(Visvader_135_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations,"/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_135_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations.txt",sep="\t",row.names=FALSE)
# Remove metadata and annotations
rm(Visvader_135_TNBC_Navin_Visvader_NORM_panCK_Reference.meta.data)
rm(Visvader_135_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations)

########## Generate Count Matrix
# Extract counts from the Seurat object
Visvader_135_TNBC_Navin_Visvader_NORM_panCK_Reference.count <- Visvader_135_TNBC_Navin_Visvader_NORM_panCK_Reference[["RNA"]]$counts
# Convert counts to a sparse matrix
Visvader_135_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as(Visvader_135_TNBC_Navin_Visvader_NORM_panCK_Reference.count, 'sparseMatrix')
rm(Visvader_135_TNBC_Navin_Visvader_NORM_panCK_Reference.count)

# Convert to requisite data format
Visvader_135_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as.data.frame(Visvader_135_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix)
# Save table
write.csv(Visvader_135_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_Count_Matrix/Visvader_135_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv")
# Remove Objects
rm(Visvader_135_TNBC_Navin_Visvader_NORM_panCK_Reference)
rm(Visvader_135_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix)

########### Run inferCNV
# Note: check.names = FALSE prevents cell names from having dashes replaced with dots; essential to match annotations file
Visvader_135_TNBC_inferCNV <- CreateInfercnvObject(raw_counts_matrix=read.csv(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_Count_Matrix/Visvader_135_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv", header = TRUE, row.names = 1,check.names = FALSE),
                                                   annotations_file=read.table(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_135_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations.txt", header = FALSE, row.names = 1),
                                                   delim="\t",
                                                   gene_order_file="/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/hg38_gencode_v27.txt",
                                                   ref_group_names=c("Normal"))


# perform infercnv operations to reveal cnv signal
Visvader_135_TNBC_inferCNV = infercnv::run(Visvader_135_TNBC_inferCNV,
                                           cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                           out_dir="/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_135_TNBC_panCK/",  # dir is auto-created for storing outputs
                                           cluster_by_groups=T,   # cluster
                                           denoise=T,
                                           HMM=T,
                                           num_ref_groups = 8, analysis_mode = "subclusters")

###########  Integrate data with Seurat object
# Read RDS and add inferCNV results to metadata
Visvader_135_TNBC_Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Visvader_135_TNBC_Navin_Visvader_NORM_panCK_Reference.rds") 
Visvader_135_TNBC_panCK_cnv <- infercnv::add_to_seurat(infercnv_output_path="/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_135_TNBC_panCK/",
                                                       seurat_obj=Visvader_135_TNBC_Navin_Visvader_NORM_panCK_Reference, # optional
                                                       top_n=10)

# Collect Tumor cells from Seurat Object
Visvader_135_TNBC_panCK_cnv <- subset(x = Visvader_135_TNBC_panCK_cnv, subset = cnv.ident == "Tumor")
# Quantify chromosomes with alterations per cell
Visvader_135_TNBC_panCK_cnv_Pos_Chromosomes <- as.data.frame(rowSums(Visvader_135_TNBC_panCK_cnv@meta.data[,startsWith(names(Visvader_135_TNBC_panCK_cnv@meta.data),"has_cnv_")] > 0 ))
Visvader_135_TNBC_panCK_cnv_Neg_Chromosomes <- as.data.frame(rowSums(Visvader_135_TNBC_panCK_cnv@meta.data[,startsWith(names(Visvader_135_TNBC_panCK_cnv@meta.data),"has_cnv_")] == 0))
Visvader_135_TNBC_panCK_cnv_Chromosome_Quant <- as.data.frame(c(Visvader_135_TNBC_panCK_cnv_Pos_Chromosomes, Visvader_135_TNBC_panCK_cnv_Neg_Chromosomes))
colnames(Visvader_135_TNBC_panCK_cnv_Chromosome_Quant) <- c("Chromosomes CNV Pos", "Chromosomes CNV Neg")
write.csv(Visvader_135_TNBC_panCK_cnv_Chromosome_Quant, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_135_TNBC_panCK/Visvader_135_TNBC_panCK_cnv_Chromosome_Quant.csv")

# Save CNV Metadata
Visvader_135_TNBC_panCK_cnv.meta.data <- as.data.frame(as.matrix(Visvader_135_TNBC_panCK_cnv@meta.data))
write.csv(Visvader_135_TNBC_panCK_cnv.meta.data, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_135_TNBC_panCK/Visvader_135_TNBC_panCK_cnv.meta.data.csv")

# Subset Cells with inferred CNVs -------------------------------------------------------------------------------------------------------------
Visvader_135_TNBC_panCK_cnv <- subset(Visvader_135_TNBC_panCK_cnv, cells=rownames(Visvader_135_TNBC_panCK_cnv@meta.data)[rowSums(Visvader_135_TNBC_panCK_cnv@meta.data[,startsWith(names(Visvader_135_TNBC_panCK_cnv@meta.data),"has_cnv_")]) > 0 ])

# Save RDS with inferredCNVs
saveRDS(Visvader_135_TNBC_panCK_cnv, file = "/R/R_Visvader/Visvader_RDS/RDS_inferCNV/Visvader_135_TNBC_panCK_cnv.rds")


# Free memory
rm(Visvader_135_TNBC_inferCNV)
rm(Visvader_135_TNBC_Navin_Visvader_NORM_panCK_Reference)
rm(Visvader_135_TNBC_panCK_cnv)
rm(Visvader_135_TNBC_panCK_cnv_Pos_Chromosomes)
rm(Visvader_135_TNBC_panCK_cnv_Neg_Chromosomes)
rm(Visvader_135_TNBC_panCK_cnv_Chromosome_Quant)
rm(Visvader_135_TNBC_panCK_cnv.meta.data)
gc()



############################################
#### InferCNV: Visvader_0177_TNBC_panCK ####
############################################

########## Merge with normal mammary glad, extract epithelium, and process via Seurat
# Load RDS
Visvader_0177_TNBC_panCK <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0177_TNBC_panCK.rds")
Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Navin_Visvader_NORM_panCK_Reference.rds")

# Set identities 
colnames(x = Visvader_0177_TNBC_panCK[[]])
Idents(object = Visvader_0177_TNBC_panCK) <- 'Tumor'
Visvader_0177_TNBC_panCK[["cnv.ident"]] <- Idents(object = Visvader_0177_TNBC_panCK)
levels(x = Visvader_0177_TNBC_panCK)

colnames(x = Navin_Visvader_NORM_panCK_Reference[[]])
Idents(object = Navin_Visvader_NORM_panCK_Reference) <- 'Normal'
Navin_Visvader_NORM_panCK_Reference[["cnv.ident"]] <- Idents(object = Navin_Visvader_NORM_panCK_Reference)
levels(x = Navin_Visvader_NORM_panCK_Reference)

# Merge subsets into recombined RDS
Visvader_0177_TNBC_Navin_Visvader_NORM_panCK_Reference <- merge(x = Visvader_0177_TNBC_panCK, y = Navin_Visvader_NORM_panCK_Reference)
# Join Layers
Visvader_0177_TNBC_Navin_Visvader_NORM_panCK_Reference <- JoinLayers(Visvader_0177_TNBC_Navin_Visvader_NORM_panCK_Reference)

# Remove subsetted objects
rm(Visvader_0177_TNBC_panCK)
rm(Navin_Visvader_NORM_panCK_Reference)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Visvader_0177_TNBC_Navin_Visvader_NORM_panCK_Reference <- FindVariableFeatures(Visvader_0177_TNBC_Navin_Visvader_NORM_panCK_Reference, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Visvader_0177_TNBC_Navin_Visvader_NORM_panCK_Reference)
Visvader_0177_TNBC_Navin_Visvader_NORM_panCK_Reference <- ScaleData(Visvader_0177_TNBC_Navin_Visvader_NORM_panCK_Reference, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Visvader_0177_TNBC_Navin_Visvader_NORM_panCK_Reference <- RunPCA(Visvader_0177_TNBC_Navin_Visvader_NORM_panCK_Reference, features = VariableFeatures(object = Visvader_0177_TNBC_Navin_Visvader_NORM_panCK_Reference))

# Examine and visualize PCA results a few different ways
print(Visvader_0177_TNBC_Navin_Visvader_NORM_panCK_Reference[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Visvader_0177_TNBC_Navin_Visvader_NORM_panCK_Reference, dims = 1:2, reduction = "pca")
DimPlot(Visvader_0177_TNBC_Navin_Visvader_NORM_panCK_Reference, reduction = "pca")

Visvader_0177_TNBC_Navin_Visvader_NORM_panCK_Reference <- FindNeighbors(Visvader_0177_TNBC_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
Visvader_0177_TNBC_Navin_Visvader_NORM_panCK_Reference <- FindClusters(Visvader_0177_TNBC_Navin_Visvader_NORM_panCK_Reference, resolution = 0.1)

# Umap clustering
Visvader_0177_TNBC_Navin_Visvader_NORM_panCK_Reference <- RunUMAP(Visvader_0177_TNBC_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
DimPlot(Visvader_0177_TNBC_Navin_Visvader_NORM_panCK_Reference, reduction = "umap")
DimPlot(Visvader_0177_TNBC_Navin_Visvader_NORM_panCK_Reference, reduction = "umap", group.by = "orig.ident")

# Save RDS
saveRDS(Visvader_0177_TNBC_Navin_Visvader_NORM_panCK_Reference, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Visvader_0177_TNBC_Navin_Visvader_NORM_panCK_Reference.rds")

########## Prepare Cell Annotations
# Extract & Save Meta Data
Visvader_0177_TNBC_Navin_Visvader_NORM_panCK_Reference.meta.data <- as.data.frame(as.matrix(Visvader_0177_TNBC_Navin_Visvader_NORM_panCK_Reference@meta.data))
# Extract pertinent columns
# Make row names first column
Visvader_0177_TNBC_Navin_Visvader_NORM_panCK_Reference.meta.data <- tibble::rownames_to_column(Visvader_0177_TNBC_Navin_Visvader_NORM_panCK_Reference.meta.data, "Cells")
# In Seurat, "X" becomes the cell names while "cnv.ident" is specified above
Visvader_0177_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations <- Visvader_0177_TNBC_Navin_Visvader_NORM_panCK_Reference.meta.data[, c("Cells", "cnv.ident")]
# Delete header
names(Visvader_0177_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations) <- NULL
# Save as table
write.table(Visvader_0177_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations,"/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_0177_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations.txt",sep="\t",row.names=FALSE)
# Remove metadata and annotations
rm(Visvader_0177_TNBC_Navin_Visvader_NORM_panCK_Reference.meta.data)
rm(Visvader_0177_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations)

########## Generate Count Matrix
# Extract counts from the Seurat object
Visvader_0177_TNBC_Navin_Visvader_NORM_panCK_Reference.count <- Visvader_0177_TNBC_Navin_Visvader_NORM_panCK_Reference[["RNA"]]$counts
# Convert counts to a sparse matrix
Visvader_0177_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as(Visvader_0177_TNBC_Navin_Visvader_NORM_panCK_Reference.count, 'sparseMatrix')
rm(Visvader_0177_TNBC_Navin_Visvader_NORM_panCK_Reference.count)

# Convert to requisite data format
Visvader_0177_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as.data.frame(Visvader_0177_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix)
# Save table
write.csv(Visvader_0177_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_Count_Matrix/Visvader_0177_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv")
# Remove Objects
rm(Visvader_0177_TNBC_Navin_Visvader_NORM_panCK_Reference)
rm(Visvader_0177_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix)

########### Run inferCNV
# Note: check.names = FALSE prevents cell names from having dashes replaced with dots; essential to match annotations file
Visvader_0177_TNBC_inferCNV <- CreateInfercnvObject(raw_counts_matrix=read.csv(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_Count_Matrix/Visvader_0177_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv", header = TRUE, row.names = 1,check.names = FALSE),
                                                    annotations_file=read.table(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_0177_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations.txt", header = FALSE, row.names = 1),
                                                    delim="\t",
                                                    gene_order_file="/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/hg38_gencode_v27.txt",
                                                    ref_group_names=c("Normal"))


# perform infercnv operations to reveal cnv signal
Visvader_0177_TNBC_inferCNV = infercnv::run(Visvader_0177_TNBC_inferCNV,
                                            cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                            out_dir="/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0177_TNBC_panCK/",  # dir is auto-created for storing outputs
                                            cluster_by_groups=T,   # cluster
                                            denoise=T,
                                            HMM=T,
                                            num_ref_groups = 8, analysis_mode = "subclusters")

###########  Integrate data with Seurat object
# Read RDS and add inferCNV results to metadata
Visvader_0177_TNBC_Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Visvader_0177_TNBC_Navin_Visvader_NORM_panCK_Reference.rds") 
Visvader_0177_TNBC_panCK_cnv <- infercnv::add_to_seurat(infercnv_output_path="/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0177_TNBC_panCK/",
                                                        seurat_obj=Visvader_0177_TNBC_Navin_Visvader_NORM_panCK_Reference, # optional
                                                        top_n=10)

# Collect Tumor cells from Seurat Object
Visvader_0177_TNBC_panCK_cnv <- subset(x = Visvader_0177_TNBC_panCK_cnv, subset = cnv.ident == "Tumor")
# Quantify chromosomes with alterations per cell
Visvader_0177_TNBC_panCK_cnv_Pos_Chromosomes <- as.data.frame(rowSums(Visvader_0177_TNBC_panCK_cnv@meta.data[,startsWith(names(Visvader_0177_TNBC_panCK_cnv@meta.data),"has_cnv_")] > 0 ))
Visvader_0177_TNBC_panCK_cnv_Neg_Chromosomes <- as.data.frame(rowSums(Visvader_0177_TNBC_panCK_cnv@meta.data[,startsWith(names(Visvader_0177_TNBC_panCK_cnv@meta.data),"has_cnv_")] == 0))
Visvader_0177_TNBC_panCK_cnv_Chromosome_Quant <- as.data.frame(c(Visvader_0177_TNBC_panCK_cnv_Pos_Chromosomes, Visvader_0177_TNBC_panCK_cnv_Neg_Chromosomes))
colnames(Visvader_0177_TNBC_panCK_cnv_Chromosome_Quant) <- c("Chromosomes CNV Pos", "Chromosomes CNV Neg")
write.csv(Visvader_0177_TNBC_panCK_cnv_Chromosome_Quant, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0177_TNBC_panCK/Visvader_0177_TNBC_panCK_cnv_Chromosome_Quant.csv")

# Save CNV Metadata
Visvader_0177_TNBC_panCK_cnv.meta.data <- as.data.frame(as.matrix(Visvader_0177_TNBC_panCK_cnv@meta.data))
write.csv(Visvader_0177_TNBC_panCK_cnv.meta.data, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0177_TNBC_panCK/Visvader_0177_TNBC_panCK_cnv.meta.data.csv")

# Subset Cells with inferred CNVs -------------------------------------------------------------------------------------------------------------
Visvader_0177_TNBC_panCK_cnv <- subset(Visvader_0177_TNBC_panCK_cnv, cells=rownames(Visvader_0177_TNBC_panCK_cnv@meta.data)[rowSums(Visvader_0177_TNBC_panCK_cnv@meta.data[,startsWith(names(Visvader_0177_TNBC_panCK_cnv@meta.data),"has_cnv_")]) > 0 ])

# Save RDS with inferredCNVs
saveRDS(Visvader_0177_TNBC_panCK_cnv, file = "/R/R_Visvader/Visvader_RDS/RDS_inferCNV/Visvader_0177_TNBC_panCK_cnv.rds")


# Free memory
rm(Visvader_0177_TNBC_inferCNV)
rm(Visvader_0177_TNBC_Navin_Visvader_NORM_panCK_Reference)
rm(Visvader_0177_TNBC_panCK_cnv)
rm(Visvader_0177_TNBC_panCK_cnv_Pos_Chromosomes)
rm(Visvader_0177_TNBC_panCK_cnv_Neg_Chromosomes)
rm(Visvader_0177_TNBC_panCK_cnv_Chromosome_Quant)
rm(Visvader_0177_TNBC_panCK_cnv.meta.data)
gc()



############################################
#### InferCNV: Visvader_0554_TNBC_panCK ####
############################################

########## Merge with normal mammary glad, extract epithelium, and process via Seurat
# Load RDS
Visvader_0554_TNBC_panCK <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0554_TNBC_panCK.rds")
Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Navin_Visvader_NORM_panCK_Reference.rds")

# Set identities 
colnames(x = Visvader_0554_TNBC_panCK[[]])
Idents(object = Visvader_0554_TNBC_panCK) <- 'Tumor'
Visvader_0554_TNBC_panCK[["cnv.ident"]] <- Idents(object = Visvader_0554_TNBC_panCK)
levels(x = Visvader_0554_TNBC_panCK)

colnames(x = Navin_Visvader_NORM_panCK_Reference[[]])
Idents(object = Navin_Visvader_NORM_panCK_Reference) <- 'Normal'
Navin_Visvader_NORM_panCK_Reference[["cnv.ident"]] <- Idents(object = Navin_Visvader_NORM_panCK_Reference)
levels(x = Navin_Visvader_NORM_panCK_Reference)

# Merge subsets into recombined RDS
Visvader_0554_TNBC_Navin_Visvader_NORM_panCK_Reference <- merge(x = Visvader_0554_TNBC_panCK, y = Navin_Visvader_NORM_panCK_Reference)
# Join Layers
Visvader_0554_TNBC_Navin_Visvader_NORM_panCK_Reference <- JoinLayers(Visvader_0554_TNBC_Navin_Visvader_NORM_panCK_Reference)

# Remove subsetted objects
rm(Visvader_0554_TNBC_panCK)
rm(Navin_Visvader_NORM_panCK_Reference)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Visvader_0554_TNBC_Navin_Visvader_NORM_panCK_Reference <- FindVariableFeatures(Visvader_0554_TNBC_Navin_Visvader_NORM_panCK_Reference, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Visvader_0554_TNBC_Navin_Visvader_NORM_panCK_Reference)
Visvader_0554_TNBC_Navin_Visvader_NORM_panCK_Reference <- ScaleData(Visvader_0554_TNBC_Navin_Visvader_NORM_panCK_Reference, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Visvader_0554_TNBC_Navin_Visvader_NORM_panCK_Reference <- RunPCA(Visvader_0554_TNBC_Navin_Visvader_NORM_panCK_Reference, features = VariableFeatures(object = Visvader_0554_TNBC_Navin_Visvader_NORM_panCK_Reference))

# Examine and visualize PCA results a few different ways
print(Visvader_0554_TNBC_Navin_Visvader_NORM_panCK_Reference[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Visvader_0554_TNBC_Navin_Visvader_NORM_panCK_Reference, dims = 1:2, reduction = "pca")
DimPlot(Visvader_0554_TNBC_Navin_Visvader_NORM_panCK_Reference, reduction = "pca")

Visvader_0554_TNBC_Navin_Visvader_NORM_panCK_Reference <- FindNeighbors(Visvader_0554_TNBC_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
Visvader_0554_TNBC_Navin_Visvader_NORM_panCK_Reference <- FindClusters(Visvader_0554_TNBC_Navin_Visvader_NORM_panCK_Reference, resolution = 0.1)

# Umap clustering
Visvader_0554_TNBC_Navin_Visvader_NORM_panCK_Reference <- RunUMAP(Visvader_0554_TNBC_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
DimPlot(Visvader_0554_TNBC_Navin_Visvader_NORM_panCK_Reference, reduction = "umap")
DimPlot(Visvader_0554_TNBC_Navin_Visvader_NORM_panCK_Reference, reduction = "umap", group.by = "orig.ident")

# Save RDS
saveRDS(Visvader_0554_TNBC_Navin_Visvader_NORM_panCK_Reference, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Visvader_0554_TNBC_Navin_Visvader_NORM_panCK_Reference.rds")

########## Prepare Cell Annotations
# Extract & Save Meta Data
Visvader_0554_TNBC_Navin_Visvader_NORM_panCK_Reference.meta.data <- as.data.frame(as.matrix(Visvader_0554_TNBC_Navin_Visvader_NORM_panCK_Reference@meta.data))
# Extract pertinent columns
# Make row names first column
Visvader_0554_TNBC_Navin_Visvader_NORM_panCK_Reference.meta.data <- tibble::rownames_to_column(Visvader_0554_TNBC_Navin_Visvader_NORM_panCK_Reference.meta.data, "Cells")
# In Seurat, "X" becomes the cell names while "cnv.ident" is specified above
Visvader_0554_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations <- Visvader_0554_TNBC_Navin_Visvader_NORM_panCK_Reference.meta.data[, c("Cells", "cnv.ident")]
# Delete header
names(Visvader_0554_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations) <- NULL
# Save as table
write.table(Visvader_0554_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations,"/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_0554_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations.txt",sep="\t",row.names=FALSE)
# Remove metadata and annotations
rm(Visvader_0554_TNBC_Navin_Visvader_NORM_panCK_Reference.meta.data)
rm(Visvader_0554_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations)

########## Generate Count Matrix
# Extract counts from the Seurat object
Visvader_0554_TNBC_Navin_Visvader_NORM_panCK_Reference.count <- Visvader_0554_TNBC_Navin_Visvader_NORM_panCK_Reference[["RNA"]]$counts
# Convert counts to a sparse matrix
Visvader_0554_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as(Visvader_0554_TNBC_Navin_Visvader_NORM_panCK_Reference.count, 'sparseMatrix')
rm(Visvader_0554_TNBC_Navin_Visvader_NORM_panCK_Reference.count)

# Convert to requisite data format
Visvader_0554_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as.data.frame(Visvader_0554_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix)
# Save table
write.csv(Visvader_0554_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_Count_Matrix/Visvader_0554_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv")
# Remove Objects
rm(Visvader_0554_TNBC_Navin_Visvader_NORM_panCK_Reference)
rm(Visvader_0554_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix)

########### Run inferCNV
# Note: check.names = FALSE prevents cell names from having dashes replaced with dots; essential to match annotations file
Visvader_0554_TNBC_inferCNV <- CreateInfercnvObject(raw_counts_matrix=read.csv(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_Count_Matrix/Visvader_0554_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv", header = TRUE, row.names = 1,check.names = FALSE),
                                                    annotations_file=read.table(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_0554_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations.txt", header = FALSE, row.names = 1),
                                                    delim="\t",
                                                    gene_order_file="/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/hg38_gencode_v27.txt",
                                                    ref_group_names=c("Normal"))


# perform infercnv operations to reveal cnv signal
Visvader_0554_TNBC_inferCNV = infercnv::run(Visvader_0554_TNBC_inferCNV,
                                            cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                            out_dir="/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0554_TNBC_panCK/",  # dir is auto-created for storing outputs
                                            cluster_by_groups=T,   # cluster
                                            denoise=T,
                                            HMM=T,
                                            num_ref_groups = 8, analysis_mode = "subclusters")

###########  Integrate data with Seurat object
# Read RDS and add inferCNV results to metadata
Visvader_0554_TNBC_Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Visvader_0554_TNBC_Navin_Visvader_NORM_panCK_Reference.rds") 
Visvader_0554_TNBC_panCK_cnv <- infercnv::add_to_seurat(infercnv_output_path="/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0554_TNBC_panCK/",
                                                        seurat_obj=Visvader_0554_TNBC_Navin_Visvader_NORM_panCK_Reference, # optional
                                                        top_n=10)

# Collect Tumor cells from Seurat Object
Visvader_0554_TNBC_panCK_cnv <- subset(x = Visvader_0554_TNBC_panCK_cnv, subset = cnv.ident == "Tumor")
# Quantify chromosomes with alterations per cell
Visvader_0554_TNBC_panCK_cnv_Pos_Chromosomes <- as.data.frame(rowSums(Visvader_0554_TNBC_panCK_cnv@meta.data[,startsWith(names(Visvader_0554_TNBC_panCK_cnv@meta.data),"has_cnv_")] > 0 ))
Visvader_0554_TNBC_panCK_cnv_Neg_Chromosomes <- as.data.frame(rowSums(Visvader_0554_TNBC_panCK_cnv@meta.data[,startsWith(names(Visvader_0554_TNBC_panCK_cnv@meta.data),"has_cnv_")] == 0))
Visvader_0554_TNBC_panCK_cnv_Chromosome_Quant <- as.data.frame(c(Visvader_0554_TNBC_panCK_cnv_Pos_Chromosomes, Visvader_0554_TNBC_panCK_cnv_Neg_Chromosomes))
colnames(Visvader_0554_TNBC_panCK_cnv_Chromosome_Quant) <- c("Chromosomes CNV Pos", "Chromosomes CNV Neg")
write.csv(Visvader_0554_TNBC_panCK_cnv_Chromosome_Quant, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0554_TNBC_panCK/Visvader_0554_TNBC_panCK_cnv_Chromosome_Quant.csv")

# Save CNV Metadata
Visvader_0554_TNBC_panCK_cnv.meta.data <- as.data.frame(as.matrix(Visvader_0554_TNBC_panCK_cnv@meta.data))
write.csv(Visvader_0554_TNBC_panCK_cnv.meta.data, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_0554_TNBC_panCK/Visvader_0554_TNBC_panCK_cnv.meta.data.csv")

# Subset Cells with inferred CNVs -------------------------------------------------------------------------------------------------------------
Visvader_0554_TNBC_panCK_cnv <- subset(Visvader_0554_TNBC_panCK_cnv, cells=rownames(Visvader_0554_TNBC_panCK_cnv@meta.data)[rowSums(Visvader_0554_TNBC_panCK_cnv@meta.data[,startsWith(names(Visvader_0554_TNBC_panCK_cnv@meta.data),"has_cnv_")]) > 0 ])

# Save RDS with inferredCNVs
saveRDS(Visvader_0554_TNBC_panCK_cnv, file = "/R/R_Visvader/Visvader_RDS/RDS_inferCNV/Visvader_0554_TNBC_panCK_cnv.rds")


# Free memory
rm(Visvader_0554_TNBC_inferCNV)
rm(Visvader_0554_TNBC_Navin_Visvader_NORM_panCK_Reference)
rm(Visvader_0554_TNBC_panCK_cnv)
rm(Visvader_0554_TNBC_panCK_cnv_Pos_Chromosomes)
rm(Visvader_0554_TNBC_panCK_cnv_Neg_Chromosomes)
rm(Visvader_0554_TNBC_panCK_cnv_Chromosome_Quant)
rm(Visvader_0554_TNBC_panCK_cnv.meta.data)
gc()



############################################
#### InferCNV: Visvader_4031_TNBC_panCK ####
############################################

########## Merge with normal mammary glad, extract epithelium, and process via Seurat
# Load RDS
Visvader_4031_TNBC_panCK <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_4031_TNBC_panCK.rds")
Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Navin_Visvader_NORM_panCK_Reference.rds")

# Set identities 
colnames(x = Visvader_4031_TNBC_panCK[[]])
Idents(object = Visvader_4031_TNBC_panCK) <- 'Tumor'
Visvader_4031_TNBC_panCK[["cnv.ident"]] <- Idents(object = Visvader_4031_TNBC_panCK)
levels(x = Visvader_4031_TNBC_panCK)

colnames(x = Navin_Visvader_NORM_panCK_Reference[[]])
Idents(object = Navin_Visvader_NORM_panCK_Reference) <- 'Normal'
Navin_Visvader_NORM_panCK_Reference[["cnv.ident"]] <- Idents(object = Navin_Visvader_NORM_panCK_Reference)
levels(x = Navin_Visvader_NORM_panCK_Reference)

# Merge subsets into recombined RDS
Visvader_4031_TNBC_Navin_Visvader_NORM_panCK_Reference <- merge(x = Visvader_4031_TNBC_panCK, y = Navin_Visvader_NORM_panCK_Reference)
# Join Layers
Visvader_4031_TNBC_Navin_Visvader_NORM_panCK_Reference <- JoinLayers(Visvader_4031_TNBC_Navin_Visvader_NORM_panCK_Reference)

# Remove subsetted objects
rm(Visvader_4031_TNBC_panCK)
rm(Navin_Visvader_NORM_panCK_Reference)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Visvader_4031_TNBC_Navin_Visvader_NORM_panCK_Reference <- FindVariableFeatures(Visvader_4031_TNBC_Navin_Visvader_NORM_panCK_Reference, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Visvader_4031_TNBC_Navin_Visvader_NORM_panCK_Reference)
Visvader_4031_TNBC_Navin_Visvader_NORM_panCK_Reference <- ScaleData(Visvader_4031_TNBC_Navin_Visvader_NORM_panCK_Reference, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Visvader_4031_TNBC_Navin_Visvader_NORM_panCK_Reference <- RunPCA(Visvader_4031_TNBC_Navin_Visvader_NORM_panCK_Reference, features = VariableFeatures(object = Visvader_4031_TNBC_Navin_Visvader_NORM_panCK_Reference))

# Examine and visualize PCA results a few different ways
print(Visvader_4031_TNBC_Navin_Visvader_NORM_panCK_Reference[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Visvader_4031_TNBC_Navin_Visvader_NORM_panCK_Reference, dims = 1:2, reduction = "pca")
DimPlot(Visvader_4031_TNBC_Navin_Visvader_NORM_panCK_Reference, reduction = "pca")

Visvader_4031_TNBC_Navin_Visvader_NORM_panCK_Reference <- FindNeighbors(Visvader_4031_TNBC_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
Visvader_4031_TNBC_Navin_Visvader_NORM_panCK_Reference <- FindClusters(Visvader_4031_TNBC_Navin_Visvader_NORM_panCK_Reference, resolution = 0.1)

# Umap clustering
Visvader_4031_TNBC_Navin_Visvader_NORM_panCK_Reference <- RunUMAP(Visvader_4031_TNBC_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
DimPlot(Visvader_4031_TNBC_Navin_Visvader_NORM_panCK_Reference, reduction = "umap")
DimPlot(Visvader_4031_TNBC_Navin_Visvader_NORM_panCK_Reference, reduction = "umap", group.by = "orig.ident")

# Save RDS
saveRDS(Visvader_4031_TNBC_Navin_Visvader_NORM_panCK_Reference, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Visvader_4031_TNBC_Navin_Visvader_NORM_panCK_Reference.rds")

########## Prepare Cell Annotations
# Extract & Save Meta Data
Visvader_4031_TNBC_Navin_Visvader_NORM_panCK_Reference.meta.data <- as.data.frame(as.matrix(Visvader_4031_TNBC_Navin_Visvader_NORM_panCK_Reference@meta.data))
# Extract pertinent columns
# Make row names first column
Visvader_4031_TNBC_Navin_Visvader_NORM_panCK_Reference.meta.data <- tibble::rownames_to_column(Visvader_4031_TNBC_Navin_Visvader_NORM_panCK_Reference.meta.data, "Cells")
# In Seurat, "X" becomes the cell names while "cnv.ident" is specified above
Visvader_4031_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations <- Visvader_4031_TNBC_Navin_Visvader_NORM_panCK_Reference.meta.data[, c("Cells", "cnv.ident")]
# Delete header
names(Visvader_4031_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations) <- NULL
# Save as table
write.table(Visvader_4031_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations,"/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_4031_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations.txt",sep="\t",row.names=FALSE)
# Remove metadata and annotations
rm(Visvader_4031_TNBC_Navin_Visvader_NORM_panCK_Reference.meta.data)
rm(Visvader_4031_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations)

########## Generate Count Matrix
# Extract counts from the Seurat object
Visvader_4031_TNBC_Navin_Visvader_NORM_panCK_Reference.count <- Visvader_4031_TNBC_Navin_Visvader_NORM_panCK_Reference[["RNA"]]$counts
# Convert counts to a sparse matrix
Visvader_4031_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as(Visvader_4031_TNBC_Navin_Visvader_NORM_panCK_Reference.count, 'sparseMatrix')
rm(Visvader_4031_TNBC_Navin_Visvader_NORM_panCK_Reference.count)

# Convert to requisite data format
Visvader_4031_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as.data.frame(Visvader_4031_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix)
# Save table
write.csv(Visvader_4031_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_Count_Matrix/Visvader_4031_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv")
# Remove Objects
rm(Visvader_4031_TNBC_Navin_Visvader_NORM_panCK_Reference)
rm(Visvader_4031_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix)

########### Run inferCNV
# Note: check.names = FALSE prevents cell names from having dashes replaced with dots; essential to match annotations file
Visvader_4031_TNBC_inferCNV <- CreateInfercnvObject(raw_counts_matrix=read.csv(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_Count_Matrix/Visvader_4031_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv", header = TRUE, row.names = 1,check.names = FALSE),
                                                    annotations_file=read.table(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_4031_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations.txt", header = FALSE, row.names = 1),
                                                    delim="\t",
                                                    gene_order_file="/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/hg38_gencode_v27.txt",
                                                    ref_group_names=c("Normal"))


# perform infercnv operations to reveal cnv signal
Visvader_4031_TNBC_inferCNV = infercnv::run(Visvader_4031_TNBC_inferCNV,
                                            cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                            out_dir="/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_4031_TNBC_panCK/",  # dir is auto-created for storing outputs
                                            cluster_by_groups=T,   # cluster
                                            denoise=T,
                                            HMM=T,
                                            num_ref_groups = 8, analysis_mode = "subclusters")

###########  Integrate data with Seurat object
# Read RDS and add inferCNV results to metadata
Visvader_4031_TNBC_Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Visvader_4031_TNBC_Navin_Visvader_NORM_panCK_Reference.rds") 
Visvader_4031_TNBC_panCK_cnv <- infercnv::add_to_seurat(infercnv_output_path="/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_4031_TNBC_panCK/",
                                                        seurat_obj=Visvader_4031_TNBC_Navin_Visvader_NORM_panCK_Reference, # optional
                                                        top_n=10)

# Collect Tumor cells from Seurat Object
Visvader_4031_TNBC_panCK_cnv <- subset(x = Visvader_4031_TNBC_panCK_cnv, subset = cnv.ident == "Tumor")
# Quantify chromosomes with alterations per cell
Visvader_4031_TNBC_panCK_cnv_Pos_Chromosomes <- as.data.frame(rowSums(Visvader_4031_TNBC_panCK_cnv@meta.data[,startsWith(names(Visvader_4031_TNBC_panCK_cnv@meta.data),"has_cnv_")] > 0 ))
Visvader_4031_TNBC_panCK_cnv_Neg_Chromosomes <- as.data.frame(rowSums(Visvader_4031_TNBC_panCK_cnv@meta.data[,startsWith(names(Visvader_4031_TNBC_panCK_cnv@meta.data),"has_cnv_")] == 0))
Visvader_4031_TNBC_panCK_cnv_Chromosome_Quant <- as.data.frame(c(Visvader_4031_TNBC_panCK_cnv_Pos_Chromosomes, Visvader_4031_TNBC_panCK_cnv_Neg_Chromosomes))
colnames(Visvader_4031_TNBC_panCK_cnv_Chromosome_Quant) <- c("Chromosomes CNV Pos", "Chromosomes CNV Neg")
write.csv(Visvader_4031_TNBC_panCK_cnv_Chromosome_Quant, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_4031_TNBC_panCK/Visvader_4031_TNBC_panCK_cnv_Chromosome_Quant.csv")

# Save CNV Metadata
Visvader_4031_TNBC_panCK_cnv.meta.data <- as.data.frame(as.matrix(Visvader_4031_TNBC_panCK_cnv@meta.data))
write.csv(Visvader_4031_TNBC_panCK_cnv.meta.data, file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Output/Visvader_4031_TNBC_panCK/Visvader_4031_TNBC_panCK_cnv.meta.data.csv")

# Subset Cells with inferred CNVs -------------------------------------------------------------------------------------------------------------
Visvader_4031_TNBC_panCK_cnv <- subset(Visvader_4031_TNBC_panCK_cnv, cells=rownames(Visvader_4031_TNBC_panCK_cnv@meta.data)[rowSums(Visvader_4031_TNBC_panCK_cnv@meta.data[,startsWith(names(Visvader_4031_TNBC_panCK_cnv@meta.data),"has_cnv_")]) > 0 ])

# Save RDS with inferredCNVs
saveRDS(Visvader_4031_TNBC_panCK_cnv, file = "/R/R_Visvader/Visvader_RDS/RDS_inferCNV/Visvader_4031_TNBC_panCK_cnv.rds")


# Free memory
rm(Visvader_4031_TNBC_inferCNV)
rm(Visvader_4031_TNBC_Navin_Visvader_NORM_panCK_Reference)
rm(Visvader_4031_TNBC_panCK_cnv)
rm(Visvader_4031_TNBC_panCK_cnv_Pos_Chromosomes)
rm(Visvader_4031_TNBC_panCK_cnv_Neg_Chromosomes)
rm(Visvader_4031_TNBC_panCK_cnv_Chromosome_Quant)
rm(Visvader_4031_TNBC_panCK_cnv.meta.data)
gc()




#############################################################################
# Step 6: Merge Neoplastic Cells from Different Tumors into Combined Object #  
#############################################################################
# Load Libraries
library(Seurat)
library(dplyr)
library(patchwork)

#########################################################
############# This is ER/PR Tumor Batch # ###############
#########################################################


#########################
# Load Individual Files #
#########################
Visvader_0001_ER_panCK_cnv <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_inferCNV/Visvader_0001_ER_panCK_cnv.rds")
Visvader_0025_ER_panCK_cnv <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_inferCNV/Visvader_0025_ER_panCK_cnv.rds")
Visvader_0029_7C_ER_panCK_cnv <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_inferCNV/Visvader_0029_7C_ER_panCK_cnv.rds")
Visvader_0029_9C_ER_panCK_cnv <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_inferCNV/Visvader_0029_9C_ER_panCK_cnv.rds")
Visvader_0032_ER_panCK_cnv <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_inferCNV/Visvader_0032_ER_panCK_cnv.rds")
Visvader_0040_ER_panCK_cnv <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_inferCNV/Visvader_0040_ER_panCK_cnv.rds")
Visvader_0042_ER_panCK_cnv <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_inferCNV/Visvader_0042_ER_panCK_cnv.rds")
Visvader_0043_ER_panCK_cnv <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_inferCNV/Visvader_0043_ER_panCK_cnv.rds")
Visvader_0056_ER_panCK_cnv <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_inferCNV/Visvader_0056_ER_panCK_cnv.rds")
Visvader_0064_ER_panCK_cnv <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_inferCNV/Visvader_0064_ER_panCK_cnv.rds")
Visvader_0068_ER_panCK_cnv <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_inferCNV/Visvader_0068_ER_panCK_cnv.rds")
Visvader_0114_ER_panCK_cnv <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_inferCNV/Visvader_0114_ER_panCK_cnv.rds")
Visvader_0125_ER_panCK_cnv <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_inferCNV/Visvader_0125_ER_panCK_cnv.rds")
Visvader_0151_ER_panCK_cnv <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_inferCNV/Visvader_0151_ER_panCK_cnv.rds")
Visvader_0163_ER_panCK_cnv <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_inferCNV/Visvader_0163_ER_panCK_cnv.rds")
Visvader_0167_ER_panCK_cnv <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_inferCNV/Visvader_0167_ER_panCK_cnv.rds")
Visvader_0173_ER_panCK_cnv <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_inferCNV/Visvader_0173_ER_panCK_cnv.rds")
Visvader_0178_ER_panCK_cnv <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_inferCNV/Visvader_0178_ER_panCK_cnv.rds")
Visvader_0360_ER_panCK_cnv <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_inferCNV/Visvader_0360_ER_panCK_cnv.rds")
Visvader_0319_PR_panCK_cnv <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_inferCNV/Visvader_0319_PR_panCK_cnv.rds")



##############################################
# Merge Individual Files into Single Dataset #
##############################################
# Merge Datasets
Visvader_merged_HR_panCK_CNV <- merge(x = Visvader_0001_ER_panCK_cnv, y = c(Visvader_0025_ER_panCK_cnv, Visvader_0029_7C_ER_panCK_cnv, Visvader_0029_9C_ER_panCK_cnv, 
                                                                            Visvader_0032_ER_panCK_cnv, Visvader_0040_ER_panCK_cnv, Visvader_0042_ER_panCK_cnv,
                                                                            Visvader_0043_ER_panCK_cnv, Visvader_0056_ER_panCK_cnv, Visvader_0064_ER_panCK_cnv,
                                                                            Visvader_0068_ER_panCK_cnv, Visvader_0114_ER_panCK_cnv, Visvader_0125_ER_panCK_cnv,
                                                                            Visvader_0151_ER_panCK_cnv, Visvader_0163_ER_panCK_cnv, Visvader_0167_ER_panCK_cnv,
                                                                            Visvader_0173_ER_panCK_cnv, Visvader_0178_ER_panCK_cnv, Visvader_0360_ER_panCK_cnv,
                                                                            Visvader_0319_PR_panCK_cnv))

# Join Layers
Visvader_merged_HR_panCK_CNV <- JoinLayers(Visvader_merged_HR_panCK_CNV)

# Remove Objects
rm(Visvader_0001_ER_panCK_cnv)
rm(Visvader_0025_ER_panCK_cnv)
rm(Visvader_0029_7C_ER_panCK_cnv)
rm(Visvader_0029_9C_ER_panCK_cnv)
rm(Visvader_0032_ER_panCK_cnv)
rm(Visvader_0040_ER_panCK_cnv)
rm(Visvader_0042_ER_panCK_cnv)
rm(Visvader_0043_ER_panCK_cnv)
rm(Visvader_0056_ER_panCK_cnv)
rm(Visvader_0064_ER_panCK_cnv)
rm(Visvader_0068_ER_panCK_cnv)
rm(Visvader_0114_ER_panCK_cnv)
rm(Visvader_0125_ER_panCK_cnv)
rm(Visvader_0151_ER_panCK_cnv)
rm(Visvader_0163_ER_panCK_cnv)
rm(Visvader_0167_ER_panCK_cnv)
rm(Visvader_0173_ER_panCK_cnv)
rm(Visvader_0178_ER_panCK_cnv)
rm(Visvader_0360_ER_panCK_cnv)
rm(Visvader_0319_PR_panCK_cnv)
gc()


####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Visvader_merged_HR_panCK_CNV <- FindVariableFeatures(Visvader_merged_HR_panCK_CNV, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Visvader_merged_HR_panCK_CNV)
Visvader_merged_HR_panCK_CNV <- ScaleData(Visvader_merged_HR_panCK_CNV, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Visvader_merged_HR_panCK_CNV <- RunPCA(Visvader_merged_HR_panCK_CNV, features = VariableFeatures(object = Visvader_merged_HR_panCK_CNV))

# Examine and visualize PCA results a few different ways
print(Visvader_merged_HR_panCK_CNV[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Visvader_merged_HR_panCK_CNV, dims = 1:2, reduction = "pca")
DimPlot(Visvader_merged_HR_panCK_CNV, reduction = "pca")

# Visualize PCs
ElbowPlot(Visvader_merged_HR_panCK_CNV)

# Choose dimensions
Visvader_merged_HR_panCK_CNV <- FindNeighbors(Visvader_merged_HR_panCK_CNV, dims = 1:40)
Visvader_merged_HR_panCK_CNV <- FindClusters(Visvader_merged_HR_panCK_CNV, resolution = 0.12)

# Umap clustering
Visvader_merged_HR_panCK_CNV <- RunUMAP(Visvader_merged_HR_panCK_CNV, dims = 1:40)
DimPlot(Visvader_merged_HR_panCK_CNV, reduction = "umap")

# Identify Samples
DimPlot(Visvader_merged_HR_panCK_CNV, reduction = "umap", group.by = "orig.ident")

# Save
saveRDS(Visvader_merged_HR_panCK_CNV, file = "/R/R_Visvader/Visvader_RDS/RDS_Merged/Visvader_merged_HR_panCK_CNV.rds")
# Clear global environment
rm(Visvader_merged_HR_panCK_CNV)
gc()


#######################################
# Merge and Subset HER2 Tumor Samples #
#######################################

#########################
# Load Individual Files #
#########################
Visvader_0031_HER2_panCK_cnv <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_inferCNV/Visvader_0031_HER2_panCK_cnv.rds")
Visvader_0069_HER2_panCK_cnv <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_inferCNV/Visvader_0069_HER2_panCK_cnv.rds")
Visvader_0161_HER2_panCK_cnv <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_inferCNV/Visvader_0161_HER2_panCK_cnv.rds")
Visvader_0176_HER2_panCK_cnv <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_inferCNV/Visvader_0176_HER2_panCK_cnv.rds")
Visvader_0308_HER2_panCK_cnv <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_inferCNV/Visvader_0308_HER2_panCK_cnv.rds")
Visvader_0337_HER2_panCK_cnv <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_inferCNV/Visvader_0337_HER2_panCK_cnv.rds")


##############################################
# Merge Individual Files into Single Dataset #
##############################################
# Merge Datasets
Visvader_merged_HER2_panCK_CNV <- merge(x = Visvader_0031_HER2_panCK_cnv, y = c(Visvader_0069_HER2_panCK_cnv, Visvader_0161_HER2_panCK_cnv, Visvader_0176_HER2_panCK_cnv,
                                                                                Visvader_0308_HER2_panCK_cnv, Visvader_0337_HER2_panCK_cnv))

# Join Layers
Visvader_merged_HER2_panCK_CNV <- JoinLayers(Visvader_merged_HER2_panCK_CNV)

# Remove Objects
rm(Visvader_0031_HER2_panCK_cnv)
rm(Visvader_0069_HER2_panCK_cnv)
rm(Visvader_0161_HER2_panCK_cnv)
rm(Visvader_0176_HER2_panCK_cnv)
rm(Visvader_0308_HER2_panCK_cnv)
rm(Visvader_0337_HER2_panCK_cnv)
gc()


####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Visvader_merged_HER2_panCK_CNV <- FindVariableFeatures(Visvader_merged_HER2_panCK_CNV, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Visvader_merged_HER2_panCK_CNV)
Visvader_merged_HER2_panCK_CNV <- ScaleData(Visvader_merged_HER2_panCK_CNV, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Visvader_merged_HER2_panCK_CNV <- RunPCA(Visvader_merged_HER2_panCK_CNV, features = VariableFeatures(object = Visvader_merged_HER2_panCK_CNV))

# Examine and visualize PCA results a few different ways
print(Visvader_merged_HER2_panCK_CNV[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Visvader_merged_HER2_panCK_CNV, dims = 1:2, reduction = "pca")
DimPlot(Visvader_merged_HER2_panCK_CNV, reduction = "pca")

# Visualize PCs
ElbowPlot(Visvader_merged_HER2_panCK_CNV)

# Choose dimensions
Visvader_merged_HER2_panCK_CNV <- FindNeighbors(Visvader_merged_HER2_panCK_CNV, dims = 1:40)
Visvader_merged_HER2_panCK_CNV <- FindClusters(Visvader_merged_HER2_panCK_CNV, resolution = 0.12)

# Umap clustering
Visvader_merged_HER2_panCK_CNV <- RunUMAP(Visvader_merged_HER2_panCK_CNV, dims = 1:40)
DimPlot(Visvader_merged_HER2_panCK_CNV, reduction = "umap")

# Identify Samples
DimPlot(Visvader_merged_HER2_panCK_CNV, reduction = "umap", group.by = "orig.ident")

# Save
saveRDS(Visvader_merged_HER2_panCK_CNV, file = "/R/R_Visvader/Visvader_RDS/RDS_Merged/Visvader_merged_HER2_panCK_CNV.rds")

# Clear global environment
rm(Visvader_merged_HER2_panCK_CNV)
gc()




#######################################
# Merge and Subset TNBC Tumor Samples #
#######################################

#########################
# Load Individual Files #
#########################
Visvader_106_TNBC_panCK_cnv <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_inferCNV/Visvader_106_TNBC_panCK_cnv.rds")
Visvader_114_TNBC_panCK_cnv <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_inferCNV/Visvader_114_TNBC_panCK_cnv.rds")
Visvader_126_TNBC_panCK_cnv <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_inferCNV/Visvader_126_TNBC_panCK_cnv.rds")
Visvader_0131_TNBC_panCK_cnv <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_inferCNV/Visvader_0131_TNBC_panCK_cnv.rds")
Visvader_135_TNBC_panCK_cnv <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_inferCNV/Visvader_135_TNBC_panCK_cnv.rds")
Visvader_0177_TNBC_panCK_cnv <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_inferCNV/Visvader_0177_TNBC_panCK_cnv.rds")
Visvader_0554_TNBC_panCK_cnv <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_inferCNV/Visvader_0554_TNBC_panCK_cnv.rds")
Visvader_4031_TNBC_panCK_cnv <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_inferCNV/Visvader_4031_TNBC_panCK_cnv.rds")


##############################################
# Merge Individual Files into Single Dataset #
##############################################
# Merge Datasets
Visvader_merged_TNBC_panCK_cnv <- merge(x = Visvader_106_TNBC_panCK_cnv, y = c(Visvader_114_TNBC_panCK_cnv, Visvader_126_TNBC_panCK_cnv, Visvader_0131_TNBC_panCK_cnv,
                                                                               Visvader_135_TNBC_panCK_cnv, Visvader_0177_TNBC_panCK_cnv, Visvader_0554_TNBC_panCK_cnv,
                                                                               Visvader_4031_TNBC_panCK_cnv))

# Join Layers
Visvader_merged_TNBC_panCK_cnv <- JoinLayers(Visvader_merged_TNBC_panCK_cnv)

# Remove Objects
rm(Visvader_106_TNBC_panCK_cnv)
rm(Visvader_114_TNBC_panCK_cnv)
rm(Visvader_126_TNBC_panCK_cnv)
rm(Visvader_0131_TNBC_panCK_cnv)
rm(Visvader_135_TNBC_panCK_cnv)
rm(Visvader_0177_TNBC_panCK_cnv)
rm(Visvader_0554_TNBC_panCK_cnv)
rm(Visvader_4031_TNBC_panCK_cnv)
gc()


####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
Visvader_merged_TNBC_panCK_cnv <- FindVariableFeatures(Visvader_merged_TNBC_panCK_cnv, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Visvader_merged_TNBC_panCK_cnv)
Visvader_merged_TNBC_panCK_cnv <- ScaleData(Visvader_merged_TNBC_panCK_cnv, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
Visvader_merged_TNBC_panCK_cnv <- RunPCA(Visvader_merged_TNBC_panCK_cnv, features = VariableFeatures(object = Visvader_merged_TNBC_panCK_cnv))

# Examine and visualize PCA results a few different ways
print(Visvader_merged_TNBC_panCK_cnv[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Visvader_merged_TNBC_panCK_cnv, dims = 1:2, reduction = "pca")
DimPlot(Visvader_merged_TNBC_panCK_cnv, reduction = "pca")

# Visualize PCs
ElbowPlot(Visvader_merged_TNBC_panCK_cnv)

# Choose dimensions
Visvader_merged_TNBC_panCK_cnv <- FindNeighbors(Visvader_merged_TNBC_panCK_cnv, dims = 1:40)
Visvader_merged_TNBC_panCK_cnv <- FindClusters(Visvader_merged_TNBC_panCK_cnv, resolution = 0.12)

# Umap clustering
Visvader_merged_TNBC_panCK_cnv <- RunUMAP(Visvader_merged_TNBC_panCK_cnv, dims = 1:40)
DimPlot(Visvader_merged_TNBC_panCK_cnv, reduction = "umap")

# Identify Samples
DimPlot(Visvader_merged_TNBC_panCK_cnv, reduction = "umap", group.by = "orig.ident")

# Save
saveRDS(Visvader_merged_TNBC_panCK_cnv, file = "/R/R_Visvader/Visvader_RDS/RDS_Merged/Visvader_merged_TNBC_panCK_cnv.rds")

# Clear global environment
rm(Visvader_merged_TNBC_panCK_cnv)
gc()



##################################################################
# Merge the panCK-subsetted files into one for combined analysis #
##################################################################

########
# LOAD #
########
Visvader_merged_HR_panCK_CNV <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Merged/Visvader_merged_HR_panCK_CNV.rds")
Visvader_merged_HER2_panCK_CNV <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Merged/Visvader_merged_HER2_panCK_CNV.rds")
Visvader_merged_TNBC_panCK_cnv <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_Merged/Visvader_merged_TNBC_panCK_cnv.rds")


##############################################
# Merge Individual Files into Single Dataset #
##############################################
# Merge Datasets & remove intermediate files
all_Visvader_breast_tumors_panCK_cnv <- merge(x = Visvader_merged_HR_panCK_CNV, y = c(Visvader_merged_HER2_panCK_CNV, Visvader_merged_TNBC_panCK_cnv))

rm(Visvader_merged_TNBC_panCK_cnv)
rm(Visvader_merged_HER2_panCK_CNV)
rm(Visvader_merged_HR_panCK_CNV)
gc()

# Join Layers
all_Visvader_breast_tumors_panCK_cnv <- JoinLayers(all_Visvader_breast_tumors_panCK_cnv)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
all_Visvader_breast_tumors_panCK_cnv <- FindVariableFeatures(all_Visvader_breast_tumors_panCK_cnv, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(all_Visvader_breast_tumors_panCK_cnv)
all_Visvader_breast_tumors_panCK_cnv <- ScaleData(all_Visvader_breast_tumors_panCK_cnv, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################reduction
all_Visvader_breast_tumors_panCK_cnv <- RunPCA(all_Visvader_breast_tumors_panCK_cnv, features = VariableFeatures(object = all_Visvader_breast_tumors_panCK_cnv))

# Examine and visualize PCA results a few different ways
print(all_Visvader_breast_tumors_panCK_cnv[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(all_Visvader_breast_tumors_panCK_cnv, dims = 1:2, reduction = "pca")
DimPlot(all_Visvader_breast_tumors_panCK_cnv, reduction = "pca")

# Visualize PCs
ElbowPlot(all_Visvader_breast_tumors_panCK_cnv, ndims = 40)

all_Visvader_breast_tumors_panCK_cnv <- FindNeighbors(all_Visvader_breast_tumors_panCK_cnv, dims = 1:43)
all_Visvader_breast_tumors_panCK_cnv <- FindClusters(all_Visvader_breast_tumors_panCK_cnv, resolution = 0.2)

# Umap clustering
all_Visvader_breast_tumors_panCK_cnv <- RunUMAP(all_Visvader_breast_tumors_panCK_cnv, dims = 1:43)

DimPlot(all_Visvader_breast_tumors_panCK_cnv, reduction = "umap", raster = FALSE, label = TRUE)
DimPlot(all_Visvader_breast_tumors_panCK_cnv, reduction = "umap", group.by = "orig.ident", raster = FALSE)


########
# SAVE #
########
saveRDS(all_Visvader_breast_tumors_panCK_cnv, file = "/R/R_Visvader/Visvader_RDS/RDS_Merged/all_Visvader_breast_tumors_panCK_cnv.rds")


