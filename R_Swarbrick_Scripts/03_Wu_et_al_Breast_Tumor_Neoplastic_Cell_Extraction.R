#####################################################################################################################
#                       Swarbrick (Wu et al) Dataset Breast Tumor Analysis Steps                                    #
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

##############################################
#### InferCNV: Swarbrick_CID3941_ER_panCK ####
##############################################

########## Merge with normal mammary glad, extract epithelium, and process via Seurat
# Load RDS
Swarbrick_CID3941_ER_panCK <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_panCK/Swarbrick_CID3941_ER_panCK.rds")
Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Navin_Visvader_NORM_panCK_Reference.rds")

# Set identities 
colnames(x = Swarbrick_CID3941_ER_panCK[[]])
Idents(object = Swarbrick_CID3941_ER_panCK) <- 'Tumor'
Swarbrick_CID3941_ER_panCK[["cnv.ident"]] <- Idents(object = Swarbrick_CID3941_ER_panCK)
levels(x = Swarbrick_CID3941_ER_panCK)

colnames(x = Navin_Visvader_NORM_panCK_Reference[[]])
Idents(object = Navin_Visvader_NORM_panCK_Reference) <- 'Normal'
Navin_Visvader_NORM_panCK_Reference[["cnv.ident"]] <- Idents(object = Navin_Visvader_NORM_panCK_Reference)
levels(x = Navin_Visvader_NORM_panCK_Reference)

# Merge subsets into recombined RDS
Swarbrick_CID3941_ER_Navin_Visvader_NORM_panCK_Reference <- merge(x = Swarbrick_CID3941_ER_panCK, y = Navin_Visvader_NORM_panCK_Reference)
# Join Layers
Swarbrick_CID3941_ER_Navin_Visvader_NORM_panCK_Reference <- JoinLayers(Swarbrick_CID3941_ER_Navin_Visvader_NORM_panCK_Reference)

# Remove subsetted objects
rm(Swarbrick_CID3941_ER_panCK)
rm(Navin_Visvader_NORM_panCK_Reference)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples before merge

# Find Variable Features
Swarbrick_CID3941_ER_Navin_Visvader_NORM_panCK_Reference <- FindVariableFeatures(Swarbrick_CID3941_ER_Navin_Visvader_NORM_panCK_Reference, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Swarbrick_CID3941_ER_Navin_Visvader_NORM_panCK_Reference)
Swarbrick_CID3941_ER_Navin_Visvader_NORM_panCK_Reference <- ScaleData(Swarbrick_CID3941_ER_Navin_Visvader_NORM_panCK_Reference, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################

Swarbrick_CID3941_ER_Navin_Visvader_NORM_panCK_Reference <- RunPCA(Swarbrick_CID3941_ER_Navin_Visvader_NORM_panCK_Reference, verbose = FALSE)

# Examine and visualize PCA results a few different ways
print(Swarbrick_CID3941_ER_Navin_Visvader_NORM_panCK_Reference[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Swarbrick_CID3941_ER_Navin_Visvader_NORM_panCK_Reference, dims = 1:2, reduction = "pca")
DimPlot(Swarbrick_CID3941_ER_Navin_Visvader_NORM_panCK_Reference, reduction = "pca")

Swarbrick_CID3941_ER_Navin_Visvader_NORM_panCK_Reference <- FindNeighbors(Swarbrick_CID3941_ER_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
Swarbrick_CID3941_ER_Navin_Visvader_NORM_panCK_Reference <- FindClusters(Swarbrick_CID3941_ER_Navin_Visvader_NORM_panCK_Reference, resolution = 0.1)

# Umap clustering
Swarbrick_CID3941_ER_Navin_Visvader_NORM_panCK_Reference <- RunUMAP(Swarbrick_CID3941_ER_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
DimPlot(Swarbrick_CID3941_ER_Navin_Visvader_NORM_panCK_Reference, reduction = "umap")
DimPlot(Swarbrick_CID3941_ER_Navin_Visvader_NORM_panCK_Reference, reduction = "umap", group.by = "orig.ident")

# Save RDS
saveRDS(Swarbrick_CID3941_ER_Navin_Visvader_NORM_panCK_Reference, file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_inferCNV_Input_RDS/Swarbrick_CID3941_ER_Navin_Visvader_NORM_panCK_Reference.rds")

########## Prepare Cell Annotations
# Extract & Save Meta Data
Swarbrick_CID3941_ER_Navin_Visvader_NORM_panCK_Reference.meta.data <- as.data.frame(as.matrix(Swarbrick_CID3941_ER_Navin_Visvader_NORM_panCK_Reference@meta.data))
# Extract pertinent columns
# Make row names first column
Swarbrick_CID3941_ER_Navin_Visvader_NORM_panCK_Reference.meta.data <- tibble::rownames_to_column(Swarbrick_CID3941_ER_Navin_Visvader_NORM_panCK_Reference.meta.data, "Cells")
# In Seurat, "X" becomes the cell names while "cnv.ident" is specified above
Swarbrick_CID3941_ER_Navin_Visvader_NORM_panCK_Reference.annotations <- Swarbrick_CID3941_ER_Navin_Visvader_NORM_panCK_Reference.meta.data[, c("Cells", "cnv.ident")]
# Delete header
names(Swarbrick_CID3941_ER_Navin_Visvader_NORM_panCK_Reference.annotations) <- NULL
# Save as table
write.table(Swarbrick_CID3941_ER_Navin_Visvader_NORM_panCK_Reference.annotations,"/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_CID3941_ER_Navin_Visvader_NORM_panCK_Reference.annotations.txt",sep="\t",row.names=FALSE)
# Remove metadata and annotations
rm(Swarbrick_CID3941_ER_Navin_Visvader_NORM_panCK_Reference.meta.data)
rm(Swarbrick_CID3941_ER_Navin_Visvader_NORM_panCK_Reference.annotations)

########## Generate Count Matrix
# Extract Count Matrix
Swarbrick_CID3941_ER_Navin_Visvader_NORM_panCK_Reference.count <- Swarbrick_CID3941_ER_Navin_Visvader_NORM_panCK_Reference[["RNA"]]$counts
# Convert counts to a sparse matrix
Swarbrick_CID3941_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as(Swarbrick_CID3941_ER_Navin_Visvader_NORM_panCK_Reference.count, 'sparseMatrix')
rm(Swarbrick_CID3941_ER_Navin_Visvader_NORM_panCK_Reference.count)
# Convert to requisite data format
Swarbrick_CID3941_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as.data.frame(Swarbrick_CID3941_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix)
# Save table
write.csv(Swarbrick_CID3941_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix, file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_inferCNV_Input_Count_Matrix/Swarbrick_CID3941_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv")
# Remove Objects
rm(Swarbrick_CID3941_ER_Navin_Visvader_NORM_panCK_Reference)
rm(Swarbrick_CID3941_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix)

########### Run inferCNV
# Note: check.names = FALSE prevents cell names from having dashes replaced with dots; essential to match annotations file
Swarbrick_CID3941_ER_inferCNV <- CreateInfercnvObject(raw_counts_matrix=read.csv(file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_inferCNV_Input_Count_Matrix/Swarbrick_CID3941_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv", header = TRUE, row.names = 1,check.names = FALSE),
                                            annotations_file=read.table(file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_CID3941_ER_Navin_Visvader_NORM_panCK_Reference.annotations.txt", header = FALSE, row.names = 1),
                                            delim="\t",
                                            gene_order_file="/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/hg38_gencode_v27.txt",
                                            ref_group_names=c("Normal"))


# perform infercnv operations to reveal cnv signal
Swarbrick_CID3941_ER_inferCNV = infercnv::run(Swarbrick_CID3941_ER_inferCNV,
                                    cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                    out_dir="/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/Swarbrick_CID3941_ER_panCK/",  # dir is auto-created for storing outputs
                                    cluster_by_groups=T,   # cluster
                                    denoise=T,
                                    HMM=T,
                                    num_ref_groups = 9, analysis_mode = "subclusters")

###########  Integrate data with Seurat object
# Read RDS and add inferCNV results to metadata
Swarbrick_CID3941_ER_Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_inferCNV_Input_RDS/Swarbrick_CID3941_ER_Navin_Visvader_NORM_panCK_Reference.rds") 
Swarbrick_CID3941_ER_panCK_cnv <- infercnv::add_to_seurat(infercnv_output_path="/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/Swarbrick_CID3941_ER_panCK/",
                                                seurat_obj=Swarbrick_CID3941_ER_Navin_Visvader_NORM_panCK_Reference, # optional
                                                top_n=10)

# Collect Tumor cells from Seurat Object
Swarbrick_CID3941_ER_panCK_cnv <- subset(x = Swarbrick_CID3941_ER_panCK_cnv, subset = cnv.ident == "Tumor")
# Quantify chromosomes with alterations per cell
Swarbrick_CID3941_ER_panCK_cnv_Pos_Chromosomes <- as.data.frame(rowSums(Swarbrick_CID3941_ER_panCK_cnv@meta.data[,startsWith(names(Swarbrick_CID3941_ER_panCK_cnv@meta.data),"has_cnv_")] > 0 ))
Swarbrick_CID3941_ER_panCK_cnv_Neg_Chromosomes <- as.data.frame(rowSums(Swarbrick_CID3941_ER_panCK_cnv@meta.data[,startsWith(names(Swarbrick_CID3941_ER_panCK_cnv@meta.data),"has_cnv_")] == 0))
Swarbrick_CID3941_ER_panCK_cnv_Chromosome_Quant <- as.data.frame(c(Swarbrick_CID3941_ER_panCK_cnv_Pos_Chromosomes, Swarbrick_CID3941_ER_panCK_cnv_Neg_Chromosomes))
colnames(Swarbrick_CID3941_ER_panCK_cnv_Chromosome_Quant) <- c("Chromosomes CNV Pos", "Chromosomes CNV Neg")
write.csv(Swarbrick_CID3941_ER_panCK_cnv_Chromosome_Quant, file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/Swarbrick_CID3941_ER_panCK/Swarbrick_CID3941_ER_panCK_cnv_Chromosome_Quant.csv")

# Save CNV Metadata
Swarbrick_CID3941_ER_panCK_cnv.meta.data <- as.data.frame(as.matrix(Swarbrick_CID3941_ER_panCK_cnv@meta.data))
write.csv(Swarbrick_CID3941_ER_panCK_cnv.meta.data, file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/Swarbrick_CID3941_ER_panCK/Swarbrick_CID3941_ER_panCK_cnv.meta.data.csv")

# Subset Cells with inferred CNVs -------------------------------------------------------------------------------------------------------------
Swarbrick_CID3941_ER_panCK_cnv <- subset(Swarbrick_CID3941_ER_panCK_cnv, cells=rownames(Swarbrick_CID3941_ER_panCK_cnv@meta.data)[rowSums(Swarbrick_CID3941_ER_panCK_cnv@meta.data[,startsWith(names(Swarbrick_CID3941_ER_panCK_cnv@meta.data),"has_cnv_")]) > 0 ])

# Save RDS with inferredCNVs
saveRDS(Swarbrick_CID3941_ER_panCK_cnv, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_inferCNV/Swarbrick_CID3941_ER_panCK_cnv.rds")


# Free memory
rm(Swarbrick_CID3941_ER_inferCNV)
rm(Swarbrick_CID3941_ER_Navin_Visvader_NORM_panCK_Reference)
rm(Swarbrick_CID3941_ER_panCK_cnv)
rm(Swarbrick_CID3941_ER_panCK_cnv_Pos_Chromosomes)
rm(Swarbrick_CID3941_ER_panCK_cnv_Neg_Chromosomes)
rm(Swarbrick_CID3941_ER_panCK_cnv_Chromosome_Quant)
rm(Swarbrick_CID3941_ER_panCK_cnv.meta.data)
gc()



##############################################
#### InferCNV: Swarbrick_CID3948_ER_panCK ####
##############################################

########## Merge with normal mammary glad, extract epithelium, and process via Seurat
# Load RDS
Swarbrick_CID3948_ER_panCK <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_panCK/Swarbrick_CID3948_ER_panCK.rds")
Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Navin_Visvader_NORM_panCK_Reference.rds")

# Set identities 
colnames(x = Swarbrick_CID3948_ER_panCK[[]])
Idents(object = Swarbrick_CID3948_ER_panCK) <- 'Tumor'
Swarbrick_CID3948_ER_panCK[["cnv.ident"]] <- Idents(object = Swarbrick_CID3948_ER_panCK)
levels(x = Swarbrick_CID3948_ER_panCK)

colnames(x = Navin_Visvader_NORM_panCK_Reference[[]])
Idents(object = Navin_Visvader_NORM_panCK_Reference) <- 'Normal'
Navin_Visvader_NORM_panCK_Reference[["cnv.ident"]] <- Idents(object = Navin_Visvader_NORM_panCK_Reference)
levels(x = Navin_Visvader_NORM_panCK_Reference)

# Merge subsets into recombined RDS
Swarbrick_CID3948_ER_Navin_Visvader_NORM_panCK_Reference <- merge(x = Swarbrick_CID3948_ER_panCK, y = Navin_Visvader_NORM_panCK_Reference)
# Join Layers
Swarbrick_CID3948_ER_Navin_Visvader_NORM_panCK_Reference <- JoinLayers(Swarbrick_CID3948_ER_Navin_Visvader_NORM_panCK_Reference)

# Remove subsetted objects
rm(Swarbrick_CID3948_ER_panCK)
rm(Navin_Visvader_NORM_panCK_Reference)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples before merge

# Find Variable Features
Swarbrick_CID3948_ER_Navin_Visvader_NORM_panCK_Reference <- FindVariableFeatures(Swarbrick_CID3948_ER_Navin_Visvader_NORM_panCK_Reference, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Swarbrick_CID3948_ER_Navin_Visvader_NORM_panCK_Reference)
Swarbrick_CID3948_ER_Navin_Visvader_NORM_panCK_Reference <- ScaleData(Swarbrick_CID3948_ER_Navin_Visvader_NORM_panCK_Reference, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################

Swarbrick_CID3948_ER_Navin_Visvader_NORM_panCK_Reference <- RunPCA(Swarbrick_CID3948_ER_Navin_Visvader_NORM_panCK_Reference, verbose = FALSE)

# Examine and visualize PCA results a few different ways
print(Swarbrick_CID3948_ER_Navin_Visvader_NORM_panCK_Reference[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Swarbrick_CID3948_ER_Navin_Visvader_NORM_panCK_Reference, dims = 1:2, reduction = "pca")
DimPlot(Swarbrick_CID3948_ER_Navin_Visvader_NORM_panCK_Reference, reduction = "pca")

Swarbrick_CID3948_ER_Navin_Visvader_NORM_panCK_Reference <- FindNeighbors(Swarbrick_CID3948_ER_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
Swarbrick_CID3948_ER_Navin_Visvader_NORM_panCK_Reference <- FindClusters(Swarbrick_CID3948_ER_Navin_Visvader_NORM_panCK_Reference, resolution = 0.1)

# Umap clustering
Swarbrick_CID3948_ER_Navin_Visvader_NORM_panCK_Reference <- RunUMAP(Swarbrick_CID3948_ER_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
DimPlot(Swarbrick_CID3948_ER_Navin_Visvader_NORM_panCK_Reference, reduction = "umap")
DimPlot(Swarbrick_CID3948_ER_Navin_Visvader_NORM_panCK_Reference, reduction = "umap", group.by = "orig.ident")

# Save RDS
saveRDS(Swarbrick_CID3948_ER_Navin_Visvader_NORM_panCK_Reference, file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_inferCNV_Input_RDS/Swarbrick_CID3948_ER_Navin_Visvader_NORM_panCK_Reference.rds")

########## Prepare Cell Annotations
# Extract & Save Meta Data
Swarbrick_CID3948_ER_Navin_Visvader_NORM_panCK_Reference.meta.data <- as.data.frame(as.matrix(Swarbrick_CID3948_ER_Navin_Visvader_NORM_panCK_Reference@meta.data))
# Extract pertinent columns
# Make row names first column
Swarbrick_CID3948_ER_Navin_Visvader_NORM_panCK_Reference.meta.data <- tibble::rownames_to_column(Swarbrick_CID3948_ER_Navin_Visvader_NORM_panCK_Reference.meta.data, "Cells")
# In Seurat, "X" becomes the cell names while "cnv.ident" is specified above
Swarbrick_CID3948_ER_Navin_Visvader_NORM_panCK_Reference.annotations <- Swarbrick_CID3948_ER_Navin_Visvader_NORM_panCK_Reference.meta.data[, c("Cells", "cnv.ident")]
# Delete header
names(Swarbrick_CID3948_ER_Navin_Visvader_NORM_panCK_Reference.annotations) <- NULL
# Save as table
write.table(Swarbrick_CID3948_ER_Navin_Visvader_NORM_panCK_Reference.annotations,"/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_CID3948_ER_Navin_Visvader_NORM_panCK_Reference.annotations.txt",sep="\t",row.names=FALSE)
# Remove metadata and annotations
rm(Swarbrick_CID3948_ER_Navin_Visvader_NORM_panCK_Reference.meta.data)
rm(Swarbrick_CID3948_ER_Navin_Visvader_NORM_panCK_Reference.annotations)

########## Generate Count Matrix
# Extract Count Matrix
Swarbrick_CID3948_ER_Navin_Visvader_NORM_panCK_Reference.count <- Swarbrick_CID3948_ER_Navin_Visvader_NORM_panCK_Reference[["RNA"]]$counts
# Convert counts to a sparse matrix
Swarbrick_CID3948_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as(Swarbrick_CID3948_ER_Navin_Visvader_NORM_panCK_Reference.count, 'sparseMatrix')
rm(Swarbrick_CID3948_ER_Navin_Visvader_NORM_panCK_Reference.count)
# Convert to requisite data format
Swarbrick_CID3948_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as.data.frame(Swarbrick_CID3948_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix)
# Save table
write.csv(Swarbrick_CID3948_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix, file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_inferCNV_Input_Count_Matrix/Swarbrick_CID3948_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv")
# Remove Objects
rm(Swarbrick_CID3948_ER_Navin_Visvader_NORM_panCK_Reference)
rm(Swarbrick_CID3948_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix)

########### Run inferCNV
# Note: check.names = FALSE prevents cell names from having dashes replaced with dots; essential to match annotations file
Swarbrick_CID3948_ER_inferCNV <- CreateInfercnvObject(raw_counts_matrix=read.csv(file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_inferCNV_Input_Count_Matrix/Swarbrick_CID3948_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv", header = TRUE, row.names = 1,check.names = FALSE),
                                            annotations_file=read.table(file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_CID3948_ER_Navin_Visvader_NORM_panCK_Reference.annotations.txt", header = FALSE, row.names = 1),
                                            delim="\t",
                                            gene_order_file="/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/hg38_gencode_v27.txt",
                                            ref_group_names=c("Normal"))


# perform infercnv operations to reveal cnv signal
Swarbrick_CID3948_ER_inferCNV = infercnv::run(Swarbrick_CID3948_ER_inferCNV,
                                    cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                    out_dir="/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/Swarbrick_CID3948_ER_panCK/",  # dir is auto-created for storing outputs
                                    cluster_by_groups=T,   # cluster
                                    denoise=T,
                                    HMM=T,
                                    num_ref_groups = 9, analysis_mode = "subclusters")

###########  Integrate data with Seurat object
# Read RDS and add inferCNV results to metadata
Swarbrick_CID3948_ER_Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_inferCNV_Input_RDS/Swarbrick_CID3948_ER_Navin_Visvader_NORM_panCK_Reference.rds") 
Swarbrick_CID3948_ER_panCK_cnv <- infercnv::add_to_seurat(infercnv_output_path="/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/Swarbrick_CID3948_ER_panCK/",
                                                seurat_obj=Swarbrick_CID3948_ER_Navin_Visvader_NORM_panCK_Reference, # optional
                                                top_n=10)

# Collect Tumor cells from Seurat Object
Swarbrick_CID3948_ER_panCK_cnv <- subset(x = Swarbrick_CID3948_ER_panCK_cnv, subset = cnv.ident == "Tumor")
# Quantify chromosomes with alterations per cell
Swarbrick_CID3948_ER_panCK_cnv_Pos_Chromosomes <- as.data.frame(rowSums(Swarbrick_CID3948_ER_panCK_cnv@meta.data[,startsWith(names(Swarbrick_CID3948_ER_panCK_cnv@meta.data),"has_cnv_")] > 0 ))
Swarbrick_CID3948_ER_panCK_cnv_Neg_Chromosomes <- as.data.frame(rowSums(Swarbrick_CID3948_ER_panCK_cnv@meta.data[,startsWith(names(Swarbrick_CID3948_ER_panCK_cnv@meta.data),"has_cnv_")] == 0))
Swarbrick_CID3948_ER_panCK_cnv_Chromosome_Quant <- as.data.frame(c(Swarbrick_CID3948_ER_panCK_cnv_Pos_Chromosomes, Swarbrick_CID3948_ER_panCK_cnv_Neg_Chromosomes))
colnames(Swarbrick_CID3948_ER_panCK_cnv_Chromosome_Quant) <- c("Chromosomes CNV Pos", "Chromosomes CNV Neg")
write.csv(Swarbrick_CID3948_ER_panCK_cnv_Chromosome_Quant, file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/Swarbrick_CID3948_ER_panCK/Swarbrick_CID3948_ER_panCK_cnv_Chromosome_Quant.csv")

# Save CNV Metadata
Swarbrick_CID3948_ER_panCK_cnv.meta.data <- as.data.frame(as.matrix(Swarbrick_CID3948_ER_panCK_cnv@meta.data))
write.csv(Swarbrick_CID3948_ER_panCK_cnv.meta.data, file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/Swarbrick_CID3948_ER_panCK/Swarbrick_CID3948_ER_panCK_cnv.meta.data.csv")

# Subset Cells with inferred CNVs -------------------------------------------------------------------------------------------------------------
Swarbrick_CID3948_ER_panCK_cnv <- subset(Swarbrick_CID3948_ER_panCK_cnv, cells=rownames(Swarbrick_CID3948_ER_panCK_cnv@meta.data)[rowSums(Swarbrick_CID3948_ER_panCK_cnv@meta.data[,startsWith(names(Swarbrick_CID3948_ER_panCK_cnv@meta.data),"has_cnv_")]) > 0 ])

# Save RDS with inferredCNVs
saveRDS(Swarbrick_CID3948_ER_panCK_cnv, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_inferCNV/Swarbrick_CID3948_ER_panCK_cnv.rds")


# Free memory
rm(Swarbrick_CID3948_ER_inferCNV)
rm(Swarbrick_CID3948_ER_Navin_Visvader_NORM_panCK_Reference)
rm(Swarbrick_CID3948_ER_panCK_cnv)
rm(Swarbrick_CID3948_ER_panCK_cnv_Pos_Chromosomes)
rm(Swarbrick_CID3948_ER_panCK_cnv_Neg_Chromosomes)
rm(Swarbrick_CID3948_ER_panCK_cnv_Chromosome_Quant)
rm(Swarbrick_CID3948_ER_panCK_cnv.meta.data)
gc()



##############################################
#### InferCNV: Swarbrick_CID4040_ER_panCK ####
##############################################

########## Merge with normal mammary glad, extract epithelium, and process via Seurat
# Load RDS
Swarbrick_CID4040_ER_panCK <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_panCK/Swarbrick_CID4040_ER_panCK.rds")
Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Navin_Visvader_NORM_panCK_Reference.rds")

# Set identities 
colnames(x = Swarbrick_CID4040_ER_panCK[[]])
Idents(object = Swarbrick_CID4040_ER_panCK) <- 'Tumor'
Swarbrick_CID4040_ER_panCK[["cnv.ident"]] <- Idents(object = Swarbrick_CID4040_ER_panCK)
levels(x = Swarbrick_CID4040_ER_panCK)

colnames(x = Navin_Visvader_NORM_panCK_Reference[[]])
Idents(object = Navin_Visvader_NORM_panCK_Reference) <- 'Normal'
Navin_Visvader_NORM_panCK_Reference[["cnv.ident"]] <- Idents(object = Navin_Visvader_NORM_panCK_Reference)
levels(x = Navin_Visvader_NORM_panCK_Reference)

# Merge subsets into recombined RDS
Swarbrick_CID4040_ER_Navin_Visvader_NORM_panCK_Reference <- merge(x = Swarbrick_CID4040_ER_panCK, y = Navin_Visvader_NORM_panCK_Reference)
# Join Layers
Swarbrick_CID4040_ER_Navin_Visvader_NORM_panCK_Reference <- JoinLayers(Swarbrick_CID4040_ER_Navin_Visvader_NORM_panCK_Reference)

# Remove subsetted objects
rm(Swarbrick_CID4040_ER_panCK)
rm(Navin_Visvader_NORM_panCK_Reference)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples before merge

# Find Variable Features
Swarbrick_CID4040_ER_Navin_Visvader_NORM_panCK_Reference <- FindVariableFeatures(Swarbrick_CID4040_ER_Navin_Visvader_NORM_panCK_Reference, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Swarbrick_CID4040_ER_Navin_Visvader_NORM_panCK_Reference)
Swarbrick_CID4040_ER_Navin_Visvader_NORM_panCK_Reference <- ScaleData(Swarbrick_CID4040_ER_Navin_Visvader_NORM_panCK_Reference, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################

Swarbrick_CID4040_ER_Navin_Visvader_NORM_panCK_Reference <- RunPCA(Swarbrick_CID4040_ER_Navin_Visvader_NORM_panCK_Reference, verbose = FALSE)

# Examine and visualize PCA results a few different ways
print(Swarbrick_CID4040_ER_Navin_Visvader_NORM_panCK_Reference[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Swarbrick_CID4040_ER_Navin_Visvader_NORM_panCK_Reference, dims = 1:2, reduction = "pca")
DimPlot(Swarbrick_CID4040_ER_Navin_Visvader_NORM_panCK_Reference, reduction = "pca")

Swarbrick_CID4040_ER_Navin_Visvader_NORM_panCK_Reference <- FindNeighbors(Swarbrick_CID4040_ER_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
Swarbrick_CID4040_ER_Navin_Visvader_NORM_panCK_Reference <- FindClusters(Swarbrick_CID4040_ER_Navin_Visvader_NORM_panCK_Reference, resolution = 0.1)

# Umap clustering
Swarbrick_CID4040_ER_Navin_Visvader_NORM_panCK_Reference <- RunUMAP(Swarbrick_CID4040_ER_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
DimPlot(Swarbrick_CID4040_ER_Navin_Visvader_NORM_panCK_Reference, reduction = "umap")
DimPlot(Swarbrick_CID4040_ER_Navin_Visvader_NORM_panCK_Reference, reduction = "umap", group.by = "orig.ident")

# Save RDS
saveRDS(Swarbrick_CID4040_ER_Navin_Visvader_NORM_panCK_Reference, file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_inferCNV_Input_RDS/Swarbrick_CID4040_ER_Navin_Visvader_NORM_panCK_Reference.rds")

########## Prepare Cell Annotations
# Extract & Save Meta Data
Swarbrick_CID4040_ER_Navin_Visvader_NORM_panCK_Reference.meta.data <- as.data.frame(as.matrix(Swarbrick_CID4040_ER_Navin_Visvader_NORM_panCK_Reference@meta.data))
# Extract pertinent columns
# Make row names first column
Swarbrick_CID4040_ER_Navin_Visvader_NORM_panCK_Reference.meta.data <- tibble::rownames_to_column(Swarbrick_CID4040_ER_Navin_Visvader_NORM_panCK_Reference.meta.data, "Cells")
# In Seurat, "X" becomes the cell names while "cnv.ident" is specified above
Swarbrick_CID4040_ER_Navin_Visvader_NORM_panCK_Reference.annotations <- Swarbrick_CID4040_ER_Navin_Visvader_NORM_panCK_Reference.meta.data[, c("Cells", "cnv.ident")]
# Delete header
names(Swarbrick_CID4040_ER_Navin_Visvader_NORM_panCK_Reference.annotations) <- NULL
# Save as table
write.table(Swarbrick_CID4040_ER_Navin_Visvader_NORM_panCK_Reference.annotations,"/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_CID4040_ER_Navin_Visvader_NORM_panCK_Reference.annotations.txt",sep="\t",row.names=FALSE)
# Remove metadata and annotations
rm(Swarbrick_CID4040_ER_Navin_Visvader_NORM_panCK_Reference.meta.data)
rm(Swarbrick_CID4040_ER_Navin_Visvader_NORM_panCK_Reference.annotations)

########## Generate Count Matrix
# Extract Count Matrix
Swarbrick_CID4040_ER_Navin_Visvader_NORM_panCK_Reference.count <- Swarbrick_CID4040_ER_Navin_Visvader_NORM_panCK_Reference[["RNA"]]$counts
# Convert counts to a sparse matrix
Swarbrick_CID4040_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as(Swarbrick_CID4040_ER_Navin_Visvader_NORM_panCK_Reference.count, 'sparseMatrix')
rm(Swarbrick_CID4040_ER_Navin_Visvader_NORM_panCK_Reference.count)
# Convert to requisite data format
Swarbrick_CID4040_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as.data.frame(Swarbrick_CID4040_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix)
# Save table
write.csv(Swarbrick_CID4040_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix, file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_inferCNV_Input_Count_Matrix/Swarbrick_CID4040_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv")
# Remove Objects
rm(Swarbrick_CID4040_ER_Navin_Visvader_NORM_panCK_Reference)
rm(Swarbrick_CID4040_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix)

########### Run inferCNV
# Note: check.names = FALSE prevents cell names from having dashes replaced with dots; essential to match annotations file
Swarbrick_CID4040_ER_inferCNV <- CreateInfercnvObject(raw_counts_matrix=read.csv(file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_inferCNV_Input_Count_Matrix/Swarbrick_CID4040_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv", header = TRUE, row.names = 1,check.names = FALSE),
                                            annotations_file=read.table(file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_CID4040_ER_Navin_Visvader_NORM_panCK_Reference.annotations.txt", header = FALSE, row.names = 1),
                                            delim="\t",
                                            gene_order_file="/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/hg38_gencode_v27.txt",
                                            ref_group_names=c("Normal"))


# perform infercnv operations to reveal cnv signal
Swarbrick_CID4040_ER_inferCNV = infercnv::run(Swarbrick_CID4040_ER_inferCNV,
                                    cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                    out_dir="/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/Swarbrick_CID4040_ER_panCK/",  # dir is auto-created for storing outputs
                                    cluster_by_groups=T,   # cluster
                                    denoise=T,
                                    HMM=T,
                                    num_ref_groups = 9, analysis_mode = "subclusters")

###########  Integrate data with Seurat object
# Read RDS and add inferCNV results to metadata
Swarbrick_CID4040_ER_Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_inferCNV_Input_RDS/Swarbrick_CID4040_ER_Navin_Visvader_NORM_panCK_Reference.rds") 
Swarbrick_CID4040_ER_panCK_cnv <- infercnv::add_to_seurat(infercnv_output_path="/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/Swarbrick_CID4040_ER_panCK/",
                                                seurat_obj=Swarbrick_CID4040_ER_Navin_Visvader_NORM_panCK_Reference, # optional
                                                top_n=10)

# Collect Tumor cells from Seurat Object
Swarbrick_CID4040_ER_panCK_cnv <- subset(x = Swarbrick_CID4040_ER_panCK_cnv, subset = cnv.ident == "Tumor")
# Quantify chromosomes with alterations per cell
Swarbrick_CID4040_ER_panCK_cnv_Pos_Chromosomes <- as.data.frame(rowSums(Swarbrick_CID4040_ER_panCK_cnv@meta.data[,startsWith(names(Swarbrick_CID4040_ER_panCK_cnv@meta.data),"has_cnv_")] > 0 ))
Swarbrick_CID4040_ER_panCK_cnv_Neg_Chromosomes <- as.data.frame(rowSums(Swarbrick_CID4040_ER_panCK_cnv@meta.data[,startsWith(names(Swarbrick_CID4040_ER_panCK_cnv@meta.data),"has_cnv_")] == 0))
Swarbrick_CID4040_ER_panCK_cnv_Chromosome_Quant <- as.data.frame(c(Swarbrick_CID4040_ER_panCK_cnv_Pos_Chromosomes, Swarbrick_CID4040_ER_panCK_cnv_Neg_Chromosomes))
colnames(Swarbrick_CID4040_ER_panCK_cnv_Chromosome_Quant) <- c("Chromosomes CNV Pos", "Chromosomes CNV Neg")
write.csv(Swarbrick_CID4040_ER_panCK_cnv_Chromosome_Quant, file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/Swarbrick_CID4040_ER_panCK/Swarbrick_CID4040_ER_panCK_cnv_Chromosome_Quant.csv")

# Save CNV Metadata
Swarbrick_CID4040_ER_panCK_cnv.meta.data <- as.data.frame(as.matrix(Swarbrick_CID4040_ER_panCK_cnv@meta.data))
write.csv(Swarbrick_CID4040_ER_panCK_cnv.meta.data, file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/Swarbrick_CID4040_ER_panCK/Swarbrick_CID4040_ER_panCK_cnv.meta.data.csv")

# Subset Cells with inferred CNVs -------------------------------------------------------------------------------------------------------------
Swarbrick_CID4040_ER_panCK_cnv <- subset(Swarbrick_CID4040_ER_panCK_cnv, cells=rownames(Swarbrick_CID4040_ER_panCK_cnv@meta.data)[rowSums(Swarbrick_CID4040_ER_panCK_cnv@meta.data[,startsWith(names(Swarbrick_CID4040_ER_panCK_cnv@meta.data),"has_cnv_")]) > 0 ])

# Save RDS with inferredCNVs
saveRDS(Swarbrick_CID4040_ER_panCK_cnv, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_inferCNV/Swarbrick_CID4040_ER_panCK_cnv.rds")


# Free memory
rm(Swarbrick_CID4040_ER_inferCNV)
rm(Swarbrick_CID4040_ER_Navin_Visvader_NORM_panCK_Reference)
rm(Swarbrick_CID4040_ER_panCK_cnv)
rm(Swarbrick_CID4040_ER_panCK_cnv_Pos_Chromosomes)
rm(Swarbrick_CID4040_ER_panCK_cnv_Neg_Chromosomes)
rm(Swarbrick_CID4040_ER_panCK_cnv_Chromosome_Quant)
rm(Swarbrick_CID4040_ER_panCK_cnv.meta.data)
gc()



##############################################
#### InferCNV: Swarbrick_CID4067_ER_panCK ####
##############################################

########## Merge with normal mammary glad, extract epithelium, and process via Seurat
# Load RDS
Swarbrick_CID4067_ER_panCK <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_panCK/Swarbrick_CID4067_ER_panCK.rds")
Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Navin_Visvader_NORM_panCK_Reference.rds")

# Set identities 
colnames(x = Swarbrick_CID4067_ER_panCK[[]])
Idents(object = Swarbrick_CID4067_ER_panCK) <- 'Tumor'
Swarbrick_CID4067_ER_panCK[["cnv.ident"]] <- Idents(object = Swarbrick_CID4067_ER_panCK)
levels(x = Swarbrick_CID4067_ER_panCK)

colnames(x = Navin_Visvader_NORM_panCK_Reference[[]])
Idents(object = Navin_Visvader_NORM_panCK_Reference) <- 'Normal'
Navin_Visvader_NORM_panCK_Reference[["cnv.ident"]] <- Idents(object = Navin_Visvader_NORM_panCK_Reference)
levels(x = Navin_Visvader_NORM_panCK_Reference)

# Merge subsets into recombined RDS
Swarbrick_CID4067_ER_Navin_Visvader_NORM_panCK_Reference <- merge(x = Swarbrick_CID4067_ER_panCK, y = Navin_Visvader_NORM_panCK_Reference)
# Join Layers
Swarbrick_CID4067_ER_Navin_Visvader_NORM_panCK_Reference <- JoinLayers(Swarbrick_CID4067_ER_Navin_Visvader_NORM_panCK_Reference)

# Remove subsetted objects
rm(Swarbrick_CID4067_ER_panCK)
rm(Navin_Visvader_NORM_panCK_Reference)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples before merge

# Find Variable Features
Swarbrick_CID4067_ER_Navin_Visvader_NORM_panCK_Reference <- FindVariableFeatures(Swarbrick_CID4067_ER_Navin_Visvader_NORM_panCK_Reference, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Swarbrick_CID4067_ER_Navin_Visvader_NORM_panCK_Reference)
Swarbrick_CID4067_ER_Navin_Visvader_NORM_panCK_Reference <- ScaleData(Swarbrick_CID4067_ER_Navin_Visvader_NORM_panCK_Reference, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################

Swarbrick_CID4067_ER_Navin_Visvader_NORM_panCK_Reference <- RunPCA(Swarbrick_CID4067_ER_Navin_Visvader_NORM_panCK_Reference, verbose = FALSE)

# Examine and visualize PCA results a few different ways
print(Swarbrick_CID4067_ER_Navin_Visvader_NORM_panCK_Reference[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Swarbrick_CID4067_ER_Navin_Visvader_NORM_panCK_Reference, dims = 1:2, reduction = "pca")
DimPlot(Swarbrick_CID4067_ER_Navin_Visvader_NORM_panCK_Reference, reduction = "pca")

Swarbrick_CID4067_ER_Navin_Visvader_NORM_panCK_Reference <- FindNeighbors(Swarbrick_CID4067_ER_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
Swarbrick_CID4067_ER_Navin_Visvader_NORM_panCK_Reference <- FindClusters(Swarbrick_CID4067_ER_Navin_Visvader_NORM_panCK_Reference, resolution = 0.1)

# Umap clustering
Swarbrick_CID4067_ER_Navin_Visvader_NORM_panCK_Reference <- RunUMAP(Swarbrick_CID4067_ER_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
DimPlot(Swarbrick_CID4067_ER_Navin_Visvader_NORM_panCK_Reference, reduction = "umap")
DimPlot(Swarbrick_CID4067_ER_Navin_Visvader_NORM_panCK_Reference, reduction = "umap", group.by = "orig.ident")

# Save RDS
saveRDS(Swarbrick_CID4067_ER_Navin_Visvader_NORM_panCK_Reference, file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_inferCNV_Input_RDS/Swarbrick_CID4067_ER_Navin_Visvader_NORM_panCK_Reference.rds")

########## Prepare Cell Annotations
# Extract & Save Meta Data
Swarbrick_CID4067_ER_Navin_Visvader_NORM_panCK_Reference.meta.data <- as.data.frame(as.matrix(Swarbrick_CID4067_ER_Navin_Visvader_NORM_panCK_Reference@meta.data))
# Extract pertinent columns
# Make row names first column
Swarbrick_CID4067_ER_Navin_Visvader_NORM_panCK_Reference.meta.data <- tibble::rownames_to_column(Swarbrick_CID4067_ER_Navin_Visvader_NORM_panCK_Reference.meta.data, "Cells")
# In Seurat, "X" becomes the cell names while "cnv.ident" is specified above
Swarbrick_CID4067_ER_Navin_Visvader_NORM_panCK_Reference.annotations <- Swarbrick_CID4067_ER_Navin_Visvader_NORM_panCK_Reference.meta.data[, c("Cells", "cnv.ident")]
# Delete header
names(Swarbrick_CID4067_ER_Navin_Visvader_NORM_panCK_Reference.annotations) <- NULL
# Save as table
write.table(Swarbrick_CID4067_ER_Navin_Visvader_NORM_panCK_Reference.annotations,"/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_CID4067_ER_Navin_Visvader_NORM_panCK_Reference.annotations.txt",sep="\t",row.names=FALSE)
# Remove metadata and annotations
rm(Swarbrick_CID4067_ER_Navin_Visvader_NORM_panCK_Reference.meta.data)
rm(Swarbrick_CID4067_ER_Navin_Visvader_NORM_panCK_Reference.annotations)

########## Generate Count Matrix
# Extract Count Matrix
Swarbrick_CID4067_ER_Navin_Visvader_NORM_panCK_Reference.count <- Swarbrick_CID4067_ER_Navin_Visvader_NORM_panCK_Reference[["RNA"]]$counts
# Convert counts to a sparse matrix
Swarbrick_CID4067_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as(Swarbrick_CID4067_ER_Navin_Visvader_NORM_panCK_Reference.count, 'sparseMatrix')
rm(Swarbrick_CID4067_ER_Navin_Visvader_NORM_panCK_Reference.count)
# Convert to requisite data format
Swarbrick_CID4067_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as.data.frame(Swarbrick_CID4067_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix)
# Save table
write.csv(Swarbrick_CID4067_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix, file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_inferCNV_Input_Count_Matrix/Swarbrick_CID4067_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv")
# Remove Objects
rm(Swarbrick_CID4067_ER_Navin_Visvader_NORM_panCK_Reference)
rm(Swarbrick_CID4067_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix)

########### Run inferCNV
# Note: check.names = FALSE prevents cell names from having dashes replaced with dots; essential to match annotations file
Swarbrick_CID4067_ER_inferCNV <- CreateInfercnvObject(raw_counts_matrix=read.csv(file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_inferCNV_Input_Count_Matrix/Swarbrick_CID4067_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv", header = TRUE, row.names = 1,check.names = FALSE),
                                            annotations_file=read.table(file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_CID4067_ER_Navin_Visvader_NORM_panCK_Reference.annotations.txt", header = FALSE, row.names = 1),
                                            delim="\t",
                                            gene_order_file="/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/hg38_gencode_v27.txt",
                                            ref_group_names=c("Normal"))


# perform infercnv operations to reveal cnv signal
Swarbrick_CID4067_ER_inferCNV = infercnv::run(Swarbrick_CID4067_ER_inferCNV,
                                    cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                    out_dir="/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/Swarbrick_CID4067_ER_panCK/",  # dir is auto-created for storing outputs
                                    cluster_by_groups=T,   # cluster
                                    denoise=T,
                                    HMM=T,
                                    num_ref_groups = 9, analysis_mode = "subclusters")

###########  Integrate data with Seurat object
# Read RDS and add inferCNV results to metadata
Swarbrick_CID4067_ER_Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_inferCNV_Input_RDS/Swarbrick_CID4067_ER_Navin_Visvader_NORM_panCK_Reference.rds") 
Swarbrick_CID4067_ER_panCK_cnv <- infercnv::add_to_seurat(infercnv_output_path="/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/Swarbrick_CID4067_ER_panCK/",
                                                seurat_obj=Swarbrick_CID4067_ER_Navin_Visvader_NORM_panCK_Reference, # optional
                                                top_n=10)

# Collect Tumor cells from Seurat Object
Swarbrick_CID4067_ER_panCK_cnv <- subset(x = Swarbrick_CID4067_ER_panCK_cnv, subset = cnv.ident == "Tumor")
# Quantify chromosomes with alterations per cell
Swarbrick_CID4067_ER_panCK_cnv_Pos_Chromosomes <- as.data.frame(rowSums(Swarbrick_CID4067_ER_panCK_cnv@meta.data[,startsWith(names(Swarbrick_CID4067_ER_panCK_cnv@meta.data),"has_cnv_")] > 0 ))
Swarbrick_CID4067_ER_panCK_cnv_Neg_Chromosomes <- as.data.frame(rowSums(Swarbrick_CID4067_ER_panCK_cnv@meta.data[,startsWith(names(Swarbrick_CID4067_ER_panCK_cnv@meta.data),"has_cnv_")] == 0))
Swarbrick_CID4067_ER_panCK_cnv_Chromosome_Quant <- as.data.frame(c(Swarbrick_CID4067_ER_panCK_cnv_Pos_Chromosomes, Swarbrick_CID4067_ER_panCK_cnv_Neg_Chromosomes))
colnames(Swarbrick_CID4067_ER_panCK_cnv_Chromosome_Quant) <- c("Chromosomes CNV Pos", "Chromosomes CNV Neg")
write.csv(Swarbrick_CID4067_ER_panCK_cnv_Chromosome_Quant, file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/Swarbrick_CID4067_ER_panCK/Swarbrick_CID4067_ER_panCK_cnv_Chromosome_Quant.csv")

# Save CNV Metadata
Swarbrick_CID4067_ER_panCK_cnv.meta.data <- as.data.frame(as.matrix(Swarbrick_CID4067_ER_panCK_cnv@meta.data))
write.csv(Swarbrick_CID4067_ER_panCK_cnv.meta.data, file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/Swarbrick_CID4067_ER_panCK/Swarbrick_CID4067_ER_panCK_cnv.meta.data.csv")

# Subset Cells with inferred CNVs -------------------------------------------------------------------------------------------------------------
Swarbrick_CID4067_ER_panCK_cnv <- subset(Swarbrick_CID4067_ER_panCK_cnv, cells=rownames(Swarbrick_CID4067_ER_panCK_cnv@meta.data)[rowSums(Swarbrick_CID4067_ER_panCK_cnv@meta.data[,startsWith(names(Swarbrick_CID4067_ER_panCK_cnv@meta.data),"has_cnv_")]) > 0 ])

# Save RDS with inferredCNVs
saveRDS(Swarbrick_CID4067_ER_panCK_cnv, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_inferCNV/Swarbrick_CID4067_ER_panCK_cnv.rds")


# Free memory
rm(Swarbrick_CID4067_ER_inferCNV)
rm(Swarbrick_CID4067_ER_Navin_Visvader_NORM_panCK_Reference)
rm(Swarbrick_CID4067_ER_panCK_cnv)
rm(Swarbrick_CID4067_ER_panCK_cnv_Pos_Chromosomes)
rm(Swarbrick_CID4067_ER_panCK_cnv_Neg_Chromosomes)
rm(Swarbrick_CID4067_ER_panCK_cnv_Chromosome_Quant)
rm(Swarbrick_CID4067_ER_panCK_cnv.meta.data)
gc()



###############################################
#### InferCNV: Swarbrick_CID4290A_ER_panCK ####
###############################################

########## Merge with normal mammary glad, extract epithelium, and process via Seurat
# Load RDS
Swarbrick_CID4290A_ER_panCK <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_panCK/Swarbrick_CID4290A_ER_panCK.rds")
Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Navin_Visvader_NORM_panCK_Reference.rds")

# Set identities 
colnames(x = Swarbrick_CID4290A_ER_panCK[[]])
Idents(object = Swarbrick_CID4290A_ER_panCK) <- 'Tumor'
Swarbrick_CID4290A_ER_panCK[["cnv.ident"]] <- Idents(object = Swarbrick_CID4290A_ER_panCK)
levels(x = Swarbrick_CID4290A_ER_panCK)

colnames(x = Navin_Visvader_NORM_panCK_Reference[[]])
Idents(object = Navin_Visvader_NORM_panCK_Reference) <- 'Normal'
Navin_Visvader_NORM_panCK_Reference[["cnv.ident"]] <- Idents(object = Navin_Visvader_NORM_panCK_Reference)
levels(x = Navin_Visvader_NORM_panCK_Reference)

# Merge subsets into recombined RDS
Swarbrick_CID4290A_ER_Navin_Visvader_NORM_panCK_Reference <- merge(x = Swarbrick_CID4290A_ER_panCK, y = Navin_Visvader_NORM_panCK_Reference)
# Join Layers
Swarbrick_CID4290A_ER_Navin_Visvader_NORM_panCK_Reference <- JoinLayers(Swarbrick_CID4290A_ER_Navin_Visvader_NORM_panCK_Reference)

# Remove subsetted objects
rm(Swarbrick_CID4290A_ER_panCK)
rm(Navin_Visvader_NORM_panCK_Reference)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples before merge

# Find Variable Features
Swarbrick_CID4290A_ER_Navin_Visvader_NORM_panCK_Reference <- FindVariableFeatures(Swarbrick_CID4290A_ER_Navin_Visvader_NORM_panCK_Reference, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Swarbrick_CID4290A_ER_Navin_Visvader_NORM_panCK_Reference)
Swarbrick_CID4290A_ER_Navin_Visvader_NORM_panCK_Reference <- ScaleData(Swarbrick_CID4290A_ER_Navin_Visvader_NORM_panCK_Reference, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################

Swarbrick_CID4290A_ER_Navin_Visvader_NORM_panCK_Reference <- RunPCA(Swarbrick_CID4290A_ER_Navin_Visvader_NORM_panCK_Reference, verbose = FALSE)

# Examine and visualize PCA results a few different ways
print(Swarbrick_CID4290A_ER_Navin_Visvader_NORM_panCK_Reference[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Swarbrick_CID4290A_ER_Navin_Visvader_NORM_panCK_Reference, dims = 1:2, reduction = "pca")
DimPlot(Swarbrick_CID4290A_ER_Navin_Visvader_NORM_panCK_Reference, reduction = "pca")

Swarbrick_CID4290A_ER_Navin_Visvader_NORM_panCK_Reference <- FindNeighbors(Swarbrick_CID4290A_ER_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
Swarbrick_CID4290A_ER_Navin_Visvader_NORM_panCK_Reference <- FindClusters(Swarbrick_CID4290A_ER_Navin_Visvader_NORM_panCK_Reference, resolution = 0.1)

# Umap clustering
Swarbrick_CID4290A_ER_Navin_Visvader_NORM_panCK_Reference <- RunUMAP(Swarbrick_CID4290A_ER_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
DimPlot(Swarbrick_CID4290A_ER_Navin_Visvader_NORM_panCK_Reference, reduction = "umap")
DimPlot(Swarbrick_CID4290A_ER_Navin_Visvader_NORM_panCK_Reference, reduction = "umap", group.by = "orig.ident")

# Save RDS
saveRDS(Swarbrick_CID4290A_ER_Navin_Visvader_NORM_panCK_Reference, file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_inferCNV_Input_RDS/Swarbrick_CID4290A_ER_Navin_Visvader_NORM_panCK_Reference.rds")

########## Prepare Cell Annotations
# Extract & Save Meta Data
Swarbrick_CID4290A_ER_Navin_Visvader_NORM_panCK_Reference.meta.data <- as.data.frame(as.matrix(Swarbrick_CID4290A_ER_Navin_Visvader_NORM_panCK_Reference@meta.data))
# Extract pertinent columns
# Make row names first column
Swarbrick_CID4290A_ER_Navin_Visvader_NORM_panCK_Reference.meta.data <- tibble::rownames_to_column(Swarbrick_CID4290A_ER_Navin_Visvader_NORM_panCK_Reference.meta.data, "Cells")
# In Seurat, "X" becomes the cell names while "cnv.ident" is specified above
Swarbrick_CID4290A_ER_Navin_Visvader_NORM_panCK_Reference.annotations <- Swarbrick_CID4290A_ER_Navin_Visvader_NORM_panCK_Reference.meta.data[, c("Cells", "cnv.ident")]
# Delete header
names(Swarbrick_CID4290A_ER_Navin_Visvader_NORM_panCK_Reference.annotations) <- NULL
# Save as table
write.table(Swarbrick_CID4290A_ER_Navin_Visvader_NORM_panCK_Reference.annotations,"/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_CID4290A_ER_Navin_Visvader_NORM_panCK_Reference.annotations.txt",sep="\t",row.names=FALSE)
# Remove metadata and annotations
rm(Swarbrick_CID4290A_ER_Navin_Visvader_NORM_panCK_Reference.meta.data)
rm(Swarbrick_CID4290A_ER_Navin_Visvader_NORM_panCK_Reference.annotations)

########## Generate Count Matrix
# Extract Count Matrix
Swarbrick_CID4290A_ER_Navin_Visvader_NORM_panCK_Reference.count <- Swarbrick_CID4290A_ER_Navin_Visvader_NORM_panCK_Reference[["RNA"]]$counts
# Convert counts to a sparse matrix
Swarbrick_CID4290A_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as(Swarbrick_CID4290A_ER_Navin_Visvader_NORM_panCK_Reference.count, 'sparseMatrix')
rm(Swarbrick_CID4290A_ER_Navin_Visvader_NORM_panCK_Reference.count)
# Convert to requisite data format
Swarbrick_CID4290A_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as.data.frame(Swarbrick_CID4290A_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix)
# Save table
write.csv(Swarbrick_CID4290A_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix, file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_inferCNV_Input_Count_Matrix/Swarbrick_CID4290A_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv")
# Remove Objects
rm(Swarbrick_CID4290A_ER_Navin_Visvader_NORM_panCK_Reference)
rm(Swarbrick_CID4290A_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix)

########### Run inferCNV
# Note: check.names = FALSE prevents cell names from having dashes replaced with dots; essential to match annotations file
Swarbrick_CID4290A_ER_inferCNV <- CreateInfercnvObject(raw_counts_matrix=read.csv(file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_inferCNV_Input_Count_Matrix/Swarbrick_CID4290A_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv", header = TRUE, row.names = 1,check.names = FALSE),
                                             annotations_file=read.table(file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_CID4290A_ER_Navin_Visvader_NORM_panCK_Reference.annotations.txt", header = FALSE, row.names = 1),
                                             delim="\t",
                                             gene_order_file="/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/hg38_gencode_v27.txt",
                                             ref_group_names=c("Normal"))


# perform infercnv operations to reveal cnv signal
Swarbrick_CID4290A_ER_inferCNV = infercnv::run(Swarbrick_CID4290A_ER_inferCNV,
                                     cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                     out_dir="/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/Swarbrick_CID4290A_ER_panCK/",  # dir is auto-created for storing outputs
                                     cluster_by_groups=T,   # cluster
                                     denoise=T,
                                     HMM=T,
                                     num_ref_groups = 9, analysis_mode = "subclusters")

###########  Integrate data with Seurat object
# Read RDS and add inferCNV results to metadata
Swarbrick_CID4290A_ER_Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_inferCNV_Input_RDS/Swarbrick_CID4290A_ER_Navin_Visvader_NORM_panCK_Reference.rds") 
Swarbrick_CID4290A_ER_panCK_cnv <- infercnv::add_to_seurat(infercnv_output_path="/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/Swarbrick_CID4290A_ER_panCK/",
                                                 seurat_obj=Swarbrick_CID4290A_ER_Navin_Visvader_NORM_panCK_Reference, # optional
                                                 top_n=10)

# Collect Tumor cells from Seurat Object
Swarbrick_CID4290A_ER_panCK_cnv <- subset(x = Swarbrick_CID4290A_ER_panCK_cnv, subset = cnv.ident == "Tumor")
# Quantify chromosomes with alterations per cell
Swarbrick_CID4290A_ER_panCK_cnv_Pos_Chromosomes <- as.data.frame(rowSums(Swarbrick_CID4290A_ER_panCK_cnv@meta.data[,startsWith(names(Swarbrick_CID4290A_ER_panCK_cnv@meta.data),"has_cnv_")] > 0 ))
Swarbrick_CID4290A_ER_panCK_cnv_Neg_Chromosomes <- as.data.frame(rowSums(Swarbrick_CID4290A_ER_panCK_cnv@meta.data[,startsWith(names(Swarbrick_CID4290A_ER_panCK_cnv@meta.data),"has_cnv_")] == 0))
Swarbrick_CID4290A_ER_panCK_cnv_Chromosome_Quant <- as.data.frame(c(Swarbrick_CID4290A_ER_panCK_cnv_Pos_Chromosomes, Swarbrick_CID4290A_ER_panCK_cnv_Neg_Chromosomes))
colnames(Swarbrick_CID4290A_ER_panCK_cnv_Chromosome_Quant) <- c("Chromosomes CNV Pos", "Chromosomes CNV Neg")
write.csv(Swarbrick_CID4290A_ER_panCK_cnv_Chromosome_Quant, file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/Swarbrick_CID4290A_ER_panCK/Swarbrick_CID4290A_ER_panCK_cnv_Chromosome_Quant.csv")

# Save CNV Metadata
Swarbrick_CID4290A_ER_panCK_cnv.meta.data <- as.data.frame(as.matrix(Swarbrick_CID4290A_ER_panCK_cnv@meta.data))
write.csv(Swarbrick_CID4290A_ER_panCK_cnv.meta.data, file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/Swarbrick_CID4290A_ER_panCK/Swarbrick_CID4290A_ER_panCK_cnv.meta.data.csv")

# Subset Cells with inferred CNVs -------------------------------------------------------------------------------------------------------------
Swarbrick_CID4290A_ER_panCK_cnv <- subset(Swarbrick_CID4290A_ER_panCK_cnv, cells=rownames(Swarbrick_CID4290A_ER_panCK_cnv@meta.data)[rowSums(Swarbrick_CID4290A_ER_panCK_cnv@meta.data[,startsWith(names(Swarbrick_CID4290A_ER_panCK_cnv@meta.data),"has_cnv_")]) > 0 ])

# Save RDS with inferredCNVs
saveRDS(Swarbrick_CID4290A_ER_panCK_cnv, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_inferCNV/Swarbrick_CID4290A_ER_panCK_cnv.rds")


# Free memory
rm(Swarbrick_CID4290A_ER_inferCNV)
rm(Swarbrick_CID4290A_ER_Navin_Visvader_NORM_panCK_Reference)
rm(Swarbrick_CID4290A_ER_panCK_cnv)
rm(Swarbrick_CID4290A_ER_panCK_cnv_Pos_Chromosomes)
rm(Swarbrick_CID4290A_ER_panCK_cnv_Neg_Chromosomes)
rm(Swarbrick_CID4290A_ER_panCK_cnv_Chromosome_Quant)
rm(Swarbrick_CID4290A_ER_panCK_cnv.meta.data)
gc()



##############################################
#### InferCNV: Swarbrick_CID4398_ER_panCK ####
##############################################

########## Merge with normal mammary glad, extract epithelium, and process via Seurat
# Load RDS
Swarbrick_CID4398_ER_panCK <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_panCK/Swarbrick_CID4398_ER_panCK.rds")
Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Navin_Visvader_NORM_panCK_Reference.rds")

# Set identities 
colnames(x = Swarbrick_CID4398_ER_panCK[[]])
Idents(object = Swarbrick_CID4398_ER_panCK) <- 'Tumor'
Swarbrick_CID4398_ER_panCK[["cnv.ident"]] <- Idents(object = Swarbrick_CID4398_ER_panCK)
levels(x = Swarbrick_CID4398_ER_panCK)

colnames(x = Navin_Visvader_NORM_panCK_Reference[[]])
Idents(object = Navin_Visvader_NORM_panCK_Reference) <- 'Normal'
Navin_Visvader_NORM_panCK_Reference[["cnv.ident"]] <- Idents(object = Navin_Visvader_NORM_panCK_Reference)
levels(x = Navin_Visvader_NORM_panCK_Reference)

# Merge subsets into recombined RDS
Swarbrick_CID4398_ER_Navin_Visvader_NORM_panCK_Reference <- merge(x = Swarbrick_CID4398_ER_panCK, y = Navin_Visvader_NORM_panCK_Reference)
# Join Layers
Swarbrick_CID4398_ER_Navin_Visvader_NORM_panCK_Reference <- JoinLayers(Swarbrick_CID4398_ER_Navin_Visvader_NORM_panCK_Reference)

# Remove subsetted objects
rm(Swarbrick_CID4398_ER_panCK)
rm(Navin_Visvader_NORM_panCK_Reference)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples before merge

# Find Variable Features
Swarbrick_CID4398_ER_Navin_Visvader_NORM_panCK_Reference <- FindVariableFeatures(Swarbrick_CID4398_ER_Navin_Visvader_NORM_panCK_Reference, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Swarbrick_CID4398_ER_Navin_Visvader_NORM_panCK_Reference)
Swarbrick_CID4398_ER_Navin_Visvader_NORM_panCK_Reference <- ScaleData(Swarbrick_CID4398_ER_Navin_Visvader_NORM_panCK_Reference, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################

Swarbrick_CID4398_ER_Navin_Visvader_NORM_panCK_Reference <- RunPCA(Swarbrick_CID4398_ER_Navin_Visvader_NORM_panCK_Reference, verbose = FALSE)

# Examine and visualize PCA results a few different ways
print(Swarbrick_CID4398_ER_Navin_Visvader_NORM_panCK_Reference[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Swarbrick_CID4398_ER_Navin_Visvader_NORM_panCK_Reference, dims = 1:2, reduction = "pca")
DimPlot(Swarbrick_CID4398_ER_Navin_Visvader_NORM_panCK_Reference, reduction = "pca")

Swarbrick_CID4398_ER_Navin_Visvader_NORM_panCK_Reference <- FindNeighbors(Swarbrick_CID4398_ER_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
Swarbrick_CID4398_ER_Navin_Visvader_NORM_panCK_Reference <- FindClusters(Swarbrick_CID4398_ER_Navin_Visvader_NORM_panCK_Reference, resolution = 0.1)

# Umap clustering
Swarbrick_CID4398_ER_Navin_Visvader_NORM_panCK_Reference <- RunUMAP(Swarbrick_CID4398_ER_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
DimPlot(Swarbrick_CID4398_ER_Navin_Visvader_NORM_panCK_Reference, reduction = "umap")
DimPlot(Swarbrick_CID4398_ER_Navin_Visvader_NORM_panCK_Reference, reduction = "umap", group.by = "orig.ident")

# Save RDS
saveRDS(Swarbrick_CID4398_ER_Navin_Visvader_NORM_panCK_Reference, file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_inferCNV_Input_RDS/Swarbrick_CID4398_ER_Navin_Visvader_NORM_panCK_Reference.rds")

########## Prepare Cell Annotations
# Extract & Save Meta Data
Swarbrick_CID4398_ER_Navin_Visvader_NORM_panCK_Reference.meta.data <- as.data.frame(as.matrix(Swarbrick_CID4398_ER_Navin_Visvader_NORM_panCK_Reference@meta.data))
# Extract pertinent columns
# Make row names first column
Swarbrick_CID4398_ER_Navin_Visvader_NORM_panCK_Reference.meta.data <- tibble::rownames_to_column(Swarbrick_CID4398_ER_Navin_Visvader_NORM_panCK_Reference.meta.data, "Cells")
# In Seurat, "X" becomes the cell names while "cnv.ident" is specified above
Swarbrick_CID4398_ER_Navin_Visvader_NORM_panCK_Reference.annotations <- Swarbrick_CID4398_ER_Navin_Visvader_NORM_panCK_Reference.meta.data[, c("Cells", "cnv.ident")]
# Delete header
names(Swarbrick_CID4398_ER_Navin_Visvader_NORM_panCK_Reference.annotations) <- NULL
# Save as table
write.table(Swarbrick_CID4398_ER_Navin_Visvader_NORM_panCK_Reference.annotations,"/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_CID4398_ER_Navin_Visvader_NORM_panCK_Reference.annotations.txt",sep="\t",row.names=FALSE)
# Remove metadata and annotations
rm(Swarbrick_CID4398_ER_Navin_Visvader_NORM_panCK_Reference.meta.data)
rm(Swarbrick_CID4398_ER_Navin_Visvader_NORM_panCK_Reference.annotations)

########## Generate Count Matrix
# Extract Count Matrix
Swarbrick_CID4398_ER_Navin_Visvader_NORM_panCK_Reference.count <- Swarbrick_CID4398_ER_Navin_Visvader_NORM_panCK_Reference[["RNA"]]$counts# Convert counts to a sparse matrix
Swarbrick_CID4398_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as(Swarbrick_CID4398_ER_Navin_Visvader_NORM_panCK_Reference.count, 'sparseMatrix')
rm(Swarbrick_CID4398_ER_Navin_Visvader_NORM_panCK_Reference.count)
# Convert to requisite data format
Swarbrick_CID4398_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as.data.frame(Swarbrick_CID4398_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix)
# Save table
write.csv(Swarbrick_CID4398_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix, file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_inferCNV_Input_Count_Matrix/Swarbrick_CID4398_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv")
# Remove Objects
rm(Swarbrick_CID4398_ER_Navin_Visvader_NORM_panCK_Reference)
rm(Swarbrick_CID4398_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix)

########### Run inferCNV
# Note: check.names = FALSE prevents cell names from having dashes replaced with dots; essential to match annotations file
Swarbrick_CID4398_ER_inferCNV <- CreateInfercnvObject(raw_counts_matrix=read.csv(file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_inferCNV_Input_Count_Matrix/Swarbrick_CID4398_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv", header = TRUE, row.names = 1,check.names = FALSE),
                                            annotations_file=read.table(file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_CID4398_ER_Navin_Visvader_NORM_panCK_Reference.annotations.txt", header = FALSE, row.names = 1),
                                            delim="\t",
                                            gene_order_file="/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/hg38_gencode_v27.txt",
                                            ref_group_names=c("Normal"))


# perform infercnv operations to reveal cnv signal
Swarbrick_CID4398_ER_inferCNV = infercnv::run(Swarbrick_CID4398_ER_inferCNV,
                                    cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                    out_dir="/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/Swarbrick_CID4398_ER_panCK/",  # dir is auto-created for storing outputs
                                    cluster_by_groups=T,   # cluster
                                    denoise=T,
                                    HMM=T,
                                    num_ref_groups = 9, analysis_mode = "subclusters")

###########  Integrate data with Seurat object
# Read RDS and add inferCNV results to metadata
Swarbrick_CID4398_ER_Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_inferCNV_Input_RDS/Swarbrick_CID4398_ER_Navin_Visvader_NORM_panCK_Reference.rds") 
Swarbrick_CID4398_ER_panCK_cnv <- infercnv::add_to_seurat(infercnv_output_path="/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/Swarbrick_CID4398_ER_panCK/",
                                                seurat_obj=Swarbrick_CID4398_ER_Navin_Visvader_NORM_panCK_Reference, # optional
                                                top_n=10)

# Collect Tumor cells from Seurat Object
Swarbrick_CID4398_ER_panCK_cnv <- subset(x = Swarbrick_CID4398_ER_panCK_cnv, subset = cnv.ident == "Tumor")
# Quantify chromosomes with alterations per cell
Swarbrick_CID4398_ER_panCK_cnv_Pos_Chromosomes <- as.data.frame(rowSums(Swarbrick_CID4398_ER_panCK_cnv@meta.data[,startsWith(names(Swarbrick_CID4398_ER_panCK_cnv@meta.data),"has_cnv_")] > 0 ))
Swarbrick_CID4398_ER_panCK_cnv_Neg_Chromosomes <- as.data.frame(rowSums(Swarbrick_CID4398_ER_panCK_cnv@meta.data[,startsWith(names(Swarbrick_CID4398_ER_panCK_cnv@meta.data),"has_cnv_")] == 0))
Swarbrick_CID4398_ER_panCK_cnv_Chromosome_Quant <- as.data.frame(c(Swarbrick_CID4398_ER_panCK_cnv_Pos_Chromosomes, Swarbrick_CID4398_ER_panCK_cnv_Neg_Chromosomes))
colnames(Swarbrick_CID4398_ER_panCK_cnv_Chromosome_Quant) <- c("Chromosomes CNV Pos", "Chromosomes CNV Neg")
write.csv(Swarbrick_CID4398_ER_panCK_cnv_Chromosome_Quant, file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/Swarbrick_CID4398_ER_panCK/Swarbrick_CID4398_ER_panCK_cnv_Chromosome_Quant.csv")

# Save CNV Metadata
Swarbrick_CID4398_ER_panCK_cnv.meta.data <- as.data.frame(as.matrix(Swarbrick_CID4398_ER_panCK_cnv@meta.data))
write.csv(Swarbrick_CID4398_ER_panCK_cnv.meta.data, file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/Swarbrick_CID4398_ER_panCK/Swarbrick_CID4398_ER_panCK_cnv.meta.data.csv")

# Subset Cells with inferred CNVs -------------------------------------------------------------------------------------------------------------
Swarbrick_CID4398_ER_panCK_cnv <- subset(Swarbrick_CID4398_ER_panCK_cnv, cells=rownames(Swarbrick_CID4398_ER_panCK_cnv@meta.data)[rowSums(Swarbrick_CID4398_ER_panCK_cnv@meta.data[,startsWith(names(Swarbrick_CID4398_ER_panCK_cnv@meta.data),"has_cnv_")]) > 0 ])

# Save RDS with inferredCNVs
saveRDS(Swarbrick_CID4398_ER_panCK_cnv, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_inferCNV/Swarbrick_CID4398_ER_panCK_cnv.rds")


# Free memory
rm(Swarbrick_CID4398_ER_inferCNV)
rm(Swarbrick_CID4398_ER_Navin_Visvader_NORM_panCK_Reference)
rm(Swarbrick_CID4398_ER_panCK_cnv)
rm(Swarbrick_CID4398_ER_panCK_cnv_Pos_Chromosomes)
rm(Swarbrick_CID4398_ER_panCK_cnv_Neg_Chromosomes)
rm(Swarbrick_CID4398_ER_panCK_cnv_Chromosome_Quant)
rm(Swarbrick_CID4398_ER_panCK_cnv.meta.data)
gc()



##############################################
#### InferCNV: Swarbrick_CID4461_ER_panCK ####
##############################################

########## Merge with normal mammary glad, extract epithelium, and process via Seurat
# Load RDS
Swarbrick_CID4461_ER_panCK <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_panCK/Swarbrick_CID4461_ER_panCK.rds")
Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Navin_Visvader_NORM_panCK_Reference.rds")

# Set identities 
colnames(x = Swarbrick_CID4461_ER_panCK[[]])
Idents(object = Swarbrick_CID4461_ER_panCK) <- 'Tumor'
Swarbrick_CID4461_ER_panCK[["cnv.ident"]] <- Idents(object = Swarbrick_CID4461_ER_panCK)
levels(x = Swarbrick_CID4461_ER_panCK)

colnames(x = Navin_Visvader_NORM_panCK_Reference[[]])
Idents(object = Navin_Visvader_NORM_panCK_Reference) <- 'Normal'
Navin_Visvader_NORM_panCK_Reference[["cnv.ident"]] <- Idents(object = Navin_Visvader_NORM_panCK_Reference)
levels(x = Navin_Visvader_NORM_panCK_Reference)

# Merge subsets into recombined RDS
Swarbrick_CID4461_ER_Navin_Visvader_NORM_panCK_Reference <- merge(x = Swarbrick_CID4461_ER_panCK, y = Navin_Visvader_NORM_panCK_Reference)
# Join Layers
Swarbrick_CID4461_ER_Navin_Visvader_NORM_panCK_Reference <- JoinLayers(Swarbrick_CID4461_ER_Navin_Visvader_NORM_panCK_Reference)

# Remove subsetted objects
rm(Swarbrick_CID4461_ER_panCK)
rm(Navin_Visvader_NORM_panCK_Reference)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples before merge

# Find Variable Features
Swarbrick_CID4461_ER_Navin_Visvader_NORM_panCK_Reference <- FindVariableFeatures(Swarbrick_CID4461_ER_Navin_Visvader_NORM_panCK_Reference, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Swarbrick_CID4461_ER_Navin_Visvader_NORM_panCK_Reference)
Swarbrick_CID4461_ER_Navin_Visvader_NORM_panCK_Reference <- ScaleData(Swarbrick_CID4461_ER_Navin_Visvader_NORM_panCK_Reference, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################

Swarbrick_CID4461_ER_Navin_Visvader_NORM_panCK_Reference <- RunPCA(Swarbrick_CID4461_ER_Navin_Visvader_NORM_panCK_Reference, verbose = FALSE)

# Examine and visualize PCA results a few different ways
print(Swarbrick_CID4461_ER_Navin_Visvader_NORM_panCK_Reference[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Swarbrick_CID4461_ER_Navin_Visvader_NORM_panCK_Reference, dims = 1:2, reduction = "pca")
DimPlot(Swarbrick_CID4461_ER_Navin_Visvader_NORM_panCK_Reference, reduction = "pca")

Swarbrick_CID4461_ER_Navin_Visvader_NORM_panCK_Reference <- FindNeighbors(Swarbrick_CID4461_ER_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
Swarbrick_CID4461_ER_Navin_Visvader_NORM_panCK_Reference <- FindClusters(Swarbrick_CID4461_ER_Navin_Visvader_NORM_panCK_Reference, resolution = 0.1)

# Umap clustering
Swarbrick_CID4461_ER_Navin_Visvader_NORM_panCK_Reference <- RunUMAP(Swarbrick_CID4461_ER_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
DimPlot(Swarbrick_CID4461_ER_Navin_Visvader_NORM_panCK_Reference, reduction = "umap")
DimPlot(Swarbrick_CID4461_ER_Navin_Visvader_NORM_panCK_Reference, reduction = "umap", group.by = "orig.ident")

# Save RDS
saveRDS(Swarbrick_CID4461_ER_Navin_Visvader_NORM_panCK_Reference, file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_inferCNV_Input_RDS/Swarbrick_CID4461_ER_Navin_Visvader_NORM_panCK_Reference.rds")

########## Prepare Cell Annotations
# Extract & Save Meta Data
Swarbrick_CID4461_ER_Navin_Visvader_NORM_panCK_Reference.meta.data <- as.data.frame(as.matrix(Swarbrick_CID4461_ER_Navin_Visvader_NORM_panCK_Reference@meta.data))
# Extract pertinent columns
# Make row names first column
Swarbrick_CID4461_ER_Navin_Visvader_NORM_panCK_Reference.meta.data <- tibble::rownames_to_column(Swarbrick_CID4461_ER_Navin_Visvader_NORM_panCK_Reference.meta.data, "Cells")
# In Seurat, "X" becomes the cell names while "cnv.ident" is specified above
Swarbrick_CID4461_ER_Navin_Visvader_NORM_panCK_Reference.annotations <- Swarbrick_CID4461_ER_Navin_Visvader_NORM_panCK_Reference.meta.data[, c("Cells", "cnv.ident")]
# Delete header
names(Swarbrick_CID4461_ER_Navin_Visvader_NORM_panCK_Reference.annotations) <- NULL
# Save as table
write.table(Swarbrick_CID4461_ER_Navin_Visvader_NORM_panCK_Reference.annotations,"/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_CID4461_ER_Navin_Visvader_NORM_panCK_Reference.annotations.txt",sep="\t",row.names=FALSE)
# Remove metadata and annotations
rm(Swarbrick_CID4461_ER_Navin_Visvader_NORM_panCK_Reference.meta.data)
rm(Swarbrick_CID4461_ER_Navin_Visvader_NORM_panCK_Reference.annotations)

########## Generate Count Matrix
# Extract Count Matrix
Swarbrick_CID4461_ER_Navin_Visvader_NORM_panCK_Reference.count <- Swarbrick_CID4461_ER_Navin_Visvader_NORM_panCK_Reference[["RNA"]]$counts# Convert counts to a sparse matrix
Swarbrick_CID4461_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as(Swarbrick_CID4461_ER_Navin_Visvader_NORM_panCK_Reference.count, 'sparseMatrix')
rm(Swarbrick_CID4461_ER_Navin_Visvader_NORM_panCK_Reference.count)
# Convert to requisite data format
Swarbrick_CID4461_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as.data.frame(Swarbrick_CID4461_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix)
# Save table
write.csv(Swarbrick_CID4461_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix, file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_inferCNV_Input_Count_Matrix/Swarbrick_CID4461_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv")
# Remove Objects
rm(Swarbrick_CID4461_ER_Navin_Visvader_NORM_panCK_Reference)
rm(Swarbrick_CID4461_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix)

########### Run inferCNV
# Note: check.names = FALSE prevents cell names from having dashes replaced with dots; essential to match annotations file
Swarbrick_CID4461_ER_inferCNV <- CreateInfercnvObject(raw_counts_matrix=read.csv(file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_inferCNV_Input_Count_Matrix/Swarbrick_CID4461_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv", header = TRUE, row.names = 1,check.names = FALSE),
                                            annotations_file=read.table(file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_CID4461_ER_Navin_Visvader_NORM_panCK_Reference.annotations.txt", header = FALSE, row.names = 1),
                                            delim="\t",
                                            gene_order_file="/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/hg38_gencode_v27.txt",
                                            ref_group_names=c("Normal"))


# perform infercnv operations to reveal cnv signal
Swarbrick_CID4461_ER_inferCNV = infercnv::run(Swarbrick_CID4461_ER_inferCNV,
                                    cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                    out_dir="/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/Swarbrick_CID4461_ER_panCK/",  # dir is auto-created for storing outputs
                                    cluster_by_groups=T,   # cluster
                                    denoise=T,
                                    HMM=T,
                                    num_ref_groups = 9, analysis_mode = "subclusters")

###########  Integrate data with Seurat object
# Read RDS and add inferCNV results to metadata
Swarbrick_CID4461_ER_Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_inferCNV_Input_RDS/Swarbrick_CID4461_ER_Navin_Visvader_NORM_panCK_Reference.rds") 
Swarbrick_CID4461_ER_panCK_cnv <- infercnv::add_to_seurat(infercnv_output_path="/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/Swarbrick_CID4461_ER_panCK/",
                                                seurat_obj=Swarbrick_CID4461_ER_Navin_Visvader_NORM_panCK_Reference, # optional
                                                top_n=10)

# Collect Tumor cells from Seurat Object
Swarbrick_CID4461_ER_panCK_cnv <- subset(x = Swarbrick_CID4461_ER_panCK_cnv, subset = cnv.ident == "Tumor")
# Quantify chromosomes with alterations per cell
Swarbrick_CID4461_ER_panCK_cnv_Pos_Chromosomes <- as.data.frame(rowSums(Swarbrick_CID4461_ER_panCK_cnv@meta.data[,startsWith(names(Swarbrick_CID4461_ER_panCK_cnv@meta.data),"has_cnv_")] > 0 ))
Swarbrick_CID4461_ER_panCK_cnv_Neg_Chromosomes <- as.data.frame(rowSums(Swarbrick_CID4461_ER_panCK_cnv@meta.data[,startsWith(names(Swarbrick_CID4461_ER_panCK_cnv@meta.data),"has_cnv_")] == 0))
Swarbrick_CID4461_ER_panCK_cnv_Chromosome_Quant <- as.data.frame(c(Swarbrick_CID4461_ER_panCK_cnv_Pos_Chromosomes, Swarbrick_CID4461_ER_panCK_cnv_Neg_Chromosomes))
colnames(Swarbrick_CID4461_ER_panCK_cnv_Chromosome_Quant) <- c("Chromosomes CNV Pos", "Chromosomes CNV Neg")
write.csv(Swarbrick_CID4461_ER_panCK_cnv_Chromosome_Quant, file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/Swarbrick_CID4461_ER_panCK/Swarbrick_CID4461_ER_panCK_cnv_Chromosome_Quant.csv")

# Save CNV Metadata
Swarbrick_CID4461_ER_panCK_cnv.meta.data <- as.data.frame(as.matrix(Swarbrick_CID4461_ER_panCK_cnv@meta.data))
write.csv(Swarbrick_CID4461_ER_panCK_cnv.meta.data, file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/Swarbrick_CID4461_ER_panCK/Swarbrick_CID4461_ER_panCK_cnv.meta.data.csv")

# Subset Cells with inferred CNVs -------------------------------------------------------------------------------------------------------------
Swarbrick_CID4461_ER_panCK_cnv <- subset(Swarbrick_CID4461_ER_panCK_cnv, cells=rownames(Swarbrick_CID4461_ER_panCK_cnv@meta.data)[rowSums(Swarbrick_CID4461_ER_panCK_cnv@meta.data[,startsWith(names(Swarbrick_CID4461_ER_panCK_cnv@meta.data),"has_cnv_")]) > 0 ])

# Save RDS with inferredCNVs
saveRDS(Swarbrick_CID4461_ER_panCK_cnv, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_inferCNV/Swarbrick_CID4461_ER_panCK_cnv.rds")


# Free memory
rm(Swarbrick_CID4461_ER_inferCNV)
rm(Swarbrick_CID4461_ER_Navin_Visvader_NORM_panCK_Reference)
rm(Swarbrick_CID4461_ER_panCK_cnv)
rm(Swarbrick_CID4461_ER_panCK_cnv_Pos_Chromosomes)
rm(Swarbrick_CID4461_ER_panCK_cnv_Neg_Chromosomes)
rm(Swarbrick_CID4461_ER_panCK_cnv_Chromosome_Quant)
rm(Swarbrick_CID4461_ER_panCK_cnv.meta.data)
gc()



##############################################
#### InferCNV: Swarbrick_CID4463_ER_panCK ####
##############################################

########## Merge with normal mammary glad, extract epithelium, and process via Seurat
# Load RDS
Swarbrick_CID4463_ER_panCK <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_panCK/Swarbrick_CID4463_ER_panCK.rds")
Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Navin_Visvader_NORM_panCK_Reference.rds")

# Set identities 
colnames(x = Swarbrick_CID4463_ER_panCK[[]])
Idents(object = Swarbrick_CID4463_ER_panCK) <- 'Tumor'
Swarbrick_CID4463_ER_panCK[["cnv.ident"]] <- Idents(object = Swarbrick_CID4463_ER_panCK)
levels(x = Swarbrick_CID4463_ER_panCK)

colnames(x = Navin_Visvader_NORM_panCK_Reference[[]])
Idents(object = Navin_Visvader_NORM_panCK_Reference) <- 'Normal'
Navin_Visvader_NORM_panCK_Reference[["cnv.ident"]] <- Idents(object = Navin_Visvader_NORM_panCK_Reference)
levels(x = Navin_Visvader_NORM_panCK_Reference)

# Merge subsets into recombined RDS
Swarbrick_CID4463_ER_Navin_Visvader_NORM_panCK_Reference <- merge(x = Swarbrick_CID4463_ER_panCK, y = Navin_Visvader_NORM_panCK_Reference)
# Join Layers
Swarbrick_CID4463_ER_Navin_Visvader_NORM_panCK_Reference <- JoinLayers(Swarbrick_CID4463_ER_Navin_Visvader_NORM_panCK_Reference)

# Remove subsetted objects
rm(Swarbrick_CID4463_ER_panCK)
rm(Navin_Visvader_NORM_panCK_Reference)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples before merge

# Find Variable Features
Swarbrick_CID4463_ER_Navin_Visvader_NORM_panCK_Reference <- FindVariableFeatures(Swarbrick_CID4463_ER_Navin_Visvader_NORM_panCK_Reference, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Swarbrick_CID4463_ER_Navin_Visvader_NORM_panCK_Reference)
Swarbrick_CID4463_ER_Navin_Visvader_NORM_panCK_Reference <- ScaleData(Swarbrick_CID4463_ER_Navin_Visvader_NORM_panCK_Reference, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################

Swarbrick_CID4463_ER_Navin_Visvader_NORM_panCK_Reference <- RunPCA(Swarbrick_CID4463_ER_Navin_Visvader_NORM_panCK_Reference, verbose = FALSE)

# Examine and visualize PCA results a few different ways
print(Swarbrick_CID4463_ER_Navin_Visvader_NORM_panCK_Reference[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Swarbrick_CID4463_ER_Navin_Visvader_NORM_panCK_Reference, dims = 1:2, reduction = "pca")
DimPlot(Swarbrick_CID4463_ER_Navin_Visvader_NORM_panCK_Reference, reduction = "pca")

Swarbrick_CID4463_ER_Navin_Visvader_NORM_panCK_Reference <- FindNeighbors(Swarbrick_CID4463_ER_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
Swarbrick_CID4463_ER_Navin_Visvader_NORM_panCK_Reference <- FindClusters(Swarbrick_CID4463_ER_Navin_Visvader_NORM_panCK_Reference, resolution = 0.1)

# Umap clustering
Swarbrick_CID4463_ER_Navin_Visvader_NORM_panCK_Reference <- RunUMAP(Swarbrick_CID4463_ER_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
DimPlot(Swarbrick_CID4463_ER_Navin_Visvader_NORM_panCK_Reference, reduction = "umap")
DimPlot(Swarbrick_CID4463_ER_Navin_Visvader_NORM_panCK_Reference, reduction = "umap", group.by = "orig.ident")

# Save RDS
saveRDS(Swarbrick_CID4463_ER_Navin_Visvader_NORM_panCK_Reference, file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_inferCNV_Input_RDS/Swarbrick_CID4463_ER_Navin_Visvader_NORM_panCK_Reference.rds")

########## Prepare Cell Annotations
# Extract & Save Meta Data
Swarbrick_CID4463_ER_Navin_Visvader_NORM_panCK_Reference.meta.data <- as.data.frame(as.matrix(Swarbrick_CID4463_ER_Navin_Visvader_NORM_panCK_Reference@meta.data))
# Extract pertinent columns
# Make row names first column
Swarbrick_CID4463_ER_Navin_Visvader_NORM_panCK_Reference.meta.data <- tibble::rownames_to_column(Swarbrick_CID4463_ER_Navin_Visvader_NORM_panCK_Reference.meta.data, "Cells")
# In Seurat, "X" becomes the cell names while "cnv.ident" is specified above
Swarbrick_CID4463_ER_Navin_Visvader_NORM_panCK_Reference.annotations <- Swarbrick_CID4463_ER_Navin_Visvader_NORM_panCK_Reference.meta.data[, c("Cells", "cnv.ident")]
# Delete header
names(Swarbrick_CID4463_ER_Navin_Visvader_NORM_panCK_Reference.annotations) <- NULL
# Save as table
write.table(Swarbrick_CID4463_ER_Navin_Visvader_NORM_panCK_Reference.annotations,"/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_CID4463_ER_Navin_Visvader_NORM_panCK_Reference.annotations.txt",sep="\t",row.names=FALSE)
# Remove metadata and annotations
rm(Swarbrick_CID4463_ER_Navin_Visvader_NORM_panCK_Reference.meta.data)
rm(Swarbrick_CID4463_ER_Navin_Visvader_NORM_panCK_Reference.annotations)

########## Generate Count Matrix
# Extract Count Matrix
Swarbrick_CID4463_ER_Navin_Visvader_NORM_panCK_Reference.count <- Swarbrick_CID4463_ER_Navin_Visvader_NORM_panCK_Reference[["RNA"]]$counts# Convert counts to a sparse matrix
Swarbrick_CID4463_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as(Swarbrick_CID4463_ER_Navin_Visvader_NORM_panCK_Reference.count, 'sparseMatrix')
rm(Swarbrick_CID4463_ER_Navin_Visvader_NORM_panCK_Reference.count)
# Convert to requisite data format
Swarbrick_CID4463_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as.data.frame(Swarbrick_CID4463_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix)
# Save table
write.csv(Swarbrick_CID4463_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix, file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_inferCNV_Input_Count_Matrix/Swarbrick_CID4463_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv")
# Remove Objects
rm(Swarbrick_CID4463_ER_Navin_Visvader_NORM_panCK_Reference)
rm(Swarbrick_CID4463_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix)

########### Run inferCNV
# Note: check.names = FALSE prevents cell names from having dashes replaced with dots; essential to match annotations file
Swarbrick_CID4463_ER_inferCNV <- CreateInfercnvObject(raw_counts_matrix=read.csv(file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_inferCNV_Input_Count_Matrix/Swarbrick_CID4463_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv", header = TRUE, row.names = 1,check.names = FALSE),
                                            annotations_file=read.table(file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_CID4463_ER_Navin_Visvader_NORM_panCK_Reference.annotations.txt", header = FALSE, row.names = 1),
                                            delim="\t",
                                            gene_order_file="/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/hg38_gencode_v27.txt",
                                            ref_group_names=c("Normal"))


# perform infercnv operations to reveal cnv signal
Swarbrick_CID4463_ER_inferCNV = infercnv::run(Swarbrick_CID4463_ER_inferCNV,
                                    cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                    out_dir="/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/Swarbrick_CID4463_ER_panCK/",  # dir is auto-created for storing outputs
                                    cluster_by_groups=T,   # cluster
                                    denoise=T,
                                    HMM=T,
                                    num_ref_groups = 9, analysis_mode = "subclusters")

###########  Integrate data with Seurat object
# Read RDS and add inferCNV results to metadata
Swarbrick_CID4463_ER_Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_inferCNV_Input_RDS/Swarbrick_CID4463_ER_Navin_Visvader_NORM_panCK_Reference.rds") 
Swarbrick_CID4463_ER_panCK_cnv <- infercnv::add_to_seurat(infercnv_output_path="/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/Swarbrick_CID4463_ER_panCK/",
                                                seurat_obj=Swarbrick_CID4463_ER_Navin_Visvader_NORM_panCK_Reference, # optional
                                                top_n=10)

# Collect Tumor cells from Seurat Object
Swarbrick_CID4463_ER_panCK_cnv <- subset(x = Swarbrick_CID4463_ER_panCK_cnv, subset = cnv.ident == "Tumor")
# Quantify chromosomes with alterations per cell
Swarbrick_CID4463_ER_panCK_cnv_Pos_Chromosomes <- as.data.frame(rowSums(Swarbrick_CID4463_ER_panCK_cnv@meta.data[,startsWith(names(Swarbrick_CID4463_ER_panCK_cnv@meta.data),"has_cnv_")] > 0 ))
Swarbrick_CID4463_ER_panCK_cnv_Neg_Chromosomes <- as.data.frame(rowSums(Swarbrick_CID4463_ER_panCK_cnv@meta.data[,startsWith(names(Swarbrick_CID4463_ER_panCK_cnv@meta.data),"has_cnv_")] == 0))
Swarbrick_CID4463_ER_panCK_cnv_Chromosome_Quant <- as.data.frame(c(Swarbrick_CID4463_ER_panCK_cnv_Pos_Chromosomes, Swarbrick_CID4463_ER_panCK_cnv_Neg_Chromosomes))
colnames(Swarbrick_CID4463_ER_panCK_cnv_Chromosome_Quant) <- c("Chromosomes CNV Pos", "Chromosomes CNV Neg")
write.csv(Swarbrick_CID4463_ER_panCK_cnv_Chromosome_Quant, file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/Swarbrick_CID4463_ER_panCK/Swarbrick_CID4463_ER_panCK_cnv_Chromosome_Quant.csv")

# Save CNV Metadata
Swarbrick_CID4463_ER_panCK_cnv.meta.data <- as.data.frame(as.matrix(Swarbrick_CID4463_ER_panCK_cnv@meta.data))
write.csv(Swarbrick_CID4463_ER_panCK_cnv.meta.data, file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/Swarbrick_CID4463_ER_panCK/Swarbrick_CID4463_ER_panCK_cnv.meta.data.csv")

# Subset Cells with inferred CNVs -------------------------------------------------------------------------------------------------------------
Swarbrick_CID4463_ER_panCK_cnv <- subset(Swarbrick_CID4463_ER_panCK_cnv, cells=rownames(Swarbrick_CID4463_ER_panCK_cnv@meta.data)[rowSums(Swarbrick_CID4463_ER_panCK_cnv@meta.data[,startsWith(names(Swarbrick_CID4463_ER_panCK_cnv@meta.data),"has_cnv_")]) > 0 ])

# Save RDS with inferredCNVs
saveRDS(Swarbrick_CID4463_ER_panCK_cnv, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_inferCNV/Swarbrick_CID4463_ER_panCK_cnv.rds")


# Free memory
rm(Swarbrick_CID4463_ER_inferCNV)
rm(Swarbrick_CID4463_ER_Navin_Visvader_NORM_panCK_Reference)
rm(Swarbrick_CID4463_ER_panCK_cnv)
rm(Swarbrick_CID4463_ER_panCK_cnv_Pos_Chromosomes)
rm(Swarbrick_CID4463_ER_panCK_cnv_Neg_Chromosomes)
rm(Swarbrick_CID4463_ER_panCK_cnv_Chromosome_Quant)
rm(Swarbrick_CID4463_ER_panCK_cnv.meta.data)
gc()



##############################################
#### InferCNV: Swarbrick_CID4471_ER_panCK ####
##############################################

########## Merge with normal mammary glad, extract epithelium, and process via Seurat
# Load RDS
Swarbrick_CID4471_ER_panCK <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_panCK/Swarbrick_CID4471_ER_panCK.rds")
Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Navin_Visvader_NORM_panCK_Reference.rds")

# Set identities 
colnames(x = Swarbrick_CID4471_ER_panCK[[]])
Idents(object = Swarbrick_CID4471_ER_panCK) <- 'Tumor'
Swarbrick_CID4471_ER_panCK[["cnv.ident"]] <- Idents(object = Swarbrick_CID4471_ER_panCK)
levels(x = Swarbrick_CID4471_ER_panCK)

colnames(x = Navin_Visvader_NORM_panCK_Reference[[]])
Idents(object = Navin_Visvader_NORM_panCK_Reference) <- 'Normal'
Navin_Visvader_NORM_panCK_Reference[["cnv.ident"]] <- Idents(object = Navin_Visvader_NORM_panCK_Reference)
levels(x = Navin_Visvader_NORM_panCK_Reference)

# Merge subsets into recombined RDS
Swarbrick_CID4471_ER_Navin_Visvader_NORM_panCK_Reference <- merge(x = Swarbrick_CID4471_ER_panCK, y = Navin_Visvader_NORM_panCK_Reference)
# Join Layers
Swarbrick_CID4471_ER_Navin_Visvader_NORM_panCK_Reference <- JoinLayers(Swarbrick_CID4471_ER_Navin_Visvader_NORM_panCK_Reference)

# Remove subsetted objects
rm(Swarbrick_CID4471_ER_panCK)
rm(Navin_Visvader_NORM_panCK_Reference)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples before merge

# Find Variable Features
Swarbrick_CID4471_ER_Navin_Visvader_NORM_panCK_Reference <- FindVariableFeatures(Swarbrick_CID4471_ER_Navin_Visvader_NORM_panCK_Reference, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Swarbrick_CID4471_ER_Navin_Visvader_NORM_panCK_Reference)
Swarbrick_CID4471_ER_Navin_Visvader_NORM_panCK_Reference <- ScaleData(Swarbrick_CID4471_ER_Navin_Visvader_NORM_panCK_Reference, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################

Swarbrick_CID4471_ER_Navin_Visvader_NORM_panCK_Reference <- RunPCA(Swarbrick_CID4471_ER_Navin_Visvader_NORM_panCK_Reference, verbose = FALSE)

# Examine and visualize PCA results a few different ways
print(Swarbrick_CID4471_ER_Navin_Visvader_NORM_panCK_Reference[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Swarbrick_CID4471_ER_Navin_Visvader_NORM_panCK_Reference, dims = 1:2, reduction = "pca")
DimPlot(Swarbrick_CID4471_ER_Navin_Visvader_NORM_panCK_Reference, reduction = "pca")

Swarbrick_CID4471_ER_Navin_Visvader_NORM_panCK_Reference <- FindNeighbors(Swarbrick_CID4471_ER_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
Swarbrick_CID4471_ER_Navin_Visvader_NORM_panCK_Reference <- FindClusters(Swarbrick_CID4471_ER_Navin_Visvader_NORM_panCK_Reference, resolution = 0.1)

# Umap clustering
Swarbrick_CID4471_ER_Navin_Visvader_NORM_panCK_Reference <- RunUMAP(Swarbrick_CID4471_ER_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
DimPlot(Swarbrick_CID4471_ER_Navin_Visvader_NORM_panCK_Reference, reduction = "umap")
DimPlot(Swarbrick_CID4471_ER_Navin_Visvader_NORM_panCK_Reference, reduction = "umap", group.by = "orig.ident")

# Save RDS
saveRDS(Swarbrick_CID4471_ER_Navin_Visvader_NORM_panCK_Reference, file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_inferCNV_Input_RDS/Swarbrick_CID4471_ER_Navin_Visvader_NORM_panCK_Reference.rds")

########## Prepare Cell Annotations
# Extract & Save Meta Data
Swarbrick_CID4471_ER_Navin_Visvader_NORM_panCK_Reference.meta.data <- as.data.frame(as.matrix(Swarbrick_CID4471_ER_Navin_Visvader_NORM_panCK_Reference@meta.data))
# Extract pertinent columns
# Make row names first column
Swarbrick_CID4471_ER_Navin_Visvader_NORM_panCK_Reference.meta.data <- tibble::rownames_to_column(Swarbrick_CID4471_ER_Navin_Visvader_NORM_panCK_Reference.meta.data, "Cells")
# In Seurat, "X" becomes the cell names while "cnv.ident" is specified above
Swarbrick_CID4471_ER_Navin_Visvader_NORM_panCK_Reference.annotations <- Swarbrick_CID4471_ER_Navin_Visvader_NORM_panCK_Reference.meta.data[, c("Cells", "cnv.ident")]
# Delete header
names(Swarbrick_CID4471_ER_Navin_Visvader_NORM_panCK_Reference.annotations) <- NULL
# Save as table
write.table(Swarbrick_CID4471_ER_Navin_Visvader_NORM_panCK_Reference.annotations,"/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_CID4471_ER_Navin_Visvader_NORM_panCK_Reference.annotations.txt",sep="\t",row.names=FALSE)
# Remove metadata and annotations
rm(Swarbrick_CID4471_ER_Navin_Visvader_NORM_panCK_Reference.meta.data)
rm(Swarbrick_CID4471_ER_Navin_Visvader_NORM_panCK_Reference.annotations)

########## Generate Count Matrix
# Extract Count Matrix
Swarbrick_CID4471_ER_Navin_Visvader_NORM_panCK_Reference.count <- Swarbrick_CID4471_ER_Navin_Visvader_NORM_panCK_Reference[["RNA"]]$counts# Convert counts to a sparse matrix
Swarbrick_CID4471_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as(Swarbrick_CID4471_ER_Navin_Visvader_NORM_panCK_Reference.count, 'sparseMatrix')
rm(Swarbrick_CID4471_ER_Navin_Visvader_NORM_panCK_Reference.count)
# Convert to requisite data format
Swarbrick_CID4471_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as.data.frame(Swarbrick_CID4471_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix)
# Save table
write.csv(Swarbrick_CID4471_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix, file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_inferCNV_Input_Count_Matrix/Swarbrick_CID4471_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv")
# Remove Objects
rm(Swarbrick_CID4471_ER_Navin_Visvader_NORM_panCK_Reference)
rm(Swarbrick_CID4471_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix)

########### Run inferCNV
# Note: check.names = FALSE prevents cell names from having dashes replaced with dots; essential to match annotations file
Swarbrick_CID4471_ER_inferCNV <- CreateInfercnvObject(raw_counts_matrix=read.csv(file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_inferCNV_Input_Count_Matrix/Swarbrick_CID4471_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv", header = TRUE, row.names = 1,check.names = FALSE),
                                            annotations_file=read.table(file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_CID4471_ER_Navin_Visvader_NORM_panCK_Reference.annotations.txt", header = FALSE, row.names = 1),
                                            delim="\t",
                                            gene_order_file="/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/hg38_gencode_v27.txt",
                                            ref_group_names=c("Normal"))


# perform infercnv operations to reveal cnv signal
Swarbrick_CID4471_ER_inferCNV = infercnv::run(Swarbrick_CID4471_ER_inferCNV,
                                    cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                    out_dir="/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/Swarbrick_CID4471_ER_panCK/",  # dir is auto-created for storing outputs
                                    cluster_by_groups=T,   # cluster
                                    denoise=T,
                                    HMM=T,
                                    num_ref_groups = 9, analysis_mode = "subclusters")

###########  Integrate data with Seurat object
# Read RDS and add inferCNV results to metadata
Swarbrick_CID4471_ER_Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_inferCNV_Input_RDS/Swarbrick_CID4471_ER_Navin_Visvader_NORM_panCK_Reference.rds") 
Swarbrick_CID4471_ER_panCK_cnv <- infercnv::add_to_seurat(infercnv_output_path="/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/Swarbrick_CID4471_ER_panCK/",
                                                seurat_obj=Swarbrick_CID4471_ER_Navin_Visvader_NORM_panCK_Reference, # optional
                                                top_n=10)

# Collect Tumor cells from Seurat Object
Swarbrick_CID4471_ER_panCK_cnv <- subset(x = Swarbrick_CID4471_ER_panCK_cnv, subset = cnv.ident == "Tumor")
# Quantify chromosomes with alterations per cell
Swarbrick_CID4471_ER_panCK_cnv_Pos_Chromosomes <- as.data.frame(rowSums(Swarbrick_CID4471_ER_panCK_cnv@meta.data[,startsWith(names(Swarbrick_CID4471_ER_panCK_cnv@meta.data),"has_cnv_")] > 0 ))
Swarbrick_CID4471_ER_panCK_cnv_Neg_Chromosomes <- as.data.frame(rowSums(Swarbrick_CID4471_ER_panCK_cnv@meta.data[,startsWith(names(Swarbrick_CID4471_ER_panCK_cnv@meta.data),"has_cnv_")] == 0))
Swarbrick_CID4471_ER_panCK_cnv_Chromosome_Quant <- as.data.frame(c(Swarbrick_CID4471_ER_panCK_cnv_Pos_Chromosomes, Swarbrick_CID4471_ER_panCK_cnv_Neg_Chromosomes))
colnames(Swarbrick_CID4471_ER_panCK_cnv_Chromosome_Quant) <- c("Chromosomes CNV Pos", "Chromosomes CNV Neg")
write.csv(Swarbrick_CID4471_ER_panCK_cnv_Chromosome_Quant, file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/Swarbrick_CID4471_ER_panCK/Swarbrick_CID4471_ER_panCK_cnv_Chromosome_Quant.csv")

# Save CNV Metadata
Swarbrick_CID4471_ER_panCK_cnv.meta.data <- as.data.frame(as.matrix(Swarbrick_CID4471_ER_panCK_cnv@meta.data))
write.csv(Swarbrick_CID4471_ER_panCK_cnv.meta.data, file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/Swarbrick_CID4471_ER_panCK/Swarbrick_CID4471_ER_panCK_cnv.meta.data.csv")

# Subset Cells with inferred CNVs -------------------------------------------------------------------------------------------------------------
Swarbrick_CID4471_ER_panCK_cnv <- subset(Swarbrick_CID4471_ER_panCK_cnv, cells=rownames(Swarbrick_CID4471_ER_panCK_cnv@meta.data)[rowSums(Swarbrick_CID4471_ER_panCK_cnv@meta.data[,startsWith(names(Swarbrick_CID4471_ER_panCK_cnv@meta.data),"has_cnv_")]) > 0 ])

# Save RDS with inferredCNVs
saveRDS(Swarbrick_CID4471_ER_panCK_cnv, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_inferCNV/Swarbrick_CID4471_ER_panCK_cnv.rds")


# Free memory
rm(Swarbrick_CID4471_ER_inferCNV)
rm(Swarbrick_CID4471_ER_Navin_Visvader_NORM_panCK_Reference)
rm(Swarbrick_CID4471_ER_panCK_cnv)
rm(Swarbrick_CID4471_ER_panCK_cnv_Pos_Chromosomes)
rm(Swarbrick_CID4471_ER_panCK_cnv_Neg_Chromosomes)
rm(Swarbrick_CID4471_ER_panCK_cnv_Chromosome_Quant)
rm(Swarbrick_CID4471_ER_panCK_cnv.meta.data)
gc()



###############################################
#### InferCNV: Swarbrick_CID4530N_ER_panCK ####
###############################################

########## Merge with normal mammary glad, extract epithelium, and process via Seurat
# Load RDS
Swarbrick_CID4530N_ER_panCK <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_panCK/Swarbrick_CID4530N_ER_panCK.rds")
Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Navin_Visvader_NORM_panCK_Reference.rds")

# Set identities 
colnames(x = Swarbrick_CID4530N_ER_panCK[[]])
Idents(object = Swarbrick_CID4530N_ER_panCK) <- 'Tumor'
Swarbrick_CID4530N_ER_panCK[["cnv.ident"]] <- Idents(object = Swarbrick_CID4530N_ER_panCK)
levels(x = Swarbrick_CID4530N_ER_panCK)

colnames(x = Navin_Visvader_NORM_panCK_Reference[[]])
Idents(object = Navin_Visvader_NORM_panCK_Reference) <- 'Normal'
Navin_Visvader_NORM_panCK_Reference[["cnv.ident"]] <- Idents(object = Navin_Visvader_NORM_panCK_Reference)
levels(x = Navin_Visvader_NORM_panCK_Reference)

# Merge subsets into recombined RDS
Swarbrick_CID4530N_ER_Navin_Visvader_NORM_panCK_Reference <- merge(x = Swarbrick_CID4530N_ER_panCK, y = Navin_Visvader_NORM_panCK_Reference)
# Join Layers
Swarbrick_CID4530N_ER_Navin_Visvader_NORM_panCK_Reference <- JoinLayers(Swarbrick_CID4530N_ER_Navin_Visvader_NORM_panCK_Reference)

# Remove subsetted objects
rm(Swarbrick_CID4530N_ER_panCK)
rm(Navin_Visvader_NORM_panCK_Reference)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples before merge

# Find Variable Features
Swarbrick_CID4530N_ER_Navin_Visvader_NORM_panCK_Reference <- FindVariableFeatures(Swarbrick_CID4530N_ER_Navin_Visvader_NORM_panCK_Reference, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Swarbrick_CID4530N_ER_Navin_Visvader_NORM_panCK_Reference)
Swarbrick_CID4530N_ER_Navin_Visvader_NORM_panCK_Reference <- ScaleData(Swarbrick_CID4530N_ER_Navin_Visvader_NORM_panCK_Reference, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################

Swarbrick_CID4530N_ER_Navin_Visvader_NORM_panCK_Reference <- RunPCA(Swarbrick_CID4530N_ER_Navin_Visvader_NORM_panCK_Reference, verbose = FALSE)

# Examine and visualize PCA results a few different ways
print(Swarbrick_CID4530N_ER_Navin_Visvader_NORM_panCK_Reference[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Swarbrick_CID4530N_ER_Navin_Visvader_NORM_panCK_Reference, dims = 1:2, reduction = "pca")
DimPlot(Swarbrick_CID4530N_ER_Navin_Visvader_NORM_panCK_Reference, reduction = "pca")

Swarbrick_CID4530N_ER_Navin_Visvader_NORM_panCK_Reference <- FindNeighbors(Swarbrick_CID4530N_ER_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
Swarbrick_CID4530N_ER_Navin_Visvader_NORM_panCK_Reference <- FindClusters(Swarbrick_CID4530N_ER_Navin_Visvader_NORM_panCK_Reference, resolution = 0.1)

# Umap clustering
Swarbrick_CID4530N_ER_Navin_Visvader_NORM_panCK_Reference <- RunUMAP(Swarbrick_CID4530N_ER_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
DimPlot(Swarbrick_CID4530N_ER_Navin_Visvader_NORM_panCK_Reference, reduction = "umap")
DimPlot(Swarbrick_CID4530N_ER_Navin_Visvader_NORM_panCK_Reference, reduction = "umap", group.by = "orig.ident")

# Save RDS
saveRDS(Swarbrick_CID4530N_ER_Navin_Visvader_NORM_panCK_Reference, file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_inferCNV_Input_RDS/Swarbrick_CID4530N_ER_Navin_Visvader_NORM_panCK_Reference.rds")

########## Prepare Cell Annotations
# Extract & Save Meta Data
Swarbrick_CID4530N_ER_Navin_Visvader_NORM_panCK_Reference.meta.data <- as.data.frame(as.matrix(Swarbrick_CID4530N_ER_Navin_Visvader_NORM_panCK_Reference@meta.data))
# Extract pertinent columns
# Make row names first column
Swarbrick_CID4530N_ER_Navin_Visvader_NORM_panCK_Reference.meta.data <- tibble::rownames_to_column(Swarbrick_CID4530N_ER_Navin_Visvader_NORM_panCK_Reference.meta.data, "Cells")
# In Seurat, "X" becomes the cell names while "cnv.ident" is specified above
Swarbrick_CID4530N_ER_Navin_Visvader_NORM_panCK_Reference.annotations <- Swarbrick_CID4530N_ER_Navin_Visvader_NORM_panCK_Reference.meta.data[, c("Cells", "cnv.ident")]
# Delete header
names(Swarbrick_CID4530N_ER_Navin_Visvader_NORM_panCK_Reference.annotations) <- NULL
# Save as table
write.table(Swarbrick_CID4530N_ER_Navin_Visvader_NORM_panCK_Reference.annotations,"/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_CID4530N_ER_Navin_Visvader_NORM_panCK_Reference.annotations.txt",sep="\t",row.names=FALSE)
# Remove metadata and annotations
rm(Swarbrick_CID4530N_ER_Navin_Visvader_NORM_panCK_Reference.meta.data)
rm(Swarbrick_CID4530N_ER_Navin_Visvader_NORM_panCK_Reference.annotations)

########## Generate Count Matrix
# Extract Count Matrix
Swarbrick_CID4530N_ER_Navin_Visvader_NORM_panCK_Reference.count <- Swarbrick_CID4530N_ER_Navin_Visvader_NORM_panCK_Reference[["RNA"]]$counts# Convert counts to a sparse matrix
Swarbrick_CID4530N_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as(Swarbrick_CID4530N_ER_Navin_Visvader_NORM_panCK_Reference.count, 'sparseMatrix')
rm(Swarbrick_CID4530N_ER_Navin_Visvader_NORM_panCK_Reference.count)
# Convert to requisite data format
Swarbrick_CID4530N_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as.data.frame(Swarbrick_CID4530N_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix)
# Save table
write.csv(Swarbrick_CID4530N_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix, file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_inferCNV_Input_Count_Matrix/Swarbrick_CID4530N_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv")
# Remove Objects
rm(Swarbrick_CID4530N_ER_Navin_Visvader_NORM_panCK_Reference)
rm(Swarbrick_CID4530N_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix)

########### Run inferCNV
# Note: check.names = FALSE prevents cell names from having dashes replaced with dots; essential to match annotations file
Swarbrick_CID4530N_ER_inferCNV <- CreateInfercnvObject(raw_counts_matrix=read.csv(file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_inferCNV_Input_Count_Matrix/Swarbrick_CID4530N_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv", header = TRUE, row.names = 1,check.names = FALSE),
                                             annotations_file=read.table(file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_CID4530N_ER_Navin_Visvader_NORM_panCK_Reference.annotations.txt", header = FALSE, row.names = 1),
                                             delim="\t",
                                             gene_order_file="/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/hg38_gencode_v27.txt",
                                             ref_group_names=c("Normal"))


# perform infercnv operations to reveal cnv signal
Swarbrick_CID4530N_ER_inferCNV = infercnv::run(Swarbrick_CID4530N_ER_inferCNV,
                                     cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                     out_dir="/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/Swarbrick_CID4530N_ER_panCK/",  # dir is auto-created for storing outputs
                                     cluster_by_groups=T,   # cluster
                                     denoise=T,
                                     HMM=T,
                                     num_ref_groups = 9, analysis_mode = "subclusters")

###########  Integrate data with Seurat object
# Read RDS and add inferCNV results to metadata
Swarbrick_CID4530N_ER_Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_inferCNV_Input_RDS/Swarbrick_CID4530N_ER_Navin_Visvader_NORM_panCK_Reference.rds") 
Swarbrick_CID4530N_ER_panCK_cnv <- infercnv::add_to_seurat(infercnv_output_path="/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/Swarbrick_CID4530N_ER_panCK/",
                                                 seurat_obj=Swarbrick_CID4530N_ER_Navin_Visvader_NORM_panCK_Reference, # optional
                                                 top_n=10)

# Collect Tumor cells from Seurat Object
Swarbrick_CID4530N_ER_panCK_cnv <- subset(x = Swarbrick_CID4530N_ER_panCK_cnv, subset = cnv.ident == "Tumor")
# Quantify chromosomes with alterations per cell
Swarbrick_CID4530N_ER_panCK_cnv_Pos_Chromosomes <- as.data.frame(rowSums(Swarbrick_CID4530N_ER_panCK_cnv@meta.data[,startsWith(names(Swarbrick_CID4530N_ER_panCK_cnv@meta.data),"has_cnv_")] > 0 ))
Swarbrick_CID4530N_ER_panCK_cnv_Neg_Chromosomes <- as.data.frame(rowSums(Swarbrick_CID4530N_ER_panCK_cnv@meta.data[,startsWith(names(Swarbrick_CID4530N_ER_panCK_cnv@meta.data),"has_cnv_")] == 0))
Swarbrick_CID4530N_ER_panCK_cnv_Chromosome_Quant <- as.data.frame(c(Swarbrick_CID4530N_ER_panCK_cnv_Pos_Chromosomes, Swarbrick_CID4530N_ER_panCK_cnv_Neg_Chromosomes))
colnames(Swarbrick_CID4530N_ER_panCK_cnv_Chromosome_Quant) <- c("Chromosomes CNV Pos", "Chromosomes CNV Neg")
write.csv(Swarbrick_CID4530N_ER_panCK_cnv_Chromosome_Quant, file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/Swarbrick_CID4530N_ER_panCK/Swarbrick_CID4530N_ER_panCK_cnv_Chromosome_Quant.csv")

# Save CNV Metadata
Swarbrick_CID4530N_ER_panCK_cnv.meta.data <- as.data.frame(as.matrix(Swarbrick_CID4530N_ER_panCK_cnv@meta.data))
write.csv(Swarbrick_CID4530N_ER_panCK_cnv.meta.data, file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/Swarbrick_CID4530N_ER_panCK/Swarbrick_CID4530N_ER_panCK_cnv.meta.data.csv")

# Subset Cells with inferred CNVs -------------------------------------------------------------------------------------------------------------
Swarbrick_CID4530N_ER_panCK_cnv <- subset(Swarbrick_CID4530N_ER_panCK_cnv, cells=rownames(Swarbrick_CID4530N_ER_panCK_cnv@meta.data)[rowSums(Swarbrick_CID4530N_ER_panCK_cnv@meta.data[,startsWith(names(Swarbrick_CID4530N_ER_panCK_cnv@meta.data),"has_cnv_")]) > 0 ])

# Save RDS with inferredCNVs
saveRDS(Swarbrick_CID4530N_ER_panCK_cnv, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_inferCNV/Swarbrick_CID4530N_ER_panCK_cnv.rds")


# Free memory
rm(Swarbrick_CID4530N_ER_inferCNV)
rm(Swarbrick_CID4530N_ER_Navin_Visvader_NORM_panCK_Reference)
rm(Swarbrick_CID4530N_ER_panCK_cnv)
rm(Swarbrick_CID4530N_ER_panCK_cnv_Pos_Chromosomes)
rm(Swarbrick_CID4530N_ER_panCK_cnv_Neg_Chromosomes)
rm(Swarbrick_CID4530N_ER_panCK_cnv_Chromosome_Quant)
rm(Swarbrick_CID4530N_ER_panCK_cnv.meta.data)
gc()



##############################################
#### InferCNV: Swarbrick_CID4535_ER_panCK ####
##############################################

########## Merge with normal mammary glad, extract epithelium, and process via Seurat
# Load RDS
Swarbrick_CID4535_ER_panCK <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_panCK/Swarbrick_CID4535_ER_panCK.rds")
Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Navin_Visvader_NORM_panCK_Reference.rds")

# Set identities 
colnames(x = Swarbrick_CID4535_ER_panCK[[]])
Idents(object = Swarbrick_CID4535_ER_panCK) <- 'Tumor'
Swarbrick_CID4535_ER_panCK[["cnv.ident"]] <- Idents(object = Swarbrick_CID4535_ER_panCK)
levels(x = Swarbrick_CID4535_ER_panCK)

colnames(x = Navin_Visvader_NORM_panCK_Reference[[]])
Idents(object = Navin_Visvader_NORM_panCK_Reference) <- 'Normal'
Navin_Visvader_NORM_panCK_Reference[["cnv.ident"]] <- Idents(object = Navin_Visvader_NORM_panCK_Reference)
levels(x = Navin_Visvader_NORM_panCK_Reference)

# Merge subsets into recombined RDS
Swarbrick_CID4535_ER_Navin_Visvader_NORM_panCK_Reference <- merge(x = Swarbrick_CID4535_ER_panCK, y = Navin_Visvader_NORM_panCK_Reference)
# Join Layers
Swarbrick_CID4535_ER_Navin_Visvader_NORM_panCK_Reference <- JoinLayers(Swarbrick_CID4535_ER_Navin_Visvader_NORM_panCK_Reference)

# Remove subsetted objects
rm(Swarbrick_CID4535_ER_panCK)
rm(Navin_Visvader_NORM_panCK_Reference)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples before merge

# Find Variable Features
Swarbrick_CID4535_ER_Navin_Visvader_NORM_panCK_Reference <- FindVariableFeatures(Swarbrick_CID4535_ER_Navin_Visvader_NORM_panCK_Reference, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Swarbrick_CID4535_ER_Navin_Visvader_NORM_panCK_Reference)
Swarbrick_CID4535_ER_Navin_Visvader_NORM_panCK_Reference <- ScaleData(Swarbrick_CID4535_ER_Navin_Visvader_NORM_panCK_Reference, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################

Swarbrick_CID4535_ER_Navin_Visvader_NORM_panCK_Reference <- RunPCA(Swarbrick_CID4535_ER_Navin_Visvader_NORM_panCK_Reference, verbose = FALSE)

# Examine and visualize PCA results a few different ways
print(Swarbrick_CID4535_ER_Navin_Visvader_NORM_panCK_Reference[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Swarbrick_CID4535_ER_Navin_Visvader_NORM_panCK_Reference, dims = 1:2, reduction = "pca")
DimPlot(Swarbrick_CID4535_ER_Navin_Visvader_NORM_panCK_Reference, reduction = "pca")

Swarbrick_CID4535_ER_Navin_Visvader_NORM_panCK_Reference <- FindNeighbors(Swarbrick_CID4535_ER_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
Swarbrick_CID4535_ER_Navin_Visvader_NORM_panCK_Reference <- FindClusters(Swarbrick_CID4535_ER_Navin_Visvader_NORM_panCK_Reference, resolution = 0.1)

# Umap clustering
Swarbrick_CID4535_ER_Navin_Visvader_NORM_panCK_Reference <- RunUMAP(Swarbrick_CID4535_ER_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
DimPlot(Swarbrick_CID4535_ER_Navin_Visvader_NORM_panCK_Reference, reduction = "umap")
DimPlot(Swarbrick_CID4535_ER_Navin_Visvader_NORM_panCK_Reference, reduction = "umap", group.by = "orig.ident")

# Save RDS
saveRDS(Swarbrick_CID4535_ER_Navin_Visvader_NORM_panCK_Reference, file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_inferCNV_Input_RDS/Swarbrick_CID4535_ER_Navin_Visvader_NORM_panCK_Reference.rds")

########## Prepare Cell Annotations
# Extract & Save Meta Data
Swarbrick_CID4535_ER_Navin_Visvader_NORM_panCK_Reference.meta.data <- as.data.frame(as.matrix(Swarbrick_CID4535_ER_Navin_Visvader_NORM_panCK_Reference@meta.data))
# Extract pertinent columns
# Make row names first column
Swarbrick_CID4535_ER_Navin_Visvader_NORM_panCK_Reference.meta.data <- tibble::rownames_to_column(Swarbrick_CID4535_ER_Navin_Visvader_NORM_panCK_Reference.meta.data, "Cells")
# In Seurat, "X" becomes the cell names while "cnv.ident" is specified above
Swarbrick_CID4535_ER_Navin_Visvader_NORM_panCK_Reference.annotations <- Swarbrick_CID4535_ER_Navin_Visvader_NORM_panCK_Reference.meta.data[, c("Cells", "cnv.ident")]
# Delete header
names(Swarbrick_CID4535_ER_Navin_Visvader_NORM_panCK_Reference.annotations) <- NULL
# Save as table
write.table(Swarbrick_CID4535_ER_Navin_Visvader_NORM_panCK_Reference.annotations,"/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_CID4535_ER_Navin_Visvader_NORM_panCK_Reference.annotations.txt",sep="\t",row.names=FALSE)
# Remove metadata and annotations
rm(Swarbrick_CID4535_ER_Navin_Visvader_NORM_panCK_Reference.meta.data)
rm(Swarbrick_CID4535_ER_Navin_Visvader_NORM_panCK_Reference.annotations)

########## Generate Count Matrix
# Extract Count Matrix
Swarbrick_CID4535_ER_Navin_Visvader_NORM_panCK_Reference.count <- Swarbrick_CID4535_ER_Navin_Visvader_NORM_panCK_Reference[["RNA"]]$counts# Convert counts to a sparse matrix
Swarbrick_CID4535_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as(Swarbrick_CID4535_ER_Navin_Visvader_NORM_panCK_Reference.count, 'sparseMatrix')
rm(Swarbrick_CID4535_ER_Navin_Visvader_NORM_panCK_Reference.count)
# Convert to requisite data format
Swarbrick_CID4535_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as.data.frame(Swarbrick_CID4535_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix)
# Save table
write.csv(Swarbrick_CID4535_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix, file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_inferCNV_Input_Count_Matrix/Swarbrick_CID4535_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv")
# Remove Objects
rm(Swarbrick_CID4535_ER_Navin_Visvader_NORM_panCK_Reference)
rm(Swarbrick_CID4535_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix)

########### Run inferCNV
# Note: check.names = FALSE prevents cell names from having dashes replaced with dots; essential to match annotations file
Swarbrick_CID4535_ER_inferCNV <- CreateInfercnvObject(raw_counts_matrix=read.csv(file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_inferCNV_Input_Count_Matrix/Swarbrick_CID4535_ER_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv", header = TRUE, row.names = 1,check.names = FALSE),
                                            annotations_file=read.table(file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_CID4535_ER_Navin_Visvader_NORM_panCK_Reference.annotations.txt", header = FALSE, row.names = 1),
                                            delim="\t",
                                            gene_order_file="/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/hg38_gencode_v27.txt",
                                            ref_group_names=c("Normal"))


# perform infercnv operations to reveal cnv signal
Swarbrick_CID4535_ER_inferCNV = infercnv::run(Swarbrick_CID4535_ER_inferCNV,
                                    cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                    out_dir="/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/Swarbrick_CID4535_ER_panCK/",  # dir is auto-created for storing outputs
                                    cluster_by_groups=T,   # cluster
                                    denoise=T,
                                    HMM=T,
                                    num_ref_groups = 9, analysis_mode = "subclusters")

###########  Integrate data with Seurat object
# Read RDS and add inferCNV results to metadata
Swarbrick_CID4535_ER_Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_inferCNV_Input_RDS/Swarbrick_CID4535_ER_Navin_Visvader_NORM_panCK_Reference.rds") 
Swarbrick_CID4535_ER_panCK_cnv <- infercnv::add_to_seurat(infercnv_output_path="/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/Swarbrick_CID4535_ER_panCK/",
                                                seurat_obj=Swarbrick_CID4535_ER_Navin_Visvader_NORM_panCK_Reference, # optional
                                                top_n=10)

# Collect Tumor cells from Seurat Object
Swarbrick_CID4535_ER_panCK_cnv <- subset(x = Swarbrick_CID4535_ER_panCK_cnv, subset = cnv.ident == "Tumor")
# Quantify chromosomes with alterations per cell
Swarbrick_CID4535_ER_panCK_cnv_Pos_Chromosomes <- as.data.frame(rowSums(Swarbrick_CID4535_ER_panCK_cnv@meta.data[,startsWith(names(Swarbrick_CID4535_ER_panCK_cnv@meta.data),"has_cnv_")] > 0 ))
Swarbrick_CID4535_ER_panCK_cnv_Neg_Chromosomes <- as.data.frame(rowSums(Swarbrick_CID4535_ER_panCK_cnv@meta.data[,startsWith(names(Swarbrick_CID4535_ER_panCK_cnv@meta.data),"has_cnv_")] == 0))
Swarbrick_CID4535_ER_panCK_cnv_Chromosome_Quant <- as.data.frame(c(Swarbrick_CID4535_ER_panCK_cnv_Pos_Chromosomes, Swarbrick_CID4535_ER_panCK_cnv_Neg_Chromosomes))
colnames(Swarbrick_CID4535_ER_panCK_cnv_Chromosome_Quant) <- c("Chromosomes CNV Pos", "Chromosomes CNV Neg")
write.csv(Swarbrick_CID4535_ER_panCK_cnv_Chromosome_Quant, file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/Swarbrick_CID4535_ER_panCK/Swarbrick_CID4535_ER_panCK_cnv_Chromosome_Quant.csv")

# Save CNV Metadata
Swarbrick_CID4535_ER_panCK_cnv.meta.data <- as.data.frame(as.matrix(Swarbrick_CID4535_ER_panCK_cnv@meta.data))
write.csv(Swarbrick_CID4535_ER_panCK_cnv.meta.data, file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/Swarbrick_CID4535_ER_panCK/Swarbrick_CID4535_ER_panCK_cnv.meta.data.csv")

# Subset Cells with inferred CNVs -------------------------------------------------------------------------------------------------------------
Swarbrick_CID4535_ER_panCK_cnv <- subset(Swarbrick_CID4535_ER_panCK_cnv, cells=rownames(Swarbrick_CID4535_ER_panCK_cnv@meta.data)[rowSums(Swarbrick_CID4535_ER_panCK_cnv@meta.data[,startsWith(names(Swarbrick_CID4535_ER_panCK_cnv@meta.data),"has_cnv_")]) > 0 ])

# Save RDS with inferredCNVs
saveRDS(Swarbrick_CID4535_ER_panCK_cnv, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_inferCNV/Swarbrick_CID4535_ER_panCK_cnv.rds")


# Free memory
rm(Swarbrick_CID4535_ER_inferCNV)
rm(Swarbrick_CID4535_ER_Navin_Visvader_NORM_panCK_Reference)
rm(Swarbrick_CID4535_ER_panCK_cnv)
rm(Swarbrick_CID4535_ER_panCK_cnv_Pos_Chromosomes)
rm(Swarbrick_CID4535_ER_panCK_cnv_Neg_Chromosomes)
rm(Swarbrick_CID4535_ER_panCK_cnv_Chromosome_Quant)
rm(Swarbrick_CID4535_ER_panCK_cnv.meta.data)
gc()



###################################################################################################
#############################     InferCNV: HER2 Tumor Specimens      #############################     
###################################################################################################

################################################
#### InferCNV: Swarbrick_CID3586_HER2_panCK ####
################################################

########## Merge with normal mammary glad, extract epithelium, and process via Seurat
# Load RDS
Swarbrick_CID3586_HER2_panCK <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_panCK/Swarbrick_CID3586_HER2_panCK.rds")
Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Navin_Visvader_NORM_panCK_Reference.rds")

# Set identities 
colnames(x = Swarbrick_CID3586_HER2_panCK[[]])
Idents(object = Swarbrick_CID3586_HER2_panCK) <- 'Tumor'
Swarbrick_CID3586_HER2_panCK[["cnv.ident"]] <- Idents(object = Swarbrick_CID3586_HER2_panCK)
levels(x = Swarbrick_CID3586_HER2_panCK)

colnames(x = Navin_Visvader_NORM_panCK_Reference[[]])
Idents(object = Navin_Visvader_NORM_panCK_Reference) <- 'Normal'
Navin_Visvader_NORM_panCK_Reference[["cnv.ident"]] <- Idents(object = Navin_Visvader_NORM_panCK_Reference)
levels(x = Navin_Visvader_NORM_panCK_Reference)

# Merge subsets into recombined RDS
Swarbrick_CID3586_HER2_Navin_Visvader_NORM_panCK_Reference <- merge(x = Swarbrick_CID3586_HER2_panCK, y = Navin_Visvader_NORM_panCK_Reference)
# Join Layers
Swarbrick_CID3586_HER2_Navin_Visvader_NORM_panCK_Reference <- JoinLayers(Swarbrick_CID3586_HER2_Navin_Visvader_NORM_panCK_Reference)

# Remove subsetted objects
rm(Swarbrick_CID3586_HER2_panCK)
rm(Navin_Visvader_NORM_panCK_Reference)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples before merge

# Find Variable Features
Swarbrick_CID3586_HER2_Navin_Visvader_NORM_panCK_Reference <- FindVariableFeatures(Swarbrick_CID3586_HER2_Navin_Visvader_NORM_panCK_Reference, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Swarbrick_CID3586_HER2_Navin_Visvader_NORM_panCK_Reference)
Swarbrick_CID3586_HER2_Navin_Visvader_NORM_panCK_Reference <- ScaleData(Swarbrick_CID3586_HER2_Navin_Visvader_NORM_panCK_Reference, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################

Swarbrick_CID3586_HER2_Navin_Visvader_NORM_panCK_Reference <- RunPCA(Swarbrick_CID3586_HER2_Navin_Visvader_NORM_panCK_Reference, verbose = FALSE)

# Examine and visualize PCA results a few different ways
print(Swarbrick_CID3586_HER2_Navin_Visvader_NORM_panCK_Reference[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Swarbrick_CID3586_HER2_Navin_Visvader_NORM_panCK_Reference, dims = 1:2, reduction = "pca")
DimPlot(Swarbrick_CID3586_HER2_Navin_Visvader_NORM_panCK_Reference, reduction = "pca")

Swarbrick_CID3586_HER2_Navin_Visvader_NORM_panCK_Reference <- FindNeighbors(Swarbrick_CID3586_HER2_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
Swarbrick_CID3586_HER2_Navin_Visvader_NORM_panCK_Reference <- FindClusters(Swarbrick_CID3586_HER2_Navin_Visvader_NORM_panCK_Reference, resolution = 0.1)

# Umap clustering
Swarbrick_CID3586_HER2_Navin_Visvader_NORM_panCK_Reference <- RunUMAP(Swarbrick_CID3586_HER2_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
DimPlot(Swarbrick_CID3586_HER2_Navin_Visvader_NORM_panCK_Reference, reduction = "umap")
DimPlot(Swarbrick_CID3586_HER2_Navin_Visvader_NORM_panCK_Reference, reduction = "umap", group.by = "orig.ident")

# Save RDS
saveRDS(Swarbrick_CID3586_HER2_Navin_Visvader_NORM_panCK_Reference, file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_inferCNV_Input_RDS/Swarbrick_CID3586_HER2_Navin_Visvader_NORM_panCK_Reference.rds")

########## Prepare Cell Annotations
# Extract & Save Meta Data
Swarbrick_CID3586_HER2_Navin_Visvader_NORM_panCK_Reference.meta.data <- as.data.frame(as.matrix(Swarbrick_CID3586_HER2_Navin_Visvader_NORM_panCK_Reference@meta.data))
# Extract pertinent columns
# Make row names first column
Swarbrick_CID3586_HER2_Navin_Visvader_NORM_panCK_Reference.meta.data <- tibble::rownames_to_column(Swarbrick_CID3586_HER2_Navin_Visvader_NORM_panCK_Reference.meta.data, "Cells")
# In Seurat, "X" becomes the cell names while "cnv.ident" is specified above
Swarbrick_CID3586_HER2_Navin_Visvader_NORM_panCK_Reference.annotations <- Swarbrick_CID3586_HER2_Navin_Visvader_NORM_panCK_Reference.meta.data[, c("Cells", "cnv.ident")]
# Delete header
names(Swarbrick_CID3586_HER2_Navin_Visvader_NORM_panCK_Reference.annotations) <- NULL
# Save as table
write.table(Swarbrick_CID3586_HER2_Navin_Visvader_NORM_panCK_Reference.annotations,"/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_CID3586_HER2_Navin_Visvader_NORM_panCK_Reference.annotations.txt",sep="\t",row.names=FALSE)
# Remove metadata and annotations
rm(Swarbrick_CID3586_HER2_Navin_Visvader_NORM_panCK_Reference.meta.data)
rm(Swarbrick_CID3586_HER2_Navin_Visvader_NORM_panCK_Reference.annotations)

########## Generate Count Matrix
# Extract Count Matrix
Swarbrick_CID3586_HER2_Navin_Visvader_NORM_panCK_Reference.count <- Swarbrick_CID3586_HER2_Navin_Visvader_NORM_panCK_Reference[["RNA"]]$counts
Swarbrick_CID3586_HER2_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as(Swarbrick_CID3586_HER2_Navin_Visvader_NORM_panCK_Reference.count, 'sparseMatrix')
rm(Swarbrick_CID3586_HER2_Navin_Visvader_NORM_panCK_Reference.count)
# Convert to requisite data format
Swarbrick_CID3586_HER2_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as.data.frame(Swarbrick_CID3586_HER2_Navin_Visvader_NORM_panCK_Reference.count.matrix)
# Save table
write.csv(Swarbrick_CID3586_HER2_Navin_Visvader_NORM_panCK_Reference.count.matrix, file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_inferCNV_Input_Count_Matrix/Swarbrick_CID3586_HER2_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv")
# Remove Objects
rm(Swarbrick_CID3586_HER2_Navin_Visvader_NORM_panCK_Reference)
rm(Swarbrick_CID3586_HER2_Navin_Visvader_NORM_panCK_Reference.count.matrix)

########### Run inferCNV
# Note: check.names = FALSE prevents cell names from having dashes replaced with dots; essential to match annotations file
Swarbrick_CID3586_HER2_inferCNV <- CreateInfercnvObject(raw_counts_matrix=read.csv(file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_inferCNV_Input_Count_Matrix/Swarbrick_CID3586_HER2_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv", header = TRUE, row.names = 1,check.names = FALSE),
                                              annotations_file=read.table(file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_CID3586_HER2_Navin_Visvader_NORM_panCK_Reference.annotations.txt", header = FALSE, row.names = 1),
                                              delim="\t",
                                              gene_order_file="/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/hg38_gencode_v27.txt",
                                              ref_group_names=c("Normal"))


# perform infercnv operations to reveal cnv signal
Swarbrick_CID3586_HER2_inferCNV = infercnv::run(Swarbrick_CID3586_HER2_inferCNV,
                                      cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                      out_dir="/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/Swarbrick_CID3586_HER2_panCK/",  # dir is auto-created for storing outputs
                                      cluster_by_groups=T,   # cluster
                                      denoise=T,
                                      HMM=T,
                                      num_ref_groups = 9, analysis_mode = "subclusters")

###########  Integrate data with Seurat object
# Read RDS and add inferCNV results to metadata
Swarbrick_CID3586_HER2_Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_inferCNV_Input_RDS/Swarbrick_CID3586_HER2_Navin_Visvader_NORM_panCK_Reference.rds") 
Swarbrick_CID3586_HER2_panCK_cnv <- infercnv::add_to_seurat(infercnv_output_path="/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/Swarbrick_CID3586_HER2_panCK/",
                                                  seurat_obj=Swarbrick_CID3586_HER2_Navin_Visvader_NORM_panCK_Reference, # optional
                                                  top_n=10)

# Collect Tumor cells from Seurat Object
Swarbrick_CID3586_HER2_panCK_cnv <- subset(x = Swarbrick_CID3586_HER2_panCK_cnv, subset = cnv.ident == "Tumor")
# Quantify chromosomes with alterations per cell
Swarbrick_CID3586_HER2_panCK_cnv_Pos_Chromosomes <- as.data.frame(rowSums(Swarbrick_CID3586_HER2_panCK_cnv@meta.data[,startsWith(names(Swarbrick_CID3586_HER2_panCK_cnv@meta.data),"has_cnv_")] > 0 ))
Swarbrick_CID3586_HER2_panCK_cnv_Neg_Chromosomes <- as.data.frame(rowSums(Swarbrick_CID3586_HER2_panCK_cnv@meta.data[,startsWith(names(Swarbrick_CID3586_HER2_panCK_cnv@meta.data),"has_cnv_")] == 0))
Swarbrick_CID3586_HER2_panCK_cnv_Chromosome_Quant <- as.data.frame(c(Swarbrick_CID3586_HER2_panCK_cnv_Pos_Chromosomes, Swarbrick_CID3586_HER2_panCK_cnv_Neg_Chromosomes))
colnames(Swarbrick_CID3586_HER2_panCK_cnv_Chromosome_Quant) <- c("Chromosomes CNV Pos", "Chromosomes CNV Neg")
write.csv(Swarbrick_CID3586_HER2_panCK_cnv_Chromosome_Quant, file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/Swarbrick_CID3586_HER2_panCK/Swarbrick_CID3586_HER2_panCK_cnv_Chromosome_Quant.csv")

# Save CNV Metadata
Swarbrick_CID3586_HER2_panCK_cnv.meta.data <- as.data.frame(as.matrix(Swarbrick_CID3586_HER2_panCK_cnv@meta.data))
write.csv(Swarbrick_CID3586_HER2_panCK_cnv.meta.data, file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/Swarbrick_CID3586_HER2_panCK/Swarbrick_CID3586_HER2_panCK_cnv.meta.data.csv")

# Subset Cells with inferred CNVs -------------------------------------------------------------------------------------------------------------
Swarbrick_CID3586_HER2_panCK_cnv <- subset(Swarbrick_CID3586_HER2_panCK_cnv, cells=rownames(Swarbrick_CID3586_HER2_panCK_cnv@meta.data)[rowSums(Swarbrick_CID3586_HER2_panCK_cnv@meta.data[,startsWith(names(Swarbrick_CID3586_HER2_panCK_cnv@meta.data),"has_cnv_")]) > 0 ])

# Save RDS with inferredCNVs
saveRDS(Swarbrick_CID3586_HER2_panCK_cnv, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_inferCNV/Swarbrick_CID3586_HER2_panCK_cnv.rds")

# Free memory
rm(Swarbrick_CID3586_HER2_inferCNV)
rm(Swarbrick_CID3586_HER2_Navin_Visvader_NORM_panCK_Reference)
rm(Swarbrick_CID3586_HER2_panCK_cnv)
rm(Swarbrick_CID3586_HER2_panCK_cnv_Pos_Chromosomes)
rm(Swarbrick_CID3586_HER2_panCK_cnv_Neg_Chromosomes)
rm(Swarbrick_CID3586_HER2_panCK_cnv_Chromosome_Quant)
rm(Swarbrick_CID3586_HER2_panCK_cnv.meta.data)
gc()


################################################
#### InferCNV: Swarbrick_CID3838_HER2_panCK ####
################################################

########## Merge with normal mammary glad, extract epithelium, and process via Seurat
# Load RDS
Swarbrick_CID3838_HER2_panCK <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_panCK/Swarbrick_CID3838_HER2_panCK.rds")
Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Navin_Visvader_NORM_panCK_Reference.rds")

# Set identities 
colnames(x = Swarbrick_CID3838_HER2_panCK[[]])
Idents(object = Swarbrick_CID3838_HER2_panCK) <- 'Tumor'
Swarbrick_CID3838_HER2_panCK[["cnv.ident"]] <- Idents(object = Swarbrick_CID3838_HER2_panCK)
levels(x = Swarbrick_CID3838_HER2_panCK)

colnames(x = Navin_Visvader_NORM_panCK_Reference[[]])
Idents(object = Navin_Visvader_NORM_panCK_Reference) <- 'Normal'
Navin_Visvader_NORM_panCK_Reference[["cnv.ident"]] <- Idents(object = Navin_Visvader_NORM_panCK_Reference)
levels(x = Navin_Visvader_NORM_panCK_Reference)

# Merge subsets into recombined RDS
Swarbrick_CID3838_HER2_Navin_Visvader_NORM_panCK_Reference <- merge(x = Swarbrick_CID3838_HER2_panCK, y = Navin_Visvader_NORM_panCK_Reference)
# Join Layers
Swarbrick_CID3838_HER2_Navin_Visvader_NORM_panCK_Reference <- JoinLayers(Swarbrick_CID3838_HER2_Navin_Visvader_NORM_panCK_Reference)

# Remove subsetted objects
rm(Swarbrick_CID3838_HER2_panCK)
rm(Navin_Visvader_NORM_panCK_Reference)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples before merge

# Find Variable Features
Swarbrick_CID3838_HER2_Navin_Visvader_NORM_panCK_Reference <- FindVariableFeatures(Swarbrick_CID3838_HER2_Navin_Visvader_NORM_panCK_Reference, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Swarbrick_CID3838_HER2_Navin_Visvader_NORM_panCK_Reference)
Swarbrick_CID3838_HER2_Navin_Visvader_NORM_panCK_Reference <- ScaleData(Swarbrick_CID3838_HER2_Navin_Visvader_NORM_panCK_Reference, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################

Swarbrick_CID3838_HER2_Navin_Visvader_NORM_panCK_Reference <- RunPCA(Swarbrick_CID3838_HER2_Navin_Visvader_NORM_panCK_Reference, verbose = FALSE)

# Examine and visualize PCA results a few different ways
print(Swarbrick_CID3838_HER2_Navin_Visvader_NORM_panCK_Reference[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Swarbrick_CID3838_HER2_Navin_Visvader_NORM_panCK_Reference, dims = 1:2, reduction = "pca")
DimPlot(Swarbrick_CID3838_HER2_Navin_Visvader_NORM_panCK_Reference, reduction = "pca")

Swarbrick_CID3838_HER2_Navin_Visvader_NORM_panCK_Reference <- FindNeighbors(Swarbrick_CID3838_HER2_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
Swarbrick_CID3838_HER2_Navin_Visvader_NORM_panCK_Reference <- FindClusters(Swarbrick_CID3838_HER2_Navin_Visvader_NORM_panCK_Reference, resolution = 0.1)

# Umap clustering
Swarbrick_CID3838_HER2_Navin_Visvader_NORM_panCK_Reference <- RunUMAP(Swarbrick_CID3838_HER2_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
DimPlot(Swarbrick_CID3838_HER2_Navin_Visvader_NORM_panCK_Reference, reduction = "umap")
DimPlot(Swarbrick_CID3838_HER2_Navin_Visvader_NORM_panCK_Reference, reduction = "umap", group.by = "orig.ident")

# Save RDS
saveRDS(Swarbrick_CID3838_HER2_Navin_Visvader_NORM_panCK_Reference, file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_inferCNV_Input_RDS/Swarbrick_CID3838_HER2_Navin_Visvader_NORM_panCK_Reference.rds")

########## Prepare Cell Annotations
# Extract & Save Meta Data
Swarbrick_CID3838_HER2_Navin_Visvader_NORM_panCK_Reference.meta.data <- as.data.frame(as.matrix(Swarbrick_CID3838_HER2_Navin_Visvader_NORM_panCK_Reference@meta.data))
# Extract pertinent columns
# Make row names first column
Swarbrick_CID3838_HER2_Navin_Visvader_NORM_panCK_Reference.meta.data <- tibble::rownames_to_column(Swarbrick_CID3838_HER2_Navin_Visvader_NORM_panCK_Reference.meta.data, "Cells")
# In Seurat, "X" becomes the cell names while "cnv.ident" is specified above
Swarbrick_CID3838_HER2_Navin_Visvader_NORM_panCK_Reference.annotations <- Swarbrick_CID3838_HER2_Navin_Visvader_NORM_panCK_Reference.meta.data[, c("Cells", "cnv.ident")]
# Delete header
names(Swarbrick_CID3838_HER2_Navin_Visvader_NORM_panCK_Reference.annotations) <- NULL
# Save as table
write.table(Swarbrick_CID3838_HER2_Navin_Visvader_NORM_panCK_Reference.annotations,"/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_CID3838_HER2_Navin_Visvader_NORM_panCK_Reference.annotations.txt",sep="\t",row.names=FALSE)
# Remove metadata and annotations
rm(Swarbrick_CID3838_HER2_Navin_Visvader_NORM_panCK_Reference.meta.data)
rm(Swarbrick_CID3838_HER2_Navin_Visvader_NORM_panCK_Reference.annotations)

########## Generate Count Matrix
# Extract Count Matrix
Swarbrick_CID3838_HER2_Navin_Visvader_NORM_panCK_Reference.count <- Swarbrick_CID3838_HER2_Navin_Visvader_NORM_panCK_Reference[["RNA"]]$counts
Swarbrick_CID3838_HER2_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as(Swarbrick_CID3838_HER2_Navin_Visvader_NORM_panCK_Reference.count, 'sparseMatrix')
rm(Swarbrick_CID3838_HER2_Navin_Visvader_NORM_panCK_Reference.count)
# Convert to requisite data format
Swarbrick_CID3838_HER2_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as.data.frame(Swarbrick_CID3838_HER2_Navin_Visvader_NORM_panCK_Reference.count.matrix)
# Save table
write.csv(Swarbrick_CID3838_HER2_Navin_Visvader_NORM_panCK_Reference.count.matrix, file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_inferCNV_Input_Count_Matrix/Swarbrick_CID3838_HER2_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv")
# Remove Objects
rm(Swarbrick_CID3838_HER2_Navin_Visvader_NORM_panCK_Reference)
rm(Swarbrick_CID3838_HER2_Navin_Visvader_NORM_panCK_Reference.count.matrix)

########### Run inferCNV
# Note: check.names = FALSE prevents cell names from having dashes replaced with dots; essential to match annotations file
Swarbrick_CID3838_HER2_inferCNV <- CreateInfercnvObject(raw_counts_matrix=read.csv(file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_inferCNV_Input_Count_Matrix/Swarbrick_CID3838_HER2_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv", header = TRUE, row.names = 1,check.names = FALSE),
                                              annotations_file=read.table(file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_CID3838_HER2_Navin_Visvader_NORM_panCK_Reference.annotations.txt", header = FALSE, row.names = 1),
                                              delim="\t",
                                              gene_order_file="/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/hg38_gencode_v27.txt",
                                              ref_group_names=c("Normal"))


# perform infercnv operations to reveal cnv signal
Swarbrick_CID3838_HER2_inferCNV = infercnv::run(Swarbrick_CID3838_HER2_inferCNV,
                                      cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                      out_dir="/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/Swarbrick_CID3838_HER2_panCK/",  # dir is auto-created for storing outputs
                                      cluster_by_groups=T,   # cluster
                                      denoise=T,
                                      HMM=T,
                                      num_ref_groups = 9, analysis_mode = "subclusters")

###########  Integrate data with Seurat object
# Read RDS and add inferCNV results to metadata
Swarbrick_CID3838_HER2_Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_inferCNV_Input_RDS/Swarbrick_CID3838_HER2_Navin_Visvader_NORM_panCK_Reference.rds") 
Swarbrick_CID3838_HER2_panCK_cnv <- infercnv::add_to_seurat(infercnv_output_path="/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/Swarbrick_CID3838_HER2_panCK/",
                                                  seurat_obj=Swarbrick_CID3838_HER2_Navin_Visvader_NORM_panCK_Reference, # optional
                                                  top_n=10)

# Collect Tumor cells from Seurat Object
Swarbrick_CID3838_HER2_panCK_cnv <- subset(x = Swarbrick_CID3838_HER2_panCK_cnv, subset = cnv.ident == "Tumor")
# Quantify chromosomes with alterations per cell
Swarbrick_CID3838_HER2_panCK_cnv_Pos_Chromosomes <- as.data.frame(rowSums(Swarbrick_CID3838_HER2_panCK_cnv@meta.data[,startsWith(names(Swarbrick_CID3838_HER2_panCK_cnv@meta.data),"has_cnv_")] > 0 ))
Swarbrick_CID3838_HER2_panCK_cnv_Neg_Chromosomes <- as.data.frame(rowSums(Swarbrick_CID3838_HER2_panCK_cnv@meta.data[,startsWith(names(Swarbrick_CID3838_HER2_panCK_cnv@meta.data),"has_cnv_")] == 0))
Swarbrick_CID3838_HER2_panCK_cnv_Chromosome_Quant <- as.data.frame(c(Swarbrick_CID3838_HER2_panCK_cnv_Pos_Chromosomes, Swarbrick_CID3838_HER2_panCK_cnv_Neg_Chromosomes))
colnames(Swarbrick_CID3838_HER2_panCK_cnv_Chromosome_Quant) <- c("Chromosomes CNV Pos", "Chromosomes CNV Neg")
write.csv(Swarbrick_CID3838_HER2_panCK_cnv_Chromosome_Quant, file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/Swarbrick_CID3838_HER2_panCK/Swarbrick_CID3838_HER2_panCK_cnv_Chromosome_Quant.csv")

# Save CNV Metadata
Swarbrick_CID3838_HER2_panCK_cnv.meta.data <- as.data.frame(as.matrix(Swarbrick_CID3838_HER2_panCK_cnv@meta.data))
write.csv(Swarbrick_CID3838_HER2_panCK_cnv.meta.data, file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/Swarbrick_CID3838_HER2_panCK/Swarbrick_CID3838_HER2_panCK_cnv.meta.data.csv")

# Subset Cells with inferred CNVs -------------------------------------------------------------------------------------------------------------
Swarbrick_CID3838_HER2_panCK_cnv <- subset(Swarbrick_CID3838_HER2_panCK_cnv, cells=rownames(Swarbrick_CID3838_HER2_panCK_cnv@meta.data)[rowSums(Swarbrick_CID3838_HER2_panCK_cnv@meta.data[,startsWith(names(Swarbrick_CID3838_HER2_panCK_cnv@meta.data),"has_cnv_")]) > 0 ])

# Save RDS with inferredCNVs
saveRDS(Swarbrick_CID3838_HER2_panCK_cnv, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_inferCNV/Swarbrick_CID3838_HER2_panCK_cnv.rds")


# Free memory
rm(Swarbrick_CID3838_HER2_inferCNV)
rm(Swarbrick_CID3838_HER2_Navin_Visvader_NORM_panCK_Reference)
rm(Swarbrick_CID3838_HER2_panCK_cnv)
rm(Swarbrick_CID3838_HER2_panCK_cnv_Pos_Chromosomes)
rm(Swarbrick_CID3838_HER2_panCK_cnv_Neg_Chromosomes)
rm(Swarbrick_CID3838_HER2_panCK_cnv_Chromosome_Quant)
rm(Swarbrick_CID3838_HER2_panCK_cnv.meta.data)
gc()


################################################
#### InferCNV: Swarbrick_CID3921_HER2_panCK ####
################################################

########## Merge with normal mammary glad, extract epithelium, and process via Seurat
# Load RDS
Swarbrick_CID3921_HER2_panCK <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_panCK/Swarbrick_CID3921_HER2_panCK.rds")
Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Navin_Visvader_NORM_panCK_Reference.rds")

# Set identities 
colnames(x = Swarbrick_CID3921_HER2_panCK[[]])
Idents(object = Swarbrick_CID3921_HER2_panCK) <- 'Tumor'
Swarbrick_CID3921_HER2_panCK[["cnv.ident"]] <- Idents(object = Swarbrick_CID3921_HER2_panCK)
levels(x = Swarbrick_CID3921_HER2_panCK)

colnames(x = Navin_Visvader_NORM_panCK_Reference[[]])
Idents(object = Navin_Visvader_NORM_panCK_Reference) <- 'Normal'
Navin_Visvader_NORM_panCK_Reference[["cnv.ident"]] <- Idents(object = Navin_Visvader_NORM_panCK_Reference)
levels(x = Navin_Visvader_NORM_panCK_Reference)

# Merge subsets into recombined RDS
Swarbrick_CID3921_HER2_Navin_Visvader_NORM_panCK_Reference <- merge(x = Swarbrick_CID3921_HER2_panCK, y = Navin_Visvader_NORM_panCK_Reference)
# Join Layers
Swarbrick_CID3921_HER2_Navin_Visvader_NORM_panCK_Reference <- JoinLayers(Swarbrick_CID3921_HER2_Navin_Visvader_NORM_panCK_Reference)

# Remove subsetted objects
rm(Swarbrick_CID3921_HER2_panCK)
rm(Navin_Visvader_NORM_panCK_Reference)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples before merge

# Find Variable Features
Swarbrick_CID3921_HER2_Navin_Visvader_NORM_panCK_Reference <- FindVariableFeatures(Swarbrick_CID3921_HER2_Navin_Visvader_NORM_panCK_Reference, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Swarbrick_CID3921_HER2_Navin_Visvader_NORM_panCK_Reference)
Swarbrick_CID3921_HER2_Navin_Visvader_NORM_panCK_Reference <- ScaleData(Swarbrick_CID3921_HER2_Navin_Visvader_NORM_panCK_Reference, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################

Swarbrick_CID3921_HER2_Navin_Visvader_NORM_panCK_Reference <- RunPCA(Swarbrick_CID3921_HER2_Navin_Visvader_NORM_panCK_Reference, verbose = FALSE)

# Examine and visualize PCA results a few different ways
print(Swarbrick_CID3921_HER2_Navin_Visvader_NORM_panCK_Reference[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Swarbrick_CID3921_HER2_Navin_Visvader_NORM_panCK_Reference, dims = 1:2, reduction = "pca")
DimPlot(Swarbrick_CID3921_HER2_Navin_Visvader_NORM_panCK_Reference, reduction = "pca")

Swarbrick_CID3921_HER2_Navin_Visvader_NORM_panCK_Reference <- FindNeighbors(Swarbrick_CID3921_HER2_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
Swarbrick_CID3921_HER2_Navin_Visvader_NORM_panCK_Reference <- FindClusters(Swarbrick_CID3921_HER2_Navin_Visvader_NORM_panCK_Reference, resolution = 0.1)

# Umap clustering
Swarbrick_CID3921_HER2_Navin_Visvader_NORM_panCK_Reference <- RunUMAP(Swarbrick_CID3921_HER2_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
DimPlot(Swarbrick_CID3921_HER2_Navin_Visvader_NORM_panCK_Reference, reduction = "umap")
DimPlot(Swarbrick_CID3921_HER2_Navin_Visvader_NORM_panCK_Reference, reduction = "umap", group.by = "orig.ident")

# Save RDS
saveRDS(Swarbrick_CID3921_HER2_Navin_Visvader_NORM_panCK_Reference, file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_inferCNV_Input_RDS/Swarbrick_CID3921_HER2_Navin_Visvader_NORM_panCK_Reference.rds")

########## Prepare Cell Annotations
# Extract & Save Meta Data
Swarbrick_CID3921_HER2_Navin_Visvader_NORM_panCK_Reference.meta.data <- as.data.frame(as.matrix(Swarbrick_CID3921_HER2_Navin_Visvader_NORM_panCK_Reference@meta.data))
# Extract pertinent columns
# Make row names first column
Swarbrick_CID3921_HER2_Navin_Visvader_NORM_panCK_Reference.meta.data <- tibble::rownames_to_column(Swarbrick_CID3921_HER2_Navin_Visvader_NORM_panCK_Reference.meta.data, "Cells")
# In Seurat, "X" becomes the cell names while "cnv.ident" is specified above
Swarbrick_CID3921_HER2_Navin_Visvader_NORM_panCK_Reference.annotations <- Swarbrick_CID3921_HER2_Navin_Visvader_NORM_panCK_Reference.meta.data[, c("Cells", "cnv.ident")]
# Delete header
names(Swarbrick_CID3921_HER2_Navin_Visvader_NORM_panCK_Reference.annotations) <- NULL
# Save as table
write.table(Swarbrick_CID3921_HER2_Navin_Visvader_NORM_panCK_Reference.annotations,"/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_CID3921_HER2_Navin_Visvader_NORM_panCK_Reference.annotations.txt",sep="\t",row.names=FALSE)
# Remove metadata and annotations
rm(Swarbrick_CID3921_HER2_Navin_Visvader_NORM_panCK_Reference.meta.data)
rm(Swarbrick_CID3921_HER2_Navin_Visvader_NORM_panCK_Reference.annotations)

########## Generate Count Matrix
# Extract Count Matrix
Swarbrick_CID3921_HER2_Navin_Visvader_NORM_panCK_Reference.count <- Swarbrick_CID3921_HER2_Navin_Visvader_NORM_panCK_Reference[["RNA"]]$counts
Swarbrick_CID3921_HER2_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as(Swarbrick_CID3921_HER2_Navin_Visvader_NORM_panCK_Reference.count, 'sparseMatrix')
rm(Swarbrick_CID3921_HER2_Navin_Visvader_NORM_panCK_Reference.count)
# Convert to requisite data format
Swarbrick_CID3921_HER2_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as.data.frame(Swarbrick_CID3921_HER2_Navin_Visvader_NORM_panCK_Reference.count.matrix)
# Save table
write.csv(Swarbrick_CID3921_HER2_Navin_Visvader_NORM_panCK_Reference.count.matrix, file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_inferCNV_Input_Count_Matrix/Swarbrick_CID3921_HER2_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv")
# Remove Objects
rm(Swarbrick_CID3921_HER2_Navin_Visvader_NORM_panCK_Reference)
rm(Swarbrick_CID3921_HER2_Navin_Visvader_NORM_panCK_Reference.count.matrix)

########### Run inferCNV
# Note: check.names = FALSE prevents cell names from having dashes replaced with dots; essential to match annotations file
Swarbrick_CID3921_HER2_inferCNV <- CreateInfercnvObject(raw_counts_matrix=read.csv(file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_inferCNV_Input_Count_Matrix/Swarbrick_CID3921_HER2_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv", header = TRUE, row.names = 1,check.names = FALSE),
                                              annotations_file=read.table(file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_CID3921_HER2_Navin_Visvader_NORM_panCK_Reference.annotations.txt", header = FALSE, row.names = 1),
                                              delim="\t",
                                              gene_order_file="/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/hg38_gencode_v27.txt",
                                              ref_group_names=c("Normal"))


# perform infercnv operations to reveal cnv signal
Swarbrick_CID3921_HER2_inferCNV = infercnv::run(Swarbrick_CID3921_HER2_inferCNV,
                                      cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                      out_dir="/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/Swarbrick_CID3921_HER2_panCK/",  # dir is auto-created for storing outputs
                                      cluster_by_groups=T,   # cluster
                                      denoise=T,
                                      HMM=T,
                                      num_ref_groups = 9, analysis_mode = "subclusters")

###########  Integrate data with Seurat object
# Read RDS and add inferCNV results to metadata
Swarbrick_CID3921_HER2_Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_inferCNV_Input_RDS/Swarbrick_CID3921_HER2_Navin_Visvader_NORM_panCK_Reference.rds") 
Swarbrick_CID3921_HER2_panCK_cnv <- infercnv::add_to_seurat(infercnv_output_path="/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/Swarbrick_CID3921_HER2_panCK/",
                                                  seurat_obj=Swarbrick_CID3921_HER2_Navin_Visvader_NORM_panCK_Reference, # optional
                                                  top_n=10)

# Collect Tumor cells from Seurat Object
Swarbrick_CID3921_HER2_panCK_cnv <- subset(x = Swarbrick_CID3921_HER2_panCK_cnv, subset = cnv.ident == "Tumor")
# Quantify chromosomes with alterations per cell
Swarbrick_CID3921_HER2_panCK_cnv_Pos_Chromosomes <- as.data.frame(rowSums(Swarbrick_CID3921_HER2_panCK_cnv@meta.data[,startsWith(names(Swarbrick_CID3921_HER2_panCK_cnv@meta.data),"has_cnv_")] > 0 ))
Swarbrick_CID3921_HER2_panCK_cnv_Neg_Chromosomes <- as.data.frame(rowSums(Swarbrick_CID3921_HER2_panCK_cnv@meta.data[,startsWith(names(Swarbrick_CID3921_HER2_panCK_cnv@meta.data),"has_cnv_")] == 0))
Swarbrick_CID3921_HER2_panCK_cnv_Chromosome_Quant <- as.data.frame(c(Swarbrick_CID3921_HER2_panCK_cnv_Pos_Chromosomes, Swarbrick_CID3921_HER2_panCK_cnv_Neg_Chromosomes))
colnames(Swarbrick_CID3921_HER2_panCK_cnv_Chromosome_Quant) <- c("Chromosomes CNV Pos", "Chromosomes CNV Neg")
write.csv(Swarbrick_CID3921_HER2_panCK_cnv_Chromosome_Quant, file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/Swarbrick_CID3921_HER2_panCK/Swarbrick_CID3921_HER2_panCK_cnv_Chromosome_Quant.csv")

# Save CNV Metadata
Swarbrick_CID3921_HER2_panCK_cnv.meta.data <- as.data.frame(as.matrix(Swarbrick_CID3921_HER2_panCK_cnv@meta.data))
write.csv(Swarbrick_CID3921_HER2_panCK_cnv.meta.data, file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/Swarbrick_CID3921_HER2_panCK/Swarbrick_CID3921_HER2_panCK_cnv.meta.data.csv")

# Subset Cells with inferred CNVs -------------------------------------------------------------------------------------------------------------
Swarbrick_CID3921_HER2_panCK_cnv <- subset(Swarbrick_CID3921_HER2_panCK_cnv, cells=rownames(Swarbrick_CID3921_HER2_panCK_cnv@meta.data)[rowSums(Swarbrick_CID3921_HER2_panCK_cnv@meta.data[,startsWith(names(Swarbrick_CID3921_HER2_panCK_cnv@meta.data),"has_cnv_")]) > 0 ])

# Save RDS with inferredCNVs
saveRDS(Swarbrick_CID3921_HER2_panCK_cnv, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_inferCNV/Swarbrick_CID3921_HER2_panCK_cnv.rds")


# Free memory
rm(Swarbrick_CID3921_HER2_inferCNV)
rm(Swarbrick_CID3921_HER2_Navin_Visvader_NORM_panCK_Reference)
rm(Swarbrick_CID3921_HER2_panCK_cnv)
rm(Swarbrick_CID3921_HER2_panCK_cnv_Pos_Chromosomes)
rm(Swarbrick_CID3921_HER2_panCK_cnv_Neg_Chromosomes)
rm(Swarbrick_CID3921_HER2_panCK_cnv_Chromosome_Quant)
rm(Swarbrick_CID3921_HER2_panCK_cnv.meta.data)
gc()


################################################
#### InferCNV: Swarbrick_CID4066_HER2_panCK ####
################################################

########## Merge with normal mammary glad, extract epithelium, and process via Seurat
# Load RDS
Swarbrick_CID4066_HER2_panCK <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_panCK/Swarbrick_CID4066_HER2_panCK.rds")
Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Navin_Visvader_NORM_panCK_Reference.rds")

# Set identities 
colnames(x = Swarbrick_CID4066_HER2_panCK[[]])
Idents(object = Swarbrick_CID4066_HER2_panCK) <- 'Tumor'
Swarbrick_CID4066_HER2_panCK[["cnv.ident"]] <- Idents(object = Swarbrick_CID4066_HER2_panCK)
levels(x = Swarbrick_CID4066_HER2_panCK)

colnames(x = Navin_Visvader_NORM_panCK_Reference[[]])
Idents(object = Navin_Visvader_NORM_panCK_Reference) <- 'Normal'
Navin_Visvader_NORM_panCK_Reference[["cnv.ident"]] <- Idents(object = Navin_Visvader_NORM_panCK_Reference)
levels(x = Navin_Visvader_NORM_panCK_Reference)

# Merge subsets into recombined RDS
Swarbrick_CID4066_HER2_Navin_Visvader_NORM_panCK_Reference <- merge(x = Swarbrick_CID4066_HER2_panCK, y = Navin_Visvader_NORM_panCK_Reference)
# Join Layers
Swarbrick_CID4066_HER2_Navin_Visvader_NORM_panCK_Reference <- JoinLayers(Swarbrick_CID4066_HER2_Navin_Visvader_NORM_panCK_Reference)

# Remove subsetted objects
rm(Swarbrick_CID4066_HER2_panCK)
rm(Navin_Visvader_NORM_panCK_Reference)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples before merge

# Find Variable Features
Swarbrick_CID4066_HER2_Navin_Visvader_NORM_panCK_Reference <- FindVariableFeatures(Swarbrick_CID4066_HER2_Navin_Visvader_NORM_panCK_Reference, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Swarbrick_CID4066_HER2_Navin_Visvader_NORM_panCK_Reference)
Swarbrick_CID4066_HER2_Navin_Visvader_NORM_panCK_Reference <- ScaleData(Swarbrick_CID4066_HER2_Navin_Visvader_NORM_panCK_Reference, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################

Swarbrick_CID4066_HER2_Navin_Visvader_NORM_panCK_Reference <- RunPCA(Swarbrick_CID4066_HER2_Navin_Visvader_NORM_panCK_Reference, verbose = FALSE)

# Examine and visualize PCA results a few different ways
print(Swarbrick_CID4066_HER2_Navin_Visvader_NORM_panCK_Reference[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Swarbrick_CID4066_HER2_Navin_Visvader_NORM_panCK_Reference, dims = 1:2, reduction = "pca")
DimPlot(Swarbrick_CID4066_HER2_Navin_Visvader_NORM_panCK_Reference, reduction = "pca")

Swarbrick_CID4066_HER2_Navin_Visvader_NORM_panCK_Reference <- FindNeighbors(Swarbrick_CID4066_HER2_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
Swarbrick_CID4066_HER2_Navin_Visvader_NORM_panCK_Reference <- FindClusters(Swarbrick_CID4066_HER2_Navin_Visvader_NORM_panCK_Reference, resolution = 0.1)

# Umap clustering
Swarbrick_CID4066_HER2_Navin_Visvader_NORM_panCK_Reference <- RunUMAP(Swarbrick_CID4066_HER2_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
DimPlot(Swarbrick_CID4066_HER2_Navin_Visvader_NORM_panCK_Reference, reduction = "umap")
DimPlot(Swarbrick_CID4066_HER2_Navin_Visvader_NORM_panCK_Reference, reduction = "umap", group.by = "orig.ident")

# Save RDS
saveRDS(Swarbrick_CID4066_HER2_Navin_Visvader_NORM_panCK_Reference, file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_inferCNV_Input_RDS/Swarbrick_CID4066_HER2_Navin_Visvader_NORM_panCK_Reference.rds")

########## Prepare Cell Annotations
# Extract & Save Meta Data
Swarbrick_CID4066_HER2_Navin_Visvader_NORM_panCK_Reference.meta.data <- as.data.frame(as.matrix(Swarbrick_CID4066_HER2_Navin_Visvader_NORM_panCK_Reference@meta.data))
# Extract pertinent columns
# Make row names first column
Swarbrick_CID4066_HER2_Navin_Visvader_NORM_panCK_Reference.meta.data <- tibble::rownames_to_column(Swarbrick_CID4066_HER2_Navin_Visvader_NORM_panCK_Reference.meta.data, "Cells")
# In Seurat, "X" becomes the cell names while "cnv.ident" is specified above
Swarbrick_CID4066_HER2_Navin_Visvader_NORM_panCK_Reference.annotations <- Swarbrick_CID4066_HER2_Navin_Visvader_NORM_panCK_Reference.meta.data[, c("Cells", "cnv.ident")]
# Delete header
names(Swarbrick_CID4066_HER2_Navin_Visvader_NORM_panCK_Reference.annotations) <- NULL
# Save as table
write.table(Swarbrick_CID4066_HER2_Navin_Visvader_NORM_panCK_Reference.annotations,"/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_CID4066_HER2_Navin_Visvader_NORM_panCK_Reference.annotations.txt",sep="\t",row.names=FALSE)
# Remove metadata and annotations
rm(Swarbrick_CID4066_HER2_Navin_Visvader_NORM_panCK_Reference.meta.data)
rm(Swarbrick_CID4066_HER2_Navin_Visvader_NORM_panCK_Reference.annotations)

########## Generate Count Matrix
# Extract Count Matrix
Swarbrick_CID4066_HER2_Navin_Visvader_NORM_panCK_Reference.count <- Swarbrick_CID4066_HER2_Navin_Visvader_NORM_panCK_Reference[["RNA"]]$counts
Swarbrick_CID4066_HER2_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as(Swarbrick_CID4066_HER2_Navin_Visvader_NORM_panCK_Reference.count, 'sparseMatrix')
rm(Swarbrick_CID4066_HER2_Navin_Visvader_NORM_panCK_Reference.count)
# Convert to requisite data format
Swarbrick_CID4066_HER2_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as.data.frame(Swarbrick_CID4066_HER2_Navin_Visvader_NORM_panCK_Reference.count.matrix)
# Save table
write.csv(Swarbrick_CID4066_HER2_Navin_Visvader_NORM_panCK_Reference.count.matrix, file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_inferCNV_Input_Count_Matrix/Swarbrick_CID4066_HER2_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv")
# Remove Objects
rm(Swarbrick_CID4066_HER2_Navin_Visvader_NORM_panCK_Reference)
rm(Swarbrick_CID4066_HER2_Navin_Visvader_NORM_panCK_Reference.count.matrix)

########### Run inferCNV
# Note: check.names = FALSE prevents cell names from having dashes replaced with dots; essential to match annotations file
Swarbrick_CID4066_HER2_inferCNV <- CreateInfercnvObject(raw_counts_matrix=read.csv(file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_inferCNV_Input_Count_Matrix/Swarbrick_CID4066_HER2_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv", header = TRUE, row.names = 1,check.names = FALSE),
                                              annotations_file=read.table(file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_CID4066_HER2_Navin_Visvader_NORM_panCK_Reference.annotations.txt", header = FALSE, row.names = 1),
                                              delim="\t",
                                              gene_order_file="/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/hg38_gencode_v27.txt",
                                              ref_group_names=c("Normal"))


# perform infercnv operations to reveal cnv signal
Swarbrick_CID4066_HER2_inferCNV = infercnv::run(Swarbrick_CID4066_HER2_inferCNV,
                                      cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                      out_dir="/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/Swarbrick_CID4066_HER2_panCK/",  # dir is auto-created for storing outputs
                                      cluster_by_groups=T,   # cluster
                                      denoise=T,
                                      HMM=T,
                                      num_ref_groups = 9, analysis_mode = "subclusters")

###########  Integrate data with Seurat object
# Read RDS and add inferCNV results to metadata
Swarbrick_CID4066_HER2_Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_inferCNV_Input_RDS/Swarbrick_CID4066_HER2_Navin_Visvader_NORM_panCK_Reference.rds") 
Swarbrick_CID4066_HER2_panCK_cnv <- infercnv::add_to_seurat(infercnv_output_path="/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/Swarbrick_CID4066_HER2_panCK/",
                                                  seurat_obj=Swarbrick_CID4066_HER2_Navin_Visvader_NORM_panCK_Reference, # optional
                                                  top_n=10)

# Collect Tumor cells from Seurat Object
Swarbrick_CID4066_HER2_panCK_cnv <- subset(x = Swarbrick_CID4066_HER2_panCK_cnv, subset = cnv.ident == "Tumor")
# Quantify chromosomes with alterations per cell
Swarbrick_CID4066_HER2_panCK_cnv_Pos_Chromosomes <- as.data.frame(rowSums(Swarbrick_CID4066_HER2_panCK_cnv@meta.data[,startsWith(names(Swarbrick_CID4066_HER2_panCK_cnv@meta.data),"has_cnv_")] > 0 ))
Swarbrick_CID4066_HER2_panCK_cnv_Neg_Chromosomes <- as.data.frame(rowSums(Swarbrick_CID4066_HER2_panCK_cnv@meta.data[,startsWith(names(Swarbrick_CID4066_HER2_panCK_cnv@meta.data),"has_cnv_")] == 0))
Swarbrick_CID4066_HER2_panCK_cnv_Chromosome_Quant <- as.data.frame(c(Swarbrick_CID4066_HER2_panCK_cnv_Pos_Chromosomes, Swarbrick_CID4066_HER2_panCK_cnv_Neg_Chromosomes))
colnames(Swarbrick_CID4066_HER2_panCK_cnv_Chromosome_Quant) <- c("Chromosomes CNV Pos", "Chromosomes CNV Neg")
write.csv(Swarbrick_CID4066_HER2_panCK_cnv_Chromosome_Quant, file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/Swarbrick_CID4066_HER2_panCK/Swarbrick_CID4066_HER2_panCK_cnv_Chromosome_Quant.csv")

# Save CNV Metadata
Swarbrick_CID4066_HER2_panCK_cnv.meta.data <- as.data.frame(as.matrix(Swarbrick_CID4066_HER2_panCK_cnv@meta.data))
write.csv(Swarbrick_CID4066_HER2_panCK_cnv.meta.data, file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/Swarbrick_CID4066_HER2_panCK/Swarbrick_CID4066_HER2_panCK_cnv.meta.data.csv")

# Subset Cells with inferred CNVs -------------------------------------------------------------------------------------------------------------
Swarbrick_CID4066_HER2_panCK_cnv <- subset(Swarbrick_CID4066_HER2_panCK_cnv, cells=rownames(Swarbrick_CID4066_HER2_panCK_cnv@meta.data)[rowSums(Swarbrick_CID4066_HER2_panCK_cnv@meta.data[,startsWith(names(Swarbrick_CID4066_HER2_panCK_cnv@meta.data),"has_cnv_")]) > 0 ])

# Save RDS with inferredCNVs
saveRDS(Swarbrick_CID4066_HER2_panCK_cnv, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_inferCNV/Swarbrick_CID4066_HER2_panCK_cnv.rds")


# Free memory
rm(Swarbrick_CID4066_HER2_inferCNV)
rm(Swarbrick_CID4066_HER2_Navin_Visvader_NORM_panCK_Reference)
rm(Swarbrick_CID4066_HER2_panCK_cnv)
rm(Swarbrick_CID4066_HER2_panCK_cnv_Pos_Chromosomes)
rm(Swarbrick_CID4066_HER2_panCK_cnv_Neg_Chromosomes)
rm(Swarbrick_CID4066_HER2_panCK_cnv_Chromosome_Quant)
rm(Swarbrick_CID4066_HER2_panCK_cnv.meta.data)
gc()


#################################################
#### InferCNV: Swarbrick_CID45171_HER2_panCK ####
#################################################

########## Merge with normal mammary glad, extract epithelium, and process via Seurat
# Load RDS
Swarbrick_CID45171_HER2_panCK <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_panCK/Swarbrick_CID45171_HER2_panCK.rds")
Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Navin_Visvader_NORM_panCK_Reference.rds")

# Set identities 
colnames(x = Swarbrick_CID45171_HER2_panCK[[]])
Idents(object = Swarbrick_CID45171_HER2_panCK) <- 'Tumor'
Swarbrick_CID45171_HER2_panCK[["cnv.ident"]] <- Idents(object = Swarbrick_CID45171_HER2_panCK)
levels(x = Swarbrick_CID45171_HER2_panCK)

colnames(x = Navin_Visvader_NORM_panCK_Reference[[]])
Idents(object = Navin_Visvader_NORM_panCK_Reference) <- 'Normal'
Navin_Visvader_NORM_panCK_Reference[["cnv.ident"]] <- Idents(object = Navin_Visvader_NORM_panCK_Reference)
levels(x = Navin_Visvader_NORM_panCK_Reference)

# Merge subsets into recombined RDS
Swarbrick_CID45171_HER2_Navin_Visvader_NORM_panCK_Reference <- merge(x = Swarbrick_CID45171_HER2_panCK, y = Navin_Visvader_NORM_panCK_Reference)
# Join Layers
Swarbrick_CID45171_HER2_Navin_Visvader_NORM_panCK_Reference <- JoinLayers(Swarbrick_CID45171_HER2_Navin_Visvader_NORM_panCK_Reference)

# Remove subsetted objects
rm(Swarbrick_CID45171_HER2_panCK)
rm(Navin_Visvader_NORM_panCK_Reference)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples before merge

# Find Variable Features
Swarbrick_CID45171_HER2_Navin_Visvader_NORM_panCK_Reference <- FindVariableFeatures(Swarbrick_CID45171_HER2_Navin_Visvader_NORM_panCK_Reference, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Swarbrick_CID45171_HER2_Navin_Visvader_NORM_panCK_Reference)
Swarbrick_CID45171_HER2_Navin_Visvader_NORM_panCK_Reference <- ScaleData(Swarbrick_CID45171_HER2_Navin_Visvader_NORM_panCK_Reference, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################

Swarbrick_CID45171_HER2_Navin_Visvader_NORM_panCK_Reference <- RunPCA(Swarbrick_CID45171_HER2_Navin_Visvader_NORM_panCK_Reference, verbose = FALSE)

# Examine and visualize PCA results a few different ways
print(Swarbrick_CID45171_HER2_Navin_Visvader_NORM_panCK_Reference[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Swarbrick_CID45171_HER2_Navin_Visvader_NORM_panCK_Reference, dims = 1:2, reduction = "pca")
DimPlot(Swarbrick_CID45171_HER2_Navin_Visvader_NORM_panCK_Reference, reduction = "pca")

Swarbrick_CID45171_HER2_Navin_Visvader_NORM_panCK_Reference <- FindNeighbors(Swarbrick_CID45171_HER2_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
Swarbrick_CID45171_HER2_Navin_Visvader_NORM_panCK_Reference <- FindClusters(Swarbrick_CID45171_HER2_Navin_Visvader_NORM_panCK_Reference, resolution = 0.1)

# Umap clustering
Swarbrick_CID45171_HER2_Navin_Visvader_NORM_panCK_Reference <- RunUMAP(Swarbrick_CID45171_HER2_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
DimPlot(Swarbrick_CID45171_HER2_Navin_Visvader_NORM_panCK_Reference, reduction = "umap")
DimPlot(Swarbrick_CID45171_HER2_Navin_Visvader_NORM_panCK_Reference, reduction = "umap", group.by = "orig.ident")

# Save RDS
saveRDS(Swarbrick_CID45171_HER2_Navin_Visvader_NORM_panCK_Reference, file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_inferCNV_Input_RDS/Swarbrick_CID45171_HER2_Navin_Visvader_NORM_panCK_Reference.rds")

########## Prepare Cell Annotations
# Extract & Save Meta Data
Swarbrick_CID45171_HER2_Navin_Visvader_NORM_panCK_Reference.meta.data <- as.data.frame(as.matrix(Swarbrick_CID45171_HER2_Navin_Visvader_NORM_panCK_Reference@meta.data))
# Extract pertinent columns
# Make row names first column
Swarbrick_CID45171_HER2_Navin_Visvader_NORM_panCK_Reference.meta.data <- tibble::rownames_to_column(Swarbrick_CID45171_HER2_Navin_Visvader_NORM_panCK_Reference.meta.data, "Cells")
# In Seurat, "X" becomes the cell names while "cnv.ident" is specified above
Swarbrick_CID45171_HER2_Navin_Visvader_NORM_panCK_Reference.annotations <- Swarbrick_CID45171_HER2_Navin_Visvader_NORM_panCK_Reference.meta.data[, c("Cells", "cnv.ident")]
# Delete header
names(Swarbrick_CID45171_HER2_Navin_Visvader_NORM_panCK_Reference.annotations) <- NULL
# Save as table
write.table(Swarbrick_CID45171_HER2_Navin_Visvader_NORM_panCK_Reference.annotations,"/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_CID45171_HER2_Navin_Visvader_NORM_panCK_Reference.annotations.txt",sep="\t",row.names=FALSE)
# Remove metadata and annotations
rm(Swarbrick_CID45171_HER2_Navin_Visvader_NORM_panCK_Reference.meta.data)
rm(Swarbrick_CID45171_HER2_Navin_Visvader_NORM_panCK_Reference.annotations)

########## Generate Count Matrix
# Extract Count Matrix
Swarbrick_CID45171_HER2_Navin_Visvader_NORM_panCK_Reference.count <- Swarbrick_CID45171_HER2_Navin_Visvader_NORM_panCK_Reference[["RNA"]]$counts
Swarbrick_CID45171_HER2_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as(Swarbrick_CID45171_HER2_Navin_Visvader_NORM_panCK_Reference.count, 'sparseMatrix')
rm(Swarbrick_CID45171_HER2_Navin_Visvader_NORM_panCK_Reference.count)
# Convert to requisite data format
Swarbrick_CID45171_HER2_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as.data.frame(Swarbrick_CID45171_HER2_Navin_Visvader_NORM_panCK_Reference.count.matrix)
# Save table
write.csv(Swarbrick_CID45171_HER2_Navin_Visvader_NORM_panCK_Reference.count.matrix, file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_inferCNV_Input_Count_Matrix/Swarbrick_CID45171_HER2_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv")
# Remove Objects
rm(Swarbrick_CID45171_HER2_Navin_Visvader_NORM_panCK_Reference)
rm(Swarbrick_CID45171_HER2_Navin_Visvader_NORM_panCK_Reference.count.matrix)

########### Run inferCNV
# Note: check.names = FALSE prevents cell names from having dashes replaced with dots; essential to match annotations file
Swarbrick_CID45171_HER2_inferCNV <- CreateInfercnvObject(raw_counts_matrix=read.csv(file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_inferCNV_Input_Count_Matrix/Swarbrick_CID45171_HER2_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv", header = TRUE, row.names = 1,check.names = FALSE),
                                               annotations_file=read.table(file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_CID45171_HER2_Navin_Visvader_NORM_panCK_Reference.annotations.txt", header = FALSE, row.names = 1),
                                               delim="\t",
                                               gene_order_file="/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/hg38_gencode_v27.txt",
                                               ref_group_names=c("Normal"))


# perform infercnv operations to reveal cnv signal
Swarbrick_CID45171_HER2_inferCNV = infercnv::run(Swarbrick_CID45171_HER2_inferCNV,
                                       cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                       out_dir="/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/Swarbrick_CID45171_HER2_panCK/",  # dir is auto-created for storing outputs
                                       cluster_by_groups=T,   # cluster
                                       denoise=T,
                                       HMM=T,
                                       num_ref_groups = 9, analysis_mode = "subclusters")

###########  Integrate data with Seurat object
# Read RDS and add inferCNV results to metadata
Swarbrick_CID45171_HER2_Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_inferCNV_Input_RDS/Swarbrick_CID45171_HER2_Navin_Visvader_NORM_panCK_Reference.rds") 
Swarbrick_CID45171_HER2_panCK_cnv <- infercnv::add_to_seurat(infercnv_output_path="/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/Swarbrick_CID45171_HER2_panCK/",
                                                   seurat_obj=Swarbrick_CID45171_HER2_Navin_Visvader_NORM_panCK_Reference, # optional
                                                   top_n=10)

# Collect Tumor cells from Seurat Object
Swarbrick_CID45171_HER2_panCK_cnv <- subset(x = Swarbrick_CID45171_HER2_panCK_cnv, subset = cnv.ident == "Tumor")
# Quantify chromosomes with alterations per cell
Swarbrick_CID45171_HER2_panCK_cnv_Pos_Chromosomes <- as.data.frame(rowSums(Swarbrick_CID45171_HER2_panCK_cnv@meta.data[,startsWith(names(Swarbrick_CID45171_HER2_panCK_cnv@meta.data),"has_cnv_")] > 0 ))
Swarbrick_CID45171_HER2_panCK_cnv_Neg_Chromosomes <- as.data.frame(rowSums(Swarbrick_CID45171_HER2_panCK_cnv@meta.data[,startsWith(names(Swarbrick_CID45171_HER2_panCK_cnv@meta.data),"has_cnv_")] == 0))
Swarbrick_CID45171_HER2_panCK_cnv_Chromosome_Quant <- as.data.frame(c(Swarbrick_CID45171_HER2_panCK_cnv_Pos_Chromosomes, Swarbrick_CID45171_HER2_panCK_cnv_Neg_Chromosomes))
colnames(Swarbrick_CID45171_HER2_panCK_cnv_Chromosome_Quant) <- c("Chromosomes CNV Pos", "Chromosomes CNV Neg")
write.csv(Swarbrick_CID45171_HER2_panCK_cnv_Chromosome_Quant, file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/Swarbrick_CID45171_HER2_panCK/Swarbrick_CID45171_HER2_panCK_cnv_Chromosome_Quant.csv")

# Save CNV Metadata
Swarbrick_CID45171_HER2_panCK_cnv.meta.data <- as.data.frame(as.matrix(Swarbrick_CID45171_HER2_panCK_cnv@meta.data))
write.csv(Swarbrick_CID45171_HER2_panCK_cnv.meta.data, file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/Swarbrick_CID45171_HER2_panCK/Swarbrick_CID45171_HER2_panCK_cnv.meta.data.csv")

# Subset Cells with inferred CNVs -------------------------------------------------------------------------------------------------------------
Swarbrick_CID45171_HER2_panCK_cnv <- subset(Swarbrick_CID45171_HER2_panCK_cnv, cells=rownames(Swarbrick_CID45171_HER2_panCK_cnv@meta.data)[rowSums(Swarbrick_CID45171_HER2_panCK_cnv@meta.data[,startsWith(names(Swarbrick_CID45171_HER2_panCK_cnv@meta.data),"has_cnv_")]) > 0 ])

# Save RDS with inferredCNVs
saveRDS(Swarbrick_CID45171_HER2_panCK_cnv, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_inferCNV/Swarbrick_CID45171_HER2_panCK_cnv.rds")


# Free memory
rm(Swarbrick_CID45171_HER2_inferCNV)
rm(Swarbrick_CID45171_HER2_Navin_Visvader_NORM_panCK_Reference)
rm(Swarbrick_CID45171_HER2_panCK_cnv)
rm(Swarbrick_CID45171_HER2_panCK_cnv_Pos_Chromosomes)
rm(Swarbrick_CID45171_HER2_panCK_cnv_Neg_Chromosomes)
rm(Swarbrick_CID45171_HER2_panCK_cnv_Chromosome_Quant)
rm(Swarbrick_CID45171_HER2_panCK_cnv.meta.data)
gc()





###################################################################################################
#############################     InferCNV: TNBC Tumor Specimens      #############################     
###################################################################################################

################################################
#### InferCNV: Swarbrick_CID3946_TNBC_panCK ####
################################################

########## Merge with normal mammary glad, extract epithelium, and process via Seurat
# Load RDS
Swarbrick_CID3946_TNBC_panCK <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_panCK/Swarbrick_CID3946_TNBC_panCK.rds")
Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Navin_Visvader_NORM_panCK_Reference.rds")

# Set identities 
colnames(x = Swarbrick_CID3946_TNBC_panCK[[]])
Idents(object = Swarbrick_CID3946_TNBC_panCK) <- 'Tumor'
Swarbrick_CID3946_TNBC_panCK[["cnv.ident"]] <- Idents(object = Swarbrick_CID3946_TNBC_panCK)
levels(x = Swarbrick_CID3946_TNBC_panCK)

colnames(x = Navin_Visvader_NORM_panCK_Reference[[]])
Idents(object = Navin_Visvader_NORM_panCK_Reference) <- 'Normal'
Navin_Visvader_NORM_panCK_Reference[["cnv.ident"]] <- Idents(object = Navin_Visvader_NORM_panCK_Reference)
levels(x = Navin_Visvader_NORM_panCK_Reference)

# Merge subsets into recombined RDS
Swarbrick_CID3946_TNBC_Navin_Visvader_NORM_panCK_Reference <- merge(x = Swarbrick_CID3946_TNBC_panCK, y = Navin_Visvader_NORM_panCK_Reference)
# Join Layers
Swarbrick_CID3946_TNBC_Navin_Visvader_NORM_panCK_Reference <- JoinLayers(Swarbrick_CID3946_TNBC_Navin_Visvader_NORM_panCK_Reference)

# Remove subsetted objects
rm(Swarbrick_CID3946_TNBC_panCK)
rm(Navin_Visvader_NORM_panCK_Reference)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples before merge

# Find Variable Features
Swarbrick_CID3946_TNBC_Navin_Visvader_NORM_panCK_Reference <- FindVariableFeatures(Swarbrick_CID3946_TNBC_Navin_Visvader_NORM_panCK_Reference, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Swarbrick_CID3946_TNBC_Navin_Visvader_NORM_panCK_Reference)
Swarbrick_CID3946_TNBC_Navin_Visvader_NORM_panCK_Reference <- ScaleData(Swarbrick_CID3946_TNBC_Navin_Visvader_NORM_panCK_Reference, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################

Swarbrick_CID3946_TNBC_Navin_Visvader_NORM_panCK_Reference <- RunPCA(Swarbrick_CID3946_TNBC_Navin_Visvader_NORM_panCK_Reference, verbose = FALSE)

# Examine and visualize PCA results a few different ways
print(Swarbrick_CID3946_TNBC_Navin_Visvader_NORM_panCK_Reference[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Swarbrick_CID3946_TNBC_Navin_Visvader_NORM_panCK_Reference, dims = 1:2, reduction = "pca")
DimPlot(Swarbrick_CID3946_TNBC_Navin_Visvader_NORM_panCK_Reference, reduction = "pca")

Swarbrick_CID3946_TNBC_Navin_Visvader_NORM_panCK_Reference <- FindNeighbors(Swarbrick_CID3946_TNBC_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
Swarbrick_CID3946_TNBC_Navin_Visvader_NORM_panCK_Reference <- FindClusters(Swarbrick_CID3946_TNBC_Navin_Visvader_NORM_panCK_Reference, resolution = 0.1)

# Umap clustering
Swarbrick_CID3946_TNBC_Navin_Visvader_NORM_panCK_Reference <- RunUMAP(Swarbrick_CID3946_TNBC_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
DimPlot(Swarbrick_CID3946_TNBC_Navin_Visvader_NORM_panCK_Reference, reduction = "umap")
DimPlot(Swarbrick_CID3946_TNBC_Navin_Visvader_NORM_panCK_Reference, reduction = "umap", group.by = "orig.ident")

# Save RDS
saveRDS(Swarbrick_CID3946_TNBC_Navin_Visvader_NORM_panCK_Reference, file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_inferCNV_Input_RDS/Swarbrick_CID3946_TNBC_Navin_Visvader_NORM_panCK_Reference.rds")

########## Prepare Cell Annotations
# Extract & Save Meta Data
Swarbrick_CID3946_TNBC_Navin_Visvader_NORM_panCK_Reference.meta.data <- as.data.frame(as.matrix(Swarbrick_CID3946_TNBC_Navin_Visvader_NORM_panCK_Reference@meta.data))
# Extract pertinent columns
# Make row names first column
Swarbrick_CID3946_TNBC_Navin_Visvader_NORM_panCK_Reference.meta.data <- tibble::rownames_to_column(Swarbrick_CID3946_TNBC_Navin_Visvader_NORM_panCK_Reference.meta.data, "Cells")
# In Seurat, "X" becomes the cell names while "cnv.ident" is specified above
Swarbrick_CID3946_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations <- Swarbrick_CID3946_TNBC_Navin_Visvader_NORM_panCK_Reference.meta.data[, c("Cells", "cnv.ident")]
# Delete header
names(Swarbrick_CID3946_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations) <- NULL
# Save as table
write.table(Swarbrick_CID3946_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations,"/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_CID3946_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations.txt",sep="\t",row.names=FALSE)
# Remove metadata and annotations
rm(Swarbrick_CID3946_TNBC_Navin_Visvader_NORM_panCK_Reference.meta.data)
rm(Swarbrick_CID3946_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations)

########## Generate Count Matrix
# Extract Count Matrix
Swarbrick_CID3946_TNBC_Navin_Visvader_NORM_panCK_Reference.count <- Swarbrick_CID3946_TNBC_Navin_Visvader_NORM_panCK_Reference[["RNA"]]$counts
Swarbrick_CID3946_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as(Swarbrick_CID3946_TNBC_Navin_Visvader_NORM_panCK_Reference.count, 'sparseMatrix')
rm(Swarbrick_CID3946_TNBC_Navin_Visvader_NORM_panCK_Reference.count)
# Convert to requisite data format
Swarbrick_CID3946_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as.data.frame(Swarbrick_CID3946_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix)
# Save table
write.csv(Swarbrick_CID3946_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix, file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_inferCNV_Input_Count_Matrix/Swarbrick_CID3946_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv")
# Remove Objects
rm(Swarbrick_CID3946_TNBC_Navin_Visvader_NORM_panCK_Reference)
rm(Swarbrick_CID3946_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix)

########### Run inferCNV
# Note: check.names = FALSE prevents cell names from having dashes replaced with dots; essential to match annotations file
Swarbrick_CID3946_TNBC_inferCNV <- CreateInfercnvObject(raw_counts_matrix=read.csv(file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_inferCNV_Input_Count_Matrix/Swarbrick_CID3946_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv", header = TRUE, row.names = 1,check.names = FALSE),
                                              annotations_file=read.table(file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_CID3946_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations.txt", header = FALSE, row.names = 1),
                                              delim="\t",
                                              gene_order_file="/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/hg38_gencode_v27.txt",
                                              ref_group_names=c("Normal"))


# perform infercnv operations to reveal cnv signal
Swarbrick_CID3946_TNBC_inferCNV = infercnv::run(Swarbrick_CID3946_TNBC_inferCNV,
                                      cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                      out_dir="/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/Swarbrick_CID3946_TNBC_panCK/",  # dir is auto-created for storing outputs
                                      cluster_by_groups=T,   # cluster
                                      denoise=T,
                                      HMM=T,
                                      num_ref_groups = 9, analysis_mode = "subclusters")

###########  Integrate data with Seurat object
# Read RDS and add inferCNV results to metadata
Swarbrick_CID3946_TNBC_Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_inferCNV_Input_RDS/Swarbrick_CID3946_TNBC_Navin_Visvader_NORM_panCK_Reference.rds") 
Swarbrick_CID3946_TNBC_panCK_cnv <- infercnv::add_to_seurat(infercnv_output_path="/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/Swarbrick_CID3946_TNBC_panCK/",
                                                  seurat_obj=Swarbrick_CID3946_TNBC_Navin_Visvader_NORM_panCK_Reference, # optional
                                                  top_n=10)

# Collect Tumor cells from Seurat Object
Swarbrick_CID3946_TNBC_panCK_cnv <- subset(x = Swarbrick_CID3946_TNBC_panCK_cnv, subset = cnv.ident == "Tumor")
# Quantify chromosomes with alterations per cell
Swarbrick_CID3946_TNBC_panCK_cnv_Pos_Chromosomes <- as.data.frame(rowSums(Swarbrick_CID3946_TNBC_panCK_cnv@meta.data[,startsWith(names(Swarbrick_CID3946_TNBC_panCK_cnv@meta.data),"has_cnv_")] > 0 ))
Swarbrick_CID3946_TNBC_panCK_cnv_Neg_Chromosomes <- as.data.frame(rowSums(Swarbrick_CID3946_TNBC_panCK_cnv@meta.data[,startsWith(names(Swarbrick_CID3946_TNBC_panCK_cnv@meta.data),"has_cnv_")] == 0))
Swarbrick_CID3946_TNBC_panCK_cnv_Chromosome_Quant <- as.data.frame(c(Swarbrick_CID3946_TNBC_panCK_cnv_Pos_Chromosomes, Swarbrick_CID3946_TNBC_panCK_cnv_Neg_Chromosomes))
colnames(Swarbrick_CID3946_TNBC_panCK_cnv_Chromosome_Quant) <- c("Chromosomes CNV Pos", "Chromosomes CNV Neg")
write.csv(Swarbrick_CID3946_TNBC_panCK_cnv_Chromosome_Quant, file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/Swarbrick_CID3946_TNBC_panCK/Swarbrick_CID3946_TNBC_panCK_cnv_Chromosome_Quant.csv")

# Save CNV Metadata
Swarbrick_CID3946_TNBC_panCK_cnv.meta.data <- as.data.frame(as.matrix(Swarbrick_CID3946_TNBC_panCK_cnv@meta.data))
write.csv(Swarbrick_CID3946_TNBC_panCK_cnv.meta.data, file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/Swarbrick_CID3946_TNBC_panCK/Swarbrick_CID3946_TNBC_panCK_cnv.meta.data.csv")

# Subset Cells with inferred CNVs -------------------------------------------------------------------------------------------------------------
Swarbrick_CID3946_TNBC_panCK_cnv <- subset(Swarbrick_CID3946_TNBC_panCK_cnv, cells=rownames(Swarbrick_CID3946_TNBC_panCK_cnv@meta.data)[rowSums(Swarbrick_CID3946_TNBC_panCK_cnv@meta.data[,startsWith(names(Swarbrick_CID3946_TNBC_panCK_cnv@meta.data),"has_cnv_")]) > 0 ])

# Save RDS with inferredCNVs
saveRDS(Swarbrick_CID3946_TNBC_panCK_cnv, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_inferCNV/Swarbrick_CID3946_TNBC_panCK_cnv.rds")


# Free memory
rm(Swarbrick_CID3946_TNBC_inferCNV)
rm(Swarbrick_CID3946_TNBC_Navin_Visvader_NORM_panCK_Reference)
rm(Swarbrick_CID3946_TNBC_panCK_cnv)
rm(Swarbrick_CID3946_TNBC_panCK_cnv_Pos_Chromosomes)
rm(Swarbrick_CID3946_TNBC_panCK_cnv_Neg_Chromosomes)
rm(Swarbrick_CID3946_TNBC_panCK_cnv_Chromosome_Quant)
rm(Swarbrick_CID3946_TNBC_panCK_cnv.meta.data)
gc()



################################################
#### InferCNV: Swarbrick_CID3963_TNBC_panCK ####
################################################

########## Merge with normal mammary glad, extract epithelium, and process via Seurat
# Load RDS
Swarbrick_CID3963_TNBC_panCK <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_panCK/Swarbrick_CID3963_TNBC_panCK.rds")
Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Navin_Visvader_NORM_panCK_Reference.rds")

# Set identities 
colnames(x = Swarbrick_CID3963_TNBC_panCK[[]])
Idents(object = Swarbrick_CID3963_TNBC_panCK) <- 'Tumor'
Swarbrick_CID3963_TNBC_panCK[["cnv.ident"]] <- Idents(object = Swarbrick_CID3963_TNBC_panCK)
levels(x = Swarbrick_CID3963_TNBC_panCK)

colnames(x = Navin_Visvader_NORM_panCK_Reference[[]])
Idents(object = Navin_Visvader_NORM_panCK_Reference) <- 'Normal'
Navin_Visvader_NORM_panCK_Reference[["cnv.ident"]] <- Idents(object = Navin_Visvader_NORM_panCK_Reference)
levels(x = Navin_Visvader_NORM_panCK_Reference)

# Merge subsets into recombined RDS
Swarbrick_CID3963_TNBC_Navin_Visvader_NORM_panCK_Reference <- merge(x = Swarbrick_CID3963_TNBC_panCK, y = Navin_Visvader_NORM_panCK_Reference)
# Join Layers
Swarbrick_CID3963_TNBC_Navin_Visvader_NORM_panCK_Reference <- JoinLayers(Swarbrick_CID3963_TNBC_Navin_Visvader_NORM_panCK_Reference)

# Remove subsetted objects
rm(Swarbrick_CID3963_TNBC_panCK)
rm(Navin_Visvader_NORM_panCK_Reference)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples before merge

# Find Variable Features
Swarbrick_CID3963_TNBC_Navin_Visvader_NORM_panCK_Reference <- FindVariableFeatures(Swarbrick_CID3963_TNBC_Navin_Visvader_NORM_panCK_Reference, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Swarbrick_CID3963_TNBC_Navin_Visvader_NORM_panCK_Reference)
Swarbrick_CID3963_TNBC_Navin_Visvader_NORM_panCK_Reference <- ScaleData(Swarbrick_CID3963_TNBC_Navin_Visvader_NORM_panCK_Reference, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################

Swarbrick_CID3963_TNBC_Navin_Visvader_NORM_panCK_Reference <- RunPCA(Swarbrick_CID3963_TNBC_Navin_Visvader_NORM_panCK_Reference, verbose = FALSE)

# Examine and visualize PCA results a few different ways
print(Swarbrick_CID3963_TNBC_Navin_Visvader_NORM_panCK_Reference[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Swarbrick_CID3963_TNBC_Navin_Visvader_NORM_panCK_Reference, dims = 1:2, reduction = "pca")
DimPlot(Swarbrick_CID3963_TNBC_Navin_Visvader_NORM_panCK_Reference, reduction = "pca")

Swarbrick_CID3963_TNBC_Navin_Visvader_NORM_panCK_Reference <- FindNeighbors(Swarbrick_CID3963_TNBC_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
Swarbrick_CID3963_TNBC_Navin_Visvader_NORM_panCK_Reference <- FindClusters(Swarbrick_CID3963_TNBC_Navin_Visvader_NORM_panCK_Reference, resolution = 0.1)

# Umap clustering
Swarbrick_CID3963_TNBC_Navin_Visvader_NORM_panCK_Reference <- RunUMAP(Swarbrick_CID3963_TNBC_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
DimPlot(Swarbrick_CID3963_TNBC_Navin_Visvader_NORM_panCK_Reference, reduction = "umap")
DimPlot(Swarbrick_CID3963_TNBC_Navin_Visvader_NORM_panCK_Reference, reduction = "umap", group.by = "orig.ident")

# Save RDS
saveRDS(Swarbrick_CID3963_TNBC_Navin_Visvader_NORM_panCK_Reference, file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_inferCNV_Input_RDS/Swarbrick_CID3963_TNBC_Navin_Visvader_NORM_panCK_Reference.rds")

########## Prepare Cell Annotations
# Extract & Save Meta Data
Swarbrick_CID3963_TNBC_Navin_Visvader_NORM_panCK_Reference.meta.data <- as.data.frame(as.matrix(Swarbrick_CID3963_TNBC_Navin_Visvader_NORM_panCK_Reference@meta.data))
# Extract pertinent columns
# Make row names first column
Swarbrick_CID3963_TNBC_Navin_Visvader_NORM_panCK_Reference.meta.data <- tibble::rownames_to_column(Swarbrick_CID3963_TNBC_Navin_Visvader_NORM_panCK_Reference.meta.data, "Cells")
# In Seurat, "X" becomes the cell names while "cnv.ident" is specified above
Swarbrick_CID3963_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations <- Swarbrick_CID3963_TNBC_Navin_Visvader_NORM_panCK_Reference.meta.data[, c("Cells", "cnv.ident")]
# Delete header
names(Swarbrick_CID3963_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations) <- NULL
# Save as table
write.table(Swarbrick_CID3963_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations,"/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_CID3963_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations.txt",sep="\t",row.names=FALSE)
# Remove metadata and annotations
rm(Swarbrick_CID3963_TNBC_Navin_Visvader_NORM_panCK_Reference.meta.data)
rm(Swarbrick_CID3963_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations)

########## Generate Count Matrix
# Extract Count Matrix
Swarbrick_CID3963_TNBC_Navin_Visvader_NORM_panCK_Reference.count <- Swarbrick_CID3963_TNBC_Navin_Visvader_NORM_panCK_Reference[["RNA"]]$counts
Swarbrick_CID3963_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as(Swarbrick_CID3963_TNBC_Navin_Visvader_NORM_panCK_Reference.count, 'sparseMatrix')
rm(Swarbrick_CID3963_TNBC_Navin_Visvader_NORM_panCK_Reference.count)
# Convert to requisite data format
Swarbrick_CID3963_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as.data.frame(Swarbrick_CID3963_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix)
# Save table
write.csv(Swarbrick_CID3963_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix, file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_inferCNV_Input_Count_Matrix/Swarbrick_CID3963_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv")
# Remove Objects
rm(Swarbrick_CID3963_TNBC_Navin_Visvader_NORM_panCK_Reference)
rm(Swarbrick_CID3963_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix)

########### Run inferCNV
# Note: check.names = FALSE prevents cell names from having dashes replaced with dots; essential to match annotations file
Swarbrick_CID3963_TNBC_inferCNV <- CreateInfercnvObject(raw_counts_matrix=read.csv(file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_inferCNV_Input_Count_Matrix/Swarbrick_CID3963_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv", header = TRUE, row.names = 1,check.names = FALSE),
                                              annotations_file=read.table(file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_CID3963_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations.txt", header = FALSE, row.names = 1),
                                              delim="\t",
                                              gene_order_file="/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/hg38_gencode_v27.txt",
                                              ref_group_names=c("Normal"))


# perform infercnv operations to reveal cnv signal
Swarbrick_CID3963_TNBC_inferCNV = infercnv::run(Swarbrick_CID3963_TNBC_inferCNV,
                                      cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                      out_dir="/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/Swarbrick_CID3963_TNBC_panCK/",  # dir is auto-created for storing outputs
                                      cluster_by_groups=T,   # cluster
                                      denoise=T,
                                      HMM=T,
                                      num_ref_groups = 9, analysis_mode = "subclusters")

###########  Integrate data with Seurat object
# Read RDS and add inferCNV results to metadata
Swarbrick_CID3963_TNBC_Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_inferCNV_Input_RDS/Swarbrick_CID3963_TNBC_Navin_Visvader_NORM_panCK_Reference.rds") 
Swarbrick_CID3963_TNBC_panCK_cnv <- infercnv::add_to_seurat(infercnv_output_path="/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/Swarbrick_CID3963_TNBC_panCK/",
                                                  seurat_obj=Swarbrick_CID3963_TNBC_Navin_Visvader_NORM_panCK_Reference, # optional
                                                  top_n=10)

# Collect Tumor cells from Seurat Object
Swarbrick_CID3963_TNBC_panCK_cnv <- subset(x = Swarbrick_CID3963_TNBC_panCK_cnv, subset = cnv.ident == "Tumor")
# Quantify chromosomes with alterations per cell
Swarbrick_CID3963_TNBC_panCK_cnv_Pos_Chromosomes <- as.data.frame(rowSums(Swarbrick_CID3963_TNBC_panCK_cnv@meta.data[,startsWith(names(Swarbrick_CID3963_TNBC_panCK_cnv@meta.data),"has_cnv_")] > 0 ))
Swarbrick_CID3963_TNBC_panCK_cnv_Neg_Chromosomes <- as.data.frame(rowSums(Swarbrick_CID3963_TNBC_panCK_cnv@meta.data[,startsWith(names(Swarbrick_CID3963_TNBC_panCK_cnv@meta.data),"has_cnv_")] == 0))
Swarbrick_CID3963_TNBC_panCK_cnv_Chromosome_Quant <- as.data.frame(c(Swarbrick_CID3963_TNBC_panCK_cnv_Pos_Chromosomes, Swarbrick_CID3963_TNBC_panCK_cnv_Neg_Chromosomes))
colnames(Swarbrick_CID3963_TNBC_panCK_cnv_Chromosome_Quant) <- c("Chromosomes CNV Pos", "Chromosomes CNV Neg")
write.csv(Swarbrick_CID3963_TNBC_panCK_cnv_Chromosome_Quant, file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/Swarbrick_CID3963_TNBC_panCK/Swarbrick_CID3963_TNBC_panCK_cnv_Chromosome_Quant.csv")

# Save CNV Metadata
Swarbrick_CID3963_TNBC_panCK_cnv.meta.data <- as.data.frame(as.matrix(Swarbrick_CID3963_TNBC_panCK_cnv@meta.data))
write.csv(Swarbrick_CID3963_TNBC_panCK_cnv.meta.data, file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/Swarbrick_CID3963_TNBC_panCK/Swarbrick_CID3963_TNBC_panCK_cnv.meta.data.csv")

# Subset Cells with inferred CNVs -------------------------------------------------------------------------------------------------------------
Swarbrick_CID3963_TNBC_panCK_cnv <- subset(Swarbrick_CID3963_TNBC_panCK_cnv, cells=rownames(Swarbrick_CID3963_TNBC_panCK_cnv@meta.data)[rowSums(Swarbrick_CID3963_TNBC_panCK_cnv@meta.data[,startsWith(names(Swarbrick_CID3963_TNBC_panCK_cnv@meta.data),"has_cnv_")]) > 0 ])

# Save RDS with inferredCNVs
saveRDS(Swarbrick_CID3963_TNBC_panCK_cnv, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_inferCNV/Swarbrick_CID3963_TNBC_panCK_cnv.rds")


# Free memory
rm(Swarbrick_CID3963_TNBC_inferCNV)
rm(Swarbrick_CID3963_TNBC_Navin_Visvader_NORM_panCK_Reference)
rm(Swarbrick_CID3963_TNBC_panCK_cnv)
rm(Swarbrick_CID3963_TNBC_panCK_cnv_Pos_Chromosomes)
rm(Swarbrick_CID3963_TNBC_panCK_cnv_Neg_Chromosomes)
rm(Swarbrick_CID3963_TNBC_panCK_cnv_Chromosome_Quant)
rm(Swarbrick_CID3963_TNBC_panCK_cnv.meta.data)
gc()



################################################
#### InferCNV: Swarbrick_CID4465_TNBC_panCK ####
################################################

########## Merge with normal mammary glad, extract epithelium, and process via Seurat
# Load RDS
Swarbrick_CID4465_TNBC_panCK <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_panCK/Swarbrick_CID4465_TNBC_panCK.rds")
Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Navin_Visvader_NORM_panCK_Reference.rds")

# Set identities 
colnames(x = Swarbrick_CID4465_TNBC_panCK[[]])
Idents(object = Swarbrick_CID4465_TNBC_panCK) <- 'Tumor'
Swarbrick_CID4465_TNBC_panCK[["cnv.ident"]] <- Idents(object = Swarbrick_CID4465_TNBC_panCK)
levels(x = Swarbrick_CID4465_TNBC_panCK)

colnames(x = Navin_Visvader_NORM_panCK_Reference[[]])
Idents(object = Navin_Visvader_NORM_panCK_Reference) <- 'Normal'
Navin_Visvader_NORM_panCK_Reference[["cnv.ident"]] <- Idents(object = Navin_Visvader_NORM_panCK_Reference)
levels(x = Navin_Visvader_NORM_panCK_Reference)

# Merge subsets into recombined RDS
Swarbrick_CID4465_TNBC_Navin_Visvader_NORM_panCK_Reference <- merge(x = Swarbrick_CID4465_TNBC_panCK, y = Navin_Visvader_NORM_panCK_Reference)
# Join Layers
Swarbrick_CID4465_TNBC_Navin_Visvader_NORM_panCK_Reference <- JoinLayers(Swarbrick_CID4465_TNBC_Navin_Visvader_NORM_panCK_Reference)

# Remove subsetted objects
rm(Swarbrick_CID4465_TNBC_panCK)
rm(Navin_Visvader_NORM_panCK_Reference)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples before merge

# Find Variable Features
Swarbrick_CID4465_TNBC_Navin_Visvader_NORM_panCK_Reference <- FindVariableFeatures(Swarbrick_CID4465_TNBC_Navin_Visvader_NORM_panCK_Reference, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Swarbrick_CID4465_TNBC_Navin_Visvader_NORM_panCK_Reference)
Swarbrick_CID4465_TNBC_Navin_Visvader_NORM_panCK_Reference <- ScaleData(Swarbrick_CID4465_TNBC_Navin_Visvader_NORM_panCK_Reference, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################

Swarbrick_CID4465_TNBC_Navin_Visvader_NORM_panCK_Reference <- RunPCA(Swarbrick_CID4465_TNBC_Navin_Visvader_NORM_panCK_Reference, verbose = FALSE)

# Examine and visualize PCA results a few different ways
print(Swarbrick_CID4465_TNBC_Navin_Visvader_NORM_panCK_Reference[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Swarbrick_CID4465_TNBC_Navin_Visvader_NORM_panCK_Reference, dims = 1:2, reduction = "pca")
DimPlot(Swarbrick_CID4465_TNBC_Navin_Visvader_NORM_panCK_Reference, reduction = "pca")

Swarbrick_CID4465_TNBC_Navin_Visvader_NORM_panCK_Reference <- FindNeighbors(Swarbrick_CID4465_TNBC_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
Swarbrick_CID4465_TNBC_Navin_Visvader_NORM_panCK_Reference <- FindClusters(Swarbrick_CID4465_TNBC_Navin_Visvader_NORM_panCK_Reference, resolution = 0.1)

# Umap clustering
Swarbrick_CID4465_TNBC_Navin_Visvader_NORM_panCK_Reference <- RunUMAP(Swarbrick_CID4465_TNBC_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
DimPlot(Swarbrick_CID4465_TNBC_Navin_Visvader_NORM_panCK_Reference, reduction = "umap")
DimPlot(Swarbrick_CID4465_TNBC_Navin_Visvader_NORM_panCK_Reference, reduction = "umap", group.by = "orig.ident")

# Save RDS
saveRDS(Swarbrick_CID4465_TNBC_Navin_Visvader_NORM_panCK_Reference, file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_inferCNV_Input_RDS/Swarbrick_CID4465_TNBC_Navin_Visvader_NORM_panCK_Reference.rds")

########## Prepare Cell Annotations
# Extract & Save Meta Data
Swarbrick_CID4465_TNBC_Navin_Visvader_NORM_panCK_Reference.meta.data <- as.data.frame(as.matrix(Swarbrick_CID4465_TNBC_Navin_Visvader_NORM_panCK_Reference@meta.data))
# Extract pertinent columns
# Make row names first column
Swarbrick_CID4465_TNBC_Navin_Visvader_NORM_panCK_Reference.meta.data <- tibble::rownames_to_column(Swarbrick_CID4465_TNBC_Navin_Visvader_NORM_panCK_Reference.meta.data, "Cells")
# In Seurat, "X" becomes the cell names while "cnv.ident" is specified above
Swarbrick_CID4465_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations <- Swarbrick_CID4465_TNBC_Navin_Visvader_NORM_panCK_Reference.meta.data[, c("Cells", "cnv.ident")]
# Delete header
names(Swarbrick_CID4465_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations) <- NULL
# Save as table
write.table(Swarbrick_CID4465_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations,"/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_CID4465_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations.txt",sep="\t",row.names=FALSE)
# Remove metadata and annotations
rm(Swarbrick_CID4465_TNBC_Navin_Visvader_NORM_panCK_Reference.meta.data)
rm(Swarbrick_CID4465_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations)

########## Generate Count Matrix
# Extract Count Matrix
Swarbrick_CID4465_TNBC_Navin_Visvader_NORM_panCK_Reference.count <- Swarbrick_CID4465_TNBC_Navin_Visvader_NORM_panCK_Reference[["RNA"]]$counts
Swarbrick_CID4465_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as(Swarbrick_CID4465_TNBC_Navin_Visvader_NORM_panCK_Reference.count, 'sparseMatrix')
rm(Swarbrick_CID4465_TNBC_Navin_Visvader_NORM_panCK_Reference.count)
# Convert to requisite data format
Swarbrick_CID4465_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as.data.frame(Swarbrick_CID4465_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix)
# Save table
write.csv(Swarbrick_CID4465_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix, file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_inferCNV_Input_Count_Matrix/Swarbrick_CID4465_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv")
# Remove Objects
rm(Swarbrick_CID4465_TNBC_Navin_Visvader_NORM_panCK_Reference)
rm(Swarbrick_CID4465_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix)

########### Run inferCNV
# Note: check.names = FALSE prevents cell names from having dashes replaced with dots; essential to match annotations file
Swarbrick_CID4465_TNBC_inferCNV <- CreateInfercnvObject(raw_counts_matrix=read.csv(file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_inferCNV_Input_Count_Matrix/Swarbrick_CID4465_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv", header = TRUE, row.names = 1,check.names = FALSE),
                                              annotations_file=read.table(file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_CID4465_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations.txt", header = FALSE, row.names = 1),
                                              delim="\t",
                                              gene_order_file="/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/hg38_gencode_v27.txt",
                                              ref_group_names=c("Normal"))


# perform infercnv operations to reveal cnv signal
Swarbrick_CID4465_TNBC_inferCNV = infercnv::run(Swarbrick_CID4465_TNBC_inferCNV,
                                      cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                      out_dir="/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/Swarbrick_CID4465_TNBC_panCK/",  # dir is auto-created for storing outputs
                                      cluster_by_groups=T,   # cluster
                                      denoise=T,
                                      HMM=T,
                                      num_ref_groups = 9, analysis_mode = "subclusters")

###########  Integrate data with Seurat object
# Read RDS and add inferCNV results to metadata
Swarbrick_CID4465_TNBC_Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_inferCNV_Input_RDS/Swarbrick_CID4465_TNBC_Navin_Visvader_NORM_panCK_Reference.rds") 
Swarbrick_CID4465_TNBC_panCK_cnv <- infercnv::add_to_seurat(infercnv_output_path="/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/Swarbrick_CID4465_TNBC_panCK/",
                                                  seurat_obj=Swarbrick_CID4465_TNBC_Navin_Visvader_NORM_panCK_Reference, # optional
                                                  top_n=10)

# Collect Tumor cells from Seurat Object
Swarbrick_CID4465_TNBC_panCK_cnv <- subset(x = Swarbrick_CID4465_TNBC_panCK_cnv, subset = cnv.ident == "Tumor")
# Quantify chromosomes with alterations per cell
Swarbrick_CID4465_TNBC_panCK_cnv_Pos_Chromosomes <- as.data.frame(rowSums(Swarbrick_CID4465_TNBC_panCK_cnv@meta.data[,startsWith(names(Swarbrick_CID4465_TNBC_panCK_cnv@meta.data),"has_cnv_")] > 0 ))
Swarbrick_CID4465_TNBC_panCK_cnv_Neg_Chromosomes <- as.data.frame(rowSums(Swarbrick_CID4465_TNBC_panCK_cnv@meta.data[,startsWith(names(Swarbrick_CID4465_TNBC_panCK_cnv@meta.data),"has_cnv_")] == 0))
Swarbrick_CID4465_TNBC_panCK_cnv_Chromosome_Quant <- as.data.frame(c(Swarbrick_CID4465_TNBC_panCK_cnv_Pos_Chromosomes, Swarbrick_CID4465_TNBC_panCK_cnv_Neg_Chromosomes))
colnames(Swarbrick_CID4465_TNBC_panCK_cnv_Chromosome_Quant) <- c("Chromosomes CNV Pos", "Chromosomes CNV Neg")
write.csv(Swarbrick_CID4465_TNBC_panCK_cnv_Chromosome_Quant, file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/Swarbrick_CID4465_TNBC_panCK/Swarbrick_CID4465_TNBC_panCK_cnv_Chromosome_Quant.csv")

# Save CNV Metadata
Swarbrick_CID4465_TNBC_panCK_cnv.meta.data <- as.data.frame(as.matrix(Swarbrick_CID4465_TNBC_panCK_cnv@meta.data))
write.csv(Swarbrick_CID4465_TNBC_panCK_cnv.meta.data, file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/Swarbrick_CID4465_TNBC_panCK/Swarbrick_CID4465_TNBC_panCK_cnv.meta.data.csv")

# Subset Cells with inferred CNVs -------------------------------------------------------------------------------------------------------------
Swarbrick_CID4465_TNBC_panCK_cnv <- subset(Swarbrick_CID4465_TNBC_panCK_cnv, cells=rownames(Swarbrick_CID4465_TNBC_panCK_cnv@meta.data)[rowSums(Swarbrick_CID4465_TNBC_panCK_cnv@meta.data[,startsWith(names(Swarbrick_CID4465_TNBC_panCK_cnv@meta.data),"has_cnv_")]) > 0 ])

# Save RDS with inferredCNVs
saveRDS(Swarbrick_CID4465_TNBC_panCK_cnv, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_inferCNV/Swarbrick_CID4465_TNBC_panCK_cnv.rds")


# Free memory
rm(Swarbrick_CID4465_TNBC_inferCNV)
rm(Swarbrick_CID4465_TNBC_Navin_Visvader_NORM_panCK_Reference)
rm(Swarbrick_CID4465_TNBC_panCK_cnv)
rm(Swarbrick_CID4465_TNBC_panCK_cnv_Pos_Chromosomes)
rm(Swarbrick_CID4465_TNBC_panCK_cnv_Neg_Chromosomes)
rm(Swarbrick_CID4465_TNBC_panCK_cnv_Chromosome_Quant)
rm(Swarbrick_CID4465_TNBC_panCK_cnv.meta.data)
gc()



################################################
#### InferCNV: Swarbrick_CID4495_TNBC_panCK ####
################################################

########## Merge with normal mammary glad, extract epithelium, and process via Seurat
# Load RDS
Swarbrick_CID4495_TNBC_panCK <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_panCK/Swarbrick_CID4495_TNBC_panCK.rds")
Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Navin_Visvader_NORM_panCK_Reference.rds")

# Set identities 
colnames(x = Swarbrick_CID4495_TNBC_panCK[[]])
Idents(object = Swarbrick_CID4495_TNBC_panCK) <- 'Tumor'
Swarbrick_CID4495_TNBC_panCK[["cnv.ident"]] <- Idents(object = Swarbrick_CID4495_TNBC_panCK)
levels(x = Swarbrick_CID4495_TNBC_panCK)

colnames(x = Navin_Visvader_NORM_panCK_Reference[[]])
Idents(object = Navin_Visvader_NORM_panCK_Reference) <- 'Normal'
Navin_Visvader_NORM_panCK_Reference[["cnv.ident"]] <- Idents(object = Navin_Visvader_NORM_panCK_Reference)
levels(x = Navin_Visvader_NORM_panCK_Reference)

# Merge subsets into recombined RDS
Swarbrick_CID4495_TNBC_Navin_Visvader_NORM_panCK_Reference <- merge(x = Swarbrick_CID4495_TNBC_panCK, y = Navin_Visvader_NORM_panCK_Reference)
# Join Layers
Swarbrick_CID4495_TNBC_Navin_Visvader_NORM_panCK_Reference <- JoinLayers(Swarbrick_CID4495_TNBC_Navin_Visvader_NORM_panCK_Reference)

# Remove subsetted objects
rm(Swarbrick_CID4495_TNBC_panCK)
rm(Navin_Visvader_NORM_panCK_Reference)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples before merge

# Find Variable Features
Swarbrick_CID4495_TNBC_Navin_Visvader_NORM_panCK_Reference <- FindVariableFeatures(Swarbrick_CID4495_TNBC_Navin_Visvader_NORM_panCK_Reference, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Swarbrick_CID4495_TNBC_Navin_Visvader_NORM_panCK_Reference)
Swarbrick_CID4495_TNBC_Navin_Visvader_NORM_panCK_Reference <- ScaleData(Swarbrick_CID4495_TNBC_Navin_Visvader_NORM_panCK_Reference, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################

Swarbrick_CID4495_TNBC_Navin_Visvader_NORM_panCK_Reference <- RunPCA(Swarbrick_CID4495_TNBC_Navin_Visvader_NORM_panCK_Reference, verbose = FALSE)

# Examine and visualize PCA results a few different ways
print(Swarbrick_CID4495_TNBC_Navin_Visvader_NORM_panCK_Reference[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Swarbrick_CID4495_TNBC_Navin_Visvader_NORM_panCK_Reference, dims = 1:2, reduction = "pca")
DimPlot(Swarbrick_CID4495_TNBC_Navin_Visvader_NORM_panCK_Reference, reduction = "pca")

Swarbrick_CID4495_TNBC_Navin_Visvader_NORM_panCK_Reference <- FindNeighbors(Swarbrick_CID4495_TNBC_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
Swarbrick_CID4495_TNBC_Navin_Visvader_NORM_panCK_Reference <- FindClusters(Swarbrick_CID4495_TNBC_Navin_Visvader_NORM_panCK_Reference, resolution = 0.1)

# Umap clustering
Swarbrick_CID4495_TNBC_Navin_Visvader_NORM_panCK_Reference <- RunUMAP(Swarbrick_CID4495_TNBC_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
DimPlot(Swarbrick_CID4495_TNBC_Navin_Visvader_NORM_panCK_Reference, reduction = "umap")
DimPlot(Swarbrick_CID4495_TNBC_Navin_Visvader_NORM_panCK_Reference, reduction = "umap", group.by = "orig.ident")

# Save RDS
saveRDS(Swarbrick_CID4495_TNBC_Navin_Visvader_NORM_panCK_Reference, file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_inferCNV_Input_RDS/Swarbrick_CID4495_TNBC_Navin_Visvader_NORM_panCK_Reference.rds")

########## Prepare Cell Annotations
# Extract & Save Meta Data
Swarbrick_CID4495_TNBC_Navin_Visvader_NORM_panCK_Reference.meta.data <- as.data.frame(as.matrix(Swarbrick_CID4495_TNBC_Navin_Visvader_NORM_panCK_Reference@meta.data))
# Extract pertinent columns
# Make row names first column
Swarbrick_CID4495_TNBC_Navin_Visvader_NORM_panCK_Reference.meta.data <- tibble::rownames_to_column(Swarbrick_CID4495_TNBC_Navin_Visvader_NORM_panCK_Reference.meta.data, "Cells")
# In Seurat, "X" becomes the cell names while "cnv.ident" is specified above
Swarbrick_CID4495_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations <- Swarbrick_CID4495_TNBC_Navin_Visvader_NORM_panCK_Reference.meta.data[, c("Cells", "cnv.ident")]
# Delete header
names(Swarbrick_CID4495_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations) <- NULL
# Save as table
write.table(Swarbrick_CID4495_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations,"/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_CID4495_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations.txt",sep="\t",row.names=FALSE)
# Remove metadata and annotations
rm(Swarbrick_CID4495_TNBC_Navin_Visvader_NORM_panCK_Reference.meta.data)
rm(Swarbrick_CID4495_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations)

########## Generate Count Matrix
# Extract Count Matrix
Swarbrick_CID4495_TNBC_Navin_Visvader_NORM_panCK_Reference.count <- Swarbrick_CID4495_TNBC_Navin_Visvader_NORM_panCK_Reference[["RNA"]]$counts
Swarbrick_CID4495_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as(Swarbrick_CID4495_TNBC_Navin_Visvader_NORM_panCK_Reference.count, 'sparseMatrix')
rm(Swarbrick_CID4495_TNBC_Navin_Visvader_NORM_panCK_Reference.count)
# Convert to requisite data format
Swarbrick_CID4495_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as.data.frame(Swarbrick_CID4495_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix)
# Save table
write.csv(Swarbrick_CID4495_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix, file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_inferCNV_Input_Count_Matrix/Swarbrick_CID4495_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv")
# Remove Objects
rm(Swarbrick_CID4495_TNBC_Navin_Visvader_NORM_panCK_Reference)
rm(Swarbrick_CID4495_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix)

########### Run inferCNV
# Note: check.names = FALSE prevents cell names from having dashes replaced with dots; essential to match annotations file
Swarbrick_CID4495_TNBC_inferCNV <- CreateInfercnvObject(raw_counts_matrix=read.csv(file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_inferCNV_Input_Count_Matrix/Swarbrick_CID4495_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv", header = TRUE, row.names = 1,check.names = FALSE),
                                              annotations_file=read.table(file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_CID4495_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations.txt", header = FALSE, row.names = 1),
                                              delim="\t",
                                              gene_order_file="/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/hg38_gencode_v27.txt",
                                              ref_group_names=c("Normal"))


# perform infercnv operations to reveal cnv signal
Swarbrick_CID4495_TNBC_inferCNV = infercnv::run(Swarbrick_CID4495_TNBC_inferCNV,
                                      cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                      out_dir="/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/Swarbrick_CID4495_TNBC_panCK/",  # dir is auto-created for storing outputs
                                      cluster_by_groups=T,   # cluster
                                      denoise=T,
                                      HMM=T,
                                      num_ref_groups = 9, analysis_mode = "subclusters")

###########  Integrate data with Seurat object
# Read RDS and add inferCNV results to metadata
Swarbrick_CID4495_TNBC_Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_inferCNV_Input_RDS/Swarbrick_CID4495_TNBC_Navin_Visvader_NORM_panCK_Reference.rds") 
Swarbrick_CID4495_TNBC_panCK_cnv <- infercnv::add_to_seurat(infercnv_output_path="/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/Swarbrick_CID4495_TNBC_panCK/",
                                                  seurat_obj=Swarbrick_CID4495_TNBC_Navin_Visvader_NORM_panCK_Reference, # optional
                                                  top_n=10)

# Collect Tumor cells from Seurat Object
Swarbrick_CID4495_TNBC_panCK_cnv <- subset(x = Swarbrick_CID4495_TNBC_panCK_cnv, subset = cnv.ident == "Tumor")
# Quantify chromosomes with alterations per cell
Swarbrick_CID4495_TNBC_panCK_cnv_Pos_Chromosomes <- as.data.frame(rowSums(Swarbrick_CID4495_TNBC_panCK_cnv@meta.data[,startsWith(names(Swarbrick_CID4495_TNBC_panCK_cnv@meta.data),"has_cnv_")] > 0 ))
Swarbrick_CID4495_TNBC_panCK_cnv_Neg_Chromosomes <- as.data.frame(rowSums(Swarbrick_CID4495_TNBC_panCK_cnv@meta.data[,startsWith(names(Swarbrick_CID4495_TNBC_panCK_cnv@meta.data),"has_cnv_")] == 0))
Swarbrick_CID4495_TNBC_panCK_cnv_Chromosome_Quant <- as.data.frame(c(Swarbrick_CID4495_TNBC_panCK_cnv_Pos_Chromosomes, Swarbrick_CID4495_TNBC_panCK_cnv_Neg_Chromosomes))
colnames(Swarbrick_CID4495_TNBC_panCK_cnv_Chromosome_Quant) <- c("Chromosomes CNV Pos", "Chromosomes CNV Neg")
write.csv(Swarbrick_CID4495_TNBC_panCK_cnv_Chromosome_Quant, file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/Swarbrick_CID4495_TNBC_panCK/Swarbrick_CID4495_TNBC_panCK_cnv_Chromosome_Quant.csv")

# Save CNV Metadata
Swarbrick_CID4495_TNBC_panCK_cnv.meta.data <- as.data.frame(as.matrix(Swarbrick_CID4495_TNBC_panCK_cnv@meta.data))
write.csv(Swarbrick_CID4495_TNBC_panCK_cnv.meta.data, file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/Swarbrick_CID4495_TNBC_panCK/Swarbrick_CID4495_TNBC_panCK_cnv.meta.data.csv")

# Subset Cells with inferred CNVs -------------------------------------------------------------------------------------------------------------
Swarbrick_CID4495_TNBC_panCK_cnv <- subset(Swarbrick_CID4495_TNBC_panCK_cnv, cells=rownames(Swarbrick_CID4495_TNBC_panCK_cnv@meta.data)[rowSums(Swarbrick_CID4495_TNBC_panCK_cnv@meta.data[,startsWith(names(Swarbrick_CID4495_TNBC_panCK_cnv@meta.data),"has_cnv_")]) > 0 ])

# Save RDS with inferredCNVs
saveRDS(Swarbrick_CID4495_TNBC_panCK_cnv, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_inferCNV/Swarbrick_CID4495_TNBC_panCK_cnv.rds")


# Free memory
rm(Swarbrick_CID4495_TNBC_inferCNV)
rm(Swarbrick_CID4495_TNBC_Navin_Visvader_NORM_panCK_Reference)
rm(Swarbrick_CID4495_TNBC_panCK_cnv)
rm(Swarbrick_CID4495_TNBC_panCK_cnv_Pos_Chromosomes)
rm(Swarbrick_CID4495_TNBC_panCK_cnv_Neg_Chromosomes)
rm(Swarbrick_CID4495_TNBC_panCK_cnv_Chromosome_Quant)
rm(Swarbrick_CID4495_TNBC_panCK_cnv.meta.data)
gc()




################################################
#### InferCNV: Swarbrick_CID4513_TNBC_panCK ####
################################################

########## Merge with normal mammary glad, extract epithelium, and process via Seurat
# Load RDS
Swarbrick_CID4513_TNBC_panCK <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_panCK/Swarbrick_CID4513_TNBC_panCK.rds")
Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Navin_Visvader_NORM_panCK_Reference.rds")

# Set identities 
colnames(x = Swarbrick_CID4513_TNBC_panCK[[]])
Idents(object = Swarbrick_CID4513_TNBC_panCK) <- 'Tumor'
Swarbrick_CID4513_TNBC_panCK[["cnv.ident"]] <- Idents(object = Swarbrick_CID4513_TNBC_panCK)
levels(x = Swarbrick_CID4513_TNBC_panCK)

colnames(x = Navin_Visvader_NORM_panCK_Reference[[]])
Idents(object = Navin_Visvader_NORM_panCK_Reference) <- 'Normal'
Navin_Visvader_NORM_panCK_Reference[["cnv.ident"]] <- Idents(object = Navin_Visvader_NORM_panCK_Reference)
levels(x = Navin_Visvader_NORM_panCK_Reference)

# Merge subsets into recombined RDS
Swarbrick_CID4513_TNBC_Navin_Visvader_NORM_panCK_Reference <- merge(x = Swarbrick_CID4513_TNBC_panCK, y = Navin_Visvader_NORM_panCK_Reference)
# Join Layers
Swarbrick_CID4513_TNBC_Navin_Visvader_NORM_panCK_Reference <- JoinLayers(Swarbrick_CID4513_TNBC_Navin_Visvader_NORM_panCK_Reference)

# Remove subsetted objects
rm(Swarbrick_CID4513_TNBC_panCK)
rm(Navin_Visvader_NORM_panCK_Reference)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples before merge

# Find Variable Features
Swarbrick_CID4513_TNBC_Navin_Visvader_NORM_panCK_Reference <- FindVariableFeatures(Swarbrick_CID4513_TNBC_Navin_Visvader_NORM_panCK_Reference, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Swarbrick_CID4513_TNBC_Navin_Visvader_NORM_panCK_Reference)
Swarbrick_CID4513_TNBC_Navin_Visvader_NORM_panCK_Reference <- ScaleData(Swarbrick_CID4513_TNBC_Navin_Visvader_NORM_panCK_Reference, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################

Swarbrick_CID4513_TNBC_Navin_Visvader_NORM_panCK_Reference <- RunPCA(Swarbrick_CID4513_TNBC_Navin_Visvader_NORM_panCK_Reference, verbose = FALSE)

# Examine and visualize PCA results a few different ways
print(Swarbrick_CID4513_TNBC_Navin_Visvader_NORM_panCK_Reference[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Swarbrick_CID4513_TNBC_Navin_Visvader_NORM_panCK_Reference, dims = 1:2, reduction = "pca")
DimPlot(Swarbrick_CID4513_TNBC_Navin_Visvader_NORM_panCK_Reference, reduction = "pca")

Swarbrick_CID4513_TNBC_Navin_Visvader_NORM_panCK_Reference <- FindNeighbors(Swarbrick_CID4513_TNBC_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
Swarbrick_CID4513_TNBC_Navin_Visvader_NORM_panCK_Reference <- FindClusters(Swarbrick_CID4513_TNBC_Navin_Visvader_NORM_panCK_Reference, resolution = 0.1)

# Umap clustering
Swarbrick_CID4513_TNBC_Navin_Visvader_NORM_panCK_Reference <- RunUMAP(Swarbrick_CID4513_TNBC_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
DimPlot(Swarbrick_CID4513_TNBC_Navin_Visvader_NORM_panCK_Reference, reduction = "umap")
DimPlot(Swarbrick_CID4513_TNBC_Navin_Visvader_NORM_panCK_Reference, reduction = "umap", group.by = "orig.ident")

# Save RDS
saveRDS(Swarbrick_CID4513_TNBC_Navin_Visvader_NORM_panCK_Reference, file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_inferCNV_Input_RDS/Swarbrick_CID4513_TNBC_Navin_Visvader_NORM_panCK_Reference.rds")

########## Prepare Cell Annotations
# Extract & Save Meta Data
Swarbrick_CID4513_TNBC_Navin_Visvader_NORM_panCK_Reference.meta.data <- as.data.frame(as.matrix(Swarbrick_CID4513_TNBC_Navin_Visvader_NORM_panCK_Reference@meta.data))
# Extract pertinent columns
# Make row names first column
Swarbrick_CID4513_TNBC_Navin_Visvader_NORM_panCK_Reference.meta.data <- tibble::rownames_to_column(Swarbrick_CID4513_TNBC_Navin_Visvader_NORM_panCK_Reference.meta.data, "Cells")
# In Seurat, "X" becomes the cell names while "cnv.ident" is specified above
Swarbrick_CID4513_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations <- Swarbrick_CID4513_TNBC_Navin_Visvader_NORM_panCK_Reference.meta.data[, c("Cells", "cnv.ident")]
# Delete header
names(Swarbrick_CID4513_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations) <- NULL
# Save as table
write.table(Swarbrick_CID4513_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations,"/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_CID4513_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations.txt",sep="\t",row.names=FALSE)
# Remove metadata and annotations
rm(Swarbrick_CID4513_TNBC_Navin_Visvader_NORM_panCK_Reference.meta.data)
rm(Swarbrick_CID4513_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations)

########## Generate Count Matrix
# Extract Count Matrix
Swarbrick_CID4513_TNBC_Navin_Visvader_NORM_panCK_Reference.count <- Swarbrick_CID4513_TNBC_Navin_Visvader_NORM_panCK_Reference[["RNA"]]$counts
Swarbrick_CID4513_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as(Swarbrick_CID4513_TNBC_Navin_Visvader_NORM_panCK_Reference.count, 'sparseMatrix')
rm(Swarbrick_CID4513_TNBC_Navin_Visvader_NORM_panCK_Reference.count)
# Convert to requisite data format
Swarbrick_CID4513_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as.data.frame(Swarbrick_CID4513_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix)
# Save table
write.csv(Swarbrick_CID4513_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix, file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_inferCNV_Input_Count_Matrix/Swarbrick_CID4513_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv")
# Remove Objects
rm(Swarbrick_CID4513_TNBC_Navin_Visvader_NORM_panCK_Reference)
rm(Swarbrick_CID4513_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix)

########### Run inferCNV
# Note: check.names = FALSE prevents cell names from having dashes replaced with dots; essential to match annotations file
Swarbrick_CID4513_TNBC_inferCNV <- CreateInfercnvObject(raw_counts_matrix=read.csv(file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_inferCNV_Input_Count_Matrix/Swarbrick_CID4513_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv", header = TRUE, row.names = 1,check.names = FALSE),
                                              annotations_file=read.table(file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_CID4513_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations.txt", header = FALSE, row.names = 1),
                                              delim="\t",
                                              gene_order_file="/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/hg38_gencode_v27.txt",
                                              ref_group_names=c("Normal"))


# perform infercnv operations to reveal cnv signal
Swarbrick_CID4513_TNBC_inferCNV = infercnv::run(Swarbrick_CID4513_TNBC_inferCNV,
                                      cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                      out_dir="/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/Swarbrick_CID4513_TNBC_panCK/",  # dir is auto-created for storing outputs
                                      cluster_by_groups=T,   # cluster
                                      denoise=T,
                                      HMM=T,
                                      num_ref_groups = 9, analysis_mode = "subclusters")

###########  Integrate data with Seurat object
# Read RDS and add inferCNV results to metadata
Swarbrick_CID4513_TNBC_Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_inferCNV_Input_RDS/Swarbrick_CID4513_TNBC_Navin_Visvader_NORM_panCK_Reference.rds") 
Swarbrick_CID4513_TNBC_panCK_cnv <- infercnv::add_to_seurat(infercnv_output_path="/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/Swarbrick_CID4513_TNBC_panCK/",
                                                  seurat_obj=Swarbrick_CID4513_TNBC_Navin_Visvader_NORM_panCK_Reference, # optional
                                                  top_n=10)

# Collect Tumor cells from Seurat Object
Swarbrick_CID4513_TNBC_panCK_cnv <- subset(x = Swarbrick_CID4513_TNBC_panCK_cnv, subset = cnv.ident == "Tumor")
# Quantify chromosomes with alterations per cell
Swarbrick_CID4513_TNBC_panCK_cnv_Pos_Chromosomes <- as.data.frame(rowSums(Swarbrick_CID4513_TNBC_panCK_cnv@meta.data[,startsWith(names(Swarbrick_CID4513_TNBC_panCK_cnv@meta.data),"has_cnv_")] > 0 ))
Swarbrick_CID4513_TNBC_panCK_cnv_Neg_Chromosomes <- as.data.frame(rowSums(Swarbrick_CID4513_TNBC_panCK_cnv@meta.data[,startsWith(names(Swarbrick_CID4513_TNBC_panCK_cnv@meta.data),"has_cnv_")] == 0))
Swarbrick_CID4513_TNBC_panCK_cnv_Chromosome_Quant <- as.data.frame(c(Swarbrick_CID4513_TNBC_panCK_cnv_Pos_Chromosomes, Swarbrick_CID4513_TNBC_panCK_cnv_Neg_Chromosomes))
colnames(Swarbrick_CID4513_TNBC_panCK_cnv_Chromosome_Quant) <- c("Chromosomes CNV Pos", "Chromosomes CNV Neg")
write.csv(Swarbrick_CID4513_TNBC_panCK_cnv_Chromosome_Quant, file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/Swarbrick_CID4513_TNBC_panCK/Swarbrick_CID4513_TNBC_panCK_cnv_Chromosome_Quant.csv")

# Save CNV Metadata
Swarbrick_CID4513_TNBC_panCK_cnv.meta.data <- as.data.frame(as.matrix(Swarbrick_CID4513_TNBC_panCK_cnv@meta.data))
write.csv(Swarbrick_CID4513_TNBC_panCK_cnv.meta.data, file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/Swarbrick_CID4513_TNBC_panCK/Swarbrick_CID4513_TNBC_panCK_cnv.meta.data.csv")

# Subset Cells with inferred CNVs -------------------------------------------------------------------------------------------------------------
Swarbrick_CID4513_TNBC_panCK_cnv <- subset(Swarbrick_CID4513_TNBC_panCK_cnv, cells=rownames(Swarbrick_CID4513_TNBC_panCK_cnv@meta.data)[rowSums(Swarbrick_CID4513_TNBC_panCK_cnv@meta.data[,startsWith(names(Swarbrick_CID4513_TNBC_panCK_cnv@meta.data),"has_cnv_")]) > 0 ])

# Save RDS with inferredCNVs
saveRDS(Swarbrick_CID4513_TNBC_panCK_cnv, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_inferCNV/Swarbrick_CID4513_TNBC_panCK_cnv.rds")


# Free memory
rm(Swarbrick_CID4513_TNBC_inferCNV)
rm(Swarbrick_CID4513_TNBC_Navin_Visvader_NORM_panCK_Reference)
rm(Swarbrick_CID4513_TNBC_panCK_cnv)
rm(Swarbrick_CID4513_TNBC_panCK_cnv_Pos_Chromosomes)
rm(Swarbrick_CID4513_TNBC_panCK_cnv_Neg_Chromosomes)
rm(Swarbrick_CID4513_TNBC_panCK_cnv_Chromosome_Quant)
rm(Swarbrick_CID4513_TNBC_panCK_cnv.meta.data)
gc()



################################################
#### InferCNV: Swarbrick_CID4515_TNBC_panCK ####
################################################

########## Merge with normal mammary glad, extract epithelium, and process via Seurat
# Load RDS
Swarbrick_CID4515_TNBC_panCK <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_panCK/Swarbrick_CID4515_TNBC_panCK.rds")
Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Navin_Visvader_NORM_panCK_Reference.rds")

# Set identities 
colnames(x = Swarbrick_CID4515_TNBC_panCK[[]])
Idents(object = Swarbrick_CID4515_TNBC_panCK) <- 'Tumor'
Swarbrick_CID4515_TNBC_panCK[["cnv.ident"]] <- Idents(object = Swarbrick_CID4515_TNBC_panCK)
levels(x = Swarbrick_CID4515_TNBC_panCK)

colnames(x = Navin_Visvader_NORM_panCK_Reference[[]])
Idents(object = Navin_Visvader_NORM_panCK_Reference) <- 'Normal'
Navin_Visvader_NORM_panCK_Reference[["cnv.ident"]] <- Idents(object = Navin_Visvader_NORM_panCK_Reference)
levels(x = Navin_Visvader_NORM_panCK_Reference)

# Merge subsets into recombined RDS
Swarbrick_CID4515_TNBC_Navin_Visvader_NORM_panCK_Reference <- merge(x = Swarbrick_CID4515_TNBC_panCK, y = Navin_Visvader_NORM_panCK_Reference)
# Join Layers
Swarbrick_CID4515_TNBC_Navin_Visvader_NORM_panCK_Reference <- JoinLayers(Swarbrick_CID4515_TNBC_Navin_Visvader_NORM_panCK_Reference)

# Remove subsetted objects
rm(Swarbrick_CID4515_TNBC_panCK)
rm(Navin_Visvader_NORM_panCK_Reference)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples before merge

# Find Variable Features
Swarbrick_CID4515_TNBC_Navin_Visvader_NORM_panCK_Reference <- FindVariableFeatures(Swarbrick_CID4515_TNBC_Navin_Visvader_NORM_panCK_Reference, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Swarbrick_CID4515_TNBC_Navin_Visvader_NORM_panCK_Reference)
Swarbrick_CID4515_TNBC_Navin_Visvader_NORM_panCK_Reference <- ScaleData(Swarbrick_CID4515_TNBC_Navin_Visvader_NORM_panCK_Reference, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################

Swarbrick_CID4515_TNBC_Navin_Visvader_NORM_panCK_Reference <- RunPCA(Swarbrick_CID4515_TNBC_Navin_Visvader_NORM_panCK_Reference, verbose = FALSE)

# Examine and visualize PCA results a few different ways
print(Swarbrick_CID4515_TNBC_Navin_Visvader_NORM_panCK_Reference[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Swarbrick_CID4515_TNBC_Navin_Visvader_NORM_panCK_Reference, dims = 1:2, reduction = "pca")
DimPlot(Swarbrick_CID4515_TNBC_Navin_Visvader_NORM_panCK_Reference, reduction = "pca")

Swarbrick_CID4515_TNBC_Navin_Visvader_NORM_panCK_Reference <- FindNeighbors(Swarbrick_CID4515_TNBC_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
Swarbrick_CID4515_TNBC_Navin_Visvader_NORM_panCK_Reference <- FindClusters(Swarbrick_CID4515_TNBC_Navin_Visvader_NORM_panCK_Reference, resolution = 0.1)

# Umap clustering
Swarbrick_CID4515_TNBC_Navin_Visvader_NORM_panCK_Reference <- RunUMAP(Swarbrick_CID4515_TNBC_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
DimPlot(Swarbrick_CID4515_TNBC_Navin_Visvader_NORM_panCK_Reference, reduction = "umap")
DimPlot(Swarbrick_CID4515_TNBC_Navin_Visvader_NORM_panCK_Reference, reduction = "umap", group.by = "orig.ident")

# Save RDS
saveRDS(Swarbrick_CID4515_TNBC_Navin_Visvader_NORM_panCK_Reference, file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_inferCNV_Input_RDS/Swarbrick_CID4515_TNBC_Navin_Visvader_NORM_panCK_Reference.rds")

########## Prepare Cell Annotations
# Extract & Save Meta Data
Swarbrick_CID4515_TNBC_Navin_Visvader_NORM_panCK_Reference.meta.data <- as.data.frame(as.matrix(Swarbrick_CID4515_TNBC_Navin_Visvader_NORM_panCK_Reference@meta.data))
# Extract pertinent columns
# Make row names first column
Swarbrick_CID4515_TNBC_Navin_Visvader_NORM_panCK_Reference.meta.data <- tibble::rownames_to_column(Swarbrick_CID4515_TNBC_Navin_Visvader_NORM_panCK_Reference.meta.data, "Cells")
# In Seurat, "X" becomes the cell names while "cnv.ident" is specified above
Swarbrick_CID4515_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations <- Swarbrick_CID4515_TNBC_Navin_Visvader_NORM_panCK_Reference.meta.data[, c("Cells", "cnv.ident")]
# Delete header
names(Swarbrick_CID4515_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations) <- NULL
# Save as table
write.table(Swarbrick_CID4515_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations,"/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_CID4515_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations.txt",sep="\t",row.names=FALSE)
# Remove metadata and annotations
rm(Swarbrick_CID4515_TNBC_Navin_Visvader_NORM_panCK_Reference.meta.data)
rm(Swarbrick_CID4515_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations)

########## Generate Count Matrix
# Extract Count Matrix
Swarbrick_CID4515_TNBC_Navin_Visvader_NORM_panCK_Reference.count <- Swarbrick_CID4515_TNBC_Navin_Visvader_NORM_panCK_Reference[["RNA"]]$counts
Swarbrick_CID4515_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as(Swarbrick_CID4515_TNBC_Navin_Visvader_NORM_panCK_Reference.count, 'sparseMatrix')
rm(Swarbrick_CID4515_TNBC_Navin_Visvader_NORM_panCK_Reference.count)
# Convert to requisite data format
Swarbrick_CID4515_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as.data.frame(Swarbrick_CID4515_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix)
# Save table
write.csv(Swarbrick_CID4515_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix, file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_inferCNV_Input_Count_Matrix/Swarbrick_CID4515_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv")
# Remove Objects
rm(Swarbrick_CID4515_TNBC_Navin_Visvader_NORM_panCK_Reference)
rm(Swarbrick_CID4515_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix)

########### Run inferCNV
# Note: check.names = FALSE prevents cell names from having dashes replaced with dots; essential to match annotations file
Swarbrick_CID4515_TNBC_inferCNV <- CreateInfercnvObject(raw_counts_matrix=read.csv(file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_inferCNV_Input_Count_Matrix/Swarbrick_CID4515_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv", header = TRUE, row.names = 1,check.names = FALSE),
                                              annotations_file=read.table(file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_CID4515_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations.txt", header = FALSE, row.names = 1),
                                              delim="\t",
                                              gene_order_file="/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/hg38_gencode_v27.txt",
                                              ref_group_names=c("Normal"))


# perform infercnv operations to reveal cnv signal
Swarbrick_CID4515_TNBC_inferCNV = infercnv::run(Swarbrick_CID4515_TNBC_inferCNV,
                                      cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                      out_dir="/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/Swarbrick_CID4515_TNBC_panCK/",  # dir is auto-created for storing outputs
                                      cluster_by_groups=T,   # cluster
                                      denoise=T,
                                      HMM=T,
                                      num_ref_groups = 9, analysis_mode = "subclusters")

###########  Integrate data with Seurat object
# Read RDS and add inferCNV results to metadata
Swarbrick_CID4515_TNBC_Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_inferCNV_Input_RDS/Swarbrick_CID4515_TNBC_Navin_Visvader_NORM_panCK_Reference.rds") 
Swarbrick_CID4515_TNBC_panCK_cnv <- infercnv::add_to_seurat(infercnv_output_path="/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/Swarbrick_CID4515_TNBC_panCK/",
                                                  seurat_obj=Swarbrick_CID4515_TNBC_Navin_Visvader_NORM_panCK_Reference, # optional
                                                  top_n=10)

# Collect Tumor cells from Seurat Object
Swarbrick_CID4515_TNBC_panCK_cnv <- subset(x = Swarbrick_CID4515_TNBC_panCK_cnv, subset = cnv.ident == "Tumor")
# Quantify chromosomes with alterations per cell
Swarbrick_CID4515_TNBC_panCK_cnv_Pos_Chromosomes <- as.data.frame(rowSums(Swarbrick_CID4515_TNBC_panCK_cnv@meta.data[,startsWith(names(Swarbrick_CID4515_TNBC_panCK_cnv@meta.data),"has_cnv_")] > 0 ))
Swarbrick_CID4515_TNBC_panCK_cnv_Neg_Chromosomes <- as.data.frame(rowSums(Swarbrick_CID4515_TNBC_panCK_cnv@meta.data[,startsWith(names(Swarbrick_CID4515_TNBC_panCK_cnv@meta.data),"has_cnv_")] == 0))
Swarbrick_CID4515_TNBC_panCK_cnv_Chromosome_Quant <- as.data.frame(c(Swarbrick_CID4515_TNBC_panCK_cnv_Pos_Chromosomes, Swarbrick_CID4515_TNBC_panCK_cnv_Neg_Chromosomes))
colnames(Swarbrick_CID4515_TNBC_panCK_cnv_Chromosome_Quant) <- c("Chromosomes CNV Pos", "Chromosomes CNV Neg")
write.csv(Swarbrick_CID4515_TNBC_panCK_cnv_Chromosome_Quant, file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/Swarbrick_CID4515_TNBC_panCK/Swarbrick_CID4515_TNBC_panCK_cnv_Chromosome_Quant.csv")

# Save CNV Metadata
Swarbrick_CID4515_TNBC_panCK_cnv.meta.data <- as.data.frame(as.matrix(Swarbrick_CID4515_TNBC_panCK_cnv@meta.data))
write.csv(Swarbrick_CID4515_TNBC_panCK_cnv.meta.data, file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/Swarbrick_CID4515_TNBC_panCK/Swarbrick_CID4515_TNBC_panCK_cnv.meta.data.csv")

# Subset Cells with inferred CNVs -------------------------------------------------------------------------------------------------------------
Swarbrick_CID4515_TNBC_panCK_cnv <- subset(Swarbrick_CID4515_TNBC_panCK_cnv, cells=rownames(Swarbrick_CID4515_TNBC_panCK_cnv@meta.data)[rowSums(Swarbrick_CID4515_TNBC_panCK_cnv@meta.data[,startsWith(names(Swarbrick_CID4515_TNBC_panCK_cnv@meta.data),"has_cnv_")]) > 0 ])

# Save RDS with inferredCNVs
saveRDS(Swarbrick_CID4515_TNBC_panCK_cnv, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_inferCNV/Swarbrick_CID4515_TNBC_panCK_cnv.rds")


# Free memory
rm(Swarbrick_CID4515_TNBC_inferCNV)
rm(Swarbrick_CID4515_TNBC_Navin_Visvader_NORM_panCK_Reference)
rm(Swarbrick_CID4515_TNBC_panCK_cnv)
rm(Swarbrick_CID4515_TNBC_panCK_cnv_Pos_Chromosomes)
rm(Swarbrick_CID4515_TNBC_panCK_cnv_Neg_Chromosomes)
rm(Swarbrick_CID4515_TNBC_panCK_cnv_Chromosome_Quant)
rm(Swarbrick_CID4515_TNBC_panCK_cnv.meta.data)
gc()



################################################
#### InferCNV: Swarbrick_CID4523_TNBC_panCK ####
################################################

########## Merge with normal mammary glad, extract epithelium, and process via Seurat
# Load RDS
Swarbrick_CID4523_TNBC_panCK <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_panCK/Swarbrick_CID4523_TNBC_panCK.rds")
Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Navin_Visvader_NORM_panCK_Reference.rds")

# Set identities 
colnames(x = Swarbrick_CID4523_TNBC_panCK[[]])
Idents(object = Swarbrick_CID4523_TNBC_panCK) <- 'Tumor'
Swarbrick_CID4523_TNBC_panCK[["cnv.ident"]] <- Idents(object = Swarbrick_CID4523_TNBC_panCK)
levels(x = Swarbrick_CID4523_TNBC_panCK)

colnames(x = Navin_Visvader_NORM_panCK_Reference[[]])
Idents(object = Navin_Visvader_NORM_panCK_Reference) <- 'Normal'
Navin_Visvader_NORM_panCK_Reference[["cnv.ident"]] <- Idents(object = Navin_Visvader_NORM_panCK_Reference)
levels(x = Navin_Visvader_NORM_panCK_Reference)

# Merge subsets into recombined RDS
Swarbrick_CID4523_TNBC_Navin_Visvader_NORM_panCK_Reference <- merge(x = Swarbrick_CID4523_TNBC_panCK, y = Navin_Visvader_NORM_panCK_Reference)
# Join Layers
Swarbrick_CID4523_TNBC_Navin_Visvader_NORM_panCK_Reference <- JoinLayers(Swarbrick_CID4523_TNBC_Navin_Visvader_NORM_panCK_Reference)

# Remove subsetted objects
rm(Swarbrick_CID4523_TNBC_panCK)
rm(Navin_Visvader_NORM_panCK_Reference)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples before merge

# Find Variable Features
Swarbrick_CID4523_TNBC_Navin_Visvader_NORM_panCK_Reference <- FindVariableFeatures(Swarbrick_CID4523_TNBC_Navin_Visvader_NORM_panCK_Reference, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Swarbrick_CID4523_TNBC_Navin_Visvader_NORM_panCK_Reference)
Swarbrick_CID4523_TNBC_Navin_Visvader_NORM_panCK_Reference <- ScaleData(Swarbrick_CID4523_TNBC_Navin_Visvader_NORM_panCK_Reference, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################

Swarbrick_CID4523_TNBC_Navin_Visvader_NORM_panCK_Reference <- RunPCA(Swarbrick_CID4523_TNBC_Navin_Visvader_NORM_panCK_Reference, verbose = FALSE)

# Examine and visualize PCA results a few different ways
print(Swarbrick_CID4523_TNBC_Navin_Visvader_NORM_panCK_Reference[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Swarbrick_CID4523_TNBC_Navin_Visvader_NORM_panCK_Reference, dims = 1:2, reduction = "pca")
DimPlot(Swarbrick_CID4523_TNBC_Navin_Visvader_NORM_panCK_Reference, reduction = "pca")

Swarbrick_CID4523_TNBC_Navin_Visvader_NORM_panCK_Reference <- FindNeighbors(Swarbrick_CID4523_TNBC_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
Swarbrick_CID4523_TNBC_Navin_Visvader_NORM_panCK_Reference <- FindClusters(Swarbrick_CID4523_TNBC_Navin_Visvader_NORM_panCK_Reference, resolution = 0.1)

# Umap clustering
Swarbrick_CID4523_TNBC_Navin_Visvader_NORM_panCK_Reference <- RunUMAP(Swarbrick_CID4523_TNBC_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
DimPlot(Swarbrick_CID4523_TNBC_Navin_Visvader_NORM_panCK_Reference, reduction = "umap")
DimPlot(Swarbrick_CID4523_TNBC_Navin_Visvader_NORM_panCK_Reference, reduction = "umap", group.by = "orig.ident")

# Save RDS
saveRDS(Swarbrick_CID4523_TNBC_Navin_Visvader_NORM_panCK_Reference, file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_inferCNV_Input_RDS/Swarbrick_CID4523_TNBC_Navin_Visvader_NORM_panCK_Reference.rds")

########## Prepare Cell Annotations
# Extract & Save Meta Data
Swarbrick_CID4523_TNBC_Navin_Visvader_NORM_panCK_Reference.meta.data <- as.data.frame(as.matrix(Swarbrick_CID4523_TNBC_Navin_Visvader_NORM_panCK_Reference@meta.data))
# Extract pertinent columns
# Make row names first column
Swarbrick_CID4523_TNBC_Navin_Visvader_NORM_panCK_Reference.meta.data <- tibble::rownames_to_column(Swarbrick_CID4523_TNBC_Navin_Visvader_NORM_panCK_Reference.meta.data, "Cells")
# In Seurat, "X" becomes the cell names while "cnv.ident" is specified above
Swarbrick_CID4523_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations <- Swarbrick_CID4523_TNBC_Navin_Visvader_NORM_panCK_Reference.meta.data[, c("Cells", "cnv.ident")]
# Delete header
names(Swarbrick_CID4523_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations) <- NULL
# Save as table
write.table(Swarbrick_CID4523_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations,"/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_CID4523_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations.txt",sep="\t",row.names=FALSE)
# Remove metadata and annotations
rm(Swarbrick_CID4523_TNBC_Navin_Visvader_NORM_panCK_Reference.meta.data)
rm(Swarbrick_CID4523_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations)

########## Generate Count Matrix
# Extract Count Matrix
Swarbrick_CID4523_TNBC_Navin_Visvader_NORM_panCK_Reference.count <- Swarbrick_CID4523_TNBC_Navin_Visvader_NORM_panCK_Reference[["RNA"]]$counts
Swarbrick_CID4523_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as(Swarbrick_CID4523_TNBC_Navin_Visvader_NORM_panCK_Reference.count, 'sparseMatrix')
rm(Swarbrick_CID4523_TNBC_Navin_Visvader_NORM_panCK_Reference.count)
# Convert to requisite data format
Swarbrick_CID4523_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as.data.frame(Swarbrick_CID4523_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix)
# Save table
write.csv(Swarbrick_CID4523_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix, file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_inferCNV_Input_Count_Matrix/Swarbrick_CID4523_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv")
# Remove Objects
rm(Swarbrick_CID4523_TNBC_Navin_Visvader_NORM_panCK_Reference)
rm(Swarbrick_CID4523_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix)

########### Run inferCNV
# Note: check.names = FALSE prevents cell names from having dashes replaced with dots; essential to match annotations file
Swarbrick_CID4523_TNBC_inferCNV <- CreateInfercnvObject(raw_counts_matrix=read.csv(file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_inferCNV_Input_Count_Matrix/Swarbrick_CID4523_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv", header = TRUE, row.names = 1,check.names = FALSE),
                                              annotations_file=read.table(file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_CID4523_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations.txt", header = FALSE, row.names = 1),
                                              delim="\t",
                                              gene_order_file="/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/hg38_gencode_v27.txt",
                                              ref_group_names=c("Normal"))


# perform infercnv operations to reveal cnv signal
Swarbrick_CID4523_TNBC_inferCNV = infercnv::run(Swarbrick_CID4523_TNBC_inferCNV,
                                      cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                      out_dir="/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/Swarbrick_CID4523_TNBC_panCK/",  # dir is auto-created for storing outputs
                                      cluster_by_groups=T,   # cluster
                                      denoise=T,
                                      HMM=T,
                                      num_ref_groups = 9, analysis_mode = "subclusters")

###########  Integrate data with Seurat object
# Read RDS and add inferCNV results to metadata
Swarbrick_CID4523_TNBC_Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_inferCNV_Input_RDS/Swarbrick_CID4523_TNBC_Navin_Visvader_NORM_panCK_Reference.rds") 
Swarbrick_CID4523_TNBC_panCK_cnv <- infercnv::add_to_seurat(infercnv_output_path="/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/Swarbrick_CID4523_TNBC_panCK/",
                                                  seurat_obj=Swarbrick_CID4523_TNBC_Navin_Visvader_NORM_panCK_Reference, # optional
                                                  top_n=10)

# Collect Tumor cells from Seurat Object
Swarbrick_CID4523_TNBC_panCK_cnv <- subset(x = Swarbrick_CID4523_TNBC_panCK_cnv, subset = cnv.ident == "Tumor")
# Quantify chromosomes with alterations per cell
Swarbrick_CID4523_TNBC_panCK_cnv_Pos_Chromosomes <- as.data.frame(rowSums(Swarbrick_CID4523_TNBC_panCK_cnv@meta.data[,startsWith(names(Swarbrick_CID4523_TNBC_panCK_cnv@meta.data),"has_cnv_")] > 0 ))
Swarbrick_CID4523_TNBC_panCK_cnv_Neg_Chromosomes <- as.data.frame(rowSums(Swarbrick_CID4523_TNBC_panCK_cnv@meta.data[,startsWith(names(Swarbrick_CID4523_TNBC_panCK_cnv@meta.data),"has_cnv_")] == 0))
Swarbrick_CID4523_TNBC_panCK_cnv_Chromosome_Quant <- as.data.frame(c(Swarbrick_CID4523_TNBC_panCK_cnv_Pos_Chromosomes, Swarbrick_CID4523_TNBC_panCK_cnv_Neg_Chromosomes))
colnames(Swarbrick_CID4523_TNBC_panCK_cnv_Chromosome_Quant) <- c("Chromosomes CNV Pos", "Chromosomes CNV Neg")
write.csv(Swarbrick_CID4523_TNBC_panCK_cnv_Chromosome_Quant, file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/Swarbrick_CID4523_TNBC_panCK/Swarbrick_CID4523_TNBC_panCK_cnv_Chromosome_Quant.csv")

# Save CNV Metadata
Swarbrick_CID4523_TNBC_panCK_cnv.meta.data <- as.data.frame(as.matrix(Swarbrick_CID4523_TNBC_panCK_cnv@meta.data))
write.csv(Swarbrick_CID4523_TNBC_panCK_cnv.meta.data, file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/Swarbrick_CID4523_TNBC_panCK/Swarbrick_CID4523_TNBC_panCK_cnv.meta.data.csv")

# Subset Cells with inferred CNVs -------------------------------------------------------------------------------------------------------------
Swarbrick_CID4523_TNBC_panCK_cnv <- subset(Swarbrick_CID4523_TNBC_panCK_cnv, cells=rownames(Swarbrick_CID4523_TNBC_panCK_cnv@meta.data)[rowSums(Swarbrick_CID4523_TNBC_panCK_cnv@meta.data[,startsWith(names(Swarbrick_CID4523_TNBC_panCK_cnv@meta.data),"has_cnv_")]) > 0 ])

# Save RDS with inferredCNVs
saveRDS(Swarbrick_CID4523_TNBC_panCK_cnv, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_inferCNV/Swarbrick_CID4523_TNBC_panCK_cnv.rds")


# Free memory
rm(Swarbrick_CID4523_TNBC_inferCNV)
rm(Swarbrick_CID4523_TNBC_Navin_Visvader_NORM_panCK_Reference)
rm(Swarbrick_CID4523_TNBC_panCK_cnv)
rm(Swarbrick_CID4523_TNBC_panCK_cnv_Pos_Chromosomes)
rm(Swarbrick_CID4523_TNBC_panCK_cnv_Neg_Chromosomes)
rm(Swarbrick_CID4523_TNBC_panCK_cnv_Chromosome_Quant)
rm(Swarbrick_CID4523_TNBC_panCK_cnv.meta.data)
gc()



#################################################
#### InferCNV: Swarbrick_CID44041_TNBC_panCK ####
#################################################

########## Merge with normal mammary glad, extract epithelium, and process via Seurat
# Load RDS
Swarbrick_CID44041_TNBC_panCK <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_panCK/Swarbrick_CID44041_TNBC_panCK.rds")
Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Navin_Visvader_NORM_panCK_Reference.rds")

# Set identities 
colnames(x = Swarbrick_CID44041_TNBC_panCK[[]])
Idents(object = Swarbrick_CID44041_TNBC_panCK) <- 'Tumor'
Swarbrick_CID44041_TNBC_panCK[["cnv.ident"]] <- Idents(object = Swarbrick_CID44041_TNBC_panCK)
levels(x = Swarbrick_CID44041_TNBC_panCK)

colnames(x = Navin_Visvader_NORM_panCK_Reference[[]])
Idents(object = Navin_Visvader_NORM_panCK_Reference) <- 'Normal'
Navin_Visvader_NORM_panCK_Reference[["cnv.ident"]] <- Idents(object = Navin_Visvader_NORM_panCK_Reference)
levels(x = Navin_Visvader_NORM_panCK_Reference)

# Merge subsets into recombined RDS
Swarbrick_CID44041_TNBC_Navin_Visvader_NORM_panCK_Reference <- merge(x = Swarbrick_CID44041_TNBC_panCK, y = Navin_Visvader_NORM_panCK_Reference)
# Join Layers
Swarbrick_CID44041_TNBC_Navin_Visvader_NORM_panCK_Reference <- JoinLayers(Swarbrick_CID44041_TNBC_Navin_Visvader_NORM_panCK_Reference)

# Remove subsetted objects
rm(Swarbrick_CID44041_TNBC_panCK)
rm(Navin_Visvader_NORM_panCK_Reference)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples before merge

# Find Variable Features
Swarbrick_CID44041_TNBC_Navin_Visvader_NORM_panCK_Reference <- FindVariableFeatures(Swarbrick_CID44041_TNBC_Navin_Visvader_NORM_panCK_Reference, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Swarbrick_CID44041_TNBC_Navin_Visvader_NORM_panCK_Reference)
Swarbrick_CID44041_TNBC_Navin_Visvader_NORM_panCK_Reference <- ScaleData(Swarbrick_CID44041_TNBC_Navin_Visvader_NORM_panCK_Reference, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################

Swarbrick_CID44041_TNBC_Navin_Visvader_NORM_panCK_Reference <- RunPCA(Swarbrick_CID44041_TNBC_Navin_Visvader_NORM_panCK_Reference, verbose = FALSE)

# Examine and visualize PCA results a few different ways
print(Swarbrick_CID44041_TNBC_Navin_Visvader_NORM_panCK_Reference[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Swarbrick_CID44041_TNBC_Navin_Visvader_NORM_panCK_Reference, dims = 1:2, reduction = "pca")
DimPlot(Swarbrick_CID44041_TNBC_Navin_Visvader_NORM_panCK_Reference, reduction = "pca")

Swarbrick_CID44041_TNBC_Navin_Visvader_NORM_panCK_Reference <- FindNeighbors(Swarbrick_CID44041_TNBC_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
Swarbrick_CID44041_TNBC_Navin_Visvader_NORM_panCK_Reference <- FindClusters(Swarbrick_CID44041_TNBC_Navin_Visvader_NORM_panCK_Reference, resolution = 0.1)

# Umap clustering
Swarbrick_CID44041_TNBC_Navin_Visvader_NORM_panCK_Reference <- RunUMAP(Swarbrick_CID44041_TNBC_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
DimPlot(Swarbrick_CID44041_TNBC_Navin_Visvader_NORM_panCK_Reference, reduction = "umap")
DimPlot(Swarbrick_CID44041_TNBC_Navin_Visvader_NORM_panCK_Reference, reduction = "umap", group.by = "orig.ident")

# Save RDS
saveRDS(Swarbrick_CID44041_TNBC_Navin_Visvader_NORM_panCK_Reference, file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_inferCNV_Input_RDS/Swarbrick_CID44041_TNBC_Navin_Visvader_NORM_panCK_Reference.rds")

########## Prepare Cell Annotations
# Extract & Save Meta Data
Swarbrick_CID44041_TNBC_Navin_Visvader_NORM_panCK_Reference.meta.data <- as.data.frame(as.matrix(Swarbrick_CID44041_TNBC_Navin_Visvader_NORM_panCK_Reference@meta.data))
# Extract pertinent columns
# Make row names first column
Swarbrick_CID44041_TNBC_Navin_Visvader_NORM_panCK_Reference.meta.data <- tibble::rownames_to_column(Swarbrick_CID44041_TNBC_Navin_Visvader_NORM_panCK_Reference.meta.data, "Cells")
# In Seurat, "X" becomes the cell names while "cnv.ident" is specified above
Swarbrick_CID44041_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations <- Swarbrick_CID44041_TNBC_Navin_Visvader_NORM_panCK_Reference.meta.data[, c("Cells", "cnv.ident")]
# Delete header
names(Swarbrick_CID44041_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations) <- NULL
# Save as table
write.table(Swarbrick_CID44041_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations,"/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_CID44041_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations.txt",sep="\t",row.names=FALSE)
# Remove metadata and annotations
rm(Swarbrick_CID44041_TNBC_Navin_Visvader_NORM_panCK_Reference.meta.data)
rm(Swarbrick_CID44041_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations)

########## Generate Count Matrix
# Extract Count Matrix
Swarbrick_CID44041_TNBC_Navin_Visvader_NORM_panCK_Reference.count <- Swarbrick_CID44041_TNBC_Navin_Visvader_NORM_panCK_Reference[["RNA"]]$counts
Swarbrick_CID44041_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as(Swarbrick_CID44041_TNBC_Navin_Visvader_NORM_panCK_Reference.count, 'sparseMatrix')
rm(Swarbrick_CID44041_TNBC_Navin_Visvader_NORM_panCK_Reference.count)
# Convert to requisite data format
Swarbrick_CID44041_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as.data.frame(Swarbrick_CID44041_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix)
# Save table
write.csv(Swarbrick_CID44041_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix, file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_inferCNV_Input_Count_Matrix/Swarbrick_CID44041_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv")
# Remove Objects
rm(Swarbrick_CID44041_TNBC_Navin_Visvader_NORM_panCK_Reference)
rm(Swarbrick_CID44041_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix)

########### Run inferCNV
# Note: check.names = FALSE prevents cell names from having dashes replaced with dots; essential to match annotations file
Swarbrick_CID44041_TNBC_inferCNV <- CreateInfercnvObject(raw_counts_matrix=read.csv(file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_inferCNV_Input_Count_Matrix/Swarbrick_CID44041_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv", header = TRUE, row.names = 1,check.names = FALSE),
                                               annotations_file=read.table(file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_CID44041_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations.txt", header = FALSE, row.names = 1),
                                               delim="\t",
                                               gene_order_file="/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/hg38_gencode_v27.txt",
                                               ref_group_names=c("Normal"))


# perform infercnv operations to reveal cnv signal
Swarbrick_CID44041_TNBC_inferCNV = infercnv::run(Swarbrick_CID44041_TNBC_inferCNV,
                                       cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                       out_dir="/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/Swarbrick_CID44041_TNBC_panCK/",  # dir is auto-created for storing outputs
                                       cluster_by_groups=T,   # cluster
                                       denoise=T,
                                       HMM=T,
                                       num_ref_groups = 9, analysis_mode = "subclusters")

###########  Integrate data with Seurat object
# Read RDS and add inferCNV results to metadata
Swarbrick_CID44041_TNBC_Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_inferCNV_Input_RDS/Swarbrick_CID44041_TNBC_Navin_Visvader_NORM_panCK_Reference.rds") 
Swarbrick_CID44041_TNBC_panCK_cnv <- infercnv::add_to_seurat(infercnv_output_path="/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/Swarbrick_CID44041_TNBC_panCK/",
                                                   seurat_obj=Swarbrick_CID44041_TNBC_Navin_Visvader_NORM_panCK_Reference, # optional
                                                   top_n=10)

# Collect Tumor cells from Seurat Object
Swarbrick_CID44041_TNBC_panCK_cnv <- subset(x = Swarbrick_CID44041_TNBC_panCK_cnv, subset = cnv.ident == "Tumor")
# Quantify chromosomes with alterations per cell
Swarbrick_CID44041_TNBC_panCK_cnv_Pos_Chromosomes <- as.data.frame(rowSums(Swarbrick_CID44041_TNBC_panCK_cnv@meta.data[,startsWith(names(Swarbrick_CID44041_TNBC_panCK_cnv@meta.data),"has_cnv_")] > 0 ))
Swarbrick_CID44041_TNBC_panCK_cnv_Neg_Chromosomes <- as.data.frame(rowSums(Swarbrick_CID44041_TNBC_panCK_cnv@meta.data[,startsWith(names(Swarbrick_CID44041_TNBC_panCK_cnv@meta.data),"has_cnv_")] == 0))
Swarbrick_CID44041_TNBC_panCK_cnv_Chromosome_Quant <- as.data.frame(c(Swarbrick_CID44041_TNBC_panCK_cnv_Pos_Chromosomes, Swarbrick_CID44041_TNBC_panCK_cnv_Neg_Chromosomes))
colnames(Swarbrick_CID44041_TNBC_panCK_cnv_Chromosome_Quant) <- c("Chromosomes CNV Pos", "Chromosomes CNV Neg")
write.csv(Swarbrick_CID44041_TNBC_panCK_cnv_Chromosome_Quant, file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/Swarbrick_CID44041_TNBC_panCK/Swarbrick_CID44041_TNBC_panCK_cnv_Chromosome_Quant.csv")

# Save CNV Metadata
Swarbrick_CID44041_TNBC_panCK_cnv.meta.data <- as.data.frame(as.matrix(Swarbrick_CID44041_TNBC_panCK_cnv@meta.data))
write.csv(Swarbrick_CID44041_TNBC_panCK_cnv.meta.data, file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/Swarbrick_CID44041_TNBC_panCK/Swarbrick_CID44041_TNBC_panCK_cnv.meta.data.csv")

# Subset Cells with inferred CNVs -------------------------------------------------------------------------------------------------------------
Swarbrick_CID44041_TNBC_panCK_cnv <- subset(Swarbrick_CID44041_TNBC_panCK_cnv, cells=rownames(Swarbrick_CID44041_TNBC_panCK_cnv@meta.data)[rowSums(Swarbrick_CID44041_TNBC_panCK_cnv@meta.data[,startsWith(names(Swarbrick_CID44041_TNBC_panCK_cnv@meta.data),"has_cnv_")]) > 0 ])

# Save RDS with inferredCNVs
saveRDS(Swarbrick_CID44041_TNBC_panCK_cnv, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_inferCNV/Swarbrick_CID44041_TNBC_panCK_cnv.rds")


# Free memory
rm(Swarbrick_CID44041_TNBC_inferCNV)
rm(Swarbrick_CID44041_TNBC_Navin_Visvader_NORM_panCK_Reference)
rm(Swarbrick_CID44041_TNBC_panCK_cnv)
rm(Swarbrick_CID44041_TNBC_panCK_cnv_Pos_Chromosomes)
rm(Swarbrick_CID44041_TNBC_panCK_cnv_Neg_Chromosomes)
rm(Swarbrick_CID44041_TNBC_panCK_cnv_Chromosome_Quant)
rm(Swarbrick_CID44041_TNBC_panCK_cnv.meta.data)
gc()



#################################################
#### InferCNV: Swarbrick_CID44971_TNBC_panCK ####
#################################################

########## Merge with normal mammary glad, extract epithelium, and process via Seurat
# Load RDS
Swarbrick_CID44971_TNBC_panCK <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_panCK/Swarbrick_CID44971_TNBC_panCK.rds")
Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Navin_Visvader_NORM_panCK_Reference.rds")

# Set identities 
colnames(x = Swarbrick_CID44971_TNBC_panCK[[]])
Idents(object = Swarbrick_CID44971_TNBC_panCK) <- 'Tumor'
Swarbrick_CID44971_TNBC_panCK[["cnv.ident"]] <- Idents(object = Swarbrick_CID44971_TNBC_panCK)
levels(x = Swarbrick_CID44971_TNBC_panCK)

colnames(x = Navin_Visvader_NORM_panCK_Reference[[]])
Idents(object = Navin_Visvader_NORM_panCK_Reference) <- 'Normal'
Navin_Visvader_NORM_panCK_Reference[["cnv.ident"]] <- Idents(object = Navin_Visvader_NORM_panCK_Reference)
levels(x = Navin_Visvader_NORM_panCK_Reference)

# Merge subsets into recombined RDS
Swarbrick_CID44971_TNBC_Navin_Visvader_NORM_panCK_Reference <- merge(x = Swarbrick_CID44971_TNBC_panCK, y = Navin_Visvader_NORM_panCK_Reference)
# Join Layers
Swarbrick_CID44971_TNBC_Navin_Visvader_NORM_panCK_Reference <- JoinLayers(Swarbrick_CID44971_TNBC_Navin_Visvader_NORM_panCK_Reference)

# Remove subsetted objects
rm(Swarbrick_CID44971_TNBC_panCK)
rm(Navin_Visvader_NORM_panCK_Reference)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples before merge

# Find Variable Features
Swarbrick_CID44971_TNBC_Navin_Visvader_NORM_panCK_Reference <- FindVariableFeatures(Swarbrick_CID44971_TNBC_Navin_Visvader_NORM_panCK_Reference, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Swarbrick_CID44971_TNBC_Navin_Visvader_NORM_panCK_Reference)
Swarbrick_CID44971_TNBC_Navin_Visvader_NORM_panCK_Reference <- ScaleData(Swarbrick_CID44971_TNBC_Navin_Visvader_NORM_panCK_Reference, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################

Swarbrick_CID44971_TNBC_Navin_Visvader_NORM_panCK_Reference <- RunPCA(Swarbrick_CID44971_TNBC_Navin_Visvader_NORM_panCK_Reference, verbose = FALSE)

# Examine and visualize PCA results a few different ways
print(Swarbrick_CID44971_TNBC_Navin_Visvader_NORM_panCK_Reference[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Swarbrick_CID44971_TNBC_Navin_Visvader_NORM_panCK_Reference, dims = 1:2, reduction = "pca")
DimPlot(Swarbrick_CID44971_TNBC_Navin_Visvader_NORM_panCK_Reference, reduction = "pca")

Swarbrick_CID44971_TNBC_Navin_Visvader_NORM_panCK_Reference <- FindNeighbors(Swarbrick_CID44971_TNBC_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
Swarbrick_CID44971_TNBC_Navin_Visvader_NORM_panCK_Reference <- FindClusters(Swarbrick_CID44971_TNBC_Navin_Visvader_NORM_panCK_Reference, resolution = 0.1)

# Umap clustering
Swarbrick_CID44971_TNBC_Navin_Visvader_NORM_panCK_Reference <- RunUMAP(Swarbrick_CID44971_TNBC_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
DimPlot(Swarbrick_CID44971_TNBC_Navin_Visvader_NORM_panCK_Reference, reduction = "umap")
DimPlot(Swarbrick_CID44971_TNBC_Navin_Visvader_NORM_panCK_Reference, reduction = "umap", group.by = "orig.ident")

# Save RDS
saveRDS(Swarbrick_CID44971_TNBC_Navin_Visvader_NORM_panCK_Reference, file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_inferCNV_Input_RDS/Swarbrick_CID44971_TNBC_Navin_Visvader_NORM_panCK_Reference.rds")

########## Prepare Cell Annotations
# Extract & Save Meta Data
Swarbrick_CID44971_TNBC_Navin_Visvader_NORM_panCK_Reference.meta.data <- as.data.frame(as.matrix(Swarbrick_CID44971_TNBC_Navin_Visvader_NORM_panCK_Reference@meta.data))
# Extract pertinent columns
# Make row names first column
Swarbrick_CID44971_TNBC_Navin_Visvader_NORM_panCK_Reference.meta.data <- tibble::rownames_to_column(Swarbrick_CID44971_TNBC_Navin_Visvader_NORM_panCK_Reference.meta.data, "Cells")
# In Seurat, "X" becomes the cell names while "cnv.ident" is specified above
Swarbrick_CID44971_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations <- Swarbrick_CID44971_TNBC_Navin_Visvader_NORM_panCK_Reference.meta.data[, c("Cells", "cnv.ident")]
# Delete header
names(Swarbrick_CID44971_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations) <- NULL
# Save as table
write.table(Swarbrick_CID44971_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations,"/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_CID44971_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations.txt",sep="\t",row.names=FALSE)
# Remove metadata and annotations
rm(Swarbrick_CID44971_TNBC_Navin_Visvader_NORM_panCK_Reference.meta.data)
rm(Swarbrick_CID44971_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations)

########## Generate Count Matrix
# Extract Count Matrix
Swarbrick_CID44971_TNBC_Navin_Visvader_NORM_panCK_Reference.count <- Swarbrick_CID44971_TNBC_Navin_Visvader_NORM_panCK_Reference[["RNA"]]$counts
Swarbrick_CID44971_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as(Swarbrick_CID44971_TNBC_Navin_Visvader_NORM_panCK_Reference.count, 'sparseMatrix')
rm(Swarbrick_CID44971_TNBC_Navin_Visvader_NORM_panCK_Reference.count)
# Convert to requisite data format
Swarbrick_CID44971_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as.data.frame(Swarbrick_CID44971_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix)
# Save table
write.csv(Swarbrick_CID44971_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix, file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_inferCNV_Input_Count_Matrix/Swarbrick_CID44971_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv")
# Remove Objects
rm(Swarbrick_CID44971_TNBC_Navin_Visvader_NORM_panCK_Reference)
rm(Swarbrick_CID44971_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix)

########### Run inferCNV
# Note: check.names = FALSE prevents cell names from having dashes replaced with dots; essential to match annotations file
Swarbrick_CID44971_TNBC_inferCNV <- CreateInfercnvObject(raw_counts_matrix=read.csv(file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_inferCNV_Input_Count_Matrix/Swarbrick_CID44971_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv", header = TRUE, row.names = 1,check.names = FALSE),
                                               annotations_file=read.table(file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_CID44971_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations.txt", header = FALSE, row.names = 1),
                                               delim="\t",
                                               gene_order_file="/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/hg38_gencode_v27.txt",
                                               ref_group_names=c("Normal"))


# perform infercnv operations to reveal cnv signal
Swarbrick_CID44971_TNBC_inferCNV = infercnv::run(Swarbrick_CID44971_TNBC_inferCNV,
                                       cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                       out_dir="/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/Swarbrick_CID44971_TNBC_panCK/",  # dir is auto-created for storing outputs
                                       cluster_by_groups=T,   # cluster
                                       denoise=T,
                                       HMM=T,
                                       num_ref_groups = 9, analysis_mode = "subclusters")

###########  Integrate data with Seurat object
# Read RDS and add inferCNV results to metadata
Swarbrick_CID44971_TNBC_Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_inferCNV_Input_RDS/Swarbrick_CID44971_TNBC_Navin_Visvader_NORM_panCK_Reference.rds") 
Swarbrick_CID44971_TNBC_panCK_cnv <- infercnv::add_to_seurat(infercnv_output_path="/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/Swarbrick_CID44971_TNBC_panCK/",
                                                   seurat_obj=Swarbrick_CID44971_TNBC_Navin_Visvader_NORM_panCK_Reference, # optional
                                                   top_n=10)

# Collect Tumor cells from Seurat Object
Swarbrick_CID44971_TNBC_panCK_cnv <- subset(x = Swarbrick_CID44971_TNBC_panCK_cnv, subset = cnv.ident == "Tumor")
# Quantify chromosomes with alterations per cell
Swarbrick_CID44971_TNBC_panCK_cnv_Pos_Chromosomes <- as.data.frame(rowSums(Swarbrick_CID44971_TNBC_panCK_cnv@meta.data[,startsWith(names(Swarbrick_CID44971_TNBC_panCK_cnv@meta.data),"has_cnv_")] > 0 ))
Swarbrick_CID44971_TNBC_panCK_cnv_Neg_Chromosomes <- as.data.frame(rowSums(Swarbrick_CID44971_TNBC_panCK_cnv@meta.data[,startsWith(names(Swarbrick_CID44971_TNBC_panCK_cnv@meta.data),"has_cnv_")] == 0))
Swarbrick_CID44971_TNBC_panCK_cnv_Chromosome_Quant <- as.data.frame(c(Swarbrick_CID44971_TNBC_panCK_cnv_Pos_Chromosomes, Swarbrick_CID44971_TNBC_panCK_cnv_Neg_Chromosomes))
colnames(Swarbrick_CID44971_TNBC_panCK_cnv_Chromosome_Quant) <- c("Chromosomes CNV Pos", "Chromosomes CNV Neg")
write.csv(Swarbrick_CID44971_TNBC_panCK_cnv_Chromosome_Quant, file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/Swarbrick_CID44971_TNBC_panCK/Swarbrick_CID44971_TNBC_panCK_cnv_Chromosome_Quant.csv")

# Save CNV Metadata
Swarbrick_CID44971_TNBC_panCK_cnv.meta.data <- as.data.frame(as.matrix(Swarbrick_CID44971_TNBC_panCK_cnv@meta.data))
write.csv(Swarbrick_CID44971_TNBC_panCK_cnv.meta.data, file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/Swarbrick_CID44971_TNBC_panCK/Swarbrick_CID44971_TNBC_panCK_cnv.meta.data.csv")

# Subset Cells with inferred CNVs -------------------------------------------------------------------------------------------------------------
Swarbrick_CID44971_TNBC_panCK_cnv <- subset(Swarbrick_CID44971_TNBC_panCK_cnv, cells=rownames(Swarbrick_CID44971_TNBC_panCK_cnv@meta.data)[rowSums(Swarbrick_CID44971_TNBC_panCK_cnv@meta.data[,startsWith(names(Swarbrick_CID44971_TNBC_panCK_cnv@meta.data),"has_cnv_")]) > 0 ])

# Save RDS with inferredCNVs
saveRDS(Swarbrick_CID44971_TNBC_panCK_cnv, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_inferCNV/Swarbrick_CID44971_TNBC_panCK_cnv.rds")


# Free memory
rm(Swarbrick_CID44971_TNBC_inferCNV)
rm(Swarbrick_CID44971_TNBC_Navin_Visvader_NORM_panCK_Reference)
rm(Swarbrick_CID44971_TNBC_panCK_cnv)
rm(Swarbrick_CID44971_TNBC_panCK_cnv_Pos_Chromosomes)
rm(Swarbrick_CID44971_TNBC_panCK_cnv_Neg_Chromosomes)
rm(Swarbrick_CID44971_TNBC_panCK_cnv_Chromosome_Quant)
rm(Swarbrick_CID44971_TNBC_panCK_cnv.meta.data)
gc()





#################################################
#### InferCNV: Swarbrick_CID44991_TNBC_panCK ####
#################################################

########## Merge with normal mammary glad, extract epithelium, and process via Seurat
# Load RDS
Swarbrick_CID44991_TNBC_panCK <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_panCK/Swarbrick_CID44991_TNBC_panCK.rds")
Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Visvader/Visvader_inferCNV/Visvader_inferCNV_Input/Visvader_inferCNV_Input_RDS/Navin_Visvader_NORM_panCK_Reference.rds")

# Set identities 
colnames(x = Swarbrick_CID44991_TNBC_panCK[[]])
Idents(object = Swarbrick_CID44991_TNBC_panCK) <- 'Tumor'
Swarbrick_CID44991_TNBC_panCK[["cnv.ident"]] <- Idents(object = Swarbrick_CID44991_TNBC_panCK)
levels(x = Swarbrick_CID44991_TNBC_panCK)

colnames(x = Navin_Visvader_NORM_panCK_Reference[[]])
Idents(object = Navin_Visvader_NORM_panCK_Reference) <- 'Normal'
Navin_Visvader_NORM_panCK_Reference[["cnv.ident"]] <- Idents(object = Navin_Visvader_NORM_panCK_Reference)
levels(x = Navin_Visvader_NORM_panCK_Reference)

# Merge subsets into recombined RDS
Swarbrick_CID44991_TNBC_Navin_Visvader_NORM_panCK_Reference <- merge(x = Swarbrick_CID44991_TNBC_panCK, y = Navin_Visvader_NORM_panCK_Reference)
# Join Layers
Swarbrick_CID44991_TNBC_Navin_Visvader_NORM_panCK_Reference <- JoinLayers(Swarbrick_CID44991_TNBC_Navin_Visvader_NORM_panCK_Reference)

# Remove subsetted objects
rm(Swarbrick_CID44991_TNBC_panCK)
rm(Navin_Visvader_NORM_panCK_Reference)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples before merge

# Find Variable Features
Swarbrick_CID44991_TNBC_Navin_Visvader_NORM_panCK_Reference <- FindVariableFeatures(Swarbrick_CID44991_TNBC_Navin_Visvader_NORM_panCK_Reference, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Swarbrick_CID44991_TNBC_Navin_Visvader_NORM_panCK_Reference)
Swarbrick_CID44991_TNBC_Navin_Visvader_NORM_panCK_Reference <- ScaleData(Swarbrick_CID44991_TNBC_Navin_Visvader_NORM_panCK_Reference, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################

Swarbrick_CID44991_TNBC_Navin_Visvader_NORM_panCK_Reference <- RunPCA(Swarbrick_CID44991_TNBC_Navin_Visvader_NORM_panCK_Reference, verbose = FALSE)

# Examine and visualize PCA results a few different ways
print(Swarbrick_CID44991_TNBC_Navin_Visvader_NORM_panCK_Reference[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Swarbrick_CID44991_TNBC_Navin_Visvader_NORM_panCK_Reference, dims = 1:2, reduction = "pca")
DimPlot(Swarbrick_CID44991_TNBC_Navin_Visvader_NORM_panCK_Reference, reduction = "pca")

Swarbrick_CID44991_TNBC_Navin_Visvader_NORM_panCK_Reference <- FindNeighbors(Swarbrick_CID44991_TNBC_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
Swarbrick_CID44991_TNBC_Navin_Visvader_NORM_panCK_Reference <- FindClusters(Swarbrick_CID44991_TNBC_Navin_Visvader_NORM_panCK_Reference, resolution = 0.1)

# Umap clustering
Swarbrick_CID44991_TNBC_Navin_Visvader_NORM_panCK_Reference <- RunUMAP(Swarbrick_CID44991_TNBC_Navin_Visvader_NORM_panCK_Reference, dims = 1:20)
DimPlot(Swarbrick_CID44991_TNBC_Navin_Visvader_NORM_panCK_Reference, reduction = "umap")
DimPlot(Swarbrick_CID44991_TNBC_Navin_Visvader_NORM_panCK_Reference, reduction = "umap", group.by = "orig.ident")

# Save RDS
saveRDS(Swarbrick_CID44991_TNBC_Navin_Visvader_NORM_panCK_Reference, file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_inferCNV_Input_RDS/Swarbrick_CID44991_TNBC_Navin_Visvader_NORM_panCK_Reference.rds")

########## Prepare Cell Annotations
# Extract & Save Meta Data
Swarbrick_CID44991_TNBC_Navin_Visvader_NORM_panCK_Reference.meta.data <- as.data.frame(as.matrix(Swarbrick_CID44991_TNBC_Navin_Visvader_NORM_panCK_Reference@meta.data))
# Extract pertinent columns
# Make row names first column
Swarbrick_CID44991_TNBC_Navin_Visvader_NORM_panCK_Reference.meta.data <- tibble::rownames_to_column(Swarbrick_CID44991_TNBC_Navin_Visvader_NORM_panCK_Reference.meta.data, "Cells")
# In Seurat, "X" becomes the cell names while "cnv.ident" is specified above
Swarbrick_CID44991_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations <- Swarbrick_CID44991_TNBC_Navin_Visvader_NORM_panCK_Reference.meta.data[, c("Cells", "cnv.ident")]
# Delete header
names(Swarbrick_CID44991_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations) <- NULL
# Save as table
write.table(Swarbrick_CID44991_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations,"/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_CID44991_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations.txt",sep="\t",row.names=FALSE)
# Remove metadata and annotations
rm(Swarbrick_CID44991_TNBC_Navin_Visvader_NORM_panCK_Reference.meta.data)
rm(Swarbrick_CID44991_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations)

########## Generate Count Matrix
# Extract Count Matrix
Swarbrick_CID44991_TNBC_Navin_Visvader_NORM_panCK_Reference.count <- Swarbrick_CID44991_TNBC_Navin_Visvader_NORM_panCK_Reference[["RNA"]]$counts
Swarbrick_CID44991_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as(Swarbrick_CID44991_TNBC_Navin_Visvader_NORM_panCK_Reference.count, 'sparseMatrix')
rm(Swarbrick_CID44991_TNBC_Navin_Visvader_NORM_panCK_Reference.count)
# Convert to requisite data format
Swarbrick_CID44991_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix <- as.data.frame(Swarbrick_CID44991_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix)
# Save table
write.csv(Swarbrick_CID44991_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix, file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_inferCNV_Input_Count_Matrix/Swarbrick_CID44991_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv")
# Remove Objects
rm(Swarbrick_CID44991_TNBC_Navin_Visvader_NORM_panCK_Reference)
rm(Swarbrick_CID44991_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix)

########### Run inferCNV
# Note: check.names = FALSE prevents cell names from having dashes replaced with dots; essential to match annotations file
Swarbrick_CID44991_TNBC_inferCNV <- CreateInfercnvObject(raw_counts_matrix=read.csv(file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_inferCNV_Input_Count_Matrix/Swarbrick_CID44991_TNBC_Navin_Visvader_NORM_panCK_Reference.count.matrix.csv", header = TRUE, row.names = 1,check.names = FALSE),
                                               annotations_file=read.table(file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_CID44991_TNBC_Navin_Visvader_NORM_panCK_Reference.annotations.txt", header = FALSE, row.names = 1),
                                               delim="\t",
                                               gene_order_file="/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/hg38_gencode_v27.txt",
                                               ref_group_names=c("Normal"))


# perform infercnv operations to reveal cnv signal
Swarbrick_CID44991_TNBC_inferCNV = infercnv::run(Swarbrick_CID44991_TNBC_inferCNV,
                                       cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                       out_dir="/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/Swarbrick_CID44991_TNBC_panCK/",  # dir is auto-created for storing outputs
                                       cluster_by_groups=T,   # cluster
                                       denoise=T,
                                       HMM=T,
                                       num_ref_groups = 9, analysis_mode = "subclusters")

###########  Integrate data with Seurat object
# Read RDS and add inferCNV results to metadata
Swarbrick_CID44991_TNBC_Navin_Visvader_NORM_panCK_Reference <- readRDS(file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Input/Swarbrick_inferCNV_Input_RDS/Swarbrick_CID44991_TNBC_Navin_Visvader_NORM_panCK_Reference.rds") 
Swarbrick_CID44991_TNBC_panCK_cnv <- infercnv::add_to_seurat(infercnv_output_path="/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/Swarbrick_CID44991_TNBC_panCK/",
                                                   seurat_obj=Swarbrick_CID44991_TNBC_Navin_Visvader_NORM_panCK_Reference, # optional
                                                   top_n=10)

# Collect Tumor cells from Seurat Object
Swarbrick_CID44991_TNBC_panCK_cnv <- subset(x = Swarbrick_CID44991_TNBC_panCK_cnv, subset = cnv.ident == "Tumor")
# Quantify chromosomes with alterations per cell
Swarbrick_CID44991_TNBC_panCK_cnv_Pos_Chromosomes <- as.data.frame(rowSums(Swarbrick_CID44991_TNBC_panCK_cnv@meta.data[,startsWith(names(Swarbrick_CID44991_TNBC_panCK_cnv@meta.data),"has_cnv_")] > 0 ))
Swarbrick_CID44991_TNBC_panCK_cnv_Neg_Chromosomes <- as.data.frame(rowSums(Swarbrick_CID44991_TNBC_panCK_cnv@meta.data[,startsWith(names(Swarbrick_CID44991_TNBC_panCK_cnv@meta.data),"has_cnv_")] == 0))
Swarbrick_CID44991_TNBC_panCK_cnv_Chromosome_Quant <- as.data.frame(c(Swarbrick_CID44991_TNBC_panCK_cnv_Pos_Chromosomes, Swarbrick_CID44991_TNBC_panCK_cnv_Neg_Chromosomes))
colnames(Swarbrick_CID44991_TNBC_panCK_cnv_Chromosome_Quant) <- c("Chromosomes CNV Pos", "Chromosomes CNV Neg")
write.csv(Swarbrick_CID44991_TNBC_panCK_cnv_Chromosome_Quant, file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/Swarbrick_CID44991_TNBC_panCK/Swarbrick_CID44991_TNBC_panCK_cnv_Chromosome_Quant.csv")

# Save CNV Metadata
Swarbrick_CID44991_TNBC_panCK_cnv.meta.data <- as.data.frame(as.matrix(Swarbrick_CID44991_TNBC_panCK_cnv@meta.data))
write.csv(Swarbrick_CID44991_TNBC_panCK_cnv.meta.data, file = "/R/R_Swarbrick/Swarbrick_inferCNV/Swarbrick_inferCNV_Output/Swarbrick_CID44991_TNBC_panCK/Swarbrick_CID44991_TNBC_panCK_cnv.meta.data.csv")

# Subset Cells with inferred CNVs -------------------------------------------------------------------------------------------------------------
Swarbrick_CID44991_TNBC_panCK_cnv <- subset(Swarbrick_CID44991_TNBC_panCK_cnv, cells=rownames(Swarbrick_CID44991_TNBC_panCK_cnv@meta.data)[rowSums(Swarbrick_CID44991_TNBC_panCK_cnv@meta.data[,startsWith(names(Swarbrick_CID44991_TNBC_panCK_cnv@meta.data),"has_cnv_")]) > 0 ])

# Save RDS with inferredCNVs
saveRDS(Swarbrick_CID44991_TNBC_panCK_cnv, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_inferCNV/Swarbrick_CID44991_TNBC_panCK_cnv.rds")


# Free memory
rm(Swarbrick_CID44991_TNBC_inferCNV)
rm(Swarbrick_CID44991_TNBC_Navin_Visvader_NORM_panCK_Reference)
rm(Swarbrick_CID44991_TNBC_panCK_cnv)
rm(Swarbrick_CID44991_TNBC_panCK_cnv_Pos_Chromosomes)
rm(Swarbrick_CID44991_TNBC_panCK_cnv_Neg_Chromosomes)
rm(Swarbrick_CID44991_TNBC_panCK_cnv_Chromosome_Quant)
rm(Swarbrick_CID44991_TNBC_panCK_cnv.meta.data)
gc()





#############################################################################
# Step 6: Merge Neoplastic Cells from Different Tumors into Combined Object #  
#############################################################################
# Load Libraries
library(Seurat)
library(dplyr)
library(patchwork)


#########################################################
############# This is ER Tumor Batch # ##################
#########################################################

#########################
# Load Individual Files #
#########################
Swarbrick_CID3941_ER_panCK_cnv <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_inferCNV/Swarbrick_CID3941_ER_panCK_cnv.rds")
Swarbrick_CID3948_ER_panCK_cnv <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_inferCNV/Swarbrick_CID3948_ER_panCK_cnv.rds")
Swarbrick_CID4040_ER_panCK_cnv <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_inferCNV/Swarbrick_CID4040_ER_panCK_cnv.rds")
Swarbrick_CID4067_ER_panCK_cnv <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_inferCNV/Swarbrick_CID4067_ER_panCK_cnv.rds")
Swarbrick_CID4290A_ER_panCK_cnv <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_inferCNV/Swarbrick_CID4290A_ER_panCK_cnv.rds")
Swarbrick_CID4398_ER_panCK_cnv <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_inferCNV/Swarbrick_CID4398_ER_panCK_cnv.rds")
Swarbrick_CID4461_ER_panCK_cnv <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_inferCNV/Swarbrick_CID4461_ER_panCK_cnv.rds")
Swarbrick_CID4463_ER_panCK_cnv <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_inferCNV/Swarbrick_CID4463_ER_panCK_cnv.rds")
Swarbrick_CID4471_ER_panCK_cnv <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_inferCNV/Swarbrick_CID4471_ER_panCK_cnv.rds")
Swarbrick_CID4530N_ER_panCK_cnv <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_inferCNV/Swarbrick_CID4530N_ER_panCK_cnv.rds")
Swarbrick_CID4535_ER_panCK_cnv <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_inferCNV/Swarbrick_CID4535_ER_panCK_cnv.rds")


##############################################
# Merge Individual Files into Single Dataset #
##############################################
# Merge Datasets
Swarbrick_merged_ER_panCK_cnv <- merge(x = Swarbrick_CID3941_ER_panCK_cnv, y = c(Swarbrick_CID3948_ER_panCK_cnv, Swarbrick_CID4040_ER_panCK_cnv, Swarbrick_CID4067_ER_panCK_cnv,
                                                                                 Swarbrick_CID4290A_ER_panCK_cnv, Swarbrick_CID4398_ER_panCK_cnv, Swarbrick_CID4461_ER_panCK_cnv,
                                                                                 Swarbrick_CID4463_ER_panCK_cnv, Swarbrick_CID4471_ER_panCK_cnv, Swarbrick_CID4530N_ER_panCK_cnv,
                                                                                 
                                                                                 Swarbrick_CID4535_ER_panCK_cnv))
# Join Layers
Swarbrick_merged_ER_panCK_cnv <- JoinLayers(Swarbrick_merged_ER_panCK_cnv)


# Remove Objects
rm(Swarbrick_CID3941_ER_panCK_cnv)
rm(Swarbrick_CID3948_ER_panCK_cnv)
rm(Swarbrick_CID4040_ER_panCK_cnv)
rm(Swarbrick_CID4067_ER_panCK_cnv)
rm(Swarbrick_CID4290A_ER_panCK_cnv)
rm(Swarbrick_CID4398_ER_panCK_cnv)
rm(Swarbrick_CID4461_ER_panCK_cnv)
rm(Swarbrick_CID4463_ER_panCK_cnv)
rm(Swarbrick_CID4471_ER_panCK_cnv)
rm(Swarbrick_CID4530N_ER_panCK_cnv)
rm(Swarbrick_CID4535_ER_panCK_cnv)
gc()



####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples before merge

# Find Variable Features
Swarbrick_merged_ER_panCK_cnv <- FindVariableFeatures(Swarbrick_merged_ER_panCK_cnv, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Swarbrick_merged_ER_panCK_cnv)
Swarbrick_merged_ER_panCK_cnv <- ScaleData(Swarbrick_merged_ER_panCK_cnv, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################

Swarbrick_merged_ER_panCK_cnv <- RunPCA(Swarbrick_merged_ER_panCK_cnv, verbose = FALSE)

# Examine and visualize PCA results a few different ways
print(Swarbrick_merged_ER_panCK_cnv[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Swarbrick_merged_ER_panCK_cnv, dims = 1:2, reduction = "pca")
DimPlot(Swarbrick_merged_ER_panCK_cnv, reduction = "pca")

# Visualize PCs
ElbowPlot(Swarbrick_merged_ER_panCK_cnv)

# Choose dimensions
Swarbrick_merged_ER_panCK_cnv <- FindNeighbors(Swarbrick_merged_ER_panCK_cnv, dims = 1:20)
Swarbrick_merged_ER_panCK_cnv <- FindClusters(Swarbrick_merged_ER_panCK_cnv, resolution = 0.05)

# Umap clustering
Swarbrick_merged_ER_panCK_cnv <- RunUMAP(Swarbrick_merged_ER_panCK_cnv, dims = 1:20)
DimPlot(Swarbrick_merged_ER_panCK_cnv, reduction = "umap")

# Identify Samples
DimPlot(Swarbrick_merged_ER_panCK_cnv, reduction = "umap", group.by = "orig.ident")

# Save
saveRDS(Swarbrick_merged_ER_panCK_cnv, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Merged/Swarbrick_merged_ER_panCK_cnv.rds")
# Clear global environment
rm(Swarbrick_merged_ER_panCK_cnv)
gc()


#######################################
# Merge and Subset HER2 Tumor Samples #
#######################################

#########################
# Load Individual Files #
#########################
Swarbrick_CID3586_HER2_panCK_cnv <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_inferCNV/Swarbrick_CID3586_HER2_panCK_cnv.rds")
Swarbrick_CID3838_HER2_panCK_cnv <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_inferCNV/Swarbrick_CID3838_HER2_panCK_cnv.rds")
Swarbrick_CID3921_HER2_panCK_cnv <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_inferCNV/Swarbrick_CID3921_HER2_panCK_cnv.rds")
Swarbrick_CID4066_HER2_panCK_cnv <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_inferCNV/Swarbrick_CID4066_HER2_panCK_cnv.rds")
Swarbrick_CID45171_HER2_panCK_cnv <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_inferCNV/Swarbrick_CID45171_HER2_panCK_cnv.rds")


##############################################
# Merge Individual Files into Single Dataset #
##############################################
# Merge Datasets
Swarbrick_merged_HER2_panCK_CNV <- merge(x = Swarbrick_CID3586_HER2_panCK_cnv, y = c(Swarbrick_CID3838_HER2_panCK_cnv, Swarbrick_CID3921_HER2_panCK_cnv, Swarbrick_CID4066_HER2_panCK_cnv,
                                                                                     Swarbrick_CID45171_HER2_panCK_cnv))
# Join Layers
Swarbrick_merged_HER2_panCK_CNV <- JoinLayers(Swarbrick_merged_HER2_panCK_CNV)


# Remove Objects
rm(Swarbrick_CID3586_HER2_panCK_cnv)
rm(Swarbrick_CID3838_HER2_panCK_cnv)
rm(Swarbrick_CID3921_HER2_panCK_cnv)
rm(Swarbrick_CID4066_HER2_panCK_cnv)
rm(Swarbrick_CID45171_HER2_panCK_cnv)
gc()



####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples before merge

# Find Variable Features
Swarbrick_merged_HER2_panCK_CNV <- FindVariableFeatures(Swarbrick_merged_HER2_panCK_CNV, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Swarbrick_merged_HER2_panCK_CNV)
Swarbrick_merged_HER2_panCK_CNV <- ScaleData(Swarbrick_merged_HER2_panCK_CNV, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################

Swarbrick_merged_HER2_panCK_CNV <- RunPCA(Swarbrick_merged_HER2_panCK_CNV, verbose = FALSE)

# Examine and visualize PCA results a few different ways
print(Swarbrick_merged_HER2_panCK_CNV[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Swarbrick_merged_HER2_panCK_CNV, dims = 1:2, reduction = "pca")
DimPlot(Swarbrick_merged_HER2_panCK_CNV, reduction = "pca")

# Visualize PCs
ElbowPlot(Swarbrick_merged_HER2_panCK_CNV)

# Choose dimensions
Swarbrick_merged_HER2_panCK_CNV <- FindNeighbors(Swarbrick_merged_HER2_panCK_CNV, dims = 1:20)
Swarbrick_merged_HER2_panCK_CNV <- FindClusters(Swarbrick_merged_HER2_panCK_CNV, resolution = 0.05)

# Umap clustering
Swarbrick_merged_HER2_panCK_CNV <- RunUMAP(Swarbrick_merged_HER2_panCK_CNV, dims = 1:20)
DimPlot(Swarbrick_merged_HER2_panCK_CNV, reduction = "umap")

# Identify Samples
DimPlot(Swarbrick_merged_HER2_panCK_CNV, reduction = "umap", group.by = "orig.ident")

# Save
saveRDS(Swarbrick_merged_HER2_panCK_CNV, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Merged/Swarbrick_merged_HER2_panCK_cnv.rds")

# Clear global environment
rm(Swarbrick_merged_HER2_panCK_CNV)
gc()





#########################################################
############# This is TNBC Tumor Batch # ################
#########################################################

#########################
# Load Individual Files #
#########################
Swarbrick_CID3946_TNBC_panCK_cnv <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_inferCNV/Swarbrick_CID3946_TNBC_panCK_cnv.rds")
Swarbrick_CID3963_TNBC_panCK_cnv <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_inferCNV/Swarbrick_CID3963_TNBC_panCK_cnv.rds")
Swarbrick_CID4465_TNBC_panCK_cnv <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_inferCNV/Swarbrick_CID4465_TNBC_panCK_cnv.rds")
Swarbrick_CID4495_TNBC_panCK_cnv <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_inferCNV/Swarbrick_CID4495_TNBC_panCK_cnv.rds")
Swarbrick_CID4513_TNBC_panCK_cnv <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_inferCNV/Swarbrick_CID4513_TNBC_panCK_cnv.rds")
Swarbrick_CID4515_TNBC_panCK_cnv <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_inferCNV/Swarbrick_CID4515_TNBC_panCK_cnv.rds")
Swarbrick_CID4523_TNBC_panCK_cnv <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_inferCNV/Swarbrick_CID4523_TNBC_panCK_cnv.rds")
Swarbrick_CID44041_TNBC_panCK_cnv <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_inferCNV/Swarbrick_CID44041_TNBC_panCK_cnv.rds")
Swarbrick_CID44971_TNBC_panCK_cnv <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_inferCNV/Swarbrick_CID44971_TNBC_panCK_cnv.rds")
Swarbrick_CID44991_TNBC_panCK_cnv <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_inferCNV/Swarbrick_CID44991_TNBC_panCK_cnv.rds")

##############################################
# Merge Individual Files into Single Dataset #
##############################################
# Merge Datasets
Swarbrick_merged_TNBC_panCK_cnv <- merge(x = Swarbrick_CID3946_TNBC_panCK_cnv, y = c(Swarbrick_CID3963_TNBC_panCK_cnv, Swarbrick_CID4465_TNBC_panCK_cnv, Swarbrick_CID4495_TNBC_panCK_cnv,
                                                                                     Swarbrick_CID4513_TNBC_panCK_cnv, Swarbrick_CID4515_TNBC_panCK_cnv, Swarbrick_CID4523_TNBC_panCK_cnv,
                                                                                     Swarbrick_CID44041_TNBC_panCK_cnv, Swarbrick_CID44971_TNBC_panCK_cnv, Swarbrick_CID44991_TNBC_panCK_cnv))


# Join Layers
Swarbrick_merged_TNBC_panCK_cnv <- JoinLayers(Swarbrick_merged_TNBC_panCK_cnv)


# Remove Objects
rm(Swarbrick_CID3946_TNBC_panCK_cnv)
rm(Swarbrick_CID3963_TNBC_panCK_cnv)
rm(Swarbrick_CID4465_TNBC_panCK_cnv)
rm(Swarbrick_CID4495_TNBC_panCK_cnv)
rm(Swarbrick_CID4513_TNBC_panCK_cnv)
rm(Swarbrick_CID4515_TNBC_panCK_cnv)
rm(Swarbrick_CID4523_TNBC_panCK_cnv)
rm(Swarbrick_CID44041_TNBC_panCK_cnv)
rm(Swarbrick_CID44971_TNBC_panCK_cnv)
rm(Swarbrick_CID44991_TNBC_panCK_cnv)
gc()


####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples before merge

# Find Variable Features
Swarbrick_merged_TNBC_panCK_cnv <- FindVariableFeatures(Swarbrick_merged_TNBC_panCK_cnv, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(Swarbrick_merged_TNBC_panCK_cnv)
Swarbrick_merged_TNBC_panCK_cnv <- ScaleData(Swarbrick_merged_TNBC_panCK_cnv, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################

Swarbrick_merged_TNBC_panCK_cnv <- RunPCA(Swarbrick_merged_TNBC_panCK_cnv, verbose = FALSE)

# Examine and visualize PCA results a few different ways
print(Swarbrick_merged_TNBC_panCK_cnv[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Swarbrick_merged_TNBC_panCK_cnv, dims = 1:2, reduction = "pca")
DimPlot(Swarbrick_merged_TNBC_panCK_cnv, reduction = "pca")

# Visualize PCs
ElbowPlot(Swarbrick_merged_TNBC_panCK_cnv)

# Choose dimensions
Swarbrick_merged_TNBC_panCK_cnv <- FindNeighbors(Swarbrick_merged_TNBC_panCK_cnv, dims = 1:20)
Swarbrick_merged_TNBC_panCK_cnv <- FindClusters(Swarbrick_merged_TNBC_panCK_cnv, resolution = 0.05)

# Umap clustering
Swarbrick_merged_TNBC_panCK_cnv <- RunUMAP(Swarbrick_merged_TNBC_panCK_cnv, dims = 1:20)
DimPlot(Swarbrick_merged_TNBC_panCK_cnv, reduction = "umap")

# Identify Samples
DimPlot(Swarbrick_merged_TNBC_panCK_cnv, reduction = "umap", group.by = "orig.ident")

# Save
saveRDS(Swarbrick_merged_TNBC_panCK_cnv, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Merged/Swarbrick_merged_TNBC_panCK_cnv.rds")
# Clear global environment
rm(Swarbrick_merged_TNBC_panCK_cnv)
gc()


##################################################################
# Merge the panCK-subsetted files into one for combined analysis #
##################################################################

########
# LOAD #
########
Swarbrick_merged_ER_panCK_cnv <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Merged/Swarbrick_merged_ER_panCK_cnv.rds")
Swarbrick_merged_HER2_panCK_cnv <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Merged/Swarbrick_merged_HER2_panCK_cnv.rds")
Swarbrick_merged_TNBC_panCK_cnv <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Merged/Swarbrick_merged_TNBC_panCK_cnv.rds")


##############################################
# Merge Individual Files into Single Dataset #
##############################################
# Merge Datasets & remove intermediate files
all_Swarbrick_breast_tumors_panCK_cnv <- merge(x = Swarbrick_merged_ER_panCK_cnv, y = c(Swarbrick_merged_HER2_panCK_cnv, Swarbrick_merged_TNBC_panCK_cnv))

# Join Layers
all_Swarbrick_breast_tumors_panCK_cnv <- JoinLayers(all_Swarbrick_breast_tumors_panCK_cnv)

# Remove Objects
rm(Swarbrick_merged_ER_panCK_cnv)
rm(Swarbrick_merged_HER2_panCK_cnv)
rm(Swarbrick_merged_TNBC_panCK_cnv)


####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples before merge

# Find Variable Features
all_Swarbrick_breast_tumors_panCK_cnv <- FindVariableFeatures(all_Swarbrick_breast_tumors_panCK_cnv, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(all_Swarbrick_breast_tumors_panCK_cnv)
all_Swarbrick_breast_tumors_panCK_cnv <- ScaleData(all_Swarbrick_breast_tumors_panCK_cnv, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################

# Reduce dimensions using cell surface gene list
all_Swarbrick_breast_tumors_panCK_cnv <- RunPCA(all_Swarbrick_breast_tumors_panCK_cnv, verbose = FALSE)

# Examine and visualize PCA results a few different ways
print(all_Swarbrick_breast_tumors_panCK_cnv[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(all_Swarbrick_breast_tumors_panCK_cnv, dims = 1:2, reduction = "pca")
DimPlot(all_Swarbrick_breast_tumors_panCK_cnv, reduction = "pca")

# Visualize PCs
ElbowPlot(all_Swarbrick_breast_tumors_panCK_cnv, ndims = 40)

all_Swarbrick_breast_tumors_panCK_cnv <- FindNeighbors(all_Swarbrick_breast_tumors_panCK_cnv, dims = 1:37)
all_Swarbrick_breast_tumors_panCK_cnv <- FindClusters(all_Swarbrick_breast_tumors_panCK_cnv, resolution = 0.12)

# Umap clustering
all_Swarbrick_breast_tumors_panCK_cnv <- RunUMAP(all_Swarbrick_breast_tumors_panCK_cnv, dims = 1:37)
DimPlot(all_Swarbrick_breast_tumors_panCK_cnv, reduction = "umap", raster = FALSE, label = TRUE)
DimPlot(all_Swarbrick_breast_tumors_panCK_cnv, reduction = "umap", group.by = "orig.ident", raster = FALSE)


########
# SAVE #
########
saveRDS(all_Swarbrick_breast_tumors_panCK_cnv, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Merged/all_Swarbrick_breast_tumors_panCK_cnv.rds")


