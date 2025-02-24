#####################################################################################################################
#                       Swarbrick (Wu et al) Dataset Breast Tumor Analysis Steps                                    #
#####################################################################################################################
# Step 7: Determine Cluster Markers & Save Summary Information                                                      #
# Step 8: Shannon-Weaver Diversity Analysis                                                                         #
# Step 9: Performing Fgsea to Identify Biological Processes Enriched in Shared Clusters                             #
# Step 10: Annotating Diverse Clusters Based on Differentially Enriched Biological Processes                        #
# Step 11: Identifying Surface Receptors Enriched in Immune-like Clusters                                           #
# Step 12: Enumerating Immune Mimicry Surface Markers                                                               #
# Step 13: Export Immune Mimicry Marker Expression for Comparison Across Subgroups                                  #
# Step 14: Quantify Immune Mimicry Marker Positivity per Patient                                                    #
# Step 15: Compare DEGs and Signatures in Immune-like vs Non-Immune-like Neoplastic Cells                           #
# Step 16: Evaluate Cell Cycle Status of CD69-pos vs CD69-neg cells                                                 #
# Step 17: Assessing Other Features of Immune-Mimicked Cells                                                        #
# Step 18: Generate Expression Matrices for Each Sample: All Genes vs Immune Mimicry Markers                        #
#####################################################################################################################

################################################################################################
################################# Main Immune Mimicry Analysis #################################
################################################################################################

################################################################
# Step 7: Determine Cluster Markers & Save Summary Information #
################################################################
# Load Libraries
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)

########
# LOAD #
########
all_Swarbrick_breast_tumors_panCK_cnv <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Merged/all_Swarbrick_breast_tumors_panCK_cnv.rds")

# Perform Marker Analysis
all_Swarbrick_breast_tumors_panCK_cnv.cluster.markers <- FindAllMarkers(all_Swarbrick_breast_tumors_panCK_cnv, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(all_Swarbrick_breast_tumors_panCK_cnv.cluster.markers, file = "/R/R_Swarbrick/Swarbrick_Output/all_Swarbrick_breast_tumors_panCK_cnv.cluster.markers.csv")

# Save Example Plots
DimPlot(all_Swarbrick_breast_tumors_panCK_cnv, reduction = "umap", raster = FALSE)
ggsave("/R/R_Swarbrick/Swarbrick_Output/all_Swarbrick_breast_tumors_panCK_cnv_clusters_UMAP.tiff", plot = last_plot(), device = "tiff",
       scale = 1, width = 16, height = 10,
       dpi = 200, limitsize = TRUE)

DimPlot(all_Swarbrick_breast_tumors_panCK_cnv, reduction = "umap", raster = FALSE, label = TRUE)
ggsave("/R/R_Swarbrick/Swarbrick_Output/all_Swarbrick_breast_tumors_panCK_cnv_clusters_UMAP_labels.tiff", plot = last_plot(), device = "tiff",
       scale = 1, width = 16, height = 10,
       dpi = 200, limitsize = TRUE)

DimPlot(all_Swarbrick_breast_tumors_panCK_cnv, reduction = "umap", group.by = "orig.ident", raster = FALSE)
ggsave("/R/R_Swarbrick/Swarbrick_Output/all_Swarbrick_breast_tumors_panCK_cnv_clusters_orig.ident_UMAP.tiff", plot = last_plot(), device = "tiff",
       scale = 1, width = 16, height = 10,
       dpi = 200, limitsize = TRUE)

FeaturePlot(all_Swarbrick_breast_tumors_panCK_cnv, c("KRT14"),  raster = FALSE, cols = c("grey90", "firebrick"))
ggsave("/R/R_Swarbrick/Swarbrick_Output/all_Swarbrick_breast_tumors_Neoplastic_KRT14_UMAP_Custom_Colors.tiff", plot = last_plot(), device = "tiff", 
       scale = 1, width = 16, height = 10,
       dpi = 200, limitsize = TRUE)

FeaturePlot(all_Swarbrick_breast_tumors_panCK_cnv, c("KRT18"),  raster = FALSE, cols = c("grey90", "firebrick"))
ggsave("/R/R_Swarbrick/Swarbrick_Output/all_Swarbrick_breast_tumors_Neoplastic_KRT18_UMAP_Custom_Colors.tiff", plot = last_plot(), device = "tiff", 
       scale = 1, width = 16, height = 10,
       dpi = 200, limitsize = TRUE)

FeaturePlot(all_Swarbrick_breast_tumors_panCK_cnv, c("KRT19"),  raster = FALSE, cols = c("grey90", "firebrick"))
ggsave("/R/R_Swarbrick/Swarbrick_Output/all_Swarbrick_breast_tumors_Neoplastic_KRT19_UMAP_Custom_Colors.tiff", plot = last_plot(), device = "tiff", 
       scale = 1, width = 16, height = 10,
       dpi = 200, limitsize = TRUE)

# Extract  Meta Data
all_Swarbrick_breast_tumors_panCK_cnv.meta.data <- as.data.frame(as.matrix(all_Swarbrick_breast_tumors_panCK_cnv@meta.data))

# Save Meta Data
write.csv(all_Swarbrick_breast_tumors_panCK_cnv.meta.data, file = "/R/R_Swarbrick/Swarbrick_Output/all_Swarbrick_breast_tumors_panCK_cnv.meta.data.csv")

# Count Total Cells 
Total <- Idents(all_Swarbrick_breast_tumors_panCK_cnv)
# Total = 31,086 cells from 26 tumors



#############################################
# Step 8: Shannon-Weaver Diversity Analysis #            
#############################################

###################################################
# Calculate Number of Cells per Ident per Cluster #
# for downstream diversity calculations ###########
###################################################
# Store cluster identities in object@meta.data$cluster.designation
all_Swarbrick_breast_tumors_panCK_cnv[["cluster.designation"]] <- Idents(object = all_Swarbrick_breast_tumors_panCK_cnv)

# Get number of cells per cluster and per sample of origin
all_Swarbrick_breast_tumors_panCK_cnv.idents.per.cluster <- table(all_Swarbrick_breast_tumors_panCK_cnv@meta.data$cluster.designation, all_Swarbrick_breast_tumors_panCK_cnv@meta.data$orig.ident)

# Save table
write.csv(all_Swarbrick_breast_tumors_panCK_cnv.idents.per.cluster, file = "/R/R_Swarbrick/Swarbrick_Output/all_Swarbrick_breast_tumors_panCK_cnv.idents.per.cluster.csv")


###########################################
# Begin Shannon-Weaver Diversity Analysis #
###########################################
library(vegan)

# Load Table
ident.cluster.table <- read.csv(file = "/R/R_Swarbrick/Swarbrick_Output/all_Swarbrick_breast_tumors_panCK_cnv.idents.per.cluster.csv", header = TRUE, row.names = 1)

# Calculate Shannon index
shannon_diversity <- diversity(ident.cluster.table)  

# Save results
write.csv(shannon_diversity, file = "/R/R_Swarbrick/Swarbrick_Output/all_Swarbrick_breast_tumors_panCK_cnv.idents.per.cluster.shannon_diversity.csv")

# Remove objects
rm(all_Swarbrick_breast_tumors_panCK_cnv.idents.per.cluster)
rm(ident.cluster.table)
rm(shannon_diversity)
gc()



#########################################################################################
# Step 9: Performing Fgsea to Identify Biological Processes Enriched in Shared Clusters #
#########################################################################################
# Load Libraries
library(tidyverse)
library(data.table)
library(fgsea)

#########################################
# Most Diverse Clusters (Shannon Index) #
# SI Cutoff > 1.5                       #
# Cluster 5 = 2.52                      #
# Cluster 14 = 2.38                     #
# Cluster 17 = 1.93                     #
# Cluster 11 = 1.83                     #
# Cluster 9 = 1.82                      #
# Cluster 20 = 1.54                     #
# Cluster 6 = 1.50                      #
#########################################


# Load Cluster Markers
all_Swarbrick_breast_tumors_panCK_cnv.cluster.markers <- read.csv(file = "/R/R_Swarbrick/Swarbrick_Output/all_Swarbrick_breast_tumors_panCK_cnv.cluster.markers.csv")

# Remove genes with adjusted p-values above <0.01
all_Swarbrick_breast_tumors_panCK_cnv.cluster.markers <- all_Swarbrick_breast_tumors_panCK_cnv.cluster.markers %>% filter(p_val_adj<0.01)

# Collect columns: gene names, avg_log2FC, and cluster; sort by descending 
all_Swarbrick_breast_tumors_panCK_cnv.cluster.markers <- all_Swarbrick_breast_tumors_panCK_cnv.cluster.markers %>% dplyr::select(gene, avg_log2FC, cluster)
all_Swarbrick_breast_tumors_panCK_cnv.cluster.markers <- all_Swarbrick_breast_tumors_panCK_cnv.cluster.markers[order(all_Swarbrick_breast_tumors_panCK_cnv.cluster.markers$avg_log2FC, decreasing = TRUE),]  

# Subset highly diverse clusters & filter on gene name 
Cluster5 <- all_Swarbrick_breast_tumors_panCK_cnv.cluster.markers %>% filter(cluster == 5) 
Cluster14 <- all_Swarbrick_breast_tumors_panCK_cnv.cluster.markers %>% filter(cluster == 14) 
Cluster17 <- all_Swarbrick_breast_tumors_panCK_cnv.cluster.markers %>% filter(cluster == 17) 
Cluster11 <- all_Swarbrick_breast_tumors_panCK_cnv.cluster.markers %>% filter(cluster == 11) 
Cluster9 <- all_Swarbrick_breast_tumors_panCK_cnv.cluster.markers %>% filter(cluster == 9)
Cluster20 <- all_Swarbrick_breast_tumors_panCK_cnv.cluster.markers %>% filter(cluster == 20)
Cluster6 <- all_Swarbrick_breast_tumors_panCK_cnv.cluster.markers %>% filter(cluster == 6)

# Remove marker table
rm(all_Swarbrick_breast_tumors_panCK_cnv.cluster.markers)

# Save Output
write.csv(Cluster5, file = "/R/R_Swarbrick/Swarbrick_fgsea/Swarbrick_fgsea_Output/Cluster5.csv")
write.csv(Cluster14, file = "/R/R_Swarbrick/Swarbrick_fgsea/Swarbrick_fgsea_Output/Cluster14.csv")
write.csv(Cluster17, file = "/R/R_Swarbrick/Swarbrick_fgsea/Swarbrick_fgsea_Output/Cluster17.csv")
write.csv(Cluster11, file = "/R/R_Swarbrick/Swarbrick_fgsea/Swarbrick_fgsea_Output/Cluster11.csv")
write.csv(Cluster9, file = "/R/R_Swarbrick/Swarbrick_fgsea/Swarbrick_fgsea_Output/Cluster9.csv")
write.csv(Cluster20, file = "/R/R_Swarbrick/Swarbrick_fgsea/Swarbrick_fgsea_Output/Cluster20.csv")
write.csv(Cluster6, file = "/R/R_Swarbrick/Swarbrick_fgsea/Swarbrick_fgsea_Output/Cluster6.csv")

# Collapse into gene name and LogFC
Cluster5 <- Cluster5 %>% dplyr::select(gene, avg_log2FC)
Cluster14 <- Cluster14 %>% dplyr::select(gene, avg_log2FC)
Cluster17 <- Cluster17 %>% dplyr::select(gene, avg_log2FC)
Cluster11 <- Cluster11 %>% dplyr::select(gene, avg_log2FC)
Cluster9 <- Cluster9 %>% dplyr::select(gene, avg_log2FC)
Cluster20 <- Cluster20 %>% dplyr::select(gene, avg_log2FC)
Cluster6 <- Cluster6 %>% dplyr::select(gene, avg_log2FC)

# Create a vector of gene ranks
Cluster5.ranks <- deframe(Cluster5)
Cluster14.ranks <- deframe(Cluster14)
Cluster17.ranks <- deframe(Cluster17)
Cluster11.ranks <- deframe(Cluster11)
Cluster9.ranks <- deframe(Cluster9)
Cluster20.ranks <- deframe(Cluster20)
Cluster6.ranks <- deframe(Cluster6)

# Remove original objects
rm(Cluster5)
rm(Cluster14)
rm(Cluster17)
rm(Cluster11)
rm(Cluster9)
rm(Cluster20)
rm(Cluster6)

# Load the pathways into a named list
pathways.GO.BP <- gmtPathways("/R/R_Swarbrick/Swarbrick_fgsea/Swarbrick_fgsea_MSigDB/c5.go.bp.v2022.1.Hs.symbols.gmt")


##########################
# Run fgsea on Cluster5  #
##########################
fgsea.Cluster5.res <- fgsea(pathways = pathways.GO.BP, 
                            stats    = Cluster5.ranks,
                            eps      = 0.0,
                            minSize  = 15,
                            maxSize  = 500, nPermSimple = 100000)

# Summarize Cluster5.results by collapsing into main pathways, to omit redundancies
Cluster5.collapsedPathways <- collapsePathways(fgsea.Cluster5.res[order(pval)][padj < 0.05], 
                                               pathways.GO.BP, Cluster5.ranks)
Cluster5.mainPathways <- fgsea.Cluster5.res[pathway %in% Cluster5.collapsedPathways$mainPathways][
  order(-NES), pathway]
plotGseaTable(pathways.GO.BP[Cluster5.mainPathways], Cluster5.ranks, fgsea.Cluster5.res, 
              gseaParam = 0.5)

# Filter Cluster5.results so they only include the main pathways
fgsea.Cluster5.res_collapsed <- filter(fgsea.Cluster5.res, pathway %in% Cluster5.mainPathways)

# Save Cluster5.results
fwrite(fgsea.Cluster5.res, file="/R/R_Swarbrick/Swarbrick_fgsea/Swarbrick_fgsea_Output/fgsea.Cluster5.res_GO.BP.csv")
fwrite(fgsea.Cluster5.res_collapsed, file="/R/R_Swarbrick/Swarbrick_fgsea/Swarbrick_fgsea_Output/fgsea.Cluster5.res_GO.BP_collapsed.csv")

# Retain pathways with adjusted p-values below <0.05 
fgsea.Cluster5.res_sig <- fgsea.Cluster5.res %>% filter(padj<0.05)

# Save Cluster5.results
fwrite(fgsea.Cluster5.res_sig, file="/R/R_Swarbrick/Swarbrick_fgsea/Swarbrick_fgsea_Output/fgsea.Cluster5.res_sig.csv")


# Remove objects
rm(fgsea.Cluster5.res)
rm(fgsea.Cluster5.res_collapsed)
rm(Cluster5.collapsedPathways)
rm(Cluster5.mainPathways)
rm(Cluster5.ranks)
rm(fgsea.Cluster5.res_sig)


###########################
# Run fgsea on Cluster14  #
###########################
fgsea.Cluster14.res <- fgsea(pathways = pathways.GO.BP, 
                             stats    = Cluster14.ranks,
                             eps      = 0.0,
                             minSize  = 15,
                             maxSize  = 500, nPermSimple = 100000)

# Summarize Cluster14.results by collapsing into main pathways, to omit redundancies
Cluster14.collapsedPathways <- collapsePathways(fgsea.Cluster14.res[order(pval)][padj < 0.05], 
                                                pathways.GO.BP, Cluster14.ranks)
Cluster14.mainPathways <- fgsea.Cluster14.res[pathway %in% Cluster14.collapsedPathways$mainPathways][
  order(-NES), pathway]
plotGseaTable(pathways.GO.BP[Cluster14.mainPathways], Cluster14.ranks, fgsea.Cluster14.res, 
              gseaParam = 0.5)

# Filter Cluster14.results so they only include the main pathways
fgsea.Cluster14.res_collapsed <- filter(fgsea.Cluster14.res, pathway %in% Cluster14.mainPathways)

# Save Cluster14.results
fwrite(fgsea.Cluster14.res, file="/R/R_Swarbrick/Swarbrick_fgsea/Swarbrick_fgsea_Output/fgsea.Cluster14.res_GO.BP.csv")
fwrite(fgsea.Cluster14.res_collapsed, file="/R/R_Swarbrick/Swarbrick_fgsea/Swarbrick_fgsea_Output/fgsea.Cluster14.res_GO.BP_collapsed.csv")

# Retain pathways with adjusted p-values below <0.05 
fgsea.Cluster14.res_sig <- fgsea.Cluster14.res %>% filter(padj<0.05)

# Save Cluster14.results
fwrite(fgsea.Cluster14.res_sig, file="/R/R_Swarbrick/Swarbrick_fgsea/Swarbrick_fgsea_Output/fgsea.Cluster14.res_sig.csv")

# Remove objects
rm(fgsea.Cluster14.res)
rm(fgsea.Cluster14.res_collapsed)
rm(Cluster14.collapsedPathways)
rm(Cluster14.mainPathways)
rm(Cluster14.ranks)
rm(fgsea.Cluster14.res_sig)


##########################
# Run fgsea on Cluster17 #
##########################
fgsea.Cluster17.res <- fgsea(pathways = pathways.GO.BP, 
                             stats    = Cluster17.ranks,
                             eps      = 0.0,
                             minSize  = 15,
                             maxSize  = 500, nPermSimple = 100000)

# Summarize Cluster17.results by collapsing into main pathways, to omit redundancies
Cluster17.collapsedPathways <- collapsePathways(fgsea.Cluster17.res[order(pval)][padj < 0.05], 
                                                pathways.GO.BP, Cluster17.ranks)
Cluster17.mainPathways <- fgsea.Cluster17.res[pathway %in% Cluster17.collapsedPathways$mainPathways][
  order(-NES), pathway]
plotGseaTable(pathways.GO.BP[Cluster17.mainPathways], Cluster17.ranks, fgsea.Cluster17.res, 
              gseaParam = 0.5)

# Filter Cluster17.results so they only include the main pathways
fgsea.Cluster17.res_collapsed <- filter(fgsea.Cluster17.res, pathway %in% Cluster17.mainPathways)

# Save Cluster17.results
fwrite(fgsea.Cluster17.res, file="/R/R_Swarbrick/Swarbrick_fgsea/Swarbrick_fgsea_Output/fgsea.Cluster17.res_GO.BP.csv")
fwrite(fgsea.Cluster17.res_collapsed, file="/R/R_Swarbrick/Swarbrick_fgsea/Swarbrick_fgsea_Output/fgsea.Cluster17.res_GO.BP_collapsed.csv")

# Retain pathways with adjusted p-values below <0.05 
fgsea.Cluster17.res_sig <- fgsea.Cluster17.res %>% filter(padj<0.05)

# Save Cluster14.results
fwrite(fgsea.Cluster17.res_sig, file="/R/R_Swarbrick/Swarbrick_fgsea/Swarbrick_fgsea_Output/fgsea.Cluster17.res_sig.csv")

# Remove objects
rm(fgsea.Cluster17.res)
rm(fgsea.Cluster17.res_collapsed)
rm(Cluster17.collapsedPathways)
rm(Cluster17.mainPathways)
rm(Cluster17.ranks)
rm(fgsea.Cluster17.res_sig)


###########################
# Run fgsea on Cluster11  #
###########################
fgsea.Cluster11.res <- fgsea(pathways = pathways.GO.BP, 
                            stats    = Cluster11.ranks,
                            eps      = 0.0,
                            minSize  = 15,
                            maxSize  = 500, nPermSimple = 100000)

# Summarize Cluster11.results by collapsing into main pathways, to omit redundancies
Cluster11.collapsedPathways <- collapsePathways(fgsea.Cluster11.res[order(pval)][padj < 0.05], 
                                               pathways.GO.BP, Cluster11.ranks)
Cluster11.mainPathways <- fgsea.Cluster11.res[pathway %in% Cluster11.collapsedPathways$mainPathways][
  order(-NES), pathway]
plotGseaTable(pathways.GO.BP[Cluster11.mainPathways], Cluster11.ranks, fgsea.Cluster11.res, 
              gseaParam = 0.5)

# Filter Cluster11.results so they only include the main pathways
fgsea.Cluster11.res_collapsed <- filter(fgsea.Cluster11.res, pathway %in% Cluster11.mainPathways)

# Save Cluster11.results
fwrite(fgsea.Cluster11.res, file="/R/R_Swarbrick/Swarbrick_fgsea/Swarbrick_fgsea_Output/fgsea.Cluster11.res_GO.BP.csv")
fwrite(fgsea.Cluster11.res_collapsed, file="/R/R_Swarbrick/Swarbrick_fgsea/Swarbrick_fgsea_Output/fgsea.Cluster11.res_GO.BP_collapsed.csv")

# Retain pathways with adjusted p-values below <0.05 
fgsea.Cluster11.res_sig <- fgsea.Cluster11.res %>% filter(padj<0.05)

# Save Cluster14.results
fwrite(fgsea.Cluster11.res_sig, file="/R/R_Swarbrick/Swarbrick_fgsea/Swarbrick_fgsea_Output/fgsea.Cluster11.res_sig.csv")

# Remove objects
rm(fgsea.Cluster11.res)
rm(fgsea.Cluster11.res_collapsed)
rm(Cluster11.collapsedPathways)
rm(Cluster11.mainPathways)
rm(Cluster11.ranks)
rm(fgsea.Cluster11.res_sig)


##########################
# Run fgsea on Cluster9  #
##########################
fgsea.Cluster9.res <- fgsea(pathways = pathways.GO.BP, 
                             stats    = Cluster9.ranks,
                             eps      = 0.0,
                             minSize  = 15,
                             maxSize  = 500, nPermSimple = 100000)

# Summarize Cluster9.results by collapsing into main pathways, to omit redundancies
Cluster9.collapsedPathways <- collapsePathways(fgsea.Cluster9.res[order(pval)][padj < 0.05], 
                                                pathways.GO.BP, Cluster9.ranks)
Cluster9.mainPathways <- fgsea.Cluster9.res[pathway %in% Cluster9.collapsedPathways$mainPathways][
  order(-NES), pathway]
plotGseaTable(pathways.GO.BP[Cluster9.mainPathways], Cluster9.ranks, fgsea.Cluster9.res, 
              gseaParam = 0.5)

# Filter Cluster9.results so they only include the main pathways
fgsea.Cluster9.res_collapsed <- filter(fgsea.Cluster9.res, pathway %in% Cluster9.mainPathways)

# Save Cluster9.results
fwrite(fgsea.Cluster9.res, file="/R/R_Swarbrick/Swarbrick_fgsea/Swarbrick_fgsea_Output/fgsea.Cluster9.res_GO.BP.csv")
fwrite(fgsea.Cluster9.res_collapsed, file="/R/R_Swarbrick/Swarbrick_fgsea/Swarbrick_fgsea_Output/fgsea.Cluster9.res_GO.BP_collapsed.csv")

# Retain pathways with adjusted p-values below <0.05 
fgsea.Cluster9.res_sig <- fgsea.Cluster9.res %>% filter(padj<0.05)

# Save Cluster14.results
fwrite(fgsea.Cluster9.res_sig, file="/R/R_Swarbrick/Swarbrick_fgsea/Swarbrick_fgsea_Output/fgsea.Cluster9.res_sig.csv")

# Remove objects
rm(fgsea.Cluster9.res)
rm(fgsea.Cluster9.res_collapsed)
rm(Cluster9.collapsedPathways)
rm(Cluster9.mainPathways)
rm(Cluster9.ranks)
rm(fgsea.Cluster9.res_sig)


##########################
# Run fgsea on Cluster20 #
##########################
fgsea.Cluster20.res <- fgsea(pathways = pathways.GO.BP, 
                            stats    = Cluster20.ranks,
                            eps      = 0.0,
                            minSize  = 15,
                            maxSize  = 500, nPermSimple = 100000)

# Summarize Cluster20.results by collapsing into main pathways, to omit redundancies
Cluster20.collapsedPathways <- collapsePathways(fgsea.Cluster20.res[order(pval)][padj < 0.05], 
                                               pathways.GO.BP, Cluster20.ranks)
Cluster20.mainPathways <- fgsea.Cluster20.res[pathway %in% Cluster20.collapsedPathways$mainPathways][
  order(-NES), pathway]
plotGseaTable(pathways.GO.BP[Cluster20.mainPathways], Cluster20.ranks, fgsea.Cluster20.res, 
              gseaParam = 0.5)

# Filter Cluster20.results so they only include the main pathways
fgsea.Cluster20.res_collapsed <- filter(fgsea.Cluster20.res, pathway %in% Cluster20.mainPathways)

# Save Cluster20.results
fwrite(fgsea.Cluster20.res, file="/R/R_Swarbrick/Swarbrick_fgsea/Swarbrick_fgsea_Output/fgsea.Cluster20.res_GO.BP.csv")
fwrite(fgsea.Cluster20.res_collapsed, file="/R/R_Swarbrick/Swarbrick_fgsea/Swarbrick_fgsea_Output/fgsea.Cluster20.res_GO.BP_collapsed.csv")

# Retain pathways with adjusted p-values below <0.05 
fgsea.Cluster20.res_sig <- fgsea.Cluster20.res %>% filter(padj<0.05)

# Save Cluster14.results
fwrite(fgsea.Cluster20.res_sig, file="/R/R_Swarbrick/Swarbrick_fgsea/Swarbrick_fgsea_Output/fgsea.Cluster20.res_sig.csv")

# Remove objects
rm(fgsea.Cluster20.res)
rm(fgsea.Cluster20.res_collapsed)
rm(Cluster20.collapsedPathways)
rm(Cluster20.mainPathways)
rm(Cluster20.ranks)
rm(fgsea.Cluster20.res_sig)



##########################
# Run fgsea on Cluster6  #
##########################
fgsea.Cluster6.res <- fgsea(pathways = pathways.GO.BP, 
                             stats    = Cluster6.ranks,
                             eps      = 0.0,
                             minSize  = 15,
                             maxSize  = 500, nPermSimple = 100000)

# Summarize Cluster6.results by collapsing into main pathways, to omit redundancies
Cluster6.collapsedPathways <- collapsePathways(fgsea.Cluster6.res[order(pval)][padj < 0.05], 
                                                pathways.GO.BP, Cluster6.ranks)
Cluster6.mainPathways <- fgsea.Cluster6.res[pathway %in% Cluster6.collapsedPathways$mainPathways][
  order(-NES), pathway]
plotGseaTable(pathways.GO.BP[Cluster6.mainPathways], Cluster6.ranks, fgsea.Cluster6.res, 
              gseaParam = 0.5)

# Filter Cluster6.results so they only include the main pathways
fgsea.Cluster6.res_collapsed <- filter(fgsea.Cluster6.res, pathway %in% Cluster6.mainPathways)

# Save Cluster6.results
fwrite(fgsea.Cluster6.res, file="/R/R_Swarbrick/Swarbrick_fgsea/Swarbrick_fgsea_Output/fgsea.Cluster6.res_GO.BP.csv")
fwrite(fgsea.Cluster6.res_collapsed, file="/R/R_Swarbrick/Swarbrick_fgsea/Swarbrick_fgsea_Output/fgsea.Cluster6.res_GO.BP_collapsed.csv")

# Retain pathways with adjusted p-values below <0.05 
fgsea.Cluster6.res_sig <- fgsea.Cluster6.res %>% filter(padj<0.05)

# Save Cluster14.results
fwrite(fgsea.Cluster6.res_sig, file="/R/R_Swarbrick/Swarbrick_fgsea/Swarbrick_fgsea_Output/fgsea.Cluster6.res_sig.csv")

# Remove objects
rm(fgsea.Cluster6.res)
rm(fgsea.Cluster6.res_collapsed)
rm(Cluster6.collapsedPathways)
rm(Cluster6.mainPathways)
rm(Cluster6.ranks)
rm(fgsea.Cluster6.res_sig)

# Remove pathway list
rm(pathways.GO.BP)



##############################################################################################
# Step 10: Annotating Diverse Clusters Based on Differentially Enriched Biological Processes #
##############################################################################################

################################################
# Most Diverse Clusters (Shannon Index)        #
# SI Cutoff > 1.5                              #
# Cluster 5 = 2.52 = Immune (Lymphoid)         #
# Cluster 14 = 2.38 = Immune (Myeloid)         #
# Cluster 17 = 1.93 = Mesenchymal              #
# Cluster 11 = 1.83 = Developmental            #
# Cluster 9 = 1.82 = Endothelial               #
# Cluster 20 = 1.54 = Immune (Myeloid)        #
# Cluster 6 = 1.50 =  Catabolic (Not annotated)#
################################################

# Load Seurat Object
all_Swarbrick_breast_tumors_panCK_cnv <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Merged/all_Swarbrick_breast_tumors_panCK_cnv.rds")

# Visualize Clusters
DimPlot(all_Swarbrick_breast_tumors_panCK_cnv, reduction = "umap", raster = FALSE)
DimPlot(all_Swarbrick_breast_tumors_panCK_cnv, reduction = "umap", group.by = 'orig.ident', raster = FALSE)

# Annotate Clusters = detailed.ident
Mimicry_IDs_detailed <- c("Neoplastic Cells", "Neoplastic Cells", "Neoplastic Cells", "Neoplastic Cells", "Neoplastic Cells",
                          "Lymphoid-like", "Neoplastic Cells", "Neoplastic Cells", "Neoplastic Cells", "Endothelial-like",
                          "Neoplastic Cells", "Developmental-like", "Neoplastic Cells", "Neoplastic Cells", "Myeloid-like",
                          "Neoplastic Cells", "Neoplastic Cells", "Mesenchymal-like", "Neoplastic Cells", "Neoplastic Cells",
                          "Myeloid-like", "Neoplastic Cells", "Neoplastic Cells", "Neoplastic Cells")

# Apply New Labels to Clusters in Already Annotated Object
names(Mimicry_IDs_detailed) <- levels(all_Swarbrick_breast_tumors_panCK_cnv)
all_Swarbrick_breast_tumors_panCK_cnv_annotated <- RenameIdents(all_Swarbrick_breast_tumors_panCK_cnv, Mimicry_IDs_detailed)

# Remove old Seurat Object and ID list
rm(Mimicry_IDs_detailed)
rm(all_Swarbrick_breast_tumors_panCK_cnv)

# Store Annotations Under Classification
all_Swarbrick_breast_tumors_panCK_cnv_annotated[["detailed.ident"]] <- Idents(object = all_Swarbrick_breast_tumors_panCK_cnv_annotated)

# Draw plot with Changed colors
DimPlot(all_Swarbrick_breast_tumors_panCK_cnv_annotated, reduction = "umap", raster = FALSE, group.by = 'detailed.ident', cols = c('Neoplastic Cells' = '#00A9FF', 'Myeloid-like' = '#F8766D', 'Lymphoid-like' = '#CD9600', 'Endothelial-like' = '#C77CFF', 
                                                                                                                                   'Mesenchymal-like' = '#7CAE00', 'Developmental-like' = '#00BE67'))

# Save Plot
ggsave("/R/R_Swarbrick/Swarbrick_Output/all_Swarbrick_breast_tumors_panCK_cnv_detailed_annotation_with_custom_colors.tiff", plot = last_plot(), device = "tiff", 
       scale = 1, width = 16, height = 10,
       dpi = 200, limitsize = TRUE)

# Check Idents
Idents(object = all_Swarbrick_breast_tumors_panCK_cnv_annotated)
table(Idents(all_Swarbrick_breast_tumors_panCK_cnv_annotated))



# Annotate Clusters = main.ident
Mimicry_IDs_main <- c("Neoplastic Cells", "Immune-like", "Endothelial-like", "Developmental-like", "Immune-like", "Mesenchymal-like")

# Apply New Labels to Clusters
names(Mimicry_IDs_main) <- levels(all_Swarbrick_breast_tumors_panCK_cnv_annotated)
all_Swarbrick_breast_tumors_panCK_cnv_annotated <- RenameIdents(all_Swarbrick_breast_tumors_panCK_cnv_annotated, Mimicry_IDs_main)

# Remove old Seurat Object and ID list
rm(Mimicry_IDs_main)

# Store Annotations Under Classification
all_Swarbrick_breast_tumors_panCK_cnv_annotated[["main.ident"]] <- Idents(object = all_Swarbrick_breast_tumors_panCK_cnv_annotated)

# Draw plot with Changed colors
DimPlot(all_Swarbrick_breast_tumors_panCK_cnv_annotated, reduction = "umap", raster = FALSE, group.by = 'main.ident', cols = c('Neoplastic Cells' = '#00A9FF', 'Immune-like' = '#F8766D', 'Endothelial-like' = '#C77CFF', 
                                                                                                                               'Mesenchymal-like' = '#7CAE00', 'Developmental-like' = '#00BE67'))

# Save Plot
ggsave("/R/R_Swarbrick/Swarbrick_Output/all_Swarbrick_breast_tumors_panCK_cnv_main_annotation_with_custom_colors.tiff", plot = last_plot(), device = "tiff", 
       scale = 1, width = 16, height = 10,
       dpi = 200, limitsize = TRUE)

# Check Idents
table(Idents(all_Swarbrick_breast_tumors_panCK_cnv_annotated))

# Set Identity for moving forward
all_Swarbrick_breast_tumors_panCK_cnv_annotated <- SetIdent(all_Swarbrick_breast_tumors_panCK_cnv_annotated, value = "detailed.ident")

# Verify Idents have changed
table(Idents(all_Swarbrick_breast_tumors_panCK_cnv_annotated))

# Save Annotated Seurat Object
saveRDS(all_Swarbrick_breast_tumors_panCK_cnv_annotated, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Annotated/all_Swarbrick_breast_tumors_panCK_cnv_annotated.rds")

# Extract  Meta Data
all_Swarbrick_breast_tumors_panCK_cnv_annotated.meta.data <- as.data.frame(as.matrix(all_Swarbrick_breast_tumors_panCK_cnv_annotated@meta.data))

# Save Meta Data
write.csv(all_Swarbrick_breast_tumors_panCK_cnv_annotated.meta.data, file = "/R/R_Swarbrick/Swarbrick_Output/all_Swarbrick_breast_tumors_panCK_cnv_annotated.meta.data.csv")

# Remove Meta Data
rm(all_Swarbrick_breast_tumors_panCK_cnv_annotated.meta.data)

#########################################
# Draw IM Gene Expression Feature Plots #
#########################################
#all_Swarbrick_breast_tumors_panCK_cnv_annotated <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Annotated/all_Swarbrick_breast_tumors_panCK_cnv_annotated.rds")

# Use patchwork to combine plots and incorporate custom colors

p1 <- FeaturePlot(all_Swarbrick_breast_tumors_panCK_cnv_annotated, c("CD3D"), cols = c("grey90", "darkgoldenrod4"))
p2 <- FeaturePlot(all_Swarbrick_breast_tumors_panCK_cnv_annotated, c("CD52"), cols = c("grey90", "darkgoldenrod4"))
p3 <- FeaturePlot(all_Swarbrick_breast_tumors_panCK_cnv_annotated, c("CD69"), cols = c("grey90", "darkgoldenrod4"))
p4 <- FeaturePlot(all_Swarbrick_breast_tumors_panCK_cnv_annotated, c("CD14"), cols = c("grey90", "dodgerblue4"))
p5 <- FeaturePlot(all_Swarbrick_breast_tumors_panCK_cnv_annotated, c("CD68"), cols = c("grey90", "dodgerblue4"))
p6 <- FeaturePlot(all_Swarbrick_breast_tumors_panCK_cnv_annotated, c("CD74"), cols = c("grey90", "dodgerblue4"))
p7 <- FeaturePlot(all_Swarbrick_breast_tumors_panCK_cnv_annotated, c("PTPRC"), cols = c("grey90", "darkgreen"))
p8 <- FeaturePlot(all_Swarbrick_breast_tumors_panCK_cnv_annotated, c("CD44"), cols = c("grey90", "darkgreen"))
p9 <- FeaturePlot(all_Swarbrick_breast_tumors_panCK_cnv_annotated, c("CXCR4"), cols = c("grey90", "darkgreen"))

# Combine plots with patchwork
(p1 | p2 | p3) / (p4 | p5 | p6) / (p7 | p8 | p9)

ggsave("/R/R_Swarbrick/Swarbrick_Output/all_Swarbrick_breast_tumors_panCK_cnv_annotated_IM_Marker_Patchwork.tiff", plot = last_plot(), device = "tiff", 
       scale = 1, width = 16, height = 10,
       dpi = 200, limitsize = TRUE)

# Remove Plots & Object
rm(p1)
rm(p2)
rm(p3)
rm(p4)
rm(p5)
rm(p6)
rm(p7)
rm(p8)
rm(p9)



###########################################################################
# Step 11: Identifying Surface Receptors Enriched in Immune-like Clusters #
###########################################################################

######################################################
# Surface Receptors Enriched in Immune-like Clusters #
######################################################
# Surface receptors associated with the immune mimicry phenotype are determined below
# List of CD receptor genes were obtained from the Uniprot database:
# https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/docs/cdlist.txt
# Gene symbols from the above table collected into a csv, with "gene" as its header; this being done outside R
# Below I extract CD DEGs specifically for the immune-mimicked clusters, and then filter them based on significance and Log2FC


###########################################
# Most Diverse Clusters (Shannon Index)   #
# SI Cutoff > 1.5                         #
# Cluster 5 = 2.52 = Immune (Lymphoid)    #
# Cluster 14 = 2.38 = Immune (Myeloid)    #
# Cluster 17 = 1.93 = Mesenchymal         #
# Cluster 11 = 1.83 = Developmental       #
# Cluster 9 = 1.82 = Endothelial          #
# Cluster 20 = 1.54 = Immune (Myeloid)    #
# Cluster 6 = 1.50 =  Catabolic           #
###########################################


########
# LOAD #
########
# Read cluster markers
all_Swarbrick_breast_tumors_panCK_cnv <- read.csv("/R/R_Swarbrick/Swarbrick_Output/all_Swarbrick_breast_tumors_panCK_cnv.cluster.markers.csv")

# Read surface protein list 
uniprot_cd_gene_list <- read.csv(file = "/R/R_Swarbrick/Swarbrick_Gene_Lists/Uniprot_CD_Gene_List_2024.csv", header = TRUE)

# Specify gene name
# Note column 8 is the gene name
colnames(all_Swarbrick_breast_tumors_panCK_cnv)[8] <- "Gene"
# Filter Immune-mimicked cells: Lymphoid = Cluster 5, Myeloid = Cluster 14, Myeloid = Cluster 20
all_Swarbrick_breast_tumors_panCK_cnv_filtered <- all_Swarbrick_breast_tumors_panCK_cnv[all_Swarbrick_breast_tumors_panCK_cnv$cluster == 5 | all_Swarbrick_breast_tumors_panCK_cnv$cluster == 14 | all_Swarbrick_breast_tumors_panCK_cnv$cluster == 20, ]

# Remove unfiltered cluster markers
rm(all_Swarbrick_breast_tumors_panCK_cnv)

# Change cluster column to mimicry status
all_Swarbrick_breast_tumors_panCK_cnv_filtered <- all_Swarbrick_breast_tumors_panCK_cnv_filtered %>% 
  mutate(cluster = ifelse(cluster == 5, "lymphoid",
                          ifelse(cluster == 13, "myeloid",
                                 ifelse(cluster == 20, "Myeloid-2", cluster))))
colnames(all_Swarbrick_breast_tumors_panCK_cnv_filtered)[7] <- "mimicry status"

# Filter DEGs based on Surface Protein List
all_Swarbrick_breast_tumors_panCK_cnv_filtered <- all_Swarbrick_breast_tumors_panCK_cnv_filtered[all_Swarbrick_breast_tumors_panCK_cnv_filtered$Gene %in% uniprot_cd_gene_list$gene,]

# Filter CD DEGs to only significant genes with fold-change > 1.5, e.g. (log2FC(1.5) = 0.58)
Swarbrick_scRNAseq_Immune_Mimicry_CD_Markers_FC1.5_Cutoff <- all_Swarbrick_breast_tumors_panCK_cnv_filtered[all_Swarbrick_breast_tumors_panCK_cnv_filtered$p_val_adj <= 0.001 & all_Swarbrick_breast_tumors_panCK_cnv_filtered$avg_log2FC >= 0.58, ] 

# Save Results
write.csv(Swarbrick_scRNAseq_Immune_Mimicry_CD_Markers_FC1.5_Cutoff, file = "/R/R_Swarbrick/Swarbrick_Output/Swarbrick_scRNAseq_Immune_Mimicry_CD_Markers_FC1.5_Cutoff.csv")

# Remove objects
rm(all_Swarbrick_breast_tumors_panCK_cnv_filtered)
rm(Swarbrick_scRNAseq_Immune_Mimicry_CD_Markers_FC1.5_Cutoff)
rm(uniprot_cd_gene_list)



#######################################################
# Step 12: Enumerating Immune Mimicry Surface Markers #
#######################################################

# Major Immune Mimicry Surface Markers
# Selecting Surface Markers Consistent in Both Visvader & Swarbrick Datasets
#FCGR3A
#CD2
#CD3E
#CD7
#CD3G
#CD3D
#IL7R
#FCGR2A
#CD68
#CD52
#CD69
#CD14
#PTPRC
#TNFSF13B
#CD37
#ITGB2
#CXCR4
#CD74
#CLEC7A
#CD53
#CD83
#PLAUR
#CD44

# Found in Swarbrick but not Visvader
# Many others not listed here for brevity


###############
# Load Object #
###############
all_Swarbrick_breast_tumors_panCK_cnv <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Merged/all_Swarbrick_breast_tumors_panCK_cnv.rds")

# Cells Expressing Marker
panCK_count <- length(WhichCells(all_Swarbrick_breast_tumors_panCK_cnv, expression = KRT14 > 1 | KRT18 > 1 | KRT19 > 1 ))
FCGR3A_count <- length(WhichCells(all_Swarbrick_breast_tumors_panCK_cnv, expression = FCGR3A > 0))
CD2_count <- length(WhichCells(all_Swarbrick_breast_tumors_panCK_cnv, expression = CD2 > 0))
CD3E_count <- length(WhichCells(all_Swarbrick_breast_tumors_panCK_cnv, expression = CD3E > 0))
CD7_count <- length(WhichCells(all_Swarbrick_breast_tumors_panCK_cnv, expression = CD7 > 0))
CD3G_count <- length(WhichCells(all_Swarbrick_breast_tumors_panCK_cnv, expression = CD3G > 0))
CD3D_count <- length(WhichCells(all_Swarbrick_breast_tumors_panCK_cnv, expression = CD3D > 0))
IL7R_count <- length(WhichCells(all_Swarbrick_breast_tumors_panCK_cnv, expression = IL7R > 0))
FCGR2A_count <- length(WhichCells(all_Swarbrick_breast_tumors_panCK_cnv, expression = FCGR2A > 0))
CD68_count <- length(WhichCells(all_Swarbrick_breast_tumors_panCK_cnv, expression = CD68 > 0))
CD52_count <- length(WhichCells(all_Swarbrick_breast_tumors_panCK_cnv, expression = CD52 > 0))
CD69_count <- length(WhichCells(all_Swarbrick_breast_tumors_panCK_cnv, expression = CD69 > 0))
CD14_count <- length(WhichCells(all_Swarbrick_breast_tumors_panCK_cnv, expression = CD14 > 0))
PTPRC_count <- length(WhichCells(all_Swarbrick_breast_tumors_panCK_cnv, expression = PTPRC > 0))
TNFSF13B_count <- length(WhichCells(all_Swarbrick_breast_tumors_panCK_cnv, expression = TNFSF13B > 0))
CD37_count <- length(WhichCells(all_Swarbrick_breast_tumors_panCK_cnv, expression = CD37 > 0))
ITGB2_count <- length(WhichCells(all_Swarbrick_breast_tumors_panCK_cnv, expression = ITGB2 > 0))
CXCR4_count <- length(WhichCells(all_Swarbrick_breast_tumors_panCK_cnv, expression = CXCR4 > 0))
CD74_count <- length(WhichCells(all_Swarbrick_breast_tumors_panCK_cnv, expression = CD74 > 0))
CLEC7A_count <- length(WhichCells(all_Swarbrick_breast_tumors_panCK_cnv, expression = CLEC7A > 0))
CD53_count <- length(WhichCells(all_Swarbrick_breast_tumors_panCK_cnv, expression = CD53 > 0))
CD83_count <- length(WhichCells(all_Swarbrick_breast_tumors_panCK_cnv, expression = CD83 > 0))
PLAUR_count <- length(WhichCells(all_Swarbrick_breast_tumors_panCK_cnv, expression = PLAUR > 0))
CD44_count <- length(WhichCells(all_Swarbrick_breast_tumors_panCK_cnv, expression = CD44 > 0))

#######################################
# Compile Counts & Save as Data Frame #
#######################################

Immune_Mimicry_count <- c(panCK_count,
                          FCGR3A_count,
                          CD2_count,
                          CD3E_count,
                          CD7_count,
                          CD3G_count,
                          CD3D_count,
                          IL7R_count,
                          FCGR2A_count,
                          CD68_count,
                          CD52_count,
                          CD69_count,
                          CD14_count,
                          PTPRC_count,
                          TNFSF13B_count,
                          CD37_count,
                          ITGB2_count,
                          CXCR4_count,
                          CD74_count,
                          CLEC7A_count,
                          CD53_count,
                          CD83_count,
                          PLAUR_count,
                          CD44_count)

# Convert to data frame
Immune_Mimicry_count <- as.data.frame(Immune_Mimicry_count)

# Generate list for row names
Immune_Mimicry_count_rn <- c("panCK_count",
                             "FCGR3A_count",
                             "CD2_count",
                             "CD3E_count",
                             "CD7_count",
                             "CD3G_count",
                             "CD3D_count",
                             "IL7R_count",
                             "FCGR2A_count",
                             "CD68_count",
                             "CD52_count",
                             "CD69_count",
                             "CD14_count",
                             "PTPRC_count",
                             "TNFSF13B_count",
                             "CD37_count",
                             "ITGB2_count",
                             "CXCR4_count",
                             "CD74_count",
                             "CLEC7A_count",
                             "CD53_count",
                             "CD83_count",
                             "PLAUR_count",
                             "CD44_count")

# Add row names to data frame
rownames(Immune_Mimicry_count) <- Immune_Mimicry_count_rn 

# Save results
write.csv(Immune_Mimicry_count, file = "/R/R_Swarbrick/Swarbrick_Output/Swarbrick_Immune_Mimicry_count.csv")

# Remove Objects
rm(panCK_count)
rm(FCGR3A_count)
rm(CD2_count)
rm(CD3E_count)
rm(CD7_count)
rm(CD3G_count)
rm(CD3D_count)
rm(IL7R_count)
rm(FCGR2A_count)
rm(CD68_count)
rm(CD52_count)
rm(CD69_count)
rm(CD14_count)
rm(PTPRC_count)
rm(TNFSF13B_count)
rm(CD37_count)
rm(ITGB2_count)
rm(CXCR4_count)
rm(CD74_count)
rm(CLEC7A_count)
rm(CD53_count)
rm(CD83_count)
rm(PLAUR_count)
rm(CD44_count)
rm(Immune_Mimicry_count_rn)
rm(Immune_Mimicry_count)
gc()



####################################################################################
# Step 13: Export Immune Mimicry Marker Expression for Comparison Across Subgroups #
####################################################################################

# Export Averaged and Scale Expression Matrices for Mimicry Markers
all_Swarbrick_breast_tumors_panCK_cnv_annotated <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Annotated/all_Swarbrick_breast_tumors_panCK_cnv_annotated.rds")
# Set Ident
table(Idents(all_Swarbrick_breast_tumors_panCK_cnv_annotated))
all_Swarbrick_breast_tumors_panCK_cnv_annotated <- SetIdent(all_Swarbrick_breast_tumors_panCK_cnv_annotated, value = "detailed.ident")
table(Idents(all_Swarbrick_breast_tumors_panCK_cnv_annotated))
# Subset mimicry subtypes from combined object
Myeloid.like <- subset(x = all_Swarbrick_breast_tumors_panCK_cnv_annotated, idents = c("Myeloid-like"))
Lymphoid.like <- subset(x = all_Swarbrick_breast_tumors_panCK_cnv_annotated, idents = c("Lymphoid-like"))
Mimicry.neg <- subset(x = all_Swarbrick_breast_tumors_panCK_cnv_annotated, idents = c("Neoplastic Cells"))

# Extract and save IM expression matrix: Myeloid-like
Myeloid.like <- Myeloid.like[["RNA"]]$scale.data
Myeloid.like <- as.matrix(Myeloid.like, 'sparseMatrix')
Myeloid.like <- subset(Myeloid.like, rownames(Myeloid.like) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Myeloid.like <- as.data.frame(rowMeans(Myeloid.like))
# Further organize and save data frame
Myeloid.like <-t(Myeloid.like)
write.csv(Myeloid.like, file = "/R/R_Swarbrick/Swarbrick_Expression_Matrices/Swarbrick_Expression_Matrices_Output/Swarbrick_Myeloid.like_marker_scaled_expression.csv")
rm(Myeloid.like)
gc()

# Extract and save IM expression matrix: Lymphoid-like
Lymphoid.like <- Lymphoid.like[["RNA"]]$scale.data
Lymphoid.like <- as.matrix(Lymphoid.like, 'sparseMatrix')
Lymphoid.like <- subset(Lymphoid.like, rownames(Lymphoid.like) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Lymphoid.like <- as.data.frame(rowMeans(Lymphoid.like))
# Further organize and save data frame
Lymphoid.like <-t(Lymphoid.like)
write.csv(Lymphoid.like, file = "/R/R_Swarbrick/Swarbrick_Expression_Matrices/Swarbrick_Expression_Matrices_Output/Swarbrick_Lymphoid.like_marker_scaled_expression.csv")
rm(Lymphoid.like)
gc()

# Extract and save IM expression matrix: Mimicry.neg
Mimicry.neg <- Mimicry.neg[["RNA"]]$scale.data
Mimicry.neg <- as.matrix(Mimicry.neg, 'sparseMatrix')
Mimicry.neg <- subset(Mimicry.neg, rownames(Mimicry.neg) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Mimicry.neg <- as.data.frame(rowMeans(Mimicry.neg))
# Further organize and save data frame
Mimicry.neg <-t(Mimicry.neg)
write.csv(Mimicry.neg, file = "/R/R_Swarbrick/Swarbrick_Expression_Matrices/Swarbrick_Expression_Matrices_Output/Swarbrick_Mimicry.neg_marker_scaled_expression.csv")
rm(Mimicry.neg)
gc()



##################################################################
# Step 14: Quantify Immune Mimicry Marker Positivity per Patient #
##################################################################
# Load Object
all_Swarbrick_breast_tumors_panCK_cnv_annotated <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Annotated/all_Swarbrick_breast_tumors_panCK_cnv_annotated.rds")


# Extract the gene expression matrix
expression_matrix <- GetAssayData(all_Swarbrick_breast_tumors_panCK_cnv_annotated, slot = "data")
# Extract the metadata to identify samples
metadata <- all_Swarbrick_breast_tumors_panCK_cnv_annotated@meta.data
# Initialize a list to store counts for each gene
gene_counts_list <- list()

# Define a function to process a gene of interest
process_gene <- function(gene_symbol, expression_matrix, metadata) {
  if (gene_symbol %in% rownames(expression_matrix)) {
    # Identify cells expressing the gene (expression > 0)
    gene_positive_cells <- expression_matrix[gene_symbol, ] > 0
    
    # Combine the gene expression data with metadata
    metadata[[paste0(gene_symbol, "_positive")]] <- gene_positive_cells
    
    # Count the number of positive cells per sample
    gene_counts_per_sample <- aggregate(gene_positive_cells ~ orig.ident, data = metadata, sum)
    names(gene_counts_per_sample)[2] <- gene_symbol  # Rename the column to gene symbol
    
    # Store gene counts in the list
    gene_counts_list[[gene_symbol]] <<- gene_counts_per_sample
    
    # Print the counts
    print(gene_counts_per_sample)
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
combined_counts <- gene_counts_list[[genes_of_interest[1]]]  # Start with the first gene
for (gene in genes_of_interest[-1]) {
  combined_counts <- merge(combined_counts, gene_counts_list[[gene]], by = "orig.ident", all = TRUE)
}

# Sort columns alphabetically by gene symbol
combined_counts <- combined_counts[, c("orig.ident", sort(genes_of_interest))]

# Save the combined counts to a single CSV file
write.csv(combined_counts, file = "/R/R_Swarbrick/Swarbrick_Output/all_Swarbrick_breast_tumors_panCK_cnv_annotated_IMBC_Genes_per_Sample.csv", row.names = FALSE)


###################
# Enumerate panCK #
###################
genes_of_interest <- c("KRT14", "KRT18", "KRT19")
if (all(genes_of_interest %in% rownames(expression_matrix))) {
  # Identify cells expressing KRT14 > 1, KRT18 > 1, or KRT19 > 1
  KRT14_positive_cells <- expression_matrix["KRT14", ] > 1
  KRT18_positive_cells <- expression_matrix["KRT18", ] > 1
  KRT19_positive_cells <- expression_matrix["KRT19", ] > 1
  
  # Combine the logical vectors using OR condition
  positive_cells <- KRT14_positive_cells | KRT18_positive_cells | KRT19_positive_cells
  
  
  # Combine the result with metadata
  metadata$positive_cells <- positive_cells
  
  # Count the number of positive cells per sample
  positive_counts_per_sample <- aggregate(positive_cells ~ orig.ident, data = metadata, sum)
  
  # Save the counts as a data frame
  # You can change the file path and name as needed
  write.csv(positive_counts_per_sample, file = "/R/R_Swarbrick/Swarbrick_Output/all_Swarbrick_breast_tumors_panCK_cnv_annotated_panCK_per_Sample.csv", row.names = FALSE)
  
  
  # Print the counts
  print(positive_counts_per_sample)
} else {
  print("CD69 gene is not present in the dataset.")
}

# Remove objects
rm(expression_matrix) 
rm(metadata)
rm(gene_counts_list)
rm(process_gene)
rm(genes_of_interest)
rm(combined_counts)
gc()



###########################################################################################
# Step 15: Compare DEGs and Signatures in Immune-like vs Non-Immune-like Neoplastic Cells #
###########################################################################################
# all_Swarbrick_breast_tumors_panCK_cnv_annotated <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Annotated/all_Swarbrick_breast_tumors_panCK_cnv_annotated.rds")

# Check Idents
table(Idents(all_Swarbrick_breast_tumors_panCK_cnv_annotated))

# Set Identity
all_Swarbrick_breast_tumors_panCK_cnv_annotated <- SetIdent(all_Swarbrick_breast_tumors_panCK_cnv_annotated, value = "detailed.ident")

# Check Idents
table(Idents(all_Swarbrick_breast_tumors_panCK_cnv_annotated))

# Perform differential gene expression analysis for Lymphoid immune mimicry
Immune_Mimicry.Lymphoid.markers <- FindMarkers(all_Swarbrick_breast_tumors_panCK_cnv_annotated, ident.1 = "Lymphoid-like", ident.2 = "Neoplastic Cells")
write.csv(Immune_Mimicry.Lymphoid.markers, file = "/R/R_Swarbrick/Swarbrick_Output/Immune_Mimicry.Lymphoid.markers.csv")
rm(Immune_Mimicry.Lymphoid.markers)

# Perform differential gene expression analysis for Myeloid immune mimicry
Immune_Mimicry.Myeloid.markers <- FindMarkers(all_Swarbrick_breast_tumors_panCK_cnv_annotated, ident.1 = "Myeloid-like", ident.2 = "Neoplastic Cells")
write.csv(Immune_Mimicry.Myeloid.markers, file = "/R/R_Swarbrick/Swarbrick_Output/Immune_Mimicry.Myeloid.markers.csv")
rm(Immune_Mimicry.Myeloid.markers)

# Set Next Identity
all_Swarbrick_breast_tumors_panCK_cnv_annotated <- SetIdent(all_Swarbrick_breast_tumors_panCK_cnv_annotated, value = "main.ident")

# Verify Idents have changed
table(Idents(all_Swarbrick_breast_tumors_panCK_cnv_annotated))


# Perform differential gene expression analysis for immune mimicry
Immune_Mimicry.de.markers <- FindMarkers(all_Swarbrick_breast_tumors_panCK_cnv_annotated, ident.1 = "Immune-like", ident.2 = "Neoplastic Cells")
write.csv(Immune_Mimicry.de.markers, file = "/R/R_Swarbrick/Swarbrick_Output/Immune_Mimicry.de.markers.csv")

# Remove Markers
rm(Immune_Mimicry.de.markers)


#########################################################
### Mammary Stem Cell & Functional Signature Analysis ###
#########################################################
# Signatures obtained from msigdb
# Published in: https://pubmed.ncbi.nlm.nih.gov/20346151/
# Lim_Mammary_Stem: http://www.gsea-msigdb.org/gsea/msigdb/human/geneset/LIM_MAMMARY_STEM_CELL_UP.html?keywords=lim%20mammary
# Lim_Luminal_Progenitor: http://www.gsea-msigdb.org/gsea/msigdb/human/geneset/LIM_MAMMARY_LUMINAL_PROGENITOR_UP.html?keywords=lim%20mammary
# Lim_Mature: http://www.gsea-msigdb.org/gsea/msigdb/human/geneset/LIM_MAMMARY_LUMINAL_MATURE_UP.html?keywords=lim%20mammary
# BIOCARTA_NFKB_PATHWAY: https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/BIOCARTA_NFKB_PATHWAY.html
# POSITIVE_REGULATION_OF_LEUKOCYTE_PROLIFERATION: https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_PROLIFERATION.html
# POSITIVE_REGULATION_OF_MAPK_CASCADE: https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/GOBP_POSITIVE_REGULATION_OF_MAPK_CASCADE.html
# POSITIVE_REGULATION_OF_CYTOKINE_PRODUCTION: https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/POSITIVE_REGULATION_OF_CYTOKINE_PRODUCTION.html
# POSITIVE_REGULATION_OF_LYMPHOCYTE_MIGRATION: https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/POSITIVE_REGULATION_OF_LYMPHOCYTE_MIGRATION.html
# POSITIVE_REGULATION_OF_MAPK_CASCADE: https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/GOBP_POSITIVE_REGULATION_OF_MAPK_CASCADE.html

# Compile Signatures
# A few genes were outdated and these were updated manually to current symbols (e.g. BLAM1 -> ARNTL)
Lim_Mammary_Stem <- list(c("ABI3BP","ABTB3","ACER2","ACTA2","ACTG2","ACVR2A","ADAMTS1","ADAMTS2","ADARB1","ADGRA2","ADGRL1","AEBP1","AGPAT4","AHI1","AKT3","ALDH1L2","AMOTL1","ANGPTL2","ANKRD1","ANTXR1","AOPEP","APOE","AQP9","ARC","ARHGAP20","ARHGAP24","ARHGAP25","ARHGEF25","ARHGEF28","ARMH4","ARSI","ASPHD2","ATP2A2","AXL","BACH1","BAG3","BCAR1","BCL2L11","BCOR","ARNTL","BMP1","BMP7","BNC1","BVES","C19orf12","C1QTNF12","C1QTNF4","CACNB4","CADM1","CALD1","CALML3","CALU","CARD10","CAV1","CAVIN2","CBLB","CCDC3","CCDC85B","CTGF","CCND2","CD36","CD70","CDC42EP2","CDH13","CDH3","CDH4","CDKN1A","CHST3","CHST7","CLIP3","CLMP","CLXN","CNN1","CNP","CNRIP1","COL12A1","COL14A1","COL16A1","COL17A1","COL18A1","COL23A1","COL4A1","COL4A2","COL5A1","COL5A2","COL7A1","COL9A2","CPNE8","CPXM1","CRISPLD1","CRLF1","CRYAB","CSDC2","CSPG4","CSRP2","CTNNAL1","CXCL14","CYGB","DCBLD2","DCHS1","DCUN1D3","DIPK2A","DKK3","DLK2","DLL1","DMWD","DOCK10","DPYSL3","DST","DUSP6","DUSP7","DZIP1L","EBF3","ECRG4","EDARADD","EDNRB","EEPD1","EFNB1","EGFR","EGR2","EGR3","EID3","ELK3","ELOVL4","ELP5","ENC1","ENPP2","EOGT","EPAS1","EPDR1","EPHB1","ERF","ETS1","EVA1A","EVC","EXT1","FABP5","FAM184A","FAM216A","FAS","FBLN7","FBXO30","FERMT2","FEZ1","FGFRL1","FGL2","FHL1","FHOD3","FJX1","FLNC","FLRT2","FMOD","FOXP1","FST","FXYD1","FZD8","GEM","GJA1","GJC1","GNAI1","GNB4","GNG11","GOLIM4","GPC3","GPR176","GPR3","GPR87","GPSM1","GSN","GYPC","HACD1","HAS2","HDAC4","HEG1","HGFAC","HRAS","HS3ST3A1","HSPB2","HSPG2","HTRA1","ICAM1","ID4","IGFBP2","IGFBP3","IGFBP4","IGFBP6","IL17B","IL17RD","IL1B","IL24","IL6","IL6ST","INKA1","MRVI1","IRX4","ISM1","ITGA1","ITGA6","ITGA9","ITGB1","ITGB4","ITM2A","JAG1","JAM2","JAM3","KANK4","KCNIP3","KCNMA1","KCNMB1","KLHL21","KLHL29","KLHL42","KRT14","KRT16","KRT5","KRT75","L3HYPDH","LAG3","LAMA1","LAMA3","LAMB1","LAMB3","LAMC1","LARGE2","LBH","LCA5","LCAT","LEP","LGALS1","LGALS7","LGR6","LHFPL6","LIFR","LIMA1","LIMS2","LMOD1","LRCH2","LRP1","LRP4","LRRC8C","LRRN1","LTBP4","LUZP1","MALT1","MAMDC2","MAOB","MAP3K7CL","MATN2","MBNL1","MCAM","MEDAG","MEF2C","MEG3","MEST","MFNG","MGARP","MIA","MICAL2","MME","MMP2","MPDZ","MRGPRF","MSRB3","MSX1","MTSS1","MXRA7","MYC","MYH11","MYL9","MYLK","MYOCD","NBL1","NDN","NECTIN3","NETO2","NGF","NGFR","NLGN2","NNMT","NPTX2","NRCAM","NRG1","NRP1","NRP2","NSG1","NT5E","NTF3","NTRK2","NUDT11","NXN","OSBPL6","OSR1","OXTR","P3H1","P3H2","PALM2AKAP2","PAMR1","PARD6G","PCBP4","PCDH18","PCDH19","PCDH7","PCDHGC3","PCOLCE","PDGFA","PDLIM4","PDLIM7","PDPN","PEG3","PGF","PHLDA3","PHLDB1","PKD1","PKNOX2","PKP1","PLA2G7","PLCH2","PLEKHA4","PLPP1","PLPP3","PLS3","PLXNA2","PODN","POGLUT2","POPDC2","POSTN","POU3F1","PPP1R14A","PPP1R16B","PPP1R18","PPP1R3C","PPP2R2B","PRDM1","PRICKLE1","PRICKLE2","PRNP","PROS1","PRRX1","PRX","PSD2","PTGS2","PTPRE","PTPRT","PXDC1","PXN","QKI","QRICH2","RAB34","RAPGEF1","RARB","RARRES2","RASIP1","RASL12","RBPMS","RCN3","RCSD1","RECK","RELN","RFLNB","RFX2","RHOJ","RND3","RNF165","RUSC2","SBSPON","SCARF2","SCHIP1","SCML2","SCN4B","SDK2","SEC24D","SEMA3C","SEMA5A","SERPINF1","SERPING1","SERPINH1","SGCB","SGIP1","SH2D5","SH3TC1","SHE","SIAH2","SIMC1","SKI","SLC12A4","SLC1A3","SLC1A5","SLC25A4","SLC27A3","SLC27A6","SLC2A3","SLC38A5","SLC4A3","SLC6A8","SLCO3A1","SLIT2","SLIT3","SMIM13","SMTN","SNAI2","SNCA","SNTB2","SOBP","SOGA1","SORBS1","SORCS1","SOX11","SPHK1","SPRED1","SRGN","SRPX","SSBP2","SSH1","STAC","STAC2","STARD8","STXBP4","SULF1","SVEP1","SYDE1","SYNM","TACC1","TAGLN","TAMALIN","TBX2","TCF4","TCF7L1","TCOF1","TENM3","TES","TGFB1I1","TGFBR3","THBS1","THSD1","THY1","TIE1","TIMP3","TINAGL1","TM7SF3","TMEM121","TMEM178A","TMEM201","TMEM204","TMEM255B","TMEM47","TMEM64","TNS1","TNS4","TOX","TP63","TPM2","TPST1","TRIM29","TRIM9","TRO","TRPC1","TSHZ2","TSHZ3","TSKU","TSPY26P","TSPYL2","TTYH2","TUBB6","TWIST2","UCN2","UNC45A","UPP1","VCAN","VGLL3","VIM","VIT","VSIR","VSNL1","WIF1","WIPF1","WTIP","YAF2","ZC3H12B","ZNF219","ZNF423"))
Lim_Luminal_Progenitor <- list(c("ACSL1","ALDH1A3","ANPEP","ASIC1","ATP6V1B1","ATP6V1C2","BBOX1","C10orf90","C1QTNF1","C3","CCDC88B","CD14","CKMT1B","CLDN1","CSN2","CSN3","CTSC","CXCR4","CYP24A1","DAPP1","ELF5","FOLR1","FOXI1","GALNT15","GGT5","GJB2","GNE","HAPLN3","HIVEP3","HSD17B12","IL15","IL4I1","ITPR2","KIT","LALBA","LBP","LPCAT1","MELTF","NCALD","NOXO1","PDZK1IP1","PIGR","PLB1","QPCT","RASAL1","RASGEF1C","RFTN2","RPS6KL1","S100A8","SECTM1","SLC13A2","SLC28A3","SLC34A2","SORBS2","TNFAIP2","TSPAN33","WFDC3","XDH"))
Lim_Mature <- list(c("ABCA7","ABCC8","ACOT11","ALCAM","ALDH3B1","ALDH3B2","ANKMY2","AQP11","ARFGEF3","ATP6V0E2","BATF","BBOF1","BTRC","C1orf210","CACNB3","CACNG4","CASZ1","CCDC92","CITED1","DNAAF3","DNAJC12","DRC3","DUSP10","EDEM1","EEF1A2","EPS8L1","ERN1","ESR1","FAAH","FBXO36","FER1L4","FGF13","FGL1","FLVCR2","FOXA1","FYCO1","G6PD","GADD45G","GALE","GMPR","GPRC5C","H4C12","H4C8","HDAC11","HES6","HID1","HMGCS2","HOXB2","HSD11B2","IL13RA1","KBTBD4","KLHL5","LAMA5","LMNTD2","LNX2","MBOAT1","MEIS1","MEIS3","MINDY1","MLPH","MYB","NECTIN4","NSD3","PAK4","PGAP6","PGR","PHKA1","PLEKHG3","PON3","PRLR","PROM2","PSD4","PTPN6","PVALB","PXYLP1","RABL3","RASEF","REEP6","SCMH1","SGMS1","SLC16A5","SLC22A18","SLC40A1","SLC44A4","SLC7A2","SLC7A4","SORT1","SPDEF","SPINK1","SPRR1A","SULT2B1","TANGO2","TBX3","TGM2","TMCO3","TMPRSS6","TNFSF11","TOX3","TP53INP2","TRIM6","TSPAN1","TSPAN13","TUBG1","VOPP1","VPS33B","WNK4","WNT4","WNT5A","WNT7B","YIPF6","ZDHHC1","ZFHX2","ZSCAN18"))
BIOCARTA_NFKB_PATHWAY <- list(c("CHUK","FADD","IKBKB","IKBKG","IL1A","IL1R1","MAP3K1","MAP3K14","MAP3K7","MYD88","NFKB1","NFKBIA","RELA","RIPK1","TAB1","TNF","TNFAIP3","TNFRSF1A","TNFRSF1B","TRADD","TRAF6"))
POSITIVE_REGULATION_OF_LEUKOCYTE_PROLIFERATION <- list(c("ADA","AGER","AIF1","ANXA1","ATAD5","BCL2","BCL2L1","BCL6","BMI1","BST1","BST2","BTK","CARD11","CCDC88B","CCL19","CCL5","CCR2","CD1D","CD209","CD24","CD274","CD276","CD28","CD320","CD38","CD3E","CD40","CD40LG","CD46","CD55","CD6","CD70","CD74","CD80","CD81","CD86","CDKN1A","CHRNB2","CLCF1","CLECL1P","CORO1A","CSF1","CSF1R","CSF2","CSF2RA","CSF2RB","DHPS","DNAJA3","EBI3","EFNB1","EPHB2","EPO","FADD","FCGR3A","FCRL3","FGF10","FOXP3","GPAM","GPR183","HAVCR2","HES1","HHLA2","HLA-A","HLA-DMB","HLA-DPA1","HLA-DPB1","HLA-E","HMGB1","ICOSLG","IGF1","IGF2","IGFBP2","IL12A","IL12B","IL12RB1","IL13","IL15","IL18","IL1A","IL1B","IL2","IL21","IL23A","IL23R","IL2RA","IL34","IL4","IL5","IL5RA","IL6","IL6ST","IL7","IRS2","JAK2","JAK3","KIT","KITLG","LEP","LGALS9","LILRB2","LYN","MAPK1","MAPK3","MEF2C","MIF","MIR181B1","MIR21","MIR30B","MPL","MYD88","NCK1","NCK2","NCKAP1L","NFATC2","NMB","NMBR","OCSTAMP","PDCD1LG2","PELI1","PNP","PPP3CA","PRKCQ","PRLR","PTH","PTK2","PTPN22","PTPRC","PYCARD","RAC2","RASAL3","RIPK2","RPS3","SASH3","SELENOK","SHH","SLAMF1","SLC39A10","SLC7A1","SPTA1","ST6GAL1","STAT5B","SYK","TACR1","TFRC","TGFBR2","TICAM1","TIRAP","TLR4","TLR9","TMIGD2","TNFRSF13C","TNFRSF4","TNFSF13B","TNFSF4","TNFSF9","TRAF6","TYK2","VAV3","VCAM1","VTCN1","WNT3A","XCL1","ZAP70","ZNF335","ZP3","ZP4"))
POSITIVE_REGULATION_OF_MAPK_CASCADE <- list(c("ABCA7","ABL1","ACKR3","ACTA2","ADAM8","ADAM9","ADCYAP1","ADORA1","ADRA1A","ADRA1B","ADRA2A","ADRA2B","ADRA2C","ADRB2","ADRB3","AGER","AJUBA","AKAP12","AKAP13","ALKAL1","ALKAL2","ALOX12B","ALOX15","ANGPT1","ANKRD6","APELA","APOE","APP","AR","ARHGAP8","ARHGEF5","ARL6IP5","ARRB1","ARRB2","AVPI1","AVPR1B","AXIN1","BANK1","BCAR3","BIRC7","BMP2","BMP4","BMPER","BRAF","C1QTNF1","C5AR1","CALCR","CARD9","CARTPT","CASR","CAV2","CAVIN3","CCL1","CCL11","CCL13","CCL14","CCL15","CCL16","CCL17","CCL18","CCL19","CCL2","CCL20","CCL21","CCL22","CCL23","CCL24","CCL25","CCL26","CCL3","CCL3L1","CCL3L3","CCL4","CCL5","CCL7","CCL8","CCN2","CCR1","CCR7","CD24","CD27","CD36","CD4","CD40","CD44","CD74","CD81","CDC42","CDH2","CDK10","CDON","CFLAR","CHI3L1","CHRNA7","CIB1","CRK","CRKL","CSF1R","CSK","CSPG4","CTNNB1","CX3CL1","CXCL17","DAB2IP","DDR1","DDR2","DDT","DENND2B","DHX33","DIRAS1","DIRAS2","DIXDC1","DKK1","DNAJC27","DOK1","DOK2","DOK3","DOK4","DOK5","DOK6","DRD2","DRD4","DSTYK","DUSP19","DUSP22","DVL2","DVL3","EDA2R","EDAR","EDN1","EDN3","EFNA1","EGF","EGFR","EIF2AK2","ELANE","EPGN","EPHA4","EPHA8","EPO","ERBB2","ERBB4","ERN1","ERN2","ERP29","EZH2","F2R","F2RL1","FBXW7","FCGR2B","FCRL3","FERMT2","FFAR4","FGA","FGB","FGD2","FGF1","FGF10","FGF18","FGF19","FGF2","FGF20","FGF21","FGF23","FGF4","FGF8","FGFR1","FGFR2","FGFR3","FGFR4","FGG","FLT1","FLT3","FLT4","FPR2","FRS2","FSHR","FZD10","FZD4","FZD5","FZD7","FZD8","GADD45A","GADD45B","GADD45G","GAREM1","GAS6","GATA4","GCG","GCNT2","GDF15","GDF6","GFRAL","GH1","GHR","GHRL","GLIPR2","GNAI2","GPBAR1","GPER1","GPNMB","GPR183","GPR37","GPR37L1","GPR55","GRM1","GRM4","GRM5","GSDME","HAND2","HAVCR2","HCRTR1","HGF","HIPK2","HLA-DRB1","HMGB1","HRAS","HTR2A","HTR2B","HTR2C","IAPP","ICAM1","IGF1","IGF1R","IGF2","IGFBP3","IGFBP4","IGFBP6","IL11","IL1A","IL1B","IL26","IL34","IL6","INAVA","INHBA","INS","INSR","IQGAP1","IQGAP3","IRAK1","ITGA1","ITGB3","JAK2","JCAD","JUN","KDR","KISS1","KIT","KITLG","KL","KLB","KLHDC10","KSR1","LAMTOR1","LAMTOR2","LAMTOR3","LAPTM5","LEP","LGALS9","LIF","LILRA5","LPAR1","LPAR2","LPAR3","LRRK2","LTBR","MADD","MAGED1","MAP2K1","MAP2K2","MAP2K3","MAP2K4","MAP2K5","MAP2K6","MAP2K7","MAP3K10","MAP3K11","MAP3K12","MAP3K13","MAP3K3","MAP3K4","MAP3K5","MAP3K7","MAP4K1","MAP4K2","MAPK3","MAPK8IP1","MAPK8IP2","MAPK8IP3","MAPKBP1","MARCO","MBIP","MEF2C","MFAP3","MFHAS1","MID1","MIF","MINK1","MIR126","MIR181A2","MIR181B1","MIR181D","MIR21","MIR221","MIR222","MIR23A","MIR24-1","MIR27A","MIR27B","MIR519D","MIR92A1","MIRLET7B","MMP8","MOS","MST1R","MT3","MTURN","MUSK","MYD88","MYDGF","NAIP","NCF1","NDRG4","NDST1","NECAB2","NEK10","NELFE","NENF","NOD1","NOD2","NODAL","NOTCH1","NOTCH2","NOX1","NOX4","NPNT","NPSR1","NPTN","NPY","NPY5R","NRG1","NRP1","NTF3","NTRK1","NTRK2","NTRK3","OPRK1","OPRM1","OR2AT4","OSM","P2RX7","P2RY1","P2RY6","PAK1","PDCD10","PDE5A","PDE6G","PDE6H","PDE8A","PDGFA","PDGFB","PDGFC","PDGFD","PDGFRA","PDGFRB","PELI2","PHB1","PHB2","PIK3CG","PIK3R5","PIK3R6","PJA2","PLA2G1B","PLA2G2A","PLA2G5","PLCB1","PLCE1","PLCG2","PPIA","PRDX2","PRKCA","PRKCE","PRKCZ","PRKD2","PRMT1","PROK1","PRXL2C","PSEN1","PTK2B","PTPN1","PTPN11","PTPN22","PTPRC","PTPRJ","PYCARD","RAF1","RAMP3","RAP1A","RAP1B","RAPGEF2","RASGRP1","RASSF2","RB1CC1","RELL1","RELL2","RET","RIPK1","RIPK2","RIT2","ROBO1","ROCK1","ROCK2","ROR1","ROR2","RPS3","RYK","S100A12","S100A7","SASH1","SCIMP","SDCBP","SEMA3A","SEMA4C","SEMA7A","SERPINF2","SH3RF1","SH3RF2","SH3RF3","SHC1","SLAMF1","SLC30A10","SOD1","SORBS3","SOX2","SPAG9","SPHK1","SPI1","SPRY2","SRC","SSTR4","STK25","STK3","STK39","SYK","SYT14P1","TAB1","TAOK1","TAOK2","TAOK3","TBX1","TDGF1","TEK","TENM1","TGFA","TGFB1","TGFB2","TGFB3","TGFBR1","THBS1","THPO","TIRAP","TLR3","TLR4","TLR6","TLR9","TMEM106A","TNF","TNFAIP8L3","TNFRSF11A","TNFRSF19","TNFSF11","TNIK","TP73","TPBG","TPD52L1","TRAF1","TRAF2","TRAF3","TRAF4","TRAF5","TRAF6","TRAF7","TREM2","TRIM5","TRPV4","UNC5CL","VEGFA","WNT16","WNT5A","WNT7A","WNT7B","WWC1","XCL1","XCL2","XDH","XIAP","ZC3H12A","ZNF622"))
POSITIVE_REGULATION_OF_CYTOKINE_PRODUCTION <- list(c("ABCC8","ABL1","ADAM17","ADAM8","ADCYAP1","ADIPOQ","ADORA2B","ADRA2A","AFAP1L2","AGER","AGPAT1","AGPAT2","AGT","AIF1","AIM2","AIRE","AKAP12","AKIRIN2","ALOX15B","ANXA1","APOA2","APP","APPL1","ARFGEF2","ARHGEF2","ARID5A","ARNT","ATF2","ATF4","ATP6AP2","AZU1","B2M","BATF","BCL10","BCL3","BMPR1A","BRCA1","BTK","BTN3A1","BTN3A2","C1QTNF3","C1QTNF4","C3","C3AR1","C5","C5AR1","CADM1","CAMK4","CARD11","CARD8","CARD9","CASP1","CASP8","CCBE1","CCDC88B","CCL1","CCL19","CCL3","CCR2","CCR7","CD14","CD160","CD2","CD200","CD226","CD244","CD274","CD276","CD28","CD34","CD36","CD3E","CD4","CD40","CD40LG","CD46","CD55","CD58","CD6","CD74","CD80","CD81","CD83","CD84","CD86","CEACAM20","CEBPB","CEBPG","CGAS","CHI3L1","CHIA","CHUK","CLEC4E","CLEC5A","CLEC6A","CLEC7A","CLEC9A","CLECL1P","CLNK","CLU","CRLF2","CRTAM","CSF1R","CSF2","CX3CL1","CXCL17","CYBA","CYBB","CYP1B1","CYRIB","DDIT3","DDT","DDX1","DDX21","DDX3X","DEFA5","DEFB124","DENND1B","DHX33","DHX36","DHX58","DHX9","DRD2","EBI3","EGR1","EIF2AK2","EIF2AK3","ELANE","EPHB2","EPX","EREG","F2R","F2RL1","F3","FADD","FCER1G","FCGR3A","FCN1","FERMT1","FFAR2","FFAR3","FGR","FLOT1","FLT4","FOXP1","FOXP3","FRMD8","FURIN","FZD5","G3BP1","GAPDH","GARIN5A","GATA3","GATA4","GATA6","GBP5","GDF2","GLMN","GPRC5B","GPSM3","GSDMD","H19","HAVCR2","HDAC2","HEG1","HGF","HHLA2","HIF1A","HILPDA","HK1","HLA-A","HLA-DPA1","HLA-DPB1","HLA-E","HLA-F","HLA-G","HMGB1","HMGB2","HMHB1","HMOX1","HPSE","HRAS","HSP90AA1","HSPA1A","HSPA1B","HSPB1","HSPD1","HTR2B","HYAL2","IDO1","IFI16","IFIH1","IFNG","IFNGR1","IFNL1","IGHD","IL10","IL12A","IL12B","IL12RB1","IL12RB2","IL13","IL15","IL16","IL17A","IL17B","IL17D","IL17F","IL17RA","IL17RC","IL18","IL18R1","IL1A","IL1B","IL1R1","IL1RAP","IL1RL1","IL1RL2","IL2","IL20RB","IL21","IL23A","IL23R","IL26","IL27","IL27RA","IL33","IL36A","IL4","IL4R","IL6","IL6R","IL6ST","IL7","IL9","INAVA","INS","IRAK1","IRAK3","IRF1","IRF3","IRF4","IRF5","IRF7","IRF8","ISG15","ISL1","ITK","JAK2","KAT2A","KIR2DL4","KIT","KLRC4-KLRK1","KLRF2","KLRK1","KPNA6","LACC1","LAMTOR5","LAPTM5","LBP","LEP","LGALS9","LILRA2","LILRA5","LILRB1","LILRB2","LPL","LRRK2","LTA","LTB","LUM","LURAP1","LY9","LY96","MALT1","MAP3K7","MAPK11","MAPK13","MAPK14","MAPKAPK2","MAVS","MBP","MCOLN2","MDK","MEFV","MIF","MIR132","MIR144","MIR145","MIR149","MIR17","MIR182","MIR206","MIR21","MIR27B","MIR324","MIR657","MIR675","MIR92A1","MMP12","MMP8","MNDA","MYD88","NAIP","NFAM1","NFATC4","NLRC4","NLRP1","NLRP10","NLRP12","NLRP2","NLRP3","NLRP9","NMB","NMBR","NOD1","NOD2","NODAL","NOS2","NOX1","NOX5","NR1H4","NR4A3","OAS1","OAS2","OAS3","ORM1","ORM2","OSM","P2RX7","PAEP","PANX1","PANX2","PANX3","PARK7","PDE4B","PDE4D","PELI1","PF4","PHB1","PIBF1","PIK3CD","PIK3CG","PIK3R1","PLA2G1B","PLA2G3","PLA2R1","PLCB1","PLCG2","PNP","POLR3A","POLR3B","POLR3C","POLR3D","POLR3F","POLR3G","POU2AF1","POU2F2","PQBP1","PRG2","PRG3","PRKCQ","PRKCZ","PRKD2","PSEN1","PTAFR","PTGER4","PTGS2","PTPN11","PTPN22","PTPRC","PTPRJ","PYCARD","PYDC1","PYHIN1","RAB1A","RAB2B","RAB7B","RAET1G","RARA","RASGRP1","RBM47","RELA","RFTN1","RGCC","RIGI","RIOK3","RIPK1","RIPK2","RNF135","ROCK2","RORA","RPS3","RSAD2","RUNX1","S100A13","SAA1","SASH3","SCAMP5","SCIMP","SCRIB","SELENOK","SEMA7A","SERPINB7","SERPINE1","SERPINF2","SETD2","SETD4","SIGLEC16","SIRT1","SLAMF1","SLAMF6","SLC11A1","SLC7A5","SMAD3","SOD1","SORL1","SPHK1","SPHK2","SPN","SPON2","SPTBN1","SRC","STAT1","STAT3","STAT5B","STING1","STMP1","STOML2","SULF1","SULF2","SYK","TBK1","TBX21","TGFB1","THBS1","TICAM1","TICAM2","TIGIT","TIRAP","TLR1","TLR2","TLR3","TLR4","TLR5","TLR6","TLR7","TLR8","TLR9","TMED10","TMEM106A","TMF1","TMIGD2","TNF","TNFRSF14","TNFRSF8","TNFSF4","TNXB","TOMM70","TRAF2","TRAF6","TREM2","TRIM15","TRIM16","TRIM27","TRIM32","TRIM56","TRIM6","TRPV4","TSLP","TUSC2","TWIST1","TXK","TYK2","TYROBP","UCN","UNC93B1","USP50","VTCN1","WNT11","WNT3A","WNT5A","XBP1","XCL1","ZBTB20","ZBTB7B","ZCCHC3","ZFPM1","ZNF580","ZP3"))
POSITIVE_REGULATION_OF_LYMPHOCYTE_MIGRATION <- list(c("ABL1","ABL2","ADAM10","ADAM17","ADAM8","AIF1","APP","CCL20","CCL21","CCL3","CCL4","CCL5","CCL7","CCR2","CD99L2","CORO1A","CXCL10","CXCL12","CXCL13","DOCK8","FADD","ITGA4","ITGB3","JAM2","MADCAM1","NEDD9","OXSR1","PTK2B","PYCARD","RHOA","S100A7","SELENOK","SPN","STK39","TMEM102","TNFRSF14","TNFSF14","WNK1","WNT5A","XCL1","XCL2"))
POSITIVE_REGULATION_OF_MAPK_CASCADE <- list(c("ABCA7","ABL1","ACKR3","ACTA2","ADAM8","ADAM9","ADCYAP1","ADORA1","ADRA1A","ADRA1B","ADRA2A","ADRA2B","ADRA2C","ADRB2","ADRB3","AGER","AJUBA","AKAP12","AKAP13","ALKAL1","ALKAL2","ALOX12B","ALOX15","ANGPT1","ANKRD6","APELA","APOE","APP","AR","ARHGAP8","ARHGEF5","ARL6IP5","ARRB1","ARRB2","AVPI1","AVPR1B","AXIN1","BANK1","BCAR3","BIRC7","BMP2","BMP4","BMPER","BRAF","C1QTNF1","C5AR1","CALCR","CARD9","CARTPT","CASR","CAV2","CAVIN3","CCL1","CCL11","CCL13","CCL14","CCL15","CCL16","CCL17","CCL18","CCL19","CCL2","CCL20","CCL21","CCL22","CCL23","CCL24","CCL25","CCL26","CCL3","CCL3L1","CCL3L3","CCL4","CCL5","CCL7","CCL8","CCN2","CCR1","CCR7","CD24","CD27","CD36","CD4","CD40","CD44","CD74","CD81","CDC42","CDH2","CDK10","CDON","CFLAR","CHI3L1","CHRNA7","CIB1","CRK","CRKL","CSF1R","CSK","CSPG4","CTNNB1","CX3CL1","CXCL17","DAB2IP","DDR1","DDR2","DDT","DENND2B","DHX33","DIRAS1","DIRAS2","DIXDC1","DKK1","DNAJC27","DOK1","DOK2","DOK3","DOK4","DOK5","DOK6","DRD2","DRD4","DSTYK","DUSP19","DUSP22","DVL2","DVL3","EDA2R","EDAR","EDN1","EDN3","EFNA1","EGF","EGFR","EIF2AK2","ELANE","EPGN","EPHA4","EPHA8","EPO","ERBB2","ERBB4","ERN1","ERN2","ERP29","EZH2","F2R","F2RL1","FBXW7","FCGR2B","FCRL3","FERMT2","FFAR4","FGA","FGB","FGD2","FGF1","FGF10","FGF18","FGF19","FGF2","FGF20","FGF21","FGF23","FGF4","FGF8","FGFR1","FGFR2","FGFR3","FGFR4","FGG","FLT1","FLT3","FLT4","FPR2","FRS2","FSHR","FZD10","FZD4","FZD5","FZD7","FZD8","GADD45A","GADD45B","GADD45G","GAREM1","GAS6","GATA4","GCG","GCNT2","GDF15","GDF6","GFRAL","GH1","GHR","GHRL","GLIPR2","GNAI2","GPBAR1","GPER1","GPNMB","GPR183","GPR37","GPR37L1","GPR55","GRM1","GRM4","GRM5","GSDME","HAND2","HAVCR2","HCRTR1","HGF","HIPK2","HLA-DRB1","HMGB1","HRAS","HTR2A","HTR2B","HTR2C","IAPP","ICAM1","IGF1","IGF1R","IGF2","IGFBP3","IGFBP4","IGFBP6","IL11","IL1A","IL1B","IL26","IL34","IL6","INAVA","INHBA","INS","INSR","IQGAP1","IQGAP3","IRAK1","ITGA1","ITGB3","JAK2","JCAD","JUN","KDR","KISS1","KIT","KITLG","KL","KLB","KLHDC10","KSR1","LAMTOR1","LAMTOR2","LAMTOR3","LAPTM5","LEP","LGALS9","LIF","LILRA5","LPAR1","LPAR2","LPAR3","LRRK2","LTBR","MADD","MAGED1","MAP2K1","MAP2K2","MAP2K3","MAP2K4","MAP2K5","MAP2K6","MAP2K7","MAP3K10","MAP3K11","MAP3K12","MAP3K13","MAP3K3","MAP3K4","MAP3K5","MAP3K7","MAP4K1","MAP4K2","MAPK3","MAPK8IP1","MAPK8IP2","MAPK8IP3","MAPKBP1","MARCO","MBIP","MEF2C","MFAP3","MFHAS1","MID1","MIF","MINK1","MIR126","MIR181A2","MIR181B1","MIR181D","MIR21","MIR221","MIR222","MIR23A","MIR24-1","MIR27A","MIR27B","MIR519D","MIR92A1","MIRLET7B","MMP8","MOS","MST1R","MT3","MTURN","MUSK","MYD88","MYDGF","NAIP","NCF1","NDRG4","NDST1","NECAB2","NEK10","NELFE","NENF","NOD1","NOD2","NODAL","NOTCH1","NOTCH2","NOX1","NOX4","NPNT","NPSR1","NPTN","NPY","NPY5R","NRG1","NRP1","NTF3","NTRK1","NTRK2","NTRK3","OPRK1","OPRM1","OR2AT4","OSM","P2RX7","P2RY1","P2RY6","PAK1","PDCD10","PDE5A","PDE6G","PDE6H","PDE8A","PDGFA","PDGFB","PDGFC","PDGFD","PDGFRA","PDGFRB","PELI2","PHB1","PHB2","PIK3CG","PIK3R5","PIK3R6","PJA2","PLA2G1B","PLA2G2A","PLA2G5","PLCB1","PLCE1","PLCG2","PPIA","PRDX2","PRKCA","PRKCE","PRKCZ","PRKD2","PRMT1","PROK1","PRXL2C","PSEN1","PTK2B","PTPN1","PTPN11","PTPN22","PTPRC","PTPRJ","PYCARD","RAF1","RAMP3","RAP1A","RAP1B","RAPGEF2","RASGRP1","RASSF2","RB1CC1","RELL1","RELL2","RET","RIPK1","RIPK2","RIT2","ROBO1","ROCK1","ROCK2","ROR1","ROR2","RPS3","RYK","S100A12","S100A7","SASH1","SCIMP","SDCBP","SEMA3A","SEMA4C","SEMA7A","SERPINF2","SH3RF1","SH3RF2","SH3RF3","SHC1","SLAMF1","SLC30A10","SOD1","SORBS3","SOX2","SPAG9","SPHK1","SPI1","SPRY2","SRC","SSTR4","STK25","STK3","STK39","SYK","SYT14P1","TAB1","TAOK1","TAOK2","TAOK3","TBX1","TDGF1","TEK","TENM1","TGFA","TGFB1","TGFB2","TGFB3","TGFBR1","THBS1","THPO","TIRAP","TLR3","TLR4","TLR6","TLR9","TMEM106A","TNF","TNFAIP8L3","TNFRSF11A","TNFRSF19","TNFSF11","TNIK","TP73","TPBG","TPD52L1","TRAF1","TRAF2","TRAF3","TRAF4","TRAF5","TRAF6","TRAF7","TREM2","TRIM5","TRPV4","UNC5CL","VEGFA","WNT16","WNT5A","WNT7A","WNT7B","WWC1","XCL1","XCL2","XDH","XIAP","ZC3H12A","ZNF622"))

# Create factor level order for VlnPlot
mimicry_order1 <- c("Neoplastic Cells","Immune-like", "Mesenchymal-like", "Endothelial-like", "Developmental-like")
all_Swarbrick_breast_tumors_panCK_cnv_annotated[["main.ident"]] <- factor(x = all_Swarbrick_breast_tumors_panCK_cnv_annotated@meta.data$main.ident, levels = mimicry_order1)

mimicry_order2 <- c("Neoplastic Cells", "Lymphoid-like", "Myeloid-like", "Mesenchymal-like", "Endothelial-like", "Developmental-like")
all_Swarbrick_breast_tumors_panCK_cnv_annotated[["detailed.ident"]] <- factor(x = all_Swarbrick_breast_tumors_panCK_cnv_annotated@meta.data$detailed.ident, levels = mimicry_order2)


# Lim_Mammary_Stem 
# Score cells based on signature expression
all_Swarbrick_breast_tumors_panCK_cnv_annotated <- AddModuleScore(all_Swarbrick_breast_tumors_panCK_cnv_annotated, features = Lim_Mammary_Stem, name = "Lim_Mammary_Stem_Up")
# Draw feature plot
FeaturePlot(all_Swarbrick_breast_tumors_panCK_cnv_annotated, c("Lim_Mammary_Stem_Up1"))
# Draw Vln plot
VlnPlot(all_Swarbrick_breast_tumors_panCK_cnv_annotated, c("Lim_Mammary_Stem_Up1"), group.by = "main.ident", pt.size = 0)
ggsave("/R/R_Swarbrick/Swarbrick_Output/all_Swarbrick_breast_tumors_panCK_cnv_annotated.VLN.main.Lim_Mammary_Stem.tiff", plot = last_plot(), device = "tiff",
       scale = 1, width = 6, height = 10,
       dpi = 200, limitsize = TRUE)
VlnPlot(all_Swarbrick_breast_tumors_panCK_cnv_annotated, c("Lim_Mammary_Stem_Up1"), group.by = "detailed.ident", pt.size = 0)
ggsave("/R/R_Swarbrick/Swarbrick_Output/all_Swarbrick_breast_tumors_panCK_cnv_annotated.VLN.detailed.ident.Lim_Mammary_Stem.tiff", plot = last_plot(), device = "tiff",
       scale = 1, width = 8, height = 10,
       dpi = 200, limitsize = TRUE)



# Lim_Luminal_Progenitor 
# Score cells based on signature expression
all_Swarbrick_breast_tumors_panCK_cnv_annotated <- AddModuleScore(all_Swarbrick_breast_tumors_panCK_cnv_annotated, features = Lim_Luminal_Progenitor, name = "Lim_Luminal_Progenitor_Up")
# Draw feature plot
FeaturePlot(all_Swarbrick_breast_tumors_panCK_cnv_annotated, c("Lim_Luminal_Progenitor_Up1"))
# Draw Vln plot
VlnPlot(all_Swarbrick_breast_tumors_panCK_cnv_annotated, c("Lim_Luminal_Progenitor_Up1"), group.by = "main.ident", pt.size = 0)
ggsave("/R/R_Swarbrick/Swarbrick_Output/all_Swarbrick_breast_tumors_panCK_cnv_annotated.VLN.main.Lim_Luminal_Progenitor.tiff", plot = last_plot(), device = "tiff",
       scale = 1, width = 6, height = 10,
       dpi = 200, limitsize = TRUE)
VlnPlot(all_Swarbrick_breast_tumors_panCK_cnv_annotated, c("Lim_Luminal_Progenitor_Up1"), group.by = "detailed.ident", pt.size = 0)
ggsave("/R/R_Swarbrick/Swarbrick_Output/all_Swarbrick_breast_tumors_panCK_cnv_annotated.VLN.detailed.ident.Lim_Luminal_Progenitor.tiff", plot = last_plot(), device = "tiff",
       scale = 1, width = 8, height = 10,
       dpi = 200, limitsize = TRUE)



# Lim_Mature 
# Score cells based on signature expression
all_Swarbrick_breast_tumors_panCK_cnv_annotated <- AddModuleScore(all_Swarbrick_breast_tumors_panCK_cnv_annotated, features = Lim_Mature, name = "Lim_Mature_Up")
# Draw feature plot
FeaturePlot(all_Swarbrick_breast_tumors_panCK_cnv_annotated, c("Lim_Mature_Up1"))
# Draw Vln plot
VlnPlot(all_Swarbrick_breast_tumors_panCK_cnv_annotated, c("Lim_Mature_Up1"), group.by = "main.ident", pt.size = 0)
ggsave("/R/R_Swarbrick/Swarbrick_Output/all_Swarbrick_breast_tumors_panCK_cnv_annotated.VLN.main.Lim_Mature.tiff", plot = last_plot(), device = "tiff",
       scale = 1, width = 6, height = 10,
       dpi = 200, limitsize = TRUE)
VlnPlot(all_Swarbrick_breast_tumors_panCK_cnv_annotated, c("Lim_Mature_Up1"), group.by = "detailed.ident", pt.size = 0)
ggsave("/R/R_Swarbrick/Swarbrick_Output/all_Swarbrick_breast_tumors_panCK_cnv_annotated.VLN.detailed.ident.Lim_Mature.tiff", plot = last_plot(), device = "tiff",
       scale = 1, width = 8, height = 10,
       dpi = 200, limitsize = TRUE)



#################################################
###      Functional Signature Analysis        ###
#################################################

# BIOCARTA_NFKB_PATHWAY 
# Score cells based on signature expression
all_Swarbrick_breast_tumors_panCK_cnv_annotated <- AddModuleScore(all_Swarbrick_breast_tumors_panCK_cnv_annotated, features = BIOCARTA_NFKB_PATHWAY, name = "BIOCARTA_NFKB_PATHWAY")
# Draw feature plot
FeaturePlot(all_Swarbrick_breast_tumors_panCK_cnv_annotated, c("BIOCARTA_NFKB_PATHWAY1"))
# Draw Vln plot
VlnPlot(all_Swarbrick_breast_tumors_panCK_cnv_annotated, c("BIOCARTA_NFKB_PATHWAY1"), group.by = "main.ident", pt.size = 0)
ggsave("/R/R_Swarbrick/Swarbrick_Output/all_Swarbrick_breast_tumors_panCK_cnv_annotated.VLN.main.ident.BIOCARTA_NFKB_PATHWAY.tiff", plot = last_plot(), device = "tiff",
       scale = 1, width = 6, height = 10,
       dpi = 200, limitsize = TRUE)
VlnPlot(all_Swarbrick_breast_tumors_panCK_cnv_annotated, c("BIOCARTA_NFKB_PATHWAY1"), group.by = "detailed.ident", pt.size = 0)
ggsave("/R/R_Swarbrick/Swarbrick_Output/all_Swarbrick_breast_tumors_panCK_cnv_annotated.VLN.detailed.ident.BIOCARTA_NFKB_PATHWAY.tiff", plot = last_plot(), device = "tiff",
       scale = 1, width = 8, height = 10,
       dpi = 200, limitsize = TRUE)


# POSITIVE_REGULATION_OF_LEUKOCYTE_PROLIFERATION
# Score cells based on signature expression
all_Swarbrick_breast_tumors_panCK_cnv_annotated <- AddModuleScore(all_Swarbrick_breast_tumors_panCK_cnv_annotated, features = POSITIVE_REGULATION_OF_LEUKOCYTE_PROLIFERATION, name = "POSITIVE_REGULATION_OF_LEUKOCYTE_PROLIFERATION")
# Draw feature plot
FeaturePlot(all_Swarbrick_breast_tumors_panCK_cnv_annotated, c("POSITIVE_REGULATION_OF_LEUKOCYTE_PROLIFERATION1"))
# Draw Vln plot
VlnPlot(all_Swarbrick_breast_tumors_panCK_cnv_annotated, c("POSITIVE_REGULATION_OF_LEUKOCYTE_PROLIFERATION1"), group.by = "main.ident", pt.size = 0)
ggsave("/R/R_Swarbrick/Swarbrick_Output/all_Swarbrick_breast_tumors_panCK_cnv_annotated.VLN.main.POSITIVE_REGULATION_OF_LEUKOCYTE_PROLIFERATION.tiff", plot = last_plot(), device = "tiff",
       scale = 1, width = 6, height = 10,
       dpi = 200, limitsize = TRUE)
VlnPlot(all_Swarbrick_breast_tumors_panCK_cnv_annotated, c("POSITIVE_REGULATION_OF_LEUKOCYTE_PROLIFERATION1"), group.by = "detailed.ident", pt.size = 0)
ggsave("/R/R_Swarbrick/Swarbrick_Output/all_Swarbrick_breast_tumors_panCK_cnv_annotated.VLN.detailed.ident.POSITIVE_REGULATION_OF_LEUKOCYTE_PROLIFERATION.tiff", plot = last_plot(), device = "tiff",
       scale = 1, width = 8, height = 10,
       dpi = 200, limitsize = TRUE)


# POSITIVE_REGULATION_OF_MAPK_CASCADE
# Score cells based on signature expression
all_Swarbrick_breast_tumors_panCK_cnv_annotated <- AddModuleScore(all_Swarbrick_breast_tumors_panCK_cnv_annotated, features = POSITIVE_REGULATION_OF_MAPK_CASCADE, name = "POSITIVE_REGULATION_OF_MAPK_CASCADE")
# Draw feature plot
FeaturePlot(all_Swarbrick_breast_tumors_panCK_cnv_annotated, c("POSITIVE_REGULATION_OF_MAPK_CASCADE1"))
# Draw Vln plot
VlnPlot(all_Swarbrick_breast_tumors_panCK_cnv_annotated, c("POSITIVE_REGULATION_OF_MAPK_CASCADE1"), group.by = "main.ident", pt.size = 0)
ggsave("/R/R_Swarbrick/Swarbrick_Output/all_Swarbrick_breast_tumors_panCK_cnv_annotated.VLN.main.POSITIVE_REGULATION_OF_MAPK_CASCADE.tiff", plot = last_plot(), device = "tiff",
       scale = 1, width = 6, height = 10,
       dpi = 200, limitsize = TRUE)
VlnPlot(all_Swarbrick_breast_tumors_panCK_cnv_annotated, c("POSITIVE_REGULATION_OF_MAPK_CASCADE1"), group.by = "detailed.ident", pt.size = 0)
ggsave("/R/R_Swarbrick/Swarbrick_Output/all_Swarbrick_breast_tumors_panCK_cnv_annotated.VLN.detailed.ident.POSITIVE_REGULATION_OF_MAPK_CASCADE.tiff", plot = last_plot(), device = "tiff",
       scale = 1, width = 8, height = 10,
       dpi = 200, limitsize = TRUE)


# POSITIVE_REGULATION_OF_CYTOKINE_PRODUCTION
# Score cells based on signature expression
all_Swarbrick_breast_tumors_panCK_cnv_annotated <- AddModuleScore(all_Swarbrick_breast_tumors_panCK_cnv_annotated, features = POSITIVE_REGULATION_OF_CYTOKINE_PRODUCTION, name = "POSITIVE_REGULATION_OF_CYTOKINE_PRODUCTION")
# Draw feature plot
FeaturePlot(all_Swarbrick_breast_tumors_panCK_cnv_annotated, c("POSITIVE_REGULATION_OF_CYTOKINE_PRODUCTION1"))
# Draw Vln plot
VlnPlot(all_Swarbrick_breast_tumors_panCK_cnv_annotated, c("POSITIVE_REGULATION_OF_CYTOKINE_PRODUCTION1"), group.by = "main.ident", pt.size = 0)
ggsave("/R/R_Swarbrick/Swarbrick_Output/all_Swarbrick_breast_tumors_panCK_cnv_annotated.VLN.main.POSITIVE_REGULATION_OF_CYTOKINE_PRODUCTION.tiff", plot = last_plot(), device = "tiff",
       scale = 1, width = 6, height = 10,
       dpi = 200, limitsize = TRUE)
VlnPlot(all_Swarbrick_breast_tumors_panCK_cnv_annotated, c("POSITIVE_REGULATION_OF_CYTOKINE_PRODUCTION1"), group.by = "detailed.ident", pt.size = 0)
ggsave("/R/R_Swarbrick/Swarbrick_Output/all_Swarbrick_breast_tumors_panCK_cnv_annotated.VLN.detailed.ident.POSITIVE_REGULATION_OF_CYTOKINE_PRODUCTION.tiff", plot = last_plot(), device = "tiff",
       scale = 1, width = 8, height = 10,
       dpi = 200, limitsize = TRUE)


# POSITIVE_REGULATION_OF_LYMPHOCYTE_MIGRATION
# Score cells based on signature expression
all_Swarbrick_breast_tumors_panCK_cnv_annotated <- AddModuleScore(all_Swarbrick_breast_tumors_panCK_cnv_annotated, features = POSITIVE_REGULATION_OF_LYMPHOCYTE_MIGRATION, name = "POSITIVE_REGULATION_OF_LYMPHOCYTE_MIGRATION")
# Draw feature plot
FeaturePlot(all_Swarbrick_breast_tumors_panCK_cnv_annotated, c("POSITIVE_REGULATION_OF_LYMPHOCYTE_MIGRATION1"))
# Draw Vln plot
VlnPlot(all_Swarbrick_breast_tumors_panCK_cnv_annotated, c("POSITIVE_REGULATION_OF_LYMPHOCYTE_MIGRATION1"), group.by = "main.ident", pt.size = 0)
ggsave("/R/R_Swarbrick/Swarbrick_Output/all_Swarbrick_breast_tumors_panCK_cnv_annotated.VLN.main.POSITIVE_REGULATION_OF_LYMPHOCYTE_MIGRATION.tiff", plot = last_plot(), device = "tiff",
       scale = 1, width = 6, height = 10,
       dpi = 200, limitsize = TRUE)
VlnPlot(all_Swarbrick_breast_tumors_panCK_cnv_annotated, c("POSITIVE_REGULATION_OF_LYMPHOCYTE_MIGRATION1"), group.by = "detailed.ident", pt.size = 0)
ggsave("/R/R_Swarbrick/Swarbrick_Output/all_Swarbrick_breast_tumors_panCK_cnv_annotated.VLN.detailed.identPOSITIVE_REGULATION_OF_LYMPHOCYTE_MIGRATION.tiff", plot = last_plot(), device = "tiff",
       scale = 1, width = 8, height = 10,
       dpi = 200, limitsize = TRUE)


# Extract  Meta Data
all_Swarbrick_breast_tumors_panCK_cnv_annotated.signatures.meta.data <- as.data.frame(as.matrix(all_Swarbrick_breast_tumors_panCK_cnv_annotated@meta.data))

# Save Meta Data
write.csv(all_Swarbrick_breast_tumors_panCK_cnv_annotated.signatures.meta.data, file = "/R/R_Swarbrick/Swarbrick_Output/all_Swarbrick_breast_tumors_panCK_cnv_annotated.signatures.meta.data.csv")

# Remove Meta Data
rm(all_Swarbrick_breast_tumors_panCK_cnv_annotated.signatures.meta.data)

# Save Annotated Seurat Object with Signature Scoring in Meta Data
saveRDS(all_Swarbrick_breast_tumors_panCK_cnv_annotated, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Annotated/all_Swarbrick_breast_tumors_panCK_cnv_annotated.rds")

# Clear memory
rm(all_Swarbrick_breast_tumors_panCK_cnv_annotated)
rm(Lim_Mammary_Stem)
rm(Lim_Luminal_Progenitor)
rm(Lim_Mature)
rm(BIOCARTA_NFKB_PATHWAY)
rm(POSITIVE_REGULATION_OF_LEUKOCYTE_PROLIFERATION)
rm(POSITIVE_REGULATION_OF_MAPK_CASCADE)
rm(POSITIVE_REGULATION_OF_CYTOKINE_PRODUCTION)
rm(POSITIVE_REGULATION_OF_LYMPHOCYTE_MIGRATION)
rm(mimicry_order1)
rm(mimicry_order2)
gc()



#####################################################################
# Step 16: Evaluate Cell Cycle Status of CD69-pos vs CD69-neg cells # 
#####################################################################

###########################################################
#### Subject CD69-pos and CD69-neg for Cell Cycle Scoring # 
###########################################################
# Load SeuratObject
all_Swarbrick_breast_tumors_panCK_cnv_annotated <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Annotated/all_Swarbrick_breast_tumors_panCK_cnv_annotated.rds")

###################
# Subset CD69-pos #
###################
# Cells greater than cutoffs for at least one of the three main mammary epithelial cell types
CD69.pos.Swarbrick_panCK_cnv <- subset(x = all_Swarbrick_breast_tumors_panCK_cnv_annotated, subset = CD69 > 1)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
CD69.pos.Swarbrick_panCK_cnv <- FindVariableFeatures(CD69.pos.Swarbrick_panCK_cnv, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(CD69.pos.Swarbrick_panCK_cnv)
CD69.pos.Swarbrick_panCK_cnv <- ScaleData(CD69.pos.Swarbrick_panCK_cnv, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
CD69.pos.Swarbrick_panCK_cnv <- RunPCA(CD69.pos.Swarbrick_panCK_cnv, features = VariableFeatures(object = CD69.pos.Swarbrick_panCK_cnv))

# Examine and visualize PCA results a few different ways
print(CD69.pos.Swarbrick_panCK_cnv[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(CD69.pos.Swarbrick_panCK_cnv, dims = 1:2, reduction = "pca")
DimPlot(CD69.pos.Swarbrick_panCK_cnv, reduction = "pca")

# Visualize PCs
ElbowPlot(CD69.pos.Swarbrick_panCK_cnv)

CD69.pos.Swarbrick_panCK_cnv <- FindNeighbors(CD69.pos.Swarbrick_panCK_cnv, dims = 1:20)
CD69.pos.Swarbrick_panCK_cnv <- FindClusters(CD69.pos.Swarbrick_panCK_cnv, resolution = 0.2)

# Umap clustering
CD69.pos.Swarbrick_panCK_cnv <- RunUMAP(CD69.pos.Swarbrick_panCK_cnv, dims = 1:20)
DimPlot(CD69.pos.Swarbrick_panCK_cnv, reduction = "umap")
DimPlot(CD69.pos.Swarbrick_panCK_cnv, reduction = "umap", group.by = "orig.ident")

########
# SAVE #
########
saveRDS(CD69.pos.Swarbrick_panCK_cnv, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Annotated/CD69.pos.Swarbrick_panCK_cnv.rds")

#################################
# Score Cell Cycle for CD69-pos #
#################################

# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

# Score Cell Cycle
CD69.pos.Swarbrick_panCK_cnv <- CellCycleScoring(CD69.pos.Swarbrick_panCK_cnv, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# View cell cycle scores and phase assignments
head(CD69.pos.Swarbrick_panCK_cnv[[]])

# Visualize the distribution of cell cycle markers across
RidgePlot(CD69.pos.Swarbrick_panCK_cnv, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)

# Draw plots
DimPlot(CD69.pos.Swarbrick_panCK_cnv, raster = FALSE)
DimPlot(CD69.pos.Swarbrick_panCK_cnv, group.by = "main.ident")

# Extract  Meta Data
CD69.pos.Swarbrick_cell_cycle <- as.data.frame(as.matrix(CD69.pos.Swarbrick_panCK_cnv@meta.data))

# Save Meta Data
write.csv(CD69.pos.Swarbrick_cell_cycle, file = "/R/R_Swarbrick/Swarbrick_Output/CD69.pos.Swarbrick_cell_cycle.csv")

# Remove Meta Data
rm(CD69.pos.Swarbrick_cell_cycle)
# Remove Object
rm(CD69.pos.Swarbrick_panCK_cnv)


###################
# Subset CD69-neg #
###################
# Cells greater than cutoffs for at least one of the three main mammary epithelial cell types
CD69.neg.Swarbrick_panCK_cnv <- subset(x = all_Swarbrick_breast_tumors_panCK_cnv_annotated, subset = CD69 == 0)

####################################
# FIND VARIABLE FEATURES & SCALING #
####################################
# Note: Normalization already performed for samples

# Find Variable Features
CD69.neg.Swarbrick_panCK_cnv <- FindVariableFeatures(CD69.neg.Swarbrick_panCK_cnv, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(CD69.neg.Swarbrick_panCK_cnv)
CD69.neg.Swarbrick_panCK_cnv <- ScaleData(CD69.neg.Swarbrick_panCK_cnv, features = all.genes)

#####################
# REDUCE DIMENSIONS #
#####################
CD69.neg.Swarbrick_panCK_cnv <- RunPCA(CD69.neg.Swarbrick_panCK_cnv, features = VariableFeatures(object = CD69.neg.Swarbrick_panCK_cnv))

# Examine and visualize PCA results a few different ways
print(CD69.neg.Swarbrick_panCK_cnv[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(CD69.neg.Swarbrick_panCK_cnv, dims = 1:2, reduction = "pca")
DimPlot(CD69.neg.Swarbrick_panCK_cnv, reduction = "pca")

# Visualize PCs
ElbowPlot(CD69.neg.Swarbrick_panCK_cnv)

CD69.neg.Swarbrick_panCK_cnv <- FindNeighbors(CD69.neg.Swarbrick_panCK_cnv, dims = 1:40)
CD69.neg.Swarbrick_panCK_cnv <- FindClusters(CD69.neg.Swarbrick_panCK_cnv, resolution = 0.08)

# Umap clustering
CD69.neg.Swarbrick_panCK_cnv <- RunUMAP(CD69.neg.Swarbrick_panCK_cnv, dims = 1:40)
DimPlot(CD69.neg.Swarbrick_panCK_cnv, reduction = "umap")
DimPlot(CD69.neg.Swarbrick_panCK_cnv, reduction = "umap", group.by = "orig.ident")

########
# SAVE #
########
saveRDS(CD69.neg.Swarbrick_panCK_cnv, file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Annotated/CD69.neg.Swarbrick_panCK_cnv.rds")

#################################
# Score Cell Cycle for CD69-neg #
#################################

# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

# Score Cell Cycle
CD69.neg.Swarbrick_panCK_cnv <- CellCycleScoring(CD69.neg.Swarbrick_panCK_cnv, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# View cell cycle scores and phase assignments
head(CD69.neg.Swarbrick_panCK_cnv[[]])

# Visualize the distribution of cell cycle markers across
RidgePlot(CD69.neg.Swarbrick_panCK_cnv, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)

# Draw plots
DimPlot(CD69.neg.Swarbrick_panCK_cnv, raster = FALSE)
DimPlot(CD69.neg.Swarbrick_panCK_cnv, group.by = "main.ident")

# Extract  Meta Data
CD69.neg.Swarbrick_cell_cycle <- as.data.frame(as.matrix(CD69.neg.Swarbrick_panCK_cnv@meta.data))

# Save Meta Data
write.csv(CD69.neg.Swarbrick_cell_cycle, file = "/R/R_Swarbrick/Swarbrick_Output/CD69.neg.Swarbrick_cell_cycle.csv")

# Remove Meta Data
rm(CD69.neg.Swarbrick_cell_cycle)
# Remove Object
rm(CD69.neg.Swarbrick_panCK_cnv)
# Remove others
rm(s.genes)
rm(g2m.genes)
rm(all_Swarbrick_breast_tumors_panCK_cnv_annotated)
gc()


#########################################################
######### Quantify Cell Cycle Status by Patient #########
#########################################################

#############################################
# Neoplastic CD69-pos Cell Cycle by Patient #
#############################################

# Read Cell Cycle Meta Data
CD69.pos.Swarbrick_cell_cycle <- read.csv(file = "/R/R_Swarbrick/Swarbrick_Output/CD69.pos.Swarbrick_cell_cycle.csv")

# Check if "phase" column exists in your dataframe
if ("Phase" %in% colnames(CD69.pos.Swarbrick_cell_cycle)) {
  
  # Count the number of individuals in each Phase per group
  phase_counts_per_group <- table(CD69.pos.Swarbrick_cell_cycle$orig.ident, CD69.pos.Swarbrick_cell_cycle$Phase)
  
  # Convert table to data frame and add column names
  phase_counts_df <- as.data.frame.matrix(phase_counts_per_group)
  colnames(phase_counts_df) <- c("G1", "G2M", "S")
  
  # Print the counts
  print(phase_counts_df)
  
  # Save the counts as a CSV file (optional)
  write.csv(phase_counts_df, file = "/R/R_Swarbrick/Swarbrick_Output/CD69.pos.Swarbrick_cell_cycle_status_by_patient.csv", row.names = TRUE)
  
} else {
  print("The 'phase' column is not present in the dataframe.")
}


#############################################
# Neoplastic CD69-neg Cell Cycle by Patient #
#############################################

# Read Cell Cycle Meta Data
CD69.neg.Swarbrick_cell_cycle <- read.csv(file = "/R/R_Swarbrick/Swarbrick_Output/CD69.neg.Swarbrick_cell_cycle.csv")

# Check if "phase" column exists in your dataframe
if ("Phase" %in% colnames(CD69.neg.Swarbrick_cell_cycle)) {
  
  # Count the number of individuals in each Phase per group
  phase_counts_per_group <- table(CD69.neg.Swarbrick_cell_cycle$orig.ident, CD69.neg.Swarbrick_cell_cycle$Phase)
  
  # Convert table to data frame and add column names
  phase_counts_df <- as.data.frame.matrix(phase_counts_per_group)
  colnames(phase_counts_df) <- c("G1", "G2M", "S")
  
  # Print the counts
  print(phase_counts_df)
  
  # Save the counts as a CSV file (optional)
  write.csv(phase_counts_df, file = "/R/R_Swarbrick/Swarbrick_Output/CD69.neg.Swarbrick_cell_cycle_status_by_patient.csv", row.names = TRUE)
  
} else {
  print("The 'phase' column is not present in the dataframe.")
}



##############################################################
# Step 17: Assessing Other Features of Immune-Mimicked Cells #
##############################################################

########################################
# Assessing Features in Mimicked Cells #
########################################
# all_Swarbrick_breast_tumors_panCK_cnv_annotated <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Annotated/all_Swarbrick_breast_tumors_panCK_cnv_annotated.rds")

# Feature Assessment
VlnPlot(all_Swarbrick_breast_tumors_panCK_cnv_annotated, c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "detailed.ident", pt.size = 0.1)
ggsave("/R/R_Swarbrick/Swarbrick_Output/all_Swarbrick_breast_tumors_panCK_cnv_annotated.General_Features.tiff", plot = last_plot(), device = "tiff",
       scale = 1, width = 12, height = 10,
       dpi = 200, limitsize = TRUE)

# Feature Assessment - No Points 
VlnPlot(all_Swarbrick_breast_tumors_panCK_cnv_annotated, c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "detailed.ident", pt.size = 0)
ggsave("/R/R_Swarbrick/Swarbrick_Output/all_Swarbrick_breast_tumors_panCK_cnv_annotated.General_Features_No_Points.tiff", plot = last_plot(), device = "tiff",
       scale = 1, width = 12, height = 10,
       dpi = 200, limitsize = TRUE)

#########################################################################
# Extract pANN Doublet Score from Immune Mimicry Annotated SeuratObject #
#########################################################################
# Load Library
library(dplyr)

# Read RDS
# all_Swarbrick_breast_tumors_panCK_cnv_annotated <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Annotated/all_Swarbrick_breast_tumors_panCK_cnv_annotated.rds")

# Extract MetaData
all_Swarbrick_breast_tumors_panCK_cnv_annotated.meta.data <- as.data.frame(as.matrix(all_Swarbrick_breast_tumors_panCK_cnv_annotated@meta.data))

# Identify columns that start with "pANN"
pANN_columns <- grep("^pANN", colnames(all_Swarbrick_breast_tumors_panCK_cnv_annotated.meta.data), value = TRUE)

# all_Swarbrick_breast_tumors_panCK_cnv_annotated.meta.data the values from the pANN columns into a single column
all_Swarbrick_breast_tumors_panCK_cnv_annotated.meta.data <- all_Swarbrick_breast_tumors_panCK_cnv_annotated.meta.data %>%
  rowwise() %>%
  mutate(DoubletFinder_pANN = coalesce(!!!syms(pANN_columns)))

# Switch back to data frame
all_Swarbrick_breast_tumors_panCK_cnv_annotated.meta.data <- as.data.frame(all_Swarbrick_breast_tumors_panCK_cnv_annotated.meta.data)

# Convert the 'DoubletFinder_pANN' column to numeric
all_Swarbrick_breast_tumors_panCK_cnv_annotated.meta.data$DoubletFinder_pANN <- as.numeric(all_Swarbrick_breast_tumors_panCK_cnv_annotated.meta.data$DoubletFinder_pANN)

# Ensure the row names of the meta data match the cell identities in the Seurat object
rownames(all_Swarbrick_breast_tumors_panCK_cnv_annotated.meta.data) <- Cells(all_Swarbrick_breast_tumors_panCK_cnv_annotated)

# Add the new metadata to the Seurat object
all_Swarbrick_breast_tumors_panCK_cnv_annotated <- AddMetaData(object = all_Swarbrick_breast_tumors_panCK_cnv_annotated, metadata = all_Swarbrick_breast_tumors_panCK_cnv_annotated.meta.data)


# Create factor level order for VlnPlot
mimicry_order2 <- c(mimicry_order2 <- c("Neoplastic Cells", "Lymphoid-like", "Myeloid-like", "Mesenchymal-like", "Endothelial-like", "Developmental-like"))
all_Swarbrick_breast_tumors_panCK_cnv_annotated[["detailed.ident"]] <- factor(x = all_Swarbrick_breast_tumors_panCK_cnv_annotated@meta.data$detailed.ident, levels = mimicry_order2)

# pANN Assessment - Transparent Points 
VlnPlot(all_Swarbrick_breast_tumors_panCK_cnv_annotated, c("DoubletFinder_pANN"), group.by = "detailed.ident", pt.size = 0) + geom_jitter(width = 0.35, alpha = 0.15, size = 0.1) 
ggsave("/R/R_Swarbrick/Swarbrick_Output/all_Swarbrick_breast_tumors_panCK_cnv_annotated_pANN_VLN_Transparent_Points.tiff", plot = last_plot(), device = "tiff",
       scale = 1, width = 12, height = 10,
       dpi = 200, limitsize = TRUE)

# Remove objects
rm(all_Swarbrick_breast_tumors_panCK_cnv_annotated)
rm(all_Swarbrick_breast_tumors_panCK_cnv_annotated.meta.data)
rm(pANN_columns)
rm(mimicry_order2)



#########################################################
# Step 18: Generate Expression Matrices for Each Sample # 
#########################################################
# Load Libraries
library(Seurat)
library(dplyr)
library(patchwork)

################################################################
# Extract Expression Matrices per Sample: Immune Mimicry Genes #
################################################################

#########################################
# all_Swarbrick_breast_tumors_panCK_cnv #
#########################################
# Load data file
all_Swarbrick_breast_tumors_panCK_cnv <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Merged/all_Swarbrick_breast_tumors_panCK_cnv.rds")
# Extract Expression Matrix
all_Swarbrick_breast_tumors_panCK_cnv <- all_Swarbrick_breast_tumors_panCK_cnv[["RNA"]]$counts
all_Swarbrick_breast_tumors_panCK_cnv <- as.matrix(all_Swarbrick_breast_tumors_panCK_cnv, 'sparseMatrix')
all_Swarbrick_breast_tumors_panCK_cnv <- subset(all_Swarbrick_breast_tumors_panCK_cnv, rownames(all_Swarbrick_breast_tumors_panCK_cnv) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
all_Swarbrick_breast_tumors_panCK_cnv <-t(all_Swarbrick_breast_tumors_panCK_cnv)
write.csv(all_Swarbrick_breast_tumors_panCK_cnv, file = "/R/R_Swarbrick/Swarbrick_Expression_Matrices/Swarbrick_Expression_Matrices_Output/all_Swarbrick_breast_tumors_panCK_cnv_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(all_Swarbrick_breast_tumors_panCK_cnv)
gc()


##################################
# Swarbrick_CID3941_ER_panCK_cnv #
##################################
# Load Data
Swarbrick_CID3941_ER_panCK_cnv <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_inferCNV/Swarbrick_CID3941_ER_panCK_cnv.rds")
# Extract Expression Matrix
Swarbrick_CID3941_ER_panCK_cnv <- Swarbrick_CID3941_ER_panCK_cnv[["RNA"]]$counts
Swarbrick_CID3941_ER_panCK_cnv <- as.matrix(Swarbrick_CID3941_ER_panCK_cnv, 'sparseMatrix')
Swarbrick_CID3941_ER_panCK_cnv <- subset(Swarbrick_CID3941_ER_panCK_cnv, rownames(Swarbrick_CID3941_ER_panCK_cnv) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Swarbrick_CID3941_ER_panCK_cnv <-t(Swarbrick_CID3941_ER_panCK_cnv)
write.csv(Swarbrick_CID3941_ER_panCK_cnv, file = "/R/R_Swarbrick/Swarbrick_Expression_Matrices/Swarbrick_Expression_Matrices_Output/Swarbrick_CID3941_ER_panCK_cnv_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Swarbrick_CID3941_ER_panCK_cnv)
gc()


##################################
# Swarbrick_CID3948_ER_panCK_cnv #
##################################
# Load Data
Swarbrick_CID3948_ER_panCK_cnv <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_inferCNV/Swarbrick_CID3948_ER_panCK_cnv.rds")
# Extract Expression Matrix
Swarbrick_CID3948_ER_panCK_cnv <- Swarbrick_CID3948_ER_panCK_cnv[["RNA"]]$counts
Swarbrick_CID3948_ER_panCK_cnv <- as.matrix(Swarbrick_CID3948_ER_panCK_cnv, 'sparseMatrix')
Swarbrick_CID3948_ER_panCK_cnv <- subset(Swarbrick_CID3948_ER_panCK_cnv, rownames(Swarbrick_CID3948_ER_panCK_cnv) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Swarbrick_CID3948_ER_panCK_cnv <-t(Swarbrick_CID3948_ER_panCK_cnv)
write.csv(Swarbrick_CID3948_ER_panCK_cnv, file = "/R/R_Swarbrick/Swarbrick_Expression_Matrices/Swarbrick_Expression_Matrices_Output/Swarbrick_CID3948_ER_panCK_cnv_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Swarbrick_CID3948_ER_panCK_cnv)
gc()


##################################
# Swarbrick_CID4040_ER_panCK_cnv #
##################################
# Load Data
Swarbrick_CID4040_ER_panCK_cnv <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_inferCNV/Swarbrick_CID4040_ER_panCK_cnv.rds")
# Extract Expression Matrix
Swarbrick_CID4040_ER_panCK_cnv <- Swarbrick_CID4040_ER_panCK_cnv[["RNA"]]$counts
Swarbrick_CID4040_ER_panCK_cnv <- as.matrix(Swarbrick_CID4040_ER_panCK_cnv, 'sparseMatrix')
Swarbrick_CID4040_ER_panCK_cnv <- subset(Swarbrick_CID4040_ER_panCK_cnv, rownames(Swarbrick_CID4040_ER_panCK_cnv) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Swarbrick_CID4040_ER_panCK_cnv <-t(Swarbrick_CID4040_ER_panCK_cnv)
write.csv(Swarbrick_CID4040_ER_panCK_cnv, file = "/R/R_Swarbrick/Swarbrick_Expression_Matrices/Swarbrick_Expression_Matrices_Output/Swarbrick_CID4040_ER_panCK_cnv_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Swarbrick_CID4040_ER_panCK_cnv)
gc()


##################################
# Swarbrick_CID4067_ER_panCK_cnv #
##################################
# Load Data
Swarbrick_CID4067_ER_panCK_cnv <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_inferCNV/Swarbrick_CID4067_ER_panCK_cnv.rds")
# Extract Expression Matrix
Swarbrick_CID4067_ER_panCK_cnv <- Swarbrick_CID4067_ER_panCK_cnv[["RNA"]]$counts
Swarbrick_CID4067_ER_panCK_cnv <- as.matrix(Swarbrick_CID4067_ER_panCK_cnv, 'sparseMatrix')
Swarbrick_CID4067_ER_panCK_cnv <- subset(Swarbrick_CID4067_ER_panCK_cnv, rownames(Swarbrick_CID4067_ER_panCK_cnv) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Swarbrick_CID4067_ER_panCK_cnv <-t(Swarbrick_CID4067_ER_panCK_cnv)
write.csv(Swarbrick_CID4067_ER_panCK_cnv, file = "/R/R_Swarbrick/Swarbrick_Expression_Matrices/Swarbrick_Expression_Matrices_Output/Swarbrick_CID4067_ER_panCK_cnv_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Swarbrick_CID4067_ER_panCK_cnv)
gc()


###################################
# Swarbrick_CID4290A_ER_panCK_cnv #
###################################
# Load Data
Swarbrick_CID4290A_ER_panCK_cnv <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_inferCNV/Swarbrick_CID4290A_ER_panCK_cnv.rds")
# Extract Expression Matrix
Swarbrick_CID4290A_ER_panCK_cnv <- Swarbrick_CID4290A_ER_panCK_cnv[["RNA"]]$counts
Swarbrick_CID4290A_ER_panCK_cnv <- as.matrix(Swarbrick_CID4290A_ER_panCK_cnv, 'sparseMatrix')
Swarbrick_CID4290A_ER_panCK_cnv <- subset(Swarbrick_CID4290A_ER_panCK_cnv, rownames(Swarbrick_CID4290A_ER_panCK_cnv) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Swarbrick_CID4290A_ER_panCK_cnv <-t(Swarbrick_CID4290A_ER_panCK_cnv)
write.csv(Swarbrick_CID4290A_ER_panCK_cnv, file = "/R/R_Swarbrick/Swarbrick_Expression_Matrices/Swarbrick_Expression_Matrices_Output/Swarbrick_CID4290A_ER_panCK_cnv_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Swarbrick_CID4290A_ER_panCK_cnv)
gc()


##################################
# Swarbrick_CID4398_ER_panCK_cnv #
##################################
# Load Data
Swarbrick_CID4398_ER_panCK_cnv <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_inferCNV/Swarbrick_CID4398_ER_panCK_cnv.rds")
# Extract Expression Matrix
Swarbrick_CID4398_ER_panCK_cnv <- Swarbrick_CID4398_ER_panCK_cnv[["RNA"]]$counts
Swarbrick_CID4398_ER_panCK_cnv <- as.matrix(Swarbrick_CID4398_ER_panCK_cnv, 'sparseMatrix')
Swarbrick_CID4398_ER_panCK_cnv <- subset(Swarbrick_CID4398_ER_panCK_cnv, rownames(Swarbrick_CID4398_ER_panCK_cnv) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Swarbrick_CID4398_ER_panCK_cnv <-t(Swarbrick_CID4398_ER_panCK_cnv)
write.csv(Swarbrick_CID4398_ER_panCK_cnv, file = "/R/R_Swarbrick/Swarbrick_Expression_Matrices/Swarbrick_Expression_Matrices_Output/Swarbrick_CID4398_ER_panCK_cnv_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Swarbrick_CID4398_ER_panCK_cnv)
gc()


##################################
# Swarbrick_CID4461_ER_panCK_cnv #
##################################
# Load Data
Swarbrick_CID4461_ER_panCK_cnv <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_inferCNV/Swarbrick_CID4461_ER_panCK_cnv.rds")
# Extract Expression Matrix
Swarbrick_CID4461_ER_panCK_cnv <- Swarbrick_CID4461_ER_panCK_cnv[["RNA"]]$counts
Swarbrick_CID4461_ER_panCK_cnv <- as.matrix(Swarbrick_CID4461_ER_panCK_cnv, 'sparseMatrix')
Swarbrick_CID4461_ER_panCK_cnv <- subset(Swarbrick_CID4461_ER_panCK_cnv, rownames(Swarbrick_CID4461_ER_panCK_cnv) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Swarbrick_CID4461_ER_panCK_cnv <-t(Swarbrick_CID4461_ER_panCK_cnv)
write.csv(Swarbrick_CID4461_ER_panCK_cnv, file = "/R/R_Swarbrick/Swarbrick_Expression_Matrices/Swarbrick_Expression_Matrices_Output/Swarbrick_CID4461_ER_panCK_cnv_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Swarbrick_CID4461_ER_panCK_cnv)
gc()


##################################
# Swarbrick_CID4463_ER_panCK_cnv #
##################################
# Load Data
Swarbrick_CID4463_ER_panCK_cnv <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_inferCNV/Swarbrick_CID4463_ER_panCK_cnv.rds")
# Extract Expression Matrix
Swarbrick_CID4463_ER_panCK_cnv <- Swarbrick_CID4463_ER_panCK_cnv[["RNA"]]$counts
Swarbrick_CID4463_ER_panCK_cnv <- as.matrix(Swarbrick_CID4463_ER_panCK_cnv, 'sparseMatrix')
Swarbrick_CID4463_ER_panCK_cnv <- subset(Swarbrick_CID4463_ER_panCK_cnv, rownames(Swarbrick_CID4463_ER_panCK_cnv) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Swarbrick_CID4463_ER_panCK_cnv <-t(Swarbrick_CID4463_ER_panCK_cnv)
write.csv(Swarbrick_CID4463_ER_panCK_cnv, file = "/R/R_Swarbrick/Swarbrick_Expression_Matrices/Swarbrick_Expression_Matrices_Output/Swarbrick_CID4463_ER_panCK_cnv_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Swarbrick_CID4463_ER_panCK_cnv)
gc()


##################################
# Swarbrick_CID4471_ER_panCK_cnv #
##################################
# Load Data
Swarbrick_CID4471_ER_panCK_cnv <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_inferCNV/Swarbrick_CID4471_ER_panCK_cnv.rds")
# Extract Expression Matrix
Swarbrick_CID4471_ER_panCK_cnv <- Swarbrick_CID4471_ER_panCK_cnv[["RNA"]]$counts
Swarbrick_CID4471_ER_panCK_cnv <- as.matrix(Swarbrick_CID4471_ER_panCK_cnv, 'sparseMatrix')
Swarbrick_CID4471_ER_panCK_cnv <- subset(Swarbrick_CID4471_ER_panCK_cnv, rownames(Swarbrick_CID4471_ER_panCK_cnv) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Swarbrick_CID4471_ER_panCK_cnv <-t(Swarbrick_CID4471_ER_panCK_cnv)
write.csv(Swarbrick_CID4471_ER_panCK_cnv, file = "/R/R_Swarbrick/Swarbrick_Expression_Matrices/Swarbrick_Expression_Matrices_Output/Swarbrick_CID4471_ER_panCK_cnv_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Swarbrick_CID4471_ER_panCK_cnv)
gc()


###################################
# Swarbrick_CID4530N_ER_panCK_cnv #
###################################
# Load Data
Swarbrick_CID4530N_ER_panCK_cnv <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_inferCNV/Swarbrick_CID4530N_ER_panCK_cnv.rds")
# Extract Expression Matrix
Swarbrick_CID4530N_ER_panCK_cnv <- Swarbrick_CID4530N_ER_panCK_cnv[["RNA"]]$counts
Swarbrick_CID4530N_ER_panCK_cnv <- as.matrix(Swarbrick_CID4530N_ER_panCK_cnv, 'sparseMatrix')
Swarbrick_CID4530N_ER_panCK_cnv <- subset(Swarbrick_CID4530N_ER_panCK_cnv, rownames(Swarbrick_CID4530N_ER_panCK_cnv) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Swarbrick_CID4530N_ER_panCK_cnv <-t(Swarbrick_CID4530N_ER_panCK_cnv)
write.csv(Swarbrick_CID4530N_ER_panCK_cnv, file = "/R/R_Swarbrick/Swarbrick_Expression_Matrices/Swarbrick_Expression_Matrices_Output/Swarbrick_CID4530N_ER_panCK_cnv_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Swarbrick_CID4530N_ER_panCK_cnv)
gc()


##################################
# Swarbrick_CID4535_ER_panCK_cnv #
##################################
# Load Data
Swarbrick_CID4535_ER_panCK_cnv <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_inferCNV/Swarbrick_CID4535_ER_panCK_cnv.rds")
# Extract Expression Matrix
Swarbrick_CID4535_ER_panCK_cnv <- Swarbrick_CID4535_ER_panCK_cnv[["RNA"]]$counts
Swarbrick_CID4535_ER_panCK_cnv <- as.matrix(Swarbrick_CID4535_ER_panCK_cnv, 'sparseMatrix')
Swarbrick_CID4535_ER_panCK_cnv <- subset(Swarbrick_CID4535_ER_panCK_cnv, rownames(Swarbrick_CID4535_ER_panCK_cnv) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Swarbrick_CID4535_ER_panCK_cnv <-t(Swarbrick_CID4535_ER_panCK_cnv)
write.csv(Swarbrick_CID4535_ER_panCK_cnv, file = "/R/R_Swarbrick/Swarbrick_Expression_Matrices/Swarbrick_Expression_Matrices_Output/Swarbrick_CID4535_ER_panCK_cnv_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Swarbrick_CID4535_ER_panCK_cnv)
gc()


####################################
# Swarbrick_CID3586_HER2_panCK_cnv #
####################################
# Load Data
Swarbrick_CID3586_HER2_panCK_cnv <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_inferCNV/Swarbrick_CID3586_HER2_panCK_cnv.rds")
# Extract Expression Matrix
Swarbrick_CID3586_HER2_panCK_cnv <- Swarbrick_CID3586_HER2_panCK_cnv[["RNA"]]$counts
Swarbrick_CID3586_HER2_panCK_cnv <- as.matrix(Swarbrick_CID3586_HER2_panCK_cnv, 'sparseMatrix')
Swarbrick_CID3586_HER2_panCK_cnv <- subset(Swarbrick_CID3586_HER2_panCK_cnv, rownames(Swarbrick_CID3586_HER2_panCK_cnv) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Swarbrick_CID3586_HER2_panCK_cnv <-t(Swarbrick_CID3586_HER2_panCK_cnv)
write.csv(Swarbrick_CID3586_HER2_panCK_cnv, file = "/R/R_Swarbrick/Swarbrick_Expression_Matrices/Swarbrick_Expression_Matrices_Output/Swarbrick_CID3586_HER2_panCK_cnv_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Swarbrick_CID3586_HER2_panCK_cnv)
gc()


####################################
# Swarbrick_CID3838_HER2_panCK_cnv #
####################################
# Load Data
Swarbrick_CID3838_HER2_panCK_cnv <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_inferCNV/Swarbrick_CID3838_HER2_panCK_cnv.rds")
# Extract Expression Matrix
Swarbrick_CID3838_HER2_panCK_cnv <- Swarbrick_CID3838_HER2_panCK_cnv[["RNA"]]$counts
Swarbrick_CID3838_HER2_panCK_cnv <- as.matrix(Swarbrick_CID3838_HER2_panCK_cnv, 'sparseMatrix')
Swarbrick_CID3838_HER2_panCK_cnv <- subset(Swarbrick_CID3838_HER2_panCK_cnv, rownames(Swarbrick_CID3838_HER2_panCK_cnv) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Swarbrick_CID3838_HER2_panCK_cnv <-t(Swarbrick_CID3838_HER2_panCK_cnv)
write.csv(Swarbrick_CID3838_HER2_panCK_cnv, file = "/R/R_Swarbrick/Swarbrick_Expression_Matrices/Swarbrick_Expression_Matrices_Output/Swarbrick_CID3838_HER2_panCK_cnv_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Swarbrick_CID3838_HER2_panCK_cnv)
gc()


####################################
# Swarbrick_CID3921_HER2_panCK_cnv #
####################################
# Load Data
Swarbrick_CID3921_HER2_panCK_cnv <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_inferCNV/Swarbrick_CID3921_HER2_panCK_cnv.rds")
# Extract Expression Matrix
Swarbrick_CID3921_HER2_panCK_cnv <- Swarbrick_CID3921_HER2_panCK_cnv[["RNA"]]$counts
Swarbrick_CID3921_HER2_panCK_cnv <- as.matrix(Swarbrick_CID3921_HER2_panCK_cnv, 'sparseMatrix')
Swarbrick_CID3921_HER2_panCK_cnv <- subset(Swarbrick_CID3921_HER2_panCK_cnv, rownames(Swarbrick_CID3921_HER2_panCK_cnv) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Swarbrick_CID3921_HER2_panCK_cnv <-t(Swarbrick_CID3921_HER2_panCK_cnv)
write.csv(Swarbrick_CID3921_HER2_panCK_cnv, file = "/R/R_Swarbrick/Swarbrick_Expression_Matrices/Swarbrick_Expression_Matrices_Output/Swarbrick_CID3921_HER2_panCK_cnv_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Swarbrick_CID3921_HER2_panCK_cnv)
gc()


####################################
# Swarbrick_CID4066_HER2_panCK_cnv #
####################################
# Load Data
Swarbrick_CID4066_HER2_panCK_cnv <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_inferCNV/Swarbrick_CID4066_HER2_panCK_cnv.rds")
# Extract Expression Matrix
Swarbrick_CID4066_HER2_panCK_cnv <- Swarbrick_CID4066_HER2_panCK_cnv[["RNA"]]$counts
Swarbrick_CID4066_HER2_panCK_cnv <- as.matrix(Swarbrick_CID4066_HER2_panCK_cnv, 'sparseMatrix')
Swarbrick_CID4066_HER2_panCK_cnv <- subset(Swarbrick_CID4066_HER2_panCK_cnv, rownames(Swarbrick_CID4066_HER2_panCK_cnv) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Swarbrick_CID4066_HER2_panCK_cnv <-t(Swarbrick_CID4066_HER2_panCK_cnv)
write.csv(Swarbrick_CID4066_HER2_panCK_cnv, file = "/R/R_Swarbrick/Swarbrick_Expression_Matrices/Swarbrick_Expression_Matrices_Output/Swarbrick_CID4066_HER2_panCK_cnv_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Swarbrick_CID4066_HER2_panCK_cnv)
gc()


#####################################
# Swarbrick_CID45171_HER2_panCK_cnv #
#####################################
# Load Data
Swarbrick_CID45171_HER2_panCK_cnv <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_inferCNV/Swarbrick_CID45171_HER2_panCK_cnv.rds")
# Extract Expression Matrix
Swarbrick_CID45171_HER2_panCK_cnv <- Swarbrick_CID45171_HER2_panCK_cnv[["RNA"]]$counts
Swarbrick_CID45171_HER2_panCK_cnv <- as.matrix(Swarbrick_CID45171_HER2_panCK_cnv, 'sparseMatrix')
Swarbrick_CID45171_HER2_panCK_cnv <- subset(Swarbrick_CID45171_HER2_panCK_cnv, rownames(Swarbrick_CID45171_HER2_panCK_cnv) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Swarbrick_CID45171_HER2_panCK_cnv <-t(Swarbrick_CID45171_HER2_panCK_cnv)
write.csv(Swarbrick_CID45171_HER2_panCK_cnv, file = "/R/R_Swarbrick/Swarbrick_Expression_Matrices/Swarbrick_Expression_Matrices_Output/Swarbrick_CID45171_HER2_panCK_cnv_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Swarbrick_CID45171_HER2_panCK_cnv)
gc()


####################################
# Swarbrick_CID3946_TNBC_panCK_cnv #
####################################
# Load Data
Swarbrick_CID3946_TNBC_panCK_cnv <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_inferCNV/Swarbrick_CID3946_TNBC_panCK_cnv.rds")
# Extract Expression Matrix
Swarbrick_CID3946_TNBC_panCK_cnv <- Swarbrick_CID3946_TNBC_panCK_cnv[["RNA"]]$counts
Swarbrick_CID3946_TNBC_panCK_cnv <- as.matrix(Swarbrick_CID3946_TNBC_panCK_cnv, 'sparseMatrix')
Swarbrick_CID3946_TNBC_panCK_cnv <- subset(Swarbrick_CID3946_TNBC_panCK_cnv, rownames(Swarbrick_CID3946_TNBC_panCK_cnv) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Swarbrick_CID3946_TNBC_panCK_cnv <-t(Swarbrick_CID3946_TNBC_panCK_cnv)
write.csv(Swarbrick_CID3946_TNBC_panCK_cnv, file = "/R/R_Swarbrick/Swarbrick_Expression_Matrices/Swarbrick_Expression_Matrices_Output/Swarbrick_CID3946_TNBC_panCK_cnv_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Swarbrick_CID3946_TNBC_panCK_cnv)
gc()


####################################
# Swarbrick_CID3963_TNBC_panCK_cnv #
####################################
# Load Data
Swarbrick_CID3963_TNBC_panCK_cnv <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_inferCNV/Swarbrick_CID3963_TNBC_panCK_cnv.rds")
# Extract Expression Matrix
Swarbrick_CID3963_TNBC_panCK_cnv <- Swarbrick_CID3963_TNBC_panCK_cnv[["RNA"]]$counts
Swarbrick_CID3963_TNBC_panCK_cnv <- as.matrix(Swarbrick_CID3963_TNBC_panCK_cnv, 'sparseMatrix')
Swarbrick_CID3963_TNBC_panCK_cnv <- subset(Swarbrick_CID3963_TNBC_panCK_cnv, rownames(Swarbrick_CID3963_TNBC_panCK_cnv) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Swarbrick_CID3963_TNBC_panCK_cnv <-t(Swarbrick_CID3963_TNBC_panCK_cnv)
write.csv(Swarbrick_CID3963_TNBC_panCK_cnv, file = "/R/R_Swarbrick/Swarbrick_Expression_Matrices/Swarbrick_Expression_Matrices_Output/Swarbrick_CID3963_TNBC_panCK_cnv_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Swarbrick_CID3963_TNBC_panCK_cnv)
gc()


####################################
# Swarbrick_CID4465_TNBC_panCK_cnv #
####################################
# Load Data
Swarbrick_CID4465_TNBC_panCK_cnv <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_inferCNV/Swarbrick_CID4465_TNBC_panCK_cnv.rds")
# Extract Expression Matrix
Swarbrick_CID4465_TNBC_panCK_cnv <- Swarbrick_CID4465_TNBC_panCK_cnv[["RNA"]]$counts
Swarbrick_CID4465_TNBC_panCK_cnv <- as.matrix(Swarbrick_CID4465_TNBC_panCK_cnv, 'sparseMatrix')
Swarbrick_CID4465_TNBC_panCK_cnv <- subset(Swarbrick_CID4465_TNBC_panCK_cnv, rownames(Swarbrick_CID4465_TNBC_panCK_cnv) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Swarbrick_CID4465_TNBC_panCK_cnv <-t(Swarbrick_CID4465_TNBC_panCK_cnv)
write.csv(Swarbrick_CID4465_TNBC_panCK_cnv, file = "/R/R_Swarbrick/Swarbrick_Expression_Matrices/Swarbrick_Expression_Matrices_Output/Swarbrick_CID4465_TNBC_panCK_cnv_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Swarbrick_CID4465_TNBC_panCK_cnv)
gc()


####################################
# Swarbrick_CID4495_TNBC_panCK_cnv #
####################################
# Load Data
Swarbrick_CID4495_TNBC_panCK_cnv <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_inferCNV/Swarbrick_CID4495_TNBC_panCK_cnv.rds")
# Extract Expression Matrix
Swarbrick_CID4495_TNBC_panCK_cnv <- Swarbrick_CID4495_TNBC_panCK_cnv[["RNA"]]$counts
Swarbrick_CID4495_TNBC_panCK_cnv <- as.matrix(Swarbrick_CID4495_TNBC_panCK_cnv, 'sparseMatrix')
Swarbrick_CID4495_TNBC_panCK_cnv <- subset(Swarbrick_CID4495_TNBC_panCK_cnv, rownames(Swarbrick_CID4495_TNBC_panCK_cnv) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Swarbrick_CID4495_TNBC_panCK_cnv <-t(Swarbrick_CID4495_TNBC_panCK_cnv)
write.csv(Swarbrick_CID4495_TNBC_panCK_cnv, file = "/R/R_Swarbrick/Swarbrick_Expression_Matrices/Swarbrick_Expression_Matrices_Output/Swarbrick_CID4495_TNBC_panCK_cnv_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Swarbrick_CID4495_TNBC_panCK_cnv)
gc()


####################################
# Swarbrick_CID4513_TNBC_panCK_cnv #
####################################
# Load Data
Swarbrick_CID4513_TNBC_panCK_cnv <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_inferCNV/Swarbrick_CID4513_TNBC_panCK_cnv.rds")
# Extract Expression Matrix
Swarbrick_CID4513_TNBC_panCK_cnv <- Swarbrick_CID4513_TNBC_panCK_cnv[["RNA"]]$counts
Swarbrick_CID4513_TNBC_panCK_cnv <- as.matrix(Swarbrick_CID4513_TNBC_panCK_cnv, 'sparseMatrix')
Swarbrick_CID4513_TNBC_panCK_cnv <- subset(Swarbrick_CID4513_TNBC_panCK_cnv, rownames(Swarbrick_CID4513_TNBC_panCK_cnv) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Swarbrick_CID4513_TNBC_panCK_cnv <-t(Swarbrick_CID4513_TNBC_panCK_cnv)
write.csv(Swarbrick_CID4513_TNBC_panCK_cnv, file = "/R/R_Swarbrick/Swarbrick_Expression_Matrices/Swarbrick_Expression_Matrices_Output/Swarbrick_CID4513_TNBC_panCK_cnv_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Swarbrick_CID4513_TNBC_panCK_cnv)
gc()


####################################
# Swarbrick_CID4515_TNBC_panCK_cnv #
####################################
# Load Data
Swarbrick_CID4515_TNBC_panCK_cnv <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_inferCNV/Swarbrick_CID4515_TNBC_panCK_cnv.rds")
# Extract Expression Matrix
Swarbrick_CID4515_TNBC_panCK_cnv <- Swarbrick_CID4515_TNBC_panCK_cnv[["RNA"]]$counts
Swarbrick_CID4515_TNBC_panCK_cnv <- as.matrix(Swarbrick_CID4515_TNBC_panCK_cnv, 'sparseMatrix')
Swarbrick_CID4515_TNBC_panCK_cnv <- subset(Swarbrick_CID4515_TNBC_panCK_cnv, rownames(Swarbrick_CID4515_TNBC_panCK_cnv) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Swarbrick_CID4515_TNBC_panCK_cnv <-t(Swarbrick_CID4515_TNBC_panCK_cnv)
write.csv(Swarbrick_CID4515_TNBC_panCK_cnv, file = "/R/R_Swarbrick/Swarbrick_Expression_Matrices/Swarbrick_Expression_Matrices_Output/Swarbrick_CID4515_TNBC_panCK_cnv_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Swarbrick_CID4515_TNBC_panCK_cnv)
gc()


####################################
# Swarbrick_CID4523_TNBC_panCK_cnv #
####################################
# Load Data
Swarbrick_CID4523_TNBC_panCK_cnv <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_inferCNV/Swarbrick_CID4523_TNBC_panCK_cnv.rds")
# Extract Expression Matrix
Swarbrick_CID4523_TNBC_panCK_cnv <- Swarbrick_CID4523_TNBC_panCK_cnv[["RNA"]]$counts
Swarbrick_CID4523_TNBC_panCK_cnv <- as.matrix(Swarbrick_CID4523_TNBC_panCK_cnv, 'sparseMatrix')
Swarbrick_CID4523_TNBC_panCK_cnv <- subset(Swarbrick_CID4523_TNBC_panCK_cnv, rownames(Swarbrick_CID4523_TNBC_panCK_cnv) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Swarbrick_CID4523_TNBC_panCK_cnv <-t(Swarbrick_CID4523_TNBC_panCK_cnv)
write.csv(Swarbrick_CID4523_TNBC_panCK_cnv, file = "/R/R_Swarbrick/Swarbrick_Expression_Matrices/Swarbrick_Expression_Matrices_Output/Swarbrick_CID4523_TNBC_panCK_cnv_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Swarbrick_CID4523_TNBC_panCK_cnv)
gc()


#####################################
# Swarbrick_CID44041_TNBC_panCK_cnv #
#####################################
# Load Data
Swarbrick_CID44041_TNBC_panCK_cnv <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_inferCNV/Swarbrick_CID44041_TNBC_panCK_cnv.rds")
# Extract Expression Matrix
Swarbrick_CID44041_TNBC_panCK_cnv <- Swarbrick_CID44041_TNBC_panCK_cnv[["RNA"]]$counts
Swarbrick_CID44041_TNBC_panCK_cnv <- as.matrix(Swarbrick_CID44041_TNBC_panCK_cnv, 'sparseMatrix')
Swarbrick_CID44041_TNBC_panCK_cnv <- subset(Swarbrick_CID44041_TNBC_panCK_cnv, rownames(Swarbrick_CID44041_TNBC_panCK_cnv) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Swarbrick_CID44041_TNBC_panCK_cnv <-t(Swarbrick_CID44041_TNBC_panCK_cnv)
write.csv(Swarbrick_CID44041_TNBC_panCK_cnv, file = "/R/R_Swarbrick/Swarbrick_Expression_Matrices/Swarbrick_Expression_Matrices_Output/Swarbrick_CID44041_TNBC_panCK_cnv_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Swarbrick_CID44041_TNBC_panCK_cnv)
gc()


#####################################
# Swarbrick_CID44971_TNBC_panCK_cnv #
#####################################
# Load Data
Swarbrick_CID44971_TNBC_panCK_cnv <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_inferCNV/Swarbrick_CID44971_TNBC_panCK_cnv.rds")
# Extract Expression Matrix
Swarbrick_CID44971_TNBC_panCK_cnv <- Swarbrick_CID44971_TNBC_panCK_cnv[["RNA"]]$counts
Swarbrick_CID44971_TNBC_panCK_cnv <- as.matrix(Swarbrick_CID44971_TNBC_panCK_cnv, 'sparseMatrix')
Swarbrick_CID44971_TNBC_panCK_cnv <- subset(Swarbrick_CID44971_TNBC_panCK_cnv, rownames(Swarbrick_CID44971_TNBC_panCK_cnv) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Swarbrick_CID44971_TNBC_panCK_cnv <-t(Swarbrick_CID44971_TNBC_panCK_cnv)
write.csv(Swarbrick_CID44971_TNBC_panCK_cnv, file = "/R/R_Swarbrick/Swarbrick_Expression_Matrices/Swarbrick_Expression_Matrices_Output/Swarbrick_CID44971_TNBC_panCK_cnv_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Swarbrick_CID44971_TNBC_panCK_cnv)
gc()


#####################################
# Swarbrick_CID44991_TNBC_panCK_cnv #
#####################################
# Load Data
Swarbrick_CID44991_TNBC_panCK_cnv <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_inferCNV/Swarbrick_CID44991_TNBC_panCK_cnv.rds")
# Extract Expression Matrix
Swarbrick_CID44991_TNBC_panCK_cnv <- Swarbrick_CID44991_TNBC_panCK_cnv[["RNA"]]$counts
Swarbrick_CID44991_TNBC_panCK_cnv <- as.matrix(Swarbrick_CID44991_TNBC_panCK_cnv, 'sparseMatrix')
Swarbrick_CID44991_TNBC_panCK_cnv <- subset(Swarbrick_CID44991_TNBC_panCK_cnv, rownames(Swarbrick_CID44991_TNBC_panCK_cnv) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Swarbrick_CID44991_TNBC_panCK_cnv <-t(Swarbrick_CID44991_TNBC_panCK_cnv)
write.csv(Swarbrick_CID44991_TNBC_panCK_cnv, file = "/R/R_Swarbrick/Swarbrick_Expression_Matrices/Swarbrick_Expression_Matrices_Output/Swarbrick_CID44991_TNBC_panCK_cnv_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Swarbrick_CID44991_TNBC_panCK_cnv)
gc()




#############################################################################
# Extract Immune Mimicry Expression Matrices for Merged Data with All Genes #
#############################################################################

#########################################
# all_Swarbrick_breast_tumors_panCK_cnv #
#########################################
# Load data file
all_Swarbrick_breast_tumors_panCK_cnv <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_Merged/all_Swarbrick_breast_tumors_panCK_cnv.rds")
# Extract Expression Matrix
all_Swarbrick_breast_tumors_panCK_cnv <- all_Swarbrick_breast_tumors_panCK_cnv[["RNA"]]$counts
all_Swarbrick_breast_tumors_panCK_cnv <- as.matrix(all_Swarbrick_breast_tumors_panCK_cnv, 'sparseMatrix')
all_Swarbrick_breast_tumors_panCK_cnv <-t(all_Swarbrick_breast_tumors_panCK_cnv)
write.csv(all_Swarbrick_breast_tumors_panCK_cnv, file = "/R/R_Swarbrick/Swarbrick_Expression_Matrices/Swarbrick_Expression_Matrices_Output/all_Swarbrick_breast_tumors_panCK_cnv_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(all_Swarbrick_breast_tumors_panCK_cnv)
gc()


#####################################################
# Extract Expression Matrices per Sample: All Genes #
#####################################################

##############################
# Swarbrick_CID3941_ER_panCK #
##############################
# Load Data
Swarbrick_CID3941_ER_panCK <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_panCK/Swarbrick_CID3941_ER_panCK.rds")
# Extract Expression Matrix
Swarbrick_CID3941_ER_panCK <- Swarbrick_CID3941_ER_panCK[["RNA"]]$counts
Swarbrick_CID3941_ER_panCK <- as.matrix(Swarbrick_CID3941_ER_panCK, 'sparseMatrix')
Swarbrick_CID3941_ER_panCK <-t(Swarbrick_CID3941_ER_panCK)
write.csv(Swarbrick_CID3941_ER_panCK, file = "/R/R_Swarbrick/Swarbrick_Expression_Matrices/Swarbrick_Expression_Matrices_Output/Swarbrick_CID3941_ER_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Swarbrick_CID3941_ER_panCK)
gc()


##############################
# Swarbrick_CID3948_ER_panCK #
##############################
# Load Data
Swarbrick_CID3948_ER_panCK <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_panCK/Swarbrick_CID3948_ER_panCK.rds")
# Extract Expression Matrix
Swarbrick_CID3948_ER_panCK <- Swarbrick_CID3948_ER_panCK[["RNA"]]$counts
Swarbrick_CID3948_ER_panCK <- as.matrix(Swarbrick_CID3948_ER_panCK, 'sparseMatrix')
Swarbrick_CID3948_ER_panCK <-t(Swarbrick_CID3948_ER_panCK)
write.csv(Swarbrick_CID3948_ER_panCK, file = "/R/R_Swarbrick/Swarbrick_Expression_Matrices/Swarbrick_Expression_Matrices_Output/Swarbrick_CID3948_ER_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Swarbrick_CID3948_ER_panCK)
gc()


##############################
# Swarbrick_CID4040_ER_panCK #
##############################
# Load Data
Swarbrick_CID4040_ER_panCK <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_panCK/Swarbrick_CID4040_ER_panCK.rds")
# Extract Expression Matrix
Swarbrick_CID4040_ER_panCK <- Swarbrick_CID4040_ER_panCK[["RNA"]]$counts
Swarbrick_CID4040_ER_panCK <- as.matrix(Swarbrick_CID4040_ER_panCK, 'sparseMatrix')
Swarbrick_CID4040_ER_panCK <-t(Swarbrick_CID4040_ER_panCK)
write.csv(Swarbrick_CID4040_ER_panCK, file = "/R/R_Swarbrick/Swarbrick_Expression_Matrices/Swarbrick_Expression_Matrices_Output/Swarbrick_CID4040_ER_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Swarbrick_CID4040_ER_panCK)
gc()


##############################
# Swarbrick_CID4067_ER_panCK #
##############################
# Load Data
Swarbrick_CID4067_ER_panCK <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_panCK/Swarbrick_CID4067_ER_panCK.rds")
# Extract Expression Matrix
Swarbrick_CID4067_ER_panCK <- Swarbrick_CID4067_ER_panCK[["RNA"]]$counts
Swarbrick_CID4067_ER_panCK <- as.matrix(Swarbrick_CID4067_ER_panCK, 'sparseMatrix')
Swarbrick_CID4067_ER_panCK <-t(Swarbrick_CID4067_ER_panCK)
write.csv(Swarbrick_CID4067_ER_panCK, file = "/R/R_Swarbrick/Swarbrick_Expression_Matrices/Swarbrick_Expression_Matrices_Output/Swarbrick_CID4067_ER_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Swarbrick_CID4067_ER_panCK)
gc()


##############################
# Swarbrick_CID4290A_ER_panCK #
##############################
# Load Data
Swarbrick_CID4290A_ER_panCK <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_panCK/Swarbrick_CID4290A_ER_panCK.rds")
# Extract Expression Matrix
Swarbrick_CID4290A_ER_panCK <- Swarbrick_CID4290A_ER_panCK[["RNA"]]$counts
Swarbrick_CID4290A_ER_panCK <- as.matrix(Swarbrick_CID4290A_ER_panCK, 'sparseMatrix')
Swarbrick_CID4290A_ER_panCK <-t(Swarbrick_CID4290A_ER_panCK)
write.csv(Swarbrick_CID4290A_ER_panCK, file = "/R/R_Swarbrick/Swarbrick_Expression_Matrices/Swarbrick_Expression_Matrices_Output/Swarbrick_CID4290A_ER_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Swarbrick_CID4290A_ER_panCK)
gc()


##############################
# Swarbrick_CID4398_ER_panCK #
##############################
# Load Data
Swarbrick_CID4398_ER_panCK <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_panCK/Swarbrick_CID4398_ER_panCK.rds")
# Extract Expression Matrix
Swarbrick_CID4398_ER_panCK <- Swarbrick_CID4398_ER_panCK[["RNA"]]$counts
Swarbrick_CID4398_ER_panCK <- as.matrix(Swarbrick_CID4398_ER_panCK, 'sparseMatrix')
Swarbrick_CID4398_ER_panCK <-t(Swarbrick_CID4398_ER_panCK)
write.csv(Swarbrick_CID4398_ER_panCK, file = "/R/R_Swarbrick/Swarbrick_Expression_Matrices/Swarbrick_Expression_Matrices_Output/Swarbrick_CID4398_ER_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Swarbrick_CID4398_ER_panCK)
gc()


##############################
# Swarbrick_CID4461_ER_panCK #
##############################
# Load Data
Swarbrick_CID4461_ER_panCK <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_panCK/Swarbrick_CID4461_ER_panCK.rds")
# Extract Expression Matrix
Swarbrick_CID4461_ER_panCK <- Swarbrick_CID4461_ER_panCK[["RNA"]]$counts
Swarbrick_CID4461_ER_panCK <- as.matrix(Swarbrick_CID4461_ER_panCK, 'sparseMatrix')
Swarbrick_CID4461_ER_panCK <-t(Swarbrick_CID4461_ER_panCK)
write.csv(Swarbrick_CID4461_ER_panCK, file = "/R/R_Swarbrick/Swarbrick_Expression_Matrices/Swarbrick_Expression_Matrices_Output/Swarbrick_CID4461_ER_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Swarbrick_CID4461_ER_panCK)
gc()


##############################
# Swarbrick_CID4463_ER_panCK #
##############################
# Load Data
Swarbrick_CID4463_ER_panCK <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_panCK/Swarbrick_CID4463_ER_panCK.rds")
# Extract Expression Matrix
Swarbrick_CID4463_ER_panCK <- Swarbrick_CID4463_ER_panCK[["RNA"]]$counts
Swarbrick_CID4463_ER_panCK <- as.matrix(Swarbrick_CID4463_ER_panCK, 'sparseMatrix')
Swarbrick_CID4463_ER_panCK <-t(Swarbrick_CID4463_ER_panCK)
write.csv(Swarbrick_CID4463_ER_panCK, file = "/R/R_Swarbrick/Swarbrick_Expression_Matrices/Swarbrick_Expression_Matrices_Output/Swarbrick_CID4463_ER_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Swarbrick_CID4463_ER_panCK)
gc()


##############################
# Swarbrick_CID4471_ER_panCK #
##############################
# Load Data
Swarbrick_CID4471_ER_panCK <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_panCK/Swarbrick_CID4471_ER_panCK.rds")
# Extract Expression Matrix
Swarbrick_CID4471_ER_panCK <- Swarbrick_CID4471_ER_panCK[["RNA"]]$counts
Swarbrick_CID4471_ER_panCK <- as.matrix(Swarbrick_CID4471_ER_panCK, 'sparseMatrix')
Swarbrick_CID4471_ER_panCK <-t(Swarbrick_CID4471_ER_panCK)
write.csv(Swarbrick_CID4471_ER_panCK, file = "/R/R_Swarbrick/Swarbrick_Expression_Matrices/Swarbrick_Expression_Matrices_Output/Swarbrick_CID4471_ER_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Swarbrick_CID4471_ER_panCK)
gc()


###############################
# Swarbrick_CID4530N_ER_panCK #
###############################
# Load Data
Swarbrick_CID4530N_ER_panCK <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_panCK/Swarbrick_CID4530N_ER_panCK.rds")
# Extract Expression Matrix
Swarbrick_CID4530N_ER_panCK <- Swarbrick_CID4530N_ER_panCK[["RNA"]]$counts
Swarbrick_CID4530N_ER_panCK <- as.matrix(Swarbrick_CID4530N_ER_panCK, 'sparseMatrix')
Swarbrick_CID4530N_ER_panCK <-t(Swarbrick_CID4530N_ER_panCK)
write.csv(Swarbrick_CID4530N_ER_panCK, file = "/R/R_Swarbrick/Swarbrick_Expression_Matrices/Swarbrick_Expression_Matrices_Output/Swarbrick_CID4530N_ER_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Swarbrick_CID4530N_ER_panCK)
gc()


##############################
# Swarbrick_CID4535_ER_panCK #
##############################
# Load Data
Swarbrick_CID4535_ER_panCK <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_panCK/Swarbrick_CID4535_ER_panCK.rds")
# Extract Expression Matrix
Swarbrick_CID4535_ER_panCK <- Swarbrick_CID4535_ER_panCK[["RNA"]]$counts
Swarbrick_CID4535_ER_panCK <- as.matrix(Swarbrick_CID4535_ER_panCK, 'sparseMatrix')
Swarbrick_CID4535_ER_panCK <-t(Swarbrick_CID4535_ER_panCK)
write.csv(Swarbrick_CID4535_ER_panCK, file = "/R/R_Swarbrick/Swarbrick_Expression_Matrices/Swarbrick_Expression_Matrices_Output/Swarbrick_CID4535_ER_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Swarbrick_CID4535_ER_panCK)
gc()


################################
# Swarbrick_CID3586_HER2_panCK #
################################
# Load Data
Swarbrick_CID3586_HER2_panCK <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_panCK/Swarbrick_CID3586_HER2_panCK.rds")
# Extract Expression Matrix
Swarbrick_CID3586_HER2_panCK <- Swarbrick_CID3586_HER2_panCK[["RNA"]]$counts
Swarbrick_CID3586_HER2_panCK <- as.matrix(Swarbrick_CID3586_HER2_panCK, 'sparseMatrix')
Swarbrick_CID3586_HER2_panCK <-t(Swarbrick_CID3586_HER2_panCK)
write.csv(Swarbrick_CID3586_HER2_panCK, file = "/R/R_Swarbrick/Swarbrick_Expression_Matrices/Swarbrick_Expression_Matrices_Output/Swarbrick_CID3586_HER2_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Swarbrick_CID3586_HER2_panCK)
gc()


################################
# Swarbrick_CID3838_HER2_panCK #
################################
# Load Data
Swarbrick_CID3838_HER2_panCK <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_panCK/Swarbrick_CID3838_HER2_panCK.rds")
# Extract Expression Matrix
Swarbrick_CID3838_HER2_panCK <- Swarbrick_CID3838_HER2_panCK[["RNA"]]$counts
Swarbrick_CID3838_HER2_panCK <- as.matrix(Swarbrick_CID3838_HER2_panCK, 'sparseMatrix')
Swarbrick_CID3838_HER2_panCK <-t(Swarbrick_CID3838_HER2_panCK)
write.csv(Swarbrick_CID3838_HER2_panCK, file = "/R/R_Swarbrick/Swarbrick_Expression_Matrices/Swarbrick_Expression_Matrices_Output/Swarbrick_CID3838_HER2_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Swarbrick_CID3838_HER2_panCK)
gc()


################################
# Swarbrick_CID3921_HER2_panCK #
################################
# Load Data
Swarbrick_CID3921_HER2_panCK <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_panCK/Swarbrick_CID3921_HER2_panCK.rds")
# Extract Expression Matrix
Swarbrick_CID3921_HER2_panCK <- Swarbrick_CID3921_HER2_panCK[["RNA"]]$counts
Swarbrick_CID3921_HER2_panCK <- as.matrix(Swarbrick_CID3921_HER2_panCK, 'sparseMatrix')
Swarbrick_CID3921_HER2_panCK <-t(Swarbrick_CID3921_HER2_panCK)
write.csv(Swarbrick_CID3921_HER2_panCK, file = "/R/R_Swarbrick/Swarbrick_Expression_Matrices/Swarbrick_Expression_Matrices_Output/Swarbrick_CID3921_HER2_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Swarbrick_CID3921_HER2_panCK)
gc()


################################
# Swarbrick_CID4066_HER2_panCK #
################################
# Load Data
Swarbrick_CID4066_HER2_panCK <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_panCK/Swarbrick_CID4066_HER2_panCK.rds")
# Extract Expression Matrix
Swarbrick_CID4066_HER2_panCK <- Swarbrick_CID4066_HER2_panCK[["RNA"]]$counts
Swarbrick_CID4066_HER2_panCK <- as.matrix(Swarbrick_CID4066_HER2_panCK, 'sparseMatrix')
Swarbrick_CID4066_HER2_panCK <-t(Swarbrick_CID4066_HER2_panCK)
write.csv(Swarbrick_CID4066_HER2_panCK, file = "/R/R_Swarbrick/Swarbrick_Expression_Matrices/Swarbrick_Expression_Matrices_Output/Swarbrick_CID4066_HER2_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Swarbrick_CID4066_HER2_panCK)
gc()


#################################
# Swarbrick_CID45171_HER2_panCK #
#################################
# Load Data
Swarbrick_CID45171_HER2_panCK <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_panCK/Swarbrick_CID45171_HER2_panCK.rds")
# Extract Expression Matrix
Swarbrick_CID45171_HER2_panCK <- Swarbrick_CID45171_HER2_panCK[["RNA"]]$counts
Swarbrick_CID45171_HER2_panCK <- as.matrix(Swarbrick_CID45171_HER2_panCK, 'sparseMatrix')
Swarbrick_CID45171_HER2_panCK <-t(Swarbrick_CID45171_HER2_panCK)
write.csv(Swarbrick_CID45171_HER2_panCK, file = "/R/R_Swarbrick/Swarbrick_Expression_Matrices/Swarbrick_Expression_Matrices_Output/Swarbrick_CID45171_HER2_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Swarbrick_CID45171_HER2_panCK)
gc()


################################
# Swarbrick_CID3946_TNBC_panCK #
################################
# Load Data
Swarbrick_CID3946_TNBC_panCK <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_panCK/Swarbrick_CID3946_TNBC_panCK.rds")
# Extract Expression Matrix
Swarbrick_CID3946_TNBC_panCK <- Swarbrick_CID3946_TNBC_panCK[["RNA"]]$counts
Swarbrick_CID3946_TNBC_panCK <- as.matrix(Swarbrick_CID3946_TNBC_panCK, 'sparseMatrix')
Swarbrick_CID3946_TNBC_panCK <-t(Swarbrick_CID3946_TNBC_panCK)
write.csv(Swarbrick_CID3946_TNBC_panCK, file = "/R/R_Swarbrick/Swarbrick_Expression_Matrices/Swarbrick_Expression_Matrices_Output/Swarbrick_CID3946_TNBC_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Swarbrick_CID3946_TNBC_panCK)
gc()


################################
# Swarbrick_CID3963_TNBC_panCK #
################################
# Load Data
Swarbrick_CID3963_TNBC_panCK <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_panCK/Swarbrick_CID3963_TNBC_panCK.rds")
# Extract Expression Matrix
Swarbrick_CID3963_TNBC_panCK <- Swarbrick_CID3963_TNBC_panCK[["RNA"]]$counts
Swarbrick_CID3963_TNBC_panCK <- as.matrix(Swarbrick_CID3963_TNBC_panCK, 'sparseMatrix')
Swarbrick_CID3963_TNBC_panCK <-t(Swarbrick_CID3963_TNBC_panCK)
write.csv(Swarbrick_CID3963_TNBC_panCK, file = "/R/R_Swarbrick/Swarbrick_Expression_Matrices/Swarbrick_Expression_Matrices_Output/Swarbrick_CID3963_TNBC_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Swarbrick_CID3963_TNBC_panCK)
gc()


################################
# Swarbrick_CID4465_TNBC_panCK #
################################
# Load Data
Swarbrick_CID4465_TNBC_panCK <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_panCK/Swarbrick_CID4465_TNBC_panCK.rds")
# Extract Expression Matrix
Swarbrick_CID4465_TNBC_panCK <- Swarbrick_CID4465_TNBC_panCK[["RNA"]]$counts
Swarbrick_CID4465_TNBC_panCK <- as.matrix(Swarbrick_CID4465_TNBC_panCK, 'sparseMatrix')
Swarbrick_CID4465_TNBC_panCK <-t(Swarbrick_CID4465_TNBC_panCK)
write.csv(Swarbrick_CID4465_TNBC_panCK, file = "/R/R_Swarbrick/Swarbrick_Expression_Matrices/Swarbrick_Expression_Matrices_Output/Swarbrick_CID4465_TNBC_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Swarbrick_CID4465_TNBC_panCK)
gc()


################################
# Swarbrick_CID4495_TNBC_panCK #
################################
# Load Data
Swarbrick_CID4495_TNBC_panCK <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_panCK/Swarbrick_CID4495_TNBC_panCK.rds")
# Extract Expression Matrix
Swarbrick_CID4495_TNBC_panCK <- Swarbrick_CID4495_TNBC_panCK[["RNA"]]$counts
Swarbrick_CID4495_TNBC_panCK <- as.matrix(Swarbrick_CID4495_TNBC_panCK, 'sparseMatrix')
Swarbrick_CID4495_TNBC_panCK <-t(Swarbrick_CID4495_TNBC_panCK)
write.csv(Swarbrick_CID4495_TNBC_panCK, file = "/R/R_Swarbrick/Swarbrick_Expression_Matrices/Swarbrick_Expression_Matrices_Output/Swarbrick_CID4495_TNBC_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Swarbrick_CID4495_TNBC_panCK)
gc()


################################
# Swarbrick_CID4513_TNBC_panCK #
################################
# Load Data
Swarbrick_CID4513_TNBC_panCK <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_panCK/Swarbrick_CID4513_TNBC_panCK.rds")
# Extract Expression Matrix
Swarbrick_CID4513_TNBC_panCK <- Swarbrick_CID4513_TNBC_panCK[["RNA"]]$counts
Swarbrick_CID4513_TNBC_panCK <- as.matrix(Swarbrick_CID4513_TNBC_panCK, 'sparseMatrix')
Swarbrick_CID4513_TNBC_panCK <-t(Swarbrick_CID4513_TNBC_panCK)
write.csv(Swarbrick_CID4513_TNBC_panCK, file = "/R/R_Swarbrick/Swarbrick_Expression_Matrices/Swarbrick_Expression_Matrices_Output/Swarbrick_CID4513_TNBC_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Swarbrick_CID4513_TNBC_panCK)
gc()


################################
# Swarbrick_CID4515_TNBC_panCK #
################################
# Load Data
Swarbrick_CID4515_TNBC_panCK <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_panCK/Swarbrick_CID4515_TNBC_panCK.rds")
# Extract Expression Matrix
Swarbrick_CID4515_TNBC_panCK <- Swarbrick_CID4515_TNBC_panCK[["RNA"]]$counts
Swarbrick_CID4515_TNBC_panCK <- as.matrix(Swarbrick_CID4515_TNBC_panCK, 'sparseMatrix')
Swarbrick_CID4515_TNBC_panCK <-t(Swarbrick_CID4515_TNBC_panCK)
write.csv(Swarbrick_CID4515_TNBC_panCK, file = "/R/R_Swarbrick/Swarbrick_Expression_Matrices/Swarbrick_Expression_Matrices_Output/Swarbrick_CID4515_TNBC_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Swarbrick_CID4515_TNBC_panCK)
gc()


################################
# Swarbrick_CID4523_TNBC_panCK #
################################
# Load Data
Swarbrick_CID4523_TNBC_panCK <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_panCK/Swarbrick_CID4523_TNBC_panCK.rds")
# Extract Expression Matrix
Swarbrick_CID4523_TNBC_panCK <- Swarbrick_CID4523_TNBC_panCK[["RNA"]]$counts
Swarbrick_CID4523_TNBC_panCK <- as.matrix(Swarbrick_CID4523_TNBC_panCK, 'sparseMatrix')
Swarbrick_CID4523_TNBC_panCK <-t(Swarbrick_CID4523_TNBC_panCK)
write.csv(Swarbrick_CID4523_TNBC_panCK, file = "/R/R_Swarbrick/Swarbrick_Expression_Matrices/Swarbrick_Expression_Matrices_Output/Swarbrick_CID4523_TNBC_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Swarbrick_CID4523_TNBC_panCK)
gc()


#################################
# Swarbrick_CID44041_TNBC_panCK #
#################################
# Load Data
Swarbrick_CID44041_TNBC_panCK <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_panCK/Swarbrick_CID44041_TNBC_panCK.rds")
# Extract Expression Matrix
Swarbrick_CID44041_TNBC_panCK <- Swarbrick_CID44041_TNBC_panCK[["RNA"]]$counts
Swarbrick_CID44041_TNBC_panCK <- as.matrix(Swarbrick_CID44041_TNBC_panCK, 'sparseMatrix')
Swarbrick_CID44041_TNBC_panCK <-t(Swarbrick_CID44041_TNBC_panCK)
write.csv(Swarbrick_CID44041_TNBC_panCK, file = "/R/R_Swarbrick/Swarbrick_Expression_Matrices/Swarbrick_Expression_Matrices_Output/Swarbrick_CID44041_TNBC_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Swarbrick_CID44041_TNBC_panCK)
gc()


#################################
# Swarbrick_CID44971_TNBC_panCK #
#################################
# Load Data
Swarbrick_CID44971_TNBC_panCK <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_panCK/Swarbrick_CID44971_TNBC_panCK.rds")
# Extract Expression Matrix
Swarbrick_CID44971_TNBC_panCK <- Swarbrick_CID44971_TNBC_panCK[["RNA"]]$counts
Swarbrick_CID44971_TNBC_panCK <- as.matrix(Swarbrick_CID44971_TNBC_panCK, 'sparseMatrix')
Swarbrick_CID44971_TNBC_panCK <-t(Swarbrick_CID44971_TNBC_panCK)
write.csv(Swarbrick_CID44971_TNBC_panCK, file = "/R/R_Swarbrick/Swarbrick_Expression_Matrices/Swarbrick_Expression_Matrices_Output/Swarbrick_CID44971_TNBC_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Swarbrick_CID44971_TNBC_panCK)
gc()


#################################
# Swarbrick_CID44991_TNBC_panCK #
#################################
# Load Data
Swarbrick_CID44991_TNBC_panCK <- readRDS(file = "/R/R_Swarbrick/Swarbrick_RDS/RDS_panCK/Swarbrick_CID44991_TNBC_panCK.rds")
# Extract Expression Matrix
Swarbrick_CID44991_TNBC_panCK <- Swarbrick_CID44991_TNBC_panCK[["RNA"]]$counts
Swarbrick_CID44991_TNBC_panCK <- as.matrix(Swarbrick_CID44991_TNBC_panCK, 'sparseMatrix')
Swarbrick_CID44991_TNBC_panCK <-t(Swarbrick_CID44991_TNBC_panCK)
write.csv(Swarbrick_CID44991_TNBC_panCK, file = "/R/R_Swarbrick/Swarbrick_Expression_Matrices/Swarbrick_Expression_Matrices_Output/Swarbrick_CID44991_TNBC_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Swarbrick_CID44991_TNBC_panCK)
gc()











