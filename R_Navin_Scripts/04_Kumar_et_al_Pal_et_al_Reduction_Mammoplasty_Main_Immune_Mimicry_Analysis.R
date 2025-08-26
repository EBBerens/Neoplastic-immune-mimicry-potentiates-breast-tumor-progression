#####################################################################################################################
#        Navin (Kumar et al) and Visvader (Pal et al) Datasets Reduction Mammoplasty Analysis Steps                 #
#####################################################################################################################
# Step 7: Determine Cluster Markers & Save Summary Information                                                      #
# Step 8: Shannon-Weaver Diversity Analysis                                                                         #
# Step 9: Performing Fgsea to Identify Biological Processes Enriched in Shared Clusters                             #
# Step 10: Annotating Diverse Clusters Based on Differentially Enriched Biological Processes                        #
# Step 11: Compare Signatures in Immune-like vs Non-Immune-like Neoplastic Cells                                    #
# Step 12: Generate Expression Matrices for Each Sample: All Genes vs Immune Mimicry Markers                        #
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
Navin_Visvader_NORM_panCK_combined <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Merged/Navin_Visvader_NORM_panCK_combined.rds")

# Perform Marker Analysis
Navin_Visvader_NORM_panCK_combined.cluster.markers <- FindAllMarkers(Navin_Visvader_NORM_panCK_combined, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(Navin_Visvader_NORM_panCK_combined.cluster.markers, file = "/R/R_Navin/Navin_Output/Navin_Visvader_NORM_panCK_combined.cluster.markers.csv")

# Save Example Plots
DimPlot(Navin_Visvader_NORM_panCK_combined, reduction = "umap", raster = FALSE)
ggsave("/R/R_Navin/Navin_Output/Navin_Visvader_Combined_all_KRTpos_cells_clusters_UMAP.tiff", plot = last_plot(), device = "tiff",
       scale = 1, width = 16, height = 10,
       dpi = 200, limitsize = TRUE)

DimPlot(Navin_Visvader_NORM_panCK_combined, reduction = "umap", raster = FALSE, label = TRUE)
ggsave("/R/R_Navin/Navin_Output/Navin_Visvader_Combined_all_KRTpos_cells_clusters_UMAP_labels.tiff", plot = last_plot(), device = "tiff",
       scale = 1, width = 16, height = 10,
       dpi = 200, limitsize = TRUE)

DimPlot(Navin_Visvader_NORM_panCK_combined, reduction = "umap", group.by = "orig.ident", raster = FALSE)
ggsave("/R/R_Navin/Navin_Output/Navin_Visvader_Combined_all_KRTpos_cells_clusters_orig.ident_UMAP.tiff", plot = last_plot(), device = "tiff",
       scale = 1, width = 16, height = 10,
       dpi = 200, limitsize = TRUE)

DimPlot(Navin_Visvader_NORM_panCK_combined, reduction = "umap", group.by = "orig.ident", raster = FALSE) + NoLegend()
ggsave("/R/R_Navin/Navin_Output/Navin_Visvader_Combined_all_KRTpos_cells_clusters_orig.ident_UMAP_No_Legend.tiff", plot = last_plot(), device = "tiff",
       scale = 1, width = 16, height = 10,
       dpi = 200, limitsize = TRUE)

FeaturePlot(Navin_Visvader_NORM_panCK_combined, c("KRT14"),  raster = FALSE, cols = c("grey90", "firebrick"))
ggsave("/R/R_Navin/Navin_Output/Navin_Visvader_Combined_breast_tumors_Neoplastic_KRT14_UMAP_Custom_Colors.tiff", plot = last_plot(), device = "tiff", 
       scale = 1, width = 16, height = 10,
       dpi = 200, limitsize = TRUE)

FeaturePlot(Navin_Visvader_NORM_panCK_combined, c("KRT18"),  raster = FALSE, cols = c("grey90", "firebrick"))
ggsave("/R/R_Navin/Navin_Output/Navin_Visvader_Combined_breast_tumors_Neoplastic_KRT18_UMAP_Custom_Colors.tiff", plot = last_plot(), device = "tiff", 
       scale = 1, width = 16, height = 10,
       dpi = 200, limitsize = TRUE)

FeaturePlot(Navin_Visvader_NORM_panCK_combined, c("KRT19"),  raster = FALSE, cols = c("grey90", "firebrick"))
ggsave("/R/R_Navin/Navin_Output/Navin_Visvader_Combined_breast_tumors_Neoplastic_KRT19_UMAP_Custom_Colors.tiff", plot = last_plot(), device = "tiff", 
       scale = 1, width = 16, height = 10,
       dpi = 200, limitsize = TRUE)

FeaturePlot(Navin_Visvader_NORM_panCK_combined, c("CD69"),  raster = FALSE, cols = c("grey90", "firebrick"))
ggsave("/R/R_Navin/Navin_Output/Navin_Visvader_Combined_breast_tumors_Neoplastic_CD69_UMAP_Custom_Colors.tiff", plot = last_plot(), device = "tiff", 
       scale = 1, width = 16, height = 10,
       dpi = 200, limitsize = TRUE)

FeaturePlot(Navin_Visvader_NORM_panCK_combined, c("PTPRC"),  raster = FALSE, cols = c("grey90", "firebrick"))
ggsave("/R/R_Navin/Navin_Output/Navin_Visvader_Combined_breast_tumors_Neoplastic_PTPRC_UMAP_Custom_Colors.tiff", plot = last_plot(), device = "tiff", 
       scale = 1, width = 16, height = 10,
       dpi = 200, limitsize = TRUE)

# Count Total Cells 
Total <- Idents(Navin_Visvader_NORM_panCK_combined)
# Total = 404,201 cells from 124 samples



#############################################
# Step 8: Shannon-Weaver Diversity Analysis #            
#############################################

###################################################
# Calculate Number of Cells per Ident per Cluster #
# for downstream diversity calculations ###########
###################################################
# Store cluster identities in object@meta.data$cluster.designation
Navin_Visvader_NORM_panCK_combined[["cluster.designation"]] <- Idents(object = Navin_Visvader_NORM_panCK_combined)

# Get number of cells per cluster and per sample of origin
Navin_Visvader_NORM_panCK_combined.idents.per.cluster <- table(Navin_Visvader_NORM_panCK_combined@meta.data$cluster.designation, Navin_Visvader_NORM_panCK_combined@meta.data$orig.ident)

# Save table
write.csv(Navin_Visvader_NORM_panCK_combined.idents.per.cluster, file = "/R/R_Navin/Navin_Output/Navin_Visvader_NORM_panCK_combined.idents.per.cluster.csv")

###########################################
# Begin Shannon-Weaver Diversity Analysis #
###########################################
library(vegan)

# Load Table
ident.cluster.table <- read.csv(file = "/R/R_Navin/Navin_Output/Navin_Visvader_NORM_panCK_combined.idents.per.cluster.csv", header = TRUE, row.names = 1)

# Calculate Shannon index
shannon_diversity <- diversity(ident.cluster.table)  

# Save results
write.csv(shannon_diversity, file = "/R/R_Navin/Navin_Output/Navin_Visvader_NORM_panCK_combined.idents.per.cluster.shannon_diversity.csv")

# Remove objects
rm(Navin_Visvader_NORM_panCK_combined.idents.per.cluster)
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

# All clusters are diverse in merged mammoplasties, so specifically running fgsea on clusters with immune receptor expression

#################################
# Clusters 7 & 9 & 12 for FGSEA #
#################################

# Load Cluster Markers
Navin_Visvader_NORM_panCK_combined.cluster.markers <- read.csv(file = "/R/R_Navin/Navin_Output/Navin_Visvader_NORM_panCK_combined.cluster.markers.csv")

# Remove genes with adjusted p-values above <0.01
Navin_Visvader_NORM_panCK_combined.cluster.markers <- Navin_Visvader_NORM_panCK_combined.cluster.markers %>% filter(p_val_adj<0.01)

# Collect columns: gene names, avg_log2FC, and cluster; sort by descending 
Navin_Visvader_NORM_panCK_combined.cluster.markers <- Navin_Visvader_NORM_panCK_combined.cluster.markers %>% dplyr::select(gene, avg_log2FC, cluster)
Navin_Visvader_NORM_panCK_combined.cluster.markers <- Navin_Visvader_NORM_panCK_combined.cluster.markers[order(Navin_Visvader_NORM_panCK_combined.cluster.markers$avg_log2FC, decreasing = TRUE),]  

# Subset highly diverse clusters & filter on gene name 
Cluster7 <- Navin_Visvader_NORM_panCK_combined.cluster.markers %>% filter(cluster == 7) 
Cluster9 <- Navin_Visvader_NORM_panCK_combined.cluster.markers %>% filter(cluster == 9) 
Cluster12 <- Navin_Visvader_NORM_panCK_combined.cluster.markers %>% filter(cluster == 12) 

# Remove marker table
rm(Navin_Visvader_NORM_panCK_combined.cluster.markers)

# Save Output
write.csv(Cluster7, file = "/R/R_Navin/Navin_fgsea/Navin_fgsea_Output/Cluster7.csv")
write.csv(Cluster9, file = "/R/R_Navin/Navin_fgsea/Navin_fgsea_Output/Cluster9.csv")
write.csv(Cluster12, file = "/R/R_Navin/Navin_fgsea/Navin_fgsea_Output/Cluster12.csv")

# Collapse into gene name and LogFC
Cluster7 <- Cluster7 %>% dplyr::select(gene, avg_log2FC)
Cluster9 <- Cluster9 %>% dplyr::select(gene, avg_log2FC)
Cluster12 <- Cluster12 %>% dplyr::select(gene, avg_log2FC)

# Create a vector of gene ranks
Cluster7.ranks <- deframe(Cluster7)
Cluster9.ranks <- deframe(Cluster9)
Cluster12.ranks <- deframe(Cluster12)

# Remove original objects
rm(Cluster7)
rm(Cluster9)
rm(Cluster12)


# Load the pathways into a named list
pathways.GO.BP <- gmtPathways("/R/R_Navin/Navin_fgsea/Navin_fgsea_MSigDB/c5.go.bp.v2022.1.Hs.symbols.gmt")


##########################
# Run fgsea on Cluster7  #
##########################
fgsea.Cluster7.res <- fgsea(pathways = pathways.GO.BP, 
                            stats    = Cluster7.ranks,
                            eps      = 0.0,
                            minSize  = 15,
                            maxSize  = 500, nPermSimple = 100000)

# Summarize Cluster7.results by collapsing into main pathways, to omit redundancies
Cluster7.collapsedPathways <- collapsePathways(fgsea.Cluster7.res[order(pval)][padj < 0.05], 
                                               pathways.GO.BP, Cluster7.ranks)
Cluster7.mainPathways <- fgsea.Cluster7.res[pathway %in% Cluster7.collapsedPathways$mainPathways][
  order(-NES), pathway]
plotGseaTable(pathways.GO.BP[Cluster7.mainPathways], Cluster7.ranks, fgsea.Cluster7.res, 
              gseaParam = 0.5)

# Filter Cluster7.results so they only include the main pathways
fgsea.Cluster7.res_collapsed <- filter(fgsea.Cluster7.res, pathway %in% Cluster7.mainPathways)

# Save Cluster7.results
fwrite(fgsea.Cluster7.res, file="/R/R_Navin/Navin_fgsea/Navin_fgsea_Output/fgsea.Cluster7.res_GO.BP.csv")
fwrite(fgsea.Cluster7.res_collapsed, file="/R/R_Navin/Navin_fgsea/Navin_fgsea_Output/fgsea.Cluster7.res_GO.BP_collapsed.csv")

# Retain pathways with adjusted p-values below <0.05 
fgsea.Cluster7.res_sig <- fgsea.Cluster7.res %>% filter(padj<0.05)

# Save Cluster7.results
fwrite(fgsea.Cluster7.res_sig, file="/R/R_Navin/Navin_fgsea/Navin_fgsea_Output/fgsea.Cluster7.res_sig.csv")


# Remove objects
rm(fgsea.Cluster7.res)
rm(fgsea.Cluster7.res_collapsed)
rm(Cluster7.collapsedPathways)
rm(Cluster7.mainPathways)
rm(Cluster7.ranks)
rm(fgsea.Cluster7.res_sig)


###########################
# Run fgsea on Cluster9   #
###########################
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
fwrite(fgsea.Cluster9.res, file="/R/R_Navin/Navin_fgsea/Navin_fgsea_Output/fgsea.Cluster9.res_GO.BP.csv")
fwrite(fgsea.Cluster9.res_collapsed, file="/R/R_Navin/Navin_fgsea/Navin_fgsea_Output/fgsea.Cluster9.res_GO.BP_collapsed.csv")

# Retain pathways with adjusted p-values below <0.05 
fgsea.Cluster9.res_sig <- fgsea.Cluster9.res %>% filter(padj<0.05)

# Save Cluster13.results
fwrite(fgsea.Cluster9.res_sig, file="/R/R_Navin/Navin_fgsea/Navin_fgsea_Output/fgsea.Cluster9.res_sig.csv")

# Remove objects
rm(fgsea.Cluster9.res)
rm(fgsea.Cluster9.res_collapsed)
rm(Cluster9.collapsedPathways)
rm(Cluster9.mainPathways)
rm(Cluster9.ranks)
rm(fgsea.Cluster9.res_sig)




############################
# Run fgsea on Cluster12   #
############################
fgsea.Cluster12.res <- fgsea(pathways = pathways.GO.BP, 
                            stats    = Cluster12.ranks,
                            eps      = 0.0,
                            minSize  = 15,
                            maxSize  = 500, nPermSimple = 100000)

# Summarize Cluster12.results by collapsing into main pathways, to omit redundancies
Cluster12.collapsedPathways <- collapsePathways(fgsea.Cluster12.res[order(pval)][padj < 0.05], 
                                               pathways.GO.BP, Cluster12.ranks)
Cluster12.mainPathways <- fgsea.Cluster12.res[pathway %in% Cluster12.collapsedPathways$mainPathways][
  order(-NES), pathway]
plotGseaTable(pathways.GO.BP[Cluster12.mainPathways], Cluster12.ranks, fgsea.Cluster12.res, 
              gseaParam = 0.5)

# Filter Cluster12.results so they only include the main pathways
fgsea.Cluster12.res_collapsed <- filter(fgsea.Cluster12.res, pathway %in% Cluster12.mainPathways)

# Save Cluster12.results
fwrite(fgsea.Cluster12.res, file="/R/R_Navin/Navin_fgsea/Navin_fgsea_Output/fgsea.Cluster12.res_GO.BP.csv")
fwrite(fgsea.Cluster12.res_collapsed, file="/R/R_Navin/Navin_fgsea/Navin_fgsea_Output/fgsea.Cluster12.res_GO.BP_collapsed.csv")

# Retain pathways with adjusted p-values below <0.05 
fgsea.Cluster12.res_sig <- fgsea.Cluster12.res %>% filter(padj<0.05)

# Save Cluster13.results
fwrite(fgsea.Cluster12.res_sig, file="/R/R_Navin/Navin_fgsea/Navin_fgsea_Output/fgsea.Cluster12.res_sig.csv")

# Remove objects
rm(fgsea.Cluster12.res)
rm(fgsea.Cluster12.res_collapsed)
rm(Cluster12.collapsedPathways)
rm(Cluster12.mainPathways)
rm(Cluster12.ranks)
rm(fgsea.Cluster12.res_sig)





##############################################################################################
# Step 10: Annotating Diverse Clusters Based on Differentially Enriched Biological Processes #
##############################################################################################  

################################################
# Mimicked Clusters                            #
# Cluster 7 = Immune (Lymphoid)                #
# Cluster 9 = Immune (Lymphoid)                #
# Cluster 12 = Immune (Lymphoid)               #
################################################

# Load Seurat Object
Navin_Visvader_NORM_panCK_combined <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Merged/Navin_Visvader_NORM_panCK_combined.rds")

# Visualize Clusters
DimPlot(Navin_Visvader_NORM_panCK_combined, reduction = "umap", raster = FALSE)
DimPlot(Navin_Visvader_NORM_panCK_combined, reduction = "umap", group.by = 'orig.ident', raster = FALSE)

# Annotate Clusters = detailed.ident
Mimicry_IDs_detailed <- c("Epithelial Cells", "Epithelial Cells", "Epithelial Cells", "Epithelial Cells", "Epithelial Cells",
                          "Epithelial Cells", "Epithelial Cells", "Lymphoid-like", "Epithelial Cells", "Lymphoid-like",
                          "Epithelial Cells", "Epithelial Cells", "Lymphoid-like", "Epithelial Cells")

# Apply New Labels to Clusters in Already Annotated Object
names(Mimicry_IDs_detailed) <- levels(Navin_Visvader_NORM_panCK_combined)
Navin_Visvader_NORM_panCK_combined_annotated <- RenameIdents(Navin_Visvader_NORM_panCK_combined, Mimicry_IDs_detailed)

# Remove old Seurat Object and ID list
rm(Mimicry_IDs_detailed)
rm(Navin_Visvader_NORM_panCK_combined)

# Store Annotations Under Classification
Navin_Visvader_NORM_panCK_combined_annotated[["detailed.ident"]] <- Idents(object = Navin_Visvader_NORM_panCK_combined_annotated)

# Draw plot with Changed colors
DimPlot(Navin_Visvader_NORM_panCK_combined_annotated, reduction = "umap", raster = FALSE, group.by = 'detailed.ident', cols = c('Epithelial Cells' = '#00A9FF', 'Lymphoid-like' = '#CD9600'))

# Save Plot
ggsave("/R/R_Navin/Navin_Output/Navin_Visvader_NORM_panCK_combined_detailed_annotation_with_custom_colors.tiff", plot = last_plot(), device = "tiff", 
       scale = 1, width = 16, height = 10,
       dpi = 200, limitsize = TRUE)

# Check Idents
Idents(object = Navin_Visvader_NORM_panCK_combined_annotated)
table(Idents(Navin_Visvader_NORM_panCK_combined_annotated))



# Annotate Clusters = main.ident
Mimicry_IDs_main <- c("Epithelial Cells", "Immune-like")

# Apply New Labels to Clusters
names(Mimicry_IDs_main) <- levels(Navin_Visvader_NORM_panCK_combined_annotated)
Navin_Visvader_NORM_panCK_combined_annotated <- RenameIdents(Navin_Visvader_NORM_panCK_combined_annotated, Mimicry_IDs_main)

# Remove old Seurat Object and ID list
rm(Mimicry_IDs_main)

# Store Annotations Under Classification
Navin_Visvader_NORM_panCK_combined_annotated[["main.ident"]] <- Idents(object = Navin_Visvader_NORM_panCK_combined_annotated)

# Draw plot with Changed colors
DimPlot(Navin_Visvader_NORM_panCK_combined_annotated, reduction = "umap", raster = FALSE, group.by = 'main.ident', cols = c('Epithelial Cells' = '#00A9FF', 'Immune-like' = '#F8766D'))

# Save Plot
ggsave("/R/R_Navin/Navin_Output/Navin_Visvader_NORM_panCK_combined_main_annotation_with_custom_colors.tiff", plot = last_plot(), device = "tiff", 
       scale = 1, width = 16, height = 10,
       dpi = 200, limitsize = TRUE)

# Check Idents
table(Idents(Navin_Visvader_NORM_panCK_combined_annotated))

# Set Identity for moving forward
Navin_Visvader_NORM_panCK_combined_annotated <- SetIdent(Navin_Visvader_NORM_panCK_combined_annotated, value = "detailed.ident")

# Verify Idents have changed
table(Idents(Navin_Visvader_NORM_panCK_combined_annotated))

# Save Annotated Seurat Object
saveRDS(Navin_Visvader_NORM_panCK_combined_annotated, file = "/R/R_Navin/Navin_RDS/RDS_Annotated/Navin_Visvader_NORM_panCK_combined_annotated.rds")

# Extract  Meta Data
Navin_Visvader_NORM_panCK_combined_annotated.meta.data <- as.data.frame(as.matrix(Navin_Visvader_NORM_panCK_combined_annotated@meta.data))

# Save Meta Data
write.csv(Navin_Visvader_NORM_panCK_combined_annotated.meta.data, file = "/R/R_Navin/Navin_Output/Navin_Visvader_NORM_panCK_combined_annotated.meta.data.csv")

# Remove Meta Data
rm(Navin_Visvader_NORM_panCK_combined_annotated.meta.data)

# Count Total Cells 
Total <- Idents(Navin_Visvader_NORM_panCK_combined_annotated)
# Total = 404,201 cells from 124 mammoplasties

#########################################
# Draw IM Gene Expression Feature Plots #
#########################################
# Navin_Visvader_NORM_panCK_combined_annotated <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_Annotated/Navin_Visvader_NORM_panCK_combined_annotated.rds")

# Use patchwork to combine plots and incorporate custom colors

p1 <- FeaturePlot(Navin_Visvader_NORM_panCK_combined_annotated, c("CD3D"), cols = c("grey90", "darkgoldenrod4"))
p2 <- FeaturePlot(Navin_Visvader_NORM_panCK_combined_annotated, c("CD52"), cols = c("grey90", "darkgoldenrod4"))
p3 <- FeaturePlot(Navin_Visvader_NORM_panCK_combined_annotated, c("CD69"), cols = c("grey90", "darkgoldenrod4"))
p4 <- FeaturePlot(Navin_Visvader_NORM_panCK_combined_annotated, c("CD14"), cols = c("grey90", "dodgerblue4"))
p5 <- FeaturePlot(Navin_Visvader_NORM_panCK_combined_annotated, c("CD68"), cols = c("grey90", "dodgerblue4"))
p6 <- FeaturePlot(Navin_Visvader_NORM_panCK_combined_annotated, c("CD74"), cols = c("grey90", "dodgerblue4"))
p7 <- FeaturePlot(Navin_Visvader_NORM_panCK_combined_annotated, c("PTPRC"), cols = c("grey90", "darkgreen"))
p8 <- FeaturePlot(Navin_Visvader_NORM_panCK_combined_annotated, c("CD44"), cols = c("grey90", "darkgreen"))
p9 <- FeaturePlot(Navin_Visvader_NORM_panCK_combined_annotated, c("CXCR4"), cols = c("grey90", "darkgreen"))

# Combine plots with patchwork
(p1 | p2 | p3) / (p4 | p5 | p6) / (p7 | p8 | p9)

ggsave("/R/R_Navin/Navin_Output/Navin_Visvader_NORM_panCK_combined_annotated_IM_Marker_Patchwork.tiff", plot = last_plot(), device = "tiff", 
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

##################################################################################
# Step 11: Compare Signatures in Immune-like vs Non-Immune-like Neoplastic Cells #
##################################################################################

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
mimicry_order1 <- c("Epithelial Cells", "Immune-like")
Navin_Visvader_NORM_panCK_combined_annotated[["main.ident"]] <- factor(x = Navin_Visvader_NORM_panCK_combined_annotated@meta.data$main.ident, levels = mimicry_order1)

mimicry_order2 <- c("Epithelial Cells", "Lymphoid-like")
Navin_Visvader_NORM_panCK_combined_annotated[["detailed.ident"]] <- factor(x = Navin_Visvader_NORM_panCK_combined_annotated@meta.data$detailed.ident, levels = mimicry_order2)

DimPlot(Navin_Visvader_NORM_panCK_combined_annotated, reduction = "umap", raster = FALSE, group.by = 'detailed.ident', cols = c('Neoplastic Cells' = '#00A9FF', 'Lymphoid-like' = '#CD9600'))


# Lim_Mammary_Stem 
# Score cells based on signature expression
Navin_Visvader_NORM_panCK_combined_annotated <- AddModuleScore(Navin_Visvader_NORM_panCK_combined_annotated, features = Lim_Mammary_Stem, name = "Lim_Mammary_Stem_Up")
# Draw feature plot
FeaturePlot(Navin_Visvader_NORM_panCK_combined_annotated, c("Lim_Mammary_Stem_Up1"))
# Draw Vln plot
VlnPlot(Navin_Visvader_NORM_panCK_combined_annotated, c("Lim_Mammary_Stem_Up1"), group.by = "main.ident", pt.size = 0)
ggsave("/R/R_Navin/Navin_Output/Navin_Visvader_NORM_panCK_combined_annotated.VLN.main.Lim_Mammary_Stem.tiff", plot = last_plot(), device = "tiff",
       scale = 1, width = 6, height = 10,
       dpi = 200, limitsize = TRUE)
VlnPlot(Navin_Visvader_NORM_panCK_combined_annotated, c("Lim_Mammary_Stem_Up1"), group.by = "detailed.ident", pt.size = 0)
ggsave("/R/R_Navin/Navin_Output/Navin_Visvader_NORM_panCK_combined_annotated.VLN.detailed.ident.Lim_Mammary_Stem.tiff", plot = last_plot(), device = "tiff",
       scale = 1, width = 8, height = 10,
       dpi = 200, limitsize = TRUE)



# Lim_Luminal_Progenitor 
# Score cells based on signature expression
Navin_Visvader_NORM_panCK_combined_annotated <- AddModuleScore(Navin_Visvader_NORM_panCK_combined_annotated, features = Lim_Luminal_Progenitor, name = "Lim_Luminal_Progenitor_Up")
# Draw feature plot
FeaturePlot(Navin_Visvader_NORM_panCK_combined_annotated, c("Lim_Luminal_Progenitor_Up1"))
# Draw Vln plot
VlnPlot(Navin_Visvader_NORM_panCK_combined_annotated, c("Lim_Luminal_Progenitor_Up1"), group.by = "main.ident", pt.size = 0)
ggsave("/R/R_Navin/Navin_Output/Navin_Visvader_NORM_panCK_combined_annotated.VLN.main.Lim_Luminal_Progenitor.tiff", plot = last_plot(), device = "tiff",
       scale = 1, width = 6, height = 10,
       dpi = 200, limitsize = TRUE)
VlnPlot(Navin_Visvader_NORM_panCK_combined_annotated, c("Lim_Luminal_Progenitor_Up1"), group.by = "detailed.ident", pt.size = 0)
ggsave("/R/R_Navin/Navin_Output/Navin_Visvader_NORM_panCK_combined_annotated.VLN.detailed.ident.Lim_Luminal_Progenitor.tiff", plot = last_plot(), device = "tiff",
       scale = 1, width = 8, height = 10,
       dpi = 200, limitsize = TRUE)



# Lim_Mature 
# Score cells based on signature expression
Navin_Visvader_NORM_panCK_combined_annotated <- AddModuleScore(Navin_Visvader_NORM_panCK_combined_annotated, features = Lim_Mature, name = "Lim_Mature_Up")
# Draw feature plot
FeaturePlot(Navin_Visvader_NORM_panCK_combined_annotated, c("Lim_Mature_Up1"))
# Draw Vln plot
VlnPlot(Navin_Visvader_NORM_panCK_combined_annotated, c("Lim_Mature_Up1"), group.by = "main.ident", pt.size = 0)
ggsave("/R/R_Navin/Navin_Output/Navin_Visvader_NORM_panCK_combined_annotated.VLN.main.Lim_Mature.tiff", plot = last_plot(), device = "tiff",
       scale = 1, width = 6, height = 10,
       dpi = 200, limitsize = TRUE)
VlnPlot(Navin_Visvader_NORM_panCK_combined_annotated, c("Lim_Mature_Up1"), group.by = "detailed.ident", pt.size = 0)
ggsave("/R/R_Navin/Navin_Output/Navin_Visvader_NORM_panCK_combined_annotated.VLN.detailed.ident.Lim_Mature.tiff", plot = last_plot(), device = "tiff",
       scale = 1, width = 8, height = 10,
       dpi = 200, limitsize = TRUE)


#################################################
###      Functional Signature Analysis        ###
#################################################

# BIOCARTA_NFKB_PATHWAY 
# Score cells based on signature expression
Navin_Visvader_NORM_panCK_combined_annotated <- AddModuleScore(Navin_Visvader_NORM_panCK_combined_annotated, features = BIOCARTA_NFKB_PATHWAY, name = "BIOCARTA_NFKB_PATHWAY")
# Draw feature plot
FeaturePlot(Navin_Visvader_NORM_panCK_combined_annotated, c("BIOCARTA_NFKB_PATHWAY1"))
# Draw Vln plot
VlnPlot(Navin_Visvader_NORM_panCK_combined_annotated, c("BIOCARTA_NFKB_PATHWAY1"), group.by = "main.ident", pt.size = 0)
ggsave("/R/R_Navin/Navin_Output/Navin_Visvader_NORM_panCK_combined_annotated.VLN.main.ident.BIOCARTA_NFKB_PATHWAY.tiff", plot = last_plot(), device = "tiff",
       scale = 1, width = 6, height = 10,
       dpi = 200, limitsize = TRUE)
VlnPlot(Navin_Visvader_NORM_panCK_combined_annotated, c("BIOCARTA_NFKB_PATHWAY1"), group.by = "detailed.ident", pt.size = 0)
ggsave("/R/R_Navin/Navin_Output/Navin_Visvader_NORM_panCK_combined_annotated.VLN.detailed.ident.BIOCARTA_NFKB_PATHWAY.tiff", plot = last_plot(), device = "tiff",
       scale = 1, width = 8, height = 10,
       dpi = 200, limitsize = TRUE)


# POSITIVE_REGULATION_OF_LEUKOCYTE_PROLIFERATION
# Score cells based on signature expression
Navin_Visvader_NORM_panCK_combined_annotated <- AddModuleScore(Navin_Visvader_NORM_panCK_combined_annotated, features = POSITIVE_REGULATION_OF_LEUKOCYTE_PROLIFERATION, name = "POSITIVE_REGULATION_OF_LEUKOCYTE_PROLIFERATION")
# Draw feature plot
FeaturePlot(Navin_Visvader_NORM_panCK_combined_annotated, c("POSITIVE_REGULATION_OF_LEUKOCYTE_PROLIFERATION1"))
# Draw Vln plot
VlnPlot(Navin_Visvader_NORM_panCK_combined_annotated, c("POSITIVE_REGULATION_OF_LEUKOCYTE_PROLIFERATION1"), group.by = "main.ident", pt.size = 0)
ggsave("/R/R_Navin/Navin_Output/Navin_Visvader_NORM_panCK_combined_annotated.VLN.main.POSITIVE_REGULATION_OF_LEUKOCYTE_PROLIFERATION.tiff", plot = last_plot(), device = "tiff",
       scale = 1, width = 6, height = 10,
       dpi = 200, limitsize = TRUE)
VlnPlot(Navin_Visvader_NORM_panCK_combined_annotated, c("POSITIVE_REGULATION_OF_LEUKOCYTE_PROLIFERATION1"), group.by = "detailed.ident", pt.size = 0)
ggsave("/R/R_Navin/Navin_Output/Navin_Visvader_NORM_panCK_combined_annotated.VLN.detailed.ident.POSITIVE_REGULATION_OF_LEUKOCYTE_PROLIFERATION.tiff", plot = last_plot(), device = "tiff",
       scale = 1, width = 8, height = 10,
       dpi = 200, limitsize = TRUE)


# POSITIVE_REGULATION_OF_MAPK_CASCADE
# Score cells based on signature expression
Navin_Visvader_NORM_panCK_combined_annotated <- AddModuleScore(Navin_Visvader_NORM_panCK_combined_annotated, features = POSITIVE_REGULATION_OF_MAPK_CASCADE, name = "POSITIVE_REGULATION_OF_MAPK_CASCADE")
# Draw feature plot
FeaturePlot(Navin_Visvader_NORM_panCK_combined_annotated, c("POSITIVE_REGULATION_OF_MAPK_CASCADE1"))
# Draw Vln plot
VlnPlot(Navin_Visvader_NORM_panCK_combined_annotated, c("POSITIVE_REGULATION_OF_MAPK_CASCADE1"), group.by = "main.ident", pt.size = 0)
ggsave("/R/R_Navin/Navin_Output/Navin_Visvader_NORM_panCK_combined_annotated.VLN.main.POSITIVE_REGULATION_OF_MAPK_CASCADE.tiff", plot = last_plot(), device = "tiff",
       scale = 1, width = 6, height = 10,
       dpi = 200, limitsize = TRUE)
VlnPlot(Navin_Visvader_NORM_panCK_combined_annotated, c("POSITIVE_REGULATION_OF_MAPK_CASCADE1"), group.by = "detailed.ident", pt.size = 0)
ggsave("/R/R_Navin/Navin_Output/Navin_Visvader_NORM_panCK_combined_annotated.VLN.detailed.ident.POSITIVE_REGULATION_OF_MAPK_CASCADE.tiff", plot = last_plot(), device = "tiff",
       scale = 1, width = 8, height = 10,
       dpi = 200, limitsize = TRUE)


# POSITIVE_REGULATION_OF_CYTOKINE_PRODUCTION
# Score cells based on signature expression
Navin_Visvader_NORM_panCK_combined_annotated <- AddModuleScore(Navin_Visvader_NORM_panCK_combined_annotated, features = POSITIVE_REGULATION_OF_CYTOKINE_PRODUCTION, name = "POSITIVE_REGULATION_OF_CYTOKINE_PRODUCTION")
# Draw feature plot
FeaturePlot(Navin_Visvader_NORM_panCK_combined_annotated, c("POSITIVE_REGULATION_OF_CYTOKINE_PRODUCTION1"))
# Draw Vln plot
VlnPlot(Navin_Visvader_NORM_panCK_combined_annotated, c("POSITIVE_REGULATION_OF_CYTOKINE_PRODUCTION1"), group.by = "main.ident", pt.size = 0)
ggsave("/R/R_Navin/Navin_Output/Navin_Visvader_NORM_panCK_combined_annotated.VLN.main.POSITIVE_REGULATION_OF_CYTOKINE_PRODUCTION.tiff", plot = last_plot(), device = "tiff",
       scale = 1, width = 6, height = 10,
       dpi = 200, limitsize = TRUE)
VlnPlot(Navin_Visvader_NORM_panCK_combined_annotated, c("POSITIVE_REGULATION_OF_CYTOKINE_PRODUCTION1"), group.by = "detailed.ident", pt.size = 0)
ggsave("/R/R_Navin/Navin_Output/Navin_Visvader_NORM_panCK_combined_annotated.VLN.detailed.ident.POSITIVE_REGULATION_OF_CYTOKINE_PRODUCTION.tiff", plot = last_plot(), device = "tiff",
       scale = 1, width = 8, height = 10,
       dpi = 200, limitsize = TRUE)


# POSITIVE_REGULATION_OF_LYMPHOCYTE_MIGRATION
# Score cells based on signature expression
Navin_Visvader_NORM_panCK_combined_annotated <- AddModuleScore(Navin_Visvader_NORM_panCK_combined_annotated, features = POSITIVE_REGULATION_OF_LYMPHOCYTE_MIGRATION, name = "POSITIVE_REGULATION_OF_LYMPHOCYTE_MIGRATION")
# Draw feature plot
FeaturePlot(Navin_Visvader_NORM_panCK_combined_annotated, c("POSITIVE_REGULATION_OF_LYMPHOCYTE_MIGRATION1"))
# Draw Vln plot
VlnPlot(Navin_Visvader_NORM_panCK_combined_annotated, c("POSITIVE_REGULATION_OF_LYMPHOCYTE_MIGRATION1"), group.by = "main.ident", pt.size = 0)
ggsave("/R/R_Navin/Navin_Output/Navin_Visvader_NORM_panCK_combined_annotated.VLN.main.POSITIVE_REGULATION_OF_LYMPHOCYTE_MIGRATION.tiff", plot = last_plot(), device = "tiff",
       scale = 1, width = 6, height = 10,
       dpi = 200, limitsize = TRUE)
VlnPlot(Navin_Visvader_NORM_panCK_combined_annotated, c("POSITIVE_REGULATION_OF_LYMPHOCYTE_MIGRATION1"), group.by = "detailed.ident", pt.size = 0)
ggsave("/R/R_Navin/Navin_Output/Navin_Visvader_NORM_panCK_combined_annotated.VLN.detailed.identPOSITIVE_REGULATION_OF_LYMPHOCYTE_MIGRATION.tiff", plot = last_plot(), device = "tiff",
       scale = 1, width = 8, height = 10,
       dpi = 200, limitsize = TRUE)


# Extract  Meta Data
Navin_Visvader_NORM_panCK_combined_annotated.signatures.meta.data <- as.data.frame(as.matrix(Navin_Visvader_NORM_panCK_combined_annotated@meta.data))

# Save Meta Data
write.csv(Navin_Visvader_NORM_panCK_combined_annotated.signatures.meta.data, file = "/R/R_Navin/Navin_Output/Navin_Visvader_NORM_panCK_combined_annotated.signatures.meta.data.csv")

# Remove Meta Data
rm(Navin_Visvader_NORM_panCK_combined_annotated.signatures.meta.data)

# Save Annotated Seurat Object with Signature Scoring in Meta Data
saveRDS(Navin_Visvader_NORM_panCK_combined_annotated, file = "/R/R_Navin/Navin_RDS/RDS_Annotated/Navin_Visvader_NORM_panCK_combined_annotated.rds")

# Clear memory
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




#########################################################
# Step 18: Generate Expression Matrices for Each Sample # 
#########################################################
# Load Libraries
library(Seurat)
library(dplyr)
library(patchwork)

######################################################
################## Navin Dataset #####################
######################################################

############################ 
# Navin_hbca_c14_panCK   #
############################ 
# Load Object
Navin_hbca_c14_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c14_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c14_panCK_All <- Navin_hbca_c14_panCK[["RNA"]]$counts
Navin_hbca_c14_panCK_All <- as.matrix(Navin_hbca_c14_panCK_All, 'sparseMatrix')
Navin_hbca_c14_panCK_All <-t(Navin_hbca_c14_panCK_All)
write.csv(Navin_hbca_c14_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c14_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c14_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c14_panCK_IM <- Navin_hbca_c14_panCK[["RNA"]]$counts
Navin_hbca_c14_panCK_IM <- as.matrix(Navin_hbca_c14_panCK_IM, 'sparseMatrix')
Navin_hbca_c14_panCK_IM <- subset(Navin_hbca_c14_panCK_IM, rownames(Navin_hbca_c14_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c14_panCK_IM <-t(Navin_hbca_c14_panCK_IM)
write.csv(Navin_hbca_c14_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c14_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c14_panCK_IM)

# Remove Object
rm(Navin_hbca_c14_panCK)
gc()



############################ 
# Navin_hbca_c15_panCK   #
############################ 
# Load Object
Navin_hbca_c15_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c15_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c15_panCK_All <- Navin_hbca_c15_panCK[["RNA"]]$counts
Navin_hbca_c15_panCK_All <- as.matrix(Navin_hbca_c15_panCK_All, 'sparseMatrix')
Navin_hbca_c15_panCK_All <-t(Navin_hbca_c15_panCK_All)
write.csv(Navin_hbca_c15_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c15_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c15_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c15_panCK_IM <- Navin_hbca_c15_panCK[["RNA"]]$counts
Navin_hbca_c15_panCK_IM <- as.matrix(Navin_hbca_c15_panCK_IM, 'sparseMatrix')
Navin_hbca_c15_panCK_IM <- subset(Navin_hbca_c15_panCK_IM, rownames(Navin_hbca_c15_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c15_panCK_IM <-t(Navin_hbca_c15_panCK_IM)
write.csv(Navin_hbca_c15_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c15_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c15_panCK_IM)

# Remove Object
rm(Navin_hbca_c15_panCK)
gc()



############################ 
# Navin_hbca_c19_panCK   #
############################ 
# Load Object
Navin_hbca_c19_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c19_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c19_panCK_All <- Navin_hbca_c19_panCK[["RNA"]]$counts
Navin_hbca_c19_panCK_All <- as.matrix(Navin_hbca_c19_panCK_All, 'sparseMatrix')
Navin_hbca_c19_panCK_All <-t(Navin_hbca_c19_panCK_All)
write.csv(Navin_hbca_c19_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c19_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c19_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c19_panCK_IM <- Navin_hbca_c19_panCK[["RNA"]]$counts
Navin_hbca_c19_panCK_IM <- as.matrix(Navin_hbca_c19_panCK_IM, 'sparseMatrix')
Navin_hbca_c19_panCK_IM <- subset(Navin_hbca_c19_panCK_IM, rownames(Navin_hbca_c19_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c19_panCK_IM <-t(Navin_hbca_c19_panCK_IM)
write.csv(Navin_hbca_c19_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c19_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c19_panCK_IM)

# Remove Object
rm(Navin_hbca_c19_panCK)
gc()



############################ 
# Navin_hbca_c20_panCK   #
############################ 
# Load Object
Navin_hbca_c20_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c20_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c20_panCK_All <- Navin_hbca_c20_panCK[["RNA"]]$counts
Navin_hbca_c20_panCK_All <- as.matrix(Navin_hbca_c20_panCK_All, 'sparseMatrix')
Navin_hbca_c20_panCK_All <-t(Navin_hbca_c20_panCK_All)
write.csv(Navin_hbca_c20_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c20_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c20_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c20_panCK_IM <- Navin_hbca_c20_panCK[["RNA"]]$counts
Navin_hbca_c20_panCK_IM <- as.matrix(Navin_hbca_c20_panCK_IM, 'sparseMatrix')
Navin_hbca_c20_panCK_IM <- subset(Navin_hbca_c20_panCK_IM, rownames(Navin_hbca_c20_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c20_panCK_IM <-t(Navin_hbca_c20_panCK_IM)
write.csv(Navin_hbca_c20_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c20_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c20_panCK_IM)

# Remove Object
rm(Navin_hbca_c20_panCK)
gc()



############################ 
# Navin_hbca_c22_panCK   #
############################ 
# Load Object
Navin_hbca_c22_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c22_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c22_panCK_All <- Navin_hbca_c22_panCK[["RNA"]]$counts
Navin_hbca_c22_panCK_All <- as.matrix(Navin_hbca_c22_panCK_All, 'sparseMatrix')
Navin_hbca_c22_panCK_All <-t(Navin_hbca_c22_panCK_All)
write.csv(Navin_hbca_c22_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c22_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c22_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c22_panCK_IM <- Navin_hbca_c22_panCK[["RNA"]]$counts
Navin_hbca_c22_panCK_IM <- as.matrix(Navin_hbca_c22_panCK_IM, 'sparseMatrix')
Navin_hbca_c22_panCK_IM <- subset(Navin_hbca_c22_panCK_IM, rownames(Navin_hbca_c22_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c22_panCK_IM <-t(Navin_hbca_c22_panCK_IM)
write.csv(Navin_hbca_c22_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c22_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c22_panCK_IM)

# Remove Object
rm(Navin_hbca_c22_panCK)
gc()



############################ 
# Navin_hbca_c23_panCK   #
############################ 
# Load Object
Navin_hbca_c23_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c23_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c23_panCK_All <- Navin_hbca_c23_panCK[["RNA"]]$counts
Navin_hbca_c23_panCK_All <- as.matrix(Navin_hbca_c23_panCK_All, 'sparseMatrix')
Navin_hbca_c23_panCK_All <-t(Navin_hbca_c23_panCK_All)
write.csv(Navin_hbca_c23_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c23_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c23_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c23_panCK_IM <- Navin_hbca_c23_panCK[["RNA"]]$counts
Navin_hbca_c23_panCK_IM <- as.matrix(Navin_hbca_c23_panCK_IM, 'sparseMatrix')
Navin_hbca_c23_panCK_IM <- subset(Navin_hbca_c23_panCK_IM, rownames(Navin_hbca_c23_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c23_panCK_IM <-t(Navin_hbca_c23_panCK_IM)
write.csv(Navin_hbca_c23_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c23_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c23_panCK_IM)

# Remove Object
rm(Navin_hbca_c23_panCK)
gc()



############################ 
# Navin_hbca_c24_panCK   #
############################ 
# Load Object
Navin_hbca_c24_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c24_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c24_panCK_All <- Navin_hbca_c24_panCK[["RNA"]]$counts
Navin_hbca_c24_panCK_All <- as.matrix(Navin_hbca_c24_panCK_All, 'sparseMatrix')
Navin_hbca_c24_panCK_All <-t(Navin_hbca_c24_panCK_All)
write.csv(Navin_hbca_c24_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c24_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c24_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c24_panCK_IM <- Navin_hbca_c24_panCK[["RNA"]]$counts
Navin_hbca_c24_panCK_IM <- as.matrix(Navin_hbca_c24_panCK_IM, 'sparseMatrix')
Navin_hbca_c24_panCK_IM <- subset(Navin_hbca_c24_panCK_IM, rownames(Navin_hbca_c24_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c24_panCK_IM <-t(Navin_hbca_c24_panCK_IM)
write.csv(Navin_hbca_c24_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c24_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c24_panCK_IM)

# Remove Object
rm(Navin_hbca_c24_panCK)
gc()



############################ 
# Navin_hbca_c25_panCK   #
############################ 
# Load Object
Navin_hbca_c25_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c25_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c25_panCK_All <- Navin_hbca_c25_panCK[["RNA"]]$counts
Navin_hbca_c25_panCK_All <- as.matrix(Navin_hbca_c25_panCK_All, 'sparseMatrix')
Navin_hbca_c25_panCK_All <-t(Navin_hbca_c25_panCK_All)
write.csv(Navin_hbca_c25_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c25_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c25_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c25_panCK_IM <- Navin_hbca_c25_panCK[["RNA"]]$counts
Navin_hbca_c25_panCK_IM <- as.matrix(Navin_hbca_c25_panCK_IM, 'sparseMatrix')
Navin_hbca_c25_panCK_IM <- subset(Navin_hbca_c25_panCK_IM, rownames(Navin_hbca_c25_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c25_panCK_IM <-t(Navin_hbca_c25_panCK_IM)
write.csv(Navin_hbca_c25_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c25_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c25_panCK_IM)

# Remove Object
rm(Navin_hbca_c25_panCK)
gc()



############################ 
# Navin_hbca_c26_panCK   #
############################ 
# Load Object
Navin_hbca_c26_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c26_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c26_panCK_All <- Navin_hbca_c26_panCK[["RNA"]]$counts
Navin_hbca_c26_panCK_All <- as.matrix(Navin_hbca_c26_panCK_All, 'sparseMatrix')
Navin_hbca_c26_panCK_All <-t(Navin_hbca_c26_panCK_All)
write.csv(Navin_hbca_c26_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c26_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c26_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c26_panCK_IM <- Navin_hbca_c26_panCK[["RNA"]]$counts
Navin_hbca_c26_panCK_IM <- as.matrix(Navin_hbca_c26_panCK_IM, 'sparseMatrix')
Navin_hbca_c26_panCK_IM <- subset(Navin_hbca_c26_panCK_IM, rownames(Navin_hbca_c26_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c26_panCK_IM <-t(Navin_hbca_c26_panCK_IM)
write.csv(Navin_hbca_c26_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c26_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c26_panCK_IM)

# Remove Object
rm(Navin_hbca_c26_panCK)
gc()



############################ 
# Navin_hbca_c31_panCK   #
############################ 
# Load Object
Navin_hbca_c31_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c31_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c31_panCK_All <- Navin_hbca_c31_panCK[["RNA"]]$counts
Navin_hbca_c31_panCK_All <- as.matrix(Navin_hbca_c31_panCK_All, 'sparseMatrix')
Navin_hbca_c31_panCK_All <-t(Navin_hbca_c31_panCK_All)
write.csv(Navin_hbca_c31_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c31_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c31_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c31_panCK_IM <- Navin_hbca_c31_panCK[["RNA"]]$counts
Navin_hbca_c31_panCK_IM <- as.matrix(Navin_hbca_c31_panCK_IM, 'sparseMatrix')
Navin_hbca_c31_panCK_IM <- subset(Navin_hbca_c31_panCK_IM, rownames(Navin_hbca_c31_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c31_panCK_IM <-t(Navin_hbca_c31_panCK_IM)
write.csv(Navin_hbca_c31_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c31_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c31_panCK_IM)

# Remove Object
rm(Navin_hbca_c31_panCK)
gc()





############################ 
# Navin_hbca_c32_panCK   #
############################ 
# Load Object
Navin_hbca_c32_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c32_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c32_panCK_All <- Navin_hbca_c32_panCK[["RNA"]]$counts
Navin_hbca_c32_panCK_All <- as.matrix(Navin_hbca_c32_panCK_All, 'sparseMatrix')
Navin_hbca_c32_panCK_All <-t(Navin_hbca_c32_panCK_All)
write.csv(Navin_hbca_c32_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c32_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c32_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c32_panCK_IM <- Navin_hbca_c32_panCK[["RNA"]]$counts
Navin_hbca_c32_panCK_IM <- as.matrix(Navin_hbca_c32_panCK_IM, 'sparseMatrix')
Navin_hbca_c32_panCK_IM <- subset(Navin_hbca_c32_panCK_IM, rownames(Navin_hbca_c32_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c32_panCK_IM <-t(Navin_hbca_c32_panCK_IM)
write.csv(Navin_hbca_c32_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c32_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c32_panCK_IM)

# Remove Object
rm(Navin_hbca_c32_panCK)
gc()



############################ 
# Navin_hbca_c50_panCK   #
############################ 
# Load Object
Navin_hbca_c50_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c50_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c50_panCK_All <- Navin_hbca_c50_panCK[["RNA"]]$counts
Navin_hbca_c50_panCK_All <- as.matrix(Navin_hbca_c50_panCK_All, 'sparseMatrix')
Navin_hbca_c50_panCK_All <-t(Navin_hbca_c50_panCK_All)
write.csv(Navin_hbca_c50_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c50_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c50_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c50_panCK_IM <- Navin_hbca_c50_panCK[["RNA"]]$counts
Navin_hbca_c50_panCK_IM <- as.matrix(Navin_hbca_c50_panCK_IM, 'sparseMatrix')
Navin_hbca_c50_panCK_IM <- subset(Navin_hbca_c50_panCK_IM, rownames(Navin_hbca_c50_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c50_panCK_IM <-t(Navin_hbca_c50_panCK_IM)
write.csv(Navin_hbca_c50_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c50_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c50_panCK_IM)

# Remove Object
rm(Navin_hbca_c50_panCK)
gc()



############################ 
# Navin_hbca_c51_panCK   #
############################ 
# Load Object
Navin_hbca_c51_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c51_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c51_panCK_All <- Navin_hbca_c51_panCK[["RNA"]]$counts
Navin_hbca_c51_panCK_All <- as.matrix(Navin_hbca_c51_panCK_All, 'sparseMatrix')
Navin_hbca_c51_panCK_All <-t(Navin_hbca_c51_panCK_All)
write.csv(Navin_hbca_c51_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c51_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c51_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c51_panCK_IM <- Navin_hbca_c51_panCK[["RNA"]]$counts
Navin_hbca_c51_panCK_IM <- as.matrix(Navin_hbca_c51_panCK_IM, 'sparseMatrix')
Navin_hbca_c51_panCK_IM <- subset(Navin_hbca_c51_panCK_IM, rownames(Navin_hbca_c51_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c51_panCK_IM <-t(Navin_hbca_c51_panCK_IM)
write.csv(Navin_hbca_c51_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c51_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c51_panCK_IM)

# Remove Object
rm(Navin_hbca_c51_panCK)
gc()



############################ 
# Navin_hbca_c52_panCK   #
############################ 
# Load Object
Navin_hbca_c52_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c52_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c52_panCK_All <- Navin_hbca_c52_panCK[["RNA"]]$counts
Navin_hbca_c52_panCK_All <- as.matrix(Navin_hbca_c52_panCK_All, 'sparseMatrix')
Navin_hbca_c52_panCK_All <-t(Navin_hbca_c52_panCK_All)
write.csv(Navin_hbca_c52_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c52_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c52_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c52_panCK_IM <- Navin_hbca_c52_panCK[["RNA"]]$counts
Navin_hbca_c52_panCK_IM <- as.matrix(Navin_hbca_c52_panCK_IM, 'sparseMatrix')
Navin_hbca_c52_panCK_IM <- subset(Navin_hbca_c52_panCK_IM, rownames(Navin_hbca_c52_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c52_panCK_IM <-t(Navin_hbca_c52_panCK_IM)
write.csv(Navin_hbca_c52_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c52_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c52_panCK_IM)

# Remove Object
rm(Navin_hbca_c52_panCK)
gc()



############################ 
# Navin_hbca_c53_panCK   #
############################ 
# Load Object
Navin_hbca_c53_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c53_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c53_panCK_All <- Navin_hbca_c53_panCK[["RNA"]]$counts
Navin_hbca_c53_panCK_All <- as.matrix(Navin_hbca_c53_panCK_All, 'sparseMatrix')
Navin_hbca_c53_panCK_All <-t(Navin_hbca_c53_panCK_All)
write.csv(Navin_hbca_c53_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c53_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c53_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c53_panCK_IM <- Navin_hbca_c53_panCK[["RNA"]]$counts
Navin_hbca_c53_panCK_IM <- as.matrix(Navin_hbca_c53_panCK_IM, 'sparseMatrix')
Navin_hbca_c53_panCK_IM <- subset(Navin_hbca_c53_panCK_IM, rownames(Navin_hbca_c53_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c53_panCK_IM <-t(Navin_hbca_c53_panCK_IM)
write.csv(Navin_hbca_c53_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c53_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c53_panCK_IM)

# Remove Object
rm(Navin_hbca_c53_panCK)
gc()



############################ 
# Navin_hbca_c54_panCK   #
############################ 
# Load Object
Navin_hbca_c54_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c54_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c54_panCK_All <- Navin_hbca_c54_panCK[["RNA"]]$counts
Navin_hbca_c54_panCK_All <- as.matrix(Navin_hbca_c54_panCK_All, 'sparseMatrix')
Navin_hbca_c54_panCK_All <-t(Navin_hbca_c54_panCK_All)
write.csv(Navin_hbca_c54_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c54_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c54_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c54_panCK_IM <- Navin_hbca_c54_panCK[["RNA"]]$counts
Navin_hbca_c54_panCK_IM <- as.matrix(Navin_hbca_c54_panCK_IM, 'sparseMatrix')
Navin_hbca_c54_panCK_IM <- subset(Navin_hbca_c54_panCK_IM, rownames(Navin_hbca_c54_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c54_panCK_IM <-t(Navin_hbca_c54_panCK_IM)
write.csv(Navin_hbca_c54_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c54_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c54_panCK_IM)

# Remove Object
rm(Navin_hbca_c54_panCK)
gc()



############################ 
# Navin_hbca_c55_panCK   #
############################ 
# Load Object
Navin_hbca_c55_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c55_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c55_panCK_All <- Navin_hbca_c55_panCK[["RNA"]]$counts
Navin_hbca_c55_panCK_All <- as.matrix(Navin_hbca_c55_panCK_All, 'sparseMatrix')
Navin_hbca_c55_panCK_All <-t(Navin_hbca_c55_panCK_All)
write.csv(Navin_hbca_c55_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c55_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c55_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c55_panCK_IM <- Navin_hbca_c55_panCK[["RNA"]]$counts
Navin_hbca_c55_panCK_IM <- as.matrix(Navin_hbca_c55_panCK_IM, 'sparseMatrix')
Navin_hbca_c55_panCK_IM <- subset(Navin_hbca_c55_panCK_IM, rownames(Navin_hbca_c55_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c55_panCK_IM <-t(Navin_hbca_c55_panCK_IM)
write.csv(Navin_hbca_c55_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c55_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c55_panCK_IM)

# Remove Object
rm(Navin_hbca_c55_panCK)
gc()



############################ 
# Navin_hbca_c56_panCK   #
############################ 
# Load Object
Navin_hbca_c56_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c56_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c56_panCK_All <- Navin_hbca_c56_panCK[["RNA"]]$counts
Navin_hbca_c56_panCK_All <- as.matrix(Navin_hbca_c56_panCK_All, 'sparseMatrix')
Navin_hbca_c56_panCK_All <-t(Navin_hbca_c56_panCK_All)
write.csv(Navin_hbca_c56_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c56_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c56_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c56_panCK_IM <- Navin_hbca_c56_panCK[["RNA"]]$counts
Navin_hbca_c56_panCK_IM <- as.matrix(Navin_hbca_c56_panCK_IM, 'sparseMatrix')
Navin_hbca_c56_panCK_IM <- subset(Navin_hbca_c56_panCK_IM, rownames(Navin_hbca_c56_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c56_panCK_IM <-t(Navin_hbca_c56_panCK_IM)
write.csv(Navin_hbca_c56_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c56_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c56_panCK_IM)

# Remove Object
rm(Navin_hbca_c56_panCK)
gc()



############################ 
# Navin_hbca_c57_panCK   #
############################ 
# Load Object
Navin_hbca_c57_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c57_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c57_panCK_All <- Navin_hbca_c57_panCK[["RNA"]]$counts
Navin_hbca_c57_panCK_All <- as.matrix(Navin_hbca_c57_panCK_All, 'sparseMatrix')
Navin_hbca_c57_panCK_All <-t(Navin_hbca_c57_panCK_All)
write.csv(Navin_hbca_c57_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c57_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c57_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c57_panCK_IM <- Navin_hbca_c57_panCK[["RNA"]]$counts
Navin_hbca_c57_panCK_IM <- as.matrix(Navin_hbca_c57_panCK_IM, 'sparseMatrix')
Navin_hbca_c57_panCK_IM <- subset(Navin_hbca_c57_panCK_IM, rownames(Navin_hbca_c57_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c57_panCK_IM <-t(Navin_hbca_c57_panCK_IM)
write.csv(Navin_hbca_c57_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c57_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c57_panCK_IM)

# Remove Object
rm(Navin_hbca_c57_panCK)
gc()



############################ 
# Navin_hbca_c58_panCK   #
############################ 
# Load Object
Navin_hbca_c58_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c58_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c58_panCK_All <- Navin_hbca_c58_panCK[["RNA"]]$counts
Navin_hbca_c58_panCK_All <- as.matrix(Navin_hbca_c58_panCK_All, 'sparseMatrix')
Navin_hbca_c58_panCK_All <-t(Navin_hbca_c58_panCK_All)
write.csv(Navin_hbca_c58_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c58_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c58_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c58_panCK_IM <- Navin_hbca_c58_panCK[["RNA"]]$counts
Navin_hbca_c58_panCK_IM <- as.matrix(Navin_hbca_c58_panCK_IM, 'sparseMatrix')
Navin_hbca_c58_panCK_IM <- subset(Navin_hbca_c58_panCK_IM, rownames(Navin_hbca_c58_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c58_panCK_IM <-t(Navin_hbca_c58_panCK_IM)
write.csv(Navin_hbca_c58_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c58_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c58_panCK_IM)

# Remove Object
rm(Navin_hbca_c58_panCK)
gc()





############################ 
# Navin_hbca_c59_panCK   #
############################ 
# Load Object
Navin_hbca_c59_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c59_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c59_panCK_All <- Navin_hbca_c59_panCK[["RNA"]]$counts
Navin_hbca_c59_panCK_All <- as.matrix(Navin_hbca_c59_panCK_All, 'sparseMatrix')
Navin_hbca_c59_panCK_All <-t(Navin_hbca_c59_panCK_All)
write.csv(Navin_hbca_c59_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c59_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c59_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c59_panCK_IM <- Navin_hbca_c59_panCK[["RNA"]]$counts
Navin_hbca_c59_panCK_IM <- as.matrix(Navin_hbca_c59_panCK_IM, 'sparseMatrix')
Navin_hbca_c59_panCK_IM <- subset(Navin_hbca_c59_panCK_IM, rownames(Navin_hbca_c59_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c59_panCK_IM <-t(Navin_hbca_c59_panCK_IM)
write.csv(Navin_hbca_c59_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c59_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c59_panCK_IM)

# Remove Object
rm(Navin_hbca_c59_panCK)
gc()



############################ 
# Navin_hbca_c60_panCK   #
############################ 
# Load Object
Navin_hbca_c60_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c60_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c60_panCK_All <- Navin_hbca_c60_panCK[["RNA"]]$counts
Navin_hbca_c60_panCK_All <- as.matrix(Navin_hbca_c60_panCK_All, 'sparseMatrix')
Navin_hbca_c60_panCK_All <-t(Navin_hbca_c60_panCK_All)
write.csv(Navin_hbca_c60_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c60_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c60_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c60_panCK_IM <- Navin_hbca_c60_panCK[["RNA"]]$counts
Navin_hbca_c60_panCK_IM <- as.matrix(Navin_hbca_c60_panCK_IM, 'sparseMatrix')
Navin_hbca_c60_panCK_IM <- subset(Navin_hbca_c60_panCK_IM, rownames(Navin_hbca_c60_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c60_panCK_IM <-t(Navin_hbca_c60_panCK_IM)
write.csv(Navin_hbca_c60_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c60_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c60_panCK_IM)

# Remove Object
rm(Navin_hbca_c60_panCK)
gc()



############################ 
# Navin_hbca_c61_panCK   #
############################ 
# Load Object
Navin_hbca_c61_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c61_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c61_panCK_All <- Navin_hbca_c61_panCK[["RNA"]]$counts
Navin_hbca_c61_panCK_All <- as.matrix(Navin_hbca_c61_panCK_All, 'sparseMatrix')
Navin_hbca_c61_panCK_All <-t(Navin_hbca_c61_panCK_All)
write.csv(Navin_hbca_c61_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c61_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c61_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c61_panCK_IM <- Navin_hbca_c61_panCK[["RNA"]]$counts
Navin_hbca_c61_panCK_IM <- as.matrix(Navin_hbca_c61_panCK_IM, 'sparseMatrix')
Navin_hbca_c61_panCK_IM <- subset(Navin_hbca_c61_panCK_IM, rownames(Navin_hbca_c61_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c61_panCK_IM <-t(Navin_hbca_c61_panCK_IM)
write.csv(Navin_hbca_c61_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c61_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c61_panCK_IM)

# Remove Object
rm(Navin_hbca_c61_panCK)
gc()



############################ 
# Navin_hbca_c62_panCK   #
############################ 
# Load Object
Navin_hbca_c62_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c62_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c62_panCK_All <- Navin_hbca_c62_panCK[["RNA"]]$counts
Navin_hbca_c62_panCK_All <- as.matrix(Navin_hbca_c62_panCK_All, 'sparseMatrix')
Navin_hbca_c62_panCK_All <-t(Navin_hbca_c62_panCK_All)
write.csv(Navin_hbca_c62_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c62_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c62_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c62_panCK_IM <- Navin_hbca_c62_panCK[["RNA"]]$counts
Navin_hbca_c62_panCK_IM <- as.matrix(Navin_hbca_c62_panCK_IM, 'sparseMatrix')
Navin_hbca_c62_panCK_IM <- subset(Navin_hbca_c62_panCK_IM, rownames(Navin_hbca_c62_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c62_panCK_IM <-t(Navin_hbca_c62_panCK_IM)
write.csv(Navin_hbca_c62_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c62_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c62_panCK_IM)

# Remove Object
rm(Navin_hbca_c62_panCK)
gc()



############################ 
# Navin_hbca_c63_panCK   #
############################ 
# Load Object
Navin_hbca_c63_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c63_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c63_panCK_All <- Navin_hbca_c63_panCK[["RNA"]]$counts
Navin_hbca_c63_panCK_All <- as.matrix(Navin_hbca_c63_panCK_All, 'sparseMatrix')
Navin_hbca_c63_panCK_All <-t(Navin_hbca_c63_panCK_All)
write.csv(Navin_hbca_c63_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c63_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c63_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c63_panCK_IM <- Navin_hbca_c63_panCK[["RNA"]]$counts
Navin_hbca_c63_panCK_IM <- as.matrix(Navin_hbca_c63_panCK_IM, 'sparseMatrix')
Navin_hbca_c63_panCK_IM <- subset(Navin_hbca_c63_panCK_IM, rownames(Navin_hbca_c63_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c63_panCK_IM <-t(Navin_hbca_c63_panCK_IM)
write.csv(Navin_hbca_c63_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c63_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c63_panCK_IM)

# Remove Object
rm(Navin_hbca_c63_panCK)
gc()



############################ 
# Navin_hbca_c64_panCK   #
############################ 
# Load Object
Navin_hbca_c64_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c64_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c64_panCK_All <- Navin_hbca_c64_panCK[["RNA"]]$counts
Navin_hbca_c64_panCK_All <- as.matrix(Navin_hbca_c64_panCK_All, 'sparseMatrix')
Navin_hbca_c64_panCK_All <-t(Navin_hbca_c64_panCK_All)
write.csv(Navin_hbca_c64_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c64_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c64_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c64_panCK_IM <- Navin_hbca_c64_panCK[["RNA"]]$counts
Navin_hbca_c64_panCK_IM <- as.matrix(Navin_hbca_c64_panCK_IM, 'sparseMatrix')
Navin_hbca_c64_panCK_IM <- subset(Navin_hbca_c64_panCK_IM, rownames(Navin_hbca_c64_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c64_panCK_IM <-t(Navin_hbca_c64_panCK_IM)
write.csv(Navin_hbca_c64_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c64_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c64_panCK_IM)

# Remove Object
rm(Navin_hbca_c64_panCK)
gc()



############################ 
# Navin_hbca_c65_panCK   #
############################ 
# Load Object
Navin_hbca_c65_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c65_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c65_panCK_All <- Navin_hbca_c65_panCK[["RNA"]]$counts
Navin_hbca_c65_panCK_All <- as.matrix(Navin_hbca_c65_panCK_All, 'sparseMatrix')
Navin_hbca_c65_panCK_All <-t(Navin_hbca_c65_panCK_All)
write.csv(Navin_hbca_c65_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c65_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c65_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c65_panCK_IM <- Navin_hbca_c65_panCK[["RNA"]]$counts
Navin_hbca_c65_panCK_IM <- as.matrix(Navin_hbca_c65_panCK_IM, 'sparseMatrix')
Navin_hbca_c65_panCK_IM <- subset(Navin_hbca_c65_panCK_IM, rownames(Navin_hbca_c65_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c65_panCK_IM <-t(Navin_hbca_c65_panCK_IM)
write.csv(Navin_hbca_c65_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c65_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c65_panCK_IM)

# Remove Object
rm(Navin_hbca_c65_panCK)
gc()



############################ 
# Navin_hbca_c66_panCK   #
############################ 
# Load Object
Navin_hbca_c66_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c66_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c66_panCK_All <- Navin_hbca_c66_panCK[["RNA"]]$counts
Navin_hbca_c66_panCK_All <- as.matrix(Navin_hbca_c66_panCK_All, 'sparseMatrix')
Navin_hbca_c66_panCK_All <-t(Navin_hbca_c66_panCK_All)
write.csv(Navin_hbca_c66_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c66_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c66_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c66_panCK_IM <- Navin_hbca_c66_panCK[["RNA"]]$counts
Navin_hbca_c66_panCK_IM <- as.matrix(Navin_hbca_c66_panCK_IM, 'sparseMatrix')
Navin_hbca_c66_panCK_IM <- subset(Navin_hbca_c66_panCK_IM, rownames(Navin_hbca_c66_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c66_panCK_IM <-t(Navin_hbca_c66_panCK_IM)
write.csv(Navin_hbca_c66_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c66_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c66_panCK_IM)

# Remove Object
rm(Navin_hbca_c66_panCK)
gc()



############################ 
# Navin_hbca_c67_panCK   #
############################ 
# Load Object
Navin_hbca_c67_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c67_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c67_panCK_All <- Navin_hbca_c67_panCK[["RNA"]]$counts
Navin_hbca_c67_panCK_All <- as.matrix(Navin_hbca_c67_panCK_All, 'sparseMatrix')
Navin_hbca_c67_panCK_All <-t(Navin_hbca_c67_panCK_All)
write.csv(Navin_hbca_c67_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c67_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c67_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c67_panCK_IM <- Navin_hbca_c67_panCK[["RNA"]]$counts
Navin_hbca_c67_panCK_IM <- as.matrix(Navin_hbca_c67_panCK_IM, 'sparseMatrix')
Navin_hbca_c67_panCK_IM <- subset(Navin_hbca_c67_panCK_IM, rownames(Navin_hbca_c67_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c67_panCK_IM <-t(Navin_hbca_c67_panCK_IM)
write.csv(Navin_hbca_c67_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c67_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c67_panCK_IM)

# Remove Object
rm(Navin_hbca_c67_panCK)
gc()



############################ 
# Navin_hbca_c68_panCK   #
############################ 
# Load Object
Navin_hbca_c68_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c68_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c68_panCK_All <- Navin_hbca_c68_panCK[["RNA"]]$counts
Navin_hbca_c68_panCK_All <- as.matrix(Navin_hbca_c68_panCK_All, 'sparseMatrix')
Navin_hbca_c68_panCK_All <-t(Navin_hbca_c68_panCK_All)
write.csv(Navin_hbca_c68_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c68_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c68_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c68_panCK_IM <- Navin_hbca_c68_panCK[["RNA"]]$counts
Navin_hbca_c68_panCK_IM <- as.matrix(Navin_hbca_c68_panCK_IM, 'sparseMatrix')
Navin_hbca_c68_panCK_IM <- subset(Navin_hbca_c68_panCK_IM, rownames(Navin_hbca_c68_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c68_panCK_IM <-t(Navin_hbca_c68_panCK_IM)
write.csv(Navin_hbca_c68_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c68_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c68_panCK_IM)

# Remove Object
rm(Navin_hbca_c68_panCK)
gc()





############################ 
# Navin_hbca_c69_panCK   #
############################ 
# Load Object
Navin_hbca_c69_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c69_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c69_panCK_All <- Navin_hbca_c69_panCK[["RNA"]]$counts
Navin_hbca_c69_panCK_All <- as.matrix(Navin_hbca_c69_panCK_All, 'sparseMatrix')
Navin_hbca_c69_panCK_All <-t(Navin_hbca_c69_panCK_All)
write.csv(Navin_hbca_c69_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c69_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c69_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c69_panCK_IM <- Navin_hbca_c69_panCK[["RNA"]]$counts
Navin_hbca_c69_panCK_IM <- as.matrix(Navin_hbca_c69_panCK_IM, 'sparseMatrix')
Navin_hbca_c69_panCK_IM <- subset(Navin_hbca_c69_panCK_IM, rownames(Navin_hbca_c69_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c69_panCK_IM <-t(Navin_hbca_c69_panCK_IM)
write.csv(Navin_hbca_c69_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c69_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c69_panCK_IM)

# Remove Object
rm(Navin_hbca_c69_panCK)
gc()



############################ 
# Navin_hbca_c70_panCK   #
############################ 
# Load Object
Navin_hbca_c70_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c70_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c70_panCK_All <- Navin_hbca_c70_panCK[["RNA"]]$counts
Navin_hbca_c70_panCK_All <- as.matrix(Navin_hbca_c70_panCK_All, 'sparseMatrix')
Navin_hbca_c70_panCK_All <-t(Navin_hbca_c70_panCK_All)
write.csv(Navin_hbca_c70_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c70_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c70_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c70_panCK_IM <- Navin_hbca_c70_panCK[["RNA"]]$counts
Navin_hbca_c70_panCK_IM <- as.matrix(Navin_hbca_c70_panCK_IM, 'sparseMatrix')
Navin_hbca_c70_panCK_IM <- subset(Navin_hbca_c70_panCK_IM, rownames(Navin_hbca_c70_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c70_panCK_IM <-t(Navin_hbca_c70_panCK_IM)
write.csv(Navin_hbca_c70_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c70_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c70_panCK_IM)

# Remove Object
rm(Navin_hbca_c70_panCK)
gc()



############################ 
# Navin_hbca_c71_panCK   #
############################ 
# Load Object
Navin_hbca_c71_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c71_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c71_panCK_All <- Navin_hbca_c71_panCK[["RNA"]]$counts
Navin_hbca_c71_panCK_All <- as.matrix(Navin_hbca_c71_panCK_All, 'sparseMatrix')
Navin_hbca_c71_panCK_All <-t(Navin_hbca_c71_panCK_All)
write.csv(Navin_hbca_c71_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c71_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c71_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c71_panCK_IM <- Navin_hbca_c71_panCK[["RNA"]]$counts
Navin_hbca_c71_panCK_IM <- as.matrix(Navin_hbca_c71_panCK_IM, 'sparseMatrix')
Navin_hbca_c71_panCK_IM <- subset(Navin_hbca_c71_panCK_IM, rownames(Navin_hbca_c71_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c71_panCK_IM <-t(Navin_hbca_c71_panCK_IM)
write.csv(Navin_hbca_c71_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c71_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c71_panCK_IM)

# Remove Object
rm(Navin_hbca_c71_panCK)
gc()



############################ 
# Navin_hbca_c72_panCK   #
############################ 
# Load Object
Navin_hbca_c72_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c72_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c72_panCK_All <- Navin_hbca_c72_panCK[["RNA"]]$counts
Navin_hbca_c72_panCK_All <- as.matrix(Navin_hbca_c72_panCK_All, 'sparseMatrix')
Navin_hbca_c72_panCK_All <-t(Navin_hbca_c72_panCK_All)
write.csv(Navin_hbca_c72_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c72_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c72_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c72_panCK_IM <- Navin_hbca_c72_panCK[["RNA"]]$counts
Navin_hbca_c72_panCK_IM <- as.matrix(Navin_hbca_c72_panCK_IM, 'sparseMatrix')
Navin_hbca_c72_panCK_IM <- subset(Navin_hbca_c72_panCK_IM, rownames(Navin_hbca_c72_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c72_panCK_IM <-t(Navin_hbca_c72_panCK_IM)
write.csv(Navin_hbca_c72_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c72_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c72_panCK_IM)

# Remove Object
rm(Navin_hbca_c72_panCK)
gc()



############################ 
# Navin_hbca_c73_panCK   #
############################ 
# Load Object
Navin_hbca_c73_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c73_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c73_panCK_All <- Navin_hbca_c73_panCK[["RNA"]]$counts
Navin_hbca_c73_panCK_All <- as.matrix(Navin_hbca_c73_panCK_All, 'sparseMatrix')
Navin_hbca_c73_panCK_All <-t(Navin_hbca_c73_panCK_All)
write.csv(Navin_hbca_c73_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c73_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c73_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c73_panCK_IM <- Navin_hbca_c73_panCK[["RNA"]]$counts
Navin_hbca_c73_panCK_IM <- as.matrix(Navin_hbca_c73_panCK_IM, 'sparseMatrix')
Navin_hbca_c73_panCK_IM <- subset(Navin_hbca_c73_panCK_IM, rownames(Navin_hbca_c73_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c73_panCK_IM <-t(Navin_hbca_c73_panCK_IM)
write.csv(Navin_hbca_c73_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c73_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c73_panCK_IM)

# Remove Object
rm(Navin_hbca_c73_panCK)
gc()



############################ 
# Navin_hbca_c74_panCK   #
############################ 
# Load Object
Navin_hbca_c74_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c74_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c74_panCK_All <- Navin_hbca_c74_panCK[["RNA"]]$counts
Navin_hbca_c74_panCK_All <- as.matrix(Navin_hbca_c74_panCK_All, 'sparseMatrix')
Navin_hbca_c74_panCK_All <-t(Navin_hbca_c74_panCK_All)
write.csv(Navin_hbca_c74_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c74_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c74_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c74_panCK_IM <- Navin_hbca_c74_panCK[["RNA"]]$counts
Navin_hbca_c74_panCK_IM <- as.matrix(Navin_hbca_c74_panCK_IM, 'sparseMatrix')
Navin_hbca_c74_panCK_IM <- subset(Navin_hbca_c74_panCK_IM, rownames(Navin_hbca_c74_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c74_panCK_IM <-t(Navin_hbca_c74_panCK_IM)
write.csv(Navin_hbca_c74_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c74_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c74_panCK_IM)

# Remove Object
rm(Navin_hbca_c74_panCK)
gc()



############################ 
# Navin_hbca_c75_panCK   #
############################ 
# Load Object
Navin_hbca_c75_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c75_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c75_panCK_All <- Navin_hbca_c75_panCK[["RNA"]]$counts
Navin_hbca_c75_panCK_All <- as.matrix(Navin_hbca_c75_panCK_All, 'sparseMatrix')
Navin_hbca_c75_panCK_All <-t(Navin_hbca_c75_panCK_All)
write.csv(Navin_hbca_c75_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c75_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c75_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c75_panCK_IM <- Navin_hbca_c75_panCK[["RNA"]]$counts
Navin_hbca_c75_panCK_IM <- as.matrix(Navin_hbca_c75_panCK_IM, 'sparseMatrix')
Navin_hbca_c75_panCK_IM <- subset(Navin_hbca_c75_panCK_IM, rownames(Navin_hbca_c75_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c75_panCK_IM <-t(Navin_hbca_c75_panCK_IM)
write.csv(Navin_hbca_c75_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c75_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c75_panCK_IM)

# Remove Object
rm(Navin_hbca_c75_panCK)
gc()



############################ 
# Navin_hbca_c76_panCK   #
############################ 
# Load Object
Navin_hbca_c76_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c76_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c76_panCK_All <- Navin_hbca_c76_panCK[["RNA"]]$counts
Navin_hbca_c76_panCK_All <- as.matrix(Navin_hbca_c76_panCK_All, 'sparseMatrix')
Navin_hbca_c76_panCK_All <-t(Navin_hbca_c76_panCK_All)
write.csv(Navin_hbca_c76_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c76_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c76_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c76_panCK_IM <- Navin_hbca_c76_panCK[["RNA"]]$counts
Navin_hbca_c76_panCK_IM <- as.matrix(Navin_hbca_c76_panCK_IM, 'sparseMatrix')
Navin_hbca_c76_panCK_IM <- subset(Navin_hbca_c76_panCK_IM, rownames(Navin_hbca_c76_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c76_panCK_IM <-t(Navin_hbca_c76_panCK_IM)
write.csv(Navin_hbca_c76_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c76_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c76_panCK_IM)

# Remove Object
rm(Navin_hbca_c76_panCK)
gc()



############################ 
# Navin_hbca_c77_panCK   #
############################ 
# Load Object
Navin_hbca_c77_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c77_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c77_panCK_All <- Navin_hbca_c77_panCK[["RNA"]]$counts
Navin_hbca_c77_panCK_All <- as.matrix(Navin_hbca_c77_panCK_All, 'sparseMatrix')
Navin_hbca_c77_panCK_All <-t(Navin_hbca_c77_panCK_All)
write.csv(Navin_hbca_c77_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c77_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c77_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c77_panCK_IM <- Navin_hbca_c77_panCK[["RNA"]]$counts
Navin_hbca_c77_panCK_IM <- as.matrix(Navin_hbca_c77_panCK_IM, 'sparseMatrix')
Navin_hbca_c77_panCK_IM <- subset(Navin_hbca_c77_panCK_IM, rownames(Navin_hbca_c77_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c77_panCK_IM <-t(Navin_hbca_c77_panCK_IM)
write.csv(Navin_hbca_c77_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c77_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c77_panCK_IM)

# Remove Object
rm(Navin_hbca_c77_panCK)
gc()



############################ 
# Navin_hbca_c78_panCK   #
############################ 
# Load Object
Navin_hbca_c78_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c78_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c78_panCK_All <- Navin_hbca_c78_panCK[["RNA"]]$counts
Navin_hbca_c78_panCK_All <- as.matrix(Navin_hbca_c78_panCK_All, 'sparseMatrix')
Navin_hbca_c78_panCK_All <-t(Navin_hbca_c78_panCK_All)
write.csv(Navin_hbca_c78_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c78_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c78_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c78_panCK_IM <- Navin_hbca_c78_panCK[["RNA"]]$counts
Navin_hbca_c78_panCK_IM <- as.matrix(Navin_hbca_c78_panCK_IM, 'sparseMatrix')
Navin_hbca_c78_panCK_IM <- subset(Navin_hbca_c78_panCK_IM, rownames(Navin_hbca_c78_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c78_panCK_IM <-t(Navin_hbca_c78_panCK_IM)
write.csv(Navin_hbca_c78_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c78_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c78_panCK_IM)

# Remove Object
rm(Navin_hbca_c78_panCK)
gc()





############################ 
# Navin_hbca_c79_panCK   #
############################ 
# Load Object
Navin_hbca_c79_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c79_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c79_panCK_All <- Navin_hbca_c79_panCK[["RNA"]]$counts
Navin_hbca_c79_panCK_All <- as.matrix(Navin_hbca_c79_panCK_All, 'sparseMatrix')
Navin_hbca_c79_panCK_All <-t(Navin_hbca_c79_panCK_All)
write.csv(Navin_hbca_c79_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c79_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c79_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c79_panCK_IM <- Navin_hbca_c79_panCK[["RNA"]]$counts
Navin_hbca_c79_panCK_IM <- as.matrix(Navin_hbca_c79_panCK_IM, 'sparseMatrix')
Navin_hbca_c79_panCK_IM <- subset(Navin_hbca_c79_panCK_IM, rownames(Navin_hbca_c79_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c79_panCK_IM <-t(Navin_hbca_c79_panCK_IM)
write.csv(Navin_hbca_c79_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c79_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c79_panCK_IM)

# Remove Object
rm(Navin_hbca_c79_panCK)
gc()



############################ 
# Navin_hbca_c80_panCK   #
############################ 
# Load Object
Navin_hbca_c80_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c80_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c80_panCK_All <- Navin_hbca_c80_panCK[["RNA"]]$counts
Navin_hbca_c80_panCK_All <- as.matrix(Navin_hbca_c80_panCK_All, 'sparseMatrix')
Navin_hbca_c80_panCK_All <-t(Navin_hbca_c80_panCK_All)
write.csv(Navin_hbca_c80_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c80_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c80_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c80_panCK_IM <- Navin_hbca_c80_panCK[["RNA"]]$counts
Navin_hbca_c80_panCK_IM <- as.matrix(Navin_hbca_c80_panCK_IM, 'sparseMatrix')
Navin_hbca_c80_panCK_IM <- subset(Navin_hbca_c80_panCK_IM, rownames(Navin_hbca_c80_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c80_panCK_IM <-t(Navin_hbca_c80_panCK_IM)
write.csv(Navin_hbca_c80_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c80_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c80_panCK_IM)

# Remove Object
rm(Navin_hbca_c80_panCK)
gc()



############################ 
# Navin_hbca_c81_panCK   #
############################ 
# Load Object
Navin_hbca_c81_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c81_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c81_panCK_All <- Navin_hbca_c81_panCK[["RNA"]]$counts
Navin_hbca_c81_panCK_All <- as.matrix(Navin_hbca_c81_panCK_All, 'sparseMatrix')
Navin_hbca_c81_panCK_All <-t(Navin_hbca_c81_panCK_All)
write.csv(Navin_hbca_c81_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c81_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c81_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c81_panCK_IM <- Navin_hbca_c81_panCK[["RNA"]]$counts
Navin_hbca_c81_panCK_IM <- as.matrix(Navin_hbca_c81_panCK_IM, 'sparseMatrix')
Navin_hbca_c81_panCK_IM <- subset(Navin_hbca_c81_panCK_IM, rownames(Navin_hbca_c81_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c81_panCK_IM <-t(Navin_hbca_c81_panCK_IM)
write.csv(Navin_hbca_c81_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c81_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c81_panCK_IM)

# Remove Object
rm(Navin_hbca_c81_panCK)
gc()



############################ 
# Navin_hbca_c82_panCK   #
############################ 
# Load Object
Navin_hbca_c82_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c82_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c82_panCK_All <- Navin_hbca_c82_panCK[["RNA"]]$counts
Navin_hbca_c82_panCK_All <- as.matrix(Navin_hbca_c82_panCK_All, 'sparseMatrix')
Navin_hbca_c82_panCK_All <-t(Navin_hbca_c82_panCK_All)
write.csv(Navin_hbca_c82_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c82_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c82_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c82_panCK_IM <- Navin_hbca_c82_panCK[["RNA"]]$counts
Navin_hbca_c82_panCK_IM <- as.matrix(Navin_hbca_c82_panCK_IM, 'sparseMatrix')
Navin_hbca_c82_panCK_IM <- subset(Navin_hbca_c82_panCK_IM, rownames(Navin_hbca_c82_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c82_panCK_IM <-t(Navin_hbca_c82_panCK_IM)
write.csv(Navin_hbca_c82_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c82_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c82_panCK_IM)

# Remove Object
rm(Navin_hbca_c82_panCK)
gc()



############################ 
# Navin_hbca_c83_panCK   #
############################ 
# Load Object
Navin_hbca_c83_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c83_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c83_panCK_All <- Navin_hbca_c83_panCK[["RNA"]]$counts
Navin_hbca_c83_panCK_All <- as.matrix(Navin_hbca_c83_panCK_All, 'sparseMatrix')
Navin_hbca_c83_panCK_All <-t(Navin_hbca_c83_panCK_All)
write.csv(Navin_hbca_c83_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c83_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c83_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c83_panCK_IM <- Navin_hbca_c83_panCK[["RNA"]]$counts
Navin_hbca_c83_panCK_IM <- as.matrix(Navin_hbca_c83_panCK_IM, 'sparseMatrix')
Navin_hbca_c83_panCK_IM <- subset(Navin_hbca_c83_panCK_IM, rownames(Navin_hbca_c83_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c83_panCK_IM <-t(Navin_hbca_c83_panCK_IM)
write.csv(Navin_hbca_c83_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c83_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c83_panCK_IM)

# Remove Object
rm(Navin_hbca_c83_panCK)
gc()



############################ 
# Navin_hbca_c84_panCK   #
############################ 
# Load Object
Navin_hbca_c84_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c84_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c84_panCK_All <- Navin_hbca_c84_panCK[["RNA"]]$counts
Navin_hbca_c84_panCK_All <- as.matrix(Navin_hbca_c84_panCK_All, 'sparseMatrix')
Navin_hbca_c84_panCK_All <-t(Navin_hbca_c84_panCK_All)
write.csv(Navin_hbca_c84_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c84_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c84_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c84_panCK_IM <- Navin_hbca_c84_panCK[["RNA"]]$counts
Navin_hbca_c84_panCK_IM <- as.matrix(Navin_hbca_c84_panCK_IM, 'sparseMatrix')
Navin_hbca_c84_panCK_IM <- subset(Navin_hbca_c84_panCK_IM, rownames(Navin_hbca_c84_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c84_panCK_IM <-t(Navin_hbca_c84_panCK_IM)
write.csv(Navin_hbca_c84_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c84_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c84_panCK_IM)

# Remove Object
rm(Navin_hbca_c84_panCK)
gc()



############################ 
# Navin_hbca_c85_panCK   #
############################ 
# Load Object
Navin_hbca_c85_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c85_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c85_panCK_All <- Navin_hbca_c85_panCK[["RNA"]]$counts
Navin_hbca_c85_panCK_All <- as.matrix(Navin_hbca_c85_panCK_All, 'sparseMatrix')
Navin_hbca_c85_panCK_All <-t(Navin_hbca_c85_panCK_All)
write.csv(Navin_hbca_c85_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c85_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c85_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c85_panCK_IM <- Navin_hbca_c85_panCK[["RNA"]]$counts
Navin_hbca_c85_panCK_IM <- as.matrix(Navin_hbca_c85_panCK_IM, 'sparseMatrix')
Navin_hbca_c85_panCK_IM <- subset(Navin_hbca_c85_panCK_IM, rownames(Navin_hbca_c85_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c85_panCK_IM <-t(Navin_hbca_c85_panCK_IM)
write.csv(Navin_hbca_c85_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c85_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c85_panCK_IM)

# Remove Object
rm(Navin_hbca_c85_panCK)
gc()



############################ 
# Navin_hbca_c86_panCK   #
############################ 
# Load Object
Navin_hbca_c86_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c86_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c86_panCK_All <- Navin_hbca_c86_panCK[["RNA"]]$counts
Navin_hbca_c86_panCK_All <- as.matrix(Navin_hbca_c86_panCK_All, 'sparseMatrix')
Navin_hbca_c86_panCK_All <-t(Navin_hbca_c86_panCK_All)
write.csv(Navin_hbca_c86_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c86_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c86_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c86_panCK_IM <- Navin_hbca_c86_panCK[["RNA"]]$counts
Navin_hbca_c86_panCK_IM <- as.matrix(Navin_hbca_c86_panCK_IM, 'sparseMatrix')
Navin_hbca_c86_panCK_IM <- subset(Navin_hbca_c86_panCK_IM, rownames(Navin_hbca_c86_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c86_panCK_IM <-t(Navin_hbca_c86_panCK_IM)
write.csv(Navin_hbca_c86_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c86_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c86_panCK_IM)

# Remove Object
rm(Navin_hbca_c86_panCK)
gc()



############################ 
# Navin_hbca_c87_panCK   #
############################ 
# Load Object
Navin_hbca_c87_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c87_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c87_panCK_All <- Navin_hbca_c87_panCK[["RNA"]]$counts
Navin_hbca_c87_panCK_All <- as.matrix(Navin_hbca_c87_panCK_All, 'sparseMatrix')
Navin_hbca_c87_panCK_All <-t(Navin_hbca_c87_panCK_All)
write.csv(Navin_hbca_c87_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c87_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c87_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c87_panCK_IM <- Navin_hbca_c87_panCK[["RNA"]]$counts
Navin_hbca_c87_panCK_IM <- as.matrix(Navin_hbca_c87_panCK_IM, 'sparseMatrix')
Navin_hbca_c87_panCK_IM <- subset(Navin_hbca_c87_panCK_IM, rownames(Navin_hbca_c87_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c87_panCK_IM <-t(Navin_hbca_c87_panCK_IM)
write.csv(Navin_hbca_c87_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c87_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c87_panCK_IM)

# Remove Object
rm(Navin_hbca_c87_panCK)
gc()



############################ 
# Navin_hbca_c88_panCK   #
############################ 
# Load Object
Navin_hbca_c88_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c88_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c88_panCK_All <- Navin_hbca_c88_panCK[["RNA"]]$counts
Navin_hbca_c88_panCK_All <- as.matrix(Navin_hbca_c88_panCK_All, 'sparseMatrix')
Navin_hbca_c88_panCK_All <-t(Navin_hbca_c88_panCK_All)
write.csv(Navin_hbca_c88_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c88_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c88_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c88_panCK_IM <- Navin_hbca_c88_panCK[["RNA"]]$counts
Navin_hbca_c88_panCK_IM <- as.matrix(Navin_hbca_c88_panCK_IM, 'sparseMatrix')
Navin_hbca_c88_panCK_IM <- subset(Navin_hbca_c88_panCK_IM, rownames(Navin_hbca_c88_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c88_panCK_IM <-t(Navin_hbca_c88_panCK_IM)
write.csv(Navin_hbca_c88_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c88_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c88_panCK_IM)

# Remove Object
rm(Navin_hbca_c88_panCK)
gc()





############################ 
# Navin_hbca_c89_panCK   #
############################ 
# Load Object
Navin_hbca_c89_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c89_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c89_panCK_All <- Navin_hbca_c89_panCK[["RNA"]]$counts
Navin_hbca_c89_panCK_All <- as.matrix(Navin_hbca_c89_panCK_All, 'sparseMatrix')
Navin_hbca_c89_panCK_All <-t(Navin_hbca_c89_panCK_All)
write.csv(Navin_hbca_c89_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c89_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c89_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c89_panCK_IM <- Navin_hbca_c89_panCK[["RNA"]]$counts
Navin_hbca_c89_panCK_IM <- as.matrix(Navin_hbca_c89_panCK_IM, 'sparseMatrix')
Navin_hbca_c89_panCK_IM <- subset(Navin_hbca_c89_panCK_IM, rownames(Navin_hbca_c89_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c89_panCK_IM <-t(Navin_hbca_c89_panCK_IM)
write.csv(Navin_hbca_c89_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c89_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c89_panCK_IM)

# Remove Object
rm(Navin_hbca_c89_panCK)
gc()



############################ 
# Navin_hbca_c90_panCK   #
############################ 
# Load Object
Navin_hbca_c90_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c90_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c90_panCK_All <- Navin_hbca_c90_panCK[["RNA"]]$counts
Navin_hbca_c90_panCK_All <- as.matrix(Navin_hbca_c90_panCK_All, 'sparseMatrix')
Navin_hbca_c90_panCK_All <-t(Navin_hbca_c90_panCK_All)
write.csv(Navin_hbca_c90_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c90_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c90_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c90_panCK_IM <- Navin_hbca_c90_panCK[["RNA"]]$counts
Navin_hbca_c90_panCK_IM <- as.matrix(Navin_hbca_c90_panCK_IM, 'sparseMatrix')
Navin_hbca_c90_panCK_IM <- subset(Navin_hbca_c90_panCK_IM, rownames(Navin_hbca_c90_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c90_panCK_IM <-t(Navin_hbca_c90_panCK_IM)
write.csv(Navin_hbca_c90_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c90_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c90_panCK_IM)

# Remove Object
rm(Navin_hbca_c90_panCK)
gc()



############################ 
# Navin_hbca_c91_panCK   #
############################ 
# Load Object
Navin_hbca_c91_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c91_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c91_panCK_All <- Navin_hbca_c91_panCK[["RNA"]]$counts
Navin_hbca_c91_panCK_All <- as.matrix(Navin_hbca_c91_panCK_All, 'sparseMatrix')
Navin_hbca_c91_panCK_All <-t(Navin_hbca_c91_panCK_All)
write.csv(Navin_hbca_c91_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c91_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c91_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c91_panCK_IM <- Navin_hbca_c91_panCK[["RNA"]]$counts
Navin_hbca_c91_panCK_IM <- as.matrix(Navin_hbca_c91_panCK_IM, 'sparseMatrix')
Navin_hbca_c91_panCK_IM <- subset(Navin_hbca_c91_panCK_IM, rownames(Navin_hbca_c91_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c91_panCK_IM <-t(Navin_hbca_c91_panCK_IM)
write.csv(Navin_hbca_c91_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c91_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c91_panCK_IM)

# Remove Object
rm(Navin_hbca_c91_panCK)
gc()



############################ 
# Navin_hbca_c92_panCK   #
############################ 
# Load Object
Navin_hbca_c92_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c92_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c92_panCK_All <- Navin_hbca_c92_panCK[["RNA"]]$counts
Navin_hbca_c92_panCK_All <- as.matrix(Navin_hbca_c92_panCK_All, 'sparseMatrix')
Navin_hbca_c92_panCK_All <-t(Navin_hbca_c92_panCK_All)
write.csv(Navin_hbca_c92_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c92_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c92_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c92_panCK_IM <- Navin_hbca_c92_panCK[["RNA"]]$counts
Navin_hbca_c92_panCK_IM <- as.matrix(Navin_hbca_c92_panCK_IM, 'sparseMatrix')
Navin_hbca_c92_panCK_IM <- subset(Navin_hbca_c92_panCK_IM, rownames(Navin_hbca_c92_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c92_panCK_IM <-t(Navin_hbca_c92_panCK_IM)
write.csv(Navin_hbca_c92_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c92_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c92_panCK_IM)

# Remove Object
rm(Navin_hbca_c92_panCK)
gc()



############################ 
# Navin_hbca_c93_panCK   #
############################ 
# Load Object
Navin_hbca_c93_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c93_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c93_panCK_All <- Navin_hbca_c93_panCK[["RNA"]]$counts
Navin_hbca_c93_panCK_All <- as.matrix(Navin_hbca_c93_panCK_All, 'sparseMatrix')
Navin_hbca_c93_panCK_All <-t(Navin_hbca_c93_panCK_All)
write.csv(Navin_hbca_c93_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c93_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c93_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c93_panCK_IM <- Navin_hbca_c93_panCK[["RNA"]]$counts
Navin_hbca_c93_panCK_IM <- as.matrix(Navin_hbca_c93_panCK_IM, 'sparseMatrix')
Navin_hbca_c93_panCK_IM <- subset(Navin_hbca_c93_panCK_IM, rownames(Navin_hbca_c93_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c93_panCK_IM <-t(Navin_hbca_c93_panCK_IM)
write.csv(Navin_hbca_c93_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c93_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c93_panCK_IM)

# Remove Object
rm(Navin_hbca_c93_panCK)
gc()



############################ 
# Navin_hbca_c94_panCK   #
############################ 
# Load Object
Navin_hbca_c94_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c94_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c94_panCK_All <- Navin_hbca_c94_panCK[["RNA"]]$counts
Navin_hbca_c94_panCK_All <- as.matrix(Navin_hbca_c94_panCK_All, 'sparseMatrix')
Navin_hbca_c94_panCK_All <-t(Navin_hbca_c94_panCK_All)
write.csv(Navin_hbca_c94_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c94_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c94_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c94_panCK_IM <- Navin_hbca_c94_panCK[["RNA"]]$counts
Navin_hbca_c94_panCK_IM <- as.matrix(Navin_hbca_c94_panCK_IM, 'sparseMatrix')
Navin_hbca_c94_panCK_IM <- subset(Navin_hbca_c94_panCK_IM, rownames(Navin_hbca_c94_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c94_panCK_IM <-t(Navin_hbca_c94_panCK_IM)
write.csv(Navin_hbca_c94_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c94_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c94_panCK_IM)

# Remove Object
rm(Navin_hbca_c94_panCK)
gc()



############################ 
# Navin_hbca_c95_panCK   #
############################ 
# Load Object
Navin_hbca_c95_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c95_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c95_panCK_All <- Navin_hbca_c95_panCK[["RNA"]]$counts
Navin_hbca_c95_panCK_All <- as.matrix(Navin_hbca_c95_panCK_All, 'sparseMatrix')
Navin_hbca_c95_panCK_All <-t(Navin_hbca_c95_panCK_All)
write.csv(Navin_hbca_c95_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c95_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c95_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c95_panCK_IM <- Navin_hbca_c95_panCK[["RNA"]]$counts
Navin_hbca_c95_panCK_IM <- as.matrix(Navin_hbca_c95_panCK_IM, 'sparseMatrix')
Navin_hbca_c95_panCK_IM <- subset(Navin_hbca_c95_panCK_IM, rownames(Navin_hbca_c95_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c95_panCK_IM <-t(Navin_hbca_c95_panCK_IM)
write.csv(Navin_hbca_c95_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c95_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c95_panCK_IM)

# Remove Object
rm(Navin_hbca_c95_panCK)
gc()



############################ 
# Navin_hbca_c96_panCK   #
############################ 
# Load Object
Navin_hbca_c96_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c96_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c96_panCK_All <- Navin_hbca_c96_panCK[["RNA"]]$counts
Navin_hbca_c96_panCK_All <- as.matrix(Navin_hbca_c96_panCK_All, 'sparseMatrix')
Navin_hbca_c96_panCK_All <-t(Navin_hbca_c96_panCK_All)
write.csv(Navin_hbca_c96_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c96_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c96_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c96_panCK_IM <- Navin_hbca_c96_panCK[["RNA"]]$counts
Navin_hbca_c96_panCK_IM <- as.matrix(Navin_hbca_c96_panCK_IM, 'sparseMatrix')
Navin_hbca_c96_panCK_IM <- subset(Navin_hbca_c96_panCK_IM, rownames(Navin_hbca_c96_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c96_panCK_IM <-t(Navin_hbca_c96_panCK_IM)
write.csv(Navin_hbca_c96_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c96_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c96_panCK_IM)

# Remove Object
rm(Navin_hbca_c96_panCK)
gc()



############################ 
# Navin_hbca_c97_panCK   #
############################ 
# Load Object
Navin_hbca_c97_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c97_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c97_panCK_All <- Navin_hbca_c97_panCK[["RNA"]]$counts
Navin_hbca_c97_panCK_All <- as.matrix(Navin_hbca_c97_panCK_All, 'sparseMatrix')
Navin_hbca_c97_panCK_All <-t(Navin_hbca_c97_panCK_All)
write.csv(Navin_hbca_c97_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c97_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c97_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c97_panCK_IM <- Navin_hbca_c97_panCK[["RNA"]]$counts
Navin_hbca_c97_panCK_IM <- as.matrix(Navin_hbca_c97_panCK_IM, 'sparseMatrix')
Navin_hbca_c97_panCK_IM <- subset(Navin_hbca_c97_panCK_IM, rownames(Navin_hbca_c97_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c97_panCK_IM <-t(Navin_hbca_c97_panCK_IM)
write.csv(Navin_hbca_c97_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c97_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c97_panCK_IM)

# Remove Object
rm(Navin_hbca_c97_panCK)
gc()



############################ 
# Navin_hbca_c98_panCK   #
############################ 
# Load Object
Navin_hbca_c98_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c98_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c98_panCK_All <- Navin_hbca_c98_panCK[["RNA"]]$counts
Navin_hbca_c98_panCK_All <- as.matrix(Navin_hbca_c98_panCK_All, 'sparseMatrix')
Navin_hbca_c98_panCK_All <-t(Navin_hbca_c98_panCK_All)
write.csv(Navin_hbca_c98_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c98_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c98_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c98_panCK_IM <- Navin_hbca_c98_panCK[["RNA"]]$counts
Navin_hbca_c98_panCK_IM <- as.matrix(Navin_hbca_c98_panCK_IM, 'sparseMatrix')
Navin_hbca_c98_panCK_IM <- subset(Navin_hbca_c98_panCK_IM, rownames(Navin_hbca_c98_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c98_panCK_IM <-t(Navin_hbca_c98_panCK_IM)
write.csv(Navin_hbca_c98_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c98_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c98_panCK_IM)

# Remove Object
rm(Navin_hbca_c98_panCK)
gc()





############################ 
# Navin_hbca_c99_panCK   #
############################ 
# Load Object
Navin_hbca_c99_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c99_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c99_panCK_All <- Navin_hbca_c99_panCK[["RNA"]]$counts
Navin_hbca_c99_panCK_All <- as.matrix(Navin_hbca_c99_panCK_All, 'sparseMatrix')
Navin_hbca_c99_panCK_All <-t(Navin_hbca_c99_panCK_All)
write.csv(Navin_hbca_c99_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c99_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c99_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c99_panCK_IM <- Navin_hbca_c99_panCK[["RNA"]]$counts
Navin_hbca_c99_panCK_IM <- as.matrix(Navin_hbca_c99_panCK_IM, 'sparseMatrix')
Navin_hbca_c99_panCK_IM <- subset(Navin_hbca_c99_panCK_IM, rownames(Navin_hbca_c99_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c99_panCK_IM <-t(Navin_hbca_c99_panCK_IM)
write.csv(Navin_hbca_c99_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c99_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c99_panCK_IM)

# Remove Object
rm(Navin_hbca_c99_panCK)
gc()



############################ 
# Navin_hbca_c100_panCK   #
############################ 
# Load Object
Navin_hbca_c100_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c100_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c100_panCK_All <- Navin_hbca_c100_panCK[["RNA"]]$counts
Navin_hbca_c100_panCK_All <- as.matrix(Navin_hbca_c100_panCK_All, 'sparseMatrix')
Navin_hbca_c100_panCK_All <-t(Navin_hbca_c100_panCK_All)
write.csv(Navin_hbca_c100_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c100_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c100_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c100_panCK_IM <- Navin_hbca_c100_panCK[["RNA"]]$counts
Navin_hbca_c100_panCK_IM <- as.matrix(Navin_hbca_c100_panCK_IM, 'sparseMatrix')
Navin_hbca_c100_panCK_IM <- subset(Navin_hbca_c100_panCK_IM, rownames(Navin_hbca_c100_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c100_panCK_IM <-t(Navin_hbca_c100_panCK_IM)
write.csv(Navin_hbca_c100_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c100_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c100_panCK_IM)

# Remove Object
rm(Navin_hbca_c100_panCK)
gc()



############################ 
# Navin_hbca_c101_panCK   #
############################ 
# Load Object
Navin_hbca_c101_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c101_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c101_panCK_All <- Navin_hbca_c101_panCK[["RNA"]]$counts
Navin_hbca_c101_panCK_All <- as.matrix(Navin_hbca_c101_panCK_All, 'sparseMatrix')
Navin_hbca_c101_panCK_All <-t(Navin_hbca_c101_panCK_All)
write.csv(Navin_hbca_c101_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c101_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c101_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c101_panCK_IM <- Navin_hbca_c101_panCK[["RNA"]]$counts
Navin_hbca_c101_panCK_IM <- as.matrix(Navin_hbca_c101_panCK_IM, 'sparseMatrix')
Navin_hbca_c101_panCK_IM <- subset(Navin_hbca_c101_panCK_IM, rownames(Navin_hbca_c101_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c101_panCK_IM <-t(Navin_hbca_c101_panCK_IM)
write.csv(Navin_hbca_c101_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c101_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c101_panCK_IM)

# Remove Object
rm(Navin_hbca_c101_panCK)
gc()



############################ 
# Navin_hbca_c102_panCK   #
############################ 
# Load Object
Navin_hbca_c102_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c102_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c102_panCK_All <- Navin_hbca_c102_panCK[["RNA"]]$counts
Navin_hbca_c102_panCK_All <- as.matrix(Navin_hbca_c102_panCK_All, 'sparseMatrix')
Navin_hbca_c102_panCK_All <-t(Navin_hbca_c102_panCK_All)
write.csv(Navin_hbca_c102_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c102_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c102_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c102_panCK_IM <- Navin_hbca_c102_panCK[["RNA"]]$counts
Navin_hbca_c102_panCK_IM <- as.matrix(Navin_hbca_c102_panCK_IM, 'sparseMatrix')
Navin_hbca_c102_panCK_IM <- subset(Navin_hbca_c102_panCK_IM, rownames(Navin_hbca_c102_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c102_panCK_IM <-t(Navin_hbca_c102_panCK_IM)
write.csv(Navin_hbca_c102_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c102_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c102_panCK_IM)

# Remove Object
rm(Navin_hbca_c102_panCK)
gc()



############################ 
# Navin_hbca_c103_panCK   #
############################ 
# Load Object
Navin_hbca_c103_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c103_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c103_panCK_All <- Navin_hbca_c103_panCK[["RNA"]]$counts
Navin_hbca_c103_panCK_All <- as.matrix(Navin_hbca_c103_panCK_All, 'sparseMatrix')
Navin_hbca_c103_panCK_All <-t(Navin_hbca_c103_panCK_All)
write.csv(Navin_hbca_c103_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c103_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c103_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c103_panCK_IM <- Navin_hbca_c103_panCK[["RNA"]]$counts
Navin_hbca_c103_panCK_IM <- as.matrix(Navin_hbca_c103_panCK_IM, 'sparseMatrix')
Navin_hbca_c103_panCK_IM <- subset(Navin_hbca_c103_panCK_IM, rownames(Navin_hbca_c103_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c103_panCK_IM <-t(Navin_hbca_c103_panCK_IM)
write.csv(Navin_hbca_c103_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c103_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c103_panCK_IM)

# Remove Object
rm(Navin_hbca_c103_panCK)
gc()



############################ 
# Navin_hbca_c104_panCK   #
############################ 
# Load Object
Navin_hbca_c104_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c104_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c104_panCK_All <- Navin_hbca_c104_panCK[["RNA"]]$counts
Navin_hbca_c104_panCK_All <- as.matrix(Navin_hbca_c104_panCK_All, 'sparseMatrix')
Navin_hbca_c104_panCK_All <-t(Navin_hbca_c104_panCK_All)
write.csv(Navin_hbca_c104_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c104_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c104_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c104_panCK_IM <- Navin_hbca_c104_panCK[["RNA"]]$counts
Navin_hbca_c104_panCK_IM <- as.matrix(Navin_hbca_c104_panCK_IM, 'sparseMatrix')
Navin_hbca_c104_panCK_IM <- subset(Navin_hbca_c104_panCK_IM, rownames(Navin_hbca_c104_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c104_panCK_IM <-t(Navin_hbca_c104_panCK_IM)
write.csv(Navin_hbca_c104_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c104_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c104_panCK_IM)

# Remove Object
rm(Navin_hbca_c104_panCK)
gc()



############################ 
# Navin_hbca_c105_panCK   #
############################ 
# Load Object
Navin_hbca_c105_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c105_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c105_panCK_All <- Navin_hbca_c105_panCK[["RNA"]]$counts
Navin_hbca_c105_panCK_All <- as.matrix(Navin_hbca_c105_panCK_All, 'sparseMatrix')
Navin_hbca_c105_panCK_All <-t(Navin_hbca_c105_panCK_All)
write.csv(Navin_hbca_c105_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c105_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c105_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c105_panCK_IM <- Navin_hbca_c105_panCK[["RNA"]]$counts
Navin_hbca_c105_panCK_IM <- as.matrix(Navin_hbca_c105_panCK_IM, 'sparseMatrix')
Navin_hbca_c105_panCK_IM <- subset(Navin_hbca_c105_panCK_IM, rownames(Navin_hbca_c105_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c105_panCK_IM <-t(Navin_hbca_c105_panCK_IM)
write.csv(Navin_hbca_c105_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c105_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c105_panCK_IM)

# Remove Object
rm(Navin_hbca_c105_panCK)
gc()



############################ 
# Navin_hbca_c106_panCK   #
############################ 
# Load Object
Navin_hbca_c106_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c106_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c106_panCK_All <- Navin_hbca_c106_panCK[["RNA"]]$counts
Navin_hbca_c106_panCK_All <- as.matrix(Navin_hbca_c106_panCK_All, 'sparseMatrix')
Navin_hbca_c106_panCK_All <-t(Navin_hbca_c106_panCK_All)
write.csv(Navin_hbca_c106_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c106_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c106_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c106_panCK_IM <- Navin_hbca_c106_panCK[["RNA"]]$counts
Navin_hbca_c106_panCK_IM <- as.matrix(Navin_hbca_c106_panCK_IM, 'sparseMatrix')
Navin_hbca_c106_panCK_IM <- subset(Navin_hbca_c106_panCK_IM, rownames(Navin_hbca_c106_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c106_panCK_IM <-t(Navin_hbca_c106_panCK_IM)
write.csv(Navin_hbca_c106_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c106_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c106_panCK_IM)

# Remove Object
rm(Navin_hbca_c106_panCK)
gc()



############################ 
# Navin_hbca_c107_panCK   #
############################ 
# Load Object
Navin_hbca_c107_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c107_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c107_panCK_All <- Navin_hbca_c107_panCK[["RNA"]]$counts
Navin_hbca_c107_panCK_All <- as.matrix(Navin_hbca_c107_panCK_All, 'sparseMatrix')
Navin_hbca_c107_panCK_All <-t(Navin_hbca_c107_panCK_All)
write.csv(Navin_hbca_c107_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c107_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c107_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c107_panCK_IM <- Navin_hbca_c107_panCK[["RNA"]]$counts
Navin_hbca_c107_panCK_IM <- as.matrix(Navin_hbca_c107_panCK_IM, 'sparseMatrix')
Navin_hbca_c107_panCK_IM <- subset(Navin_hbca_c107_panCK_IM, rownames(Navin_hbca_c107_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c107_panCK_IM <-t(Navin_hbca_c107_panCK_IM)
write.csv(Navin_hbca_c107_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c107_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c107_panCK_IM)

# Remove Object
rm(Navin_hbca_c107_panCK)
gc()



############################ 
# Navin_hbca_c108_panCK   #
############################ 
# Load Object
Navin_hbca_c108_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c108_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c108_panCK_All <- Navin_hbca_c108_panCK[["RNA"]]$counts
Navin_hbca_c108_panCK_All <- as.matrix(Navin_hbca_c108_panCK_All, 'sparseMatrix')
Navin_hbca_c108_panCK_All <-t(Navin_hbca_c108_panCK_All)
write.csv(Navin_hbca_c108_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c108_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c108_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c108_panCK_IM <- Navin_hbca_c108_panCK[["RNA"]]$counts
Navin_hbca_c108_panCK_IM <- as.matrix(Navin_hbca_c108_panCK_IM, 'sparseMatrix')
Navin_hbca_c108_panCK_IM <- subset(Navin_hbca_c108_panCK_IM, rownames(Navin_hbca_c108_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c108_panCK_IM <-t(Navin_hbca_c108_panCK_IM)
write.csv(Navin_hbca_c108_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c108_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c108_panCK_IM)

# Remove Object
rm(Navin_hbca_c108_panCK)
gc()





############################ 
# Navin_hbca_c109_panCK   #
############################ 
# Load Object
Navin_hbca_c109_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c109_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c109_panCK_All <- Navin_hbca_c109_panCK[["RNA"]]$counts
Navin_hbca_c109_panCK_All <- as.matrix(Navin_hbca_c109_panCK_All, 'sparseMatrix')
Navin_hbca_c109_panCK_All <-t(Navin_hbca_c109_panCK_All)
write.csv(Navin_hbca_c109_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c109_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c109_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c109_panCK_IM <- Navin_hbca_c109_panCK[["RNA"]]$counts
Navin_hbca_c109_panCK_IM <- as.matrix(Navin_hbca_c109_panCK_IM, 'sparseMatrix')
Navin_hbca_c109_panCK_IM <- subset(Navin_hbca_c109_panCK_IM, rownames(Navin_hbca_c109_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c109_panCK_IM <-t(Navin_hbca_c109_panCK_IM)
write.csv(Navin_hbca_c109_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c109_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c109_panCK_IM)

# Remove Object
rm(Navin_hbca_c109_panCK)
gc()



############################ 
# Navin_hbca_c110_panCK   #
############################ 
# Load Object
Navin_hbca_c110_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c110_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c110_panCK_All <- Navin_hbca_c110_panCK[["RNA"]]$counts
Navin_hbca_c110_panCK_All <- as.matrix(Navin_hbca_c110_panCK_All, 'sparseMatrix')
Navin_hbca_c110_panCK_All <-t(Navin_hbca_c110_panCK_All)
write.csv(Navin_hbca_c110_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c110_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c110_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c110_panCK_IM <- Navin_hbca_c110_panCK[["RNA"]]$counts
Navin_hbca_c110_panCK_IM <- as.matrix(Navin_hbca_c110_panCK_IM, 'sparseMatrix')
Navin_hbca_c110_panCK_IM <- subset(Navin_hbca_c110_panCK_IM, rownames(Navin_hbca_c110_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c110_panCK_IM <-t(Navin_hbca_c110_panCK_IM)
write.csv(Navin_hbca_c110_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c110_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c110_panCK_IM)

# Remove Object
rm(Navin_hbca_c110_panCK)
gc()



############################ 
# Navin_hbca_c111_panCK   #
############################ 
# Load Object
Navin_hbca_c111_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c111_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c111_panCK_All <- Navin_hbca_c111_panCK[["RNA"]]$counts
Navin_hbca_c111_panCK_All <- as.matrix(Navin_hbca_c111_panCK_All, 'sparseMatrix')
Navin_hbca_c111_panCK_All <-t(Navin_hbca_c111_panCK_All)
write.csv(Navin_hbca_c111_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c111_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c111_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c111_panCK_IM <- Navin_hbca_c111_panCK[["RNA"]]$counts
Navin_hbca_c111_panCK_IM <- as.matrix(Navin_hbca_c111_panCK_IM, 'sparseMatrix')
Navin_hbca_c111_panCK_IM <- subset(Navin_hbca_c111_panCK_IM, rownames(Navin_hbca_c111_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c111_panCK_IM <-t(Navin_hbca_c111_panCK_IM)
write.csv(Navin_hbca_c111_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c111_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c111_panCK_IM)

# Remove Object
rm(Navin_hbca_c111_panCK)
gc()



############################ 
# Navin_hbca_c113_panCK   #
############################ 
# Load Object
Navin_hbca_c113_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c113_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c113_panCK_All <- Navin_hbca_c113_panCK[["RNA"]]$counts
Navin_hbca_c113_panCK_All <- as.matrix(Navin_hbca_c113_panCK_All, 'sparseMatrix')
Navin_hbca_c113_panCK_All <-t(Navin_hbca_c113_panCK_All)
write.csv(Navin_hbca_c113_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c113_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c113_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c113_panCK_IM <- Navin_hbca_c113_panCK[["RNA"]]$counts
Navin_hbca_c113_panCK_IM <- as.matrix(Navin_hbca_c113_panCK_IM, 'sparseMatrix')
Navin_hbca_c113_panCK_IM <- subset(Navin_hbca_c113_panCK_IM, rownames(Navin_hbca_c113_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c113_panCK_IM <-t(Navin_hbca_c113_panCK_IM)
write.csv(Navin_hbca_c113_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c113_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c113_panCK_IM)

# Remove Object
rm(Navin_hbca_c113_panCK)
gc()



############################ 
# Navin_hbca_c114_panCK   #
############################ 
# Load Object
Navin_hbca_c114_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c114_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c114_panCK_All <- Navin_hbca_c114_panCK[["RNA"]]$counts
Navin_hbca_c114_panCK_All <- as.matrix(Navin_hbca_c114_panCK_All, 'sparseMatrix')
Navin_hbca_c114_panCK_All <-t(Navin_hbca_c114_panCK_All)
write.csv(Navin_hbca_c114_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c114_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c114_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c114_panCK_IM <- Navin_hbca_c114_panCK[["RNA"]]$counts
Navin_hbca_c114_panCK_IM <- as.matrix(Navin_hbca_c114_panCK_IM, 'sparseMatrix')
Navin_hbca_c114_panCK_IM <- subset(Navin_hbca_c114_panCK_IM, rownames(Navin_hbca_c114_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c114_panCK_IM <-t(Navin_hbca_c114_panCK_IM)
write.csv(Navin_hbca_c114_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c114_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c114_panCK_IM)

# Remove Object
rm(Navin_hbca_c114_panCK)
gc()



############################ 
# Navin_hbca_c115_panCK   #
############################ 
# Load Object
Navin_hbca_c115_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c115_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c115_panCK_All <- Navin_hbca_c115_panCK[["RNA"]]$counts
Navin_hbca_c115_panCK_All <- as.matrix(Navin_hbca_c115_panCK_All, 'sparseMatrix')
Navin_hbca_c115_panCK_All <-t(Navin_hbca_c115_panCK_All)
write.csv(Navin_hbca_c115_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c115_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c115_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c115_panCK_IM <- Navin_hbca_c115_panCK[["RNA"]]$counts
Navin_hbca_c115_panCK_IM <- as.matrix(Navin_hbca_c115_panCK_IM, 'sparseMatrix')
Navin_hbca_c115_panCK_IM <- subset(Navin_hbca_c115_panCK_IM, rownames(Navin_hbca_c115_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c115_panCK_IM <-t(Navin_hbca_c115_panCK_IM)
write.csv(Navin_hbca_c115_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c115_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c115_panCK_IM)

# Remove Object
rm(Navin_hbca_c115_panCK)
gc()



############################ 
# Navin_hbca_c118_panCK   #
############################ 
# Load Object
Navin_hbca_c118_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c118_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c118_panCK_All <- Navin_hbca_c118_panCK[["RNA"]]$counts
Navin_hbca_c118_panCK_All <- as.matrix(Navin_hbca_c118_panCK_All, 'sparseMatrix')
Navin_hbca_c118_panCK_All <-t(Navin_hbca_c118_panCK_All)
write.csv(Navin_hbca_c118_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c118_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c118_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c118_panCK_IM <- Navin_hbca_c118_panCK[["RNA"]]$counts
Navin_hbca_c118_panCK_IM <- as.matrix(Navin_hbca_c118_panCK_IM, 'sparseMatrix')
Navin_hbca_c118_panCK_IM <- subset(Navin_hbca_c118_panCK_IM, rownames(Navin_hbca_c118_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c118_panCK_IM <-t(Navin_hbca_c118_panCK_IM)
write.csv(Navin_hbca_c118_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c118_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c118_panCK_IM)

# Remove Object
rm(Navin_hbca_c118_panCK)
gc()



############################ 
# Navin_hbca_c119_panCK   #
############################ 
# Load Object
Navin_hbca_c119_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c119_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c119_panCK_All <- Navin_hbca_c119_panCK[["RNA"]]$counts
Navin_hbca_c119_panCK_All <- as.matrix(Navin_hbca_c119_panCK_All, 'sparseMatrix')
Navin_hbca_c119_panCK_All <-t(Navin_hbca_c119_panCK_All)
write.csv(Navin_hbca_c119_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c119_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c119_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c119_panCK_IM <- Navin_hbca_c119_panCK[["RNA"]]$counts
Navin_hbca_c119_panCK_IM <- as.matrix(Navin_hbca_c119_panCK_IM, 'sparseMatrix')
Navin_hbca_c119_panCK_IM <- subset(Navin_hbca_c119_panCK_IM, rownames(Navin_hbca_c119_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c119_panCK_IM <-t(Navin_hbca_c119_panCK_IM)
write.csv(Navin_hbca_c119_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c119_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c119_panCK_IM)

# Remove Object
rm(Navin_hbca_c119_panCK)
gc()



############################ 
# Navin_hbca_c121_panCK   #
############################ 
# Load Object
Navin_hbca_c121_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c121_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c121_panCK_All <- Navin_hbca_c121_panCK[["RNA"]]$counts
Navin_hbca_c121_panCK_All <- as.matrix(Navin_hbca_c121_panCK_All, 'sparseMatrix')
Navin_hbca_c121_panCK_All <-t(Navin_hbca_c121_panCK_All)
write.csv(Navin_hbca_c121_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c121_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c121_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c121_panCK_IM <- Navin_hbca_c121_panCK[["RNA"]]$counts
Navin_hbca_c121_panCK_IM <- as.matrix(Navin_hbca_c121_panCK_IM, 'sparseMatrix')
Navin_hbca_c121_panCK_IM <- subset(Navin_hbca_c121_panCK_IM, rownames(Navin_hbca_c121_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c121_panCK_IM <-t(Navin_hbca_c121_panCK_IM)
write.csv(Navin_hbca_c121_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c121_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c121_panCK_IM)

# Remove Object
rm(Navin_hbca_c121_panCK)
gc()



############################ 
# Navin_hbca_c122_panCK   #
############################ 
# Load Object
Navin_hbca_c122_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c122_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c122_panCK_All <- Navin_hbca_c122_panCK[["RNA"]]$counts
Navin_hbca_c122_panCK_All <- as.matrix(Navin_hbca_c122_panCK_All, 'sparseMatrix')
Navin_hbca_c122_panCK_All <-t(Navin_hbca_c122_panCK_All)
write.csv(Navin_hbca_c122_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c122_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c122_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c122_panCK_IM <- Navin_hbca_c122_panCK[["RNA"]]$counts
Navin_hbca_c122_panCK_IM <- as.matrix(Navin_hbca_c122_panCK_IM, 'sparseMatrix')
Navin_hbca_c122_panCK_IM <- subset(Navin_hbca_c122_panCK_IM, rownames(Navin_hbca_c122_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c122_panCK_IM <-t(Navin_hbca_c122_panCK_IM)
write.csv(Navin_hbca_c122_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c122_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c122_panCK_IM)

# Remove Object
rm(Navin_hbca_c122_panCK)
gc()





############################ 
# Navin_hbca_c123_panCK   #
############################ 
# Load Object
Navin_hbca_c123_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c123_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c123_panCK_All <- Navin_hbca_c123_panCK[["RNA"]]$counts
Navin_hbca_c123_panCK_All <- as.matrix(Navin_hbca_c123_panCK_All, 'sparseMatrix')
Navin_hbca_c123_panCK_All <-t(Navin_hbca_c123_panCK_All)
write.csv(Navin_hbca_c123_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c123_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c123_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c123_panCK_IM <- Navin_hbca_c123_panCK[["RNA"]]$counts
Navin_hbca_c123_panCK_IM <- as.matrix(Navin_hbca_c123_panCK_IM, 'sparseMatrix')
Navin_hbca_c123_panCK_IM <- subset(Navin_hbca_c123_panCK_IM, rownames(Navin_hbca_c123_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c123_panCK_IM <-t(Navin_hbca_c123_panCK_IM)
write.csv(Navin_hbca_c123_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c123_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c123_panCK_IM)

# Remove Object
rm(Navin_hbca_c123_panCK)
gc()



############################ 
# Navin_hbca_c124_panCK   #
############################ 
# Load Object
Navin_hbca_c124_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c124_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c124_panCK_All <- Navin_hbca_c124_panCK[["RNA"]]$counts
Navin_hbca_c124_panCK_All <- as.matrix(Navin_hbca_c124_panCK_All, 'sparseMatrix')
Navin_hbca_c124_panCK_All <-t(Navin_hbca_c124_panCK_All)
write.csv(Navin_hbca_c124_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c124_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c124_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c124_panCK_IM <- Navin_hbca_c124_panCK[["RNA"]]$counts
Navin_hbca_c124_panCK_IM <- as.matrix(Navin_hbca_c124_panCK_IM, 'sparseMatrix')
Navin_hbca_c124_panCK_IM <- subset(Navin_hbca_c124_panCK_IM, rownames(Navin_hbca_c124_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c124_panCK_IM <-t(Navin_hbca_c124_panCK_IM)
write.csv(Navin_hbca_c124_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c124_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c124_panCK_IM)

# Remove Object
rm(Navin_hbca_c124_panCK)
gc()



############################ 
# Navin_hbca_c125_panCK   #
############################ 
# Load Object
Navin_hbca_c125_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c125_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c125_panCK_All <- Navin_hbca_c125_panCK[["RNA"]]$counts
Navin_hbca_c125_panCK_All <- as.matrix(Navin_hbca_c125_panCK_All, 'sparseMatrix')
Navin_hbca_c125_panCK_All <-t(Navin_hbca_c125_panCK_All)
write.csv(Navin_hbca_c125_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c125_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c125_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c125_panCK_IM <- Navin_hbca_c125_panCK[["RNA"]]$counts
Navin_hbca_c125_panCK_IM <- as.matrix(Navin_hbca_c125_panCK_IM, 'sparseMatrix')
Navin_hbca_c125_panCK_IM <- subset(Navin_hbca_c125_panCK_IM, rownames(Navin_hbca_c125_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c125_panCK_IM <-t(Navin_hbca_c125_panCK_IM)
write.csv(Navin_hbca_c125_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c125_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c125_panCK_IM)

# Remove Object
rm(Navin_hbca_c125_panCK)
gc()



############################ 
# Navin_hbca_c126_panCK   #
############################ 
# Load Object
Navin_hbca_c126_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c126_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c126_panCK_All <- Navin_hbca_c126_panCK[["RNA"]]$counts
Navin_hbca_c126_panCK_All <- as.matrix(Navin_hbca_c126_panCK_All, 'sparseMatrix')
Navin_hbca_c126_panCK_All <-t(Navin_hbca_c126_panCK_All)
write.csv(Navin_hbca_c126_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c126_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c126_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c126_panCK_IM <- Navin_hbca_c126_panCK[["RNA"]]$counts
Navin_hbca_c126_panCK_IM <- as.matrix(Navin_hbca_c126_panCK_IM, 'sparseMatrix')
Navin_hbca_c126_panCK_IM <- subset(Navin_hbca_c126_panCK_IM, rownames(Navin_hbca_c126_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c126_panCK_IM <-t(Navin_hbca_c126_panCK_IM)
write.csv(Navin_hbca_c126_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c126_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c126_panCK_IM)

# Remove Object
rm(Navin_hbca_c126_panCK)
gc()



############################# 
# Navin_hbca_c127_panCK   #
############################# 
# Load Object
Navin_hbca_c127_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c127_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c127_panCK_All <- Navin_hbca_c127_panCK[["RNA"]]$counts
Navin_hbca_c127_panCK_All <- as.matrix(Navin_hbca_c127_panCK_All, 'sparseMatrix')
Navin_hbca_c127_panCK_All <-t(Navin_hbca_c127_panCK_All)
write.csv(Navin_hbca_c127_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c127_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c127_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c127_panCK_IM <- Navin_hbca_c127_panCK[["RNA"]]$counts
Navin_hbca_c127_panCK_IM <- as.matrix(Navin_hbca_c127_panCK_IM, 'sparseMatrix')
Navin_hbca_c127_panCK_IM <- subset(Navin_hbca_c127_panCK_IM, rownames(Navin_hbca_c127_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c127_panCK_IM <-t(Navin_hbca_c127_panCK_IM)
write.csv(Navin_hbca_c127_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c127_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c127_panCK_IM)

# Remove Object
rm(Navin_hbca_c127_panCK)
gc()



############################ 
# Navin_hbca_c128_panCK   #
############################ 
# Load Object
Navin_hbca_c128_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c128_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c128_panCK_All <- Navin_hbca_c128_panCK[["RNA"]]$counts
Navin_hbca_c128_panCK_All <- as.matrix(Navin_hbca_c128_panCK_All, 'sparseMatrix')
Navin_hbca_c128_panCK_All <-t(Navin_hbca_c128_panCK_All)
write.csv(Navin_hbca_c128_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c128_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c128_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c128_panCK_IM <- Navin_hbca_c128_panCK[["RNA"]]$counts
Navin_hbca_c128_panCK_IM <- as.matrix(Navin_hbca_c128_panCK_IM, 'sparseMatrix')
Navin_hbca_c128_panCK_IM <- subset(Navin_hbca_c128_panCK_IM, rownames(Navin_hbca_c128_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c128_panCK_IM <-t(Navin_hbca_c128_panCK_IM)
write.csv(Navin_hbca_c128_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c128_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c128_panCK_IM)

# Remove Object
rm(Navin_hbca_c128_panCK)
gc()



############################ 
# Navin_hbca_c129_panCK   #
############################ 
# Load Object
Navin_hbca_c129_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c129_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c129_panCK_All <- Navin_hbca_c129_panCK[["RNA"]]$counts
Navin_hbca_c129_panCK_All <- as.matrix(Navin_hbca_c129_panCK_All, 'sparseMatrix')
Navin_hbca_c129_panCK_All <-t(Navin_hbca_c129_panCK_All)
write.csv(Navin_hbca_c129_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c129_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c129_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c129_panCK_IM <- Navin_hbca_c129_panCK[["RNA"]]$counts
Navin_hbca_c129_panCK_IM <- as.matrix(Navin_hbca_c129_panCK_IM, 'sparseMatrix')
Navin_hbca_c129_panCK_IM <- subset(Navin_hbca_c129_panCK_IM, rownames(Navin_hbca_c129_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c129_panCK_IM <-t(Navin_hbca_c129_panCK_IM)
write.csv(Navin_hbca_c129_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c129_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c129_panCK_IM)

# Remove Object
rm(Navin_hbca_c129_panCK)
gc()



############################ 
# Navin_hbca_c130_panCK   #
############################ 
# Load Object
Navin_hbca_c130_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c130_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c130_panCK_All <- Navin_hbca_c130_panCK[["RNA"]]$counts
Navin_hbca_c130_panCK_All <- as.matrix(Navin_hbca_c130_panCK_All, 'sparseMatrix')
Navin_hbca_c130_panCK_All <-t(Navin_hbca_c130_panCK_All)
write.csv(Navin_hbca_c130_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c130_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c130_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c130_panCK_IM <- Navin_hbca_c130_panCK[["RNA"]]$counts
Navin_hbca_c130_panCK_IM <- as.matrix(Navin_hbca_c130_panCK_IM, 'sparseMatrix')
Navin_hbca_c130_panCK_IM <- subset(Navin_hbca_c130_panCK_IM, rownames(Navin_hbca_c130_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c130_panCK_IM <-t(Navin_hbca_c130_panCK_IM)
write.csv(Navin_hbca_c130_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c130_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c130_panCK_IM)

# Remove Object
rm(Navin_hbca_c130_panCK)
gc()



############################ 
# Navin_hbca_c131_panCK   #
############################ 
# Load Object
Navin_hbca_c131_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c131_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c131_panCK_All <- Navin_hbca_c131_panCK[["RNA"]]$counts
Navin_hbca_c131_panCK_All <- as.matrix(Navin_hbca_c131_panCK_All, 'sparseMatrix')
Navin_hbca_c131_panCK_All <-t(Navin_hbca_c131_panCK_All)
write.csv(Navin_hbca_c131_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c131_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c131_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c131_panCK_IM <- Navin_hbca_c131_panCK[["RNA"]]$counts
Navin_hbca_c131_panCK_IM <- as.matrix(Navin_hbca_c131_panCK_IM, 'sparseMatrix')
Navin_hbca_c131_panCK_IM <- subset(Navin_hbca_c131_panCK_IM, rownames(Navin_hbca_c131_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c131_panCK_IM <-t(Navin_hbca_c131_panCK_IM)
write.csv(Navin_hbca_c131_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c131_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c131_panCK_IM)

# Remove Object
rm(Navin_hbca_c131_panCK)
gc()



############################ 
# Navin_hbca_c132_panCK   #
############################ 
# Load Object
Navin_hbca_c132_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c132_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c132_panCK_All <- Navin_hbca_c132_panCK[["RNA"]]$counts
Navin_hbca_c132_panCK_All <- as.matrix(Navin_hbca_c132_panCK_All, 'sparseMatrix')
Navin_hbca_c132_panCK_All <-t(Navin_hbca_c132_panCK_All)
write.csv(Navin_hbca_c132_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c132_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c132_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c132_panCK_IM <- Navin_hbca_c132_panCK[["RNA"]]$counts
Navin_hbca_c132_panCK_IM <- as.matrix(Navin_hbca_c132_panCK_IM, 'sparseMatrix')
Navin_hbca_c132_panCK_IM <- subset(Navin_hbca_c132_panCK_IM, rownames(Navin_hbca_c132_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c132_panCK_IM <-t(Navin_hbca_c132_panCK_IM)
write.csv(Navin_hbca_c132_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c132_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c132_panCK_IM)

# Remove Object
rm(Navin_hbca_c132_panCK)
gc()





############################ 
# Navin_hbca_c133_panCK   #
############################ 
# Load Object
Navin_hbca_c133_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c133_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c133_panCK_All <- Navin_hbca_c133_panCK[["RNA"]]$counts
Navin_hbca_c133_panCK_All <- as.matrix(Navin_hbca_c133_panCK_All, 'sparseMatrix')
Navin_hbca_c133_panCK_All <-t(Navin_hbca_c133_panCK_All)
write.csv(Navin_hbca_c133_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c133_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c133_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c133_panCK_IM <- Navin_hbca_c133_panCK[["RNA"]]$counts
Navin_hbca_c133_panCK_IM <- as.matrix(Navin_hbca_c133_panCK_IM, 'sparseMatrix')
Navin_hbca_c133_panCK_IM <- subset(Navin_hbca_c133_panCK_IM, rownames(Navin_hbca_c133_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c133_panCK_IM <-t(Navin_hbca_c133_panCK_IM)
write.csv(Navin_hbca_c133_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c133_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c133_panCK_IM)

# Remove Object
rm(Navin_hbca_c133_panCK)
gc()



############################ 
# Navin_hbca_c134_panCK   #
############################ 
# Load Object
Navin_hbca_c134_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c134_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c134_panCK_All <- Navin_hbca_c134_panCK[["RNA"]]$counts
Navin_hbca_c134_panCK_All <- as.matrix(Navin_hbca_c134_panCK_All, 'sparseMatrix')
Navin_hbca_c134_panCK_All <-t(Navin_hbca_c134_panCK_All)
write.csv(Navin_hbca_c134_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c134_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c134_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c134_panCK_IM <- Navin_hbca_c134_panCK[["RNA"]]$counts
Navin_hbca_c134_panCK_IM <- as.matrix(Navin_hbca_c134_panCK_IM, 'sparseMatrix')
Navin_hbca_c134_panCK_IM <- subset(Navin_hbca_c134_panCK_IM, rownames(Navin_hbca_c134_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c134_panCK_IM <-t(Navin_hbca_c134_panCK_IM)
write.csv(Navin_hbca_c134_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c134_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c134_panCK_IM)

# Remove Object
rm(Navin_hbca_c134_panCK)
gc()



############################ 
# Navin_hbca_c135_panCK   #
############################ 
# Load Object
Navin_hbca_c135_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c135_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c135_panCK_All <- Navin_hbca_c135_panCK[["RNA"]]$counts
Navin_hbca_c135_panCK_All <- as.matrix(Navin_hbca_c135_panCK_All, 'sparseMatrix')
Navin_hbca_c135_panCK_All <-t(Navin_hbca_c135_panCK_All)
write.csv(Navin_hbca_c135_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c135_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c135_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c135_panCK_IM <- Navin_hbca_c135_panCK[["RNA"]]$counts
Navin_hbca_c135_panCK_IM <- as.matrix(Navin_hbca_c135_panCK_IM, 'sparseMatrix')
Navin_hbca_c135_panCK_IM <- subset(Navin_hbca_c135_panCK_IM, rownames(Navin_hbca_c135_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c135_panCK_IM <-t(Navin_hbca_c135_panCK_IM)
write.csv(Navin_hbca_c135_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c135_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c135_panCK_IM)

# Remove Object
rm(Navin_hbca_c135_panCK)
gc()



############################ 
# Navin_hbca_c136_panCK   #
############################ 
# Load Object
Navin_hbca_c136_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c136_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c136_panCK_All <- Navin_hbca_c136_panCK[["RNA"]]$counts
Navin_hbca_c136_panCK_All <- as.matrix(Navin_hbca_c136_panCK_All, 'sparseMatrix')
Navin_hbca_c136_panCK_All <-t(Navin_hbca_c136_panCK_All)
write.csv(Navin_hbca_c136_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c136_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c136_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c136_panCK_IM <- Navin_hbca_c136_panCK[["RNA"]]$counts
Navin_hbca_c136_panCK_IM <- as.matrix(Navin_hbca_c136_panCK_IM, 'sparseMatrix')
Navin_hbca_c136_panCK_IM <- subset(Navin_hbca_c136_panCK_IM, rownames(Navin_hbca_c136_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c136_panCK_IM <-t(Navin_hbca_c136_panCK_IM)
write.csv(Navin_hbca_c136_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c136_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c136_panCK_IM)

# Remove Object
rm(Navin_hbca_c136_panCK)
gc()



############################ 
# Navin_hbca_c137_panCK   #
############################ 
# Load Object
Navin_hbca_c137_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c137_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c137_panCK_All <- Navin_hbca_c137_panCK[["RNA"]]$counts
Navin_hbca_c137_panCK_All <- as.matrix(Navin_hbca_c137_panCK_All, 'sparseMatrix')
Navin_hbca_c137_panCK_All <-t(Navin_hbca_c137_panCK_All)
write.csv(Navin_hbca_c137_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c137_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c137_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c137_panCK_IM <- Navin_hbca_c137_panCK[["RNA"]]$counts
Navin_hbca_c137_panCK_IM <- as.matrix(Navin_hbca_c137_panCK_IM, 'sparseMatrix')
Navin_hbca_c137_panCK_IM <- subset(Navin_hbca_c137_panCK_IM, rownames(Navin_hbca_c137_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c137_panCK_IM <-t(Navin_hbca_c137_panCK_IM)
write.csv(Navin_hbca_c137_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c137_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c137_panCK_IM)

# Remove Object
rm(Navin_hbca_c137_panCK)
gc()



############################ 
# Navin_hbca_c138_panCK   #
############################ 
# Load Object
Navin_hbca_c138_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c138_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c138_panCK_All <- Navin_hbca_c138_panCK[["RNA"]]$counts
Navin_hbca_c138_panCK_All <- as.matrix(Navin_hbca_c138_panCK_All, 'sparseMatrix')
Navin_hbca_c138_panCK_All <-t(Navin_hbca_c138_panCK_All)
write.csv(Navin_hbca_c138_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c138_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c138_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c138_panCK_IM <- Navin_hbca_c138_panCK[["RNA"]]$counts
Navin_hbca_c138_panCK_IM <- as.matrix(Navin_hbca_c138_panCK_IM, 'sparseMatrix')
Navin_hbca_c138_panCK_IM <- subset(Navin_hbca_c138_panCK_IM, rownames(Navin_hbca_c138_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c138_panCK_IM <-t(Navin_hbca_c138_panCK_IM)
write.csv(Navin_hbca_c138_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c138_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c138_panCK_IM)

# Remove Object
rm(Navin_hbca_c138_panCK)
gc()



############################ 
# Navin_hbca_c139_panCK   #
############################ 
# Load Object
Navin_hbca_c139_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c139_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c139_panCK_All <- Navin_hbca_c139_panCK[["RNA"]]$counts
Navin_hbca_c139_panCK_All <- as.matrix(Navin_hbca_c139_panCK_All, 'sparseMatrix')
Navin_hbca_c139_panCK_All <-t(Navin_hbca_c139_panCK_All)
write.csv(Navin_hbca_c139_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c139_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c139_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c139_panCK_IM <- Navin_hbca_c139_panCK[["RNA"]]$counts
Navin_hbca_c139_panCK_IM <- as.matrix(Navin_hbca_c139_panCK_IM, 'sparseMatrix')
Navin_hbca_c139_panCK_IM <- subset(Navin_hbca_c139_panCK_IM, rownames(Navin_hbca_c139_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c139_panCK_IM <-t(Navin_hbca_c139_panCK_IM)
write.csv(Navin_hbca_c139_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c139_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c139_panCK_IM)

# Remove Object
rm(Navin_hbca_c139_panCK)
gc()



############################ 
# Navin_hbca_c140_panCK   #
############################ 
# Load Object
Navin_hbca_c140_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c140_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c140_panCK_All <- Navin_hbca_c140_panCK[["RNA"]]$counts
Navin_hbca_c140_panCK_All <- as.matrix(Navin_hbca_c140_panCK_All, 'sparseMatrix')
Navin_hbca_c140_panCK_All <-t(Navin_hbca_c140_panCK_All)
write.csv(Navin_hbca_c140_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c140_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c140_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c140_panCK_IM <- Navin_hbca_c140_panCK[["RNA"]]$counts
Navin_hbca_c140_panCK_IM <- as.matrix(Navin_hbca_c140_panCK_IM, 'sparseMatrix')
Navin_hbca_c140_panCK_IM <- subset(Navin_hbca_c140_panCK_IM, rownames(Navin_hbca_c140_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c140_panCK_IM <-t(Navin_hbca_c140_panCK_IM)
write.csv(Navin_hbca_c140_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c140_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c140_panCK_IM)

# Remove Object
rm(Navin_hbca_c140_panCK)
gc()



############################ 
# Navin_hbca_c141_panCK   #
############################ 
# Load Object
Navin_hbca_c141_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c141_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c141_panCK_All <- Navin_hbca_c141_panCK[["RNA"]]$counts
Navin_hbca_c141_panCK_All <- as.matrix(Navin_hbca_c141_panCK_All, 'sparseMatrix')
Navin_hbca_c141_panCK_All <-t(Navin_hbca_c141_panCK_All)
write.csv(Navin_hbca_c141_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c141_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c141_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c141_panCK_IM <- Navin_hbca_c141_panCK[["RNA"]]$counts
Navin_hbca_c141_panCK_IM <- as.matrix(Navin_hbca_c141_panCK_IM, 'sparseMatrix')
Navin_hbca_c141_panCK_IM <- subset(Navin_hbca_c141_panCK_IM, rownames(Navin_hbca_c141_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c141_panCK_IM <-t(Navin_hbca_c141_panCK_IM)
write.csv(Navin_hbca_c141_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c141_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c141_panCK_IM)

# Remove Object
rm(Navin_hbca_c141_panCK)
gc()



############################ 
# Navin_hbca_c142_panCK   #
############################ 
# Load Object
Navin_hbca_c142_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c142_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c142_panCK_All <- Navin_hbca_c142_panCK[["RNA"]]$counts
Navin_hbca_c142_panCK_All <- as.matrix(Navin_hbca_c142_panCK_All, 'sparseMatrix')
Navin_hbca_c142_panCK_All <-t(Navin_hbca_c142_panCK_All)
write.csv(Navin_hbca_c142_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c142_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c142_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c142_panCK_IM <- Navin_hbca_c142_panCK[["RNA"]]$counts
Navin_hbca_c142_panCK_IM <- as.matrix(Navin_hbca_c142_panCK_IM, 'sparseMatrix')
Navin_hbca_c142_panCK_IM <- subset(Navin_hbca_c142_panCK_IM, rownames(Navin_hbca_c142_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c142_panCK_IM <-t(Navin_hbca_c142_panCK_IM)
write.csv(Navin_hbca_c142_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c142_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c142_panCK_IM)

# Remove Object
rm(Navin_hbca_c142_panCK)
gc()



############################ 
# Navin_hbca_c143_panCK   #
############################ 
# Load Object
Navin_hbca_c143_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c143_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c143_panCK_All <- Navin_hbca_c143_panCK[["RNA"]]$counts
Navin_hbca_c143_panCK_All <- as.matrix(Navin_hbca_c143_panCK_All, 'sparseMatrix')
Navin_hbca_c143_panCK_All <-t(Navin_hbca_c143_panCK_All)
write.csv(Navin_hbca_c143_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c143_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c143_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c143_panCK_IM <- Navin_hbca_c143_panCK[["RNA"]]$counts
Navin_hbca_c143_panCK_IM <- as.matrix(Navin_hbca_c143_panCK_IM, 'sparseMatrix')
Navin_hbca_c143_panCK_IM <- subset(Navin_hbca_c143_panCK_IM, rownames(Navin_hbca_c143_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c143_panCK_IM <-t(Navin_hbca_c143_panCK_IM)
write.csv(Navin_hbca_c143_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c143_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c143_panCK_IM)

# Remove Object
rm(Navin_hbca_c143_panCK)
gc()




############################ 
# Navin_hbca_c144_panCK   #
############################ 
# Load Object
Navin_hbca_c144_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c144_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c144_panCK_All <- Navin_hbca_c144_panCK[["RNA"]]$counts
Navin_hbca_c144_panCK_All <- as.matrix(Navin_hbca_c144_panCK_All, 'sparseMatrix')
Navin_hbca_c144_panCK_All <-t(Navin_hbca_c144_panCK_All)
write.csv(Navin_hbca_c144_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c144_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c144_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c144_panCK_IM <- Navin_hbca_c144_panCK[["RNA"]]$counts
Navin_hbca_c144_panCK_IM <- as.matrix(Navin_hbca_c144_panCK_IM, 'sparseMatrix')
Navin_hbca_c144_panCK_IM <- subset(Navin_hbca_c144_panCK_IM, rownames(Navin_hbca_c144_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c144_panCK_IM <-t(Navin_hbca_c144_panCK_IM)
write.csv(Navin_hbca_c144_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c144_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c144_panCK_IM)

# Remove Object
rm(Navin_hbca_c144_panCK)
gc()




############################ 
# Navin_hbca_c145_panCK   #
############################ 
# Load Object
Navin_hbca_c145_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c145_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c145_panCK_All <- Navin_hbca_c145_panCK[["RNA"]]$counts
Navin_hbca_c145_panCK_All <- as.matrix(Navin_hbca_c145_panCK_All, 'sparseMatrix')
Navin_hbca_c145_panCK_All <-t(Navin_hbca_c145_panCK_All)
write.csv(Navin_hbca_c145_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c145_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c145_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c145_panCK_IM <- Navin_hbca_c145_panCK[["RNA"]]$counts
Navin_hbca_c145_panCK_IM <- as.matrix(Navin_hbca_c145_panCK_IM, 'sparseMatrix')
Navin_hbca_c145_panCK_IM <- subset(Navin_hbca_c145_panCK_IM, rownames(Navin_hbca_c145_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c145_panCK_IM <-t(Navin_hbca_c145_panCK_IM)
write.csv(Navin_hbca_c145_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c145_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c145_panCK_IM)

# Remove Object
rm(Navin_hbca_c145_panCK)
gc()




############################ 
# Navin_hbca_c146_panCK   #
############################ 
# Load Object
Navin_hbca_c146_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c146_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c146_panCK_All <- Navin_hbca_c146_panCK[["RNA"]]$counts
Navin_hbca_c146_panCK_All <- as.matrix(Navin_hbca_c146_panCK_All, 'sparseMatrix')
Navin_hbca_c146_panCK_All <-t(Navin_hbca_c146_panCK_All)
write.csv(Navin_hbca_c146_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c146_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c146_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c146_panCK_IM <- Navin_hbca_c146_panCK[["RNA"]]$counts
Navin_hbca_c146_panCK_IM <- as.matrix(Navin_hbca_c146_panCK_IM, 'sparseMatrix')
Navin_hbca_c146_panCK_IM <- subset(Navin_hbca_c146_panCK_IM, rownames(Navin_hbca_c146_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c146_panCK_IM <-t(Navin_hbca_c146_panCK_IM)
write.csv(Navin_hbca_c146_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c146_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c146_panCK_IM)

# Remove Object
rm(Navin_hbca_c146_panCK)
gc()




############################ 
# Navin_hbca_c147_panCK   #
############################ 
# Load Object
Navin_hbca_c147_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c147_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c147_panCK_All <- Navin_hbca_c147_panCK[["RNA"]]$counts
Navin_hbca_c147_panCK_All <- as.matrix(Navin_hbca_c147_panCK_All, 'sparseMatrix')
Navin_hbca_c147_panCK_All <-t(Navin_hbca_c147_panCK_All)
write.csv(Navin_hbca_c147_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c147_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c147_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c147_panCK_IM <- Navin_hbca_c147_panCK[["RNA"]]$counts
Navin_hbca_c147_panCK_IM <- as.matrix(Navin_hbca_c147_panCK_IM, 'sparseMatrix')
Navin_hbca_c147_panCK_IM <- subset(Navin_hbca_c147_panCK_IM, rownames(Navin_hbca_c147_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c147_panCK_IM <-t(Navin_hbca_c147_panCK_IM)
write.csv(Navin_hbca_c147_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c147_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c147_panCK_IM)

# Remove Object
rm(Navin_hbca_c147_panCK)
gc()




############################ 
# Navin_hbca_c148_panCK   #
############################ 
# Load Object
Navin_hbca_c148_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c148_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c148_panCK_All <- Navin_hbca_c148_panCK[["RNA"]]$counts
Navin_hbca_c148_panCK_All <- as.matrix(Navin_hbca_c148_panCK_All, 'sparseMatrix')
Navin_hbca_c148_panCK_All <-t(Navin_hbca_c148_panCK_All)
write.csv(Navin_hbca_c148_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c148_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c148_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c148_panCK_IM <- Navin_hbca_c148_panCK[["RNA"]]$counts
Navin_hbca_c148_panCK_IM <- as.matrix(Navin_hbca_c148_panCK_IM, 'sparseMatrix')
Navin_hbca_c148_panCK_IM <- subset(Navin_hbca_c148_panCK_IM, rownames(Navin_hbca_c148_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c148_panCK_IM <-t(Navin_hbca_c148_panCK_IM)
write.csv(Navin_hbca_c148_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c148_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c148_panCK_IM)

# Remove Object
rm(Navin_hbca_c148_panCK)
gc()




############################ 
# Navin_hbca_c149_panCK   #
############################ 
# Load Object
Navin_hbca_c149_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c149_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c149_panCK_All <- Navin_hbca_c149_panCK[["RNA"]]$counts
Navin_hbca_c149_panCK_All <- as.matrix(Navin_hbca_c149_panCK_All, 'sparseMatrix')
Navin_hbca_c149_panCK_All <-t(Navin_hbca_c149_panCK_All)
write.csv(Navin_hbca_c149_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c149_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c149_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c149_panCK_IM <- Navin_hbca_c149_panCK[["RNA"]]$counts
Navin_hbca_c149_panCK_IM <- as.matrix(Navin_hbca_c149_panCK_IM, 'sparseMatrix')
Navin_hbca_c149_panCK_IM <- subset(Navin_hbca_c149_panCK_IM, rownames(Navin_hbca_c149_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c149_panCK_IM <-t(Navin_hbca_c149_panCK_IM)
write.csv(Navin_hbca_c149_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c149_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c149_panCK_IM)

# Remove Object
rm(Navin_hbca_c149_panCK)
gc()




############################ 
# Navin_hbca_c150_panCK   #
############################ 
# Load Object
Navin_hbca_c150_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c150_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c150_panCK_All <- Navin_hbca_c150_panCK[["RNA"]]$counts
Navin_hbca_c150_panCK_All <- as.matrix(Navin_hbca_c150_panCK_All, 'sparseMatrix')
Navin_hbca_c150_panCK_All <-t(Navin_hbca_c150_panCK_All)
write.csv(Navin_hbca_c150_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c150_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c150_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c150_panCK_IM <- Navin_hbca_c150_panCK[["RNA"]]$counts
Navin_hbca_c150_panCK_IM <- as.matrix(Navin_hbca_c150_panCK_IM, 'sparseMatrix')
Navin_hbca_c150_panCK_IM <- subset(Navin_hbca_c150_panCK_IM, rownames(Navin_hbca_c150_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c150_panCK_IM <-t(Navin_hbca_c150_panCK_IM)
write.csv(Navin_hbca_c150_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c150_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c150_panCK_IM)

# Remove Object
rm(Navin_hbca_c150_panCK)
gc()




############################ 
# Navin_hbca_c151_panCK   #
############################ 
# Load Object
Navin_hbca_c151_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c151_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c151_panCK_All <- Navin_hbca_c151_panCK[["RNA"]]$counts
Navin_hbca_c151_panCK_All <- as.matrix(Navin_hbca_c151_panCK_All, 'sparseMatrix')
Navin_hbca_c151_panCK_All <-t(Navin_hbca_c151_panCK_All)
write.csv(Navin_hbca_c151_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c151_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c151_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c151_panCK_IM <- Navin_hbca_c151_panCK[["RNA"]]$counts
Navin_hbca_c151_panCK_IM <- as.matrix(Navin_hbca_c151_panCK_IM, 'sparseMatrix')
Navin_hbca_c151_panCK_IM <- subset(Navin_hbca_c151_panCK_IM, rownames(Navin_hbca_c151_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c151_panCK_IM <-t(Navin_hbca_c151_panCK_IM)
write.csv(Navin_hbca_c151_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c151_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c151_panCK_IM)

# Remove Object
rm(Navin_hbca_c151_panCK)
gc()




############################ 
# Navin_hbca_c152_panCK   #
############################ 
# Load Object
Navin_hbca_c152_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c152_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c152_panCK_All <- Navin_hbca_c152_panCK[["RNA"]]$counts
Navin_hbca_c152_panCK_All <- as.matrix(Navin_hbca_c152_panCK_All, 'sparseMatrix')
Navin_hbca_c152_panCK_All <-t(Navin_hbca_c152_panCK_All)
write.csv(Navin_hbca_c152_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c152_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c152_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c152_panCK_IM <- Navin_hbca_c152_panCK[["RNA"]]$counts
Navin_hbca_c152_panCK_IM <- as.matrix(Navin_hbca_c152_panCK_IM, 'sparseMatrix')
Navin_hbca_c152_panCK_IM <- subset(Navin_hbca_c152_panCK_IM, rownames(Navin_hbca_c152_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c152_panCK_IM <-t(Navin_hbca_c152_panCK_IM)
write.csv(Navin_hbca_c152_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c152_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c152_panCK_IM)

# Remove Object
rm(Navin_hbca_c152_panCK)
gc()




############################ 
# Navin_hbca_c153_panCK   #
############################ 
# Load Object
Navin_hbca_c153_panCK  <- readRDS(file = "/R/R_Navin/Navin_RDS/RDS_panCK/Navin_hbca_c153_panCK .rds")

# Extract Entire Expression Matrix
Navin_hbca_c153_panCK_All <- Navin_hbca_c153_panCK[["RNA"]]$counts
Navin_hbca_c153_panCK_All <- as.matrix(Navin_hbca_c153_panCK_All, 'sparseMatrix')
Navin_hbca_c153_panCK_All <-t(Navin_hbca_c153_panCK_All)
write.csv(Navin_hbca_c153_panCK_All, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c153_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Navin_hbca_c153_panCK_All)

# Extract IM Expression Matrix
Navin_hbca_c153_panCK_IM <- Navin_hbca_c153_panCK[["RNA"]]$counts
Navin_hbca_c153_panCK_IM <- as.matrix(Navin_hbca_c153_panCK_IM, 'sparseMatrix')
Navin_hbca_c153_panCK_IM <- subset(Navin_hbca_c153_panCK_IM, rownames(Navin_hbca_c153_panCK_IM) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Navin_hbca_c153_panCK_IM <-t(Navin_hbca_c153_panCK_IM)
write.csv(Navin_hbca_c153_panCK_IM, file = "/R/R_Navin/Navin_Expression_Matrices/Navin_Expression_Matrices_Output/Navin_hbca_c153_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Navin_hbca_c153_panCK_IM)

# Remove Object
rm(Navin_hbca_c153_panCK)
gc()




######################################################
################## Visvader Dataset ##################
######################################################

#############
# All Genes #
#############
############################
# Visvader_0019_NORM_panCK #
############################
# Load Data
Visvader_0019_NORM_panCK <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0019_NORM_panCK.rds")
# Extract Expression Matrix
Visvader_0019_NORM_panCK <- Visvader_0019_NORM_panCK[["RNA"]]$counts
Visvader_0019_NORM_panCK <- as.matrix(Visvader_0019_NORM_panCK, 'sparseMatrix')
Visvader_0019_NORM_panCK <-t(Visvader_0019_NORM_panCK)
write.csv(Visvader_0019_NORM_panCK, file = "/R/R_Visvader/Visvader_Expression_Matrices/Visvader_Expression_Matrices_Output/Visvader_0019_NORM_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Visvader_0019_NORM_panCK)
gc()



############################
# Visvader_0021_NORM_panCK #
############################
# Load Data
Visvader_0021_NORM_panCK <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0021_NORM_panCK.rds")
# Extract Expression Matrix
Visvader_0021_NORM_panCK <- Visvader_0021_NORM_panCK[["RNA"]]$counts
Visvader_0021_NORM_panCK <- as.matrix(Visvader_0021_NORM_panCK, 'sparseMatrix')
Visvader_0021_NORM_panCK <-t(Visvader_0021_NORM_panCK)
write.csv(Visvader_0021_NORM_panCK, file = "/R/R_Visvader/Visvader_Expression_Matrices/Visvader_Expression_Matrices_Output/Visvader_0021_NORM_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Visvader_0021_NORM_panCK)
gc()



############################
# Visvader_0064_NORM_panCK #
############################
# Load Data
Visvader_0064_NORM_panCK <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0064_NORM_panCK.rds")
# Extract Expression Matrix
Visvader_0064_NORM_panCK <- Visvader_0064_NORM_panCK[["RNA"]]$counts
Visvader_0064_NORM_panCK <- as.matrix(Visvader_0064_NORM_panCK, 'sparseMatrix')
Visvader_0064_NORM_panCK <-t(Visvader_0064_NORM_panCK)
write.csv(Visvader_0064_NORM_panCK, file = "/R/R_Visvader/Visvader_Expression_Matrices/Visvader_Expression_Matrices_Output/Visvader_0064_NORM_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Visvader_0064_NORM_panCK)
gc()



############################
# Visvader_0092_NORM_panCK #
############################
# Load Data
Visvader_0092_NORM_panCK <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0092_NORM_panCK.rds")
# Extract Expression Matrix
Visvader_0092_NORM_panCK <- Visvader_0092_NORM_panCK[["RNA"]]$counts
Visvader_0092_NORM_panCK <- as.matrix(Visvader_0092_NORM_panCK, 'sparseMatrix')
Visvader_0092_NORM_panCK <-t(Visvader_0092_NORM_panCK)
write.csv(Visvader_0092_NORM_panCK, file = "/R/R_Visvader/Visvader_Expression_Matrices/Visvader_Expression_Matrices_Output/Visvader_0092_NORM_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Visvader_0092_NORM_panCK)
gc()



############################
# Visvader_0093_NORM_panCK #
############################
# Load Data
Visvader_0093_NORM_panCK <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0093_NORM_panCK.rds")
# Extract Expression Matrix
Visvader_0093_NORM_panCK <- Visvader_0093_NORM_panCK[["RNA"]]$counts
Visvader_0093_NORM_panCK <- as.matrix(Visvader_0093_NORM_panCK, 'sparseMatrix')
Visvader_0093_NORM_panCK <-t(Visvader_0093_NORM_panCK)
write.csv(Visvader_0093_NORM_panCK, file = "/R/R_Visvader/Visvader_Expression_Matrices/Visvader_Expression_Matrices_Output/Visvader_0093_NORM_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Visvader_0093_NORM_panCK)
gc()


############################
# Visvader_0123_NORM_panCK #
############################
# Load Data
Visvader_0123_NORM_panCK <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0123_NORM_panCK.rds")
# Extract Expression Matrix
Visvader_0123_NORM_panCK <- Visvader_0123_NORM_panCK[["RNA"]]$counts
Visvader_0123_NORM_panCK <- as.matrix(Visvader_0123_NORM_panCK, 'sparseMatrix')
Visvader_0123_NORM_panCK <-t(Visvader_0123_NORM_panCK)
write.csv(Visvader_0123_NORM_panCK, file = "/R/R_Visvader/Visvader_Expression_Matrices/Visvader_Expression_Matrices_Output/Visvader_0123_NORM_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Visvader_0123_NORM_panCK)
gc()



############################
# Visvader_0169_NORM_panCK #
############################
# Load Data
Visvader_0169_NORM_panCK <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0169_NORM_panCK.rds")
# Extract Expression Matrix
Visvader_0169_NORM_panCK <- Visvader_0169_NORM_panCK[["RNA"]]$counts
Visvader_0169_NORM_panCK <- as.matrix(Visvader_0169_NORM_panCK, 'sparseMatrix')
Visvader_0169_NORM_panCK <-t(Visvader_0169_NORM_panCK)
write.csv(Visvader_0169_NORM_panCK, file = "/R/R_Visvader/Visvader_Expression_Matrices/Visvader_Expression_Matrices_Output/Visvader_0169_NORM_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Visvader_0169_NORM_panCK)
gc()


############################
# Visvader_0230_NORM_panCK #
############################
# Load Data
Visvader_0230_NORM_panCK <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0230_NORM_panCK.rds")
# Extract Expression Matrix
Visvader_0230_NORM_panCK <- Visvader_0230_NORM_panCK[["RNA"]]$counts
Visvader_0230_NORM_panCK <- as.matrix(Visvader_0230_NORM_panCK, 'sparseMatrix')
Visvader_0230_NORM_panCK <-t(Visvader_0230_NORM_panCK)
write.csv(Visvader_0230_NORM_panCK, file = "/R/R_Visvader/Visvader_Expression_Matrices/Visvader_Expression_Matrices_Output/Visvader_0230_NORM_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Visvader_0230_NORM_panCK)
gc()


############################
# Visvader_0233_NORM_panCK #
############################
# Load Data
Visvader_0233_NORM_panCK <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0233_NORM_panCK.rds")
# Extract Expression Matrix
Visvader_0233_NORM_panCK <- Visvader_0233_NORM_panCK[["RNA"]]$counts
Visvader_0233_NORM_panCK <- as.matrix(Visvader_0233_NORM_panCK, 'sparseMatrix')
Visvader_0233_NORM_panCK <-t(Visvader_0233_NORM_panCK)
write.csv(Visvader_0233_NORM_panCK, file = "/R/R_Visvader/Visvader_Expression_Matrices/Visvader_Expression_Matrices_Output/Visvader_0233_NORM_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Visvader_0233_NORM_panCK)
gc()


############################
# Visvader_0275_NORM_panCK #
############################
# Load Data
Visvader_0275_NORM_panCK <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0275_NORM_panCK.rds")
# Extract Expression Matrix
Visvader_0275_NORM_panCK <- Visvader_0275_NORM_panCK[["RNA"]]$counts
Visvader_0275_NORM_panCK <- as.matrix(Visvader_0275_NORM_panCK, 'sparseMatrix')
Visvader_0275_NORM_panCK <-t(Visvader_0275_NORM_panCK)
write.csv(Visvader_0275_NORM_panCK, file = "/R/R_Visvader/Visvader_Expression_Matrices/Visvader_Expression_Matrices_Output/Visvader_0275_NORM_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Visvader_0275_NORM_panCK)
gc()

############################
# Visvader_0288_NORM_panCK #
############################
# Load Data
Visvader_0288_NORM_panCK <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0288_NORM_panCK.rds")
# Extract Expression Matrix
Visvader_0288_NORM_panCK <- Visvader_0288_NORM_panCK[["RNA"]]$counts
Visvader_0288_NORM_panCK <- as.matrix(Visvader_0288_NORM_panCK, 'sparseMatrix')
Visvader_0288_NORM_panCK <-t(Visvader_0288_NORM_panCK)
write.csv(Visvader_0288_NORM_panCK, file = "/R/R_Visvader/Visvader_Expression_Matrices/Visvader_Expression_Matrices_Output/Visvader_0288_NORM_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Visvader_0288_NORM_panCK)
gc()


############################
# Visvader_0342_NORM_panCK #
############################
# Load Data
Visvader_0342_NORM_panCK <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0342_NORM_panCK.rds")
# Extract Expression Matrix
Visvader_0342_NORM_panCK <- Visvader_0342_NORM_panCK[["RNA"]]$counts
Visvader_0342_NORM_panCK <- as.matrix(Visvader_0342_NORM_panCK, 'sparseMatrix')
Visvader_0342_NORM_panCK <-t(Visvader_0342_NORM_panCK)
write.csv(Visvader_0342_NORM_panCK, file = "/R/R_Visvader/Visvader_Expression_Matrices/Visvader_Expression_Matrices_Output/Visvader_0342_NORM_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Visvader_0342_NORM_panCK)
gc()


############################
# Visvader_0372_NORM_panCK #
############################
# Load Data
Visvader_0372_NORM_panCK <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0372_NORM_panCK.rds")
# Extract Expression Matrix
Visvader_0372_NORM_panCK <- Visvader_0372_NORM_panCK[["RNA"]]$counts
Visvader_0372_NORM_panCK <- as.matrix(Visvader_0372_NORM_panCK, 'sparseMatrix')
Visvader_0372_NORM_panCK <-t(Visvader_0372_NORM_panCK)
write.csv(Visvader_0372_NORM_panCK, file = "/R/R_Visvader/Visvader_Expression_Matrices/Visvader_Expression_Matrices_Output/Visvader_0372_NORM_panCK_Immune_Mimicry_Cell_Expression_Matrix_All_Genes.csv")
rm(Visvader_0372_NORM_panCK)
gc()



##############
# IMBC Genes #
##############
############################
# Visvader_0019_NORM_panCK #
############################
# Load Data
Visvader_0019_NORM_panCK <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0019_NORM_panCK.rds")
# Extract Expression Matrix
Visvader_0019_NORM_panCK <- Visvader_0019_NORM_panCK[["RNA"]]$counts
Visvader_0019_NORM_panCK <- as.matrix(Visvader_0019_NORM_panCK, 'sparseMatrix')
Visvader_0019_NORM_panCK <- subset(Visvader_0019_NORM_panCK, rownames(Visvader_0019_NORM_panCK) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Visvader_0019_NORM_panCK <-t(Visvader_0019_NORM_panCK)
write.csv(Visvader_0019_NORM_panCK, file = "/R/R_Visvader/Visvader_Expression_Matrices/Visvader_Expression_Matrices_Output/Visvader_0019_NORM_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Visvader_0019_NORM_panCK)
gc()



############################
# Visvader_0021_NORM_panCK #
############################
# Load Data
Visvader_0021_NORM_panCK <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0021_NORM_panCK.rds")
# Extract Expression Matrix
Visvader_0021_NORM_panCK <- Visvader_0021_NORM_panCK[["RNA"]]$counts
Visvader_0021_NORM_panCK <- as.matrix(Visvader_0021_NORM_panCK, 'sparseMatrix')
Visvader_0021_NORM_panCK <- subset(Visvader_0021_NORM_panCK, rownames(Visvader_0021_NORM_panCK) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Visvader_0021_NORM_panCK <-t(Visvader_0021_NORM_panCK)
write.csv(Visvader_0021_NORM_panCK, file = "/R/R_Visvader/Visvader_Expression_Matrices/Visvader_Expression_Matrices_Output/Visvader_0021_NORM_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Visvader_0021_NORM_panCK)
gc()



############################
# Visvader_0064_NORM_panCK #
############################
# Load Data
Visvader_0064_NORM_panCK <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0064_NORM_panCK.rds")
# Extract Expression Matrix
Visvader_0064_NORM_panCK <- Visvader_0064_NORM_panCK[["RNA"]]$counts
Visvader_0064_NORM_panCK <- as.matrix(Visvader_0064_NORM_panCK, 'sparseMatrix')
Visvader_0064_NORM_panCK <- subset(Visvader_0064_NORM_panCK, rownames(Visvader_0064_NORM_panCK) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Visvader_0064_NORM_panCK <-t(Visvader_0064_NORM_panCK)
write.csv(Visvader_0064_NORM_panCK, file = "/R/R_Visvader/Visvader_Expression_Matrices/Visvader_Expression_Matrices_Output/Visvader_0064_NORM_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Visvader_0064_NORM_panCK)
gc()



############################
# Visvader_0092_NORM_panCK #
############################
# Load Data
Visvader_0092_NORM_panCK <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0092_NORM_panCK.rds")
# Extract Expression Matrix
Visvader_0092_NORM_panCK <- Visvader_0092_NORM_panCK[["RNA"]]$counts
Visvader_0092_NORM_panCK <- as.matrix(Visvader_0092_NORM_panCK, 'sparseMatrix')
Visvader_0092_NORM_panCK <- subset(Visvader_0092_NORM_panCK, rownames(Visvader_0092_NORM_panCK) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Visvader_0092_NORM_panCK <-t(Visvader_0092_NORM_panCK)
write.csv(Visvader_0092_NORM_panCK, file = "/R/R_Visvader/Visvader_Expression_Matrices/Visvader_Expression_Matrices_Output/Visvader_0092_NORM_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Visvader_0092_NORM_panCK)
gc()



############################
# Visvader_0093_NORM_panCK #
############################
# Load Data
Visvader_0093_NORM_panCK <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0093_NORM_panCK.rds")
# Extract Expression Matrix
Visvader_0093_NORM_panCK <- Visvader_0093_NORM_panCK[["RNA"]]$counts
Visvader_0093_NORM_panCK <- as.matrix(Visvader_0093_NORM_panCK, 'sparseMatrix')
Visvader_0093_NORM_panCK <- subset(Visvader_0093_NORM_panCK, rownames(Visvader_0093_NORM_panCK) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Visvader_0093_NORM_panCK <-t(Visvader_0093_NORM_panCK)
write.csv(Visvader_0093_NORM_panCK, file = "/R/R_Visvader/Visvader_Expression_Matrices/Visvader_Expression_Matrices_Output/Visvader_0093_NORM_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Visvader_0093_NORM_panCK)
gc()


############################
# Visvader_0123_NORM_panCK #
############################
# Load Data
Visvader_0123_NORM_panCK <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0123_NORM_panCK.rds")
# Extract Expression Matrix
Visvader_0123_NORM_panCK <- Visvader_0123_NORM_panCK[["RNA"]]$counts
Visvader_0123_NORM_panCK <- as.matrix(Visvader_0123_NORM_panCK, 'sparseMatrix')
Visvader_0123_NORM_panCK <- subset(Visvader_0123_NORM_panCK, rownames(Visvader_0123_NORM_panCK) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Visvader_0123_NORM_panCK <-t(Visvader_0123_NORM_panCK)
write.csv(Visvader_0123_NORM_panCK, file = "/R/R_Visvader/Visvader_Expression_Matrices/Visvader_Expression_Matrices_Output/Visvader_0123_NORM_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Visvader_0123_NORM_panCK)
gc()



############################
# Visvader_0169_NORM_panCK #
############################
# Load Data
Visvader_0169_NORM_panCK <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0169_NORM_panCK.rds")
# Extract Expression Matrix
Visvader_0169_NORM_panCK <- Visvader_0169_NORM_panCK[["RNA"]]$counts
Visvader_0169_NORM_panCK <- as.matrix(Visvader_0169_NORM_panCK, 'sparseMatrix')
Visvader_0169_NORM_panCK <- subset(Visvader_0169_NORM_panCK, rownames(Visvader_0169_NORM_panCK) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Visvader_0169_NORM_panCK <-t(Visvader_0169_NORM_panCK)
write.csv(Visvader_0169_NORM_panCK, file = "/R/R_Visvader/Visvader_Expression_Matrices/Visvader_Expression_Matrices_Output/Visvader_0169_NORM_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Visvader_0169_NORM_panCK)
gc()


############################
# Visvader_0230_NORM_panCK #
############################
# Load Data
Visvader_0230_NORM_panCK <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0230_NORM_panCK.rds")
# Extract Expression Matrix
Visvader_0230_NORM_panCK <- Visvader_0230_NORM_panCK[["RNA"]]$counts
Visvader_0230_NORM_panCK <- as.matrix(Visvader_0230_NORM_panCK, 'sparseMatrix')
Visvader_0230_NORM_panCK <- subset(Visvader_0230_NORM_panCK, rownames(Visvader_0230_NORM_panCK) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Visvader_0230_NORM_panCK <-t(Visvader_0230_NORM_panCK)
write.csv(Visvader_0230_NORM_panCK, file = "/R/R_Visvader/Visvader_Expression_Matrices/Visvader_Expression_Matrices_Output/Visvader_0230_NORM_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Visvader_0230_NORM_panCK)
gc()


############################
# Visvader_0233_NORM_panCK #
############################
# Load Data
Visvader_0233_NORM_panCK <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0233_NORM_panCK.rds")
# Extract Expression Matrix
Visvader_0233_NORM_panCK <- Visvader_0233_NORM_panCK[["RNA"]]$counts
Visvader_0233_NORM_panCK <- as.matrix(Visvader_0233_NORM_panCK, 'sparseMatrix')
Visvader_0233_NORM_panCK <- subset(Visvader_0233_NORM_panCK, rownames(Visvader_0233_NORM_panCK) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Visvader_0233_NORM_panCK <-t(Visvader_0233_NORM_panCK)
write.csv(Visvader_0233_NORM_panCK, file = "/R/R_Visvader/Visvader_Expression_Matrices/Visvader_Expression_Matrices_Output/Visvader_0233_NORM_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Visvader_0233_NORM_panCK)
gc()


############################
# Visvader_0275_NORM_panCK #
############################
# Load Data
Visvader_0275_NORM_panCK <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0275_NORM_panCK.rds")
# Extract Expression Matrix
Visvader_0275_NORM_panCK <- Visvader_0275_NORM_panCK[["RNA"]]$counts
Visvader_0275_NORM_panCK <- as.matrix(Visvader_0275_NORM_panCK, 'sparseMatrix')
Visvader_0275_NORM_panCK <- subset(Visvader_0275_NORM_panCK, rownames(Visvader_0275_NORM_panCK) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Visvader_0275_NORM_panCK <-t(Visvader_0275_NORM_panCK)
write.csv(Visvader_0275_NORM_panCK, file = "/R/R_Visvader/Visvader_Expression_Matrices/Visvader_Expression_Matrices_Output/Visvader_0275_NORM_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Visvader_0275_NORM_panCK)
gc()

############################
# Visvader_0288_NORM_panCK #
############################
# Load Data
Visvader_0288_NORM_panCK <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0288_NORM_panCK.rds")
# Extract Expression Matrix
Visvader_0288_NORM_panCK <- Visvader_0288_NORM_panCK[["RNA"]]$counts
Visvader_0288_NORM_panCK <- as.matrix(Visvader_0288_NORM_panCK, 'sparseMatrix')
Visvader_0288_NORM_panCK <- subset(Visvader_0288_NORM_panCK, rownames(Visvader_0288_NORM_panCK) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Visvader_0288_NORM_panCK <-t(Visvader_0288_NORM_panCK)
write.csv(Visvader_0288_NORM_panCK, file = "/R/R_Visvader/Visvader_Expression_Matrices/Visvader_Expression_Matrices_Output/Visvader_0288_NORM_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Visvader_0288_NORM_panCK)
gc()


############################
# Visvader_0342_NORM_panCK #
############################
# Load Data
Visvader_0342_NORM_panCK <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0342_NORM_panCK.rds")
# Extract Expression Matrix
Visvader_0342_NORM_panCK <- Visvader_0342_NORM_panCK[["RNA"]]$counts
Visvader_0342_NORM_panCK <- as.matrix(Visvader_0342_NORM_panCK, 'sparseMatrix')
Visvader_0342_NORM_panCK <- subset(Visvader_0342_NORM_panCK, rownames(Visvader_0342_NORM_panCK) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Visvader_0342_NORM_panCK <-t(Visvader_0342_NORM_panCK)
write.csv(Visvader_0342_NORM_panCK, file = "/R/R_Visvader/Visvader_Expression_Matrices/Visvader_Expression_Matrices_Output/Visvader_0342_NORM_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Visvader_0342_NORM_panCK)
gc()


############################
# Visvader_0372_NORM_panCK #
############################
# Load Data
Visvader_0372_NORM_panCK <- readRDS(file = "/R/R_Visvader/Visvader_RDS/RDS_panCK/Visvader_0372_NORM_panCK.rds")
# Extract Expression Matrix
Visvader_0372_NORM_panCK <- Visvader_0372_NORM_panCK[["RNA"]]$counts
Visvader_0372_NORM_panCK <- as.matrix(Visvader_0372_NORM_panCK, 'sparseMatrix')
Visvader_0372_NORM_panCK <- subset(Visvader_0372_NORM_panCK, rownames(Visvader_0372_NORM_panCK) %in% c("FCGR3A", "CD2", "CD3E", "CD7", "CD3G", "CD3D", "IL7R", "FCGR2A", "CD68", "CD52", "CD69", "CD14", "PTPRC", "TNFSF13B", "CD37", "ITGB2", "CXCR4", "CD74", "CLEC7A", "CD53", "CD83", "PLAUR", "CD44"))
Visvader_0372_NORM_panCK <-t(Visvader_0372_NORM_panCK)
write.csv(Visvader_0372_NORM_panCK, file = "/R/R_Visvader/Visvader_Expression_Matrices/Visvader_Expression_Matrices_Output/Visvader_0372_NORM_panCK_Immune_Mimicry_Cell_Expression_Matrix.csv")
rm(Visvader_0372_NORM_panCK)
gc()










