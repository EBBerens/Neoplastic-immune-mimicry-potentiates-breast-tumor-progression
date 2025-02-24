# **Immune Mimicry scRNA-seq R Analysis Setup** #

### **Folder Structure** ###

Running the immune mimicry single-cell RNA-sequencing R analysis as described will require the folder structure shown below. The folder structure depicted for the R_Visvader breast tumor dataset mirrors organization for the R_Swarbrick and R_PANNTHR datasets, whereas the R_Navin, R_All_Breast_Tumors, and R_PMBC do not use every folder as they are peripheral to the main analysis. Input data were downloaded from the referenced repository and organized into their corresponding input folder; note that here breast tumors are organized by disease subtype and reduction mammoplasties go into the NORM folder.  

R \
... R_Visvader  \
... Visvader_Doublet_Tables  \
... Visvader_Expression_Matrices  \
... ... Visvader_Expression_Matrices_Input  \
... ... Visvader_Expression_Matrices_Output  \
... Visvader_fgsea  \
... ... Visvader_fgsea_MSigDB  \
... ... Visvader_fgsea_Output  \
... Visvader_Gene_Lists  \
... Visvader_inferCNV  \
... ... Visvader_inferCNV_Input  \
... ... ... Visvader_inferCNV_Input_Count_Matrix  \
... ... ... Visvader_inferCNV_Input_RDS  \
... ... Visvader_inferCNV_Output  \
... Visvader_Input  \
... ... ER  \
... ... HER2  \
... ... NORM  \
... ... PR  \
... ... TNBC  \
... Visvader_Output  \
... Visvader_RDS  \
... ... RDS_Annotated  \
... ... RDS_inferCNV  \
... ... RDS_Leukocytes  \
... ... RDS_Merged  \
... ... RDS_panCK  \
... ... RDS_Total  \
... ... RDS_Total_Singlets  \
... Visvader_Scripts \
\
R \
... R_Visvader  \
... R_Swarbrick  \
... R_PANNTHR  \
... R_Navin  \
... R_All_Breast_Tumors  \
... R_PBMC

### **After the folders are created, the following files need to be downloaded:** ###

**1.**	Download the following file from the Broad Institute and place into the “[Dataset]_inferCNV_Input” folder: hg38_gencode_v27

* This is a gene/chromosome positions file for inferCNV.

* https://data.broadinstitute.org/Trinity/CTAT/cnv/

**2.**	Download the following file from MSigDB courtesy of the Broad Institute and place into the “[Dataset]_ fgsea_MSigDB” folder: c5.go.bp.v2022.1.Hs.symbols.gmt

* These are Gene Ontology Biological Processes used for pathway analysis.

* https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2022.1.Hs/

**3.**	Download the following file from Uniprot place into the “[Dataset]_Gene_List” folder: cdlist.txt

* This file contains gene names for known CD markers.

* Note: cdlist.txt was converted to Uniprot_CD_Gene_List_2024.csv

* https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/docs/cdlist.txt

### **Recommended R Script Order** ###

If reproducing the entire analysis, I recommend executing the R scripts in the order shown below. Keys points are that: (1) the reduction mammoplasty reference must be generated for inferCNV to run properly on breast tumor datasets and (2) analyses with all breast tumors cannot proceed until the main immune mimicry analysis has been completed. 

##### **Preprocess Breast Tumors and Extract Epithelial Cells** #####
<pre>
01_Pal_et_al_Breast_Tumor_Preprocessing.R
02_Pal_et_al_Breast_Tumor_Epithelial_Cell_Extraction.R
01_Wu_et_al_Breast_Tumor_Preprocessing.R
02_Wu_et_al_Breast_Tumor_Epithelial_Cell_Extraction.R
01_PANNTHR_Breast_Tumor_Preprocessing.R
02_PANNTHR_Breast_Tumor_Epithelial_Cell_Extraction.R
</pre>

##### **Preprocess Reduction Mammoplasties, Extract Epithelial Cells, and Create Epithelial Reference** #####
<pre>
01_Pal_et_al_Reduction_Mammoplasty_Preprocessing.R
02_Pal_et_al_Reduction_Mammoplasty_Epithelial_Cell_Extraction.R
01_Kumar_et_al_Reduction_Mammoplasty_Preprocessing.R
02_Kumar_et_al_Reduction_Mammoplasty_Epithelial_Cell_Extraction.R
03_Pal_et_al_and_Kumar_et_al_Reduction_Mammoplasty_Reference_Creation.R
</pre>

##### **Extract and Merge Neoplastic Cells from Breast Tumors** #####
<pre>
03_Pal_et_al_Breast_Tumor_Neoplastic_Cell_Extraction.R
03_Wu_et_al_Breast_Tumor_Neoplastic_Cell_Extraction.R
03_PANNTHR_Breast_Tumor_Neoplastic_Cell_Extraction.R
</pre>

##### **Conduct Main Immune Mimicry Analysis on Breat Tumors** #####
<pre>
04_Pal_et_al_Breast_Tumor_Main_Immune_Mimicry_Analysis.R
04_Wu_et_al_Breast_Tumor_Main_Immune_Mimicry_Analysis.R
04_PANNTHR_Breast_Tumor_Main_Immune_Mimicry_Analysis.R
05_Additional_Analyses_All_Breast_Tumors.R
</pre>

##### **Evaluate Immune Mimicry in Reduction Mammoplastiess** #####
<pre>
04_Kumar_et_al_Pal_et_al_Reduction_Mammoplasty_Main_Immune_Mimicry_Analysis.R
</pre>

##### **Assess Mammary Keratins in Peripheral Blood Cells** #####
<pre>
01_Terekhova_et_al_PMBC_Mammary_Keratin_Quantification.R
</pre>

