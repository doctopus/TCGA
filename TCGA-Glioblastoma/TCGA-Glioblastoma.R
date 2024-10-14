#To Plot the survival curve of males and females in TCGA-GBM database

# Install and load required packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("TCGAbiolinks", "SummarizedExperiment", "biomaRt", "survival", "survminer", "ggplot2", "pheatmap", "reshape2"))

library(TCGAbiolinks)
library(SummarizedExperiment)
library(biomaRt)
library(dplyr)
library(survival)
library(survminer)
library(ggplot2)
library(pheatmap)
library(reshape2)

#Find Human Ortholog Genes from Mouse Genes----
# List of mouse genes
mouse_genes <- c("Sox2", "Adam11", "Dkk4", "Wnt5b", "Wnt11", "Snai2", "Ngf", "Rspo3", "Adamts5",
                 "Wnt3", "Wnt6", "Fzd10", "Wnt9b", "Sox17",
                 "Adrb1", "Adrb2", "Adra1a", "Adra1b", "Adra1d", "Adra2a", "Adra2b", "Adra2c")
genes_of_interest <- toupper(mouse_genes)
# Set up Ensembl connections
ensembl_mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/") #Newer version will cause errors in getLDS step
ensembl_human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/") #Newer version will cause errors in getLDS step

# Get mouse gene symbol and Ensemble IDs
mouse_genes_info <- getBM(attributes = c("mgi_symbol", "ensembl_gene_id"),
                          filters = "mgi_symbol",
                          values = mouse_genes,
                          mart = ensembl_mouse)

# Get human orthologs from the ensemble IDs of the mouse genes
human_orthologs <- getLDS(attributes = c("ensembl_gene_id", "mgi_symbol"),
                          filters = "ensembl_gene_id",
                          values = mouse_genes_info$ensembl_gene_id,
                          mart = ensembl_mouse,
                          attributesL = c("hgnc_symbol", "ensembl_gene_id"),
                          martL = ensembl_human)

# Merge to get a clean list of mouse to human gene symbols
gene_mapping <- merge(mouse_genes_info, human_orthologs, by.x = "ensembl_gene_id", by.y = "Gene.stable.ID")
gene_mapping <- gene_mapping[, c("mgi_symbol", "HGNC.symbol")]
colnames(gene_mapping) <- c("mouse_gene", "human_gene")

print(gene_mapping)

#SKIP-START IF NOT SAVING DATA (Not functional)----
#Symlink GDCdata folder to store data downloaded from TCGAbiolinks outside the 
# first link is where data is, second link is where symlink is going to be
system("ln -s /Users/i/Documents/Clinic3.0/TCGA/TCGA-Glioblastoma/GDCdata /Users/i/Dropbox/Clinic3.0/Developer/RStudio/TCGA/TCGA-Glioblastoma/GDCdata")
#Test if it works by adding a file to the folderr where the data will go through the symlink
list.files("/Users/i/Dropbox/Clinic3.0/Developer/RStudio/TCGA/TCGA-Glioblastoma/GDCdata")
 
#Define Target directory for TCGAbiolinks to download data
tcga_download_directory <- "/Users/i/Dropbox/Clinic3.0/Developer/RStudio/TCGA/TCGA-Glioblastoma/GDCdata"
# tcga_download_directory <- "/Users/i/Documents/Clinic3.0/TCGA/TCGA-BreastCa/GDCdata"

# Download TCGA BRCA gene expression data
query <- GDCquery(project = "TCGA-GBM", 
                  data.category = "Clinical")

GDCdownload(query, directory = tcga_download_directory)
View(getResults(query))
#SKIP-END IF NOT SAVING DATA----

#SURVIVAL PLOT-START-Download and plot the clinical data for Survival Plot ----
clinical_data <- GDCquery_clinic(project = "TCGA-GBM",
                                 type = "clinical",
                                 save.csv = FALSE)

# Filter columns for survival plot
survival_data <- clinical_data[, c("bcr_patient_barcode", "submitter_id", "primary_diagnosis", "gender", "vital_status", "days_to_death", "days_to_last_follow_up")]

# survival_data <- survival_data
# Create a new column for survival time
survival_data$survival_time <- ifelse(survival_data$vital_status == "Dead", survival_data$days_to_death, survival_data$days_to_last_follow_up)

# Create a new column for censoring status
survival_data$censor <- ifelse(survival_data$vital_status == "Alive", 0, 1)

# Convert survival time and censor to numeric
survival_data$survival_time <- as.numeric(survival_data$survival_time)
survival_data$censor <- as.numeric(survival_data$censor)

# Create a new column for gender
# survival_data$gender <- substr(survival_data$bcr_patient_barcode, 14, 14)

# Plot Kaplan-Meier curve comparing survival between male and female patients
fit <- survfit(Surv(survival_time, censor) ~ gender, data = survival_data)
ggsurvplot(fit, data = survival_data, pval = TRUE)

#SURVIVAL PLOT-END ----

