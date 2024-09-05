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

# List of mouse genes
mouse_genes <- c("Sox2", "Adam11", "Dkk4", "Wnt5b", "Wnt11", "Snai2", "Ngf", "Rspo3", "Adamts5",
                 "Wnt3", "Wnt6", "Fzd10", "Wnt9b", "Sox17")

# Set up Ensembl connections
ensembl_mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/") #Newer version will cause errors in getLDS step
ensembl_human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/") #Newer version will cause errors in getLDS step

# Get mouse gene IDs
mouse_genes_info <- getBM(attributes = c("mgi_symbol", "ensembl_gene_id"),
                          filters = "mgi_symbol",
                          values = mouse_genes,
                          mart = ensembl_mouse)

# Get human orthologs
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


#Symlink GDCdata folder to store data downloaded from TCGAbiolinks outside the 
# first link is where data is, second link is where symlink is
system("ln -s /Users/i/Documents/Clinic3.0/TCGA/TCGA-BreastCa/GDCdata /Users/i/Dropbox/Clinic3.0/Developer/RStudio/TCGA/TCGA-BreastCa/GDCdata")
#Test if it works
list.files("/Users/i/Dropbox/Clinic3.0/Developer/RStudio/TCGA/TCGA-BreastCa/GDCdata")

#Define Target directory for TCGAbiolinks to download data
tcga_download_directory <- "/Users/i/Dropbox/Clinic3.0/Developer/RStudio/TCGA/TCGA-BreastCa/GDCdata"

# Download TCGA BRCA data
query <- GDCquery(project = "TCGA-BRCA",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "STAR - Counts")

GDCdownload(query, directory = tcga_download_directory)

data <- GDCprepare(query)

# Filter for TNBC samples (you may need to adjust this based on the specific clinical data available)
clinical <- colData(data)

# tnbc_samples <- clinical$primary_diagnosis == c("Infiltrating duct carcinoma, NOS", "Infiltrating duct mixed with other types of carcinoma")
# data_tnbc <- data[, tnbc_samples]


tnbc_samples <- clinical$primary_diagnosis %in% c("Infiltrating duct carcinoma, NOS", "Infiltrating duct mixed with other types of carcinoma")
data_tnbc <- data[, tnbc_samples]

# unique_primary_diagnosis <- unique(clinical$primary_diagnosis)
# print(unique_primary_diagnosis)

# Extract expression data for human orthologs
human_genes <- gene_mapping$human_gene
expr_data <- assay(data_tnbc)[rowData(data_tnbc)$gene_name %in% human_genes, ]

# If some genes are missing, print a warning
missing_genes <- human_genes[!human_genes %in% rowData(data_tnbc)$gene_name]
if (length(missing_genes) > 0) {
  warning(paste("The following genes were not found in the TCGA data:", paste(missing_genes, collapse = ", ")))
}

# Get survival data
survival_data <- data.frame(
  time = clinical$days_to_last_follow_up[tnbc_samples],
  status = clinical$vital_status[tnbc_samples],
  stringsAsFactors = FALSE
)
survival_data$status <- ifelse(survival_data$status == "Alive", 0, 1)

# Save the results for use in visualization
save(expr_data, survival_data, gene_mapping, file = "tcga_analysis_data.RData")

print("Analysis complete. Data saved in tcga_analysis_data.RData")

#### Optional ----
genes_of_interest = c("ADAMTS5", "WNT3", "WNT9B", "WNT11")
# Perform survival analysis for each gene
for (gene in genes_of_interest) {
  fit <- survfit(Surv(time, status) ~ ifelse(expr_data[gene, ] > median(expr_data[gene, ]), "High", "Low"), data = survival_data)
  plot <- ggsurvplot(fit, data = survival_data, pval = TRUE, risk.table = TRUE, title = paste("Survival analysis for", gene))
  print(plot)
}

# Logistic regression for metastasis (assuming you have metastasis data)
# You'll need to create a metastasis variable based on available clinical data
metastasis <- clinical$metastasis_status  # Replace with actual metastasis data

for (gene in genes_of_interest) {
  model <- glm(metastasis ~ expr_data[gene, ], family = binomial())
  print(summary(model))
}

####Plot Data----
#Inputs: expr_data, survival_data, and genes_of_interest
library(ggplot2)
library(survminer)
library(pheatmap)
library(reshape2)
# Load the data we saved earlier
load("tcga_analysis_data.RData")

# 1. Heatmap of gene expression
pheatmap(expr_data, 
         scale = "row", 
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         main = "Gene Expression Heatmap",
         fontsize_row = 10,
         fontsize_col = 8)


# 2. Box plot of gene expression
expr_data_long <- reshape2::melt(expr_data)
colnames(expr_data_long) <- c("Gene", "Sample", "Expression")

ggplot(expr_data_long, aes(x = Gene, y = Expression)) +
  geom_boxplot() +
  theme_bw() +
  ggtitle("Gene Expression Distribution") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 3. Correlation plot between genes
cor_matrix <- cor(t(expr_data))
pheatmap(cor_matrix, 
         main = "Gene Expression Correlation",
         fontsize = 10)


# 4. Survival curves for high/low expression of each gene
for (gene in rownames(expr_data)) {
  expr_median <- median(expr_data[gene, ])
  survival_data$expr_group <- ifelse(expr_data[gene, ] > expr_median, "High", "Low")
  
  fit <- survfit(Surv(time, status) ~ expr_group, data = survival_data)
  
  p <- ggsurvplot(fit, 
                  data = survival_data, 
                  pval = TRUE, 
                  conf.int = TRUE,
                  risk.table = TRUE, 
                  risk.table.col = "strata",
                  linetype = "strata",
                  surv.median.line = "hv",
                  ggtheme = theme_bw(),
                  palette = c("#E7B800", "#2E9FDF"),
                  title = paste("Survival analysis for", gene))
  
  print(p)
}

# 5. Differential Expression Analysis and Volcano Plot
# Create a DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = round(expr_data),
                              colData = data.frame(condition = survival_data$status),
                              design = ~ condition)

# Run DESeq2
dds <- DESeq(dds)
res <- results(dds)

# Create Volcano plot
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Differential Expression',
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 3.0,
                labSize = 6.0)


# 6. PCA plot
vst_data <- vst(dds, blind = FALSE)
plotPCA(vst_data, intgroup = "condition") +
  ggtitle("PCA of Gene Expression Data")

print("Visualization complete. Please check the generated plots.")