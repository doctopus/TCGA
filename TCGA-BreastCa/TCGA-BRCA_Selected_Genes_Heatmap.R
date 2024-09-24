# Load required libraries
library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(tibble)
library(ComplexHeatmap)
library(circlize)

# Define target directory for TCGAbiolinks to download data
tcga_download_directory <- "/Users/i/Dropbox/Clinic3.0/Developer/RStudio/TCGA/TCGA-BreastCa/GDCdata"

# List of mouse genes
mouse_genes <- c("Wnt5b", "Wnt11", "Snai2", "Ngf", "Rspo3", "Adamts5",
                 "Wnt3", "Wnt6", "Fzd10", "Wnt9b", "Sox17",
                 "Adrb1", "Adrb2")
genes_of_interest <- toupper(mouse_genes)

# Download TCGA BRCA gene expression data
query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

GDCdownload(query, directory = tcga_download_directory)
exp_data <- GDCprepare(query)

# Extract gene expression data
gene_expression <- assay(exp_data, "tpm_unstrand")  # Using TPM values

# Get gene symbols
gene_symbols <- rowData(exp_data)$gene_name

# Create a data frame with gene symbols and expression data
expression_df <- data.frame(
  gene_symbol = gene_symbols,
  gene_expression,
  stringsAsFactors = FALSE,
  row.names = NULL  # Ensure no row names are present
)

# Filter for genes of interest
filtered_expression <- expression_df %>%
  filter(gene_symbol %in% genes_of_interest) %>%
  
  mutate(gene_symbol = factor(gene_symbol, levels = genes_of_interest)) %>%
  arrange(gene_symbol) %>%
  
  tibble::column_to_rownames("gene_symbol")

# Convert filtered_expression to a matrix
filtered_expression_matrix <- as.matrix(filtered_expression)

# Check if the matrix is numeric
if (!is.numeric(filtered_expression_matrix)) {
  stop("The expression data is not numeric. Please check your input data.")
}

# Replace any non-finite values with NA
filtered_expression_matrix[!is.finite(filtered_expression_matrix)] <- NA

# Log2 transform the expression values, adding a small constant to avoid log(0)
log2_expression <- log2(filtered_expression_matrix + 1)

# Remove any rows or columns that are all NA
log2_expression <- log2_expression[rowSums(is.na(log2_expression)) != ncol(log2_expression), 
                                   colSums(is.na(log2_expression)) != nrow(log2_expression)]

# Calculate color breaks, excluding NA values
col_breaks <- c(min(log2_expression, na.rm = TRUE), 
                mean(log2_expression, na.rm = TRUE), 
                max(log2_expression, na.rm = TRUE))

# Create the heatmap
heatmap <- Heatmap(
  log2_expression,
  name = "Log2(TPM + 1)",
  cluster_rows = FALSE,  # Ensure rows are not clustered
  cluster_columns = TRUE,
  show_column_names = FALSE,
  show_row_names = TRUE,
  row_order = rownames(log2_expression),  # Explicitly set row order
  column_title = "TCGA-BRCA Breast Cancer Samples",
  row_title = "Genes of Interest",
  col = colorRamp2(col_breaks, c("blue", "white", "red")),
  na_col = "grey"  # Color for NA values
)
# Draw the heatmap
pdf("tcga_brca_heatmap_4.pdf", width = 12, height = 8)
draw(heatmap)
dev.off()


# Print the order of genes in the final data
print(rownames(log2_expression))

# Print summary of the data to check for any issues
print(summary(log2_expression))

# Print the structure of the data
print(str(log2_expression))
