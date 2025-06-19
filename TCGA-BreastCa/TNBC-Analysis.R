#   Alternative Methods to Identify TNBC in TCGA Data
# 1. Using PAM50 Molecular Subtypes----
# The basal-like subtype is strongly associated with TNBC:
# 
# Check available clinical variables
colnames(clinical_data)

# Look for PAM50 or other molecular subtype classifications
subtype_cols <- grep("subtype|PAM50|molecular_subtype", colnames(clinical_data), 
                     ignore.case = TRUE, value = TRUE)
print(subtype_cols)

# If PAM50 is available, use basal subtype as proxy for TNBC
if(length(subtype_cols) > 0) {
  # Choose the most appropriate column
  subtype_col <- subtype_cols[1]  # Adjust if needed
  
  # Print distribution of subtypes
  print(table(clinical_data[[subtype_col]], useNA = "ifany"))
  
  # Identify TNBC as basal-like
  clinical_data$TNBC_PAM50 <- grepl("Basal", clinical_data[[subtype_col]], ignore.case = TRUE)
  
  # Print TNBC count
  print(paste("TNBC samples (basal):", sum(clinical_data$TNBC_PAM50, na.rm = TRUE)))
}


print(paste("TNBC samples (PAM50-inferred):", sum(clinical_data$TNBC_PAM50, na.rm = TRUE)))

clinical_data$NON_TNBC_PAM50 <- !is.na(clinical_data$paper_BRCA_Subtype_PAM50) & clinical_data$paper_BRCA_Subtype_PAM50 != "Basal"
# clinical_data <- clinical_data %>% mutate(NON_TNBC_PAM50 = !is.na(paper_BRCA_Subtype_PAM50) & paper_BRCA_Subtype_PAM50 != "Basal")
print(paste("Non-TNBC samples (PAM50-inferred):", sum(clinical_data$NON_TNBC_PAM50, na.rm = TRUE)))



# 2. Using RNA-seq Data to Infer Receptor Status----
# Map Ensembl IDs to gene symbols: Create a complete mapping first
# Function to Map
create_gene_mapping <- function(expression_data) {
  # Get the expression data
  expr_data <- assay(expression_data, "tpm_unstrand")
  
  # Extract Ensembl IDs (remove version numbers if present)
  ensembl_ids <- rownames(expr_data)
  ensembl_ids_clean <- gsub("\\..*", "", ensembl_ids)
  
  # Set up biomaRt connection
  mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  
  # Map all Ensembl IDs to gene symbols
  gene_mapping <- getBM(
    attributes = c("ensembl_gene_id", "external_gene_name"),
    filters = "ensembl_gene_id",
    values = unique(ensembl_ids_clean),
    mart = mart
  )
  
  # Create a named vector for quick lookup
  mapping_vector <- setNames(gene_mapping$external_gene_name, gene_mapping$ensembl_gene_id)
  
  # Add gene symbols to expression data
  gene_symbols <- mapping_vector[ensembl_ids_clean]
  
  # Create a data frame with mapping information
  mapping_df <- data.frame(
    ensembl_id = ensembl_ids,
    ensembl_id_clean = ensembl_ids_clean,
    gene_symbol = gene_symbols,
    stringsAsFactors = FALSE
  )
  
  return(mapping_df)
}

# Create complete mapping
gene_mapping_df <- create_gene_mapping(BRCA_data)

# Function to extract expression by gene symbol
extract_by_symbol <- function(gene_symbol, mapping_df, expression_matrix) {
  # Find rows matching the gene symbol
  matching_rows <- which(mapping_df$gene_symbol == gene_symbol)
  
  if (length(matching_rows) == 0) {
    warning(paste("Gene", gene_symbol, "not found"))
    return(NULL)
  }
  
  if (length(matching_rows) > 1) {
    warning(paste("Multiple entries found for", gene_symbol, ". Using the first one."))
    print(paste("Ensembl IDs found:", paste(mapping_df$ensembl_id[matching_rows], collapse = ", ")))
    matching_rows <- matching_rows[1]
  }
  
  # Get the corresponding Ensembl ID
  ensembl_id <- mapping_df$ensembl_id[matching_rows]
  
  # Extract expression
  return(expression_matrix[ensembl_id, ])
}


# Extract expression for receptor genes
print("Extracting expression for receptor genes...")

esr1_expr <- extract_by_symbol("ESR1", gene_mapping_df, expression_data)
pgr_expr <- extract_by_symbol("PGR", gene_mapping_df, expression_data)
erbb2_expr <- extract_by_symbol("ERBB2", gene_mapping_df, expression_data)


# Check if all genes were found
if (is.null(esr1_expr)) stop("ESR1 not found")
if (is.null(pgr_expr)) stop("PGR not found")
if (is.null(erbb2_expr)) stop("ERBB2 not found")

print("Successfully extracted expression for all receptor genes")

# Plot distributions to determine thresholds
par(mfrow=c(3,1))
hist(log2(esr1_expr+1), breaks=50, main="ESR1 Expression", xlab="log2(TPM+1)")
abline(v=median(log2(esr1_expr+1)), col="red", lty=2)

hist(log2(pgr_expr+1), breaks=50, main="PGR Expression", xlab="log2(TPM+1)")
abline(v=median(log2(pgr_expr+1)), col="red", lty=2)

hist(log2(erbb2_expr+1), breaks=50, main="ERBB2 Expression", xlab="log2(TPM+1)")
abline(v=median(log2(erbb2_expr+1)), col="red", lty=2)

# Print summary statistics to help set thresholds
print("Expression summary statistics:")
print(paste("ESR1 - Median:", round(median(log2(esr1_expr+1)), 2), 
            "Mean:", round(mean(log2(esr1_expr+1)), 2),
            "75th percentile:", round(quantile(log2(esr1_expr+1), 0.75), 2)))
print(paste("PGR - Median:", round(median(log2(pgr_expr+1)), 2), 
            "Mean:", round(mean(log2(pgr_expr+1)), 2),
            "75th percentile:", round(quantile(log2(pgr_expr+1), 0.75), 2)))
print(paste("ERBB2 - Median:", round(median(log2(erbb2_expr+1)), 2), 
            "Mean:", round(mean(log2(erbb2_expr+1)), 2),
            "75th percentile:", round(quantile(log2(erbb2_expr+1), 0.75), 2)))

# Set thresholds based on distributions (adjust these based on your histograms)
# Some approaches
# 1. Literature-based thresholds
# 2. Median split
# 3. Quantile-based (e.g., 75th percentile)
# 4. Mixture model-based thresholds
# USING MEDIAN Approach
esr1_threshold <- median(log2(esr1_expr+1))
pgr_threshold <- median(log2(pgr_expr+1))
erbb2_threshold <- median(log2(erbb2_expr+1))

# OR USING THRESHOLD Approach (Custom Defined)
# esr1_threshold <- 1.5  # log2(TPM+1)
# pgr_threshold <- 1.5   # log2(TPM+1)
# erbb2_threshold <- 4   # log2(TPM+1)

# OR more stringent thresholds (75th percentile)
# esr1_threshold <- quantile(log2(esr1_expr+1), 0.75)
# pgr_threshold <- quantile(log2(pgr_expr+1), 0.75)
# erbb2_threshold <- quantile(log2(erbb2_expr+1), 0.75)

# Define threshold for Non-TNBC (25th percentile)
esr1_threshold_control <- quantile(log2(esr1_expr+1), 0.75)
pgr_threshold_control <- quantile(log2(pgr_expr+1), 0.75)
erbb2_threshold_control <- quantile(log2(erbb2_expr+1), 0.75)



print(paste("Using thresholds - ESR1:", round(esr1_threshold, 2), 
            "PGR:", round(pgr_threshold, 2), 
            "ERBB2:", round(erbb2_threshold, 2)))

# Classify Receptor Positive samples
clinical_data$ER_status_RNA <- log2(esr1_expr+1) > esr1_threshold
clinical_data$PR_status_RNA <- log2(pgr_expr+1) > pgr_threshold
clinical_data$HER2_status_RNA <- log2(erbb2_expr+1) > erbb2_threshold

# Classify Control with High Receptor Positive samples
clinical_data$ER_status_RNA_Control <- log2(esr1_expr+1) > esr1_threshold_control
clinical_data$PR_status_RNA_Control <- log2(pgr_expr+1) > pgr_threshold_control
clinical_data$HER2_status_RNA_Control <- log2(erbb2_expr+1) > erbb2_threshold_control

# Print classification summary
print("Receptor status classification:")
print(paste("ER positive:", sum(clinical_data$ER_status_RNA, na.rm = TRUE)))
print(paste("PR positive:", sum(clinical_data$PR_status_RNA, na.rm = TRUE)))
print(paste("HER2 positive:", sum(clinical_data$HER2_status_RNA, na.rm = TRUE)))

# Define TNBC as triple-negative based on RNA
clinical_data$TNBC_RNA <- !clinical_data$ER_status_RNA & 
  !clinical_data$PR_status_RNA & 
  !clinical_data$HER2_status_RNA

# Define Non-TNBC as Receptor expression Above 75 quartile threshold based on RNA (Use OR)
clinical_data$NON_TNBC_RNA <- clinical_data$ER_status_RNA_Control | 
  clinical_data$PR_status_RNA_Control | 
  clinical_data$HER2_status_RNA_Control

# Print TNBC count
print(paste("TNBC samples (RNA-inferred):", sum(clinical_data$TNBC_RNA, na.rm = TRUE)))
print(paste("Non-TNBC samples (RNA-inferred):", sum(clinical_data$NON_TNBC_RNA, na.rm = TRUE)))

# Save the gene mapping for future use
# write.csv(gene_mapping_df, "gene_mapping_ensembl_to_symbol.csv", row.names = FALSE)


# 3. Using Published TNBC Classifications----
# Several papers have published curated TNBC classifications for TCGA:

# Check if TCGA biospecimen barcode contains sample type information
clinical_data$patient_id <- substr(rownames(clinical_data), 1, 12)

# Load published TNBC classifications (this is an example - you would need to get this data)
# You could download this from published supplementary materials
published_tnbc <- read.csv("published_tnbc_classifications.csv", stringsAsFactors = FALSE)

# Merge with your clinical data
clinical_data$TNBC_published <- clinical_data$patient_id %in% published_tnbc$patient_id

# 4. Using The Cancer Proteome Atlas (TCPA) Data----
# TCPA has protein expression data including ER/PR/HER2:
# 

# Install TCPA package if needed
if (!requireNamespace("TCPA", quietly = TRUE)) {
  BiocManager::install("TCPA")
}
library(TCPA)

# Load TCPA data
tcpa_data <- TCPA.loadRPPA("BRCA")

# Extract protein expression for receptors
er_protein <- tcpa_data$data["ER-alpha",]
pr_protein <- tcpa_data$data["PR",]  # May need to check exact names
her2_protein <- tcpa_data$data["HER2",]

# Set threshold for protein expression
er_thresh <- median(er_protein)
pr_thresh <- median(pr_protein)
her2_thresh <- median(her2_protein)

# Match samples and assign TNBC status
# [Code to match TCPA sample IDs with your clinical data]


# 5. Combining Multiple Methods----
# Select Superset of TNBC_PAM50 and TNBC_RNA and Overlap of NON_TNBC_PAM50 and NON_TNBC_RNA

#Create TNBC Column Mark TNBC as TRUE and The NON_TNBC_RNA as FLASE and rest as NA; Comparison will be between TNBC TRUE vs FALSE
clinical_data$TNBC <- NA
clinical_data$TNBC[clinical_data$TNBC_PAM50 == TRUE | clinical_data$TNBC_RNA == TRUE] <- TRUE
clinical_data$TNBC[clinical_data$NON_TNBC_PAM50 == TRUE & clinical_data$NON_TNBC_RNA == TRUE] <- FALSE

table(clinical_data$TNBC, useNA = "always")

# Filter out samples with unknown TNBC status
clinical_filtered <- clinical_data[!is.na(clinical_data$TNBC), ]
expression_filtered <- expression_data[, rownames(clinical_filtered)]

# Print final counts
print(paste("Total BRCA samples in TCGA:", nrow(clinical_data)))
print(paste("TNBC Analysis samples TNBC + Non_TNBC:", nrow(clinical_filtered)))

# 3. DIFFERENTIAL EXPRESSION ANALYSIS -------------------------------------
# Function to perform differential expression analysis
perform_deg_analysis <- function(expression_data, clinical_data, group_column) {
  # Create DGEList object (for edgeR analysis)
  dge <- DGEList(counts = expression_data)
  
  # Filter low-expression genes
  keep <- filterByExpr(dge, group = clinical_data[[group_column]])
  dge <- dge[keep, ]
  
  # Normalize the data
  dge <- calcNormFactors(dge)
  
  # Create design matrix dynamically
  design_formula <- as.formula(paste("~", group_column))
  design <- model.matrix(design_formula, data = clinical_data)
  rownames(design) <- rownames(clinical_data)
  
  # Estimate dispersion
  dge <- estimateDisp(dge, design)
  
  # Fit model
  fit <- glmQLFit(dge, design)
  
  # Test for differential expression
  qlf <- glmQLFTest(fit, coef = 2)  # Testing group vs control
  
  # Get results table
  de_results <- topTags(qlf, n = nrow(dge$counts), sort.by = "PValue")
  de_table <- as.data.frame(de_results)
  
  return(de_table)
}

# Function to add gene symbols to DE results using existing mapping
add_gene_symbols_to_de <- function(de_table, gene_mapping_df) {
  # Add ensembl_id column (clean version without .version)
  de_table$ensembl_id_clean <- gsub("\\..*", "", rownames(de_table))
  
  # Merge with gene mapping
  de_table_with_symbols <- merge(
    de_table, 
    gene_mapping_df[, c("ensembl_id_clean", "gene_symbol")], 
    by = "ensembl_id_clean", 
    all.x = TRUE
  )
  
  # Add back the original ensembl_id with version
  de_table_with_symbols$ensembl_id_original <- rownames(de_table)[match(de_table_with_symbols$ensembl_id_clean, gsub("\\..*", "", rownames(de_table)))]
  
  # Reorder columns for better readability
  col_order <- c("ensembl_id_original", "ensembl_id_clean", "gene_symbol", 
                 "logFC", "logCPM", "F", "PValue", "FDR")
  de_table_with_symbols <- de_table_with_symbols[, col_order]
  
  return(de_table_with_symbols)
}

# Function to filter DE results by custom gene list
filter_de_by_genes <- function(de_table_with_symbols, gene_list, filter_type = "symbol") {
  if (filter_type == "symbol") {
    # Filter by gene symbols
    filtered_results <- de_table_with_symbols[de_table_with_symbols$gene_symbol %in% gene_list, ]
  } else if (filter_type == "ensembl") {
    # Filter by Ensembl IDs
    filtered_results <- de_table_with_symbols[de_table_with_symbols$ensembl_id_clean %in% gene_list, ]
  } else {
    stop("filter_type must be either 'symbol' or 'ensembl'")
  }
  
  # Remove rows where gene_symbol is NA (if filtering by symbol)
  if (filter_type == "symbol") {
    filtered_results <- filtered_results[!is.na(filtered_results$gene_symbol), ]
  }
  
  # Sort by p-value
  filtered_results <- filtered_results[order(filtered_results$PValue), ]
  
  return(filtered_results)
}

# Function to get top N DE genes
get_top_de_genes <- function(de_table_with_symbols, n = 20, criteria = "PValue") {
  # Remove rows with NA gene symbols
  de_clean <- de_table_with_symbols[!is.na(de_table_with_symbols$gene_symbol), ]
  
  # Sort by specified criteria
  if (criteria == "PValue") {
    top_genes <- de_clean[order(de_clean$PValue), ]
  } else if (criteria == "FDR") {
    top_genes <- de_clean[order(de_clean$FDR), ]
  } else if (criteria == "logFC") {
    top_genes <- de_clean[order(abs(de_clean$logFC), decreasing = TRUE), ]
  } else {
    stop("criteria must be 'PValue', 'FDR', or 'logFC'")
  }
  
  return(head(top_genes, n))
}

# Function to create volcano plot with highlighted genes
create_volcano_plot <- function(de_table_with_symbols, highlight_genes = NULL, 
                                pval_threshold = 0.05, fc_threshold = 1,
                                title = "Volcano Plot") {
  library(ggplot2)
  
  # Prepare data for plotting
  plot_data <- de_table_with_symbols
  plot_data$significant <- (plot_data$PValue < pval_threshold) & (abs(plot_data$logFC) > fc_threshold)
  plot_data$neg_log10_pval <- -log10(plot_data$PValue)
  
  # Create basic plot
  p <- ggplot(plot_data, aes(x = logFC, y = neg_log10_pval)) +
    geom_point(aes(color = significant), alpha = 0.6) +
    scale_color_manual(values = c("FALSE" = "grey", "TRUE" = "red")) +
    geom_vline(xintercept = c(-fc_threshold, fc_threshold), linetype = "dashed", alpha = 0.5) +
    geom_hline(yintercept = -log10(pval_threshold), linetype = "dashed", alpha = 0.5) +
    labs(
      title = title,
      x = "Log2 Fold Change",
      y = "-Log10 P-value",
      color = "Significant"
    ) +
    theme_minimal()
  
  # Highlight specific genes if provided
  if (!is.null(highlight_genes)) {
    highlight_data <- plot_data[plot_data$gene_symbol %in% highlight_genes, ]
    if (nrow(highlight_data) > 0) {
      p <- p + 
        geom_point(data = highlight_data, aes(x = logFC, y = neg_log10_pval), 
                   color = "blue", size = 3, alpha = 0.8) +
        ggrepel::geom_text_repel(data = highlight_data, 
                                 aes(label = gene_symbol),
                                 size = 3, max.overlaps = Inf)
    }
  }
  
  return(p)
}
#Better Volcano Plot---================
# Enhanced volcano plot function with group labels (CORRECTED)
create_volcano_plot <- function(de_table_with_symbols, highlight_genes = NULL, 
                                pval_threshold = 0.05, fc_threshold = 1,
                                title = "Volcano Plot",
                                group1_name = "TNBC", group2_name = "non-TNBC") {
  library(ggplot2)
  library(ggrepel)
  
  # Prepare data for plotting
  plot_data <- de_table_with_symbols
  plot_data$significant <- (plot_data$PValue < pval_threshold) & (abs(plot_data$logFC) > fc_threshold)
  plot_data$neg_log10_pval <- -log10(plot_data$PValue)
  
  # Create significance categories for better coloring
  plot_data$category <- "Not Significant"
  plot_data$category[plot_data$logFC > fc_threshold & plot_data$PValue < pval_threshold] <- "Higher_in_group1"
  plot_data$category[plot_data$logFC < -fc_threshold & plot_data$PValue < pval_threshold] <- "Higher_in_group2"
  
  # Define colors using separate variables
  higher_group1_label <- paste("Higher in", group1_name)
  higher_group2_label <- paste("Higher in", group2_name)
  
  colors <- c("Not Significant" = "grey70",
              "Higher_in_group1" = "#d73027",    # Red for group1 (TNBC)
              "Higher_in_group2" = "#1a9850")    # Green for group2 (non-TNBC)
  
  # Create labels for legend
  legend_labels <- c("Not Significant" = "Not Significant",
                     "Higher_in_group1" = higher_group1_label,
                     "Higher_in_group2" = higher_group2_label)
  
  # Create basic plot
  p <- ggplot(plot_data, aes(x = logFC, y = neg_log10_pval)) +
    geom_point(aes(color = category), alpha = 0.7, size = 1) +
    scale_color_manual(values = colors, labels = legend_labels, name = "Expression") +
    geom_vline(xintercept = c(-fc_threshold, fc_threshold), 
               linetype = "dashed", alpha = 0.5, color = "black") +
    geom_hline(yintercept = -log10(pval_threshold), 
               linetype = "dashed", alpha = 0.5, color = "black") +
    labs(
      title = title,
      subtitle = paste("Comparison:", group1_name, "vs", group2_name),
      x = paste("Log2 Fold Change\n←", group2_name, "higher | ", group1_name, "higher →"),
      y = "-Log10 P-value"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      axis.title.x = element_text(size = 11),
      axis.title.y = element_text(size = 11),
      legend.position = "bottom",
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 10)
    )
  
  # Add count annotations
  n_higher_group1 <- sum(plot_data$category == "Higher_in_group1")
  n_higher_group2 <- sum(plot_data$category == "Higher_in_group2")
  
  # Add text annotations showing counts
  max_y <- max(plot_data$neg_log10_pval, na.rm = TRUE)
  max_x <- max(abs(plot_data$logFC), na.rm = TRUE)
  
  p <- p + 
    annotate("text", x = max_x * 0.7, y = max_y * 0.95, 
             label = paste("n =", n_higher_group1), 
             color = "#d73027", 
             size = 4, fontface = "bold") +
    annotate("text", x = -max_x * 0.7, y = max_y * 0.95, 
             label = paste("n =", n_higher_group2), 
             color = "#1a9850", 
             size = 4, fontface = "bold")
  
  # Highlight specific genes if provided
  if (!is.null(highlight_genes)) {
    highlight_data <- plot_data[plot_data$gene_symbol %in% highlight_genes & 
                                  !is.na(plot_data$gene_symbol), ]
    if (nrow(highlight_data) > 0) {
      p <- p + 
        geom_point(data = highlight_data, 
                   aes(x = logFC, y = neg_log10_pval), 
                   color = "blue", size = 3, alpha = 0.9, shape = 21, 
                   fill = "lightblue", stroke = 1.5) +
        geom_text_repel(data = highlight_data, 
                        aes(label = gene_symbol),
                        size = 3.5, 
                        max.overlaps = Inf,
                        box.padding = 0.5,
                        point.padding = 0.3,
                        segment.color = "blue",
                        segment.alpha = 0.7,
                        fontface = "bold")
    }
  }
  
  return(p)
}

# Summary function
print_de_summary <- function(de_table_with_symbols, group1_name = "TNBC", group2_name = "non-TNBC",
                             pval_threshold = 0.05, fc_threshold = 1) {
  
  # Calculate summary statistics
  total_genes <- nrow(de_table_with_symbols)
  sig_genes <- sum(de_table_with_symbols$PValue < pval_threshold, na.rm = TRUE)
  sig_up_group1 <- sum(de_table_with_symbols$logFC > fc_threshold & 
                         de_table_with_symbols$PValue < pval_threshold, na.rm = TRUE)
  sig_down_group1 <- sum(de_table_with_symbols$logFC < -fc_threshold & 
                           de_table_with_symbols$PValue < pval_threshold, na.rm = TRUE)
  
  cat("=== DIFFERENTIAL EXPRESSION SUMMARY ===\n")
  cat("Comparison:", group1_name, "vs", group2_name, "\n")
  cat("Total genes tested:", total_genes, "\n")
  cat("Significance threshold: p <", pval_threshold, ", |logFC| >", fc_threshold, "\n")
  cat("Total significant genes:", sig_genes, "\n")
  cat("Higher in", group1_name, "(logFC >", fc_threshold, "):", sig_up_group1, "\n")
  cat("Higher in", group2_name, "(logFC <", -fc_threshold, "):", sig_down_group1, "\n")
  cat("========================================\n\n")
}

# Corrected group-specific genes function
get_group_specific_genes <- function(de_table_with_symbols, group1_name = "TNBC", group2_name = "non-TNBC",
                                     pval_threshold = 0.05, fc_threshold = 1, top_n = 10) {
  
  # Remove NA gene symbols
  de_clean <- de_table_with_symbols[!is.na(de_table_with_symbols$gene_symbol), ]
  
  # Genes higher in group1 (positive logFC)
  higher_group1 <- de_clean[de_clean$logFC > fc_threshold & 
                              de_clean$PValue < pval_threshold, ]
  higher_group1 <- higher_group1[order(higher_group1$logFC, decreasing = TRUE), ]
  
  # Genes higher in group2 (negative logFC)
  higher_group2 <- de_clean[de_clean$logFC < -fc_threshold & 
                              de_clean$PValue < pval_threshold, ]
  higher_group2 <- higher_group2[order(higher_group2$logFC), ]
  
  result <- list(
    higher_in_group1 = head(higher_group1, top_n),
    higher_in_group2 = head(higher_group2, top_n),
    group1_name = group1_name,
    group2_name = group2_name
  )
  
  return(result)
}



# MAIN ANALYSIS EXECUTION -------------------------------------------------
library(readxl)
library(tidyverse)
patients_tnbc_nature <- read_excel("TCGA_TNBC_Nature.xlsx") 
clinical_filtered_nature <- clinical_filtered %>%
  filter(TNBC == FALSE | 
           (TNBC == TRUE & patient %in% patients_tnbc_nature$bcr_patient_barcode))

expression_filtered_nature <- expression_filtered %>% select(all_of(clinical_filtered_nature$barcode)) #Does not work as it is matrix
expression_filtered_nature <- expression_filtered[,colnames(expression_filtered) %in% clinical_filtered_nature$barcode]

clinical_filtered %>% count(TNBC)
clinical_filtered_nature %>% count(TNBC)




# Perform differential expression analysis
print("Performing differential expression analysis...")
de_results <- perform_deg_analysis(expression_filtered, clinical_filtered, "TNBC")
de_results_nature <- perform_deg_analysis(expression_filtered_nature, clinical_filtered_nature, "TNBC")
# Add gene symbols using existing mapping
print("Adding gene symbols...")
de_results_with_symbols <- add_gene_symbols_to_de(de_results, gene_mapping_df)
de_results_with_symbols_nature <- add_gene_symbols_to_de(de_results_nature, gene_mapping_df)

# Print summary
print(paste("Total genes tested:", nrow(de_results_with_symbols)))
print(paste("Genes with symbols:", sum(!is.na(de_results_with_symbols$gene_symbol))))
print(paste("Significant genes (p < 0.05):", sum(de_results_with_symbols$PValue < 0.05, na.rm = TRUE)))
print(paste("Significant genes (FDR < 0.05):", sum(de_results_with_symbols$FDR < 0.05, na.rm = TRUE)))

# OPTION 1: Get top N differentially expressed genes
print("\n=== TOP 20 DIFFERENTIALLY EXPRESSED GENES ===")
top_20_genes <- get_top_de_genes(de_results_with_symbols, n = 20, criteria = "PValue")
print(top_20_genes[, c("gene_symbol", "logFC", "PValue", "FDR")])

# OPTION 2: Filter by custom gene list 
print("\n=== NEURO RECEPTOR GENES DIFFERENTIAL EXPRESSION ===")
de_neuro <- filter_de_by_genes(de_results_with_symbols, neuro_receptor_genes, filter_type = "symbol")
print("Differentially expressed neuro receptor genes:")
if (nrow(de_neuro) > 0) {
  print(de_neuro[, c("gene_symbol", "logFC", "PValue", "FDR")])
} else {
  print("No neuro receptor genes found in the differential expression results.")
  print("Available neuro receptor genes in dataset:")
  available_neuro <- intersect(neuro_receptor_genes, de_results_with_symbols$gene_symbol)
  print(available_neuro)
}

# OPTION 3: Custom gene list example (you can replace this with any gene list)
custom_genes <- c("FZD10","WNT3","WNT6","WNT5B","WNT9B","WNT11","RSPO3","DKK4", "DRAXIN", "NGF","SNAI2","SOX2","SOX17","ADAMTS5","ADAM11")
print("\n=== CUSTOM GENE LIST (Cancer-related genes) ===")
de_custom <- filter_de_by_genes(de_results_with_symbols, custom_genes, filter_type = "symbol") 
print("Differentially expressed cancer-related genes:")
if (nrow(de_custom) > 0) {
  print(de_custom[, c("gene_symbol", "logFC", "PValue", "FDR")])
}

#START Define Custom Gene-Sets----
# Define the adrenergic and cholinergic receptor genes of interest
adrenergic_genes <- c(
  # Alpha-adrenergic receptors
  "ADRA1A", "ADRA1B", "ADRA1D", "ADRA2A", "ADRA2B", "ADRA2C",
  # Beta-adrenergic receptors
  "ADRB1", "ADRB2", "ADRB3"
)

cholinergic_genes <- c(
  # Muscarinic acetylcholine receptors
  "CHRM1", "CHRM2", "CHRM3", "CHRM4", "CHRM5",
  # Nicotinic acetylcholine receptors
  "CHRNA1", "CHRNA2", "CHRNA3", "CHRNA4", "CHRNA5", "CHRNA6", "CHRNA7",
  "CHRNA9", "CHRNA10", "CHRNB1", "CHRNB2", "CHRNB3", "CHRNB4", "CHRND",
  "CHRNE", "CHRNG"
)

# Combined list of neuro receptor genes
neuro_receptor_genes <- c(adrenergic_genes, cholinergic_genes)

# Comprehensive WNT signaling pathway genes
wnt_pathway_genes <- c(
  # WNT ligands (19 human WNTs)
  "WNT1", "WNT2", "WNT2B", "WNT3", "WNT3A", "WNT4", "WNT5A", "WNT5B", 
  "WNT6", "WNT7A", "WNT7B", "WNT8A", "WNT8B", "WNT9A", "WNT9B", "WNT10A", 
  "WNT10B", "WNT11", "WNT16",
  
  # Frizzled receptors (10 FZDs)
  "FZD1", "FZD2", "FZD3", "FZD4", "FZD5", "FZD6", "FZD7", "FZD8", "FZD9", "FZD10",
  
  # LRP co-receptors
  "LRP5", "LRP6",
  
  # Beta-catenin and destruction complex
  "CTNNB1", "APC", "AXIN1", "AXIN2", "GSK3A", "GSK3B", "CK1A1", "CSNK1E",
  
  # TCF/LEF transcription factors
  "TCF7", "TCF7L1", "TCF7L2", "LEF1",
  
  # WNT target genes (key ones)
  "MYC", "CCND1", "CD44", "AXIN2", "DKK1", "NOTUM", "SP5", "ASCL2",
  
  # WNT inhibitors
  "DKK1", "DKK2", "DKK3", "DKK4", "SFRP1", "SFRP2", "SFRP4", "SFRP5", 
  "WIF1", "NOTUM", "TIKI1", "TIKI2",
  
  # R-spondins (WNT enhancers)
  "RSPO1", "RSPO2", "RSPO3", "RSPO4",
  
  # Planar cell polarity pathway
  "VANGL1", "VANGL2", "DVL1", "DVL2", "DVL3", "PRICKLE1", "PRICKLE2",
  "CELSR1", "CELSR2", "CELSR3",
  
  # WNT/Ca2+ pathway
  "CAMK2A", "CAMK2B", "PRKCA", "PRKCB", "NFATC1", "NFATC2", "NFATC3", "NFATC4",
  
  # Additional regulatory components
  "KREMEN1", "KREMEN2", "RNF43", "ZNRF3", "LGR4", "LGR5", "LGR6",
  "PORCN", "WLS", "GPR177",
  
  # Transcriptional regulators often associated with WNT
  "SOX2", "SOX17", "SNAI1", "SNAI2", "TWIST1", "TWIST2",
  
  # Matrix metalloproteinases related to WNT
  "ADAMTS5", "ADAM10", "ADAM17"
)

# Remove duplicates and sort
wnt_pathway_genes <- sort(unique(wnt_pathway_genes))
print(paste("Total WNT pathway genes:", length(wnt_pathway_genes)))


# WNT genes particularly relevant in breast cancer
breast_cancer_wnt_genes <- c(
  # Frequently dysregulated in breast cancer
  "WNT1", "WNT3", "WNT3A", "WNT5A", "WNT5B", "WNT10B", "WNT11",
  "CTNNB1", "APC", "AXIN1", "AXIN2", "GSK3B", 
  "TCF7L2", "LEF1", "MYC", "CCND1", "CD44",
  "DKK1", "DKK3", "SFRP1", "SFRP2", "SFRP4", "WIF1",
  "FZD1", "FZD2", "FZD7", "FZD8", "LRP6",
  "DVL1", "DVL3", "SNAI1", "SNAI2", "TWIST1"
)
#END   Define Custom Gene-Sets----
# VISUALIZATION -----------------------------------------------------------

# Create volcano plot highlighting neuro receptor genes
library(ggplot2)
library(ggrepel)

# Volcano plot with neuro receptor genes highlighted
volcano_neuro <- create_volcano_plot(
  de_results_with_symbols, 
  highlight_genes = neuro_receptor_genes,
  title = "TNBC vs non-TNBC: Neuro Receptor Genes Highlighted"
)
print(volcano_neuro)


# Volcano plot with top 10 genes highlighted
top_10_gene_names <- get_top_de_genes(de_results_with_symbols, n = 10)$gene_symbol
volcano_top10 <- create_volcano_plot(
  de_results_with_symbols, 
  highlight_genes = top_10_gene_names,
  title = "TNBC vs non-TNBC: Top 10 DE Genes Highlighted"
)
print(volcano_top10)

# SAVE RESULTS ------------------------------------------------------------

# Save all results
write.csv(de_results_with_symbols, "differential_expression_results_all.csv", row.names = FALSE)
write.csv(top_20_genes, "top_20_DE_genes.csv", row.names = FALSE)
write.csv(de_neuro, "neuro_receptor_DE_genes.csv", row.names = FALSE)

print("\nAnalysis complete! Results saved to CSV files.")


#Better Visualization ==========

# ENHANCED VISUALIZATION with group labels -----------------------------------

# ENHANCED VISUALIZATION with group labels (CORRECTED) ----------------------

# Print summary
print_de_summary(de_results_with_symbols, "TNBC", "non-TNBC")
print_de_summary(de_results_with_symbols_nature, "TNBC", "non-TNBC")

# Get top genes in each group
group_genes <- get_group_specific_genes(de_results_with_symbols, "TNBC", "non-TNBC", top_n = 15)

# Print top genes higher in each group
cat("=== TOP GENES HIGHER IN TNBC ===\n")
print(group_genes$higher_in_group1[, c("gene_symbol", "logFC", "PValue", "FDR")])

cat("\n=== TOP GENES HIGHER IN NON-TNBC ===\n")
print(group_genes$higher_in_group2[, c("gene_symbol", "logFC", "PValue", "FDR")])

# Create enhanced volcano plots
volcano_neuro_enhanced <- create_volcano_plot(
  de_results_with_symbols, 
  highlight_genes = neuro_receptor_genes,
  title = "Neuro Receptor Genes in TNBC vs non-TNBC",
  group1_name = "TNBC",
  group2_name = "non-TNBC"
)
print(volcano_neuro_enhanced)


# Volcano plot with Cutom (WNT Signaling) genes highlighted
volcano_custom <- create_volcano_plot(
  de_results_with_symbols, 
  highlight_genes = wnt_pathway_genes, #breast_cancer_wnt_genes, #wnt_pathway_genes,
  title = "TNBC vs non-TNBC: WNT Signaling Genes Highlighted (Published Dataset of TNBC)",
  group1_name = "TNBC",
  group2_name = "Non-TNBC"
)
print(volcano_custom)

# To Prevent multiple ENSEMLE IDs showing uphighlighted in the plot, filter one gene
de_results_filtered <- de_results_with_symbols[order(de_results_with_symbols$PValue), ]
de_results_filtered <- de_results_filtered[!duplicated(de_results_filtered$gene_symbol) | is.na(de_results_filtered$gene_symbol), ]


# Volcano plot with top genes from each group highlighted
top_genes_both_groups <- c(
  head(group_genes$higher_in_group1$gene_symbol, 25),
  head(group_genes$higher_in_group2$gene_symbol, 25)
)

volcano_top_both <- create_volcano_plot(
  de_results_filtered, #de_results_with_symbols, 
  highlight_genes = top_genes_both_groups,
  title = "Top DE Genes in Each Group",
  group1_name = "TNBC",
  group2_name = "non-TNBC"
)
print(volcano_top_both)

# Create a focused plot for significantly DE neuro receptor genes only
if (nrow(de_neuro) > 0) {
  # Filter for significant neuro genes
  sig_neuro <- de_neuro[de_neuro$PValue < 0.05, ]
  
  if (nrow(sig_neuro) > 0) {
    volcano_sig_neuro <- create_volcano_plot(
      de_results_with_symbols, 
      highlight_genes = de_neuro$gene_symbol, #sig_neuro$gene_symbol,
      title = "Significantly DE Neuro Receptor Genes",
      group1_name = "TNBC",
      group2_name = "non-TNBC",
      pval_threshold = 0.05,
      fc_threshold = 0.5  # Lower threshold to see more neuro genes
    )
    print(volcano_sig_neuro)
    
    cat("\n=== SIGNIFICANT NEURO RECEPTOR GENES ===\n")
    for (i in 1:nrow(sig_neuro)) {
      gene <- sig_neuro[i, ]
      direction <- ifelse(gene$logFC > 0, "higher in TNBC", "higher in non-TNBC")
      cat(sprintf("%s: logFC = %.2f, p = %.2e (%s)\n", 
                  gene$gene_symbol, gene$logFC, gene$PValue, direction))
    }
  }
}


## TESTING IF GENE NAME MATCHES MULTIPLE ENSEMBLE IDS
# Check how many rows match your gene list
matching_rows <- de_results_with_symbols[de_results_with_symbols$gene_symbol %in% top_genes_both_groups & 
                                           !is.na(de_results_with_symbols$gene_symbol), ]

print(paste("Number of gene symbols in top_genes_both_groups:", length(top_genes_both_groups)))
print(paste("Number of unique gene symbols in top_genes_both_groups:", length(unique(top_genes_both_groups))))
print(paste("Number of rows found in data:", nrow(matching_rows)))
print(paste("Number of unique gene symbols found:", length(unique(matching_rows$gene_symbol))))

# See which genes have multiple entries
gene_counts <- table(matching_rows$gene_symbol)
multiple_entries <- gene_counts[gene_counts > 1]
if(length(multiple_entries) > 0) {
  print("Genes with multiple Ensembl IDs:")
  print(multiple_entries)
}
