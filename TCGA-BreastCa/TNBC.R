### TCGA BRCA Data Acquisition ###
# This script performs:
# 1. Data acquisition from TCGA
# 2. TNBC subtype identification
# 3. Differential expression analysis of adrenergic/cholinergic receptor genes
# 4. Survival analysis correlated with receptor expression
# 5. Visualization with volcano plots and other relevant graphics

## Install & Load Packages ----
install_and_load_packages <- function(cran_packages, bioc_packages) {
# Install missing CRAN packages
new_packages_cran <- cran_packages[!(cran_packages %in% installed.packages()[, "Package"])]
if (length(new_packages_cran) > 0) {install.packages(new_packages_cran)}
# Install missing Bioconductor packages
new_packages_bioc <- bioc_packages[!(bioc_packages %in% installed.packages()[, "Package"])]
if (length(new_packages_bioc) > 0) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {install.packages("BiocManager")}
  BiocManager::install(new_packages_bioc, update = FALSE)
}
# Load all packages
all_packages <- c(cran_packages, bioc_packages)
sapply(all_packages, require, character.only = TRUE)
}
cran_packages <- c("BiocManager", "dplyr", "ggplot2", "pheatmap", "survival", "survminer", "tidyr")
bioc_packages <- c("biomaRt", "edgeR", "EnhancedVolcano", "SummarizedExperiment", "TCGAbiolinks")
install_and_load_packages(cran_packages, bioc_packages)

## DATA ACQUISITION -----
# Define location to store TCGA Data out of DropBox
actual_data_location <- "/Users/i/Documents/Clinic3.0/TCGA/TCGA-BreastCa/GDCdata"

#Symlink GDCdata folder to store data downloaded from TCGAbiolinks outside the 
# first link is where data is, second link is where symlink is
# Create symlink if it doesn't already exist
if (!file.exists("GDCdata")) {
  system("ln -s /Users/i/Documents/Clinic3.0/TCGA/TCGA-BreastCa/GDCdata /Users/i/Dropbox/Clinic3.0/Developer/RStudio/TCGA/TCGA-BreastCa/GDCdata")
}

#Verify Symlink works
list.files("/Users/i/Dropbox/Clinic3.0/Developer/RStudio/TCGA/TCGA-BreastCa/GDCdata")

# Or simply use:
list.files("GDCdata")

#Define Target directory for TCGAbiolinks to download data
tcga_download_directory <- "/Users/i/Dropbox/Clinic3.0/Developer/RStudio/TCGA/TCGA-BreastCa/GDCdata"
#Either is fine, to download to dropbox symlink or to the target directory directly (actual_data_location)
#tcga_download_directory <- "/Users/i/Documents/Clinic3.0/TCGA/TCGA-BreastCa/GDCdata"


# Query and download TCGA breast cancer data
query_BRCA <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

# Download data (this may take time depending on your internet connection)
GDCdownload(query_BRCA, directory = tcga_download_directory) #or actual_data_location is fine too

# Prepare expression data
BRCA_data <- GDCprepare(query_BRCA, directory = tcga_download_directory) #or actual_data_location is fine too

# Extract expression and clinical data
expression_data <- assay(BRCA_data, "tpm_unstrand")
clinical_data <- colData(BRCA_data) %>% as.data.frame()

# 4. SURVIVAL ANALYSIS ----------------------------------------------------

# Prepare survival data
clinical_surv <- clinical_filtered
clinical_surv$overall_survival <- clinical_surv$days_to_death
clinical_surv$overall_survival[is.na(clinical_surv$overall_survival)] <- 
  clinical_surv$days_to_last_follow_up[is.na(clinical_surv$overall_survival)]
clinical_surv$death_event <- ifelse(is.na(clinical_surv$days_to_death), 0, 1)

# Function to perform survival analysis for a gene
analyze_gene_survival <- function(gene_id, expr_data, clinical_data) {
  # Get expression data for the gene
  if(gene_id %in% rownames(expr_data)) {
    gene_expr <- expr_data[gene_id, ]
  } else {
    return(NULL)  # Gene not found
  }
  
  # Create data frame with expression and survival data
  surv_data <- data.frame(
    patient_id = colnames(expr_data),
    gene_expression = as.numeric(gene_expr)
  )
  
  # Merge with clinical data
  surv_data <- merge(surv_data, clinical_data[, c("barcode", "overall_survival", "death_event")], 
                     by.x = "patient_id", by.y = "barcode")
  
  # Divide into high/low expression groups based on median
  surv_data$expr_group <- ifelse(surv_data$gene_expression > median(surv_data$gene_expression), 
                                 "High", "Low")
  
  # Perform log-rank test
  surv_formula <- Surv(overall_survival, death_event) ~ expr_group
  surv_fit <- survfit(surv_formula, data = surv_data)
  surv_diff <- survdiff(surv_formula, data = surv_data)
  
  # Calculate p-value
  p_value <- 1 - pchisq(surv_diff$chisq, df = 1)
  
  # Create survival plot
  surv_plot <- ggsurvplot(
    surv_fit,
    data = surv_data,
    pval = TRUE,
    risk.table = TRUE,
    title = paste0("Survival based on ", gene_id, " expression"),
    legend.labs = c("High expression", "Low expression"),
    palette = c("#E7B800", "#2E9FDF")
  )
  
  # Return results
  return(list(
    gene = gene_id,
    p_value = p_value,
    plot = surv_plot,
    data = surv_data
  ))
}

# Analyze survival for all wnt pathway genes
survival_results <- list()
for(gene in wnt_pathway_genes) {
  if(gene %in% rownames(expression_filtered)) {
    result <- analyze_gene_survival(gene, expression_filtered, clinical_filtered)
    if(!is.null(result)) {
      survival_results[[gene]] <- result
      print(paste("Survival analysis for", gene, "- p-value:", result$p_value))
    }
  }
}

# 5. VISUALIZATION --------------------------------------------------------

# Volcano plot for all genes with highlighted neuro receptor genes
volcano_plot <- EnhancedVolcano(
  de_table,
  lab = de_table$external_gene_name,
  x = 'logFC',
  y = 'PValue',
  title = 'TNBC vs non-TNBC',
  subtitle = 'Differential Expression',
  pointSize = 3.0,
  labSize = 3.5,
  colAlpha = 0.4,
  legend = c("NS", "Log2FC", "P-value", "P-value & Log2FC"),
  legendPosition = "right",
  legendLabSize = 12,
  legendIconSize = 4.0,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colConnectors = 'grey30',
  selectLab = neuro_receptor_genes[neuro_receptor_genes %in% de_table$external_gene_name]
)

# Save the volcano plot
pdf("volcano_plot_neuro_receptors.pdf", width = 12, height = 10)
print(volcano_plot)
dev.off()

# Heatmap of neuro receptor gene expression in TNBC vs non-TNBC
if(nrow(de_neuro) > 0) {
  # Extract expression data for neuro receptor genes
  neuro_expr <- expression_filtered[rownames(expression_filtered) %in% rownames(de_neuro), ]
  
  # Normalize data for heatmap
  neuro_expr_norm <- t(scale(t(neuro_expr)))
  
  # Create annotation data
  anno_col <- data.frame(
    TNBC_Status = ifelse(clinical_filtered$TNBC, "TNBC", "non-TNBC"),
    row.names = colnames(neuro_expr)
  )
  
  # Define colors
  anno_colors <- list(
    TNBC_Status = c(TNBC = "#FF0000", `non-TNBC` = "#0000FF")
  )
  
  # Generate heatmap
  heatmap <- pheatmap(
    neuro_expr_norm,
    annotation_col = anno_col,
    annotation_colors = anno_colors,
    show_colnames = FALSE,
    cluster_cols = TRUE,
    cluster_rows = TRUE,
    main = "Expression of Adrenergic/Cholinergic Receptor Genes",
    fontsize_row = 8,
    filename = "neuro_receptor_heatmap.pdf",
    width = 12,
    height = 10
  )
}

# Box plots comparing expression between TNBC and non-TNBC for significant genes
# Select top significant genes from de_neuro
if(nrow(de_neuro) > 0) {
  sig_genes <- de_neuro[de_neuro$PValue < 0.05, ]
  
  if(nrow(sig_genes) > 0) {
    top_genes <- rownames(sig_genes)[order(sig_genes$PValue)][1:min(5, nrow(sig_genes))]
    
    for(gene in top_genes) {
      expr_values <- expression_filtered[gene, ]
      plot_data <- data.frame(
        Expression = expr_values,
        Group = ifelse(clinical_filtered$TNBC, "TNBC", "non-TNBC")
      )
      
      # Create boxplot
      p <- ggplot(plot_data, aes(x = Group, y = Expression, fill = Group)) +
        geom_boxplot() +
        geom_jitter(width = 0.2, alpha = 0.5) +
        scale_fill_manual(values = c("TNBC" = "#FF0000", "non-TNBC" = "#0000FF")) +
        labs(title = paste("Expression of", gene), 
             y = "Expression (TPM)", 
             x = "") +
        theme_bw()
      
      # Save plot
      ggsave(paste0("boxplot_", gene, ".pdf"), p, width = 8, height = 6)
    }
  }
}

# Forest plot for survival analysis results
if(length(survival_results) > 0) {
  # Extract p-values and calculate hazard ratios
  forest_data <- data.frame(
    Gene = character(),
    HR = numeric(),
    lower_CI = numeric(),
    upper_CI = numeric(),
    P_Value = numeric(),
    stringsAsFactors = FALSE
  )
  
  for(gene in names(survival_results)) {
    # Get survival data
    surv_data <- survival_results[[gene]]$data
    
    # Fit Cox proportional hazards model
    cox_model <- coxph(Surv(overall_survival, death_event) ~ gene_expression, 
                       data = surv_data)
    
    # Extract hazard ratio and confidence interval
    hr <- exp(coef(cox_model))[1]
    ci <- exp(confint(cox_model))
    
    # Add to forest_data
    forest_data <- rbind(forest_data, data.frame(
      Gene = gene,
      HR = hr,
      lower_CI = ci[1],
      upper_CI = ci[3],
      P_Value = survival_results[[gene]]$p_value
    ))
  }
  
  # Order by P-value
  forest_data <- forest_data[order(forest_data$P_Value), ]
  
  # Create forest plot
  forest_plot <- ggplot(forest_data, aes(x = HR, y = Gene)) +
    geom_point(size = 3, aes(color = P_Value < 0.05)) +
    geom_errorbarh(aes(xmin = lower_CI, xmax = upper_CI, color = P_Value < 0.05), height = 0.2) +
    geom_vline(xintercept = 1, linetype = "dashed") +
    scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
    labs(title = "Hazard Ratios for Neuro Receptor Genes",
         x = "Hazard Ratio (95% CI)",
         y = "") +
    theme_bw() +
    theme(legend.position = "none")
  
  # Save forest plot
  # ggsave("forest_plot_survival.pdf", forest_plot, width = 10, height = 8)
}

# Correlation matrix between adrenergic and cholinergic genes
# Extract expression of neuro receptor genes
neuro_expr_matrix <- expression_filtered[rownames(expression_filtered) %in% neuro_receptor_genes, ]
wnt_expr_matrix <- expression_filtered[rownames(expression_filtered) %in% wnt_pathway_genes, ]

if(nrow(wnt_expr_matrix) > 1) {
  # Calculate correlation matrix
  cor_matrix <- cor(t(wnt_expr_matrix), method = "spearman")
  
  # Generate correlation heatmap
  pheatmap(
    cor_matrix,
    main = "Correlation Between Adrenergic and Cholinergic Receptor Genes",
    color = colorRampPalette(c("blue", "white", "red"))(100),
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    fontsize_row = 8,
    fontsize_col = 8,
    filename = "correlation_heatmap.pdf",
    width = 10,
    height = 10
  )
}

# Create plots comparing adrenergic vs cholinergic expression in TNBC
# Calculate mean expression for adrenergic and cholinergic genes by group
adrenergic_expr <- expression_filtered[rownames(expression_filtered) %in% adrenergic_genes, ]
cholinergic_expr <- expression_filtered[rownames(expression_filtered) %in% cholinergic_genes, ]

if(nrow(adrenergic_expr) > 0 && nrow(cholinergic_expr) > 0) {
  # Calculate mean expression for each sample
  adrenergic_mean <- colMeans(adrenergic_expr)
  cholinergic_mean <- colMeans(cholinergic_expr)
  
  # Create data frame for plotting
  plot_data <- data.frame(
    Sample = colnames(expression_filtered),
    TNBC = clinical_filtered$TNBC,
    Adrenergic = adrenergic_mean,
    Cholinergic = cholinergic_mean
  )
  
  # Convert to long format for plotting
  plot_data_long <- tidyr::pivot_longer(
    plot_data,
    cols = c(Adrenergic, Cholinergic),
    names_to = "ReceptorType",
    values_to = "Expression"
  )
  
  # Create boxplot comparing expression by group and receptor type
  comparison_plot <- ggplot(plot_data_long, aes(x = ReceptorType, y = Expression, fill = TNBC)) +
    geom_boxplot() +
    scale_fill_manual(values = c("FALSE" = "blue", "TRUE" = "red"),
                      labels = c("FALSE" = "non-TNBC", "TRUE" = "TNBC")) +
    labs(title = "Comparison of Adrenergic vs Cholinergic Receptor Expression",
         x = "Receptor Type",
         y = "Mean Expression",
         fill = "Group") +
    theme_bw()
  
  # Save plot
  ggsave("adrenergic_vs_cholinergic.pdf", comparison_plot, width = 10, height = 8)
  
  # Create scatter plot of adrenergic vs cholinergic expression
  scatter_plot <- ggplot(plot_data, aes(x = Adrenergic, y = Cholinergic, color = TNBC)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_smooth(method = "lm", se = TRUE) +
    scale_color_manual(values = c("FALSE" = "blue", "TRUE" = "red"),
                       labels = c("FALSE" = "non-TNBC", "TRUE" = "TNBC")) +
    labs(title = "Correlation between Adrenergic and Cholinergic Expression",
         x = "Mean Adrenergic Receptor Expression",
         y = "Mean Cholinergic Receptor Expression",
         color = "Group") +
    theme_bw()
  
  # Save plot
  ggsave("adrenergic_vs_cholinergic_scatter.pdf", scatter_plot, width = 10, height = 8)
}

# Save key results
# Save DE results for neuro receptor genes
write.csv(de_neuro, "de_neuro_receptors.csv")

# Save survival analysis summary
if(length(survival_results) > 0) {
  surv_summary <- data.frame(
    Gene = names(survival_results),
    P_Value = sapply(survival_results, function(x) x$p_value)
  )
  write.csv(surv_summary, "survival_analysis_summary.csv")
}

print("Analysis complete!")