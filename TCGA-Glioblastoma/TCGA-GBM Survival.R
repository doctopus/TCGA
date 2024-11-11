# Install and load required packages
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(c("TCGAbiolinks", "SummarizedExperiment", "dplyr", "survival"))

library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(survival)
library(survminer)
library(mclust)


#Symlink GDCdata folder to store data downloaded from TCGAbiolinks outside the 
# first link is where data is, second link is where symlink is going to be
system("ln -s /Users/i/Documents/Clinic3.0/TCGA/TCGA-Glioblastoma/GDCdata /Users/i/Dropbox/Clinic3.0/Developer/RStudio/TCGA/TCGA-Glioblastoma/GDCdata")
#Test if it works by adding a file to the folderr where the data will go through the symlink
list.files("/Users/i/Dropbox/Clinic3.0/Developer/RStudio/TCGA/TCGA-Glioblastoma/GDCdata")

#Define Target directory for TCGAbiolinks to download data
tcga_gbm_download_directory <- "/Users/i/Dropbox/Clinic3.0/Developer/RStudio/TCGA/TCGA-Glioblastoma/GDCdata"

#Get TCGA Project Summary
TCGAbiolinks:::getProjectSummary('TCGA-GBM')

#-------------↓
# Query and download clinical data
query_clinical <- GDCquery(
  project = "TCGA-GBM",
  data.category = "Clinical",
  data.type = "Clinical Supplement",
  data.format = "BCR Biotab"
)
GDCdownload(query_clinical, directory = tcga_gbm_download_directory)
clinical_gbm <- GDCprepare(query_clinical)



#-------------↑

#-------------↓
# Query and download Transcriptome Data
query_transcriptome <- GDCquery(
  project = "TCGA-GBM",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

GDCdownload(query_transcriptome, directory = tcga_gbm_download_directory)
transcriptome_gbm <- GDCprepare(query_transcriptome)
#Subset the SummarizedExperiment object to keep only IDH WT
# Subset the SummarizedExperiment object
transcriptome_gbm_subset <- transcriptome_gbm[, colData(transcriptome_gbm)$paper_IDH.status == "WT" & !is.na(colData(transcriptome_gbm)$paper_IDH.status)]
# Extract the barcodes from transcriptome_gbm_subset

patients_to_keep <- colData(transcriptome_gbm_subset)$patient
# Subset clinical_data to keep only rows with matching barcodes
clinical_data_subset <- survival_data[survival_data$bcr_patient_barcode %in% patients_to_keep, ]

#Plot Start - Survival Plot ----
survival_data <- clinical_data_subset
# Create a new column for survival time
# survival_data$survival_time <- ifelse(survival_data$vital_status == "Dead", survival_data$days_to_death, survival_data$days_to_last_follow_up)

# Create a new column for censoring status
# survival_data$censor <- ifelse(survival_data$vital_status == "Alive", 0, 1)

# Convert survival time and censor to numeric
# survival_data$survival_time <- as.numeric(survival_data$survival_time)
# survival_data$censor <- as.numeric(survival_data$censor)

# Create a new column for gender
# survival_data$gender <- substr(survival_data$bcr_patient_barcode, 14, 14)

# Plot Kaplan-Meier curve comparing survival between male and female patients
# fit <- survfit(Surv(survival_time, censor) ~ gender, data = survival_data)
# ggsurvplot(fit, data = survival_data, pval = TRUE)

survival_data <- clinical_data_subset %>%
  mutate(
    survival_time = as.numeric(ifelse(vital_status == "Alive", days_to_last_follow_up, days_to_death)),
    survival_time = survival_time / 30.417,  # Divide by 365.25 to Convert days to years OR divide by 30.417 to convet to Month
    prognosis = ifelse(survival_time > 18, "Good", "Poor")
  )
survival_data <- survival_data %>% 
  mutate_at(c('age_at_initial_pathologic_diagnosis', 'death_days_to', 'last_contact_days_to'), as.numeric) %>% 
  mutate_at('prognosis', as.factor)

# Create a new column for censoring status
survival_data$censor <- ifelse(survival_data$vital_status == "Alive", 0, 1)

fit <- survfit(Surv(survival_time, censor) ~ gender, data = survival_data)
ggsurvplot(fit, data = survival_data, title="TCGA-GBM Survival of Male & Female",
           legend =c(0.8, 0.8), palette = c("#FC6BF4", "#007DEF"), legend.title="Sex",
           risk.table =TRUE, risk.table.col="gender",
           legend.labs=c("Female", "Male"),
           ggtheme=theme_survminer() + theme(plot.title=element_text(hjust = 0.5, face = "plain" )),
           xlab="Time in Months", size=1,
           break.x.by=12,
           pval = FALSE, conf.int = FALSE,
           surv.median.line = "v",
           surv.scale="percent",
           axes.offset=TRUE)

ggsurvplot(fit, data = survival_data, title="TCGA-GBM Survival of Male & Female",           
           facet.by = c("prognosis"))
#Plot End -Survival Plot ----
#-------------↑

#-------------↓
# Query and download biospecimen data
query_biospecimen <- GDCquery(
  project = "TCGA-GBM",
  data.category = "Biospecimen",
  data.type = "Biospecimen Supplement",
  data.format = "BCR Biotab"
)
GDCdownload(query_biospecimen, directory = tcga_gbm_download_directory)
biospecimen_gbm <- GDCprepare(query_biospecimen)

#-------------↑
# Extract relevant clinical information
clinical_gbm_data <- clinical_gbm$clinical_patient_gbm %>%
  slice(3:n()) %>%
  select(bcr_patient_barcode, age_at_initial_pathologic_diagnosis, 
         death_days_to, last_contact_days_to, vital_status, gender,radiation_treatment_adjuvant,
         tumor_status)

# Calculate survival time and create prognosis groups
clinical_gbm_data <- clinical_gbm_data %>%
  mutate(
    survival_time = as.numeric(ifelse(vital_status == "Alive", last_contact_days_to, death_days_to)),
    survival_time = survival_time / 30.417,  # Divide by 365.25 to Convert days to years OR divide by 30.417 to convet to Month
    prognosis = ifelse(survival_time > 18, "Good", "Poor")
  )
clinical_gbm_data <- clinical_gbm_data %>% 
  mutate_at(c('age_at_initial_pathologic_diagnosis', 'death_days_to', 'last_contact_days_to'), as.numeric) %>% 
  mutate_at('prognosis', as.factor)

# Create a new column for censoring status
clinical_gbm_data$censor <- ifelse(clinical_gbm_data$vital_status == "Alive", 0, 1)



# Plot Kaplan-Meier curve comparing survival between male and female patients
fit <- survfit(Surv(survival_time, censor) ~ gender, data = clinical_gbm_data)
ggsurvplot(fit, data = clinical_gbm_data, title="TCGA-GBM Survival of Male & Female",
           legend =c(0.8, 0.8), palette = c("#FC6BF4", "#007DEF"), legend.title="Sex",
           risk.table =TRUE, risk.table.col="gender",
           legend.labs=c("Female", "Male"),
           ggtheme=theme_survminer() + theme(plot.title=element_text(hjust = 0.5, face = "plain" )),
           xlab="Time in Months", size=1,
           break.x.by=12,
           pval = FALSE, conf.int = FALSE,
           surv.median.line = "v",
           surv.scale="percent",
           axes.offset=TRUE)

ggsurvplot(fit, data = clinical_gbm_data, title="TCGA-GBM Survival of Male & Female",           
           facet.by = c("prognosis"))


###-----------------Gaussian Mixure Modelling via mclust packaghe, to identify a threshold for IDH1 mutant
#IDH1 mutation should typically determined by DNA sequencing or immunohistochemistry. Inferring from gene expression levels is a simplification.
# Gaussian mixture model assumes a bimodal distribution of IDH1 expression, which may not always be the case. Need to visually inspect the resulting plot to ensure this assumption holds.
# This approach may not be as accurate as direct genetic testing for IDH1 mutations. It is a proxy method and need to be validated against known mutation data if possible.
# The Ensembl ID for IDH1 (ENSG00000138413) is used to extract its expression data.
# Query and download gene expression data
query_expression <- GDCquery(
  project = "TCGA-GBM", 
  data.category = "Transcriptome Profiling", 
  data.type = "Gene Expression Quantification", 
  workflow.type = "STAR - Counts"
)
GDCdownload(query_expression)
expression <- GDCprepare(query_expression)

# Extract IDH1 expression data
idh1_expression <- assay(expression)["ENSG00000138413",]  # IDH1 Ensembl ID
idh1_data <- data.frame(
  bcr_patient_barcode = colData(expression)$barcode,
  idh1_expression = idh1_expression
)

# Fit Gaussian mixture model
gmm_fit <- Mclust(idh1_data$idh1_expression, G=2)

# Determine threshold
threshold <- mean(gmm_fit$parameters$mean)

# Classify samples
idh1_data$idh1_status <- ifelse(idh1_data$idh1_expression > threshold, "Mutant", "Wild-type")

# Merge clinical and IDH1 status data
clinical_data <- merge(clinical_data, idh1_data, by="bcr_patient_barcode")

# Summarize clinical data
summary_age <- summary(clinical_data$age_at_initial_pathologic_diagnosis)
summary_stage <- table(clinical_data$tumor_stage)
summary_idh1 <- table(clinical_data$idh1_status)

print(summary_age)
print(summary_stage)
print(summary_idh1)

# Display the first few rows of the processed clinical data
head(clinical_data)

# Plot IDH1 expression distribution
library(ggplot2)
ggplot(clinical_data, aes(x=idh1_expression, fill=idh1_status)) +
  geom_density(alpha=0.5) +
  geom_vline(xintercept=threshold, linetype="dashed") +
  labs(title="Distribution of IDH1 Expression", x="IDH1 Expression", y="Density") +
  theme_minimal()