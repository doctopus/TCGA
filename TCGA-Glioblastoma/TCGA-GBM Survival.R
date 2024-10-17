# Install and load required packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("TCGAbiolinks", "SummarizedExperiment", "dplyr", "survival"))

library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(survival)
library(survminer)

#Symlink GDCdata folder to store data downloaded from TCGAbiolinks outside the 
# first link is where data is, second link is where symlink is going to be
system("ln -s /Users/i/Documents/Clinic3.0/TCGA/TCGA-Glioblastoma/GDCdata /Users/i/Dropbox/Clinic3.0/Developer/RStudio/TCGA/TCGA-Glioblastoma/GDCdata")
#Test if it works by adding a file to the folderr where the data will go through the symlink
list.files("/Users/i/Dropbox/Clinic3.0/Developer/RStudio/TCGA/TCGA-Glioblastoma/GDCdata")

#Define Target directory for TCGAbiolinks to download data
tcga_gbm_download_directory <- "/Users/i/Dropbox/Clinic3.0/Developer/RStudio/TCGA/TCGA-Glioblastoma/GDCdata"



# Query and download clinical data
query <- GDCquery(
  project = "TCGA-GBM",
  data.category = "Clinical",
  data.type = "Clinical Supplement",
  data.format = "BCR Biotab"
)

GDCdownload(query, directory = tcga_gbm_download_directory)

clinical_gbm <- GDCprepare(query)

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
