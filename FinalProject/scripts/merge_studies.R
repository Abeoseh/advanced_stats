suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyverse))

getwd()
setwd("C:/Users/brean/Downloads/masters/advanced_stats/FinalProject")
getwd()

alal_rna <- read.csv("alal_target_gdc/data_mrna_seq_read_counts.txt", sep="\t", comment = "#", check.names = FALSE)
alal_meta <- read.csv("alal_target_gdc/data_clinical_sample.txt", sep="\t", comment = "#")


# laml_rna <- read.csv("laml_tcga/data_mrna_seq_read_counts.txt", sep="\t", comment = "#", check.names = FALSE)
# laml_meta <- read.csv("laml_tcga/data_clinical_sample.txt", sep="\t", comment = "#")



#### merge metadata with samples ####

## alal ##
t_alal_rna <- as.data.frame(t(alal_rna))
colnames(t_alal_rna) <- t_alal_rna[1,]
t_alal_rna <- t_alal_rna[-c(1),]
t_alal_rna$SAMPLE_ID <- row.names(t_alal_rna)

alal <- alal_meta %>% filter(DISEASE_TYPE != "Leukemia, NOS") %>%
  select(SAMPLE_ID, PATIENT_ID, DISEASE_TYPE) %>% 
  merge(t_alal_rna, by = "SAMPLE_ID")

alal$Study_ID = "alal"

## laml ##
# t_laml_rna <- as.data.frame(t(laml_rna))
# colnames(t_laml_rna) <- t_laml_rna[1,]
# t_laml_rna <- t_laml_rna[-c(1),]
# t_laml_rna$SAMPLE_ID <- row.names(t_laml_rna)
# 
# laml <- laml_meta %>% select(SAMPLE_ID, PATIENT_ID, DISEASE_TYPE) %>%
#   merge(t_laml_rna, by = "SAMPLE_ID")
# 
# laml$Study_ID = "laml"


## both merged ##
# FALSE %in% (colnames(laml) == colnames(alal)) ## they both have the same columns
# laml_alal <- cbind(laml, alal)

alal <- relocate(alal, "Study_ID")
colnames(alal)[1:10]

write.csv(alal, "csvs/rna_alal.csv", row.names = F)

#### Make Kaplan Meier Dataset ####
## Needed columns for data_clinical_patient.txt:
## Study ID,	Patient ID,	OS_STATUS,	OS_MONTHS


