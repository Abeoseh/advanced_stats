suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyverse))

getwd()
setwd("C:/Users/brean/Downloads/masters/advanced_stats/FinalProject")
getwd()

alal_rna <- read.csv("alal_target_gdc/data_mrna_seq_read_counts.txt", sep="\t", comment = "#", check.names = FALSE)
alal_meta <- read.csv("alal_target_gdc/data_clinical_sample.txt", sep="\t", comment = "#")


# aml_rna <- read.csv("aml_tcga_gdc/data_mrna_seq_read_counts.txt", sep="\t", comment = "#", check.names = FALSE)
# aml_meta <- read.csv("aml_tcga_gdc/data_clinical_sample.txt", sep="\t", comment = "#")



#### merge metadata with samples ####

## alal ##
t_alal_rna <- as.data.frame(t(alal_rna))
colnames(t_alal_rna) <- t_alal_rna[1,]
t_alal_rna <- t_alal_rna[-c(1),]
t_alal_rna$SAMPLE_ID <- row.names(t_alal_rna)

alal <- alal_meta %>% filter(DISEASE_TYPE != "Leukemia, NOS") %>%
  select(SAMPLE_ID, PATIENT_ID, DISEASE_TYPE) %>% 
  merge(t_alal_rna, by = "SAMPLE_ID")

# alal$Study_ID = "alal"
# 
## aml ##
# t_aml_rna <- as.data.frame(t(aml_rna))
# colnames(t_aml_rna) <- t_aml_rna[1,]
# t_aml_rna <- t_aml_rna[-c(1),]
# t_aml_rna$SAMPLE_ID <- row.names(t_aml_rna)
# 
# colnames(aml_meta)[8] = "DISEASE_TYPE"
# 
# aml <- aml_meta %>% select(SAMPLE_ID, PATIENT_ID, DISEASE_TYPE) %>%
#   merge(t_aml_rna, by = "SAMPLE_ID")
# 
# aml$Study_ID = "aml"


## both merged ##
# FALSE %in% (colnames(aml) == colnames(alal)) ## they both have the same columns so it prints FALSE
# aml_alal <- rbind(aml, alal)
# 
# aml_alal <- relocate(aml_alal, "Study_ID")
# aml_alal$DISEASE_TYPE = gsub("Acute Myeloid Leukemia", "Myeloid leukemia" , gsub(", NOS", "", aml_alal$DISEASE_TYPE))
# colnames(aml_alal)[1:10]

write.csv(alal, "csvs/rna_aml_alal.csv", row.names = F)



