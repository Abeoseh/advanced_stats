install.packages("survival")
library(survival)
suppressPackageStartupMessages(library(dplyr))

alal_meta <- read.csv("alal_target_gdc/data_clinical_patient.txt", sep="\t", comment = "#")
alal_meta <- dplyr::select(alal_meta, PATIENT_ID,	OS_STATUS, OS_MONTHS, PRIMARY_DIAGNOSIS) %>%
  dplyr::filter(PRIMARY_DIAGNOSIS == "Lymphoid leukemia, NOS" | PRIMARY_DIAGNOSIS == "Myeloid leukemia, NOS")

alal_meta$OS_STATUS <- as.numeric(sub(":.*","", alal_meta$OS_STATUS))

alal_meta <- alal_meta[complete.cases(alal_meta), ]

alal_meta$PRIMARY_DIAGNOSIS <- as.factor(alal_meta$PRIMARY_DIAGNOSIS)

# other <- read.csv("KM_Plot__Overall_(months).txt", sep = "\t")

# setdiff(alal_meta$PATIENT_ID, other$Patient.ID) |> length()

# filter(alal_meta, PATIENT_ID %in% setdiff(alal_meta$PATIENT_ID, other$Patient.ID)) %>% View()



## Surv creates survival object which is the response variable
Y = Surv(time = alal_meta$OS_MONTHS, event = alal_meta$OS_STATUS == 1)



## Stratify by primary diagnosis variable:

kmfit = survfit(Y ~ alal_meta$PRIMARY_DIAGNOSIS)

summary(kmfit, times = c(seq(0, 100, by = 10)))

## Actual Plot:
plot(kmfit, lty = c("solid", "dashed"), col = c("black", "grey"), xlab = "Survival Time In Months", ylab = "Survival Probabilities")
legend("topright", c("Lymphoid leukemia, NOS", "Myeloid leukemia, NOS"), lty = c("solid", "dashed"), col = c("black", "grey"))
