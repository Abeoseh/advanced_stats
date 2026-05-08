library(tidyverse)
library(ggplot2)
suppressPackageStartupMessages(library(vegan))

getwd()
setwd("C:/Users/brean/Downloads/masters/advanced_stats/FinalProject")
getwd()

#### read in data ####
alal <- read.csv("csvs/rna_alal.csv", check.names = F)

numeric_alal <- alal[,5:length(alal)]
zero_cols <- names(numeric_alal)[apply(numeric_alal, 2, function(x) sum(x) == 0)]


#### PCA ####

results <- prcomp(numeric_alal[,!names(numeric_alal) %in% zero_cols], scale = T)

# results$x %>% View()

## porportion of variation explained
var_explained = summary(results)$importance[2,]



#### PERMANOVA ####
perm_df <- adonis2(numeric_alal[,!names(numeric_alal) %in% zero_cols] ~ DISEASE_TYPE, method = "euclidean", data = alal) %>% as.data.frame()



#### plot ####
png("output/PCA.png")
p <- ggplot(cbind( alal[, c("DISEASE_TYPE"),drop = FALSE], (results$x[,c("PC1", "PC2")])), aes(PC1, PC2, color = DISEASE_TYPE)) + 
  geom_point() +
  theme(plot.caption = element_text(hjust = 0)) +
  labs(x = paste("PC1 ", var_explained[1]*100, "%", sep = ""), y = paste("PC2 ", var_explained[2]*100, "%", sep = ""),
  caption = paste("R^2:", round(perm_df[1,3], 4), "p-value:", round(perm_df[1,5], 2)))
  
print(p)
dev.off()

