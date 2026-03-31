.libPaths( c( .libPaths(), "~/my_R_libs") )

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stringr))
# PCoA
suppressPackageStartupMessages(library(ecodist))
suppressPackageStartupMessages(library(vegan))
# PCA
suppressPackageStartupMessages(library(factoextra))
suppressPackageStartupMessages(library(lemon))
suppressPackageStartupMessages(library(ggplot2))

### IF statement explaination:
# the images are arranged in 2x2 grids

args = commandArgs(trailingOnly = TRUE)
input = args[1]
output = args[2]
pheno1 =  args[3]
pheno2 = args[4]

output = strsplit(output, ".png")[[1]]

PCOA_df <- read.csv(input)
Studies = read.csv("./csv_files/phenotypes.csv")


#### PERMANOVA ####

print("Starting PERMANOVAs")



PCOA_df = merge(PCOA_df, Studies, by.x = "Study_ID", by.y = "ID", all.x = T)
PCOA_df = relocate(PCOA_df, Author)
colnames(PCOA_df)[1:6]

IDs = unique(PCOA_df$Author)
print("IDs:")
print(IDs)

data_cols <- filter(PCOA_df, Author == IDs[1])
data_cols <- data_cols[,-c(1:6)]

dim(data_cols)
factors <- PCOA_df[,c(1,2,4)]
dim(factors)


perm_df <- data.frame()




#### PERMANOVA of Study_ID
print("starting PERMANOVA of Study ID")

data_cols <- PCOA_df[,-c(1:7)]
meta <- PCOA_df[,c(1,2,4)]


perm_df <- adonis2(data_cols ~ Study_ID, method = "bray", meta) %>% as.data.frame()


write.csv(perm_df, paste(output,"_StudyID_permanova.csv",sep=""))


print("PERMANOVA finished, starting POCAs")

############ POCAs ############
IDs = unique(PCOA_df$Study_ID)

group.colors <- c(setNames("#880808",pheno1), setNames("#333BFF", pheno2))

#### PCOA of all data combine colored by Study ID ####
perm_df <- read.csv(paste(output,"_StudyID_permanova.csv",sep=""), row.names=1, check.names=F)
perm_df <- perm_df[1,]

group.colors <- c(`Hospital: Lax et al. 2017` = "#880808", `Air Force: Sharma et al. 2019` = "#333BFF", 
                  `Dorm: Richardson et al. 2019` = "#32a848", `House: Lax et al. 2014` = "#a832a8" ) # #8a8328 burnt yellow color


rownames(PCOA_df) <- make.names(PCOA_df$Study_ID, unique = TRUE)
rownames(PCOA_df) <- gsub("^X", "", rownames(PCOA_df))

numeric_df <- PCOA_df[,7:length(PCOA_df)]
row.names(numeric_df) <- rownames(PCOA_df)

bray <- vegdist(numeric_df, method = "bray")
pcoa_val <- pco(bray, negvals = "zero", dround = 0)

pco.labels = lapply(row.names(pcoa_val$vectors), function(x) unlist(strsplit(x, split = ".", fixed=TRUE))[1]) # remove .numbers from the end of the row names

unique(pco.labels)

eigenvalues = pcoa_val$values


pcoa_val.df = data.frame(Study_ID = unlist(pco.labels),
                         PCoA1 = pcoa_val$vectors[,1],
                         PCoA2 = pcoa_val$vectors[,2])

print(colnames(pcoa_val.df))


pcoa_val.df = merge(pcoa_val.df, Studies, by.x = "Study_ID", by.y = "ID", all.x=TRUE)

print(dim(pcoa_val.df))

pco.plot = ggplot(data = pcoa_val.df, mapping = aes(x = PCoA1, y = PCoA2)) +
  geom_point(aes(col = as.factor(Author)), alpha = 0.7) +
  # scale_color_brewer(name = "phenotype", palette = "Accent") +
  stat_ellipse(level = 0.95, aes(group = Author, color = Author)) +
  scale_color_manual(values = group.colors) +
  theme(plot.margin = margin(10, 10, 20, 10), plot.caption = element_text(hjust = 0)) +
  labs(title = "PCoA of Count data", x = paste("PCo1 (", round((eigenvalues[1] / sum(eigenvalues)) * 100, 2), "%)",sep=""), 
       y = paste("PCo2 (", round((eigenvalues[2] / sum(eigenvalues)) * 100, 2),"%)",sep=""), color="Study",
       caption = paste("PERMANOVA p-value =", signif(perm_df[["Pr(>F)"]],2), "R2 =", round(perm_df$R2, 2)) )

png(paste(output,"_combine_Study_ID.png",sep=""))
print(pco.plot)
dev.off()

print("script complete")
