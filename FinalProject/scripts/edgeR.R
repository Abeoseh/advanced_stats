# BiocManager::install("edgeR")
library(edgeR)

getwd()
setwd("C:/Users/brean/Downloads/masters/advanced_stats/FinalProject")
getwd()

# https://olvtools.com/en/documents/edger


## read in the data
alal <- read.csv("csvs/rna_alal.csv", check.names = F)


counts <- alal[,-c(1:4)] |> t() |> as.data.frame()
group <- factor(alal$DISEASE_TYPE)
dge <- DGEList(counts=counts, group=group)

## filter lowly expressed genes (0 inflated data)
keep <- filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes=FALSE]


## Estimate dispersion
dge <- estimateDisp(dge, design=model.matrix(~group))

## Fit GLM
fit <- glmQLFit(dge, design = model.matrix(~group))
qlf <- glmQLFTest(fit, coef=2)

## saving results
topTags(qlf)

result <- as.data.frame(topTags(qlf, n=nrow(dge)))
result$higher_in <- ifelse(result$logFC > 0, "Myeloid leukemia, NOS", "Lymphoid leukemia, NOS")

result$gene <- row.names(result)

write.csv(result, "./csvs/edgeR_results.csv", row.names = F)
