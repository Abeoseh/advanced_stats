library(edgeR)

# https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html

getwd()
setwd("C:/Users/brean/Downloads/masters/advanced_stats/FinalProject")
getwd()

# https://olvtools.com/en/documents/edger

alal <- read.csv("csvs/rna_alal.csv", check.names = F)


counts <- alal[,-c(1:4)] |> t() |> as.data.frame()

group <- factor(gsub(" ", "_", gsub(", NOS", "", alal$DISEASE_TYPE)))
dge <- DGEList(counts=counts, group=group)

## filter lowly expressed genes (0 inflated data)
keep <- filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes=FALSE]

## 3. Voom transformation and calculation of variance weights
mm <- model.matrix(~0 + group)


y <- voom(dge, mm, plot = T)

## fit the model
fit <- lmFit(y, mm)
head(coef(fit))

## comparison between lymphoid leukemia and meyloid leukemia
contr <- makeContrasts(groupLymphoid_leukemia - groupMyeloid_leukemia, levels = colnames(coef(fit)))

## estimate contrasts for each gene
tmp <- contrasts.fit(fit, contr)

## Empirical Bayes smoothing of standard errors
tmp <- eBayes(tmp)

## Write results out
top.table <- topTable(tmp, sort.by = "P", n = Inf)
top.table$higher_in <- ifelse(top.table$logFC > 0, "Myeloid leukemia, NOS", "Lymphoid leukemia, NOS")
top.table$gene <- row.names(top.table)

write.csv(top.table, "./csvs/limma-voom_results.csv", row.names = F)
