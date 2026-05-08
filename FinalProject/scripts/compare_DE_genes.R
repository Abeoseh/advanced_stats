suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))


getwd()
setwd("C:/Users/brean/Downloads/masters/advanced_stats/FinalProject")
getwd()

mine <- read.csv("csvs/my_DE.csv") #%>% filter(padj < 0.05)
limma.voom <- read.csv("csvs/limma-voom_results.csv") #%>% filter(adj.P.Val < 0.05)
edger <- read.csv("csvs/edgeR_results.csv") #%>% filter(FDR < 0.05)


intersect_me_limma <- intersect(mine$gene, limma.voom$gene)
intersect_me_edger <- intersect(mine$gene, edger$gene)
intersect_limma_edger <- intersect(edger$gene, limma.voom$gene)

all_in_common <- intersect(intersect_limma_edger, mine$gene)

length(intersect_me_limma)
length(intersect_me_edger)
length(intersect_limma_edger)
length(all_in_common)


#### compare logFC ####
compare_logFC <- function(df1, df2, intersect_vec, logFC1, logFC2,xtitle, ytitle, title){
  
  df1 <- filter(df1, gene %in% intersect_vec) %>% dplyr::select(gene, {{logFC1}})
  plot_df <- filter(df2, gene %in% intersect_vec) %>% dplyr::select(gene, {{logFC2}}) %>% merge(df1, by = "gene")
  

  
  plot_df$colors <- ifelse(plot_df[[logFC1]] > 0 & plot_df[[logFC2]] > 0, "Myeloid leukemia, NOS", 
                   ifelse(plot_df[[logFC1]] < 0 & plot_df[[logFC2]] < 0, "Lymphoid leukemia, NOS",
                          "disagreement"))
  
  group.colors <- c(`disagreement` = "#880808", `Lymphoid leukemia, NOS` = "#333BFF", 
                    `Myeloid leukemia, NOS` = "#32a848" )
  # print(plot_df$colors[1:10])
  # print(length(plot_df$colors))
  # print(head(plot_df))

  logFC1 <- sym(logFC1)
  logFC2 <- sym(logFC2)  
  
  p <- ggplot(plot_df, aes(!!logFC1, !!logFC2, colour = colors)) +
    geom_point() +
    geom_abline(intercept = 0, slope = 1, color = "red") +
    scale_color_manual(values = group.colors, name = "group") +
    labs(x = xtitle, y = ytitle, title = title)
  
  
  png(paste("./output/", title, ".png", sep = ""))
  print(p)
  dev.off()
}


## log2FC comparison of all genes
colnames(edger)[1] = "logFC_edger" ## edgeR and limma.voom had the same colname for logFC
compare_logFC(df1 = mine, df2 = limma.voom, intersect_vec = intersect_me_limma, logFC1 = "log2FC", logFC2 = "logFC", xtitle = "mine", ytitle = "limma", 
              title = "log2FC comparison of all genes between mine and limma-voom")

compare_logFC(df1 = mine, df2 = edger, intersect_vec = intersect_me_edger, logFC1 = "log2FC", logFC2 = "logFC_edger", xtitle = "mine", ytitle = "edgeR", 
              title = "log2FC comparison of all genes between mine and edgeR")


compare_logFC(df1 = edger, df2 = limma.voom, intersect_vec = intersect_limma_edger, logFC1 = "logFC_edger", logFC2 = "logFC", xtitle = "edgeR", ytitle = "limma", 
              title = "log2FC comparison of all genes between limma-voom and edgeR")



### log2FC comparison of sig. genes
compare_logFC(df1 = filter(mine, padj < 0.05), df2 = filter(limma.voom, adj.P.Val < 0.05), intersect_vec = intersect_me_limma, logFC1 = "log2FC", logFC2 = "logFC", xtitle = "mine", ytitle = "limma", 
              title = "log2FC comparison of sig. genes between mine and limma-voom")

compare_logFC(df1 = filter(mine, padj < 0.05), df2 = filter(edger, FDR < 0.05), intersect_vec = intersect_me_edger, logFC1 = "log2FC", logFC2 = "logFC_edger", xtitle = "mine", ytitle = "edgeR", 
              title = "log2FC comparison of sig. genes between mine and edgeR")


compare_logFC(df1 = filter(edger, FDR < 0.05), df2 = filter(limma.voom, adj.P.Val < 0.05), intersect_vec = intersect_limma_edger, logFC1 = "logFC_edger", logFC2 = "logFC", xtitle = "edgeR", ytitle = "limma", 
              title = "log2FC comparison of sig. genes between limma-voom and edgeR")





#### p-value ####
compare_pvalue <- function(df1, df2, logtransform = FALSE, filter_pvals=FALSE, intersect_vec, pval1, pval2, xtitle, ytitle, title){
  
  df1 <- filter(df1, gene %in% intersect_vec) %>% dplyr::select(gene, {{pval1}})
  
  if (logtransform){
    df1[[pval1]] = -log10(df1[[pval1]])
    df2[[pval2]] = -log10(df2[[pval2]])
  }

  plot_df <- filter(df2, gene %in% intersect_vec) %>% dplyr::select(gene, {{pval2}}) %>% merge(df1, by = "gene")
  plot_df <- plot_df[!is.infinite(rowSums(plot_df)),]
  
  if (filter_pvals){
    ## remove p-values 10 sd from the variance from the p-val vector with the largest range
    if (range(plot_df[[pval1]])[2] - range(plot_df[[pval1]])[1] > range(plot_df[[pval2]])[2] - range(plot_df[[pval2]])[1]){plot_df <- plot_df[plot_df[[pval1]] < var(plot_df[[pval1]]) * 10,]}
    else(plot_df <- plot_df[plot_df[[pval2]] < var(plot_df[[pval2]]) * 10,])
  }
  
  colors <- ifelse(plot_df[[pval1]] >= 0.05 & plot_df[[pval2]] >= 0.05, "Agreement, N.S", 
                       ifelse(plot_df[[pval1]] < 0.05 & plot_df[[pval2]] < 0.05, "Agreement, Sig",
                              "disagreement"))
  
  group.colors <- c(`disagreement` = "#880808", `Agreement, N.S` = "#333BFF", 
                    `Agreement, Sig` = "#32a848" )
  
  # print(plot_df$colors[1:10])
  # print(length(plot_df$colors))
  # print(head(plot_df))
  
  pval1 <- sym(pval1)
  pval2 <- sym(pval2)  
  
  p <- ggplot(plot_df, aes(!!pval1, !!pval2, colour = colors)) +
    geom_point() +
    geom_abline(intercept = 0, slope = 1, color = "red") +
    scale_color_manual(values = group.colors, name = "group") +
    labs(x = xtitle, y = ytitle, title = title)
  
  
  png(paste("./output/", title, ".png", sep = ""))
  print(p)
  dev.off()
}


# mine$padj <- -log10(mine$padj)
# limma.voom$adj.P.Val <- -log10(limma.voom$adj.P.Val)

## No transformations/manupliations 
compare_pvalue(df1 = mine, df2 = limma.voom, logtransform = FALSE, filter_pvals = FALSE, intersect_vec = intersect_me_limma, pval1 = "padj", pval2 = "adj.P.Val", xtitle = "mine", ytitle = "limma", 
               title = "p-value comparison of all genes between mine and limma-voom")


compare_pvalue(df1 = mine, df2 = edger, logtransform = FALSE, filter_pvals = FALSE, intersect_vec = intersect_me_edger, pval1 = "padj", pval2 = "FDR", xtitle = "mine", ytitle = "edgeR", 
               title = "p-value comparison of all genes between mine and edgeR")


compare_pvalue(df1 = edger, df2 = limma.voom, logtransform = FALSE, filter_pvals = FALSE, intersect_vec = intersect_limma_edger, pval1 = "FDR", pval2 = "adj.P.Val", xtitle = "edgeR", ytitle = "limma", 
               title = "p-value comparison of all genes between limma-voom and edgeR")



## log 10 ##
compare_pvalue(df1 = mine, df2 = limma.voom, logtransform = TRUE, filter_pvals = FALSE, intersect_vec = intersect_me_limma, pval1 = "padj", pval2 = "adj.P.Val", xtitle = "mine", ytitle = "limma", 
              title = "-log10 p-value comparison of all genes between mine and limma-voom")


compare_pvalue(df1 = mine, df2 = edger, logtransform = TRUE, filter_pvals = FALSE, intersect_vec = intersect_me_edger, pval1 = "padj", pval2 = "FDR", xtitle = "mine", ytitle = "edgeR", 
              title = "-log10 p-value comparison of all genes between mine and edgeR")


compare_pvalue(df1 = edger, df2 = limma.voom, logtransform = TRUE, filter_pvals = FALSE, intersect_vec = intersect_limma_edger, pval1 = "FDR", pval2 = "adj.P.Val", xtitle = "edgeR", ytitle = "limma", 
              title = "-log10 p-value comparison of all genes between limma-voom and edgeR")



## filtered ##
compare_pvalue(df1 = mine, df2 = limma.voom, logtransform = TRUE, filter_pvals = TRUE, intersect_vec = intersect_me_limma, pval1 = "padj", pval2 = "adj.P.Val", xtitle = "mine", ytitle = "limma", 
               title = "filtered -log10 p-value comparison of all genes between mine and limma-voom")


compare_pvalue(df1 = mine, df2 = edger, logtransform = TRUE, filter_pvals = TRUE, intersect_vec = intersect_me_edger, pval1 = "padj", pval2 = "FDR", xtitle = "mine", ytitle = "edgeR", 
               title = "filtered -log10 p-value comparison of all genes between mine and edgeR")


compare_pvalue(df1 = edger, df2 = limma.voom, logtransform = TRUE, filter_pvals = TRUE, intersect_vec = intersect_limma_edger, pval1 = "FDR", pval2 = "adj.P.Val", xtitle = "edgeR", ytitle = "limma", 
               title = "filtered -log10 p-value comparison of all genes between limma-voom and edgeR")



#### Histogram of p-values as a sanity check ####

png("./output/hist of my p-values.png")
hist(mine$pvalue, xlab = "p-value",
     main = "Histogram of my p-values")
dev.off()


png("./output/hist of limma-voom p-values.png")
hist(limma.voom$P.Value, xlab = "p-value",
     main = "Histogram of limma-voom p-values")
dev.off()



png("./output/hist of edgeR p-values.png")
hist(edger$PValue, xlab = "p-value",
     main = "Histogram of edgeR p-values")
dev.off()




png("./output/hist of p-values.png")

p <- ggplot() +
  geom_histogram(data = mine, aes(pvalue, fill = "Mine")) +
  geom_histogram(data = edger, aes(PValue, fill = "edgeR"), alpha = 0.7) +
  geom_histogram(data = limma.voom, aes(P.Value, fill = "limma-voom"), alpha = 0.6) +
  geom_vline(xintercept = 0.05) +
  labs(title = "Histogram of p-values", x = "p-value") + 
  scale_fill_manual(name = "Model",
                     breaks = c("Mine", "edgeR", "limma-voom"),
                     values = c("Mine" = "grey", "edgeR" = "cornflowerblue", "limma-voom" = "pink"))

print(p)



dev.off()


