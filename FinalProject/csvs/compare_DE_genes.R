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


#### compare mine and limma ####
## logFC ##
compare_logFC <- function(df1, df2, intersect_vec, logFC1, logFC2,xtitle, ytitle, title){
  
  df1 <- filter(df1, gene %in% intersect_vec) %>% dplyr::select(gene, {{logFC1}})
  plot_df <- filter(df2, gene %in% intersect_vec) %>% dplyr::select(gene, {{logFC2}}) %>% merge(df1, by = "gene")
  

  
  colors <- ifelse(plot_df[[logFC1]] > 0 & plot_df[[logFC2]] > 0, "Myeloid leukemia, NOS", 
                   ifelse(plot_df[[logFC1]] < 0 & plot_df[[logFC2]] < 0, "Lymphoid leukemia, NOS",
                          "disagreement"))
  
  # print(plot_df$colors[1:10])
  # print(length(plot_df$colors))
  # print(head(plot_df))

  logFC1 <- sym(logFC1)
  logFC2 <- sym(logFC2)  
  
  p <- ggplot(plot_df, aes(!!logFC1, !!logFC2, colour = colors)) +
    geom_point() +
    labs(x = xtitle, y = ytitle, title = title)
  
  
  png(paste("./output/", title, ".png", sep = ""))
  print(p)
  dev.off()
}


## log2FC comparison of all genes
colnames(edger)[1] = "logFC_edger" ## edgeR and limma.voom had the same colname for logFC
compare_logFC(df1 = mine, df2 = limma.voom, intersect_vec = intersect_me_limma, logFC1 = "log2FC", logFC2 = "logFC", xtitle = "limma", ytitle = "mine", 
              title = "log2FC comparison of all genes between mine and limma-voom")

compare_logFC(df1 = mine, df2 = edger, intersect_vec = intersect_me_edger, logFC1 = "log2FC", logFC2 = "logFC_edger", xtitle = "edgeR", ytitle = "mine", 
              title = "log2FC comparison of all genes between mine and edgeR")


compare_logFC(df1 = edger, df2 = limma.voom, intersect_vec = intersect_limma_edger, logFC1 = "logFC_edger", logFC2 = "logFC", xtitle = "limma", ytitle = "edgeR", 
              title = "log2FC comparison of all genes between limma-voom and edgeR")



### log2FC comparison of sig. genes
compare_logFC(df1 = filter(mine, padj < 0.05), df2 = filter(limma.voom, adj.P.Val < 0.05), intersect_vec = intersect_me_limma, logFC1 = "log2FC", logFC2 = "logFC", xtitle = "limma", ytitle = "mine", 
              title = "log2FC comparison of sig. genes between mine and limma-voom")

compare_logFC(df1 = filter(mine, padj < 0.05), df2 = filter(edger, FDR < 0.05), intersect_vec = intersect_me_edger, logFC1 = "log2FC", logFC2 = "logFC_edger", xtitle = "edgeR", ytitle = "mine", 
              title = "log2FC comparison of sig. genes between mine and edgeR")


compare_logFC(df1 = filter(edger, FDR < 0.05), df2 = filter(limma.voom, adj.P.Val < 0.05), intersect_vec = intersect_limma_edger, logFC1 = "logFC_edger", logFC2 = "logFC", xtitle = "limma", ytitle = "edgeR", 
              title = "log2FC comparison of sig. genes between limma-voom and edgeR")





#### p-value ####
compare_pvalue <- function(df1, df2, intersect_vec, pval1, pval2, xtitle, ytitle, title){
  
  df1 <- filter(df1, gene %in% intersect_vec) %>% dplyr::select(gene, {{pval1}})
  plot_df <- filter(df2, gene %in% intersect_vec) %>% dplyr::select(gene, {{pval2}}) %>% merge(df1, by = "gene")
  
  
  
  colors <- ifelse(plot_df[[pval1]] >= 0.05 & plot_df[[pval2]] >= 0.05, "Agreement, N.S", 
                       ifelse(plot_df[[pval1]] < 0.05 & plot_df[[pval2]] < 0.05, "Agreement, Sig",
                              "disagreement"))
  # print(plot_df$colors[1:10])
  # print(length(plot_df$colors))
  # print(head(plot_df))
  
  pval1 <- sym(pval1)
  pval2 <- sym(pval2)  
  
  p <- ggplot(plot_df, aes(!!pval1, !!pval2, colour = colors)) +
    geom_point() +
    labs(x = xtitle, y = ytitle, title = title)
  
  
  png(paste("./output/", title, ".png", sep = ""))
  print(p)
  dev.off()
}



compare_pvalue(df1 = mine, df2 = limma.voom, intersect_vec = intersect_me_limma, pval1 = "padj", pval2 = "adj.P.Val", xtitle = "limma", ytitle = "mine", 
              title = "p-value comparison of all genes between mine and limma-voom")

compare_pvalue(df1 = mine, df2 = edger, intersect_vec = intersect_me_edger, pval1 = "padj", pval2 = "FDR", xtitle = "edgeR", ytitle = "mine", 
              title = "p-value comparison of all genes between mine and edgeR")


compare_pvalue(df1 = edger, df2 = limma.voom, intersect_vec = intersect_limma_edger, pval1 = "FDR", pval2 = "adj.P.Val", xtitle = "limma", ytitle = "edgeR", 
              title = "p-value comparison of all genes between limma-voom and edgeR")




