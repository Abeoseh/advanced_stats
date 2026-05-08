# install.packages("BiocManager")
# BiocManager::install("org.Hs.eg.db")
# BiocManager::install("annotate")
library(org.Hs.eg.db)
library(annotate)
library(ggplot2)
library(dplyr)


getwd()
setwd("C:/Users/brean/Downloads/masters/advanced_stats/FinalProject")
getwd()


mine <- read.csv("csvs/my_DE.csv") %>% filter(padj < 0.05)
limma.voom <- read.csv("csvs/limma-voom_results.csv") %>% filter(adj.P.Val < 0.05)
edger <- read.csv("csvs/edgeR_results.csv") %>% filter(FDR < 0.05)

colnames(edger)[1] = "edgeR_log2FC"
colnames(limma.voom)[1] = "limma_log2FC"







#### plot top and bottom 10 log2 FC ####

log2FC_plot <- function(df, lof2FC_column, group1, group2, title){
  
  lof2FC_column_var <- sym(lof2FC_column)
  
  plot_df <- rbind({{df}} |> arrange(desc(!!lof2FC_column_var)) |> slice(1:10), {{df}} |> arrange((!!lof2FC_column_var)) |> slice(1:10)) #%>% View()
  
  wilcox_plot = ggplot(
    rbind({{df}} |> arrange(desc(!!lof2FC_column_var)) |> slice(1:10), {{df}} |> arrange((!!lof2FC_column_var)) |> slice(1:10)), 
    aes(.data[[lof2FC_column]], y = reorder(gene_names, -.data[[lof2FC_column]]), color = higher_in)) + 
  geom_point() +
  scale_color_hue(name = "Disease type", labels = c(group1, group2)) +
  geom_vline(xintercept = 0) +
  labs(title = title, x = "log2FC", y = "Gene")

  png(paste("./output/", title, ".png", sep=""), width = 700)
  print(wilcox_plot)
  dev.off()
  
  return(plot_df)
  
}

## Mine
mine$gene_names <- getSYMBOL(as.character(mine$gene), data='org.Hs.eg')
my_top_bottom_10 <- log2FC_plot(mine, lof2FC_column = "log2FC", group1 = "Lymphoid leukemia, NOS", group2 = "Myeloid leukemia, NOS",
            title = "Top and bottom 10 log2FC for Mine")


## limma
limma.voom$gene_names <- getSYMBOL(as.character(limma.voom$gene), data='org.Hs.eg')
limma_top_bottom_10 <- log2FC_plot(limma.voom, lof2FC_column = "limma_log2FC", group1 = "Lymphoid leukemia, NOS", group2 = "Myeloid leukemia, NOS",
            title = "Top and bottom 10 log2FC for limma-voom")


## edgeR
edger$gene_names <- getSYMBOL(as.character(edger$gene), data='org.Hs.eg')
edger_top_bottom_10 <- log2FC_plot(edger, lof2FC_column = "edgeR_log2FC", group1 = "Lymphoid leukemia, NOS", group2 = "Myeloid leukemia, NOS",
            title = "Top and bottom 10 log2FC for edgeR")



intersect(edger_top_bottom_10$gene_names, my_top_bottom_10$gene_names)
# [1] "OR5P4P"
intersect(edger_top_bottom_10$gene_names, limma_top_bottom_10$gene_names)
# [1] "MSLN"
intersect(my_top_bottom_10$gene_names, limma_top_bottom_10$gene_names)
# character(0)




## All 3 combine
plot_df <- merge(mine, edger, by = "gene") %>% 
  merge(select(limma.voom, -c(AveExpr, t, B)), by = "gene")

plot_df$gene_names <- getSYMBOL(as.character(plot_df$gene), data='org.Hs.eg')

## normalize log2FC, (log2FC_i - u)/sd ... mean = 0, sd = 1
plot_df <- plot_df %>% mutate(mine_std_log2fc = ( (log2FC - mean(log2FC)) / sd(log2FC) ),
                   edgeR_std_log2FC = ( (edgeR_log2FC - mean(edgeR_log2FC)) / sd(edgeR_log2FC) ),
                   limma_std_log2FC = ( (limma_log2FC - mean(limma_log2FC)) / sd(limma_log2FC) ),
                   combine_std_log2FC = (mine_std_log2fc + edgeR_std_log2FC + limma_std_log2FC)) 
                  


## 
common_top_bottom_mine <- log2FC_plot(plot_df, lof2FC_column = "log2FC", group1 = "Lymphoid leukemia, NOS", group2 = "Myeloid leukemia, NOS",
                                   title = "Common top and bottom 10 log2FC for Mine")

common_top_bottom_edgeR <- log2FC_plot(plot_df, lof2FC_column = "edgeR_log2FC", group1 = "Lymphoid leukemia, NOS", group2 = "Myeloid leukemia, NOS",
                           title = "Common top and bottom 10 log2FC for edgeR")

common_top_bottom_limma <- log2FC_plot(plot_df, lof2FC_column = "limma_log2FC", group1 = "Lymphoid leukemia, NOS", group2 = "Myeloid leukemia, NOS",
                           title = "Common top and bottom 10 log2FC for limma-voom")

throwawaydf <- log2FC_plot(plot_df, lof2FC_column = "combine_std_log2FC", group1 = "Lymphoid leukemia, NOS", group2 = "Myeloid leukemia, NOS",
                           title = "Top and bottom 10 log2FC for all three")


intersect(common_top_bottom_mine$gene_names, common_top_bottom_edgeR$gene_names)
#  [1] "IGSF1"     "MSLN"      "HMX3"      "PDCD6IPP1" "LINC02600" "ADGRG6"    "ROBO1"    
# [8] "HIF3A"     "DES"       "LINC01749" "MIR663AHG" "CNN2P6"    "H3P11"     "B4GALNT4" 
# [15] "BHLHE23"   "MIR466"    "MYT1L"  
intersect(common_top_bottom_edgeR$gene_names, common_top_bottom_limma$gene_names)
# [1] "MSLN"      "MIR663AHG" "CNN2P6"    "FFAR1" 
intersect(common_top_bottom_mine$gene_names, common_top_bottom_limma$gene_names)
# "MSLN" "MIR663AHG" "CNN2P6" 




