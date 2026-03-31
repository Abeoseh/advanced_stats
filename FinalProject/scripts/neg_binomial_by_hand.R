library("gap")
library(MASS)

getwd()
setwd("C:/Users/brean/Downloads/masters/advanced_stats/FinalProject")
getwd()

alal <- read.csv("csvs/rna_alal.csv", check.names = F)
alal$DISEASE_TYPE <- factor(alal$DISEASE_TYPE)
countsTable <- alal[,-c(1:4)]
results <- list()

geometricMean <- geometricMean <- apply(countsTable, 2, function(x) {
  if(all(x == 0)) return(NA)
  exp(mean(log(x[x > 0])))
})

ratios <- sweep(countsTable, 2, geometricMean, "/")
size_factors <- apply(ratios, 1, median, na.rm = TRUE) ## median ratio

for( gene in colnames(countsTable) ){
  
  y <- countsTable[, gene]
  
  ## estimate dispersion
  mu_hat <- mean(y)
  if(mu_hat < 1e-6) next
  
  var_hat <- var(y)
  phi <- max((var_hat - mu_hat) / mu_hat^2, 1e-3) ## estimate dispersion... using 1e-8 since 0 later I do 1/phi and I can't do 1/0. Even though I want max(x,0) I use 1e-3 since 1/1e-3 is only 1k which is reasonable
  ## fit full model while considering disease type (DISEASE_TYPE)
  beta_full <- try(glm.nb(y ~ alal$DISEASE_TYPE + offset(log(size_factors)), init.theta = 1/phi,
    link = log), silent = TRUE)
  if(inherits(beta_full, "try-error")) next
  
  log_likelihood_full <- logLik(beta_full)
  
  ## fit reduced model... or model with no condition ... in this case there is only one group
  ## so this is the null that there is no difference between groups
  beta_null <- try(glm.nb(y ~ 1 + offset(log(size_factors)), init.theta = 1/phi,
    link = log), silent = TRUE)
  if(inherits(beta_null, "try-error")) next
  
  log_likelihood_null <- logLik(beta_null)
  
  ## likelihood ratio test
  stat <- 2 * (log_likelihood_full - log_likelihood_null)
  pval <- pchisq(stat, df=1, lower.tail=FALSE)
    
  results[[gene]] <- list(
    log2FC = beta_full$coefficients[2] / log(2),
    pvalue = pval[1]
  )
}


## multiple testing correction
pvals <- sapply(results, function(x) x$pvalue)
padj <- p.adjust(pvals, method="BH")

names(padj) <- names(results)

res_df <- data.frame(
  gene = names(results),
  log2FC = sapply(results, function(x) x$log2FC),
  pvalue = sapply(results, function(x) x$pvalue),
  padj = padj
)

res_df$higher_in <- ifelse(res_df$log2FC > 0, "Myeloid leukemia, NOS", "Lymphoid leukemia, NOS")

write.csv(res_df, "csvs/my_DE.csv", row.names = F)
