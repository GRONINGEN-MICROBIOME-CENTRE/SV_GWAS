args <- commandArgs(trailingOnly = TRUE)
library(ieugwasr)

res_path <- args[1]
d <- read.delim(res_path, header = T,  sep = "\t", as.is = T, check.names = F)

res2 <- data.frame()

colnames(d) <- gsub("SV", "id", colnames(d))
colnames(d) <- gsub("Pvalue", "pval", colnames(d))
res2 <-  ld_clump(d, clump_r2 = 0.1, clump_kb = 250000)

write.table(res2, file = paste0(res_path, ".clumped_0.1.txt"), sep = "\t", quote = F, row.names = F)
