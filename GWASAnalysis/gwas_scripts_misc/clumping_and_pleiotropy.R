args <- commandArgs(trailingOnly = TRUE)
library(ieugwasr)

#setwd("C://Users/Dasha/work/UMCG/data/SV_GWAS/v4")
#res_path <- "vSV_meta_all_5e-08.sorted.rsids.genes.gwas_annot.txt"
res_path <- args[1]
d <- read.delim(res_path, header = T,  sep = "\t", as.is = T, check.names = F)
svs <- unique(d$ProbeName)
res2 <- data.frame()
cnt <- 1
for (sv in svs){
  subs <- d[d$ProbeName == sv,]
  colnames(subs)[1] <- "pval"
  colnames(subs)[2] <- "rsid"
  if (cnt == 1) {
    res2 <-  ld_clump(subs, clump_r2 = 0.1, clump_kb = 250000)
    cnt <- cnt + 1
  } else {
    res2 <- rbind(res2, ld_clump(subs))
  }
}
write.table(res2, file = paste0(res_path, ".clumped_0.1.txt"), sep = "\t", quote = F, col.names = NA)

