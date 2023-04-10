setwd("~/work/UMCG/data/SV_GWAS/v3")
abund <- read.delim("heritability_abund.txt", header = T, as.is = T, check.names = F, sep = "\t")
dsv <- read.delim("heritability_dsv.txt", header = T, as.is = T, check.names = F, sep = "\t")
vsv <- read.delim("heritability_vsv.txt", header = T, as.is = T, check.names = F, sep = "\t")

abund$lower <- abund$h2 - 1.96*abund$SE
abund$upper <- abund$h2 + 1.96*abund$SE
abund$signif <- as.factor(ifelse(abund$P < 0.05, 1, 0))
abund$type <- "abundance"

dsv$lower <- dsv$h2 - 1.96*dsv$SE
dsv$upper <- dsv$h2 + 1.96*dsv$SE
dsv$signif <- as.factor(ifelse(dsv$P < 0.05, 1, 0))
dsv$type <- "dSV"

vsv$lower <- vsv$h2 - 1.96*vsv$SE
vsv$upper <- vsv$h2 + 1.96*vsv$SE
vsv$signif <- as.factor(ifelse(vsv$P < 0.05, 1, 0))
vsv$type <- "vSV"


merged <- rbind(abund, dsv[,-1], vsv[,-1])
merged$p_bin <- cut(-merged$P, breaks = -c(0, 0.00001, 0.0001, 0.001, 0.01, 0.05, 1), labels = c(10,8,5,3,2,1))

pdf("heritability_cmp.pdf")
ggplot(merged, aes(x = SP, y = h2, color = type, size = - log10(P))) + geom_point(aes(alpha = - log10(P))) +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("Species") + labs(size = "-log10(P)")
dev.off() 

ggplot (merged, aes(x = type, y = h2)) + geom_boxplot()


sp="F.prausnitzii"
for (sp in abund$SP){
  
  lower <- abund[abund$SP == sp,"lower"]
  upper <- abund[abund$SP == sp,"upper"]

  cat (sp, nrow(dsv[dsv$SP == sp & dsv$lower > upper,]), nrow(dsv[dsv$SP == sp & dsv$upper < lower,]), "\n")
}
