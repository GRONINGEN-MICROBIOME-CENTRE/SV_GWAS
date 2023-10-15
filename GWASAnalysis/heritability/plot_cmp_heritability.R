setwd("~/work/UMCG/data/SV_GWAS/v4")
abund <- read.delim("heritability_abund_Ranko.all_sp.txt", header = T, as.is = T, check.names = F, sep = "\t")
dsv <- read.delim("heritability_dsv.txt", header = T, as.is = T, check.names = F, sep = "\t")
vsv <- read.delim("heritability_vsv.txt", header = T, as.is = T, check.names = F, sep = "\t")

#abund$lower <- abund$h2 - 1.96*abund$SE
#abund$upper <- abund$h2 + 1.96*abund$SE
abund$signif <- as.factor(ifelse(abund$P < 0.05, 1, 0))
abund$type <- "abundance"

dsv$lower <- dsv$h2 - 1.96*dsv$SE
dsv$upper <- dsv$h2 + 1.96*dsv$SE
dsv[dsv$lower < 0, "lower"] <- 0
dsv$signif <- as.factor(ifelse(dsv$P < 0.05, 1, 0))
dsv$type <- "dSV"
dsv$SV <- paste0("dSV_", dsv$SV)

vsv$lower <- vsv$h2 - 1.96*vsv$SE
vsv$upper <- vsv$h2 + 1.96*vsv$SE
vsv[vsv$lower < 0, "lower"] <- 0
vsv$signif <- as.factor(ifelse(vsv$P < 0.05, 1, 0))
vsv$type <- "vSV"
vsv$SV <- paste0("vSV_", vsv$SV)

# merged <- rbind(abund, dsv[,-1], vsv[,-1])
# merged[merged$P == 0, "P"] <- 1e-4
# merged$p_bin <- cut(-merged$P, breaks = -c(-0, 0.00001, 0.0001, 0.001, 0.01, 0.05, 1), labels = c(10,8,5,3,2,1))
# 
# 
# species_heritable <- abund[abund$signif == 1, "SP"]
# pdf("heritability_cmp_v2.pdf")
# ggplot(merged[merged$signif == 1 & merged$SP %in% species_heritable,], aes(x = SP, y = h2, color = type, size = - log10(P))) + geom_point(aes(alpha = - log10(P))) +
#   theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("Species") + labs(size = "-log10(P)")
# dev.off() 
# 
# ggplot (merged, aes(x = type, y = h2)) + geom_boxplot()


merged_sv <- rbind(dsv, vsv)
abund2 <- abund
row.names(abund2) <- abund2$SP
merged_sv$abundance_h2 <- abund2[merged_sv$SP, "h2"]
merged_sv$abundance_h2_lower <- abund2[merged_sv$SP, "lower"]
merged_sv$abundance_h2_upper <- abund2[merged_sv$SP, "upper"]

merged_sv <- merged_sv[order(merged_sv$h2, decreasing = T),]

# ggplot(merged_sv[merged_sv$signif == 1,], aes (x = reorder(SV, -h2), y = h2, color = type)) + 
#   geom_bar(fill = "white", position = "dodge", stat = "identity" ) + 
#   geom_errorbar(color = "grey", aes(ymin=lower, ymax=upper), width=.2,
#                 position="dodge") +
#   theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "none") +
#   facet_grid(.~SP, scale="free_x", space="free") + 
#   xlab("Species") + ylab ("heritability")

pdf("heritability_cmp_v4.pdf", width = 6, height = 9)
colors = c("#3abdaf", "#345d88")
ggplot(merged_sv[merged_sv$signif == 1,], aes (x = reorder(SV, h2), y = h2, fill = type)) + 
  geom_bar( position = "dodge", stat = "identity" , width=.7) + 
  geom_errorbar(color = "#727271", aes(ymin=lower, ymax=upper), width=.5, linewidth = 0.3,
                position="dodge") +
  geom_hline(aes(yintercept = abundance_h2), color = "red") +
  geom_rect(aes(ymin = abundance_h2_lower, ymax = abundance_h2_upper, xmin = -Inf, xmax = Inf), fill = "red", alpha = 0.07) +
  theme(axis.text.y = element_blank(), axis.ticks.y=element_blank(), 
         legend.position = "none", strip.text.y.left = element_text(angle = 0),
         panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6), panel.spacing = unit(.5, "lines")) +
  xlab("Species") + ylab ("heritability") + coord_flip() +
  facet_grid(SP~., scale="free", space="free",switch = "y") +
  scale_fill_manual(values = colors) 
dev.off()
