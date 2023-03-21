setwd("/Users/Dasha/work/UMCG/data/SV_GWAS/v3")
svtype="vSV"

sv = "B.obeum:22"
snp = "20:60298955"

c="500FG"
print(c)
pheno <- read.delim(paste0(c, ".", svtype, ".filtered.txt"), header = T, sep = "\t", as.is = T, check.names = F, row.names = 2)
covar <- read.delim(paste0( c, ".covariates.txt"), header = T, sep = "\t", as.is = T, check.names = F, row.names = 2)
if (svtype == "dSV") pheno[,sv] <- as.factor(pheno[,sv])
geno <- as.data.frame(t(read.delim(paste0(c, ".", snp, ".genotypes.txt"), header = T, sep = "\t", as.is = T, check.names = F, row.names = 1)))
geno[,1] <- as.factor(geno[,1])

sex <- read.table(paste0( c , "_gender.txt"), header = F,  as.is = T, check.names = F, row.names = 2)

ids <- intersect(row.names(pheno), row.names(covar))
m <- cbind(pheno[ids,sv], covar[ids, c(gsub(":[0-9]+","", sv), "age", "read_number")], sex[ids,2], geno[ids,])
colnames(m) <- c("sv", "abundance", "age", "read_number", "sex", "geno_factor")
m <- na.omit(m)

covar_formula <- as.formula("sv ~ abundance + age + read_number + sex")
fit_covar <- lm(covar_formula, data = m)
m$sv_corrected <- residuals(fit_covar)


colors = c("#225EA8", "#CC4C02")
ggplot(m, aes(x = geno_factor, y = sv_corrected)) +
  geom_violin(trim=FALSE, color = colors[1]) +
  geom_boxplot(width=0.5, color = colors[1]) +
  labs(title=paste0(sv, " : ", snp), x = "genotype", y = paste0(svtype, " corrected abundance")) +
  theme_classic() + theme(legend.position = "none") +
  scale_color_manual(values = c(colors[1]))

