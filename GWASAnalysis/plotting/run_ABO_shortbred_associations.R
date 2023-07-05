# ml RPlus/4.0.3-foss-2018c-v21.12.10
library(kinship2, lib.loc="/groups/umcg-lld/tmp01/other-users/umcg-dzhernakova/R/x86_64-pc-linux-gnu-library/4.0.3")
library(lme4qtl, lib.loc="/groups/umcg-lld/tmp01/other-users/umcg-dzhernakova/R/x86_64-pc-linux-gnu-library/4.0.3")
library(lme4)
library(ggplot2)
library(patchwork)
library('ggpubr')
library(cowplot)

fut2_snp <- "19:49206674"
recode_genotypes <- function(geno, a1, a2){
  geno$ac <- geno[,1]
  geno$ac <- sub(paste0(a2, "/", a2), 0, geno$ac)
  geno$ac <- sub(paste0(a2, "/", a1), 1, geno$ac)
  geno$ac <- sub(paste0(a1, "/", a2), 1, geno$ac)
  geno$ac <- sub(paste0(a1, "/", a1), 2, geno$ac)
  geno$ac <- as.numeric(geno$ac)
  geno[,1] <- as.factor(geno[,1])
  return(geno)
}

get_fut2_secretor_status <- function(fut2_geno){
  fut2_geno$secr <- "secretor"
  fut2_geno$secr[fut2_geno[,"19:49206674"] == "A/A"] <- "non-secretor"
  fut2_geno$secr <- as.factor(fut2_geno$secr)
  fut2_geno$secr <- as.factor(fut2_geno$secr)
  return (fut2_geno)
}

format_abo_bloodgroup <- function(abo){
  abo$ABO_O <- "A/B"
  abo[abo$Bloodtype == "O", "ABO_O"] <- "O"
  abo$ABO_O <- as.factor(abo$ABO_O)
  abo$Bloodtype <- as.factor(abo$Bloodtype)
  abo$ABO_A <- "B/O"
  abo[abo$Bloodtype %in% c("AB", "A"), "ABO_A"] <- "A/AB"
  abo$ABO_A <- as.factor(abo$ABO_A)
  abo$A_geno <- 0
  abo[abo$Blood_genotype %in% c("A", "AB"), "A_geno"] <- 1
  abo[abo$Blood_genotype == "AA", "A_geno"] <- 2
  abo$ABO_blood_group <- abo$Bloodtype
  return(abo)
}


rename_genes <- function(pheno){
  id_conv <- read.delim("/groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/umcg-dzhernakova/SV_GWAS/results/shortbred/id_convertion.txt", header = T, sep = "\t", as.is = T, check.names = F, row.names = 1)
  colnames(pheno) <- id_conv[colnames(pheno),"new_id"]
  return(pheno)
}

fit_lmm <- function(dat, abo_group){
  covars = "age + read_number + sex + abundance + sv1 +sv2 + sv3 + sv4 + sv5 + (1 | ID)"
  
  lmm <- relmatLmer(as.formula(paste0("gene ~ ", abo_group, " + ", covars)), data = dat, relmat = list(ID = kinMat))

  fixedEffectsPv <- car::Anova(lmm)
  fixedEffT2 <- summary(lmm)[[10]]
  rownames(fixedEffT2)[2] <- abo_group
  fixedEffTable <- merge(fixedEffT2,fixedEffectsPv,by='row.names')
  row.names(fixedEffTable) <- fixedEffTable$Row.names
  fixedEffTable$Row.names <- NULL
  coeff <- fixedEffTable[,c(1,2,6)]
  coeff$N <- nobs(lmm)
  return (coeff[abo_group,])
}

make_plots_non_corrected <- function(m, gene, p_sec, p_non){
  #colors = c("#1D91C0", "#EC7014")
  colors = c("#B3B4B6","#225EA8")
  
  pvals <- data.frame(matrix(nrow = 2, ncol = 3))
  colnames(pvals) <- c("group1", "group2", "p")
  pvals[1,] <- c("secretor.A/AB", "secretor.B/O", changeSciNot(formatC(p_sec, digits = 3)))
  pvals[2,] <- c("non-secretor.A/AB", "non-secretor.B/O", changeSciNot(formatC(p_non, digits = 3)))
  
  p <- ggplot(m, aes(x = ABO_A_FUT2, y = gene, color = FUT2_status)) +
    geom_violin(trim=FALSE) +
    geom_boxplot(width=0.5) +
    labs( x = "", y = bquote("ln ("~italic(.(gene))~" RPKM )")) +
    ggtitle(bquote(italic(.(gene)))) +
    theme_classic() + theme(legend.position = "none") +
    scale_color_manual(values = colors) + scale_x_discrete(labels = c("A/AB\nFUT2 secretor", "B/O", "A/AB\nFUT2 non-secretor", "B/O")) +
    stat_pvalue_manual(pvals, label = "p",  y.position = max(m$gene)-1) 

 return(p)
}

make_plots_corrected <- function(m, gene, p_sec, p_non){
  #colors = c("#1D91C0", "#EC7014")
  colors = c("#B3B4B6","#225EA8")
  
  pvals <- data.frame(matrix(nrow = 2, ncol = 3))
  colnames(pvals) <- c("group1", "group2", "p")
  pvals[1,] <- c("secretor.A/AB", "secretor.B/O", changeSciNot(formatC(p_sec, digits = 3)))
  pvals[2,] <- c("non-secretor.A/AB", "non-secretor.B/O", changeSciNot(formatC(p_non, digits = 3)))
  
  p <- ggplot(m, aes(x = ABO_A_FUT2, y = gene_resid, color = FUT2_status)) +
    geom_violin(trim=FALSE) +
    geom_boxplot(width=0.5) +
    labs( x = "", y = bquote("ln ("~italic(.(gene))~" RPKM)")) +
    theme_classic() + theme(legend.position = "none") +
    ggtitle(bquote(italic(.(gene)))) +
    scale_color_manual(values = colors) + scale_x_discrete(labels = c("A/AB\nFUT2 secretor", "B/O", "A/AB\nFUT2 non-secretor", "B/O")) +
    stat_pvalue_manual(pvals, label = "p",  y.position = max(m$gene)-1) 
  
  return(p)
}


make_plots_galnac <- function(m, gene, pval){
  #colors = c("#1D91C0", "#EC7014")
  colors = c("#225EA8", "#B3B4B6")
  
  pvals <- data.frame(matrix(nrow = 1, ncol = 3))
  colnames(pvals) <- c("group1", "group2", "p")
  pvals[1,] <- c("GalNAc-depleted", "GalNAc-enriched", formatC(pval, digits = 3))

  
  p <- ggplot(m, aes(x = ABO_A_vs_rest, y = gene)) +
    geom_violin(trim=FALSE, color = "#225EA8") +
    geom_boxplot(width=0.5, color = "#225EA8") +
    labs( x = "", y = bquote("ln ("~italic(.(gene))~" RPKM)")) +
    theme_classic() + theme(legend.position = "none") +
    scale_x_discrete(labels = c("GalNAc-low", "GalNAc-high")) +
    stat_pvalue_manual(pvals, label = "p",  y.position = max(m$gene)-1) 
  
  return(p)
}
make_plots_galnac_corrected <- function(m, gene, pval){
  #colors = c("#1D91C0", "#EC7014")
  colors = c("#225EA8", "#B3B4B6")
  
  pvals <- data.frame(matrix(nrow = 1, ncol = 3))
  colnames(pvals) <- c("group1", "group2", "p")
  pvals[1,] <- c("GalNAc-depleted", "GalNAc-enriched", changeSciNot(formatC(pval, digits = 3)))
  
  
  p <- ggplot(m, aes(x = ABO_A_vs_rest, y = gene_resid)) +
    geom_violin(trim=FALSE, color = "#225EA8") +
    geom_boxplot(width=0.5, color = "#225EA8") +
    labs( x = "", y = bquote("ln ("~italic(.(gene))~" RPKM)")) +
    theme_classic() + theme(legend.position = "none") +
    scale_x_discrete(labels = c("GalNAc-low", "GalNAc-high")) +
    stat_pvalue_manual(pvals, label = "p",  y.position = max(m$gene) -1) 
  
  return(p)
}

changeSciNot <- function(n) {
  output <- format(n, scientific = TRUE) #Transforms the number into scientific notation even if small
  output <- sub("e", "x10", output) #Replace e with 10^
  output <- sub("\\+0?", "", output) #Remove + symbol and leading zeros on expoent, if > 1
  output <- sub("-0?", "-", output) #Leaves - symbol but removes leading zeros on expoent, if < 1
  output
}

d <- "/groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/umcg-dzhernakova/SV_GWAS/v2/"
c <- "DAG3"

pheno <- as.data.frame(t(read.delim(paste0(d, "/data_fastGWA/shortbred_final.res"), header = T, sep = "\t", as.is = T, check.names = F, row.names = 1)))
pheno <- log(pheno + 1)
pheno <- rename_genes(pheno)
covar <- read.delim(paste0(d, "/data_fastGWA/", c, ".covariates.with_svs.txt"), header = T, sep = "\t", as.is = T, check.names = F, row.names = 1)

fut2_geno <- as.data.frame(t(read.delim(paste0(d, "/genotypes/", c, "/with_relatives/text_genotypes/", c, ".", fut2_snp, ".genotypes.txt"), header = T, sep = "\t", as.is = T, check.names = F, row.names = 1)))
fut2_geno <- get_fut2_secretor_status(fut2_geno)

sex <- read.table(paste0(d, "/genotypes/", c, "/with_relatives/", c , "_filtered_withrel.fam"), header = F,  as.is = T, check.names = F)
row.names(sex) <- sex[,2]

abo <- read.delim(paste0(d, "/data_fastGWA/pheno/", c, ".abo_blood_group.txt"), header = T, sep = "\t", as.is = T, check.names = F, row.names = 1)
abo <- format_abo_bloodgroup(abo)

ids <- intersect(row.names(pheno), row.names(covar))

kinMat <- readRDS(paste0(d, "/genotypes/", c, "/with_relatives/GCTA/GRM_", c, ".text.grm.ibd_matrix.RDS"))*2

gene_list <- c("agaF", "agaV", "agaC", "agaD", "nagA", "agaS", "lacC", "gatY-kbaY",  "gatZ-kbaZ")
cnt <- 1
all_plots <- list()
all_plots_adj <- list()
galnac_plots_adj <- list()
galnac_plots <- list()
assoc_res <- data.frame(matrix(nrow = length(gene_list), ncol = 1+4*3))
assoc_res_galnac <- data.frame(matrix(nrow = length(gene_list), ncol = 5))

for (gene in gene_list){
  
  print(gene)    
  
  m <- cbind(pheno[ids,gene], covar[ids, c("F.prausnitzii:101","F.prausnitzii:102","F.prausnitzii:9","F.prausnitzii:69", "F.prausnitzii:33","age", "read_number", "F.prausnitzii")], sex[ids,5], fut2_geno[ids, "secr"], abo[ids, ])
  colnames(m) <- c("gene", "sv1", "sv2", "sv3", "sv4", "sv5","age", "read_number", "abundance", "sex", "FUT2_status", colnames(abo))
  m <- na.omit(m)
  m$ID <- row.names(m)
  #m$read_number <- scale(m$read_number)
  m$age <- scale(m$age)
  #m$abundance <- scale(m$abundance)
  m$sex <- as.factor(m$sex)
  
  m$ABO_A_FUT2 <- interaction(m$FUT2_status, m$ABO_A)
  m$ABO_A_FUT2 <- factor(m$ABO_A_FUT2, levels = c("secretor.A/AB", "secretor.B/O", "non-secretor.A/AB", "non-secretor.B/O"))
  
  m$ABO_A_vs_rest <- "GalNAc-depleted"
  m[m$ABO_A == "A/AB" & m$FUT2_status == "secretor","ABO_A_vs_rest"] <- "GalNAc-enriched"
  m$ABO_A_vs_rest <- factor(m$ABO_A_vs_rest, levels = c("GalNAc-depleted", "GalNAc-enriched"))
  
  
  sec <- m[m$FUT2_status == "secretor",]
  non <- m[m$FUT2_status == "non-secretor",]
  
  assoc_res[cnt,] <- cbind(gene,fit_lmm(m, "ABO_A"),  fit_lmm(sec, "ABO_A"),  fit_lmm(non, "ABO_A"))
  assoc_res_galnac[cnt,] <- cbind(gene, fit_lmm(m, "ABO_A_vs_rest"))
  
  covars = "age + read_number + sex + abundance + sv1 +sv2 + sv3 + sv4 + sv5 + (1 | ID)"
  m$gene_resid <- residuals(relmatLmer(as.formula(paste0("gene ~ ", covars)), data = m, relmat = list(ID = kinMat)))
  
  all_plots[[cnt]] <- make_plots_non_corrected(m, gene, assoc_res[cnt,8], assoc_res[cnt,12])
  all_plots_adj[[cnt]] <- make_plots_corrected(m, gene, assoc_res[cnt,8], assoc_res[cnt,12])
  
  galnac_plots[[cnt]] <- make_plots_galnac(m, gene, assoc_res_galnac[cnt,4])
  galnac_plots_adj[[cnt]] <- make_plots_galnac_corrected(m, gene, assoc_res_galnac[cnt,4])
  
  cnt <- cnt + 1
}
write.table(assoc_res, file = "assoc/shortbred_assoc_ABO_A.txt", sep = "\t", quote = F, col.names = NA)
write.table(assoc_res_galnac, file = "assoc/shortbred_assoc_galnac.txt", sep = "\t", quote = F, col.names = NA)

pdf("shortbred.all_combined_noadj.pdf", width = 8.3, height = 11.7, useDingbats = F)
plot_grid(all_plots[[1]], all_plots[[2]], all_plots[[3]], all_plots[[4]], all_plots[[5]], all_plots[[6]], all_plots[[7]], all_plots[[8]], all_plots[[9]], nrow = 3, ncol = 3)
dev.off()

pdf("shortbred.all_combined_adj.pdf", width = 8.3, height = 11.7, useDingbats = F)
plot_grid(all_plots_adj[[1]], all_plots_adj[[2]], all_plots_adj[[3]], all_plots_adj[[4]], all_plots_adj[[5]], all_plots_adj[[6]], all_plots_adj[[7]], all_plots_adj[[8]], all_plots_adj[[9]], nrow = 3, ncol = 3)
dev.off()

pdf("shortbred.all_combined_galnac_noadj.pdf", width = 8.3, height = 11.7, useDingbats = F)
plot_grid(galnac_plots[[1]], galnac_plots[[2]], galnac_plots[[3]], galnac_plots[[4]], galnac_plots[[5]], galnac_plots[[6]], galnac_plots[[7]], galnac_plots[[8]], galnac_plots[[9]], nrow = 3, ncol = 3)
dev.off()

pdf("shortbred.all_combined_galnac_adj.pdf", width = 8.3, height = 11.7, useDingbats = F)
plot_grid(galnac_plots_adj[[1]], galnac_plots_adj[[2]], galnac_plots_adj[[3]], galnac_plots_adj[[4]], galnac_plots_adj[[5]], galnac_plots_adj[[6]], galnac_plots_adj[[7]], galnac_plots_adj[[8]], galnac_plots_adj[[9]], nrow = 3, ncol = 3)
dev.off()







