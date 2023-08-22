# ml RPlus/4.0.3-foss-2018c-v21.12.10
library(kinship2, lib.loc="/groups/umcg-lld/tmp01/other-users/umcg-dzhernakova/R/x86_64-pc-linux-gnu-library/4.0.3")
library(lme4qtl, lib.loc="/groups/umcg-lld/tmp01/other-users/umcg-dzhernakova/R/x86_64-pc-linux-gnu-library/4.0.3")
library(lme4)
library(ggplot2)
library(patchwork)
library('ggpubr')

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
  return (fut2_geno)
}

format_abo_bloodgroup <- function(abo){
  abo$ABO_O <- "A/B"
  abo[abo$Bloodtype == "O", "ABO_O"] <- "O"
  abo$ABO_O <- as.factor(abo$ABO_O)
  abo$Bloodtype <- as.factor(abo$Bloodtype)
  abo$ABO_A <- "O/B"
  abo[abo$Bloodtype %in% c("AB", "A"), "ABO_A"] <- "A"
  abo$A_geno <- 0
  abo[abo$Blood_genotype %in% c("A", "AB"), "A_geno"] <- 1
  abo[abo$Blood_genotype == "AA", "A_geno"] <- 2
  abo$ABO_blood_group <- abo$Bloodtype
  return(abo)
}
changeSciNot <- function(n) {
  output <- format(n, scientific = TRUE) #Transforms the number into scientific notation even if small
  output <- sub("e", "x10", output) #Replace e with 10^
  output <- sub("\\+0?", "", output) #Remove + symbol and leading zeros on expoent, if > 1
  output <- sub("-0?", "-", output) #Leaves - symbol but removes leading zeros on expoent, if < 1
  output
}

fit_lmm <- function(dat, abo_group, svtype){
  covars = paste0(abo_group, " + age + read_number + sex + abundance + (1 | ID)")
  dat$age <- scale(dat$age)
  dat$abundance <- scale(dat$abundance)
  if (svtype == "vSV"){
    lmm <- relmatLmer(as.formula(paste0("sv ~ ", covars)), data = dat, relmat = list(ID = kinMat), control = lmerControl(optimizer = "bobyqa"))
    fixedEffectsPv <- car::Anova(lmm)
    fixedEffT2 <- summary(lmm)[[10]]
    rownames(fixedEffT2) <- gsub('sex2','sex',row.names(fixedEffT2))
    rownames(fixedEffT2)[2] <- abo_group
    fixedEffTable <- merge(fixedEffT2,fixedEffectsPv,by='row.names')
    row.names(fixedEffTable) <- fixedEffTable$Row.names
    fixedEffTable$Row.names <- NULL
    coeff <- fixedEffTable[,c(1,2,6)]
    coeff$N <- nobs(lmm)
    
  } else if (svtype == "dSV"){
    lmm <- relmatGlmer(as.formula(paste0("sv ~ ", covars)), data = dat, relmat = list(ID = kinMat), family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
    coeff <- as.data.frame(summary(lmm)$coefficients[,c(1,2,4)])
    #if (coeff[2,3] == 0){
    #  lmm <- relmatGlmer(as.formula(paste0("sv ~ ", covars)), data = dat, relmat = list(ID = kinMat), family = binomial)
    #  coeff <- as.data.frame(summary(lmm)$coefficients[,c(1,2,4)])
    #}
    coeff$N <- nobs(lmm)
  }
  return (coeff[abo_group,])
}

make_plots <- function(m, svtype, sv_name, p_sec, p_non){
  m$ABO_A_FUT2 <- interaction(m$FUT2_status, m$ABO_A)
  m$ABO_A_FUT2 <- factor(m$ABO_A_FUT2, levels = c("secretor.A", "secretor.O/B", "non-secretor.A", "non-secretor.O/B"))
  
  pvals <- data.frame(matrix(nrow = 2, ncol = 3))
  colnames(pvals) <- c("group1", "group2", "p")
  pvals[1,] <- c("secretor.A/AB", "secretor.B/O", changeSciNot(formatC(p_sec, digits = 3)))
  pvals[2,] <- c("non-secretor.A/AB", "non-secretor.B/O", changeSciNot(formatC(p_non, digits = 3)))
  
  
  if (svtype == "vSV"){
    #colors = c("#1D91C0", "#EC7014")
    colors = c("#B3B4B6", "#225EA8")
    
    p <-ggplot(m, aes(x = ABO_A_FUT2, y = sv, color = FUT2_status)) +
      geom_violin(trim=FALSE) +
      geom_boxplot(width=0.5) +
      labs(x = "", y = paste0(svtype, " abundance")) +
      theme_classic() + theme(legend.position = "none") +
      scale_color_manual(values = colors) + 
      scale_x_discrete(labels = c("A/AB\nFUT2 secretors", "B/O", "A/AB\nFUT2 non-secretors", "B/O"))
    
  } else if (svtype == "dSV"){
    colors = c("white", "#B3B4B6", "#225EA8")
    m[m$sv == 0, "sv_text"] <- "deletion"
    m[m$sv == 1, "sv_text"] <- "no deletion"
    m$sv_text <- as.factor(m$sv_text)
    
    p <- ggplot(m, aes(x = ABO_A_FUT2)) +
      geom_bar(aes(fill = interaction(sv_text, FUT2_status)), position = "fill") +
      scale_fill_manual(sv, values = c(colors[1], colors[2], colors[1], colors[3])) +
      theme_classic() + labs(x = "", y = paste0("fraction of samples with ", svtype)) +
      theme(legend.position = "none") +
      scale_x_discrete(labels=c("A/AB\nFUT2 secretors", "B/O", "A/AB\nFUT2 non-secretors", "B/O")) 
    
  }
  return(p)
}

sv <- "F.prausnitzii:102"
svtype <- "dSV"
d <- "/groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/umcg-dzhernakova/SV_GWAS/v3/"
setwd(paste0(d, "/plots/"))

svs <- read.delim("svs_to_plot.txt", sep = "\t", as.is = T)
plots <- list()
cnt <- 1
for (sv_num in 1:nrow(svs)){
  sv <- svs[sv_num, "sv"]
  cohorts <- unlist(strsplit(svs[sv_num, "cohorts"], ",", fixed = T))
  svtype <- svs[sv_num, "svtype"]
  coord <- svs[sv_num, "coord"]
  
  print(sv)
  for (c in cohorts){
    print(c)
  
    pheno <- read.delim(paste0(d, "/data_fastGWA/", c, ".", svtype, ".filtered.txt"), header = T, sep = "\t", as.is = T, check.names = F, row.names = 2)
    covar <- read.delim(paste0(d, "/data_fastGWA/", c, ".covariates.txt"), header = T, sep = "\t", as.is = T, check.names = F, row.names = 2)
    if (svtype == "dSV") pheno[,sv] <- as.factor(pheno[,sv])
    fut2_geno <- as.data.frame(t(read.delim(paste0(d, "/genotypes/", c, "/with_relatives/text_genotypes/", c, ".", fut2_snp, ".genotypes.txt"), header = T, sep = "\t", as.is = T, check.names = F, row.names = 1)))
    fut2_geno <- get_fut2_secretor_status(fut2_geno)
    
    sex <- read.table(paste0(d, "/genotypes/", c, "/with_relatives/", c , "_filtered_withrel.fam"), header = F,  as.is = T, check.names = F)
    row.names(sex) <- sex[,2]
    
    abo <- read.delim(paste0(d, "/data_fastGWA/pheno/", c, ".abo_blood_group.txt"), header = T, sep = "\t", as.is = T, check.names = F, row.names = 1)
    abo <- format_abo_bloodgroup(abo)
    
    ids <- intersect(row.names(pheno), row.names(covar))
    
    m <- cbind(pheno[ids,sv], covar[ids, c(gsub(":[0-9]+","", sv), "age", "read_number")], sex[ids,5], fut2_geno[ids, "secr"], abo[ids, ])
    colnames(m) <- c("sv", "abundance", "age", "read_number", "sex", "FUT2_status", colnames(abo))
    m <- na.omit(m)
    m$sex <- as.factor(m$sex)
    m$ID <- row.names(m)
    #m$age <- scale(m$age)
    #m$abundance <- scale(m$abundance)
    
    kinMat <- readRDS(paste0(d, "/genotypes/", c, "/with_relatives/GCTA/GRM_", c, ".text.grm.ibd_matrix.RDS"))*2
    
    sec <- m[m$FUT2_status == "secretor",]
    non <- m[m$FUT2_status == "non-secretor",]
    
    res <- data.frame(matrix(nrow = 6, ncol = 7))
    colnames(res) <- c("ID", "A1", "A2", "b", "se", "P", "N")
    
    res[1,] <- c("ABO_A", "B/O", "A/AB", fit_lmm(m, "ABO_A", svtype))
    res[2,] <- c("ABO_A_secretors", "B/O", "A/AB", fit_lmm(sec, "ABO_A", svtype))
    res[3,] <- c("ABO_A_nonsecretors", "B/O", "A/AB", fit_lmm(non, "ABO_A", svtype))
    res[4,] <- c("ABO_O", "O", "A/AB/B", fit_lmm(m, "ABO_O", svtype))
    res[5,] <- c("ABO_O_secretors", "O", "A/AB/B", fit_lmm(sec, "ABO_O", svtype))
    res[6,] <- c("ABO_O_nonsecretors", "O", "A/AB/B", fit_lmm(non, "ABO_O", svtype))
    write.table(res, paste0("assoc/", c, ".", sv, "-ABO_association.txt"), sep = "\t", quote = F, row.names = F)
    
    plots[[cnt]] <- make_plots(m, svtype, coord, res[2,6], res[3,6])
    cnt <- cnt + 1
  }
}

pdf("all_combined.pdf", width = 8.3, height = 11.7, useDingbats = F)
(plots[[1]] + plots[[2]] + plots[[3]]) /  (plots[[4]] + plots[[5]] + plot_spacer()) / (plots[[6]] + plots[[7]] + plot_spacer()) / (plots[[8]] + plots[[9]] + plots[[10]]) / (plots[[11]] + plots[[12]] + plots[[13]])
dev.off()



