# ml RPlus/4.0.3-foss-2018c-v21.12.10
library(kinship2, lib.loc="/groups/umcg-lld/tmp01/other-users/umcg-dzhernakova/R/x86_64-pc-linux-gnu-library/4.0.3")
library(lme4qtl, lib.loc="/groups/umcg-lld/tmp01/other-users/umcg-dzhernakova/R/x86_64-pc-linux-gnu-library/4.0.3")
library(lme4)
library(ggplot2)
library(patchwork)
library(plyr)

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

sv <- "F.prausnitzii:102"
svtype = "dSV"

d <- "/groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/umcg-dzhernakova/SV_GWAS/v3/"
setwd(paste0(d, "/plots/"))
c="DAG3"
print(c)
pheno <- read.delim(paste0(d, "/data_fastGWA/", c, ".", svtype, ".filtered.txt"), header = T, sep = "\t", as.is = T, check.names = F, row.names = 2)
if (svtype == "dSV") pheno[,sv] <- as.factor(pheno[,sv])

covar <- read.delim(paste0(d, "/data_fastGWA/", c, ".covariates.txt"), header = T, sep = "\t", as.is = T, check.names = F, row.names = 2)
fut2_geno <- as.data.frame(t(read.delim(paste0(d, "/genotypes/", c, "/text_genotypes/", c, ".", fut2_snp, ".genotypes.txt"), header = T, sep = "\t", as.is = T, check.names = F, row.names = 1)))
fut2_geno <- get_fut2_secretor_status(fut2_geno)

sex <- read.table(paste0(d, "/genotypes/", c, "/with_relatives/", c , "_filtered_withrel.fam"), header = F,  as.is = T, check.names = F)
row.names(sex) <- sex[,2]

abo <- read.delim(paste0(d, "/data_fastGWA/pheno/", c, ".abo_blood_group.txt"), header = T, sep = "\t", as.is = T, check.names = F, row.names = 1)
abo <- format_abo_bloodgroup(abo)

ids <- intersect(row.names(sex), row.names(covar))

m <- cbind(pheno[ids ,sv], covar[ids, c(gsub(":[0-9]+","", sv), "age", "read_number")], sex[ids,5], fut2_geno[ids, "secr"], abo[ids, ])
colnames(m) <- c("sv","abundance", "age", "read_number", "sex", "FUT2_status", colnames(abo))
m <- na.omit(m)
m$sex <- as.factor(m$sex)
m$ID <- row.names(m)
m$age <- scale(m$age)

kinMat <- readRDS(paste0(d, "/genotypes/", c, "/with_relatives/GCTA/GRM_", c, ".text.grm.ibd_matrix.RDS"))*2


m$ABO_A_FUT2 <- interaction(m$FUT2_status, m$ABO_A)
m$ABO_A_FUT2 <- factor(m$ABO_A_FUT2, levels = c("secretor.A", "secretor.O/B", "non-secretor.A", "non-secretor.O/B"))


# Fprau abundance vs ABO A
# Associations
sec <- m[m$FUT2_status == "secretor",]
non <- m[m$FUT2_status == "non-secretor",]

covars = "ABO_A + age + read_number + sex + (1 | ID)"
lmm_a <- relmatLmer(as.formula(paste0("abundance ~ ", covars)), data = m, relmat = list(ID = kinMat), control = lmerControl(optimizer = "bobyqa"))
Pv_a <- car::Anova(lmm_a)
lmm_s <- relmatLmer(as.formula(paste0("abundance ~ ", covars)), data = sec, relmat = list(ID = kinMat), control = lmerControl(optimizer = "bobyqa"))
Pv_s <- car::Anova(lmm_s)
lmm_n <- relmatLmer(as.formula(paste0("abundance ~ ", covars)), data = non, relmat = list(ID = kinMat), control = lmerControl(optimizer = "bobyqa"))
Pv_n <- car::Anova(lmm_n)
paste("ABO A/AB vs B/O vs F.prausnitzii abundance: all ", Pv_a['ABO_A', 3], "; FUT2 secretors: ", Pv_s['ABO_A', 3], "; FUT2 non-secretors: ", Pv_n['ABO_A', 3])

# Make violin plot

pdf("Fprau_abundance_vs_ABO_A.pdf", width = 6, height = 4, useDingbats = F)
colors = c("#B3B4B6", "#3077B9")
ggplot(m, aes(x = ABO_A_FUT2, y = abundance, color = FUT2_status)) +
      geom_violin(trim=FALSE) +
      geom_boxplot(width=0.5) +
      labs(x = "", y =  "F. prausnitzii abundance") +
      theme_classic() + theme(legend.position = "none") +
      scale_color_manual(values = colors) + scale_x_discrete(labels = c("A/AB", "B/O", "A/AB", "B/O"))
dev.off()


# Interaction analysis
lmm_a_inter <- relmatLmer(abundance ~ sv + ABO_A + ABO_A*sv + age + read_number + sex + (1 | ID), data = m, relmat = list(ID = kinMat), control = lmerControl(optimizer = "bobyqa"))
pv_a <- car::Anova(lmm_a_inter)['sv:ABO_A',3]
lmm_s_inter <- relmatLmer(abundance ~ sv + ABO_A + ABO_A*sv + age + read_number + sex + (1 | ID), data = sec, relmat = list(ID = kinMat), control = lmerControl(optimizer = "bobyqa"))
pv_s <- car::Anova(lmm_s_inter)['sv:ABO_A',3]
lmm_n_inter <- relmatLmer(abundance ~ sv + ABO_A + ABO_A*sv + age + read_number + sex + (1 | ID), data = non, relmat = list(ID = kinMat), control = lmerControl(optimizer = "bobyqa"))
pv_n <- car::Anova(lmm_n_inter)['sv:ABO_A',3]

paste("ABO_A * SV interaction P-value: all: ", pv_a, "FUT2 secretors: ", pv_s, "FUT2 non-secretors: ", pv_n)


