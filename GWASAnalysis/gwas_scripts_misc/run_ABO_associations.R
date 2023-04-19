# ml RPlus/4.0.3-foss-2018c-v21.12.10
library(kinship2, lib.loc="/groups/umcg-lld/tmp01/other-users/umcg-dzhernakova/R/x86_64-pc-linux-gnu-library/4.0.3")
library(lme4qtl, lib.loc="/groups/umcg-lld/tmp01/other-users/umcg-dzhernakova/R/x86_64-pc-linux-gnu-library/4.0.3")



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
svtype <- "dSV"
d <- "/groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/umcg-dzhernakova/SV_GWAS/v2/"

c = "DAG3"

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
m$ID <- row.names(m)
#m$read_number <- scale(m$read_number)
m$age <- scale(m$age)

kinMat <- readRDS(paste0(d, "/genotypes/", c, "/with_relatives/GCTA/GRM_", c, ".text.grm.ibd_matrix.RDS"))*2


sec <- m[m$FUT2_status == "secretor",]
non <- m[m$FUT2_status == "non-secretor",]

res <- data.frame(matrix(nrow = 6, ncol = 7))
colnames(res) <- c("ID", "A1", "A2", "N", "b", "se", "P")
covars = "age + read_number + sex + abundance + (1 | ID)"

if (svtype == "vSV"){
  lmm_a <- relmatLmer(as.formula(paste0("sv ~ ABO_A + ", covars)), data = m, relmat = list(ID = kinMat))
  lmm_a_s <- relmatLmer(as.formula(paste0("sv ~ ABO_A + ", covars)), data = sec, relmat = list(ID = kinMat))
  lmm_a_n <- relmatLmer(as.formula(paste0("sv ~ ABO_A + ", covars)), data = non, relmat = list(ID = kinMat))
  
  lmm_o <- relmatLmer(as.formula(paste0("sv ~ ABO_O + ", covars)), data = m , relmat = list(ID = kinMat))
  lmm_o_s <- relmatLmer(as.formula(paste0("sv ~ ABO_O + ", covars)), data = sec, relmat = list(ID = kinMat))
  lmm_o_n <- relmatLmer(as.formula(paste0("sv ~ ABO_O + ", covars)), data = non, relmat = list(ID = kinMat))
} else if (svtype == "dSV"){
  lmm_a <- relmatGlmer(as.formula(paste0("sv ~ ABO_A + ", covars)), data = m, relmat = list(ID = kinMat), family = binomial)
  lmm_a_s <- relmatGlmer(as.formula(paste0("sv ~ ABO_A + ", covars)), data = sec, relmat = list(ID = kinMat), family = binomial)
  lmm_a_n <- relmatGlmer(as.formula(paste0("sv ~ ABO_A + ", covars)), data = non, relmat = list(ID = kinMat), family = binomial, calc.derivs = FALSE)
  # lmm_a_n <- relmatGlmer(as.formula(paste0("sv ~ ABO_A + ", covars)), data = non, relmat = list(ID = kinMat), family = binomial, calc.derivs = FALSE)
  
  lmm_o <- relmatGlmer(as.formula(paste0("sv ~ ABO_O + ", covars)), data = m , relmat = list(ID = kinMat), family = binomial, calc.derivs = FALSE)
  lmm_o_s <- relmatGlmer(as.formula(paste0("sv ~ ABO_O + ", covars)), data = sec, relmat = list(ID = kinMat), family = binomial)
  lmm_o_n <- relmatGlmer(as.formula(paste0("sv ~ ABO_O + ", covars)), data = non, relmat = list(ID = kinMat), family = binomial , calc.derivs = FALSE)
}


#lmm_abo_a <- relmatGlmer(sv ~ ABO_A + age + read_number + sex + abundance + (1 | ID), data = m, relmat = list(ID = kinMat), family = binomial)
#p <- summary(mod_abo_a)$coefficients["ABO_AO/B", 4]
# 1.175974e-32

#lm_abo_a <- glm(sv ~ ABO_A + age + read_number + sex + abundance, data = m, family = binomial)
#lm_p <- summary(lm_abo_a)$coefficients["ABO_AO/B", 4]
# 3.844771e-41

res[1,] <- c("ABO_A", "B/O", "A/AB", nobs(lmm_a), summary(lmm_a)$coefficients["ABO_AO/B",1], summary(lmm_a)$coefficients["ABO_AO/B",2], summary(lmm_a)$coefficients["ABO_AO/B",4])
res[2,] <- c("ABO_A_secretors", "B/O", "A/AB", nobs(lmm_a_s), summary(lmm_a_s)$coefficients[2,1], summary(lmm_a_s)$coefficients[2,2], summary(lmm_a_s)$coefficients[2,4])
res[3,] <- c("ABO_A_nonsecretors", "B/O", "A/AB", nobs(lmm_a_n), summary(lmm_a_n)$coefficients[2,1], summary(lmm_a_n)$coefficients[2,2], summary(lmm_a_n)$coefficients[2,4])
res[4,] <- c("ABO_O", "O", "A/AB/B", nobs(lmm_o), summary(lmm_o)$coefficients[2,1], summary(lmm_o)$coefficients[2,2], summary(lmm_o)$coefficients[2,4])
res[5,] <- c("ABO_O_secretors", "O", "A/AB/B", nobs(lmm_o_s), summary(lmm_o_s)$coefficients[2,1], summary(lmm_o_s)$coefficients[2,2], summary(lmm_o_s)$coefficients[2,4])
res[6,] <- c("ABO_O_nonsecretors", "O", "A/AB/B", nobs(lmm_o_n), summary(lmm_o_n)$coefficients[2,1], summary(lmm_o_n)$coefficients[2,2], summary(lmm_o_n)$coefficients[2,4])
write.table(res, paste0("assoc/", c, ".", sv, "-ABO_association.txt"), sep = "\t", quote = F, row.names = F)









