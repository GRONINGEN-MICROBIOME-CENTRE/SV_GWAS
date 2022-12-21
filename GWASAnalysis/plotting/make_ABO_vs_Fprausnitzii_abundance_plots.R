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


cohorts <- c("LLD", "500FG", "DAG3")
snp <- "9:136141870"
rsid <- "rs2519093"
a1 <- "T"
a2 <- "C"
sv <- "F.prausnitzii:102"
sv_name <- "Faecalibacterium cf. prausnitzii KLE1255:577_579"
svtype <- "dSV"

d <- "/groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/umcg-dzhernakova/SV_GWAS/"
c="DAG3"
#for (c in cohorts){
    print(c)
    pheno <- read.delim(paste0(d, "/data/", c, ".", svtype, ".filtered.txt"), header = T, sep = "\t", as.is = T, check.names = F, row.names = 1)
    covar <- read.delim(paste0(d, "/data/", c, ".covariates.txt"), header = T, sep = "\t", as.is = T, check.names = F, row.names = 1)
    if (svtype == "dSV") pheno[,sv] <- as.factor(pheno[,sv])
    geno <- as.data.frame(t(read.delim(paste0(d, "/genotypes/", c, "/text_genotypes/", c, ".", snp, ".genotypes.txt"), header = T, sep = "\t", as.is = T, check.names = F, row.names = 1)))
    geno <- recode_genotypes(geno, a1, a2)
    fut2_geno <- as.data.frame(t(read.delim(paste0(d, "/genotypes/", c, "/text_genotypes/", c, ".", fut2_snp, ".genotypes.txt"), header = T, sep = "\t", as.is = T, check.names = F, row.names = 1)))
    fut2_geno <- get_fut2_secretor_status(fut2_geno)

    sex <- read.table(paste0(d, "/genotypes/", c, "/", c , "_filtered.fam"), header = F,  as.is = T, check.names = F)
    row.names(sex) <- sex[,2]

    abo <- read.delim(paste0(d, "/data/pheno/", c, ".abo_blood_group.txt"), header = T, sep = "\t", as.is = T, check.names = F, row.names = 1)
    abo <- format_abo_bloodgroup(abo)

    ids <- intersect(row.names(pheno), row.names(covar))

    if (c != "DAG3"){
        m <- cbind(pheno[ids,sv], covar[ids, c(gsub(":[0-9]+","", sv), "PC1", "PC2", "age", "read_number")], sex[ids,5], geno[ids,], fut2_geno[ids, "secr"], abo[ids, ])
        colnames(m) <- c("sv", "abundance", "PC1", "PC2", "age", "read_number", "sex", "geno_factor", "genotype", "FUT2_status", colnames(abo))

        # SV vs genotype
        full_formula <- as.formula("sv ~ genotype + abundance + PC1 + PC2 + age + read_number + sex")
        covar_formula <- as.formula("sv ~ abundance + PC1 + PC2 + age + read_number + sex")
        covar_formula_abundance <- as.formula("abundance ~  PC1 + PC2  + age + read_number + sex")
    } else {
        m <- cbind(pheno[ids,sv], covar[ids, c(gsub(":[0-9]+","", sv), "PC1", "PC2", "PC3", "PC4", "PC5", "age", "read_number")], sex[ids,5], geno[ids,], fut2_geno[ids, "secr"], abo[ids, ])
        colnames(m) <- c("sv", "abundance", "PC1", "PC2",  "PC3", "PC4", "PC5", "age", "read_number", "sex", "geno_factor", "genotype", "FUT2_status", colnames(abo))

        # SV vs genotype
        full_formula <- as.formula("sv ~ genotype + abundance + PC1 + PC2 + PC3 + PC4 + PC5 + age + read_number + sex")
        covar_formula <- as.formula("sv ~ abundance + PC1 + PC2 + PC3 + PC4 + PC5 + age + read_number + sex")
        covar_formula_abundance <- as.formula("abundance ~ PC1 + PC2 + PC3 + PC4 + PC5 + age + read_number + sex")
    }
    m <- na.omit(m)

    m$ABO_A_FUT2 <- interaction(m$FUT2_status, m$ABO_A)
    m$ABO_A_FUT2 <- factor(m$ABO_A_FUT2, levels = c("secretor.A", "secretor.O/B", "non-secretor.A", "non-secretor.O/B"))

    
    if (svtype == "dSV"){
        fit <- glm(full_formula, data = m, family = binomial(link = "logit"))
        fit_covar <- glm(covar_formula, data = m, family = binomial(link = "logit"))
        fit_covar_abund <- lm(covar_formula_abundance, data = m)
        sv_counts <- rbind(sv_counts, cbind(as.data.frame(table(m[,c("sv", "FUT2_status","ABO_A")])), c))

    } else if (svtype == "vSV"){
        fit <- lm(full_formula, data = m, )
        fit_covar <- lm(covar_formula, data = m)
        fit_covar_abund <- lm(covar_formula_abundance, data = m)
    }

    print_sv_summary(m, sv, svtype)
    print_association_summary(m, svtype, c)
    m$sv_corrected <- residuals(fit_covar)
    m$abundance_corrected <- residuals(fit_covar_abund)


	# Fprau abundance vs ABO A
	# Associations
	sec <- m[m$FUT2_status == "secretor",]
    non <- m[m$FUT2_status == "non-secretor",]

	covars = "PC1 + PC2 + age + read_number + sex"
    if (c == "DAG3"){
        covars = "PC1 + PC2 + PC3 + PC4 + PC5 + age + read_number + sex"
    }
	lm_a <- lm(as.formula(paste0("abundance ~ ABO_A + ", covars)), d = m)
    lm_a_s <- lm(as.formula(paste0("abundance ~ ABO_A + ", covars)), d = sec)
    lm_a_n <- lm(as.formula(paste0("abundance ~ ABO_A + ", covars)), d = non)
    
	paste("ABO A/AB vs B/O vs F.prausnitzii abundance: all ", summary(lm_a)$coefficients[2,4], "; FUT2 secretors: ", summary(lm_a_s)$coefficients[2,4], "; FUT2 non-secretors: ", summary(lm_a_n)$coefficients[2,4])

	# Make violin plot
	m$abundance_corrected <- residuals(lm(covar_formula_abundance, data = m))
	pdf("Fprau_abundance_vs_ABO_A.pdf", width = 6, height = 4, useDingbats = F)
	colors = c("#B3B4B6", "#3077B9")
	ggplot(m, aes(x = ABO_A_FUT2, y = abundance_corrected, color = FUT2_status)) +
        geom_violin(trim=FALSE) +
        geom_boxplot(width=0.5) +
        labs(x = "", y =  "F. prausnitzii adjusted abundance") +
        theme_classic() + theme(legend.position = "none") +
        scale_color_manual(values = colors) + scale_x_discrete(labels = c("A/AB", "B/O", "A/AB", "B/O"))
	dev.off()
	
	
	# Interaction analysis
	lm_a_inter <- lm(as.formula(paste0("abundance ~ ABO_A + sv + ABO_A*sv + ", covars)), d = m)
    lm_a_s_inter <- lm(as.formula(paste0("abundance ~ ABO_A + sv + ABO_A*sv + ", covars)), d = sec)
    lm_a_n_inter <- lm(as.formula(paste0("abundance ~ ABO_A + sv + ABO_A*sv + ", covars)), d = non)
    paste("ABO_A * SV interaction P-value: all: ", summary(lm_a_inter)$coefficients["ABO_AO/B:sv2",4], "FUT2 secretors: ", summary(lm_a_s_inter)$coefficients["ABO_AO/B:sv2",4], "FUT2 non-secretors: ", summary(lm_a_n_inter)$coefficients["ABO_AO/B:sv2",4])
#}

