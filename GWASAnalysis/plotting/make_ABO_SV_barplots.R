library(ggplot2)
library(patchwork)
library(plyr)


cohorts <- c("LLD", "500FG", "DAG3")
snp <- "9:136141870"
rsid <- "rs2519093"
a1 <- "T"
a2 <- "C"

sv <- "F.prausnitzii:102"
sv_name <- "Faecalibacterium cf. prausnitzii KLE1255:577_579"
svtype <- "dSV"

cohorts <- c("LLD", "DAG3")
snp <- "9:136146597"
rsid <- "rs550057"
a1 <- "T"
a2 <- "C"

sv <- "F.prausnitzii:9"
sv_name <- "Faecalibacterium cf. prausnitzii KLE1255:1154_1155"
svtype <- "dSV"

cohorts <- c("LLD", "DAG3")
snp <- "9:136146597"
rsid <- "rs550057"
a1 <- "T"
a2 <- "C"
sv <- "F.prausnitzii:101"
sv_name <- "Faecalibacterium cf. prausnitzii KLE1255:575_577"
svtype <- "dSV"

 cohorts <- c("LLD", "500FG", "DAG3")
 snp <- "9:136146597"
 rsid <- "rs550057"
 a1 <- "T"
 a2 <- "C"
 sv <- "F.prausnitzii:69"
 sv_name <- "Faecalibacterium cf. prausnitzii KLE1255:2910_2911"
 svtype <- "dSV"



cohorts <- c("LLD", "500FG", "DAG3")
snp <- "9:136149399"
rsid <- "rs507666"
a1 <- "A"
a2 <- "G"
sv <- "F.prausnitzii:33"
sv_name <- "Faecalibacterium cf. prausnitzii KLE1255:885_887"
svtype <- "vSV"


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

print_sv_summary <- function(m, sv, svtype){
    if (svtype == "dSV"){

        #print(table(m[,c("sv", "geno_factor")]))
        #print(table(m[m$FUT2_status == "secretor", c("sv", "geno_factor")]))
        #print(table(m[m$FUT2_status == "non-secretor", c("sv", "geno_factor")]))
        
        print(table(m[,c("sv", "ABO_blood_group")]))
        print(table(m[m$FUT2_status == "secretor", c("sv", "ABO_blood_group")]))
        print(table(m[m$FUT2_status == "non-secretor", c("sv", "ABO_blood_group")]))
        

    } else if (svtype == "vSV"){
        sv_mean <- ddply(m, "SNP_FUT2", summarise, grp.mean=mean(sv))
        sv_median <- ddply(m, "SNP_FUT2", summarise, grp.median=median(sv))
        sv_sd <- ddply(m, "SNP_FUT2", summarise, grp.sd=sd(sv))
        print(merge(merge(sv_mean, sv_median,  by = "SNP_FUT2"), sv_sd, by = "SNP_FUT2"))
    }
}


print_association_summary <- function(m, svtype, c){

    sec <- m[m$FUT2_status == "secretor",]
    non <- m[m$FUT2_status == "non-secretor",]

    res <- data.frame(matrix(nrow = 6, ncol = 7))
    colnames(res) <- c("ID", "A1", "A2", "N", "b", "se", "P")
    covars = "abundance + PC1 + PC2 + age + read_number + sex"
    if (c == "DAG3"){
        covars = "abundance + PC1 + PC2 + PC3 + PC4 + PC5 + age + read_number + sex"
    }

    if (svtype == "vSV"){
        #lm_g <- lm(as.formula(paste0("sv ~ genotype + ", covars)), d = m)
        #lm_g_s <- lm(as.formula(paste0("sv ~ genotype + ", covars)), d = sec)
        #lm_g_n <- lm(as.formula(paste0("sv ~ genotype + ", covars)), d = non)

        lm_a <- lm(as.formula(paste0("sv ~ ABO_A + ", covars)), d = m)
        lm_a_s <- lm(as.formula(paste0("sv ~ ABO_A + ", covars)), d = sec)
        lm_a_n <- lm(as.formula(paste0("sv ~ ABO_A + ", covars)), d = non)

        lm_o <- lm(as.formula(paste0("sv ~ ABO_O + ", covars)), d = m)
        lm_o_s <- lm(as.formula(paste0("sv ~ ABO_O + ", covars)), d = sec)
        lm_o_n <- lm(as.formula(paste0("sv ~ ABO_O + ", covars)), d = non)

        res[1,] <- c("ABO_A", "B/O", "A/AB", nobs(lm_a), summary(lm_a)$coefficients[2,1], summary(lm_a)$coefficients[2,2], summary(lm_a)$coefficients[2,4])
        res[2,] <- c("ABO_A_secretors", "B/O", "A/AB", nobs(lm_a_s), summary(lm_a_s)$coefficients[2,1], summary(lm_a_s)$coefficients[2,2], summary(lm_a_s)$coefficients[2,4])
        res[3,] <- c("ABO_A_nonsecretors", "B/O", "A/AB", nobs(lm_a_n), summary(lm_a_n)$coefficients[2,1], summary(lm_a_n)$coefficients[2,2], summary(lm_a_n)$coefficients[2,4])
        res[4,] <- c("ABO_O", "O", "A/AB/B", nobs(lm_o), summary(lm_o)$coefficients[2,1], summary(lm_o)$coefficients[2,2], summary(lm_o)$coefficients[2,4])
        res[5,] <- c("ABO_O_secretors", "O", "A/AB/B", nobs(lm_o_s), summary(lm_o_s)$coefficients[2,1], summary(lm_o_s)$coefficients[2,2], summary(lm_o_s)$coefficients[2,4])
        res[6,] <- c("ABO_O_nonsecretors", "O", "A/AB/B", nobs(lm_o_n), summary(lm_o_n)$coefficients[2,1], summary(lm_o_n)$coefficients[2,2], summary(lm_o_n)$coefficients[2,4])
        write.table(res, paste0("assoc/", c, ".", sv, "-ABO_association.txt"), sep = "\t", quote = F, row.names = F)

    } else if (svtype == "dSV"){
        #lm_g <- glm(as.formula(paste0("sv ~ genotype + ", covars)), d = m, family = binomial(link = "logit"))
        #lm_g_s <- lm(as.formula(paste0("sv ~ genotype + ", covars)), d = sec, family = binomial(link = "logit"))
        #summary(glm(as.formula(paste0("sv ~ genotype + ", covars)), d = non, family = binomial(link = "logit"))

        lm_a <- glm(as.formula(paste0("sv ~ ABO_A + ", covars)), d = m, family = binomial(link = "logit"))
        lm_a_s <- glm(as.formula(paste0("sv ~ ABO_A + ", covars)), d = sec, family = binomial(link = "logit"))
        lm_a_n <- glm(as.formula(paste0("sv ~ ABO_A + ", covars)), d = non, family = binomial(link = "logit"))

        lm_o <- glm(as.formula(paste0("sv ~ ABO_O + ", covars)), d = m, family = binomial(link = "logit"))
        lm_o_s <- glm(as.formula(paste0("sv ~ ABO_O + ", covars)), d = sec, family = binomial(link = "logit"))
        lm_o_n <- glm(as.formula(paste0("sv ~ ABO_O + ", covars)), d = non, family = binomial(link = "logit"))

        res[1,] <- c("ABO_A", "B/O", "A/AB", nobs(lm_a), summary(lm_a)$coefficients[2,1], summary(lm_a)$coefficients[2,2], summary(lm_a)$coefficients[2,4])
        res[2,] <- c("ABO_A_secretors", "B/O", "A/AB", nobs(lm_a_s), summary(lm_a_s)$coefficients[2,1], summary(lm_a_s)$coefficients[2,2], summary(lm_a_s)$coefficients[2,4])
        res[3,] <- c("ABO_A_nonsecretors", "B/O", "A/AB", nobs(lm_a_n), summary(lm_a_n)$coefficients[2,1], summary(lm_a_n)$coefficients[2,2], summary(lm_a_n)$coefficients[2,4])
       res[4,] <- c("ABO_O", "O", "A/AB/B", nobs(lm_o), summary(lm_o)$coefficients[2,1], summary(lm_o)$coefficients[2,2], summary(lm_o)$coefficients[2,4])
        res[5,] <- c("ABO_O_secretors", "O", "A/AB/B", nobs(lm_o_s), summary(lm_o_s)$coefficients[2,1], summary(lm_o_s)$coefficients[2,2], summary(lm_o_s)$coefficients[2,4])
        res[6,] <- c("ABO_O_nonsecretors", "O", "A/AB/B", nobs(lm_o_n), summary(lm_o_n)$coefficients[2,1], summary(lm_o_n)$coefficients[2,2], summary(lm_o_n)$coefficients[2,4])
        write.table(res, paste0("assoc/", c, ".", sv, "-ABO_association.txt"), sep = "\t", quote = F, row.names = F)
    }
}



make_plots_vSV <- function(m, out_fname, svtype, sv_name){
    #colors = c("#1D91C0", "#EC7014")
    colors = c("#225EA8", "#CC4C02")

    p1 <- ggplot(m, aes(x = geno_factor, y = sv_corrected)) +
        geom_violin(trim=FALSE, color = colors[1]) +
        geom_boxplot(width=0.5, color = colors[1]) +
        labs(title=paste0(sv, " : ", rsid), x = "genotype", y = paste0(svtype, " corrected abundance")) +
        theme_classic() + theme(legend.position = "none") +
        scale_color_manual(values = c(colors[1]))

    p2 <- ggplot(m, aes(x = FUT2_status, y = sv_corrected)) +
        geom_violin(trim=FALSE, color = colors[1]) +
        geom_boxplot(width=0.5, color = colors[1]) +
        labs(title=paste0(sv, " : FUT2 status"), x = "FUT2 status", y = paste0(svtype, " corrected abundance")) +
        theme_classic() + theme(legend.position = "none") +
        scale_color_manual(values = colors)

    p3 <- ggplot(m, aes(x = ABO_blood_group, y = sv_corrected)) +
        geom_violin(trim=FALSE, color = colors[1]) +
        geom_boxplot(width=0.5, color = colors[1]) +
        labs(title=paste0(sv, " : ABO blood group"), x = "ABO blood group", y = paste0(svtype, " corrected abundance")) +
        theme_classic() + theme(legend.position = "none") +
        scale_color_manual(values = colors)

    p4 <- ggplot(m, aes(x = Blood_genotype, y = sv_corrected)) +
        geom_violin(trim=FALSE, color = colors[1]) +
        geom_boxplot(width=0.5, color = colors[1]) +
        labs(title=paste0(sv, " : ABO genotype"), x = "ABO blood group geno", y = paste0(svtype, " corrected abundance")) +
        theme_classic() + theme(legend.position = "none") +
        scale_color_manual(values = colors)

    p5 <-ggplot(m, aes(x = ABO_O_FUT2, y = sv_corrected, color = FUT2_status)) +
        geom_violin(trim=FALSE) +
        geom_boxplot(width=0.5) +
        labs(title=paste0(sv, " : ABO O vs FUT2"), x = "", y = paste0(svtype, " corrected abundance")) +
        theme_classic() + theme(legend.position = "none") +
        scale_color_manual(values = colors) + scale_x_discrete(labels = c("O", "A/AB/B", "O", "A/AB/B"))

    p6 <-ggplot(m, aes(x = ABO_A_FUT2, y = sv_corrected, color = FUT2_status)) +
        geom_violin(trim=FALSE) +
        geom_boxplot(width=0.5) +
        labs(title=paste0(sv, " : ABO A vs FUT2"), x = "", y = paste0(svtype, " corrected abundance")) +
        theme_classic() + theme(legend.position = "none") +
        scale_color_manual(values = colors) + scale_x_discrete(labels = c("A/AB", "B/O", "A/AB", "B/O"))

    p7 <-ggplot(m, aes(x = ABO_FUT2, y = sv_corrected, color = FUT2_status)) +
        geom_violin(trim=FALSE) +
        geom_boxplot(width=0.3) +
        labs(title=paste0(sv, " : ABO-FUT2"), x = "", y = paste0(svtype, " corrected abundance")) +
        theme_classic() + theme(axis.text.x = element_text(angle = 45,hjust=0.95,vjust=0.1), legend.position = "none") +
        scale_color_manual(values = colors)

    p8 <-ggplot(m, aes(x = SNP_FUT2, y = sv_corrected, color = FUT2_status)) +
        geom_violin(trim=FALSE) +
        geom_boxplot(width=0.3) +
        labs(title=paste0(sv, " : genotype vs FUT2"), x = "", y = paste0(svtype, " corrected abundance")) +
        theme_classic() + theme(axis.text.x = element_text(angle = 45,hjust=0.95,vjust=0.1), legend.position = "none") +
        scale_color_manual(values = colors)

    pdf(out_fname, width = 12, height = 30)
    par(mfrow=c(4, 2))

    print ((p1 + p2) / (p3 + p4) / (p5 + p6) / p7 / p8)
    dev.off()
    return(list(p3,p4,p6))
}

make_plots_dSV <- function(m, out_fname, svtype, sv_name){
    colors = c("#d2d6db", "#225EA8", "#CC4C02")
    p1 <- ggplot(m, aes(x = geno_factor)) +
        geom_bar(aes(fill = sv_text), position = "fill") +
        scale_fill_manual(sv, values = alpha( c(colors), 0.5)) +
        theme_classic() + labs( x = rsid, y = paste0("fraction of samples with ", svtype))

    p2 <- ggplot(m, aes(x = ABO_blood_group)) +
        geom_bar(aes(fill = sv_text), position = "fill") +
        scale_fill_manual(sv, values = alpha( c(colors), 0.5)) +
        theme_classic() + labs(title =  "SV vs ABO blood group", y = paste0("fraction of samples with ", svtype)) +
        theme(legend.position = "none")

    p3 <- ggplot(m, aes(x = Blood_genotype)) +
        geom_bar(aes(fill = sv_text), position = "fill") +
        scale_fill_manual(sv, values = alpha( c(colors), 0.5)) +
        theme_classic() + labs(title =  "SV vs ABO blood group genotype", y = paste0("fraction of samples with ", svtype)) +
        theme(legend.position = "none")


    p4 <- ggplot(m, aes(x = SNP_FUT2)) +
        geom_bar(aes(fill = interaction(sv_text, FUT2_status)), position = "fill") +
        scale_fill_manual(sv, values = alpha( c(colors[1], colors[2], colors[1], colors[3]), 0.5)) +
        theme_classic() + labs(title = "SV vs genotype by FUT2 status", x = "", y = paste0("fraction of samples with ", svtype)) +
        theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0))


    p5 <- ggplot(m, aes(x = ABO_FUT2)) +
        geom_bar(aes(fill = interaction(sv_text, FUT2_status)), position = "fill") +
        scale_fill_manual(sv, values = alpha( c(colors[1], colors[2], colors[1], colors[3]), 0.5)) +
        theme_classic() + labs(title = "SV vs ABO blood group by FUT2 status", x = "", y = paste0("fraction of samples with ", svtype)) +
        theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0), legend.position = "none")


    p6 <- ggplot(m, aes(x = ABO_O_FUT2)) +
        geom_bar(aes(fill = interaction(sv_text, FUT2_status)), position = "fill") +
        scale_fill_manual(sv, values = alpha( c(colors[1], colors[2], colors[1], colors[3]), 0.5)) +
        theme_classic() + labs(title =  "SV vs blood group O by FUT2 status", x = "", y = paste0("fraction of samples with ", svtype)) +
        theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0), legend.position = "none")

    p7 <- ggplot(m, aes(x = ABO_A_FUT2)) +
        geom_bar(aes(fill = interaction(sv_text, FUT2_status)), position = "fill") +
        scale_fill_manual(sv, values = alpha( c(colors[1], colors[2], colors[1], colors[3]), 0.5)) +
        theme_classic() + labs(title =  "SV vs blood group A by FUT2 status", x = "", y = paste0("fraction of samples with ", svtype)) +
        theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0), legend.position = "none")

    pdf(out_fname, width = 8, height = 19)

    ps <- (p1 + p2) / p3 / p4 / p5 / (p6 + p7)
    print(ps + plot_annotation(title = paste0(svtype, " ", sv_name)))
    dev.off()
    return(list(p2,p3,p7))

}
abo_geno_plots <- list()
abo_plots <- list()
abo_a_plots <- list()
sv_counts <- data.frame()

d <- "/groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/umcg-dzhernakova/SV_GWAS/"
for (c in cohorts){
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

    m$ABO_FUT2 <- interaction(m$FUT2_status, m$Bloodtype)
    m$ABO_FUT2 <- factor(m$ABO_FUT2, levels = c("secretor.A",  "secretor.AB", "secretor.B", "secretor.O", "non-secretor.A", "non-secretor.AB", "non-secretor.B", "non-secretor.O"))
    m$ABO_O_FUT2 <- interaction(m$FUT2_status, m$ABO_O)
    m$ABO_O_FUT2 <- factor(m$ABO_O_FUT2, levels = c("secretor.O", "secretor.A/B", "non-secretor.O", "non-secretor.A/B"))
    m$ABO_A_FUT2 <- interaction(m$FUT2_status, m$ABO_A)
    m$ABO_A_FUT2 <- factor(m$ABO_A_FUT2, levels = c("secretor.A", "secretor.O/B", "non-secretor.A", "non-secretor.O/B"))

    m$SNP_FUT2 <- interaction(m$FUT2_status, m$geno_factor)
    m$SNP_FUT2 <- factor(m$SNP_FUT2, levels = c(paste0("secretor", ".", levels(m$geno_factor)), paste0("non-secretor", ".", levels(m$geno_factor))))

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

	if (svtype == "vSV"){
		plots <- make_plots_vSV(m, paste0(c, ".ABO-", rsid, ".", sv, ".", svtype, ".pdf"), svtype, sv_name)
		abo_plots[[c]] <- plots[[1]]
		abo_geno_plots[[c]] <- plots[[2]]
		abo_a_plots[[c]] <- plots[[3]]

	} else if (svtype == "dSV"){
		# add sdev
		m[m$sv == 1, "sv_text"] <- "deletion"
		m[m$sv == 2, "sv_text"] <- "no deletion"
		m$sv_text <- as.factor(m$sv_text)
		plots <- make_plots_dSV(m, paste0(c, ".ABO-", rsid, ".", sv, ".", svtype,".pdf"), svtype, sv_name)
		abo_plots[[c]] <- plots[[1]]
		abo_geno_plots[[c]] <- plots[[2]]
		abo_a_plots[[c]] <- plots[[3]]
		as.data.frame(table(m[,c("sv", "geno_factor", "FUT2_status")]))
    }
}

if (svtype == "dSV"){
    print(sv_counts)
}



library(patchwork)
pdf("all_plots_combined.pdf", width = 15, height = 20, useDingbats = F)
(plots_fig[[1]] + plots_fig[[2]] + plots_fig[[3]]) /  (plots_fig[[4]] + plots_fig[[5]] + plot_spacer()) / (plots_fig[[6]] + plots_fig[[7]] + plot_spacer()) / (plots_fig[[8]] + plots_fig[[9]] + plots_fig[[10]]) / (plots_fig[[11]] + plots_fig[[12]] + plots_fig[[13]])
dev.off()