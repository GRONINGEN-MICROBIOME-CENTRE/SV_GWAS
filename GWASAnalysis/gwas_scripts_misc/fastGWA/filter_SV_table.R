#
# Filter the SV table to remove SVs with low call rate and low MAF and remove samples with low call rate.
# Write filtered SV tables per cohort and a file showing for each SV the cohorts that has it
#

args <- commandArgs(trailingOnly = TRUE)
library(UpSetR)
library(tibble)
library(dplyr)

infile <- args[1] # Original raw SV table
cohort_table <- args[2] # The cohort information for each sample
sv_type <- args[3]

# Filter the per-cohort SV table by call rate (@param cr), MAF (5%) and number of cases (>80) and controls (>80);
# Draw diagnostic plots, save to @param outpath
run_qc_per_dsv <- function(d, cohort_name, outpath, cr = 0.1){
  qc <- as.data.frame(t(apply(d, 2, function(x) table(factor(x, levels = c(0,1,NA)), exclude = NULL))))
  
  qc$num_called <- qc[,1] + qc[,2]
  qc$presence_rate <- qc[,2]/qc$num_called
  qc$call_rate <- qc$num_called / (qc[,1] + qc[,2] + qc[,3])
  
  write.table(qc, file = paste0(outpath, ".qc_per_dSV.txt"),sep = "\t", quote = F, col.names = NA)
  
  # Make plots:
  pdf(paste0(outpath, ".qc_per_dSV.CR",cr, ".pdf"), useDingbats = F)
  par(mfrow=c(2,2))
  hist(qc$num_called, breaks = 100, main = "Number of called dSV", col = "black")
  hist(qc$call_rate, breaks = 100, xlim = c(0,1),  main = "dSV call rate", col = "black")
  abline(v=cr, col = "red")
  plot(qc$presence_rate, qc$call_rate, pch = 16, cex = 0.5, col = "black", main = "Fraction of 1s vs call rate")
  abline(v=0.05, col = "red")
  abline(v=0.95, col = "red")
  abline(h=cr, col = "red")
  hist(qc$presence_rate, breaks = 100, xlim = c(0,1), col = "black", main = "dSV fraction of 1")
  abline(v=0.05, col = "red")
  abline(v=0.95, col = "red")
  dev.off()
  
  # Get the filtered dSVs
  qc_05 <- qc[qc[,1] > 80 & qc[,2] > 80 & qc$presence_rate > 0.05 & qc$presence_rate < 0.95 & qc$call_rate > cr,]
  
  cat(cohort_name, ": number dSV before QC:", nrow(qc), "number dSV after filtering:", nrow(qc_05), "\n", sep = " ")
  
  d <- d[colnames(d) %in% row.names(qc_05)]
  return(d)
}

# Filter the per-cohort SV table by call rate (@param cr), MAF (5%) and number of cases (>80) and controls (>80);
# Draw diagnostic plots, save to @param outpath
run_qc_per_vsv <- function(d, cohort_name, outpath, cr = 0.1){
  qc <- as.data.frame(sapply(d, function(y) length(which(! is.na(y)))))
  
  n=nrow(d)
  nrow(qc)
  colnames(qc) <- "num_called"
  qc$call_rate <- qc[,1]/n
  
  # Make plots:
  pdf(paste0(outpath, ".qc_per_vSV.CR",cr, ".pdf"), useDingbats = F)
  par(mfrow=c(1,2))
  hist(qc$num_called, breaks = 100, main = "Number of called vSV", col = "black")
  hist(qc$call_rate, breaks = 100, xlim = c(0,1),  main = "vSV call rate", col = "black")
  abline(v=cr, col = "red")
  
  dev.off()
  
  # Get the filtered vSVs
  qc_05 <- qc[qc$call_rate > cr,]
  
  cat(cohort_name, ": number vSV before QC:", nrow(qc), "number vSV after filtering:", nrow(qc_05), "\n", sep = " ")
  d <- d[colnames(d) %in% row.names(qc_05)]
  
  return(d)
}


# Remove the samples with low call rate (@param cr)
run_qc_per_sample <- function(d, cohort_name, outpath, cr = 0.05){
  qc2 <-as.data.frame(apply(d, 1, function(y) length(which(! is.na(y)))))
  n2 <- ncol(d)
  colnames(qc2) <- "num_called"
  
  qc2$call_rate <- qc2$num_called/n2
  
  # Make plots
  pdf(paste0(outpath, ".qc_per_sample.pdf"), useDingbats = F)
  hist(qc2$call_rate, breaks = 100, main = "Call rate per sample", col = "black")
  abline(v=cr, col = "red")
  dev.off()
  
  # Get the filtered samples
  qc2_05 <- qc2[qc2$call_rate > cr,]
  cat(cohort_name, ": number samples before QC:", nrow(qc2), "number samples after filtering:", nrow(qc2_05), "\n", sep = " ")
  
  d <- d[row.names(d) %in% row.names(qc2_05),]
  return(d)
}

#
# Main
#


# Raw SV table:
d <- read.delim(infile, header = T, sep = "\t", as.is = T, check.names = F, row.names = NULL)

# SV name conversion table
conv <- read.delim(paste0("data_fastGWA/", sv_type, "_name_conversion_table.txt"), header = T,  sep = "\t", as.is = T, check.names = F)

# Table with cohort info for each sample id
cohorts <- read.delim(cohort_table, header = T, sep = "\t", as.is = T, check.names = F)

cs <- c("LLD1", "300-OB", "500FG", "DAG3")
cohorts <- cohorts[cohorts$Cohort %in%  cs,]

sv_per_cohort <- matrix(nrow = ncol(d), ncol = length(cs))
row.names(sv_per_cohort) <- colnames(d)
colnames(sv_per_cohort) <- c("LLD", "300OB", "500FG", "DAG3")



for (c in cs){
  #per cohort SV table
  d_cohort <- d[d[,1] %in% cohorts[cohorts$Cohort == c, 1],]
  row.names(d_cohort) <- d_cohort[,1]
  d_cohort <- d_cohort[,2:ncol(d_cohort)]
  
  if (c == "LLD1") c = "LLD"
  if (c == "300-OB") c = "300OB"
  
  # get samples with genotype and covariate data
  geno <- read.delim(paste0("data_fastGWA/",c, ".covariates.txt") , header = T, sep = "\t", as.is = T, check.names = F, row.names = 1)
  cat("Number of samples with dSVs: ", nrow(d_cohort), "\n")
  cat("Number of samples with genotypes and covariates: ", nrow(geno), "\n")
  sample_overlap <- row.names(d_cohort)[row.names(d_cohort) %in% row.names(geno)]
  cat ("Number overlapping samples: ", length(sample_overlap), "\n")
  d_cohort <- d_cohort[row.names(d_cohort) %in% sample_overlap,]
  
  #filter SVs and samples
  if (sv_type == "dSV"){
    d_flt <- run_qc_per_dsv(d_cohort, c, paste0("data_fastGWA/QC/", c, ".dSV_filtering"))
    d_flt <- run_qc_per_sample(d_flt, c, paste0("data_fastGWA/QC/", c, ".dSV_filtering"))
  
    sv_per_cohort[row.names(sv_per_cohort) %in% colnames(d_flt), c] <- 1
  
    #change 0/1 into 1/2
    #d_flt[] <- lapply(as.data.frame(d_flt), function(x) sub(1,2,x))
    #d_flt[] <- lapply(as.data.frame(d_flt), function(x) sub(0,1,x))
  } else if (sv_type == "vSV"){
    d_flt <- run_qc_per_vsv(d_cohort, c, paste0("data_fastGWA/QC/", c, ".vSV_filtering"))
    d_flt <- run_qc_per_sample(d_flt, c, paste0("data_fastGWA/QC/", c, ".vSV_filtering"))
    
    sv_per_cohort[row.names(sv_per_cohort) %in% colnames(d_flt), c] <- 1
    
    # normalize: apply INT
    d_norm <- apply(d_flt, 2, function(x) qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x))))
    
  }
  # rename SVs
  new_sv_ids <-(conv[match(colnames(d_flt), conv$sv_id, nomatch = 0),"new_sv_id"])
  length(new_sv_ids) == ncol(d_flt)
  colnames(d_flt) <- new_sv_ids
  
  d_flt <- d_flt %>% 
    rownames_to_column(var = "#IID")
  #write filtered table
  write.table(d_flt, file = paste0("data_fastGWA/", c, ".", sv_type,".filtered.txt"), sep = "\t", quote = F, row.names = F) 
}


# plot number of SVs and their overlap between cohorts
sv_per_cohort <- as.data.frame(sv_per_cohort)
sv_per_cohort[is.na(sv_per_cohort)] <- 0
pdf(paste0("data_fastGWA/QC/", sv_type, "_overlap_v3.pdf"),  useDingbats = F)
upset(sv_per_cohort, order.by = "freq")
dev.off()


# Write SV in any number of cohorts
sv_per_cohort <- sv_per_cohort[2:nrow(sv_per_cohort),]
new_sv_ids3 <-(conv[match(row.names(sv_per_cohort), conv$sv_id, nomatch = 0),"new_sv_id"])
cat("name conversion check:", length(new_sv_ids3) == nrow(sv_per_cohort), "\n")
row.names(sv_per_cohort) <- new_sv_ids3
sv_per_cohort$cohorts <- NA
sv_per_cohort$cohorts <- apply(sv_per_cohort, 1, function(x) {paste(colnames(sv_per_cohort)[which(x==1)], collapse = ",")})
sv_per_cohort <- sv_per_cohort[sv_per_cohort$cohorts != "",]

write.table(sv_per_cohort, file = paste0("data_fastGWA/", sv_type, "_per_cohort_all_cohorts.txt"), sep = "\t", quote = F, col.names = NA) 


# Write SVs present in > 1 cohort
sv_per_cohort2 <- sv_per_cohort[rowSums(sv_per_cohort[,1:4]) > 1,]
write.table(sv_per_cohort2, file = paste0("data_fastGWA/", sv_type, "_per_cohort.txt"), sep = "\t", quote = F, col.names = NA)



