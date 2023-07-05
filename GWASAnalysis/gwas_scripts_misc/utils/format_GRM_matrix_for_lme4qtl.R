# ml RPlus/4.0.3-foss-2018c-v21.12.10
args <- commandArgs(trailingOnly = TRUE)
library(kinship2, lib.loc="/groups/umcg-lld/tmp01/other-users/umcg-dzhernakova/R/x86_64-pc-linux-gnu-library/4.0.3")

file_prefix <- args[1]
grm <- read.delim(gzfile(paste0(file_prefix,".gz")), as.is = T, sep = "\t", check.names = F, header = F)
ids <- read.delim(paste0(file_prefix,".id"), as.is = T, sep = "\t", check.names = F, header = F)
grm$id1 <- ids[grm$V1,"V2"]
grm$id2 <- ids[grm$V2,"V2"]

grm[grm$V4 < 0.05, "V4"] <- 0

ibd_matrix <- ibdMatrix(grm$id1, grm$id2, grm$V4)
saveRDS(ibd_matrix, file = paste0(file_prefix, ".ibd_matrix.RDS"))