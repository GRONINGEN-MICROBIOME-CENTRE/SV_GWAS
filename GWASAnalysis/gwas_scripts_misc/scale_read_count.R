args <- commandArgs(trailingOnly = TRUE)

f <- read.delim(args[1], header = T, sep = "\t", as.is = T, check.names = F)
f$read_number <- scale(f$read_number)
write.table(f, file = args[2], sep = "\t", quote = F, col.names = NA)