args <- commandArgs(trailingOnly = TRUE)
library(qqman)

fname = args[1]
d <- read.delim(fname, header = F, as.is = T, check.names = F, sep = "\t")

chisq <- qchisq(d[,7], 1, lower.tail = F)
lambda <- median(chisq) / qchisq(0.5,1)
cat("lambda =", lambda, "\n")
png(paste0(fname,".qqplot.png"), type = "cairo", width = 7, height = 7, units = 'in', res = 400)
qq(d[,7])
dev.off()