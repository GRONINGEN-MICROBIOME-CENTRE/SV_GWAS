args <- commandArgs(trailingOnly = TRUE)

f <- read.delim(args[1], header = T,  sep = "\t", as.is = T, check.names = F, row.names = 1)
f[] <- lapply(as.data.frame(f), function(x) sub(1,0,x))
f[] <- lapply(as.data.frame(f), function(x) sub(2,1,x))
f2 <- cbind(a = 0, rownames(f), data.frame(f, row.names=NULL))
colnames(f2)[1] = "#FID"
colnames(f2)[2] = "IID"

write.table(f2, file = args[2], sep = "\t", quote = F, row.names = F) 