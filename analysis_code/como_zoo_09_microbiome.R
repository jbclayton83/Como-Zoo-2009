#Como zoo analysis
require(tibble)
otus <- read.delim("data/como_aligned_otutable.txt",header = 1,row.names = 1,check.names = FALSE)
szscids= c("SZ.SC1.S137.L001.R1",'SZ.SC2.S138.L001.R1',"SZ.SC3.S139.L001.R1","SZ.SC7.S143.L001.R1","SZ.SC4.S140.L001.R1","SZ.SC6.S142.L001.R1","SZ.SC5.S141.L001.R1")
otus <- otus[, !(names(otus) %in% szscids)]
# transpose so that rows are samples and columns are OTUs
otus <- t(otus)
depths <- rowSums(otus)
hist(depths,breaks=30)
otu.counts <- colSums(otus > 0)
hist(otu.counts,breaks=30)
otus <- otus[,colMeans(otus > 0) >= .02]
depths <- rowSums(otus)
dim(otus)
otu.counts <- colSums(otus > 0)
hist(otu.counts,breaks=30)
sort(depths)[1:10]
otus <- otus[depths >= 8900,]
dim(otus)
otus <- t(otus)
otus <- data.frame(otus)
otus <- tibble::rownames_to_column(otus, var = "#OTU ID")
write.table(otus, file = "como_aligned_otutable_filtered2.txt",sep = "\t", quote = FALSE, row.names = FALSE)
