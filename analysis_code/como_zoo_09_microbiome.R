#Como zoo analysis - unstitched
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

#Como zoo analysis - stitched
require(tibble)
otus <- read.delim("data/stitched/como_stitched_otu.txt",header = 1,row.names = 1,check.names = FALSE)
szscids= c("SZ.SC1.S137.L001.",'SZ.SC2.S138.L001.',"SZ.SC3.S139.L001.","SZ.SC7.S143.L001.","SZ.SC4.S140.L001.","SZ.SC6.S142.L001.","SZ.SC5.S141.L001.")
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
otus <- otus[depths >= 9500,]
dim(otus)
otus <- t(otus)
otus <- data.frame(otus)
otus <- tibble::rownames_to_column(otus, var = "#OTU ID")
write.table(otus, file = "data/stitched/como_stitched_otu_filtered.txt",sep = "\t", quote = FALSE, row.names = FALSE)

# Rename the columns to match the mapping file
otus <- read.delim("data/stitched/como_stitched_otu_filtered.txt",header = 1,row.names = 1,check.names = FALSE)
samples <- colnames(otus)
map_friendly <- lapply(samples, function(x) paste(strsplit(x, split = '\\.')[[1]][1:2], collapse = '.'))
colnames(otus) <- map_friendly
colnames(otus)
otus <- tibble::rownames_to_column(otus, var = "#OTU ID")
write.table(otus, file = "data/stitched/como_stitched_otu_filtered_final.txt",sep = "\t", quote = FALSE, row.names = FALSE)


# Do the same as above for the stitched taxa table
taxa <- read.delim("**taxatablefile**",header = 1,row.names = 1,check.names = FALSE)
szscids= c("SZ.SC1.S137.L001.",'SZ.SC2.S138.L001.',"SZ.SC3.S139.L001.","SZ.SC7.S143.L001.","SZ.SC4.S140.L001.","SZ.SC6.S142.L001.","SZ.SC5.S141.L001.")
taxa <- taxa[, !(names(taxa) %in% szscids)]
samples <- colnames(taxa)
map_friendly <- lapply(samples, function(x) paste(strsplit(x, split = '\\.')[[1]][1:2], collapse = '.'))
colnames(taxa) <- map_friendly
colnames(taxa)
# transpose so that rows are samples and columns are OTUs
taxa <- t(taxa)
depths <- rowSums(taxa)
hist(depths,breaks=30)
taxa.counts <- colSums(taxa > 0)
hist(taxa.counts,breaks=30)
taxa <- taxa[,colMeans(taxa > 0) >= .02]
depths <- rowSums(taxa)
dim(taxa)
taxa.counts <- colSums(taxa > 0)
hist(taxa.counts,breaks=30)
sort(depths)[1:10]
taxa <- taxa[depths >= 9500,]
dim(taxa)
taxa <- t(taxa)
taxa <- data.frame(taxa)
taxa <- tibble::rownames_to_column(taxa, var = "#OTU ID")
write.table(taxa, file = "data/stitched/como_stitched_taxa_filtered_final.txt",sep = "\t", quote = FALSE, row.names = FALSE)

# Format the references to match the tree for Refseq from markerGMG.
otus <- read.delim('data/stitched/como_stitched_otu_filtered_final.txt', header=1, row.names = 1, check.names = F, sep = '\t')
refs <- rownames(otus)
refs <- lapply(refs, function(x) paste(strsplit(x, split = ' ')[[1]][1]))
refs <- lapply(refs, function(x) paste(strsplit(x, split = '_')[[1]][1:2], collapse = ' '))
refs[1]
rownames(otus) <- refs
otus <- tibble::rownames_to_column(otus, var = '#OTU ID')
write.table(otus, file = 'data/stitched/como_stitched_otu_filtered_finalref.txt', sep = '\t', quote = F, row.names = F)











