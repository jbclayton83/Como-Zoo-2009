## R stats analysis of the Como Zoo NHP data, working with Bryan Featherstone from Drexel
## Begun on 2016/12/21

install.packages('biom', repo='http://cran.wustl.edu')
library(biom)
setwd('~/Box Sync/knights_box/pmp_abx/')
cz.rec.otus.biom = read_biom('./ninja97/ninja97_json_sub_otutable.biom')
cz.rec.otus = as.matrix(biom_data(cz.rec.otus.biom))
cz.rec.otus = t(cz.rec.otus)
depths = rowSums(cz.rec.otus)
hist(depths,breaks=30)
otu.counts = colSums(cz.rec.otus > 0)
hist(otu.counts, breaks=30)

# Remove OTUs in less than 10% of samples
cz.rec.otus = cz.rec.otus[,colMeans(cz.rec.otus > 0) >= 0.1]
depths = rowSums(cz.rec.otus)
dim(cz.rec.otus)
otu.counts = colSums(cz.rec.otus > 0)
hist(otu.counts, breaks=30)

# Remove samples with low depth
sort(depths)[1:10]
cz.rec.otus = cz.rec.otus[depths >= 1000,]
dim(cz.rec.otus)

# Load mapping file
cz.rec.map = read.table('map_cz_recent.txt', sep = '\t', header = T, row=1, check.names = F, comment.char = '')
sample.ids = intersect(rownames(cz.rec.otus), rownames(cz.rec.map))
sample.ids = sort(sample.ids)
cz.rec.otus = cz.rec.otus[sample.ids,]
cz.rec.map = cz.rec.map[sample.ids,]
dim(cz.rec.map)
dim(cz.rec.otus)

## STATISTICAL ANALYSIS
library(vegan)
library(biom)
setwd('~/Box Sync/knights_box/pmp_abx/')
genus.biom = read_biom('ninja97/ninja97_json_sub_otutable_L6.biom')
genus = as.matrix(biom_data(genus.biom))
genus = t(genus)
map = read.table('map_cz_recent.txt', sep = '\t', comment.char = '', header = T, row.names = 1)
common.ids = intersect(rownames(map), rownames(genus))
genus = genus[common.ids,]
map = map[common.ids,]
dim(map)
# Drop genera present in <10% of samples
genus = genus[,colMeans(genus > 0) >= 0.1]
dim(genus)
colnames(genus)[1:5]
colnames(map)
table(map$SpeciesCommonNameGeneral)
# What column is prevotella in?
grep('Prevotella$',colnames(genus))
prevotella = genus[,grep('Prevotella$',colnames(genus))]
hist(prevotella)
# Pearson correlation of Prevotella and a continuous variable
cor.test(prevotella, map$Total_number_of_courses_of_antibiotics)
fit = lm(prevotella ~ map$Total_number_of_courses_of_antibiotics)
summary(fit)
pval = anova(fit)['map$Total_number_of_courses_of_antibiotics','Pr(>F)']
pval
qqnorm(rstudent(fit)); abline(0,1)
ks.test(rstudent(fit), pnorm, mean=mean(rstudent(fit)), sd=sd(rstudent(fit)))
# Run linear model controlling for confounders
fit = lm(prevotella ~ map$Total_number_of_courses_of_antibiotics + map$SpeciesCommonName_Original)
summary(fit)
fit = lm(prevotella ~ map$Total_number_of_courses_of_antibiotics + map$DigestivePhysiology)
summary(fit)

# TEST all genera for correlation with abx exposure
# Make a vector of zeros for each genus in the table
pvals = numeric(ncol(genus))
names(pvals) = colnames(genus)
# Loop through genera and test each one
for(i in 1:ncol(genus)) {
  fit = lm(genus[,i] ~ map$Total_number_of_courses_of_antibiotics + map$Sex)
  pvals[i] = anova(fit)['map$Total_number_of_courses_of_antibiotics','Pr(>F)']
}
sort(pvals)[1:10]
# Multiple hypothesis testing with False Discovery Rate
qvals = p.adjust(pvals, 'fdr')
sort(qvals)[1:10]

# Look at the normality of all genera with KS
ks.pvals = numeric(ncol(genus))
names(ks.pvals) = colnames(genus)
options(warn = -1)
for(i in 1:ncol(genus)) {
  fit = lm(genus[,i] ~ map$Total_number_of_courses_of_antibiotics + map$SpeciesCommonName_Original)
  ks.pvals[i] = ks.test(rstudent(fit), pnorm, mean=mean(rstudent(fit)), sd=sd(rstudent(fit)), exact=F)$p.value
}
options(warn = 0)
ks.qvals = p.adjust(ks.pvals, 'fdr')
sort(ks.qvals)[1:10]
plot(-log10(seq(0,1,(dim(genus)[2])+1))[-1]), -log10(sort(ks.pvals))); abline(0,1)
# Not a good set of residuals!

# Install edgeR
source("https://bioconductor.org/biocLite.R")
biocLite('edgeR')
source("~/bin/wrap.edge.R")
genus.biom = read_biom('ninja_genus_absolute/ninja97_otutable_L6_json.biom')
genus.a = as.matrix(biom_data(genus.biom))
genus.a = t(genus.a)
genus.a = genus.a[common.ids,]
genus.a = genus.a[,colMeans(genus.a > 0) >= 0.1]
# Run the Generalized Linear Model on the matrix
result = glm.edgeR(x=map$Total_number_of_courses_of_antibiotics, Y=genus.a)
topTags(result)
pvals = topTags(result,n=Inf)$table[,'PValue']
plot(-log10(seq(0,1,length=((dim(genus.a)[2])+1))[-1]), -log10(sort(pvals))); abline(0,1)

# Now include covariates
result = glm.edgeR(x=map$Total_number_of_courses_of_antibiotics, Y=genus.a, covariates = map$SpeciesCommonName_Original + map$Age_2009)
topTags(result)

result = glm.edgeR(x=map$Total_number_of_courses_of_antibiotics, Y=genus.a, covariates = map$DigestivePhysiology)
topTags(result)





