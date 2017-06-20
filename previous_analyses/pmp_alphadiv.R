# Alpha Diversity analysis of PMP-abx data with Bryan Featherstone from Drexel
# Begun on 12/28/2016

library(ggplot2)
library(dplyr)

#
# Chao1 diversity analysis

setwd('~/Box Sync/knights_box/pmp_abx/alphadiv/')
chao1.data = read.delim('chao1_full.txt', header = T, check.names = F)
colnames(chao1.data) = c('sample_id', 'chao1_distance')
map = read.delim('../map_cz_full.txt', header = T, row.names = 1, check.names = F)
sample.ids = intersect(rownames(map), chao1.data$sample_id)
sample.ids = sort(sample.ids)
chao1.data = chao1.data[order(chao1.data$sample_id),]
map = map[sample.ids,]

table(map$SpeciesCommonName_Original, map$DigestivePhysiology)

p = ggplot(chao1.data, aes(x=map$Antibiotics_ever_used, y=chao1.data$chao1)) + geom_boxplot(fill="#009E73", width = 0.6, outlier.size = 1.1)
p = p + labs(x="Antibiotics ever used?", y="chao1")
p = p + theme(panel.grid.major = element_blank()) + theme(panel.grid.minor = element_blank())
p = p + theme(axis.title.x = element_text(size= rel(1.2)), axis.title.y = element_text(size = rel(1.2)))
p = p + theme(axis.text.y = element_text(size = rel(1.3), colour = 'black'), axis.text.x = element_text(size = rel(1.3), color = 'black'))
p = p + theme(panel.background = element_rect(color="grey10"))
p

p = ggplot(chao1.data, aes(x=map$DigestivePhysiology, y=chao1.data$chao1_distance)) + geom_boxplot(fill="light blue", width = 0.6, outlier.size = 1.1)
p = p + labs(x="digestive physiology", y="chao1")
p = p + theme(panel.grid.major = element_blank()) + theme(panel.grid.minor = element_blank())
p = p + theme(axis.title.x = element_text(size= rel(1.2)), axis.title.y = element_text(size = rel(1.2)))
p = p + theme(axis.text.y = element_text(size = rel(1.3), colour = 'black'), axis.text.x = element_text(size = rel(1.3), color = 'black'))
p = p + theme(panel.background = element_rect(color="grey10"))
p

p = ggplot(chao1.data, aes(x=map$SpeciesCommonName_Original, y=chao1.data$chao1, fill=map$DigestivePhysiology)) + geom_boxplot(aes(width = 0.6, outlier.size = 1.1))
p = p + labs(x="primate species", y="chao1") + xlim("Black-handed spider monkey", "Blue-eyed black lemur", "De Brazzas monkey", "Emperor tamarin", "Geoffroys tamarin", "White-faced saki", "Orangutan", "Western lowland gorilla")
p = p + theme(panel.grid.major = element_blank()) + theme(panel.grid.minor = element_blank())
p = p + theme(axis.title.x = element_text(size= rel(1.2)), axis.title.y = element_text(size = rel(1.2)))
p = p + theme(axis.text.y = element_text(size = rel(1.3), colour = 'black'), axis.text.x = element_text(size = rel(1.3), color = 'black', angle = 45, vjust = 1.01, hjust = 1))
p = p + theme(panel.background = element_rect(color="grey10"))
p = p + scale_fill_discrete(name="digestive\nphysiology")+ theme(legend.position="top")
p = p + theme(plot.margin = unit(c(0.2,0.2,0.2,1.3), 'cm'))
p

chao1.hindgut = filter(chao1.data, map$DigestivePhysiology == "Hindgut-fermenter")
hindgut.ids = intersect(rownames(map), chao1.hindgut$sample_id)
hindgut.ids = sort(hindgut.ids)
map.hindgut = map[hindgut.ids,]

chao1.simple = filter(chao1.data, map$DigestivePhysiology == "Simple-gut")
simple.ids = intersect(rownames(map), chao1.simple$sample_id)
simple.ids = sort(simple.ids)
map.simple = map[simple.ids,]

p = ggplot(chao1.data, aes(x=map$DateCollection, y=chao1.data$chao1_distance, group=map$AnimalID, color=map$Antibiotics_ever_used)) + geom_line()
p = p + labs(x="collection date", y="chao1")
p = p + scale_color_discrete(name="antibiotic\nexposure")+ theme(legend.position="right")
p

p = ggplot(chao1.hindgut, aes(x=map.hindgut$DateCollection, y=chao1.hindgut$chao1_distance, group=map.hindgut$AnimalID, color=map.hindgut$Antibiotics_ever_used)) + geom_line()
p = p + labs(x="collection date", y="chao1")
p = p + scale_color_discrete(name="antibiotic\nexposure")+ theme(legend.position="right")
p

p = ggplot(chao1.simple, aes(x=map.simple$DateCollection, y=chao1.simple$chao1_distance, group=map.simple$AnimalID, color=map.simple$Antibiotics_ever_used)) + geom_line()
p = p + labs(x="collection date", y="chao1")
p = p + scale_color_discrete(name="antibiotic\nexposure")+ theme(legend.position="right")
p

p = ggplot(chao1.data, aes(x=map$DateCollection, y=chao1.data$chao1_distance, group=map$AnimalID, color=map$DigestivePhysiology)) + geom_line()
p = p + labs(x="collection date", y="chao1")
p = p + scale_color_discrete(name="digestion")+ theme(legend.position="right")
p



#
#
#
# Shannon diversity
# Note the use of a rarefied dataset here, at 16,086 seqs/sample.

setwd('~/Box Sync/knights_box/pmp_abx/alphadiv/')
shannon.data = read.delim('shannon_rare_16k.txt', header = T, check.names = F)
colnames(shannon.data) = c('sample_id', 'shannon_diversity')
map = read.delim('../map_cz_full.txt', header = T, row.names = 1, check.names = F)
sample.ids = intersect(rownames(map), shannon.data$sample_id)
sample.ids = sort(sample.ids)
shannon.data = shannon.data[order(shannon.data$sample_id),]
map = map[sample.ids,]


p = ggplot(shannon.data, aes(x=map$Antibiotics_ever_used, y=shannon.data$shannon_diversity)) + geom_boxplot(fill="#009E73", width = 0.6, outlier.size = 1.1)
p = p + labs(x="Antibiotics ever used?", y="Shannon diversity")
p = p + theme(panel.grid.major = element_blank()) + theme(panel.grid.minor = element_blank())
p = p + theme(axis.title.x = element_text(size= rel(1.2)), axis.title.y = element_text(size = rel(1.2)))
p = p + theme(axis.text.y = element_text(size = rel(1.3), colour = 'black'), axis.text.x = element_text(size = rel(1.3), color = 'black'))
p = p + theme(panel.background = element_rect(color="grey10"))
p

p = ggplot(shannon.data, aes(x=map$DigestivePhysiology, y=shannon.data$shannon_diversity)) + geom_boxplot(fill="light blue", width = 0.6, outlier.size = 1.1)
p = p + labs(x="digestive physiology", y="Shannon diversity")
p = p + theme(panel.grid.major = element_blank()) + theme(panel.grid.minor = element_blank())
p = p + theme(axis.title.x = element_text(size= rel(1.2)), axis.title.y = element_text(size = rel(1.2)))
p = p + theme(axis.text.y = element_text(size = rel(1.3), colour = 'black'), axis.text.x = element_text(size = rel(1.3), color = 'black'))
p = p + theme(panel.background = element_rect(color="grey10"))
p

p = ggplot(shannon.data, aes(x=map$SpeciesCommonName_Original, y=shannon.data$shannon_diversity, fill=map$DigestivePhysiology)) + geom_boxplot(aes(width = 0.6, outlier.size = 1.1))
p = p + labs(x="primate species", y="Shannon diversity") + xlim("Black-handed spider monkey", "Blue-eyed black lemur", "De Brazzas monkey", "Emperor tamarin", "Geoffroys tamarin", "White-faced saki", "Orangutan", "Western lowland gorilla")
p = p + theme(panel.grid.major = element_blank()) + theme(panel.grid.minor = element_blank())
p = p + theme(axis.title.x = element_text(size= rel(1.2)), axis.title.y = element_text(size = rel(1.2)))
p = p + theme(axis.text.y = element_text(size = rel(1.3), colour = 'black'), axis.text.x = element_text(size = rel(1.3), color = 'black', angle = 45, vjust = 1.01, hjust = 1))
p = p + theme(panel.background = element_rect(color="grey10"))
p = p + scale_fill_discrete(name="digestive\nphysiology")+ theme(legend.position="top")
p = p + theme(plot.margin = unit(c(0.1,0.2,0.2,1.8), 'cm'))
p


#
# PD whole distance
## Note the use of a rarefied dataset here, at 16,086 seqs/sample.

setwd('~/Box Sync/knights_box/pmp_abx/alphadiv/')
PDwhole.data = read.delim('pdwholetree_rare_16k.txt', header = T, check.names = F)
colnames(PDwhole.data) = c('sample_id', 'PD_whole_tree')
map = read.delim('../map_cz_adult.txt', header = T, row.names = 1, check.names = F)
sample.ids = intersect(rownames(map), PDwhole.data$sample_id)
sample.ids = sort(sample.ids)
PDwhole.data = PDwhole.data[order(PDwhole.data$sample_id),]
rownames(PDwhole.data) = PDwhole.data[,1]
#PDwhole.data[,1] = NULL
PDwhole.data = PDwhole.data[sample.ids,]
PDwhole.data[,1] = NULL
map = map[sample.ids,]
dim(PDwhole.data)
dim(map)

p = ggplot(PDwhole.data, aes(x=map$Antibiotics_ever_used, y=PDwhole.data$PD_whole_tree)) + geom_boxplot(fill="#009E73", width = 0.6, outlier.size = 1.1)
p = p + labs(x="Antibiotics ever used?", y="phylogenetic distance")
p = p + theme(panel.grid.major = element_blank()) + theme(panel.grid.minor = element_blank())
p = p + theme(axis.title.x = element_text(size= rel(1.2)), axis.title.y = element_text(size = rel(1.2)))
p = p + theme(axis.text.y = element_text(size = rel(1.3), colour = 'black'), axis.text.x = element_text(size = rel(1.3), color = 'black'))
p = p + theme(panel.background = element_rect(color="grey10"))
p

plot(PDwhole.data$PD_whole_tree)
ks.test(PDwhole.data$PD_whole_tree, pnorm, mean=mean(PDwhole.data$PD_whole_tree, sd=sd(PDwhole.data$PD_whole_tree)))
# Definitely not a normal distribution, using non-parametric test
PD.abx.mannpval = wilcox.test(PD_whole_tree ~ map$Antibiotics_ever_used, data=PDwhole.data)
PD.abx.mannpval

p = ggplot(PDwhole.data, aes(x=map$DateCollection, y=PDwhole.data$PD_whole_tree, group=map$AnimalID, color=map$Antibiotics_ever_used)) + geom_line()
p = p + labs(x="collection date", y="phylogenetic distance")
p = p + scale_color_discrete(name="antibiotic\nexposure")+ theme(legend.position="right")
p

# Sex and antibiotics
p = ggplot(map, aes(x=map$Sex, y=map$Total_number_of_courses_of_antibiotics)) + geom_boxplot(fill="light yellow", width = 0.6, outlier.size = 1.1)
p = p + labs(x="sex", y="total number of antibiotics courses")
p = p + theme(panel.grid.major = element_blank()) + theme(panel.grid.minor = element_blank())
p = p + theme(axis.title.x = element_text(size= rel(1.2)), axis.title.y = element_text(size = rel(1.2)))
p = p + theme(axis.text.y = element_text(size = rel(1.3), colour = 'black'), axis.text.x = element_text(size = rel(1.3), color = 'black'))
p = p + theme(panel.background = element_rect(color="grey10"))
p

sex.abx.mannpval = wilcox.test(map$Total_number_of_courses_of_antibiotics ~ map$Sex, data=map)
sex.abx.mannpval

table(map$FecalConsistency, map$Total_number_of_courses_of_antibiotics)

p = ggplot(PDwhole.data, aes(x=map$DigestivePhysiology, y=PDwhole.data$PD_whole_tree)) + geom_boxplot(fill="light blue", width = 0.6, outlier.size = 1.1)
p = p + labs(x="digestive physiology", y="phylogenetic distance")
p = p + theme(panel.grid.major = element_blank()) + theme(panel.grid.minor = element_blank())
p = p + theme(axis.title.x = element_text(size= rel(1.2)), axis.title.y = element_text(size = rel(1.2)))
p = p + theme(axis.text.y = element_text(size = rel(1.3), colour = 'black'), axis.text.x = element_text(size = rel(1.3), color = 'black'))
p = p + theme(panel.background = element_rect(color="grey10"))
p

p = ggplot(PDwhole.data, aes(x=map$Sex, y=PDwhole.data$PD_whole_tree)) + geom_boxplot(fill="light blue", width = 0.6, outlier.size = 1.1)
p = p + labs(x="sex", y="phylogenetic distance")
p = p + theme(panel.grid.major = element_blank()) + theme(panel.grid.minor = element_blank())
p = p + theme(axis.title.x = element_text(size= rel(1.2)), axis.title.y = element_text(size = rel(1.2)))
p = p + theme(axis.text.y = element_text(size = rel(1.3), colour = 'black'), axis.text.x = element_text(size = rel(1.3), color = 'black'))
p = p + theme(panel.background = element_rect(color="grey10"))
p

p = ggplot(PDwhole.data, aes(x=map$SpeciesCommonName_Original, y=PDwhole.data$PD_whole_tree, fill=map$DigestivePhysiology)) + geom_boxplot(aes(width = 0.6, outlier.size = 1.1))
p = p + labs(x="primate species", y="phylogenetic distance") + xlim("Black-handed spider monkey", "Blue-eyed black lemur", "De Brazzas monkey", "Emperor tamarin", "Geoffroys tamarin", "White-faced saki", "Orangutan", "Western lowland gorilla")
p = p + theme(panel.grid.major = element_blank()) + theme(panel.grid.minor = element_blank())
p = p + theme(axis.title.x = element_text(size= rel(1.2)), axis.title.y = element_text(size = rel(1.2)))
p = p + theme(axis.text.y = element_text(size = rel(1.3), colour = 'black'), axis.text.x = element_text(size = rel(1.3), color = 'black', angle = 45, vjust = 1.01, hjust = 1))
p = p + theme(panel.background = element_rect(color="grey10"))
p = p + scale_fill_discrete(name="digestive\nphysiology")+ theme(legend.position="top")
p = p + theme(plot.margin = unit(c(0.1,0.2,0.2,1.8), 'cm'))
p


#
#
#
## Now same as above using only the most recent data timepoint from each monkey-- only adults again
#

setwd('~/Box Sync/knights_box/pmp_abx/alphadiv/')
PDwhole.data = read.delim('pdwholetree_rare_16k.txt', header = T, check.names = F)
colnames(PDwhole.data) = c('sample_id', 'PD_whole_tree')
map = read.delim('../map_cz_recent_adult.txt', header = T, row.names = 1, check.names = F)
sample.ids = intersect(rownames(map), PDwhole.data$sample_id)
sample.ids = sort(sample.ids)
PDwhole.data = PDwhole.data[order(PDwhole.data$sample_id),]
rownames(PDwhole.data) = PDwhole.data[,1]
#PDwhole.data[,1] = NULL
PDwhole.data = PDwhole.data[sample.ids,]
PDwhole.data[,1] = NULL
map = map[sample.ids,]
dim(PDwhole.data)
dim(map)

p = ggplot(PDwhole.data, aes(x=map$Antibiotics_ever_used, y=PDwhole.data$PD_whole_tree)) + geom_boxplot(fill="#009E73", width = 0.6, outlier.size = 1.1)
p = p + labs(x="Antibiotics ever used?", y="phylogenetic distance")
p = p + theme(panel.grid.major = element_blank()) + theme(panel.grid.minor = element_blank())
p = p + theme(axis.title.x = element_text(size= rel(1.2)), axis.title.y = element_text(size = rel(1.2)))
p = p + theme(axis.text.y = element_text(size = rel(1.3), colour = 'black'), axis.text.x = element_text(size = rel(1.3), color = 'black'))
p = p + theme(panel.background = element_rect(color="grey10"))
p

plot(PDwhole.data$PD_whole_tree)
ks.test(PDwhole.data$PD_whole_tree, pnorm, mean=mean(PDwhole.data$PD_whole_tree, sd=sd(PDwhole.data$PD_whole_tree)))
# Definitely not a normal distribution, using non-parametric test
PD.abx.mannpval = wilcox.test(PDwhole.data$PD_whole_tree ~ map$Antibiotics_ever_used, data=PDwhole.data)
PD.abx.mannpval


# Sex and antibiotics
p = ggplot(map, aes(x=map$Sex, y=map$Total_number_of_courses_of_antibiotics)) + geom_boxplot(fill="light yellow", width = 0.6, outlier.size = 1.1)
p = p + labs(x="sex", y="total number of antibiotics courses")
p = p + theme(panel.grid.major = element_blank()) + theme(panel.grid.minor = element_blank())
p = p + theme(axis.title.x = element_text(size= rel(1.2)), axis.title.y = element_text(size = rel(1.2)))
p = p + theme(axis.text.y = element_text(size = rel(1.3), colour = 'black'), axis.text.x = element_text(size = rel(1.3), color = 'black'))
p = p + theme(panel.background = element_rect(color="grey10"))
p

sex.abx.mannpval = wilcox.test(map$Total_number_of_courses_of_antibiotics ~ map$Sex, data=map)
sex.abx.mannpval

table(map$FecalConsistency, map$Total_number_of_courses_of_antibiotics)

p = ggplot(PDwhole.data, aes(x=map$DigestivePhysiology, y=PDwhole.data$PD_whole_tree)) + geom_boxplot(fill="light blue", width = 0.6, outlier.size = 1.1)
p = p + labs(x="digestive physiology", y="phylogenetic distance")
p = p + theme(panel.grid.major = element_blank()) + theme(panel.grid.minor = element_blank())
p = p + theme(axis.title.x = element_text(size= rel(1.2)), axis.title.y = element_text(size = rel(1.2)))
p = p + theme(axis.text.y = element_text(size = rel(1.3), colour = 'black'), axis.text.x = element_text(size = rel(1.3), color = 'black'))
p = p + theme(panel.background = element_rect(color="grey10"))
p

p = ggplot(PDwhole.data, aes(x=map$Sex, y=PDwhole.data$PD_whole_tree)) + geom_boxplot(fill="light blue", width = 0.6, outlier.size = 1.1)
p = p + labs(x="sex", y="phylogenetic distance")
p = p + theme(panel.grid.major = element_blank()) + theme(panel.grid.minor = element_blank())
p = p + theme(axis.title.x = element_text(size= rel(1.2)), axis.title.y = element_text(size = rel(1.2)))
p = p + theme(axis.text.y = element_text(size = rel(1.3), colour = 'black'), axis.text.x = element_text(size = rel(1.3), color = 'black'))
p = p + theme(panel.background = element_rect(color="grey10"))
p

p = ggplot(PDwhole.data, aes(x=map$SpeciesCommonName_Original, y=PDwhole.data$PD_whole_tree, fill=map$DigestivePhysiology)) + geom_boxplot(aes(width = 0.6, outlier.size = 1.1))
p = p + labs(x="primate species", y="phylogenetic distance") + xlim("Black-handed spider monkey", "Blue-eyed black lemur", "De Brazzas monkey", "Emperor tamarin", "Geoffroys tamarin", "White-faced saki", "Orangutan", "Western lowland gorilla")
p = p + theme(panel.grid.major = element_blank()) + theme(panel.grid.minor = element_blank())
p = p + theme(axis.title.x = element_text(size= rel(1.2)), axis.title.y = element_text(size = rel(1.2)))
p = p + theme(axis.text.y = element_text(size = rel(1.3), colour = 'black'), axis.text.x = element_text(size = rel(1.3), color = 'black', angle = 45, vjust = 1.01, hjust = 1))
p = p + theme(panel.background = element_rect(color="grey10"))
p = p + scale_fill_discrete(name="digestive\nphysiology")+ theme(legend.position="top")
p = p + theme(plot.margin = unit(c(0.1,0.2,0.2,1.8), 'cm'))
p


#
#
## Now looking at the alpha metrics averaged across the timepoints for each monkey.
#

setwd('~/Box Sync/knights_box/pmp_abx/alphadiv/')
alphameans.data = read.delim('alphadiv_adultmeans.txt', header = T, check.names = F)
colnames(alphameans.data) = c('sample_id', 'chao1', 'shannon','PD_whole_tree')
map = read.delim('../map_cz_alpha_adultmean.txt', header = T, row.names = 1, check.names = F)
sample.ids = intersect(rownames(map), alphameans.data$sample_id)
sample.ids = sort(sample.ids)
alphameans.data = alphameans.data[order(alphameans.data$sample_id),]
rownames(alphameans.data) = alphameans.data[,1]
#PDwhole.data[,1] = NULL
alphameans.data = alphameans.data[sample.ids,]
alphameans.data[,1] = NULL
map = map[sample.ids,]
dim(alphameans.data)
dim(map)

p = ggplot(alphameans.data, aes(x=map$Antibiotics_ever_used, y=alphameans.data$PD_whole_tree)) + geom_boxplot(fill="#009E73", width = 0.6, outlier.size = 1.1)
p = p + labs(x="Antibiotics ever used?", y="phylogenetic distance")
p = p + theme(panel.grid.major = element_blank()) + theme(panel.grid.minor = element_blank())
p = p + theme(axis.title.x = element_text(size= rel(1.2)), axis.title.y = element_text(size = rel(1.2)))
p = p + theme(axis.text.y = element_text(size = rel(1.3), colour = 'black'), axis.text.x = element_text(size = rel(1.3), color = 'black'))
p = p + theme(panel.background = element_rect(color="grey10"))
p

plot(alphameans.data$PD_whole_tree)
ks.test(alphameans.data$PD_whole_tree, pnorm, mean=mean(PDwhole.data$PD_whole_tree, sd=sd(PDwhole.data$PD_whole_tree)))
# Definitely not a normal distribution, using non-parametric test
PD.abx.mannpval = wilcox.test(alphameans.data$PD_whole_tree ~ map$Antibiotics_ever_used, data=alphameans.data)
PD.abx.mannpval


p = ggplot(alphameans.data, aes(x=map$DigestivePhysiology, y=alphameans.data$PD_whole_tree)) + geom_boxplot(fill="light blue", width = 0.6, outlier.size = 1.1)
p = p + labs(x="digestive physiology", y="mean phylogenetic distance")
p = p + theme(panel.grid.major = element_blank()) + theme(panel.grid.minor = element_blank())
p = p + theme(axis.title.x = element_text(size= rel(1.2)), axis.title.y = element_text(size = rel(1.2)))
p = p + theme(axis.text.y = element_text(size = rel(1.3), colour = 'black'), axis.text.x = element_text(size = rel(1.3), color = 'black'))
p = p + theme(panel.background = element_rect(color="grey10"))
p

p = ggplot(PDwhole.data, aes(x=map$Sex, y=alphameans.data$PD_whole_tree)) + geom_boxplot(fill="light blue", width = 0.6, outlier.size = 1.1)
p = p + labs(x="sex", y="mean phylogenetic distance")
p = p + theme(panel.grid.major = element_blank()) + theme(panel.grid.minor = element_blank())
p = p + theme(axis.title.x = element_text(size= rel(1.2)), axis.title.y = element_text(size = rel(1.2)))
p = p + theme(axis.text.y = element_text(size = rel(1.3), colour = 'black'), axis.text.x = element_text(size = rel(1.3), color = 'black'))
p = p + theme(panel.background = element_rect(color="grey10"))
p

p = ggplot(alphameans.data, aes(x=map$SpeciesCommonName_Original, y=alphameans.data$PD_whole_tree, fill=map$DigestivePhysiology)) + geom_boxplot(aes(width = 0.6, outlier.size = 1.1))
p = p + labs(x="primate species", y="phylogenetic distance") + xlim("Black-handed spider monkey", "Blue-eyed black lemur", "De Brazzas monkey", "Emperor tamarin", "Geoffroys tamarin", "White-faced saki", "Orangutan", "Western lowland gorilla")
p = p + theme(panel.grid.major = element_blank()) + theme(panel.grid.minor = element_blank())
p = p + theme(axis.title.x = element_text(size= rel(1.2)), axis.title.y = element_text(size = rel(1.2)))
p = p + theme(axis.text.y = element_text(size = rel(1.3), colour = 'black'), axis.text.x = element_text(size = rel(1.3), color = 'black', angle = 45, vjust = 1.01, hjust = 1))
p = p + theme(panel.background = element_rect(color="grey10"))
p = p + scale_fill_discrete(name="digestive\nphysiology")+ theme(legend.position="top")
p = p + theme(plot.margin = unit(c(0.1,0.2,0.2,1.8), 'cm'))
p


## Mixed model




