#Redo TUD distance export for the 60 phage in the Hatful paper
#read TUD for all 
all <- as.matrix(read.table('GitHub/tango/data/with_reverse_complement/all_phages_TUD_4_RC.tsv'))
sixty <- as.matrix(read.table('GitHub/tango/data/hatful60phageNames.txt'))
allnames <- sapply(rownames(all), function(x) strsplit(x, split="\\(")[[1]][1])
sixtynames <- sapply(sixty[,1], function(x) strsplit(x, split="\\(")[[1]][1])
sixtynames[33] <- 'PLot'
all_subset <- all[allnames %in% sixtynames, ]

d <- as.matrix(dist(all_subset))
d_all <- as.matrix(dist(all))
write.table(d, file='GitHub/tango/data/with_reverse_complement/Hatful_60_TUD_distance.tsv', sep='\t', quote=F)
write.table(all_subset, file='GitHub/tango/data/with_reverse_complement/Hatful_60_TUD.tsv', sep='\t', quote=F)
write.table(d_all, file='GitHub/tango/data/with_reverse_complement/all_TUD_distance.tsv', sep='\t', quote=F)


# get 60 subset for 3mers,5mers,6mers
dev2 <- as.matrix(read.table('~/GitHub/tango/data/kmer_counts/all_dev_k2.tsv'))
dev3 <- as.matrix(read.table('~/GitHub/tango/data/kmer_counts/all_dev_k3.tsv'))
dev4 <- as.matrix(read.table('~/GitHub/tango/data/kmer_counts/all_dev_k4.tsv'))
dev5 <- as.matrix(read.table('~/GitHub/tango/data/kmer_counts/all_dev_k5.tsv'))
dev6 <- as.matrix(read.table('~/GitHub/tango/data/kmer_counts/all_dev_k6.tsv'))

#write deviations
write.table(dev2[allnames %in% sixtynames,],file='~/GitHub/tango/data/kmer_counts/Hatful_60_subset/Hatful_60_dev_k2.tsv',sep='\t',quote=F)
write.table(dev3[allnames %in% sixtynames,],file='~/GitHub/tango/data/kmer_counts/Hatful_60_subset/Hatful_60_dev_k3.tsv',sep='\t',quote=F)
write.table(dev4[allnames %in% sixtynames,],file='~/GitHub/tango/data/kmer_counts/Hatful_60_subset/Hatful_60_dev_k4.tsv',sep='\t',quote=F)
write.table(dev5[allnames %in% sixtynames,],file='~/GitHub/tango/data/kmer_counts/Hatful_60_subset/Hatful_60_dev_k5.tsv',sep='\t',quote=F)
write.table(dev6[allnames %in% sixtynames,],file='~/GitHub/tango/data/kmer_counts/Hatful_60_subset/Hatful_60_dev_k5.tsv',sep='\t',quote=F)
write.table(dev6[allnames %in% sixtynames,],file='~/GitHub/tango/data/kmer_counts/Hatful_60_subset/Hatful_60_dev_k6.tsv',sep='\t',quote=F)
#write distances
write.table(as.matrix(dist(dev2[allnames %in% sixtynames,])),file='~/GitHub/tango/data/kmer_counts/Hatful_60_subset/Hatful_60_dist_k2.tsv',sep='\t',quote=F)
write.table(as.matrix(dist(dev3[allnames %in% sixtynames,])),file='~/GitHub/tango/data/kmer_counts/Hatful_60_subset/Hatful_60_dist_k3.tsv',sep='\t',quote=F)
write.table(as.matrix(dist(dev4[allnames %in% sixtynames,])),file='~/GitHub/tango/data/kmer_counts/Hatful_60_subset/Hatful_60_dist_k4.tsv',sep='\t',quote=F)
write.table(as.matrix(dist(dev5[allnames %in% sixtynames,])),file='~/GitHub/tango/data/kmer_counts/Hatful_60_subset/Hatful_60_dist_k5.tsv',sep='\t',quote=F)
write.table(as.matrix(dist(dev6[allnames %in% sixtynames,])),file='~/GitHub/tango/data/kmer_counts/Hatful_60_subset/Hatful_60_dist_k5.tsv',sep='\t',quote=F)
write.table(as.matrix(dist(dev6[allnames %in% sixtynames,])),file='~/GitHub/tango/data/kmer_counts/Hatful_60_subset/Hatful_60_dist_k6.tsv',sep='\t',quote=F)
