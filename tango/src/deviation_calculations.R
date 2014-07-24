# analysis of palindromes in phage genomes 
#get phage data 
pdata <- as.data.frame(read.table('~/GitHub/tango/data/phagesDB/sequenced_phage_metadata_simple.txt',row.names=1, header=T))
# read data for various nucleotides
mers1 <- as.matrix(read.table('~/GitHub/tango/data/kmer_counts/all_freq_k1.tsv'))
mers2 <- as.matrix(read.table('~/GitHub/tango/data/kmer_counts/all_freq_k2.tsv'))
mers3 <- as.matrix(read.table('~/GitHub/tango/data/kmer_counts/all_freq_k3.tsv'))
mers4 <- as.matrix(read.table('~/GitHub/tango/data/kmer_counts/all_freq_k4.tsv'))
mers5 <- as.matrix(read.table('~/GitHub/tango/data/kmer_counts/all_freq_k5.tsv'))
mers6 <- as.matrix(read.table('~/GitHub/tango/data/kmer_counts/all_freq_k6.tsv'))
#normalize each to fraction
mers1p <- mers1/rowSums(mers1)
mers2p <- mers2/rowSums(mers2)
mers3p <- mers3/rowSums(mers3)
mers4p <- mers4/rowSums(mers4)
mers5p <- mers5/rowSums(mers5)
mers6p <- mers6/rowSums(mers6)

#normalize to zero order markov
usageDeviation <- function (observedCounts, baseFracitons){
  expected <- observedCounts
  for (i in 1:dim(observedCounts)[1]){
    size = sum(observedCounts[i,])
    for (j in 1:dim(observedCounts)[2]){
      bases <- table(strsplit(colnames(observedCounts)[j], split=''))
      #print(bases)
      #print(sapply(1:length(bases), function(x) baseFractions[i,names(bases)[x]]**bases[x]))
      #print(size)
      exp <- prod(sapply(1:length(bases), function(x) baseFractions[i,names(bases)[x]]**bases[x]))*size
      #print(exp)
      expected[i,j] <- exp
    }
  }
  #print(head(expected))
  return(observedCounts/expected)
}

dev2 <- usageDeviation(mers2, mers1p)
dev3 <- usageDeviation(mers3, mers1p)
dev4 <- usageDeviation(mers4, mers1p)
dev5 <- usageDeviation(mers5, mers1p)
dev6 <- usageDeviation(mers6, mers1p)

write.table(dev2,file='~/GitHub/tango/data/kmer_counts/all_dev_k2.tsv',sep='\t',quote=F)
write.table(dev3,file='~/GitHub/tango/data/kmer_counts/all_dev_k3.tsv',sep='\t',quote=F)
write.table(dev4,file='~/GitHub/tango/data/kmer_counts/all_dev_k4.tsv',sep='\t',quote=F)
write.table(dev5,file='~/GitHub/tango/data/kmer_counts/all_dev_k5.tsv',sep='\t',quote=F)
write.table(dev6,file='~/GitHub/tango/data/kmer_counts/all_dev_k6.tsv',sep='\t',quote=F)

#5mers with "GATC"
gatc5 <- mers5p[,grep("GATC", colnames(mers5p))]
#GGATC and GATCC are most descriminatory
#6mers with "GATC"
gatc6 <- mers6p[,grep("GATC", colnames(mers6p))]
# GGATCC are most different
ggatcc <- sort(mers6p[,"GGATCC"], decreasing=T)
hist(mers6p[,"GGATCC"],breaks=30,main="GGATCC counts",xlab="Frequency in genome")
#plot 6mers with GATC
pdf(file='~/GitHub/tango/figures/with_reverse_complement/motif_plots/GATC6mers.pdf')
par(mfrow=c(4,4),mai=c(0.4,0.25,0.5,0.1))
for (i in 1:48){
  hist(gatc6[,i],breaks=40,main=colnames(gatc6)[i])
}
dev.off()