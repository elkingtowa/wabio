# doing Naive Bayes on tetranucleotide probabilities with multinomial distributions. 
# NEW CODE 2014-07-01
setwd("~/GitHub/tango")
#read phage names and clusters, tetranucleotides for all
phageData <- as.data.frame(read.table('kmer_analysis/examples/PhagesDB_Data.txt',skip=1))
allK4 <- as.matrix(read.table('data/kmer_counts/all_phage/all_freq_k4.csv',sep=',',header=T))
allK4p <- allK4 / rowSums(allK4)

# try for a cluster O - catdawg
baseName = 'data/kmer_counts/indiviudal_windows/'
ending = '_500_250_k4.csv'

catdawg = as.matrix(read.table(paste(baseName, 'Catdawg-O', ending, sep=''), sep=',', header=T))
selfSimilarity <- apply(windows, 1, function(x) dmultinom(x,prob=allK4p['Catdawg-O',], log=T))
plot(1:288, selfSimilarity, type='o')
plot(250:288, selfSimilarity[250:288], type='o')

#self similarity is lowest in window 4000:4500
compare4k <- sort(apply(allK4p, 1, function(x) dmultinom(windows['4000:4500',], prob=x, log=T)), decreasing=T)

#patience-U has high similarity with window 4000:4500. look for similarity across this 
patienceSimilarity <- apply(windows, 1, function(x) dmultinom(x,prob=allK4p['Patience-U',], log=T))

plot(1:50, selfSimilarity[1:50], type='o')
lines(patienceSimilarity[1:50], col='red', type='o')
lines(apply(windows, 1, function(x) dmultinom(x,prob=allK4p['Gaia-Singleton',], log=T))[1:50], col='blue',type='o')
lines(apply(windows, 1, function(x) dmultinom(x,prob=allK4p['Hawkeye-D2',], log=T))[1:50], col='green',type='o')
lines(apply(windows, 1, function(x) dmultinom(x,prob=allK4p['Nilo-R',], log=T))[1:50], col='orange',type='o')

#look at patience?
patience = as.matrix(read.table(paste(baseName, 'Patience-U', ending, sep=''), sep=',', header=T))
plot(apply(windows, 1, function(x) dmultinom(x,prob=allK4p['Patience-U',], log=T))[210:240],type='o')

# NEW CODE 2014-05-18
# calclulations similar to TDI. In sliding window, get probability sequence came from cluster. 
#get phage data 
pdata <- as.data.frame(read.table('~/GitHub/tango/data/phagesDB/sequenced_phage_metadata_simple.txt',row.names=1, header=T))
# read data for various nucleotides
mers1 <- as.matrix(read.table('~/GitHub/tango/data/kmer_counts/all_freq_k1.tsv'))
mers4 <- as.matrix(read.table('~/GitHub/tango/data/kmer_counts/all_freq_k4.tsv'))
#normalize each to fraction
mers1p <- mers1/rowSums(mers1)
mers4p <- mers4/rowSums(mers4)

#get data for a sliding window I have.
jd <- as.matrix(read.table('~/GitHub/tango/data/kmer_counts/indiviudal/JoeDirt(L1)_2000_2000_k4_freq.txt'))
#sums of counts in each cluster
subclusters <- levels(pdata[,'subcluster'])
subclusterNums <- table(pdata[,'subcluster'])
subclustermers4 <- matrix(ncol=dim(mers4)[2])
for (subcluster in subclusters){
  if (subclusterNums[subcluster] > 1){
    subclustermers4 <- rbind(subclustermers4, colSums(mers4[pdata[,'subcluster']==subcluster,]))
  }
  else subclustermers4 <- rbind(subclustermers4, mers4[pdata[,'subcluster']==subcluster,])
}
subclustermers4 <- subclustermers4[2:(length(subclusters)+1),]
rownames(subclustermers4) <- subclusters
subclustermers4p <- subclustermers4 / rowSums(subclustermers4)

#Likelihood of each section of each cluster being the origin
likelihood.cluster = apply(jd, 1, function(y) apply(subclustermers4p, 1, function(x) dmultinom(y, size=sum(y),prob=x,log=T)))
#plot likelihoods of various clusters
num = 48
cp <- rainbow(num)
plot(likelihood.cluster[1,],type='o',col=cp[1],ylim=c(-1500,max(likelihood.cluster)))
for (i in 1:48){
  lines(likelihood.cluster[i,],type='o',col=cp[i])
}
legend(1, 0.15, rownames(likelihood.cluster), cex=0.5, col=cp, lty=1)
# find max in each window
rownames(likelihood.cluster)[apply(likelihood.cluster,2,which.max)]

# nicer figure to use in report
plot((1:37)*2000, likelihood.cluster["L1",1:37], ,ylim=c(-950,-600),type='o', col='blue',xlab="Genomic position", ylab="log-likelihood", main="4-mer usage in JoeDirt compared to cluster averages, 2000bp sliding window")
lines( (1:37)*2000, likelihood.cluster["J",1:37], type='o', col='red')
legend(1, -710, c("p(cluster L1|JoeDirt)", "p(cluster J|JoeDirt)") , cex=0.8, col=c("blue","red"), lty=1)

# REPEAT FOR BONGO(M1)
#get data for a sliding window I have.
bongo <- as.matrix(read.table('~/GitHub/tango/data/kmer_counts/indiviudal/Bongo(M1)_5000_2500_k4_freq.txt'))

#Likelihood of each section of each cluster being the origin
likelihood.cluster.bongo = apply(bongo, 1, function(y) apply(subclustermers4p, 1, function(x) dmultinom(y, size=sum(y),prob=x,log=T)))
max.likelihood.bongo <- apply(likelihood.cluster.bongo,2,which.max)
rownames(likelihood.cluster.bongo)[max.likelihood.bongo]

#plot likelihoods of various clusters
num = 48
cp <- rainbow(num)
plot(likelihood.cluster.bongo[1,],type='o',col=cp[1],ylim=c(-1500,max(likelihood.cluster)))
for (i in 1:48){
  lines(likelihood.cluster.bongo[i,],type='o',col=cp[i])
}
#legend(1, 0.15, rownames(likelihood.cluster), cex=0.5, col=cp, lty=1)



#read in 4-mer probabilities and frequencies for all phage
probs <- as.data.frame(read.table('~/GitHub/tango/data/with_reverse_complement/FCGR_all_probability.tsv'))
freqs <- as.data.frame(read.table('~/GitHub/tango/data/with_reverse_complement/FCGR_all_frequency.tsv'))

#have to compute with integers (freqs)
#compute with JoeDirt against all others
compare <- freqs["JoeDirt(L1)",]
compared <- apply(probs, 1, function(x) dmultinom(compare, size=sum(compare),prob=x,log=T))
compared <- sort(compared,decreasing=T)
head(compared)
# Other L1s, then L2s, then L3 are the most likely. makes sense. 

#compare the last section of JoeDirt agains all others 
jd <- as.matrix(read.table('~/GitHub/tango/data/window_kmer_count//JoeDirt(L1)_2000_2000_k4_freq.txt'))
compare_start = apply(jd[1:32,], 2, sum)
compare_end = apply(jd[36:38,], 2, sum)
#compare_end=jd[38,]
compared_start = apply(probs, 1, function(x) dmultinom(compare_start, size=sum(compare_start),prob=x,log=T))
compared_end = apply(probs, 1, function(x) dmultinom(compare_end, size=sum(compare_end),prob=x,log=T))
compared_start <- sort(compared_start,decreasing=T)
compared_end <- sort(compared_end,decreasing=T)
head(compared_start,n=15)
head(compared_end,n=30)

#compare last of jd against all windows
probs <- as.data.frame(read.table('~/GitHub/tango/data/window_kmer_count/all_prob_5000_k4.tsv'))
freqs <- as.data.frame(read.table('~/GitHub/tango/data/window_kmer_count//all_freq_5000_k4.tsv'))
compare_end = freqs["JoeDirt(L1)_70000:75000",]
compared_end = apply(probs, 1, function(x) dmultinom(compare_end, size=sum(compare_end),prob=x,log=T))
compared_end <- sort(compared_end,decreasing=T)
head(compared_end,n=100)
