# plot counts of GATC in a sliding window in B1 B2 B3 phage 
install.packages('extrafont')
library(extrafont)
font_import()
# read data
b1 <- as.matrix(read.table('GitHub/tango/data/without_reverse_complement/GATC_motif_B1.tsv',row.names=1, header=F, sep=','))
b2 <- as.matrix(read.table('GitHub/tango/data/without_reverse_complement/GATC_motif_B2.tsv',row.names=1, header=F, sep=','))
b3 <- as.matrix(read.table('GitHub/tango/data/without_reverse_complement/GATC_motif_B3.tsv',row.names=1, header=F, sep=','))

b1.mean <- apply(b1[2:83,], 2, function(x) mean(x,na.rm=T))
b1.var <- apply(b1[2:83,], 2, function(x) var(x,na.rm=T))
b2.mean <- apply(b2[2:11,],2,function(x) mean(x,na.rm=T))
b2.var <- apply(b2[2:11,], 2, function(x) var(x,na.rm=T))
b3.mean <- apply(b3[2:14,], 2, function(x) mean(x,na.rm=T))
b3.var <- apply(b3[2:14,], 2, function(x) var(x,na.rm=T))

col1=rgb(247/255,144/255,41/255,1)
col2=rgb(181/255,207/255,114/255,1)
col3=rgb(67/255,105/255,189/255,1)
col3a=rgb(67/255,105/255,189/255,0.5)
plot(b1[1,], b1.mean, type='l', ylim=c(0, max(b3.mean)+10), lwd=3, col=col1,
     main="GATC counts in 5kb sliding window", xlab="Genomic position", ylab="GATC counts",
     )
lines(b2[1,], b2.mean, type='l', lwd=3, col=col2)
lines(b3[1,], b3.mean, type='l', lwd=3, col=col3)
legend(8000, 83, c("B1 mean","B2 mean","B3 mean"), cex=1, col=c(col1, col2, col3), lty=1, lwd=3,horiz=T)


  # SD BARS LOOK LIKE SHIT
# segments(b1[1,],b1.mean-b1.var,b1[1,],b1.mean+b1.var,col='gray50')
# epsilon <- 300
# segments(b1[1,]-epsilon,b1.mean-b1.var,b1[1,]+epsilon,b1.mean-b1.var,col='gray50')
# segments(b1[1,]-epsilon,b1.mean+b1.var,b1[1,]+epsilon,b1.mean+b1.var,col='gray50')
# 
# segments(b2[1,],b2.mean-b2.var,b2[1,],b2.mean+b2.var,col='gray50')
# epsilon <- 300
# segments(b2[1,]-epsilon,b2.mean-b2.var,b2[1,]+epsilon,b2.mean-b2.var,col='gray50')
# segments(b2[1,]-epsilon,b2.mean+b2.var,b2[1,]+epsilon,b2.mean+b2.var,col='gray50')
# 
# segments(b3[1,],b3.mean-b3.var,b3[1,],b3.mean+b3.var,col='gray50')
# epsilon <- 300
# segments(b3[1,]-epsilon,b3.mean-b3.var,b3[1,]+epsilon,b3.mean-b3.var,col='gray50')
# segments(b3[1,]-epsilon,b3.mean+b3.var,b3[1,]+epsilon,b3.mean+b3.var,col='gray50')

# plot new Histogram of GATC and GGATCC for presentation
# THIS IS WITHOUT REVERSE COMPLEMENT
dev4 <- as.data.frame(read.table('~/GitHub/tango/data/kmer_counts/all_dev_k4.tsv'))
dev6 <- as.data.frame(read.table('~/GitHub/tango/data/kmer_counts/all_dev_k6.tsv'))
clusters <- sapply(rownames(dev4), function(i) strsplit(x=strsplit(i,split="\\(")[[1]][2], split='\\)')[[1]][1])
hist(dev4[,"GATC"], breaks=50)
hist(dev6[clusters!='B3',"GGATCC"], breaks=50)
library(ggplot2)
#plot GATC
ggplot(as.data.frame(dev4[,'GATC']), aes(x=dev4[,'GATC'])) + geom_histogram(binwidth=.1, colour="black", fill=col3a) +
  ggtitle("GATC usage deviation") +xlab('Tetranucleotide usage deviation') + ylab('Number of phage')
#plot GGATCC
ggplot(as.data.frame(dev6[,'GGATCC']), aes(x=dev6[,'GGATCC'])) + geom_histogram(binwidth=.1, colour="black", fill=col3a) +
  ggtitle("GGATCC usage deviation") +xlab('Hexanucleotide usage deviation') + ylab('Number of phage')
