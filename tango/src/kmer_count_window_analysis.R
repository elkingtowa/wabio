# analysis of tetranucleotide probability in windows across the genome
#look at 5kb windows in L1 genomes

jd <- as.matrix(read.table('~/GitHub/tango/data/window_kmer_count//JoeDirt(L1)_5000_2500.txt'))
jd2 <- as.matrix(read.table('~/GitHub/tango/data/window_kmer_count//JoeDirt(L1)_5000_2500_k2.txt'))
jd3 <- as.matrix(read.table('~/GitHub/tango/data/window_kmer_count//JoeDirt(L1)_5000_2500_k3.txt'))


#compute averages and variances
jdMean <- sort(apply(jd, 2, mean),decreasing=T)
jdVar <- sort(apply(jd[20:29,], 2, var),decreasing=T)
jdTopEnd <- jd[,order(jd[29,],decreasing=T)]

#plot first few by mean and var
num = 15
cp <- rainbow(num)
#by mean
#plot(jd[,names(jdMean[1])], type='o', col=cp[1], ylim=c(0,0.02))
#for (i in 2:num){
#  lines(jd[,names(jdMean[i])], type='o', col=cp[i])
#}
#legend(1, 0.02, names(jdMean[1:num]), cex=0.3, col=cp, lty=1)
#by var
plot(x=(2:30 *2500)[20:29], y=jd[20:29,names(jdVar[1])], type='o', col=cp[1], ylim=c(0,0.02),xlab='genomic position',ylab='relative frequency',main='most variable 4-mers JoeDirt')
for (i in 2:num){
  lines(x=(2:30 *2500)[20:29], y=jd[20:29,names(jdVar[i])], type='o', col=cp[i])
}
legend(52500, 0.02, names(jdVar[1:num]), cex=0.5, col=cp, lty=1)
#by top at last point
plot(jdTopEnd[20:29,1], type='o', col=cp[1], ylim=c(0,0.02))
for (i in 2:num){
  lines(jdTopEnd[20:29,i], type='o', col=cp[i])
}
legend(1, 0.02, colnames(jdTopEnd[,1:num]), cex=0.5, col=cp, lty=1)

for (i in 1:256){
  plot(jdTopEnd[,i],type='o',col='red',main=colnames(jdTopEnd)[i])
  par(ask=T)
}

#plot dimers
num = 16
cp <- rainbow(num)
plot(jd2[,1],type='o',col=cp[1],ylim=c(0,0.15))
for (i in 2:16){
  lines(jd2[,i],type='o',col=cp[i])
}
legend(1, 0.15, colnames(jd2), cex=0.5, col=cp, lty=1)

for (i in 1:16){
  plot(jd2[,i],type='o',col='red',main=colnames(jd2)[i])
  par(ask=T)
}

#plot trimers
num = 16
cp <- rainbow(num)
for (i in 1:64){
  plot(jd3[,i],type='o',col='red',main=colnames(jd3)[i])
  par(ask=T)
}