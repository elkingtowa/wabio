# repeat some of the analysis on the data with reverse complements added
tud <- as.data.frame(read.table('~/GitHub/tango/data/all_phages_TUD.tsv'))
tudrc <- as.data.frame(read.table('~/GitHub/tango/data/with_reverse_complement/all_phages_TUD_4_RC.tsv'))
#tudrc <- as.data.frame(read.table('~/projects/tango/data/with_reverse_complement/all_phages_TUD_4_RC.tsv'))

d<- dist(tudrc)
plot(hist(d))
#holy shit it looks even better this way. 

clusters <- sapply(rownames(tudrc), function(i) strsplit(x=strsplit(i,split="\\(")[[1]][2], split='\\)')[[1]][1])
cf <- factor(clusters)
tudrc["cluster"] <- clusters

tudVars <- sort(apply(tudrc[1:256],2,var), decreasing=T)
#TCGA has the most variance and is also very bimodal. Try split on this
TCGA1 <- tudrc[tudrc[,'TCGA']<2.2,"cluster"]
TCGA2 <- tudrc[tudrc[,'TCGA']>2.2,"cluster"]
# Literally a perfect split between clusters.
table(TCGA1)
table(TCGA2)
# only singletons are in both tables
sum(names(table(TCGA1)) %in% names(table(TCGA2))) + sum(names(table(TCGA2)) %in% names(table(TCGA1)))
#look at top part of TCGA
TCGA3 <- tudrc[tudrc[,'TCGA']<3.8,"cluster"]
TCGA4 <- tudrc[tudrc[,'TCGA']>3.8,"cluster"]
# Literally a perfect split between clusters.
table(TCGA3)
table(TCGA4)


#B3 cluster can be differentiated just on the basis of TUD in GATC
length(tudrc[tudrc[,'GATC']>3.5,"cluster"])
sum(tudrc[,"cluster"]=="B3")

#Try GAAG
plot(hist(tudrc[,"GAAG"],breaks=50))
# GAAG is characteristic of B1 phage
plot(hist(tudrc[tudrc[,"cluster"]=="B1","GAAG"],breaks=50))
GAAG1 <- tudrc[tudrc[,'GAAG']<0.49,"cluster"]
GAAG2 <- tudrc[tudrc[,'GAAG']>0.49,"cluster"]
#only the Ps contaminate this split, otherwise perfect
table(GAAG1)
table(GAAG2)

#Try CCTA
plot(hist(tudrc[,"CCTA"],breaks=50))
# CCTA 
#plot(hist(tudrc[tudrc[,"cluster"]=="B1","GAAG"],breaks=50))
CCTA1 <- tudrc[tudrc[,'CCTA']<0.55,"cluster"]
CCTA2 <- tudrc[tudrc[,'CCTA']>0.55,"cluster"]
CCTA3 <- tudrc[tudrc[,'CCTA']>0.7,"cluster"]
table(CCTA1)
table(CCTA2) #Ls and Rs come out here
table(CCTA3) #Rs isolated further with this split

#Try CTAA
plot(hist(tudrc[,"CTAA"],breaks=50)) 
#plot(hist(tudrc[tudrc[,"cluster"]=="B1","GAAG"],breaks=50))
CCTA1 <- tudrc[tudrc[,'CTAA']<1.0,"cluster"]
CCTA2 <- tudrc[tudrc[,'CTAA']>1.0,"cluster"]
table(CCTA1)
table(CCTA2) #Ls and come out here

# are tetras normally distributed within a given cluster?
# check 
qqnorm(tudrc[tudrc[,"cluster"]=="B1",1])
plot(hist(tudrc[tudrc[,"cluster"]=="B1",1]))

#For each cluster, for each TUD, how many other clusters are in the range of those values
clusterMaxVals <- t(sapply(levels(cf), function(x) apply(tudrc[tudrc[,"cluster"]==x,1:256],2,max)))
clusterMinVals <- t(sapply(levels(cf), function(x) apply(tudrc[tudrc[,"cluster"]==x,1:256],2,min)))

#dip statistic
library(diptest)
tudDip <- sort(apply(tudrc[1:256], 2, function(x) dip(x)),decreasing=T)
tudDipO <- order(apply(tudrc[1:256], 2, function(x) dip(x)),decreasing=T)

#sort dataframe starting with largest variance 
tudVarsO <- order(apply(tudrc[1:256],2,var), decreasing=T)
tudrc2<- tudrc[,tudVarsO]
tudrc3<- tudrc[,tudDipO]
# try plotting lots of histograms on same plot
pdf(file='figures/with_reverse_complement/motif_plots/all_motif_histograms_by_var.pdf')
par(mfrow=c(4,4),mai=c(0.4,0.25,0.5,0.1))
for (j in 1:16){
    for (i in ((j-1)*16+1):(j*16)){
    hist(tudrc2[,i],breaks=40,main=colnames(tudrc2[i]),xlab='TUD')
  }
}
dev.off()

pdf(file='figures/with_reverse_complement/motif_plots/all_motif_histograms_by_dip.pdf')
par(mfrow=c(4,4),mai=c(0.4,0.25,0.5,0.1))
for (j in 1:16){
  for (i in ((j-1)*16+1):(j*16)){
    hist(tudrc3[,i],breaks=40,main=colnames(tudrc3[i]),xlab='TUD')
  }
}
dev.off()
