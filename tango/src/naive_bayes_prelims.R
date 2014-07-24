# Building things for the naive bayes classifier - Cluster version. 
tudrc <- as.data.frame(read.table('~/GitHub/tango/data/with_reverse_complement//all_phages_TUD_4_RC.tsv'))
# renormalize everything to Log2 so values aren't constrained
tudrcLog <- data.frame(apply(tudrc, 1:2, log2))
clusters <- sapply(rownames(tudrc), function(i) strsplit(x=strsplit(i,split="\\(")[[1]][2], split='\\)')[[1]][1])
cf <- factor(clusters)
tudrcLog["cluster"] <- clusters


#mean and variance of each tetra within a cluster
clusterMean1 <- t(sapply(levels(cf), function(x) apply(tudrc[tudrc[,"cluster"]==x,1:256],2,mean)))
clusterVar1 <- t(sapply(levels(cf), function(x) apply(tudrc[tudrc[,"cluster"]==x,1:256],2,var)))
# remove singletons
clusterMean = clusterMean1[rownames(clusterVar1)!="Singleton",]
clusterVar = clusterVar1[rownames(clusterVar1)!="Singleton",]

write.table(clusterMean,file='~/GitHub/tango/data/with_reverse_complement/NB_clusterMean.tsv', sep='\t',quote=F)
write.table(clusterVar,file='~/GitHub/tango/data/with_reverse_complement/NB_clusterVar.tsv', sep='\t',quote=F)
write.table(table(clusters[clusters!="Singleton"]),file='~/GitHub/tango/data/with_reverse_complement/NB_clusterNums.tsv', sep='\t',quote=F,row.names=F)
