#read data table
match_table <- as.matrix(read.table('blast/blast_all_filtered.csv', sep=',', header=T))
seq_phage_clusters <- as.matrix(read.table('phages/sequenced_phage_clusters.csv', sep=',',header=T))
all_phage_clusters <- as.matrix(read.table('phages/all_phage_clusters.csv', sep=',',header=T))

#plot of distribution of clusters across significant matches
cluster_nums <- table(match_table[,'cluster'])
png('figures/match_cluster_distribution.png',height=600, width=800)
barplot(cluster_nums, main='Distribution of all CRISPR matches',xlab='cluster',ylab='# matches')
dev.off()

#null distribution of clusters - sequenced phages
seq_cluster_nums <- table(seq_phage_clusters[,2])
png('figures/sequenced_phage_cluster_distribution.png',height=600, width=800)
barplot(seq_cluster_nums, main='Distribution of clusters in sequenced phages',xlab='cluster',ylab='# matches')
dev.off()
#all phages
all_cluster_nums <- table(all_phage_clusters[,2])[c((1:41),(43:49))]
png('figures/all_phage_cluster_distribution.png',height=600, width=800)
barplot(all_cluster_nums, main='Distribution of clusters in all phages (not none, unclustered)',xlab='cluster',ylab='# matches')
dev.off()

#plot of distribution of crispr match sequences
myco_nums <- table(match_table[,'qseqid'])
names(myco_nums) <- 1:34
png('figures/match_qseqid_distribution.png',height=600, width=800)
barplot(myco_nums, main='Distribution of all CRISPR matches',xlab='qseqid',ylab='# matches', )
dev.off()

#comparison to theroretical null distribution
cluster_names= names(all_cluster_nums)
m<- matrix(0, nrow=48, ncol=1)
rownames(m) <- cluster_names
cluster_freq = cluster_nums/sum(cluster_nums)
for (i in 1:48){
  if (cluster_names[i] %in% names(cluster_freq)){
    m[i] <- cluster_freq[rownames(m)[i]]
  }
}
cluster_mat= cbind(m,seq_cluster_nums/(sum(seq_cluster_nums)),all_cluster_nums/sum(all_cluster_nums))
colnames(cluster_mat) <- c('CRISPR matches', 'sequenced phages','all phages')

png('figures/cluster_distribution_comparison.png',width=2000, height=800)
barplot(t(cluster_mat), col=c("green","red","darkblue"),legend=colnames(cluster_mat), main='Distribution of clusters', xlab='Cluster',ylab='Frequency',beside=T)
dev.off()
