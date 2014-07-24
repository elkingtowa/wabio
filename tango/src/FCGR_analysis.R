# analysis of tetra frequencies as FCGR
# import matrix
freqs <- as.data.frame(read.table('~/GitHub/tango/data/with_reverse_complement/FCGR_all_probability.tsv'))

clusters <- sapply(rownames(freqs), function(i) strsplit(x=strsplit(i,split="\\(")[[1]][2], split='\\)')[[1]][1])
cf <- factor(clusters)

freqs["cluster"] <- clusters
#sort by cluster
freqs = freqs[with(freqs,order(freqs[,'cluster'])),]

mycolor = colorRampPalette(c('white', 'red'))(n=3000)
heatmap(matrix(as.numeric(freqs[8,1:256]),nrow=16), Rowv=NA, Colv=NA,revC=T,col=myc, main=rownames(freqs[8,]))

for (i in 1:dim(freqs)[1]){
  png(file=paste('~/GitHub/tango/figures/with_reverse_complement/FCGR_heatmaps/',freqs[i,257],'_',rownames(freqs)[i],'.png',sep=''))
  heatmap(matrix(as.numeric(freqs[i,1:256]),nrow=16), Rowv=NA, Colv=NA,revC=T,col=myc, main=rownames(freqs[i,]))
  dev.off()
}

nameMat = matrix(colnames(freqs[,1:256]),nrow=16)
write.table(nameMat,file='~/GitHub/tango/figures/with_reverse_complement/FCGR_heatmaps/heatmap_legend.tsv',sep='\t')

  