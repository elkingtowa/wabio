# CODE FOR EDDIE
tud <- as.matrix(read.table('~/GitHub/tango/data/mycobacteria_TUD.tsv'))
d <- dist(tud, method='euclidean',)
fit1 <- hclust(d,method='ward.D')
plot(fit1)

#clustering 1-similarity matrix from pairwise alignment
aln <- as.matrix(read.table('~/GitHub/tango/data/pairwise_alignment_distance_60.tsv'))
d <- dist(aln, method='euclidean')
fit1 <- hclust(d,method='ward.D')
plot(fit1)

# clustering on TUD distance 
tud <- as.matrix(read.table('~/GitHub/tango/data/all_phages_TUD_3.tsv'))
d <- dist(tud, method='euclidean',)
# try different clustering methods
fit1 <- hclust(d,method='ward.D')
fit2 <- hclust(d,method='ward.D2')
fit3 <- hclust(d,method='complete')
fit4 <- hclust(d,method='single')
par(mfrow=c(2,2))
plot(fit1)
plot(fit2)
plot(fit3)
plot(fit4)

# working with ape for tree building? 
library(ape)
njtree <- nj(d)

# testing things like PCA and MDS to dimension reduce/cluster
fit <- cmdscale(d, eig=T, k=2)
fit
# plot solution 
x <- fit$points[,1]
y <- fit$points[,2]

clusters2 <- sapply(rownames(tud), function (x) strsplit(x,split="\\(")[[1]][2])
clusters <- sapply(clusters2, function (x) strsplit(x,split="\\)")[[1]][1])

plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", 
     main="Metric  MDS",   type="n")
text(x, y, labels = clusters, cex=.7)

#nonmetric MDS
library(MASS)
fit<-isoMDS(d,k=2)
x <- fit$points[,1]
y <- fit$points[,2]

clusters2 <- sapply(rownames(tud), function (x) strsplit(x,split="\\(")[[1]][2])
clusters <- sapply(clusters2, function (x) strsplit(x,split="\\)")[[1]][1])

plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", 
     main="Metric  MDS",   type="n")
text(x, y, labels = clusters, cex=.7)

#trying different distance metrics
#manhattan distance
d <- dist(tud, method='manhattan')
fit1 <- hclust(d,method='ward.D')
fit2 <- hclust(d,method='ward.D2')
fit3 <- hclust(d,method='complete')
fit4 <- hclust(d,method='single')
par(mfrow=c(2,2))
plot(fit1)
plot(fit2)
plot(fit3)
plot(fit4)
#dev.off()
#MDS
fit <- cmdscale(d, eig=T, k=2)
# plot solution 
x <- fit$points[,1]
y <- fit$points[,2]
clusters2 <- sapply(rownames(tud), function (x) strsplit(x,split="\\(")[[1]][2])
clusters <- sapply(clusters2, function (x) strsplit(x,split="\\)")[[1]][1])
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", 
     main="Metric  MDS",   type="n")
text(x, y, labels = clusters, cex=.7)
#manhattan verdict:
# hclustering seems to perform better. MDS is not better or even worse. 
# plot of clustering saved, MDS no. 

#minkowski distance
d <- dist(tud, method='minkowski')
fit1 <- hclust(d,method='ward.D')
fit2 <- hclust(d,method='ward.D2')
fit3 <- hclust(d,method='complete')
fit4 <- hclust(d,method='single')
par(mfrow=c(2,2))
plot(fit1)
plot(fit2)
plot(fit3)
plot(fit4)
dev.off()
#MDS
fit <- cmdscale(d, eig=T, k=2)
# plot solution 
x <- fit$points[,1]
y <- fit$points[,2]
clusters2 <- sapply(rownames(tud), function (x) strsplit(x,split="\\(")[[1]][2])
clusters <- sapply(clusters2, function (x) strsplit(x,split="\\)")[[1]][1])
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", 
     main="Metric  MDS",   type="n")
text(x, y, labels = clusters, cex=.7)
# verdict minkowski 
# hclust just as good as manhattan. MDS better than manhattan. 
# plot saved for both. 

#correlation as distance (pearson)
d <- cor(t(tud))
fit1 <- hclust(d,method='ward.D')
fit2 <- hclust(d,method='ward.D2')
fit3 <- hclust(d,method='complete')
fit4 <- hclust(d,method='single')
par(mfrow=c(2,2))
plot(fit1)
plot(fit2)
plot(fit3)
plot(fit4)
dev.off()
#MDS
fit <- cmdscale(d, eig=T, k=2)
# plot solution 
x <- fit$points[,1]
y <- fit$points[,2]
clusters2 <- sapply(rownames(tud), function (x) strsplit(x,split="\\(")[[1]][2])
clusters <- sapply(clusters2, function (x) strsplit(x,split="\\)")[[1]][1])
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", 
     main="Metric  MDS",   type="n")
text(x, y, labels = clusters, cex=.7)
