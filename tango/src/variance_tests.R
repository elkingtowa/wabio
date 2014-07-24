# TUD variation tests - within cluster variation of certain signals, etc
tud <- as.matrix(read.table('~/GitHub/tango/data/all_phages_TUD.tsv'))
tudf <- as.data.frame(read.table('~/GitHub/tango/data/all_phages_TUD.tsv'))
strsplit("224(E)",split="\\(")[[1]][2]

clusters <- sapply(rownames(tud), function(i) strsplit(x=strsplit(i,split="\\(")[[1]][2], split='\\)')[[1]][1])
cf <- factor(clusters)

tudf["cluster"] <- clusters

cnums <- sapply(levels(cf), function(x) sum(tudf["cluster"]==x))
b1var <- apply(tudf[tudf["cluster"]=="B1",1:256], 2, var)
b1mean <- apply(tudf[tudf["cluster"]=="B1",1:256], 2, mean)
b1median <- apply(tudf[tudf["cluster"]=="B1",1:256], 2, median)

a1var <- apply(tudf[tudf["cluster"]=="A1",1:256], 2, var)
a1mean <- apply(tudf[tudf["cluster"]=="A1",1:256], 2, mean)
a1median <- apply(tudf[tudf["cluster"]=="A1",1:256], 2, median)

allvar <- apply(tudf[,1:256], 2, var)

tcor <- cor(tud)
samplecor <- cor(t(tud))
pc <- prcomp(tud)
predict1 <- predict(pc)[,1]

# look at signals of individual tetras. Many appear bimodal. Interesting... 
library(diptest)
tudDip <- apply(tudf[1:256], 2, function(x) dip(x))
tudDipSorted <- sort(tudDip, decreasing=T)
#largest dip value is most bimodal - AGCC
plot(hist(tudf[,names(tudDipSorted[1])]))
#try splitting into two classes based on modes and see how clusters come out. 
#middle of AGCC looks to be at 0.8
AGCC1 <- tudf[tudf[,'AGCC']<0.8,"cluster"]
AGCC2 <- tudf[tudf[,'AGCC']>0.8,"cluster"]
# This makes an incredible split across clusters on just one value!!!
table(AGCC1)
table(AGCC2)

# look at top 10 most bimodal tetras and see if they can hcluster well.
tudSubset <- tudf[, names(tudDipSorted[1:10])]
plot(hclust(dist(tudSubset[1:100,])))

#is most bimodal the one with the most variance
tudVars <- sort(apply(tudf[1:256],2,var), decreasing=T)
#TCGA has the most variance and is also very bimodal. Try split on this
TCGA1 <- tudf[tudf[,'TCGA']<2.2,"cluster"]
TCGA2 <- tudf[tudf[,'TCGA']>2.2,"cluster"]
# Literally a perfect split between clusters.
table(TCGA1)
table(TCGA2)
# only singletons are in both tables
sum(names(table(TCGA1)) %in% names(table(TCGA2))) + sum(names(table(TCGA2)) %in% names(table(TCGA1)))

#B3 cluster can be differentiated just on the basis of TUD in GATC
length(tudf[tudf[,'GATC']>3.5,"cluster"])
sum(tudf[,"cluster"]=="B3")

#Try GAAG
plot(hist(tudf[,"GAAG"],breaks=50))
# GAAG is characteristic of B1 phage
plot(hist(tudf[tudf[,"cluster"]=="B1","GAAG"],breaks=50))

#look at each cluster and see what tetra has the least variance within that cluster.
clusterVar1 <- t(sapply(levels(cf), function(x) apply(tudf[tudf[,"cluster"]==x,1:256],2,var)))
# remove singletons
clusterVar = clusterVar1[rownames(clusterVar1)!="Singleton",]
# maximillly/minimally variant tetra in each cluster
maxVarCluster <- apply(clusterVar, 1, which.max)
minVarCluster <- apply(clusterVar, 1, which.min)
