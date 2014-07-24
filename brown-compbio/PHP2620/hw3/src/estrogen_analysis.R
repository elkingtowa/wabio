#processing of estrogen data
library(limma)
library(estrogen)
library(affy)
library(ggplot2)

#loading data
datadir <- file.path(find.package("estrogen"),"extdata")
targets<-readTargets("phenoData.txt",path=datadir,sep="",row.names="filename")
ab<-ReadAffy(filenames=targets$filename,celfile.path=datadir,phenoData=targets)

#visualization of raw images 
for (i in 1:length(ab@phenoData@data$filename)){
  png(paste(ab@phenoData@data$filename[i],".png", sep=''))
  image(ab[,i])
  dev.off()
}

#pre norm MvA plot
exp <- exprs(ab)
#control
png("mva1.png")
ma.plot(rowMeans(log2(exp[,c(1,2)])),log2(exp[,1])-log2(exp[,2]),cex=1)
title("Pre-norm replicate, sample 1v2")
dev.off()
#different samples
png("mva2.png")
ma.plot(rowMeans(log2(exp[,c(1,8)])),log2(exp[,1])-log2(exp[,8]),cex=1)
title("Pre-norm difference, sample 1v8")
dev.off()

#boxplot of intensities across samples
png('raw_boxplot.png', height=350, width=525)
boxplot(ab, names=c('l.10.1','l.10.2','h.10.1','h.10.2','l.48.1','l.48.2','h.48.1','h.48.2'),
        main='Raw probe intensity across samples',
        xlab='Sample',ylab='log2(intensity)', col='red')
dev.off()

#do RMA and convert to expression set
eset<-rma(ab)
x <- exprs(eset)

#redo mVa plot and check differences
#control
png("mva3.png")
ma.plot(rowMeans(log2(x[,c(1,2)])),log2(x[,1])-log2(x[,2]),cex=1)
title("Post-norm replicate, sample 1v2")
dev.off()
#different samples
png("mva4.png")
ma.plot(rowMeans(log2(x[,c(1,8)])),log2(x[,1])-log2(x[,8]),cex=1)
title("Post-norm difference, sample 1v8")
dev.off()

#redo boxplot and check if equal
png('eset_boxplot.png', height=350, width=525)
boxplot(x, names=c('l.10.1','l.10.2','h.10.1','h.10.2','l.48.1','l.48.2','h.48.1','h.48.2'),
        main='Normalized expression across samples',
        xlab='Sample',ylab='rma normalized expression', col='blue')
dev.off()

#make design matrix from information in targets
TS <- paste(targets[,2],targets[,3],sep="")
TS <- factor(TS)
design <- model.matrix(~TS-1)
colnames(design)<-levels(TS)

#lmFit to the data. 
fit = lmFit(x,design)

#Contrast matrix for each comparison
cont.matrix = makeContrasts(ESTROGEN10=present10-absent10, 
                            ESTROGEN48=present48-absent48,
                            TIME.CONTROL=absent48-absent10,
                            TIME.ESTROGEN=present48-present10,
                            INTERACTION=(present48-present10)-(absent48-absent10),
                            levels=design)
fit1 = contrasts.fit(fit, cont.matrix)
fit1 = eBayes(fit1)

#get top 20 genes for each contrast. Default sorting is by B-statistic.
top1 = topTable(fit1, coef=1,number=20,lfc=1.5,sort.by="P")
top2 = topTable(fit1, coef=2,number=20,lfc=1.5,sort.by="P")
top3 = topTable(fit1, coef=3,number=20,lfc=1.5,sort.by="P")
top4 = topTable(fit1, coef=4,number=20,lfc=1.5,sort.by="P")
top5 = topTable(fit1, coef=5,number=20,lfc=1.5,sort.by="P")
#top2 in top1?
ov1 = rownames(top2)[rownames(top2) %in% rownames(top1)]
#top4 in top3?
ov2 = rownames(top4)[rownames(top4) %in% rownames(top3)]

#decide test using default parameters
results = decideTests(fit1)
#vennDiagram(results)

#how many control genes?
controls <- rownames(fit1$coefficients[grep("AFFX", rownames(fit1$coefficients)),])
control_expression <- x[rownames(x)%in%controls,]
control_variance <- apply(control_expression, 1, var)
#how many with variance below 1?
dim(control_expression[control_variance<1,])

#how many below given FDR?
for (val in c(0.001, 0.01, 0.05, 0.1, 0.2)){
  print(val)
  print(paste("coef 1: ", dim(topTable(fit1, coef=1,number=Inf,p.value=val))[1]))
  print(paste("coef 2: ", dim(topTable(fit1, coef=2,number=Inf,p.value=val))[1]))
  print(paste("coef 3: ", dim(topTable(fit1, coef=3,number=Inf,p.value=val))[1]))
  print(paste("coef 4: ", dim(topTable(fit1, coef=4,number=Inf,p.value=val))[1]))
  print(paste("coef 5: ", dim(topTable(fit1, coef=5,number=Inf,p.value=val))[1]))
}

#how many genes in FDR<0.01 are control?
top1a = topTable(fit1, coef=1,number=Inf,p.value=0.01, lfc=1.5)
top2a = topTable(fit1, coef=2,number=Inf,p.value=0.01, lfc=1.5)
top3a = topTable(fit1, coef=3,number=Inf,p.value=0.01, lfc=1.5)
top4a = topTable(fit1, coef=4,number=Inf,p.value=0.01, lfc=1.5)
top5a = topTable(fit1, coef=5,number=Inf,p.value=0.01, lfc=1.5)
dim(top1a[grep("AFFX",rownames(top1a)),])
dim(top2a[grep("AFFX",rownames(top2a)),])
dim(top3a[grep("AFFX",rownames(top3a)),])
dim(top4a[grep("AFFX",rownames(top4a)),])

#make volcano plots for each contrast 
top1b = topTable(fit1, coef=1,number=Inf)
top2b = topTable(fit1, coef=2,number=Inf)
top3b = topTable(fit1, coef=3,number=Inf)
top4b = topTable(fit1, coef=4,number=Inf)
top5b = topTable(fit1, coef=5,number=Inf)
#threshold for first 4 contrasts
top1b$threshold = as.factor(abs(top1b$logFC) > 1.5 & top1b$adj.P.Val < 0.01)
top2b$threshold = as.factor(abs(top2b$logFC) > 1.5 & top2b$adj.P.Val < 0.01)
top3b$threshold = as.factor(abs(top3b$logFC) > 1.5 & top3b$adj.P.Val < 0.01)
top4b$threshold = as.factor(abs(top4b$logFC) > 1.5 & top4b$adj.P.Val < 0.01)

#treatment 10h
png("vol1.png")
g = ggplot(data=top1b, aes(x=logFC, y=-log10(adj.P.Val), colour=threshold)) +
  geom_point(alpha=0.4, size=2.25) +
  xlim(c(-5, 5)) + ylim(c(0,5)) +
  opts(legend.position = "none") +
  xlab("log2 fold change") + ylab("-log10 FDR") +ggtitle("Treatment 10h")
g
dev.off()
#treatment 48h
png("vol2.png")
g = ggplot(data=top2b, aes(x=logFC, y=-log10(adj.P.Val), colour=threshold)) +
  geom_point(alpha=0.4, size=2.25) +
  xlim(c(-5, 5)) + ylim(c(0,5)) +
  opts(legend.position = "none") +
  xlab("log2 fold change") + ylab("-log10 FDR") +ggtitle("Treatment 48h")
g
dev.off()
#time effect control
png("vol3.png")
g = ggplot(data=top3b, aes(x=logFC, y=-log10(adj.P.Val), colour=threshold)) +
  geom_point(alpha=0.4, size=2.25) +
  xlim(c(-5, 5)) + ylim(c(0,5)) +
  opts(legend.position = "none") +
  xlab("log2 fold change") + ylab("-log10 FDR") +ggtitle("Time effect control")
g
dev.off()
#time effect estrogen
png("vol4.png")
g = ggplot(data=top4b, aes(x=logFC, y=-log10(adj.P.Val), colour=threshold)) +
  geom_point(alpha=0.4, size=2.25) +
  xlim(c(-5, 5)) + ylim(c(0,5)) +
  opts(legend.position = "none") +
  xlab("log2 fold change") + ylab("-log10 FDR") +ggtitle("Time effect estrogen")
g
dev.off()


