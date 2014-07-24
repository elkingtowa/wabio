#REPEATING GDS2938 analysis with pre-normalized data
library(GEOquery)
library(ggplot2)
library(limma)
library(affy)

gds2938 <- getGEO(filename='GDS2938.soft.gz')
eset <- GDS2eSet(gds2938)
x <- exprs(eset)

#redo boxplot and check if equal
png('eset_boxplot.png', height=350, width=525)
boxplot(x, names=1:12,
        main='Normalized expression across samples',
        xlab='Sample',ylab='pre-normalized expression', col='blue')
dev.off()

#make design matrix from information in targets - read phenodata from previous sample
datadir = "C:/GSE5054"
targets<-readTargets("phenodata.txt",path=datadir,sep="",row.names="filename")
TS <- paste(targets[,2],targets[,3],sep="")
TS <- factor(TS)
design <- model.matrix(~TS-1)
colnames(design)<-c('control','IFN.gamma','IL1.beta','both')

#lmFit to the data. 
fit = lmFit(x,design)
cont.matrix = makeContrasts(IFN=IFN.gamma-control,
                            IL1=IL1.beta-control,
                            BOTH=both-control,
                            levels=design)
fit1 = contrasts.fit(fit, cont.matrix)
fit1 = eBayes(fit1)

#get top 20 genes for each contrast. Default sorting is by B-statistic.
top1 = topTable(fit1, coef=1,number=20)
top2 = topTable(fit1, coef=2,number=20)
top3 = topTable(fit1, coef=3,number=20)
top1a = topTable(fit1, coef=1,number=20,lfc=1.5)
top2a = topTable(fit1, coef=2,number=20,lfc=1.5)
top3a = topTable(fit1, coef=3,number=Inf,lfc=1.5)              

#make volcano plots for each contrast 
top1b = topTable(fit1, coef=1,number=Inf)
top2b = topTable(fit1, coef=2,number=Inf)
top3b = topTable(fit1, coef=3,number=Inf)
#threshold for contrasts
top1b$threshold = as.factor(top1b$adj.P.Val < 0.01)
top2b$threshold = as.factor(top2b$adj.P.Val < 0.01)
top3b$threshold = as.factor(top3b$adj.P.Val < 0.01)

#IFN-gamma treatment
png("vol1.png")
g = ggplot(data=top1b, aes(x=logFC, y=-log10(adj.P.Val), colour=threshold)) +
  geom_point(alpha=0.4, size=2.25) +
  xlim(c(-5, 5)) + ylim(c(0,5)) +
  opts(legend.position = "none") +
  xlab("log2 fold change") + ylab("-log10 FDR") +ggtitle("IFN-gamma treatment")
g
dev.off()
#IL1-beta treatment
png("vol2.png")
g = ggplot(data=top2b, aes(x=logFC, y=-log10(adj.P.Val), colour=threshold)) +
  geom_point(alpha=0.4, size=2.25) +
  xlim(c(-5, 5)) + ylim(c(0,5)) +
  opts(legend.position = "none") +
  xlab("log2 fold change") + ylab("-log10 FDR") +ggtitle("IL1-beta treatment")
g
dev.off()
#both treatments 
png("vol3.png")
g = ggplot(data=top3b, aes(x=logFC, y=-log10(adj.P.Val), colour=threshold)) +
  geom_point(alpha=0.4, size=2.25) +
  xlim(c(-5, 5)) + ylim(c(0,5)) +
  opts(legend.position = "none") +
  xlab("log2 fold change") + ylab("-log10 FDR") +ggtitle("both treatments")
g
dev.off()
