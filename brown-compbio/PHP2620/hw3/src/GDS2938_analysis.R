#analysis of GSE5054 data
library(affy)
library(ggplot2)
library(limma)

#load raw CEL files
datadir = "C:/GSE5054"
targets<-readTargets("phenodata.txt",path=datadir,sep="",row.names="filename")
ab<-ReadAffy(filenames=targets$filename,celfile.path=datadir,phenoData=targets)

#image of each array
for (i in 1:length(ab@phenoData@data$filename)){
  png(paste(ab@phenoData@data$filename[i],".png", sep=''))
  image(ab[,i])
  dev.off()
}

#boxplot of raw intensities
png('raw_boxplot.png', height=350, width=525)
boxplot(ab, names=1:12,
        main='Raw probe intensity across samples',
        xlab='Sample',ylab='log2(intensity)', col='red')
dev.off()

#do RMA and convert to expression set
eset<-rma(ab)
x <- exprs(eset)

#redo boxplot and check if equal
png('eset_boxplot.png', height=350, width=525)
boxplot(x, names=1:12,
        main='Normalized expression across samples',
        xlab='Sample',ylab='rma normalized expression', col='blue')
dev.off()

#make design matrix from information in targets
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
top3a = topTable(fit1, coef=3,number=20,lfc=1.5)
top1b = topTable(fit1, coef=1,number=Inf)
top2b = topTable(fit1, coef=2,number=Inf)
top3b = topTable(fit1, coef=3,number=Inf)                            

#volcano plots
#IFN-gamma
png("vol1.png")
g = ggplot(data=top1b, aes(x=logFC, y=-log10(P.Value))) +
  geom_point(alpha=0.4, size=2.25) +
  xlim(c(-5, 5)) + ylim(c(0,5)) +
  opts(legend.position = "none") +
  xlab("log2 fold change") + ylab("-log10 p-value") +ggtitle("IFN-gamma treatment")
g
dev.off()

#IL1-beta   
png("vol2.png")
g = ggplot(data=top2b, aes(x=logFC, y=-log10(P.Value))) +
  geom_point(alpha=0.4, size=2.25) +
  xlim(c(-5, 5)) + ylim(c(0,5)) +
  opts(legend.position = "none") +
  xlab("log2 fold change") + ylab("-log10 p-value") +ggtitle("IL1-beta treatment")
g
dev.off()

#both   
png("vol3.png")
g = ggplot(data=top3b, aes(x=logFC, y=-log10(P.Value))) +
  geom_point(alpha=0.4, size=2.25) +
  xlim(c(-5, 5)) + ylim(c(0,5)) +
  opts(legend.position = "none") +
  xlab("log2 fold change") + ylab("-log10 p-value") +ggtitle("Both treatments")
g
dev.off()