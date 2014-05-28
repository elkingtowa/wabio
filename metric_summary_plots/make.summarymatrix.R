make.summarymatrix<-function(files,ancestry) {
##This function generates a matrixed datasset of the summary files
##INPUT: files: file listing the absolute path for the summary files

	summary.mat<-NULL
	for (f in 1:length(files)) {
		cat('processing file #: ',f,'...\n')
		f1<-read.table(files[f],fill=T,sep='\t',header=T)
		if (f==1) 
			summary.mat<-f1
		if (f>1)
			summary.mat<-cbind(summary.mat,f1[,2])
		if (f==length(files)) {
			rownames(summary.mat)<-f1[,1]
			summary.mat<-summary.mat[,-1]
			colnames(summary.mat)<-ancestry[,1]	
		}
	}
	return(summary.mat)
}