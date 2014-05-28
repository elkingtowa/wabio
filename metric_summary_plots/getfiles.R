getfiles<-function(ancestry,ver=37) {
##this file reads in a flat file named 'summary_list' which lists the absolute path for all relevant summary files
##INPUTS: 
##	ancestry
##	ver=36 or 37(default), 

	if (ver == 36 || ver ==37) {
	print(ver)
		#################read in the filelist for summary file locations#################
		files<-as.character(unlist(read.table('summary_list')))

		#################reorder fileorder based on ancestry annotation#################
		file.reorder.id<-match(ancestry[,1],sub('GS','NA',sub(paste('-1100-',ver,'-ASM.tsv',sep=''),'',sub('summary-','',basename(files)))))
		tmp<- files[file.reorder.id]
		files<-tmp

		return(files)
	} else {

	stop('Incorrect genome version')
	}	


}