hrv1A = strsplit(paste(read.table('hrv1A.fasta',comment.char=">",stringsAsFactors=F)$V1, collapse=''),split='')[[1]]
hrv16 = strsplit(paste(read.table('hrv16.fasta',comment.char=">",stringsAsFactors=F)$V1, collapse=''),split='')[[1]]

mat1 = matrix(hrv1A, nrow=length(hrv1A), ncol=length(hrv1A))
mat2 = matrix(hrv16, nrow=length(hrv16), ncol=length(hrv16), byrow=T)

dot_mat = 1*(mat1==mat2)
image(dot_mat, col=c("white","black"))

dot <- function(seq1, seq2, split){
  if (split) {
    seq1 <- strsplit(seq1,split="")[[1]]
    seq2 <- strsplit(seq2,split="")[[1]]
  }

  mat1 = matrix(seq1, nrow=length(seq1), ncol=length(seq1))
  mat2 = matrix(seq2, nrow=length(seq2), ncol=length(seq2), byrow=T)
  dot_mat = 1*matrix(mapply(identical, mat1,mat2),ncol=length(seq1), nrow=length(seq2))
  image(dot_mat, col=c("white","black"))
}

# try on big seq
h1mb = strsplit(paste(read.table('human_1mb.fasta',comment.char=">",stringsAsFactors=F)$V1, collapse=''),split='')[[1]]
m1mb = strsplit(paste(read.table('mouse_1mb.fasta',comment.char=">",stringsAsFactors=F)$V1, collapse=''),split='')[[1]]

