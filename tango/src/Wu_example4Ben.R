HIV="ACTCTTTGGCAACGACCCCTCGTCACAATAAAGATAGGGGGGCAACTAAAGGAAGCTCTATTAGATACAGGAGCAGATGATACAGTATTAGAAGAAATGAGTTTGCCAGGAAGATGGAAACCAAAAATGATAGGGGGAATTGGAGGTTTTATCAAAGTAAGACAGTATGATCAGATACTCATAGAAATCTGTGGACATAAAGCTATAGGTACAGTATTAGTAGGACCTACACCTGTCAACATAATTGGAAGAAATCTGTTGACTCAGATTGGTTGCACTTTAAATTTTATGGAAAAGGAAGGGAAAATTTCAAAAATTGGGCCTGAAAATCCATACAATACTCCAGTATTTGCCATAAAGAAAAAAGACAGTACTAAATGGAGAAAATTAGTAGATTTCAGAGAACTTAATAAGAGAACTAAGACTTCTGGGAAGTTCAATTAGGAATACCACATCCCGCAGGGTTAAAAAAGAAAAAATCAGTAACAGTACTGGATGTGGGTGATGCATATTTTTCAGTTCCCTTAGATGAAGACTTCAGGAAGTATACTGCATTTACCATACCTAGTATAAACAATGAGACACCAGGGATTAGATATCAGTACAATGTGCTTCCACAGGGATGGAAAGGATCACCAGCAATATTCCAAAGTAGCATGACAAAAATCTTAGAGCCTTTTAGAAAACAAAATCCAGACATAGTTATCTATCAATACATGGATGATTTGTATGTAGGATCTGACTTAGAAATAGGGCAGCATAGAACAAAAATAGAGGAGCTGAGACAACATCTGTTGAGGTGGGGACTTACCACACCAGACAAAAAACATCAGAAAGAACCTCCATTCCTTTGGATGGGTTATGAACTCCATCCTGATAAATGGACA"
#### example data from http://www.ncbi.nlm.nih.gov/nuccore/574960918?report=fasta
FLU="ATGGACTCCAACACCATGTCAAGCTTTCAGGTAGACTGTTTCCTTTGGCATATCCGCAAGCGATTTGCAGACAATGGATTGGGTGATGCCCCATTCCTTGATCGGCTCCGCCGAGATCAAAAGTCCTTAAAAGGAAGAGGCAACACCCTTGGCCTCGATATCGAAACAGCCACTCTTGTTGGGAAACAAATCGTGGAATGGATCTTGAAAGAGGAATCCAGCGAGACACTTAGAATGACAATTGCATCTGTACCTACTTCGCGCTACCTTTCTGACATGACCCTCGAGGAAATGTCACGAGACTGGTTCATGCTCATGCCTAGGCAAAAGATAATAGGCCCTCTTTGCGTGCGATTGGACCAGGCGATCATGGAAAAGAACATAGTACTGAAAGCGAACTTCAGTGTAATCTTTAACCGATTAGAGACCTTGATACTACTAAGGGCTTTCACTGAGGAGGGAGCAATAGTTGGAGAAATTTCACCATTACCTTCTCTTCCAGGACATACTTATGAGGATGTCAAAAATGCAGTTGGGGTCCTCATCGGAGGACTTGAATGGAATGGTAACACGGTTCGAGTCTCTGAAAATATACAGAGATTCGCTTGGAGAAACTGTGATGAGAATGGGAGACCTTCACTACCTCCAGAGCAGAAATGAAAAGTGGCGAGAGCAATTGGGACAGAAATTTGAGGAAATAAGGTGGTTAATTGAAGAAATGCGGCACAGATTGAAAGCGACAGAGAATAGTTTCGAACAAATAACATTTATGCAAGCCTTACAACTACTGCTTGAAGTAGAACAAGAGATAAGAGCTTTCTCGTTTCAGCTTATTTAA"
#http://www.ncbi.nlm.nih.gov/nuccore/227809838?report=fasta

n1=nchar(HIV);n2=nchar(FLU)
x3.hiv=substring(HIV,1:(n1-2),3:n1)
x3.flu=substring(FLU,1:(n2-2),3:n2)

T3.hiv=table(x3.hiv)### empirical counts of all three letter word
T3.hiv### you will see ther are missing data so some words are "nAA" which is not a real word.
T3.flu=table(x3.flu)

### creating the three letter word
a1=c("A","T","G","C")
a2=expand.grid(a1,a1)
a2=apply(as.matrix(a2),1,paste,collapse="")
a3=expand.grid(a2,a1)
a3=apply(as.matrix(a3),1,paste,collapse="") ## all 4^3=64 of these
### you also can create the 4-letter words this way. Same of doing loops but easier to write.

##let's keep only the ones of the real 3-letter words
## also make sure if a word in the dictionary does not appear in a table, 
## it returns 0 instead of NA
clean.table=function(x,dict){
  new.table=rep(0,length(dict))
  new.table=x[dict]
  new.table[is.na(new.table)]=0
  names(new.table)=dict
  new.table
}

T3.hiv2=clean.table(T3.hiv,dict=a3)
T3.flu2=clean.table(T3.flu,dict=a3)

##now make the counts into relative frequency
T3.hiv2=T3.hiv2/sum(T3.hiv2)
T3.flu2=T3.flu2/sum(T3.flu2)

plot(T3.hiv2,T3.flu2)
abline(0,1)
 
k=identify(T3.hiv2,T3.flu2)## i highlited a few, you could highlight others

#> k
#[1]  1 12 18 24 38
text(T3.hiv2[k],T3.flu2[k],a3[k])

## so if the two viruses are generated from a common word frequency, 
## i would estimate it jointly
T0=(T3.hiv2+T3.flu2)/2

### suppose we run into a virus sequence,
## since flu is more more common than HIV, I give a prior probability of 0.9 to be
## flu, after we get a piece of its sequence, we might compute a posterior prob

virusX=substr(HIV,100,150)
virusX=table(substring(virusX,1:(nchar(virusX)-2),3:nchar(virusX)))
virusX=clean.table(virusX,dict=a3)

pX.given.hiv=dmultinom(virusX,size=sum(virusX),prob=T3.hiv2)
pX.given.flu=dmultinom(virusX,size=sum(virusX),prob=T3.flu2)
phiv.given.X=(0.10*pX.given.hiv)/(0.10*pX.given.hiv+0.90*pX.given.flu)
phiv.given.X

## this example shows that even if I only sequenced from 100 to 150 , 
##that is 50nt of the HIV sequence, I'd be able to have close to 100% posterior probability that it is HIV.
## but of course this is this particular example. It does not mean that based on that 50nt piece one
## has confidence it is HIV. IT says, if you already know it is either HIV or FLU,
## then you are pretty sure it is HIV instead of FLU.
 

