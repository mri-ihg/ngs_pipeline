# build the exomeCountObject for ExomeDepth 
#source('getBamCounts.R')

# output TestCount

library('ExomeDepth')
myargs      <-commandArgs(trailingOnly = TRUE)
bamFile     <-myargs[1]
target.file <-myargs[2]
refFile     <-myargs[3]
outFile     <-myargs[4]
#cat (bamFile)
#cat (target.file)

#if (assay == 'SureSelect50Mbv5') {
#	target.file<-'/usr/local/packages/seq/tools/solexa/exome50MbV5/ontarget.merged.coding_exons.bed'
#}
#if (assay == 'SureSelect50Mbv4') {
#	target.file<-'/usr/local/packages/seq/tools/solexa/exome50MbV4/ontarget.merged.coding_exons.bed'
#}

#bam.files muss vector sein
#bam.files <- read.table(file=bamFile, header=F)
#bam.files <- as.vector(as.matrix(bam.files))

bam.files <- strsplit(bamFile,",")[[1]]


#bed.frame muss data.frame sein
tmp       <- read.table(file= target.file, header=F)
#tmp       <- as.data.frame(tmp)
#tmp       <- tmp[1:10,,drop=F]
#bed.frame <- t(apply(tmp,1,FUN=function(y) {(unlist(strsplit(as.character(y),"[:-]")))} ))
bed.frame <- as.data.frame(tmp)
bed.frame[,2] <- as.numeric(as.character(bed.frame[,2]))
bed.frame[,3] <- as.numeric(as.character(bed.frame[,3]))

TestCount<-getBamCounts(
bed.frame = bed.frame,
bam.files = bam.files,
min.mapq = 40,
read.width = 300,
include.chr = FALSE,
referenceFasta = refFile)

save(TestCount,file=outFile)

