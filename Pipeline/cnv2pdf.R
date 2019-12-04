library('ExomeDepth')

myargs <- commandArgs(trailingOnly = TRUE)
sample        <-myargs[1]
ExomeCountFile<-myargs[2]
outfile       <-myargs[3]

#outfile <- 'cnv.pdf'
#ExomeCountFile <- 'cnv.csv'

ExomeCount<-read.csv(ExomeCountFile,header=TRUE)
#sample<-colnames(ExomeCount[6])
#sample<-sub('Exome','',sample)
ExomeCount[,'space']<-paste('chr', ExomeCount[,'space'], sep='')

#################################
# print graphic
#################################
chr<-c('chr1','chr2','chr3','chr4','chr5',
'chr6','chr7','chr8','chr9','chr10','chr11',
'chr12','chr13','chr14','chr15','chr16',
'chr17','chr18','chr19','chr20','chr21','chr22')


pdf(file=outfile,paper='A4')
par(oma=c(1,1,2,1))
par(mar=c(4,4,3,2))
par(mfrow=c(2,1))

for(chrom in chr) {
	sl <-ExomeCount[,'space'] == chrom
	tmp<-ExomeCount[sl,]
	plot(
	tmp[,'start'],
	tmp[,6],
	pch=20,cex=0.8,
	xlim=c(0,tmp[length(tmp[,1]),'end']),
	ylim=c(0,5),
	xlab=chrom,
	ylab='Ratio',
	axes=FALSE,
	main=chrom
	)
	xlabels<-seq(0,3e2,1e1)
	axis(side=1,at=seq(0,3e8,1e7),labels=xlabels,pos=c(-0.1,0))
	ylabels<-seq(0,150,1)
	axis(side=2,pos=0,at=seq(0,150,1),labels=ylabels)
	axis(side=3,labels=F,lwd.ticks=-1,at=c(0,3e8))
	axis(side=4,labels=F,lwd.ticks=-1,at=c(-0.1,105))
	abline(h=1,col='red')
        if (chrom=='chr1') {
                tmp.header<-paste(date(),'Sample:',sample)
                mtext(tmp.header,outer=TRUE,adj=0,padj=-1,cex=0.8) 
        }
}

dev.off()

