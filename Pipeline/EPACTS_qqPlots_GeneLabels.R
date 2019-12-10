#create qq plots
#R
library(Cairo)
#library(ggplot2)
#library(ggrepel)

options(echo=FALSE) # don't see commands in output file
args <- commandArgs(trailingOnly = TRUE)
argslen<-length(args);

path<-'./';

if (argslen>0)
{
	path<-args[1]
}

print(path);

minVal<-2.76
if (argslen>1)
{
	minVal<-args[2];
}

print(minVal);

files <- system(paste('ls ', path, '/*.epacts | grep -v case.ctl',sep=""), intern=T);

for(i in 1:length(files)){
	pvals <- read.table(files[i],header=T)[,10]
	genes <- read.table(files[i],header=T)[,4]

	pg = data.frame(pval=pvals,gene=gsub("(^.*_)","", genes) )
	pg<-pg[order(pg$pval),];

	#pvals <- read.table(files[i],header=T)[,9]
	#pvals <- read.table(files[i],header=T)[,11]

	N <- length(pg$pval)
	ci <- 0.975
	df = data.frame(observed=-log10(pg$pval),
                expected=-log10(1:N / N),
                cupper=-log10(qbeta(ci, 1:N, N - (1:N) + 1)),
                clower=-log10(qbeta(1 - ci, 1:N, N - (1:N) + 1)),
		genename=pg$gene)



	log10Pe = expression(paste("Expected -log"[10], plain(P)))
	log10Po = expression(paste("Observed -log"[10], plain(P)))

	Cairo(file=paste(files[i],"LABELS.png",sep="."), width=500, height=500, bg="white", canvas="white", type="png")

	#LASTEDIT:
	#plot(df$expected,df$observed,xlab=log10Pe,ylab=log10Po, type='n',cex.axis=1.2,cex.lab=1.2,ylim=c(0,9.5))
	plot(df$expected,df$observed,xlab=log10Pe,ylab=log10Po, type='n',cex.axis=1.2,cex.lab=1.2)
	#plot(0,0);

	#p<-qplot(expected,observed, data=df, xlab=log10Pe,ylab=log10Po, ylim=c(0,11))
	#plot(p);

	topgenes<-(df$observed>minVal)
	polygon(c(df$expected,rev(df$expected)),c(df$cupper,rev(df$clower) ),col='lightgrey',border=NA)
	abline(0,1)
	points(df$expected,df$observed,pch=20,cex=1)

	text(df$expected[topgenes],df$observed[topgenes], labels=df$genename[topgenes], cex=0.8, pos=2); #sample(1:4,length(df$genename),replace=T) )


	dev.off()

}
