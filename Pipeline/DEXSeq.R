library("DEXSeq")
library("DESeq2")
library("Cairo")
library("RColorBrewer")
library("pvclust")
library("ggplot2")
library("gplots")
library("gdata")
library("genefilter")
library("matrixStats")

condition.map <- function(cond.biol) { if (cond.biol == 1) "case" else "control" }

myargs      <- commandArgs(trailingOnly = TRUE)
caseFile    <- myargs[1]
controlFile <- myargs[2]
outputFolder<- myargs[3]
projectName <- myargs[4]

caseTable    <- read.table(file=caseFile, sep="\t")
controlTable <- read.table(file=controlFile, sep="\t")

countFileVector <- c()
nameVector  <- c()
conditionVector <- c()

for (i in 1:nrow(caseTable)) {
	countFileVector <- c(countFileVector, paste(caseTable[i,2], ".dexseq.counts", sep=""))
	nameVector      <- c(nameVector, paste(caseTable[i,1]))
	conditionVector <- c(conditionVector, "case") 
}
for (i in 1:nrow(controlTable)) {
	countFileVector <- c(countFileVector, paste(controlTable[i,2], ".dexseq.counts", sep=""))
	nameVector      <- c(nameVector, paste(controlTable[i,1]))
	conditionVector <- c(conditionVector, "control") 
}

flattenedFile = "/data/isilon/seq/analysis/exomehg19/IDG_Floss_RNAseq/RNASeq/DEU/Mus_musculus.NCBIM37.67.gff.chrtransformed"
save.image()
sampleTable = data.frame(
		row.names=nameVector,
		condition=conditionVector
)
dxd = DEXSeqDataSetFromHTSeq(
		countFileVector,
		sampleData=sampleTable,
		design= ~ sample + exon + condition:exon,
		flattenedfile=flattenedFile
)

dxd = estimateSizeFactors( dxd )
dxd = estimateDispersions( dxd )
dxd = testForDEU( dxd )
dxd = estimateExonFoldChanges( dxd, fitExpToVar="condition" )
res =  DEXSeqResults( dxd )

Cairo(file=paste(outputFolder, "REST.exons.png", sep="/"), widht=800, height=480, bg="white", canvas="white", type="png")
#ENSMUSG00000029249 = REST gene
plotDEXSeq( res, "ENSMUSG00000029249", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2)
dev.off()

#save this session data
save.image(paste(outputFolder, "session.RData", sep="/"))
