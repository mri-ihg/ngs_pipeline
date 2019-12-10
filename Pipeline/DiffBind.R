library(DiffBind)
library(Cairo)

myargs       <- commandArgs(trailingOnly = TRUE)
outputFolder <- myargs[1]
sampleSheet  <- myargs[2]
species      <- myargs[3]

db <- dba(sampleSheet=sampleSheet)
db <- dba.count(db, summits=250)
db <- dba.contrast(db, categories=DBA_CONDITION, minMembers=2)
db <- dba.analyze(db)
db.rep <- dba.report(db, th=1)
db.df <- as.data.frame(db.rep)

Cairo(file=paste(outputFolder,"PCA_ob.png",sep="/"), widht=640, height=640, bg="white", canvas="white", type="png")
dba.plotPCA(db, label=DBA_ID)
dev.off()
Cairo(file=paste(outputFolder, "Heatmap_ob.png", sep="/"), widht=640, height=640, bg="white", canvas="white", type="png")
plot(db)
dev.off()

write.table(db.df[order(db.df$FDR, db.df$p.value),], file=paste(outputFolder,"DEPeaks.tsv",sep="/"), row.names=F, sep="\t",quote = FALSE)


## peak annotation ##
# from: https://www.bioconductor.org/packages/devel/bioc/vignettes/ChIPseeker/inst/doc/ChIPseeker.html
library(ChIPseeker)
library(clusterProfiler)
library(ChIPpeakAnno)

if (species == "human") {
	library(TxDb.Hsapiens.UCSC.hg19.knownGene)
	library(EnsDb.Hsapiens.v75) ##(hg19)
	library(org.Hs.eg.db)
	txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
	ensdb <- EnsDb.Hsapiens.v75
	eg.db <- "org.Hs.eg.db"
} else if (species == "mouse") {
	library(TxDb.Mmusculus.UCSC.mm10.knownGene)
	library(EnsDb.Mmusculus.v79)
	library(org.Mm.eg.db)
	txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
	ensdb <- EnsDb.Mmusculus.v79
	eg.db <- "org.Mm.eg.db"
} else {
	txdb <- ""
	ensdb <- ""
	eg.db <- ""
	print("Unknown species specified. Annotation not possible.")
}


peak <- db.rep
Cairo(file=paste(outputFolder, "Peak_overview.png", sep="/"), width=1280, height=720, bg="white", canvas="white", type="png")
covplot(peak, weightCol="Fold")
dev.off()

peakAnno <- annotatePeak(db.rep, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb=eg.db)
Cairo(file=paste(outputFolder, "Peak_annotation_pie.png", sep="/"), width=640, height=640, bg="white", canvas="white", type="png")
plotAnnoPie(peakAnno)
dev.off()

Cairo(file=paste(outputFolder, "Peak_annotation_bar.png", sep="/"), width=640, height=640, bg="white", canvas="white", type="png")
plotAnnoBar(peakAnno)
dev.off()

Cairo(file=paste(outputFolder, "Peak_upsetplot.png", sep="/"), width=720, height=480, bg="white", canvas="white", type="png")
upsetplot(peakAnno)
dev.off()

Cairo(file=paste(outputFolder, "Peak_upsetplot_pie.png", sep="/"), width=720, height=480, bg="white", canvas="white", type="png")
upsetplot(peakAnno, vennpie=TRUE)
dev.off()

Cairo(file=paste(outputFolder, "Peak_distance2TSS.png", sep="/"), width=720, height=480, bg="white", canvas="white", type="png")
plotDistToTSS(peakAnno, title="Distribution of transcription factor-binding loci\nrelative to TSS")
dev.off()

#annoData <- toGRanges(TxDb.Hsapiens.UCSC.hg19.knownGene, feature="CDS")
annoData <- toGRanges(ensdb, feature="gene")

aCR<-assignChromosomeRegion(db.rep, nucleotideLevel=FALSE, precedence=c("Promoters", "immediateDownstream", "fiveUTRs", "threeUTRs", "Exons", "Introns"),	TxDb=txdb)

Cairo(file=paste(outputFolder, "Peak_featureBarplot.png", sep="/"), width=720, height=480, bg="white", canvas="white", type="png")
barplot(aCR$percentage, las=3)
dev.off()

overlaps.anno <- annotatePeakInBatch(db.rep, AnnotationData=annoData,	output="nearestBiDirectionalPromoters",	bindingRegion=c(-2000, 2000))
library(org.Hs.eg.db)
overlaps.anno <- addGeneIDs(overlaps.anno, eg.db, IDs2Add = c("genename", "refseq", "symbol", "entrez_id"))
#head(overlaps.anno)

write.csv(as.data.frame(unname(overlaps.anno)), file=paste(outputFolder, "anno.csv", sep="/"))

Cairo(file=paste(outputFolder, "Peak_featurePie.png", sep="/"), width=720, height=480, bg="white", canvas="white", type="png")
pie1(table(overlaps.anno$insideFeature))
dev.off()
