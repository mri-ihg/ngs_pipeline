library(gage)
library(GenomicFeatures)
library(Rsamtools)
library(annotate)
library(pathview)

require(gage)
require(gageData)

myargs       <- commandArgs(trailingOnly = TRUE)
settings    <- myargs[1]	#mm9, mm10 or hg19
resFile <- myargs[2]		#file where the DE-results are stored (with logFC)
outputFolder <- myargs[3]
experiment <- myargs[4]

pathviewSpecies = ""
keggSet = ""


#read result table in res.gage
res.gage <- read.table(resFile, header=T, row.names=1, check.names=F, quote="\"", sep=",")

if (startsWith(settings, "hg19")) {
	data(kegg.sets.hs)
	library(org.Hs.eg.db)
	
	dict <- id2eg(rownames(res.gage), category="SYMBOL", org="Hs")
	
	pathviewSpecies = "hsa"
	keggSet <- kegg.sets.hs
} else if (settings == "mm10" || settings == "mm9" || settings == "mouse")  {
	data(kegg.sets.mm)
	library(org.Mm.eg.db)
	
	dict <- id2eg(rownames(res.gage), category="SYMBOL", org="Mm")
	
	pathviewSpecies <- "mmu"
	keggSet <- kegg.sets.mm
}

res.gage.df <- data.frame(res.gage, entrezid=dict[which(dict[,1] == row.names(res.gage)), 2])
res.gage.df <- res.gage.df[!is.na(res.gage.df$entrezid),]

fc.gage <- res.gage.df$log2FoldChange
names(fc.gage) <- res.gage.df$entrezid

fc.gage[fc.gage>10] = 10
fc.gage[fc.gage<-10] = -10

out.suffix = "DESeq2"
fc.kegg.p <- gage(fc.gage, gsets=keggSet, ref=NULL, samp=NULL)

setwd(outputFolder)
sel <- fc.kegg.p$greater[, "q.val"] < 0.1 & !is.na(fc.kegg.p$greater[,"q.val"])
path.ids <- rownames(fc.kegg.p$greater)[sel]
path.ids2 <- substr(path.ids, 1, 8)
pv.out.list <- sapply(path.ids2, function(pid) pathview(gene.data = fc.gage, pathway.id = pid,species = pathviewSpecies))

sel <- fc.kegg.p$less[, "q.val"] < 0.1 & !is.na(fc.kegg.p$less[,"q.val"])
path.ids <- rownames(fc.kegg.p$less)[sel]
path.ids2 <- substr(path.ids, 1, 8)
pv.out.list <- sapply(path.ids2, function(pid) pathview(gene.data = fc.gage, pathway.id = pid,species = pathviewSpecies))

fc.kegg.greater <- na.omit(fc.kegg.p$greater)
fc.kegg.less    <- na.omit(fc.kegg.p$less)

fc.kegg.greater.df <- data.frame(pathway=rownames(fc.kegg.greater), fc.kegg.greater[,3:4])
fc.kegg.less.df    <- data.frame(pathway=rownames(fc.kegg.less), fc.kegg.less[,3:4])

write.table(fc.kegg.greater.df, file=paste(outputFolder, "/", experiment, ".", "upregulated.pathways.tsv",   sep=""), sep="\t", row.names=F)
write.table(fc.kegg.less.df,    file=paste(outputFolder, "/", experiment, ".", "downregulated.pathways.tsv", sep=""), sep="\t", row.names=F)

save.image(paste(outputFolder, "session.RData", sep="/"))
