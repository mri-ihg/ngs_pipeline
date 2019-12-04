#run all neccessary steps for RNA-seq analysis
#main component of this script is the DESeq2 package which handles most of the RNA-seq analysis
######

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

myargs       <- commandArgs(trailingOnly = TRUE)
countFile    <- myargs[1]
conditionFile <- myargs[2]
outputFolder <- myargs[3]
projectName  <- myargs[4]
fpkmFile <- myargs[5]
noCreateHeatMap <- myargs[6]

countTable = read.table(countFile, header=T, row.names=1, check.names=F, quote="\"")
countTable = countTable[1:(nrow(countTable)-5),] 			#delete last 5 lines; includes unneeded information
conditionTable = read.table(conditionFile, header=F, check.names=F, row.names=1)
names(conditionTable)[1] <- paste("condition")
conditions = c()
for (i in 1:ncol(countTable)) {
	conditions <- c(conditions, conditionTable[which(rownames(conditionTable)==colnames(countTable)[i]), 1])
}
formatted.conditionTable <- transform(conditionTable, condition=as.factor(condition))

#DESeq2: create a new DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = countTable, colData=formatted.conditionTable, design= ~ condition)
dds$condition <- relevel(dds$condition, "0")
dds <- DESeq(dds, betaPrior=F)
res <- results(dds, cooksCutoff=FALSE, independentFiltering=FALSE)

rld <- rlog(dds)
vsd <- varianceStabilizingTransformation(dds)
rlogMat <- assay(rld)
vstMat <- assay(vsd)

select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:30]
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)

distsRL <- dist(t(assay(rld)))
mat <- as.matrix(distsRL)
rownames(mat) <- colnames(mat) <- with(colData(dds), {x<-row.names(colData(dds)); paste(x, colData(dds)$condition, sep=":")})
Cairo(file=paste(outputFolder, "heatmap.png", sep="/"), width=800, height=800, bg="white", canvas="white", type="png")
heatmap.2(mat, trace="none", col = rev(hmcol), margin=c(13, 13))
dev.off()

#taken and modified from DESeq2::plotPCA
intgroup = "condition"
ntop = 500
rv <- rowVars(assay(rld))
select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
pca <- prcomp(t(assay(rld)[select, ]))
percentVar <- pca$sdev^2/sum(pca$sdev^2)
if (!all(intgroup %in% names(colData(rld)))) {
	stop("the argument 'intgroup' should specify columns of colData(dds)")
}
intgroup.df <- as.data.frame(colData(rld)[, intgroup, drop = FALSE])
d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], group = sapply(unlist(intgroup.df), condition.map), intgroup.df, names = colnames(rld))
Cairo(file=paste(outputFolder, "PCA.png", sep="/"), width=640, height=480, bg="white", canvas="white", type="png")
ggplot(data = d, aes_string(x = "PC1", y = "PC2", color = "group")) + 
			geom_point(size = 3) + 
			geom_text(aes(label=rownames(d)), hjust=0, vjust=0, show_guide=F, size=4) +
			xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) +
			ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) +
			coord_cartesian(xlim=c(min((pca$x)[,"PC1"]) - (2/100*(max((pca$x)[,"PC1"]) - min((pca$x)[,"PC1"]))), max((pca$x)[,"PC1"])+ (20/100*(max((pca$x)[,"PC1"]) - min((pca$x)[,"PC1"])))))	+
			ggtitle(paste("Project: ", projectName, sep=""))
dev.off()

#plotMA
Cairo(file=paste(outputFolder, "MAPlot.png", sep="/"), widht=800, height=480, bg="white", canvas="white", type="png")
plotMA(res, main=paste("Project: ", projectName, sep=""))
dev.off()


Cairo(file=paste(outputFolder, "EstimatedDispersion.png", sep="/"), widht=800, height=480, bg="white", canvas="white", type="png")
plotDispEsts(dds, main=paste("Project: ", projectName, sep=""))
dev.off()


#plot heatmap of top 100 rowMeans Variance Transcripts
distCor <- function(x) as.dist(1-cor(t(x)))				#for heatmap.2: distfun=distCor
hclustAvg <- function(x) hclust(x, method="average")	#for heatmap.2: hclustfun=hclustAvg

select = order(rowMeans(assay(vsd)), decreasing=T)[1:100]
res.hm = assay(vsd)[select,]

column.conditions <- unlist(lapply(conditions, condition.map))
column.names <- paste(column.conditions, colnames(res.hm), sep=":")
colnames(res.hm) <- column.names
color.map <- function(mol.biol) { if (startsWith(mol.biol, "case", trim=T, ignore.case=T)) "grey50" else "grey80" }
columncolors <- unlist(lapply(colnames(res.hm), color.map))
Cairo(file=paste(outputFolder, "top1000_rowMeans_HeatMap.png", sep="/"), width=800, height=1200, bg="white", canvas="white", type="png")
heatmap.2(res.hm, col=greenred(200), key=T, keysize=1, density.info="none", trace="none",cexCol=2, cexRow=1, na.rm=T, main=paste("Project: ", projectName, sep=""), margin=c(15,15), ColSideColors=columncolors)
dev.off()

#plot heatmap of top 100 differentially expressed genes
top100id = order(res$padj, res$pvalue)[1:100]
top100 = res[top100id,]
select = which(rownames(assay(vsd)) %in% rownames(top100))
res.hm = assay(vsd)[select,]
column.conditions <- unlist(lapply(conditions, condition.map))
column.names <- paste(column.conditions, colnames(res.hm), sep=":")
colnames(res.hm) <- column.names
color.map <- function(mol.biol) { if (startsWith(mol.biol, "case", trim=T, ignore.case=T)) "grey50" else "grey80" }
columncolors <- unlist(lapply(colnames(res.hm), color.map))
Cairo(file=paste(outputFolder, "top100_de_HeatMap.png", sep="/"), width=800, height=1200, bg="white", canvas="white", type="png")
heatmap.2(res.hm, col=greenred(200), key=T, keysize=1, density.info="none", trace="none",cexCol=2, cexRow=1, na.rm=T, main=paste("Project: ", projectName, sep=""), margin=c(15,15), ColSideColors=columncolors)
dev.off()

#plot histogramm of pval's
Cairo(file=paste(outputFolder, "pValueHistogram.png", sep="/"), widht=640, height=480, bg="white", canvas="white", type="png")
hist(res$pvalue, breaks=100, col="skyblue", border="slateblue", main=paste("Project: ", projectName, sep=""), xlab="p-value")
dev.off()
Cairo(file=paste(outputFolder, "pAdjHistogram.png", sep="/"), widht=640, height=480, bg="white", canvas="white", type="png")
hist(res$padj, breaks=100, col="skyblue", border="slateblue", main=paste("Project: ", projectName, sep=""), xlab="p-adjusted")
dev.off()

#volcano plot (taken from Tutorial: Statistical analysis of RNAseq data)
dds.volc = DESeq(dds, betaPrior=F)
res.volc = results(dds.volc)
tab = data.frame(logFC=res.volc$log2FoldChange, negLogPval = -log10(res.volc$pvalue), rownames=rownames(res.volc) )
Cairo(file=paste(outputFolder, "volcanoPlot.png", sep="/"), widht=640, height=1080, bg="white", canvas="white", type="png")
par(mar = c(5, 4, 4, 4))
plot(tab$logFC, tab$negLogPval, pch = 16, cex = 0.6, xlab = expression(log[2]~fold~change), ylab = expression(-log[10]~pvalue))
lfc=2
pval=0.01
signGenes = (abs(tab$logFC) > lfc & tab$negLogPval > -log10(pval))

points(tab[signGenes, ], pch = 16, cex = 0.8, col = "red")
abline(h = -log10(pval), col = "green3", lty = 2)
abline(v = c(-lfc, lfc), col = "blue", lty = 2)
mtext(paste("pval =", pval), side = 4, at = -log10(pval), cex = 0.8, line = 0.5, las = 1)
mtext(c(paste("-", lfc, "fold"), paste("+", lfc, "fold")), side = 3, at = c(-lfc, lfc),	cex = 0.8, line = 0.5)
text(tab[signGenes,"logFC"], tab[signGenes, "negLogPval"], labels=tab[signGenes, "rownames"], cex= 0.7, pos=3, col="red")
dev.off()

#volcano plot (taken from Tutorial: Statistical analysis of RNAseq data)
dds.volc = dds
res.volc = res
tab = data.frame(logFC=res.volc$log2FoldChange, negLogPval = -log10(res.volc$pvalue), rownames=rownames(res.volc) )
Cairo(file=paste(outputFolder, "volcanoPlot2.png", sep="/"), widht=640, height=1080, bg="white", canvas="white", type="png")
par(mar = c(5, 4, 4, 4))
plot(tab$logFC, tab$negLogPval, pch = 16, cex = 0.6, xlab = expression(log[2]~fold~change), ylab = expression(-log[10]~pvalue))
lfc=2
pval=0.01
signGenes = (abs(tab$logFC) > lfc & tab$negLogPval > -log10(pval))
points(tab[signGenes, ], pch = 16, cex = 0.8, col = "red")
abline(h = -log10(pval), col = "green3", lty = 2)
abline(v = c(-lfc, lfc), col = "blue", lty = 2)
mtext(paste("pval =", pval), side = 4, at = -log10(pval), cex = 0.8, line = 0.5, las = 1)
mtext(c(paste("-", lfc, "fold"), paste("+", lfc, "fold")), side = 3, at = c(-lfc, lfc),	cex = 0.8, line = 0.5)
text(tab[signGenes,"logFC"], tab[signGenes, "negLogPval"], labels=tab[signGenes, "rownames"], cex= 0.7, pos=3, col="red")
dev.off()

#create heatmap of all reads
if (noCreateHeatMap == 0) {
	tryCatch({
		mat=as.table(counts(dds))
		#get rid of NAs
		mat[is.na(mat)] <- 0
		min<-apply(mat, 1,min)
		max<-apply(mat,1,max)
		smat <- mat[order(max - min, decreasing=TRUE),][1:20000,]
		#cluster on correlation
		cdist <- as.dist(1-cor(t(smat)))
		hc <- hclust(cdist, "average")
		color_palette <- colorRampPalette(c("blue", "white", "red"))(n = 300)
		Cairo(file=paste(outputFolder, "heatmap_all_rawreads.png", sep="/"), width=800, height=1200, bg="white", canvas="white", type="png")
		heatmap.2(as.matrix(smat),
						Rowv=as.dendrogram(hc),
						Colv=T,
						col=color_palette, 
						scale="row", 
						dendrogram="both",
						key=T, 
						keysize=.5, 
						density.info="none", 
						trace="none",
						na.rm=T, 
						cexRow=1,
						cexCol=1,
						labRow=c(""))
		dev.off()
				
	},  interrupt = function(ex) {
		write.csv(mat,file=paste(outputFolder, "heatmapCreationInterruption.ERR", sep="/"), sep="\t");
	}, error = function(ex) {
		write.csv(mat,file=paste(outputFolder, "heatmapCreationError.ERR", sep="/"), sep="\t");
	})
}

#write DE results to the file system
res[which(res$log2FoldChange<0), "regulation"] <- paste("down")
res[which(res$log2FoldChange>0), "regulation"] <- paste("up")
res[which(res$log2FolgChange==0), "regulation"] <- paste("none")
write.csv(na.omit(res[order(res$padj, res$pvalue),c(1,2,3,5,6,7)]), file=paste(outputFolder, paste(projectName, "csv", sep="."), sep="/"), sep="\t")
write(row.names(res[which(res$padj<0.01),]), file=paste(outputFolder, "sigGenes.out", sep="/"));

#save this session data
save.image(paste(outputFolder, "session.RData", sep="/"))
