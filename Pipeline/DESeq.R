#run all neccessary steps for RNA-seq analysis
#main component of this script is the DESeq package which handles most of the RNA-seq analysis

#DEPRECATED!

library("DESeq")
library("Cairo")
library("RColorBrewer")
library("pvclust")
library("ggplot2")
library("gplots")
library("gdata")
library("genefilter")
library("scatterplot3d")


myargs       <- commandArgs(trailingOnly = TRUE)
countFile    <- myargs[1]
conditionFile <- myargs[2]
outputFolder <- myargs[3]
projectName  <- myargs[4]
fpkmFile <- myargs[5]

countTable = read.table(countFile, header=T, row.names=1, check.names=F, quote="\"")
fpkmTable = read.table(fpkmFile, header=T, row.names=1, check.names=F, quote="\"")
countTable.waste = countTable[nrow(countTable)-4:nrow(countTable),]
countTable = countTable[1:(nrow(countTable)-5),] 			#delete last 5 lines; includes unneeded information
conditionTable = read.table(conditionFile, header=F, check.names=F)
conditions = c()
for (i in 1:ncol(countTable)) {
	conditions <- c(conditions, conditionTable[which(conditionTable[,1]==colnames(countTable)[i]), 2])
}
save.image()
cds = newCountDataSet(countTable, conditions)
cds = estimateSizeFactors(cds)

if (length(which(conditions == "0")) < 2 || length(which(conditions == "1")) < 2) {			#check if there are enough sample for cases and controls, otherwise estimateDispersions() throws an error without the blind argument
	#it's possible that "cds = estimateDispersions(cds, method='blind')" throws an error, then one has to use fitType="local" as well
	tryCatch({
				cds = estimateDispersions(cds, method='blind')
			},
			error = function (ex) {
				cat("OUTPUT: error while cds = estimateDispersions(cds, method='blind')\n")
				cat("OUTPUT: trying to add fitType=local\n")
			},
			finally = {
				cds = estimateDispersions(cds, method='blind', fitType='local')
			})
	
} else {
	cds = estimateDispersions(cds)
}
Cairo(file=paste(outputFolder, "EstimatedDispersion.png", sep="/"), widht=800, height=480, bg="white", canvas="white", type="png")
plotDispEsts(cds, main=paste("Project: ", projectName, sep=""))
dev.off()

#DESeq
vsd = varianceStabilizingTransformation(cds)

#plot heatmap of top 100 rowMeans Variance Transcripts
distCor <- function(x) as.dist(1-cor(t(x)))				#for heatmap.2: distfun=distCor
hclustAvg <- function(x) hclust(x, method="average")	#for heatmap.2: hclustfun=hclustAvg

select = order(rowMeans(exprs(vsd)), decreasing=T)[1:100]
res.hm = exprs(vsd)[select,]
condition.map <- function(cond.biol) { if (cond.biol == 1) "case" else "control" }
column.conditions <- unlist(lapply(conditions, condition.map))
column.names <- paste(column.conditions, colnames(res.hm), sep=":")
colnames(res.hm) <- column.names
color.map <- function(mol.biol) { if (startsWith(mol.biol, "case", trim=T, ignore.case=T)) "grey50" else "grey80" }
columncolors <- unlist(lapply(colnames(res.hm), color.map))
Cairo(file=paste(outputFolder, "top1000_rowMeans_HeatMap.png", sep="/"), width=800, height=1200, bg="white", canvas="white", type="png")
heatmap.2(res.hm, col=greenred(200), key=T, keysize=1, density.info="none", trace="none",cexCol=2, cexRow=1, na.rm=T, main=paste("Project: ", projectName, sep=""), margin=c(15,15), ColSideColors=columncolors)
dev.off()


#heatmap (sample-to-sample)
dists= dist(t(exprs(vsd)))
mat = as.matrix(dists)
rownames(mat) = colnames(mat) = paste(rownames(pData(cds)), (pData(cds))$condition, sep=":")
Cairo(file=paste(outputFolder, "Sample2Sample.png", sep="/"), width=800, height=800, bg="white", canvas="white", type="png")
heatmap.2(mat, trace="none", col=heat.colors, margin=c(15,15), main=paste("Project: ", projectName, sep=""))
dev.off()

#Principal component analysis
pca = prcomp(t(exprs(vsd)[,]))
pcax = as.data.frame(pca$x)
pcax[which(conditions(cds)==1), "gn"] = "case"
pcax[which(conditions(cds)==0), "gn"] = "control"

Cairo(file=paste(outputFolder, "PCA.png", sep="/"), width=640, height=480, bg="white", canvas="white", type="png")
ggplot(as.data.frame(pca$x), aes(PC1, PC2, color=pcax$gn)) 	+ 
		scale_color_discrete(name="Samples") + 
		geom_point(size=3) +
		geom_text(aes(label=rownames(pca$x)), hjust=0, vjust=0, show_guide=F, size=4) +
		coord_cartesian(xlim=c(min((pca$x)[,"PC1"]) - (2/100*(max((pca$x)[,"PC1"]) - min((pca$x)[,"PC1"]))), max((pca$x)[,"PC1"])+ (20/100*(max((pca$x)[,"PC1"]) - min((pca$x)[,"PC1"])))))	+
		ggtitle(paste("Project: ", projectName, sep=""))
dev.off()

#DESeq: DE analysis
res = nbinomTest(cds, "0", "1")
#plotMA
Cairo(file=paste(outputFolder, "MAPlot.png", sep="/"), widht=800, height=480, bg="white", canvas="white", type="png")
plotMA(res, main=paste("Project: ", projectName, sep=""))
dev.off()

#plot heatmap of top 100 differentially expressed genes
top100id = order(res$padj, res$pval)[1:100]
top100 = res[top100id,]
select = which(rownames(exprs(vsd)) %in% top100$id)
res.hm = exprs(vsd)[select,]
condition.map <- function(cond.biol) { if (cond.biol == 1) "case" else "control" }
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
hist(res$pval, breaks=100, col="skyblue", border="slateblue", main=paste("Project: ", projectName, sep=""), xlab="p-value")
dev.off()

#volcano plot (taken from Tutorial: Statistical analysis of RNAseq data)
tab = data.frame(logFC=res$log2FoldChange, negLogPval = -log10(res$pval), rownames=res$id)
Cairo(file=paste(outputFolder, "volcanoPlot.png", sep="/"), widht=640, height=1080, bg="white", canvas="white", type="png")
par(mar = c(5, 4, 4, 4))
plot(tab$logFC, tab$negLogPval, pch = 16, cex = 0.6, xlab = expression(log[2]~fold~change),ylab = expression(-log[10]~pvalue))
lfc=1.5
pval=0.01
signGenes = (abs(tab$logFC) > lfc & tab$negLogPval > -log10(pval))
points(tab[signGenes, ], pch = 16, cex = 0.8, col = "red")
abline(h = -log10(pval), col = "green", lty = 1)
abline(v = c(-lfc, lfc), col = "blue", lty = 1)
mtext(paste("pval =", pval), side = 4, at = -log10(pval), cex = 0.8, line = 0.5, las = 1)
mtext(c(paste("-", lfc, "fold"), paste("+", lfc, "fold")), side = 3, at = c(-lfc, lfc),	cex = 0.8, line = 0.5)
text(tab[signGenes,"logFC"], tab[signGenes, "negLogPval"], labels=tab[signGenes, "rownames"], cex= 0.7, pos=3, col="red")
dev.off()



resSig = na.omit(res[res$padj < 0.1,])
if (dim(resSig)[1] > 1) {  #check if generally significant expression is present
	if (dim(resSig)[1] > 100) {resSigSelect = head(resSig[order(resSig$pval),], n=100L)} else {resSigSelect = head(resSig[order(resSig$pval),], n=(dim(resSig)[1]))}
	vsd=varianceStabilizingTransformation(cds)
	select = which(rownames(counts(cds)) %in% resSigSelect$id)
	gexprs=exprs(vsd)
	Cairo(file=paste(outputFolder, "heatmap.2.pval_z_outliers_excluded.png", sep="/"), width=580, height=840, bg="white", canvas="white", type="png")
	heatmap.2(gexprs[which(rownames(gexprs) %in% resSigSelect$id),], col="heat.colors", trace="none", margin=c(10,6), main=paste("Project: ", "projectName", sep=""))
	dev.off()
}

#write DE results to the file system
res[which(res$log2FoldChange<0), "regulation"] <- paste("down")
res[which(res$log2FoldChange>0), "regulation"] <- paste("up")
res[which(res$log2FolgChange==0), "regulation"] <- paste("none")
write.csv(na.omit(res[order(res$padj, res$pval),]), file=paste(outputFolder, paste(projectName, "csv", sep="."), sep="/"))


#save this session data
save.image(paste(outputFolder, "session.RData", sep="/"))



countTable.stack <- stack(as.data.frame(counts(cds)))
countTable.stack[,"log2values"] <- log2(countTable.stack$values)
Cairo(file=paste(outputFolder, "density_unnormalized.png", sep="/"), width=840, height=840, bg="white", canvas="white", type="png")
ggplot(countTable.stack, aes(x=log2values) ) + geom_density(aes(group=ind, colour=ind, fill=ind), alpha=.3) + xlim(-1, 20)
dev.off()

cds.norm.stack <- stack(as.data.frame(counts(cds, normalize=T)))
cds.norm.stack[,"log2values"] <- log2(cds.norm.stack$values)

Cairo(file=paste(outputFolder, "density_normalized.png", sep="/"), width=840, height=840, bg="white", canvas="white", type="png")
ggplot(cds.norm.stack, aes(x=log2values) ) + geom_density(aes(group=ind, colour=ind, fill=ind), alpha=.3) + xlim(-1, 20)
dev.off()

save.image(paste(outputFolder, "session.RData", sep="/"))