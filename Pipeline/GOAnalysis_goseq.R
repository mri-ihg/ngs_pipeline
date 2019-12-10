library(goseq)
library(GO.db)
library(DESeq2)

myargs       <- commandArgs(trailingOnly = TRUE)
goseqOutputFolder <- myargs[1]
sessionFile  <- myargs[2]
settings <- myargs[3]

load(sessionFile)

assayed.genes <- rownames(na.omit(res))
de.genes.up   <- rownames(na.omit(res[which(res$padj < 0.1 & res$log2FoldChange > 0),]))
de.genes.down <- rownames(na.omit(res[which(res$padj < 0.1 & res$log2FoldChange < 0),]))

gene.vector.up   <- as.integer(assayed.genes %in% de.genes.up)
gene.vector.down <- as.integer(assayed.genes %in% de.genes.down)
	
names(gene.vector.up)   <- assayed.genes
names(gene.vector.down) <- assayed.genes

if (length(grep("^hg19", settings))>0) {
	pwf.up   <- nullp(gene.vector.up, "hg19", "geneSymbol")
	pwf.down <- nullp(gene.vector.down, "hg19", "geneSymbol")
	
	GO.wall.up = goseq(pwf.up, "hg19", "geneSymbol")
	GO.samp.up = goseq(pwf.up, "hg19", "geneSymbol", method="Sampling", repcnt=1000)
	
	GO.wall.down = goseq(pwf.down, "hg19", "geneSymbol")
	GO.samp.down = goseq(pwf.down, "hg19", "geneSymbol", method="Sampling", repcnt=1000)
} else {
	pwf.up   <- nullp(gene.vector.up, "mm9", "geneSymbol")
	pwf.down <- nullp(gene.vector.down, "mm9", "geneSymbol")
		
	GO.wall.up = goseq(pwf.up, "mm9", "geneSymbol")
	GO.samp.up = goseq(pwf.up, "mm9", "geneSymbol", method="Sampling", repcnt=1000)
	
	GO.wall.down = goseq(pwf.down, "mm9", "geneSymbol")
	GO.samp.down = goseq(pwf.down, "mm9", "geneSymbol", method="Sampling", repcnt=1000)
}

#correct for multiple testing
GO.wall.up[,"padj"]   <- p.adjust(GO.wall.up$over_represented_pvalue, method="BH")
GO.samp.up[,"padj"]   <- p.adjust(GO.samp.up$over_represented_pvalue, method="BH")

GO.wall.down[,"padj"] <- p.adjust(GO.wall.down$over_represented_pvalue, method="BH")
GO.samp.down[,"padj"] <- p.adjust(GO.samp.down$over_represented_pvalue, method="BH")


paste(outputFolder, paste(projectName, "downregulated","tsv", sep="."), sep="/")

write.table(na.omit(GO.wall.up[order(GO.wall.up$padj, GO.wall.up$over_represented_pvalue),c(1,6,7,2,8)]), file=paste(goseqOutputFolder, paste(projectName, "upregulated","tsv", sep="."), sep="/"), sep="\t")
write.table(na.omit(GO.samp.up[order(GO.samp.up$padj, GO.samp.up$over_represented_pvalue),c(1,6,7,2,8)]), file=paste(goseqOutputFolder, paste(projectName, "upregulated","samp","tsv", sep="."), sep="/"), sep="\t")

write.table(na.omit(GO.wall.down[order(GO.wall.down$padj, GO.wall.down$over_represented_pvalue),c(1,6,7,2,8)]), file=paste(goseqOutputFolder, paste(projectName, "downregulated","tsv", sep="."), sep="/"), sep="\t")
write.table(na.omit(GO.samp.down[order(GO.samp.down$padj, GO.samp.down$over_represented_pvalue),c(1,6,7,2,8)]), file=paste(goseqOutputFolder, paste(projectName, "downregulated","samp","tsv", sep="."), sep="/"), sep="\t")


save.image(paste(goseqOutputFolder, "session.RData", sep="/"))