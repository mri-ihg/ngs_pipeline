library("edgeR")
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
doRuv <- myargs[6]

countTable = read.table(countFile, header=T, row.names=1, check.names=F, quote="\"")
fpkmTable = read.table(fpkmFile, header=T, row.names=1, check.names=F, quote="\"")
countTable = countTable[1:(nrow(countTable)-5),] 			#delete last 5 lines; includes unneeded information
conditionTable = read.table(conditionFile, header=F, check.names=F, row.names=1)
names(conditionTable)[1] <- paste("condition")
conditions = c()
for (i in 1:ncol(countTable)) {
	conditions <- c(conditions, conditionTable[which(rownames(conditionTable)==colnames(countTable)[i]), 1])
}
formatted.conditionTable <- transform(conditionTable, condition=as.factor(condition))
save.image(paste(outputFolder, "session.RData", sep="/"))

#edgeR stuff
group = formatted.conditionTable$condition
keep <- rowSums(cpm(countTable)>1) >= 3
y <- countTable[keep,]
y.countTable <- y

if (doRuv == 0) {
	design <- model.matrix(~group)
	y <- DGEList(counts=y)
	y <- estimateGLMCommonDisp(y, design)
	y <- estimateGLMTrendedDisp(y, design)
	y <- estimateGLMTagwiseDisp(y, design)
	fit <- glmFit(y, design)
	lrt <- glmLRT(fit, coef=2)
	tt <- topTags(lrt, nrow(lrt))
	save.image()
} else {
	library(RUVSeq)
	set <- newSeqExpressionSet(as.matrix(y.countTable),	phenoData = data.frame(group, row.names=colnames(y.countTable)))
	design <- model.matrix(~group, data=pData(set))
	
	y.ruv <- DGEList(counts=counts(set), group=group)
	y.ruv <- calcNormFactors(y.ruv, method="upperquartile")
	y.ruv <- estimateGLMCommonDisp(y.ruv, design)
	y.ruv <- estimateGLMTagwiseDisp(y.ruv, design)
	fit.ruv <- glmFit(y.ruv, design)
	res.ruv <- residuals(fit.ruv, type="deviance")
	set.ruv <- RUVr(set, row.names(y.countTable), k=1, res.ruv)

	design <- model.matrix( ~ group + W_1, data=pData(set.ruv))
	y <- DGEList(counts=counts(set.ruv), group=group)
	y <- calcNormFactors(y, method="upperquartile")
	y <- estimateGLMCommonDisp(y, design)
	y <- estimateGLMTagwiseDisp(y, design)
	fit <- glmFit(y, design)
	lrt <- glmLRT(fit, coef=2)
	tt <- topTags(lrt, nrow(lrt))
	save.image()
}

###PCA####
Cairo(file=paste(outputFolder, "PCA.png", sep="/"), width=800, height=800, bg="white", canvas="white", type="png")
plotMDS(y)
dev.off()
###BCV####
Cairo(file=paste(outputFolder, "BCV.png", sep="/"), width=800, height=800, bg="white", canvas="white", type="png")
plotBCV(y)
dev.off()

###logFC###
de <- decideTestsDGE(lrt)
detags <- rownames(y)[as.logical(de)]
Cairo(file=paste(outputFolder, "logFC.png", sep="/"), width=800, height=800, bg="white", canvas="white", type="png")
plotSmear(lrt, de.tags=detags)
abline(h=c(-1, 1), col="blue")
dev.off()

#write result table to csv
tt.df <- as.data.frame(tt)
tt.df[which(tt.df$logFC<0), "regulation"] <- paste("down")
tt.df[which(tt.df$logFC>0), "regulation"] <- paste("up")
tt.df[which(tt.df$logFC==0), "regulation"] <- paste("none")
write.csv(tt.df, file=paste(outputFolder, paste(projectName, "csv", sep="."), sep="/"), sep="\t")

#save this session data
save.image(paste(outputFolder, "session.RData", sep="/"))