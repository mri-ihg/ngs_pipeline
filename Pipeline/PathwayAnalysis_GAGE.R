#######
##gage
#######
library(gage)
library(gageData)
library(pathview)
data(kegg.gs)
data(go.sets.mm)
data(go.subs.mm)



#TRANSLATE ENTREZ2SYMBOL
library(mygene)
xli <- rownames(res)
xli.translated <- queryMany(xli, scopes="symbol", fields="entrezgene", species="human")

xli.transl.df <- as.data.frame(xli.translated)

#remove rows where "entrezid" == NA
xli.transl.df.rect <- xli.transl.df[complete.cases(xli.transl.df[,"entrezgene"]),]

res.df <- as.data.frame(res)

res.transl <- merge(res.df, xli.transl.df.rect, by.x="row.names", by.y="query")
res.transl.rectified <- res.transl[which(!is.na(res.transl$entrezgene)),]
deseq2.fc <- res.transl.rectified$log2FoldChange
names(deseq2.fc) <- res.transl.rectified$entrezgene
exp.fc=deseq2.fc
out.suffix="deseq2"

##KEGG Pathway
fc.kegg.p <- gage(exp.fc, gsets = kegg.gs, ref=NULL, samp=NULL, same.dir=F)
sel.g <- fc.kegg.p$greater[, "q.val"] < 0.1 & !is.na(fc.kegg.p$greater[, "q.val"])
path.ids.g <- rownames(fc.kegg.p$greater)[sel.g]
sel.l <- fc.kegg.p$less[, "q.val"] < 0.1 & !is.na(fc.kegg.p$less[,"q.val"])
path.ids.l <- rownames(fc.kegg.p$less)[sel.l]

path.ids.g.final <- substr(path.ids.g, 1, 8)
path.ids.l.final <- substr(path.ids.l, 1, 8)

pv.out.list.g <- sapply(path.ids.g.final, function(pid) pathview(gene.data = exp.fc, pathway.id = pid,species = "hsa", out.suffix="deseq2.greater"))
pv.out.list.l <- sapply(path.ids.l.final, function(pid) pathview(gene.data = exp.fc, pathway.id = pid,species = "hsa", out.suffix="deseq2.less"))


###with counts
###BETTER!!!
cnts <- counts(dds, normalize=T)
cnts.df <- as.data.frame(cnts)

cnts.transl <- merge(cnts.df, xli.transl.df.rect, by.x="row.names", by.y="query")
cnts.transl.rectified <- cnts.transl[which(!is.na(cnts.transl$entrezgene)),]
cnts <- cnts.transl.rectified
cnts <- cnts[!duplicated(cnts$entrezgene),]

#TODO: repair this duplicate issue in gtf file
rownames(cnts) <- cnts$entrezgene
cnts.mtrx <- as.matrix(cnts[,c(2:7)])
cntrl.ids <- c(2,4,6)
case.ids <- c(1,3,5)
cnts.kegg.p <- gage(cnts.mtrx, gsets=kegg.gs, ref=cntrl.ids, samp=case.ids, compare="paired")

sel.g <- cnts.kegg.p$greater[, "q.val"] < 0.1 & !is.na(cnts.kegg.p$greater[, "q.val"])
path.ids.g <- rownames(cnts.kegg.p$greater)[sel.g]
sel.l <- cnts.kegg.p$less[, "q.val"] < 0.1 & !is.na(cnts.kegg.p$less[,"q.val"])
path.ids.l <- rownames(cnts.kegg.p$less)[sel.l]

path.ids.g.final <- substr(path.ids.g, 1, 8)
path.ids.l.final <- substr(path.ids.l, 1, 8)

pv.out.list.up <- sapply(path.ids.g.final, function(pid) pathview(gene.data = exp.fc, pathway.id = pid,species = "hsa", out.suffix="deseq2.up"))
pv.out.list.down <- sapply(path.ids.l.final, function(pid) pathview(gene.data = exp.fc, pathway.id = pid,species = "hsa", out.suffix="deseq2.down"))

#write results
write.table(cnts.kegg.p$greater[complete.cases(cnts.kegg.p$greater[,"q.val"]),c(3:4)], file="KEGG.pathways.up.csv", sep="\t")
write.table(cnts.kegg.p$less[complete.cases(cnts.kegg.p$less[,"q.val"]),c(3:4)], file="KEGG.pathways.down.csv", sep="\t")


