setRepositories(1,2)
install.packages("GEOquery")
install.packages("limma")
install.packages("umap")
install.packages("pheatmap")
install.packages("ggplot2")
install.packages("gplots")
install.packages("reshape2")
install.packages("plyr")
install.packages("stringi")
install.packages("Biobase")
install.packages('MASS')
install.packages('umap')
install.packages('qkerntool')
install.packages('M3C')
install.packages('gridExtra')
install.packages("annotate")
install.packages("hugene10stprobeset.db")

library(GEOquery)
library(limma)
library(umap)
library(maptools)  
library(pheatmap)
library(gplots)
library(ggplot2)
library(reshape2)
library(plyr)
library(Biobase)
library(MASS)
library(qkerntool)
library(umap)
library(M3C)
library(dplyr)
library(gridExtra)
library(annotate)
library(hugene10stprobeset.db)

setwd(getwd())

# load series and platform data from GEO
series = "GSE48558"
platform = "GPL6244"
gset <- getGEO(series, GSEMatrix =TRUE, getGPL=TRUE,destdir = "data")

if (length(gset) > 1) {
   idx <- grep(platform, attr(gset, "names"))
}else {
  idx <- 1}


gset <- gset[[idx]]


source_name <- pData(phenoData(gset))$source_name_ch1
source_name.aml = grepl("AML Patient",source_name)
source_name.raw = gset[,source_name.aml]

ph <-  pData(phenoData(gset))$characteristics_ch1
ph.normal <-  grepl("Normal",ph)
ph.raw = gset[,ph.normal]

ex <- exprs(gset)
title <- paste (series, "/", platform, sep ="")

source_name.aml = exprs(source_name.raw)
ph.normal = exprs(ph.raw)


colors = c(rep("green",ncol(ph.normal)),rep("red",ncol(source_name.aml)))
group = c(rep("healthy",ncol(ph.normal)),rep("aml",ncol(source_name.aml)))


data.scaled = t(scale(t(data.normalized)))
data = cbind(ph.normal,source_name.aml)
data.normalized = normalizeQuantiles(data)



pdf("result/boxplot_health.pdf",width = 32,height = 16)
boxplot(data, boxwex=0.7, notch=T, main=title, outline=FALSE, las=2,col=colors)
dev.off()

pdf("result/boxplot_normalized.pdf",width = 32,height = 16)
boxplot(data.normalized, boxwex=0.7, notch=T, main=title, outline=FALSE, las=2,col=colors)
dev.off()


newColName = function(names,matrix){
  paste(names,seq(1,ncol(matrix)),"")
}

healthy.new_name = exprs(ph.raw)
colnames(healthy.new_name) = newColName(ph.raw$source_name_ch1,healthy.new_name)
healthy.new_name = healthy.new_name[, sort(colnames(healthy.new_name))]


aml.new_name = exprs(source_name.raw)
colnames(aml.new_name) = rep("aml",ncol(aml.new_name))

data.cor = cbind(healthy.new_name,aml.new_name)
data.cor = cor(data.cor)


pdf("result/corheatmap_data.pdf",width = 32,height = 32)
pheatmap(data.cor,main=title, display_numbers = TRUE,border_color = "black")
dev.off()

pdf("result/corheatmap_data_no_cluster.pdf",width = 32,height = 32)
pheatmap(data.cor,main=title, display_numbers = TRUE,border_color = "black",cluster_rows = F,cluster_cols = F)
dev.off()


data.umap = umap( data.scaled)
data.umap_col = data.umap[["data"]]
v1 = data.umap_col[,1]
v2 = data.umap_col[,2]

pdf("result/umap.pdf",width = 8,height = 6)
plot(v1,v2, main = paste(title,"/","umap",sep = ""), col = colors, pch = 19)
legend("topright", legend=c("healthy", "AML"),
       col=c("green", "red"), pch= 19, cex=0.8)
dev.off()

data.pca = prcomp(data.scaled,scale. = F)

pcr = data.frame(data.pca$rotation[,1:3], group = c(rep("normal",ncol(ph.normal)),rep("aml",ncol(source_name.aml))))

pdf("result/geom_point_pca12.pdf",width = 6,height = 4)
ggplot(pcr,aes(PC1,PC2,color = group)) + geom_point(size = 1) + theme_bw()
dev.off()

pdf("result/geom_point_pca13.pdf",width = 6,height = 4)
ggplot(pcr,aes(PC1,PC3,color = group)) + geom_point(size = 1) + theme_bw()
dev.off()

pdf("result/geom_point_pca23.pdf",width = 6,height = 4)
ggplot(pcr,aes(PC2,PC3,color = group)) + geom_point(size = 1) + theme_bw()
dev.off()

data.distance = dist(t(data.scaled))
data.mds_classical  = cmdscale(data.distance)
v1 = data.mds_classical[,1]
v2 = data.mds_classical[,2]

pdf("result/mds_classical.pdf",width = 8,height = 6)
plot(v1,v2, main = paste (title, "/", "mds", sep =""), col = colors, pch = 19)
legend("topleft", legend=c("healthy", "AML"),
       col=c("green", "red"),pch = 19 )
dev.off()


pdf("result/tSNE_Samples.pdf",width = 8,height = 6)
tsne(data,labels=group)
dev.off()

#phase 2

gset <- getGEO(series, GSEMatrix =TRUE, getGPL=TRUE,destdir = "data")

if (length(gset) > 1) {
  idx <- grep(platform, attr(gset, "names"))
}else {
  idx <- 1}


gset <- gset[[idx]]

source_name <- pData(phenoData(gset))$source_name_ch1
source_name.CD34 = grepl("CD34",source_name)
source_name.CD34_raw = gset[,source_name.CD34]

source_name.aml = grepl("AML Patient",source_name)

aml = exprs(source_name.raw)
CD34 = exprs(source_name.CD34_raw)
fac = function(cd34,aml)
{
  ifelse(cd34,"CD34",ifelse(aml,"AML","none"))
}
f <- factor(fac(c(source_name.CD34),c(source_name.aml)),levels=c("AML","CD34","none"))
gset$description <- f
design<- model.matrix(~ description + 0, gset)
colnames(design) <- levels(f)
fit = lmFit(gset, design) 
cont.matrix = makeContrasts(AML - CD34, levels = design) 
fit2 = contrasts.fit(fit, cont.matrix) 
fit2 = eBayes(fit2, 0.01)
compare_data =topTable(fit2,number=Inf,adjust="fdr")
compare_data = subset(compare_data,select=c("ID","P.Value","adj.P.Val","logFC","mrna_assignment"))

inc =  filter(compare_data, logFC > 1 & adj.P.Val < 0.05) 
inc = inc[order(inc$logFC),]
n_inc = unique( as.character(strsplit2( (inc$mrna_assignmentl),"///")))
write.table(rownames(inc), file = "result/inc.txt",col.names = F,row.names = F,quote = F)


dec =  filter(compare_data, logFC < -1 & adj.P.Val < 0.05) 
dec = dec[order(dec$logFC,decreasing = T),]
write.table(rownames(dec), file = "result/dec.txt",col.names = F,row.names = F,quote = F)







