setRepositories(1,2)

library(GEOquery)
library(limma)
library(Biobase)
library(pheatmap)
library(dplyr)
library(magrittr)
library(plyr)

series = "GSE48558"
platform = "GPL6244"
gset <- getGEO(series, GSEMatrix =TRUE,  AnnotGPL=TRUE , destdir = "data")
     
if (length(gset) > 1) {
  idx <- grep(platform, attr(gset, "names"))
}else {
  idx <- 1}


gset <- gset[[idx]]
source_name <- pData(phenoData(gset))$source_name_ch1

source_name.CD34 = grepl("CD34",source_name)
source_name.aml = grepl("AML Patient",source_name)


fac = function(cd34,aml,name)
{
 ifelse(cd34,"CD34",ifelse(aml,"AML",name))
}
f <- factor(fac(c(source_name.CD34),c(source_name.aml),source_name),levels=c("AML","CD34",levels(factor(source_name))))

which_select = which(f == "B Cells" | f == "T Cells"
                     | f == "AML"| f == "CD34"|
                       f == "Granulocytes" | f == "Monocytes")

gset <- gset[,which_select ]
gset$description <- as.character(factor(f[which_select]))
gset$description <- ifelse(gset$description == "B Cells","B_Cells",gset$description)
gset$description <- ifelse(gset$description == "T Cells","T_Cells",gset$description)
gset$description <- factor(gset$description )
design<- model.matrix(~ description + 0, gset)
colnames(design) <- levels(gset$description)
fit = lmFit(gset, design) 
cont.matrix = makeContrasts(AML - CD34, levels = design) 
fit2 = contrasts.fit(fit, cont.matrix) 
fit2 = eBayes(fit2, 0.01)
compare_data =topTable(fit2,number=Inf,adjust="fdr")

compare_data = subset(compare_data,select=c("P.Value","adj.P.Val","logFC","Gene.symbol"))

inc =  compare_data %>% dplyr::filter( logFC > 1 & adj.P.Val < 0.05 ) 
inc = inc[order(inc$logFC),]
n_inc = unique( as.character(strsplit2( (inc$Gene.symbol[(!is.na(inc$Gene.symbol) | inc$Gene.symbol != "")]),"///")))
write.table(n_inc, file = "result/inc.txt",col.names = F,row.names = F,quote = F)


dec =  compare_data %>% dplyr::filter( logFC < -1 & adj.P.Val < 0.05 & (!is.na(Gene.symbol) | Gene.symbol != "")) 
dec = dec[order(dec$logFC,decreasing = T),]
n_dec = unique( as.character(strsplit2( (dec$Gene.symbol[(!is.na(dec$Gene.symbol) | dec$Gene.symbol != "")]),"///")))
write.table(n_dec, file = "result/dec.txt",col.names = F,row.names = F,quote = F)


