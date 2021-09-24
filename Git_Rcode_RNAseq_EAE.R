
library(DBI)
library(RODBC)
library(lattice)
library(latticeExtra)
library(ggplot2) 
library(gplots) 
library(RColorBrewer) 
library(xlsx)
library(devtools)
library(factoextra)
library(caret)  
library(Rtsne)
library(reshape2)
library(DESeq2)
library(BiocStyle)
library(rmarkdown)
library(geneplotter)
library(plyr)
library(LSD)
library(stringr)
library(genefilter)
library(biomaRt)
library(dplyr)
library(EDASeq)
library(fdrtool)
library(org.Mm.eg.db)
library(limma)
library(edgeR)
library(ggrepel)
library(FactoMineR)
library(ggbiplot)
library(gridExtra)
library(plot3D)
library(EnhancedVolcano)
library(readxl)


close(con)  # close the connection
rm(list = ls(all.names = TRUE)) # clear workspace
dev.off() # clear all plots

########## Build connection to MouseDB
con = odbcDriverConnect('driver={SQL Server};server=localhost\\SQLEXPRESS;database=MouseDB;trusted_connection=true')


################### Retrieve the raw data from SQL MouseDB database

data_raw <- sqlQuery(con,  "   select * from omics_results_raw where omics_name='RNA-seq_Martina_all_tissue' and sample_name like  '%Bonemarrow%'")

data_raw=dcast (data_raw, gene_id_name ~ sample_name, value.var='raw_value' )
cts=as.matrix(data_raw[,-1])

rownames(cts) <- data_raw$gene_id_name
cts[is.na(cts)] <- 0

fdr = 5
min_counts=1
bcv=0.1

keep <- filterByExpr(cts,min.count = 5, min.total.count = 4)
res <- cts[keep,]
cond <- factor(c("CE","CE","CE","RT","RT","RT"))


y <- DGEList(res, group= cond,  genes= row.names(res) )

y <- calcNormFactors(y)
options(digits=3)

y= estimateCommonDisp(y,robust=TRUE)

#########Plot MDS
mds=plotMDS(y)

toplot <- data.frame(Dim1 = mds$x, Dim2 = mds$y)
type=rownames(toplot)


#########Plot MD
forMD <- exactTest(y,pair =list("RT","CE")  )
forMD$table["Padj"] =p.adjust(forMD$table$PValue, method = "BH")

plotMD(forMD, values=c(1,-1), col=c("red","blue"), 
       legend="topright")

library("xlsx")
write.xlsx(forMD$table, file = "forMD_BM_second.xls")    

######## for Volcano plot
# object construction
dds <- DESeqDataSetFromMatrix(res, DataFrame(cond), ~ cond)

# standard analysis
dds <- DESeq(dds)
count<-counts(dds,normalize=TRUE)

dds <- estimateSizeFactors(dds)
idx <- rowSums( counts(dds, normalized=TRUE) >= 5 ) >= 4

dds <- dds[idx,]
dds <- DESeq(dds)



for_volcano <- results(dds, contrast=c("cond","CE","RT"))

EnhancedVolcano(for_volcano,
                lab = rownames(for_volcano),
                title = 'CE vs RT',
                x = 'log2FoldChange',
                y = 'pvalue',
                xlim = c(-3, 3),
                ylim = c(0, 10),
                pCutoff = 10e-2,
                FCcutoff = 1,
                labSize= 5.0,
                pointSize = 3.0,
                colAlpha = 1)
##### PCA and MDS ########################################################################################################################################

######################## Retrieve the normalized data from SQL MouseDB database

data_norm <- sqlQuery(con,  "   select * from omics_results_norm where omics_name='RNA-seq_Martina_all_tissue' and sample_name like  '%Bonemarrow%' ")


#reshape the molt data to a dataframe
data_norm=dcast (data_norm, gene_id_name ~ sample_name, value.var='norm_value' )


res=as.matrix(data_norm[,-1])

rownames(res) <- data_norm$gene_id_name

#filter normalized counts
res[is.na(res)] <- 0

keep <- filterByExpr(res,min.count = 0.4, min.total.count = 1)
res <- res[keep,]


pca.normcount <- prcomp(t(res), center=TRUE, scale=TRUE)
pca_data<- data.frame(type = rownames(pca.normcount$x), PC1 = pca.normcount$x[,1], PC2=pca.normcount$x[,2])
summary(pca.normcount)


percentage <- round((pca.normcount$sdev^2) / (sum(pca.normcount$sdev^2))*100)
percentage <- paste( colnames(pca_data)[2:3], "(", paste( as.character(percentage), "%", ")", sep="") )

theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"),legend.position = "none")

# Plotting
p1=ggplot(toplot, aes(Dim1, Dim2,col=type) ) + geom_point(size=4) + geom_text_repel(aes(label = type), hjust=0, vjust=0,size=6,parse=TRUE)+theme
p2=ggplot(pca_data,aes(x=PC1, y=PC2, col=type)) + geom_point(size=4) + geom_text_repel(aes(label = type), hjust=0, vjust=0,size=6,parse=TRUE)+theme+  xlab(percentage[1]) + ylab(percentage[2])

grid.arrange(p1 + ggtitle("MDS") +   xlab("Leading logFC dim 1") + ylab("Leading logFC dim 2"), p2+ggtitle("PCA"), ncol=2)



###### Plotting the scatter of logFC############################################################################################################################

data_scatter <- sqlQuery(con, "with Table1 as (select * from omics_results_Fc
where omics_name = 'RNA-seq_Martina_SC' and comparison = 'D20RT VS D0RT' and logCPM > 0 and logFC  NOT BETWEEN -1 AND 1),
Table2 as (select * from omics_results_Fc
where omics_name = 'RNA-seq_Martina_SC' and comparison = 'D20CE VS D0CE' )
select Table1.gene_id_name,
Table1.logFC,
Table2.logFC
from Table1
join Table2
on Table1.gene_id_name = Table2.gene_id_name")




x.tick.number <- nrow(data_scatter)
at <- seq(1, nrow(data_scatter), length.out=x.tick.number)
y.tick.number <- 10
lablist<-data_scatter$gene_id_name
names(data_scatter)[names(data_scatter) == "gene_id_name"] <- "x"
names(data_scatter)[names(data_scatter) == "logFC"] <- "SC_EAED20_D0_RT"
names(data_scatter)[names(data_scatter) == "logFC.1"] <- "SC_EAED20_D0_COLD"

xyplot(SC_EAED20_D0_RT+SC_EAED20_D0_COLD~reorder(x, -SC_EAED20_D0_RT),data_scatter , auto.key = list(columns = 2,cex=1,title=""),main = "",
       xlab =list(label=  "gene index",cex=1.9),scales = list(x= list(tick.number=20),y=list(tick.number=10)),
       ylab = list(label=  "log2FC",cex=1.5), par.settings = list(superpose.symbol = list(pch = 16 ,cex = 1,alpha=0.5,
                                                                                          col = c("red", "blue"))),grid = TRUE,xlim=c(1, nrow(data_scatter)))

write.xlsx(data_scatter, file = "MC_EAE_RT_Cold.xlsx")    

################Pathwya Analysis
DE <- enrichPathway(unique(bitr(data_scatter$x,fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")$ENTREZID),organism = "mouse", pvalueCutoff = 0.05, qvalueCutoff = 1)
write.xlsx(DE, file = "Reactome_data_toPlot_MC_EAE_RT_Cold.xlsx",row.names=FALSE)   

rm(list = ls(all.names = TRUE)) # clear workspace
dev.off() # clear all plots

con = odbcDriverConnect('driver={SQL Server};server=localhost\\SQLEXPRESS;database=MouseDB;trusted_connection=true')

inputDir="C:/Users/hadadi/Dropbox/UniGen/database/RNAseq_noonoo/DB/with_rep/"
FileNames <- list.files(path=inputDir)

setwd("C:/Users/hadadi/Dropbox/UniGen/Martina/")


comparisons= list(list("Spleen_RT_rep","Spleen_CE_rep"))
cond= factor(c("BR", "BR", "BR", "BC", "BC", "BC" ,"BC"))

f = 'RNA-seq_Martina_MC_B.xlsx'
omics_name = 'new_FC_allTissue_3removed_SC_BM'
data = read_excel(paste(inputDir,f,sep=""))
raw_counts=as.matrix(data[,-1])
rownames(raw_counts) <- data$gene_id_name   
keep <- filterByExpr(raw_counts,min.count = 5, min.total.count = 4)
cts <- raw_counts[keep,]

y <- DGEList(cts, group= cond,  genes= row.names(cts) )
y <- calcNormFactors(y)
options(digits=3)
y= estimateCommonDisp(y,robust=TRUE)


CPM <- as.data.frame(cpm(y, normalized.lib.sizes = T, log = T))
write.xlsx(CPM, file = "RNA-seq_Martina_BM_secondAnalysis_normalized.xlsx",row.names=TRUE)    



for ( i in 1:Ncomparisons) 
{
  
  comparison0 = comparisons[[1]]
  comparison = paste0(comparison0[2]," VS ",comparison0[1])
  
  res <- exactTest(y ,  pair =  comparison0)$table
  res["Padj"] =p.adjust(res$PValue, method = "BH")
  gene_id_name = rownames(res)
  rownames(res) <- NULL
  res = cbind(omics_name,comparison,gene_id_name,res)
  sqlSave(con, res, tablename = "new_FC_allTissue_3removed_SC_BM",rownames=FALSE, append = TRUE)
  
  write.xlsx(res[,-c(1:2)], paste(omics_name,comparison,"_FC.xlsx"), row.names = FALSE)
  
}



mds=plotMDS(y)

toplot <- data.frame(Dim1 = mds$x, Dim2 = mds$y)
write.xlsx(toplot, file = "MDS_data_toPlot.xlsx",row.names=TRUE)    


CPM <- as.data.frame(cpm(y, normalized.lib.sizes = T, log = T))
pca.normcount <- prcomp(t(CPM), center=TRUE, scale=TRUE)
pca_data<- data.frame(type = rownames(pca.normcount$x), PC1 = pca.normcount$x[,1], PC2=pca.normcount$x[,2])

write.xlsx(pca_data, file = "PCA_data_toPlot.xlsx",row.names=FALSE)    


########Volcano
dds <- DESeqDataSetFromMatrix(cts, DataFrame(cond), ~ cond)


dds <- DESeq(dds)
count<-counts(dds,normalize=TRUE)

dds <- estimateSizeFactors(dds)

for_volcano <- results(dds, contrast=c("cond","BC","BR"))
write.xlsx(for_volcano, file = "FC_VolcanoPlot_MC_B.xlsx")    


EnhancedVolcano(for_volcano,
                lab = rownames(for_volcano),
                x = 'log2FoldChange',
                y = 'pvalue',
                ylim = c(0, 35),
                pCutoff = 10e-2,
                FCcutoff = 1,
                labSize= 5.0,
                pointSize = 3.0,
                colAlpha = 1)

ggsave("VolcanoPlot_SC_D20CE_RT.eps")
ggsave("VolcanoPlot_SC_D20CE_RT.pdf")


####Reactome#############################
library(ReactomePA)
library(org.Mm.eg.db)
library(clusterProfiler)
library(xlsx)

res <- exactTest(y,   list("D0RT","D0CE"))
res$table["Padj"] =p.adjust(res$table$PValue, method = "BH")


Dereg2= as.data.frame(res$table)

DE <- enrichPathway(unique(bitr(rownames(Dereg[ Dereg$PValue <= 0.05 , ]),fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")$ENTREZID),organism = "mouse", pvalueCutoff = 0.05, qvalueCutoff = 1)
write.xlsx(DE, file = "Reactome_data_toPlot_D0CE_D0RT.xlsx",row.names=FALSE)    

###### TopGO analysis ############################################################################################################################

library(topGO)
library(Rgraphviz)

Dereg$de <- Dereg$PValue < 0.05
de_results=Dereg
gene_universe <- as.numeric(de_results$de)
gene_universe <- factor(gene_universe)

names(gene_universe) <- rownames(de_results)

myGenes=de_results[de_results$de == TRUE,]
myGenes=rownames(myGenes) 



go_data_BP <- new("topGOdata",
                  ontology = "BP",
                  allGenes = gene_universe,
                  nodeSize = 5,
                  annotationFun = annFUN.org,
                  mapping = "org.Mm.eg",
                  ID = "SYMBOL")

go_data_MF <- new("topGOdata",
                  ontology = "MF",
                  allGenes = gene_universe,
                  nodeSize = 5,
                  annotationFun = annFUN.org,
                  mapping = "org.Mm.eg",
                  ID = "SYMBOL")

go_data_CC <- new("topGOdata",
                  ontology = "CC",
                  allGenes = gene_universe,
                  nodeSize = 5,
                  annotationFun = annFUN.org,
                  mapping = "org.Mm.eg",
                  ID = "SYMBOL")

go_test_BP <- runTest(go_data_BP, algorithm = "weight01", statistic = "fisher")
go_table_BP <- GenTable(go_data_BP, weightFisher = go_test_BP,
                        orderBy = "weightFisher", ranksOf = "weightFisher",
                        topNodes = sum(score(go_test_BP) < .01),numChar = 500)





go_table_BP$genes <- sapply(go_table_BP$GO.ID, function(x)
{
  genes<-genesInTerm(go_data_BP, x)
  genes[[1]][genes[[1]] %in% myGenes] # myGenes is the queried gene list
})
go_table_BP$genes[which(go_table_BP$weightFisher<0.05)] # print those only with p-value < 0.05

go_table_BP$genes <-vapply(go_table_BP$genes, paste, collapse = ",", character(1L))



write.xlsx(go_table_BP,file ="go_table_BP.xlsx",row.names=FALSE)

printGraph(go_data_BP, go_test_BP, firstSigNodes = 5, fn.prefix = "GO_BP", useInfo = "all", pdfSW = TRUE)



go_test_MF <- runTest(go_data_MF, algorithm = "weight01", statistic = "fisher")
go_table_MF <- GenTable(go_data_MF, weightFisher = go_test_MF,
                        orderBy = "weightFisher", ranksOf = "weightFisher",
                        topNodes = sum(score(go_test_MF) < .01),numChar = 500)

go_table_MF$genes <- sapply(go_table_MF$GO.ID, function(x)
{
  genes<-genesInTerm(go_data_MF, x)
  genes[[1]][genes[[1]] %in% myGenes] # myGenes is the queried gene list
})
go_table_MF$genes[which(go_table_MF$weightFisher<0.05)] # print those only with p-value < 0.05

go_table_MF$genes <-vapply(go_table_MF$genes, paste, collapse = ",", character(1L))



write.xlsx(go_table_MF,file ="go_table_MF.xlsx",row.names=FALSE)

printGraph(go_data_MF, go_test_MF, firstSigNodes = 5, fn.prefix = "GO_MF", useInfo = "all", pdfSW = TRUE)



go_test_CC <- runTest(go_data_CC, algorithm = "weight01", statistic = "fisher")
go_table_CC <- GenTable(go_data_CC, weightFisher = go_test_CC,
                        orderBy = "weightFisher", ranksOf = "weightFisher",
                        topNodes = sum(score(go_test_CC) < .01),numChar = 500)

go_table_CC$genes <- sapply(go_table_CC$GO.ID, function(x)
{
  genes<-genesInTerm(go_data_CC, x)
  genes[[1]][genes[[1]] %in% myGenes] # myGenes is the queried gene list
})
go_table_CC$genes[which(go_table_CC$weightFisher<0.05)] # print those only with p-value < 0.05

go_table_CC$genes <-vapply(go_table_CC$genes, paste, collapse = ",", character(1L))


write.xlsx(go_table_CC,file ="go_table_CC.xlsx",row.names=FALSE)

printGraph(go_data_CC, go_test_CC, firstSigNodes = 5, fn.prefix = "GO_CC", useInfo = "all", pdfSW = TRUE)




