## This script will take data from callpeaks files in order to  
## determine expression peaks and the genes involved.
## 
## 
## 
## 
## 
## 

params<- commandArgs(trailingOnly=T)
print("First argument:")
print(as.numeric(params[[1]]))

print("Second argument:")
print(params[[2]])

print("Second argument:")
print(params[[3]])

##Parameters are read
promoter_length<-as.numeric(params[[1]])
peaks <- params[[2]]
peaks_summits <- params[[3]]

## Loading packages...
library(ChIPseeker)
library(TxDb.Athaliana.BioMart.plantsmart28)
library(org.At.tair.db)
library(DO.db)
library(clusterProfiler)
library(pathview)


txdb <- TxDb.Athaliana.BioMart.plantsmart28


jpeg("covplot_peaks.jpg")
covplot(peaks, weightCol = "V5")
dev.off()

## Promotor is defined
promoter <- getPromoters(TxDb=txdb, upstream=promoter_length, downstream=promoter_length)

##Peaks are annoted using annotatePeak function

peak_annotation <- annotatePeak(peak = peaks, tssRegion = c(-promoter_length, promoter_length), TxDb = txdb, annoDb = org.At.tair)

jpeg("annotation_pie.jpg")
plotAnnoPie(peak_annotation)
dev.off()

jpeg("annotation_bar.jpg")
plotAnnoBar(peak_annotation)
dev.off()

##Target genes are obtained, considering them as genes which promotor is recognised by the transcription factor
peak_annotation_dataframe <- as.data.frame(peak_annotation)

target_genes <- peak_annotation_dataframe$geneId[peak_annotation_dataframe$annotation == "Promoter"]
target_genes<-unique(target_genes)
write(x = target_genes,file = "target_genes.txt")

#GO terms enrichment
peak_annotation_summits <- annotatePeak(peak = peaks_summits, tssRegion=c(-promoter_length, promoter_length), TxDb=txdb)
peak_annotation_summits_dataframe <- as.data.frame(peak_annotation_summits)
target_genes_summits <- subset(peak_annotation_summits_dataframe, annotation == "Promoter")
TF_regulome<-target_genes_summits$geneId

txdb_genes<-as.data.frame(genes(txdb))
atha_genes<-txdb_genes$gene_id #universe


eGO_BP <- enrichGO(gene = TF_regulome, OrgDb = org.At.tair.db, 
                universe = atha_genes, keyType = "TAIR", ont = "BP")

eGO_MF <- enrichGO(gene = TF_regulome, OrgDb = org.At.tair.db, 
                   universe = atha_genes, keyType = "TAIR", ont = "MF")

eGO_CC <- enrichGO(gene = TF_regulome, OrgDb = org.At.tair.db, 
                   universe = atha_genes, keyType = "TAIR", ont = "CC")


#Plots
jpeg("GO_BiologicalProcess_Barplot.jpg")
barplot(eGO_BP, showCategory = 20)
dev.off()

jpeg("GO_MolecularFunction_Barplot.jpg")
barplot(eGO_MF, showCategory = 20)
dev.off()

jpeg("GO_CelularComponent_Barplot.jpg")
barplot(eGO_CC, showCategory = 20)
dev.off()


#KEGG Pathways.

#Photosintesis transduction pathway
keggpathway<-pathview(gene.data=target_genes, species = "ath", pathway.id = "ath00195", gene.idtype = "TAIR")

#DNA replication pathway
keggpathway<-pathview(gene.data=target_genes, species = "ath", pathway.id = "ath03030", gene.idtype = "TAIR")