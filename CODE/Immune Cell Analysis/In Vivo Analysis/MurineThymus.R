# Author: Marissa McDonald
# 440 Project Code
# Thymus analysis, GSE18388

# setting my directory
setwd("~/MEMP/440/Project/MurineThymus")

#update R if you need to, version used to execute code is 4.2.2
#install.packages("installr")
#library(installr)
#updateR()

#loading data

data =  load("GLDS-4_array_normalized-annotated.rda")

#converting RefSeq ids to gene symbols

#BiocManager::install("biomaRt")
library(biomaRt)
mart<- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))

#selecting for coding mrna, removing nc rna
refSeq <- rownames(eset)
refSeq = refSeq[grepl("NM",refSeq)]
esetNew <- eset[rownames(eset) %in% refSeq,]
geneInfo <- getBM(filters="refseq_mrna", attributes=c("refseq_mrna","external_gene_name"), values=refSeq, mart=mart)
#remove duplicates in geneInfo, this is when genes correspond to the same ref seq id? 
geneInfo <- geneInfo[!duplicated(geneInfo$external_gene_name),]
#now renaming rownames of dataset
newRownames <- list()
for (i in rownames(esetNew)){
  if (i %in% geneInfo$refseq_mrna){
    geneSym <- geneInfo[geneInfo$refseq_mrna == i,2]
  } else {
    geneSym <- "NA"
  }
  #print(geneSym)
  newRownames <- append(newRownames,geneSym)
}

#adding column
esetNew$GeneSymbol <- t(as.data.frame(newRownames))


#delete rows with Gene Symbol = NA, these are likely controls
esetNew <- esetNew[!(esetNew$GeneSymbol=="NA"),]


#remove duplicates from data & removing unnecessary columns
rownames(esetNew) <- esetNew$GeneSymbol
esetNew$GeneSymbol <- NULL


#running ImmuCellAI for t cell subsets
# microarray data must be log2-transformed signal
#install.packages("devtools")
# library(devtools)
# BiocManager::install("GSVA")
# install_github("lydiaMyr/ImmuCellAI@main")
library(ImmuCellAI)
#The prediction process of ImmuCellAI is updated, which simulated the flow cytometry process to predict cell type abundance by hierarchical strategy. All 24 cell types were divided into two layers, layer1:DC, B cell, Monocyte, Macrophage, NK, Neutrophil, CD4 T, CD8 T, NKT, Tgd; layer2:CD4 naive, CD8 naive, Tc, Tex, Tr1, nTreg, iTreg, Th1, Th2, Th17, Tfh, Tcm, Tem, MAIT.
#Parameters
#sample_expression: Sample expression profile in FPKM, TPM format by RNA-seq or log2-transformed signal by microarray.
#datatype: One of "rnaseq" and "microarray"
#group_tag: One of 0 and 1, if there is the need to perform the comparision between different groups. If the value is 1, users need to add a group tag row in the input epxression matrix to explain the group of each sample.
#response_tag: One of 0 and 1, if there is the need to predict the ICB response of each sample.
#rownames aka gene symbols must be in all caps!!


#adding group tag row, needs to be first row
groupTag <- data.frame('Flight','Flight','Flight','Flight','Ground Control','Ground Control','Ground Control','Ground Control')
colnames(groupTag) <- colnames(esetNew)
rownames(groupTag) <- "Group Tag"
esetImmu <- rbind(groupTag, esetNew)
rownames(esetImmu) <- toupper(rownames(esetImmu))

Subsets <- ImmuCellAI_new(esetImmu,"microarray",1,0,0)
subsetAbundance <- data.frame(Subsets[["Sample_abundance"]])
groupResult <- data.frame(Subsets$Group_result)
#flight group had less central memory t cells, only sig diff
library(openxlsx)
write.xlsx(groupResult, file='cellTypeDiffs.xlsx',rowNames=TRUE,colnames=TRUE)

#DE analysis
#testing to see if there are differences
#install.packages("pheatmap")
library(pheatmap)
library(ggplot2)

#BiocManager::install("limma")
library(limma)
groups <- c('Flight','Flight','Flight','Flight','GroundControl','GroundControl','GroundControl','GroundControl')
design <- model.matrix(~0 + groups)
fit<- lmFit(esetNew, design)
contr <- makeContrasts(groupsFlight - groupsGroundControl, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit,contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp,sort.by="P",n=Inf,p.value=0.05)
write.xlsx(top.table, file='DEgenes.xlsx',rowNames=TRUE,colnames=TRUE)


#filtering for heatmap
esetNewfilt <- esetNew[rownames(esetNew) %in% rownames(top.table),]
esetNewfilt$Flight <- rowMeans(esetNewfilt[,1:4])
esetNewfilt$GroundControl <- rowMeans(esetNewfilt[,5:8])
esetNewfilt[,1:8] <- NULL

matDataNew <- as.matrix(esetNewfilt)
pheatmap(matDataNew, fontsize_row=8, angle_col=0)

#enrichment for chronic inflammation, CVD
library(enrichR)
#finding upregulated and downregulated genes in spaceflight 
upSpaceFlight <- rownames(top.table[top.table$logFC > 0,])
downSpaceFlight <- rownames(top.table[top.table$logFC < 0,])

upSpaceFlightHum <- toupper(upSpaceFlight)
enrSpaceFlightDisease <- enrichr(upSpaceFlightHum, "DisGeNET")

enrSpaceFlightDisease[[1]]$P.value <- enrSpaceFlightDisease[[1]]$Adjusted.P.value
enrSpaceFlightDisease[[1]] <- enrSpaceFlightDisease$DisGeNET[enrSpaceFlightDisease$DisGeNET$Adjusted.P.value <0.05,]
plotEnrich(enrSpaceFlightDisease[[1]],showTerms = 47,numChar = 100,
           y = "Count", orderBy = "P.value", xlab = "DiseaseProfile", ylab = NULL, 
           title = "Enriched Diseases of Upregulated Genes")

#one: neurogenic inflamm

downSpaceFlightHum <- toupper(downSpaceFlight)
enrSpaceFlightDiseaseDown <- enrichr(downSpaceFlightHum, "DisGeNET")

enrSpaceFlightDiseaseDown[[1]]$P.value <- enrSpaceFlightDiseaseDown[[1]]$Adjusted.P.value
enrSpaceFlightDiseaseDown[[1]] <- enrSpaceFlightDiseaseDown$DisGeNET[enrSpaceFlightDiseaseDown$DisGeNET$Adjusted.P.value <0.05,]
plotEnrich(enrSpaceFlightDiseaseDown[[1]],showTerms = 47,numChar = 100,
           y = "Count", orderBy = "P.value", xlab = "DiseaseProfile", ylab = NULL, 
           title = "Enriched Diseases of Downregulated Genes")


enrGroundBio <- enrichr(downSpaceFlightHum, "GO_Biological_Process_2021")
enrGroundBio[[1]]$P.value <- enrGroundBio[[1]]$Adjusted.P.value

for (i in 1:length(enrGroundBio[[1]]$Term)) {
  s <- enrGroundBio[[1]]$Term[i]
  enrGroundBio[[1]]$Term[i] <- unlist(strsplit(s, split='(', fixed=TRUE))[1]
}


enrGroundBio$GO_Biological_Process_2021 <- enrGroundBio$GO_Biological_Process_2021[enrGroundBio$GO_Biological_Process_2021$Adjusted.P.value <0.05,]
write.xlsx(enrGroundBio[[1]], file='processesDown.xlsx',rowNames=TRUE,colnames=TRUE)

filter <- c("regulation of fever generation ",
            "regulation of endothelial tube morphogenesis ",
            "positive regulation of fever generation ",
            "regulation of NK T cell activation ",
            "positive regulation of NK T cell activation ",
            "T cell chemotaxis ",
            "phasic smooth muscle contraction ",
            "negative regulation of secretion ",
            "positive regulation of lymphocyte migration ",
            "positive regulation of acute inflammatory response ",
            "regulation of T cell chemotaxis ",
            "inflammatory response ")

enrGroundBioFilt <- enrGroundBio$GO_Biological_Process_2021[enrGroundBio$GO_Biological_Process_2021$Term %in% filter,]
plotEnrich(enrGroundBioFilt,showTerms = 47,numChar = 100,
           y = "Count", orderBy = "P.value", xlab = "Biological Process", ylab = NULL, 
           title = "Downregulated Pathways of Mouth Thymus in Space")

enrFlightBio <- enrichr(upSpaceFlightHum, "GO_Biological_Process_2021")
enrFlightBio$GO_Biological_Process_2021 <- enrFlightBio$GO_Biological_Process_2021[enrFlightBio$GO_Biological_Process_2021$Adjusted.P.value <0.05,]

enrFlightBio[[1]]$P.value <- enrFlightBio[[1]]$Adjusted.P.value
for (i in 1:length(enrFlightBio[[1]]$Term)) {
  s <- enrFlightBio[[1]]$Term[i]
  enrFlightBio[[1]]$Term[i] <- unlist(strsplit(s, split='(', fixed=TRUE))[1]
}
write.xlsx(enrFlightBio[[1]], file='processesUp.xlsx',rowNames=TRUE,colnames=TRUE)

plotEnrich(enrFlightBio[[1]],showTerms = 20,numChar = 100,
           y = "Count", orderBy = "P.value", xlab = "Biological Process", ylab = NULL, 
           title = "Upregulated Pathways of Mouse Thymus in Space")