# Author: Marissa McDonald
# 440 Project Code
# Tcell analysis, E-GEOD-4209

# setting my directory
setwd("~/MEMP/440/Project/Tcells/E-GEOD-4209")

#update R if you need to, version used to execute code is 4.2.2
#install.packages("installr")
#library(installr)
#updateR()

#loading data
#install.packages("readxl")
library("readxl")
expData <- read_excel("E-GEOD-4209-processed-data-1519310732 4.xlsx")

#converting affymatrix to gene symbols with hgu133a.db from bioconductor
#if (!require("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")

#BiocManager::install("hgu133a.db")
library(hgu133a.db)

#getting probe IDs and annotations
probeIDs <- expData$`Composite Element REF`
annotLookup <- select(hgu133a.db, keys = probeIDs,
                      columns = c('PROBEID', 'ENSEMBL', 'SYMBOL'))

#there are duplicates in the annotation key so remove these
annotLookupNew <- annotLookup[!duplicated(annotLookup$PROBEID), ]

#empty df for storage
geneSym = data.frame(matrix(nrow = 0, ncol = 1)) 
#adding to df based on corresponding gene symbol
for (i in probeIDs) {
  if (i %in% annotLookupNew$PROBEID) {
    geneSym <- rbind(geneSym, annotLookupNew$SYMBOL[annotLookupNew$PROBEID==i])
  } else {
    geneSym <- rbind(geneSym, "NA")
  }
}
colnames(geneSym) = "Symbol"

#rename expData with gene Symbol
expDataNew <- expData
expDataNew$`Composite Element REF`<- geneSym$Symbol
colnames(expDataNew)[1] <- 'Gene Symbol'

#delete rows with Gene Symbol = NA, these are likely controls
completeData <- expDataNew[complete.cases(expDataNew),]

#for duplicated rows, keep row with greatest average difference between groups
n_occur <- data.frame(table(completeData$`Gene Symbol`))
repeats <- completeData[completeData$`Gene Symbol` %in% n_occur$Var1[n_occur$Freq > 1],]
repeats$AverageFlight <- rowMeans(repeats[,c(2,5,6)])
repeats$AverageControl <- rowMeans(repeats[,c(3,4,7)])
repeats$Difference <- abs(repeats$AverageFlight - repeats$AverageControl)

for (gene in n_occur$Var1[n_occur$Freq > 1]){
  #print(gene)
  subset <- repeats[repeats$`Gene Symbol`== gene,]
  repeats[repeats$`Gene Symbol`== gene,] <- subset[which.max(subset$Difference),]
}

#remove duplicates from data & removing unnecessary columns
completeDataNew <- repeats[!duplicated(repeats$`Gene Symbol`),]
completeDataNew <- as.data.frame(completeDataNew)
rownames(completeDataNew) <- completeDataNew$`Gene Symbol`
completeDataNew$`Gene Symbol` <- NULL
completeDataNew[,7:9] <- NULL

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

#adding group tag row, needs to be first row
groupTag <- data.frame('MMG','static','static','MMG','MMG','static')
colnames(groupTag) <- colnames(completeDataNew)
completeDataImmu <- rbind(groupTag, completeDataNew)

Subsets <- ImmuCellAI_new(completeDataImmu,"microarray",1,0,0)
subsetAbundance <- data.frame(Subsets[["Sample_abundance"]])
groupResult <- data.frame(Subsets$Group_result)
#No significant differences 

#DE analysis
#install.packages("pheatmap")
library(pheatmap)
library(ggplot2)

#preliminary heatmap
#BiocManager::install("limma")
library(limma)
groups <- c('MMG','static','static','MMG','MMG','static')
design <- model.matrix(~0 + groups)
fit<- lmFit(completeDataNew, design)
contr <- makeContrasts(groupsMMG - groupsstatic, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit,contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp,sort.by="P",n=Inf)
view <- top.table[1:5,]

#No significantly DE genes--> stopping analysis here for the paper  


# Optional Continued Analysis ---------------------------------------------

#enrichment for chronic inflammation, CVD

library(enrichR)


#finding upregulated genes in spaceflight with DisGeNET db
upSpaceFlight <- rownames(view[view$logFC > 0,])

#finding library of gene list for enrichR
dbs <- listEnrichrDbs()
dbs <- c("DisGeNET") #this is a gene set related to human disease

enrSpaceFlight <- enrichr(upSpaceFlight, dbs[1])

enrSpaceFlight[[1]]$P.value <- enrSpaceFlight[[1]]$Adjusted.P.value
#plotEnrich(enrSpaceFlight[["DisGeNET"]],showTerms = 20,numChar = 40,
           #y = "Count", orderBy = "P.value", xlab = "Enriched Disease Profile", ylab = NULL, 
           #title = "Enrichment Analysis of T-cells in Spaceflight")
enrSpaceFlight$Filtered <- enrSpaceFlight$DisGeNET[enrSpaceFlight$DisGeNET$Term %in% c('Carotid artery occlusion', 'Carotid Stenosis', 'Myocardial Ischemia', 'Amyloidosis', 'Cardiovascular Diseases','Vascular Diseases' ),]
plotEnrich(enrSpaceFlight[["Filtered"]],showTerms = 20,numChar = 40,
           y = "Count", orderBy = "P.value", xlab = "Disease Profile", ylab = NULL, 
           title = "Enrichment Analysis of T-cells in Simulated Microgravity")

enrSpaceFlightBio <- enrichr(upSpaceFlight, "GO_Biological_Process_2021")

enrSpaceFlightBio[[1]]$P.value <- enrSpaceFlightBio[[1]]$Adjusted.P.value
for (i in 1:length(enrSpaceFlightBio[[1]]$Term)) {
  s <- enrSpaceFlightBio[[1]]$Term[i]
  enrSpaceFlightBio[[1]]$Term[i] <- unlist(strsplit(s, split='(', fixed=TRUE))[1]
}

enrSpaceFlightBio$Filtered <- enrSpaceFlightBio$GO_Biological_Process_2021[enrSpaceFlightBio$GO_Biological_Process_2021$Term %in% c('positive regulation of myeloid cell differentiation ',
                                                                                                                                    'positive regulation of interleukin-8 production ','negative regulation of transforming growth factor beta receptor signaling pathway ','positive regulation of NF-kappaB transcription factor activity ',
                                                                                                                                    'positive regulation of cytokine production ','neutrophil degranulation ','neutrophil activation involved in immune response ','neutrophil mediated immunity '),]
plotEnrich(enrSpaceFlightBio[[2]],showTerms = 20,numChar = 60,
           y = "Count", orderBy = "P.value", xlab = "Biological Process", ylab = NULL, 
           title = "Enrichment Analysis of T-cells in Simulated Microgravity")
processes <- as.data.frame(enrSpaceFlightBio$GO_Biological_Process_2021$Term)
processes$Ad.PValue <- enrSpaceFlightBio$GO_Biological_Process_2021$Adjusted.P.value
processes <- processes[processes$Ad.PValue < 0.05,]


