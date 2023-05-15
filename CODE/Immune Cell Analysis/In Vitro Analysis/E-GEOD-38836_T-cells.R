# Author: Marissa McDonald
# 440 Project Code
# Early activation T-cell analysis, E-GEOD-38836

# setting my directory
setwd("~/MEMP/440/Project/Tcells/E-GEOD-38836")

#loading data
data =  load("GLDS-13_array_normalized-annotated.rda")

#annotating exp dat
library(biomaRt)
mart<- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

#selecting for coding mrna, removing non coding rna
refSeq <- rownames(eset)
refSeq = refSeq[grepl("NM",refSeq)]
esetNew <- eset[rownames(eset) %in% refSeq,]

#getting gene info 
geneInfo <- getBM(filters="refseq_mrna", attributes=c("refseq_mrna","external_gene_name"), values=refSeq, mart=mart)
#remove duplicates in geneInfo, this is when genes correspond to the same ref seq id? 
geneInfo <- geneInfo[!duplicated(geneInfo$external_gene_name),]
geneInfo <- geneInfo[!duplicated(geneInfo$refseq_mrna),]

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

#adding gene symbol column
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
esetExp <- esetNew[,1:6]
groupTag <- data.frame('1G-Activated','1G-Activated','1G-Activated','uG-Activated','uG-Activated','uG-Activated')
colnames(groupTag) <- colnames(esetExp)
esetImmu <- rbind(groupTag, esetExp)

Subsets <- ImmuCellAI_new(esetImmu,"microarray",1,0,0)
subsetAbundance <- data.frame(Subsets[["Sample_abundance"]])
groupResult <- data.frame(Subsets$Group_result)
#No significant cell type diffs

#DE analysis
#testing to see if there are differences
#BiocManager::install("limma")
library(limma)
groups <- c('1GActivated','1GActivated','1GActivated','uGActivated','uGActivated','uGActivated','uGNonActivated','uGNonActivated','uGNonActivated')
design <- model.matrix(~0 + groups)
fit<- lmFit(esetNew, design)
contruG21G <- makeContrasts(groupsuGActivated - groups1GActivated, levels = colnames(coef(fit)))
tmpuG21G <- contrasts.fit(fit,contruG21G)
tmpuG21G <- eBayes(tmpuG21G)
top.tableUG21G <- topTable(tmpuG21G,sort.by="P",n=Inf,p.value=0.05)
#save table
library(openxlsx)
write.xlsx(top.tableUG21G,file='DEGenes.xlsx',colNames=TRUE, rowNames=TRUE)

#showing as heatmap
#install.packages("pheatmap")
library(pheatmap)
library(ggplot2)

#filtering exp data to include only significant DE genes
#rownames(esetNew) <- str_to_title(rownames(esetNew))
esetNewfilt <- esetNew[rownames(esetNew) %in% rownames(top.tableUG21G),]
esetNewfilt$Activated1G <- rowMeans(esetNewfilt[,1:3])
esetNewfilt$ActivateduG <- rowMeans(esetNewfilt[,4:6])
esetNewfilt$NonActivateduG <- rowMeans(esetNewfilt[,7:9])
esetNewfilt[,1:9] <- NULL

matDataNew <- as.matrix(esetNewfilt)
pheatmap(matDataNew, fontsize_row=8, angle_col=0)

#enrichment for chronic inflammation, CVD
library(enrichR)
#finding upregulated and downregulated genes and processes in spaceflight
upSpaceFlight <- rownames(top.tableUG21G[top.tableUG21G$logFC > 0,])
downSpaceFlight <- rownames(top.tableUG21G[top.tableUG21G$logFC < 0,])

#No upreg. genes, only downreg
enr1Gact <- enrichr(downSpaceFlight,c("DisGeNET","GO_Biological_Process_2021"))

enr1Gact[[1]]$P.value <- enr1Gact[[1]]$Adjusted.P.value
enr1Gact[[1]] <- enr1Gact$DisGeNET[enr1Gact$DisGeNET$Adjusted.P.value <0.05,]
write.xlsx(enr1Gact[[1]],file="Diseases.xlsx",rowNames=TRUE,colNames=TRUE)

filter <- c("Vascular Diseases", "Myocardial Ischemia", "Arteriosclerosis", "Atherosclerosis", "Vascular inflammations", "Inflammation",
            "Hypertrophic Cardiomyopathy","Vascular lesions", "Coronary Arteriosclerosis","Disseminated Intravascular Coagulation",
            "Cardiovascular Diseases", "Stunned Myocardium","Cardiac fibrosis")
enr1GactDiseaseFilt <- enr1Gact$DisGeNET[enr1Gact$DisGeNET$Term %in% filter,]
plotEnrich(enr1GactDiseaseFilt,showTerms = 20,numChar = 40,
           y = "Count", orderBy = "P.value", xlab = "Disease Profile", ylab = NULL, 
           title = "Enriched Diseases of Activated Human T Cells in 1G")

enr1Gact[[2]]$P.value <- enr1Gact[[2]]$Adjusted.P.value
enr1Gact[[2]] <- enr1Gact$GO_Biological_Process_2021[enr1Gact$GO_Biological_Process_2021$Adjusted.P.value <0.05,]
write.xlsx(enr1Gact[[2]],file="processes.xlsx",rowNames=TRUE,colNames=TRUE)

for (i in 1:length(enr1Gact[[2]]$Term)) {
  s <- enr1Gact[[2]]$Term[i]
  enr1Gact[[2]]$Term[i] <- unlist(strsplit(s, split='(', fixed=TRUE))[1]
}

filter <- c("cytokine-mediated signaling pathway ",
            "cellular response to cytokine stimulus ",
            "fat cell differentiation ",
            "cellular response to tumor necrosis factor ",
            "cellular response to vascular endothelial growth factor stimulus ",
            "cellular response to granulocyte macrophage colony-stimulating factor stimulus ",
            "endothelial cell chemotaxis ",
            "mast cell mediated immunity ",
            "myeloid dendritic cell activation ",
            "regulation of monocyte differentiation ",
            "regulation of fat cell differentiation ",
            "negative regulation of myeloid cell differentiation ",
            "blood vessel endothelial cell migration ",
            "tumor necrosis factor-mediated signaling pathway ",
            "negative regulation of interleukin-6 production ",
            "positive regulation of leukocyte cell-cell adhesion ",
            "positive regulation of myeloid leukocyte differentiation ",
            "endothelial cell migration ",
            "regulation of chemokine production ",
            "positive regulation of smooth muscle cell proliferation ",
            "regulation of interleukin-2 production ",
            "regulation of smooth muscle cell proliferation ",
            "positive regulation of cytokine production ",
            "positive regulation of fat cell differentiation ",
            "positive regulation of chemokine production ",
            "positive regulation of NIK/NF-kappaB signaling ",
            "positive regulation of interleukin-1 beta production ",
            "negative regulation of cytokine production ",
            "positive regulation of interleukin-1 production ",
            "cellular response to type I interferon ",
            "type I interferon signaling pathway ",
            "regulation of myeloid cell differentiation ")

enr1GactProcessesFiltered <- enr1Gact$GO_Biological_Process_2021[enr1Gact$GO_Biological_Process_2021$Term %in% filter,]
                                                                                                                                    
plotEnrich(enr1GactProcessesFiltered,showTerms = 32,numChar = 100,
           y = "Count", orderBy = "P.value", xlab = "Biological Process", ylab = NULL, 
           title = "Enriched Pathways of Activated Human T Cells in 1G")
