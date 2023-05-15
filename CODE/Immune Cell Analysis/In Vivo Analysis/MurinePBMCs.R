# Author: Marissa McDonald
# 440 Project Code
# Mouse PBMC analysis

# setting my directory
setwd("~/MEMP/440/Project/PBMCs")

library(readr)
library(Matrix)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(purrr)
library(tibble)
library(openxlsx)
library(HGNChelper)
#library(metap)


#create individidual seurats, filtering out cells with < 200 genes and genes detected in <3 cells
setwd("~/MEMP/440/Project/PBMCs/GC_LAR_YNG_GY9")
matrixGC_YNG_GY9 = readMM('matrix.mtx')
rownames(matrixGC_YNG_GY9) = read_tsv("features.tsv",col_names=FALSE)[, 2, drop=TRUE]
colnames(matrixGC_YNG_GY9) = read_tsv("barcodes.tsv",col_names=FALSE)[, 1, drop=TRUE]
GC_YNG_GY9 <- CreateSeuratObject(counts=matrixGC_YNG_GY9,min.cells = 3, min.features = 200)

setwd("~/MEMP/440/Project/PBMCs/GC_LAR_YNG_GY10")
matrixGC_YNG_GY10 = readMM('matrix.mtx')
rownames(matrixGC_YNG_GY10) = read_tsv("features.tsv",col_names=FALSE)[, 2, drop=TRUE]
colnames(matrixGC_YNG_GY10) = read_tsv("barcodes.tsv",col_names=FALSE)[, 1, drop=TRUE]
GC_YNG_GY10 <- CreateSeuratObject(counts=matrixGC_YNG_GY10,min.cells = 3, min.features = 200)

setwd("~/MEMP/440/Project/PBMCs/GC_LAR_YNG_GY6")
matrixGC_YNG_GY6 = readMM('matrix.mtx')
rownames(matrixGC_YNG_GY6) = read_tsv("features.tsv",col_names=FALSE)[, 2, drop=TRUE]
colnames(matrixGC_YNG_GY6) = read_tsv("barcodes.tsv",col_names=FALSE)[, 1, drop=TRUE]
GC_YNG_GY6 <- CreateSeuratObject(counts=matrixGC_YNG_GY6,min.cells = 3, min.features = 200)

setwd("~/MEMP/440/Project/PBMCs/GC_LAR_YNG_GY4")
matrixGC_YNG_GY4 = readMM('matrix.mtx')
rownames(matrixGC_YNG_GY4) = read_tsv("features.tsv",col_names=FALSE)[, 2, drop=TRUE]
colnames(matrixGC_YNG_GY4) = read_tsv("barcodes.tsv",col_names=FALSE)[, 1, drop=TRUE]
GC_YNG_GY4 <- CreateSeuratObject(counts=matrixGC_YNG_GY4,min.cells = 3, min.features = 200)

setwd("~/MEMP/440/Project/PBMCs/GC_LAR_OLD_GO19")
matrixGC_OLD_GO19 = readMM('matrix.mtx')
rownames(matrixGC_OLD_GO19) = read_tsv("features.tsv",col_names=FALSE)[, 2, drop=TRUE]
colnames(matrixGC_OLD_GO19) = read_tsv("barcodes.tsv",col_names=FALSE)[, 1, drop=TRUE]
GC_OLD_GO19 <- CreateSeuratObject(counts=matrixGC_OLD_GO19,min.cells = 3, min.features = 200)

setwd("~/MEMP/440/Project/PBMCs/GC_LAR_OLD_GO20")
matrixGC_OLD_GO20 = readMM('matrix.mtx')
rownames(matrixGC_OLD_GO20) = read_tsv("features.tsv",col_names=FALSE)[, 2, drop=TRUE]
colnames(matrixGC_OLD_GO20) = read_tsv("barcodes.tsv",col_names=FALSE)[, 1, drop=TRUE]
GC_OLD_GO20 <- CreateSeuratObject(counts=matrixGC_OLD_GO20,min.cells = 3, min.features = 200)

setwd("~/MEMP/440/Project/PBMCs/GC_LAR_OLD_GO16")
matrixGC_OLD_GO16 = readMM('matrix.mtx')
rownames(matrixGC_OLD_GO16) = read_tsv("features.tsv",col_names=FALSE)[, 2, drop=TRUE]
colnames(matrixGC_OLD_GO16) = read_tsv("barcodes.tsv",col_names=FALSE)[, 1, drop=TRUE]
GC_OLD_GO16 <- CreateSeuratObject(counts=matrixGC_OLD_GO16,min.cells = 3, min.features = 200)

setwd("~/MEMP/440/Project/PBMCs/GC_LAR_OLD_GO14")
matrixGC_OLD_GO14 = readMM('matrix.mtx')
rownames(matrixGC_OLD_GO14) = read_tsv("features.tsv",col_names=FALSE)[, 2, drop=TRUE]
colnames(matrixGC_OLD_GO14) = read_tsv("barcodes.tsv",col_names=FALSE)[, 1, drop=TRUE]
GC_OLD_GO14 <- CreateSeuratObject(counts=matrixGC_OLD_GO14,min.cells = 3, min.features = 200)

setwd("~/MEMP/440/Project/PBMCs/FLT_LAR_YNG_FY15")
matrixFLT_LAR_YNG_FY15 = readMM('matrix.mtx')
rownames(matrixFLT_LAR_YNG_FY15) = read_tsv("features.tsv",col_names=FALSE)[, 2, drop=TRUE]
colnames(matrixFLT_LAR_YNG_FY15) = read_tsv("barcodes.tsv",col_names=FALSE)[, 1, drop=TRUE]
FLT_LAR_YNG_FY15 <- CreateSeuratObject(counts=matrixFLT_LAR_YNG_FY15,min.cells = 3, min.features = 200)

setwd("~/MEMP/440/Project/PBMCs/FLT_LAR_YNG_FY14")
matrixFLT_LAR_YNG_FY14 = readMM('matrix.mtx')
rownames(matrixFLT_LAR_YNG_FY14) = read_tsv("features.tsv",col_names=FALSE)[, 2, drop=TRUE]
colnames(matrixFLT_LAR_YNG_FY14) = read_tsv("barcodes.tsv",col_names=FALSE)[, 1, drop=TRUE]
FLT_LAR_YNG_FY14 <- CreateSeuratObject(counts=matrixFLT_LAR_YNG_FY14,min.cells = 3, min.features = 200)

setwd("~/MEMP/440/Project/PBMCs/FLT_LAR_YNG_FY11")
matrixFLT_LAR_YNG_FY11 = readMM('matrix.mtx')
rownames(matrixFLT_LAR_YNG_FY11) = read_tsv("features.tsv",col_names=FALSE)[, 2, drop=TRUE]
colnames(matrixFLT_LAR_YNG_FY11) = read_tsv("barcodes.tsv",col_names=FALSE)[, 1, drop=TRUE]
FLT_LAR_YNG_FY11 <- CreateSeuratObject(counts=matrixFLT_LAR_YNG_FY11,min.cells = 3, min.features = 200)

setwd("~/MEMP/440/Project/PBMCs/FLT_LAR_YNG_FY9")
matrixFLT_LAR_YNG_FY9 = readMM('matrix.mtx')
rownames(matrixFLT_LAR_YNG_FY9) = read_tsv("features.tsv",col_names=FALSE)[, 2, drop=TRUE]
colnames(matrixFLT_LAR_YNG_FY9) = read_tsv("barcodes.tsv",col_names=FALSE)[, 1, drop=TRUE]
FLT_LAR_YNG_FY9 <- CreateSeuratObject(counts=matrixFLT_LAR_YNG_FY9,min.cells = 3, min.features = 200)

setwd("~/MEMP/440/Project/PBMCs/FLT_LAR_OLD_FO20")
matrixFLT_LAR_OLD_FO20 = readMM('matrix.mtx')
rownames(matrixFLT_LAR_OLD_FO20) = read_tsv("features.tsv",col_names=FALSE)[, 2, drop=TRUE]
colnames(matrixFLT_LAR_OLD_FO20) = read_tsv("barcodes.tsv",col_names=FALSE)[, 1, drop=TRUE]
FLT_LAR_OLD_FO20 <- CreateSeuratObject(counts=matrixFLT_LAR_OLD_FO20,min.cells = 3, min.features = 200)

setwd("~/MEMP/440/Project/PBMCs/FLT_LAR_OLD_FO19")
matrixFLT_LAR_OLD_FO19 = readMM('matrix.mtx')
rownames(matrixFLT_LAR_OLD_FO19) = read_tsv("features.tsv",col_names=FALSE)[, 2, drop=TRUE]
colnames(matrixFLT_LAR_OLD_FO19) = read_tsv("barcodes.tsv",col_names=FALSE)[, 1, drop=TRUE]
FLT_LAR_OLD_FO19 <- CreateSeuratObject(counts=matrixFLT_LAR_OLD_FO19,min.cells = 3, min.features = 200)

setwd("~/MEMP/440/Project/PBMCs/FLT_LAR_OLD_FO16")
matrixFLT_LAR_OLD_FO16 = readMM('matrix.mtx')
rownames(matrixFLT_LAR_OLD_FO16) = read_tsv("features.tsv",col_names=FALSE)[, 2, drop=TRUE]
colnames(matrixFLT_LAR_OLD_FO16) = read_tsv("barcodes.tsv",col_names=FALSE)[, 1, drop=TRUE]
FLT_LAR_OLD_FO16 <- CreateSeuratObject(counts=matrixFLT_LAR_OLD_FO16,min.cells = 3, min.features = 200)

setwd("~/MEMP/440/Project/PBMCs/FLT_LAR_OLD_FO4")
matrixFLT_LAR_OLD_FO4 = readMM('matrix.mtx')
rownames(matrixFLT_LAR_OLD_FO4) = read_tsv("features.tsv",col_names=FALSE)[, 2, drop=TRUE]
colnames(matrixFLT_LAR_OLD_FO4) = read_tsv("barcodes.tsv",col_names=FALSE)[, 1, drop=TRUE]
FLT_LAR_OLD_FO4 <- CreateSeuratObject(counts=matrixFLT_LAR_OLD_FO4,min.cells = 3, min.features = 200)

#normalize data independently, 10000 reads/cell, natural log transformation 
GC_YNG_GY9 <- NormalizeData(GC_YNG_GY9)
GC_YNG_GY10 <- NormalizeData(GC_YNG_GY10)
GC_YNG_GY6 <- NormalizeData(GC_YNG_GY6)
GC_YNG_GY4 <- NormalizeData(GC_YNG_GY4)
GC_OLD_GO20 <- NormalizeData(GC_OLD_GO20)
GC_OLD_GO19 <- NormalizeData(GC_OLD_GO19)
GC_OLD_GO14 <- NormalizeData(GC_OLD_GO14)
GC_OLD_GO16 <- NormalizeData(GC_OLD_GO16)
FLT_LAR_YNG_FY15 <- NormalizeData(FLT_LAR_YNG_FY15)
FLT_LAR_YNG_FY14 <- NormalizeData(FLT_LAR_YNG_FY14)
FLT_LAR_YNG_FY11 <- NormalizeData(FLT_LAR_YNG_FY11)
FLT_LAR_YNG_FY9 <- NormalizeData(FLT_LAR_YNG_FY9)
FLT_LAR_OLD_FO20 <- NormalizeData(FLT_LAR_OLD_FO20)
FLT_LAR_OLD_FO19 <- NormalizeData(FLT_LAR_OLD_FO19)
FLT_LAR_OLD_FO16 <- NormalizeData(FLT_LAR_OLD_FO16)
FLT_LAR_OLD_FO4 <- NormalizeData(FLT_LAR_OLD_FO4)

#identify variable features independently, 2000 highest variable genes
GC_YNG_GY9 <- FindVariableFeatures(GC_YNG_GY9,selection.method = "vst", nfeatures = 2000)
GC_YNG_GY10 <- FindVariableFeatures(GC_YNG_GY10,selection.method = "vst", nfeatures = 2000)
GC_YNG_GY6 <- FindVariableFeatures(GC_YNG_GY6,selection.method = "vst", nfeatures = 2000)
GC_YNG_GY4 <- FindVariableFeatures(GC_YNG_GY4,selection.method = "vst", nfeatures = 2000)
GC_OLD_GO20 <- FindVariableFeatures(GC_OLD_GO20,selection.method = "vst", nfeatures = 2000)
GC_OLD_GO19 <- FindVariableFeatures(GC_OLD_GO19,selection.method = "vst", nfeatures = 2000)
GC_OLD_GO14 <- FindVariableFeatures(GC_OLD_GO14,selection.method = "vst", nfeatures = 2000)
GC_OLD_GO16 <- FindVariableFeatures(GC_OLD_GO16,selection.method = "vst", nfeatures = 2000)
FLT_LAR_YNG_FY15 <- FindVariableFeatures(FLT_LAR_YNG_FY15,selection.method = "vst", nfeatures = 2000)
FLT_LAR_YNG_FY14 <- FindVariableFeatures(FLT_LAR_YNG_FY14,selection.method = "vst", nfeatures = 2000)
FLT_LAR_YNG_FY11 <- FindVariableFeatures(FLT_LAR_YNG_FY11,selection.method = "vst", nfeatures = 2000)
FLT_LAR_YNG_FY9 <- FindVariableFeatures(FLT_LAR_YNG_FY9,selection.method = "vst", nfeatures = 2000)
FLT_LAR_OLD_FO20 <- FindVariableFeatures(FLT_LAR_OLD_FO20,selection.method = "vst", nfeatures = 2000)
FLT_LAR_OLD_FO19 <- FindVariableFeatures(FLT_LAR_OLD_FO19,selection.method = "vst", nfeatures = 2000)
FLT_LAR_OLD_FO16 <- FindVariableFeatures(FLT_LAR_OLD_FO16,selection.method = "vst", nfeatures = 2000)
FLT_LAR_OLD_FO4 <- FindVariableFeatures(FLT_LAR_OLD_FO4,selection.method = "vst", nfeatures = 2000)

# select features that are repeatedly variable across datasets for integration
dataList = c(GC_YNG_GY9,GC_YNG_GY10,GC_YNG_GY6,GC_YNG_GY4,GC_OLD_GO20,GC_OLD_GO19,GC_OLD_GO14,
             GC_OLD_GO16,FLT_LAR_YNG_FY15,FLT_LAR_YNG_FY14,FLT_LAR_YNG_FY11,FLT_LAR_YNG_FY9,
             FLT_LAR_OLD_FO20,FLT_LAR_OLD_FO19,FLT_LAR_OLD_FO16,FLT_LAR_OLD_FO4)
features <- SelectIntegrationFeatures(object.list = dataList, nfeatures=2000)

dataList <- lapply(X=dataList,FUN =function(x) {
  x <- ScaleData(x,features=features,verbose=FALSE)
  x<- RunPCA(x, features=features,verbose=FALSE)
})


#perform integration with rpca, chose 1 reference for each condition
immune.anchors <- FindIntegrationAnchors(object.list = dataList, reference = c(1,5,9,13),
                                         reduction = "rpca",dims=1:50)

# this command creates an 'integrated' data assay
immune.combined <- IntegrateData(anchorset = immune.anchors, dims=1:50)

#integrated analysis
#workflow for visualization and clustering
immune.combined <- ScaleData(immune.combined,verbose=FALSE)
immune.combined <- RunPCA(immune.combined,verbose=FALSE)
#t-sne and clustering
immune.combined <- RunUMAP(immune.combined,dims=1:50)
immune.combined <- FindNeighbors(immune.combined,reduction='pca',dims=1:50)
immune.combined <- FindClusters(immune.combined,resolution = 0.5)

#cell type assignment with sctype 
# load gene set preparation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

#load cell marker DB and gene sets
# DB file
db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
tissue = "Immune system" 

# prepare gene sets
gs_list = gene_sets_prepare(db_, tissue)

# get cell-type by cell matrix
es.max = sctype_score(scRNAseqData = immune.combined[["integrated"]]@scale.data, scaled = TRUE, 
                      gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 

# NOTE: scRNAseqData parameter should correspond to your input scRNA-seq matrix. 
# In case Seurat is used, it is either pbmc[["RNA"]]@scale.data (default), pbmc[["SCT"]]@scale.data, in case sctransform is used for normalization,
# or pbmc[["integrated"]]@scale.data, in case a joint analysis of multiple single-cell datasets is performed.

# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(immune.combined@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(immune.combined@meta.data[immune.combined@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(immune.combined@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
print(sctype_scores[,1:3])

#re visualize umap with new clusters
immune.combined@meta.data$customclassif = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  immune.combined@meta.data$customclassif[immune.combined@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

DimPlot(immune.combined, reduction = "umap", repel = TRUE, group.by = 'customclassif') +ggtitle("")     


#visualize differences in cell type 
#get cell type counts for each condition
#install.packages("remotes")
#remotes::install_github("genomics-kl/seurathelpeR")
library(seurathelpeR)

#relabel identification as sctype classification & get proportions
Idents(immune.combined) <- immune.combined$customclassif
counts = count_cells(immune.combined,'orig.ident')
pt <- table(Idents(immune.combined), immune.combined$orig.ident)
pt <- as.data.frame(pt)
pt$CellType <- as.character(pt$Var1)
pt$Total <- rep(c(counts[counts$orig.ident=='GC_Old',2],counts[counts$orig.ident=='GC_Young',2],
                  counts[counts$orig.ident=='SF_Old',2],counts[counts$orig.ident=='SF_Young',2]),
                each=length(levels(immune.combined)))
pt$Total <- as.numeric(pt$Total)
pt$Proportion <- pt$Freq / pt$Total
library(RColorBrewer)
ggplot(pt) + 
  list(
    geom_bar(aes(x=Var2, y=Proportion, fill=CellType), 
             colour="black", stat="identity") ,
    xlab(""),
    ylab("Proportion") ,
    labs(fill="Cell Type"))

#add significance testing on cell type prop
library(stats)
ptOld = subset(pt, subset=(Var2 == "GC_Old" | Var2=="SF_Old"))
ptYoung = subset(pt, subset=(Var2 == "GC_Young" | Var2=="SF_Young"))
ptOld$pVal = NA
for (type in unique(pt$Var1)){
  n = c(ptOld[ptOld$Var1== type & ptOld$Var2 == "SF_Old", "Total"],ptOld[ptOld$Var1== type & ptOld$Var2 == "GC_Old", "Total"])
  x = c(ptOld[ptOld$Var1== type & ptOld$Var2 == "SF_Old", "Freq"],ptOld[ptOld$Var1== type & ptOld$Var2 == "GC_Old", "Freq"])
  pval = prop.test(x,n)$p.value
  ptOld[ptOld$Var1==type,"pVal"] = pval
}
  
for (type in unique(pt$Var1)){
  n = c(ptYoung[ptYoung$Var1== type & ptYoung$Var2 == "SF_Young", "Total"],ptYoung[ptYoung$Var1== type & ptYoung$Var2 == "GC_Young", "Total"])
  x = c(ptYoung[ptYoung$Var1== type & ptYoung$Var2 == "SF_Young", "Freq"],ptYoung[ptYoung$Var1== type & ptYoung$Var2 == "GC_Young", "Freq"])
  pval = prop.test(x,n)$p.value
  ptYoung[ptYoung$Var1==type,"pVal"] = pval
}  

sigDiffCellsOld = ptOld[ptOld$pVal < .05,]
sigDiffCellsYoung = ptYoung[ptYoung$pVal < .05,]
sigDiffCellsOld[,1] = NULL
colnames(sigDiffCellsOld)[1] <- "Group"
sigDiffCellsYoung[,1] = NULL
colnames(sigDiffCellsYoung)[1] <- "Group"
write.xlsx(sigDiffCellsOld,file="CellTypeSigPropOld.xlsx", colNames=TRUE)
write.xlsx(sigDiffCellsYoung,file="CellTypeSigPropYoung.xlsx", colNames=TRUE)

#diff exp analysis between conditions
# Switch the identity of your cells to condition
Idents(immune.combined) <- 'orig.ident'
#find markers diff exp in spaceflight conditions
sfResponse_young <- FindMarkers(immune.combined, ident.1 = "SF_Young", ident.2 = "GC_Young", verbose = FALSE)
sfResponse_old <- FindMarkers(immune.combined, ident.1 = "SF_Old", ident.2 = "GC_Old", verbose = FALSE)
#filtering for significance
sfResponse_young = sfResponse_young[sfResponse_young$p_val_adj < 0.05,]
sfResponse_old = sfResponse_old[sfResponse_old$p_val_adj < 0.05,]
young_genes <- rownames(sfResponse_young)
youngSubset <- subset(immune.combined,subset = (orig.ident == "GC_Young" | orig.ident == "SF_Young"))

old_genes <- rownames(sfResponse_old)
oldSubset <- subset(immune.combined,subset = (orig.ident == "GC_Old" | orig.ident == "SF_Old"))

write.xlsx(sfResponse_young, file="DEgenesYoung.xlsx", rowNames=TRUE,colNames=TRUE)
write.xlsx(sfResponse_old, file="DEgenesOld.xlsx", rowNames=TRUE,colNames=TRUE)

#heatmap, look at cell type specific expression in shared genes 
sharedGenes <- young_genes[young_genes %in% old_genes]
DoHeatmap(oldSubset, features = sharedGenes, group.by = "customclassif", angle=90, size=4) + guides(color="none")
DoHeatmap(youngSubset, features = sharedGenes, group.by = "customclassif", label=FALSE) + guides(color="none")
DoHeatmap(youngSubset,features = rownames(sfResponse_young),group.by="orig.ident",label=FALSE, size=3) + guides(color="none")
DoHeatmap(oldSubset,features = rownames(sfResponse_old),group.by="orig.ident", size=3) + guides(color="none")

#enrichment
library(enrichR)
#disease enrichment for young mice
SpaceFlightHum <- toupper(young_genes)
enrSpaceFlightDisease <- enrichr(SpaceFlightHum, "DisGeNET")

enrSpaceFlightDisease[[1]]$P.value <- enrSpaceFlightDisease[[1]]$Adjusted.P.value
enrSpaceFlightDisease[[1]] <- enrSpaceFlightDisease$DisGeNET[enrSpaceFlightDisease$DisGeNET$Adjusted.P.value <0.05,]

plotEnrich(enrSpaceFlightDisease[[1]],showTerms = 20,numChar = 40,
           y = "Count", orderBy = "P.value", xlab = "Disease Profile", ylab = NULL, 
           title = "Young Mouse PBMCs in Spaceflight")
write.xlsx(enrSpaceFlightDisease[[1]],file="youngDisease.xlsx", colNames=TRUE)

#process enrichment in young mice
enrSpaceFlightProcess <- enrichr(SpaceFlightHum, "GO_Biological_Process_2021")

enrSpaceFlightProcess[[1]]$P.value <- enrSpaceFlightProcess[[1]]$Adjusted.P.value
enrSpaceFlightProcess[[1]] <- enrSpaceFlightProcess$GO_Biological_Process_2021[enrSpaceFlightProcess$GO_Biological_Process_2021$Adjusted.P.value <0.05,]

for (i in 1:length(enrSpaceFlightProcess[[1]]$Term)) {
  s <- enrSpaceFlightProcess[[1]]$Term[i]
  enrSpaceFlightProcess[[1]]$Term[i] <- unlist(strsplit(s, split='(', fixed=TRUE))[1]
}

interestingTerms <- c("innate immune response ","regulation of tumor necrosis factor production ",
                      "neutrophil mediated immunity ","response to interferon-gamma ","response to interferon-alpha ",
                      "negative regulation of myeloid cell differentiation ","response to cytokine ",
                      "negative regulation of cytokine production ","regulation of megakaryocyte differentiation ")
intEnrSFProcess <- enrSpaceFlightProcess$GO_Biological_Process_2021[is.element(enrSpaceFlightProcess$GO_Biological_Process_2021$Term,interestingTerms),]

plotEnrich(intEnrSFProcess,showTerms = 20,numChar = 100,
           y = "Count", orderBy = "P.value", xlab = "Biological Process", ylab = NULL, 
           title = "Young Mouse PBMCs in Spaceflight")
write.xlsx(enrSpaceFlightProcess[[1]],file="youngProcess.xlsx", colNames=TRUE)

#disease enrichment for old mice
SpaceFlightHumOld <- toupper(old_genes)
enrSpaceFlightDiseaseOld <- enrichr(SpaceFlightHumOld, "DisGeNET")

enrSpaceFlightDiseaseOld[[1]]$P.value <- enrSpaceFlightDiseaseOld[[1]]$Adjusted.P.value
enrSpaceFlightDiseaseOld[[1]] <- enrSpaceFlightDiseaseOld$DisGeNET[enrSpaceFlightDiseaseOld$DisGeNET$Adjusted.P.value <0.05,]

diseaseOI = c("Myocarditis","Arteriosclerosis", "Atherosclerosis","Endothelial dysfunction","Cardiovascular Diseases",
              "Coronary Artery Disease","External Carotid Artery Diseases","Internal Carotid Artery Diseases", 
              "Atherosclerotic occlusive disease","Acute myocardial infarction", "Acute Coronary Syndrome",
              "Coronary Arteriosclerosis", "Coronary heart disease", "Atherosclerosis of aorta")
enrSFDiseaseOld_Filt <- enrSpaceFlightDiseaseOld$DisGeNET[is.element(enrSpaceFlightDiseaseOld$DisGeNET$Term, diseaseOI),]
plotEnrich(enrSFDiseaseOld_Filt,showTerms = 20,numChar = 100,
           y = "Count", orderBy = "P.value", xlab = "Disease Profile", ylab = NULL, 
           title = "Old Mouse PBMCs in Spaceflight")
write.xlsx(enrSpaceFlightDiseaseOld[[1]],file="oldDisease.xlsx", colNames=TRUE)

#process enrichment in old mice
enrSpaceFlightProcessOld <- enrichr(SpaceFlightHumOld, "GO_Biological_Process_2021")

enrSpaceFlightProcessOld[[1]]$P.value <- enrSpaceFlightProcessOld[[1]]$Adjusted.P.value
enrSpaceFlightProcessOld[[1]] <- enrSpaceFlightProcessOld$GO_Biological_Process_2021[enrSpaceFlightProcessOld$GO_Biological_Process_2021$Adjusted.P.value <0.05,]

for (i in 1:length(enrSpaceFlightProcessOld[[1]]$Term)) {
  s <- enrSpaceFlightProcessOld[[1]]$Term[i]
  enrSpaceFlightProcessOld[[1]]$Term[i] <- unlist(strsplit(s, split='(', fixed=TRUE))[1]
}

processesOI = c("neutrophil degranulation ", "neutrophil mediated immunity ","regulation of inflammatory response ",
"innate immune response ", "alpha-beta T cell differentiation ", "negative regulation of cholesterol metabolic process ", 
"positive regulation of prostaglandin biosynthetic process ","regulation of low-density lipoprotein particle receptor catabolic process ",
"regulation of cytokine production ","positive regulation of NF-kappaB transcription factor activity ",
"regulation of T-helper 1 cell differentiation ","negative regulation of type 2 immune response ",
"regulation of T-helper 2 cell differentiation ","alpha-beta T cell activation ",
"leukocyte aggregation ", "positive regulation of cholesterol esterification ",
"positive regulation of cholesterol esterification ", "cytokine-mediated signaling pathway ",
"regulation of interleukin-1 production ", "negative regulation of interleukin-8 production ",
"positive regulation of interleukin-2 production ","negative regulation of blood vessel endothelial cell migration ",
"cellular response to vascular endothelial growth factor stimulus ", "negative regulation of blood coagulation ", 
"monocyte chemotaxis")

enrSFProcessOld_Filt <- enrSpaceFlightProcessOld$GO_Biological_Process_2021[is.element(enrSpaceFlightProcessOld$GO_Biological_Process_2021$Term, processesOI),]
plotEnrich(enrSFProcessOld_Filt,showTerms = 25,numChar = 100,
           y = "Count", orderBy = "P.value", xlab = "Biological Process", ylab = NULL, 
           title = "Old Mouse PBMCs in Spaceflight")
write.xlsx(enrSpaceFlightProcessOld[[1]],file="oldProcess.xlsx", colNames=TRUE)

#run NicheNet
setwd("~/MEMP/440/Project/PBMCs")
library(nichenetr)
library(tidyverse)
url = "https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"
download.file(url,"lig_targ_mat.rds", method="curl")
ligand_target_matrix = readRDS("lig_targ_mat.rds")

#first do with all features
#merge young SF mice expressed features
senderFeatures <- c(rownames(FLT_LAR_YNG_FY9),rownames(FLT_LAR_YNG_FY11),rownames(FLT_LAR_YNG_FY14),rownames(FLT_LAR_YNG_FY15))
#remove duplicates
senderFeatures <- senderFeatures[!duplicated(senderFeatures)]
senderFeatures <- toupper(senderFeatures)

#define sender, receiver, and interesting genes
expressed_genes_sender = senderFeatures
expressed_genes_receiver = read.table("heart_expressed_genes_caps.txt", sep='\t')
expressed_genes_receiver = expressed_genes_receiver$V1
geneset_oi = read.table("geneset_oi.txt",sep="\t")
geneset_oi = geneset_oi$V1

#define background genes
background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

#define potential ligands
lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
ligands = lr_network %>% pull(from) %>% unique()
expressed_ligands = intersect(ligands,expressed_genes_sender)

receptors = lr_network %>% pull(to) %>% unique()
expressed_receptors = intersect(receptors,expressed_genes_receiver) 

lr_network_expressed = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) 
potential_ligands = lr_network_expressed %>% pull(from) %>% unique()

#ligand activity analysis
ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
#rank based on activity 
best_upstream_ligands = ligand_activities %>% top_n(20, pearson) %>% arrange(-pearson) %>% pull(test_ligand)


# show histogram of ligand activity scores
p_hist_lig_activity = ggplot(ligand_activities, aes(x=pearson)) + 
  geom_histogram(color="black", fill="darkorange")  + 
  # geom_density(alpha=.1, fill="orange") +
  geom_vline(aes(xintercept=min(ligand_activities %>% top_n(20, pearson) %>% pull(pearson))), color="red", linetype="dashed", size=1) + 
  labs(x="Ligand Activity (PCC)", y = "# Ligands") +
  theme_classic()
p_hist_lig_activity

#infer target genes
active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 250) %>% bind_rows()
active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.25)
order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets = active_ligand_target_links_df$target %>% unique()
vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()

p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot("Prioritized PBMC-Ligands","DE & CVD Genes in Heart Cells", color = "purple",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential") + scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.005,0.01)) + theme(axis.text.x = element_text(face = "italic"))

p_ligand_target_network

#network inference 
# get the ligand-receptor network of the top-ranked ligands
lr_network_top = lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()

# get the weights of the ligand-receptor interactions as used in the NicheNet model
url ="https://zenodo.org/record/3260758/files/weighted_networks.rds"
download.file(url,"weightedNetworks.rds", method="curl")
weighted_networks = readRDS("weightedNetworks.rds")
lr_network_top_df = weighted_networks$lr_sig %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)

# convert to a matrix
lr_network_top_df = lr_network_top_df %>% spread("from","weight",fill = 0)
lr_network_top_matrix = lr_network_top_df %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)

# perform hierarchical clustering to order the ligands and receptors
dist_receptors = dist(lr_network_top_matrix, method = "binary")
hclust_receptors = hclust(dist_receptors, method = "ward.D2")
order_receptors = hclust_receptors$labels[hclust_receptors$order]

dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
hclust_ligands = hclust(dist_ligands, method = "ward.D2")
order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]

vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]
#subset to only include rows and columns with potential >.5
Filtvis_ligand_receptor_network <- vis_ligand_receptor_network
#filter for ligands in gene target plot 
Filtvis_ligand_receptor_network = Filtvis_ligand_receptor_network[,is.element(colnames(Filtvis_ligand_receptor_network),rownames(vis_ligand_target))]
Filtvis_ligand_receptor_network =  Filtvis_ligand_receptor_network[rowSums(Filtvis_ligand_receptor_network > 0.5) >= 1,]
Filtvis_ligand_receptor_network2 = Filtvis_ligand_receptor_network[,rownames(vis_ligand_target)]

p_ligand_receptor_network = Filtvis_ligand_receptor_network2 %>% t() %>% make_heatmap_ggplot("Prioritized PBMC-Ligands","Receptors Expressed by Heart Cells", color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential") +
            scale_fill_gradient2(low = "whitesmoke",  high = "mediumvioletred") +
             theme(axis.text = element_text(size = 6)) + guides(fill = guide_colourbar(
              label.theme = element_text(size=8)))
p_ligand_receptor_network

library(RColorBrewer)
library(cowplot)
library(ggpubr)
ligand_pearson_matrix = ligand_activities %>% select(pearson) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities$test_ligand)
vis_ligand_pearson = ligand_pearson_matrix[order_ligands, ] %>% as.matrix(ncol = 1) %>% magrittr::set_colnames("Pearson")
p_ligand_pearson = vis_ligand_pearson %>% make_heatmap_ggplot("Prioritized PBMC-ligands","Ligand activity", color = "darkorange",legend_position = "top", x_axis_position = "top", legend_title = "Pearson correlation coefficient\ntarget gene prediction ability)")+
   guides(fill = guide_colourbar(
    label.theme = element_text(size=8))) + theme(axis.text.x = element_blank())
p_ligand_pearson

#repeat for old mice
senderFeatures <- c(rownames(FLT_LAR_OLD_FO16),rownames(FLT_LAR_OLD_FO19),rownames(FLT_LAR_OLD_FO20),rownames(FLT_LAR_OLD_FO4))
#remove duplicates
senderFeatures <- senderFeatures[!duplicated(senderFeatures)]
senderFeatures <- toupper(senderFeatures)

#define sender, receiver, and interesting genes
expressed_genes_sender = senderFeatures
expressed_ligands = intersect(ligands,expressed_genes_sender)
lr_network_expressed = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) 
potential_ligands = lr_network_expressed %>% pull(from) %>% unique()

#ligand activity analysis
ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
#rank based on activity 
best_upstream_ligands = ligand_activities %>% top_n(20, pearson) %>% arrange(-pearson) %>% pull(test_ligand)

# show histogram of ligand activity scores
p_hist_lig_activity = ggplot(ligand_activities, aes(x=pearson)) + 
  geom_histogram(color="black", fill="darkorange")  + 
  # geom_density(alpha=.1, fill="orange") +
  geom_vline(aes(xintercept=min(ligand_activities %>% top_n(20, pearson) %>% pull(pearson))), color="red", linetype="dashed", size=1) + 
  labs(x="Ligand Activity (PCC)", y = "# Ligands") +
  theme_classic()
p_hist_lig_activity

#infer target genes
active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 250) %>% bind_rows()
active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.25)
order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets = active_ligand_target_links_df$target %>% unique()
vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()

p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot("Prioritized PBMC-Ligands","DE & CVD Genes in Heart Cells", color = "purple",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential") + scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.005,0.01)) + theme(axis.text.x = element_text(face = "italic"))

p_ligand_target_network

#network inference 
# get the ligand-receptor network of the top-ranked ligands
lr_network_top = lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()

# get the weights of the ligand-receptor interactions as used in the NicheNet model
lr_network_top_df = weighted_networks$lr_sig %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)

# convert to a matrix
lr_network_top_df = lr_network_top_df %>% spread("from","weight",fill = 0)
lr_network_top_matrix = lr_network_top_df %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)

# perform hierarchical clustering to order the ligands and receptors
dist_receptors = dist(lr_network_top_matrix, method = "binary")
hclust_receptors = hclust(dist_receptors, method = "ward.D2")
order_receptors = hclust_receptors$labels[hclust_receptors$order]

dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
hclust_ligands = hclust(dist_ligands, method = "ward.D2")
order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]

vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]
#subset to only include rows and columns with potential >.5
Filtvis_ligand_receptor_network <- vis_ligand_receptor_network
Filtvis_ligand_receptor_network =  Filtvis_ligand_receptor_network[rowSums(Filtvis_ligand_receptor_network > 0.5) >= 1, colSums(Filtvis_ligand_receptor_network > 0.5) >= 1]
#Filtvis_ligand_receptor_network[Filtvis_ligand_receptor_network < 0.5] <- NA
p_ligand_receptor_network = Filtvis_ligand_receptor_network %>% t() %>% make_heatmap_ggplot("Prioritized PBMC-Ligands","Receptors Expressed by Heart Cells", color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential") +
  theme(axis.text = element_text(size = 6)) + guides(fill = guide_colourbar(
    label.theme = element_text(size=8)))
p_ligand_receptor_network

library(RColorBrewer)
library(cowplot)
library(ggpubr)
ligand_pearson_matrix = ligand_activities %>% select(pearson) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities$test_ligand)
vis_ligand_pearson = ligand_pearson_matrix[order_ligands, ] %>% as.matrix(ncol = 1) %>% magrittr::set_colnames("Pearson")
p_ligand_pearson = vis_ligand_pearson %>% make_heatmap_ggplot("Prioritized CAF-ligands","Ligand activity", color = "darkorange",legend_position = "top", x_axis_position = "top", legend_title = "Pearson correlation coefficient\ntarget gene prediction ability)")+
  guides(fill = guide_colourbar(
    label.theme = element_text(size=8))) + theme(axis.text.x = element_blank())
p_ligand_pearson