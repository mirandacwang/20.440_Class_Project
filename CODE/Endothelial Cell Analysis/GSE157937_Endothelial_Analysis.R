
#data is from a raw counts table obtained through GREIN
#data is processed using EdgeR as detailed here: https://www.nathalievialaneix.eu/doc/html/solution_edgeR-tomato-withcode.html#starting-from-count-table
BiocManager::install("mixOmics")
library(edgeR)
library(limma)
library(RColorBrewer)
library(mixOmics)
library(HTSFilter)
#read in the raw counts table
raw_counts_hmec <- read.csv('~/Desktop/GREIN/GSE157937_GeneLevel_Raw_data.csv')

#remove extra gene symbol column and space 1g columns
raw_counts_hmec <- subset(raw_counts_hmec, select = -c(GSM4781325,GSM4781326, gene_symbol))
#turn ensembl IDs into rownames
raw_counts_hmec <- raw_counts_hmec[!duplicated(raw_counts_hmec$X),]
rownames(raw_counts_hmec) <- raw_counts_hmec$X
raw_counts_hmec$X <- NULL
head(raw_counts_hmec)

#read in the design table
sample_info_hmec <- read.table('~/Desktop/GREIN/design_table.csv', header=TRUE, sep=",", row.names=1)
#create a DGEList data object
dgeFull_hmec <- DGEList(raw_counts_hmec, group=sample_info_hmec$condition)

#data exploration and quality assessment
pseudoCounts_hmec <- log2(dgeFull_hmec$counts+1)
head(pseudoCounts_hmec)
hist(pseudoCounts_hmec[,"GSM4781322"])
boxplot(pseudoCounts_hmec, col="gray", las=3)
plotMDS(pseudoCounts_hmec)
sampleDists_hmec <- as.matrix(dist(t(pseudoCounts_hmec)))
cimColor <- colorRampPalette(rev(brewer.pal(9, "Blues")))(16)
cim(sampleDists_hmec, color=cimColor, symkey=FALSE)

#DEanalysis
dgeFull_hmec <- DGEList(dgeFull_hmec$counts[apply(dgeFull_hmec$counts, 1, sum) !=0, ], group=dgeFull_hmec$samples$group)
head(dgeFull_hmec$counts)
dgeFull_hmec <- calcNormFactors(dgeFull_hmec, method='TMM')
eff.lib.size_hmec <- dgeFull_hmec$samples$lib.size*dgeFull_hmec$samples$norm.factors
normCounts_hmec <- cpm(dgeFull_hmec)
pseudoNormCounts_hmec <- log2(normCounts_hmec + 1)
boxplot(pseudoNormCounts_hmec, col='gray', las=3)
plotMDS(pseudoNormCounts_hmec)
dgeFull_hmec <- estimateCommonDisp(dgeFull_hmec)
dgeFull_hmec <- estimateTagwiseDisp(dgeFull_hmec)
dgeTest_hmec <- exactTest(dgeFull_hmec)
filtData_hmec <- HTSFilter(dgeFull_hmec)$filteredData
dgeTestFilt_hmec <- exactTest(filtData_hmec)
resNoFilt_hmec <- topTags(dgeTest_hmec, n=nrow(dgeTest_hmec$table))
resFilt_hmec <- topTags(dgeTestFilt_hmec, n=nrow(dgeTest_hmec$table))
sigDownReg_hmec <- resFilt_hmec$table[resFilt_hmec$table$FDR<.05,]
sigDownReg_hmec <- sigDownReg_hmec[order(sigDownReg_hmec$logFC),]
sigUpReg_hmec <- sigDownReg_hmec[order(sigDownReg_hmec$logFC, decreasing=TRUE),]

#plotting DEanalysis
plotSmear(dgeTestFilt_hmec,
          de.tags = rownames(resFilt_hmec$table)[which(resFilt_hmec$table$FDR<0.01)])

volcanoData <- cbind(resFilt_hmec$table$logFC, -log10(resFilt_hmec$table$FDR))
colnames(volcanoData) <- c("logFC", "negLogPval")
head(volcanoData)
plot(volcanoData, pch=19)


y_hmec <- cpm(dgeFull_hmec, log=TRUE, prior.count=1)
selY_hmec <- y_hmec[rownames(resFilt_hmec$table)[resFilt_hmec$table$FDR<.05],]
cimColor <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)[255:1]
finalHM_hmec <- cim(t(selY_hmec), color=cimColor, symkey=FALSE)

#GSEA analysis
#uses the package clusterProfiler and the visualization enrichPlot
#tutorial: http://yulab-smu.top/biomedical-knowledge-mining-book/universal-api.html
#another good tutorial: https://alexslemonade.github.io/refinebio-examples/03-rnaseq/pathway-analysis_rnaseq_02_gsea.html
gene_list_hmec <- resFilt_hmec$table
gene_list_hmec <- gene_list_hmec$logFC
names(gene_list_hmec) <- row.names(resFilt_hmec$table)
gene_list_hmec = sort(gene_list_hmec, decreasing = TRUE)


library(msigdbr)
library(clusterProfiler)
library(enrichplot)
hallmark <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, ensembl_gene)
gse_result_hmec <- GSEA(gene_list_hmec, TERM2GENE=hallmark, pvalueCutoff = .05, pAdjustMethod="BH")
head(gse_result_hmec)
require(DOSE)
pdf(file = '/Users/wang/Desktop/GSEA_microvascular.pdf',
    width =10, 
    height=5)
dotplot(gse_result_hmec, title = "Space vs Ground GSEA results", showCategory=10, split='.sign') + facet_grid(.~.sign)
dev.off()

dotplot(gse_result_hmec)
