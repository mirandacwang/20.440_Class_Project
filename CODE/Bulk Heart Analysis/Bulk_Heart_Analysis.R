#data analyzed on 04/28/23
#data is from a raw counts table provided on the 
#nasa open data portal from the RR-3 mission GLDS-270
#data is processed using EdgeR as detailed here: https://www.nathalievialaneix.eu/doc/html/solution_edgeR-tomato-withcode.html#starting-from-count-table
#the raw counts table and design table are provided in the github repository
BiocManager::install("mixOmics")
library(edgeR)
library(limma)
library(RColorBrewer)
library(mixOmics)
library(HTSFilter)
#read in the raw counts table
raw_counts_hrt <- read.csv('~/Desktop/RR3/GLDS-270_rna_seq_RSEM_Unnormalized_Counts.csv')

#make gene ensemble IDs the row names
rownames(raw_counts_hrt) <- raw_counts_hrt$X
raw_counts_hrt$X <- NULL
head(raw_counts_hrt)

#read in the design table
sample_info_hrt <- read.csv('~/Desktop/RR3/design_RR3.csv')
#create a DGEList data object
dgeFull_hrt <- DGEList(raw_counts_hrt, group=sample_info_hrt$condition)

#data exploration and quality assessment
pseudoCounts_hrt <- log2(dgeFull_hrt$counts+1)
head(pseudoCounts_hrt)
hist(pseudoCounts_hrt[,"RR3_HRT_FLT_F1"])
boxplot(pseudoCounts_hrt, col="gray", las=3)
plotMDS(pseudoCounts_hrt)
sampleDists_hrt <- as.matrix(dist(t(pseudoCounts_hrt)))
cimColor <- colorRampPalette(rev(brewer.pal(9, "Blues")))(16)
cim(sampleDists_hrt, color=cimColor, symkey=FALSE)

#DEanalysis
dgeFull_hrt <- DGEList(dgeFull_hrt$counts[apply(dgeFull_hrt$counts, 1, sum) !=0, ], group=dgeFull_hrt$samples$group)
head(dgeFull_hrt$counts)
dgeFull_hrt <- calcNormFactors(dgeFull_hrt, method='TMM')
eff.lib.size_hrt <- dgeFull_hrt$samples$lib.size*dgeFull_hrt$samples$norm.factors
normCounts_hrt <- cpm(dgeFull_hrt)
pseudoNormCounts_hrt <- log2(normCounts_hrt + 1)
boxplot(pseudoNormCounts_hrt, col='gray', las=3)
plotMDS(pseudoNormCounts_hrt)
dgeFull_hrt <- estimateCommonDisp(dgeFull_hrt)
dgeFull_hrt <- estimateTagwiseDisp(dgeFull_hrt)
dgeTest_hrt <- exactTest(dgeFull_hrt)
filtData_hrt <- HTSFilter(dgeFull_hrt)$filteredData
dgeTestFilt_hrt <- exactTest(filtData_hrt)
resNoFilt_hrt <- topTags(dgeTest_hrt, n=nrow(dgeTest_hrt$table))
resFilt_hrt <- topTags(dgeTestFilt_hrt, n=nrow(dgeTest_hrt$table))
sigDownReg_hrt <- resFilt_hrt$table[resFilt_hrt$table$FDR<.05,]
sigDownReg_hrt <- sigDownReg_hrt[order(sigDownReg_hrt$logFC),]
sigUpReg_hrt <- sigDownReg_hrt[order(sigDownReg_hrt$logFC, decreasing=TRUE),]

#plotting DEanalysis
plotSmear(dgeTestFilt_hrt,
          de.tags = rownames(resFilt_hrt$table)[which(resFilt_hrt$table$FDR<0.01)])

volcanoData <- cbind(resFilt_hrt$table$logFC, -log10(resFilt_hrt$table$FDR))
colnames(volcanoData) <- c("logFC", "negLogPval")
head(volcanoData)
plot(volcanoData, pch=19)


y_hrt <- cpm(dgeFull_hrt, log=TRUE, prior.count=1)
selY_hrt <- y_hrt[rownames(resFilt_hrt$table)[resFilt_hrt$table$FDR<.05],]
cimColor <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)[255:1]
finalHM_hrt <- cim(t(selY_hrt), color=cimColor, symkey=FALSE)

#GSEA analysis
#uses the package clusterProfiler and the visualization enrichPlot
#tutorial: http://yulab-smu.top/biomedical-knowledge-mining-book/universal-api.html
#another good tutorial: https://alexslemonade.github.io/refinebio-examples/03-rnaseq/pathway-analysis_rnaseq_02_gsea.html
gene_list_hrt <- resFilt_hrt$table
gene_list_hrt <- gene_list_hrt$logFC
names(gene_list_hrt) <- row.names(resFilt_hrt$table)
gene_list_hrt = sort(gene_list_hrt, decreasing = TRUE)
library(msigdbr)
library(clusterProfiler)
library(enrichplot)
hallmark_mouse <- msigdbr(species = "Mus musculus", category = "H") %>%
  dplyr::select(gs_name, ensembl_gene)
gse_result_hrt <- GSEA(gene_list_hrt, TERM2GENE=hallmark_mouse, pvalueCutoff = .05, pAdjustMethod="BH")
head(gse_result_hrt)
require(DOSE)
pdf(file = '/Users/wang/Desktop/Figures/heart_GSEA.pdf',
    width =10, 
    height=6)
dotplot(gse_result_hrt, title = "Space vs Ground GSEA results", showCategory=10, split='.sign') + facet_grid(.~.sign)
dev.off()

#KEGG Enrichment Analysis
#kegg tutorial https://learn.gencore.bio.nyu.edu/rna-seq-analysis/gene-set-enrichment-analysis/
#prepare for KEGG enrichment analysis
de_table <- resFilt_hrt$table
de_table$ENSEMBL <- rownames(de_table)
de_gene_list <- de_table$logFC
names(de_gene_list) <- de_table$ENSEMBL
de_gene_list = sort(de_gene_list, decreasing = TRUE)
ids <- bitr(names(de_gene_list), fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
dedup_ids = ids[!duplicated(ids[c("ENSEMBL")]),]
mapped_df = de_table[de_table$ENSEMBL %in% dedup_ids$ENSEMBL,]
mapped_df$ENTREZID = dedup_ids$ENTREZID
kegg_gene_list <- mapped_df$logFC
names(kegg_gene_list) <- mapped_df$ENTREZID
kegg_gene_list <- na.omit(kegg_gene_list)
kegg_gene_list = sort(kegg_gene_list, decreasing=TRUE)
kegg_organism = "mmu"
kk2 <- gseKEGG(geneList = kegg_gene_list, organism=kegg_organism, pvalueCutoff = .05, pAdjustMethod = 'none', keyType = "ncbi-geneid")
library(pathview)
#create images of enriched KEGG pathways
#b cell signaling
dme <- pathview(gene.data=kegg_gene_list, pathway.id="mmu04662", species=kegg_organism)
knitr::include_graphics("mmu04662.pathview.png")
kegg_data <- data.frame(kk2)
#cell adhesion
dme <- pathview(gene.data=kegg_gene_list, pathway.id="mmu04514", species=kegg_organism)
knitr::include_graphics("mmu04514.pathview.png")
#cell senescence
dme <- pathview(gene.data=kegg_gene_list, pathway.id="mmu04218", species=kegg_organism)
knitr::include_graphics("mmu04218.pathview.png")

#cell-ECm interactions
dme <- pathview(gene.data=kegg_gene_list, pathway.id="mmu04512", species=kegg_organism)
knitr::include_graphics("mmu04512.pathview.png")

#IL-17 signaling pathway
dme <- pathview(gene.data=kegg_gene_list, pathway.id="mmu04657", species=kegg_organism)
knitr::include_graphics("mmu04657.pathview.png")

#PD-L1
dme <- pathview(gene.data=kegg_gene_list, pathway.id="mmu05235", species=kegg_organism)
knitr::include_graphics("mmu05235.pathview.png")



#reactome analysis
BiocManager::install("org.Mm.eg.db")
library(org.Mm.eg.db)
#add a column of gene short names
infile_hrt <- rownames(sigDownReg_hrt)
data_hrt = as.vector(infile_hrt)
annots_hrt <- select(org.Mm.eg.db, keys=data_hrt, columns="SYMBOL", keytype="ENSEMBL")
entrez_hrt <- select(org.Mm.eg.db, keys=data_hrt, columns="ENTREZID", keytype="ENSEMBL")
gene_symbols_hrt <- annots_hrt$SYMBOL
entrez_id_hrt <- entrez_hrt$ENTREZID
sigDownReg_hrt$SYMBOL <- gene_symbols_hrt
sigDownReg_hrt$ENTREZID <- entrez_id_hrt
sigDownReg_hrt <- na.omit(sigDownReg_hrt)


library(ReactomePA)

reactome_gene_list <- sigDownReg_hrt$logFC
names(reactome_gene_list) <- sigDownReg_hrt$ENTREZID
reactome_gene_list <- na.omit(reactome_gene_list)
reactome_gene_list = sort(reactome_gene_list, decreasing=TRUE)
de_reactome <- names(reactome_gene_list)
reactome_results <- enrichPathway(gene=de_reactome, pvalueCutoff=.05, readable=TRUE, organism="mouse")
viewPathway("Cellular responses to stimuli", readable=TRUE)


#generate table for cibersort
hrt_counts_ciber <- raw_counts_hrt[rowSums(raw_counts_hrt[])>0,]
hrt_ensembls <- rownames(hrt_counts_ciber)
hrt_counts_ciber$ENSEMBL = hrt_ensembls
library(clusterProfiler)
hrt_symbols <- bitr(hrt_ensembls, fromType="ENSEMBL", toType="SYMBOL", OrgDb="org.Mm.eg.db")
dedup_symbols = hrt_symbols[!duplicated(hrt_symbols[c("ENSEMBL")]),]
mapped_ciber = hrt_counts_ciber[hrt_counts_ciber$ENSEMBL %in% dedup_symbols$ENSEMBL,]
mapped_ciber$gene_symbol = dedup_symbols$SYMBOL
mapped_ciber <- subset(mapped_ciber, select = -c(ENSEMBL))
cibersort_df <- mapped_ciber %>% select('gene_symbol', everything())
cibersort <- cibersort_df[!duplicated(cibersort_df$gene_symbol),]
write.table(cibersort, file = "heart_expression.txt", sep = "\t",
            row.names = FALSE, quote=FALSE)


#list of expressed genes for NicheNet Analysis
expressed_genes <- rownames(filtData_hrt$counts)
expr_genes = as.vector(expressed_genes)
library(org.Mm.eg.db)
annot_expr <- select(org.Mm.eg.db, keys=expr_genes, columns="SYMBOL", keytype="ENSEMBL")
expr_genes_symbols <- annot_expr$SYMBOL
expr_genes_symbols <- unique(expr_genes_symbols)
write.table(expr_genes_symbols, file = "reciever_expressed_genes.txt", sep = "\t",
            row.names = FALSE,)
write.table(gene_symbols_hrt, file = "diff_exp_genes.txt", sep = "\t",
            row.names = FALSE,)
sigDownReg_hrt$SYMBOL <- gene_symbols_hrt


infile_hrt <- rownames(sigDownReg_hrt)
data_hrt = as.vector(infile_hrt)
annots_hrt <- select(org.Mm.eg.db, keys=data_hrt, columns="SYMBOL", keytype="ENSEMBL")
entrez_hrt <- select(org.Mm.eg.db, keys=data_hrt, columns="ENTREZID", keytype="ENSEMBL")
gene_symbols_hrt <- annots_hrt$SYMBOL
entrez_id_hrt <- entrez_hrt$ENTREZID
sigDownReg_hrt$SYMBOL <- gene_symbols_hrt
sigDownReg_hrt$ENTREZID <- entrez_id_hrt
sigDownReg_hrt <- na.omit(sigDownReg_hrt)

library(clusterProfiler)
library(org.Mm.eg.db)
BiocManager::install('EnhancedVolcano')

#creates volcanoplot of differentially expressed genes
volcano_table <- resFilt_hrt$table
volcano_table$ENSEMBL <- rownames(volcano_table)
gene_list <- volcano_table$ENSEMBL
symbols <- bitr(gene_list, fromType="ENSEMBL", toType="SYMBOL", OrgDb = 'org.Mm.eg.db')
dedup_symbols = symbols[!duplicated(symbols[c("ENSEMBL")]),]
mapped_volcano = volcano_table[volcano_table$ENSEMBL %in% dedup_symbols$ENSEMBL,]
mapped_volcano$SYMBOL = dedup_symbols$SYMBOL
library(EnhancedVolcano)
EnhancedVolcano(mapped_volcano, title = "Space vs Ground",legendLabels=c( 'adjPvalue','Not sig.','adjPvalue','adjPvalue'), lab=mapped_volcano$SYMBOL,selectLab =c('Dnajb1', "Ifi44", 'Stat1', 'Irf7', 'Cd79a', 'Cldn5', 'Fasn', 'Adamtsl2', 'Fam177a2', 'Alb', 'Car3', 'Mat1a', 'Myl4', 'Myl7'),  x='logFC', y='FDR', FCcutoff = 0, pCutoff = .05,drawConnectors = TRUE, max.overlaps = Inf)
options(ggrepel.max.overlaps = Inf)

pdf(file = '/Users/wang/Desktop/Figures/heart_volcano.pdf',
    width =10, 
    height=10)
EnhancedVolcano(mapped_volcano, title = "Space vs Ground",legendLabels=c( 'adjPvalue','Not sig.','adjPvalue','adjPvalue'), lab=mapped_volcano$SYMBOL,selectLab =c('Dnajb1', 'Stat1', 'Irf7', 'Adamtsl2', 'Fam177a2', 'Alb', 'Car3', 'Mat1a', 'Myl4', 'Myl7'),  x='logFC', y='FDR', FCcutoff = 0, pCutoff = .05,drawConnectors = TRUE, max.overlaps = Inf, labSize=7)
dev.off()

#generates visual representation of enriched KEGG pathways
library(enrichplot)
require(DOSE)
dotplot(kk2, title = "Space vs Ground KEGG Pathways", showCategory=10, split='.sign') + facet_grid(.~.sign)
de_table$ENSEMBL <- rownames(de_table)
de_gene_list <- de_table$logFC
names(de_gene_list) <- de_table$ENSEMBL
de_gene_list = sort(de_gene_list, decreasing = TRUE)
ids <- bitr(names(de_gene_list), fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
dedup_ids = ids[!duplicated(ids[c("ENSEMBL")]),]
mapped_df = de_table[de_table$ENSEMBL %in% dedup_ids$ENSEMBL,]
mapped_df$ENTREZID = dedup_ids$ENTREZID
keg_cols <- c('ID', 'Description', 'enrichmentScore', 'p.adjust')
kegg_pretty <- kegg_data[, colnames(kegg_data) %in% keg_cols]
write.csv(kegg_pretty, file = "kegg_heart.csv",
            row.names = FALSE, quote=FALSE)
library(enrichplot)
pdf(file = '/Users/wang/Desktop/Figures/heart_kegg.pdf',
    width =7, 
    height=6)
dotplot(kk2, title = "Enriched KEGG Pathways", showCategory=10, split='.sign', font.size = 8) + facet_grid(.~.sign)
dev.off()
write.csv(sigDownReg_hrt, 'Heart_DE_genes.csv', row.names=FALSE)


