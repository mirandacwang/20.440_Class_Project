#sample analysis of provided normalized sphinx data on april first
#returns no significant results with correction
BiocManager::install("GO.db")
BiocManager::install("annotate")
BiocManager::install("clusterProfiler")
BiocManager::install("pathview")
BiocManager::install("enrichplot")
BiocManager::install("msigdbr")
library(msigdbr)
library(enrichplot)
library(ggplot2)
library(clusterProfiler)
library(GeneAnswers)
library(annotate)
library(affy)
library(Biobase)
library(limma)
library(readxl)
library(edgeR)
library(Glimma)
library(org.Hs.eg.db)
library(RColorBrewer)
library(gplots)
library(GO.db)
library(pheatmap)
#generating basic heatmap
pheatmap(eset)
object <- new("ExpressionSet", exprs=as.matrix(eset))
expression_set <- object
normalized_data <- as.matrix(eset)
#colnames(normalized_data) <- c("Flight Rep 1", "Flight Rep 2", "Flight Rep 3", "Ground Rep 1", "Ground Rep 2", "Ground Rep 3")
head(normalized_data)
var_genes <- apply(normalized_data, 2, var)
head(var_genes)
select_var <- names(sort(var_genes, decreasing=TRUE))[1:100]
head(select_var)
highly_variable_lcpm <- normalized_data[select_var,]
dim(highly_variable_lcpm)
head(highly_variable_lcpm)

pdf(file = '/Users/wang/Desktop/Figures/highly_variable_heatmap1.pdf',
    width =10, 
    height=10)
pheatmap(highly_variable_lcpm, fontsize_row = 4, main = "Top 100 Variable Genes")
dev.off()


#generating differential expression data
design <- model.matrix(~0+factor(c("flight", "flight", "flight", "ground", "ground", "ground")))
colnames(design) <- c("flightgroup", "groundgroup")
fit <- lmFit(expression_set, design)
contrast.matrix <- makeContrasts(flightgroup - groundgroup, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
res_table_sphinx <- topTable(fit2, coef="flightgroup - groundgroup", adjust='BH', num=Inf)

results <- decideTests(fit2, adjust.method = "BH", p.value = .05, lfc=0)
summary(results)
vennDiagram(results)
summary(results)
par(mfrow=c(1,2))
plotMD(fit2, coef=1, status=results[,"flightgroup - groundgroup"], values = c(-1,1), h1.col=c("blue", "red"), main = "Flight vs Ground")
volcanoplot(fit2, coef=1)
go <- goana(fit2, coef="flightgroup - groundgroup", species = "Hs")
topGO(go, n=10)

refseq_data <- read.csv("~/Desktop/SPHINX DATA/refseq_dataframe.csv")
refseq_dataframe <- data.frame(refseq_data)
head(refseq_dataframe)
for(i in 18735:nrow(refseq_dataframe)){
  entrezID = convert2EntrezID(IDs=refseq_dataframe[i,'ID'], orgAnn='org.Hs.eg.db', ID_type='refseq_id')
  if (!is.null(entrezID))
    refseq_dataframe[i,'ID'] <- entrezID
}

entrez_ID_dataframe = subset(refseq_dataframe, select= -c(X, X.1, X.2, X.3, X.4))
unique_ID_dataframe <- entrez_ID_dataframe[!duplicated(entrez_ID_dataframe$ID), ]
rownames(unique_ID_dataframe) <- unique_ID_dataframe$ID
unique_ID_dataframe$ID <- NULL
object <- new("ExpressionSet", exprs=as.matrix(unique_ID_dataframe))
expression_set <- object

groups <- c("flight", "flight", "flight", "ground", "ground", "ground")
design <- model.matrix(~0 + groups)
head(design)
fit <- lmFit(expression_set, design)
head(coef(fit))
contr <- makeContrasts(groupsflight - groupsground, levels = colnames(coef(fit)))
contr
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)






design <- model.matrix(~0+factor(c("flight", "flight", "flight", "ground", "ground", "ground")))
colnames(design) <- c("flightgroup", "groundgroup")






fit <- lmFit(expression_set, design)
contrast.matrix <- makeContrasts(flightgroup - groundgroup, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
res_table <- topTable(fit2, coef="flightgroup - groundgroup", adjust='BH', num=Inf)
results <- decideTests(fit2, adjust.method = "BH", p.value = .05, lfc=0)
gene_list <- res_table$logFC
names(gene_list) <- row.names(res_table)
gene_list = sort(gene_list, decreasing = TRUE)
hallmark <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, entrez_gene)
head(hallmark)
gse_result <- GSEA(gene_list, TERM2GENE=hallmark, pvalueCutoff = .05, pAdjustMethod="BH")
head(gse_result)
require(DOSE)
dotplot(gse_result, title = "SPHINX Space vs Ground GSEA Results", showCategory=10, split='.sign') + facet_grid(.~.sign)
dotplot(gse_result)
data_gsea <- data.frame(gse_result, num=Inf)

head(unique_ID_dataframe)
object <- new("ExpressionSet", exprs=as.matrix(unique_ID_dataframe))
expression_set <- object
write.table(expression_set, "gsea_dataframe.txt")
new_dataframe <- unique_ID_dataframe
new_dataframe$ID <- NULL
data_matrix <- as.matrix(new_dataframe)
pheatmap(new_dataframe)


library(ChIPpeakAnno)

athero_genes <- c('JAM2', 'JAM3', 'CCL2', 'CCL8', 'CCL7', 'CCL3', 'CCL5', 'CCL4', 'CRP', 'SELE', 'SELP', 'CD40LG', 'CCR2', 'CD81', 'IL6', 'IL6ST', 'IL6R', 'IFNG', 'TNF', 'ICAM1', 'ICAM2', 'PECAM1', 'CD99', 'PVR', 'CD47', 'VCAM1', 'MCP1', 'IL4', 'IL33', 'KLF2')
athero_IDs <- convert2EntrezID(IDs=athero_genes, orgAnn='org.Hs.eg.db', ID_type='gene_symbol')
athero_data <- unique_ID_dataframe[rownames(unique_ID_dataframe) %in% athero_IDs, ]
athero_symbols <- getSYMBOL(athero_data$ID, 'org.Hs.eg.db')
athero_data$ID <- athero_symbols
rownames(athero_data) <- athero_data$ID
athero_data$ID <- NULL
new_object <- new("ExpressionSet", exprs=as.matrix(athero_data))
athero_expression_set <- new_object
design <- model.matrix(~0+factor(c("flight", "flight", "flight", "ground", "ground", "ground")))
colnames(design) <- c("flightgroup", "groundgroup")
fit <- lmFit(athero_expression_set, design)
contrast.matrix <- makeContrasts(flightgroup - groundgroup, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
res_table <- topTable(fit2, coef="flightgroup - groundgroup", adjust='BH', num=Inf)
write.csv(res_table, "Differential Expression of Atherosclerosis Genes.csv")
results <- decideTests(fit2, adjust.method = "BH", p.value = .05, lfc=0)
summary(results)
vennDiagram(results)
summary(results)
par(mfrow=c(1,2))
plotMD(fit2, coef=1, status=results[,"flightgroup - groundgroup"], values = c(-1,1), h1.col=c("blue", "red"), main = "Atherosclerosis Gene Expression: Flight vs Ground")
volcanoplot(fit2, coef=1)
go <- goana(fit2, coef="flightgroup - groundgroup", species = "Hs")
topGO(go, n=10)

row.names(res_table_sphinx)[row.names(res_table_sphinx)=="NR_003125"] <- "SNORD14E"
pdf(file = '/Users/wang/Desktop/EC2_volcano.pdf',
    width =10, 
    height=5)
EnhancedVolcano(res_table_sphinx, title = "SPHINX Space vs Ground",legendLabels=c( 'adjPvalue','Not sig.','adjPvalue','adjPvalue'), lab=rownames(res_table_sphinx),  x='logFC', y='adj.P.Val', FCcutoff = 0, pCutoff = .05,drawConnectors = TRUE, max.overlaps = Inf)
dev.off()

pdf(file = '/Users/wang/Desktop/GSEA_sphinx.pdf',
    width =10, 
    height=5)
dotplot(gse_result, title = "SPHINX Space vs Ground GSEA Results", showCategory=10, split='.sign') + facet_grid(.~.sign)
dev.off()
