#cibersort signature matrix generation from single-cell reference

heart_gene_counts <- read.csv('~/Desktop/CiberSort/Heart-counts.csv')
heart_annots <- read.csv('~/Desktop/CiberSort/annotations_facs.csv')
only_heart_annots <- heart_annots[heart_annots$tissue == 'Heart',]
rownames(only_heart_annots) <- only_heart_annots$cell
cell_list <- colnames(heart_gene_counts)
named_list <- list('X')
for (i in cell_list){
  named_list <- c(named_list, only_heart_annots[i, 'cell_ontology_class'])
}
library(gtools)
named_list <- na.replace(named_list, 0)
named_list <- named_list[-2]
colnames(heart_gene_counts) <- c(named_list)
#get rid of empty rows

library(dplyr)
for (i in colnames(heart_gene_counts)){
  if (i == "X"){
    names(heart_gene_counts)[names(heart_gene_counts) == 'X'] <- "gene_symbol"
  }
  new_name <- only_heart_annots[i, 'cell_ontology_class']
  if (is.na(new_name)){
    names(heart_gene_counts)[names(heart_gene_counts) == i] <- '0'
  }
  names(heart_gene_counts)[names(heart_gene_counts) == i] <- new_name
}

heart_gene_counts <- heart_gene_counts[,!names(heart_gene_counts) %in% c('0')]
heart_gene_counts <- heart_gene_counts[,!names(heart_gene_counts) %in% c("")]
colnames(heart_gene_counts) <- gsub("[0-9]+","", colnames(heart_gene_counts))
colnames(heart_gene_counts) <- gsub("[.]", "", colnames(heart_gene_counts))
install.packages("writexl")
library("writexl")
write_xlsx(heart_gene_counts,"~/Desktop/CiberSort/signature_matrix.xlsx")
write.table(heart_gene_counts, file = "signature_matrix.txt", sep = "\t",
            row.names = FALSE, quote=FALSE)

library(readxl)
help <- read_excel("~/Desktop/signature_matrix_edit.xlsx")
colnames(help) <- gsub("[0-9]+","", colnames(help))
colnames(help) <- gsub("[.]", "", colnames(help))
write.table(help, file = "signature_matrix.txt", sep = "\t",
            row.names = FALSE, quote=FALSE)

cols_keep <- c("gene_symbol", "fibroblast", "endothelial cell", "leukocyte", "cardiac muscle cell")
small_sig <- help[, colnames(help) %in% cols_keep] 
write.table(small_sig, file = "signature_small.txt", sep = "\t",
            row.names = FALSE, quote=FALSE)
cols_keep_cm <- c("gene_symbol", "cardiac muscle cell")
small_matrix <- help[,colnames(help)%in% cols_keep_cm]
write.table(small_matrix, file = "signature_tiny.txt", sep = "\t",
            row.names = FALSE, quote=FALSE)
