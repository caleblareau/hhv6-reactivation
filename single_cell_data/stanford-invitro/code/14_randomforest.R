library(Seurat)
library(dplyr)
library(data.table)
library(randomForest)
library(Matrix)
source("00_functions.R")

lib1 <- "ALLO_Sample98"
process_rf <- function(lib1){
  
  # Import and munge counts
  mat <- Read10X_h5(paste0("../data/gexp/",lib1,"_rna_filtered_feature_bc_matrix.h5"))
  
  nHHV6 <- import_from_kite(unique(fread(paste0("../data/feature_counts/",lib1,"_HHV6b.kb.txt.gz"))),
                            gsub("-1", "", colnames(mat)), logme = FALSE)
  se <- mat[,as.numeric(nHHV6 )> 9] 
  n_cells_neg <- 3000
  set.seed(42)
  neg <- mat[,sample(which(as.numeric(nHHV6 )==0), n_cells_neg)] 
  
  # Set up features
  n_features <- 3000
  so <- CreateSeuratObject(
    cbind(se, neg)
  ) %>% NormalizeData() %>% FindVariableFeatures(selection.method = "vst", nfeatures = n_features)
  rfo <- randomForest::randomForest(
    t(data.matrix(so@assays$RNA@data[so@assays$RNA@var.features,])),
    c(rep(0, dim(se)[2]), rep(1, n_cells_neg))
  )
  
  # Create summary statistics
  cpm_norm <- function(mat){
    x<- rowSums(mat)
    round(x/sum(x)*1000000,1)
  }
  genes <- rownames(rfo$importance)
  data.frame(
    gene = genes,
    importance = unname(rfo$importance[,1]),
    cpm_se = cpm_norm(se)[genes],
    cpm_neg = cpm_norm(neg)[genes]
  ) %>% arrange(desc(importance)) %>% filter(importance > 0.1) -> odf
  write.table(odf, file = paste0("../output/rf_output/",lib1,".tsv"), 
              sep = "\t", quote= FALSE, row.names = FALSE, col.names = TRUE)
  lib1
  
}

#process_rf("ALLO_Sample98")
process_rf("ALLO_Sample34")

