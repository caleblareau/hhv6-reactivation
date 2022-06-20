library(SingleCellExperiment)
library(monocle)
library(Seurat)
library(org.Hs.eg.db)
library(garnett)

process_me <- function(x){
  ri <- Read10X_h5(paste0("../data/rnaseq/Z",x,"_filtered_feature_bc_matrix.h5"))
  
  pdata <- data.frame(
    cell = colnames(ri) 
  )
  rownames(pdata) <- pdata[[1]]
  fdata <- data.frame(
    gene = rownames(ri)
  )
  rownames(fdata) <- fdata[[1]]
  
  # create a new CDS object
  pd <- new("AnnotatedDataFrame", data = pdata)
  fd <- new("AnnotatedDataFrame", data = fdata)
  sc <- newCellDataSet(ri, 
                       phenoData = pd,
                       featureData = fd)
  sc <- estimateSizeFactors(sc)
  
  
  t_classifier <- train_cell_classifier(cds = sc,
                                        marker_file =  "../data/hsTcell.txt",
                                        db=org.Hs.eg.db,
                                        cds_gene_id_type = "SYMBOL",
                                        num_unknown = 50,
                                        marker_file_gene_id_type = "SYMBOL")
  
  t_cds <- classify_cells(sc, t_classifier,
                          db = org.Hs.eg.db,
                          cluster_extend = TRUE, rank_prob_ratio = 1.5,
                          cds_gene_id_type = "SYMBOL")
  df <- data.frame(pData(t_cds)[,c(1,4)])
  df$cd8markers <- colSums(ri[c("CD8A", "CD8B"),])
  df$cd4markers <- colSums(ri[c("CD4", "FOXP3"),])
  write.table(df, file = paste0("../output/garnett_classifier_Z", x,".tsv"), 
              row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
}

lapply(as.character(1:6), process_me)


