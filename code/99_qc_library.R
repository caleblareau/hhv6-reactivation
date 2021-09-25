library(data.table)
library(dplyr)
library(Seurat)

process_qc <- function(library_go){
  print(library_go)
  h5_file = paste0("../data/gexp/",library_go,"_rna_filtered_feature_bc_matrix.h5")
  input <- Read10X_h5(h5_file)
  if(is.list(input)){
    counts <- input[[1]]
  } else {
    counts <- input
  }
  og_counts <- dim(counts)[2]
  # Create seurat object for easy QC
  so <- CreateSeuratObject(counts = counts, project = "cart", min.cells = 1, min.features = 1)
  so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^MT-")
  so[["percent.ribo"]] <- PercentageFeatureSet(so, pattern = "^RPL|^RPS")
  summary(so$nCount_RNA)
  summary(so$nFeature_RNA)
  summary(so$percent.ribo)
  summary(so$percent.mt)
  so_filt <- subset(so, percent.mt < 20 & percent.ribo < 40 & nCount_RNA > 1000 & nCount_RNA < 30000 & 
                      nFeature_RNA > 500 & nFeature_RNA < 5000)
  data.frame(
    library_go,
    cell_n = dim(so_filt)[2],
    cell_prop = round(dim(so_filt)[2]/og_counts,2),
    mean_ribo = round(mean(so_filt$percent.ribo),2),
    mean_MT = round(mean(so_filt$percent.mt),2),
    mean_UMI = round(mean(so_filt$nCount_RNA),1),
    mean_Genes = round(mean(so_filt$nFeature_RNA),1)
  )
}

qc_df <- rbind(process_qc("tenxPublic"),
      process_qc("ALLO_Sample34"),
      process_qc("ALLO_Sample38"),
      process_qc("ALLO_Sample97"),
      process_qc("ALLO_Sample98"))
write.table(qc_df, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
