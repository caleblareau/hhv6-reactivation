library(data.table)
library(dplyr)
library(Seurat)
library(BuenColors)
source("00_functions.R")
library(Matrix)

dim(Read10X_h5("../data/gexp/ALLO_Sample61-Day7_rna_filtered_feature_bc_matrix.h5"))
    
make_superexpressor_plot <- function(donor){
  
  # Import data
  h5_file = paste0("../data/gexp/ALLO_Sample",donor,"_rna_filtered_feature_bc_matrix.h5")
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
  
  # Append kite counts
  so$WPRE <- import_from_kite(unique(fread(paste0("../data/feature_counts/ALLO_Sample",donor,"_HHV6b.wpre.txt.gz"))),
                                   gsub("-1", "", colnames(so)), logme = FALSE)
  so$HHV6 <- import_from_kite(unique(fread(paste0("../data/feature_counts/ALLO_Sample",donor,"_HHV6b.kb.txt.gz"))),
                                   gsub("-1", "", colnames(so)), logme = FALSE)
  
  so_filt <- subset(so,  HHV6 >= 1)
  data.frame(
    HHV6 = so_filt$HHV6,
    t(so_filt@assays$RNA@counts[c("CD3E", "CD4", "CD8A", "CD8B", "CCR7","IL7R", "CREM"),])
  ) %>% arrange(desc(HHV6))
}

make_superexpressor_plot("38")
