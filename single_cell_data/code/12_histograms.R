library(data.table)
library(dplyr)
library(Seurat)
library(BuenColors)
source("00_functions.R")

make_density_plot <- function(donor){
  
  # Import data
  h5_file = paste0("../data/gexp/ALLO_",donor,"_rna_filtered_feature_bc_matrix.h5")
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
  
  so_filt <- subset(so, percent.mt < 20 & percent.ribo < 40 & nCount_RNA > 1000 & nCount_RNA < 30000 & 
                      nFeature_RNA > 500 & nFeature_RNA < 5000)
  so_filt <- so 
  # Append kite counts
  so_filt$HHV6 <- import_from_kite(unique(fread(paste0("../data/feature_counts/ALLO_",donor,"_HHV6b.kb.txt.gz"))),
                                   gsub("-1", "", colnames(so_filt)), logme = FALSE)
  
  pq <- ggplot(so_filt@meta.data,aes(x = so_filt$HHV6 + 1)) +
    geom_histogram(color = "black", fill = "lightgrey", bins = 20) + 
    pretty_plot(fontsize = 7) + L_border() +
    scale_x_log10() + scale_y_log10(expand = c(0,0)) + 
    labs(x = "HHV6 expression (log10 scaled)", y = "Number of cells")
  cowplot::ggsave2(pq, file = paste0("../plots/histogram_", donor, ".pdf"), width = 2.8, height = 1.4)
  # Now plot the distribution
  print(table(so_filt$HHV6 > 10))
  print(sum(so_filt$HHV6 [so_filt$HHV6 > 10])/sum(so_filt$HHV6 ))
  
}

lapply(c("Sample61-Day7", "Sample34-Day5", "Sample34-Day7", "GMP6-Day6", "D5720-Day15"), make_density_plot)
lapply(c("Sample34", "Sample38", "Sample97", "Sample98"), make_density_plot)