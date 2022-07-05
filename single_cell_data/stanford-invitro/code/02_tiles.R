library(ggplot2)
library(data.table)
library(dplyr)
library(Seurat)
source("00_functions.R")

process_tiles <- function(donor){
  
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
  
  so_filt <- subset(so, percent.mt < 20 & percent.ribo < 40 & nCount_RNA > 1000 & nCount_RNA < 30000 & 
                      nFeature_RNA > 500 & nFeature_RNA < 5000)
  
  # Append kite counts
  so_filt$WPRE <- import_from_kite(unique(fread(paste0("../data/feature_counts/ALLO_Sample",donor,"_HHV6b.wpre.txt.gz"))),
                                   gsub("-1", "", colnames(so_filt)))
  so_filt$HHV6 <- import_from_kite(unique(fread(paste0("../data/feature_counts/ALLO_Sample",donor,"_HHV6b.kb.txt.gz"))),
                                   gsub("-1", "", colnames(so_filt)))
  so_filt <- so_filt %>%
    FindVariableFeatures() %>% NormalizeData() %>%
    ScaleData(vars.to.regress = c("nFeature_RNA", "percent.mt")) %>% RunPCA() %>%
    FindNeighbors() %>% FindClusters(resolution = 0.5) %>% RunUMAP(dims = 1:20)
  
  ####
  p1 <- DimPlot(so_filt, reduction = "umap", label = TRUE) + NoLegend()
  p2 <- FeaturePlot(so_filt, features = c( "WPRE"),sort.cell = TRUE, max.cutoff = "q99")+ theme_void()
  p3 <- FeaturePlot(so_filt, features = c( "HHV6"), sort.cell = TRUE, max.cutoff = "q99")+ theme_void()
  p4 <- FeaturePlot(so_filt, features = c( "nFeature_RNA"), sort.cell = TRUE)+ theme_void()
  p5 <- FeaturePlot(so_filt, features = c( "percent.mt"), sort.cell = TRUE)+ theme_void()
  p6 <- FeaturePlot(so_filt, features = c( "percent.ribo"), sort.cell = TRUE)+ theme_void()
  
  p7 <- FeaturePlot(so_filt, features = c( "CD3E"),sort.cell = TRUE)+ theme_void()
  p8 <- FeaturePlot(so_filt, features = c( "TNFRSF4"),sort.cell = TRUE)+ theme_void()
  p9 <- FeaturePlot(so_filt, features = c( "CD8A"), sort.cell = TRUE)+ theme_void()
  p10 <- FeaturePlot(so_filt, features = c( "CD4"),sort.cell = TRUE)+ theme_void()
  p11 <- FeaturePlot(so_filt, features = c( "MKI67"),sort.cell = TRUE)+ theme_void()
  p12 <- FeaturePlot(so_filt, features = c( "LEF1"),sort.cell = TRUE)+ theme_void()
  p13 <- FeaturePlot(so_filt, features = c( "CCR7"), sort.cell = TRUE)+ theme_void()
  p14 <- FeaturePlot(so_filt, features = c( "S100A4"),sort.cell = TRUE)+ theme_void()
  p15 <- FeaturePlot(so_filt, features = c( "LTA"),sort.cell = TRUE)+ theme_void()
  p16 <- FeaturePlot(so_filt, features = c( "LTB"),sort.cell = TRUE) + theme_void()
  
  cowplot::ggsave2(cowplot::plot_grid(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, nrow = 4), 
                   width = 15, height = 12, filename = paste0("../plots/Sample", donor, ".png"))
  
  cowplot::ggsave2(cowplot::plot_grid(p9, p10, p3, p11, nrow = 1), 
                   width = 12, height = 2.5, filename = paste0("../plots/simpleSample", donor, ".png"), dpi = 400)
}

process_tiles("34-Day5")
process_tiles("34-Day7")
process_tiles("61-Day7")

process_tiles("34")
process_tiles("38")
process_tiles("97")
process_tiles("98")
