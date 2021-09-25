library(data.table)
library(dplyr)
library(Seurat)
library(BuenColors)
source("00_functions.R")

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
  
  so_filt <- subset(so, percent.mt < 20 & percent.ribo < 40 & nCount_RNA > 1000 & nCount_RNA < 30000 & 
                      nFeature_RNA > 500 & nFeature_RNA < 5000)
  so_filt <- so 
  # Append kite counts
  so_filt$WPRE <- import_from_kite(unique(fread(paste0("../data/feature_counts/ALLO_Sample",donor,"_HHV6b.wpre.txt.gz"))),
                                   gsub("-1", "", colnames(so_filt)), logme = FALSE)
  so_filt$HHV6 <- import_from_kite(unique(fread(paste0("../data/feature_counts/ALLO_Sample",donor,"_HHV6b.kb.txt.gz"))),
                                   gsub("-1", "", colnames(so_filt)), logme = FALSE)
  
  # Examine the "jackpot" cells
  total_HHV6 <- sum(so_filt@meta.data$HHV6)
  total_RNAs <- sum(so_filt@meta.data$nCount_RNA)
  
  viral_abundance_df <- data.frame(
    null = sort(total_HHV6 * (so_filt@meta.data$nCount_RNA / total_RNAs), decreasing = TRUE),
    observed = sort(so_filt@meta.data$HHV6, decreasing = TRUE),
    rank = c(rep(1:dim(so_filt)[2],2))
  )
  
  # Count super expressors too
  n_SE <- sum(so_filt@meta.data$HHV6 >= 5)
  n_cells <- dim(so_filt)[2]
  title = paste0("S", donor, " n_super=", as.character(n_SE), "; n=", as.character(n_cells), "; HHV6max=", as.character(max(so_filt@meta.data$HHV6)))
  viral_abundance_df %>% reshape2::melt(id.vars = "rank") %>%
    ggplot(aes(x = rank, y = value , color = variable)) + 
    geom_point(size = 0.5) + pretty_plot() + scale_color_manual(values = c("grey", "firebrick")) + L_border() +
    labs(x = "Rank ordered cells", color = "Data", y = "# of HHV6 UMIs (log scaled counts)") +
    geom_hline(yintercept = 0, linetype = 2)+ scale_y_log10() + 
    geom_hline(yintercept = 4, linetype = 2, color = "darkblue") + ggtitle(title) -> plot1
  cowplot::ggsave2(plot1, file = paste0("../plots/nofilter_SuperExpressor_Sample", donor, ".png"), 
                   width = 4, height = 3)
  donor
}
lapply(c("34", "38", "97", "98"), make_superexpressor_plot)