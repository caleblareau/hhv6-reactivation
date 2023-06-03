library(data.table)
library(dplyr)
library(Seurat)
source("00_functions.R")


# Import data
process_numer_assoc <- function(ss = "Sample34-Day7"){
  h5_file = paste0("../data/gexp/ALLO_",ss,"_rna_filtered_feature_bc_matrix.h5")
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
  so_filt$HHV6 <- import_from_kite(unique(fread(paste0("../data/feature_counts/ALLO_",ss,"_HHV6b.kb.txt.gz"))),
                                   gsub("-1", "", colnames(so_filt)), logme = FALSE)
  so_filt$pct_HHV6 <- so_filt$HHV6/(so_filt$HHV6 + so_filt$nCount_RNA)*100
  so_filt <- so_filt %>%
    FindVariableFeatures() %>% NormalizeData() %>%
    ScaleData(vars.to.regress = c("nFeature_RNA", "percent.mt")) 
  
  # Compute correlations
  cor_vec <- cor(t(data.matrix(so_filt@assays$RNA@data)), so_filt@meta.data$pct_HHV6)
  set.seed(42)
  cor_vec_perm <- cor(t(data.matrix(so_filt@assays$RNA@data)), sample(so_filt@meta.data$pct_HHV6))
  
  df <- data.frame(
    cor = c(round(cor_vec[,1], 3)),
    gene = c(rownames(cor_vec))
  ) %>% arrange(desc(cor))
  df$rank <- 1:length(cor_vec[,1])
  permdf <- data.frame(
    cor = c(round(cor_vec_perm[,1], 3)),
    gene = rep("perm", dim(cor_vec_perm)[1])
  ) %>% arrange(desc(cor))
  permdf$rank <- 1:length(cor_vec_perm[,1])
  
  df <- df[!is.na(df$cor),]
  permdf <- permdf[!is.na(permdf$cor),]
  cutoff = quantile(sort(abs(permdf$cor)), c(0.99))
  data.frame(
    ss, 
    cutoff = as.numeric(cutoff), 
    n_hits = sum(abs(df$cor) > cutoff)
  )
}
process_numer_assoc("Sample34-Day7")
process_numer_assoc("Sample34-Day5")
process_numer_assoc("Sample61-Day7")

