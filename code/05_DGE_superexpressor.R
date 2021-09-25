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
  
  # Append kite counts
  so$WPRE <- import_from_kite(unique(fread(paste0("../data/feature_counts/ALLO_Sample",donor,"_HHV6b.wpre.txt.gz"))),
                                   gsub("-1", "", colnames(so)), logme = FALSE)
  so$HHV6 <- import_from_kite(unique(fread(paste0("../data/feature_counts/ALLO_Sample",donor,"_HHV6b.kb.txt.gz"))),
                                   gsub("-1", "", colnames(so)), logme = FALSE)
  
  so_filt <- subset(so, (percent.mt < 20 & percent.ribo < 40 & nCount_RNA > 1000 & nCount_RNA < 30000 & 
                      nFeature_RNA > 500 & nFeature_RNA < 5000) | HHV6 >= 5)

  so_filt$super_expressor <- ifelse(so_filt$HHV6 >= 2, "SE", "other")
  
  so_filt2 <- subset(so_filt, (HHV6 >= 1) | (HHV6 == 0))
  so_filt2 <- so_filt2 %>%
    FindVariableFeatures() %>% NormalizeData() %>%
    ScaleData()
  so_filt2
}

fm98 <- FindMarkers(make_superexpressor_plot("98"), min.pct = 0.01, logfc.threshold = 0.01,
                    group.by = "super_expressor",
                    ident.1 = "SE", ident.2 = "other", max.cells.per.ident =1000)

fm34 <- FindMarkers(make_superexpressor_plot("34"), min.pct = 0.01, logfc.threshold = 0.01,
                    group.by = "super_expressor",
                    ident.1 = "SE", ident.2 = "other", max.cells.per.ident =1000)

fm97 <- FindMarkers(make_superexpressor_plot("97"), min.pct = 0.01, logfc.threshold = 0.01,
                    group.by = "super_expressor",
                    ident.1 = "SE", ident.2 = "other", max.cells.per.ident =1000)
head(fm97)
head(fm98)
head(fm34)

fm97[c("LTA", "LTB", "IFITM1", "TMEM97"),]
fm34$stat_34 <- -log10(fm34$p_val) * fm34$avg_log2FC
fm98$stat_98 <- -log10(fm98$p_val) * fm98$avg_log2FC
mdf <- merge(fm34, fm98, by = "row.names")

p1 <- ggplot(mdf, aes(x = stat_98, y = stat_34, label = Row.names)) +
  geom_text(size = 2) +
  pretty_plot(fontsize = 8) + L_border() +
  ggtitle("Statistic = -log10pval * logFC HHV6+/-") +
  labs(x = "Donor 98 Statistic", y = "Donor 34 Statistic")
cowplot::ggsave2(p1, file = "../plots/foldchange_statistic.png", width = 4, height = 4)
mdf %>% filter(stat_98 > 4 | stat_34 > 4)
