library(data.table)
library(dplyr)
library(Seurat)
source("00_functions.R")


# Import data
h5_file = paste0("../data/gexp/ALLO_D5720-Day15_rna_filtered_feature_bc_matrix.h5")
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
so_filt$HHV6 <- import_from_kite(unique(fread(paste0("../data/feature_counts/ALLO_D5720-Day15_HHV6b.kb.txt.gz"))),
                                 gsub("-1", "", colnames(so_filt)), logme = FALSE)
so_filt$pct_HHV6 <- so_filt$HHV6/(so_filt$HHV6 + so_filt$nCount_RNA)*100
so_filt <- so_filt %>%
  FindVariableFeatures() %>% NormalizeData() %>%
  ScaleData(vars.to.regress = c("nFeature_RNA", "percent.mt")) %>% RunPCA() %>%
  FindNeighbors() %>% FindClusters(resolution = 0.5) %>% RunUMAP(dims = 1:30)
p1 <- DimPlot(so_filt, reduction = "umap", label = TRUE) + NoLegend()
p2 <- FeaturePlot(so_filt, features = c( "pct_HHV6"), sort.cell = TRUE, max.cutoff = "q99")
p3 <- FeaturePlot(so_filt, features = c( "HHV6"), sort.cell = TRUE, max.cutoff = "q99")
p4 <- FeaturePlot(so_filt, features = c( "nFeature_RNA"), sort.cell = TRUE)
p5 <- FeaturePlot(so_filt, features = c( "percent.mt"), sort.cell = TRUE)
p6 <- FeaturePlot(so_filt, features = c( "percent.ribo"), sort.cell = TRUE)

p7 <- FeaturePlot(so_filt, features = c( "CD3E"),sort.cell = TRUE)
p8 <- FeaturePlot(so_filt, features = c( "S100A4"),sort.cell = TRUE)
p9 <- FeaturePlot(so_filt, features = c( "CD8A"), sort.cell = TRUE)
p10 <- FeaturePlot(so_filt, features = c( "CD4"),sort.cell = TRUE)
p11 <- FeaturePlot(so_filt, features = c( "MS4A1"),sort.cell = TRUE)
p12 <- FeaturePlot(so_filt, features = c( "GNLY"),sort.cell = TRUE)
p13 <- FeaturePlot(so_filt, features = c( "LTB"), sort.cell = TRUE)
p14 <- FeaturePlot(so_filt, features = c( "POLQ"),sort.cell = TRUE)
p15 <- FeaturePlot(so_filt, features = c( "HSPA5"),sort.cell = TRUE)
p16 <- FeaturePlot(so_filt, features = c( "BRCA1"),sort.cell = TRUE)
cowplot::ggsave2(cowplot::plot_grid(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, nrow = 4), 
                 width = 15, height = 12, filename = paste0("../plots/PBMC-D15-infection.png"))


cor_vec <- cor(t(data.matrix(so_filt@assays$RNA@data)), so_filt@meta.data$pct_HHV6)
df <- data.frame(
  cor = round(cor_vec[,1], 3),
  gene = rownames(cor_vec)
) %>% arrange(desc(cor))
df <- df[!is.na(df$cor),]

write.table(df, file = "../output/PBMC-Day15-correlation.tsv", 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


head(df, 15)
