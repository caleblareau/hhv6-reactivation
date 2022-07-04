library(data.table)
library(dplyr)
library(Seurat)
library(BuenColors)
library(DESeq2)
source("00_functions.R")

getPseudoBulks <- function(donor){
  
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
  so <- CreateSeuratObject(counts = counts, project = "cart", min.cells = 0, min.features = 0)
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
  
  so_filt2 <- subset(so_filt, (HHV6 >= 3) | (HHV6 == 0))
  so_filt2 <- so_filt2 %>%
    FindVariableFeatures() %>% NormalizeData() %>%
    ScaleData()
  
  df <- data.frame(
    high = Matrix::rowSums(so_filt2@assays$RNA@counts[,so_filt2$HHV6 >= 2]),
    low = Matrix::rowSums(so_filt2@assays$RNA@counts[,so_filt2$HHV6 == 0])
  )
  df
}

dd <- cbind(getPseudoBulks("34"),
            getPseudoBulks("98"))
d <- data.matrix(dd)
cpm <- round(t(t(d)/colSums(d)) * 1000000, 0)
gene <- rownames(dd); colnames(dd) <- c("SE34", "low34", "SE38", "low38")
colData <- data.frame(
  sample = colnames(dd),
  condition = c("high", "low", "high", "low")
)
dds <- DESeqDataSetFromMatrix(countData = dd,
                              colData = colData,
                              design = ~ condition)
keep <- rowSums(counts(dds)) >= 10

dds <- dds[keep,]
dds <- DESeq(dds)
res <- results(dds)
data.frame(res) %>% arrange((padj))

cpm[c("CD3D","CD3E", "TRAC","NPM3","CD4", "CD8A", "CD52", "CCL4", "GNLY"),]


library(EnhancedVolcano)
res$gene <- rownames(res)
EnhancedVolcano(res, x = "log2FoldChange", y = "padj", lab =  rownames(res), pCutoff = 0.0005, FCcutoff = 0.8)
