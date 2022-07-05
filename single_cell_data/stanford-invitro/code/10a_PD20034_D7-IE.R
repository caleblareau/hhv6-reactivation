library(data.table)
library(dplyr)
library(Seurat)
source("00_functions.R")

# Import data
sampleID <- "ALLO_Sample34-Day7"
sampleID <- "ALLO_Sample61-Day7"


h5_file = paste0("../data/gexp/",sampleID,"_rna_filtered_feature_bc_matrix.h5")
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

# Append kite counts
hhv6_mat <- import_from_kite_matrix(unique(fread(paste0("../data/feature_counts/",sampleID,"_HHV6b.kb.txt.gz"))),
                                    gsub("-1", "", colnames(so_filt)))
hhv6_annotation <- fread("../../../hhv6-reference/HHV6b_expression_annotations.tsv", header = FALSE)
process_hhv6_annotations_to_pct <- function(mat, annotation_df){
  mat <-mat[,colnames(mat) %in% annotation_df[["V1"]]]
  rs <- rowSums(mat) 
  possible <- unique(annotation_df[[2]])
  mat_pct <- sapply(possible, function(x){
    genes_keep <- annotation_df[x == annotation_df[[2]],][[1]]
    rowSums(mat[,colnames(mat) %in% genes_keep])/rs
  })
  mat_pct
}
hhv6_pct_mat <- process_hhv6_annotations_to_pct(hhv6_mat, hhv6_annotation)

# Append to seurat transcriptome
so_filt$immediate_early_pct <- hhv6_pct_mat[,"immediate_early"]
so_filt$intermediate_early_early_pct <- hhv6_pct_mat[,"intermediate_early_early"]

so_filt$early_pct <- hhv6_pct_mat[,"early"]
so_filt$late_pct <- hhv6_pct_mat[,"late"]
so_filt$nHHV6 <- rowSums(hhv6_mat)



ggplot(so_filt@meta.data %>% filter(nHHV6 > 10), aes(x = nHHV6, y = immediate_early_pct)) +
  geom_point() + scale_x_log10()


# Do more seurat things
so_filt <- so_filt %>%
  FindVariableFeatures() %>% NormalizeData() 

ct <- cor.test(so_filt@assays$RNA@data["TNFRSF4",],hhv6_pct_mat[,"immediate_early"])

cor_mat <- cor(t(data.matrix(so_filt@assays$RNA@data)), hhv6_pct_mat, use = "pairwise.complete")
df <- data.frame(
  cor = round(cor_mat, 3),
  gene = rownames(cor_mat)
)
cbind(data.frame(sampleID,
  pvalue_ie = ct$p.value), (cor_mat["TNFRSF4",,drop = FALSE]))
