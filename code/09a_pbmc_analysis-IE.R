library(data.table)
library(dplyr)
library(BuenColors)
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
hhv6_mat <- import_from_kite_matrix(unique(fread(paste0("../data/feature_counts/ALLO_D5720-Day15_HHV6b.kb.txt.gz"))),
                                 gsub("-1", "", colnames(so_filt)))
hhv6_annotation <- fread("../data/reference/HHV6b_expression_annotations.tsv", header = FALSE)
process_hhv6_annotations_to_pct <- function(mat, annotation_df){
  rs <- rowSums(mat)
  possible <- unique(annotation_df[[2]])
  mat_pct <- sapply(possible, function(x){
    genes_keep <- annotation_df[x == annotation_df[[2]],][[1]]
    rowSums(mat[,colnames(mat) %in% genes_keep])/rs
  })
}
hhv6_pct_mat <- process_hhv6_annotations_to_pct(hhv6_mat, hhv6_annotation)

# Append to seurat transcriptome
so_filt$immediate_early_pct <- hhv6_pct_mat[,"immediate_early"]
so_filt$intermediate_early_early_pct <- hhv6_pct_mat[,"intermediate_early_early"]
so_filt$early_pct <- hhv6_pct_mat[,"early"]
so_filt$late_pct <- hhv6_pct_mat[,"late"]
so_filt <- so_filt %>%
  FindVariableFeatures() %>% NormalizeData() %>%
  ScaleData(vars.to.regress = c("nFeature_RNA", "percent.mt")) %>% RunPCA() %>%
  FindNeighbors() %>% FindClusters(resolution = 0.5) %>% RunUMAP(dims = 1:30)

cor_mat <- cor(t(data.matrix(so_filt@assays$RNA@data)), hhv6_pct_mat)
df <- data.frame(
  cor = round(cor_mat, 3),
  gene = rownames(cor_mat)
)
df <- df[!(is.na(df[,1]) | is.na(df[,2]) | is.na(df[,3]) | is.na(df[,4])), ]

df <- df %>%
  mutate(label_rep = case_when(
    cor.late < -0.15 ~ gene,
    cor.late > 0.39  ~ gene,
    TRUE ~ ""
  ))

library(ggrepel)
ggplot(df, aes(x = cor.late, y = cor.early, label = label_rep)) +
  geom_point(aes( color = label_rep != "")) + 
  geom_text_repel(max.overlaps = Inf, min.segment.length = 0, seed = 42, box.padding = 0.5) + 
  scale_color_manual(values = c("black", "red")) + 
  scale_x_continuous(limits = c(-0.25, 0.5)) +
  scale_y_continuous(limits = c(-0.5, 0.3)) +
  
  pretty_plot(fontsize = 20) + L_border()+  theme(legend.position = "none") +
  labs(x = "Cor. w/ HHV6 late genes", 
       y = "Cor. w/ HHV6 early genes")



df %>% arrange(desc(cor.immediate_early)) %>% tail(20)
