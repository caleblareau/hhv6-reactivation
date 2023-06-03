library(data.table)
library(dplyr)
library(Seurat)
source("../../stanford-invitro/code/00_functions.R")
library(ggplot2)
library(BuenColors)
library(viridis)

counts <- Read10X_h5("../data/tenx_h5/CD7_Day19_raw_feature_bc_matrix.h5")
# Create seurat object for easy QC
so <- CreateSeuratObject(counts = counts, project = "cart", min.cells = 1, min.features = 1)
so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^MT-")
so[["percent.ribo"]] <- PercentageFeatureSet(so, pattern = "^RPL|^RPS")

so_filt <- subset(so, percent.mt < 20 & percent.ribo < 40 & nCount_RNA > 1000 & nCount_RNA < 30000 & 
                    nFeature_RNA > 500 & nFeature_RNA < 5000)

# Append kite counts
hhv6_mat <- import_from_kite_matrix(unique(fread(paste0("../data/hhv6/fastq_CD7_Day19_S3_L001_R1_001_fastq_gz_HHV6b.quant.txt"))),
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
# CAGAGAGAGAATTGTG-1 and TTGGCAAGTAGCCTAT-1 are the CAR T products
# kinda random on the early - late axis
ggplot(so_filt@meta.data %>% filter(nHHV6 > 5), aes(x = late_pct*100, y = early_pct*100, color = immediate_early_pct)) +
  geom_point(size = 0.5) + 
  scale_color_viridis(limits = c(0,1)) + scale_y_continuous(limits = c(0,100)) + 
  scale_x_continuous(limits = c(0,100)) + 
  pretty_plot(fontsize = 8) + L_border() + labs(x = "%HHV6 UMIs - Late", y = "%HHV6 UMIs - Early") +
  theme(legend.position = "none") -> px
px
cowplot::ggsave2(px, file = paste0("../plots/",sampleID,"scatter_hhv6signature.pdf"), width = 2, height = 2)


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
