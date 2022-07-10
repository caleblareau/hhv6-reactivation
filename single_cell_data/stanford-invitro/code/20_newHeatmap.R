library(data.table)
library(dplyr)
library(BuenColors)
library(Seurat)
library(ComplexHeatmap)

gene_idx <- fread("../../../hhv6-reference/HHV6b_only.index.txt", header = FALSE)

pp <- function(donor){
  kite_df <- fread(paste0("../data/feature_counts/ALLO_",donor,"_HHV6b.kb.txt.gz"))
  vec <- gene_idx[["V2"]]
  ncount <- kite_df %>% filter(V3 != "120" & V3 != "6") %>% group_by(V1,V3) %>%
    summarize(count = n()) %>% mutate(gene = vec[as.numeric(V3) + 1], sample = donor) %>%
    ungroup() %>%group_by(V1) %>% filter(sum(count) > 100)
  ncount[,c("sample", "V1", "gene", "count")]
  }
rbdf <- rbind(pp("Sample34"), pp("Sample97"), pp("Sample98"))
refmat <- rbdf[complete.cases(rbdf),] %>%
  reshape2::dcast(sample + V1 ~ gene, fill = 0, value.var = "count")
donor <- refmat[[1]]
hhv6_mat <- t(data.matrix(data.frame(refmat[,c(-1, -2)])))

# Set up gene annotations
anno_df <- fread("../../../hhv6-reference/HHV6b_expression_annotations.tsv", header = FALSE)
anno_df <- anno_df %>% 
  mutate(simple = case_when(
    anno_df$V2 == "late" ~ "late", 
    anno_df$V2 == "early" ~ "early", 
    TRUE ~ "aImmediate_early"
  ))

vec <- anno_df[[3]]; names(vec) <- anno_df[[1]]

anno_df <- data.frame(gene = rownames(hhv6_mat), anno = vec[as.character(rownames(hhv6_mat))])
adf_go <- anno_df[complete.cases(anno_df),] %>%
  arrange((anno))


# Summarize in heatmap annotation
ha <- HeatmapAnnotation(
  gene_annotation = adf_go$anno,
  col = list(
    gene_annotation = c("early" = "dodgerblue3",
                        "late" = "dodgerblue4",
                        "aImmediate_early" = "dodgerblue")
  ),
  gp = gpar(col = "black"),
  show_legend = FALSE, annotation_label = ""
)

# Now plot data
hm <- Heatmap(t(log1p(hhv6_mat[adf_go$gene,])),
              cluster_rows = TRUE, cluster_columns = FALSE,
        col = jdb_palette("solar_rojos"), top_annotation =  ha, 
        column_names_gp = grid::gpar(fontsize = 4),
        column_split  = adf_go$anno,
        row_split = donor,
        row_names_gp = grid::gpar(fontsize = 0), show_heatmap_legend = FALSE)

pdf(file="../plots/HHV6expressionHeatmap-viralgenes.pdf", width = 5, height = 2)  
par(cex.main=0.8,mar=c(1,1,1,1))
hm
dev.off()
