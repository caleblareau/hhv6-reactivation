library(data.table)
library(dplyr)
library(BuenColors)
library(Seurat)
library(ComplexHeatmap)

# Helper function to annotate genes
import_from_kite_matrix <- function(kite_df){
  gene_idx <- fread("../../../hhv6-reference/HHV6b_only.index.txt", header = FALSE)
  vec <- gene_idx[["V2"]]
  ncount <- kite_df %>% filter(V3 != "120" & V3 != "6") %>% group_by(V1,V3) %>%
    summarize(count = n()) %>% mutate(gene = vec[as.numeric(V3) + 1]) %>%
    data.table
  d <- (data.frame(dcast.data.table(ncount, V1 ~ gene, fill = 0, value.var = "count", fun.aggregate = sum)))
  rownames(d) <- d[[1]]
  mat <- data.matrix(d[,-c(1,2)])
  return(mat)
}

# Import counts matrix
import_hhv6_mat <- function(){
  hhv6_mat <- import_from_kite_matrix(
    rbind(unique(fread(paste0("../data/new-data/BBB7_HHV6b.quantFinal.txt")))
    ))
  rownames(hhv6_mat) <- paste0(rownames(hhv6_mat), "-1")
  t(hhv6_mat)
}

# Subset to cells expressing hhv6
hhv6_mat <- import_hhv6_mat()

anno_df <- fread("../../../hhv6-reference/HHV6b_expression_annotations.tsv", header = FALSE)
anno_df <- anno_df %>% 
  mutate(simple = case_when(
    anno_df$anno == "late" ~ "late", 
    anno_df$anno == "early" ~ "early", 
    TRUE ~ "aImmediate_early"
  ))

vec <- anno_df[[3]]; names(vec) <- anno_df[[1]]

anno_df <- data.frame(gene = rownames(hhv6_mat), anno = vec[as.character(rownames(hhv6_mat))])
adf_go <- anno_df[complete.cases(anno_df),] %>%
  arrange((anno))

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

hm <- Heatmap(t(log1p(hhv6_mat[adf_go$gene,])), cluster_rows = FALSE, cluster_columns = FALSE,
        col = jdb_palette("solar_rojos"), top_annotation =  ha, 
        column_names_gp = grid::gpar(fontsize = 4),
        column_split  = adf_go$anno,
        row_names_gp = grid::gpar(fontsize = 0), show_heatmap_legend = FALSE)

pdf(file="../plots/HHV6expressionHeatmap-day14invivo.pdf", width = 3, height = 2)  
par(cex.main=0.8,mar=c(1,1,1,1))
hm
dev.off()
