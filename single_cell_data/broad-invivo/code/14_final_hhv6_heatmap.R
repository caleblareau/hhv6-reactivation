library(data.table)
library(dplyr)
library(BuenColors)
library(Seurat)
library(ComplexHeatmap)

gene_idx <- fread("../../../hhv6-reference/HHV6b_only.index.txt", header = FALSE)
vec <- gene_idx[["V2"]]

# Helper function to annotate genes
import_from_kite_matrix <- function(kite_df){
  ncount <- kite_df %>% filter(V3 != "120" & V3 != "6") %>% group_by(V1,V3) %>%
    summarize(count = n()) %>% mutate(gene = vec[as.numeric(V3) + 1]) %>%
    mutate(barcode = paste0(V1, "-1")) %>%data.frame()
    return(ncount)
}

# Import counts matrix
import_hhv6_mat <- function(){
  hhv6_mat <- import_from_kite_matrix(
    rbind(unique(fread(paste0("../data/E1_1_HHV6b.quantFinal.txt"))),
          unique(fread(paste0("../data/E1_2_HHV6b.quantFinal.txt"))),
          unique(fread(paste0("../data/F1_HHV6b.quantFinal.txt")))
    ))
  rownames(hhv6_mat) <- paste0(rownames(hhv6_mat), "-1")
  return(hhv6_mat)
}

# Subset to cells expressing hhv6
edf <- read.table("../output/HHV6-positive-cells.tsv")
hhv6_df_dfci <- import_hhv6_mat()

# Now add st jude
count_df2 <- rbind(
  fread("../../stjude-invivo/data/hhv6-quant/1878509_JCC212_SJCAR19-09_WK2_PB_Gex_S10_L001_HHV6b_nocorrection.txt"),
  fread("../../stjude-invivo/data/hhv6-quant/1878509_JCC212_SJCAR19-09_WK2_PB_Gex_S10_L002_HHV6b_nocorrection.txt"),
  fread("../../stjude-invivo/data/hhv6-quant/1878509_JCC212_SJCAR19-09_WK2_PB_Gex_S10_L003_HHV6b_nocorrection.txt"),
  fread("../../stjude-invivo/data/hhv6-quant/1878509_JCC212_SJCAR19-09_WK2_PB_Gex_S10_L004_HHV6b_nocorrection.txt")
) %>%
  filter(!(V3 %in% c(6, 120))) %>%
  distinct() %>% group_by(V1,V3) %>%
  summarize(count = n()) %>% mutate(gene = vec[as.numeric(V3) + 1]) %>%
  mutate(barcode = paste0(V1, "-2")) %>%data.frame()

count_df3 <- rbind(
  fread("../../stjude-invivo/data/hhv6-quant/1894736_JCC212_SJCAR19-09_Wk3_PB_Gex_S1_L001_HHV6b_nocorrection.txt"),
  fread("../../stjude-invivo/data/hhv6-quant/1894736_JCC212_SJCAR19-09_Wk3_PB_Gex_S1_L002_HHV6b_nocorrection.txt")
) %>%
  filter(!(V3 %in% c(6, 120))) %>%
  distinct() %>% group_by(V1,V3) %>%
  summarize(count = n()) %>% mutate(gene = vec[as.numeric(V3) + 1]) %>%
  mutate(barcode = paste0(V1, "-3")) %>% data.frame()


count_dfreseq <- read.table("../data/new-data/features/BBB7_HHV6b.quantFinal.txt") %>%
  filter(!(V3 %in% c(6, 120))) %>%
  distinct() %>% group_by(V1,V3) %>%
  summarize(count = n()) %>% mutate(gene = vec[as.numeric(V3) + 1]) %>%
  mutate(barcode = paste0(V1, "-4")) %>% data.frame()
reseq_positive_cells <- gsub("-1", "-4", fread("../output/HHV6-positive-cells-reseq.tsv") %>% data.frame() %>% pull(V1))

# import st jude gexp data
exp_wk2 <- Read10X_h5("../../stjude-invivo/data/gexp/JCC212_SJCAR19-09_Wk2_PB_short_raw_feature_bc_matrix.h5")
exp_wk3 <- Read10X_h5("../../stjude-invivo/data/gexp/JCC212_SJCAR19-09_Wk3_PB_short_raw_feature_bc_matrix.h5")
keep_week2 <- gsub("-1", "-2", names(colSums(exp_wk2))[(colSums(exp_wk2)) > 1000])
keep_week3 <- gsub("-1", "-3", names(colSums(exp_wk3))[(colSums(exp_wk3)) > 1000])

d <- (data.frame(dcast.data.table(data.table(rbind(
  hhv6_df_dfci %>%filter(barcode %in% rownames(edf)),
  count_dfreseq %>% filter(barcode %in% reseq_positive_cells),
  count_df2%>%filter(barcode %in% keep_week2),
  count_df3%>%filter(barcode %in% keep_week3) )), barcode ~ gene, fill = 0, value.var = "count", fun.aggregate = sum)))
rownames(d) <- d[[1]]
mat <- data.matrix(d[,-c(1,2)])
substr(rownames(mat), 18, 18)


# import ordering of cells from supp table
supp_table_df <- fread("../output/supp-invivo-table-reformatted.txt")

#########
# Import HHV-6 annotations'
####

anno_df <- fread("../../../hhv6-reference/HHV6b_expression_annotations.tsv", header = FALSE)
anno_df <- anno_df %>% 
  mutate(simple = case_when(
    anno_df$V2 == "late" ~ "late", 
    anno_df$V2 == "early" ~ "early", 
    TRUE ~ "aImmediate_early"
  ))

vec2 <- anno_df[[3]]; names(vec2) <- anno_df[[1]]

#-----
# Add st jude data
#-----

# Annotate for a combined heatmap
hhv6_mat <- t(mat)
anno_df <- data.frame(gene = rownames(hhv6_mat), anno = vec2[as.character(rownames(hhv6_mat))])
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

hhv6_mat[adf_go$gene,supp_table_df$barcode]+1

hm <- Heatmap(t(log1p(hhv6_mat[adf_go$gene,supp_table_df$barcode]+1)), cluster_rows = FALSE, cluster_columns = FALSE,
        col = jdb_palette("solar_rojos"), top_annotation =  ha, 
        column_names_gp = grid::gpar(fontsize = 4),
        column_split  = adf_go$anno,
        row_names_gp = grid::gpar(fontsize = 0), show_heatmap_legend = FALSE)

pdf(file="../plots/HHV6expressionHeatmap-wstj-wreseq-revision.pdf", width = 3, height = 3)  
par(cex.main=0.8,mar=c(1,1,1,1))
hm
dev.off()


### Now make the panel with gene expression / etc.

library(BuenColors)
library(forcats)

plot_df <- data.frame(
  barcode = supp_table_df$barcode,
  hhv6 = "hhv6",
  car = as.character(supp_table_df$n_car!= 0),
  tcr = ifelse(supp_table_df$channel == "Day7reseq" , "NA", (supp_table_df$nTCRA + supp_table_df$nTCRB) > 0),
  cd3 = as.character((supp_table_df$CD3D + supp_table_df$CD3E) > 0)
)

p1 <- reshape2::melt(plot_df, id.vars = "barcode") %>%
  mutate(name = factor(as.character(barcode), levels = rev(as.character( supp_table_df$barcode)))) %>%
  ggplot(aes(x = variable, y = name, fill = value)) + 
  geom_tile(color = "black") +
  pretty_plot(fontsize = 8) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) + labs(x = "", y = "") +
  scale_fill_manual(values = c("FALSE" = "white", "hhv6" = "dodgerblue3", "NA" = "darkgrey", "TRUE" = "darkgreen")) +
  theme(legend.position = "none")
p1
cowplot::ggsave2(p1, file = "../plots/grid.pdf", width = 2, height = 1.7)


