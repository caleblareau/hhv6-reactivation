library(data.table)
library(dplyr)
library(pheatmap)
library(BuenColors)
# Import data
gene_vec <- fread("../data/reference/HHV6b_only.index.txt", header = FALSE)[[2]]
fread(paste0("../data/feature_counts/ALLO_D5720-Day15_HHV6b.kb.txt.gz")) %>% 
  filter(V3 != 120) %>%
  mutate(idx = V3 +1) %>% 
  filter(idx <= 103) %>% unique() -> udf
udf$gene <-gene_vec[as.numeric(as.character(udf$idx))]
udf$barcode <- udf$V1
udf %>%
  group_by(barcode, gene) %>% summarize(count = n()) %>% 
  reshape2::dcast(barcode ~ gene, fill = 0, value.var = "count") -> count_df

barcodes <- colnames(Seurat::Read10X_h5("../data/gexp/ALLO_D5720-Day15_rna_filtered_feature_bc_matrix.h5"))
count_df <- count_df %>% filter(barcode %in% gsub("-1", "", barcodes))

# Make correlation plot
pheatmap(cor((data.matrix(count_df[,-1]))))


cor(data.matrix(count_df[,-1]), rowSums(data.matrix(count_df[,-1])))

qplot(log10(rowSums(data.matrix(count_df[,-1]))))
highmat <- count_df[rowSums(data.matrix(count_df[,-1])) > 500,-1]
lowmat <- count_df[rowSums(data.matrix(count_df[,-1])) < 100,-1]
df <- data.frame(
  high_pct = colSums(highmat)/sum(highmat)*100,
  low_pct = colSums(lowmat)/sum(lowmat)*100,
  gene = colnames(highmat)
)
ggplot(df, aes(x = low_pct, y = high_pct, label = gene)) +
  geom_text() +
  scale_x_log10() + scale_y_log10() + 
  pretty_plot(fontsize = 16) + L_border() +
  labs(x = "% HHV6 RNA - Low HHV6 Expressing Cells", 
       y = "% HHV6 RNA - High HHV6 Expressing Cells")
