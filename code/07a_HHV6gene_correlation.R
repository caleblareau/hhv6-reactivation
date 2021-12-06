library(data.table)
library(dplyr)
library(pheatmap)

# Import data
gene_vec <- fread("../data/reference/HHV6b_only.index.txt", header = FALSE)[[2]]
fread(paste0("../data/feature_counts/ALLO_D5720-Day15_HHV6b.kb.txt.gz")) %>% 
  mutate(idx = V3 +1) %>% 
  filter(idx <= 103) %>% unique() -> udf
udf$gene <-gene_vec[as.numeric(as.character(udf$idx))]
udf$barcode <- udf$V1
udf %>%
  group_by(barcode, gene) %>% summarize(count = n()) %>% 
  reshape2::dcast(barcode ~ gene, fill = 0, value.var = "count") -> count_df

# Make correlation plot
pheatmap(cor((data.matrix(count_df[,-1]))))
