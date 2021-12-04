library(data.table)
library(dplyr)

gene_vec <- fread("../data/reference/HHV6b_only.index.txt", header = FALSE)[[2]]
lapply(c("34", "38", "97", "98"), function(ss){
  fread(paste0("../data/feature_counts/ALLO_Sample",ss,"_HHV6b.kb.txt.gz")) %>% 
    mutate(idx = V3 +1) %>% 
    filter(idx <= 103) %>% unique() -> udf
  udf$gene <-gene_vec[as.numeric(as.character(udf$idx))]
  udf$barcode <- udf$V1
  udf %>%
    group_by(barcode, gene) %>% summarize(count = n()) %>% 
    reshape2::dcast(barcode ~ gene, fill = 0, value.var = "count") -> count_df
  write.table(count_df, file = paste0("../countdf/HHV6b_counts_ALLO", ss, ".tsv"), 
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
})

