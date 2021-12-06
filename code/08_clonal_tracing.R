library(data.table)
library(dplyr)

combine_hhv6_clones <- function(ss){
  fread(paste0("../data/feature_counts/ALLO_",ss,"_HHV6b.kb.txt.gz")) %>% 
    mutate(idx = V3 +1) %>% 
    filter(idx <= 103) %>% unique() %>% group_by(V1) %>% summarize(count = n()) -> udf
  vec <- udf$count; names(vec) <- paste0(udf$V1, "-1")
  tcrdf <- fread(paste0("../data/tcrs/all/ALLO-TCRs_",ss,"_VDJ_all_contig_annotations.csv.gz"))
  tcrdf$count <- vec[as.character(tcrdf$barcode)]
  tcrdf %>% arrange(desc(count))
}

s34 <- combine_hhv6_clones("Sample34")
write.table(s34, file = paste0("../output/",ss, "_tcr_hhv6.tsv"),
            sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)