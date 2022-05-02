library(data.table)
library(dplyr)

combine_hhv6_clones <- function(ss){
  fread(paste0("../data/feature_counts/ALLO_",ss,"_HHV6b.kb.txt.gz")) %>% 
    filter(V3!= 120) %>% 
    mutate(idx = V3 +1) %>% 
    filter(idx <= 103) %>% unique() %>% group_by(V1) %>% summarize(count = n()) -> udf
  vec <- udf$count; names(vec) <- paste0(udf$V1, "-1")
  tcrdf <- fread(paste0("../data/tcrs/all/ALLO-TCRs_",ss,"_VDJ_all_contig_annotations.csv.gz"))
  tcrdf$hhv6_count <- vec[as.character(tcrdf$barcode)]
  tcrdf %>% arrange(desc(hhv6_count))
}

s34 <- combine_hhv6_clones("Sample61")
s34 %>% filter(chain == "TRA") %>% 
  arrange(desc(hhv6_count))

s34 %>% filter(hhv6_count > 10) %>%
  select(c(raw_clonotype_id, barcode)) %>%
  distinct() %>% 
  group_by(raw_clonotype_id) %>% 
  summarize(count = n()) %>%
  arrange(desc(count))
