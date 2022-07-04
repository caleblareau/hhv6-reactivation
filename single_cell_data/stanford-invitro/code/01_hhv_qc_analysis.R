library(data.table)
library(dplyr)

lapply(c("34", "38", "97", "98"), function(ss){
  fread(paste0("../data/feature_counts/ALLO_Sample",ss,"_both.kb.txt.gz")) %>% 
    filter(V3 < 88) %>% unique() %>% dim()
})

lapply(c("34", "38", "97", "98"), function(ss){
  fread(paste0("../data/feature_counts/ALLO_Sample",ss,"_HHV6b.kb.txt.gz")) %>% 
    arrange(V1) %>% 
    filter(V3 != 120) %>% unique() %>% dim()
})