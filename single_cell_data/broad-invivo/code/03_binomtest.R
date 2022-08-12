library(data.table)
library(dplyr)
library(Seurat)

fread("../data/fromNick/cell-annotation-metadata.txt", header = FALSE) %>%
  filter(V3 == "Axi-R-15" ) %>% dim()

fread("../data/fromNick/CART_global_obs.csv.gz", header = TRUE) %>%
  filter(barcode == "Axi-R-15" & timepoint_fine == "Infusion") %>% pull(timepoint_fine) %>% table()


# 1068 cells 
# we have 4 SEs
fisher.test(matrix(c(4,4731,0,1068), nrow = 2)) %>% str

# with resequencing, we got more data
#Day 0: 0/10171
# Day 7: 6/5026
# Day 14: 0/12039
# Day 21: 0/4717

fisher.test(matrix(c(4+6,4731+5026,0,1068+10171), nrow = 2)) 

10/(4731+5026)*100
fisher.test(matrix(c(4+6,4731+5026,0,12039), nrow = 2)) 
