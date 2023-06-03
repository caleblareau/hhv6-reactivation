library(data.table)
library(Seurat)
library(dplyr)
library(BuenColors)
library(ggbeeswarm)
library(viridis)

#  Define hhv6 expressing cells from broad script
experiment <- "CD7_Day14"
experimentX <- paste0(experiment, "_S2")
exp <- Read10X_h5(paste0("../data/tenx_h5/",experiment,"_raw_feature_bc_matrix.h5"))
count_df2 <- fread(paste0("../data/hhv6/fastq_",experimentX,"_L001_R1_001_fastq_gz_HHV6b.quant.txt")) %>% 
  filter(!(V3 %in% c(6, 120))) %>%
  distinct() %>%
  group_by(V1) %>% summarize(tc = n()) %>%
  mutate(barcode = paste0(V1, "-1")) %>% arrange(desc(tc)) %>% data.frame()
count_df2$count_UMI <- colSums(exp)[count_df2$barcode]
count_df2
count_df2 %>% filter(barcode %in% colnames(exp)) %>% dim()
