library(data.table)
library(Seurat)
library(dplyr)


#  Define hhv6 expressing cells from broad script
hits <- c("AACCGCGGTAACGCGA-3", "ACTGTCCGTTTGGGCC-3", "GCATACATCAGCTCGG-3", "GCGACCACACGTAAGG-3", "GGGCATCAGAACTCGG-3", "TGACTAGCACTTGGAT-2")

count_df2 <- rbind(
  fread("../data/hhv6-quant/1878509_JCC212_SJCAR19-09_WK2_PB_Gex_S10_L001_HHV6b_nocorrection.txt"),
  fread("../data/hhv6-quant/1878509_JCC212_SJCAR19-09_WK2_PB_Gex_S10_L002_HHV6b_nocorrection.txt"),
  fread("../data/hhv6-quant/1878509_JCC212_SJCAR19-09_WK2_PB_Gex_S10_L003_HHV6b_nocorrection.txt"),
  fread("../data/hhv6-quant/1878509_JCC212_SJCAR19-09_WK2_PB_Gex_S10_L004_HHV6b_nocorrection.txt")
) %>%
  filter(!(V3 %in% c(6, 120))) %>%
  distinct() %>%
  group_by(V1) %>% summarize(tc = n()) %>%
  mutate(barcode = paste0(V1, "-2")) %>% arrange(desc(tc)) %>% data.frame()

week2_positives <- count_df2 %>% pull(barcode)

exp_wk2 <- Read10X_h5("../data/gexp/JCC212_SJCAR19-09_Wk2_PB_short_raw_feature_bc_matrix.h5")
colnames(exp_wk2) <- gsub("-1", "-2", colnames(exp_wk2))
tail(exp)
rownames(exp)

genes <- c("JCC_SJCAR19short", "CD3D","CD3E", "CD4", "CD8A", "CD8B")
exp_wk2[genes,colnames(exp_wk2) %in% hits]

fread("../data/tcrs/1879720_JCC212_SJCAR19-09_Wk2_PB_TCR_all_contig_annotations.csv") %>%
  filter(gsub("-1", "-2", barcode) %in% hits) %>% group_by(barcode, chain) %>% summarize(sum(umis))

count_df2$totalUMIs <- colSums(exp_wk2)[week2_positives]
count_df2 %>%
  filter(barcode %in% hits)

#################

count_df3 <- rbind(
  fread("../data/hhv6-quant/1894736_JCC212_SJCAR19-09_Wk3_PB_Gex_S1_L001_HHV6b_nocorrection.txt"),
  fread("../data/hhv6-quant/1894736_JCC212_SJCAR19-09_Wk3_PB_Gex_S1_L002_HHV6b_nocorrection.txt")
) %>%
  filter(!(V3 %in% c(6, 120))) %>%
  distinct() %>%
  group_by(V1) %>% summarize(tc = n()) %>%
  mutate(barcode = paste0(V1, "-3")) %>% arrange(desc(tc)) %>% data.frame()
week3_positives <- count_df3 %>% pull(barcode)

exp_wk3 <- Read10X_h5("../data/gexp/JCC212_SJCAR19-09_Wk3_PB_short_raw_feature_bc_matrix.h5")
colnames(exp_wk3) <- gsub("-1", "-3", colnames(exp_wk3))

tail(exp)
rownames(exp)
genes <- c("JCC_SJCAR19short", "CD3D","CD3E", "CD4", "CD8A","CD8B" )
fdf_week3 <- data.frame(
  t(data.matrix(exp_wk3[genes,(colnames(exp_wk3) %in% hits)])
  ))

fread("../data/tcrs/1972991_JCC212_SJCAR19-09_Wk3_PB_TCR_all_contig_annotations.csv") %>%
  filter(gsub("-1", "-3", barcode) %in% week3_positives) %>% group_by(barcode, chain) %>% summarize(sum(umis))


count_df3$totalUMIs <- colSums(exp_wk3)[hits]

