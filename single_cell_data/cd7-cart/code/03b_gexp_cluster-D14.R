library(data.table)
library(Seurat)
library(dplyr)
library(BuenColors)

#  Define hhv6 expressing cells from broad script
  count_df2 <- fread("../data/hhv6/fastq_CD7_Day14_S2_L001_R1_001_fastq_gz_HHV6b.quant.txt") %>% 
  filter(!(V3 %in% c(6, 120))) %>%
  distinct() %>%
  group_by(V1) %>% summarize(tc = n()) %>%
  mutate(barcode = paste0(V1, "-1")) %>% arrange(desc(tc)) %>% data.frame()
count_df2

car_vec <- car_df <- fread("../data/car/CD7_Day14_CAR.quant.txt") %>%
  mutate(barcode = paste0(V1, "-1")) %>% pull(barcode) %>% table()
count_df2$car_copyN <- ifelse(is.na(car_vec[count_df2$barcode]), 0, car_vec[count_df2$barcode])

exp <- Read10X_h5("../data/tenx_h5/CD7_Day14_filtered_feature_bc_matrix.h5")
count_df2$called_cell <- count_df2$barcode %in% colnames(exp)
count_df2$valid_barcode <- count_df2$barcode %in% colnames(Read10X_h5("../data/tenx_h5/CD7_Day14_raw_feature_bc_matrix.h5"))

count_df2$totalUMIs <- colSums(Read10X_h5("../data/tenx_h5/CD7_Day14_raw_feature_bc_matrix.h5"))[count_df2$barcode]
head(count_df2)

# Look at snps
snp_df <- fread("../data/snp_assignment/Day14_mitoSNPassignments.tsv")
soup_df <- fread("../data/snp_assignment/CD7_Day14_clusters.tsv")
merge(snp_df, soup_df, by.x = "cell_id", by.y = "barcode") %>%
  group_by(assignment, assign, status) %>%
  summarize(count = n())

mdf <- merge(soup_df, count_df2, all.x = TRUE, by.x = "barcode", by.y = "barcode") %>% # 
  arrange(desc(tc))
mdf$CD7 <- exp["CD7", mdf$barcode]
mdf$CD7cp10k <- exp["CD7", mdf$barcode]/colSums(exp[, mdf$barcode])*10000
mdf$HHV6cp10k <- mdf$tc/(colSums(exp[, mdf$barcode])+ mdf$tc)*10000
mdf$HHV6log <- log1p(mdf$tc/(colSums(exp[, mdf$barcode])+ mdf$tc)*10000)

mdf %>% filter(tc > 0)

# SEurat things
mdf <- data.frame(mdf)
rownames(mdf) <- mdf$barcode
so <- CreateSeuratObject(counts = exp[, mdf$barcode], meta.data = mdf)
so <- NormalizeData(so) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% 
  FindNeighbors() %>%  FindClusters(resolution = 0.2) %>% RunUMAP(dims = 1:30)
DimPlot(so, group.by = c( "assignment", "seurat_clusters"), label = FALSE)
FeaturePlot(so, features = c("CD7", "GZMB", "TREM1", "CD3D", "CD4", "CD8A","GZMK"), sort = TRUE)
FeaturePlot(so, features = c("HHV6log", "MKI67"), sort = TRUE)


fm <- FindMarkers(so, "0", "1")
head(fm)
