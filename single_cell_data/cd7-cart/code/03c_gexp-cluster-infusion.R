library(data.table)
library(Seurat)
library(dplyr)
library(BuenColors)

#  Define hhv6 expressing cells from broad script
  count_df2 <- fread("../data/hhv6/fastq_CD7_Infusion_S1_L001_R1_001_fastq_gz_HHV6b.quant.txt") %>% 
  filter(!(V3 %in% c(6, 120))) %>%
  distinct() %>%
  group_by(V1) %>% summarize(tc = n()) %>%
  mutate(barcode = paste0(V1, "-1")) %>% arrange(desc(tc)) %>% data.frame()
count_df2

car_vec <- car_df <- fread("../data/car/CD7_Infusion_CAR.quant.txt") %>%
  mutate(barcode = paste0(V1, "-1")) %>% pull(barcode) %>% table()
count_df2$car_copyN <- ifelse(is.na(car_vec[count_df2$barcode]), 0, car_vec[count_df2$barcode])

exp <- Read10X_h5("../data/tenx_h5/CD7_Infusion_filtered_feature_bc_matrix.h5")
count_df2$called_cell <- count_df2$barcode %in% colnames(exp)
count_df2$valid_barcode <- count_df2$barcode %in% colnames(Read10X_h5("../data/tenx_h5/CD7_Infusion_raw_feature_bc_matrix.h5"))

count_df2$totalUMIs <- colSums(Read10X_h5("../data/tenx_h5/CD7_Infusion_raw_feature_bc_matrix.h5"))[count_df2$barcode]
head(count_df2)


# SEurat things

so <- CreateSeuratObject(counts = exp)

# dimension reduction
so <- NormalizeData(so) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% 
  FindNeighbors() %>%  FindClusters() %>% RunUMAP(dims = 1:30)
DimPlot(so, group.by = c("seurat_clusters"), label = FALSE)
FeaturePlot(so, features = c("CD7", "GZMB", "TREM1", "CD3D", "CD4", "CD8A"), sort = TRUE)
FeaturePlot(so, features = c("HHV6log", "MKI67"), sort = TRUE)


fm <- FindMarkers(so, "9")
head()