library(data.table)
library(Seurat)
library(dplyr)
library(BuenColors)
source("00_blanco.R")

#  Define hhv6 expressing cells from broad script
count_df2 <- fread("../data/hhv6/fastq_CD7_Day19_S3_L001_R1_001_fastq_gz_HHV6b.quant.txt") %>% 
  filter(!(V3 %in% c(6, 120))) %>%
  distinct() %>%
  group_by(V1) %>% summarize(tc = n()) %>%
  mutate(barcode = paste0(V1, "-1")) %>% arrange(desc(tc)) %>% data.frame()

car_vec <- car_df <- fread("../data/car/CD7_Day19_CAR.quant.txt") %>%
  mutate(barcode = paste0(V1, "-1")) %>% pull(barcode) %>% table()
count_df2$car_copyN <- ifelse(is.na(car_vec[count_df2$barcode]), 0, car_vec[count_df2$barcode])

exp <- Read10X_h5("../data/tenx_h5/CD7_Day19_filtered_feature_bc_matrix.h5")
count_df2$called_cell <- count_df2$barcode %in% colnames(exp)
count_df2$valid_barcode <- count_df2$barcode %in% colnames(Read10X_h5("../data/tenx_h5/CD7_Day19_raw_feature_bc_matrix.h5"))

count_df2$totalUMIs <- colSums(Read10X_h5("../data/tenx_h5/CD7_Day19_raw_feature_bc_matrix.h5"))[count_df2$barcode]
head(count_df2)

# Look at snps
snp_df <- fread("../data/snp_assignment/Day19_mitoSNPassignments.tsv")
soup_df <- fread("../data/snp_assignment/CD7_Day19_clusters.tsv")
soup_mito <- merge(snp_df, soup_df, by.x = "cell_id", by.y = "barcode") 
soup_mito$barcode <- soup_mito$cell_id
mdf <- merge(soup_df, count_df2, all.x = TRUE, by.x = "barcode", by.y = "barcode") %>% # 
  arrange(desc(tc))
mdf$CD7 <- exp["CD7", mdf$barcode]
mdf$CD7cp10k <- exp["CD7", mdf$barcode]/colSums(exp[, mdf$barcode])*10000
mdf$HHV6cp10k <- mdf$tc/(colSums(exp[, mdf$barcode])+ mdf$tc)*10000
mdf$HHV6log <- log10(mdf$tc/(colSums(exp[, mdf$barcode])+ mdf$tc)*10000 + 1)
head(mdf)


# SEurat things
mdf <- data.frame(mdf)
rownames(mdf) <- mdf$barcode
so <- CreateSeuratObject(counts = exp[, mdf$barcode], meta.data = mdf)
so <- so[,c(so$assignment %in% c("0", "1"))]
so$HHV6log <- ifelse(is.na(so$HHV6log), 0, so$HHV6log)
so <- NormalizeData(so) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% 
  FindNeighbors() %>%  FindClusters(resolution = 0.3) %>% RunUMAP(dims = 1:25)
DimPlot(so, group.by = c( "assignment", "seurat_clusters"), label = TRUE) 

so@meta.data %>%
  group_by(seurat_clusters) %>%
  summarize(count = sum(tc, na.rm = TRUE))
write.table(rownames(so@meta.data %>% filter(seurat_clusters %in% c(5))), row.names = FALSE, col.names = FALSE, quote = FALSE,
            file = "pot_hhv_for_cnv.tsv")


FeaturePlot(so, features = c("CD7", "GZMB", "TREM1", "CD3D", "CD4", "CD8A", "MS4A1", "MKI67", "HHV6log", "ALAS2", "PDCD1", "MYC", "NOTCH1"), sort = TRUE)


###3
assign_plot <- DimPlot(so, group.by = c( "assignment"), label = FALSE, pt.size = 1) +
  theme_blank() +
  scale_color_manual(values = jdb_palette("corona")[c(2,1)]) +
  theme(legend.position = "none")

cowplot:::ggsave2(assign_plot, 
                  filename = "../plots/assign_cd7_plot.png", width = 6, height = 6, dpi = 500)

assign_plot_z <- DimPlot(so, group.by = c( "assignment"), label = FALSE, pt.size = 2) +
  theme_blank() +coord_cartesian(xlim = c(4,8.8), ylim = c(2,6.8)) +
  scale_color_manual(values = jdb_palette("corona")[c(2,1)]) +
  theme(legend.position = "none")
cowplot:::ggsave2(assign_plot_z,
                  filename = "../plots/ZOOMEDassign_cd7_plot.png", width = 2, height = 2, dpi = 900)


#######
cowplot:::ggsave2(FeaturePlot(so, features = c("HHV6log"), sort = TRUE, pt.size = 1) + scale_color_viridis() + theme_blank() + theme(legend.position = "none"),
                  filename = paste0("../plots/hhv6log_viridis_plot.png"), width = 6, height = 6, dpi = 500)
# now plot markers
make_plot <- function(one){
  cowplot:::ggsave2(FeaturePlot(so, one, sort.cell = TRUE, pt.size = 1) + theme_blank() + theme(legend.position = "none"),
                    filename = paste0("../plots/",one,"_plot.png"), width = 6, height = 6, dpi = 500)
  one
}
lapply(c("CD7", "GZMB", "TREM1", "CD3D", "CD4", "CD8A", "MS4A1", "MKI67", "HHV6log", "ALAS2", "PDCD1", "MYC", "CD8B"), make_plot)

make_plot_zoom <- function(one){
  cowplot:::ggsave2(FeaturePlot(so, one, sort.cell = TRUE, pt.size = 2) + theme_blank() + theme(legend.position = "none") +
                      coord_cartesian(xlim = c(4,8.5), ylim = c(2,6.5)),
                    filename = paste0("../plots/zoomp/",one,"_zoom_plot.png"), width = 2, height = 2, dpi = 900)
  one
}
lapply(c("CD7", "GZMB", "TREM1", "CD3D", "CD4", "CD8A", "MS4A1", "MKI67", "HHV6log", "ALAS2", "PDCD1", "MYC", "CD8B"), make_plot_zoom)


