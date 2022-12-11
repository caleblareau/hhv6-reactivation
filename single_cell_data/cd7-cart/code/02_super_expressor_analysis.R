library(data.table)
library(Seurat)
library(dplyr)
library(BuenColors)

#  Define hhv6 expressing cells from broad script
count_df2 <- fread("../data/hhv6/fastq_CD7_Day19_S3_L001_R1_001_fastq_gz_HHV6b.quant.txt") %>% 
  filter(!(V3 %in% c(6, 120))) %>%
  distinct() %>%
  group_by(V1) %>% summarize(tc = n()) %>%
  mutate(barcode = paste0(V1, "-1")) %>% arrange(desc(tc)) %>% data.frame()

car_vec <- car_df <- fread("../data/car/CD7_D19_CARquant.txt") %>%
  mutate(barcode = paste0(V1, "-1")) %>% pull(barcode) %>% table()
count_df2$car_copyN <- ifelse(is.na(car_vec[count_df2$barcode]), 0, car_vec[count_df2$barcode])

exp <- Read10X_h5("../data/tenx_h5/CD7_Day19_filtered_feature_bc_matrix.h5")
count_df2$called_cell <- count_df2$barcode %in% colnames(exp)
count_df2$valid_barcode <- count_df2$barcode %in% colnames(Read10X_h5("../data/tenx_h5/CD7_Day19_raw_feature_bc_matrix.h5"))

count_df2$totalUMIs <- colSums(Read10X_h5("../data/tenx_h5/CD7_Day19_raw_feature_bc_matrix.h5"))[count_df2$barcode]
head(count_df2)

# Look at snps
snp_df <- fread("../data/snp_assignment/Day19_SNPassignments.tsv")
mdf <- merge(snp_df, count_df2, all.x = TRUE, by.x = "cell_id", by.y = "barcode") %>%
  arrange(desc(tc))
mdf$CD7 <- exp["CD7", mdf$cell_id]
mdf$CD7cp10k <- exp["CD7", mdf$cell_id]/colSums(exp[, mdf$cell_id])*10000
mdf$HHV6cp10k <- mdf$tc/(colSums(exp[, mdf$cell_id])+ mdf$tc)*10000

table(mdf$assign)
library(ggbeeswarm)
library(viridis)
P1 <- ggplot(mdf %>% filter(called_cell & assign %in% c("D1", "D2")), aes(x = assign, y = HHV6cp10k)) + 
  geom_quasirandom(size = 0.5)  + pretty_plot(fontsize = 7) + L_border() + labs(x = "", y = "HHV6 UMIs per 10k")
cowplot::ggsave2(P1, file = "../plots/hhv6umis.pdf", width = 1.8, height = 1.2)

hhv6cp10k <- mdf %>% filter(called_cell & assign %in% c("D1", "D2")) %>% pull(HHV6cp10k)
CD7cp10k <- mdf %>% filter(called_cell & assign %in% c("D1", "D2")) %>% pull(CD7cp10k)
donor <- mdf %>% filter(called_cell & assign %in% c("D1", "D2")) %>% pull(assign)
str(wilcox.test( hhv6cp10k[donor == "D1"],hhv6cp10k[donor == "D2"]))
str(wilcox.test( CD7cp10k[donor == "D1"],CD7cp10k[donor == "D2"]))

P2 <- ggplot(mdf %>% filter(called_cell & assign %in% c("D1", "D2")), aes(x = assign, y = CD7cp10k)) + 
  geom_quasirandom(size = 0.5)  + pretty_plot(fontsize = 7) + L_border() + labs(x = "", y = "CD7 UMIs per 10k")
cowplot::ggsave2(P2, file = "../plots/cd7umis.pdf", width = 1.8, height = 1.2)

ggplot(mdf %>% filter(called_cell) %>% arrange(CD7cp10k), aes(x = D1_sum, y = D2_sum, color = CD7cp10k)) + 
  geom_quasirandom()  + pretty_plot() + scale_color_viridis()


mdf %>% 
  head(count_df2, 50)

mdf %>%
  filter(called_cell) %>%
  group_by(assign) %>% summarize(sum(tc >= 10))

tail(exp)
rownames(exp)

genes <- c("CD3D","CD3E", "CD4", "CD8A", "CD8B", "CD7")
exp[genes,colnames(exp) %in% (count_df2 %>% filter(tc > 23) %>% pull(barcode))]

mdf %>% group_by(assign) %>% summarize(mean(CD7))


####


#  Define infusion

#  Define hhv6 expressing cells from broad script
count_df2 <- fread("../data/hhv6/fastq_CD7_Infusion_S1_L001_R1_001_fastq_gz_HHV6b.quant.txt") %>% 
  filter(!(V3 %in% c(6, 120))) %>%
  distinct() %>%
  group_by(V1) %>% summarize(tc = n()) %>%
  mutate(barcode = paste0(V1, "-1")) %>% arrange(desc(tc)) %>% data.frame()

car_vec <- car_df <- fread("../data/car/CD7_Infusion_CARquant.txt") %>%
  mutate(barcode = paste0(V1, "-1")) %>% pull(barcode) %>% table()
count_df2$car_copyN <- ifelse(is.na(car_vec[count_df2$barcode]), 0, car_vec[count_df2$barcode])

exp <- Read10X_h5("../data/tenx_h5/CD7_Infusion_filtered_feature_bc_matrix.h5")
count_df2$called_cell <- count_df2$barcode %in% colnames(exp)
count_df2$valid_barcode <- count_df2$barcode %in% colnames(Read10X_h5("../data/tenx_h5/CD7_Infusion_raw_feature_bc_matrix.h5"))

count_df2$totalUMIs <- colSums(Read10X_h5("../data/tenx_h5/CD7_Infusion_raw_feature_bc_matrix.h5"))[count_df2$barcode]
count_df2


#  Define hhv6 expressing cells from broad script
count_df2 <- fread("../data/hhv6/fastq_CD7_Day14_S2_L001_R1_001_fastq_gz_HHV6b.quant.txt") %>% 
  filter(!(V3 %in% c(6, 120))) %>%
  distinct() %>%
  group_by(V1) %>% summarize(tc = n()) %>%
  mutate(barcode = paste0(V1, "-1")) %>% arrange(desc(tc)) %>% data.frame()

car_vec <- car_df <- fread("../data/car/CD7_D14_CARquant.txt") %>%
  mutate(barcode = paste0(V1, "-1")) %>% pull(barcode) %>% table()
count_df2$car_copyN <- ifelse(is.na(car_vec[count_df2$barcode]), 0, car_vec[count_df2$barcode])

exp <- Read10X_h5("../data/tenx_h5/CD7_Day14_filtered_feature_bc_matrix.h5")
count_df2$called_cell <- count_df2$barcode %in% colnames(exp)
count_df2$valid_barcode <- count_df2$barcode %in% colnames(Read10X_h5("../data/tenx_h5/CD7_Day14_raw_feature_bc_matrix.h5"))

count_df2$totalUMIs <- colSums(Read10X_h5("../data/tenx_h5/CD7_Day14_raw_feature_bc_matrix.h5"))[count_df2$barcode]

# Look at snps
snp_df <- fread("../data/snp_assignment/Day14_SNPassignments.tsv")
mdf <- merge(snp_df, count_df2, all.x = TRUE, by.x = "cell_id", by.y = "barcode") %>%
  arrange(desc(tc))
mdf$CD7 <- exp["CD7", mdf$cell_id]

table(mdf$assign)
ggplot(mdf, aes(x = assign, y = tc)) + geom_violin() 
mdf %>% 
  head(, 50)
tail(exp)
rownames(exp)

genes <- c("CD3D","CD3E", "CD4", "CD8A", "CD8B", "CD7")
exp[genes,colnames(exp) %in% (count_df2 %>% filter(tc > 23) %>% pull(barcode))]

mdf %>% group_by(assign) %>% summarize(mean(CD7))

