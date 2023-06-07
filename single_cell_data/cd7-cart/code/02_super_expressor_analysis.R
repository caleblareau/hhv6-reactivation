library(data.table)
library(Seurat)
library(dplyr)
library(BuenColors)
library(ggbeeswarm)
library(viridis)

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
merge(snp_df, soup_df, by.x = "cell_id", by.y = "barcode") %>%
  group_by(assignment, assign, status) %>%
  summarize(count = n())

mdf <- merge(soup_df, count_df2, all.x = TRUE, by.x = "barcode", by.y = "barcode") %>% # 
  arrange(desc(tc))
mdf$CD7 <- exp["CD7", mdf$barcode]
mdf$CD7cp10k <- exp["CD7", mdf$barcode]/colSums(exp[, mdf$barcode])*10000
mdf$GZMBcp10k <- exp["GZMB", mdf$barcode]/colSums(exp[, mdf$barcode])*10000

mdf$HHV6cp10k <- mdf$tc/(colSums(exp[, mdf$barcode])+ mdf$tc)*10000

mdf$called_cell <- ifelse(is.na(mdf$called_cell), TRUE, mdf$called_cell)
mdf$tc <- ifelse(is.na(mdf$tc), 0, mdf$tc)
mdf$car_copyN <- ifelse(is.na(mdf$car_copyN), 0, mdf$called_cell)

mdf$totalUMIs <- ifelse(is.na(mdf$totalUMIs), colSums(exp)[mdf$barcode], mdf$totalUMIs)

mdf %>% filter(called_cell & assignment %in% c("0", "1")) %>%
  group_by(assignment) %>%
  summarize(count = n(), total_umis = sum(totalUMIs),
            pct_positive = sum(tc > 0)/count,
            total_hhv6 = sum(tc), n_SE = sum(tc >= 10)) %>%
  mutate(total_hhv6/(total_umis + total_hhv6)*10000)
mdf %>% filter(called_cell & assignment %in% c("0", "1")) %>%
  filter(tc >0 )%>%
  group_by(assignment) %>% summarize(count = n())
P1 <- ggplot(mdf %>% filter(called_cell & assignment %in% c("0", "1")), aes(x =assignment, y = HHV6cp10k, color = tc >= 10)) + 
  geom_quasirandom(size = 0.5)  + pretty_plot(fontsize = 7) + L_border() + labs(x = "Genetic demultiplex annotation", y = "#HHV-6 RNAs per 10k UMIs") +
  scale_color_manual(values = c("black", "firebrick")) +
  scale_y_log10() + theme(legend.position = "none")
P1
cowplot::ggsave2(P1, file = "../plots/hhv6umis.pdf", width = 1.1, height = 1.5)

mdf %>% filter(called_cell & assignment %in% c("0", "1")) %>%
  filter(HHV6cp10k > 0) %>%
  group_by(assignment) %>% 
  summarize(nHigh = sum(HHV6cp10k > 20), max(HHV6cp10k))


### Other analyses
ggplot(mdf %>% filter(called_cell & assignment %in% c("0", "1")), aes(x = as.character(abs(as.numeric(assignment)-1)), y = GZMBcp10k)) + 
  geom_quasirandom(size = 0.5)  + pretty_plot(fontsize = 7) + L_border() + labs(x = "", y = "GZMB UMIs per 10k")

hhv6cp10k <- mdf %>% filter(called_cell & assignment %in% c("0", "1")) %>% pull(HHV6cp10k)
CD7cp10k <- mdf %>% filter(called_cell & assignment %in% c("0", "1")) %>% pull(CD7cp10k)
donor <- mdf %>% filter(called_cell & assignment %in% c("0", "1")) %>% pull(assignment)
str(wilcox.test( hhv6cp10k[donor == "0"],hhv6cp10k[donor == "1"]))
str(wilcox.test( CD7cp10k[donor == "0"],CD7cp10k[donor == "1"]))

P2 <- ggplot(mdf %>% filter(called_cell & assignment %in% c("0", "1")), aes(x = as.character(abs(as.numeric(assignment)-1)), y = log10(CD7cp10k + 1))) + 
  geom_quasirandom(size = 0.5)  + pretty_plot(fontsize = 7) + L_border() + labs(x = "", y = "CD7 UMIs per 10k")
P2
cowplot::ggsave2(P2, file = "../plots/cd7umis.pdf", width = 1.5, height = 1.2)

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

