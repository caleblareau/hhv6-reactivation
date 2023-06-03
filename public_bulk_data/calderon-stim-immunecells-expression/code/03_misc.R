library(data.table)
library(dplyr)
library(BuenColors)
library(annotables)
msy_genes <- c("RPS4Y1", "ZFY", "TBL1Y", "USP9Y", "DDX3Y", "UTY", "TMSB4Y", "NLGN4Y", "KDM5D", "EIF1AY")

ens <- grch38 %>% filter(chr == "Y" & symbol %in% msy_genes) %>% pull(ensgene)
ens_vec <- grch38 %>% filter(chr == "Y" & symbol %in% msy_genes) %>% pull(symbol); names(ens_vec) <- ens
dt <- fread("../data/GSE118165_RNA_gene_abundance.txt")

dt2 <- dt %>% filter(V1 %in% ens)
mat_big <- data.matrix(data.frame(dt[,-1]))
mat_big_cpm <- t(t(mat_big)/colSums(mat_big)*1000000)
mat <- mat_big_cpm[dt$V1 %in% ens,]
rownames(mat) <- ens_vec[as.character(dt2$V1)]
mat
pheatmap(log1p(mat[,grepl("^X1002", colnames(mat))]))

# Append meta data
md <- fread("../data/calderon-metadata.txt")
modf <- merge(order_df, md, by.x = "sn", by.y = "ID")
modf$stimstim <- ifelse(modf$treatment == "no_treament", "none", "astimulation")
modf <- modf %>%
  group_by(stimstim) %>%
  arrange(desc(CD21)) %>%
  mutate(rank2 = 1:n())

p1 <- ggplot(modf %>% arrange(rank), aes(x = rank2, y = CD21, color = lineage)) + 
  geom_point(size = 0.5) + scale_color_manual(values = jdb_palette("corona")[c(3,1,2,4,5,6)]) +
  pretty_plot(fontsize = 8) + labs(x = "Rank ordered expression", y = "CD21 transcripts per million") +
  facet_wrap(~stimstim)
p1
cowplot::ggsave2(p1, file = "../output/tpm_CD21.pdf", width = 6, height = 2)
