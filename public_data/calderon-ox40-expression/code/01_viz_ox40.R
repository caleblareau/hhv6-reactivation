library(data.table)
library(dplyr)
library(BuenColors)

dt <- fread("../data/GSE118165_RNA_gene_abundance.txt")
mat <- data.matrix(data.frame(dt[,-1]))
rownames(mat) <- dt[["V1"]]
order_df <- data.frame(sn = colnames(dt)[-1], OX40 = round(mat["ENSG00000186827",]/colSums(mat) * 1000000, 1)) %>%
  arrange(desc(OX40)) %>%
  mutate(rank = 1:n())

# Append meta data
md <- fread("../data/calderon-metadata.txt")
modf <- merge(order_df, md, by.x = "sn", by.y = "ID")

p1 <- ggplot(modf %>% arrange(rank), aes(x = rank, y = OX40, color = lineage)) + 
  geom_point(size = 0.5) + scale_color_manual(values = jdb_palette("corona")[c(3,1,2,4,5,6)]) +
  pretty_plot() + L_border() + labs(x = "Rank ordered expression", y = "OX40 transcripts per million")
cowplot::ggsave2(p1, file = "../output/tpm_ox40.pdf", width = 6, height = 2)
