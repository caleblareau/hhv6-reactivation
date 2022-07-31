library(data.table)
library(dplyr)
library(BuenColors)

dt <- fread("../data/GSE118165_RNA_gene_abundance.txt")
mat <- data.matrix(data.frame(dt[,-1]))
rownames(mat) <- dt[["V1"]]
gene = "ENSG00000104490"
order_df <- data.frame(sn = colnames(dt)[-1], gexp = round(mat[gene,]/colSums(mat) * 1000000, 1)) %>%
  arrange(desc(gexp)) %>%
  mutate(rank = 1:n())
order_df



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
