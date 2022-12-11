library(data.table)
library(dplyr)
library(BuenColors)
library(annotables)

dt <- fread("../data/GSE118165_RNA_gene_abundance.txt")
mat <- data.matrix(data.frame(dt[,-1]))
rownames(mat) <- dt[["V1"]]
cpm <- round(t(t(mat)/colSums(mat)*1000000),1)

vec <- grch38$symbol; names(vec) <- grch38$ensgene
vec2 <- grch38$ensgene; names(vec2) <- grch38$symbol

robbie_list <- vec2[c("CD247", "LCK", "ZAP70", "LCP2", "PLCG1", "NCK1", "SH2D1A", "LAT", "VAV1", "SYK", "SLAMF1", "MAPK3", "MAPK1")]
robbie_mat <- t(cpm[unname(robbie_list),]); colnames(robbie_mat) <- names(robbie_list)

order_df <- data.frame(sn = colnames(dt)[-1], robbie_mat) 

# Append meta data
md <- fread("../data/calderon-metadata.txt")
modf <- merge(order_df, md, by.x = "sn", by.y = "ID")
modf$stimstim <- ifelse(modf$treatment == "no_treament", "none", "stimulation")

melt_df <- modf %>%
  reshape2::melt(id.vars = c("Experiment", "GEO", "lineage", "treatment", "Cell_type", "Donor", "stimstim", "sn"))

library(viridis)
ggplot(melt_df %>% group_by(stimstim, lineage, variable, Cell_type) %>% top_n(1),
       aes(x = variable, y = Cell_type, fill = log10(value), label = value)) +
  geom_tile() + geom_text() + facet_grid(lineage~stimstim, scales = "free", space = "free") +
  scale_fill_viridis() + theme_bw() +
  theme(legend.position = "bottom") + labs(x = "Gene", y = "Celltype", color = "log10 counts per million") +
  ggtitle("Heatmap of RNA expression (counts per million); color is log10; text is value")

p1
cowplot::ggsave2(p1, file = "../output/tpm_ox40.pdf", width = 6, height = 2)
