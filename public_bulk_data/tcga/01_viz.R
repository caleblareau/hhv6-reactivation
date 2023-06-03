library(data.table)
library(BuenColors)
library(dplyr)
library(ggbeeswarm)

fread('ox40_tcga.tsv') %>%
  mutate(tcga = gsub("TCGA-", "", tcgaid)) %>%
  mutate(ox40cpm = OX40/total*1000000) %>%
  ggplot(aes(x = reorder(tcga, -1*ox40cpm, FUN = median, ), y = ox40cpm)) + 
  geom_boxplot(outlier.shape = NA) + 
  scale_y_log10() +
  coord_cartesian(ylim= c(5e-2, 2e2)) +
  pretty_plot(fontsize = 8) + L_border() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  labs(x = "TCGA tumor type (rank ordered)", y = "OX40 counts per million") -> p1
cowplot::ggsave2(p1, file = "ox40_tcga.pdf", width = 6.5, height = 1.8)
