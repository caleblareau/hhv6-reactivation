library(data.table)
library(dplyr)
library(BuenColors)
library(cowplot)

fread("data.tsv") %>%
  mutate(conditions = c("aUntreated", "IFN-g", "TNF-a", "zIFN-g + TNF-a")) %>%
  ggplot(aes(x = conditions, y= TPM)) +
  geom_bar(stat = "identity", color = "black", fill = "lightgrey", width = 0.8) +
  pretty_plot(fontsize = 8) + L_border() + labs(x = "", y = "OX40 expression (transcripts per million)") +
  scale_y_continuous(expand = c(0,0)) -> plotME

ggsave2(plotME, filename = "CRS_endothelial.pdf", width = 2.5, height = 1.8)
