library(data.table)
library(dplyr)
library(BuenColors)
library(ggbeeswarm)
library(stringr)

dt <- fread("../data/Qu2017_summary_statistics.txt") 
dt$donor <- str_split_fixed(dt$source_name, "_", 3)[,2]
dt$donor <- factor(as.character(dt$donor ), levels = rev(unique(dt$donor )))
p1 <- dt %>% filter(donor != "") %>%
  ggplot(aes(x = nHHV6reads, y = donor, color = donor)) +
  geom_quasirandom() + scale_color_manual(values = jdb_palette("corona")) + 
  pretty_plot(fontsize = 8) + L_border() + theme(legend.position = "none") + 
  labs(x = "# of ATAC-seq reads mapping to HHV-6B genome", y = "")
cowplot::ggsave2(p1, file = "../plots/ctcl_atac.pdf", width = 3, height = 3)
