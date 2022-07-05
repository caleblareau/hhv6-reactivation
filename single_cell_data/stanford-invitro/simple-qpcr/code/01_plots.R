library(data.table)
library(dplyr)
library(BuenColors)


fread("../data/first_qPCR.txt") %>%
  ggplot(aes(x = Day, y = qPCR, color = Donor)) + 
  geom_point() + geom_line()  +
  scale_y_log10() +
  scale_color_manual(values = jdb_palette("corona")[c(2,5)]) +
  pretty_plot(fontsize = 7) + L_border()  + labs(x = "Days in Culture", y = "HHV-6B U31 qPCR") -> plong1
cowplot::ggsave2(plong1, file = "../plots/first_timecourse.pdf", width = 2.8, height = 1.2)


fread("../data/percell_qPCR.txt") %>%
  arrange((qPCR_U31_pCell)) %>%
  mutate(Donor = factor(as.character(Donor), levels = as.character(Donor)), rank = 1:n()) %>%
  data.frame() %>%
  ggplot(aes(x = Donor, y = qPCR_U31_pCell, fill = Donor)) + 
  geom_bar(stat = "identity", color = "black", width = 0.8) +
  scale_fill_manual(values = jdb_palette("corona")[4:1]) +
  pretty_plot(fontsize = 7) + L_border()  +
  scale_y_continuous(expand = c(0,0)) + theme(legend.position = "none") -> pbar
cowplot::ggsave2(pbar, file = "../plots/percell.pdf", width = 1.8, height = 1.2)


fread("../data/extended_culture.txt") %>%
  ggplot(aes(x = Day, y = U31_qPCR, color = Donor)) + 
  geom_point() + geom_line()  +
  scale_y_log10() +
  scale_color_manual(values = jdb_palette("corona")[c(2,10)]) +
  pretty_plot(fontsize = 7) + L_border()  + labs(x = "Days in Culture", y = "HHV-6B U31 qPCR")
