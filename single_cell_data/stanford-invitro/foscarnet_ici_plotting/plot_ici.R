library(data.table)
library(dplyr)
library(BuenColors)

fread("allogene_ici.csv") %>%
  reshape2::melt(id.vars = c("Day")) %>%
  ggplot(aes(x = Day, y = value, color = variable)) + 
  geom_hline(yintercept = 4.73, linetype = 2) + 
  geom_point() + scale_y_continuous(limits = c(3, 6)) +
  pretty_plot(fontsize = 8) + L_border() + 
  scale_color_manual(values = jdb_palette("corona")) +
  labs(x = "Day", y = "HHV-6 Genome equivalents") -> p1
cowplot::ggsave2(p1, file = "ici_hhv6.pdf", width = 3.5, height = 2)
