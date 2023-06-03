library(data.table)
library(dplyr)
library(viridis)
library(BuenColors)

data.frame(
  donor = c("D34", "D38", "D97", "D98"), 
  qPCR = c(0.37, 0.0015, 0.059, 0.90),
  pct_SE = c(0.14, 0, 0.01, 0.28)
) %>%
  ggplot(aes(x = pct_SE, y = qPCR)) + 
  geom_point() +
  geom_smooth(method='lm', formula= y~x, se = FALSE, linetype = 2) + 
  pretty_plot(fontsize = 7) + L_border() + 
  labs(x = "% of scRNA-seq cells with high HHV-6", y = "Sample U31 qPCR / cell") -> p1
cowplot::ggsave2(p1, file = "../plots/bulk_plot.pdf", width = 1.5, height = 1.5)
  