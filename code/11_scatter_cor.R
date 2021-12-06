library(dplyr)
library(data.table)
library(BuenColors)
library(ggrepel)

pbmc <- fread("../output/PBMC-Day15-correlation.tsv")
s34_d7 <- fread("../output/Sample34-Day7-correlation.tsv")
mdf <- merge(pbmc, s34_d7, by = "gene")

mdf <- mdf %>%
  mutate(label_rep = case_when(
    cor.x < -0.3 ~ gene,
    cor.x > 0.3  ~ gene,
    cor.y > 0.37  ~ gene,
    cor.y < -0.35  ~ gene,
    
    TRUE ~ ""
  ))
cor(mdf[,2], mdf[,3])
ggplot(mdf, aes(x = cor.x, y = cor.y, label = label_rep)) +
  geom_point(aes( color = label_rep != "")) + 
  
  geom_text_repel(max.overlaps = Inf, min.segment.length = 0, seed = 42, box.padding = 0.5,
                  fill = "white", xlim = c(-Inf, Inf), ylim = c(-Inf, Inf)) + 
  scale_color_manual(values = c("black", "red")) + 
  pretty_plot(fontsize = 20) + L_border()+  theme(legend.position = "none") +
  labs(x = "Cor. w/ HHV6 in PBMC infection (Day 15)", 
       y = "Cor. w/ HHV6 in PPD20034 (Day 7)")

