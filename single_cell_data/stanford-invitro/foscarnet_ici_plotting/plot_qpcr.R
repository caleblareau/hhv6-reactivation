library(data.table)
library(dplyr)
library(BuenColors)
library(ggbeeswarm)

qpcr_df <- fread("bence_qpcr.tsv") 
sem_df <- qpcr_df %>%
  group_by(V1) %>%
  summarize(sem = sqrt(var(V2))/sqrt(n()), V2 = mean(V2))
qpcr_df %>%
  ggplot(aes(x = V1, y = V2)) + 
  geom_bar(data = sem_df, stat = "identity", fill = "lightgrey", color = "black") +
  geom_errorbar(data = sem_df, aes(ymin=V2-sem, ymax=V2+sem), width=.2,
                position=position_dodge(.9))+
  geom_quasirandom() +
  pretty_plot(fontsize = 8) + L_border() + 
  labs(x = "Foscarnet treatment", y = "U31 RT-qPCR")  +
  scale_y_continuous(expand = c(0,0)) -> p1


cowplot::ggsave2(p1, file = "foscarnet_qpcr_hhv6.pdf", width = 2, height = 2)

# do some statistics tho we won't report them
ctr <- qpcr_df$V2[1:3]
mm25 <- qpcr_df$V2[4:6]
mm5 <- qpcr_df$V2[7:8]

t.test(ctr, mm25)
t.test(ctr, mm5)
