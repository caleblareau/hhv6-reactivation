library(data.table)
library(dplyr)
library(RcppRoll)

dv <- fread(paste0("cat ../data/qu2027_all_alignments.txt | awk '{print $1 , $2 , $3 , $4, $5, $6}'"), skip = 2, fill = TRUE, select = 1:6)
dv <- dv[dv$V5 >29,]

library(GenomicRanges)
# Expand based on the length of the read molecule
dd <- data.frame(
  chr = "chr1",
  start = ifelse(dv$V2 ==0, dv$V4, dv$V4 - 100),
  end = ifelse(dv$V2 ==0, dv$V4 + 100, dv$V4 )
) %>% makeGRangesFromDataFrame()
cov <- as.integer(coverage(dd)[["chr1"]])
ddf <- data.frame(
  index = roll_mean(1:length(cov), n = 1),
  coverage = roll_mean(cov, n = 1)
)

library(BuenColors)
p1 <- ggplot(ddf[seq(1, 162075, 20),], aes(x = index, y = coverage+1)) +
  geom_line(color = "dodgerblue2") + 
  scale_y_log10(limits = c(1, 100)) +
  pretty_plot(fontsize = 7) + L_border() +
  labs (y = "total coverage", x = "Position in HHVb reference genome") 
p1

cowplot::ggsave2(p1, file = "../plots/coverage_across_virus-quATAC-nosmooth.pdf", height = 1.5, width = 2.3)
