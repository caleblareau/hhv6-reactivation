library(data.table)
library(dplyr)
ff <- list.files("data/")
md <- fread("chipseq-metadata.txt")

counts <- sapply(md$SRR, function(i){
  print(i)
  # Written in a way to facilitate splits
  files = list.files("data/", pattern = i, full.names = TRUE)
  cmds = paste0("zcat < ",files," | awk '{print $1 , $2 , $3 , $4, $5, $6}'")
  
  # Filter for non DR segments
  #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC112820/#:~:text=As%20sequenced%2C%20the%20HHV%2D6B,present%20only%20in%20HHV%2D6B.
  
  totals <- sapply(cmds, function(cmd){
    dv <- fread(cmd, skip = 2, fill = TRUE, select = 1:6) %>%
      filter(V4 > 8793 & V4 < ( 8793 + 144528) & V5 > 29)  %>%
      dim()
    dv[1]
  })
  sum(totals)
})
md$counts <- counts
md$pct <- md$counts / md$nReads *100
md$DayStimCharacter <- factor(as.character(md$DayStim), levels = c("0", "1", "5", "14"))
mds <- md %>% group_by(DayStimCharacter, DonorID, CellInput) %>%
  summarize(nc = sum(counts), nr = sum(nReads)) %>%
  mutate(pct = nc/nr*100)

p2 <- ggplot(mds %>% filter(DayStimCharacter != "1") , aes(x = as.numeric(as.character(DayStimCharacter)), y = pct, color = as.character(DonorID), group = as.character(DonorID)))+
  facet_wrap(~CellInput)  +
  geom_point() + geom_line() +
  scale_y_log10() +
  scale_color_manual(values = c( "purple", "dodgerblue", "dodgerblue4", "purple4")) +
  pretty_plot(fontsize = 7) + L_border() +
  scale_x_continuous(limits = c(0,15))+
  theme(legend.position = "none") + labs(x = "Day in culture", y = "% of DNA from HHV6")
cowplot::ggsave2(p2, file = "plots/scripps_chipseq_viz.pdf", width = 4, height = 1.65 )

mds %>% arrange(desc(pct))
#----------------------------
# Plot the chipseq input library with several files

# Make a visualization across the reference genome
dvA <- fread(paste0("zcat < data/SRR2453798a.alignments.txt.gz | awk '{print $1 , $2 , $3 , $4, $5, $6}'"), skip = 2, fill = TRUE, select = 1:6)
dvB <- fread(paste0("zcat < data/SRR2453798b.alignments.txt.gz | awk '{print $1 , $2 , $3 , $4, $5, $6}'"), skip = 2, fill = TRUE, select = 1:6)
dv <- rbind(dvA, dvB)
dv <- dv[dv$V5 >29,]

dv <- fread(paste0("zcat < data/SRR2453739.alignments.txt.gz | awk '{print $1 , $2 , $3 , $4, $5, $6}'"), skip = 2, fill = TRUE, select = 1:6)
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
  index = 1:length(cov), 
  coverage = cov
)

library(BuenColors)
p1 <- ggplot(ddf[seq(1, 162075, 20),], aes(x = index, y = coverage+1)) +
  geom_line(color = "dodgerblue2") + 
  scale_y_log10(limits = c(1, 2000)) +
  pretty_plot(fontsize = 7) + L_border() +
  labs (y = "total coverage", x = "Position in HHVb reference genome") 
p1

cowplot::ggsave2(p1, file = "plots/coverage_across_virus-SRR2453784.pdf", height = 1.5, width = 2.3)
