library(data.table)
library(dplyr)
library(stringr)

# Import serratus data
sra <- fread("../data/hu_SraRunInfo.csv.gz")
virus <- fread(paste0("../data/datapull-19February2022/SerratusMatches-HHV6B.csv.gz"))
human_reactivation <- merge(virus, sra, by.x = "run_id", by.y = "Run")

sra %>% filter(Run %in% c("SRR11684512"))

# Import meta data
gse12_raw <- fread("../data/tcell-followup/GSE127468-counts.txt")
heidelberg_df <- data.frame(
  id = gse12_raw$V1, 
  count = gse12_raw$V2,
  srr = gse12_raw$V4,
  str_split_fixed(gse12_raw$V3, "_", 3)
)

gse73_raw <- fread("../data/tcell-followup/GSE73213-counts.txt")
scripps_df <- data.frame(
  id = gse73_raw$V1, 
  count = gse73_raw$V2,
  srr = gse73_raw$V4,
  str_split_fixed(gse73_raw$V3, "-", 3)
)

# Merge with Serratus data
scripps_df_full <- merge(scripps_df, human_reactivation[,c("run_id", "n_reads", "score")] %>%
  filter(run_id %in% paste0("SRR24540", 50:81)), by.x = "srr", by.y = "run_id", all = TRUE)
scripps_df_full$n_reads <- ifelse(is.na(scripps_df_full$n_reads ), 0, scripps_df_full$n_reads )
scripps_df_full$pct <- round(scripps_df_full$n_reads/scripps_df_full$count *100,4)
scripps_df_full$day <- as.numeric(gsub("T", "", scripps_df_full$X3 ))/24
  
heidelberg_df_full <- merge(heidelberg_df, human_reactivation[,c("run_id", "n_reads", "score")] %>%
                           filter(run_id %in% paste0("SRR86461", 51:74)), by.x = "srr", by.y = "run_id", all = TRUE)
heidelberg_df_full$n_reads <- ifelse(is.na(heidelberg_df_full$n_reads ), 0, heidelberg_df_full$n_reads )
heidelberg_df_full$pct <- round(heidelberg_df_full$n_reads/heidelberg_df_full$count *100,4)
heidelberg_df_full$day <- as.numeric(gsub("day", "", heidelberg_df_full$X3 ))
  
library(BuenColors)
library(ggbeeswarm)
p1 <- ggplot(heidelberg_df_full, aes(x = day, y = pct, color = X1, shape = X2)) +
  geom_point() +  scale_y_log10() + geom_line() +
  scale_color_manual(values = c( "firebrick","blue", "black")) +
  pretty_plot(fontsize = 7) +
  scale_x_continuous(limits = c(3,15)) + L_border() +
  theme(legend.position = "bottom") + labs(x = "Day in culture", y = "% of RNA from HHV6")
cowplot::ggsave2(p1, file = "../plots/heidelberg_viz.pdf", width = 2, height = 2)

p2 <- ggplot(scripps_df_full, aes(x = day, y =  pct, color = as.character(X1))) +
  geom_point() + geom_line() + scale_y_log10() +
  facet_wrap(~X2) +
  scale_color_manual(values = c( "purple", "dodgerblue", "dodgerblue4", "purple4")) +
  pretty_plot(fontsize = 7) +
  scale_x_continuous(limits = c(0,15)) +
  theme(legend.position = "bottom") + labs(x = "Day in culture", y = "% of RNA from HHV6")
cowplot::ggsave2(p2, file = "../plots/scripps_viz.pdf", width = 4, height = 2)
