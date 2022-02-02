library(data.table)
library(dplyr)

sra <- fread("hu_SraRunInfo.csv.gz")
hhv6 <- fread("hhv6b-serratus.csv")
human_hhv6b_reactivation <- merge(hhv6, sra, by.x = "run_id", by.y = "Run")

human_hhv6b_reactivation_filt <- human_hhv6b_reactivation %>%
  filter(n_reads >= 3) %>%
  arrange(desc(n_reads))

human_hhv6b_reactivation_filt[1:10,c("run_id", "n_reads", "score","SampleName")]

human_hhv6b_reactivation[,c("run_id", "n_reads", "score","SampleName")] %>%
  filter(SampleName %in% paste0("GSM1888", 809:840))

human_hhv6b_reactivation %>%
  filter(SampleName == "SAMD00009187")

human_hhv6b_reactivation %>% 
  filter(percent_identity >= 90 & n_reads >= 2) %>%
  group_by(SRAStudy) %>% summarize(n = n(), total = sum(n_reads)) %>%
  arrange(desc(total)) %>%
  data.frame()

human_hhv6b_reactivation_filt[,c("run_id", "n_reads", "score","SampleName", "SRAStudy", "percent_identity")] %>%
  filter(SRAStudy == "SRP089788")

