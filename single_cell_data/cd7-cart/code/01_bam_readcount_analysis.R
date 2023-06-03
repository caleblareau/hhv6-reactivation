library(data.table)
library(dplyr)
library(matrixStats)


###

all_possible_muts <- function(lib){
  dt <- fread(lib) 
  mat <- (data.frame(dt[,c("A", "C", "G", "T", "TOTAL", "BP", "REF")]))
  reshape2::melt(mat, id.vars = c("BP", "TOTAL", "REF")) %>%
    mutate(var = paste0(REF, BP, variable), af = value/ TOTAL)
}
mdf10 <- merge(all_possible_muts("../data/bam-readcount/SRR2454080.hhv6b.se.bam.br.txt"),
               all_possible_muts("../data/bam-readcount/SRR8646162.hhv6b.paired.bam.br.txt"),
               by = "var", all = TRUE) %>%
  filter(TOTAL.x >= 10 & TOTAL.y >=10 )

# This has the scRNA-seq coverage
well_covered <- fread("../data/bam-readcount/D19.hhv6b.paired.bam.br.txt") %>% filter(TOTAL > 100) %>% pull(BP)
mdf10 %>% filter(af.x > 0.95 & af.y > 0.95 & REF.y != variable.y) %>% filter(BP.x %in% well_covered)
mdf10 %>% filter(af.x > 0.95 & af.y < 0.05 & REF.y != variable.y) %>% filter(BP.x %in% well_covered)
mdf10 %>% filter(af.x < 0.05 & af.y > 0.95 & REF.y != variable.y) %>% filter(BP.x %in% well_covered)

####
mdf10 <- merge(all_possible_muts("../data/bam-readcount/D19.hhv6b.se.bam.br.txt"),
               all_possible_muts("../data/bam-readcount/D0infusion.hhv6b.se.bam.br.txt"),
               by = "var", all = TRUE) %>%
  filter(TOTAL.x >= 5 & TOTAL.y >= 5)
mdf10 %>% filter(af.x > 0.95 & af.y > 0.95 & REF.y != variable.y)
mdf10 %>% filter(af.x > 0.95 & af.y < 0.05 & REF.y != variable.y) 
mdf10 %>% filter(af.x < 0.05 & af.y > 0.95 & REF.y != variable.y)


