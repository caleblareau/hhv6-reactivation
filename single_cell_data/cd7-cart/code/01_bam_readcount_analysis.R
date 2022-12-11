library(data.table)
library(dplyr)
library(matrixStats)

process_lib <- function(lib){
  dt <- fread(lib) %>%
    filter(TOTAL >= 1) 
  mat <- data.matrix(data.frame(dt[,c("A", "C", "G", "T")]))
  dt$max_allele <-c("A", "C", "G", "T")[max.col(mat)]
  dt$AF <- round(100*(rowMaxs(mat)/ dt$TOTAL),0)
  dt %>% mutate(mutation = max_allele != REF) %>%
    filter(BP > 8788 & BP < 153321)
  
}
d19 <- process_lib("../data/bam-readcount/D19.hhv6b.se.bam.br.txt")
d14 <- process_lib("../data/bam-readcount/D14.hhv6b.se.bam.br.txt")
infusion <- process_lib("../data/bam-readcount/D0infusion.hhv6b.se.bam.br.txt")
pub1 <- process_lib("../data/bam-readcount/SRR2454080.hhv6b.se.bam.br.txt")
pub2 <- process_lib("../data/bam-readcount/SRR8646162.hhv6b.paired.bam.br.txt")

muts <- infusion %>% filter(AF > 90 & TOTAL >= 2) %>%
  filter(max_allele != REF) %>%
  pull(BP)
d19 %>% filter(BP %in% muts)  #%>% filter(!mutation)
d14 %>% filter(BP %in% muts)  #%>% filter(!mutation)

pub1 %>% filter(BP %in% muts)%>% filter(!mutation)
pub2 %>% filter(BP %in% muts) %>% filter(!mutation) %>%
  data.frame()
