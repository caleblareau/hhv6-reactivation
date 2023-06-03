library(data.table)
library(dplyr)
wl <- fread("../data/737K-august-2016.txt.gz")

whichone <- c("7")
bc <- fread(paste0("../data/barcodes/StJ",whichone,"_rna_barcodes.tsv.gz"))[[1]]
hhv6df <- fread(paste0("../data/hhv6/fastqs_StJ",whichone,"_S",whichone,"_R1_001_fastq_gz_HHV6b.quant.txt")) %>%
  filter(!(V3 %in% c(6,120))) %>% mutate(barcode = paste0(V1, "-1")) %>%
  mutate(is_cell = barcode %in% bc)
hhv6df

lapply(c("5", "7", "10", "14"), function(whichone){
  fread(paste0("../data/hhv6/fastqs_StJ",whichone,"_S",whichone,"_R1_001_fastq_gz_HHV6b.quant.txt")) %>%
    filter(!(V3 %in% c(6,120))) %>% mutate(barcode = paste0(V1, "-1")) %>%
    mutate(is_cell = barcode %in% bc) %>%pull(V3)
}) %>% reshape2::melt() %>% pull(value) %>% unique() %>% length()
