library(data.table)
library(stringr)
library(Seurat)
library(dplyr)

# Import data
high <- fread("../data/ID34high_HHV6b.kb.txt") %>% filter(V3 != 120)
txt <- data.frame(fread("../data/transcripts.txt", header = FALSE))
mat <- Read10X_h5("../data/ID34high.filtered_feature_bc_matrix.h5")
umi_vec <- colSums(mat)
ngenes_vec <- colSums(mat > 0)
ddf <- data.frame(table(paste0(high[[1]], "-1")))
ddf$ngenes <- ngenes_vec[as.character(ddf$Var1)]; ddf$ngenes <- ifelse(is.na(ddf$ngenes), 0, ddf$ngenes)
ddf$numis <- umi_vec[as.character(ddf$Var1)]; ddf$numis <- ifelse(is.na(ddf$numis), 0, ddf$numis)

ddf %>% arrange(desc(Freq)) %>% head(20)
