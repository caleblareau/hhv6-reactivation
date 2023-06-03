library(data.table)
library(dplyr)
library(Seurat)

library(annotables)
gene_vec <- grch38$symbol; names(gene_vec) <- grch38$ensgene

process_channel <- function(channel){

  # Import kallisto dat
  car <- fread(paste0("../data/new-data/features/",channel,"_CAR.quantFinal.txt")) %>%
    mutate(bc = paste0(V1, "-1"))
  hhv6 <- fread(paste0("../data/new-data/features/",channel,"_HHV6b.quantFinal.txt")) %>%
    filter(!(V3 %in% c(6, 120))) %>%
    mutate(bc = paste0(V1, "-1"))
  print(dim(hhv6))
  hhv6p_cells <- unique(hhv6$bc)
  print(length(hhv6p_cells))
  hhv6p_count <- table(hhv6[["bc"]])
  
  # Import Gex data
  bcs <- paste0(fread(paste0("../data/new-data/kallisto-scrna/",channel,"_counts.barcodes.txt.gz"), header = FALSE)[[1]], "-1")
  genes <- fread(paste0("../data/new-data/kallisto-scrna/",channel,"_counts.genes.txt.gz"), header = FALSE)[[1]]
  genes_common <- make.unique(unname(gene_vec[stringr::str_split_fixed(genes, "[.]", 2)[,1]]))
  genes_common[is.na(genes_common)] <- "IDK"
  mtx <- fread(paste0("../data/new-data/kallisto-scrna/",channel,"_counts.mtx.gz"), skip = 4)
  library(Matrix)
  gex <- Matrix::sparseMatrix(
      i = c(mtx[[1]],length(bcs)),
      j = c(mtx[[2]],length(genes_common)),
      x = c(mtx[[3]],1)
  ) %>% t()
  print(sum(colSums(gex) > 800  ))
  colnames(gex) <- bcs; rownames(gex) <- genes_common
  gex_filtered <- gex[, colnames(gex) %in% hhv6p_cells]
  
  
  car_seq_count <- table(car[["bc"]])[colnames(gex_filtered)]
  
  
  # Compbine it all
  data.frame(
    n_car = as.numeric(car_seq_count),
    n_hhv6 = as.numeric(hhv6p_count[colnames(gex_filtered)]),
    t(data.matrix(gex_filtered[c("CD3D","CD3E", "CD4", "CD8A", "CD8B"),])), # "CCL5", "MS4A1", "TREM1"
    nGExpUMIs = colSums(gex_filtered), 
    pctMito = colSums(gex_filtered[grepl("^MT", rownames(gex_filtered)),])/colSums(gex_filtered)*100,
    channel = channel
  ) %>% arrange(desc(n_hhv6))
}

df <- process_channel(channel = "BBB7")
filt_df <- df %>% filter(nGExpUMIs > 800) %>% arrange(desc(n_hhv6))

write.table(filt_df, file = "../output/HHV6-positive-cells-reseq.tsv", 
            sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

fix <- function(vec){
  return(ifelse(is.na(vec), 0, vec))
}

plot_df <- data.frame(
  patient.unique  = make.unique(filt_df$patient),
  n_hhv6 = filt_df$n_hhv6,
  CAR = case_when(filt_df$channel %in% c("F1") ~ "F", 
                  filt_df$n_car > 0 ~ "T", 
                  TRUE ~ ""),
  aTCRab = fix(filt_df$nTCRA) + fix(filt_df$nTCRB),
  CD3 = filt_df$CD3D + filt_df$CD3E ,
  cytotoxic = filt_df$GZMK  + filt_df$GNLY+  filt_df$KLRG1 + filt_df$ZEB2 + filt_df$NKG7,
  CD8_4 = filt_df$CD8A + filt_df$CD8B +  filt_df$CD4
)

ptn <- as.character(plot_df$patient.unique)
library(BuenColors)
library(forcats)

p1 <- reshape2::melt(plot_df, id.vars = "patient.unique") %>%
  mutate(name = factor(as.character(patient.unique), levels = rev(ptn))) %>%
  ggplot(aes(x = variable, y = name, label = value)) + 
  geom_text(size = 2) + geom_tile(fill = NA, color = "black") +
  pretty_plot(fontsize = 8) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) + labs(x = "", y = "")

cowplot::ggsave2(p1, file = "../plots/grid.pdf", width = 2, height = 1.7)
