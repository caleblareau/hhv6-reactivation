library(data.table)
library(dplyr)
library(Seurat)

md <- fread("../data/cell-annotation-metadata.txt", header = FALSE)


process_channel <- function(channel){
  # Subset meta data for specific channel
  mdc <- md %>% filter(V2 %in% channel)
  patient_vec <- as.character(mdc$V3); names(patient_vec) <- substr(mdc$V1, 1, 18)
  
  # Import kallisto dat
  car <- fread(paste0("../data/",channel,"_CAR.quantFinal.txt")) %>%
    mutate(bc = paste0(V1, "-1"))
  hhv6 <- fread(paste0("../data/",channel,"_HHV6b.quantFinal.txt")) %>%
    filter(!(V3 %in% c(6, 120))) %>%
    mutate(bc = paste0(V1, "-1"))
  print(dim(hhv6))
  hhv6p_cells <- unique(hhv6$bc)
  hhv6p_count <- table(hhv6[["bc"]])
  
  # Import Gex data
  gex <- Read10X_h5(paste0("../data/",channel,"_raw_feature_bc_matrix.h5"))
  gex_filtered <- gex[, hhv6p_cells]
  car_seq_count <- table(car[["bc"]])[hhv6p_cells]
  
  # Import TCR data
  tcr <- fread(paste0("../data/",channel,"_TCR_tcr_all_contig_annotations.csv"), sep = ",")
  tcra <- tcr %>% filter(chain == "TRA") %>% group_by(barcode) %>% summarize(count = sum(umis))
  tcrb <- tcr %>% filter(chain == "TRB") %>% group_by(barcode) %>% summarize(count = sum(umis))
  tcravec <- tcra$count; names(tcravec) <- tcra$barcode
  tcrbvec <- tcrb$count; names(tcrbvec) <- tcrb$barcode
  
  # Compbine it al 
  data.frame(
    n_car = as.numeric(car_seq_count),
    n_hhv6 = as.numeric(hhv6p_count[hhv6p_cells]),
    t(data.matrix(gex_filtered[c("CD3D","CD3E", "CD4", "CD8A", "CD8B","ZAP70",  "LCK", "TNF", "GZMK", "GNLY", "KLRG1", "ZEB2", "NKG7"),])), # "CCL5", "MS4A1", "TREM1"
    nTCRA = tcravec[hhv6p_cells],
    nTCRB = tcrbvec[hhv6p_cells],
    nGExpUMIs = colSums(gex_filtered), 
    pctMito = colSums(gex_filtered[grepl("^MT", rownames(gex_filtered)),])/colSums(gex_filtered)*100,
    patient = patient_vec[hhv6p_cells],
    channel = channel
  ) %>% arrange(desc(nGExpUMIs))
}

rbind(process_channel("E1_1"), # 2-3 cells
process_channel(channel = "E1_2"), # 2-3 cells
process_channel(channel = "F1") # 1 cell
) %>% filter(n_hhv6 > 0) -> df

filt_df <- df %>% filter(!is.na(patient)) %>% arrange(desc(n_hhv6))
write.table(filt_df, file = "../output/HHV6-positive-cells.tsv", 
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
