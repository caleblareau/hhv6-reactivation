library(data.table)
library(dplyr)
library(viridis)
library(Seurat)
library(BuenColors)
source("00a_multiome_helper_functions.R")
source("00_functions.R")
x <- Read10X_h5("../../../../hhv6-large-data-files/invitro_cultures/ALLO98_ARC_v2.feature_mat.h5")
donorz <- "D98"
bcs <- fread(paste0("../data/multiome/98ARC_per_barcode_metrics.csv.gz")) %>%
  filter(is_cell == 1) %>%
  mutate(barcode = paste0(substr(barcode, 1, 16)))
cdf <- bcs[,c("barcode", "atac_fragments", "gex_exonic_umis", "gex_intronic_umis")]

# Import RNA
lapply(donorz, function(d){
  dt <- fread(paste0("../data/multiome/98GEX_HHV6b.quant.txt"))
  dt %>% filter(!(V3 %in% c(6, 120))) %>%
    filter(V1 %in% unname(bc_translate)) %>%
    mutate(barcode = paste0( V1)) %>% 
    group_by(barcode) %>% summarize(RNA_count = n()) %>%
    arrange(desc(RNA_count)) %>%
    mutate(Donor = d, rna_rank = 1:n())
}) %>% rbindlist() %>% data.frame() -> rna_df
colnames(rna_df) <- c("barcode", "RNA_HHV6_count", "Donor", "rna_rank")

crap_positions <- c(68576:69516)

# Append kite counts
hhv6_mat <- import_from_kite_matrix(unique(fread(paste0("../data/multiome/98GEX_HHV6b.quant.txt"))),rna_df$barcode)
                                   
hhv6_annotation <- fread("../../../hhv6-reference/HHV6b_expression_annotations.tsv", header = FALSE)
process_hhv6_annotations_to_pct <- function(mat, annotation_df){
  mat <-mat[,colnames(mat) %in% annotation_df[["V1"]]]
  rs <- rowSums(mat) 
  possible <- unique(annotation_df[[2]])
  mat_pct <- sapply(possible, function(x){
    genes_keep <- annotation_df[x == annotation_df[[2]],][[1]]
    rowSums(mat[,colnames(mat) %in% genes_keep])/rs
  })
  mat_pct
}
hhv6_pct_mat <- process_hhv6_annotations_to_pct(hhv6_mat, hhv6_annotation)
rna_df <- data.frame(rna_df, hhv6_pct_mat)

# Import ATAC
lapply(donorz, function(d){
  adt <- fread(paste0("../data/multiome/98ATAC_HHV6b.bed"))
  adt$barcode <- unname(bc_translate[revComp(as.character(adt$V4))])
  
  adt[complete.cases(adt),] %>%
    mutate(barcode = paste0(barcode)) %>%
    filter(!(V2 %in% crap_positions)) %>%
    filter(!(V3 %in% crap_positions)) %>%
    group_by(barcode) %>% summarize(ATAC_HHV6_count = n()) %>%
    arrange(desc(ATAC_HHV6_count)) %>%
    mutate(Donor = d, atac_rank = 1:n())
}) %>% rbindlist() %>% data.frame() -> atac_df

# Now merge them and filter for called cells
mdf <- merge(atac_df, rna_df, by = c("barcode", "Donor"), all = TRUE)
mdf <- mdf[mdf$barcode %in% cdf$barcode,]

mdf2 <- merge(cdf, mdf, all.x = TRUE, by.x ="barcode", by.y = "barcode")
mdf2[is.na(mdf2)] <- 0

mdf2 %>% arrange(desc(RNA_HHV6_count))
ps <- ggplot(shuf(mdf2), aes(x = ATAC_HHV6_count + 1, y = RNA_HHV6_count + 1)) + 
  geom_point(size = 0.5) + scale_x_log10() + scale_y_log10() +
  pretty_plot(fontsize = 7) + labs( x= "log10(#ATAC Frags. + 1)", y =  "log10(#RNA UMIs + 1)") + L_border()
cowplot::ggsave2(ps, file = paste0("../plots/multiome_scatter_both_abundance.pdf"), width = 1.8, height = 1.8)

mdf3 <- mdf2 %>% filter(ATAC_HHV6_count >0 | RNA_HHV6_count>0)
cor.test(log1p(mdf2$ATAC_HHV6_count), log1p(mdf2$RNA_HHV6_count), method = "pearson") %>%
  str()
  

saveRDS(mdf2, file = "../output/all_meta_data_multiome.rds")


px <- ggplot(mdf2 %>% filter(RNA_HHV6_count >= 10), aes(x = late*100, y = early*100, color = immediate_early)) +
  geom_point(size = 0.5) + 
  scale_color_viridis(limits = c(0,1)) + scale_y_continuous(limits = c(0,100)) + 
  scale_x_continuous(limits = c(0,100)) + 
  pretty_plot(fontsize = 8) + L_border() + labs(x = "%HHV6 UMIs - Late", y = "%HHV6 UMIs - Early") +
  theme(legend.position = "none") 

pa <- ggplot(mdf2 %>% filter(RNA_HHV6_count >= 10), aes(x = late*100, y = early*100, color = log1p(ATAC_HHV6_count))) +
  geom_point(size = 0.5) + 
  scale_color_viridis() + scale_y_continuous(limits = c(0,100)) + 
  scale_x_continuous(limits = c(0,100)) + 
  pretty_plot(fontsize = 8) + L_border() + labs(x = "%HHV6 UMIs - Late", y = "%HHV6 UMIs - Early") +
  theme(legend.position = "none") 


cowplot::ggsave2(px, file = paste0("../plots/multiome_scatter_hhv6signature.pdf"), width = 1.8, height = 1.8)
cowplot::ggsave2(pa, file = paste0("../plots/multiome_scatter_atac_abundance.pdf"), width = 1.8, height = 1.8)


mdf2$logatac <- log1p(mdf2$ATAC_HHV6_count)


for_cor <- data.matrix(data.frame(mdf2 %>% filter(RNA_HHV6_count >= 10))[,c("immediate_early","intermediate_early_early","early","late", "logatac")]) 
cor(for_cor)
cor.test(for_cor[,"logatac"], for_cor[,"intermediate_early_early"])
cor.test(for_cor[,"logatac"], for_cor[,"immediate_early"])
cor.test(for_cor[,"logatac"], for_cor[,"late"])
cor.test(for_cor[,"logatac"], for_cor[,"early"])

qplot(for_cor[,"logatac"], for_cor[,"early"])

data.frame(value = cor(for_cor)[5,-5]) %>%
  mutate(variable = rownames(.)) %>%
  ggplot(aes(x = variable, y = value)) + 
  geom_bar(stat = "identity", position='dodge', color = "black", fill = "lightgrey", width = 0.5) +
  pretty_plot(fontsize = 7) + L_border() + geom_hline(yintercept = 0,) + 
  labs(x = "Gene signature", y = "Pearson correlation", fill = "") -> p1

cowplot::ggsave2(p1, file = "../plots/bar_gene_sig_hhv6_correlation.pdf", width = 2, height = 1.4)

