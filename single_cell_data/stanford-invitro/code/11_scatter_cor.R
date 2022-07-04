library(dplyr)
library(data.table)
library(BuenColors)
library(ggrepel)

pbmc <- fread("../output/PBMC-Day15-correlation.tsv")
s34_d5 <- fread("../output/Sample34-Day5-correlation.tsv")
s34_d5permuted <- fread("../output/Sample34-Day5-permuted.tsv")

s34_d7 <- fread("../output/Sample34-Day7-correlation.tsv")
mdf <- merge(pbmc, s34_d7, by = "gene")

# make a plot to call out genes
s34_d5$what <- "obs"
s34_d5permuted$what <- "perm"
s34_d5_all <- rbind(s34_d5, s34_d5permuted)
p1 <- ggplot(s34_d5_all, aes(x = rank, y = cor, color = what)) +
  geom_point(size = 0.5) +
  scale_color_manual(values = c("black", "lightgrey")) +
  pretty_plot(fontsize  = 7) + L_border() + 
  labs(x = "Rank ordered correlation", y = "Correlation with HHV-6B expression") +
  theme_void() + theme(legend.position = "none")
cowplot::ggsave2(p1, file = "../plots/correlation_rank.png", width = 2*2, height = 1.8*2, dpi = 300)


library(fgsea)
library(msigdbr)

gene_list <- s34_d7$cor; names(gene_list) <- s34_d7$gene

h_df_hallmark = msigdbr(species = "Homo sapiens") %>% dplyr::filter(gs_cat == "H")
pathway_list_hallmark = h_df_hallmark %>% split(x = .$gene_symbol, f = .$gs_name)
names(pathway_list_hallmark) <- gsub("HALLMARK_", "", names(pathway_list_hallmark))

fgseaRes_hallmark <- fgseaMultilevel(pathways = pathway_list_hallmark, 
                                     stats    = gene_list,
                                     minSize  = 15, eps = 0,
                                     maxSize  = 500)

p1 <- fgseaRes_hallmark[,c(1:7)] %>% arrange(pval) %>%
  ggplot(aes(x = NES, y = -log10(padj))) +
  geom_point() +
  pretty_plot(fontsize = 7) + L_border() +
  labs(x = "Normalized Enrichment Score")
cowplot::ggsave2(p1, file = "../plots/gsea_day5.pdf", width = 2, height = 1.8)


library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(org.Hs.eg.db)
organism = "org.Hs.eg.db"
gene_list <- s34_d7$cor; names(gene_list) <- s34_d7$gene

gse <- gseGO(geneList=gene_list, 
             ont ="ALL", 
             keyType = "SYMBOL", 
             nPerm = 100000, 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = org.Hs.eg.db, 
             pAdjustMethod = "none")
data.frame(gse)[,c(1:11)] %>% arrange(desc(enrichmentScore))


mdf <- mdf %>%
  mutate(label_rep = case_when(
    cor.x < -0.3 ~ gene,
    cor.x > 0.3  ~ gene,
    cor.y > 0.37  ~ gene,
    cor.y < -0.35  ~ gene,
    
    TRUE ~ ""
  ))
cor(mdf[,2], mdf[,3])
ggplot(mdf, aes(x = cor.x, y = cor.y, label = label_rep)) +
  geom_point(aes( color = label_rep != "")) + 
  
  geom_text_repel(max.overlaps = Inf, min.segment.length = 0, seed = 42, box.padding = 0.5,
                  fill = "white", xlim = c(-Inf, Inf), ylim = c(-Inf, Inf)) + 
  scale_color_manual(values = c("black", "red")) + 
  pretty_plot(fontsize = 20) + L_border()+  theme(legend.position = "none") +
  labs(x = "Cor. w/ HHV6 in PBMC infection (Day 15)", 
       y = "Cor. w/ HHV6 in PPD20034 (Day 7)")

