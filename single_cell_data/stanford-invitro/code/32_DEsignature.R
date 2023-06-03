library(data.table)
library(dplyr)
library(Matrix)
library(BuenColors)
library(Seurat)

# Annotate stuff
D34C <- Read10X_h5("../data/foscarnet/D34C_rna.h5"); colnames(D34C) <- paste0("D34C-", colnames(D34C))
D34T <-Read10X_h5("../data/foscarnet/D34T_rna.h5"); colnames(D34T) <- paste0("D34T-", colnames(D34T))

common <- intersect(rownames(D34C), rownames(D34T))
d <- cbind(D34C[common,],D34T[common,])

# Filter object
so <- CreateSeuratObject(d)
so <- subset(so, nFeature_RNA > 1000 )
so$channel <- substr(colnames(so), 1, 4)
so <- so %>% NormalizeData() %>%
  FindVariableFeatures %>%
  ScaleData %>% RunPCA %>% RunUMAP(dims = 1:30) %>%
  FindNeighbors(dims = 1:20)

DimPlot(so, group.by = c("channel"), shuffle = TRUE) +
  theme_void() + ggtitle("") + 
  scale_color_manual(values = c("dodgerblue3", "firebrick")) +
  theme(legend.position = "none") -> p1
#cowplot::ggsave2(p1, file = "../plots/umap_foscarnet.png", width = 3, height = 3, dpi = 400)

d34_diff <- FindMarkers(so, group.by = "channel", ident.1 = "D34C", ident.2 = "D34T", logfc.threshold = 0.05)
head(d34_diff)
dim(d34_diff)
p2 <- ggplot(d34_diff, aes(x = avg_log2FC, y = -log10(p_val_adj)))+
  geom_point() +
  scale_x_continuous(limits = c(-1.5, 1.5)) +
  pretty_plot(fontsize = 8) + L_border()
cowplot::ggsave2(p2, file = "../plots/foscarnet_dge.pdf", width = 1.7, height = 1.7)

FeaturePlot(so, c("CD8A", "CD8B", "CD4", "CD3E", "nFeature_RNA"), sort.cell = TRUE)
