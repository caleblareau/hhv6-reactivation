library(data.table)
library(dplyr)
library(Seurat)
library(annotables)
gene_vec <- grch38$symbol; names(gene_vec) <- grch38$ensgene

md <- fread("../data/fromNick/cell-annotation-metadata.txt", header = FALSE)
e11_bcs <- paste0(md %>% filter(V2 == "E1_1" & V3 == "Axi-R-15") %>% pull(V1)%>% substr(1, 16), "-2")
e12_bcs <- paste0(md %>% filter(V2 == "E1_2" & V3 == "Axi-R-15") %>% pull(V1) %>% substr(1, 16), "-3")

# Process h5s
channel <- "BBB7"
bbb7 <- Read10X_h5("../data/BBB7_raw_feature_bc_matrix.h5")
e11 <- Read10X_h5("../data/E1_1_raw_feature_bc_matrix.h5"); colnames(e11) <- gsub("-1", "-2", colnames(e11))
e12 <- Read10X_h5("../data/E1_2_raw_feature_bc_matrix.h5"); colnames(e12) <- gsub("-1", "-3", colnames(e12))
gex <- cbind(bbb7[,colSums(bbb7>0) >= 500],
             e11[,colSums(e11>0) >= 500 & colnames(e11) %in% e11_bcs],
             e12[,colSums(e12>0) >= 500 & colnames(e12) %in% e12_bcs])
dim(gex)
# Import kallisto dat
car <- rbind(
  fread(paste0("../data/new-data/features/BBB7_CAR.quantFinal.txt")) %>%
  mutate(bc = paste0(V1, "-1")),
  fread(paste0("../data/E1_1_CAR.quantFinal.txt")) %>%
    mutate(bc = paste0(V1, "-2")),
  fread(paste0("../data/E1_2_CAR.quantFinal.txt")) %>%
    mutate(bc = paste0(V1, "-3"))
)

hhv6 <- rbind(
  fread(paste0("../data/new-data/features/BBB7_HHV6b.quantFinal.txt")) %>%
    mutate(bc = paste0(V1, "-1")),
  fread(paste0("../data/E1_1_HHV6b.quantFinal.txt")) %>%
    mutate(bc = paste0(V1, "-2")),
  fread(paste0("../data/E1_2_HHV6b.quantFinal.txt")) %>%
    mutate(bc = paste0(V1, "-3"))
) %>% filter(!(V3 %in% c(6, 120))) 
  
# Summarize more readily
print(dim(hhv6))
hhv6p_cells <- unique(hhv6$bc)
print(length(hhv6p_cells))
hhv6p_count <- table(hhv6[["bc"]])
car_count <- table(car[["bc"]])

# Now subset for gene expression
gex2 <- gex[,colSums(gex) > 1000]
dim(gex2)
mdf <- data.frame(
  row.names = colnames(gex2),
  CAR = unname(car_count[colnames(gex2)]),
  hhv6count = unname(hhv6p_count[colnames(gex2)])
)[,c(2,4)]
mdf[is.na(mdf)] <- 0
mdf %>% filter(hhv6count.Freq > 0)

# Set up seurat analyses
so <- CreateSeuratObject(counts = gex2, meta.data = mdf)
so <- NormalizeData(so) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% 
  FindNeighbors() %>%  FindClusters() %>% RunUMAP(dims = 1:30)
so$HHV6log <- log1p(so$hhv6count.Freq)
so$carlog <- log1p(so$CAR.Freq)
so$channel <- substr(colnames(so), 18, 18)
DimPlot(so, group.by  = c("channel", "seurat_clusters"), label = TRUE)
FeaturePlot(so, features = c("HHV6log", "CAR.Freq", "CD69", "CD3D", "CD3E", "CD4", "CD8A", "CD8B", "CD14", "MS4A1", "NCAM1", "LAMP1"), sort = TRUE)

# Remove the junk (nk / mono)
monocyte_clusters <- c(7,8,13)
nkclusters <- c(6,11,12)
so2 <- so[,!c(so$seurat_clusters %in% c(monocyte_clusters, nkclusters))]
so2$any_hhv6 <- so2$hhv6count.Freq > 0
so2 <- NormalizeData(so2) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% 
  FindNeighbors() %>%  FindClusters() %>% RunUMAP(dims = 1:30)
FeaturePlot(so2, features = c("HHV6log", "hhv6count.Freq", "carlog", "CAR.Freq", "CD69", "CD3D", "CD3E", "CD4", "CD8A", "CD8B", "CD14", "MS4A1", "NCAM1", "GZMB"), sort = TRUE)
DimPlot(so2, group.by  = c("any_hhv6"), label = TRUE)

library(Nebulosa)
plot_density(so2, c("CD4", "carlog"))

