library(data.table)
library(stringr)
library(Seurat)
library(dplyr)
library(BuenColors)

# Import data
high_raw <- fread("../data/ID34high_HHV6b.kb.txt") %>% filter(V3 != 120) # 120 is the bustools set for the DR1 gene
high <- high_raw %>% group_by(V1, V2, V3) %>% summarize(count = n())
dim(high); dim(high_raw)

# Import transcript maps
txt <- data.frame(fread("../data/transcripts.txt", header = FALSE, sep = " "))
txt[as.numeric(names(tail(sort(table(high$V3)))))+1,] # the plus one is critical here for the 0-based indexing

# Now process with human scRNA
mat <- Read10X_h5("../data/ID34high_feature_bc_matrix.h5")
umi_vec <- colSums(mat)
ngenes_vec <- colSums(mat > 0)
ddf <- data.frame(table(paste0(high[[1]], "-1")))
ddf$ngenes <- ngenes_vec[as.character(ddf$Var1)]; ddf$ngenes <- ifelse(is.na(ddf$ngenes), 0, ddf$ngenes)
ddf$numis <- umi_vec[as.character(ddf$Var1)]; ddf$numis <- ifelse(is.na(ddf$numis), 0, ddf$numis)

ddf %>% arrange(desc(Freq)) %>% head(20)

# Now do the total viral count
viral_count <- ddf$Freq; names(viral_count) <- ddf$Var1

# Do some sort of filtering
mat_use <- mat[,colSums(mat>0)> 500]

# Do seurat things
rna_object <- CreateSeuratObject(counts = mat_use)
rna_object$pct.mt <- PercentageFeatureSet(rna_object, pattern = "MT-")
rna_object <- subset(rna_object, pct.mt < 10)
rna_object <- NormalizeData(rna_object)
rna_object <- ScaleData(rna_object)
rna_object <- FindVariableFeatures(rna_object)
rna_object <- RunPCA(rna_object)
rna_object$HHV6count <- ifelse(unname(is.na(viral_count[colnames(rna_object)])), 0, viral_count[colnames(rna_object)])

# Correlation comparison with permuted
sc <- cor(rna_object$HHV6count/rna_object$nCount_RNA,
          t(rna_object@assays$RNA@scale.data[rna_object@assays$RNA@var.features,]),
          method = "spearman")
data.frame(spearman_cor = sort(sc[1,], decreasing = TRUE)) %>% head(20)


set.seed(1)
scperm <- cor(sample(rna_object$HHV6count/rna_object$nCount_RNA),
          t(rna_object@assays$RNA@scale.data[rna_object@assays$RNA@var.features,]),
          method = "spearman")

plot_df <- data.frame(spearman_cor = c(sort(sc[1,], decreasing = TRUE), sort(scperm[1,], decreasing = TRUE)),
           rank = c(1:2000, 1:2000), 
           color = c(rep("observed", 2000), rep("permuted", 2000)))

ggplot(plot_df, aes(x = rank, y = spearman_cor, color = color)) +
  geom_point(size = 0.5) + pretty_plot() + scale_color_manual(values = c("firebrick", "grey")) + L_border() +
  labs(x = "2000 variable genes (Rank ordered)", color = "", y = "Spearman correlation with HHV6 #") +
  geom_hline(yintercept = 0, linetype = 2)


# Examine the "jackpot" cells
total_HHV6 <- sum(rna_object@meta.data$HHV6count)
total_RNAs <- sum(rna_object@meta.data$nCount_RNA)

viral_abundance_df <- data.frame(
  null = sort(total_HHV6 * (rna_object@meta.data$nCount_RNA / total_RNAs), decreasing = TRUE),
  observed = sort(rna_object@meta.data$HHV6count, decreasing = TRUE),
  rank = c(rep(1:dim(rna_object)[2],2))
)
viral_abundance_df %>% reshape2::melt(id.vars = "rank") %>%
  ggplot(aes(x = rank, y = value , color = variable)) + 
  geom_point(size = 0.5) + pretty_plot() + scale_color_manual(values = c("grey", "firebrick")) + L_border() +
  labs(x = "Rank ordered cells", color = "", y = "Number of HHV6 UMIs (log10)") +
  geom_hline(yintercept = 0, linetype = 2)+ scale_y_log10()

sum(rna_object@meta.data$HHV6count >= 1)
sum(rna_object@meta.data$HHV6count >= 2)
sum(rna_object@meta.data$HHV6count >= 5)

# Try classic DGE
rna_object$HHV6_on <- ifelse(rna_object$HHV6count>0, "yes", "no")
fm <- FindMarkers(rna_object, ident.1 = "yes", ident.2 = "no", group.by = "HHV6_on")
