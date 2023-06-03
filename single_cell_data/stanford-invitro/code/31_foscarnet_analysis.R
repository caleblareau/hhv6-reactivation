library(data.table)
library(dplyr)
library(Matrix)
library(Seurat)

# Define donors
ll <- c("D34C" ,"D34T")


hhv6files <- paste0("../data/foscarnet/",ll,"_HHV6b.quant.txt")
sapply(hhv6files, function(x){
  (fread(x) %>%
    group_by(V1) %>%
    summarize(count = n()) %>%
    filter(count >= 10) %>% dim())[1]
}) -> ses
data.frame(
  ses, 
  ll
)

lapply(ll, function(x){
  # Import vec
  mat <- Read10X_h5(paste0("../data/foscarnet/",x,"_rna.h5"))
  mat_filt <-mat[,colSums(mat >0) >= 1000]
  
  # Import hhv6
  fread(paste0("../data/foscarnet/",x,"_HHV6b.quant.txt")) %>%
    filter(!c(V3 %in% c(6, 120))) %>%
    group_by(V1) %>%
    summarize(count = n()) %>% mutate(barcode = paste0(V1, "-1"))  -> df
  
  hhv6vec <- df$count; names(hhv6vec) <- df$barcode
  cv <- hhv6vec[colnames(mat_filt)]
  cv <- ifelse(is.na(cv), 0, cv)
  data.frame(
    cell = colnames(mat_filt), 
    count = cv,
    exp = x
  )
}) %>% rbindlist() %>% data.frame() -> odf

library(BuenColors)
ggplot(odf, aes(x = count+1)) + 
  scale_x_log10(limits = c(0.5, 1000)) +
  scale_y_log10() +
  geom_histogram() + facet_wrap(~exp, scales = "free_y")  + labs(x = "# HHV6 reads")
odf %>%
  group_by(exp) %>%
  summarize(n_cells = n(), n_hhv6 = sum(count > 0), n_se = sum(count >= 10)) %>%
  mutate(pct_0hhv = (n_cells - n_hhv6)/n_cells*100) -> df

fisher.test(matrix(c(df$n_hhv6, df$n_cells-df$n_hhv6), nrow = 2)) %>% str()

odf %>%
  group_by(exp) %>%
  summarize(n_cells = n(), n_hhv6 = sum(count > 0), n_se = sum(count >= 10)) %>%
  mutate(pct_0hhv = (n_cells - n_hhv6)/n_cells*100) %>%
  ggplot(aes(x = exp, y = pct_0hhv)) + 
  geom_bar(stat = "identity", fill = "lightgrey", color = "black", width = 0.5) +
  pretty_plot(fontsize = 8) + scale_y_continuous(expand = c(0,0)) + L_border() +
  labs(x ="", y = "% cells with 0 HHV-6 UMIs") -> p1
p1
cowplot::ggsave2(p1, file = "../plots/bar_pct.pdf", width = 1.8, height = 1.8)
