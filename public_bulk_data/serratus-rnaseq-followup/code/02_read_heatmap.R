library(data.table)
library(dplyr)
library(BuenColors)

samples <- gsub("_abundance.tsv", "", list.files("../data/tcell-followup/kallisto/"))
annotations <- fread('../../../hhv6-reference/HHV6b_only.index.txt', header = FALSE)
lapply(samples, function(x){
  (fread(paste0("../data/tcell-followup/kallisto/", x, "_abundance.tsv")) %>%
    mutate(gene = annotations[[2]], sample = x) %>%
    filter(gene != "DR1") %>% head(40) %>%
    mutate(count = round(est_counts)))[,c("gene", "count", "sample")]
}) %>% rbindlist() %>% data.frame() -> out_df

p1 <- ggplot(out_df, aes(y = sample, x = gene, label = count, fill = log2(count+1))) +
  geom_tile() + geom_text(size = 1) +
  pretty_plot(fontsize = 7) + pretty_plot(fontsize = 7) + L_border() +
  theme(legend.position = "none") +
  scale_fill_gradientn(colors = jdb_palette("solar_rojos"))
cowplot::ggsave2(p1, file = "../plots/read_heatmap.pdf", width = 7, height = 2)
