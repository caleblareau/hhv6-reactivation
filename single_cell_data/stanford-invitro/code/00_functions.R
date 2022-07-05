library(dplyr)
import_from_kite <- function(kite_df, barcodes, logme= TRUE){
  ncount <- kite_df %>% filter(V3 != "120") %>% group_by(V1) %>%
    summarize(count = n())
  vec <- ncount[[2]]; names(vec) <- ncount[[1]]
  try <- vec[barcodes]
  if(logme){
    out <- log1p(ifelse(is.na(try), 0, try))
  } else {
    ifelse(is.na(try), 0, try)
  }
}


import_from_kite_matrix <- function(kite_df, barcodes){
  gene_idx <- fread("../../../hhv6-reference/HHV6b_only.index.txt", header = FALSE)
  vec <- gene_idx[["V2"]]
  ncount <- kite_df %>% filter(V3 != "120") %>% group_by(V1,V3) %>%
    summarize(count = n()) %>% mutate(gene = vec[as.numeric(V3) + 1]) %>%
    data.table
  d <- (data.frame(dcast.data.table(ncount, V1 ~ gene, fill = 0, value.var = "count", fun.aggregate = sum)))
  rownames(d) <- d[[1]]
  mat <- data.matrix(d[,-c(1,2)])
  missing_barcodes <- barcodes[!(barcodes %in% rownames(mat))]
  if(length(missing_barcodes) > 0){
    zero_mat <- matrix(0, nrow = length(missing_barcodes), ncol = dim(mat)[2])
    rownames(zero_mat) <- missing_barcodes
    colnames(zero_mat) <- colnames(mat)
    out_mat <- rbind(zero_mat, mat[rownames(mat) %in% barcodes,])[barcodes,]
  } else {
    out_mat <- mat[barcodes,]
  }
  return(out_mat)
  
}