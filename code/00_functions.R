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