library(dplyr)

# Function for reverse complement
revComp <- function(invec) {
  chartr(old = "AaTtGgCc", new = "TTAACCGG",  stringi::stri_reverse(invec))
}

# Set up method for barcode translation
bc_translate <- fread("../data/multiome-whitelists/RNA-737K-arc-v1.txt.gz", header = FALSE)[[1]]
names(bc_translate) <- fread("../data/multiome-whitelists/ATAC-737K-arc-v1.txt.gz", header = FALSE)[[1]]
