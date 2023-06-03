theme_blank <- function(...) {
  ret <- theme_bw(...)
  ret$line <- element_blank()
  ret$rect <- element_blank()
  ret$strip.text <- element_blank()
  ret$axis.text <- element_blank()
  ret$plot.title <- element_blank()
  ret$axis.title <- element_blank()
  ret$plot.margin <- structure(c(0, 0, -1, -1), unit = "lines", valid.unit = 3L, class = "unit")
  ret
}