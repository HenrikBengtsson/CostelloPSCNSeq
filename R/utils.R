seqToHumanReadable2 <- function(idxs, collapse=", ") {
  x <- suppressWarnings(as.numeric(idxs))
  nas <- is.na(x)
  idxs_numeric <- idxs[!nas]
  x_numeric <- seqToHumanReadable(idxs_numeric)
  x_non_numeric <- idxs[nas]
  paste(c(x_numeric, x_non_numeric), collapse = collapse)
}
