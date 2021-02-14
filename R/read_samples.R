#' @importFrom R.filesets readDataFrame
read_samples <- function(pathname) {
  samples <- readDataFrame(pathname, fill=TRUE)
  str(samples)
  stopifnot(
    "Patient_ID" %in% colnames(samples),
    "Sample_ID"  %in% colnames(samples)
  )
  o <- order(samples$Patient_ID, samples$Sample_ID)
  samples <- samples[o,]
  samples
}
