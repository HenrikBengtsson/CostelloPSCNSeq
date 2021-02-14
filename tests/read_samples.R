library(CostelloPSCNSeq)

read_samples <- CostelloPSCNSeq:::read_samples
pairs_from_samples <- CostelloPSCNSeq:::pairs_from_samples

pathname <- system.file(package = "CostelloPSCNSeq", "samples.tsv", mustWork = TRUE)
samples <- read_samples(pathname)
print(samples)
stopifnot(nrow(samples) == 3)

pairs <- pairs_from_samples(samples)
print(pairs)
stopifnot(
  unique(pairs$PatientID) == unique(samples$PatientID),
  nrow(pairs) == 2
)