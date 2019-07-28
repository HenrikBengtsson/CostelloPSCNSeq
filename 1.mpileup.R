## Example:
## qcmd --exec Rscript 1.mpileup.R --config=config.yml --samples=sampleData/20161014_samplesforPSCN.txt

library(R.utils)
library(CostelloPSCNSeq)

options("R.filesets::onRemapping"="ignore")

if (!interactive()) {
  mprint(future::sessionDetails())
  mprint(findSamtools())
}

message("* Loading configuration")
config <- cmdArg(config = "config.yml")
config_data <- yaml::yaml.load_file(config)
str(config_data)

message("* Assertions")
assert %<-% {
  ver <- attr(findSamtools(), "version")
  print(ver)
  stopifnot(ver < "1.4")
} %label% "samtools-version"
print(assert)

message("* Loading configuration")
config <- cmdArg(config = "config.yml")
config_data <- yaml::yaml.load_file(config)
str(config_data)

dataset <- cmdArg(dataset = config_data$dataset)
organism <- cmdArg(organism = config_data$organism)
chrs <- cmdArg(chrs = eval(parse(text = config_data$chromosomes)))
samples <- cmdArg(samples = config_data$samples)
fasta <- cmdArg(fasta = config_data$fasta)
gcbase <- cmdArg(gcbase = config_data$gcbase)
bam_pattern <- config_data$bam_pattern

mps <- pscnseq_mpileup(dataset, organism = organism, chrs = chrs, samples = samples, fasta = fasta, gcbase = gcbase, bam_pattern = bam_pattern, verbose = TRUE)
print(mps)

if (!interactive()) {
  mprint(future::sessionDetails())
  mprint(findSamtools())
}



