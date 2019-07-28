## Example:
## qcmd --exec Rscript 1.mpileup.R --config=config.yml --samples=sampleData/20161014_samplesforPSCN.txt

library("aroma.seq")
source("R/pairs_from_samples.R")
source("R/pscnseq_mpileup.R")
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

mps <- pscnseq_mpileup(dataset, organism = organism, chrs = chrs, samples = samples)

if (!interactive()) {
  mprint(future::sessionDetails())
  mprint(findSamtools())
}



