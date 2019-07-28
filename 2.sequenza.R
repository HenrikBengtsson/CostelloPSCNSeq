## Example:
## qcmd --exec Rscript 2.sequenza.R --config=config.yml --samples=sampleData/20161014_samplesforPSCN.txt

## https://github.com/HenrikBengtsson/Costello-PSCN-Seq/issues/25
stopifnot(packageVersion("sequenza") <= "2.1.2")

library(R.utils)
library(CostelloPSCNSeq)

options("R.filesets::onRemapping"="ignore")

if (!interactive()) mprint(future::sessionDetails())

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

seqzList <- pscnseq_sequenza(dataset, organism = organism, chrs = chrs, samples = samples, fasta = fasta, gcbase = gcbase, verbose = TRUE)
print(seqzList)

if (!interactive()) mprint(future::sessionDetails())
