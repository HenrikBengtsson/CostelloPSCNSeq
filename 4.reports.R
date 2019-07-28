## USAGE:
## qcmd --exec Rscript 4.reports.R --config=config.yml --samples=sampleData/20161014_samplesforPSCN.txt

library(R.utils)
library(CostelloPSCNSeq)

if (!interactive()) mprint(future::sessionDetails())

mprintf("Script: 4.reports.R ...\n")

message("* Loading configuration")
config <- cmdArg(config = "config.yml")
config_data <- yaml::yaml.load_file(config)
str(config_data)

dataset <- cmdArg(dataset = config_data$dataset)
organism <- cmdArg(organism = config_data$organism)
chrs <- cmdArg(chrs = eval(parse(text = config_data$chromosomes)))
samples <- cmdArg(samples = config_data$samples)

reports <- pscnseq_reports(dataset, organism = organism, chrs = chrs, samples = samples, verbose = TRUE)
print(reports)

if (!interactive()) mprint(future::sessionDetails())

