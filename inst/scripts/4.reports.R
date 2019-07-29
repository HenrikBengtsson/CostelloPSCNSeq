## USAGE:
## qcmd --exec Rscript 4.reports.R --config=config.yml --samples=sampleData/20161014_samplesforPSCN.txt

res <- CostelloPSCNSeq::pscnseq(what = "reports", verbose = TRUE)
print(res)
