## USAGE:
## qcmd --exec Rscript 3.pscbs.R --config=config.yml --samples=sampleData/20161014_samplesforPSCN.txt

res <- CostelloPSCNSeq::pscnseq(what = "pscbs", verbose = TRUE)
print(res)
