## Example:
## qcmd --exec Rscript 2.sequenza.R --config=config.yml --samples=sampleData/20161014_samplesforPSCN.txt

res <- CostelloPSCNSeq::pscnseq(what = "sequenza", verbose = TRUE)
print(res)
