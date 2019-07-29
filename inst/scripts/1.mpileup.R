## Example:
## qcmd --exec Rscript 1.mpileup.R --config=config.yml --samples=sampleData/20161014_samplesforPSCN.txt

res <- CostelloPSCNSeq::pscnseq(what = "mpileup", verbose = TRUE)
print(res)



