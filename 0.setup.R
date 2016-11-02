library("R.utils")
createLink(target = "/home/shared/cbc/annotationData")

path <- "bamData/CostelloP_2015-Exome,bwa,realigned,rmDups,recal/"
mkdirs(path)
createLink(link = file.path(path, "Homo_sapiens"), target = "/costellolab/jocostello/LG3/exomes_recal")

mkdirs("sampleData")
mkdirs("seqzData")

if (!file_test("-f", ".future.R")) {
cat(file = ".future.R", '
library("future.BatchJobs")
plan(list(
  samples     = tweak(batchjobs_torque, resources=list(vmem="1gb")),
  chromosomes = tweak(batchjobs_torque, resources=list(vmem="5gb"))
))
')
}
