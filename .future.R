## Future toplogy configuration
library("future.BatchJobs")

plan(list(
  samples     = tweak(batchjobs_torque, label = "sample", resources=list(vmem="1gb")),
  chromosomes = tweak(batchjobs_torque, label = "chr",    resources=list(vmem="5gb"))
))

R.utils::mprintf("Using future plan:\n")
R.utils::mprint(attr(plan(), "call"))


