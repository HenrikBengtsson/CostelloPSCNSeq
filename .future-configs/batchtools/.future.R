## Future toplogy configuration
library("future.batchtools")

plan(list(
  samples     = tweak(batchtools_torque, label = "sample", resources=list(vmem = "3gb")),
  chromosomes = tweak(batchtools_torque, label = "chr",    resources=list(vmem = "5gb"))
))

R.utils::mprintf("Using future plan:\n")
R.utils::mprint(attr(plan(), "call"))
