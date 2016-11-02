# CostelloJ-PSCN-Seq

Parent-specific copy number estimation pipeline.

The following scripts should be run in order:

* `Rscript 1.mpileup.R`
* `Rscript 2.sequenza.R`
* `Rscript 3.pscbs.R`
* `Rscript 4.reports.R`

The all use `./config.yml` for their default setup,
which can be overridden by individual command-line options.
