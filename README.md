# CostelloJ-PSCN-Seq

Parent-specific copy number estimation pipeline.

## Setup (once)

1. Run `Rscript 0.setup.R` once. This will setup links to shared annotation data sets and lab data files on the TIPCC compute cluster.
2. Make sure `./config.yml` is correct.  It specify the default analysis settings.  The individual entries can be overridden by individual command-line options to the below `Rscript` calls.


## Data processing

The following scripts should be run in order:

* `Rscript 1.mpileup.R`
* `Rscript 2.sequenza.R`
* `Rscript 3.pscbs.R`
* `Rscript 4.reports.R`

You may adjust `./config.yml` when new files come in, or you can overridden the defaults individually via command-line options.
