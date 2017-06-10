# CostelloJ-PSCN-Seq

Parent-specific copy-number estimation pipeline.


## Requirements

This pipeline is implemented in [R] and requires the R package [aroma.seq].  To install the latter and all of the required dependencies, call the following from R:
```r
> source("http://callr.org/install#HenrikBengtsson/aroma.seq")
```
To run the pipeline on a computer cluster, the [future.batchtools] package is also needed, which can be installed directly from CRAN as:
```r
> install.packages("future.batchtools")
```


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


[R]: https://www.r-project.org/
[aroma.seq]: https://github.com/HenrikBengtsson/aroma.seq/
[future.batchtools]: https://cloud.r-project.org/package=future.batchtools
