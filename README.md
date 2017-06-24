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

We've updated to future.batchtools which looks for ~/.batchtools.torque.tmpl; just copy Henrik's:

```
$ cp /home/henrik/.batchtools.torque.tmpl ~
```

You can verify that it works but trying the following in the project directory:

```r
> library("future")
Using future plan:
plan(list(samples = tweak(batchtools_torque, label = "sample", 
    resources = list(vmem = "1gb")), chromosomes = tweak(batchtools_torque, 
    label = "chr", resources = list(vmem = "5gb"))))
> x %<-% Sys.info()[["nodename"]]
> x  
[1] "n17"
> 
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

You may want to adjust [`./config.yml`](https://github.com/HenrikBengtsson/Costello-PSCN-Seq/blob/master/config.yml) to process other data set. Alternatively, you can specify another file that this default via command-line option `--config`, e.g. `Rscript 1.mpileup.R --config=config_set_a.yml`.


[R]: https://www.r-project.org/
[aroma.seq]: https://github.com/HenrikBengtsson/aroma.seq/
[future.batchtools]: https://cran.r-project.org/package=future.batchtools
