# CostelloJ-PSCN-Seq

Parent-specific copy-number estimation pipeline.


## Requirements

### Required software
This pipeline is implemented in [R] and requires R packages [aroma.seq] and [sequenza].  To install these packages and all of their dependencies, call the following from R:
```r
> source("http://callr.org/install#sequenza")
> source("http://callr.org/install#HenrikBengtsson/aroma.seq")
```
In addition to the above R dependencies, the pipeline requires that [samtools] is on the `PATH`.


## Setup (once)

1. Run `Rscript 0.setup.R` once. This will setup links to shared annotation data sets and lab data files on the TIPCC compute cluster.

2. Make sure `./config.yml` is correct.  It specify the default analysis settings.  The individual entries can be overridden by individual command-line options to the below `Rscript` calls.

3. Configure parallel processing following the instructions in Section 'Configure parallel processing' below.


## Data processing

The following scripts should be run in order:

* `Rscript 1.mpileup.R`
* `Rscript 2.sequenza.R`
* `Rscript 3.pscbs.R`
* `Rscript 4.reports.R`

You may want to adjust [`./config.yml`](https://github.com/HenrikBengtsson/Costello-PSCN-Seq/blob/master/config.yml) to process other data set. Alternatively, you can specify another file that this default via command-line option `--config`, e.g. `Rscript 1.mpileup.R --config=config_set_a.yml`.


### Data processing via scheduler

To process the above four steps via the Torque/PBS scheduler, use:

```sh
$ qsub -d $(pwd) 1-4.all.pbs
```

This will in turn submit the corresponding PBS scripts `1.mpileup.pbs`, `2.sequenza.pbs`, `3.pscbs.pbs`, and `4.reports.pbs` to the scheduler.  Those PBS scripts "freeze" software versions to R 3.4.4 and samtools 1.3.1.



## Configure parallel processing

The pipeline supports both sequential and parallel processing on a large number of backends and compute resources.  By default the pipeline is configured to process the data sequentially on the current machine, but this can easily be changed to run in parallel, say, on a compute cluster.  In order not to clutter up the analysis scripts, these settings are preferably done in a separate `.future.R` (loaded automatically by the [future] framework) in the project root directory.

To process data via a TORQUE / PBS job scheduler using the [future.batchtools] package, try with the configuration that we use for our UCSF TIPCC cluster;
```
# Copy to project directory
$ cp .future-configs/batchtools/.future.R .

# Copy to home directory
$ cp .future-configs/batchtools/.batchtools.torque.tmpl ~

# Install the future.batchtools package
$ Rscript -e "install.packages('future.batchtools')"
```
These should be generic enough to also run on other TORQUE / PBS systems.  If not, see the [batchtools] package for how to configure the template file.

You can verify that it works by trying the following in the project directory:
```r
> library("future")
Using future plan:
plan(list(samples = tweak(batchtools_torque, label = "sample", 
    resources = list(vmem = "2gb")), chromosomes = tweak(batchtools_torque, 
    label = "chr", resources = list(vmem = "5gb"))))
```
This confirms that as soon as the [future] package is loaded, it will source the `.future.R` script which in turn will setup the parallel settings.  It is `.future.R` that reports on the future plan used.

Next, we can try to submit a job to the scheduler using these settings by:
```r
> x %<-% Sys.info()[["nodename"]]
```
In this step, future.batchtools will import the `.batchtools.torque.tmpl` file.  If it fails to locate that file, there will be an error.  If it succeeds, a batchtools job will be submitted to the job scheduler, cf. `qstat -u $USER`.

Finally, if we try to look at the value of `x`;
```r
> x
[1] "n17"
> 
```
it will block until the job is finished and then its value will be printed. Here we see that the job was running on compute node n17.


[R]: https://www.r-project.org/
[samtools]: http://www.htslib.org/
[aroma.seq]: https://github.com/HenrikBengtsson/aroma.seq/
[sequenza]: https://cran.r-project.org/package=sequenza
[batchtools]: https://cran.r-project.org/package=batchtools
[future]: https://cran.r-project.org/package=future
[future.batchtools]: https://cran.r-project.org/package=future.batchtools
