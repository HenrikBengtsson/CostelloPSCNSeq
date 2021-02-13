#' Calling the Parent-Specific Copy-Number Pipeline Step by Step
#'
#' @param what (character) The step to be performed; in order, one of
#' `"mpileup"`, `"sequenza"`, `"pscbs"`, or `"reports"`.
#'
#' @param dataset (character) The name of the dataset as on file.
#'
#' @param organism (character) The name of the organism as on file.
#'
#' @param chrs (character vector) The name of the chromosomes to be processed,
#' e.g. `c("1", "2", "X")`.
#'
#' @param samples (character) Pathname to a tab-delimited sample specification
#' file, typically named \file{*.tsv}, e.g. \file{samples.tsv}.
#'
#' @param fasta (character) The pathname to the FASTA reference file,
#' typically named \file{*.fa} or \file{*.fasta}, e.g. \file{hg19.fa}.
#'
#' @param gcbase (character) The pathname to the FASTA reference file,
#' typically named \file{*.txt.gz}, e.g. \file{hg19.gc50Base.txt.gz}.
#'
#' @param bam_pattern (character; optional) Regular expression to identify
#' subset of BAM files to be processed.  If NULL (default), then BAM files
#' matching `.bwa.realigned.rmDups(|.recal)(|.bam)$` are included.
#'
#' @param binSize (integer or numeric) The bin size (in basepairs) used for
#' binning reads into bins that then are passed to the segmentation method.
#'
#' @param config (character) Pathname to YAML configuration file.
#' If NULL, then the configuration file is skipped.
#'
#' @param session_details (logical) If TRUE, session details are reported
#' before starting the processing and after it completed.
#'
#' @param verbose (logical) If TRUE, then verbose output is produced,
#' otherwise not.
#'
#' @return Returns what the called `pscnseq_nnn()` function returns, i.e.
#' [pscnseq_mpileup()], [pscnseq_sequenza()], [pscnseq_pscbs()], or
#' [pscnseq_reports()].
#'
#' @section Configuration File:
#' The default arguments can be set in an YAML-formatted configuration file
#' as given by argument `config`.  The default is to look for a file named
#' \file{config.yml} in the current directory.  To skip this file, specify
#' `config = NULL`.  An example of such a file is:
#'
#' ```yaml
#' organism: Homo_sapiens
#' chromosomes: c(1:22, "X", "Y", "M")
#' fasta: annotationData/organisms/Homo_sapiens/GRCh37,hg19/UCSC/hg19.fa
#' gcbase: annotationData/organisms/Homo_sapiens/GRCh37,hg19/UCSC/hg19.gc50Base.txt.gz
#' dataset: CostelloP_2015-Exome,bwa,realigned,rmDups,recal
#' binsize: 100e3
#' samples: sampleData/samples.tsv
#' ```
#'
#' @section Specifying arguments via command-line options:
#' The arguments can be overridden by command-line options, e.g.
#' `--organism=Homo_sapiens` will take precedence of argument `organism`,
#' which in turn will take precedent of what is specified in the configuration
#' file.
#'
#' @section How to call pipeline from the command line:
#' Below is how you could run the pipeline step by step:
#'
#' ```sh
#' Rscript -e "CostelloPSCNSeq::pscnseq(what='mpileup', verbose=TRUE)"  # ~25 min
#' Rscript -e "CostelloPSCNSeq::pscnseq(what='sequenza', verbose=TRUE)" # ~60 min
#' Rscript -e "CostelloPSCNSeq::pscnseq(what='pscbs', verbose=TRUE)"    #  ~5 min
#' Rscript -e "CostelloPSCNSeq::pscnseq(what='reports', verbose=TRUE)"  #  ~2 min
#' ```
#'
#' @importFrom R.utils mprint cmdArg
#' @importFrom future %<-% %label% sessionDetails
#' @importFrom yaml yaml.load_file
#' @importFrom utils str
#' @importFrom aroma.seq findSamtools
#' @export
pscnseq <- function(what = c("mpileup", "sequenza", "pscbs", "reports"), dataset = NULL, organism = NULL, chrs = NULL, samples = NULL, fasta = NULL, gcbase = NULL, bam_pattern = NULL, binSize = NULL, config = "config.yml", session_details = interactive(), verbose = TRUE) {
  assert <- NULL  ## To please R CMD check
  what <- match.arg(what)

  verbose <- Arguments$getVerbose(verbose)

  if (verbose) {
    enter(verbose, "pscnseq() ...")
    on.exit(exit(verbose))
  }

  oopts <- options("R.filesets::onRemapping"="ignore")
  on.exit(oopts, add = TRUE)

  if (session_details) {
    mprint(sessionDetails())
    mprint(findSamtools())
  }

  # Arguments
  args <- list(
    dataset = dataset,
    organism = organism,
    chrs = chrs,
    samples = samples,
    fasta = fasta,
    gcbase = gcbase,
    bam_pattern = bam_pattern,
    binSize = binSize
  )
  
  config <- cmdArg(config = config)
  if (!is.null(config)) {
    message("* Loading configuration")  
    config_data <- yaml.load_file(config)
  } else {
    config_data <- list()
  }
  str(config_data)

  for (name in names(args)) {
    ## Priority 1: Use command-line argument, if set
    value <- cmdArg(name = name, default = NULL)
    if (!is.null(value)) {
      args[[name]] <- value
      next
    }

    ## Priority 2: Use function argument, if set
    if (!is.null(args[[name]])) next
    
    ## Priority 3: Use configure-file settings, if set
    value <- config_data[[name]]
    if (!is.null(value)) {
      args[[name]] <- value
      next
    }
  }

  verbose && cat(verbose, "Arguments:")
  verbose && str(verbose, args)

  res <- list()
  if (what == "mpileup") {
    message("* Assertions")
    assert %<-% {
      ver <- attr(findSamtools(), "version")
      print(ver)
      stopifnot(ver < "1.4")
    } %label% "samtools-version"
    print(assert)
    
    mps <- pscnseq_mpileup(dataset, organism = organism, chrs = chrs, samples = samples, fasta = fasta, gcbase = gcbase, bam_pattern = bam_pattern, verbose = verbose)
    print(mps)
    res[[what]] <- mps
  } else if (what == "sequenza") {
    ## https://github.com/HenrikBengtsson/Costello-PSCN-Seq/issues/25
    stopifnot(packageVersion("sequenza") <= "2.1.2")

    seqzList <- pscnseq_sequenza(dataset, organism = organism, chrs = chrs, samples = samples, fasta = fasta, gcbase = gcbase, verbose = verbose)
    print(seqzList)
    res[[what]] <- seqzList
  } else if (what == "pscbs") {
    fitList <- pscnseq_pscbs(dataset, organism = organism, chrs = chrs, samples = samples, binSize = binSize, verbose = TRUE)
    print(fitList)
    res[[what]] <- fitList
  } else if (what == "reports") {
    reports <- pscnseq_reports(dataset, organism = organism, chrs = chrs, samples = samples, verbose = TRUE)
    print(reports)
    res[[what]] <- reports
  }
  
  if (session_details) {
    mprint(sessionDetails())
  }

  res
}
