#' Calling the Parent-Specific Copy-Number Pipeline Step by Step
#'
#' @param what (character) The step to be performed; in order, one of
#' `"mpileup"`, `"sequenza"`, `"pscbs"`, or `"reports"`.
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
#' The pipeline looks for a YAML-formatted configuration file in the current
#' directory named \file{config.yml}.  An example of is:
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
#' This file is optional.  If not found, or incomplete, then the missing fields
#' have to be specified as command-line options, e.g.
#' `--organism=Homo_sapiens`.  If a field is specified in both, the
#' command-line version will take precedence.
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
pscnseq <- function(what = c("mpileup", "sequenza", "pscbs", "reports"), session_details = interactive(), verbose = TRUE) {
  assert <- NULL  ## To please R CMD check
  what <- match.arg(what)
  
  oopts <- options("R.filesets::onRemapping"="ignore")
  on.exit(oopts)

  if (session_details) {
    mprint(sessionDetails())
    mprint(findSamtools())
  }

  message("* Loading configuration")
  config <- cmdArg(config = "config.yml")
  config_data <- yaml.load_file(config)
  str(config_data)

  dataset <- cmdArg(dataset = config_data$dataset)
  organism <- cmdArg(organism = config_data$organism)
  chrs <- cmdArg(chrs = eval(parse(text = config_data$chromosomes)))
  samples <- cmdArg(samples = config_data$samples)
  fasta <- cmdArg(fasta = config_data$fasta)
  gcbase <- cmdArg(gcbase = config_data$gcbase)
  bam_pattern <- config_data$bam_pattern
  binSize <- cmdArg(binsize = eval(parse(text = config_data$binsize)))

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
