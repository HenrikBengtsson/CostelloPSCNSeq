#' @importFrom R.utils mprint cmdArg
#' @importFrom future %<-% %label% sessionDetails
#' @importFrom yaml yaml.load_file
#' @importFrom utils str
#' @importFrom aroma.seq findSamtools
#' @export
pscnseq <- function(what = c("mpileup", "sequenza"), session_details = interactive(), verbose = TRUE) {
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

  invisible(res)
}
