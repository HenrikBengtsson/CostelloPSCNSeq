#' Produces "pileup" Output Files from Read Alignments (BAM Files)
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
#' @param verbose (logical) If TRUE, then verbose output is produced,
#' otherwise not.
#'
#' @return A [aroma.seq::MPileupFileSet].
#'
#' @seealso
#' This function uses `aroma.seq:mpileup()`.
#'
#' @importFrom R.utils Arguments mprintf mstr hpaste isGzipped less
#' @importFrom R.filesets getNames getTags readDataFrame fullname getFullNames setFullNamesTranslator getPathnames
#' @importFrom aroma.core isCompatibleWith
#' @importFrom aroma.seq directoryStructure<- findSamtools BamDataSet hasIndex MPileupFileSet mpileup MPileupFileSet getSeqNames
#' @importFrom utils str
#'
#' @export
pscnseq_mpileup <- function(dataset, organism, chrs, samples, fasta, gcbase, bam_pattern = NULL, verbose = FALSE) {
  assert <- NULL ## To please R CMD check
  
  verbose <- Arguments$getVerbose(verbose)
  
  message("* Assertions")
  assert %<-% {
    ver <- attr(findSamtools(), "version")
    print(ver)
    stopifnot(ver < "1.4")
  } %label% "samtools-version"
  print(assert)
  
  
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Sample data
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - -
  samples <- readDataFrame(samples, fill=TRUE)
  o <- order(samples$Patient_ID, samples$Sample_ID)
  samples <- samples[o,]
  str(samples)
  
  
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Annotation data
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - -
  message("* Loading annotation data files ...")
  fa <- FastaReferenceFile(fasta)
  print(fa)
  stopifnot(!isGzipped(fa))
  gc <- GcBaseFile(gcbase)
  print(gc)
  
  ## IMPORTANT: Sequenza requires that chromosome names in GC file
  ## and the FASTA file (hg19.fa) need to be identical and in the
  ## same order.
  ## PS. It is ok that BAM files are in a different order.
  stopifnot(isCompatibleWith(gc, fa))
  stopifnot(all(getSeqNames(gc) == getSeqNames(fa)))
  
  message("* Loading annotation data files ... DONE")
  
  
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Sequence read data
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - -
  message("* Loading all BAM files ...")
  
  if (is.null(bam_pattern)) bam_pattern <- ".bwa.realigned.rmDups(|.recal)(|.bam)$"
  path <- file.path("bamData", dataset, organism)
  bams <- BamDataSet$byPath(path, recursive=TRUE, pattern=bam_pattern)
  stopifnot(length(bams) > 0)
  bams <- bams[grep("old", getPathnames(bams), invert=TRUE)]
  stopifnot(length(bams) > 0)
  bams <- setFullNamesTranslator(bams, function(name, ...) {
    name <- gsub(".bwa.realigned.rmDups.recal.bam", "", name, fixed=TRUE)
    name <- gsub(".bwa.realigned.rmDups.bam", "", name, fixed=TRUE)
    name <- gsub(".bam", "", name, fixed=TRUE)
    name <- gsub("_merged", "", name, fixed=TRUE)
    name <- gsub("_[ACGT]{8}_L[0-9]{3}", "", name)
    name
  })
  print(bams)
  
  directoryStructure(bams) <- list(
    pattern="([^/]*)/([^/]*)/([^/]*)/([^/]*)/([^/]*)",
    replacement=c(rootpath="\\1", dataset="\\2", organism="\\3", sample="\\4,\\5")
  )
  
  
  ## Keep patients of interest
  names <- getNames(bams)
  str(names)
  tags <- sapply(bams, FUN=function(bam) getTags(bam)[1])
  str(tags)
  keep <- which(paste(names, tags, sep=",") %in% paste(samples$Patient_ID, samples$A0, sep=","))
  str(keep)
  stopifnot(length(keep) > 0)
  bams <- bams[keep]
  print(bams)
  stopifnot(length(bams) > 0)
  
  ## Assert that all BAMs have been indexed
  stopifnot(all(vapply(bams, FUN = hasIndex, FUN.VALUE = FALSE)))
  
  message("* Loading all BAM files ... DONE")
  
  
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Process
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (interactive()) readline("Press ENTER to start processing of data: ")
  
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Count nucleotides at every(!) genomic position
  ## using 'samtools mpileup'
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - -
  chrLabels <- sprintf("chr%s", chrs)
  
  ## Note that the generated *.mpileup files are very large.
  res <- mpileup(bams, fa=fa, chromosomes=chrLabels, verbose=less(verbose, 20))
  print(res)
  
  mps <- MPileupFileSet(res)
  mps
}

