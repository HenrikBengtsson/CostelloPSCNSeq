#' @importFrom future %<-% %label% resolve
#' @importFrom listenv listenv
#' @importFrom R.utils Arguments isFile saveObject mprintf mstr hpaste seqToHumanReadable cat enter enterf exit
#' @importFrom utils file_test
#' @importFrom R.filesets readDataFrame sortBy
#' @importFrom aroma.seq SeqzFileSet
#' @importFrom PSCBS segmentByPairedPSCBS normalizeTotalCNs
#'
#' @export
pscnseq_pscbs <- function(dataset, organism, chrs, samples, binSize) {
  `%<-%` <- future::`%<-%`
  `%label%` <- future::`%label%`
  resolve <- future::resolve
  listenv <- listenv::listenv
  readDataFrame <- R.filesets::readDataFrame
  sortBy <- R.filesets::sortBy
  library(R.utils)
  library(utils)
  library(aroma.seq)
  library(PSCBS)
  
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Sample data
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (samples != "*") {
    ## Tumor-normal pairs from sample data
    data <- readDataFrame(samples, fill=TRUE)
    pairs <- pairs_from_samples(data)
    mstr(pairs)
  }
  
  
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Overview
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - -
  dataset <- config_data$dataset
  organism <- config_data$organism
  chrs[chrs == "X"] <- 23
  chrs[chrs == "Y"] <- 24
  chrs[chrs == "M"] <- 25
  chrsTag <- sprintf("chrs=%s", seqToHumanReadable(chrs, collapse = ","))
  binSizeTag <- sprintf("%gkb", binSize/1000)
  
  mprintf("Dataset: %s\n", dataset)
  mprintf("Organism: %s\n", organism)
  mprintf("Chromosomes: %s\n", seqToHumanReadable(chrs))
  mprintf("Bin size: %s (%d bp)\n", binSizeTag, binSize)
  
  
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Identify files to be processed
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Identify availble Sequenza files
  datasetS <- paste(c(dataset, "seqz"), collapse = ",")
  pathS <- file.path("seqzData", datasetS, organism)
  pathS <- Arguments$getReadablePath(pathS)
  filenamesS <- dir(path=pathS, pattern="[.]seqz(|[.]gz)$")
  ## Sanity check
  if (length(filenamesS) == 0) {
    stop("Failed to locate any *.seqz(.gz) files: ", sQuote(pathS))
  }
  pathnamesS <- file.path(pathS, filenamesS)
  fullnamesS <- gsub("[.]seqz(|[.]gz)$", "", filenamesS)
  
  datasetD <- fullname(datasetS, binSizeTag, "tcn=2")
  pathD <- file.path("pscbsData", datasetD, organism)
  pathD <- Arguments$getWritablePath(pathD)
  fullnamesS2 <- gsub("chrX", "chr23", fullnamesS)
  fullnamesS2 <- gsub("chrY", "chr24", fullnamesS2)
  fullnamesS2 <- gsub("chrM", "chr25", fullnamesS2)
  fullnamesD <- gsub(",chr=chr", ",chr=", fullnamesS2)
  filenamesD <- sprintf("%s,PairedPSCBS.xdr", fullnamesD)
  pathnamesD <- file.path(pathD, filenamesD)
  
  todo <- which(!file_test("-f", pathnamesD))
  if (length(todo) == 0) {
    mprintf("All %d Sequenza files have been PSCBS segmented: %s", length(pathnamesS), pathD)
  }
  
  ## Keep non-processed files
  if (length(todo) != length(pathnamesD)) {
    pathnamesS <- pathnamesS[todo]
    pathnamesD <- pathnamesD[todo]
    fullnamesS <- fullnamesS[todo]
  }
  
  mprintf("Sequenza files located: [n=%d] %s\n", length(fullnamesS), hpaste(sQuote(fullnamesS)))
  all_samples <- unique(gsub(",chr=.*", "", fullnamesS))
  
  
  if (FALSE && samples != "*") {
    stop("Not yet implemented.")
    sample_pairs <- sprintf("%s,%s_vs_%s", pairs$Patient_ID, pairs$T.A0, pairs$N.A0)
    keep <- (all_samples %in% sample_pairs)
    all_samples <- all_samples[keep]
    pathnamesS <- pathnamesS[keep]
    pathnamesD <- pathnamesD[keep]
    fullnamesS <- fullnamesS[keep]
  
    todo <- !(sample_pairs %in% unique(all_samples))
  }
  
  mprintf("Sequenza files to process: [n=%d] %s\n", length(fullnamesS), hpaste(sQuote(fullnamesS)))
  mprintf("Tumor-normal samples: [n=%d] %s\n", length(all_samples), hpaste(sQuote(all_samples)))
  
  
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Process
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (interactive()) readline("Press ENTER to start processing of data: ")
  
  fitList <- listenv()
  for (ii in seq_along(all_samples)) {
    sample <- all_samples[ii]
    verbose && enterf(verbose, "Sample %d ('%s') of %d ...\n", ii, sample, length(all_samples))
  
    ## Future label
    label <- sprintf("sample_%s", ii)
  
    filenameF <- sprintf("%s,PairedPSCBS.rds", fullname(sample, chrsTag))
    pathnameF <- file.path(pathD, filenameF)
    if (isFile(pathnameF)) {
      verbose && cat(verbose, "Already processed. Skipping")
      fitList[[ii]] <- pathnameF
    } else {
      fitList[[ii]] %<-% {
        pattern <- sprintf("%s,chr=.*[.]seqz(|[.]gz)$", sample)
        seqz <- SeqzFileSet$byPath(pathS, pattern=pattern, depth=1L)
        seqz <- sortBy(seqz, by="mixeddecimal") # chr1, chr2, ..., chr10, ...
        verbose && print(verbose, seqz)
        verbose && printf(verbose, "Tumor-normal samples: [n=%d] %s\n",
                          length(seqz), hpaste(sQuote(getFullNames(seqz))))
    
        fits <- segmentByPairedPSCBS(seqz, binSize=binSize, verbose=TRUE)
        verbose && cat(verbose, "PSCBS segmentation:")
        verbose && print(verbose, fits)
      
        ## Since 'fits' consists of futures, we'll collect here
        fits <- resolve(fits, value = TRUE)
        ## Coerce to list
        fits <- as.list(fits)
        fit <- do.call(c, fits)
        verbose && print(verbose, fit)
        fits <- NULL ## Not needed anymore
        
        ## Assert all chromosomes have been processed
        missing <- setdiff(chrs, getChromosomes(fit))
        if (length(missing) > 0) {
          throw(sprintf("Sanity check failure: Some chromosomes were not processed for sample %s: %s", sQuote(sample), paste(missing, collapse=", ")))
        }
  
        verbose && cat(verbose, "PSCBS total copy-number normalization:")
        fit <- normalizeTotalCNs(fit)
        verbose && print(verbose, fit)
        
        ## Save to file (atomically)
        saveObject(fit, pathnameF)
        
        pathnameF
      } %label% label
    } ## if (isFile(pathnameF))
  
    ## Have at most 10 jobs on the cluster at any time
    if (ii %% 10 == 0) resolve(fitList, value = TRUE)
    verbose && exit(verbose)
  } ## for (ii ...)
  
  ## Resolve futures and gather their values
  fitList <- resolve(fitList, value = TRUE)
  
  ## Vectorize
  fitList <- unlist(fitList)
  
  fitList
}