## USAGE:
## qcmd --exec Rscript 2b.pscbs.R --config=config.yml --samples=sampleData/20161014_samplesforPSCN.txt

library("aroma.seq")
mprint(sessionDetails())
library("listenv")
source("R/pairs_from_samples.R")

mprintf("Script: 2.pscbs.R ...\n")


message("* Loading configuration")
config <- cmdArg(config = "config.yml")
config_data <- yaml::yaml.load_file(config)
str(config_data)

dataset <- cmdArg(dataset = config_data$dataset)
organism <- cmdArg(organism = config_data$organism)
chrs <- cmdArg(chrs = eval(parse(text = config_data$chromosomes)))
samples <- cmdArg(samples = config_data$samples)
binSize <- cmdArg(binsize = eval(parse(text = config_data$binsize)))


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
chrsTag <- sprintf("chrs=%s", seqToHumanReadable(chrs))
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
pathnamesS <- file.path(pathS, filenamesS)
fullnamesS <- gsub("[.]seqz(|[.]gz)$", "", filenamesS)

datasetD <- fullname(datasetS, binSizeTag, "tcn=2")
pathD <- file.path("pscbsData", datasetD, organism)
pathD <- Arguments$getWritablePath(pathD)
fullnamesD <- gsub(",chr=chr", ",chr=", fullnamesS)
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
      fit <- Reduce(PSCBS::append, fits)
      verbose && print(verbose, fit)
      fits <- NULL ## Not needed anymore
      
      ## Assert all chromosomes have been processed
      missing <- setdiff(chrs, getChromosomes(fit))
      if (length(missing) > 0) {
        throw(sprintf("Sanity check failure: Some chromosomes were not processed for sample %s: %s", sQuote(sample), paste(chrs, collapse=", ")))
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
print(fitList)

mprint(sessionDetails())
