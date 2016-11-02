## Example:
## qcmd --exec Rscript 1b.sequenza.R --config=config.yml --samples=sampleData/20161014_samplesforPSCN.txt

library("aroma.seq")
mprint(sessionDetails())
options("R.filesets::onRemapping"="ignore")
library("listenv")
source("R/pairs_from_samples.R")

message("* Loading configuration")
config <- cmdArg(config = "config.yml")
config_data <- yaml::yaml.load_file(config)
str(config_data)

dataset <- cmdArg(dataset = config_data$dataset)
organism <- cmdArg(organism = config_data$organism)
chrs <- cmdArg(chrs = eval(parse(text = config_data$chromosomes)))
samples <- cmdArg(samples = config_data$samples)


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
fa <- FastaReferenceFile(config_data$fasta)
print(fa)
gc <- GcBaseFile(config_data$gcbase)
print(gc)

## IMPORTANT: Sequenza requires that chromosome names in GC file
## and the FASTA file (hg19.fa) need to be identical and in the
## same order.
## PS. It is ok that BAM files are in a different order.
stopifnot(isCompatibleWith(gc, fa))
stopifnot(all(getSeqNames(gc) == getSeqNames(fa)))


## - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Pileup data
## - - - - - - - - - - - - - - - - - - - - - - - - - - -
pathMP <- file.path("seqzData", fullname(dataset, "mpileup"), organism)
pathMP <- Arguments$getReadablePath(pathMP)
filenamesMP <- dir(path=pathMP, pattern="[.]mpileup(|[.]gz)$")

pathD <- file.path("seqzData", fullname(dataset, "seqz"), organism)
pathD <- Arguments$getWritablePath(pathD)


## - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Process
## - - - - - - - - - - - - - - - - - - - - - - - - - - -
if (interactive()) readline("Press ENTER to start processing of data: ")


## - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Sequenza
## - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Tumor-normal pairs from sample data
pairs <- pairs_from_samples(samples)
mstr(pairs)
npairs <- nrow(pairs)

chrTags <- sprintf("chr=chr%s", chrs)
seqzList <- listenv()
for (ii in seq_len(npairs)) {
  pair <- pairs[ii,]
  patient <- pair$Patient_ID
  tumor <- pair$T.A0
  normal <- pair$N.A0
  sampleName <- fullname(patient, paste(c(tumor, normal), collapse="_vs_"))
  if (is.na(tumor)) {
    message(sprintf("Sample #%d (%s) of %d: Skipping; non-existing tumor", ii, sampleName, npairs))
    next
  }
  if (is.na(normal)) {
    message(sprintf("Sample #%d (%s) of %d: Skipping; non-existing normal", ii, sampleName, npairs))
    next
  }
  stopifnot(nzchar(patient), nzchar(tumor), nzchar(normal),
            !is.na(patient), !is.na(tumor), !is.na(normal))

  patternD <- sprintf("%s,chr=chr.*[.]seqz(|[.]gz)$", sampleName)
  seqz <- SeqzFileSet$byPath(pathD, pattern=patternD)
  
  namesD <- getFullNames(seqz)
  chrTags_ii <- gsub(".*,chr=chr", "", namesD)
  todo <- setdiff(chrs, chrTags_ii)

  ## Already done
  if (length(todo) == 0) {
    seqzList[[ii]] <- seqz
    message(sprintf("Sample #%d (%s) %d: Already processed", ii, sampleName, npairs))
    next
  }

  message(sprintf("Sample #%d (%s) %d: Processing %d chromosomes %s",
          ii, sampleName, npairs, length(todo), hpaste(todo)))

  ## Future label
  label <- sprintf("sample_%s", ii)
  
  seqzList[[ii]] %<-% {
    patternT <- sprintf("^%s,%s,.*[.]mpileup(|[.]gz)$", patient, tumor)
    patternN <- sprintf("^%s,%s,.*[.]mpileup(|[.]gz)$", patient, normal)
    pusT <- MPileupFileSet$byPath(pathMP, pattern=patternT)
    pusN <- MPileupFileSet$byPath(pathMP, pattern=patternN)
    namesT <- getFullNames(pusT)
    namesN <- getFullNames(pusN)

    ## Pair up tumor and normal for each chromosome
    idxsT <- idxsN <- NULL
    for (cc in seq_along(todo)) {
      chr <- todo[cc]
      pattern <- sprintf(",chr=%s$", chr)
      idxT <- grep(pattern, namesT)
      idxN <- grep(pattern, namesN)
      if (length(idxT) == 0) {
        stop(sprintf("Sample #%d (%s) of %d: Missing tumor %s for given pattern: %s",
	             ii, sampleName, npairs, sQuote(tumor), sQuote(pattern)))
      } else if (length(idxN) == 0) {
        stop(sprintf("Sample #%d (%s) of %d: Missing normal %sor sample %s for given pattern: %s", sQuote(normal),
	             sQuote(patient), sQuote(pattern)))
      }
      stopifnot(length(idxT) == 1, length(idxN) == 1)
      idxsT <- c(idxsT, idxT)
      idxsN <- c(idxsN, idxN)
    }
    pusT <- as.list(pusT[idxsT])
    pusN <- as.list(pusN[idxsN])
    stopifnot(length(pusT) == length(pusN), length(pusN) == length(todo))
    ## Normal always comes before tumor!
    pus <- list(pusN, pusT)
    seqz <- pileup2seqz(pus, gc=gc, sampleName=sampleName,
                        dataset=dataset, organism=organism, verbose=-10)
    print(seqz)
    seqz
  } %label% label
} ## for (ii in ...)

print(Sys.time())

seqzList <- resolve(seqzList, value=TRUE)
print(Sys.time())

for (ii in seq_along(seqzList)) {
  try({
    print(list(ii=ii, seqz=seqzList[[ii]]))
  })
}

mprint(sessionDetails())
