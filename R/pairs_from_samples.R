pairs_from_samples <- function(samples) {
  mixedorder <- gtools::mixedorder
  
  patients <- sort(unique(samples$Patient_ID))
##  mprintf("Patients: [%d] %s\n", length(patients), hpaste(sQuote(patients)))
  
  pairs <- NULL
  
  for (ii in seq_along(patients)) {
    patient <- patients[[ii]]
##    mprintf("Patient #%d (%s) ...\n", ii, sQuote(patient))
    dataII <- samples[samples$Patient_ID == patient,]
  
    ## Sanity check on assay kits
    kits <- sort(unique(dataII$Kit))
##    mprintf("Assay kit(s) used: [%d] %s\n", length(kits), hpaste(sQuote(kits)))
    if (length(kits) > 1) {
      warning(sprintf("More than one kit was used for patient %s: %s", sQuote(patient), hpaste(sQuote(kits))))
    }
  
    normal <- (dataII$Sample_ID == "Normal")
    dataN <- dataII[normal,]
    dataT <- dataII[!normal,]
  
    if (nrow(dataN) == 0) {
      warning(sprintf("Patient %s has no normal sample", sQuote(patient)))
    } else if (nrow(dataN) > 1) {
      warning(sprintf("Patient %s has more than one normal sample", sQuote(patient)))
    } else {
      dataP <- data.frame(Patient_ID=patient, T=dataT[,c("Sample_ID", "A0", "SF")], N=dataN[,c("Sample_ID", "A0", "SF")], stringsAsFactors=FALSE)
      npairs <- nrow(dataP)
      labels <- sprintf("%s,%s_vs_%s", patient, dataP$T.Sample_ID, dataP$N.Sample_ID)
##      mprintf("Tumor-normal pairs: [%d] %s\n", npairs, hpaste(sQuote(labels)))
    }
  
    pairs <- rbind(pairs, dataP)
  
##    mprintf("Patient #%d (%s) ... DONE\n", ii, sQuote(patient))
  } ## for (ii ...)
  rownames(pairs) <- NULL
  
  ## FIXME: For some reason we might end up with duplicated pairs
  pairs <- unique(pairs)
  
  ## Sort in human readable order
  o <- order(pairs$Patient_ID, pairs$T.Sample_ID, pairs$N.Sample_ID)
  pairs <- pairs[o,]
  o <- mixedorder(pairs$Patient_ID)
  pairs <- pairs[o,]

  pairs
} ## pairs_from_sample()
