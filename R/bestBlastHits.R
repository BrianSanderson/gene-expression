# Identify the best homologs based on BLAST results optimizing for alignment quality and length
#
# Author: Brian J. Sanderson <brian.sanderson@ttu.edu>
#
# <blastResults>: A text file of BLAST results containing the following fields: 
#                 "qseqid", "qlen", "sseqid", "length", "pident", "evalue", "bitscore"
#

bestBlastHits <- function(fullResults) {
    colnames(fullResults) <- c("qseqid", "qlen", "sseqid", "length", "pident", "evalue", "bitscore")
    fullResults$transcript <- NA

    uniqFrame <- as.data.frame(rbind(rep(NA, ncol(fullResults))), stringsAsFactors=F)
    colnames(uniqFrame) <- colnames(fullResults)
    for(i in 1:length(unique(fullResults$qseqid))) {
      subFrame <- fullResults[fullResults$qseqid==unique(fullResults$qseqid)[i],]
      cand <- subFrame[subFrame$bitscore==max(subFrame$bitscore),]
      if(nrow(cand) > 1) {
        cand2 <- cand[cand$transcript==min(cand$transcript),]
        if(nrow(cand2) > 1) {
          cand3 <- cand2[cand2$length==max(cand2$length),]
          if(nrow(cand3) > 1) {
            cand3$rand <- rnorm(nrow(cand3))
            cand4 <- cand3[cand3$rand==max(cand3$rand),]
            if(nrow(cand4) > 1) {
              print("WTF")
            } else{uniqFrame <- rbind(uniqFrame, cand4[,1:length(colnames(cand4))-1])}
          } else{uniqFrame <- rbind(uniqFrame, cand3[,1:length(colnames(cand3))-1])}
        } else{uniqFrame <- rbind(uniqFrame, cand2)}
      } else{uniqFrame <- rbind(uniqFrame, cand)}
    }

    uniqFrame <- uniqFrame[!is.na(uniqFrame$qseqid),]

    return(uniqFrame)
}
