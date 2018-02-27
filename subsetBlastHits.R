#!/usr/local/bin/Rscript
args <- commandArgs(TRUE)
dater <- read.delim(args[1], stringsAsFactors=F)
colnames(dater) <- c("qseqid", "qlen", "sseqid", "length", "pident", "evalue", "bitscore")
dater$transcript <- NA

uniqFrame <- as.data.frame(rbind(rep(NA, ncol(dater))), stringsAsFactors=F)
colnames(uniqFrame) <- colnames(dater)
for(i in 1:length(unique(dater$qseqid))) {
  subFrame <- dater[dater$qseqid==unique(dater$qseqid)[i],]
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

write.table(uniqFrame, args[2], quote=FALSE, 
            sep="\t", eol="\n", row.names=F, col.names=T)
