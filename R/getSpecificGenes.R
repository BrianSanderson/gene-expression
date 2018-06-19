# Sex-specific Genes Function
# Author: Brian J. Sanderson <brian.sanderson@ttu.edu>
#
# Quantify genes that exhibit sex- or tissue-limited expression, estimated by a total count threshold or 
# a CPM threshold. 
# 
# <dgeObj>: A DGElist object containing read counts for all libraries, with group names corresponding to 
#           male leaves "Ml", female leaves "Fl", male flowers "Mf", and female flowers "Ff"
#
# <incLib>: The number of libraries that are required to pass the threshold for inclusion.
#
# <excLib>: The number of libraries that are required to pass the threshold for exclusion.
# 
# <thresh>: The threshold for including libraries, either in CPM or total reads.
# 
# <type>: Determines whether the threshold is based on total reads (counts), or counts per million (CPM). 
#         CPM is the default.
# 
# Returns a data frame of genes, with an indication of in which tissues the genes are expressed (Y/N).
#
# Note: The letters A, B, C, & D correspond to different sex-tissue combinations. 
# A: male flowers, B: male leaves, c: female flowers, d: female leaves. 
# This was convenient for grouping (e.g. male flowers and male leaves is AB).


getSpecificGenes <- function(dgeObj, incLib=5, excLib=5, thresh=0.1, type="CPM") {

    if(type=="CPM") {
        A <- dgeObj[rowSums(cpm(dgeObj[,dgeObj$samples$group=="Ml"]) < thresh) >= excLib &
                    rowSums(cpm(dgeObj[,dgeObj$samples$group=="Fl"]) < thresh) >= excLib &
                    rowSums(cpm(dgeObj[,dgeObj$samples$group=="Mf"]) >= thresh) >= incLib &
                    rowSums(cpm(dgeObj[,dgeObj$samples$group=="Ff"]) < thresh) >= excLib,
                   ,
                   keep.lib.sizes=FALSE]

        A <- as.data.frame(cbind(as.character(names(rowSums(A$counts)))),
                           stringsAsFactors=F)
        colnames(A) <- c("gene")
        if(nrow(A)>0) {
            A$maleFlowers <- "Y"
            A$femaleFlowers <- "N"
            A$maleLeaves <- "N"
            A$femaleLeaves <- "N"
        } else {
            A$maleFlowers <- character(0)
            A$femaleFlowers <- character(0)
            A$maleLeaves <- character(0)
            A$femaleLeaves <- character(0)
        }
        
  
        B <- dgeObj[rowSums(cpm(dgeObj[,dgeObj$samples$group=="Ml"]) >= thresh) >= incLib &        
                    rowSums(cpm(dgeObj[,dgeObj$samples$group=="Fl"]) < thresh) >= excLib &                 
                    rowSums(cpm(dgeObj[,dgeObj$samples$group=="Mf"]) < thresh) >= excLib &                 
                    rowSums(cpm(dgeObj[,dgeObj$samples$group=="Ff"]) < thresh) >= excLib,               
                    ,               
                    keep.lib.sizes=FALSE]
  
        B <- as.data.frame(cbind(as.character(names(rowSums(B$counts)))),                   
                           stringsAsFactors=F)
  
        colnames(B) <- c("gene")
        if(nrow(B)>0) {
            B$maleFlowers <- "N"
            B$femaleFlowers <- "N"
            B$maleLeaves <- "Y"
            B$femaleLeaves <- "N"
        } else {
            B$maleFlowers <- character(0)
            B$femaleFlowers <- character(0)
            B$maleLeaves <- character(0)
            B$femaleLeaves <- character(0)
        }


  
        C <- dgeObj[rowSums(cpm(dgeObj[,dgeObj$samples$group=="Ml"]) < thresh) >= excLib &
                    rowSums(cpm(dgeObj[,dgeObj$samples$group=="Fl"]) < thresh) >= excLib &                 
                    rowSums(cpm(dgeObj[,dgeObj$samples$group=="Mf"]) < thresh) >= excLib &                 
                    rowSums(cpm(dgeObj[,dgeObj$samples$group=="Ff"]) >= thresh) >= incLib,               
                    ,               
                    keep.lib.sizes=FALSE]
  
        C <- as.data.frame(cbind(as.character(names(rowSums(C$counts)))),                   
                           stringsAsFactors=F)
  
        colnames(C) <- c("gene")
        if(nrow(C)>0) {
            C$maleFlowers <- "N"
            C$femaleFlowers <- "Y"
            C$maleLeaves <- "N"
            C$femaleLeaves <- "N"
        } else {
            C$maleFlowers <- character(0)
            C$femaleFlowers <- character(0)
            C$maleLeaves <- character(0)
            C$femaleLeaves <- character(0)
        }
  
        D <- dgeObj[rowSums(cpm(dgeObj[,dgeObj$samples$group=="Ml"]) < thresh) >= excLib &       
                    rowSums(cpm(dgeObj[,dgeObj$samples$group=="Fl"]) >= thresh) >= incLib &                 
                    rowSums(cpm(dgeObj[,dgeObj$samples$group=="Mf"]) < thresh) >= excLib &                 
                    rowSums(cpm(dgeObj[,dgeObj$samples$group=="Ff"]) < thresh) >= excLib,               
                    ,               
                    keep.lib.sizes=FALSE]
  
        D <- as.data.frame(cbind(as.character(names(rowSums(D$counts)))),                   
                           stringsAsFactors=F)
  
        colnames(D) <- c("gene")
        if(nrow(D)>0) {
            D$maleFlowers <- "N"
            D$femaleFlowers <- "N"
            D$maleLeaves <- "N"
            D$femaleLeaves <- "Y"
        } else {
            D$maleFlowers <- character(0)
            D$femaleFlowers <- character(0)
            D$maleLeaves <- character(0)
            D$femaleLeaves <- character(0)
        }

        AB <- dgeObj[rowSums(cpm(dgeObj[,dgeObj$samples$group=="Ml"]) >= thresh) >= incLib &        
                     rowSums(cpm(dgeObj[,dgeObj$samples$group=="Fl"]) < thresh) >= excLib &                  
                     rowSums(cpm(dgeObj[,dgeObj$samples$group=="Mf"]) >= thresh) >= incLib &                  
                     rowSums(cpm(dgeObj[,dgeObj$samples$group=="Ff"]) < thresh) >= excLib,                
                     ,                
                     keep.lib.sizes=FALSE]
  
        AB <- as.data.frame(cbind(as.character(names(rowSums(AB$counts)))),                   
                            stringsAsFactors=F)

        colnames(AB) <- c("gene")
        if(nrow(AB)>0) {
            AB$maleFlowers <- "Y"
            AB$femaleFlowers <- "N"
            AB$maleLeaves <- "Y"
            AB$femaleLeaves <- "N"
        } else {
            AB$maleFlowers <- character(0)
            AB$femaleFlowers <- character(0)
            AB$maleLeaves <- character(0)
            AB$femaleLeaves <- character(0)
        }
  
        AC <- dgeObj[rowSums(cpm(dgeObj[,dgeObj$samples$group=="Ml"]) < thresh) >= excLib &        
                     rowSums(cpm(dgeObj[,dgeObj$samples$group=="Fl"]) < thresh) >= excLib &                  
                     rowSums(cpm(dgeObj[,dgeObj$samples$group=="Mf"]) >= thresh) >= incLib &                  
                     rowSums(cpm(dgeObj[,dgeObj$samples$group=="Ff"]) >= thresh) >= incLib,                
                     ,      
                     keep.lib.sizes=FALSE]

 
        AC <- as.data.frame(cbind(as.character(names(rowSums(AC$counts)))), 
                            stringsAsFactors=F)

        colnames(AC) <- c("gene")
        if(nrow(AC)>0) {
            AC$maleFlowers <- "Y"
            AC$femaleFlowers <- "Y"
            AC$maleLeaves <- "N"
            AC$femaleLeaves <- "N"
        } else {
            AC$maleFlowers <- character(0)
            AC$femaleFlowers <- character(0)
            AC$maleLeaves <- character(0)
            AC$femaleLeaves <- character(0)
        }

        AD <- dgeObj[rowSums(cpm(dgeObj[,dgeObj$samples$group=="Ml"]) < thresh) >= excLib &  
                     rowSums(cpm(dgeObj[,dgeObj$samples$group=="Fl"]) >= thresh) >= incLib &
                     rowSums(cpm(dgeObj[,dgeObj$samples$group=="Mf"]) >= thresh) >= incLib &   
                     rowSums(cpm(dgeObj[,dgeObj$samples$group=="Ff"]) < thresh) >= excLib,
                     ,
                     keep.lib.sizes=FALSE]

        AD <- as.data.frame(cbind(as.character(names(rowSums(AD$counts)))),          
                            stringsAsFactors=F)

  
        colnames(AD) <- c("gene")
        if(nrow(AD)>0) {
            AD$maleFlowers <- "Y"
            AD$femaleFlowers <- "N"
            AD$maleLeaves <- "N"
            AD$femaleLeaves <- "Y"
        } else {
            AD$maleFlowers <- character(0)
            AD$femaleFlowers <- character(0)
            AD$maleLeaves <- character(0)
            AD$femaleLeaves <- character(0)
        }
  
        BC <- dgeObj[rowSums(cpm(dgeObj[,dgeObj$samples$group=="Ml"]) >= thresh) >= incLib &        
                     rowSums(cpm(dgeObj[,dgeObj$samples$group=="Fl"]) < thresh) >= excLib &                  
                     rowSums(cpm(dgeObj[,dgeObj$samples$group=="Mf"]) < thresh) >= excLib &                 
                     rowSums(cpm(dgeObj[,dgeObj$samples$group=="Ff"]) >= thresh) >= incLib,                
                     ,                
                     keep.lib.sizes=FALSE]
  
        BC <- as.data.frame(cbind(as.character(names(rowSums(BC$counts)))),                    
                            stringsAsFactors=F)

        colnames(BC) <- c("gene")
        if(nrow(BC)>0) {
            BC$maleFlowers <- "N"
            BC$femaleFlowers <- "Y"
            BC$maleLeaves <- "Y"
            BC$femaleLeaves <- "N"
        } else {
            BC$maleFlowers <- character(0)
            BC$femaleFlowers <- character(0)
            BC$maleLeaves <- character(0)
            BC$femaleLeaves <- character(0)
        }

  
        BD <- dgeObj[rowSums(cpm(dgeObj[,dgeObj$samples$group=="Ml"]) >= thresh) >= incLib &              
                     rowSums(cpm(dgeObj[,dgeObj$samples$group=="Fl"]) >= thresh) >= incLib &                            
                     rowSums(cpm(dgeObj[,dgeObj$samples$group=="Mf"]) < thresh) >= excLib &               
                     rowSums(cpm(dgeObj[,dgeObj$samples$group=="Ff"]) < thresh) >= excLib,                     
                     ,
                     keep.lib.sizes=FALSE]
 
        BD <- as.data.frame(cbind(as.character(names(rowSums(BD$counts)))),                   
                            stringsAsFactors=F)


        colnames(BD) <- c("gene")
        if(nrow(BD)>0) {
            BD$maleFlowers <- "N"
            BD$femaleFlowers <- "N"
            BD$maleLeaves <- "Y"
            BD$femaleLeaves <- "Y"
        } else {
            BD$maleFlowers <- character(0)
            BD$femaleFlowers <- character(0)
            BD$maleLeaves <- character(0)
            BD$femaleLeaves <- character(0)
        }

  
        CD <- dgeObj[rowSums(cpm(dgeObj[,dgeObj$samples$group=="Ml"]) < thresh) >= excLib &
                     rowSums(cpm(dgeObj[,dgeObj$samples$group=="Fl"]) >= thresh) >= incLib &                  
                     rowSums(cpm(dgeObj[,dgeObj$samples$group=="Mf"]) < thresh) >= excLib &                  
                     rowSums(cpm(dgeObj[,dgeObj$samples$group=="Ff"]) >= thresh) >= incLib,                
                     ,                
                     keep.lib.sizes=FALSE]
  
        CD <- as.data.frame(cbind(as.character(names(rowSums(CD$counts)))),                    
                            stringsAsFactors=F)
  
        colnames(CD) <- c("gene")
        if(nrow(CD)>0) {
            CD$maleFlowers <- "N"
            CD$femaleFlowers <- "Y"
            CD$maleLeaves <- "N"
            CD$femaleLeaves <- "Y"
        } else {
            CD$maleFlowers <- character(0)
            CD$femaleFlowers <- character(0)
            CD$maleLeaves <- character(0)
            CD$femaleLeaves <- character(0)
        }
 
        ABC <- dgeObj[rowSums(cpm(dgeObj[,dgeObj$samples$group=="Ml"]) >= thresh) >= incLib &        
                      rowSums(cpm(dgeObj[,dgeObj$samples$group=="Fl"]) < thresh) >= excLib &                  
                      rowSums(cpm(dgeObj[,dgeObj$samples$group=="Mf"]) >= thresh) >= incLib &                  
                      rowSums(cpm(dgeObj[,dgeObj$samples$group=="Ff"]) >= thresh) >= incLib,                
                      ,                 
                      keep.lib.sizes=FALSE]
  
        ABC <- as.data.frame(cbind(as.character(names(rowSums(ABC$counts)))),                    
                             stringsAsFactors=F)
  
        colnames(ABC) <- c("gene")
        if(nrow(ABC)>0) {
            ABC$maleFlowers <- "Y"
            ABC$femaleFlowers <- "Y"
            ABC$maleLeaves <- "Y"
            ABC$femaleLeaves <- "N"
        } else {
            ABC$maleFlowers <- character(0)
            ABC$femaleFlowers <- character(0)
            ABC$maleLeaves <- character(0)
            ABC$femaleLeaves <- character(0)
        }
  
        ABD <- dgeObj[rowSums(cpm(dgeObj[,dgeObj$samples$group=="Ml"]) >= thresh) >= incLib &        
                      rowSums(cpm(dgeObj[,dgeObj$samples$group=="Fl"]) >= thresh) >= incLib &                   
                      rowSums(cpm(dgeObj[,dgeObj$samples$group=="Mf"]) >= thresh) >= incLib &                   
                      rowSums(cpm(dgeObj[,dgeObj$samples$group=="Ff"]) < thresh) >= excLib,                 
                      ,                 
                      keep.lib.sizes=FALSE]
  
        ABD <- as.data.frame(cbind(as.character(names(rowSums(ABD$counts)))),                     
                             stringsAsFactors=F)
  
        colnames(ABD) <- c("gene")
        if(nrow(ABD)>0) {
            ABD$maleFlowers <- "Y"
            ABD$femaleFlowers <- "N"
            ABD$maleLeaves <- "Y"
            ABD$femaleLeaves <- "Y"
        } else {
            ABD$maleFlowers <- character(0)
            ABD$femaleFlowers <- character(0)
            ABD$maleLeaves <- character(0)
            ABD$femaleLeaves <- character(0)
        }
  
        ACD <- dgeObj[rowSums(cpm(dgeObj[,dgeObj$samples$group=="Ml"]) < thresh) >= excLib &        
                      rowSums(cpm(dgeObj[,dgeObj$samples$group=="Fl"]) >= thresh) >= incLib &                   
                      rowSums(cpm(dgeObj[,dgeObj$samples$group=="Mf"]) >= thresh) >= incLib &                   
                      rowSums(cpm(dgeObj[,dgeObj$samples$group=="Ff"]) >= thresh) >= incLib,                 
                      ,                 
                      keep.lib.sizes=FALSE]
  
        ACD <- as.data.frame(cbind(as.character(names(rowSums(ACD$counts)))),       
                             stringsAsFactors=F)
  
        colnames(ACD) <- c("gene")
        if(nrow(ACD)>0) {
            ACD$maleFlowers <- "Y"
            ACD$femaleFlowers <- "Y"
            ACD$maleLeaves <- "N"
            ACD$femaleLeaves <- "Y"
        } else {
            ACD$maleFlowers <- character(0)
            ACD$femaleFlowers <- character(0)
            ACD$maleLeaves <- character(0)
            ACD$femaleLeaves <- character(0)
        }
  
        BCD <- dgeObj[rowSums(cpm(dgeObj[,dgeObj$samples$group=="Ml"]) >= thresh) >= incLib &        
                      rowSums(cpm(dgeObj[,dgeObj$samples$group=="Fl"]) >= thresh) >= incLib &                  
                      rowSums(cpm(dgeObj[,dgeObj$samples$group=="Mf"]) < thresh) >= excLib &                  
                      rowSums(cpm(dgeObj[,dgeObj$samples$group=="Ff"]) >= thresh) >= incLib,                 
                      ,                 
                      keep.lib.sizes=FALSE]
  
        BCD <- as.data.frame(cbind(as.character(names(rowSums(BCD$counts)))),                     
                             stringsAsFactors=F)
  
        colnames(BCD) <- c("gene")
        if(nrow(BCD)>0) {
            BCD$maleFlowers <- "N"
            BCD$femaleFlowers <- "Y"
            BCD$maleLeaves <- "Y"
            BCD$femaleLeaves <- "Y"
        } else {
            BCD$maleFlowers <- character(0)
            BCD$femaleFlowers <- character(0)
            BCD$maleLeaves <- character(0)
            BCD$femaleLeaves <- character(0)
        }
  
        ABCD <- dgeObj[rowSums(cpm(dgeObj[,dgeObj$samples$group=="Ml"]) >= thresh) >= incLib &       
                       rowSums(cpm(dgeObj[,dgeObj$samples$group=="Fl"]) >= thresh) >= incLib &                  
                       rowSums(cpm(dgeObj[,dgeObj$samples$group=="Mf"]) >= thresh) >= incLib &                   
                       rowSums(cpm(dgeObj[,dgeObj$samples$group=="Ff"]) >= thresh) >= incLib,                  
                       ,                  
                       keep.lib.sizes=FALSE]
 
        ABCD <- as.data.frame(cbind(as.character(names(rowSums(ABCD$counts)))),                     
                              stringsAsFactors=F)
 
        colnames(ABCD) <- c("gene")
        if(nrow(ABCD)>0) {
            ABCD$maleFlowers <- "Y"
            ABCD$femaleFlowers <- "Y"
            ABCD$maleLeaves <- "Y"
            ABCD$femaleLeaves <- "Y"
        } else {
            ABCD$maleFlowers <- character(0)
            ABCD$femaleFlowers <- character(0)
            ABCD$maleLeaves <- character(0)
            ABCD$femaleLeaves <- character(0)
        }
  

        allSpecificGenes <- as.data.frame(rbind(A, B, C, D,
                                                AB, AC, AD, BC, BD, CD,
                                                ABC, ABD, ACD, BCD,
                                                ABCD), stringsAsFactors=F)
        
        return(allSpecificGenes)  
   
    }
    
    if(type=="counts") {
        A <- dgeObj[rowSums(dgeObj[,dgeObj$samples$group=="Ml"]$counts < thresh) >= excLib &
                    rowSums(dgeObj[,dgeObj$samples$group=="Fl"]$counts < thresh) >= excLib &
                    rowSums(dgeObj[,dgeObj$samples$group=="Mf"]$counts >= thresh) >= incLib &
                    rowSums(dgeObj[,dgeObj$samples$group=="Ff"]$counts < thresh) >= excLib,
                   ,
                   keep.lib.sizes=FALSE]

        A <- as.data.frame(cbind(as.character(names(rowSums(A$counts)))),
                           stringsAsFactors=F)
        colnames(A) <- c("gene")
        if(nrow(A)>0) {
            A$maleFlowers <- "Y"
            A$femaleFlowers <- "N"
            A$maleLeaves <- "N"
            A$femaleLeaves <- "N"
        } else {
            A$maleFlowers <- character(0)
            A$femaleFlowers <- character(0)
            A$maleLeaves <- character(0)
            A$femaleLeaves <- character(0)
        }
        
  
        B <- dgeObj[rowSums(dgeObj[,dgeObj$samples$group=="Ml"]$counts >= thresh) >= incLib &        
                    rowSums(dgeObj[,dgeObj$samples$group=="Fl"]$counts < thresh) >= excLib &                 
                    rowSums(dgeObj[,dgeObj$samples$group=="Mf"]$counts < thresh) >= excLib &                 
                    rowSums(dgeObj[,dgeObj$samples$group=="Ff"]$counts < thresh) >= excLib,               
                    ,               
                    keep.lib.sizes=FALSE]
  
        B <- as.data.frame(cbind(as.character(names(rowSums(B$counts)))),                   
                           stringsAsFactors=F)
  
        colnames(B) <- c("gene")
        if(nrow(B)>0) {
            B$maleFlowers <- "N"
            B$femaleFlowers <- "N"
            B$maleLeaves <- "Y"
            B$femaleLeaves <- "N"
        } else {
            B$maleFlowers <- character(0)
            B$femaleFlowers <- character(0)
            B$maleLeaves <- character(0)
            B$femaleLeaves <- character(0)
        }


  
        C <- dgeObj[rowSums(dgeObj[,dgeObj$samples$group=="Ml"]$counts < thresh) >= excLib &
                    rowSums(dgeObj[,dgeObj$samples$group=="Fl"]$counts < thresh) >= excLib &                 
                    rowSums(dgeObj[,dgeObj$samples$group=="Mf"]$counts < thresh) >= excLib &                 
                    rowSums(dgeObj[,dgeObj$samples$group=="Ff"]$counts >= thresh) >= incLib,               
                    ,               
                    keep.lib.sizes=FALSE]
  
        C <- as.data.frame(cbind(as.character(names(rowSums(C$counts)))),                   
                           stringsAsFactors=F)
  
        colnames(C) <- c("gene")
        if(nrow(C)>0) {
            C$maleFlowers <- "N"
            C$femaleFlowers <- "Y"
            C$maleLeaves <- "N"
            C$femaleLeaves <- "N"
        } else {
            C$maleFlowers <- character(0)
            C$femaleFlowers <- character(0)
            C$maleLeaves <- character(0)
            C$femaleLeaves <- character(0)
        }
  
        D <- dgeObj[rowSums(dgeObj[,dgeObj$samples$group=="Ml"]$counts < thresh) >= excLib &       
                    rowSums(dgeObj[,dgeObj$samples$group=="Fl"]$counts >= thresh) >= incLib &                 
                    rowSums(dgeObj[,dgeObj$samples$group=="Mf"]$counts < thresh) >= excLib &                 
                    rowSums(dgeObj[,dgeObj$samples$group=="Ff"]$counts < thresh) >= excLib,               
                    ,               
                    keep.lib.sizes=FALSE]
  
        D <- as.data.frame(cbind(as.character(names(rowSums(D$counts)))),                   
                           stringsAsFactors=F)
  
        colnames(D) <- c("gene")
        if(nrow(D)>0) {
            D$maleFlowers <- "N"
            D$femaleFlowers <- "N"
            D$maleLeaves <- "N"
            D$femaleLeaves <- "Y"
        } else {
            D$maleFlowers <- character(0)
            D$femaleFlowers <- character(0)
            D$maleLeaves <- character(0)
            D$femaleLeaves <- character(0)
        }

        AB <- dgeObj[rowSums(dgeObj[,dgeObj$samples$group=="Ml"]$counts >= thresh) >= incLib &        
                     rowSums(dgeObj[,dgeObj$samples$group=="Fl"]$counts < thresh) >= excLib &                  
                     rowSums(dgeObj[,dgeObj$samples$group=="Mf"]$counts >= thresh) >= incLib &                  
                     rowSums(dgeObj[,dgeObj$samples$group=="Ff"]$counts < thresh) >= excLib,                
                     ,                
                     keep.lib.sizes=FALSE]
  
        AB <- as.data.frame(cbind(as.character(names(rowSums(AB$counts)))),                   
                            stringsAsFactors=F)

        colnames(AB) <- c("gene")
        if(nrow(AB)>0) {
            AB$maleFlowers <- "Y"
            AB$femaleFlowers <- "N"
            AB$maleLeaves <- "Y"
            AB$femaleLeaves <- "N"
        } else {
            AB$maleFlowers <- character(0)
            AB$femaleFlowers <- character(0)
            AB$maleLeaves <- character(0)
            AB$femaleLeaves <- character(0)
        }
  
        AC <- dgeObj[rowSums(dgeObj[,dgeObj$samples$group=="Ml"]$counts < thresh) >= excLib &        
                     rowSums(dgeObj[,dgeObj$samples$group=="Fl"]$counts < thresh) >= excLib &                  
                     rowSums(dgeObj[,dgeObj$samples$group=="Mf"]$counts >= thresh) >= incLib &                  
                     rowSums(dgeObj[,dgeObj$samples$group=="Ff"]$counts >= thresh) >= incLib,                
                     ,      
                     keep.lib.sizes=FALSE]

 
        AC <- as.data.frame(cbind(as.character(names(rowSums(AC$counts)))), 
                            stringsAsFactors=F)

        colnames(AC) <- c("gene")
        if(nrow(AC)>0) {
            AC$maleFlowers <- "Y"
            AC$femaleFlowers <- "Y"
            AC$maleLeaves <- "N"
            AC$femaleLeaves <- "N"
        } else {
            AC$maleFlowers <- character(0)
            AC$femaleFlowers <- character(0)
            AC$maleLeaves <- character(0)
            AC$femaleLeaves <- character(0)
        }

        AD <- dgeObj[rowSums(dgeObj[,dgeObj$samples$group=="Ml"]$counts < thresh) >= excLib &  
                     rowSums(dgeObj[,dgeObj$samples$group=="Fl"]$counts >= thresh) >= incLib &
                     rowSums(dgeObj[,dgeObj$samples$group=="Mf"]$counts >= thresh) >= incLib &   
                     rowSums(dgeObj[,dgeObj$samples$group=="Ff"]$counts < thresh) >= excLib,
                     ,
                     keep.lib.sizes=FALSE]

        AD <- as.data.frame(cbind(as.character(names(rowSums(AD$counts)))),          
                            stringsAsFactors=F)

  
        colnames(AD) <- c("gene")
        if(nrow(AD)>0) {
            AD$maleFlowers <- "Y"
            AD$femaleFlowers <- "N"
            AD$maleLeaves <- "N"
            AD$femaleLeaves <- "Y"
        } else {
            AD$maleFlowers <- character(0)
            AD$femaleFlowers <- character(0)
            AD$maleLeaves <- character(0)
            AD$femaleLeaves <- character(0)
        }
  
        BC <- dgeObj[rowSums(dgeObj[,dgeObj$samples$group=="Ml"]$counts >= thresh) >= incLib &        
                     rowSums(dgeObj[,dgeObj$samples$group=="Fl"]$counts < thresh) >= excLib &                  
                     rowSums(dgeObj[,dgeObj$samples$group=="Mf"]$counts < thresh) >= excLib &                 
                     rowSums(dgeObj[,dgeObj$samples$group=="Ff"]$counts >= thresh) >= incLib,                
                     ,                
                     keep.lib.sizes=FALSE]
  
        BC <- as.data.frame(cbind(as.character(names(rowSums(BC$counts)))),                    
                            stringsAsFactors=F)

        colnames(BC) <- c("gene")
        if(nrow(BC)>0) {
            BC$maleFlowers <- "N"
            BC$femaleFlowers <- "Y"
            BC$maleLeaves <- "Y"
            BC$femaleLeaves <- "N"
        } else {
            BC$maleFlowers <- character(0)
            BC$femaleFlowers <- character(0)
            BC$maleLeaves <- character(0)
            BC$femaleLeaves <- character(0)
        }

  
        BD <- dgeObj[rowSums(dgeObj[,dgeObj$samples$group=="Ml"]$counts >= thresh) >= incLib &              
                     rowSums(dgeObj[,dgeObj$samples$group=="Fl"]$counts >= thresh) >= incLib &                            
                     rowSums(dgeObj[,dgeObj$samples$group=="Mf"]$counts < thresh) >= excLib &               
                     rowSums(dgeObj[,dgeObj$samples$group=="Ff"]$counts < thresh) >= excLib,                     
                     ,
                     keep.lib.sizes=FALSE]
 
        BD <- as.data.frame(cbind(as.character(names(rowSums(BD$counts)))),                   
                            stringsAsFactors=F)


        colnames(BD) <- c("gene")
        if(nrow(BD)>0) {
            BD$maleFlowers <- "N"
            BD$femaleFlowers <- "N"
            BD$maleLeaves <- "Y"
            BD$femaleLeaves <- "Y"
        } else {
            BD$maleFlowers <- character(0)
            BD$femaleFlowers <- character(0)
            BD$maleLeaves <- character(0)
            BD$femaleLeaves <- character(0)
        }

  
        CD <- dgeObj[rowSums(dgeObj[,dgeObj$samples$group=="Ml"]$counts < thresh) >= excLib &
                     rowSums(dgeObj[,dgeObj$samples$group=="Fl"]$counts >= thresh) >= incLib &                  
                     rowSums(dgeObj[,dgeObj$samples$group=="Mf"]$counts < thresh) >= excLib &                  
                     rowSums(dgeObj[,dgeObj$samples$group=="Ff"]$counts >= thresh) >= incLib,                
                     ,                
                     keep.lib.sizes=FALSE]
  
        CD <- as.data.frame(cbind(as.character(names(rowSums(CD$counts)))),                    
                            stringsAsFactors=F)
  
        colnames(CD) <- c("gene")
        if(nrow(CD)>0) {
            CD$maleFlowers <- "N"
            CD$femaleFlowers <- "Y"
            CD$maleLeaves <- "N"
            CD$femaleLeaves <- "Y"
        } else {
            CD$maleFlowers <- character(0)
            CD$femaleFlowers <- character(0)
            CD$maleLeaves <- character(0)
            CD$femaleLeaves <- character(0)
        }
 
        ABC <- dgeObj[rowSums(dgeObj[,dgeObj$samples$group=="Ml"]$counts >= thresh) >= incLib &        
                      rowSums(dgeObj[,dgeObj$samples$group=="Fl"]$counts < thresh) >= excLib &                  
                      rowSums(dgeObj[,dgeObj$samples$group=="Mf"]$counts >= thresh) >= incLib &                  
                      rowSums(dgeObj[,dgeObj$samples$group=="Ff"]$counts >= thresh) >= incLib,                
                      ,                 
                      keep.lib.sizes=FALSE]
  
        ABC <- as.data.frame(cbind(as.character(names(rowSums(ABC$counts)))),                    
                             stringsAsFactors=F)
  
        colnames(ABC) <- c("gene")
        if(nrow(ABC)>0) {
            ABC$maleFlowers <- "Y"
            ABC$femaleFlowers <- "Y"
            ABC$maleLeaves <- "Y"
            ABC$femaleLeaves <- "N"
        } else {
            ABC$maleFlowers <- character(0)
            ABC$femaleFlowers <- character(0)
            ABC$maleLeaves <- character(0)
            ABC$femaleLeaves <- character(0)
        }
  
        ABD <- dgeObj[rowSums(dgeObj[,dgeObj$samples$group=="Ml"]$counts >= thresh) >= incLib &        
                      rowSums(dgeObj[,dgeObj$samples$group=="Fl"]$counts >= thresh) >= incLib &                   
                      rowSums(dgeObj[,dgeObj$samples$group=="Mf"]$counts >= thresh) >= incLib &                   
                      rowSums(dgeObj[,dgeObj$samples$group=="Ff"]$counts < thresh) >= excLib,                 
                      ,                 
                      keep.lib.sizes=FALSE]
  
        ABD <- as.data.frame(cbind(as.character(names(rowSums(ABD$counts)))),                     
                             stringsAsFactors=F)
  
        colnames(ABD) <- c("gene")
        if(nrow(ABD)>0) {
            ABD$maleFlowers <- "Y"
            ABD$femaleFlowers <- "N"
            ABD$maleLeaves <- "Y"
            ABD$femaleLeaves <- "Y"
        } else {
            ABD$maleFlowers <- character(0)
            ABD$femaleFlowers <- character(0)
            ABD$maleLeaves <- character(0)
            ABD$femaleLeaves <- character(0)
        }
  
        ACD <- dgeObj[rowSums(dgeObj[,dgeObj$samples$group=="Ml"]$counts < thresh) >= excLib &        
                      rowSums(dgeObj[,dgeObj$samples$group=="Fl"]$counts >= thresh) >= incLib &                   
                      rowSums(dgeObj[,dgeObj$samples$group=="Mf"]$counts >= thresh) >= incLib &                   
                      rowSums(dgeObj[,dgeObj$samples$group=="Ff"]$counts >= thresh) >= incLib,                 
                      ,                 
                      keep.lib.sizes=FALSE]
  
        ACD <- as.data.frame(cbind(as.character(names(rowSums(ACD$counts)))),       
                             stringsAsFactors=F)
  
        colnames(ACD) <- c("gene")
        if(nrow(ACD)>0) {
            ACD$maleFlowers <- "Y"
            ACD$femaleFlowers <- "Y"
            ACD$maleLeaves <- "N"
            ACD$femaleLeaves <- "Y"
        } else {
            ACD$maleFlowers <- character(0)
            ACD$femaleFlowers <- character(0)
            ACD$maleLeaves <- character(0)
            ACD$femaleLeaves <- character(0)
        }
  
        BCD <- dgeObj[rowSums(dgeObj[,dgeObj$samples$group=="Ml"]$counts >= thresh) >= incLib &        
                      rowSums(dgeObj[,dgeObj$samples$group=="Fl"]$counts >= thresh) >= incLib &                  
                      rowSums(dgeObj[,dgeObj$samples$group=="Mf"]$counts < thresh) >= excLib &                  
                      rowSums(dgeObj[,dgeObj$samples$group=="Ff"]$counts >= thresh) >= incLib,                 
                      ,                 
                      keep.lib.sizes=FALSE]
  
        BCD <- as.data.frame(cbind(as.character(names(rowSums(BCD$counts)))),                     
                             stringsAsFactors=F)
  
        colnames(BCD) <- c("gene")
        if(nrow(BCD)>0) {
            BCD$maleFlowers <- "N"
            BCD$femaleFlowers <- "Y"
            BCD$maleLeaves <- "Y"
            BCD$femaleLeaves <- "Y"
        } else {
            BCD$maleFlowers <- character(0)
            BCD$femaleFlowers <- character(0)
            BCD$maleLeaves <- character(0)
            BCD$femaleLeaves <- character(0)
        }
  
        ABCD <- dgeObj[rowSums(dgeObj[,dgeObj$samples$group=="Ml"]$counts >= thresh) >= incLib &       
                       rowSums(dgeObj[,dgeObj$samples$group=="Fl"]$counts >= thresh) >= incLib &                  
                       rowSums(dgeObj[,dgeObj$samples$group=="Mf"]$counts >= thresh) >= incLib &                   
                       rowSums(dgeObj[,dgeObj$samples$group=="Ff"]$counts >= thresh) >= incLib,                  
                       ,                  
                       keep.lib.sizes=FALSE]
 
        ABCD <- as.data.frame(cbind(as.character(names(rowSums(ABCD$counts)))),                     
                              stringsAsFactors=F)
 
        colnames(ABCD) <- c("gene")
        if(nrow(ABCD)>0) {
            ABCD$maleFlowers <- "Y"
            ABCD$femaleFlowers <- "Y"
            ABCD$maleLeaves <- "Y"
            ABCD$femaleLeaves <- "Y"
        } else {
            ABCD$maleFlowers <- character(0)
            ABCD$femaleFlowers <- character(0)
            ABCD$maleLeaves <- character(0)
            ABCD$femaleLeaves <- character(0)
        }
  

        allSpecificGenes <- as.data.frame(rbind(A, B, C, D,
                                                AB, AC, AD, BC, BD, CD,
                                                ABC, ABD, ACD, BCD,
                                                ABCD), stringsAsFactors=F)
   
    }
}