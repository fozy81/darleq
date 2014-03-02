### Based on WFD method statement  
### This function requires input data.frame with columns containing: #Date #Code #sampleID #siteID #Abundance #Taxa #Alk
### This function should work with S plus or R languages

darleqFunc <- function(diatomTDI){  # create function called darleqFunc

  lengthTDI <- length(diatomTDI)
  taxaList <- read.csv("DARLEQ2_TAXON_LIST.csv")  # get DARLEQ taxa list and scores
  diatomTDI <- merge(diatomTDI,taxaList,by.x="Code", by.y="TaxonId") # merge scores with taxa uploaded/input
  
  ## format data & create reference values:
  
  diatomTDI$Alk[diatomTDI$Alk > 250] <- 250  # Alkalinity capped at 250 according to method statement?
  diatomTDI$score <- diatomTDI$TDI4 * diatomTDI$Abundance # create scores for TDI4 'as' value (abundance * TDI3 score)
  diatomTDI$score3 <- diatomTDI$TDI3 * diatomTDI$Abundance # create scores for TDI3 'as' value  (abundance * TDI4 score)
  diatomTDI$scoreL4 <- diatomTDI$LTDI2 * diatomTDI$Abundance # create scores for LTDI2 'as' value
  diatomTDI$scoreL3 <- diatomTDI$LTDI1 * diatomTDI$Abundance # create scores for LTDI2 'as' value + for NEMS TDi3 use $LTDLR for mastertaxonlist use $LTDI1
  diatomTDI$Date <- as.Date(diatomTDI$Date, "%d-%b-%Y")
  diatomTDI$Date <- as.character(diatomTDI$Date)
  diatomTDI$month <- as.numeric(substring(diatomTDI$Date, 6,7)) ### Dares season for TDI3 calculation creates months
  diatomTDI$DaresSeason[diatomTDI$month >= 7] = 1 # create DARES season (is sample in first or second half of year?)
  diatomTDI$DaresSeason[diatomTDI$month <= 6] = 0
  diatomTDI$LochRefValue[diatomTDI$Alk < 10] = 20 ## Create LTDI reference/expected diatom index value based on two Alkalinity bands
  diatomTDI$LochRefValue[diatomTDI$Alk >= 10] = 25
  diatomTDI$eLTDI2 <- ifelse(diatomTDI$Alk >= 10 & diatomTDI$Alk <= 50,35, diatomTDI$Alk)
  diatomTDI$eLTDI2[diatomTDI$Alk < 10] = 22 ## Create LTDI2 reference/expected diatom index value based on three Alkalinity bands
  diatomTDI$eLTDI2[diatomTDI$Alk > 50] = 42
  
  dataTDI <- lapply(split(diatomTDI, diatomTDI$SampleID), function(TDI){ #creates a new function which applies a series of commands on the split data by sample ID  
    
    outDF <- 0 # creates new output data.frame for data to be put into
    
    ## Abundances & Sums (ignores any 'NAs' missing data i.e. na.rm=TRUE):
  
    outDF$sumAbund <- sum(TDI$Abundance[TDI$Planktic != TRUE],na.rm=TRUE) # sum abundance but exclude planktonic
    outDF$tdi3SumAbund <- sum(TDI$Abundance[TDI$TDI3 > 0],na.rm=TRUE) # sum abundance of scoring TDI3 taxa
    outDF$tdi4SumAbund <- sum(TDI$Abundance[TDI$TDI4 > 0],na.rm=TRUE) # sum abundance of scoring TDI4 taxa
    outDF$allTaxaAbund <- nrow(TDI) # count all taxa (the number of rows i.e. should be a row per taxa)
    outDF$allTaxaAbund3 <- sum(TDI$Abundance[TDI$TDI3 >= 0],na.rm=TRUE) # count all taxa abundance checking TDI.x field is not null/na 
    outDF$allTaxaAbund4 <-  sum(TDI$Abundance[TDI$TDI4 >= 0],na.rm=TRUE) # count all taxa abundance checking TDI.x field is not null/na 
    outDF$notPlankticAbund <- sum(TDI$Abundance[TDI$Planktic == FALSE],na.rm=TRUE) # count non-planktonic taxa to fit with NEMS field: 'Abundance of non-planktonic taxa'
    
    #### percentages TDI4:
    
    outDF$plankticAbund <- sum(TDI$Abundance[TDI$Planktic == TRUE],na.rm=TRUE) # count planktonic taxa
    outDF$plankticPercent <- (outDF$plankticAbund / outDF$allTaxaAbund4) * 100 # percentage planktonic using TDI4
    outDF$motileAbund <- sum(TDI$Abundance[TDI$Motile == TRUE],na.rm=TRUE) # sum abundance of motile taxa
    outDF$motilePercent <- (outDF$motileAbund / outDF$allTaxaAbund4) * 100 # percentage motile using TDI4
    outDF$organicAbund <- sum(TDI$Abundance[TDI$OrganicTolerant == TRUE],na.rm=TRUE)  # sum abundance of organic taxa
    outDF$organicPercent <- (outDF$organicAbund / outDF$allTaxaAbund4) * 100 # percentage organic using TDI4
    
    #### TDI3 & TDI4 & eTDI4 & eTDI3 scores:
    
    outDF$as <- sum(TDI$score,na.rm=TRUE)
    outDF$as3 <- sum(TDI$score3,na.rm=TRUE)
    outDF$TDI4 <- ((outDF$as / outDF$tdi4SumAbund)*25) - 25 
    outDF$TDI3 <- ((outDF$as3 / outDF$tdi3SumAbund)*25) - 25
    outDF$eTDI <-  9.933*exp(log10(unique(TDI$Alk))*0.81)
    outDF$eTDI3 <- -25.36+(56.83*(log10(unique(TDI$Alk))))-(12.96*(log10(unique(TDI$Alk))*log10(unique(TDI$Alk))))+(3.21*(unique(TDI$DaresSeason)))
    outDF$EQR <- (100-outDF$TDI4)/(100-outDF$eTDI)
    outDF$tdi3EQR <-  (100-outDF$TDI3)/(100-outDF$eTDI3)
    outDF$SiteID <- unique(TDI$SiteID)
    outDF$Date <- as.character(unique(TDI$Date))
    outDF$SampleID <- unique(TDI$SampleID)
      
    ### Loch LTDI2 scores:
   
    outDF$totalabundsum <- sum(TDI$Abundance,na.rm=TRUE)
    outDF$sumLTDI2 <- sum(TDI$Abundance[TDI$LTDI2 > 0],na.rm=TRUE) #  should this exclude zero scoring taxa? -
    outDF$LTDI4SumAbund <- sum(TDI$scoreL4,na.rm=TRUE) # sum of 
    outDF$w <- outDF$LTDI4SumAbund / outDF$sumLTDI2
    outDF$LTDI2 <- (outDF$w * 25) - 25
    outDF$'EQR LTDI2' <- (100 - (outDF$LTDI2)) / (100 - (unique(TDI$eLTDI2)))
    outDF$eLTDI2 <- unique(TDI$eLTDI2)
    
    ### Loch LTDI1 scores:
    
    outDF$sumLTDI1 <- sum(TDI$Abundance[TDI$LTDI1 > 0],na.rm=TRUE) # should this exclude zero scoring taxa? -
    outDF$LTDI3SumAbund <- sum(TDI$scoreL3,na.rm=TRUE) # issue with duplicate taxa in NEMS
    outDF$w2 <- outDF$LTDI3SumAbund / outDF$sumLTDI1
    outDF$lochTDI3 <- (outDF$w2 * 25) - 25
    outDF$'EQR LTDI1' <- (100 - (outDF$lochTDI3)) / (100 - (unique(TDI$LochRefValue)))
     
    lengthTDI <- length(TDI) # how many columns in input data
    if (lengthTDI > 7){      # does input data include extra column for waterbodyID
      outDF$WaterbodyID <- unique(TDI$WaterbodyID)  # add waterbodyID to output data
    }
   
      
   return(outDF) ##runs output and stores under function name i.e. creates 'splitTDI2' as a list
    
  })
  
 dataTDI <- do.call(rbind, lapply(dataTDI, data.frame, stringsAsFactors=FALSE,check.names=F))
  row.names(dataTDI) <- NULL  # remove row names not required for display
 dataTDI[,1] <- NULL # removes empty column created when outDF was created as start of function (i.e. outDF <- 0)
 
 if (lengthTDI > 7){   # check if waterbodyID in data.frame
wbEQR <- lapply(split(dataTDI, dataTDI$WaterbodyID), function(EQR){ # split by waterbody
 Eqr <- 0
  Eqr$WBEQR <- mean(as.numeric(EQR$'EQR LTDI2')) # create mean waterbody EQR
  Eqr$Waterbody <- unique(EQR$WaterbodyID)
 return(Eqr)  }) 

wbEQR <- do.call(rbind, lapply(wbEQR, data.frame, stringsAsFactors=FALSE,check.names=F)) create data.frame of wb eqr
#dataTDI$wbEQR[] <- wbEQR$WBEQR[as.numeric(dataTDI$WaterbodyID) == as.numeric(wbEQR$Waterbody)] - not working but easier way to merge?
dataTDI <- merge(dataTDI,wbEQR,by.x="WaterbodyID",by.y="Waterbody",all.x=TRUE) # merge wb eqr with existing output data.frame

 }

 return(dataTDI) 
  
}