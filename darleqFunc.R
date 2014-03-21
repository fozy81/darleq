### README
# Based on standalone desktop tool: http://www.staff.ncl.ac.uk/staff/stephen.juggins/DARLEQII.htm 
# + WFD Method statements for Rivers and Lake Phytobenthos
# WFD UKTAG document: http://www.wfduk.org/sites/default/files/Media/Environmental%20standards/Annex%202%20Rivers%20Macrophytes%20%26%20Phytobenthos%20DARLEQ.pdf
# WFD UKTAG document: http://www.wfduk.org/sites/default/files/Media/Environmental%20standards/Annex%2010%20Lakes%20Macrophytes%20and%20Phytobenthos%20DARLEQ.pdf
# This function requires input data.frame with columns containing: #Date #Code #SampleID #siteID #Abundance #Taxa #Alk
# Optional #WaterbodyID column will create mean eqr for waterbodyID
# This function should work with S plus but will require some changes for importing csv and Time/Date function

darleqFunc <- function(diatomTDI){  # create function called darleqFunc

  lengthTDI <- length(diatomTDI) # see if optional waterbody ID included
  taxaList <- read.csv("DARLEQ2_TAXON_LIST.csv")  # get DARLEQ taxa list and scores
  diatomTDI <- merge(diatomTDI,taxaList,by.x="Code", by.y="TaxonId") # merge scores with taxa uploaded/input
  
  ## format data & create reference values:
  
  diatomTDI$Alk[diatomTDI$Alk > 250] <- 250  # Alkalinity capped at 250 according to method statement?
  diatomTDI$score <- diatomTDI$TDI4 * diatomTDI$Abundance # create scores for TDI4 'as' value (abundance * TDI3 score)
  diatomTDI$score3 <- diatomTDI$TDI3 * diatomTDI$Abundance # create scores for TDI3 'as' value  (abundance * TDI4 score)
  diatomTDI$scoreL3 <- diatomTDI$LTDI2 * diatomTDI$Abundance # create scores for LTDI2 'as' value
  diatomTDI$scoreL1 <- diatomTDI$LTDI1 * diatomTDI$Abundance # create scores for LTDI1 'as' value + for NEMS TDi3 use $LTDLR for mastertaxonlist use $LTDI1
  diatomTDI$Date <- as.Date(diatomTDI$Date, "%d-%b-%Y")
  diatomTDI$Date <- as.character(diatomTDI$Date)
  diatomTDI$month <- as.numeric(substring(diatomTDI$Date, 6,7)) ### Dares season for TDI3 calculation creates months
  diatomTDI$DaresSeason[diatomTDI$month >= 7] = 1 # create DARES season (is sample in first or second half of year?)
  diatomTDI$DaresSeason[diatomTDI$month <= 6] = 0
  diatomTDI$eLTDI1[diatomTDI$Alk < 10] = 20 ## Create LTDI reference/expected diatom index value based on two Alkalinity bands
  diatomTDI$eLTDI1[diatomTDI$Alk >= 10] = 25
  diatomTDI$eLTDI2 <- ifelse(diatomTDI$Alk >= 10 & diatomTDI$Alk <= 50,35, diatomTDI$Alk)
  diatomTDI$eLTDI2[diatomTDI$Alk < 10] = 22 ## Create LTDI2 reference/expected diatom index value based on three Alkalinity bands
  diatomTDI$eLTDI2[diatomTDI$Alk > 50] = 42
  
  dataTDI <- lapply(split(diatomTDI, diatomTDI$SampleID), function(TDI){ #creates a new function which applies a series of commands on the split data by sample ID  
    
    outDF <- 0 # creates new output data.frame for data to be put into
    
    ## TDI Abundances & Sums (ignores any 'NAs' missing data i.e. na.rm=TRUE):
     
    outDF$'RIVER TDI3 SumAbund' <- sum(TDI$Abundance[TDI$TDI3 > 0],na.rm=TRUE) # sum abundance of scoring TDI3 taxa
    outDF$'RIVER TDI4 SumAbund' <- sum(TDI$Abundance[TDI$TDI4 > 0],na.rm=TRUE) # sum abundance of scoring TDI4 taxa
    outDF$'LAKE LTDI1 SumAbund' <- sum(TDI$Abundance[TDI$LTDI1 > 0],na.rm=TRUE) # sum abundance of scoring TDI3 taxa
    outDF$'LAKE LTDI2 SumAbund' <- sum(TDI$Abundance[TDI$LTDI2 > 0],na.rm=TRUE) # sum abundance of scoring TDI4 taxa
    outDF$'SAMPLE Total Abundance' <- sum(TDI$Abundance, na.rm=TRUE) # sum abundance - all taxa
    
    #### Plantic/Organic/Motile percentages:
    
  outDF$'SAMPLE plankticAbund' <- sum(TDI$Abundance[TDI$Planktic == TRUE],na.rm=TRUE) # count planktonic taxa
  outDF$'SAMPLE motileAbund' <- sum(TDI$Abundance[TDI$Motile == TRUE],na.rm=TRUE) # sum abundance of motile taxa
  outDF$'SAMPLE organicAbund' <- sum(TDI$Abundance[TDI$OrganicTolerant == TRUE],na.rm=TRUE)  # sum abundance of organic taxa
 
  outDF$'SAMPLE plankticPercent' <- (outDF$'SAMPLE plankticAbund' / outDF$'SAMPLE Total Abundance') * 100 # percentage planktonic 
  outDF$'SAMPLE motilePercent' <- (outDF$'SAMPLE motileAbund' / outDF$'SAMPLE Total Abundance') * 100 # percentage motile 
  outDF$'SAMPLE organicPercent' <- (outDF$'SAMPLE organicAbund' / outDF$'SAMPLE Total Abundance') * 100 # percentage organic 
  
  #### TDI3 & TDI4 & eTDI4 & eTDI3 scores:
    
  outDF$'RIVER as TDI4' <- sum(TDI$score,na.rm=TRUE)
  outDF$'RIVER as3 TDI3' <- sum(TDI$score3,na.rm=TRUE)
  outDF$'RIVER TDI4' <- ((outDF$'RIVER as TDI4' / outDF$'RIVER TDI4 SumAbund')*25) - 25 
  outDF$'RIVER TDI3' <- ((outDF$'RIVER as3 TDI3' / outDF$'RIVER TDI3 SumAbund')*25) - 25
  outDF$'RIVER eTDI4' <-  9.933*exp(log10(unique(TDI$Alk))*0.81)
  outDF$'RIVER eTDI3' <- -25.36+(56.83*(log10(unique(TDI$Alk))))-(12.96*(log10(unique(TDI$Alk))*log10(unique(TDI$Alk))))+(3.21*(unique(TDI$DaresSeason)))
  outDF$'RIVER TDI4 EQR' <- (100-outDF$'RIVER TDI4')/(100-outDF$'RIVER eTDI4')
  outDF$'RIVER TDI3 EQR' <-  (100-outDF$'RIVER TDI3')/(100-outDF$'RIVER eTDI3')
    
  ### Date, Site, Sample
    
  outDF$'SAMPLE SiteID' <- unique(TDI$SiteID)
  outDF$'SAMPLE Date' <- as.character(unique(TDI$Date))
  outDF$'SAMPLE ID' <- unique(TDI$SampleID)
        
  ### LAKE LTDI2 scores:
   
  outDF$'LAKE totalabundsum' <- sum(TDI$Abundance,na.rm=TRUE)
  outDF$'LAKE LTDI2 SumScore' <- sum(TDI$scoreL3,na.rm=TRUE) # sum of 
  outDF$'LAKE w' <- outDF$'LAKE LTDI2 SumScore' / outDF$'LAKE LTDI2 SumAbund'
  outDF$'LAKE LTDI2' <- (outDF$'LAKE w'  * 25) - 25
  outDF$'LAKE EQR LTDI2' <- (100 - (outDF$'LAKE LTDI2')) / (100 - (unique(TDI$eLTDI2)))
  outDF$'LAKE eLTDI2' <- unique(TDI$eLTDI2)
    
  ### LAKE LTDI1 scores:
    
  outDF$'LAKE LTDI1 SumScore' <- sum(TDI$scoreL1,na.rm=TRUE) # issue with duplicate taxa in NEMS
  outDF$'LAKE w2' <- outDF$'LAKE LTDI1 SumScore' / outDF$'LAKE LTDI1 SumAbund'
  outDF$'LAKE LTDI1' <- (outDF$'LAKE w2' * 25) - 25
  outDF$'LAKE EQR LTDI1' <- (100 - (outDF$'LAKE LTDI1')) / (100 - (unique(TDI$eLTDI1)))
 # outDF$'LAKE EQR LTDI1 capped at 1.0' <- 1[outDF$'LAKE EQR LTDI1' >= 1]
 # outDF$'LAKE EQR LTDI1 capped at 1.0' <- outDF$'LAKE EQR LTDI1'[outDF$'LAKE EQR LTDI1' < 1]  
  
 lengthTDI <- length(TDI) # how many columns in input data
   if (lengthTDI > 7){      # does input data include extra column for waterbodyID
   outDF$'SAMPLE WaterbodyID' <- unique(TDI$WaterbodyID)  # add optional waterbodyID to output data
    }
     return(outDF) ##runs output and stores under function name i.e. creates 'splitTDI2' as a list
   })
  
 dataTDI <- do.call(rbind, lapply(dataTDI, data.frame, stringsAsFactors=FALSE,check.names=F))
 
 dataTDI[,1] <- NULL # removes empty column created when outDF was created as start of function (i.e. outDF <- 0)
 
 if (lengthTDI > 7){   # check if waterbodyID in data.frame
wbEQR <- lapply(split(dataTDI, dataTDI$'SAMPLE WaterbodyID'), function(EQR){ # split by waterbody
 Eqr <- 0
 std <- function(x) sd(x)/sqrt(length(x)) # standard error function
 Eqr$'LAKE WB STANDARD ERROR LTDI2' <- std(as.numeric(EQR$'LAKE EQR LTDI2')) # create SE for EQR
 Eqr$'LAKE WB STANDARD ERROR LTDI1' <- std(as.numeric(EQR$'LAKE EQR LTDI1')) 
 Eqr$'RIVER WB STANDARD ERROR TDI3' <- std(as.numeric(EQR$'RIVER TDI3 EQR')) 
 Eqr$'RIVER WB STANDARD ERROR TDI4' <- std(as.numeric(EQR$'RIVER TDI4 EQR')) 
 Eqr$'LAKE WB EQR LTDI2'  <- mean(as.numeric(EQR$'LAKE EQR LTDI2')) # create mean waterbody EQR
 Eqr$'LAKE WB EQR LTDI1'  <- mean(as.numeric(EQR$'LAKE EQR LTDI1')) 
 Eqr$'RIVER WB EQR TDI3' <- mean(as.numeric(EQR$'RIVER TDI3 EQR')) 
 Eqr$'RIVER WB EQR TDI4' <- mean(as.numeric(EQR$'RIVER TDI4 EQR')) 
 Eqr$Waterbody <- unique(EQR$'SAMPLE WaterbodyID')
 Eqr$numberOfSamplesInWaterBody <- length(unique(EQR$'SAMPLE ID'))
 return(Eqr)  }) 

wbEQR <- do.call(rbind, lapply(wbEQR, data.frame, stringsAsFactors=FALSE,check.names=F)) #create data.frame of wb eqr
wbEQR[,1] <- NULL  # removes empty column created when outDF was created as start of function (i.e. Eqr <- 0)
dataTDI <- merge(dataTDI,wbEQR,by.x="SAMPLE WaterbodyID",by.y="Waterbody",all.x=TRUE) # merge wb eqr with existing output data.frame

 }
 return(dataTDI) 
}