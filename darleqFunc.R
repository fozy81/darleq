### function requires #Data #Code #sampleID #siteID #Abundance #Taxa

darleqFunc <- function(diatomTDI2){
  
  require(lubridate)
  
  taxaList <- read.csv("DARLEQ2_TAXON_LIST.csv")  # get taxa list and scores
  diatomTDI2 <- merge(diatomTDI2,taxaList,by.x="Code", by.y="TaxonId") # merge scores with taxa uploaded/input
  
  ## tidy up data:
  
  diatomTDI2$Alk <- as.numeric(diatomTDI2$Alk)
  diatomTDI2$Alk[diatomTDI2$Alk > 250] <- 250
  diatomTDI2$TDI4 <-  as.numeric(diatomTDI2$TDI4)
  diatomTDI2$TDI3 <-  as.numeric(diatomTDI2$TDI3)
  diatomTDI2$Abundance <- as.numeric(diatomTDI2$Abundance)
  diatomTDI2$score <- diatomTDI2$TDI4 * diatomTDI2$Abundance # create scores for TDI4 'as' value
  diatomTDI2$score3 <- diatomTDI2$TDI3 * diatomTDI2$Abundance # create scores for TDI3 'as' value
  diatomTDI2$scoreL4 <- diatomTDI2$LTDI2 * diatomTDI2$Abundance # create scores for LTDI2 'as' value
  diatomTDI2$scoreL3 <- diatomTDI2$LTDI1 * diatomTDI2$Abundance # create scores for LTDI2 'as' value + for NEMS TDi3 use $LTDLR for mastertaxonlist use $LTDI1
  diatomTDI2$Date <- as.Date(diatomTDI2$Date, "%d-%b-%Y")
  diatomTDI2$month <- as.numeric(month(diatomTDI2$Date)) ### Dares season for TDI3 calculation creates months
  diatomTDI2$DaresSeason[diatomTDI2$month >= 7] = 1 # create DARES season (is sample in first or second half of year?)
  diatomTDI2$DaresSeason[diatomTDI2$month <= 6] = 0
  diatomTDI2$LochRefValue[diatomTDI2$Alk < 10] = 20 ## Create LTDI reference/expected diatom index value based on two Alkalinity bands
  diatomTDI2$LochRefValue[diatomTDI2$Alk >= 10] = 25
  diatomTDI2$LochRefValue2 <- ifelse(diatomTDI2$Alk >= 10 & diatomTDI2$Alk <= 50,35, diatomTDI2$Alk)
  diatomTDI2$LochRefValue2[diatomTDI2$Alk < 10] = 22 ## Create LTDI2 reference/expected diatom index value based on three Alkalinity bands
  diatomTDI2$LochRefValue2[diatomTDI2$Alk > 50] = 42
  
 dataf <- lapply(split(diatomTDI2, diatomTDI2$SampleID), function(TDI){ #creates a new function which applies a series of commands on the split data by sample ID  

   ## Abundances & Sums:
  
  sumAbund <- sum(TDI$Abundance[TDI$Planktic != TRUE],na.rm=TRUE) # sum abundance but exclude planktonic
  tdi3SumAbund <- sum(TDI$Abundance[TDI$TDI3 > 0],na.rm=TRUE)
  tdi4SumAbund <- sum(TDI$Abundance[TDI$TDI4 > 0],na.rm=TRUE)# sum abundance of scoring TDI4 taxa
  allTaxaAbund <- nrow(TDI) # count all taxa (the number of rows i.e. should be a row per taxa)
  allTaxaAbund3 <- sum(TDI$Abundance[TDI$TDI3 >= 0],na.rm=TRUE) # count all taxa abundance checking TDI.x field is not null/na 
  allTaxaAbund4 <-  sum(TDI$Abundance[TDI$TDI4 >= 0],na.rm=TRUE) # count all taxa abundance checking TDI.x field is not null/na 
  notPlankticAbund <- sum(TDI$Abundance[TDI$Planktic == FALSE],na.rm=TRUE) # count non-planktonic taxa to fit with NEMS field: 'Abundance of non-planktonic taxa'
  
  #### percentages TDI4:
  
  plankticAbund <- sum(TDI$Abundance[TDI$Planktic == TRUE],na.rm=TRUE) # count planktonic taxa
  plankticPercent <- (plankticAbund / allTaxaAbund4) * 100 # percentage planktonic using TDI4
  motileAbund <- sum(TDI$Abundance[TDI$Motile == TRUE],na.rm=TRUE) # sum abundance of motile taxa
  motilePercent <- (motileAbund / allTaxaAbund4) * 100 # percentage motile using TDI4
  organicAbund <- sum(TDI$Abundance[TDI$OrganicTolerant == TRUE],na.rm=TRUE)  # sum abundance of organic taxa
  organicPercent <- (organicAbund / allTaxaAbund4) * 100 # percentage organic using TDI4
  
  #### percentages TDI3:
  
  plankticAbund3 <- sum(TDI$Abundance[TDI$Planktic == TRUE],na.rm=TRUE) # count planktonic taxa
  plankticPercent3 <- (plankticAbund3 / allTaxaAbund3) * 100 # percentage planktonic using TDI3
  motileAbund3 <- sum(TDI$Abundance[TDI$Motile == TRUE],na.rm=TRUE) # sum abundance of motile taxa
  motilePercent3 <- (motileAbund3 / allTaxaAbund3) * 100 # percentage motile using TDI3
  organicAbund3 <- sum(TDI$Abundance[TDI$OrganicTolerant == TRUE],na.rm=TRUE)  # sum abundance of organic taxa
  organicPercent3 <- (organicAbund3 / allTaxaAbund3) * 100 # percentage organic using TDI3
  
  #### TDI3 & TDI4 & eTDI4 & eTDI3 scores:
  
  as <- sum(TDI$score,na.rm=TRUE)
  as3 <- sum(TDI$score3,na.rm=TRUE)
  TDIscore <- ((as / tdi4SumAbund)*25) - 25 
  TDIscore3 <- ((as3 / tdi3SumAbund)*25) - 25
  eTDI <-  9.933*exp(log10(unique(TDI$Alk))*0.81)
  eTDI3 <- -25.36+(56.83*(log10(unique(TDI$Alk))))-(12.96*(log10(unique(TDI$Alk))*log10(unique(TDI$Alk))))+(3.21*(unique(TDI$DaresSeason)))
  EQR <- (100-TDIscore)/(100-eTDI)
  DaresEQR <-  (100-TDIscore3)/(100-eTDI3)
  desc <- unique(TDI$SiteID)
  date <- unique(TDI$Date)
  samplenumber <- unique(TDI$SampleID)
  #loc <- unique(TDI$Loc.x)
  
  ### Loch LTDI2 scores:
  
  lochAbundsum <- sum(TDI$Abundance,na.rm=TRUE) # [TDI$LTDI2 > 0]) should this exclude zero scoring taxa? -
  LTDI4SumAbund <- sum(TDI$scoreL4,na.rm=TRUE)
  w <- LTDI4SumAbund / lochAbundsum
  lochTDI4 <- (w * 25) - 25
  lochEQR4 <- (100 - (lochTDI4)) / (100 - (unique(TDI$LochRefValue2)))
  
  ### Loch LTDI1 scores:
  
  totalabundsum <- sum(TDI$Abundance,na.rm=TRUE)
  lochAbundsum2 <- sum(TDI$Abundance[TDI$LTDI1 > 0],na.rm=TRUE) # should this exclude zero scoring taxa? -
  LTDI3SumAbund <- sum(TDI$scoreL3,na.rm=TRUE) # issue with duplicate taxa in NEMS
  w2 <- LTDI3SumAbund / lochAbundsum2
  lochTDI3 <- (w2 * 25) - 25
  lochEQR3 <- (100 - (lochTDI3)) / (100 - (unique(TDI$LochRefValue)))
  
  outDF <- data.frame(sumAbund =sumAbund,
                      tdi3SumAbund=tdi3SumAbund,
                      tdi4SumAbund=tdi4SumAbund,
                      allTaxaAbund=allTaxaAbund,
                      allTaxaAbund3= allTaxaAbund3,
                      allTaxaAbund4= allTaxaAbund4,
                      notPlankticAbund=notPlankticAbund,
                      plankticAbund=plankticAbund,
                      plankticPercent=plankticPercent,
                      motileAbund=motileAbund,
                      motilePercent=motilePercent,
                      organicAbund=organicAbund,
                      organicPercent=organicPercent,
                      plankticAbund3=plankticAbund3,
                      plankticPercent3=plankticPercent3,
                      motileAbund3=motileAbund3,
                      motilePercent3=motilePercent3,
                      organicAbund3=organicAbund3,
                      organicPercent3=organicPercent3,
                      as=as,
                      as3=as3,
                      TDIscore=TDIscore,
                      TDIscore3= TDIscore3,
                      eTDI=eTDI,
                      eTDI3=eTDI3, 
                      EQR=EQR,
                      DaresEQR=DaresEQR,
                      desc=desc,
                      date=date,
                      samplenumber=samplenumber,
                      lochAbundsum=lochAbundsum,
                      LTDI4SumAbund=LTDI4SumAbund,
                      w=w,
                      lochTDI4=lochTDI4,
                      lochEQR4=lochEQR4,
                      totalabundsum=totalabundsum,
                      lochAbundsum2=lochAbundsum2,
                      LTDI3SumAbund,
                      w2=w2,
                      lochTDI3=lochTDI3,
                      lochEQR3=lochEQR3
                      )
  
  return(outDF) ##runs output and stores under function name i.e. creates 'splitTDI2' as a list
  
 })
  
  dataf <- do.call("rbind",dataf) # combines lits into dataframe using useful column names we have created
  return(dataf) # pri
  
  }
  
