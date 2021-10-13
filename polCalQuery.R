#-------- Parse arguments
parseArg <- function( args ){
	argNum <- length(args)
	argList <- list(as.character(as.POSIXct(Sys.time())), 3600, 3, 0.05, 100.0, FALSE)
	names(argList) <- c('startTime', 'execDuration', 'Band', 'threshFlux', 'refFreq', 'load')
	for( index in 1:argNum ){
		if(substr(args[index], 1,2) == "-s"){ argList$startTime <- as.character(substring(args[index], 3)) }
		if(substr(args[index], 1,2) == "-d"){ argList$execDuration <- as.numeric(substring(args[index], 3)) }
		if(substr(args[index], 1,2) == "-b"){ argList$Band <- as.integer(substring(args[index], 3)) }
		if(substr(args[index], 1,2) == "-t"){ argList$threshFlux <- as.numeric(substring(args[index], 3))}
		if(substr(args[index], 1,2) == "-#"){ argList$refFreq <- as.numeric(substring(args[index], 3))}
		if(substr(args[index], 1,2) == "-L"){ argList$load <- TRUE}
	}
	if( !file.exists("Flux.Rdata") ){ argList$load <- TRUE }
	return(argList)
}
#-------- Time Constants
SEC_PER_DAY <- 86400
MJD_1901 <- 15384
DAY_PER_YEAR <- 365
DAY_PER_4YEAR <- 1461
DOY_MON <- c(0, 31, 59, 90, 120, 151, 181, 212, 242, 273, 303, 334)
#-------- FE-specific PA
#         Band1      2     3      4     5      6     7      8      9   10
BandPA <- c(45.0, -45.0, 80.0, -80.0, 45.0, -45.0, 36.45, 90.0, -90.0, 0.0)*pi/180
BandFreq <- c(43.0, 75.0, 97.5, 132.0, 183.0, 233.0, 343.5, 400.0, 650.0, 800.0)
ALMA_LAT <- -23.029* pi/180.0
ALMA_LONG <- -67.755* pi/180.0
cos_phi <- cos(ALMA_LAT)
sin_phi <- sin(ALMA_LAT)
maxSinEL <- sin(86/180*pi)
minSinEL <- sin(20/180*pi)

#-------- Calculate Day of year from Month and date
md2doy <- function(year, month, date){
	is_leap <- ((year%%4 == 0) && (month > 2))	# Leap Year Frag
	DOY_MON[month] + date + is_leap
}

#-------- Calculate (fractional) Modified Julian Date in unit of second. This unit is used in CASA
doy2mjdSec <- function(year, doy, hour, min, sec){
	if((year < 1901) || (year > 2099)){ return(-1)}
	year <- year - 1901
	mjd <- year %/%4 * DAY_PER_4YEAR + year%%4 * DAY_PER_YEAR + doy + MJD_1901
	sod <- (hour*60 + min)*60 + sec
	return(mjd* SEC_PER_DAY + sod)
}

#-------- Convert ISO8601 format string into mjdSec
ISO8601mjdSec <- function( timeString ){				# Input string YYYY-MM-DDTHH:MM:SS.S
	year <- as.integer(substring(timeString, 1, 4))
	month <- as.integer(substring(timeString, 6, 7))
	day <- as.integer(substring(timeString, 9, 10))
	UT_hour <- as.integer(substring(timeString, 12, 13))
	minute <- as.integer(substring(timeString, 15, 16))
	second <- as.integer(substring(timeString, 18, 19))
	return(doy2mjdSec(year, md2doy(year, month, day), UT_hour, minute, second))
}

#-------- MJD to Greenwich mean sidereal time
mjd2gmst <- function(mjd, ut1utc = 0){
	# mjd2gmst : Calculate Greenwidge Mean Sidereal Time
	# mjd : Modified Julian Date
	# ut1utc : UT1 - UTC [sec]
	
	FACT <- c(24110.54841, 8640184.812866, 0.093104, 0.0000062)
	MJD_EPOCH <- 51544.5  # MJD at 2000 1/1 12:00:00 UT
	TU_UNIT <- 36525.0
	SEC_PER_DAY   <- 86400
	
	tu <- (mjd - MJD_EPOCH) / TU_UNIT
	ut1 <- (mjd%%1)* SEC_PER_DAY + ut1utc
	gmst <- (ut1 + FACT[1] + ((FACT[4]* tu + FACT[3])* tu + FACT[2])* tu) / SEC_PER_DAY
	return(2* pi* (gmst %% 1))
}


#-------- # cos(hour angle) when it passes the given EL
EL_HA <- function(sinEL, dec){
	cosHA <- (sinEL - sin_phi* sin(dec)) / (cos_phi* cos(dec))
	return(cosHA)
}

#-------- # filtering by source polarization
sourceDataFrame <- function(DF, refFreq=100.0, refDate=Sys.Date()){
    DateRange <- 30    # 30-day window
    #-------- Filter by observing date
	DF <- DF[abs(as.Date(DF$Date) - as.Date(refDate)) < DateRange,]
	filterDF <- DF[((DF$Band == 7) | (DF$Band == 6)),]
    sourceList <- unique(filterDF$Src)
    sourceList <- sourceList[grep('^J[0-9]', sourceList)]  # Filter SSOs out
    #-------- filter by observed band coverage
    flagSoruceList <- sourceList
    for(src in sourceList){
        srcDF <- DF[DF$Src == src,]
        bandList <- unique(srcDF$Band)
        B3obsNum <- nrow(srcDF[srcDF$Band == 3,])
        B7obsNum <- nrow(srcDF[((srcDF$Band == 7) | (srcDF$Band == 6)),])
        if( (length(bandList) < 2) | (B3obsNum < 3) | (B7obsNum < 3) ){
            flagSoruceList <- flagSoruceList[-which(flagSoruceList == src)]
		}
	}
    sourceList <- flagSoruceList
    numSrc <- length(sourceList)
    #-------- Source properties
    RAList <- (60.0* as.numeric(substring(sourceList, 2,3)) + as.numeric(substring(sourceList, 4,5))) / 720 * pi # RA in [rad]
    DecList<- as.numeric(substring(sourceList, 6,8))
    DecList<- DecList + sign(as.numeric(substring(sourceList, 6,10)))* as.numeric(substring(sourceList, 9,10))/60.0
    DecList<- DecList / 180 * pi # DEC in [rad]
	sourceDF <- data.frame(Src=sourceList, RA=RAList, Dec=DecList, I=numeric(numSrc), Q=numeric(numSrc), U=numeric(numSrc), V=numeric(numSrc), P=numeric(numSrc), EVPA=numeric(numSrc), RM=numeric(numSrc), stringsAsFactors=F)
	for(src in sourceList){
		index <- which(sourceDF$Src == src)
		srcDF <- DF[((DF$Src == src) & (abs(as.Date(DF$Date) - as.Date(refDate)) < DateRange)),] 
		srcDF$relTime <-as.numeric(srcDF$Date) - as.numeric(as.POSIXct(refDate))  # Relative time since reference [sec]
		srcDF$relFreq <- srcDF$Freq / refFreq                                	  # Relative frequency ratio
		fitI <- lm(formula=log(I) ~ relTime + log(relFreq), data=srcDF, weight=(srcDF$I/srcDF$eI)^2)
		fitP <- lm(formula=log(P) ~ relTime + log(relFreq), data=srcDF, weight=(srcDF$P/(srcDF$eP_upper - srcDF$eP_lower))^2)
		B3DF <- srcDF[srcDF$Band == 3,]
		B7DF <- srcDF[((srcDF$Band == 7) | (srcDF$Band == 6)),]
		B7DF[B7DF$Band == 6,]$I <- B7DF[B7DF$Band == 6,]$I* (BandFreq[7]/BandFreq[6])^coef(fitI)['log(relFreq)'][[1]]
		B7DF[B7DF$Band == 6,]$P <- B7DF[B7DF$Band == 6,]$P* (BandFreq[7]/BandFreq[6])^coef(fitP)['log(relFreq)'][[1]]
		B7DF[B7DF$Band == 6,]$Q <- B7DF[B7DF$Band == 6,]$Q* (BandFreq[7]/BandFreq[6])^coef(fitP)['log(relFreq)'][[1]]
		B7DF[B7DF$Band == 6,]$U <- B7DF[B7DF$Band == 6,]$U* (BandFreq[7]/BandFreq[6])^coef(fitP)['log(relFreq)'][[1]]
		B7DF[B7DF$Band == 6,]$V <- B7DF[B7DF$Band == 6,]$V* (BandFreq[7]/BandFreq[6])^coef(fitP)['log(relFreq)'][[1]]		
		fitB3Q <- lm(formula = Q ~ relTime, data=B3DF, weight=(B3DF$P / B3DF$eQ)^2)
		fitB3U <- lm(formula = U ~ relTime, data=B3DF, weight=(B3DF$P / B3DF$eU)^2)
		fitB3V <- lm(formula = V ~ relTime, data=B3DF, weight=(B3DF$P / B3DF$eV)^2)
		fitB7Q <- lm(formula = Q ~ relTime, data=B7DF, weight=(B7DF$P / B7DF$eQ)^2)
		fitB7U <- lm(formula = U ~ relTime, data=B7DF, weight=(B7DF$P / B7DF$eU)^2)
		fitB7V <- lm(formula = V ~ relTime, data=B7DF, weight=(B7DF$P / B7DF$eV)^2)
		sourceDF[index,]$V <- (coef(fitB3V)['(Intercept)'][[1]]*(BandFreq[7] - refFreq) + coef(fitB7V)['(Intercept)'][[1]]*(refFreq - BandFreq[3])) / (BandFreq[7] - BandFreq[3])		
		B3EVPA <- 0.5* atan2( coef(fitB3U)['(Intercept)'][[1]], coef(fitB3Q)['(Intercept)'][[1]] )
		B7EVPA <- 0.5* atan2( coef(fitB7U)['(Intercept)'][[1]], coef(fitB7Q)['(Intercept)'][[1]] )
		lambdaSQ <- (0.299792458 / c(BandFreq[3], BandFreq[7], refFreq))^2
		if(B3EVPA - B7EVPA >  pi/2){ B3EVPA <- B3EVPA - pi }
		if(B3EVPA - B7EVPA < -pi/2){ B3EVPA <- B3EVPA + pi }
		sourceDF[index,]$RM <- (B7EVPA - B3EVPA) / (lambdaSQ[2] - lambdaSQ[1])
		sourceDF[index,]$EVPA <- (lambdaSQ[2]* B3EVPA - lambdaSQ[1]* B7EVPA) / (lambdaSQ[2] - lambdaSQ[1]) + sourceDF[index,]$RM* lambdaSQ[3]
		sourceDF[index,]$I <- exp(coef(fitI)['(Intercept)'][[1]])
		sourceDF[index,]$P <- exp(coef(fitP)['(Intercept)'][[1]])
	}
	sourceDF$Q <- sourceDF$P * cos(2.0* sourceDF$EVPA)
	sourceDF$U <- sourceDF$P * sin(2.0* sourceDF$EVPA)
	return( sourceDF )
}

#-------- Input parameters for debugging
#Arguments <- "-s2019-12-13T03:10:21"
Arguments <- commandArgs(trailingOnly = T)
setwd('./')
argList <- parseArg(Arguments)
startmjdSec <- ISO8601mjdSec(argList$startTime)
endmjdSec <- startmjdSec + argList$execDuration
startLST <- mjd2gmst(startmjdSec/SEC_PER_DAY) + ALMA_LONG
endLST   <- mjd2gmst(endmjdSec/SEC_PER_DAY) + ALMA_LONG
if( endLST < startLST){ endLST <- endLST + 2*pi}

#-------- Load Flux.Rdata from web
if( argList$load ){
	cat('--- Loading Flux.Rdata from the web\n')
	FluxDataURL <- "https://www.alma.cl/~skameno/Grid/Stokes/"
	load(url(paste(FluxDataURL, "Flux.Rdata", sep='')))     # Data frame of FLDF
    save(FLDF, file='Flux.Rdata')
} else {
	load('Flux.Rdata')
}
#-------- Filter out single-sideband Band-3 frequency 
pos <- regexpr("RB",FLDF$File)
FLDF$Band <- as.integer(substr(FLDF$File, pos+3, pos+4))
FLDF$BandPA <- BandPA[FLDF$Band]
FLDF <- FLDF[-which((FLDF$Band == 3) & (abs(FLDF$Freq - 97.45) > 1.0)),]

SDF <- sourceDataFrame( FLDF, argList$refFreq, argList$startTime)

#---- Filter by P > 0.05 Jy
SDF <- SDF[SDF$P > argList$threshFlux,]

#---- Filter by m > 3%
SDF <- SDF[SDF$P / SDF$I > 0.03,]

#---- Filter by EL
SDF$startHA <- startLST - SDF$RA
SDF$endHA   <- endLST - SDF$RA
SDF <- SDF[ which(cos(SDF$startHA) > EL_HA(minSinEL, SDF$Dec)), ]	# EL >20º at the beggining
SDF <- SDF[ which(cos(SDF$endHA)   > EL_HA(minSinEL, SDF$Dec)), ] # EL >20º at the end
SDF <- SDF[ which(cos(SDF$startHA) < EL_HA(maxSinEL, SDF$Dec)), ] # EL < 86º at the beggining
SDF <- SDF[ which(cos(SDF$endHA)   < EL_HA(maxSinEL, SDF$Dec)), ] # EL < 86º at the end
SDF <- SDF[ which( (EL_HA(maxSinEL, SDF$Dec) > 1.0) | (sin(SDF$startHA)* sin(SDF$endHA) > 0)),] # transit for EL>86º

#---- Calculate feed-EVPA angle
sourceList <- SDF$Src
SDF$QC_H <- SDF$QC_L <- SDF$QC_0 <- SDF$UC_H <- SDF$UC_L <- SDF$UC_0 <- numeric(length(sourceList))
# H <- L <- numeric(length(sourceList))
for(src in sourceList){
	index <- which(SDF$Src == src)
	HA <- seq(SDF$startHA[index], SDF$endHA[index], by=0.01)
	cos_HA  <- cos(HA)
	sin_HA  <- sin(HA)
	cos_dec <- cos(SDF$Dec[index])
	sin_dec <- sin(SDF$Dec[index])
	PA <- atan2(sin_HA, sin_phi*cos_dec/cos_phi - sin_dec* cos_HA) + BandPA[argList$Band]
	QCpUS <- SDF$Q[index] * cos(2.0* PA) + SDF$U[index] * sin(2.0* PA)
	UCmQS <- SDF$U[index] * cos(2.0* PA) - SDF$Q[index] * sin(2.0* PA)
	SDF$UC_0[index] <- HA[which.min( abs(UCmQS) )] + SDF$RA[index]
	SDF$QC_0[index] <- HA[which.min( abs(QCpUS) )] + SDF$RA[index]
	SDF$QC_H[index] <- max(QCpUS)
	SDF$QC_L[index] <- min(QCpUS)
	SDF$UC_H[index] <- max(UCmQS)
	SDF$UC_L[index] <- min(UCmQS)
}

calDF <- SDF[which.max(SDF$UC_H - SDF$UC_L),]		# Primary calibrator for max XY
SDF <- SDF[-(SDF$Src == calDF$Src),]
if( calDF$UC_H* calDF$UC_L > 0 ){					# require 2nd calibrator
	calDF[2,] <- SDF[which.min(SDF$UC_H* SDF$UC_L),]	# 2nd calibrator for XY intercept
}
if( min(calDF$QC_L* calDF$QC_H) > 0 ){				# 3rd calibrator for gain equalization (QCpUS = 0)
	index <- which(SDF$QC_L * SDF$QC_H < 0)
	if(length(index) > 0){
		calDF[nrow(calDF)+1,] <- SDF[which.max(SDF$P),]
	}
}
write.table(na.omit(calDF[,c('Src', 'I', 'P', 'EVPA', 'UC_L', 'UC_H', 'UC_0', 'QC_0')]), file="CalQU.data", append=F, quote=F, col.names=T, row.name=F)
