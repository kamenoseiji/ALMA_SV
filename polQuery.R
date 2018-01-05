Arguments <- commandArgs(trailingOnly = T)
timeWindow <- 60	# Days
library(RCurl)
#-------- Function to return residuals in RM fit
residEVPA <- function(x, y, w){	# Optimization function for RM fit
	return(function(para){
		RM <- para[1]
		EVPA <- para[2]
		return( sum(w* sin(y - RM*x - EVPA)^2 ))
	})
}

#-------- Parse arguments
parseArg <- function( args ){
	srcNum <- argNum <- length(args)
	for( index in 1:argNum ){
		if(substr(args[index], 1,2) == "-D"){ refDate <- as.Date(substring(args[index], 3));    srcNum <- srcNum - 1 }
		if(substr(args[index], 1,2) == "-F"){ refFreq <- as.numeric(substring(args[index], 3)); srcNum <- srcNum - 1}
	}
	srcList = args[(argNum - srcNum + 1):argNum]
	return(list(refDate = refDate, refFreq = refFreq, srcList = srcList[grep('^J[0-9]',srcList )]))
}

argList <- parseArg(Arguments)
refDate <- argList$refDate
refFreq <- argList$refFreq
srcList <- argList$srcList

#-------- Load Flux.Rdata from web
load(url("http://www.alma.cl/~skameno/Grid/Stokes/Flux.Rdata")) 
FLDF$timeDiff <- as.numeric(difftime(FLDF$Date, refDate, units='days'))

#-------- For each source
IatRef <- QatRef <- UatRef <- numeric(0)
for(sourceName in srcList){
	srcDF <- FLDF[((FLDF$Src == sourceName) & (abs(FLDF$timeDiff) < timeWindow)),]
	if(nrow(srcDF) < 3){ next }
	freqList <- as.numeric(unique(srcDF$Freq))
	freqNum <- length(freqList)
	freqList <- freqList[order(freqList)]
	estI <- errI <- estQ <- errQ <- estU <- errU <- numeric(freqNum)
	#-------- For each frequency
	for(freq_index in 1:freqNum){
		srcFreqDF <- srcDF[srcDF$Freq == freqList[freq_index],]
		if(nrow(srcFreqDF) < 3){
			estI[freq_index] <- median(srcFreqDF$I); errI[freq_index] <- median(srcFreqDF$eI) * 10.0
			estQ[freq_index] <- median(srcFreqDF$Q); errQ[freq_index] <- median(srcFreqDF$eQ) * 10.0
			estU[freq_index] <- median(srcFreqDF$U); errU[freq_index] <- median(srcFreqDF$eU) * 10.0
		} else {
			fit <- lm(data=srcFreqDF, formula=I ~ timeDiff, weights=1.0 / eI^2 / abs(timeDiff + 1))
			estI[freq_index] <- summary(fit)$coefficients[1,'Estimate'];  errI[freq_index] <- summary(fit)$coefficients[1,'Std. Error']
			fit <- lm(data=srcFreqDF, formula=Q ~ timeDiff, weights=1.0 / eQ^2 / abs(timeDiff + 1))
			estQ[freq_index] <- summary(fit)$coefficients[1,'Estimate'];  errQ[freq_index] <- summary(fit)$coefficients[1,'Std. Error']
			fit <- lm(data=srcFreqDF, formula=U ~ timeDiff, weights=1.0 / eU^2 / abs(timeDiff + 1))
			estU[freq_index] <- summary(fit)$coefficients[1,'Estimate'];  errU[freq_index] <- summary(fit)$coefficients[1,'Std. Error']
		}
	}
	lambdaSQ <- (0.299792458 / freqList)^2; lambdasqSpan <- diff(range(lambdaSQ))
	estP <- sqrt(estQ^2 + estU^2); errP <- sqrt(errQ^2 + errU^2); estEVPA <- 0.5*atan2(estU, estQ)
	fit <- lm(log(estI) ~ log(freqList/100.0), weights=1.0/errI^2); I100 <- exp(as.numeric(coef(fit)[1])); spixI <- as.numeric(coef(fit)[2])
	fit <- lm(log(estP) ~ log(freqList/100.0), weights=1.0/errP^2); P100 <- exp(as.numeric(coef(fit)[1])); spixP <- as.numeric(coef(fit)[2])
	estEVPAend <- estEVPA[freqNum]
	if(estEVPAend - estEVPA[1] >  pi/2){ estEVPAend <- estEVPAend - pi }
	if(estEVPAend - estEVPA[1] < -pi/2){ estEVPAend <- estEVPAend + pi }
	RMinit <- (estEVPAend - estEVPA[1]) / (lambdaSQ[freqNum] - lambdaSQ[1])
	fit <- optim(par=c(RMinit, estEVPA[freqNum]), fn=residEVPA(lambdaSQ, estEVPA, 1.0/errP^2), method='Nelder-Mead')
	RM <- fit$par[1]; EVPAintercept <- fit$par[2]
	PatFreqRef <- P100*(refFreq/100)^spixP
	EVPAatFreqRef <- RM* (0.299792458 / refFreq)^2 + EVPAintercept
	IatRef <- append(IatRef, I100*(refFreq/100)^spixI)
	QatRef <- append(QatRef, PatFreqRef * cos(2.0* EVPAatFreqRef))
	UatRef <- append(UatRef, PatFreqRef * sin(2.0* EVPAatFreqRef))
}
DF <- data.frame(Src=srcList, I=IatRef, Q=QatRef, U=UatRef)
write.table(DF, file="CalQU.data", append=F, quote=F, col.names=F, row.name=F)
