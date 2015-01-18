require(plyr)
require(longitudinalData)
sac.process2 = function(pathway = "./", Pop = "NA"){
	sac = c()
	for (i in unique(Pop)){
		for (j in c("Pre","Targ","Rew")){# Something is off with the reward section ,"Rew")){
		file.list <- list.files(path = pathway,full.names = T,pattern = paste(i,".*","Sac",j, sep = ""))
		print(file.list)
		sac.temp <- read.delim(file.list, header = T)
		sac.temp$Period = ifelse(j == "Targ","Naming",j)
		sac <- rbind(sac, sac.temp)	
		}
	}
	names(sac)[names(sac) == "RECORDING_SESSION_LABEL"] <- "Subj"
	names(sac)[names(sac) == "CURRENT_SAC_NEAREST_START_INTEREST_AREA_LABEL"] <- "SacStart"
	names(sac)[names(sac) == "CURRENT_SAC_NEAREST_END_INTEREST_AREA_LABEL"] <- "SacEnd"
	names(sac)[names(sac) == "CURRENT_SAC_END_TIME"] <- "SacEndTime"
	sac$Sac <- 1
	sac$SacSwitch <- ifelse(sac$SacStart == sac$SacEnd,0,1)
	sac$SacToTarg <- 0
	sac$SacFromTarg <- 0
	sac[sac$SacStart == "Pre_Targ " & sac$SacEnd == "Pre_D1 ",]$SacFromTarg <- 1
	sac[sac$SacStart == "Pre_D1 " & sac$SacEnd == "Pre_Targ ",]$SacToTarg <- 1
	sac$SacTarg <- 0
	sac[sac$SacStart == "Pre_Targ " & sac$SacEnd == "Pre_D1 ",]$SacTarg <- 1
	sac[sac$SacStart == "Pre_D1 " & sac$SacEnd == "Pre_Targ ",]$SacTarg <- 1

	sac$SacDist1 <- 0
	sac[sac$SacStart == "Pre_D1 " & sac$SacEnd == "Pre_D2 ",]$SacDist1 <- 1
	sac[sac$SacStart == "Pre_D2 " & sac$SacEnd == "Pre_D1 ",]$SacDist1 <- 1
	sac$SacDist2 <- 0
	sac[sac$SacStart == "Pre_Targ " & sac$SacEnd == "Pre_D2 ",]$SacDist2 <- 1
	sac[sac$SacStart == "Pre_D2 " & sac$SacEnd == "Pre_Targ ",]$SacDist2 <- 1
	sac <- sac[sac$type != "Filler",]
	sac$cond <- factor(sac$cond)
  contrasts(sac$cond)[1] <- -1
	
	return(sac)
}

ad.sac <- sac.process2("./EyeData/",c("AdSil"))
ad.sac$order <- 1:length(ad.sac$SacTarg)
names <- read.delim("./EyeData/AdultWriteNames.txt")
ad.sac <- merge(ad.sac,names, by = c("Subj","trialnum"), all.x = TRUE)
ad.sac<- ad.sac[order(ad.sac$order),]

ad.sac <- ddply(ad.sac, .(Subj,trialnum,Period), transform, CumTarg = cumsum(SacTarg),CumD1 = cumsum(SacDist1),CumD2 = cumsum(SacDist2), SacTime = SacEndTime - min(SacEndTime), SacBin = round((SacEndTime - min(SacEndTime))/100))

ad.sac$Item <- as.factor(sapply(strsplit(as.character(ad.sac$targname),"[.12]"), "[", 1))


summaryBy(SacDist1+SacDist2+SacTarg ~Subj+Item+cond+type+Period, data = ad.sac, keep.names = T) -> Sac.sum
Sac.sum$SacOnTrial = ifelse(Sac.sum$SacTarg > 0,1,0)
na.omit(summaryBy(SacTarg +SacDist2 + SacOnTrial~Period+type+cond, data = Sac.sum, keep.names = T, FUN = c(mean,sd)))

summary(lmer(SacTarg~cond + (1+cond|Subj)+(1+cond|Item), data = subset(Sac.sum,  Period == "Pre")))