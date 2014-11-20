require(lme4)
require(boot)
require(doBy)
require(plyr)
require(longitudinalData)
sac.process2 = function(pathway = "./", Pop = "NA"){
	sac = c()
	for (i in unique(Pop)){
		for (j in c("Pre","Targ","Rew")){# Something is off with the reward section ,"Rew")){
		file.list <- list.files(path = pathway,full.names = T,pattern = paste(i,".*","Sac",j, sep = ""))
		sac.temp <- read.delim(file.list, header = T)
		sac.temp$Period = ifelse(j == "Targ","Naming",j)
		sac <- rbind(sac, sac.temp)	
		}
	}
	names(sac)[names(sac) == "RECORDING_SESSION_LABEL"] <- "Subj"
	names(sac)[names(sac) == "CURRENT_SAC_NEAREST_START_INTEREST_AREA_LABEL"] <- "SacStart"
	names(sac)[names(sac) == "CURRENT_SAC_NEAREST_END_INTEREST_AREA_LABEL"] <- "SacEnd"
	names(sac)[names(sac) == "CURRENT_SAC_END_TIME"] <- "SacEndTime"
	
	sac$Subj <- as.factor(paste(sac$Subj,".edf", sep = ""))
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
	#sac <- summaryBy(SacTarg + SacSwitch + Sac + SacDist1 + SacDist2 ~ Subj + trialnum + cond+Period, data = sac, keep.names = T)
	sac <- sac[sac$type != "Filler",]
	return(sac)
}

kid.sac <- sac.process2("./Eyedata/","Kids")
kid.sac$order <- 1:length(kid.sac$SacToTarg)
names <- read.delim("./Eyedata/WriteNames-Times.txt")
#Groups = read.delim("SubjNames.txt", header = T)
#merge(kid.sac,Groups, all = T) -> kid.sac
#kid.sac[is.na(kid.sac$AgeGroup),]$AgeGroup <- "Old"
kid.sac <- merge(kid.sac,names, by = c("Subj","trialnum"), all.x = TRUE)
kid.sac<- kid.sac[order(kid.sac$order),]
names(kid.sac)[names(kid.sac) == "StartTime..ms."] <- "StartTime"
kid.sac <- ddply(kid.sac, .(Subj,trialnum,Period), transform, CumFromTarg = cumsum(SacFromTarg),CumToTarg = cumsum(SacToTarg),CumTarg = cumsum(SacTarg),CumD1 = cumsum(SacDist1),CumD2 = cumsum(SacDist2), SacTime = SacEndTime - min(SacEndTime), SacBin = round((SacEndTime - min(SacEndTime))/100))
save(kid.sac,file = "kidmonitor.RDATA")

load('kidmonitor.RDATA')
# Create LabelCond variable; used for data analysis
kid.sac$LabelCond <- NA
kid.sac[kid.sac$cond == "Control" & kid.sac$Label %in% c(1,0),]$LabelCond <- "Control"
kid.sac[kid.sac$cond == "Ambig" & kid.sac$Label %in% c(1),]$LabelCond <- "Test-1Hit"
kid.sac[kid.sac$cond == "Ambig" & kid.sac$Label %in% c(0),]$LabelCond <- "Test-0Miss"
kid.sac$LabelCond <- as.factor(kid.sac$LabelCond)



ddply(kid.sac, .(Subj,trialnum,cond,Period,targname,LabelCond,Label,StartTime), summarize, SacBin = c(0:181)) -> kid.sac.s
ddply(kid.sac[kid.sac$SacBin <182,], .(Subj,trialnum,cond,Period,targname,LabelCond,Label,SacBin), summarize, CumTarg = mean(CumTarg),CumD1 = mean(CumD1),CumD2 = mean(CumD2)) -> kid.sac.sum
kid.sac.bin <- merge(kid.sac.s, kid.sac.sum, all = TRUE)
kid.sac.bin$CumTarg <- t(imputation(matrix(kid.sac.bin$CumTarg, nrow = 1),method = "locf"))
kid.sac.bin$CumD1 <- t(imputation(matrix(kid.sac.bin$CumD1, nrow = 1),method = "locf"))
kid.sac.bin$CumD2 <- t(imputation(matrix(kid.sac.bin$CumD2, nrow = 1),method = "locf"))
kid.sac.bin$Item <- as.factor(sapply(strsplit(as.character(kid.sac.bin$targname),"[.12]"), "[", 1))

kid.sac.bin$SubjTrial <- paste(kid.sac.bin$Subj,kid.sac.bin$Item, sep = "") 
TrialTime <- summaryBy(TRIAL_DWELL_TIME~SubjTrial, data = subset(kid.ref.t, subset = Period == "Pre" & Picture == "Targ"), keep.names = T, na.rm = T)
kid.sac.bin <- kid.sac.bin[kid.sac.bin$SubjTrial %in% TrialTime[TrialTime$TRIAL_DWELL_TIME > 1500,]$SubjTrial,]



kid.sac.bin$StartTime <- kid.sac.bin$StartTime/100
kid.sac.bin$StartTime <- round(kid.sac.bin$StartTime)
kid.sac.bin$SacBin2 <- kid.sac.bin$SacBin - kid.sac.bin$StartTime
#Change the SacBin2 above to SacBin to see this timelocked to naming onset.

kid.sac.bin.naming <- subset(kid.sac.bin, SacBin2 >= -10 & SacBin2 <= 50)
kid.sac.bin.naming$SacBin <- kid.sac.bin.naming$SacBin2
kid.sac.bin.naming <- ddply(kid.sac.bin.naming, .(Subj,trialnum, Period), transform, CumTarg = CumTarg - min(CumTarg),CumD1 = CumD1 - min(CumD1),CumD2 = CumD2 - min(CumD2))

# We don't have naming times currently
summary(lmer(SacTarg~LabelCond+ (1|Subj), data = subset(Sac.sum,  Period == "Pre")))
	 
summaryBy(SacDist1+SacDist2+SacFromTarg+SacToTarg +SacTarg ~Subj+trialnum+LabelCond+Period, data = kid.sac, keep.names = T) -> Sac.sum
na.omit(summaryBy(SacTarg +SacDist2 +SacDist1+SacFromTarg+SacToTarg ~Period+LabelCond, data = Sac.sum, keep.names = T))


na.omit(summaryBy(SacTarg~Period+LabelCond+Subj, data = Sac.sum[Sac.sum$Period != "Rew",], FUN = c(mean,sd), keep.names = T)) -> Sac.graph
na.omit(summaryBy(SacTarg.mean~Period+LabelCond, data = Sac.graph, FUN = c(mean,sd))) -> Sac.graph
Sac.graph$Period <- factor(Sac.graph$Period, levels = c("Pre","Naming"),labels = c("Pre","Naming"), ordered = T)
Sac.graph$SE = Sac.graph$SacTarg.mean.sd/sqrt(length(unique(Sac.sum$Subj)))

Sac.graph$Time = "Naming"
Sac.graph[Sac.graph$Period == "Pre" , ]$Time = "Preview"
Sac.graph$Time <- factor(Sac.graph$Time, levels = c("Preview","Naming"),labels = c("Preview","Naming"), ordered = T)

tapply(Sac.graph$SacTarg.mean.mean, list(Sac.graph$LabelCond,Sac.graph$Time), FUN = mean) -> o
tapply(Sac.graph$SE, list(Sac.graph$LabelCond,Sac.graph$Time), FUN = mean) -> se

barplot(o, beside =T , ylim = c(0,0.40),col = "white",  border = NA, ylab = "Critical Saccades", names.arg = c("Preview", "Naming"))
 legend(1.2,0.15, legend = c("Control", "Not Detecting Ambiguity", "Detecting Ambiguity"), bty = "n", col = c("blue","grey","red"), pch = 20)
 points(c(1.5,6), o[1,], pch = 20, cex = 2, col = "blue")
 points(c(2.5,6.8), o[2,], pch = 20, cex = 2, col = "grey")
  points(c(3.5,7.6), o[3,], pch = 20, cex = 2, col = "red")
 grid(nx = NA, ny = NULL, col = "gray", lty = "dotted",lwd = par("lwd"), equilogs = TRUE)
abline(v = c(4.5,8.5), col = "grey", lty = "dashed")
  arrows(c(1.5,2.5,3.5,6,6.8,7.6), (c(o) + c(se)+0.01), c(1.5,2.5,3.5,6,6.8,7.6), (c(o) - c(se)-0.01), code = 0)	 