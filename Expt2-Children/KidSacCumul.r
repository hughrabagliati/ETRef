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

kid.sac <- sac.process2("./","Kids")
kid.sac$order <- 1:length(kid.sac$SacToTarg)
names <- read.delim("WriteNames-Times.txt")
Groups = read.delim("SubjNames.txt", header = T)
merge(kid.sac,Groups, all = T) -> kid.sac
#kid.sac[is.na(kid.sac$AgeGroup),]$AgeGroup <- "Old"
kid.sac <- merge(kid.sac,names, by = c("Subj","trialnum"), all.x = TRUE)
kid.sac<- kid.sac[order(kid.sac$order),]
names(kid.sac)[names(kid.sac) == "StartTime..ms."] <- "StartTime"
kid.sac <- ddply(kid.sac, .(Subj,trialnum,Period), transform, CumFromTarg = cumsum(SacFromTarg),CumToTarg = cumsum(SacToTarg),CumTarg = cumsum(SacTarg),CumD1 = cumsum(SacDist1),CumD2 = cumsum(SacDist2), SacTime = SacEndTime - min(SacEndTime), SacBin = round((SacEndTime - min(SacEndTime))/100))


# Create LabelCond variable; used for data analysis
kid.sac$LabelCond <- NA
kid.sac[kid.sac$cond == "Control" & kid.sac$Label %in% c(1,0),]$LabelCond <- "Control"
kid.sac[kid.sac$cond == "Ambig" & kid.sac$Label %in% c(1),]$LabelCond <- "Test-Ambig"
kid.sac[kid.sac$cond == "Ambig" & kid.sac$Label %in% c(0),]$LabelCond <- "Test-Unambig"
kid.sac$LabelCond <- as.factor(kid.sac$LabelCond)



ddply(kid.sac, .(Subj,trialnum,cond,Period,targname,LabelCond,Label,StartTime,AgeGroup,Lang), summarize, SacBin = c(0:181)) -> kid.sac.s
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


plot.sac.cum <- function(sac.bin, Peri){
summaryBy(CumTarg+CumD1+CumD2~LabelCond+SacBin, data = subset(sac.bin, Period == Peri), keep.names = T) -> k 
plot(k[k$LabelCond == "Control",]$SacBin, k[k$LabelCond == "Control",]$CumTarg, ylim = c(0,2), xlim = c(-10,40), pch = 21, col = "darkgray",bg = "darkgrey")
points(k[k$LabelCond == "Test-Ambig",]$SacBin, k[k$LabelCond == "Test-Ambig",]$CumTarg, ylim = c(0,2), pch = 21, col = "blue",bg = "blue")
points(k[k$LabelCond == "Test-Unambig",]$SacBin, k[k$LabelCond == "Test-Unambig",]$CumTarg, ylim = c(0,2), pch = 21, col = "red",bg = "red")
}
plot.sac.cum(subset(kid.sac.bin.naming,AgeGroup == "Young" & Lang == "Mon"),"Naming")
abline(v=(median(kid.sac$StartTime, na.rm = T)/100))

plot.sac.cum(kid.sac.bin,"Pre")
summary(lmer(CumTarg~LabelCond*SacBin + (1+LabelCond|Subj)+(1+LabelCond|Item), data = subset(kid.sac.bin.naming,Period == "Naming" & SacBin <= 20 & AgeGroup == "Young" & Lang == "Mon")))
 summary(lmer(CumTarg~LabelCond*SacBin + (1+LabelCond|Subj)+(1+LabelCond|Item), data = kid.sac.bin[kid.sac.bin$Period == "Pre",]))
 
 kid.sac$Start <- "Before"
 kid.sac[kid.sac$Period == "Naming",]$Start <- ifelse(kid.sac[kid.sac$Period == "Naming",]$SacTime < kid.sac[kid.sac$Period == "Naming",]$StartTime,"Before","After")
 #kid.sac[kid.sac$Period == "Pre",]$Start <- ifelse(kid.sac[kid.sac$Period == "Pre",]$SacTime < 2250 ,"Before","After")
  kid.sac[kid.sac$Period == "Rew",]$Start <- "After"
	 summaryBy(SacDist1+SacDist2+SacFromTarg+SacToTarg +SacTarg ~AgeGroup+Subj+trialnum+LabelCond+Period+Start, data = kid.sac, keep.names = T) -> Sac.sum
	 na.omit(summaryBy(SacTarg +SacDist2 +SacDist1+SacFromTarg+SacToTarg ~AgeGroup+Period+Start+LabelCond, data = Sac.sum, keep.names = T))
	 
	 summary(lmer(SacTarg~LabelCond+ (1|Subj), data = subset(Sac.sum, AgeGroup !="Excl" & Period == "Pre")))
	 	 summary(lmer(SacTarg~LabelCond+ (1|Subj), data = subset(Sac.sum, AgeGroup !="Excl" & Period == "Naming" & Start == "Before")))
	 	 summary(lmer(SacTarg~LabelCond+ (1|Subj), data = subset(Sac.sum, AgeGroup !="Excl" & Period == "Naming" & Start == "After")))	 	 
	
	 
	 
 	 summaryBy(SacDist1+SacDist2+SacFromTarg+SacToTarg +SacTarg ~AgeGroup+Subj+trialnum+LabelCond+Period+Start, data = kid.sac[kid.sac$SacSwitch ==1,], keep.names = T, FUN = c(sum)) -> Sac.sum2
 	 Sac.sum2$SacProp <- Sac.sum2$SacTarg/(Sac.sum2$SacTarg + Sac.sum2$SacDist2)
	 na.omit(summaryBy(SacTarg +SacDist2 +SacDist1+SacFromTarg+SacToTarg ~AgeGroup+Period+Start+LabelCond, data = Sac.sum2, keep.names = T))
	 
na.omit(summaryBy(SacTarg~Period+Start+LabelCond+Subj, data = Sac.sum[Sac.sum$Period != "Rew",], FUN = c(mean,sd), keep.names = T)) -> Sac.graph
na.omit(summaryBy(SacTarg.mean~Period+Start+LabelCond, data = Sac.graph, FUN = c(mean,sd))) -> Sac.graph
Sac.graph$Period <- factor(Sac.graph$Period, levels = c("Pre","Naming"),labels = c("Pre","Naming"), ordered = T)
Sac.graph$Start <- factor(Sac.graph$Start, levels = c("Before","After"),labels = c("Before","After"), ordered = T)
Sac.graph$SE = Sac.graph$SacTarg.mean.sd/sqrt(length(unique(Sac.sum$Subj)))

Sac.graph$Time = "Pre-Naming"
Sac.graph[Sac.graph$Start == "After" , ]$Time = "Post-Naming"
Sac.graph[Sac.graph$Period == "Pre" , ]$Time = "Preview"
Sac.graph$Time <- factor(Sac.graph$Time, levels = c("Preview","Pre-Naming","Post-Naming"),labels = c("Preview","Pre-Naming","Post-Naming"), ordered = T)

tapply(Sac.graph$SacTarg.mean.mean, list(Sac.graph$Time,Sac.graph$LabelCond), FUN = mean) -> o
tapply(Sac.graph$SE, list(Sac.graph$Time,Sac.graph$LabelCond), FUN = mean) -> se

barplot(o, beside =T , ylim = c(0,0.25),col = "white",  border = NA, ylab = "Proportion actions using Instrument", names.arg = c("Preview", "Pre-Naming","Post-Naming"))
 legend(1.2,0.10, legend = c("Control", "Ambiguous Description", "Unambiguous Description"), bty = "n", col = c("blue","grey","red"), pch = 20)
 points(c(1.5,6,10), o[1,], pch = 20, cex = 2, col = "blue")
 points(c(2.5,6.8,11), o[2,], pch = 20, cex = 2, col = "grey")
  points(c(3.5,7.6,12), o[3,], pch = 20, cex = 2, col = "red")
 grid(nx = NA, ny = NULL, col = "gray", lty = "dotted",lwd = par("lwd"), equilogs = TRUE)
abline(v = c(4.5,8.5), col = "grey", lty = "dashed")
  arrows(c(1.5,2.5,3.5,6,6.8,7.6,10,11,12), (c(o) + c(se)+0.01), c(1.5,2.5,3.5,6,6.8,7.6,10,11,12), (c(o) - c(se)-0.01), code = 0)	 
