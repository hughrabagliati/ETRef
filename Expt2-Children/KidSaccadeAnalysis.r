require(lme4)
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
CreateNewDataFrame <- ginput("Press n to load the current dataset, or y to recreate the dataframe: ")

if (CreateNewDataFrame == "y"){
kid.sac <- sac.process2("./EyeData/","Kids")
  kid.sac$order <- 1:length(kid.sac$SacToTarg)
  names <- read.delim("./EyeData/WriteNames-Times.txt")
  Groups = read.delim("./EyeData/SubjNames.txt", header = T)
  merge(kid.sac,Groups, all = T) -> kid.sac
  #kid.sac[is.na(kid.sac$AgeGroup),]$AgeGroup <- "Old"
  kid.sac <- merge(kid.sac,names, by = c("Subj","trialnum"), all.x = TRUE)
  kid.sac<- kid.sac[order(kid.sac$order),]
  names(kid.sac)[names(kid.sac) == "StartTime..ms."] <- "StartTime"
  kid.sac <- ddply(kid.sac, .(Subj,trialnum,Period), transform, CumFromTarg = cumsum(SacFromTarg),CumToTarg = cumsum(SacToTarg),CumTarg = cumsum(SacTarg),CumD1 = cumsum(SacDist1),CumD2 = cumsum(SacDist2), SacTime = SacEndTime - min(SacEndTime), SacBin = round((SacEndTime - min(SacEndTime))/100))
  save(kid.sac,"kid.sac.RDATA")
}
load("kidsac.RDATA")

# Create LabelCond variable; used for data analysis
kid.sac$LabelCond <- NA
kid.sac[kid.sac$cond == "Control" & kid.sac$Label %in% c(1,0),]$LabelCond <- "Control"
kid.sac[kid.sac$cond == "Ambig" & kid.sac$Label %in% c(1),]$LabelCond <- "Test-Ambig"
kid.sac[kid.sac$cond == "Ambig" & kid.sac$Label %in% c(0),]$LabelCond <- "Test-Unambig"
kid.sac$LabelCond <- as.factor(kid.sac$LabelCond)

# Exclude kids with no utterance data
kid.sac <- kid.sac[!is.na(kid.sac$LabelCond),]
kid.sac <- kid.sac[!kid.sac$AgeGroup == "Excl",]
kid.sac <- kid.sac[!kid.sac$Lang == "Exc",]

# Discern when the name was said.
 kid.sac$Start <- "Before"
 kid.sac[kid.sac$Period == "Naming",]$Start <- ifelse(kid.sac[kid.sac$Period == "Naming",]$SacTime < kid.sac[kid.sac$Period == "Naming",]$StartTime,"Before","After")
 #kid.sac[kid.sac$Period == "Pre",]$Start <- ifelse(kid.sac[kid.sac$Period == "Pre",]$SacTime < 2250 ,"Before","After")
  kid.sac[kid.sac$Period == "Rew",]$Start <- "After"
	 summaryBy(SacDist1+SacDist2+SacFromTarg+SacToTarg +SacTarg ~AgeGroup+Subj+trialnum+LabelCond+Period+Start, data = kid.sac, keep.names = T) -> Sac.sum
	 na.omit(summaryBy(SacTarg +SacDist2 +SacDist1+SacFromTarg+SacToTarg ~AgeGroup+Period+Start+LabelCond, data = Sac.sum, keep.names = T))
	 
summary(lmer(SacTarg~LabelCond+ (1|Subj), data = subset(Sac.sum, AgeGroup !="Excl" & Period == "Pre")))
summary(lmer(SacTarg~LabelCond+ (1|Subj), data = subset(Sac.sum, AgeGroup !="Excl" & Period == "Naming" & Start == "Before")))
summary(lmer(SacTarg~LabelCond+ (1|Subj), data = subset(Sac.sum, AgeGroup !="Excl" & Period == "Naming" & Start == "After")))	 	 
	
	 
	 	 
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

barplot(o, beside =T , ylim = c(0,0.25),col = "white",  border = NA, ylab = "Proportion Critical Saccades", names.arg = c("Preview", "Pre-Naming","Post-Naming"))
legend(1.2,0.10, legend = c("Control", "Ambiguous Description", "Unambiguous Description"), bty = "n", col = c("blue","grey","red"), pch = 20)
points(c(1.5,6,10), o[,1], pch = 20, cex = 2, col = "blue")
points(c(2.5,6.8,11), o[,2], pch = 20, cex = 2, col = "grey")
points(c(3.5,7.6,12), o[,3], pch = 20, cex = 2, col = "red")
grid(nx = NA, ny = NULL, col = "gray", lty = "dotted",lwd = par("lwd"), equilogs = TRUE)
abline(v = c(4.5,8.5), col = "grey", lty = "dashed")
arrows(c(1.5,6,10,2.5,6.8,11,3.5,7.6,12), (c(o) + c(se)+0.01), c(1.5,6,10,2.5,6.8,11,3.5,7.6,12), (c(o) - c(se)-0.01), code = 0)