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
  sac$type <- factor(sac$type)
	sac$cond <- factor(sac$cond)
  contrasts(sac$cond)[1] <- -1
	contrasts(sac$type)[1] <- -1
	return(sac)
}

ad.sac <- sac.process2("./EyeData/",c("AdSame","Homoph"))
ad.sac$order <- 1:length(ad.sac$SacTarg)
names <- read.delim("./EyeData/AdultWriteNames.txt")
ad.sac <- merge(ad.sac,names, by = c("Subj","trialnum"), all.x = TRUE)
ad.sac<- ad.sac[order(ad.sac$order),]

ad.sac$dupic_sac <- NA
for (i in unique(ad.sac$Subj)){
  for (j in unique(subset(ad.sac, Subj == i)$trialnum)){
    ad.sac[ad.sac$Subj == i & ad.sac$trialnum ==j,]$dupic_sac <- duplicated(subset(ad.sac, Subj == i & trialnum ==j)$SacEndTime,fromLast = TRUE)
  }
}
ad.sac <- ad.sac[ad.sac$dupic_sac == FALSE,]
ad.sac <- ddply(ad.sac, .(Subj,trialnum,Period), transform, CumTarg = cumsum(SacTarg),CumD1 = cumsum(SacDist1),CumD2 = cumsum(SacDist2), SacTime = SacEndTime - min(SacEndTime), SacBin = round((SacEndTime - min(SacEndTime))/100))

ad.sac$Item <- as.factor(sapply(strsplit(as.character(ad.sac$targname),"[.12]"), "[", 1))
# Create LabelCond variable; used for data analysis
ad.sac$LabelCond <- NA
ad.sac[ad.sac$cond == "Control" & ad.sac$Label %in% c(1,0),]$LabelCond <- "Control"
ad.sac[ad.sac$cond == "Ambig" & ad.sac$Label %in% c(1),]$LabelCond <- "Test-Ambig"
ad.sac[ad.sac$cond == "Ambig" & ad.sac$Label %in% c(0),]$LabelCond <- "Test-Unambig"
ad.sac$LabelCond <- as.factor(ad.sac$LabelCond)


ad.sac$Start <- "Preview"
ad.sac[ad.sac$Period == "Naming",]$Start <- ifelse(ad.sac[ad.sac$Period == "Naming",]$SacTime < ad.sac[ad.sac$Period == "Naming",]$StartTime,"Before","After")
ad.sac[ad.sac$Period == "Rew",]$Start <- "After"

summaryBy(SacDist1+SacDist2+SacTarg ~Subj+Item+cond+type+Start, data = ad.sac, keep.names = T) -> Sac.sum
Sac.sum$SacOnTrial = ifelse(Sac.sum$SacTarg > 0,1,0)
na.omit(summaryBy(SacTarg +SacDist2 + SacOnTrial~Start+type+cond, data = Sac.sum, keep.names = T, FUN = c(mean,sd)))

summary(lmer(SacTarg~cond*type + (1+cond|Subj)+(1+cond|Item), data = subset(Sac.sum,  Start == "Preview")))
summary(lmer(SacTarg~cond*type + (1+cond|Subj)+(1+cond|Item), data = subset(Sac.sum,   Start == "Before")))
summary(lmer(SacTarg~cond*type + (1+cond|Subj)+(1+cond|Item), data = subset(Sac.sum,  Start == "After")))        

summary(lmer(SacTarg~cond + (1+cond|Subj)+(1+cond|Item), data = subset(Sac.sum,  Start == "Preview" & type == "Homoph")))
summary(lmer(SacTarg~cond + (1+cond|Subj)+(1+cond|Item), data = subset(Sac.sum,  Start == "Preview" & type == "Same")))


summary(lmer(SacOnTrial~cond*type + (1+cond|Subj)+(1+cond|Item), data = subset(Sac.sum,  Start == "Preview"), family = "binomial"))
summary(lmer(SacOnTrial~cond*type + (1+cond|Subj)+(1+cond|Item), data = subset(Sac.sum,   Start == "Before"), family = "binomial"))
summary(lmer(SacOnTrial~cond*type + (1+cond|Subj)+(1+cond|Item), data = subset(Sac.sum,  Start == "After"), family = "binomial"))        


SacGraph <- function(Sac.graph){
  na.omit(summaryBy(SacTarg.mean~Start+cond, data = Sac.graph, FUN = c(mean,sd))) -> Sac.graph
  Sac.graph$Start <- factor(Sac.graph$Start, levels = c("Preview","Before","After"),labels = c("Preview","Before","After"), ordered = T)
  Sac.graph$SE = Sac.graph$SacTarg.mean.sd/sqrt(length(unique(Sac.sum$Subj)))
  
  Sac.graph$Time = "Pre-Naming"
  Sac.graph[Sac.graph$Start == "After" , ]$Time = "Post-Naming"
  Sac.graph[Sac.graph$Start == "Preview" , ]$Time = "Preview"
  Sac.graph$Time <- factor(Sac.graph$Time, levels = c("Preview","Pre-Naming","Post-Naming"),labels = c("Preview","Pre-Naming","Post-Naming"), ordered = T)
  Sac.graph$cond <- factor(Sac.graph$cond, levels = c("Control","Ambig"),labels = c("Control","Ambig"), ordered = T)
  tapply(Sac.graph$SacTarg.mean.mean, list(Sac.graph$cond,Sac.graph$Time), FUN = mean) -> o
  tapply(Sac.graph$SE, list(Sac.graph$cond,Sac.graph$Time), FUN = mean) -> se
  
  barplot(o, beside =T , ylim = c(0,0.3),col = "white",  border = NA, ylab = "Proportion Critical Saccades", names.arg = c("Preview", "Pre-Naming","Post-Naming"))
  legend(1.2,0.10, legend = c("Control Scene", "Ambiguous Scene"), bty = "n", col = c("blue","red"), pch = 20)
  points(c(1.5,4.3,7.3), o[1,], pch = 20, cex = 2, col = "blue")
  points(c(2.5,5.5,8.5), o[2,], pch = 20, cex = 2, col = "red")
  grid(nx = NA, ny = NULL, col = "gray", lty = "dotted",lwd = par("lwd"), equilogs = TRUE)
  abline(v = c(3.4,6.4), col = "grey", lty = "dashed")
  arrows(c(1.5,2.5,4.3,5.5,7.3,8.5), (c(o) + c(se)+0.01), c(1.5,2.5,4.3,5.5,7.3,8.5), (c(o) - c(se)-0.01), code = 0)
}

# Homophones
SacGraph(na.omit(summaryBy(SacTarg~Start+cond+Subj, data = Sac.sum[ Sac.sum$type == "Homoph",], FUN = c(mean,sd), keep.names = T)))
# Same Cat
SacGraph(na.omit(summaryBy(SacTarg~Start+cond+Subj, data = Sac.sum[ Sac.sum$type == "Same",], FUN = c(mean,sd), keep.names = T)))
