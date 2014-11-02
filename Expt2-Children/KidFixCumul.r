		require(plyr)
	require(longitudinalData)
	fix.process2 = function(pathway = "./", Pop = "NA"){
		sac = c()
		for (i in unique(Pop)){
			for (j in c("Pre","Targ")){# Something is off with the reward section ,"Rew")){
			file.list <- list.files(path = pathway,full.names = T,pattern = paste(i,".*","FiRep",j, sep = ""))
	 		print(file.list)
	 		sac.temp <- read.delim(file.list, header = T)
			sac.temp$Period = ifelse(j == "Targ","Naming",j)
			sac <- rbind(sac, sac.temp)	
			}
		}
		
		names(sac)[names(sac) == "RECORDING_SESSION_LABEL"] <- "Subj"
		names(sac)[names(sac) == "CURRENT_FIX_NEAREST_INTEREST_AREA_LABEL"] <- "Picture"
		names(sac)[names(sac) == "CURRENT_FIX_DURATION"] <- "Duration"
		names(sac)[names(sac) == "CURRENT_FIX_END"] <- "FixEndTime"
		
		levels(sac$Picture) <- list(Foil="Pre_D1 ", Dist="Pre_D2 ", Targ="Pre_Targ ")
	
		
		sac$Subj <- as.factor(paste(sac$Subj,".edf", sep = ""))
		sac$FixTarg <- 0
		sac[sac$Picture == "Targ",]$FixTarg <- 1
		sac$FixFoil <- 0
		sac[sac$Picture == "Foil",]$FixFoil <- 1
		sac$FixDist <- 0
		sac[sac$Picture == "Dist",]$FixDist <- 1
		sac <- sac[sac$type != "Filler",]
		return(sac)
	}
	
	kid.fix <- fix.process2("./","Kids")
	kid.fix$order <- 1:length(kid.fix$FixTarg)
	Groups = read.delim("SubjNames.txt", header = T)
	merge(kid.fix,Groups, all = T) -> kid.fix
	kid.fix[is.na(kid.fix$AgeGroup),]$AgeGroup <- "Old"
	names <- read.delim("WriteNames-Times.txt")
	kid.fix <- merge(kid.fix,names, by = c("Subj","trialnum"), all.x = TRUE)
	kid.fix<- kid.fix[order(kid.fix$order),]
	#names(kid.fix)[names(kid.fix) == "StartTime..ms."] <- "StartTime"
	
	# The next line adjusts for fixations that straddle the boundary between preview and naming
	kid.fix[kid.fix$Period == "Pre" & kid.fix$FixEndTime > 4250,]$Duration <- kid.fix[kid.fix$Period == "Pre" & kid.fix$FixEndTime > 4250,]$Duration - (kid.fix[kid.fix$Period == "Pre" & kid.fix$FixEndTime > 4250,]$FixEndTime - 4250)
	kid.fix[kid.fix$Period == "Pre" & kid.fix$FixEndTime > 4250,]$FixEndTime <- 4250
	
	kid.fix$FixTime = kid.fix$FixEndTime
	# Can also make fixtime relative to the onset of the namen
	kid.fix[kid.fix$Period == "Naming",]$FixTime <- kid.fix[kid.fix$Period == "Naming",]$FixTime - kid.fix[kid.fix$Period == "Naming",]$StartTime
	kid.fix <- ddply(kid.fix, .(Subj,trialnum,Period), transform, CumTarg = cumsum(FixTarg*Duration),CumFoil = cumsum(FixFoil*Duration),CumDist = cumsum(FixDist*Duration), FixBin = round((FixTime)/100))
	
	
	# Create LabelCond variable; used for data analysis
	kid.fix$LabelCond <- NA
	kid.fix[kid.fix$cond == "Control" & kid.fix$Label %in% c(1,0),]$LabelCond <- "Control"
	kid.fix[kid.fix$cond == "Ambig" & kid.fix$Label %in% c(1),]$LabelCond <- "Test-Ambig"
	kid.fix[kid.fix$cond == "Ambig" & kid.fix$Label %in% c(0),]$LabelCond <- "Test-Unambig"
	kid.fix$LabelCond <- as.factor(kid.fix$LabelCond)
	
	ddply(kid.fix, .(Subj,trialnum,cond,Period,targname,LabelCond,Label, StartTime,AgeGroup,Lang), summarize, FixBin = c(0:81)) -> kid.fix.s
	ddply(kid.fix[kid.fix$FixBin <82,], .(Subj,trialnum,cond,Period,targname,LabelCond,Label,FixBin), summarize, CumTarg = max(CumTarg),CumFoil = max(CumFoil),CumDist = max(CumDist)) -> kid.fix.sum
	
	kid.fix.bin <- merge(kid.fix.s, kid.fix.sum, all = TRUE)
	kid.fix.bin[kid.fix.bin$FixBin == 0,][,c("CumTarg","CumFoil","CumDist")] <- 0
	kid.fix.bin$CumTarg <- t(imputation(matrix(kid.fix.bin$CumTarg, nrow = 1),method = "locf"))
	kid.fix.bin$CumFoil <- t(imputation(matrix(kid.fix.bin$CumFoil, nrow = 1),method = "locf"))
	kid.fix.bin$CumDist <- t(imputation(matrix(kid.fix.bin$CumDist, nrow = 1),method = "locf"))
	kid.fix.bin$Item <- as.factor(sapply(strsplit(as.character(kid.fix.bin$targname),"[.12]"), "[", 1))
	
	kid.fix.bin$SubjTrial <- paste(kid.fix.bin$Subj,kid.fix.bin$Item, sep = "") 
	TrialTime <- summaryBy(TRIAL_DWELL_TIME~SubjTrial, data = subset(kid.ref.t, subset = Period == "Pre" & Picture == "Targ"), keep.names = T, na.rm = T)
	kid.fix.bin <- kid.fix.bin[kid.fix.bin$SubjTrial %in% TrialTime[TrialTime$TRIAL_DWELL_TIME > 1500,]$SubjTrial,]
	
	kid.fix.bin$StartTime <- kid.fix.bin$StartTime/100
	kid.fix.bin$StartTime <- round(kid.fix.bin$StartTime)
	kid.fix.bin$FixBin2 <- kid.fix.bin$FixBin - kid.fix.bin$StartTime
	
	plot.fix.cum <- function(sac.bin, Peri){
	summaryBy(CumTarg+CumFoil+CumDist~LabelCond+FixBin, data = subset(sac.bin, Period == Peri), keep.names = T) -> k 
	plot(k[k$LabelCond == "Control",]$FixBin, k[k$LabelCond == "Control",]$CumTarg, ylim = c(0,2400), pch = 21, col = "darkgray",bg = "darkgrey")
	points(k[k$LabelCond == "Test-Ambig",]$FixBin, k[k$LabelCond == "Test-Ambig",]$CumTarg, ylim = c(0,2), pch = 21, col = "blue",bg = "blue")
	points(k[k$LabelCond == "Test-Unambig",]$FixBin, k[k$LabelCond == "Test-Unambig",]$CumTarg, ylim = c(0,2), pch = 21, col = "red",bg = "red")
	}
	plot.fix.cum(kid.fix.bin,"Naming")
	abline(v=(median(kid.fix$StartTime, na.rm = T)/100))
	plot.fix.cum(kid.fix.bin,"Pre")
	summary(lmer(CumTarg~LabelCond*FixBin + (1+LabelCond|Subj)+(1+LabelCond|Item), data = kid.fix.bin[kid.fix.bin$Period == "Naming",]))
	 summary(lmer(CumTarg~LabelCond*FixBin + (1+LabelCond|Subj)+(1+LabelCond|Item), data = subset(kid.fix.bin,Period == "Pre" &FixBin <43 & AgeGroup == "Young" & Lang == "Mon")))
	 
	 kid.fix$Start <- "Before"
 	 kid.fix[kid.fix$Period == "Naming",]$Start <- ifelse(kid.fix[kid.fix$Period == "Naming",]$FixTime < kid.fix[kid.fix$Period == "Naming",]$StartTime,"Before","After")

	 summaryBy(FixFoil+FixTarg+FixDist ~AgeGroup+Subj+trialnum+cond+Period+Start, data = kid.fix, keep.names = T) -> fix.sum
	 fix.sum$FixDist <- ifelse(fix.sum$FixDist > 0,1,0)
	 fix.sum$FixFoil <- ifelse(fix.sum$FixFoil > 0,1,0)
	 fix.sum$FixTarg <- ifelse(fix.sum$FixTarg > 0,1,0)
	 na.omit(summaryBy(FixFoil +FixDist+FixTarg ~AgeGroup+Period+Start+cond, data = fix.sum, keep.names = T))
	 
	 