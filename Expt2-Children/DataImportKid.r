require(doBy)
require(lme4)
require(boot)


# Import the data and tidy up the column names etc
data_import <- function(path_name,type = "Kids"){ 
	print(paste("Importing", type))
	data.set = c()
	#Import Homoph data
	for (i in c("Kids")){
		for (j in c("Pre","Targ")){# Something is off with the reward section ,"Rew")){
		file.list <- list.files(path = path_name,full.names = T,pattern = paste(i,".*","Fix",j, sep = ""))
		print(file.list)
		data.temp = read.delim(file = file.list, header = T)
		data.temp$Period = ifelse(j == "Targ","Naming",j)
		data.set = rbind(data.set,data.temp)
		}
	}
	data.set$Period = ordered(data.set$Period, levels = c("Pre", "Naming","Rew"))	
	data.set$cond = ordered(data.set$cond, levels = c("Control", "Ambig"))	
	
	data.set$IA_FIXATION_.[data.set$IA_FIXATION_. == "."] <- NA
	data.set$IA_DWELL_TIME_.[data.set$IA_DWELL_TIME_. == "."] <- NA

	data.set$IA_FIXATION_. <- as.numeric(as.character(data.set$IA_FIXATION_.))
	data.set$IA_DWELL_TIME_. <- as.numeric(as.character(data.set$IA_DWELL_TIME_.))
	
	names(data.set)[names(data.set) == "IA_LABEL"] <- "Picture"
	names(data.set)[names(data.set) == "IA_DWELL_TIME"] <- "DwellTime"
	names(data.set)[names(data.set) == "IA_FIXATION_COUNT"] <- "FixCount"
	names(data.set)[names(data.set) == "DATA_FILE"] <- "Subj"
	names(data.set)[names(data.set) == "IA_FIXATION_."] <- "PropFix"
	names(data.set)[names(data.set) == "IA_DWELL_TIME_."] <- "PropDwell"
	data.set$PropFix = as.numeric(data.set$PropFix)
	data.set$PropDwell = as.numeric(data.set$PropDwell)
	levels(data.set$Picture) <- list(Foil="Pre_D1 ", Dist="Pre_D2 ", Targ="Pre_Targ ")
	return(data.set)
	}

#Plot the data
BarPlotGaze <- function(solo_sum,IV, ylabel){
tapply(solo_sum[,IV],list(solo_sum$LabelCond),mean,na.rm = T) -> solo_plot

tapply(solo_sum[,IV], solo_sum[,c("LabelCond")], FUN = boot, statistic = boot.mean, R = 10000) -> b.solo
sapply(lapply(b.solo,boot.ci),"[[","normal")[2:3,] -> b.ci

print(solo_plot)
barplot(solo_plot, beside = T, col = "white",border = NA, ylim = c(0, max(b.ci)+0.05),ylab = ylabel, legend = F,args.legend = list(x = 1.5, bty = "n"),xlab = "Ambiguity Type") -> solo_bplot
grid(nx = NA, ny = NULL, col = "gray", lty = "dotted",
     lwd = par("lwd"), equilogs = TRUE)
points(solo_bplot, solo_plot, pch = 15, cex = 5, col = c("red","blue","blue"))
arrows(solo_bplot,(b.ci[1,]),solo_bplot,(b.ci[2,]),code = 0, length = 0.1, angle = 90)
}

# For bootstrap
boot.mean <- function(x, ind){mean(x[ind],na.rm = T, trim =0)}

# Process saccades
sac.process = function(pathway = "./", Pop = "NA"){
	sac = c()
	for (i in unique(Pop)){
		for (j in c("Pre","Targ")){# Something is off with the reward section ,"Rew")){
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
	sac$Subj <- as.factor(paste(sac$Subj,".edf", sep = ""))
	sac$Sac <- 1
	sac$SacSwitch <- ifelse(sac$SacStart == sac$SacEnd,0,1)
	sac$SacTarg <- 0
	sac[sac$SacStart == "Pre_Targ " & sac$SacEnd == "Pre_D1 ",]$SacTarg <- 1
	sac[sac$SacStart == "Pre_D1 " & sac$SacEnd == "Pre_Targ ",]$SacTarg <- 1
	sac$SacDist1 <- 0
	sac[sac$SacStart == "Pre_D1 " & sac$SacEnd == "Pre_D2 ",]$SacDist1 <- 1
	sac[sac$SacStart == "Pre_D2 " & sac$SacEnd == "Pre_D1 ",]$SacDist1 <- 1
	sac$SacDist2 <- 0
	sac[sac$SacStart == "Pre_Targ " & sac$SacEnd == "Pre_D2 ",]$SacDist2 <- 1
	sac[sac$SacStart == "Pre_D2 " & sac$SacEnd == "Pre_Targ ",]$SacDist2 <- 1
	sac <- summaryBy(SacTarg + SacSwitch + Sac + SacDist1 + SacDist2 ~ Subj + trialnum + cond+Period, data = sac, keep.names = T)
	return(sac)
}



# Import and plot the data,
kid.ref <- data_import("./EyeData/")
kid.ref.t <- kid.ref[kid.ref$type != "Filler",]
kid.ref.t$type <- kid.ref.t$type[drop= TRUE]
kid.ref.t$cond <- kid.ref.t$cond[drop= TRUE]
kid.ref.t <- kid.ref.t[!is.na(kid.ref.t$cond),]

Groups = read.delim("./EyeData/SubjNames.txt", header = T)
merge(kid.ref.t,Groups, all = T) -> kid.ref.t
kid.ref.t[is.na(kid.ref.t$AgeGroup),]$AgeGroup <- "Old"
#kid.ref.t$Lang <- "Mon"

kid.sac = sac.process("./EyeData/","Kids")
kid.ref.t <- merge(kid.ref.t, kid.sac, by = c("Subj", "trialnum","cond","Period"), all.x = TRUE)


# If statements tidy up the importing of sound coded trials -- writes a file where sound trials can be imported into.
if (file.exists("./EyeData/WriteNames-Times.txt") == TRUE){
	names = read.delim(./EyeData/"WriteNames-Times.txt")
	if(length(names$FixCount > 0)){names$FixCount <- NULL}
	kid.ref.t <- merge(kid.ref.t,names, by = c("Subj","trialnum"), all.x = TRUE)
	prntout <- data.frame(summaryBy(Label+StartTime~Subj+trialnum, data = kid.ref.t, FUN = mean, na.rm=T, keep.names = T))
	write.table(prntout[order(prntout$Subj),], file = "./EyeData/WriteNames-Times.txt", sep = "\t", row.names = F)
	}else{
			prntout <- data.frame(cbind(aggregate(FixCount~Subj+trialnum, data = kid.ref.t, sum),Label = NA,StartTime = NA))
			write.table(prntout[order(prntout$Subj),], file = "./EyeData/WriteNames-Times.txt", sep = "\t", row.names = F)
			}

# Create a measure of difference in labeling
#kid.ref.t <- merge(kid.ref.t,data.frame(Subj = unique(kid.ref.t$Subj), Sub = (summaryBy(Label~Subj, data = subset(kid.ref.t, cond == "Control"),na.rm = T)$Label.mean  - summaryBy(Label~Subj, data = subset(kid.ref.t, cond == "Ambig"),na.rm = T)$Label.mean)))

# Create LabelCond variable; used for data analysis
kid.ref.t$LabelCond <- NA
kid.ref.t[kid.ref.t$cond == "Control" & kid.ref.t$Label %in% c(1,0),]$LabelCond <- "Control"
kid.ref.t[kid.ref.t$cond == "Ambig" & kid.ref.t$Label %in% c(1),]$LabelCond <- "Test-Ambig"
kid.ref.t[kid.ref.t$cond == "Ambig" & kid.ref.t$Label %in% c(0),]$LabelCond <- "Test-Unambig"
kid.ref.t$LabelCond <- as.factor(kid.ref.t$LabelCond)

#Create Item labels
kid.ref.t$Item <- as.factor(sapply(strsplit(as.character(kid.ref.t$targname),"[.12]"), "[", 1))

kid.ref.t$SubjTrial <- paste(kid.ref.t$Subj,kid.ref.t$Item, sep = "") 
TrialTime <- summaryBy(TRIAL_DWELL_TIME~SubjTrial, data = subset(kid.ref.t, subset = Period == "Pre" & Picture == "Targ"), keep.names = T, na.rm = T)
kid.ref.t <- kid.ref.t[kid.ref.t$SubjTrial %in% TrialTime[TrialTime$TRIAL_DWELL_TIME > 1500,]$SubjTrial,]


summary(lmer(PropDwell~LabelCond + (1+LabelCond|Subj)+(1|Item), data = subset(kid.ref.t, subset = Period == "Pre" & Picture == "Targ" & AgeGroup == "Young" & Lang == "Mon")))
summary(lmer(PropDwell~LabelCond + (1+LabelCond|Subj)+(1|Item), data = subset(kid.ref.t, subset = Period == "Pre" & Picture == "Targ" & AgeGroup == "Old" & Lang == "Mon")))
summary(lmer(PropDwell~LabelCond + (1+LabelCond|Subj)+(1|Item), data = subset(kid.ref.t, subset = Period == "Pre" & Picture == "Dist" & AgeGroup == "Young" & Lang == "Mon")))
summary(lmer(PropDwell~LabelCond + (1+LabelCond|Subj)+(1|Item), data = subset(kid.ref.t, subset = Period == "Pre" & Picture == "Dist" & AgeGroup == "Old" & Lang == "Mon")))

summary(glmer(Label~PropDwell + (1+PropDwell|Subj)+(1|Item), data = subset(kid.ref.t, subset = Period == "Pre" & Picture == "Targ" & cond == "Ambig" & AgeGroup == "Young" & Lang == "Mon") , family = "binomial"))
summary(glmer(Label~PropDwell + (1+PropDwell|Subj)+(1|Item), data = subset(kid.ref.t, subset = Period == "Pre" & Picture == "Targ" & cond == "Ambig" & AgeGroup == "Old" & Lang == "Mon") , family = "binomial"))
summary(glmer(Label~PropDwell + (1+PropDwell|Subj)+(1|Item), data = subset(kid.ref.t, subset = Period == "Pre" & Picture == "Dist" & cond == "Ambig" & AgeGroup == "Young" & Lang == "Mon") , family = "binomial"))
summary(glmer(Label~PropDwell + (1+PropDwell|Subj)+(1|Item), data = subset(kid.ref.t, subset = Period == "Pre" & Picture == "Dist" & cond == "Ambig" & AgeGroup == "Old" & Lang == "Mon") , family = "binomial"))


kid.ref.t.s <- kid.ref.t[kid.ref.t$Period == "Pre" & kid.ref.t$Picture == "Targ",]
kid.ref.t.s$Label2 <- 1- kid.ref.t.s$Label

BarPlotGaze(kid.ref.t.s, "PropDwell","PropDwell")
