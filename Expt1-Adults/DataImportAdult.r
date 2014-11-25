require(doBy)
require(lme4)
require(boot)

# Import the data and tidy up the column names etc
data_import <- function(path_name,type = "Homoph"){ 
	print(paste("Importing", type))
	data.set = c()
	#Import Homoph data
	for (i in c("Homoph","AdSame")){
		for (j in c("Pre","Targ","Rew")){# Something is off with the reward section ,"Rew")){
		file.list <- list.files(path = path_name,full.names = T,pattern = paste(i,".*","Fix",j, sep = ""))
		print(file.list)
		data.temp = read.delim(file = file.list, header = T)
		data.temp$Period = ifelse(j == "Targ","Naming",j)
		data.set = rbind(data.set,data.temp)
		}
	}
	data.set$Period = ordered(data.set$Period, levels = c("Pre", "Naming","Rew"))	
	data.set$cond = ordered(data.set$cond, levels = c("Control", "Ambig"))	

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

# Plot the data
BarPlotGaze <- function(solo_sum,IV,ylabel){
tapply(solo_sum[,IV],list(solo_sum$cond,solo_sum$type),mean, na.rm = T) -> solo_plot

tapply(solo_sum[,IV], solo_sum[,c("cond","type")], FUN = boot, statistic = boot.mean, R = 10000) -> b.solo
sapply(lapply(b.solo,boot.ci),"[[","normal")[2:3,] -> b.ci

print(solo_plot)
barplot(solo_plot, beside = T, col = "white",border = NA, ylim = c(0, max(b.ci)+0.1),ylab = ylabel, legend = F,xlab = "Ambiguity Type") -> solo_bplot
legend(3.5,0.25, legend = c("Unambig", "Ambig"), bty = "n", fill = c("blue","red"))
grid(nx = NA, ny = NULL, col = "gray", lty = "dotted",
     lwd = par("lwd"), equilogs = TRUE)
points(solo_bplot, solo_plot, pch = 15, cex = 5, col = c("blue","red"))
arrows(solo_bplot,(b.ci[1,]),solo_bplot,(b.ci[2,]),code = 0, length = 0.1, angle = 90)
}

# For bootstrap
boot.mean <- function(x, ind){mean(x[ind],na.rm = T, trim =0)}

# Process saccades
sac.process = function(pathway = "./", Pop = "NA"){
	sac = c()
	for (i in unique(Pop)){
		file.list <- list.files(path = pathway,full.names = T,pattern = paste(i,".*","SacPre", sep = ""))
		sac.temp <- read.delim(file.list, header = T)
		sac <- rbind(sac, sac.temp)	
	}
	names(sac)[names(sac) == "RECORDING_SESSION_LABEL"] <- "Subj"
	names(sac)[names(sac) == "CURRENT_SAC_START_INTEREST_AREA_LABEL"] <- "SacStart"
	names(sac)[names(sac) == "CURRENT_SAC_END_INTEREST_AREA_LABEL"] <- "SacEnd"
	#sac$Subj <- as.factor(paste(sac$Subj,".edf", sep = ""))
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
	sac <- summaryBy(SacTarg + SacSwitch + Sac + SacDist1 + SacDist2 ~ Subj + trialnum + cond, data = sac, keep.names = T)
	return(sac)
}



# Import and plot the data,
adult.ref <- data_import("./EyeData/")
adult.ref.t <- adult.ref[adult.ref$type != "Filler",]
adult.ref.t$type <- adult.ref.t$type[drop= TRUE]
adult.ref.t$cond <- adult.ref.t$cond[drop= TRUE]

adult.sac = sac.process("./EyeData/",Pop = c("AdSame","Homoph"))
adult.ref.t <- merge(adult.ref.t, adult.sac, by = c("Subj", "trialnum","cond"), all.x = TRUE)

# If statements tidy up the importing of sound coded trials -- writes a file where sound trials can be imported into.
if (file.exists("./EyeData/AdultWriteNames.txt") == TRUE){
	names = read.delim("./EyeData/AdultWriteNames.txt")
	if(length(names$FixCount > 0)){names$FixCount <- NULL}
	adult.ref.t <- merge(adult.ref.t,names, by = c("Subj","trialnum"), all.x = TRUE)
	prntout <- data.frame(aggregate(Label~Subj+trialnum, data = adult.ref.t, mean, na.action = na.pass))
	write.table(prntout[order(prntout$Subj),], file = "AdultWriteNames.txt", sep = "\t", row.names = F)
	}else{
			prntout <- data.frame(cbind(aggregate(FixCount~Subj+trialnum, data = adult.ref.t, sum),Label = NA))
			write.table(prntout[order(prntout$Subj),], file = "AdultWriteNames.txt", sep = "\t", row.names = F)
			}




adult.ref.t$Item <- as.factor(sapply(strsplit(as.character(adult.ref.t$targname),"[.12]"), "[", 1))

summary(lmer(SacTarg~cond*type + (1|Subj)+(1|Item), data = subset(adult.ref.t, subset = Period == "Pre" & Picture == "Targ" & TRIAL_DWELL_TIME > 1500)))


adult.ref.t.s <- adult.ref.t[adult.ref.t$Period == "Pre" & adult.ref.t$Picture == "Targ",]
adult.ref.t.s$Label2 <- 1 - adult.ref.t.s$Label
BarPlotGaze(adult.ref.t.s, "SacTarg", "Proportion critical saccades")
BarPlotGaze(adult.ref.t.s, "Label","Proportion ambiguous labels")
