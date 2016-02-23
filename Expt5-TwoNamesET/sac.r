library(plyr)
library(doBy)
library(lme4)


# Initial anlyses -- how does informativity of response vary across response order and ambiguity condition?
resp <- read.csv("./data/responses_ambig_and_cont.csv")
resp <- reshape(resp, 
  varying = c("informative_target", "informative_foil"), 
  v.names = "Informative",
  timevar = "Response", 
  times = c("target", "foil"),
  direction = "long")

resp$Response <- as.factor(resp$Response)
contrasts(resp$Response)[1] <- -1
summaryBy(Informative ~ Response  + condition, data = resp, na.rm = T)
summary(glmer(Informative ~ condition * as.factor(Response) + (1+condition|Subject) + (1|Trial_ID), data = resp, family = "binomial"))


# How do saccades during Preview and 500ms ISI predict informative responses?  
# Note that we are now using an SMI system for data collection, and the data output is 
# v different. 
sac <- read.csv("./data/transition_data.csv")
sac$Targ <- ifelse(sac$informative_target == "True", 1, ifelse(sac$informative_target == "False", 0, NA))
sac$Foil <- ifelse(sac$informative_foil == "True", 1, ifelse(sac$informative_foil == "False", 0, NA))

sac$TotalSaccades <- sac$target_to_distractor + sac$target_to_foil + sac$distractor_to_target + sac$distractor_to_foil + sac$foil_to_target + sac$foil_to_distractor
sac$CriticalSaccades <- sac$target_to_foil + sac$foil_to_target
sac$PropSac <- sac$CriticalSaccades/sac$TotalSaccades
sac$PropSac[is.na(sac$PropSac)] <- 0


# How do fixations during 500ms ISI predict informative responses?  
# Note that we are now using an SMI system for data collection, and the data output is 
# v different. 
fix <- read.csv("./data/all_fixation_durations.csv")
fix$Targ <- ifelse(fix$informative_target == "TRUE", 1, ifelse(fix$informative_target == "FALSE", 0, NA))
fix$Foil <- ifelse(fix$informative_foil == "TRUE", 1, ifelse(fix$informative_foil == "FALSE", 0, NA))

fix$foil_duration_R <- (fix$foil_sum_duration_Right - mean(fix$foil_sum_duration_Right))/sd(fix$foil_sum_duration_Right)
fix$target_duration_R <- (fix$target_sum_duration_Right - mean(fix$target_sum_duration_Right))/sd(fix$target_sum_duration_Right)

fix$TargetLabelCond <- "Control"
fix[fix$condition %in% "ambig" & fix$Targ %in% NA,]$TargetLabelCond <- NA
fix[fix$condition %in% "ambig" & fix$Targ %in% 1,]$TargetLabelCond <- "Ambig - Informative"
fix[fix$condition %in% "ambig" & fix$Targ %in% 0,]$TargetLabelCond <- "Ambig - Uninformative"
fix$TargetLabelCond <- ordered(fix$TargetLabelCond, levels = c("Control", "Ambig - Uninformative", "Ambig - Informative"))

fix$FoilLabelCond <- "Control"
fix[fix$condition %in% "ambig" & fix$Foil %in% NA,]$FoilLabelCond <- NA
fix[fix$condition %in% "ambig" & fix$Foil %in% 1,]$FoilLabelCond <- "Ambig - Informative"
fix[fix$condition %in% "ambig" & fix$Foil %in% 0,]$FoilLabelCond <- "Ambig - Uninformative"
fix$FoilLabelCond <- ordered(fix$FoilLabelCond, levels = c("Control", "Ambig - Uninformative", "Ambig - Informative"))


summaryBy(target_duration_R + foil_duration_R ~ phase + condition + TargetLabelCond, data = fix)
summaryBy(target_duration_R + foil_duration_R ~ phase + condition + FoilLabelCond, data = fix)