library(plyr)
library(doBy)
library(lme4)
library(ggplot2)


# Initial anlyses -- does fixation duration during 500ms ISI predict informative response  
isi <- read.csv("./data/500ms_Ambig.csv")

isi$foil_duration_R <- (isi$foil_duration_R - mean(isi$foil_duration_R))/sd(isi$foil_duration_R)
isi$target_duration_R <- (isi$target_duration_R - mean(isi$target_duration_R))/sd(isi$target_duration_R)

summary(glmer(informative_foil ~ foil_duration_R + (1|participant_number), data = isi, family = "binomial"))
summary(glmer(informative_foil ~ target_duration_R + (1|participant_number), data = isi, family = "binomial"))

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
# Note that model does not converge with response as a random effect for Subject or Trial_ID, nor with condition as random effect for
# Trial_ID 

resp.graph <- summaryBy(Informative ~ condition + Response + Subject, data = resp, keep.names = T ,na.rm = T)
resp.graph <- summaryBy(Informative ~ condition + Response , data = resp.graph, FUN = c(mean,sd) )
resp.graph$SE <- resp.graph$Informative.sd/sqrt(length(unique(resp$Subject)))
resp.graph$condition <- ordered(resp.graph$condition, levels = c("control","ambiguous"), labels = c("Control","Ambiguous"))
resp.graph$Response <- ordered(resp.graph$Response, levels = c("target","foil"), labels = c("Target","Foil"))

dodge <- position_dodge(width=0.9)
qplot(resp.graph$condition, resp.graph$Informative.mean, fill = resp.graph$Response, ylab = "Proportion Informative Trials", xlab = "Trial Type", ylim = c(0,0.5)) +  geom_errorbar(aes(ymax = resp.graph$Informative.mean + resp.graph$SE, ymin = resp.graph$Informative.mean - resp.graph$SE), width=0.25, position = dodge) + labs(fill = "Response Order") + theme(axis.text.x = element_text(colour = "black", size = 12)) + geom_bar(stat = "identity", position = "dodge")

ggplot(isi, aes(x= foil_duration_R, y= informative_foil)) +
    geom_point() +    # Use hollow circles
    geom_smooth(   # Add linear regression line
                )  + labs(y = "Informativeness of Response", x = "Standardized Gaze Time to Foil") + geom_jitter(width = 0.05, height = 0.02)