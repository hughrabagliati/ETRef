library(plyr)
library(doBy)
library(lme4)


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