library(lme4)
library(doBy)
TwoNames <- read.csv("TwoNames.csv", header = T)
summaryBy(Label+Age~Condition+Order, data = TwoNames)

summary(glmer(Label~Order + (1|Participant.Name) + (1|Target.Name), data = TwoNames[TwoNames$Condition != "Control",], family = "binomial"))

# Now run the analyses based on the hypothesis that Unambig (second name) should have less modification than
# Ambiguous (first name), which should have less modification than Ambiguous (second name)

TwoNames$Ordering <- 1
TwoNames[TwoNames$Condition == "Control" & TwoNames$Order == 1,]$Ordering <- NA
TwoNames[TwoNames$Condition == "Ambig" & TwoNames$Order == 1,]$Ordering <- 2
TwoNames[TwoNames$Condition == "Ambig" & TwoNames$Order == 2,]$Ordering <- 3

# Function to write Forward Difference Coding, comparing each level against the previous.
forward_diff <- function(k){ # k = number of levels of the factor
j <- k-1
fd <- matrix(NA, nrow =k, ncol = j , dimnames = list(1:k,1:j))

r = 1
for (h in 1:j){
	for (i in 1:k){
		if (i <= r){
			fd[i,h] <- (k-r)/k
				}else{
			fd[i,h] <- -r/k
			}
		}
	r <- r+1
}
return(fd)
}

TwoNames$Ordering <- as.factor(TwoNames$Ordering)
contrasts(TwoNames$Ordering) <- forward_diff(length(levels(TwoNames$Ordering)))

summary(glmer(Label~as.factor(Ordering) + (1|Participant.Name) + (1|Target.Name), data = TwoNames, family = "binomial"))
