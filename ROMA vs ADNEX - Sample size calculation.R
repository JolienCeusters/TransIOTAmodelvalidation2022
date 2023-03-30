#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%# 
#  ROMA vs ADNEX: Sample size  #
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Load packages
library(sn)
library(fGarch)
library(boot)
library(ggplot2)
library(readxl)
library(caret)

#### 1. O/E ####

# target confidence interval width of 1 for O/E; assuming O/E is 1 then this corresponds to SE(lnOE) of 0.245
selnoe = c(0.026, 0.038, 0.051, 0.077, 0.100)
outcome_prop = 0.49
sampsize_OE = (1 - outcome_prop)/(outcome_prop * selnoe^2)
sampsize_OE

#### 2. Calibration slope ####
set.seed(1212)

#### 2.1 ADNEX ####
load("C:/Users/u0123496/OneDrive - KU Leuven/IOTA/IOTA5/IOTA5 - Datasets.RData")

## Previous study
Iota5ValNI$LP <- logit(Iota5ValNI$pmalw)
ggplot(Iota5ValNI, aes(x = LP, group = binaryCDcorrect)) +
  geom_density()
summary(Iota5ValNI$LP[Iota5ValNI$binaryCDcorrect == 0])
sd(Iota5ValNI$LP[Iota5ValNI$binaryCDcorrect == 0])
summary(Iota5ValNI$LP[Iota5ValNI$binaryCDcorrect == 1])
sd(Iota5ValNI$LP[Iota5ValNI$binaryCDcorrect == 1])

set.seed(1212)
LPben <- rsnorm(n = 0.51 * 5000000, mean = -3.5, sd = 1, xi = 1.2)
summary(LPben)
plot(density(Iota5ValNI$LP[Iota5ValNI$binaryCDcorrect == 0]), lwd = 2, col = "red", main = "LP for Benign")
lines(density(DataCC$LP[DataCC$OutcomeBin == 0]), lwd = 2, col = "black")
lines(density(DataImp$LP[DataImp$OutcomeBin == 0]), lwd = 2, col = "grey")
lines(density(LPben), col = "blue", lwd = 2)
legend("topright", legend = c("Observed", "Estimated"),
       col=c("red", "blue"), lwd = 2)
# 550 x 550

set.seed(1212)
LPmal <- rsnorm(n = 0.49 * 5000000, mean = 1, sd = 1.9, xi = -1)
summary(LPmal)
plot(density(Iota5ValNI$LP[Iota5ValNI$binaryCDcorrect == 1]), lwd = 2, col = "red", main = "LP for Malignant")
lines(density(DataCC$LP[DataCC$OutcomeBin == 1]), lwd = 2, col = "black")
lines(density(DataImp$LP[DataImp$OutcomeBin == 1]), lwd = 2, col = "grey")
lines(density(LPmal), col = "blue", lwd = 2)
legend("topright", legend = c("Observed", "Estimated"),
       col=c("red", "blue"), lwd = 2)
# 550 x 550

## Current study
DataCC$LP <- logit(DataCC$pmalw)
DataImp$LP <- logit(DataImp$pmalw)

set.seed(1212)
LPben <- rsnorm(n = 0.51 * 5000000, mean = -2.75, sd = 1.7, xi = 1.2) # Mine
summary(LPben)
plot(density(DataCC$LP[DataCC$OutcomeBin == 0]), lwd = 2, col = "red", main = "LP for Benign")
lines(density(LPben), col = "blue", lwd = 2)
legend("topright", legend = c("Observed", "Estimated"),
       col=c("red", "blue"), lwd = 2)
# 550 x 550

set.seed(1212)
LPmal <- rsnorm(n = 0.49 * 5000000, mean = 1.4, sd = 2.5, xi = 1.2)
summary(LPmal)
plot(density(DataCC$LP[DataCC$OutcomeBin == 1]), lwd = 2, col = "red", main = "LP for Malignant")
lines(density(LPmal), col = "blue", lwd = 2)
legend("topright", legend = c("Observed", "Estimated"),
       col=c("red", "blue"), lwd = 2)
# 550 x 550

LP <- c(LPben, LPmal)

# Input assumed parameters of calibration model: 'weak' level calibration is obtained, i.e. intercept of 0 and slope of 1
beta0 = 0
beta1 = 1
# calculate elements of I matrix
Borenstein_00 = exp(beta0 + (beta1 * LP))/((1 + exp(beta0 + (beta1 * LP)))^2)
Borenstein_01 = LP * exp(beta0 + (beta1 * LP))/((1 + exp(beta0 + (beta1 * LP)))^2)
Borenstein_11 = LP^2 * exp(beta0 + (beta1 * LP))/((1 + exp(beta0 + (beta1 * LP)))^2)

I_00 = mean(Borenstein_00)
I_01 = mean(Borenstein_01)
I_11 = mean(Borenstein_11)
# input desired SE: target a CI width of 0.2, which corresponds to a SE of 0.051
seslope = c(0.026, 0.038, 0.051, 0.064, 0.077, 0.089, 0.100)
# alternatively: selnoe = CIwidth/(2*1.96)
outcome_prop = 0.49
# calculate sample size
sampsize_slope = (I_00/(seslope^2 * ((I_00 * I_11) - (I_01^2))))
sampsize_slope

#### 2.2 ROMA ####
ROMA <- read_excel("C:/Users/u0123496/OneDrive - KU Leuven/IOTA/ROMA vs ADNEX/Sample size calculation/roma risks Kaijser et al 2013.xlsx")

## Previous study
ROMA$LP <- logit(ROMA$ROMAfuji_risk)
ggplot(ROMA, aes(x = LP, group = outcome1)) +
  geom_density()
summary(ROMA$LP[ROMA$outcome1 == 0])
sd(ROMA$LP[ROMA$outcome1 == 0])
summary(ROMA$LP[ROMA$outcome1 == 1])
sd(ROMA$LP[ROMA$outcome1 == 1])

set.seed(1212)
LPben <- rsnorm(n = 0.51 * 5000000, mean = -2.5, sd = 0.9, xi = 1.3)
summary(LPben)
plot(density(ROMA$LP[ROMA$outcome1 == 0]), lwd = 2, col = "red", main = "LP for Benign")
lines(density(LPben), col = "blue", lwd = 2)
legend("topright", legend = c("Observed", "Estimated"),
       col=c("red", "blue"), lwd = 2)
# 550 x 550

set.seed(1212)
LPmal <- rsnorm(n = 0.49 * 5000000, mean = 1.5, sd = 3.3, xi = 1.55)
summary(LPmal)
plot(density(ROMA$LP[ROMA$outcome1 == 1]), lwd = 2, col = "red", main = "LP for Malignant")
lines(density(LPmal), col = "blue", lwd = 2)
legend("topright", legend = c("Observed", "Estimated"),
       col=c("red", "blue"), lwd = 2)
# 550 x 550

## Current study
DataCC$LP <- logit(DataCC$ROMA)
DataImp$LP <- logit(DataImp$ROMA)

set.seed(1212)
LPben <- rsnorm(n = 0.51 * 5000000, mean = -2.2, sd = 0.7, xi = 1.1) # Mine
summary(LPben)
plot(density(DataCC$LP[DataCC$OutcomeBin == 0]), lwd = 2, col = "red", main = "LP for Benign")
lines(density(LPben), col = "blue", lwd = 2)
legend("topright", legend = c("Observed", "Estimated"),
       col=c("red", "blue"), lwd = 2)
# 550 x 550

set.seed(1212)
LPmal <- rsnorm(n = 0.49 * 5000000, mean = 0.3, sd = 2.6, xi = 1.6)
summary(LPmal)
plot(density(DataCC$LP[DataCC$OutcomeBin == 1]), lwd = 2, col = "red", main = "LP for Malignant")
lines(density(LPmal), col = "blue", lwd = 2)
legend("topright", legend = c("Observed", "Estimated"),
       col=c("red", "blue"), lwd = 2)
# 550 x 550

LP <- c(LPben, LPmal)

# Input assumed parameters of calibration model: 'weak' level calibration is obtained, i.e. intercept of 0 and slope of 1
beta0 = 0
beta1 = 1
# calculate elements of I matrix
Borenstein_00 = exp(beta0 + (beta1 * LP))/((1 + exp(beta0 + (beta1 * LP)))^2)
Borenstein_01 = LP * exp(beta0 + (beta1 * LP))/((1 + exp(beta0 + (beta1 * LP)))^2)
Borenstein_11 = LP^2 * exp(beta0 + (beta1 * LP))/((1 + exp(beta0 + (beta1 * LP)))^2)

I_00 = mean(Borenstein_00)
I_01 = mean(Borenstein_01)
I_11 = mean(Borenstein_11)
# input desired SE: target a CI width of 0.2, which corresponds to a SE of 0.051
seslope = c(0.026, 0.038, 0.051, 0.064, 0.077, 0.089, 0.100)
# alternatively: selnoe = CIwidth/(2*1.96)
outcome_prop = 0.49
# calculate sample size
sampsize_slope = (I_00/(seslope^2 * ((I_00 * I_11) - (I_01^2))))
sampsize_slope

#### 3. C-statistic ####

outcome_prop = 0.49

#### 3.1 ROMA ####
Cstat = 0.85 # Werk met laagste AUC-waarde: sample size voor AUC is meestal groter bij de veronderstelling van een lage AUC
# No closed form solution -> an iterative or deductive approach
# Calculate the SE for a range of N
N <- seq(1, 1000000, 1)
seCstatsq = Cstat * (1 - Cstat) * (1 + (((N/2) - 1) * ((1 - Cstat)/(2 - Cstat))) + ((((N/2) - 1) * Cstat)/(1 + Cstat)))/(N^2 * outcome_prop * (1 - outcome_prop))
seCstat = sqrt(seCstatsq)
CIwidth = 2*1.96*seCstat 
# Identify the sample sizes that give a CIwidth no wider than the desired value
Result <- data.frame(cbind(N, CIwidth))
min(Result$N[Result$CIwidth <= 0.050])
min(Result$N[Result$CIwidth <= 0.060])
min(Result$N[Result$CIwidth <= 0.070])
min(Result$N[Result$CIwidth <= 0.080])
min(Result$N[Result$CIwidth <= 0.090])
min(Result$N[Result$CIwidth <= 0.100])

#### 3.2 ADNEX ####
Cstat = 0.90 # Werk met laagste AUC-waarde: sample size voor AUC is meestal groter bij de veronderstelling van een lage AUC
# No closed form solution -> an iterative or deductive approach
# Calculate the SE for a range of N
N <- seq(1, 1000000, 1)
seCstatsq = Cstat * (1 - Cstat) * (1 + (((N/2) - 1) * ((1 - Cstat)/(2 - Cstat))) + ((((N/2) - 1) * Cstat)/(1 + Cstat)))/(N^2 * outcome_prop * (1 - outcome_prop))
seCstat = sqrt(seCstatsq)
CIwidth = 2*1.96*seCstat 
# Identify the sample sizes that give a CIwidth no wider than the desired value
Result <- data.frame(cbind(N, CIwidth))
min(Result$N[Result$CIwidth <= 0.050])
min(Result$N[Result$CIwidth <= 0.060])
min(Result$N[Result$CIwidth <= 0.070])
min(Result$N[Result$CIwidth <= 0.080])
min(Result$N[Result$CIwidth <= 0.090])
min(Result$N[Result$CIwidth <= 0.100])


#### 4. Net benefit ####

#### 4.1 ADNEX ####
outcome_prop = 0.49
sens = c(0.93, 0.91, 0.87, 0.83, 0.79, 0.77, 0.72, 0.66)
spec = c(0.74, 0.84, 0.88, 0.90, 0.92, 0.93, 0.95, 0.96)
threshold = c(0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50)
NB = (sens * outcome_prop) - ((1 - spec) * (1 - outcome_prop) * (threshold/(1 - threshold)))
sNB = NB/outcome_prop
NB
sNB

w = ((1 - outcome_prop)/outcome_prop) * (threshold/(1 - threshold))
# target CI width for sNB of 0.2, which corresponds to SE of 0.051
sesNB = 0.026
sampsize_sNB = (1/(sesNB^2)) * ((sens * (1 - sens)/outcome_prop) + (w^2 * spec * (1 - spec)/(1 - outcome_prop)) + (w^2 * (1 - spec)^2/(outcome_prop * (1 - outcome_prop))))
sampsize_sNB

sesNB = 0.038
sampsize_sNB = (1/(sesNB^2)) * ((sens * (1 - sens)/outcome_prop) + (w^2 * spec * (1 - spec)/(1 - outcome_prop)) + (w^2 * (1 - spec)^2/(outcome_prop * (1 - outcome_prop))))
sampsize_sNB

sesNB = 0.051
sampsize_sNB = (1/(sesNB^2)) * ((sens * (1 - sens)/outcome_prop) + (w^2 * spec * (1 - spec)/(1 - outcome_prop)) + (w^2 * (1 - spec)^2/(outcome_prop * (1 - outcome_prop))))
sampsize_sNB


#### 4.2 ROMA ####
ROMA <- read_excel("C:/Users/u0123496/OneDrive - KU Leuven/IOTA/ROMA vs ADNEX/Sample size calculation/roma risks Kaijser et al 2013.xlsx")

# Calculate sensitivity and specificity based on previous data
threshold = 0.50
predicted_values<-ifelse(ROMA$ROMAfuji_risk > threshold, 1, 0)
actual_values <- ROMA$outcome1
conf_matrix <- table(predicted_values, actual_values)
conf_matrix
sensitivity(conf_matrix, positive = "1")
specificity(conf_matrix, negative = "0")

outcome_prop = 0.49
sens = c(0.96, 0.88, 0.79, 0.72, 0.69, 0.66, 0.59, 0.53)
spec = c(0.16, 0.51, 0.75, 0.83, 0.88, 0.92, 0.95, 0.97)
threshold = c(0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50)
NB = (sens * outcome_prop) - ((1 - spec) * (1 - outcome_prop) * (threshold/(1 - threshold)))
sNB = NB/outcome_prop
NB
sNB

w = ((1 - outcome_prop)/outcome_prop) * (threshold/(1 - threshold))
# target CI width for sNB of 0.2, which corresponds to SE of 0.051
sesNB = 0.026
sampsize_sNB = (1/(sesNB^2)) * ((sens * (1 - sens)/outcome_prop) + (w^2 * spec * (1 - spec)/(1 - outcome_prop)) + (w^2 * (1 - spec)^2/(outcome_prop * (1 - outcome_prop))))
sampsize_sNB

sesNB = 0.038
sampsize_sNB = (1/(sesNB^2)) * ((sens * (1 - sens)/outcome_prop) + (w^2 * spec * (1 - spec)/(1 - outcome_prop)) + (w^2 * (1 - spec)^2/(outcome_prop * (1 - outcome_prop))))
sampsize_sNB

sesNB = 0.051
sampsize_sNB = (1/(sesNB^2)) * ((sens * (1 - sens)/outcome_prop) + (w^2 * spec * (1 - spec)/(1 - outcome_prop)) + (w^2 * (1 - spec)^2/(outcome_prop * (1 - outcome_prop))))
sampsize_sNB

