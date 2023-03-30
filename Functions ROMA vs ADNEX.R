#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#   Functions ROMA vs ADNEX analyses    #
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

#### Packages
library(car)
library(metafor)
library(metamisc)
library(boot)
library(auRoc)
library(data.table)
library(SDMTools)
library(doBy)
library(matrixStats)
library(Metatron)
library(Hmisc)
library(plotrix)
library(rms)


###########################
### 1. Function for AUC ###
###########################

AUC.IOTA <- function(pred, outcome, center, data, method.MA = "BAYES", titleGraph = "AUC per center"){
  
  arguments <- as.list(match.call())[-1]
  pred = eval(arguments$pred, data)
  outcome <- eval(arguments$outcome, data)
  center <- eval(arguments$center, data)
  
  Df = data.frame(p = pred, y = outcome, center = center, stringsAsFactors = F)
  
  centers <- unique(Df$center)
  
  AUCcenter <- matrix(ncol = 6, nrow = length(centers))
  colnames(AUCcenter) <- c('Center', 'SampleSize', 'Prevalence', 'AUC', 'LL', 'UL')
  AUCcenter <- data.frame(AUCcenter)
  AUCcenter$Center <- centers
  
  # AUC per center
  for(i in seq_along(centers)){
    AUCcenter[i, 4:6] <- auc.nonpara.mw(Df$p[Df$center == centers[i] & Df$y == 1], Df$p[Df$center == centers[i] & Df$y == 0], method = "pepe")
    AUCcenter[i, 2]   <- nrow(Df[Df$center == centers[i],])
    AUCcenter[i, 3]   <- round(nrow(Df[Df$y == 1 & Df$center == centers[i],])/nrow(Df[Df$center == centers[i],])*100)
    
    ## Additional part for AUCs of 1
    if(AUCcenter[i, 4] == 1){
      AUCcenter[i, 4:6] <- auc.nonpara.mw(Df$p[Df$center == centers[i] & Df$y == 1], Df$p[Df$center == centers[i] & Df$y == 0], method = "newcombe") # Newcombe ipv pepe
    } else{
      AUCcenter[i, 4:6] <- auc.nonpara.mw(Df$p[Df$center == centers[i] & Df$y == 1], Df$p[Df$center == centers[i] & Df$y == 0], method = "pepe")
    }
    
    if(AUCcenter$AUC[i] != 1){
      AUCcenter$logit.AUC[i] <- logit(AUCcenter$AUC[i])
      AUCcenter$logit.se[i]  <- (logit(AUCcenter$AUC[i]) - logit(AUCcenter$LL[i]))/1.96
      AUCcenter$logit.var[i] <- AUCcenter$logit.se[i]^2
    } else{
      AUCcenter$logit.AUC[i] <- logit(0.999)
      AUCcenter$logit.se[i]  <- (logit(0.999) - logit(AUCcenter$LL[i]))/1.96
      AUCcenter$logit.var[i] <- AUCcenter$logit.se[i]^2
    }
  }
  
  
  AUCcenter <- AUCcenter[order(AUCcenter$SampleSize, decreasing = TRUE),]
  
  AUCoverall <- matrix(nrow = 2, ncol = 6)
  colnames(AUCoverall) <- c('Center', 'SampleSize', 'Prevalence', 'AUC', 'LL', 'UL')
  AUCoverall <- data.frame(AUCoverall)
  AUCoverall$Center <- c("Meta-analysis", "95% Prediction interval")
  AUCoverall$SampleSize <- nrow(Df)
  AUCoverall$Prevalence <- round(nrow(Df[Df$y == 1,])/nrow(Df)*100)
  
  # Meta-analyse for overall estimate
  fit.RE = uvmeta(r = AUCcenter$logit.AUC, r.se = AUCcenter$logit.se, method = method.MA, labels = AUCcenter$Center)
  
  AUCoverall$AUC[1] <- inv.logit(fit.RE$est)
  AUCoverall$LL[1]  <- inv.logit(fit.RE$ci.lb)
  AUCoverall$UL[1]  <- inv.logit(fit.RE$ci.ub)
  AUCoverall$AUC[2] <- inv.logit(fit.RE$est)
  AUCoverall$LL[2]  <- inv.logit(fit.RE$pi.lb)
  AUCoverall$UL[2]  <- inv.logit(fit.RE$pi.ub)
  AUCoverall
  
  NAforest <- matrix(nrow = 1, ncol = 6)
  colnames(NAforest) <- c('Center', 'SampleSize', 'Prevalence', 'AUC', 'LL', 'UL')
  NAforest <- data.frame(NAforest)
  NAforest <- NA
  
  AUC <- rbind(AUCcenter[, 1:6], NAforest, NAforest, AUCoverall)
  
  # Layout forestplot
  RRcenter <- c('Other', 'Meta-analysis')
  nrobs <- nrow(AUC)
  RRcenter <- 1:nrobs
  for(i in 1:nrobs){
    RRcenter[i] <- AUC$Center[i]
  }
  RRauc <- 1:nrobs
  for(i in 1:nrobs){
    RRauc[i] <- paste(format(round(AUC$AUC[i], 2), nsmall = 2), " (", format(round(AUC$LL[i], 2), nsmall = 2), " to ", format(round(AUC$UL[i], 2), nsmall = 2), ")", sep = "")
  }
  RRprev <- 1:nrobs
  for(i in 1:nrobs){
    RRprev[i] <- AUC$SampleSize[i]
  }
  
  Labels <- c('Centre', 'AUC (95% CI)', 'N')
  Combined <- cbind(RRcenter, RRauc, RRprev)
  Tabletext <- rbind(Labels, NAforest, Combined)
  
  AUCna <- rbind(NAforest, NAforest, AUC)
  
  return(structure(list(Performance = AUCoverall, ModelFit = fit.RE, AUCcenters = AUCcenter, IncludedCenters = centers, dataPlot = AUCna, Plot = Tabletext, data = Df)))
  
}

#####################################################
#### 2. Function for sensitivity and specificity ####
#####################################################

SensSpec <- function(pred, outcome, threshold, center, data, method.MA = "REML", c = 0.5){
  
  arguments <- as.list(match.call())[-1]
  pred = eval(arguments$pred, data)
  outcome <- eval(arguments$outcome, data)
  center <- eval(arguments$center, data)
  threshold = threshold
  
  Df = data.frame(p = pred, y = outcome, center = center, stringsAsFactors = F)
  
  centers <- unique(Df$center)
  
  # Sensitivity and specificity per center
  Senscenter <- matrix(ncol = 7, nrow = length(centers))
  colnames(Senscenter) <- c("Center", "TP", "TN", "FP", "FN", "Sensitivity", "Specificity")
  Senscenter <- data.frame(Senscenter)
  Senscenter$Center <- centers
  
  for(i in seq_along(centers)){
    predicted_values <- ifelse(Df$p[Df$center == centers[i]] >= threshold, 1, 0)
    actual_values <- Df$y[Df$center == centers[i]]
    conf_matrix <- table(factor(predicted_values, levels = c(0, 1)), factor(actual_values, levels = c(0, 1)))
    Senscenter$TN[i]     <- conf_matrix[1,1]
    Senscenter$TP[i]     <- conf_matrix[2,2]
    Senscenter$FN[i]     <- conf_matrix[1,2]
    Senscenter$FP[i]     <- conf_matrix[2,1]
    
    if(any(Senscenter[i, 2:5] == 0)){
      Senscenter$TN[i]     <- conf_matrix[1,1] + c
      Senscenter$TP[i]     <- conf_matrix[2,2] + c
      Senscenter$FN[i]     <- conf_matrix[1,2] + c
      Senscenter$FP[i]     <- conf_matrix[2,1] + c
    }
    
  }
  
  Senscenter$Sensitivity  <- Senscenter$TP / (Senscenter$TP + Senscenter$FN)
  Senscenter$se_sens      <- sqrt((Senscenter$Sensitivity * (1 - Senscenter$Sensitivity))/(Senscenter$TP + Senscenter$FN))
  Senscenter$Specificity  <- Senscenter$TN / (Senscenter$TN + Senscenter$FP)
  Senscenter$se_spec      <- sqrt((Senscenter$Specificity * (1 - Senscenter$Specificity))/(Senscenter$TN + Senscenter$FP))
  
  # Logit transformation
  Senscenter$logit.sens     <- logit(Senscenter$Sensitivity)
  Senscenter$logit.se.sens  <- sqrt(1/Senscenter$TP + 1/Senscenter$FN)
  Senscenter$logit.var.sens <- Senscenter$logit.se.sens^2
  Senscenter$logit.LL.sens  <- Senscenter$logit.sens - 1.96*Senscenter$logit.se.sens
  Senscenter$logit.UL.sens  <- Senscenter$logit.sens + 1.96*Senscenter$logit.se.sens
  
  Senscenter$logit.spec     <- logit(Senscenter$Specificity)
  Senscenter$logit.se.spec  <- sqrt(1/Senscenter$TN + 1/Senscenter$FP)
  Senscenter$logit.var.spec <- Senscenter$logit.se.spec^2
  Senscenter$logit.LL.spec  <- Senscenter$logit.spec - 1.96*Senscenter$logit.se.spec
  Senscenter$logit.UL.spec  <- Senscenter$logit.spec + 1.96*Senscenter$logit.se.spec
  
  Sensoverall <- matrix(nrow = 1, ncol = 7)
  colnames(Sensoverall) <- c('Center', 'Sens', 'LL.sens', 'UL.sens', 'Spec', 'LL.spec', 'UL.spec')
  Sensoverall <- data.frame(Sensoverall)
  Sensoverall$Center <- "Overall"
  
  # Meta-analyse for overall estimate: Bivariate random-effects model
  p1 <- 2*length(centers)
  p2 <- 4*length(centers)
  Sensbi <- Senscenter[, c("Center", "logit.sens", "logit.spec", "logit.se.sens", "logit.se.spec")]
  Sensbi.long <- melt(Sensbi, id.vars = "Center")
  Sensbi2 <- cbind(Sensbi.long[1:p1,], Sensbi.long[(p1 + 1):p2,])
  colnames(Sensbi2) <- c("Center", "SensSpec", "Value", "Centre", "Pooled", "SE")
  Sensbi2$VAR <- Sensbi2$SE^2
  Sensbi2$SensSpec <- factor(Sensbi2$SensSpec)
  fit.RE = rma.mv(yi = Value, V = VAR, mods = ~ SensSpec-1, random = ~ SensSpec-1 | Center, struct="UN", data=Sensbi2)
  Sensoverall$Sens <- inv.logit(coef(fit.RE)[1])
  Sensoverall$LL.sens <- inv.logit(fit.RE$ci.lb[1])
  Sensoverall$UL.sens <- inv.logit(fit.RE$ci.ub[1])
  Sensoverall$Spec <- inv.logit(coef(fit.RE)[2])
  Sensoverall$LL.spec <- inv.logit(fit.RE$ci.lb[2])
  Sensoverall$UL.spec <- inv.logit(fit.RE$ci.ub[2])
  
  Sensoverall$Threshold   <- threshold
  Sensoverall$Sensitivity <- paste0(format(round(Sensoverall$Sens, 3), nsmall = 3), " (", format(round(Sensoverall$LL.sens, 3), nsmall = 3), " - ", format(round(Sensoverall$UL.sens, 3), nsmall = 3), ")")
  Sensoverall$Specificity <- paste0(format(round(Sensoverall$Spec, 3), nsmall = 3), " (", format(round(Sensoverall$LL.spec, 3), nsmall = 3), " - ", format(round(Sensoverall$UL.spec, 3), nsmall = 3), ")")
  
  return(structure(list(OverallPer = Sensoverall[, 8:10], CenterPer = Senscenter, IncludedCenter = centers, data = Df)))
}

## Fixed sensitivity
FixSens <- function(pred, outcome, Sensitivity, center, data, method.MA = "REML", c = 0.5){
  
  arguments <- as.list(match.call())[-1]
  pred = eval(arguments$pred, data)
  outcome <- eval(arguments$outcome, data)
  center <- eval(arguments$center, data)
  Sensitivity = Sensitivity
  
  Df = data.frame(p = pred, y = outcome, center = center, stringsAsFactors = F)
  
  centers <- unique(Df$center)
  
  # Sensitivity and specificity per center
  Senscenter <- matrix(ncol = 3, nrow = length(centers))
  colnames(Senscenter) <- c("Center", "Sensitivity", "Specificity")
  Senscenter <- data.frame(Senscenter)
  Senscenter$Center <- centers
  Senscenter$Sensitivity <- Sensitivity
  
  for(i in seq_along(centers)){
    threshold <- quantile(Df$p[Df$center == centers[i] & Df$y == 1], probs = 1 - Sensitivity)
    
    predicted_values <- ifelse(Df$p[Df$center == centers[i]] >= threshold, 1, 0)
    actual_values <- Df$y[Df$center == centers[i]]
    conf_matrix <- table(factor(predicted_values, levels = c(0, 1)), factor(actual_values, levels = c(0, 1)))
    Senscenter$TN[i]     <- conf_matrix[1,1]
    Senscenter$TP[i]     <- conf_matrix[2,2]
    Senscenter$FN[i]     <- conf_matrix[1,2]
    Senscenter$FP[i]     <- conf_matrix[2,1]
    
    if(any(Senscenter[i, 4:7] == 0)){
      Senscenter$TN[i]     <- conf_matrix[1,1] + c
      Senscenter$TP[i]     <- conf_matrix[2,2] + c
      Senscenter$FN[i]     <- conf_matrix[1,2] + c
      Senscenter$FP[i]     <- conf_matrix[2,1] + c
    }
    
  }
  
  Senscenter$Sensitivity  <- Senscenter$TP / (Senscenter$TP + Senscenter$FN)
  Senscenter$se_sens      <- sqrt((Senscenter$Sensitivity * (1 - Senscenter$Sensitivity))/(Senscenter$TP + Senscenter$FN))
  Senscenter$Specificity  <- Senscenter$TN / (Senscenter$TN + Senscenter$FP)
  Senscenter$se_spec      <- sqrt((Senscenter$Specificity * (1 - Senscenter$Specificity))/(Senscenter$TN + Senscenter$FP))
  
  # Logit transformation
  Senscenter$logit.sens     <- logit(Senscenter$Sensitivity)
  Senscenter$logit.se.sens  <- sqrt(1/Senscenter$TP + 1/Senscenter$FN)
  Senscenter$logit.var.sens <- Senscenter$logit.se.sens^2
  Senscenter$logit.LL.sens  <- Senscenter$logit.sens - 1.96*Senscenter$logit.se.sens
  Senscenter$logit.UL.sens  <- Senscenter$logit.sens + 1.96*Senscenter$logit.se.sens
  
  Senscenter$logit.spec     <- logit(Senscenter$Specificity)
  Senscenter$logit.se.spec  <- sqrt(1/Senscenter$TN + 1/Senscenter$FP)
  Senscenter$logit.var.spec <- Senscenter$logit.se.spec^2
  Senscenter$logit.LL.spec  <- Senscenter$logit.spec - 1.96*Senscenter$logit.se.spec
  Senscenter$logit.UL.spec  <- Senscenter$logit.spec + 1.96*Senscenter$logit.se.spec
  
  Sensoverall <- matrix(nrow = 1, ncol = 7)
  colnames(Sensoverall) <- c('Center', 'Sens', 'LL.sens', 'UL.sens', 'Spec', 'LL.spec', 'UL.spec')
  Sensoverall <- data.frame(Sensoverall)
  Sensoverall$Center <- "Overall"
  
  # Meta-analyse for overall estimate: Bivariate random-effects model
  p1 <- 2*length(centers)
  p2 <- 4*length(centers)
  Sensbi <- Senscenter[, c("Center", "logit.sens", "logit.spec", "logit.se.sens", "logit.se.spec")]
  Sensbi.long <- melt(Sensbi, id.vars = "Center")
  Sensbi2 <- cbind(Sensbi.long[1:p1,], Sensbi.long[(p1 + 1):p2,])
  colnames(Sensbi2) <- c("Center", "SensSpec", "Value", "Centre", "Pooled", "SE")
  Sensbi2$VAR <- Sensbi2$SE^2
  Sensbi2$SensSpec <- factor(Sensbi2$SensSpec)
  fit.RE = rma.mv(yi = Value, V = VAR, mods = ~ SensSpec-1, random = ~ SensSpec | Center, struct="UN", data=Sensbi2)
  Sensoverall$Sens <- inv.logit(coef(fit.RE)[1])
  Sensoverall$LL.sens <- inv.logit(fit.RE$ci.lb[1])
  Sensoverall$UL.sens <- inv.logit(fit.RE$ci.ub[1])
  Sensoverall$Spec <- inv.logit(coef(fit.RE)[2])
  Sensoverall$LL.spec <- inv.logit(fit.RE$ci.lb[2])
  Sensoverall$UL.spec <- inv.logit(fit.RE$ci.ub[2])
  
  Sensoverall$Threshold   <- threshold
  Sensoverall$Sensitivity <- paste0(format(round(Sensoverall$Sens, 3), nsmall = 3), " (", format(round(Sensoverall$LL.sens, 3), nsmall = 3), " - ", format(round(Sensoverall$UL.sens, 3), nsmall = 3), ")")
  Sensoverall$Specificity <- paste0(format(round(Sensoverall$Spec, 3), nsmall = 3), " (", format(round(Sensoverall$LL.spec, 3), nsmall = 3), " - ", format(round(Sensoverall$UL.spec, 3), nsmall = 3), ")")
  
  return(structure(list(OverallPer = Sensoverall[, 8:10], CenterPer = Senscenter, IncludedCenter = centers, data = Df)))
}

## Fixed specificity
FixSpec <- function(pred, outcome, Specificity, center, data, method.MA = "REML", c = 0.5){
  
  arguments <- as.list(match.call())[-1]
  pred = eval(arguments$pred, data)
  outcome <- eval(arguments$outcome, data)
  center <- eval(arguments$center, data)
  Specificity = Specificity
  
  Df = data.frame(p = pred, y = outcome, center = center, stringsAsFactors = F)
  
  centers <- unique(Df$center)
  
  # Sensitivity and specificity per center
  Senscenter <- matrix(ncol = 3, nrow = length(centers))
  colnames(Senscenter) <- c("Center", "Sensitivity", "Specificity")
  Senscenter <- data.frame(Senscenter)
  Senscenter$Center <- centers
  Senscenter$Specificity <- Specificity
  
  for(i in seq_along(centers)){
    threshold <- quantile(Df$p[Df$center == centers[i] & Df$y == 0], probs = Specificity)
    
    predicted_values <- ifelse(Df$p[Df$center == centers[i]] >= threshold, 1, 0)
    actual_values <- Df$y[Df$center == centers[i]]
    conf_matrix <- table(factor(predicted_values, levels = c(0, 1)), factor(actual_values, levels = c(0, 1)))
    Senscenter$TN[i]     <- conf_matrix[1,1]
    Senscenter$TP[i]     <- conf_matrix[2,2]
    Senscenter$FN[i]     <- conf_matrix[1,2]
    Senscenter$FP[i]     <- conf_matrix[2,1]
    
    if(any(Senscenter[i, 4:7] == 0)){
      Senscenter$TN[i]     <- conf_matrix[1,1] + c
      Senscenter$TP[i]     <- conf_matrix[2,2] + c
      Senscenter$FN[i]     <- conf_matrix[1,2] + c
      Senscenter$FP[i]     <- conf_matrix[2,1] + c
    }
    
  }
  
  Senscenter$Sensitivity  <- Senscenter$TP / (Senscenter$TP + Senscenter$FN)
  Senscenter$se_sens      <- sqrt((Senscenter$Sensitivity * (1 - Senscenter$Sensitivity))/(Senscenter$TP + Senscenter$FN))
  Senscenter$Specificity  <- Senscenter$TN / (Senscenter$TN + Senscenter$FP)
  Senscenter$se_spec      <- sqrt((Senscenter$Specificity * (1 - Senscenter$Specificity))/(Senscenter$TN + Senscenter$FP))
  
  # Logit transformation
  Senscenter$logit.sens     <- logit(Senscenter$Sensitivity)
  Senscenter$logit.se.sens  <- sqrt(1/Senscenter$TP + 1/Senscenter$FN)
  Senscenter$logit.var.sens <- Senscenter$logit.se.sens^2
  Senscenter$logit.LL.sens  <- Senscenter$logit.sens - 1.96*Senscenter$logit.se.sens
  Senscenter$logit.UL.sens  <- Senscenter$logit.sens + 1.96*Senscenter$logit.se.sens
  
  Senscenter$logit.spec     <- logit(Senscenter$Specificity)
  Senscenter$logit.se.spec  <- sqrt(1/Senscenter$TN + 1/Senscenter$FP)
  Senscenter$logit.var.spec <- Senscenter$logit.se.spec^2
  Senscenter$logit.LL.spec  <- Senscenter$logit.spec - 1.96*Senscenter$logit.se.spec
  Senscenter$logit.UL.spec  <- Senscenter$logit.spec + 1.96*Senscenter$logit.se.spec
  
  Sensoverall <- matrix(nrow = 1, ncol = 7)
  colnames(Sensoverall) <- c('Center', 'Sens', 'LL.sens', 'UL.sens', 'Spec', 'LL.spec', 'UL.spec')
  Sensoverall <- data.frame(Sensoverall)
  Sensoverall$Center <- "Overall"
  
  # Meta-analyse for overall estimate: Bivariate random-effects model
  p1 <- 2*length(centers)
  p2 <- 4*length(centers)
  Sensbi <- Senscenter[, c("Center", "logit.sens", "logit.spec", "logit.se.sens", "logit.se.spec")]
  Sensbi.long <- melt(Sensbi, id.vars = "Center")
  Sensbi2 <- cbind(Sensbi.long[1:p1,], Sensbi.long[(p1 + 1):p2,])
  colnames(Sensbi2) <- c("Center", "SensSpec", "Value", "Centre", "Pooled", "SE")
  Sensbi2$VAR <- Sensbi2$SE^2
  Sensbi2$SensSpec <- factor(Sensbi2$SensSpec)
  fit.RE = rma.mv(yi = Value, V = VAR, mods = ~ SensSpec-1, random = ~ SensSpec | Center, struct="UN", data=Sensbi2)
  Sensoverall$Sens <- inv.logit(coef(fit.RE)[1])
  Sensoverall$LL.sens <- inv.logit(fit.RE$ci.lb[1])
  Sensoverall$UL.sens <- inv.logit(fit.RE$ci.ub[1])
  Sensoverall$Spec <- inv.logit(coef(fit.RE)[2])
  Sensoverall$LL.spec <- inv.logit(fit.RE$ci.lb[2])
  Sensoverall$UL.spec <- inv.logit(fit.RE$ci.ub[2])
  
  Sensoverall$Threshold   <- threshold
  Sensoverall$Sensitivity <- paste0(format(round(Sensoverall$Sens, 3), nsmall = 3), " (", format(round(Sensoverall$LL.sens, 3), nsmall = 3), " - ", format(round(Sensoverall$UL.sens, 3), nsmall = 3), ")")
  Sensoverall$Specificity <- paste0(format(round(Sensoverall$Spec, 3), nsmall = 3), " (", format(round(Sensoverall$LL.spec, 3), nsmall = 3), " - ", format(round(Sensoverall$UL.spec, 3), nsmall = 3), ")")
  
  return(structure(list(OverallPer = Sensoverall[, 8:10], CenterPer = Senscenter[, c(1, 10:15)], IncludedCenter = centers, data = Df)))
}


#####################################
#### 3. Function for calibration ####
#####################################

Calibration <- function (p, y, center, data, CalibrLines = c("overall", "centers", 
                                                             "both"), LogCal = T, flexible = F, dostats = T, statloc = c(0, 
                                                                                                                         0.85), legendloc = c(0.5, 0.27), MethodCL = c("profile", 
                                                                                                                                                                       "Wald", "boot"), LevelCL = 0.95, roundstats = 3, cex = 0.75, 
                         cex.leg = 0.75, ncol.leg = 1, lty.overall = 1, lwd.overall = 2, 
                         col.overall = "red", RMprompt = F, xlab = "Estimated risk", nr.knots = 5,
                         ylab = "Observed proportion", xlim = c(0, 1), ylim = c(0, 
                                                                                1), d0lab = "0", d1lab = "1", cex.d01 = 0.7, dist.label = 0.04, 
                         line.bins = -0.05, dist.label2 = 0.03, las = 1, length.seg = 1, 
                         y.intersp = 1, lty.ideal = 1, col.ideal = "black", lwd.ideal = 1.75, 
                         lty.centers = NULL, lwd.centers = NULL, col.centers = NULL, 
                         EpsGrad = 0.001, fNonConv = c("warning", "stop"), Controlglmer = glmerControl(optimizer = "bobyqa"), 
                         ...) 
{
  Argz = as.list(match.call())[-1]
  CalibrLines = match.arg(CalibrLines)
  MethCL = match.arg(MethodCL)
  p = eval(Argz$p, data)
  y = eval(Argz$y, data)
  center = eval(Argz$center, data)
  LP = Logit(p)
  if (is.factor(center)) 
    center = as.character(center)
  if (length(unique(center)) == 1) 
    stop("There is only one center, hence a random effects model should not be used here.")
  if (length(unique(center)) < 5) 
    warning("There are less than 5 centers, consider using a different method.", 
            immediate. = T)
  if (LogCal & flexible) 
    stop("LogCal and flexible cannot both be true.")
  Df = data.frame(y = y, p = p, LP = LP, center = center)
  Df = Df[with(Df, order(center, p)), ]
  RmCenter = dlply(Df, .(center), function(x) if (sum(x$y == 
                                                      1) < 10 | sum(x$y == 0) < 10) 
    unique(x$center)
    else NULL)
  RmCenter = unname(unlist(RmCenter))
  if (!is.null(RmCenter)) {
    if (RMprompt) {
      Cmmnd = readline(paste0("The center(s) ", paste0(RmCenter, 
                                                       collapse = ", "), " have less than 10 (non-)events and these will be", 
                              " removed. Do you want to continue? (y/n)   "))
      if (Cmmnd == "n") 
        stop("Function was stopped by the user.")
    }
    Df = Df[!(Df$center %in% RmCenter), ]
  }
  IncludedCenters = unique(Df$center)
  centers = IncludedCenters
  nCenters = length(IncludedCenters)
  if (CalibrLines != "overall") {
    if (all(sapply(list(lty.centers, lwd.centers, col.centers), 
                   is.null))) {
      lty.centers = rep(1, nCenters)
      lwd.centers = rep(1, nCenters)
      col.centers = seq_len(nCenters)
    }
    else if (sum(!sapply(list(lty.centers, lwd.centers, col.centers), 
                         is.null)) != 3) {
      stop("If you specify one of the arguments, please also specify the others.")
    }
    else {
      if (any(sapply(list(lty.centers, lwd.centers, col.centers), 
                     function(x) length(x) < nCenters))) 
        stop("The vector length of lty.centers, lwd.centers and col.centers is less than the number of centers.")
      FixL = function(x) x[1:nCenters]
      lty.centers = FixL(lty.centers)
      col.centers = FixL(col.centers)
      lwd.centers = FixL(lwd.centers)
    }
  }
  Df$center = factor(Df$center)
  contrasts(Df$center) = "contr.sum"
  LogMM = glmer(y ~ LP + (LP | center), data = Df, family = "binomial", 
                control = Controlglmer)
  CalSlope = fixef(LogMM)[2]
  cat("\n\nComputing confidence interval for calibration slope...\n\n")
  CalSlope = c(CalSlope, confint.merMod(LogMM, parm = "LP", 
                                        quiet = T, method = MethodCL, level = LevelCL))
  LogMM2 = glmer(y ~ 1 + (1 | center), data = Df, family = "binomial", 
                 offset = LP, control = Controlglmer)
  CalInterc = fixef(LogMM2)
  cat("\n\nComputing confidence interval for calibration intercept...\n\n")
  CalInterc = c(CalInterc, confint(LogMM2, parm = "(Intercept)", 
                                   quiet = T, method = MethodCL, level = LevelCL))

  if(CalibrLines != "overall"){
    par(mar = c(5,5,1,13), xpd=TRUE, pty = 's')
  } else{
    par(pty = 's')
  }
  
  plot(0.5, 0.5, xlim = xlim, ylim = ylim, type = "n", xlab = xlab, 
       ylab = ylab, las = las, ...)
  clip(0, 1, 0, 1)
  abline(0, 1, lty = lty.ideal, col = col.ideal, lwd = lwd.ideal)
  do.call("clip", as.list(par()$usr))
  lt = lty.ideal
  lw.d = lwd.ideal
  all.col = col.ideal
  leg = "Ideal"
  NewX = data.frame(LP = seq(min(Df$LP), max(Df$LP), length = 500))
  p = Ilogit(NewX$LP)
  X = cbind(1, Logit(p))
  FE = fixef(LogMM)
  if (CalibrLines != "overall") {
    PerC = data.frame(LP = rep(NewX$LP, length(IncludedCenters)), 
                      center = sort(rep(IncludedCenters, 500)))
    EstC = predict(LogMM, newdata = PerC, re.form = ~(LP | 
                                                        center), allow.new.levels = T, type = "response")
    ResultsC = cbind.data.frame(x = Ilogit(PerC$LP), y = EstC, 
                                center = PerC$center)
    ResultsC = split(ResultsC, ResultsC$center)
  }
  if (CalibrLines == "centers") {
    for (i in IncludedCenters) lines(ResultsC[[i]], col = col.centers[which(IncludedCenters == 
                                                                              i)], lty = lty.centers[which(IncludedCenters == i)], 
                                     lwd = lwd.centers[which(IncludedCenters == i)])
    lt = c(lt, lty.centers)
    lw.d = c(lw.d, lwd.centers)
    all.col = c(all.col, col.centers)
    leg = c(leg, as.character(IncludedCenters))
  }
  else {
    p = Ilogit(NewX$LP)
    X = cbind(1, Logit(p))
    FE = fixef(LogMM)
    OverallCal = Ilogit(X[order(p), ] %*% FE)
    
    options(datadist=NULL)
    
    nr.knots <<- nr.knots
    f <- ols(OverallCal ~ rcs(p, nr.knots)) #https://stackoverflow.com/questions/35331184/how-to-plot-ols-with-r-c-splines
    ddist <<- datadist(p) # https://stackoverflow.com/questions/30634766/nomogram-of-rms-not-working-in-a-function
    options(datadist = "ddist")
    PredictAverage <- Predict(f)

    lines(PredictAverage$p, PredictAverage$yhat, lwd = lwd.overall, lty = lty.overall, 
          col = col.overall)

    lt = c(lt, lty.overall)
    lw.d = c(lw.d, lwd.overall)
    all.col = c(all.col, col.overall)
    leg = c(leg, "Overall")
    if (CalibrLines == "both") {
      for (i in IncludedCenters) lines(ResultsC[[i]], col = col.centers[which(IncludedCenters == 
                                                                                i)], lty = lty.centers[which(IncludedCenters == 
                                                                                                               i)], lwd = lwd.centers[which(IncludedCenters == 
                                                                                                                                              i)])
      lt = c(lt, lty.centers)
      lw.d = c(lw.d, lwd.centers)
      all.col = c(all.col, col.centers)
      leg = c(leg, as.character(IncludedCenters))
    }
  }
  lp = legendloc
  lp = list(x = lp[1], y = lp[2])
  legend("topright", leg, lty = lt, cex = cex.leg, bty = "n", lwd = lw.d, 
         col = all.col, y.intersp = y.intersp, ncol = ncol.leg, inset = c(-.68, 0), xpd=NA) 
  if (dostats) {
    if (CalibrLines == "centers"){
      stats.2 <- paste("")
      text(statloc[1], statloc[2], stats.2, pos = 4, cex = 0.75)
    } else{
      stats.2 <- matrix(ncol = 2, nrow = 2)
      colnames(stats.2) <- c("", "Estimate (95% CI)")
      stats.2[1, ] <- c("Intercept", paste0(format(round(CalInterc[1], 2), nsmall = 2), " (", format(round(CalInterc[2], 2), nsmall = 2), " to ", format(round(CalInterc[3], 2), nsmall = 2), ")"))
      stats.2[2, ] <- c("Slope", paste0(format(round(CalSlope[1], 2), nsmall = 2), " (", format(round(CalSlope[2], 2), nsmall = 2), " to ", format(round(CalSlope[3], 2), nsmall = 2), ")"))
      
      addtable2plot(x = statloc[1], y = statloc[2], table = stats.2, display.colnames = TRUE, cex = 0.75)
      
    }
  }
  Performance = rbind(CalInterc, CalSlope) 
  rownames(Performance) = c("Calibration intercept", "Calibration slope") 
  colnames(Performance) = c("Point estimate", "LCL", "UCL")
  return(structure(list(Performance = Performance, included = unique(Df$center), FitCalIntercept = LogMM2, 
                        FitCalSlope = LogMM, ConfLevel = LevelCL, Prediction = PredictAverage), class = "RE_ValProb"))
}

########################################
#### 4. Function for decision curve ####
########################################

DataWinBugs <- function(pred, outcome, center, data, 
                        sequence = c(0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5)){
  
  arguments <- as.list(match.call())[-1]
  pred = eval(arguments$pred, data)
  outcome <- eval(arguments$outcome, data)
  center <- eval(arguments$center, data)
  
  Df = data.frame(p = pred, y = outcome, center = center, stringsAsFactors = F)
  
  centers <- unique(Df$center)
  
  ConfusionList <- list()
  
  for(i in 1:length(sequence)){
    ConfusionList[[i]] <- list()
  }
  
  
  for(i in 1:length(sequence)){
    threshold <- sequence[i]
    
    Confusion <- matrix(nrow = length(centers), ncol = 8)
    Confusion <- data.frame(Confusion)
    colnames(Confusion) <- c('Center', 'CutOff', 'TN', 'TP', 'FP', 'FN', 'cases', 'controls')
    Confusion$CutOff <- threshold
    
    for(j in seq_along(centers)){
      
      Confusion$Center[j] <- centers[j]
      
      CM <- confusion.matrix(obs = Df$y[Df$center == centers[j]], pred = Df$p[Df$center == centers[j]], threshold = threshold)
      Confusion$TN[j] <- CM[1,1]
      Confusion$TP[j] <- CM[2,2]
      Confusion$FP[j] <- CM[2,1]
      Confusion$FN[j] <- CM[1,2]
      
      Confusion$cases <- Confusion$TP + Confusion$FN
      Confusion$controls <- Confusion$TN + Confusion$FP
      
      Confusion$n <- Confusion$cases + Confusion$controls
      Confusion$NB <- Confusion$TP / Confusion$n - Confusion$FP / Confusion$n * (threshold / (1 - threshold))
      
    }
    
    ConfusionList[[i]] <- Confusion
    
    
  }
  return(structure(list(Results = ConfusionList)))
}


#############################
#### 5. Function for PDI ####
#############################

pdI.extented <- 
  function (y, d, method = "multinom", k = 3, ...){
    num = k
    option = method
    if (num == 3) {							# For three categories
      y = as.numeric(y)
      d = data.matrix(d)
      n1 = which(y == 1)
      n2 = which(y == 2)
      n3 = which(y == 3)
      if (option == "multinom") {
        fit <- nnet::multinom(y ~ d, ...)
        predict.test.probs <- predict(fit, type = "probs")
        predict.test.df <- data.frame(predict.test.probs)
        pp = predict.test.df
      }
      else if (option == "tree") {
        y <- as.factor(y)
        fit <- rpart::rpart(y ~ d, ...)
        predict.test.probs <- predict(fit, type = "prob")
        predict.test.df <- data.frame(predict.test.probs)
        pp = predict.test.df
      }
      else if (option == "svm") {
        y <- as.factor(y)
        fit <- e1071::svm(y ~ d, ..., probability = T)
        predict.test <- predict(fit, d, probability = T)
        predict.test <- attr(predict.test, "probabilities")
        predict.test.df <- data.frame(predict.test)
        pp = predict.test.df[c("X1", "X2", "X3")]
      }
      else if (option == "lda") {
        fit <- MASS::lda(y ~ d, ...)
        predict.test.probs <- predict(fit, type = "probs")
        predict.test.fit <- predict(fit)
        predict.test <- predict.test.fit$posterior
        predict.test.df <- data.frame(predict.test)
        pp = predict.test.df
      }
      else if (option == "prob") {
        pp_sum <- apply(d, 1, sum)
        a <- pp_sum < 0.999 | pp_sum > 1.001
        b <- sum(a)
        if (b != 0) {
          cat("ERROR: The input value \"d\" should be a probability matrix.")
          return(NULL)
        }
        pp = d
      }
      pv = pp
      pv1 = pv[n1, ]
      pv2 = pv[n2, ]
      pv3 = pv[n3, ]
      pdi1 <- 0
      pdi2 <- 0
      pdi3 <- 0
      for (i in 1:length(n1)) {
        pdi1 = pdi1 + sum(pv1[i, 1] > pv2[, 1]) * sum(pv1[i, 
                                                          1] > pv3[, 1])
      }
      for (i in 1:length(n2)) {
        pdi2 = pdi2 + sum(pv2[i, 2] > pv1[, 2]) * sum(pv2[i, 
                                                          2] > pv3[, 2])
      }
      for (i in 1:length(n3)) {
        pdi3 = pdi3 + sum(pv3[i, 3] > pv1[, 3]) * sum(pv3[i, 
                                                          3] > pv2[, 3])
      }
      pdi <- (pdi1 + pdi2 + pdi3)/(3 * length(n1) * length(n2) * 
                                     length(n3))
      return(pdi)
    }
    else if (num == 4) {							# For four categories
      y = as.numeric(y)
      d = data.matrix(d)
      n1 = which(y == 1)
      n2 = which(y == 2)
      n3 = which(y == 3)
      n4 = which(y == 4)
      if (option == "multinom") {
        fit <- nnet::multinom(y ~ d, ...)
        predict.test.probs <- predict(fit, type = "probs")
        predict.test.df <- data.frame(predict.test.probs)
        pp = predict.test.df
      }
      else if (option == "tree") {
        y <- as.factor(y)
        fit <- rpart::rpart(y ~ d, ...)
        predict.test.probs <- predict(fit, type = "prob")
        predict.test.df <- data.frame(predict.test.probs)
        pp = predict.test.df
      }
      else if (option == "svm") {
        y <- as.factor(y)
        fit <- e1071::svm(y ~ d, ..., probability = T)
        predict.test <- predict(fit, d, probability = T)
        predict.test <- attr(predict.test, "probabilities")
        predict.test.df <- data.frame(predict.test)
        pp = predict.test.df[c("X1", "X2", "X3", "X4")]
      }
      else if (option == "lda") {
        fit <- MASS::lda(y ~ d, ...)
        predict.test.probs <- predict(fit, type = "probs")
        predict.test.fit <- predict(fit)
        predict.test <- predict.test.fit$posterior
        predict.test.df <- data.frame(predict.test)
        pp = predict.test.df
      }
      else if (option == "prob") {
        pp_sum <- apply(d, 1, sum)
        a <- pp_sum < 0.999 | pp_sum > 1.001
        b <- sum(a)
        if (b != 0) {
          cat("ERROR: The input value \"d\" should be a probability matrix.")
          return(NULL)
        }
        pp = d
      }
      pv = pp
      pv1 = pv[n1, ]
      pv2 = pv[n2, ]
      pv3 = pv[n3, ]
      pv4 = pv[n4, ]
      pdi1 <- 0
      pdi2 <- 0
      pdi3 <- 0
      pdi4 <- 0
      for (i in 1:length(n1)) {
        pdi1 = pdi1 + sum(pv1[i, 1] > pv2[, 1]) * sum(pv1[i, 
                                                          1] > pv3[, 1]) * sum(pv1[i, 1] > pv4[, 1])
      }
      for (i in 1:length(n2)) {
        pdi2 = pdi2 + sum(pv2[i, 2] > pv1[, 2]) * sum(pv2[i, 
                                                          2] > pv3[, 2]) * sum(pv2[i, 2] > pv4[, 2])
      }
      for (i in 1:length(n3)) {
        pdi3 = pdi3 + sum(pv3[i, 3] > pv1[, 3]) * sum(pv3[i, 
                                                          3] > pv2[, 3]) * sum(pv3[i, 3] > pv4[, 3])
      }
      for (i in 1:length(n4)) {
        pdi4 = pdi4 + sum(pv4[i, 4] > pv1[, 4]) * sum(pv4[i, 
                                                          4] > pv2[, 4]) * sum(pv4[i, 4] > pv3[, 4])
      }
      pdi <- (pdi1 + pdi2 + pdi3 + pdi4)/(4 * length(n1) * 
                                            length(n2) * length(n3) * length(n4))
      return(pdi)
    }
    else if (num == 5) {							# For five categories
      y = as.numeric(y)
      d = data.matrix(d)
      n1 = which(y == 1)
      n2 = which(y == 2)
      n3 = which(y == 3)
      n4 = which(y == 4)
      n5 = which(y == 5)
      if (option == "multinom") {
        fit <- nnet::multinom(y ~ d, ...)
        predict.test.probs <- predict(fit, type = "probs")
        predict.test.df <- data.frame(predict.test.probs)
        pp = predict.test.df
      }
      else if (option == "tree") {
        y <- as.factor(y)
        fit <- rpart::rpart(y ~ d, ...)
        predict.test.probs <- predict(fit, type = "prob")
        predict.test.df <- data.frame(predict.test.probs)
        pp = predict.test.df
      }
      else if (option == "svm") {
        y <- as.factor(y)
        fit <- e1071::svm(y ~ d, ..., probability = T)
        predict.test <- predict(fit, d, probability = T)
        predict.test <- attr(predict.test, "probabilities")
        predict.test.df <- data.frame(predict.test)
        pp = predict.test.df[c("X1", "X2", "X3", "X4", "X5")]
      }
      else if (option == "lda") {
        fit <- MASS::lda(y ~ d, ...)
        predict.test.probs <- predict(fit, type = "probs")
        predict.test.fit <- predict(fit)
        predict.test <- predict.test.fit$posterior
        predict.test.df <- data.frame(predict.test)
        pp = predict.test.df
      }
      else if (option == "prob") {
        pp_sum <- apply(d, 1, sum)
        a <- pp_sum < 0.999 | pp_sum > 1.001
        b <- sum(a)
        if (b != 0) {
          cat("ERROR: The input value \"d\" should be a probability matrix.")
          return(NULL)
        }
        pp = d
      }
      pv = pp
      pv1 = pv[n1, ]
      pv2 = pv[n2, ]
      pv3 = pv[n3, ]
      pv4 = pv[n4, ]
      pv5 = pv[n5, ]
      pdi1 <- 0
      pdi2 <- 0
      pdi3 <- 0
      pdi4 <- 0
      pdi5 <- 0
      for (i in 1:length(n1)) {  
        pdi1 = pdi1 + sum(pv1[i, 1] > pv2[, 1]) * sum(pv1[i, 
                                                          1] > pv3[, 1]) * sum(pv1[i, 1] > pv4[, 1]) * as.double(sum(pv1[i, 1] > pv5[, 1]))
      }
      for (i in 1:length(n2)) {
        pdi2 = pdi2 + sum(pv2[i, 2] > pv1[, 2]) * sum(pv2[i, 
                                                          2] > pv3[, 2]) * sum(pv2[i, 2] > pv4[, 2]) * as.double(sum(pv2[i, 2] > pv5[, 2]))
      }
      for (i in 1:length(n3)) {
        pdi3 = pdi3 + sum(pv3[i, 3] > pv1[, 3]) * sum(pv3[i, 
                                                          3] > pv2[, 3]) * sum(pv3[i, 3] > pv4[, 3]) * as.double(sum(pv3[i, 3] > pv5[, 3]))
      }
      for (i in 1:length(n4)) {
        pdi4 = pdi4 + sum(pv4[i, 4] > pv1[, 4]) * sum(pv4[i, 
                                                          4] > pv2[, 4]) * sum(pv4[i, 4] > pv3[, 4]) * as.double(sum(pv4[i, 4] > pv5[, 4]))
      }
      for (i in 1:length(n5)) {
        pdi5 = pdi5 + sum(pv5[i, 5] > pv1[, 5]) * sum(pv5[i, 
                                                          5] > pv2[, 5]) * sum(pv5[i, 5] > pv3[, 5]) * as.double(sum(pv5[i, 5] > pv4[, 5]))
      }
      pdi <- (pdi1 + pdi2 + pdi3 + pdi4 + pdi5)/(5 * length(n1) * 
                                                   length(n2) * length(n3) * length(n4) * length(n5))
      return(pdi)
    }
    
    else if (num == 2) {							# For two categories
      y = as.numeric(y)
      d = data.matrix(d)
      n1 = which(y == 1)
      n2 = which(y == 2)
      if (option == "multinom") {
        fit <- nnet::multinom(y ~ d, ...)
        predict.test.probs <- predict(fit, type = "probs")
        predict.test.df <- data.frame(predict.test.probs)
        pp = predict.test.df
        pp <- data.frame(1 - pp, pp)
      }
      else if (option == "tree") {
        y <- as.factor(y)
        fit <- rpart::rpart(y ~ d, ...)
        predict.test.probs <- predict(fit, type = "prob")
        predict.test.df <- data.frame(predict.test.probs)
        pp = predict.test.df
      }
      else if (option == "svm") {
        y <- as.factor(y)
        fit <- e1071::svm(y ~ d, ..., probability = T)
        predict.test <- predict(fit, d, probability = T)
        predict.test <- attr(predict.test, "probabilities")
        predict.test.df <- data.frame(predict.test)
        pp = predict.test.df[c("X1", "X2")]
      }
      else if (option == "lda") {
        fit <- MASS::lda(y ~ d, ...)
        predict.test.probs <- predict(fit, type = "probs")
        predict.test.fit <- predict(fit)
        predict.test <- predict.test.fit$posterior
        predict.test.df <- data.frame(predict.test)
        pp = predict.test.df
      }
      else if (option == "prob") {
        pp_sum <- apply(d, 1, sum)
        a <- pp_sum < 0.999 | pp_sum > 1.001
        b <- sum(a)
        if (b != 0) {
          cat("ERROR: The input value \"d\" should be a probability matrix.")
          return(NULL)
        }
        pp = d
      }
      pv = pp
      pv1 = pv[n1, ]
      pv2 = pv[n2, ]
      pdi1 <- 0
      pdi2 <- 0
      for (i in 1:length(n1)) {
        pdi1 = pdi1 + sum(pv1[i, 1] > pv2[, 1]) + sum(pv1[i, 
                                                          1] == pv2[, 1])
      }
      for (i in 1:length(n2)) {
        pdi2 = pdi2 + sum(pv2[i, 2] > pv1[, 2])
      }
      pdi <- (pdi1 + pdi2)/(2 * length(n1) * length(n2))
      return(pdi)
    }
  }


ests.extended <- function (y, d, acc = "hum", level = 0.95, method = "multinom", 
                           k = 3, B = 250, balance = FALSE, ...) 
{
  series = numeric()
  if (acc == "hum") {
    if (balance == FALSE) {
      for (b in 1:B) {
        nn <- length(y)
        id <- sample(1:nn, nn, replace = T)
        while (length(unique(y[id])) < k) {
          id <- sample(1:nn, nn, replace = T)
        }
        while (min(table(y[id])) < 2) {
          id <- sample(1:nn, nn, replace = T)
        }
        if (class(d) == "numeric") {
          series[b] <- hum(y = y[id], d = d[id], method = method, 
                           k = k, ...)
        }
        else {
          series[b] <- hum(y = y[id], d = d[id, ], method = method, 
                           k = k, ...)
        }
      }
    }
    if (balance == TRUE) {
      for (b in 1:B) {
        id <- unlist(caret::createResample(y, times = 1))
        if (class(d) == "numeric") {
          series[b] <- hum(y = y[id], d = d[id], method = method, 
                           k = k, ...)
        }
        else {
          series[b] <- hum(y = y[id], d = d[id, ], method = method, 
                           k = k, ...)
        }
      }
    }
    series.sort <- sort(series)
    return(list(value = hum(y = y, d = d, method = method, 
                            k = k, ...), se = sd(series), interval = c(series.sort[ifelse(B * 
                                                                                            (0.5 - level/2) < 1, 1, B * (0.5 - level/2))], series.sort[B * 
                                                                                                                                                         (0.5 + level/2)])))
  }
  if (acc == "pdi") {
    if (balance == FALSE) {
      for (b in 1:B) {
        nn <- length(y)
        id <- sample(1:nn, nn, replace = T)
        while (length(unique(y[id])) < k) {
          id <- sample(1:nn, nn, replace = T)
        }
        while (min(table(y[id])) < 2) {
          id <- sample(1:nn, nn, replace = T)
        }
        if (class(d) == "numeric") {
          series[b] <- pdI.extented(y = y[id], d = d[id], method = method, 
                                    k = k, ...)
        }
        else {
          series[b] <- pdI.extented(y = y[id], d = d[id, ], method = method, 
                                    k = k, ...)
        }
      }
    }
    if (balance == TRUE) {
      for (b in 1:B) {
        id <- unlist(caret::createResample(y, times = 1))
        if (class(d) == "numeric") {
          series[b] <- pdI.extented(y = y[id], d = d[id], method = method, 
                                    k = k, ...)
        }
        else {
          series[b] <- pdI.extented(y = y[id], d = d[id, ], method = method, 
                                    k = k, ...)
        }
      }
    }
    series.sort <- sort(series)
    return(list(value = pdI.extented(y = y, d = d, method = method, 
                                     k = k, ...), se = sd(series), interval = c(series.sort[ifelse(B * 
                                                                                                     (0.5 - level/2) < 1, 1, B * (0.5 - level/2))], series.sort[B * 
                                                                                                                                                                  (0.5 + level/2)])))
  }
  if (acc == "ccp") {
    if (balance == FALSE) {
      for (b in 1:B) {
        nn <- length(y)
        id <- sample(1:nn, nn, replace = T)
        while (length(unique(y[id])) < k) {
          id <- sample(1:nn, nn, replace = T)
        }
        while (min(table(y[id])) < 2) {
          id <- sample(1:nn, nn, replace = T)
        }
        if (class(d) == "numeric") {
          series[b] <- ccp(y = y[id], d = d[id], method = method, 
                           k = k, ...)
        }
        else {
          series[b] <- ccp(y = y[id], d = d[id, ], method = method, 
                           k = k, ...)
        }
      }
    }
    if (balance == TRUE) {
      for (b in 1:B) {
        id <- unlist(caret::createResample(y, times = 1))
        if (class(d) == "numeric") {
          series[b] <- ccp(y = y[id], d = d[id], method = method, 
                           k = k, ...)
        }
        else {
          series[b] <- ccp(y = y[id], d = d[id, ], method = method, 
                           k = k, ...)
        }
      }
    }
    series.sort <- sort(series)
    return(list(value = ccp(y = y, d = d, method = method, 
                            k = k, ...), se = sd(series), interval = c(series.sort[ifelse(B * 
                                                                                            (0.5 - level/2) < 1, 1, B * (0.5 - level/2))], series.sort[B * 
                                                                                                                                                         (0.5 + level/2)])))
  }
  if (acc == "rsq") {
    if (balance == FALSE) {
      for (b in 1:B) {
        nn <- length(y)
        id <- sample(1:nn, nn, replace = T)
        while (length(unique(y[id])) < k) {
          id <- sample(1:nn, nn, replace = T)
        }
        while (min(table(y[id])) < 2) {
          id <- sample(1:nn, nn, replace = T)
        }
        if (class(d) == "numeric") {
          series[b] <- rsq(y = y[id], d = d[id], method = method, 
                           k = k, ...)
        }
        else {
          series[b] <- rsq(y = y[id], d = d[id, ], method = method, 
                           k = k, ...)
        }
      }
    }
    if (balance == TRUE) {
      for (b in 1:B) {
        id <- unlist(caret::createResample(y, times = 1))
        if (class(d) == "numeric") {
          series[b] <- rsq(y = y[id], d = d[id], method = method, 
                           k = k, ...)
        }
        else {
          series[b] <- rsq(y = y[id], d = d[id, ], method = method, 
                           k = k, ...)
        }
      }
    }
    series.sort <- sort(series)
    return(list(value = rsq(y = y, d = d, method = method, 
                            k = k, ...), se = sd(series), interval = c(series.sort[ifelse(B * 
                                                                                            (0.5 - level/2) < 1, 1, B * (0.5 - level/2))], series.sort[B * 
                                                                                                                                                         (0.5 + level/2)])))
  }
}
