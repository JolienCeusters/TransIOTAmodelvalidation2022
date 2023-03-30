#%%%%%%%%%%%%%%%#
# ROMA vs ADNEX #
#%%%%%%%%%%%%%%%#


## Packages
library(reshape2)
library(readxl)
library(data.table)
library(auRoc)
library(caret)
library(car)
library(boot)
library(plyr)
library(forestplot)
library(REMA)
library(mice)


#-------------------------#
##### 1. Load dataset #####
#-------------------------#

setwd("C:/Users/u0123496/OneDrive - KU Leuven/IOTA/ROMA vs ADNEX")
load("TransIOTA - data.RData")

ToExclude <- read_excel("c:/Users/u0123496/OneDrive - KU Leuven/Trans-IOTA/Final analyses/To exclude from TransIOTA - Proteins.xlsx")
IDsExclude <- ToExclude$`IOTA ID`

HistSpecial <- read_excel_allsheets("c:/Users/u0123496/OneDrive - KU Leuven/Trans-IOTA/Final analyses/Histologies TransIOTA - Chiara_CL.xlsx")
HistSpecialLong <- rbindlist(HistSpecial, fill = TRUE)
IDsHis <- HistSpecialLong$`__metadata__.__caseId__`

##### 1.1 Check/Apply exclusion criteria #####

# 1. Denial or withdrawal of written informed consent
table(Trans$iota7_basic.scan_IC) ## No patients withdrew consent

# 2. Patients younger than 18 years
summary(Trans$iota7_basic.age) ## Lowest age in dataset is 16
IDminor <- Trans$`__metadata__.__caseId__`[Trans$iota7_basic.age < 18] ## 4 minors

# 3. Patients operated later than 120 days after recruitment
IntervalEventsTr <- Trans[as.Date(Trans$iota7_basic.surgery_date, "%d/%m/%Y") - as.Date(Trans$iota7_basic.scan_date, "%d/%m/%Y") > 120,
                          c("__metadata__.__userName__", "__metadata__.__realmId__", "__metadata__.__caseId__", "iota7_basic.month_of_birth", "iota7_basic.year_of_birth", "iota7_basic.ultrasound_examiner",
                            "iota7_basic.surgery_date", "iota7_basic.scan_date")]
IntervalEventsTr <- IntervalEventsTr[!is.na(IntervalEventsTr$iota7_basic.year_of_birth), ]
IntervalEventsTr ## 12 patients have a longer interval

# 4. Pregnancy
table(Trans$iota7_basic.pregnancy) # 6 patients are pregnant

# 5. Refusal of preoperative transvaginal ultrasonography and/or blood sample
table(Trans$iota_trans_local.def3)
table(Trans$iota_trans_local.serum)
IDnoBlood <- Trans$`__metadata__.__caseId__`[Trans$iota_trans_local.def2 != "true"]

# 6. Known simultaneous and/or previous malignancies within five years

# 7. simultaneous autoimmune disease and/or treatment with immune modulating drugs
table(Trans$iota_trans_local.personal_history_immunologicaldisease)
table(Trans$iota_trans_local.personal_history_immunologicaldisease_specify)
IDai <- Trans$`__metadata__.__caseId__`[Trans$iota_trans_local.personal_history_immunologicaldisease == 1]

# 8. Infectious serology (i.e. HIV, Hepatitis B, Hepatitis C)


## Apply exclusion criteria
TransIOTAanalyses <- TransIOTA[!TransIOTA$`Patient ID` %in% IDsExclude, ]


##### 1.2 Create variables #####

## Make variable ADNEX groups
table(TransIOTAanalyses$`FOR JOLIEN`)
table(TransIOTAanalyses$`Mass Outcome`[is.na(TransIOTAanalyses$`FOR JOLIEN`)])
table(TransIOTAanalyses$`Mass Outcome`[is.na(TransIOTAanalyses$ADNEXgroups)])

TransIOTAanalyses = ddply(TransIOTAanalyses, .(`Patient ID`),
                          function(x){
                            ADNEXgroups =
                              if(!is.na(x$`FOR JOLIEN`)){
                                if(x$`FOR JOLIEN` == "Benign"){
                                  "Benign"
                                } # Einde benign
                                else{
                                  if(x$`FOR JOLIEN` == "Borderline"){
                                    "Borderline"
                                  } # Borderline
                                  else{
                                    if(x$`FOR JOLIEN` == "Metastasis"){
                                      "Metastatic"
                                    } # Metastatic
                                    else{
                                      if(x$`FOR JOLIEN` == "Stage I"){
                                        "Stage I invasive"
                                      } # Stage I
                                      else{
                                        if(x$`FOR JOLIEN` == "Stage II-IV"){
                                          "Stage II-IV invasive"
                                        } # Stage II - IV
                                        else{
                                          "Unknown"
                                        } # Einde not Stage II - IV
                                      } # Einde not Stage I
                                    } # Einde not metastatic
                                  } # Einde not borderline
                                } # Einde not benign
                              } # Einde if For Jolien is not missing
                            else{
                              if(x$`Patient ID` %in% IDsHis){
                                HistSpecialLong$Group[HistSpecialLong$`__metadata__.__caseId__` == x$`Patient ID`]
                              } # If in list Chiara
                              else{
                                if(x$`Mass Outcome` == "benign"){
                                  "Benign"
                                } # Benign
                                else{
                                  if(x$`Mass Outcome` == "borderline"){
                                    "Borderline"
                                  } # Borderline
                                  else{
                                    if(x$`Mass Outcome` == "infectious_acute_chronic"){
                                      "Benign"
                                    } # Infectious acute chronic
                                    else{
                                      if(x$`Mass Outcome` == "malignant"){
                                        if(x$`FIGO stage.x` == 1){
                                          "Stage I invasive"
                                        } # FIGO stage I
                                        else{
                                          if(x$`FIGO stage.x` == 2 | x$`FIGO stage.x` == 3 | x$`FIGO stage.x` == 4){
                                            "Stage II-IV invasive"
                                          } # FIGO stage II or III or IV
                                          else{
                                            "Invasive FIGO stage unknown"
                                          } # Not FIGO stage II or III or IV
                                        } # Not FIGO stage I
                                      } # Malignant
                                      else{
                                        if(x$`Mass Outcome` == "metastatic"){
                                          "Metastatic"
                                        } # Metastatic
                                        else{
                                          if(x$`Mass Outcome` == "rare_benign"){
                                            "Benign"
                                          } # Rare benign
                                          else{
                                            "Unknown"
                                          } # Not rare benign
                                        } # Not metastatic
                                      } # Not Malignant
                                    } # Not infectious acute chronic
                                  } # Not borderline
                                } # Not benign
                              } # Not in list Chiara
                            } # For Jolien is missing
                            
                            x$ADNEXgroups = ADNEXgroups
                            return(x)
                          }# Einde functie
)
table(TransIOTAanalyses$ADNEXgroups)

## Make variable Binary outcome
TransIOTAanalyses$OutcomeBin <- ifelse(TransIOTAanalyses$ADNEXgroups == "Benign", 0, 1)

## Make variable for center
TransIOTAanalyses$CenterAUC <- ifelse(grepl("Genk", TransIOTAanalyses$center), "Genk",
                                      ifelse(grepl("Prague", TransIOTAanalyses$center), "Prague",
                                      ifelse(grepl("Milan", TransIOTAanalyses$center), "Milan",
                                      ifelse(grepl("Rome", TransIOTAanalyses$center), "Rome",
                                      ifelse(grepl("London", TransIOTAanalyses$center), "London", "Leuven")))))

## Make variable Oncocenter
TransIOTAanalyses$oncocenter <- ifelse(grepl("Genk", TransIOTAanalyses$center), 0, 1) 

## Set working directory
setwd("c:/Users/u0123496/OneDrive - KU Leuven/IOTA/ROMA vs ADNEX/Results")


#-----------------------------------#
##### 2. Complete case analysis #####
#-----------------------------------#

DataCC <- TransIOTAanalyses[!is.na(TransIOTAanalyses$`CA125 (kU/L)`) & TransIOTAanalyses$`CA125 (kU/L)` != "/",]

DataCC$`CA125 (kU/L)` <- as.numeric(DataCC$`CA125 (kU/L)`)
DataCC$`HE4 (pmol/L)` <- as.numeric(DataCC$`HE4 (pmol/L)`)

##### 2.1 Descriptive statistics #####

## Age
median(DataCC$`Patient age`) # 53
quantile(DataCC$`Patient age`, probs = 0.25) # 42
quantile(DataCC$`Patient age`, probs = 0.75) # 64
min(DataCC$`Patient age`) # 18
max(DataCC$`Patient age`) # 88

## Menopausal status
table(DataCC$Postmenopausal2) # No: 387; Yes: 507

## Presence of solid components
table(DataCC$solidbin) # 0: 284; 1: 610

## Observed CA125
median(DataCC$`CA125 (kU/L)`) # 27.93
quantile(DataCC$`CA125 (kU/L)`, probs = 0.25) # 13
quantile(DataCC$`CA125 (kU/L)`, probs = 0.75) # 135
min(DataCC$`CA125 (kU/L)`) # 3
max(DataCC$`CA125 (kU/L)`) # 24137

mean(DataCC$`CA125 (kU/L)`)

## Observed HE4
median(DataCC$`HE4 (pmol/L)`) # 64
quantile(DataCC$`HE4 (pmol/L)`, probs = 0.25) # 49
quantile(DataCC$`HE4 (pmol/L)`, probs = 0.75) # 127.825
min(DataCC$`HE4 (pmol/L)`) # 26
max(DataCC$`HE4 (pmol/L)`) # 19005

## Maximum diameter of lesion (mm)
median(DataCC$`Lesion largest diameter`) # 74
quantile(DataCC$`Lesion largest diameter`, probs = 0.25) # 48.25
quantile(DataCC$`Lesion largest diameter`, probs = 0.75) # 117
min(DataCC$`Lesion largest diameter`) # 7
max(DataCC$`Lesion largest diameter`) # 459

## Proportion solid tissue
median(DataCC$propsol) # 0.38
quantile(DataCC$propsol, probs = 0.25) # 0
quantile(DataCC$propsol, probs = 0.75) # 1
min(DataCC$propsol) # 0
max(DataCC$propsol) # 1

## Number of papillary projections
table(DataCC$papnr)

## Bilateral masses
table(DataCC$bilatbin) # 0: 661; 1: 233

## More than 10 cyst locules
table(DataCC$loc10) # 0: 759; 1: 135

## Acoustic shadows
table(DataCC$shadowsbin) # 0: 568; 1: 326

table(DataCC$shadowsbin, DataCC$CenterAUC)
table(DataCC$shadowsbin, DataCC$ADNEXgroups)

## Ascites
table(DataCC$ascitesbin) # 0: 761; 1: 133

## Multilocular cysts
table(DataCC$multibin) # 0: 506; 1: 388

## Metastases 
table(DataCC$metasbin) #0: 765; 1: 129

## Ultrasound examiner's subjective impression
table(DataCC$Certainty)

## Outcome
table(DataCC$ADNEXgroups)

## Histology
View(DataCC[, c("Histology.x", "Histology.y", "ADNEXgroups")])
table(DataCC$Histology.x)
table(DataCC$Histology.x, DataCC$ADNEXgroups)


##### 2.2 Validation #####

##### 2.2.1 Model predictions #####

source("C:/Users/u0123496/OneDrive - KU Leuven/IOTA/IOTA5/Scripts/FunctionsAllPredModels2.R")

##### ADNEX #####

## With CA125
DataCC = ADNEX(Age = `Patient age`, Oncocenter = oncocenter, lesdmax = `Lesion largest diameter`, soldmaxorig = soldmax, loc10 = loc10, Papnr = papnr, Shadows = shadowsbin, Ascites = ascitesbin, ca125 = `CA125 (kU/L)`, wica = T, woca = F, data = DataCC)

## Histogram: back-to-back
ADNEXw.hist <- ggplot(DataCC, aes(x=pmalw)) +
               geom_histogram(data=subset(DataCC,OutcomeBin == '0'),aes(y = ..count.., fill = "Benign"), binwidth = 0.01, alpha=.75) +
               geom_histogram(data=subset(DataCC,OutcomeBin == '1'),aes(y = -..count.., fill = "Malignant"), binwidth = 0.01, alpha=.75) +
               theme_minimal() +
               labs(x = "Estimated risk",
                    y = "Frequency",
                    fill = "Outcome") +
               scale_fill_manual(name = "Outcome", values = c(Benign = "#00BFC4", Malignant = "#F8766D")) +
               scale_y_continuous(breaks = c(-50, -25, 0, 25, 50, 75), labels = c(50, 25, 0, 25, 50, 75))
ADNEXw.hist
# 1000 x 400

## Without CA125
DataCC = ADNEX(Age = `Patient age`, Oncocenter = oncocenter, lesdmax = `Lesion largest diameter`, soldmaxorig = soldmax, loc10 = loc10, Papnr = papnr, Shadows = shadowsbin, Ascites = ascitesbin, ca125 = `CA125 (kU/L)`, wica = F, woca = T, data = DataCC)

## Histogram: back-to-back
ADNEXwo.hist <- ggplot(DataCC, aes(x=pmalwo)) +
                geom_histogram(data=subset(DataCC,OutcomeBin == '0'),aes(y = ..count.., fill = "Benign"), binwidth = 0.01, alpha=.75) +
                geom_histogram(data=subset(DataCC,OutcomeBin == '1'),aes(y = -..count.., fill = "Malignant"), binwidth = 0.01, alpha=.75) +
                theme_minimal() +
                labs(x = "Estimated risk",
                     y = "Frequency",
                     fill = "Outcome") +
                scale_fill_manual(name = "Outcome", values = c(Benign = "#00BFC4", Malignant = "#F8766D")) +
                scale_y_continuous(limits = c(-25, 70), breaks = c(-25, 0, 25, 50, 75), labels = c(25, 0, 25, 50, 75))
ADNEXwo.hist
# 1000 x 400

##### ROMA #####

## HE4 expressed in pmol/L
## CA125 expressed in U/mL (= kU/L)

DataCC <- ddply(.data = DataCC, .(`Patient ID`),
                function(x){
                  x$ROMA =
                    if(x$Postmenopausal2 == "yes"){
                      z = -8.09 + (1.04 * log(x$`HE4 (pmol/L)`)) + (0.732 * log(x$`CA125 (kU/L)`))
                      p = 1/(1 + exp(-z))
                      p
                    }else{
                      z = -12.0 + (2.38 * log(x$`HE4 (pmol/L)`)) + (0.0626 * log(x$`CA125 (kU/L)`))
                      p = 1/(1 + exp(-z))
                      p
                    }
                  return(x)
                })

## Histogram: back-to-back
ROMA.hist <- ggplot(DataCC, aes(x=ROMA)) +
             geom_histogram(data=subset(DataCC,OutcomeBin == '0'),aes(y = ..count.., fill = "Benign"), binwidth = 0.01, alpha=.75) +
             geom_histogram(data=subset(DataCC,OutcomeBin == '1'),aes(y = -..count.., fill = "Malignant"), binwidth = 0.01, alpha=.75) +
             theme_minimal() +
             labs(x = "Estimated risk",
                  y = "Frequency",
                  fill = "Outcome") +
             scale_fill_manual(name = "Outcome", values = c(Benign = "#00BFC4", Malignant = "#F8766D")) +
             scale_y_continuous(breaks = c(-25, 0, 25, 50), labels = c(25, 0, 25, 50))
ROMA.hist
# 1000 x 400


##### 2.2.2 Discrimination between benign and malignant tumors #####

##### ADNEX with CA125 #####

## AUC
AUC.ADNEXw <- AUC.IOTA(pred = pmalw, outcome = OutcomeBin, center = CenterAUC, data = DataCC, method.MA = "BAYES", titleGraph = "AUC per center")

AUC.ADNEXw$Plot[9, "RRauc"] <- "   "
AUC.ADNEXw$Plot[10, "RRauc"] <- "   "
AUC.ADNEXw$Plot[11, "RRcenter"] <- "AUC (95% CI)"
AUC.ADNEXw$Plot[12, "RRauc"] <- "        (0.83 to 0.97)"
AUC.ADNEXw$Plot[12, "RRprev"] <- "   "

forestplot(AUC.ADNEXw$Plot,
           align = c("l", "c", "c"),
           mean = AUC.ADNEXw$dataPlot$AUC,
           lower = AUC.ADNEXw$dataPlot$LL,
           upper = AUC.ADNEXw$dataPlot$UL,
           is.summary = c(TRUE, FALSE, rep(FALSE, length(AUC.ADNEXw$IncludedCenters)), FALSE, FALSE, TRUE, TRUE),
           xlab = "AUC", 
           xticks = c(0.70, 0.80, 0.9, 1),
           clip = c(0.71, 1),
           txt_gp = fpTxtGp(summary = gpar(cex = 1.55), xlab = gpar(cex = 1.5), label = gpar(cex = 1.5),
                            ticks = gpar(cex = 1.5), title = gpar(cex = 1.75)),
           graphwidth = unit(6, "cm"),
           col = fpColors(box = "gray", lines = "black", summary = "black"),
           fn.ci_sum = c(as.list(rep("fpDrawSummaryCI", length(AUC.ADNEXw$IncludedCenters) + 4)),
                         fpDrawSummaryCI,
                         fpDrawBarCI),
           zero = AUC.ADNEXw$Performance$AUC)
# 900 x 450

## Sensitivity and specificity
SensSpec(pred = pmalw, outcome = OutcomeBin, threshold = 0.01,  center = CenterAUC, data = DataCC)$OverallPer
SensSpec(pred = pmalw, outcome = OutcomeBin, threshold = 0.03,  center = CenterAUC, data = DataCC)$OverallPer
SensSpec(pred = pmalw, outcome = OutcomeBin, threshold = 0.05,  center = CenterAUC, data = DataCC)$OverallPer
SensSpec(pred = pmalw, outcome = OutcomeBin, threshold = 0.10,  center = CenterAUC, data = DataCC)$OverallPer
SensSpec(pred = pmalw, outcome = OutcomeBin, threshold = 0.15,  center = CenterAUC, data = DataCC)$OverallPer
SensSpec(pred = pmalw, outcome = OutcomeBin, threshold = 0.20,  center = CenterAUC, data = DataCC)$OverallPer
SensSpec(pred = pmalw, outcome = OutcomeBin, threshold = 0.25,  center = CenterAUC, data = DataCC)$OverallPer
SensSpec(pred = pmalw, outcome = OutcomeBin, threshold = 0.30,  center = CenterAUC, data = DataCC)$OverallPer
SensSpec(pred = pmalw, outcome = OutcomeBin, threshold = 0.40,  center = CenterAUC, data = DataCC)$OverallPer
SensSpec(pred = pmalw, outcome = OutcomeBin, threshold = 0.50,  center = CenterAUC, data = DataCC)$OverallPer

## Sensitivity at fixed specificity
FixSpec(pred = pmalw, outcome = OutcomeBin, Specificity = 0.90,  center = CenterAUC, data = DataCC)$OverallPer

## Specificity at fixed sensitivity
FixSens(pred = pmalw, outcome = OutcomeBin, Sensitivity = 0.90,  center = CenterAUC, data = DataCC)$OverallPer

##### ADNEX without CA125 #####

## AUC
AUC.ADNEXwo <- AUC.IOTA(pred = pmalwo, outcome = OutcomeBin, center = CenterAUC, data = DataCC, method.MA = "BAYES", titleGraph = "AUC per center")

AUC.ADNEXwo$Plot[9, "RRauc"] <- "   "
AUC.ADNEXwo$Plot[10, "RRauc"] <- "   "
AUC.ADNEXwo$Plot[11, "RRcenter"] <- "AUC (95% CI)"
AUC.ADNEXwo$Plot[12, "RRauc"] <- "        (0.70 to 0.97)"
AUC.ADNEXwo$Plot[12, "RRprev"] <- "   "

forestplot(AUC.ADNEXwo$Plot,
           align = c("l", "c", "c"),
           mean = AUC.ADNEXwo$dataPlot$AUC,
           lower = AUC.ADNEXwo$dataPlot$LL,
           upper = AUC.ADNEXwo$dataPlot$UL,
           is.summary = c(TRUE, FALSE, rep(FALSE, length(AUC.ADNEXwo$IncludedCenters)), FALSE, FALSE, TRUE, TRUE),
           xlab = "AUC", 
           xticks = c(0.70, 0.80, 0.9, 1),
           clip = c(0.71, 1),
           txt_gp = fpTxtGp(summary = gpar(cex = 1.55), xlab = gpar(cex = 1.5), label = gpar(cex = 1.5),
                            ticks = gpar(cex = 1.5), title = gpar(cex = 1.75)),
           graphwidth = unit(6, "cm"),
           col = fpColors(box = "gray", lines = "black", summary = "black"),
           fn.ci_sum = c(as.list(rep("fpDrawSummaryCI", length(AUC.ADNEXwo$IncludedCenters) + 4)),
                         fpDrawSummaryCI,
                         fpDrawBarCI),
           zero = AUC.ADNEXwo$Performance$AUC)
# 900 x 450

## Sensitivity and specificity
SensSpec(pred = pmalwo, outcome = OutcomeBin, threshold = 0.01, center = CenterAUC, data = DataCC)$OverallPer
SensSpec(pred = pmalwo, outcome = OutcomeBin, threshold = 0.03, center = CenterAUC, data = DataCC)$OverallPer
SensSpec(pred = pmalwo, outcome = OutcomeBin, threshold = 0.05, center = CenterAUC, data = DataCC)$OverallPer
SensSpec(pred = pmalwo, outcome = OutcomeBin, threshold = 0.10, center = CenterAUC, data = DataCC)$OverallPer
SensSpec(pred = pmalwo, outcome = OutcomeBin, threshold = 0.15, center = CenterAUC, data = DataCC)$OverallPer
SensSpec(pred = pmalwo, outcome = OutcomeBin, threshold = 0.20, center = CenterAUC, data = DataCC)$OverallPer
SensSpec(pred = pmalwo, outcome = OutcomeBin, threshold = 0.25, center = CenterAUC, data = DataCC)$OverallPer
SensSpec(pred = pmalwo, outcome = OutcomeBin, threshold = 0.30, center = CenterAUC, data = DataCC)$OverallPer
SensSpec(pred = pmalwo, outcome = OutcomeBin, threshold = 0.40, center = CenterAUC, data = DataCC)$OverallPer
SensSpec(pred = pmalwo, outcome = OutcomeBin, threshold = 0.50, center = CenterAUC, data = DataCC)$OverallPer

## Sensitivity at fixed specificity
FixSpec(pred = pmalwo, outcome = OutcomeBin, Specificity = 0.90, center = CenterAUC, data = DataCC)$OverallPer

## Specificity at fixed sensitivity
FixSens(pred = pmalwo, outcome = OutcomeBin, Sensitivity = 0.90, center = CenterAUC, data = DataCC)$OverallPer

##### ROMA #####

## AUC
AUC.ROMA <- AUC.IOTA(pred = ROMA, outcome = OutcomeBin, center = CenterAUC, data = DataCC, method.MA = "BAYES", titleGraph = "AUC per center")

AUC.ROMA$Plot[9, "RRauc"] <- "   "
AUC.ROMA$Plot[10, "RRauc"] <- "   "
AUC.ROMA$Plot[11, "RRcenter"] <- "AUC (95% CI)"
AUC.ROMA$Plot[12, "RRauc"] <- "        (0.73 to 0.92)"
AUC.ROMA$Plot[12, "RRprev"] <- "   "

forestplot(AUC.ROMA$Plot,
           align = c("l", "c", "c"),
           mean = AUC.ROMA$dataPlot$AUC,
           lower = AUC.ROMA$dataPlot$LL,
           upper = AUC.ROMA$dataPlot$UL,
           is.summary = c(TRUE, FALSE, rep(FALSE, length(AUC.ROMA$IncludedCenters)), FALSE, FALSE, TRUE, TRUE),
           xlab = "AUC", 
           xticks = c(0.70, 0.80, 0.9, 1),
           clip = c(0.71, 1),
           txt_gp = fpTxtGp(summary = gpar(cex = 1.55), xlab = gpar(cex = 1.5), label = gpar(cex = 1.5),
                            ticks = gpar(cex = 1.5), title = gpar(cex = 1.75)),
           graphwidth = unit(6, "cm"),
           col = fpColors(box = "gray", lines = "black", summary = "black"),
           fn.ci_sum = c(as.list(rep("fpDrawSummaryCI", length(AUC.ROMA$IncludedCenters) + 4)),
                         fpDrawSummaryCI,
                         fpDrawBarCI),
           zero = AUC.ROMA$Performance$AUC)
# 900 x 450

## Sensitivity and specificity
SensSpec(pred = ROMA, outcome = OutcomeBin, threshold = 0.01, center = CenterAUC, data = DataCC)$OverallPer
SensSpec(pred = ROMA, outcome = OutcomeBin, threshold = 0.03, center = CenterAUC, data = DataCC)$OverallPer
SensSpec(pred = ROMA, outcome = OutcomeBin, threshold = 0.05, center = CenterAUC, data = DataCC)$OverallPer
SensSpec(pred = ROMA, outcome = OutcomeBin, threshold = 0.10, center = CenterAUC, data = DataCC)$OverallPer
SensSpec(pred = ROMA, outcome = OutcomeBin, threshold = 0.15, center = CenterAUC, data = DataCC)$OverallPer
SensSpec(pred = ROMA, outcome = OutcomeBin, threshold = 0.20, center = CenterAUC, data = DataCC)$OverallPer
SensSpec(pred = ROMA, outcome = OutcomeBin, threshold = 0.25, center = CenterAUC, data = DataCC)$OverallPer
SensSpec(pred = ROMA, outcome = OutcomeBin, threshold = 0.30, center = CenterAUC, data = DataCC)$OverallPer
SensSpec(pred = ROMA, outcome = OutcomeBin, threshold = 0.40, center = CenterAUC, data = DataCC)$OverallPer
SensSpec(pred = ROMA, outcome = OutcomeBin, threshold = 0.50, center = CenterAUC, data = DataCC)$OverallPer

## ROMA cut-offs

# Postmenopausal
SensSpec(pred = ROMA, outcome = OutcomeBin, threshold = 0.299, center = CenterAUC, data = DataCC[DataCC$Postmenopausal2 == "yes",])$OverallPer

ROMApost <- SensSpec(pred = ROMA, outcome = OutcomeBin, threshold = 0.299, center = CenterAUC, data = DataCC[DataCC$Postmenopausal2 == "yes",])$CenterPer

# Premenopausal
SensSpec(pred = ROMA, outcome = OutcomeBin, threshold = 0.114, center = CenterAUC, data = DataCC[DataCC$Postmenopausal2 == "no",])$OverallPer

ROMApre <- SensSpec(pred = ROMA, outcome = OutcomeBin, threshold = 0.114, center = CenterAUC, data = DataCC[DataCC$Postmenopausal2 == "no",])$CenterPer

# Combined
SensSpecROMA <- matrix(nrow = 6, ncol = 5)
colnames(SensSpecROMA) <- c("Center", "TP", "TN", "FP", "FN")
SensSpecROMA <- data.frame(SensSpecROMA)

SensSpecROMA$Center <- c("Leuven", "Prague", "London", "Genk", "Milan", "Rome")

SensSpecROMA$TN     <- ROMApost$TN + ROMApre$TN
SensSpecROMA$TP     <- ROMApost$TP + ROMApre$TP
SensSpecROMA$FN     <- ROMApost$FN + ROMApre$FN
SensSpecROMA$FP     <- ROMApost$FP + ROMApre$FP

SensSpecROMA$Sensitivity  <- SensSpecROMA$TP / (SensSpecROMA$TP + SensSpecROMA$FN)
SensSpecROMA$se_sens      <- sqrt((SensSpecROMA$Sensitivity * (1 - SensSpecROMA$Sensitivity))/(SensSpecROMA$TP + SensSpecROMA$FN))
SensSpecROMA$Specificity  <- SensSpecROMA$TN / (SensSpecROMA$TN + SensSpecROMA$FP)
SensSpecROMA$se_spec      <- sqrt((SensSpecROMA$Specificity * (1 - SensSpecROMA$Specificity))/(SensSpecROMA$TN + SensSpecROMA$FP))

# Logit transformation
SensSpecROMA$logit.sens     <- logit(SensSpecROMA$Sensitivity)
SensSpecROMA$logit.se.sens  <- sqrt(1/SensSpecROMA$TP + 1/SensSpecROMA$FN)
SensSpecROMA$logit.LL.sens  <- SensSpecROMA$logit.sens - 1.96*SensSpecROMA$logit.se.sens
SensSpecROMA$logit.UL.sens  <- SensSpecROMA$logit.sens + 1.96*SensSpecROMA$logit.se.sens

SensSpecROMA$logit.spec     <- logit(SensSpecROMA$Specificity)
SensSpecROMA$logit.se.spec  <- sqrt(1/SensSpecROMA$TN + 1/SensSpecROMA$FP)
SensSpecROMA$logit.LL.spec  <- SensSpecROMA$logit.spec - 1.96*SensSpecROMA$logit.se.spec
SensSpecROMA$logit.UL.spec  <- SensSpecROMA$logit.spec + 1.96*SensSpecROMA$logit.se.spec

Sensoverall <- matrix(nrow = 1, ncol = 7)
colnames(Sensoverall) <- c('Center', 'Sens', 'LL.sens', 'UL.sens', 'Spec', 'LL.spec', 'UL.spec')
Sensoverall <- data.frame(Sensoverall)
Sensoverall$Center <- "Overall"

# Meta-analysis for overall estimate: Bivariate random-effects model
p1 <- 2*6
p2 <- 4*6
Sensbi <- SensSpecROMA[, c("Center", "logit.sens", "logit.spec", "logit.se.sens", "logit.se.spec")]
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

## Sensitivity at fixed specificity
FixSpec(pred = ROMA, outcome = OutcomeBin, Specificity = 0.90, center = CenterAUC, data = DataCC)$OverallPer

## Specificity at fixed sensitivity
FixSens(pred = ROMA, outcome = OutcomeBin, Sensitivity = 0.90, center = CenterAUC, data = DataCC)$OverallPer

##### Summary plot for AUC #####

NA.forest <- AUC.ROMA
NA.forest <- NA

Summary.AUC <- rbind(NA.forest, AUC.ROMA$Performance[1,], AUC.ADNEXwo$Performance[1,], AUC.ADNEXw$Performance[1,])
Summary.AUC$Model <- c('', 'ROMA', 'ADNEX without CA125', 'ADNEX with CA125')

Summary.PI <- rbind(NA.forest, AUC.ROMA$Performance[2,], AUC.ADNEXwo$Performance[2,], AUC.ADNEXw$Performance[2,])
Summary.PI$Model <- c('', 'ROMA', 'ADNEX without CA125', 'ADNEX with CA125')

tabletext <- cbind(
  c('Model', 'ROMA', 'ADNEX without CA125', 'ADNEX with CA125'),
  c('AUC (95% CI)', 
    paste(format(round(AUC.ROMA$Performance$AUC[1], 2), nsmall = 2), " (", format(round(AUC.ROMA$Performance$LL[1], 2), nsmall = 2), "; ", format(round(AUC.ROMA$Performance$UL[1], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(AUC.ADNEXwo$Performance$AUC[1], 2), nsmall = 2), " (", format(round(AUC.ADNEXwo$Performance$LL[1], 2), nsmall = 2), "; ", format(round(AUC.ADNEXwo$Performance$UL[1], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(AUC.ADNEXw$Performance$AUC[1], 2), nsmall = 2), " (", format(round(AUC.ADNEXw$Performance$LL[1], 2), nsmall = 2), "; ", format(round(AUC.ADNEXw$Performance$UL[1], 2), nsmall = 2), ")", sep = "")),
  c('95% PI', 
    paste("(", format(round(AUC.ROMA$Performance$LL[2], 2), nsmall = 2), "; ", format(round(AUC.ROMA$Performance$UL[2], 2), nsmall = 2), ")", sep = ""),
    paste("(", format(round(AUC.ADNEXwo$Performance$LL[2], 2), nsmall = 2), "; ", format(round(AUC.ADNEXwo$Performance$UL[2], 2), nsmall = 2), ")", sep = ""),
    paste("(", format(round(AUC.ADNEXw$Performance$LL[2], 2), nsmall = 2), "; ", format(round(AUC.ADNEXw$Performance$UL[2], 2), nsmall = 2), ")", sep = ""))
)

forestplot(labeltext = tabletext,
           title = "",
           mean = round(Summary.AUC$AUC, 3),
           lower = round(Summary.AUC$LL, 3),
           upper = round(Summary.AUC$UL, 3),
           is.summary = c(FALSE, TRUE, TRUE, TRUE),
           xlab = "AUC (95% CI)",
           boxsize = .5,
           txt_gp = fpTxtGp(summary = gpar(cex = 1.55, fontface = "plain"), xlab = gpar(cex = 1.5, fontface = "plain"), label = gpar(cex = 1.5, fontface = "bold"),
                            ticks = gpar(cex = 1.5, fontface = "plain"), title = gpar(cex = 1.75)),
           graphwidth = unit(9, "cm"),
           graph.pos = 3,
           xticks = c(0.80, 0.85, 0.9, 0.95, 1), xlog = TRUE, clip = c(0.80, 1))
# 1050 x 400

##### Difference in AUC #####

## Calculate the difference per center
centers <- unique(DataCC$CenterAUC)
NRcent <- length(centers)

AUCcenter <- matrix(ncol = 3, nrow = length(centers))
colnames(AUCcenter) <- c('Center', 'SampleSize', 'Prevalence')
AUCcenter <- data.frame(AUCcenter)
AUCcenter$Center <- centers

for(i in seq_along(centers)){
  
  AUCcenter[i, 2]   <- nrow(DataCC[DataCC$CenterAUC == centers[i],])
  AUCcenter[i, 3]   <- round(nrow(DataCC[DataCC$OutcomeBin == 1 & DataCC$CenterAUC == centers[i],])/nrow(DataCC[DataCC$CenterAUC == centers[i],])*100)
  
  # ROMA
  roc1 <- roc(DataCC$OutcomeBin[DataCC$CenterAUC == centers[i]], DataCC$ROMA[DataCC$CenterAUC == centers[i]])
  # ADNEXwo
  roc2 <- roc(DataCC$OutcomeBin[DataCC$CenterAUC == centers[i]], DataCC$pmalwo[DataCC$CenterAUC == centers[i]])
  # ADNEXw
  roc3 <- roc(DataCC$OutcomeBin[DataCC$CenterAUC == centers[i]], DataCC$pmalw[DataCC$CenterAUC == centers[i]])
  
  v1 <- roc.test(roc2, roc1, paired = TRUE, method = "delong") # ADNEXwo vs ROMA
  v2 <- roc.test(roc3, roc1, paired = TRUE, method = "delong") # ADNEXw vs ROMA
  
  AUCcenter$AUC.ROMA[i] <-v1$estimate[2]
  AUCcenter$AUC.ADNEXwo[i] <- v1$estimate[1]
  AUCcenter$AUC.ADNEXw[i] <- v2$estimate[1]
  
  AUCcenter$ROMAvsADNEXwo[i] <- round(v1$estimate[1] - v1$estimate[2], 8)
  AUCcenter$ROMAvsADNEXwo.Var[i] <- ((AUCcenter$ROMAvsADNEXwo[i])/v1$statistic)^2
  AUCcenter$ROMAvsADNEXwo.p[i] <- v1$p.value
  
  AUCcenter$ROMAvsADNEXw[i] <- round(v2$estimate[1] - v2$estimate[2], 8)
  AUCcenter$ROMAvsADNEXw.Var[i] <- ((AUCcenter$ROMAvsADNEXw[i])/v2$statistic)^2
  AUCcenter$ROMAvsADNEXw.p[i] <- v2$p.value
}

## Meta-analysis for overall estimate
# ADNEXwo vs ROMA
fit.RE = rma.uni(AUCcenter$ROMAvsADNEXwo, sei = sqrt(AUCcenter$ROMAvsADNEXwo.Var), method = "REML")
coef(fit.RE)
fit.RE$ci.lb
fit.RE$ci.ub

# ADNEXw vs ROMA
fit.RE = rma.uni(AUCcenter$ROMAvsADNEXw, sei = sqrt(AUCcenter$ROMAvsADNEXw.Var), method = "REML")
coef(fit.RE)
fit.RE$ci.lb
fit.RE$ci.ub


##### 2.2.3 Calibration of the risk of malignancy #####

##### ADNEX with CA125 #####

Cal.ADNEXw <- Calibration(p = pmalw, y = OutcomeBin, center = CenterAUC, data = DataCC, CalibrLines = "overall", MethodCL = "Wald", nr.knots = 3) 

##### ADNEX without CA125 #####

Cal.ADNEXwo <- Calibration(p = pmalwo, y = OutcomeBin, center = CenterAUC, data = DataCC, CalibrLines = "overall", MethodCL = "Wald")

##### ROMA #####

Cal.ROMA <- Calibration(p = ROMA, y = OutcomeBin, center = CenterAUC, data = DataCC, CalibrLines = "overall", MethodCL = "Wald") 

##### Summary plot for Calibration #####

table <- matrix(ncol = 3, nrow = 3)
colnames(table) <- c('Model', 'Intercept (95% CI)', 'Slope (95% CI)')
table[, 1] <- c('ROMA', 'ADNEX without CA125', 'ADNEX with CA125')
table[1, 2:3] <- c(paste0(format(round(Cal.ROMA$Performance[1,1], 2), nsmall = 2), " (", format(round(Cal.ROMA$Performance[1,2], 2), nsmall = 2), " to ", format(round(Cal.ROMA$Performance[1,3], 2), nsmall = 2), ")"), paste0(format(round(Cal.ROMA$Performance[2,1], 2), nsmall = 2), " (", format(round(Cal.ROMA$Performance[2,2], 2), nsmall = 2), " to ", format(round(Cal.ROMA$Performance[2,3], 2), nsmall = 2), ")"))
table[3, 2:3] <- c(paste0(format(round(Cal.ADNEXw$Performance[1,1], 2), nsmall = 2), " (", format(round(Cal.ADNEXw$Performance[1,2], 2), nsmall = 2), " to ", format(round(Cal.ADNEXw$Performance[1,3], 2), nsmall = 2), ")"), paste0(format(round(Cal.ADNEXw$Performance[2,1], 2), nsmall = 2), " (", format(round(Cal.ADNEXw$Performance[2,2], 2), nsmall = 2), " to ", format(round(Cal.ADNEXw$Performance[2,3], 2), nsmall = 2), ")"))
table[2, 2:3] <- c(paste0(format(round(Cal.ADNEXwo$Performance[1,1], 2), nsmall = 2), " (", format(round(Cal.ADNEXwo$Performance[1,2], 2), nsmall = 2), " to ", format(round(Cal.ADNEXwo$Performance[1,3], 2), nsmall = 2), ")"), paste0(format(round(Cal.ADNEXwo$Performance[2,1], 2), nsmall = 2), " (", format(round(Cal.ADNEXwo$Performance[2,2], 2), nsmall = 2), " to ", format(round(Cal.ADNEXwo$Performance[2,3], 2), nsmall = 2), ")"))

library(plotrix)

## Make graph
x = seq(0, 1, by = 0.05)
y = seq(0, 1, by = 0.05)
tiff("Summary calibration - CC.tiff", width = 14, height = 14, units = "cm", res = 300)
plot(x, y, xlim = c(0,1), ylim = c(0,1), type = "l", col = "gray50", lwd = 2, lty = 2, 
     xlab = "Estimated risk", ylab = "Observed proportion",
     main = "", cex.lab = 1, cex.axis = 1, las = 1) 
lines(Cal.ROMA$Prediction$p, Cal.ROMA$Prediction$yhat, lwd = 2, col = "maroon")
lines(Cal.ADNEXwo$Prediction$p, Cal.ADNEXwo$Prediction$yhat, lwd = 2, col = "darkgreen")
lines(Cal.ADNEXw$Prediction$p, Cal.ADNEXw$Prediction$yhat, lwd = 2, col = "darkorange")
legend(x = -0.035, y = 1, legend = c( "Ideal", "ROMA", "ADNEX without CA125", "ADNEX with CA125"),
       col = c("gray50", "maroon", "darkgreen", "darkorange"), lty = c(2,1,1,1), lwd = 2, cex = 0.7, bty = "n", ncol = 1)
addtable2plot(x = 0.17, y = 0.02, table = table, display.colnames= TRUE, cex = 0.7)
dev.off()


##### 2.2.4 Clinical utility #####

##### ADNEX with CA125 #####

Dec.ADNEXw <- DataWinBugs(pred = pmalw, outcome = OutcomeBin, center = CenterAUC, data = DataCC)$Results

Dec.ADNEXwlong <- rbindlist(Dec.ADNEXw)
Dec.ADNEXwlong$prev <- (Dec.ADNEXwlong$TP + Dec.ADNEXwlong$FN) / Dec.ADNEXwlong$n
Dec.ADNEXwlong$TA <- Dec.ADNEXwlong$prev - (1 - Dec.ADNEXwlong$prev) * (Dec.ADNEXwlong$CutOff/(1 - Dec.ADNEXwlong$CutOff))

ggplot(data = Dec.ADNEXwlong, aes(x = CutOff, y = NB, colour = Center)) + 
  geom_line() +
  geom_line(aes(x = CutOff, y = TA))

##### ADNEX without CA125 #####

Dec.ADNEXwo <- DataWinBugs(pred = pmalwo, outcome = OutcomeBin, center = CenterAUC, data = DataCC)$Results

Dec.ADNEXwolong <- rbindlist(Dec.ADNEXwo)
Dec.ADNEXwolong$prev <- (Dec.ADNEXwolong$TP + Dec.ADNEXwolong$FN) / Dec.ADNEXwolong$n
Dec.ADNEXwolong$TA <- Dec.ADNEXwolong$prev - (1 - Dec.ADNEXwolong$prev) * (Dec.ADNEXwolong$CutOff/(1 - Dec.ADNEXwolong$CutOff))

ggplot(data = Dec.ADNEXwolong, aes(x = CutOff, y = NB, colour = Center)) + 
  geom_line() +
  geom_line(aes(x = CutOff, y = TA))

##### ROMA #####

Dec.ROMA <- DataWinBugs(pred = ROMA, outcome = OutcomeBin, center = CenterAUC, data = DataCC)$Results

Dec.ROMAlong <- rbindlist(Dec.ROMA)
Dec.ROMAlong$prev <- (Dec.ROMAlong$TP + Dec.ROMAlong$FN) / Dec.ROMAlong$n
Dec.ROMAlong$TA <- Dec.ROMAlong$prev - (1 - Dec.ROMAlong$prev) * (Dec.ROMAlong$CutOff/(1 - Dec.ROMAlong$CutOff))

ggplot(data = Dec.ROMAlong, aes(x = CutOff, y = NB, colour = Center)) + 
  geom_line() +
  geom_line(aes(x = CutOff, y = TA))

##### Summary plot for Clinical utility #####

Overview <- read_excel("Results WinBugs - CC.xlsx", sheet = "Overview")

plot(Overview$`Risk threshold`, Overview$`Treat all`, type = "l", col = "grey", xlab = "Risk threshold", ylab = "Net benefit", xlim = c(0.055, 0.5), ylim = c(-0.03, 0.55), lwd = 2, cex.lab = 1, cex.axis = 1, las = 1)
points(Overview$`Risk threshold`, Overview$`NB ROMA`, type = 'l', col = "maroon", lwd = 2)
points(Overview$`Risk threshold`, Overview$`NB ADNEX wo`, type = 'l', col = "darkgreen", lwd = 2)
points(Overview$`Risk threshold`, Overview$`NB ADNEXw`, type = 'l', col = "darkorange", lwd = 2)
abline(a = 0, b = 0, lty = 2, lwd = 2)
legend(x = 0.04, y = 0.55, legend = c("ROMA", "ADNEX without CA125", "ADNEX with CA125", "Treat all", "Treat none"), ncol = 3, #cex = 1.2, #bty = "n", #cex = 0.7, pt.cex = 1,
       col = c("maroon", "darkgreen", "darkorange", "grey", 'black'), lty = c(1,1,1,1,2), lwd = 2, cex = 0.7, bty = "n")


##### 2.2.5 Multinomial discrimination #####

#### ADNEX with CA125 ####

## 5 groups
DataCC$pbenw
DataCC$pborw
DataCC$pst1w
DataCC$pst2_4w
DataCC$pmetaw

## The 10 pairs
# Pair 1: benign vs borderline
DataCC$pp1 <- DataCC$pbenw / (DataCC$pbenw + DataCC$pborw)
round(auc.nonpara.mw(DataCC[DataCC$ADNEXgroups == "Benign", "pp1"], DataCC[DataCC$ADNEXgroups == "Borderline", "pp1"], method = "pepe"), 2)
# AUC = 0.86

# Pair 2: benign vs stage I
DataCC$pp2 <- DataCC$pbenw / (DataCC$pbenw + DataCC$pst1w)
round(auc.nonpara.mw(DataCC[DataCC$ADNEXgroups == "Benign", "pp2"], DataCC[DataCC$ADNEXgroups == "Stage I invasive", "pp2"], method = "pepe"), 2)
# AUC = 0.93

# Pair 3: benign vs stage II-IV
DataCC$pp3 <- DataCC$pbenw / (DataCC$pbenw + DataCC$pst2_4w)
round(auc.nonpara.mw(DataCC[DataCC$ADNEXgroups == "Benign", "pp3"], DataCC[DataCC$ADNEXgroups == "Stage II-IV invasive", "pp3"], method = "pepe"), 2)
# AUC = 0.98

# Pair 4: benign vs metastatic
DataCC$pp4 <- DataCC$pbenw / (DataCC$pbenw + DataCC$pmetaw)
round(auc.nonpara.mw(DataCC[DataCC$ADNEXgroups == "Benign", "pp4"], DataCC[DataCC$ADNEXgroups == "Metastatic", "pp4"], method = "pepe"), 2)
# AUC = 0.93

# Pair 5: borderline vs stage I
DataCC$pp5 <- DataCC$pborw / (DataCC$pborw + DataCC$pst1w)
DataCC$binbo <- ifelse(DataCC$ADNEXgroups == "Borderline", 0, 1)
round(auc.nonpara.mw(DataCC[DataCC$ADNEXgroups == "Borderline", "pp5"], DataCC[DataCC$ADNEXgroups == "Stage I invasive", "pp5"], method = "pepe"), 2)
# AUC = 0.69

# Pair 6: borderline vs stage II-IV
DataCC$pp6 <- DataCC$pborw / (DataCC$pborw + DataCC$pst2_4w)
round(auc.nonpara.mw(DataCC[DataCC$ADNEXgroups == "Borderline", "pp6"], DataCC[DataCC$ADNEXgroups == "Stage II-IV invasive", "pp6"], method = "pepe"), 2)
# AUC = 0.92

# Pair 7: borderline vs metastatic
DataCC$pp7 <- DataCC$pborw / (DataCC$pborw + DataCC$pmetaw)
round(auc.nonpara.mw(DataCC[DataCC$ADNEXgroups == "Borderline", "pp7"], DataCC[DataCC$ADNEXgroups == "Metastatic", "pp7"], method = "pepe"), 2)
# AUC = 0.86

# Pair 8: stage I vs stage II-IV
DataCC$pp8 <- DataCC$pst1w / (DataCC$pst1w + DataCC$pst2_4w)
DataCC$bini <- ifelse(DataCC$ADNEXgroups == "Stage I invasive", 0, 1)
round(auc.nonpara.mw(DataCC[DataCC$ADNEXgroups == "Stage I invasive", "pp8"], DataCC[DataCC$ADNEXgroups == "Stage II-IV invasive", "pp8"], method = "pepe"), 2)
# AUC = 0.86

# Pair 9: stage I vs metastatic
DataCC$pp9 <- DataCC$pst1w / (DataCC$pst1w + DataCC$pmetaw)
round(auc.nonpara.mw(DataCC[DataCC$ADNEXgroups == "Stage I invasive", "pp9"], DataCC[DataCC$ADNEXgroups == "Metastatic", "pp9"], method = "pepe"), 2)
# AUC = 0.73

# Pair 10: stage II-IV vs metastatic
DataCC$pp10 <- DataCC$pst2_4w / (DataCC$pst2_4w + DataCC$pmetaw)
DataCC$bin2 <- ifelse(DataCC$ADNEXgroups == "Stage II-IV invasive", 0, 1)
round(auc.nonpara.mw(DataCC[DataCC$ADNEXgroups == "Stage II-IV invasive", "pp10"], DataCC[DataCC$ADNEXgroups == "Metastatic", "pp10"], method = "pepe"), 2)
# AUC = 0.85

# Polytomous Discrimination Index (PDI)
library(mcca)

DataCC$Groups <- ifelse(DataCC$ADNEXgroups == "Benign", '1',
                        ifelse(DataCC$ADNEXgroups == "Borderline", '2',
                        ifelse(DataCC$ADNEXgroups == "Stage I invasive", '3',
                        ifelse(DataCC$ADNEXgroups == "Stage II-IV invasive", '4', '5'))))


y = outcome = DataCC$Groups # y
d = response = as.matrix(DataCC[, c('pbenw', 'pborw', 'pst1w', 'pst2_4w', 'pmetaw')]) # d

PDI <- pdI.extented(y = outcome, d = response, k = 5, method = "prob") # 0.55
PDI.se <- ests.extended(y = outcome, d = response, acc = "pdi", level = 0.95, method="prob", k = 5, B = 250, balance = TRUE)$se

LL <- PDI - 1.96*PDI.se
round(LL, 2) # 0.51
UL <- PDI + 1.96*PDI.se
round(UL, 2) # 0.59


#### ADNEX with CA125 ####

## 5 groups
DataCC$pbenwo
DataCC$pborwo
DataCC$pst1wo
DataCC$pst2_4wo
DataCC$pmetawo

## The 10 pairs
# Pair 1: benign vs borderline
DataCC$pp1wo <- DataCC$pbenwo / (DataCC$pbenwo + DataCC$pborwo)
round(auc.nonpara.mw(DataCC[DataCC$ADNEXgroups == "Benign", "pp1wo"], DataCC[DataCC$ADNEXgroups == "Borderline", "pp1wo"], method = "pepe"), 2)
# AUC = 0.85

# Pair 2: benign vs stage I
DataCC$pp2wo <- DataCC$pbenwo / (DataCC$pbenwo + DataCC$pst1wo)
round(auc.nonpara.mw(DataCC[DataCC$ADNEXgroups == "Benign", "pp2wo"], DataCC[DataCC$ADNEXgroups == "Stage I invasive", "pp2wo"], method = "pepe"), 2)
# AUC = 0.92

# Pair 3: benign vs stage II-IV
DataCC$pp3wo <- DataCC$pbenwo / (DataCC$pbenwo + DataCC$pst2_4wo)
round(auc.nonpara.mw(DataCC[DataCC$ADNEXgroups == "Benign", "pp3wo"], DataCC[DataCC$ADNEXgroups == "Stage II-IV invasive", "pp3wo"], method = "pepe"), 2)
# AUC = 0.96

# Pair 4: benign vs metastatic
DataCC$pp4wo <- DataCC$pbenwo / (DataCC$pbenwo + DataCC$pmetawo)
round(auc.nonpara.mw(DataCC[DataCC$ADNEXgroups == "Benign", "pp4wo"], DataCC[DataCC$ADNEXgroups == "Metastatic", "pp4wo"], method = "pepe"), 2)
# AUC = 0.93

# Pair 5: borderline vs stage I
DataCC$pp5wo <- DataCC$pborwo / (DataCC$pborwo + DataCC$pst1wo)
DataCC$binbo <- ifelse(DataCC$ADNEXgroups == "Borderline", 0, 1)
round(auc.nonpara.mw(DataCC[DataCC$ADNEXgroups == "Borderline", "pp5wo"], DataCC[DataCC$ADNEXgroups == "Stage I invasive", "pp5wo"], method = "pepe"), 2)
# AUC = 0.71

# Pair 6: borderline vs stage II-IV
DataCC$pp6wo <- DataCC$pborwo / (DataCC$pborwo + DataCC$pst2_4wo)
round(auc.nonpara.mw(DataCC[DataCC$ADNEXgroups == "Borderline", "pp6wo"], DataCC[DataCC$ADNEXgroups == "Stage II-IV invasive", "pp6wo"], method = "pepe"), 2)
# AUC = 0.91

# Pair 7: borderline vs metastatic
DataCC$pp7wo <- DataCC$pborwo / (DataCC$pborwo + DataCC$pmetawo)
round(auc.nonpara.mw(DataCC[DataCC$ADNEXgroups == "Borderline", "pp7wo"], DataCC[DataCC$ADNEXgroups == "Metastatic", "pp7wo"], method = "pepe"), 2)
# AUC = 0.88

# Pair 8: stage I vs stage II-IV
DataCC$pp8wo <- DataCC$pst1wo / (DataCC$pst1wo + DataCC$pst2_4wo)
DataCC$bini <- ifelse(DataCC$ADNEXgroups == "Stage I invasive", 0, 1)
round(auc.nonpara.mw(DataCC[DataCC$ADNEXgroups == "Stage I invasive", "pp8wo"], DataCC[DataCC$ADNEXgroups == "Stage II-IV invasive", "pp8wo"], method = "pepe"), 2)
# AUC = 0.81

# Pair 9: stage I vs metastatic
DataCC$pp9wo <- DataCC$pst1wo / (DataCC$pst1wo + DataCC$pmetawo)
round(auc.nonpara.mw(DataCC[DataCC$ADNEXgroups == "Stage I invasive", "pp9wo"], DataCC[DataCC$ADNEXgroups == "Metastatic", "pp9wo"], method = "pepe"), 2)
# AUC = 0.74

# Pair 10: stage II-IV vs metastatic
DataCC$pp10wo <- DataCC$pst2_4wo / (DataCC$pst2_4wo + DataCC$pmetawo)
DataCC$bin2 <- ifelse(DataCC$ADNEXgroups == "Stage II-IV invasive", 0, 1)
round(auc.nonpara.mw(DataCC[DataCC$ADNEXgroups == "Stage II-IV invasive", "pp10wo"], DataCC[DataCC$ADNEXgroups == "Metastatic", "pp10wo"], method = "pepe"), 2)
# AUC = 0.59

# Polytomous Discrimination Index (PDI)
library(mcca)

DataCC$Groups <- ifelse(DataCC$ADNEXgroups == "Benign", '1',
                        ifelse(DataCC$ADNEXgroups == "Borderline", '2',
                        ifelse(DataCC$ADNEXgroups == "Stage I invasive", '3',
                        ifelse(DataCC$ADNEXgroups == "Stage II-IV invasive", '4', '5'))))


y = outcome = DataCC$Groups # y
d = response = as.matrix(DataCC[, c('pbenwo', 'pborwo', 'pst1wo', 'pst2_4wo', 'pmetawo')]) # d

PDI <- pdI.extented(y = outcome, d = response, k = 5, method = "prob") # 0.49
PDI.se <- ests.extended(y = outcome, d = response, acc = "pdi", level = 0.95, method="prob", k = 5, B = 250, balance = TRUE)$se

LL <- PDI - 1.96*PDI.se
round(LL, 2) # 0.45
UL <- PDI + 1.96*PDI.se
round(UL, 2) # 0.52


#-----------------------------------------------#
##### 3. Sensitivity analysis: All patients #####
#-----------------------------------------------#

DataVal <- TransIOTAanalyses

##### 3.1 Imputation #####

## Variable indicating whether CA125 is missing
DataVal$CA125missing = factor(sapply(DataVal$`CA125 (kU/L)`, function(x) as.numeric(is.na(x))),
                              levels = 0:1, labels = c("CA125 present", "CA125 missing"))
DataVal[DataVal$CA125missing == "CA125 missing", c("CA125 (kU/L)", "CA125missing")]

## Variable indicating whether HE4 is missing
DataVal$HE4missing = factor(sapply(DataVal$`HE4 (pmol/L)`, function(x) as.numeric(is.na(x))),
                            levels = 0:1, labels = c("HE4 present", "HE4 missing"))
DataVal[DataVal$HE4missing == "HE4 missing", c("HE4 (pmol/L)", "HE4missing")]

## Factor for Nlocules
DataVal$Nlocules <- factor(DataVal$Nlocules, levels = c("1", "2 - 10", "> 10", "Other"))

## Log transformation of 'maximum diameter of lesion'
DataVal$Llesdmax <- log2(DataVal$`Lesion largest diameter`)

## Quadratic term for 'Proportion of solid tissue'
DataVal$Qpropsol <- DataVal$propsol^2

## Factor for CD_5groups
DataVal$ADNEXgroups <- factor(DataVal$ADNEXgroups, levels = c('Benign', 'Borderline', 'Stage I invasive', 'Stage II-IV invasive', 'Metastatic'),
                              labels = c('Benign', 'Borderline', 'Stage I invasive', 'Stage II-IV invasive', 'Metastatic'))

## Transform CA125 and HE4
DataVal$`CA125 (kU/L)` <- as.numeric(DataVal$`CA125 (kU/L)`)
DataVal$`HE4 (pmol/L)` <- as.numeric(DataVal$`HE4 (pmol/L)`)

hist(DataVal$`CA125 (kU/L)`)
hist(DataVal$`HE4 (pmol/L)`)

hist(log(log(DataVal$`CA125 (kU/L)` + 1)))
hist(log(log(DataVal$`HE4 (pmol/L)` + 1)))

DataVal$T_CA125 = log(log(DataVal$`CA125 (kU/L)` + 1))   # Transform the original data

hist(DataVal$T_CA125, prob = TRUE)
lines(density(DataVal$T_CA125[!is.na(DataVal$T_CA125)]))

library(MASS)
Box <- boxcox(DataVal$`HE4 (pmol/L)` ~ 1, plotit = FALSE) # lambda = -0.8
Cox = data.frame(Box$x, Box$y)            # Create a data frame with the results
Cox2 = Cox[with(Cox, order(-Cox$Box.y)),] # Order the new data frame by decreasing y
Cox2[1,]                                  # Display the lambda with the greatest log likelihood
lambda_HE = Cox2[1, "Box.x"]                 # Extract that lambda
DataVal$T_HE4 = (DataVal$`HE4 (pmol/L)` ^ lambda_HE - 1)/lambda_HE   # Transform the original data

hist(DataVal$T_HE4, prob = TRUE)
lines(density(DataVal$T_HE4[!is.na(DataVal$T_HE4)]))

## Variables needed in imputation model
colnames(DataVal)
VarsImp = c('Patient age', 'Llesdmax', 'propsol', 'Qpropsol', 'papnr', 'Nlocules',
            'ascitesbin', 'shadowsbin', 'bilatbin', 'metasbin', 'ADNEXgroups')
VarsAdd = c("T_CA125", "T_HE4", "Patient ID", "Lesion largest diameter", 
            "Type", "Histology.x", "FIGO stage.x", "FIGO a.x", "Additiol details (histology)", 'loc10',
            "Certainty", "soldmax", "multibin", "solidbin", "Postmenopausal2", "Mass Outcome", "CenterAUC", 
            "Histology.y", "FIGO stage.y", "FIGO a.y", "oncocenter", "OutcomeBin", 
            "CA125missing", "HE4missing")
AllVars <- c(VarsImp, VarsAdd)

iotaImp <- DataVal[, AllVars]

## Preparation imputation
PredMatr = matrix(0, length(AllVars), length(AllVars))
colnames(PredMatr) <- AllVars
rownames(PredMatr) <- AllVars

# PredMatr: A value of 1 indicates that the column variable is used as a predictor to impute the target (row) variable, and a 0 means that it is not used

# Prediction Matrix for CA125: the first row needs to be 1, except for ll_CA125. The rest of the matrix is zero
PredMatr['T_CA125', ] <- 1
# Variables present in the dataset that are not necessary for the imputation set to zero
PredMatr[, VarsAdd] <- 0
PredMatr['T_CA125', 'T_HE4'] <- 1
PredMatr

# Prediction Matrix for HE4
PredMatr['T_HE4', ] <- 1
PredMatr['T_HE4', VarsAdd] <- 0
PredMatr['T_HE4', 'T_CA125'] <- 1
PredMatr


colnames(iotaImp) = gsub(" ","\\.",colnames(iotaImp))
rownames(PredMatr)   = colnames(PredMatr) = gsub(" ","\\.",colnames(PredMatr))

Init = mice(iotaImp, m = 1, maxit = 0, predictorMatrix = PredMatr) #, threshold = 1.0)
Init$loggedEvents # 8 logged events
for(i in Init$loggedEvents$out)
  iotaImp[, i] %<>% as.factor
Init$method

Method = sapply(colnames(iotaImp),
                function(x){
                  if(x == "T_CA125")
                    "pmm"
                  else if(x == "T_HE4")
                    "pmm"
                  else
                    ""
                })
Method

## Imputation
iotaSI <- mice(iotaImp, m = 1, me = Method,
               predictorMatrix = PredMatr, maxit = 100,
               seed = 1213, vis = "monotone", ridge = 1e-3)
iotaSI$loggedEvents

## Check convergence
plot(iotaSI)

## Density plots

densityplot(iotaSI, ~ T_CA125| .imp)

densityplot(iotaSI, ~ T_CA125)

densityplot(iotaSI, ~ T_HE4| .imp)

densityplot(iotaSI, ~ T_HE4)

## Boxplot CA125 & HE4
DataImp                = mice::complete(iotaSI, "long")
DataImp$CA125 = exp(exp(DataImp$T_CA125)) - 1
DataImp$HE4 = (lambda_HE * DataImp$T_HE4 + 1)^(1/lambda_HE)

boxplot(T_CA125 ~ .imp, DataImp[DataImp$CA125missing == "CA125 missing",])


boxplot(T_HE4 ~ .imp, DataImp[DataImp$HE4missing == "HE4 missing",])


##### 3.2 Descriptive statistics #####

## Age
median(DataImp$Patient.age, na.rm = T) # 54
quantile(DataImp$Patient.age, na.rm = T, probs = 0.25) # 42
quantile(DataImp$Patient.age, na.rm = T, probs = 0.75) # 64.25
min(DataImp$Patient.age, na.rm = T) # 18
max(DataImp$Patient.age, na.rm = T) # 88

mean(DataImp$Patient.age)

## Menopausal status
table(DataImp$Postmenopausal2) # No: 396; Yes: 536

## Presence of solid components
table(DataImp$solidbin) # 0: 296; 1: 636

## Observed CA125
median(DataVal$`CA125 (kU/L)`) # 26.12
quantile(DataVal$`CA125 (kU/L)`, probs = 0.25) # 14.51
quantile(DataVal$`CA125 (kU/L)`, probs = 0.75) # 133
min(DataVal$`CA125 (kU/L)`) # 3.2
max(DataVal$`CA125 (kU/L)`) # 4798

## Observed HE4
median(DataVal$`HE4 (pmol/L)`) # 63.42
quantile(DataVal$`HE4 (pmol/L)`, probs = 0.25) # 49.78
quantile(DataVal$`HE4 (pmol/L)`, probs = 0.75) # 131.8
min(DataVal$`HE4 (pmol/L)`) # 26.82
max(DataVal$`HE4 (pmol/L)`) # 8260

## Maximum diameter of lesion (mm)
median(DataImp$Lesion.largest.diameter) # 75
quantile(DataImp$Lesion.largest.diameter, probs = 0.25) # 49
quantile(DataImp$Lesion.largest.diameter, probs = 0.75) # 117
min(DataImp$Lesion.largest.diameter) # 7
max(DataImp$Lesion.largest.diameter) # 459

## Proportion solid tissue
median(DataImp$propsol) # 0.38
quantile(DataImp$propsol, probs = 0.25) # 0
quantile(DataImp$propsol, probs = 0.75) # 1
min(DataImp$propsol) # 0
max(DataImp$propsol) # 1

## Number of papillary projections
table(DataImp$papnr)

## Bilateral masses
table(DataImp$bilatbin) # 0: 691; 1: 241

## More than 10 cyst locules
table(DataImp$loc10) # 0: 789; 1: 143

## Acoustic shadows
table(DataImp$shadowsbin) # 0: 594; 1: 338

table(DataImp$shadowsbin, DataImp$CenterAUC)
table(DataImp$shadowsbin, DataImp$ADNEXgroups)

table(DataImp$CenterAUC[DataImp$shadowsbin == 1], DataImp$ADNEXgroups[DataImp$shadowsbin == 1])

## Ascites
table(DataImp$ascitesbin) # 0: 791; 1: 141

## Multilocular cysts
table(DataImp$multibin) # 0: 527; 1: 405

## Metastases 
table(DataImp$metasbin) #0: 798; 1: 134

## Ultrasound examiner's subjective impression
table(DataImp$Certainty)

## Outcome
table(DataImp$ADNEXgroups)


##### 3.3 Validation #####

##### 3.3.1 Model predictions #####

##### ADNEX #####

## With CA125
DataImp = ADNEX(Age = Patient.age, Oncocenter = oncocenter, lesdmax = Lesion.largest.diameter, soldmaxorig = soldmax, loc10 = loc10, Papnr = papnr, Shadows = shadowsbin, Ascites = ascitesbin, ca125 = CA125, wica = T, woca = F, data = DataImp)

## Without CA125
DataImp = ADNEX(Age = Patient.age, Oncocenter = oncocenter, lesdmax = Lesion.largest.diameter, soldmaxorig = soldmax, loc10 = loc10, Papnr = papnr, Shadows = shadowsbin, Ascites = ascitesbin, ca125 = CA125, wica = F, woca = T, data = DataImp)

##### ROMA #####

## HE4 expressed in pmol/L
## CA125 expressed in U/mL (= kU/L)

DataImp <- ddply(.data = DataImp, .(Patient.ID),
                 function(x){
                   x$ROMA =
                     if(x$Postmenopausal2 == "yes"){
                       z = -8.09 + (1.04 * log(x$HE4)) + (0.732 * log(x$CA125))
                       p = 1/(1 + exp(-z))
                       p
                     }else{
                       z = -12.0 + (2.38 * log(x$HE4)) + (0.0626 * log(x$CA125))
                       p = 1/(1 + exp(-z))
                       p
                     }
                   return(x)
                 })

##### 3.3.2 Discrimination between benign and malignant tumors #####

##### ADNEX with CA125 #####

## AUC
AUC.ADNEXwSens <- AUC.IOTA(pred = pmalw, outcome = OutcomeBin, center = CenterAUC, data = DataImp, method.MA = "BAYES", titleGraph = "AUC per center")

AUC.ADNEXwSens$Plot[9, "RRauc"] <- "   "
AUC.ADNEXwSens$Plot[10, "RRauc"] <- "   "
AUC.ADNEXwSens$Plot[11, "RRcenter"] <- "AUC (95% CI)"
AUC.ADNEXwSens$Plot[12, "RRauc"] <- "        (0.82 to 0.97)"
AUC.ADNEXwSens$Plot[12, "RRprev"] <- "   "

forestplot(AUC.ADNEXwSens$Plot,
           align = c("l", "c", "c"),
           mean = AUC.ADNEXwSens$dataPlot$AUC,
           lower = AUC.ADNEXwSens$dataPlot$LL,
           upper = AUC.ADNEXwSens$dataPlot$UL,
           is.summary = c(TRUE, FALSE, rep(FALSE, length(AUC.ADNEXwSens$IncludedCenters)), FALSE, FALSE, TRUE, TRUE),
           xlab = "AUC", 
           xticks = c(0.70, 0.80, 0.9, 1),
           clip = c(0.71, 1),
           txt_gp = fpTxtGp(summary = gpar(cex = 1.55), xlab = gpar(cex = 1.5), label = gpar(cex = 1.5),
                            ticks = gpar(cex = 1.5), title = gpar(cex = 1.75)),
           graphwidth = unit(6, "cm"),
           col = fpColors(box = "gray", lines = "black", summary = "black"),
           fn.ci_sum = c(as.list(rep("fpDrawSummaryCI", length(AUC.ADNEXwSens$IncludedCenters) + 4)),
                         fpDrawSummaryCI,
                         fpDrawBarCI),
           zero = AUC.ADNEXwSens$Performance$AUC)
# 900 x 450

## Sensitivity and specificity
SensSpec(pred = pmalw, outcome = OutcomeBin, threshold = 0.01,  center = CenterAUC, data = DataImp)$OverallPer
SensSpec(pred = pmalw, outcome = OutcomeBin, threshold = 0.03,  center = CenterAUC, data = DataImp)$OverallPer
SensSpec(pred = pmalw, outcome = OutcomeBin, threshold = 0.05,  center = CenterAUC, data = DataImp)$OverallPer
SensSpec(pred = pmalw, outcome = OutcomeBin, threshold = 0.10,  center = CenterAUC, data = DataImp)$OverallPer
SensSpec(pred = pmalw, outcome = OutcomeBin, threshold = 0.15,  center = CenterAUC, data = DataImp)$OverallPer
SensSpec(pred = pmalw, outcome = OutcomeBin, threshold = 0.20,  center = CenterAUC, data = DataImp)$OverallPer
SensSpec(pred = pmalw, outcome = OutcomeBin, threshold = 0.25,  center = CenterAUC, data = DataImp)$OverallPer
SensSpec(pred = pmalw, outcome = OutcomeBin, threshold = 0.30,  center = CenterAUC, data = DataImp)$OverallPer
SensSpec(pred = pmalw, outcome = OutcomeBin, threshold = 0.40,  center = CenterAUC, data = DataImp)$OverallPer
SensSpec(pred = pmalw, outcome = OutcomeBin, threshold = 0.50,  center = CenterAUC, data = DataImp)$OverallPer

## Sensitivity at fixed specificity
FixSpec(pred = pmalw, outcome = OutcomeBin, Specificity = 0.90, center = CenterAUC, data = DataImp)$OverallPer

## Specificity at fixed sensitivity
FixSens(pred = pmalw, outcome = OutcomeBin, Sensitivity = 0.90, center = CenterAUC, data = DataImp)$OverallPer

##### ADNEX without CA125 #####

## AUC
AUC.ADNEXwoSens <- AUC.IOTA(pred = pmalwo, outcome = OutcomeBin, center = CenterAUC, data = DataImp, method.MA = "BAYES", titleGraph = "AUC per center")

AUC.ADNEXwoSens$Plot[9, "RRauc"] <- "   "
AUC.ADNEXwoSens$Plot[10, "RRauc"] <- "   "
AUC.ADNEXwoSens$Plot[11, "RRcenter"] <- "AUC (95% CI)"
AUC.ADNEXwoSens$Plot[12, "RRauc"] <- "        (0.70 to 0.97)"
AUC.ADNEXwoSens$Plot[12, "RRprev"] <- "   "

forestplot(AUC.ADNEXwoSens$Plot,
           align = c("l", "c", "c"),
           mean = AUC.ADNEXwoSens$dataPlot$AUC,
           lower = AUC.ADNEXwoSens$dataPlot$LL,
           upper = AUC.ADNEXwoSens$dataPlot$UL,
           is.summary = c(TRUE, FALSE, rep(FALSE, length(AUC.ADNEXwoSens$IncludedCenters)), FALSE, FALSE, TRUE, TRUE),
           xlab = "AUC", 
           xticks = c(0.70, 0.80, 0.9, 1),
           clip = c(0.71, 1),
           txt_gp = fpTxtGp(summary = gpar(cex = 1.55), xlab = gpar(cex = 1.5), label = gpar(cex = 1.5),
                            ticks = gpar(cex = 1.5), title = gpar(cex = 1.75)),
           graphwidth = unit(6, "cm"),
           col = fpColors(box = "gray", lines = "black", summary = "black"),
           fn.ci_sum = c(as.list(rep("fpDrawSummaryCI", length(AUC.ADNEXwoSens$IncludedCenters) + 4)),
                         fpDrawSummaryCI,
                         fpDrawBarCI),
           zero = AUC.ADNEXwoSens$Performance$AUC)
# 900 x 450

## Sensitivity and specificity
SensSpec(pred = pmalwo, outcome = OutcomeBin, threshold = 0.01, center = CenterAUC, data = DataImp)$OverallPer
SensSpec(pred = pmalwo, outcome = OutcomeBin, threshold = 0.03, center = CenterAUC, data = DataImp)$OverallPer
SensSpec(pred = pmalwo, outcome = OutcomeBin, threshold = 0.05, center = CenterAUC, data = DataImp)$OverallPer
SensSpec(pred = pmalwo, outcome = OutcomeBin, threshold = 0.10, center = CenterAUC, data = DataImp)$OverallPer
SensSpec(pred = pmalwo, outcome = OutcomeBin, threshold = 0.15, center = CenterAUC, data = DataImp)$OverallPer
SensSpec(pred = pmalwo, outcome = OutcomeBin, threshold = 0.20, center = CenterAUC, data = DataImp)$OverallPer
SensSpec(pred = pmalwo, outcome = OutcomeBin, threshold = 0.25, center = CenterAUC, data = DataImp)$OverallPer
SensSpec(pred = pmalwo, outcome = OutcomeBin, threshold = 0.30, center = CenterAUC, data = DataImp)$OverallPer
SensSpec(pred = pmalwo, outcome = OutcomeBin, threshold = 0.40, center = CenterAUC, data = DataImp)$OverallPer
SensSpec(pred = pmalwo, outcome = OutcomeBin, threshold = 0.50, center = CenterAUC, data = DataImp)$OverallPer

## Sensitivity at fixed specificity
FixSpec(pred = pmalwo, outcome = OutcomeBin, Specificity = 0.90, center = CenterAUC, data = DataImp)$OverallPer

## Specificity at fixed sensitivity
FixSens(pred = pmalwo, outcome = OutcomeBin, Sensitivity = 0.90, center = CenterAUC, data = DataImp)$OverallPer

##### ROMA #####

## AUC
AUC.ROMAsens <- AUC.IOTA(pred = ROMA, outcome = OutcomeBin, center = CenterAUC, data = DataImp, method.MA = "BAYES", titleGraph = "AUC per center")

AUC.ROMAsens$Plot[9, "RRauc"] <- "   "
AUC.ROMAsens$Plot[10, "RRauc"] <- "   "
AUC.ROMAsens$Plot[11, "RRcenter"] <- "AUC (95% CI)"
AUC.ROMAsens$Plot[12, "RRauc"] <- "        (0.72 to 0.92)"
AUC.ROMAsens$Plot[12, "RRprev"] <- "   "

forestplot(AUC.ROMAsens$Plot,
           align = c("l", "c", "c"),
           mean = AUC.ROMAsens$dataPlot$AUC,
           lower = AUC.ROMAsens$dataPlot$LL,
           upper = AUC.ROMAsens$dataPlot$UL,
           is.summary = c(TRUE, FALSE, rep(FALSE, length(AUC.ROMAsens$IncludedCenters)), FALSE, FALSE, TRUE, TRUE),
           xlab = "AUC",
           xticks = c(0.70, 0.80, 0.9, 1),
           clip = c(0.71, 1),
           txt_gp = fpTxtGp(summary = gpar(cex = 1.55), xlab = gpar(cex = 1.5), label = gpar(cex = 1.5),
                            ticks = gpar(cex = 1.5), title = gpar(cex = 1.75)),
           graphwidth = unit(6, "cm"),
           col = fpColors(box = "gray", lines = "black", summary = "black"),
           fn.ci_sum = c(as.list(rep("fpDrawSummaryCI", length(AUC.ROMAsens$IncludedCenters) + 4)),
                         fpDrawSummaryCI,
                         fpDrawBarCI),
           zero = AUC.ROMAsens$Performance$AUC)
# 900 x 450

## Sensitivity and specificity
SensSpec(pred = ROMA, outcome = OutcomeBin, threshold = 0.01, center = CenterAUC, data = DataImp)$OverallPer
SensSpec(pred = ROMA, outcome = OutcomeBin, threshold = 0.03, center = CenterAUC, data = DataImp)$OverallPer
SensSpec(pred = ROMA, outcome = OutcomeBin, threshold = 0.05, center = CenterAUC, data = DataImp)$OverallPer
SensSpec(pred = ROMA, outcome = OutcomeBin, threshold = 0.10, center = CenterAUC, data = DataImp)$OverallPer
SensSpec(pred = ROMA, outcome = OutcomeBin, threshold = 0.15, center = CenterAUC, data = DataImp)$OverallPer
SensSpec(pred = ROMA, outcome = OutcomeBin, threshold = 0.20, center = CenterAUC, data = DataImp)$OverallPer
SensSpec(pred = ROMA, outcome = OutcomeBin, threshold = 0.25, center = CenterAUC, data = DataImp)$OverallPer
SensSpec(pred = ROMA, outcome = OutcomeBin, threshold = 0.30, center = CenterAUC, data = DataImp)$OverallPer
SensSpec(pred = ROMA, outcome = OutcomeBin, threshold = 0.40, center = CenterAUC, data = DataImp)$OverallPer
SensSpec(pred = ROMA, outcome = OutcomeBin, threshold = 0.50, center = CenterAUC, data = DataImp)$OverallPer

## ROMA cut-offs

# Postmenopausal
SensSpec(pred = ROMA, outcome = OutcomeBin, threshold = 0.299, center = CenterAUC, data = DataImp[DataImp$Postmenopausal2 == "yes",])$OverallPer

ROMApostSens <- SensSpec(pred = ROMA, outcome = OutcomeBin, threshold = 0.299, center = CenterAUC, data = DataImp[DataImp$Postmenopausal2 == "yes",])$CenterPer

# Premenopausal
SensSpec(pred = ROMA, outcome = OutcomeBin, threshold = 0.114, center = CenterAUC, data = DataImp[DataImp$Postmenopausal2 == "no",])$OverallPer

ROMApreSens <- SensSpec(pred = ROMA, outcome = OutcomeBin, threshold = 0.114, center = CenterAUC, data = DataImp[DataImp$Postmenopausal2 == "no",])$CenterPer

# Combined
SensSpecROMA <- matrix(nrow = 6, ncol = 5)
colnames(SensSpecROMA) <- c("Center", "TP", "TN", "FP", "FN")
SensSpecROMA <- data.frame(SensSpecROMA)

SensSpecROMA$Center <- c("Leuven", "Prague", "London", "Genk", "Milan", "Rome")

SensSpecROMA$TN     <- ROMApostSens$TN + ROMApreSens$TN
SensSpecROMA$TP     <- ROMApostSens$TP + ROMApreSens$TP
SensSpecROMA$FN     <- ROMApostSens$FN + ROMApreSens$FN
SensSpecROMA$FP     <- ROMApostSens$FP + ROMApreSens$FP

SensSpecROMA$Sensitivity  <- SensSpecROMA$TP / (SensSpecROMA$TP + SensSpecROMA$FN)
SensSpecROMA$se_sens      <- sqrt((SensSpecROMA$Sensitivity * (1 - SensSpecROMA$Sensitivity))/(SensSpecROMA$TP + SensSpecROMA$FN))
SensSpecROMA$Specificity  <- SensSpecROMA$TN / (SensSpecROMA$TN + SensSpecROMA$FP)
SensSpecROMA$se_spec      <- sqrt((SensSpecROMA$Specificity * (1 - SensSpecROMA$Specificity))/(SensSpecROMA$TN + SensSpecROMA$FP))

# Logit transformation
SensSpecROMA$logit.sens     <- logit(SensSpecROMA$Sensitivity)
SensSpecROMA$logit.se.sens  <- sqrt(1/SensSpecROMA$TP + 1/SensSpecROMA$FN)
SensSpecROMA$logit.LL.sens  <- SensSpecROMA$logit.sens - 1.96*SensSpecROMA$logit.se.sens
SensSpecROMA$logit.UL.sens  <- SensSpecROMA$logit.sens + 1.96*SensSpecROMA$logit.se.sens

SensSpecROMA$logit.spec     <- logit(SensSpecROMA$Specificity)
SensSpecROMA$logit.se.spec  <- sqrt(1/SensSpecROMA$TN + 1/SensSpecROMA$FP)
SensSpecROMA$logit.LL.spec  <- SensSpecROMA$logit.spec - 1.96*SensSpecROMA$logit.se.spec
SensSpecROMA$logit.UL.spec  <- SensSpecROMA$logit.spec + 1.96*SensSpecROMA$logit.se.spec

Sensoverall <- matrix(nrow = 1, ncol = 7)
colnames(Sensoverall) <- c('Center', 'Sens', 'LL.sens', 'UL.sens', 'Spec', 'LL.spec', 'UL.spec')
Sensoverall <- data.frame(Sensoverall)
Sensoverall$Center <- "Overall"

# Meta-analysis for overall estimate: Bivariate random-effects model
p1 <- 2*6
p2 <- 4*6
Sensbi <- SensSpecROMA[, c("Center", "logit.sens", "logit.spec", "logit.se.sens", "logit.se.spec")]
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

## Sensitivity at fixed specificity
FixSpec(pred = ROMA, outcome = OutcomeBin, Specificity = 0.90, center = CenterAUC, data = DataImp)$OverallPer

## Specificity at fixed sensitivity
FixSens(pred = ROMA, outcome = OutcomeBin, Sensitivity = 0.90, center = CenterAUC, data = DataImp)$OverallPer

##### Summary plot for AUC #####

NA.forest <- AUC.ROMAsens
NA.forest <- NA

Summary.AUC <- rbind(NA.forest, AUC.ROMAsens$Performance[1,], AUC.ADNEXwoSens$Performance[1,], AUC.ADNEXwSens$Performance[1,])
Summary.AUC$Model <- c('', 'ROMA', 'ADNEX without CA125', 'ADNEX with CA125')

Summary.PI <- rbind(NA.forest, AUC.ROMAsens$Performance[2,], AUC.ADNEXwoSens$Performance[2,], AUC.ADNEXwSens$Performance[2,])
Summary.PI$Model <- c('', 'ROMA', 'ADNEX without CA125', 'ADNEX with CA125')

tabletext <- cbind(
  c('Model', 'ROMA', 'ADNEX without CA125', 'ADNEX with CA125'),
  c('AUC (95% CI)', 
    paste(format(round(AUC.ROMAsens$Performance$AUC[1], 2), nsmall = 2), " (", format(round(AUC.ROMAsens$Performance$LL[1], 2), nsmall = 2), "; ", format(round(AUC.ROMAsens$Performance$UL[1], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(AUC.ADNEXwoSens$Performance$AUC[1], 2), nsmall = 2), " (", format(round(AUC.ADNEXwoSens$Performance$LL[1], 2), nsmall = 2), "; ", format(round(AUC.ADNEXwoSens$Performance$UL[1], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(AUC.ADNEXwSens$Performance$AUC[1], 2), nsmall = 2), " (", format(round(AUC.ADNEXwSens$Performance$LL[1], 2), nsmall = 2), "; ", format(round(AUC.ADNEXwSens$Performance$UL[1], 2), nsmall = 2), ")", sep = "")),
  c('95% PI', 
    paste("(", format(round(AUC.ROMAsens$Performance$LL[2], 2), nsmall = 2), "; ", format(round(AUC.ROMAsens$Performance$UL[2], 2), nsmall = 2), ")", sep = ""),
    paste("(", format(round(AUC.ADNEXwoSens$Performance$LL[2], 2), nsmall = 2), "; ", format(round(AUC.ADNEXwoSens$Performance$UL[2], 2), nsmall = 2), ")", sep = ""),
    paste("(", format(round(AUC.ADNEXwSens$Performance$LL[2], 2), nsmall = 2), "; ", format(round(AUC.ADNEXwSens$Performance$UL[2], 2), nsmall = 2), ")", sep = ""))
)

forestplot(labeltext = tabletext,
           title = "",
           mean = round(Summary.AUC$AUC, 3),
           lower = round(Summary.AUC$LL, 3),
           upper = round(Summary.AUC$UL, 3),
           is.summary = c(FALSE, TRUE, TRUE, TRUE),
           xlab = "AUC (95% CI)",
           boxsize = .5,
           txt_gp = fpTxtGp(summary = gpar(cex = 1.55, fontface = "plain"), xlab = gpar(cex = 1.5, fontface = "plain"), label = gpar(cex = 1.5, fontface = "bold"),
                            ticks = gpar(cex = 1.5, fontface = "plain"), title = gpar(cex = 1.75)),
           graphwidth = unit(9, "cm"),
           graph.pos = 3,
           xticks = c(0.80, 0.85, 0.9, 0.95, 1), xlog = TRUE, clip = c(0.80, 1))
# 1050 x 400


##### 3.3.3 Calibration of the risk of malignancy #####

##### ADNEX with CA125 #####

Cal.ADNEXwSens <- Calibration(p = pmalw, y = OutcomeBin, center = CenterAUC, data = DataImp, CalibrLines = "overall", MethodCL = "Wald", nr.knots = 3) 

##### ADNEX without CA125 #####

Cal.ADNEXwoSens <- Calibration(p = pmalwo, y = OutcomeBin, center = CenterAUC, data = DataImp, CalibrLines = "overall", MethodCL = "Wald") 

##### ROMA #####

Cal.ROMAsens <- Calibration(p = ROMA, y = OutcomeBin, center = CenterAUC, data = DataImp, CalibrLines = "overall", MethodCL = "Wald") 

##### Summary plot for Calibration #####

table <- matrix(ncol = 3, nrow = 3)
colnames(table) <- c('Model', 'Intercept (95% CI)', 'Slope (95% CI)')
table[, 1] <- c('ROMA', 'ADNEX without CA125', 'ADNEX with CA125')
table[1, 2:3] <- c(paste0(format(round(Cal.ROMAsens$Performance[1,1], 2), nsmall = 2), " (", format(round(Cal.ROMAsens$Performance[1,2], 2), nsmall = 2), " to ", format(round(Cal.ROMAsens$Performance[1,3], 2), nsmall = 2), ")"), paste0(format(round(Cal.ROMAsens$Performance[2,1], 2), nsmall = 2), " (", format(round(Cal.ROMAsens$Performance[2,2], 2), nsmall = 2), " to ", format(round(Cal.ROMAsens$Performance[2,3], 2), nsmall = 2), ")"))
table[3, 2:3] <- c(paste0(format(round(Cal.ADNEXwSens$Performance[1,1], 2), nsmall = 2), " (", format(round(Cal.ADNEXwSens$Performance[1,2], 2), nsmall = 2), " to ", format(round(Cal.ADNEXwSens$Performance[1,3], 2), nsmall = 2), ")"), paste0(format(round(Cal.ADNEXwSens$Performance[2,1], 2), nsmall = 2), " (", format(round(Cal.ADNEXwSens$Performance[2,2], 2), nsmall = 2), " to ", format(round(Cal.ADNEXwSens$Performance[2,3], 2), nsmall = 2), ")"))
table[2, 2:3] <- c(paste0(format(round(Cal.ADNEXwoSens$Performance[1,1], 2), nsmall = 2), " (", format(round(Cal.ADNEXwoSens$Performance[1,2], 2), nsmall = 2), " to ", format(round(Cal.ADNEXwoSens$Performance[1,3], 2), nsmall = 2), ")"), paste0(format(round(Cal.ADNEXwoSens$Performance[2,1], 2), nsmall = 2), " (", format(round(Cal.ADNEXwoSens$Performance[2,2], 2), nsmall = 2), " to ", format(round(Cal.ADNEXwoSens$Performance[2,3], 2), nsmall = 2), ")"))

library(plotrix)

## Make graph
x = seq(0, 1, by = 0.05)
y = seq(0, 1, by = 0.05)
tiff("Summary calibration - Sensitivity.tiff", width = 14, height = 14, units = "cm", res = 300)
plot(x, y, xlim = c(0,1), ylim = c(0,1), type = "l", col = "gray50", lwd = 2, lty = 2, 
     xlab = "Estimated risk", ylab = "Observed proportion",
     main = "", cex.lab = 1, cex.axis = 1, las = 1) 
lines(Cal.ROMAsens$Prediction$p, Cal.ROMAsens$Prediction$yhat, lwd = 2, col = "maroon")
lines(Cal.ADNEXwoSens$Prediction$p, Cal.ADNEXwoSens$Prediction$yhat, lwd = 2, col = "darkgreen")
lines(Cal.ADNEXwSens$Prediction$p, Cal.ADNEXwSens$Prediction$yhat, lwd = 2, col = "darkorange")
legend(x = -0.035, y = 1, legend = c( "Ideal", "ROMA", "ADNEX without CA125", "ADNEX with CA125"),
       col = c("gray50", "maroon", "darkgreen", "darkorange"), lty = c(2,1,1,1), lwd = 2, cex = 0.7, bty = "n", ncol = 1)
addtable2plot(x = 0.17, y = 0.02, table = table, display.colnames= TRUE, cex = 0.7)
dev.off()


##### 3.3.4 Clinical utility #####

##### ADNEX with CA125 #####

Dec.ADNEXwSens <- DataWinBugs(pred = pmalw, outcome = OutcomeBin, center = CenterAUC, data = DataImp)$Results

Dec.ADNEXwSenslong <- rbindlist(Dec.ADNEXwSens)
Dec.ADNEXwSenslong$prev <- (Dec.ADNEXwSenslong$TP + Dec.ADNEXwSenslong$FN) / Dec.ADNEXwSenslong$n
Dec.ADNEXwSenslong$TA <- Dec.ADNEXwSenslong$prev - (1 - Dec.ADNEXwSenslong$prev) * (Dec.ADNEXwSenslong$CutOff/(1 - Dec.ADNEXwSenslong$CutOff))

ggplot(data = Dec.ADNEXwSenslong, aes(x = CutOff, y = NB, colour = Center)) + 
  geom_line() +
  geom_line(aes(x = CutOff, y = TA))

##### ADNEX without CA125 #####

Dec.ADNEXwoSens <- DataWinBugs(pred = pmalwo, outcome = OutcomeBin, center = CenterAUC, data = DataImp)$Results

Dec.ADNEXwoSenslong <- rbindlist(Dec.ADNEXwoSens)
Dec.ADNEXwoSenslong$prev <- (Dec.ADNEXwoSenslong$TP + Dec.ADNEXwoSenslong$FN) / Dec.ADNEXwoSenslong$n
Dec.ADNEXwoSenslong$TA <- Dec.ADNEXwoSenslong$prev - (1 - Dec.ADNEXwoSenslong$prev) * (Dec.ADNEXwoSenslong$CutOff/(1 - Dec.ADNEXwoSenslong$CutOff))

ggplot(data = Dec.ADNEXwoSenslong, aes(x = CutOff, y = NB, colour = Center)) + 
  geom_line() +
  geom_line(aes(x = CutOff, y = TA))

##### ROMA #####

Dec.ROMAsens <- DataWinBugs(pred = ROMA, outcome = OutcomeBin, center = CenterAUC, data = DataImp)$Results

Dec.ROMAsenslong <- rbindlist(Dec.ROMAsens)
Dec.ROMAsenslong$prev <- (Dec.ROMAsenslong$TP + Dec.ROMAsenslong$FN) / Dec.ROMAsenslong$n
Dec.ROMAsenslong$TA <- Dec.ROMAsenslong$prev - (1 - Dec.ROMAsenslong$prev) * (Dec.ROMAsenslong$CutOff/(1 - Dec.ROMAsenslong$CutOff))

ggplot(data = Dec.ROMAsenslong, aes(x = CutOff, y = NB, colour = Center)) + 
  geom_line() +
  geom_line(aes(x = CutOff, y = TA))

##### Summary plot for Clinical utility #####

Overview <- read_excel("Results WinBugs - Imp.xlsx", sheet = "Overview")

tiff("Summary Clinical utility - Sensitivity.tiff", width = 14, height = 14, units = "cm", res = 300)
plot(Overview$`Risk threshold`, Overview$`Treat all`, type = "l", col = "grey", xlab = "Risk threshold", ylab = "Net benefit", xlim = c(0.055, 0.5), ylim = c(-0.03, 0.55), lwd = 2, cex.lab = 1, cex.axis = 1, las = 1)
points(Overview$`Risk threshold`, Overview$`NB ROMA`, type = 'l', col = "maroon", lwd = 2)
points(Overview$`Risk threshold`, Overview$`NB ADNEX wo`, type = 'l', col = "darkgreen", lwd = 2)
points(Overview$`Risk threshold`, Overview$`NB ADNEXw`, type = 'l', col = "darkorange", lwd = 2)
abline(a = 0, b = 0, lty = 2, lwd = 2)
legend(x = 0.04, y = 0.55, legend = c("ROMA", "ADNEX without CA125", "ADNEX with CA125", "Treat all", "Treat none"), ncol = 3, #cex = 1.2, #bty = "n", #cex = 0.7, pt.cex = 1,
       col = c("maroon", "darkgreen", "darkorange", "grey", 'black'), lty = c(1,1,1,1,2), lwd = 2, cex = 0.7, bty = "n")
dev.off()


##### 3.3.5 Multinomial discrimination #####

#### ADNEX with CA125 ####

## 5 groups
DataImp$pbenw
DataImp$pborw
DataImp$pst1w
DataImp$pst2_4w
DataImp$pmetaw

## The 10 pairs
# Pair 1: benign vs borderline
DataImp$pp1 <- DataImp$pbenw / (DataImp$pbenw + DataImp$pborw)
round(auc.nonpara.mw(DataImp[DataImp$ADNEXgroups == "Benign", "pp1"], DataImp[DataImp$ADNEXgroups == "Borderline", "pp1"], method = "pepe"), 2)
# AUC = 0.86

# Pair 2: benign vs stage I
DataImp$pp2 <- DataImp$pbenw / (DataImp$pbenw + DataImp$pst1w)
round(auc.nonpara.mw(DataImp[DataImp$ADNEXgroups == "Benign", "pp2"], DataImp[DataImp$ADNEXgroups == "Stage I invasive", "pp2"], method = "pepe"), 2)
# AUC = 0.93

# Pair 3: benign vs stage II-IV
DataImp$pp3 <- DataImp$pbenw / (DataImp$pbenw + DataImp$pst2_4w)
round(auc.nonpara.mw(DataImp[DataImp$ADNEXgroups == "Benign", "pp3"], DataImp[DataImp$ADNEXgroups == "Stage II-IV invasive", "pp3"], method = "pepe"), 2)
# AUC = 0.98

# Pair 4: benign vs metastatic
DataImp$pp4 <- DataImp$pbenw / (DataImp$pbenw + DataImp$pmetaw)
round(auc.nonpara.mw(DataImp[DataImp$ADNEXgroups == "Benign", "pp4"], DataImp[DataImp$ADNEXgroups == "Metastatic", "pp4"], method = "pepe"), 2)
# AUC = 0.93

# Pair 5: borderline vs stage I
DataImp$pp5 <- DataImp$pborw / (DataImp$pborw + DataImp$pst1w)
DataImp$binbo <- ifelse(DataImp$ADNEXgroups == "Borderline", 0, 1)
round(auc.nonpara.mw(DataImp[DataImp$ADNEXgroups == "Borderline", "pp5"], DataImp[DataImp$ADNEXgroups == "Stage I invasive", "pp5"], method = "pepe"), 2)
# AUC = 0.71

# Pair 6: borderline vs stage II-IV
DataImp$pp6 <- DataImp$pborw / (DataImp$pborw + DataImp$pst2_4w)
round(auc.nonpara.mw(DataImp[DataImp$ADNEXgroups == "Borderline", "pp6"], DataImp[DataImp$ADNEXgroups == "Stage II-IV invasive", "pp6"], method = "pepe"), 2)
# AUC = 0.93

# Pair 7: borderline vs metastatic
DataImp$pp7 <- DataImp$pborw / (DataImp$pborw + DataImp$pmetaw)
round(auc.nonpara.mw(DataImp[DataImp$ADNEXgroups == "Borderline", "pp7"], DataImp[DataImp$ADNEXgroups == "Metastatic", "pp7"], method = "pepe"), 2)
# AUC = 0.86

# Pair 8: stage I vs stage II-IV
DataImp$pp8 <- DataImp$pst1w / (DataImp$pst1w + DataImp$pst2_4w)
DataImp$bini <- ifelse(DataImp$ADNEXgroups == "Stage I invasive", 0, 1)
round(auc.nonpara.mw(DataImp[DataImp$ADNEXgroups == "Stage I invasive", "pp8"], DataImp[DataImp$ADNEXgroups == "Stage II-IV invasive", "pp8"], method = "pepe"), 2)
# AUC = 0.87

# Pair 9: stage I vs metastatic
DataImp$pp9 <- DataImp$pst1w / (DataImp$pst1w + DataImp$pmetaw)
round(auc.nonpara.mw(DataImp[DataImp$ADNEXgroups == "Stage I invasive", "pp9"], DataImp[DataImp$ADNEXgroups == "Metastatic", "pp9"], method = "pepe"), 2)
# AUC = 0.73

# Pair 10: stage II-IV vs metastatic
DataImp$pp10 <- DataImp$pst2_4w / (DataImp$pst2_4w + DataImp$pmetaw)
DataImp$bin2 <- ifelse(DataImp$ADNEXgroups == "Stage II-IV invasive", 0, 1)
round(auc.nonpara.mw(DataImp[DataImp$ADNEXgroups == "Stage II-IV invasive", "pp10"], DataImp[DataImp$ADNEXgroups == "Metastatic", "pp10"], method = "pepe"), 2)
# AUC = 0.86

# Polytomous Discrimination Index (PDI)
library(mcca)

DataImp$Groups <- ifelse(DataImp$ADNEXgroups == "Benign", '1',
                         ifelse(DataImp$ADNEXgroups == "Borderline", '2',
                         ifelse(DataImp$ADNEXgroups == "Stage I invasive", '3',
                         ifelse(DataImp$ADNEXgroups == "Stage II-IV invasive", '4', '5'))))

y = outcome = DataImp$Groups # y
d = response = as.matrix(DataImp[, c('pbenw', 'pborw', 'pst1w', 'pst2_4w', 'pmetaw')]) # d

PDI <- pdI.extented(y = outcome, d = response, k = 5, method = "prob") # 0.56
PDI.se <- ests.extended(y = outcome, d = response, acc = "pdi", level = 0.95, method="prob", k = 5, B = 250, balance = TRUE)$se

LL <- PDI - 1.96*PDI.se
round(LL, 2) # 0.52
UL <- PDI + 1.96*PDI.se
round(UL, 2) # 0.60


#### ADNEX with CA125 ####

## 5 groups
DataImp$pbenwo
DataImp$pborwo
DataImp$pst1wo
DataImp$pst2_4wo
DataImp$pmetawo

## The 10 pairs
# Pair 1: benign vs borderline
DataImp$pp1wo <- DataImp$pbenwo / (DataImp$pbenwo + DataImp$pborwo)
round(auc.nonpara.mw(DataImp[DataImp$ADNEXgroups == "Benign", "pp1wo"], DataImp[DataImp$ADNEXgroups == "Borderline", "pp1wo"], method = "pepe"), 2)
# AUC = 0.85

# Pair 2: benign vs stage I
DataImp$pp2wo <- DataImp$pbenwo / (DataImp$pbenwo + DataImp$pst1wo)
round(auc.nonpara.mw(DataImp[DataImp$ADNEXgroups == "Benign", "pp2wo"], DataImp[DataImp$ADNEXgroups == "Stage I invasive", "pp2wo"], method = "pepe"), 2)
# AUC = 0.92

# Pair 3: benign vs stage II-IV
DataImp$pp3wo <- DataImp$pbenwo / (DataImp$pbenwo + DataImp$pst2_4wo)
round(auc.nonpara.mw(DataImp[DataImp$ADNEXgroups == "Benign", "pp3wo"], DataImp[DataImp$ADNEXgroups == "Stage II-IV invasive", "pp3wo"], method = "pepe"), 2)
# AUC = 0.96

# Pair 4: benign vs metastatic
DataImp$pp4wo <- DataImp$pbenwo / (DataImp$pbenwo + DataImp$pmetawo)
round(auc.nonpara.mw(DataImp[DataImp$ADNEXgroups == "Benign", "pp4wo"], DataImp[DataImp$ADNEXgroups == "Metastatic", "pp4wo"], method = "pepe"), 2)
# AUC = 0.92

# Pair 5: borderline vs stage I
DataImp$pp5wo <- DataImp$pborwo / (DataImp$pborwo + DataImp$pst1wo)
DataImp$binbo <- ifelse(DataImp$ADNEXgroups == "Borderline", 0, 1)
round(auc.nonpara.mw(DataImp[DataImp$ADNEXgroups == "Borderline", "pp5wo"], DataImp[DataImp$ADNEXgroups == "Stage I invasive", "pp5wo"], method = "pepe"), 2)
# AUC = 0.72

# Pair 6: borderline vs stage II-IV
DataImp$pp6wo <- DataImp$pborwo / (DataImp$pborwo + DataImp$pst2_4wo)
round(auc.nonpara.mw(DataImp[DataImp$ADNEXgroups == "Borderline", "pp6wo"], DataImp[DataImp$ADNEXgroups == "Stage II-IV invasive", "pp6wo"], method = "pepe"), 2)
# AUC = 0.92

# Pair 7: borderline vs metastatic
DataImp$pp7wo <- DataImp$pborwo / (DataImp$pborwo + DataImp$pmetawo)
round(auc.nonpara.mw(DataImp[DataImp$ADNEXgroups == "Borderline", "pp7wo"], DataImp[DataImp$ADNEXgroups == "Metastatic", "pp7wo"], method = "pepe"), 2)
# AUC = 0.88

# Pair 8: stage I vs stage II-IV
DataImp$pp8wo <- DataImp$pst1wo / (DataImp$pst1wo + DataImp$pst2_4wo)
DataImp$bini <- ifelse(DataImp$ADNEXgroups == "Stage I invasive", 0, 1)
round(auc.nonpara.mw(DataImp[DataImp$ADNEXgroups == "Stage I invasive", "pp8wo"], DataImp[DataImp$ADNEXgroups == "Stage II-IV invasive", "pp8wo"], method = "pepe"), 2)
# AUC = 0.81

# Pair 9: stage I vs metastatic
DataImp$pp9wo <- DataImp$pst1wo / (DataImp$pst1wo + DataImp$pmetawo)
round(auc.nonpara.mw(DataImp[DataImp$ADNEXgroups == "Stage I invasive", "pp9wo"], DataImp[DataImp$ADNEXgroups == "Metastatic", "pp9wo"], method = "pepe"), 2)
# AUC = 0.75

# Pair 10: stage II-IV vs metastatic
DataImp$pp10wo <- DataImp$pst2_4wo / (DataImp$pst2_4wo + DataImp$pmetawo)
DataImp$bin2 <- ifelse(DataImp$ADNEXgroups == "Stage II-IV invasive", 0, 1)
round(auc.nonpara.mw(DataImp[DataImp$ADNEXgroups == "Stage II-IV invasive", "pp10wo"], DataImp[DataImp$ADNEXgroups == "Metastatic", "pp10wo"], method = "pepe"), 2)
# AUC = 0.62

# Polytomous Discrimination Index (PDI)
library(mcca)

DataImp$Groups <- ifelse(DataImp$ADNEXgroups == "Benign", '1',
                         ifelse(DataImp$ADNEXgroups == "Borderline", '2',
                         ifelse(DataImp$ADNEXgroups == "Stage I invasive", '3',
                         ifelse(DataImp$ADNEXgroups == "Stage II-IV invasive", '4', '5'))))


y = outcome = DataImp$Groups # y
d = response = as.matrix(DataImp[, c('pbenwo', 'pborwo', 'pst1wo', 'pst2_4wo', 'pmetawo')]) # d

PDI <- pdI.extented(y = outcome, d = response, k = 5, method = "prob") # 0.49
PDI.se <- ests.extended(y = outcome, d = response, acc = "pdi", level = 0.95, method="prob", k = 5, B = 250, balance = TRUE)$se

LL <- PDI - 1.96*PDI.se
round(LL, 2) # 0.46
UL <- PDI + 1.96*PDI.se
round(UL, 2) # 0.52


#------------------------------#
##### 4. Subgroup analysis #####
#------------------------------#


##### 4.1 Premenopausal patients #####

PreCC <- DataCC[DataCC$Postmenopausal2 == "no",]
PreImp <- DataImp[DataImp$Postmenopausal2 == "no",]

DataPre <- PreCC

table(DataPre$CenterAUC, DataPre$OutcomeBin)

##### 4.1.1 Discrimination between benign and malignant tumors #####

##### ADNEX with CA125 #####

## AUC
AUC.ADNEXwPre <- AUC.IOTA(pred = pmalw, outcome = OutcomeBin, center = CenterAUC, data = DataPre, method.MA = "BAYES", titleGraph = "AUC per center")

AUC.ADNEXwPre$Plot[9, "RRauc"] <- "   "
AUC.ADNEXwPre$Plot[10, "RRauc"] <- "   "
AUC.ADNEXwPre$Plot[11, "RRcenter"] <- "AUC (95% CI)"
AUC.ADNEXwPre$Plot[12, "RRauc"] <- "        (0.40 to 0.99)"
AUC.ADNEXwPre$Plot[12, "RRprev"] <- "   "

forestplot(AUC.ADNEXwPre$Plot,
           align = c("l", "c", "c"),
           mean = AUC.ADNEXwPre$dataPlot$AUC,
           lower = AUC.ADNEXwPre$dataPlot$LL,
           upper = AUC.ADNEXwPre$dataPlot$UL,
           is.summary = c(TRUE, FALSE, rep(FALSE, length(AUC.ADNEXwPre$IncludedCenters)), FALSE, FALSE, TRUE, TRUE),
           xlab = "AUC", 
           xticks = c(0.60, 0.70, 0.80, 0.9, 1),
           clip = c(0.61, 1),
           txt_gp = fpTxtGp(summary = gpar(cex = 1.55), xlab = gpar(cex = 1.5), label = gpar(cex = 1.5),
                            ticks = gpar(cex = 1.5), title = gpar(cex = 1.75)),
           graphwidth = unit(6, "cm"),
           col = fpColors(box = "gray", lines = "black", summary = "black"),
           fn.ci_sum = c(as.list(rep("fpDrawSummaryCI", length(AUC.ADNEXwPre$IncludedCenters) + 4)),
                         fpDrawSummaryCI,
                         fpDrawBarCI),
           zero = AUC.ADNEXwPre$Performance$AUC)
# 900 x 450

## Sensitivity and specificity
SensSpec(pred = pmalw, outcome = OutcomeBin, threshold = 0.01,  center = CenterAUC, data = DataPre)$OverallPer
SensSpec(pred = pmalw, outcome = OutcomeBin, threshold = 0.03,  center = CenterAUC, data = DataPre)$OverallPer
SensSpec(pred = pmalw, outcome = OutcomeBin, threshold = 0.05,  center = CenterAUC, data = DataPre)$OverallPer
SensSpec(pred = pmalw, outcome = OutcomeBin, threshold = 0.10,  center = CenterAUC, data = DataPre)$OverallPer
SensSpec(pred = pmalw, outcome = OutcomeBin, threshold = 0.15,  center = CenterAUC, data = DataPre)$OverallPer
SensSpec(pred = pmalw, outcome = OutcomeBin, threshold = 0.20,  center = CenterAUC, data = DataPre)$OverallPer
SensSpec(pred = pmalw, outcome = OutcomeBin, threshold = 0.25,  center = CenterAUC, data = DataPre)$OverallPer
SensSpec(pred = pmalw, outcome = OutcomeBin, threshold = 0.30,  center = CenterAUC, data = DataPre)$OverallPer
SensSpec(pred = pmalw, outcome = OutcomeBin, threshold = 0.40,  center = CenterAUC, data = DataPre)$OverallPer
SensSpec(pred = pmalw, outcome = OutcomeBin, threshold = 0.50,  center = CenterAUC, data = DataPre)$OverallPer

## Sensitivity at fixed specificity
FixSpec(pred = pmalw, outcome = OutcomeBin, Specificity = 0.90,  center = CenterAUC, data = DataPre)$OverallPer

## Specificity at fixed sensitivity
FixSens(pred = pmalw, outcome = OutcomeBin, Sensitivity = 0.90,  center = CenterAUC, data = DataPre)$OverallPer

##### ADNEX without CA125 #####

## AUC
AUC.ADNEXwoPre <- AUC.IOTA(pred = pmalwo, outcome = OutcomeBin, center = CenterAUC, data = DataPre, method.MA = "BAYES", titleGraph = "AUC per center")

AUC.ADNEXwoPre$Plot[9, "RRauc"] <- "   "
AUC.ADNEXwoPre$Plot[10, "RRauc"] <- "   "
AUC.ADNEXwoPre$Plot[11, "RRcenter"] <- "AUC (95% CI)"
AUC.ADNEXwoPre$Plot[12, "RRauc"] <- "        (0.46 to 0.99)"
AUC.ADNEXwoPre$Plot[12, "RRprev"] <- "   "

forestplot(AUC.ADNEXwoPre$Plot,
           align = c("l", "c", "c"),
           mean = AUC.ADNEXwoPre$dataPlot$AUC,
           lower = AUC.ADNEXwoPre$dataPlot$LL,
           upper = AUC.ADNEXwoPre$dataPlot$UL,
           is.summary = c(TRUE, FALSE, rep(FALSE, length(AUC.ADNEXwoPre$IncludedCenters)), FALSE, FALSE, TRUE, TRUE),
           xlab = "AUC",
           xticks = c(0.60, 0.70, 0.80, 0.9, 1),
           clip = c(0.61, 1),
           txt_gp = fpTxtGp(summary = gpar(cex = 1.55), xlab = gpar(cex = 1.5), label = gpar(cex = 1.5),
                            ticks = gpar(cex = 1.5), title = gpar(cex = 1.75)),
           graphwidth = unit(6, "cm"),
           col = fpColors(box = "gray", lines = "black", summary = "black"),
           fn.ci_sum = c(as.list(rep("fpDrawSummaryCI", length(AUC.ADNEXwoPre$IncludedCenters) + 4)),
                         fpDrawSummaryCI,
                         fpDrawBarCI),
           zero = AUC.ADNEXwoPre$Performance$AUC)
# 900 x 450

## Sensitivity and specificity
SensSpec(pred = pmalwo, outcome = OutcomeBin, threshold = 0.01, center = CenterAUC, data = DataPre)$OverallPer
SensSpec(pred = pmalwo, outcome = OutcomeBin, threshold = 0.03, center = CenterAUC, data = DataPre)$OverallPer
SensSpec(pred = pmalwo, outcome = OutcomeBin, threshold = 0.05, center = CenterAUC, data = DataPre)$OverallPer
SensSpec(pred = pmalwo, outcome = OutcomeBin, threshold = 0.10, center = CenterAUC, data = DataPre)$OverallPer
SensSpec(pred = pmalwo, outcome = OutcomeBin, threshold = 0.15, center = CenterAUC, data = DataPre)$OverallPer
SensSpec(pred = pmalwo, outcome = OutcomeBin, threshold = 0.20, center = CenterAUC, data = DataPre)$OverallPer
SensSpec(pred = pmalwo, outcome = OutcomeBin, threshold = 0.25, center = CenterAUC, data = DataPre)$OverallPer
SensSpec(pred = pmalwo, outcome = OutcomeBin, threshold = 0.30, center = CenterAUC, data = DataPre)$OverallPer
SensSpec(pred = pmalwo, outcome = OutcomeBin, threshold = 0.40, center = CenterAUC, data = DataPre)$OverallPer
SensSpec(pred = pmalwo, outcome = OutcomeBin, threshold = 0.50, center = CenterAUC, data = DataPre)$OverallPer

## Sensitivity at fixed specificity
FixSpec(pred = pmalwo, outcome = OutcomeBin, Specificity = 0.90, center = CenterAUC, data = DataPre)$OverallPer

## Specificity at fixed sensitivity
FixSens(pred = pmalwo, outcome = OutcomeBin, Sensitivity = 0.90, center = CenterAUC, data = DataPre)$OverallPer

##### ROMA #####

## AUC
AUC.ROMApre <- AUC.IOTA(pred = ROMA, outcome = OutcomeBin, center = CenterAUC, data = DataPre, method.MA = "BAYES", titleGraph = "AUC per center")

AUC.ROMApre$Plot[9, "RRauc"] <- "   "
AUC.ROMApre$Plot[10, "RRauc"] <- "   "
AUC.ROMApre$Plot[11, "RRcenter"] <- "AUC (95% CI)"
AUC.ROMApre$Plot[12, "RRauc"] <- "        (0.42 to 0.94)"
AUC.ROMApre$Plot[12, "RRprev"] <- "   "

forestplot(AUC.ROMApre$Plot,
           align = c("l", "c", "c"),
           mean = AUC.ROMApre$dataPlot$AUC,
           lower = AUC.ROMApre$dataPlot$LL,
           upper = AUC.ROMApre$dataPlot$UL,
           is.summary = c(TRUE, FALSE, rep(FALSE, length(AUC.ROMApre$IncludedCenters)), FALSE, FALSE, TRUE, TRUE),
           xlab = "AUC", 
           xticks = c(0.60, 0.70, 0.80, 0.9, 1),
           clip = c(0.61, 1),
           txt_gp = fpTxtGp(summary = gpar(cex = 1.55), xlab = gpar(cex = 1.5), label = gpar(cex = 1.5),
                            ticks = gpar(cex = 1.5), title = gpar(cex = 1.75)),
           graphwidth = unit(6, "cm"),
           col = fpColors(box = "gray", lines = "black", summary = "black"),
           fn.ci_sum = c(as.list(rep("fpDrawSummaryCI", length(AUC.ROMApre$IncludedCenters) + 4)),
                         fpDrawSummaryCI,
                         fpDrawBarCI),
           zero = AUC.ROMApre$Performance$AUC)
# 900 x 450

## Sensitivity and specificity
SensSpec(pred = ROMA, outcome = OutcomeBin, threshold = 0.01, center = CenterAUC, data = DataPre)$OverallPer
SensSpec(pred = ROMA, outcome = OutcomeBin, threshold = 0.03, center = CenterAUC, data = DataPre)$OverallPer
SensSpec(pred = ROMA, outcome = OutcomeBin, threshold = 0.05, center = CenterAUC, data = DataPre)$OverallPer
SensSpec(pred = ROMA, outcome = OutcomeBin, threshold = 0.10, center = CenterAUC, data = DataPre)$OverallPer
SensSpec(pred = ROMA, outcome = OutcomeBin, threshold = 0.15, center = CenterAUC, data = DataPre)$OverallPer
SensSpec(pred = ROMA, outcome = OutcomeBin, threshold = 0.20, center = CenterAUC, data = DataPre)$OverallPer
SensSpec(pred = ROMA, outcome = OutcomeBin, threshold = 0.25, center = CenterAUC, data = DataPre)$OverallPer
SensSpec(pred = ROMA, outcome = OutcomeBin, threshold = 0.30, center = CenterAUC, data = DataPre)$OverallPer
SensSpec(pred = ROMA, outcome = OutcomeBin, threshold = 0.40, center = CenterAUC, data = DataPre)$OverallPer
SensSpec(pred = ROMA, outcome = OutcomeBin, threshold = 0.50, center = CenterAUC, data = DataPre)$OverallPer

## Sensitivity at fixed specificity
FixSpec(pred = ROMA, outcome = OutcomeBin, Specificity = 0.90, center = CenterAUC, data = DataPre)$OverallPer

## Specificity at fixed sensitivity
FixSens(pred = ROMA, outcome = OutcomeBin, Sensitivity = 0.90, center = CenterAUC, data = DataPre)$OverallPer

##### Summary plot for AUC #####

NA.forest <- AUC.ROMApre
NA.forest <- NA

Summary.AUC <- rbind(NA.forest, AUC.ROMApre$Performance[1,], AUC.ADNEXwoPre$Performance[1,], AUC.ADNEXwPre$Performance[1,])
Summary.AUC$Model <- c('', 'ROMA', 'ADNEX w/o CA125', 'ADNEX w/ CA125')

Summary.PI <- rbind(NA.forest, AUC.ROMApre$Performance[2,], AUC.ADNEXwoPre$Performance[2,], AUC.ADNEXwPre$Performance[2,])
Summary.PI$Model <- c('', 'ROMA', 'ADNEX w/o CA125', 'ADNEX w/ CA125')

tabletext <- cbind(
  c('Model', 'ROMA', 'ADNEX w/o CA125', 'ADNEX w/ CA125'),
  c('AUC (95% CI)', 
    paste(format(round(AUC.ROMApre$Performance$AUC[1], 2), nsmall = 2), " (", format(round(AUC.ROMApre$Performance$LL[1], 2), nsmall = 2), "; ", format(round(AUC.ROMApre$Performance$UL[1], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(AUC.ADNEXwoPre$Performance$AUC[1], 2), nsmall = 2), " (", format(round(AUC.ADNEXwoPre$Performance$LL[1], 2), nsmall = 2), "; ", format(round(AUC.ADNEXwoPre$Performance$UL[1], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(AUC.ADNEXwPre$Performance$AUC[1], 2), nsmall = 2), " (", format(round(AUC.ADNEXwPre$Performance$LL[1], 2), nsmall = 2), "; ", format(round(AUC.ADNEXwPre$Performance$UL[1], 2), nsmall = 2), ")", sep = "")),
  c('95% PI', 
    paste("(", format(round(AUC.ROMApre$Performance$LL[2], 2), nsmall = 2), "; ", format(round(AUC.ROMApre$Performance$UL[2], 2), nsmall = 2), ")", sep = ""),
    paste("(", format(round(AUC.ADNEXwoPre$Performance$LL[2], 2), nsmall = 2), "; ", format(round(AUC.ADNEXwoPre$Performance$UL[2], 2), nsmall = 2), ")", sep = ""),
    paste("(", format(round(AUC.ADNEXwPre$Performance$LL[2], 2), nsmall = 2), "; ", format(round(AUC.ADNEXwPre$Performance$UL[2], 2), nsmall = 2), ")", sep = ""))
)

forestplot(labeltext = tabletext,
           title = "",
           mean = round(Summary.AUC$AUC, 3),
           lower = round(Summary.AUC$LL, 3),
           upper = round(Summary.AUC$UL, 3),
           is.summary = c(FALSE, TRUE, TRUE, TRUE),
           xlab = "AUC (95% CI)",
           boxsize = .5,
           txt_gp = fpTxtGp(summary = gpar(cex = 1.55, fontface = "plain"), xlab = gpar(cex = 1.5, fontface = "plain"), label = gpar(cex = 1.5, fontface = "bold"),
                            ticks = gpar(cex = 1.5, fontface = "plain"), title = gpar(cex = 1.75)),
           graphwidth = unit(9, "cm"),
           graph.pos = 3,
           xticks = c(0.60, 0.70, 0.80, 0.90, 1), xlog = TRUE, clip = c(0.60, 1))
# 1050 x 400


##### 4.2 Postmenopausal patients #####

PostCC <- DataCC[DataCC$Postmenopausal2 == "yes",]
PostImp <- DataImp[DataImp$Postmenopausal2 == "yes",]

DataPost <- PostCC

table(DataPost$CenterAUC, DataPost$OutcomeBin)

##### 4.2.1 Discrimination between benign and malignant tumors #####

##### ADNEX with CA125 #####

## AUC
AUC.ADNEXwPost <- AUC.IOTA(pred = pmalw, outcome = OutcomeBin, center = CenterAUC, data = DataPost, method.MA = "BAYES", titleGraph = "AUC per center")

AUC.ADNEXwPost$Plot[9, "RRauc"] <- "   "
AUC.ADNEXwPost$Plot[10, "RRauc"] <- "   "
AUC.ADNEXwPost$Plot[11, "RRcenter"] <- "AUC (95% CI)"
AUC.ADNEXwPost$Plot[12, "RRauc"] <- "        (0.77 to 0.98)"
AUC.ADNEXwPost$Plot[12, "RRprev"] <- "   "

forestplot(AUC.ADNEXwPost$Plot,
           align = c("l", "c", "c"),
           mean = AUC.ADNEXwPost$dataPlot$AUC,
           lower = AUC.ADNEXwPost$dataPlot$LL,
           upper = AUC.ADNEXwPost$dataPlot$UL,
           is.summary = c(TRUE, FALSE, rep(FALSE, length(AUC.ADNEXwPost$IncludedCenters)), FALSE, FALSE, TRUE, TRUE),
           xlab = "AUC",
           xticks = c(0.70, 0.80, 0.9, 1),
           clip = c(0.71, 1),
           txt_gp = fpTxtGp(summary = gpar(cex = 1.55), xlab = gpar(cex = 1.5), label = gpar(cex = 1.5),
                            ticks = gpar(cex = 1.5), title = gpar(cex = 1.75)),
           graphwidth = unit(6, "cm"),
           col = fpColors(box = "gray", lines = "black", summary = "black"),
           fn.ci_sum = c(as.list(rep("fpDrawSummaryCI", length(AUC.ADNEXwPost$IncludedCenters) + 4)),
                         fpDrawSummaryCI,
                         fpDrawBarCI),
           zero = AUC.ADNEXwPost$Performance$AUC)
# 900 x 450

## Sensitivity and specificity
SensSpec(pred = pmalw, outcome = OutcomeBin, threshold = 0.01,  center = CenterAUC, data = DataPost)$OverallPer
SensSpec(pred = pmalw, outcome = OutcomeBin, threshold = 0.03,  center = CenterAUC, data = DataPost)$OverallPer
SensSpec(pred = pmalw, outcome = OutcomeBin, threshold = 0.05,  center = CenterAUC, data = DataPost)$OverallPer
SensSpec(pred = pmalw, outcome = OutcomeBin, threshold = 0.10,  center = CenterAUC, data = DataPost)$OverallPer
SensSpec(pred = pmalw, outcome = OutcomeBin, threshold = 0.15,  center = CenterAUC, data = DataPost)$OverallPer
SensSpec(pred = pmalw, outcome = OutcomeBin, threshold = 0.20,  center = CenterAUC, data = DataPost)$OverallPer
SensSpec(pred = pmalw, outcome = OutcomeBin, threshold = 0.25,  center = CenterAUC, data = DataPost)$OverallPer
SensSpec(pred = pmalw, outcome = OutcomeBin, threshold = 0.30,  center = CenterAUC, data = DataPost)$OverallPer
SensSpec(pred = pmalw, outcome = OutcomeBin, threshold = 0.40,  center = CenterAUC, data = DataPost)$OverallPer
SensSpec(pred = pmalw, outcome = OutcomeBin, threshold = 0.50,  center = CenterAUC, data = DataPost)$OverallPer

## Sensitivity at fixed specificity
FixSpec(pred = pmalw, outcome = OutcomeBin, Specificity = 0.90,  center = CenterAUC, data = DataPost)$OverallPer

## Specificity at fixed sensitivity
FixSens(pred = pmalw, outcome = OutcomeBin, Sensitivity = 0.90,  center = CenterAUC, data = DataPost)$OverallPer

##### ADNEX without CA125 #####

## AUC
AUC.ADNEXwoPost <- AUC.IOTA(pred = pmalwo, outcome = OutcomeBin, center = CenterAUC, data = DataPost, method.MA = "BAYES", titleGraph = "AUC per center")

AUC.ADNEXwoPost$Plot[9, "RRauc"] <- "   "
AUC.ADNEXwoPost$Plot[10, "RRauc"] <- "   "
AUC.ADNEXwoPost$Plot[11, "RRcenter"] <- "AUC (95% CI)"
AUC.ADNEXwoPost$Plot[12, "RRauc"] <- "        (0.71 to 0.97)"
AUC.ADNEXwoPost$Plot[12, "RRprev"] <- "   "

forestplot(AUC.ADNEXwoPost$Plot,
           align = c("l", "c", "c"),
           mean = AUC.ADNEXwoPost$dataPlot$AUC,
           lower = AUC.ADNEXwoPost$dataPlot$LL,
           upper = AUC.ADNEXwoPost$dataPlot$UL,
           is.summary = c(TRUE, FALSE, rep(FALSE, length(AUC.ADNEXwoPost$IncludedCenters)), FALSE, FALSE, TRUE, TRUE),
           xlab = "AUC", 
           xticks = c(0.70, 0.80, 0.9, 1),
           clip = c(0.71, 1),
           txt_gp = fpTxtGp(summary = gpar(cex = 1.55), xlab = gpar(cex = 1.5), label = gpar(cex = 1.5),
                            ticks = gpar(cex = 1.5), title = gpar(cex = 1.75)),
           graphwidth = unit(6, "cm"),
           col = fpColors(box = "gray", lines = "black", summary = "black"),
           fn.ci_sum = c(as.list(rep("fpDrawSummaryCI", length(AUC.ADNEXwoPost$IncludedCenters) + 4)),
                         fpDrawSummaryCI,
                         fpDrawBarCI),
           zero = AUC.ADNEXwoPost$Performance$AUC)
# 900 x 450

## Sensitivity and specificity
SensSpec(pred = pmalwo, outcome = OutcomeBin, threshold = 0.01, center = CenterAUC, data = DataPost)$OverallPer
SensSpec(pred = pmalwo, outcome = OutcomeBin, threshold = 0.03, center = CenterAUC, data = DataPost)$OverallPer
SensSpec(pred = pmalwo, outcome = OutcomeBin, threshold = 0.05, center = CenterAUC, data = DataPost)$OverallPer
SensSpec(pred = pmalwo, outcome = OutcomeBin, threshold = 0.10, center = CenterAUC, data = DataPost)$OverallPer
SensSpec(pred = pmalwo, outcome = OutcomeBin, threshold = 0.15, center = CenterAUC, data = DataPost)$OverallPer
SensSpec(pred = pmalwo, outcome = OutcomeBin, threshold = 0.20, center = CenterAUC, data = DataPost)$OverallPer
SensSpec(pred = pmalwo, outcome = OutcomeBin, threshold = 0.25, center = CenterAUC, data = DataPost)$OverallPer
SensSpec(pred = pmalwo, outcome = OutcomeBin, threshold = 0.30, center = CenterAUC, data = DataPost)$OverallPer
SensSpec(pred = pmalwo, outcome = OutcomeBin, threshold = 0.40, center = CenterAUC, data = DataPost)$OverallPer
SensSpec(pred = pmalwo, outcome = OutcomeBin, threshold = 0.50, center = CenterAUC, data = DataPost)$OverallPer

## Sensitivity at fixed specificity
FixSpec(pred = pmalwo, outcome = OutcomeBin, Specificity = 0.90, center = CenterAUC, data = DataPost)$OverallPer

## Specificity at fixed sensitivity
FixSens(pred = pmalwo, outcome = OutcomeBin, Sensitivity = 0.90, center = CenterAUC, data = DataPost)$OverallPer

##### ROMA #####

## AUC
AUC.ROMApost <- AUC.IOTA(pred = ROMA, outcome = OutcomeBin, center = CenterAUC, data = DataPost, method.MA = "BAYES", titleGraph = "AUC per center")

AUC.ROMApost$Plot[9, "RRauc"] <- "   "
AUC.ROMApost$Plot[10, "RRauc"] <- "   "
AUC.ROMApost$Plot[11, "RRcenter"] <- "AUC (95% CI)"
AUC.ROMApost$Plot[12, "RRauc"] <- "        (0.74 to 0.96)"
AUC.ROMApost$Plot[12, "RRprev"] <- "   "

forestplot(AUC.ROMApost$Plot,
           align = c("l", "c", "c"),
           mean = AUC.ROMApost$dataPlot$AUC,
           lower = AUC.ROMApost$dataPlot$LL,
           upper = AUC.ROMApost$dataPlot$UL,
           is.summary = c(TRUE, FALSE, rep(FALSE, length(AUC.ROMApost$IncludedCenters)), FALSE, FALSE, TRUE, TRUE),
           xlab = "AUC", 
           xticks = c(0.70, 0.80, 0.9, 1),
           clip = c(0.71, 1),
           txt_gp = fpTxtGp(summary = gpar(cex = 1.55), xlab = gpar(cex = 1.5), label = gpar(cex = 1.5),
                            ticks = gpar(cex = 1.5), title = gpar(cex = 1.75)),
           graphwidth = unit(6, "cm"),
           col = fpColors(box = "gray", lines = "black", summary = "black"),
           fn.ci_sum = c(as.list(rep("fpDrawSummaryCI", length(AUC.ROMApost$IncludedCenters) + 4)),
                         fpDrawSummaryCI,
                         fpDrawBarCI),
           zero = AUC.ROMApost$Performance$AUC)
# 900 x 450

## Sensitivity and specificity
SensSpec(pred = ROMA, outcome = OutcomeBin, threshold = 0.01, center = CenterAUC, data = DataPost)$OverallPer
SensSpec(pred = ROMA, outcome = OutcomeBin, threshold = 0.03, center = CenterAUC, data = DataPost)$OverallPer
SensSpec(pred = ROMA, outcome = OutcomeBin, threshold = 0.05, center = CenterAUC, data = DataPost)$OverallPer
SensSpec(pred = ROMA, outcome = OutcomeBin, threshold = 0.10, center = CenterAUC, data = DataPost)$OverallPer
SensSpec(pred = ROMA, outcome = OutcomeBin, threshold = 0.15, center = CenterAUC, data = DataPost)$OverallPer
SensSpec(pred = ROMA, outcome = OutcomeBin, threshold = 0.20, center = CenterAUC, data = DataPost)$OverallPer
SensSpec(pred = ROMA, outcome = OutcomeBin, threshold = 0.25, center = CenterAUC, data = DataPost)$OverallPer
SensSpec(pred = ROMA, outcome = OutcomeBin, threshold = 0.30, center = CenterAUC, data = DataPost)$OverallPer
SensSpec(pred = ROMA, outcome = OutcomeBin, threshold = 0.40, center = CenterAUC, data = DataPost)$OverallPer
SensSpec(pred = ROMA, outcome = OutcomeBin, threshold = 0.50, center = CenterAUC, data = DataPost)$OverallPer

## Sensitivity at fixed specificity
FixSpec(pred = ROMA, outcome = OutcomeBin, Specificity = 0.90, center = CenterAUC, data = DataPost)$OverallPer

## Specificity at fixed sensitivity
FixSens(pred = ROMA, outcome = OutcomeBin, Sensitivity = 0.90, center = CenterAUC, data = DataPost)$OverallPer

##### Summary plot for AUC #####

NA.forest <- AUC.ROMApost
NA.forest <- NA

Summary.AUC <- rbind(NA.forest, AUC.ROMApost$Performance[1,], AUC.ADNEXwoPost$Performance[1,], AUC.ADNEXwPost$Performance[1,])
Summary.AUC$Model <- c('', 'ROMA', 'ADNEX w/o CA125', 'ADNEX w/ CA125')

Summary.PI <- rbind(NA.forest, AUC.ROMApost$Performance[1,], AUC.ADNEXwoPost$Performance[1,], AUC.ADNEXwPost$Performance[1,])
Summary.PI$Model <- c('', 'ROMA', 'ADNEX w/o CA125', 'ADNEX w/ CA125')

tabletext <- cbind(
  c('Model', 'ROMA', 'ADNEX w/o CA125', 'ADNEX w/ CA125'),
  c('AUC (95% CI)', 
    paste(format(round(AUC.ROMApost$Performance$AUC[1], 2), nsmall = 2), " (", format(round(AUC.ROMApost$Performance$LL[1], 2), nsmall = 2), "; ", format(round(AUC.ROMApost$Performance$UL[1], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(AUC.ADNEXwoPost$Performance$AUC[1], 2), nsmall = 2), " (", format(round(AUC.ADNEXwoPost$Performance$LL[1], 2), nsmall = 2), "; ", format(round(AUC.ADNEXwoPost$Performance$UL[1], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(AUC.ADNEXwPost$Performance$AUC[1], 2), nsmall = 2), " (", format(round(AUC.ADNEXwPost$Performance$LL[1], 2), nsmall = 2), "; ", format(round(AUC.ADNEXwPost$Performance$UL[1], 2), nsmall = 2), ")", sep = "")),
  c('95% PI', 
    paste("(", format(round(AUC.ROMApost$Performance$LL[2], 2), nsmall = 2), "; ", format(round(AUC.ROMApost$Performance$UL[2], 2), nsmall = 2), ")", sep = ""),
    paste("(", format(round(AUC.ADNEXwoPost$Performance$LL[2], 2), nsmall = 2), "; ", format(round(AUC.ADNEXwoPost$Performance$UL[2], 2), nsmall = 2), ")", sep = ""),
    paste("(", format(round(AUC.ADNEXwPost$Performance$LL[2], 2), nsmall = 2), "; ", format(round(AUC.ADNEXwPost$Performance$UL[2], 2), nsmall = 2), ")", sep = ""))
)

forestplot(labeltext = tabletext,
           title = "",
           mean = round(Summary.AUC$AUC, 3),
           lower = round(Summary.AUC$LL, 3),
           upper = round(Summary.AUC$UL, 3),
           is.summary = c(FALSE, TRUE, TRUE, TRUE),
           xlab = "AUC (95% CI)",
           boxsize = .5,
           txt_gp = fpTxtGp(summary = gpar(cex = 1.55, fontface = "plain"), xlab = gpar(cex = 1.5, fontface = "plain"), label = gpar(cex = 1.5, fontface = "bold"),
                            ticks = gpar(cex = 1.5, fontface = "plain"), title = gpar(cex = 1.75)),
           graphwidth = unit(9, "cm"),
           graph.pos = 3,
           xticks = c(0.80, 0.85, 0.9, 0.95, 1), xlog = TRUE, clip = c(0.80, 1))
# 1050 x 400


##### 4.3 Summary plot for subgroup analyses #####

NA.forest <- AUC.ROMApre
NA.forest <- NA

Summary.AUC <- rbind(NA.forest, AUC.ROMApre$Performance[1,], AUC.ADNEXwoPre$Performance[1,], AUC.ADNEXwPre$Performance[1,], NA.forest, NA.forest, AUC.ROMApost$Performance[1,], AUC.ADNEXwoPost$Performance[1,], AUC.ADNEXwPost$Performance[1,])
Summary.AUC$Model <- c('', 'ROMA', 'ADNEX without CA125', 'ADNEX with CA125', '', '', 'ROMA', 'ADNEX without CA125', 'ADNEX with CA125')

Summary.PI <- rbind(NA.forest, AUC.ROMApre$Performance[2,], AUC.ADNEXwoPre$Performance[2,], AUC.ADNEXwPre$Performance[2,], NA.forest, NA.forest, AUC.ROMApost$Performance[1,], AUC.ADNEXwoPost$Performance[1,], AUC.ADNEXwPost$Performance[1,])
Summary.PI$Model <- c('', 'ROMA', 'ADNEX without CA125', 'ADNEX with CA125', '', '', 'ROMA', 'ADNEX without CA125', 'ADNEX with CA125')

tabletext <- cbind(
  c('Premenopausal patients', 'ROMA', 'ADNEX without CA125', 'ADNEX with CA125', '', 'Postmenopausal patients', 'ROMA', 'ADNEX without CA125', 'ADNEX with CA125'),
  c('AUC (95% CI)', 
    paste(format(round(AUC.ROMApre$Performance$AUC[1], 2), nsmall = 2), " (", format(round(AUC.ROMApre$Performance$LL[1], 2), nsmall = 2), "; ", format(round(AUC.ROMApre$Performance$UL[1], 2), nsmall = 2), ")", sep = ""),
    paste("0.90 (", format(round(AUC.ADNEXwoPre$Performance$LL[1], 2), nsmall = 2), "; ", format(round(AUC.ADNEXwoPre$Performance$UL[1], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(AUC.ADNEXwPre$Performance$AUC[1], 2), nsmall = 2), " (", format(round(AUC.ADNEXwPre$Performance$LL[1], 2), nsmall = 2), "; ", format(round(AUC.ADNEXwPre$Performance$UL[1], 2), nsmall = 2), ")", sep = ""),
    '',
    'AUC (95% CI)', 
    paste(format(round(AUC.ROMApost$Performance$AUC[1], 2), nsmall = 2), " (", format(round(AUC.ROMApost$Performance$LL[1], 2), nsmall = 2), "; ", format(round(AUC.ROMApost$Performance$UL[1], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(AUC.ADNEXwoPost$Performance$AUC[1], 2), nsmall = 2), " (", format(round(AUC.ADNEXwoPost$Performance$LL[1], 2), nsmall = 2), "; ", format(round(AUC.ADNEXwoPost$Performance$UL[1], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(AUC.ADNEXwPost$Performance$AUC[1], 2), nsmall = 2), " (", format(round(AUC.ADNEXwPost$Performance$LL[1], 2), nsmall = 2), "; ", format(round(AUC.ADNEXwPost$Performance$UL[1], 2), nsmall = 2), ")", sep = "")),
  c('95% PI', 
    paste("(", format(round(AUC.ROMApre$Performance$LL[2], 2), nsmall = 2), "; ", format(round(AUC.ROMApre$Performance$UL[2], 2), nsmall = 2), ")", sep = ""),
    paste("(", format(round(AUC.ADNEXwoPre$Performance$LL[2], 2), nsmall = 2), "; ", format(round(AUC.ADNEXwoPre$Performance$UL[2], 2), nsmall = 2), ")", sep = ""),
    paste("(", format(round(AUC.ADNEXwPre$Performance$LL[2], 2), nsmall = 2), "; ", format(round(AUC.ADNEXwPre$Performance$UL[2], 2), nsmall = 2), ")", sep = ""),
    '',
    '95% PI', 
    paste("(", format(round(AUC.ROMApost$Performance$LL[2], 2), nsmall = 2), "; ", format(round(AUC.ROMApost$Performance$UL[2], 2), nsmall = 2), ")", sep = ""),
    paste("(", format(round(AUC.ADNEXwoPost$Performance$LL[2], 2), nsmall = 2), "; ", format(round(AUC.ADNEXwoPost$Performance$UL[2], 2), nsmall = 2), ")", sep = ""),
    paste("(", format(round(AUC.ADNEXwPost$Performance$LL[2], 2), nsmall = 2), "; ", format(round(AUC.ADNEXwPost$Performance$UL[2], 2), nsmall = 2), ")", sep = ""))
)

forestplot(labeltext = tabletext,
           title = "",
           mean = round(Summary.AUC$AUC, 3),
           lower = round(Summary.AUC$LL, 3),
           upper = round(Summary.AUC$UL, 3),
           is.summary = c(FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, TRUE),
           xlab = "AUC (95% CI)",
           boxsize = .5,
           txt_gp = fpTxtGp(summary = gpar(cex = 1.55, fontface = "plain"), xlab = gpar(cex = 1.5, fontface = "plain"), label = gpar(cex = 1.5, fontface = "bold"),
                            ticks = gpar(cex = 1.5, fontface = "plain"), title = gpar(cex = 1.75)),
           graphwidth = unit(9, "cm"),
           graph.pos = 3,
           xticks = c(0.60, 0.70, 0.80, 0.90, 1), xlog = TRUE, clip = c(0.60, 1))
# 1100 x 775


save.image(file = "C:/Users/u0123496/OneDrive - KU Leuven/IOTA/ROMA vs ADNEX/ROMA vs ADNEX - 112021.RData") # Save the entire workspace


