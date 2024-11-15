library(dplyr)
source("~/Downloads/SPIRIT/SPIRIT/fun_pval.R")
data = read.csv("Spirit_metadata.csv")
data = data %>% filter(TreatmentAssignment %in% c("metformin", "self-directed"))
treatment = data$TreatmentAssignment
X = factor(ifelse(treatment== "metformin", "metformin", "control"))
#Inflammation markers
Y = list()
Y[[1]] = data$IL_6_pg_ml_UMD_12mth
Y[[2]] = data$hsCRP_UMD_12mth
#SCFAs
M = list()
M[[1]] = data$AceticAcid_6V_serum
M[[2]] = data$PropionicAcid_6V_serum
M[[3]] = data$IsobutyricAcid_6V_serum
M[[4]] = data$ButyricAcid_6V_serum
M[[5]] = data$MethylbutyricAcid_6V_serum
M[[6]] = data$IsovalericAcid_6V_serum
M[[7]] = data$ValericAcid_6V_serum
M[[8]] = data$HexanoicAcid_6V_serum

#confounders
age = data$Age_RZ
race = as.factor(data$Race)
BMI = data$BMI_RZ
sex = as.factor(data$SEX)
trt = data$Yrs_since_treatment
other_covariates = cbind(age, race, BMI, sex, trt)
pvalue_mat = matrix(nrow = 16, ncol = 4)
counter = 1
#set.seed(333)
for(inflam in 1:2)
{
  for(scfa in 1:8)
  {
    outcome = scale(Y[[inflam]])
    med = scale(M[[scfa]])
    pvalue_mat[counter,] = fun.pval(X, med, outcome, other_covariates = other_covariates)
      counter = counter+1
      print(counter)
  }
}

colnames(pvalue_mat) = c("Subbel", "ABtest", "MaxP", "Sobel")
#write.csv(pvalue_mat, "Spirit-data-pvalues.csv")

Treatment = rep("metformin/control", 16)
Outcome = rep(c("IL-6", "hs-CRP"), each = 8)
Mediator = rep(c("Acetic Acid", "Propionic Acid", "Isobutyric Acid", "Butyric Acid", "Methylbutyric Acid", "Isovaleric Acid", "Valeric Acid", "Hexanoic Acid"), times = 2)
pvalue_complete = data.frame(Treatment, Outcome, Mediator, pvalue_mat)


###Repeat everything with Treatment = coach-directed weight loss vs self directed 
data = read.csv("Spirit_metadata.csv")
data = data %>% filter(TreatmentAssignment %in% c("coach-directed", "self-directed"))
treatment = data$TreatmentAssignment
X = factor(ifelse(treatment== "coach-directed", "lifestyle-intervention", "control"))
#Inflammation markers
Y = list()
Y[[1]] = data$IL_6_pg_ml_UMD_12mth
Y[[2]] = data$hsCRP_UMD_12mth
#SCFAs
M = list()
M[[1]] = data$AceticAcid_6V_serum
M[[2]] = data$PropionicAcid_6V_serum
M[[3]] = data$IsobutyricAcid_6V_serum
M[[4]] = data$ButyricAcid_6V_serum
M[[5]] = data$MethylbutyricAcid_6V_serum
M[[6]] = data$IsovalericAcid_6V_serum
M[[7]] = data$ValericAcid_6V_serum
M[[8]] = data$HexanoicAcid_6V_serum
#confounders
age = data$Age_RZ
race = as.factor(data$Race)
BMI = data$BMI_RZ
sex = as.factor(data$SEX)
trt = data$Yrs_since_treatment
other_covariates = cbind(age, race, BMI, sex, trt)
pvalue_mat2 = matrix(nrow = 16, ncol = 4)
counter = 1
#set.seed(333)
for(inflam in 1:2)
{
  for(scfa in 1:8)
  {
    outcome = scale(Y[[inflam]])
    med = scale(M[[scfa]])
    pvalue_mat2[counter,] = fun.pval(X, med, outcome, other_covariates = other_covariates)
    counter = counter+1
    print(counter)
  }
}
colnames(pvalue_mat2) = c("Subbel", "ABtest", "MaxP", "Sobel")
Treatment = rep("lifestyle/control", 16)
Outcome = rep(c("IL-6", "hs-CRP"), each = 8)
Mediator = rep(c("Acetic Acid", "Propionic Acid", "Isobutyric Acid", "Butyric Acid", "Methylbutyric Acid", "Isovaleric Acid", "Valeric Acid", "Hexanoic Acid"), times = 2)
pvalue_complete2 = data.frame(Treatment, Outcome, Mediator, pvalue_mat2)
pvalue_complete2
p_all = rbind(pvalue_complete, pvalue_complete2)
p_all[,4:7] = round(p_all[,4:7], digits = 2)
write.csv(p_all, "SPIRIT-output.csv")
