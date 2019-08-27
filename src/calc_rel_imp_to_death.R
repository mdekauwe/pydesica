#!/usr/bin/Rscript

library(visreg)
library(ggplot2)
library(ppcor)
library(relaimpo)
library(tidyverse)


setwd("/Users/mdekauwe/src/python/pydesica/outputs")

pft = "rf"
fname <- sprintf("%s_trait_sens_OAT.csv", pft)
df <- read.csv(fname)
fit <- lm(day_of_death ~ gmin + lai + p50 + depth + Cl + Cs, data=df)
rel_imp <- calc.relimp(fit, type=c("lmg"), rela=TRUE)
#plot(rel_imp)
rf_imp <- rel_imp$lmg

pft = "wsf"
fname <- sprintf("%s_trait_sens_OAT.csv", pft)
df <- read.csv(fname)
fit <- lm(day_of_death ~ gmin + lai + p50 + depth + Cl + Cs, data=df)
rel_imp <- calc.relimp(fit, type=c("lmg","last", "first"), rela=TRUE)
wsf_imp <- rel_imp$lmg

pft = "dsf"
fname <- sprintf("%s_trait_sens_OAT.csv", pft)
df <- read.csv(fname)
fit <- lm(day_of_death ~ gmin + lai + p50 + depth + Cl + Cs, data=df)
rel_imp <- calc.relimp(fit, type=c("lmg","last", "first"), rela=TRUE)
dsf_imp <- rel_imp$lmg

pft = "grw"
fname <- sprintf("%s_trait_sens_OAT.csv", pft)
df <- read.csv(fname)
fit <- lm(day_of_death ~ gmin + lai + p50 + depth + Cl + Cs, data=df)
rel_imp <- calc.relimp(fit, type=c("lmg","last", "first"), rela=TRUE)
grw_imp <- rel_imp$lmg

pft = "saw"
fname <- sprintf("%s_trait_sens_OAT.csv", pft)
df <- read.csv(fname)
fit <- lm(day_of_death ~ gmin + lai + p50 + depth + Cl + Cs, data=df)
rel_imp <- calc.relimp(fit, type=c("lmg","last", "first"), rela=TRUE)
saw_imp <- rel_imp$lmg

# Clean up and save the output
all <- rbind(rf_imp, wsf_imp, dsf_imp, grw_imp, saw_imp)
all <- cbind(c("RF", "WSF", "DSF", "GRW", "SAW"), all)
colnames(all) <- c("pft","gmin","lai","p50","depth","Cl","Cs")
write.csv(all,"rel_imp.csv", row.names=FALSE)

