#!/usr/bin/Rscript

library(visreg)
library(ggplot2)
library(ppcor)
library(relaimpo)
library(tidyverse)


setwd("/Users/mdekauwe/src/python/pydesica/outputs")


pfts <- c("rf", "wsf", "dsf", "grw", "saw")
vars <- c("gmin","lai", "p50", "depth", "Cl", "Cs", "b", "psi_e")
all <- matrix(nrow=length(pfts), ncol=length(vars))
for (i in 1:length(pfts)) {
  #fname <- sprintf("%s_trait_sens_OAT.csv", pfts[i])
  fname <- sprintf("%s_trait_sensitivity_all.csv", pfts[i])
  df <- read.csv(fname)
  df <- df[df$day_of_death> 0, ]
  fit <- lm(day_of_death ~ gmin + lai + p50 + depth + Cl + Cs + b + psi_e, data=df)
  rel_imp <- calc.relimp(fit, type=c("lmg"), rela=TRUE)
  all[i,] <-  rel_imp$lmg

}

all <- cbind.data.frame(pfts, all)
colnames(all) <- c("pft","gmin","lai","p50","depth","Cl","Cs","b","psi_e")
write.csv(all,"rel_imp.csv", row.names=FALSE)

