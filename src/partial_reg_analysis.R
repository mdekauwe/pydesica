#!/usr/bin/Rscript

library(visreg)
library(ggplot2)
library(ppcor)
library(relaimpo)
library(tidyverse)


setwd("/Users/mdekauwe/src/python/pydesica/outputs")

#pft = "wsf"
#fname <- sprintf("%s_trait_sensitivity.csv", pft)
#df <- read.csv(fname)

pft = "grw"
fname <- sprintf("%s_trait_sens_OAT.csv", pft)
df <- read.csv(fname)

df <- df[df$day_of_death> 0, ]

#gmin <- df %>% distinct(gmin, .keep_all=TRUE)
#lai <- df %>% distinct(lai, .keep_all=TRUE)
#p50 <- df %>% distinct(p50, .keep_all=TRUE)
#depth <- df %>% distinct(depth, .keep_all=TRUE)
#Cl <- df %>% distinct(Cl, .keep_all=TRUE)
#Cs <- df %>% distinct(Cs, .keep_all=TRUE)
#dfx <- rbind(gmin, lai, p50, depth, Cl, Cs)


#df <- df[!duplicated(df), ]
df <- unique(df)
fit <- lm(day_of_death ~ gmin + lai + p50 + depth + Cl + Cs, data=df)
visreg(fit)



#x <- pcor(df[,c("gmin","lai","p50","depth","Cl","Cs")])
#sens_2 <- sort(abs(x$estimate[1,-1]),T)

#visreg(fit)
par(mfrow=c(3,3))
visreg(fit, "gmin")
visreg(fit, "lai")
visreg(fit, "p50")
visreg(fit, "Cl")
visreg(fit, "Cs")
visreg(fit, "depth")

summary(fit)


rel_imp <- calc.relimp(fit, type=c("lmg"), rela=TRUE)
plot(rel_imp)

pft = "rf"
fname <- sprintf("%s_trait_sens_OAT.csv", pft)
df <- read.csv(fname)
fit <- lm(day_of_death ~ gmin + lai + p50 + depth + Cl + Cs, data=df)
rel_imp <- calc.relimp(fit, type=c("lmg","last", "first"), rela=TRUE)
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




boot <- boot.relimp(fit, b=1000, type=c("lmg","last", "first"),
                    rank=TRUE, diff=TRUE, rela=TRUE)
#booteval.relimp(boot)
plot(booteval.relimp(boot, sort=TRUE))


df_dead <- df %>% filter(day_of_death > 0)
dev.off()
ggplot(df_dead, aes(x=day_of_death)) + geom_histogram()

fit <- lm(day_of_death ~ gmin + p50 + Cl + Cs + depth, data=df_dead)

#visreg(fit)
par(mfrow=c(2,3))
visreg(fit, "gmin")
visreg(fit, "p50")
visreg(fit, "Cl")
visreg(fit, "Cs")
visreg(fit, "depth")

rel_imp <- calc.relimp(fit, type=c("lmg","last", "first"), rela=TRUE)
plot(rel_imp)

boot <- boot.relimp(fit, b=1000, type=c("lmg","last", "first"),
                    rank=TRUE, diff=TRUE, rela=TRUE)
#booteval.relimp(boot)
plot(booteval.relimp(boot, sort=TRUE))

