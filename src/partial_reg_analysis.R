#!/usr/bin/Rscript

library(visreg)
library(ggplot2)
library(ppcor)
library(relaimpo) 

setwd("/Users/mdekauwe/src/python/pydesica/outputs")

pft = "dsf"
fname <- sprintf("%s_trait_sensitivity.csv", pft)
df <- read.csv(fname)
#head(df)

#fit <- lm(plc ~ gmin + AL + p50 + Cl + Cs + Tmax + Dmax, data=df)
fit <- lm(plc ~ gmin + AL + p50 + Cl + Cs, data=df)

#visreg(fit)
par(mfrow=c(2,4))
visreg(fit, "Tmax")
visreg(fit, "Dmax")
visreg(fit, "gmin")
visreg(fit, "AL")
visreg(fit, "p50")
visreg(fit, "Cl")
visreg(fit, "Cs")

summary(fit)

 
rel_imp <- calc.relimp(fit, type=c("lmg","last", "first"), rela=TRUE)
plot(rel_imp)

boot <- boot.relimp(fit, b=1000, type=c("lmg","last", "first"),
                    rank=TRUE, diff=TRUE, rela=TRUE)
#booteval.relimp(boot) 
plot(booteval.relimp(boot, sort=TRUE))


