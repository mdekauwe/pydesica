#!/usr/bin/Rscript

library(visreg)
library(ggplot2)

setwd("/Users/mdekauwe/src/python/pydesica/outputs")

pft = "dsf"
fname <- sprintf("%s_trait_sensitivity.csv", pft)
df <- read.csv(fname)
#head(df)

fit <- lm(plc ~ gmin + AL + p50 + Cl + Cs + Tmax + Dmax, data=df)
visreg(fit)
summary(fit)

#par(mfrow=c(2,4))
#visreg(fit, "Tmax")
#visreg(fit, "Dmax")
#visreg(fit, "gmin")
#visreg(fit, "AL")
#visreg(fit, "p50")
#visreg(fit, "Cl")
#visreg(fit, "Cs")

