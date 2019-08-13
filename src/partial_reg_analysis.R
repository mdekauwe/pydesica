#!/usr/bin/Rscript

library(visreg)
library(ggplot2)
library(ppcor)
library(relaimpo)
library(tidyverse)


setwd("/Users/mdekauwe/src/python/pydesica/outputs")

pft = "saw"
fname <- sprintf("%s_trait_sensitivity.csv", pft)
df <- read.csv(fname)
#head(df)

fit <- lm(plc ~ gmin + p50 + Cl + Cs + depth, data=df)

#visreg(fit)
par(mfrow=c(2,3))
visreg(fit, "gmin")
visreg(fit, "p50")
visreg(fit, "Cl")
visreg(fit, "Cs")
visreg(fit, "depth")

summary(fit)


rel_imp <- calc.relimp(fit, type=c("lmg","last", "first"), rela=TRUE)
plot(rel_imp)

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
