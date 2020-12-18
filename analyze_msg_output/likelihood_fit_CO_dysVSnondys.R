####################################################################################
#### Fitting likelihood models to determine effects on crossover number
####################################################################################


#### Load libraries

library(lme4)
library(lsmeans)

#### Import data

brks <- read.csv("Hemmer_etal_2020_MobileDNA/data/CO_sum_per_parent_FileS8.csv", header=TRUE)

brks$par <- as.character(brks$par)
brks$batch <- as.character(brks$batch)

#### analyze crossover number per parent

## student t-test first between dysgenic and non-dysgenic progeny

t.test(brks$CO.sum[brks$dys=="dys"],brks$CO.sum[brks$dys=="nondys"]) # p = 0.4278, non-significant

## fitting a comprehesive likelihood model, poisson distribution, for batch, dysgenic / non-dysgenic, fecundity
# parent, fecundity are treated as random variables

fit.full <- glmer(CO.sum~batch*dys.nondys*(1|dys.nondys:parent)*(1|dys.nondys:parent:all.progeny.count),data=brks,family=poisson(link=log))

summary(fit.full)

## fecundity doesn't seem to make a difference, neither does parent as random effects, p > 0.05
# can remove random effects from likelihood model

## fit model with just effects of dysgenesis

fit.dys <- glm(CO.sum~dys.nondys,data=brks,family=poisson(link=log))

summary(fit.dys) # p = 0.523, no difference between dysgenic and non-dysgenic


#### compare model 

anova(fit.full,fit.dys) ## p = 0.5155, keep simplier model


#### moral of the story, random effects have no effect, full model is not preferred over (almost) bare minimum

#### almost a difference between dys and nondys, believe a lot of that is driven by parent 7

#### lets see if there is a difference between dysgenic flies and set, nondysgenic flies and set, no fecundity 

## select for dysgenic only

brks.dys <- brks[brks$dys=="dys",]
fit.dys.set <- glmer(CO.sum~batch*(1|parent),data=brks.dys,family=poisson(link=log))

summary(fit.dys.set) ## effect of batch is insignificant p = 0.291
car::Anova(fit.dys.set, type=3)

## parent and fecundity do not make much of a difference

## non-dysgenic now

brks.non <- brks[brks$dys=="nondys",]
fit.non.set <- glmer(sum~batch*(1|parent),data=brks.non,family=poisson(link=log))

summary(fit.non.set) # batch no effect, p = 0.2383
car::Anova(fit.non.set, type=3)

## parent really has no effect, no difference between batch 


#### what about low fecund dysgenic
brks.low <- brks[brks$dys=="dys" & brks$fecund=="low",]
fit.dys.low <- glm(sum~par,data=brks.low,family=poisson(link=log))


summary(fit.dys.low)
car::Anova(fit.dys.low, type=3)

## nothing different between low 

#### high fecund dysgenic

brks.high <- brks[brks$dys=="dys" & brks$fecund=="high",]
fit.dys.high <- glm(sum~par,data=brks.high,family=poisson(link=log))

summary(fit.dys.high) # par701 p = 0.0234, others p > 0.05
car::Anova(fit.dys.high, type=3) ## p = 0.01707

## parent has a signficant effect, specifically from parent 701 

means.fit <- lsmeans(fit.dys.high,pairwise~par)
contrast(means.fit)

#### parent 701 is different from the rest, also effect of 4029 in means.contrast

## lets see what is going on with 701
brks.701.ab <- brks[brks$par!="701",]
fit.701.ab <- glmer(sum~batch*dys.nondys*(1|dys.nondys:parent)*(1|dys.nondys:parent:all.progeny.count),data=brks.701.ab,family=poisson(link=log))

summary(fit.701.ab)
car::Anova(fit.701.ab, type=3)

#### p values increased, crossovers in progeny of 701 driving signal of near signficance

fit.701.ab.dys <- glmer(sum~dys.nondys*(1|dys.nondys:par)*(1|dys.nondys:par:all.progeny.count),data=brks.701.ab,family=poisson(link=log))
summary(fit.701.ab.dys)
car::Anova(fit.701.ab.dys, type=3)

################################
#### no difference in crossover number between dysgenic and non-dysgenic progeny
#### progeny of 701 have slightly higher crossover numbers among dysgenic flies, but does little to the effect of dysgenesis
################################

