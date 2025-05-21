#Scritpt - JEcol 2025 - Reassembly dynamics of tropical secondary succession: evidence from Atlantic Forests
#Toledo RM; Pivello VR; Vangansbeke P; Jakovac C; Verheyen K; Verdade LM; Lima RAF.

library(vegan)
library(reshape2)
library(ggplot2)
library(nlme)
library(sciplot)
library(dplyr)
library(indicspecies)
library(mvabund)
library(tidyr)
library(ggplot2)
library(MuMIn)
library(rstatix)

#A1- First, composition similarity between and within AFSubtype and succession. Then,  additional effects of methods and GEO-context  

dtaS1 <- read.csv("T1_CommunityMatrix_SIM.csv", header = T)
M <- dtaS1[, c(10:2726)]
dst.env <-  dtaS1[, c(1:9)]
dstS1 <- vegdist(M, method="bray")
braymat <- as.matrix(dstS1)
write.csv(braymat,"braymat.csv")

#A1 data for figure2

braylong <- melt(braymat)
write.csv(braylong,"braylong.csv")

#A1a - AFSubtype and succession
P.div1 <- adonis2(M ~ succession*af.SubType, data = dst.env, permutations = 999, by = "terms", method="bray")
P.div1
#A1b - Adding survey effort (ha) and inclusion DBH (cm)
P.div2 <- adonis2(M ~ succession*af.SubType*effort*c.o.DBHcm, data = dst.env, permutations = 999, by = "terms", method="bray")
P.div2
#A1b - Adding soil and climate
P.div3 <- adonis2(M ~ succession*af.SubType*climate*soil, data = dst.env, permutations = 999, by = "terms", method="bray")
P.div3

#####
#### Indicator species for Early secondary forests and Old Growth

options(max.print=90000)

#ARAUCARIA
dta_ARA <- read.csv("T2_IndicSPP_ARA.csv", header = T, sep = ",")
COMAT <- dta_ARA[,c(4:2720)]
groups_ARA = c(rep(1, 30), rep(2, 28))
indval_ARA = multipatt(COMAT, groups_ARA, control = how(nperm=999))
summary(indval_ARA)
sink("INDICSPP_ARA_OUT.txt")
summary(indval_ARA, indvalcomp=TRUE, alpha=1)
sink()

#RAINFOREST
dta_RAI <- read.csv("T3_IndicSPP_RAI.csv", header = T, sep = ",")
COMAT <- dta_RAI[,c(4:2752)]
groups_RAI = c(rep(1, 118), rep(2, 61))
indval_RAI = multipatt(COMAT, groups_RAI, control = how(nperm=999))
sink("INDICSPP_RAI_OUT.txt")
summary(indval_RAI, indvalcomp=TRUE, alpha=1)
sink()

#SEASONAL
dta_SEA <- read.csv("T4_IndicSPP_SEA.csv", header = T, sep = ",")
COMAT <- dta_SEA[,c(4:2752)]
groups_SEA = c(rep(1, 168), rep(2, 78))
indval_SEA = multipatt(COMAT, groups_SEA, control = how(nperm=999))
sink("INDICSPP_SEA_OUT.txt")
summary(indval_SEA, indvalcomp=TRUE, alpha=1)
sink()

#####
#### MVABUND - species responses to succession accounting for AF-subtype and survey methods 

dta <- read.csv("T5_CommunityMatrix_mvabund.csv", header = T)
#summarize predictors at site level##
env <- dta %>%
  select(Survey:SuccessionIndex) %>%
  group_by(Survey) %>%
  summarize(DBH = (DBH), Stage =(Stage), SuccessionIndex = (SuccessionIndex), Subtype = (Subtype), effort = (effort))

# Mvabund prep objects #

treedata <- select(dta, Survey:Survey, A.bra:Z.lat)
treedat_mva<-mvabund(treedata[, -1])

# Relative abundance as a function of the successional index (EarlyIndicators:OldGrowthIndicators) ####  

glm_treea <- manyglm(treedat_mva~ Subtype + DBH + effort + SuccessionIndex+SuccessionIndex:Subtype, data = env,  family = "binomial")
newdata <- data.frame(expand.grid(SuccessionIndex = seq(-0.8, 0.8, 1.6), DBH = 5, effort  = 1, Subtype = factor(c("rainforest", "seasonal", "araucaria"))))
preds <- predict(glm_treea, newdata = newdata, type = "response")
preds_env <- data.frame(newdata, preds)
Apspp <- preds_env %>%
  gather(Species, Abundance, A.bra:Z.lat)
write.csv(Apspp,"Preds_BINOMIAL.csv")
sink("MVABUND_output.txt")
anova(glm_treea, nBoot = 999, p.uni = "unadjusted", show.time ="total")
sink()

### tests figue 3 - Differences of species groups abundance along the successional gradient

dtRAG <- read.csv("T6_GroupsRelativeAbundance_and_stages.csv", header = T)
fig3 <- aov(dtRAG$RelativAbund ~ dtRAG$Stage, data = dtRAG)
summary(fig3)
TukeyHSD(fig3)

#####
### Functional patterns 

## SPECIES LEVEL - Successional groups and traits (Figure 4)

tdta <- read.csv("T7_ARAUCARIA_traits_standartized.csv", header = T)
tdtr <- read.csv("T8_RAINFOREST_traits_standartized.csv", header = T)
tdts <- read.csv("T9_SEAZONAL_traits_standartized.csv", header = T)


#log transformed Specific Leaf Area SLA

Lma1 <- aov(tdta$L10SLA ~ tdta$ARA, data = tdta)
summary(Lma1)

Lmc1 <- aov(tdtr$L10SLA ~ tdtr$RAI, data = tdtr)
summary(Lmc1)

Lmc1 <- aov(tdts$L10SLA ~ tdts$SEA, data = tdts)
summary(Lmc1)
TukeyHSD(Lmc1)


#log transformed Seed Mass

Lma4 <- lm((L10SM)~(ARA), data = tdta)
summary(Lma4)

Lmb4 <- lm((L10SM)~(RAI), data = tdtr)
summary(Lmb4)

Lmc4 <- lm((L10SM)~(SEA), data = tdts)
summary(Lmc4)

kruskal.test(tdta$L10SM ~ tdta$ARA, data = tdta)
dunn_test(L10SM ~ ARA, data = tdta, p.adjust.method = "bonferroni")
e2smar <- kruskal.test(tdta$L10SM ~ tdta$ARA, data = tdta)
H <- e2smar$statistic
k <- length(unique(tdta$ARA))
n <- length(tdta$L10SM)
epsilon_squared <- (H - k + 1) / (n - k)
epsilon_squared

kruskal.test(tdtr$L10SM ~ tdtr$RAI, data = tdtr)
dunn_test(L10SM ~ RAI, data = tdtr, p.adjust.method = "bonferroni")
e2smra <- kruskal.test(tdtr$L10SM ~ tdtr$RAI, data = tdtr)
H <- e2mara$statistic
k <- length(unique(tdtr$RAI))
n <- length(tdtr$L10SM)
epsilon_squared <- (H - k + 1) / (n - k)
epsilon_squared

kruskal.test(tdts$L10SM ~ tdts$SEA, data = tdts)
dunn_test(L10SM ~ SEA, data = tdts, p.adjust.method = "bonferroni")
e2smse <- kruskal.test(tdts$L10SM ~ tdts$SEA, data = tdts)
H <- e2smse$statistic
k <- length(unique(tdts$SEA))
n <- length(tdts$L10SM)
epsilon_squared <- (H - k + 1) / (n - k)
epsilon_squared


#Wood specific gravity

Lma2 <- aov((tdta$WSG)~(tdta$ARA), data = tdta)
summary(Lma2)

Lmb2 <- aov((tdtr$WSG)~(tdtr$RAI), data = tdtr)
summary(Lmb2)
TukeyHSD(Lmb2)

Lmc2 <- aov((tdts$WSG)~(tdts$SEA), data = tdts)
summary(Lmc2)

#Maximum Height
Lma2 <- aov((tdta$WSG)~(tdta$ARA), data = tdta)
summary(Lma2)

Lmb2 <- aov((tdtr$WSG)~(tdtr$RAI), data = tdtr)
summary(Lmb2)
TukeyHSD(Lmb2)

Lmc2 <- aov((tdts$WSG)~(tdts$SEA), data = tdts)
summary(Lmc2)

kruskal.test(tdta$MaxH ~ tdta$ARA, data = tdta)
e2maar <- kruskal.test(tdta$MaxH ~ tdta$ARA, data = tdta)
H <- e2maar$statistic
k <- length(unique(tdta$ARA))
n <- length(tdta$MaxH)
epsilon_squared <- (H - k + 1) / (n - k)
epsilon_squared

dunn_test(MaxH ~ ARA, data = tdta, p.adjust.method = "bonferroni")


kruskal.test(tdtr$MaxH ~ tdtr$RAI, data = tdtr)
dunn_test(MaxH ~ RAI, data = tdtr, p.adjust.method = "bonferroni")
e2mara <- kruskal.test(tdtr$MaxH ~ tdtr$RAI, data = tdtr)
H <- e2mara$statistic
k <- length(unique(tdtr$RAI))
n <- length(tdtr$MaxH)
epsilon_squared <- (H - k + 1) / (n - k)
epsilon_squared

kruskal.test(tdts$MaxH ~ tdts$SEA, data = tdts)
dunn_test(MaxH ~ SEA, data = tdts, p.adjust.method = "bonferroni")
e2mase <- kruskal.test(tdts$MaxH ~ tdts$SEA, data = tdts)
H <- e2mase$statistic
k <- length(unique(tdts$SEA))
n <- length(tdts$MaxH)
epsilon_squared <- (H - k + 1) / (n - k)
epsilon_squared

#####
#SITE LEVEL - Linear mixed models of traits CWM as a function of Successional Stage
# see FIGURE 5
# considering surveys minimum DBH as a random factor

tdtb <- read.csv("T10_stdCWMsites.csv", header = T)

#sla
LME1<- lme(wgsSLA~Stage, random =~ 1 | catDBH, data = tdtb, na.action = na.omit)
summary(LME1)
r.squaredGLMM(LME1, null, envir = parent.frame(), pj2014 = FALSE)

#wood 
LME2<- lme(wgsWD~Stage, random =~ 1 | catDBH, data = tdtb, na.action = na.omit)
summary(LME2)
r.squaredGLMM(LME2, null, envir = parent.frame(), pj2014 = FALSE)

#seed 
LME3<- lme(wgsSEED~Stage, random =~ 1 | catDBH, data = tdtb, na.action = na.omit)
summary(LME3)
r.squaredGLMM(LME3, null, envir = parent.frame(), pj2014 = FALSE)

#height
LME4<- lme(wgsHEIGHT~Stage, random =~ 1 | catDBH, data = tdtb, na.action = na.omit)
summary(LME4)
r.squaredGLMM(LME4, null, envir = parent.frame(), pj2014 = FALSE)

# also considering survey effort as fixed factor surveys minimum DBH as a random factor

LME1b<- lme(wgsSLA~effort+Stage, random =~ 1 | catDBH, data = tdtb, na.action = na.omit)
summary(LME1b)
r.squaredGLMM(LME1b, null, envir = parent.frame(), pj2014 = FALSE)

LME2b<- lme(wgsWD~effort+Stage, random =~ 1 | catDBH, data = tdtb, na.action = na.omit)
summary(LME2b)
r.squaredGLMM(LME2b, null, envir = parent.frame(), pj2014 = FALSE)

LME3b<- lme(wgsSEED~effort+Stage, random =~ 1 | catDBH, data = tdtb, na.action = na.omit)
summary(LME3b)
r.squaredGLMM(LME3b, null, envir = parent.frame(), pj2014 = FALSE)

LME4b<- lme(wgsHEIGHT~effort+Stage, random =~ 1 | catDBH, data = tdtb, na.action = na.omit)
summary(LME4b)
r.squaredGLMM(LME4b, null, envir = parent.frame(), pj2014 = FALSE)