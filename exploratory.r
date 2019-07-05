#library(parallel)
library(splines)
library(reshape2)
library(tidyverse)
library(lmtest)
#cl <- makeCluster(3)

source('functions.r')
source('../generalCode/estimationFunctions.r')
source('~/Box Sync/rcode/binnedplot.r')

load('vetSurvData.RData')

dat <- filter(dat,agep>=25)

dat <- filter(dat,(dratx==1)|(deaf=='hearing'))

print(tab(~deaf+dratx+diss,dat))


## deafness by age
dat%>%
  mutate(deafBin=ifelse(deaf=='deaf',1,0))%>%
  group_by(agep)%>%
  summarize(perDeaf=svmean(deafBin,pwgtp),perSevere=svmean(drat==5,pwgtp))%>%
  gather(pp,percent,-agep)%>%
ggplot(aes(agep,percent,color=pp,group=pp))+
  geom_point()+
  geom_smooth()+
  xlab('Age')+
  scale_y_continuous(name='Percent Deaf',labels=scales::percent)

## enrollment proportions
enrModD <- glm(enrolled~ns(agep,5),family=binomial,data=filter(dat,deaf=='deaf'),weights=pwgtp)
enrModH <- glm(enrolled~ns(agep,5),family=binomial,data=filter(dat,deaf=='hearing'),weights=pwgtp)
pd <- data.frame(agep=seq(25,55))
pdD <- pd%>%mutate(pred=predict(enrModD,pd,type='response'),deaf='deaf')
pdH <- pd%>%mutate(pred=predict(enrModH,pd,type='response'),deaf='hearing')
pd <- rbind(pdD,pdH)
ggplot(pd,aes(agep,pred,group=deaf,color=deaf))+geom_line()+xlab('Proportion Enrolled')+ylab('Age')

enrProp <- dat%>%group_by(deaf,agep)%>%summarize(propEnr=svmean(enrolled,pwgtp))

ggplot(full_join(enrProp,pd),aes(agep,propEnr,color=deaf,group=deaf))+geom_point()+geom_line(aes(y=pred))+ylab('Proportion Enrolled')+xlab('Age')
ggsave('enrolledProportionByAge.jpg')

## transform data
ldat <- transformData(dat)


## now just deaf
modDeaf <- model(data=filter(ldat,deaf=='deaf'),surveySEs=FALSE)

## now both, just with interactions
modInt <- model(event~level*deaf,data=ldat,surveySEs=FALSE)
plotInt(modInt,ldat)
ggsave('noCovariates.pdf',width=10,height=5)
plotDiff(modInt,ldat)
ggsave('noCoveraitesDiff.pdf')


### now control for stuff
ldat$female <- ldat$sex=='Female'
ldat$raceEth <- relevel(factor(ldat$raceEth),ref='White')
ldat$ageCentered <- ldat$agep-35
#modFull <- model(event~level*deaf+female+ns(ageCentered,df=5)*deaf+raceEth+nativity,data=filter(ldat,level!='postBA'),surveySEs=FALSE)
modFull <- model(event~level*deaf+female+ns(ageCentered,df=5)*deaf+raceEth+nativity+recentVet,data=ldat,surveySEs=FALSE)
modFull2 <- model(event~level*deaf+female+ns(ageCentered,df=5)+raceEth+nativity,data=ldat,surveySEs=FALSE)
save(modFull,file='output/fullMod.RData')
plotInt(modFull,ldat)
ggsave('full.pdf',width=10,height=6)
plotInt(modFull,ldat,reverse=TRUE)
ggsave('fullReverse.pdf',width=10,height=6)

advFull <- shTab(modFull,ldat,'advanceProb')

fullCoef <- coeftest(modFull,vcov.=modFull$vcov)
hazardOddsRatios(modFull,ldat,age=35)
ggsave('oddRatiosFull35.pdf')

modFullRecent <-  model(event~level*deaf+female+ns(ageCentered,df=5)*deaf+raceEth+nativity,data=filter(ldat,recentVet),surveySEs=FALSE)
save(modFullRecent,file='output/fullModRecent.RData')
plotInt(modFullRecent,ldat)+ggtitle('Post-9/11 Vets')
ggsave('recent.pdf',width=10,height=6)

######### what about drat==5? (i.e. disability rating .70, 80, 90, or 100 percent)

ldatSev <- filter(ldat,deaf=='hearing'|drat==5)

### now control for stuff
modFullSev <- model(event~level*deaf+female+ns(ageCentered,df=5)*deaf+raceEth+nativity+recentVet,data=ldatSev,surveySEs=FALSE)
save(modFullSev,file='fullModSev.RData')
arm::binnedplot(predict(modFullSev),resid(modFullSev))
mf <- model.frame(modFullSev)
par(mfrow=c(3,3))
walk(mf%>%select(-event,-starts_with('ns'))%>%mutate(age=ldatSev$ageCentered),
  ~arm::binnedplot(.,resid(modFullSev)))

plotInt(modFullSev,ldat=ldatSev)+ggtitle('70%+ Disability Rating')
ggsave('fullSev.pdf',width=10,height=6)
plotInt(modFullSev,ldat=ldatSev,reverse=TRUE)+ggtitle('70%+ Disability Rating')
ggsave('fullReverseSev.pdf',width=10,height=6)
survFullSev <- survivalHazard(modFullSev,ldat=ldatSev)%>%select(level,deaf,hazard,overallAttainment=attainment)%>%
  mutate(advanceProb=1-hazard)%>%select(hazard,advanceProb,overallAttainment,everything())%>%
  melt()%>%dcast(level~variable+deaf)
fullCoefSev <- coeftest(modFullSev,vcov.=modFullSev$vcov)
hazardOddsRatios(modFullSev,ldatSev,age=35)+ggtitle('70%+ Disability Rating')
ggsave('oddRatiosFull35Sev.pdf')

modFullRecentSev <-  model(event~level*deaf+female+ns(ageCentered,df=5)*deaf+raceEth+nativity,data=filter(ldatSev,recentVet),surveySEs=FALSE)
save(modFullRecentSev,file='fullModRecentSev.RData')
plotInt(modFullRecentSev,ldatSev,reverse=TRUE)+ggtitle('Post-9/11 Vets, 70%+ Disability Rating')
ggsave('recentSev.pdf',width=10,height=6)


###################################################################################
#### ages 25-45
###################################################################################
dat3 <- filter(dat2,agep<46)


## transform data
ldat3 <- transformData(dat3)


## now just deaf
modDeaf3 <- model(data=filter(ldat3,deaf=='deaf'),surveySEs=FALSE)

## now both, just with interactions
modInt3 <- model(event~level*deaf,data=ldat3,surveySEs=FALSE)
plotInt(modInt3,ldat3)
ggsave('noCovariates3.pdf',width=10,height=5)
plotDiff(modInt3,ldat3)
ggsave('noCoveraitesDiff3.pdf')


### now control for stuff
ldat3$female <- ldat3$sex=='Female'
ldat3$raceEth <- relevel(factor(ldat3$raceEth),ref='White')
ldat3$ageCentered <- ldat3$agep-35
modFull3 <- model(event~level*deaf+female+ns(ageCentered,df=5)*deaf+raceEth+nativity,data=ldat3,surveySEs=TRUE)
save(modFull3,file='output/fullMod3.RData')
plotInt(modFull3,ldat=ldat3)+ggtitle('Ages 25-45')
ggsave('full3.pdf',width=10,height=6)
plotInt(modFull3,ldat=ldat3,reverse=TRUE)+ggtitle('Ages 25-45')
ggsave('fullReverse3.pdf',width=10,height=6)
survFull3 <- survivalHazard(modFull3,ldat=ldat3)%>%select(level,deaf,hazard,overallAttainment=attainment)%>%
  mutate(advanceProb=1-hazard)%>%select(hazard,advanceProb,overallAttainment,everything())%>%
  melt()%>%dcast(level~variable+deaf)
fullCoef3 <- coeftest(modFull3,vcov.=modFull3$vcov)
hazardOddsRatios(modFull3,ldat3,age=35)+ggtitle('Ages 25-45')
ggsave('oddRatiosFull353.pdf')

modFullRecent3 <-  model(event~level*deaf+female+ns(ageCentered,df=5)*deaf+raceEth+nativity,data=filter(ldat3,recentVet),surveySEs=TRUE)
save(modFullRecent3,file='output/fullModRecent3.RData')
plotInt(modFullRecent3,ldat3)+ggtitle('Post-9/11 Vets, 70%+ Disability Rating')
ggsave('recent3.pdf',width=10,height=6)



### what about when drat=1 (ie 0%)
ldat0 <- transformData(filter(dat,(deaf=='hearing')|(drat==1)))
ldat0$female <- ldat0$sex=='Female'
ldat0$raceEth <- relevel(factor(ldat0$raceEth),ref='White')
ldat0$ageCentered <- ldat0$agep-35
modFull0 <- model(event~level*deaf+female+ns(ageCentered,df=5)*deaf+raceEth+nativity,data=ldat0,surveySEs=FALSE)
save(modFull0,file='output/fullMod0.RData')
plotInt(modFull0,ldat=ldat0)+ggtitle('0% Disability Rating')
ggsave('full0.pdf',width=10,height=6)
plotInt(modFull0,ldat=ldat0,reverse=TRUE)+ggtitle('0% Disability Rating')
ggsave('fullReverse0.pdf',width=10,height=6)

hazardOddsRatios(modFull0,ldat0,age=35)+ggtitle('0% Disability Rating')
ggsave('oddRatiosFull350.pdf')
