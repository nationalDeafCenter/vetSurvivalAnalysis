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

dat <- dat%>%mutate(foreignBorn=nativity==2)

print(tab(~deaf+dratx+diss,dat))

full <- function(dat,survSE=TRUE,centeringAge=35){
  crossTabs <- addmargins(xtabs(
    if(any(dat$diss=='disabled')) ~deaf+drat+recentVet+diss else ~deaf+drat+recentVet,
    dat,addNA=TRUE
  ))

  rates <- dat%>%
    mutate(Nsamp=n(),Npop=sum(pwgtp))%>%
    bind_rows(filter(.,drat==5)%>%mutate(deaf=ifelse(deaf=='deaf','deaf: severe','hearing:severe')))%>%
    group_by(deaf)%>%
    mutate(n=n(),propSamp=n/Nsamp*100,propPop=sum(pwgtp)/Npop*100)%>%
    group_map(
      ~bind_rows(
        tibble(
          level=levels(dat$postSec),#levels=c(levels(dat$postSec)),
          percent=sapply(level,function(x) svmean(.x$postSec>=x,.x$pwgtp)*100)
        ),
        tibble(level=c('n','propSamp','propPop'),percent=unlist(.x[1,level]))
      )%>%
        mutate(level=factor(level,levels=level))
    )%>%
    spread(deaf,percent)

  ## transform data
  ldat <- transformData(dat)

  ### now control for stuff
  ldat$female <- ldat$sex=='Female'
  ldat$raceEth <- relevel(factor(ldat$raceEth),ref='White')
  ldat$ageCentered <- ldat$agep-centeringAge

  form <- event~level*deaf+female+ns(ageCentered,df=5)*deaf+raceEth+foreignBorn
  if(between(mean(dat$recentVet),0.01,0.99)) form <- update(form,.~.+recentVet)

  mod <- model(form,data=ldat,surveySEs=survSE)

  list(
    raw=rates,
    sample=crossTabs,
    mod=mod,
    tables=shTab(mod,recent=TRUE,centeringAge=centeringAge)
  )
}

results <-
  list(
    total=full(dat),
    recent=full(filter(dat,recentVet),FALSE),
    severe=full(filter(dat,(drat==5)|(deaf=='hearing')),FALSE)
  )

system('mv fittedModels/*.RData fittedModels/old')

save(results,file=paste0('fittedModels/results',Sys.Date(),'.RData'))

makeSheet <- function(x){
  full_join(
    x$raw[-1,]%>%rename(deaf_raw=deaf,hearing_raw=hearing,severe_deaf_raw="deaf: severe",severe_hearing_raw='hearing:severe'),
    x$tables%>%mutate_if(is.numeric,~.*100))%>%
    mutate(`hazard odds ratio`=hazard_deaf*(100-hazard_hearing)/(hazard_hearing*(100-hazard_deaf)))%>%
    add_case(level='no other ACS disabilities; hearing=w/ and w/o rating; deaf= w/rating')%>%
    add_case(level=c(toupper('model-based results for:'),apply(attr(x$tables,"forWhom"),1,paste,collapse='=')))

}


tables <- map(results,makeSheet)
tables$recent <- add_case(tables$recent, level='SAMPLE: Post-9/11 vets')
tables$severe <- add_case(tables$severe,level=c('DEFINITION: deaf= only 70+ rating'))

system('mv output/*.xlsx output/old')

openxlsx::write.xlsx(tables,paste0('output/resultsTables',Sys.Date(),'.xlsx'))

system('mv output/*.pdf output/old')

plotInt(results$total$mod,reverse=TRUE,recent=TRUE)
ggsave(paste0('output/full',Sys.Date(),'.pdf'),width=10,height=6)

plotInt(results$recent$mod,reverse=TRUE)+ggtitle('Post 9/11 vets')
ggsave(paste0('output/recent',Sys.Date(),'.pdf'),width=10,height=6)

plotInt(results$severe$mod,reverse=TRUE,recent=TRUE)+ggtitle('70+ Disability Rating')
ggsave(paste0('output/severe',Sys.Date(),'.pdf'),width=10,height=6)


hazardOddsRatios(results$total$mod,age=35)
ggsave(paste0('output/hazardRatiosFull',Sys.Date(),'.pdf'))

map(quantile(dat$agep), ~hazardOddsRatios(results$total$mod,age=.,returnData=TRUE))%>%#mutate(Age=.))%>%
  bind_rows()%>%
  mutate(Age=as.factor(Age))%>%

  hazAges%>%ggplot(aes(level,hazardRatio,group=gr))+
  geom_point()+geom_line()+
  geom_errorbar(aes(ymin=ciL,ymax=ciH),width=0.1)+
  scale_y_continuous(breaks=function(y) seq(0,ceiling(max(y)),.25),limits=c(0,NA))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab(NULL)+ylab('Hazard Odds Ratio deaf/hearing')+
  geom_hline(yintercept=1)+
  facet_wrap(~Age,labeller=function(variable,value) paste('Age',value))

ggsave('hazardRatiosByAge.jpg')


### what if we remove AA?

## option 1: group AA with previous level:

resultsNoAA <- list(
  `AA as some college`=full(
    within(
      dat,
      postSec[postSec=='Associates degree'] <- '1 or more years of college credit, no degree'
    ),
    survSE=FALSE
  ),
  `undergrad degree`=full(
    mutate(dat,postSec=fct_collapse(postSec,`Undergrad Degree`=c("Associates degree","Bachelors degree"))),
    survSE=FALSE
  )
)
save(resultsNoAA,file=paste0('fittedModels/resultsNoAA',Sys.Date(),'.RData'))

openxlsx::write.xlsx(map(resultsNoAA,makeSheet),paste0('output/resultsNoAA',Sys.Date(),'.xlsx'))
