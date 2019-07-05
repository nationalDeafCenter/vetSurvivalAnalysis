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

full <- function(dat,survSE=TRUE,centeringAge=35){
  crossTabs <- addmargins(xtabs(
    if(any(dat$diss=='disabled')) ~deaf+drat+recentVet+diss else ~deaf+drat+recentVet,
    dat,addNA=TRUE
  ))

  rates <- dat%>%
    mutate(Nsamp=n(),Npop=sum(pwgtp))%>%
    bind_rows(filter(.,drat==5)%>%mutate(deaf='severe'))%>%
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

  form <- event~level*deaf+female+ns(ageCentered,df=5)*deaf+raceEth+nativity
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

system('mv fittedModels/results* fittedModels/old')

save(results,file=paste0('fittedModels/results',Sys.Date(),'.RData'))

makeSheet <- function(x){
  full_join(
    x$raw[-1,]%>%rename(deaf_raw=deaf,hearing_raw=hearing,severe_raw=severe),
    x$tables%>%mutate_if(is.numeric,~.*100))%>%
    mutate(`hazard odds ratio`=hazard_deaf*(100-hazard_hearing)/(hazard_hearing*(100-hazard_deaf)))%>%
    add_case(level='no other ACS disabilities; hearing=w/ and w/o rating; deaf= w/rating')%>%
    add_case(level=c(toupper('model-based results for:'),apply(attr(x$tables,"forWhom"),1,paste,collapse='=')))

}


tables <- map(results,makeSheet)
tables$recent <- add_case(tables$recent, level='SAMPLE: Post-9/11 vets')
tables$severe <- add_case(tables$severe,level=c('DEFINITION: deaf= only 70+ rating'))

system('mv output/* output/old')

openxlsx::write.xlsx(tables,paste0('output/resultsTables',Sys.Date(),'.xlsx'))

plotInt(results$total$mod,reverse=TRUE,recent=TRUE)
ggsave(paste0('output/full',Sys.Date(),'.pdf'),width=10,height=6)

plotInt(results$recent$mod,reverse=TRUE)+ggtitle('Post 9/11 vets')
ggsave(paste0('output/recent',Sys.Date(),'.pdf'),width=10,height=6)

plotInt(results$severe$mod,reverse=TRUE,recent=TRUE)+ggtitle('70+ Disability Rating')
ggsave(paste0('output/severe',Sys.Date(),'.pdf'),width=10,height=6)


hazardOddsRatios(results$total$mod,age=35)
ggsave(paste0('output/hazardRatiosFull',Sys.Date(),'.pdf'))

