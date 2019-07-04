
## for suvival model: (nothing in stone here)
## exclude disabled
## include only dratx=1
## use drat as cont covariate (?)
## compare deaf nondisabled dratx=1 to hearing nondisabled dratx=0 vets
## period of service: pre/post 911
library(tidyverse)
library(openxlsx)
source('../generalCode/estimationFunctions.r')
states <- read.csv('../../../data/states.csv')

#1) a simple breakdown of current enrollment, and completion data, across type of institution (4 year colleges, community colleges, etc) using all the 'type of institution' data we have, so that would give us some nice descriptives and allow us to make a final decision on how we want to categorize 'community colleges and 2-year institutions'




varNames <- c('ST','AGEP','DDRS','DEAR','DEYE','DOUT','DPHY','DRATX','DREM','DRAT','FDEARP','ESR','SCHL','SCH','SCHG','RAC1P','HISP','SEX','PERNP','PINCP','SSIP','WKHP','WKW','ADJINC','PWGTP','RELP','DRAT','MIL','NATIVITY','POBP','QTRBIR','VPS',
  'WAOB',
              paste0('MLP',c('A','B','CD','E','FG',LETTERS[8:11])),
              paste0('PWGTP',1:80))

ctypes <- rep('i',length(varNames))
names(ctypes) <- varNames
ctypes$.default <- '_'
colTypes <- do.call('cols',as.list(ctypes))

dat <- read_csv('../../../data/acs5yr2017/psam_pusa.csv',col_types=colTypes)
dim(dat)

setdiff(names(dat),varNames)

for(ll in c('b','c','d'))
  dat <- bind_rows(dat,
    read_csv(paste0('../../../data/acs5yr2017/psam_pus',ll,'.csv'),col_types=colTypes))
gc()

names(dat) <- tolower(names(dat))

dat$state <- states$abb[match(dat$st,states$x)]
gc()
### let's figure out what the definition of "veteran" is
### https://www.census.gov/acs/www/about/why-we-ask-each-question/veterans/ says 8%
### https://factfinder.census.gov/faces/tableservices/jsf/pages/productview.xhtml?src=bkmk "civilian veterans" 8.0%
## factorProps('mil',dat,cum=FALSE)
##     % 1            1 SE             % 2            2 SE             % 3
## 0.403647        0.002822        7.963292        0.009436        1.291155
##     3 SE             % 4            4 SE               n
## 0.004261       90.341905        0.011969 12364760.000000
## well, it looks like mil==2 does it.
##  from data dictionary (MIL):   2    .On active duty in the past, but not now


edlevs <- c(
    '<Grade 10',
    'Grade 10',
    'Grade 11',
    '12th grade - no diploma',
    'Regular high school diploma',
    'GED or alternative credential',
    'Some college, but less than 1 year',
    '1 or more years of college credit, no degree',
    'Associates degree',
    'Bachelors degree',
    'Masters degree',
    'Professional degree beyond a bachelors degree',
    'Doctorate degree')

dat$attain <- ifelse(dat$schl<13,1,dat$schl-11)
dat$attain <- factor(edlevs[dat$attain],levels=edlevs,ordered=TRUE)

gc()

dat <- dat%>%filter(agep>17,agep<55,mil==2)%>% ## mil==2 => only vets in dataset
  mutate(
    currentMil=factor(ifelse(esr%in%4:5,'Currently serving','Not currently serving')),
    recentVet=!is.na(mlpa)&(mlpa==1), ## since 9/11 not currently active duty (except reservs,ng?)
    deaf=factor(ifelse(dear==1,'deaf','hearing')),
    diss=ifelse(ddrs==1|deye==1|dout==1|dphy==1|drem==1,'disabled','nondisabled'),
    Age=ordered(
      ifelse(agep<25,'18-24',
        ifelse(agep<35,'25-34',
          ifelse(agep<45,'35-44',
            ifelse(agep<55,'45-54','55-64'))))),
    attainCum=ordered(
      ifelse(attain<'Regular high school diploma','No HS',
        ifelse(attain<'Some college, but less than 1 year','HS Diploma',
          ifelse(attain<'Associates degree','Some College',
            ifelse(attain<'Bachelors degree','Associates',
              ifelse(attain<'Masters degree','Bachelors','Post-Graduate'))))),
      levels=c('No HS','HS Diploma','Some College','Associates','Bachelors','Post-Graduate')),

    postSec= ordered(
      ifelse(attain<'Regular high school diploma','No HS',
            ifelse(attain<'Some college, but less than 1 year','HS',as.character(attain))),
      levels=c('No HS','HS','Some college, but less than 1 year',
        '1 or more years of college credit, no degree',
        'Associates degree',
        'Bachelors degree',
        'Masters degree',
        'Professional degree beyond a bachelors degree',
        'Doctorate degree')),
    postSec=ordered(
      ifelse(postSec<'Masters degree',as.character(postSec),'postBA'),
      levels=c('No HS','HS','Some college, but less than 1 year',
        '1 or more years of college credit, no degree',
        'Associates degree',
        'Bachelors degree',
        'postBA')),

    enrollment=ifelse(sch==1,'not enrolled',
      ifelse(schg<15,'k-12',
        ifelse(schg==15,'undergrad','grad school'))),

    enrolled=sch>1,

    nativeBorn=nativity==1,


    employment=factor(ifelse(esr%in%c(1,2,4,5),'Employed',
      ifelse(esr==3,'Unemployed','Not In Labor Force'))),


    fulltime=(employment=='Employed')&(wkw==1 & wkhp>=35),

    raceEth=ifelse(hisp>1,"Hispanic",
      ifelse(rac1p==2,"African American",
        ifelse(rac1p==6| rac1p==7,"Asian/PacIsl",
          ifelse(rac1p%in%c(3,4,5),'American Indian',
            ifelse(rac1p==1,"White","Other"))))),

    sex=ifelse(sex==1,'Male','Female'))
gc()

dat <- dat%>%filter(diss=='nondisabled')#,(deaf=='deaf'&dratx==1)|(deaf=='hearing'&dratx==2))

save(dat,file='vetSurvData.RData')
gc()
