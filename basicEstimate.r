load('../vetData.RData') # dat

options(stringsAsFactors=FALSE)
library(tidyverse)
## if(!exists("dat")){
##     if('vetData.RData'%in%list.files()){
##         load('vetData.RData')
##     } else source('makeData.r')
## }


source('../../generalCode/estimationFunctions.r')
source('../../generalCode/median.r')
source('../functions.r')

agerange <- paste(min(dat$agep),'-',max(dat$agep))
totaln <- nrow(dat)

everyone <- overall('postSecEnrolled',factorProps,dat)

dat <- filter(dat,vet=='vet')
overall18 <- overall('postSecEnrolled',factorProps,dat)
byDratx18 <- bysub('postSecEnrolled',factorProps,dat,dratx)
rates18 <- as.data.frame(t(byDratx18[c(1,6),-c(1:3)]))
names(rates18) <- paste(byDratx18[c(1,6),1],byDratx18[c(1,6),3],sep='DRATX=')
rates18s <- rates18[-grep(' SE',rownames(rates18)),]

overall25 <- overall('postSecEnrolled',factorProps,filter(dat,agep>24))
byDratx25 <- bysub('postSecEnrolled',factorProps,filter(dat,agep>24),dratx)
rates25 <- as.data.frame(t(byDratx25[c(1,6),-c(1:3)]))
names(rates25) <- paste(byDratx25[c(1,6),1],byDratx25[c(1,6),3],sep='DRATX=')
rates25s <- rates25[-grep(' SE',rownames(rates25)),]

save(everyone,overall18,byDratx18,overall25,byDratx25,file='proportions.RData')

openxlsx::write.xlsx(list(`18-54`=rates18s,`25-54`=rates25s),'proportions.xlsx',row.names=TRUE,col.names=TRUE)
