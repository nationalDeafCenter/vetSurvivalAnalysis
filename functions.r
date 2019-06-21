library(lmtest)
library(tidyverse)
transformData <- function(dat){
  dat$postSec <- droplevels(dat$postSec)
  dat$id <- 1:nrow(dat)
    ## time = stages of education: HS, begin college, 1 year college, AA, BA, ma, doctorate/prof
    ## currently enrolled=censored
    ## fixed effects for time
    ## yij=0 if attain_i>=j
    ## yij=1 if attain_i=j-1, not censored
  ## yij=0 for all j if censored
  levs <- levels(dat$postSec)
  nlevs <- nlevels(dat$postSec)
  for(ll in 2:nlevs){
    lev <- levs[ll]
    dat[[lev]] <- 0
    dat[[lev]][dat$postSec<levs[ll-1]] <- NA
    dat[[lev]][dat$postSec==levs[ll-1]&!dat$enrolled] <- 1
    dat[[lev]][dat$postSec==levs[ll-1]&dat$enrolled] <- NA
  }
  ldat <- gather(dat,"level","event",!!levs[2]:!!levs[nlevs],
    factor_key=TRUE)
  ldat <- filter(ldat,!is.na(event))
  ldat
}

## fit a DTSM
model <- function(formula=event~level-1,data,surveySEs=FALSE,returnReps=FALSE,cl=NULL){
  #mf <- model.frame(formula,data=data)
  mod <- glm(formula,data=data,#mf,
    family=binomial,weights=pwgtp)

  if(surveySEs){
    if(is.null(cl)){
      betas <- matrix(ncol=length(coef(mod)),nrow=80)
      for(i in 1:80){
        data$ww <- pmax(data[[paste0('pwgtp',i)]],0)
        betas[i,] <- coef(update(mod,weights=ww))
        if(i %in% seq(10,70,10)) cat(i,'/80 ',sep='')
      }
      cat('\n')
    } else{
      fn <- function(ww){
        mf$ww <- pmax(ww,0)
        coef(update(mod,weights=ww))
      }
      clusterExport(cl,c('mod','mf','fn'),envir=environment())
      betas <- parLapply(cl,data[,paste0('pwgtp',1:80)],fn)
      betas <- do.call('rbind',betas)

    }
    vvv <- crossprod(sweep(betas,2,coef(mod),'-')/20)
  } else vvv <- vcov(mod)

  mod$vcov <- vvv

  mod
}


plotBase <- function(mod,ldat){
  pdat <-
    tibble(
      level=factor(levels(ldat$level),
        levels=levels(ldat$level)),
      logitHazard=coef(mod)[paste0('level',level)],
      g=1)%>%
    mutate(hazard=plogis(logitHazard))%>%
    select(-logitHazard)
  pdat$attainment=sapply(1:nrow(pdat),function(i) prod(1-pdat$hazard[1:i]))
  pdat <- gather(pdat,"measure","est",-level,-g)

  ggplot(ungroup(pdat),aes(x=level,y=est))+geom_point()+geom_line(group=pdat$g)+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    facet_wrap(~measure,scales="free_y")
}

survivalHazard <- function(mod,ldat){
 pdat <-
    expand.grid(
      level=factor(levels(ldat$level),
        levels=levels(ldat$level)),
      deaf=factor(levels(ldat$deaf)))

  moreVars <- setdiff(names(get_all_vars(mod$formula,mod$data)),c('event',names(pdat)))
  if(length(moreVars)){
    for(vv in moreVars) pdat[[vv]] <-
                          if(is.factor(mod$data[[vv]])){
                            levels(mod$data[[vv]])[1]
                          } else if(is.character(mod$data[[vv]])){
                            sort(unique(mod$data[[vv]]))[1]
                          } else if(is.logical(mod$data[[vv]])) FALSE else 0
  }
  pdat$logitHazard <- predict(mod,pdat)#,type='link')
  pdat$hazard <- plogis(pdat$logitHazard)
  pdat%>%group_by(deaf)%>%
    mutate(attainment=sapply(1:n(),function(i) prod(1-hazard[1:i])))%>%
    select(-logitHazard)%>%ungroup()
}

plotInt <- function(mod,ldat,reverse=FALSE){

  pdat <- survivalHazard(mod,ldat)

  if(reverse){
    pdat$advanceProb <- 1-pdat$hazard
    pdat <- gather(pdat,"measure","est",advanceProb,attainment)
  } else pdat <- gather(pdat,"measure","est",hazard,attainment)

  ggplot(pdat,aes(x=level,y=est,group=deaf,color=deaf))+geom_point()+geom_line()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    facet_wrap(~measure,scales="free_y")
}




plotDiff <- function(mod,ldat){
  pdat <-
    expand.grid(
      level=factor(levels(ldat$level),
        levels=levels(ldat$level)),
      deaf=factor(levels(ldat$deaf)))
  moreVars <- setdiff(names(mod$data),c('event',names(pdat)))
  if(length(moreVars)){
    for(vv in moreVars) pdat[[vv]] <-
                          if(is.factor(mod$data[[vv]])){
                            levels(mod$data[[vv]])[1]
                          } else if(is.character(mod$data[[vv]])){
                            sort(unique(mod$data[[vv]]))[1]
                          } else 0
  }
  pdat$logitHazard <- predict(mod,pdat,type='link')
  pdat$hazard <- plogis(pdat$logitHazard)
  pdat <- pdat%>%group_by(level)%>%
    summarize(diffHazard=hazard[deaf=='deaf']-hazard[deaf=='hearing'],g=1)

  ggplot(pdat,aes(level,diffHazard,group=g))+geom_point()+geom_line()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
}

varSum <- function(index,vcv){
  index <- ifelse(seq(nrow(vcv))%in%index,1,0)
  (t(index)%*%vcv%*%index)[1,1]
}

hazardOddsRatios <- function(mod,ldat,age=0){
  rownames(mod$vcov) <- colnames(mod$vcov) <- names(coef(mod))
  diffs <- list()
  vars <- list()
  mm <- model.matrix(mod)
  agePred <- mm[mod$dat$ageCentered==age-35,grep('deafhearing:ns',colnames(mm),fixed=TRUE)][1,]
  ageNum <- grep('deafhearing:ns',colnames(mm),fixed=TRUE)

  vcv <- mod$vcov

  for(i in 1:length(ageNum)){
    vcv[ageNum[i],] <- vcv[ageNum[i],]*agePred[i]
    vcv[,ageNum[i]] <- vcv[,ageNum[i]]*agePred[i]
  }

  dh <- which(names(coef(mod))=='deafhearing')
  diffs[[levels(ldat$level)[1]]] <- (coef(mod)[dh]+crossprod(agePred,coef(mod)[names(agePred)]))[1,1]
  vars[[levels(ldat$level)[1]]] <- varSum(c(dh,ageNum),vcv)

  for(ll in levels(ldat$level)[-1]){
    nn <- paste0('level',ll,':deafhearing')
    diffs[[ll]] <- diffs[[1]]+coef(mod)[nn]
    vars[[ll]] <- varSum(c(dh,ageNum,which(names(coef(mod))==nn)),vcv)
  }

  cis <- mapply(FUN=function(d,v) exp(-d+c(2,-2)*sqrt(v)),diffs,vars)
  diffs <- sapply(-unlist(diffs),exp)


  pdat <- data.frame(
    level=factor(levels(ldat$level),levels=levels(ldat$level)),
    hazardRatio=diffs,
    ciL=cis[1,],
    ciH=cis[2,],
    gr=1)
  ggplot(pdat,aes(level,hazardRatio,group=gr))+
    geom_point()+geom_line()+
    geom_errorbar(aes(ymin=ciL,ymax=ciH),width=0.1)+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    xlab(NULL)+ylab('Hazard Odds Ratio deaf/hearing')+
    geom_hline(yintercept=1)
}



plotAge <- function(mod){
  pdat <- expand.grid(
    ageCentered=seq(min(mod$data$ageCentered),max(mod$data$ageCentered),.1),
    deaf=levels(mod$data$deaf)
  )
  moreVars <- setdiff(names(get_all_vars(mod$formula,mod$data)),c('event',names(pdat)))
  if(length(moreVars)){
    for(vv in moreVars) pdat[[vv]] <-
                          if(is.factor(mod$data[[vv]])){
                            levels(mod$data[[vv]])[1]
                          } else if(is.character(mod$data[[vv]])){
                            sort(unique(mod$data[[vv]]))[1]
                          } else if(is.logical(mod$data[[vv]])) FALSE else 0
  }
  pdat$ageEff <- predict(mod,pdat)
  pdat$ageEff[pdat$deaf=='hearing'] <- pdat$ageEff[pdat$deaf=='hearing']-coef(mod)['deafhearing']
  pdat$age=pdat$ageCentered+35
  ggplot(pdat,aes(age,ageEff,color=deaf))+geom_line()
}
