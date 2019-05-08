output <- function(dat=tags,res,ressum=F,pltagger=T,anon=F,pltagtype=T,pldiag=F,plxprt=T,plcradle=T,condor=T,pltagger.type=F,
                plOTC=T,plcond=T,plqual=T,plsize=T,plassist=F,plevent=F,do_ci=T,prefix_nm="",io_agg=T,tagprog="iotc",ymaxtg=1.5,ymaxas=1.5,pltype="png") {
  # pltagtype = plot tag types
  # pldiag = plot diagnostics
  # pl xprt=0; ;1 ;2 
  # condor = T only if running on condor
  # tagger type if 
  
  
  if(ressum) a <- res else a <-summary(res)
  attributes(a)
  #a$coefficients[grep("tag_type",dimnames(a$coefficients)[[1]]),]
  a$coefficients[grep("tagger",dimnames(a$coefficients)[[1]]),]
  a$coefficients[grep("len5",dimnames(a$coefficients)[[1]]),]
  a$coefficients[grep("Cond",dimnames(a$coefficients)[[1]]),]
  a$coefficients[grep("Qual",dimnames(a$coefficients)[[1]]),]
  a$coefficients[grep("cradle",dimnames(a$coefficients)[[1]]),]
  a$coefficients[grep("OTC",dimnames(a$coefficients)[[1]]),]
  a$coefficients[grep("sp_id",dimnames(a$coefficients)[[1]]),]
  a$coefficients[grep("tag_type",dimnames(a$coefficients)[[1]]),]
  a$coefficients[grep("xprt",dimnames(a$coefficients)[[1]]),]
  a$coefficients[grep("xptype",dimnames(a$coefficients)[[1]]),]
  
  fn <- function(x) length(grep(x,names(res$model)))>0  
  inmodel <- data.frame(t(sapply(c("tagger","len5","Cond2","Qual","OTC","cradle","sp","tag_type","xprt","xptype","tagger.tagtype"),fn)))
  
  if(io_agg) tmp <- tapply(dat$wts,list(dat$tagger,dat$Qual),sum) else tmp <- tapply(dat$relrec[,1]+dat$relrec[,2],list(dat$tagger,dat$Qual),sum)
  tmp[is.na(tmp)==T] <- 0
  tmp <- tmp[apply(tmp,1,sum)!=0,]
  listqual <- cbind(apply(tmp,1,sum),tmp/apply(tmp,1,sum))
  dimnames(listqual)[[2]][1] <- "Tagged"
  write.csv(listqual,paste(prefix_nm,"listqual.csv",sep="_"))
  
  if(io_agg) tmp <- tapply(dat$wts,list(dat$tagger,dat$Cond),sum) else tmp <- tapply(dat$relrec[,1]+dat$relrec[,2],list(dat$tagger,dat$Cond),sum)
  tmp[is.na(tmp)==T] <- 0
  tmp <- tmp[apply(tmp,1,sum)!=0,]
  listcond <- cbind(apply(tmp,1,sum),tmp/apply(tmp,1,sum))
  dimnames(listcond)[[2]][1] <- "Tagged"
  write.csv(listcond,paste(prefix_nm,"listcond.csv",sep="_"))
  
  if(io_agg) tmp <- tapply(dat$wts,list(dat$tagger,dat$sp),sum) else tmp <-  cbind(tapply(dat$relrec[,1]+dat$relrec[,2],list(dat$tagger,dat$sp),sum))
  tmp[is.na(tmp)==T] <- 0
  tmp <- tmp[apply(tmp,1,sum)!=0,]
  listsp <- data.frame(cbind(apply(tmp,1,sum),tmp/apply(tmp,1,sum)))
  names(listsp)[1] <- c("Tagged")
  write.csv(listsp,paste(prefix_nm,"listsp.csv",sep="_"))

  if(io_agg) tmp <- tapply(dat$wts,list(dat$tagger,dat$cradle2),sum) else tmp <-  cbind(tapply(dat$relrec[,1]+dat$relrec[,2],list(dat$tagger,dat$cradle2),sum))
  tmp[is.na(tmp)==T] <- 0
  tmp <- tmp[apply(tmp,1,sum)!=0,]
  listsp <- data.frame(cbind(apply(tmp,1,sum),tmp/apply(tmp,1,sum)))
  names(listsp)[1] <- c("Tagged")
  write.csv(listsp,paste(prefix_nm,"list_cradle.csv",sep="_"))

  if(inmodel$tag_type) {
    if(io_agg) tmp <- tapply(dat$wts,list(dat$tagger,dat$tag_type),sum) else tmp <- tapply(dat$relrec[,1]+dat$relrec[,2],list(dat$tagger,dat$tag_type),sum)
    tmp[is.na(tmp)==T] <- 0
    tmp <- tmp[apply(tmp,1,sum)!=0,]
    list_tag_type <- data.frame(cbind(apply(tmp,1,sum),tmp/apply(tmp,1,sum)))
    names(list_tag_type)[1] <- c("Tagged")
    write.csv(list_tag_type,paste(prefix_nm,"list_tag_type.csv",sep="_"))
  
    if(io_agg) list_tag_type2 <- tapply(dat$wts,list(dat$tagger,dat$tag_type),sum) else list_tag_type2 <- tapply(dat$relrec[,1]+dat$relrec[,2],dat$tag_type,sum)
    list_tag_type2[is.na(list_tag_type2)==T] <- 0
    #names(list_tag_type2) <- c("Y13","Y11","ONG","GRN")
    write.csv(list_tag_type2,paste(prefix_nm,"list_tag_type2.csv",sep="_"))
  }
  
  # diagnostics
#  x11()
  if(pldiag==T) {
    par(mfrow=c(2,2))
    plot(res)
    if(condor==F) savePlot(paste(prefix_nm,"diagnostics",type=pltype))
  }
  
if(inmodel$tagger) basetagger=levels(dat$tagger)[1]
if(inmodel$tagger.tagtype & is.factor(dat$tagger.tagtype)) basetagger=levels(dat$tagger.tagtype)[1]
if(inmodel$tagger.tagtype & !is.factor(dat$tagger.tagtype)) basetagger="ADL.conv"
basett=levels(dat$tag_type)[1]
par(mfrow=c(1,1))
if(pltagger & !pltagger.type) plotit(a,prefix_nm=prefix_nm,parname="tagger",basenm=basetagger,xlb="Tagger",anon=T,condor=condor,ymax=ymaxtg,pltype=pltype)
if(pltagger & !pltagger.type) plotit(a,prefix_nm=prefix_nm,parname="tagger",basenm=basetagger,xlb="Tagger",anon=F,condor=condor,ymax=ymaxtg,pltype=pltype)
if(pltagger.type) plotit(a,prefix_nm=prefix_nm,parname="tagger.tagtype",basenm=basetagger,xlb="Tagger",anon=T,condor=condor,ymax=ymaxtg,pltype=pltype)  # will need to be checked and fixed
if(pltagger.type) plotit(a,prefix_nm=prefix_nm,parname="tagger.tagtype",basenm=basetagger,xlb="Tagger",anon=F,condor=condor,ymax=ymaxtg,pltype=pltype)  # will need to be checked and fixed
if(pltagtype) plotit(a,prefix_nm=prefix_nm,parname="tag_type",basenm=basett,xlb="Tag type",condor=condor,pltype=pltype)
if(plcond) {
  if(length(grep("Cond2",dimnames(a$coefficients)[[1]])) > 0) pn <- "Cond2" else pn <- "Cond"
  plotit(a,prefix_nm=prefix_nm,parname=pn,basenm="Good",xlb="Fish condition on release",condor=condor,pltype=pltype)
  }
if(plassist) plotit(a,prefix_nm=prefix_nm,parname="assist",basenm="CRW",xlb="Tagging assistant",anon=T,condor=condor,ymax=ymaxas,pltype=pltype)
if(plassist) plotit(a,prefix_nm=prefix_nm,parname="assist",basenm="CRW",xlb="Tagging assistant",anon=F,condor=condor,ymax=ymaxas,pltype=pltype)
if(plqual) plotit(a,prefix_nm=prefix_nm,parname="Qual",basenm="Good",xlb="Tag placement quality",condor=condor,pltype=pltype)
if(plOTC) plotit(a,prefix_nm=prefix_nm,parname="OTC",basenm="no OTC",xlb="OTC used",condor=condor,pltype=pltype)
if(plsize) plotsz(res=res,dat=dat,prefix_nm=prefix_nm,condor=condor,basenm=basetagger,do_ci=do_ci,pltype=pltype)
if(plxprt) plot_xprt(res=res,dat=dat,prefix_nm=prefix_nm,plxprt=plxprt,inmodel=inmodel,condor=condor,pltype=pltype)
if(plevent) plot_event(res=res,dat=dat,prefix_nm=prefix_nm,condor=condor,basenm=basetagger,pltype=pltype)
if(plcradle==T & inmodel$cradle) plotit(a,prefix_nm=prefix_nm,parname="cradle2",basenm="Port BOW",xlb="Cradle position",condor=condor,pltype=pltype)

  
  # plot for cradles
  if(plcradle==T & inmodel$cradle) {
    param <- "cradle2"; label <- "Tagging station"; baselev <- "Bow"
    cf.param <- data.frame(a$coefficients[grep(param,dimnames(a$coefficients)[[1]]),])
    brpos <- grep("2",unlist(strsplit(dimnames(cf.param)[[1]][1],NULL))) + 1
    dimnames(cf.param)[[1]] <- substr(dimnames(cf.param)[[1]],brpos,brpos+2)
    b <- cf.param[cf.param$Std..Error<2,]
    L <- b[,1]-1.96*b[,2]; U <- b[,1]+1.96*b[,2]
    labs <- dimnames(b)[[1]]
    par(mfrow=c(1,1))
    plot(b[,1],ylim=c(min(L),max(U)),xlab=label,ylab=paste("logit return rate (vs",baselev," at 0)"))
    X <- seq(1,length(b[,1]))
    segments(X,L,X,U)
    abline(h=0)
    text(X,b[,1]-0.1,labels=labs,cex=0.8)
    if(condor==F) savePlot(paste(prefix_nm,"cradle",sep="_"),type=pltype)
    pr <- (exp(b[,1])/(1+exp(b[,1])))/0.5
    Lp <- (exp(L)/(1+exp(L)))/0.5
    Up <- (exp(U)/(1+exp(U)))/0.5
    par(mfrow=c(1,1))
    plot(pr,xaxt='n',xlim=c(0.6,4.4),ylim=c(0,1.2),xlab=label,ylab=paste("Relative return rate (vs", baselev," at 1)"))
    text(X,pr-0.03,labels=labs,cex=0.75)
    segments(X,Lp,X,Up)
    abline(h=1)
    if(condor==F) savePlot(paste(prefix_nm,"cradle_probs",sep="_"),type=pltype)
  }
}

