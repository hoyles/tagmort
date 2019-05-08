plotit <- function(x,prefix_nm,parname,basenm,xlb="Tagger",anon=F,condor,ymax=1.5,pltype="png") {
  #if(parname=="cradle2") browser()
  pos <- grep(parname,dimnames(x$coefficients)[[1]])
  namelen <- nchar(parname)
  if(length(pos) >1) {
    cf <- data.frame(x$coefficients[pos,])
    dimnames(cf)[[1]] <- substr(dimnames(cf)[[1]],namelen+1,namelen+20)
    } else {
    cf <- data.frame(rbind(x$coefficients[pos,],c(0,2,0,0)))
    dimnames(cf)[[1]][1] <- substr(dimnames(x$coefficients)[[1]][pos],namelen+1,namelen+20)
    }
  b <- cf[cf$Std..Error<2,]
  L <- b[,1]-1.96*b[,2]; U <- b[,1]+1.96*b[,2]
  labs <- dimnames(b)[[1]]
  if(anon) basenm <- "Reference"
  if(anon & parname%in% c("tagger","assist")) labs <- seq(2:(length(dimnames(b)[[1]])+1))
  if(anon & parname=="tagger.tagtype") {
    labs <- dimnames(b)[[1]]
    k <- 0
    for(i in 1:length(labs)) {
      if(substring(labs[i],5,8)=="conv") k <- k + 1
      labs[i] <- paste(k,substring(labs[i],5,10),sep=".")
      }
    }
  #if(io_agg) b_tagger <- "JDD" else b_tagger <- "Tony Lewis"
  if(parname %in% c("tagger.tagtype","assist")) windows(width=9,height=7) else windows()
  plot(b[,1],ylim=c(min(L),max(U)),xlab=xlb,ylab=paste("logit Return rate (vs",basenm,"at 0)"))
  X <- seq(1,length(b[,1]))
  segments(X,L,X,U)
  abline(h=0)
  text(X,b[,1]-0.1,labels=labs,cex=0.8)
  if(anon) prefix_nm <- paste(prefix_nm,"_anon",sep="")
  if(condor==F) savePlot(paste(prefix_nm,parname,"logit",sep="_"),type=pltype)
  pr <- (exp(b[,1])/(1+exp(b[,1])))/0.5
  L <- (exp(L)/(1+exp(L)))/0.5
  U <- (exp(U)/(1+exp(U)))/0.5
  plot(pr,xlim=c(0.5,length(pr)+0.5),ylim=c(0,ymax),xaxt="n",xlab=xlb,ylab=paste("Relative Return rate (vs",basenm,"at 1)"))
  text(X,pr-0.03,labels=labs,cex=0.75)
  segments(X,L,X,U)
  abline(h=1)
  pars <- cbind(b,pr,L,U)
  if(anon & parname %in% c("tagger","tagger.tagtype","assist")) rownames(pars) <- labs
  write.csv(pars,file=paste(prefix_nm,parname,"parests.csv",sep="_"))
  if(condor==F) savePlot(paste(prefix_nm,parname,"probs",sep="_"),type=pltype)
  }

plot_tagger.tagtype <- function(x,prefix_nm,parname,basenm,xlb="Tagger",anon=F,condor,ymax=1.5,pltype="png") {
  pos <- grep(parname,dimnames(x$coefficients)[[1]])
  namelen <- nchar(parname)
  cf <- data.frame(x$coefficients[pos,])
  dimnames(cf)[[1]] <- substr(dimnames(cf)[[1]],namelen+1,namelen+20)
  b <- cf[cf$Std..Error<2,]
  L1a <- b[,1]-1.96*b[,2]; L2a <- b[,1] - 0.01
  U1a <- b[,1]+1.96*b[,2]; U2a <- b[,1] + 0.01
  labs <- dimnames(b)[[1]]
  k <- 0
  for(i in 1:length(labs)) {
    if(substring(labs[i],5,8)=="conv") k <- k + 1
    labs[i] <- paste(k,substring(labs[i],5,10),sep=".")
    }
  numlabs <- gsub(".intern","",gsub(".conv","",labs))
  windows(width=9,height=7)
  ttype <- rep(7,length(labs))
  ttype[grep("conv",labs)] <- 1
  ltype <- ttype
  ltype[ltype==7] <- 2
  X <- seq(1,length(b[,1]))
  if(anon) prefix_nm <- paste(prefix_nm,"_anon",sep="")
  pr <- (exp(b[,1])/(1+exp(b[,1])))/0.5
  L1 <- (exp(L1a)/(1+exp(L1a)))/0.5; L2 <- pr - 0.01
  U1 <- (exp(U1a)/(1+exp(U1a)))/0.5; U2 <- pr + 0.01
  plot(pr,xlim=c(0.5,length(pr)+0.5),ylim=c(0,ymax),xaxt="n",xlab=xlb,ylab=paste("Relative Return rate (vs Reference at 1)"),pch=ttype,cex=1)
  text(X,L1-0.03,labels=numlabs,cex=1)
  segments(X,U1,X,U2,lty=ltype)
  segments(X,L1,X,L2,lty=ltype)
  abline(h=1)
  legend("topleft",legend=c("Conventional tag","Internal tag"),pch=c(1,7),lty=c(1,2))
  }

plot_tagger <- function(x,prefix_nm,parname,basenm,xlb="Tagger",anon=F,condor,ymax=1.5,pltype="png") {
  pos <- grep(parname,dimnames(x$coefficients)[[1]])
  namelen <- nchar(parname)
  cf <- data.frame(x$coefficients[pos,])
  dimnames(cf)[[1]] <- substr(dimnames(cf)[[1]],namelen+1,namelen+20)
  b <- cf[cf$Std..Error<2,]
  L1a <- b[,1]-1.96*b[,2]; L2a <- b[,1] - 0.01
  U1a <- b[,1]+1.96*b[,2]; U2a <- b[,1] + 0.01
  labs <- dimnames(b)[[1]]
  k <- 0
  for(i in 1:length(labs)) {
    labs[i] <- i
    }
  numlabs <- gsub(".intern","",gsub(".conv","",labs))
  windows()
  ttype <- rep(2,length(labs))
  ttype[grep("conv",labs)] <- 2
  ltype <- 1
  ltype[ltype==7] <- 2
  X <- seq(1,length(b[,1]))
  if(anon) prefix_nm <- paste(prefix_nm,"_anon",sep="")
  pr <- (exp(b[,1])/(1+exp(b[,1])))/0.5
  L1 <- (exp(L1a)/(1+exp(L1a)))/0.5; L2 <- pr - 0.01
  U1 <- (exp(U1a)/(1+exp(U1a)))/0.5; U2 <- pr + 0.01
  plot(pr,xlim=c(0.5,length(pr)+0.5),ylim=c(0,ymax),xaxt="n",xlab=xlb,ylab=paste("Relative Return rate (vs Reference at 1)"),pch=ttype,cex=1)
  text(X,L1-0.03,labels=numlabs,cex=1)
  segments(X,U1,X,U2,lty=ltype)
  segments(X,L1,X,L2,lty=ltype)
  abline(h=1)
  }

plot_event <- function(res,dat=indat,prefix_nm,condor,basenm,pltype="png") {
  nms <- names(dat)
  nms <- nms[!nms %in% c("sp_id","len5","wts")]
  nmlist <- as.list(dat[1,nms])
  if(length(grep("xprt",nms))>0) nmlist$xprt <- 2500
  if(length(grep("assist",nms))>0) nmlist$assist <- levels(dat$assist)[1]
  if(length(grep("cradle2",nms))>0) nmlist$cradle2 <- levels(dat$cradle2)[1]
  if(length(grep("tagger.tagtype",nms))>0) {
    nmlist$tagger.tagtype <- basenm
    } else  nmlist$tagger <- basenm
  nmlist$tag_sch_id <- unique(dat$tag_sch_id)
  nmlist$len5 <- 50
  nmlist$Qual <- levels(dat$Qual)[1]
  spp <- sort(unique(dat$sp_id))
  if(length(spp) > 1) nmlist$sp_id <- spp[2] else nmlist$sp_id <- spp
  pd <- expand.grid(nmlist)
  p <- predict(res,newdata=pd,type="response",se.fit=F)
  windows()
  hist(p,breaks=seq(0,1,0.01),xlim=c(0,1),xlab="Return rate",ylab="Frequency",main="")
  if(condor==F) savePlot(paste(prefix_nm,"event",sep="_"),type=pltype)
  }


plotsz <- function(res,dat=indat,prefix_nm,xlb="Length",condor,basenm,do_ci=T,pltype="png") {
if(length(unique(dat$sp_id))>1) {
  nms <- names(dat)
  nms <- nms[!nms %in% c("sp_id","len5","wts")]
  nmlist <- as.list(dat[1,nms])
  if(length(grep("tagger.tagtype",nms))>0) {
    nmlist$tagger.tagtype <- basenm
    } else nmlist$tagger <- basenm
  if(length(grep("assist",nms))>0) nmlist$assist <- levels(dat$assist)[1]
  if(length(grep("OTC",nms))>0) nmlist$OTC <- levels(dat$OTC)[1]
  if(length(grep("xprt",nms))>0) nmlist$xprt <- 2500
  if(length(grep("cradle2",nms))>0) nmlist$cradle2 <- levels(dat$cradle2)[1]
  if(nmlist$cradle2=="ARC") nmlist$cradle2 <- levels(dat$cradle2)[2]

  nmlist$len5 <- sort(unique(dat$len5))
  nmlist$Qual <- levels(dat$Qual)[1]
  nmlist$sp_id <- unique(dat$sp_id)
  pd <- expand.grid(nmlist)
  pd <- pd[paste(pd$sp_id,pd$len5) %in% unique(paste(dat$sp_id,dat$len5)),]
  par(mfrow=c(2,2))
  p <- predict(res,newdata=pd,type="link",se.fit=T)
  pr <- inv.logit(cbind(p$fit,p$fit + 2*p$se.fit,p$fit - 2*p$se.fit))
  if(do_ci) yl <- c(0,max(pr)) else yl <- c(0,max(pr[,1]))
  if(length(grep("S",dat$sp_id)) > 0) {
    plot(pd[pd$sp_id=="S",]$len5,pr[pd$sp_id=="S",1],ylim=yl,type="l",xlab="Skipjack length (cm)",ylab="Return rate")
    if(do_ci) lines(pd[pd$sp_id=="S",]$len5,pr[pd$sp_id=="S",2])
    if(do_ci) lines(pd[pd$sp_id=="S",]$len5,pr[pd$sp_id=="S",3])
    }
  if(length(grep("B",dat$sp_id)) > 0) {
    plot(pd[pd$sp_id=="B",]$len5,pr[pd$sp_id=="B",1],ylim=yl,type="l",xlab="Bigeye length (cm)",ylab="Return rate")
    if(do_ci) lines(pd[pd$sp_id=="B",]$len5,pr[pd$sp_id=="B",2])
    if(do_ci) lines(pd[pd$sp_id=="B",]$len5,pr[pd$sp_id=="B",3])
  }
  if(length(grep("Y",dat$sp_id)) > 0) {
    plot(pd[pd$sp_id=="Y",]$len5,pr[pd$sp_id=="Y",1],ylim=yl,type="l",xlab="Yellowfin length (cm)",ylab="Return rate")
    if(do_ci) lines(pd[pd$sp_id=="Y",]$len5,pr[pd$sp_id=="Y",2])
    if(do_ci) lines(pd[pd$sp_id=="Y",]$len5,pr[pd$sp_id=="Y",3])
  }
  if(condor==F) savePlot(paste(prefix_nm,"length",sep="_"),type=pltype)
  } else {
  
  nms <- names(dat)
  nms <- nms[!nms %in% c("len5","wts")]
  nmlist <- as.list(dat[1,nms])
  nmlist$len5 <- sort(unique(dat$len5))
  if(length(grep("tagger.tagtype",nms))>0) {
    nmlist$tagger.tagtype <- basenm
    } else nmlist$tagger <- basenm
  if(length(grep("assist",nms))>0) nmlist$assist <- levels(dat$assist)[1]
  if(length(grep("OTC",nms))>0) nmlist$OTC <- levels(dat$OTC)[1]
  if(length(grep("xprt",nms))>0) nmlist$xprt <- 2500
  if(length(grep("cradle2",nms))>0) nmlist$cradle2 <- levels(dat$cradle2)[1]
  nmlist$Qual <- levels(dat$Qual)[1]
  spp <- sort(unique(dat$sp_id))
  if(length(spp) > 1) nmlist$sp_id <- spp[2] else nmlist$sp_id <- spp

  pd <- expand.grid(nmlist)
  p <- predict(res,newdata=pd,type="response",se.fit=T)
  p2 <- predict(res,newdata=pd,type="terms",se.fit=T)
  plot(pd$len5,p$fit,ylim=c(0,max(p$fit+2*p$se.fit)),type="b",xlab="Length (cm)",ylab="Return rate")
  sd2 <- p2$se.fit[,2]
#  lines(pd$len5,p$fit+2*p$se.fit)
#  lines(pd$len5,p$fit-2*p$se.fit)
  lines(pd$len5,inv.logit(logit(p$fit) - 1.96*sd2))
  lines(pd$len5,inv.logit(logit(p$fit) + 1.96*sd2))
  if(condor==F) savePlot(paste(prefix_nm,"length",sep="_"),type=pltype)
  }}
  
plotsz2 <- function(res,dat=indat,prefix_nm,xlb="Length",condor,pltype="png") {
if(length(unique(dat$sp_id))>1) {
  nms <- names(dat)
  nms <- nms[!nms %in% c("sp_id","len5","wts")]
  nmlist <- as.list(dat[1,nms])
  nmlist$len5 <- sort(unique(dat$len5))
  nmlist$sp_id <- unique(dat$sp_id)
  pd <- expand.grid(nmlist)
  pd <- pd[paste(pd$sp_id,pd$len5) %in% unique(paste(dat$sp_id,dat$len5)),]
  par(mfrow=c(2,2))
  p <- predict(res,newdata=pd,type="response",se.fit=F)
  pter <- predict(res,newdata=pd,type="terms",se.fit=T)
  loc <- grep(":sp_id",colnames(pter$fit))
  plgt <- logit(p)
  pr_adj <- inv.logit(cbind(plgt,plgt + 2*pter$se.fit[,loc],plgt - 2*pter$se.fit[,loc]))
  plot(pd[pd$sp_id=="S",]$len5,pr_adj[pd$sp_id=="S",1],ylim=c(0,max(pr_adj)),type="b",xlab="Skipjack length (cm)",ylab="Return rate")
  lines(pd[pd$sp_id=="S",]$len5,pr_adj[pd$sp_id=="S",2])
  lines(pd[pd$sp_id=="S",]$len5,pr_adj[pd$sp_id=="S",3])
  plot(pd[pd$sp_id=="B",]$len5,pr_adj[pd$sp_id=="B",1],ylim=c(0,max(pr_adj)),type="b",xlab="Bigeye length (cm)",ylab="Return rate")
  lines(pd[pd$sp_id=="B",]$len5,pr_adj[pd$sp_id=="B",2])
  lines(pd[pd$sp_id=="B",]$len5,pr_adj[pd$sp_id=="B",3])
  plot(pd[pd$sp_id=="Y",]$len5,pr_adj[pd$sp_id=="Y",1],ylim=c(0,max(pr_adj)),type="b",xlab="Yellowfin length (cm)",ylab="Return rate")
  lines(pd[pd$sp_id=="Y",]$len5,pr_adj[pd$sp_id=="Y",2])
  lines(pd[pd$sp_id=="Y",]$len5,pr_adj[pd$sp_id=="Y",3])
  if(condor==F) savePlot(paste(prefix_nm,"length",sep="_"),type=pltype)
  } else {
  nms <- names(dat)
  nms <- nms[!nms %in% c("len5","wts")]
  nmlist <- as.list(dat[1,nms])
  nmlist$len5 <- sort(unique(dat$len5))
  pd <- expand.grid(nmlist)
  p <- predict(res,newdata=pd,type="response",se.fit=T)
  plot(pd$len5,p$fit,ylim=c(0,max(p$fit+2*p$se.fit)),type="b",xlab="Length (cm)",ylab="Return rate")
  lines(pd$len5,p$fit+2*p$se.fit)
  lines(pd$len5,p$fit-2*p$se.fit)
  if(condor==F) savePlot(paste(prefix_nm,"length",sep="_"),type=pltype)
  }}
  
plot_xprt<- function(res,dat=indat,prefix_nm,plxprt,inmodel,condor,pltype="png") {
  if(plxprt & inmodel$xprt & inmodel$xptype) {
    nms <- names(dat)
    nms <- nms[!nms %in% c("xprt","xptype","wts")]
    nmlist <- as.list(dat[1,nms])
    nmlist$xprt <- sort(unique(dat$xprt))
    nmlist$xptype <- unique(dat$xptype)
    pd <- expand.grid(nmlist)
    pd <- pd[paste(pd$xptype,pd$xprt) %in% unique(paste(dat$xptype,dat$xprt)),]
    par(mfrow=c(2,2))
    p <- predict(res,newdata=pd,type="response",se.fit=F)
    pter <- predict(res,newdata=pd,type="terms",se.fit=T)
    loc <- grep("xptype:",colnames(pter$fit))
    plgt <- logit(p)
    maxp <- max(plgt)
    pr_adj <- inv.logit(cbind(plgt,plgt + 2*pter$se.fit[,loc],plgt - 2*pter$se.fit[,loc]))
    pr_adj[pd$xptype=="N",] <- pr_adj[pd$xptype=="N",] / max(pr_adj[pd$xptype=="N",1])
    pr_adj[pd$xptype=="O",] <- pr_adj[pd$xptype=="O",] / max(pr_adj[pd$xptype=="O",1])
    pr_adj[pd$xptype=="T",] <- pr_adj[pd$xptype=="T",] / max(pr_adj[pd$xptype=="T",1])
    plot(pd[pd$xptype=="N",]$xprt,pr_adj[pd$xptype=="N",1],ylim=c(0,max(pr_adj)),type="b",xlab="Tuna tagged",ylab="Relative return rate")
    lines(pd[pd$xptype=="N",]$xprt,pr_adj[pd$xptype=="N",2])
    lines(pd[pd$xptype=="N",]$xprt,pr_adj[pd$xptype=="N",3])
    plot(pd[pd$xptype=="O",]$xprt,pr_adj[pd$xptype=="O",1],ylim=c(0,max(pr_adj)),type="b",xlab="Tuna tagged",ylab="Relative return rate")
    lines(pd[pd$xptype=="O",]$xprt,pr_adj[pd$xptype=="O",2])
    lines(pd[pd$xptype=="O",]$xprt,pr_adj[pd$xptype=="O",3])
    plot(pd[pd$xptype=="T",]$xprt,pr_adj[pd$xptype=="T",1],ylim=c(0,max(pr_adj)),type="b",xlab="Tuna tagged",ylab="Relative return rate")
    lines(pd[pd$xptype=="T",]$xprt,pr_adj[pd$xptype=="T",2])
    lines(pd[pd$xptype=="T",]$xprt,pr_adj[pd$xptype=="T",3])
    if(condor==F) savePlot(paste(prefix_nm,"xprt_xptype",sep="_"),type=pltype)
    } else if(plxprt & inmodel$xprt & !inmodel$xptype) {
    nms <- names(dat)
    nms <- nms[!nms %in% c("xprt","wts")]
    nmlist <- as.list(dat[1,nms])
    nmlist$xprt <- sort(unique(dat$xprt))
    pd <- expand.grid(nmlist)
    pd <- pd[paste(pd$xprt) %in% unique(paste(dat$xprt)),]
    par(mfrow=c(1,1))
    p <- predict(res,newdata=pd,type="response",se.fit=F)
    pter <- predict(res,newdata=pd,type="terms",se.fit=T)
    loc <- grep("xprt",colnames(pter$fit))
    plgt <- logit(p)
    pr_adj <- inv.logit(cbind(plgt,plgt + 2*pter$se.fit[,loc],plgt - 2*pter$se.fit[,loc]))
    pr_adj <- pr_adj / max(pr_adj[,1])
    plot(pd$xprt,pr_adj[,1],ylim=c(0,max(pr_adj)),type="b",xlab="Tuna tagged",ylab="Relative return rate")
    lines(pd$xprt,pr_adj[,2])
    lines(pd$xprt,pr_adj[,3])
    if(condor==F) savePlot(paste(prefix_nm,"xprt",sep="_"),type=pltype)
   } else if(plxprt==2) {
    cf.xprt <- data.frame(a$coefficients[grep("xprt",dimnames(a$coefficients)[[1]]),])
    dimnames(cf.xprt)[[1]] <- substr(dimnames(cf.xprt)[[1]],16,44)
    b <- cf.xprt[cf.xprt$Std..Error<2,]
    L <- b[,1]-1.96*b[,2]; U <- b[,1]+1.96*b[,2]
    labs <- dimnames(b)[[1]]
    plot(b[,1],xlim=c(0,3),ylim=c(min(L),max(c(0,max(U)))),xlab="Expertise",ylab="logit Return rate (vs Good at 0)")
    X <- seq(1,length(b[,1]))
    segments(X,L,X,U)
    abline(h=0)
    text(X,b[,1]-0.1,labels=labs,cex=0.8)
    }
  if(condor==F) savePlot(paste(prefix_nm,"xprt",sep="_"),type=pltype)
  }


aggdata <- function(dat=tag_all,minTperSch=50,archival=T,xptype=F,tagger.tagtype=F,assist=F,cond=F,cradle=F) {
  if(!cond & !cradle) {
    if(!tagger.tagtype & !xptype & !assist) tag_agg <- aggregate(recap ~ eval(recap==T) + tagger + len5 + sp_id + Qual + Cond + Cond2 + OTC + tag_type + xprt + tag_sch_id,data=dat,length)
    if(!tagger.tagtype & xptype & !assist) tag_agg <- aggregate(recap ~ eval(recap==T) + tagger + len5 + sp_id + Qual + Cond + Cond2 + OTC + tag_type + xprt + tag_sch_id + xptype,data=dat,FUN=length)
    if(tagger.tagtype & !xptype & !assist) tag_agg <- aggregate(recap ~ eval(recap==T) + tagger.tagtype + len5 + sp_id + Qual + Cond + Cond2 + OTC + xprt + tag_sch_id,data=dat,length)
    if(tagger.tagtype & xptype & !assist) tag_agg <- aggregate(recap ~ eval(recap==T) + tagger.tagtype + len5 + sp_id + Qual + Cond + Cond2 + OTC + xprt + tag_sch_id + xptype,data=dat,FUN=length)
    if(!tagger.tagtype & !xptype & assist) tag_agg <- aggregate(recap ~ eval(recap==T) + tagger + len5 + sp_id + Qual + Cond + Cond2 + OTC + tag_type + xprt + assist + tag_sch_id,data=dat,length)
    if(!tagger.tagtype & xptype & assist) tag_agg <- aggregate(recap ~ eval(recap==T) + tagger + len5 + sp_id + Qual + Cond + Cond2 + OTC + tag_type + xprt + assist + tag_sch_id + xptype,data=dat,FUN=length)
    if(tagger.tagtype & !xptype & assist) tag_agg <- aggregate(recap ~ eval(recap==T) + tagger.tagtype + len5 + sp_id + Qual + Cond + Cond2 + OTC + xprt + assist + tag_sch_id,data=dat,length)
    if(tagger.tagtype & xptype & assist) tag_agg <- aggregate(recap ~ eval(recap==T) + tagger.tagtype + len5 + sp_id + Qual + Cond + Cond2 + OTC + xprt + assist + tag_sch_id + xptype,data=dat,FUN=length)
    }
  if(cond==1 & !cradle) {
    if(tagger.tagtype & !xptype & assist) tag_agg <- aggregate(recap ~ eval(recap==T) + tagger.tagtype + len5 + sp_id + Qual + Cond + OTC + xprt + assist + tag_sch_id,data=dat,length)
    if(!tagger.tagtype & !xptype & assist) tag_agg <- aggregate(recap ~ eval(recap==T) + tagger + len5 + sp_id + Qual + Cond + OTC + tag_type + xprt + assist + tag_sch_id,data=dat,length)
    }
  if(cond==1 & cradle) {
    if(tagger.tagtype & !xptype & assist) tag_agg <- aggregate(recap ~ eval(recap==T) + tagger.tagtype + len5 + sp_id + Qual + Cond + cradle2 + xprt + assist + tag_sch_id,data=dat,length)
    if(!tagger.tagtype & !xptype & assist) tag_agg <- aggregate(recap ~ eval(recap==T) + tagger + len5 + sp_id + Qual + Cond + cradle2 + tag_type + xprt + assist + tag_sch_id,data=dat,length)
    }
  if(!cond & cradle) {
    if(tagger.tagtype & !xptype & !assist) tag_agg <- aggregate(recap ~ eval(recap==T) + tagger.tagtype + len5 + sp_id + Qual + Cond + Cond2 + OTC + xprt + cradle2 + tag_sch_id,data=dat,length)
    if(!tagger.tagtype & !xptype & !assist) tag_agg <- aggregate(recap ~ eval(recap==T) + tagger + len5 + sp_id + Qual + Cond + Cond2 + OTC + tag_type + xprt + cradle2 + tag_sch_id,data=dat,length)
    }

  names(tag_agg)[c(1,3)] <- c("recap","len5")
  names(tag_agg)[grep("recap",names(tag_agg))[2]] <- c("wts")
  if(!archival) tag_agg <- tag_agg[tag_agg$tag_type %in% c("ST","DT","a.Y13","b.Y11"),]
  a <- tapply(tag_agg$wts,tag_agg$tag_sch_id,sum,na.rm=T)
  a <- a[a > minTperSch]
  tagsub <- tag_agg[tag_agg$tag_sch_id %in% names(a),]
  tagsub$tag_sch_id <- as.factor(as.character(tagsub$tag_sch_id))
  tagsub <- tagsub[tagsub$sp_id != "U",]
  tagsub <- tagsub[tagsub$Qual != "Unknown",]
  tagsub$Qual <- as.factor(as.character(tagsub$Qual))
  tagsub$Qual <- relevel(tagsub$Qual,ref="Good")
  tagsub$sp_id <- as.factor(as.character(tagsub$sp_id))
  if(!tagger.tagtype) {
    tagsub$tag_type <- as.factor(as.character(tagsub$tag_type))
    if(!is.na(match("ST",levels(tagsub$tag_type)))) tagsub$tag_type <- relevel(tagsub$tag_type,ref="ST") else tagsub$tag_type <- relevel(tagsub$tag_type,ref="a.Y13")
    }
  if(length(grep("OTC",names(tagsub))) > 0) tagsub$OTC <- as.factor(tagsub$OTC)
  tagsub$recap <- na.omit(tagsub$recap)
  return(tagsub)
  }

plot_assocs <- function(model) {
  a <- summary(model)@coefs
  coefs <- a[grep("assoc",rownames(a)),]
  passc <- exp(coefs)
  nass <- dim(coefs)[1]+1
  assc <- factor(1:nass,levels=1:nass,labels=c("Unassociated","Log","Anchored FAD","Drifting FAD")[1:nass])
  plotCI(1:nass,exp(c(0,coefs[,1])),ui=exp(c(0,coefs[,1]+2*coefs[,2])), li=exp(c(0,coefs[,1]-2*coefs[,2])),sfrac=0.01, col=1,lwd=1,xaxt="n",xlab="",
      ylab="Relative recapture rate",ylim=c(0,max(c(exp(coefs[,1]+2*coefs[,2])),1)))
  abline(h=1,lty=2)
  axis(1, at = 1:nass, labels = assc)
  }

