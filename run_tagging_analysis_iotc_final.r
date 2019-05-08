#dir <- "Desktop/tagging"
dir<- "~/../SkyDrive/Work/Papers_finished/Tagger_effects/iotc"
setwd(dir)
require(splines)
require(lme4)
require(pkpkg)
require(plotrix)
require(boot)
require(car)
require(MASS)
require(R4MFCL)
ddrive <- "D:/assessments/tagging/tagger_fx/TagFX_paper/"
#ddrive <- "P:/tagging/tagger_fx/TagFX_paper/"
dd_iotc <- ppath(ddrive,"iotc")

#"I:\assessments/tagging/Rep22_Q1_by_sch_simonh (02_12_2008).csv"

source("../R_scripts/load_iotc.r")
source("../R_scripts/prepare_data_both.r")
source("../R_scripts/output.r")
source("../R_scripts/support_functions.r")

infile <- "TAG_REC_AB3.txt"
iotags <- load_iotc(infile)
iotag_all <- prepare_data(iotags,maxprt=8000,xpf="TaggerExp.csv",reg="iotc")
dat <- iotagsub <- aggdata(dat=iotag_all,minTperSch=20,archival=T,xptype=F,tagger.tagtype=F,cond=F,cradle=T); str(iotagsub)
dat$cradle2[dat$cradle2=="ARC"] <- unique(dat$cradle2)[1]
dat$cradle2 <- as.factor(as.character(dat$cradle2))
levels(dat$cradle2)
tapply(iotagsub$wts,iotagsub$cradle2,sum,na.rm=T)
tapply(iotagsub$wts,iotagsub$OTC,sum,na.rm=T)
tapply(iotag_all$OTC,list(iotag_all$OTC,iotag_all$sp_id),length)
###################

AIClist <- vector()
mod_att <- glm(recap ~ tagger + ns(len5,df=4)*sp_id + Qual + Cond2 + OTC + tag_type + tag_sch_id + cradle2,data=dat,weights=wts,family=binomial("logit")) # originally mod_noa (wrong name)
summary(mod_att)
save(mod_att,file=ppath(ddrive,"iotc/iotc_mod2.RData"))

mod_noa_B <- glm(recap ~ tagger + ns(len5,df=4) + Qual + Cond2 + OTC + tag_type + tag_sch_id + cradle2,data=dat[dat$sp_id=="B",],weights=wts,family=binomial("logit"))
save(mod_noa_B,file=ppath(ddrive,"/iotc/iotc_mod2_B.RData")); rm(mod_noa_B); gc()
mod_noa_S <- glm(recap ~ tagger + ns(len5,df=4) + Qual + Cond2 + OTC + tag_type + tag_sch_id + cradle2,data=dat[dat$sp_id=="S",],weights=wts,family=binomial("logit"))
save(mod_noa_S,file=ppath(ddrive,"/iotc/iotc_mod2_S.RData")); rm(mod_noa_S); gc()
mod_noa_Y <- glm(recap ~ tagger + ns(len5,df=4) + Qual + Cond2 + OTC + tag_type + tag_sch_id + cradle2,data=dat[dat$sp_id=="Y",],weights=wts,family=binomial("logit"))
save(mod_noa_Y,file=ppath(ddrive,"/iotc/iotc_mod2_Y.RData")); rm(mod_noa_Y); gc()
load(ppath(ddrive,"/iotc/iotc_mod2_B.RData"))
load(ppath(ddrive,"/iotc/iotc_mod2_Y.RData"))
load(ppath(ddrive,"/iotc/iotc_mod2_S.RData"))
load(ppath(ddrive,"/iotc/iotc_mod2.RData"))
mod_noa_B$formula
mod_noa_Y$formula
mod_noa_S$formula
mod_att$formula
sum(AIC(mod_noa_S,mod_noa_B,mod_noa_Y))  # mod_noa has OTC but not cradle
AIC(mod_att)

# Load mod_att and the species-level mods with all tag types


# Plot figures for best model
load(ppath(ddrive,"/iotc/iotc_mod2.RData"))
output(dat=dat,                 res=mod_att,  pltagger=T,anon=T,pltagtype=T,pldiag=F,plxprt=F,plcradle=T,condor=F,pltagger.type=F,plOTC=T,plcond=T,plqual=T,plsize=T,plevent=T,do_ci=F,prefix_nm="iotc_att2",io_agg=T,pltype="pdf");graphics.off()
output(dat=dat[dat$sp_id=="B",],res=mod_noa_B,pltagger=T,anon=T,pltagtype=T,pldiag=F,plxprt=F,plcradle=T,condor=F,pltagger.type=F,plOTC=T,plcond=T,plqual=T,plsize=T,plevent=T,do_ci=T,prefix_nm="iotc_att2_B",io_agg=T,pltype="pdf");graphics.off()
output(dat=dat[dat$sp_id=="Y",],res=mod_noa_Y,pltagger=T,anon=T,pltagtype=T,pldiag=F,plxprt=F,plcradle=T,condor=F,pltagger.type=F,plOTC=T,plcond=T,plqual=T,plsize=T,plevent=T,do_ci=T,prefix_nm="iotc_att2_Y",io_agg=T,pltype="pdf");graphics.off()
output(dat=dat[dat$sp_id=="S",],res=mod_noa_S,pltagger=T,anon=T,pltagtype=T,pldiag=F,plxprt=F,plcradle=T,condor=F,pltagger.type=F,plOTC=T,plcond=T,plqual=T,plsize=T,plevent=T,do_ci=T,prefix_nm="iotc_att2_S",io_agg=T,pltype="pdf");graphics.off()


AIC_mod_att <- drop1(mod_att)
write.csv(AIC_mod_att,file="ioAIC_mod_att2.csv")

load(ppath(ddrive,"/iotc/iotc_mod2_B.RData"))
AIC_mod_noa_B <- drop1(mod_noa_B); rm(mod_noa_B)
save(AIC_mod_noa_B,file="ioAIC_mod2_B_noa.RData")  # Old version file="ioAIC_mod_noa.RData" lacked cradle2 and used Cond instead of cond2
load(ppath(ddrive,"/iotc/iotc_mod2_Y.RData"))
AIC_mod_noa_Y <- drop1(mod_noa_Y); rm(mod_noa_Y);gc()
save(AIC_mod_noa_Y,file="ioAIC_mod2_Y_noa.RData")  # Old version file="ioAIC_mod_noa.RData" lacked cradle2 and used Cond instead of cond2
load(ppath(ddrive,"/iotc/iotc_mod2_S.RData"))
AIC_mod_noa_S <- drop1(mod_noa_S); rm(mod_noa_S);gc()
save(AIC_mod_noa_S,file="ioAIC_mod2_S_noa.RData")  # Old version file="ioAIC_mod_noa.RData" lacked cradle2 and used Cond instead of cond2

write.csv(AIClist,file="AIClist.csv")


load(file=ppath(ddrive,"paper_xs/iotc/imod_noa.RData"))
AIC(mod_noa)
load(file="imod_noa_B.RData")
load(file="imod_noa_S.RData")
load(file="imod_noa_Y.RData")
sum(AIC(mod_noa_B,mod_noa_S,mod_noa_Y))

# Summarise for mod_att by species, with archival tags "recap ~ tagger.tagtype + ns(len5,df=4) + Qual + Cond + OTC + tag_sch_id"
load(file="ioAIC_mod2_B_noa.RData")
load(file="ioAIC_mod2_Y_noa.RData")
load(file="ioAIC_mod2_S_noa.RData")
#AIC_mod_NOAx <- read.csv(file="ioAIC_mod_noa.csv")
AIC_mod_NOAx <- read.csv(file="ioAIC_mod_att2.csv")
require(xtable)
rownms <- c("Full model","Tag type","Tagger","Tag placement quality","Tagging cradle","Fish condition","OTC","Length","Length:Species",                  "Tagging event")
lrows <- c("<none>",     "tag_type", "tagger","Qual",                "cradle2",       "Cond",          "OTC","ns(len5, df = 4)","ns(len5, df = 4):sp_id","tag_sch_id")
lrows2 <- c("<none>",     "tag_type", "tagger","Qual",                "cradle2",       "Cond2",          "OTC","ns(len5, df = 4)","ns(len5, df = 4):sp_id","tag_sch_id")
AICsummaryIO <- data.frame(mod=rownms,
            BET_YFT_SKJ= c(AIC_mod_NOAx[1,4],  AIC_mod_NOAx  [match(lrows2,AIC_mod_NOAx[,1]),4][-1] -         AIC_mod_NOAx[1,4]),
            BET=         c(AIC_mod_noa_B[1,3], AIC_mod_noa_B [match(lrows2,rownames( AIC_mod_noa_B)),3][-1] - AIC_mod_noa_B[1,3]),
            YFT=         c(AIC_mod_noa_Y[1,3], AIC_mod_noa_Y [match(lrows2,rownames( AIC_mod_noa_Y)),3][-1] - AIC_mod_noa_Y[1,3]),
            SKJ=         c(AIC_mod_noa_S[1,3], AIC_mod_noa_S [match(lrows2,rownames( AIC_mod_noa_S)),3][-1] - AIC_mod_noa_S[1,3]))
xtAICsummaryIO <- xtable(AICsummaryIO,digits=1)
print.xtable(xtAICsummaryIO,type="html",file="AICsummaryIO.html")
AICsummaryIO

# Table 3: Table of parameter estimates
# parameter estimates in the file
aBYS <- rbind(read.csv("iotc_att2_Qual_parests.csv"),read.csv("iotc_att2_cradle2_parests.csv"),read.csv("iotc_att2_OTC_parests.csv"))
aB <- rbind(read.csv("iotc_att2_B_Qual_parests.csv"),read.csv("iotc_att2_B_cradle2_parests.csv"),read.csv("iotc_att2_B_OTC_parests.csv"))
aS <- rbind(read.csv("iotc_att2_S_Qual_parests.csv"),read.csv("iotc_att2_S_cradle2_parests.csv"),read.csv("iotc_att2_S_OTC_parests.csv"))
aY <- rbind(read.csv("iotc_att2_Y_Qual_parests.csv"),read.csv("iotc_att2_Y_cradle2_parests.csv"),read.csv("iotc_att2_Y_OTC_parests.csv"))
cols <- c(6:8)
tab3 <- xtable(cbind(aBYS[,c(1,cols)],aB[,cols],aY[,cols],aS[,cols]),digits=2)
print.xtable(tab3,type="html",file="tab3.html")


AICsummaryALL <- data.frame(mod=AIC_mod_att[,1],deltaAIC_all=AIC_mod_att[,3]-AIC_mod_att[1,3])                           # Full model
AICsummarySPP <- data.frame(mod=rownames(AIC_mod_att_B_obj),deltaAIC_B=AIC_mod_att_B_obj[,3]-AIC_mod_att_B_obj[1,3],   #single sp models
deltaAIC_S=AIC_mod_att_S_obj[,3]-AIC_mod_att_S_obj[1,3],
deltaAIC_Y=AIC_mod_att_Y_obj[,3]-AIC_mod_att_Y_obj[1,3])
AICsummaryBY <- data.frame(mod=rownames(AIC_mod_att_BY_obj),deltaAIC_BY=AIC_mod_att_BY_obj[,3]-AIC_mod_att_BY_obj[1,3]) # BY model
AICsummary <- cbind(AICsummarySPP,deltaAIC_all=AICsummaryALL[c(1,2,1,3:8),2],deltaAIC_BY=AICsummaryBY[c(1,2,1,3:8),2])

# Table 2. estimate relative tag recovery
############### effective tags
# predict recoveries using existing model, then as if everyone was JPH and condition etc were good. Divide releases through by the ratio.
dir<- "C:/Users/simonh/Dropbox/Papers/Tagger_effects/iotc"
setwd(dir)
require(splines)
require(lme4)
require(pkpkg)
require(plotrix)
require(boot)
require(car)
require(MASS)
require(R4MFCL)
ddrive <- "D:/assessments/tagging/tagger_fx/TagFX_paper/"
dd_iotc <- ppath(ddrive,"iotc")
source("../R_scripts/load_iotc.r")
source("../R_scripts/prepare_data_both.r")
source("../R_scripts/output.r")
source("../R_scripts/support_functions.r")
infile <- "TAG_REC_AB3.txt"
iotags <- load_iotc(infile)
iotag_all <- prepare_data(iotags,maxprt=8000,xpf="TaggerExp.csv",reg="iotc")
dat <- iotagsub <- aggdata(dat=iotag_all,minTperSch=20,archival=T,xptype=F,tagger.tagtype=F,cond=F,cradle=T); str(iotagsub)
dat$cradle2[dat$cradle2=="ARC"] <- unique(dat$cradle2)[1]
dat$cradle2 <- as.factor(as.character(dat$cradle2))
levels(dat$cradle2)
##
load(ppath(ddrive,"/iotc/iotc_mod.RData"))
load(ppath(ddrive,"/iotc/iotc_mod2.RData"))
pfull <- predict.glm(mod_att,newdata=dat,type="response")
tag.pred <- dat
tag.pred$tagger <- dat$tagger[grep("JPH",dat$tagger)][1]
tag.pred$tag_type <- levels(dat$tag_type)[1]
tag.pred$cradle2 <- dat$cradle2[grep("Front",dat$cradle2)][1]
tag.pred$Cond <- levels(dat$Cond)[1]
tag.pred$Qual <- levels(dat$Qual)[1]
tag.pred$OTC <- levels(dat$OTC)[1]
ppred <- predict.glm(mod_att,newdata=tag.pred,type="response")
sum(dat$wts*pfull)/sum(dat$wts*ppred)
#for(ss in c("B","S","Y")) {
#  tag.ss <- iotagsub[iotagsub$sp_id == ss,]
#  pss <- cbind(tag.ss,pfull[iotagsub$sp_id == ss],ppred[iotagsub$sp_id == ss])
#  print(sum(pss$wts*pss[,14]) / sum(pss$wts*pss[,15]))
#  }
ntp=7
effectlist <- data.frame(tp=rep("x",ntp),ratio=rep(0,ntp),stringsAsFactors = F)
pfull <- predict.glm(mod_att,newdata=dat,type="response")
for(x in 1:ntp) {
  tag.pred <- dat
  if(x==1) { tag.pred$tagger <- dat$tagger[grep("JPH",dat$tagger)][1]; effectlist[x,]$tp <- "tagger" }
  if(x==2) { tag.pred$tag_type <- levels(dat$tag_type)[1]; effectlist[x,]$tp <- "tag_type" }
  if(x==3) { tag.pred$cradle2 <- dat$cradle2[grep("Front",dat$cradle2)][1]; effectlist[x,]$tp <- "cradle2" }
  if(x==4) { tag.pred$Cond <- levels(dat$Cond)[1]; effectlist[x,]$tp <- "Cond" }
  if(x==5) { tag.pred$Qual <- levels(dat$Qual)[1]; effectlist[x,]$tp <- "Qual" }
  if(x==6) { tag.pred$OTC <- levels(dat$OTC)[1]; effectlist[x,]$tp <- "OTC" }
  if(x==7) {
    tag.pred$tagger  <- dat$tagger[grep("JPH",dat$tagger)][1]
    tag.pred$OTC  <- levels(dat$OTC)[1]
    tag.pred$cradle2 <- dat$cradle2[grep("Front",dat$cradle2)][1]
    tag.pred$Cond <- levels(dat$Cond)[1]
    tag.pred$Qual <- levels(dat$Qual)[1]
    tag.pred$tag_type<- levels(dat$tag_type)[1]
    effectlist[x,]$tp <- "All"
    }
  ppred <- predict.glm(mod_att,newdata=tag.pred,type="response")
  tag.ss <- dat
  pss <- cbind(tag.ss,pfull,ppred)
  effectlist[x,]$ratio <- sum(pss$wts*pss$pfull) / sum(pss$wts*pss$ppred)
  }
effectlist
xteffectsIO <- xtable(effectlist,digits=3)
print.xtable(xteffectsIO,type="html",file="Tab2_IOeffects.html")

#### Table of values
a <- tapply(dat$wts,dat$tagger,sum);rbind(c("tagger",sum(a)),cbind(names(a),formatC(100*a/sum(a),digits=2,format="f")))
a <- tapply(dat$wts,dat$tag_type,sum);rbind(c("tag_type",sum(a)),cbind(names(a),formatC(100*a/sum(a),digits=2,format="f")))
a <- tapply(dat$wts,dat$cradle2,sum);rbind(c("cradle2",sum(a,na.rm=T)),cbind(names(a),formatC(100*a/sum(a,na.rm=T),digits=2,format="f")))
a <- tapply(dat$wts,dat$Cond,sum);rbind(c("Cond",sum(a,na.rm=T)),cbind(names(a),formatC(100*a/sum(a,na.rm=T),digits=2,format="f")))
a <- tapply(dat$wts,dat$Qual,sum);rbind(c("Qual",sum(a)),cbind(names(a),formatC(100*a/sum(a,na.rm=T),digits=2,format="f")))
a <- tapply(dat$wts,dat$OTC,sum);rbind(c("OTC",sum(a)),cbind(names(a),formatC(100*a/sum(a),digits=2,format="f")))

a <- tapply(dat$wts,dat$tagger,sum,na.rm=T)
names(a) <- seq(0,(length(a)-1),1)
barplot(a/sum(a),xlab="Tagger",ylab="Relative frequency", main="IO-RTTP")
savePlot("IO_Tagger_frequency_distribution",type="png")
savePlot("Figure_2_IO_Tagger_frequency_distribution",type="pdf")

