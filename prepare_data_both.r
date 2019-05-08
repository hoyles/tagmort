# prepare data

prepare_data <- function(tags,maxprt=8000,first1000=F,xpf=NA,reg="iotc") {

a <- table(tags$tagger)
a <- names(a[a>100])
tags <- tags[tags$tagger %in% a,]
 
# Expertise
a <- tags[order(tags$tagger,tags$tag_sch_id),]
taggers <- sort(unique(as.character(a$tagger)))
cumtags <- NULL
for (tg in taggers) {
  ct <- cumsum(a[a$tagger==tg,]$tagger==tg)
  ct <- c(0,ct[1:(length(ct)-1)])
  cumtags <- c(cumtags,ct)
  }
a.tg.sc <- paste(a$tagger,a$tag_sch_id)
min.tg.sc <- tapply(cumtags,a.tg.sc,min)   # get the minimum no. tags for the school
a$xprt <- min.tg.sc[match(a.tg.sc, names(min.tg.sc))]
a$xprt <- as.numeric(a$xprt)+1
#a[a$tagger %in% c("ADL","WJH","DGI","PGW","RNP","PBS"),]$xprt <- a[a$tagger %in% c("ADL","WJH","DGI","PGW","RNP","PBS"),]$xprt + floor(0.75 * maxprt)
a[a$xprt>maxprt,]$xprt <- maxprt
if(first1000==T){
  a[a$xprt<=200,]$xprt <- 0
  a[a$xprt>200 & a$xprt<=500,]$xprt <- 1
  a[a$xprt>500,]$xprt <- 2
  }
tags<-a
tags$tagger <- as.factor(as.character(tags$tagger))
if(reg=="iotc") tags$tagger <- relevel(tags$tagger, ref="JUD")
if(reg=="pttp") tags$tagger <- relevel(tags$tagger, ref="ADL")
if(reg=="iotc") {
  tags$assist <- NA
  crloc <- grep("cradle",names(tags))
  assloc <- grep("assist",names(tags))
  tags <- tags[,c(1:crloc,assloc,c((crloc+1):(assloc-1)))]
  }
if(reg=="pttp") tags$assist <- relevel(tags$assist, ref="CRW")

if(!is.na(xpf)) {
  xplist <- read.csv(xpf,stringsAsFactors=F,col.names=c("tg_code","tg","tuna","other","none"))
  xplist$xp[xplist$tuna=="x"] <- "T"
  xplist$xp[xplist$other=="x"] <- "O"
  xplist$xp[xplist$none=="x"] <- "N"
  tags$xptype <- xplist$xp[match(tags$tagger,xplist$tg_code)]
  } else tags$xptype <- NA

tags$xptype <- as.factor(tags$xptype)
if(reg=="iotc") tags$rel_date <- as.Date(tags$SIG_Date,format="%d/%m/%Y") # Release date

tty <- rep("conv",length(tags$tag_type))
tty[tags$tag_type %in% c("ET","SO","c.archival","d.sonic")] <- "internal"
tags$tagger.tagtype <- as.factor(paste(tags$tagger,tty,sep="."))
if(reg=="iotc") tags$tagger.tagtype <- relevel(tags$tagger.tagtype, ref="JUD.conv")
if(reg=="pttp") tags$tagger.tagtype <- relevel(tags$tagger.tagtype, ref="ADL.conv")

if(reg=="iotc") tags$tag_type <- relevel(tags$tag_type,ref="ST")
if(reg=="pttp") tags$tag_type <- relevel(tags$tag_type,ref="a.Y13")

if(reg=="iotc") {
  tmp <- tags$Cond
  tags[grepl("BLOOD",tags$TAG_Comments) & tmp == "Good" ,]$Cond <- "Bleeding"
  tags[grepl("lood",tags$TAG_Comments) & tmp == "Good" ,]$Cond <- "Bleeding"
  tags[grepl("leed",tags$TAG_Comments) & tmp == "Good" ,]$Cond <- "Bleeding"
}
tags$Cond2 <- as.character(tags$Cond)
tags$Cond2[tags$Cond2 %in% c("Dropped on deck","Hit side of boat")] <- "Impact"
tags$Cond2[tags$Cond2 %in% c("Bleeding","Mouth damage","Tail damage","Eye damage")] <- "Damaged during capture"
table(tags$Cond2)
if(reg=="iotc") tags$Cond2 <- factor(tags$Cond2,levels=c("Good","Impact","Damaged during capture","Shark bite","Too slow / other","Unknown"))
if(reg=="pttp") tags$Cond2 <- factor(tags$Cond2,levels=c("Damaged during capture","Good","Impact","Long time on hook","Minor palate damage","Shark bite","Unknown"))
tags$Cond2 <- relevel(tags$Cond2,ref="Good")
tags <- tags[!tags$Cond %in% c("Unknown","Minor palate damage"),]

if(reg=="iotc") {
  tmp <- tags$Qual
  tags[grepl(" low",tags$TAG_Comments) & tmp == "Good" & !grepl("lower",tags$TAG_Comments) & !grepl("vcr",tags$TAG_Comments),]$Qual <- "Badly placed"
  tags[grepl(" high",tags$TAG_Comments) & tmp == "Good" & !grepl("deck",tags$TAG_Comments),]$Qual <- "Badly placed"
  tags$Qual <- relevel(tags$Qual,ref="Good")
  }
tags <- tags[!tags$Qual %in% c("Tag lost","Tag rejected","Unknown"),]

table(tags$cradle)
if(reg=="iotc") {
  tags$cradle2 <- NA
  cradles <- cbind(c("ARC","BPC","BSC","FPC","FSC","MAT","MPC","MSC","SPC","SSC"),c("ARC","Bow","Bow","Front","Front","X","Mid","Mid","Stern","Stern"))
  tags$cradle2 <- as.factor(cradles[match(tags$cradle,cradles[,1]),2])
  }
if(reg=="pttp") {
  tags$cradle2 <- as.character(tags$cradle)
  tags[tags$cradle %in% c("PB"),]$cradle2 <- "BOW"
  tags[tags$cradle %in% c("SB","AB"),]$cradle2 <- "SABOW"
  tags[tags$cradle %in% c("PM","SM","MS","BM"),]$cradle2 <- "MID"
  tags[tags$cradle %in% c("PS","SS","SR","AS","TS"),]$cradle2 <- "STN"
  tags$cradle2 <- relevel(as.factor(tags$cradle2),ref="BOW")
  }

tags$assoc_id <- as.factor(as.character(tags$assoc_id))
if(reg=="iotc") tags$len5 <- 5*floor(tags$length/5)
tags <- tags[tags$len5>15,]

if(reg=="iotc") {
  tags$recap <- !(is.na(tags$REC_Tag1) & is.na(tags$REC_Tag2))
  tags$recap_tag1 <- !is.na(tags$REC_Tag1)
  tags$TAG_Comments <- as.character(tags$TAG_Comments)
  }
if(reg=="pttp") {
  tags$recap <- (tags$tag_recovered=="Y" | tags$dbl_tag_recovered=="Y") & !is.na(tags$tag_recovered)
  tags$recap_tag1 <- tags$tag_recovered=="Y" & !is.na(tags$tag_recovered)
  }

tags$tag_sch_id <- as.factor(tags$tag_sch_id)
tags$TAG_SchoolNb <- NULL

if(reg=="iotc") names(tags) <- c("proj_id","vessel","tag_no","sp_id","length","Qual","Cond","tag_type","tagger","cradle","assist","assoc_id","SIG_SchoolType",
    "SIG_Date","SIG_NS","SIG_LatDegree","SIG_LatMinute","SIG_LonDegree","SIG_LonMinute","REC_Tag1","REC_Tag2","TAG_Comments","tag_sch_id","OTC","rel_lat","rel_lon",
    "xprt","xptype","rel_date","tagger.tagtype","Cond2","cradle2","len5","recap","recap_tag1")
if(reg=="pttp") names(tags) <- c("proj_id", "vessel", "tag_sch_id", "tag_no", "tagger", "tag_type", "assoc_id", "Qual", "Cond", "cradle", "assist", "len5",
    "sp_id", "tag_rel_id", "Arc_tag_no", "Sonic_tag_no", "OTC", "Dbl_tag_no", "recovery", "rel_date", "best_catch_date", "dbl_tag", "rel_lat",
    "rel_lon", "tag_recovered", "dbl_tag_recovered", "xprt", "xptype", "tagger.tagtype", "Cond2", "cradle2", "recap", "recap_tag1")

tags <- tags[,c("proj_id","vessel","tag_sch_id","tag_no","tagger","tag_type","assoc_id","Qual","Cond","cradle","assist","sp_id","len5",    "OTC","rel_lat","rel_lon","xprt","xptype","rel_date","tagger.tagtype","Cond2","cradle2","recap","recap_tag1")]

tags$rel_lat5long <- paste(5*floor(tags$rel_lat/5),5*floor(tags$rel_lon/5),sep="_")
tags$rel_lat2long <- paste(2*floor(tags$rel_lat/2),2*floor(tags$rel_lon/2),sep="_")
tags$rel_lat1long <- paste(1*floor(tags$rel_lat/1),1*floor(tags$rel_lon/1),sep="_")


return(tags)
}
