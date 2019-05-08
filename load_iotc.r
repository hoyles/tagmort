#load_iotc <- function(infile="C:/Users/simonh/Dropbox/Papers/Tagger_effects/TAG_REC_AB_FLAG2.txt") {
load_iotc <- function(infile="C:/Users/simonh/Dropbox/Papers/Tagger_effects/TAG_REC_AB3.txt") {
dat <- read.csv(infile)  #

dat <- dat[dat$TAG_Project!="SS",]    # small-scale project
dat$TAG_Vessel <- as.factor(as.character(dat$TAG_Vessel))

dat <- dat[,c("TAG_Project","TAG_Vessel","TAG_SchoolNb","TAG_Tag1","TAG_Sp","TAG_Length","TAG_Tag1Rel",
"TAG_FishRel","TAG_Type","TAG_Tagger","TAG_Cradle","SIG_SchoolAssociation","SIG_SchoolType","SIG_Date","SIG_NS","SIG_LatDegree","SIG_LatMinute","SIG_LonDegree","SIG_LonMinute",
"REC_Tag1","REC_Tag2","TAG_Comments")]

names(dat) <- c("proj_id","vessel","TAG_SchoolNb","tag_no","sp_id","length","Qual","Cond","tag_type","tagger","cradle","assoc_id","SIG_SchoolType",
    "SIG_Date","SIG_NS","SIG_LatDegree","SIG_LatMinute","SIG_LonDegree","SIG_LonMinute","REC_Tag1","REC_Tag2","TAG_Comments")

dat$tag_sch_id <- paste(dat$TAG_SchoolNb,dat$vessel)

dat$tag_type[dat$tag_type=="St"] <- "ST"
dat$tag_type <- as.character(dat$tag_type)
dat$OTC <- F
dat$OTC[dat$tag_type %in% c("OT","OTS")] <- T
dat$tag_type[dat$tag_type=="OTS"] <- "SO"
dat$tag_type[dat$tag_type=="OT"] <- "ST"
dat$tag_type[dat$tag_type==""] <- "UN"
dat$tag_type <- as.factor(dat$tag_type)

dat <- dat[dat$Qual %in% c(1,2,7),]  # Tag Reliability for the 1st tag: 1 = Good; 2 = Badly placed; 3 = rejected; 4 = lost; 5 = seeded; 6 = doubts on the tag number; 7 = Unknown.
dat$Qual <- as.factor(as.character(dat$Qual))
levels(dat$Qual) <- c("Good","Badly placed","Unknown")
dat$Qual <- relevel(dat$Qual,ref="Good")

dat$Cond <- as.factor(as.character(dat$Cond))
levels(dat$Cond) <- c("Good","Bleeding","Tail damage","Mouth damage","Dropped on deck","Hit side of boat","Shark bite","Too slow / other","Unknown")
dat$Cond <- relevel(dat$Cond,ref="Good")

dat$cradle <- as.factor(as.character(dat$cradle))

conv_lat <- function(NS,degree,minute) {
  outlat <- degree + minute / 60
  outlat <- ifelse(NS == "S", -1 * outlat, outlat)
  return(outlat)
}
conv_lon <- function(EW,degree,minute) {
  outlon <- degree + (minute / 60)
  outlon <- ifelse(EW == "W", 360 - outlon, outlon)
  return(outlon)
  }
conv_lon <- function(degree,minute) {
  outlon <- degree + (minute / 60)
  return(outlon)
  }
conv_locs <- function(a) {
  a$rel_lat <- conv_lat(a$SIG_NS,a$SIG_LatDegree,a$SIG_LatMinute)
  a$rel_lon <- conv_lon(a$SIG_LonDegree,a$SIG_LonMinute)
  return(a)
  }
dat <- conv_locs(dat)


return(dat)
}
