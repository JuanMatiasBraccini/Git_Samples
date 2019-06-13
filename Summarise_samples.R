a=read.csv("H:/Matias WA Fisheries/Analyses/Samples/Book2.csv")

b=read.csv("H:/Matias WA Fisheries/Data/Species.code.csv")
library(RODBC)
setwd("M:/Fisheries Research/Production Databases/Shark")  # working directory
channel <- odbcConnectAccess2007("Sharks.mdb")      
Boat_bio=sqlFetch(channel, "Boat_bio", colnames = F) 
close(channel)


a$sample.type=as.character(a$sample.type)
a$sample.type=with(a,ifelse(
  sample.type=="blood serum?","blood serum",
  ifelse(sample.type%in%c("Fin","fins/vert","vert/ fin"),"fin clip",
         ifelse(sample.type%in%c("Jaw/vert"),"jaw",
      ifelse(sample.type=="MUSCLE","muscle",
    ifelse(sample.type%in%c("tissues","?","","Heart, liver, vertebrae, muscle",
                            "Liver, spiral valve, muscle, vertebrae, reproductive tissue"),"other tissues",
      ifelse(sample.type%in%c("vertebrae","vert/tissue"),"Vertebrae",
             sample.type)))))))
a=subset(a,!(sample.type=="na"))

a$ssp=as.character(a$ssp)

a$Bag=as.character(a$"Bag..")
Boat_bio$BAG_NO=as.character(Boat_bio$"BAG NO")
Boat_bio=Boat_bio[Boat_bio$BAG_NO%in%unique(a$"Bag"),]
Boat_bio$SPECIES2=as.character(Boat_bio$SPECIES)
Boat_bio=subset(Boat_bio,select=c(SPECIES2,BAG_NO))
a=merge(a,Boat_bio,by.x="Bag",by.y="BAG_NO",all.x=T)


a$ssp=with(a,ifelse(is.na(ssp)|ssp%in%c("?","","TK?","PE?","na","LP (ZE)","BW?"),SPECIES2,ssp))



b=subset(b,Species%in%unique(a$ssp))
a=merge(a,b,by.x="ssp",by.y="Species",all.x=T)

a=subset(a,!is.na(ssp))
a=subset(a,!ssp=="XX")


A=table(a$ssp,a$sample.type,useNA = 'ifany')
B=data.frame(SPECIES=row.names(A))
B=merge(B,subset(b,select=c(Species,COMMON_NAME)),
        by.x="SPECIES",by.y="Species",all.x=T)
B$COMMON_NAME=as.character(B$COMMON_NAME)
B$COMMON_NAME=with(B,ifelse(is.na(COMMON_NAME),as.character(SPECIES),COMMON_NAME))
D=as.matrix(A) 
row.names(D)=B$COMMON_NAME
write.csv(D,"H:/Matias WA Fisheries/Analyses/Samples/sumary.csv")
