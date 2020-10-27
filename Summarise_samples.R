#Script for summarising available samples


# Data section ------------------------------------------------------------

#Sharks data base 
User="Matias"
source("C:/Matias/Analyses/SOURCE_SCRIPTS/Git_other/Source_Shark_bio.R")


#Biological storage room
library(tidyverse)
library(dplyr)
library(readxl)
Biol.storage.room=read_excel("C:/Matias/Data/Shark_bio/Biol store samples.xlsx")

#Vertebrae in freezer
Vertebrae.frozen=read.csv('U:/Shark/Fish_processor_age_sampling.csv',stringsAsFactors = F)


#Dried vertebrae
Vertebrae.dried_boat=read_excel("C:/Matias/Data/Shark_bio/Santi_historic.vertebrae_General sample data.xlsx",
                           sheet='Boat')
Vertebrae.dried_new=read_excel("C:/Matias/Data/Shark_bio/Santi_historic.vertebrae_General sample data.xlsx",
                                sheet='NEW')



# Genetic samples ------------------------------------------------------------
DATA=DATA%>%mutate(BAG_NO=tolower(BAG_NO),
                   BAG_NO.match=paste(BAG_NO,SPECIES),
                   BAG_NO=ifelse(BAG_NO=='',NA,BAG_NO),
                   FinClipFlag=tolower(FinClipFlag),
                   MuscleFlag=tolower(MuscleFlag))

#Genetic samples in biological store room
Genetic.biol.store.room=Biol.storage.room%>%
                          data.frame%>%
                          mutate(DATA=tolower(DATA))%>%
                          filter(DATA=='yes')%>%
                          rename(BAG_NO=SAMPLE..)%>%
                          mutate(BAG_NO=tolower(BAG_NO),
                                 BAG_NO.match=paste(BAG_NO,SSP))
Genetic.biol.store.room=Genetic.biol.store.room%>%
                          left_join(DATA%>%
                                      filter(BAG_NO%in%unique(Genetic.biol.store.room$BAG_NO))%>%
                                      dplyr::select(BAG_NO.match,SPECIES,COMMON_NAME,SCIENTIFIC_NAME,SHEET_NO,
                                                    date,Mid.Lat,Mid.Long,FL,TL,SEX),by='BAG_NO.match')%>%
                          mutate(Sp.match=ifelse(SSP==SPECIES,"yes","no"))
Bag.number.not.matching=Genetic.biol.store.room%>%
                            filter(Sp.match!='yes')
Genetic.biol.store.room=Genetic.biol.store.room%>%
                            filter(Sp.match=='yes')%>%
                            mutate(Data.set="Biol.store.room")

#Genetics from PA project (only need to accout for tissue from vertebare, muscle and fin clips already in Biol storage)
Genetic.PA.samples=DATA[grep("PA", DATA$SHEET_NO), ]%>%
                      mutate(VERT_SAMPL=tolower(VERT_SAMPL))%>%
                      filter(year>=2020 & VERT_SAMPL=='yes')%>%
                      mutate(Data.set="ParksAustralia.vertebrae")

#Genetics from tissue from frozen vertebrae
SP.names=DATA%>%distinct(SPECIES,.keep_all = T)%>%
  dplyr::select(SPECIES,COMMON_NAME,SCIENTIFIC_NAME)
Genetic.Vertebrae.frozen=Vertebrae.frozen %>%
            mutate(Species=substr(Specimen.No.,1,2),
                   Date=as.POSIXct(Sampling.date,format="%d/%m/%Y"),
                   Year=year(Date),
                   Month=month(Date),
                   Sex=ifelse(Sex=='',NA,Sex),
                   Area=ifelse(grepl("Fatal Attraction",Comments),"Esperance","Albany"),
                   Gear=ifelse(grepl("LONGLINE",Comments),"Longline","Net"))%>%
  filter(!Species=="" & !Recorder=="")%>%
  left_join(SP.names,by=c('Species'='SPECIES'))

  #conversions from PA analysis (Missing: issues with GM, dummy for TK)
IDL.to.TL_FL=data.frame(Species=c('BW','GM','WH','TK'),
                     slope.TL=c(2.6203,1.3001,2.0467,2.6),
                     inter.TL=c(9.8941,49.586,19.151,9),
                     slope.FL=c(2.0313,1.1966,1.8562,2.1),
                     inter.FL=c(9.3876,39.876,13.5,9))

Genetic.Vertebrae.frozen=Genetic.Vertebrae.frozen%>%
                    left_join(IDL.to.TL_FL,by='Species')%>%
                    mutate(TL=ifelse(is.na(TL),IDL*slope.TL+inter.TL,TL),
                           FL=ifelse(is.na(FL),IDL*slope.FL+inter.FL,FL),
                           Mid.Lat=ifelse(Area=='Albany',-35,ifelse(Area=='Esperance',-33.8,NA)),
                           Mid.Long=ifelse(Area=='Albany',117.9,ifelse(Area=='Esperance',121.9,NA)),
                           Data.set="Frozen.vertebrae")%>%
                    rename(date=Date,
                           SEX=Sex,
                           BAG_NO=Specimen.No.)



this=c('Data.set','date','COMMON_NAME','SCIENTIFIC_NAME','FL','Mid.Lat','Mid.Long','SEX','BAG_NO')
Dat=rbind(Genetic.biol.store.room[,this],
          Genetic.PA.samples[,this],
          Genetic.Vertebrae.frozen[,this])
  
table(Dat$COMMON_NAME)

tab=Dat%>%
    filter(COMMON_NAME%in%c('Whiskery shark','Dusky shark','Gummy Shark','Sandbar shark'))%>%
    group_by(COMMON_NAME)%>%
    summarise(min.size=min(FL,na.rm=T),
         max.size=max(FL,na.rm=T),
         mean.size=mean(FL,na.rm=T),
         sd.size=sd(FL,na.rm=T))%>%
  data.frame()

Dat%>%
  filter(!is.na(SEX) & COMMON_NAME%in%c('Whiskery shark','Dusky shark','Gummy Shark','Sandbar shark'))%>%
  ggplot(aes(x=FL, fill=COMMON_NAME)) +
  geom_density(alpha=0.4) +
  facet_grid(. ~ SEX)
ggsave("C:/Matias/Analyses/Samples/gen_size.dist.tiff", width = 12,height = 8, dpi = 300,compression = "lzw")



#Table of genetic samples by data set and species
T1=Dat%>%
  group_by(Data.set,COMMON_NAME)%>%
  tally()%>%
  spread(Data.set,n,fill=0)%>%
  data.frame
write.csv(T1,'C:/Matias/Analyses/Samples/Gen.samples.in.storage.csv',row.names=F)


# Vertebrae samples ------------------------------------------------------------
#Historic dried vertebrae
Vertebrae.dried_boat=Vertebrae.dried_boat%>%
              rename(BAG_NO="BAG NO",
                     Box="Box Location")%>%
              mutate(Available=ifelse(Available=="???",NA,Available),
                     Available=tolower(Available),
                     SEX=tolower(SEX),
                     BAG_NO=tolower(BAG_NO),
                     Data.set='historic_boat')%>%
              filter(Available=='yes')

Vertebrae.dried_new=Vertebrae.dried_new%>%
              rename(BAG_NO="BAG NO",
                     Box="Box Location")%>%
              mutate(Available=ifelse(Available=="???",NA,Available),
                     Available=tolower(Available),
                     SEX=tolower(SEX),
                     BAG_NO=tolower(BAG_NO),
                     Data.set='historic_new')%>%
              filter(Available=='yes')

this.var=c('SHEET_NO','LINE_NO','UNIQUE_ID','SPECIES','TL','FL','SEX','BAG_NO','Box')
Vertebrae.dried=rbind(Vertebrae.dried_boat%>%dplyr::select(all_of(this.var)),
                      Vertebrae.dried_new%>%dplyr::select(all_of(this.var)))%>%
                  mutate(BAG_NO.match=paste(BAG_NO,SPECIES))%>%
                  left_join(SP.names,by='SPECIES')
Vertebrae.dried=Vertebrae.dried%>%
                  left_join(DATA%>%
                              filter(BAG_NO%in%unique(Vertebrae.dried$BAG_NO) & !is.na(BAG_NO))%>%
                              dplyr::select(BAG_NO.match,date,year,Mid.Lat,Mid.Long,zone),by='BAG_NO.match')
  

Vertebrae.dried%>%
  filter(!is.na(SEX) & COMMON_NAME%in%c('Whiskery shark','Dusky shark','Gummy Shark','Sandbar shark'))%>%
  filter(SEX!='?')%>%
  ggplot(aes(x=FL, fill=COMMON_NAME)) +
  geom_density(alpha=0.4) +
  facet_grid(. ~ SEX)
ggsave("C:/Matias/Analyses/Samples/vertebrae_size.dist.tiff", width = 12,height = 8, dpi = 300,compression = "lzw")

sort(table(Vertebrae.dried$COMMON_NAME))

Tab.age.obs_yr=Vertebrae.dried%>%
  filter(COMMON_NAME%in%c('Whiskery shark','Dusky shark','Gummy Shark','Sandbar shark'))%>%
  group_by(COMMON_NAME,year)%>%
  tally%>%
  spread(COMMON_NAME,n,fill=0)

Tab.age.obs_yr.zone=Vertebrae.dried%>%
  filter(COMMON_NAME%in%c('Whiskery shark','Dusky shark','Gummy Shark','Sandbar shark'))%>%
  group_by(COMMON_NAME,year,zone)%>%
  tally%>%
  spread(COMMON_NAME,n,fill=0)


#these are not in Boat bio but have a location and date
Missing=Vertebrae.dried_new%>%
  filter(SPECIES%in%c("WH","GM","TK","BW"))%>%
  dplyr::select(Box,SHEET_NO,LINE_NO,UNIQUE_ID,ID2,ID3,
                SPECIES,TL,FL,SEX,BAG_NO,Vessel,Location,Date)