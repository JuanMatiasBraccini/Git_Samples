#Script for summarising available samples

#note: vertebrae_processor.R summarises vertebrae collected from southern seafood producers in 2019

library(tidyverse)
library(dplyr)
library(readxl)

# Data section ------------------------------------------------------------
if(!exists('handl_OneDrive')) source('C:/Users/myb/OneDrive - Department of Primary Industries and Regional Development/Matias/Analyses/SOURCE_SCRIPTS/Git_other/handl_OneDrive.R')

#1. Sharks data base 
User="Matias"
source(handl_OneDrive("Analyses/SOURCE_SCRIPTS/Git_other/Source_Shark_bio.R"))


#2. Biological storage room
Biol.storage.room=read_excel(handl_OneDrive("Data/Shark_bio/Biol store samples.xlsx"))
Stocktake_Abbey=read_excel(handl_OneDrive('Data/Shark_bio/Shakr genetics Stocktake.xlsx'),sheet='Good Samples',skip=1) #done in Nov 2022
TG.in.storage=read.csv(handl_OneDrive('Data/Shark_bio/TG tissue sample in stock_2024_02.csv'))

#3. Vertebrae in freezer from processor sampling (Southern Seafood Producers)
Vertebrae.frozen=read.csv('M:/Production Databases/Shark/Fish_processor_age_sampling.csv',stringsAsFactors = F)

  #modify species for sandbar samples genetically identified as whiskery shark by Brenton Pember
Vertebrae.frozen=Vertebrae.frozen%>%
  mutate(Species=substr(Specimen.No.,1,2),
         Species=ifelse(Specimen.No.%in%c("TK190315-28","TK190409-8","TK190429-06"),"WH",
                        Species))%>%
  filter(!Specimen.No.=="")

  #stocktake of freezer vertebrae done by Matt C. May 2024 (this has PA_2019 and Albany Fish Processor samples)
NMs=c('Sheet1','PA incom.','GM incom.','BW incom.','TK incom.','WH incom.','N incom.','R incom.','S incom.')
Stocktake_Matt=vector('list',length(NMs))
names(Stocktake_Matt)=NMs
for(s in 1:length(Stocktake_Matt))
{
  NN=names(Stocktake_Matt)[s]
  SKIP=0
  if(NN=='Sheet1') SKIP=1
  Stocktake_Matt[[s]]=read_excel(handl_OneDrive('Data/Shark_bio/Vertebrae tracking_1_Mat C_2024_05.xlsx'),sheet=NN,skip=SKIP)
}
Stocktake_Matt$Sheet1=Stocktake_Matt$Sheet1%>%
  rename(Proccesed.by.Genetics="Processed by...2",
         Date.genetics="Date...3",
         Proccesed.by.cut.up="Processed by...4",
         Date.cut.up="Date...5")
for(s in 2:length(Stocktake_Matt)) #note: 'Dried': vertebrae cleaned; 'Genetics': muscle collected
{
  Stocktake_Matt[[s]]=Stocktake_Matt[[s]]%>%
                        rename(ID='Frozen Samples',
                               No.record="No record")
}

ID.samples.not.in.Sheet1=unique(do.call(rbind,Stocktake_Matt[-1])%>%pull(ID))
ID.samples.not.in.Sheet1=ID.samples.not.in.Sheet1[which(!ID.samples.not.in.Sheet1 %in% unique(Stocktake_Matt$Sheet1$ID))]
dumi=Stocktake_Matt$Sheet1[1:length(ID.samples.not.in.Sheet1),]
dumi[,]=NA
Available.freezer.samples=rbind(Stocktake_Matt$Sheet1,
                                dumi%>%mutate(ID=ID.samples.not.in.Sheet1))%>%
                          distinct(ID,.keep_all = TRUE)


#4. Dried vertebrae
  #Santiago Bertacca revised by Matt Cieslak
#note from Matt: The complete datasheet for the stored vertebrae samples. 
#                 Since I added extra information such as condition and secondary box location there are 
#                 some samples without this. You will now be able to determine what is in the storage as 
#                 all samples have been processed with YES if available and blank if not. 
#                 If we decide to use these samples I would suggest that they are all reorganized with 
#                 appropriate secondary boxes, labels and a more succinct directory than using this database

Vertebrae.dried_boat=read_excel(handl_OneDrive("Data/Shark_bio/Santi_historic.vertebrae_revised by Mat.xlsx"),
                           sheet='Boat')
Vertebrae.dried_new=read_excel(handl_OneDrive("Data/Shark_bio/Santi_historic.vertebrae_revised by Mat.xlsx"),
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
TG.in.storage=TG.in.storage%>%
            mutate(Sample=tolower(Sample))%>%
            filter(!Sample%in%c("","bag # i?"))
add.dis.tiger=TG.in.storage$Sample[which(!TG.in.storage$Sample%in%Genetic.biol.store.room$BAG_NO)]
add.dis.tiger=data.frame(BAG_NO=add.dis.tiger,
                         SSP='TG')%>%
              mutate(DATA=NA,Notes=NA,BAG_NO.match=paste(BAG_NO,SSP))

Genetic.biol.store.room=rbind(Genetic.biol.store.room,add.dis.tiger)

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

#Genetics from PA project (only need to account for tissue from vertebrae, muscle and fin clips already in Biol storage)
Genetic.PA.samples=DATA[grep("PA", DATA$SHEET_NO), ]%>%
                      mutate(VERT_SAMPL=tolower(VERT_SAMPL))%>%
                      filter(year>=2020 & VERT_SAMPL=='yes')%>%
                      mutate(Data.set="ParksAustralia.vertebrae")

#Genetics from tissue from frozen vertebrae
SP.names=DATA%>%distinct(SPECIES,.keep_all = T)%>%
  dplyr::select(SPECIES,COMMON_NAME,SCIENTIFIC_NAME)
Genetic.Vertebrae.frozen=Vertebrae.frozen %>%
            mutate(Date=as.POSIXct(Sampling.date,format="%d/%m/%Y"),
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
ggsave(handl_OneDrive("Analyses/Samples/gen_size.dist_indicator species.tiff"), width = 12,height = 8, dpi = 300,compression = "lzw")



#Table of genetic samples by data set and species
T1=Dat%>%
  group_by(Data.set,COMMON_NAME)%>%
  tally()%>%
  spread(Data.set,n,fill=0)%>%
  data.frame
write.csv(T1,handl_OneDrive('Analyses/Samples/Gen.samples.in.storage.csv'),row.names=F)



# Stock take done by Matt C for Kyle Zawada-------------------------------------------------------------------------
Out.Kyle=rbind(
  DATA%>%
    mutate(BAG_NO=tolower(BAG_NO))%>%
    filter(BAG_NO%in%tolower(Available.freezer.samples$ID))%>%
    dplyr::select(SHEET_NO,BAG_NO,date,Mid.Lat,Mid.Long,COMMON_NAME,SEX,FL),
  Dat%>%
    distinct(BAG_NO,.keep_all = T)%>%
    left_join(DATA%>%dplyr::select(SHEET_NO,BAG_NO),by='BAG_NO')%>%
    mutate(BAG_NO=tolower(BAG_NO),
           BAG_NO=ifelse(grepl('-',BAG_NO) & nchar(BAG_NO)==10,paste0(str_extract(BAG_NO, "[^-]+"),'-',0,str_extract(BAG_NO, "(?<=-)[[:digit:]]+")),
                         BAG_NO))%>%
    filter(BAG_NO%in%tolower(Available.freezer.samples$ID))%>%
    dplyr::select(SHEET_NO,BAG_NO,date,Mid.Lat,Mid.Long,COMMON_NAME,SEX,FL))%>%
  arrange(COMMON_NAME)%>%
  distinct(BAG_NO,.keep_all=T)
Out.Kyle=Out.Kyle%>%
  left_join(Available.freezer.samples%>%
              mutate(BAG_NO=tolower(ID),
                     Sample.status=ifelse(!is.na(Date.genetics) | !is.na(Date.cut.up),'dried vertebrae and muscle in ethanol',
                                          'frozen vertebrae with muscle'))%>%
              dplyr::select(BAG_NO,Sample.status),
            by='BAG_NO')
write.csv(Out.Kyle,handl_OneDrive("Analyses/Samples/Out.Kyle.Zawada_freezer.samples.csv"),row.names = F)

in.freezer.but.no.metadata=Available.freezer.samples[which(!tolower(Available.freezer.samples$ID)%in%Out.Kyle$BAG_NO),]

# Vertebrae in freezer (including processor samples)------------------------------------------------------------
Out.Kyle%>%
  ggplot(aes(FL))+
  geom_histogram()+
  facet_wrap(~COMMON_NAME, scales = 'free')
ggsave(handl_OneDrive("Analyses/Samples/vertebrae_freezer(all, including processor)_size.dist.tiff"), width = 12,height = 8, dpi = 300,compression = "lzw")


write.csv(Out.Kyle%>%
            mutate(Year=year(date))%>%
            group_by(COMMON_NAME,Year)%>%
            tally()%>%
            spread(Year,n,fill=0),
          handl_OneDrive("Analyses/Samples/vertebrae_freezer(all, including processor)_samples by species and year.csv"),row.names = F)


# Historic dried vertebrae ------------------------------------------------------------
Vertebrae.dried_boat=Vertebrae.dried_boat%>%
              rename(BAG_NO="BAG NO",
                     Box="Box no.")%>%
              mutate(Available=ifelse(Available=="???",NA,Available),
                     Available=tolower(Available),
                     SEX=tolower(SEX),
                     BAG_NO=tolower(BAG_NO),
                     Data.set='historic_boat')%>%
              filter(Available=='yes')

Vertebrae.dried_new=Vertebrae.dried_new%>%
              rename(BAG_NO="BAG NO",
                     Box="Box no.")%>%
              mutate(Available='yes',
                     #Available=ifelse(Available=="???",NA,Available),
                     #Available=tolower(Available),
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
                              dplyr::select(BAG_NO.match,date,year,Mid.Lat,Mid.Long,zone),by='BAG_NO.match')%>%
  filter(!SEX%in%c("j143",'?'))
  

Vertebrae.dried%>%
  filter(!is.na(SEX) & COMMON_NAME%in%c('Whiskery shark','Dusky shark','Gummy Shark','Sandbar shark'))%>%
  ggplot(aes(x=FL, fill=COMMON_NAME)) +
  geom_density(alpha=0.4) +
  facet_grid(. ~ SEX)
ggsave(handl_OneDrive("Analyses/Samples/vertebrae_historic_size.dist.tiff"), width = 12,height = 8, dpi = 300,compression = "lzw")


write.csv(sort(table(Vertebrae.dried$COMMON_NAME)),
          handl_OneDrive("Analyses/Samples/vertebrae_historic_samples by species.csv"),row.names = F)

Tab.age.obs_yr=Vertebrae.dried%>%
  filter(COMMON_NAME%in%c('Whiskery shark','Dusky shark','Gummy shark','Sandbar shark'))%>%
  group_by(COMMON_NAME,year)%>%
  tally%>%
  spread(COMMON_NAME,n,fill=0)
write.csv(Tab.age.obs_yr,
          handl_OneDrive("Analyses/Samples/vertebrae_historic_samples by year_indicators.csv"),row.names = F)


Tab.age.obs_yr.zone=Vertebrae.dried%>%
  filter(COMMON_NAME%in%c('Whiskery shark','Dusky shark','Gummy shark','Sandbar shark'))%>%
  group_by(COMMON_NAME,year,zone)%>%
  tally%>%
  spread(COMMON_NAME,n,fill=0)
write.csv(Tab.age.obs_yr.zone,
          handl_OneDrive("Analyses/Samples/vertebrae_historic_samples by year_zone_indicators.csv"),row.names = F)


#these are not in Boat bio but have a location and date
Missing=Vertebrae.dried_new%>%
  filter(SPECIES%in%c("WH","GM","TK","BW"))%>%
  dplyr::select(Box,SHEET_NO,LINE_NO,UNIQUE_ID,ID2,ID3,
                SPECIES,TL,FL,SEX,BAG_NO,Vessel,Location,Date)



# Export all bag numbers-------------------------------------------------------------------------

Out.bag.number=DATA%>%
  filter(!is.na(BAG_NO))%>%
  dplyr::select(BAG_NO,SHEET_NO,LINE_NO,date,SPECIES,COMMON_NAME,FinClipFlag,MuscleFlag,VERT_SAMPL)%>%
  arrange(BAG_NO)

Dat.bagnumbers=sort(unique(Dat$BAG_NO))
A=Dat.bagnumbers[which(!Dat.bagnumbers%in%unique(Out.bag.number$BAG_NO))]
Out.bag.number.not.in.DATA=Dat%>%
  filter(BAG_NO%in%A)

write.csv(Out.bag.number,handl_OneDrive("Analyses/Samples/Out.bag.number.csv"),row.names = F)
write.csv(Out.bag.number.not.in.DATA,handl_OneDrive("Analyses/Samples/Out.bag.number.not.in.DATA.csv"),row.names = F)


