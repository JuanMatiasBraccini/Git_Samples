#Script for tallying collection of vertebrae from fish processors

#best way of identifying logbook trip:
# Identify trips from combining logbook data and VMS
#Nick: "Re the trip info, are you meaning the logbook returns which we could get from Anja (fishing block, weights etc), or VMS, where we would put in a data request form for each trip. That would give us gps of where they went fishing.
# Sampling date from the datasheet would be a day or 2 ahead of the trip finish date"
#Tim: "logbook return daily sheets (with block numbers) will be ahead of unload dates , by many days in some instances. Available from Anja.
#   VMS data won't distinguish which fish caught when, but in combination with logbook might narrow location down if different spp. caught. Will show discrete areas fished and allow us to remove areas where just transiting. Available by request (or from friendly FMO with a few minutes to spare?"



library(dplyr)
library(lubridate)
if(!exists('handl_OneDrive')) source('C:/Users/myb/OneDrive - Department of Primary Industries and Regional Development/Matias/Analyses/SOURCE_SCRIPTS/Git_other/handl_OneDrive.R')

# Read in data-------------------------------------------------------------------------
dat=read.csv('M:/Production Databases/Shark/Fish_processor_age_sampling.csv',stringsAsFactors = F)
#dat=read.csv(handl_OneDrive('/Data/Shark_bio/vertebrae collection_processor/Fish_processor_age_sampling.csv'),stringsAsFactors = F)  #offline

# Input parameters-------------------------------------------------------------------------
IDL.to.FL.conversion=data.frame(Species=c("BW","GM","WH","TK"),
                                intercept=c(10.604,16.369,20.284,14.999),
                                slope=c(1.9881,1.7181,1.7207,1.6541))
# Manipulate data-------------------------------------------------------------------------
dat=dat %>% mutate(Species=substr(Specimen.No.,1,2),
                   Date=as.POSIXct(Sampling.date,format="%d/%m/%Y"),
                   Year=year(Date),
                   Month=month(Date),
                   Sex=ifelse(Sex=='',NA,Sex),
                   Area=ifelse(grepl("Fatal Attraction",Comments),"Esperance","Albany"),
                   Gear=ifelse(grepl("LONGLINE",Comments),"Longline","Net"))%>%
            filter(!Species=="" & !Recorder=="")%>%
            left_join(IDL.to.FL.conversion,by="Species")%>%
            mutate(FL=ifelse(is.na(FL),IDL*slope+intercept,FL))



# Plots-------------------------------------------------------------------------

#reconstructed FL
n.tab=dat%>%group_by(Species)%>%tally()%>%rename(N=n)
dat%>%
  mutate(FL.5=5*round(FL/5))%>%
  group_by(Species,FL.5)%>%
  tally()%>%
  left_join(n.tab,by='Species')%>%
  mutate(LBL=paste( Species," (n=",N,")",sep=""))%>%
  ggplot(aes(x=FL.5,y=n,fill=Species)) + 
  geom_bar(position="dodge", stat="identity")+
  facet_wrap(~LBL)
ggsave(handl_OneDrive("Analyses/Samples/vertebrae_processor_Reconstructed_FL.tiff"),width = 10,height = 8,compression = "lzw")


#sample summaries
Tab.age.obs_yr=dat%>%
  group_by(Species,Year)%>%
  tally%>%
  spread(Species,n,fill=0)
write.csv(Tab.age.obs_yr,
          handl_OneDrive("Analyses/Samples/vertebrae_processor_samples by year_indicators.csv"),row.names = F)



fun.plt=function(d)
{
  d=d[order(d$Date),]
  attach(d)
  par(mfcol=c(3,2),mar=c(2,2,1,1),oma=c(1,1,1,1),mgp=c(1,.5,0))
  barplot(table(Vessel,useNA = 'ifany'),main='vessels',col=2)
  barplot(table(Date,useNA = 'ifany'),main='Date',col=2)
  barplot(table(Sex,useNA = 'ifany'),main='Sex',col=2)
  barplot(table(10*round(IDL/10),useNA = 'ifany'),main=paste('IDL (n=',nrow(d),')',sep=""),col=2)
  barplot(table(Area,useNA = 'ifany'),main='Area',col=2)
  barplot(table(Gear,useNA = 'ifany'),main='Gear',col=2)
  mtext("Frequency",2,outer=T,cex=1.5,line=-0.35)
  mtext(unique(Species),3,outer=T,cex=1.5,line=-1)
  detach(d)
}

sp=unique(dat$Species)
pdf(handl_OneDrive("Analyses/Samples/vertebrae_processor_summary.pdf"))
for(s in 1:length(sp)) fun.plt(d=subset(dat,Species==sp[s]))
dev.off()