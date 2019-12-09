#Script for tallying collection of vertebrae from fish processors

#best way of identifying logbook trip:
#Nick: "Re the trip info, are you meaning the logbook returns which we could get from Anja (fishing block, weights etc), or VMS, where we would put in a data request form for each trip. That would give us gps of where they went fishing.
# Sampling date from the datasheet would be a day or 2 ahead of the trip finish date"
#Tim: "logbook return daily sheets (with block numbers) will be ahead of unload dates , by many days in some instances. Available from Anja.
#   VMS data won't distinguish which fish caught when, but in combination with logbook might narrow location down if different spp. caught. Will show discrete areas fished and allow us to remove areas where just transiting. Available by request (or from friendly FMO with a few minutes to spare?"



library(dplyr)
library(lubridate)

#Read in data
setwd('C:\\Matias\\Analyses\\Samples\\vertebrae collection')
dat=read.csv('U:/Shark/Fish_processor_age_sampling.csv',stringsAsFactors = F)
#dat=read.csv('Fish_processor_age_sampling.csv',stringsAsFactors = F)  #offline

#Manipulate data
dat=dat %>% mutate(Species=substr(Specimen.No.,1,2),
                   Date=as.POSIXct(Sampling.date,format="%d/%m/%Y"),
                   Year=year(Date),
                   Month=month(Date),
                   Sex=ifelse(Sex=='',NA,Sex),
                   Area=ifelse(grepl("Fatal Attraction",Comments),"Esperance","Albany"),
                   Gear=ifelse(grepl("LONGLINE",Comments),"Longline","Net"))%>%
            filter(!Species=="" & !Recorder=="")

#plots
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
pdf("summary.pdf")
for(s in 1:length(sp)) fun.plt(d=subset(dat,Species==sp[s]))
dev.off()