library(tidyverse)
library(lubridate)

mydat=gather(read_csv("dropoutmodel.csv"), "nhs_name", "sgtf", -date)
phedat=read_csv("sgtf-2021-01-18.csv")
phebrief = rename(read_csv("PHEbriefing7.csv") %>% mutate(date=dmy(`week`)), nhs_name=Region) %>% mutate(nhs_name=replace(nhs_name,nhs_name=="Yorkshire and Humber","Yorkshire and The Humber"))

onsdat=read_csv("ons_ct.csv") %>% mutate(date=dmy(`Week started`))# Only use this to read the last date
lastdate=max(onsdat$date)

mindate="2020-10-01"
mydat=mydat %>% filter(date >= mindate)
phedat=phedat %>% filter(date >= mindate)
phebrief=phebrief %>% filter(date >= mindate)

# mydat, phebrief uses seven NHS regions: https://www.arcgis.com/sharing/rest/content/items/e4dd34a2d07b43088e17ed73e1bd3357/data
# phedat uses nine regions: https://www.arcgis.com/sharing/rest/content/items/738566d958834b70a809dc72a154eaf9/data
# To compare, restrict to a common set
# Note that "North West" is a bit different between the two despite having the same name, but I'm including it anyway
common=intersect(unique(pull(phedat,"nhs_name")), unique(pull(mydat,"nhs_name")))

ggplot()+geom_line(data=phedat %>% filter(sgtf+other>=200 & nhs_name %in% common), aes(x=date,y=log(sgtf/other),colour="PHE Pillar 2 via LSHTM",group=nhs_name))+geom_line(data=mydat %>% filter(nhs_name %in% common),aes(x=date,y=log(sgtf/(1-sgtf)),colour="Modelled from ONS survey data",group=nhs_name))+geom_line(data=phebrief %>% filter(nhs_name %in% common),aes(x=date,y=log(`n_Cases with confirmed SGTF`/`n_Cases with confirmed S-gene`),colour="PHE technical briefing 7",group=nhs_name))+labs(x="Date",y="Log odds estimated B.1.1.7")+scale_y_continuous()+coord_cartesian(ylim=c(-5,NA))+geom_vline(xintercept=lastdate,linetype=3)+theme(legend.title=element_blank())+facet_wrap(~nhs_name)
ggsave("sgtf.compare3.logodds.png",units="in",width=10,height=7,dpi=120)

ggplot()+geom_line(data=mydat,aes(x=date,y=log(sgtf/(1-sgtf)),colour="Modelled from ONS survey data",group=nhs_name))+geom_line(data=phebrief,aes(x=date,y=log(`n_Cases with confirmed SGTF`/`n_Cases with confirmed S-gene`),colour="PHE technical briefing 7",group=nhs_name))+labs(x="Date",y="Log odds estimated B.1.1.7")+scale_y_continuous()+coord_cartesian(ylim=c(-5,NA))+geom_vline(xintercept=lastdate,linetype=3)+theme(legend.title=element_blank())+facet_wrap(~nhs_name)
ggsave("sgtf.compare2.logodds.png",units="in",width=10,height=7,dpi=120)


ggplot()+geom_line(data=phedat %>% filter(sgtf+other>=200 & nhs_name %in% common), aes(x=date,y=sgtf/(sgtf+other),colour="PHE Pillar 2 via LSHTM",group=nhs_name))+geom_line(data=mydat %>% filter(nhs_name %in% common),aes(x=date,y=sgtf,colour="Modelled from ONS survey data",group=nhs_name))+geom_line(data=phebrief %>% filter(nhs_name %in% common),aes(x=date,y=`percent_SGTF`/100,colour="PHE technical briefing 7",group=nhs_name))+labs(x="Date",y="Estimated percentage B.1.1.7")+scale_y_continuous(labels=scales::percent)+geom_vline(xintercept=lastdate,linetype=3)+theme(legend.title=element_blank())+facet_wrap(~nhs_name)
ggsave("sgtf.compare3.percent.png",units="in",width=10,height=7,dpi=120)

ggplot()+geom_line(data=mydat,aes(x=date,y=sgtf,colour="Modelled from ONS survey data",group=nhs_name))+geom_line(data=phebrief,aes(x=date,y=`percent_SGTF`/100,colour="PHE technical briefing 7",group=nhs_name))+labs(x="Date",y="Estimated percentage B.1.1.7")+scale_y_continuous(labels=scales::percent)+geom_vline(xintercept=lastdate,linetype=3)+theme(legend.title=element_blank())+facet_wrap(~nhs_name)
ggsave("sgtf.compare2.percent.png",units="in",width=10,height=7,dpi=120)
