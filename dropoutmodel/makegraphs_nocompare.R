library(tidyverse)
library(lubridate)

mydat=gather(read_csv("dropoutmodel.csv"), "place", "sgtf", -date)

onsdat=read_csv("ons_ct.csv") %>% mutate(date=dmy(`Week started`))# Only use this to read the last date
lastdate=max(onsdat$date)

mindate="2020-10-01"
mydat=mydat %>% filter(date >= mindate)

ggplot()+geom_line(data=mydat,aes(x=date,y=log(sgtf/(1-sgtf)),colour=place,group=place))+labs(x="Date",y="Log odds estimated B.1.1.7")+scale_x_date(date_breaks="months",date_labels="%b %y")+scale_y_continuous()+coord_cartesian(ylim=c(-5,NA))+geom_vline(xintercept=lastdate,linetype=3)+theme(legend.title=element_blank())
ggsave("sgtf.logodds.png",units="in",width=10,height=7,dpi=180)

ggplot()+geom_line(data=mydat,aes(x=date,y=sgtf,colour=place,group=place))+scale_x_date(date_breaks="months",date_labels="%b %y")+scale_y_continuous(labels=scales::percent)+geom_vline(xintercept=lastdate,linetype=3)+theme(legend.title=element_blank())
ggsave("sgtf.percent.png",units="in",width=10,height=7,dpi=180)
