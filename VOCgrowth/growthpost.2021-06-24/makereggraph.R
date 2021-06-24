library(tidyverse)
library(lubridate)
library(plyr)

lad2reg=read.csv('Local_Authority_District_to_Region__December_2019__Lookup_in_England.csv')
sanger=read.csv('lineages_by_ltla_and_week.tsv',sep='\t')

a=merge(sanger, lad2reg, by.x="LTLA", by.y="LAD19CD")# Could also use dplyr::left_join. See https://stackoverflow.com/questions/21888910/how-to-specify-names-of-columns-for-x-and-y-when-joining-in-dplyr
b=a %>% filter(WeekEndDate>="2021-04-01") %>% filter(Lineage=="B.1.1.7" | Lineage=="B.1.617.2")
# Using || instead of | in the above wouldn't work, perhaps because || does something odd with vectors. See: x=c(1:10);x>8||x<5;x>8|x<5. https://stackoverflow.com/a/22251971/7881139
d=aggregate(b$Count, by=list(Variant=b$Lineage, Date=ymd(b$WeekEndDate), Region=b$RGN19NM), FUN=sum)
e=spread(d,Variant,x)
names(e)[names(e) == "B.1.1.7"] <- "Alpha"
names(e)[names(e) == "B.1.617.2"] <- "Delta"
minfitdate=ymd("2021-05-15")
lme=function(df){
  m=lm(log(`Delta`/`Alpha`) ~ Date, df %>% filter(Date>=minfitdate), weights=1/(1/`Alpha`+1/`Delta`))
  a=coef(m)[1];b=coef(m)[2]
  data.frame(intercept=a,slope=b,crossover=as.Date("1970-01-01")-a/b+.5,instgrowth=format(b,digits=3),fivedaygrowth=format((exp(5*b)-1)*100,digits=2))
}
fit <- ddply(e, .(Region), lme)

ggplot(data = e, aes(x = Date, y = log(`Delta`/`Alpha`))) + geom_point(na.rm=T) + geom_abline(data=fit,aes(intercept=intercept,slope=slope)) + geom_text(data=fit,x=ymd("2021-04-23"),y=2.4,label=paste("log(Delta/Alpha) using Sanger data\nLine fits to values from ",minfitdate,"\n","Continuous growth rate: ",fit$instgrowth,"/day\nFive day growth ~= Rt(δ)/Rt(α): ",fit$fivedaygrowth,"%",sep="")) + facet_wrap(~Region)
f=0.75
ggsave("growthratesbyregion.png",units="in",width=20*f,height=16*f,dpi=80/f)
